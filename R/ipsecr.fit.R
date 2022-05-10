###############################################################################
## package 'ipsecr'
## ipsecr.fit.R
## 2022-04-01,18,19, 2022-05-08
###############################################################################

ipsecr.fit <- function (
    capthist, 
    proxyfn = proxyfn1, 
    model = list(D~1, g0~1, sigma~1), 
    mask = NULL,
    buffer = 100, 
    detectfn = 'HN', 
    binomN = NULL,
    start = NULL, 
    link = list(),
    fixed = list(),
    timecov = NULL,
    details = list(), 
    verify = TRUE, 
    verbose = TRUE, 
    ncores = NULL, 
    seed = NULL,
    ...) {
    
    ## ... passed to proxyfn
    ## boxsize may be vector of length NP
    ## proxyfn1 is default defined separately
    
    if (ms(capthist)) {
        stop ("'ipsecr.fit' not implemented for multi-session 'capthist'")
    }
    
    ptm  <- proc.time()
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    cl <- match.call(expand.dots = TRUE)
    cl <- paste(names(cl)[-1],cl[-1], sep=' = ', collapse=', ' )
    cl <- paste('ipsecr.fit(', cl, ')')
    code <- 0  # exit code
    
    #################################################
    # optionally replace detector type
    if (!is.null(details$newdetector)) {
        warning("replacement detector type specified by user")
        detector(traps(capthist)) <- details$newdetector       
    }
    #################################################
    
    
    #################################################
    ## inputs
    #################################################
    
    traps <- traps(capthist)
    noccasions <- ncol(capthist)
    ncores <- setNumThreads(ncores)
    if (is.null(mask)) {
        mask <- make.mask(traps, buffer = buffer, nx = 64)
    }
    
    #################################################
    ## detection function
    #################################################
    
    if (is.character(detectfn)) {
        detectfn <- detectionfunctionnumber(detectfn)
    }
    if (!(detectfn %in% c(0,2,4, 14,16))) {
        stop (detectionfunctionname(detectfn),
            " detection function not implemented in ipsecr.fit")
    }
    
    #################################################
    
    #################################################
    defaultdetails <- list(
        boxsize      = 0.2, 
        boxsize2     = 0.05, 
        boxtype      = 'absolute',
        centre       = 3,
        dev.max      = 0.002, 
        var.nsim     = 2000,
        min.nsim     = 20,
        max.nsim     = 2000, 
        max.nbox     = 5, 
        max.ntries   = 2,
        distribution = 'poisson',
        even         = FALSE,
        binomN       = 0,               ## Poisson counts
        param        = 0,
        ignoreusage  = FALSE,
        ignorenontarget = FALSE,
        debug        = FALSE,
        savecall     = TRUE,
        newdetector  = NULL,
        contrasts    = NULL,
        CHmethod     = 'internal',
        popmethod    = 'internal',
        Nmax         = 1e4,
        factorial    = 'full',
        FrF2args     = NULL
    )
    if (any(!(names(details) %in% names(defaultdetails)))) {
        stop ("details list includes invalid name")
    }
    details <- replace (defaultdetails, names(details), details)
    if (details$max.nbox<2) stop("ipsecr.fit details$max.nbox >= 2")
    details$distribution <- match.arg(details$distribution, choices = c('poisson','binomial','even'))
    details$boxtype <- match.arg(details$boxtype, choices = c('absolute','relative'))
    details$popmethod <- match.arg(details$popmethod, choices = c('internal','sim.popn'))
    if (!is.function(details$CHmethod)) {
        details$CHmethod <- match.arg(details$CHmethod, choices = c('internal','sim.capthist'))
    }
    # choices for factorial depend on FrF2
    details$factorial<- match.arg(details$factorial, choices=c('full','fractional'))
    details$verbose <- verbose
    
    #################################################
    ## optional data check
    #################################################
    if (verify) {
        memo ('Checking data', verbose)
        test <- verify(capthist, report = 1)
        if (test$errors)
            stop ("'verify' found errors in 'capthist' argument")
        
        if (!is.null(mask)) {
            notOK <- verify(mask, report = 1)$errors
            if (notOK)
                stop ("'verify' found errors in 'mask' argument")
        }
    }
    
    #################################################
    ## nontarget interference 
    #################################################
    
    nontarget <- attr(capthist, 'nontarget', exact = TRUE)
    modelnontarget <- !is.null(nontarget) && !details$ignorenontarget
    ## check
    if (modelnontarget) {
        # this duplicates a check in secr::verify from 4.5.5
        ks <- apply(abs(capthist), 2:3, sum)
        if (any(t(apply(capthist,2:3,sum)) & nontarget)) {
            stop ("nontarget conflicts with capture at least once")
        }
        if (isTRUE(all.equal(proxyfn, proxyfn1))) {
            proxyfn <- proxyfn2
            warning("replacing default proxy function with proxyfn2 for nontarget model")
        }
    }
    
    #################################################
    ## standardize user model and parameterisation
    #################################################
    
    if ('formula' %in% class(model)) model <- list(model)
    model <- stdform (model)  ## named, no LHS
    model <- updatemodel(model, detectfn, 14:20, 'g0', 'lambda0')
    
    if (any(sapply(model[-1], "!=", ~1))) stop ("not ready for varying detection")
    
    fnames <- names(fixed)
    
    #################################################
    ## build default model and update with user input
    #################################################
    
    defaultmodel <- list(D=~1, g0=~1, lambda0=~1, sigma=~1, z=~1, w=~1, lambdak=~1)
    defaultmodel <- replace (defaultmodel, names(model), model)
    
    #################################################
    ## parameter names
    #################################################
    
    pnames <- valid.pnames (details, FALSE, detectfn, FALSE, FALSE, 1)
    if (modelnontarget) pnames <- c(pnames, 'lambdak')
    
    #################################################
    ## test for irrelevant parameters in user's model
    #################################################
    
    OK <- names(model) %in% pnames
    if (any(!OK))
        stop ("parameters in model not consistent with detectfn etc. : ",
            paste(names(model)[!OK], collapse = ', '))
    OK <- fnames %in% pnames
    if (any(!OK))
        stop ("attempt to fix parameters not in model : ",
            paste(fnames[!OK], collapse = ', '))
    
    #################################################
    ## finalise model
    #################################################
    
    pnames <- pnames[!(pnames %in% fnames)]   ## drop fixed real parameters
    model <- defaultmodel[pnames]             ## select real parameters
    # valid.model(model, FALSE, detectfn, NULL, NULL, '')
    vars <-  unlist(lapply(model, all.vars))
    
    #################################################
    # Link functions (model-specific)
    #################################################
    
    defaultlink <- list(D = 'log', g0 = 'logit', lambda0 = 'log', sigma = 'log', lambdak = 'log')
    link <- replace (defaultlink, names(link), link)
    link[!(names(link) %in% pnames)] <- NULL
    
    ##############################################
    # Prepare detection design matrices and lookup
    ##############################################
    # memo ('Preparing detection design matrices', verbose)
    design <- secr.design.MS (capthist, model, timecov, NULL, NULL, NULL,
        NULL, ignoreusage = details$ignoreusage, naive = FALSE,
        CL = FALSE, contrasts = details$contrasts)
    design0 <- design   # for now
    # design0 <- secr.design.MS (capthist, model, timecov, NULL, NULL, NULL,
    #     NULL, ignoreusage = details$ignoreusage, naive = TRUE,
    #     CL = FALSE, contrasts = details$contrasts)
    
    ############################################
    # Prepare density design matrix
    ############################################
    
    D.modelled <- is.null(fixed$D)
    if (!D.modelled) {
        designD <- matrix(nrow = 0, ncol = 0)
        grouplevels <- 1    ## was NULL
        attr(designD, 'dimD') <- NA
        nDensityParameters <- integer(0)
    }
    else {
        memo ('Preparing density design matrix', verbose)
        temp <- D.designdata( mask, model$D, 1, 1, NULL)
        designD <- general.model.matrix(model$D, data = temp, gamsmth = NULL,
            contrasts = details$contrasts)
        attr(designD, 'dimD') <- attr(temp, 'dimD')
        Dnames <- colnames(designD)
        nDensityParameters <- length(Dnames)
    }

    ############################################
    # Parameter mapping (general)
    ############################################
    nDetectionParameters <- sapply(design$designMatrices, ncol)
    np <- c(D = nDensityParameters, nDetectionParameters)
    NP <- sum(np)
    parindx <- split(1:NP, rep(1:length(np), np))
    names(parindx) <- names(np)[np>0]
    
    ###########################################
    # setup boxes
    ###########################################
    
    if (length(details$boxsize)==1) boxsize1 <- rep(details$boxsize, NP)
    else if (length(boxsize1) != NP)
        stop ("invalid boxsize vector")
    else boxsize1 <- details$boxsize
    if (length(details$boxsize2)==1) boxsize2 <- rep(details$boxsize2, NP)
    else if (length(details$boxsize2) != NP)
        stop ("invalid boxsize vector")
    else boxsize2 <- details$boxsize2

    ###########################################
    # cluster for parallel processing
    ###########################################
    
    if (ncores > 1) {
        memo ('Preparing cluster for parallel processing', verbose)
        if(.Platform$OS.type == "unix") {
            clust <- makeCluster(ncores, type = "FORK")
        }
        else {
            clust <- makeCluster(ncores, type = "PSOCK")
            clusterExport(clust, c(
                "mask", "link", "fixed", "details", "traps",
                "detectfn", "noccasions", "proxyfn",
                "parindx", "designD", "getD", "simCH",
                "modelnontarget"),
                environment())
        }
        on.exit(stopCluster(clust))
        clusterSetRNGStream(clust, seed)
    }
    else {
        clust <- NULL
        set.seed(seed)
    }
    
    #---------------------------------------------------------------------------
    # to simulate one realization
    #---------------------------------------------------------------------------
    
    simfn <- function (beta, distribution = 'binomial', ...) {
        attempts <- 0
        allOK <- FALSE
        D <- getD(designD, beta, mask, parindx, link, fixed)
        N <- sum(D) * attr(mask, 'area')
        Ndist <- switch(distribution, poisson = 'poisson', binomial = 'fixed', even = 'fixed')    
        N <- switch(tolower(Ndist), poisson = rpois(1, N), fixed = round(N), NA)
        # cannot simulate zero animals, so return NA for predicted
        if (is.na(N) || N<=0 || N>details$Nmax) return(rep(NA, NP))
        
        if (details$popmethod == 'internal') {
            prob <- D/sum(D)
        }
        else if (details$popmethod == 'sim.popn') {
            covariates(mask) <- data.frame(D = D)
        }
        
        # detectpar is list of named parameters of detectfn
        # not yet modelled using design, design0: assumes 1:1
        betalist <- lapply(parindx[-1], function(x) unname(beta[x]))
        detectpar <- mapply(untransform, betalist, link[-1], SIMPLIFY = FALSE)
        
        repeat {
            #------------------------------------------------
            # generate population
            if (details$popmethod == 'internal') {   # C++
                if (distribution == 'even') {
                    bounds <- apply(mask,2,range)
                    popn <- popevencpp(as.matrix(bounds), as.integer(N))
                }
                else {
                    # cat('N ', N, ' sum(prob) ', sum(prob), '\n')
                    popn <- popcpp(as.matrix(mask), as.double(prob), 
                        as.double(spacing(mask)/100), as.integer(N))
                }
            }
            else if (details$popmethod == "sim.popn") {
                popn <- sim.popn (D = 'D', core = mask, model2D = 'IHP', Ndist = Ndist)
            }
            # abort if no animals to sample
            if (length(popn)<2) {
                warning ("ipsecr.fit: no animals in simulated population", call. = FALSE)
                return(rep(NA, NP))
            }
            
            #------------------------------------------------
            # sample from population
            if (is.function(details$CHmethod)) {   # user
                ch <- details$CHmethod(traps, popn, detectfn, detectpar, noccasions)
            }
            else if (details$CHmethod == 'internal') {   # C++
                ch <- simCH(traps, popn, detectfn, detectpar, noccasions)
            }
            else if (details$CHmethod == 'sim.capthist') {   
                if (modelnontarget) stop ("sim.capthist does not simulate nontarget detections")
                ch <- sim.capthist(traps, popn, detectfn, detectpar, noccasions)
            }
            
            # check valid
            if (nrow(ch)==0) {
                warning ("ipsecr.fit: no captures in simulation", call. = FALSE)
            }
            else {
                if ((sum(abs(ch)>0) - nrow(ch)) < 1)
                    warning ("ipsecr.fit: no re-captures in simulation", call. = FALSE)
            }
            
            #------------------------------------------------
            # predict values
            predicted <- try (proxyfn (ch, ...))
            if (inherits(predicted, 'try-error')) {
                predicted <- rep(NA, NP)
            }
            if (details$debug) {
                cat ('saving ch to ch.RDS\n')
                saveRDS(ch, file='ch.RDS')
                cat('N ', N, '\n')
                cat('detectpar', unlist(detectpar), '\n')
                cat('detectfn', detectfn, '\n')
                print(summary(ch))
                cat('\nproxyfn\n')
                print(proxyfn)
                cat('\npredicted', predicted, '\n')
                stop()
            }
            attempts <- attempts+1
            
            #------------------------------------------------
            ## exit loop if exceeded max.ntries or all OK
            allOK <- !any(is.na(predicted)) && all(is.finite(predicted))
            if (attempts >= details$max.ntries || allOK)
                break
        }
        #----------------------------------------------------
        
        if (!allOK) {
            warning ("ipsecr.fit: no successful simulation after ", details$max.ntries,
                " attempts", call. = FALSE)
            return (rep(NA, NP))
        }
        else {
            if (attempts > 1)
                warning ("ipsecr.fit: simulation repeated", call. = FALSE)
            return(predicted)
        }
    }   # end of simfn()
    
    ##########################################
    ## function to test if current solution is within box (vector result)
    ##########################################
    within <- function (i) {
        beta[i] >= vertices[[i]][1] & 
            beta[i] <= vertices[[i]][2]
    }
    
    ##########################################
    ## target values of predictor
    ##########################################
    y <- proxyfn(capthist, ...)
    if (length(y) != NP)
        stop ("need one proxy for each parameter ",
            paste(pnames, collapse=" "))
    
    ##########################################
    ## starting values
    ##########################################
    ## ad hoc exclusion of lambdak 2022-05-06
    details$trace <- verbose  # as used by makeStart
    start <- makeStart(start, parindx[names(parindx) != 'lambdak'], 
        capthist, mask, detectfn, link, details, fixed)
    if (modelnontarget) {
        start <- c(start, y[4])
    }
    
    ############################################
    # Fixed beta parameters
    ############################################
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        if (!(length(fb)== NP))
            stop ("invalid fixed beta - require NP-vector")
        if (sum(is.na(fb))==0)
            stop ("cannot fix all beta parameters")
        ## drop unwanted betas; remember later to adjust parameter count
        start <- start[is.na(fb)]
        NP <- length(start)
    }
    
    ############################################
    # Variable names (general)
    ############################################
    betanames <- unlist(sapply(design$designMatrices, colnames))
    names(betanames) <- NULL
    realnames <- names(model)
    ## coefficients for D precede all others
    betanames <- c(paste('D', Dnames, sep='.'), betanames)
    betanames <- sub('..(Intercept))','',betanames)
    
    ####################################################################
    ## Loop over boxes
    ## keep trying until estimates within box or number exceeds max.nbox
    ####################################################################
    
    beta <- start
    
    for (m in 1:details$max.nbox) {
        if (verbose) {
            cat('\nFitting box', m, '...  \n')
        }
        boxsize <- if (m == 1) boxsize1 else boxsize2
        if (details$boxtype == 'relative') {
            vertices <- sweep (1 + outer(c(-1,1), boxsize), MARGIN = 2,
                FUN = '*', STATS = beta)
        }
        else {
            vertices <- sweep (outer(c(-1,1), boxsize), MARGIN = 2,
                FUN = '+', STATS = beta)
        }
        vertices <- data.frame(vertices)
        names(vertices) <- betanames
        rownames(vertices) <- c('min','max')
        if (verbose) {
            print(vertices)
            cat('\n')
            flush.console()
        }
        
        if (details$factorial == 'full') {
            # full factorial
            designbeta <- as.matrix(expand.grid (as.list(vertices)))
            centrepoints <- matrix(beta, nrow = details$centre, ncol = NP, byrow = T)
            designbeta <- rbind(designbeta, centrepoints)
        }
        else {
            if (!requireNamespace('FrF2'))
                stop ('package FrF2 required for fractional factorial')
            if (is.null(details$FrF2args)) {
                details$FrF2args <- list(nruns=2^(NP-1), nfactors=NP, ncenter = details$centre)
            }
            Fr <- do.call(FrF2::FrF2, details$FrF2args)
            # recast factors as numeric
            base <- sapply(Fr, function(x) as.numeric(as.character(x)))
            base <- sweep(base, MARGIN = 2, STATS = boxsize, FUN = '*')
            if (details$boxtype == "absolute") {
                designbeta <- sweep(base, MARGIN = 2, STATS = beta, FUN = '+')
            }
            else {
                designbeta <- sweep(base, MARGIN = 2, STATS = beta, FUN = '*')
            }
        }
        
        designpoints <- nrow(designbeta)
        sim <- NULL
        alldesignbeta <- NULL
        tempdistn <- if (details$even) 'even' else 'binomial'
        basedesign <- designbeta[rep(1:nrow(designbeta), details$min.nsim),]
        # accumulate simulations until reach precision target or exceed max.nsim
        tries <- 0
        repeat {
            dev <- 0   # criterion for precision (SE, RSE) before set
            if (ncores > 1) {
                list(...) # evaluate any promises cf boot
                newsim <- parRapply(clust, basedesign, simfn, distribution = tempdistn, ...)
                newsim <- t(matrix(newsim, ncol = nrow(basedesign)))
            }
            else {
                newsim <- t(apply(basedesign,1,simfn, distribution = tempdistn, ...))
            }
            OK <- apply(!apply(newsim,1, is.na), 2, all)
            if (sum(OK) == 0) {
                code <- 3
                break
            }
            sim <- rbind(sim, newsim[OK,])
            alldesignbeta <- rbind(alldesignbeta,basedesign[OK,])
            if ((nrow(alldesignbeta) > details$max.nsim)) {
                break
            }
            sim.lm <- lm ( sim ~ alldesignbeta )   # proxy ~ link(param)
            sum.sim.lm <- summary(sim.lm)
            code <- 0
            if (details$boxtype == 'absolute') {
                dev <- sapply(sum.sim.lm, function(x) x$sigma) / sqrt(nrow(sim))
            }
            else {
                dev <- sapply(sum.sim.lm, function(x) x$sigma) / y / sqrt(nrow(sim))
            }
            # break if have achieved target precision
            if (!is.null(dev) && !any(is.na(dev)) && all(dev <= details$dev.max)) break
        }
        
        if (any(dev > details$dev.max)) {
            crit <- switch(details$boxtype, absolute = 'SE', relative = 'RSE')
            memo(paste0("simulations for box ", m, " did not reach target for proxy ", 
                crit, " ", details$dev.max), verbose)
        }
        
        if (code>2) {
            beta <- rep(NA, NP)
            if (code == 3) warning ("no valid simulations")
            # if (code == 4) warning ("exceeded maximum allowable replicates ",
            #     "without achieving precision better than 'dev.max'")
            # if (code == 5) warning ("ipsecr.fit: invalid lm")
        }
        else {
            B <- coef(sim.lm)[-1,]
            B <- solve(t(B))  ## invert
            lambda <- coef(sim.lm)[1,]   ## intercepts
            beta <- as.numeric(B %*% matrix((y - lambda), ncol = 1))
            ## only break on second or later box if differ boxsize
            if (all(sapply(1:NP, within)) && (all(boxsize == boxsize2) || (m>1))) break
        }
    }    # end of loop over boxes
    ####################################################################
    
    if (code == 0) {
        if (!all(sapply(1:NP, within))) {
            warning ("solution not found after ", details$max.nbox, " attempts")
            code <- 2
        }
        else {
            code <- 1
        }
    }
    ####################################################################
    
    if (details$var.nsim>1 && code == 1) {
        if (verbose) {
            cat('Simulating for variance ...\n')
            flush.console()
            cat('\n')
        }
        
        vardesign <- matrix(beta, nrow = details$var.nsim, ncol = NP, byrow = T)
        colnames(vardesign) <- betanames
        if (ncores > 1) {
            list(...) # evaluate any promises cf boot
            newsim <- parRapply(clust, vardesign, simfn, distribution = details$distribution, ...)
            newsim <- t(matrix(newsim, ncol = nrow(vardesign)))
        }
        else {
            newsim <- t(apply(vardesign,1,simfn, distribution = details$distribution, ...))
        }
        OK <- apply(!apply(newsim,1, is.na), 2, all)
        memo(paste(sum(OK), " variance simulations successful\n"), verbose)
        if (sum(OK) != details$var.nsim) {
            warning(paste(details$var.nsim-sum(OK), "of", details$var.nsim, "variance simulations failed"))
        }
        newsim <- newsim[OK,]
        V <- var(newsim)  ## additional simulations for var-covar matrix
        vcov <- B %*% V %*% t(B)
        
        ## compare estimates to parametric bootstrap
        n <- apply(newsim, 2, function(x) sum(!is.na(x)))
        ymean <- apply(newsim, 2, mean, na.rm=T)
        yse <- apply(newsim, 2, function(x) sd(x, na.rm=T) / sum(!is.na(x)))
        
        bootstrap <- data.frame (target = y, nsim = n, simulated = ymean,
            SE.simulated = yse)
        
        ## biasest not reported, yet
        yest <- as.numeric(B %*% matrix((ymean - lambda), ncol = 1))
        biasest <- data.frame (estimate = 100 * (beta - yest) / yest,
            SE = 100 * (beta - yest) / yest)
        
    }
    else {
        vcov <- matrix(nrow = NP, ncol = NP)
        bootstrap <- NA
    }
    dimnames(vcov) <- list(betanames, betanames)
    
    desc <- packageDescription("ipsecr")  ## for version number
    
    output <- list(
        call = cl,
        capthist = capthist,
        mask = mask,
        detectfn = detectfn,
        timecov = timecov,
        start = start,
        link = link,
        fixed = fixed,
        model = model,
        details = details,
        designD = designD,        
        design = design,
        design0 = design0,
        parindx = parindx,
        vars = vars,
        betanames = betanames,
        realnames = realnames,
        code = code,                # success of fit
        beta = beta,
        beta.vcv = vcov,
        designbeta = designbeta,   # last 
        ip.nsim = nrow(sim),
        var.nsim.OK = if(details$var.nsim>1) sum(OK) else NA,
        variance.bootstrap = bootstrap,
        version = desc$Version,
        starttime = starttime,
        proctime = as.numeric((proc.time() - ptm)[3])
    )
    class(output) <- 'ipsecr'
    memo(paste('Completed in ', round(output$proctime,2), ' seconds at ',
        format(Sys.time(), "%H:%M:%S %d %b %Y")), verbose)
    output
}
##################################################

