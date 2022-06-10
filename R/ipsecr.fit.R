###############################################################################
## package 'ipsecr'
## ipsecr.fit.R
## 2022-04-01,18,19, 2022-05-08, 2022-06-07
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
    sessioncov = NULL,
    details = list(), 
    verify = TRUE, 
    verbose = TRUE, 
    ncores = NULL, 
    seed = NULL,
    ...) {
    
    ## ... passed to proxyfn
    ## boxsize may be vector of length NP
    ## proxyfn1 is default defined separately

    ptm  <- proc.time()
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    cl <- match.call(expand.dots = TRUE)
    cl <- paste(names(cl)[-1],cl[-1], sep=' = ', collapse=', ' )
    cl <- paste('ipsecr.fit(', cl, ')')
    code <- 0  # exit code

    #################################################
    ## inputs
    #################################################
    
    traps <- traps(capthist)
    #------------------------------------------
    # optionally replace detector type
    if (!is.null(details$newdetector)) {
        warning("replacement detector type specified by user")
        detector(traps) <- details$newdetector       
    }
    #------------------------------------------
    ncores <- setNumThreads(ncores)
    if (is.null(mask)) {
        mask <- make.mask(traps, buffer = buffer, nx = 64)
    }
    sessionlevels <- session(capthist)
    if (ms(capthist)) {
        nsessions <- length(capthist)
        noccasions <- sapply(capthist, ncol)
        if (details$popmethod == 'internal') {
            stop ("'internal' method does not simulate multi-session populations; try 'sim.popn'")
        }
        if (details$popmethod == 'sim.pop' && packageVersion('secr') < '4.5.6') {
            stop ("simulation of multi-session population with sim.pop requires secr >= 4.5.6")
        }
        if (!is.null(mask) && !ms(mask)) {
            ## inefficiently replicate mask for each session!
            mask <- lapply(sessionlevels, function(x) mask)
            class (mask) <- c('mask', 'list')
            names(mask) <- sessionlevels
        }
        cellarea <- sapply(mask, attr, 'area')
    }
    else {
        nsessions <- 1
        noccasions <- ncol(capthist)
        cellarea <- attr(mask, 'area')
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
        min.nbox     = 2, 
        max.nbox     = 5, 
        max.ntries   = 2,
        distribution = 'poisson',
        even         = FALSE,
        binomN       = 0,               ## Poisson counts
        param        = 0,
        ignoreusage  = FALSE,
        ignorenontarget = FALSE,
        nontargettype = 'exclusive',
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
    # nontarget model 
    #################################################

    validnontargettype <- c('exclusive', 'truncated','erased','independent')
    details$nontargettype <- match.arg(details$nontargettype, choices = validnontargettype)
    if (any(detector(traps) %in% c('multi', 'proximity', 'count'))) {
        if (details$nontargettype == 'exclusive') {
            details$nontargettype <- 'truncated'   # downgrade for these detectors
        }
    }
    
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
        if (isTRUE(all.equal(proxyfn, proxyfn1))) {
            proxyfn <- proxy.nt
            warning("replacing default proxy function with proxy.nt for nontarget model")
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
    design <- secr.design.MS (capthist, model, timecov, sessioncov, NULL, NULL,
        NULL, ignoreusage = details$ignoreusage, naive = FALSE,
        CL = FALSE, contrasts = details$contrasts)
    design0 <- design   # for now
    # design0 <- secr.design.MS (capthist, model, timecov, sessioncov, NULL, NULL,
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
        temp <- D.designdata( mask, model$D, 1, sessionlevels, sessioncov)
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

    #---------------------------------------------------------------------------
    # to simulate one realization
    #---------------------------------------------------------------------------
    
    simfn <- function (beta, distribution = 'binomial', ...) {
        NP <- length(beta)
        attempts <- 0
        allOK <- FALSE
        D <- getD(designD, beta, mask, parindx, link, fixed, nsessions)
        N <- apply(D,2,sum) * cellarea
        # cat('step 1 \n', apply(D,2,sum) , '\n', cellarea, '\n', distribution, '\n', N, '\n')        
        Ndist <- switch(distribution, poisson = 'poisson', binomial = 'fixed', even = 'fixed')    
        N <- switch(tolower(Ndist), poisson = rpois(1, N), fixed = round(N), NA)
        
        # cannot simulate zero animals, so return NA for predicted

        if (is.na(N) || any(N<=0) || any(N>details$Nmax)) return(rep(NA, NP))
        
        if (details$popmethod == 'internal') {
            prob <- sweep(D, STATS = apply(D,2,sum), FUN = '/', MARGIN = 2)
        }
        # else if (details$popmethod == 'sim.popn') {
        #     # covariates(mask) <- as.data.frame(D)
        # }
        
        # detectpar is list of named parameters of detectfn
        # not yet modelled using design, design0: assumes 1:1
        betalist <- lapply(parindx[-1], function(x) unname(beta[x]))
        detectpar <- mapply(untransform, betalist, link[-1], SIMPLIFY = FALSE)

        repeat {
            
            #------------------------------------------------
            # generate population
            if (is.function(details$popmethod)) {   # user
                popn <- details$popmethod(mask, prob, N)
            }
            else if (details$popmethod == 'internal') {   # C++
                if (nsessions>1) {
                    stop("internal not ready for multisession")
                }
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
                popn <- sim.popn (D = as.data.frame(D), core = mask, 
                    model2D = 'IHP', Ndist = Ndist, nsessions = nsessions)
            }
            # abort if no animals to sample
            if (length(popn)<2) {
                warning ("ipsecr.fit: no animals in simulated population", call. = FALSE)
                return(rep(NA, NP))
            }

            #------------------------------------------------
            # sample from population
            if (is.function(details$CHmethod)) {   # user
                ch <- details$CHmethod(traps, popn, detectfn, detectpar, noccasions, details)
            }
            else if (details$CHmethod == 'internal') {   # C++
                ch <- simCH(traps, popn, detectfn, detectpar, noccasions, details)
            }
            else if (details$CHmethod == 'sim.capthist') {   
                if (modelnontarget) stop ("sim.capthist does not simulate nontarget detections")
                ch <- sim.capthist(traps, popn, detectfn, detectpar, noccasions, nsessions, 
                    renumber = FALSE)
            }
            # check valid
            n <- if (nsessions == 1) nrow(ch) else sapply(ch,nrow)
            r <- if (nsessions == 1) sum(abs(ch)>0) else sapply(ch, function(x) sum(abs(x)>0))
            if (any(n==0)) {
                warning ("ipsecr.fit: no captures in simulation", call. = FALSE)
            }
            else {
                if (any((r - n) < 1))
                    warning ("ipsecr.fit: no re-captures in simulation", call. = FALSE)
            }

            #------------------------------------------------
            # predict values
            predicted <- try (proxyfn (ch, model = model, ...))
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
    y <- proxyfn(capthist, model = model, ...)
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

    ###########################################
    # cluster for parallel processing
    ###########################################
    
    if (ncores > 1) {
        memo ('Preparing cluster for parallel processing', verbose)
        if(.Platform$OS.type == "unix") {
            # clust <- makeCluster(ncores, type = "FORK", outfile = "clusterlogfile.txt")
            clust <- makeCluster(ncores, type = "FORK")
        }
        else {
            # clust <- makeCluster(ncores, type = "PSOCK", outfile = "clusterlogfile.txt")
            clust <- makeCluster(ncores, type = "PSOCK")
            clusterExport(clust, c(
                "mask", "link", "fixed", "details", "traps",
                "detectfn", "noccasions", "nsessions", "proxyfn",
                "model", "parindx", "designD", "getD", "simCH",
                "modelnontarget", "cellarea"), 
                environment())
        }
        on.exit(stopCluster(clust))
        clusterSetRNGStream(clust, seed)
    }
    else {
        clust <- NULL
        set.seed(seed)
    }

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
            if (!requireNamespace('FrF2', quietly = TRUE))
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
        
        # check number of simulations
        if (2 * designpoints * details$min.nsim > details$max.nsim) {
            stop("max.nsim = ", details$max.nsim, 
                " does not allow two boxes with designpoints = ", designpoints, 
                " and min.nsim = ", details$min.nsim)
        }
        
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
                newsim <- parRapply(clust, basedesign, simfn, 
                    distribution = tempdistn, ...)
                newsim <- t(matrix(newsim, ncol = nrow(basedesign)))
            }
            else {
                newsim <- t(apply(basedesign, 1, simfn, 
                    distribution = tempdistn, ...))
            }
            OK <- apply(!apply(newsim,1, is.na), 2, all)
            if (sum(OK) == 0) {
                code <- 3
                break
            }
            sim <- rbind(sim, newsim[OK,])
            alldesignbeta <- rbind(alldesignbeta,basedesign[OK,])
            sim.lm <- lm ( sim ~ alldesignbeta )   # proxy ~ link(param)
            sum.sim.lm <- summary(sim.lm)
            code <- 0
            if (details$boxtype == 'absolute') {
                dev <- sapply(sum.sim.lm, function(x) x$sigma) / sqrt(nrow(sim))
            }
            else {
                dev <- sapply(sum.sim.lm, function(x) x$sigma) / y / sqrt(nrow(sim))
            }
            # break if have exceeded allowed simulations
            if ((nrow(alldesignbeta) > details$max.nsim)) break
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
            # B <- solve(t(B))  ## invert
            B <- try(MASS::ginv(t(B)), silent = TRUE)   ## 2022-05-15
            if (inherits(B, 'try-error')) {
                code <- 5
                beta <- rep(NA, NP)
            }
            else {
                lambda <- coef(sim.lm)[1,]   ## intercepts
                beta <- as.numeric(B %*% matrix((y - lambda), ncol = 1))
                ## only break on second or later box if differ boxsize
                # if (all(sapply(1:NP, within)) && (all(boxsize == boxsize2) || (m>1))) break
                if (all(sapply(1:NP, within)) && (m >= details$min.nbox)) break
            }
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
            newsim <- parRapply(clust, vardesign, simfn, 
                distribution = details$distribution, ...)
            newsim <- t(matrix(newsim, ncol = nrow(vardesign)))
        }
        else {
            newsim <- t(apply(vardesign,1,simfn, 
                distribution = details$distribution, ...))
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
        proxyfn = proxyfn,
        model = model,
        mask = mask,
        detectfn = detectfn,
        start = start,
        link = link,
        fixed = fixed,
        timecov = timecov,
        sessioncov = sessioncov,
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

