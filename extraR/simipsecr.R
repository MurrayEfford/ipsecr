############################################################################################
## package 'ipsecr'
## simipsecr.R
## uses Dsurface.R, utility.R
## 2022-04-02 UNUSED
############################################################################################

simpop <- function (designD, beta, mask, parindx, link, fixed,
    vars, model, details, nsim = 1, seed = NULL, ...)
{
    ##  check input
    if (any(c("bn", "bkn", "bkc", "Bkc") %in% tolower(vars)))
        stop ("simulate works only with binary behavioural responses")
    
    D <- getD(designD, beta, mask, parindx, link, fixed)
    
    ##################
    ## set random seed
    ## copied from simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    ##################
    
    ## loop over replicates
    runone <- function(r) {
        ## argument r is unused dummy
        
        popn <- list()
        if (model$D == ~1) {
            density <- D[1]         ## scalar
            mod2D <- 'poisson'      ## homogeneous
        }
        else {
            density <- D            ## vector
            mod2D <- 'IHP'          ## inhomogeneous
        }
        
        Ndist <- switch (details$distribution,
            binomial = 'fixed',
            poisson = 'poisson',
            'poisson')
        sim.popn (density, core = mask, model2D = mod2D, Ndist = Ndist)
    }        
    poplist <- lapply(1:nsim, runone)
    
    attr(poplist,'seed') <- RNGstate   ## save random seed
    class(poplist) <- c('secrdata', 'list')
    poplist
}
############################################################################################

simCHslow <- function (
    beta, mask, parindx, model, link, fixed, design, design0, 
    traps, noccasions, details, detectfn, vars, popnlist, timecov,
    maxperpoint = 100) {

    ## popnlist is always a list of popn objects
    
    ## we use fake CH to extract parameter value dependent on previous capt
    ## BUT this is slow and clunky...
    
    dummycapthist <- function (pop, fillvalue = 1) {
        if (length(pop)>1) {
            output <- list()
            for (i in 1:length(pop))
                output[[i]] <- dummycapthist (pop = pop[i], fillvalue = fillvalue)
            class(output) <- c('capthist', 'list')
            session(output) <- 1:length(pop)
            output
        }
        else {
            
            newdim <- c(nrow(pop[[1]]), noccasions, nrow(traps))
            output <- array(fillvalue, dim = newdim)
            ## CAPTURE ON LAST OCCASION
            ## trick to keep array valid without misleading
            output[,noccasions,] <- 1   ###???
            class(output) <- 'capthist'
            traps(output) <- traps
            output
        }
    }
    ## --------------------------------------------------------------------
    ## process behavioural responses
    Markov <- any(c('B','Bk','K') %in% vars)
    btype <- which (c("b", "bk", "k") %in% tolower(vars))
    if (length(btype) > 1)
        stop ("all behavioural responses must be of same type in sim.detect")
    if (length(btype) == 0)
        btype <- 0
    ## --------------------------------------------------------------------
    ## setup
    if (is.null(details$ignoreusage)) {
        details$ignoreusage <- FALSE
    }
    N <- sapply(popnlist, nrow)
    S <- noccasions
    K <- nrow(traps)
    ncores <- setNumThreads()
    grain <- if (ncores==1) 0 else 1
    
    ##------------------------------------------
    ## detection parameters for each session, animal, occasion, detector
    ## realparval is lookup table of values,
    ## design$PIA is index to lookup
    
    ## real parameter values for naive animals
    
    ## using existing design to save time in secr.design.MS...
    if (length(unique(design0$PIA)) == 1) {
        dim0 <- dim(design0$PIA)
        design0$PIA <- array (1, dim = c(dim0[1], max(N), dim0[3:5]))
        dummyCH <- NULL
    }
    else {
        dummyCH <- dummycapthist(popnlist, fillvalue = 1)
        design0 <- secr.design.MS (dummyCH, model, timecov, NULL, NULL, NULL,
            NULL, ignoreusage = details$ignoreusage, naive = TRUE,
            CL = FALSE, contrasts = details$contrasts)
    }
    
    ## uniform over individuals
    ## C++ code uses individual-specific PIA even in full-likelihood case
    realparval0 <- makerealparameters (
        design0, 
        beta, 
        parindx, 
        link,
        fixed)
    
    ## real parameter  values for 'caughtbefore'  animals or detectors
    ## -- this  definition of  design1 differs  from that in secr.fit()
    ## where  parameter  values  are  appropriate to  the  particular
    ## (realised) detection histories
    
    if (btype > 0) {
        if (is.null(dummyCH))
            dummyCH <- dummycapthist(capthist, popnlist, fillvalue = 1)
        
        design1 <- secr.design.MS (dummyCH, model, timecov, NULL, NULL, NULL,
            NULL, ignoreusage = details$ignoreusage, naive = TRUE,
            CL = FALSE, contrasts = details$contrasts)
        realparval1 <- makerealparameters (design1, beta, parindx, link,
            fixed)  # caught before
    }
    else {   ## faster to just copy if no behavioural response
        design1 <- design0
        realparval1 <- realparval0
    }
    ##------------------------------------------
    
    # D <- getD (designD, beta, mask, parindx, link, fixed)
    dettype <- detectorcode(traps, MLonly = FALSE, noccasions = S)
    binomN <- expandbinomN(details$binomN, dettype)
    usge <- usage(traps)
    if (is.null(usge) | details$ignoreusage) {
        usge <- matrix(1, nrow = K, ncol = S)
    }
    output <- list()
    
    #------------------------------------------
    ## simulate this population...
    
    for (popnum in 1:length(N)) {
        
        NR <- N[popnum]
        
        if (all(detector(traps) %in% .localstuff$exclusivedetectors)) {
            maxdet <- NR * S
        }
        else {
            # safety margin : average detections per animal per detector per occasion
            maxdet <- NR * S * K * maxperpoly
        }
        
        if (all(dettype %in% c(-1,0,1,2,5,8))) {
            ## pre-compute distances from detectors to animals
            animals <- popnlist[[popnum]]
            distmat2 <- edist(as.matrix(traps), as.matrix(animals))^2
            
            ## precompute gk, hk for point detectors
            ## uses RcppParallel 
            
            gkhk0 <- makegkPointcpp (
                as.integer(detectfn),
                as.integer(grain),
                as.integer(ncores),
                as.matrix(realparval0),
                as.matrix(distmat2),
                as.double(0))
            
            gkhk <- makegkPointcpp (
                as.integer(detectfn),
                as.integer(grain),
                as.integer(ncores),
                as.matrix(realparval1),
                as.matrix(distmat2),
                as.double(0))
            
            if (any(dettype == 8)) {
                stop("sim.secr not ready for capped detectors")
                #     gkhk <- cappedgkhkcpp (as.integer(nrow(Xrealparval1)), as.integer(nrow(traps)),
                #                            as.double(getcellsize(mask)), as.double(density[,1]), 
                #                            as.double(gkhk$gk), as.double(gkhk$hk))  
            }
        }
        
        if (all(dettype %in% c(-1,0,1,2,8))) {
            temp <- simdetectpointcpp (dettype[1],      ## detector -1 single, 0 multi, 1 proximity, 2 count,... 
                as.integer(NR), 
                as.integer(nrow(realparval1)),
                as.double(gkhk0$gk), 
                as.double(gkhk$gk), 
                as.double(gkhk0$hk), 
                as.double(gkhk$hk), 
                as.integer(design0$PIA[1,1:NR,1:S,1:K,]),       ## naive
                as.integer(design1$PIA[1,1:NR,1:S,1:K,]),       
                as.matrix(usge),
                as.integer(btype),       ## code for behavioural response  0 none etc. 
                as.integer(Markov),      ## learned vs transient behavioural response 0 learned 1 Markov 
                as.integer(binomN)      ## number of trials for 'count' detector modelled with binomial 
            )
        }
        else {
            stop ("point detectors only")
        }
        if (temp$resultcode != 0) {
            stop ("simulated detection failed, code ", temp$resultcode)
        }
        w <- array(temp$value, dim=c(S, K, NR), dimnames = list(1:S,NULL,NULL))
        w <- aperm(w, c(3,1,2))
        w <- w[apply(w,1,sum)>0,,, drop = FALSE] 
        class(w)   <- 'capthist'    # NOT data.frame
        traps(w)   <- traps
        attr(w, 'detectedXY') <- NULL
        if (temp$n>0) {
            rownames(w) <- 1:temp$n
        }
        output[[popnum]] <- w
    }
    
    output
}
############################################################################################