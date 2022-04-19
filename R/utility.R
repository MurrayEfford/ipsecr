#######################################################################################
## ipsecr/R/utility.R
#######################################################################################

## 2022-04-04
#######################################################################################

# Global variables in namespace
#
## define a local environment for temporary variables e.g. iter
## e.g. Roger Peng https://stat.ethz.ch/pipermail/r-devel/2009-March/052883.html

.localstuff <- new.env()

.localstuff$packageType <- ' pre-release'
# .localstuff$packageType <- ''

.localstuff$validdetectors <- c('single','multi','proximity','count', 
    'polygonX', 'transectX', 'signal', 'polygon', 'transect', 
    'capped', 'null','null','null','null', 'telemetry', 'signalnoise')
.localstuff$simpledetectors <- c('single','multi','proximity','count', 'capped')
.localstuff$individualdetectors <- c('single','multi','proximity','count',
    'polygonX', 'transectX', 'signal', 'signalnoise', 'polygon', 'transect',
                                     'telemetry', 'capped')
.localstuff$pointdetectors <- c('single','multi','proximity','count',
    'signal', 'signalnoise', 'unmarked','presence','capped')
.localstuff$polydetectors <- c('polygon','transect','polygonX','transectX')
.localstuff$exclusivedetectors <- c('single','multi','polygonX','transectX')
.localstuff$countdetectors <- c('count','polygon','transect','unmarked','telemetry')
.localstuff$detectionfunctions <-
  c('halfnormal',
    'hazard rate',
    'exponential',
    'compound halfnormal',
    'uniform',
    'w exponential',
    'annular normal',
    'cumulative lognormal',
    'cumulative gamma',
    'binary signal strength',
    'signal strength',
    'signal strength spherical',
    'signal-noise',
    'signal-noise spherical',
    'hazard halfnormal',
    'hazard hazard rate',
    'hazard exponential',
    'hazard annular normal',
    'hazard cumulative gamma',
    'hazard variable power',
    'hazard pixelar')

.localstuff$DFN <- c('HN', 'HR', 'EX', 'CHN', 'UN', 'WEX', 'ANN', 'CLN', 'CG',
  'BSS', 'SS', 'SSS', 'SN', 'SNS',
  'HHN', 'HHR', 'HEX', 'HAN', 'HCG', 'HVP','HPX')

.localstuff$learnedresponses <- c('b', 'bk', 'B', 'k', 'Bk')   ## Bk added 2020-02-26

detectionfunctionname <- function (fn) {
    .localstuff$detectionfunctions[fn+1]
}

detectionfunctionnumber <- function (detname) {
    dfn <- match (toupper(detname), .localstuff$DFN)
    if (is.na(dfn))
        dfn <- match (tolower(detname), .localstuff$detectionfunctions)
    if (is.na(dfn))
        stop ("unrecognised detection function ", detname)
    dfn-1
}
parnames <- function (detectfn) {
    parnames <- switch (detectfn+1,
        c('g0','sigma'),   ## 0
        c('g0','sigma','z'),
        c('g0','sigma'),
        c('g0','sigma','z'),
        c('g0','sigma'),
        c('g0','sigma','w'),
        c('g0','sigma','w'),
        c('g0','sigma','z'),
        c('g0','sigma','z'),
        c('b0','b1'),
        c('beta0','beta1', 'sdS'),    ## include cutval?
        c('beta0','beta1', 'sdS'),    ## include cutval?
        c('beta0','beta1', 'sdS','muN','sdN'),
        c('beta0','beta1', 'sdS','muN','sdN'),
        c('lambda0','sigma'),
        c('lambda0','sigma','z'),
        c('lambda0','sigma'),
        c('lambda0','sigma','w'),
        c('lambda0','sigma','z'),
        c('lambda0','sigma','z'),
        c('lambda0')    ## 20
    )
}
getdfn <- function (detectfn) {
    switch (detectfn+1, HN, HR, EX, CHN, UN, WEX, ANN, CLN, CG, BSS, SS, SSS,
                       SN, SNS, HHN, HHR, HEX, HAN, HCG, HVP, HPX)
}

valid.detectfn <- function (detectfn, valid = c(0:3,5:19, 20)) {
# exclude 4 uniform: too numerically flakey
    if (is.null(detectfn))
        stop ("requires 'detectfn'")
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if (!(detectfn %in% valid))
        stop ("invalid detection function")
    detectfn
}

valid.detectpar <- function (detectpar, detectfn) {
    if (is.null(detectpar) | is.null(detectfn))
        stop ("requires valid 'detectpar' and 'detectfn'")
    if (!all(parnames(detectfn) %in% names(detectpar)))
        stop ("requires 'detectpar' ", paste(parnames(detectfn), collapse=','),
            " for ", detectionfunctionname(detectfn), " detectfn")
    detectpar[parnames(detectfn)]
}

valid.pnames <- function (details, CL, detectfn, alltelem, sighting, nmix) {
    ## modelled parameters
    pnames <- switch (detectfn+1,
        c('g0','sigma'),           # 0 halfnormal
        c('g0','sigma','z'),       # 1 hazard rate
        c('g0','sigma'),           # 2 exponential
        c('g0','sigma','z'),       # 3
        c('g0','sigma'),           # 4
        c('g0','sigma','w'),       # 5
        c('g0','sigma','w'),       # 6
        c('g0','sigma','z'),       # 7
        c('g0','sigma','z'),       # 8
        c('b0','b1'),              # 9
        c('beta0','beta1','sdS'),  # 10
        c('beta0','beta1','sdS'),  # 11
        c('beta0','beta1','sdS'),  # 12  cf parnames() in utility.R: muN, sdN?
        c('beta0','beta1','sdS'),  # 13  cf parnames() in utility.R: muN, sdN?
        c('lambda0','sigma'),      # 14 hazard halfnormal
        c('lambda0','sigma','z'),  # 15 hazard hazard rate
        c('lambda0','sigma'),      # 16 hazard exponential
        c('lambda0','sigma','w'),  # 17
        c('lambda0','sigma','z'),  # 18
        c('lambda0','sigma','z'),  # 19
        c('lambda0','sigma'))      # 20 hazard pixelar 2021-03-25    

    if (details$param %in% c(2,6))
        pnames[1] <- 'esa'
    if (details$param %in% c(3,5))
        pnames[1] <- 'a0'
    if (details$param %in% 4:6) {
        pnames[2] <- 'sigmak'
        pnames <- c(pnames, 'c')
        pnames <- c(pnames, 'd')
    }
    c('D', pnames)
}
#-------------------------------------------------------------------------------

expanddet <- function(CH) {
    trps <- traps(CH)
    if (is.null(trps))
        return ('nonspatial')
    else {
        det <- detector(trps)
        if (length(det)<ncol(CH))
            rep(det[1], ncol(CH))
        else det
    }
}

ndetectpar <- function (detectfn) {
    length(parnames(detectfn))
}

replacedefaults <- function (default, user) replace(default, names(user), user)

discreteN <- function (n, N) {
    tN <- trunc(N)
    if (N != tN) tN + sample (x = c(1,0), prob = c(N-tN, 1-(N-tN)),
        replace = T, size = n)
    else rep(tN,n)
}

ndetector <- function (traps) {
    if (is.null(traps))
        return(1)
    else if (all(detector(traps) %in% .localstuff$polydetectors))
        length(levels(polyID(traps)))
    else
        nrow(traps)
}

memo <- function (text, trace) {
    ## could use message(text), but does not immediately flush console
    if (trace) { cat (text, '\n')
    flush.console() }
}

## regularize a list of formulae
stdform <- function (flist) {
    LHS <- function (form) {
        trms <- as.character (form)
        if (length(trms)==2) '' else trms[2]
    }
    RHS <- function (form) {
        trms <- as.character (form)
        ## 2020-05-14 for compatibility with R 4.0
        if (length(trms)==3) as.formula(paste(trms[c(1,3)], collapse = " ")) else form
    }
    lhs <- sapply(flist, LHS)
    temp <- lapply(flist, RHS)
    if (is.null(names(flist))) names(temp) <- lhs
    else names(temp) <- ifelse(names(flist) == '', lhs, names(flist))
    temp
}

## Start of miscellaneous functions

invlogit <- function (y) 1/(1+exp(-y))   # plogis(y)
logit    <- function (x) log(x/(1-x))    # qlogis(x), except for invalid argument
sine     <- function (x) asin (x*2-1)
invsine  <- function (y) (sin(y)+1) / 2
odds     <- function (x) x / (1-x)
invodds  <- function (y) y / (1+y)

lnbinomial <- function (x,size,prob) {
  lgamma (size+1) - lgamma (size-x+1) - lgamma (x+1) +
      x * log(prob) + (size-x) * log (1-prob)
}
############################################################################################
## moved from methods.r 2012-10-28

model.string <- function (model, userDfn) {
    if (!is.null(userDfn)) {
        if (!is.null(model$D))
            model$D <- paste('~userD', userDfn('name'), sep='.')
    }
    temp <- paste (names(model), as.character(model), collapse=' ', sep='')
    temp
}
fixed.string <- function (fixed) {
    if (is.null(fixed) | length(fixed)==0) 'none'
    else paste (names(fixed), as.character(fixed), collapse=', ', sep=' = ')
}
############################################################################################

var.in.model <- function(v,m) v %in% unlist(lapply(m, all.vars))

############################################################################################

add.cl <- function (df, alpha, loginterval, lowerbound = 0) {

## add lognormal or standard Wald intervals to dataframe with columns
## 'estimate' and 'SE.estimate'
## lowerbound added 2011-07-15
    z <- abs(qnorm(1-alpha/2))
    if (loginterval) {
        delta <- df$estimate - lowerbound
        df$lcl <- delta / exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
        df$ucl <- delta * exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
    }
    else {
        df$lcl <- pmax(lowerbound, df$estimate - z * df$SE.estimate)
        df$ucl <- df$estimate + z * df$SE.estimate
    }
    df
}

###############################################################################

spatialscale <- function (object, detectfn, session = '') {
    if (inherits(object, 'secr')) {
        if (ms(object))
            detpar <- detectpar(object)[[session]]
        else
            detpar <- detectpar(object)
        cutval <- object$details$cutval
    }
    else {
        detpar <- object
        cutval <- object$cutval
    }
    if (!is.null(detpar$sigma)) detpar$sigma
    else if (detectfn == 10) {
        (cutval - detpar$beta0) / detpar$beta1
    }
    else if (detectfn == 11) {
        d11 <- function(d, beta0, beta1, c) beta0 +
            beta1 * (d-1) - 10 * log10(d^2) - c
        interval <- c(0,10 * (cutval - detpar$beta0) / detpar$beta1)
        uniroot (d11, interval, detpar$beta0, detpar$beta1, cutval)$root
    }
    else if (detectfn == 9) {
#        (0.5 - detpar$b0) / detpar$b1
        - 1 / detpar$b1   ## 2010-11-01
    }
    else stop ("unrecognised detectfn")
}

###############################################################################
## mean and SD if x numeric
## was statfn 2011-11-08
getMeanSD <- function(xy) {
    MeanSD <- function (x) {
        if (is.numeric(x))
            c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
        else
            c(NA,NA)
    }
    as.data.frame (apply(xy, 2, MeanSD))
}
###############################################################################

## clunky but effective re-write 2012-09-04, improved 2016-02-20, 2016-05-10
leadingzero <- function (x) {
    xc <- as.character(x)
    w <- max(nchar(xc))
    n0 <- function(n) paste(rep('0',n), collapse='')
    paste(sapply(w-nchar(xc), n0), x, sep='')

    ## or, 2016-01-15, 2016-02-20 BUT DOESN'T HANDLE NON-INTEGER 2016-05-10
    #     if (is.character(x)) x <- as.numeric(x)
    #     sprintf(paste("%0", w, "d", sep = ""), x)
}

###############################################################################

## Detection functions

HN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]
    g0 * exp (-r^2 / 2 / sigma^2)
}
HR <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    g0 * (1 - exp (-(r / sigma)^-z))
}
EX <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]
    g0 * exp (-r / sigma)
}
UN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]
    ifelse (r<=sigma, g0, 0)
}
CHN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    g0 * ( 1 - (1 - exp (-r^2 / 2 / sigma^2)) ^ z )
}
WEX <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
    ifelse(r<=w, g0, g0*exp (-(r-w) / sigma))
}
ANN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
    g0 * exp (-(r-w)^2 / 2 / sigma^2)
}
CLN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    CV2 <- (z/sigma)^2
    sdlog <- log(1 + CV2)^0.5
    meanlog <- log(sigma) - sdlog^2/2
    g0 * plnorm(r, meanlog, sdlog, lower.tail = FALSE)
}
CG <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    g0 * pgamma(r, shape=z, scale=sigma/z, lower.tail = FALSE)
}
CN <- function (r, pars, cutval) {
    g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    x <- z * (r - sigma)
    g0 * (1 + (1 - exp(x)) / (1 + exp(x)))/2
}
BSS <- function (r, pars, cutval) {
    b0 <- pars[1]; b1 <- pars[2]
    gam <- -(b0 + b1 * r);
    pnorm (gam, mean=0, sd=1, lower.tail=FALSE)
}
SS <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
    if (is.null(cutval))
        stop ("require 'details$cutval' for signal strength plot")
    mu <- beta0 + beta1 * r
    1 - pnorm (q=cutval, mean=mu, sd=sdS)
}
SSS <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
    if (is.null(cutval))
        stop ("require 'details$cutval' for signal strength plot")
    ## spherical so assume distance r measured from 1 m
    mu <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
    mu[r<1] <- beta0
    1 - pnorm (q=cutval, mean=mu, sd=sdS)
}
SN <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3];
    muN <- pars[4]; sdN <- pars[5]
    muS <- beta0 + beta1 * r
    1 - pnorm (q=cutval, mean=muS-muN, sd=sqrt(sdS^2+sdN^2))
}
SNS <- function (r, pars, cutval) {
    beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3];
    muN <- pars[4]; sdN <- pars[5]
    ## spherical so assume distance r measured from 1 m
    muS <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
    muS[r<1] <- beta0
    1 - pnorm (q=cutval, mean=muS-muN, sd=sqrt(sdS^2+sdN^2))
}
HHN <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]
    1 - exp(-lambda0 * exp (-r^2 / 2 / sigma^2))
}
HHR <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    1 - exp(-lambda0 * ( 1 - exp (-(r / sigma)^-z)))
}
HEX <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]
    1 - exp(-lambda0 * exp (-r / sigma))
}
HAN <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
    1 - exp(-lambda0 * exp (-(r-w)^2 / 2 / sigma^2))
}
HCG <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    lambda0 * pgamma(r, shape=z, scale=sigma/z, lower.tail = FALSE)
}
HVP <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    1 - exp(-lambda0 * exp(-(r/sigma)^z))
}
HPX <- function (r, pars, cutval) {
    g0 <- 1-exp(-pars[1])
    radius <- pars[2]
    ifelse (r<=radius, g0, 0)  # circular, not square! crude approx
}

############################################################################################

gradient <- function (pars, fun, eps=0.001, ...)
## quick & dirty 2009 09 14
## used by plot.secr for delta method limits
{
  est <- pars
  g   <- pars
  for (i in 1:length(est))
  {
      temp     <- est[i]
      if (temp != 0.0) delta <- eps * abs(temp)
      else             delta <- eps
      est[i]  <- temp - delta
      fminus  <- fun (est, ...)
      est[i]  <- temp + delta
      fplus   <- fun (est, ...)
      g[i]    <- (fplus - fminus) / (2.0 * delta)
      est[i]  <- temp;
  }
  g
}
############################################################################################

#-------------------------------------------------------------------------------

# transformation tidy up 2021-12-16
# arbitrary link function specified with functions X, invX, se.invX

transform <- function (x, link) {
    switch (link,
        identity = x,
        i1000 = x * 1000,
        log = log(x),
        neglog = log(-x),
        logit = logit(x),
        odds = odds(x),
        sin = sine(x),
        do.call(link, list(x))
    )
}
#-------------------------------------------------------------------------------

# used only in model.average, modelAverage
se.transform <- function (real, sereal, link) {
    switch (link,
        identity = sereal,
        i1000 = sereal / 1000,
        log = log((sereal/real)^2 + 1)^0.5,
        neglog = log((sereal/-real)^2 + 1)^0.5,
        logit = sereal / real / (1 - real),
        sin = NA,
        do.call(paste0('se.',link), list(real, sereal) )
    )
}
#-------------------------------------------------------------------------------

untransform <- function (beta, link) {
    switch (link,
        identity = beta,
        i1000 = beta / 1000,
        log = exp(beta),
        neglog = -exp(beta),
        logit = invlogit(beta),
        odds = invodds(beta),
        sin = invsine(beta),
        do.call(paste0('inv',link), list(beta))
    )
}
#-------------------------------------------------------------------------------

se.untransform <- function (beta, sebeta, link) {
    # Approximate translation of SE to untransformed scale
    # Delta method cf Lebreton et al 1992 p 77
    switch (link,
        identity = sebeta,
        i1000 = sebeta / 1000,
        log = exp(beta) * sqrt(exp(sebeta^2)-1),
        neglog = exp(beta) * sqrt(exp(sebeta^2)-1),
        logit = invlogit(beta) * (1-invlogit(beta)) * sebeta,
        sin = NA,                ####!!!!
        do.call(paste0('se.inv', link), list(beta=beta, sebeta=sebeta))
    )
}
#-------------------------------------------------------------------------------

# vectorized transformations
Xtransform <- function (real, linkfn, varnames) {
    mapply(transform, real, linkfn[varnames])
}
se.Xtransform <- function (real, sereal, linkfn, varnames) {
    mapply(se.transform, real, sereal, linkfn[varnames])
}
Xuntransform <- function (beta, linkfn, varnames) {
    mapply(untransform, beta, linkfn[varnames])
}
se.Xuntransform <- function (beta, sebeta, linkfn, varnames)
{
    if (length(beta)!=length(sebeta))
        stop ("'beta' and 'sebeta' do not match")
    if (!all(varnames %in% names(linkfn)))
        stop ("'linkfn' component missing for at least one real variable")
    mapply(se.untransform, beta, sebeta, linkfn[varnames])
}
#-------------------------------------------------------------------------------

mlogit.untransform <- function (beta, latentmodel) {
    if (!missing(latentmodel)) {
        for (i in unique(latentmodel))
            beta[latentmodel==i] <- mlogit.untransform(beta[latentmodel==i])
        beta
    }
    else {
        ## beta should include values for all classes (mixture components)
        nmix <- length(beta)
        if (sum(is.na(beta)) != 1) {
            ## require NA for a single reference class
            rep(NA, length(beta))
        }
        else {
            nonreference <- !is.na(beta)   # not reference class
            b <- beta[nonreference]
            pmix <- numeric(nmix)
            pmix[nonreference] <- exp(b) / (1+sum(exp(b)))
            pmix[!nonreference] <- 1 - sum(pmix[nonreference])
            pmix
        }
    }
}

clean.mlogit <- function(x) {
    ## 2014-08-19 for robustness...
    if (is.na(x[2])) x[2] <- 1-x[1]
    x[1] <- NA   ## assumed reference class
    logit(mlogit.untransform(x))
}

mlogit <- function (x) {
    ## return the mlogit of an unscaled vector of positive values
    ## 2013-04-14
    logit(x/sum(x))
}

## End of miscellaneous functions

############################################################################################

###############################################################################

## expand beta parameter vector using template of 'fixed beta'
## fixed beta fb input is missing (NA) for estimated beta parameters
fullbeta <- function (beta, fb) {
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta  ## partial beta (varying only)
        beta <- fb             ## complete beta
    }
    beta
}
###############################################################################

## moved from pdot.R 2013-11-09
## scalar 2016-10-14
getbinomN <- function (binomN, detectr) {
    if (any(detectr %in% .localstuff$countdetectors)) {
        if (is.null(binomN))
            return(0)
        else if (binomN == 'usage')
            return(1)
        else
            return(binomN)
    }
    else
        return(1)
}
###############################################################################

complete.beta <- function (object) {
    fb <- object$details$fixedbeta
    if (!is.null(fb)) {
        nbeta <- length(fb)
        fb[is.na(fb)] <- object$beta
        beta <- fb
    }
    else {
        beta <- object$beta
    }
    beta
}
###############################################################################

complete.beta.vcv <- function (object) {
    fb <- object$details$fixedbeta
    if (!is.null(fb)) {
        nbeta <- length(fb)
        beta.vcv <- matrix(NA, nrow = nbeta, ncol = nbeta)
        beta.vcv[is.na(fb[row(beta.vcv)]) & is.na(fb[col(beta.vcv)])] <- object$beta.vcv
    }
    else {
        beta.vcv <- object$beta.vcv
    }
    beta.vcv
}
###############################################################################

general.model.matrix <- function (formula, data, gamsmth = NULL, 
    contrasts = NULL, ...) {

    ## A function to compute the design matrix for the model in
    ## 'formula' given the data in 'data'. This is merely the result
    ## of model.matrix() unless 'formula' includes smooth terms -- s()
    ## or te() as described in mgcv ?formula.gam.

    ## smooth terms blocked in ipsecr
    dots <- list(...)

    ## model.matrix(formula, data, ...)
    mat <- model.matrix(formula, data = data, contrasts.arg = contrasts)
    rownames (mat) <- NULL
    mat
}
###############################################################################

## shifted from secrloglik 2016-10-16
makerealparameters <- function (design, beta, parindx, link, fixed) {
    modelfn <- function(i) {
        ## linear predictor for real parameter i
        Yp <- design$designMatrices[[i]] %*% beta[parindx[[i]]]
        if (names(link)[i] == 'pmix') {
            ## 2013-04-14 index of class groups (pmix sum to 1.0 within latentmodel)
            cols <- dimnames(design$designMatrices[[i]])[[2]]
            h2 <- grep('.h2', cols, fixed=T)
            h3 <- grep('.h3', cols, fixed=T)
            h2c <- grep(':h2', cols, fixed=T)
            h3c <- grep(':h3', cols, fixed=T)
            h.cols <- c(h2,h3,h2c,h3c)
            tmp <- design$designMatrices[[i]][,-h.cols, drop = FALSE]
            tmph <- design$designMatrices[[i]][,h.cols, drop = FALSE]
            ## 2018-02-23 why as.numeric()? 
            latentmodel <- as.numeric(factor(apply(tmp,1,paste, collapse='')))
            refclass <- apply(tmph,1,sum) == 0
            Yp[refclass] <- NA
            Yp <- mlogit.untransform(Yp, latentmodel)
            Yp[design$parameterTable[,i]]
        }
        else {
            Yp <- untransform(Yp, link[[i]])
            Yp[design$parameterTable[,i]]   ## replicate as required
        }
    }
    ## construct matrix of detection parameters
    nrealpar  <- length(design$designMatrices)
    parindx$D <- NULL ## detection parameters only
    link$D    <- NULL ## detection parameters only
    parindx$noneuc <- NULL ## detection parameters only
    link$noneuc    <- NULL ## detection parameters only
    detectionparameters <- names(link)
    fixed.dp <- fixed[detectionparameters[detectionparameters %in% names(fixed)]]
    
    if (length(fixed.dp)>0)
        for (a in names(fixed.dp))  ## bug fixed by adding this line 2011-09-28
            link[[a]] <- NULL
    if (length(link) != nrealpar)
        stop ("number of links does not match design matrices")
    
    if (nrealpar == 0) {
        return(matrix(unlist(fixed.dp),nrow = 1))
    }
    
    temp <- sapply (1:nrealpar, modelfn)
    if (nrow(design$parameterTable)==1) temp <- t(temp)
    nrw <- nrow(temp)
    ## make new matrix and insert columns in right place
    temp2 <- as.data.frame(matrix(nrow = nrw, ncol = length(detectionparameters)))
    names(temp2) <- detectionparameters
    temp2[ , names(design$designMatrices)] <- temp          ## modelled
    if (!is.null(fixed.dp) & length(fixed.dp)>0)
        temp2[ , names(fixed.dp)] <- sapply(fixed.dp, rep, nrw)    ## fixed
    as.matrix(temp2)
    
}
############################################################################################

secr.lpredictor <- function (formula, newdata, indx, beta, field, beta.vcv=NULL,
    smoothsetup = NULL, contrasts = NULL, f = NULL) {
    ## form linear predictor for a single 'real' parameter
    ## smoothsetup should be provided whenever newdata differs from
    ## data used to fit model and the model includes smooths from gam
    vars <- all.vars(formula)
    OK <- vars %in% names(newdata)
    if (any(!OK)) {
        missingvars <- paste(vars[!OK], collapse = ', ')
        if (sum(!OK) == 1)
            stop ("model covariate ", missingvars, " not found in 'newdata'")
        else
            stop ("model covariates ", missingvars, " not found in 'newdata'")
    }
    newdata <- as.data.frame(newdata)
    lpred <- matrix(ncol = 2, nrow = nrow(newdata), dimnames = list(NULL,c('estimate','se')))

    if (!is.null(f) && field == 'D') {
       Yp <- f(newdata[,vars[1]], beta = beta[indx]) 
       mat <- as.matrix(newdata[,vars[1], drop = FALSE])
    }
    else {
        
        mat <- general.model.matrix(formula, data = newdata, gamsmth = NULL, 
            contrasts = contrasts)
        if (nrow(mat) < nrow(newdata))
            warning ("missing values in predictors?")
        
        nmix <- 1
        if (field=='pmix') {
            ## drop pmix beta0 column from design matrix (always zero)
            mat <- mat[,-1,drop=FALSE]
            if ('h2' %in% names(newdata)) nmix <- 2
            if ('h3' %in% names(newdata)) nmix <- 3
            mixfield <- c('h2','h3')[nmix-1]
        }
        
        ###############################
        Yp <- mat %*% beta[indx]
        ###############################
        
        ## A latent model comprises one row for each latent class.
        ## Back transformation of pmix in mlogit.untransform() requires all rows of 
        ## each latent model. That function splits vector Yp by latent model.
        
        if (field == 'pmix') {
            nonh <- newdata[, names(newdata) != mixfield, drop = FALSE]
            latentmodel <- factor(apply(nonh, 1, paste, collapse = ''))
            refclass <- as.numeric(newdata[, mixfield]) == 1
            Yp[refclass] <- NA
            Yp <- mlogit.untransform(Yp, latentmodel)
            Yp <- logit(Yp)  # return to logit scale for later untransform!
            if (nmix==2) {
                h2.1 <- as.numeric(newdata$h2)==1
                h2.2 <- as.numeric(newdata$h2)==2
            }
        }
    }

    lpred[,1] <- Yp
    if (is.null(beta.vcv) || (any(is.na(beta[indx])))) return ( cbind(newdata,lpred) )
    else {
        if (is.null(f) || field != 'D') {
            vcv <- beta.vcv[indx,indx, drop = FALSE]
            vcv[is.na(vcv)] <- 0
            nrw <- nrow(mat)
            vcv <- apply(expand.grid(1:nrw, 1:nrw), 1, function(ij)
                mat[ij[1],, drop=F] %*% vcv %*% t(mat[ij[2],, drop=F])) 
            
            vcv <- matrix (vcv, nrow = nrw)
            if (field=='pmix') {
                if (nmix==2)
                    vcv[h2.1,h2.1] <- vcv[h2.2,h2.2]
                else
                    vcv[,] <- NA
            }
            lpred[,2] <- diag(vcv)^0.5
        }
        else {
            vcv <- NULL
        }
        
        temp <- cbind(newdata,lpred)
        attr(temp, 'vcv') <- vcv
        return(temp)
    }
}

############################################################################################

## 2014-10-25, 2017-01-24
## intercept and fix certain ,models with bad defaults
updatemodel <- function (model, detectfn, detectfns, oldvar, newvar, warn = FALSE) {
    if (detectfn %in% detectfns) {
        for (i in 1:length(oldvar)) {
            if (oldvar[i] %in% names(model)) {
                names(model)[names(model) == oldvar[i]] <- newvar[i]
                if (warn)
                    warning ("replacing ", oldvar[i], " by ", newvar[i],
                             " in model for detectfn ", detectfn)
            }
        }
    }
    model
}
############################################################################################

## Manually remove some mask points
# simplified 2022-02-03

deleteMaskPoints <- function (mask, onebyone = TRUE, add = FALSE, poly = NULL,
                              poly.habitat = FALSE, ...) {
    ## interface does not work properly in RStudio

    if (ms(mask)) {         ## a list of mask objects
        if (inherits(poly, 'list') & (!is.data.frame(poly)))
            stop ("lists of polygons not implemented in 'make.mask'")
        temp <- lapply (mask, deleteMaskPoints, onebyone = onebyone, add = add,
                        poly = poly, poly.habitat = poly.habitat, ...)
        class (temp) <- c('mask', 'list')
        temp
    }
    else {
        plot(mask, add = add, ...)
        if (!is.null(poly)) {
            if (poly.habitat)
                pointstodrop <- (1:nrow(mask))[!pointsInPolygon(mask, poly)]
            else
                pointstodrop <- (1:nrow(mask))[pointsInPolygon(mask, poly)]
        }
        else if (onebyone) {
            cat ('Click to select points; right-click to stop\n')
            flush.console()
            xy <- locator(type = 'p', pch=1, col='red')
            pointstodrop <- if (length(xy$x)==0)
                numeric(0)
            else
                nearesttrap(xy, mask)
        }
        else {
            cat ('Click to select polygon vertices; right-click to stop\n')
            flush.console()
            xy <- locator(type = 'l', col='red')
            xy <- as.data.frame(xy)
            xy <- rbind(xy, xy[1,])
            if (poly.habitat)
                pointstodrop <- (1:nrow(mask))[!pointsInPolygon(mask, xy)]
            else
                pointstodrop <- (1:nrow(mask))[pointsInPolygon(mask, xy)]
        }
        npts <- length(pointstodrop)
        if (npts>0) {
            points(mask[pointstodrop,], pch = 16, col = 'red')
            if(.Platform$OS.type == "windows") {
                pl <- if (npts>1) 's' else ''
                msg <- paste ('Delete ', npts, ' red point',pl, '?', sep='')
                response <-  utils::winDialog(type = "okcancel", msg)
            } else {
                response <- 'OK'
            }
            if (response == 'OK') {
                mask <- subset(mask, -pointstodrop)
            if (npts==1)
                message("1 point deleted")
            else
                message(npts, " points deleted")
            }
        else
            message ("point(s) not deleted")
        }
        else
            message ("no points to delete")
        plot(mask, col='green')
        mask
    }
}
############################################################################################

# number of beta parameters
nparameters <- function (object) {
    Npar <- max(unlist(object$parindx))
    ## allow for fixed beta parameters
    if (!is.null(object$details$fixedbeta))
        Npar <- Npar - sum(!is.na(object$details$fixedbeta))
    Npar
}
############################################################################################

mapbeta <- function (parindx0, parindx1, beta0, betaindex)

    ## Extend beta vector from simple model (beta0) to a more complex (i.e. general)
    ## model, inserting neutral values (zero) as required.
    ## For each real parameter, a 1:1 match is assumed between
    ## beta values until all beta values from the simpler model are
    ## used up. THIS ASSUMPTION MAY NOT BE JUSTIFIED.
    ## betaindex is a user-controlled alternative.

{
    ## list of zeroed vectors, one per real parameter
    beta1 <- lapply(parindx1, function (x) {x[]<-0; x})

    if (!is.null(betaindex)) {
        beta1 <- unlist(beta1)
        if (sum(betaindex>0) != length(beta0))
            stop ("invalid 'betaindex'")
        beta1[betaindex] <- beta0
        beta1
    }
    else {
        ## indx is within-parameter rather than absolute index
        ## for each _original_ real parameter
        indx <- lapply(parindx0, function(x) x-x[1]+1)
        ## for (j in 1:length(beta1))
        ## improved replace by name2015-11-17
        for (j in names(beta1)) {
            if (j %in% names(beta0))
                beta1[[j]][indx[[j]]] <- beta0[parindx0[[j]]]
        }
        unlist(beta1)
    }
}
############################################################################################

detectorcode <- function (object, MLonly = TRUE, noccasions = NULL) {
  ## numeric detector code from a traps object
  detcode <- sapply(detector(object), switch,
    single      = -1,
    multi       = 0,
    proximity   = 1,
    count       = 2,
    polygonX    = 3,
    transectX   = 4,
    signal      = 5,
    polygon     = 6,
    transect    = 7,
    capped      = 8,
    unmarked    = 10,
    presence    = 11,
    signalnoise = 12,
    telemetry   = 13,
    -2)
  
  if (MLonly) {
    detcode <- ifelse (detcode==-1, rep(0,length(detcode)), detcode)
    if (any(detcode<0))
      stop ("Unrecognised detector type")
  }
  
  if (!is.null(noccasions) & (length(detcode)==1))
    detcode <- rep(detcode, noccasions)
  detcode
}
############################################################################################

expandbinomN <- function (binomN, detectorcodes) {
    # assumes detectorcodes is a vector of length = noccasions
    binomN <- ifelse (detectorcodes %in% c(2,6,7), binomN, 1)
    if (any(is.na(binomN))) stop ("NA value in binomN")
    binomN
}
############################################################################################


newstr <-function (strings) {
    ## compress a character vector
    ## use run length encoding function
    rl <- rle(strings)
    st <- rl$values
    le <- paste0(' (',as.character(rl$lengths), ')')
    le[le==' (1)'] <- ''
    paste(paste0(st, le), collapse = ', ')
}
# newstr(c("single", rep("proximity",4)))
############################################################################################

## Based on Charles C. Berry on R-help 2008-01-13
## drop = FALSE 2018-11-22

## used in secr.design.MS

n.unique.rows <- function(x) {
  order.x <- do.call(order, as.data.frame(x))
  equal.to.previous <- rowSums(x[tail(order.x,-1),,drop = FALSE] != 
                                 x[head(order.x,-1),,drop = FALSE])==0 
  1 + sum(!equal.to.previous)
}

individualcovariates <- function (PIA) {
  pia <- matrix(aperm(PIA, c(2:5,1)), nrow = dim(PIA)[2])
  n.unique.rows(pia) > 1
}
##############################################################################

# modified from secr
getD <- function (designD, beta, mask, parindx, link, fixed) {
  if (!is.function(designD)) {
    if ((is.null(designD) || nrow(designD)==0) && (is.null(fixed[['D']]))) return(NULL)
  }
  nmask <- nrow(mask)
  D <- numeric(nmask)
  if (!is.null(fixed[['D']])) {
    D[] <- fixed[['D']]
  }
  else {
    beta <- beta[parindx[['D']]]
    D[] <- designD %*% beta
    D[] <- untransform (D, link[['D']])
    # silently truncate D at zero
    D[D<0] <- 0
  }
  D
}
##############################################################################

# copied from secr
mapbeta <- function (parindx0, parindx1, beta0, betaindex)
  
  ## Extend beta vector from simple model (beta0) to a more complex (i.e. general)
  ## model, inserting neutral values (zero) as required.
  ## For each real parameter, a 1:1 match is assumed between
  ## beta values until all beta values from the simpler model are
  ## used up. THIS ASSUMPTION MAY NOT BE JUSTIFIED.
  ## betaindex is a user-controlled alternative.
  
{
  ## list of zeroed vectors, one per real parameter
  beta1 <- lapply(parindx1, function (x) {x[]<-0; x})
  
  if (!is.null(betaindex)) {
    beta1 <- unlist(beta1)
    if (sum(betaindex>0) != length(beta0))
      stop ("invalid 'betaindex'")
    beta1[betaindex] <- beta0
    beta1
  }
  else {
    ## indx is within-parameter rather than absolute index
    ## for each _original_ real parameter
    indx <- lapply(parindx0, function(x) x-x[1]+1)
    ## for (j in 1:length(beta1))
    ## improved replace by name2015-11-17
    for (j in names(beta1)) {
      if (j %in% names(beta0))
        beta1[[j]][indx[[j]]] <- beta0[parindx0[[j]]]
    }
    unlist(beta1)
  }
}
############################################################################################

