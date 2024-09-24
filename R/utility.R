#######################################################################################
## ipsecr/R/utility.R
#######################################################################################

## 2022-04-04, 2022-05-04, 2022-05-09, 2022-12-26, 2024-01-15
#######################################################################################

# Global variables in namespace
#
## define a local environment for temporary variables e.g. iter

.localstuff <- new.env()

.localstuff$packageType <- ' pre-release'
##.localstuff$packageType <- ''

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
getdfn <- function (detectfn) {
    switch (detectfn+1, HN, HR, EX, CHN, UN, WEX, ANN, CLN, CG, BSS, SS, SSS,
                       SN, SNS, HHN, HHR, HEX, HAN, HCG, HVP, HPX)
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
        c('beta0','beta1','sdS'),  # 12  
        c('beta0','beta1','sdS'),  # 13 
        c('lambda0','sigma'),      # 14 hazard halfnormal
        c('lambda0','sigma','z'),  # 15 hazard hazard rate
        c('lambda0','sigma'),      # 16 hazard exponential
        c('lambda0','sigma','w'),  # 17
        c('lambda0','sigma','z'),  # 18
        c('lambda0','sigma','z'),  # 19
        c('lambda0','sigma'))      # 20 hazard pixelar 2021-03-25    

    # 2022-12-21 extended to allow user-defined parameters
    if (!is.null(details$extraparam)) {
        pnames <- c(pnames, names(details$extraparam))
    }
    if (details$param>0) warning ("details$param>0 is not recognised by ipsecr.fit")
    c('D', pnames)
}
#-------------------------------------------------------------------------------

memo <- function (text, verbose) {
    ## could use message(text), but does not immediately flush console
    if (verbose) { cat (text, '\n')
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

############################################################################################

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
}
SS <- function (r, pars, cutval) {
}
SSS <- function (r, pars, cutval) {
}
SN <- function (r, pars, cutval) {
}
SNS <- function (r, pars, cutval) {
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

################################################################################

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
################################################################################

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

## End of miscellaneous functions

################################################################################

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


ipsecr.lpredictor <- function (formula, newdata, indx, beta, field, beta.vcv=NULL,
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
        
        mat <- model.matrix(formula, data = newdata, contrasts = contrasts)
        rownames(mat) <- NULL
        if (nrow(mat) < nrow(newdata))
            warning ("missing values in predictors?")
        
        nmix <- 1

        ###############################
        Yp <- mat %*% beta[indx]
        ###############################
        
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

################################################################################

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
################################################################################

# number of beta parameters
nparameters <- function (object) {
    Npar <- max(unlist(object$parindx))
    ## allow for fixed beta parameters
    if (!is.null(object$details$fixedbeta))
        Npar <- Npar - sum(!is.na(object$details$fixedbeta))
    Npar
}
################################################################################

# modified from secr
getD <- function (parm = 'D', designD, beta, mask, parindx, link, fixed, nsessions) {
    if (is.function(designD)) {
        stop ("designD cannot be a function in ipsecr")
    }
    else {
        if ((is.null(designD) || nrow(designD)==0) && (is.null(fixed[[parm]]))) return(NULL)
    }
    if (nsessions>1)
        nmask <- max(sapply(mask, nrow))
    else
        nmask <- nrow(mask)
    D <- matrix(nrow = nmask, ncol = nsessions)
    
    if (!is.null(fixed[[parm]])) {
        D[] <- fixed[[parm]]
    }
    else {
        beta <- beta[parindx[[parm]]]
        D[] <- designD %*% beta
        D[] <- untransform (D, link[[parm]])
        # silently truncate D at zero
        D[D<0] <- 0
    }
    dimnames(D)[[2]] <- paste0(parm, 1:nsessions)
    D
}
##############################################################################

getDetDesignData <- function(popn, model, session, sessionlevels) {
    designdata <- popn   ## x,y
    if (!is.null(covariates(popn))) {
        designdata <- cbind(designdata, covariates(popn))
    }
    vars <- unlist(sapply(model, all.vars))
    if ('random' %in% vars) {
        designdata$random <- rnorm(nrow(popn))   # hold space
    }
    if ('session' %in% vars && !is.null(session)) {
        designdata$session <- rep(factor(session, levels = sessionlevels), nrow(popn))
    }
    found <- vars %in% names(designdata)
    if (sum(!found)>0) stop('detection predictors not found: ', vars[!found])
    designdata
}
################################################################################

getDetParMat <- function (popn, model, detectfn, beta, parindx, link, fixed, 
    details, sessionlevels, session = NULL) {
    # if ((is.null(design) || nrow(design)==0) && (is.null(fixed))) return(NULL)
    if (ms(popn)) {
        # for each session
        out <- mapply(getDetParMat, popn = popn, session = sessionlevels,
            MoreArgs = list(model, detectfn, beta, parindx, link, fixed, 
                details, sessionlevels), SIMPLIFY = FALSE)
        out
    }
    else {
        detectparnames <- secr::parnames(detectfn)
        npop <- nrow(popn)
        detparmat <- matrix(nrow = npop, ncol = length(detectparnames), 
            dimnames =list(NULL, detectparnames))
        designdata <- getDetDesignData(popn, model, session, sessionlevels)
        for (parm in detectparnames) {
            if (!is.null(fixed[[parm]])) {
                detparmat[,parm] <- fixed[[parm]]
            }
            else {
                if ('random' %in% all.vars(model[[parm]])) {
                    designdata$random <- rnorm(npop)
                }
                design <- model.matrix(model[[parm]], data = designdata, 
                    contrasts.arg = details$contrasts)                   
                detparmat[,parm] <- design %*% beta[parindx[[parm]]]
                detparmat[,parm] <- untransform (detparmat[,parm], link[[parm]])
            }
        }
        # silently truncate at zero
        detparmat[detparmat<0] <- 0
        detparmat
    }
}
##############################################################################

# return names of coefficients for each parameter
# requires prior simulation of popn, including any individual covariates
detBetaNames <- function(popn, model, detectfn, sessionlevels, fixed = NULL, 
    details = NULL) {
    if (ms(popn)) popn <- popn[[1]]
    detectparnames <- secr::parnames(detectfn)
    detparmat <- matrix(nrow = nrow(popn), ncol = length(detectparnames), 
        dimnames =list(NULL, detectparnames))
    designdata <- getDetDesignData(popn, model, sessionlevels[1], sessionlevels)
    nb <- function (parm) {
        out <- if (!is.null(fixed[[parm]])) character(0)
        else colnames(model.matrix(model[[parm]], data = designdata, 
            contrasts.arg = details$contrasts)) 
        out[out=='(Intercept)'] <- parm
        out
    }
    out <- lapply(detectparnames, nb)
    names(out) <- detectparnames
    out
    
}
################################################################################

# return names of extra parameters 2022-12-21
extraParNames <- function (details, fixed) {
    if (!is.null(details$extraparam)) {
        out <- names(details$extraparam)
        out <- out[!(out %in% names(fixed))]
        out
    }
    else character(0)
}
################################################################################

rpsv <- function (capthist)
{
    if (inherits (capthist, 'list')) {
        lapply(capthist, rpsv)   ## recursive
    }
    else {
        if (nrow(capthist) < 1) return(NA)
        trm <- as.matrix(traps(capthist))
        temp <- apply(abs(capthist), 1, rpsvcpp, trm)
        temp <- apply(temp,1,sum, na.rm = TRUE)
        if (any(is.na(temp) | temp<0)) {
            temp <- NA  
        }
        else {
            temp <- sqrt((temp[2]+temp[3]) / (2 * temp[1]))
        }
        temp
    }
}
rpsvi <- function (capthist)
{
    if (inherits (capthist, 'list')) {
        lapply(capthist, rpsvi)   ## recursive
    }
    else {
        onedxy <- function (dxy) {
            if (dxy[1]==0) NA else sqrt((dxy[2]+dxy[3]) / (2 * dxy[1]))
        }
        if (nrow(capthist) < 1) return(NA)
        trm <- as.matrix(traps(capthist))
        temp <- apply(abs(capthist), 1, rpsvcpp, trm)
        unname(apply(temp,2,onedxy))
    }
}

################################################################################
replacedefaults <- function (default, user) replace(default, names(user), user)
