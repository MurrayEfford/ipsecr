###############################################################################
## package 'ipsecr'
## proxyfn.R
## 2022-05-08
## 2022-06-11 proxy.ms
## 2022-06-15 proxy.ms extended for count detectors
## 2022-08-24 fixed major bug in zippin
## 2022-08-24 defunct proxyfn0
## 2022-08-27 proxy.ms glm log(si+1)

###############################################################################

proxyfn0 <- function (capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...) {
    .Defunct('proxyfn1', package = 'ipsecr') ## 2022-08-24
}

##################################################

proxyfn1 <- function (capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...) {
    N.estimator <- tolower(N.estimator)
    N.estimator <- match.arg(N.estimator)
    ## capthist single-session only; ignoring losses
    n <- nrow(capthist)         ## number of individuals
    ch <- apply(abs(capthist), 1:2, sum) > 0
    
    nocc <- ncol(capthist)      ## number of occasions
    ni <- apply(ch, 2, sum)     ## individuals on each occasion
    estimates <- {
        if (N.estimator == "n")
            c(n, sum(ni)/n/nocc)
        else if (N.estimator == "null")
            M0(c(sum(ni), n, nocc))
        else if (N.estimator == "zippin") {
            tempx2 <- apply(ch, 1, function(x) cumsum(abs(x))>0)
            Mt1 <- apply(tempx2,1,sum)
            ui <- c(ni[1], diff(Mt1))
            Mb(ni, ui)
        }
        else if (N.estimator == "jackknife") {
            fi <- tabulate(apply(ch,1,sum), nbins = nocc)
            Mh(fi)
        }
    }
    c(
        logN = log(estimates[1]), 
        cloglogp = log(-log(1-estimates[2])), 
        logrpsv= log(rpsv(capthist))
    )
}
##################################################

detectionDesignData <- function (capthist, byoccasion = FALSE, ...) {
    # start with centroids
    getxys <- function(ch) {
        df <- as.data.frame(centroids(ch))
        names(df) <- c('x','y')
        df
    }
    xys <- lapply(capthist, getxys)
    
    # add covariates from spatial data source if requested
    if ('spatialdata' %in% names(list(...))) {
        # force class for addCovariates secr <= 4.5.6
        class(xys) <- c('traps','list')  
        tmp <- covariates(addCovariates(xys, ...))
        xys <- mapply(cbind, xys, tmp, SIMPLIFY = FALSE)
    }
    
    # add individual covariates if present
    covar <- covariates(capthist)
    if (!is.null(covar) && all(sapply(covar,nrow)>0)) {
        xys <- mapply(cbind, xys, covar, SIMPLIFY = FALSE)
    }
    
    # add session
    addsess <- function(df, fact) {
        df$session <- fact; df
    }
    xys <- mapply(addsess, xys, factor(1:length(capthist)), SIMPLIFY = FALSE)
    
    # replicate by nocc within session, only if byoccasion
    if (byoccasion) {
        addocc <- function(df, occ) {
            df$occasion <- occ; df
        }
        n    <- sapply(capthist, nrow)
        nocc <- sapply(capthist, ncol)
        id   <-  mapply(rep, lapply(n,seq_len), nocc, SIMPLIFY = FALSE)
        occ  <- mapply(rep, lapply(nocc,seq_len), each = n, SIMPLIFY = FALSE)
        xys  <- mapply(function(x,i) x[i,], xys, id, SIMPLIFY = FALSE)
        xys  <- mapply(addocc, xys, occ, SIMPLIFY = FALSE)
    }
    
    # form single dataframe from list of dataframes
    designdata <- do.call(rbind, xys)
    
    designdata
}
##################################################

proxy.ms <- function (capthist, model = NULL, trapdesigndata = NULL, ...) {
    
    ## force list for simplicity
    ## -------------------------
    
    if (!secr::ms(capthist)) {   
        capthist <- list(capthist)
    }
    ## basics
    ## ------
    
    ni <- function (chi) {
        if (binary) {
            sum(abs(chi)>0)
        }
        else {
            sum(abs(chi))
        }
    }
    nim <- function (chi) {
        # average over occasions of sum over traps
        if (binary) {
            mean(apply(abs(chi)>0,1,sum))
        }
        else {
            mean(apply(abs(chi),1,sum))
        }
    }
    binary <- detector(traps(capthist[[1]]))[1] %in% c('single','multi','proximity')
    binom <- detector(traps(capthist[[1]]))[1] %in% c('single','multi')
    n    <- sapply(capthist, nrow)         ## number of individuals per session
    nocc <- sapply(capthist, ncol)         ## number of occasions
    K    <- sapply(traps(capthist), nrow)  ## detectors per session
    
    defaultmodel <- list(D = ~1, g0 = ~1, lambda0 = ~1, sigma = ~1, NT = ~1)
    model <- replacedefaults(defaultmodel, model)
    pmodel <- model$g0
    if ('lambda0' %in% names(model)) {
        pmodel <- model$lambda0
    }
    smodel <- model$sigma
    if (pmodel != ~1 ) {
        animaldesigndata <- detectionDesignData(capthist, byoccasion = TRUE, ...)
    }    
    if (smodel != ~1) {
        animaldesigndata.s <- detectionDesignData(capthist, byoccasion = FALSE, ...)
        freq <- unlist(sapply(capthist,function(x) apply(abs(x)>0,1,sum)-1))
        animaldesigndata.s$freq <- unname(as.numeric(freq))
    }    
    
    if (binary) {
        if (pmodel == ~1) {
            p    <- lapply(capthist, function(x) apply(x,1,nim))
            pterms <- c(cloglogp = log(-log(1-mean(unlist(p)))))
        }
        else {
            # binary animal x occasion data
            ni <- lapply(capthist, function(x) apply(abs(x),1:2,sum))
            
            animaldesigndata$ni <- unlist(lapply(ni, as.numeric)) 
            pmodel <- update(pmodel, ni ~ .)  ## ni on LHS
            if (binom) {
                # glmfit <- glm(pmodel, data = animaldesigndata, family = binomial(link = "logit"))
                glmfit <- glm(pmodel, data = animaldesigndata, family = binomial(link = "cloglog"))  # 1.4.0
            }
            else {
                # using Poisson approximation
                glmfit <- glm(pmodel, data = animaldesigndata, family = poisson(link = "log"))
            }
            pterms <- coef(glmfit)    
        }
        if (smodel == ~1) {
            sterms <- c(logRPSV = log(mean(unlist(rpsv(capthist)))))
        }
        else {
            getsi <- function (capthist) log( as.numeric(rpsvi(capthist)) + 1 )
            animaldesigndata.s$si <- unlist(lapply(capthist, getsi))
            animaldesigndata.s <- animaldesigndata.s[animaldesigndata.s$freq>0,]

            smodel <- update(smodel, si ~ .)  ## si on LHS
            glmfit <- glm(smodel, data = animaldesigndata.s, family = gaussian(link = "identity"),
                na.action = na.omit, weights = freq)
            sterms <- coef(glmfit)
                   
            # experimental Tobit model 2023-01-05
            # if (!requireNamespace("survival")) stop ("sigma model requires survival package; please install")
            # smodel <- update(smodel, survival::Surv(si, si>0, type = 'left') ~ .)  ## si on LHS
            # tobitfit <- survival::survreg(smodel, data = animaldesigndata.s,
            #     na.action = na.omit, weights = freq, dist = "gaussian", y = FALSE)
            # sterms <- coef(tobitfit)
            
        }
    }
    else {
        lambda  <- sapply(capthist, ni) / n / nocc
        pterms <- c(logL = log(mean(lambda)))
        sterms <- c(logRPSV = log(mean(unlist(rpsv(capthist)))))
    }
    
    ## Optional density model
    ## ----------------------
    if (model$D == ~1) {
        Dterms <- c(logn = log(sum(n)))
    }
    else {
        nk <- mapply(function(x,Kj) tabulate(trap(x, names = FALSE), Kj), capthist, K)
        trapdesigndata$nk <- as.numeric(nk)
        model$D <- update(model$D, nk ~ .)  ## nk on LHS
        glmfit <- glm(model$D, data = trapdesigndata, family = poisson(link = "log"))
        Dterms <- coef(glmfit)    
    }
    
    ## Optional model of nontarget data
    ## --------------------------------
    nontarget <- lapply(capthist, attr, which = 'nontarget', exact = TRUE)
    usenontarget <- !any(sapply(nontarget, is.null)) && !is.null(model$NT)
    if (usenontarget) {
        if (any(sapply(nontarget,nrow)!= K | sapply(nontarget, ncol) != nocc)) {
            stop ("invalid nontarget data in proxy.ms")
        }
        if (model$NT == ~1) {
            pdisturb <- sapply(nontarget, mean)
            if (binary) {
                pdisturb[pdisturb>=1] <- 1 - 1e-4  ## arbitrary to dodge log(0)
                NTterms <-  c(cloglogNT = log(-log(1-mean(pdisturb))))
            }
            else {
                NTterms <- c(logNT = log(mean(pdisturb)))
            }
        }
        else {
            NTk <- lapply(nontarget, apply, 1, mean)   # by detector
            trapdesigndata$NTk <- as.numeric(unlist(NTk))
            trapdesigndata$nocc <- rep(nocc, each = max(K))
            model$NT <- update(model$NT, NTk ~ .)  ## NTk on LHS
            if (binary) {
                glmfitNT <- glm(model$NT, data = trapdesigndata, weights = nocc, 
                    family = binomial(link = "cloglog"))
            }
            else {
                glmfitNT <- glm(model$NT, data = trapdesigndata, weights = nocc, 
                    family = poisson(link = "log"))
            }
            NTterms <- coef(glmfitNT)    
            names(NTterms) <- paste('NT', names(NTterms), sep = '.')
            names(NTterms) <- sub('..(Intercept))', '', names(NTterms))
            
        }
    }    
    else {
        NTterms <- NULL
    }
    
    ## compile output vector
    ## ---------------------
    c(
        Dterms,
        pterms,
        sterms,
        NTterms
    )
}
##################################################

## beta binomial
## Dorazio & Royle
## based in part on S+ code of Shirley Pledger 24/4/98

# proxy.Mhbeta <- function (capthist, ...) {
#     loglik <- function (pr) {
#         pr <- exp(pr)   ## all on log scale
#         N <- Mt1 + pr[1]
#         rat   <- pr[2] * (1- pr[2])/ pr[3]
#         if (is.na(rat) || (rat<1) || (N>maxN))  return (1e10)
#         alpha <- pr[2] * (rat-1)
#         beta  <- (1 - pr[2]) * (rat-1)
#         i <- 1:tt
#         terms <-  lgamma(alpha+i) + lgamma(beta+tt-i)
#         LL <- lgamma (N+1) - lgamma(N-Mt1+1) - lgamma(Mt1+1) +
#             N * (lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) -
#                     lgamma(alpha + beta + tt)) +
#             (N-Mt1) * (lgamma(alpha) + lgamma(beta+tt)) +
#             sum(fi*terms)
#         -LL
#     }
#     if (ms(capthist)) stop("proxy.Mhbeta is for single session only")
#     maxN <- 1e7  # little risk in hard-wiring this
#     tt <- ncol(capthist)
#     Mt1  <- nrow(capthist)
#     fi <- tabulate(apply(apply(abs(capthist),1:2,sum)>0,1,sum), nbins=tt)
#     start <- log(c(10, 1/tt, 0.2 * 1/tt * (1 - 1/tt) ))
#     fit <- nlm (p = start, f = loglik, hessian = TRUE)
#     Nhat <- exp(fit$estimate[1]) + Mt1
#     phat <- sum(abs(capthist)) / tt / Nhat
#     pr <- exp(fit$estimate) # all log scale
#     rat <- pr[2] * (1- pr[2])/ pr[3]
#     alpha <- pr[2] * (rat-1)
#     beta  <- (1 - pr[2]) * (rat-1)
#     mean <- alpha / (alpha+beta)    
#     var <- alpha * beta / (alpha+beta)^2 / (alpha+beta+1)
#     CV <- sqrt(var)/mean
#     c(
#         logN     = log(Nhat), 
#         cloglogp = log(-log(1-phat)), 
#         logCV    = log(CV), 
#         logRPSV  = log(rpsv(capthist))
#     )
# }
##################################################
