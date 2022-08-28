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
        # prepare animaldesigndata
        getxy <- function(ch,nocc) as.data.frame(centroids(ch))[rep(1:nrow(ch), nocc),]
        animaldesigndata <- do.call(rbind, mapply(getxy, capthist, nocc, SIMPLIFY = FALSE))
        names(animaldesigndata) <- c('x','y')
        animaldesigndata$session <- factor(rep(1:length(capthist), n*nocc))
        animaldesigndata$nocc <- rep(nocc, n*nocc)
        covar1 <- covariates(capthist[[1]])
        if (!is.null(covar1) && nrow(covar1)>0) {
            getcov <- function(ch,nocc) covariates(ch)[rep(1:nrow(ch), nocc),, drop = FALSE]
            animalcovlist <- mapply(getcov, capthist, nocc, SIMPLIFY = FALSE)
            animalcov <- do.call(rbind, animalcovlist)
            animaldesigndata <- cbind(animaldesigndata, animalcov)
        }
    }    
    if (smodel != ~1) {
        # prepare animaldesigndata.s
        getxys <- function(ch) as.data.frame(centroids(ch))
        xys <- lapply(capthist, getxys)
        class(xys) <- c('traps','list')  # fool addCovariates secr <= 4.5.6
        if ('spatialdata' %in% names(list(...))) {
            tmp <- covariates(addCovariates(xys, ...))
            xys <- mapply(cbind, xys, tmp, SIMPLIFY = FALSE)
        }
        animaldesigndata.s <- do.call(rbind, xys)
        names(animaldesigndata.s)[1:2] <- c('x','y')   # replaces meanx, meany
        animaldesigndata.s$session <- factor(rep(1:length(capthist), n))
        covar1 <- covariates(capthist[[1]])
        if (!is.null(covar1) && nrow(covar1)>0) {
            animalcov <- do.call(rbind, lapply(capthist, covariates))
            animaldesigndata.s <- cbind(animaldesigndata.s, animalcov)
        }
        animaldesigndata.s$freq <- unlist(sapply(capthist,function(x) apply(abs(x)>0,1,sum)-1))
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
                glmfit <- glm(pmodel, data = animaldesigndata, family = binomial())
            }
            else {
                # ntraps <- sapply(capthist, function(x) dim(x)[3])
                # animaldesigndata$ntraps <- rep(ntraps, n*nocc)
                # glmfit <- glm(pmodel, data = animaldesigndata, family = binomial(), weights = ntraps)
                # using Poisson approximation
                glmfit <- glm(pmodel, data = animaldesigndata, family = poisson())
            }
            pterms <- coef(glmfit)    
        }
        if (smodel == ~1) {
            sterms <- c(logRPSV = log(mean(unlist(rpsv(capthist)))))
        }
        else {
            # si <- lapply(capthist, rpsvi)
            # animaldesigndata.s$si <- unlist(lapply(si, as.numeric)) 
            getsi <- function (capthist) log( as.numeric(rpsvi(capthist)) + 1 )
            animaldesigndata.s$si <- unlist(lapply(capthist, getsi))
            smodel <- update(smodel, si ~ .)  ## si on LHS
            freq <- animaldesigndata.s$freq
            glmfit <- glm(smodel, data = animaldesigndata.s, 
                na.action = na.omit, weights = freq)
            sterms <- coef(glmfit)
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
        glmfit <- glm(model$D, data = trapdesigndata, family = poisson())
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
                    family = poisson())
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
