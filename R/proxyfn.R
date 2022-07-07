###############################################################################
## package 'ipsecr'
## proxyfn.R
## 2022-05-08
## 2022-06-11 proxy.ms
## 2022-06-15 proxy.ms extended for count detectors
###############################################################################

proxyfn0 <- function (capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...) {
    N.estimator <- tolower(N.estimator)
    N.estimator <- match.arg(N.estimator)
    ## capthist single-session only; ignoring losses
    n <- nrow(capthist)         ## number of individuals
    ch <- abs(capthist)>0
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
            Mb(ui)
        }
        else if (N.estimator == "jackknife") {
            fi <- tabulate(apply(ch,1,sum), nbins = nocc)
            Mh(fi)
        }
    }
    c(N=estimates[1], oddsp = odds(estimates[2]), rpsv=rpsv(capthist))
}
##################################################

proxyfn1 <- function (capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...) {
    N.estimator <- tolower(N.estimator)
    N.estimator <- match.arg(N.estimator)
    ## capthist single-session only; ignoring losses
    n <- nrow(capthist)         ## number of individuals
    ch <- abs(capthist)>0
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
            Mb(ui)
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

proxy.ms <- function (capthist, model = NULL, trapdesigndata = NULL, 
    animaldesigndata = NULL, ...) {
    
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
    if (is.null(pmodel)) {
        pmodel <- model$lambda0
    }
    
    if (binary) {
        # p    <- sapply(capthist, ni) / n / nocc
        if (pmodel == ~1) {
            p    <- lapply(capthist, function(x) apply(x,1,nim))
            pterms <- c(cloglogp = log(-log(1-mean(unlist(p)))))
        }
        else {
            getxy <- function(ch,nocc) as.data.frame(centroids(ch))[rep(1:nrow(ch), nocc),]
            animaldesigndata <- do.call(rbind, mapply(getxy, capthist, nocc, SIMPLIFY = FALSE))
            names(animaldesigndata) <- c('x','y')
            animaldesigndata$session <- rep(1:length(capthist), each = n*nocc)
            animaldesigndata$nocc <- rep(nocc, each = n*nocc)
            covar1 <- covariates(capthist[[1]])
            if (!is.null(covar1) && nrow(covar1>0)) {
                getcov <- function(ch,nocc) covariates(ch)[rep(1:nrow(ch), nocc),]
                animalcov <- do.call(rbind, mapply(getcov, capthist, nocc, SIMPLIFY = FALSE))
                animaldesigndata <- cbind(animaldesigndata, animalcov)
            }
            # binary animal x occasion data
            ni <- lapply(capthist, function(x) apply(x>0,1:2,sum))

            animaldesigndata$ni <- unlist(lapply(ni, as.numeric)) 
            pmodel <- update(pmodel, ni ~ .)  ## ni on LHS
            if (binom) {
                glmfit <- glm(pmodel, data = animaldesigndata, family = binomial())
            }
            else {
                stop("not ready for non-binom")
            }
            pterms <- coef(glmfit)    
        }
    }
    else {
        lambda  <- sapply(capthist, ni) / n / nocc
        pterms <- c(logL = log(mean(lambda)))
    }
    
    rpsv <- unlist(rpsv(capthist))
    
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
        logRPSV = log(mean(rpsv)),
        NTterms
    )
}
##################################################

## beta binomial
## Dorazio & Royle
## based in part on S+ code of Shirley Pledger 24/4/98

proxy.Mhbeta <- function (capthist, ...) {
    loglik <- function (pr) {
        pr <- exp(pr)   ## all on log scale
        N <- Mt1 + pr[1]
        rat   <- pr[2] * (1- pr[2])/ pr[3]
        if (is.na(rat) || (rat<1) || (N>maxN))  return (1e10)
        alpha <- pr[2] * (rat-1)
        beta  <- (1 - pr[2]) * (rat-1)
        i <- 1:tt
        terms <-  lgamma(alpha+i) + lgamma(beta+tt-i)
        LL <- lgamma (N+1) - lgamma(N-Mt1+1) - lgamma(Mt1+1) +
            N * (lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) -
                    lgamma(alpha + beta + tt)) +
            (N-Mt1) * (lgamma(alpha) + lgamma(beta+tt)) +
            sum(fi*terms)
        -LL
    }
    if (ms(capthist)) stop("proxy.Mhbeta is for single session only")
    maxN <- 1e7  # little risk in hard-wiring this
    tt <- ncol(capthist)
    Mt1  <- nrow(capthist)
    fi <- tabulate(apply(apply(abs(capthist),1:2,sum)>0,1,sum), nbins=tt)
    start <- log(c(10, 1/tt, 0.2 * 1/tt * (1 - 1/tt) ))
    fit <- nlm (p = start, f = loglik, hessian = TRUE)
    Nhat <- exp(fit$estimate[1]) + Mt1
    phat <- sum(abs(capthist)) / tt / Nhat
    pr <- exp(fit$estimate) # all log scale
    rat <- pr[2] * (1- pr[2])/ pr[3]
    alpha <- pr[2] * (rat-1)
    beta  <- (1 - pr[2]) * (rat-1)
    mean <- alpha / (alpha+beta)    
    var <- alpha * beta / (alpha+beta)^2 / (alpha+beta+1)
    CV <- sqrt(var)/mean
    c(
        logN     = log(Nhat), 
        cloglogp = log(-log(1-phat)), 
        logCV    = log(CV), 
        logRPSV  = log(rpsv(capthist))
    )
}
##################################################

# 
# ms 
# meanSD <- lapply(traps(capthist), getMeanSD)
# not ms 
# meanSD <- getMeanSD(traps(capthist))
# trapdesigndata <- D.designdata(traps(capthist), model$D, 1, sessionlevels, sessioncov, meanSD)

# system.time(proxy.ms(ovenCHp, model=list(D= ~y)))
# ipsecr.fit(ovenCHp, mask=msk, model=list(D~y), proxyfn = proxy.ms, ncores=4)
