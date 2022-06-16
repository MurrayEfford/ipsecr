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

proxy.ms <- function (capthist, model = list(D = ~1, NT = ~1), trapdesigndata = NULL, ...) {
    
    ## force list for simplicity
    ## -------------------------
    
    if (!secr::ms(capthist)) {   
        capthist <- list(capthist)
    }
    
    ## basics
    ## ------
    
    ni <- function (chi) {
        # ch <- abs(chi) > 0
        # sum(apply(ch, 2, sum))
        if (binary) {
            sum(abs(chi)>0)
        }
        else {
            sum(abs(chi))
        }
    }
    binary <- detector(traps(capthist[[1]]))[1] %in% c('single','multi','proximity')
    n    <- sapply(capthist, nrow)         ## number of individuals per session
    nocc <- sapply(capthist, ncol)         ## number of occasions
    K    <- sapply(traps(capthist), nrow)  ## detectors per session
    
    if (binary) {
        p    <- sapply(capthist, ni) / n / nocc
        pterms <- c(cloglogp = log(-log(1-mean(p))))
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
    usenontarget <- !any(sapply(nontarget, is.null))
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
# 
# ms 
# meanSD <- lapply(traps(capthist), getMeanSD)
# not ms 
# meanSD <- getMeanSD(traps(capthist))
# trapdesigndata <- D.designdata(traps(capthist), model$D, 1, sessionlevels, sessioncov, meanSD)

# system.time(proxy.ms(ovenCHp, model=list(D= ~y)))
# ipsecr.fit(ovenCHp, mask=msk, model=list(D~y), proxyfn = proxy.ms, ncores=4)
