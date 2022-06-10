###############################################################################
## package 'ipsecr'
## proxyfn.R
## 2022-05-08
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
    c(N=estimates[1], oddsp = odds(estimates[2]), rpsv=RPSV(capthist, CC = TRUE))
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
        logrpsv= log(RPSV(capthist, CC = TRUE))
    )
}
##################################################

proxy.nt <- function (capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...) {
    N.estimator <- tolower(N.estimator)
    N.estimator <- match.arg(N.estimator)
    ## capthist single-session only; ignoring losses
    n <- nrow(capthist)         ## number of individuals
    ch <- abs(capthist)>0
    nocc <- ncol(capthist)      ## number of occasions
    ni <- apply(ch, 2, sum)     ## individuals on each occasion
    nontarget <- attr(capthist, 'nontarget', exact = TRUE)
    if (is.null(nontarget)) {
        stop ("proxy.nt requires nontarget")
    }
    if (nrow(nontarget)!= dim(capthist)[3] || ncol(nontarget) != nocc) {
        stop ("invalid nontarget data in proxy.nt")
    }
    pdisturb <- mean(nontarget)
    if (pdisturb>=1) {
        pdisturb <- 1 - 1e-4  ## arbitrary to dodge log(0)
    }
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
        logrpsv = log(RPSV(capthist, CC = TRUE)),
        cloglogNT = log(-log(1-pdisturb))
    )
}
##################################################

proxy.ms <- function (capthist, model, ...) {
    if (!inherits(capthist, 'list')) stop ("proxy.s expects a multi-session capthist")
    n <- sapply(capthist, nrow) ## number of individuals per session
    ni <- function (chi) {
        ch <- abs(chi) > 0
        sum(apply(ch, 2, sum))
    }
    nocc <- sapply(capthist, ncol)      ## number of occasions
    nsess <- length(n)
    p <- sapply(capthist, ni) / n / nocc
    rpsv <- unlist(secr::RPSV(capthist, CC = TRUE))
    D0 <- log(n[1])
    if (model$D == ~1) {  # constant density
        Dterms <- NULL
    }
    else if (model$D == ~session) {   # session-specific density
        Dterms <- log(n[-1]/n[1])
    }
    else if (model$D == ~Session) {   # log-linear trend in density
        Dlm <- coef(lm(log(n)~I(0:(nsess-1))))
        D0 <- Dlm[1]
        Dterms <- Dlm[2]
    }
    else stop ("proxyfn.ms does not accept model ", model$D)
    c(
        D0,
        Dterms,
        logp = log(-log(1-mean(p))),
        logRPSV = log(mean(rpsv))
    )
}

