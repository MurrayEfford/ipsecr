###############################################################################
## package 'ipsecr'
## simCH.R
## 2022-05-10
###############################################################################

# function simCH is used by ipsecr.fit for CHmethod 'internal'

simCH <- function (traps, popn, detectfn, detectpar, noccasions) {
    
    K <- nrow(traps)
    usge <- usage(traps)
    if (is.null(usge)) {
        usge <- matrix(1, K, noccasions)
    }
    detectcode <- switch(detector(traps)[1], single = -1, 
        multi = 0, proximity = 1, count = 2, capped = 8, 9)
    if (detectcode == 9) stop ("unsupported detector type")
    # optional nontarget rate
    lambdak <- detectpar[['lambdak']][1]
    if (is.null(lambdak)) lambdak <- -1 
    
    temp <- CHcpp(
        as.matrix(popn), 
        as.matrix(traps), 
        as.matrix(usge),
        as.integer(detectfn), 
        as.integer(detectcode), 
        unlist(detectpar[parnames(detectfn)]),  # robust to order of detectpar
        as.double(lambdak),
        0, 0, 0)
    
    if (temp$resultcode != 0) {
        stop ("simulated detection failed, code ", temp$resultcode)
    }
    npop <- nrow(popn)
    w <- array(temp$value, dim = c(noccasions, K, npop), 
        dimnames = list(1:noccasions, NULL, 1:npop))
    w <- aperm(w, c(3,1,2))
    if (lambdak > 0) {
        # retrieve nontarget from last row
        nontarget <- temp$nontarget
    }
    else {
        nontarget <- NULL
    }
    w <- w[apply(w,1,sum)>0,,, drop = FALSE] 
    class(w)   <- 'capthist'
    attr(w, 'nontarget') <- nontarget
    traps(w)   <- traps
    w
}
