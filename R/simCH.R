###############################################################################
## package 'ipsecr'
## simCH.R
## 2022-05-10
###############################################################################

# function simCH is used by ipsecr.fit for CHmethod 'internal'

simCH <- function (traps, popn, detectfn, detectpar, noccasions, details = list()) {
    K <- nrow(traps)
    usge <- usage(traps)
    if (is.null(usge)) {
        usge <- matrix(1, K, noccasions)
    }
    detectcode <- switch(detector(traps)[1], single = -1, 
        multi = 0, proximity = 1, count = 2, capped = 8, 9)
    if (detectcode == 9) stop ("unsupported detector type")
    
    detpar <- unlist(detectpar[parnames(detectfn)])  # robust to order of detectpar
    # optional nontarget rate
    lambdak <- detectpar[['lambdak']][1]
    if (is.null(lambdak)) lambdak <- -1 
    
    if (lambdak>0) {
        validnontargettype <- c('exclusive', 'truncated','erased','independent')
        # nontargettype defaults to 'exclusive'
        details$nontargettype <- match.arg(details$nontargettype, validnontargettype)
        nontargetcode <- match(details$nontargettype, validnontargettype)
        if (detectcode %in% c(0,1,2) && nontargetcode == 1) {
            warning("exclusive interference not possible with detector type, assuming 'truncated'")
            nontargetcode <- 2  # truncated
        }
    }
    else {
        nontargetcode <- 0
    }
    
    temp <- CHcpp(
        as.matrix(popn), 
        as.matrix(traps), 
        as.matrix(usge),
        as.double(detpar), 
        as.integer(detectfn), 
        as.integer(detectcode), 
        as.double(lambdak),
        as.integer(nontargetcode),
        0, 0, 0)
    
    if (temp$resultcode != 0) {
        stop ("simulated detection failed, code ", temp$resultcode)
    }
    npop <- nrow(popn)
    # w <- array(temp$CH, dim = c(noccasions, K, npop), 
    #     dimnames = list(1:noccasions, NULL, rownames(pop)))
    # w <- aperm(w, c(3,1,2))
    w <- array(temp$CH, dim = c(npop, noccasions, K),
        dimnames = list(rownames(pop), 1:noccasions, NULL))
    w <- w[apply(w,1,sum)>0,,, drop = FALSE] 
    class(w)   <- 'capthist'
    if (lambdak > 0) {
        nontarget <- temp$nontarget
    }
    else {
        nontarget <- NULL
    }
    attr(w, 'nontarget') <- nontarget
    traps(w)   <- traps
    w
}
