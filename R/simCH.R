###############################################################################
## package 'ipsecr'
## simCH.R
## 2022-05-10, 2022-06-12, 2022-06-14, 2022-06-16
###############################################################################

# function simCH is used by ipsecr.fit for CHmethod 'internal'

simCH <- function (traps, popn, detectfn, detectpar, NT, noccasions, details = list()) {
    if (ms(traps)) {
        # if detectpar vary across sessions then should be list of lists
        if (!is.list(detectpar[[1]])) detectpar <- list(detectpar)
        if (length(NT) == 0 || length(unlist(sapply(NT, length))) < length(traps)) NT <- 0
        tmp <- mapply(simCH, 
            traps = traps, 
            popn = popn, 
            detectpar = detectpar, 
            NT = as.data.frame(NT),  # detector x session 
            noccasions = noccasions,
            MoreArgs = list(detectfn = detectfn, details = details), 
            SIMPLIFY = FALSE)
        MS.capthist(tmp)
    }
    else {
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
        if (is.null(NT) || all(NT<=0)) {
            nontargetcode <- 0
        }
        else {
            validnontargettype <- c('exclusive', 'truncated','erased','independent','dependent')
            # nontargettype defaults to 'exclusive'
            details$nontargettype <- match.arg(details$nontargettype, validnontargettype)
            nontargetcode <- match(details$nontargettype, validnontargettype)
            if (detectcode %in% c(0,1,2) && nontargetcode == 1) {
                warning("exclusive interference not possible with detector type, assuming 'truncated'")
                nontargetcode <- 2  # truncated
            }
            NT <- pmax(NT,0)    # ensure non-negative
            NT <- rep(NT,length.out = K)
        }

        temp <- CHcpp(
            as.matrix(popn), 
            as.matrix(traps), 
            as.matrix(usge),
            as.double(detpar), 
            as.double(NT),   # vector 2022-06-13
            as.integer(detectfn), 
            as.integer(detectcode), 
            as.integer(nontargetcode),
            0, 0, rep(0, noccasions))
     
        if (temp$resultcode != 0) {
            stop ("simulated detection failed, code ", temp$resultcode)
        }
        w <- array(temp$CH, dim = c(nrow(popn), noccasions, K),
            dimnames = list(1:nrow(popn), 1:noccasions, NULL))
        w <- w[apply(w,1,sum)>0,,, drop = FALSE] 
        class(w)   <- 'capthist'
        if (nontargetcode > 0) {
            attr(w, 'nontarget') <- temp$nontarget
        }
        traps(w)   <- traps
        w
    }
}
