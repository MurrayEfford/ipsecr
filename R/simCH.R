###############################################################################
## package 'ipsecr'
## simCH.R
## 2022-05-10, 2022-06-12, 2022-06-14, 2022-06-16
## 2022-09-07 RcppArmadillo armaCHcpp()
###############################################################################

# function simCH is used by ipsecr.fit for CHmethod 'internal'

simCH <- function (traps, popn, detectfn, detparmat, noccasions, NT = NULL, 
    details = list()) {
    if (ms(traps)) {
        # detparmat should be list of matrices
        if (!is.list(detparmat)) detparmat <- list(detparmat)
        if (length(NT) == 0 || length(unlist(sapply(NT, length))) < length(traps)) NT <- 0
        tmp <- mapply(simCH, 
            traps = traps, 
            popn = popn, 
            detparmat = detparmat,     # automatically replicated if single component
            noccasions = noccasions,
            NT = as.data.frame(NT),  # detector x session 
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
        
        if (!is.matrix(detparmat)) {
            detparmat <- matrix(unlist(detparmat), byrow = TRUE, 
                nrow = nrow(popn), ncol = length(unlist(detparmat)))
        }
        # optional nontarget rate
        if (is.null(NT) || all(NT<=0)) {
            if (is.null(NT)) NT <- 0
            NT <- rep(NT, length.out = K)
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
            NT <- pmax(NT, 0)    # ensure non-negative
            NT <- rep(NT, length.out = K)
        }

        binomN <- 0  # Poisson count if count detector
        binomN <- rep(binomN, length.out = noccasions)
        
        #-----------------------------------------------------------
        # new code 1.3.0 2022-09-07
        
        w <- armaCHcpp(
            as.matrix(edist(popn, traps)),
            as.matrix(usge),
            as.matrix(detparmat),
            as.double(NT),
            as.integer(binomN),
            as.integer(detectfn),
            as.integer(detectcode),
            as.integer(nontargetcode))
        
        dimnames(w) <- list(1:nrow(w), 1:noccasions, NULL)
        if (nrow(w)>0) {
            ## strip nontarget detections from last row of array
            if (nontargetcode>0) {
                nontarget <- t(array(w[nrow(w),,], dim=dim(w)[2:3]))
                w <- w[-nrow(w),,, drop = FALSE]
            }
            else {
                nontarget <- NULL
            }
        }
        #-----------------------------------------------------------
        
        ## drop empty capture histories
        w <- w[apply(w,1,sum)>0,,, drop = FALSE]
        
        class(w)   <- 'capthist'
        if (nontargetcode > 0) {
            attr(w, 'nontarget') <- nontarget
        }
        traps(w)   <- traps
        w
      
    }
}
