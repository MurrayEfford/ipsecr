#############################################################################
## package 'secr'
## predict.ipsecr.R
#############################################################################

## 2022-03-31 based on predict.secr

############################################################################################
predict.ipsecr <- function (object, newdata = NULL, type = c("response", "link"), se.fit = TRUE,
                          alpha = 0.05, savenew = FALSE, ...) {

    type <- match.arg(type)

    if (is.null(newdata)) {
      newdata <- makeNewData (object, ...)
    }
    
    parindices <- object$parindx
    models <- object$model
    realnames <- names(models)
    
    ## drop unused columns
    vars <- unlist(lapply(models, all.vars))
    newdata <- newdata[,names(newdata) %in% vars, drop = FALSE]

    ## allow for fixed beta parameters 
    beta <- complete.beta(object)
    beta.vcv <- complete.beta.vcv(object)

    getfield <- function (x) {
      secr.lpredictor (
        formula = models[[x]], 
        newdata = newdata,
        indx = parindices[[x]], 
        beta = beta, 
        field = x,
        beta.vcv = beta.vcv, 
        contrasts = object$details$contrasts
      )
    }
    
    predict <- sapply (realnames, getfield, simplify = FALSE)
  
    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard z!
    if (se.fit)  out <- list(nrow(newdata))
    else {
        out <- newdata
        ## add columns for real parameter estimates
        for (varname in realnames)
            out[,varname] <- rep(NA,nrow(out))
    }
  
    for (new in 1:nrow(newdata)) {
        lpred  <- sapply (predict, function(x) x[new,'estimate'])
        if (type == "response")
            Xlpred <- Xuntransform(lpred, object$link[realnames], realnames)

        # n.mash <- attr (object$capthist, 'n.mash')
        # n.clust <- length(n.mash)
        # if (unmash) {
        #     n.clust <- 1
        # }

        if (se.fit) {
            selpred <- sapply (predict,function(x) x[new,'se'])
            if (type == "response") {
                temp <- data.frame (
                    row.names = realnames,
                    link = unlist(object$link[realnames]),
                    estimate = Xlpred,
                    SE.estimate = se.Xuntransform (lpred, selpred, object$link[realnames], realnames),
                    lcl = Xuntransform(lpred-z*selpred, object$link[realnames], realnames),
                    ucl = Xuntransform(lpred+z*selpred, object$link[realnames], realnames)
                )
                ## truncate density at zero; adjust for mash(); adjust for telemetry
                if ('D' %in% row.names(temp)) {
                    temp['D', -1][temp['D',-1]<0] <- 0
                    # if (!is.null(n.mash)) {
                    #     temp['D', -1] <- temp['D', -1] / n.clust
                    # }
                }
            }
            else {
                temp <- data.frame (
                    row.names = realnames,
                    link = unlist(object$link[realnames]),
                    estimate = lpred,
                    SE.estimate = selpred,
                    lcl = lpred-z*selpred,
                    ucl = lpred+z*selpred
                )
            }

            det <- detector(traps(object$capthist))

            if (nrow(newdata)==1) out <- temp
            else {
                out[[new]] <- temp
                names(out)[new] <- paste (
                        paste(names(newdata),'=', unlist(lapply(newdata[new,],as.character)),
                        sep=' ',collapse=', '),
                    sep=',')
            }
        }
        else { # no SE; terse format
            if (type == "link") {
                out[new, (ncol(newdata)+1) : ncol(out)] <- lpred
            }
            else {
                if ('D' %in% names(Xlpred)) {
                    Xlpred['D'] <- ifelse (Xlpred['D']<0, 0, Xlpred['D'])
                    # if (!is.null(n.mash)) {
                    #     Xlpred['D'] <- Xlpred['D'] / n.clust
                    # }
                }
                out[new, (ncol(newdata)+1) : ncol(out)] <- Xlpred
            }
        }
    }
    if (savenew) attr(out, 'newdata') <- newdata
    out
}
