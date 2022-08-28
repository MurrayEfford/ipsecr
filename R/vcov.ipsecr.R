############################################################################################
## vcov.ipsecr.R
## S3 method cf vcov.ipsecr
## 2022-08-24  bug fixed; newdata not fully tested
############################################################################################

vcov.ipsecr <- function (object, realnames = NULL, newdata = NULL, byrow = FALSE, ...) {
    ## return either the beta-parameter variance-covariance matrix
    ## or vcv each real parameters between points given by newdata (byrow = TRUE)
    ## or vcv for real parameters at points given by newdata (byrow = TRUE)
    
    if (is.null(dimnames(object$beta.vcv)))
        dimnames(object$beta.vcv) <- list(object$betanames, object$betanames)
    
    if (is.null(realnames))
        ## average beta parameters
        return( object$beta.vcv )
    else {
        if (is.null(newdata)) {
            newdata <- makeNewData (object)
        }
        ## average real parameters
        ## vcv among multiple rows
        if (byrow) {
            ## need delta-method variance of reals given object$beta.vcv & newdata
            nreal <- length(realnames)
            nbeta <- length(object$fit$par)
            
            rowi <- function (newdatai) {
                reali <- function (beta, rn) {
                    ## real from all beta pars eval at newdata[i,]
                    par.rn <- object$parindx[[rn]]
                    mat <- model.matrix(object$model[[rn]], data = newdatai, 
                        contrasts = object$details$contrasts)
                    lp <- mat %*% matrix(beta[par.rn], ncol = 1)
                    untransform (lp, object$link[[rn]])
                }
                grad <- matrix(nrow = nreal, ncol = nbeta)
                dimnames(grad) <- list(realnames, object$betanames)
                for (rn in realnames)
                    grad[rn,] <- fdHess (pars = object$beta, fun = reali, rn = rn)$gradient
                vcv <- grad %*% object$beta.vcv %*% t(grad)
                vcv
            }
            
            vcvlist <- list(nrow(newdata))
            for (i in 1:nrow(newdata)) vcvlist[[i]] <- rowi(newdata[i,])
            if (length(vcvlist) == 1) vcvlist <- vcvlist[[1]]
            return(vcvlist)
        }
        else {
            newdata <- as.data.frame(newdata)
            rownames <- apply(newdata, 1, function(x) paste(names(newdata), '=', x, sep='',
                collapse=','))
            vcvlist <- list()
            for (rn in realnames) {
                par.rn <- object$parindx[[rn]]
                mat <- model.matrix(
                    object$model[[rn]], 
                    data = newdata,
                    contrasts = object$details$contrasts)
                lp <- mat %*% matrix(object$beta[par.rn], ncol = 1)
                real <- untransform (lp, object$link[[rn]])
                real <- as.vector(real)
                ## from Jeff Laake's 'compute.real' in RMark...
                deriv.real <- switch(object$link[[rn]],
                    logit = mat * real * (1-real),
                    log = mat * real,
                    identity = mat,
                    sin = mat * cos(asin(2*real-1))/2)
                vcvlist[[rn]] <- deriv.real %*% object$beta.vcv[par.rn, par.rn] %*% t(deriv.real)
                dimnames(vcvlist[[rn]]) <- list(rownames, rownames)
            }
            names (vcvlist) <- realnames
            return (vcvlist)
        }
        ## DIFFERENT VARIANCE TO ipsecr.lpredictor for sigma because there use se.Xuntransfom
    }
}
