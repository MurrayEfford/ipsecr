coef.ipsecr <- function (object, alpha=0.05, ...) {
    beta   <- object$beta
    if (!is.null(object$beta.vcv))
        sebeta <- suppressWarnings(sqrt(diag(object$beta.vcv)))
    else sebeta <- rep(NA, length(beta))
    z <- abs(qnorm(1-alpha/2))
    temp <- data.frame(
        row.names = object$betanames,
        beta    = beta,
        SE.beta = sebeta,
        lcl = beta - z*sebeta,
        ucl = beta + z*sebeta
    )
    attr(temp, 'alpha') <- alpha
    temp
}
############################################################################################

