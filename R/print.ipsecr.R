print.ipsecr <- function (x, newdata = NULL, alpha = 0.05, deriv = FALSE, call = TRUE, ...) {
    
    if (!is.null(x$call) & call) {
        cat ('\n')
        print(x$call)
    }
    cat ('ipsecr ', x$version, ', ', x$starttime, '\n', sep='')
    cat ('\n')
    
    print(summary(traps(x$capthist)), terse=TRUE)
    
    cat ('\n')
    ###################
    ## Data description
    
    if (ms(x$capthist)) {
        print (summary(x$capthist, terse = TRUE))
        det <- detector(traps(x$capthist)[[1]])
    } 
    else {
        det <- detector(traps(x$capthist))
        n  <- nrow(x$capthist)     # number caught
        if (length(dim(x$capthist))>2)
            ncapt <- sum(abs(x$capthist))
        else
            ncapt <- sum(abs(x$capthist)>0)
        cat ('N animals       : ', n, '\n')
        cat ('N detections    : ', ncapt, '\n')
        cat ('N occasions     : ', ncol(x$capthist), '\n')
    }
    if (any(det %in% .localstuff$countdetectors)) {
        cat ('Count model     :  ')
        if (x$details$binomN == 0) cat ('Poisson \n')
        else if (x$details$binomN == 1) cat ('Binomial, size from usage\n')
        else if (x$details$binomN < 0) cat ('Negative binomial k = ', abs(x$details$binomN), '\n')
        else if (x$details$binomN > 1) cat('Binomial', x$details$binomN, '\n')
    }
    
    if (!ms(x$capthist)) {
        cat ('Mask area       : ', maskarea(x$mask), 'ha \n')
    }
    
    ####################
    ## Model description
    
    Npar <- nparameters(x)   ## see utility.R
    cat ('\n')
    cat ('Model           : ', model.string(x$model, x$details$userDfn), '\n')
    
    cat ('Fixed (real)    : ', fixed.string(x$fixed), '\n')
    cat ('Detection fn    : ', detectionfunctionname(x$detectfn), '\n')
    cat ('Distribution    : ', x$details$distribution, '\n')
    cat ('N parameters    : ', Npar, '\n')
    cat ('\n')
    cat ('Beta parameters (coefficients)', '\n')
    
    print(coef(x), ...)
    
    cat ('\n')
    cat ('Variance-covariance matrix of beta parameters', '\n')
    print (x$beta.vcv, ...)
    
    cat ('\n')
    cat ('Variance bootstrap \n')
    print(x$variance.bootstrap)

    # scale newdata covariates... NOT FINISHED 10 05 08
    meanSD <- attr(x$mask,'meanSD',exact = TRUE)
    if (!is.null(newdata)) {
        for (i in 1:length(newdata)) {
            ind <- match (names(newdata[i]),names(meanSD))
            if (ind>0 & !is.na(meanSD[1,ind]))
                newdata[[i]] <- (newdata[[i]] - meanSD[1,ind]) / meanSD[2,ind]
        }
    }
    
    cat ('\n')
    cat ('Fitted (real) parameters evaluated at base levels of covariates', '\n')
    
    temp <- predict (x, newdata, type = "response", alpha = alpha, 
        se.fit = x$details$var.nsim>1)
    print(temp, ...)
    
    cat ('\n')
}
#################################################################################
