############################################################################################
## plot.ipsecr.R
## Method for plotting detection function from fitted ipsecr object
## 2022-04-01
############################################################################################

plot.ipsecr <- function (x, newdata = NULL, add = FALSE,
    sigmatick = FALSE, rgr = FALSE, limits = FALSE, alpha = 0.05, # bootstrap = FALSE, nboot = 10000, 
    xval = 0:200, ylim = NULL, xlab = NULL, ylab = NULL, ...)
{
    gline <- function (predicted, rowi = 1, eps = 1e-10) {
        ## eps is used to limit y to range where gradient() works
        ## may need later adjustment
        if (!is.data.frame(predicted)) {
            out <- list()
            for (i in 1:length(predicted))
                out[[i]] <- gline(predicted[[i]], i)
            names(out) <- names(predicted)
            out
        } 
        else {
            pars <- predicted[secr:::parnames(x$detectfn),'estimate']
            pars[is.na(pars)] <- unlist(x$fixed)
            dfn <- getdfn(x$detectfn)
            if (sigmatick) {
              sigma <- pars[2]
              y <- dfn(sigma, pars, x$details$cutval)
              dy <- par()$usr[4]/20
              segments (sigma, y-dy, sigma, y+dy)
            }

            y <- dfn(xval, pars, x$details$cutval)
            if (rgr) {
              y <- xval * y
              ymax <- par()$usr[4]
              lines (xval, y * 0.8 * ymax / max(y), lty = 2, ...)
            }
            else lines (xval, y, ...)

            if (limits & !rgr) {
                ## delta method variance of g()

                grad <- matrix(nrow = length(xval), ncol = length(x$fit$par))  ## beta pars
                if (is.null(newdata)) {
                    newdata <- makeNewData(x)
                }
                parnamvec <- secr:::parnames(x$detectfn)
                ## 2019-01-25 added beta0
                if (!parnamvec[1] %in% c('g0','lambda0','beta0'))
                    stop ("first detection parameter not g0 or lambda0")
                
                lkdfn <- function (beta, r) {
                    ## real g() from all beta pars and model.matrix
                    real <- numeric(length(parnamvec))
                    names(real) <- parnamvec
                    
                    for (rn in parnamvec) {
                        par.rn <- x$parindx[[rn]]
                        mat <- model.matrix(
                            x$model[[rn]], 
                            data = newdata[rowi,,drop=FALSE], 
                            contrasts = x$details$contrasts)
                        lp <- mat %*% matrix(beta[par.rn], ncol = 1)
                        real[rn] <- untransform (lp, x$link[[rn]])
                    }
                    gr <- dfn(r, real, x$details$cutval)
                    if (parnamvec[1] == 'lambda0')
                        log(-log(1-gr))
                    else 
                        logit(gr) 
                }
                
                for (i in 1:length(xval))
                    grad[i,] <- gradient (pars = x$fit$par, fun = lkdfn, r = xval[i])  ## see 'utility.R'
                
                vc <- vcov (x)
                gfn <- function(gg) {
                    gg <- matrix(gg, nrow = 1)
                    gg %*% vc %*% t(gg)
                }
                se <- apply(grad, 1, gfn)^0.5
                
                ## limits on g(r) scale 2017-10-16
                if (parnamvec[1] == 'lambda0') {
                    lcl <- ifelse ((y>eps) & (y<(1-eps)), 1 - exp(-exp(log(-log(1-y)) - z*se)), NA)
                    ucl <- ifelse ((y>eps) & (y<(1-eps)), 1 - exp(-exp(log(-log(1-y)) + z*se)), NA)
                } 
                else {
                    lcl <- ifelse ((y>eps) & (y<(1-eps)), invlogit(logit(y) - z*se), NA)
                    ucl <- ifelse ((y>eps) & (y<(1-eps)), invlogit(logit(y) + z*se), NA)
                }
                lines (xval, lcl, lty=2, ...)
                lines (xval, ucl, lty=2, ...)
            }
            
            if (limits & !rgr)
                data.frame(x=xval, y=y, lcl = lcl, ucl = ucl)
            else
                data.frame(x=xval, y=y)
        }
    }
    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard-rate z!
    temp <- predict (x, newdata)
    if (is.null(ylim)) {
        if (x$detectfn %in% c(9,10,11,12,13)) {      ## included 9 2010-11-01
            ylim <- c(0, 1)
        } 
        else {
            if (x$detectfn %in% 14:19)
                yname <- 'lambda0'
            else
                yname <- 'g0'
            getmax <- function(x) {
                g0 <- x[yname,'estimate']
                se.g0 <- x[yname,'SE.estimate']
                if (limits & is.finite(se.g0))  ## is.finite 2010-10-10
                    min(1, g0 + z * se.g0)
                else
                    g0
            }
            if (is.data.frame(temp)) maxg0 <- getmax(temp)
            else maxg0 <- max(sapply (temp, getmax))

            if (is.na(maxg0)) maxg0 <- x$fixed[[yname]]
            if (maxg0 > 0.75) maxg0 <- max(maxg0,1)
            ylim <- c(0, maxg0)
        }
    }
    if (!add) {
        if (is.null(xlab))
            xlab <- 'Distance  (m)'
        if (is.null(ylab)) {
           binomN <- ifelse(is.null(x$details$binomN),0,x$details$binomN)
           ylab <- 'Detection probability'
        }
        plot (type ='n', 0,0, xlim=range(xval), ylim=ylim,
            xlab=xlab, ylab=ylab,...)
    }
    invisible(gline(temp))
}
############################################################################################
