###############################################################################
## package 'ipsecr'
## plotProxy.R
# 2022-07-06
###############################################################################

plotProxy <- function (
    parameter  = 'sigma', 
    proxyfn    = proxy.ms, 
    traps,
    mask,
    detectfn   = 'HHN',
    noccasions = 5,
    basepar    = list(),
    xvals      = NULL, 
    nrepl      = 100, 
    add        = FALSE,
    trend      = TRUE, 
    points     = FALSE,
    boxplot    = TRUE, 
    boxplotargs = list(),
    link       = 'log',
    details    = NULL, 
    ...
) {

    onecombo <- function (parval) {
        # parval is (usually) vector of D, lambda0, sigma
        out <- numeric(nrepl)
        N <- parval[1] * maskarea(mask)   # expected number in popn
        for (r in 1:nrepl) {
            pop <- simpop(mask, D = D, N = N)
            detparmat <- matrix(parval[-1], byrow = TRUE, nrow = nrow(pop), 
                ncol = length(parval)-1)
            ch <- simCH(traps, pop, detectfn, detparmat, noccasions)
            out[r] <- proxyfn(ch)[proxy]
        }
        out
    }
    detectfn <- detectionfunctionnumber(detectfn)
    D <- rep(1,nrow(mask))  # relative not absolute
    if (is.null(xvals)) {
        xvals <- seq(0.8,1.2,0.1) * basepar[[parameter]]
    }
    detpar <- if(detectfn<14) 'g0' else 'lambda0'
    parameters <- c('D',detpar,'sigma','z')
    if (!parameter %in% parameters) stop ("parameter not recognised")
    proxy <- names(proxyfn(secr::captdata))[match(parameter, parameters)]
    basepar[[parameter]] <- xvals
    combo <- expand.grid(basepar)
    tmp <- apply(combo, 1, onecombo)
    tmp[!is.finite(tmp)] <- NA
    dimnames(tmp) <- list(NULL, xvals)
    if (boxplot || points) {
        trxvals <- transform(xvals, link)
        if (!add) {
            xlab <- if (link=='identity') parameter else 
                paste0(parameter, ' (', link, ' scale)')
            xlim <- range(trxvals, na.rm = TRUE)
            xlim <- xlim + diff(xlim)/20 * c(-1,1) # extend by 5%
            plot(
                range(trxvals, na.rm = TRUE),
                range(tmp, na.rm = TRUE), 
                xlim = xlim,
                type = 'n', 
                xlab = xlab, 
                ylab = proxy, 
                axes = FALSE, 
                ...)
            axis(1, at = trxvals, labels = xvals)
            axis(2)
            box()
        }
        if (boxplot) {
            boxwidth <- diff(par()$usr[1:2])/30
            boxplotargs <- replacedefaults(list(x = tmp, add = TRUE, 
                at = trxvals, col = NA, pars = list(boxwex = boxwidth)), 
                boxplotargs)
            do.call('boxplot', boxplotargs)
        }
        if (points) {
            apply(tmp, 1, points, x = trxvals)
        }
        if (trend) {
            x <- rep(trxvals, each = nrepl)
            y <- as.numeric(tmp)
            abline(lm(y~x))
        }
    }
    invisible(tmp)
}