###############################################################################
## package 'ipsecr'
## simpop.R
## 2022-06-13, 2022-07-04
###############################################################################

# function simpop is used by ipsecr.fit for popmethod 'internal'
# 2022-07-04 distribution replaced by details$distribution

simpop <- function (mask, D, N, details = list()) {
    if (ms(mask)) {
        tmp <- mapply(simpop, 
            mask = mask, 
            D = as.data.frame(D),    # each column of matrix
            N = N, 
            MoreArgs = list(details = details), 
            SIMPLIFY = FALSE
        )
        class(tmp) <- c('popn','list')
        tmp
    }
    else {
        if (is.null(details$distribution) || details$distribution == 'even') {
            popcpp(
                as.matrix(mask), 
                as.double(D/sum(D)), 
                as.double(spacing(mask)/100), 
                as.integer(N)
            )
        }
        else {  # details$distribution == 'even' 
            bounds <- apply(mask,2,range)
            popevencpp(
                as.matrix(bounds), 
                as.integer(N)
            )
        }
    }
}