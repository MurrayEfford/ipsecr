###############################################################################
## package 'ipsecr'
## simpop.R
## 2022-06-13
###############################################################################

# function simpop is used by ipsecr.fit for popmethod 'internal'

simpop <- function (mask, D, N, distribution) {
    if (ms(mask)) {
        tmp <- mapply(simpop, 
            mask = mask, 
            D = as.data.frame(D),    # each column of matrix
            N = N, 
            MoreArgs = list(distribution = distribution), 
            SIMPLIFY = FALSE
        )
        class(tmp) <- c('popn','list')
        tmp
    }
    else {
        if (distribution == 'even') {
            bounds <- apply(mask,2,range)
            popevencpp(
                as.matrix(bounds), 
                as.integer(N)
            )
        }
        else {
            popcpp(
                as.matrix(mask), 
                as.double(D/sum(D)), 
                as.double(spacing(mask)/100), 
                as.integer(N)
            )
        }
    }
}
