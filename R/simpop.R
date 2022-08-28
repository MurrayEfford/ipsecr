###############################################################################
## package 'ipsecr'
## simpop.R
## 2022-06-13, 2022-07-04
###############################################################################

# function simpop is used by ipsecr.fit for popmethod 'internal'
# 2022-07-04 distribution replaced by details$distribution
# 2022-07-04 class popn

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
        if (is.null(details$distribution) || 
                tolower(details$distribution) == 'poisson') {
            D <- rep(D, length.out = nrow(mask))  # 2022-07-04
            pop <- popcpp(
                as.matrix(mask), 
                as.double(D/sum(D)), 
                as.double(spacing(mask)/100), 
                as.integer(N)
            )
        }
        else {  # details$distribution == 'even' 
            bounds <- apply(mask,2,range)
            pop <- popevencpp(
                as.matrix(bounds), 
                as.integer(N)
            )
            pop <- pop[!is.na(pop[,1]),]
        }
        dimnames(pop)[[2]] <- c('x','y')
        pop <- as.data.frame(pop)
        attr(pop, 'boundingbox') <- attr(mask, 'boundingbox')
        class(pop) <- c('popn','data.frame')
        
        # add mask covariates (could make conditional on details argument)
        if (!is.null(covariates(mask))) {
            pop <- addCovariates(pop, mask)   # add all
        }
        
        pop
    }
}
