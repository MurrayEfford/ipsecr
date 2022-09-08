## 2022-08-24 miscellaneous
## 2022-08-30 explicit RNGkind in set.seed

library(ipsecr)
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

###############################################################################

set.seed(1235)
msk <- make.mask(traps(subset(ovenCH,sessions=1:2)), buffer = 100, nx = 16)

pop1 <- simpop(
    mask = msk, 
    D = 1,
    N = 100, 
    details = list(distribution = 'poisson'))

pop2 <- simpop(
    mask = msk, 
    D = 1,
    N = 100, 
    details = list(distribution = 'even'))

# summarise nearest neighbour distances of simulated populations 
dfn <- function (pop) {
    if (ms(pop)) {
        unlist(lapply(pop, dfn))
    }
    else {
        d <- as.matrix(dist(pop))
        diag(d) <- NA
        mind <- apply(d,1,min, na.rm = TRUE)
        c(mean = mean(mind), sd = sd(mind))
    }
}

test_that("simulated Poisson and even populations have expected spacing", {
    # Poisson
    expect_equal(dfn(pop1), 
        c(24.401273, 14.376827, 25.382389, 14.218940),
        tolerance = 1e-4, check.attributes = FALSE)
    # even
    expect_equal(dfn(pop2), 
        c(45.77484, 0.00000, 45.77484, 0.00000),
        tolerance = 1e-4, check.attributes = FALSE)
})
###############################################################################

