## 2022-04-17 start
## 2022-06-13 new proxy.ms
## 2022-06-21 set RCPP_PARALLEL_BACKEND
## 2022-08-24 uses summary.ipsecr
## 2022-08-24 vcov.ipsecr
## 2022-08-24 slimmed down test fit, better random seed handling
## 2022-08-30 explicit RNGkind 
## 2022-09-01 test RNG

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

###############################################################################
test_that("random number generator stable", {
    set.seed(1235)
    expect_equal(rnorm(1), -0.6979879, tolerance = 1e-6)
})
###############################################################################

library(ipsecr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

###############################################################################
set.seed(1235)

ch <- subset(captdata, traps = 1:40)
msk <- make.mask(traps(ch), buffer = 100, nx = 32)
fit2 <- ipsecr.fit(ch, detectfn = 'HHN', mask = msk, ncores = 1, verbose = TRUE,
    details=list(var.nsim = 500, boxsize1 = 0.3, dev.max = c(0.005, 0.005)))

vcv <- vcov(fit2, realnames = c('D', 'lambda0', 'sigma'))
fitsum <- summary(fit2)
###############################################################################

test_that("correct single-catch estimate", {
    expect_equal(fitsum$predicted[,'estimate'], 
        # c(6.051439, 0.247498, 28.834703),  # 1.2.0
        c(6.0684585, 0.2425415, 29.0726786),  # 1.2.1
        tolerance = 1e-4, check.attributes = FALSE)
})


###############################################################################

test_that("correct single-catch SE", {
    expect_equal(fitsum$predicted[,'SE.estimate'], 
        # c(1.14083246, 0.04887866, 2.73538888), # 1.2.0
        c(1.29664150, 0.06362387, 3.56967573), # 1.2.1
        tolerance = 1e-4, check.attributes = FALSE)
})
###############################################################################

test_that("print.ipsecr no warnings", {
    # see https://stackoverflow.com/questions/22003306/is-there-something-in-testthat-like-expect-no-warnings
    expect_warning(print(fit2), regexp = NA)
})
###############################################################################

test_that("vcov matches expectation", {
    expect_equal(unlist(vcv), 
        # c(1.27890439, 0.00234371, 7.44888511 ),   # 1.2.0
        c(1.644029768, 0.003914797, 12.647485557),  # 1.2.1
        tolerance = 1e-4, check.attributes = FALSE)
})
###############################################################################
