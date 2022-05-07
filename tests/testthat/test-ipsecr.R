## Started 2022-04-17,19,20

library(ipsecr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
## Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

# create small working datasets

###############################################################################
set.seed(123)
setNumThreads(2)
fit <- ipsecr.fit(captdata, buffer = 100, detectfn = 'HHN', verbose = FALSE)
pred <- predict(fit)

test_that("correct single-catch estimate", {
    expect_equal(pred[,'estimate'], c(5.6772568, 0.4365298, 28.2901364), 
        tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct single-catch SE", {
    expect_equal(pred[,'SE.estimate'], c(0.81674255, 0.06611273, 1.40548990), 
        tolerance = 1e-4, check.attributes = FALSE)
})

###############################################################################
