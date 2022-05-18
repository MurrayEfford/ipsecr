## 2022-04-17 start
## 2022-05-08 new proxyfn1

library(ipsecr)
# library(testthat)

## Not needed as RcppParallel not used, but keep as a reminder
## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
## Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

###############################################################################
set.seed(123)
setNumThreads(2)
fit <- ipsecr.fit(captdata, buffer = 100, detectfn = 'HHN', proxyfn = proxyfn1, 
    verbose = FALSE)
pred <- predict(fit)

test_that("correct single-catch estimate", {
    expect_equal(pred[,'estimate'], c(5.6231340, 0.4393444, 28.3127091), 
        tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct single-catch SE", {
    expect_equal(pred[,'SE.estimate'], c(0.5954782, 0.0658304, 1.4052508), 
        tolerance = 1e-4, check.attributes = FALSE)
})
###############################################################################
