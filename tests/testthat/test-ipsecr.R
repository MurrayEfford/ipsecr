## 2022-04-17 start
## 2022-06-13 new proxy.ms
## 2022-06-21 set RCPP_PARALLEL_BACKEND

library(ipsecr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

###############################################################################
set.seed(123)
setNumThreads(2)
fit <- ipsecr.fit(captdata, buffer = 100, detectfn = 'HHN', proxyfn = proxy.ms, 
    verbose = FALSE)
pred <- predict(fit)

test_that("correct single-catch estimate", {
    expect_equal(pred[,'estimate'], c(5.642231, 0.439468, 28.232423), 
        tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct single-catch SE", {
    expect_equal(pred[,'SE.estimate'], c( 0.660144, 0.062686, 1.362865), 
        tolerance = 1e-4, check.attributes = FALSE)
})
###############################################################################
