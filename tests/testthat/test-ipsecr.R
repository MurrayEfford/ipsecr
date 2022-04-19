## Started 2022-04-17,19

library(ipsecr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
# Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

# create small working datasets

###############################################################################
set.seed(1237)
setNumThreads(2)
fit <- ipsecr.fit(captdata, buffer = 100, detectfn = 'HHN')
pred <- predict(fit)

test_that("correct single-catch estimate", {
    expect_equal(pred[,'estimate'], c(5.6225459,  0.4374263, 28.2450044), 
        tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct single-catch SE", {
    expect_equal(pred[,'SE.estimate'], c(0.65039471, 0.07577939, 1.33063897), 
        tolerance = 1e-4, check.attributes = FALSE)
})

###############################################################################
