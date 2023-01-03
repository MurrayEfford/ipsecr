## 2023-01-01 1.4.0

# test simulations of without non-target interference

library(ipsecr)
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

#------------------------------------------------------------------
set.seed(123)
trs <- make.grid(8, 8, spacing = 30, detector = 'single')
pop <- sim.popn(10, trs, 100)
detparmat <- matrix(c(0.2,20), byrow = TRUE, nrow = nrow(pop), ncol = 2)

chs <- simCH(
    traps = trs, 
    popn = pop, 
    detectfn = 14, 
    detparmat = detparmat,
    noccasions = 5, 
    NT = 0)

test_that("RPSV of single-catch simulations, detectfn HHN", {
    # expect_equal(rpsv(chs), 18.853879,    # 1.3.0
    expect_equal(rpsv(chs), 16.2889591,    # 1.4.0
        tolerance = 1e-4, check.attributes = FALSE)
})

#------------------------------------------------------------------

set.seed(123)
trm <- make.grid(8, 8, spacing = 30, detector = 'multi')
pop <- sim.popn(10, trm, 100)
detparmat <- matrix(c(0.2,20), byrow = TRUE, nrow = nrow(pop), ncol = 2)
chm <- simCH(
    traps = trm, 
    popn = pop, 
    detectfn = 14, 
    detparmat = detparmat,
    noccasions = 5, 
    NT = 0)

test_that("RPSV of multi-catch simulations", {
    # expect_equal(rpsv(chm), 19.348499,    # 1.3.0
    expect_equal(rpsv(chm), 16.81577445,    # 1.4.0
        tolerance = 1e-4, check.attributes = FALSE)
})
#------------------------------------------------------------------

set.seed(123)
trp <- make.grid(8, 8, spacing = 30, detector = 'proximity')
pop <- sim.popn(10, trp, 100)
chp <- simCH(
    traps = trp, 
    popn = pop, 
    detectfn = 14, 
    detparmat = list(lambda0 = 0.2, sigma = 20),
    noccasions = 5, 
    NT = 0)

test_that("RPSV of proximity simulations", {
    # expect_equal(rpsv(chp), 20.269715,    # 1.3.0
    expect_equal(rpsv(chp), 17.56461286,    # 1.4.0
        tolerance = 1e-4, check.attributes = FALSE)
})
#------------------------------------------------------------------

# count detectors HN detectfn

set.seed(123)
trC <- make.grid(8, 8, spacing = 30, detector = 'count')
pop <- sim.popn(10, trC, 100)
chC <- simCH(
    traps = trC, 
    popn = pop, 
    detectfn = 0, 
    detparmat = list(lambda0 = 0.2, sigma = 20),
    noccasions = 5, 
    NT = 0)

test_that("RPSV of count simulations detectfn HN", {
    # expect_equal(rpsv(chC),  18.234897,    # 1.3.0
    expect_equal(rpsv(chC),  17.84443647,    # 1.4.0
        tolerance = 1e-4, check.attributes = FALSE)
})
#------------------------------------------------------------------


