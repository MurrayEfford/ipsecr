# 2022-05-13

# test simulations of non-target interference

library(ipsecr)

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
    NT = 0.5,
    details = list(nontargettype='exclusive'))
summ <- summary(chs)
nt <- as.numeric(summ$nontarget)[1:5]
nd <- as.numeric(summ$counts['detections',1:5])

test_that("single-catch nontarget simulations", {
    expect_equal(nt, c(25, 24, 19, 15, 17), 
        tolerance = 1e-4, check.attributes = FALSE)
})
test_that("single-catch simulated detections in presence of interference", {
    expect_equal(nd, c(16, 17, 17, 22, 21), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(summ$counts$Total, c(93, 54, 54, 54, 0, 93, 93, 320), 
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
    NT = 0.5, 
    details = list(nontargettype='truncated'))
summ <- summary(chm)
nt <- as.numeric(summ$nontarget)[1:5]
nd <- as.numeric(summ$counts['detections',1:5])

test_that("multi-catch simulated frequency of interference", {
    expect_equal(nt, c(23, 28, 26, 18, 29), 
        tolerance = 1e-4, check.attributes = FALSE)
})
test_that("multi-catch simulated detections in presence of interference", {
    expect_equal(nd, c(16, 26, 24, 19, 21), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(summ$counts$Total, c(106, 64, 64, 64, 0, 106, 84, 320), 
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
    NT = 0.5,
    details = list(nontargettype = 'truncated'))
summ <- summary(chp)
nt <- as.numeric(summ$nontarget)[1:5]
nd <- as.numeric(summ$counts['detections',1:5])

test_that("proximity simulated frequency of interference", {
    expect_equal(nt, c(23, 28, 26, 18, 29), 
        tolerance = 1e-4, check.attributes = FALSE)
})
test_that("proximity simulated detections in presence of interference", {
    expect_equal(nd, c(15, 24, 28, 26, 15), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(summ$counts$Total, c(98, 64, 64, 64, 0, 108, 92, 320), 
        tolerance = 1e-4, check.attributes = FALSE)
})
#------------------------------------------------------------------

set.seed(123)
trc <- make.grid(8, 8, spacing = 30, detector = 'capped')
pop <- sim.popn(10, trc, 100)
chc <- simCH(
    traps = trc, 
    popn = pop, 
    detectfn = 14, 
    detparmat = list(lambda0 = 0.2, sigma = 20),
    noccasions = 5, 
    NT = 0.5,
    details = list(nontargettype = 'exclusive'))
summ <- summary(chc)
nt <- as.numeric(summ$nontarget)[1:5]
nd <- as.numeric(summ$counts['detections',1:5])

test_that("capped nontarget simulations", {
    expect_equal(nt, c(22, 18, 22, 14, 24), 
        tolerance = 1e-4, check.attributes = FALSE)
})
test_that("capped simulated detections in presence of interference", {
    expect_equal(nd, c(13, 21, 23, 22, 13), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(summ$counts$Total, c(84, 58, 58, 58,  0, 92, 92, 320), 
        tolerance = 1e-4, check.attributes = FALSE)
})
#------------------------------------------------------------------
