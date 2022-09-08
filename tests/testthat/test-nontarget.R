## 2022-05-13
## 2022-08-30 explicit RNGkind in set.seed
## 2022-09-07 1.3.0

# test simulations of non-target interference

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
    NT = 0.5,
    details = list(nontargettype = 'exclusive'))
summ <- summary(chs)
nt <- as.numeric(summ$nontarget)[1:5]
nd <- as.numeric(summ$counts['detections',1:5])

test_that("single-catch nontarget simulations", {
    expect_equal(nt, 
        # c(25, 24, 19, 15, 17),  # 1.2.0
        c(14, 23, 25, 21, 21),    # 1.3.0
        tolerance = 1e-4, check.attributes = FALSE)
})

test_that("single-catch simulated detections in presence of interference", {
    expect_equal(nd, 
        # c(16, 17, 17, 22, 21),    # 1.2.0
        c(22, 14, 21, 19, 17),    # 1.3.0
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(summ$counts$Total, 
        # c(93, 54, 54, 54, 0, 93, 93, 320),    # 1.2.0
        c(93, 56, 56, 56,  0, 93, 93, 320),     # 1.3.0
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
    expect_equal(nt, c(17, 25, 31, 25, 25), 
        tolerance = 1e-4, check.attributes = FALSE)
})
test_that("multi-catch simulated detections in presence of interference", {
    expect_equal(nd, c(27, 19, 31, 28, 23), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(summ$counts$Total, c(128, 65, 65, 65, 0, 128, 107, 320), 
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
    expect_equal(nt, c(17, 25, 31, 25, 25), 
        tolerance = 1e-4, check.attributes = FALSE)
})
test_that("proximity simulated detections in presence of interference", {
    expect_equal(nd, c(34, 20, 35, 32, 29), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(summ$counts$Total, c(128, 65, 65, 65, 0, 150, 119, 320), 
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
    expect_equal(nt, c(14, 23, 25, 21, 20), 
        tolerance = 1e-4, check.attributes = FALSE)
})
test_that("capped simulated detections in presence of interference", {
    expect_equal(nd, c(26, 15, 21, 20, 20), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(summ$counts$Total, c(91, 56, 56, 56, 0, 102, 102, 320), 
        tolerance = 1e-4, check.attributes = FALSE)
})
#------------------------------------------------------------------

set.seed(123)
trC <- make.grid(8, 8, spacing = 30, detector = 'count')
pop <- sim.popn(10, trC, 100)
chc <- simCH(
    traps = trC, 
    popn = pop, 
    detectfn = 1, 
    detparmat = list(lambda0 = 0.2, sigma = 20, z = 2),
    noccasions = 5, 
    NT = 0.5,
    details = list(nontargettype = 'dependent'))

summ <- summary(chc)
nt <- as.numeric(summ$nontarget)[1:5]
nd <- as.numeric(summ$counts['detections',1:5])
chv <- function(v) paste0('c(', paste(v, collapse=', '), ')')

# chv(nt)
# chv(nd)
# chv(summ$counts$Total)

test_that("count nontarget simulations", {
    expect_equal(nt, c(36, 50, 37, 42, 34), 
        tolerance = 1e-4, check.attributes = FALSE)
})
test_that("count simulated detections in presence of interference", {
    expect_equal(nd, c(42, 47, 43, 46, 47), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(summ$counts$Total, c(184, 101, 101, 101, 0, 225, 164, 320), 
        tolerance = 1e-4, check.attributes = FALSE)
})
#------------------------------------------------------------------


