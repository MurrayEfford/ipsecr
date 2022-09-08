## 2022-08-24 1.2.0 start
## 2022-08-30 explicit RNGkind in set.seed

library(ipsecr)
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

###############################################################################

test_that("proxyfn1 correct", {
    p1 <- proxyfn1(captdata, 'n')
    p2 <- proxyfn1(captdata, 'null')
    p3 <- proxyfn1(captdata, 'zippin')
    p4 <- proxyfn1(captdata, 'jackknife')
    expect_equal(p1, c(4.33073334, -0.03724765, 3.24371994), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(p2, c(4.33163254, -0.03875965, 3.24371994 ),  
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(p3, c(4.336191, 1.478975, 3.243720), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(p4, c(4.5089284, -0.3164427, 3.2437199), 
        tolerance = 1e-4, check.attributes = FALSE)
})
###############################################################################

test_that("proxy.ms correct", {
    tdd <- data.frame(x=traps(captdata)$x)
    m1 <- list(D=~1, lambda0=~1, sigma=~1)
    m2 <- list(D=~x, lambda0=~1, sigma=~1)
    m3 <- list(D=~1, lambda0=~y, sigma=~1)
    m4 <- list(D=~1, lambda0=~1, sigma=~x)
    p1 <- proxy.ms(captdata, model = m1)  # same as proxyfn1()
    p2 <- proxy.ms(captdata, model = m2, trapdesigndata = tdd)
    p3 <- proxy.ms(captdata, model = m3)
    p4 <- proxy.ms(captdata, model = m4)
    expect_equal(p1, c(4.33073334, -0.03724765, 3.24371994), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(p2, c(0.338209883, 0.001024622, -0.037247650, 3.243719942), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(p3, c(4.330733e+00, 4.942692e-01, -2.281577e-05, 3.243720e+00), 
        tolerance = 1e-4, check.attributes = FALSE)
    expect_equal(p4, c(4.330733e+00, -3.724765e-02, 3.159766e+00, 1.023099e-06), 
        tolerance = 1e-4, check.attributes = FALSE)
})
###############################################################################

test_that("proxy.ms spatial sigma OK for Feb96 possum data", {
    skip_if (!requireNamespace("sf"), "spatial sigma test skipped: sf unavailable")
    Feb96 <- OVpossumCH[[1]] # one trapping session (February 1996)
    boundary <- system.file("extdata/OVforest.shp", package = "secr")
    OVforest <- sf::st_read(boundary)
    msk <- make.mask(traps(Feb96), buffer = 120, type = "trapbuffer",
        poly = OVforest[1:2,], spacing = 15, keep.poly = FALSE)
    msk <- suppressWarnings(addCovariates(msk, OVforest[1:2,]))
    m5 <- list(D=~1, lambda0=~1, sigma=~forest)
    p5 <- proxy.ms(Feb96, model = m5, spatialdata = msk)
    expect_equal(p5, c(5.40717177, -0.66005551, 3.08557700, -0.09749582), 
        tolerance = 1e-4, check.attributes = FALSE)
    
})
###############################################################################


test_that("plotProxy means match", {
    set.seed(123)
    trps <- traps(captdata)
    msk <- make.mask(trps, buffer = 100)
    base <- list(D = 5, lambda0 = 0.2, sigma = 25)
    out <- plotProxy (parameter = 'D', traps = trps, mask = msk,
        basepar = base, points = FALSE, boxplot = FALSE, nrepl = 20)
    expect_equal(apply(out,2,mean), 
        # c(3.784703, 3.900069, 4.016776, 4.078801, 4.194038),  # 1.2.0
        c(3.797947196, 3.929653612, 3.998393785, 4.073763702, 4.225303122), # 1.3.0
        tolerance = 1e-4, check.attributes = FALSE)
})
###############################################################################

    