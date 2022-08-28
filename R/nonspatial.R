## closed population size likelihoods
negloglikM0   <- function (theta, n) {
    NN   <- exp(theta)
    ndot <- n[1]  # total captures
    Mt1  <- n[2]  # total animals
    K    <- n[3]  # number of occasions
    LL   <- lgamma (NN+1) - lgamma(NN-Mt1+1) + ndot * log(ndot) +
        (K*NN - ndot) * log(K*NN - ndot) - K*NN*log(K*NN)
    if (is.finite(LL)) -LL
    else 1e10
}

negloglikMbold   <- function (theta, u) {
    p <- invlogit(theta)
    K <- length(u)
    n  <- sum(u)
    -(n * log(p)+ sum(u *(0:(K-1))) * log (1-p) - n * log(1 - (1-p)^K))
}

#     tt <- length(ut)
#     m <- sum(nt-ut)
#     Mj <- c(0,cumsum(ut))
#     M.  <- sum(Mj[2:tt])
#     Mt1 <- Mj[tt+1]
    
#     ## OK <- sum( (tt+1-2*(1:tt)) * ut) > 0
#     fit <- optimize(f = loglik, interval = c(Mt1, maxN),  maximum = TRUE)
#     nhat <- fit$maximum
#     phat <- Mt1 / (tt*nhat - M.)
#     senhat <- sqrt((nhat * (1-phat)^tt * (1 - (1-phat)^tt)) /
#             ((1 - (1-phat)^tt)^2 - tt^2*phat^2*(1-phat)^(tt-1)))
#     c(Mt1 = Mt1, Nhat = nhat, seNhat = senhat, npar = 3, LL = fit$objective)
# }

##################################################

jack.est <- function (inp, deads = 0, full = F)
{
    
    # Calculate Burnham & Overton's jackknife estimate for closed populations
    #
    # inp may be a vector of capture frequencies, or
    #            a list comprising such a vector as its first element and
    #                              the scalar number of 'deads' as its second
    #
    # 'deads' are assumed to be additional to the tabulated capture frequencies
    # They are added to the calculated population size (cf Otis et al. 1978)
    #
    # MODIFIED
    # 30/3/95 Fix bug when mt=5 (last jacknife selected) - length(test) s/b 5
    # 30/3/95 Optional full output
    # 30/3/95 Add confIDence limits and deads to 'short' output
    # 3/4/95  Implement unconditional se using method of K.P.Burnham
    # 3/4/95  Minor changes to full output
    #
    
    first	<- function(vec) match(1, vec)
    
    jack.fill <- function(tt)
    {
        
        # Input:  tt is the number of capture occasions (e.g., days)
        # Output: matrix of jackknife coefficients for Burnham & Overton estimator
        
        # Murray Efford 30/3/95
        
        T1 <- tt - 1
        T2 <- tt - 2
        T3 <- tt - 3
        T4 <- tt - 4
        T5 <- tt - 5
        occ5 <- min(tt, 5)
        fcoeff <- matrix(data = 0, ncol = 5, nrow = 5)
        fcoeff[1, 1] <- T1/tt
        fcoeff[2, 1] <- (2 * tt - 3)/tt
        fcoeff[2, 2] <-  - T2^2/(tt * T1)
        fcoeff[3, 1] <- (3 * tt - 6)/tt
        fcoeff[3, 2] <-  - (3 * tt^2 - 15 * tt + 19)/(tt * T1)
        fcoeff[3, 3] <- T3^3/(tt * T1 * T2)
        fcoeff[4, 1] <- (4 * tt - 10)/tt
        fcoeff[4, 2] <-  - (6 * tt^2 - 36 * tt + 55)/(tt * T1)
        fcoeff[4, 3] <- (4 * tt^3 - 42 * tt^2 + 148 * tt - 175)/(tt * T1 * T2)
        fcoeff[4, 4] <-  - T4^4/(tt * T1 * T2 * T3)
        fcoeff[5, 1] <- (5 * tt - 15)/tt
        fcoeff[5, 2] <-  - (10 * tt^2 - 70 * tt + 125)/(tt * T1)
        fcoeff[5, 3] <- (10 * tt^3 - 120 * tt^2 + 485 * tt - 660)/(tt * T1 * T2)
        fcoeff[5, 4] <-  - (T4^5 - T5^5)/(tt * T1 * T2 * T3)
        fcoeff[5, 5] <- T5^5/(tt * T1 * T2 * T3 * T4)
        fcoeff <- fcoeff[1:occ5, 1:occ5]	# Use sub-matrix if tt < 5
        if(tt > 5) {
            fcoeff <- cbind(fcoeff, matrix(data = 0, nrow = 5, ncol = tt - 5))
        }
        fcoeff + 1	# Adding one allows: Nj<-fcoeff%*%fi
    }
    
    if(is.list(inp)) {
        fi 	<- inp[[1]]
        deads	<- inp[[2]]
    }
    else fi	<- inp
    
    S	<- sum(fi)
    occ	<- length(fi)
    occ5	<- min(5, occ)
    aki	<- jack.fill(occ)
    nj		<- aki %*% fi
    varnj 	<- (aki^2 %*% fi) - nj
    difnk	<- nj[2:occ5] - nj[1:(occ5 - 1)]
    b2f	<- ((aki[2:occ5,  ] - aki[1:(occ5 - 1),  ])^2) %*% fi
    test	<- (difnk/sqrt(S/(S - 1) * (b2f - difnk^2/S)))^2
    test[is.na(test)] <- 0
    test	<- rbind(test, 0)
    pk		<- 1 - pchisq((difnk/sqrt(S/(S - 1) * (b2f - difnk^2/S)))^2, 1)
    pk[difnk == 0] <- 1
    
    # Select jacknife on basis of test results, and interpolate estimate
    # The conditional variance calculations commented out here are superceded by the
    # unconditional calculations below
    
    mt		<- first(test < 3.84) #
    if(mt == 1) {
        xtest	<- (nj[1] * fi[1])/(nj[1] - fi[1])
        if(xtest > 3.84) {
            alpha	<- (xtest - 3.84)/(xtest - test[1])
            beta	<- 1 - alpha
            z	<- aki[1,  ] * alpha + beta
            N <- z %*% fi  # varN <- (z * z) %*% fi - N
        } 
        else {
            N	<- nj[1]
            varN	<- varnj[1]
        }
    } 
    else {
        alpha	<- (test[mt - 1] - 3.84)/(test[mt - 1] - test[mt])
        beta	<- 1 - alpha
        z 	<- aki[mt,  ] * alpha + aki[mt - 1,  ] * beta
        N 	<- z %*% fi  # varN <- ((z * z) %*% fi) - N
    }
    k1		<- occ5 - 1
    Pi		<- rep(1, occ5)
    Beta 	<- pchisq(rep(3.8415, occ5), 1, pmax(test - 1, 0))
    for(i in 2:occ5) Pi[i] <- Pi[i - 1] * (1 - Beta[i - 1])
    Pi[1:k1]	<- Pi[1:k1] * Beta[1:k1]
    if(occ5 < 5) Pi[occ5] <- 1 - sum(Pi[1:k1])
    varN 	<- sum(Pi * varnj) + sum(Pi * (nj - sum(Pi * nj))^2)
    sejack <- sqrt(varN)
    
    # Confidence limits as in Rexstad & Burnham 1991 CAPTURE Users' Guide
    
    f0		<- N - S
    if (f0>0)
        cc1 	<- exp(1.96 * sqrt(log(1 + varN/f0^2)))
    else
        cc1    <- NA
    jacklcl	<- S + f0 / cc1
    jackucl	<- S + f0 * cc1
    
    # Output
    if(full) list(fi = fi, jack = N + deads, deads = deads, sejack = sejack, jacklcl
        = jacklcl + deads, jackucl = jackucl + deads, mt = mt, aki = aki,
        nj = data.frame(nj, se = sqrt(varnj), X2 = test, Pi = Pi, Beta = Beta))
    else c(N + deads, sqrt(varN), jacklcl + deads, jackucl + deads, deads)
}
##################################################

M0 <- function (ni) {
    ## data vector ni is c(total captures, total animals, n occasions)
    Mt1 <- ni[2]
    if (Mt1==0) {
        c(0,NA,NA)
    } 
    else {
        theta <- optimize (f = negloglikM0, lower = log(Mt1), upper = log(1000*Mt1),
            n = ni)$minimum
        c(exp(theta), ni[1]/exp(theta)/ni[3])
    }
}
##################################################

Mbold <- function (ui) {
    ## data vector ui is number of new animals on each occasion
    n <- sum(ui)
    if (n==0) {
        c(0,NA)
    } 
    else {
        theta <- optimize (f = negloglikMb, lower = logit(0.0001), upper = logit(0.9999),
            u = ui)$minimum
        p <- invlogit(theta)
        p <- 1 - (1-p)^length(ui)
        c(n/p, p)  ## Nhat, p
    }
}
##################################################


negloglikMb   <- function (theta, tt, m, M., Mt1) {
    N <- exp(theta)
    cterm <- ifelse (m>0, m*log(m/M.) + (M.-m) * log(1-m/M.), 0)
    lgamma(N+1) - lgamma(N-Mt1+1) - lgamma(Mt1+1) +
        Mt1*log(Mt1) +
        (tt*N - M. - Mt1) * log(tt*N - M. - Mt1) -
        (tt*N - M.) * log(tt*N - M.) +
        cterm
}

Mb <- function (ni, ui) {
    ## data vector ni is number of new animals on each occasion
    ## data vector ui is number of new animals on each occasion
    n <- sum(ui)
    if (n==0) {
        c(0,NA)
    } 
    else {
        tt <- length(ui)
        m <- sum(ni-ui)
        Mj <- c(0,cumsum(ui))
        M.  <- sum(Mj[2:tt])
        Mt1 <- Mj[tt+1]
        theta <- optimize (f = negloglikMb, lower = log(Mt1), upper = log(1000*Mt1),
            tt = tt, m = m, M. = M., Mt1 = Mt1, maximum = TRUE)$maximum
        Nhat <- exp(theta)
        phat <- Mt1 / (tt * Nhat - M.)
        c(Nhat, 1 - (1-phat)^tt)
    }
}
##################################################

Mh <- function (fi) {
    ## data vector fi is number of animals caught on exactly i occasions
    temp <- jack.est(fi)
    nocc <- length(fi)
    p <- sum(fi*(1:nocc)) / temp[1] / nocc
    c(temp[1], p)
}
##################################################
