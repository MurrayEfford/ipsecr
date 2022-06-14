
onecombo <- function (D, lambda0, sigma, nrepl = 50) {
    out <- matrix(nrow = nrepl, ncol=3)
    for (r in 1:nrepl) {
        pop <- sim.popn(D = D, core = tr, buffer = 120)
        ch <- sim.capthist(traps = tr, popn = pop, detectfn = 'HHN', 
            detectpar = list(lambda0=lambda0, sigma=sigma), noccasions = 5)
        out[r,] <- proxyfn1(ch)
    }
    out
}

library(ipsecr)
tr <- traps(captdata)
Dval <- seq(3,7,0.5)
lambda0val <- seq(0.1,0.3,0.025)
sigmaval <- seq(20,40,2.5)

Dlist <- lapply(Dval, onecombo, lambda0 = 0.2, sigma = 30)
lamlist <- lapply(lambda0val, onecombo, D = 5, sigma = 30)
siglist <- lapply(sigmaval, onecombo, D = 5, lambda0 = 0.2)
Dproxy <- sapply(Dlist, '[',,1)
lamproxy <- sapply(lamlist, '[',,2)
sigproxy <- sapply(siglist, '[',,3)

png('d:/density secr 4.5/ipsecr/vignettes/linearityplot.png', width=900, height=320)

par(mfrow=c(1,3), mgp=c(2.4,0.7,0), mar=c(4,4,2,2), cex=1.2)
plot(1,1,type='n', xlim=c(2.8,8), ylim=c(3,4.8),
    xlab = 'Density / ha', ylab = 'log(n)', log='x')
apply(Dproxy,1,points, x=Dval)
points(Dval, apply(Dproxy,2,mean), pch = 16, type='b', col='red', lwd=1.5)

plot(1,1,type='n', xlim=c(0.09,0.35), ylim=c(-1,0.1),
    xlab = 'lambda0', ylab = 'cloglog(p)', log='x')
apply(lamproxy, 1 ,points, x=lambda0val)
points(lambda0val, apply(lamproxy,2,mean), pch = 16, type='b', col='red', lwd=1.5)

plot(1,1,type='n', xlim = c(18,46), ylim=c(2.5, 3.7),
    xlab = 'sigma m', ylab = 'log(RPSV)', log='x')
apply(sigproxy,1,points, x=sigmaval)
points(sigmaval, apply(sigproxy,2,mean), pch = 16, type='b', col='red', lwd=1.5)
dev.off()
