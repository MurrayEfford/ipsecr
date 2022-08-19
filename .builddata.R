# .builddata

library(ipsecr)
ipsecrdemo <- ipsecr.fit(captdata, ncores = 1, detectfn = 'HHN', seed = 1237,
    details = list(keep.sim = TRUE))
save(ipsecrdemo, file = "d:/density secr 4.5/ipsecr/data/ipsecrdemo.RData")
tools::resaveRdaFiles("d:/density secr 4.5/ipsecr/data")

