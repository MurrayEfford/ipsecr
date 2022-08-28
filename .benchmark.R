# 2022-08-25

library(rbenchmark)

# in progress

library(ipsecr)
set.seed(1235)
benchmark(
    ncores1 = ipsecr.fit(captdata, detectfn = 'HHN', ncores = 1, verbose = FALSE),
    ncores2 = ipsecr.fit(captdata, detectfn = 'HHN', ncores = 2, verbose = FALSE),
    ncores7 = ipsecr.fit(captdata, detectfn = 'HHN', ncores = 7, verbose = FALSE),
    replications = 10)

#      test replications elapsed relative user.self sys.self user.child sys.child
# 1 ncores1           10 1377.41    1.444   1375.93     1.02         NA        NA
# 2 ncores2           10  953.68    1.000      3.83     0.82         NA        NA
# 3 ncores7           10 1037.17    1.088      9.41     2.87         NA        NA
