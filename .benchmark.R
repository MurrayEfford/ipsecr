# 2022-08-25

library(rbenchmark)

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

# 2022-09-02 ipsecr 1.2.1
# test replications elapsed relative user.self sys.self user.child sys.child
# 1 ncores1           10 1245.25    1.340   1220.50     7.34         NA        NA
# 2 ncores2           10  929.62    1.000      4.54     2.28         NA        NA
# 3 ncores7           10  999.08    1.075     10.20     8.33         NA        NA

# 2022-09-07 ipsecr 1.3.0

#      test replications elapsed relative user.self sys.self user.child sys.child
# 1 ncores1           10 1112.97    1.292   1236.41    12.84         NA        NA
# 2 ncores2           10  861.53    1.000      4.36     2.35         NA        NA
# 3 ncores7           10 1030.96    1.197     10.05     8.40         NA        NA