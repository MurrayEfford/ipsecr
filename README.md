<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/ipsecr)](https://cran.r-project.org/package=ipsecr)
[![Downloads (total)](https://cranlogs.r-pkg.org/badges/grand-total/ipsecr)](https://www.r-pkg.org/pkg/ipsecr)
[![R-CMD-check](https://github.com/MurrayEfford/ipsecr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MurrayEfford/ipsecr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
  
# ipsecr
Spatially Explicit Capture-Recapture by Inverse Prediction

This package implements a simulation-based 'inverse prediction' algorithm for fitting spatially explicit capture-recapture models to data from tricky detectors such as single-catch traps (Efford, 2004; Efford, Dawson and Robbins, 2004; Efford, 2023). 

**ipsecr** depends on **secr** 4.5.8 or later. It improves on the functionality of `ip.secr` in earlier versions of **secr**. Many functions mirror those in **secr**.

**ipsecr** is available on CRAN.

The code here is under development. It may be installed using
```
devtools::install_github("MurrayEfford/ipsecr")
```

Compilation of C++ code is required.

Alternatively, the development version may be installed with
```
install.packages("ipsecr", repos = "https://MurrayEfford.r-universe.dev")
```

Please report problems as Issues on GitHub.

Efford, M. G. (2004) Density estimation in live-trapping studies. *Oikos* **106**, 598--610.

Efford, M. G. (2023) ipsecr: An R package for awkward spatial capture-recapture data. 
*Methods in Ecology and Evolution* **14**, 1182--1189.

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software
for analysing capture-recapture data from passive detector arrays.
*Animal Biodiversity and Conservation* **27**, 217--228.
