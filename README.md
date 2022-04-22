# ipsecr
Spatially Explicit Capture-Recapture by Inverse Prediction

This package implements a simulation-based 'inverse prediction' algorithm for fitting spatially explicit capture-recapture models to data from tricky detectors such as single-catch traps (Efford, 2004; Efford, Dawson and Robbins, 2004). 

**ipsecr** depends on **secr** 4.5.4 or later. It improves on the functionality of `ip.secr` in earlier versions of **secr**. Many functions mirror those in **secr**.

The code here is under development. It may be installed using
```
devtools::install_github("MurrayEfford/ipsecr")
```

Compilation of C++ code is required.

Please report problems as Issues on GitHub.

Efford, M. G. (2004) Density estimation in live-trapping studies. *Oikos* **106**, 598--610.

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software
for analysing capture-recapture data from passive detector arrays.
*Animal Biodiversity and Conservation* **27**, 217--228.