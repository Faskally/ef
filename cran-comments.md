## Update

This is an update of ef from version 1.2.0 to 2.0.0.

* The likelihood is now computed using 'RTMB', so the package no longer
  contains compiled code and no longer links to 'TMB' and 'RcppEigen'.
* `overdispersion()` and `BICadj()` have been removed, along with the
  `sample_re` and `init` arguments of `efp()`. They are replaced by random
  effects specified in the model formula, see NEWS.md.

## Test environments

* local: Ubuntu 24.04, R 4.5.3

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.
