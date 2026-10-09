# ef 1.3.0 (2026-10-09)

* move likelihood computation to RTMB, the package no longer has compiled code
* random effects can be included using `s(f, bs = "re")` terms in the model
  formula, for example for between visit and between sample variation
* remove `overdispersion()` and `BICadj()`, these are replaced by visit and
  sample random effects
* remove the `sample_re` argument of `efp()`

# ef 1.2.0 (2020-10-18)

* move likelihood computation to TMB
* introduce random effects at the sample level

# ef 1.1.0 (2016-08-01)

* initial package developed using STAN for likelihood computations
