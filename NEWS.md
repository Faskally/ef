# ef 2.0.0 (2026-10-09)

* move likelihood computation to RTMB, the package no longer has compiled code
* random effects can be included using `s(f, bs = "re")` terms in the model
  formula, for example for between visit and between sample variation
* remove `overdispersion()` and `BICadj()`, these are replaced by visit and
  sample random effects
* remove the `sample_re` argument of `efp()`
* add a `predict()` method for `efp` fits, with predictions for new data,
  standard errors and the option to exclude random effects
* `fitted(type = "pi")` now returns the probability of capture in each pass
* add examples of smooth terms to the help and the README (#1)
* add `summary()`, `vcov()`, `confint()`, `logLik()` and `formula()` methods
  for `efp` fits. `residuals()`, `anova()` and `plot()` now give an
  informative error, rather than failing in the methods inherited from `glm`
* remove the unused `init` argument of `efp()`

# ef 1.2.0 (2020-10-18)

* move likelihood computation to TMB
* introduce random effects at the sample level

# ef 1.1.0 (2016-08-01)

* initial package developed using STAN for likelihood computations
