# ef 2.0.0 (2026-10-09)

* move likelihood computation to RTMB, the package no longer has compiled code
* random effects can be included using `s(f, bs = "re")` terms in the model
  formula, for example for between visit and between sample variation
* remove `overdispersion()` and `BICadj()`, these are replaced by visit and
  sample random effects
* remove the `sample_re` argument of `efp()`
* smooth terms are now penalised and the amount of smoothing is estimated,
  by fitting smooths as random effects. Previously smooths were unpenalised
  regression splines, which are still available using `fx = TRUE`, for
  example `s(x, k = 3, fx = TRUE)`
* the recommended smooths are first and second order random walks, using
  `s(f, bs = "mrf", xt = list(penalty = rw_penalty(levels(f))))`. The new
  function `rw_penalty()` creates the penalty, and these are fitted using
  sparse matrices so they are quick to fit even with many levels. Other
  smooths from `mgcv` can be used but are slower to fit
* add a `predict()` method for `efp` fits, with predictions for new data,
  standard errors and the option to exclude random effects
* `fitted(type = "pi")` now returns the probability of capture in each pass
* add day of the year (`doy`) to the `ef_data` example dataset, a continuous
  covariate which can be used to fit smooth terms
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
