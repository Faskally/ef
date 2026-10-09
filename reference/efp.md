# Estimate capture probabilities from electrofishing data

This function uses the marginal likelihood of capture probabilities to
estimate model parameters

## Usage

``` r
efp(
  formula,
  data = NULL,
  pass = pass,
  id = id,
  offset = NULL,
  verbose = FALSE,
  hessian = TRUE,
  fit = TRUE,
  control = list(),
  sparse = TRUE
)
```

## Arguments

- formula:

  a formula object. Random effects are specified using `s(f, bs = "re")`
  terms, where `f` is a factor, and smooths using other `s()` terms, see
  Details.

- data:

  a data.frame containing all relevant info

- pass:

  a vector of integers giving the pass number of the observation

- id:

  a vector of integers identifying an observation (a set of
  electrofishing passes)

- offset:

  an optional offset for the linear predictor of capture probability

- verbose:

  if TRUE optimiser messages are printed to the screen

- hessian:

  if TRUE the hessian is computed and the covariance matrix of the
  parameters is returned via Vb

- fit:

  if TRUE model is fitted if FALSE the data that would be passed to the
  optimiser is returned

- control:

  a list of control commands to be passed to
  [`optim`](https://rdrr.io/r/stats/optim.html)

- sparse:

  if TRUE (the default) penalised smooths are set up so that sparse
  bases and penalties stay sparse, see Details. If FALSE the constrained
  basis and penalties from `mgcv` are used, which are dense and so
  slower to fit.

## Value

If `fit = TRUE`, an object of class `efp`, a list which inherits from
`glm`, containing among other things the fixed effect estimates
(`coefficients`), the random effect estimates (`random`), the standard
deviations of the random effects and smooths (`sigma`), the fitted
capture probabilities (`fitted`), the log likelihood (`llik`) and the
variance matrix of the fixed effects (`Vb`). See
[`efp-methods`](https://faskally.github.io/ef/reference/efp-methods.md)
for the methods available. If `fit = FALSE`, a list containing the data
that would be passed to the optimiser.

## Details

Random effects on the logit capture probability are included by adding
`s(f, bs = "re")` terms to the formula, where `f` is a factor in `data`.
For example, `s(visitID, bs = "re")` adds between visit variation and
`s(sampleID, bs = "re")` adds between sample variation. Each term has
its own standard deviation, and the random effects are integrated out of
the likelihood using the Laplace approximation, so that `logLik` and
`AIC` are based on the marginal likelihood. Fitted values are
conditional on the estimated random effects.

Smooth terms can also be included in the formula. These are penalised,
and the amount of smoothing is estimated: a smooth is fitted as a random
effect, in which the part of the smooth that is not penalised (for
example a linear trend) is a fixed effect, and the standard deviation of
the random part controls how smooth the term is. A standard deviation
close to zero means that only the unpenalised part of the smooth is
needed. As for other random effects, `logLik` and `AIC` are based on the
marginal likelihood.

The recommended smooths are first and second order random walks, using
the Markov random field basis of `mgcv` with a penalty from
[`rw_penalty`](https://faskally.github.io/ef/reference/rw_penalty.md):
`s(f, bs = "mrf", xt = list(penalty = rw_penalty(levels(f))))`, where
`f` is a factor whose levels are the steps of the walk. The basis and
the penalty are sparse, and are kept sparse when fitting, so these
models are quick to fit even when there are many levels. The standard
deviation is that of the steps of the random walk, and the starting
values of the walk (the first level for a first order walk, the first
two for a second order walk) are fixed effects: this means that the
intercept is the logit capture probability at the first level of the
walk, and that a second order walk has an additional fixed effect for
its linear trend, which is on an arbitrary scale.

Any other penalised smooth from `mgcv` can be used, such as `s(x)` or
`te(x, z)`, but their penalties are dense so they will take longer to
fit, particularly if the basis dimension `k` is large. To fit an
unpenalised regression spline, whose flexibility is set entirely by `k`,
use `fx = TRUE`, for example `s(x, k = 3, fx = TRUE)`.

The fixed effects are estimated by maximum likelihood, so the likelihood
depends slightly on how the unpenalised part of a smooth is separated
from the penalised part. With `sparse = TRUE` one coefficient of the
smooth is fixed for each function that is not penalised, as described
above for random walks. With `sparse = FALSE`, and for smooths with `by`
variables or models without an intercept, the identifiability
constraints of `mgcv` are applied to the smooth first, which makes the
basis and the penalty dense. The two approaches give very similar, but
not identical, fits.

## References

Millar, C. P., Fryer, R. J., Millidine, K. J. and Malcolm, I. A. (2016)
Modelling capture probability of Atlantic salmon (Salmo salar) from a
diverse national electrofishing dataset: Implications for the estimation
of abundance. Fisheries Research, 177, 1-12.
[doi:10.1016/j.fishres.2016.01.001](https://doi.org/10.1016/j.fishres.2016.01.001)

## Examples

``` r
# create two electrofishing site visits with 3 and 4 passes and 2 lifestages
ef_data <- data.frame(n      = c(100, 53, 24, 50, 26, 12,
                                 100, 53, 24, 50, 26, 12),
                      pass   = c(  1,  2,  3,  1,  2,  3,
                                   1,  2,  3,  1,  2,  3),
                      stage  = c(  1,  1,  1,  2,  2,  2,
                                   1,  1,  1,  2,  2,  2),
                      sample = c(  1,  1,  1,  2,  2,  2,
                                   3,  3,  3,  4,  4,  4))

ef_data2 <- data.frame(n      = c(100, 53, 24, 50, 26, 12,
                                  100, 53, 24, 12, 50, 26, 12, 6),
                       pass   = c(  1,  2,  3,  1,  2,  3,
                                    1,  2,  3,  4,  1,  2,  3, 4),
                       stage  = c(  1,  1,  1,  2,  2,  2,
                                    1,  1,  1,  1,  2,  2,  2, 2),
                       sample = c(  1,  1,  1,  2,  2,  2,
                                    3,  3,  3,  3,  4,  4,  4, 4))

ef_data3 <- data.frame(n      = c(100, 53, 24, 50, 26, 12, 40,
                                  100, 53, 24, 12, 50, 26, 12, 6, 40),
                       pass   = c(  1,  2,  3,  1,  2,  3, 1,
                                    1,  2,  3,  4,  1,  2,  3, 4, 1),
                       stage  = c(  1,  1,  1,  2,  2,  2, 1,
                                    1,  1,  1,  1,  2,  2,  2, 2, 2),
                       sample = c(  1,  1,  1,  2,  2,  2, 5,
                                    3,  3,  3,  3,  4,  4,  4, 4, 6))

# Fit a simple model
m2 <- efp(n ~ 1 + factor(stage), data = ef_data, pass = pass, id = sample)
cbind(ef_data, fit = fitted(m2))
#>      n pass stage sample       fit
#> 1  100    1     1      1 0.5007600
#> 2   53    2     1      1 0.5007600
#> 3   24    3     1      1 0.5007600
#> 4   50    1     2      2 0.5030559
#> 5   26    2     2      2 0.5030559
#> 6   12    3     2      2 0.5030559
#> 7  100    1     1      3 0.5007600
#> 8   53    2     1      3 0.5007600
#> 9   24    3     1      3 0.5007600
#> 10  50    1     2      4 0.5030559
#> 11  26    2     2      4 0.5030559
#> 12  12    3     2      4 0.5030559
m3 <- efp(n ~ 1 + factor(stage), data = ef_data2, pass = pass, id = sample)
cbind(ef_data2, fit = fitted(m3))
#>      n pass stage sample       fit
#> 1  100    1     1      1 0.5033905
#> 2   53    2     1      1 0.5033905
#> 3   24    3     1      1 0.5033905
#> 4   50    1     2      2 0.5047679
#> 5   26    2     2      2 0.5047679
#> 6   12    3     2      2 0.5047679
#> 7  100    1     1      3 0.5033905
#> 8   53    2     1      3 0.5033905
#> 9   24    3     1      3 0.5033905
#> 10  12    4     1      3 0.5033905
#> 11  50    1     2      4 0.5047679
#> 12  26    2     2      4 0.5047679
#> 13  12    3     2      4 0.5047679
#> 14   6    4     2      4 0.5047679
m4 <- efp(n ~ 1 + factor(stage), data = ef_data3, pass = pass, id = sample)
cbind(ef_data3, fit = fitted(m4))
#>      n pass stage sample       fit
#> 1  100    1     1      1 0.5033905
#> 2   53    2     1      1 0.5033905
#> 3   24    3     1      1 0.5033905
#> 4   50    1     2      2 0.5047679
#> 5   26    2     2      2 0.5047679
#> 6   12    3     2      2 0.5047679
#> 7   40    1     1      5 0.5033905
#> 8  100    1     1      3 0.5033905
#> 9   53    2     1      3 0.5033905
#> 10  24    3     1      3 0.5033905
#> 11  12    4     1      3 0.5033905
#> 12  50    1     2      4 0.5047679
#> 13  26    2     2      4 0.5047679
#> 14  12    3     2      4 0.5047679
#> 15   6    4     2      4 0.5047679
#> 16  40    1     2      6 0.5047679

# create two electrofishing site visits with 3 and 4 passes and 2 lifestages
ef_data <- data.frame(n      = c(200, 53, 24, 100, 26, 12,
                                 200, 53, 24, 100, 26, 12),
                      pass   = c(  1,  2,  3,  1,  2,  3,
                                   1,  2,  3,  1,  2,  3),
                      stage  = c(  1,  1,  1,  2,  2,  2,
                                   1,  1,  1,  2,  2,  2),
                      sample = c(  1,  1,  1,  2,  2,  2,
                                   3,  3,  3,  4,  4,  4))
# Fit a simple model
m2 <- efp(n ~ 1 + factor(stage) + factor(replace(pass, pass> 2, 2)),
          data = ef_data, pass = pass, id = sample)
out <- cbind(ef_data, p = fitted(m2, type = "p"))
out
#>      n pass stage sample         p
#> 1  200    1     1      1 0.6731839
#> 2   53    2     1      1 0.5434508
#> 3   24    3     1      1 0.5434508
#> 4  100    1     2      2 0.6754774
#> 5   26    2     2      2 0.5460407
#> 6   12    3     2      2 0.5460407
#> 7  200    1     1      3 0.6731839
#> 8   53    2     1      3 0.5434508
#> 9   24    3     1      3 0.5434508
#> 10 100    1     2      4 0.6754774
#> 11  26    2     2      4 0.5460407
#> 12  12    3     2      4 0.5460407

# a smooth effect of a continuous covariate, day of the year, using a
# second order random walk
data(ef_data)
ef_data$sampleID <- factor(paste(ef_data$date, ef_data$siteID,
                                 ef_data$lifestage, ef_data$species))
#> Error in `$<-.data.frame`(`*tmp*`, sampleID, value = structure(integer(0), levels = character(0), class = "factor")): replacement has 0 rows, data has 12
ef_data$pass12 <- factor(replace(ef_data$pass, ef_data$pass > 2, 2))
days <- min(ef_data$doy):max(ef_data$doy)
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Error in min(ef_data$doy):max(ef_data$doy): result would be too long a vector
ef_data$fdoy <- factor(ef_data$doy, levels = days)
#> Error: object 'days' not found
rw2 <- rw_penalty(days, order = 2)
#> Error: object 'days' not found
m5 <- efp(count ~ siteID + pass12 + lifestage +
            s(fdoy, bs = "mrf", xt = list(penalty = rw2)),
          data = ef_data, pass = pass, id = sampleID)
#> Error in eval(id, data, environment(formula)): object 'sampleID' not found
m5
#> Error: object 'm5' not found

# predicted capture probability over the smooth, with standard errors
newdata <- data.frame(siteID = "Site1", pass12 = "1", lifestage = "Fry",
                      fdoy = factor(days, levels = days))
#> Error: object 'days' not found
pred <- predict(m5, newdata = newdata, type = "response", se.fit = TRUE)
#> Error: object 'm5' not found
cbind(newdata, fit = pred$fit, se = pred$se.fit)
#> Error: object 'newdata' not found

# an unpenalised regression spline
m6 <- efp(count ~ siteID + pass12 + lifestage + s(doy, k = 3, fx = TRUE),
          data = ef_data, pass = pass, id = sampleID)
#> Error in eval(id, data, environment(formula)): object 'sampleID' not found
m6
#> Error: object 'm6' not found
```
