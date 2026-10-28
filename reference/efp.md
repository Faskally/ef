# Estimate capture probabilites from electrofishing data

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
  init = "0",
  hessian = TRUE,
  fit = TRUE,
  sample_re = FALSE,
  control = list()
)
```

## Arguments

- formula:

  a formula object

- data:

  a data.frame containing all relavent info

- pass:

  a vector of integers giving the pass number of the observation

- id:

  a vector of integers identifying an observation (a set of
  electrofishing passes)

- offset:

  an possible offset for the linear predictor of capture probability

- verbose:

  if TRUE optimiser messages are printed to the screen

- init:

  should initialisatiom be random?

- hessian:

  if TRUE the hessian is computed and the covariance matrix of the
  parameters is returned via Vb

- fit:

  if TRUE model is fitted if FALSE the data that would be passed to the
  optimiser is returned

- sample_re:

  should sample random effects be included

- control:

  a list of control commands to be passed to optim

## Value

glm type object

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
#> 'as(<matrix>, "dgTMatrix")' is deprecated.
#> Use 'as(as(as(., "dMatrix"), "generalMatrix"), "TsparseMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").
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
```
