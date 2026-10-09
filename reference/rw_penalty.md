# Random walk penalties for smooth terms

Creates a first or second order random walk penalty matrix for the
levels of an ordered covariate, for use with the Markov random field
(mrf) basis of `mgcv`. This is the recommended way to include smooth
terms in [`efp`](https://faskally.github.io/ef/reference/efp.md): the
basis and the penalty are both sparse, so models are quick to fit even
when there are many levels.

## Usage

``` r
rw_penalty(levels, order = 1)
```

## Arguments

- levels:

  the levels of the covariate, in order, and assumed to be equally
  spaced. A level is needed for every step of the random walk, including
  any for which there are no data.

- order:

  the order of the random walk, 1 or 2

## Value

a square matrix with one row and column for each level, with the levels
as dimnames. The matrix has an attribute, `null.space`, which is a basis
for the functions which are not penalised.

## Details

A first order random walk (`order = 1`) penalises the differences
between neighbouring levels, and shrinks towards a constant. A second
order random walk (`order = 2`) penalises the second differences, and
shrinks towards a linear trend. In both cases the standard deviation
estimated by [`efp`](https://faskally.github.io/ef/reference/efp.md) is
the standard deviation of the steps of the random walk, on the logit
scale.

To use the penalty, the covariate must be a factor with the same levels
as the penalty, see the examples.

## Examples

``` r
rw_penalty(1:5)
#>    1  2  3  4  5
#> 1  1 -1  0  0  0
#> 2 -1  2 -1  0  0
#> 3  0 -1  2 -1  0
#> 4  0  0 -1  2 -1
#> 5  0  0  0 -1  1
#> attr(,"null.space")
#>            [,1]
#> [1,] -0.4472136
#> [2,] -0.4472136
#> [3,] -0.4472136
#> [4,] -0.4472136
#> [5,] -0.4472136
rw_penalty(1:5, order = 2)
#>    1  2  3  4  5
#> 1  1 -2  1  0  0
#> 2 -2  5 -4  1  0
#> 3  1 -4  6 -4  1
#> 4  0  1 -4  5 -2
#> 5  0  0  1 -2  1
#> attr(,"null.space")
#>            [,1]          [,2]
#> [1,] -0.4472136 -6.324555e-01
#> [2,] -0.4472136 -3.162278e-01
#> [3,] -0.4472136  8.326673e-17
#> [4,] -0.4472136  3.162278e-01
#> [5,] -0.4472136  6.324555e-01

# a first order random walk over day of the year
data(ef_data)
ef_data$sampleID <- factor(paste(ef_data$date, ef_data$siteID,
                                 ef_data$lifestage, ef_data$species))
days <- min(ef_data$doy):max(ef_data$doy)
ef_data$fdoy <- factor(ef_data$doy, levels = days)
rw1 <- rw_penalty(days)

m1 <- efp(count ~ lifestage + s(fdoy, bs = "mrf", xt = list(penalty = rw1)),
          data = ef_data, pass = pass, id = sampleID)
m1
#> 
#> Call:
#> efp(formula = count ~ lifestage + s(fdoy, bs = "mrf", xt = list(penalty = rw1)), 
#>     data = ef_data, pass = pass, id = sampleID)
#> 
#> Formula:
#> count ~ lifestage + s(fdoy, bs = "mrf", xt = list(penalty = rw1))
#> <environment: 0x55888bb86b60>
#> 
#> Random effects:
#>               Estimate Std. Error
#> sigma s(fdoy) 0.001289    0.01585
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z^2|)
#> (Intercept)     0.6820    0.07131   9.564  1.138e-21
#> lifestageParr   0.6212    0.10292   6.036  1.577e-09
#> 
#> Degrees of Freedom: 90 Total (i.e. Null);  87 Residual
#> AIC: 3194
```
