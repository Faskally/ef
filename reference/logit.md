# Logistic transformation

This function transforms values on the probability scale to values on
the logistic scale.

## Usage

``` r
logit(p)
```

## Arguments

- p:

  a numeric vector with values between 0 and 1

## Value

vector of values between -Inf and Inf

## Examples

``` r
logit(c(0.1, 0.5, 0.9))
#> [1] -2.197225  0.000000  2.197225
```
