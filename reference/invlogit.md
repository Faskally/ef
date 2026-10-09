# Inverse logistic transformation

This function transforms values on the logistic scale to values on the
probability scale.

## Usage

``` r
invlogit(x)
```

## Arguments

- x:

  a numeric vector

## Value

vector of values between 0 and 1

## Examples

``` r
invlogit(c(-2, 0, 2))
#> [1] 0.1192029 0.5000000 0.8807971
```
