# Utility function to convert parameters to probabilities

The matrix G should be of dimension n x p, and the parameter vector
should be length p

## Usage

``` r
transpar(par, G)
```

## Arguments

- par:

  fitted model parameters

- G:

  The design matrix for a model

## Value

a vector of probabilities of length n

## Examples

``` r
G <- cbind(1, c(0, 0, 1, 1))
transpar(c(0.5, -0.2), G)
#> [1] 0.6224593 0.6224593 0.5744425 0.5744425
```
