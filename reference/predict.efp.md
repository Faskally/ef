# Predict capture probabilities from an efp fit

Predictions of capture probability from a fitted `efp` model, for the
data used to fit the model or for new data.

## Usage

``` r
# S3 method for class 'efp'
predict(
  object,
  newdata = NULL,
  type = c("link", "response", "lpmatrix"),
  se.fit = FALSE,
  exclude = NULL,
  ...
)
```

## Arguments

- object:

  an efp fitted model

- newdata:

  a data.frame containing the covariates in the model formula. If NULL,
  predictions are for the data used to fit the model.

- type:

  the type of prediction: `"link"` for the logit capture probability,
  `"response"` for the capture probability and `"lpmatrix"` for the
  design matrix which gives the logit capture probability when
  multiplied by the model parameters.

- se.fit:

  if TRUE standard errors of the predictions are also returned

- exclude:

  a character vector of terms to leave out of the prediction, named as
  in the model output, for example `"s(visitID)"`. This is the way to
  predict for levels of a random effect that were not in the data used
  to fit the model. Variables that are only used by excluded smooth or
  random effect terms need not be supplied in `newdata`.

- ...:

  further arguments passed to
  [`predict.gam`](https://rdrr.io/pkg/mgcv/man/predict.gam.html) when
  building the design matrix

## Value

a vector of predictions, or if `se.fit = TRUE` a list with elements
`fit` and `se.fit`. If `type = "lpmatrix"` a matrix.

## Details

Predictions are conditional on the estimated random effects, unless
these are removed using `exclude`. For a random walk smooth, predictions
can be made for any level in the penalty, including levels with no data,
which are interpolated. Standard errors on the response scale are
calculated using the delta method.

The `offset` supplied to
[`efp`](https://faskally.github.io/ef/reference/efp.md) is included in
the predictions only when `newdata` is NULL.

## Examples

``` r
data(ef_data)
ef_data$sampleID <- factor(paste(ef_data$date, ef_data$siteID,
                                 ef_data$lifestage, ef_data$species))
m1 <- efp(count ~ lifestage, data = ef_data, pass = pass, id = sampleID)

newdata <- data.frame(lifestage = c("Fry", "Parr"))
predict(m1, newdata = newdata, type = "response", se.fit = TRUE)
#> $fit
#>         1         2 
#> 0.6641846 0.7863670 
#> 
#> $se.fit
#>          1          2 
#> 0.01587840 0.01248373 
#> 
```
