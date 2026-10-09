# Conversion of an efp fit to a gam fit for plotting

Converts an efp fit to a gam object, so that the plotting and prediction
tools in the `mgcv` package can be used. The coefficients and their
variances are those estimated by `efp`, but the smoothing parameters and
effective degrees of freedom in the gam object are not.

## Usage

``` r
as.gam(object)
```

## Arguments

- object:

  an efp fitted model

## Value

an object of class `gam`, see
[`gamObject`](https://rdrr.io/pkg/mgcv/man/gamObject.html)

## Examples

``` r
data(ef_data)
ef_data$sampleID <- factor(paste(ef_data$date, ef_data$siteID,
                                 ef_data$lifestage, ef_data$species))
m1 <- efp(count ~ lifestage + s(doy, k = 3),
          data = ef_data, pass = pass, id = sampleID)

m1_gam <- as.gam(m1)
#> Note that, although predictions and plots will be correct
#>  the smoothing parameters and degrees of freedom in the gam object are not those estimated by efp!
plot(m1_gam)

```
