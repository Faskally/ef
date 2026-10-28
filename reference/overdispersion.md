# Estimating overdispersion

Complete function for returning overdispersion estimates

## Usage

``` r
overdispersion(
  data,
  visitID = "visitID",
  count = "count",
  pass = "pass",
  sampleID = "sampleID",
  largemodel,
  control = "control"
)
```

## Arguments

- data:

  dataframe containing EF data

- visitID:

  a number identifying each unique site visit

- count:

  the number of fish caught for a particular combination of site visit,
  species, lifestage and pass (defaults to "count")

- pass:

  EF pass number e.g. 1,2,3,4 (defaults to "pass")

- sampleID:

  sample ID i.e. unique combinations of site visit, species & lifestage

- largemodel:

  a large model that captures most of the systematic variation in the
  data - this is specified before running the overdispersion function

- control:

  passes control information to optimiser

## Value

a data.frame summarising overdispersion

## Note

ensure column names in function call are in inverted commas
