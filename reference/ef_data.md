# Counts of salmon fry and parr from electrofishing.

A dataset containing counts from multipass electrofishing samples over
three sites and 5 years.

## Usage

``` r
ef_data
```

## Format

A data frame with 90 rows and 10 variables:

- siteID:

  Unique identifier for each electrofishing site

- year:

  year of data collection

- date:

  date of data collections

- doy:

  day of the year of data collection, a continuous covariate which can
  be used to fit smooth terms

- EFpasscount:

  the total number of electrofishing passes

- pass:

  the electrofishing pass on which the fish were caught

- species:

  fish species, Salmon

- lifestage:

  fish lifestage, Fry or Parr

- count:

  the number of fish caught per pass for each site visit and species,
  etc.

- area:

  the area of river fished

## Source

<https://www.gov.scot/publications/scotland-river-temperature-monitoring-network-srtmn/>
