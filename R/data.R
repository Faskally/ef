#' Counts of salmon fry and parr from electrofishing.
#'
#' A dataset containing counts from multipass electrofishing
#' samples over three sites and 5 years.
#'
#' @format A data frame with 90 rows and 9 variables:
#' \describe{
#'   \item{siteID}{Unique identifier for each elecrofishing site}
#'   \item{year}{year of data collection}
#'   \item{date}{date of data collections}
#'   \item{EFpasscount:}{the total number of electrofishing passes}
#'   \item{pass}{the electrofishing pass on which the fish were caught}
#'   \item{species}{fish species, Salmon}
#'   \item{lifestage}{fish lifestage, Fry or Parr}
#'   \item{count}{the number of fish caught per pass for each site visit and species, etc.}
#'   \item{area}{the area of river fished}
#' }
#'
#' @source \url{https://www2.gov.scot/Topics/marine/Salmon-Trout-Coarse/Freshwater/Monitoring/temperature}
"ef_data"
