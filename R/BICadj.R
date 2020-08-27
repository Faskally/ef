#' Adjusted BIC function
#'
#' Complete function for returning overdispersion estimates
#'
#'
#' @param model a fitted ef model
#' @param data the data used to fit the model
#' @param overdispersion.output output form the function \link{overdispersion}
#' @return the adjusted BIC
#'
#' @note When calling the function, you need to specify the data source for
#'   the model so that the number of site visits can be determined.
#'   You also need to specify the output from the overdispersion model
#'   to get the measure
#'
#' @export
BICadj <- function(model, data, overdispersion.output) {

  # The adjusted BIC value is calculated in a similar way to BIC, however the
  # loglikelihood of the model is divided by the measure of between sample
  # overdispersion, which makes it more stringent when adding or removing terms
  # from the model.

  -2 * model$llik / overdispersion.output[3, 5]
  + log(length(unique(paste(data$siteID, data$date))))
      * length(model$coefficients)
}
