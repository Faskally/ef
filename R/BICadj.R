#' Adjusted BIC function
#'
#' Complete function for returning overdispersion estimates
#'
#'
#' @param model a fitted ef model
#' @param overdispersion.output output from the function \link{overdispersion}
#' @return the adjusted BIC
#'
#' @note VisitID calculated from the data contained in the model object
#'
#' @export
BICadj <- function(model, overdispersion.output) {

  # The adjusted BIC value is calculated in a similar way to BIC, however the
  # loglikelihood of the model is divided by the greatest measure of 
  # overdispersion (within or between visits), which makes it more stringent 
  # when adding or removing terms from the model.

  -2 * model$llik / max(overdispersion.output[c("sitevisit", "large"), "disp"]) +
   log(length(unique(model$data$visitID))) * length(model$coefficients)
}
