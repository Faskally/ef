

#' Inverse logistic transformation
#'
#' This function transforms values on the logistic scale to values on
#' the probability scale.
#' 
#' @param x a numeric vector
#' @return vector of values between 0 and 1
#' @export
invlogit <- function(x) {
  1/(1 + exp(-x))
}
