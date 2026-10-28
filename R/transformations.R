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



#' Logistic transformation
#'
#' This function transforms values on the probability scale to values on
#' the logistic scale.
#'
#' @param p a numeric vector with values between 0 and 1
#' @return vector of values between -Inf and Inf
#' @export
logit <- function(p) {
  log(p / (1-p))
}



#' Utility function to convert parameters to probabilities
#'
#' The matrix G should be of dimension n x p,
#' and the parameter vector should be length p
#'
#' @param par fitted model parameters
#' @param G The design matrix for a model
#' @return a data frame
#' @export
transpar <- function(par, G) {
  1/(1 + exp(-c(G %*% par)))
}
