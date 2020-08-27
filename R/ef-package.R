#' Modelling Framework for the Estimation of Salmonid Abundance in Scottish Rivers
#'
#' A set of functions to estimate capture probabilities and densities
#' from multipass pass removal data.
#'
#' @docType package
#' @name ef-package
#' @aliases ef
#' @useDynLib ef, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @importFrom utils capture.output
#' @importFrom stats simulate AIC fitted as.formula formula binomial deviance
#' @importFrom stats update.formula drop.scope model.matrix
#' @importFrom stats aggregate pchisq qlogis
#' @importFrom mgcv gam smooth.construct Predict.matrix
#' @importFrom rstan optimizing
#' @importFrom dplyr left_join
#' @importFrom Matrix sparse.model.matrix colSums
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'

NULL
