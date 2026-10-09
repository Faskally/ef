#' Conversion of an efp fit to a gam fit for plotting
#'
#' Converts an efp fit to a gam object, so that the plotting and prediction
#' tools in the \code{mgcv} package can be used. The coefficients and their
#' variances are those estimated by \code{efp}, but the smoothing parameters
#' and effective degrees of freedom in the gam object are not.
#'
#'
#' @param object an efp fitted model
#' @return an object of class \code{gam}, see \code{\link[mgcv]{gamObject}}
#'
#' @examples
#' data(ef_data)
#' ef_data$sampleID <- factor(paste(ef_data$date, ef_data$siteID,
#'                                  ef_data$lifestage, ef_data$species))
#' m1 <- efp(count ~ lifestage + s(doy, k = 3),
#'           data = ef_data, pass = pass, id = sampleID)
#'
#' m1_gam <- as.gam(m1)
#' plot(m1_gam)
#'
#' @importFrom mgcv gam
#' @importFrom stats binomial
#' @importFrom methods is
#'
#' @export
as.gam <- function(object) {
  stopifnot(is(object)[1] == "efp")
  g1 <- efp_gam(object)

  # columns of the design matrix that were estimated
  whichkeep <- object$gam_cols
  redundant <- names(g1$coefficients[-1 * whichkeep])
  if (length(redundant)) {
    warning(
      "there is some redundancy in the model specification.",
      "This may be fine, but best to check."
    )
  }

  message("Note that, although predictions and plots will be correct\n the smoothing parameters and degrees of freedom in the gam object are not those estimated by efp!")
  g1
}

# gam container holding the efp estimates, used for prediction and plotting
efp_gam <- function(object) {
  # make a gam container from the original setup
  g1 <- gam(G = object $ Gsetup)

  # columns of the design matrix that were estimated
  whichkeep <- object$gam_cols

  # assign coefficients to gam container
  g1$coefficients[] <- object$gam_coefficients

  # assign variance matrix to gam container
  g1$Vp[] <- 0
  diag(g1$Vp[]) <- 1e-5
  if (is.matrix(object$Vp)) {
    g1$Vp[whichkeep, whichkeep] <- object$Vp
  }

  # assign binomial family to gam container
  g1$family <- binomial()

  g1
}
