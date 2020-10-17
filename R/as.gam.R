#' Coversion of an efp fit to a gam fit for plotting
#'
#' Coversion of an efp fit to a gam fit for plotting
#'
#'
#' @param object an efp fitted model
#' @return gam type object
#'
#' @importFrom mgcv gam
#' @importFrom stats binomial
#' @importFrom methods is
#'
#' @export
as.gam <- function(object) {
  stopifnot(is(object)[1] == "efp")
  # make a gam container from the originoal setup
  g1 <- gam(G = object $ Gsetup)

  # recheck for rank deficiency
  qr.G <- qr(object$G)
  rank.deficient <- qr.G$pivot[abs(diag(qr.G$qr)) < 1e-7]
  whichkeep <- -rank.deficient
  if (!length(whichkeep)) {
    whichkeep <- 1:length(object$ coefficients)
  }
  redudant <- names(g1$coefficients[-1 * whichkeep])
  if (length(redudant)) {
    warning(
      "there is some reduandancy in the model specification.",
      "This may be fine, but best to check."
    )
  }

  # assign coefficients to gam container
  g1$coefficients[] <- 0
  g1$coefficients[whichkeep] <- object $ coefficients

  # assign variance matrix to gam container
  g1$Vp[] <- 0
  diag(g1$Vp[]) <- 1e-5
  g1$Vp[whichkeep, whichkeep] <- object$Vb

  # assign binomial family to gam container
  g1$family <- binomial()

  message("Note that, although predictions and plots will be correct\n smoothing parameter estimates are all wrong as no smoothing took place!")
  g1
}
