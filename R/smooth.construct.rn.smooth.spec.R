#' @export
smooth.construct.rn.smooth.spec <- function (object, data, knots)
{
  if (!is.null(object$id))
    stop("rn smoothers don't work with ids.")
  # set basis dimension
  if (object$bs.dim == -1) {
    object$X <- object$xt$X
  } else {
    if (object$bs.dim > ncol(object$xt$X) + 1) {
      warning("dimension of smoother too large for supplied X matrix, reducing to k = ", ncol(object$xt$X) + 1)
      object$bs.dim <- ncol(object$xt$X) + 1
    }
    object$X <- object$xt$X[,ncol(object$xt$X) - object$bs.dim:1 + 1]
  }
  object$bs.dim <- ncol(object$X)
  # get correct rows
  object$X <- object$X[data[[object$term]],]
  # set penalty matrix
  if (is.null(object$xt$S)) {
    object$S <- list(diag(object$bs.dim))
  } else {
    object$S <- list(object$xt$S)
  }

  # set smoother properties
  object$rank <- object$bs.dim
  object$null.space.dim <- 0
  object$C <- matrix(0, 0, ncol(object$X))
  object$side.constrain <- FALSE
  object$plot.me <- FALSE
  object$te.ok <- if (inherits(object, "tensor.smooth.spec")) 0 else 2
  object$random <- TRUE
  class(object) <- "rn.smooth"
  object
}
