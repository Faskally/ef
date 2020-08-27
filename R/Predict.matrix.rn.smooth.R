#' @export
Predict.matrix.rn.smooth <- function(object, data) {
  X <- object$xt$X[, ncol(object$xt$X) - object$bs.dim:1 + 1]

  X[data[[object$term]], ]
}
