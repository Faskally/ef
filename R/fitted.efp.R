#' @export
fitted.efp <- function(object, type = "p", ...) {
  type <- match.arg(type, c("p", "pi", "both", "lpmatrix"))
  if (type == "p") {
    object$p
  } else if (type == "pi") {
    object$pi
  } else if (type == "both") {
    object[c("p", "pi")]
  } else if (type == "lpmatrix") {
    object$G
  }
}
