

#' @export
invlogit <- function(x) {
  1/(1 + exp(-x))
}
