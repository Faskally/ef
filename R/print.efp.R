#' @rdname efp-methods
#' @export
print.efp <- function(x, digits = max(3L, getOption("digits") - 3L),  ...) {
  print(summary(x), digits = digits, ...)
  invisible(x)
}
