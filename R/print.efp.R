#' @rdname efp-methods
#' @exportS3Method print efp
print.efp <- function(x, digits = max(3L, getOption("digits") - 3L),  ...) {
  print(summary(x), digits = digits, ...)
  invisible(x)
}
