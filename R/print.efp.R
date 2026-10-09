#' @importFrom stats formula coef
#' @export
print.efp <- function(x, digits = max(3L, getOption("digits") - 3L),  ...) {

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  cat("Formula:\n")
  print(x$formula)
  cat("\n")

  if (!is.null(x$sdrep$par.random)) {
    cat("Random effects:\n")
    re <- summary(x$sdrep, "report", p.value = FALSE)
    re <- re[rownames(re) == "sigma", , drop = FALSE]
    rownames(re) <- paste("sigma", x$re_labels)
    print(re, digits = digits, ...)
    cat("\n")
  }
  else {
    cat("No random effects\n\n")
  }

  if (!is.null(x$sdrep$par.fixed)) {
    cat("Coefficients:\n")
    fx <- summary(x$sdrep, "fixed", p.value = TRUE)
    fx <- fx[!grepl("log_", rownames(fx)), ,drop = FALSE]
    rownames(fx) <- names(coef(x))
    print(fx, digits = digits, ...)
  }
  else {
    cat("No coefficients\n\n")
  }

  cat(
    "\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
    x$df.residual, "Residual\n"
  )

  cat("AIC:", format(signif(x$aic, digits)))

  cat("\n")
  invisible(x)
}