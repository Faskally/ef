#' @export
print.efp <- function(x, digits = max(3L, getOption("digits") - 3L),  ...) {

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  cat("Formula:\n")
  print(formula(x))
  cat("\n")

  if (!is.null(x$sdrep$par.random)) {
    cat("Random effects:\n")
    print(summary(x$sdrep, "report", p.value = FALSE), digits = digits, ...)
    cat("\n")
  }
  else {
    cat("No random effects\n\n")
  }

  if (length(coef(x))) {
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

  if (nzchar(mess <- naprint(x$na.action))) {
    cat("  (", mess, ")\n", sep = "")
  }

  cat("AIC:", format(signif(x$aic, digits)))

  cat("\n")
  invisible(x)
}