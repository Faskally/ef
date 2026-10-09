#' Methods for efp fits
#'
#' Methods to summarise and extract information from a model fitted
#' using \code{\link{efp}}.
#'
#' @param object,x an efp fitted model
#' @param type the type of fitted value: \code{"p"} for the capture
#'   probability in each pass, \code{"pi"} for the probability that a fish
#'   is first caught in each pass, \code{"both"} for a list containing both
#'   and \code{"lpmatrix"} for the design matrix.
#' @param parm the parameters to compute confidence intervals for, either a
#'   vector of numbers or names. The default is all fixed effects.
#' @param level the confidence level required
#' @param digits the number of significant digits to print
#' @param ... further arguments, ignored unless stated
#'
#' @details
#' The variance matrix returned by \code{vcov} is for the fixed effects
#' only, the joint variance matrix of the fixed and random effects is
#' available as \code{object$Vp}. Confidence intervals are Wald intervals
#' based on this variance matrix.
#'
#' The log likelihood is the marginal log likelihood when the model
#' contains random effects, and the degrees of freedom are the number of
#' fixed effects plus the number of random effect standard deviations.
#'
#' Residuals, analysis of deviance tables and diagnostic plots are not
#' available for efp fits, and the corresponding methods return an error.
#'
#' @return \code{fitted} returns a vector, list or matrix depending on
#'   \code{type}; \code{model.matrix} a matrix; \code{vcov} a matrix;
#'   \code{confint} a matrix with columns giving the lower and upper
#'   confidence limits; \code{logLik} an object of class \code{logLik};
#'   \code{formula} the model formula; \code{summary} an object of class
#'   \code{summary.efp}, a list containing the coefficient table
#'   (\code{coefficients}) and the table of random effect standard
#'   deviations (\code{sigma}). The \code{print} methods return their
#'   argument invisibly.
#'
#' @examples
#' data(ef_data)
#' ef_data$sampleID <- factor(paste(ef_data$date, ef_data$siteID,
#'                                  ef_data$lifestage, ef_data$species))
#' m1 <- efp(count ~ lifestage, data = ef_data, pass = pass, id = sampleID)
#'
#' summary(m1)
#' vcov(m1)
#' confint(m1)
#' logLik(m1)
#' head(fitted(m1))
#'
#' @name efp-methods
NULL

#' @rdname efp-methods
#' @export
vcov.efp <- function(object, ...) {
  if (!is.matrix(object$Vb)) {
    stop("the variance matrix is not available: refit the model with hessian = TRUE")
  }
  object$Vb
}

#' @rdname efp-methods
#' @importFrom stats confint.default logLik
#' @export
confint.efp <- function(object, parm, level = 0.95, ...) {
  confint.default(object, parm, level = level, ...)
}

#' @rdname efp-methods
#' @export
logLik.efp <- function(object, ...) {
  structure(
    object$llik,
    df = object$rank, nobs = nrow(object$G), class = "logLik"
  )
}

#' @rdname efp-methods
#' @export
formula.efp <- function(x, ...) {
  x$formula
}

#' @rdname efp-methods
#' @exportS3Method summary efp
summary.efp <- function(object, ...) {
  out <- object[c("call", "formula", "aic", "df.null", "df.residual", "convergence")]
  out$logLik <- logLik(object)

  fx <- summary(object$sdrep, "fixed", p.value = TRUE)
  fx <- fx[rownames(fx) == "beta", , drop = FALSE]
  rownames(fx) <- names(object$coefficients)
  out$coefficients <- fx

  if (length(object$sigma)) {
    re <- summary(object$sdrep, "report", p.value = FALSE)
    re <- re[rownames(re) == "sigma", , drop = FALSE]
    rownames(re) <- paste("sigma", object$re_labels)
    out$sigma <- re
  }

  class(out) <- "summary.efp"
  out
}

#' @rdname efp-methods
#' @exportS3Method print summary.efp
print.summary.efp <- function(x, digits = max(3L, getOption("digits") - 3L),  ...) {

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  cat("Formula:\n")
  print(x$formula)
  cat("\n")

  if (!is.null(x$sigma)) {
    cat("Random effects:\n")
    print(x$sigma, digits = digits, ...)
    cat("\n")
  }
  else {
    cat("No random effects\n\n")
  }

  if (nrow(x$coefficients)) {
    cat("Coefficients:\n")
    print(x$coefficients, digits = digits, ...)
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

#' @rdname efp-methods
#' @export
residuals.efp <- function(object, ...) {
  stop("residuals are not available for efp fits")
}

#' @rdname efp-methods
#' @export
anova.efp <- function(object, ...) {
  stop("anova is not available for efp fits, models can be compared using AIC")
}

#' @rdname efp-methods
#' @export
plot.efp <- function(x, ...) {
  stop("plot is not available for efp fits, see as.gam for plotting model terms")
}
