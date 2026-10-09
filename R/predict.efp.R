#' Predict capture probabilities from an efp fit
#'
#' Predictions of capture probability from a fitted \code{efp} model, for
#' the data used to fit the model or for new data.
#'
#' @param object an efp fitted model
#' @param newdata a data.frame containing the covariates in the model
#'   formula. If NULL, predictions are for the data used to fit the model.
#' @param type the type of prediction: \code{"link"} for the logit capture
#'   probability, \code{"response"} for the capture probability and
#'   \code{"lpmatrix"} for the design matrix which gives the logit capture
#'   probability when multiplied by the model parameters.
#' @param se.fit if TRUE standard errors of the predictions are also returned
#' @param exclude a character vector of terms to leave out of the prediction,
#'   named as in the model output, for example \code{"s(visitID)"}. This is
#'   the way to predict for levels of a random effect that were not in the
#'   data used to fit the model. Variables that are only used by excluded
#'   smooth or random effect terms need not be supplied in \code{newdata}.
#' @param ... further arguments passed to \code{\link[mgcv]{predict.gam}}
#'   when building the design matrix
#' @return a vector of predictions, or if \code{se.fit = TRUE} a list with
#'   elements \code{fit} and \code{se.fit}. If \code{type = "lpmatrix"} a
#'   matrix.
#'
#' @details
#' Predictions are conditional on the estimated random effects, unless these
#' are removed using \code{exclude}. Standard errors on the response scale are
#' calculated using the delta method.
#'
#' The \code{offset} supplied to \code{\link{efp}} is included in the
#' predictions only when \code{newdata} is NULL.
#'
#' @examples
#' data(ef_data)
#' ef_data$sampleID <- factor(paste(ef_data$date, ef_data$siteID,
#'                                  ef_data$lifestage, ef_data$species))
#' m1 <- efp(count ~ lifestage, data = ef_data, pass = pass, id = sampleID)
#'
#' newdata <- data.frame(lifestage = c("Fry", "Parr"))
#' predict(m1, newdata = newdata, type = "response", se.fit = TRUE)
#'
#' @importFrom stats predict
#' @export
predict.efp <- function(object, newdata = NULL,
  type = c("link", "response", "lpmatrix"), se.fit = FALSE,
  exclude = NULL, ...) {

  type <- match.arg(type)

  if (is.null(newdata) && is.null(exclude)) {
    X <- object$G
  } else {
    g1 <- efp_gam(object)
    if (is.null(newdata)) {
      X <- predict(g1, type = "lpmatrix", exclude = exclude, ...)
    } else {
      # variables only used by excluded smooths do not affect the prediction
      # and may be missing, or have levels that were not in the fitted data
      labels <- vapply(g1$smooth, function(x) x$label, "")
      smooth_vars <- lapply(g1$smooth, function(x) setdiff(c(x$term, x$by), "NA"))
      used <- c(all.vars(g1$pterms), unlist(smooth_vars[!labels %in% exclude]))
      unused <- setdiff(unlist(smooth_vars[labels %in% exclude]), used)
      for (v in unused) {
        newdata[[v]] <- rep(object$data[[v]][1], nrow(newdata))
      }
      X <- predict(g1, newdata = newdata, type = "lpmatrix", exclude = exclude, ...)
    }
  }

  if (type == "lpmatrix") {
    return(X)
  }

  if (is.null(newdata) && !is.null(object$offset)) {
    offset <- object$offset
  } else {
    offset <- 0
  }

  # redundant parameters were not estimated and are zero
  b <- numeric(ncol(X))
  b[object$fx_cols] <- object$coefficients
  b[object$re_cols] <- object$random

  fit <- drop(X %*% b) + offset

  if (se.fit) {
    if (!is.matrix(object$Vp)) {
      stop("standard errors are not available: refit the model with hessian = TRUE")
    }
    whichkeep <- sort(c(object$fx_cols, object$re_cols))
    Xfit <- X[, whichkeep, drop = FALSE]
    se <- sqrt(pmax(0, rowSums((Xfit %*% object$Vp) * Xfit)))
  }

  if (type == "response") {
    fit <- invlogit(fit)
    if (se.fit) {
      se <- se * fit * (1 - fit)
    }
  }

  if (se.fit) {
    list(fit = fit, se.fit = se)
  } else {
    fit
  }
}
