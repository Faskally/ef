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
#' @param ... further arguments, ignored
#' @return a vector of predictions, or if \code{se.fit = TRUE} a list with
#'   elements \code{fit} and \code{se.fit}. If \code{type = "lpmatrix"} a
#'   matrix.
#'
#' @details
#' Predictions are conditional on the estimated random effects, unless these
#' are removed using \code{exclude}. For a random walk smooth, predictions can
#' be made for any level in the penalty, including levels with no data, which
#' are interpolated. Standard errors on the response scale are
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
#' @importFrom stats predict model.frame model.matrix delete.response
#' @importFrom mgcv PredictMat
#' @export
predict.efp <- function(object, newdata = NULL,
  type = c("link", "response", "lpmatrix"), se.fit = FALSE,
  exclude = NULL, ...) {

  type <- match.arg(type)

  X <- efp_lpmatrix(object, newdata, exclude)

  if (type == "lpmatrix") {
    return(X)
  }

  if (is.null(newdata) && !is.null(object$offset)) {
    offset <- object$offset
  } else {
    offset <- 0
  }

  fit <- drop(X %*% object$gam_coefficients) + offset

  if (se.fit) {
    if (!is.matrix(object$Vp)) {
      stop("standard errors are not available: refit the model with hessian = TRUE")
    }
    Xfit <- X[, object$gam_cols, drop = FALSE]
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

# the design matrix for new data, with the columns of excluded terms set to
# zero
efp_lpmatrix <- function(object, newdata = NULL, exclude = NULL) {
  Gsetup <- object$Gsetup
  labels <- vapply(Gsetup$smooth, function(x) x$label, "")
  excluded <- labels %in% exclude

  if (is.null(newdata)) {
    X <- object$G
  } else {
    newdata <- as.data.frame(newdata)

    # variables only used by excluded smooths do not affect the prediction
    # and may be missing, or have levels that were not in the fitted data
    smooth_vars <- lapply(Gsetup$smooth, function(x) setdiff(c(x$term, x$by), "NA"))
    used <- c(all.vars(Gsetup$pterms), unlist(smooth_vars[!excluded]))
    unused <- setdiff(unlist(smooth_vars[excluded]), used)
    for (v in unused) {
      newdata[[v]] <- rep(Gsetup$mf[[v]][1], nrow(newdata))
    }

    # factors need the levels used in the fit
    for (v in intersect(names(newdata), names(Gsetup$mf))) {
      if (is.factor(Gsetup$mf[[v]])) {
        x <- factor(newdata[[v]], levels = levels(Gsetup$mf[[v]]))
        if (any(is.na(x) & !is.na(newdata[[v]]))) {
          stop("newdata has levels of ", v, " which were not in the fitted data")
        }
        newdata[[v]] <- x
      }
    }

    X <- matrix(0, nrow(newdata), ncol(object$G))
    colnames(X) <- colnames(object$G)

    # parametric terms
    if (Gsetup$nsdf > 0) {
      pterms <- delete.response(Gsetup$pterms)
      mf <- model.frame(pterms, newdata, xlev = Gsetup$xlevels, na.action = stats::na.pass)
      X[, seq_len(Gsetup$nsdf)] <- model.matrix(pterms, mf, contrasts.arg = Gsetup$contrasts)
    }

    # smooth terms
    for (sm in Gsetup$smooth[!excluded]) {
      X[, sm$first.para:sm$last.para] <- PredictMat(sm, newdata)
    }
  }

  # set excluded terms to zero
  for (sm in Gsetup$smooth[excluded]) {
    X[, sm$first.para:sm$last.para] <- 0
  }
  pexclude <- which(attr(Gsetup$pterms, "term.labels") %in% exclude)
  if (length(pexclude)) {
    X[, which(Gsetup$assign %in% pexclude)] <- 0
  }

  X
}
