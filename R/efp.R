#' Estimate capture probabilities from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters
#'
#'
#' @param formula a formula object. Random effects are specified using
#'   \code{s(f, bs = "re")} terms, where \code{f} is a factor, see Details.
#' @param data a data.frame containing all relevant info
#' @param pass a vector of integers giving the pass number of the observation
#' @param id a vector of integers identifying an observation (a set of electrofishing passes)
#' @param offset an optional offset for the linear predictor of capture probability
#' @param verbose if TRUE optimiser messages are printed to the screen
#' @param hessian if TRUE the hessian is computed and the covariance matrix of the parameters is returned via Vb
#' @param fit if TRUE model is fitted if FALSE the data that would be passed to the optimiser is returned
#' @param control a list of control commands to be passed to \code{\link[stats]{optim}}
#' @param sparse if TRUE (the default) penalised smooths are set up so that
#'   sparse bases and penalties stay sparse, see Details. If FALSE the
#'   constrained basis and penalties from \code{mgcv} are used, which are
#'   dense. Both give the same model, so this is mainly useful for testing.
#' @return If \code{fit = TRUE}, an object of class \code{efp}, a list which
#'   inherits from \code{glm}, containing among other things the fixed effect
#'   estimates (\code{coefficients}), the random effect estimates
#'   (\code{random}), the random effect standard deviations (\code{sigma}),
#'   the fitted capture probabilities (\code{fitted}), the log likelihood
#'   (\code{llik}) and the variance matrix of the fixed effects (\code{Vb}).
#'   See \code{\link{efp-methods}} for the methods available. If
#'   \code{fit = FALSE}, a list containing the data that would be passed to
#'   the optimiser.
#'
#' @references
#' Millar, C. P., Fryer, R. J., Millidine, K. J. and Malcolm, I. A. (2016)
#' Modelling capture probability of Atlantic salmon (Salmo salar) from a
#' diverse national electrofishing dataset: Implications for the estimation
#' of abundance. Fisheries Research, 177, 1-12.
#' \doi{10.1016/j.fishres.2016.01.001}
#'
#' @details
#' Random effects on the logit capture probability are included by adding
#' \code{s(f, bs = "re")} terms to the formula, where \code{f} is a factor
#' in \code{data}. For example, \code{s(visitID, bs = "re")} adds between
#' visit variation and \code{s(sampleID, bs = "re")} adds between sample
#' variation. Each term has its own standard deviation, and the
#' random effects are integrated out of the likelihood using the Laplace
#' approximation, so that \code{logLik} and \code{AIC} are based on the marginal
#' likelihood. Fitted values are conditional on the estimated random effects.
#'
#' Smooth terms, such as \code{s(x, k = 3)}, can also be included in the
#' formula. These are fitted as unpenalised regression splines: no smoothing
#' parameters are estimated, so the flexibility of a smooth is set entirely
#' by its basis dimension \code{k}, which should be kept small (typically
#' \code{k = 3}, i.e. 2 degrees of freedom).
#'
#' @examples
#' # create two electrofishing site visits with 3 and 4 passes and 2 lifestages
#' ef_data <- data.frame(n      = c(100, 53, 24, 50, 26, 12,
#'                                  100, 53, 24, 50, 26, 12),
#'                       pass   = c(  1,  2,  3,  1,  2,  3,
#'                                    1,  2,  3,  1,  2,  3),
#'                       stage  = c(  1,  1,  1,  2,  2,  2,
#'                                    1,  1,  1,  2,  2,  2),
#'                       sample = c(  1,  1,  1,  2,  2,  2,
#'                                    3,  3,  3,  4,  4,  4))
#'
#' ef_data2 <- data.frame(n      = c(100, 53, 24, 50, 26, 12,
#'                                   100, 53, 24, 12, 50, 26, 12, 6),
#'                        pass   = c(  1,  2,  3,  1,  2,  3,
#'                                     1,  2,  3,  4,  1,  2,  3, 4),
#'                        stage  = c(  1,  1,  1,  2,  2,  2,
#'                                     1,  1,  1,  1,  2,  2,  2, 2),
#'                        sample = c(  1,  1,  1,  2,  2,  2,
#'                                     3,  3,  3,  3,  4,  4,  4, 4))
#'
#' ef_data3 <- data.frame(n      = c(100, 53, 24, 50, 26, 12, 40,
#'                                   100, 53, 24, 12, 50, 26, 12, 6, 40),
#'                        pass   = c(  1,  2,  3,  1,  2,  3, 1,
#'                                     1,  2,  3,  4,  1,  2,  3, 4, 1),
#'                        stage  = c(  1,  1,  1,  2,  2,  2, 1,
#'                                     1,  1,  1,  1,  2,  2,  2, 2, 2),
#'                        sample = c(  1,  1,  1,  2,  2,  2, 5,
#'                                     3,  3,  3,  3,  4,  4,  4, 4, 6))
#'
#' # Fit a simple model
#' m2 <- efp(n ~ 1 + factor(stage), data = ef_data, pass = pass, id = sample)
#' cbind(ef_data, fit = fitted(m2))
#' m3 <- efp(n ~ 1 + factor(stage), data = ef_data2, pass = pass, id = sample)
#' cbind(ef_data2, fit = fitted(m3))
#' m4 <- efp(n ~ 1 + factor(stage), data = ef_data3, pass = pass, id = sample)
#' cbind(ef_data3, fit = fitted(m4))
#'
#' # create two electrofishing site visits with 3 and 4 passes and 2 lifestages
#' ef_data <- data.frame(n      = c(200, 53, 24, 100, 26, 12,
#'                                  200, 53, 24, 100, 26, 12),
#'                       pass   = c(  1,  2,  3,  1,  2,  3,
#'                                    1,  2,  3,  1,  2,  3),
#'                       stage  = c(  1,  1,  1,  2,  2,  2,
#'                                    1,  1,  1,  2,  2,  2),
#'                       sample = c(  1,  1,  1,  2,  2,  2,
#'                                    3,  3,  3,  4,  4,  4))
#' # Fit a simple model
#' m2 <- efp(n ~ 1 + factor(stage) + factor(replace(pass, pass> 2, 2)),
#'           data = ef_data, pass = pass, id = sample)
#' out <- cbind(ef_data, p = fitted(m2, type = "p"))
#' out
#'
#' # a smooth effect of a continuous covariate, day of the year. Smooths are unpenalised, so
#' # the flexibility is set by k
#' data(ef_data)
#' ef_data$sampleID <- factor(paste(ef_data$date, ef_data$siteID,
#'                                  ef_data$lifestage, ef_data$species))
#' ef_data$pass12 <- factor(replace(ef_data$pass, ef_data$pass > 2, 2))
#' m5 <- efp(count ~ siteID + pass12 + lifestage + s(doy, k = 3),
#'           data = ef_data, pass = pass, id = sampleID)
#' m5
#'
#' # predicted capture probability over the smooth, with standard errors
#' newdata <- data.frame(siteID = "Site1", pass12 = "1", lifestage = "Fry",
#'                       doy = seq(195, 215, by = 5))
#' pred <- predict(m5, newdata = newdata, type = "response", se.fit = TRUE)
#' cbind(newdata, fit = pred$fit, se = pred$se.fit)
#'
#' @importFrom RTMB MakeADFun sdreport ADREPORT REPORT dnorm dgmrf
#' @importFrom mgcv gam
#' @importFrom methods as
#' @import Matrix
#' @importFrom stats binomial
#' @importFrom stats optim
#' @importFrom stats ave
#'
#' @export
efp <- function(formula, data = NULL, pass = pass, id = id,
  offset = NULL, verbose = FALSE, hessian = TRUE,
  fit = TRUE, control = list(), sparse = TRUE) {

  # some checks
  if (is.null(data)) stop("must supply data")
  
  # get pass information
  pass <- substitute(pass)
  pass <- eval(pass, data, environment(formula))
  # get id
  id <- substitute(id)
  id <- eval(id, data, environment(formula))
  
  # this is necessary to ensure that passes are ordered within samples so that
  # the order of input data and model outputs match. Likelihood function 
  # requires data ordered by pass.
  if (any(tapply(pass, id, is.unsorted))) {
    stop("data must be sorted by pass")
  }

  # set up offset
  if (is.null(offset)) {
    offset <- rep(0, nrow(data))
  }

  # set up model
  Gsetup <- gam(formula, data = data, fit = FALSE, drop.unused.levels = FALSE)
  G <- Gsetup$X
  colnames(G) <- Gsetup$term.names
  if (nrow(G) != nrow(data)) stop("I think there are NAs in your data, please check and remove them before rerunning.")

  # identify random effects: these are the s(..., bs = "re") terms
  is_re <- vapply(Gsetup$smooth, inherits, TRUE, what = "random.effect")
  re_cols <- lapply(Gsetup$smooth[is_re], function(x) x$first.para:x$last.para)
  re_cols <- as.integer(unlist(re_cols))

  # set up the penalised terms: random effects and penalised smooths. Each
  # is split into fixed effects, for the part of the term that is not
  # penalised, and random effects
  blocks <- smooth_setup(Gsetup, formula, sparse = sparse)
  pen_cols <- as.integer(unlist(lapply(blocks, function(x) x$cols)))
  fx_cols <- setdiff(seq_len(ncol(G)), pen_cols)
  nsigma <- vapply(blocks, function(x) x$npen, 0)
  re_labels <-
    unlist(lapply(blocks, function(x) {
      if (x$npen == 1) x$label else paste0(x$label, seq_len(x$npen))
    }))
  if (is.null(re_labels)) re_labels <- character(0)

  # the fixed effects design matrix
  X <- do.call(cbind, c(list(G[, fx_cols, drop = FALSE]), lapply(blocks, function(x) x$Xf)))
  Z <- do.call(cbind, c(list(Matrix(0, nrow(G), 0)), lapply(blocks, function(x) x$Zu)))

  # remove redundant / collinear parameters (fixed effects only)
  qr.X <- qr(X)
  rank.deficient <- qr.X $ pivot[abs(diag(qr.X $ qr)) < 1e-7]
  if (length(rank.deficient)) {
    droppar <- paste(colnames(X)[rank.deficient], collapse = "\n\t")
    warning("*** Model has ", length(rank.deficient)," too many parameter(s)!!\n    i will remove the redundant ones:\n\t", droppar, call. = FALSE)
  }
  fx_keep <- setdiff(seq_len(ncol(X)), rank.deficient)
  fx_cols <- setdiff(fx_cols, fx_cols[rank.deficient[rank.deficient <= length(fx_cols)]])

  # the matrix which converts the model parameters (fixed then random
  # effects) to the coefficients of the columns of G
  nfx <- ncol(G) - length(pen_cols)
  nX <- ncol(X)
  to_gam <- Matrix(0, ncol(G), nX + ncol(Z), sparse = TRUE)
  to_gam[cbind(setdiff(seq_len(ncol(G)), pen_cols), seq_len(nfx))] <- 1
  at_g <- nfx
  at_u <- nX
  for (i in seq_along(blocks)) {
    block <- blocks[[i]]
    ipar <- c(at_g + seq_len(ncol(block$Xf)), at_u + seq_len(ncol(block$Zu)))
    to_gam[block$cols, ipar] <- block$TG
    if (!is.null(block$Tint)) {
      to_gam[which(colnames(G) == "(Intercept)"), ipar] <- block$Tint
    }
    # where the random effects and standard deviations are
    blocks[[i]]$idx <- at_u - nX + seq_len(ncol(block$Zu))
    blocks[[i]]$sig <- sum(nsigma[seq_len(i - 1)]) + seq_len(block$npen)
    at_g <- at_g + ncol(block$Xf)
    at_u <- at_u + ncol(block$Zu)
  }
  to_gam <- to_gam[, c(fx_keep, nX + seq_len(ncol(Z))), drop = FALSE]
  gam_cols <- which(rowSums(abs(to_gam)) > 0)
  Gfit <- G[, gam_cols, drop = FALSE]

  # get data for likelihood
  ord_output <- order(id, pass)

  efdat <-
    list(
      y = Gsetup$y,
      offset = offset,
      sample_id = factor(id),
      X = as(as(X[, fx_keep, drop = FALSE], "generalMatrix"), "TsparseMatrix"),
      Z = as(as(Z, "generalMatrix"), "TsparseMatrix"),
      re_id = as.integer(unlist(lapply(blocks, function(x) rep(x$sig[1], length(x$idx))))),
      blocks = lapply(blocks, function(x) x[c("idx", "sig", "S", "logdet")])
    )

  if (!fit) {
    return(efdat)
  }

  params <-
    list(
      alpha = rep(0, ncol(efdat$Z)),
      beta = rep(0, ncol(efdat$X)),
      log_sigma = rep(0, length(re_labels))
    )

  if (length(blocks)) {
    rand <- "alpha"
  } else {
    rand <- NULL
  }

  obj <-
    MakeADFun(
      function(params) ef_nll(params, efdat),
      params,
      random = rand,
      silent = !verbose
    )

  opt <- optim(obj$par, obj$fn, obj$gr, method = "BFGS", control = control)
  
  # check convergence
  
  if (opt$convergence != 0) {
    warning("model has not converged: optim convergence code = ", opt$convergence)
  }
  
  # set up output object
  out <- list()

  # evaluate at the optimum, this also sets the random effects to their
  # estimates at the optimum
  out$llik <- -1 * as.numeric(obj$fn(opt$par))
  par_full <- obj$env$last.par
  out$rep <- obj$report(par_full)

  # extract transformed parameters
  out$p <- unname(unlist(out$rep$ps))
  # Re-order to match input data
  out$p <- out$p[order(ord_output)]
  # probability of capture in each pass, allowing for the fish removed in
  # previous passes
  out$pi <- as.vector(
    ave(out$p, factor(id), FUN = function(x) x * c(1, cumprod(1 - x)[-length(x)]))
  )
  
  fx_pars <- grep("beta", names(opt$par))
  out$sdrep <- sdreport(obj, getJointPrecision = hessian && length(blocks) > 0)
  # keep data for re-fitting
  # efdat needs to be reordered to be compatible with other structures 
  out$data <- data
  out$offset <- offset
  out$formula <- formula # for printing and summary
  out$terms <- Gsetup$terms
  out$call <- match.call()
  out$rank <- length(opt$par)
  out$aic <- -2 * out$llik + 2 * out$rank
  out$G <- G
  out$Gfit <- Gfit
  out$fx_cols <- fx_cols
  out$re_cols <- re_cols
  out$gam_cols <- gam_cols
  out$re_labels <- re_labels
  out$coefficients <- opt$par[fx_pars]
  names(out$coefficients) <- colnames(X)[fx_keep]
  out$random <- par_full[names(par_full) == "alpha"]
  names(out$random) <- colnames(Z)
  out$sigma <- exp(opt$par[names(opt$par) == "log_sigma"])
  names(out$sigma) <- re_labels
  # coefficients for the columns of G, used for prediction
  out$gam_coefficients <- drop(as.matrix(to_gam %*% c(out$coefficients, out$random)))
  names(out$gam_coefficients) <- colnames(G)
  out$df.null <- nrow(G)
  out$df.residual <- nrow(G) - out$rank
  out$fitted <- out$p

  out$residuals <- rep(0, nrow(data))
  out$null.deviance <- NA
  out$deviance <- NA
  out$family <- binomial()
  if (!hessian) {
    out$hessian <- NULL
    out$Vb <- NULL
    out$Vp <- NULL
  } else if (length(blocks)) {
    # covariance of fixed and random effects
    out$hessian <- NULL
    V <- try(as.matrix(solve(out$sdrep$jointPrecision)))
    if (inherits(V, "try-error")) {
      out$Vb <- out$Vp <- V
    } else {
      nms <- rownames(out$sdrep$jointPrecision)
      which_b <- which(nms == "beta")
      which_a <- which(nms == "alpha")
      out$Vb <- V[which_b, which_b, drop = FALSE]
      rownames(out$Vb) <- colnames(out$Vb) <- names(out$coefficients)
      # covariance of the coefficients of the columns of Gfit
      keep <- c(which_b, which_a)
      to_gam_fit <- to_gam[gam_cols, , drop = FALSE]
      out$Vp <- as.matrix(to_gam_fit %*% V[keep, keep, drop = FALSE] %*% t(to_gam_fit))
      rownames(out$Vp) <- colnames(out$Vp) <- colnames(Gfit)
    }
  } else {
    out$hessian <- obj$he(opt$par)[fx_pars, fx_pars, drop = FALSE]
    rownames(out$hessian) <- colnames(out$hessian) <- colnames(Gfit)
    out$Vb <- try(solve(out$hessian))
    out$Vp <- out$Vb
  }

  out$Gsetup <- Gsetup
  out$convergence <- opt$convergence

  class(out) <- c("efp", "glm", "lm")
  out
}
