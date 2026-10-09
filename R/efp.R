#' Estimate capture probabilites from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters
#'
#'
#' @param formula a formula object. Random effects are specified using
#'   \code{s(f, bs = "re")} terms, where \code{f} is a factor, see Details.
#' @param data a data.frame containing all relavent info
#' @param pass a vector of integers giving the pass number of the observation
#' @param id a vector of integers identifying an observation (a set of electrofishing passes)
#' @param offset an possible offset for the linear predictor of capture probability
#' @param verbose if TRUE optimiser messages are printed to the screen
#' @param init should initialisatiom be random?
#' @param hessian if TRUE the hessian is computed and the covariance matrix of the parameters is returned via Vb
#' @param fit if TRUE model is fitted if FALSE the data that would be passed to the optimiser is returned
#' @param control a list of control commands to be passed to optim
#' @return glm type object
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
#' # a smooth effect of a continuous covariate. Smooths are unpenalised, so
#' # the flexibility is set by k
#' data(ef_data)
#' ef_data$sampleID <- factor(paste(ef_data$date, ef_data$siteID,
#'                                  ef_data$lifestage, ef_data$species))
#' ef_data$pass12 <- factor(replace(ef_data$pass, ef_data$pass > 2, 2))
#' ef_data$nyear <- as.numeric(as.character(ef_data$year))
#' m5 <- efp(count ~ siteID + pass12 + lifestage + s(nyear, k = 3),
#'           data = ef_data, pass = pass, id = sampleID)
#' m5
#'
#' # predicted capture probability over the smooth, with standard errors
#' newdata <- data.frame(siteID = "Site1", pass12 = "1", lifestage = "Fry",
#'                       nyear = seq(2011, 2015, by = 0.5))
#' pred <- predict(m5, newdata = newdata, type = "response", se.fit = TRUE)
#' cbind(newdata, fit = pred$fit, se = pred$se.fit)
#'
#' @importFrom RTMB MakeADFun sdreport ADREPORT REPORT dnorm
#' @importFrom mgcv gam
#' @importFrom methods as
#' @importClassesFrom Matrix TsparseMatrix
#' @importFrom stats binomial
#' @importFrom stats optim
#' @importFrom stats ave
#'
#' @export
efp <- function(formula, data = NULL, pass = pass, id = id,
  offset = NULL, verbose = FALSE, init = "0", hessian = TRUE,
  fit = TRUE, control = list()) {

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
  re_labels <- vapply(Gsetup$smooth[is_re], function(x) x$label, "")
  re_cols <- lapply(Gsetup$smooth[is_re], function(x) x$first.para:x$last.para)
  re_id <- rep(seq_along(re_cols), lengths(re_cols))
  re_cols <- as.integer(unlist(re_cols))
  fx_cols <- setdiff(seq_len(ncol(G)), re_cols)

  # remove redundant / collinear parameters (fixed effects only)
  qr.G <- qr(G[, fx_cols, drop = FALSE])
  rank.deficient <- fx_cols[qr.G $ pivot[abs(diag(qr.G $ qr)) < 1e-7]]
  if (length(rank.deficient)) {
    droppar <- paste(colnames(G)[rank.deficient], collapse = "\n\t")
    warning("*** Model has ", length(rank.deficient)," too many parameter(s)!!\n    i will remove the redundant ones:\n\t", droppar, call. = FALSE)
    fx_cols <- setdiff(fx_cols, rank.deficient)
  }
  Gfit <- G[, sort(c(fx_cols, re_cols)), drop = FALSE]

  # get data for likelihood
  ord_output <- order(id, pass)

  efdat <-
    list(
      y = Gsetup$y,
      offset = offset,
      sample_id = factor(id),
      X = as(G[, fx_cols, drop = FALSE], "TsparseMatrix"),
      Z = as(G[, re_cols, drop = FALSE], "TsparseMatrix"),
      re_id = re_id
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

  if (length(re_cols)) {
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

#  if (verbose) {
#    control = list(trace = 1)
#  } else {
#    control = list()
#  }

  #opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
  
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
  out$p <- unlist(out$rep$ps)
  # Re-order to match input data
  out$p <- out$p[order(ord_output)]
  # probability of capture in each pass, allowing for the fish removed in
  # previous passes
  out$pi <- as.vector(
    ave(out$p, factor(id), FUN = function(x) x * c(1, cumprod(1 - x)[-length(x)]))
  )
  
  fx_pars <- grep("beta", names(opt$par))
  out$sdrep <- sdreport(obj, getJointPrecision = hessian && length(re_cols) > 0)
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
  out$re_labels <- re_labels
  out$coefficients <- opt$par[fx_pars]
  names(out$coefficients) <- colnames(G)[fx_cols]
  out$random <- par_full[names(par_full) == "alpha"]
  names(out$random) <- colnames(G)[re_cols]
  out$sigma <- exp(opt$par[names(opt$par) == "log_sigma"])
  names(out$sigma) <- re_labels
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
  } else if (length(re_cols)) {
    # covariance of fixed and random effects, in the column order of Gfit
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
      keep <- c(which_b, which_a)[order(c(fx_cols, re_cols))]
      out$Vp <- V[keep, keep, drop = FALSE]
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
