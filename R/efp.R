#' Estimate capture probabilites from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters
#'
#'
#' @param formula a formula object
#' @param data a data.frame containing all relavent info
#' @param pass a vector of integers giving the pass number of the observation
#' @param id a vector of integers identifying an observation (a set of electrofishing passes)
#' @param offset an possible offset for the linear predictor of capture probability
#' @param verbose if TRUE stan optimiser messages are printed to the screen
#' @param init should initialisatiom be random?
#' @param hessian if TRUE the hessian is computed and the covariance matrix of the parameters is returned via Vb
#' @param fit if TRUE model is fitted if FALSE the data that would be passed to the optimiser is returned
#' @param sample_re should sample random effects be included
#' @return glm type object
#' @examples
#' \dontrun{
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
#' out <- cbind(ef_data, p = fitted(m2, type = "p"), pi = fitted(m2, type = "pi"))
#' out
#'
#' # to get offset for fitting density do this
#' m2$piT
#'
#' msim2 <- simulate(m2, nsim = 1)
#' }
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb
#' @importFrom Matrix sparse.model.matrix
#' @importFrom mgcv gam
#' @importFrom methods as
#'
#' @export
efp <- function(formula, data = NULL, pass = pass, id = id,
  offset = NULL, verbose = FALSE, init = "0", hessian = TRUE,
  fit = TRUE, sample_re = FALSE) {

  # some checks
  if (is.null(data)) stop("must supply data")

  # get pass information
  pass <- substitute(pass)
  pass <- eval(pass, data, environment(formula))
  # get id
  id <- substitute(id)
  id <- eval(id, data, environment(formula))

  # set up offset
  if (is.null(offset)) {
    offset <- rep(0, nrow(data))
  }

  # set up model
  Gsetup <- gam(formula, data = data, fit = FALSE, drop.unused.levels = FALSE)
  G <- Gsetup$X
  colnames(G) <- Gsetup$term.names
  if (nrow(G) != nrow(data)) stop("I think there are NAs in your data, please check and remove them before rerunning.")

  # remove redundant / collinear parameters
  qr.G <- qr(G)
  rank.deficient <- qr.G $ pivot[abs(diag(qr.G $ qr)) < 1e-7]
  if (length(rank.deficient)) {
    droppar <- paste(colnames(G)[rank.deficient], collapse = "\n\t")
    warning("*** Model has ", length(rank.deficient)," too many parameter(s)!!\n    i will remove the redundant ones:\n\t", droppar, call. = FALSE)
    Gfit <- G[,-rank.deficient]
  } else {
    Gfit <- G
  }

  # get data for likelihood
  ord <- order(id, pass)

  X <- Gfit[ord, , drop = FALSE]
  X <- as(X, "dgTMatrix")
  if (sample_re) {
    Z <- sparse.model.matrix(~ factor(id[ord]) - 1)
  } else {
    Z <- matrix(1, nrow(X), 1)
  }
  Z <- as(Z, "dgTMatrix")

  efdat <-
    list(
      y = Gsetup$y[ord],
      offset = offset[ord],
      sample_id = factor(id[ord]),
      X = X,
      Z = Z,
      sample_re = as.integer(sample_re)
    )

  if (!fit) {
    return(efdat)
  }

  params <-
    list(
      alpha = rep(0, ncol(efdat$Z)),
      beta = rep(0, ncol(efdat$X)),
      log_sigma = 0
    )

  if (sample_re == 0L) {
    rand <- NULL
    map <- list(alpha = factor(NA), log_sigma = factor(NA))
  } else {
    rand <- "alpha"
    map <- list()
  }

  obj <-
    MakeADFun(
      efdat,
      params,
      random = rand,
      DLL = "ef",
      map = map,
      hessian = hessian,
      silent = !verbose
    )

  if (verbose) {
    control = list(trace = 1)
  } else {
    control = list()
  }

  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)

  # set up output object
  out <- list()

  # extract transformed parameters
  out$p <- unlist(obj$report()$ps)

  fx_pars <- grep("beta", names(opt$par))
  out$sdrep <- sdreport(obj)
  out$rep <- obj$report()
  # keep data for reffiting
  out$data <- data
  out$formula <- formula # for printing and summary
  out$llik <- -1 * obj$fn(opt$par)
  out$terms <- Gsetup$terms
  out$call <- match.call()
  out$aic <- -2 * out$llik + 2*ncol(Gfit)
  out$G <- G
  out$Gfit <- Gfit
  out$coefficients <- opt$par[fx_pars]
  names(out$coefficients) <- colnames(Gfit)
  out$df.null <- nrow(G)
  out$df.residual <- nrow(G) - ncol(Gfit)
  out$rank <- ncol(Gfit)
  out$fitted <- out$p

  out$residuals <- rep(0, nrow(data))
  out$null.deviance <- NA
  out$deviance <- NA
  out$family <- binomial()
  if (sample_re | !hessian) {
    out$hessian <- NULL
    out$Vb <- NULL
  } else {
    out$hessian <- obj$he(opt$par)[fx_pars, fx_pars, drop = FALSE]
    rownames(out$hessian) <- colnames(out$hessian) <- colnames(Gfit)
    out$Vb <- try(solve(out$hessian))
  }

  out$Gsetup <- Gsetup

  class(out) <- c("efp", "glm", "lm")
  out
}
