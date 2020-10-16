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
#' library(rstan)
#' plot(msim2, pars = "p")
#' }
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb
#' @importFrom Matrix sparse.model.matrix
#' @importFrom mgcv gam
#'
#' @export
efp <- function(formula, data = NULL, pass = pass, id = id, offset = NULL,
                verbose = FALSE, init = "0", hessian = TRUE, fit = TRUE) {

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

  efdat <-
    list(
      y = Gsetup$y[ord],
      offset = offset[ord],
      sample_id = factor(id[ord]),
      A = Gfit[ord, , drop = FALSE]
    )

  if (!fit) {
    return(efdat)
  }

  obj <-
    MakeADFun(
      efdat,
      list(alpha = rep(0, ncol(Gfit))),
      DLL = "ef",
      # map = map,
      inner.control = list(maxit = 500, trace = TRUE)
    )

  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))

  # set up output object
  out <- list()

  out$tmb_report <- obj$report()

  # extract transformed parameters
  #p <- opt$par[ grep("^p[[][0-9]*[,][0-9]*[]]", names(opt$par)) ]
  #pi <- opt$par[ grep("^pi[[][0-9]*[,][0-9]*[]]", names(opt$par)) ]
  #piT <- opt$par[ grep("^piT[[][0-9]*[]]", names(opt$par)) ]

  # convert to same order as input data
  #ordering <- as.numeric(factor(pass)) + (as.numeric(factor(id)) - 1) * s
  #out$p <- p[ordering]
  #out$pi <- pi[ordering]
  #out$piT <- data.frame(id = 1:N, piT = piT, total = y_tot)

  # keep data for reffiting
  out$data <- data
  out$formula <- formula # for printing and summary
  out$llik <- -1 * obj$fn()
  out$terms <- Gsetup$terms
  out$call <- match.call()
  out$aic <- -2 * out$llik + 2*ncol(Gfit)
  out$G <- G
  out$Gfit <- Gfit
  out$coefficients <- opt$par
  names(out$coefficients) <- colnames(Gfit)
  out$df.null <- nrow(G)
  out$df.residual <- nrow(G) - ncol(Gfit)
  out$rank <- ncol(Gfit)
  out$fitted <- p <- transpar(out$coefficients, Gfit)

  out$residuals <- rep(0, nrow(data))
  out$null.deviance <- NA
  out$deviance <- NA
  out$family <- binomial()
  out$hessian <- obj$he()
  if (hessian) rownames(out$hessian) <- colnames(out$hessian) <- colnames(Gfit)
  out$Vb <- if (hessian) try(solve(-1 * out$hessian)) else NULL
  out$Gsetup <- Gsetup

  class(out) <- c("efp", "glm", "lm")
  out
}
