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

  # for model matrices
  if (length(unique(id)) > 1) {
    Xs <- sparse.model.matrix(~ factor(id) - 1)
  } else {
    Xs <- sparse.model.matrix(~ id - 1)
  }
  Xp <- sparse.model.matrix(~ factor(pass) - 1)

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

  # define inputs for likelihood
  npasses <- Matrix::colSums(Xs)
  N <- ncol(Xs)
  K <- ncol(Gfit)
  s <- max(npasses)

  # get data in the correct order
  mat <-
    dplyr::left_join(
      expand.grid(pass = 1:s, id = sort(unique(id))),
      data.frame(
        pass = pass, id = id,
        y = Gsetup$y,
        offset = offset,
        irow = 1:length(pass)
      ),
      by = c("pass", "id")
    )
  mat[is.na(mat$y), c("y", "offset")] <- 0

  # we want observations in columns, passes in rows
  y <- matrix(mat$y, nrow = s, ncol = N)
  y_tot <- colSums(y)

  # same for offset
  offset <- matrix(mat$offset, nrow = s, ncol = N)

  # same for design matrix, but this is a bit trickier
  A <- array(0, c(s, N, K))
  for (i in 1:N) {
    iid <- sort(unique(id))[i]
    .irow <- subset(mat, id == iid)$irow
    A[!is.na(.irow),i,] <- Gfit[.irow[!is.na(.irow)],]
  }

  efdat <-
    list(
      N = N,
      K = K,
      S = s,
      npasses = npasses,
      y = y,
      yT = y_tot,
      offset = offset,
      A = A
    )

  if (!fit) {
    return(efdat)
  }

  obj <-
    MakeADFun(
      efdat,
      list(alpha = rep(0, K)),
      DLL = "ef",
      # map = map,
      inner.control = list(maxit = 500, trace = TRUE)
    )

  opt <-
    nlminb(
      obj$par,
      obj$fn,
      obj$gr,
      control = list(trace = 1)
    )


  # extract transformed parameters
  p <- opt$par[ grep("^p[[][0-9]*[,][0-9]*[]]", names(opt$par)) ]
  pi <- opt$par[ grep("^pi[[][0-9]*[,][0-9]*[]]", names(opt$par)) ]
  piT <- opt$par[ grep("^piT[[][0-9]*[]]", names(opt$par)) ]

  # convert to same order as input data
  #Xps <- Matrix::sparse.model.matrix( ~ factor(pass):factor(id) - 1)
  #opt $ p <- as.vector(Xps %*% p)
  #opt $ pi <- as.vector(Xps %*% pi)
  #opt $ piT <- data.frame(id = 1:N, piT = piT, total = y_tot)
  ordering <- as.numeric(factor(pass)) + (as.numeric(factor(id)) - 1) * s
  opt $ p <- p[ordering]
  opt $ pi <- pi[ordering]
  opt $ piT <- data.frame(id = 1:N, piT = piT, total = y_tot)

  # keep data for reffiting
  opt $ standat <- standat

  opt $ formula <- formula # for printing and summary
  opt $ llik <- opt $ value
  opt $ terms <- Gsetup $ terms
  opt $ call <- match.call()
  opt $ aic <- -2 * opt $ llik + 2 * ncol(Gfit)
  opt $ G <- G
  opt $ Gfit <- Gfit
  opt $ coefficients <- opt $ par[grep("alpha", names(opt$par))]
  names(opt $ coefficients) <- colnames(Gfit)
  opt $ df.null <- nrow(G)
  opt $ df.residual <- nrow(G) - ncol(Gfit)
  opt $ rank <- ncol(Gfit)
  opt $ fitted <- p <- transpar(opt $ coefficients, Gfit)

  opt $ residuals <- rep(0, nrow(data))
  opt $ null.deviance <- NA
  opt $ deviance <- NA
  opt $ family <- binomial()
  if (hessian) rownames(opt $ hessian) <- colnames(opt $ hessian) <- colnames(Gfit)
  opt $ Vb <- if (hessian) try(solve(-1 * opt $ hessian)) else NULL
  opt $ Gsetup <- Gsetup

  # get a gam container
  # g1 <- gam(G = Gsetup)
  # g1 $ coefficients[] <- opt $ par
  # g1 $ Vp[] <- opt $ Vb
  # g1 $ family <- binomial()
  # X <- predict(g1, type = "lpmatrix")
  # g1 $ linear.predictors <-  c(X %*% g1 $ coef)
  # g1 $ fitted.values <- c(1/(1 + exp(-g1 $ linear.predictors)))
  # g1 $ aic <- opt $ aic

  class(opt) <- c("efp", "glm", "lm")
  opt
}
