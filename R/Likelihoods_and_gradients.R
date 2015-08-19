
#
# Purpose: fit capture probabilities and abundances to electrofishing data
# 
# author: CP Millar, millarc@marlab.ac.uk
# origional date: 12/2014
#

#' Estimate capture probabilites from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters 
#' 
#'
#' @param formula a formula object
#' @param data a data.frame containing all relavent info
#' @param pass a vector of integers giving the pass number of the observation
#' @param verbose if TRUE stan optimiser messages are printed to the screen
#' @param init should initialisatiom be random?
#' @param hessian if TRUE the hessian is computed and the covariance matrix of the parameters is returned via Vb
#' @return glm type object
#' @export
efp <- function(formula, data = NULL, pass, verbose=TRUE, init = "0", hessian = FALSE) {

  # some checks
  if (is.null(data)) stop("must supply data")

  # get pass information
  if (missing(pass)) stop("must supply pass number")
  pass <- substitute(pass)
  pass <- eval(pass, data, environment(formula))
  pass <- as.integer(pass)
  if (length(unique(pass)) != 3 || !all(sort(unique(pass)) == 1:3)) stop("There should only be 3 passes and they should be numbered 1 to 3")
  # the within sample structure is then,
  X <- model.matrix(~ factor(pass) - 1)

  # set up model
  Gsetup <- gam(formula, data = data, fit = FALSE)
  G <- Gsetup $ X
 
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
  # need to check that data has rows of multiples of 3
  N <- as.integer(nrow(Gfit) / 3)
  K <- ncol(Gfit)
  # get data in the correct order
  y <- Gsetup $ y[row(X)[as.logical(X)]]
  dim(y) <- c(N, 3)
  y <- aperm(y, c(2,1))
  # same for design matrix, but this is a bit trickier
  A <- t(Gfit[row(X)[as.logical(X)],])
  dim(A) <- c(K, N, 3)
  A <- aperm(A, c(3,2,1))

  standat <- 
    list(N = N,  
         K = K,
         y = y,
         A = A)

  if (!verbose) {
    tmp <- 
      capture.output(
        opt <- rstan::optimizing(.efEnv $ stanmod, data = standat, algorith = "BFGS", hessian = hessian, verbose = verbose, init = init)
      )
  } else {
    opt <- rstan::optimizing(.efEnv $ stanmod, data = standat, algorith = "BFGS", hessian = hessian, verbose = verbose, init = init)
  } 

  opt $ formula <- formula # for printing and summary
  opt $ llik <- opt $ value
  opt $ terms <- Gsetup $ terms
  opt $ call <- match.call()
  opt $ aic <- -2 * opt $ llik + 2 * ncol(Gfit)
  opt $ G <- G
  opt $ Gfit <- Gfit
  opt $ coefficients <- opt $ par
  names(opt $ coefficients) <- colnames(Gfit)
  opt $ df.null <- nrow(G)
  opt $ df.residual <- nrow(G) - ncol(Gfit)
  opt $ rank <- ncol(Gfit)
  opt $ fitted <- p <- transpar(opt $ par, Gfit)
  opt $ residuals <- rep(0, nrow(data))
  opt $ null.deviance <- NA
  opt $ deviance <- NA 
  opt $ family <- binomial()
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





#' Utility function to convert parameters to probabilities 
#'
#' The matrix G shoudl be of dimension n x p, 
#' and the parameter vector should be lenght p
#'
#' @param par fitted model parameters
#' @param G The design matrix for a model 
#' @return a data frame
#' @export
transpar <- function(par, G) {
   1/(1 + exp(-c(G %*% par)))
}

