
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
#' @return glm type object
#' @export
efp <- function(formula, data = NULL, passes = NULL, verbose=TRUE, init = "0", hessian = FALSE) {

  if (!exists("stanmod")) {
    message("Building optimiser for first use...")
    stanmod <- rstan::stan_model(model_code = "
      data {
        int<lower=0> N; // number of observations
        int<lower=0> K; // number of parameters
        real S[N]; // the number of fishing passes
        real R[N]; // Zippins R (see seber p 312, eq 7.22)
        real T[N]; // total catches
        matrix[N,K] A;
      }
      parameters {
        vector[K] alpha;
      } 
      model {
        vector[N] expeta;
        expeta <- exp(A * alpha);
        for (i in 1:N) {
          real p;
          p <- expeta[i]/(1.0 + expeta[i]);
          increment_log_prob(T[i] * log(p));
          increment_log_prob(T[i] * R[i] * log(1-p));
          increment_log_prob(-T[i] * log(1 - (1-p)^S[i]) );
        }
      }")
    assign("stanmod", stanmod, .GlobalEnv)
  }

  if (is.null(data)) stop("must supply data")
  data0 <- subset(data, T > 0)

  if (is.null(passes)) stop("must supply the number of fishing runs")
  data0 $ S <- data0[[passes]]

  # set up model
  if (nrow(data0) == 1) {
    G <- matrix(1, 1, 1)
    Gsetup <- list()
  } else {
    Gsetup <- gam(formula, data = data0, fit = FALSE)
    G <- Gsetup $ X
  }

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


  S = data0 $ S
  T = data0 $ T
  R = with(data0, S - 1 - Z)

  if(nrow(data0) == 1) {
    dim(S) <- dim(T) <- dim(R) <- 1
  }

  standat <- 
    list(N = nrow(Gfit), K = ncol(Gfit), 
         S = S, T = T, R = R,
         A = Gfit)
  if (!verbose) {
    tmp <- 
      capture.output(
        opt <- rstan::optimizing(stanmod, data = standat, algorith = "BFGS", hessian = hessian, verbose = verbose, init = init)
      )
  } else {
    opt <- rstan::optimizing(stanmod, data = standat, algorith = "BFGS", hessian = hessian, verbose = verbose, init = init)
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
  opt $ residuals <- rep(0, nrow(data0))
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






#' Calculate the relavent statistics for fitting electrifishing models
#'
#' Basically, only three bits of info are needed from every electrofishing 
#' event: The total catch, T, and a kind of cumulative catch measure termed X
#' in Seber (1978)
#'
#' The data.frame must have the number of passes stored in a colum called 'Runs'
#' And fishing area must be in a row called 'Area'
#'
#' @param data a data.frame containing all relavent info.  
#' @param passnames a character vector giving the names of the columns in which
#'                  the catch data reside.  These must be ordered so that the 
#'                  first comes first, etc. 
#' @return a data frame
#' @export
getData <- function(data, passnames = paste0("S0_R", 1:6)) {
  catch <- as.matrix(data[passnames])
  rownames(catch) <- NULL
  s <- data $ Runs

  T <- sapply(1:nrow(catch), function(i) sum(catch[i,1:s[i]]))
  maxs <- ncol(catch)
  X <- colSums(sapply(s, function(s) c(s - 1:s, rep(0, maxs - s))) * t(catch), na.rm = TRUE)
  Z <- X / T
  phi <- 1 - Z/(s-1)

  data.frame(
      s = s,
      T = T,
      X = X,
      Z = ifelse(T > 0, Z, 0),
      phi = ifelse(T > 0, phi, 0)
    )
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




#' @export
summaryMods <- function(lst, m0 = NULL) {
  #aics <- sapply(lst, AIC)
  aics <- sapply(lst, BIC)

  tab <- 
   data.frame(
    forms = sapply(lst, function(x) paste(deparse(x$formula, width.cutoff = 500L))),
    aic = aics
    )

  if (!is.null(m0)) tab $ Daic <- tab $ aic - AIC(m0)
  tab <- tab[order(aics),]

  unique(tab)  
}

#' @export
getModels <- function(vars, n) {
  if (n > length(vars)) stop("n too big for number of variable")
  out <- do.call(expand.grid, lapply(1:n, function(i) 1:length(vars)))
  out <- unique(t(apply(out, 1, sort)))
  out <- out[apply(out, 1, function(x) !any(table(x)>1)),,drop=FALSE]
  if (n > 1) {
    apply(out, 1, function(x) paste(vars[x], collapse = " + "))
  } else {
    vars[out]
  }
}

