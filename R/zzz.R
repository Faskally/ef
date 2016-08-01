


.efEnv <- new.env()


.buildOptimiser <- function(force = FALSE) {

  if (exists("stanmod", envir = .efEnv) & !force) return (invisible(FALSE))

  message("Building optimiser for first use...")

  stanmod <- rstan::stan_model(model_code = "
      data {
        int<lower=0> N; // number of observations
        int<lower=0> K; // number of parameters
        row_vector[N] y[3]; // data - 3 passes
        vector[N] offset[3]; // offset - 3 passes
        matrix[N,K] A[3]; // the design matrices - 3 passes
      }
      parameters {
        vector[K] alpha;
      }
      model {
        vector[N] p[3]; // calculate all the probs required
        for (s in 1:3) {
          p[s] = 1.0 ./ (1.0 + exp(-1.0 * A[s] * alpha - offset[s]));
        }
        target += y[1] * log(p[1]);
        target += y[2] * log((1-p[1]) .* p[2]);
        target += y[3] * log((1-p[1]) .* (1-p[2]) .* p[3]);
        target += -1.0 * (y[1] + y[2] + y[3]) * log(p[1] + (1-p[1]) .* p[2] + (1-p[1]) .* (1-p[2]) .* p[3]);
      }")

  assign("stanmod", stanmod, envir = .efEnv)

  invisible(TRUE)
}


.buildOptimiser2 <- function(force = FALSE) {

  if (exists("stanmod2", envir = .efEnv) & !force) return (invisible(FALSE))

  message("Building optimiser for first use...")

  stanmod2 <- rstan::stan_model(model_code = "
    data {
      int<lower=0> N; // number of multipass observations
      int<lower=0> K; // number of parameters
      int<lower=0> s; // max number of passes
      int<lower=0> npasses[N]; // number of passes per observation - vectors of s passes
      row_vector[N] y[s]; // data - vectors of s passes
      row_vector[N] yT; // data - total catch
      vector[N] offset[s]; // offset - s passes
      matrix[N,K] A[s]; // the design matrices - one for each pass
    }
    parameters {
       vector[K] alpha;
    }
    transformed parameters {
      vector[N] p[s];
      vector[N] pi[s];
      vector[N] piT;
      // calculate all the within pass probs required
      for (i in 1:s) {
        p[i] = 1.0 ./ (1.0 + exp(-1.0 * A[i] * alpha - offset[i]));
      }
      // calculate all the marginal(?) probs required
      pi[1] = p[1];
      for (i in 2:s) {
        pi[i] = p[i];
        for (j in 1:(i-1)) {
          pi[i] = pi[i] .* (1-p[j]);
        }
      }
      // calculate the probability of capture - needs number of passes
      piT = pi[1];
      for (i in 1:N) {
        for (j in 2:npasses[i]) {
          piT[i] = piT[i] + pi[j,i];
        }
      }
    }
    model {
      // calculate the probability of capture - needs number of passes
      for (i in 1:s) {
        target += y[i] * log(pi[i]);
      }
      target += -1.0 * yT * log(piT);
    }")

  assign("stanmod2", stanmod2, envir = .efEnv)

  invisible(TRUE)
}





.buildOptimiser3 <- function(force = FALSE) {

  if (exists("stanmod3", envir = .efEnv) & !force) return (invisible(FALSE))

  message("Building optimiser 3 for first use...")

  code <-
  "data {
    int<lower=0> N; // number of multipass observations
    int<lower=0> K; // number of parameters
    int<lower=0> s; // max number of passes
    int<lower=0> npasses[N]; // number of passes per observation - vectors of s passes
    row_vector[N] y[s]; // data - vectors of s passes
    row_vector[N] yT; // data - total catch
    vector[N] offset[s]; // offset - s passes
    matrix[N,K] A[s]; // the design matrices - one for each pass
    matrix[K,K] Q; // penalty matrix
  }
  parameters {
    vector[K] alpha;
  }
  transformed parameters {
    vector[N] p[s];
    vector[N] pi[s];
    vector[N] piT;
    // calculate all the within pass probs required
    for (i in 1:s) {
      p[i] <- 1.0 ./ (1.0 + exp(-1.0 * A[i] * alpha - offset[i]));
    }
    // calculate all the marginal(?) probs required
    pi[1] = p[1];
    for (i in 2:s) {
      pi[i] = p[i];
      for (j in 1:(i-1)) {
        pi[i] = pi[i] .* (1-p[j]);
      }
    }
    // calculate the probability of capture - needs number of passes
    piT = pi[1];
    for (i in 1:N) {
      for (j in 2:npasses[i]) {
        piT[i] = piT[i] + pi[j,i];
      }
    }
  }
  model {
    // calculate the probability of capture - needs number of passes
    for (i in 1:s) {
      target += y[i] * log(pi[i]);
    }
    target += -1.0 * yT * log(piT);
    // add in penalty
    target += -1.0 * quad_form(Q, alpha);
  }"

  stanmod3 <- rstan::stan_model(model_code = code)

  assign("stanmod3", stanmod3, envir = .efEnv)

  invisible(TRUE)
}





## these optimisers us tmb

.buildTMB_dll <- function() {


  message("Building optimiser for first use...")



  invisible(TRUE)
}











##.onLoad <- function(libname, pkgname)
##{
##  #packageStartupMessage("Building optimiser for first use...")
##
##}

##.onUnload <- function(libpath)
##    library.dynam.unload("lattice", libpath)

