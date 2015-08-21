


.efEnv <- new.env()

.onLoad <- function(libname, pkgname) 
{
  packageStartupMessage("Building optimiser for first use...")

  tmp <- capture.output(
    stanmod <- rstan::stan_model(model_code = "
        data {
          int<lower=0> N; // number of observations
          int<lower=0> K; // number of parameters
          row_vector[N] y[3]; // data - 3 passes 
          matrix[N,K] A[3]; // the design matrices - 3 passes 
        }
        parameters {
          vector[K] alpha;
        } 
        model {
          vector[N] p[3]; // calculate all the probs required
          for (k in 1:3) {
            p[k] <- 1.0 ./ (1.0 + exp(-1.0 * A[k] * alpha)); 
          }
          increment_log_prob(y[1] * log(p[1])); 
          increment_log_prob(y[2] * log((1-p[1]) .* p[2]));
          increment_log_prob(y[3] * log((1-p[1]) .* (1-p[2]) .* p[3]));
          increment_log_prob(-1.0 * (y[1] + y[2] + y[3]) * log(p[1] + (1-p[1]) .* p[2] + (1-p[1]) .* (1-p[2]) .* p[3]));
        }")
  )
  assign("stanmod", stanmod, .efEnv)
}

##.onUnload <- function(libpath)
##    library.dynam.unload("lattice", libpath)

