#' @export
simulate.efp <- function(object, nsim = 1000, seed = NULL, ...) {
  if (nsim < 1000) message("Increasing nsim to 1000.")
  nsim <- max(nsim, 1000)
  rstan::sampling(.efEnv$stanmod2, data = object$standat, iter = 1000 + nsim, warmup = 1000, chains = 1, ...)
}
