# log likelihood of a single multipass electrofishing sample
ll_sample <- function(p, y) {
  npasses <- length(y)
  pi <- p
  for (i in seq_len(npasses)[-1]) {
    for (j in 1:(i - 1)) {
      pi[i] <- pi[i] * (1 - p[j])
    }
  }
  pos <- which(y > 0)
  ll <- sum(y[pos] * log(pi[pos]))
  ll - sum(y) * log(sum(pi))
}

# negative log likelihood, the RTMB objective function
ef_nll <- function(params, efdat) {
  alpha <- params$alpha
  beta <- params$beta
  log_sigma <- params$log_sigma

  eta <- drop(efdat$X %*% beta) + efdat$offset
  if (length(alpha)) {
    # add random effects
    eta <- eta + drop(efdat$Z %*% alpha)
  }

  p <- 1 / (1 + exp(-eta))
  idx <- split(seq_along(efdat$y), efdat$sample_id)
  ps <- unname(lapply(idx, function(i) p[i]))
  ys <- lapply(idx, function(i) efdat$y[i])

  nll <- 0 # negative log likelihood

  sigma <- exp(log_sigma)

  if (length(alpha)) {
    # add random effects, one sigma per random effect term
    nll <- nll - sum(dnorm(alpha, 0, sigma[efdat$re_id], log = TRUE))
  }

  for (i in seq_along(idx)) {
    nll <- nll - ll_sample(ps[[i]], ys[[i]])
  }

  # delta method for sigma
  if (length(sigma)) {
    ADREPORT(sigma)
  }
  REPORT(sigma)

  REPORT(ps)

  nll
}
