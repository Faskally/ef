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

  # add random effects and penalised smooths
  for (block in efdat$blocks) {
    u <- alpha[block$idx]
    if (is.null(block$S)) {
      # independent random effects
      nll <- nll - sum(dnorm(u, 0, sigma[block$sig], log = TRUE))
    } else {
      # gaussian with sparse precision matrix sum_j S_j / sigma_j^2, with
      # one term for each penalty
      Q <- block$S[[1]] / sigma[block$sig[1]]^2
      for (j in seq_along(block$S)[-1]) {
        Q <- Q + block$S[[j]] / sigma[block$sig[j]]^2
      }
      nll <- nll - dgmrf(u, 0, Q, log = TRUE)
    }
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
