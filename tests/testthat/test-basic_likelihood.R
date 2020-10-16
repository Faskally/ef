
nloglik <- function(efdat, b) {
  p <- 1 / (1 + exp(-efdat$A %*% b - efdat$offset))

  dat <- data.frame(y = efdat$y, p)
  ll <-
    by(dat, efdat$sample_id, function(x) {
      pi <- x$p * c(1, cumprod(1 - x$p)[-nrow(x)])
      sum(x$y * log(pi)) - sum(x$y) * log(sum(pi))
    })
  -1 * sum(unlist(ll))
}

test_that("single sample, 3 pass, constant p", {
  # data with p = 0.5
  data <-
    data.frame(
      n = c(100, 50, 25),
      pass = 1:3,
      id = 1
    )

  formula <- n ~ 1

  # the data
  efdat <-
    efp(
      formula,
      data = data,
      fit = FALSE
    )

  obj <-
    MakeADFun(
      efdat,
      list(alpha = rep(0, ncol(efdat$A))),
      DLL = "ef",
      # map = map,
      inner.control = list(maxit = 500, trace = TRUE)
    )

  b <- 0
  tmb_nll <- obj$fn(b)
  r_nll <- nloglik(efdat, b)
  expect_equal(obj$fn(b), nloglik(efdat, b))
})
