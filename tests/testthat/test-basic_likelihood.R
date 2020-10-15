
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
      data = ef_data,
      fit = FALSE
    )

  obj <-
    MakeADFun(
      efdat,
      list(alpha = rep(0, K)),
      DLL = "ef",
      # map = map,
      inner.control = list(maxit = 500, trace = TRUE)
    )

  llik_tmb <- obj$fn(0)
  llik <- sum(c(efdat$y) * log(0.5)^(1:3)) - efdat$yT * log(sum(0.5^(1:3)))
  testthat::expect_equal(llik_tmb, llik)
})
