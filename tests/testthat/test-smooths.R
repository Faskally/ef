# simulate removal data with a smooth trend in capture probability
sim_trend <- function(nyear = 12, nsite = 4, seed = 1) {
  set.seed(seed)
  data <- expand.grid(pass = 1:3, site = factor(1:nsite), year = 1:nyear)
  data$visit <- factor(paste(data$site, data$year))
  data$fyear <- factor(data$year, levels = 1:nyear)
  data$eta <- 0.3 + 0.8 * sin(2 * pi * data$year / nyear)
  p <- 1 / (1 + exp(-data$eta))
  N <- rep(150, nrow(data) / 3)
  data$n <- NA
  for (i in 1:3) {
    ipass <- which(data$pass == i)
    data$n[ipass] <- rbinom(length(ipass), N, p[ipass])
    N <- N - data$n[ipass]
  }
  data
}

test_that("random walk penalties", {
  rw1 <- rw_penalty(1:5)
  rw2 <- rw_penalty(letters[1:5], order = 2)

  expect_equal(dim(rw1), c(5, 5))
  expect_equal(rownames(rw2), letters[1:5])
  expect_equal(unname(rw1[1, ]), c(1, -1, 0, 0, 0))
  expect_equal(unname(rw2[3, ]), c(1, -4, 6, -4, 1))

  # constant and linear functions are not penalised
  expect_equal(drop(rw1 %*% rep(1, 5)), rep(0, 5), check.attributes = FALSE)
  expect_equal(drop(rw2 %*% 1:5), rep(0, 5), check.attributes = FALSE)
  expect_equal(max(abs(rw1 %*% attr(rw1, "null.space"))), 0)
  expect_equal(max(abs(rw2 %*% attr(rw2, "null.space"))), 0)
  expect_equal(dim(attr(rw2, "null.space")), c(5, 2))

  expect_error(rw_penalty(1:5, order = 3), "order")
  expect_error(rw_penalty(c(1, 1, 2)), "unique")
  expect_error(rw_penalty(1:2, order = 2), "more than")
})

test_that("random walk likelihood", {
  # a random walk is the cumulative sum of independent steps, so the model
  # can also be written using independent random effects, with a dense
  # design matrix. The marginal likelihoods should be identical.
  data <- sim_trend()

  for (order in 1:2) {
    rw <- rw_penalty(1:12, order = order)
    efdat <-
      efp(
        n ~ s(fyear, bs = "mrf", xt = list(penalty = rw)),
        data = data, pass = pass, id = visit, fit = FALSE
      )

    # the starting values of the walk are fixed effects
    expect_equal(ncol(efdat$X), order)
    expect_equal(ncol(efdat$Z), 12 - order)
    expect_equal(colnames(efdat$Z), paste0("s(fyear).", (order + 1):12))
    # the penalty is sparse
    expect_s4_class(efdat$blocks[[1]]$S[[1]], "sparseMatrix")
    expect_lte(length(efdat$blocks[[1]]$S[[1]]@x), (2 * order + 1) * 12)
    expect_equal(efdat$blocks[[1]]$logdet, 0)

    k <- ncol(efdat$Z)
    L <- 1 * lower.tri(diag(k), diag = TRUE)
    if (order == 2) L <- L %*% L
    efdat_steps <- efdat
    efdat_steps$Z <- as(as.matrix(efdat$Z) %*% L, "TsparseMatrix")
    efdat_steps$blocks[[1]]["S"] <- list(NULL)

    # joint likelihood
    set.seed(order)
    steps <- rnorm(k, 0, 0.2)
    beta <- rnorm(order, 0, 0.1)
    expect_equal(
      ef:::ef_nll(list(alpha = drop(L %*% steps), beta = beta, log_sigma = -1.5), efdat),
      ef:::ef_nll(list(alpha = steps, beta = beta, log_sigma = -1.5), efdat_steps)
    )

    # marginal likelihood
    params <- list(alpha = rep(0, k), beta = beta, log_sigma = -1.5)
    obj <-
      RTMB::MakeADFun(
        function(params) ef:::ef_nll(params, efdat), params,
        random = "alpha", silent = TRUE
      )
    obj_steps <-
      RTMB::MakeADFun(
        function(params) ef:::ef_nll(params, efdat_steps), params,
        random = "alpha", silent = TRUE
      )
    expect_equal(obj$fn(obj$par), obj_steps$fn(obj_steps$par), tolerance = 1e-8)
  }
})

test_that("random walk smooths fit", {
  data <- sim_trend()
  rw1 <- rw_penalty(1:12)
  rw2 <- rw_penalty(1:12, order = 2)

  m0 <- efp(n ~ 1, data = data, pass = pass, id = visit)
  m1 <- efp(n ~ s(fyear, bs = "mrf", xt = list(penalty = rw1)),
    data = data, pass = pass, id = visit
  )
  m2 <- efp(n ~ s(fyear, bs = "mrf", xt = list(penalty = rw2)),
    data = data, pass = pass, id = visit
  )

  expect_equal(m1$convergence, 0)
  expect_equal(m2$convergence, 0)

  # the unpenalised part of the smooth is a fixed effect
  expect_equal(names(coef(m1)), "(Intercept)")
  expect_equal(names(coef(m2)), c("(Intercept)", "s(fyear).null1"))
  expect_equal(names(m1$sigma), "s(fyear)")
  expect_true(all(c(m1$sigma, m2$sigma) > 0))
  expect_equal(length(m1$random), 11)
  expect_equal(length(m2$random), 10)
  expect_equal(attr(logLik(m1), "df"), 2)
  expect_equal(attr(logLik(m2), "df"), 3)

  # there is a strong trend in the data
  expect_gt(as.numeric(logLik(m1)), as.numeric(logLik(m0)) + 10)
  expect_lt(sqrt(mean((predict(m1) - data$eta)^2)), 0.15)
  expect_lt(sqrt(mean((predict(m2) - data$eta)^2)), 0.15)

  # predictions agree with fitted values
  for (m in list(m1, m2)) {
    pred <- predict(m, newdata = data, type = "response", se.fit = TRUE)
    expect_equal(as.vector(pred$fit), fitted(m))
    expect_true(all(pred$se.fit > 0 & pred$se.fit < 0.1))
    expect_equal(dim(m$Vp), rep(ncol(m$Gfit), 2))

    m_gam <- suppressMessages(as.gam(m))
    expect_equal(
      as.vector(predict(m_gam, newdata = data, type = "response")),
      fitted(m)
    )
  }

  expect_output(print(m1), "sigma s\\(fyear\\)")
})

test_that("random walk smooths with levels in any order", {
  # the walk follows the order of the levels in the penalty, which is not
  # the alphabetical order used by mgcv. Year 5 has no data.
  data <- sim_trend()
  data <- data[data$year != 5, ]
  data$visit <- droplevels(data$visit)
  data$fyear <- factor(data$year + 5, levels = 6:17)
  rw1 <- rw_penalty(6:17)

  efdat <-
    efp(
      n ~ s(fyear, bs = "mrf", xt = list(penalty = rw1)),
      data = data, pass = pass, id = visit, fit = FALSE
    )
  expect_equal(colnames(efdat$Z), paste0("s(fyear).", 7:17))
  expect_equal(
    unname(as.matrix(efdat$blocks[[1]]$S[[1]])),
    unname(rw1[-1, -1]),
    check.attributes = FALSE
  )

  m1 <- efp(n ~ s(fyear, bs = "mrf", xt = list(penalty = rw1)),
    data = data, pass = pass, id = visit
  )
  expect_equal(m1$convergence, 0)
  expect_equal(
    as.vector(predict(m1, newdata = data, type = "response")),
    fitted(m1)
  )

  # the level with no data is interpolated
  newdata <- data.frame(fyear = factor(9:11, levels = 6:17))
  pred <- predict(m1, newdata = newdata, se.fit = TRUE)
  expect_equal(unname(pred$fit[2]), unname(mean(pred$fit[c(1, 3)])), tolerance = 1e-6)
  expect_gt(pred$se.fit[2], max(pred$se.fit[c(1, 3)]))
})

test_that("other penalised smooths fit", {
  data <- sim_trend()
  set.seed(2)
  data$x <- rep(runif(nlevels(data$visit)), each = 3)

  # thin plate spline, sparse and dense set up
  m1 <- efp(n ~ s(year, k = 6), data = data, pass = pass, id = visit)
  m1_dense <- efp(n ~ s(year, k = 6),
    data = data, pass = pass, id = visit, sparse = FALSE
  )
  expect_equal(m1$convergence, 0)
  expect_equal(names(coef(m1)), c("(Intercept)", "s(year).null1"))
  expect_equal(as.vector(predict(m1, newdata = data, type = "response")), fitted(m1))
  expect_equal(
    as.vector(predict(m1_dense, newdata = data, type = "response")),
    fitted(m1_dense)
  )
  # these are not identical, as the unpenalised part of the smooth is
  # separated in different ways, but should be very similar
  expect_equal(fitted(m1), fitted(m1_dense), tolerance = 0.01)
  expect_equal(unname(m1$sigma), unname(m1_dense$sigma), tolerance = 0.1)

  # tensor product smooth, which has two penalties
  m2 <- efp(n ~ te(year, x, k = 3), data = data, pass = pass, id = visit)
  expect_equal(m2$convergence, 0)
  expect_equal(names(m2$sigma), c("te(year,x)1", "te(year,x)2"))
  expect_equal(as.vector(predict(m2, newdata = data, type = "response")), fitted(m2))

  # smooths and random effects together
  m3 <- efp(n ~ s(year, bs = "cr", k = 6) + s(site, bs = "re"),
    data = data, pass = pass, id = visit
  )
  expect_equal(names(m3$sigma), c("s(year)", "s(site)"))
  expect_equal(as.vector(predict(m3, newdata = data, type = "response")), fitted(m3))

  # by variable smooths use the mgcv constrained basis
  m4 <- efp(n ~ site + s(year, by = site, k = 4),
    data = data, pass = pass, id = visit
  )
  expect_equal(length(m4$sigma), nlevels(data$site))
  expect_equal(as.vector(predict(m4, newdata = data, type = "response")), fitted(m4))
})
