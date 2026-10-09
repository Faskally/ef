library(testthat)

nloglik <- function(efdat, b) {
  p <- 1 / (1 + exp(-efdat$X %*% b - efdat$offset))

  dat <- data.frame(y = efdat$y, p[,1])
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

  params <-
    list(
      alpha = rep(0, ncol(efdat$Z)),
      beta = rep(0, ncol(efdat$X)),
      log_sigma = rep(0, max(c(0, efdat$re_id)))
    )

  obj <-
    RTMB::MakeADFun(
      function(params) ef:::ef_nll(params, efdat),
      params,
      # map = map,
      inner.control = list(maxit = 500, trace = TRUE)
    )

  tmb_nll <- obj$fn(obj$par)
  r_nll <- nloglik(efdat, obj$par["beta"])
  expect_equal(tmb_nll, r_nll)
})


test_that("simple model fits", {

  data(ef_data)

  ef_data$pass12 <- factor(replace(ef_data$pass, ef_data$pass > 2, 2))

  ef_data$sampleID <-
    as.numeric(factor(paste(
      ef_data$date, ef_data$siteID,
      ef_data$lifestage, ef_data$species
    )))

  ef_data$visitID <- as.numeric(factor(paste(ef_data$date, ef_data$siteID)))


  m1 <- efp(count ~ 1,
    data = ef_data, pass = pass, id = sampleID
  )

  m2 <- efp(count ~ lifestage + pass12,
    data = ef_data, pass = pass, id = sampleID
  )

  m3 <- efp(count ~ lifestage + pass + year,
    data = ef_data, pass = pass, id = sampleID
  )

  expect_equivalent(coef(m1), 0.9734342, tolerance = 1e-5)
  expect_equivalent(as.numeric(logLik(m1)), -1612.81836111119, tolerance = 1e-9)

  expect_equivalent(coef(m2), c(0.6666140, 0.6450585, -0.2563096), tolerance = 1e-4)
  expect_equivalent(as.numeric(logLik(m2)), -1592.30666991703, tolerance = 1e-9)

  expect_equivalent(
    coef(m3),
    c(0.67592094, 1.07051329, -0.80160463, -0.39612386, 0.19011401, 0.67211445, 0.02518319),
    tolerance = 1e-4
  )
  expect_equivalent(as.numeric(logLik(m3)), -1584.7000365589, tolerance = 1e-9)
})

test_that("random effects likelihood", {
  data(ef_data)
  ef_data$sampleID <-
    factor(paste(
      ef_data$date, ef_data$siteID,
      ef_data$lifestage, ef_data$species
    ))
  ef_data$visitID <- factor(paste(ef_data$date, ef_data$siteID))

  efdat <-
    efp(
      count ~ lifestage + s(visitID, bs = "re") + s(sampleID, bs = "re"),
      data = ef_data, pass = pass, id = sampleID,
      fit = FALSE
    )

  expect_equal(ncol(efdat$X), 2)
  expect_equal(ncol(efdat$Z), nlevels(ef_data$visitID) + nlevels(ef_data$sampleID))
  expect_equal(
    efdat$re_id,
    rep(1:2, c(nlevels(ef_data$visitID), nlevels(ef_data$sampleID)))
  )

  set.seed(1)
  params <-
    list(
      alpha = rnorm(ncol(efdat$Z)),
      beta = c(0.5, -0.2),
      log_sigma = c(-1, 0.5)
    )

  # joint likelihood of data and random effects
  nll <- ef:::ef_nll(params, efdat)

  efdat_fixed <- efdat
  efdat_fixed$offset <- as.vector(efdat$Z %*% params$alpha)
  r_nll <-
    nloglik(efdat_fixed, params$beta) -
    sum(dnorm(params$alpha, 0, exp(params$log_sigma)[efdat$re_id], log = TRUE))

  expect_equal(nll, r_nll)
})


test_that("visit and sample random effects fit", {
  data(ef_data)
  ef_data$sampleID <-
    factor(paste(
      ef_data$date, ef_data$siteID,
      ef_data$lifestage, ef_data$species
    ))
  ef_data$visitID <- factor(paste(ef_data$date, ef_data$siteID))

  m0 <- efp(count ~ lifestage,
    data = ef_data, pass = pass, id = sampleID
  )

  m1 <- efp(count ~ lifestage + s(visitID, bs = "re"),
    data = ef_data, pass = pass, id = sampleID
  )

  m2 <- efp(count ~ lifestage + s(visitID, bs = "re") + s(sampleID, bs = "re"),
    data = ef_data, pass = pass, id = sampleID
  )

  expect_equal(m1$convergence, 0)
  expect_equal(m2$convergence, 0)

  # fixed effects only
  expect_equal(names(coef(m1)), names(coef(m0)))
  expect_equal(names(coef(m2)), names(coef(m0)))

  # one sigma per random effect term
  expect_equal(names(m1$sigma), "s(visitID)")
  expect_equal(names(m2$sigma), c("s(visitID)", "s(sampleID)"))
  expect_true(all(m2$sigma > 0))
  expect_equal(length(m1$random), nlevels(ef_data$visitID))
  expect_equal(
    length(m2$random),
    nlevels(ef_data$visitID) + nlevels(ef_data$sampleID)
  )

  # models are nested, so the marginal likelihood can not decrease
  expect_gte(as.numeric(logLik(m1)), as.numeric(logLik(m0)) - 1e-6)
  expect_gte(as.numeric(logLik(m2)), as.numeric(logLik(m1)) - 1e-6)
  expect_equal(as.numeric(logLik(m2)), m2$llik)
  expect_equal(AIC(m2), -2 * m2$llik + 2 * 4)

  # fitted values are conditional on the random effects and are in the
  # same order as the data
  eta <-
    m2$G[, m2$fx_cols, drop = FALSE] %*% coef(m2) +
    m2$G[, m2$re_cols, drop = FALSE] %*% m2$random
  expect_equal(fitted(m2), as.vector(1 / (1 + exp(-eta))))

  # variance of fixed and random effects
  expect_equal(dim(m2$Vb), c(2, 2))
  expect_equal(dim(m2$Vp), rep(ncol(m2$Gfit), 2))
  expect_true(all(diag(m2$Vp) > 0))

  expect_output(print(m2), "sigma s\\(visitID\\)")
})

test_that("gam conversion", {
  data(ef_data)
  ef_data$sampleID <-
    factor(paste(
      ef_data$date, ef_data$siteID,
      ef_data$lifestage, ef_data$species
    ))
  ef_data$visitID <- factor(paste(ef_data$date, ef_data$siteID))

  m1 <- efp(count ~ lifestage,
    data = ef_data, pass = pass, id = sampleID
  )
  m2 <- efp(count ~ lifestage + s(visitID, bs = "re"),
    data = ef_data, pass = pass, id = sampleID
  )

  m1_gam <- suppressMessages(as.gam(m1))
  m2_gam <- suppressMessages(as.gam(m2))

  expect_equal(
    as.vector(predict(m1_gam, newdata = ef_data, type = "response")),
    fitted(m1)
  )
  expect_equal(
    as.vector(predict(m2_gam, newdata = ef_data, type = "response")),
    fitted(m2)
  )
})
