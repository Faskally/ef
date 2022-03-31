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
      log_sigma = 0
    )

  obj <-
    MakeADFun(
      efdat,
      params,
      DLL = "ef",
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

  expect_equivalent(coef(m1), 0.9734365, tolerance = 1e-5)
  expect_equivalent(as.numeric(logLik(m1)), -1612.81836111119, tolerance = 1e-9)

  expect_equivalent(coef(m2), c(0.6667471, 0.6449035, -0.2562596), tolerance = 1e-4)
  expect_equivalent(as.numeric(logLik(m2)), -1592.30666991703, tolerance = 1e-9)

  expect_equivalent(
    coef(m3),
    c(0.67588338, 1.07064066, -0.80165087, -0.39616011, 0.19015039, 0.67198771, 0.02512267),
    tolerance = 1e-4
  )
  expect_equivalent(as.numeric(logLik(m3)), -1584.7000365589, tolerance = 1e-9)
})

test_that("gam conversion", {
  # m1_gam <- as.gam(m1)
  # m2_gam <- as.gam(m2)
  # m3_gam <- as.gam(m2)
})

test_that("overdispersion estimate", {

  data(ef_data)

  ef_data$pass12 <- factor(replace(ef_data$pass, ef_data$pass > 2, 2))

  ef_data$sampleID <-
    as.numeric(factor(paste(
      ef_data$date, ef_data$siteID,
      ef_data$lifestage, ef_data$species
    )))

  ef_data$visitID <- as.numeric(factor(paste(ef_data$date, ef_data$siteID)))

  large_model <- efp(count ~ pass12 + lifestage + siteID + year,
    data = ef_data, pass = pass, id = sampleID
  )

  od_estimate <- overdispersion(
    data = ef_data, visitID = "visitID",
    siteID = "siteID", sampleID = "sampleID",
    largemodel = large_model
  )

  mfull <- efp(count ~ pass12 + lifestage + siteID + year,
    data = ef_data, pass = pass, id = sampleID
  )

  bica <- BICadj(mfull, ef_data, od_estimate)
  expect_equivalent(bica, -3129.33815406195, tol = 1e-9)
})