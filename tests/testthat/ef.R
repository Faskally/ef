context("fitted values")
library(ef)
library(mgcv)

# do a model fit, ef_data

data(ef_data)

ef_data $ pass12 <- factor(replace(ef_data $ pass, ef_data $ pass > 2, 2))
ef_data $ sampleID <- as.numeric(factor(paste(ef_data $ date, ef_data $ siteID,
                                              ef_data $ lifestage, ef_data $ species)))
ef_data $ visitID <- as.numeric(factor(paste(ef_data $ date, ef_data $ siteID)))
m1 <- efp(count ~ pass12 + lifestage + siteID + year,
          data = ef_data, pass = pass, id = sampleID)


test_that("model_output", {
#  expect_equal(str_length("a"), 1)
#  expect_equal(str_length("ab"), 2)
#  expect_equal(str_length("abc"), 3)
})


# Get Predicted Values for _P_
m1_gam <- as.gam(m1)
ef_data $ prob <- predict.gam(m1_gam, newdata=ef_data, type = "response")

test_that("as.gam", {
  #  expect_equal(str_length("a"), 1)
  #  expect_equal(str_length("ab"), 2)
  #  expect_equal(str_length("abc"), 3)
})


# Comparing Models

if (FALSE) {
  large_model <- efp(count ~ pass12 + lifestage + siteID + year,
    data = ef_data, pass = pass, id = sampleID
  )

  od_estimate <- overdispersion(
    data = ef_data, visitID = "visitID",
    siteID = "siteID", sampleID = "sampleID",
    largemodel = large_model
  )


  test_that("od_estimate", {
    #  expect_equal(str_length("a"), 1)
    #  expect_equal(str_length("ab"), 2)
    #  expect_equal(str_length("abc"), 3)
  })

  mfull <- efp(count ~ pass12 + lifestage + siteID + year,
    data = ef_data, pass = pass, id = sampleID
  )

  BICadj(mfull, ef_data, od_estimate)

  test_that("BICadj", {
    #  expect_equal(str_length("a"), 1)
    #  expect_equal(str_length("ab"), 2)
    #  expect_equal(str_length("abc"), 3)
  })


  mlife <- efp(count ~ pass12 + siteID + year,
    data = ef_data, pass = pass, id = sampleID
  )

  BICadj(mlife, ef_data, od_estimate)

  test_that("BICadj", {
    #  expect_equal(str_length("a"), 1)
    #  expect_equal(str_length("ab"), 2)
    #  expect_equal(str_length("abc"), 3)
  })
}
