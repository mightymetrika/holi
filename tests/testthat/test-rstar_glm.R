test_that("rstar_glm works for logistic", {

  # Test with Seatbelts data
  rs_seatbelts <- rstar_glm(law ~ DriversKilled + VanKilled + drivers + kms,
                            .data = Seatbelts,
                            .model = "logistic") |> suppressWarnings()

  expect_s3_class(rs_seatbelts, "rstar_glm_result")
})


test_that("rstar_glm works for linear", {

  # Linear regression example
  rs_linear <- rstar_glm(mpg ~ wt + hp,
                         .data = mtcars,
                         .model = "linear") |> suppressWarnings()

  expect_s3_class(rs_linear, "rstar_glm_result")
})

test_that("rstar_glm works for poisson", {

  # Poisson regression example
  rs_poisson <- rstar_glm(count ~ spray,
                          .data = InsectSprays,
                          .model = "poisson") |> suppressWarnings()

  expect_s3_class(rs_poisson, "rstar_glm_result")
})
