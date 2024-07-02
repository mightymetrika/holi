test_that("rstar_glm works for logistic", {

  # Test with Seatbelts data (no confidence intervals)
  rs_seatbelts <- rstar_glm(law ~ DriversKilled + VanKilled + drivers + kms,
                            .data = Seatbelts,
                            .model = "logistic") |> suppressWarnings()

  expect_s3_class(rs_seatbelts, "rstar_glm_result")
  expect_null(rs_seatbelts$rs_ci)

  # Test with Seatbelts data (with confidence intervals)
  rs_seatbelts <- rstar_glm(law ~ DriversKilled + VanKilled + drivers + kms,
                            .data = Seatbelts,
                            .model = "logistic",
                            .rstar.ci = TRUE) |> suppressWarnings()

  expect_s3_class(rs_seatbelts, "rstar_glm_result")
  expect_s3_class(rs_seatbelts$rs_ci, "rstarci")
})


test_that("rstar_glm works for linear", {

  # Linear regression example (no confidence intervals)
  rs_linear <- rstar_glm(mpg ~ wt + hp,
                         .data = mtcars,
                         .model = "linear") |> suppressWarnings()

  expect_s3_class(rs_linear, "rstar_glm_result")
  expect_null(rs_linear$rs_ci)

  # Linear regression example (with confidence intervals)
  rs_linear <- rstar_glm(mpg ~ wt + hp,
                         .data = mtcars,
                         .model = "linear",
                         .rstar.ci = TRUE) |> suppressWarnings()

  expect_s3_class(rs_linear, "rstar_glm_result")
  expect_s3_class(rs_linear$rs_ci, "rstarci")
})

test_that("rstar_glm works for poisson", {

  # Poisson regression example (without confidence intervals)
  rs_poisson <- rstar_glm(count ~ spray,
                          .data = InsectSprays,
                          .model = "poisson") |> suppressWarnings()

  expect_s3_class(rs_poisson, "rstar_glm_result")
  expect_null(rs_poisson$rs_ci)

  # Poisson regression example (with confidence intervals)
  rs_poisson <- rstar_glm(count ~ spray,
                          .data = InsectSprays,
                          .model = "poisson",
                          .rstar.ci = TRUE) |> suppressWarnings()

  expect_s3_class(rs_poisson, "rstar_glm_result")
  expect_s3_class(rs_poisson$rs_ci, "rstarci")

})
