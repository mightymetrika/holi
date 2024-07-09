test_that("run_sim_rstar_glm works with logistic model", {

  # Test with 3 skewed covariatesovariates
  sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  run_simulation_result1 <- run_sim_rstar_glm(n_sims = 2,
                                             n_main = 20, n_control = 20, n_covariates = 3,
                                     true_coef_main = c(.5, -0.3, 0.2),
                                     true_coef_control = c(.5, -0.3, 0.2),
                                     treatment_effect = 0.5,
                                     model = "logistic",
                                     Sigma_main = sig,
                                     Sigma_control = sig,
                                     skewness_main = c(0.08, 0.08, 0.08), skewness_control = c(0.02, 0.02, 0.02)) |>
    suppressWarnings()

  expect_equal(length(run_simulation_result1), 2)

  # Test with 0 covariates
  run_simulation_result2 <- run_sim_rstar_glm(n_sims = 2,
                                             n_main = 15, n_control = 15, n_covariates = 0,
                                             true_coef_main = c(0),
                                             true_coef_control = c(0),
                                             treatment_effect = 0.5,
                                             model = "logistic") |>
    suppressWarnings()

  expect_equal(length(run_simulation_result2), 2)

  # Test with 3 non-skewed covariatesovariates
  sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  run_simulation_result3 <- run_sim_rstar_glm(n_sims = 2,
                                              n_main = 20, n_control = 20, n_covariates = 3,
                                              true_coef_main = c(.5, -0.3, 0.2),
                                              true_coef_control = c(.5, -0.3, 0.2),
                                              treatment_effect = 0.5,
                                              model = "logistic",
                                              Sigma_main = sig,
                                              Sigma_control = sig) |>
    suppressWarnings()

  expect_equal(length(run_simulation_result3), 2)

})


test_that("run_sim_rstar_glm works with linear model", {

  # Test with 3 skewed covariatesovariates
  sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  run_simulation_result1 <- run_sim_rstar_glm(n_sims = 2,
                                              n_main = 20, n_control = 20, n_covariates = 3,
                                              true_coef_main = c(.5, -0.3, 0.2),
                                              true_coef_control = c(.5, -0.3, 0.2),
                                              treatment_effect = 0.5,
                                              model = "linear",
                                              Sigma_main = sig,
                                              Sigma_control = sig,
                                              skewness_main = c(0.08, 0.08, 0.08), skewness_control = c(0.02, 0.02, 0.02)) |>
    suppressWarnings()

  expect_equal(length(run_simulation_result1), 2)

  # Test with 0 covariates
  run_simulation_result2 <- run_sim_rstar_glm(n_sims = 2,
                                              n_main = 15, n_control = 15, n_covariates = 0,
                                              true_coef_main = c(0),
                                              true_coef_control = c(0),
                                              treatment_effect = 0.5,
                                              model = "linear") |>
    suppressWarnings()

  expect_equal(length(run_simulation_result2), 2)

  # Test with 3 non-skewed covariatesovariates
  sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  run_simulation_result3 <- run_sim_rstar_glm(n_sims = 2,
                                              n_main = 20, n_control = 20, n_covariates = 3,
                                              true_coef_main = c(.5, -0.3, 0.2),
                                              true_coef_control = c(.5, -0.3, 0.2),
                                              treatment_effect = 0.5,
                                              model = "linear",
                                              Sigma_main = sig,
                                              Sigma_control = sig) |>
    suppressWarnings()

  expect_equal(length(run_simulation_result3), 2)

})


test_that("run_sim_rstar_glm works with poisson model", {

  # Test with 3 skewed covariatesovariates
  sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  run_simulation_result1 <- run_sim_rstar_glm(n_sims = 2,
                                              n_main = 20, n_control = 20, n_covariates = 3,
                                              true_coef_main = c(.5, -0.3, 0.2),
                                              true_coef_control = c(.5, -0.3, 0.2),
                                              treatment_effect = 0.5,
                                              model = "poisson",
                                              Sigma_main = sig,
                                              Sigma_control = sig,
                                              skewness_main = c(0.08, 0.08, 0.08), skewness_control = c(0.02, 0.02, 0.02)) |>
    suppressWarnings()

  expect_equal(length(run_simulation_result1), 2)

  # Test with 0 covariates
  run_simulation_result2 <- run_sim_rstar_glm(n_sims = 2,
                                              n_main = 15, n_control = 15, n_covariates = 0,
                                              true_coef_main = c(0),
                                              true_coef_control = c(0),
                                              treatment_effect = 0.5,
                                              model = "poisson") |>
    suppressWarnings()

  expect_equal(length(run_simulation_result2), 2)

  # Test with 3 non-skewed covariatesovariates
  sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  run_simulation_result3 <- run_sim_rstar_glm(n_sims = 2,
                                              n_main = 20, n_control = 20, n_covariates = 3,
                                              true_coef_main = c(.5, -0.3, 0.2),
                                              true_coef_control = c(.5, -0.3, 0.2),
                                              treatment_effect = 0.5,
                                              model = "poisson",
                                              Sigma_main = sig,
                                              Sigma_control = sig) |>
    suppressWarnings()

  expect_equal(length(run_simulation_result3), 2)

})
