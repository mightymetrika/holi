test_that("run_sim_rstar_glm works with logistic model", {

  # Test with 3 skewed covariatesovariates
  # set.seed(140)
  # sig <- rWishart(1, 20, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  sig <- matrix(c(
    6.548557, 5.179998, 4.568447,
    5.179998, 6.615009, 6.790185,
    4.568447, 6.790185, 9.491619
  ), nrow = 3, ncol = 3, byrow = TRUE)
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
  # set.seed(249)
  # sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  sig <- matrix(c(
    1.0700681, 0.9252148, 0.9915194,
    0.9252148, 1.7532156, 1.9882893,
    0.9915194, 1.9882893, 3.5348444
  ), nrow = 3, ncol = 3, byrow = TRUE)
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
  # set.seed(482)
  # sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  sig <- matrix(c(
    1.574233, 1.647567, 0.423899,
    1.647567, 3.666186, 2.042598,
    0.423899, 2.042598, 2.232371
  ), nrow = 3, ncol = 3, byrow = TRUE)
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
  # set.seed(1124)
  # sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  sig <- matrix(c(
    2.063991, 1.509164, 1.454282,
    1.509164, 3.057952, 3.344217,
    1.454282, 3.344217, 4.698425
  ), nrow = 3, ncol = 3, byrow = TRUE)
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
  # set.seed(3874)
  # sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  sig <- matrix(c(
    0.7724809, 0.7610845, 0.2756609,
    0.7610845, 2.2112480, 1.2285775,
    0.2756609, 1.2285775, 1.3178022
  ), nrow = 3, ncol = 3, byrow = TRUE)
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
  # set.seed(34729)
  # sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  sig <- matrix(c(
    2.1768436,  0.6693993, -0.2795894,
    0.6693993,  0.9687892,  0.3094327,
    -0.2795894,  0.3094327,  0.8243939
  ), nrow = 3, ncol = 3, byrow = TRUE)
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
