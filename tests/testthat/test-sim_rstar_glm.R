test_that("multiplication works", {


  # simulation_result <- sim_rstar_glm(n_main = 100, n_control = 100, n_covariates = 5,
  #                                    true_coef_main = c(1, -0.5, 0.3, 0.7, -0.2),
  #                                    true_coef_control = c(1, -0.5, 0.3, 0.7, -0.2),
  #                                    treatment_effect = 0.5,
  #                                    model = "logistic", Sigma_main = toeplitz((5:1)/(5*5)),
  #                                    Sigma_control = toeplitz((5:1)/(5*5)))
  #
  # simulation_result <- sim_rstar_glm(n_main = 100, n_control = 100, n_covariates = 5,
  #                                    true_coef_main = c(1, -0.5, 0.3, 0.7, -0.2),
  #                                    true_coef_control = c(1, -0.5, 0.3, 0.7, -0.2),
  #                                    treatment_effect = 0,
  #                                    model = "logistic",
  #                                    Sigma_main = rWishart(1, 10, Sigma = toeplitz((5:1)/(5*5)))[,,1],
  #                                    Sigma_control = rWishart(1, 10, Sigma = toeplitz((5:1)/(5*5)))[,,1])

  # Test with 3 covariates
  sig <- rWishart(1, 10, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  run_simulation_result1 <- run_sim_rstar_glm(n_sims = 10,
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
  run_simulation_result2 <- run_sim_rstar_glm(n_sims = 10,
                                             n_main = 15, n_control = 15, n_covariates = 0,
                                             true_coef_main = c(0),
                                             true_coef_control = c(0),
                                             treatment_effect = 0.5,
                                             model = "logistic") |>
    suppressWarnings()

  expect_equal(length(run_simulation_result2), 2)

  # sig <- rWishart(1, 5, Sigma = toeplitz((1)/(1*5)))[,,1]
  # run_simulation_result <- run_sim_rstar_glm(n_sims = 10,
  #                                            n_main = 20, n_control = 20, n_covariates = 1,
  #                                            true_coef_main = c(0),
  #                                            true_coef_control = c(0),
  #                                            treatment_effect = 0.5,
  #                                            model = "logistic",
  #                                            Sigma_main = diag(1),
  #                                            Sigma_control = diag(1)) |>
  #   suppressWarnings()
  #
  #
  #
  #
  #
  #
  # sig <- rWishart(1, 4, Sigma = toeplitz((3:1)/3))[,,1]
  # run_simulation_result <- run_sim_rstar_glm(n_sims = 10,
  #                                            n_main = 7, n_control = 7, n_covariates = 3,
  #                                            true_coef_main = c(.5, -0.3, 0.2),
  #                                            true_coef_control = c(.5, -0.3, 0.2),
  #                                            treatment_effect = 0,
  #                                            model = "logistic",
  #                                            Sigma_main = sig,
  #                                            Sigma_control = sig)
  #
  #
  #
  #
  # sig <- rWishart(1, 7, Sigma = toeplitz((3:1)/(3*5)))[,,1]
  # run_simulation_result <- run_sim_rstar_glm(n_sims = 10,
  #                                            n_main = 20, n_control = 20, n_covariates = 0,
  #                                            true_coef_main = NULL,
  #                                            true_coef_control = NULL,
  #                                            treatment_effect = 0.5,
  #                                            model = "logistic",
  #                                            Sigma_main = NULL,
  #                                            Sigma_control = NULL)

})
