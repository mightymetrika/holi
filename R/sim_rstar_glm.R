sim_rstar_glm <- function(n_main, n_covariates, true_coef_main,
                          n_control = NULL, true_coef_control = NULL,
                          treatment_effect = NULL,
                          model = c("logistic", "linear", "poisson"),
                          skewness_main = NULL, skewness_control = NULL,
                          Sigma_main = NULL, Sigma_control = NULL, ...) {
  model <- match.arg(model)

  # Helper function to generate covariate data
  generate_covariate_data <- function(n, n_covariates, Sigma, skewness = NULL) {
    if (n_covariates == 0) return(matrix(1, nrow = n, ncol = 1))  # Intercept only

    if (is.null(Sigma)) Sigma <- diag(n_covariates)

    if (is.null(skewness)) {
      X <- tryCatch({
        MASS::mvrnorm(n, mu = rep(0, n_covariates), Sigma = Sigma)
      }, error = function(e) {
        warning("Covariate data generation failed: ", e$message)
        return(NULL)
      })
    } else {
      cpM <- list(mean = rep(0, n_covariates), var.cov = Sigma, gamma1 = skewness)
      dpM <- sn::cp2dp(cpM, family = "SN")
      X <- tryCatch({
        sn::rmsn(n, dp = dpM)
      }, error = function(e) {
        warning("Covariate data generation with skewness failed: ", e$message)
        return(NULL)
      })
    }
    #return(cbind(1, X))  # Add intercept
    return(X)
  }

  # Helper function to generate outcome variable
  generate_outcome <- function(X, true_coef, model, treatment_effect = NULL, group = NULL) {
    if (is.null(X)) return(NULL)
    eta <- X %*% true_coef
    if (!is.null(treatment_effect) & !is.null(group)) {
      eta <- eta + group * treatment_effect
    }
    y <- tryCatch({
      if (model == "logistic") {
        p <- stats::plogis(eta)
        stats::rbinom(nrow(X), size = 1, prob = p)
      } else if (model == "linear") {
        eta + stats::rnorm(nrow(X), mean = 0, sd = 1)
      } else if (model == "poisson") {
        mu <- exp(eta)
        stats::rpois(nrow(X), lambda = mu)
      }
    }, error = function(e) {
      warning("Outcome data generation failed: ", e$message)
      return(NULL)
    })
    return(y)
  }

  # Generate data for main group
  X_main <- generate_covariate_data(n_main, n_covariates, Sigma_main, skewness_main)
  y_main <- generate_outcome(X_main, true_coef_main, model, treatment_effect, group = 1)
  if (is.null(y_main)) return(NULL)
  data_main <- data.frame(y = y_main, X_main)
  if(n_covariates == 0){
    colnames(data_main) <- c("y", paste0("X", 0))
  } else{
    colnames(data_main) <- c("y", paste0("X", 1:n_covariates))
  }


  # Check if control group data is needed
  if (!is.null(n_control) && !is.null(true_coef_control)) {
    X_control <- generate_covariate_data(n_control, n_covariates, Sigma_control, skewness_control)
    y_control <- generate_outcome(X_control, true_coef_control, model, treatment_effect, group = 0)
    if (is.null(y_control)) return(NULL)
    data_control <- data.frame(y = y_control, X_control)
    if(n_covariates == 0){
      colnames(data_control) <- c("y", paste0("X", 0))
    } else{
      colnames(data_control) <- c("y", paste0("X", 1:n_covariates))
    }

    # Combine main and control data, add group indicator
    data_control$group <- 0
    data_main$group <- 1
    data <- rbind(data_main, data_control)

    # Update formula to include group indicator
    if (n_covariates > 0) {
      formula <- stats::as.formula(paste("y ~ group +", paste(colnames(data)[-(1:2)], collapse = " + ")))
    } else {
      formula <- stats::as.formula("y ~ group")
    }
  } else {
    data <- data_main
    if (n_covariates > 0) {
      formula <- stats::as.formula(paste("y ~", paste(colnames(data)[-1], collapse = " + ")))
    } else {
      formula <- stats::as.formula("y ~ 1")
    }
  }

  # Fit the models
  fit_glm <- tryCatch({
    stats::glm(formula = formula, family = if (model == "logistic") stats::binomial() else if (model == "linear") stats::gaussian() else stats::poisson(), data = data)
  }, error = function(e) {
    warning("GLM model fitting failed: ", e$message)
    return(NULL)
  })

  fit_rstar <- tryCatch({
    rstar_glm(.formula = formula, .data = data, .model = model, ...)
  }, error = function(e) {
    warning("rstar model fitting failed: ", e$message)
    return(NULL)
  })

  if (is.null(fit_glm) || is.null(fit_rstar)) return(NULL)

  return(list(rstar = fit_rstar, fit_glm = fit_glm, data = data))
}

run_sim_rstar_glm <- function(n_sims, alpha_level = 0.05,
                              n_main, n_covariates, true_coef_main,
                              n_control = NULL, true_coef_control = NULL,
                              treatment_effect = NULL,
                              model = c("logistic", "linear", "poisson"),
                              skewness_main = NULL, skewness_control = NULL,
                              Sigma_main = NULL, Sigma_control = NULL, ...) {

  results <- replicate(n_sims, {
    sim_result <- sim_rstar_glm(n_main, n_covariates, true_coef_main,
                                n_control, true_coef_control,
                                treatment_effect, model,
                                skewness_main, skewness_control,
                                Sigma_main, Sigma_control, ...)
    if (is.null(sim_result)) return(rep(NA, 6))

    # Adjust for the case when there's no group variable (i.e., no control group)
    if ("group" %in% names(sim_result$fit_glm$coefficients)) {
      glm_coef <- sim_result$fit_glm$coefficients["group"]
      glm_se <- summary(sim_result$fit_glm)$coefficients["group", "Std. Error"]
      glm_p_value <- summary(sim_result$fit_glm)$coefficients["group", "Pr(>|z|)"]
      rs_estimate <- sim_result$rstar$rs$theta.hat["group"]
    } else {
      glm_coef <- sim_result$fit_glm$coefficients["(Intercept)"]
      glm_se <- summary(sim_result$fit_glm)$coefficients["(Intercept)", "Std. Error"]
      glm_p_value <- summary(sim_result$fit_glm)$coefficients["(Intercept)", "Pr(>|z|)"]
      rs_estimate <- sim_result$rstar$rs$theta.hat["(Intercept)"]
    }

    c(
      glm_estimate = glm_coef,
      glm_p_value = glm_p_value,
      rs_estimate = rs_estimate,
      rs_p_value = stats::pnorm(sim_result$rstar$rs$rs),
      glm_se = glm_se
    )
  }, simplify = "matrix")

  results_df <- as.data.frame(t(results))
  converged <- stats::complete.cases(results_df)
  results_df <- results_df[converged, ]

  # Calculate summary statistics
  n_converged <- sum(converged)
  rejection_rate_glm <- mean(results_df$glm_p_value < alpha_level, na.rm = TRUE)
  rejection_rate_rs <- mean(results_df$rs_p_value < alpha_level, na.rm = TRUE)

  # Standard error of rejection rate (using binomial proportion standard error)
  se_rejection_rate_glm <- sqrt(rejection_rate_glm * (1 - rejection_rate_glm) / n_converged)
  se_rejection_rate_rs <- sqrt(rejection_rate_rs * (1 - rejection_rate_rs) / n_converged)

  # Bias
  bias_glm <- mean(results_df$glm_estimate - treatment_effect, na.rm = TRUE)
  bias_rs <- mean(results_df$rs_estimate - treatment_effect, na.rm = TRUE)

  # Empirical SE
  empse_glm <- sqrt(stats::var(results_df$glm_estimate, na.rm = TRUE))
  empse_rs <- sqrt(stats::var(results_df$rs_estimate, na.rm = TRUE))

  # MSE
  mse_glm <- mean((results_df$glm_estimate - treatment_effect)^2, na.rm = TRUE)
  mse_rs <- mean((results_df$rs_estimate - treatment_effect)^2, na.rm = TRUE)

  # RMSE
  rmse_glm <- sqrt(mse_glm)
  rmse_rs <- sqrt(mse_rs)

  # Format summary statistics as a data frame
  summary_df <- data.frame(
    Model = c("GLM", "RS"),
    Rejection_Rate = c(rejection_rate_glm, rejection_rate_rs),
    Rejection_Rate_SE = c(se_rejection_rate_glm, se_rejection_rate_rs),
    Bias = c(bias_glm, bias_rs),
    EmpSE = c(empse_glm, empse_rs),
    MSE = c(mse_glm, mse_rs),
    RMSE = c(rmse_glm, rmse_rs),
    Converged_Proportion = rep(n_converged / n_sims, 2),
    Iter = rep(n_sims, 2),
    Alpha = rep(alpha_level, 2),
    n_covs = rep(n_covariates, 2),
    n_main = rep(n_main, 2),
    n_cont = rep(n_control, 2),
    true_coef_main = rep(paste(true_coef_main, collapse = ","), 2),
    true_coef_control = rep(paste(true_coef_control, collapse = ","), 2),
    treatment_effect = rep(treatment_effect, 2),
    mod = rep(model, 2),
    skewness_main = rep(paste(skewness_main, collapse = ","), 2),
    skewness_control = rep(paste(skewness_control, collapse = ","), 2),
    sigma_main = rep(paste(Sigma_main, collapse = ","), 2),
    sigma_control = rep(paste(Sigma_control, collapse = ","), 2)
  )

  # Round numeric values to 4 decimal places
  numeric_columns <- sapply(summary_df, is.numeric)
  summary_df[, numeric_columns] <- round(summary_df[, numeric_columns], 4)

  return(list(results = results_df, summary = summary_df))
}










# sim_rstar_glm <- function(n_main, n_covariates, true_coef_main,
#                           n_control = NULL, true_coef_control = NULL,
#                           treatment_effect = NULL,
#                           model = c("logistic", "linear", "poisson"),
#                           skewness_main = NULL, skewness_control = NULL,
#                           Sigma_main, Sigma_control, ...) {
#   model <- match.arg(model)
#
#   # Helper function to generate covariate data
#   generate_covariate_data <- function(n, n_covariates, Sigma, skewness = NULL) {
#     if (is.null(skewness)) {
#       X <- tryCatch({
#         MASS::mvrnorm(n, mu = rep(0, n_covariates), Sigma = Sigma)
#       }, error = function(e) {
#         warning("Covariate data generation failed: ", e$message)
#         return(NULL)
#       })
#     } else {
#       cpM <- list(mean = rep(0, n_covariates), var.cov = Sigma, gamma1 = skewness)
#       dpM <- sn::cp2dp(cpM, family = "SN")
#       X <- tryCatch({
#         sn::rmsn(n, dp = dpM)
#       }, error = function(e) {
#         warning("Covariate data generation with skewness failed: ", e$message)
#         return(NULL)
#       })
#     }
#     return(X)
#   }
#
#   # Helper function to generate outcome variable
#   generate_outcome <- function(X, true_coef, model, treatment_effect = NULL, group = NULL) {
#     if (is.null(X)) return(NULL)
#     eta <- X %*% true_coef
#     if (!is.null(treatment_effect) & !is.null(group)) {
#       eta <- eta + group * treatment_effect
#     }
#     y <- tryCatch({
#       if (model == "logistic") {
#         p <- stats::plogis(eta)
#         stats::rbinom(nrow(X), size = 1, prob = p)
#       } else if (model == "linear") {
#         eta + stats::rnorm(nrow(X), mean = 0, sd = 1)
#       } else if (model == "poisson") {
#         mu <- exp(eta)
#         stats::rpois(nrow(X), lambda = mu)
#       }
#     }, error = function(e) {
#       warning("Outcome data generation failed: ", e$message)
#       return(NULL)
#     })
#     return(y)
#   }
#
#   # Generate data for main group
#   X_main <- generate_covariate_data(n_main, n_covariates, Sigma_main, skewness_main)
#   y_main <- generate_outcome(X_main, true_coef_main, model, treatment_effect, group = 1)
#   if (is.null(y_main)) return(NULL)
#   data_main <- data.frame(y = y_main, X_main)
#   colnames(data_main) <- c("y", paste0("X", 1:n_covariates))
#
#   # Check if control group data is needed
#   if (!is.null(n_control) && !is.null(true_coef_control)) {
#     X_control <- generate_covariate_data(n_control, n_covariates, Sigma_control, skewness_control)
#     y_control <- generate_outcome(X_control, true_coef_control, model, treatment_effect, group = 0)
#     if (is.null(y_control)) return(NULL)
#     data_control <- data.frame(y = y_control, X_control)
#     colnames(data_control) <- c("y", paste0("X", 1:n_covariates))
#
#     # Combine main and control data, add group indicator
#     data_control$group <- 0
#     data_main$group <- 1
#     data <- rbind(data_main, data_control)
#
#     # Update formula to include group indicator
#     formula <- as.formula(paste("y ~ group +", paste(colnames(data)[-c(1, ncol(data))], collapse = " + ")))
#   } else {
#     data <- data_main
#     formula <- as.formula(paste("y ~", paste(colnames(data)[-1], collapse = " + ")))
#   }
#
#   # Fit the models
#   fit_glm <- tryCatch({
#     stats::glm(formula = formula, family = if (model == "logistic") stats::binomial() else if (model == "linear") stats::gaussian() else stats::poisson(), data = data)
#   }, error = function(e) {
#     warning("GLM model fitting failed: ", e$message)
#     return(NULL)
#   })
#
#   fit_rstar <- tryCatch({
#     rstar_glm(.formula = formula, .data = data, .model = model, ...)
#   }, error = function(e) {
#     warning("rstar model fitting failed: ", e$message)
#     return(NULL)
#   })
#
#   if (is.null(fit_glm) || is.null(fit_rstar)) return(NULL)
#
#   return(list(rstar = fit_rstar, fit_glm = fit_glm, data = data))
# }
#
# run_sim_rstar_glm <- function(n_sims, alpha_level = 0.05,
#                               n_main, n_covariates, true_coef_main,
#                               n_control = NULL, true_coef_control = NULL,
#                               treatment_effect = NULL,
#                               model = c("logistic", "linear", "poisson"),
#                               skewness_main = NULL, skewness_control = NULL,
#                               Sigma_main, Sigma_control, ...) {
#
#   results <- replicate(n_sims, {
#     sim_result <- sim_rstar_glm(n_main, n_covariates, true_coef_main,
#                                 n_control, true_coef_control,
#                                 treatment_effect, model,
#                                 skewness_main, skewness_control,
#                                 Sigma_main, Sigma_control, ...)
#     if (is.null(sim_result)) return(rep(NA, 6))
#     glm_se <- summary(sim_result$fit_glm)$coefficients["group", "Std. Error"]
#     c(
#       glm_estimate = sim_result$fit_glm$coefficients["group"],
#       glm_p_value = summary(sim_result$fit_glm)$coefficients["group", "Pr(>|z|)"],
#       rs_estimate = sim_result$rstar$rs$theta.hat["group"],
#       #r_p_value = stats::pnorm(sim_result$rstar$rs$r),
#       rs_p_value = stats::pnorm(sim_result$rstar$rs$rs),
#       glm_se = glm_se
#     )
#   }, simplify = "matrix")
#
#   results_df <- as.data.frame(t(results))
#   converged <- complete.cases(results_df)
#   results_df <- results_df[converged, ]
#
#   # Calculate summary statistics
#   n_converged <- sum(converged)
#   rejection_rate_glm <- mean(results_df$glm_p_value < alpha_level, na.rm = TRUE)
#   rejection_rate_rs <- mean(results_df$rs_p_value < alpha_level, na.rm = TRUE)
#
#   # Standard error of rejection rate (using binomial proportion standard error)
#   se_rejection_rate_glm <- sqrt(rejection_rate_glm * (1 - rejection_rate_glm) / n_converged)
#   se_rejection_rate_rs <- sqrt(rejection_rate_rs * (1 - rejection_rate_rs) / n_converged)
#
#   # Bias
#   bias_glm <- mean(results_df$glm_estimate - treatment_effect, na.rm = TRUE)
#   bias_rs <- mean(results_df$rs_estimate - treatment_effect, na.rm = TRUE)
#
#   # Empirical SE
#   empse_glm <- sqrt(var(results_df$glm_estimate, na.rm = TRUE))
#   empse_rs <- sqrt(var(results_df$rs_estimate, na.rm = TRUE))
#
#   # MSE
#   mse_glm <- mean((results_df$glm_estimate - treatment_effect)^2, na.rm = TRUE)
#   mse_rs <- mean((results_df$rs_estimate - treatment_effect)^2, na.rm = TRUE)
#
#   # RMSE
#   rmse_glm <- sqrt(mse_glm)
#   rmse_rs <- sqrt(mse_rs)
#
#   # Format summary statistics as a data frame
#   summary_df <- data.frame(
#     Model = c("GLM", "RS"),
#     Rejection_Rate = c(rejection_rate_glm, rejection_rate_rs),
#     Rejection_Rate_SE = c(se_rejection_rate_glm, se_rejection_rate_rs),
#     Bias = c(bias_glm, bias_rs),
#     EmpSE = c(empse_glm, empse_rs),
#     MSE = c(mse_glm, mse_rs),
#     RMSE = c(rmse_glm, rmse_rs),
#     Converged_Proportion = rep(n_converged / n_sims, 2)
#   )
#
#   # Round numeric values to 4 decimal places
#   summary_df[, -1] <- round(summary_df[, -1], 4)
#
#   return(list(results = results_df, summary = summary_df))
# }



# sim_rstar_glm <- function(n_main, n_covariates, true_coef_main,
#                           n_control = NULL, true_coef_control = NULL,
#                           treatment_effect = NULL,
#                           model = c("logistic", "linear", "poisson"),
#                           skewness_main = NULL, skewness_control = NULL,
#                           Sigma_main = NULL, Sigma_control = NULL, ...) {
#   model <- match.arg(model)
#
#   # Default Sigma to identity matrix if not provided
#   if (is.null(Sigma_main)) {
#     Sigma_main <- diag(n_covariates)
#   }
#   if (is.null(Sigma_control)) {
#     Sigma_control <- diag(n_covariates)
#   }
#
#   # Helper function to generate covariate data
#   generate_covariate_data <- function(n, n_covariates, Sigma, skewness = NULL) {
#     if (n_covariates == 0) {
#       return(matrix(nrow = n, ncol = 0))
#     }
#     if (is.null(skewness)) {
#       X <- tryCatch({
#         MASS::mvrnorm(n, mu = rep(0, n_covariates), Sigma = Sigma)
#       }, error = function(e) {
#         warning("Covariate data generation failed: ", e$message)
#         return(NULL)
#       })
#     } else {
#       cpM <- list(mean = rep(0, n_covariates), var.cov = Sigma, gamma1 = skewness)
#       dpM <- sn::cp2dp(cpM, family = "SN")
#       X <- tryCatch({
#         sn::rmsn(n, dp = dpM)
#       }, error = function(e) {
#         warning("Covariate data generation with skewness failed: ", e$message)
#         return(NULL)
#       })
#     }
#     return(X)
#   }
#
#   # Helper function to generate outcome variable
#   generate_outcome <- function(X, true_coef, model, treatment_effect = NULL, group = NULL) {
#     if (is.null(X)) return(NULL)
#     if (ncol(X) == 0) {
#       eta <- rep(0, nrow(X))
#     } else {
#       eta <- X %*% true_coef
#     }
#     if (!is.null(treatment_effect) & !is.null(group)) {
#       eta <- eta + group * treatment_effect
#     }
#     y <- tryCatch({
#       if (model == "logistic") {
#         p <- stats::plogis(eta)
#         stats::rbinom(nrow(X), size = 1, prob = p)
#       } else if (model == "linear") {
#         eta + stats::rnorm(nrow(X), mean = 0, sd = 1)
#       } else if (model == "poisson") {
#         mu <- exp(eta)
#         stats::rpois(nrow(X), lambda = mu)
#       }
#     }, error = function(e) {
#       warning("Outcome data generation failed: ", e$message)
#       return(NULL)
#     })
#     return(y)
#   }
#
#   # Generate data for main group
#   X_main <- generate_covariate_data(n_main, n_covariates, Sigma_main, skewness_main)
#   y_main <- generate_outcome(X_main, true_coef_main, model, treatment_effect, group = 1)
#   if (is.null(y_main)) return(NULL)
#   data_main <- data.frame(y = y_main, X_main)
#   if (n_covariates > 0) {
#     colnames(data_main) <- c("y", paste0("X", 1:n_covariates))
#   } else {
#     colnames(data_main) <- "y"
#   }
#
#   # Check if control group data is needed
#   if (!is.null(n_control) && !is.null(true_coef_control)) {
#     X_control <- generate_covariate_data(n_control, n_covariates, Sigma_control, skewness_control)
#     y_control <- generate_outcome(X_control, true_coef_control, model, treatment_effect, group = 0)
#     if (is.null(y_control)) return(NULL)
#     data_control <- data.frame(y = y_control, X_control)
#     if (n_covariates > 0) {
#       colnames(data_control) <- c("y", paste0("X", 1:n_covariates))
#     } else {
#       colnames(data_control) <- "y"
#     }
#
#     # Combine main and control data, add group indicator
#     data_control$group <- 0
#     data_main$group <- 1
#     data <- rbind(data_main, data_control)
#
#     # Update formula to include group indicator
#     if (n_covariates > 0) {
#       formula <- as.formula(paste("y ~ group +", paste(colnames(data)[-c(1, ncol(data))], collapse = " + ")))
#     } else {
#       formula <- as.formula("y ~ group")
#     }
#   } else {
#     data <- data_main
#     if (n_covariates > 0) {
#       formula <- as.formula(paste("y ~", paste(colnames(data)[-1], collapse = " + ")))
#     } else {
#       formula <- as.formula("y ~ 1")
#     }
#   }
#
#   # Fit the models
#   fit_glm <- tryCatch({
#     stats::glm(formula = formula, family = if (model == "logistic") stats::binomial() else if (model == "linear") stats::gaussian() else stats::poisson(), data = data)
#   }, error = function(e) {
#     warning("GLM model fitting failed: ", e$message)
#     return(NULL)
#   })
#
#   fit_rstar <- tryCatch({
#     rstar_glm(.formula = formula, .data = data, .model = model, ...)
#   }, error = function(e) {
#     warning("rstar model fitting failed: ", e$message)
#     return(NULL)
#   })
#
#   if (is.null(fit_glm) || is.null(fit_rstar)) return(NULL)
#
#   return(list(rstar = fit_rstar, fit_glm = fit_glm, data = data))
# }
#
# run_sim_rstar_glm <- function(n_sims, alpha_level = 0.05,
#                               n_main, n_covariates, true_coef_main,
#                               n_control = NULL, true_coef_control = NULL,
#                               treatment_effect = NULL,
#                               model = c("logistic", "linear", "poisson"),
#                               skewness_main = NULL, skewness_control = NULL,
#                               Sigma_main = NULL, Sigma_control = NULL, ...) {
#
#   # Default Sigma to identity matrix if not provided
#   if (is.null(Sigma_main)) {
#     Sigma_main <- diag(n_covariates)
#   }
#   if (is.null(Sigma_control)) {
#     Sigma_control <- diag(n_covariates)
#   }
#
#   results <- replicate(n_sims, {
#     sim_result <- sim_rstar_glm(n_main, n_covariates, true_coef_main,
#                                 n_control, true_coef_control,
#                                 treatment_effect, model,
#                                 skewness_main, skewness_control,
#                                 Sigma_main, Sigma_control, ...)
#     if (is.null(sim_result)) return(rep(NA, 6))
#     glm_se <- summary(sim_result$fit_glm)$coefficients["group", "Std. Error"]
#     c(
#       glm_estimate = sim_result$fit_glm$coefficients["group"],
#       glm_p_value = summary(sim_result$fit_glm)$coefficients["group", "Pr(>|z|)"],
#       rs_estimate = sim_result$rstar$rs$theta.hat["group"],
#       r_p_value = stats::pnorm(sim_result$rstar$rs$r),
#       rs_p_value = stats::pnorm(sim_result$rstar$rs$rs),
#       glm_se = glm_se
#     )
#   }, simplify = "matrix")
#
#   results_df <- as.data.frame(t(results))
#   converged <- complete.cases(results_df)
#   results_df <- results_df[converged, ]
#
#   # Calculate summary statistics
#   n_converged <- sum(converged)
#   rejection_rate_glm <- mean(results_df$glm_p_value < alpha_level, na.rm = TRUE)
#   rejection_rate_r <- mean(results_df$r_p_value < alpha_level, na.rm = TRUE)
#   rejection_rate_rs <- mean(results_df$rs_p_value < alpha_level, na.rm = TRUE)
#
#   # Standard error of rejection rate (using binomial proportion standard error)
#   se_rejection_rate_glm <- sqrt(rejection_rate_glm * (1 - rejection_rate_glm) / n_converged)
#   se_rejection_rate_r <- sqrt(rejection_rate_r * (1 - rejection_rate_r) / n_converged)
#   se_rejection_rate_rs <- sqrt(rejection_rate_rs * (1 - rejection_rate_rs) / n_converged)
#
#   # Bias
#   bias_glm <- mean(results_df$glm_estimate - treatment_effect, na.rm = TRUE)
#   bias_rs <- mean(results_df$rs_estimate - treatment_effect, na.rm = TRUE)
#
#   # Empirical SE
#   empse_glm <- sqrt(var(results_df$glm_estimate, na.rm = TRUE))
#   empse_rs <- sqrt(var(results_df$rs_estimate, na.rm = TRUE))
#
#   # MSE
#   mse_glm <- mean((results_df$glm_estimate - treatment_effect)^2, na.rm = TRUE)
#   mse_rs <- mean((results_df$rs_estimate - treatment_effect)^2, na.rm = TRUE)
#
#   # RMSE
#   rmse_glm <- sqrt(mse_glm)
#   rmse_rs <- sqrt(mse_rs)
#
#   # Format summary statistics as a data frame
#   summary_df <- data.frame(
#     Model = c("GLM", "R", "RS"),
#     Rejection_Rate = c(rejection_rate_glm, rejection_rate_r, rejection_rate_rs),
#     Rejection_Rate_SE = c(se_rejection_rate_glm, se_rejection_rate_r, se_rejection_rate_rs),
#     Bias = c(bias_glm, NA, bias_rs),
#     EmpSE = c(empse_glm, NA, empse_rs),
#     MSE = c(mse_glm, NA, mse_rs),
#     RMSE = c(rmse_glm, NA, rmse_rs),
#     Converged_Proportion = rep(n_converged / n_sims, 3)
#   )
#
#   # Round numeric values to 4 decimal places
#   summary_df[, -1] <- round(summary_df[, -1], 4)
#
#   return(list(results = results_df, summary = summary_df))
# }



# sim_rstar_glm <- function(n_main, n_covariates, true_coef_main,
#                           n_control = NULL, true_coef_control = NULL,
#                           treatment_effect = NULL,
#                           model = c("logistic", "linear", "poisson"),
#                           skewness_main = NULL, skewness_control = NULL,
#                           Sigma_main = NULL, Sigma_control = NULL, ...) {
#   model <- match.arg(model)
#
#   # Default Sigma to identity matrix if not provided
#   if (is.null(Sigma_main)) {
#     Sigma_main <- diag(n_covariates)
#   }
#   if (is.null(Sigma_control)) {
#     Sigma_control <- diag(n_covariates)
#   }
#
#   # Helper function to generate covariate data
#   generate_covariate_data <- function(n, n_covariates, Sigma, skewness = NULL) {
#     if (n_covariates == 0) {
#       return(matrix(nrow = n, ncol = 0))
#     }
#     if (is.null(skewness)) {
#       X <- tryCatch({
#         MASS::mvrnorm(n, mu = rep(0, n_covariates), Sigma = Sigma)
#       }, error = function(e) {
#         warning("Covariate data generation failed: ", e$message)
#         return(NULL)
#       })
#     } else {
#       cpM <- list(mean = rep(0, n_covariates), var.cov = Sigma, gamma1 = skewness)
#       dpM <- sn::cp2dp(cpM, family = "SN")
#       X <- tryCatch({
#         sn::rmsn(n, dp = dpM)
#       }, error = function(e) {
#         warning("Covariate data generation with skewness failed: ", e$message)
#         return(NULL)
#       })
#     }
#     return(X)
#   }
#
#   # Helper function to generate outcome variable
#   generate_outcome <- function(X, true_coef, model, treatment_effect = NULL, group = NULL) {
#     if (is.null(X)) return(NULL)
#     if (ncol(X) == 0) {
#       eta <- rep(0, nrow(X))
#     } else {
#       eta <- X %*% true_coef
#     }
#     if (!is.null(treatment_effect) & !is.null(group)) {
#       eta <- eta + group * treatment_effect
#     }
#     y <- tryCatch({
#       if (model == "logistic") {
#         p <- stats::plogis(eta)
#         stats::rbinom(nrow(X), size = 1, prob = p)
#       } else if (model == "linear") {
#         eta + stats::rnorm(nrow(X), mean = 0, sd = 1)
#       } else if (model == "poisson") {
#         mu <- exp(eta)
#         stats::rpois(nrow(X), lambda = mu)
#       }
#     }, error = function(e) {
#       warning("Outcome data generation failed: ", e$message)
#       return(NULL)
#     })
#     return(y)
#   }
#
#   # Generate data for main group
#   X_main <- generate_covariate_data(n_main, n_covariates, Sigma_main, skewness_main)
#   y_main <- generate_outcome(X_main, true_coef_main, model, treatment_effect, group = 1)
#   if (is.null(y_main)) return(NULL)
#   data_main <- data.frame(y = y_main, X_main)
#   colnames(data_main) <- c("y", paste0("X", 1:n_covariates))
#
#   # Check if control group data is needed
#   if (!is.null(n_control) && !is.null(true_coef_control)) {
#     X_control <- generate_covariate_data(n_control, n_covariates, Sigma_control, skewness_control)
#     y_control <- generate_outcome(X_control, true_coef_control, model, treatment_effect, group = 0)
#     if (is.null(y_control)) return(NULL)
#     data_control <- data.frame(y = y_control, X_control)
#     colnames(data_control) <- c("y", paste0("X", 1:n_covariates))
#
#     # Combine main and control data, add group indicator
#     data_control$group <- 0
#     data_main$group <- 1
#     data <- rbind(data_main, data_control)
#
#     # Update formula to include group indicator
#     formula <- as.formula(paste("y ~ group +", paste(colnames(data)[-c(1, ncol(data))], collapse = " + ")))
#   } else {
#     data <- data_main
#     formula <- as.formula(paste("y ~", paste(colnames(data)[-1], collapse = " + ")))
#   }
#
#   # Fit the models
#   fit_glm <- tryCatch({
#     stats::glm(formula = formula, family = if (model == "logistic") stats::binomial() else if (model == "linear") stats::gaussian() else stats::poisson(), data = data)
#   }, error = function(e) {
#     warning("GLM model fitting failed: ", e$message)
#     return(NULL)
#   })
#
#   fit_rstar <- tryCatch({
#     rstar_glm(.formula = formula, .data = data, .model = model, ...)
#   }, error = function(e) {
#     warning("rstar model fitting failed: ", e$message)
#     return(NULL)
#   })
#
#   if (is.null(fit_glm) || is.null(fit_rstar)) return(NULL)
#
#   return(list(rstar = fit_rstar, fit_glm = fit_glm, data = data))
# }
#
# run_sim_rstar_glm <- function(n_sims, alpha_level = 0.05,
#                               n_main, n_covariates, true_coef_main,
#                               n_control = NULL, true_coef_control = NULL,
#                               treatment_effect = NULL,
#                               model = c("logistic", "linear", "poisson"),
#                               skewness_main = NULL, skewness_control = NULL,
#                               Sigma_main = NULL, Sigma_control = NULL, ...) {
#
#   # Default Sigma to identity matrix if not provided
#   if (is.null(Sigma_main)) {
#     Sigma_main <- diag(n_covariates)
#   }
#   if (is.null(Sigma_control)) {
#     Sigma_control <- diag(n_covariates)
#   }
#
#   results <- replicate(n_sims, {
#     sim_result <- sim_rstar_glm(n_main, n_covariates, true_coef_main,
#                                 n_control, true_coef_control,
#                                 treatment_effect, model,
#                                 skewness_main, skewness_control,
#                                 Sigma_main, Sigma_control, ...)
#     if (is.null(sim_result)) return(rep(NA, 6))
#     glm_se <- summary(sim_result$fit_glm)$coefficients["group", "Std. Error"]
#     c(
#       glm_estimate = sim_result$fit_glm$coefficients["group"],
#       glm_p_value = summary(sim_result$fit_glm)$coefficients["group", "Pr(>|z|)"],
#       rs_estimate = sim_result$rstar$rs$theta.hat["group"],
#       r_p_value = stats::pnorm(sim_result$rstar$rs$r),
#       rs_p_value = stats::pnorm(sim_result$rstar$rs$rs),
#       glm_se = glm_se
#     )
#   }, simplify = "matrix")
#
#   results_df <- as.data.frame(t(results))
#   converged <- complete.cases(results_df)
#   results_df <- results_df[converged, ]
#
#   # Calculate summary statistics
#   n_converged <- sum(converged)
#   rejection_rate_glm <- mean(results_df$glm_p_value < alpha_level, na.rm = TRUE)
#   rejection_rate_r <- mean(results_df$r_p_value < alpha_level, na.rm = TRUE)
#   rejection_rate_rs <- mean(results_df$rs_p_value < alpha_level, na.rm = TRUE)
#
#   # Standard error of rejection rate (using binomial proportion standard error)
#   se_rejection_rate_glm <- sqrt(rejection_rate_glm * (1 - rejection_rate_glm) / n_converged)
#   se_rejection_rate_r <- sqrt(rejection_rate_r * (1 - rejection_rate_r) / n_converged)
#   se_rejection_rate_rs <- sqrt(rejection_rate_rs * (1 - rejection_rate_rs) / n_converged)
#
#   # Bias
#   bias_glm <- mean(results_df$glm_estimate - treatment_effect, na.rm = TRUE)
#   bias_rs <- mean(results_df$rs_estimate - treatment_effect, na.rm = TRUE)
#
#   # Empirical SE
#   empse_glm <- sqrt(var(results_df$glm_estimate, na.rm = TRUE))
#   empse_rs <- sqrt(var(results_df$rs_estimate, na.rm = TRUE))
#
#   # MSE
#   mse_glm <- mean((results_df$glm_estimate - treatment_effect)^2, na.rm = TRUE)
#   mse_rs <- mean((results_df$rs_estimate - treatment_effect)^2, na.rm = TRUE)
#
#   # RMSE
#   rmse_glm <- sqrt(mse_glm)
#   rmse_rs <- sqrt(mse_rs)
#
#   # Format summary statistics as a data frame
#   summary_df <- data.frame(
#     Model = c("GLM", "R", "RS"),
#     Rejection_Rate = c(rejection_rate_glm, rejection_rate_r, rejection_rate_rs),
#     Rejection_Rate_SE = c(se_rejection_rate_glm, se_rejection_rate_r, se_rejection_rate_rs),
#     Bias = c(bias_glm, NA, bias_rs),
#     EmpSE = c(empse_glm, NA, empse_rs),
#     MSE = c(mse_glm, NA, mse_rs),
#     RMSE = c(rmse_glm, NA, rmse_rs),
#     Converged_Proportion = rep(n_converged / n_sims, 3)
#   )
#
#   # Round numeric values to 4 decimal places
#   summary_df[, -1] <- round(summary_df[, -1], 4)
#
#   return(list(results = results_df, summary = summary_df))
# }












# sim_rstar_glm <- function(n_main, n_covariates, true_coef_main,
#                           n_control = NULL, true_coef_control = NULL,
#                           treatment_effect = NULL,
#                           model = c("logistic", "linear", "poisson"),
#                           skewness_main = NULL, skewness_control = NULL,
#                           Sigma_main, Sigma_control, ...) {
#   model <- match.arg(model)
#
#   # Helper function to generate covariate data
#   generate_covariate_data <- function(n, n_covariates, Sigma, skewness = NULL) {
#     if (is.null(skewness)) {
#       X <- tryCatch({
#         MASS::mvrnorm(n, mu = rep(0, n_covariates), Sigma = Sigma)
#       }, error = function(e) {
#         warning("Covariate data generation failed: ", e$message)
#         return(NULL)
#       })
#     } else {
#       cpM <- list(mean = rep(0, n_covariates), var.cov = Sigma, gamma1 = skewness)
#       dpM <- sn::cp2dp(cpM, family = "SN")
#       X <- tryCatch({
#         sn::rmsn(n, dp = dpM)
#       }, error = function(e) {
#         warning("Covariate data generation with skewness failed: ", e$message)
#         return(NULL)
#       })
#     }
#     return(X)
#   }
#
#   # Helper function to generate outcome variable
#   generate_outcome <- function(X, true_coef, model, treatment_effect = NULL, group = NULL) {
#     if (is.null(X)) return(NULL)
#     eta <- X %*% true_coef
#     if (!is.null(treatment_effect) & !is.null(group)) {
#       eta <- eta + group * treatment_effect
#     }
#     y <- tryCatch({
#       if (model == "logistic") {
#         p <- stats::plogis(eta)
#         stats::rbinom(nrow(X), size = 1, prob = p)
#       } else if (model == "linear") {
#         eta + stats::rnorm(nrow(X), mean = 0, sd = 1)
#       } else if (model == "poisson") {
#         mu <- exp(eta)
#         stats::rpois(nrow(X), lambda = mu)
#       }
#     }, error = function(e) {
#       warning("Outcome data generation failed: ", e$message)
#       return(NULL)
#     })
#     return(y)
#   }
#
#   # Generate data for main group
#   X_main <- generate_covariate_data(n_main, n_covariates, Sigma_main, skewness_main)
#   y_main <- generate_outcome(X_main, true_coef_main, model, treatment_effect, group = 1)
#   if (is.null(y_main)) return(NULL)
#   data_main <- data.frame(y = y_main, X_main)
#   colnames(data_main) <- c("y", paste0("X", 1:n_covariates))
#
#   # Check if control group data is needed
#   if (!is.null(n_control) && !is.null(true_coef_control)) {
#     X_control <- generate_covariate_data(n_control, n_covariates, Sigma_control, skewness_control)
#     y_control <- generate_outcome(X_control, true_coef_control, model, treatment_effect, group = 0)
#     if (is.null(y_control)) return(NULL)
#     data_control <- data.frame(y = y_control, X_control)
#     colnames(data_control) <- c("y", paste0("X", 1:n_covariates))
#
#     # Combine main and control data, add group indicator
#     data_control$group <- 0
#     data_main$group <- 1
#     data <- rbind(data_main, data_control)
#
#     # Update formula to include group indicator
#     formula <- as.formula(paste("y ~ group +", paste(colnames(data)[-c(1, ncol(data))], collapse = " + ")))
#   } else {
#     data <- data_main
#     formula <- as.formula(paste("y ~", paste(colnames(data)[-1], collapse = " + ")))
#   }
#
#   # Fit the models
#   fit_glm <- tryCatch({
#     stats::glm(formula = formula, family = if (model == "logistic") stats::binomial() else if (model == "linear") stats::gaussian() else stats::poisson(), data = data)
#   }, error = function(e) {
#     warning("GLM model fitting failed: ", e$message)
#     return(NULL)
#   })
#
#   fit_rstar <- tryCatch({
#     rstar_glm(.formula = formula, .data = data, .model = model, ...)
#   }, error = function(e) {
#     warning("rstar model fitting failed: ", e$message)
#     return(NULL)
#   })
#
#   if (is.null(fit_glm) || is.null(fit_rstar)) return(NULL)
#
#   return(list(rstar = fit_rstar, fit_glm = fit_glm, data = data))
# }
#
# run_sim_rstar_glm <- function(n_sims, alpha_level = 0.05,
#                               n_main, n_covariates, true_coef_main,
#                               n_control = NULL, true_coef_control = NULL,
#                               treatment_effect = NULL,
#                               model = c("logistic", "linear", "poisson"),
#                               skewness_main = NULL, skewness_control = NULL,
#                               Sigma_main, Sigma_control, ...) {
#
#   results <- replicate(n_sims, {
#     sim_result <- sim_rstar_glm(n_main, n_covariates, true_coef_main,
#                                 n_control, true_coef_control,
#                                 treatment_effect, model,
#                                 skewness_main, skewness_control,
#                                 Sigma_main, Sigma_control, ...)
#     if (is.null(sim_result)) return(rep(NA, 6))
#     glm_se <- summary(sim_result$fit_glm)$coefficients["group", "Std. Error"]
#     c(
#       glm_estimate = sim_result$fit_glm$coefficients["group"],
#       glm_p_value = summary(sim_result$fit_glm)$coefficients["group", "Pr(>|z|)"],
#       rs_estimate = sim_result$rstar$rs$theta.hat["group"],
#       r_p_value = stats::pnorm(sim_result$rstar$rs$r),
#       rs_p_value = stats::pnorm(sim_result$rstar$rs$rs),
#       glm_se = glm_se
#     )
#   }, simplify = "matrix")
#
#   results_df <- as.data.frame(t(results))
#   converged <- complete.cases(results_df)
#   results_df <- results_df[converged, ]
#
#   # Calculate summary statistics
#   n_converged <- sum(converged)
#   rejection_rate_glm <- mean(results_df$glm_p_value < alpha_level, na.rm = TRUE)
#   rejection_rate_r <- mean(results_df$r_p_value < alpha_level, na.rm = TRUE)
#   rejection_rate_rs <- mean(results_df$rs_p_value < alpha_level, na.rm = TRUE)
#
#   # Standard error of rejection rate (using binomial proportion standard error)
#   se_rejection_rate_glm <- sqrt(rejection_rate_glm * (1 - rejection_rate_glm) / n_converged)
#   se_rejection_rate_r <- sqrt(rejection_rate_r * (1 - rejection_rate_r) / n_converged)
#   se_rejection_rate_rs <- sqrt(rejection_rate_rs * (1 - rejection_rate_rs) / n_converged)
#
#   # Bias
#   bias_glm <- mean(results_df$glm_estimate - treatment_effect, na.rm = TRUE)
#   bias_rs <- mean(results_df$rs_estimate - treatment_effect, na.rm = TRUE)
#
#   # Empirical SE
#   empse_glm <- sqrt(var(results_df$glm_estimate, na.rm = TRUE))
#   empse_rs <- sqrt(var(results_df$rs_estimate, na.rm = TRUE))
#
#   # MSE
#   mse_glm <- mean((results_df$glm_estimate - treatment_effect)^2, na.rm = TRUE)
#   mse_rs <- mean((results_df$rs_estimate - treatment_effect)^2, na.rm = TRUE)
#
#   # RMSE
#   rmse_glm <- sqrt(mse_glm)
#   rmse_rs <- sqrt(mse_rs)
#
#   # Format summary statistics as a data frame
#   summary_df <- data.frame(
#     Model = c("GLM", "R", "RS"),
#     Rejection_Rate = c(rejection_rate_glm, rejection_rate_r, rejection_rate_rs),
#     Rejection_Rate_SE = c(se_rejection_rate_glm, se_rejection_rate_r, se_rejection_rate_rs),
#     Bias = c(bias_glm, NA, bias_rs),
#     EmpSE = c(empse_glm, NA, empse_rs),
#     MSE = c(mse_glm, NA, mse_rs),
#     RMSE = c(rmse_glm, NA, rmse_rs),
#     Converged_Proportion = rep(n_converged / n_sims, 3)
#   )
#
#   # Round numeric values to 4 decimal places
#   summary_df[, -1] <- round(summary_df[, -1], 4)
#
#   return(list(results = results_df, summary = summary_df))
# }









# run_sim_rstar_glm <- function(n_sims, alpha_level = 0.05,
#                               n_main, n_covariates, true_coef_main,
#                               n_control = NULL, true_coef_control = NULL,
#                               treatment_effect = NULL,
#                               model = c("logistic", "linear", "poisson"),
#                               skewness_main = NULL, skewness_control = NULL,
#                               Sigma_main, Sigma_control, ...) {
#
#   results <- replicate(n_sims, {
#     sim_result <- sim_rstar_glm(n_main, n_covariates, true_coef_main,
#                                 n_control, true_coef_control,
#                                 treatment_effect, model,
#                                 skewness_main, skewness_control,
#                                 Sigma_main, Sigma_control, ...)
#     if (is.null(sim_result)) return(rep(NA, 6))
#     glm_se <- summary(sim_result$fit_glm)$coefficients["group", "Std. Error"]
#     c(
#       glm_estimate = sim_result$fit_glm$coefficients["group"],
#       glm_p_value = summary(sim_result$fit_glm)$coefficients["group", "Pr(>|z|)"],
#       rs_estimate = sim_result$rstar$rs$theta.hat["group"],
#       r_p_value = stats::pnorm(sim_result$rstar$rs$r),
#       rs_p_value = stats::pnorm(sim_result$rstar$rs$rs),
#       glm_se = glm_se
#     )
#   }, simplify = "matrix")
#
#   results_df <- as.data.frame(t(results))
#   converged <- complete.cases(results_df)
#   results_df <- results_df[converged, ]
#
#   # Calculate summary statistics
#   n_converged <- sum(converged)
#   rejection_rate_glm <- mean(results_df$glm_p_value < alpha_level, na.rm = TRUE)
#   rejection_rate_r <- mean(results_df$r_p_value < alpha_level, na.rm = TRUE)
#   rejection_rate_rs <- mean(results_df$rs_p_value < alpha_level, na.rm = TRUE)
#
#   # Standard error of rejection rate (using binomial proportion standard error)
#   se_rejection_rate_glm <- sqrt(rejection_rate_glm * (1 - rejection_rate_glm) / n_converged)
#   se_rejection_rate_r <- sqrt(rejection_rate_r * (1 - rejection_rate_r) / n_converged)
#   se_rejection_rate_rs <- sqrt(rejection_rate_rs * (1 - rejection_rate_rs) / n_converged)
#
#   # Bias
#   bias_glm <- mean(results_df$glm_estimate - treatment_effect, na.rm = TRUE)
#   bias_rs <- mean(results_df$rs_estimate - treatment_effect, na.rm = TRUE)
#
#   # Empirical SE
#   empse_glm <- sqrt(var(results_df$glm_estimate, na.rm = TRUE))
#   empse_rs <- sqrt(var(results_df$rs_estimate, na.rm = TRUE))
#
#   # MSE
#   mse_glm <- mean((results_df$glm_estimate - treatment_effect)^2, na.rm = TRUE)
#   mse_rs <- mean((results_df$rs_estimate - treatment_effect)^2, na.rm = TRUE)
#
#   # RMSE
#   rmse_glm <- sqrt(mse_glm)
#   rmse_rs <- sqrt(mse_rs)
#
#   # Format summary statistics as a data frame
#   summary_df <- data.frame(
#     Statistic = c("Rejection Rate", "Rejection Rate SE", "Bias", "EmpSE", "MSE", "RMSE", "Converged Proportion"),
#     GLM = c(rejection_rate_glm, se_rejection_rate_glm, bias_glm, empse_glm, mse_glm, rmse_glm, n_converged / n_sims),
#     R = c(rejection_rate_r, se_rejection_rate_r, NA, NA, NA, NA, n_converged / n_sims),
#     RS = c(rejection_rate_rs, se_rejection_rate_rs, bias_rs, empse_rs, mse_rs, rmse_rs, n_converged / n_sims)
#   )
#
#   # Round numeric values to 4 decimal places
#   summary_df[, 2:4] <- round(summary_df[, 2:4], 4)
#
#   return(list(results = results_df, summary = summary_df))
# }








# sim_rstar_glm <- function(n_main, n_covariates, true_coef_main,
#                           n_control = NULL, true_coef_control = NULL,
#                           treatment_effect = NULL,
#                           model = c("logistic", "linear", "poisson"),
#                           skewness_main = NULL, skewness_control = NULL,
#                           Sigma_main, Sigma_control, ...) {
#   model <- match.arg(model)
#
#   # Helper function to generate covariate data
#   generate_covariate_data <- function(n, n_covariates, Sigma, skewness = NULL) {
#     if (is.null(skewness)) {
#       X <- MASS::mvrnorm(n, mu = rep(0, n_covariates), Sigma = Sigma)
#     } else {
#       cpM <- list(mean = rep(0, n_covariates), var.cov = Sigma, gamma1 = skewness)
#       dpM <- sn::cp2dp(cpM, family = "SN")
#       X <- sn::rmsn(n, dp = dpM)
#     }
#     return(X)
#   }
#
#   # Helper function to generate outcome variable
#   generate_outcome <- function(X, true_coef, model, treatment_effect = NULL, group = NULL) {
#     eta <- X %*% true_coef
#     if (!is.null(treatment_effect) & !is.null(group)) {
#       eta <- eta + group * treatment_effect
#     }
#     if (model == "logistic") {
#       p <- stats::plogis(eta)
#       y <- stats::rbinom(nrow(X), size = 1, prob = p)
#     } else if (model == "linear") {
#       y <- eta + stats::rnorm(nrow(X), mean = 0, sd = 1)
#     } else if (model == "poisson") {
#       mu <- exp(eta)
#       y <- stats::rpois(nrow(X), lambda = mu)
#     }
#     return(y)
#   }
#
#   # Generate data for main group
#   X_main <- generate_covariate_data(n_main, n_covariates, Sigma_main, skewness_main)
#   y_main <- generate_outcome(X_main, true_coef_main, model, treatment_effect, group = 1)
#   data_main <- data.frame(y = y_main, X_main)
#   colnames(data_main) <- c("y", paste0("X", 1:n_covariates))
#
#   # Check if control group data is needed
#   if (!is.null(n_control) && !is.null(true_coef_control)) {
#     X_control <- generate_covariate_data(n_control, n_covariates, Sigma_control, skewness_control)
#     y_control <- generate_outcome(X_control, true_coef_control, model, treatment_effect, group = 0)
#     data_control <- data.frame(y = y_control, X_control)
#     colnames(data_control) <- c("y", paste0("X", 1:n_covariates))
#
#     # Combine main and control data, add group indicator
#     data_control$group <- 0
#     data_main$group <- 1
#     data <- rbind(data_main, data_control)
#
#     # Update formula to include group indicator
#     formula <- as.formula(paste("y ~ group +", paste(colnames(data)[-c(1, ncol(data))], collapse = " + ")))
#   } else {
#     data <- data_main
#     formula <- as.formula(paste("y ~", paste(colnames(data)[-1], collapse = " + ")))
#   }
#
#   # Fit the models
#   fit_glm <- stats::glm(formula = formula, family = if (model == "logistic") stats::binomial() else if (model == "linear") stats::gaussian() else stats::poisson(), data = data)
#   fit_rstar <- rstar_glm(.formula = formula, .data = data, .model = model, ...)
#
#   return(list(rstar = fit_rstar, data = data))
# }
#
# run_sim_rstar_glm <- function(n_sims, alpha_level = 0.05,
#                               n_main, n_covariates, true_coef_main,
#                               n_control = NULL, true_coef_control = NULL,
#                               treatment_effect = NULL,
#                               model = c("logistic", "linear", "poisson"),
#                               skewness_main = NULL, skewness_control = NULL,
#                               Sigma_main, Sigma_control, ...) {
#
#   results <- replicate(n_sims, {
#     sim_result <- sim_rstar_glm(n_main, n_covariates, true_coef_main,
#                                 n_control, true_coef_control,
#                                 treatment_effect, model,
#                                 skewness_main, skewness_control,
#                                 Sigma_main, Sigma_control, ...)
#     glm_se <- summary(sim_result$rstar$fit_glm)$coefficients["group", "Std. Error"]
#     c(
#       glm_estimate = sim_result$rstar$fit_glm$coefficients["group"],
#       glm_p_value = summary(sim_result$rstar$fit_glm)$coefficients["group", "Pr(>|z|)"],
#       rs_estimate = sim_result$rstar$rs$theta.hat["group"],
#       r_p_value = stats::pnorm(sim_result$rstar$rs$r),
#       rs_p_value = stats::pnorm(sim_result$rstar$rs$rs),
#       glm_se = glm_se
#     )
#   }, simplify = "matrix")
#
#   results_df <- as.data.frame(t(results))
#
#   # Calculate summary statistics
#   rejection_rate_glm <- mean(results_df$glm_p_value < alpha_level)
#   rejection_rate_r <- mean(results_df$r_p_value < alpha_level)
#   rejection_rate_rs <- mean(results_df$rs_p_value < alpha_level)
#
#   # Standard error of rejection rate (using binomial proportion standard error)
#   se_rejection_rate_glm <- sqrt(rejection_rate_glm * (1 - rejection_rate_glm) / n_sims)
#   se_rejection_rate_r <- sqrt(rejection_rate_r * (1 - rejection_rate_r) / n_sims)
#   se_rejection_rate_rs <- sqrt(rejection_rate_rs * (1 - rejection_rate_rs) / n_sims)
#
#   # Bias
#   bias_glm <- mean(results_df$glm_estimate - treatment_effect)
#   bias_rs <- mean(results_df$rs_estimate - treatment_effect)
#
#   # Empirical SE
#   empse_glm <- sqrt(var(results_df$glm_estimate))
#   empse_rs <- sqrt(var(results_df$rs_estimate))
#
#   # Relative % increase in precision
#   relative_precision <- 100 * ((empse_glm^2 / empse_rs^2) - 1)
#
#   # MSE
#   mse_glm <- mean((results_df$glm_estimate - treatment_effect)^2)
#   mse_rs <- mean((results_df$rs_estimate - treatment_effect)^2)
#
#   # RMSE
#   rmse_glm <- sqrt(mse_glm)
#   rmse_rs <- sqrt(mse_rs)
#
#   # Relative RMSE (relative to the true treatment effect)
#   rel_rmse_glm <- rmse_glm / abs(treatment_effect)
#   rel_rmse_rs <- rmse_rs / abs(treatment_effect)
#
#   # Average ModSE (using GLM model-based SE)
#   modse_glm <- mean(results_df$glm_se)
#
#   # Relative % error in ModSE
#   relative_modse_error <- 100 * ((modse_glm / empse_glm) - 1)
#
#   # Format summary statistics as a data frame
#   summary_df <- data.frame(
#     Statistic = c("Rejection Rate", "Rejection Rate SE", "Bias", "EmpSE", "Relative Precision", "MSE", "RMSE", "Relative RMSE", "ModSE", "Relative ModSE Error"),
#     GLM = c(rejection_rate_glm, se_rejection_rate_glm, bias_glm, empse_glm, relative_precision, mse_glm, rmse_glm, rel_rmse_glm, modse_glm, relative_modse_error),
#     R = c(rejection_rate_r, se_rejection_rate_r, NA, NA, NA, NA, NA, NA, NA, NA),
#     RS = c(rejection_rate_rs, se_rejection_rate_rs, bias_rs, empse_rs, relative_precision, mse_rs, rmse_rs, rel_rmse_rs, NA, NA)
#   )
#
#   # Round numeric values to 4 decimal places
#   summary_df[, 2:4] <- round(summary_df[, 2:4], 4)
#
#   return(list(results = results_df, summary = summary_df))
# }










#
# run_sim_rstar_glm <- function(n_sims, alpha_level = 0.05,
#                               n_main, n_covariates, true_coef_main,
#                               n_control = NULL, true_coef_control = NULL,
#                               treatment_effect = NULL,
#                               model = c("logistic", "linear", "poisson"),
#                               skewness_main = NULL, skewness_control = NULL,
#                               Sigma_main, Sigma_control, ...) {
#
#   results <- replicate(n_sims, {
#     sim_result <- sim_rstar_glm(n_main, n_covariates, true_coef_main,
#                                 n_control, true_coef_control,
#                                 treatment_effect, model,
#                                 skewness_main, skewness_control,
#                                 Sigma_main, Sigma_control, ...)
#     c(
#       glm_estimate = sim_result$rstar$fit_glm$coefficients["group"],
#       glm_p_value = summary(sim_result$rstar$fit_glm)$coefficients["group", "Pr(>|z|)"],
#       rs_estimate = sim_result$rstar$rs$theta.hat["group"],
#       r_p_value = stats::pnorm(sim_result$rstar$rs$r),
#       rs_p_value = stats::pnorm(sim_result$rstar$rs$rs)
#     )
#   }, simplify = "matrix")
#
#   results_df <- as.data.frame(t(results))
#
#   # Calculate summary statistics
#   rejection_rate_glm <- mean(results_df$glm_p_value < alpha_level)
#   rejection_rate_r <- mean(results_df$r_p_value < alpha_level)
#   rejection_rate_rs <- mean(results_df$rs_p_value < alpha_level)
#
#   # Standard error of rejection rate (using binomial proportion standard error)
#   se_rejection_rate_glm <- sqrt(rejection_rate_glm * (1 - rejection_rate_glm) / n_sims)
#   se_rejection_rate_r <- sqrt(rejection_rate_r * (1 - rejection_rate_r) / n_sims)
#   se_rejection_rate_rs <- sqrt(rejection_rate_rs * (1 - rejection_rate_rs) / n_sims)
#
#   # Bias
#   bias_glm <- mean(results_df$glm_estimate - treatment_effect)
#   bias_rs <- mean(results_df$rs_estimate - treatment_effect)
#
#   # RMSE
#   rmse_glm <- sqrt(mean((results_df$glm_estimate - treatment_effect)^2))
#   rmse_rs <- sqrt(mean((results_df$rs_estimate - treatment_effect)^2))
#
#   # Relative RMSE (relative to the true treatment effect)
#   rel_rmse_glm <- rmse_glm / abs(treatment_effect)
#   rel_rmse_rs <- rmse_rs / abs(treatment_effect)
#
#   # Format summary statistics as a data frame
#   summary_df <- data.frame(
#     Statistic = c("Rejection Rate", "Rejection Rate SE", "Bias", "RMSE", "Relative RMSE"),
#     GLM = c(rejection_rate_glm, se_rejection_rate_glm, bias_glm, rmse_glm, rel_rmse_glm),
#     R = c(rejection_rate_r, se_rejection_rate_r, NA, NA, NA),
#     RS = c(rejection_rate_rs, se_rejection_rate_rs, bias_rs, rmse_rs, rel_rmse_rs)
#   )
#
#   # Round numeric values to 4 decimal places
#   summary_df[, 2:4] <- round(summary_df[, 2:4], 4)
#
#   return(list(results = results_df, summary = summary_df))
# }



