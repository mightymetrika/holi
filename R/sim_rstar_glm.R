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

  # Generate data for control group
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


  # Generate a unique code for the simulation run
  timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
  random_string <- paste(sample(c(letters, LETTERS, 0:9), 5, replace = TRUE), collapse = "")
  run_code <- paste0(timestamp, "_", random_string)

  results <- replicate(n_sims, {
    sim_result <- sim_rstar_glm(n_main, n_covariates, true_coef_main,
                                n_control, true_coef_control,
                                treatment_effect, model,
                                skewness_main, skewness_control,
                                Sigma_main, Sigma_control, ...)
    if (is.null(sim_result)) return(rep(NA, 6))

    glm_coef <- sim_result$fit_glm$coefficients["group"]
    glm_se <- summary(sim_result$fit_glm)$coefficients["group", "Std. Error"]
    if (model %in% c("logistic", "poisson")){
      glm_p_value <- summary(sim_result$fit_glm)$coefficients["group", "Pr(>|z|)"]
      } else if (model == "linear"){
        glm_p_value <- summary(sim_result$fit_glm)$coefficients["group", "Pr(>|t|)"]
      }
    rs_estimate <- sim_result$rstar$rs$theta.hat["group"]

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
    sigma_control = rep(paste(Sigma_control, collapse = ","), 2),
    run_code = c(run_code, run_code)
  )

  # Round numeric values to 4 decimal places
  numeric_columns <- sapply(summary_df, is.numeric)
  summary_df[, numeric_columns] <- round(summary_df[, numeric_columns], 4)

  return(list(results = results_df, summary = summary_df))
}
