rstar_glm <- function(.formula, .data, .model = "logistic", ...) {
  UseMethod("rstar_glm")
}

rstar_glm.logistic <- function(.formula, .data, .psidesc = "Coefficient of Interest",
                               .psival = 0, .fpsi = 2, ...) {

  # Fit logistic regression model with binary outcome
  fit_glm <- stats::glm(formula = .formula, family = stats::binomial, data = .data)

  # Build data object
  data_obj <- list(y = fit_glm$y, X = stats::model.matrix(fit_glm))

  # Create loglikelihood function
  loglik_lgr <- function(theta, data) {

    # Extract response variable and covariates
    y <- data$y
    X <- data$X

    # Calculate linear predictor (eta)
    eta <- X %*% theta

    # Apply logistic transformation to get probabilities (p)
    p <- stats::plogis(eta) # p = exp(eta) / (1 + exp(eta))

    # Compute log likelihood
    # For binary outcomes, the log likelihood is:
    l <- sum(y * log(p) + (1 - y) * log(1 - p))

    # Return the log likelihood
    return(l)
  }

  # Create gradient function
  grad_lgr <- function(theta, data) {

    # Extract response variable and covariates
    y <- data$y
    X <- data$X

    # Calculate linear predictor (eta)
    eta <- X %*% theta

    # Apply logistic transformation to get probabilities (p)
    p <- stats::plogis(eta) # p = exp(eta) / (1 + exp(eta))

    # Compute the gradient (score function)
    # The gradient is the derivative of the log likelihood with respect to the parameters
    # For binary outcomes, the gradient is:
    out <- t(y - p) %*% X

    # Return the gradient, flattened to a vector
    return(drop(out))
  }

  # Create data generation function
  sim_lgr <- function(theta, data) {

    # Extract covariates
    X <- data$X

    # Calculate linear predictor (eta)
    eta <- X %*% theta

    # Apply logistic transformation to get probabilities (p)
    p <- stats::plogis(eta) # p = exp(eta) / (1 + exp(eta))

    # Generate new response variable (y) from a binomial distribution
    # Each observation's response is generated as a binomial random variable
    # with probability p (since size is 1, it's equivalent to a binary outcome)
    out <- data
    out$y <- stats::rbinom(length(data$y), size = 1, prob = p)

    # Return the new dataset
    return(out)
  }

  # Compute the r* statistic for the hypothesis test on the specified parameter
  rs <- likelihoodAsy::rstar(data = data_obj,
                             thetainit = stats::coef(fit_glm),
                             floglik = loglik_lgr,
                             fpsi = function(theta) theta[.fpsi],
                             psival = .psival,
                             fscore = grad_lgr,
                             datagen = sim_lgr,
                             psidesc = .psidesc,
                             ...)

  # Build return object
  ret <- list(rs = rs, fit_glm = fit_glm)
  class(ret) <- "rstar_glm_result"
  return(ret)
}

rstar_glm.linear <- function(.formula, .data, .psidesc = "Coefficient of Interest",
                             .psival = 0, .fpsi = 2, ...) {
  # Fit linear regression model with Gaussian link
  fit_glm <- stats::glm(formula = .formula, family = stats::gaussian, data = .data)

  # Build data object
  data_obj <- list(y = fit_glm$y, X = stats::model.matrix(fit_glm))

  loglik_lin <- function(theta, data) {
    y <- data$y
    X <- data$X
    eta <- X %*% theta
    l <- -0.5 * sum((y - eta)^2)
    return(l)
  }

  grad_lin <- function(theta, data) {
    y <- data$y
    X <- data$X
    eta <- X %*% theta
    out <- t(y - eta) %*% X
    return(drop(out))
  }

  sim_lin <- function(theta, data) {
    X <- data$X
    eta <- X %*% theta
    out <- data
    out$y <- stats::rnorm(length(data$y), mean = eta, sd = sqrt(stats::var(data$y)))
    return(out)
  }

#   loglik_lin <- function(theta, data) {
#     y <- data$y
#     X <- data$X
#     eta <- X %*% theta
#     N <- length(y)
#     sigma2 <- var(y - eta)  # Estimate of sigma^2
#
#     # Compute log-likelihood
#     log_likelihood <- -0.5 * N * log(2 * pi) - 0.5 * N * log(sigma2) - 0.5 / sigma2 * sum((y - eta)^2)
#
#     return(log_likelihood)
#   }
#
#   grad_lin <- function(theta, data) {
#     y <- data$y
#     X <- data$X
#     eta <- X %*% theta
#     sigma2 <- var(y - eta)  # Estimate of sigma^2
#     residuals <- as.matrix(y - eta)  # Ensure residuals is a column vector
#     out <- (1 / sigma2) * t(residuals) %*% X
#     return(drop(out))
#   }
#
#   sim_lin <- function(theta, data) {
#     X <- data$X
#     eta <- X %*% theta
#     sigma2 <- var(data$y - eta)  # Estimate of sigma^2
#     out <- data
#     out$y <- rnorm(length(data$y), mean = eta, sd = sqrt(sigma2))
#     return(out)
#   }

  # Compute the r* statistic for the hypothesis test on the specified parameter
  rs <- likelihoodAsy::rstar(data = data_obj,
                             thetainit = stats::coef(fit_glm),
                             floglik = loglik_lin,
                             fpsi = function(theta) theta[.fpsi],
                             psival = .psival,
                             fscore = grad_lin,
                             datagen = sim_lin,
                             psidesc = .psidesc,
                             ...)

  # Build return object
  ret <- list(rs = rs, fit_glm = fit_glm)
  class(ret) <- "rstar_glm_result"
  return(ret)
}

rstar_glm.default <- function(.formula, .data, .model = "logistic", ...) {
  stop("Unsupported model type")
}

rstar_glm <- function(.formula, .data, .model = "logistic", ...) {
  method <- paste("rstar_glm", .model, sep = ".")
  if (!exists(method, mode = "function")) {
    method <- "rstar_glm.default"
  }
  do.call(method, list(.formula = .formula, .data = .data, ...))
}

