#' Compute r* Statistics for Generalized Linear Models
#'
#' The `rstar_glm` function computes r* statistics for hypothesis testing
#' on coefficients of interest in generalized linear models (GLMs).
#' It supports logistic, linear, and Poisson regression models. For logistic
#' models, the outcome must be binary.
#'
#' @param .formula A formula specifying the model.
#' @param .data A data frame containing the variables in the model.
#' @param .model The type of GLM model: "logistic", "linear", or "poisson".
#' @param .psidesc A description of the parameter of interest.
#' @param .psival The value of the parameter of interest under the null hypothesis.
#' @param .fpsi The index of the parameter of interest.
#' @param .rstar.ci Logical; if TRUE, compute confidence intervals for r*.
#' @param ... Additional arguments passed to the likelihoodAsy functions.
#'
#' @return A list with the object returned from likelihoodAsy::rstar (`rs`),
#' the object returned from likelihoodAsy::rstar.ci (`rs_ci`), and the object
#' returned from stats::glm (`fit_glm`).
#'
#' @examples
#'
#' # Logistic model
#' rstar_glm(law ~ DriversKilled + VanKilled + drivers + kms,
#'           .data = Seatbelts,
#'           .model = "logistic") |> suppressWarnings()
#'
#' # Poisson model
#' rstar_glm(count ~ spray,
#'           .data = InsectSprays,
#'           .model = "poisson") |> suppressWarnings()
#'
#' # Linear model
#' rstar_glm(mpg ~ wt + hp,
#'           .data = mtcars,
#'           .model = "linear") |> suppressWarnings()
#'
#' @references
#' Pierce, D. A., & Bellio, R. (2017). Modern Likelihood-Frequentist Inference.
#' International Statistical Review / Revue Internationale de Statistique, 85(3),
#' 519–541. <doi:10.1111/insr.12232>
#'
#' Bellio R, Pierce D (2020). likelihoodAsy: Functions for Likelihood Asymptotics.
#' R package version 0.51, \url{https://CRAN.R-project.org/package=likelihoodAsy}.
#'
#' @export
rstar_glm <- function(.formula, .data, .model = c("logistic", "linear", "poisson"),
                      .psidesc = "Coefficient of Interest", .psival = 0, .fpsi = 2,
                      .rstar.ci = FALSE, ...) {
  UseMethod("rstar_glm", .model)
}

#' @rdname rstar_glm
#' @export
rstar_glm.logistic <- function(.formula, .data, .model = c("logistic", "linear", "poisson"),
                               .psidesc = "Coefficient of Interest", .psival = 0, .fpsi = 2,
                               .rstar.ci = FALSE, ...) {

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
  # Get confidence intervals for r*
  if (.rstar.ci){
    rs_ci <- likelihoodAsy::rstar.ci(data = data_obj,
                                     thetainit = stats::coef(fit_glm),
                                     floglik = loglik_lgr,
                                     fpsi = function(theta) theta[.fpsi],
                                     fscore = grad_lgr,
                                     datagen = sim_lgr,
                                     psidesc = .psidesc,
                                     ...)

  } else {
    rs_ci <- NULL
  }


  # Build return object
  ret <- list(rs = rs, rs_ci = rs_ci, fit_glm = fit_glm)
  class(ret) <- "rstar_glm_result"
  return(ret)
}

#' @rdname rstar_glm
#' @export
rstar_glm.linear <- function(.formula, .data, .model = c("logistic", "linear", "poisson"),
                             .psidesc = "Coefficient of Interest", .psival = 0, .fpsi = 2,
                             .rstar.ci = FALSE, ...) {
  # Fit linear regression model with Gaussian link
  fit_glm <- stats::glm(formula = .formula, family = stats::gaussian, data = .data)

  # Build data object
  data_obj <- list(y = fit_glm$y, X = stats::model.matrix(fit_glm))

  loglik_lin <- function(theta, data) {
    y <- data$y
    X <- data$X
    N <- length(y)
    eta <- X %*% theta
    sigma2 <- mean((y - eta)^2)  # MLE of sigma^2

    log_likelihood <- -0.5 * N * log(2 * pi) - 0.5 * N * log(sigma2) - 0.5 / sigma2 * sum((y - eta)^2)

    return(log_likelihood)
  }

  grad_lin <- function(theta, data) {
    y <- data$y
    X <- data$X
    eta <- X %*% theta
    sigma2 <- mean((y - eta)^2)  # MLE of sigma^2

    out <- (1 / sigma2) * t(X) %*% (y - eta)
    return(drop(out))
  }

  sim_lin <- function(theta, data) {
    X <- data$X
    eta <- X %*% theta
    sigma2 <- mean((data$y - eta)^2)  # MLE of sigma^2
    out <- data
    out$y <- stats::rnorm(length(data$y), mean = eta, sd = sqrt(sigma2))
    return(out)
  }

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

  # Get confidence intervals for r*
  if (.rstar.ci){
    rs_ci <- likelihoodAsy::rstar.ci(data = data_obj,
                                     thetainit = stats::coef(fit_glm),
                                     floglik = loglik_lin,
                                     fpsi = function(theta) theta[.fpsi],
                                     fscore = grad_lin,
                                     datagen = sim_lin,
                                     psidesc = .psidesc,
                                     ...)

  } else {
    rs_ci <- NULL
  }

  # Build return object
  ret <- list(rs = rs, rs_ci = rs_ci, fit_glm = fit_glm)
  class(ret) <- "rstar_glm_result"
  return(ret)
}

#' @rdname rstar_glm
#' @export
rstar_glm.poisson <- function(.formula, .data, .model = c("logistic", "linear", "poisson"),
                              .psidesc = "Coefficient of Interest", .psival = 0, .fpsi = 2,
                              .rstar.ci = FALSE, ...) {
  # Fit Poisson regression model
  fit_glm <- stats::glm(formula = .formula, family = stats::poisson, data = .data)

  # Build data object
  data_obj <- list(y = fit_glm$y, X = stats::model.matrix(fit_glm))

  # Create log-likelihood function
  loglik_pois <- function(theta, data) {
    y <- data$y
    X <- data$X
    eta <- X %*% theta
    mu <- exp(eta)
    l <- sum(y * log(mu) - mu - log(factorial(y)))
    return(l)
  }

  # Create gradient function (score function)
  grad_pois <- function(theta, data) {
    y <- data$y
    X <- data$X
    eta <- X %*% theta
    mu <- exp(eta)
    out <- t(X) %*% (y - mu)
    return(drop(out))
  }

  # Create data generation function
  sim_pois <- function(theta, data) {
    X <- data$X
    eta <- X %*% theta
    mu <- exp(eta)
    out <- data
    out$y <- stats::rpois(length(data$y), lambda = mu)
    return(out)
  }

  # Compute the r* statistic for the hypothesis test on the specified parameter
  rs <- likelihoodAsy::rstar(data = data_obj,
                             thetainit = stats::coef(fit_glm),
                             floglik = loglik_pois,
                             fpsi = function(theta) theta[.fpsi],
                             psival = .psival,
                             fscore = grad_pois,
                             datagen = sim_pois,
                             psidesc = .psidesc,
                             ...)

  # Get confidence intervals for r*
  if (.rstar.ci){
    rs_ci <- likelihoodAsy::rstar.ci(data = data_obj,
                                     thetainit = stats::coef(fit_glm),
                                     floglik = loglik_pois,
                                     fpsi = function(theta) theta[.fpsi],
                                     fscore = grad_pois,
                                     datagen = sim_pois,
                                     psidesc = .psidesc,
                                     ...)

  } else {
    rs_ci <- NULL
  }

  # Build return object
  ret <- list(rs = rs, rs_ci = rs_ci, fit_glm = fit_glm)
  class(ret) <- "rstar_glm_result"
  return(ret)
}

#' @rdname rstar_glm
#' @export
rstar_glm.default <- function(.formula, .data, .model = c("logistic", "linear", "poisson"),
                              .psidesc = "Coefficient of Interest", .psival = 0, .fpsi = 2,
                              .rstar.ci = FALSE, ...) {
    .model <- match.arg(.model)
    method <- paste("rstar_glm", .model, sep = ".")
    if (!exists(method, mode = "function")) {
      stop("Unsupported model type")
    }
    # do.call(method, list(.formula = .formula, .data = .data, ...))
    do.call(method, list(.formula = .formula, .data = .data, .model = .model,
                         .psidesc = .psidesc, .psival = .psival, .fpsi = .fpsi,
                         .rstar.ci = .rstar.ci, ...))
}
