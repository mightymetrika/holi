rstar_lr <- function(.formula, .data, .psidesc="Coefficient of Interest",
                     .psival = 0, .fpsi = 3){

  # Fit logistic regression model with binary outcome
  .mod_glm <- stats::glm(formula = .formula,
                         family = stats::binomial,
                         data = .data)

  print(summary(.mod_glm))

  # Build data object
  data_obj <- list(#y = .data[,all.vars(.formula)[1]],
                   y = .mod_glm$y,
                   X = stats::model.matrix(.mod_glm))


  loglik_lgr <- function(theta, data) {
    # Extract response variable and covariates
    y <- data$y
    print(paste("y: ", y))
    X <- data$X

    # Calculate linear predictor (eta)
    eta <- X %*% theta

    # Apply logistic transformation to get probabilities (p)
    p <- stats::plogis(eta)  # p = exp(eta) / (1 + exp(eta))

    print(paste("p: ", p))

    # Compute log likelihood
    # For binary outcomes, the log likelihood is:
    l <- sum(y * log(p) + (1 - y) * log(1 - p))

    # Return the log likelihood
    return(l)
  }

  grad_lgr <- function(theta, data) {
    # Extract response variable and covariates
    y <- data$y
    X <- data$X

    # Calculate linear predictor (eta)
    eta <- X %*% theta

    # Apply logistic transformation to get probabilities (p)
    p <- stats::plogis(eta)  # p = exp(eta) / (1 + exp(eta))

    # Compute the gradient (score function)
    # The gradient is the derivative of the log likelihood with respect to the parameters
    # For binary outcomes, the gradient is:
    #out <- t(y - p) %*% X
    out <- t(y - p) %*% X

    # Return the gradient, flattened to a vector
    return(drop(out))
  }

  sim_lgr <- function(theta, data) {
    # Extract covariates
    X <- data$X

    # Calculate linear predictor (eta)
    eta <- X %*% theta

    # Apply logistic transformation to get probabilities (p)
    p <- stats::plogis(eta)  # p = exp(eta) / (1 + exp(eta))

    # Generate new response variable (y) from a binomial distribution
    # Each observation's response is generated as a binomial random variable
    # with probability p (since size is 1, it's equivalent to a binary outcome)
    out <- data
    out$y <- stats::rbinom(length(data$y), size = 1, prob = p)

    # Return the new dataset
    return(out)
  }

  # Compute the r* statistic for the hypothesis test on the 19th parameter
  rs <- likelihoodAsy::rstar(data=data_obj,
                             thetainit = stats::coef(.mod_glm),
                             floglik = loglik_lgr,
                             fpsi = function(theta) theta[2],
                             psival = 1,
                             #fscore=grad_lgr,
                             datagen=sim_lgr,
                             trace=FALSE,
                             psidesc=.psidesc)

}

