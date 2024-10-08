% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstar_glm.R
\name{rstar_glm}
\alias{rstar_glm}
\alias{rstar_glm.logistic}
\alias{rstar_glm.linear}
\alias{rstar_glm.poisson}
\alias{rstar_glm.default}
\title{Compute r* Statistics for Generalized Linear Models}
\usage{
rstar_glm(
  .formula,
  .data,
  .model = c("logistic", "linear", "poisson"),
  .psidesc = "Coefficient of Interest",
  .psival = 0,
  .fpsi = 2,
  .rstar.ci = FALSE,
  trace = FALSE,
  ...
)

\method{rstar_glm}{logistic}(
  .formula,
  .data,
  .model = c("logistic", "linear", "poisson"),
  .psidesc = "Coefficient of Interest",
  .psival = 0,
  .fpsi = 2,
  .rstar.ci = FALSE,
  trace = FALSE,
  ...
)

\method{rstar_glm}{linear}(
  .formula,
  .data,
  .model = c("logistic", "linear", "poisson"),
  .psidesc = "Coefficient of Interest",
  .psival = 0,
  .fpsi = 2,
  .rstar.ci = FALSE,
  trace = FALSE,
  ...
)

\method{rstar_glm}{poisson}(
  .formula,
  .data,
  .model = c("logistic", "linear", "poisson"),
  .psidesc = "Coefficient of Interest",
  .psival = 0,
  .fpsi = 2,
  .rstar.ci = FALSE,
  trace = FALSE,
  ...
)

\method{rstar_glm}{default}(
  .formula,
  .data,
  .model = c("logistic", "linear", "poisson"),
  .psidesc = "Coefficient of Interest",
  .psival = 0,
  .fpsi = 2,
  .rstar.ci = FALSE,
  trace = FALSE,
  ...
)
}
\arguments{
\item{.formula}{A formula specifying the model.}

\item{.data}{A data frame containing the variables in the model.}

\item{.model}{The type of GLM model: "logistic", "linear", or "poisson".}

\item{.psidesc}{A description of the parameter of interest.}

\item{.psival}{The value of the parameter of interest under the null hypothesis.}

\item{.fpsi}{The index of the parameter of interest.}

\item{.rstar.ci}{Logical; if TRUE, compute confidence intervals for r*.}

\item{trace}{Logical; if TRUE, print information about computation. (Default is FALSE)}

\item{...}{Additional arguments passed to the likelihoodAsy functions.}
}
\value{
A list with the object returned from likelihoodAsy::rstar (\code{rs}),
the object returned from likelihoodAsy::rstar.ci (\code{rs_ci}), and the object
returned from stats::glm (\code{fit_glm}).
}
\description{
The \code{rstar_glm} function computes r* statistics for hypothesis testing
on coefficients of interest in generalized linear models (GLMs).
It supports logistic, linear, and Poisson regression models. For logistic
models, the outcome must be binary.
}
\examples{

# Logistic model
rstar_glm(law ~ DriversKilled + VanKilled + drivers + kms,
          .data = Seatbelts,
          .model = "logistic") |> suppressWarnings()

# Poisson model
rstar_glm(count ~ spray,
          .data = InsectSprays,
          .model = "poisson") |> suppressWarnings()

# Linear model
rstar_glm(mpg ~ wt + hp,
          .data = mtcars,
          .model = "linear") |> suppressWarnings()

}
\references{
Pierce, D. A., & Bellio, R. (2017). Modern Likelihood-Frequentist Inference.
International Statistical Review / Revue Internationale de Statistique, 85(3),
519–541. \url{doi:10.1111/insr.12232}

Bellio R, Pierce D (2020). likelihoodAsy: Functions for Likelihood Asymptotics.
R package version 0.51, \url{https://CRAN.R-project.org/package=likelihoodAsy}.
}
