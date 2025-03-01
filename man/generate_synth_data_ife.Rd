% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synthetic_data.R
\name{generate_synth_data_ife}
\alias{generate_synth_data_ife}
\title{Synthetic Data with Interactive Fixed Effect Model}
\usage{
generate_synth_data_ife(
  t0 = 99,
  n = 20,
  tau = 10,
  phi = NULL,
  mu = NULL,
  alpha = NULL,
  nu = NULL,
  X = NULL,
  beta = NULL,
  epsilon_sd = 1,
  J = NULL
)
}
\arguments{
\item{t0}{Numeric, the number of pre-treatment time periods (observations).}

\item{n}{Numeric, the number of donor units (variables).}

\item{tau}{Numeric, the value of ATT (Average Treatment Effect on the Treated).}

\item{phi}{Numeric matrix, \eqn{((n+1) \times J)} matrix of unit fixed effects.}

\item{mu}{Numeric matrix, \eqn{(J \times (t0+1))} matrix of time fixed effects.}

\item{alpha}{Numeric vector, \eqn{(n+1)}-length vector of unit fixed effects.}

\item{nu}{Numeric vector, \eqn{(t0+1)}-length vector of time fixed effects.}

\item{X}{Numeric array, \eqn{((n+1) \times (t0+1) \times k)} \eqn{k}-dimensional time varying covariates (optional).}

\item{beta}{Numeric vector, \eqn{k}-length vector of coefficients for time varying covariates.}

\item{epsilon_sd}{Numeric, standard deviation of error term.}

\item{J}{Numeric, dimension of factor (interacted fixed effects).}
}
\value{
A list containing the following elements:
\itemize{
\item \code{df_pre_full}: A tibble of complete pre-treatment data \eqn{(t0 \times (n+1))}. The last column is the observation for treated unit (\eqn{Y}).
\item \code{df_prepost_full}: A tibble of complete pre- and post-treatment data \eqn{((t0+1) \times (n+1))}. The second to last column is the observation for treated unit (\eqn{Y}), and the last column is the treatment indicator (\eqn{D}). The last row is the post-treatment observation.
\item \code{tau}: A numeric value of ATT
\item \code{params}: A list of parameters used to generate the data.
\item \code{.call}: The matched call.
}
}
\description{
This function generates synthetic data based on an interactive fixed effect model.
}
\examples{
set.seed(123)
synth <- generate_synth_data_ife(
  t0 = 99,
  n = 20,
  tau = 10,
  phi = NULL,
  mu = NULL,
  alpha = NULL,
  nu = NULL,
  X = NULL,
  beta = NULL,
  epsilon_sd = 1,
  J = 2
)

}
