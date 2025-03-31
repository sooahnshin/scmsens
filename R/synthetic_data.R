#' Synthetic Data with Vanilla Vertical Regression
#'
#' This function generates synthetic data based on a vanilla vertical regression model given a set of parameters.
#' This method is designed to analyze data with a specified number of pre-treatment time periods and donor units,
#' incorporating an ATT (Average Treatment Effect on the Treated) value and parameters for non-overlap.
#'
#' @param t0 Numeric, the number of pre-treatment time periods (observations).
#' @param n Numeric, the number of donor units (variables).
#' @param tau Numeric, the value of ATT (Average Treatment Effect on the Treated).
#' @param rho Numeric, a parameter for non-overlap.
#' @param mu Numeric, another parameter for non-overlap.
#' @param beta Numeric vector, weights.
#'
#' @return A list containing the following elements:
#'    - `df_pre_full`: A tibble of complete pre-treatment data \eqn{(t0 \times (n+1))}. The last column is the observation for treated unit (\eqn{Y}).
#'    - `df_prepost_full`: A tibble of complete pre- and post-treatment data \eqn{((t0+1) \times (n+1))}. The second to last column is the observation for treated unit (\eqn{Y}), and the last column is the treatment indicator (\eqn{D}). The last row is the post-treatment observation.
#'    - `beta`: A numeric vector of weights.
#'    - `tau`: A numeric value of ATT
#'    - `.call`: The matched call.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
#' @importFrom tibble as_tibble
#'
#' @examples
#' set.seed(123)
#' synth <- generate_synth_data_vertreg(t0 = 99, n = 20, tau = 10, rho = 0.04, mu = NULL, beta = NULL)
#'
#' @export

generate_synth_data_vertreg <- function(t0 = 99,
                                        n = 20,
                                        tau = 10,
                                        rho = 0.04,
                                        mu = NULL,
                                        beta = NULL) {
  # observed donors
  A <- matrix(rnorm(n^2, 0, rho), ncol = n)
  # to induce correlation among donors
  # -> this may have an impact on the non-overlap term
  sigma <- t(A) %*% A
  diag(sigma) <- 1
  if (is.null(mu)) {
    mu <- rnorm(n, sd = 2)
  }
  # donors are generated from different mean values
  # -> this may have an impact on the non-overlap term
  X0 <- MASS::mvrnorm(n = t0, mu = mu, Sigma = sigma)
  # treatment
  D <- rep(0, t0)
  D <- c(D, 1)

  # True coefficients -> Note that these will not be used in subsequent parts;
  # they are usued solely for data generation to control the distribution of the weights term
  if (is.null(beta)) {
    beta <- rnorm(n, sd = 2)
  }
  epsilon_full <- rnorm(t0)

  # pretreatment outcome
  Y0 <- X0 %*% beta + epsilon_full
  # posttreatment outcome
  x_T <- rnorm(n)
  epsilon_T <- rnorm(1)
  y_T <- beta %*% x_T + tau + epsilon_T

  # labels
  xlab <- paste0("X", 1:n)
  colnames(X0) <- xlab

  df_pre_full <- tibble::as_tibble(cbind(X0, Y0), .name_repair = ~{c(xlab, "Y")})

  df_prepost_full <- df_pre_full |>
    dplyr::mutate(D = 0)
  df_prepost_full <- rbind(
    df_prepost_full,
    c(x_T, y_T, 1)
  )

  res <- list(
    df_pre_full = df_pre_full,
    df_prepost_full = df_prepost_full,
    beta = beta,
    tau = tau,
    .call = match.call()
  )

  return(res)
}

#' Simulate Time Series from AR(l) Model with Specific Initial Values
#'
#' This function simulates a time series from an autoregressive model of order \eqn{l} (AR(\eqn{l}))
#' with specified initial values for each lag.
#'
#' @param phi A numeric vector of length \code{l} containing the autoregressive coefficients
#'   for each lag. The AR(1) coefficient corresponds to \code{phi[1]}, AR(2) to \code{phi[2]},
#'   and so on up to AR(\code{l}).
#' @param y0 A numeric vector of length \code{l} specifying the initial values for each lag
#'   of the time series. The value \code{y0[i]} represents the initial value for lag \eqn{i}.
#' @param n An integer specifying the number of observations to simulate.
#' @param sd_eps A numeric value representing the standard deviation of the innovations.
#'   The innovations are assumed to be independent and identically distributed (iid) Gaussian white noise
#'   with mean 0 and standard deviation \code{sd_eps}.
#'
#' @return A numeric vector containing the simulated time series.
#'
#' @examples
#' # Simulate an AR(2) series with specific initial values
#' phi <- c(0.3, 0.2) # AR(2) coefficients
#' y0 <- c(10, 15) # Initial values for each lag
#' n <- 100 # Number of observations
#' sd_eps <- 1 # Standard deviation of innovations
#' simulated_series <- simulate_ar_l_with_initial(phi, y0, n, sd_eps)
#'
#' # Plot the simulated series
#' plot(simulated_series,
#'   type = "l", main = "Simulated AR(2) Series with Specific Initial Values",
#'   xlab = "Time", ylab = "Value"
#' )
#'
#' @export
#' @seealso [generate_synth_data_ar()]
simulate_ar_l_with_initial <- function(phi, y0, n, sd_eps) {
  l <- length(phi)
  y <- numeric(n)
  y[1:l] <- y0
  eps <- rnorm(n, mean = 0, sd = sd_eps)
  for (t in (l + 1):n) {
    y[t] <- sum(phi * y[(t - l):(t - 1)]) + eps[t]
  }
  return(y)
}


#' Synthetic Data with Autoregressive Model AR(l)
#'
#' This function generates synthetic data based on an autoregressive model of order \eqn{l} (AR(\eqn{l})).
#'
#' @param t0 Numeric, the number of pre-treatment time periods (observations).
#' @param n Numeric, the number of donor units (variables).
#' @param tau Numeric, the value of ATT (Average Treatment Effect on the Treated).
#' @param ar Numeric vector of length \eqn{l}, coefficients for AR model.
#' @param init \eqn{(n+1)} by \eqn{l} matrix, initial values for AR model (optional).
#'
#' @return A list containing the following elements:
#'    - `df_pre_full`: A tibble of complete pre-treatment data \eqn{(t0 \times (n+1))}. The last column is the observation for treated unit (\eqn{Y}).
#'    - `df_prepost_full`: A tibble of complete pre- and post-treatment data \eqn{((t0+1) \times (n+1))}. The second to last column is the observation for treated unit (\eqn{Y}), and the last column is the treatment indicator (\eqn{D}). The last row is the post-treatment observation.
#'    - `tau`: A numeric value of ATT
#'    - `.call`: The matched call.
#'
#' @examples
#' set.seed(123)
#' synth <- generate_synth_data_ar(t0 = 99, n = 20, tau = 10, ar = c(0.7, 0.2, -0.1, -0.3))
#'
#' @importFrom stats arima.sim
#' @importFrom dplyr mutate
#'
#' @export
#' @seealso [simulate_ar_l_with_initial()]
generate_synth_data_ar <- function(
    t0 = 99,
    n = 20,
    tau = 10,
    ar = c(0.7, 0.2, -0.1, -0.3),
    init = NULL) {
  dat <- matrix(NA, nrow = n + 1, ncol = t0 + 1)

  if (!is.null(init)) {
    for (i in 1:(n + 1)) {
      dat[i, ] <- simulate_ar_l_with_initial(ar, init[i, ], t0 + 1, 1)
    }
  } else {
    dat <- t(replicate((n + 1), arima.sim(n = (t0 + 1), list(ar = ar), sd = 1)))
  }
  # labels
  xlab <- paste0("X", 1:n)
  df_pre_full <- tibble::as_tibble(t(dat), .name_repair = ~{c(xlab, "Y")})

  df_prepost_full <- df_pre_full |>
    dplyr::mutate(D = 0)
  df_prepost_full[(t0 + 1), "D"] <- 1
  df_prepost_full[(t0 + 1), "Y"] <- df_prepost_full[(t0 + 1), "Y"] + tau

  df_pre_full <- df_pre_full[-(t0 + 1), ]

  res <- list(
    df_pre_full = df_pre_full,
    df_prepost_full = df_prepost_full,
    tau = tau,
    .call = match.call()
  )

  return(res)
}

#' Synthetic Data with Interactive Fixed Effect Model
#'
#' This function generates synthetic data based on an interactive fixed effect model.
#'
#' @param t0 Numeric, the number of pre-treatment time periods (observations).
#' @param n Numeric, the number of donor units (variables).
#' @param tau Numeric, the value of ATT (Average Treatment Effect on the Treated).
#' @param phi Numeric matrix, \eqn{((n+1) \times J)} matrix of unit fixed effects.
#' @param mu Numeric matrix, \eqn{(J \times (t0+1))} matrix of time fixed effects.
#' @param alpha Numeric vector, \eqn{(n+1)}-length vector of unit fixed effects.
#' @param nu Numeric vector, \eqn{(t0+1)}-length vector of time fixed effects.
#' @param X Numeric array, \eqn{((n+1) \times (t0+1) \times k)} \eqn{k}-dimensional time varying covariates (optional).
#' @param beta Numeric vector, \eqn{k}-length vector of coefficients for time varying covariates.
#' @param epsilon_sd Numeric, standard deviation of error term.
#' @param J Numeric, dimension of factor (interacted fixed effects).
#'
#' @return A list containing the following elements:
#'    - `df_pre_full`: A tibble of complete pre-treatment data \eqn{(t0 \times (n+1))}. The last column is the observation for treated unit (\eqn{Y}).
#'    - `df_prepost_full`: A tibble of complete pre- and post-treatment data \eqn{((t0+1) \times (n+1))}. The second to last column is the observation for treated unit (\eqn{Y}), and the last column is the treatment indicator (\eqn{D}). The last row is the post-treatment observation.
#'    - `tau`: A numeric value of ATT
#'    - `params`: A list of parameters used to generate the data.
#'    - `.call`: The matched call.
#'
#' @importFrom dplyr slice
#'
#' @examples
#' set.seed(123)
#' synth <- generate_synth_data_ife(
#'   t0 = 99,
#'   n = 20,
#'   tau = 10,
#'   phi = NULL,
#'   mu = NULL,
#'   alpha = NULL,
#'   nu = NULL,
#'   X = NULL,
#'   beta = NULL,
#'   epsilon_sd = 1,
#'   J = 2
#' )
#'
#' @export
generate_synth_data_ife <- function(
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
    J = NULL) {
  if (is.null(phi) & is.null(J)) {
    stop("Either phi or J should be specified.")
  }
  if (is.null(mu) & is.null(J)) {
    stop("Either mu or J should be specified.")
  }
  if (!is.null(phi)) {
    J <- ncol(phi)
  }
  if (is.null(phi)) {
    phi <- matrix(rnorm((n + 1) * J), ncol = J)
  }
  if (is.null(mu)) {
    mu <- matrix(rnorm(J * (t0 + 1)), ncol = (t0 + 1))
  }
  if (is.null(alpha)) {
    alpha <- rep(0, (n + 1))
  }
  if (is.null(nu)) {
    nu <- rep(0, (t0 + 1))
  }
  if (is.null(X)) {
    X <- array(0, dim = c((n + 1), (t0 + 1), 3))
  }
  if (is.null(beta)) {
    beta <- matrix(rnorm(3), ncol = 1)
  }

  ## error term
  epsilon <- rnorm((n + 1) * (t0 + 1), sd = epsilon_sd)

  ## covariates
  Xbeta <- matrix(0, nrow = n + 1, ncol = t0 + 1)
  for (i in 1:(n + 1)) {
    for (j in 1:(t0 + 1)) {
      Xbeta[i, j] <- sum(X[i, j, ] * beta)
    }
  }

  ## fixed effects
  alpha.matrix <- replicate(t0 + 1, alpha)
  nu.matrix <- t(replicate(n + 1, nu))

  Y <- matrix(0, nrow = (n + 1), ncol = (t0 + 1))
  Y <- Xbeta + alpha.matrix + nu.matrix + phi %*% mu + epsilon

  df_pre_full <- t(Y)
  colnames(df_pre_full) <- c(paste0("X", 1:n), "Y")
  df_pre_full <- tibble::as_tibble(df_pre_full, .name_repair = "check_unique")

  df_prepost_full <- df_pre_full |>
    dplyr::mutate(D = 0)
  df_prepost_full[(t0 + 1), "D"] <- 1
  df_prepost_full[(t0 + 1), "Y"] <- df_prepost_full[(t0 + 1), "Y"] + tau

  params <- list(
    alpha = alpha,
    nu = nu,
    phi = phi,
    mu = mu,
    beta = beta,
    epsilon = epsilon
  )

  res <- list(
    df_pre_full = df_pre_full |> dplyr::slice(1:(t0)),
    df_prepost_full = df_prepost_full,
    tau = tau,
    params = params,
    .call = match.call()
  )
  class(res) <- "synth_ife"
  return(res)
}
