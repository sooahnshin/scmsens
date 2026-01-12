#' Generate Synthetic Data with Vertical Regression Model
#'
#' Generate synthetic panel data using a basic vertical regression model for
#' simulation studies. The data generation process creates correlated control units
#' with specified mean differences, which can be used to test sensitivity analysis
#' methods under controlled conditions.
#'
#' The model generates outcomes as: \eqn{Y = X \beta + \epsilon}, where \eqn{X} are
#' the control units drawn from a multivariate normal distribution.
#'
#' @param t0 Integer specifying the number of pre-treatment time periods.
#'   Default is 99.
#' @param n Integer specifying the number of control (donor) units. Default is 20.
#' @param tau Numeric value specifying the true Average Treatment Effect on the
#'   Treated (ATT). Default is 10.
#' @param rho Numeric parameter controlling the correlation among donor units.
#'   Larger values create stronger correlations. Default is 0.04.
#' @param mu Numeric vector of length \code{n} specifying the mean values for each
#'   donor unit. If \code{NULL} (default), randomly generated from N(0, 2).
#' @param beta Numeric vector of length \code{n} specifying the true weights (coefficients)
#'   for each donor unit. If \code{NULL} (default), randomly generated from N(0, 2).
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{df_pre_full}{A tibble of dimension (t0 x (n+1)) containing the pre-treatment
#'     data. Columns X1, X2, ..., Xn are the control units, and Y is the treated unit.}
#'   \item{df_prepost_full}{A tibble of dimension ((t0+1) x (n+2)) containing both
#'     pre- and post-treatment data. Includes a treatment indicator column D.
#'     The last row is the post-treatment observation.}
#'   \item{beta}{The numeric vector of weights used in data generation.}
#'   \item{tau}{The true ATT value.}
#'   \item{.call}{The matched call used to create this object.}
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
#' @importFrom tibble as_tibble
#'
#' @examples
#' # Generate synthetic data
#' set.seed(123)
#' synth <- generate_synth_data_vertreg(t0 = 99, n = 20, tau = 10, rho = 0.04)
#'
#' # Access the pre-treatment data
#' head(synth$df_pre_full)
#'
#' # The true effect is stored in the output
#' synth$tau
#'
#' @seealso \code{\link{generate_synth_data_ar}} for AR model-based data,
#'   \code{\link{generate_synth_data_ife}} for interactive fixed effects model
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


#' Generate Synthetic Data with Autoregressive Model
#'
#' Generate synthetic panel data where all units (treated and controls) follow an
#' autoregressive process of order l. This is useful for simulation studies where
#' temporal dynamics are important.
#'
#' Each unit is generated independently from an AR(l) process with the same
#' coefficients but potentially different initial values.
#'
#' @param t0 Integer specifying the number of pre-treatment time periods.
#'   Default is 99.
#' @param n Integer specifying the number of control (donor) units. Default is 20.
#' @param tau Numeric value specifying the true Average Treatment Effect on the
#'   Treated (ATT). This is added to the treated unit's outcome in the post-treatment
#'   period. Default is 10.
#' @param ar Numeric vector of AR coefficients. The length determines the order of
#'   the AR process. Default is \code{c(0.7, 0.2, -0.1, -0.3)} for an AR(4) process.
#' @param init Optional matrix of dimension (n+1) x length(ar) specifying initial
#'   values for each unit. If \code{NULL} (default), random initial values are used.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{df_pre_full}{A tibble of dimension (t0 x (n+1)) containing the pre-treatment
#'     data. Columns X1, X2, ..., Xn are the control units, and Y is the treated unit.}
#'   \item{df_prepost_full}{A tibble of dimension ((t0+1) x (n+2)) containing both
#'     pre- and post-treatment data. Includes a treatment indicator column D.
#'     The last row is the post-treatment observation.}
#'   \item{tau}{The true ATT value.}
#'   \item{.call}{The matched call used to create this object.}
#' }
#'
#' @examples
#' # Generate AR(4) synthetic data
#' set.seed(123)
#' synth <- generate_synth_data_ar(t0 = 99, n = 20, tau = 10)
#'
#' # Access the data
#' head(synth$df_pre_full)
#'
#' # Use with custom AR coefficients
#' synth_ar2 <- generate_synth_data_ar(t0 = 50, n = 10, tau = 5, ar = c(0.5, 0.3))
#'
#' @importFrom stats arima.sim
#' @importFrom dplyr mutate
#'
#' @export
#' @seealso \code{\link{simulate_ar_l_with_initial}} for simulating individual AR series,
#'   \code{\link{generate_synth_data_vertreg}}, \code{\link{generate_synth_data_ife}}
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

#' Generate Synthetic Data with Interactive Fixed Effects Model
#'
#' Generate synthetic panel data using an interactive fixed effects (IFE) model,
#' which is commonly used in synthetic control applications. The IFE model allows
#' for unit-specific factor loadings interacting with time-varying factors.
#'
#' The model generates outcomes as:
#' \deqn{Y_{it} = X_{it}'\beta + \alpha_i + \nu_t + \phi_i'\mu_t + \epsilon_{it}}
#'
#' where \eqn{\phi_i} are unit-specific factor loadings and \eqn{\mu_t} are
#' time-varying factors.
#'
#' @param t0 Integer specifying the number of pre-treatment time periods.
#'   Default is 99.
#' @param n Integer specifying the number of control (donor) units. Default is 20.
#' @param tau Numeric value specifying the true Average Treatment Effect on the
#'   Treated (ATT). Default is 10.
#' @param phi Optional matrix of dimension ((n+1) x J) containing unit-specific
#'   factor loadings. If \code{NULL}, generated randomly from N(0, 1).
#' @param mu Optional matrix of dimension (J x (t0+1)) containing time-varying
#'   factors. If \code{NULL}, generated randomly from N(0, 1).
#' @param alpha Optional vector of length (n+1) containing additive unit fixed effects.
#'   If \code{NULL} (default), set to zero.
#' @param nu Optional vector of length (t0+1) containing additive time fixed effects.
#'   If \code{NULL} (default), set to zero.
#' @param X Optional array of dimension ((n+1) x (t0+1) x k) containing time-varying
#'   covariates. If \code{NULL} (default), no covariates are included.
#' @param beta Optional vector of length k containing coefficients for time-varying
#'   covariates. If \code{NULL} (default), generated randomly from N(0, 1).
#' @param epsilon_sd Numeric value specifying the standard deviation of the
#'   error term. Default is 1.
#' @param J Integer specifying the number of factors (dimension of interactive
#'   fixed effects). Required if \code{phi} and \code{mu} are not provided.
#'
#' @return A list of class \code{"synth_ife"} containing:
#' \describe{
#'   \item{df_pre_full}{A tibble of dimension (t0 x (n+1)) containing the pre-treatment
#'     data. Columns X1, X2, ..., Xn are the control units, and Y is the treated unit.}
#'   \item{df_prepost_full}{A tibble of dimension ((t0+1) x (n+2)) containing both
#'     pre- and post-treatment data. Includes a treatment indicator column D.}
#'   \item{tau}{The true ATT value.}
#'   \item{params}{A list containing the parameters used for data generation:
#'     alpha, nu, phi, mu, beta, and epsilon.}
#'   \item{.call}{The matched call used to create this object.}
#' }
#'
#' @importFrom dplyr slice
#'
#' @examples
#' # Generate IFE synthetic data with 2 factors
#' set.seed(123)
#' synth <- generate_synth_data_ife(t0 = 99, n = 20, tau = 10, J = 2)
#'
#' # Access the generated data
#' head(synth$df_pre_full)
#'
#' # Access the parameters used for generation
#' dim(synth$params$phi)  # Factor loadings
#' dim(synth$params$mu)   # Time factors
#'
#' # Generate with custom factor loadings
#' custom_phi <- matrix(rnorm(42), ncol = 2)  # 21 units x 2 factors
#' custom_mu <- matrix(rnorm(200), nrow = 2)  # 2 factors x 100 time periods
#' synth_custom <- generate_synth_data_ife(
#'   t0 = 99, n = 20, tau = 5,
#'   phi = custom_phi, mu = custom_mu
#' )
#'
#' @seealso \code{\link{generate_synth_data_vertreg}}, \code{\link{generate_synth_data_ar}}
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
    stop("Either 'phi' (factor loadings matrix) or 'J' (number of factors) must be specified.",
         call. = FALSE)
  }
  if (is.null(mu) & is.null(J)) {
    stop("Either 'mu' (time factors matrix) or 'J' (number of factors) must be specified.",
         call. = FALSE)
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
