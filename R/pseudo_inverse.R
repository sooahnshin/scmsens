#' Fitting Linear Models Using Pseudo-inverse
#'
#' Function to fit linear models using the pseudo-inverse.
#'
#' @param formula A string to specify the linear model.
#' @param data A data.frame containing the data.
#' @return A list containing the estimated coefficients.
#'
#' @importFrom stats model.matrix model.response model.frame
#' @importFrom MASS ginv
#' @importFrom tibble tibble
#'
#' @export
#'
#' @seealso [summary.lm_pseudo()], [predict.lm_pseudo()], [coef.lm_pseudo()], [tidy.lm_pseudo()]
lm_pseudo <- function(formula, data) {
  X <- model.matrix(formula, data)
  y <- model.response(model.frame(formula, data))
  # Compute pseudoinverse of X
  X_pinv <- MASS::ginv(X)

  # Compute coefficients using pseudoinverse
  beta_hat <- X_pinv %*% y

  res <- list(
    summary = tibble(term = colnames(X), estimate = as.numeric(beta_hat)),
    formula = formula,
    data = data
  )

  class(res) <- "lm_pseudo"

  return(res)
}

#' Summary Method for lm_pseudo Objects
#'
#' Extracts the summary of an `lm_pseudo` object, returning a tibble with estimated coefficients.
#'
#' @param object An object of class `lm_pseudo`.
#' @param ... Additional arguments (ignored).
#' @return A tibble containing the estimated coefficients with their corresponding terms.
#' @export
#'
#' @seealso [lm_pseudo()], [predict.lm_pseudo()], [coef.lm_pseudo()], [tidy.lm_pseudo()]
summary.lm_pseudo <- function(object, ...) {
  return(object$summary)
}

#' Predict Method for lm_pseudo Objects
#'
#' Generates predictions from an `lm_pseudo` object using new data.
#'
#' @param object An object of class `lm_pseudo`.
#' @param newdata A data frame containing new predictor values.
#' @param ... Additional arguments (ignored, for consistency with generic predict()).
#' @return A numeric vector of predicted values.
#' @importFrom stats model.matrix
#' @export
#'
#' @seealso [lm_pseudo()], [summary.lm_pseudo()], [coef.lm_pseudo()], [tidy.lm_pseudo()]
predict.lm_pseudo <- function(object, newdata, ...) {
  X <- model.matrix(object$formula, newdata)
  y_hat <- X %*% object$summary$estimate
  return(y_hat)
}

#' Extract Coefficients from lm_pseudo Objects
#'
#' Retrieves the estimated coefficients from an `lm_pseudo` object.
#'
#' @param object An object of class `lm_pseudo`.
#' @param ... Additional arguments (ignored, for consistency with generic `coef()`).
#' @return A named numeric vector of estimated coefficients.
#' @export
coef.lm_pseudo <- function(object, ...) {
  beta <- object$summary$estimate
  names(beta) <- object$summary$term
  return(beta)
}


#' Tidy Method for lm_pseudo Objects
#'
#' @param x An object of class `lm_pseudo`.
#' @param ... Additional arguments (ignored).
#' @return A tibble with columns: term, estimate, std.error, statistic, p.value.
#' @importFrom stats model.matrix model.response model.frame pt
#' @importFrom tibble tibble
#' @export
#'
#' @seealso [lm_pseudo()], [summary.lm_pseudo()], [predict.lm_pseudo()], [coef.lm_pseudo()]
tidy.lm_pseudo <- function(x, ...) {
  X <- model.matrix(x$formula, x$data)
  y <- model.response(model.frame(x$formula, x$data))

  # Compute residuals
  y_hat <- X %*% x$summary$estimate
  residuals <- y - y_hat
  n <- length(y)
  p <- ncol(X)

  # Estimate residual standard error
  sigma_sq <- sum(residuals^2) / (n - p)

  # Compute variance-covariance matrix using SVD (more stable)
  svd_X <- svd(X)
  X_pinv_svd <- svd_X$v %*% diag(1 / svd_X$d, length(svd_X$d)) %*% t(svd_X$u)  # Stable pseudo-inverse
  cov_matrix <- sigma_sq * (X_pinv_svd %*% t(X_pinv_svd))

  # Extract standard errors
  std_error <- sqrt(diag(cov_matrix))

  # Compute t-statistics and p-values
  t_stat <- x$summary$estimate / std_error
  p_values <- 2 * pt(-abs(t_stat), df = n - p)

  # Return a tibble
  tibble(
    term = x$summary$term,
    estimate = x$summary$estimate,
    std.error = std_error,
    statistic = t_stat,
    p.value = p_values
  )
}
