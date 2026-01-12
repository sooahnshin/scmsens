#' Leave-one-out Sensitivity Analysis
#'
#' Perform sensitivity analysis by systematically dropping one control unit at a time
#' to estimate bias parameters. This function implements the leave-one-out approach
#' described in Liu, Shin, and Yamauchi (2024) for assessing the robustness of
#' synthetic control estimates to potential unobserved donor units.
#'
#' @param var_y_name A character string specifying the name of the treated unit column
#'   in the data. This should match exactly with a column name in both \code{data_pre}
#'   and \code{data_post}.
#' @param var_x_name A character vector specifying the names of the control unit columns.
#'   Each name should match exactly with column names in both \code{data_pre} and
#'   \code{data_post}. Typically created with \code{paste0("X", 1:n)} for n control units.
#' @param data_pre A data.frame, tibble, or matrix containing the pre-treatment period data.
#'   Rows represent time points and columns represent units (both treated and controls).
#'   Must contain all columns specified in \code{var_y_name} and \code{var_x_name}.
#'   See \code{\link{synth_pre}} for an example of the expected format.
#' @param data_post A data.frame, tibble, or matrix containing the post-treatment period data.
#'   Must have the same column structure as \code{data_pre}. Can contain one or more
#'   post-treatment time periods. See \code{\link{synth_post}} for an example.
#' @param pseudo_inverse Logical indicating whether to use the Moore-Penrose pseudo-inverse
#'   to fit the regression model. Set to \code{TRUE} when dealing with collinear or
#'   near-collinear control units. Default is \code{FALSE}.
#'
#' @return A data.frame (tibble) with the following columns:
#'    - `dropped_unit`: Name of the control unit that is dropped.
#'    - `time_from_treatment`: Time from the treatment.
#'    - `estimate`: Estimated treatment effect.
#'    - `std.error`: Standard error of the treatment effect (based on the vertical regression).
#'    - `statistic`: t-statistic of the treatment effect.
#'    - `p.value`: p-value of the treatment effect.
#'    - `gamma`: Weight of the control unit.
#'    - `delta`: Imbalance of the control unit.
#'    - `bias`: Bias of the treatment effect.
#'    - `tau`: Estimated treatment effect with the entire data.
#'    - `tau_plus_bias`: Sum of the estimated treatment effect and the bias. Must be close to the estimate (= estimate + noise).
#'    - `r2_Y_Z`: Partial R squared of the control unit with the outcome.
#'    - `r2_D_Z`: Partial R squared of the control unit with the treatment indicator.
#'    - `tau_se`: Standard error of the treatment effect with the entire data.
#'    - `t_stat`: t-Statistic of the treatment effect with the entire data.
#'
#'
#' @examples
#' estimate_params(
#'   var_y_name = "Y",
#'   var_x_name = paste0("X", 1:16),
#'   data_pre = synth_pre,
#'   data_post = synth_post,
#'   pseudo_inverse = FALSE
#' )
#'
#' @export
#' @importFrom dplyr pull
#' @importFrom sensemakr partial_r2
#'
#' @seealso [estimate_params_partial()], [estimate_params_partial_multi()]
#'
#' @references
#' Liu, N., Shin, S., & Yamauchi, S. (2024). Synthetic Control Method with
#' Missing Pre-treatment Outcomes. Working Paper.
#'
estimate_params <- function(var_y_name,
                            var_x_name,
                            data_pre,
                            data_post,
                            pseudo_inverse = FALSE) {
  # Input validation

  if (!is.character(var_y_name) || length(var_y_name) != 1) {
    stop("'var_y_name' must be a single character string specifying the treated unit column name.",
         call. = FALSE)
  }
  if (!is.character(var_x_name) || length(var_x_name) < 1) {
    stop("'var_x_name' must be a character vector with at least one control unit name.",
         call. = FALSE)
  }
  if (is.null(data_pre) || (nrow(data_pre) < 1)) {
    stop("'data_pre' must be a non-empty data.frame, tibble, or matrix.",
         call. = FALSE)
  }
  if (is.null(data_post) || (nrow(data_post) < 1)) {
    stop("'data_post' must be a non-empty data.frame, tibble, or matrix.",
         call. = FALSE)
  }
  
  if (is.matrix(data_pre) | is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) | is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
  }
  
  # Validate column names exist in data
  all_vars <- c(var_y_name, var_x_name)
  missing_pre <- setdiff(all_vars, names(data_pre))
  if (length(missing_pre) > 0) {
    stop("The following columns are missing from 'data_pre': ",
         paste(missing_pre, collapse = ", "), ".",
         call. = FALSE)
  }
  missing_post <- setdiff(all_vars, names(data_post))
  if (length(missing_post) > 0) {
    stop("The following columns are missing from 'data_post': ",
         paste(missing_post, collapse = ", "), ".",
         call. = FALSE)
  }

  data_prepost <- dplyr::bind_rows(
    dplyr::mutate(data_pre, D = 0, t = -rev((1:dplyr::n()) - 1)),
    dplyr::mutate(data_post, D = 1, t = 1:dplyr::n())
  )
  # fit VR with entire data to get weights
  formula_full <- as.formula(paste(var_y_name, "~ -1 + ", paste(var_x_name, collapse = " + ")))
  fit_full <- lm(formula_full, data = data_pre)
  weights_df <- broom::tidy(fit_full)

  # fit VR (stacked) with entire data to get the treatment effect and its standard error
  fit_complete <- vertreg_stacked(
    formula = formula_full,
    data_pre = data_pre,
    data_post = data_post,
    pseudo_inverse = pseudo_inverse
  )

  res <- purrr::map(var_x_name, function(z) {
    formula_drop <- paste(
      var_y_name,
      "~ -1 + ",
      paste(var_x_name |> setdiff(z), collapse = " + ")
    )
    # estimate treatment effect
    tau <- vertreg_stacked(
      formula = formula_drop,
      data_pre = data_pre,
      data_post = data_post,
      pseudo_inverse = pseudo_inverse
    ) |>
      dplyr::select(-df)
    # estimate gamma
    weight <- weights_df |>
      dplyr::filter(term == z) |>
      dplyr::pull(estimate)
    # estimate delta
    fm_z_on_x <- paste(z, "~ -1 +", paste(var_x_name |> setdiff(z), collapse = " + "))
    fit_zx <- lm(as.formula(fm_z_on_x), data = data_pre)
    pred_zT <- predict(fit_zx, newdata = data_post)
    obs_zT <- data_post |>
      dplyr::select(z) |>
      dplyr::pull()
    imbalance <- obs_zT - pred_zT
    # partial r squared
    fm_z_on_x_and_d <- as.formula(paste(z, "~ -1 +", paste(var_x_name |> setdiff(z), collapse = " + "), " + D"))
    r2_Y_Z <- sensemakr::partial_r2(fit_full, z)
    r2_D_Z <- purrr::map(1:nrow(data_post), function(x) {
      fit_zxd <- lm(fm_z_on_x_and_d, data = data_prepost |> dplyr::filter(t %in% c(-((1:nrow(data_pre)) - 1), x)))
      sensemakr::partial_r2(fit_zxd, "D")
    }) |>
      dplyr::bind_rows() |>
      dplyr::pull(D)
    tau |>
      dplyr::mutate(gamma = weight, delta = imbalance) |>
      dplyr::mutate(
        bias = gamma * delta,
        tau = fit_complete$estimate,
        tau_plus_bias = bias + tau # = estimate + noise
      ) |>
      dplyr::mutate(
        r2_Y_Z = r2_Y_Z, r2_D_Z = r2_D_Z,
        tau_se = std.error * sqrt(((1 - r2_Y_Z) / (1 - r2_D_Z)) * ((fit_complete$df + 1) / (fit_complete$df))),
        t_stat = tau / tau_se
      ) |>
      dplyr::mutate(dropped_unit = z) |>
      dplyr::relocate(dropped_unit)
  }) |>
    dplyr::bind_rows()

  return(res)
}

#' Estimate Sensitivity Parameters with Partially Observed Data (Single Unit)
#'
#' Estimate the bias parameters (gamma, delta, and bias) when a single control unit
#' has partially observed data (i.e., missing values in some pre-treatment periods).
#' This function uses the observed portion of the data to estimate how much bias
#' would result from excluding this control unit.
#'
#' For multiple control units with missing data, use \code{\link{estimate_params_partial_multi}}.
#'
#' @param fm_z_on_x A character string specifying the formula for regressing the
#'   partially observed control unit (Z) on the fully observed control units (X).
#'   For example, \code{"X1 ~ X2 + X3 + X4"} if X1 is the partially observed unit.
#'   The intercept is automatically removed.
#' @param fm_y_on_z_and_x A character string specifying the vertical regression formula
#'   with the treated unit as the outcome and all control units (including the
#'   partially observed one) as predictors. For example, \code{"Y ~ X1 + X2 + X3 + X4"}.
#'   The intercept is automatically removed.
#' @param data_pre A data.frame, tibble, or matrix containing the pre-treatment period data.
#'   Should only contain complete cases for the partially observed unit (i.e., rows
#'   where the partially observed unit is not NA). See \code{\link{synth_pre}} for format.
#' @param data_post A data.frame, tibble, or matrix containing the post-treatment period data.
#'   Currently only supports a single post-treatment period (one row).
#'   See \code{\link{synth_post}} for format.
#' @param pseudo_inverse Logical indicating whether to use the Moore-Penrose pseudo-inverse
#'   to fit the regression model. Default is \code{FALSE}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{gamma}{The weight (coefficient) of the partially observed control unit
#'     in the vertical regression.}
#'   \item{delta}{The imbalance (prediction error) for the partially observed unit
#'     in the post-treatment period.}
#'   \item{bias}{The estimated bias from excluding this control unit, calculated
#'     as \code{gamma * delta}.}
#' }
#'
#' @importFrom stats coef
#' @importFrom dplyr row_number
#'
#' @export
#'
#' @examples
#' # Suppose X1 is partially observed
#' synth_pre_partial <- synth_pre |>
#'   dplyr::mutate(X1 = ifelse(dplyr::row_number() < dplyr::n() / 2, NA, X1)) |>
#'   dplyr::filter(!is.na(X1))
#' estimate_params_partial(
#'   fm_z_on_x = paste(
#'     "X1 ~ -1 + ",
#'     paste(paste0("X", 2:16), collapse = " + ")
#'   ),
#'   fm_y_on_z_and_x = paste(
#'     "Y ~ -1 + ",
#'     paste(paste0("X", 1:16), collapse = " + ")
#'   ),
#'   data_pre = synth_pre_partial,
#'   data_post = synth_post,
#'   pseudo_inverse = FALSE
#' )
#'
#' @seealso [estimate_params()], [estimate_params_partial_multi()]
estimate_params_partial <- function(
    fm_z_on_x,
    fm_y_on_z_and_x,
    data_pre,
    data_post,
    pseudo_inverse = FALSE) {
  # Input validation
  if (!is.character(fm_z_on_x) || length(fm_z_on_x) != 1) {
    stop("'fm_z_on_x' must be a single character string specifying the formula.",
         call. = FALSE)
  }
  if (!is.character(fm_y_on_z_and_x) || length(fm_y_on_z_and_x) != 1) {
    stop("'fm_y_on_z_and_x' must be a single character string specifying the formula.",
         call. = FALSE)
  }
  if (is.null(data_pre) || (nrow(data_pre) < 1)) {
    stop("'data_pre' must be a non-empty data.frame, tibble, or matrix.",
         call. = FALSE)
  }
  if (is.null(data_post) || (nrow(data_post) < 1)) {
    stop("'data_post' must be a non-empty data.frame, tibble, or matrix.",
         call. = FALSE)
  }

  if (is.matrix(data_pre) | is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) | is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
  }

  # Check for single post-period with more informative message
  if (nrow(data_post) != 1) {
    stop("'data_post' must contain exactly one row (one post-treatment period). ",
         "You provided ", nrow(data_post), " rows. ",
         "For multiple post-treatment periods, run the analysis separately for each period.",
         call. = FALSE)
  }

  # variable with missing
  var_z_name <- all.vars(as.formula(fm_z_on_x))[1]

  # convert to formula
  fm_z_on_x <- as.formula(fm_z_on_x)
  # remove the intercept
  fm_z_on_x <- update(fm_z_on_x, . ~ . - 1)

  # estimate Z - ((X'X)^{-1}X'Z)'X using complete observations
  if (isTRUE(pseudo_inverse)) {
    fit_zx <- lm_pseudo(
      fm_z_on_x,
      data = data_pre
    )
  } else {
    fit_zx <- lm(
      fm_z_on_x,
      data = data_pre
    )
  }

  pred_zT <- predict(fit_zx, newdata = data_post)
  obs_zT <- data_post[1, var_z_name]
  delta <- obs_zT - pred_zT

  # convert to formula
  fm_y_on_z_and_x <- as.formula(fm_y_on_z_and_x)
  # remove the intercept
  fm_y_on_z_and_x <- update(fm_y_on_z_and_x, . ~ . - 1)

  # estimate \gamma coefficient using complete observations
  if (isTRUE(pseudo_inverse)) {
    fit_yzx <- lm_pseudo(
      fm_y_on_z_and_x,
      data = data_pre
    )
  } else {
    fit_yzx <- lm(
      fm_y_on_z_and_x,
      data = data_pre
    )
  }

  gamma <- coef(fit_yzx)[names(coef(fit_yzx)) == var_z_name]

  return(
    list(
      gamma = as.vector(gamma),
      delta = as.vector(delta),
      bias = as.vector(gamma * delta)
    )
  )
}

#' Estimate Sensitivity Parameters with Partially Observed Data (Multiple Units)
#'
#' Estimate the combined bias parameters when multiple control units have partially
#' observed data. This function creates a weighted combination of the partially
#' observed units and estimates the overall bias from excluding these units.
#'
#' For a single control unit with missing data, use \code{\link{estimate_params_partial}}.
#'
#' @param var_z_name A character vector specifying the names of the partially observed
#'   control unit columns. For example, \code{c("X1", "X2")} if both X1 and X2 have
#'   missing values.
#' @param var_x_name A character vector specifying the names of the fully observed
#'   control unit columns. These units should have complete data for all time periods.
#' @param var_y_name A character string specifying the name of the treated unit column.
#' @param data_pre A data.frame, tibble, or matrix containing the pre-treatment period data.
#'   Should contain complete cases where all partially observed units have values.
#'   See \code{\link{synth_pre}} for format.
#' @param data_post A data.frame, tibble, or matrix containing the post-treatment period data.
#'   Currently only supports a single post-treatment period (one row).
#'   See \code{\link{synth_post}} for format.
#' @param pseudo_inverse Logical indicating whether to use the Moore-Penrose pseudo-inverse
#'   to fit the regression model. Default is \code{FALSE}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{gamma}{Always returns 1, as the combined effect is captured in the imbalance.}
#'   \item{imbalance}{The combined imbalance (prediction error) for all partially
#'     observed units, weighted by their coefficients.}
#'   \item{bias}{The estimated bias from excluding these control units.}
#' }
#'
#' @importFrom tidyselect all_of
#'
#' @export
#'
#' @examples
#' # Suppose X1 and X2 are partially observed
#' synth_pre_partial <- synth_pre |>
#'   dplyr::mutate(
#'     X1 = ifelse(dplyr::row_number() < dplyr::n() / 2, NA, X1),
#'     X2 = ifelse(dplyr::row_number() < dplyr::n() / 2, NA, X2)
#'   ) |>
#'   dplyr::filter(!is.na(X1) | is.na(X1))
#' estimate_params_partial_multi(
#'   var_z_name = c("X1", "X2"),
#'   var_x_name = paste0("X", 3:16),
#'   var_y_name = "Y",
#'   data_pre = synth_pre_partial,
#'   data_post = synth_post,
#'   pseudo_inverse = FALSE
#' )
#'
#' @seealso [estimate_params()], [estimate_params_partial()]
estimate_params_partial_multi <- function(
    var_z_name,
    var_x_name,
    var_y_name,
    data_pre,
    data_post,
    pseudo_inverse = FALSE) {
  # Input validation
  if (!is.character(var_z_name) || length(var_z_name) < 1) {
    stop("'var_z_name' must be a character vector with at least one partially observed unit name.",
         call. = FALSE)
  }
  if (!is.character(var_x_name) || length(var_x_name) < 1) {
    stop("'var_x_name' must be a character vector with at least one fully observed unit name.",
         call. = FALSE)
  }
  if (!is.character(var_y_name) || length(var_y_name) != 1) {
    stop("'var_y_name' must be a single character string specifying the treated unit column name.",
         call. = FALSE)
  }
  if (is.null(data_pre) || (nrow(data_pre) < 1)) {
    stop("'data_pre' must be a non-empty data.frame, tibble, or matrix.",
         call. = FALSE)
  }
  if (is.null(data_post) || (nrow(data_post) < 1)) {
    stop("'data_post' must be a non-empty data.frame, tibble, or matrix.",
         call. = FALSE)
  }
  
  if (is.matrix(data_pre) | is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) | is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
  }
  
  # Check for single post-period with more informative message
  if (nrow(data_post) != 1) {
    stop("'data_post' must contain exactly one row (one post-treatment period). ",
         "You provided ", nrow(data_post), " rows. ",
         "For multiple post-treatment periods, run the analysis separately for each period.",
         call. = FALSE)
  }
  
  # Validate column names exist in data
  all_vars <- c(var_y_name, var_z_name, var_x_name)
  missing_pre <- setdiff(all_vars, names(data_pre))
  if (length(missing_pre) > 0) {
    stop("The following columns are missing from 'data_pre': ",
         paste(missing_pre, collapse = ", "), ".",
         call. = FALSE)
  }
  missing_post <- setdiff(all_vars, names(data_post))
  if (length(missing_post) > 0) {
    stop("The following columns are missing from 'data_post': ",
         paste(missing_post, collapse = ", "), ".",
         call. = FALSE)
  }

  fm_y_on_z_and_x <- paste(var_y_name, "~ -1 +", paste(var_z_name, collapse = "+"), "+", paste(var_x_name, collapse = "+"))
  fm_y_on_z_and_x <- as.formula(fm_y_on_z_and_x)
  # estimate \gamma coefficient using complete observations
  if (isTRUE(pseudo_inverse)) {
    fit_yzx <- lm_pseudo(
      fm_y_on_z_and_x,
      data = data_pre
    )
  } else {
    fit_yzx <- lm(
      fm_y_on_z_and_x,
      data = data_pre
    )
  }
  gamma <- coef(fit_yzx)[names(coef(fit_yzx)) %in% var_z_name]

  # make new single Z variable that is a linear combination of missing variables with weights being gamma
  data_pre_aug <- data_pre
  data_pre_aug$comb_Z <- as.matrix(select(data_pre, all_of(var_z_name))) %*% gamma
  data_post_aug <- data_post
  data_post_aug$comb_Z <- as.matrix(select(data_post, all_of(var_z_name))) %*% gamma

  # estimate Z - ((X'X)^{-1}X'Z)'X using complete observations
  fm_z_on_x <- paste("comb_Z ~ -1 + ", paste(var_x_name, collapse = "+"))
  fm_z_on_x <- as.formula(fm_z_on_x)

  if (isTRUE(pseudo_inverse)) {
    fit_zx <- lm_pseudo(
      fm_z_on_x,
      data = data_pre_aug
    )
  } else {
    fit_zx <- lm(
      fm_z_on_x,
      data = data_pre_aug
    )
  }

  pred_zT <- predict(fit_zx, newdata = data_post_aug)
  obs_zT <- data_post_aug[1, "comb_Z"]
  imbalance <- obs_zT - pred_zT

  return(
    list(
      gamma = 1,
      imbalance = as.numeric(imbalance),
      bias = as.numeric(imbalance)
    )
  )
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "D", "bias", "delta", "df", "dropped_unit", "estimate", "std.error", "tau_se", "term"
  ))
}
