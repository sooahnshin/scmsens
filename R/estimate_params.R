#' Leave-one-out Analysis
#'
#' Function to estimate bias parameters with dropping one control unit at a time
#'
#' @param var_y_name Name of the treated unit.
#' @param var_x_name Names of the control units.
#' @param data_pre Pre-period data where the rows are time points and the columns include treated/control units.
#' @param data_post Post-period data. Should have identical columns to data_pre.
#' @param pseudo_inverse A logical to indicate whether to use the pseudo-inverse to fit the model.
#' @return A data.frame with the following columns:
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
#' df_pre_full <- synth_data |>
#'   dplyr::filter(Dt == 0) |>
#'   dplyr::select(-year, -Dt)
#' df_post_full <- synth_data |>
#'   dplyr::filter(Dt == 1) |>
#'   dplyr::select(-year, -Dt)
#' formula_full <- as.formula(
#'   paste("Y ~ -1 + ", paste(paste0("X", 1:16), collapse = " + "))
#' )
#' estimate_params(var_y_name = "Y",
#'                var_x_name = paste0("X", 1:16),
#'                data_pre = df_pre_full,
#'                data_post = df_post_full,
#'                pseudo_inverse = FALSE)
#'
#' @export
#' @importFrom dplyr pull
#' @importFrom sensemakr partial_r2
estimate_params <- function(var_y_name,
                            var_x_name,
                            data_pre,
                            data_post,
                            pseudo_inverse = FALSE) {
  if (is.matrix(data_pre) | is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) | is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
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

#' Estimate Sensitivity Parameters
#'
#' Function to use partially observed data to estimate gamma and imbalance.
#'
#' @param fm_z_on_x A string specifying the regression of Z on X in the formula
#'   format (e.g., "Z ~ X1 + X2").
#' @param fm_y_on_z_and_x A string specifying the vertical regression with
#'   missing obs in the formula format (e.g, "Y ~ X1 + X2 + Z").
#' @param data_pre Pre-period data where the rows are time points and the columns include treated/control units.
#' @param data_post Post-period data. Should have identical columns to data_pre.
#' @param pseudo_inverse A logical to indicate whether to use the pseudo-inverse to fit the model.
#' @return A list of `gamma`, `imbalance`, and `bias`.
#'
#' @importFrom stats coef
#'
#' @export
#'
#' @examples
#' df_pre_full <- synth_data |>
#'   dplyr::filter(Dt == 0) |>
#'   dplyr::select(-year, -Dt)
#' df_post_full <- synth_data |>
#'   dplyr::filter(Dt == 1) |>
#'   dplyr::select(-year, -Dt)
#' # Suppose X1 is partially observed
#' df_pre_full <- df_pre_full |>
#'   dplyr::mutate(X1 = ifelse(rep(c(FALSE, TRUE), length.out = dplyr::n()), NA, X1))
#' df_pre_partial <- df_pre_full |>
#'  dplyr::filter(!is.na(X1))
#' fm_z_on_x <- as.formula(
#'   paste("X1 ~ -1 + ", paste(paste0("X", 2:16), collapse = " + "))
#' )
#' fm_y_on_z_and_x <- as.formula(
#'  paste("Y ~ -1 + ", paste(paste0("X", 1:16), collapse = " + "))
#'  )
#' estimate_params_partial(fm_z_on_x = fm_z_on_x,
#'                        fm_y_on_z_and_x = fm_y_on_z_and_x,
#'                        data_pre = df_pre_partial,
#'                        data_post = df_post_full,
#'                        pseudo_inverse = FALSE)
estimate_params_partial <- function(
    fm_z_on_x,
    fm_y_on_z_and_x,
    data_pre,
    data_post,
    pseudo_inverse = FALSE
) {

  if (is.matrix(data_pre) | is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) | is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
  }

  # assume T = T0 + 1 (only one post-period)
  if (nrow(data_post) != 1) stop("Only one post_period is allowed")

  # variable with missing
  var_z_name <- all.vars(as.formula(fm_z_on_x))[1]

  # convert to formula
  fm_z_on_x <- as.formula(fm_z_on_x)
  # remove the intercept
  fm_z_on_x <- update(fm_z_on_x, . ~ . -1)

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
  fm_y_on_z_and_x <- update(fm_y_on_z_and_x, . ~ . -1)

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

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "D", "bias", "delta", "df", "dropped_unit", "estimate", "std.error", "tau_se", "term"
  ))
}
