#' Vertical Regression
#'
#' Function to run a vertical regression.
#'
#' @param formula A string to specify the vertical regression model (e.g., "Y ~ X1 + X2" where "Y" is treated unit, "X1" and "X2" are control units).
#' @param data_pre Pre-period data where the rows are time points and the columns include treated/control units. See `synth_pre` for example.
#' @param data_post Post-period data. Should have identical columns to data_pre. See `synth_post` for example.
#' @param pseudo_inverse A logical to indicate whether to use the pseudo-inverse to fit the model.
#' @return A vector of treatment effects.
#'
#' @seealso [vertreg_stacked()]
#'
#' @importFrom stats as.formula lm predict update
#' @importFrom tibble is_tibble
#'
#' @examples
#' vertreg(
#'   formula = paste("Y ~ -1 + ", paste(paste0("X", 1:16), collapse = " + ")),
#'   data_pre = synth_pre,
#'   data_post = synth_post
#' )
#'
#' @export
vertreg <- function(formula, data_pre, data_post, pseudo_inverse = FALSE) {
  if (is.matrix(data_pre) | is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) | is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
  }
  # get the response variable and remove intercept
  y_var_name <- all.vars(as.formula(formula))[1]
  fm <- as.formula(formula)
  fm <- update(fm, . ~ . - 1)

  # fit the vertical regression
  if (isTRUE(pseudo_inverse)) {
    fit <- lm_pseudo(fm, data_pre)
  } else {
    fit <- lm(fm, data = data_pre)
  }

  # predict Y(0) for the post-period
  Y0 <- predict(fit, newdata = data_post)

  # get Y(1)
  Y1 <- data_post[, y_var_name]

  # treatment effect
  tau <- as.vector(Y1 - Y0)

  return(tau)
}

#' Vertical Regression (With Treatment Indicator)
#'
#' Function to run a vertical regression in stacked form with a treatment indicator.
#'
#' @param formula A string to specify the vertical regression model (e.g., "Y ~ X1 + X2" where "Y" is treated unit, "X1" and "X2" are control units).
#' @param data_pre Pre-period data where the rows are time points and the columns include treated/control units. See `synth_pre` for example.
#' @param data_post Post-period data. Should have identical columns to data_pre. See `synth_post` for example.
#' @param pseudo_inverse A logical to indicate whether to use the pseudo-inverse to fit the model.
#' @return A data.frame with the treatment effect and its standard error based on the vertical regression:
#'  - `time_from_treatment`: the time from the treatment
#'  - `estimate`: the estimated treatment effect
#'  - `std.error`: the standard error of the estimated treatment effect
#'  - `statistic`: the t-statistic of the estimated treatment effect
#'  - `p.value`: the p-value of the estimated treatment effect
#'  - `df`: the degree of freedom of the regression
#'
#'  @seealso [vertreg()]
#'
#' @importFrom dplyr bind_rows filter mutate relocate select n
#' @importFrom purrr map
#' @importFrom broom tidy
#'
#' @examples
#' vertreg_stacked(
#'   formula = paste("Y ~ -1 + ", paste(paste0("X", 1:16), collapse = " + ")),
#'   data_pre = synth_pre,
#'   data_post = synth_post
#' )
#'
#' @export
vertreg_stacked <- function(formula, data_pre, data_post, pseudo_inverse = FALSE) {
  if (is.matrix(data_pre) | is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) | is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
  }

  # combine data and add post treatment time indicator
  data_prepost <- dplyr::bind_rows(
    dplyr::mutate(data_pre, D = 0, t = -rev((1:dplyr::n()) - 1)),
    dplyr::mutate(data_post, D = 1, t = 1:dplyr::n())
  )

  # remove intercept
  fm <- as.formula(formula) |>
    update(. ~ . - 1 + D)

  # fit the vertical regression for each post treatment period
  if (isTRUE(pseudo_inverse)) {
    # for computing the number of control units
    X <- model.matrix(fm, data_prepost)
    res <- purrr::map(1:nrow(data_post), ~ {
      fit <- lm_pseudo(fm, data = data_prepost |> dplyr::filter(t <= 0 | t == .x))
      tidy.lm_pseudo(fit) |>
        dplyr::filter(term == "D") |>
        dplyr::select(-term) |>
        dplyr::mutate(time_from_treatment = .x) |>
        dplyr::relocate(time_from_treatment) |>
        # degree of freedom of restricted regression: # of pretreatment periods + 1 - # of control units - 1
        dplyr::mutate(df = nrow(data_pre) - (ncol(X) - 1))
    }) |>
      dplyr::bind_rows()
  } else {
    res <- purrr::map(1:nrow(data_post), ~ {
      fit <- lm(fm, data = data_prepost |> dplyr::filter(t <= 0 | t == .x))
      broom::tidy(fit) |>
        dplyr::filter(term == "D") |>
        dplyr::select(-term) |>
        dplyr::mutate(time_from_treatment = .x) |>
        dplyr::relocate(time_from_treatment) |>
        dplyr::mutate(df = fit$df)
    }) |>
      dplyr::bind_rows()
  }

  return(res)
}
