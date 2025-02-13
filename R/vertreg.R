#' Function to run a vertical regression.
#'
#' @param formula A string to specify the vertical regression model (e.g., "Y ~ X1 + X2" where "Y" is treated unit, "X1" and "X2" are control units).
#' @param data_pre Pre-period data where the rows are time points and the columns include treated/control units.
#' @param data_post Post-period data. Should have identical columns to data_pre.
#' @return A vector of treatment effects.
#'
#' @importFrom stats as.formula lm predict update
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
#' vertreg(
#'   formula = formula_full,
#'   data_pre = df_pre_full,
#'   data_post = df_post_full
#' )
#'
#' @export
vertreg <- function(
    formula, data_pre, data_post) {
  # get the response variable and remove intercept
  y_var_name <- all.vars(as.formula(formula))[1]
  fm <- as.formula(formula)
  fm <- update(fm, . ~ . - 1)

  # fit the vertical regression
  fit <- lm(fm, data = data_pre)

  # predict Y(0) for the post-period
  Y0 <- predict(fit, newdata = data_post)

  # get Y(1)
  Y1 <- data_post[, y_var_name]

  # treatment effect
  tau <- as.vector(Y1 - Y0)

  return(tau)
}
