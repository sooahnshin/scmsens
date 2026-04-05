#' Vertical Regression for Synthetic Control
#'
#' Run a vertical regression to estimate treatment effects in the Synthetic Control
#' Method framework. This function fits a regression of the treated unit's outcomes
#' on the control units' outcomes in the pre-treatment period, then uses the fitted
#' model to predict the counterfactual outcome in the post-treatment period.
#'
#' The treatment effect is calculated as the difference between the observed
#' outcome and the predicted counterfactual outcome.
#'
#' @param formula A character string or formula object specifying the vertical
#'   regression model. The left-hand side should be the treated unit and the
#'   right-hand side should include all control units. For example,
#'   \code{"Y ~ X1 + X2 + X3"} or \code{"Y ~ -1 + X1 + X2 + X3"} (without intercept).
#'   The intercept is automatically removed regardless of the formula specification.
#' @param data_pre A data.frame, tibble, or matrix containing the pre-treatment period data.
#'   Rows represent time points and columns represent units. Must contain all
#'   variables specified in the formula. See \code{\link{synth_pre}} for an example.
#' @param data_post A data.frame, tibble, or matrix containing the post-treatment period data.
#'   Must have the same column structure as \code{data_pre}. Can contain one or more
#'   post-treatment time periods. See \code{\link{synth_post}} for an example.
#' @param pseudo_inverse Logical indicating whether to use the Moore-Penrose pseudo-inverse
#'   to fit the regression model. Set to \code{TRUE} when dealing with collinear or
#'   near-collinear control units. Default is \code{FALSE}.
#'
#' @return A numeric vector of estimated treatment effects, one for each post-treatment
#'   time period (row in \code{data_post}).
#'
#' @seealso \code{\link{vertreg_stacked}} for inference with standard errors,
#'   \code{\link{estimate_params}} for sensitivity analysis
#'
#' @importFrom stats as.formula lm predict update
#' @importFrom tibble is_tibble
#'
#' @examples
#' # Basic usage with synthetic data
#' vertreg(
#'   formula = paste("Y ~ -1 + ", paste(paste0("X", 1:16), collapse = " + ")),
#'   data_pre = synth_pre,
#'   data_post = synth_post
#' )
#'
#' # With pseudo-inverse for potentially collinear data
#' vertreg(
#'   formula = "Y ~ X1 + X2 + X3",
#'   data_pre = synth_pre,
#'   data_post = synth_post,
#'   pseudo_inverse = TRUE
#' )
#'
#' @export
vertreg <- function(formula, data_pre, data_post, pseudo_inverse = FALSE) {
  # Input validation
  if (is.null(formula) || (!is.character(formula) && !inherits(formula, "formula"))) {
    stop("'formula' must be a character string or formula object.",
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
#' Run a vertical regression in stacked form with a treatment indicator D,
#' as described in Liu, Shin, and Yamauchi (2024). Algebraically equivalent
#' to \code{\link{vertreg}} but connects to the omitted variable bias framework
#' and enables inference under additional assumptions.
#'
#' @param formula A character string or formula object specifying the vertical
#'   regression model. The left-hand side should be the treated unit and the
#'   right-hand side should include all control units. For example,
#'   \code{"Y ~ X1 + X2 + X3"}. The treatment indicator (D) is automatically added.
#' @param data_pre A data.frame, tibble, or matrix containing the pre-treatment period data.
#'   Rows represent time points and columns represent units. Must contain all
#'   variables specified in the formula. See \code{\link{synth_pre}} for an example.
#' @param data_post A data.frame, tibble, or matrix containing the post-treatment period data.
#'   Must have the same column structure as \code{data_pre}. Can contain one or more
#'   post-treatment time periods. See \code{\link{synth_post}} for an example.
#' @param pseudo_inverse Logical indicating whether to use the Moore-Penrose pseudo-inverse
#'   to fit the regression model. Default is \code{FALSE}.
#'
#' @return A data.frame (tibble) with the following columns:
#' \describe{
#'   \item{time_from_treatment}{Integer indicating the post-treatment period (1, 2, 3, ...)}
#'   \item{estimate}{The estimated treatment effect for that period}
#'   \item{std.error}{Standard error of the treatment effect estimate}
#'   \item{statistic}{The t-statistic for testing if the effect differs from zero}
#'   \item{p.value}{Two-sided p-value for the t-test}
#'   \item{df}{Degrees of freedom used in the regression}
#' }
#'
#' @seealso \code{\link{vertreg}} for estimates without standard errors,
#'   \code{\link{estimate_params}} for sensitivity analysis
#'
#' @importFrom dplyr bind_rows filter mutate relocate select n
#' @importFrom purrr map
#' @importFrom broom tidy
#'
#' @examples
#' # Basic usage
#' vertreg_stacked(
#'   formula = paste("Y ~ -1 + ", paste(paste0("X", 1:16), collapse = " + ")),
#'   data_pre = synth_pre,
#'   data_post = synth_post
#' )
#'
#' # Access specific results
#' result <- vertreg_stacked(
#'   formula = "Y ~ X1 + X2 + X3 + X4",
#'   data_pre = synth_pre,
#'   data_post = synth_post
#' )
#' result$estimate    # Treatment effect
#' result$std.error   # Standard error
#' result$p.value     # P-value
#'
#' @export
vertreg_stacked <- function(formula, data_pre, data_post, pseudo_inverse = FALSE) {
  # Input validation
  if (is.null(formula) || (!is.character(formula) && !inherits(formula, "formula"))) {
    stop("'formula' must be a character string or formula object.",
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
