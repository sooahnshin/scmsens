#' Sensitivity Contour Plot (Weight and Imbalance)
#'
#' Generate a contour plot to visualize how the treatment effect estimate changes
#' based on the weight (gamma) and imbalance (delta) parameters of potential
#' unobserved donor units. This is a key visualization for sensitivity analysis
#' in the Synthetic Control Method framework.
#'
#' The contour lines show combinations of weight and imbalance that would result
#' in the same adjusted treatment effect. Points representing dropped control units
#' (from leave-one-out analysis) are plotted to show their contribution to potential bias.
#'
#' @param tau_complete Numeric value of the estimated treatment effect using
#'   all control units (the unadjusted estimate).
#' @param var_x_gamma Numeric vector of weights (gamma) for each control unit or
#'   benchmark point. These typically come from leave-one-out analysis via
#'   \code{\link{estimate_params}}.
#' @param var_x_delta Numeric vector of imbalances (delta) for each control unit or
#'   benchmark point. Same length as \code{var_x_gamma}.
#' @param var_x_name Character vector of labels for each point. Same length as
#'   \code{var_x_gamma}. Often includes the effect estimate, e.g.,
#'   \code{"Drop X1 (-2.3)"}.
#' @param title Character string for the plot title. Default is \code{NULL} (no title).
#' @param gamma_seq Numeric vector specifying the range of gamma values for the
#'   contour lines. If \code{NULL} (default), automatically determined from data.
#' @param delta_seq Numeric vector specifying the range of delta values for the
#'   contour lines. If \code{NULL} (default), automatically determined from data.
#' @param text_size Numeric value for the size of text labels. Default is 3.
#' @param nudge_x Numeric value for horizontal offset of text labels. Default is 0.02.
#' @param nudge_y Numeric value for vertical offset of text labels. Default is 0.
#' @param repel Logical indicating whether to use \code{ggrepel::geom_text_repel()}
#'   for automatic label positioning to avoid overlaps. Default is \code{FALSE}.
#' @param tau_round Integer specifying the number of decimal places for rounding
#'   the treatment effect in the legend. Default is 1.
#' @param nudge_x2 Numeric value for horizontal offset of the unadjusted estimate label.
#'   If \code{NULL} (default), uses \code{nudge_x}.
#' @param nudge_y2 Numeric value for vertical offset of the unadjusted estimate label.
#'   If \code{NULL} (default), uses \code{nudge_y}.
#'
#' @return A ggplot2 object that can be further customized or printed.
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#'
#' @examples
#' # Perform leave-one-out sensitivity analysis
#' res <- estimate_params(var_y_name = "Y",
#'                        var_x_name = paste0("X", 1:16),
#'                        data_pre = synth_pre,
#'                        data_post = synth_post)
#'
#' # Extract parameters for plotting
#' tau_complete <- res |> dplyr::pull(tau) |> unique()
#' var_x_gamma <- res |> dplyr::pull(gamma)
#' var_x_delta <- res |> dplyr::pull(delta)
#' var_x_tau <- res |> dplyr::pull(estimate)
#' var_x_name <- paste0("Drop ", paste0("X", 1:16), " (", round(var_x_tau, 1), ")")
#'
#' # Create sensitivity plot
#' plot_sensitivity(tau_complete = tau_complete,
#'                  var_x_gamma = var_x_gamma,
#'                  var_x_delta = var_x_delta,
#'                  var_x_name = var_x_name,
#'                  title = "Sensitivity of the SCM Estimate",
#'                  text_size = 4)
#'
#' # Use repel for better label placement with many points
#' plot_sensitivity(tau_complete = tau_complete,
#'                  var_x_gamma = var_x_gamma,
#'                  var_x_delta = var_x_delta,
#'                  var_x_name = var_x_name,
#'                  repel = TRUE)
#'
#' @seealso [plot_sensitivity_r2()] for R-squared based sensitivity plots,
#'   [estimate_params()] for generating the input data
plot_sensitivity <- function(
    tau_complete, var_x_gamma, var_x_delta, var_x_name, title = NULL,
    gamma_seq = NULL, delta_seq = NULL, text_size = 3, nudge_x = 0.02, nudge_y = 0,
    repel = FALSE, tau_round = 1, nudge_x2 = NULL, nudge_y2 = NULL) {
  # Input validation
  if (!is.numeric(tau_complete) || length(tau_complete) != 1) {
    stop("'tau_complete' must be a single numeric value.", call. = FALSE)
  }
  if (!is.numeric(var_x_gamma) || length(var_x_gamma) < 1) {
    stop("'var_x_gamma' must be a numeric vector with at least one element.", call. = FALSE)
  }
  if (!is.numeric(var_x_delta) || length(var_x_delta) < 1) {
    stop("'var_x_delta' must be a numeric vector with at least one element.", call. = FALSE)
  }
  if (length(var_x_gamma) != length(var_x_delta)) {
    stop("'var_x_gamma' and 'var_x_delta' must have the same length.", call. = FALSE)
  }
  if (length(var_x_gamma) != length(var_x_name)) {
    stop("'var_x_name' must have the same length as 'var_x_gamma' and 'var_x_delta'.", call. = FALSE)
  }
  
  if (is.null(nudge_x2)) nudge_x2 <- nudge_x
  if (is.null(nudge_y2)) nudge_y2 <- nudge_y
  if (is.null(gamma_seq)) {
    gamma_seq <- seq(
      round(min(var_x_gamma), 0) - 0.5,
      round(max(var_x_gamma), 0) + 0.5,
      length.out = 100
    )
  }
  if (is.null(delta_seq)) {
    delta_seq <- seq(
      round(min(var_x_delta), 0) - 0.5,
      round(max(var_x_delta), 0) + 0.5,
      length.out = 100
    )
  }
  out <- tibble(
    var_name = var_x_name,
    gamma = var_x_gamma,
    delta = var_x_delta
  )
  tau <- tibble(
    label = paste0("Unadjusted (", round(tau_complete, tau_round), ")"),
    gamma = 0, delta = 0
  )

  sens_df <- expand.grid(gamma_seq, delta_seq) |>
    as.data.frame() |>
    dplyr::mutate(z = tau_complete + Var1 * Var2)

  if (!repel) {
    p <- sens_df |>
      ggplot() +
      geom_vline(xintercept = 0, color = "black", alpha = 0.1) +
      geom_hline(yintercept = 0, color = "black", alpha = 0.1) +
      geom_contour(aes(x = Var1, y = Var2, z = z), color = "black", alpha = 0.5) +
      # metR::geom_text_contour(aes(x = Var1, y = Var2, z = z), size = 3, skip = 0, alpha = 0.5, nudge_x = 0.05) +
      geom_point(data = out, aes(x = gamma, y = delta)) +
      geom_text(
        data = out, aes(x = gamma, y = delta, label = var_name),
        hjust = 0, nudge_x = nudge_x, nudge_y = nudge_y, size = text_size
      ) +
      geom_point(data = tau, aes(x = gamma, y = delta), shape = 17, color = "#990000") +
      geom_text(
        data = tau, aes(x = gamma, y = delta, label = label), color = "#990000",
        hjust = 0, nudge_x = nudge_x2, nudge_y = nudge_y2, size = text_size
      ) +
      labs(
        x = expression(hat(gamma) ~ "(Weight)"),
        y = expression(hat(delta) ~ "(Imbalance)"),
        title = title
      ) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  } else {
    p <- sens_df |>
      ggplot() +
      geom_vline(xintercept = 0, color = "black", alpha = 0.1) +
      geom_hline(yintercept = 0, color = "black", alpha = 0.1) +
      geom_contour(aes(x = Var1, y = Var2, z = z), color = "black", alpha = 0.5) +
      # metR::geom_text_contour(aes(x = Var1, y = Var2, z = z), size = 3, skip = 0, alpha = 0.5, nudge_x = 0.05) +
      geom_point(data = out, aes(x = gamma, y = delta)) +
      ggrepel::geom_text_repel(
        data = out, aes(x = gamma, y = delta, label = var_name),
        hjust = 0, nudge_x = nudge_x, nudge_y = nudge_y, size = text_size
      ) +
      geom_point(data = tau, aes(x = gamma, y = delta), shape = 17, color = "#990000") +
      geom_text(
        data = tau, aes(x = gamma, y = delta, label = label), color = "#990000",
        hjust = 0, nudge_x = nudge_x2, nudge_y = nudge_y2, size = text_size
      ) +
      labs(
        x = expression(hat(gamma) ~ "(Weight)"),
        y = expression(hat(delta) ~ "(Imbalance)"),
        title = title
      ) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  }
  return(p)
}

#' Sensitivity Contour Plot (Partial R-squared)
#'
#' Generate a sensitivity contour plot using partial R-squared values to assess
#' how unobserved confounders with certain explanatory power could affect the
#' treatment effect estimate. This approach follows the framework of Cinelli and
#' Hazlett (2020).
#'
#' The x-axis shows the partial R-squared of the outcome regressed on the
#' confounder (controlling for treatment and observed controls), while the y-axis
#' shows the partial R-squared of the treatment indicator regressed on the
#' confounder (controlling for observed controls).
#'
#' @param tau_complete Numeric value of the estimated treatment effect using
#'   all control units.
#' @param df Integer specifying the degrees of freedom from the vertical regression.
#'   Can be obtained from the output of \code{\link{vertreg_stacked}}.
#' @param tau_se Numeric value of the standard error of the treatment effect.
#'   Can be obtained from the output of \code{\link{vertreg_stacked}}.
#' @param sign_bias Integer value (1 or -1) indicating the sign of the bias to consider.
#'   Use 1 for positive bias scenarios and -1 for negative bias scenarios.
#' @param var_x_r2_Y_Z Numeric vector of partial R-squared values for the outcome
#'   regression. Typically from \code{\link{estimate_params}}.
#' @param var_x_r2_D_Z Numeric vector of partial R-squared values for the treatment
#'   indicator regression. Same length as \code{var_x_r2_Y_Z}.
#' @param var_x_name Character vector of labels for each benchmark point.
#'   Same length as \code{var_x_r2_Y_Z}.
#' @param title Character string for the plot title. Default is \code{NULL} (no title).
#' @param r2_Y_Z_seq Numeric vector specifying the range of R-squared values for the
#'   x-axis. Default is \code{seq(0, 1, length.out = 100)}.
#' @param r2_D_Z_seq Numeric vector specifying the range of R-squared values for the
#'   y-axis. Default is \code{seq(0, 1, length.out = 100)}.
#' @param text_size Numeric value for the size of text labels. Default is 3.
#' @param plot_t_stat Logical indicating whether to plot t-statistics on contours
#'   instead of treatment effects. If \code{TRUE}, the \code{t_stat} argument should
#'   also be provided. Default is \code{FALSE}.
#' @param t_stat Numeric value of the t-statistic for the treatment effect (required
#'   only if \code{plot_t_stat = TRUE}).
#' @param critical_value Numeric value specifying a critical t-value for highlighting
#'   a significance threshold contour (e.g., 1.96 for 95% confidence).
#'   If \code{NULL} (default), no critical value line is drawn.
#'
#' @return A ggplot2 object that can be further customized or printed.
#'
#' @references
#' Cinelli, C., & Hazlett, C. (
#'   2020). Making sense of sensitivity: Extending omitted variable bias.
#' \emph{Journal of the Royal Statistical Society: Series B}, 82(1), 39-67.
#'
#' @import ggplot2
#' @importFrom metR geom_text_contour
#' @importFrom grDevices contourLines
#'
#' @export
#' @seealso [plot_sensitivity()] for weight/imbalance plots,
#'   [estimate_params()] for generating the input data
#'
#' @examples
#' # First, fit the model to get df and standard error
#' fit <- vertreg_stacked(formula = paste("Y ~ -1 + ",
#'                                        paste(paste0("X", 1:16),
#'                                              collapse = " + ")),
#'                        data_pre = synth_pre,
#'                        data_post = synth_post,
#'                        pseudo_inverse = FALSE)
#' df <- fit$df
#' tau_se <- fit$std.error
#'
#' # Run leave-one-out analysis
#' res <- estimate_params(var_y_name = "Y",
#'                        var_x_name = paste0("X", 1:16),
#'                        data_pre = synth_pre,
#'                        data_post = synth_post)
#'
#' # Plot only top 5 units with positive bias
#' res <- res[order(-res$bias), ][1:5, ]
#' tau_complete <- res |> dplyr::pull(tau) |> unique()
#' var_x_r2_Y_Z <- res |> dplyr::pull(r2_Y_Z)
#' var_x_r2_D_Z <- res |> dplyr::pull(r2_D_Z)
#' var_x_name <- paste0("Drop ", paste0("X", 1:5))
#' plot_sensitivity_r2(tau_complete = tau_complete, df = df, tau_se = tau_se,
#'                     sign_bias = 1,
#'                     var_x_r2_Y_Z = var_x_r2_Y_Z, var_x_r2_D_Z = var_x_r2_D_Z,
#'                     var_x_name = var_x_name,
#'                     title = "Sensitivity of the SCM Estimate (Positive Bias)",
#'                     r2_Y_Z_seq = seq(0, 0.1, length.out = 100),
#'                     r2_D_Z_seq = seq(0, 0.1, length.out = 100),
#'                     text_size = 5)
#'
plot_sensitivity_r2 <- function(
    tau_complete,
    df, tau_se, sign_bias,
    var_x_r2_Y_Z, var_x_r2_D_Z, var_x_name, title = NULL,
    r2_Y_Z_seq = seq(0, 1, length.out = 100), r2_D_Z_seq = seq(0, 1, length.out = 100), text_size = 3,
    plot_t_stat = FALSE, t_stat = NULL, critical_value = NULL) {
  # Input validation
  if (!is.numeric(tau_complete) || length(tau_complete) != 1) {
    stop("'tau_complete' must be a single numeric value.", call. = FALSE)
  }
  if (!is.numeric(df) || length(df) != 1 || df < 1) {
    stop("'df' must be a single positive integer (degrees of freedom).", call. = FALSE)
  }
  if (!is.numeric(tau_se) || length(tau_se) != 1 || tau_se <= 0) {
    stop("'tau_se' must be a single positive numeric value (standard error).", call. = FALSE)
  }
  if (!sign_bias %in% c(-1, 1)) {
    stop("'sign_bias' must be either 1 (positive bias) or -1 (negative bias).", call. = FALSE)
  }
  if (!is.numeric(var_x_r2_Y_Z) || length(var_x_r2_Y_Z) < 1) {
    stop("'var_x_r2_Y_Z' must be a numeric vector with at least one element.", call. = FALSE)
  }
  if (!is.numeric(var_x_r2_D_Z) || length(var_x_r2_D_Z) < 1) {
    stop("'var_x_r2_D_Z' must be a numeric vector with at least one element.", call. = FALSE)
  }
  if (length(var_x_r2_Y_Z) != length(var_x_r2_D_Z)) {
    stop("'var_x_r2_Y_Z' and 'var_x_r2_D_Z' must have the same length.", call. = FALSE)
  }
  if (length(var_x_r2_Y_Z) != length(var_x_name)) {
    stop("'var_x_name' must have the same length as 'var_x_r2_Y_Z' and 'var_x_r2_D_Z'.", call. = FALSE)
  }
  if (plot_t_stat && is.null(t_stat)) {
    stop("'t_stat' must be provided when 'plot_t_stat = TRUE'.", call. = FALSE)
  }
  
  out <- tibble(
    var_name = var_x_name,
    x = var_x_r2_Y_Z,
    y = var_x_r2_D_Z
  )
  tau <- tibble(
    label = paste0("Unadjusted (", round(tau_complete, 1), ")"),
    x = 0, y = 0
  )

  sens_df <- expand.grid(r2_Y_Z_seq, r2_D_Z_seq) |>
    as.data.frame() |>
    mutate(z = tau_complete + sign_bias * (tau_se * sqrt(df * (Var1 * Var2) / (1 - Var2))))

  if (plot_t_stat) {
    sens_df <- sens_df |>
      mutate(
        tau_se_star = tau_se * sqrt(((1 - Var1) / (1 - Var2)) * (df / (df - 1))),
        t_stat = z / tau_se_star
      ) |>
      mutate(z = t_stat)
    tau <- tibble(
      label = paste0("Unadjusted (", round(t_stat, 1), ")"),
      x = 0, y = 0
    )
  }

  p <- sens_df |>
    ggplot() +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1) +
    geom_hline(yintercept = 0, color = "black", alpha = 0.1) +
    geom_contour(aes(x = Var1, y = Var2, z = z), color = "black", alpha = 0.5) +
    metR::geom_text_contour(aes(x = Var1, y = Var2, z = z), size = 3, skip = 0, alpha = 0.5, nudge_x = 0.01, nudge_y = 0.01) +
    geom_point(data = out, aes(x = x, y = y)) +
    geom_text(
      data = out, aes(x = x, y = y, label = var_name),
      hjust = 0, nudge_x = 0.02, size = text_size
    ) +
    geom_point(data = tau, aes(x = x, y = y), shape = 17, color = "#990000") +
    geom_text(
      data = tau, aes(x = x, y = y, label = label), color = "#990000",
      hjust = 0, nudge_x = 0.02, size = text_size
    ) +
    labs(
      x = expression(Partial ~ R^2 ~ of ~ Y ~ "~" ~ Z * "|" * D ~ "," ~ X),
      y = expression(Partial ~ R^2 ~ of ~ Z ~ "~" ~ D * "|" * X),
      title = title
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  if (!is.null(critical_value)) {
    # Generate contour lines using base R
    contour_lines <- contourLines(
      x = r2_Y_Z_seq,
      y = r2_D_Z_seq,
      z = matrix(sens_df$z, nrow = length(r2_Y_Z_seq), ncol = length(r2_D_Z_seq), byrow = FALSE),
      levels = critical_value
    )
    # Convert contour lines to a data frame
    contour_df <- do.call(rbind, lapply(contour_lines, function(contour) {
      data.frame(x = contour$x, y = contour$y, z = critical_value)
    }))
    p <- p + geom_path(data = contour_df, aes(x = x, y = y, group = z), color = "red", size = 1, linetype = "dashed")
  }
  return(p)
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Var1", "Var2", "label", "tau_se_star", "var_name", "x", "y", "z"
  ))
}
