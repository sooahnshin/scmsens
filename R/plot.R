#' Sensitivity Contour Plot
#'
#' Generate a sensitivity contour plot with the estimated treatment effect, weights, and imbalance.
#'
#' @param tau_complete Estimated treatment effect with the entire data.
#' @param var_x_gamma Weights of the control units.
#' @param var_x_delta Imbalance of the control units.
#' @param var_x_name Names of the control units.
#' @param title Title of the plot
#' @param gamma_seq Sequence of gamma values to plot.
#' @param delta_seq Sequence of delta values to plot.
#' @param text_size Size of the text
#' @param nudge_x Nudge value for the x-axis text
#' @param nudge_y Nudge value for the y-axis text
#' @param repel Whether to use `ggrepel::geom_text_repel()` for text labels.
#' @param tau_round Number of decimal places to round the treatment effect.
#' @param nudge_x2 Nudge value for the x-axis text of the treatment effect.
#' @param nudge_y2 Nudge value for the y-axis text of the treatment effect.
#' @return ggplot2 object
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#'
#' @examples
#' # leave-one-out sensitivity analysis
#' res <- estimate_params(var_y_name = "Y",
#'                        var_x_name = paste0("X", 1:16),
#'                        data_pre = synth_pre,
#'                        data_post = synth_post)
#' tau_complete <- res |> dplyr::pull(tau) |> unique()
#' var_x_gamma <- res |> dplyr::pull(gamma)
#' var_x_delta <- res |> dplyr::pull(delta)
#' var_x_tau <- res |> dplyr::pull(estimate)
#' var_x_name <- paste0("Drop ", paste0("X", 1:16), " (", round(var_x_tau, 1), ")")
#' plot_sensitivity(tau_complete = tau_complete,
#'                  var_x_gamma = var_x_gamma, var_x_delta = var_x_delta, var_x_name = var_x_name,
#'                  title = "Sensitivity of the SCM Estimate", text_size = 4)
#'
#'
#' @seealso [plot_sensitivity_r2()]
plot_sensitivity <- function(
    tau_complete, var_x_gamma, var_x_delta, var_x_name, title = NULL,
    gamma_seq = NULL, delta_seq = NULL, text_size = 3, nudge_x = 0.02, nudge_y = 0,
    repel = FALSE, tau_round = 1, nudge_x2 = NULL, nudge_y2 = NULL) {
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

#' Sensitivity Contour Plot with Partial R Squared
#'
#' Generate a sensitivity contour plot of partial R squared with the estimated treatment effect, weights, and imbalance.
#'
#' @param tau_complete Estimated treatment effect with the entire data.
#' @param df Degree of freedom of the vertical regression with the entire data.
#' @param tau_se Standard error of the treatment effect with the entire data.
#' @param sign_bias Sign of the bias.
#' @param var_x_r2_Y_Z Partial R squared of the control units with the outcome.
#' @param var_x_r2_D_Z Partial R squared of the control units with the treatment indicator.
#' @param var_x_name Names of the control units.
#' @param title Title of the plot
#' @param r2_Y_Z_seq Sequence of x values to plot.
#' @param r2_D_Z_seq Sequence of y values to plot.
#' @param text_size Size of the text
#' @param plot_t_stat Whether to plot t-statistic instead of bias. If true, `t_stat` argument should be provided.
#' @param t_stat t-statistic of the treatment effect with the entire data.
#' @param critical_value Critical value for the t-statistic contour line.
#' @return ggplot2 object
#'
#' @import ggplot2
#' @importFrom metR geom_text_contour
#' @importFrom grDevices contourLines
#'
#' @export
#' @seealso [plot_sensitivity()]
#'
#' @examples
#' # leave-one-out sensitivity analysis
#' fit <- vertreg_stacked(formula = paste("Y ~ -1 + ",
#'                                        paste(paste0("X", 1:16),
#'                                              collapse = " + ")),
#'                        data_pre = synth_pre,
#'                        data_post = synth_post,
#'                        pseudo_inverse = FALSE)
#' df <- fit$df
#' tau_se <- fit$std.error
#' res <- estimate_params(var_y_name = "Y",
#'                        var_x_name = paste0("X", 1:16),
#'                        data_pre = synth_pre,
#'                        data_post = synth_post)
#' # plot only top 5 units with positive bias
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
