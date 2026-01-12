#' @keywords internal
"_PACKAGE"

#' scmsens: Sensitivity Analysis Tools for Synthetic Control Method
#'
#' The scmsens package provides sensitivity analysis tools for the Synthetic
#' Control Method (SCM) framework. It helps researchers test the robustness
#' of SCM results to potential unobserved donor units, extending the
#' traditional omitted variable bias framework developed by Cinelli and
#' Hazlett (2020) for sensitivity analysis.
#'
#' @section Main Functions:
#' The package provides the following main functions:
#' \itemize{
#'   \item \code{\link{estimate_params}}: Perform leave-one-out sensitivity analysis
#'     by dropping one control unit at a time to estimate bias parameters.
#'   \item \code{\link{estimate_params_partial}}: Estimate sensitivity parameters
#'     with a single control unit that has partially observed data.
#'   \item \code{\link{estimate_params_partial_multi}}: Estimate sensitivity parameters
#'     with multiple control units that have partially observed data.
#'   \item \code{\link{vertreg}}: Run a vertical regression to estimate
#'     treatment effects in SCM settings.
#'   \item \code{\link{vertreg_stacked}}: Run a stacked vertical regression
#'     with a treatment indicator for inference.
#'   \item \code{\link{plot_sensitivity}}: Generate sensitivity contour plots
#'     showing how treatment effects vary with weight and imbalance parameters.
#'   \item \code{\link{plot_sensitivity_r2}}: Generate sensitivity contour plots
#'     using partial R-squared metrics.
#' }
#'
#' @section Data Generation Functions:
#' For simulation studies and testing:
#' \itemize{
#'   \item \code{\link{generate_synth_data_vertreg}}: Generate synthetic data
#'     using a vertical regression model.
#'   \item \code{\link{generate_synth_data_ar}}: Generate synthetic data using
#'     an autoregressive model.
#'   \item \code{\link{generate_synth_data_ife}}: Generate synthetic data using
#'     an interactive fixed effects model.
#' }
#'
#' @section Included Datasets:
#' \itemize{
#'   \item \code{\link{synth_data}}: Synthetic data generated with an interactive
#'     fixed effect model.
#'   \item \code{\link{synth_pre}}, \code{\link{synth_post}}: Pre-period and
#'     post-period data derived from \code{synth_data}.
#'   \item \code{\link{wgermany}}: West German reunification data from
#'     Abadie, Diamond, and Hainmueller (2015).
#'   \item \code{\link{taiwan}}: Taiwan's expulsion from the IMF data from
#'     Lipscy and Lee (2019).
#' }
#'
#' @section Getting Started:
#' A typical workflow involves:
#' \enumerate{
#'   \item Prepare your data in wide format with columns for treated and control units.
#'   \item Split your data into pre-period (\code{data_pre}) and post-period (\code{data_post}).
#'   \item Use \code{\link{estimate_params}} for leave-one-out sensitivity analysis.
#'   \item Visualize results with \code{\link{plot_sensitivity}} or \code{\link{plot_sensitivity_r2}}.
#' }
#'
#' See the package vignette for detailed examples: \code{vignette("simulation_study", package = "scmsens")}
#'
#' @references
#' Liu, N., Shin, S., & Yamauchi, S. (2024). Synthetic Control Method with
#' Missing Pre-treatment Outcomes. Working Paper.
#' \url{https://sooahnshin.com/SCM_Missing.pdf}
#'
#' Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity: Extending
#' omitted variable bias. \emph{Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology)}, 82(1), 39-67.
#' \doi{10.1111/rssb.12348}
#'
#' @docType package
#' @name scmsens-package
NULL
