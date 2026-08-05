test_that("simplex_lm_fit returns feasible simplex coefficients", {
  dat <- data.frame(
    Y = c(1, 2, 3, 4),
    X1 = c(1, 2, 3, 4),
    X2 = c(4, 3, 2, 1)
  )

  fit <- simplex_lm_fit(Y ~ X1 + X2, dat)
  beta <- coef(fit)

  expect_true(fit$converged)
  expect_equal(sum(beta), 1, tolerance = 1e-8)
  expect_true(all(beta >= -1e-8))
  expect_lt(fit$kkt_residual_scaled, 1e-5)
})

test_that("projected-gradient and quadprog simplex fits agree", {
  set.seed(20260719)
  x <- matrix(stats::rnorm(240), nrow = 60, ncol = 4)
  colnames(x) <- paste0("X", 1:4)
  beta <- c(0.1, 0.2, 0.3, 0.4)
  dat <- data.frame(Y = as.vector(x %*% beta), x)

  fit_pg <- simplex_lm_fit(
    Y ~ X1 + X2 + X3 + X4,
    dat,
    solver = "projected_gradient"
  )
  fit_ref <- simplex_lm_fit(
    Y ~ X1 + X2 + X3 + X4,
    dat,
    solver = "quadprog",
    tol = 1e-10
  )

  expect_true(fit_pg$converged)
  expect_true(fit_ref$converged)
  expect_equal(coef(fit_pg), coef(fit_ref), tolerance = 1e-4)
  expect_lt(abs(fit_pg$objective - fit_ref$objective), 1e-7)
  expect_lt(fit_pg$kkt_residual_scaled, 1e-5)
  expect_lt(fit_ref$kkt_residual_scaled, 1e-5)
})

test_that("quadprog ridge fallback preserves fitted values under rank deficiency", {
  dat <- data.frame(
    Y = c(1, 2, 3, 4),
    X1 = c(1, 2, 3, 4),
    X2 = c(1, 2, 3, 4),
    X3 = c(4, 3, 2, 1)
  )

  fit_pg <- simplex_lm_fit(
    Y ~ X1 + X2 + X3,
    dat,
    solver = "projected_gradient"
  )
  fit_ref <- simplex_lm_fit(
    Y ~ X1 + X2 + X3,
    dat,
    solver = "quadprog"
  )

  expect_equal(fit_pg$fitted.values, fit_ref$fitted.values, tolerance = 1e-4)
  expect_equal(fit_pg$objective, fit_ref$objective, tolerance = 1e-7)
  expect_gte(fit_ref$solver_ridge, 0)
})

test_that("quadprog is invariant to Germany-scale outcomes", {
  dat <- data.frame(
    Y = c(20465, 21602, 22154, 21878, 22371, 23035),
    X1 = c(20471, 21351, 22085, 22540, 23453, 24270),
    X2 = c(19000, 20000, 21000, 22000, 23000, 24000),
    X3 = c(22000, 21800, 22400, 22900, 23100, 23800)
  )

  fit <- simplex_lm_fit(Y ~ X1 + X2 + X3, dat)

  expect_true(fit$converged)
  expect_equal(sum(coef(fit)), 1, tolerance = 1e-8)
  expect_true(all(coef(fit) >= -1e-8))
  expect_lt(fit$kkt_residual_scaled, 1e-5)
})

test_that("Q inverse uses the squared design-matrix tolerance", {
  set.seed(20260720)
  u <- qr.Q(qr(matrix(stats::rnorm(90), nrow = 30, ncol = 3)))
  v <- qr.Q(qr(matrix(stats::rnorm(9), nrow = 3, ncol = 3)))
  x <- u %*% diag(c(2e4, 20, 16)) %*% t(v)
  q <- crossprod(x)
  x_post <- c(1, -0.5, 0.25)

  expected <- sqrt(as.numeric(t(x_post) %*% solve(q, x_post)))

  expect_gt(min(svd(x)$d) / max(svd(x)$d), get_pseudo_inverse_tol())
  expect_lt(min(svd(q)$d) / max(svd(q)$d), get_pseudo_inverse_tol())
  expect_equal(compute_lx(x_post, q), expected, tolerance = 1e-6)
})

test_that("vertreg_simplex computes a treatment effect from simplex weights", {
  data_pre <- data.frame(
    Y = 0.25 * (1:5) + 0.75 * (5:1),
    X1 = 1:5,
    X2 = 5:1
  )
  data_post <- data.frame(Y = 20, X1 = 10, X2 = 20)

  tau <- vertreg_simplex(Y ~ X1 + X2, data_pre, data_post)

  expect_equal(as.numeric(tau), 2.5, tolerance = 1e-3)
})

test_that("partial simplex bound matches a hand-checkable example", {
  data_pre <- data.frame(
    Y = c(0.2, 0.3, 0.5),
    X1 = c(1, 0, 0),
    X2 = c(0, 1, 0),
    Z = c(0, 0, 1)
  )
  data_post <- data.frame(Y = 0, X1 = 1, X2 = 1, Z = 2)

  res <- estimate_params_partial_simplex(
    fm_z_on_x = Z ~ X1 + X2,
    fm_y_on_z_and_x = Y ~ Z + X1 + X2,
    data_pre = data_pre,
    data_post = data_post,
    data_pre_restricted = data_pre,
    pseudo_inverse = TRUE
  )

  expect_equal(res$gamma_const, 0.5, tolerance = 1e-5)
  expect_equal(res$delta, 2, tolerance = 1e-8)
  expect_equal(res$main_term, 1, tolerance = 1e-5)
  expect_equal(res$q_distance, sqrt(0.125), tolerance = 1e-5)
  expect_equal(res$L_X, sqrt(2), tolerance = 1e-8)
  expect_equal(res$simplex_penalty, 0.5, tolerance = 1e-5)
  expect_equal(res$bound_deterministic, 1.5, tolerance = 1e-5)
  expect_equal(res$product_signed, 1, tolerance = 1e-5)
  expect_equal(res$projection_signed, -0.5, tolerance = 1e-5)
  expect_equal(res$projection_abs, 0.5, tolerance = 1e-5)
  expect_equal(res$stochastic_signed, 0, tolerance = 1e-8)
  expect_equal(res$cs_ratio, 1, tolerance = 1e-5)
  expect_lt(res$full_kkt_residual_scaled, 1e-5)
  expect_lt(res$restricted_kkt_residual_scaled, 1e-5)

  tau_complete <- vertreg_simplex(
    Y ~ Z + X1 + X2,
    data_pre,
    data_post
  )
  tau_deletion <- vertreg_simplex(
    Y ~ X1 + X2,
    data_pre,
    data_post
  )
  expect_equal(
    as.numeric(tau_deletion - tau_complete),
    res$product_signed + res$projection_signed + res$stochastic_signed,
    tolerance = 1e-8
  )
})

test_that("multi-unit simplex collapsed product equals vector product", {
  data_pre <- data.frame(
    Y = c(0.1, 0.2, 0.3, 0.4),
    X1 = c(1, 0, 0, 0),
    X2 = c(0, 1, 0, 0),
    Z1 = c(0, 0, 1, 0),
    Z2 = c(0, 0, 0, 1)
  )
  data_post <- data.frame(Y = 0, X1 = 1, X2 = 1, Z1 = 2, Z2 = 3)

  res <- estimate_params_partial_multi_simplex(
    var_z_name = c("Z1", "Z2"),
    var_x_name = c("X1", "X2"),
    var_y_name = "Y",
    data_pre = data_pre,
    data_post = data_post,
    data_pre_restricted = data_pre,
    pseudo_inverse = TRUE
  )

  direct_product <- sum(res$component_gamma_const * res$component_delta)

  expect_equal(res$gamma_total * res$delta_collapsed, direct_product, tolerance = 1e-8)
  expect_equal(res$product_signed, direct_product, tolerance = 1e-8)
  expect_equal(res$main_term, abs(direct_product), tolerance = 1e-8)
})

test_that("plot_sensitivity_simplex returns a ggplot object", {
  plot <- plot_sensitivity_simplex(
    var_x_gamma_const = c(0.1, 0.4),
    var_x_delta = c(-2, 3),
    var_x_name = c("A", "B"),
    simplex_penalty = 0,
    var_x_bound = c(1.24, 5.67),
    bound_digits = 1
  )

  expect_s3_class(plot, "ggplot")
  benchmark_layers <- Filter(
    function(layer) "var_name" %in% names(layer$data),
    plot$layers
  )
  expect_true(length(benchmark_layers) >= 1)
  expect_equal(
    benchmark_layers[[1]]$data$var_name,
    c("A (1.2)", "B (5.7)")
  )
  expect_equal(plot$data$bound, abs(plot$data$Var1 * plot$data$Var2))
})

test_that("plot_sensitivity_simplex validates benchmark-specific bounds", {
  expect_error(
    plot_sensitivity_simplex(
      var_x_gamma_const = c(0.1, 0.4),
      var_x_delta = c(-2, 3),
      var_x_bound = 1
    ),
    "one nonnegative finite value or NA per point"
  )
  expect_error(
    plot_sensitivity_simplex(
      var_x_gamma_const = c(0.1, 0.4),
      var_x_delta = c(-2, 3),
      var_x_bound = c(1, -1)
    ),
    "one nonnegative finite value or NA per point"
  )
})
