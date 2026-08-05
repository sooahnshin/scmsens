#' Project a Vector onto the Simplex
#'
#' @param x Numeric vector.
#' @param z Simplex radius. Defaults to 1.
#'
#' @return Numeric vector on the simplex.
simplex_project <- function(x, z = 1) {
  if (!is.numeric(x) || length(x) < 1 || anyNA(x)) {
    stop("'x' must be a non-empty numeric vector without missing values.",
         call. = FALSE)
  }
  if (!is.numeric(z) || length(z) != 1 || is.na(z) || z <= 0) {
    stop("'z' must be a positive number.", call. = FALSE)
  }

  if (length(x) == 1) {
    return(z)
  }

  u <- sort(x, decreasing = TRUE)
  cssv <- cumsum(u) - z
  rho <- max(which(u - cssv / seq_along(u) > 0))
  theta <- cssv[rho] / rho
  pmax(x - theta, 0)
}

matrix_pseudo_inverse <- function(x, tol = get_pseudo_inverse_tol()) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (nrow(x) == 0 || ncol(x) == 0) {
    stop("'x' must be a non-empty matrix.", call. = FALSE)
  }
  svd_x <- svd(x)
  if (length(svd_x$d) == 0 || max(svd_x$d) == 0) {
    return(matrix(0, ncol = nrow(x), nrow = ncol(x)))
  }
  keep <- svd_x$d > tol * max(svd_x$d)
  out <- matrix(0, nrow = ncol(x), ncol = nrow(x))
  if (any(keep)) {
    out <- svd_x$v[, keep, drop = FALSE] %*%
      (diag(1 / svd_x$d[keep], nrow = sum(keep)) %*%
         t(svd_x$u[, keep, drop = FALSE]))
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

unrestricted_lm_coef <- function(formula, data, pseudo_inverse = FALSE) {
  formula <- update(as.formula(formula), . ~ . - 1)
  if (isTRUE(pseudo_inverse)) {
    coef(lm_pseudo(formula, data = data))
  } else {
    coef(lm(formula, data = data))
  }
}

simplex_kkt_diagnostics <- function(x, y, beta, active_tol = 1e-8) {
  beta <- as.numeric(beta)
  gradient <- as.numeric(crossprod(x, as.vector(x %*% beta) - y))
  active <- beta > active_tol
  if (!any(active)) {
    active[which.max(beta)] <- TRUE
  }
  lambda <- -mean(gradient[active])
  reduced_gradient <- gradient + lambda
  inactive <- !active
  gradient_scale <- 1 + max(abs(gradient))

  equality_violation <- abs(sum(beta) - 1)
  nonnegativity_violation <- max(c(0, -beta))
  active_stationarity <- max(abs(reduced_gradient[active]))
  dual_violation <- if (any(inactive)) {
    max(c(0, -reduced_gradient[inactive]))
  } else {
    0
  }
  complementarity <- max(abs(beta * reduced_gradient))
  singular_values <- svd(x, nu = 0, nv = 0)$d
  lipschitz <- if (length(singular_values) == 0 || max(singular_values) == 0) {
    1
  } else {
    max(singular_values)^2
  }
  projected_gradient <- max(abs(
    beta - simplex_project(beta - gradient / lipschitz)
  ))
  kkt_residual_scaled <- max(
    equality_violation,
    nonnegativity_violation,
    active_stationarity / gradient_scale,
    dual_violation / gradient_scale,
    complementarity / gradient_scale,
    projected_gradient
  )

  list(
    lambda = as.numeric(lambda),
    equality_violation = as.numeric(equality_violation),
    nonnegativity_violation = as.numeric(nonnegativity_violation),
    active_stationarity = as.numeric(active_stationarity),
    dual_violation = as.numeric(dual_violation),
    complementarity = as.numeric(complementarity),
    projected_gradient = as.numeric(projected_gradient),
    kkt_residual_scaled = as.numeric(kkt_residual_scaled),
    active_size = sum(active)
  )
}

solve_simplex_active_set <- function(x, y, tol = 1e-10, maxit = 1000) {
  p <- ncol(x)
  q <- crossprod(x)
  cxy <- as.numeric(crossprod(x, y))

  solve_on_active <- function(active) {
    k <- length(active)
    kkt <- rbind(
      cbind(q[active, active, drop = FALSE], rep(1, k)),
      c(rep(1, k), 0)
    )
    rhs <- c(cxy[active], 1)
    sol <- as.numeric(matrix_pseudo_inverse(kkt, tol = tol) %*% rhs)
    beta <- rep(0, p)
    beta[active] <- sol[seq_len(k)]
    list(beta = beta, lambda = sol[k + 1])
  }

  vertex_objective <- vapply(seq_len(p), function(j) {
    sum((y - x[, j])^2)
  }, numeric(1))
  beta <- rep(0, p)
  beta[which.min(vertex_objective)] <- 1
  active <- which(beta > tol)

  for (iter in seq_len(maxit)) {
    candidate <- solve_on_active(active)
    beta_candidate <- candidate$beta

    while (any(beta_candidate[active] < -tol) && length(active) > 1) {
      negative <- active[beta_candidate[active] <= tol]
      alpha <- min(beta[negative] / (beta[negative] - beta_candidate[negative]))
      beta <- beta + alpha * (beta_candidate - beta)
      beta[abs(beta) < tol] <- 0
      active <- which(beta > tol)
      candidate <- solve_on_active(active)
      beta_candidate <- candidate$beta
    }

    beta <- pmax(beta_candidate, 0)
    if (sum(beta) <= 0) {
      beta[active] <- 1 / length(active)
    } else {
      beta <- beta / sum(beta)
    }
    active <- which(beta > tol)
    candidate <- solve_on_active(active)
    lambda <- candidate$lambda

    inactive <- setdiff(seq_len(p), active)
    reduced_cost <- if (length(inactive) > 0) {
      as.numeric(q[inactive, , drop = FALSE] %*% beta - cxy[inactive] + lambda)
    } else {
      numeric()
    }

    if (length(inactive) == 0 || all(reduced_cost >= -tol)) {
      return(list(beta = beta, converged = TRUE, iterations = iter, solver = "active_set"))
    }

    active <- sort(c(active, inactive[which.min(reduced_cost)]))
  }

  list(beta = beta, converged = FALSE, iterations = maxit, solver = "active_set")
}

solve_simplex_constr_optim <- function(x, y, tol = 1e-10, maxit = 5000) {
  p <- ncol(x)
  if (p == 1) {
    return(list(beta = 1, converged = TRUE, iterations = 0L, solver = "constrOptim"))
  }
  x_base <- x[, p, drop = FALSE]
  x_free <- x[, seq_len(p - 1), drop = FALSE] - as.vector(x_base)
  objective <- function(theta) {
    residual <- y - as.vector(x_base) - as.vector(x_free %*% theta)
    sum(residual^2)
  }
  gradient <- function(theta) {
    residual <- y - as.vector(x_base) - as.vector(x_free %*% theta)
    as.numeric(-2 * crossprod(x_free, residual))
  }
  start <- rep(1 / p, p - 1)
  ui <- rbind(diag(p - 1), -rep(1, p - 1))
  ci <- c(rep(0, p - 1), -1)
  opt <- tryCatch(
    stats::constrOptim(
      theta = start,
      f = objective,
      grad = gradient,
      ui = ui,
      ci = ci,
      control = list(maxit = as.integer(maxit), reltol = tol)
    ),
    error = function(e) e
  )
  if (inherits(opt, "error")) {
    return(list(beta = simplex_project(rep(1 / p, p)), converged = FALSE, iterations = maxit, solver = "constrOptim"))
  }
  beta <- simplex_project(c(opt$par, 1 - sum(opt$par)))
  list(
    beta = beta,
    converged = is.null(opt$convergence) || identical(opt$convergence, 0L),
    iterations = opt$outer.iterations %||% opt$counts[["function"]] %||% NA_integer_,
    solver = "constrOptim"
  )
}

solve_simplex_quadprog <- function(x, y, tol = 1e-10) {
  p <- ncol(x)
  if (p == 1) {
    return(list(
      beta = 1,
      converged = TRUE,
      iterations = 0L,
      solver = "quadprog",
      ridge = 0
    ))
  }

  x_base <- x[, p]
  x_free <- x[, seq_len(p - 1), drop = FALSE] - x_base
  y_free <- y - x_base
  q <- crossprod(x_free)
  cxy <- as.numeric(crossprod(x_free, y_free))
  objective_scale <- max(c(abs(q), 1))
  q_scaled <- q / objective_scale
  cxy_scaled <- cxy / objective_scale
  ridge_candidates <- c(
    0,
    .Machine$double.eps,
    1e-14,
    1e-12,
    1e-10,
    tol,
    sqrt(.Machine$double.eps)
  )
  ridge_candidates <- unique(sort(ridge_candidates))
  constraints <- cbind(diag(p - 1), -rep(1, p - 1))
  rhs <- c(rep(0, p - 1), -1)

  for (ridge in ridge_candidates) {
    fit <- tryCatch(
      quadprog::solve.QP(
        Dmat = 2 * (q_scaled + diag(ridge, p - 1)),
        dvec = 2 * cxy_scaled,
        Amat = constraints,
        bvec = rhs,
        meq = 0
      ),
      error = function(e) e
    )
    if (!inherits(fit, "error")) {
      beta <- c(as.numeric(fit$solution), 1 - sum(fit$solution))
      beta[abs(beta) < tol] <- 0
      beta <- pmax(beta, 0)
      beta <- beta / sum(beta)
      return(list(
        beta = beta,
        converged = TRUE,
        iterations = as.integer(fit$iterations[[1]]),
        solver = "quadprog",
        ridge = as.numeric(ridge * objective_scale)
      ))
    }
  }

  stop(
    "quadprog could not solve the simplex regression, even after numerical ridge stabilization.",
    call. = FALSE
  )
}

#' Fit Linear Regression with Simplex-Constrained Coefficients
#'
#' @param formula A model formula. The intercept is removed before fitting.
#' @param data A data.frame containing the variables in `formula`.
#' @param tol Numerical tolerance used by the selected solver.
#' @param maxit Maximum iterations for projected gradient or `constrOptim`.
#' @param step_scale Multiplicative factor for the inverse-Lipschitz step size.
#' @param solver Optimization method. The default uses `"quadprog"`;
#'   `"projected_gradient"` and `"constrOptim"` are retained as diagnostic
#'   alternatives.
#'
#' @return A `simplex_lm` object.
#'
#' @export
simplex_lm_fit <- function(formula, data, tol = 1e-8, maxit = 5000,
                           step_scale = 1,
                           solver = c("quadprog", "projected_gradient", "constrOptim")) {
  solver <- match.arg(solver)
  if (is.null(formula) || (!is.character(formula) && !inherits(formula, "formula"))) {
    stop("'formula' must be a character string or formula object.",
         call. = FALSE)
  }
  if (is.null(data) || nrow(data) < 1) {
    stop("'data' must be a non-empty data.frame, tibble, or matrix.",
         call. = FALSE)
  }
  if (is.matrix(data) || tibble::is_tibble(data)) {
    data <- as.data.frame(data)
  }
  if (!is.numeric(tol) || length(tol) != 1 || is.na(tol) || tol <= 0) {
    stop("'tol' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || is.na(maxit) || maxit < 1) {
    stop("'maxit' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(step_scale) || length(step_scale) != 1 ||
      is.na(step_scale) || step_scale <= 0) {
    stop("'step_scale' must be a positive number.", call. = FALSE)
  }

  fm <- update(as.formula(formula), . ~ . - 1)
  mf <- model.frame(fm, data)
  x <- model.matrix(fm, mf)
  y <- model.response(mf)
  p <- ncol(x)
  if (p < 1) {
    stop("The formula must contain at least one predictor.", call. = FALSE)
  }
  if (anyNA(x) || anyNA(y)) {
    stop("'formula' variables must not contain missing values.", call. = FALSE)
  }

  sv <- svd(x, nu = 0, nv = 0)$d
  beta <- rep(0, p)
  converged <- FALSE
  iter <- 0L
  step_size <- NA_real_
  solver_used <- solver
  solver_ridge <- 0

  if (p == 1) {
    beta <- 1
    converged <- TRUE
    solver_used <- "closed_form"
  } else if (identical(solver, "constrOptim")) {
    solver_fit <- solve_simplex_constr_optim(x, y, tol = tol, maxit = maxit)
    beta <- solver_fit$beta
    converged <- solver_fit$converged
    iter <- solver_fit$iterations
    solver_used <- solver_fit$solver
  } else if (identical(solver, "quadprog")) {
    solver_fit <- solve_simplex_quadprog(x, y, tol = tol)
    beta <- solver_fit$beta
    converged <- solver_fit$converged
    iter <- solver_fit$iterations
    solver_used <- solver_fit$solver
    solver_ridge <- solver_fit$ridge
  } else {
    lipschitz <- if (length(sv) == 0 || max(sv) == 0) 1 else max(sv)^2
    step_size <- step_scale / lipschitz
    beta <- rep(1 / p, p)
    yk <- beta
    tk <- 1
    previous_objective <- Inf
    for (i in seq_len(as.integer(maxit))) {
      grad <- as.vector(crossprod(x, as.vector(x %*% yk) - y))
      beta_new <- simplex_project(yk - step_size * grad)
      objective_new <- sum((y - as.vector(x %*% beta_new))^2)
      tk_new <- (1 + sqrt(1 + 4 * tk^2)) / 2
      yk <- beta_new + ((tk - 1) / tk_new) * (beta_new - beta)
      iter <- i
      if (is.finite(previous_objective) &&
          abs(previous_objective - objective_new) <= tol * (1 + abs(previous_objective))) {
        candidate_kkt <- simplex_kkt_diagnostics(
          x,
          y,
          beta_new,
          active_tol = max(sqrt(.Machine$double.eps), tol)
        )
        kkt_tol <- max(1e-8, sqrt(tol) / 10)
        if (candidate_kkt$kkt_residual_scaled <= kkt_tol) {
          beta <- beta_new
          converged <- TRUE
          break
        }
      }
      beta <- beta_new
      tk <- tk_new
      previous_objective <- objective_new
    }
  }

  names(beta) <- colnames(x)
  fitted_values <- as.vector(x %*% beta)
  residuals <- y - fitted_values
  objective <- sum(residuals^2)
  rank <- if (length(sv) == 0 || max(sv) == 0) 0L else sum(sv > get_pseudo_inverse_tol() * max(sv))
  kkt <- simplex_kkt_diagnostics(
    x,
    y,
    beta,
    active_tol = max(sqrt(.Machine$double.eps), tol)
  )

  res <- list(
    coefficients = beta,
    fitted.values = fitted_values,
    residuals = residuals,
    objective = objective,
    formula = fm,
    data = data,
    model_matrix = x,
    response = y,
    converged = converged,
    iterations = iter,
    tol = tol,
    maxit = as.integer(maxit),
    step_size = step_size,
    solver = solver_used,
    solver_ridge = solver_ridge,
    singular_values = sv,
    rank = rank,
    simplex_sum = sum(beta),
    simplex_min = min(beta),
    kkt = kkt,
    kkt_residual_scaled = kkt$kkt_residual_scaled,
    projected_gradient_residual = kkt$projected_gradient
  )
  class(res) <- "simplex_lm"
  res
}

#' @export
coef.simplex_lm <- function(object, ...) {
  object$coefficients
}

#' @export
predict.simplex_lm <- function(object, newdata, ...) {
  if (is.matrix(newdata) || tibble::is_tibble(newdata)) {
    newdata <- as.data.frame(newdata)
  }
  x <- model.matrix(object$formula, newdata)
  as.vector(x %*% object$coefficients)
}

#' @export
summary.simplex_lm <- function(object, ...) {
  tibble::tibble(
    term = names(object$coefficients),
    estimate = as.numeric(object$coefficients)
  )
}

#' Vertical Regression with Simplex-Constrained Weights
#'
#' @param formula A vertical-regression formula. The intercept is removed before fitting.
#' @param data_pre Pre-treatment data.
#' @param data_post Post-treatment data.
#' @param ... Passed to [simplex_lm_fit()].
#'
#' @return Numeric treatment-effect estimates, one for each post-treatment row.
#'
#' @export
vertreg_simplex <- function(formula, data_pre, data_post, ...) {
  if (is.matrix(data_pre) || tibble::is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) || tibble::is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
  }
  y_var_name <- all.vars(as.formula(formula))[1]
  fm <- update(as.formula(formula), . ~ . - 1)
  fit <- simplex_lm_fit(fm, data_pre, ...)
  y0 <- predict(fit, newdata = data_post)
  y1 <- data_post[, y_var_name]
  as.vector(y1 - y0)
}

make_rhs_formula <- function(lhs, rhs) {
  as.formula(paste(lhs, "~ -1 +", paste(rhs, collapse = " + ")))
}

get_formula_lhs <- function(formula) {
  all.vars(as.formula(formula))[1]
}

get_formula_rhs <- function(formula) {
  attr(stats::terms(update(as.formula(formula), . ~ . - 1)), "term.labels")
}

q_distance_to_simplex <- function(beta_res, beta_simplex, q) {
  beta_res <- as.numeric(beta_res)
  beta_simplex <- as.numeric(beta_simplex)
  if (length(beta_res) != length(beta_simplex)) {
    stop("'beta_res' and 'beta_simplex' must have the same length.",
         call. = FALSE)
  }
  diff <- beta_res - beta_simplex
  sqrt(max(0, as.numeric(t(diff) %*% q %*% diff)))
}

compute_lx <- function(x_post, q, tol = get_pseudo_inverse_tol()^2) {
  x_post <- as.numeric(x_post)
  q_pinv <- matrix_pseudo_inverse(q, tol = tol)
  sqrt(max(0, as.numeric(t(x_post) %*% q_pinv %*% x_post)))
}

compute_simplex_bound <- function(gamma_const, delta, L_X, q_distance) {
  main_term <- abs(gamma_const * delta)
  simplex_penalty <- L_X * q_distance
  list(
    main_term = as.numeric(main_term),
    simplex_penalty = as.numeric(simplex_penalty),
    bound_deterministic = as.numeric(main_term + simplex_penalty)
  )
}

compute_simplex_stochastic_term <- function(x_pre, x_post, q,
                                            fit_full_const,
                                            tol = get_pseudo_inverse_tol()^2) {
  full_x <- fit_full_const$model_matrix
  same_rows <- nrow(full_x) == nrow(x_pre) &&
    identical(rownames(full_x), rownames(x_pre))
  if (!same_rows) {
    return(NA_real_)
  }

  adjustment <- matrix_pseudo_inverse(q, tol = tol) %*%
    crossprod(x_pre, fit_full_const$residuals)
  -as.numeric(sum(x_post * adjustment))
}

#' Estimate Simplex Sensitivity Parameters for One Partially Observed Unit
#'
#' @param fm_z_on_x Formula for projecting the partially observed unit on observed units.
#' @param fm_y_on_z_and_x Formula for the full vertical regression.
#' @param data_pre Rows used for regressions involving the partially observed unit.
#' @param data_post One-row post-treatment data.
#' @param data_pre_restricted Optional rows used for restricted `Y ~ X` quantities.
#' @param pseudo_inverse Logical; use pseudo-inverse for unrestricted projections.
#' @param ... Passed to [simplex_lm_fit()].
#'
#' @return A named list with simplex bound components.
#'
#' @export
estimate_params_partial_simplex <- function(
    fm_z_on_x,
    fm_y_on_z_and_x,
    data_pre,
    data_post,
    data_pre_restricted = data_pre,
    pseudo_inverse = FALSE,
    ...) {
  if (is.matrix(data_pre) || tibble::is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) || tibble::is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
  }
  if (is.matrix(data_pre_restricted) || tibble::is_tibble(data_pre_restricted)) {
    data_pre_restricted <- as.data.frame(data_pre_restricted)
  }
  if (nrow(data_post) != 1) {
    stop("'data_post' must contain exactly one row.", call. = FALSE)
  }

  var_z_name <- get_formula_lhs(fm_z_on_x)
  var_x_name <- get_formula_rhs(fm_z_on_x)
  var_y_name <- get_formula_lhs(fm_y_on_z_and_x)

  fm_z_on_x <- update(as.formula(fm_z_on_x), . ~ . - 1)
  fit_zx <- if (isTRUE(pseudo_inverse)) {
    lm_pseudo(fm_z_on_x, data = data_pre)
  } else {
    lm(fm_z_on_x, data = data_pre)
  }
  pred_zT <- predict(fit_zx, newdata = data_post)
  delta <- as.numeric(data_post[1, var_z_name] - pred_zT)

  fit_full_const <- simplex_lm_fit(
    update(as.formula(fm_y_on_z_and_x), . ~ . - 1),
    data = data_pre,
    ...
  )
  gamma_const <- coef(fit_full_const)[var_z_name]

  fm_y_on_x <- make_rhs_formula(var_y_name, var_x_name)
  fit_res_const <- simplex_lm_fit(fm_y_on_x, data = data_pre_restricted, ...)
  beta_res_const <- coef(fit_res_const)
  beta_res <- unrestricted_lm_coef(
    fm_y_on_x,
    data = data_pre_restricted,
    pseudo_inverse = pseudo_inverse
  )
  beta_res <- beta_res[names(beta_res_const)]

  x_pre <- model.matrix(fm_y_on_x, data_pre_restricted)
  x_post <- model.matrix(fm_y_on_x, data_post)[1, , drop = TRUE]
  q <- crossprod(x_pre)
  q_distance <- q_distance_to_simplex(beta_res, beta_res_const, q)
  L_X <- compute_lx(x_post, q)
  bound <- compute_simplex_bound(gamma_const, delta, L_X, q_distance)
  product_signed <- as.numeric(gamma_const * delta)
  projection_signed <- as.numeric(
    sum(x_post * (as.numeric(beta_res) - as.numeric(beta_res_const)))
  )
  projection_abs <- abs(projection_signed)
  stochastic_signed <- compute_simplex_stochastic_term(
    x_pre,
    x_post,
    q,
    fit_full_const
  )
  cs_ratio <- if (bound$simplex_penalty > sqrt(.Machine$double.eps)) {
    projection_abs / bound$simplex_penalty
  } else if (projection_abs <= sqrt(.Machine$double.eps)) {
    0
  } else {
    NA_real_
  }

  list(
    gamma_const = as.numeric(gamma_const),
    gamma = as.numeric(gamma_const),
    delta = as.numeric(delta),
    product_signed = product_signed,
    main_term = bound$main_term,
    simplex_penalty = bound$simplex_penalty,
    bound_deterministic = bound$bound_deterministic,
    projection_signed = projection_signed,
    projection_abs = projection_abs,
    stochastic_signed = stochastic_signed,
    cs_ratio = as.numeric(cs_ratio),
    q_distance = as.numeric(q_distance),
    L_X = as.numeric(L_X),
    beta_res = beta_res,
    beta_res_const = beta_res_const,
    beta_res_distance_l2 = sqrt(sum((as.numeric(beta_res) - as.numeric(beta_res_const))^2)),
    simplex_solver_status = if (isTRUE(fit_res_const$converged) && isTRUE(fit_full_const$converged)) "converged" else "maxit",
    full_kkt_residual_scaled = fit_full_const$kkt_residual_scaled,
    restricted_kkt_residual_scaled = fit_res_const$kkt_residual_scaled,
    simplex_solver_ridge = max(
      fit_res_const$solver_ridge,
      fit_full_const$solver_ridge
    ),
    full_fit = fit_full_const,
    restricted_fit = fit_res_const
  )
}

#' Estimate Simplex Sensitivity Parameters for Multiple Partially Observed Units
#'
#' @inheritParams estimate_params_partial_multi
#' @param data_pre_restricted Optional rows used for restricted `Y ~ X` quantities.
#' @param pseudo_inverse Logical; use pseudo-inverse for unrestricted projections.
#' @param ... Passed to [simplex_lm_fit()].
#'
#' @return A named list with collapsed-block simplex bound components.
#'
#' @export
estimate_params_partial_multi_simplex <- function(
    var_z_name,
    var_x_name,
    var_y_name,
    data_pre,
    data_post,
    data_pre_restricted = data_pre,
    pseudo_inverse = FALSE,
    ...) {
  if (is.matrix(data_pre) || tibble::is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) || tibble::is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
  }
  if (is.matrix(data_pre_restricted) || tibble::is_tibble(data_pre_restricted)) {
    data_pre_restricted <- as.data.frame(data_pre_restricted)
  }
  if (nrow(data_post) != 1) {
    stop("'data_post' must contain exactly one row.", call. = FALSE)
  }

  fm_y_on_z_and_x <- make_rhs_formula(var_y_name, c(var_z_name, var_x_name))
  fit_full_const <- simplex_lm_fit(fm_y_on_z_and_x, data = data_pre, ...)
  component_gamma <- coef(fit_full_const)[var_z_name]

  component_delta <- purrr::map_dbl(var_z_name, function(z_name) {
    fm_z_on_x <- make_rhs_formula(z_name, var_x_name)
    fit_zx <- if (isTRUE(pseudo_inverse)) {
      lm_pseudo(fm_z_on_x, data = data_pre)
    } else {
      lm(fm_z_on_x, data = data_pre)
    }
    as.numeric(data_post[1, z_name] - predict(fit_zx, newdata = data_post))
  })
  names(component_delta) <- var_z_name
  component_bias <- component_gamma * component_delta
  gamma_total <- sum(component_gamma)
  delta_collapsed <- if (abs(gamma_total) > sqrt(.Machine$double.eps)) {
    sum(component_bias) / gamma_total
  } else {
    NA_real_
  }

  fm_y_on_x <- make_rhs_formula(var_y_name, var_x_name)
  fit_res_const <- simplex_lm_fit(fm_y_on_x, data = data_pre_restricted, ...)
  beta_res_const <- coef(fit_res_const)
  beta_res <- unrestricted_lm_coef(
    fm_y_on_x,
    data = data_pre_restricted,
    pseudo_inverse = pseudo_inverse
  )
  beta_res <- beta_res[names(beta_res_const)]

  x_pre <- model.matrix(fm_y_on_x, data_pre_restricted)
  x_post <- model.matrix(fm_y_on_x, data_post)[1, , drop = TRUE]
  q <- crossprod(x_pre)
  q_distance <- q_distance_to_simplex(beta_res, beta_res_const, q)
  L_X <- compute_lx(x_post, q)
  bound <- compute_simplex_bound(1, sum(component_bias), L_X, q_distance)
  product_signed <- as.numeric(sum(component_bias))
  projection_signed <- as.numeric(
    sum(x_post * (as.numeric(beta_res) - as.numeric(beta_res_const)))
  )
  projection_abs <- abs(projection_signed)
  stochastic_signed <- compute_simplex_stochastic_term(
    x_pre,
    x_post,
    q,
    fit_full_const
  )
  cs_ratio <- if (bound$simplex_penalty > sqrt(.Machine$double.eps)) {
    projection_abs / bound$simplex_penalty
  } else if (projection_abs <= sqrt(.Machine$double.eps)) {
    0
  } else {
    NA_real_
  }

  list(
    gamma_const = as.numeric(gamma_total),
    gamma = as.numeric(gamma_total),
    delta = as.numeric(delta_collapsed),
    product_signed = product_signed,
    main_term = abs(sum(component_bias)),
    simplex_penalty = bound$simplex_penalty,
    bound_deterministic = abs(sum(component_bias)) + bound$simplex_penalty,
    projection_signed = projection_signed,
    projection_abs = projection_abs,
    stochastic_signed = stochastic_signed,
    cs_ratio = as.numeric(cs_ratio),
    q_distance = as.numeric(q_distance),
    L_X = as.numeric(L_X),
    beta_res = beta_res,
    beta_res_const = beta_res_const,
    beta_res_distance_l2 = sqrt(sum((as.numeric(beta_res) - as.numeric(beta_res_const))^2)),
    simplex_solver_status = if (isTRUE(fit_res_const$converged) && isTRUE(fit_full_const$converged)) "converged" else "maxit",
    full_kkt_residual_scaled = fit_full_const$kkt_residual_scaled,
    restricted_kkt_residual_scaled = fit_res_const$kkt_residual_scaled,
    simplex_solver_ridge = max(
      fit_res_const$solver_ridge,
      fit_full_const$solver_ridge
    ),
    component_gamma_const = component_gamma,
    component_delta = component_delta,
    component_bias = component_bias,
    gamma_total = as.numeric(gamma_total),
    delta_collapsed = as.numeric(delta_collapsed),
    full_fit = fit_full_const,
    restricted_fit = fit_res_const
  )
}

#' Leave-One-Out Simplex Sensitivity Analysis
#'
#' @inheritParams estimate_params
#' @param pseudo_inverse Logical; use pseudo-inverse for unrestricted projections.
#' @param ... Passed to [simplex_lm_fit()].
#'
#' @return A tibble with simplex sensitivity quantities for each dropped unit.
#'
#' @export
estimate_params_simplex <- function(var_y_name,
                                    var_x_name,
                                    data_pre,
                                    data_post,
                                    pseudo_inverse = FALSE,
                                    ...) {
  if (is.matrix(data_pre) || tibble::is_tibble(data_pre)) {
    data_pre <- as.data.frame(data_pre)
  }
  if (is.matrix(data_post) || tibble::is_tibble(data_post)) {
    data_post <- as.data.frame(data_post)
  }
  tau_complete <- vertreg_simplex(
    formula = make_rhs_formula(var_y_name, var_x_name),
    data_pre = data_pre,
    data_post = data_post,
    ...
  )
  purrr::map_dfr(var_x_name, function(z) {
    x_names <- setdiff(var_x_name, z)
    res <- estimate_params_partial_simplex(
      fm_z_on_x = make_rhs_formula(z, x_names),
      fm_y_on_z_and_x = make_rhs_formula(var_y_name, c(z, x_names)),
      data_pre = data_pre,
      data_post = data_post,
      data_pre_restricted = data_pre,
      pseudo_inverse = pseudo_inverse,
      ...
    )
    tau_deletion <- vertreg_simplex(
      formula = make_rhs_formula(var_y_name, x_names),
      data_pre = data_pre,
      data_post = data_post,
      ...
    )
    tibble::tibble(
      dropped_unit = z,
      tau_complete = as.numeric(tau_complete),
      estimate = as.numeric(tau_deletion),
      gamma_const = res$gamma_const,
      delta = res$delta,
      main_term = res$main_term,
      simplex_penalty = res$simplex_penalty,
      bound_deterministic = res$bound_deterministic,
      q_distance = res$q_distance,
      L_X = res$L_X
    )
  })
}

#' Plot Simplex Sensitivity Upper Bounds
#'
#' @param var_x_gamma_const Numeric vector of simplex weights.
#' @param var_x_delta Numeric vector of imbalance values.
#' @param var_x_name Character labels for points.
#' @param simplex_penalty Nonnegative scalar added to `abs(gamma * delta)`.
#' @param var_x_bound Optional numeric vector of benchmark-specific deterministic
#'   bounds. When supplied, each value is appended to its point label in
#'   parentheses; it does not change the common contour surface.
#' @param bound_digits Number of decimal places used for parenthetical bounds.
#' @param gamma_seq Optional grid for gamma.
#' @param delta_seq Optional grid for delta.
#' @param title Plot title.
#' @param text_size Label text size.
#' @param repel Use `ggrepel` for labels.
#'
#' @return A ggplot object.
#'
#' @export
plot_sensitivity_simplex <- function(
    var_x_gamma_const,
    var_x_delta,
    var_x_name = NULL,
    simplex_penalty = 0,
    var_x_bound = NULL,
    bound_digits = 1,
    gamma_seq = NULL,
    delta_seq = NULL,
    title = NULL,
    text_size = 3,
    repel = FALSE) {
  if (!is.numeric(var_x_gamma_const) || !is.numeric(var_x_delta)) {
    stop("'var_x_gamma_const' and 'var_x_delta' must be numeric vectors.",
         call. = FALSE)
  }
  if (length(var_x_gamma_const) != length(var_x_delta)) {
    stop("'var_x_gamma_const' and 'var_x_delta' must have the same length.",
         call. = FALSE)
  }
  if (is.null(var_x_name)) {
    var_x_name <- seq_along(var_x_gamma_const)
  }
  if (length(var_x_name) != length(var_x_gamma_const)) {
    stop("'var_x_name' must have the same length as the numeric vectors.",
         call. = FALSE)
  }
  if (!is.numeric(simplex_penalty) || length(simplex_penalty) != 1 ||
      is.na(simplex_penalty) || simplex_penalty < 0) {
    stop("'simplex_penalty' must be a nonnegative scalar.", call. = FALSE)
  }
  if (!is.null(var_x_bound)) {
    if (!is.numeric(var_x_bound) ||
        length(var_x_bound) != length(var_x_gamma_const) ||
        any(!is.na(var_x_bound) & (!is.finite(var_x_bound) | var_x_bound < 0))) {
      stop(
        "'var_x_bound' must contain one nonnegative finite value or NA per point.",
        call. = FALSE
      )
    }
  }
  if (!is.numeric(bound_digits) || length(bound_digits) != 1 ||
      is.na(bound_digits) || !is.finite(bound_digits) ||
      bound_digits < 0 || bound_digits != as.integer(bound_digits)) {
    stop("'bound_digits' must be a nonnegative integer.", call. = FALSE)
  }

  if (is.null(gamma_seq)) {
    gamma_seq <- seq(
      max(0, min(var_x_gamma_const, na.rm = TRUE) - 0.05),
      max(var_x_gamma_const, na.rm = TRUE) + 0.05,
      length.out = 100
    )
  }
  if (is.null(delta_seq)) {
    pad <- diff(range(var_x_delta, na.rm = TRUE)) * 0.1
    if (!is.finite(pad) || pad == 0) {
      pad <- 1
    }
    delta_seq <- seq(
      min(var_x_delta, na.rm = TRUE) - pad,
      max(var_x_delta, na.rm = TRUE) + pad,
      length.out = 100
    )
  }

  point_labels <- as.character(var_x_name)
  if (!is.null(var_x_bound)) {
    has_bound <- is.finite(var_x_bound)
    point_labels[has_bound] <- paste0(
      point_labels[has_bound],
      " (",
      formatC(
        var_x_bound[has_bound],
        format = "f",
        digits = as.integer(bound_digits)
      ),
      ")"
    )
  }

  out <- tibble::tibble(
    var_name = point_labels,
    gamma_const = var_x_gamma_const,
    delta = var_x_delta,
    bound_deterministic = if (is.null(var_x_bound)) NA_real_ else var_x_bound
  )
  sens_df <- dplyr::mutate(
    as.data.frame(expand.grid(gamma_seq, delta_seq)),
    bound = abs(Var1 * Var2) + simplex_penalty
  )

  p <- ggplot2::ggplot(sens_df, ggplot2::aes(x = Var1, y = Var2)) +
    ggplot2::geom_contour(ggplot2::aes(z = bound), color = "black", alpha = 0.55) +
    metR::geom_text_contour(
      ggplot2::aes(z = bound),
      color = "black",
      alpha = 0.7,
      size = 3,
      skip = 0,
      stroke = 0.15
    ) +
    ggplot2::geom_point(
      data = out,
      ggplot2::aes(x = gamma_const, y = delta),
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      x = expression(hat(gamma)[const] ~ "(Simplex weight)"),
      y = expression(hat(delta) ~ "(Imbalance)"),
      title = title
    )

  if (isTRUE(repel)) {
    p <- p + ggrepel::geom_text_repel(
      data = out,
      ggplot2::aes(x = gamma_const, y = delta, label = var_name),
      inherit.aes = FALSE,
      size = text_size,
      max.overlaps = Inf
    )
  } else {
    p <- p + ggplot2::geom_text(
      data = out,
      ggplot2::aes(x = gamma_const, y = delta, label = var_name),
      inherit.aes = FALSE,
      hjust = 0,
      nudge_x = 0.02,
      size = text_size
    )
  }

  p
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("bound", "gamma_const"))
}
