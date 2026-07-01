test_that("vertreg uses pseudo-inverse for singular designs", {
  data_pre <- data.frame(
    Y = c(5, 10, 15),
    X1 = c(1, 2, 3),
    X2 = c(2, 4, 6)
  )
  data_post <- data.frame(Y = 30, X1 = 4, X2 = 8)

  expect_equal(
    vertreg(
      formula = "Y ~ -1 + X1 + X2",
      data_pre = data_pre,
      data_post = data_post,
      pseudo_inverse = TRUE
    ),
    10,
    tolerance = 1e-12
  )
})

test_that("estimate_params_partial returns known pseudo-inverse values", {
  data_pre <- data.frame(
    Y = c(1.25, -1, 2, 0, 2.25),
    Z = c(0, 0, 1, 0, 1),
    X1 = c(1, 0, 0, 0, 1),
    X2 = c(0, 1, 0, 0, 1)
  )
  data_post <- data.frame(Y = 0, Z = 4, X1 = 2, X2 = 1)

  res <- estimate_params_partial(
    fm_z_on_x = "Z ~ -1 + X1 + X2",
    fm_y_on_z_and_x = "Y ~ -1 + Z + X1 + X2",
    data_pre = data_pre,
    data_post = data_post,
    pseudo_inverse = TRUE
  )

  expect_equal(res$gamma, 2, tolerance = 1e-12)
  expect_equal(res$delta, 3, tolerance = 1e-12)
  expect_equal(res$bias, 6, tolerance = 1e-12)
})

test_that("estimate_params_partial validates formula columns in pre and post data", {
  data_pre <- data.frame(Y = 1:4, X1 = 2:5)
  data_post <- data.frame(Y = 5, X1 = 6)

  expect_error(
    estimate_params_partial(
      fm_z_on_x = "Z ~ -1 + X1",
      fm_y_on_z_and_x = "Y ~ -1 + Z + X1",
      data_pre = data_pre,
      data_post = data_post
    ),
    "missing from 'data_pre'"
  )

  expect_error(
    estimate_params_partial(
      fm_z_on_x = "Z ~ -1 + X1",
      fm_y_on_z_and_x = "Y ~ -1 + Z + X1",
      data_pre = transform(data_pre, Z = 1:4),
      data_post = data_post
    ),
    "missing from 'data_post'"
  )
})

test_that("estimate_params_partial_multi returns known component decomposition", {
  data_pre <- data.frame(
    Y = c(1.25, -1, 2, -0.5, 1.75),
    Z1 = c(0, 0, 1, 0, 1),
    Z2 = c(0, 0, 0, 1, 1),
    X1 = c(1, 0, 0, 0, 1),
    X2 = c(0, 1, 0, 0, 1)
  )
  data_post <- data.frame(
    Y = 0,
    Z1 = 4,
    Z2 = -1,
    X1 = 2,
    X2 = 1
  )

  res <- estimate_params_partial_multi(
    var_z_name = c("Z1", "Z2"),
    var_x_name = c("X1", "X2"),
    var_y_name = "Y",
    data_pre = data_pre,
    data_post = data_post,
    pseudo_inverse = TRUE
  )

  expect_equal(res$gamma, 1)
  expect_named(res$component_gamma, c("Z1", "Z2"))
  expect_named(res$component_delta, c("Z1", "Z2"))
  expect_named(res$component_bias, c("Z1", "Z2"))
  expect_equal(res$component_gamma, c(Z1 = 2, Z2 = -0.5), tolerance = 1e-12)
  expect_equal(res$component_delta, c(Z1 = 3, Z2 = -2), tolerance = 1e-12)
  expect_equal(res$component_bias, c(Z1 = 6, Z2 = 1), tolerance = 1e-12)
  expect_equal(res$imbalance, 7, tolerance = 1e-12)
  expect_equal(res$bias, 7, tolerance = 1e-12)
})
