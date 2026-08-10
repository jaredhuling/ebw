# Formula vs matrix interface equivalence and backward-compatible arguments.

test_that("formula and matrix interfaces give identical weights", {
  set.seed(20260809L)
  n <- 300; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  trt <- rbinom(n, 1, plogis(X[, 1] - X[, 2] + 0.5 * X[, 3]))
  df <- as.data.frame(X)
  names(df) <- paste0("V", seq_len(p))
  df$trt <- trt

  # same seed before each fit so the FISTA power-iteration draws match
  set.seed(11L)
  fit_m <- energy_balance(x = X, treatment = trt, estimand = "ATE")
  set.seed(11L)
  fit_f <- energy_balance(trt ~ ., data = df, estimand = "ATE")

  expect_equal(unname(fit_m$weights), unname(fit_f$weights))
})

test_that("old-style energy_balance(trt=, x=, method=) still works", {
  set.seed(20260809L)
  n <- 200; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  trt <- rbinom(n, 1, plogis(X[, 1] - X[, 2]))

  fit_ate <- energy_balance(trt = trt, x = X, method = "ATE", standardize = TRUE)
  expect_s3_class(fit_ate, "ebw")
  expect_equal(fit_ate$treatment_type, "binary")
  expect_equal(fit_ate$estimand, "ATE")
  expect_length(fit_ate$weights, n)

  fit_att <- energy_balance(trt = trt, x = X, method = "ATT", standardize = FALSE)
  expect_s3_class(fit_att, "ebw")
  expect_equal(fit_att$estimand, "ATT")
  expect_equal(fit_att$scaling$method, "none")
})
