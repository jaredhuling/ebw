# Binary energy balancing: dual-representation self-consistency and sanity.

test_that("predict at training treated rows reproduces the fitted weights", {
  set.seed(20260809L)
  n <- 400; p <- 6
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(1, -1, 0.5, -1, 0, 0.8)
  pscore <- plogis(drop(X %*% beta))
  trt <- rbinom(n, 1, pscore)

  fit <- energy_balance(x = X, treatment = trt, estimand = "ATE")
  expect_s3_class(fit, "ebw")
  expect_equal(fit$treatment_type, "binary")

  # fitted weights for the treated (primary) group
  treated_level <- fit$engine$primary_level
  idx_t <- which(trt == treated_level)
  w_t <- fit$weights[idx_t]

  # predict at the training treated rows -> should reproduce fitted weights
  pred <- predict(fit, newdata = X[idx_t, , drop = FALSE])
  expect_lt(max(abs(pred - w_t)), 1e-3)

  # weights nonnegative
  expect_true(all(fit$weights >= -1e-10))

  # ESS sensible: strictly positive and no larger than n
  w <- fit$weights
  ess <- sum(w)^2 / sum(w^2)
  expect_gt(ess, 0)
  expect_lte(ess, n + 1e-8)
})

test_that("predict() with no newdata returns the fitted weights", {
  set.seed(20260809L)
  n <- 200; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  trt <- rbinom(n, 1, plogis(X[, 1] - X[, 2]))
  fit <- energy_balance(x = X, treatment = trt)
  expect_identical(predict(fit), weights(fit))
})
