# plot.ebw(): returns a ggplot object for the Love plot and the weight
# distribution, for both binary and continuous treatments.

test_that("plot.ebw returns ggplot objects (binary)", {
  skip_if_not_installed("ggplot2")
  set.seed(20260809L)
  dat <- sim_confounded_data(n.obs = 200, n.vars = 6)
  fit <- energy_balance(x = dat$x, treatment = dat$trt)

  p_love <- plot(fit, type = "balance")
  expect_s3_class(p_love, "ggplot")

  p_w <- plot(fit, type = "weights")
  expect_s3_class(p_w, "ggplot")
})

test_that("plot.ebw returns ggplot objects (multi-category K=3)", {
  skip_if_not_installed("ggplot2")
  set.seed(20260810L)
  n <- 240; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  tm <- sample(1:3, n, replace = TRUE)
  fit <- energy_balance(x = X, treatment = tm)

  p_love <- plot(fit, type = "balance")
  expect_s3_class(p_love, "ggplot")

  p_w <- plot(fit, type = "weights")
  expect_s3_class(p_w, "ggplot")
})

test_that("plot.ebw returns ggplot objects (continuous)", {
  skip_if_not_installed("ggplot2")
  set.seed(20260809L)
  dat <- sim_continuous_data(n.obs = 200, n.vars = 4, alpha.norm2 = 0.5)
  fit <- independence_weights(A = dat$a, X = dat$x)

  p_love <- plot(fit, type = "balance")
  expect_s3_class(p_love, "ggplot")

  p_w <- plot(fit, type = "weights")
  expect_s3_class(p_w, "ggplot")
})
