# balance(): runs for binary and continuous, in sample and out of sample, and
# returns the documented structure with a sensible effective sample size.

test_that("balance() works for a binary treatment, in and out of sample", {
  set.seed(20260809L)
  dat <- sim_confounded_data(n.obs = 300, n.vars = 6)
  fit <- energy_balance(x = dat$x, treatment = dat$trt, estimand = "ATE")

  b <- balance(fit)
  expect_s3_class(b, "balance.ebw")
  expect_equal(b$treatment_type, "binary")
  expect_true(is.data.frame(b$per_covariate))
  expect_named(b$per_covariate, c("variable", "unadjusted", "adjusted"))
  expect_equal(nrow(b$per_covariate), ncol(dat$x))
  expect_true(all(c("unadjusted", "adjusted") %in% names(b$energy_distance)))
  expect_gt(b$ess, 0)
  expect_lte(b$ess, b$n + 1e-8)
  # weighting should not worsen aggregate imbalance
  expect_lte(b$energy_distance[["adjusted"]], b$energy_distance[["unadjusted"]] + 1e-8)

  # out of sample on a fresh draw
  set.seed(20260811L)
  nd <- sim_confounded_data(n.obs = 200, n.vars = 6)
  bo <- balance(fit, newdata = list(x = nd$x, trt = nd$trt))
  expect_s3_class(bo, "balance.ebw")
  expect_equal(bo$sample, "out-of-sample")
  expect_gt(bo$ess, 0)
  expect_lte(bo$ess, bo$n + 1e-8)
})

test_that("balance() works for a continuous treatment, in and out of sample", {
  set.seed(20260809L)
  dat <- sim_continuous_data(n.obs = 300, n.vars = 4, alpha.norm2 = 0.5)
  fit <- independence_weights(A = dat$a, X = dat$x)

  b <- balance(fit)
  expect_s3_class(b, "balance.ebw")
  expect_equal(b$treatment_type, "continuous")
  expect_named(b$per_covariate, c("variable", "unadjusted", "adjusted"))
  expect_true(all(c("unadjusted", "adjusted") %in% names(b$dcov)))
  expect_true(is.data.frame(b$means))
  expect_gt(b$ess, 0)
  expect_lte(b$ess, b$n + 1e-8)
  # weighting should reduce the overall dependence between X and A
  expect_lt(b$dcov[["adjusted"]], b$dcov[["unadjusted"]])

  set.seed(20260811L)
  nd <- sim_continuous_data(n.obs = 200, n.vars = 4, alpha.norm2 = 0.5)
  bo <- balance(fit, newdata = list(x = nd$x, a = nd$a))
  expect_s3_class(bo, "balance.ebw")
  expect_equal(bo$sample, "out-of-sample")
  expect_gt(bo$ess, 0)
  expect_lte(bo$ess, bo$n + 1e-8)
})
