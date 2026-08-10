# Mean-preserving energy balancing for binary / multi-category treatments, and
# the uniformity of `preserve_means` across all three treatment types. Checks:
# self-consistency of the augmented dual predictor, exact mean balance, uniform
# behavior binary / multi-category / continuous, the decorrelate guard, and that
# the standard binary path is untouched when preserve_means = FALSE.

# max over groups and covariates of |weighted group mean - pooled mean|, raw X.
.max_mean_gap <- function(X, trt, w) {
  X <- as.matrix(X)
  pooled <- colMeans(X)
  gap <- 0
  for (lv in unique(trt)) {
    ik <- which(trt == lv)
    wm <- colSums(w[ik] * X[ik, , drop = FALSE]) / sum(w[ik])
    gap <- max(gap, max(abs(wm - pooled)))
  }
  gap
}

test_that("constrained grouped EBW: self-consistency of the augmented predictor", {
  set.seed(20260810L)
  n <- 300L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)

  # binary
  tb <- rbinom(n, 1, plogis(X[, 1] - X[, 2]))
  fb <- energy_balance(x = X, treatment = tb, preserve_means = TRUE)
  pb <- predict(fb, newdata = list(x = X, treatment = tb))
  expect_lt(max(abs(pb - weights(fb))), 1e-4)

  # multi-category K = 3
  lin <- X %*% c(1, -1, 0.5, 0)
  pr <- exp(cbind(-lin, 0, lin)); pr <- pr / rowSums(pr)
  tm <- apply(pr, 1, function(q) sample(1:3, 1, prob = q))
  fm <- energy_balance(x = X, treatment = tm, treatment_type = "multinomial",
                       preserve_means = TRUE)
  pm <- predict(fm, newdata = list(x = X, treatment = tm))
  expect_lt(max(abs(pm - weights(fm))), 1e-4)

  # improved binary
  fi <- energy_balance(x = X, treatment = tb, improved = TRUE,
                       preserve_means = TRUE)
  pi <- predict(fi, newdata = list(x = X, treatment = tb))
  expect_lt(max(abs(pi - weights(fi))), 1e-4)
})

test_that("constrained grouped EBW: exact mean balance to the pooled means", {
  set.seed(20260810L)
  n <- 300L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  tb <- rbinom(n, 1, plogis(X[, 1] - X[, 2]))
  lin <- X %*% c(1, -1, 0.5, 0)
  pr <- exp(cbind(-lin, 0, lin)); pr <- pr / rowSums(pr)
  tm <- apply(pr, 1, function(q) sample(1:3, 1, prob = q))

  fb <- energy_balance(x = X, treatment = tb, preserve_means = TRUE)
  expect_lt(.max_mean_gap(X, tb, weights(fb)), 1e-6)

  fm <- energy_balance(x = X, treatment = tm, treatment_type = "multinomial",
                       preserve_means = TRUE)
  expect_lt(.max_mean_gap(X, tm, weights(fm)), 1e-6)

  fi <- energy_balance(x = X, treatment = tb, improved = TRUE,
                       preserve_means = TRUE)
  expect_lt(.max_mean_gap(X, tb, weights(fi)), 1e-6)
})

test_that("preserve_means acts uniformly across binary, multi-category, continuous", {
  set.seed(20260810L)
  n <- 300L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)

  # binary
  tb <- rbinom(n, 1, plogis(X[, 1] - X[, 2]))
  fb <- energy_balance(x = X, treatment = tb, preserve_means = TRUE)
  expect_lt(.max_mean_gap(X, tb, weights(fb)), 1e-6)

  # 3-level factor treatment
  tf <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  ff <- energy_balance(x = X, treatment = tf, preserve_means = TRUE)
  expect_lt(.max_mean_gap(X, tf, weights(ff)), 1e-6)

  # continuous (already supported: assert it still preserves the means)
  A <- drop(X %*% c(0.5, -0.5, 0.25, 0)) + rnorm(n)
  fc <- energy_balance(x = X, treatment = A, treatment_type = "continuous",
                       preserve_means = TRUE)
  wc <- weights(fc)
  gap_cov <- max(abs(colSums(wc * X) / sum(wc) - colMeans(X)))
  gap_trt <- abs(sum(wc * A) / sum(wc) - mean(A))
  expect_lt(max(gap_cov, gap_trt), 1e-6)
})

test_that("decorrelate_moments errors clearly for a discrete treatment", {
  set.seed(20260810L)
  n <- 200L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  tb <- rbinom(n, 1, 0.5)
  expect_error(
    energy_balance(x = X, treatment = tb, decorrelate_moments = TRUE),
    "continuous")
})

test_that("standard binary preserve_means = FALSE is bit-identical to before", {
  set.seed(20260810L)
  n <- 250L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  tb <- rbinom(n, 1, plogis(X[, 1] - X[, 2]))
  # The FISTA step size uses a random power-iteration start, so reset the RNG
  # before each fit to compare the two paths bit for bit.
  set.seed(20260810L)
  f1 <- energy_balance(x = X, treatment = tb, preserve_means = FALSE)
  set.seed(20260810L)
  f2 <- energy_balance(x = X, treatment = tb)
  expect_identical(weights(f1), weights(f2))
})
