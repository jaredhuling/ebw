# Continuous independence weights: self-consistency at training points and
# out-of-sample recovery of the true GPS density ratio at genuinely new points.

test_that("continuous IW self-consistency and out-of-sample GPS recovery", {
  set.seed(20260809L)
  n <- 400; p <- 4
  alpha <- rep(sqrt(0.25 / p), p)      # ||alpha||^2 = 0.25 -> sigma^2 = 1.25
  X <- matrix(rnorm(n * p), n, p)
  A <- drop(X %*% alpha) + rnorm(n)
  sigma2 <- sum(alpha^2) + 1

  fit <- energy_balance(x = X, treatment = A, treatment_type = "continuous",
                        lambda = sqrt(n), scaling = "none")
  expect_s3_class(fit, "ebw")
  expect_equal(fit$treatment_type, "continuous")

  # (1) self-consistency: predict at the training (X_i, A_i) reproduces w
  w_pred <- predict(fit, newdata = list(x = X, a = A))
  expect_lt(max(abs(w_pred - fit$weights)), 1e-4)

  # closed-form GPS ratio w*(x,a) = f_A(a) / f_{A|X}(a|x)
  gps_ratio <- function(x, a) {
    mu <- drop(x %*% alpha)
    dnorm(a, mean = 0, sd = sqrt(sigma2)) / dnorm(a, mean = mu, sd = 1)
  }

  # (2) OUT-OF-SAMPLE recovery on genuinely new (x, a) draws via predict()
  set.seed(20260810L)
  Xh <- matrix(rnorm(2000 * p), 2000, p)
  Ah <- drop(Xh %*% alpha) + rnorm(2000)
  w_h <- predict(fit, newdata = list(x = Xh, a = Ah))
  gps_h <- gps_ratio(Xh, Ah)
  pos <- w_h > 0
  rho <- cor(w_h[pos], gps_h[pos], method = "spearman")
  expect_gt(rho, 0.8)
})
