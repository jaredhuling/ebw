# Constrained independence weights (preserve_means / decorrelate_moments):
# self-consistency of the augmented out-of-sample predictor, satisfaction of
# the moment constraints on the training sample, and out-of-sample recovery of
# the true GPS density ratio.

make_iw_data <- function(n, p, seed) {
  set.seed(seed)
  alpha <- rep(sqrt(0.25 / p), p)          # ||alpha||^2 = 0.25 -> sigma^2 = 1.25
  X <- matrix(rnorm(n * p), n, p)
  A <- drop(X %*% alpha) + rnorm(n)
  list(X = X, A = A, alpha = alpha, sigma2 = sum(alpha^2) + 1)
}

test_that("constrained IW: self-consistency, constraints, and OOS GPS recovery", {
  d <- make_iw_data(400, 4, 20260809L)
  X <- d$X; A <- d$A; alpha <- d$alpha; sigma2 <- d$sigma2

  set.seed(20260810L)
  Xh <- matrix(rnorm(2000 * 4), 2000, 4)
  Ah <- drop(Xh %*% alpha) + rnorm(2000)
  gps_ratio <- function(x, a)
    dnorm(a, 0, sqrt(sigma2)) / dnorm(a, drop(x %*% alpha), 1)
  gps_h <- gps_ratio(Xh, Ah)

  configs <- list(
    list(pm = TRUE,  dc = FALSE, lab = "preserve_means"),
    list(pm = FALSE, dc = TRUE,  lab = "decorrelate_moments"),
    list(pm = TRUE,  dc = TRUE,  lab = "both"))

  for (cfg in configs) {
    fit <- independence_weights(A = A, X = X, lambda = sqrt(nrow(X)),
                                preserve_means = cfg$pm,
                                decorrelate_moments = cfg$dc,
                                scaling = "none")

    # (i) self-consistency: predict at training points reproduces QP weights
    w_pred <- predict(fit, newdata = list(x = X, a = A))
    sc <- max(abs(w_pred - fit$weights))
    expect_lt(sc, 1e-6)

    # (ii) preserve_means: weighted mean equals unweighted mean
    if (cfg$pm) {
      w <- fit$weights
      csat <- max(abs(colSums(w * X) / sum(w) - colMeans(X)))
      expect_lt(csat, 1e-8)
    }

    # (iii) out-of-sample recovery of the GPS density ratio
    w_h <- predict(fit, newdata = list(x = Xh, a = Ah))
    pos <- w_h > 0
    rho <- cor(w_h[pos], gps_h[pos], method = "spearman")
    expect_gt(rho, 0.8)
  }
})
