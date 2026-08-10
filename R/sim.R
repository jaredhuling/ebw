# Self-contained data simulators used in examples and tests. Base R only, so
# they carry no hard dependency; MASS is not required.

# Draw n rows of a mean-zero multivariate normal with AR(1) covariance
# rho^|i-j|, using a Cholesky factor (no MASS).
#' @keywords internal
#' @noRd
.rmvnorm_ar1 <- function(n, p, rho = 0) {
  Sig <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  L <- chol(Sig)                       # upper triangular, Sig = t(L) %*% L
  matrix(stats::rnorm(n * p), n, p) %*% L
}

#' Simulate confounded data with a binary treatment
#'
#' Generates covariates, a confounded binary treatment, and a continuous
#' outcome for illustrating and testing energy balancing weights. The treatment
#' assignment and the outcome both depend on the covariates, so naive
#' comparisons are confounded.
#'
#' @param n.obs number of observations.
#' @param n.vars number of covariates (at least 4).
#' @param AR.cor AR(1) correlation of the covariates (absolute value below 1).
#' @param beta coefficient vector for the propensity (log-odds) model; recycled
#'   or truncated to `n.vars`. Defaults to a decaying alternating-sign pattern.
#' @param outcome.coef coefficient vector for the outcome model.
#' @param tau true (constant) average treatment effect.
#' @param sd.error outcome noise standard deviation.
#'
#' @return A list with covariate matrix `x`, treatment `trt`, outcome `y`,
#'   propensity `prob.trt`, and the true treatment effect `trt.eff`.
#'
#' @examples
#' set.seed(20260809)
#' dat <- sim_confounded_data(n.obs = 200, n.vars = 6, AR.cor = 0.3)
#' str(dat)
#'
#' @export
sim_confounded_data <- function(n.obs, n.vars = 10, AR.cor = 0,
                                beta = NULL, outcome.coef = NULL,
                                tau = 1, sd.error = 1) {
  if (n.vars < 4) stop("'n.vars' must be 4 or greater")
  stopifnot(abs(AR.cor) < 1)
  x <- .rmvnorm_ar1(n.obs, n.vars, AR.cor)
  if (is.null(beta)) {
    beta <- (-1)^(seq_len(n.vars) + 1) / sqrt(seq_len(n.vars))
  } else {
    beta <- rep_len(beta, n.vars)
  }
  if (is.null(outcome.coef)) {
    outcome.coef <- rep_len(c(1, -1, 0.5, 0.25), n.vars)
  } else {
    outcome.coef <- rep_len(outcome.coef, n.vars)
  }
  lp <- drop(x %*% beta)
  prob_trt <- 1 / (1 + exp(-lp))
  trt <- stats::rbinom(n.obs, 1, prob_trt)
  mu <- drop(x %*% outcome.coef)
  y <- mu + tau * trt + stats::rnorm(n.obs, sd = sd.error)
  list(x = x, trt = trt, y = y, prob.trt = prob_trt, trt.eff = tau)
}

#' Simulate confounded data with a continuous treatment
#'
#' Generates covariates and a continuous treatment that depends linearly on
#' them (a linear-Gaussian design), for illustrating and testing independence
#' weights. With standard normal covariates and treatment noise variance 1, the
#' treatment marginal variance is `1 + alpha.norm2`.
#'
#' @param n.obs number of observations.
#' @param n.vars number of covariates.
#' @param alpha.norm2 squared Euclidean norm of the covariate-to-treatment
#'   coefficient vector (confounding strength). The coefficients are set to
#'   `sqrt(alpha.norm2 / n.vars)` in every coordinate.
#' @param AR.cor AR(1) correlation of the covariates.
#' @param outcome.coef coefficient vector for the outcome model.
#' @param tau slope of the (linear) treatment effect on the outcome.
#' @param sd.error outcome noise standard deviation.
#'
#' @return A list with covariate matrix `x`, treatment `a`, outcome `y`, the
#'   coefficient vector `alpha`, and the treatment marginal variance
#'   `sigma2.A`.
#'
#' @examples
#' set.seed(20260809)
#' dat <- sim_continuous_data(n.obs = 200, n.vars = 4, alpha.norm2 = 0.25)
#' str(dat)
#'
#' @export
sim_continuous_data <- function(n.obs, n.vars = 4, alpha.norm2 = 0.25,
                                AR.cor = 0, outcome.coef = NULL,
                                tau = 1, sd.error = 1) {
  stopifnot(abs(AR.cor) < 1, alpha.norm2 >= 0)
  x <- .rmvnorm_ar1(n.obs, n.vars, AR.cor)
  alpha <- rep(sqrt(alpha.norm2 / n.vars), n.vars)
  a <- drop(x %*% alpha) + stats::rnorm(n.obs)
  if (is.null(outcome.coef)) {
    outcome.coef <- rep_len(c(1, -1, 0.5, 0.25), n.vars)
  } else {
    outcome.coef <- rep_len(outcome.coef, n.vars)
  }
  mu <- drop(x %*% outcome.coef)
  y <- mu + tau * a + stats::rnorm(n.obs, sd = sd.error)
  list(x = x, a = a, y = y, alpha = alpha, sigma2.A = 1 + alpha.norm2)
}
