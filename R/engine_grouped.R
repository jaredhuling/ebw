# Internal engine for grouped energy balancing weights (EBW), the unified dual
# representation that covers standard binary EBW (K = 2, beta = 0), improved
# binary EBW (K = 2, beta = 1), and multi-category EBW (K >= 2). Ported from the
# validated `grouped_ebw_fit` / `predict_grouped` prototype. Expects an
# already-scaled covariate matrix; scaling is handled by the caller so predict
# can reuse it.
#
# The problem places a per-group simplex on the weights (each group's weights
# sum to one) and minimizes, over the product simplex,
#   sum_k alpha_k * ED^2(mu_{k,w}, mu_0)
#     + sum_{j < k} beta_{jk} * ED^2(mu_{j,w}, mu_{k,w}) + eps ||w||^2,
# whose stationarity conditions give the dual weight representation
#   w^{(k)}_i = ( (nu_k - 2 h_k(X_i)) / (2 eps) )_+,
#   h_k(x) = alpha_k (mu_{k,w}(x) - mu_0(x))
#            + sum_{j != k} beta_{jk} (mu_{k,w}(x) - mu_{j,w}(x)),
# with mu_0 the pooled (full-sample) target embedding. This extends to new
# covariate values, so a fitted grouped weighting predicts out of sample.

#' Fit grouped energy balancing weights
#'
#' @param X scaled numeric covariate matrix (n x p).
#' @param groups integer vector in `1..K` giving each unit's group.
#' @param alpha length-K vector of per-group target weights (default all one),
#'   matching each group to the pooled target.
#' @param beta K x K symmetric matrix of cross-group coupling weights; zero gives
#'   standard (uncoupled) balancing, one adds the between-group energy distance
#'   (the improved variant).
#' @param eps ridge parameter.
#' @param kernel covariate kernel function.
#' @param max_iter,tol FISTA controls.
#'
#' @return A list carrying the fitted length-n weights `w` and the
#'   dual-representation ingredients needed by the predictor.
#'
#' @keywords internal
#' @noRd
.ebw_grouped <- function(X, groups, alpha = NULL, beta = NULL, eps = 1e-2,
                         kernel = energy_kernel_X,
                         max_iter = 30000L, tol = 1e-12) {
  X <- as.matrix(X)
  n <- nrow(X)
  groups <- as.integer(groups)
  K <- max(groups)
  if (is.null(alpha)) alpha <- rep(1, K)
  if (is.null(beta)) beta <- matrix(0, K, K)

  Kmat <- kernel(X, X)                       # n x n Gram
  tvec <- rep(1 / n, n)                       # pooled target measure
  idxs <- lapply(seq_len(K), function(k) which(groups == k))

  # a_k: length-n vector, equal to w on group k, zero elsewhere
  avec <- function(w, k) { z <- numeric(n); z[idxs[[k]]] <- w[idxs[[k]]]; z }

  # gradient: grad_i = 2 h_{g(i)}(X_i) + 2 eps w_i  (only defined on each group)
  grad_fn <- function(w) {
    A_k <- lapply(seq_len(K), function(k) avec(w, k))
    g <- numeric(n)
    for (k in seq_len(K)) {
      resid <- alpha[k] * (A_k[[k]] - tvec)
      for (j in setdiff(seq_len(K), k))
        resid <- resid + beta[j, k] * (A_k[[k]] - A_k[[j]])
      hk <- drop(Kmat %*% resid)             # h_k evaluated at all points
      ii <- idxs[[k]]; g[ii] <- 2 * hk[ii] + 2 * eps * w[ii]
    }
    g
  }
  obj_fn <- function(w) {
    A_k <- lapply(seq_len(K), function(k) avec(w, k)); val <- eps * sum(w^2)
    for (k in seq_len(K)) {
      d <- A_k[[k]] - tvec; val <- val + alpha[k] * drop(d %*% Kmat %*% d)
    }
    for (k in seq_len(K)) for (j in seq_len(K)) if (j < k) {
      d <- A_k[[j]] - A_k[[k]]; val <- val + beta[j, k] * drop(d %*% Kmat %*% d)
    }
    val
  }
  proj_all <- function(w) {
    for (k in seq_len(K)) w[idxs[[k]]] <- proj_simplex(w[idxs[[k]]], 1)
    w
  }

  # step size via power iteration on the Hessian action (grad is affine)
  g0 <- grad_fn(numeric(n))                  # gradient at 0 = linear part
  Hmv <- function(v) grad_fn(v) - g0         # Hessian * v
  v <- stats::rnorm(n); v <- v / sqrt(sum(v^2))
  for (i in 1:100) { Hv <- Hmv(v); nv <- sqrt(sum(Hv^2)); if (nv == 0) break; v <- Hv / nv }
  L <- sum(v * Hmv(v)); step <- 0.9 / L

  w <- numeric(n)
  for (k in seq_len(K)) w[idxs[[k]]] <- 1 / length(idxs[[k]])
  w_prev <- w; tk <- 1; obj_prev <- Inf
  for (it in seq_len(max_iter)) {
    mom <- (tk - 1) / (tk + 2)
    y <- w + mom * (w - w_prev); y[y < 0] <- 0
    w_new <- proj_all(y - step * grad_fn(y))
    if (it %% 200 == 0) {
      ob <- obj_fn(w_new)
      if (abs(ob - obj_prev) / (abs(ob) + 1e-14) < tol) { w <- w_new; break }
      obj_prev <- ob
    }
    w_prev <- w; w <- w_new; tk <- (1 + sqrt(1 + 4 * tk^2)) / 2
  }

  # recover per-group multipliers nu_k from active-set stationarity
  gr <- grad_fn(w)                            # = 2 h_k + 2 eps w on each group
  nu <- numeric(K)
  for (k in seq_len(K)) {
    ii <- idxs[[k]]; act <- ii[w[ii] > 1e-8 * max(w[ii])]
    nu[k] <- if (length(act)) stats::median(gr[act]) else NA_real_
  }

  list(kind = "grouped", X = X, groups = groups, w = w, eps = eps,
       alpha = alpha, beta = beta, nu = nu, n = n, K = K, idxs = idxs,
       kernel = kernel)
}

#' Predict grouped EBW weights out of sample
#'
#' Evaluates the dual-representation weight predictor at new covariate rows,
#' each carrying its group index. Does not refit.
#'
#' @param fit the list returned by `.ebw_grouped`.
#' @param Xnew already-scaled new covariate matrix (m x p).
#' @param gnew integer vector in `1..K` giving each new unit's group.
#' @param type `"weights"` for predicted weights or `"balancing"` for the
#'   balancing function \eqn{h_k(x)}.
#'
#' @return A numeric vector of predicted weights (or balancing-function values).
#'
#' @keywords internal
#' @noRd
.predict_grouped <- function(fit, Xnew, gnew, type = "weights") {
  Xnew <- as.matrix(Xnew)
  gnew <- as.integer(gnew)
  m <- nrow(Xnew); K <- fit$K
  KxX <- fit$kernel(Xnew, fit$X)             # m x n
  mu0 <- rowSums(KxX) / fit$n                 # mu_0(x)
  muk <- sapply(seq_len(K), function(k) {    # mu_{k,w}(x), m x K
    z <- numeric(fit$n); z[fit$idxs[[k]]] <- fit$w[fit$idxs[[k]]]
    drop(KxX %*% z)
  })
  if (m == 1L) muk <- matrix(muk, nrow = 1L)
  h <- numeric(m)
  for (i in seq_len(m)) {
    k <- gnew[i]
    hk <- fit$alpha[k] * (muk[i, k] - mu0[i])
    for (j in setdiff(seq_len(K), k))
      hk <- hk + fit$beta[j, k] * (muk[i, k] - muk[i, j])
    h[i] <- hk
  }
  if (identical(type, "balancing")) return(h)
  out <- numeric(m)
  for (i in seq_len(m)) out[i] <- max(0, (fit$nu[gnew[i]] - 2 * h[i]) / (2 * fit$eps))
  out
}
