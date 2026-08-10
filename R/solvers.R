# Optimization primitives shared by the binary and continuous engines:
# the exact Euclidean projection onto a scaled simplex and a FISTA loop.

#' Euclidean projection onto a scaled simplex
#'
#' Projects a vector onto \eqn{\{w \ge 0, \sum_i w_i = s\}} in Euclidean
#' distance, using the exact sort-based algorithm.
#'
#' @param v numeric vector to project.
#' @param s target sum of the projected vector (default 1).
#'
#' @return The projection of `v` onto the simplex of scale `s`.
#'
#' @keywords internal
#' @noRd
proj_simplex <- function(v, s = 1) {
  n <- length(v)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u) - s
  ind <- seq_len(n)
  cond <- u - cssv / ind > 0
  rho <- max(which(cond))
  theta <- cssv[rho] / rho
  pmax(v - theta, 0)
}

#' FISTA loop for a simplex-constrained quadratic
#'
#' Minimizes \eqn{w' Q w + a' w} over the scaled simplex
#' \eqn{\{w \ge 0, \sum w = s\}} by accelerated projected gradient (FISTA).
#' The quadratic is supplied implicitly through `Qmv`, which must return
#' \eqn{Q w} (that is, half the gradient of the quadratic part).
#'
#' @param Qmv function returning the matrix-vector product `Q %*% w`.
#' @param a_vec linear coefficient vector.
#' @param n_w length of the weight vector.
#' @param s simplex scale (sum of the weights).
#' @param w0 optional starting point; defaults to the uniform vector `s / n_w`.
#' @param project_momentum logical; if `TRUE`, clip the momentum iterate at 0
#'   before the gradient step (matches the continuous-engine convention).
#' @param max_iter maximum number of iterations.
#' @param tol relative objective tolerance, checked every 100 iterations.
#'
#' @return A list with the solution `w` and the Lipschitz estimate `L`.
#'
#' @keywords internal
#' @noRd
.fista_simplex <- function(Qmv, a_vec, n_w, s = 1, w0 = NULL,
                           project_momentum = FALSE,
                           max_iter = 20000L, tol = 1e-11) {
  # Lipschitz constant of the gradient 2 Q via power iteration.
  v <- stats::rnorm(n_w); v <- v / sqrt(sum(v^2))
  for (i in 1:100) { Qv <- Qmv(v); nv <- sqrt(sum(Qv^2)); if (nv == 0) break; v <- Qv / nv }
  L_est <- 2 * sum(v * Qmv(v)); step <- 0.99 / L_est

  w <- if (is.null(w0)) rep(s / n_w, n_w) else w0
  w_prev <- w; tk <- 1; obj_prev <- Inf
  for (it in seq_len(max_iter)) {
    mom <- (tk - 1) / (tk + 2)
    y <- w + mom * (w - w_prev)
    if (project_momentum) y <- pmax(y, 0)
    grad <- 2 * Qmv(y) + a_vec
    w_new <- proj_simplex(y - step * grad, s)
    if (it %% 100 == 0) {
      obj <- sum(w_new * Qmv(w_new)) + sum(a_vec * w_new)
      if (abs(obj - obj_prev) / (abs(obj) + 1e-14) < tol) { w <- w_new; break }
      obj_prev <- obj
    }
    w_prev <- w; w <- w_new; tk <- (1 + sqrt(1 + 4 * tk^2)) / 2
  }
  list(w = w, L = L_est)
}
