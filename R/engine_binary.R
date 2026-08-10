# Internal engine for binary energy balancing weights (EBW). Ports ebw_fit
# from the validated dual-representation prototype. Expects an already-scaled
# covariate matrix; scaling is handled by the caller so predict can reuse it.

#' Fit binary energy balancing weights
#'
#' Balances the covariate distribution of one treatment group to the full
#' sample by solving the simplex-constrained ridge QP
#' \eqn{\min_w w'(K_{aa} + \epsilon I) w + a'w} and returns the fitted weights
#' together with the dual-representation ingredients needed for out-of-sample
#' prediction.
#'
#' @param X scaled numeric covariate matrix (n x p).
#' @param trt vector giving the (binary) treatment.
#' @param level the value of `trt` identifying the group to be balanced.
#' @param target_X scaled covariate matrix defining the target empirical
#'   measure the group is balanced toward. Defaults to the full sample `X`,
#'   which is the ATE building block; passing the treated subsample gives the
#'   ATT target. The same formulas are used, only the target measure changes.
#' @param eps ridge parameter.
#' @param kernel covariate kernel function.
#' @param max_iter,tol FISTA controls.
#'
#' @return A list carrying the fitted weights `w`, the group covariate matrix
#'   `Xa`, the target covariate matrix `Xall`, the multiplier `nu`, and the
#'   scalars needed by the predictor.
#'
#' @keywords internal
#' @noRd
.ebw_binary <- function(X, trt, level = 1, target_X = NULL, eps = 1e-2,
                        kernel = energy_kernel_X,
                        max_iter = 20000L, tol = 1e-11) {
  X <- as.matrix(X)
  if (is.null(target_X)) target_X <- X
  target_X <- as.matrix(target_X)
  n <- nrow(target_X)
  idx_a <- which(trt == level)
  Xa <- X[idx_a, , drop = FALSE]
  na <- nrow(Xa)

  Kaa <- kernel(Xa, Xa)                 # na x na, within group a
  KaX <- kernel(Xa, target_X)           # na x n_target, group a vs target
  a_vec <- -2 * rowSums(KaX) / n        # linear term a_i = -(2/n) sum_j k(X_i,T_j)

  Qmv <- function(w) drop(Kaa %*% w) + eps * w
  sol <- .fista_simplex(Qmv, a_vec, n_w = na, s = 1, w0 = rep(1 / na, na),
                        project_momentum = FALSE, max_iter = max_iter, tol = tol)
  w <- sol$w

  # Recover nu from active-set stationarity 2[Q_eps w]_i + a_i = nu, using
  # well-interior active units so boundary units near 0 do not bias it.
  stat <- 2 * Qmv(w) + a_vec
  active <- w > 1e-6 * max(w)
  nu <- stats::median(stat[active])

  list(Xa = Xa, Xall = target_X, w = w, eps = eps, nu = nu, n = n, na = na,
       kernel = kernel, level = level, idx_a = idx_a)
}

#' Predict binary EBW weights out of sample
#'
#' Evaluates the dual-representation weight predictor
#' \eqn{w^*(x) = (2\epsilon)^{-1} (\nu - 2 h^*(x))_+} at new covariate rows,
#' where \eqn{h^*} is the balancing function. Does not refit.
#'
#' @param fit the list returned by `.ebw_binary`.
#' @param Xnew already-scaled new covariate matrix.
#' @param type `"weights"` for predicted weights or `"balancing"` for the
#'   balancing function \eqn{h^*(x)}.
#'
#' @return A numeric vector of predicted weights (or balancing-function values).
#'
#' @keywords internal
#' @noRd
.predict_binary <- function(fit, Xnew, type = "weights") {
  Xnew <- as.matrix(Xnew)
  Kxa <- fit$kernel(Xnew, fit$Xa)       # m x na
  Kxall <- fit$kernel(Xnew, fit$Xall)   # m x n
  h <- drop(Kxa %*% fit$w) - rowSums(Kxall) / fit$n
  if (type == "balancing") return(h)
  pmax(fit$nu - 2 * h, 0) / (2 * fit$eps)
}
