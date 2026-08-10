# Internal engine for continuous-treatment independence weights (IW). Ports
# iw_fit from the validated dual-representation prototype. Expects an
# already-scaled covariate matrix; scaling is handled by the caller.

#' Fit continuous-treatment independence weights
#'
#' Solves the simplex-constrained ridge QP that drives the weighted joint
#' distribution of covariates and treatment toward independence, and returns
#' the fitted weights with the double-centering constants and multiplier needed
#' for out-of-sample prediction.
#'
#' @param X scaled numeric covariate matrix (n x p).
#' @param A numeric treatment vector.
#' @param lambda ridge parameter; defaults to `sqrt(n)` when `NULL`.
#' @param cX,cA marginal-preservation weights; default to the dimension-adjusted
#'   constants when `NULL`.
#' @param kernel_X,kernel_A covariate and treatment kernels.
#' @param dimension_adj logical; if `FALSE`, use equal marginal-preservation
#'   weights `cX = cA = 0.5` instead of the dimension-adjusted default.
#' @param preserve_means logical; add linear equality constraints forcing the
#'   weighted mean of every covariate and of the treatment to equal the
#'   unweighted (sample) mean. Solved with `osqp` (see Details).
#' @param decorrelate_moments logical; add linear equality constraints forcing
#'   the weighted covariance between the treatment and every covariate to zero.
#'   Solved with `osqp`.
#' @param max_iter,tol FISTA controls.
#'
#' @details When either `preserve_means` or `decorrelate_moments` is `TRUE` the
#'   QP is solved with equality constraints \eqn{B w = b} through `osqp`, and the
#'   out-of-sample predictor gains a dual term linear in the constraint features:
#'   \eqn{w^*(x,a) = (n^2/2\lambda)(\nu + \sum_m \zeta_m \phi_m(x,a) - g(x,a))_+},
#'   where \eqn{(\nu, \zeta)} are recovered as the `osqp` constraint duals. The
#'   feature maps are \eqn{\phi(x,a) = x_k} for a preserved covariate mean,
#'   \eqn{\phi(x,a) = a} for the preserved treatment mean, and
#'   \eqn{\phi(x,a) = (a-\bar A)(x_k-\bar X_k)} for a decorrelation constraint.
#'   The right-hand side preserves the WEIGHTED mean under \eqn{\sum_i w_i = n},
#'   i.e. \eqn{\sum_i w_i X_{ik} = \sum_i X_{ik}}.
#'
#' @return A list carrying the fitted weights `w`, the training covariates `X`
#'   and treatment `A`, the centering constants, the multiplier `nu`, and the
#'   scalars needed by the predictor. When constrained, it also carries the
#'   constraint-feature specification `feat`, the constraint duals `zeta`, and
#'   the training centers `Xbar`, `Abar`.
#'
#' @keywords internal
#' @noRd
.ebw_continuous <- function(X, A, lambda = NULL, cX = NULL, cA = NULL,
                            kernel_X = energy_kernel_X, kernel_A = energy_kernel_A,
                            dimension_adj = TRUE,
                            preserve_means = FALSE, decorrelate_moments = FALSE,
                            max_iter = 20000L, tol = 1e-11) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  A <- as.numeric(A)
  if (is.null(cA) || is.null(cX)) {
    if (dimension_adj) {
      cA <- (1 / sqrt(p)) / (1 / sqrt(p) + 1)
      cX <- 1 / (1 / sqrt(p) + 1)
    } else {
      cA <- 0.5; cX <- 0.5
    }
  }
  if (is.null(lambda)) lambda <- sqrt(n)

  K_X <- kernel_X(X, X); K_A <- kernel_A(A, A)
  muX <- rowMeans(K_X); gmX <- mean(muX)
  muA <- rowMeans(K_A); gmA <- mean(muA)

  KXc <- K_X - outer(muX, muX, "+") + gmX     # double-centered K_X
  KAc <- K_A - outer(muA, muA, "+") + gmA
  M <- KXc * KAc + cA * K_A + cX * K_X
  q_vec <- -2 * (cA * rowSums(K_A) + cX * rowSums(K_X)) / n^2

  Xbar <- colMeans(X); Abar <- mean(A)
  constrained <- isTRUE(preserve_means) || isTRUE(decorrelate_moments)

  base_fit <- list(X = X, A = A, lambda = lambda, cX = cX, cA = cA, n = n,
                   kernel_X = kernel_X, kernel_A = kernel_A,
                   muX = muX, gmX = gmX, muA = muA, gmA = gmA,
                   Xbar = Xbar, Abar = Abar,
                   feat = list(), zeta = numeric(0), constrained = FALSE)

  if (!constrained) {
    Qmv <- function(w) (drop(M %*% w) + lambda * w) / n^2
    sol <- .fista_simplex(Qmv, q_vec, n_w = n, s = n, w0 = rep(1, n),
                          project_momentum = TRUE, max_iter = max_iter, tol = tol)
    w <- sol$w
    # Balancing function g_i; stationarity on the active set:
    # g_i + (2 lambda/n^2) w_i = nu.
    g_train <- 2 * drop(M %*% w) / n^2 + q_vec
    active <- w > 1e-8
    nu <- mean((g_train + (2 * lambda / n^2) * w)[active])
    return(c(base_fit, list(w = w, nu = nu)))
  }

  # ---- constrained path (osqp) ---------------------------------------------
  if (!requireNamespace("osqp", quietly = TRUE))
    stop("`preserve_means` / `decorrelate_moments` require the 'osqp' package. ",
         "Install it with install.packages('osqp').", call. = FALSE)

  # Assemble the constraint-feature matrix Phi (n x M) and RHS b, plus the
  # feature specification reused by the out-of-sample predictor.
  Phi <- NULL; b <- NULL; feat <- list()
  if (isTRUE(preserve_means)) {
    Phi <- cbind(Phi, X, A)
    b <- c(b, colSums(X), sum(A))            # preserve the WEIGHTED mean at sum(w)=n
    feat <- c(feat,
              lapply(seq_len(p), function(k) list(type = "meanX", k = k)),
              list(list(type = "meanA")))
  }
  if (isTRUE(decorrelate_moments)) {
    Cm <- (A - Abar) * sweep(X, 2, Xbar, "-")
    Phi <- cbind(Phi, Cm)
    b <- c(b, rep(0, p))
    feat <- c(feat, lapply(seq_len(p), function(k) list(type = "dec", k = k)))
  }

  # osqp form: min 1/2 w' Posqp w + q'w  s.t. l <= Aosqp w <= u
  Posqp <- 2 * (M + lambda * diag(n)) / n^2
  Aosqp <- rbind(diag(n), rep(1, n), t(Phi))
  lvec <- c(rep(0, n), n, b)
  uvec <- c(rep(Inf, n), n, b)

  sol <- osqp::solve_osqp(Posqp, q = q_vec, A = Aosqp, l = lvec, u = uvec,
                          pars = osqp::osqpSettings(max_iter = 2e5L,
                                                    eps_abs = 1e-9, eps_rel = 1e-9,
                                                    verbose = FALSE))
  w <- pmax(sol$x, 0)
  y <- sol$y                                 # constraint duals, same row order
  nu <- -y[n + 1]                            # dual of the sum-to-n row
  zeta <- -y[(n + 2):(n + 1 + ncol(Phi))]    # duals of the B rows

  base_fit$feat <- feat
  base_fit$zeta <- zeta
  base_fit$constrained <- TRUE
  c(base_fit, list(w = w, nu = nu))
}

# Evaluate the constraint features phi_m at new (x, a). Xnew is in the fitted
# (scaled) covariate space; Anew is the raw treatment. Returns an m x M matrix.
#' @keywords internal
#' @noRd
.phi_eval <- function(feat, Xnew, Anew, Xbar, Abar) {
  if (length(feat) == 0L) return(matrix(0, nrow(Xnew), 0))
  cols <- lapply(feat, function(f) {
    if (f$type == "meanX") Xnew[, f$k]
    else if (f$type == "meanA") Anew
    else (Anew - Abar) * (Xnew[, f$k] - Xbar[f$k])   # dec
  })
  matrix(unlist(cols), nrow = nrow(Xnew), ncol = length(feat))
}

#' Predict continuous IW weights out of sample
#'
#' Evaluates the dual-representation weight predictor
#' \eqn{w^*(x,a) = (n^2 / 2\lambda)(\nu - g(x,a))_+} at new (covariate,
#' treatment) pairs. Row centering is recomputed for the new points while the
#' column and grand means are the training constants stored in `fit`. Does not
#' refit.
#'
#' @param fit the list returned by `.ebw_continuous`.
#' @param Xnew already-scaled new covariate matrix.
#' @param Anew new treatment vector.
#' @param type `"weights"` for predicted weights or `"balancing"` for the
#'   balancing function \eqn{g(x,a)}.
#'
#' @return A numeric vector of predicted weights (or balancing-function values).
#'
#' @keywords internal
#' @noRd
.predict_continuous <- function(fit, Xnew, Anew, type = "weights") {
  Xnew <- as.matrix(Xnew)
  Anew <- as.numeric(Anew)
  m <- nrow(Xnew)
  KxX <- fit$kernel_X(Xnew, fit$X)      # m x n
  KaA <- fit$kernel_A(Anew, fit$A)      # m x n
  muX_new <- rowMeans(KxX)              # length m
  muA_new <- rowMeans(KaA)
  KxXc <- KxX - muX_new - outer(rep(1, m), fit$muX) + fit$gmX
  KaAc <- KaA - muA_new - outer(rep(1, m), fit$muA) + fit$gmA
  Mrow <- KxXc * KaAc + fit$cA * KaA + fit$cX * KxX
  Mw <- drop(Mrow %*% fit$w)
  q_new <- -2 * (fit$cA * rowSums(KaA) + fit$cX * rowSums(KxX)) / fit$n^2
  g_new <- 2 * Mw / fit$n^2 + q_new
  if (type == "balancing") return(g_new)
  # Augmented dual term, linear in the constraint features (zero when the fit is
  # unconstrained). Evaluates the same phi_m used at fit time, out of sample.
  lin <- 0
  if (length(fit$feat)) {
    Phi_new <- .phi_eval(fit$feat, Xnew, Anew, fit$Xbar, fit$Abar)
    lin <- drop(Phi_new %*% fit$zeta)
  }
  pmax(fit$nu + lin - g_new, 0) * (fit$n^2 / (2 * fit$lambda))
}
