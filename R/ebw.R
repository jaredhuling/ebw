# Main user-facing function. It handles both binary
# (energy balancing) and continuous (independence weights) treatments, through
# either a formula or a matrix interface, and stores everything needed for
# out-of-sample weight prediction.

#' Energy balancing and independence weights
#'
#' Fits covariate balancing weights for either a binary treatment (energy
#' balancing weights) or a continuous treatment (independence weights), and
#' returns an object from which the weight a new unit would receive can be
#' predicted out of sample via [predict.ebw()].
#'
#' The function accepts two calling styles. In the formula style, supply a
#' `formula` whose left-hand side is the treatment and whose right-hand side
#' gives the covariates, together with a `data` frame; the design matrix is
#' built with [stats::model.matrix()], so factors and interactions are
#' supported. In the matrix style, supply the covariate matrix `x` and the
#' treatment vector `treatment` directly. For backward compatibility the older
#' argument names `trt`, `x`, `method`, and `standardize` are also accepted.
#'
#' The treatment type is detected automatically: a treatment that is a factor,
#' is logical, or takes exactly two distinct values is treated as binary;
#' a numeric treatment with more than two distinct values is treated as
#' continuous. Use `treatment_type` to override the detection.
#'
#' @param formula a model formula `treatment ~ covariates` (formula interface).
#' @param data a data frame containing the variables in `formula`.
#' @param x covariate matrix (matrix interface).
#' @param treatment treatment vector (matrix interface).
#' @param estimand target estimand for a binary treatment, either `"ATE"`
#'   (each group balanced to the pooled sample) or `"ATT"` (the untreated group
#'   balanced to the treated sample). Ignored for continuous treatments.
#' @param treatment_type one of `"auto"`, `"binary"`, or `"continuous"`.
#' @param scaling covariate scaling rule (`.scale_covariates`): one of
#'   `"std"`, `"mad"`, `"range"`, or `"none"`.
#' @param kernel covariate kernel, either `"energy"` (default) or `"gaussian"`.
#' @param eps ridge parameter for the binary energy-balancing QP.
#' @param lambda ridge parameter for the continuous independence-weights QP;
#'   defaults to `sqrt(n)` when `NULL`.
#' @param dimension_adj logical; use dimension-adjusted marginal-preservation
#'   weights for the continuous engine (default `TRUE`).
#' @param preserve_means logical (continuous treatment only); add equality
#'   constraints forcing the weighted mean of every covariate and of the
#'   treatment to equal the unweighted sample mean. The constrained QP is solved
#'   with `osqp` and the out-of-sample predictor gains a dual term linear in the
#'   constraint features.
#' @param decorrelate_moments logical (continuous treatment only); add equality
#'   constraints forcing the weighted treatment-covariate covariances to zero.
#'   Also solved with `osqp`.
#' @param level optional value of `treatment` identifying the treated group for
#'   a binary treatment; defaults to the larger / second level.
#' @param max_iter,tol FISTA controls.
#' @param trt,method,standardize deprecated argument names retained for
#'   backward compatibility. `trt` maps to `treatment`, `method` (`"ATE"`,
#'   `"ATE.3"`, `"ATT"`) maps to `estimand`, and `standardize` (`TRUE`/`FALSE`)
#'   maps to `scaling` (`"std"`/`"none"`).
#' @param ... additional arguments (currently unused).
#'
#' @return An object of class `"ebw"`: a list with the fitted `weights`, the
#'   detected `treatment_type`, the `estimand`, the fitted `engine` (carrying
#'   the dual-representation ingredients), the `scaling` information, the
#'   `call`, the sample size `n`, and, for the formula interface, the `terms`
#'   and `xlevels` needed to rebuild the design matrix at prediction time.
#'
#' @seealso [predict.ebw()] for out-of-sample weight prediction.
#'
#' @examples
#' set.seed(20260809)
#' n <- 200; p <- 4
#' X <- matrix(rnorm(n * p), n, p)
#' a <- rbinom(n, 1, plogis(X[, 1] - X[, 2]))
#' fit <- energy_balance(x = X, treatment = a)
#' w <- weights(fit)
#'
#' @export
energy_balance <- function(formula, data,
                           x = NULL, treatment = NULL,
                           estimand = c("ATE", "ATT"),
                           treatment_type = c("auto", "binary", "continuous"),
                           scaling = c("std", "mad", "range", "none"),
                           kernel = c("energy", "gaussian"),
                           eps = 1e-2, lambda = NULL,
                           dimension_adj = TRUE,
                           preserve_means = FALSE, decorrelate_moments = FALSE,
                           level = NULL,
                           max_iter = 20000L, tol = 1e-11,
                           trt = NULL, method = NULL, standardize = NULL,
                           ...) {
  cl <- match.call()
  treatment_type <- match.arg(treatment_type)
  kernel <- match.arg(kernel)

  # ---- backward-compatible argument mapping --------------------------------
  if (!is.null(trt) && is.null(treatment)) treatment <- trt
  if (!is.null(method)) {
    method <- match.arg(method, c("ATE.3", "ATE", "ATT"))
    estimand <- if (method == "ATT") "ATT" else "ATE"
  } else {
    estimand <- match.arg(estimand)
  }
  if (!is.null(standardize)) {
    scaling <- if (isTRUE(standardize)) "std" else "none"
  } else {
    scaling <- match.arg(scaling)
  }

  # ---- assemble design matrix and treatment --------------------------------
  terms_obj <- NULL; xlevels <- NULL; fml <- NULL
  if (!missing(formula) && inherits(formula, "formula")) {
    if (missing(data)) data <- environment(formula)
    mf <- stats::model.frame(formula, data = data)
    treatment <- stats::model.response(mf)
    terms_obj <- stats::terms(mf)
    # drop the response, keep an intercept-free covariate design matrix
    rhs_terms <- stats::delete.response(terms_obj)
    X <- stats::model.matrix(rhs_terms, data = mf)
    intercept_col <- which(colnames(X) == "(Intercept)")
    if (length(intercept_col)) X <- X[, -intercept_col, drop = FALSE]
    xlevels <- stats::.getXlevels(rhs_terms, mf)
    fml <- formula
    terms_obj <- rhs_terms
  } else {
    if (is.null(x) || is.null(treatment))
      stop("Supply either a `formula` (with `data`) or both `x` and `treatment`.")
    X <- as.matrix(x)
  }
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  n <- nrow(X)

  # ---- detect treatment type -----------------------------------------------
  if (treatment_type == "auto") {
    if (is.factor(treatment) || is.logical(treatment) ||
        length(unique(treatment)) == 2L) {
      treatment_type <- "binary"
    } else {
      treatment_type <- "continuous"
    }
  }

  # ---- scale covariates (store for predict) --------------------------------
  sc <- .scale_covariates(X, method = scaling)
  Xs <- sc$X
  scaling_info <- list(method = scaling, center = sc$center, scale = sc$scale)

  # ---- covariate kernel (fixed bandwidth for gaussian) ---------------------
  if (kernel == "gaussian") {
    n1 <- rowSums(Xs^2)
    dself <- outer(n1, n1, "+") - 2 * Xs %*% t(Xs)
    dself <- sqrt(pmax(dself, 0))
    sig <- stats::median(dself[lower.tri(dself)])
    if (!is.finite(sig) || sig <= 0) sig <- 1
    kern_X <- function(Z1, Z2) gaussian_kernel_X(Z1, Z2, sigma = sig)
  } else {
    kern_X <- energy_kernel_X
  }

  # ---- fit ------------------------------------------------------------------
  if (treatment_type == "binary") {
    engine <- .fit_binary(Xs, treatment, estimand = estimand, level = level,
                          eps = eps, kernel = kern_X,
                          max_iter = max_iter, tol = tol)
    weights <- engine$weights
  } else {
    fit_c <- .ebw_continuous(Xs, treatment, lambda = lambda,
                             kernel_X = kern_X, kernel_A = energy_kernel_A,
                             dimension_adj = dimension_adj,
                             preserve_means = preserve_means,
                             decorrelate_moments = decorrelate_moments,
                             max_iter = max_iter, tol = tol)
    engine <- list(kind = "continuous", fit = fit_c)
    weights <- fit_c$w
  }

  structure(list(
    weights = weights,
    treatment_type = treatment_type,
    estimand = if (treatment_type == "binary") estimand else NA_character_,
    engine = engine,
    scaling = scaling_info,
    kernel = kernel,
    call = cl,
    n = n,
    x = X,                       # raw (unscaled) design matrix, for diagnostics
    treatment = treatment,       # raw treatment vector, for diagnostics
    terms = terms_obj,
    xlevels = xlevels,
    formula = fml
  ), class = "ebw")
}

# Assemble one or more binary engines and the length-n weight vector.
#' @keywords internal
#' @noRd
.fit_binary <- function(Xs, treatment, estimand, level, eps, kernel,
                        max_iter, tol) {
  levs <- sort(unique(treatment))
  if (length(levs) != 2L)
    stop("A binary treatment must take exactly two distinct values.")
  primary_level <- if (is.null(level)) levs[2] else level
  control_level <- levs[levs != primary_level][1]

  n <- nrow(Xs)
  w_full <- numeric(n)
  engines <- list()

  if (estimand == "ATT") {
    idx_trt <- which(treatment == primary_level)
    Xtrt <- Xs[idx_trt, , drop = FALSE]
    eng_ctrl <- .ebw_binary(Xs, treatment, level = control_level,
                            target_X = Xtrt, eps = eps, kernel = kernel,
                            max_iter = max_iter, tol = tol)
    w_full[idx_trt] <- 1 / length(idx_trt)
    w_full[eng_ctrl$idx_a] <- eng_ctrl$w
    engines[[as.character(control_level)]] <- eng_ctrl
    primary <- eng_ctrl
  } else {
    for (lv in levs) {
      eng <- .ebw_binary(Xs, treatment, level = lv, target_X = Xs,
                         eps = eps, kernel = kernel,
                         max_iter = max_iter, tol = tol)
      engines[[as.character(lv)]] <- eng
      w_full[eng$idx_a] <- eng$w
    }
    primary <- engines[[as.character(primary_level)]]
  }

  list(kind = "binary", weights = w_full, engines = engines,
       primary = primary, primary_level = primary_level,
       control_level = control_level, levels = levs, estimand = estimand)
}

#' Alias for [energy_balance()]
#'
#' `ebw()` is a thin wrapper around [energy_balance()] with identical behavior.
#'
#' @param ... arguments passed to [energy_balance()].
#'
#' @return An object of class `"ebw"`.
#'
#' @export
ebw <- function(...) {
  cl <- match.call()
  out <- energy_balance(...)
  out$call <- cl
  out
}

#' Independence weights for a continuous treatment
#'
#' Convenience wrapper that fits independence weights for a continuous
#' treatment, using the argument names of the mature CRAN interface. Returns an
#' `"ebw"` object so that [predict.ebw()] and the other methods apply.
#'
#' @param A numeric treatment vector.
#' @param X covariate matrix.
#' @param lambda ridge parameter; defaults to `sqrt(n)` internally when 0.
#' @param decorrelate_moments logical; add equality constraints forcing the
#'   weighted treatment-covariate covariances to zero (solved with `osqp`).
#' @param preserve_means logical; add equality constraints forcing the weighted
#'   mean of every covariate and of the treatment to equal the unweighted sample
#'   mean (solved with `osqp`). Unlike the CRAN `independenceWeights` default,
#'   the constraint preserves the WEIGHTED mean under `sum(w) = n`, that is
#'   `sum(w_i X_ik) = sum(X_ik)`.
#' @param dimension_adj logical; use dimension-adjusted marginal-preservation
#'   weights (default `TRUE`).
#' @param scaling covariate scaling rule (default `"none"`, matching the
#'   validated continuous prototype).
#' @param ... further arguments passed to [energy_balance()].
#'
#' @return An object of class `"ebw"`.
#'
#' @export
independence_weights <- function(A, X, lambda = 0,
                                 decorrelate_moments = FALSE,
                                 preserve_means = FALSE,
                                 dimension_adj = TRUE,
                                 scaling = "none", ...) {
  cl <- match.call()
  lam <- if (is.null(lambda) || lambda == 0) NULL else lambda
  out <- energy_balance(x = X, treatment = A, treatment_type = "continuous",
                        lambda = lam, dimension_adj = dimension_adj,
                        preserve_means = preserve_means,
                        decorrelate_moments = decorrelate_moments,
                        scaling = scaling, ...)
  out$call <- cl
  out
}
