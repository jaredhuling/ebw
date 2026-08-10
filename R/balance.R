# Balance diagnostics for a fitted "ebw" weighting. Works both in sample (on
# the training data stored in the fit) and out of sample (on fresh data under
# the predicted weights), so a user can check whether a fitted weighting still
# balances a new sample.

# Null-coalescing helper.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- distance-matrix primitives -------------------------------------------

#' @keywords internal
#' @noRd
.dist_euclid <- function(X) {
  X <- as.matrix(X)
  ss <- rowSums(X^2)
  d2 <- outer(ss, ss, "+") - 2 * tcrossprod(X)
  sqrt(pmax(d2, 0))
}

#' @keywords internal
#' @noRd
.cross_euclid <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B)
  d2 <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * tcrossprod(A, B)
  sqrt(pmax(d2, 0))
}

#' @keywords internal
#' @noRd
.dist_abs <- function(a) abs(outer(as.numeric(a), as.numeric(a), "-"))

# Weighted (biased) distance covariance / correlation from two distance
# matrices, with probability weights p = w / sum(w). Uniform w recovers the
# usual unweighted statistic.
#' @keywords internal
#' @noRd
.wtd_dcov_stats <- function(Dx, Dy, w) {
  p <- w / sum(w)
  ax <- drop(Dx %*% p); axbar <- sum(p * ax)
  ay <- drop(Dy %*% p); aybar <- sum(p * ay)
  Ax <- Dx - outer(ax, ax, "+") + axbar
  Ay <- Dy - outer(ay, ay, "+") + aybar
  dcov2  <- drop(p %*% (Ax * Ay) %*% p)
  dvarx2 <- drop(p %*% (Ax * Ax) %*% p)
  dvary2 <- drop(p %*% (Ay * Ay) %*% p)
  denom <- sqrt(dvarx2 * dvary2)
  dcor <- if (denom > 0) sqrt(max(dcov2 / denom, 0)) else 0
  list(dcov = sqrt(max(dcov2, 0)), dcor = dcor)
}

# Energy distance between two weighted empirical measures.
# ED = sqrt( 2 E||Xa-Xb|| - E||Xa-Xa'|| - E||Xb-Xb'|| ).
#' @keywords internal
#' @noRd
.energy_distance <- function(Xa, wa, Xb, wb) {
  Dab <- .cross_euclid(Xa, Xb)
  Daa <- .dist_euclid(Xa)
  Dbb <- .dist_euclid(Xb)
  cross <- drop(wa %*% Dab %*% wb)
  aa <- drop(wa %*% Daa %*% wa)
  bb <- drop(wb %*% Dbb %*% wb)
  sqrt(max(2 * cross - aa - bb, 0))
}

# ---- data extraction (in sample or out of sample) -------------------------

# Return the raw covariate matrix and treatment vector on which balance is to
# be assessed. With newdata NULL, the training data stored in the fit. With
# newdata, the fresh covariates (unscaled back to the raw space) and treatment.
#' @keywords internal
#' @noRd
.balance_data <- function(object, newdata) {
  if (is.null(newdata)) {
    return(list(X = object$x, treatment = object$treatment, newdata = FALSE))
  }
  nd <- .build_newdata(object, newdata)          # scaled X (+ A for continuous)
  Xraw <- sweep(sweep(nd$X, 2, object$scaling$scale, "*"),
                2, object$scaling$center, "+")
  colnames(Xraw) <- colnames(object$x)
  trt <- nd$A
  if (is.null(trt)) {
    if (!is.null(object$formula) && is.data.frame(newdata)) {
      rn <- all.vars(object$formula[[2]])[1]
      trt <- newdata[[rn]]
    } else if (is.list(newdata) && !is.data.frame(newdata)) {
      trt <- newdata$trt %||% newdata$a %||% newdata$treatment
    }
  }
  list(X = Xraw, treatment = trt, newdata = TRUE)
}

# ---- balance() ------------------------------------------------------------

#' Balance diagnostics for a fitted weighting
#'
#' Reports how well a fitted [energy_balance()] weighting balances the covariate
#' distribution, either on the training sample or, when `newdata` is supplied,
#' on a fresh sample under the weights the fitted balancing function predicts
#' for it. The latter checks whether a weighting learned on one sample still
#' balances another, which is the out-of-sample analogue of the usual balance
#' table.
#'
#' For a binary treatment the diagnostics are per-covariate standardized mean
#' differences (SMD) of the balanced group against the target measure (the
#' pooled sample for the ATE, the treated sample for the ATT), both before and
#' after weighting, together with the energy distance between the weighted
#' balanced group and the target. For a continuous treatment the diagnostics are
#' per-covariate weighted distance correlations between each covariate and the
#' treatment, the overall weighted distance covariance and correlation, and the
#' weighted versus unweighted marginal means. Every report includes the
#' effective sample size \eqn{\mathrm{ESS} = (\sum_i w_i)^2 / \sum_i w_i^2}.
#'
#' @param object a fitted `"ebw"` object.
#' @param newdata optional new data on which to assess balance. If omitted,
#'   balance is assessed on the training sample. Accepts the same forms as
#'   [predict.ebw()]: a data frame for a formula fit, or a list with elements
#'   `x` and (for a continuous treatment) `a` for the matrix interface. For a
#'   binary out-of-sample assessment the treatment must also be available (the
#'   response column of the data frame, or `trt`/`a` in the list).
#' @param weights optional numeric weight vector to use instead of the fitted or
#'   predicted weights. Must have one entry per row of the assessed sample.
#'
#' @return A list of class `"balance.ebw"` with a `per_covariate` data frame and
#'   the effective sample size `ess`, plus, for a binary treatment, the
#'   `energy_distance` before and after weighting and, for a continuous
#'   treatment, the overall distance `dcov`/`dcor` and a `means` data frame.
#'
#' @seealso [plot.ebw()] for a Love plot of the same quantities.
#'
#' @examples
#' set.seed(20260809)
#' dat <- sim_confounded_data(n.obs = 200, n.vars = 6)
#' fit <- energy_balance(x = dat$x, treatment = dat$trt)
#' balance(fit)
#'
#' @export
balance <- function(object, newdata = NULL, weights = NULL) {
  if (!inherits(object, "ebw")) stop("`object` must be an 'ebw' fit.")
  dat <- .balance_data(object, newdata)
  X <- as.matrix(dat$X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

  w <- weights
  if (is.null(w)) {
    w <- if (dat$newdata) predict(object, newdata = newdata) else object$weights
  }
  w <- as.numeric(w)
  sample_label <- if (dat$newdata) "out-of-sample" else "in-sample"

  kind <- object$engine$kind
  if (kind == "continuous") {
    .balance_continuous(X, dat$treatment, w, sample_label)
  } else if (kind == "binary") {
    .balance_binary(X, dat$treatment, w, object, sample_label)
  } else {
    .balance_grouped(X, dat$treatment, w, object, sample_label)
  }
}

# Pairwise weighted energy distances between the groups. `groups` is an integer
# vector in 1..K, `w` the (per-group-normalized-internally) weights. Returns a
# data frame with one row per pair and the unadjusted / adjusted energy distance.
#' @keywords internal
#' @noRd
.between_group_ed <- function(X, groups, w, level_names = NULL) {
  X <- as.matrix(X)
  K <- max(groups)
  idxs <- lapply(seq_len(K), function(k) which(groups == k))
  rows <- list()
  for (k in seq_len(K)) for (j in seq_len(K)) if (j < k) {
    ij <- idxs[[j]]; ik <- idxs[[k]]
    if (!length(ij) || !length(ik)) next
    Xj <- X[ij, , drop = FALSE]; Xk <- X[ik, , drop = FALSE]
    wj <- w[ij]; wk <- w[ik]
    nj <- length(ij); nk <- length(ik)
    un <- .energy_distance(Xj, rep(1 / nj, nj), Xk, rep(1 / nk, nk))
    ad <- .energy_distance(Xj, wj / sum(wj),    Xk, wk / sum(wk))
    lj <- if (is.null(level_names)) j else level_names[j]
    lk <- if (is.null(level_names)) k else level_names[k]
    rows[[length(rows) + 1L]] <- data.frame(
      group_j = lj, group_k = lk, unadjusted = un, adjusted = ad,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' @keywords internal
#' @noRd
.balance_binary <- function(X, trt, w, object, sample_label) {
  if (is.null(trt))
    stop("Binary balance requires the treatment labels for the assessed sample.")
  eng <- object$engine
  estimand <- object$estimand
  primary_level <- eng$primary_level
  control_level <- eng$control_level

  if (identical(estimand, "ATT")) {
    balanced_level <- control_level
    tgt_idx <- which(trt == primary_level)      # target = treated sample
  } else {
    balanced_level <- primary_level
    tgt_idx <- seq_len(nrow(X))                 # target = pooled sample
  }
  bidx <- which(trt == balanced_level)
  if (!length(bidx))
    stop("No units of the balanced treatment group are present in the sample.")
  Xb <- X[bidx, , drop = FALSE]; wb <- w[bidx]
  Xt <- X[tgt_idx, , drop = FALSE]

  tgt_mean <- colMeans(Xt)
  tgt_sd <- apply(Xt, 2, stats::sd)
  tgt_sd[!is.finite(tgt_sd) | tgt_sd == 0] <- 1
  raw_mean <- colMeans(Xb)
  wtd_mean <- colSums(wb * Xb) / sum(wb)

  per_cov <- data.frame(
    variable   = colnames(X),
    unadjusted = (raw_mean - tgt_mean) / tgt_sd,
    adjusted   = (wtd_mean - tgt_mean) / tgt_sd,
    row.names  = NULL, stringsAsFactors = FALSE)

  nb <- nrow(Xb); nt <- nrow(Xt)
  ed_unadj <- .energy_distance(Xb, rep(1 / nb, nb), Xt, rep(1 / nt, nt))
  ed_adj   <- .energy_distance(Xb, wb / sum(wb),   Xt, rep(1 / nt, nt))

  # between-group energy distance (the two treatment groups against each other)
  levs <- eng$levels
  groups2 <- match(trt, levs)
  btw <- .between_group_ed(X, groups2, w, level_names = levs)

  structure(list(
    treatment_type = "binary",
    estimand = estimand,
    metric = "standardized mean difference (SMD)",
    sample = sample_label,
    balanced_group = balanced_level,
    per_covariate = per_cov,
    energy_distance = c(unadjusted = ed_unadj, adjusted = ed_adj),
    between_group = btw,
    ess = sum(wb)^2 / sum(wb^2),
    n = nb
  ), class = "balance.ebw")
}

# Balance diagnostics for a grouped fit (improved binary or multi-category):
# each group's covariate means against the pooled sample, before and after
# weighting, plus the pairwise between-group energy distances.
#' @keywords internal
#' @noRd
.balance_grouped <- function(X, trt, w, object, sample_label) {
  if (is.null(trt))
    stop("Grouped balance requires the treatment labels for the assessed sample.")
  X <- as.matrix(X)
  levs <- object$treatment_levels
  groups <- match(trt, levs)
  if (anyNA(groups))
    stop("Assessed-sample treatment values must be among the fitted levels.")
  K <- length(levs)
  n <- nrow(X); p <- ncol(X)

  pooled_mean <- colMeans(X)
  pooled_sd <- apply(X, 2, stats::sd)
  pooled_sd[!is.finite(pooled_sd) | pooled_sd == 0] <- 1

  pc_list <- list(); ess <- numeric(K)
  for (k in seq_len(K)) {
    ik <- which(groups == k)
    Xk <- X[ik, , drop = FALSE]; wk <- w[ik]
    raw_mean <- colMeans(Xk)
    wtd_mean <- colSums(wk * Xk) / sum(wk)
    pc_list[[k]] <- data.frame(
      variable   = colnames(X),
      group      = rep(levs[k], p),
      unadjusted = (raw_mean - pooled_mean) / pooled_sd,
      adjusted   = (wtd_mean - pooled_mean) / pooled_sd,
      row.names  = NULL, stringsAsFactors = FALSE)
    ess[k] <- sum(wk)^2 / sum(wk^2)
  }
  per_cov <- do.call(rbind, pc_list)
  names(ess) <- as.character(levs)
  btw <- .between_group_ed(X, groups, w, level_names = levs)

  structure(list(
    treatment_type = object$treatment_type,
    estimand = "ATE",
    improved = isTRUE(object$improved),
    metric = "standardized mean difference (SMD) vs pooled sample",
    sample = sample_label,
    n_groups = K,
    group_levels = levs,
    per_covariate = per_cov,
    between_group = btw,
    ess = ess,
    n = n
  ), class = "balance.ebw")
}

#' @keywords internal
#' @noRd
.balance_continuous <- function(X, A, w, sample_label) {
  if (is.null(A))
    stop("Continuous balance requires the treatment values for the assessed sample.")
  A <- as.numeric(A)
  n <- nrow(X); p <- ncol(X)
  unif <- rep(1, n)                             # uniform (unweighted)

  Dy <- .dist_abs(A)
  dcor_un <- dcor_ad <- numeric(p)
  for (k in seq_len(p)) {
    Dxk <- .dist_abs(X[, k])
    dcor_un[k] <- .wtd_dcov_stats(Dxk, Dy, unif)$dcor
    dcor_ad[k] <- .wtd_dcov_stats(Dxk, Dy, w)$dcor
  }
  per_cov <- data.frame(
    variable   = colnames(X),
    unadjusted = dcor_un,
    adjusted   = dcor_ad,
    row.names  = NULL, stringsAsFactors = FALSE)

  Dx <- .dist_euclid(X)
  ov_un <- .wtd_dcov_stats(Dx, Dy, unif)
  ov_ad <- .wtd_dcov_stats(Dx, Dy, w)

  means <- data.frame(
    variable   = c(colnames(X), "(treatment)"),
    unweighted = c(colMeans(X), mean(A)),
    weighted   = c(colSums(w * X) / sum(w), sum(w * A) / sum(w)),
    row.names  = NULL, stringsAsFactors = FALSE)

  structure(list(
    treatment_type = "continuous",
    metric = "distance correlation with treatment",
    sample = sample_label,
    per_covariate = per_cov,
    dcov = c(unadjusted = ov_un$dcov, adjusted = ov_ad$dcov),
    dcor = c(unadjusted = ov_un$dcor, adjusted = ov_ad$dcor),
    means = means,
    ess = sum(w)^2 / sum(w^2),
    n = n
  ), class = "balance.ebw")
}

#' Print balance diagnostics
#'
#' @param x a `"balance.ebw"` object.
#' @param digits number of significant digits.
#' @param ... unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.balance.ebw <- function(x, digits = 3, ...) {
  grouped <- !is.null(x$n_groups)
  cat("Balance diagnostics (", x$sample, ")\n", sep = "")
  cat("  treatment type:", x$treatment_type, "\n")
  if (grouped) {
    cat("  estimand:      ", x$estimand, "\n")
    cat("  groups:        ", x$n_groups,
        if (isTRUE(x$improved)) "  (improved: cross-group term on)" else "", "\n")
    cat("  per-group ESS: ",
        paste(sprintf("%s=%.1f", names(x$ess), x$ess), collapse = ", "), "\n")
  } else if (identical(x$treatment_type, "binary")) {
    cat("  estimand:      ", x$estimand, "\n")
    cat("  balanced group:", x$balanced_group, "\n")
    cat(sprintf("  ESS:            %.1f  (n = %d)\n", x$ess, x$n))
  } else {
    cat(sprintf("  ESS:            %.1f  (n = %d)\n", x$ess, x$n))
  }
  cat("\n  per covariate (", x$metric, "):\n", sep = "")
  pc <- x$per_covariate
  pc$unadjusted <- signif(pc$unadjusted, digits)
  pc$adjusted <- signif(pc$adjusted, digits)
  print(pc, row.names = FALSE)
  if (!is.null(x$between_group)) {
    cat("\n  between-group energy distance:\n")
    bg <- x$between_group
    bg$unadjusted <- signif(bg$unadjusted, digits)
    bg$adjusted <- signif(bg$adjusted, digits)
    print(bg, row.names = FALSE)
  }
  if (identical(x$treatment_type, "binary") && !is.null(x$energy_distance)) {
    cat(sprintf("\n  energy distance to target: unadjusted %.4g, adjusted %.4g\n",
                x$energy_distance[["unadjusted"]], x$energy_distance[["adjusted"]]))
  } else if (identical(x$treatment_type, "continuous")) {
    cat(sprintf("\n  overall distance covariance: unadjusted %.4g, adjusted %.4g\n",
                x$dcov[["unadjusted"]], x$dcov[["adjusted"]]))
    cat(sprintf("  overall distance correlation: unadjusted %.4g, adjusted %.4g\n",
                x$dcor[["unadjusted"]], x$dcor[["adjusted"]]))
  }
  invisible(x)
}
