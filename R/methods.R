# S3 methods for the "ebw" class.

#' Extract fitted balancing weights
#'
#' Returns a plain numeric vector of the fitted weights, one per training unit,
#' suitable for passing directly to downstream balance tools such as
#' `cobalt::bal.tab()` or `WeightIt`.
#'
#' @param object a fitted `"ebw"` object.
#' @param ... unused.
#'
#' @return The numeric vector of fitted weights.
#'
#' @examples
#' set.seed(20260809)
#' dat <- sim_confounded_data(n.obs = 200, n.vars = 6)
#' fit <- energy_balance(x = dat$x, treatment = dat$trt)
#' w <- weights(fit)
#' \dontrun{
#' # the weights slot in directly to cobalt's balance table
#' if (requireNamespace("cobalt", quietly = TRUE))
#'   cobalt::bal.tab(dat$x, treat = dat$trt, weights = weights(fit))
#' }
#'
#' @export
weights.ebw <- function(object, ...) {
  object$weights
}

# Effective sample sizes per group (integer group index 1..K) for a grouped or
# binary fit, named by the treatment levels.
#' @keywords internal
#' @noRd
.group_ess <- function(x) {
  if (is.null(x$treatment_levels)) return(NULL)
  levs <- x$treatment_levels
  groups <- match(x$treatment, levs)
  ess <- vapply(seq_along(levs), function(k) {
    wk <- x$weights[groups == k]
    sum(wk)^2 / sum(wk^2)
  }, numeric(1))
  stats::setNames(ess, as.character(levs))
}

#' Print a fitted energy-balancing object
#'
#' @param x a fitted `"ebw"` object.
#' @param ... unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.ebw <- function(x, ...) {
  cat("Energy balancing / independence weights (class 'ebw')\n")
  cat("  treatment type:", x$treatment_type, "\n")
  if (!is.na(x$n_groups))
    cat("  groups (K):    ", x$n_groups, "\n")
  if (!identical(x$treatment_type, "continuous")) {
    cat("  estimand:      ", x$estimand, "\n")
    if (isTRUE(x$improved))
      cat("  improved:       yes (between-group energy distance included)\n")
  }
  cat("  kernel:        ", x$kernel, "\n")
  cat("  scaling:       ", x$scaling$method, "\n")
  cat("  n:             ", x$n, "\n")
  w <- x$weights
  ess <- sum(w)^2 / sum(w^2)
  cat(sprintf("  ESS:            %.1f  (%.1f%% of n)\n", ess, 100 * ess / x$n))
  gess <- .group_ess(x)
  if (!is.null(gess))
    cat("  per-group ESS: ",
        paste(sprintf("%s=%.1f", names(gess), gess), collapse = ", "), "\n")
  cat(sprintf("  weight range:   [%.4g, %.4g]\n", min(w), max(w)))
  invisible(x)
}

#' Summarize a fitted energy-balancing object
#'
#' @param object a fitted `"ebw"` object.
#' @param ... unused.
#'
#' @return A list of class `"summary.ebw"` with the sample size, treatment
#'   type, estimand, effective sample size, and weight range.
#'
#' @export
summary.ebw <- function(object, ...) {
  w <- object$weights
  ess <- sum(w)^2 / sum(w^2)
  out <- list(
    n = object$n,
    treatment_type = object$treatment_type,
    estimand = object$estimand,
    n_groups = object$n_groups,
    improved = isTRUE(object$improved),
    group_ess = .group_ess(object),
    kernel = object$kernel,
    ess = ess,
    weight_range = range(w)
  )
  class(out) <- "summary.ebw"
  out
}

#' @rdname summary.ebw
#' @param x a `"summary.ebw"` object.
#' @export
print.summary.ebw <- function(x, ...) {
  cat("Summary of energy balancing / independence weights\n")
  cat("  n:             ", x$n, "\n")
  cat("  treatment type:", x$treatment_type, "\n")
  if (!is.null(x$n_groups) && !is.na(x$n_groups))
    cat("  groups (K):    ", x$n_groups, "\n")
  if (!identical(x$treatment_type, "continuous")) {
    cat("  estimand:      ", x$estimand, "\n")
    if (isTRUE(x$improved))
      cat("  improved:       yes (between-group energy distance included)\n")
  }
  cat("  kernel:        ", x$kernel, "\n")
  cat(sprintf("  ESS:            %.1f  (%.1f%% of n)\n",
              x$ess, 100 * x$ess / x$n))
  if (!is.null(x$group_ess))
    cat("  per-group ESS: ",
        paste(sprintf("%s=%.1f", names(x$group_ess), x$group_ess),
              collapse = ", "), "\n")
  cat(sprintf("  weight range:   [%.4g, %.4g]\n",
              x$weight_range[1], x$weight_range[2]))
  invisible(x)
}
