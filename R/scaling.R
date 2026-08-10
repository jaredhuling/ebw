# Covariate scaling helpers. The fit stores the center/scale so that predict
# can re-apply exactly the same transformation to new data.

#' Scale a covariate matrix
#'
#' Centers and scales the columns of a covariate matrix by one of several
#' rules, returning the transformed matrix together with the center and scale
#' vectors so that the same transformation can later be re-applied to new data.
#'
#' @param X numeric matrix of covariates.
#' @param method one of `"std"` (mean/standard deviation), `"mad"`
#'   (median/median absolute deviation), `"range"` (min/range), or `"none"`
#'   (no transformation).
#'
#' @return A list with elements `X` (the scaled matrix), `center`, and `scale`.
#'
#' @keywords internal
#' @noRd
.scale_covariates <- function(X, method = c("std", "mad", "range", "none")) {
  method <- match.arg(method)
  X <- as.matrix(X)
  p <- ncol(X)
  if (method == "none") {
    center <- rep(0, p); scl <- rep(1, p)
  } else if (method == "std") {
    center <- colMeans(X)
    scl <- apply(X, 2, stats::sd)
  } else if (method == "mad") {
    center <- apply(X, 2, stats::median)
    scl <- apply(X, 2, stats::mad)
  } else if (method == "range") {
    center <- apply(X, 2, min)
    scl <- apply(X, 2, function(z) diff(range(z)))
  }
  scl[!is.finite(scl) | scl == 0] <- 1
  Xs <- .apply_scaling(X, center, scl)
  list(X = Xs, center = center, scale = scl)
}

#' Re-apply a stored center/scale to new data
#'
#' @param X numeric matrix of covariates.
#' @param center numeric vector of column centers.
#' @param scale numeric vector of column scales.
#'
#' @return The scaled matrix `(X - center) / scale`, applied column-wise.
#'
#' @keywords internal
#' @noRd
.apply_scaling <- function(X, center, scale) {
  X <- as.matrix(X)
  sweep(sweep(X, 2, center, "-"), 2, scale, "/")
}
