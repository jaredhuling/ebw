# Kernel functions for energy balancing and independence weights.
# Internal helpers; not exported. Each returns a cross-kernel matrix between
# the rows of its two arguments.

#' Energy kernel for covariates
#'
#' Computes the energy (distance-induced) kernel between rows of two matrices,
#' \eqn{\kappa(x, x') = \|x\| + \|x'\| - \|x - x'\|}.
#'
#' @param Z1 numeric matrix (m1 x p) of covariate rows.
#' @param Z2 numeric matrix (m2 x p) of covariate rows.
#'
#' @return An m1 x m2 matrix of cross-kernel evaluations.
#'
#' @keywords internal
#' @noRd
energy_kernel_X <- function(Z1, Z2) {
  Z1 <- as.matrix(Z1); Z2 <- as.matrix(Z2)
  n1 <- sqrt(rowSums(Z1^2)); n2 <- sqrt(rowSums(Z2^2))
  cross <- Z1 %*% t(Z2)
  d2 <- outer(n1^2, n2^2, "+") - 2 * cross
  D <- sqrt(pmax(d2, 0))
  outer(n1, n2, "+") - D
}

#' Energy kernel for a scalar treatment
#'
#' Computes the energy kernel between two numeric vectors of scalar treatment
#' values, \eqn{\kappa(a, a') = |a| + |a'| - |a - a'|}.
#'
#' @param z1 numeric vector of length m1.
#' @param z2 numeric vector of length m2.
#'
#' @return An m1 x m2 matrix of cross-kernel evaluations.
#'
#' @keywords internal
#' @noRd
energy_kernel_A <- function(z1, z2) {
  z1 <- as.numeric(z1); z2 <- as.numeric(z2)
  outer(abs(z1), abs(z2), "+") - abs(outer(z1, z2, "-"))
}

#' Gaussian kernel for covariates
#'
#' Computes a Gaussian (radial basis) kernel between rows of two matrices,
#' \eqn{\kappa(x, x') = \exp(-\|x - x'\|^2 / (2 \sigma^2))}. When `sigma` is
#' `NULL` it defaults to the median pairwise distance of the first argument
#' (the classic median heuristic), which keeps the kernel scale-adaptive.
#'
#' @param Z1 numeric matrix (m1 x p) of covariate rows.
#' @param Z2 numeric matrix (m2 x p) of covariate rows.
#' @param sigma positive bandwidth. If `NULL`, uses the median-distance
#'   heuristic computed from `Z1`.
#'
#' @return An m1 x m2 matrix of cross-kernel evaluations.
#'
#' @keywords internal
#' @noRd
gaussian_kernel_X <- function(Z1, Z2, sigma = NULL) {
  Z1 <- as.matrix(Z1); Z2 <- as.matrix(Z2)
  n1 <- rowSums(Z1^2); n2 <- rowSums(Z2^2)
  d2 <- outer(n1, n2, "+") - 2 * Z1 %*% t(Z2)
  d2 <- pmax(d2, 0)
  if (is.null(sigma)) {
    dself <- outer(n1, n1, "+") - 2 * Z1 %*% t(Z1)
    dself <- sqrt(pmax(dself, 0))
    sigma <- stats::median(dself[lower.tri(dself)])
    if (!is.finite(sigma) || sigma <= 0) sigma <- 1
  }
  exp(-d2 / (2 * sigma^2))
}
