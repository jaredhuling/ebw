#' ebw: Energy Balancing and Independence Weights
#'
#' Estimation of covariate balancing weights for causal inference with either a
#' binary treatment (energy balancing weights) or a continuous treatment
#' (independence weights), through the single function [energy_balance()]. A
#' fitted weighting carries a dual-representation balancing function that
#' extends to new covariate values, so [predict.ebw()] returns the weight a new
#' unit would receive out of sample.
#'
#' @keywords internal
"_PACKAGE"
