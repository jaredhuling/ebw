# Out-of-sample weight prediction: the headline feature. The fitted weighting
# carries a dual-representation balancing function that extends to any new
# covariate value (binary) or covariate/treatment pair (continuous), so a
# fitted object yields a function predicting the weight a new unit receives.

# Rebuild and scale a new design matrix (and, for continuous fits, extract the
# new treatment values) using only the ingredients stored at fit time.
#' @keywords internal
#' @noRd
.build_newdata <- function(object, newdata) {
  A <- NULL
  if (!is.null(object$terms)) {
    # formula fit: rebuild the design matrix via the stored terms/xlevels. When
    # the treatment column is present in `newdata` we also recover it, which is
    # what the grouped / binary predictor needs to map each unit to its group.
    if (!is.null(object$formula)) {
      rname <- all.vars(object$formula[[2]])[1]
      if (is.data.frame(newdata) && !is.null(newdata[[rname]]))
        A <- newdata[[rname]]
    }
    mf <- stats::model.frame(object$terms, data = newdata,
                             xlev = object$xlevels)
    X <- stats::model.matrix(object$terms, mf)
    ic <- which(colnames(X) == "(Intercept)")
    if (length(ic)) X <- X[, -ic, drop = FALSE]
  } else if (is.list(newdata) && !is.data.frame(newdata) &&
             !is.null(newdata$x)) {
    X <- as.matrix(newdata$x)
    A <- newdata$treatment %||% newdata$trt %||% newdata$a
  } else {
    X <- as.matrix(newdata)
  }
  Xs <- .apply_scaling(X, object$scaling$center, object$scaling$scale)
  list(X = Xs, A = A)
}

#' Predict out-of-sample balancing weights
#'
#' Predicts the weight a new unit would receive under a fitted [energy_balance()]
#' weighting, using only the stored dual-representation ingredients (the
#' balancing function \eqn{h^*(x)} for a binary treatment, or \eqn{g(x, a)} for
#' a continuous treatment). No refitting is performed.
#'
#' With no `newdata`, the fitted in-sample weights are returned (identical to
#' `weights(object)`). With `newdata`, weights are predicted for units that need
#' not appear in the training sample. For the formula interface, `newdata` is a
#' data frame from which the design matrix is rebuilt through the stored terms
#' and factor levels; for the matrix interface, `newdata` may be a matrix of
#' covariates or a list with elements `x` and (for a continuous treatment) `a`.
#' A continuous-treatment prediction requires the new treatment values, since
#' the balancing function depends on both the covariates and the treatment.
#'
#' For a binary or multi-category treatment the balancing function is specific
#' to a group, so out-of-sample prediction needs each new unit's treatment
#' level. Supply it through the response column of `newdata` (formula fits), the
#' `treatment` element of a list, or the separate `treatment` argument. For a
#' standard binary fit, if no treatment is supplied the primary (treated) group's
#' predictor is used, matching earlier releases.
#'
#' @param object a fitted `"ebw"` object.
#' @param newdata new data for prediction (see Details). If omitted, the fitted
#'   weights are returned.
#' @param type `"weights"` (default) to return predicted weights, or
#'   `"balancing"` to return the balancing function value at each new point.
#' @param treatment optional treatment vector for the new units (binary or
#'   multi-category fits), used to map each unit to its group when it is not
#'   already carried in `newdata`.
#' @param ... unused.
#'
#' @return A numeric vector of predicted weights (or balancing-function values).
#'
#' @seealso [energy_balance()]
#'
#' @examples
#' set.seed(20260809)
#' n <- 200; p <- 4
#' X <- matrix(rnorm(n * p), n, p)
#' a <- rbinom(n, 1, plogis(X[, 1] - X[, 2]))
#' fit <- energy_balance(x = X, treatment = a)
#' Xnew <- matrix(rnorm(20 * p), 20, p)
#' predict(fit, newdata = Xnew)
#'
#' @export
predict.ebw <- function(object, newdata, type = c("weights", "balancing"),
                        treatment = NULL, ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    return(object$weights)
  }
  nd <- .build_newdata(object, newdata)
  kind <- object$engine$kind
  A <- if (!is.null(treatment)) treatment else nd$A

  if (kind == "continuous") {
    if (is.null(A))
      stop("Continuous-treatment prediction requires new treatment values ",
           "(supply them in `newdata`, e.g. list(x = ..., a = ...)).")
    return(.predict_continuous(object$engine$fit, nd$X, A, type = type))
  }

  if (kind == "binary") {
    # Group-aware prediction when the treatment is available under the ATE;
    # otherwise fall back to the primary group's predictor (backward compatible).
    if (is.null(A) || !identical(object$estimand, "ATE"))
      return(.predict_binary(object$engine$primary, nd$X, type = type))
    return(.predict_binary_grouped(object, nd$X, A, type = type))
  }

  # grouped (improved binary or multi-category)
  if (is.null(A))
    stop("Prediction for this fit requires the new units' treatment level ",
         "(supply it in `newdata` or via the `treatment` argument).")
  gnew <- match(A, object$treatment_levels)
  if (anyNA(gnew))
    stop("New treatment values must be among the fitted levels: ",
         paste(object$treatment_levels, collapse = ", "), ".")
  .predict_grouped(object$engine, nd$X, gnew, type = type)
}

# Predict standard-binary EBW weights per unit, routing each unit through its
# own group's engine (ATE), so a mixed-treatment newdata is handled correctly.
#' @keywords internal
#' @noRd
.predict_binary_grouped <- function(object, Xnew, A, type = "weights") {
  eng <- object$engine
  out <- numeric(nrow(Xnew))
  for (lv in unique(A)) {
    e <- eng$engines[[as.character(lv)]]
    if (is.null(e))
      stop("Treatment level ", lv, " is not among the fitted levels.")
    ii <- which(A == lv)
    out[ii] <- .predict_binary(e, Xnew[ii, , drop = FALSE], type = type)
  }
  out
}
