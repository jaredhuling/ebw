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
    # formula fit: rebuild the design matrix via the stored terms/xlevels
    mf <- stats::model.frame(object$terms, data = newdata,
                             xlev = object$xlevels)
    X <- stats::model.matrix(object$terms, mf)
    ic <- which(colnames(X) == "(Intercept)")
    if (length(ic)) X <- X[, -ic, drop = FALSE]
    if (object$treatment_type == "continuous" &&
        !is.null(object$formula)) {
      rname <- all.vars(object$formula[[2]])[1]
      if (!is.null(newdata[[rname]])) A <- newdata[[rname]]
    }
  } else if (is.list(newdata) && !is.data.frame(newdata) &&
             !is.null(newdata$x)) {
    X <- as.matrix(newdata$x)
    if (!is.null(newdata$a)) A <- newdata$a
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
#' @param object a fitted `"ebw"` object.
#' @param newdata new data for prediction (see Details). If omitted, the fitted
#'   weights are returned.
#' @param type `"weights"` (default) to return predicted weights, or
#'   `"balancing"` to return the balancing function value at each new point.
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
                        ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    return(object$weights)
  }
  nd <- .build_newdata(object, newdata)
  if (object$treatment_type == "binary") {
    .predict_binary(object$engine$primary, nd$X, type = type)
  } else {
    if (is.null(nd$A))
      stop("Continuous-treatment prediction requires new treatment values ",
           "(supply them in `newdata`, e.g. list(x = ..., a = ...)).")
    .predict_continuous(object$engine$fit, nd$X, nd$A, type = type)
  }
}
