

#' Printing results for estimated energy balancing weights
#'
#' @description Prints results for energy balancing weights
#'
#' @param x a fitted object from \code{\link[ebw]{energy_balance}}
#' @param digits minimal number of significant digits to print.
#' @param ... further arguments passed to or from \code{\link[base]{print.default}}.
#' @seealso \code{\link[ebw]{energy_balance}} for function which produces energy balancing weights
#' @rdname print
#' @importFrom stats quantile
#' @export
print.energy_balancing_weights <- function(x, digits = max(getOption('digits')-3, 3), ...)
{
    if (x$method == "ATE")
    {
        method <- "Two way balancing"
    } else if (x$method == "ATE.3")
    {
        method <- "Three way balancing"
    } else
    {
        method <- "Two way balancing"
    }


    cat("Estimand: ", x$estimand, "\n")
    if (x$method != "ATT") cat("Method:   ", method, "\n\n")

    trt <- x$treat

    trt.levels <- levels(trt)
    K          <- length(trt.levels)

    cat("Unweighted energy distance: ", round(x$energy_dist_unweighted, digits),
        "\nOptimized energy distance:   ", round(x$energy_dist_optimized, digits), "\n\n")

    cat("Weight ranges:\n")
    for (t in 1:K)
    {
        wts_cur <- x$weights[trt == trt.levels[t]]
        cat("                 Treatment:", trt.levels[t], "\n")
        print(summary(wts_cur, digits = digits), digits = digits)
        if (t < K) cat("\n")
    }
}
