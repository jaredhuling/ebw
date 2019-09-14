

#'
#' @export
sim_confounded_data <- function(n.obs, n.vars = 10,
                                y.model          = c("A", "B", "C"),
                                propensity.model = c("I", "II", "III"))
{
    y.model          <- match.arg(y.model)
    propensity.model <- match.arg(propensity.model)

    if (n.vars < 10)
    {
        stop("'n.vars' must be 10 or greater")
    }
}
