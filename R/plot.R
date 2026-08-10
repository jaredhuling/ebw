# Diagnostic plots for a fitted "ebw" weighting: a Love-style balance plot and
# a distribution-of-weights plot, both built with ggplot2.

#' Diagnostic plots for a fitted weighting
#'
#' Produces a Love plot of covariate balance before and after weighting, or a
#' plot of the distribution of the fitted weights.
#'
#' For `type = "balance"` (the default) each covariate contributes a horizontal
#' pair of points, its imbalance before and after weighting, joined by a
#' segment. For a binary treatment the imbalance is the absolute standardized
#' mean difference against the target; for a continuous treatment it is the
#' absolute distance correlation between the covariate and the treatment. A
#' dashed reference line is drawn at 0.1. For `type = "weights"` the fitted (or
#' predicted) weights are shown as a histogram.
#'
#' @param x a fitted `"ebw"` object.
#' @param type `"balance"` for a Love plot or `"weights"` for the weight
#'   distribution.
#' @param newdata optional new data; when supplied, the plot reflects the
#'   out-of-sample balance (or predicted weights) on that sample, as in
#'   [balance()].
#' @param weights optional weight vector to use instead of the fitted or
#'   predicted weights.
#' @param bins number of histogram bins for `type = "weights"`.
#' @param ... unused.
#'
#' @return A `ggplot` object.
#'
#' @seealso [balance()]
#'
#' @examples
#' set.seed(20260809)
#' dat <- sim_confounded_data(n.obs = 200, n.vars = 6)
#' fit <- energy_balance(x = dat$x, treatment = dat$trt)
#' plot(fit)                       # Love plot
#' plot(fit, type = "weights")     # weight distribution
#'
#' @export
plot.ebw <- function(x, type = c("balance", "weights"),
                     newdata = NULL, weights = NULL, bins = 30, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("`plot.ebw` requires the 'ggplot2' package.", call. = FALSE)
  type <- match.arg(type)

  if (type == "weights") {
    w <- weights
    if (is.null(w)) {
      w <- if (is.null(newdata)) x$weights else predict(x, newdata = newdata)
    }
    df <- data.frame(w = as.numeric(w))
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = w)) +
        ggplot2::geom_histogram(bins = bins, fill = "#3366aa",
                                colour = "white", alpha = 0.9) +
        ggplot2::labs(x = "weight", y = "count",
                      title = "Distribution of fitted weights") +
        ggplot2::theme_bw()
    )
  }

  # ---- balance (Love) plot -------------------------------------------------
  b <- balance(x, newdata = newdata, weights = weights)
  pc <- b$per_covariate
  grouped <- !is.null(pc$group)
  metric_lab <- if (identical(b$treatment_type, "continuous"))
    "absolute distance correlation"
  else if (grouped)
    "absolute standardized mean difference vs pooled sample"
  else "absolute standardized mean difference"

  # For a grouped fit each covariate appears once per group; order covariates by
  # their largest unadjusted imbalance across groups and facet by group.
  if (grouped) {
    agg <- tapply(abs(pc$unadjusted), pc$variable, max)
    ord <- names(sort(agg))
    grp_key <- as.character(pc$group)
  } else {
    ord <- pc$variable[order(abs(pc$unadjusted))]
    grp_key <- rep("", nrow(pc))
  }

  long <- rbind(
    data.frame(variable = pc$variable, group = grp_key, value = abs(pc$unadjusted),
               adjustment = "Unadjusted", stringsAsFactors = FALSE),
    data.frame(variable = pc$variable, group = grp_key, value = abs(pc$adjusted),
               adjustment = "Adjusted", stringsAsFactors = FALSE))
  long$adjustment <- factor(long$adjustment, levels = c("Unadjusted", "Adjusted"))
  long$variable <- factor(long$variable, levels = ord)

  seg <- data.frame(variable = factor(pc$variable, levels = ord),
                    group = grp_key,
                    xstart = abs(pc$unadjusted), xend = abs(pc$adjusted),
                    stringsAsFactors = FALSE)

  g <- ggplot2::ggplot(long, ggplot2::aes(x = value, y = variable)) +
    ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = xstart, xend = xend,
                                       y = variable, yend = variable),
                          inherit.aes = FALSE, colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(colour = adjustment, shape = adjustment),
                        size = 2.6) +
    ggplot2::scale_colour_manual(values = c(Unadjusted = "#bbbbbb",
                                            Adjusted = "#cc3311")) +
    ggplot2::labs(x = metric_lab, y = NULL, colour = NULL, shape = NULL,
                  title = paste0("Covariate balance (", b$sample, ")")) +
    ggplot2::theme_bw()
  if (grouped)
    g <- g + ggplot2::facet_wrap(~ group, labeller = ggplot2::label_both)
  g
}
