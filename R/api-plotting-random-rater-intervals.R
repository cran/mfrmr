#' @rdname mfrm_random_rater_intervals
#' @details `plot()` method: Compare bootstrap and ordinary random-rater intervals
#' @inheritParams plot.mfrm_testlet_scores
#' @seealso [plot.mfrm_testlet_scores()] for display controls and accessible
#'   output. Bootstrap comparisons retain the interval view.
#' @param comparison Show ordinary normal intervals from the same source fit.
#' @param draw `FALSE` returns plotted data without opening a device.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object with exact bounds, availability,
#'   settings and the source estimates. Use [plot_data()] for custom graphics.
#' @details For `plot()`: Arrows identify unbounded bootstrap endpoints. Their displayed ends
#'   are plot limits, not finite interval limits. Ordinary intervals never
#'   replace an unbounded bootstrap interval. These are pointwise prediction
#'   intervals for realized rater effects, not a rater-quality classification.
#' @export
plot.mfrm_random_rater_intervals <- function(x, level = x$settings$level,
    method = c("studentized", "error"), comparison = TRUE, draw = TRUE,
    sort = c("input", "estimate", "uncertainty"), palette = c("accessible", "mono"),
    title = NULL, caption = NULL, show_title = FALSE, show_notes = FALSE,
    show_labels = TRUE, reference = 0, text_scale = 1, point_size = 2.5, ...) {
  rlang::check_dots_empty(); method <- match.arg(method)
  if (any(!vapply(list(comparison, draw), function(z) is.logical(z) && length(z) == 1L && !is.na(z), logical(1)))) {
    stop("`comparison` and `draw` must be TRUE or FALSE.", call. = FALSE)
  }
  ci <- confint(x, level = level, method = method); source <- x$source$raters
  tab <- data.frame(Rater = source$Rater, Estimate = source$Estimate,
    Lower = ci[, 1L], Upper = ci[, 2L],
    OrdinaryLower = source$Estimate - stats::qnorm((1 - level) / 2, lower.tail = FALSE) * source$PredictionSE,
    OrdinaryUpper = source$Estimate + stats::qnorm((1 - level) / 2, lower.tail = FALSE) * source$PredictionSE,
    Unresolved = attr(ci, "availability")$Unresolved, row.names = NULL)
  note <- sprintf("Nominal %.0f%% pointwise | %d planned refits | arrows: unbounded", 100 * level, nrow(x$trials))
  if (comparison) note <- paste(note, "Solid: bootstrap; dashed: ordinary normal.", sep = "\n")
  mfrm_extended_estimate_plot(tab, tab$Rater, "random_rater_interval_bootstrap",
    "Observed-rater prediction intervals", "Rater severity (logits)", note, draw,
    c(x$settings, list(method = method, plotted_level = level, comparison = comparison)),
    sort = sort, palette = palette, title = title, caption = caption,
    show_title = show_title, show_notes = show_notes, show_labels = show_labels,
    reference = reference, text_scale = text_scale, point_size = point_size, extra = list(availability = attr(ci, "availability")))
}
