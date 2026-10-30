#' @rdname mfrm_population_intervals
#' @details `plot()` method: Plot saved latent-population coefficients or contrasts
#'
#' Read point estimates and their experimental pointwise intervals from an
#' existing population analysis, without fitting or recalculating uncertainty.
#' @param type Show `"coefficients"` (default, including the intercept) or
#'   the prespecified `"contrasts"` stored in `x`. Residual variance has different
#'   units and remains in `summary(x)` rather than sharing this axis.
#' @param draw Logical; `FALSE` returns the plot data without drawing.
#' @param preset A package plot preset, including `"publication"` and
#'   `"monochrome"`. The session default is used when omitted.
#' @param reference Optional finite vertical reference, default zero; `NULL`
#'   omits it. This is not a practical-importance or automatic decision threshold.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object. `table` preserves selected
#'   source rows and adds `PlotStatus`, `Label` and `Y`; `intervals` retains the
#'   complete original table. `population_spec` preserves the formula, trained
#'   terms, factor levels and coding; `population_coding` is a readable coding
#'   table. `contrast_weights`, source identity, settings, numerical checks and
#'   interpretation notes accompany the figure. Use [plot_data()] to inspect
#'   these and [as_ggplot()] for an editable figure with the same saved values.
#' @details For `plot()`: Rows retain their saved order. Filled circles and lines show numerically
#' available intervals. Crosses with a `CI unavailable` label retain finite
#' points without eligible finite bounds. A missing point keeps its labeled row
#' without a marker; it is never plotted at zero. Titles and explanatory notes
#' are returned, not drawn. No significance colors or automatic classifications
#' are added. Plotting does not requalify a saved numerical result.
#'
#' Coefficients describe the conditional population mean on the fitted native
#' ability scale. A continuous coefficient is per unit of its fitted predictor;
#' the intercept depends on predictor zero values and factor coding. Read
#' interactions, transformations and categorical comparisons using the saved
#' coding, not an assumed reference group. Coefficients on different predictor
#' scales are not directly comparable measures of importance. Higher ability
#' predicts higher coded ratings in these RSM/PCM models; whether that is better
#' performance depends on the rubric. Coefficients are not standardized effects.
#'
#' Contrasts use the stored linear weights and may describe a conditional mean
#' or a difference between means. They are not individual ability differences,
#' causal effects, simultaneous intervals or prediction intervals. Numerical
#' availability does not establish nominal coverage: the source's experimental
#' status, failed checks and qualifications remain in `notes` and the tables.
#' A wide interval remains visible. Missing bounds are not infinite endpoints.
#' When source-anchor uncertainty was requested, the saved intervals include
#' that contribution. The returned subtitle, settings and `anchor_uncertainty`
#' distinguish it from intervals conditional on exact supplied anchors;
#' conditional bounds remain in the returned table.
#'
#' For an executable bundled-data example, run the example on
#' [mfrm_population_intervals()], which creates `ci` and both views. Then use
#' `p <- plot(ci, type = "contrasts", draw = FALSE)`,
#' `plot_data(p, "contrast_weights")` and `plot_data(p, "notes")` to read it.
#' `as_ggplot(ci, type = "contrasts")` converts the same saved result.
#' @seealso [mfrm_population_intervals()], [mfrmr_interval_guide()],
#'   [plot.mfrm_facet_intervals()]
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_population_intervals <- function(x, type = c("coefficients", "contrasts"),
                                          draw = TRUE, preset = "standard",
                                          reference = 0, ...) {
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  style <- resolve_plot_preset(preset)
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw))
    stop("`draw` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(reference) && (!is.numeric(reference) || is.complex(reference) ||
      length(reference) != 1L || !is.finite(reference)))
    stop("`reference` must be NULL or one finite number.", call. = FALSE)
  required <- c("Target", "Type", "Term", "Estimate", "Lower", "Upper",
    "CIEligible", "Status", "InferenceReview")
  if (!is.data.frame(x$table) || !all(required %in% names(x$table)))
    stop("The saved population result needs its interval table and qualification columns.", call. = FALSE)
  target_type <- if (type == "coefficients") "coefficient" else "coefficient_contrast"
  tab <- x$table[x$table$Type %in% target_type, , drop = FALSE]
  if (!nrow(tab)) stop("No saved ", type, " to plot. For contrasts, supply weights to mfrm_population_intervals() first.", call. = FALSE)
  available <- tab$CIEligible %in% TRUE & tab$Status %in% "available" &
    is.finite(tab$Estimate) & is.finite(tab$Lower) & is.finite(tab$Upper) &
    tab$Lower <= tab$Upper
  tab$PlotStatus <- ifelse(!is.finite(tab$Estimate), "Estimate unavailable",
    ifelse(available, "Interval available", "CI unavailable"))
  labels <- paste0(tab$Term, ifelse(tab$PlotStatus == "Interval available", "",
    paste0(" [", tab$PlotStatus, "]")))
  tab$Label <- vapply(labels, function(z) paste(strwrap(z, width = 38), collapse = "\n"), character(1))
  heights <- lengths(strsplit(tab$Label, "\n", fixed = TRUE))
  tab$Y <- rev(cumsum(rev(pmax(1, heights)))) - (pmax(1, heights) - 1) / 2
  limits <- c(reference, tab$Estimate, tab$Lower[available], tab$Upper[available])
  limits <- if (any(is.finite(limits))) range(limits, finite = TRUE) else c(-1, 1)
  if (diff(limits) == 0) limits <- limits + c(-.1, .1)
  population <- x$fit$population %||% list()
  formula <- if (is.null(population$formula)) "not recorded" else
    paste(deparse(population$formula), collapse = " ")
  cautions <- unique(c(x$cautions,
    "Experimental pointwise population intervals: numerical availability does not establish sampling coverage. No multiplicity adjustment, causal inference or individual prediction interval is supplied.",
    "The plot uses saved values and numerical checks without recomputing or requalifying them. Missing bounds are unavailable, not infinite; crosses retain finite estimates without eligible intervals.",
    "Native ability units are not population-SD-standardized. Predictor units, transformations, interactions and factor coding determine coefficient meaning; magnitudes are not comparable measures of predictor importance.",
    "Higher conditional population ability predicts higher coded ratings; whether higher ratings are favorable depends on the rubric. The zero line is a reference, not a practical-importance threshold.",
    "Contrasts use the saved weights and can represent conditional means or differences; the weights determine the interpretation.",
    paste0("Population formula: ", formula, ". Inspect population_spec and population_coding for fitted coding."),
    "Residual variance is retained in intervals but not drawn on the coefficient/contrast axis."))
  out <- new_mfrm_plot_data("population_intervals", list(table = tab,
    intervals = x$table, type = type, settings = x$settings,
    anchor_uncertainty = x$anchor_uncertainty,
    population_spec = population[intersect(c("formula", "terms", "xlevels", "contrasts", "design_columns"), names(population))],
    population_coding = population_coding_summary_table(population),
    contrast_weights = x$contrasts, source = x$source,
    checks = x$checks, numerical_checks = x$numerical_checks,
    information_review = x$information_review,
    comparison_information_review = x$comparison_information_review,
    cautions = cautions, reference = reference, preset = style$name,
    x_range = limits, y_range = c(.5, sum(pmax(1, heights)) + .5),
    title = paste("Conditional population", type),
    subtitle = paste0(format(100 * x$settings$level, trim = TRUE),
      "% experimental pointwise intervals",
      if (identical(x$settings$uncertainty, "independent_anchor_delta"))
        " | Independent anchor-source uncertainty included"),
    xlab = if (type == "coefficients") "Coefficient (native ability scale)" else
      "Contrast (native ability scale)"))
  if (!draw) return(invisible(out))
  apply_plot_preset(style)
  old <- graphics::par(mar = c(4.5, 10, 1.5, 2))
  on.exit(graphics::par(old), add = TRUE)
  plot(NA_real_, NA_real_, type = "n", xlim = limits,
    ylim = out$data$y_range, yaxt = "n",
    xlab = out$data$xlab, ylab = "", main = "")
  if (!is.null(reference)) abline(v = reference, col = style$grid, lty = 2)
  segments(tab$Lower[available], tab$Y[available],
    tab$Upper[available], tab$Y[available], col = style$accent_primary, lwd = 2)
  finite <- is.finite(tab$Estimate)
  points(tab$Estimate[finite], tab$Y[finite],
    pch = ifelse(available[finite], 16, 4), col = style$accent_primary)
  axis(2, at = tab$Y, labels = tab$Label, las = 1, cex.axis = .8)
  invisible(out)
}

.mfrm_gg_population_intervals <- function(x) {
  .require_mfrmr_ggplot2()
  d <- x$data; tab <- d$table
  style <- resolve_plot_preset(d$preset)
  available <- tab[tab$PlotStatus == "Interval available", , drop = FALSE]
  finite <- tab[is.finite(tab$Estimate), , drop = FALSE]
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = available, ggplot2::aes(x = .data$Lower,
      xend = .data$Upper, y = .data$Y, yend = .data$Y),
      colour = style$accent_primary, linewidth = .8) +
    ggplot2::geom_point(data = finite, ggplot2::aes(x = .data$Estimate,
      y = .data$Y, shape = .data$PlotStatus), colour = style$accent_primary, size = 2.8) +
    ggplot2::scale_shape_manual(values = c("Interval available" = 16, "CI unavailable" = 4), guide = "none") +
    ggplot2::scale_y_continuous(breaks = tab$Y, labels = tab$Label) +
    ggplot2::coord_cartesian(xlim = d$x_range, ylim = d$y_range) +
    ggplot2::labs(alt = "Conditional population estimates and experimental pointwise intervals. Crosses and row labels identify unavailable intervals; missing estimates retain a labeled row without a marker. Numerical availability does not establish coverage.") +
    .mfrmr_gg_theme()
  if (!is.null(d$reference)) p <- p + ggplot2::geom_vline(xintercept = d$reference,
    linetype = "dashed", colour = "grey60")
  p <- .mfrmr_gg_labs(p, d, x = d$xlab, y = NULL)
  attr(p, "mfrmr_plot_data") <- x
  p
}
