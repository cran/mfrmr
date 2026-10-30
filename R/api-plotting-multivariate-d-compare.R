#' @rdname mfrm_multivariate_d_compare
#' @details `plot()` method: Plot prespecified D-study plans and their differences
#'
#' Show G/Phi or SEM for every plan, including the reference. Use
#' `view = "differences"` for changes from the reference with approximate
#' pointwise intervals. In that view, the vertical zero line represents no
#' change; an interval crossing it does not
#' establish equivalence. Positive G/Phi differences or negative SEM differences
#' favor the comparison plan. The method requires normal random effects.
#' @param type `"coefficients"` for G/Phi or `"sem"` for SEM.
#' @param draw Draw the figure; `FALSE` returns its data only.
#' @param preset Plot style: `"standard"`, `"publication"`, `"compact"`, or
#'   `"monochrome"`.
#' @param view `"plans"` (default) shows point projections on their original
#'   scale, with a diamond for the reference. `"differences"` shows paired
#'   changes and their approximate intervals.
#' @details For `plot()`: The plan view connects supplied plans in their original row order.
#'   Both facet counts appear on the horizontal axis; lines are visual guides,
#'   not an interpolated response to changing one count. For conventional D-study
#'   curves that vary one count while holding the other constant, use
#'   [plot.mfrm_multivariate_d_study()]. The plan view has no sampling intervals:
#'   adding the reference estimate to a difference interval would not produce
#'   a confidence interval for the individual plan.
#'
#'   Plans and score weights must have been specified before inspecting
#'   the results. In the difference view, intervals are not simultaneous over
#'   plans or metrics. A point
#'   without an interval is retained as an open circle; an asterisk marks rows
#'   with unavailable intervals. Missing values are never replaced by zero.
#'   Each difference panel has its own horizontal scale.
#'   Original score units apply to SEM differences. [plot_data()] identifies the
#'   score/composite, its weights and the reference counts. Explanatory notes and
#'   a suggested title are returned rather than drawn. All plotted plans
#'   are future complete crossed plans, even when the source design is incomplete.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object. [plot_data()] extracts the exact
#'   comparison `table`, point-projection `series` (including the reference),
#'   `unavailable` rows for the selected view, `interval_unavailable` rows,
#'   `design_grid`, `reference`, `weights`, title and labels.
#'   Base graphics are supported; automatic ggplot conversion
#'   is not provided for this plot.
#' @seealso [mfrm_multivariate_d_compare()]
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_multivariate_d_comparison <- function(x, type = c("coefficients", "sem"),
                                                draw = TRUE, preset = "standard",
                                                view = c("plans", "differences"), ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  rlang::check_dots_empty()
  type <- match.arg(type)
  view <- match.arg(view)
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw)) stop("`draw` must be TRUE or FALSE.", call. = FALSE)
  style <- resolve_plot_preset(preset)
  metrics <- if (type == "coefficients") c("G", "Phi") else c("RelativeSEM", "AbsoluteSEM")
  table <- x$comparisons[x$comparisons$Metric %in% metrics, , drop = FALSE]
  plan_labels <- vapply(seq_len(nrow(x$design_grid)), function(i)
    paste(paste(names(x$design_grid), unlist(x$design_grid[i, ]), sep = " = "), collapse = ", "), character(1))
  table$Plan <- plan_labels[table$Scenario]
  series <- do.call(rbind, lapply(metrics, function(metric) {
    rows <- table[table$Metric == metric, , drop = FALSE]
    values <- data.frame(Scenario = c(x$reference, rows$Scenario), Metric = metric,
      Value = c(rows$ReferenceValue[1L], rows$Value))
    values[order(values$Scenario), , drop = FALSE]
  }))
  rownames(series) <- NULL
  series$Plan <- plan_labels[series$Scenario]
  series$IsReference <- series$Scenario == x$reference
  series$Status <- ifelse(is.finite(series$Value), "Available", "Point projection unavailable")
  selection <- if (x$kind == "Composite") paste0("Composite",
    if (x$score != "Composite") paste0(" ", x$score), ": ",
    paste(names(x$weights), format(x$weights, trim = TRUE), sep = " = ", collapse = ", ")) else
    paste("Score:", x$score)
  title <- paste(selection, "| Reference:", plan_labels[x$reference])
  subtitle <- paste0("Future complete crossed plans\n", if (view == "plans")
    "Point projections | Diamond: reference | No confidence intervals" else
    paste0(format(100 * x$level), "% approximate pointwise intervals | Normal random effects"))
  interval_unavailable <- table[table$Status != "Available", , drop = FALSE]
  unavailable <- if (view == "plans") series[series$Status != "Available", , drop = FALSE] else
    interval_unavailable
  note <- if (view == "plans") "Points are supplied plans in row order; lines are visual guides." else
    "Intervals crossing zero do not establish equivalence."
  if (nrow(unavailable)) note <- paste0(note, if (view == "plans")
    paste0(" Unavailable point projections: ", nrow(unavailable), ".") else
    paste0(" * = interval unavailable (", nrow(unavailable), "); open circles = point estimates only."))
  if (any(!x$component_diagnostics$PositiveSemidefinite) || !x$sampling_covariance_psd) {
    note <- paste(note, "Non-PSD covariance estimates: review diagnostics.")
  }
  out <- new_mfrm_plot_data("multivariate_d_comparison", list(table = table,
    series = series, view = view, unavailable = unavailable, interval_unavailable = interval_unavailable,
    design_grid = x$design_grid, reference = x$reference,
    weights = x$weights, title = title, subtitle = subtitle, note = note, metrics = metrics))
  if (!draw) return(invisible(out))
  apply_plot_preset(style)
  old <- graphics::par()[c("mfrow", "mar", "oma", "cex", "mex")]
  on.exit(graphics::par(old), add = TRUE)
  labels <- if (type == "coefficients") c("G: relative ordering", "Phi: absolute score levels") else
    c("Relative SEM: ordering", "Absolute SEM: score levels")
  # Keep margin line spacing compact so recorded panels also fit landscape devices.
  graphics::par(mfrow = c(2, 1), mex = .65,
    mar = if (view == "plans") c(6, 6, 2.5, 1) else c(5, 12, 2.5, 1),
    oma = c(.2, 0, .2, 0))
  for (j in seq_along(metrics)) {
    if (view == "plans") {
      tab <- series[series$Metric == metrics[j], , drop = FALSE]
      limits <- if (type == "coefficients") c(0, 1) else
        range(c(0, series$Value[is.finite(series$Value)]))
      if (diff(limits) == 0) limits <- c(0, 1)
      axis_labels <- vapply(seq_len(nrow(x$design_grid)), function(i)
        paste(paste(names(x$design_grid), unlist(x$design_grid[i, ]), sep = " = "), collapse = "\n"), character(1))
      axis_labels[x$reference] <- paste0(axis_labels[x$reference], "\n(reference)")
      plot(tab$Scenario, tab$Value, type = "n", ylim = limits, xaxt = "n", yaxt = "n",
        xlab = "", ylab = "", main = labels[j])
      mtext(if (type == "coefficients") "Dependability" else "SEM (score units)",
        side = 2, line = 4.2)
      axis(1, at = tab$Scenario, labels = axis_labels, cex.axis = .75, padj = 1)
      if (type == "coefficients") axis(2, at = seq(0, 1, .25), las = 1) else axis(2, las = 1)
      graphics::grid(nx = NA, ny = NULL, col = style$grid)
      lines(tab$Scenario, tab$Value, lwd = 2)
      points(tab$Scenario, tab$Value, pch = ifelse(tab$IsReference, 18, 16), cex = 1.2)
      if (!any(is.finite(tab$Value))) text(mean(tab$Scenario), mean(limits), "Point projections unavailable")
      next
    }
    tab <- table[table$Metric == metrics[j], , drop = FALSE]
    position <- rev(seq_len(nrow(tab)))
    extent <- unlist(tab[c("Difference", "Lower", "Upper")], use.names = FALSE)
    limits <- range(c(0, extent[is.finite(extent)]))
    if (diff(limits) == 0) limits <- limits + c(-.05, .05)
    direction <- if (type == "coefficients") "Higher favors comparison" else "Lower favors comparison"
    plot(limits, c(.5, nrow(tab) + .5), type = "n", yaxt = "n", ylab = "",
      xlab = paste("Difference from reference |", direction), main = labels[j])
    axis(2, at = position, labels = ifelse(tab$Status == "Available", tab$Plan,
      paste0(tab$Plan, " *")), las = 1, cex.axis = .8)
    abline(v = 0, lty = 2, col = "grey50")
    ok <- tab$Status == "Available"
    segments(tab$Lower[ok], position[ok], tab$Upper[ok], position[ok], lwd = 2)
    finite <- is.finite(tab$Difference)
    points(tab$Difference[finite], position[finite], pch = ifelse(ok[finite], 16, 1))
    if (!any(finite)) text(mean(limits), mean(position), "Point projections unavailable")
  }
  invisible(out)
}
