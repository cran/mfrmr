#' @rdname read_mfrm_posterior
#' @details `plot()` method: Plot saved posterior estimates and credible intervals
#' @param x Output from [read_mfrm_posterior()].
#'
#' Also for `print()`: A saved `mfrm_posterior_summary`.
#' @param quantity One of `"location"` (default), `"ability"`, `"step_offset"`,
#'   `"slope"`, `"log_slope"`, `"population_mean"` or `"population_sd"`.
#'   Population quantities require an estimated-population export.
#'   Quantities with different units are not mixed.
#' @param facet Optional facet/owner names to retain; `NULL` retains all for
#'   the selected quantity. Step rows use `facet = "Step owner"` and retain
#'   their owner labels and category transitions in the plotted table.
#' @param draw Logical; `FALSE` returns plot data without opening a device.
#' @param reference Optional finite vertical reference; `NULL` omits it.
#'   Zero for locations and one for relative slopes are references, not
#'   standards of performance or rater quality.
#' @param preset A package plot preset, such as `"standard"` or `"monochrome"`.
#' @param show_legend Logical; show interval/fixed-value symbols.
#' @param title,subtitle,caption Optional returned text. Defaults omit it.
#'   Text is not drawn by default; retrieve it with [plot_data()].
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` with the exact long `table`, source
#'   settings, parameter identities, interval statuses and interpretation notes.
#'   [as_ggplot()] returns an editable ggplot; [plot_data()] supplies the table
#'   for custom graphics. No posterior draws, densities or raincloud samples
#'   are reconstructed from marginal summaries.
#' @details For `plot()`: These are saved equal-tailed 95% posterior credible intervals,
#'   not frequentist confidence intervals. Dashed intervals and open circles
#'   mark estimates whose sampling screen requires review; solid intervals
#'   mean that no sampling flag was detected, not that inference is validated.
#'   Open diamonds mark fixed anchors/constraint-determined values without an
#'   interval. Crosses mark estimated points without an available interval.
#'   Missing points retain their labels and are never replaced by zero.
#'   Original row order, labels and endpoints are preserved. The plot never
#'   refits, changes the level, selects significant effects or validates priors.
#' @seealso [read_mfrm_posterior()], [mfrm_report()]
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @examples posterior <- read_mfrm_posterior(system.file(
#'   "extdata", "posterior-example", package = "mfrmr"))
#' p <- plot(posterior, quantity = "location", draw = FALSE)
#' plot_data(p, "table")
#' if (requireNamespace("ggplot2", quietly = TRUE)) as_ggplot(p)
#' @export
plot.mfrm_posterior_summary <- function(x,
    quantity = c("location", "ability", "step_offset", "slope", "log_slope",
                 "population_mean", "population_sd"),
    facet = NULL, draw = TRUE, reference = NULL, preset = "standard", show_legend = TRUE,
    title = NULL, subtitle = NULL, caption = NULL, ...) {
  rlang::check_dots_empty()
  .mfrm_validate_posterior(x)
  quantity <- match.arg(quantity)
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  style <- resolve_plot_preset(preset)
  if (!all(vapply(list(draw, show_legend), function(v) is.logical(v) && length(v) == 1L && !is.na(v), logical(1))))
    stop("`draw` and `show_legend` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(reference) && (!is.numeric(reference) || length(reference) != 1L || !is.finite(reference)))
    stop("`reference` must be NULL or one finite number.", call. = FALSE)
  if (!all(vapply(list(title, subtitle, caption), function(v) is.null(v) ||
      (is.character(v) && length(v) == 1L && !is.na(v)), logical(1))))
    stop("Plot text must be NULL or one nonmissing string.", call. = FALSE)
  tab <- x$parameters[x$parameters$quantity == quantity, , drop = FALSE]
  if (!is.null(facet)) {
    if (!is.character(facet) || !length(facet) || anyNA(facet) || !all(facet %in% tab$facet))
      stop("`facet` must select existing owners for the selected quantity.", call. = FALSE)
    tab <- tab[tab$facet %in% facet, , drop = FALSE]
  }
  if (!nrow(tab)) stop("No saved parameters match this plot selection.", call. = FALSE)
  tab$Display <- ifelse(tab$value_role %in% c("fixed_anchor", "fixed_by_constraint"), "Fixed value",
    ifelse(!tab$interval_computed, "Interval unavailable", ifelse(
      tab$sampling_review == "sampling_review_required", "Sampling review required", "No sampling flag")))
  tab$Label <- paste(tab$facet, tab$level, sep = ": ")
  steps <- !is.na(tab$transition)
  tab$Label[steps] <- paste(tab$Label[steps], tab$transition[steps], sep = " | ")
  tab$Label[!is.finite(tab$mean)] <- paste(tab$Label[!is.finite(tab$mean)], "(not estimated)")
  tab$Y <- rev(seq_len(nrow(tab)))
  xlab <- switch(quantity, location = "Location (fixed N(0,1) ability units)",
    ability = "Ability (fixed N(0,1) units)", step_offset = "Step offset (fixed N(0,1) ability units)",
    slope = "Component discrimination", log_slope = "Log component discrimination",
    population_mean = "Population mean (native relative-slope metric)",
    population_sd = "Population SD (native relative-slope metric)")
  if (!is.null(x$population)) xlab <- sub("fixed N(0,1) ability units",
    "native relative-slope metric", sub("fixed N(0,1) units", "native relative-slope metric", xlab, fixed = TRUE), fixed = TRUE)
  out <- new_mfrm_plot_data("posterior_intervals", list(table = tab, source = x$source,
    checks = x$checks, decision = x$decision, notes = x$notes, level = x$level,
    interval_type = x$interval, title = title, subtitle = subtitle, caption = caption,
    xlab = xlab, reference = reference, preset = style$name, show_legend = show_legend))
  if (!draw) return(invisible(out))
  apply_plot_preset(style)
  old <- graphics::par(mar = c(4.5, 10, if (show_legend) 4 else 1.5, 2)); on.exit(graphics::par(old), add = TRUE)
  values <- c(reference, tab$mean, tab$lower, tab$upper)
  limits <- if (any(is.finite(values))) range(values, finite = TRUE) else c(-1, 1)
  if (diff(limits) == 0) limits <- limits + c(-.1, .1)
  graphics::plot(tab$mean, tab$Y, type = "n", xlim = limits, ylim = c(.5, nrow(tab) + .5),
    yaxt = "n", xlab = xlab, ylab = "", main = "")
  if (!is.null(reference)) graphics::abline(v = reference, lty = 3, col = style$grid)
  graphics::segments(tab$lower, tab$Y, tab$upper, tab$Y, col = style$accent_primary,
    lty = ifelse(tab$Display == "Sampling review required", 2, 1))
  symbols <- c("Fixed value" = 5, "Interval unavailable" = 4,
    "Sampling review required" = 1, "No sampling flag" = 16)
  graphics::points(tab$mean, tab$Y, pch = unname(symbols[tab$Display]), col = style$accent_primary)
  graphics::axis(2, at = tab$Y, labels = tab$Label, las = 1)
  if (show_legend) {
    used <- unique(tab$Display)
    graphics::legend(x = mean(graphics::par("usr")[1:2]),
      y = graphics::par("usr")[4] + .6 * graphics::strheight("M"),
      xjust = .5, yjust = 0, horiz = TRUE, xpd = NA,
      title = "95% credible interval", legend = used, pch = unname(symbols[used]),
      col = style$accent_primary, bty = "n", cex = .75)
  }
  invisible(out)
}

.mfrm_gg_posterior_intervals <- function(x) {
  .require_mfrmr_ggplot2()
  d <- x$data; tab <- d$table
  style <- resolve_plot_preset(d$preset)
  finite <- tab[is.finite(tab$lower) & is.finite(tab$upper), , drop = FALSE]
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$mean, y = .data$Y)) +
    ggplot2::geom_segment(data = finite, ggplot2::aes(x = .data$lower, xend = .data$upper,
      yend = .data$Y, linetype = .data$Display), color = style$accent_primary, linewidth = .8) +
    ggplot2::geom_point(data = tab[is.finite(tab$mean), , drop = FALSE],
      ggplot2::aes(shape = .data$Display), color = style$accent_primary, size = 2.8) +
    ggplot2::scale_shape_manual(limits = unique(tab$Display), values = c("Fixed value" = 5, "Interval unavailable" = 4,
      "Sampling review required" = 1, "No sampling flag" = 16)) +
    ggplot2::scale_linetype_manual(limits = unique(tab$Display), values = c("Fixed value" = "blank", "Interval unavailable" = "blank",
      "Sampling review required" = "dashed", "No sampling flag" = "solid")) +
    ggplot2::scale_y_continuous(breaks = tab$Y, labels = tab$Label, limits = c(.5, nrow(tab) + .5)) +
    ggplot2::labs(shape = "95% credible interval", linetype = "95% credible interval",
      alt = "Saved posterior means and 95% credible intervals. Diamonds indicate fixed values; dashed lines and open circles indicate sampling review. Missing estimates are not replaced by zero.") + .mfrmr_gg_theme()
  if (!is.null(d$reference)) p <- p + ggplot2::geom_vline(xintercept = d$reference, linetype = "dotted", color = "grey60")
  if (!any(is.finite(c(d$reference, tab$mean, tab$lower, tab$upper))))
    p <- p + ggplot2::scale_x_continuous(limits = c(-1, 1))
  if (!d$show_legend) p <- p + ggplot2::theme(legend.position = "none")
  p <- .mfrmr_gg_labs(p, d, x = d$xlab)
  attr(p, "mfrmr_plot_data") <- x
  p
}
