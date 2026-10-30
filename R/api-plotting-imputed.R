#' @rdname pool_mfrm_imputed
#' @details `plot()` method: Plot pooled facet estimates or prespecified contrasts
#' @param draw Logical. With `FALSE`, return the exact plot data without
#'   opening a graphics device.
#' @param preset A package plotting preset, such as `"standard"` or
#'   `"monochrome"`.
#' @param title,subtitle Optional text; `NULL` omits it. The default subtitle
#'   retains a notice when complete-data information is weak. Full cautions
#'   remain in the returned plot data even with custom or omitted text.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object containing the plotted `table`,
#'   exact contrast coefficients and inference settings. Use [plot_data()] for
#'   custom graphics. [as_ggplot()] uses the saved pointwise MI interval
#'   endpoints; it does not substitute normal intervals from the standard errors.
#'   Default and `component = "table"` retain the complete interval view.
#'   Target order, contrasts, degrees of freedom, confidence level and
#'   complete-data information cautions remain attached through [plot_data()].
#'   Numerical interval-availability reasons and saved inference notes are
#'   retained. Both base and ggplot views show a cross for a missing interval.
#'   Fixed targets have open diamonds without intervals. In ggplot conversion,
#'   infinite endpoints have arrows (their tips are display limits, not finite
#'   bounds); missing intervals have crosses at the saved estimate. Missing
#'   point estimates cause an error rather than being placed at zero.
#'   Use source `title`/`subtitle` arguments or `ggplot2::labs()` to change text;
#'   hiding headings does not remove the attached inference cautions.
#' @details For `plot()`: Plots use the stored pointwise multiple-imputation intervals and
#'   preserve target order. Open diamonds identify fixed targets without an
#'   inferential interval. A zero reference line does not define a practical
#'   importance threshold. The plot does not refit, select significant raters,
#'   adjust for multiplicity or turn severity differences into rater quality.
#' @seealso [pool_mfrm_imputed()]
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_pooled <- function(x, draw = TRUE, preset = "standard",
    title = paste("Pooled", x$settings$facet, "estimates and contrasts"),
    subtitle = sprintf("%d imputations | %.0f%% pointwise MI intervals",
      x$settings$imputations, 100 * x$settings$ci_level), ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  rlang::check_dots_empty()
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw)) {
    stop("`draw` must be TRUE or FALSE.", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  tab <- x$table
  if (missing(subtitle) && length(x$cautions)) subtitle <- paste(subtitle,
    "Weak information: review completed-data intervals.", sep = "\n")
  if (!all(vapply(list(title, subtitle), function(a)
      is.null(a) || (is.character(a) && length(a) == 1L && !is.na(a)), logical(1)))) {
    stop("Plot text must be NULL or one nonmissing character string.", call. = FALSE)
  }
  out <- new_mfrm_plot_data("pooled_facet_intervals", list(table = tab,
    contrasts = x$contrasts, settings = x$settings, title = title,
    cautions = x$cautions, notes = x$notes, information_review = x$information_review,
    subtitle = subtitle, xlab = "Estimate (logits)", preset = style$name,
    caption = "Open diamonds: fixed targets without intervals; crosses: unavailable intervals; arrows in ggplot: infinite endpoints, not finite bounds."))
  if (!draw) return(invisible(out))
  apply_plot_preset(style)
  old <- graphics::par(mar = c(4.5, 10, 1.5, 2))
  on.exit(graphics::par(old), add = TRUE)
  y <- rev(seq_len(nrow(tab)))
  limits <- range(c(0, tab$Estimate, tab$Lower, tab$Upper), finite = TRUE)
  if (diff(limits) == 0) limits <- limits + c(-0.1, 0.1)
  fixed <- tab$Status == "fixed"
  plot(tab$Estimate, y, xlim = limits, ylim = c(0.5, nrow(tab) + 0.5),
    pch = ifelse(fixed, 5, ifelse(is.na(tab$Lower) | is.na(tab$Upper), 4, 16)),
    col = style$accent_primary, yaxt = "n",
    xlab = "Estimate (logits)", ylab = "", main = "")
  abline(v = 0, col = style$grid, lty = 2)
  segments(tab$Lower[!fixed], y[!fixed], tab$Upper[!fixed], y[!fixed],
    col = style$accent_primary, lwd = 2)
  axis(2, at = y, labels = tab$Target, las = 1)
  invisible(out)
}

.mfrmr_gg_pooled <- function(x) {
  .require_mfrmr_ggplot2()
  d <- x$data
  tab <- d$table
  required <- c("Target", "Estimate", "Lower", "Upper", "DF", "Status")
  if (!is.data.frame(tab) || !nrow(tab) || !all(required %in% names(tab)) ||
      any(!is.finite(tab$Estimate)) || anyNA(tab$Target) || anyNA(tab$Status) ||
      !all(tab$Status %in% c("pooled", "fixed"))) {
    stop("Pooled plot data lack finite estimates, saved interval fields or supported target statuses. Recreate with plot(pooled, draw = FALSE).", call. = FALSE)
  }
  free <- tab$Status == "pooled"
  complete <- free & !is.na(tab$Lower) & !is.na(tab$Upper)
  if (any(tab$Lower[complete] > tab$Upper[complete]) ||
      any(tab$Lower[complete] == Inf | tab$Upper[complete] == -Inf)) {
    stop("Saved pooled interval endpoints are not ordered lower to upper.", call. = FALSE)
  }
  tab$Y <- rev(seq_len(nrow(tab)))
  tab$Display <- ifelse(!free, "Fixed", ifelse(!complete, "Unavailable interval",
    ifelse(is.infinite(tab$Lower) | is.infinite(tab$Upper), "Unbounded interval", "Pooled interval")))
  style <- resolve_plot_preset(d$preset %||% "standard")
  endpoints <- c(0, tab$Estimate, tab$Lower[complete], tab$Upper[complete])
  limits <- range(endpoints[is.finite(endpoints)])
  padding <- if (diff(limits) > 0) .08 * diff(limits) else .1
  limits <- limits + c(-padding, padding)
  intervals <- tab[complete, , drop = FALSE]
  intervals$Left <- pmax(limits[1], intervals$Lower)
  intervals$Right <- pmin(limits[2], intervals$Upper)
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$Estimate, y = .data$Y)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = style$grid) +
    ggplot2::geom_segment(data = intervals, ggplot2::aes(x = .data$Left,
      xend = .data$Right, yend = .data$Y), colour = style$accent_primary, linewidth = .7) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$Display), colour = style$accent_primary, size = 2.8) +
    ggplot2::scale_shape_manual(values = c("Pooled interval" = 16, Fixed = 5,
      "Unavailable interval" = 4, "Unbounded interval" = 1), name = NULL) +
    ggplot2::scale_y_continuous(breaks = rev(seq_len(nrow(tab))), labels = tab$Target,
      limits = c(.5, nrow(tab) + .5)) + .mfrmr_gg_theme()
  # Arrow tips are drawing limits, never substituted confidence endpoints.
  for (side in c("Lower", "Upper")) {
    arrows <- intervals[is.infinite(intervals[[side]]), , drop = FALSE]
    if (!nrow(arrows)) next
    left <- side == "Lower"
    arrows$End <- limits[if (left) 1 else 2]
    arrows$Start <- arrows$End + if (left) padding / 2 else -padding / 2
    p <- p + ggplot2::geom_segment(data = arrows,
      ggplot2::aes(x = .data$Start, xend = .data$End, yend = .data$Y),
      colour = style$accent_primary, linewidth = .7,
      arrow = ggplot2::arrow(length = ggplot2::unit(.08, "inches")))
  }
  d$caption <- paste(c("Saved pointwise MI intervals; zero is not a practical-importance threshold.",
    if (any(!free)) "Open diamonds: fixed targets without inferential intervals.",
    if (any(tab$Display == "Unavailable interval")) "Crosses: unavailable intervals; inspect the saved table.",
    if (any(tab$Display == "Unbounded interval")) "Arrows: infinite endpoints, not finite bounds."), collapse = "\n")
  p <- .mfrmr_gg_labs(p, d, x = d$xlab, y = NULL) +
    ggplot2::labs(alt = "Pooled facet estimates or prespecified contrasts in saved target order, with stored pointwise multiple-imputation t intervals. Open diamonds mark fixed targets, crosses mark unavailable intervals, and arrows mark infinite endpoints. No multiplicity adjustment or rater-quality classification is added.")
  attr(p, "mfrmr_plot_data") <- x
  p
}
