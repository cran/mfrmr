#' Compare D-study confidence interval methods
#'
#' Compare future plans on one shared coefficient or SEM axis. Each plan has adjacent
#' method rows, distinguished by color, point shape and line type. Method
#' agreement in a single dataset is not evidence of coverage.
#' @param x A result from [confint.mfrm_multivariate_d_study()].
#' @param draw Logical; `FALSE` returns the plot data without opening a device.
#' @param preset Package plot style, including `"standard"` (color) and
#'   `"monochrome"`. Methods also have distinct point shapes.
#' @param ncol Positive integer number of columns for `layout = "panels"`.
#'   Unused by the default combined view.
#' @param layout `"combined"` (default) groups methods within each plan on
#'   one axis, with a common method legend. `"panels"` retains one panel per
#'   future plan, with at most `ncol` columns.
#' @param ... Reserved; additional arguments are rejected.
#' @return Invisibly, an `mfrm_plot_data` object retaining the exact `table`,
#'   `unavailable` rows, source `intervals`, `design_grid`, `weights` and `notes`.
#'   Use [plot_data()] to inspect them and [as_ggplot()] on the returned object
#'   for an editable plot. No title or explanatory footnote is drawn.
#' @details All methods concern the same score/composite, metric and
#'   source data. A filled point has an interval; an open point is an available
#'   estimate without an interval. Unavailable intervals are explicitly labeled
#'   in the combined view, or marked by an asterisk in the panel view, even
#'   if their point estimates are missing. An interval
#'   can exist without an admissible ANOVA/MINQUE point estimate. Missing values
#'   remain missing; no substitute interval or zero is drawn.
#'   The shared axis includes `[0,1]` for coefficients and `[-1,1]` for G/Phi
#'   differences, expanding to retain all finite estimates and accepted bounds.
#'   Delta normal-reference bounds outside `[0,1]` remain visible; they flag
#'   limitations of that approximation, not a larger population parameter space.
#'   SEM and SEM-difference axes instead span the estimates and intervals in
#'   score units, including zero. Negative SEM differences favor the comparison
#'   plan. An SEM interval is uncertainty in measurement error, not an interval
#'   for a Person's score. Chi-square and MLS labels retain the calculation
#'   used for each plan; their assumptions and coverage claims differ.
#'   Method limitations, confidence level, normal-effect assumptions and
#'   negative components remain in the returned notes and source object.
#'   Plotting neither refits a model nor recalculates an interval. Save the
#'   complete interval or plot object with `saveRDS()` for later use.
#'   Nested plan labels show child counts per parent, rather than total child
#'   pool sizes; the reference plan uses the same count convention.
#'   A mean-square-region difference is a conservative enclosure from shared
#'   source uncertainty, not a tight joint image or a paired-delta interval.
#'   Its returned source retains the individual projections and reference.
#'   REML profile and residual-contrast rows retain the original source points;
#'   the separate `IntervalCenter` and method diagnostics are available in the source
#'   interval object. A point can lie outside this separately calculated interval.
#'   For GENOVA-adjusted results, points use the selected component treatment;
#'   intervals retain their original-source population basis. `IntervalEstimate`
#'   and `IntervalBasis` in the source table distinguish them. Drawing the two
#'   together does not make the original sampling SE an adjusted-estimator SE.
#'   Unresolved numerical candidates are retained in the source object but
#'   are not drawn as accepted confidence bounds.
#'   For score and MLS intervals, an empty confidence set is labeled
#'   `"Empty set"`; it is distinct from an estimated coefficient of zero.
#'   MLS rows with `FormulaComplete = FALSE` retain candidate values where
#'   a squared width is undefined. Their conservative bounds are drawn, while
#'   the saved interval object and notes retain the reason and affected ranges.
#' @seealso [confint.mfrm_multivariate_d_study()]
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_dstudy_intervals", run.dontrun = TRUE)
#' \dontrun{
#' tasks <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
#' g <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
#' d <- mfrm_d_study(g, data.frame(Tasks = c(6, 12)))
#' ci <- confint(d, assumption = "normal")
#' p <- plot(ci, draw = FALSE)
#' plot_data(p, "notes")
#' if (requireNamespace("ggplot2", quietly = TRUE)) as_ggplot(p)
#' }
#' @export
plot.mfrm_dstudy_intervals <- function(x, draw = TRUE, preset = "standard",
                                      ncol = 2L, layout = c("combined", "panels"), ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  rlang::check_dots_empty()
  layout <- match.arg(layout)
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw))
    stop("`draw` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(ncol) || is.complex(ncol) || length(ncol) != 1L ||
      !is.finite(ncol) || ncol < 1 || ncol != floor(ncol))
    stop("`ncol` must be a positive integer.", call. = FALSE)
  if (!nrow(x)) stop("No interval rows to plot.", call. = FALSE)
  style <- resolve_plot_preset(preset)
  tab <- as.data.frame(x); grid <- attr(x, "design_grid")
  count_labels <- names(grid)
  design <- attr(x, "source_summary")$design
  if (!is.null(design$nesting)) {
    child_column <- design$count_columns[[names(design$nesting)]]
    count_labels[count_labels == child_column] <- paste(child_column, "per", unname(design$nesting))
  }
  labels <- vapply(seq_len(nrow(grid)), function(i)
    paste(paste(count_labels, unlist(grid[i, ]), sep = " = "),
      collapse = if (is.null(design$nesting)) ", " else "\n"), character(1))
  if (anyDuplicated(labels)) labels <- paste0(seq_along(labels), ": ", labels)
  plans <- unique(tab$Scenario); methods <- unique(tab$Method)
  tab$Plan <- factor(labels[tab$Scenario], levels = unique(labels[plans]))
  tab$Row <- length(methods) + 1L - match(tab$Method, methods)
  tab$Available <- is.finite(tab$Lower) & is.finite(tab$Upper)
  tab$MethodLabel <- paste0(tab$Method, ifelse(tab$Available, "", " *"))
  palette <- rep(c(style$accent_primary, style$accent_secondary), length.out = length(methods))
  if (style$name == "monochrome") palette[] <- "black"
  tab$Colour <- palette[match(tab$Method, methods)]
  tab$Shape <- ifelse(tab$Available, 15L + match(tab$Method, methods), match(tab$Method, methods))
  tab$LineType <- match(tab$Method, methods)
  centers <- length(plans) + 1L - seq_along(plans)
  offsets <- if (length(methods) == 1L) 0 else seq(.16, -.16, length.out = length(methods))
  tab$Y <- centers[match(tab$Scenario, plans)] + offsets[match(tab$Method, methods)]
  sem <- any(grepl("SEM", tab$Metric, fixed = TRUE))
  domain <- if (sem) 0 else if (any(tab$Metric %in% c("G difference", "Phi difference"))) c(-1, 1) else c(0, 1)
  values <- c(domain, tab$Estimate, tab$Lower[tab$Available], tab$Upper[tab$Available])
  xrange <- range(values[is.finite(values)])
  if (diff(xrange) == 0) xrange <- c(0, 1)
  # Place failure labels away from any available estimate, never as a point at zero.
  tab$StatusHjust <- as.numeric(!is.finite(tab$Estimate) | tab$Estimate <= mean(xrange))
  tab$StatusX <- ifelse(tab$StatusHjust == 1, xrange[2] - .02*diff(xrange),
    xrange[1] + .02*diff(xrange))
  tab$StatusLabel <- ifelse(is.finite(tab$Estimate), "CI unavailable",
    paste0(tab$Method, ": unavailable"))
  if ("Empty" %in% names(tab))
    tab$StatusLabel[!is.na(tab$Empty) & tab$Empty & !tab$Available] <- "Empty set"
  notes <- c(attr(x, "notes"),
    if (!is.null(attr(x, "reference"))) paste0("Reference plan: ",
      labels[attr(x, "reference")$scenario], "; differences are comparison minus reference."),
    paste0(format(100*unique(x$Level), trim = TRUE), "% intervals; metric: ",
      paste(unique(x$Metric), collapse = ", "), "; score: ", paste(unique(x$Score), collapse = ", "), "."),
    "Filled points: estimates with intervals. Open points: estimates without intervals. Unavailable intervals are labeled in the combined view or marked * in panels; missing estimates are never placed at zero. Inspect unavailable rows for reasons.",
    if (any(x$NegativeComponent)) "Negative ANOVA/MINQUE component estimates are present; review method restrictions.")
  out <- new_mfrm_plot_data("dstudy_intervals", list(table = tab,
    intervals = x, unavailable = tab[!tab$Available, , drop = FALSE],
    design_grid = grid, weights = attr(x, "weights"),
    notes = data.frame(Type = "Interpretation", Text = notes),
    preset = style$name, ncol = min(ncol, length(plans)), layout = layout,
    plans = data.frame(Scenario = plans, Label = labels[plans], Y = centers),
    methods = methods, palette = palette, xrange = xrange,
    xlab = if (sem) paste(sub("SEM", " SEM", unique(x$Metric), fixed = TRUE), "(score units)") else unique(x$Metric)))
  if (!draw) return(invisible(out))
  apply_plot_preset(style)
  old <- graphics::par()[c("mfrow", "mar", "oma", "cex", "mex")]
  on.exit(graphics::par(old), add = TRUE)
  combined <- layout == "combined"
  graphics::par(mfrow = if (combined) c(1, 1) else
      c(ceiling(length(plans)/out$data$ncol), out$data$ncol),
    mex = .7, mar = c(if (combined) 5 else 4, if (combined) 11 else 10,
      if (combined && length(methods) > 2L) 1.5*length(methods)+1 else if (combined) 3.5 else 2, 1),
    oma = c(0, 0, 0, 0))
  if (combined) {
    margin <- graphics::par("mar")
    margin[2] <- max(margin[2],
      max(graphics::strwidth(labels[plans], units = "inches", cex = .8))/
        (graphics::par("csi")*graphics::par("mex"))+2)
    graphics::par(mar = margin)
  }
  groups <- if (combined) list(plans) else as.list(plans)
  for (plan in groups) {
    z <- tab[tab$Scenario %in% plan, , drop = FALSE]
    y <- if (combined) z$Y else z$Row
    plot(NA_real_, NA_real_, xlim = xrange,
      ylim = c(.5, if (combined) length(plans)+.5 else length(methods)+.5),
      xlab = out$data$xlab, ylab = "", yaxt = "n", main = "")
    if (combined) {
      if (length(plans) > 1L) abline(h = centers[-1L]+.5, col = "grey90")
      axis(2, at = centers, labels = labels[plans], las = 1, cex.axis = .8, tick = FALSE)
      legend(x = mean(graphics::par("usr")[1:2]),
        y = graphics::par("usr")[4] + graphics::strheight("M")*
          if (length(methods) > 2L) .25 else 1.5,
        xjust = .5, yjust = if (length(methods) > 2L) 0 else .5,
        xpd = NA, bty = "n", horiz = length(methods) <= 2L,
        legend = methods, col = palette, pch = 15L+seq_along(methods),
        lty = seq_along(methods), lwd = 2, cex = .8)
      for (j in which(!z$Available)) text(z$StatusX[j], y[j], z$StatusLabel[j],
        adj = z$StatusHjust[j], col = z$Colour[j], cex = .7)
    } else {
      mtext(labels[plan], side = 3, line = .5, cex = .85)
      axis(2, at = z$Row, labels = z$MethodLabel, las = 1, cex.axis = .8)
    }
    segments(z$Lower, y, z$Upper, y, col = z$Colour, lwd = 2, lty = z$LineType)
    points(z$Estimate, y, col = z$Colour, pch = z$Shape)
  }
  invisible(out)
}

.mfrm_gg_dstudy_intervals <- function(x) {
  .require_mfrmr_ggplot2()
  payload <- x$data; tab <- payload$table
  if (identical(payload$layout, "combined")) {
    methods <- payload$methods
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$Estimate, y = .data$Y,
      colour = .data$Method)) +
      ggplot2::geom_blank(ggplot2::aes(x = 0)) +
      ggplot2::geom_hline(yintercept = payload$plans$Y[-1L]+.5, colour = "grey90") +
      ggplot2::geom_segment(data = tab[tab$Available, , drop = FALSE],
        ggplot2::aes(x = .data$Lower, xend = .data$Upper, yend = .data$Y,
          linetype = .data$Method), linewidth = .7) +
      ggplot2::geom_point(data = tab[is.finite(tab$Estimate), , drop = FALSE],
        ggplot2::aes(shape = .data$Shape), size = 2.5) +
      ggplot2::geom_text(data = tab[!tab$Available, , drop = FALSE],
        ggplot2::aes(x = .data$StatusX, label = .data$StatusLabel,
          hjust = .data$StatusHjust), size = 3, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = setNames(payload$palette, methods),
        limits = methods, drop = FALSE) +
      ggplot2::scale_linetype_manual(values = setNames(seq_along(methods), methods)) +
      ggplot2::scale_shape_identity() +
      ggplot2::scale_x_continuous(limits = payload$xrange) +
      ggplot2::scale_y_continuous(breaks = payload$plans$Y, labels = payload$plans$Label,
        limits = c(.5, nrow(payload$plans)+.5)) +
      ggplot2::guides(linetype = "none", colour = ggplot2::guide_legend(title = NULL,
        nrow = if (length(methods) > 2L) length(methods) else 1L,
        override.aes = list(shape = 15L+seq_along(methods), linetype = seq_along(methods)))) +
      .mfrmr_gg_theme() + ggplot2::theme(legend.position = "top",
        panel.grid.major.y = ggplot2::element_blank())
    return(.mfrmr_gg_labs(p, payload, x = payload$xlab, y = NULL))
  }
  # Facet-specific row labels keep missing intervals visible without assigning
  # an artificial coefficient to a missing estimate.
  order <- unlist(lapply(unique(tab$Method), function(m) unique(tab$MethodLabel[tab$Method == m])))
  tab$Label <- factor(tab$MethodLabel, levels = rev(order))
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$Estimate, y = .data$Label)) +
    ggplot2::geom_blank(ggplot2::aes(x = 0)) +
    ggplot2::geom_segment(data = tab[tab$Available, , drop = FALSE],
      ggplot2::aes(x = .data$Lower, xend = .data$Upper, yend = .data$Label,
        colour = .data$Colour), linewidth = .7) +
    ggplot2::geom_point(data = tab[is.finite(tab$Estimate), , drop = FALSE],
      ggplot2::aes(colour = .data$Colour, shape = .data$Shape), size = 2.5) +
    ggplot2::scale_colour_identity() + ggplot2::scale_shape_identity() +
    ggplot2::scale_x_continuous(limits = payload$xrange) +
    ggplot2::facet_wrap(~Plan, ncol = payload$ncol, scales = "free_y") +
    .mfrmr_gg_theme()
  .mfrmr_gg_labs(p, payload, x = payload$xlab, y = NULL)
}
