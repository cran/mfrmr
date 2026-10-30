#' @rdname estimate_all_bias
#' @details `plot()` method: Explore bias across facet pairs
#' @param x A result from [estimate_all_bias()].
#' @param pair An exact `Interaction` label from `summary(x)$overview`, or its
#'   row number. `NULL` (default) shows all requested pairs.
#' @param type `"overview"` shows signed cell estimates grouped by pair.
#'   With `pair` supplied, the default becomes `"heatmap"`; `"scatter"`,
#'   `"ranked"`, `"abs_t_hist"` and `"facet_profile"` select the corresponding
#'   [plot.mfrm_bias()] view.
#' @param draw Draw the figure. `FALSE` returns its data without opening a device.
#' @param preset A package plot preset, for example `"standard"` or `"monochrome"`.
#' @param ... For `plot()`: With `pair`, additional arguments to [plot.mfrm_bias()], such as
#'   `top_n` or `sort_by`. The overview rejects additional arguments.
#'
#' For `summary()`: Reserved; additional arguments are rejected.
#' @details For `plot()`: The overview plots every finite retained bias estimate, without ranking or
#' a significance cutoff. Positive shifts increase expected coded scores under
#' the fitted model; score meaning depends on the rubric. Vertical spreading
#' only separates cells visually and has no statistical meaning.
#'
#' Filled circles have a usable conditional screening statistic; open circles
#' retain a point without one. These symbols indicate availability, not whether
#' a cell is flagged. Pair labels report usable statistics / retained cells and
#' the range of observed row counts per cell (`n`), not independent Person
#' counts. Pairs share observations. Empty, failed and unavailable pairs remain
#' labeled, without substituting a point at zero. The vertical zero line means
#' no conditional shift, not a threshold for practical importance.
#'
#' No model or diagnostics are recalculated. The overview supplies no intervals,
#' multiplicity correction or formal bias test. Review the individual cell
#' counts, optimizer status and subject-matter explanations; selecting a pair
#' after this overview does not make its subsequent analysis confirmatory.
#' Titles and interpretation notes are returned outside the figure.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object. For the overview,
#'   [plot_data()] provides `table` (all combined cells, including unavailable
#'   rows), `overview`, `errors`, `optimization_failures`, `fit_identity`,
#'   `settings`, axis labels and notes. [as_ggplot()] reproduces the overview.
#'   A selected pair returns the ordinary bias plot with `collection_pair`
#'   and `fit_identity` provenance added.
#' @seealso [estimate_all_bias()], [summary.mfrm_bias_collection()],
#'   [plot.mfrm_bias()]
#' @examples # Long-running workflow; run explicitly when needed.
#' # example("plot.mfrm_bias_collection", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML")
#' screens <- estimate_all_bias(fit)
#' p <- plot(screens)
#' plot_data(p, "overview")
#' plot(screens, pair = "Rater x Criterion", type = "ranked")
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_bias_collection <- function(x, pair = NULL,
    type = c("overview", "scatter", "ranked", "heatmap", "abs_t_hist", "facet_profile"),
    draw = TRUE, preset = "standard", ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw))
    stop("`draw` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(pair) && missing(type)) type <- "heatmap"
  type <- match.arg(type)
  if (!is.null(pair)) {
    labels <- x$summary$Interaction
    if (is.numeric(pair) && length(pair) == 1L && is.finite(pair) &&
        pair == floor(pair) && pair >= 1 && pair <= length(labels)) pair <- labels[pair]
    if (!is.character(pair) || length(pair) != 1L || is.na(pair) || !pair %in% labels)
      stop("`pair` must be an exact Interaction label or its row number in summary(x)$overview.", call. = FALSE)
    if (type == "overview") stop("Use `pair = NULL` for the overview.", call. = FALSE)
    if (is.null(x$by_pair[[pair]]$table) || !nrow(x$by_pair[[pair]]$table)) {
      reason <- x$summary$Error[match(pair, labels)]
      stop("No retained cells for ", pair, if (!is.na(reason) && nzchar(reason)) paste0(": ", reason),
        ". Inspect summary(x)$overview.", call. = FALSE)
    }
    if (!is.null(x$fit_identity)) mfrm_validate_result_identity(x$by_pair[[pair]],
      x$fit_identity, paste0("by_pair[[", pair, "]]"), "estimate_bias()", "plot.mfrm_bias_collection()")
    out <- plot(x$by_pair[[pair]], type = type, draw = draw, preset = preset, ...)
    out$data$collection_pair <- pair
    out$data$fit_identity <- x$fit_identity
    return(invisible(out))
  }
  if (type != "overview") stop("Select `pair` for an individual bias view.", call. = FALSE)
  rlang::check_dots_empty()
  style <- resolve_plot_preset(preset)
  review <- summary(x)
  overview <- review$overview
  table <- review$cells
  if (!"Bias Size" %in% names(table)) table[["Bias Size"]] <- numeric(nrow(table))
  table$PlotY <- rep(NA_real_, nrow(table))
  labels <- character(nrow(overview))
  for (i in seq_len(nrow(overview))) {
    at <- which(table$Interaction == overview$Interaction[i])
    table$PlotY[at] <- nrow(overview) + 1L - i +
      if (length(at) > 1L) seq(-.22, .22, length.out = length(at)) else 0
    counts <- if ("ObsN" %in% names(table)) table$ObsN[at] else
      if ("Observd Count" %in% names(table)) table[["Observd Count"]][at] else numeric()
    counts <- counts[is.finite(counts)]
    support <- if (length(counts)) paste0("n = ", paste(unique(range(counts)), collapse = "-")) else "n unavailable"
    labels[i] <- paste0(overview$Interaction[i], "\n",
      if (overview$Status[i] %in% c("available", "partial", "unavailable"))
        paste0(overview$ScreenAvailable[i], "/", overview$Rows[i], " screens; ", support) else
        overview$Status[i])
  }
  shown <- table$PointAvailable
  limits <- range(c(0, table[["Bias Size"]][shown]))
  if (diff(limits) == 0) limits <- limits + c(-.5, .5)
  out <- new_mfrm_plot_data("bias_collection", list(table = table,
    overview = overview, errors = review$errors,
    optimization_failures = review$optimization_failures,
    fit_identity = review$fit_identity, settings = review$settings,
    pair_labels = labels, xlim = limits, preset = style$name,
    xlab = "Conditional bias shift (native ability units)", ylab = "",
    title = "Exploratory bias across facet pairs",
    notes = data.frame(Type = "interpretation", Text = c(review$notes,
      "Filled circles: screening statistic available; open circles: point only. Vertical spreading is for visibility.",
      "n is the range of observed row counts per cell, not independent Person counts. No intervals or multiplicity correction are drawn."))))
  if (!draw) return(invisible(out))
  apply_plot_preset(style)
  old <- graphics::par(c("mar", "mex"))
  on.exit(graphics::par(old), add = TRUE)
  label_width <- max(c(0, graphics::strwidth(labels, units = "inches", cex = .8)))
  if (label_width + .2 + 1.25 >= graphics::par("fin")[1])
    stop("Pair labels need a wider graphics device; use `draw = FALSE` to inspect the table.", call. = FALSE)
  graphics::par(mex = .8)
  left <- label_width / (graphics::par("csi") * graphics::par("mex")) + 1.2
  graphics::par(mar = c(4.5, max(4, left), 1, 1))
  plot(NA_real_, NA_real_, xlim = limits,
    ylim = c(.5, max(1, nrow(overview)) + .5), yaxt = "n",
    xlab = out$data$xlab, ylab = "")
  abline(v = 0, col = style$neutral, lty = 2)
  axis(2, at = rev(seq_len(nrow(overview))), labels = labels,
    las = 1, tick = FALSE, cex.axis = .8)
  points(table[["Bias Size"]][shown], table$PlotY[shown],
    pch = ifelse(table$ScreenAvailable[shown], 16, 1), cex = .9, col = style$axis)
  invisible(out)
}

.mfrmr_gg_bias_collection <- function(payload) {
  style <- resolve_plot_preset(payload$preset)
  tab <- payload$table
  tab <- tab[tab$PointAvailable, , drop = FALSE]
  ggplot2::ggplot(tab, ggplot2::aes(x = .data[["Bias Size"]], y = .data$PlotY)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = style$neutral) +
    ggplot2::geom_point(shape = ifelse(tab$ScreenAvailable, 16, 1), size = 2, colour = style$axis) +
    ggplot2::scale_x_continuous(limits = payload$xlim) +
    ggplot2::scale_y_continuous(breaks = rev(seq_len(nrow(payload$overview))),
      labels = payload$pair_labels, limits = c(.5, max(1, nrow(payload$overview)) + .5)) +
    ggplot2::labs(x = payload$xlab, y = NULL) + .mfrmr_gg_theme()
}
