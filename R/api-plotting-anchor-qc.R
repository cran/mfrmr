# ============================================================================
# Anchor Drift & Equating Chain Plots
# ============================================================================

# --- Internal plot helpers (not exported) ------------------------------------

.plot_drift_dot <- function(dt, config, draw = TRUE, style = resolve_plot_preset("standard"), ...) {
  out <- new_mfrm_plot_data(
    "anchor_drift",
    list(
      plot = "drift",
      table = dt,
      title = "Anchor drift",
      subtitle = paste0("Review threshold: |drift| >= ", format(config$drift_threshold)),
      legend = new_plot_legend(
        label = c("Within review band", "Flagged drift"),
        role = c("status", "status"),
        aesthetic = c("point", "point"),
        value = c(style$accent_primary, style$warn)
      ),
      reference_lines = new_reference_lines(
        axis = c("v", "v", "v"),
        value = c(-config$drift_threshold, 0, config$drift_threshold),
        label = c("Drift review threshold", "Centered drift reference", "Drift review threshold"),
        linetype = c("dotted", "dashed", "dotted"),
        role = c("threshold", "reference", "threshold")
      ),
      preset = style$name
    )
  )
  if (!draw) return(invisible(out))

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)

  dt <- dt |> dplyr::arrange(dplyr::desc(abs(.data$Drift)))
  labels <- paste0(dt$Facet, ":", dt$Level)
  n <- nrow(dt)

  max_abs <- max(abs(dt$Drift), na.rm = TRUE) * 1.2

  graphics::par(mar = c(4, 8, 3, 1))
  graphics::plot(dt$Drift, seq_len(n), xlim = c(-max_abs, max_abs),
                 yaxt = "n", xlab = "Drift (logits)", ylab = "",
                 main = "Anchor Drift", pch = 19,
                 col = ifelse(dt$Flag, style$warn, style$accent_primary), ...)
  graphics::axis(2, at = seq_len(n), labels = labels, las = 1, cex.axis = 0.7)
  graphics::abline(v = pretty(c(-max_abs, max_abs), n = 5), col = grDevices::adjustcolor(style$grid, alpha.f = 0.85), lty = 1)
  graphics::abline(v = 0, lty = 2, col = grDevices::adjustcolor(style$foreground, alpha.f = 0.7))
  graphics::abline(v = c(-config$drift_threshold, config$drift_threshold),
                   lty = 3, col = grDevices::adjustcolor(style$warn, alpha.f = 0.9))

  invisible(out)
}

.plot_drift_heatmap <- function(dt, config, draw = TRUE, style = resolve_plot_preset("standard"), ...) {
  dt_wide <- dt |>
    dplyr::mutate(Element = paste0(.data$Facet, ":", .data$Level)) |>
    dplyr::select("Element", "Wave", "Drift")

  mat <- tryCatch({
    tidyr::pivot_wider(dt_wide, names_from = "Wave",
                        values_from = "Drift") |>
      tibble::column_to_rownames("Element") |>
      as.matrix()
  }, error = function(e) NULL)

  if (is.null(mat) || nrow(mat) == 0) {
    if (draw) message("Insufficient data for heatmap.")
    return(invisible(NULL))
  }

  out <- new_mfrm_plot_data(
    "anchor_drift",
    list(
      plot = "heatmap",
      table = dt,
      matrix = mat,
      title = "Anchor drift heatmap",
      subtitle = paste0("Wave-by-element drift; review threshold = ", format(config$drift_threshold)),
      legend = new_plot_legend(
        label = c("Negative drift", "Positive drift"),
        role = c("drift", "drift"),
        aesthetic = c("heatmap", "heatmap"),
        value = c(style$accent_secondary, style$warn)
      ),
      reference_lines = new_reference_lines(),
      preset = style$name
    )
  )

  if (!draw) return(invisible(out))

  max_abs <- max(abs(mat), na.rm = TRUE)
  n_colors <- 21
  breaks <- seq(-max_abs, max_abs, length.out = n_colors + 1)
  blues <- grDevices::colorRampPalette(c(style$accent_secondary, "white", style$warn))(n_colors)

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mar = c(5, 8, 3, 2))

  graphics::image(t(mat[nrow(mat):1, , drop = FALSE]),
                  axes = FALSE, col = blues, breaks = breaks,
                  main = "Anchor Drift Heatmap", ...)
  graphics::axis(1, at = seq(0, 1, length.out = ncol(mat)),
                 labels = colnames(mat), las = 2, cex.axis = 0.8)
  graphics::axis(2, at = seq(0, 1, length.out = nrow(mat)),
                 labels = rev(rownames(mat)), las = 1, cex.axis = 0.7)

  invisible(out)
}

.plot_equating_chain <- function(x, draw = TRUE, style = resolve_plot_preset("standard"), ...) {
  cum <- x$cumulative
  out <- new_mfrm_plot_data(
    "anchor_drift",
    list(
      plot = "chain",
      table = cum,
      links = x$links,
      title = "Equating chain",
      subtitle = "Cumulative offsets across linked calibration waves",
      legend = new_plot_legend(
        label = c("Cumulative offset", "Centered chain reference"),
        role = c("offset", "reference"),
        aesthetic = c("line-point", "line"),
        value = c(style$accent_primary, style$foreground)
      ),
      reference_lines = new_reference_lines("h", 0, "Centered chain reference", "dashed", "reference"),
      preset = style$name
    )
  )
  if (!draw) return(invisible(out))

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mar = c(5, 4, 3, 1))

  n <- nrow(cum)
  graphics::plot(seq_len(n), cum$Cumulative_Offset, type = "b",
                 pch = 19, col = style$accent_primary, lwd = 2,
                 xaxt = "n", xlab = "", ylab = "Cumulative Offset (logits)",
                 main = "Equating Chain", ...)
  graphics::axis(1, at = seq_len(n), labels = cum$Wave, las = 2, cex.axis = 0.8)
  graphics::abline(h = pretty(cum$Cumulative_Offset, n = 5), col = grDevices::adjustcolor(style$grid, alpha.f = 0.85), lty = 1)
  graphics::abline(h = 0, lty = 2, col = grDevices::adjustcolor(style$foreground, alpha.f = 0.7))

  links <- x$links
  for (i in seq_len(nrow(links))) {
    mid_x <- i + 0.5
    mid_y <- (cum$Cumulative_Offset[i] + cum$Cumulative_Offset[i + 1]) / 2
    graphics::text(mid_x, mid_y, sprintf("n=%d", links$N_Common[i]),
                   cex = 0.7, col = grDevices::adjustcolor(style$foreground, alpha.f = 0.82))
  }

  invisible(out)
}

# --- Exported plot function --------------------------------------------------

#' Plot anchor drift or a screened linking chain
#'
#' Creates base-R plots for inspecting anchor drift across calibration waves
#' or visualising the cumulative offset in a screened linking chain.
#'
#' @param x An `mfrm_anchor_drift` or `mfrm_equating_chain` object.
#' @param type Plot type: `"drift"` (dot plot of element drift),
#'   `"chain"` (cumulative offset line plot), `"heatmap"`
#'   (wave-by-element drift heatmap), or `"forest"` (per-(Facet, Level,
#'   Wave) anchor estimate with `+/- z * SE` whiskers; requires
#'   `mfrm_anchor_drift`).
#' @param facet Optional character vector to filter drift plots to specific
#'   facets.
#' @param ci_level Confidence level used by `type = "forest"` for the
#'   anchor-estimate whiskers (default `0.95`). Ignored for other
#'   plot types.
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`, or `"monochrome"`).
#' @param draw If `FALSE`, return the plot data invisibly without drawing.
#' @param ... Additional graphical parameters passed to base plotting
#'   functions.
#'
#' @details
#' Three plot types are supported:
#'
#' - **`"drift"`** (for `mfrm_anchor_drift` objects): A dot plot of each
#'   element's drift value, grouped by facet.  Horizontal reference lines
#'   mark the drift threshold.  Red points indicate flagged elements.
#' - **`"heatmap"`** (for `mfrm_anchor_drift` objects): A wave-by-element
#'   heat matrix showing drift magnitude.  Darker cells represent larger
#'   absolute drift.  Useful for spotting systematic patterns (e.g., all
#'   criteria shifting in the same direction).
#' - **`"chain"`** (for `mfrm_equating_chain` objects): A line plot of
#'   cumulative offsets across the screened linking chain. A flatter line
#'   indicates smaller between-wave shifts; steep segments suggest larger
#'   link offsets that deserve review.
#'
#' @section Which plot should I use?:
#' - Use `type = "drift"` with an `mfrm_anchor_drift` object to review flagged
#'   elements directly.
#' - Use `type = "heatmap"` with an `mfrm_anchor_drift` object to spot
#'   wave-by-element patterns.
#' - Use `type = "chain"` with an `mfrm_equating_chain` object after
#'   [build_equating_chain()] to inspect cumulative offsets across waves.
#'
#' @section Interpreting plots:
#' **Drift** is the change in an element's estimated measure between
#' calibration waves, after accounting for the screened common-element link
#' offset. An
#' element is flagged when its absolute drift exceeds a threshold
#' (typically 0.5 logits) **and** the drift-to-SE ratio exceeds a
#' secondary criterion (typically 2.0), ensuring that only
#' practically noticeable and relatively precise shifts are flagged.
#'
#' - In drift and heatmap plots, red or dark-shaded elements exceed
#'   both thresholds.  Common causes include rater drift over time,
#'   item exposure effects, or curriculum changes.
#' - In chain plots, uneven spacing between waves suggests differential
#'   shifts in the screened linking offsets. The \eqn{y}-axis shows cumulative
#'   logit-scale offsets; flatter segments indicate more stable adjacent links.
#'   Steep segments should be checked alongside `LinkSupportAdequate` and the
#'   retained common-element counts before making longitudinal claims.
#' - For drift objects, it is usually best to read `summary(x)` first
#'   and then use the plot to see where the flagged values sit.
#'
#' @section Typical workflow:
#' 1. Build a drift or screened-linking object with [detect_anchor_drift()] or
#'    [build_equating_chain()].
#' 2. Start with `draw = FALSE` if you want the plotting data for custom
#'    reporting.
#' 3. Use the base-R plot for quick screening and then inspect the underlying
#'    tables for exact values.
#'
#' @section Further guidance:
#' For a plot-selection guide and a longer walkthrough, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @return A plotting-data object of class `mfrm_plot_data`. With
#'   `draw = FALSE`, `result$data$table` contains the filtered drift or chain
#'   table, `result$data$matrix` contains the heatmap matrix when requested,
#'   and the returned plot data includes package-native `title`, `subtitle`,
#'   `legend`, and `reference_lines`.
#'
#' @seealso [detect_anchor_drift()], [build_equating_chain()],
#'   [plot_dif_heatmap()], [plot_bubble()], [mfrmr_visual_diagnostics]
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept linking
#' @export
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' people <- unique(toy$Person)
#' d1 <- toy[toy$Person %in% people[1:12], , drop = FALSE]
#' d2 <- toy[toy$Person %in% people[13:24], , drop = FALSE]
#' fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 30)
#' fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 30)
#' drift <- detect_anchor_drift(list(W1 = fit1, W2 = fit2))
#' drift_plot <- plot_anchor_drift(drift, type = "drift", draw = FALSE)
#' class(drift_plot)
#' names(drift_plot$data)
#' chain <- build_equating_chain(list(F1 = fit1, F2 = fit2))
#' chain_plot <- plot_anchor_drift(chain, type = "chain", draw = FALSE)
#' head(chain_plot$data$table)
#' if (interactive()) {
#'   plot_anchor_drift(drift, type = "heatmap", preset = "publication")
#' }
#' }
plot_anchor_drift <- function(x, type = c("drift", "chain", "heatmap", "forest"),
                              facet = NULL,
                              ci_level = 0.95,
                              preset = c("standard", "publication", "compact", "monochrome"),
                              draw = TRUE, ...) {
  type <- match.arg(type)
  style <- resolve_plot_preset(preset)

  if (inherits(x, "mfrm_equating_chain")) {
    if (type == "chain") {
      if (isTRUE(draw)) apply_plot_preset(style)
      return(.plot_equating_chain(x, draw = draw, style = style, ...))
    }
  }

  if (inherits(x, "mfrm_anchor_drift")) {
    dt <- x$drift_table
    if (!is.null(facet)) dt <- dt |> dplyr::filter(.data$Facet %in% facet)

    if (nrow(dt) == 0) {
      if (draw) message("No drift data to plot.")
      return(invisible(NULL))
    }

    if (type == "drift") {
      if (isTRUE(draw)) apply_plot_preset(style)
      return(.plot_drift_dot(dt, x$config, draw = draw, style = style, ...))
    } else if (type == "heatmap") {
      if (isTRUE(draw)) apply_plot_preset(style)
      return(.plot_drift_heatmap(dt, x$config, draw = draw, style = style, ...))
    } else if (type == "forest") {
      if (isTRUE(draw)) apply_plot_preset(style)
      return(.plot_drift_forest(dt, x$config, ci_level = ci_level,
                                 draw = draw, style = style, ...))
    }
  }

  stop("Unsupported object class or plot type combination.", call. = FALSE)
}

# Forest-style anchor-drift visualisation. Each (Facet, Level, Wave)
# row is rendered as a horizontal CI whisker around the per-wave
# anchor estimate; rows are grouped by Facet level so reviewers see
# one ladder of CI bands per anchored element across waves.
.plot_drift_forest <- function(drift_tbl, config, ci_level = 0.95,
                               draw = TRUE, style = NULL, ...) {
  if (is.null(style)) style <- resolve_plot_preset("standard")
  dt <- as.data.frame(drift_tbl, stringsAsFactors = FALSE)
  est_col <- if ("Estimate" %in% names(dt)) "Estimate" else if ("Wave_Est" %in% names(dt)) "Wave_Est" else NA_character_
  if (is.na(est_col) || !all(c("Facet", "Level", "Wave") %in% names(dt))) {
    stop("Forest plot requires Facet, Level, Wave, and Estimate (or Wave_Est) columns.",
         call. = FALSE)
  }
  dt$Estimate <- dt[[est_col]]
  se_col <- if ("SE_Wave" %in% names(dt)) "SE_Wave" else if ("SE" %in% names(dt)) "SE" else if ("ModelSE" %in% names(dt)) "ModelSE" else NA_character_
  z_ci <- stats::qnorm(1 - (1 - ci_level) / 2)
  if (!is.na(se_col)) {
    dt$CI_Lower <- dt$Estimate - z_ci * dt[[se_col]]
    dt$CI_Upper <- dt$Estimate + z_ci * dt[[se_col]]
  } else {
    dt$CI_Lower <- NA_real_
    dt$CI_Upper <- NA_real_
  }
  dt$Label <- paste0(dt$Facet, ":", dt$Level, " @ ", dt$Wave)
  dt <- dt[order(dt$Facet, dt$Level, dt$Wave), , drop = FALSE]
  out <- new_mfrm_plot_data(
    "anchor_drift_forest",
    list(
      data = dt,
      ci_level = ci_level,
      title = "Anchor drift forest",
      subtitle = sprintf("Per-wave anchor estimates with %g%% CI",
                         round(100 * ci_level)),
      preset = style$name,
      legend = new_plot_legend(
        label = c("Anchor estimate", "CI whisker"),
        role = c("location", "uncertainty"),
        aesthetic = c("point", "segment"),
        value = c(style$accent_primary, style$accent_primary)
      ),
      reference_lines = new_reference_lines("v", 0,
                                             "Centred logit reference",
                                             "dashed", "reference")
    )
  )
  if (isTRUE(draw)) {
    y <- seq_len(nrow(dt))
    xrange <- range(c(dt$Estimate, dt$CI_Lower, dt$CI_Upper, 0),
                    finite = TRUE, na.rm = TRUE)
    graphics::plot(
      x = dt$Estimate, y = y, type = "n",
      xlim = xrange, yaxt = "n",
      xlab = "Estimate (logit)", ylab = "",
      main = "Anchor drift forest"
    )
    graphics::title(sub = sprintf("Per-wave anchor estimates with %g%% CI",
                                   round(100 * ci_level)),
                    line = 2.2, cex.sub = 0.9)
    graphics::abline(v = 0, lty = 2, col = style$neutral)
    valid <- is.finite(dt$CI_Lower) & is.finite(dt$CI_Upper)
    if (any(valid)) {
      graphics::segments(dt$CI_Lower[valid], y[valid],
                         dt$CI_Upper[valid], y[valid],
                         col = style$accent_primary, lwd = 2)
    }
    graphics::points(dt$Estimate, y, pch = 19, col = style$accent_primary)
    graphics::axis(2, at = y, labels = dt$Label, las = 1, cex.axis = 0.7)
  }
  invisible(out)
}

# ============================================================================
# QC Pipeline Plot
# ============================================================================

#' Plot QC pipeline results
#'
#' Visualizes the output from [run_qc_pipeline()] as either a traffic-light
#' bar chart or a detail panel showing values versus thresholds.
#'
#' @param x Output from [run_qc_pipeline()].
#' @param type Plot type: `"traffic_light"` (default) or `"detail"`.
#' @param draw If `FALSE`, return plot data invisibly without drawing.
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @details
#' Two plot types are provided for visual triage of QC results:
#'
#' - **`"traffic_light"`** (default): A horizontal bar chart with one row
#'   per QC check.  Bars are coloured green (Pass), amber (Warn), or red
#'   (Fail).  Provides an at-a-glance summary of the current QC review state.
#' - **`"detail"`**: A panel showing each check's observed value and its
#'   pass/warn/fail thresholds.  Useful for understanding how close a
#'   borderline result is to the next verdict level.
#'
#' @section QC checks performed:
#' The pipeline evaluates up to 10 checks (depending on available
#' diagnostics):
#' 1. **Convergence**: did the optimizer converge?
#' 2. **Overall Infit**: global information-weighted mean-square
#' 3. **Overall Outfit**: global unweighted mean-square
#' 4. **Misfit rate**: proportion of elements with \eqn{|\mathrm{ZSTD}| > 2}
#' 5. **Category usage**: minimum observations per score category
#' 6. **Disordered steps**: whether threshold estimates are monotonic
#' 7. **Separation** (per facet): element discrimination adequacy
#' 8. **Residual PCA eigenvalue**: first-component eigenvalue (if computed)
#' 9. **Displacement**: maximum absolute displacement across elements
#' 10. **Inter-rater agreement**: minimum pairwise exact agreement
#'
#' @section Interpreting plots:
#' - **Green** (Pass): the check meets the current threshold-profile criteria.
#' - **Amber** (Warn): borderline---monitor but not necessarily
#'   disqualifying.  Review the detail panel to see how close the value
#'   is to the fail threshold.
#' - **Red** (Fail): requires investigation before strong operational or
#'   interpretive claims are made from the current run. Common remedies include collapsing categories
#'   (for disordered steps), removing outlier raters (for misfit), or
#'   increasing sample size (for low separation).
#' - The detail view shows numeric values, making it easy to communicate
#'   exact results to stakeholders.
#'
#' @return Invisible verdicts tibble from the QC pipeline.
#'
#' @seealso [run_qc_pipeline()], [plot_qc_dashboard()],
#'   [build_visual_summaries()], [mfrmr_visual_diagnostics]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("study1")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' qc <- run_qc_pipeline(fit)
#' plot_qc_pipeline(qc, draw = FALSE)
#' }
#' @export
plot_qc_pipeline <- function(x, type = c("traffic_light", "detail"),
                             draw = TRUE, ...) {
  type <- match.arg(type)
  stopifnot(inherits(x, "mfrm_qc_pipeline"))

  vt <- x$verdicts
  if (!draw) return(invisible(vt))

  n <- nrow(vt)
  cols <- ifelse(vt$Verdict == "Pass", "#2ca02c",
                 ifelse(vt$Verdict == "Warn", "#ff7f0e",
                        ifelse(vt$Verdict == "Fail", "#d62728", "#999999")))

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)

  if (type == "traffic_light") {
    graphics::par(mar = c(3, 14, 3, 4))

    graphics::plot(NULL, xlim = c(0, 1), ylim = c(0.5, n + 0.5),
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                   main = paste("QC Pipeline:", x$overall), ...)

    for (i in seq_len(n)) {
      graphics::rect(0, i - 0.4, 1, i + 0.4, col = cols[i], border = NA)
      graphics::text(0.5, i, vt$Verdict[i], col = "white", font = 2, cex = 0.9)
    }

    graphics::axis(2, at = seq_len(n), labels = vt$Check,
                   las = 1, cex.axis = 0.8, tick = FALSE)

    for (i in seq_len(n)) {
      graphics::mtext(vt$Value[i], side = 4, at = i,
                      las = 1, cex = 0.6, line = 0.5)
    }
  } else {
    graphics::par(mar = c(3, 14, 3, 8))

    graphics::plot(NULL, xlim = c(0, 1), ylim = c(0.5, n + 0.5),
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                   main = paste("QC Pipeline Detail:", x$overall), ...)

    for (i in seq_len(n)) {
      graphics::rect(0, i - 0.4, 0.15, i + 0.4, col = cols[i], border = NA)
      graphics::text(0.075, i, substr(vt$Verdict[i], 1, 1),
                     col = "white", font = 2, cex = 0.8)
      graphics::text(0.2, i, vt$Detail[i], adj = 0, cex = 0.7)
    }

    graphics::axis(2, at = seq_len(n), labels = vt$Check,
                   las = 1, cex.axis = 0.8, tick = FALSE)
  }

  invisible(vt)
}
