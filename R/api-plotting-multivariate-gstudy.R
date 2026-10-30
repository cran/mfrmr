#' Inspect sources of observed-score variation
#'
#' Plot raw variance components for selected scores or one prespecified weighted
#' score, or inspect the component covariance matrices. Then use [mfrm_d_study()] to assess
#' how changing the measurement conditions affects G, Phi and SEM.
#'
#' @param x A result from [mfrm_g_study()].
#' @param score Original score names in display order. Multiple scores produce
#'   separate component panels. The default selects the first score for
#'   `type = "components"`, or all scores for `type = "covariance"`.
#'   Do not supply both `score` and `weights`.
#' @param weights Optional named numeric vector defining one composite. Name
#'   every score exactly once, with finite weights and at least one nonzero
#'   value. Names determine order. Signed weights permit difference scores;
#'   weights are not normalized. For a saved D-study, use its selected weight
#'   vector, for example `plot(d$gstudy, weights = d$weights)`.
#' @param draw Draw the figure; `FALSE` only returns its data.
#' @param preset Plot style: `"standard"`, `"publication"`, `"compact"`, or
#'   `"monochrome"`. Omission uses the session plot default.
#' @param type `"components"` for variance-component points, or `"covariance"`
#'   for one matrix heatmap per source. With `intervals = TRUE`, covariance
#'   instead displays one signed interval plot for exactly two selected scores.
#'   Covariance plots do not accept `weights`.
#' @param ncol Positive number of panel columns. Defaults to at most two.
#' @param intervals Logical; `FALSE` (default) shows point estimates. `TRUE`
#'   adds component confidence intervals in the same call; also specify
#'   `assumption = "normal"`. Available for complete balanced ANOVA component
#'   plots, including supported nesting. Covariance intervals require an
#'   explicit method and exactly two scores. Explicit `method = "reml_profile"`
#'   also provides variance and two-score covariance intervals for identifiable incomplete sources,
#'   fitting a separate REML model while preserving the original points.
#'   `method = "residual_contrast"` uses retained ratings without likelihood
#'   fitting, but only Residual receives an interval. Other points remain
#'   visible with unavailable intervals and their reasons.
#'   Other methods use saved ANOVA summaries without refitting.
#' @param assumption With `intervals = TRUE`, explicitly specify `"normal"`
#'   for independent normal random effects with constant component variances.
#'   This assumption is not tested.
#' @param level Pointwise confidence level for `intervals = TRUE`, default 0.95.
#' @param method Component interval method passed to
#'   [confint.mfrm_multivariate_gstudy()]. Single-score/composite variances
#'   use `"auto"` or `"mls"` for balanced ANOVA, or explicit experimental
#'   `"reml_profile"` for balanced/identifiable incomplete sources. Profile
#'   intervals have approximate calibration and no accompanying sampling SE.
#'   Covariance intervals require explicit `"gci"`, `"regression_region"` or
#'   `"reml_profile"`; none is an automatic default.
#'   Explicit `"residual_contrast"` also supports incomplete sources for
#'   Residual only: exact normal-model variance limits or conservative
#'   covariance limits. It retains original points and a separate interval
#'   center without supplying the original estimator's sampling SE.
#' @param draws,seed Numerical settings for covariance `method = "gci"`.
#'   Defaults are 20000 pivotal draws and seed 20261016. See the confidence
#'   interval help for requirements and numerical endpoint bands.
#' @param ... Reserved; additional arguments are rejected.
#'
#' @details Each point is a variance-component estimate in squared score units.
#'   For weights \eqn{w}, the component is \eqn{w' \Sigma_s w}, including
#'   between-score covariances. Multiplying all weights by a constant multiplies
#'   each component by its square. Components describe the source model at the
#'   individual rating level; they have not been divided by future facet counts.
#'   They are not percentages, measurement-error shares for a future mean score,
#'   or MFRM severity estimates. A larger Person component describes greater
#'   person variation, not a defect in the assessment. A large rater component
#'   does not by itself identify which rater should receive feedback.
#'
#'   Negative estimates remain to the left of zero as open circles; nonnegative
#'   estimates use filled circles. An unconstrained negative estimate is not a
#'   negative population variance. Inspect the data/design and component
#'   diagnostics before interpreting it. No component is clipped, normalized
#'   or removed. A covariance matrix can be non-positive-semidefinite even when
#'   its projection for the selected score is positive; this is retained in the
#'   returned table and notes.
#'
#'   With `intervals = TRUE`, horizontal bars replace the stems from zero.
#'   They show the nonnegative-space intervals from
#'   [confint.mfrm_multivariate_gstudy()]: approximate Ting MLS intervals for
#'   components and chi-square intervals for the combined residual. They are
#'   pointwise, not simultaneous across components or scores. Negative point
#'   estimates can lie outside their nonnegative-space intervals. Square
#'   outlines mark estimates without an interval; consult `IntervalStatus`
#'   and `unavailable_intervals` for the reason. An empty interval or degenerate
#'   mean squares do not mean zero uncertainty. Raw endpoints, SEs and the full
#'   interval tables remain in the returned data. These intervals concern
#'   source variances, not future mean-score errors, G/Phi or off-diagonal
#'   covariance elements. With `intervals = FALSE`, stems run from zero to
#'   the estimate and do not represent uncertainty.
#'
#'   With `method = "reml_profile"`, the bars instead show approximate
#'   restricted-likelihood component intervals. The points remain the original
#'   ANOVA/MINQUE estimates; the separate `IntervalCenter` is retained in the
#'   table. Profile `SE` and inapplicable MLS `RawLower/RawUpper` columns are
#'   `NA`; `SEStatus` explains this. Numerical availability does not establish
#'   coverage, especially near zero variances. The interval surrounds the REML
#'   center but does not enumerate remote disconnected accepted regions.
#'
#'   With `method = "residual_contrast"`, only Residual receives a bar, based
#'   on a projection that removes all nonresidual effects. Points still show
#'   original ANOVA/MINQUE components; `IntervalCenter` records the projected
#'   residual estimate. No remaining residual degrees of freedom or numerically
#'   unresolved residual variation means unavailable limits. Exact variance
#'   coverage and conservative covariance coverage require the normal residual
#'   model at a fixed observed assignment; they do not extend to other sources
#'   or nonrandom omissions. See *Residual contrasts* in the interval help.
#'
#'   Multiple-score component panels share axes. Compare their magnitudes only
#'   when the score units are comparable. Covariance heatmaps show the selected
#'   submatrix for every source, with variances on the diagonal and covariances
#'   off the diagonal. Covariance units are the product of the two score units;
#'   these are not correlations. All panels share a color scale centered at zero
#'   and include signed numeric labels. Color presets use a diverging palette;
#'   `"monochrome"` uses a gray scale. Values are neither standardized nor
#'   converted to correlations. A negative covariance can be legitimate; a
#'   negative variance is an unconstrained estimate, not a population variance.
#'   `ComponentPSD` and returned notes describe the full source matrix even if
#'   only a subset of scores is shown. A large covariance alone does not diagnose
#'   halo. Cell labels are rounded for display; the table keeps exact values.
#'
#'   With `type = "covariance", intervals = TRUE`, select exactly two scores
#'   and an explicit method. A single horizontal interval plot across sources
#'   replaces the heatmaps. Covariance endpoints remain signed, including
#'   negative values. `"gci"` has approximate coverage and separate Monte Carlo
#'   endpoint bands; `"regression_region"` has conservative normal-model
#'   coverage but may be very wide. Square outlines retain unavailable intervals.
#'   `SEStatus` is separate: a one-df block can prevent an SE while still
#'   permitting a conservative interval. Open circles mark negative covariances,
#'   which can be legitimate. Explicit `"reml_profile"` instead fits a joint
#'   normal model, also for incomplete sources, and retains its separate
#'   `IntervalCenter` and full numerical diagnostics. Its chi-square calibration
#'   is experimental and can under-cover with few facet levels; no sampling SE
#'   is supplied. None of these methods is a correlation interval.
#'
#'   All estimated interactions are shown. The table preserves the original
#'   `Source` identifiers and adds a readable `Label`. For a crossed Person,
#'   Rater and Task design, the last label is "Person x Rater x Task + residual".
#'   Its `Source` remains `Residual`: with one observation per cell, the
#'   highest-order interaction and other residual variation cannot be estimated
#'   separately. A one-facet design instead shows "Person x Task + residual"
#'   (or the declared facet name); a nested design shows, for example,
#'   "Person x Rater(Task) + residual". No extra independent error component
#'   should be added to these combined estimates.
#'
#'   Measurement error is not synonymous with the component called `Residual`.
#'   For the supported random-facet mean-score projections, relative error
#'   includes the person-by-condition interactions and combined residual;
#'   absolute error additionally includes condition main effects and interactions
#'   without Person. The D-study divides each component by its applicable future
#'   counts, then adds these terms. Use `plot(d, type = "sem")` to view the square
#'   roots of the two error variances. Neither SEM nor a residual-component point
#'   is a sampling confidence interval. See [mfrm_d_study()] for
#'   the formulas and assumptions, and `vignette("mfrmr-gtheory")` for a worked
#'   source-to-error example.
#'
#'   This plot reads stored component estimates and preserves missing-data rules.
#'   Explicit REML intervals fit a separate model without replacing those points;
#'   other interval methods use the saved summaries.
#'   Titles and explanations are returned separately from the figure. Use
#'   `plot_data(p, "notes")` and keep them when sharing a figure. Both base R
#'   graphics and optional [as_ggplot()] conversion are supported.
#'
#' @return Invisibly, an `mfrm_plot_data` object with `table` (source, display label, selected
#'   estimate, sign status and full component-matrix PSD diagnostic),
#'   `weights`, `score`, `kind`, `source_design`, `source_rows`, axis labels and
#'   notes. `plot`, `panels` and `ncol` describe the view. Component tables include
#'   `Score`; covariance tables include `Row`, `Column` and `Element` (variance
#'   or covariance), with a shared `fill_limits` and `heat_colors` palette.
#'   With intervals, `table` also includes `SE`, `Lower`, `Upper`, `RawLower`,
#'   `RawUpper`, `Level`, `Method`, `IntervalStatus` and `IntervalAvailable`.
#'   `intervals` retains one complete `confint()` table per selected score or
#'   composite, with its assumptions and source context; `unavailable_intervals`
#'   retains rows without a drawable interval. `Status` still describes the
#'   sign of the point estimate.
#'   REML component tables additionally retain `IntervalCenter`, `SEStatus`,
#'   candidate bounds and `Validation`; full numerical profiles are retained
#'   in `intervals`. Saving the plot preserves them for `as_ggplot()` reuse.
#'   Covariance interval plots have `kind = "Covariance"`, one forest panel,
#'   and the original pair in `score`. The full interval table retains covariance
#'   units, assumptions, numerical bands and SE/interval diagnostics.
#'   Save the result with `saveRDS()` for later `as_ggplot()` conversion,
#'   without recalculating intervals, or save the G-study to redraw it with
#'   base graphics.
#' @seealso [mfrm_g_study()], [plot.mfrm_multivariate_d_study()]
#' @md
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_multivariate_gstudy", run.dontrun = TRUE)
#' \dontrun{
#' tasks <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
#' g <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
#' p <- plot(g, weights = c(V = .5, W = .5), preset = "monochrome",
#'           intervals = TRUE, assumption = "normal")
#' plot_data(p, "table")
#' plot_data(p, "notes")
#' plot(g, score = c("V", "W"))
#' covariance_plot <- plot(g, type = "covariance")
#' plot_data(covariance_plot, "table")
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   as_ggplot(p)
#' }
#' }
#' @export
plot.mfrm_multivariate_gstudy <- function(x, score = NULL, weights = NULL,
                                         draw = TRUE, preset = "standard",
                                         type = c("components", "covariance"),
                                         ncol = NULL, intervals = FALSE,
                                         assumption = NULL, level = .95,
                                         method = "auto", draws = 20000L,
                                         seed = 20261016L, ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  rlang::check_dots_empty()
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw))
    stop("`draw` must be TRUE or FALSE.", call. = FALSE)
  style <- resolve_plot_preset(preset)
  type <- match.arg(type)
  if (!is.logical(intervals) || length(intervals) != 1L || is.na(intervals))
    stop("`intervals` must be TRUE or FALSE.", call. = FALSE)
  if (!intervals && (!missing(assumption) || !missing(level) || !missing(method) || !missing(draws) || !missing(seed)))
    stop("`assumption`, `level`, `method`, `draws` and `seed` apply only when `intervals = TRUE`.", call. = FALSE)
  if (intervals && !identical(method,"gci") && (!missing(draws) || !missing(seed)))
    stop("`draws` and `seed` apply only to `method = 'gci'`.", call. = FALSE)
  covariance_intervals <- intervals && type == "covariance"
  if (!is.null(ncol) && (!is.numeric(ncol) || is.object(ncol) || is.complex(ncol) ||
      !is.null(dim(ncol)) || length(ncol) != 1L || !is.finite(ncol) ||
      ncol < 1 || ncol != floor(ncol)))
    stop("`ncol` must be a positive integer.", call. = FALSE)
  scores <- x$design$scores
  components <- x$components
  if (!inherits(x, "mfrm_multivariate_gstudy") || !is.character(scores) ||
      !length(scores) || anyNA(scores) || anyDuplicated(scores) ||
      !is.list(components) || !length(components) || is.null(names(components)) ||
      anyNA(names(components)) || anyDuplicated(names(components)) ||
      any(!nzchar(names(components))) ||
      !all(vapply(components, function(a) is.matrix(a) && is.numeric(a) &&
        !is.complex(a) && identical(dimnames(a), list(scores, scores)) &&
        all(is.finite(a)) && identical(a, t(a)), logical(1)))) {
    stop("The G-study covariance matrices or score identities are incomplete or altered; recreate the result.", call. = FALSE)
  }
  factors <- names(x$design$counts)
  if (!identical(factors[1L], "Person") || !.mfrm_mvgt_valid_facets(factors[-1L]))
    stop("The G-study factor identities are incomplete or altered; recreate the result.", call. = FALSE)
  nesting <- .mfrm_mvgt_nesting(x$design$nesting, factors[-1L])
  if (!identical(names(components), .mfrm_mvgt_sources(factors, nesting)))
    stop("The G-study component identities do not match its design; recreate the result.", call. = FALSE)
  if (!is.null(score) && !is.null(weights))
    stop("Choose either `score` or `weights`, not both.", call. = FALSE)
  if (type == "covariance" && !is.null(weights))
    stop("Covariance plots use original scores; omit `weights`.", call. = FALSE)
  kind <- if (is.null(weights)) "Score" else "Composite"
  if (is.null(weights)) {
    if (is.null(score)) score <- if (type == "covariance") scores else scores[1L]
    if (!is.character(score) || !is.null(dim(score)) || !length(score) ||
        anyNA(score) || anyDuplicated(score) || !all(score %in% scores))
      stop("`score` must name distinct original scores.", call. = FALSE)
    weights <- diag(length(scores))[, match(score, scores), drop = FALSE]
    dimnames(weights) <- list(scores, score)
    title <- paste("Score:", paste(score, collapse = ", "))
  } else {
    if (!is.numeric(weights) || is.object(weights) || is.complex(weights) ||
        !is.null(dim(weights)) || length(weights) != length(scores) ||
        is.null(names(weights)) || anyNA(names(weights)) || anyDuplicated(names(weights)) ||
        !setequal(names(weights), scores) || any(!is.finite(weights)) || !any(weights != 0))
      stop("`weights` must name every score once with finite values and at least one nonzero value; select one vector for the plot.", call. = FALSE)
    weights <- weights[scores]
    score <- "Composite"
    title <- paste("Composite:", paste(names(weights), format(weights, trim = TRUE), sep = " = ", collapse = ", "))
  }
  diagnostics <- x$component_diagnostics
  psd <- diagnostics$PositiveSemidefinite[match(names(components), diagnostics$Source)]
  if (!is.logical(psd) || length(psd) != length(components) || anyNA(psd))
    stop("The G-study component diagnostics are incomplete; recreate the result.", call. = FALSE)
  source_labels <- gsub(":", " x ", names(components), fixed = TRUE)
  highest <- if (is.null(nesting)) paste(factors, collapse = " x ") else
    paste0("Person x ", names(nesting), "(", nesting, ")")
  source_labels[names(components) == "Residual"] <- paste(highest, "+ residual")
  if (covariance_intervals) {
    if (length(score) != 2L)
      stop("Covariance intervals require `score` to select exactly two original scores.",call.=FALSE)
    args <- list(object=x,score=score,assumption=assumption,level=level,method=method)
    if (identical(method,"gci")) args <- c(args,list(draws=draws,seed=seed))
    ci <- list(do.call(confint.mfrm_multivariate_gstudy,args))
    names(ci) <- paste(score,collapse=" / ")
    table <- data.frame(Source=names(components),Estimate=ci[[1]]$Estimate,
      ComponentPSD=psd,Score=names(ci),Panel=1L)
    panels <- data.frame(Panel=1L,Label=names(ci))
    weights <- NULL; kind <- "Covariance"
    type <- "components" # Reuse the component forest; heatmaps cannot show bounds clearly.
  } else if (type == "covariance") {
    table <- do.call(rbind, lapply(seq_along(components), function(i) {
      cells <- expand.grid(Row = score, Column = score, stringsAsFactors = FALSE)
      cells$Source <- names(components)[i]
      cells$Estimate <- as.vector(components[[i]][score, score, drop = FALSE])
      cells$Element <- ifelse(cells$Row == cells$Column, "Variance", "Covariance")
      cells$ComponentPSD <- psd[i]
      cells$Panel <- i
      cells
    }))
    panels <- data.frame(Panel = seq_along(components), Label = source_labels)
    weights <- NULL
  } else {
    if (is.null(dim(weights))) weights <- matrix(weights, ncol = 1L,
      dimnames = list(names(weights), score))
    table <- do.call(rbind, lapply(seq_along(score), function(i) {
      w <- weights[, i]
      scale <- max(abs(w))
      if (!is.finite(scale^2) || scale^2 == 0)
        stop("Composite variances exceed numeric range; rescale the weights.", call. = FALSE)
      w <- w / scale
      values <- vapply(components, function(a) as.numeric(crossprod(w, a %*% w)) * scale^2, numeric(1))
      if (any(!is.finite(values))) stop("Projected variances exceed numeric range; rescale the scores or weights.", call. = FALSE)
      data.frame(Source = names(components), Estimate = unname(values),
        ComponentPSD = psd, Score = score[i], Panel = i, row.names = NULL)
    }))
    panels <- data.frame(Panel = seq_along(score), Label = score)
    if (ncol(weights) == 1L) weights <- weights[, 1L]
  }
  table$Status <- ifelse(table$Estimate < 0, "Negative estimate", "Nonnegative estimate")
  table$Label <- source_labels[match(table$Source, names(components))]
  if (intervals) {
    if (!covariance_intervals) ci <- lapply(score, function(s) {
      if (kind == "Composite") confint.mfrm_multivariate_gstudy(x,
        weights = weights, assumption = assumption, level = level, method=method) else
        confint.mfrm_multivariate_gstudy(x, score = s,
          assumption = assumption, level = level, method=method)
    })
    if (!covariance_intervals) names(ci) <- score
    bounds <- do.call(rbind, lapply(ci, function(z) {
      columns <- intersect(c("SE", "Lower", "Upper", "RawLower", "RawUpper",
        "Level", "Method", "Status", "SEStatus", "IntervalCenter",
        "CandidateLower", "CandidateUpper", "Validation"), names(z))
      z <- z[match(names(components), z$Source), columns]
      names(z)[names(z) == "Status"] <- "IntervalStatus"
      z
    }))
    rownames(bounds) <- NULL
    table <- cbind(table, bounds)
    table$IntervalAvailable <- is.finite(table$Lower) & is.finite(table$Upper)
  }
  ncol <- min(ncol %||% 2L, nrow(panels))
  note <- c("Raw source-model variance components in squared score units; no sampling intervals. Not future mean-score error contributions or percentages.",
    paste0("Residual = ", highest, " + residual: these sources cannot be separated with one observation per person-by-condition cell."),
    "Person variation is not measurement error. Relative error combines person-by-condition interactions and residual; absolute error also includes condition main effects and interactions without Person. Divide components by their applicable future counts before summing; use the D-study SEM plot for error in mean-score units.",
    "Filled circles: nonnegative estimates; open circles: negative estimates. Weights are used as supplied.")
  if (type == "covariance") {
    note[1L] <- "Raw source-model matrices: diagonal variances in squared score units, off-diagonal covariances in products of the two score units. Not correlations or future mean-score error contributions; no sampling intervals."
    note[4L] <- "All panels share a zero-centered color scale. Cell labels are rounded; the table retains exact signed values. Negative covariance can be legitimate; a negative variance is not a population variance. No standardization or matrix repair is performed."
    note <- c(note, "PSD diagnostics refer to each full source matrix, including scores not selected for display. Large covariance alone does not diagnose halo.")
  }
  if (intervals) {
    note[1L] <- sprintf("Raw source-model variance components with %s%% pointwise normal-random-effects intervals in squared score units. Not future mean-score error contributions or percentages.",
      format(100 * level, trim = TRUE))
    note <- c(note,
      if (identical(method, "reml_profile"))
        "Horizontal bars show explicit approximate REML component intervals. Points remain original ANOVA/MINQUE estimates, including negative estimates outside their intervals; IntervalCenter records the separate REML estimate. SE is not supplied by profile inversion. Intervals are not simultaneous across sources or scores." else
      if (identical(method, "residual_contrast"))
        "Only Residual has a residual-contrast interval. Points retain original ANOVA/MINQUE estimates; IntervalCenter records the separate projected estimate. Other sources remain visible without intervals. SE is not supplied for the original estimator." else
        "Horizontal bars show nonnegative-space intervals: approximate Ting MLS for components; chi-square for the combined residual. They are not simultaneous across sources or scores. Negative points are retained even when outside those intervals.",
      "Square outlines mark estimates without intervals; inspect IntervalStatus and unavailable_intervals. No bar does not mean zero uncertainty. Raw endpoints remain in the table.",
      unique(unlist(lapply(ci, attr, which = "notes"), use.names = FALSE)))
  }
  note <- c(note, "Compare magnitudes only for comparable score units. Shared plot scales do not make different score units equivalent.")
  if (covariance_intervals) note <- c(attr(ci[[1]],"notes"),
    "Horizontal bars show signed covariance intervals for the selected pair. Open circles are negative covariance estimates, not invalid variances; square outlines mark unavailable intervals. Inspect IntervalStatus and SEStatus separately.",
    paste0("Residual = ",highest," + residual: these sources cannot be separated with one observation per cell."))
  negative_variance <- table$Estimate < 0
  if (type == "covariance") negative_variance <- negative_variance & table$Element == "Variance"
  if (covariance_intervals) negative_variance[] <- FALSE
  if (any(negative_variance)) note <- c(note,
    "Negative raw estimates are retained; they are not negative population variances. Review the design and component diagnostics.")
  if (any(!psd)) note <- c(note,
    "Non-PSD component covariance matrices are retained; a positive selected variance does not establish matrix validity.")
  labels <- vapply(source_labels, function(s) paste(strwrap(s, width = 28), collapse = "\n"), character(1))
  heat_colors <- if (style$name == "monochrome") c("#303030", "#B8B8B8", "#FFFFFF") else
    c("#B35806", "#F7F7F7", "#2166AC")
  bound <- max(abs(table$Estimate))
  if (bound == 0) bound <- 1
  out <- new_mfrm_plot_data("multivariate_gstudy", list(table = table, score = score,
    plot = type, panels = panels, ncol = ncol,
    heat_colors = heat_colors, fill_limits = c(-bound, bound),
    kind = kind, weights = weights, source_design = x$design[c("method", "counts", "rows",
      "complete", "balanced", "nesting", "score_convention")], source_rows = x$data_usage$counts,
    title = title, notes = data.frame(Type = "note", Text = note, stringsAsFactors = FALSE),
    labels = labels, x_label = if(covariance_intervals) "Covariance (product of score units)" else "Variance (squared score units)",
    preset = style$name))
  if (intervals) {
    out$data$intervals <- ci
    out$data$unavailable_intervals <- table[!table$IntervalAvailable, , drop = FALSE]
  }
  if (!draw) return(invisible(out))
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  apply_plot_preset(style)
  nr <- ceiling(nrow(panels) / ncol)
  layout <- matrix(c(seq_len(nrow(panels)), rep(0L, nr * ncol - nrow(panels))), ncol = ncol, byrow = TRUE)
  if (type == "covariance") layout <- rbind(layout, rep(nrow(panels) + 1L, ncol))
  # The color key needs room for its axis even when many panels share a device.
  graphics::layout(layout, heights = if (type == "covariance") c(rep(1, nr), graphics::lcm(1.8)) else rep(1, nr))
  graphics::par(cex = 1, mex = 1)
  colors <- grDevices::colorRampPalette(heat_colors)(101)
  limits <- range(c(0, table$Estimate))
  if (intervals) limits <- range(c(limits, table$Lower, table$Upper), na.rm = TRUE)
  if (diff(limits) == 0) limits <- c(-.05, .05)
  for (i in seq_len(nrow(panels))) {
    part <- table[table$Panel == i, , drop = FALSE]
    if (type == "components") {
      left <- max(4, max(graphics::strwidth(labels, units = "inches", cex = .75)) /
        graphics::par("csi") + 2)
      graphics::par(mar = c(4.5, left, if (nrow(panels) > 1L) 2.5 else 1, .5))
      position <- rev(seq_len(nrow(part)))
      plot(limits, c(.5, nrow(part) + .5), type = "n", yaxt = "n",
        xlab = out$data$x_label, ylab = "", cex.lab = .8, cex.axis = .8)
      axis(2, at = position, labels = labels, las = 1, cex.axis = .75, tick = FALSE)
      abline(v = 0, lty = 2, col = style$grid)
      if (intervals) {
        ok <- part$IntervalAvailable
        segments(part$Lower[ok], position[ok], part$Upper[ok], position[ok], col = style$foreground)
        for (endpoint in c("Lower", "Upper")) segments(part[[endpoint]][ok], position[ok] - .09,
          part[[endpoint]][ok], position[ok] + .09, col = style$foreground)
        points(part$Estimate[!ok], position[!ok], pch = 0, cex = 1.5, col = style$foreground)
      } else segments(0, position, part$Estimate, position, col = style$foreground)
      points(part$Estimate, position, pch = ifelse(part$Estimate < 0, 1, 16), col = style$foreground)
    } else {
      graphics::par(mar = c(3.3, 5.5, 3.5, .5))
      q <- length(score)
      plot(c(.5, q + .5), c(.5, q + .5), type = "n", axes = FALSE,
        xlab = "", ylab = "", xaxs = "i", yaxs = "i")
      xx <- match(part$Column, score); yy <- q + 1 - match(part$Row, score)
      index <- round((part$Estimate / bound + 1) * 50) + 1L
      rect(xx - .5, yy - .5, xx + .5, yy + .5, col = colors[index], border = "white")
      ink <- ifelse(colMeans(grDevices::col2rgb(colors[index])) < 140, "white", "black")
      text(xx, yy, trimws(formatC(part$Estimate, format = "g", digits = 3)), col = ink, cex = .8)
      score_labels <- vapply(score, function(s) paste(strwrap(s, 14), collapse = "\n"), character(1))
      axis(1, at = seq_len(q), labels = score_labels, tick = FALSE, cex.axis = .75, gap.axis = -1)
      axis(2, at = rev(seq_len(q)), labels = score_labels, las = 1, tick = FALSE, cex.axis = .75)
    }
    if (nrow(panels) > 1L || type == "covariance") mtext(
      paste(strwrap(panels$Label[i], 28), collapse = "\n"), side = 3, line = .5, cex = .85)
  }
  if (type == "covariance") {
    graphics::par(mar = c(2.7, 5.5, .2, 1), pty = "m")
    image(seq(-1, 1, length.out = 101), 1, matrix(seq_len(101), ncol = 1),
      col = colors, axes = FALSE, xlab = "", ylab = "")
    axis(1, at = c(-1, 0, 1), labels = format(signif(c(-bound, 0, bound), 3), trim = TRUE), cex.axis = .8)
    mtext("Variance / covariance", side = 1, line = 1.6, cex = .8)
  }
  invisible(out)
}
