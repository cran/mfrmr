#' Compare measurement-condition counts in a multivariate D-study
#'
#' See how a planned change in the number of tasks, raters, or other conditions affects the
#' dependability of a mean score. Plot an existing D-study result; no model
#' is refitted. Start with `plot(d)` for G and Phi, then use
#' `plot(d, type = "sem")` to examine error in score units. If the result has
#' several composites, select them explicitly with `composite`.
#'
#' @param x A result from [mfrm_d_study()].
#' @param type `"coefficients"` shows G and Phi in separate panels;
#'   `"sem"` shows relative and absolute standard errors of measurement;
#'   `"error_components"` shows each source's contribution to planned error
#'   variance, with separate panels for selected scores and error definitions.
#' @param score One or more distinct original score names, in legend order.
#'   An explicit name always selects
#'   an original score, even if a composite has the same name. If both `score`
#'   and `composite` are `NULL`, select the sole composite, or the first score
#'   when there are no composites. Several composites require an explicit
#'   selection. [plot_data()] identifies the selection and composite weights;
#'   a suggested title is returned rather than drawn on the figure.
#' @param x_var A count column in `x$design_grid`: `"Tasks"`/`"Raters"` for
#'   the task/rater interface, or a facet label supplied through `facets`.
#'   By default use the last varying count column, or the last column if
#'   none varies. A one-facet result permits only its included facet.
#'   Each line holds the other facet count constant.
#' @param draw Draw the plot when `TRUE`; `FALSE` only returns its data.
#' @param preset Plot style: `"standard"`, `"publication"`, `"compact"`, or
#'   `"monochrome"`. The first three use distinct curve colors; monochrome
#'   uses one dark color. Point shapes and line types also distinguish curves.
#'   An explicit preset overrides `options(mfrmr.plot_preset = ...)`.
#'   [as_ggplot()] retains colors saved in a plot payload.
#' @param composite One or more distinct composite names from `x$coefficients`, as defined
#'   by a weight-matrix column, or `"Composite"` for vector weights. Supply
#'   either `score` or `composite`, not both.
#' @param metric Optional single metric compatible with `type`: `"G"` or
#'   `"Phi"` for coefficients; `"RelativeSEM"` or `"AbsoluteSEM"` for SEM.
#'   For error components, select `"RelativeErrorVariance"` or
#'   `"AbsoluteErrorVariance"`.
#'   `NULL` shows both metrics in separate panels.
#' @param reference Optional numeric reference values: between 0 and 1 for
#'   coefficients, nonnegative for SEM. `NULL` draws no reference lines.
#'   The same values apply to each panel. These are user-chosen guides, not
#'   universal acceptance rules or confidence bounds.
#' @param panel_by Optional non-horizontal count column in `x$design_grid`,
#'   such as `"Raters"` when `x_var = "Tasks"`. Each count gets a separate
#'   panel; selected scores/composites keep the same color, line type and point
#'   shape across panels. `NULL` retains the overlaid curves.
#' @param ncol Optional positive integer number of panel columns. Requires
#'   `panel_by` for curves; by default use up to two columns.
#'   Error-component plots use one column by default and do not use `panel_by`.
#'   With both curve metrics selected,
#'   panels run through all counts for the first metric, then the second.
#' @param intervals `FALSE` (default) draws point curves. `TRUE` computes
#'   intervals and adds ribbons in one call; also specify
#'   `assumption = "normal"` and select one `metric`: `"G"`, `"Phi"`, or,
#'   with `type = "sem"`, `"RelativeSEM"` or `"AbsoluteSEM"`.
#'   Alternatively, supply a saved `confint(x)` result, or a list with one
#'   result per selected score/composite, to avoid recalculation. Saved results
#'   must use one method/level and match the source, weights and design grid;
#'   retain unavailable rows. `NULL` is equivalent to `FALSE`.
#'   Use `plot(ci)` separately to compare several interval methods.
#' @param assumption With `intervals = TRUE`, explicitly specify `"normal"`
#'   for independent normal random effects; this assumption is not tested.
#' @param method With `intervals = TRUE`, the interval method passed to
#'   [confint.mfrm_multivariate_d_study()]. `"auto"` selects exact F for
#'   one-facet G and generalized CI for one-facet Phi, also for nested plans
#'   with source child counts per parent unchanged. For two crossed random
#'   facets, `"auto"` selects approximate `"mls"` null inversion for G/Phi.
#'   Undefined MLS widths are conservatively retained and flagged in the saved
#'   interval result. For alternatives or changed nested child counts, select `"approximate_f"` or
#'   `"mean_square_region"`. For SEM, `"auto"` uses `"mls"` for supported
#'   balanced crossed/nested sources: chi-square for a single mean square,
#'   Graybill-Wang for positive sums, Ting for signed sums. Two-facet SEM also
#'   accepts `"mean_square_region"` for conservative normal-model intervals.
#'   `"reml_profile"` explicitly fits approximate G/Phi intervals for balanced
#'   or identifiable incomplete sources; lme4 and minqa are required for fitting.
#'   It preserves the original curve points and reports separate REML estimates
#'   in the saved intervals. Its finite-sample/boundary coverage is not guaranteed.
#'   For one-facet `RelativeSEM`, explicit `"residual_contrast"` accepts
#'   retained ANOVA/MINQUE ratings, including incomplete sources. It transforms
#'   a residual chi-square interval without likelihood fitting, retaining
#'   original curve points and a separate projected center in the saved table.
#'   For one-facet G, explicit `"spectral_f"` accepts retained ANOVA/MINQUE
#'   ratings, including incomplete sources, using a normal-model variance-ratio
#'   F pivot. It keeps the original points; sparse sources can yield wide bands.
#'   Balanced crossed G also accepts explicit `"score"`, an experimental
#'   extended-nuisance score interval calculated from stored mean squares.
#'   It preserves the source points; empty sets and numerical failures remain
#'   missing ribbon segments, with their reasons in the saved interval object.
#'   No method is
#'   silently substituted after failure.
#' @param level Confidence level for `intervals = TRUE`, default 0.95.
#' @param ... Reserved for future use; additional arguments are rejected.
#'
#' @details
#' G describes consistency of relative ordering: for example, ranking
#' examinees. Phi also includes shifts caused by easier tasks or more lenient
#' raters and concerns absolute score levels. Higher coefficients indicate
#' greater dependability under the fitted model; neither gives the probability
#' of a correct pass/fail decision. No universal acceptable cutoff is drawn.
#'
#' SEM is the square root of error variance. Lower values indicate less error
#' in the units of the selected mean score or composite. Relative SEM concerns
#' ordering; absolute SEM also includes condition-wide shifts. SEM is not a
#' confidence interval for G or Phi. The two SEM panels share a scale, but
#' different scores or differently scaled composites need not have comparable
#' units. Weights are not normalized. Changing weights can also change the
#' meaning of the score; a higher G or Phi alone does not justify a new set
#' of weights.
#'
#' With `type = "error_components"`, compare how complete plans change the
#' sources of error. Each point is a projected variance component divided by
#' its planned replication count. Relative error includes Person interactions
#' and residual; absolute error also includes condition-wide shifts.
#' For two crossed facets, Person-by-rater is divided by the rater count,
#' Person-by-task by the task count, and the highest interaction/residual by
#' their product. In a nested design, child counts are per parent; the parent
#' and Person-by-parent terms use the parent count, and child-within-parent
#' and residual terms use the product. See [mfrm_d_study()] for the formulas.
#' Sources absent from an error definition are omitted, including the Person
#' universe component.
#'
#' Contributions add to error **variance**, in squared score units, not to
#' SEM. Composite projections use the full covariance matrices and the saved
#' weights. Negative component estimates remain visible even when the total
#' error variance is negative and SEM is unavailable; they do not represent
#' negative population variances. A negative contribution is not evidence
#' that adding observations causes more population error. The plot compares
#' point projections, not sampling uncertainty or an optimal allocation.
#' Do not pass curve-only `x_var`, `panel_by`, `reference`, or `intervals`.
#' Use a small set of relevant plans, and increase figure height for many
#' scores or plans. Stored components and plans are checked against the saved
#' D-study estimates without refitting the G-study.
#'
#' Points represent only the supplied scenarios. Lines connect them as visual
#' guides for each selected score at a fixed count of the other facet, without fitting a curve or
#' evaluating intermediate designs. For example, with tasks on the horizontal
#' axis, compare points along one line to change tasks while keeping raters
#' constant. For overlaid curves, each legend entry identifies the score and
#' the other facet count. With `panel_by`, the panel identifies the count and
#' the legend identifies the score. Compare rater counts within the same score.
#'
#' Multivariate means several score criteria, not a requirement for a 3D plot.
#' Use these 2D panels to compare counts while holding the other facet fixed.
#' With many lines, use `panel_by` to separate the other facet count, or
#' specify a smaller set of relevant plans when creating the D-study. Panel
#' axes share limits; only requested scenarios are drawn, even in a
#' nonrectangular grid. Panel labels identify counts (and metrics when needed),
#' while the figure title and cautions remain in the returned data. Increase
#' figure height for many panels or select one `metric`. This method does not draw
#' a heatmap or 3D surface. For uncertainty in a prespecified change, use
#' [mfrm_multivariate_d_compare()] and its plot method; an SEM curve is not
#' a confidence band. [plot_data()] supplies exact values, including unavailable
#' estimates, and `preset = "monochrome"` retains line/point distinctions.
#'
#' Coefficients also depend on variation among persons. If universe-score
#' variance increases while error variance stays constant, G and Phi increase
#' but SEM does not decrease. Therefore, comparing coefficients across different
#' populations does not isolate improvement in the assessment. Read the SEMs,
#' score units and target population alongside the coefficients.
#'
#' Scenarios retain the G-study's crossed or nested structure, with complete
#' balanced future conditions shared by every person and score. For a nested
#' facet, the axis or legend says "per" to identify its count within each parent,
#' for example "Raters per Task"; this is not the total rater pool size.
#' Even when the G-study used incomplete
#' data, the figure does not describe reliability of that sparse roster.
#' Point projections use the estimated covariance components. By default no
#' sampling intervals are shown. `intervals = TRUE` uses
#' [confint.mfrm_multivariate_d_study()] for each selected score/composite,
#' without changing the G-study. `method = "reml_profile"` fits a separate
#' likelihood model once per selected score/composite. `"residual_contrast"`
#' projects the retained ratings for one-facet relative SEM, and `"spectral_f"`
#' does so for one-facet G; other methods
#' use saved ANOVA summaries. To customize other inference settings such
#' as pivotal `draws` or `seed`, pass saved `confint()` results instead.
#' This alternative does not recalculate intervals. Source summaries, score weights, plans and stored point estimates
#' are checked before overlaying. Older interval objects without source summaries
#' must be recreated with `confint(x)`; their forest plots still work.
#' Ribbons connect the supplied endpoints; intermediate plans are not evaluated.
#' Delta normal-reference ribbons may extend outside `[0,1]`; both renderers
#' expand the axis rather than clipping the saved limits. This indicates an
#' approximation limitation, not a population coefficient outside its range.
#' Missing intervals break the ribbon and leave open points when estimates exist.
#' An isolated interval is a vertical segment. An interval can remain visible
#' without a point estimate. Boundary line types distinguish overlapping bands
#' in monochrome. Confidence levels, method limitations, numerical diagnostics
#' and unavailable reasons remain in the returned object, not figure footnotes.
#'
#' Ribbons do not establish simultaneous coverage across scores or metrics.
#' Joint coverage across future plans, if any, depends on the chosen method
#' and its model assumptions; drawing a band supplies no additional guarantee.
#' Overlap is not a test of score or plan differences. Measurement SEM and
#' paired-difference intervals cannot be used as coefficient ribbons.
#' SEM ribbons show sampling uncertainty in the population SEM, not intervals
#' for Person scores. Their axes retain score units and include the entire
#' interval. MLS bands are approximate; do not read a narrow ribbon as a
#' guarantee of coverage. See [confint.mfrm_multivariate_d_study()] for the
#' equations, small-sample evaluation and boundary limitations.
#' A flattening curve can indicate limited gains from
#' adding one facet alone; it does not identify a cost-effective optimum.
#' The highest point alone does not account for uncertainty in the estimated
#' components. Read the size of differences as well as their ordering. Equal products
#' of rater and task counts give equal rating counts, not necessarily equal
#' examinee burden or total cost. Unavailable candidates prevent a complete
#' comparison for that metric.
#'
#' Unavailable estimates stay missing. A panel with no estimates states the
#' reason; otherwise returned notes report unavailable scenarios. The returned
#' `unavailable` table preserves their counts, metrics, and `Status` values.
#' Each panel uses its metric's status, so an unavailable Phi does not hide
#' an available G. Non-PSD component matrices are described in the returned
#' notes; inspect `x$component_diagnostics` before using raw projections.
#' Increasing planned counts does not repair those component estimates.
#' The saved `negative_components` convention is preserved in curves and
#' error-contribution plots. Notes identify raw, GENOVA ALGORITHM or GENOVA EMS
#' projections; original G-study diagnostics remain available. GENOVA-adjusted
#' results accept original-source population intervals calculated for that
#' treatment by `confint(x)`. Bands retain their original bounds and may not
#' contain adjusted points; notes and saved intervals retain the distinction.
#' A saved interval with a different component treatment is rejected.
#' Missing values are never plotted as zero. Previously saved D-study objects
#' retain their recorded values and availability: rerun
#' [mfrm_d_study()] on the saved G-study to apply current rules.
#'
#' @return Invisibly, an `mfrm_plot_data` object. Use [plot_data()] to extract
#'   `table` (the selected scores/composites), `series` (one row per scenario and
#'   plotted metric, including missing values), `unavailable`, `design_grid`,
#'   `weights`, score identity, axis/group names, labels, `component_note`,
#'   and legend settings. In `series`, `Status` is specific to `Metric`;
#'   `table` retains the original row status and metric-specific columns.
#'   `weights` contains a named vector for one composite, a named matrix for
#'   several composites, or `NULL` for original scores. `reference_lines` and
#'   `notes` retain the chosen guides and their interpretation. When `panel_by`
#'   is supplied, `panels` maps `series$Panel` to its metric, count and label;
#'   `ncol` records the resolved layout for base and ggplot rendering.
#'   With `intervals`, `series` also retains `Lower`, `Upper`, `IntervalStatus`,
#'   `IntervalAvailable` and `IntervalRun` (connected ribbon segments).
#'   `intervals` retains the complete input objects, including diagnostics;
#'   `unavailable_intervals` retains missing bounds separately from missing
#'   point estimates. Save the complete plot object with `saveRDS()` to replay.
#'   For `"error_components"`, `table` and `series` contain the same long
#'   contribution table: `Source`, planned `Multiplier`, projected
#'   `ComponentVariance`, `Contribution`, `DiagonalContribution` and
#'   `CovarianceContribution`. The latter two separate within-score variance
#'   and cross-score covariance terms; they sum to `Contribution`.
#'   `TotalErrorVariance` and `SEMStatus` retain the planned total and its
#'   availability; `totals` retains all original D-study metrics and statuses
#'   for selected scores. `panels`, `legend`, `weights`, `design_grid`,
#'   `component_diagnostics` and `notes` preserve the interpretation.
#'   [as_ggplot()] preserves these comparisons for editing or export with
#'   the optional ggplot2 package; for example, `as_ggplot(d, type = "sem")`.
#' @seealso [mfrm_d_study()], [mfrm_g_study()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_multivariate_d_study", run.dontrun = TRUE)
#' \dontrun{
#' # Question: how much would doubling the common tasks change dependability?
#' tasks <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
#' g <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
#' d <- mfrm_d_study(g, data.frame(Tasks = c(3, 6, 12)),
#'   weights = c(V = -1, W = 1))
#' plot(d) # Default: the difference W minus V; higher G/Phi is better.
#' plot(d, type = "sem") # Error in difference-score units; lower is better.
#' plot(d, score = "V") # Inspect an original score separately.
#' # Add intervals for both criteria in the same call.
#' ribbons <- plot(d, score = c("V", "W"), metric = "Phi", reference = 0.8,
#'   intervals = TRUE, assumption = "normal")
#' plot_data(ribbons, "intervals") # Saved intervals and their diagnostics.
#' # 0.8 above is an illustrative user criterion, not a recommended cutoff.
#' values <- plot_data(plot(d, draw = FALSE))
#' values$table # Exact values and availability, rather than reading off a line.
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   p <- as_ggplot(d, type = "sem")
#'   print(p)
#'   # ggplot2::ggsave("d-study-sem.png", p, width = 7, height = 7, dpi = 300)
#' }
#'
#' # A weight matrix names several choices; explicitly select the plotted one.
#' alternatives <- mfrm_d_study(g, data.frame(Tasks = c(3, 6, 12)),
#'   weights = cbind(Equal = c(V = 0.5, W = 0.5), Difference = c(V = -1, W = 1)))
#' plot(alternatives, composite = "Equal")
#' plot(alternatives, composite = "Difference", type = "sem")
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   print(as_ggplot(alternatives, composite = "Equal"))
#' }
#' # Which sources explain error in the difference score?
#' sources <- plot(alternatives, composite = "Difference", type = "error_components")
#' plot_data(sources, "table") # Includes signed cross-score covariance terms.
#' plot_data(sources, "totals") # Error totals, SEM, G/Phi and availability.
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @md
#' @export
plot.mfrm_multivariate_d_study <- function(x, type = c("coefficients", "sem", "error_components"),
                                          score = NULL, x_var = NULL, draw = TRUE,
                                          preset = "standard", composite = NULL,
                                          metric = NULL, reference = NULL,
                                          panel_by = NULL, ncol = NULL,
                                          intervals = FALSE, assumption = NULL,
                                          method = "auto", level = .95, ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (!identical(.mfrm_mvds_policy(x), "raw")) x <- .mfrm_mvds_replay(x)
  if (!isTRUE(intervals) && (!missing(assumption) || !missing(method) || !missing(level)))
    stop("`assumption`, `method` and `level` apply only when `intervals = TRUE`; saved intervals retain their own settings.", call. = FALSE)
  if (isFALSE(intervals)) intervals <- NULL
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw)) {
    stop("`draw` must be TRUE or FALSE.", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  tab <- x$coefficients
  if (!is.null(score) && !is.null(composite)) {
    stop("Choose either `score` or `composite`, not both.", call. = FALSE)
  }
  if (is.null(score) && is.null(composite)) {
    composites <- unique(tab$Score[tab$Kind == "Composite"])
    if (length(composites) > 1L) {
      stop("Several composites are available; select one with `composite`, or an original score with `score`.", call. = FALSE)
    }
    if (length(composites)) composite <- composites else score <- x$gstudy$design$scores[1L]
  }
  if (!is.null(composite)) {
    if (!is.character(composite) || !length(composite) || anyNA(composite) ||
        anyDuplicated(composite) || !all(composite %in% tab$Score[tab$Kind == "Composite"])) {
      stop("`composite` must name one composite or several distinct composites in the D-study result.", call. = FALSE)
    }
    kind <- "Composite"
    score <- composite
  } else {
    kind <- "Score"
    if (!is.character(score) || !length(score) || anyNA(score) ||
        anyDuplicated(score) || !all(score %in% tab$Score[tab$Kind == kind])) {
      stop("`score` must name one original score or several distinct original scores; use NULL for the default selection.", call. = FALSE)
    }
  }
  selected_weights <- if (kind == "Composite") {
    if (is.matrix(x$weights)) {
      if (length(score) == 1L) setNames(x$weights[, score], rownames(x$weights)) else
        x$weights[, score, drop = FALSE]
    } else x$weights
  } else NULL
  tab <- tab[tab$Kind == kind & tab$Score %in% score, , drop = FALSE]
  if (type == "error_components") {
    if (!is.null(x_var) || !is.null(panel_by) || !is.null(reference) || !is.null(intervals))
      stop("`error_components` uses complete plans; omit `x_var`, `panel_by`, `reference` and `intervals`.", call. = FALSE)
    return(.mfrm_mvds_error_plot(x, score, kind, selected_weights, metric, ncol, style, draw))
  }
  facets <- names(x$design_grid)
  if (is.null(x_var)) {
    varying <- facets[vapply(tab[facets], function(n) length(unique(n)) > 1L, logical(1))]
    x_var <- tail(if (length(varying)) varying else facets, 1L)
  }
  if (!is.character(x_var) || length(x_var) != 1L || is.na(x_var) || !x_var %in% facets) {
    stop("`x_var` must name a facet count present in the D-study: ",
      paste(facets, collapse = ", "), ".", call. = FALSE)
  }
  x_label <- if (x_var %in% c("Tasks", "Raters")) paste("Number of", tolower(x_var)) else
    paste(x_var, "count")
  nesting <- x$gstudy$design$nesting
  count_labels <- setNames(facets, facets)
  if (!is.null(nesting)) {
    child_column <- x$gstudy$design$count_columns[[names(nesting)]]
    count_labels[child_column] <- paste(child_column, "per", unname(nesting))
    if (x_var == child_column) x_label <- paste(x_label, "per", unname(nesting))
  }
  group_var <- setdiff(facets, x_var)
  faceted <- !is.null(panel_by)
  if (faceted && (!is.character(panel_by) || length(panel_by) != 1L ||
      is.na(panel_by) || !panel_by %in% group_var)) {
    stop("`panel_by` must name a non-horizontal facet count in the D-study; ",
      "it cannot equal `x_var`.", call. = FALSE)
  }
  if (!is.null(ncol) && (!faceted || !is.numeric(ncol) || is.complex(ncol) ||
      length(ncol) != 1L || !is.finite(ncol) || ncol < 1 || ncol != floor(ncol))) {
    stop("`ncol` must be a positive integer and requires `panel_by`.", call. = FALSE)
  }
  metrics <- if (type == "coefficients") c("G", "Phi") else c("RelativeSEM", "AbsoluteSEM")
  labels <- if (type == "coefficients") c("G: relative ordering", "Phi: absolute score levels") else
    c("Relative SEM: ordering", "Absolute SEM: score levels")
  if (!is.null(metric)) {
    if (!is.character(metric) || length(metric) != 1L || is.na(metric) || !metric %in% metrics) {
      stop("`metric` must select one metric compatible with `type`: ",
        paste(metrics, collapse = ", "), ".", call. = FALSE)
    }
    labels <- labels[match(metric, metrics)]
    metrics <- metric
  }
  references <- .mfrm_d_study_references(reference, type == "coefficients")
  counts <- if (length(group_var)) sort(unique(tab[[group_var]])) else 1
  count_groups <- if (length(group_var)) paste(count_labels[[group_var]], "=", counts) else
    paste(x_var, "only")
  multiple <- length(score) > 1L
  groups <- if (!multiple) count_groups else if (!length(group_var)) score else
    as.vector(t(outer(score, count_groups, paste, sep = " | ")))
  groups <- make.unique(groups)
  group_index <- (match(tab$Score, score) - 1L) * length(counts) +
    if (length(group_var)) match(tab[[group_var]], counts) else 1L
  if (faceted) {
    groups <- score
    group_index <- match(tab$Score, score)
  }
  show_legend <- multiple || (!faceted && length(group_var) > 0L)
  panels <- data.frame(Panel = seq_len(length(metrics)), Metric = metrics,
    Label = labels, stringsAsFactors = FALSE)
  if (faceted) {
    panels <- panels[rep(seq_along(metrics), each = length(counts)), , drop = FALSE]
    rownames(panels) <- NULL
    panels$Panel <- seq_len(nrow(panels))
    panels$Count <- rep(counts, times = length(metrics))
    panels$Label <- if (length(metrics) == 1L) count_groups else
      paste(panels$Label, rep(count_groups, times = length(metrics)), sep = "\n")
  }
  panel_columns <- if (faceted) min(ncol %||% 2, nrow(panels)) else 1L
  y_label <- if (length(metrics) == 1L) {
    if (type == "sem") sub(":", " (score units):", labels) else labels
  } else if (type == "coefficients") "Dependability (higher is better)" else
    "SEM in score units (lower is better)"
  series <- do.call(rbind, lapply(metrics, function(metric) {
    rows <- data.frame(tab, Metric = metric, X = tab[[x_var]],
      Group = groups[group_index],
      Value = tab[[metric]], row.names = NULL, check.names = FALSE)
    if (faceted) rows$Panel <- (match(metric, metrics) - 1L) * length(counts) +
      match(tab[[panel_by]], counts)
    # Preserve the recorded policy for older D-study objects; do not recalculate.
    if (!is.null(tab[[paste0(metric, "Status")]])) rows$Status <- tab[[paste0(metric, "Status")]]
    rows
  }))
  if (isTRUE(intervals)) {
    if (length(metrics) != 1L)
      stop("With `intervals`, select one coefficient or SEM using `metric`.", call. = FALSE)
    if (!identical(assumption, "normal"))
      stop("Interval computation requires `assumption = 'normal'`; the assumption is not tested.", call. = FALSE)
    if (length(method) != 1L)
      stop("Select one interval method for curves; use plot(ci) for a method comparison.", call. = FALSE)
    intervals <- lapply(score, function(s) {
      args <- list(object = x, parm = metrics, assumption = assumption, method = method, level = level)
      args[[if (kind == "Score") "score" else "composite"]] <- s
      do.call(stats::confint, args)
    })
  }
  bands <- if (!is.null(intervals))
    .mfrm_mvds_curve_intervals(x, series, intervals, metrics, kind, score) else NULL
  if (!is.null(bands)) series <- bands$series
  unavailable <- series[series$Status != "Available" | !is.finite(series$Value), , drop = FALSE]
  cols <- if (style$name == "monochrome") rep(style$foreground, length(groups)) else
    grDevices::hcl.colors(max(3L, length(groups)), "Dark 3")[seq_along(groups)]
  line_types <- ((seq_along(groups) - 1L) %/% 6L + seq_along(groups) - 1L) %% 6L + 1L
  shapes <- rep(c(16, 17, 15, 18, 1, 2), length.out = length(groups))
  title <- if (multiple) paste(paste0(kind, "s:"), paste(score, collapse = ", ")) else if (kind == "Composite") paste0("Composite",
    if (is.matrix(x$weights)) paste0(" ", score), ": ", paste(names(selected_weights),
    format(selected_weights, trim = TRUE), sep = " = ", collapse = ", ")) else paste("Score:", score)
  design_label <- if (is.null(nesting)) "Future complete crossed designs" else
    paste0("Future balanced nested designs: ", names(nesting), " within ", unname(nesting))
  subtitle <- paste(design_label, "|", if (is.null(bands))
    "Estimated components held fixed | No confidence intervals" else
    "Saved coefficient intervals | Source-study uncertainty")
  component_note <- c(if (!identical(.mfrm_mvds_policy(x), "raw")) .mfrm_mvds_component_note(x),
    if (any(!x$component_diagnostics$PositiveSemidefinite))
      "Non-PSD covariance components in the original G-study; review component diagnostics and the selected treatment.")
  out <- new_mfrm_plot_data("multivariate_d_study", list(table = tab, series = series,
    unavailable = unavailable, design_grid = x$design_grid, weights = selected_weights,
    negative_components = .mfrm_mvds_policy(x),
    nesting = nesting, panel_by = panel_by,
    panels = if (faceted) panels else NULL, ncol = if (faceted) panel_columns else NULL,
    score = score, kind = kind, plot = type, x_var = x_var, x_label = x_label,
    group_var = group_var, metric = metrics, show_legend = show_legend, y_label = y_label,
    cautions = if (multiple) paste("Coefficients do not establish equivalent score meaning.",
      "Compare SEMs only when score units are comparable; weights are not normalized."),
    reference_lines = references, interpretation_note = .mfrm_d_study_reference_note(references),
    title = title, subtitle = subtitle, component_note = component_note, panel_labels = labels,
    caption = paste("Points are requested scenarios; lines are guides.",
      if (nrow(unavailable)) sprintf("Unavailable estimates in %d scenario(s); inspect the D-study table for reasons.",
        length(unique(unavailable$Scenario)))),
    legend = new_plot_legend(groups, rep(if (faceted) "score" else if (multiple) "score_condition" else "condition_count", length(groups)),
      rep("color", length(groups)), cols),
    line_types = line_types, point_shapes = shapes, preset = style$name))
  if (!is.null(bands)) {
    out$data$intervals <- bands$intervals
    out$data$unavailable_intervals <- series[!series$IntervalAvailable, , drop = FALSE]
    out$data$notes <- rbind(out$data$notes,
      data.frame(Type = "uncertainty", Text = bands$notes))
  }
  if (!draw) return(invisible(out))
  old <- graphics::par()[c("mfrow", "mar", "oma", "cex", "mex")]
  on.exit(graphics::par(old), add = TRUE)
  apply_plot_preset(style)
  graphics::par(mfrow = c(ceiling(nrow(panels) / panel_columns), panel_columns),
    mar = if (faceted) c(4.3, 5.5, 3.2, 1) else c(4.7, 5.3, 2.3, if (show_legend) 8 else 1),
    oma = c(if (faceted && show_legend) 1 + ceiling(length(groups) / 2) else .2, 0, .2, 0))
  if (show_legend && !faceted) {
    margin <- graphics::par("mar")
    legend_width <- max(graphics::strwidth(paste("MMMM", groups), units = "inches", cex = 0.8))
    margin[4L] <- max(margin[4L], legend_width / (graphics::par("csi") * graphics::par("mex")) + 1)
    graphics::par(mar = margin)
  }
  y_lim <- .mfrm_mvds_y_limits(series, type == "coefficients", references$value)
  for (i in seq_len(nrow(panels))) {
    s <- series[if (faceted) series$Panel == panels$Panel[i] else
      series$Metric == panels$Metric[i], , drop = FALSE]
    ok <- s$Status == "Available" & is.finite(s$Value)
    x_lim <- if (faceted) range(series$X) else range(s$X)
    plot(x_lim, y_lim, type = "n", ylim = y_lim,
      xlab = x_label,
      ylab = if (length(metrics) == 1L) sub(": ", "\n", y_label) else if (type == "coefficients") "Dependability\nHigher is better" else "SEM (score units)\nLower is better",
      main = if (faceted || length(metrics) > 1L) panels$Label[i] else NULL, xaxt = "n", yaxt = "n")
    ticks <- sort(unique(if (faceted) series$X else s$X))
    if (length(ticks) > 15L) ticks <- pretty(range(ticks))
    axis(1, at = ticks[ticks >= x_lim[1] & ticks <= x_lim[2] & ticks == floor(ticks)])
    has_intervals <- !is.null(bands) && any(s$IntervalAvailable)
    if (any(ok) || has_intervals) {
      if (type == "coefficients" && identical(y_lim, c(0, 1)))
        axis(2, at = seq(0, 1, 0.2), las = 1) else axis(2, las = 1)
      graphics::grid(nx = NA, ny = NULL, col = style$grid)
      if (has_intervals) for (run in unique(stats::na.omit(s$IntervalRun))) {
        z <- s[!is.na(s$IntervalRun) & s$IntervalRun == run, , drop = FALSE]
        z <- z[order(z$X, z$Scenario), , drop = FALSE]
        j <- match(z$Group[1L], groups)
        if (length(unique(z$X)) == 1L) {
          segments(z$X, z$Lower, z$X, z$Upper, col = cols[j], lty = line_types[j])
        } else {
          graphics::polygon(c(z$X, rev(z$X)), c(z$Lower, rev(z$Upper)),
            col = grDevices::adjustcolor(cols[j], alpha.f = .14), border = NA)
          lines(z$X, z$Lower, col = cols[j], lty = line_types[j], lwd = .7)
          lines(z$X, z$Upper, col = cols[j], lty = line_types[j], lwd = .7)
        }
      }
      for (j in seq_along(groups)) {
        z <- s[s$Group == groups[j], , drop = FALSE]
        z <- z[order(z$X, z$Scenario), , drop = FALSE]
        z$Value[z$Status != "Available" | !is.finite(z$Value)] <- NA_real_
        lines(z$X, z$Value, col = cols[j], lty = line_types[j], lwd = 2)
        pch <- rep(shapes[j], nrow(z))
        if (!is.null(bands)) pch[!z$IntervalAvailable] <- c(1, 2, 0, 5, 1, 2)[(j-1L) %% 6L+1L]
        points(z$X, z$Value, col = cols[j], pch = pch)
      }
    } else {
      text(mean(x_lim), mean(y_lim), .mfrm_mvds_unavailable_message(s$Status, compact = faceted),
        cex = 0.8, col = style$foreground)
    }
    if (nrow(references)) abline(h = references$value,
      lty = references$linetype, col = style$neutral)
    if (show_legend && !faceted) {
      .mfrm_mvds_legend(groups, cols, line_types, shapes, faceted = FALSE)
    }
  }
  if (faceted && show_legend) {
    .mfrm_mvds_legend(groups, cols, line_types, shapes, faceted = TRUE)
  }
  invisible(out)
}

.mfrm_mvds_legend <- function(groups, cols, line_types, shapes, faceted) {
  # Recompute physical text spacing when pkgdown or a saved plot changes devices.
  grDevices::recordGraphics({
    if (faceted) {
      legend(graphics::grconvertX(.5, from = "ndc"),
        graphics::grconvertY(.015, from = "ndc"), legend = groups,
        col = cols, lty = line_types, pch = shapes, ncol = min(2, length(groups)),
        bty = "n", xpd = NA, xjust = .5, yjust = 0, cex = .8)
    } else {
      usr <- graphics::par("usr")
      legend(usr[2] + 0.03 * diff(usr[1:2]), usr[4], legend = groups,
        col = cols, lty = line_types, pch = shapes, bty = "n", xpd = NA, cex = .8)
    }
  }, list(groups = groups, cols = cols, line_types = line_types,
    shapes = shapes, faceted = faceted), asNamespace("mfrmr"))
}

.mfrm_mvds_y_limits <- function(series, coefficients, references) {
  values <- c(0, if (coefficients) 1, references,
    series$Value[series$Status == "Available"],
    series$Lower[series$IntervalAvailable], series$Upper[series$IntervalAvailable])
  limits <- range(values[is.finite(values)])
  if (diff(limits) == 0) c(0, 1) else limits
}

.mfrm_mvds_curve_intervals <- function(x, series, intervals, metrics, kind, scores) {
  if (length(metrics) != 1L)
    stop("With `intervals`, select one coefficient or SEM using `metric`.", call. = FALSE)
  if (inherits(intervals, "mfrm_dstudy_intervals")) intervals <- list(intervals)
  if (!is.list(intervals) || !length(intervals) ||
      !all(vapply(intervals, inherits, logical(1), "mfrm_dstudy_intervals")))
    stop("`intervals` must be a confint(d) result or a list of those results.", call. = FALSE)
  selected <- character(length(intervals))
  series$Lower <- series$Upper <- NA_real_
  series$IntervalStatus <- NA_character_
  series$IntervalAvailable <- FALSE
  for (i in seq_along(intervals)) {
    ci <- intervals[[i]]
    .mfrm_mvds_interval_source(ci, x)
    if (!isTRUE(all.equal(attr(ci, "design_grid"), x$design_grid, tolerance = 0)) ||
        !all(ci$Metric == metrics) || !all(ci$Kind == kind) ||
        length(unique(ci$Score)) != 1L || !ci$Score[1L] %in% scores)
      stop("Intervals must match the plotted plans, metric and score/composite.", call. = FALSE)
    if (length(unique(ci$MethodCode)) != 1L)
      stop("Select one interval method for curves; use plot(ci) for a method comparison.", call. = FALSE)
    if (anyDuplicated(ci$Scenario) ||
        !setequal(ci$Scenario, seq_len(nrow(x$design_grid))))
      stop("Retain one interval row per plan, including unavailable rows.", call. = FALSE)
    selected[i] <- ci$Score[1L]
    w <- if (kind == "Score") setNames(as.numeric(x$gstudy$design$scores == selected[i]), x$gstudy$design$scores) else
      if (is.matrix(x$weights)) x$weights[, selected[i]] else x$weights
    if (!isTRUE(all.equal(attr(ci, "weights"), w, tolerance = 0)))
      stop("Interval weights do not match the plotted score/composite.", call. = FALSE)
    idx <- which(series$Score == selected[i])
    ci <- ci[match(series$Scenario[idx], ci$Scenario), , drop = FALSE]
    if (!isTRUE(all.equal(ci$Estimate, series$Value[idx], tolerance = 0)) ||
        !identical(ci$PointStatus, series$Status[idx]))
      stop("Interval point estimates/statuses differ from the stored D-study; recreate the D-study and intervals.", call. = FALSE)
    available <- is.finite(ci$Lower) & is.finite(ci$Upper)
    bounded <- ci$MethodCode[available] != "delta" | !metrics %in% c("G", "Phi")
    if (!identical(ci$Available, available) ||
        any(ci$Lower[available] > ci$Upper[available] | (bounded & ci$Lower[available] < 0)) ||
        (metrics %in% c("G", "Phi") && any(bounded & ci$Upper[available] > 1)))
      stop("Inconsistent interval endpoints or availability; recreate intervals with confint(d).", call. = FALSE)
    series$Lower[idx] <- ci$Lower; series$Upper[idx] <- ci$Upper
    series$IntervalStatus[idx] <- ci$Status
    series$IntervalAvailable[idx] <- available
  }
  if (anyDuplicated(selected) || !setequal(selected, scores))
    stop("Supply exactly one interval result for every selected score/composite.", call. = FALSE)
  methods <- unique(unlist(lapply(intervals, function(ci) ci$MethodCode)))
  levels <- unique(unlist(lapply(intervals, function(ci) ci$Level)))
  if (length(methods) != 1L || length(levels) != 1L || !is.finite(levels) || levels <= 0 || levels >= 1)
    stop("All curve intervals must use the same method and confidence level.", call. = FALSE)
  # Split at unavailable rows before omitting them, so neither renderer bridges gaps.
  series$IntervalRun <- NA_integer_
  groups <- interaction(series$Group, if (is.null(series$Panel)) 1 else series$Panel, drop = TRUE)
  for (idx in split(seq_len(nrow(series)), groups)) {
    idx <- idx[order(series$X[idx], series$Scenario[idx])]
    starts <- c(TRUE, diff(as.integer(series$IntervalAvailable[idx])) != 0)
    run <- rep(idx[starts], diff(c(which(starts), length(idx)+1L)))
    series$IntervalRun[idx] <- ifelse(series$IntervalAvailable[idx], run, NA_integer_)
  }
  notes <- unique(c(sprintf("Nominal %s%% intervals; method: %s.",
      format(100*levels, trim = TRUE), paste(unique(unlist(lapply(intervals, function(ci) ci$Method))), collapse = "; ")),
    "Ribbons connect saved endpoints at requested counts; intermediate plans are not evaluated. Isolated intervals use vertical segments. Open points have no interval; gaps are not interpolated. Inspect unavailable_intervals and the saved interval objects for reasons and diagnostics.",
    "Ribbons do not establish simultaneous coverage across scores or metrics. Any joint coverage across plans depends on the interval method and its stated model; drawing a band adds no coverage guarantee. Interval overlap does not test a score or plan difference.",
    unlist(lapply(intervals, function(ci) attr(ci, "notes")), use.names = FALSE)))
  list(series = series, intervals = intervals, notes = notes)
}

.mfrm_mvds_unavailable_message <- function(status, compact = FALSE) {
  reasons <- unique(status)
  reasons[reasons == "Non-PSD component estimates"] <- "Covariance components failed validity checks."
  reasons[reasons == "Available"] <- "No finite estimates in the stored result."
  if (compact) return(paste(c("Estimates unavailable",
    strwrap(paste(reasons, collapse = "; "), width = 28)), collapse = "\n"))
  paste(c("Estimates unavailable", reasons,
    "Inspect the D-study table and component diagnostics."), collapse = "\n")
}

.mfrm_d_study_references <- function(reference, coefficients) {
  if (is.null(reference)) return(new_reference_lines())
  if (!is.numeric(reference) || is.complex(reference) || !is.null(dim(reference)) ||
      !length(reference) || any(!is.finite(reference)) || anyDuplicated(reference) ||
      any(reference < 0) || (coefficients && any(reference > 1))) {
    stop("`reference` must contain distinct finite numeric values ",
      if (coefficients) "between 0 and 1." else "greater than or equal to zero.", call. = FALSE)
  }
  new_reference_lines(axis = rep("y", length(reference)), value = reference,
    label = paste("User reference:", format(reference, trim = TRUE)),
    linetype = rep(c("dotted", "dashed"), length.out = length(reference)),
    role = rep("user_reference", length(reference)))
}

.mfrm_d_study_reference_note <- function(references) {
  if (!nrow(references)) return(NULL)
  paste0("User-chosen reference value(s): ", paste(format(references$value, trim = TRUE), collapse = ", "),
    ". Guides apply to every plotted metric; they are not universal acceptance rules, ",
    "confidence bounds or guarantees of achieved performance.")
}
