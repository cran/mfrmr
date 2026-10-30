# Public plot methods reuse the existing numerical and drawing implementations.
.mfrm_plot_call <- function(fun, x, args = list(), additional_allowed = character()) {
  nms <- names(args)
  if (length(args) && (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms)))
    stop("Name each plot argument once.", call. = FALSE)
  allowed <- setdiff(names(formals(fun)), c(names(formals(fun))[1L], "..."))
  allowed <- union(allowed, additional_allowed)
  if ("..." %in% names(formals(fun))) allowed <- c(allowed,
    "col", "pch", "lty", "lwd", "cex", "xlim", "ylim", "xlab", "ylab",
    "main", "sub", "cex.axis", "cex.lab", "cex.main", "las", "bty", "axes", "xaxt", "yaxt")
  unused <- setdiff(nms, allowed)
  if (length(unused)) stop("Arguments not supported for this plot: ",
    paste(unused, collapse = ", "), call. = FALSE)
  do.call(fun, c(list(x), args))
}

.mfrm_plot_selector <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
    stop("`", name, "` must be a single nonempty string.", call. = FALSE)
  x
}

.mfrm_plot_rename <- function(args, from, to) {
  if (!from %in% names(args)) return(args)
  if (to %in% names(args)) stop("Supply only one of `", from, "` and `", to, "`.", call. = FALSE)
  names(args)[names(args) == from] <- to
  args
}

.mfrm_fit_plot_route <- function(type, view = NULL) {
  routes <- c(apa_summary = "plot_apa_figure_one", bubble = "plot_bubble",
    displacement = "plot_displacement", facets_chisq = "plot_facets_chisq",
    fair_average = "plot_fair_average", guttman = "plot_guttman_scalogram",
    local_dependence = "plot_local_dependence_heatmap", marginal_fit = "plot_marginal_fit",
    marginal_pairwise = "plot_marginal_pairwise", person_fit = "plot_person_fit",
    qc = "plot_qc_dashboard", rater_severity = "plot_rater_severity_profile",
    reliability = "plot_reliability_snapshot", residual_matrix = "plot_residual_matrix",
    residual_qq = "plot_residual_qq", thresholds = "plot_threshold_ladder",
    unexpected = "plot_unexpected", facet_equivalence = "plot_facet_equivalence",
    facet_dashboard = "plot_facet_quality_dashboard", bias = "plot_bias_interaction",
    residual_pca = "plot_residual_pca")
  if (type == "agreement") return(if (identical(view, "heatmap"))
    "plot_rater_agreement_heatmap" else "plot_interrater_agreement")
  if (type == "wright" && identical(view, "unified")) return("plot_wright_unified")
  if (type == "shrinkage" && identical(view, "funnel")) return("plot_shrinkage_funnel")
  unname(routes[type])
}

.mfrm_dispatch_fit_plot <- function(x, type, args) {
  view <- args$view
  if (!is.null(view)) .mfrm_plot_selector(view, "view")
  if (!is.null(args$comparison)) {
    comparison <- args$comparison
    args$comparison <- NULL
    if (type == "trajectory") {
      if (inherits(comparison, "mfrm_fit")) comparison <- list(comparison)
      if (!is.list(comparison) || !length(comparison) ||
          !all(vapply(comparison, inherits, logical(1), "mfrm_fit")))
        stop("Trajectory `comparison` must contain one or more fitted models in wave order.", call. = FALSE)
      fits <- c(list(x), comparison)
      labels <- args$labels
      if (is.null(labels)) labels <- c("Reference", names(comparison))
      if (!is.character(labels) || length(labels) != length(fits) || anyNA(labels) ||
          any(!nzchar(labels)) || anyDuplicated(labels))
        stop("Supply distinct `labels` for the reference and every comparison wave, in order.", call. = FALSE)
      names(fits) <- labels
      args$labels <- NULL
      return(.mfrm_plot_call(plot_rater_trajectory, fits, args))
    }
    if (!type %in% c("wright", "ccc"))
      stop("`comparison` is supported only for Wright, CCC or trajectory plots.", call. = FALSE)
    args$comparison <- comparison
    args$type <- type
    return(.mfrm_plot_call(plot_compare_mfrm, x, args))
  }
  if (type == "trajectory") stop("Trajectory plots require `comparison` and ordered wave labels.", call. = FALSE)
  name <- .mfrm_fit_plot_route(type, view)
  if (is.na(name)) stop("No consolidated route for this plot type.", call. = FALSE)
  if (type == "agreement") {
    if (!is.null(view) && !view %in% c("pairs", "heatmap"))
      stop("Agreement `view` must be 'pairs' or 'heatmap'.", call. = FALSE)
    args$view <- NULL
    if (!identical(view, "heatmap")) {
      args <- .mfrm_plot_rename(args, "metric", "plot_type")
      if (identical(args$plot_type, "correlation")) args$plot_type <- "corr"
    }
  } else if (type %in% c("wright", "shrinkage")) {
    args$view <- NULL
  } else if (name %in% c("plot_displacement", "plot_facets_chisq", "plot_fair_average",
                         "plot_marginal_fit", "plot_unexpected", "plot_facet_quality_dashboard",
                         "plot_residual_pca")) {
    args <- .mfrm_plot_rename(args, "view", "plot_type")
  } else if (type == "bias") {
    args <- .mfrm_plot_rename(args, "view", "plot")
  } else if (type == "facet_equivalence") {
    args <- .mfrm_plot_rename(args, "view", "type")
  }
  additional_allowed <- if (identical(name, "plot_fair_average")) {
    setdiff(names(formals(fair_average_table)), c("fit", "diagnostics", "..."))
  } else character()
  .mfrm_plot_call(get(name, envir = environment(.mfrm_dispatch_fit_plot)), x, args,
    additional_allowed = additional_allowed)
}

#' @rdname compute_information
#' @details `plot()` method: Plot design-weighted information or conditional precision
#'
#' Plot a saved [compute_information()] result to see where the observed
#' assessment design provides more or less precision across ability.
#' @inheritParams plot_information
#' @param ... Additional base-graphics arguments for the selected view.
#' @inherit plot_information return
#' @inheritSection plot_information Plot types
#' @inheritSection plot_information Which type should I use?
#' @inheritSection plot_information Interpreting output
#' @inheritSection plot_information Returned data when draw = FALSE
#' @inheritSection plot_information Typical workflow
#' @seealso [compute_information()], [mfrmr_plot_migration]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_information", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_operational")
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "MML", model = "RSM",
#'                  quad_points = 7, maxit = 30)
#' info <- compute_information(fit)
#' tif_data <- plot(info, type = "tif", draw = FALSE)
#' head(tif_data$data$plot)
#' iif_data <- plot(info, type = "iif", facet = "Rater", draw = FALSE)
#' head(iif_data$data$plot)
#' }
#' @export
#' @return For \code{plot()}:
#' Invisibly, an \code{mfrm_plot_data} object.
plot.mfrm_information <- function(x, type = c("tif", "iif", "se", "sem", "csem", "both"),
                                  facet = NULL, draw = TRUE, ...) {
  .mfrm_plot_call(plot_information, x, c(list(type = match.arg(type), facet = facet, draw = draw), list(...)))
}

#' Plot differential facet functioning
#'
#' Compare saved group contrasts as ordered bars or a heatmap. Read the
#' analysis method and scale before interpreting the size of a difference.
#' @param x A result from [analyze_dff()] or [analyze_dif()].
#' @param type `"summary"` for ordered effects or `"heatmap"` for the interaction grid.
#' @param ... Named controls described below. Both views accept `draw = FALSE`
#'   to return plot data without drawing. Controls for the other view are rejected.
#' @section Summary controls:
#' - `top_n = 30` limits displayed contrasts; `sort_by` is `"abs_effect"`
#'   (default), `"effect"`, or `"classification"`.
#' - `ci_level = NULL` omits intervals. A level such as `0.95` requests
#'   approximate normal intervals from available linked-refit standard errors;
#'   residual comparisons require `NULL`.
#' - `effect_thresholds` adds symmetric effect-size guide lines, not
#'   classification boundaries. `effect_axis_label` overrides the axis label.
#' - `preset` is `"standard"`, `"publication"`, `"compact"`, or `"monochrome"`.
#' @section Heatmap controls:
#' - `metric = "obs_exp"` shows cell mean observed-minus-expected scores;
#'   `"t"` shows a descriptive scaled residual, not a t-test. These require
#'   the cell table from `method = "residual"`. `"contrast"` shows the saved
#'   pairwise group contrasts, on the analysis method's scale.
#' @section Heatmap appearance:
#' - `show_values = TRUE` labels cells; `value_digits = 2` controls rounding.
#' - `flag_threshold = NULL` adds no flags. A non-negative value outlines
#'   cells at or above that absolute value; `flag_color = "black"` sets the border.
#' - `scale_limit = NULL` uses the observed range. A positive value fixes
#'   symmetric color limits for comparing plots. It does not change saved values.
#' @inheritSection plot_dif_summary Interpreting output
#' @section Read the figure with its source table:
#' Summary labels follow `Level | Group1 | Group2`. Both summary bars and
#' contrast heatmaps use **Group1 minus Group2**. A positive residual contrast
#' means higher mean residual scores in Group1; a positive refit contrast means
#' a larger facet location in Group1. With the usual difficulty/severity
#' orientation, the latter means harder for Group1. Check the model's facet
#' signs and rubric coding before interpreting either as favorable treatment.
#'
#' `plot_data(p, "pairs")` retains the source rows, counts and analysis conditions.
#' For summaries, `PlotStatus` distinguishes displayed rows, finite contrasts
#' omitted by `top_n`, and unavailable contrasts. The last two do not mean zero
#' or no differential functioning. If all contrasts are unavailable, the
#' summary cannot draw points; inspect `x$dif_table` or a contrast heatmap.
#' Heatmap cells with missing values remain missing. An observed-minus-expected
#' cell mean can still be shown when its cell is below `min_obs`; inspect
#' `sparse`. The scaled-residual and pairwise-contrast routes withhold such values.
#'
#' `N_Group1` and `N_Group2` count retained response rows. Residual results also
#' give distinct-Person counts, which need not be independent observations.
#' `plot_data(p, "cases")` preserves residual response matches and exclusion
#' reasons; it is `NULL` for refit results, which have no such response audit.
#' Refit source rows instead report `ScaleLinkStatus`, `ContrastComparable`,
#' fixed-target flags and subgroup linking details. A finite point or an
#' adequate linking screen does not establish formal DFF inference. See
#' [mfrmr_linking_and_dff] for the different linking screens.
#' @return An `mfrm_plot_data` object. The summary's `data$data` holds the
#'   displayed contrasts; the heatmap's `data$matrix` holds the plotted matrix.
#'   Both retain source rows in `data$pairs`, response matches (when available)
#'   in `data$cases`, and source settings in `data$analysis_config`.
#'   When supplied to [analyze_dff()], the separate prospective interaction
#'   review is retained in `data$comparison_review`. Its states do not classify
#'   the displayed residual/refit effects; inspect its own model and weights.
#'   Retain the returned notes and
#'   unavailable values when exporting; a plot does not qualify new inference.
#' @seealso [analyze_dff()], [dif_interaction_table()], [mfrmr_plot_migration]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_dif", run.dontrun = TRUE)
#' \dontrun{
#' # Saved synthetic JML/RSM analysis of example_bias.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "bias-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' diag <- example_fit$diagnostics
#' # Group labels from the same synthetic data used for the saved fit.
#' ratings <- load_mfrmr_data(example_fit$dataset)
#' dff <- analyze_dff(fit, diagnostics = diag,
#'                    facet = "Rater", group = "Group", data = ratings)
#' unique(dff$dif_table$ClassificationSystem)
#' p <- plot(dff, type = "summary", draw = FALSE)
#' comparisons <- plot_data(p, "pairs")
#' comparisons[, c("Level", "Group1", "Group2", "Contrast",
#'                  "N_Group1", "N_Group2", "PlotStatus")]
#' table(plot_data(p, "cases")$Reason)
#' plot_data(p, "notes")
#' }
#' @export
plot.mfrm_dif <- function(x, type = c("summary", "heatmap"), ...) {
  type <- match.arg(type)
  .mfrm_plot_call(if (type == "summary") plot_dif_summary else plot_dif_heatmap, x, list(...))
}

#' @rdname dif_interaction_table
#' @details `plot()` method: Plot a differential-functioning interaction grid
#'
#' Display saved residual means by facet level and group. This is a descriptive
#' screen for follow-up; residual differences do not isolate differential functioning.
#' @param x A result from [dif_interaction_table()].
#' @param type The sole view, `"heatmap"`.
#' @param ... Heatmap controls below, plus `draw = FALSE` to return data
#'   without drawing. This input accepts `metric = "obs_exp"` or `"t"`;
#'   `"contrast"` requires a result from [analyze_dff()].
#' @inheritSection plot.mfrm_dif Heatmap appearance
#' @section Interpretation:
#' For `plot()`:
#'
#' Positive `"obs_exp"` values mean observed scores exceed model expectations;
#' negative values mean they fall below expectations. These values are in
#' score units, not logits. The `"t"` metric is a scaled residual without a
#' calibrated test or p-value. Color and optional threshold outlines are
#' descriptive; an unavailable cell is not a zero residual.
#' @return For `plot()`: An `mfrm_plot_data` object with `data$matrix`, `data$pairs`, the
#'   selected metric and retained interpretation notes. Use
#'   `plot_data(p, "matrix")` to extract the matrix from a saved plot `p`.
#' @seealso [dif_interaction_table()], [plot.mfrm_dif()], [mfrmr_plot_migration]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_dif_interaction", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_bias")
#'
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", model = "RSM", maxit = 300)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' int <- dif_interaction_table(fit, diag, facet = "Rater",
#'                              group = "Group", data = ratings, min_obs = 2)
#' heat <- plot(int, type = "heatmap", metric = "obs_exp", draw = FALSE)
#' dim(heat$data$matrix)
#' # Obs-Exp is in score units, not logits. A scaled residual
#' # (metric = "t") is descriptive, not a calibrated t-test.
#' }
#' @export
plot.mfrm_dif_interaction <- function(x, type = "heatmap", ...) {
  if (!identical(type, "heatmap")) stop("Interaction results support type = 'heatmap'.", call. = FALSE)
  .mfrm_plot_call(plot_dif_heatmap, x, list(...))
}

#' @rdname detect_anchor_drift
#' @details `plot()` method: Plot a saved anchor-drift review
#'
#' Inspect relative changes after the common-element offset used by
#' [detect_anchor_drift()]. Read `summary(x)` and link support alongside the plot.
#' @param type `"drift"` for relative changes with screening flags,
#'   `"heatmap"` for wave-by-element patterns, or `"forest"` for per-wave
#'   estimates with approximate normal whiskers. A linking chain is a different
#'   result: use `plot(chain, type = "chain")` after [build_equating_chain()].
#' @section Interpretation:
#' For `plot()`:
#'
#' Drift and heatmap views show the saved relative drift. Flags use the
#' screening rules chosen during analysis; changing the display does not
#' recalibrate these rules. Subtracting the common-element offset can remove
#' a panel-wide shift, so unflagged elements do not establish absence of annual drift.
#'
#' Forest whiskers describe individual wave estimates, not uncertainty in
#' their differences. They do not include linking-offset uncertainty or
#' cross-wave covariance, and their overlap is not a test of change. Missing
#' standard errors produce no whisker, not a zero-width interval.
#' @return For `plot()`: An `mfrm_plot_data` object, or `NULL` when no drift rows remain.
#'   Drift/heatmap data use `data$table` (and `data$matrix` for a heatmap);
#'   forest rows use `data$data`, including `CI_Lower` and `CI_Upper`.
#' @seealso [detect_anchor_drift()], [mfrmr_plot_migration]
#' @examples # Long-running workflow; run explicitly when needed.
#' # example("plot.mfrm_anchor_drift", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_core")
#' people <- unique(ratings$Person)
#' d1 <- ratings[ratings$Person %in% people[1:12], , drop = FALSE]
#' d2 <- ratings[ratings$Person %in% people[13:24], , drop = FALSE]
#' fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML")
#' fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML")
#' # Review convergence and common-element support before interpreting changes.
#' summary(fit1)
#' summary(fit2)
#' drift <- detect_anchor_drift(list(W1 = fit1, W2 = fit2))
#' drift_plot <- plot(drift, type = "drift", draw = FALSE)
#' class(drift_plot)
#' names(drift_plot$data)
#' chain <- build_equating_chain(list(F1 = fit1, F2 = fit2))
#' chain_plot <- plot(chain, type = "chain", draw = FALSE)
#' head(chain_plot$data$table)
#' if (interactive()) {
#'   plot(drift, type = "heatmap", preset = "publication")
#' }
#' }
#' @export
plot.mfrm_anchor_drift <- function(x, type = c("drift", "heatmap", "forest"), ...) {
  .mfrm_plot_call(plot_anchor_drift, x, c(list(type = match.arg(type)), list(...)))
}

#' @rdname run_qc_pipeline
#' @details `plot()` method: Plot a saved QC pipeline
#'
#' Display the checks already returned by [run_qc_pipeline()]. Use the colored
#' verdicts for triage, then read the returned detail for the reason and next review.
#' @inheritParams plot_qc_pipeline
#' @param ... Additional base-graphics arguments for the selected view.
#' @inherit plot_qc_pipeline details return
#' @inheritSection plot_qc_pipeline Interpretation
#' @seealso [run_qc_pipeline()], [mfrmr_plot_migration]
#' @examples # Long-running workflow; run explicitly when needed.
#' # example("plot.mfrm_qc_pipeline", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("study1")
#' # These observed data trigger a category-support warning: the QC plot
#' # displays that review requirement; convergence does not remove it.
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 300)
#' qc <- run_qc_pipeline(fit)
#' verdicts <- plot(qc, draw = FALSE)
#' head(verdicts)
#' }
#' @export
#' @details For \code{plot()}:
#' Two plot types are provided for visual triage of QC results:
#' \itemize{
#' \item \strong{\code{"traffic_light"}} (default): A horizontal bar chart with one row
#' per QC check.  Bars are coloured green (Pass), amber (Warn), or red
#' (Fail).  Provides an at-a-glance summary of the current QC review state.
#' \item \strong{\code{"detail"}}: A panel showing the saved \code{Detail} text beside each
#' verdict. Read this explanation and the returned table to understand
#' why a check was flagged or unavailable.
#' }
#' @return For \code{plot()}:
#' Invisible verdicts tibble from the QC pipeline.
plot.mfrm_qc_pipeline <- function(x, type = c("traffic_light", "detail"), draw = TRUE, ...) {
  .mfrm_plot_call(plot_qc_pipeline, x, c(list(type = match.arg(type), draw = draw), list(...)))
}

#' @rdname estimate_bias
#' @details `plot()` method: Plot saved bias-interaction results
#'
#' Inspect the pattern of conditional interaction estimates from [estimate_bias()].
#' Start with ranked interactions, then use other views to understand the pattern.
#' @param x A result from [estimate_bias()].
#' @param type `"scatter"`, `"ranked"`, `"heatmap"`, `"abs_t_hist"`, or `"facet_profile"`.
#' @param ... For `plot()`: Named controls below. `draw = FALSE` returns the plot data
#'   without drawing. Controls do not change the fitted calibration.
#'
#' For `summary()`: Reserved for generic compatibility.
#' @section Controls:
#' For `plot()`:
#'
#' - `facet_a`, `facet_b`, or `interaction_facets` identify the interaction;
#'   `top_n = 40` and `sort_by = "abs_t"` (`"abs_bias"` or `"prob"` also
#'   accepted) select ranked interactions.
#' - `abs_t_warn = 2`, `abs_bias_warn = 0.5`, and `p_max = 0.05` set
#'   screening rules when preparing a report from the bias result. They are
#'   review aids, not guarantees of a calibrated error rate.
#' - `show_ci = FALSE` omits whiskers; `ci_level = 0.95` chooses their level
#'   for scatter/ranked views. GPCM rows retain conditional profile bounds,
#'   with dashed whiskers for search-limited or unverified bounds. Other rows
#'   use conditional per-cell SEs. Missing profile bounds are not replaced
#'   with another method. These intervals omit calibration uncertainty;
#'   endpoint completion is not evidence of coverage.
#' - `preset` selects `"standard"`, `"publication"`, `"compact"`, or
#'   `"monochrome"`; `palette` overrides named colors (`normal`, `flag`,
#'   `hist`, `profile`); `label_angle` controls label rotation. `title` sets
#'   title metadata (`main` is a compatibility alias; supply only one).
#' @inherit plot_bias_interaction details
#' @inheritSection plot_bias_interaction Plot types
#' @inheritSection plot_bias_interaction Interpreting output
#' @return For `plot()`: An `mfrm_plot_data` object with view-specific tables, retained
#'   screening rules and interpretation notes. Plot agreement across views
#'   is not independent evidence of bias.
#'   Use `print(p)` or `plot_data(p, "notes")` for the interaction definition,
#'   score-direction caveat, magnitude-only facet-profile interpretation and
#'   conditional uncertainty qualifications outside the figure.
#' @seealso [estimate_bias()], [bias_interaction_report()], [mfrmr_plot_migration],
#'   [mfrmr_visual_diagnostics]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_bias", run.dontrun = TRUE)
#' \dontrun{
#' # Saved synthetic JML/RSM analysis of example_bias.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "bias-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' diagnostics <- example_fit$diagnostics
#' bias <- estimate_bias(fit, diagnostics = diagnostics,
#'   facet_a = "Rater", facet_b = "Criterion")
#' p <- plot(bias, type = "ranked", preset = "publication", draw = FALSE)
#' head(plot_data(p, "ranked_table"))
#' }
#' @export
#' @details For \code{plot()}:
#' Visualization front-end for \code{\link[=bias_interaction_report]{bias_interaction_report()}} with multiple views.
#' With \code{draw = FALSE}, the returned plot data include \code{plot_long},
#' \code{plot_annotations}, \code{flag_summary}, and \code{plot_settings} in addition to the
#' view-specific \code{ranked_table}, \code{scatter_data}, \code{facet_profile}, and heatmap
#' components. Use these fields when rebuilding the same screening view in
#' ggplot2, plotly, Quarto, or a dashboard.
#' Saved report thresholds determine flags, reference lines and returned notes;
#' recreate \code{\link[=bias_interaction_report]{bias_interaction_report()}} from the saved bias object to change
#' them. Scatter and ranked views support \code{\link[=as_ggplot]{as_ggplot()}}; the other views require
#' their base renderer or custom graphics. Conditional bounds and their status
#' remain separate from calibration-inclusive uncertainty.
plot.mfrm_bias <- function(x, type = c("scatter", "ranked", "heatmap", "abs_t_hist", "facet_profile"), ...) {
  args <- list(...)
  if (identical(type, "bias")) {
    args <- .mfrm_plot_rename(args, "view", "plot")
  } else args <- c(list(plot = match.arg(type)), args)
  .mfrm_plot_call(plot_bias_interaction, x, args)
}

#' @rdname analyze_residual_pca
#' @details `plot()` method: Plot saved residual PCA
#'
#' Display a saved residual PCA or its saved parallel-analysis comparison.
#' Plotting does not generate a new parallel reference.
#' @param x A result from [analyze_residual_pca()].
#' @param type `"scree"`, `"parallel_scree"`, `"parallel_excess"`, or `"loadings"`.
#'   Existing `"overall_*"` and `"facet_*"` spellings remain supported.
#' @param ... `facet` selects a facet for `mode = "facet"` (default: the first
#'   available facet); `component = 1` selects a loading component; `top_n = 20`
#'   limits displayed loadings. `preset` is `"standard"`, `"publication"`,
#'   `"compact"`, or `"monochrome"`; `draw = FALSE` returns data without drawing.
#' @section Interpretation:
#' For `plot()`:
#'
#' - `"scree"` shows component eigenvalues. A dominant early component is
#'   exploratory evidence of residual structure, not a standalone test of
#'   unidimensionality.
#' - `"parallel_scree"` shows the observed eigenvalues and saved reference
#'   mean/cutoff; `"parallel_excess"` shows the observed-minus-cutoff difference.
#'   Compute these references with `parallel = TRUE` in [analyze_residual_pca()]
#'   first. Exceeding a cutoff suggests follow-up, not proof of multidimensionality.
#'   Inspect `parallel_status` and `bootstrap_trials` for unavailable comparisons
#'   and distinguish residual permutation from refitted model bootstrap.
#' - `"loadings"` shows the largest absolute loadings; inspect both signs
#'   and magnitudes to identify the contributing variables or elements.
#'
#' Use overall mode before locating patterns within individual facets. Repeated
#' plots of a saved result reuse the same PCA and do not provide independent evidence.
#' @inherit plot_residual_pca return
#' @seealso [analyze_residual_pca()], [mfrmr_plot_migration]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_residual_pca", run.dontrun = TRUE)
#' \dontrun{
#' ratings_all <- load_mfrmr_data("example_core")
#' selected_people <- unique(ratings_all$Person)[1:24]
#' ratings <- ratings_all[
#'   ratings_all$Person %in% selected_people, , drop = FALSE]
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML")
#' diag <- diagnose_mfrm(fit, residual_pca = "overall")
#' pca <- analyze_residual_pca(diag, mode = "overall")
#' plt <- plot(pca, type = "scree", draw = FALSE)
#' head(plt$data)
#' # Ten replicates illustrate the display, not a stable reference cutoff.
#' pca_pa <- analyze_residual_pca(diag, mode = "overall", parallel = TRUE,
#'                                parallel_reps = 10, seed = 2026)
#' pa <- plot(pca_pa, type = "parallel_scree", draw = FALSE)
#' head(pa$data)
#' plt_load <- plot(pca, type = "loadings", component = 1, draw = FALSE)
#' head(plt_load$data)
#' if (interactive()) {
#'   plot(pca, type = "scree", preset = "publication")
#' }
#' }
#' @export
#' @return For \code{plot()}:
#' A named list of plotting data (class \code{mfrm_plot_data}) with:
#' \itemize{
#' \item \code{plot}: \code{"scree"}, \code{"parallel_scree"}, \code{"parallel_excess"}, or
#' \code{"loadings"}
#' \item \code{mode}: \code{"overall"} or \code{"facet"}
#' \item \code{facet}: facet name (or \code{NULL})
#' \item \code{title}: plot title text
#' \item \code{data}: underlying table used for plotting
#' \item \code{InferenceTier}, \code{SupportsFormalInference},
#' \code{PrimaryReportingEligible}, \code{ReportingUse}, and \code{DecisionUse}:
#' machine-readable exploratory-screening guards
#' }
plot.mfrm_residual_pca <- function(x, type = "scree", mode = "overall", ...) {
  type <- .mfrm_plot_selector(type, "type")
  args <- list(...)
  if (type == "residual_pca") {
    args <- .mfrm_plot_rename(args, "view", "plot_type")
    return(.mfrm_plot_call(plot_residual_pca, x, c(list(mode = mode), args)))
  }
  if (grepl("^(overall|facet)_", type)) {
    if (!missing(mode)) stop("Choose either a prefixed PCA type or `mode`.", call. = FALSE)
    mode <- sub("_.*$", "", type)
    type <- sub("^(overall|facet)_", "", type)
  }
  .mfrm_plot_call(plot_residual_pca, x, c(list(plot_type = type, mode = mode), args))
}

#' @rdname diagnose_mfrm
#' @details `plot()` method: Plot saved ordinary-model diagnostics
#'
#' Inspect an existing diagnostic result to locate cells or elements needing
#' follow-up. The selected display retains unavailable results and screening limits.
#' @param x An `mfrm_diagnostics` result from [diagnose_mfrm()].
#' @param type `"marginal_fit"`, `"marginal_pairwise"`, `"bubble"`, or `"residual_pca"`.
#' @param view For marginal fit, `"std_residual"` (default) or `"prop_diff"`;
#'   for bubbles, `"measure"` (default) or `"infit_outfit"`; for PCA,
#'   `"scree"`, `"parallel_scree"`, `"parallel_excess"`, or `"loadings"`.
#'   Omit it for `"marginal_pairwise"` and use `metric` instead.
#' @param ... All views accept `draw = FALSE`, `top_n`, and `preset`.
#'   Marginal views accept `facet`; pairwise plots additionally accept
#'   `metric = "exact"` or `"adjacent"`. Bubbles accept `fit_stat` (`"Infit"`
#'   or `"Outfit"`), `bubble_size` (`"SE"`, `"N"`, or `"equal"`), `facets`,
#'   `include_person = FALSE`, and the heuristic `fit_range = c(0.5, 1.5)`.
#'   PCA accepts `mode`, `facet`, and `component`, as described in
#'   [plot.mfrm_residual_pca()]. Controls for other views are rejected.
#' @section Choosing a view:
#' For `plot()`:
#'
#' Marginal fit and pairwise plots require strict marginal diagnostics from
#' [diagnose_mfrm()] with `diagnostic_mode = "both"` or `"marginal_fit"` for
#' supported MML RSM/PCM models. Positive residuals indicate more observed
#' category usage or agreement than expected; negative values indicate less.
#' Returned `full_table` and `retention` distinguish omitted and unavailable
#' rows. A flag is a descriptive review rule, not a calibrated test.
#'
#' Bubble plots place measures against mean-square fit, or Infit against Outfit.
#' Bubble size can represent SE or observation count; neither a review band
#' nor bubble size alone justifies deleting an element.
#'
#' The PCA convenience route can derive PCA from the diagnostic result.
#' For repeated PCA displays, first save [analyze_residual_pca()] and use its
#' `plot()` method. Parallel-reference generation remains an analysis step.
#' @return For `plot()`: An `mfrm_plot_data` object. Use [plot_data()] to inspect its tables
#'   and notes. The plot does not add unavailable diagnostics or qualify intervals.
#' @seealso [diagnose_mfrm()], [plot.mfrm_residual_pca()], [mfrmr_plot_migration]
#' @export
plot.mfrm_diagnostics <- function(x, type = c("marginal_fit", "marginal_pairwise", "bubble", "residual_pca"),
                                  view = NULL, ...) {
  args <- list(...)
  if (!is.null(view)) args$view <- view
  .mfrm_dispatch_fit_plot(x, match.arg(type), args)
}
