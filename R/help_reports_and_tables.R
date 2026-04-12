#' mfrmr Reports and Tables Map
#'
#' @description
#' Quick guide to choosing the right report or table helper in `mfrmr`.
#' Use this page when you know the reporting question but have not yet decided
#' which bundle, table, or reporting helper to call.
#'
#' @section Start with the question:
#' - "How should I document the model setup and run settings?"
#'   Use [specifications_report()].
#' - "Was data filtered, dropped, or mapped in unexpected ways?"
#'   Use [data_quality_report()] and [describe_mfrm_data()].
#' - "Did estimation converge cleanly and how formal is the precision layer?"
#'   Use [estimation_iteration_report()] and [precision_audit_report()].
#' - "Which facets are measurable, variable, or weakly separated?"
#'   Use [facet_statistics_report()], [measurable_summary_table()], and
#'   [facets_chisq_table()].
#' - "Are score categories functioning in a usable sequence?"
#'   Use [rating_scale_table()], [category_structure_report()], and
#'   [category_curves_report()].
#' - "Is the design linked well enough across subsets, forms, or waves?"
#'   Use [subset_connectivity_report()] and [plot_anchor_drift()].
#' - "What should go into the manuscript text and tables?"
#'   For `RSM` / `PCM`, use [reporting_checklist()], [build_apa_outputs()],
#'   and [build_summary_table_bundle()] or [export_summary_appendix()]. For
#'   bounded `GPCM`, stay on [reporting_checklist()], direct table/plot helpers,
#'   and summary-table appendix export; [build_apa_outputs()] and
#'   [export_mfrm_bundle()] remain out of scope.
#'
#' @section Recommended report route:
#' 1. Start with [specifications_report()] and [data_quality_report()] to
#'    document the run and confirm usable data.
#' 2. Continue with [estimation_iteration_report()] and
#'    [precision_audit_report()] to judge convergence and inferential strength.
#' 3. Use [facet_statistics_report()] and [subset_connectivity_report()] to
#'    describe spread, linkage, and measurability.
#' 4. Add [rating_scale_table()], [category_structure_report()], and
#'    [category_curves_report()] to document scale functioning.
#' 5. For `RSM` / `PCM`, finish with [reporting_checklist()] and
#'    [build_apa_outputs()] for manuscript-oriented output, then
#'    [build_summary_table_bundle()] for reusable handoff tables or
#'    [export_summary_appendix()] for direct appendix export. For bounded
#'    `GPCM`, skip [build_apa_outputs()] and [export_mfrm_bundle()]; use
#'    [reporting_checklist()], direct summaries/plots, and the summary-table
#'    appendix route only.
#'
#' @section Which output answers which question:
#' \describe{
#'   \item{[specifications_report()]}{Documents model type, estimation method,
#'   anchors, and core run settings. Best for method sections and audit trails.}
#'   \item{[data_quality_report()]}{Summarizes retained and dropped rows,
#'   missingness, and unknown elements. Best for data cleaning narratives.}
#'   \item{[estimation_iteration_report()]}{Shows replayed convergence
#'   trajectories. Best for diagnosing slow or unstable estimation.}
#'   \item{[precision_audit_report()]}{Summarizes whether `SE`, `CI`, and
#'   reliability indices are model-based, hybrid, or exploratory. Best for
#'   deciding how strongly to phrase inferential claims.}
#'   \item{[facet_statistics_report()]}{Bundles facet summaries, precision
#'   summaries, and variability tests. Best for facet-level reporting.}
#'   \item{[subset_connectivity_report()]}{Summarizes disconnected subsets and
#'   coverage bottlenecks. Best for linking and anchor strategy review.}
#'   \item{[rating_scale_table()]}{Gives category counts, average measures, and
#'   threshold diagnostics. Best for first-pass category evaluation.}
#'   \item{[category_structure_report()]}{Adds transition points and compact
#'   category warnings. Best for category-order interpretation.}
#'   \item{[category_curves_report()]}{Returns category-probability curve
#'   coordinates and summaries. Best for downstream graphics and report drafts.}
#'   \item{[reporting_checklist()]}{Turns analysis status into an action list
#'   with priorities and next steps. Best for closing reporting gaps.}
#'   \item{[build_apa_outputs()]}{Creates manuscript-draft text, notes,
#'   captions, and section maps from a shared reporting contract.}
#'   \item{[build_summary_table_bundle()]}{Converts supported `summary()`
#'   outputs into named `data.frame` tables with a compact index for appendix
#'   or manuscript handoff, and now supports bundle-level `summary()` /
#'   `plot()` for QC before export.}
#'   \item{[export_summary_appendix()]}{Exports those validated summary-table
#'   bundles as CSV and optional HTML appendix artifacts without requiring the
#'   broader fit-based export bundle.}
#'   \item{[apa_table()]}{Can now take those summary-table bundles directly,
#'   so a selected component can move from `summary()` to a formatted handoff
#'   table without rebuilding the analysis object path.}
#' }
#'
#' @section Practical interpretation rules:
#' - Use bundle summaries first, then drill down into component tables.
#' - Treat [precision_audit_report()] as the gatekeeper for formal inference.
#' - Treat category and bias outputs as complementary layers rather than
#'   substitutes for overall fit review.
#' - Treat zero-count score categories as scale-functioning caveats. Boundary
#'   zero-count categories can be retained with explicit `rating_min` /
#'   `rating_max`; intermediate zero-count categories require
#'   `keep_original = TRUE` and make adjacent thresholds weakly identified.
#'   `summary(describe_mfrm_data(...))` exposes these in `Notes`, printed
#'   `Caveats`, and `$caveats`; `summary(fit)` carries full structured caveats
#'   into printed `Caveats` and `$caveats`, with `Key warnings` as a short
#'   triage subset. Summary-table exports use `score_category_caveats` and
#'   `analysis_caveats`.
#' - Use [reporting_checklist()] before [build_apa_outputs()] when a report
#'   still needs missing diagnostics or clearer caveats.
#'
#' @section Typical workflow:
#' - Run documentation:
#'   [fit_mfrm()] -> [specifications_report()] -> [data_quality_report()].
#' - Precision and facet review:
#'   [diagnose_mfrm()] -> [precision_audit_report()] ->
#'   [facet_statistics_report()].
#' - Scale review:
#'   [rating_scale_table()] -> [category_structure_report()] ->
#'   [category_curves_report()].
#' - Manuscript handoff (`RSM` / `PCM`):
#'   [reporting_checklist()] -> [build_apa_outputs()] ->
#'   [build_summary_table_bundle()] -> `summary()` / `plot()` -> [apa_table()]
#'   or [export_summary_appendix()] /
#'   [export_mfrm_bundle()](include = "summary_tables").
#' - Bounded `GPCM` handoff:
#'   [reporting_checklist()] -> direct summaries/plots ->
#'   [build_summary_table_bundle()] -> [export_summary_appendix()].
#'
#' @section Companion guides:
#' - For visual follow-up, see [mfrmr_visual_diagnostics].
#' - For one-shot analysis routes, see [mfrmr_workflow_methods].
#' - For manuscript assembly, see [mfrmr_reporting_and_apa].
#' - For linking and DFF review, see [mfrmr_linking_and_dff].
#' - For legacy-compatible wrappers and exports, see [mfrmr_compatibility_layer].
#'
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' toy_small <- toy[toy$Person %in% unique(toy$Person)[1:12], , drop = FALSE]
#' fit <- fit_mfrm(
#'   toy_small,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   maxit = 200
#' )
#' diag <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "both")
#'
#' spec <- specifications_report(fit)
#' summary(spec)$overview
#'
#' prec <- precision_audit_report(fit, diagnostics = diag)
#' summary(prec)$checks
#'
#' checklist <- reporting_checklist(fit, diagnostics = diag)
#' subset(checklist$checklist, Section == "Visual Displays", c("Item", "NextAction"))
#'
#' apa <- build_apa_outputs(fit, diagnostics = diag)
#' apa$section_map[, c("Heading", "Available")]
#' bundle <- build_summary_table_bundle(checklist)
#' bundle$table_index
#' }
#'
#' @name mfrmr_reports_and_tables
NULL
