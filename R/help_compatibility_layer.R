#' mfrmr Compatibility Layer Map
#'
#' @description
#' Guide to the legacy-compatible wrappers and text/file exports in `mfrmr`.
#' Use this page when you need continuity with older compatibility-oriented workflows,
#' fixed-width reports, or graph/score file style outputs.
#'
#' @section When to use this layer:
#' - You are reproducing an older workflow that expects one-shot wrappers.
#' - You need fixed-width text blocks for console, logs, or archival handoff.
#' - You need graphfile or scorefile style outputs for downstream legacy tools.
#' - You are checking column/metric parity against a compatibility contract.
#'
#' @section When not to use this layer:
#' - For standard estimation, use [fit_mfrm()] plus [diagnose_mfrm()].
#' - For report bundles, use [mfrmr_reports_and_tables].
#' - For manuscript text, use [build_apa_outputs()] and [reporting_checklist()].
#' - For visual follow-up, use [mfrmr_visual_diagnostics].
#'
#' @section Compatibility map:
#' \describe{
#'   \item{[run_mfrm_facets()]}{One-shot legacy-compatible wrapper that fits,
#'   diagnoses, and returns key tables in one object.}
#'   \item{[mfrmRFacets()]}{Alias for [run_mfrm_facets()] kept for continuity.}
#'   \item{[build_fixed_reports()]}{Fixed-width interaction and pairwise text
#'   blocks. Best when a text-only compatibility artifact is required.}
#'   \item{[facets_output_file_bundle()]}{Graphfile/scorefile style CSV and
#'   fixed-width exports for legacy pipelines.}
#'   \item{[facets_parity_report()]}{Column and metric contract audit against
#'   the compatibility specification. Use only when explicit parity checking is
#'   part of the task.}
#' }
#'
#' @section Preferred replacements:
#' - Instead of [run_mfrm_facets()], prefer:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [reporting_checklist()].
#' - Instead of [build_fixed_reports()], prefer:
#'   [bias_interaction_report()] -> [build_apa_outputs()].
#' - Instead of [facets_output_file_bundle()], prefer:
#'   [category_curves_report()] or [category_structure_report()] plus
#'   [export_mfrm_bundle()].
#' - Instead of [facets_parity_report()] for routine QA, prefer:
#'   [reference_case_audit()] for package-native completeness auditing or
#'   [reference_case_benchmark()] for internal benchmark cases.
#'
#' @section Practical migration rules:
#' - Keep compatibility wrappers only where a downstream consumer truly needs
#'   the old layout or fixed-width format.
#' - For new scripts, start from package-native bundles and add compatibility
#'   outputs only at the export boundary.
#' - Treat compatibility outputs as presentation contracts, not as the primary
#'   analysis objects.
#'
#' @section Typical workflow:
#' - Legacy handoff:
#'   [run_mfrm_facets()] -> [build_fixed_reports()] ->
#'   [facets_output_file_bundle()].
#' - Mixed workflow:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [build_apa_outputs()] ->
#'   compatibility export only if required.
#' - Contract audit:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [facets_parity_report()].
#'
#' @section Companion guides:
#' - For standard reports/tables, see [mfrmr_reports_and_tables].
#' - For manuscript-draft reporting, see [mfrmr_reporting_and_apa].
#' - For visual diagnostics, see [mfrmr_visual_diagnostics].
#' - For linking and DFF workflows, see [mfrmr_linking_and_dff].
#' - For end-to-end routes, see [mfrmr_workflow_methods].
#'
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' toy_small <- toy[toy$Person %in% unique(toy$Person)[1:12], , drop = FALSE]
#'
#' run <- run_mfrm_facets(
#'   data = toy_small,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   maxit = 10
#' )
#' summary(run)
#'
#' fixed <- build_fixed_reports(
#'   estimate_bias(
#'     run$fit,
#'     run$diagnostics,
#'     facet_a = "Rater",
#'     facet_b = "Criterion",
#'     max_iter = 1
#'   ),
#'   branch = "original"
#' )
#' names(fixed)
#' }
#'
#' @name mfrmr_compatibility_layer
NULL
