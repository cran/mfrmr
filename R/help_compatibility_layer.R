#' mfrmr Compatibility Layer Map
#'
#' @description
#' Guide to the legacy-compatible wrappers and text/file exports in `mfrmr`.
#' Use this page when you need continuity with older compatibility-oriented workflows,
#' fixed-width reports, or graph/score file style outputs.
#'
#' This compatibility layer currently applies mainly to diagnostics-based
#' `RSM` / `PCM` workflows. First-release `GPCM` fits now also support
#' graph-only compatibility-style exports, while scorefile and
#' diagnostics-driven compatibility outputs remain limited to `RSM` / `PCM`.
#' Treat this layer as a presentation/contract surface, not as a claim of
#' FACETS or ConQuest numerical equivalence.
#'
#' SPSS is treated differently from FACETS and ConQuest: `mfrmr` currently
#' supports table/data-frame/CSV handoff for SPSS-oriented reporting workflows,
#' but it does not generate SPSS syntax, write native SPSS system files, execute
#' SPSS estimators, or claim SPSS numerical parity.
#'
#' @section When to use this layer:
#' - You are reproducing an older workflow that expects one-shot wrappers.
#' - You need fixed-width text blocks for console, logs, or archival handoff.
#' - You need graphfile or scorefile style outputs for downstream legacy tools.
#' - You are checking column coverage and metric consistency against a
#'   compatibility contract.
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
#'   the compatibility specification. Use only when an explicit compatibility
#'   contract audit is part of the task; the function name is historical and does
#'   not by itself imply external FACETS equivalence.}
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
#' - Use [compatibility_alias_table()] when you need to check which aliases are
#'   still retained and which package-native names should be used in new code.
#' - Use `reporting_checklist(fit)$software_scope` to review the current
#'   FACETS, ConQuest, and SPSS relationship wording for a fitted analysis.
#'
#' @section Typical workflow:
#' - Legacy handoff:
#'   [run_mfrm_facets()] -> [build_fixed_reports()] ->
#'   [facets_output_file_bundle()].
#' - Mixed workflow:
#'   `RSM` / `PCM`:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [build_apa_outputs()] ->
#'   compatibility export only if required.
#'   bounded `GPCM`:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [reporting_checklist()] ->
#'   graph-only compatibility export only when a legacy handoff truly requires
#'   it.
#' - Compatibility-contract audit:
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
#' compatibility_alias_table("functions")
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

#' List retained compatibility aliases and preferred names
#'
#' @param scope Which alias surface to return: `"all"`, `"functions"`,
#'   `"arguments"`, `"columns"`, or `"plot_metrics"`.
#'
#' @details
#' This helper is a compact public registry of the compatibility aliases that
#' `mfrmr` intentionally keeps visible for older scripts and downstream
#' handoffs. It is meant to answer two questions quickly:
#' 1. Which old names are still accepted?
#' 2. Which package-native names should new code use instead?
#'
#' Internal soft-deprecated helpers are deliberately excluded here. This table
#' is only for retained user-facing aliases that remain part of the public
#' surface.
#'
#' @return A data.frame with one row per retained alias and columns:
#' - `Alias`
#' - `PreferredName`
#' - `Surface`
#' - `Lifecycle`
#' - `RetainedFor`
#' - `Notes`
#'
#' @section Typical workflow:
#' 1. Call `compatibility_alias_table()` when reading older scripts or reports.
#' 2. Use `PreferredName` when writing new analysis code.
#' 3. Keep the alias only when an older workflow or external handoff requires it.
#'
#' @seealso [mfrmr_compatibility_layer], [run_mfrm_facets()], [analyze_dff()],
#'   [reporting_checklist()], [fair_average_table()], [plot_fair_average()]
#' @examples
#' compatibility_alias_table()
#' compatibility_alias_table("functions")
#' compatibility_alias_table("columns")
#' @export
compatibility_alias_table <- function(scope = c("all", "functions", "arguments", "columns", "plot_metrics")) {
  scope <- match.arg(scope)

  out <- data.frame(
    Alias = c(
      "mfrmRFacets",
      "analyze_dif",
      "JMLE",
      "ReadyForAPA",
      "SE",
      "Fair(M) Average",
      "Fair(Z) Average",
      "FairM",
      "FairZ"
    ),
    PreferredName = c(
      "run_mfrm_facets",
      "analyze_dff",
      "JML",
      "DraftReady",
      "ModelSE",
      "AdjustedAverage",
      "StandardizedAdjustedAverage",
      "AdjustedAverage",
      "StandardizedAdjustedAverage"
    ),
    Surface = c(
      "function",
      "function",
      "argument",
      "column",
      "column",
      "column",
      "column",
      "plot_metric",
      "plot_metric"
    ),
    Lifecycle = c(
      "retained_alias",
      "retained_alias",
      "retained_alias",
      "retained_alias",
      "retained_alias",
      "retained_alias",
      "retained_alias",
      "retained_alias",
      "retained_alias"
    ),
    RetainedFor = c(
      "older workflow scripts",
      "earlier DIF-oriented package code",
      "historical method labels",
      "older reporting scripts",
      "older measure-table scripts",
      "FACETS-style table continuity",
      "FACETS-style table continuity",
      "legacy plot metric shortcuts",
      "legacy plot metric shortcuts"
    ),
    Notes = c(
      "Compatibility wrapper for the legacy-compatible one-shot workflow.",
      "DFF naming is preferred for many-facet workflows; the older DIF name is still accepted.",
      "Accepted by fit wrappers, but user-facing summaries and docs use JML.",
      "Backward-compatible reporting flag; values match DraftReady exactly.",
      "Backward-compatible standard-error column; ModelSE is the preferred label.",
      "Legacy adjusted-average column retained alongside the native label.",
      "Legacy standardized adjusted-average column retained alongside the native label.",
      "Accepted by plot_fair_average() as a shortcut for AdjustedAverage.",
      "Accepted by plot_fair_average() as a shortcut for StandardizedAdjustedAverage."
    ),
    stringsAsFactors = FALSE
  )

  if (identical(scope, "all")) {
    return(out)
  }

  scope_map <- c(
    functions = "function",
    arguments = "argument",
    columns = "column",
    plot_metrics = "plot_metric"
  )
  out[out$Surface %in% scope_map[[scope]], , drop = FALSE]
}
