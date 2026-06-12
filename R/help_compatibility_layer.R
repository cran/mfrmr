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
#' SPSS estimators, or claim SPSS numerical equivalence.
#'
#' @section When to use this layer:
#' - You are reproducing an older workflow that expects one-shot wrappers.
#' - You need fixed-width text blocks for console, logs, or archival handoff.
#' - You need graphfile or scorefile style outputs for downstream legacy tools.
#' - You are checking column coverage and metric consistency against a
#'   FACETS-style output contract.
#'
#' @section When not to use this layer:
#' - For standard estimation, use [fit_mfrm()] plus [diagnose_mfrm()].
#' - For report bundles, use [mfrmr_reports_and_tables].
#' - For manuscript text, use [build_apa_outputs()] and [reporting_checklist()].
#' - For visual follow-up, use [mfrmr_visual_diagnostics].
#'
#' @section Compatibility map:
#' \describe{
#'   \item{[facets_positioning_guide()]}{User-facing wording for the package's
#'   relationship to FACETS. Use before describing compatibility outputs in a
#'   report or migration note.}
#'   \item{[run_mfrm_facets()]}{One-shot legacy-compatible wrapper that fits,
#'   diagnoses, and returns key tables in one object.}
#'   \item{[mfrmRFacets()]}{Alias for [run_mfrm_facets()] kept for continuity.}
#'   \item{[build_fixed_reports()]}{Fixed-width interaction and pairwise text
#'   blocks. Best when a text-only compatibility artifact is required.}
#'   \item{[facets_output_file_bundle()]}{Graphfile/scorefile style CSV and
#'   fixed-width exports for legacy pipelines.}
#'   \item{[write_mfrm_residual_file()]}{Package-native observation-level
#'   residual CSV/TSV export for reviewer, spreadsheet, or external QC handoff.}
#'   \item{[write_mfrm_subset_file()]}{Package-native connected-subset summary
#'   and node-membership CSV/TSV export for linking review handoff.}
#'   \item{[facets_output_contract_review()]}{Column and metric review against
#'   the FACETS-style output-contract specification. Use only when an explicit
#'   output-contract review is part of the task; it checks the package output contract
#'   and does not imply external FACETS equivalence.}
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
#' - For residual or subset handoff, prefer [write_mfrm_residual_file()] and
#'   [write_mfrm_subset_file()] over reusing graph/score compatibility exports.
#' - Instead of [facets_output_contract_review()] for routine QA, prefer:
#'   [reference_case_review()] for package-native completeness review or
#'   [reference_case_benchmark()] for packaged benchmark cases.
#'
#' @section Practical migration rules:
#' - Start FACETS-facing reports with [facets_positioning_guide()] when readers
#'   might otherwise assume FACETS numerical reproduction.
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
#' @section Retained table-field names:
#' - `row_review` is the data-quality table field used to document row
#'   filtering. FACETS-style column and metric contract results are exposed as
#'   `column_review` and `metric_checks`.
#' - Prediction outputs expose row-preparation and person-omission
#'   traceability as `row_review` and `population_review`.
#' - `hierarchical_review`, `shrinkage_review`, and `nesting_review` are
#'   manifest or model-comparison traceability fields. They are not callable
#'   helper names; user-facing helper names use review/check terminology.
#'
#' @section Typical workflow:
#' - Legacy handoff:
#'   [run_mfrm_facets()] -> [build_fixed_reports()] ->
#'   [facets_output_file_bundle()] plus [write_mfrm_residual_file()] /
#'   [write_mfrm_subset_file()] only when those standalone files are needed.
#' - Mixed workflow:
#'   `RSM` / `PCM`:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [build_apa_outputs()] ->
#'   compatibility export only if required.
#'   bounded `GPCM`:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [reporting_checklist()] ->
#'   graph-only compatibility export only when a legacy handoff truly requires
#'   it.
#' - FACETS output-contract review:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [facets_output_contract_review()].
#'
#' @section Companion guides:
#' - For FACETS coverage and boundary wording, see
#'   [facets_positioning_guide()] and [facets_feature_coverage()].
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
#'   maxit = 30
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
#'   `"arguments"`, `"fields"`, `"columns"`, or `"plot_metrics"`.
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
#' - `RemovalPlan`
#' - `Notes`
#'
#' @section Typical workflow:
#' 1. Call `compatibility_alias_table()` when reading older scripts or reports.
#' 2. Use `PreferredName` when updating older analysis code.
#' 3. Prefer the package-native name in all new outputs and scripts.
#'
#' @seealso [mfrmr_compatibility_layer], [run_mfrm_facets()], [analyze_dff()],
#'   [reporting_checklist()], [fair_average_table()], [plot_fair_average()]
#' @examples
#' compatibility_alias_table()
#' compatibility_alias_table("functions")
#' compatibility_alias_table("fields")
#' compatibility_alias_table("columns")
#' @export
compatibility_alias_table <- function(scope = c("all", "functions", "arguments", "fields", "columns", "plot_metrics")) {
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
    RemovalPlan = c(
      "No scheduled removal.",
      "No scheduled removal.",
      "No scheduled removal.",
      "No scheduled removal.",
      "No scheduled removal.",
      "No scheduled removal.",
      "No scheduled removal.",
      "No scheduled removal.",
      "No scheduled removal."
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
    fields = "field",
    columns = "column",
    plot_metrics = "plot_metric"
  )
  out[out$Surface %in% scope_map[[scope]], , drop = FALSE]
}

#' Extract canonical review components
#'
#' @param x A fitted `mfrm_fit`, diagnostics object, summary object, or
#'   compatible list containing the requested review component.
#' @param required Logical. If `TRUE`, error when no compatible review component
#'   is found. If `FALSE`, return `NULL` instead.
#'
#' @details
#' `anchor_review()` returns the fitted object's `config$anchor_review`
#' component. `precision_review()` returns the diagnostics
#' `precision_review` table. These helpers intentionally do not search older
#' field names.
#'
#' @return
#' `anchor_review()` returns an `mfrm_anchor_review`-like object when available.
#' `precision_review()` returns the precision-review table when available.
#' If `required = FALSE` and no component is available, both helpers return
#' `NULL`.
#'
#' @examples
#' fit_like <- list(config = list(
#'   anchor_review = structure(list(issue_counts = data.frame()),
#'                             class = "mfrm_anchor_review")
#' ))
#' anchor_review(fit_like)
#' diag_like <- list(
#'   precision_review = data.frame(Check = "SE", Status = "ok", Detail = "Model-based")
#' )
#' precision_review(diag_like)
#' @name review_accessors
#' @export
anchor_review <- function(x, required = TRUE) {
  required <- isTRUE(required)

  if (inherits(x, "mfrm_anchor_review")) {
    return(x)
  }
  if (!is.list(x)) {
    stop("`x` must be an `mfrm_fit`, anchor-review object, or compatible list.", call. = FALSE)
  }

  config <- x$config %||% x
  out <- config$anchor_review
  if (is.null(out)) {
    if (!required) return(NULL)
    stop(
      "No anchor-review component was found. Use `fit$config$anchor_review`, ",
      "or run `review_mfrm_anchors()`.",
      call. = FALSE
    )
  }
  out
}

#' @rdname review_accessors
#' @export
precision_review <- function(x, required = TRUE) {
  required <- isTRUE(required)

  if (inherits(x, "mfrm_precision_review")) {
    return(x)
  }
  if (!is.list(x)) {
    stop("`x` must be an `mfrm_diagnostics`, precision-review object, or compatible list.", call. = FALSE)
  }

  out <- x$precision_review
  if (is.null(out)) {
    if (!required) return(NULL)
    stop(
      "No precision-review component was found. Run `diagnose_mfrm()` and use ",
      "`diagnostics$precision_review`.",
      call. = FALSE
    )
  }
  out
}
