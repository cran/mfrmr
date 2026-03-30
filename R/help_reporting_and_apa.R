#' mfrmr Reporting and APA Guide
#'
#' @description
#' Package-native guide to moving from fitted model objects to
#' manuscript-draft text, tables, notes, and revision checklists in `mfrmr`.
#'
#' @section Start with the reporting question:
#' - "Which parts of this run are draft-complete, and with what caveats?"
#'   Use [reporting_checklist()].
#' - "How should I phrase the model, fit, and precision sections?"
#'   Use [build_apa_outputs()].
#' - "Which tables should I hand off to a manuscript or appendix?"
#'   Use [apa_table()] and [facet_statistics_report()].
#' - "How do I explain model-based vs exploratory precision?"
#'   Use [precision_audit_report()] and `summary(diagnose_mfrm(...))`.
#' - "Which caveats need to appear in the write-up?"
#'   Use [reporting_checklist()] first, then [build_apa_outputs()].
#'
#' @section Recommended reporting route:
#' 1. Fit with [fit_mfrm()].
#' 2. Build diagnostics with [diagnose_mfrm()].
#' 3. Review precision strength with [precision_audit_report()] when
#'    inferential language matters.
#' 4. Run [reporting_checklist()] to identify missing sections, caveats, and
#'    next actions.
#' 5. Create manuscript-draft prose and metadata with [build_apa_outputs()].
#' 6. Convert specific components to handoff tables with [apa_table()].
#'
#' @section Which helper answers which task:
#' \describe{
#'   \item{[reporting_checklist()]}{Turns current analysis objects into a
#'   prioritized revision guide with `DraftReady`, `Priority`, and
#'   `NextAction`. `DraftReady` means "ready to draft with the documented
#'   caveats"; `ReadyForAPA` is retained as a backward-compatible alias, and
#'   neither field means "formal inference is automatically justified".}
#'   \item{[build_apa_outputs()]}{Builds shared-contract prose, table notes,
#'   captions, and a section map from the current fit and diagnostics.}
#'   \item{[apa_table()]}{Produces reproducible base-R tables with APA-oriented
#'   labels, notes, and captions.}
#'   \item{[precision_audit_report()]}{Summarizes whether precision claims are
#'   model-based, hybrid, or exploratory.}
#'   \item{[facet_statistics_report()]}{Provides facet-level summaries that
#'   often feed result tables and appendix material.}
#'   \item{[build_visual_summaries()]}{Prepares publication-oriented figure
#'   payloads that can be cited from the report text.}
#' }
#'
#' @section Practical reporting rules:
#' - Treat [reporting_checklist()] as the gap finder and
#'   [build_apa_outputs()] as the writing engine.
#' - Phrase formal inferential claims only when the precision tier is
#'   model-based.
#' - Keep bias and differential-functioning outputs in screening language
#'   unless the current precision layer and linking evidence justify stronger
#'   claims.
#' - Treat `DraftReady` (and the legacy alias `ReadyForAPA`) as a
#'   drafting-readiness flag, not as a substitute for methodological review.
#' - Rebuild APA outputs after major model changes instead of editing old text
#'   by hand.
#'
#' @section Typical workflow:
#' - Manuscript-first route:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [reporting_checklist()] ->
#'   [build_apa_outputs()] -> [apa_table()].
#' - Appendix-first route:
#'   [facet_statistics_report()] -> [apa_table()] ->
#'   [build_visual_summaries()] -> [build_apa_outputs()].
#' - Precision-sensitive route:
#'   [diagnose_mfrm()] -> [precision_audit_report()] ->
#'   [reporting_checklist()] -> [build_apa_outputs()].
#'
#' @section Companion guides:
#' - For report/table selection, see [mfrmr_reports_and_tables].
#' - For end-to-end analysis routes, see [mfrmr_workflow_methods].
#' - For visual follow-up, see [mfrmr_visual_diagnostics].
#' - For a longer walkthrough, see
#'   `vignette("mfrmr-reporting-and-apa", package = "mfrmr")`.
#'
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(
#'   toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   maxit = 25
#' )
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#'
#' checklist <- reporting_checklist(fit, diagnostics = diag)
#' head(checklist$checklist[, c("Section", "Item", "DraftReady", "NextAction")])
#'
#' apa <- build_apa_outputs(fit, diagnostics = diag)
#' names(apa$section_map)
#'
#' tbl <- apa_table(fit, which = "summary")
#' tbl$caption
#' }
#'
#' @name mfrmr_reporting_and_apa
NULL
