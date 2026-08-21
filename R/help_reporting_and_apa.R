#' mfrmr Reporting and APA Guide
#'
#' @description
#' Package-native guide to moving from fitted model objects to
#' manuscript-draft text, tables, notes, and revision checklists in `mfrmr`.
#'
#' This guide currently applies fully to diagnostics-based `RSM` / `PCM`
#' workflows. Bounded `GPCM` fits support [reporting_checklist()],
#' [precision_review_report()], direct curve/graph and residual table helpers,
#' and caveated APA/QC/export bundles. Use [gpcm_capability_matrix()] when you
#' need the formal boundary for the current `GPCM` reporting path.
#'
#' In particular, bounded `GPCM` [build_apa_outputs()],
#' [build_visual_summaries()], [run_qc_pipeline()],
#' [build_mfrm_manifest()], [build_mfrm_replay_script()], and
#' [export_mfrm_bundle()] outputs include explicit `gpcm_boundary` caveats.
#' Full FACETS-style score-side contract review remains blocked. Scorefile
#' export, design forecasting, diagnostic/signal-detection screening, and
#' linking synthesis use their own caveated `GPCM` routes and should not be
#' treated as automatic operational-scoring evidence.
#'
#' @section Start with the reporting question:
#' - "Which parts of this run are ready to draft, and with what caveats?"
#'   Use [reporting_checklist()].
#' - "How should I phrase the model, fit, and precision sections?"
#'   For `RSM` / `PCM`, use [build_apa_outputs()].
#' - "Which tables should I hand off to a manuscript or appendix?"
#'   Use [build_summary_table_bundle()], [export_summary_appendix()],
#'   [apa_table()], and
#'   [facet_statistics_report()].
#' - "How do I explain model-based vs exploratory precision?"
#'   Use [precision_review_report()] and `summary(diagnose_mfrm(...))`.
#' - "Which caveats need to appear in the write-up?"
#'   Use [reporting_checklist()] first, then [build_apa_outputs()].
#' - "How should I report candidate-model comparisons?"
#'   Use [compare_mfrm()] for the same-data comparison table, then
#'   [build_model_choice_review()] and [build_summary_table_bundle()] for
#'   cautious model-role, route-boundary, and wording tables.
#' - "How should I start figure captions or visual-results wording?"
#'   Use [visual_reporting_template()] for conservative caption and results
#'   sentence starters, then verify availability with
#'   `reporting_checklist()$visual_scope`.
#'
#' @section Recommended reporting route:
#' 1. Fit with [fit_mfrm()].
#' 2. Build diagnostics with [diagnose_mfrm()].
#' 3. Review precision strength with [precision_review_report()] when
#'    inferential language matters.
#' 4. Run [reporting_checklist()] to identify missing sections, caveats, and
#'    next actions. Use the `"Visual Displays"` rows as the figure-routing
#'    layer for the current run.
#' 5. When strict marginal rows are available, follow up with
#'    [plot_marginal_fit()] and [plot_marginal_pairwise()] before finalizing
#'    the narrative around local misfit.
#' 6. Create manuscript-draft prose and metadata with [build_apa_outputs()].
#'    For bounded `GPCM`, treat the APA/QC/export stack as caveated
#'    sensitivity-reporting output and keep its `gpcm_boundary` visible.
#' 7. Convert summary outputs to reusable table bundles with
#'    [build_summary_table_bundle()], review the bundle with `summary()` /
#'    `plot()`, then convert specific components to handoff tables with
#'    [apa_table()] or export them directly with [export_summary_appendix()].
#' 8. When candidate models are compared, keep the comparison as a reporting
#'    review: [compare_mfrm()] -> [build_model_choice_review()] ->
#'    [build_summary_table_bundle()]. Treat bounded `GPCM` as a slope-aware
#'    sensitivity route unless the study design explicitly justifies
#'    discrimination-based operational scoring.
#'
#' @section Model-comparison reporting route:
#' Use [compare_mfrm()] to build the candidate-model table and inspect
#' `ICComparable`, `ComparisonBasis`, and any nesting warnings before reading
#' information criteria. Automatic ranking requires the current contract and
#' a shared q>=31 grid; q<31 retains raw screening/review criteria only, and a
#' close decision still needs a prespecified denser common-grid check. Then use
#' [build_model_choice_review()] to attach the
#' comparison to explicit model roles, downstream-route boundaries, wording
#' templates, and optional [build_weighting_review()] output. Convert that
#' review with [build_summary_table_bundle()] when a manuscript appendix,
#' coauthor handoff, or HTML export needs stable table names.
#'
#' A conservative bounded-`GPCM` reporting sequence is:
#' [fit_mfrm()] for the equal-weighting `RSM` / `PCM` reference,
#' [fit_mfrm()] for the bounded `GPCM` sensitivity fit,
#' [compare_mfrm()], [build_model_choice_review()],
#' [build_summary_table_bundle()], then [export_summary_appendix()] or
#' [export_mfrm_bundle()]. Do not use `AIC`, `BIC`, or log-likelihood alone as
#' an automatic operational-scoring decision.
#'
#' @section Latent-regression reporting route:
#' Active latent-regression fits expose their reporting surface through
#' `summary(fit)$population_overview`,
#' `summary(fit)$population_coefficients`,
#' `summary(fit)$population_coding`, and fit-level `caveats`. Report those
#' coefficients as conditional-normal population-model parameters, not as a
#' post-hoc regression on EAP or MLE scores. Also report the
#' `population_formula`, coding/contrast information, `population_policy`, and
#' omitted-person or omitted-row counts when complete-case handling was used.
#'
#' Prediction-side helpers [predict_mfrm_units()] and
#' [sample_mfrm_plausible_values()] can carry the fitted population model into
#' future-unit scoring and plausible-value draws. The supported route is
#' one-dimensional `MML` for `RSM` / `PCM`; avoid stronger
#' claims about multidimensional latent regression, Wald tests, posterior
#' predictive checking, or full external-engine equivalence unless those checks
#' were performed outside this helper family.
#'
#' @section Publication-readiness boundary:
#' `mfrmr` can provide a defensible measurement-output trail for a manuscript:
#' fitted model summaries, diagnostic tables, precision review, report
#' templates, APA table metadata, figure-routing guidance, and reproducible
#' exports. It does not decide whether a specific journal claim is warranted.
#' For high-stakes or selective journals, use the package outputs together
#' with the study design, measurement rationale, primary citations, sensitivity
#' checks, and substantive argument for the target field.
#'
#' Treat `DraftReady`, `ReadyForAPA`, `ClaimStrength`, and report-template rows
#' as drafting and caveat-routing aids. They are not formal acceptance rules,
#' proof of validity, or a substitute for peer-review judgment. Before copying
#' text, inspect `mfrm_report(res, style = "apa")$first_screen`,
#' `$claim_readiness`, `$report_gaps`, and `$template_index`.
#'
#' @section Standards basis and boundary:
#' The manuscript helpers are APA-oriented drafting aids informed by the
#' *Publication Manual of the American Psychological Association* (7th ed.)
#' and the quantitative Journal Article Reporting Standards (JARS-Quant;
#' Appelbaum et al., 2018). The MFRM-specific prompts also draw on the model and
#' diagnostic sources returned by `reporting_checklist(...,
#' include_references = TRUE)`, including Eckes, Myford and Wolfe, Linacre,
#' Wright and Masters, and Muraki for bounded `GPCM`.
#'
#' The package only knows the fitted measurement objects and context supplied
#' by the analyst. It therefore cannot certify research-level JARS completeness
#' for hypotheses and their confirmatory/exploratory status, recruitment and
#' participant characteristics, ethics, sample-size rationale, missing-data
#' mechanism and exclusions, multiplicity or deviations from plan, or complete
#' data/code availability statements. Those study-level fields, effect-size and
#' uncertainty choices, statistic-specific rounding, and journal typography
#' must be reviewed and completed outside the generated template.
#'
#' Appelbaum, M., Cooper, H., Kline, R. B., Mayo-Wilson, E., Nezu, A. M., and
#' Rao, S. M. (2018). Journal article reporting standards for quantitative
#' research in psychology. *American Psychologist, 73*(1), 3-25.
#' \doi{10.1037/amp0000191}
#'
#' @section Which helper answers which task:
#' \describe{
#'   \item{[reporting_checklist()]}{Turns current analysis objects into a
#'   prioritized revision guide with `DraftReady`, `Priority`, and
#'   `NextAction`. `DraftReady` means "ready to draft with the documented
#'   caveats"; `ReadyForAPA` is retained as a backward-compatible alias, and
#'   neither field means "formal inference is automatically justified". The
#'   `"Visual Displays"` rows also mirror the public plot family, so the
#'   checklist doubles as a figure-routing surface.}
#'   \item{[build_apa_outputs()]}{Builds shared-contract prose, table notes,
#'   captions, and a section map from the current fit and diagnostics.}
#'   \item{[build_summary_table_bundle()]}{Turns supported `summary()` outputs
#'   into named `data.frame` tables plus an index for manuscript or appendix
#'   handoff, and now also supports bundle-level `summary()` / `plot()` for
#'   role coverage and numeric QC.}
#'   \item{[export_summary_appendix()]}{Writes those documented summary-table
#'   bundles to CSV and optional HTML appendix artifacts without requiring a
#'   full fit-based export bundle.}
#'   \item{[apa_table()]}{Produces reproducible base-R tables with APA-oriented
#'   labels, notes, and captions.}
#'   \item{[precision_review_report()]}{Summarizes whether precision claims are
#'   model-based, hybrid, or exploratory.}
#'   \item{[facet_statistics_report()]}{Provides facet-level summaries that
#'   often feed result tables and appendix material.}
#'   \item{[build_visual_summaries()]}{Prepares publication-oriented figure
#'   data that can be cited from the report text.}
#'   \item{[visual_reporting_template()]}{Provides conservative figure
#'   placement, caption-starter, results-wording, and overclaim-avoidance
#'   guidance for public visual helpers.}
#' }
#'
#' @section Practical reporting rules:
#' - Treat [reporting_checklist()] as the gap finder and
#'   [build_apa_outputs()] as the writing engine.
#' - Use the checklist's `"Visual Displays"` rows to decide whether the next
#'   follow-up should be [plot_qc_dashboard()], [plot_marginal_fit()],
#'   [plot_residual_pca()], [plot_bias_interaction()], or another public plot.
#' - Use [visual_reporting_template()] to draft visual captions and
#'   results-sentence starters, but do not paste the skeletons without checking
#'   the actual fit, diagnostics, and study context.
#' - Phrase formal inferential claims only when the precision tier is
#'   model-based.
#' - Keep bias and differential-functioning outputs in screening language
#'   unless the current precision layer and linking evidence justify stronger
#'   claims.
#' - Treat `DraftReady` (and the legacy alias `ReadyForAPA`) as a
#'   drafting-readiness flag, not as a substitute for methodological review.
#' - Rebuild APA outputs after major model changes instead of editing old text
#'   by hand.
#' - For bounded `GPCM`, use APA/QC/export helpers only as caveated
#'   sensitivity-reporting surfaces and keep full FACETS-style score-side
#'   review outside this route.
#'
#' @section Fit-to-HTML reporting bundle:
#' When the user already has a fitted object and wants a local report folder in
#' one call, use [export_mfrm_bundle()] directly:
#' `export_mfrm_bundle(fit, include = c("core_tables", "checklist",
#' "dashboard", "apa", "summary_tables", "manifest", "script", "html"))`.
#' This route computes missing diagnostics, writes CSV/text/replay artifacts,
#' and creates a lightweight HTML summary without requiring a prior
#' [mfrm_results()] object. Use [mfrm_results()] and [mfrm_report()] first when
#' the goal is interactive triage or report-readiness review; use
#' [export_mfrm_bundle()] when the goal is a file bundle for a project folder,
#' coauthor handoff, or supplementary-methods archive. The bundle is not
#' deidentified; review every file under the study's data-handling policy before
#' any handoff.
#'
#' @section Typical workflow:
#' - Manuscript-first route:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [reporting_checklist()] ->
#'   [build_apa_outputs()] -> [build_summary_table_bundle()] -> `summary()` /
#'   `plot()` -> [apa_table()], [export_summary_appendix()], or
#'   [export_mfrm_bundle()](include = c("summary_tables", "html")).
#'   For `RSM` / `PCM` final reports, prefer `method = "MML"` and
#'   `diagnostic_mode = "both"` in the diagnostics step.
#'   For bounded `GPCM`, use the same fit-based reporting/export family only
#'   as caveated sensitivity-reporting output and inspect its `gpcm_boundary`
#'   rows before writing claims.
#' - Appendix-first route:
#'   [facet_statistics_report()] -> [apa_table()] ->
#'   [build_visual_summaries()] -> [build_apa_outputs()].
#' - Precision-sensitive route:
#'   [diagnose_mfrm()] -> [precision_review_report()] ->
#'   [reporting_checklist()] -> [build_apa_outputs()].
#' - bounded `GPCM` route:
#'   [diagnose_mfrm()] -> [precision_review_report()] ->
#'   [reporting_checklist()] -> direct residual/category/information helpers ->
#'   caveated [build_apa_outputs()], [build_visual_summaries()],
#'   [run_qc_pipeline()], or [export_mfrm_bundle()] as needed.
#' - Model-comparison route:
#'   [compare_mfrm()] -> [build_model_choice_review()] ->
#'   [build_summary_table_bundle()] -> [export_summary_appendix()] or
#'   [export_mfrm_bundle()](include = c("summary_tables", "html")).
#'
#' @section Companion guides:
#' - For report/table selection, see [mfrmr_reports_and_tables].
#' - For end-to-end analysis routes, see [mfrmr_workflow_methods].
#' - For visual follow-up, see [mfrmr_visual_diagnostics].
#' - For the bounded `GPCM` support statement, see [gpcm_capability_matrix].
#' - For a longer walkthrough, see
#'   `vignette("mfrmr-reporting-and-apa", package = "mfrmr")`.
#'
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' # A balanced slice retains every Rater and Criterion while running quickly.
#' toy <- toy[toy$Person %in% unique(toy$Person)[1:12], , drop = FALSE]
#' fit <- fit_mfrm(
#'   toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   quad_points = 7,
#'   maxit = 30
#' )
#' diag <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "both")
#'
#' checklist <- reporting_checklist(fit, diagnostics = diag)
#' visual_reporting_template("manuscript")[, c("FigureFamily", "CaptionSkeleton")]
#' head(checklist$checklist[, c("Section", "Item", "DraftReady", "NextAction")])
#' subset(
#'   checklist$checklist,
#'   Section == "Visual Displays",
#'   c("Item", "Available", "NextAction")
#' )
#'
#' apa <- build_apa_outputs(fit, diagnostics = diag)
#' apa$section_map[, c("SectionId", "Available")]
#'
#' tbl <- apa_table(fit, which = "summary")
#' tbl$caption
#' bundle <- build_summary_table_bundle(checklist)
#' bundle$table_index
#' apa_from_bundle <- apa_table(bundle, which = "section_summary")
#' apa_from_bundle$caption
#'
#' report_bundle <- export_mfrm_bundle(
#'   fit,
#'   diagnostics = diag,
#'   output_dir = tempdir(),
#'   prefix = "mfrmr_report_bundle",
#'   include = c("core_tables", "checklist", "apa", "summary_tables", "html"),
#'   overwrite = TRUE,
#'   acknowledge_sensitive = TRUE
#' )
#' report_bundle$summary[, c("FilesWritten", "HtmlWritten")]
#' }
#'
#' @name mfrmr_reporting_and_apa
NULL
