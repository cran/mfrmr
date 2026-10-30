#' mfrmr Reporting and APA Guide
#'
#' @description
#' Package-native guide to moving from fitted model objects to
#' manuscript-draft text, tables, notes, and revision checklists in `mfrmr`.
#'
#' This guide applies to ordinary MML/JML `RSM` / `PCM` reporting, with the
#' estimator-specific limits on structural inference. Native CML and
#' comparison-only fits do not supply this consolidated reporting workflow.
#' For CML, use `summary(fit)` and [mfrm_facet_intervals()], then `summary(ci)`,
#' `plot(ci)`, [apa_table()] and `saveRDS(ci, file)`; see
#' `vignette("mfrmr-facet-intervals", package = "mfrmr")`.
#' Ordinary one-family `GPCM` fits support [reporting_checklist()],
#' [precision_review_report()], direct curve/graph and residual table helpers,
#' and caveated APA/QC/export bundles. Use [gpcm_capability_matrix()] when you
#' need the formal boundary for the current `GPCM` reporting path.
#'
#' In particular, `GPCM` [build_apa_outputs()],
#' [build_visual_summaries()], [run_qc_pipeline()],
#' [build_mfrm_manifest()], [build_mfrm_replay_script()], and
#' [export_mfrm_bundle()] outputs include explicit `gpcm_boundary` caveats.
#' Full FACETS-style score-side contract review remains blocked. Scorefile
#' export, design forecasting, diagnostic/signal-detection screening, and
#' linking synthesis use their own caveated `GPCM` routes and should not be
#' treated as automatic operational-scoring evidence.
#'
#' See `vignette("mfrmr-networks-and-reporting")` for Lamprianou-inspired network
#' views, the FACETS/ConQuest reporting crosswalk and speaking/writing examples
#' combining MFRM with G/D studies. [plot_data()] supports long-format extraction
#' and original-value cluster profiles for custom distribution plots.
#'
#' @section Start with the reporting question:
#' - "Where should I start with an existing fit?"
#'   Use [mfrm_results()] to retain the fit and matching outputs, inspect
#'   `summary(res)`, then use [mfrm_report()] for a report or
#'   [export_mfrm_results()] for an analyst archive.
#' - "What can I give to an individual rater?"
#'   Use `mfrm_report(res, style = "rater", facet = "Rater", rater = "R01",
#'   output = "html")` for native additive RSM/PCM or experimental two-family
#'   GPCM MML results. Review and retain
#'   the standalone HTML file. The complete analysis archive serves the analyst
#'   and retains source data and identifiers.
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
#' @section Editable one-family GPCM report:
#' Copy the installed `templates/gpcm-report.Rmd` using
#' `system.file("templates", "gpcm-report.Rmd", package = "mfrmr")`, then render
#' it with [rmarkdown::render()]. The default synthetic example produces an
#' HTML report, separation/Strata/reliability and fit tables, a Wright map,
#' category-probability curves, ordinal Krippendorff's alpha, CSVs and an RDS
#' analysis. Set `ratings_csv`, `rating_min`, `rating_max` and `output_dir`
#' through the template's parameters for your data. It expects columns Person,
#' Rater, Criterion and Score, with one rating per cell; repeated occasions
#' require an explicit measurement unit and agreement keys.
#' The model is one-family GPCM MML with criterion-owned steps and slopes.
#' Alpha is calculated separately from observed ratings; [mfrm_report()] does
#' not automatically include it. Fit review states and precision limitations
#' remain visible. Complete the study-context and interpretation prompts before
#' treating the rendered document as a manuscript. See
#' `vignette("mfrmr-reporting-and-apa")` for the output map and full RSM lesson.
#'
#' @section Recommended reporting route:
#' 1. Fit with [fit_mfrm()] and inspect `summary(fit)` before interpreting it.
#' 2. Calculate the diagnostics or intervals needed for the question. For a
#'    fixed-rater RSM/PCM feedback sheet, use [diagnose_mfrm()] and, for
#'    eligible MML fits, [mfrm_facet_intervals()]. A rater coefficient and a
#'    difference between raters are distinct interval targets.
#' 3. Retain them with `res <- mfrm_results(fit, diagnostics = diagnostics,
#'    intervals = list(raters = ci), compute = "never")` when these objects
#'    have been calculated. Omit attachments you do not need. Here
#'    `compute = "never"` avoids filling absent diagnostics automatically;
#'    missing sections remain explicit.
#' 4. Inspect `summary(res)$triage` and `summary(res)$next_actions`. Choose a
#'    figure through `summary(res)$plot_map` or `mfrmr_output_guide("plots")`.
#'    Use `plot(res, ...)` or [as_ggplot()] for supported views. An attached
#'    rater interval is displayed by `type = "facet_raters"`; an ordinary
#'    Wright map does not acquire that interval method automatically.
#' 5. Use `mfrm_report(res)` for the analyst's quality-control report, or
#'    `style = "rater"` with explicit `facet` and `rater` for one recipient.
#'    `output = "html"` creates a temporary file; copy it to keep it.
#' 6. Save `res` with [base::saveRDS()] for later use, or use
#'    [export_mfrm_results()] for CSVs, HTML, RDS and replay files.
#'    `preset = "starter"` also requests reports and available figures.
#'    Review `written_files` and `plot_errors` before treating export as complete.
#'
#' These functions serve different purposes; figures and reports are optional
#' branches from the saved results, not compulsory steps before saving them.
#' The worked example in `vignette("mfrmr-facet-intervals")` follows this route
#' through individual feedback, interval figures and save/reopen operations.
#'
#' @section Specialist reporting tools:
#' Use [reporting_checklist()] to inspect manuscript sections and
#' [precision_review_report()] when the strength of precision claims matters.
#' [build_apa_outputs()] supplies manuscript-draft prose and metadata;
#' [build_summary_table_bundle()] collects selected tables for [apa_table()]
#' or [export_summary_appendix()]. These narrower tools remain useful and
#' are not replaced by the general report. For `GPCM`, keep the
#' `gpcm_boundary` caveats with APA/QC/export output.
#' When strict marginal rows are available, [plot(x, type = "marginal_fit")][plot.mfrm_fit] and
#' [plot(x, type = "marginal_pairwise")][plot.mfrm_fit] support local-misfit follow-up.
#' When candidate models are compared, keep the comparison as a reporting
#'    review: [compare_mfrm()] -> [build_model_choice_review()] ->
#'    [build_summary_table_bundle()]. Treat `GPCM` as a slope-aware
#'    sensitivity route unless the study design explicitly justifies
#'    discrimination-based operational scoring.
#'
#' @section Keep each fit with its own diagnostics:
#' After changing the data, model, estimator, or fitting settings, compute
#' `diagnostics <- diagnose_mfrm(fit)` again for the new fit. Supplying the
#' previous model's diagnostics can mix its standard errors and precision
#' status with the new model's estimates. Fit summaries, precision reviews,
#' APA and visual reporting, fit plots, and fit-level export helpers reject
#' mismatched or outdated diagnostic readiness records. Matching saved
#' diagnostics remain reusable; changing a table caption or note does not
#' bypass this check.
#' For saved analyses affected by the 0.2.4 calculation changes, follow
#' "Updating saved analyses for 0.2.4" in [mfrmr_workflow_methods] before
#' rebuilding reports. Matching the source fit alone does not update old
#' diagnostic calculations.
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
#' A conservative `GPCM` reporting sequence is:
#' [fit_mfrm()] for the equal-weighting `RSM` / `PCM` reference,
#' [fit_mfrm()] for the `GPCM` sensitivity fit,
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
#' Wright and Masters, and Muraki for `GPCM`.
#'
#' The package only knows the fitted measurement objects and context supplied
#' by the analyst. It therefore cannot certify research-level JARS completeness
#' for hypotheses and their confirmatory/exploratory status, recruitment and
#' participant characteristics, ethics, sample-size rationale, missing-data
#' mechanism and exclusions, multiplicity or deviations from plan, or complete
#' data/code availability statements. Those study-level fields, effect-size and
#' uncertainty choices, statistic-specific rounding, and journal typography
#' must be reviewed and completed outside the generated template.
#' The "Manuscript coverage map" in
#' `vignette("mfrmr-reporting-and-apa", package = "mfrmr")` connects each
#' reporting topic to the relevant output and the information the author must
#' supply. It includes the rating assignment and training, missingness,
#' estimation settings, uncertainty, category functioning, and the distinct
#' meanings of separation reliability and observed agreement. The vignette
#' also demonstrates a question-to-result explanation with actual estimates.
#' Its short Methods and Results example includes rating-row denominators,
#' residual counts, step uncertainty, and the difference between
#' separation-based reliability and posterior-variance EAP reliability.
#' The observed-versus-fair-score example states its reference environment
#' and distinguishes logit-measure uncertainty from fair-score uncertainty.
#' Verify author-provided context labels, such as `scale_desc`, against the
#' study record; generated-content checks do not establish their accuracy.
#'
#' Appelbaum, M., Cooper, H., Kline, R. B., Mayo-Wilson, E., Nezu, A. M., and
#' Rao, S. M. (2018). Journal article reporting standards for quantitative
#' research in psychology: The APA Publications and Communications Board task
#' force report. *American Psychologist, 73*(1), 3-25.
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
#'   follow-up should be [plot(x, type = "qc")][plot.mfrm_fit], [plot(x, type = "marginal_fit")][plot.mfrm_fit],
#'   [plot(x, type = "residual_pca")][plot.mfrm_fit], [plot(x, type = "bias")][plot.mfrm_fit], or another public plot.
#' - Use [visual_reporting_template()] to draft visual captions and
#'   results-sentence starters, but do not paste the skeletons without checking
#'   the actual fit, diagnostics, and study context.
#' - Before formal inferential claims, review the fit decision, precision
#'   checks, and restrictions for the particular statistic. A `model_based`
#'   tier alone does not establish that the claim is supported.
#' - Keep bias and differential-functioning outputs in screening language
#'   unless the current precision layer and linking evidence justify stronger
#'   claims.
#' - Treat `DraftReady` (and the legacy alias `ReadyForAPA`) as a
#'   drafting-readiness flag, not as a substitute for methodological review.
#' - Rebuild APA outputs after major model changes instead of editing old text
#'   by hand.
#' - For `GPCM`, use APA/QC/export helpers only as caveated
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
#' @section Additional task-specific workflows:
#' Use these when the reporting question calls for extra manuscript or
#' appendix assembly, after choosing the supported route above.
#' - Manuscript-first route:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [reporting_checklist()] ->
#'   [build_apa_outputs()] -> [build_summary_table_bundle()] -> `summary()` /
#'   `plot()` -> [apa_table()], [export_summary_appendix()], or
#'   [export_mfrm_bundle()](include = c("summary_tables", "html")).
#'   For `RSM` / `PCM` final reports, prefer `method = "MML"` and
#'   `diagnostic_mode = "both"` in the diagnostics step.
#'   For `GPCM`, use the same fit-based reporting/export family only
#'   as caveated sensitivity-reporting output and inspect its `gpcm_boundary`
#'   rows before writing claims.
#' - Appendix-first route:
#'   [facet_statistics_report()] -> [apa_table()] ->
#'   [build_visual_summaries()] -> [build_apa_outputs()].
#' - Precision-sensitive route:
#'   [diagnose_mfrm()] -> [precision_review_report()] ->
#'   [reporting_checklist()] -> [build_apa_outputs()].
#' - `GPCM` route:
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
#' - For the `GPCM` support statement, see [gpcm_capability_matrix].
#' - For a longer walkthrough, see
#'   `vignette("mfrmr-reporting-and-apa", package = "mfrmr")`.
#'
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("mfrmr_reporting_and_apa", run.dontrun = TRUE)
#' \dontrun{
#' # Load the package and example ratings
#' library(mfrmr)
#' toy <- load_mfrmr_data("example_operational")
#'
#' # Fit the model
#' fit <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   model = "RSM"
#' )
#'
#' # Check which report sections have supporting evidence
#' diagnostics <- diagnose_mfrm(fit)
#' checklist <- reporting_checklist(fit, diagnostics = diagnostics)
#' checklist$section_summary
#' subset(checklist$checklist, !DraftReady,
#'        c("Section", "Item", "NextAction"))
#'
#' # Format the per-facet distribution summary as a table
#' results <- summary(fit, diagnostics = diagnostics)
#' tbl <- apa_table(results, which = "facet_overview")
#' tbl # Prints the table, caption, and note; review them before using in a paper
#'
#' # For individual estimates, see as.data.frame(fit)
#' }
#' @name mfrmr_reporting_and_apa
NULL
