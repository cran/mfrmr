#' mfrmr Reports and Tables Map
#'
#' @description
#' Quick guide to choosing the right report or table helper in `mfrmr`.
#' Use this page when you know the reporting question but have not yet decided
#' which bundle, table, or reporting helper to call.
#'
#' @section Recommended report route:
#' For ordinary MML/JML `mfrm_fit` results, start here. The Examples section
#' follows the same route with bundled data.
#' 1. Inspect `summary(fit)$decision` before interpreting the estimates.
#' 2. Collect `res <- mfrm_results(fit)`. By default it computes missing
#'    requested diagnostics. Pass existing matching checks through
#'    `diagnostics = diagnostics` to reuse them; `compute = "never"` leaves
#'    missing diagnostics unavailable. Saved intervals are explicit attachments.
#' 3. Read `summary(res)$triage` and `$next_actions`; use `$plot_map` to find
#'    supported figures. A failed or unavailable component needs its own
#'    explanation, not an automatically substituted analysis.
#' 4. Choose `mfrm_report(res)` for an analyst report, `plot(res, ...)` for
#'    a supported figure, or [export_mfrm_results()] for tables and an archive.
#'    Save the complete object with `saveRDS(res, file)` for later R use.
#'    Reopen with `readRDS(file)`; reading or displaying saved outputs does not
#'    refit the model. Reports and figures are optional, not required to save.
#'
#' Native CML and comparison-only fits do not supply ordinary consolidated
#' reports. Read `summary(fit)`; for CML use [mfrm_facet_intervals()], then
#' `summary(ci)`, `plot(ci)`, [apa_table()] and `saveRDS(ci, file)`.
#' G/D-study and observed-agreement results also have their own summary and
#' saving routes. For GPCM or extended measurement models, [mfrm_results()]
#' states the supported sections and attachments; model-specific uncertainty
#' limits remain in the collected results.
#'
#' Individual helpers below answer additional questions. [build_apa_outputs()]
#' drafts manuscript text; [build_summary_table_bundle()] collects selected
#' tables for [export_summary_appendix()]. Neither is a prerequisite for
#' [mfrm_report()]. Use `mfrmr_output_guide("beginner")` for a compact
#' question-to-function map, or the full guide for specialist routes.
#'
#' @section Choose what to save:
#' - [export_mfrm_tables()]: core CSV tables from a fitted model. Supplied
#'   matching diagnostics enrich the tables; this call does not compute them.
#' - [export_mfrm_results()]: selected tables/reports and, by default, the
#'   complete results RDS plus a script that reloads it without refitting.
#' - [export_mfrm_bundle()]: fit-level tables, optional diagnostics and a
#'   refitting recipe. Its `"script"` component reruns the analysis; it does
#'   not preserve a fitted R object. Supply the original data for a colocated
#'   replay CSV, or replace the script's data-path placeholder yourself.
#' - [export_summary_appendix()]: selected summary tables for a manuscript.
#'   Display rounding follows `digits`; retain the source R object for reuse.
#'
#' No export certifies estimates as usable or changes their uncertainty target.
#' Inspect availability and caveat fields along with estimates. For G/D-study
#' or native CML objects, use `saveRDS()` for the complete object and `write.csv()`
#' for a selected supported summary table; the MFRM exporters do not widen their
#' accepted model classes.
#'
#' @section Start with the question:
#' - "How can I read assessment-plan estimates and uncertainty together?"
#'   Use [summary.mfrm_multivariate_d_study()] with `intervals = TRUE` and
#'   `assumption = "normal"`, or supply saved intervals. Its `table`,
#'   `intervals` and `notes` retain methods, design and interpretation limits.
#' - "How should I document the model setup and run settings?"
#'   Use [specifications_report()].
#' - "Was data filtered, dropped, or mapped in unexpected ways?"
#'   Use [data_quality_report()] and [describe_mfrm_data()].
#' - "Did estimation converge cleanly and how formal is the precision layer?"
#'   Use [estimation_iteration_report()] and [precision_review_report()].
#' - "Which facets are measurable, variable, or weakly separated?"
#'   Use [facet_statistics_report()], [measurable_summary_table()], and
#'   [facets_chisq_table()].
#' - "Are score categories functioning in a usable sequence?"
#'   Use [rating_scale_table()], [category_structure_report()], and
#'   [category_curves_report()].
#' - "Is the design linked well enough across subsets, forms, or waves?"
#'   Use [subset_connectivity_report()] and [plot(x)][plot.mfrm_anchor_drift].
#' - "What should go into the manuscript text and tables?"
#'   For ordinary MML/JML `RSM` / `PCM`, use [reporting_checklist()], [build_apa_outputs()],
#'   and [build_summary_table_bundle()] or [export_summary_appendix()]. For
#'   `GPCM`, use the same route only where
#'   [gpcm_capability_matrix()] marks it as `supported_with_caveat`: direct
#'   table/plot helpers, summary-table appendix export, caveated
#'   [build_apa_outputs()], and caveated [export_mfrm_bundle()] are available
#'   with a `gpcm_boundary`; score-side exports and design-forecasting evidence
#'   use their own caveated or blocked `GPCM` routes.
#' - "How can I estimate my model in Stan or brms?"
#'   [export_mfrm_model()] writes code, data, mappings and explicit Bayesian
#'   priors for supported ordinary MML specifications. Export does not run a
#'   sampler; posterior credible intervals are distinct from MML intervals.
#'   After sampling, use [read_mfrm_posterior()] to read the saved summary,
#'   inspect its sampling review, and continue with `plot()`, `confint()` and
#'   [mfrm_report()]. The summary does not contain joint posterior draws.
#' - "Did a simulation recover the known generating parameters well enough?"
#'   Use [simulate_mfrm_recovery()] for the recovery study,
#'   [assess_mfrm_recovery()] for the adequacy checklist, and then
#'   [build_summary_table_bundle()] or [export_summary_appendix()] for the
#'   appendix handoff.
#'
#' @section Which output answers which question:
#' \describe{
#'   \item{[specifications_report()]}{Documents model type, estimation method,
#'   anchors, and core run settings. Best for method sections and
#'   reproducibility records.}
#'   \item{[data_quality_report()]}{Summarizes retained and dropped rows,
#'   missingness, and unknown elements. Best for data cleaning narratives.}
#'   \item{[estimation_iteration_report()]}{Shows replayed convergence
#'   trajectories. Best for diagnosing slow or unstable estimation.}
#'   \item{[precision_review_report()]}{Summarizes whether `SE`, `CI`, and
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
#'   \item{[category_curves_report()]}{Returns category-probability,
#'   cumulative-probability, expected-ogive, total-information, and
#'   category-specific information coordinates. Best for downstream graphics
#'   and report drafts.}
#'   \item{[write_mfrm_residual_file()]}{Writes an observation-level residual
#'   file, optionally with modeled category probabilities. Best for external
#'   case review or reproducible handoff.}
#'   \item{[write_mfrm_subset_file()]}{Writes connected-subset summary and
#'   node-membership files. Best for scale-linking review outside R.}
#'   \item{[reporting_checklist()]}{Turns analysis status into an action list
#'   with priorities and next steps. Best for closing reporting gaps.}
#'   \item{[build_apa_outputs()]}{Creates manuscript-draft text, notes,
#'   captions, and section maps from a shared reporting contract.}
#'   \item{[build_summary_table_bundle()]}{Converts supported `summary()`
#'   outputs into named `data.frame` tables with a compact index for appendix
#'   or manuscript handoff, including recovery simulation and recovery
#'   assessment outputs. It also supports bundle-level `summary()` / `plot()`
#'   for QC before export.}
#'   \item{[export_summary_appendix()]}{Exports those documented summary-table
#'   bundles as CSV and optional HTML appendix artifacts without requiring the
#'   broader fit-based export bundle. This is the preferred export route for
#'   recovery simulation evidence.}
#'   \item{[apa_table()]}{Can now take those summary-table bundles directly,
#'   so a selected component can move from `summary()` to a formatted handoff
#'   table without rebuilding the analysis object path.}
#' }
#'
#' @section Read availability and method limits separately:
#' An available table or a passed numerical check does not establish that its
#' intervals have accurate repeated-sample coverage. For implemented experimental
#' methods, incomplete qualification means that a stated statistical property
#' remains unestablished; inspect the specific limitation before interpreting it.
#' Corrected JML retains this explanation in `interpretation_notes`, two-family
#' GPCM in `gpcm_interpretation_notes`, and attached population intervals in
#' `population_<name>_interpretation_notes`. These tables accompany the other
#' result tables in reports and CSV exports; keep them with the estimates.
#'
#' An `IntervalStatus` of `Incomplete` has a different meaning: only one endpoint
#' was computed. It is neither an infinite bound nor a qualified one-sided
#' interval. `Unavailable` means the requested result was not supplied; inspect
#' its reason. `CIEligible`, where present, concerns that interval procedure's
#' checks, not a universal coverage guarantee. Keep numerical status, method
#' limitations and unavailable results when saving or sharing an analysis.
#'
#' @section Practical interpretation rules:
#' - Use bundle summaries first, then drill down into component tables.
#' - Use [precision_review_report()] to determine whether formal inference is
#'   supported for the fitted result.
#' - Treat category and bias outputs as complementary layers rather than
#'   substitutes for overall fit review.
#' - Treat zero-count score categories as scale-functioning caveats. Boundary
#'   zero-count categories can be retained with explicit `rating_min` /
#'   `rating_max`; retaining an intermediate zero-count category with
#'   `keep_original = TRUE` creates an unsupported adjacent-step contrast in a
#'   polytomous fitted ladder, so new fits stop before optimization.
#'   `summary(describe_mfrm_data(...))` exposes these in `Notes`, printed
#'   `Caveats`, and `$caveats`; `summary(fit)` carries full structured caveats
#'   into printed `Caveats` and `$caveats`, with `Key warnings` as a short
#'   triage subset. Summary-table exports use `score_category_caveats` and
#'   `analysis_caveats`.
#' - Use [reporting_checklist()] before [build_apa_outputs()] when a report
#'   still needs missing diagnostics or clearer caveats.
#'
#' @section Additional task-specific workflows:
#' These extend the basic route for eligible ordinary fits; they are not
#' extra steps required for every report and do not apply to native CML.
#' - Run documentation:
#'   [fit_mfrm()] -> [specifications_report()] -> [data_quality_report()].
#' - Precision and facet review:
#'   [diagnose_mfrm()] -> [precision_review_report()] ->
#'   [facet_statistics_report()].
#' - Scale review:
#'   [rating_scale_table()] -> [category_structure_report()] ->
#'   [category_curves_report()].
#' - Manuscript handoff (`RSM` / `PCM`):
#'   [reporting_checklist()] -> [build_apa_outputs()] ->
#'   [build_summary_table_bundle()] -> `summary()` / `plot()` -> [apa_table()]
#'   or [export_summary_appendix()] /
#'   [export_mfrm_bundle()](include = "summary_tables").
#' - `GPCM` handoff:
#'   [reporting_checklist()] -> direct summaries/plots ->
#'   [build_apa_outputs()] or [build_summary_table_bundle()] ->
#'   [export_summary_appendix()] or caveated [export_mfrm_bundle()], with
#'   `gpcm_boundary` retained in report/export objects.
#' - Recovery simulation handoff:
#'   [simulate_mfrm_recovery()] -> [plot()] / [assess_mfrm_recovery()] ->
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
#' # Long-running workflow; run explicitly when needed.
#' # example("mfrmr_reports_and_tables", run.dontrun = TRUE)
#' \dontrun{
#' # Load the package and example ratings
#' library(mfrmr)
#' ratings <- load_mfrmr_data("example_operational")
#'
#' # Fit the model
#' fit <- fit_mfrm(
#'   data = ratings,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   rating_min = 1, rating_max = 4, category_policy = "preserve",
#'   method = "MML",
#'   model = "RSM"
#' )
#'
#' summary(fit)$decision # Review the fit before interpreting it
#'
#' # Build the comprehensive results and report
#' res <- mfrm_results(fit)
#' report <- mfrm_report(res)
#' summary(report, view = "reader")
#'
#' # Find report sections that need attention
#' report$first_screen[, c("Area", "Status", "MainIssue", "NextAction")]
#'
#' # Extract individual person, rater, and criterion estimates
#' estimates <- as.data.frame(fit)
#' head(subset(estimates, Facet == "Person")) # First six persons
#' subset(estimates, Facet == "Rater")       # All raters
#' subset(estimates, Facet == "Criterion")   # All criteria
#'
#' # Keep the full result, including methods and unavailable-output reasons
#' file <- tempfile(fileext = ".rds") # Choose a lasting path for your analysis
#' saveRDS(res, file)
#' saved <- readRDS(file)
#' summary(saved)$triage # No refitting or repeated diagnostic calculation
#' }
#' @name mfrmr_reports_and_tables
NULL

#' Choose an mfrmr output helper by user goal
#'
#' @description
#' Find a function for the next question in your analysis. Start with
#' `mfrmr_output_guide("beginner")` and read its `Question` and `MainFunction`
#' columns. With no argument the function returns the full specialist catalogue,
#' which includes tables, reports, reviews and exports. It does not analyze data.
#'
#' @param scope Which rows to return. `"all"` returns the full route catalogue.
#'   `"plots"` returns a separate purpose-based plot capability table (see below).
#'   `"public"` returns the canonical six-step route for most users;
#'   `"beginner"` returns the same compact route rather than combining every
#'   beginner-labelled specialist row.
#'   `"entry"` returns the recommended first-screen routes. `"viewer"` returns
#'   local-viewer routes built around `mfrm_results(include = ...)`. `"binary"`
#'   returns the two-category person-item Rasch route and checks. Other values
#'   filter to one output family or to `GPCM`-relevant routes.
#'   `"linking"` returns anchor, drift, equating and DIF sensitivity routes.
#'   `"calibration"` returns the portable fixed-calibration lifecycle and
#'   artifact-only scoring route, conditional decisions and fixed-plan
#'   classification accuracy/consistency.
#'   `"simulation"` and `"network"` return advanced design-review rows.
#'   `"response_time"` returns descriptive response-time QC rows.
#'   `"models"` compares fixed-facet, shared-rater and Person-specific testlet
#'   workflows, including their distinct prediction and reporting boundaries.
#'   `"feedback"` distinguishes fixed-rater uncertainty, unexpected rating
#'   patterns, shared-rater uncertainty, observed agreement, centrality,
#'   criterion-referenced accuracy and screening accuracy with known truth.
#'   `"features"`, `"imputation"` and `"gtheory"` show exploratory attributes,
#'   assigned-score multiple imputation and observed-score G/D-study workflows,
#'   including their own table, plot and saving routes.
#'   `"facets"`, `"conquest"`, and `"r"` return user-pathway rows for people
#'   arriving from those workflows.
#'
#' @details
#' Naming convention used by the guide:
#' - `*_table`: focused table or table-like result for one evidence source
#' - `*_report`: multi-table evidence bundle for a reporting question
#' - `*_review`: status, interpretation, or decision-support object
#' - `*_bundle`: reusable collection of tables/metadata for handoff
#' - `export_*`: writes files or appendix artifacts
#'
#' @section First-screen route:
#' Use `mfrmr_output_guide("public")` or
#' `mfrmr_output_guide("beginner")` for the shortest top-level API map:
#' an explicit [describe_mfrm_data()] check and [fit_mfrm()] MML fit, the lightweight
#' fit summary, `res <- mfrm_results(fit)` and its triage, optional figures,
#' and report/export or RDS saving. Plotting is not required before saving.
#' Earlier scripts may obtain the same result bundle from
#' `summary(fit, profile = "facets")$results`; no FACETS software is required.
#' Use `mfrmr_output_guide("entry")` when you
#' specifically need alternative first-screen creation routes, including
#' existing result objects, the optional viewer, or interactive console work.
#' After creating `res`, use `summary(res)$next_actions` to choose a
#' purpose-specific specialist helper.
#' Use `mfrmr_output_guide("viewer")` when the next step is the optional local
#' Shiny reader; it shows which `include` preset to use before calling
#' [launch_mfrmr_viewer()].
#' Use `mfrmr_output_guide("feedback")` when preparing rater feedback. Start
#' with the question and the fitted model: severity, response misfit and the
#' accuracy of a warning rule are different quantities. A severe rater need
#' not misfit, and an observed flag does not establish poor rater quality.
#' For observed agreement, choose [pairwise_kappa_table()],
#' [krippendorff_alpha_table()], [gwet_ac2_table()] or [pairwise_icc_table()]
#' according to the score scale, rater target and missingness assumptions.
#' [score_comparison_table()] compares a score source with reference ratings;
#' [rater_centrality_table()] describes residual and threshold-spread indices;
#' [fit_rater_accuracy()] models exact matches to fixed reference ratings.
#' High rater separation reliability indicates distinguishable severities,
#' not high agreement. These quantities do not replace one another.
#' For a fixed future plan, use [classification_table()] under
#' `scope = "calibration"`; its model-based accuracy/consistency differs from
#' observed reference agreement and posterior ability-band probabilities.
#' For DIF sensitivity, `scope = "linking"` points to
#' [mantel_haenszel_table()], [ordinal_logistic_table()] and the shared-roster
#' example in `vignette("mfrmr-linking-and-dff")`.
#' For an individual native additive RSM/PCM or experimental two-family GPCM sheet, use
#' `mfrm_report(res, style = "rater", facet = "Rater", rater = "R01")`.
#' The selected sheet omits source identifiers; the comprehensive result and
#' export bundle retain the original analysis.
#' Use `mfrmr_output_guide("psychometric")` for the technical table,
#' review, and reporting routes whose interpretation boundaries should be
#' checked before manuscript use.
#'
#' @section How to use this guide:
#' For route scopes, read `Question` first, then open the help for a function
#' in `MainFunction`. For `"plots"`, use the figure instructions below.
#' `UseWhen` describes its inputs and purpose; `NextStep` explains what to inspect.
#' Cells containing `...` are outlines, not complete scripts to paste and run.
#' Use [mfrmr_workflow_methods] for a runnable introduction and an explanation
#' of function names. Inspect `DecisionBoundary` before interpreting a result.
#' For construct-focused research, the *Use results as evidence about a
#' construct* section of [mfrmr_methods_reference] connects a claim to observable
#' patterns, competing explanations and the limits of a conclusion.
#' For `GPCM`, use `scope = "gpcm"` to find both the support matrix
#' and the table that explains how out-of-scope routes are handled.
#' Use `scope = "exports"` to distinguish saved-analysis archives from
#' [export_mfrm_model()], which generates a Stan/brms model for a separate
#' Bayesian analysis with explicit priors. That route returns to mfrmr through
#' [read_mfrm_posterior()]: inspect `summary(x)` before using `plot(x)`,
#' `confint(x)` or [mfrm_report()]. Other interval levels, contrasts and
#' predictions require the backend's joint draws.
#'
#' @section Check the model before following a route:
#' The six-step beginner route starts with ordinary RSM/PCM. The result class
#' `mfrm_fit` alone does not establish which diagnostics or figures are valid:
#' one-family GPCM, two-family GPCM and corrected JML share that class.
#' Read the fitted model and estimation method in `summary(fit)` first.
#' For GPCM, consult [gpcm_capability_matrix()] and
#' `vignette("mfrmr-gpcm-scope", package = "mfrmr")`.
#' Two-family GPCM and corrected JML have separate result/report and descriptive
#' response-diagnostic routes; ordinary fit tests and Wright/Pathway maps do
#' not apply. Corrected JML supplies point plots for slopes, locations and steps,
#' without structural confidence intervals.
#'
#' Use `scope = "models"` for fixed-facet, shared-rater and testlet roles,
#' `scope = "feedback"` for the distinct rater-review targets, and
#' `scope = "calibration"` for scoring or fixed-plan decisions.
#' For exploratory bias review across eligible facet pairs, use
#' [estimate_all_bias()], followed by `summary(screens)` and `plot(screens)`.
#' Use [estimate_bias()] for a prespecified pair. Retain failed pairs and
#' unavailable cells; the collection does not supply multiplicity-controlled
#' inference or evidence that no bias exists. Save the complete collection.
#'
#' Use `scope = "gtheory"` for observed-score assessment planning. Start with
#' the G-study component or covariance view, then compare D-study curves while
#' keeping other facet counts fixed or in separate panels. Use
#' `plot(d, type = "error_components")` to inspect each plan's error sources,
#' including signed covariance contributions, or `type = "sem"` for error in
#' score units. Error-source contributions are point estimates; SEM uncertainty
#' requires the separately supported interval calculation. These projections
#' and paired plan-difference intervals answer different questions. Component
#' SEs/intervals are available for the designs and methods documented in
#' [confint.mfrm_multivariate_gstudy()] and [vcov.mfrm_multivariate_gstudy()].
#' For coefficient sampling SEs, use [confint.mfrm_multivariate_d_study()]
#' with `method = "delta"` and `assumption = "normal"`. These differ from
#' measurement SEM. Its normal-reference intervals and plots retain bounds
#' outside `[0,1]` with approximation qualifications; defaults are unchanged.
#' General distribution-robust comparison intervals are not supplied.
#' Save each complete result with [base::saveRDS()] and reopen it with
#' [base::readRDS()]; a CSV alone does not retain its design or qualifications.
#' Fixed-facet output also depends on the estimator: native CML uses dedicated
#' summaries and [mfrm_facet_intervals()] plots/tables, without ordinary Person
#' scoring or consolidated reports. [mfrmr_interval_guide()] separates that
#' route from MML interaction contrasts and their population-SD standardization.
#'
#' @section Choosing a figure:
#' Use `mfrmr_output_guide("plots")` to find selected plots by purpose, input
#' class and explicit plot call. This is a curated map, not an exhaustive list
#' of every plot method, view or component. Unlisted routes are not necessarily
#' unsupported. Existing scopes, including `"all"`, keep their route-table format.
#'
#' For linked figure previews and runnable examples, open
#' `vignette("mfrmr-visual-diagnostics")`.
#' `ResultFunction` names a help page for creating the required result.
#' In `PlotCall`, replace `x` with the indicated result object and save the
#' returned object as `p`. Calls use `draw = FALSE`; change it to `TRUE` to
#' display the original figure. Then follow `GGPlotCall`, if available.
#' `DataComponent` is for [plot_data()] extraction; it is not automatically a
#' valid or equivalent `component` argument to [as_ggplot()].
#'
#' `GGPlot` distinguishes `"dedicated"` (a converter for this plot), `"native"`
#' (the plot method already returns ggplot), `"generic"` (a table-based view
#' that need not reproduce the original figure), and `"unavailable"` (use the
#' original plot or build a custom graphic from its data). ggplot routes need
#' the optional ggplot2 package. Display support does not broaden the fitted
#' model's statistical scope or establish interval coverage. Read `Notes` and
#' the source function's help before interpretation.
#'
#' @return For `scope = "plots"`, a data.frame with `Question`, `InputClass`,
#'   `PlotCall`, `PlotName` (the name in the saved plot data), `DataComponent`,
#'   `GGPlot`, `GGPlotCall` (`NA` when unavailable), `Notes`, `ResultFunction`
#'   and `NextStep`.
#'   All other scopes return a data.frame with one row per route and columns:
#' - `Scope`
#' - `Question`
#' - `OutputFamily`
#' - `Lifecycle`
#' - `UserLevel`
#' - `APILayer`
#' - `ObjectRole`
#' - `DecisionBoundary`
#' - `RecommendedEntry`
#' - `MainFunction`
#' - `UseWhen`
#' - `TypicalInput`
#' - `NextStep`
#' - `GPCMStatus`
#' - `Notes`
#'
#' @examples
#' beginner <- mfrmr_output_guide("beginner")
#' beginner[, c("Question", "MainFunction", "NextStep")]
#'
#' # Ask for a specialist map only when that question arises.
#' linking <- mfrmr_output_guide("linking")
#' linking[, c("Question", "MainFunction", "UseWhen")]
#'
#' feedback <- mfrmr_output_guide("feedback")
#' feedback[, c("Question", "MainFunction", "DecisionBoundary")]
#'
#' figures <- mfrmr_output_guide("plots")
#' figures[, c("Question", "ResultFunction", "GGPlot")]
#' @concept reporting workflow
#' @concept route selection
#' @concept GPCM boundaries
#' @export
mfrmr_output_guide <- function(scope = c("all", "public", "beginner", "psychometric",
                                         "entry", "viewer", "binary", "tables", "reports", "reviews",
                                         "bundles", "exports", "compatibility",
                                         "gpcm", "calibration", "simulation", "linking", "network",
                                         "response_time", "facets", "conquest", "r", "models",
                                         "features", "imputation", "gtheory", "feedback", "plots")) {
  scope <- match.arg(scope)
  if (identical(scope, "plots")) return(.mfrmr_plot_guide())

  out <- .mfrmr_table_mfrmr_output_guide_out

  public_rows <- .mfrmr_table_mfrmr_output_guide_public_rows

  calibration_rows <- .mfrmr_table_mfrmr_output_guide_calibration_rows

  entry_rows <- .mfrmr_table_mfrmr_output_guide_entry_rows

  viewer_rows <- .mfrmr_table_mfrmr_output_guide_viewer_rows

  binary_rows <- .mfrmr_table_mfrmr_output_guide_binary_rows

  linking_rows <- .mfrmr_table_mfrmr_output_guide_linking_rows

  advanced_rows <- .mfrmr_table_mfrmr_output_guide_advanced_rows

  gpcm_rows <- .mfrmr_table_mfrmr_output_guide_gpcm_rows

  migration_entry_rows <- .mfrmr_table_mfrmr_output_guide_migration_entry_rows

  user_pathway_rows <- .mfrmr_table_mfrmr_output_guide_user_pathway_rows

  model_rows <- .mfrmr_table_mfrmr_output_guide_model_rows

  workflow_rows <- .mfrmr_table_mfrmr_output_guide_workflow_rows

  feedback_rows <- .mfrmr_table_mfrmr_output_guide_feedback_rows

  feedback_rows <- rbind(feedback_rows, .mfrmr_table_mfrmr_output_guide_data)

  calibration_rows <- rbind(calibration_rows, .mfrmr_table_mfrmr_output_guide_data_1)

  rating_evaluation_rows <- .mfrmr_table_mfrmr_output_guide_rating_evaluation_rows

  model_export_row <- .mfrmr_table_mfrmr_output_guide_model_export_row

  workflow_rows <- rbind(workflow_rows, .mfrmr_table_mfrmr_output_guide_data_2)

  network_rows <- .mfrmr_table_mfrmr_output_guide_network_rows

  out <- rbind(
    public_rows, calibration_rows, entry_rows, viewer_rows, binary_rows, out, linking_rows,
    advanced_rows, gpcm_rows, migration_entry_rows, user_pathway_rows, model_rows, workflow_rows,
    feedback_rows, rating_evaluation_rows, model_export_row, network_rows
  )
  out$Lifecycle <- "stable"
  out$Lifecycle[out$OutputFamily %in% "compatibility" | out$Scope %in% c("compatibility", "facets", "conquest")] <- "compatibility"
  out$Lifecycle[out$Scope %in% c("simulation", "network")] <- "advanced"
  out$UserLevel <- "intermediate"
  out$UserLevel[out$Scope %in% c("public", "calibration", "entry", "viewer", "binary")] <- "beginner"
  out$UserLevel[out$Scope %in% c("simulation", "network", "reviews")] <- "advanced"
  out$UserLevel[out$Scope %in% c("compatibility", "facets", "conquest")] <- "migration"
  out$APILayer <- "specialist_component"
  out$APILayer[out$Scope %in% "public"] <- "top_level_public_surface"
  out$APILayer[out$Scope %in% c("calibration", "entry", "viewer", "binary")] <- "recommended_entry_route"
  out$APILayer[out$Scope %in% c("reports", "reviews", "tables", "bundles", "exports", "linking", "gpcm")] <- "specialist_followup"
  out$APILayer[out$Scope %in% c("simulation", "network")] <- "advanced_design_review"
  out$APILayer[out$Scope %in% c("compatibility", "facets", "conquest", "r")] <- "migration_or_integration"
  out$ObjectRole <- "specialist evidence component"
  out$ObjectRole[out$OutputFamily %in% "table"] <- "focused evidence table"
  out$ObjectRole[out$OutputFamily %in% "report"] <- "reporting evidence bundle"
  out$ObjectRole[out$OutputFamily %in% "review"] <- "review and follow-up surface"
  out$ObjectRole[out$OutputFamily %in% "bundle"] <- "handoff table bundle"
  out$ObjectRole[out$OutputFamily %in% "export"] <- "file export surface"
  out$ObjectRole[out$OutputFamily %in% "viewer"] <- "local reader over existing results"
  out$ObjectRole[out$OutputFamily %in% "guide"] <- "route-selection guide"
  out$ObjectRole[out$OutputFamily %in% "compatibility"] <- "compatibility presentation contract"
  out$ObjectRole[grepl("mfrm_results(", out$MainFunction, fixed = TRUE)] <- "comprehensive result object"
  out$ObjectRole[grepl("fit_mfrm", out$MainFunction, fixed = TRUE)] <- "model estimation and result-object entry"
  out$ObjectRole[grepl("mfrm_report", out$MainFunction, fixed = TRUE)] <- "report-readiness surface"
  out$ObjectRole[grepl("launch_mfrmr_viewer", out$MainFunction, fixed = TRUE)] <- "local reader over existing results"
  out$ObjectRole[grepl("export_", out$MainFunction, fixed = TRUE)] <- "file export surface"
  out$ObjectRole[grepl("mfrmr_output_guide", out$MainFunction, fixed = TRUE)] <- "route-selection guide"
  out$ObjectRole[grepl("view = 'unavailable'", out$MainFunction, fixed = TRUE)] <- "unavailable-route alternatives"
  out$ObjectRole[grepl("mfrm_calibration_capabilities", out$MainFunction, fixed = TRUE)] <- "route-selection guide"
  out$ObjectRole[grepl("extract_mfrm_calibration", out$MainFunction, fixed = TRUE)] <- "portable calibration lifecycle"
  out$ObjectRole[grepl("save_mfrm_calibration", out$MainFunction, fixed = TRUE)] <- "calibration persistence"
  out$ObjectRole[grepl("score_mfrm_calibration", out$MainFunction, fixed = TRUE)] <- "artifact-only posterior scoring"
  out$ObjectRole[grepl("summary(scores)", out$MainFunction, fixed = TRUE)] <- "portable score-batch review"
  out$ObjectRole[grepl("with cut_scores", out$MainFunction, fixed = TRUE)] <- "conditional ability decisions"
  out$ObjectRole[grepl("mfrm_results_interactive", out$MainFunction, fixed = TRUE)] <- "explicit opt-in interactive entry"
  out$ObjectRole[grepl('summary(fit, profile = "fit"', out$MainFunction, fixed = TRUE)] <- "fit summary surface"
  out$ObjectRole[grepl('summary(fit, profile = "facets"', out$MainFunction, fixed = TRUE)] <- "comprehensive result object"

  out$DecisionBoundary <- "Specialist follow-up: inspect the source object and help page before treating output as report evidence."
  out$DecisionBoundary[out$ObjectRole == "conditional ability decisions"] <-
    "Posterior event probabilities conditional on point calibration, scoring prior and supplied cuts; no automatic decision rule or calibration uncertainty."
  out$DecisionBoundary[out$ObjectRole %in% "model estimation and result-object entry"] <-
    "Fits model parameters first; diagnostics, reporting, and validity claims are separate follow-up steps."
  out$DecisionBoundary[out$ObjectRole %in% "fit summary surface"] <-
    "Summarizes fit status and settings without computing comprehensive diagnostics or making a reporting decision."
  out$DecisionBoundary[out$ObjectRole %in% "comprehensive result object"] <-
    "Collects existing fit, diagnostic, table, reporting, and plot surfaces; it is not a new estimator or pass/fail decision."
  out$DecisionBoundary[out$ObjectRole %in% "report-readiness surface"] <-
    "Summarizes report readiness and wording routes; it does not recompute diagnostics or create an acceptance rule."
  out$DecisionBoundary[out$ObjectRole %in% "local reader over existing results"] <-
    "Displays an existing mfrm_results object; it does not estimate models, load external web apps, or alter diagnostics."
  out$DecisionBoundary[out$ObjectRole %in% "file export surface"] <-
    "Writes already-created evidence to disk; it does not recompute or strengthen evidence."
  out$DecisionBoundary[out$ObjectRole %in% "route-selection guide"] <-
    "Chooses the next help route; it is not analysis evidence."
  out$DecisionBoundary[out$ObjectRole %in% "out-of-scope route-status table"] <-
    "Explains supported-with-caveat, blocked, and deferred GPCM route handling; it does not broaden any route beyond its current capability row."
  out$DecisionBoundary[out$ObjectRole %in% "explicit opt-in interactive entry"] <-
    "Collects column choices interactively; move replay code into an explicit script before reporting."
  out$DecisionBoundary[out$ObjectRole %in% "portable calibration lifecycle"] <-
    "Creates fixed-standard-normal RSM/PCM MML, conditional estimated-normal one-family or experimental fixed-standard-normal two-family GPCM MML, or reference-prior RSM/PCM and shared-owner GPCM JML artifacts within their documented scope; incomplete GPCM global audits remain recorded. Every refusal must be resolved before freezing."
  out$DecisionBoundary[out$ObjectRole %in% "calibration persistence"] <-
    "Preserves and validates a calibration lifecycle object; loading is not authentication of an untrusted file."
  out$DecisionBoundary[out$ObjectRole %in% "artifact-only posterior scoring"] <-
    paste(
      "Returns posterior EAP uncertainty conditional on the frozen point",
      "calibration and recorded prior, excludes calibration-parameter",
      "uncertainty, and performs no refit."
    )
  out$DecisionBoundary[out$ObjectRole %in% "portable score-batch review"] <-
    paste(
      "Reviews score dispositions and conditional uncertainty; it does not",
      "establish source-fit quality, reliability, transport, or validity."
    )
  out$DecisionBoundary[grepl("precision", out$MainFunction, ignore.case = TRUE)] <-
    "Precision and separation evidence are not inter-rater agreement or standalone validity proof."
  out$DecisionBoundary[grepl("response_time", out$MainFunction, fixed = TRUE)] <-
    "Response-time outputs are descriptive QC context; they do not alter MFRM estimates, fit a joint speed-accuracy model, or create automatic exclusion rules."
  out$DecisionBoundary[grepl("bias|DFF|dff", out$MainFunction, ignore.case = TRUE)] <-
    "Bias and DFF rows are screening prompts; they are not final fairness conclusions without follow-up design and substantive review."
  out$DecisionBoundary[grepl("anchor|drift|linking|equating", out$MainFunction, ignore.case = TRUE)] <-
    "Anchor and linking evidence support scale-maintenance review; drift and equating claims require explicit multi-fit wave or form designs."
  out$DecisionBoundary[out$Scope %in% c("simulation", "network")] <-
    "Design, network, or operating-characteristic evidence informs planning and review; it is not an observed-data validity decision by itself."
  out$DecisionBoundary[out$Scope %in% c("compatibility", "facets", "conquest")] <-
    "Compatibility rows are presentation or migration contracts, not numerical equivalence claims unless external outputs are explicitly compared."

  out$RecommendedEntry <- out$Scope %in% c("public", "calibration", "entry", "viewer", "binary")
  dedicated <- out$Scope %in% c("features", "imputation", "gtheory")
  rating_evaluation <- out$MainFunction %in% rating_evaluation_rows$MainFunction
  targeted <- dedicated | rating_evaluation | out$Scope %in% c("models", "feedback")
  out$DecisionBoundary[targeted] <- out$Notes[targeted]
  batch_bias <- grepl("estimate_all_bias()", out$MainFunction, fixed = TRUE)
  out$DecisionBoundary[batch_bias] <- out$Notes[batch_bias]
  model_export <- out$MainFunction == model_export_row$MainFunction
  out$ObjectRole[model_export] <- "external model and saved posterior summary workflow"
  out$DecisionBoundary[model_export] <- out$Notes[model_export]
  out$ObjectRole[dedicated] <- "dedicated analysis and saved-output workflow"
  out$ObjectRole[out$Scope == "feedback"] <- "rater-feedback route selection"
  out$APILayer[out$Scope == "feedback"] <- "specialist_followup"
  out$Lifecycle[dedicated] <- "advanced"
  out$Lifecycle[out$Scope == "models" & out$MainFunction != "fit_mfrm()"] <- "advanced"
  out <- out[, c(
    "Scope", "Question", "OutputFamily", "Lifecycle", "UserLevel",
    "APILayer", "ObjectRole", "DecisionBoundary", "RecommendedEntry",
    "MainFunction", "UseWhen", "TypicalInput", "NextStep", "GPCMStatus",
    "Notes"
  ), drop = FALSE]

  if (identical(scope, "all")) {
    return(out)
  }
  if (identical(scope, "beginner")) {
    return(out[out$Scope %in% "public", , drop = FALSE])
  }
  if (identical(scope, "psychometric")) {
    keep <- out$Scope %in% c("tables", "reports", "reviews", "linking") &
      out$ObjectRole %in% c(
        "focused evidence table", "reporting evidence bundle",
        "review and follow-up surface", "specialist evidence component"
      )
    return(out[keep | rating_evaluation, , drop = FALSE])
  }
  if (identical(scope, "gpcm")) {
    keep <- !out$GPCMStatus %in% c("supported", "not_applicable")
    return(out[keep, , drop = FALSE])
  }
  if (identical(scope, "response_time")) {
    haystack <- paste(out$Question, out$MainFunction, out$UseWhen,
                      out$NextStep, out$Notes)
    keep <- grepl("response[-_ ]time|timing", haystack, ignore.case = TRUE)
    return(out[keep, , drop = FALSE])
  }
  out[out$Scope == scope, , drop = FALSE]
}
