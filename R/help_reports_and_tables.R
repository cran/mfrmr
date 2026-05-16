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
#'   Use [estimation_iteration_report()] and [precision_review_report()].
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
#' - "Did a simulation recover the known generating parameters well enough?"
#'   Use [evaluate_mfrm_recovery()] for the recovery study,
#'   [assess_mfrm_recovery()] for the adequacy checklist, and then
#'   [build_summary_table_bundle()] or [export_summary_appendix()] for the
#'   appendix handoff.
#'
#' @section Recommended report route:
#' 1. Start with [specifications_report()] and [data_quality_report()] to
#'    document the run and confirm usable data.
#' 2. Continue with [estimation_iteration_report()] and
#'    [precision_review_report()] to judge convergence and inferential strength.
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
#' If you are unsure which helper to call, start with
#' [mfrmr_output_guide()]. It returns a compact purpose-to-helper table that
#' separates `*_table`, `*_report`, `*_review`, `*_bundle`, `export_*`, and
#' compatibility routes.
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
#'   \item{[export_summary_appendix()]}{Exports those validated summary-table
#'   bundles as CSV and optional HTML appendix artifacts without requiring the
#'   broader fit-based export bundle. This is the preferred export route for
#'   recovery simulation evidence.}
#'   \item{[apa_table()]}{Can now take those summary-table bundles directly,
#'   so a selected component can move from `summary()` to a formatted handoff
#'   table without rebuilding the analysis object path.}
#' }
#'
#' @section Practical interpretation rules:
#' - Use bundle summaries first, then drill down into component tables.
#' - Treat [precision_review_report()] as the gatekeeper for formal inference.
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
#' - Bounded `GPCM` handoff:
#'   [reporting_checklist()] -> direct summaries/plots ->
#'   [build_summary_table_bundle()] -> [export_summary_appendix()].
#' - Recovery simulation handoff:
#'   [evaluate_mfrm_recovery()] -> [plot()] / [assess_mfrm_recovery()] ->
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
#'   quad_points = 7,
#'   maxit = 30
#' )
#' diag <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "both")
#'
#' spec <- specifications_report(fit)
#' summary(spec)$overview
#'
#' prec <- precision_review_report(fit, diagnostics = diag)
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

#' Choose an mfrmr output helper by user goal
#'
#' @description
#' `mfrmr_output_guide()` returns a compact table for choosing among the main
#' table, report, review, bundle, export, and compatibility helpers. It is a
#' user-facing map, not an analysis result.
#'
#' @param scope Which rows to return. `"all"` returns the full guide.
#'   Other values filter to one output family or to bounded-`GPCM`-relevant
#'   routes. `"facets"`, `"conquest"`, and `"r"` return user-pathway rows for
#'   people arriving from those workflows.
#'
#' @details
#' Naming convention used by the guide:
#' - `*_table`: focused table or table-like result for one evidence source
#' - `*_report`: multi-table evidence bundle for a reporting question
#' - `*_review`: status, interpretation, or decision-support object
#' - `*_bundle`: reusable collection of tables/metadata for handoff
#' - `export_*`: writes files or appendix artifacts
#'
#' @return A data.frame with one row per recommended route and columns:
#' - `Scope`
#' - `Question`
#' - `OutputFamily`
#' - `MainFunction`
#' - `UseWhen`
#' - `TypicalInput`
#' - `NextStep`
#' - `GPCMStatus`
#' - `Notes`
#'
#' @examples
#' mfrmr_output_guide()
#' mfrmr_output_guide("reviews")
#' mfrmr_output_guide("gpcm")
#' mfrmr_output_guide("facets")
#' mfrmr_output_guide("conquest")
#' mfrmr_output_guide("r")
#' @export
mfrmr_output_guide <- function(scope = c("all", "tables", "reports", "reviews",
                                         "bundles", "exports", "compatibility",
                                         "gpcm", "facets", "conquest", "r")) {
  scope <- match.arg(scope)

  out <- data.frame(
    Scope = c(
      "reports", "reports", "reviews", "reports", "tables", "reviews",
      "reviews", "reviews", "bundles", "bundles", "exports", "compatibility"
    ),
    Question = c(
      "Document the model setup and run settings",
      "Check whether data were filtered, dropped, or remapped",
      "Decide how strongly precision claims can be phrased",
      "Summarize facet variability, separation, and measurability",
      "Review category functioning and expected-score curves",
      "Screen bias, DFF, or interaction evidence",
      "Review anchors, drift, and linking readiness",
      "Compare equal-weighting and bounded-GPCM routes",
      "Turn summaries into reusable appendix tables",
      "Assemble manuscript-oriented narrative output",
      "Write files for appendix, replay, or handoff",
      "Serve a legacy-compatible downstream layout"
    ),
    OutputFamily = c(
      "report", "report", "review", "report", "table", "review",
      "review", "review", "bundle", "bundle", "export", "compatibility"
    ),
    MainFunction = c(
      "specifications_report()",
      "data_quality_report()",
      "precision_review_report()",
      "facet_statistics_report()",
      "rating_scale_table(); category_structure_report(); category_curves_report()",
      "estimate_bias(); analyze_dff(); bias_interaction_report()",
      "review_mfrm_anchors(); detect_anchor_drift(); build_linking_review()",
      "build_model_choice_review(); build_weighting_review(); compare_mfrm()",
      "build_summary_table_bundle()",
      "reporting_checklist(); build_apa_outputs()",
      "export_summary_appendix(); export_mfrm_bundle(); build_mfrm_manifest()",
      "run_mfrm_facets(); facets_output_file_bundle(); facets_output_contract_review()"
    ),
    UseWhen = c(
      "You need methods-section settings or reproducibility context.",
      "You need to explain retained rows, missing values, or unknown elements.",
      "You need to separate model-based, hybrid, and exploratory precision evidence.",
      "You need facet-level reporting before writing interpretation.",
      "You need category diagnostics or curve data for figures.",
      "You need screening evidence for follow-up fairness or interaction review.",
      "You need operational scale-maintenance checks for RSM/PCM workflows.",
      "You need to review whether discrimination-based reweighting changes conclusions.",
      "You need named data.frame outputs for appendix or handoff.",
      "You need manuscript text, captions, notes, and section maps.",
      "You need files on disk rather than R objects.",
      "You need a fixed-width, graphfile, scorefile, or contract-check layout."
    ),
    TypicalInput = c(
      "mfrm_fit",
      "mfrm_fit plus optional raw data",
      "mfrm_fit plus mfrm_diagnostics",
      "mfrm_fit plus mfrm_diagnostics",
      "mfrm_fit",
      "mfrm_fit plus optional diagnostics/grouping",
      "anchor tables, drift results, or equating chain objects",
      "two or more fitted candidate models",
      "summary(), review, recovery, prediction, or checklist objects",
      "mfrm_fit plus diagnostics/checklist evidence",
      "summary-table bundle, manifest, or fitted analysis bundle",
      "mfrm_fit or run_mfrm_facets() output"
    ),
    NextStep = c(
      "Read summary rows before adding interpretation.",
      "Inspect row-level reasons before reporting exclusions.",
      "Use the status mix to decide wording strength.",
      "Use plots or summary-table bundles for presentation.",
      "Use plot helpers or export curve tables.",
      "Treat positive screens as follow-up prompts, not final decisions.",
      "Use linking synthesis for RSM/PCM only; keep GPCM anchor checks exploratory.",
      "Report GPCM fit gains as sensitivity evidence, not automatic scoring policy.",
      "Use apa_table() or export_summary_appendix().",
      "Use build_summary_table_bundle() for appendix tables.",
      "Check written_files and available_outputs before sharing.",
      "Prefer package-native routes unless a downstream tool requires this layout."
    ),
    GPCMStatus = c(
      "supported",
      "supported",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "exploratory_only",
      "supported_with_caveat",
      "supported_for_direct_outputs",
      "blocked_for_build_apa_outputs",
      "summary_appendix_supported; fit_bundle_blocked",
      "graph_only_or_blocked_by_score_side_semantics"
    ),
    Notes = c(
      "Best first report object when methods documentation is the goal.",
      "Pairs naturally with describe_mfrm_data().",
      "Use before formal confidence or reliability language.",
      "Good bridge from diagnostics to report text.",
      "Category curves carry reusable plot data.",
      "Bias outputs are screening layers.",
      "Formal linking synthesis remains RSM/PCM-scoped.",
      "Use with gpcm_capability_matrix().",
      "Works for several object classes; inspect table_index.",
      "RSM/PCM manuscript route; GPCM stays on checklist/direct outputs.",
      "Choose the narrowest export route that satisfies the handoff.",
      "Compatibility is a presentation contract, not numerical equivalence."
    ),
    stringsAsFactors = FALSE
  )

  user_pathway_rows <- data.frame(
    Scope = c(
      "facets", "facets", "facets", "facets", "facets", "facets", "facets",
      "facets", "facets", "facets", "facets",
      "conquest", "conquest", "conquest",
      "r", "r", "r"
    ),
    Question = c(
      "State the FACETS relationship before using FACETS-style routes",
      "Translate FACETS direct and group anchor blocks",
      "Review anchor drift across forms, raters, or waves",
      "List fit measures and misfit flags by facet",
      "Explain FACETS df and ZSTD conversion",
      "Bring an external FACETS fit table into the review",
      "Review rating-scale categories, fair averages, and expected curves",
      "Review FACETS Table 14-style bias and interaction signals",
      "Draw a Wright map / variable map on the common logit scale",
      "Check score support and rater response patterns before fitting claims",
      "Write residual and subset files for external review",
      "Prepare a scoped ConQuest overlap case",
      "Compare extracted ConQuest tables after the external run",
      "State where mfrmr is less free than ConQuest",
      "Reuse plotting data for custom graphics",
      "Start from data-quality plot data for QC dashboards",
      "Combine tables and plot data for reports"
    ),
    OutputFamily = c(
      "compatibility", "review", "review", "table", "review", "review",
      "table", "review", "plot-data", "report",
      "export",
      "compatibility", "review", "compatibility",
      "plot-data", "plot-data", "bundle"
    ),
    MainFunction = c(
      "facets_positioning_guide(); facets_feature_coverage(); run_mfrm_facets(); mfrmRFacets(); facets_output_file_bundle()",
      "review_mfrm_anchors(); make_anchor_table(); fit_mfrm(anchors = ..., group_anchors = ...)",
      "anchor_to_baseline(); detect_anchor_drift(); build_equating_chain(); plot_anchor_drift()",
      "fit_measures_table(); facets_chisq_table(); displacement_table()",
      "facets_fit_df_guide(); diagnose_mfrm(fit_df_method = \"both\")",
      "read_facets_fit_table(); facets_fit_review(); plot(..., type = \"df_sensitivity\")",
      "rating_scale_table(); category_structure_report(); category_curves_report(); fair_average_table(); plot_fair_average()",
      "estimate_bias(); bias_interaction_report(); bias_pairwise_report(); plot_bias_interaction()",
      "plot(fit, type = \"wright\"); plot_wright_unified(); plot_data(type = \"wright\")",
      "data_quality_report(); plot(..., type = \"dashboard\")",
      "write_mfrm_residual_file(); write_mfrm_subset_file(); facets_output_file_bundle()",
      "build_conquest_overlap_bundle()",
      "normalize_conquest_overlap_files(); review_conquest_overlap()",
      "reporting_checklist(); reference_case_benchmark()",
      "plot_data_components(); plot_data(); plot(..., draw = FALSE); plot(fit, type = \"pathway\", draw = FALSE); plot_data(category_curves_report(...), component = \"plot_long\"); plot_bias_interaction(..., draw = FALSE); plot_information(..., draw = FALSE)",
      "plot_data(data_quality_report(...), type = \"dashboard\")",
      "build_summary_table_bundle(); build_visual_summaries()"
    ),
    UseWhen = c(
      "You want to prevent FACETS-complete-reproduction wording while still using familiar names, coverage maps, fixed files, or graphfile-style handoff.",
      "You are translating FACETS D/A or D/G anchor blocks before fitting.",
      "You need to screen common-element drift across administrations, forms, raters, or waves.",
      "You want a Table 7-like view of measures, standard errors, and fit status.",
      "You need to show whether a ZSTD difference is caused by df convention rather than MnSq.",
      "You have already extracted a FACETS fit table and want a df/ZSTD-sensitive side-by-side review.",
      "You want FACETS Table 8/category and Table 12 fair-average style evidence.",
      "You want a FACETS Table 14-like route for local interaction-bias screening.",
      "You want a variable-map view of persons, facet levels, and step thresholds.",
      "You want a pre-interpretation QC screen for zero categories, dominant raters, or dropped rows.",
      "You need standalone residual rows or connected-subset membership files for a reviewer, spreadsheet, or external QC workflow.",
      "You need the narrow documented MML latent-regression overlap bundle.",
      "You have external ConQuest output tables already extracted to rectangular files.",
      "You need to explain why a ConQuest command-file workflow is not interchangeable with this package.",
      "You want to draw your own ggplot2/base/plotly figure from package plot data.",
      "You want the dashboard tables without accepting the package's default base-R drawing.",
      "You need stable data frames and reusable plot-data objects for an R Markdown or Quarto report."
    ),
    TypicalInput = c(
      "mfrm_fit or run_mfrm_facets() output",
      "raw data plus anchor and/or group-anchor tables",
      "baseline anchors, baseline fit, new data, or a list of fitted waves",
      "mfrm_fit plus optional diagnostics",
      "mfrm_fit, mfrm_diagnostics, or a fit-measures/review table",
      "mfrm_fit plus extracted FACETS table",
      "mfrm_fit plus optional diagnostics",
      "mfrm_fit plus diagnostics",
      "mfrm_fit",
      "mfrm_fit plus optional raw data",
      "mfrm_fit plus optional diagnostics",
      "eligible MML RSM/PCM fit",
      "ConQuest population, item, and case tables plus the overlap bundle",
      "mfrm_fit plus optional diagnostics or benchmark request",
      "mfrm_plot_data or any mfrmr object with a plot method",
      "mfrm_data_quality bundle",
      "summary/checklist/report/review objects"
    ),
    NextStep = c(
      "Start with facets_positioning_guide(), then use native mfrmr objects for inference and FACETS-style files for handoff layout.",
      "Review issue_counts and recommendations before choosing anchor_policy for fit_mfrm().",
      "Check drift flags and link residuals before treating linked measures as comparable.",
      "Use threshold profiles to show how underfit/overfit rates depend on the chosen band.",
      "Compare MnSq first, then df, then ZSTD; label convention-sensitive flags explicitly.",
      "Normalize the external table first, then inspect df_sensitivity, df_sensitive, and external_comparison.",
      "Check zero/sparse categories and threshold ordering before reporting fair averages.",
      "Use estimate_bias() for screening, bias_interaction_report() for ranked cells, and bias_pairwise_report() for pairwise contrasts.",
      "Use draw = FALSE or plot_data(type = \"wright\") when you need reusable coordinates for custom graphics.",
      "Resolve severe QC flags before interpreting fit or fairness outputs.",
      "Keep residual, subset, and graph/score files distinct; use the native writers unless a legacy graph/score contract is required.",
      "Run ConQuest externally and keep the generated command text as a conservative template.",
      "Review missing/non-numeric/status rows before treating differences as numerical evidence.",
      "Document that mfrmr currently exposes selected overlap checks, not the full ConQuest modeling language.",
      "Start with plot_data_components(); for pathway maps use pathway_long/pathway_annotations/fit_measures; for category, bias, and information plots use plot_long plus annotations/settings.",
      "Use component names such as quality_flags, category_usage_by_facet, and facet_response_patterns.",
      "Keep table objects and plot-data objects separate so custom reporting remains reproducible."
    ),
    GPCMStatus = c(
      "graph_only_or_blocked_by_score_side_semantics",
      "supported_with_caveat",
      "exploratory_for_gpcm; rsm_pcm_linking_route",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported",
      "supported",
      "blocked_for_gpcm; rsm_pcm_overlap_only",
      "blocked_for_gpcm; rsm_pcm_overlap_only",
      "blocked_for_gpcm; rsm_pcm_overlap_only",
      "supported",
      "supported",
      "supported_for_direct_outputs"
    ),
    Notes = c(
      "mfrmr is not a FACETS numerical clone; familiar names help transition, but estimates remain package-native unless external output is supplied.",
      "Direct anchors fix element logits; group anchors constrain a group mean, with direct anchors taking precedence.",
      "Formal linking synthesis is RSM/PCM-scoped; bounded GPCM drift review remains exploratory.",
      "Closest current route for FACETS users who expect fit measures in one table.",
      "The guide explains engine df, FACETS-style df, Wilson-Hilferty ZSTD, and WHEXACT caveats.",
      "This is a review of supplied rectangular output, not raw FACETS text parsing; use df_sensitive for convention-sensitive rows.",
      "GPCM fair averages are slope-aware direct outputs, not FACETS score-side equivalence.",
      "Bias outputs are conditional screening layers; use substantive review before fairness conclusions.",
      "Wright maps visualize fitted scale locations; under bounded GPCM, interpret step/threshold locations with slope-aware caveats.",
      "Designed to catch user-visible data problems before fit interpretation.",
      "Residual and subset writers are package-native CSV/TSV handoff routes, not exact FACETS fixed-field command-file clones.",
      "ConQuest support is intentionally scoped to a documented comparison case.",
      "The helper reads extracted tables, not raw ConQuest report text.",
      "Use this row when reporting scope differences to ConQuest users.",
      "Best R-first route for users who want full control over graphics, fit-annotated pathway maps, bias-screening views, and TIF/conditional-SEM plots.",
      "The accessor avoids depending on list internals in user scripts.",
      "Use this path for Quarto/R Markdown pipelines and downstream custom tables."
    ),
    stringsAsFactors = FALSE
  )

  out <- rbind(out, user_pathway_rows)

  if (identical(scope, "all")) {
    return(out)
  }
  if (identical(scope, "gpcm")) {
    keep <- out$GPCMStatus != "supported" |
      grepl("GPCM|gpcm", out$Question, ignore.case = TRUE) |
      grepl("GPCM|gpcm", out$Notes, ignore.case = TRUE)
    return(out[keep, , drop = FALSE])
  }
  out[out$Scope == scope, , drop = FALSE]
}
