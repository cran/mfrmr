#' mfrmr Workflow and Method Map
#'
#' @description
#' Quick reference for end-to-end `mfrmr` analysis and for checking which
#' output objects support `summary()` and `plot()`.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#'    For final reporting, prefer `method = "MML"` unless you explicitly want
#'    a fast exploratory JML pass.
#' 2. (Optional) Use [run_mfrm_facets()] or [mfrmRFacets()] for a
#'    legacy-compatible one-shot workflow wrapper.
#' 3. Build diagnostics with [diagnose_mfrm()].
#' 4. (Optional) Estimate interaction bias with [estimate_bias()].
#' 5. Generate reporting bundles:
#'    [apa_table()], [build_fixed_reports()], [build_visual_summaries()].
#' 6. (Optional) Audit report completeness with [reference_case_audit()].
#'    Use `facets_parity_report()` only when you explicitly need the
#'    compatibility layer.
#' 7. (Optional) Benchmark packaged reference cases with
#'    [reference_case_benchmark()] when you want an internal package-native
#'    benchmark/audit run.
#' 8. (Optional) For design planning or future scoring, move to the
#'    simulation/prediction layer:
#'    [build_mfrm_sim_spec()] / [extract_mfrm_sim_spec()] ->
#'    [evaluate_mfrm_design()] / [predict_mfrm_population()] ->
#'    [predict_mfrm_units()] / [sample_mfrm_plausible_values()].
#'    Fixed-calibration unit scoring currently requires an `MML` fit, and
#'    prediction export requires actual prediction objects in addition to
#'    `include = "predictions"`.
#' 9. Use `summary()` for compact text checks and `plot()` (or dedicated plot
#'    helpers) for base-R visual diagnostics.
#'
#' @section Three practical routes:
#' - Quick first pass:
#'   [fit_mfrm()] -> [diagnose_mfrm()] -> [plot_qc_dashboard()].
#' - Linking and coverage review:
#'   [subset_connectivity_report()] -> `plot(..., type = "design_matrix")` ->
#'   [plot_wright_unified()].
#' - Manuscript prep:
#'   [reporting_checklist()] -> [build_apa_outputs()] -> [apa_table()].
#' - Design planning and forecasting:
#'   [build_mfrm_sim_spec()] or [extract_mfrm_sim_spec()] ->
#'   [evaluate_mfrm_design()] -> [predict_mfrm_population()] ->
#'   [predict_mfrm_units()] or [sample_mfrm_plausible_values()] from an `MML`
#'   fit -> `export_mfrm_bundle(population_prediction = ..., unit_prediction = ...,
#'   plausible_values = ..., include = "predictions", ...)`.
#'
#' @section Interpreting output:
#' This help page is a map, not an estimator:
#' - use it to decide function order,
#' - confirm which objects have `summary()`/`plot()` defaults,
#' - and identify when dedicated helper functions are needed.
#'
#' @section Objects with default `summary()` and `plot()` routes:
#' - `mfrm_fit`: `summary(fit)` and `plot(fit, ...)`.
#' - `mfrm_diagnostics`: `summary(diag)`; plotting via dedicated helpers
#'   such as [plot_unexpected()], [plot_displacement()], [plot_qc_dashboard()].
#' - `mfrm_bias`: `summary(bias)` and [plot_bias_interaction()].
#' - `mfrm_data_description`: `summary(ds)` and `plot(ds, ...)`.
#' - `mfrm_anchor_audit`: `summary(aud)` and `plot(aud, ...)`.
#' - `mfrm_facets_run`: `summary(run)` and `plot(run, type = c("fit", "qc"), ...)`.
#' - `apa_table`: `summary(tbl)` and `plot(tbl, ...)`.
#' - `mfrm_apa_outputs`: `summary(apa)` for compact diagnostics of report text.
#' - `mfrm_threshold_profiles`: `summary(profiles)` for preset threshold grids.
#' - `mfrm_population_prediction`: `summary(pred)` for design-level forecast
#'   tables.
#' - `mfrm_unit_prediction`: `summary(pred)` for fixed-calibration unit-level
#'   posterior summaries.
#' - `mfrm_plausible_values`: `summary(pv)` for draw-level uncertainty
#'   summaries.
#' - `mfrm_bundle` families:
#'   `summary()` and class-aware `plot(bundle, ...)`.
#'   Key bundle classes now also use class-aware `summary(bundle)`:
#'   `mfrm_unexpected`, `mfrm_fair_average`, `mfrm_displacement`,
#'   `mfrm_interrater`, `mfrm_facets_chisq`, `mfrm_bias_interaction`,
#'   `mfrm_rating_scale`, `mfrm_category_structure`, `mfrm_category_curves`,
#'   `mfrm_measurable`, `mfrm_unexpected_after_bias`, `mfrm_output_bundle`,
#'   `mfrm_residual_pca`, `mfrm_specifications`, `mfrm_data_quality`,
#'   `mfrm_iteration_report`, `mfrm_subset_connectivity`,
#'   `mfrm_facet_statistics`, `mfrm_parity_report`, `mfrm_reference_audit`,
#'   `mfrm_reference_benchmark`.
#'
#' @section `plot.mfrm_bundle()` coverage:
#' Default dispatch now covers:
#' - `mfrm_unexpected`, `mfrm_fair_average`, `mfrm_displacement`
#' - `mfrm_interrater`, `mfrm_facets_chisq`, `mfrm_bias_interaction`
#' - `mfrm_bias_count`, `mfrm_fixed_reports`, `mfrm_visual_summaries`
#' - `mfrm_category_structure`, `mfrm_category_curves`, `mfrm_rating_scale`
#' - `mfrm_measurable`, `mfrm_unexpected_after_bias`, `mfrm_output_bundle`
#' - `mfrm_residual_pca`, `mfrm_specifications`, `mfrm_data_quality`
#' - `mfrm_iteration_report`, `mfrm_subset_connectivity`, `mfrm_facet_statistics`
#' - `mfrm_parity_report`, `mfrm_reference_audit`, `mfrm_reference_benchmark`
#'
#' For unknown bundle classes, use dedicated plotting helpers or custom base-R
#' plots from component tables.
#'
#' @seealso [fit_mfrm()], [run_mfrm_facets()], [mfrmRFacets()],
#'   [diagnose_mfrm()], [estimate_bias()], [mfrmr_visual_diagnostics],
#'   [mfrmr_reports_and_tables], [mfrmr_reporting_and_apa],
#'   [mfrmr_linking_and_dff], [mfrmr_compatibility_layer],
#'   [summary.mfrm_fit()], `summary(diag)`,
#'   `summary()`, [plot.mfrm_fit()], `plot()`
#'
#' @examples
#' \donttest{
#' toy_full <- load_mfrmr_data("example_core")
#' keep_people <- unique(toy_full$Person)[1:12]
#' toy <- toy_full[toy_full$Person %in% keep_people, , drop = FALSE]
#'
#' fit <- suppressWarnings(fit_mfrm(
#'   toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   maxit = 10
#' ))
#' class(summary(fit))
#'
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' class(summary(diag))
#'
#' t4 <- unexpected_response_table(fit, diagnostics = diag, top_n = 5)
#' class(summary(t4))
#' p <- plot(t4, draw = FALSE)
#'
#' sc <- subset_connectivity_report(fit, diagnostics = diag)
#' p_design <- plot(sc, type = "design_matrix", draw = FALSE, preset = "publication")
#' class(p_design)
#'
#' chk <- reporting_checklist(fit, diagnostics = diag)
#' head(chk$checklist[, c("Section", "Item", "DraftReady", "NextAction")])
#'
#' sim_spec <- build_mfrm_sim_spec(
#'   n_person = 30,
#'   n_rater = 4,
#'   n_criterion = 4,
#'   raters_per_person = 2,
#'   assignment = "rotating"
#' )
#' pred_pop <- predict_mfrm_population(sim_spec = sim_spec, reps = 2, maxit = 10, seed = 1)
#' summary(pred_pop)$forecast[, c("Facet", "MeanSeparation", "McseSeparation")]
#' }
#'
#' @name mfrmr_workflow_methods
NULL
