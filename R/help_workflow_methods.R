#' mfrmr Workflow and Method Map
#'
#' @description
#' Quick reference for end-to-end `mfrmr` analysis and for checking which
#' output objects support `summary()` and `plot()`.
#'
#' @section Canonical reporting route:
#' For the clearest default route in `RSM` / `PCM`, use
#' [fit_mfrm()] with `method = "MML"` ->
#' [diagnose_mfrm()] with `diagnostic_mode = "both"` ->
#' [reporting_checklist()] ->
#' [plot_qc_dashboard()] and, when flagged, [plot_marginal_fit()] /
#' [plot_marginal_pairwise()] ->
#' [build_apa_outputs()] ->
#' [build_summary_table_bundle()] -> [apa_table()] or
#' [export_summary_appendix()].
#'
#' Use `JML` only when you explicitly want a faster exploratory pass and are
#' willing to defer strict marginal follow-up and formal precision language to
#' a later `MML` run.
#'
#' @section Canonical operational review route:
#' When the main question is scale maintenance rather than manuscript reporting,
#' branch after [diagnose_mfrm()] into:
#' [audit_mfrm_anchors()] and/or [detect_anchor_drift()] ->
#' [build_equating_chain()] when adjacent-link review is needed ->
#' [build_linking_review()] ->
#' inspect `review$group_view_index` for stable wave / link / facet rollups and
#' `summary(review)$plot_routes` for the next plot helper ->
#' [plot_anchor_drift()] or `plot(anchor_audit, ...)` for the specific flagged
#' evidence family.
#'
#' For bounded `GPCM`, keep anchor/drift helpers as direct exploratory support
#' only. [build_linking_review()] remains outside the current formal `GPCM`
#' route.
#'
#' @section Canonical misfit case-review route:
#' When the main question is which observations, facet levels, or pairwise
#' structures deserve follow-up, branch after [diagnose_mfrm()] into:
#' [build_misfit_casebook()] ->
#' inspect `casebook$group_view_index`, `casebook$group_views`, and
#' `summary(casebook)$plot_routes` for stable person / facet / wave rollups and
#' the next plot helper ->
#' [plot_unexpected()], [plot_displacement()], [plot_marginal_fit()], or
#' [plot_marginal_pairwise()] according to `casebook$plot_map` ->
#' [build_summary_table_bundle()] / [export_summary_appendix()] when the
#' flagged cases need appendix-style reporting support.
#'
#' `build_misfit_casebook()` can still be used for bounded `GPCM`, but it
#' should be read as an operational exploratory screen rather than as a strict
#' Rasch-style invariance report.
#'
#' @section Latent-regression route:
#' When the fit uses `population_formula = ...`, keep the distinction between
#' the estimator and the forecast helpers explicit:
#' - [fit_mfrm()] estimates the current narrow latent-regression `MML` branch.
#'   In the returned fit object, `fit$population$person_table` is the
#'   complete-case estimation table, while
#'   `fit$population$person_table_replay` retains the observed-person-aligned
#'   pre-omit background-data table for replay/export provenance.
#' - [predict_mfrm_units()] and [sample_mfrm_plausible_values()] can then score
#'   under the fitted population model when scored units also supply
#'   one-row-per-person background data. That scoring-time `person_data`
#'   contract remains separate from the fit object's stored replay table.
#' - [predict_mfrm_population()] remains a scenario-level simulation/refit
#'   helper rather than the latent-regression estimator itself.
#'
#' @section Score-category support:
#' If the intended rating scale includes categories not observed in the current
#' data, make that support explicit. For example, use
#' `rating_min = 1, rating_max = 5` for a 1-5 scale with only 2-5 observed.
#' If an intermediate category is unobserved (for example 1, 2, 4, 5 with no
#' 3), also set `keep_original = TRUE` if the zero-count category should remain
#' in the fitted support. `summary(describe_mfrm_data(...))` reports retained
#' zero-count categories in `Notes`, printed `Caveats`, and `$caveats`;
#' `summary(fit)` carries full structured rows into printed `Caveats` and
#' `$caveats`, with `Key warnings` as a short triage subset. Summary-table
#' exports route those rows through `score_category_caveats` or
#' `analysis_caveats`. Adjacent threshold estimates should still be treated as
#' weakly identified when an intermediate category is unobserved.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#'    For final reporting, prefer `method = "MML"` unless you explicitly want
#'    a fast exploratory JML pass.
#' 2. (Optional) Use [run_mfrm_facets()] or [mfrmRFacets()] for a
#'    legacy-compatible one-shot workflow wrapper.
#' 3. For `RSM` / `PCM`, build diagnostics with [diagnose_mfrm()].
#'    For final reporting, prefer `diagnostic_mode = "both"` so the legacy
#'    residual path and the strict marginal screen remain visible side by side.
#'    For bounded `GPCM`, diagnostics are now available through
#'    [diagnose_mfrm()] together with [analyze_residual_pca()],
#'    [interrater_agreement_table()], [unexpected_response_table()],
#'    [displacement_table()], [measurable_summary_table()],
#'    [rating_scale_table()], [facet_quality_dashboard()],
#'    [reporting_checklist()], and [plot_qc_dashboard()] with its
#'    fair-average panel retained as an explicit unavailable placeholder.
#'    Treat those residual-based summaries as
#'    exploratory screens because the discrimination parameter is free.
#'    FACETS-style fair averages are Rasch-family measure-to-score
#'    transformations, so the score-side fair-average semantics remain blocked
#'    for bounded `GPCM`.
#'    Posterior scoring with [predict_mfrm_units()] /
#'    [sample_mfrm_plausible_values()], design-weighted information via
#'    [compute_information()] / [plot_information()], Wright/pathway/CCC plots
#'    via [plot.mfrm_fit()], direct category reports via
#'    [category_structure_report()] / [category_curves_report()], and direct
#'    data generation through [build_mfrm_sim_spec()], [extract_mfrm_sim_spec()],
#'    and [simulate_mfrm_data()] are also available when the simulation
#'    specification stores both thresholds and slopes. Fair-average,
#'    planning/forecasting, and APA/QC pipelines remain outside the validated
#'    `GPCM` boundary. Use [gpcm_capability_matrix()] as the formal capability
#'    map before branching into less common helpers.
#' 4. (Optional, `RSM` / `PCM`) Estimate interaction bias with [estimate_bias()].
#' 5. (Optional, `RSM` / `PCM`) Choose a downstream branch:
#'    [reporting_checklist()] for manuscript/report preparation, or
#'    [build_weighting_audit()] for Rasch-versus-bounded-`GPCM`
#'    weighting review, or [build_misfit_casebook()] / [build_linking_review()]
#'    for operational case review.
#' 6. (Optional, `RSM` / `PCM`) Generate reporting bundles:
#'    [build_summary_table_bundle()], [apa_table()],
#'    [export_summary_appendix()], [build_fixed_reports()],
#'    [build_visual_summaries()]. Weighting-review surfaces can also be routed
#'    through [build_summary_table_bundle()] -> [apa_table()] /
#'    [export_summary_appendix()]. Misfit-case review surfaces now use the same
#'    bundle/export handoff after [build_misfit_casebook()].
#' 7. (Optional, `RSM` / `PCM`) Audit report completeness with
#'    [reference_case_audit()]. Use `facets_parity_report()` only when you
#'    explicitly need the compatibility layer.
#' 8. (Optional, `RSM` / `PCM`) For operational linking follow-up, combine
#'    [audit_mfrm_anchors()], [detect_anchor_drift()], and
#'    [build_equating_chain()] inside [build_linking_review()] before
#'    exporting appendix-style tables.
#' 9. (Optional) Check packaged reference cases with
#'    [reference_case_benchmark()] when you want package-side reference checks.
#' 10. (Optional) For design planning or future scoring, move to the
#'    simulation/prediction layer:
#'    [build_mfrm_sim_spec()] / [extract_mfrm_sim_spec()] ->
#'    [evaluate_mfrm_design()] / [predict_mfrm_population()] ->
#'    [predict_mfrm_units()] / [sample_mfrm_plausible_values()]. Current
#'    fit-derived simulation specs include direct `GPCM` data generation, but
#'    design-evaluation / forecasting helpers still remain `RSM` / `PCM` only
#'    and still target the role-based person x rater-like x criterion-like
#'    contract.
#'    Unit scoring can use an ordinary `MML` fit directly, a latent-regression
#'    `MML` fit when you also supply one-row-per-person background data for the
#'    scored units, or a `JML` fit when a post hoc reference-prior EAP layer is
#'    acceptable. Intercept-only latent-regression fits
#'    (`population_formula = ~ 1`) can reconstruct that minimal person table
#'    from the scored person IDs. Keep `predict_mfrm_population()`
#'    conceptually separate from that scoring layer: it is a simulation-based
#'    scenario forecast helper, not the latent-regression estimator itself.
#'    Prediction export still requires actual prediction objects in addition to
#'    `include = "predictions"`.
#' 10. Use `summary()` for compact text checks and `plot()` (or dedicated plot
#'    helpers) for base-R visual diagnostics.
#'
#' @section Three practical routes:
#' - Quick first pass:
#'   `RSM` / `PCM`: [fit_mfrm()] -> [diagnose_mfrm()] -> [plot_qc_dashboard()] ->
#'   [reporting_checklist()] when you want the package to route the next figures.
#'   bounded `GPCM`: [fit_mfrm()] -> [diagnose_mfrm()] ->
#'   [plot_qc_dashboard()] / [unexpected_response_table()] ->
#'   [rating_scale_table()] ->
#'   [compute_information()] -> [plot_information()] ->
#'   [plot.mfrm_fit()] / [category_curves_report()]. Keep bounded `GPCM`
#'   routes on the direct table/plot side; the fit-based export family
#'   ([build_mfrm_manifest()], [build_mfrm_replay_script()],
#'   [export_mfrm_bundle()]) remains outside the formal `GPCM` boundary.
#' - Linking and coverage review:
#'   [subset_connectivity_report()] -> `plot(..., type = "design_matrix")` ->
#'   [plot_wright_unified()].
#' - Manuscript prep:
#'   `RSM` / `PCM`:
#'   [reporting_checklist()] -> inspect the `"Visual Displays"` and
#'   `"Method Section"` rows -> [build_apa_outputs()] ->
#'   [build_summary_table_bundle()] -> [apa_table()] or
#'   [export_summary_appendix()].
#'   First-release `GPCM`:
#'   [reporting_checklist()] -> direct table/plot helpers while the APA writer
#'   remains outside scope.
#' - Weighting-policy review:
#'   [compare_mfrm()] -> [build_weighting_audit()] ->
#'   [compute_information()] / [plot_information()] when you want to inspect
#'   whether bounded `GPCM` is introducing substantively acceptable
#'   discrimination-based reweighting relative to the Rasch-family reference.
#' - Design planning and forecasting:
#'   [build_mfrm_sim_spec()] or [extract_mfrm_sim_spec()] ->
#'   [evaluate_mfrm_design()] -> [predict_mfrm_population()] ->
#'   [predict_mfrm_units()] or [sample_mfrm_plausible_values()] under the
#'   fitted scoring basis (ordinary `MML`, latent-regression `MML` with
#'   person-level background data, or `JML` with the documented post hoc EAP
#'   approximation). Here again, [predict_mfrm_population()] is the
#'   scenario-level forecast helper, whereas [predict_mfrm_units()] /
#'   [sample_mfrm_plausible_values()] are the scoring layer. Prediction export
#'   requires actual prediction objects. First-release `GPCM` now supports
#'   direct data generation via
#'   [build_mfrm_sim_spec()], [extract_mfrm_sim_spec()], and
#'   [simulate_mfrm_data()], residual diagnostics, and direct curve/report
#'   helpers, but still stops before planning/forecasting helpers. The current
#'   planning layer remains role-based for two non-person facets even though
#'   estimation itself supports arbitrary facet counts; future arbitrary-facet
#'   planning fields should be treated as design metadata rather than finished
#'   public behavior.
#'
#' @section Interpreting output:
#' This help page is a map, not an estimator:
#' - use it to decide function order,
#' - confirm which objects have `summary()`/`plot()` defaults,
#' - identify when dedicated helper functions are needed,
#' - and treat [reporting_checklist()] as the package's readiness router for
#'   plot and report follow-up.
#'
#' @section Objects with default `summary()` and `plot()` routes:
#' - `mfrm_fit`: `summary(fit)` and `plot(fit, ...)`.
#' - `mfrm_diagnostics`: `summary(diag)`; plotting via dedicated helpers
#'   such as [plot_unexpected()], [plot_displacement()], [plot_qc_dashboard()].
#' - `mfrm_bias`: `summary(bias)` and [plot_bias_interaction()].
#' - `mfrm_data_description`: `summary(ds)` and `plot(ds, ...)`.
#' - `mfrm_anchor_audit`: `summary(aud)` and `plot(aud, ...)`.
#' - `mfrm_misfit_casebook`: `summary(casebook)` and `print(casebook)`, with
#'   grouping views available through `casebook$group_view_index` and
#'   `casebook$group_views`, source-specific plotting routed through
#'   `summary(casebook)$plot_routes` and `casebook$plot_map`, and
#'   appendix/report handoff available through
#'   [build_summary_table_bundle()] and [export_summary_appendix()].
#' - `mfrm_weighting_audit`: `summary(audit)` and `print(audit)`, with
#'   information follow-up routed through [compute_information()] and
#'   [plot_information()] according to `audit$plot_map`, and appendix/report
#'   handoff available through [build_summary_table_bundle()] and
#'   [export_summary_appendix()].
#' - `mfrm_linking_review`: `summary(review)` and `print(review)`, with
#'   grouping views available through `review$group_view_index` and
#'   `review$group_views`, and plotting routed through `summary(review)$plot_routes`,
#'   [plot_anchor_drift()], and `plot(anchor_audit, ...)` according to
#'   `review$plot_map`.
#' - `mfrm_facets_run`: `summary(run)` and `plot(run, type = c("fit", "qc"), ...)`.
#' - `apa_table`: `summary(tbl)` and `plot(tbl, ...)`.
#' - `mfrm_apa_outputs`: `summary(apa)` for compact diagnostics of report text.
#' - `mfrm_summary_table_bundle`: `print(bundle)` for manuscript-oriented table
#'   index plus named tables from supported `summary()` outputs,
#'   `summary(bundle)` for table-role/numeric coverage, and `plot(bundle, ...)`
#'   for table-size or numeric-column QC.
#' - `mfrm_threshold_profiles`: `summary(profiles)` for preset threshold grids.
#' - `mfrm_population_prediction`: `summary(pred)` for design-level forecast
#'   tables.
#' - `mfrm_unit_prediction`: `summary(pred)` for unit-level posterior summaries
#'   under the fitted scoring basis.
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
#'   [gpcm_capability_matrix], [mfrmr_linking_and_dff],
#'   [mfrmr_compatibility_layer],
#'   [summary.mfrm_fit()], `summary(diag)`,
#'   `summary()`, [plot.mfrm_fit()], `plot()`
#'
#' @examples
#' \donttest{
#' toy_full <- load_mfrmr_data("example_core")
#' keep_people <- unique(toy_full$Person)[1:12]
#' toy <- toy_full[toy_full$Person %in% keep_people, , drop = FALSE]
#'
#' fit <- fit_mfrm(
#'   toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   maxit = 200
#' )
#' summary(fit)$next_actions
#'
#' diag <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "both")
#' summary(diag)$next_actions
#'
#' chk <- reporting_checklist(fit, diagnostics = diag)
#' subset(
#'   chk$checklist,
#'   Section == "Visual Displays",
#'   c("Item", "DraftReady", "NextAction")
#' )
#'
#' qc <- plot_qc_dashboard(fit, diagnostics = diag, draw = FALSE, preset = "publication")
#' qc$data$preset
#' p_marg <- plot_marginal_fit(diag, draw = FALSE, preset = "publication")
#' p_marg$data$preset
#'
#' sc <- subset_connectivity_report(fit, diagnostics = diag)
#' p_design <- plot(sc, type = "design_matrix", draw = FALSE, preset = "publication")
#' p_design$data$plot
#'
#' bundle <- build_summary_table_bundle(chk, appendix_preset = "recommended")
#' summary(bundle)$role_summary
#' plot(bundle, type = "appendix_presets", draw = FALSE)$data$plot
#' }
#'
#' @name mfrmr_workflow_methods
NULL
