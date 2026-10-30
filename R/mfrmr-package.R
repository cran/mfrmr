#' mfrmr: Many-Facet Ordered-Response Modeling in R
#'
#' @description
#' `mfrmr` provides estimation, diagnostics, and reporting utilities for
#' many-facet ordered-response measurement models: the Rasch-family `RSM` /
#' `PCM` route and a `GPCM` extension in which one selected facet supplies
#' level-specific discriminations. MML permits a different facet to supply
#' category steps; JML requires the same facet for both roles. [confint.mfrm_fit()]
#' with `parm = "slopes"`
#' supplies approximate relative-slope intervals for eligible MML fits; model
#' ranking and matched PCM/GPCM tests use separate checks in [compare_mfrm()].
#' [gpcm_capability_matrix()] explains the supported uses.
#' An experimental two-family `GPCM` route also provides estimation and
#' conditional scoring, with explicit uncertainty limitations.
#' Native RSM/PCM conditional maximum likelihood (`CML`) conditions on Person
#' totals without assuming an ability distribution. Its separate result class
#' supplies structural contrast intervals through [mfrm_facet_intervals()],
#' without Person abilities, a population SD or ordinary fit diagnostics.
#' See [fit_mfrm()] for the supported conditional estimation scope.
#'
#' For assessment planning from observed scores, [mfrm_g_study()]
#' and [mfrm_d_study()] separate sources of variation and project
#' dependability for specified future rating plans. For rater feedback and
#' assessment decisions, [mfrm_report()], [interrater_agreement_table()] and
#' [classification_table()] provide their respective reports, agreement
#' measures and fixed-plan classification summaries. These analyses have
#' different assumptions and uncertainty targets; follow the linked help for
#' the question you want to answer.
#'
#' @useDynLib mfrmr, .registration = TRUE
#'
#' @details
#' Start with the complete script in the Examples section:
#'
#' 1. Load the package and example ratings with [load_mfrmr_data()].
#' 2. Fit a model with [fit_mfrm()].
#' 3. Draw the Wright map with `plot(fit)`.
#' 4. Save `results <- summary(fit)` and inspect `results$person_overview`
#'    and `results$facet_overview`.
#'
#' Each data row is one rating event. The `person`, `facets`, and `score`
#' arguments name columns. The example uses marginal maximum likelihood
#' (`MML`) and a rating-scale model (`RSM`) with shared category thresholds.
#' `head(ratings)` displays the first six rows. Only `Person`, `Rater`, `Criterion`,
#' and `Score` are used by this model; the example's `Study` and `Group` columns
#' are additional labels. `<-` saves an object, and `$` selects a named part:
#' `results$person_overview` displays one table from the saved summary.
#' The overview tables summarize distributions; `as.data.frame(fit)` returns
#' individual person, rater, and criterion estimates.
#' The assignment `ratings <- load_mfrmr_data("example_operational")` reads
#' an existing synthetic dataset supplied with mfrmr. You do not need to create
#' `ratings` first; the name on the left stores what the function returns.
#'
#' Before interpreting or reporting estimates, read `results$decision` and
#' follow its `NextAction`. The default summary does not compute diagnostics.
#' For your own data, first inspect the design and score categories with
#' [describe_mfrm_data()].
#'
#' @section Optional packages:
#' The basic example below uses packages installed with mfrmr. Its base-R
#' plots do not need ggplot2. Additional packages support particular features;
#' you do not need to install all suggested packages to get started.
#'
#' - [as_ggplot()] needs ggplot2 for ggplot2 rendering.
#' - [as_kable()] needs knitr and [as_flextable()] needs flextable for their
#'   respective table formats.
#' - [fit_mfrm_random_rater()] needs RTMB version 2.0 or later for shared
#'   random-rater estimation.
#' - Corrected JML in [fit_mfrm()], requested with `jml_correction_order`,
#'   needs nleqslv to solve the adjusted equations.
#' - [ordinal_logistic_table()] needs ordinal for ordinal-logistic DIF.
#'
#' These are common examples; other optional analyses name their dependencies
#' in their own help. A missing-package message identifies what to install.
#' For example, use `install.packages("ggplot2")` if you want [as_ggplot()].
#'
#' A saved [mfrm_results()] bundle can be reopened with [base::readRDS()] and
#' passed to [mfrm_report()] without refitting or installing optional estimation
#' packages. Reporting reuses stored results and retains unavailable diagnostics
#' as unavailable. Updating an analysis or computing new diagnostics is a
#' separate operation and may need the original dependencies. Save the complete
#' bundle with [base::saveRDS()] and retain its matching package version/source;
#' see "Updating saved analyses for 0.2.4" in [mfrmr_workflow_methods] before
#' reusing older saved analyses.
#'
#' @section Measurement scale:
#' Native response models use the logistic metric (L-metric), with scaling
#' constant \eqn{D = 1}; no extra 1.7 normal-ogive approximation factor is used.
#' This is distinct from choosing a normal ability distribution. See the
#' "Logistic metric (D = 1)" section of [fit_mfrm()] for logit interpretation
#' with fixed or varying slopes and comparison with other software.
#'
#' @section Where to go next:
#' - [mfrmr_methods_reference] for the question-to-calculation map, plain-language
#'   explanations and links to the methodological literature.
#' - [mfrmr_workflow_methods] and `vignette("mfrmr-workflow", package = "mfrmr")`
#'   for your own CSV, column mapping, input troubleshooting, and the full
#'   analysis workflow. If vignettes are not installed, start with
#'   [describe_mfrm_data()] and [recode_missing_codes()].
#' - [mfrmr_visual_diagnostics] for choosing follow-up figures.
#' - [mfrmr_reporting_and_apa] and [mfrmr_reports_and_tables] for reporting.
#' - [mfrmr_linking_and_dff] for linking and differential facet functioning.
#' - [classification_table()] for fixed-plan future classification accuracy,
#'   consistency and conditional precision under ordinary RSM/PCM calibration.
#' - [gpcm_capability_matrix] for the `GPCM` extension.
#' - [mfrmr_output_guide()] for the broader purpose-to-function map.
#'
#' A printable reference card is available at
#' `system.file("cheatsheet", "mfrmr-cheatsheet.pdf", package = "mfrmr")`.
#'
#' @section Portable fixed calibration:
#' A saved, versioned calibration artifact can be created from an eligible
#' one-scale `RSM` or `PCM` MML fit under the fixed standard-normal scoring
#' basis. Use [mfrm_calibration_capabilities()] before extraction, then follow
#' [mfrm_calibration_workflow] to review, validate, freeze, save, load, and
#' score the artifact. Review the returned batch with
#' [mfrm_calibration_score_methods] before using its estimates. One-family GPCM
#' MML also supports an estimated intercept-only normal population, shared or
#' separate slope/step owners, unit weights and no anchors/interactions, subject
#' to conditional source and per-batch numerical scoring checks. Estimated-
#' population RSM/PCM and latent-regression MML remain fitted-object-only.
#' RSM/PCM JML supports portable post-hoc EAP with a standard-normal reference
#' prior, unit weights, no anchors/interactions and passing finite identified
#' source checks. Shared-owner GPCM JML also supports conditional portable
#' EAP after its joint-likelihood/curvature checks, preserving incomplete global
#' audits. Explicit corrected GPCM JML instead uses adjusted-equation checks
#' and preserves the correction order in format 5; residual calibration bias
#' may remain. That reference prior is not estimated by JML. Experimental
#' two-family GPCM MML uses format 6 with both ordered slope owners and a fixed
#' N(0,1) prior, separate source/batch checks and no prior overrides.
#' Artifact
#' score uncertainty is conditional on the frozen point calibration and its
#' recorded prior; loading validates consistency but does not authenticate an
#' untrusted file.
#'
#' @section Advanced scope:
#' After the basic route above:
#' - the package supports latent-regression `MML`
#'   for ordered-response `RSM` / `PCM` models with a one-dimensional
#'   conditional-normal population model and explicit one-row-per-person
#'   covariates expanded through `stats::model.matrix()`
#' - `GPCM` support is summarized by [gpcm_capability_matrix()]
#' - `GPCM` supports the core fit/summary/scoring/information
#'   path, direct Wright/pathway/CCC plots, residual-PCA follow-up, and the
#'   residual-based diagnostics tables/plots as exploratory tools
#' - posterior-predictive checks and `MCMC` estimation are not available for
#'   `GPCM`; use external Bayesian software when they are required
#' - direct `GPCM` data generation through [build_mfrm_sim_spec()],
#'   [extract_mfrm_sim_spec()], and [simulate_mfrm_data()] is available when
#'   the specification carries both thresholds and slopes with the same owner;
#'   general simulation/design workflows do not support separate owners
#' - slope-aware [fair_average_table()] and [estimate_bias()] are available for
#'   `GPCM` with explicit caveats; [build_apa_outputs()],
#'   [build_visual_summaries()], [run_qc_pipeline()],
#'   [build_mfrm_manifest()], [build_mfrm_replay_script()], and
#'   [export_mfrm_bundle()] are available as caveated partial reporting/export
#'   surfaces; score-side FACETS compatibility remains limited, and fully
#'   arbitrary-facet planning is outside the documented planning scope
#' - `predict_mfrm_population()` remains a scenario-level forecast helper and
#'   should not be described as the latent-regression estimator itself
#' - the current simulation/planning layer remains role-based for two
#'   non-person facets rather than fully arbitrary-facet planning, with boundaries
#'   exposed through planner metadata such as `planning_scope`,
#'   `planning_constraints`, and `planning_schema`
#' - latent-class mixture models and response-time / careless-rating
#'   adjustment are not estimated by mfrmr; use residual, person-fit,
#'   local-dependence, and rater-drift diagnostics as screening layers rather
#'   than as mixture-model substitutes
#'
#' @section Equal weighting versus GPCM:
#' The package's operational reference route is the Rasch-family
#' `RSM` / `PCM` branch. That route enforces fixed discrimination and therefore
#' preserves an equal-weighting scoring interpretation across observed ratings.
#'
#' `GPCM` is supported because some users want a slope-aware model-
#' comparison or sensitivity layer inside the same many-facet workflow. However,
#' the package does not treat `GPCM` as a universal replacement for the
#' Rasch-family route. A better fit under `GPCM` should be read as evidence
#' about discrimination-based reweighting, not as an automatic reason to
#' discard the equal-weighting model.
#'
#' Observation weights are a different concept again. Optional `Weight`
#' columns change how observed rating events enter estimation and summaries, but
#' they do not create a free-form facet-weighting scheme and do not alter the
#' fixed-discrimination meaning of `RSM` / `PCM`.
#'
#' Public entry map:
#' - First-screen results: [mfrm_results()], `summary(res)$next_actions`,
#'   and `mfrmr_output_guide("entry")`
#' - Interactive exploration: [mfrm_results_interactive()] only when prompts
#'   are explicitly wanted at the console
#'
#' Function families:
#' - Model fitting: [fit_mfrm()], [summary.mfrm_fit()], [plot.mfrm_fit()]
#' - Legacy-compatible workflow wrapper: [run_mfrm_facets()], [mfrmRFacets()]
#' - Diagnostics: [diagnose_mfrm()], `summary(diag)`,
#'   [analyze_residual_pca()], [plot(x, type = "residual_pca")][plot.mfrm_fit]
#' - Bias and interaction: [estimate_bias()], [estimate_all_bias()],
#'   `summary(bias)`, [bias_interaction_report()], [plot(x, type = "bias")][plot.mfrm_fit]
#' - Differential functioning: [analyze_dff()], [analyze_dif()],
#'   [mantel_haenszel_table()] for separately matched binary responses,
#'   [ordinal_logistic_table()] for matched cumulative-logit comparisons,
#'   [dif_interaction_table()], [plot(x, type = "heatmap")][plot.mfrm_dif], [dif_report()]
#' - Design simulation: [build_mfrm_sim_spec()], [extract_mfrm_sim_spec()],
#'   [simulate_mfrm_data()], [simulate_mfrm_recovery()],
#'   [assess_mfrm_recovery()], [evaluate_mfrm_design()],
#'   [evaluate_mfrm_signal_detection()], [predict_mfrm_population()],
#'   [predict_mfrm_units()], [sample_mfrm_plausible_values()] (including
#'   fit-derived empirical / resampled / skeleton-based simulation
#'   specifications; fitted-object posterior scoring supports `MML` fits
#'   directly, latent-regression `MML` fits through the fitted population
#'   model when scored units also provide one-row-per-person background data,
#'   and `JML` fits through a post hoc reference-prior EAP layer;
#'   fit-derived simulation specifications also support direct
#'   `GPCM` data generation, recovery checks, role-based design evaluation,
#'   population forecasting, diagnostic-screening, and signal-detection
#'   helpers with documented caveats; curve reports, graph-only exports,
#'   fair-average tables, and bias screening are also available for
#'   `GPCM` with documented caveats)
#' - Reporting: [build_apa_outputs()], [build_visual_summaries()],
#'   [reporting_checklist()], [apa_table()] for the full `RSM` / `PCM` route;
#'   `GPCM` uses these as caveated partial surfaces and retains direct
#'   table, plot, checklist, and summary-appendix routes
#' - Weighting review: [compare_mfrm()], [build_weighting_review()],
#'   [build_model_choice_review()], [compute_information()], [plot(x)][plot.mfrm_information]
#' - Case review: [build_misfit_casebook()], [plot(x, type = "unexpected")][plot.mfrm_fit],
#'   [plot(x, type = "displacement")][plot.mfrm_fit], [plot(x, type = "marginal_fit")][plot.mfrm_fit], [plot(x, type = "marginal_pairwise")][plot.mfrm_fit]
#' - Linking and scale maintenance: [review_mfrm_anchors()],
#'   [detect_anchor_drift()], [build_equating_chain()],
#'   [build_linking_review()], [plot(x)][plot.mfrm_anchor_drift]
#' - Dashboards: [facet_quality_dashboard()], [plot(x, type = "facet_dashboard")][plot.mfrm_fit]
#' - Export / reproducibility: [build_mfrm_manifest()], [build_mfrm_replay_script()],
#'   [build_conquest_overlap_bundle()], `normalize_conquest_overlap_exports()`,
#'   [normalize_conquest_overlap_files()],
#'   [normalize_conquest_overlap_tables()],
#'   [review_conquest_overlap()],
#'   [export_mfrm_bundle()] for the diagnostics-compatible Rasch-family route;
#'   `GPCM` supports caveated partial manifest, replay, and bundle
#'   output as well as [export_summary_appendix()]
#' - Equivalence: [analyze_facet_equivalence()], [plot(x, type = "facet_equivalence")][plot.mfrm_fit]
#' - Data and anchors: [describe_mfrm_data()], [review_mfrm_anchors()],
#'   [make_anchor_table()], [load_mfrmr_data()]
#'
#' Data interface:
#' - Input analysis data is long format (one row per observed rating).
#' - Required columns are one person column, one ordered score column, and one
#'   or more non-person facet columns named in `facets = c(...)`.
#' - Score values should be ordered integer categories. Binary `0/1` or `1/2`
#'   input is supported as the two-category Rasch-family special case; by
#'   contrast, fractional score values should be recoded before fitting rather
#'   than relying on automatic coercion.
#' - If `keep_original = FALSE`, unused intermediate categories are collapsed
#'   to a contiguous internal scale and the mapping is stored in
#'   `fit$prep$score_map`.
#' - If the intended scale has unused boundary categories, such as a 1-5 scale
#'   with only 2-5 observed, set `rating_min = 1, rating_max = 5` so the
#'   zero-count boundary category remains explicit in the data-support review.
#'   A missing boundary remains review evidence for the separate element-
#'   boundary contract. If unused intermediate categories should remain in the
#'   original scale, set `keep_original = TRUE`; a retained zero-count internal
#'   category in a polytomous fitted ladder creates an unsupported adjacent-
#'   step contrast, and [fit_mfrm()] returns a structured pre-optimization error.
#' - `summary(describe_mfrm_data(...))` reports retained zero-count categories
#'   in `Notes`, printed `Caveats`, and `$caveats`; `summary(fit)` carries full
#'   structured rows into printed `Caveats` and `$caveats`, with `Key warnings`
#'   as a short triage subset. Summary-table exports route those rows through
#'   `score_category_caveats` or `analysis_caveats`.
#' - A category can be observed globally but unused by one rater or other facet
#'   level. Review that local pattern with [data_quality_report()] rather than
#'   treating it as a globally missing step or an automatic reason to select
#'   `GPCM`.
#' - Optional columns such as `Subset`, `Weight`, and `Group` support linking,
#'   weighted analysis, and fairness-focused follow-up workflows.
#' - Packaged synthetic data is available via [load_mfrmr_data()] or `data()`.
#'   Use `list_mfrmr_data(details = TRUE)` to distinguish the applied teaching,
#'   idealized, planted-effect, and larger sparse-design datasets.
#'
#' @section Interpreting output:
#' Core object classes are:
#' - `mfrm_fit`: fitted model parameters and metadata.
#' - `mfrm_diagnostics`: fit, facet-level reliability, and flag diagnostics,
#'   plus inter-rater agreement when one facet is treated as a rater facet.
#' - `mfrm_bias`: interaction bias estimates.
#' - `mfrm_dff` / `mfrm_dif`: differential-functioning contrasts and screening summaries.
#' - `mfrm_population_prediction`: scenario-level forecast summaries for one
#'   future design.
#' - `mfrm_unit_prediction`: posterior summaries for future or partially
#'   observed persons under the fitted scoring basis.
#' - `mfrm_plausible_values`: posterior draws for future or partially observed
#'   persons under the fitted scoring basis.
#' - `mfrm_bundle` families: summary/report bundles and draw-free plot data.
#'
#' @section Typical workflow:
#' 1. Prepare long-format data.
#' 2. Fit with [fit_mfrm()].
#' 3. For `RSM` / `PCM`, diagnose with [diagnose_mfrm()] and prefer
#'    `diagnostic_mode = "both"` for final `MML` runs.
#' 4. For `RSM` / `PCM`, run [analyze_dff()] or [estimate_bias()] when
#'    fairness or interaction questions matter; `GPCM` also supports
#'    [estimate_bias()] as a conditional screening review.
#' 5. For `RSM` / `PCM`, report with [build_apa_outputs()] and
#'    [build_visual_summaries()].
#' 6. For design planning, move to [build_mfrm_sim_spec()],
#'    [evaluate_mfrm_design()], [mfrm_generalizability()], [mfrm_d_study()],
#'    and [predict_mfrm_population()].
#'    `GPCM` also supports direct simulation via
#'    [extract_mfrm_sim_spec()] / [simulate_mfrm_data()] and caveated role-based
#'    [evaluate_mfrm_design()] / [predict_mfrm_population()] routes. Those
#'    helpers assume two non-person facet roles
#'    even though the estimation core supports arbitrary facet counts. Treat
#'    [evaluate_mfrm_design()] as Monte Carlo design evaluation, and use
#'    [mfrm_d_study()] for analytic G/Phi design projections. Always read
#'    the `IdentificationStatus`, `GStatus`, and `PhiStatus` columns before
#'    reporting those projections; boundary or singular mixed-model fits are
#'    design-identification warnings, not high-stakes-ready reliability
#'    evidence.
#'    `predict_mfrm_population()` remains the scenario-level forecast helper,
#'    not the latent-regression estimator.
#' 7. For future-unit scoring, retain an `MML` calibration when you want the
#'    fitted marginal model directly, use an active latent-regression `MML`
#'    fit when scored units also provide one-row-per-person background data, or
#'    use a `JML` calibration when a post hoc fitted-object EAP layer is
#'    acceptable; then score with
#'    [predict_mfrm_units()] or [sample_mfrm_plausible_values()].
#' 8. For `GPCM`, use [summary.mfrm_fit()],
#'    [diagnose_mfrm()], [analyze_residual_pca()],
#'    [predict_mfrm_units()], [sample_mfrm_plausible_values()],
#'    [compute_information()], [plot(x, type = "qc")][plot.mfrm_fit], [plot.mfrm_fit()],
#'    [category_structure_report()], [category_curves_report()],
#'    graph-only [facets_output_file_bundle()], direct simulation-spec
#'    generation/data generation, recovery checks, [fair_average_table()],
#'    [estimate_bias()], and the residual-based table helpers with their
#'    documented caveats. Caveated APA/QC/export bundles and exploratory
#'    linking review plus role-based design evaluation, population
#'    forecasting, diagnostic-screening, and signal-detection helpers are
#'    available, while full score-side FACETS review, posterior-predictive
#'    checks, and `MCMC` estimation are not available for `GPCM`. Use
#'    [gpcm_capability_matrix()] as the formal boundary statement.
#'
#' @section Model formulation:
#' The Rasch-family branch used by `RSM` and `PCM` extends the basic Rasch
#' model by incorporating multiple measurement facets into a single additive
#' linear predictor on the log-odds scale.
#'
#' **RSM/PCM adjacent-category equation**
#'
#' For an observation where person \eqn{n} with ability \eqn{\theta_n} is
#' rated by rater \eqn{j} with severity \eqn{\delta_j} on criterion \eqn{i}
#' with difficulty \eqn{\beta_i}, the probability of observing category
#' \eqn{k} (out of \eqn{K} ordered categories) is:
#'
#' \deqn{P(X_{nij} = k \mid \theta_n, \delta_j, \beta_i, \tau) =
#'   \frac{\exp\bigl[\sum_{s=1}^{k}(\theta_n - \delta_j - \beta_i - \tau_s)\bigr]}
#'        {\sum_{c=0}^{K}\exp\bigl[\sum_{s=1}^{c}(\theta_n - \delta_j - \beta_i - \tau_s)\bigr]}}
#'
#' where \eqn{\tau_s} are the Rasch-Andrich threshold (step) parameters in the
#' `RSM` reference case and
#' \eqn{\sum_{s=1}^{0}(\cdot) \equiv 0} by convention.  Additional facets
#' enter as additive terms in the linear predictor
#' \eqn{\eta = \theta_n - \delta_j - \beta_i - \ldots}.
#'
#' This additive predictor generalises to any number of facets; the
#' `facets` argument to [fit_mfrm()] accepts an arbitrary-length
#' character vector.
#'
#' **Rating Scale Model (RSM)**
#'
#' Under the RSM (Andrich, 1978), all levels of the step facet share a
#' single set of threshold parameters \eqn{\tau_1, \ldots, \tau_K}.
#'
#' **Partial Credit Model (PCM)**
#'
#' Under the PCM (Masters, 1982), each level of the designated `step_facet`
#' has its own threshold vector on the package's common observed score scale.
#' In the current implementation, threshold locations may vary by step-facet
#' level, but the fitted score range is defined by one global category
#' set taken from the observed data.
#'
#' **Generalized Partial Credit Model (GPCM)**
#'
#' With one slope family under `GPCM` (Muraki, 1992), the adjacent-category partial-credit
#' kernel is multiplied by a positive slope \eqn{\alpha_g} for the designated
#' slope-facet level \eqn{g}:
#'
#' \deqn{\ln\frac{P(X_{nij} = k)}{P(X_{nij} = k-1)} =
#'   \alpha_g(\theta_n - \delta_j - \beta_i - \tau_{hk}).}
#'
#' MML allows a slope owner \eqn{g} distinct from the step owner \eqn{h};
#' JML requires `slope_facet == step_facet`. Slopes are identified on the log
#' scale with geometric mean 1. This makes
#' `GPCM` a slope-aware sensitivity/extension route, not a replacement for the
#' equal-weighting `RSM`/`PCM` interpretation.
#' A separate provisional two-family MML route fits the multiplicative
#' Uto--Ueno equation with explicit fixed-N(0,1) identification, exactly two
#' facets and second-owner steps. The first family's slopes have geometric mean
#' one and the second family's slopes are free. Fixed-grid EM and adaptive
#' direct MML supply summaries and conditional fitted curves without intervals.
#' Both support separately checked experimental component log-Wald and
#' location/contrast intervals and descriptive posterior residuals; component
#' profiles preserve the fitted fixed-grid or adaptive integration and require
#' separate numerical checks. [predict_mfrm_units()] and the portable
#' calibration workflow supply experimental conditional new-Person EAP and
#' posterior intervals, excluding calibration uncertainty. Coverage and
#' population transport remain unqualified. The one-family diagnostic and
#' inference routes below do not apply to it. See the
#' Two slope families section of [fit_mfrm()]. Unit slopes reduce to the
#' equal-discrimination PCM kernel.
#' Under default one-family MML, an intercept-only person distribution
#' \eqn{N(\beta_0,\sigma^2)} is estimated. The geometric-mean-one slopes are
#' relative discriminations and \eqn{\sigma\alpha_g} are their equivalent
#' fixed-latent-standard-deviation values, so the conventional common
#' discrimination degree of freedom is retained. The legacy
#' `gpcm_mml_identification = "fixed_standard_normal"` mode fixes both
#' \eqn{\sigma=1} and the slope geometric mean to one and is therefore a
#' narrower relative-discrimination model. Under JML, geometric-mean-one is
#' required to resolve the scale of jointly estimated person coordinates.
#' The former "bounded" label described limited model/workflow scope, not
#' finite optimizer box constraints. Use GPCM with explicit scope restrictions.
#' The GPCM JML objective is unpenalized.
#' Certified extreme-person or facet recession is reported through typed
#' primary boundary results; a finite optimizer iterate is retained only as a
#' numerical trace and is not a finite JML maximum.
#'
#' **Ordered-response scope**
#'
#' The implemented response-model scope is ordered categorical only.
#' Binary responses are the \eqn{K = 1} special case of the same formulation,
#' so they are handled through the ordinary ordered-score interface. This means
#' `mfrmr` supports ordered binary and ordered polytomous data under `RSM` and
#' `PCM`, plus `GPCM` with one designated `slope_facet` and, under MML,
#' a possibly different `step_facet`. Unordered
#' nominal/multinomial response models are outside the documented model scope,
#' as are Poisson, negative-binomial, and grouped binomial-trial count-response
#' families. A positive observation `weight` weights the conditional
#' ordered-rating contribution; it is not a general collapsed-person frequency
#' table and does not change the response family or model dependence among
#' replicated ratings. Under MML, powering conditional responses within one
#' Person is not the same as replicating a complete Person pattern after
#' marginalization. Consequently, integer count values supplied as scores are
#' modeled as ordered category codes, not as counts from a Poisson or related
#' distribution.
#'
#' @section Estimation methods:
#' **Marginal Maximum Likelihood (MML)**
#'
#' MML integrates over the person ability distribution using Gauss-Hermite
#' quadrature, in the broader marginal-likelihood framework introduced by
#' Bock & Aitkin (1981) for IRT:
#'
#' \deqn{L = \prod_{n} \int P(\mathbf{X}_n \mid \theta, \boldsymbol{\delta})
#'   \, \phi(\theta) \, d\theta
#'   \approx \prod_{n} \sum_{q=1}^{Q} w_q \,
#'   P(\mathbf{X}_n \mid \theta_q, \boldsymbol{\delta})}
#'
#' where \eqn{\phi(\theta)} is the assumed normal prior and
#' \eqn{(\theta_q, w_q)} are quadrature nodes and weights.  Person
#' estimates are obtained post-hoc via Expected A Posteriori (EAP):
#'
#' \deqn{\hat{\theta}_n^{\mathrm{EAP}} =
#'   \frac{\sum_q \theta_q \, w_q \, L(\mathbf{X}_n \mid \theta_q)}
#'        {\sum_q w_q \, L(\mathbf{X}_n \mid \theta_q)}}
#'
#' MML does not estimate one fixed person parameter per respondent and therefore
#' avoids that JML incidental-parameter mechanism. Its structural inference is
#' nevertheless conditional on the specified population distribution, response
#' model, quadrature approximation, and regularity conditions; it is not an
#' automatic small-sample guarantee.
#'
#' Note: Bock & Aitkin (1981) is the canonical citation for the
#' Gauss-Hermite-quadrature MML *framework*. The default mfrmr engine
#' (`mml_engine = "direct"`) optimises this marginal log-likelihood by
#' direct gradient methods (BFGS / L-BFGS-B), not by Bock & Aitkin's
#' signature EM algorithm. The `"em"` and `"hybrid"` engines do follow
#' the EM template but use the selected gradient-based M-step rather than
#' B&A's probit IRLS,
#' because the target is the polytomous Rasch family rather than B&A's
#' 2PL probit model.
#'
#' **Conditional Maximum Likelihood (CML)**
#'
#' For supported RSM/PCM models, CML conditions on each Person's total and
#' observed assignments, removing ability parameters without specifying their
#' population distribution. A generalized elementary-symmetric-function
#' recursion evaluates the conditional likelihood and information. Identification
#' and finite-MLE checks precede structural contrast intervals; their Wald
#' interpretation is asymptotic, not a general finite-sample coverage guarantee.
#' The separate `mfrm_cml_fit` and comparison-only results supply no Person
#' abilities or population SD. Use [mfrm_facet_intervals()] for their
#' native-scale intervals. CML does not support GPCM or ordinary diagnostics.
#'
#' **Joint Maximum Likelihood (JML)**
#'
#' JML estimates all person and facet parameters simultaneously as fixed
#' effects by maximising the joint log-likelihood
#' \eqn{\ell(\boldsymbol{\theta}, \boldsymbol{\delta} \mid \mathbf{X})}
#' directly.  It does not assume a parametric person distribution, which
#' can be advantageous when the population shape is strongly non-normal.
#' Its incidental-parameter concern is not simply a small-person-sample issue:
#' structural-parameter bias can persist as the number of persons increases
#' when the number of ratings or items per person remains fixed (Neyman &
#' Scott, 1948).
#' The package still accepts `"JMLE"` as a backward-compatible alias, but
#' user-facing summaries and documentation use `"JML"` as the public label.
#'
#' See [fit_mfrm()] for practical guidance on choosing an estimator. General
#' GPCM inference uses MML; JML remains useful for descriptive point estimates
#' and FACETS comparisons.
#'
#' @section Strict marginal diagnostics:
#' For `RSM` / `PCM`, `diagnose_mfrm(..., diagnostic_mode = "both")`
#' returns two complementary targets: the `legacy` residual / EAP
#' diagnostics and a `marginal_fit` layer whose expected counts and
#' pairwise summaries are integrated over the posterior quadrature
#' bundle rather than plugged in at the EAP point. The screen is
#' structured as limited-information evidence (Orlando & Thissen,
#' 2000; Haberman & Sinharay, 2013; Sinharay & Monroe, 2025), not as
#' an omnibus accept / reject test, and it complements rather than
#' replaces separation / reliability and inter-rater agreement
#' summaries. The full derivation, with notation and pairwise
#' local-dependence events, lives in
#' `vignette("mfrmr-mml-and-marginal-fit", package = "mfrmr")`.
#'
#' @section Statistical background:
#' Key statistics reported throughout the package:
#'
#' **Infit (Information-Weighted Mean Square)**
#'
#' Weighted average of squared standardized residuals, where weights are the
#' model-based variance of each observation:
#'
#' \deqn{\mathrm{Infit}_j = \frac{\sum_i Z_{ij}^2 \, \mathrm{Var}_i \, w_i}
#'                               {\sum_i \mathrm{Var}_i \, w_i}}
#'
#' Expected value is 1.0 under model fit.  Values below 0.5 suggest overfit
#' (Mead-style responses); values above 1.5 suggest underfit (noise or
#' misfit).  Infit is most sensitive to unexpected patterns among on-target
#' observations (Wright & Masters, 1982).
#'
#' Note: The 0.5--1.5 range is the general "productive for measurement"
#' band given by Linacre (2002, *RMT* 16(2), 878). Context-specific bands
#' come from Wright & Linacre (1994, *RMT* 8(3), 370): 0.8--1.2 for
#' high-stakes MCQ, 0.7--1.3 for run-of-the-mill MCQ, 0.6--1.4 for
#' rating-scale surveys, 0.5--1.7 for clinical observation, and 0.4--1.2
#' for judged performance. See also Bond & Fox (2015) for textbook
#' summaries of these conventions.
#'
#' **Outfit (Unweighted Mean Square)**
#'
#' Simple average of squared standardized residuals:
#'
#' \deqn{\mathrm{Outfit}_j = \frac{\sum_i Z_{ij}^2 \, w_i}{\sum_i w_i}}
#'
#' Same expected value and flagging thresholds as Infit, but more sensitive
#' to extreme off-target outliers (e.g., a high-ability person scoring the
#' lowest category).
#'
#' **ZSTD (Standardized Fit Statistic)**
#'
#' Wilson-Hilferty (1931) cube-root transformation that converts the
#' mean-square chi-square ratio to an approximate standard normal deviate:
#'
#' \deqn{\mathrm{ZSTD} = \frac{\mathrm{MnSq}^{1/3} - (1 - 2/(9\,\mathit{df}))}
#'                            {\sqrt{2/(9\,\mathit{df})}}}
#'
#' Values near 0 indicate expected fit. The conventional
#' \eqn{|\mathrm{ZSTD}| \ge 2} and \eqn{|\mathrm{ZSTD}| \ge 3} cutoffs are heuristic
#' two- and three-standard-deviation reference bands (Wright & Linacre, 1994;
#' see also Wilson & Hilferty, 1931), not calibrated 5\% and 1\% hypothesis
#' tests. Parameter estimation, sparse cells, the selected df convention, and
#' repeated screening across elements all affect that interpretation.
#' Diagnostics default to FACETS fourth-moment df, admitting df > 0 and
#' capping ZSTD at +/-9. The explicit engine convention withholds df < 1.
#' `whexact = FALSE` uses the expression above (FACETS YES); `TRUE` removes
#' its additive q/3 term (FACETS NO). Both use a cube root. The former linear
#' TRUE transformation is available only as `whexact = "legacy_linear"`.
#'
#' **Residual basis under MML vs JMLE engines**
#'
#' For `method = "MML"` fits, the standardized residuals behind Infit,
#' Outfit, and ZSTD are evaluated at EAP person measures, which are
#' shrunken toward the population mean. JMLE engines such as FACETS
#' evaluate the same formulas at unshrunken JMLE estimates, so MnSq and
#' ZSTD values are not numerically interchangeable across the two residual
#' bases, most visibly for extreme-scoring persons. Use `method = "JML"`
#' when an external FACETS fit comparison requires a JMLE-style residual
#' basis, and see [facets_fit_df_guide()] for the separate
#' standardization-side df/ZSTD conventions.
#'
#' **PTMEA (Point-Measure Correlation)**
#'
#' Pearson correlation between observed scores and estimated person measures
#' within each facet level. A positive value is directionally consistent with
#' the fitted orientation, but it is not a confirmatory test of alignment.
#' A negative value prompts review of orientation, response patterns, and model
#' specification; it does not by itself establish a scoring error.
#'
#' **Separation**
#'
#' Package-reported separation is the ratio of adjusted true standard deviation
#' to root-mean-square measurement error:
#'
#' \deqn{G = \frac{\mathrm{SD}_{\mathrm{adj}}}{\mathrm{RMSE}}}
#'
#' where \eqn{\mathrm{SD}_{\mathrm{adj}} =
#' \sqrt{\mathrm{ObservedVariance} - \mathrm{ErrorVariance}}}. Higher values
#' indicate the facet discriminates more statistically distinct levels along the
#' measured variable. In `mfrmr`, `Separation` is the model-based value and
#' `RealSeparation` provides a more conservative companion based on `RealSE`.
#'
#' **Reliability**
#'
#' \deqn{R = \frac{G^2}{1 + G^2}}
#'
#' Analogous to Cronbach's alpha or KR-20 for the reproducibility of element
#' ordering. In `mfrmr`, `Reliability` is the model-based value and
#' `RealReliability` gives the conservative companion based on `RealSE`. For
#' `MML`, these are anchored to observed-information `ModelSE`
#' estimates for non-person facets; `JML` keeps them as exploratory summaries.
#'
#' The facet's role determines the interpretation. High Person separation
#' reliability indicates distinguishable abilities. High rater separation
#' reliability indicates distinguishable severities, **not agreement** between
#' raters. When raters are intended to be interchangeable, reliable severity
#' differences merit review rather than a high-reliability quality label.
#' Low rater separation reliability can also result from imprecise estimates;
#' it does not establish agreement or equivalence. The index neither tests
#' every rater pair nor quantifies the practical size of severity differences.
#' Its Person interpretation is analogous to test reliability, but it is not
#' numerically identical to Cronbach's alpha.
#'
#' For the person facet under `MML`, the same \eqn{G} and \eqn{R} formulas
#' are applied to EAP person measures with posterior SDs in the error slot.
#' EAP measures are shrunken, so their observed variance is already deflated
#' (approximately the true variance times the reliability), and subtracting
#' the mean posterior variance deflates it again. The reported MML person
#' separation/reliability is therefore a conservative summary: it is
#' systematically lower than the IRT empirical-reliability convention
#' \eqn{\mathrm{Var}(\mathrm{EAP}) / (\mathrm{Var}(\mathrm{EAP}) +
#' \overline{\mathrm{PSD}^2})} and is not numerically comparable to
#' JMLE-based person separation reliability from FACETS. The gap is small
#' when measurement is precise and grows as precision drops. Person rows can
#' still carry the model-based precision tier because posterior SDs are
#' model-based quantities; that tier describes the SE source, not FACETS
#' comparability. Use `method = "JML"` when a FACETS-style person separation
#' table is required, and treat MML person rows as conservative summaries.
#'
#' This is a Rasch/FACETS-style separation reliability on the fitted logit
#' scale, not an intra-class correlation. Use [compute_facet_icc()] only when
#' you want the complementary random-effects variance-share view on the
#' observed-score scale; for non-person facets, large ICC values indicate
#' systematic facet variance rather than desirable measurement reliability.
#'
#' **Strata**
#'
#' Number of statistically distinguishable groups of elements:
#'
#' \deqn{H = \frac{4G + 1}{3}}
#'
#' Three or more strata are commonly used as a practical target
#' (Wright & Masters, 1982), but in this package the estimate inherits the
#' same approximation limits as the separation index.
#'
#' @section Key references:
#' - Andrich, D. (1978). A rating formulation for ordered response
#'   categories. *Psychometrika*, 43, 561--573.
#' - Bond, T. G., & Fox, C. M. (2015). *Applying the Rasch model* (3rd
#'   ed.). Routledge.
#' - Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation
#'   of item parameters: Application of an EM algorithm. *Psychometrika*, 46,
#'   443--459.
#' - Burnham, K. P., & Anderson, D. R. (2002). *Model selection and
#'   multimodel inference: A practical information-theoretic approach*
#'   (2nd ed.). Springer. (AIC / BIC weights and Delta-IC bands used
#'   by `compare_mfrm()`.)
#' - Drasgow, F., Levine, M. V., & Williams, E. A. (1985).
#'   Appropriateness measurement with polychotomous item response
#'   models and standardized indices. *British Journal of Mathematical
#'   and Statistical Psychology*, 38(1), 67--86. (Source for the `lz`
#'   person-fit statistic implemented in `compute_person_fit_indices()`.)
#' - Haberman, S. J., & Sinharay, S. (2013). Generalized residuals for general
#'   models for contingency tables with application to item response theory.
#'   *Journal of the American Statistical Association*, 108, 1435--1444.
#' - Eckes, T. (2005). Examining rater effects in TestDaF writing and speaking
#'   performance assessments: A many-facet Rasch analysis. *Language Assessment
#'   Quarterly*, 2, 197--221.
#' - Muraki, E. (1992). A generalized partial credit model: Application of
#'   an EM algorithm. *Applied Psychological Measurement*, 16(2),
#'   159--176. (Source for the GPCM response-probability kernel used by the
#'   package's bounded `fit_mfrm(model = "GPCM")` route. The slope-aware
#'   `fair_average_table()` and `estimate_bias()` extensions are package-specific,
#'   caveated adaptations built from that kernel, not procedures proposed in
#'   Muraki 1992.)
#' - Muraki, E. (1993). Information functions of the generalized partial
#'   credit model. *Applied Psychological Measurement*, 17(4),
#'   351--363. (Companion paper to Muraki 1992 that derives the GPCM
#'   item information identity \eqn{I_j(\theta) = D^2 a_j^2
#'   \mathrm{Var}(T \mid \theta)} via Samejima's (1974) polytomous
#'   information formula. This is the canonical reference for
#'   `compute_information()` under `GPCM`.)
#' - Samejima, F. (1974). Normal ogive model on the continuous response
#'   level in the multidimensional latent space. *Psychometrika*, 39,
#'   111--121. (General polytomous information formula that Muraki
#'   1993 specializes to the GPCM.)
#' - Snijders, T. A. B. (2001). Asymptotic null distribution of person
#'   fit statistics with estimated person parameter. *Psychometrika*,
#'   66(3), 331--342. (Source for the `lz_star` correction in
#'   `compute_person_fit_indices()` when person estimates come from the
#'   JML/fixed-effect route. MML/EAP person scores are left uncorrected
#'   because EAP does not satisfy the Snijders estimating-equation setup.)
#' - Linacre, J. M. (1989). *Many-facet Rasch measurement*. MESA Press.
#' - Linacre, J. M. (2002). What do Infit and Outfit, mean-square and
#'   standardized mean? *Rasch Measurement Transactions*, 16(2), 878.
#' - Masters, G. N. (1982). A Rasch model for partial credit scoring.
#'   *Psychometrika*, 47, 149--174.
#' - Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for
#'   dichotomous item response theory models. *Applied Psychological
#'   Measurement*, 24, 50--64.
#' - Orlando, M., & Thissen, D. (2003). Further investigation of the
#'   performance of S-X2: An item fit index for use with dichotomous item
#'   response theory models. *Applied Psychological Measurement*, 27, 289--298.
#' - Sinharay, S., Johnson, M. S., & Stern, H. S. (2006). Posterior predictive
#'   assessment of item response theory models. *Applied Psychological
#'   Measurement*, 30, 298--321.
#' - Sinharay, S., & Monroe, S. (2025). Assessment of fit of item response
#'   theory models: A critical review of the status quo and some future
#'   directions. *British Journal of Mathematical and Statistical Psychology*,
#'   78, 711--733.
#' - Wright, B. D., & Masters, G. N. (1982). *Rating scale analysis*.
#'   MESA Press.
#' - Wright, B. D., & Linacre, J. M. (1994). Reasonable mean-square fit
#'   values. *Rasch Measurement Transactions*, 8(3), 370.
#' - Wilson, E. B., & Hilferty, M. M. (1931). The distribution of
#'   chi-square. *Proceedings of the National Academy of Sciences of the
#'   United States of America*, 17(12), 684-688.
#'
#' @section Model selection:
#' **RSM vs PCM**
#'
#' The Rating Scale Model (RSM; Andrich, 1978) assumes all levels of the
#' step facet share identical threshold parameters.  The Partial Credit
#' Model (PCM; Masters, 1982) allows each level of the `step_facet` to have
#' its own set of thresholds on the package's shared observed score scale.
#' Use RSM when the rating rubric is identical across all items/criteria;
#' use PCM when category boundaries are expected to vary by item or criterion.
#' In the current implementation, PCM assumes one common observed score
#' support across the fitted data, so it should not be described as a fully
#' mixed-category model with arbitrary item-specific category counts.
#'
#' **MML, CML and JML**
#'
#' Marginal Maximum Likelihood (MML) integrates over the person ability
#' distribution using Gauss-Hermite quadrature and does not directly estimate
#' person parameters; person estimates are computed post-hoc via Expected A
#' Posteriori (EAP).  Joint Maximum Likelihood (JML) estimates all person
#' and facet parameters simultaneously as fixed effects; `"JMLE"` remains a
#' backward-compatible alias.
#'
#' Estimator choice should follow the inferential target and assumptions. MML
#' avoids estimating one fixed effect per person, but its structural inference
#' depends on the specified population distribution and response model. JML
#' does not assume a parametric person distribution and can be lighter
#' computationally, but structural-parameter bias can persist as the number of
#' persons increases when the number of ratings or items per person remains
#' fixed. Sample size alone therefore does not determine the preferred route.
#' For supported RSM/PCM structural contrasts, CML instead conditions out Person
#' ability without a population model. It estimates neither abilities nor a
#' population SD; same-fit population-SD standardization requires MML.
#'
#' See [fit_mfrm()] for usage.
#'
#' **Fixed-calibration scoring after fitting**
#'
#' [predict_mfrm_units()] and [sample_mfrm_plausible_values()] score future or
#' partially observed persons on a quadrature grid under the fitted scoring
#' basis. For ordinary `MML` fits, these summaries inherit the fitted marginal
#' calibration directly. For latent-regression `MML` fits, they use the fitted
#' one-dimensional conditional normal population model and therefore require
#' one-row-per-person background data for the scored units when the fitted
#' population model includes covariates. Intercept-only latent-regression fits
#' (`population_formula = ~ 1`) can reconstruct that minimal person table from
#' the scored person IDs. For `JML` fits, `mfrmr` uses the fitted facet and
#' step parameters together with a standard normal reference prior introduced
#' only for the post hoc scoring layer. This is useful for practical
#' fixed-scale scoring, but it should still be described as a limited
#' approximation rather than as full ConQuest-style population modeling.
#'
#' **ConQuest overlap**
#'
#' The package supports a scoped latent-regression `MML` model, but
#' the overlap with ConQuest should still be described conservatively. The
#' documented comparison is binary `RSM` / `PCM`, one latent dimension,
#' exactly one item facet, complete rectangular responses, and exactly one
#' numeric person covariate beyond the intercept in a conditional-normal person
#' population model. This is a scoped file-based comparison, not a claim of
#' broad ConQuest numerical equivalence for polytomous, multidimensional,
#' incomplete, or arbitrary imported designs.
#'
#' @examples
#' # Load the package
#' library(mfrmr)
#'
#' # Load example ratings and look at the first six rows
#' ratings <- load_mfrmr_data("example_operational")
#' head(ratings)
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
#' # Plot the results (Wright map)
#' plot(fit)
#'
#' # Save the summary, then display its tables
#' results <- summary(fit)
#' results$person_overview # One row summarizing person ability estimates
#' results$facet_overview  # One row per facet: number of levels, mean, SD, range
#'
#' # Check the interpretation status and recommended next step
#' results$decision
#' @importFrom dplyr across all_of any_of arrange bind_cols bind_rows
#' @importFrom dplyr case_when coalesce count desc distinct everything filter group_by
#' @importFrom dplyr if_all inner_join lag last left_join mutate n n_distinct pull rename
#' @importFrom dplyr row_number rowwise select slice slice_head slice_tail summarise
#' @importFrom dplyr starts_with summarize transmute ungroup
#' @importFrom tidyr drop_na expand_grid pivot_wider replace_na unite
#' @importFrom tibble as_tibble column_to_rownames tibble
#' @importFrom utils combn data globalVariables head modifyList packageVersion
#' @importFrom utils read.csv tail write.csv
#' @importFrom graphics abline arrows axis barplot dotchart hist image legend lines mtext plot plot.new points rect segments text title
#' @importFrom stats cutree hclust kmeans nlminb optim optimHess optimize qt
#' @importFrom utils zip
#' @importFrom cli cli_progress_bar cli_progress_done cli_progress_update
#' @importFrom rlang .data
#' @importFrom stats ave confint na.omit optim p.adjust pchisq predict pt qnorm rnorm sd setNames uniroot
#'
#' @name mfrmr-package
"_PACKAGE"
