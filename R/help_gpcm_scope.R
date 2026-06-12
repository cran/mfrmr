#' Bounded GPCM Support Matrix
#'
#' @description
#' Public capability map for the current `GPCM` scope in `mfrmr`.
#'
#' Use this helper when you need to answer a practical question quickly:
#' which `GPCM` workflows are supported in this release, which are available
#' only with explicit caveats, and which helpers remain blocked or deferred,
#' plus the route to use instead when the requested helper is outside the
#' current boundary.
#'
#' The matrix is intentionally conservative. It is a release-scope statement,
#' not a promise that every lower-level helper can be combined with `GPCM`.
#' If a helper is not yet covered by the current validation boundary, it is
#' listed as `blocked` or `deferred` even when related components already
#' exist.
#'
#' @param status Which rows to return: `"all"` (default), `"supported"`,
#'   `"supported_with_caveat"`, `"blocked"`, or `"deferred"`.
#'
#' @details
#' The current release treats `GPCM` as a bounded supported scope inside the
#' core R package:
#'
#' - fitting and core summaries are supported,
#' - posterior-scoring and information helpers are supported,
#' - residual-based diagnostics and strict marginal follow-up are supported as
#'   exploratory screens,
#' - direct slope-aware simulation-spec generation and parameter-recovery
#'   simulation are supported with caveats,
#' - `fair_average_table()` is supported with an explicit slope-aware
#'   element-conditional caveat,
#' - `estimate_bias()` is supported as conditional screening evidence with
#'   slope-aware information and profile-likelihood follow-up columns,
#' - summary-table appendix export is available for supported direct outputs,
#' - APA writer, visual summaries, QC pipelines, manifests, replay scripts,
#'   and fit-based export bundles are available only as caveated
#'   sensitivity-reporting surfaces with an explicit `gpcm_boundary`,
#' - package-native scorefile export is available with score-side caveats,
#' - role-based design evaluation and population forecasting are available as
#'   caveated bounded-`GPCM` sensitivity evidence,
#' - role-based diagnostic and signal-detection design screening helpers are
#'   available as caveated bounded-`GPCM` sensitivity evidence,
#' - full FACETS output-contract score-side review remains outside the
#'   validated `GPCM` boundary.
#'
#' Why some helpers remain blocked:
#'
#' - full FACETS output-contract score-side review depends on Rasch-family
#'   measure-to-score semantics plus delta-method SE machinery that are not
#'   yet generalized to the free-discrimination `GPCM` branch;
#' - APA writer, fit-based report/export bundles, visual summaries, and QC
#'   pipelines stay caveated because they must not turn unsupported score-side
#'   semantics into narrative or pass/fail outputs;
#' - diagnostic, signal-detection, design-forecast, and linking helpers stay
#'   caveated because their simulation/refit summaries must not become
#'   operational screening, scoring, or arbitrary-facet planning claims.
#'
#' This boundary is aligned with the package's current validation evidence,
#' including the targeted `GPCM` recovery snapshot and the public workflow
#' checks.
#'
#' @return A data.frame with one row per public helper family and columns:
#' - `Area`
#' - `Helpers`
#' - `Status`
#' - `PrimaryUse`
#' - `Boundary`
#' - `Evidence`
#' - `RecommendedRoute`
#' - `NextValidationStep`
#'
#' @section Typical workflow:
#' 1. Call `gpcm_capability_matrix()` before using `GPCM` in a new workflow.
#' 2. Stay on rows marked `supported` or `supported_with_caveat` for the
#'    current release.
#' 3. For `blocked` and `deferred` rows, read `RecommendedRoute` before choosing
#'    a substitute workflow.
#' 4. Treat `blocked` rows as explicit non-support, not as temporary omissions.
#' 5. Treat `deferred` rows as future-extension targets rather than part of the
#'    current user-facing support.
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [compute_information()],
#'   [predict_mfrm_units()], [sample_mfrm_plausible_values()],
#'   [reporting_checklist()], [mfrmr_workflow_methods], [mfrmr-package]
#' @examples
#' gpcm_capability_matrix()
#' gpcm_capability_matrix("supported")
#' gpcm_capability_matrix("blocked")
#' @concept GPCM boundaries
#' @concept route selection
#' @export
gpcm_capability_matrix <- function(status = c("all", "supported", "supported_with_caveat", "blocked", "deferred")) {
  status <- match.arg(status)

  out <- data.frame(
    Area = c(
      "Core fitting and summaries",
      "Exploratory diagnostics and residual follow-up",
      "Fixed-calibration scoring and information",
      "Core curve and category views",
      "Checklist and summary-table appendix route",
      "Operational misfit casebook",
      "Weighting review and model-choice review",
      "Operational linking synthesis",
      "Direct simulation-spec generation and recovery",
      "APA writer and fit-based export bundles",
      "Fair-average semantics under bounded GPCM (slope-aware)",
      "Design evaluation and population forecasting under bounded GPCM",
      "Diagnostic and signal-detection design screening under bounded GPCM",
      "Differential facet functioning screening under bounded GPCM",
      "MCMC and heavy-backend extensions",
      "Residual-bias screening under bounded GPCM",
      "Score-side scorefile export under bounded GPCM",
      "FACETS output-contract score-side review"
    ),
    Helpers = c(
      "fit_mfrm(model = \"GPCM\"); summary(); print()",
      paste(
        "diagnose_mfrm(); analyze_residual_pca(); unexpected_response_table();",
        "displacement_table(); measurable_summary_table();",
        "rating_scale_table(); interrater_agreement_table();",
        "facet_quality_dashboard(); plot_qc_dashboard();",
        "plot_marginal_fit(); plot_marginal_pairwise()"
      ),
      "predict_mfrm_units(); sample_mfrm_plausible_values(); compute_information(); plot_information()",
      paste(
        "plot(fit, type = c(\"wright\", \"pathway\", \"ccc\", \"ccc_surface\"));",
        "category_structure_report(); category_curves_report();",
        "facets_output_file_bundle(include = \"graph\")"
      ),
      "reporting_checklist(); precision_review_report(); build_summary_table_bundle(); export_summary_appendix()",
      "build_misfit_casebook()",
      "compare_mfrm(); build_model_choice_review(); build_weighting_review(); compute_information(); plot_information(); build_summary_table_bundle(); export_summary_appendix()",
      "build_linking_review()",
      "build_mfrm_sim_spec(); extract_mfrm_sim_spec(); simulate_mfrm_data(); evaluate_mfrm_recovery(); assess_mfrm_recovery()",
      "build_apa_outputs(); build_visual_summaries(); run_qc_pipeline(); build_mfrm_manifest(); build_mfrm_replay_script(); export_mfrm_bundle()",
      "fair_average_table()",
      "evaluate_mfrm_design(); predict_mfrm_population()",
      "evaluate_mfrm_diagnostic_screening(); evaluate_mfrm_signal_detection()",
      "analyze_dff(); analyze_dif(); dif_interaction_table(); dif_report(); plot_dif_heatmap(); plot_dif_summary()",
      "cpp11 backend promotion; posterior predictive computation; MCMC engine; Docker-based advanced runtime",
      "estimate_bias()",
      "facets_output_file_bundle(include = \"score\")",
      "facets_output_contract_review()"
    ),
    Status = c(
      "supported",
      "supported_with_caveat",
      "supported",
      "supported",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "deferred",
      "supported_with_caveat",
      "supported_with_caveat",
      "blocked"
    ),
    PrimaryUse = c(
      "Estimate bounded GPCM models and inspect convergence, steps, and slope summaries.",
      "Screen local misfit, residual structure, and agreement patterns after fitting.",
      "Score new units or review design-weighted precision under the fitted GPCM calibration.",
      "Inspect targeting, category progression, and category-probability behavior under the generalized kernel.",
      "Check which direct tables and plots are draft-ready and export their summary tables.",
      "Combine residual, strict marginal, unexpected-response, and displacement screens into one review queue.",
      "Review whether bounded GPCM is introducing substantively acceptable discrimination-based reweighting.",
      "Synthesize anchor, drift, and chain evidence into one exploratory review surface.",
      "Generate or extract slope-aware simulation specifications, sample responses, and run direct parameter-recovery checks.",
      "Produce caveated manuscript-draft prose, fit-based report bundles, manifests, replay scripts, or full export bundles.",
      "Compute slope-aware element-conditional fair-average score adjustments for reporting tables.",
      "Evaluate role-based future designs or forecast one future administration with repeated bounded-GPCM simulation/refit runs.",
      "Run diagnostic-screening or signal-detection operating-characteristic studies.",
      "Review group-by-facet differential-functioning patterns as slope-aware screening evidence.",
      "Move beyond the current core-package release boundary.",
      "Screen residual two-way interaction-bias cells under bounded GPCM at the screening tier.",
      "Export observation-level slope-aware expected score, residual, probability, and slope fields for bounded GPCM.",
      "Run the full FACETS-style output-contract score-side review that depends on validated free-discrimination score metrics."
    ),
    Boundary = c(
      paste(
        "Requires an explicit step facet and currently keeps",
        "`slope_facet == step_facet`; MML direct is the validated default,",
        "and EM/hybrid fall back to direct."
      ),
      paste(
        "Residual-based mean-square and strict-marginal outputs remain",
        "exploratory screening tools because discrimination is free.",
        "The dashboard's fair-average panel reports an explicit",
        "unavailability status when the underlying fit is GPCM."
      ),
      paste(
        "Covers fixed-calibration posterior scoring and information only;",
        "population forecasting is a separate layer outside this row."
      ),
      paste(
        "Limited to the slope-aware probability kernel that is already",
        "generalized for the current bounded GPCM branch."
      ),
      paste(
        "Routes users to supported direct tables and plots. Caveated",
        "manuscript-draft APA and fit-based export bundles are governed by",
        "their separate capability row and carry `gpcm_boundary`."
      ),
      paste(
        "Supported with caveat for bounded GPCM because the casebook inherits",
        "exploratory screening semantics from its underlying sources."
      ),
      paste(
        "Supported with caveat because the helper is an operational review of",
        "Rasch-family equal weighting versus bounded GPCM reweighting, not an automatic model-selection rule."
      ),
      paste(
        "Supported with caveat as an exploratory synthesis over already-built",
        "anchor review, drift, and chain objects. It does not establish an",
        "operational GPCM linking decision, anchor-drift absence claim, or",
        "equating-chain adequacy claim by itself."
      ),
      paste(
        "Requires explicit slope-aware specifications and keeps the current",
        "bounded branch's facet-role restrictions. Recovery checks are direct",
        "simulation/refit summaries, not design-planning or forecasting claims.",
        "`assess_mfrm_recovery()` requires user-supplied practical thresholds",
        "before RMSE or bias can be interpreted as adequate."
      ),
      paste(
        "Supported with caveat as a partial reporting/export bundle over",
        "already-supported GPCM diagnostics, direct tables, plots, manifests,",
        "and replay scripts. Full FACETS-style score-side contract review,",
        "design forecasting, and automatic operational scoring claims remain",
        "outside this route."
      ),
      paste(
        "Slope-aware element-conditional construction: slope-facet element rows",
        "use that level's own slope; non-slope-facet rows (Person, Rater, ...)",
        "use the geometric-mean-one slope by identification convention.",
        "The historical SE columns in the output are scaled facet-measure SEs,",
        "not fair-average SEs. Use `fair_se = TRUE` to request structural",
        "delta-method fair-average SEs for non-person rows when the MML",
        "observed-information Hessian is available."
      ),
      paste(
        "Supported with caveat as a role-based person x rater-like x",
        "criterion-like Monte Carlo simulation/refit route. It uses the",
        "bounded-GPCM generator and refits bounded GPCM with the supplied or",
        "fit-derived step/slope facet contract, but it reports design-level",
        "operating characteristics only. Slope-recovery adequacy, diagnostic",
        "screening operating characteristics, signal detection, and arbitrary-",
        "facet planning remain separate routes."
      ),
      paste(
        "Supported with caveat as slope-aware repeated simulation/refit",
        "screening evidence for the current role-based person x rater-like x",
        "criterion-like design layer. The summaries are Type I proxy,",
        "sensitivity proxy, DIF target-flag, and bias-screening readouts,",
        "not calibrated inferential tests, operational screening gates, or",
        "arbitrary-facet planning validation."
      ),
      paste(
        "Supported with caveat as direct DFF/DIF screening over the fitted",
        "bounded-GPCM expected-score and residual scale. Residual-method",
        "contrasts and interaction cells remain screening evidence; refit",
        "contrasts must retain explicit linking and precision gates before",
        "any stronger subgroup-comparison wording is used."
      ),
      paste(
        "Future extensions, listed for transparency. Out of scope for",
        "the current bounded GPCM branch."
      ),
      paste(
        "Bias point estimates use the slope-aware GPCM kernel: the bias",
        "parameter is the additive shift on the linear predictor that",
        "maximises the per-cell GPCM log-likelihood. `LR ChiSq`,",
        "`LR Prob.`, and profile-CI columns compare that fitted shift",
        "with zero by conditional profile likelihood. SE / t / Prob",
        "columns use conditional plug-in information at the bias point",
        "estimate. All quantities hold theta, steps, slopes, and other",
        "facet estimates fixed, so they support screening and follow-up",
        "review rather than standalone fairness claims."
      ),
      paste(
        "Supported with caveat for package-native scorefile export only.",
        "Rows carry fitted expected score, residual, standardized residual,",
        "observed-category probability, score slope, native structural",
        "expected-score uncertainty, selectable score-side delta SEs, and",
        "explicit caveat fields when the required MML diagnostics are",
        "available. The route does not export FACETS-equivalent score-side",
        "SEs or establish operational score-scale equivalence."
      ),
      paste(
        "Not yet generalized to the full FACETS-style output-contract review.",
        "Direct scorefile export is available with caveats, but contract-wide",
        "coverage and metric claims still require a broader free-discrimination",
        "score-side review contract."
      )
    ),
    Evidence = c(
      "covered by estimation and output-stability checks",
      "covered by diagnostic and marginal-plot checks",
      "covered by scoring and information checks",
      "covered by curve, plot, and information checks",
      "covered by reporting-route and summary-appendix export checks",
      "covered by misfit-casebook and diagnostic checks",
      "covered by weighting-review and information checks",
      "covered by exploratory linking-review guardrail tests",
      "covered by slope-aware simulation and recovery checks",
      "covered by partial-reporting and export-bundle GPCM tests",
      "covered by reduction-to-PCM and worked-example numerical-agreement tests",
      "covered by caveated GPCM design-evaluation and forecast tests",
      "covered by caveated GPCM diagnostic and signal-detection screening tests",
      "covered by caveated GPCM DFF summary, report, and plot-payload tests",
      "future extension",
      "covered by an end-to-end test on a fitted GPCM example",
      "covered by GPCM scorefile export, native uncertainty, and guardrail tests",
      "not yet validated for free-discrimination score semantics"
    ),
    RecommendedRoute = c(
      paste(
        "Use `fit_mfrm(..., model = \"GPCM\", step_facet = ...,",
        "slope_facet = step_facet)` and inspect `summary(fit)`."
      ),
      paste(
        "Use diagnostics as screening evidence and return to direct residual,",
        "unexpected-response, displacement, and category tables before writing claims."
      ),
      paste(
        "Use fixed-calibration scoring and `compute_information()` /",
        "`plot_information()`; keep population forecasting on a separate route."
      ),
      paste(
        "Use draw-free plot objects and category reports for GPCM sensitivity",
        "figures and appendix tables."
      ),
      paste(
        "Export direct supported tables; use the separate APA/QC/export",
        "bundle row for caveated GPCM sensitivity prose and manifests."
      ),
      paste(
        "Use the casebook as a review queue, then confirm flagged rows with",
        "the underlying direct tables."
      ),
      paste(
        "Compare against an equal-weighting `RSM` / `PCM` reference and report",
        "reweighting as sensitivity evidence."
      ),
      paste(
        "Use `build_linking_review()` as a reader-facing index over direct",
        "anchor, drift, or chain outputs; write any GPCM linking language as",
        "exploratory and source-specific."
      ),
      paste(
        "Use ADEMP-style direct recovery checks with explicit practical RMSE,",
        "bias, and uncertainty thresholds."
      ),
      paste(
        "Use the APA/QC/export bundle for caveated GPCM sensitivity reporting;",
        "use package-native scorefile export, design forecasting, and full",
        "FACETS score-side review only through their separate caveated or",
        "blocked rows."
      ),
      paste(
        "Use `fair_average_table(fair_se = TRUE)` when structural fair-average",
        "SEs are required, and label outputs as slope-aware element-conditional."
      ),
      paste(
        "Use `evaluate_mfrm_design(..., model = \"GPCM\", step_facet = ...,",
        "slope_facet = step_facet)` or `predict_mfrm_population()` for",
        "caveated design-level operating-characteristic review; inspect",
        "`gpcm_boundary` and keep slope-recovery adequacy on",
        "`evaluate_mfrm_recovery()`."
      ),
      paste(
        "Use `evaluate_mfrm_diagnostic_screening(..., model = \"GPCM\",",
        "step_facet = ..., slope_facet = step_facet)` or",
        "`evaluate_mfrm_signal_detection()` for caveated slope-aware screening",
        "operating-characteristic review; inspect `gpcm_boundary` and keep",
        "operational screening decisions outside this route."
      ),
      paste(
        "Use `analyze_dff()` / `analyze_dif()` and `dif_interaction_table()`",
        "as screening surfaces, then carry `gpcm_boundary` through",
        "`summary()`, `dif_report()`, `plot_dif_heatmap()`, and",
        "`plot_dif_summary()` before writing claims."
      ),
      paste(
        "Keep this outside the current public GPCM route and track it as",
        "future-extension scope."
      ),
      paste(
        "Use `estimate_bias()` as screening evidence and follow up with explicit",
        "facet-pair review or external validation before fairness language."
      ),
      paste(
        "Use `facets_output_file_bundle(include = \"score\")` for a",
        "package-native bounded-GPCM scorefile with explicit caveat columns;",
        "inspect `gpcm_score_side_contract()`, and do not treat it as",
        "FACETS score-side equivalence."
      ),
      paste(
        "Use direct fair-average tables and graph-only compatibility outputs;",
        "use `gpcm_score_side_contract()` to inspect the unblock criteria,",
        "and keep full FACETS output-contract reviews on the `RSM` / `PCM` route."
      )
    ),
    NextValidationStep = c(
      paste(
        "Add identification and recovery tests before broadening beyond",
        "`slope_facet == step_facet`."
      ),
      paste(
        "Collect free-discrimination diagnostic comparison fixtures before",
        "promoting residual screens to confirmatory wording."
      ),
      paste(
        "Validate latent-regression and population-forecast semantics separately",
        "from fixed-calibration scoring."
      ),
      "Keep kernel-reduction, curve-shape, and draw-free plot-data tests current.",
      paste(
        "Keep direct-output and report-bundle caveats synchronized before",
        "broadening toward operational wording."
      ),
      "Add larger operational case-review examples if external GPCM fixtures become available.",
      paste(
        "Define a model-choice decision policy before turning reweighting review",
        "into recommendation language."
      ),
      paste(
        "Validate larger multi-wave GPCM anchor/drift and equating-chain",
        "fixtures before upgrading exploratory wording to operational linking",
        "claims."
      ),
      paste(
        "Expand multi-seed recovery coverage across slope regimes and sparse",
        "score-category support."
      ),
      paste(
        "Keep partial-section availability tests synchronized with score-side",
        "and design-forecasting guards."
      ),
      "Add external or simulation-backed checks for structural fair-average SEs.",
      paste(
        "Expand multi-seed fixtures across slope regimes, sparse score support,",
        "and fit-derived specifications before using this route for stronger",
        "operational design recommendations."
      ),
      paste(
        "Expand multi-seed fixtures across slope regimes, local-dependence,",
        "step/slope-facet misspecification, sparse score support, and",
        "fit-derived specifications before using this route for stronger",
        "screening recommendations."
      ),
      paste(
        "Add larger subgroup fixtures and simulation operating-characteristic",
        "checks before using DFF rows as fairness, invariance, or bias claims."
      ),
      paste(
        "Decide posterior-predictive, MCMC, and backend scope only after the",
        "score-side contract is stable."
      ),
      paste(
        "Add simulation operating-characteristic or external fixture evidence",
        "before using screening rows as fairness claims."
      ),
      paste(
        "Add unit-slope PCM reduction checks and slope-variation fixtures before",
        "broadening the scorefile route from package-native delta SEs toward",
        "FACETS-style score-side SEs."
      ),
      paste(
        "Complete the `gpcm_score_side_contract()` requirements, including",
        "a FACETS-compatible free-discrimination score metric and output",
        "contract, before enabling full score-side contract review."
      )
    ),
    stringsAsFactors = FALSE
  )

  if (!identical(status, "all")) {
    out <- out[out$Status == status, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

#' Bounded GPCM Score-Side Export Contract
#'
#' @description
#' Minimal contract table for the caveated bounded-`GPCM` scorefile route and
#' the still-blocked full FACETS-style score-side review route.
#'
#' @details
#' This helper does not enable full FACETS-style score-side review. It records
#' the requirements that separate the current caveated
#' `facets_output_file_bundle(include = "score")` route from a future
#' `facets_output_contract_review()` route for bounded `GPCM`.
#'
#' Use it as a release-maintenance checklist. Rows marked
#' `implemented_with_caveat` support the current package-native bounded-`GPCM`
#' scorefile route. Rows marked `required_for_full_facets_review` are still
#' blockers for full FACETS-style output-contract review. Rows marked
#' `validated_dependency` are already available in the package but are not
#' sufficient by themselves to justify full FACETS score-side equivalence.
#'
#' @param status Which rows to return: `"all"` (default),
#'   `"implemented_with_caveat"`, `"required_for_full_facets_review"`, or
#'   `"validated_dependency"`.
#'
#' @return A data.frame with columns:
#' - `ContractArea`
#' - `Requirement`
#' - `CurrentStatus`
#' - `ReleaseBoundary`
#' - `ValidationTarget`
#' - `ExitCriterion`
#'
#' @seealso [gpcm_capability_matrix()], [gpcm_runtime_guard_coverage()],
#'   [facets_output_contract_review()], [facets_output_file_bundle()]
#' @examples
#' gpcm_score_side_contract()
#' gpcm_score_side_contract("implemented_with_caveat")
#' @concept GPCM boundaries
#' @concept FACETS compatibility
#' @export
gpcm_score_side_contract <- function(status = c("all", "implemented_with_caveat", "required_for_full_facets_review", "validated_dependency")) {
  status <- match.arg(status)

  out <- data.frame(
    ContractArea = c(
      "score_estimand",
      "measure_to_score_metric",
      "score_uncertainty",
      "facets_score_uncertainty_contract",
      "structural_fair_average_se",
      "pcm_reduction",
      "export_schema",
      "runtime_guard",
      "release_wording"
    ),
    Requirement = c(
      "Define the bounded-GPCM score-side estimand separately from Rasch-family measure-to-score semantics.",
      "Specify how free-discrimination slopes enter expected-score summaries, residual score-side fields, and caveat columns.",
      "Define native observation-level expected-score uncertainty and selectable score-side delta SEs under free discrimination before exporting bounded-GPCM score files.",
      "Define the FACETS-compatible score-side uncertainty contract before enabling full output-contract review.",
      "Use structural fair-average SEs where available and document when Hessian-based SEs are unavailable.",
      "Preserve unit-slope bounded-GPCM reduction tests against the PCM route before any score-side export is advertised.",
      "Map each scorefile column to a bounded-GPCM source, caveat, or explicit unavailable status.",
      "Keep full FACETS output-contract review blocked until all required_for_full_facets_review rows are satisfied.",
      "Keep sensitivity-model output separate from operational scoring and FACETS equivalence claims."
    ),
    CurrentStatus = c(
      rep("implemented_with_caveat", 3L),
      "required_for_full_facets_review",
      "validated_dependency",
      "validated_dependency",
      "implemented_with_caveat",
      "validated_dependency",
      "implemented_with_caveat"
    ),
    ReleaseBoundary = c(
      "scorefile_supported_with_caveat",
      "scorefile_supported_with_caveat",
      "scorefile_supported_with_caveat",
      "full_facets_review_blocked",
      "available as fair_average_table(fair_se = TRUE), not as scorefile support",
      "available as reduction evidence, not as scorefile support",
      "scorefile_supported_with_caveat",
      "active guard for full FACETS review",
      "scorefile_supported_with_caveat"
    ),
    ValidationTarget = c(
      "A named estimand and interpretation note for every exported bounded-GPCM scorefile quantity.",
      "A deterministic scorefile contract with slope handling and identification conventions.",
      "Native delta-method expected-score SEs and score-side delta SEs where MML diagnostics are available, with explicit not_requested/unavailable status otherwise.",
      "A FACETS-compatible free-discrimination score metric plus uncertainty policy for contract-wide review fields.",
      "Agreement checks that structural fair-average SE columns are present only when supported by the fitted object.",
      "Numerical agreement checks showing bounded-GPCM unit-slope score-side quantities reduce to the PCM route.",
      "A column contract that separates available, caveated, and unavailable bounded-GPCM scorefile fields.",
      "Structured mfrmr_gpcm_scope_error before full FACETS output-contract review work begins.",
      "NEWS, README, help pages, and validation artifacts that prevent operational scoring overclaims."
    ),
    ExitCriterion = c(
      "Scorefile help pages can name the bounded-GPCM estimand without borrowing Rasch-family wording.",
      "Tests cover slope variation, slope_facet identification, expected-score conversion, and boundary categories.",
      "Tests cover finite native expected-score and score-side delta SEs where available, plus explicit not_requested/unavailable status where not available.",
      "facets_output_contract_review() can report bounded-GPCM score-side uncertainty without borrowing Rasch-family SE semantics.",
      "The fair-average SE route remains traceable and does not imply FACETS score-side equivalence.",
      "Unit-slope bounded-GPCM fixtures match PCM score-side outputs within stated tolerance.",
      "facets_output_contract_review() can report bounded-GPCM score rows without silently emitting unsupported fields.",
      "gpcm_runtime_guard_coverage() and score-side helper errors remain synchronized with gpcm_capability_matrix().",
      "Release wording states whether the route is supported, supported_with_caveat, or still blocked."
    ),
    stringsAsFactors = FALSE
  )

  if (!identical(status, "all")) {
    out <- out[out$CurrentStatus == status, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

#' Bounded GPCM Route-Boundary Coverage
#'
#' @description
#' Public table showing how blocked or deferred bounded-`GPCM` capability rows
#' are handled by the current release.
#'
#' @details
#' `gpcm_capability_matrix()` is the user-facing support matrix. This helper
#' records which public helpers stop with `mfrmr_gpcm_scope_error` when called
#' on a bounded `GPCM` path and which capability rows have no public route yet
#' and are therefore documented as future-extension scope.
#'
#' Package checks use this table to keep out-of-scope `GPCM` behavior aligned
#' with the capability matrix. A row with `GuardMode = "runtime_error"` should
#' have `ExpectedConditionClass = "mfrmr_gpcm_scope_error"`. A row with
#' `GuardMode = "roadmap_only"` records a documented future-extension target
#' with no public helper to call in the current release.
#'
#' @return A data.frame with columns:
#' - `Area`
#' - `Helper`
#' - `Status`
#' - `GuardMode`
#' - `ExpectedConditionClass`
#' - `RecommendedRoute`
#' - `NextValidationStep`
#' - `TestRoute`
#' - `Notes`
#'
#' @seealso [gpcm_capability_matrix()], [mfrmr_workflow_methods],
#'   [mfrmr-package]
#' @examples
#' gpcm_runtime_guard_coverage()
#' @export
gpcm_runtime_guard_coverage <- function() {
  matrix <- gpcm_capability_matrix()
  guard <- data.frame(
    Area = c(
      "FACETS output-contract score-side review",
      "MCMC and heavy-backend extensions"
    ),
    Helper = c(
      "facets_output_contract_review()",
      NA_character_
    ),
    GuardMode = c(
      "runtime_error",
      "roadmap_only"
    ),
    ExpectedConditionClass = c(
      "mfrmr_gpcm_scope_error",
      NA_character_
    ),
    TestRoute = c(
      "minimal mfrm_fit",
      "no public runtime helper in 0.2.1"
    ),
    Notes = c(
      "Full FACETS score-side contract review is intentionally unavailable for bounded GPCM; see gpcm_score_side_contract().",
      "Documented as future-extension scope until a public backend/MCMC helper is exposed."
    ),
    stringsAsFactors = FALSE
  )

  idx <- match(guard$Area, matrix$Area)
  guard$Status <- matrix$Status[idx]
  guard$RecommendedRoute <- matrix$RecommendedRoute[idx]
  guard$NextValidationStep <- matrix$NextValidationStep[idx]
  guard <- guard[, c(
    "Area", "Helper", "Status", "GuardMode", "ExpectedConditionClass",
    "RecommendedRoute", "NextValidationStep", "TestRoute", "Notes"
  )]
  rownames(guard) <- NULL
  guard
}
