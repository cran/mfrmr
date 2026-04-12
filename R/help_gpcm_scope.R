#' Bounded GPCM Support Matrix
#'
#' @description
#' Public capability map for the current `GPCM` scope in `mfrmr`.
#'
#' Use this helper when you need to answer a practical question quickly:
#' which `GPCM` workflows are formally supported in the current core package,
#' which are available only with explicit caveats, and which helpers remain
#' blocked or deferred.
#'
#' The matrix is intentionally conservative. It is a release-scope statement,
#' not a list of every internal code path that happens to run. If a helper is
#' not yet covered by the current validation boundary, it is listed as
#' `blocked` or `deferred` even when some lower-level components already exist.
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
#' - direct slope-aware simulation-spec generation is supported,
#' - APA writer, broader export bundles, fair-average semantics, and planning /
#'   forecasting helpers remain outside the validated `GPCM` boundary.
#'
#' Why some helpers remain blocked:
#'
#' - fair-average, score-side export, and FACETS compatibility-contract outputs depend on
#'   Rasch-family measure-to-score semantics that are not yet generalized to
#'   the free-discrimination `GPCM` branch;
#' - APA writer, visual summaries, and QC pipelines remain blocked because they
#'   would turn those still-unvalidated score-side semantics into narrative or
#'   pass/fail outputs;
#' - planning and forecasting remain deferred because the current design layer
#'   is still validated only for the role-based `RSM` / `PCM` planner.
#'
#' This boundary is aligned with the package's current validation evidence,
#' including the targeted `GPCM` recovery snapshot and the public-workflow
#' regression tests.
#'
#' @return A data.frame with one row per public helper family and columns:
#' - `Area`
#' - `Helpers`
#' - `Status`
#' - `PrimaryUse`
#' - `Boundary`
#' - `Evidence`
#'
#' @section Typical workflow:
#' 1. Call `gpcm_capability_matrix()` before using `GPCM` in a new workflow.
#' 2. Stay on rows marked `supported` or `supported_with_caveat` for the
#'    current release.
#' 3. Treat `blocked` rows as explicit non-support, not as temporary omissions.
#' 4. Treat `deferred` rows as future-extension targets rather than part of the
#'    current package promise.
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [compute_information()],
#'   [predict_mfrm_units()], [sample_mfrm_plausible_values()],
#'   [reporting_checklist()], [mfrmr_workflow_methods], [mfrmr-package]
#' @examples
#' gpcm_capability_matrix()
#' gpcm_capability_matrix("supported")
#' gpcm_capability_matrix("blocked")
#' @export
gpcm_capability_matrix <- function(status = c("all", "supported", "supported_with_caveat", "blocked", "deferred")) {
  status <- match.arg(status)

  out <- data.frame(
    Area = c(
      "Core fitting and summaries",
      "Exploratory diagnostics and residual follow-up",
      "Fixed-calibration scoring and information",
      "Core curve and category views",
      "Checklist and direct reporting route",
      "Operational misfit casebook",
      "Weighting audit and model-choice review",
      "Operational linking synthesis",
      "Direct simulation-spec generation",
      "APA writer and bundle-style reporting",
      "Fair-average, bias, and compatibility-contract layer",
      "Design planning and forecasting",
      "MCMC and heavy-backend extensions"
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
      "reporting_checklist(); precision_audit_report()",
      "build_misfit_casebook()",
      "compare_mfrm(); build_weighting_audit(); compute_information(); plot_information(); build_summary_table_bundle(); export_summary_appendix()",
      "build_linking_review()",
      "build_mfrm_sim_spec(); extract_mfrm_sim_spec(); simulate_mfrm_data()",
      "build_apa_outputs(); build_visual_summaries(); run_qc_pipeline(); build_mfrm_manifest(); build_mfrm_replay_script(); export_mfrm_bundle()",
      "fair_average_table(); estimate_bias(); facets_parity_report(); facets_output_file_bundle(include = \"score\")",
      "evaluate_mfrm_design(); evaluate_mfrm_diagnostic_screening(); evaluate_mfrm_signal_detection(); predict_mfrm_population()",
      "cpp11 backend promotion; posterior predictive computation; MCMC engine; Docker-based advanced runtime"
    ),
    Status = c(
      "supported",
      "supported_with_caveat",
      "supported",
      "supported",
      "supported_with_caveat",
      "supported_with_caveat",
      "supported_with_caveat",
      "deferred",
      "supported_with_caveat",
      "blocked",
      "blocked",
      "deferred",
      "deferred"
    ),
    PrimaryUse = c(
      "Estimate bounded GPCM models and inspect convergence, steps, and slope summaries.",
      "Screen local misfit, residual structure, and agreement patterns after fitting.",
      "Score new units or review design-weighted precision under the fitted GPCM calibration.",
      "Inspect targeting, category progression, and category-probability behavior under the generalized kernel.",
      "Check which direct tables and plots are draft-ready without invoking the broader APA writer.",
      "Combine residual, strict marginal, unexpected-response, and displacement screens into one review queue.",
      "Review whether bounded GPCM is introducing substantively acceptable discrimination-based reweighting.",
      "Synthesize anchor, drift, and chain evidence into one operational review surface.",
      "Generate or extract slope-aware simulation specifications and sample responses directly from them.",
      "Produce manuscript-draft prose, bundled reporting payloads, or full export bundles.",
      "Use Rasch-family score semantics, interaction-bias workflows, or compatibility-contract outputs.",
      "Evaluate designs, forecast future administrations, or run screening-design studies.",
      "Move beyond the current core-package release boundary."
    ),
    Boundary = c(
      paste(
        "Requires an explicit step facet and currently keeps",
        "`slope_facet == step_facet`; MML direct is the validated default,",
        "and EM/hybrid fall back to direct."
      ),
      paste(
        "Residual-based mean-square and strict-marginal outputs remain",
        "exploratory screening tools because discrimination is free;",
        "the dashboard's fair-average panel stays as an explicit unavailable placeholder."
      ),
      paste(
        "Covers fixed-calibration posterior scoring and information only;",
        "population forecasting is a separate layer and remains out of scope."
      ),
      paste(
        "Limited to the slope-aware probability kernel that is already",
        "generalized for the current bounded GPCM branch."
      ),
      paste(
        "Routes users to supported direct tables and plots, but does not imply",
        "that the manuscript-draft APA writer is available for GPCM."
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
        "Not yet validated for bounded GPCM. Use the underlying anchor-audit,",
        "drift, or chain helpers directly and keep the result outside the",
        "current formal GPCM route."
      ),
      paste(
        "Requires explicit slope-aware specifications and keeps the current",
        "bounded branch's facet-role restrictions."
      ),
      paste(
        "Not yet validated for GPCM because the broader reporting/export stack",
        "would convert still-blocked score-side semantics into narrative,",
        "bundle, or QC claims."
      ),
      paste(
        "Not yet generalized to free-discrimination score semantics or to the",
        "current GPCM compatibility boundary."
      ),
      paste(
        "Still validated only for the role-based RSM/PCM planning layer, not",
        "for the bounded GPCM branch."
      ),
      paste(
        "These are future extensions. They are not required to justify the",
        "current GPCM release scope."
      )
    ),
    Evidence = c(
      "covered by estimation and output-stability checks",
      "covered by diagnostic and marginal-plot checks",
      "covered by scoring and information checks",
      "covered by curve, plot, and information checks",
      "covered by reporting-route checks",
      "covered by misfit-casebook and diagnostic checks",
      "covered by weighting-review and information checks",
      "not yet validated for the bounded GPCM route",
      "covered by slope-aware simulation checks",
      "not yet validated for the bounded GPCM route",
      "not yet validated for free-discrimination score semantics",
      "not yet validated for the bounded GPCM route",
      "future extension"
    ),
    stringsAsFactors = FALSE
  )

  if (!identical(status, "all")) {
    out <- out[out$Status == status, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}
