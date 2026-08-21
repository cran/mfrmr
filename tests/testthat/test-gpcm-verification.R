# Verification tests for the bounded-GPCM workflow declared in
# `gpcm_capability_matrix()`. Each test exercises one row of the
# matrix that is `"supported"` or `"supported_with_caveat"` and
# asserts the corresponding helper returns the documented shape.
# `"blocked"` and `"deferred"` rows have negative tests where the
# helper should refuse to run (or run with an explicit caveat).

skip_if_no_lme4 <- function() {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    skip("`lme4` (Suggests) not installed; skipping GPCM verification.")
  }
}

local({
  .toy_gpcm <<- load_mfrmr_data("example_core")
  .gpcm_fit <<- suppressMessages(suppressWarnings(
    fit_mfrm(.toy_gpcm, "Person", c("Rater", "Criterion"), "Score",
             method = "MML", model = "GPCM",
             step_facet = "Criterion",
             slope_facet = "Criterion",
             quad_points = 7L, maxit = 25L)
  ))
})

test_that("GPCM core fit returns a populated mfrm_fit", {
  expect_s3_class(.gpcm_fit, "mfrm_fit")
  expect_identical(as.character(.gpcm_fit$config$model), "GPCM")
  expect_true(nrow(.gpcm_fit$summary) > 0L)
  expect_true(nrow(.gpcm_fit$facets$person) > 0L)
  expect_true(.gpcm_fit$population$estimation_converged)
  expect_false(.gpcm_fit$population$inference_ready)
  expect_identical(
    .gpcm_fit$population$converged_basis,
    "legacy_alias_of_inference_ready"
  )
  expect_identical(
    .gpcm_fit$population$converged,
    .gpcm_fit$population$inference_ready
  )
})

test_that("GPCM print / summary do not error", {
  expect_no_error(invisible(utils::capture.output(print(.gpcm_fit))))
  fit_summary <- summary(.gpcm_fit)
  expect_no_error(invisible(utils::capture.output(print(fit_summary))))
  expect_identical(
    fit_summary$decision$Interpretation,
    "Review before reporting or inference"
  )
  expect_identical(fit_summary$decision$FormalInference, "No")
  expect_match(
    fit_summary$decision$Why,
    "Identifiability evidence is incomplete",
    fixed = TRUE
  )
  expect_match(
    fit_summary$decision$Why,
    "Boundary evidence is incomplete",
    fixed = TRUE
  )

  evidence <- as.data.frame(fit_summary$inference_evidence)
  expect_identical(
    evidence$EvidenceArea,
    c(
      "optimizer_stationarity", "local_estimability",
      "observed_information_curvature", "slope_boundary_screen",
      "fit_readiness"
    )
  )
  expect_identical(
    evidence$State[evidence$EvidenceArea == "optimizer_stationarity"],
    "stationary_candidate"
  )
  expect_identical(
    evidence$State[evidence$EvidenceArea == "local_estimability"],
    "locally_full_rank_sufficient"
  )
  expect_identical(
    evidence$State[
      evidence$EvidenceArea == "observed_information_curvature"
    ],
    "locally_positive_definite_at_recorded_tolerances"
  )
  expect_match(
    evidence$State[evidence$EvidenceArea == "slope_boundary_screen"],
    "^none_certified"
  )
  expect_false(.gpcm_fit$readiness$fit$InferenceReady)
  expect_true(all(!.gpcm_fit$slopes$SEEligible))
  expect_true(all(!.gpcm_fit$slopes$CIEligible))
  expect_identical(fit_summary$slope_overview$SEEligible, 0L)
  expect_identical(fit_summary$slope_overview$CIEligible, 0L)

  console <- utils::capture.output(print(fit_summary))
  expect_true(any(grepl("GPCM-MML inference evidence", console, fixed = TRUE)))
  expect_true(any(grepl(
    "Resolve stored readiness reviews before reporting",
    console,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "Resolve numerical warnings before reporting",
    console,
    fixed = TRUE
  )))
})

test_that("q31 and q41 separate stable local evidence from formal readiness", {
  fit_at <- function(q) {
    suppressMessages(fit_mfrm(
      .toy_gpcm,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      method = "MML",
      model = "GPCM",
      step_facet = "Criterion",
      slope_facet = "Criterion",
      quad_points = q,
      maxit = 100L,
      reltol = 1e-9
    ))
  }
  fits <- list(q31 = fit_at(31L), q41 = fit_at(41L))
  evidence <- lapply(fits, function(fit) {
    as.data.frame(summary(fit)$inference_evidence)
  })
  uncertainty <- lapply(fits, function(fit) {
    covariance <- mfrmr:::compute_mml_parameter_covariance(fit)
    mfrmr:::compute_mml_structural_parameter_se(
      fit,
      covariance = covariance
    )$slopes
  })

  for (name in names(fits)) {
    fit <- fits[[name]]
    table <- evidence[[name]]
    curvature <- table[
      table$EvidenceArea == "observed_information_curvature", , drop = FALSE
    ]
    expect_identical(
      table$State[table$EvidenceArea == "optimizer_stationarity"],
      "stationary_candidate",
      info = name
    )
    expect_identical(
      table$State[table$EvidenceArea == "local_estimability"],
      "locally_full_rank_sufficient",
      info = name
    )
    expect_identical(
      curvature$State,
      "locally_positive_definite_at_recorded_tolerances",
      info = name
    )
    expect_equal(curvature$Rank, curvature$Dimension, info = name)
    expect_true(curvature$MinimumEigenvalue > 0, info = name)
    expect_false(fit$readiness$fit$InferenceReady, info = name)
    expect_match(
      fit$readiness$fit$ReasonCodes,
      "mml_gpcm_slope_boundary_not_evaluated",
      fixed = TRUE,
      info = name
    )
    expect_true(any(is.finite(uncertainty[[name]]$OptimizerSE)), info = name)
    expect_true(all(is.na(uncertainty[[name]]$SE)), info = name)
  }

  # These are regression tolerances for this fixed teaching fixture, not
  # general weak-identification or readiness thresholds.
  expect_lt(max(abs(fits$q41$slopes$Estimate - fits$q31$slopes$Estimate)),
            0.002)
  expect_lt(abs(fits$q41$population$sigma2 - fits$q31$population$sigma2),
            0.01)
  expect_lt(max(abs(
    uncertainty$q41$OptimizerSE - uncertainty$q31$OptimizerSE
  )), 0.001)
})

test_that("GPCM iteration report is an explicitly caveated replay", {
  replay <- suppressWarnings(estimation_iteration_report(
    .gpcm_fit,
    max_iter = 2L
  ))

  expect_s3_class(replay, "mfrm_iteration_report")
  expect_gt(nrow(replay$table), 0L)
  expect_true(any(is.finite(replay$table$Objective)))
  expect_equal(nrow(replay$gpcm_boundary), 1L)
  expect_identical(
    replay$gpcm_boundary$Area[1],
    "Replayed optimization diagnostics under bounded GPCM"
  )
  expect_identical(
    replay$gpcm_boundary$Status[1],
    "supported_with_caveat"
  )
  expect_match(replay$gpcm_boundary$Boundary[1], "replay")
})

test_that("GPCM diagnose_mfrm returns measures with caveat status", {
  diag <- suppressMessages(suppressWarnings(
    diagnose_mfrm(.gpcm_fit, residual_pca = "none", diagnostic_mode = "legacy")
  ))
  expect_true(is.data.frame(diag$measures))
  expect_true(nrow(diag$measures) > 0L)
  source_readiness <- mfrmr_get_readiness_record(.gpcm_fit)
  expect_equal(diag$fit_readiness, source_readiness$fit)
  expect_equal(diag$fit_readiness_components, source_readiness$components)
  expect_equal(diag$fit_readiness_parameters, source_readiness$parameters)
  diag_summary <- summary(diag)
  expect_equal(
    as.data.frame(diag_summary$fit_readiness),
    as.data.frame(source_readiness$fit)
  )
  expect_true(any(grepl(
    "source fit is review and is not inference-ready",
    diag_summary$key_warnings,
    fixed = TRUE
  )))
  expect_identical(
    diag_summary$overview$FairAverage[1],
    "available_direct_only"
  )
  expect_true(any(
    diag_summary$reporting_map$Area == "Fair averages / adjusted-score review" &
      diag_summary$reporting_map$CoveredHere == "direct_only" &
      grepl(
        "fair_average_table(fit, diagnostics = diagnostics)",
        diag_summary$reporting_map$CompanionOutput,
        fixed = TRUE
      )
  ))
  diag_console <- capture.output(print(diag_summary))
  expect_true(any(grepl("Source fit readiness", diag_console, fixed = TRUE)))
  expect_true(any(grepl("review-only", diag_console, fixed = TRUE)))
  expect_true(any(grepl(
    "Fair average: Available via direct table only",
    diag_console,
    fixed = TRUE
  )))
  # The embedded dashboard panel remains disabled under GPCM; the diagnostics
  # object points explicitly to the supported direct fair_average_table().
  if (!is.null(diag$fair_average)) {
    expect_false(diag$fair_average$available)
    expect_identical(diag$fair_average$status, "available_direct_only")
    expect_identical(
      diag$fair_average$direct_route,
      "fair_average_table(fit, diagnostics = diagnostics)"
    )
  }
})

test_that("GPCM compute_information + plot_information work", {
  info <- compute_information(.gpcm_fit)
  expect_true(is.list(info))
  expect_true("tif" %in% names(info))
  p <- plot_information(info, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

test_that("GPCM CCC / pathway / Wright plots return mfrm_plot_data", {
  for (type in c("wright", "pathway", "ccc")) {
    p <- plot(.gpcm_fit, type = type, draw = FALSE)
    expect_s3_class(p, "mfrm_plot_data")
  }
})

test_that("GPCM Wright intervals remain visibly screening-only", {
  p <- suppressWarnings(plot(
    .gpcm_fit, type = "wright", show_ci = TRUE, draw = FALSE
  ))
  locations <- as.data.frame(p$data$locations, stringsAsFactors = FALSE)
  facet_rows <- locations$PlotType == "Facet level" &
    is.finite(locations$CI_Lower) & is.finite(locations$CI_Upper)

  expect_true(any(facet_rows))
  expect_identical(p$data$interpretation_status, "review_only")
  expect_true(all(!locations$SupportsFormalInference[facet_rows]))
  expect_true(all(!locations$CIEligible[facet_rows]))
  expect_true(all(locations$CIUse[facet_rows] == "screening_only"))
  expect_true(all(grepl(
    "screening only", locations$CILabel[facet_rows], fixed = TRUE
  )))
})

test_that("GPCM capability matrix is consistent with the helper", {
  m <- gpcm_capability_matrix()
  expect_true(is.data.frame(m))
  expect_identical(
    names(m),
    c("Area", "Helpers", "Status", "Boundary", "RecommendedRoute")
  )
  expect_true(all(m$Status %in%
                    c("supported", "supported_with_caveat", "blocked", "deferred")))
})

test_that("GPCM APA/QC reporting bundle returns with explicit caveats", {
  diag <- suppressMessages(suppressWarnings(
    diagnose_mfrm(.gpcm_fit, residual_pca = "none", diagnostic_mode = "legacy")
  ))
  apa <- suppressMessages(build_apa_outputs(.gpcm_fit, diag))
  expect_s3_class(apa, "mfrm_apa_outputs")
  expect_true(nrow(apa$gpcm_boundary) > 0)
  expect_true(grepl("Bounded\\s+GPCM note", apa$report_text))

  qc <- run_qc_pipeline(.gpcm_fit, diag)
  expect_s3_class(qc, "mfrm_qc_pipeline")
  expect_true(nrow(qc$gpcm_boundary) > 0)
})

test_that("GPCM public results, APA report, export, and replay stay connected", {
  skip_on_cran()

  results <- mfrm_results(.gpcm_fit, include = "gpcm_review")
  expect_s3_class(results, "mfrm_results")
  connected_sections <- c(
    "diagnostics", "fit_summary", "diagnostics_summary", "fair_average",
    "rating_scale", "reporting_checklist"
  )
  expect_identical(
    results$status$Status[match(connected_sections, results$status$Section)],
    rep("ok", length(connected_sections))
  )

  report <- mfrm_report(results, style = "apa", output = "object")
  expect_s3_class(report, "mfrm_report")
  expect_match(report$markdown, "GPCM helper coverage")
  expect_match(report$markdown, "gpcm_capability_matrix()", fixed = TRUE)

  results_dir <- tempfile("mfrmr_gpcm_results_")
  on.exit(unlink(results_dir, recursive = TRUE), add = TRUE)
  exported_results <- export_mfrm_results(
    results,
    output_dir = results_dir,
    prefix = "gpcm_results",
    include = c("summary", "manifest", "replay"),
    overwrite = TRUE,
    acknowledge_sensitive = TRUE
  )
  expect_s3_class(exported_results, "mfrm_results_export")
  expect_true(all(file.exists(exported_results$written_files$Path)))

  bundle_dir <- tempfile("mfrmr_gpcm_replay_")
  dir.create(bundle_dir)
  on.exit(unlink(bundle_dir, recursive = TRUE), add = TRUE)
  export_mfrm_bundle(
    .gpcm_fit,
    output_dir = bundle_dir,
    prefix = "gpcm_bundle",
    include = c("core_tables", "manifest", "script"),
    data = .toy_gpcm,
    acknowledge_sensitive = TRUE
  )
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(bundle_dir)
  replay <- new.env(parent = asNamespace("mfrmr"))
  .mfrmr_muffle_expected_warnings(
    sys.source("gpcm_bundle_replay.R", envir = replay),
    patterns = c(
      "Optimization convergence review did not produce",
      "This fit-level bundle is an analysis archive"
    )
  )

  expect_s3_class(replay$fit, "mfrm_fit")
  expect_identical(as.character(replay$fit$config$model), "GPCM")
  expect_identical(replay$fit$config$step_facet, "Criterion")
  expect_identical(replay$fit$config$slope_facet, "Criterion")
  expect_identical(
    replay$fit$config$gpcm_mml_identification,
    .gpcm_fit$config$gpcm_mml_identification
  )
  expect_equal(
    replay$fit$summary$LogLik, .gpcm_fit$summary$LogLik,
    tolerance = 1e-6
  )
  slope_columns <- intersect(
    c(
      "SlopeFacet", "OptimizerEstimate",
      "FixedLatentSDOptimizerEstimate", "PopulationSD"
    ),
    intersect(names(.gpcm_fit$slopes), names(replay$fit$slopes))
  )
  expect_true(all(c(
    "SlopeFacet", "OptimizerEstimate",
    "FixedLatentSDOptimizerEstimate", "PopulationSD"
  ) %in% slope_columns))
  expect_gt(nrow(replay$fit$slopes), 0L)
  expect_equal(
    replay$fit$slopes[, slope_columns, drop = FALSE],
    .gpcm_fit$slopes[, slope_columns, drop = FALSE],
    tolerance = 1e-6
  )
  expect_equal(
    replay$fit$population$coefficients,
    .gpcm_fit$population$coefficients,
    tolerance = 1e-6
  )
  expect_equal(
    replay$fit$population$sigma2,
    .gpcm_fit$population$sigma2,
    tolerance = 1e-6
  )
})

test_that("five-category all-maximum persons and raters retain distinct contracts", {
  skip_if_not_installed("lpSolve")
  extreme <- expand.grid(
    Person = sprintf("P%02d", 1:20),
    Rater = paste0("R", 1:4),
    Criterion = paste0("C", 1:4),
    stringsAsFactors = FALSE
  )
  extreme$Score <- (
    as.integer(sub("P", "", extreme$Person)) +
      2L * as.integer(sub("R", "", extreme$Rater)) +
      3L * as.integer(sub("C", "", extreme$Criterion))
  ) %% 5L + 1L
  extreme$Score[extreme$Person == "P01"] <- 5L
  extreme$Score[extreme$Rater == "R1"] <- 5L

  jml <- suppressMessages(suppressWarnings(fit_mfrm(
    extreme, "Person", c("Rater", "Criterion"), "Score",
    method = "JML", model = "GPCM",
    step_facet = "Criterion", slope_facet = "Criterion",
    maxit = 80L, reltol = 1e-9
  )))
  expect_identical(jml$config$gpcm_estimator_family,
                   "unpenalized_identified_jml")
  expect_identical(jml$config$gpcm_statistical_penalty, "none")
  expect_false(jml$config$gpcm_finite_parameter_box)
  expect_identical(
    jml$config$gpcm_extreme_person_policy,
    "extended_real_primary_when_certified"
  )
  jml_person <- jml$facets$person[jml$facets$person$Person == "P01", , drop = FALSE]
  expect_identical(jml_person$PrimaryEstimate, Inf)
  expect_true(is.finite(jml_person$OptimizerEstimate))
  expect_identical(jml_person$OptimizerEstimateUse, "numerical_trace_only")
  expect_identical(jml_person$ParameterStatus, "unbounded_high")

  jml_summary <- summary(jml, include_person = TRUE)
  expect_identical(
    jml_summary$settings_overview$GpcmEstimatorFamily[1],
    "unpenalized_identified_jml"
  )
  expect_false(jml_summary$settings_overview$GpcmFiniteParameterBox[1])
  expect_identical(
    jml_summary$settings_overview$GpcmExtremePersonPolicy[1],
    "extended_real_primary_when_certified"
  )
  expect_identical(jml_summary$facet_support_boundaries$Level, "R1")
  expect_true(jml_summary$facet_support_boundaries$BoundaryConstant)
  rater_recession <- jml_summary$facet_recession_review[
    jml_summary$facet_recession_review$Facet == "Rater", , drop = FALSE
  ]
  expect_identical(rater_recession$Level, paste0("R", 1:4))
  expect_identical(
    rater_recession$BoundaryStatus,
    c("unbounded_low", rep("unbounded_high", 3L))
  )
  expect_true(all(rater_recession$AuditScope == "joint_person_structural"))
  expect_false(any(rater_recession$AuditComplete))
  jml_console <- capture.output(print(jml_summary))
  expect_true(any(grepl("Observed boundary-constant facet support", jml_console)))
  expect_true(any(grepl("Certified JML facet recession directions", jml_console)))
  expect_true(any(grepl("numerical traces, not finite JML maxima", jml_console)))
  expect_true(any(grepl(
    "GPCM estimator: unpenalized_identified_jml",
    jml_console,
    fixed = TRUE
  )))
  jml_fit_console <- capture.output(print(jml))
  expect_true(any(grepl(
    "Statistical penalty: none", jml_fit_console, fixed = TRUE
  )))
  expect_true(any(grepl(
    "Finite parameter box: no", jml_fit_console, fixed = TRUE
  )))
  jml_plot <- suppressWarnings(plot(jml, type = "wright", draw = FALSE))
  expect_identical(
    jml_plot$data$scale_contract$GpcmEstimatorFamily[1],
    "unpenalized_identified_jml"
  )

  mml <- suppressMessages(suppressWarnings(fit_mfrm(
    extreme, "Person", c("Rater", "Criterion"), "Score",
    method = "MML", model = "GPCM",
    step_facet = "Criterion", slope_facet = "Criterion",
    quad_points = 7L, maxit = 80L, reltol = 1e-9
  )))
  mml_person <- mml$facets$person[mml$facets$person$Person == "P01", , drop = FALSE]
  expect_true(is.finite(mml_person$PrimaryEstimate))
  expect_true(is.finite(mml_person$PosteriorSD))
  expect_identical(mml_person$PrimaryEstimateBasis, "posterior_eap")
  expect_identical(
    mml_person$ReasonCodes,
    "mml_extreme_response_prior_regularized"
  )
  expect_identical(mml_person$SourceFitReadiness, "blocked")
  expect_false(mml_person$SourceInferenceReady)
  expect_identical(
    mml_person$EstimateUse,
    "review_only_blocked_source_fit"
  )
  research_only_jml_audits <- c(
    "gpcm_fixed_objective_classification",
    "gpcm_asymptotically_affine_transport",
    "gpcm_boundary_compactification",
    "gpcm_rate_hierarchy",
    "gpcm_lexicographic_limit",
    "gpcm_parameter_path_reachability",
    "gpcm_optimizer_sequence_diagnostic",
    "gpcm_parameter_sequence_flag",
    "gpcm_global_existence",
    "gpcm_exponential_balance",
    "gpcm_terminal_gradient_stability"
  )
  expect_length(
    intersect(names(mml$config$boundary_audit), research_only_jml_audits),
    0L
  )
  mml_summary <- summary(mml, include_person = TRUE)
  expect_identical(mml_summary$facet_support_boundaries$Level, "R1")
  expect_equal(nrow(mml_summary$facet_recession_review), 0L)
  expect_identical(mml_summary$decision$FitReadiness, "blocked")
  expect_identical(mml_summary$decision$FormalInference, "No")
  expect_identical(mml_summary$person_overview$Persons, 20L)
  expect_identical(mml_summary$person_overview$DistributionN, 19L)
  expect_identical(
    mml_summary$person_overview$ReviewExcludedExtremeEAPs,
    1L
  )
  expect_identical(
    mml_summary$person_overview$EstimateUse,
    "review_only_blocked_extreme_eap_excluded"
  )
  expect_lt(
    abs(mml_summary$person_overview$Mean),
    abs(mml_person$PrimaryEstimate)
  )
  expect_true(any(grepl(
    "Person distribution and targeting summaries omit 1",
    mml_summary$key_warnings,
    fixed = TRUE
  )))
  mml_plot <- suppressWarnings(plot(
    mml, type = "wright", draw = FALSE
  ))
  expect_false("P01" %in% mml_plot$data$person$Person)
  expect_identical(mml_plot$data$person_exclusions$Person, "P01")
  expect_identical(
    mml_plot$data$person_exclusions$PlotExclusionReason,
    "prior_regularized_extreme_eap_from_blocked_source_fit"
  )
  expect_identical(mml_plot$data$person_stats$ReviewExcludedN, 1L)
  expect_lt(
    max(abs(mml_plot$data$y_range)),
    abs(mml_person$PrimaryEstimate)
  )
  expect_match(
    mml_plot$data$retention_note,
    "omitted from the plotted scale",
    fixed = TRUE
  )
})
