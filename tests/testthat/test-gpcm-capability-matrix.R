gpcm_runtime_guard_inputs <- function() {
  list(
    fit = structure(
      list(
        config = list(model = "GPCM"),
        summary = data.frame(Model = "GPCM")
      ),
      class = "mfrm_fit"
    ),
    diagnostics = structure(list(), class = "mfrm_diagnostics"),
    sim_spec = structure(list(model = "GPCM"), class = "mfrm_sim_spec"),
    drift = structure(
      list(config = list(models = "GPCM")),
      class = "mfrm_anchor_drift"
    )
  )
}

gpcm_call_runtime_guard <- function(helper, inputs) {
  calls <- list(
    "build_apa_outputs()" = function() {
      build_apa_outputs(inputs$fit, inputs$diagnostics)
    },
    "build_visual_summaries()" = function() {
      build_visual_summaries(inputs$fit, inputs$diagnostics)
    },
    "run_qc_pipeline()" = function() {
      run_qc_pipeline(inputs$fit)
    },
    "build_mfrm_manifest()" = function() {
      build_mfrm_manifest(inputs$fit)
    },
    "build_mfrm_replay_script()" = function() {
      build_mfrm_replay_script(inputs$fit)
    },
    "export_mfrm_bundle()" = function() {
      export_mfrm_bundle(
        inputs$fit,
        output_dir = tempdir(),
        include = "core_tables",
        overwrite = TRUE,
        acknowledge_sensitive = TRUE
      )
    },
    "facets_output_contract_review()" = function() {
      facets_output_contract_review(inputs$fit)
    },
    "facets_output_file_bundle(include = \"score\")" = function() {
      facets_output_file_bundle(inputs$fit, include = "score")
    },
    "evaluate_mfrm_design()" = function() {
      evaluate_mfrm_design(
        design = list(person = 12, rater = 2, criterion = 2, assignment = 2),
        reps = 1,
        model = "GPCM",
        seed = 1
      )
    },
    "evaluate_mfrm_diagnostic_screening()" = function() {
      evaluate_mfrm_diagnostic_screening(
        design = list(person = 12, rater = 2, criterion = 2, assignment = 2),
        reps = 1,
        model = "GPCM",
        seed = 1
      )
    },
    "evaluate_mfrm_signal_detection()" = function() {
      evaluate_mfrm_signal_detection(
        design = list(person = 12, rater = 2, criterion = 2, assignment = 2),
        reps = 1,
        model = "GPCM",
        seed = 1
      )
    },
    "predict_mfrm_population()" = function() {
      predict_mfrm_population(
        sim_spec = inputs$sim_spec,
        reps = 1,
        seed = 1
      )
    },
    "build_linking_review()" = function() {
      build_linking_review(drift = inputs$drift)
    }
  )

  if (!helper %in% names(calls)) {
    stop("No runtime guard input for helper: ", helper, call. = FALSE)
  }
  calls[[helper]]()
}

test_that("gpcm_capability_matrix exposes only user-facing route guidance", {
  tbl <- gpcm_capability_matrix()

  expect_s3_class(tbl, "data.frame")
  expect_identical(
    names(tbl),
    c("Area", "Helpers", "Status", "Boundary", "RecommendedRoute")
  )
  expect_true(all(tbl$Status %in% c(
    "supported", "supported_with_caveat", "blocked", "deferred"
  )))

  expect_true(any(
    tbl$Area == "Core fitting and summaries" &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    tbl$Area == "Exploratory diagnostics and residual follow-up" &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    grepl("build_apa_outputs\\(\\)", tbl$Helpers) &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    tbl$Area == "Design evaluation and population forecasting under bounded GPCM" &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    tbl$Area == "Diagnostic and signal-detection design screening under bounded GPCM" &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    tbl$Area == "Differential facet functioning screening under bounded GPCM" &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    grepl("build_misfit_casebook\\(\\)", tbl$Helpers) &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    grepl("build_weighting_review\\(\\)", tbl$Helpers) &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    grepl("build_linking_review\\(\\)", tbl$Helpers) &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(all(nzchar(tbl$Boundary)))
  expect_true(all(nzchar(tbl$RecommendedRoute)))
  expect_true(any(
    tbl$Area == "APA writer and fit-based export bundles" &
      tbl$Status == "supported_with_caveat" &
      grepl("caveated GPCM sensitivity reporting", tbl$RecommendedRoute, fixed = TRUE) &
      grepl("mfrm_results", tbl$Helpers, fixed = TRUE) &
      grepl("mfrm_report", tbl$Helpers, fixed = TRUE) &
      grepl("export_mfrm_results", tbl$Helpers, fixed = TRUE) &
      grepl("connected public reporting/export route", tbl$Boundary, fixed = TRUE)
  ))
  expect_true(any(
    tbl$Area == "Score-side scorefile export under bounded GPCM" &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    tbl$Area == "FACETS output-contract score-side review" &
      tbl$Status == "blocked"
  ))
  expect_true(any(
    grepl("unexpected_after_bias_table()", tbl$Helpers, fixed = TRUE) &
      tbl$Area == "Residual-bias screening under bounded GPCM" &
      tbl$Status == "supported_with_caveat",
    na.rm = TRUE
  ))
  expect_true(any(
    tbl$Helpers == "estimation_iteration_report()" &
      tbl$Area == "Replayed optimization diagnostics under bounded GPCM" &
      tbl$Status == "supported_with_caveat",
    na.rm = TRUE
  ))
  expect_true(any(
    tbl$Area == "Diagnostic and signal-detection design screening under bounded GPCM" &
      tbl$Status == "supported_with_caveat" &
      grepl("screening evidence", tbl$Boundary, fixed = TRUE)
  ))
})

test_that("fit-specific GPCM boundaries expose owner and primary-slope status", {
  fit <- structure(
    list(
      config = list(
        model = "GPCM",
        slope_facet = "Rater",
        step_facet = "Rater"
      ),
      summary = data.frame(Model = "GPCM"),
      slopes = data.frame(
        SlopeFacet = c("R1", "R2"),
        ParameterStatus = c("not_evaluated", "not_evaluated"),
        SEEligible = c(FALSE, FALSE),
        stringsAsFactors = FALSE
      )
    ),
    class = "mfrm_fit"
  )

  tbl <- gpcm_capability_boundary_table(
    fit,
    area = "Core fitting and summaries"
  )

  expect_equal(nrow(tbl), 1L)
  expect_identical(tbl$SlopeOwner, "Rater")
  expect_identical(tbl$StepOwner, "Rater")
  expect_identical(
    tbl$ModelFamily,
    "aligned_single_owner_relative_slope_gpcm"
  )
  expect_identical(tbl$SlopeAction, "complete_adjacent_predictor")
  expect_identical(tbl$SlopeComposition, "single_owner_relative_gm1")
  expect_identical(tbl$LatentDimensionCount, 1L)
  expect_identical(
    tbl$PrimarySlopeStatus,
    "primary_slope_suppressed_optimizer_trace_only"
  )
  expect_identical(
    tbl$SlopeUncertaintyStatus,
    "not_eligible_per_parameter_readiness"
  )
  expect_identical(
    tbl$OwnerInterpretation,
    "model_conditional_discrimination_not_rater_consistency"
  )
})

test_that("gpcm_capability_matrix filters by status", {
  full_tbl <- gpcm_capability_matrix()
  blocked_tbl <- gpcm_capability_matrix("blocked")
  supported_tbl <- gpcm_capability_matrix("supported")

  expect_true(nrow(blocked_tbl) > 0)
  expect_true(nrow(supported_tbl) > 0)
  expect_true(all(blocked_tbl$Status == "blocked"))
  expect_true(all(supported_tbl$Status == "supported"))
  expect_equal(nrow(blocked_tbl), sum(full_tbl$Status == "blocked"))
  expect_equal(nrow(supported_tbl), sum(full_tbl$Status == "supported"))
})

test_that("gpcm_capability_matrix prints a compact first screen", {
  tbl <- gpcm_capability_matrix()
  expect_s3_class(tbl, "mfrmr_gpcm_capabilities")

  printed <- capture.output(print(tbl))
  expect_match(paste(printed, collapse = "\n"), "workflow availability")
  expect_match(paste(printed, collapse = "\n"), "Route preview")
  expect_match(paste(printed, collapse = "\n"), "Filter by status")
  expect_lt(length(printed), 40L)
})

test_that("gpcm_score_side_contract exposes user-facing availability only", {
  contract <- gpcm_score_side_contract()
  available <- gpcm_score_side_contract("available_with_caveat")
  supporting <- gpcm_score_side_contract("supporting_route")
  unavailable <- gpcm_score_side_contract("unavailable")
  registry <- .gpcm_score_side_contract_registry()
  matrix <- gpcm_capability_matrix()
  score_row <- matrix[matrix$Area == "Score-side scorefile export under bounded GPCM", , drop = FALSE]
  review_row <- matrix[matrix$Area == "FACETS output-contract score-side review", , drop = FALSE]

  expect_s3_class(contract, "data.frame")
  expect_identical(
    names(contract),
    c("Capability", "Status", "Limitation", "Alternative")
  )
  expect_true(nrow(contract) >= 8L)
  expect_true(nrow(available) > 0L)
  expect_true(nrow(supporting) > 0L)
  expect_true(nrow(unavailable) > 0L)
  expect_true(all(available$Status == "available_with_caveat"))
  expect_true(all(supporting$Status == "supporting_route"))
  expect_true(all(unavailable$Status == "unavailable"))
  expect_true(any(contract$Capability == "Package-native score estimand"))
  expect_true(any(contract$Capability == "Native score uncertainty"))
  expect_true(any(contract$Capability == "Full FACETS score-side review" &
                    contract$Status == "unavailable"))
  expect_true(all(nzchar(contract$Limitation)))
  expect_true(all(nzchar(contract$Alternative)))

  expect_identical(
    names(registry),
    c(
      "CompatibilityStatus", "Capability", "PublicStatus", "Limitation",
      "Alternative"
    )
  )

  legacy <- gpcm_score_side_contract("implemented_with_caveat")
  expect_true(nrow(legacy) > 0L)
  expect_true(all(legacy$Status == "available_with_caveat"))
  expect_identical(score_row$Status[1], "supported_with_caveat")
  expect_identical(review_row$Status[1], "blocked")
  expect_true(grepl("gpcm_score_side_contract", score_row$RecommendedRoute[1], fixed = TRUE))
  expect_true(grepl("full FACETS output-contract reviews", review_row$RecommendedRoute[1], fixed = TRUE))
})

test_that("private GPCM route data contain only runtime product fields", {
  public <- gpcm_capability_matrix()
  registry <- .gpcm_capability_registry()
  guard_registry <- .gpcm_runtime_guard_registry()

  expect_identical(
    names(registry),
    c(
      "CapabilityID", "Area", "Helpers", "Status", "Boundary",
      "RecommendedRoute"
    )
  )
  expect_identical(
    names(guard_registry),
    c(
      "CapabilityID", "Area", "Helper", "Status", "AvailabilityMode",
      "ConditionClass", "RecommendedRoute"
    )
  )

  idx <- match(public$Area, registry$Area)
  expect_false(anyNA(idx))
  expect_equal(public$Status, registry$Status[idx])
  expect_equal(public$Boundary, registry$Boundary[idx])
  expect_equal(public$RecommendedRoute, registry$RecommendedRoute[idx])
  expect_identical(anyDuplicated(registry$CapabilityID), 0L)
  expect_true(all(grepl("^[a-z0-9_]+$", registry$CapabilityID)))

  outstanding <- registry[
    registry$Status %in% c("blocked", "deferred"), , drop = FALSE
  ]
  expect_true("CapabilityID" %in% names(guard_registry))
  expect_false(anyNA(match(guard_registry$CapabilityID, registry$CapabilityID)))
  expect_setequal(guard_registry$Area, outstanding$Area)
})

test_that("GPCM unavailable-route table stays user-facing", {
  matrix <- gpcm_capability_matrix()
  coverage <- gpcm_runtime_guard_coverage()
  outstanding <- matrix[matrix$Status %in% c("blocked", "deferred"), , drop = FALSE]

  expect_s3_class(coverage, "data.frame")
  expect_identical(
    names(coverage),
    c("Area", "Helper", "Status", "Boundary", "RecommendedRoute")
  )
  expect_true(nrow(coverage) > 0L)
  expect_equal(sort(unique(coverage$Area)), sort(outstanding$Area))
  expect_true(all(coverage$Status %in% c("blocked", "deferred")))

  idx <- match(coverage$Area, matrix$Area)
  expect_equal(coverage$Status, matrix$Status[idx])
  expect_equal(coverage$Boundary, matrix$Boundary[idx])
  expect_equal(coverage$RecommendedRoute, matrix$RecommendedRoute[idx])
})

test_that("GPCM runtime guards return matrix-synchronized conditions", {
  coverage <- .gpcm_runtime_guard_registry()
  runtime_rows <- coverage[
    coverage$AvailabilityMode == "structured_error", , drop = FALSE
  ]
  inputs <- gpcm_runtime_guard_inputs()

  for (i in seq_len(nrow(runtime_rows))) {
    row <- runtime_rows[i, , drop = FALSE]
    err <- tryCatch(
      gpcm_call_runtime_guard(row$Helper, inputs),
      mfrmr_gpcm_scope_error = function(e) e,
      error = function(e) e
    )

    expect_true(inherits(err, "mfrmr_gpcm_scope_error"),
                info = paste("Expected structured GPCM scope error:", row$Helper))
    expect_identical(err$helper, row$Helper)
    expect_identical(err$area, row$Area)
    expect_identical(err$status, row$Status)
    expect_identical(err$recommended_route, row$RecommendedRoute)
    expect_identical(err$follow_up, row$RecommendedRoute)
    expect_identical(err$next_validation_step, err$follow_up)
  }
})

test_that("GPCM scorefile export is caveated while full contract review stays blocked", {
  gpcm_fit <- structure(
    list(
      config = list(model = "GPCM"),
      summary = data.frame(Model = "GPCM")
    ),
    class = "mfrm_fit"
  )
  gpcm_diag <- structure(list(), class = "mfrm_diagnostics")

  err <- tryCatch(facets_output_contract_review(gpcm_fit), error = function(e) e)
  expect_true(inherits(err, "mfrmr_gpcm_scope_error"))
  expect_true(nzchar(err$helper))
  expect_identical(err$status, "blocked")
  expect_identical(err$area, "FACETS output-contract score-side review")
  expect_true(nzchar(err$recommended_route))
  expect_identical(err$follow_up, err$recommended_route)
  expect_identical(err$next_validation_step, err$follow_up)
  msg <- conditionMessage(err)
  expect_match(msg, "does not support `GPCM` fits", fixed = TRUE)
  expect_match(msg, "GPCM route:", fixed = TRUE)
  expect_match(msg, "Constraint:", fixed = TRUE)
  expect_match(msg, "Recommended route:", fixed = TRUE)
  expect_false(grepl("Next validation step:", msg, fixed = TRUE))
  expect_match(msg, "gpcm_capability_matrix()", fixed = TRUE)
  expect_match(msg, "gpcm_score_side_contract()", fixed = TRUE)
})

test_that("GPCM partial report, QC, export, and linking helpers return caveated objects", {
  skip_on_cran()
  toy <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    model = "GPCM",
    method = "MML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 5,
    maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))

  apa <- build_apa_outputs(fit, diag)
  expect_s3_class(apa, "mfrm_apa_outputs")
  expect_true(nrow(apa$gpcm_boundary) > 0)
  expect_true(any(apa$gpcm_boundary$Area == "APA writer and fit-based export bundles"))
  expect_true(grepl("Bounded\\s+GPCM note", apa$report_text))

  visual <- build_visual_summaries(fit, diag)
  expect_s3_class(visual, "mfrm_visual_summaries")
  expect_true(nrow(visual$gpcm_boundary) > 0)

  qc <- run_qc_pipeline(fit, diag)
  expect_s3_class(qc, "mfrm_qc_pipeline")
  expect_equal(nrow(qc$verdicts), 10)
  expect_true(nrow(qc$gpcm_boundary) > 0)

  manifest <- build_mfrm_manifest(fit, diagnostics = diag)
  expect_s3_class(manifest, "mfrm_manifest")
  expect_true(nrow(manifest$gpcm_boundary) > 0)

  score <- facets_output_file_bundle(fit, diagnostics = diag, include = "score")
  expect_s3_class(score, "mfrm_output_bundle")
  expect_true(all(c("scorefile", "gpcm_score_side_contract", "gpcm_boundary") %in% names(score)))
  expect_true(all(c(
    "ScoreSlope", "ObsProb", "ScoreSideStatus", "ScoreUncertaintyStatus",
    "ScoreUncertaintyMethod", "ExpectedScoreSE", "ExpectedScoreCI_Lower",
    "ExpectedScoreCI_Upper", "ResidualSE", "ScoreSideSE",
    "ScoreSideCI_Lower", "ScoreSideCI_Upper", "ScoreSideSE_Status",
    "ScoreSideSE_Method", "ScoreSideCaveat"
  ) %in% names(score$scorefile)))
  expect_identical(score$settings$score_se_method, "both")
  expect_true(any(score$scorefile$ScoreSideStatus == "supported_with_caveat"))
  expect_false(any(score$scorefile$ScoreUncertaintyStatus == "score_side_se_not_exported"))
  expect_true(any(score$scorefile$ScoreUncertaintyStatus %in% c("ok", "regularized")))
  expect_true(any(is.finite(score$scorefile$ExpectedScoreSE)))
  expect_true(any(score$scorefile$ScoreSideSE_Status == "ok"))
  expect_true(any(is.finite(score$scorefile$ScoreSideSE)))
  ok_score_side <- score$scorefile$ScoreSideSE_Status == "ok"
  expected_score_side_se <-
    abs(as.numeric(score$scorefile$ScoreSlope[ok_score_side])) *
    as.numeric(score$scorefile$Var[ok_score_side]) *
    as.numeric(score$scorefile$ScoreSideLogitSE[ok_score_side])
  expect_equal(
    as.numeric(score$scorefile$ScoreSideSE[ok_score_side]),
    expected_score_side_se,
    tolerance = 1e-3
  )
  expect_true(any(grepl("Native structural delta method",
                        score$scorefile$ScoreUncertaintyMethod,
                        fixed = TRUE)))
  expect_true(any(grepl("Score-side delta method",
                        score$scorefile$ScoreSideSE_Method,
                        fixed = TRUE)))
  expect_true(any(grepl("ScoreSlope x Var",
                        score$scorefile$ScoreSideSE_Method,
                        fixed = TRUE)))
  score_side_only <- facets_output_file_bundle(
    fit,
    diagnostics = diag,
    include = "score",
    score_se_method = "score_side"
  )
  expect_identical(score_side_only$settings$score_se_method, "score_side")
  expect_true(all(score_side_only$scorefile$ScoreUncertaintyStatus == "not_requested"))
  expect_true(any(score_side_only$scorefile$ScoreSideSE_Status == "ok"))
  expect_true(any(is.finite(score_side_only$scorefile$ScoreSideSE)))

  no_score_se <- facets_output_file_bundle(
    fit,
    diagnostics = diag,
    include = "score",
    score_se_method = "none"
  )
  expect_identical(no_score_se$settings$score_se_method, "none")
  expect_true(all(no_score_se$scorefile$ScoreUncertaintyStatus == "not_requested"))
  expect_true(all(no_score_se$scorefile$ScoreSideSE_Status == "not_requested"))
  expect_false(any(is.finite(no_score_se$scorefile$ExpectedScoreSE)))
  expect_false(any(is.finite(no_score_se$scorefile$ScoreSideSE)))

  p_score_se <- plot(score, type = "score_se", draw = FALSE)
  expect_s3_class(p_score_se, "mfrm_plot_data")
  expect_identical(p_score_se$data$se_column, "ScoreSideSE")
  expect_true(any(score$gpcm_boundary$Area == "Score-side scorefile export under bounded GPCM"))

  drift <- structure(
    list(config = list(models = "GPCM")),
    class = "mfrm_anchor_drift"
  )
  link <- build_linking_review(drift = drift)
  expect_s3_class(link, "mfrm_linking_review")
  expect_identical(link$settings$intended_use, "exploratory_gpcm_linking_review")
  expect_true(nrow(link$gpcm_boundary) > 0)

  dff <- analyze_dff(
    fit,
    diag,
    facet = "Criterion",
    group = "Group",
    data = toy,
    method = "residual",
    min_obs = 1
  )
  expect_s3_class(dff, "mfrm_dff")
  expect_true(any(dff$gpcm_boundary$Area == "Differential facet functioning screening under bounded GPCM"))
  dff_summary <- summary(dff)
  expect_true(nrow(dff_summary$gpcm_boundary) > 0)

  refit_calls <- list()
  dff_refit <- local({
    testthat::local_mocked_bindings(
      fit_mfrm = function(...) {
        refit_calls[[length(refit_calls) + 1L]] <<- list(...)
        stop("captured bounded-GPCM subgroup replay")
      },
      .package = "mfrmr"
    )
    analyze_dff(
      fit,
      diag,
      facet = "Criterion",
      group = "Group",
      data = toy,
      method = "refit",
      min_obs = 1
    )
  })
  expect_gt(length(refit_calls), 0L)
  expect_true(all(vapply(
    refit_calls,
    function(args) {
      identical(args$model, "GPCM") &&
        identical(args$step_facet, "Criterion") &&
        identical(args$slope_facet, "Criterion") &&
        isTRUE(all.equal(args$rating_min, fit$config$rating_min)) &&
        isTRUE(all.equal(args$rating_max, fit$config$rating_max)) &&
        identical(args$quad_points, fit$config$replay_inputs$quad_points) &&
        identical(args$maxit, fit$config$replay_inputs$maxit) &&
        identical(args$mml_engine, fit$config$replay_inputs$mml_engine)
    },
    logical(1)
  )))
  expect_identical(
    dff_refit$config$refit_identity_contract,
    "baseline_response_and_estimation_controls_v1"
  )
  expect_identical(dff_refit$config$refit_slope_facet, "Criterion")
  expect_equal(
    dff_refit$config$refit_rating_range,
    c(fit$config$rating_min, fit$config$rating_max)
  )

  dff_report <- dif_report(dff)
  expect_s3_class(dff_report, "mfrm_dif_report")
  expect_true(nrow(dff_report$gpcm_boundary) > 0)
  expect_match(dff_report$narrative, "Bounded GPCM note", fixed = TRUE)
  expect_true(nrow(summary(dff_report)$gpcm_boundary) > 0)

  dff_heatmap <- plot_dif_heatmap(dff, metric = "contrast", draw = FALSE)
  expect_true(nrow(dff_heatmap$data$gpcm_boundary) > 0)
  dff_summary_plot <- plot_dif_summary(dff, draw = FALSE)
  expect_true(nrow(dff_summary_plot$data$gpcm_boundary) > 0)

  dff_interaction <- dif_interaction_table(
    fit,
    diag,
    facet = "Criterion",
    group = "Group",
    data = toy,
    min_obs = 1
  )
  expect_s3_class(dff_interaction, "mfrm_dif_interaction")
  expect_true(nrow(dff_interaction$gpcm_boundary) > 0)
  expect_true(nrow(summary(dff_interaction)$gpcm_boundary) > 0)
  dff_interaction_report <- dif_report(dff_interaction)
  expect_true(nrow(dff_interaction_report$gpcm_boundary) > 0)
})

test_that("GPCM score-side delta SE stays on the expected-score scale", {
  fn <- getFromNamespace("add_gpcm_score_side_delta_se", "mfrmr")
  scorefile <- data.frame(
    Person = "P1",
    Rater = "R1",
    Criterion = "C1",
    Observed = 3,
    Expected = 2,
    Residual = 1,
    ScoreSlope = 1.25,
    Var = 0.40,
    stringsAsFactors = FALSE
  )
  res <- list(
    config = list(model = "GPCM", facet_names = c("Rater", "Criterion")),
    prep = list(rating_min = 1, rating_max = 4)
  )
  diagnostics <- list(
    measures = data.frame(
      Facet = c("Person", "Rater", "Criterion"),
      Level = c("P1", "R1", "C1"),
      ModelSE = c(0.30, 0.40, 0.50),
      stringsAsFactors = FALSE
    )
  )

  out <- fn(scorefile, res = res, diagnostics = diagnostics)
  eta_se <- sqrt(0.30^2 + 0.40^2 + 0.50^2)
  expected_score_se <- 1.25 * 0.40 * eta_se

  expect_equal(out$ScoreSideLogitSE, eta_se, tolerance = 1e-12)
  expect_equal(out$ScoreSideSE, expected_score_se, tolerance = 1e-12)
  expect_equal(out$ScoreSideResidualSE, expected_score_se, tolerance = 1e-12)
  expect_true(grepl("ScoreSlope x Var", out$ScoreSideSE_Method, fixed = TRUE))
})

test_that("GPCM design evaluation and population forecasts return caveated objects", {
  skip_on_cran()
  gpcm_spec <- build_mfrm_sim_spec(
    design = list(person = 10, rater = 2, criterion = 2, assignment = 2),
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    assignment = "rotating"
  )

  design <- suppressWarnings(evaluate_mfrm_design(
    sim_spec = gpcm_spec,
    n_person = gpcm_spec$n_person,
    n_rater = gpcm_spec$n_rater,
    n_criterion = gpcm_spec$n_criterion,
    raters_per_person = gpcm_spec$raters_per_person,
    reps = 1,
    maxit = 20,
    seed = 1,
    progress = FALSE
  ))
  expect_s3_class(design, "mfrm_design_evaluation")
  expect_identical(design$settings$model, "GPCM")
  expect_identical(design$settings$gpcm_design_status, "supported_with_caveat")
  expect_true(nrow(design$gpcm_boundary) > 0L)
  expect_true(any(design$gpcm_boundary$Area == "Design evaluation and population forecasting under bounded GPCM"))
  expect_true(any(design$results$FitModel == "GPCM"))
  expect_true(any(grepl("design-level sensitivity evidence", design$notes, fixed = TRUE)))
  design_summary <- summary(design)
  expect_true(nrow(design_summary$gpcm_boundary) > 0L)
  expect_true(any(design_summary$design_summary$RecoveryComparableRate == 1))

  forecast <- suppressWarnings(predict_mfrm_population(
    sim_spec = gpcm_spec,
    reps = 1,
    maxit = 20,
    seed = 2
  ))
  expect_s3_class(forecast, "mfrm_population_prediction")
  expect_identical(forecast$settings$model, "GPCM")
  expect_identical(forecast$settings$gpcm_design_status, "supported_with_caveat")
  expect_true(nrow(forecast$gpcm_boundary) > 0L)
  expect_true(any(forecast$forecast$RecoveryComparableRate == 1))
  forecast_summary <- summary(forecast)
  expect_true(nrow(forecast_summary$gpcm_boundary) > 0L)
  expect_true(any(grepl("design-level sensitivity evidence", forecast_summary$notes, fixed = TRUE)))
})

test_that("GPCM diagnostic and signal-detection screening return caveated objects", {
  skip_on_cran()
  gpcm_spec <- build_mfrm_sim_spec(
    design = list(person = 10, rater = 2, criterion = 2, assignment = 2),
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    slopes = c(C01 = 0.8, C02 = 1.25),
    assignment = "rotating",
    group_levels = c("A", "B")
  )

  diag_eval <- suppressWarnings(evaluate_mfrm_diagnostic_screening(
    sim_spec = gpcm_spec,
    n_person = gpcm_spec$n_person,
    n_rater = gpcm_spec$n_rater,
    n_criterion = gpcm_spec$n_criterion,
    raters_per_person = gpcm_spec$raters_per_person,
    reps = 1,
    scenarios = c("well_specified", "local_dependence"),
    maxit = 10,
    seed = 3
  ))
  expect_s3_class(diag_eval, "mfrm_diagnostic_screening")
  expect_identical(diag_eval$settings$model, "GPCM")
  expect_identical(diag_eval$settings$gpcm_screening_status, "supported_with_caveat")
  expect_true(nrow(diag_eval$gpcm_boundary) > 0L)
  expect_true(any(diag_eval$results$Model == "GPCM"))
  expect_true(any(grepl("slope-aware operating-characteristic", diag_eval$notes, fixed = TRUE)))
  diag_summary <- summary(diag_eval)
  expect_s3_class(diag_summary, "summary.mfrm_diagnostic_screening")
  expect_true(nrow(diag_summary$gpcm_boundary) > 0L)
  expect_true(any(grepl("slope-aware operating-characteristic", diag_summary$notes, fixed = TRUE)))
  diag_plot <- plot(diag_eval, draw = FALSE)
  expect_s3_class(diag_plot, "mfrm_plot_data")
  expect_true(nrow(diag_plot$data$gpcm_boundary) > 0L)

  signal_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
    sim_spec = gpcm_spec,
    n_person = gpcm_spec$n_person,
    n_rater = gpcm_spec$n_rater,
    n_criterion = gpcm_spec$n_criterion,
    raters_per_person = gpcm_spec$raters_per_person,
    reps = 1,
    fit_method = "MML",
    maxit = 10,
    quad_points = 5,
    bias_max_iter = 1,
    dif_min_obs = 1,
    seed = 4
  ))
  expect_s3_class(signal_eval, "mfrm_signal_detection")
  expect_identical(signal_eval$settings$model, "GPCM")
  expect_identical(signal_eval$settings$gpcm_screening_status, "supported_with_caveat")
  expect_true(nrow(signal_eval$gpcm_boundary) > 0L)
  signal_summary <- summary(signal_eval)
  expect_s3_class(signal_summary, "summary.mfrm_signal_detection")
  expect_true(nrow(signal_summary$gpcm_boundary) > 0L)
  expect_true(any(grepl("slope-aware operating-characteristic", signal_summary$notes, fixed = TRUE)))
  signal_plot <- plot(signal_eval, signal = "dif", metric = "power", draw = FALSE)
  expect_true(is.list(signal_plot))
  expect_true(nrow(signal_plot$gpcm_boundary) > 0L)

  sparse_spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "sparse_linked",
    sparse_controls = list(link_persons = 2, link_raters_per_person = 3),
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    slopes = c(C01 = 0.8, C02 = 1.25),
    group_levels = c("A", "B")
  )
  sparse_diag <- suppressWarnings(evaluate_mfrm_diagnostic_screening(
    sim_spec = sparse_spec,
    n_person = sparse_spec$n_person,
    n_rater = sparse_spec$n_rater,
    n_criterion = sparse_spec$n_criterion,
    raters_per_person = sparse_spec$raters_per_person,
    reps = 1,
    scenarios = "well_specified",
    maxit = 8,
    seed = 5
  ))
  expect_s3_class(sparse_diag, "mfrm_diagnostic_screening")
  expect_identical(sparse_diag$settings$generator_assignment, "sparse_linked")
  expect_true(nrow(sparse_diag$gpcm_boundary) > 0L)

  sparse_signal <- suppressWarnings(evaluate_mfrm_signal_detection(
    sim_spec = sparse_spec,
    n_person = sparse_spec$n_person,
    n_rater = sparse_spec$n_rater,
    n_criterion = sparse_spec$n_criterion,
    raters_per_person = sparse_spec$raters_per_person,
    reps = 1,
    fit_method = "MML",
    maxit = 8,
    quad_points = 5,
    bias_max_iter = 1,
    dif_min_obs = 1,
    seed = 6
  ))
  expect_s3_class(sparse_signal, "mfrm_signal_detection")
  expect_identical(sparse_signal$settings$generator_assignment, "sparse_linked")
  expect_true(nrow(sparse_signal$gpcm_boundary) > 0L)
})

test_that("GPCM scope errors can be caught by class", {
  gpcm_fit <- structure(
    list(
      config = list(model = "GPCM"),
      summary = data.frame(Model = "GPCM")
    ),
    class = "mfrm_fit"
  )

  err <- tryCatch(
    facets_output_contract_review(gpcm_fit),
    mfrmr_gpcm_scope_error = function(e) e
  )

  expect_true(inherits(err, "mfrmr_gpcm_scope_error"))
  expect_identical(err$helper, "facets_output_contract_review()")
  expect_identical(err$area, "FACETS output-contract score-side review")
  expect_identical(err$status, "blocked")
  expect_match(err$recommended_route, "direct fair-average tables", fixed = TRUE)
  expect_identical(err$follow_up, err$recommended_route)
  expect_identical(err$next_validation_step, err$follow_up)
})
