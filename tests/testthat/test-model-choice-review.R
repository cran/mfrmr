test_that("build_model_choice_review bundles comparison and user guidance", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:8]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  rsm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 18
  ))
  gpcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 5,
    maxit = 18
  ))

  review <- build_model_choice_review(RSM = rsm_fit, GPCM = gpcm_fit)

  expect_s3_class(review, "mfrm_model_choice_review")
  expect_s3_class(review, "mfrm_bundle")
  expect_true(all(c(
    "overview", "comparison", "comparison_table", "comparison_warnings",
    "model_roles",
    "downstream_routes", "report_templates", "route_map", "support_status",
    "weighting_review_status", "weighting_review", "notes", "settings"
  ) %in% names(review)))
  expect_true(isTRUE(review$overview$HasBoundedGPCM[1]))
  expect_identical(review$overview$OperationalReference[1], "RSM")
  expect_identical(review$overview$SensitivityModel[1], "GPCM")
  expect_true(any(review$model_roles$RecommendedRole == "equal_weighting_reference"))
  expect_true(any(review$model_roles$RecommendedRole == "slope_aware_sensitivity"))
  expect_true(any(grepl("slope_facet == step_facet", review$model_roles$ScoreContract, fixed = TRUE)))
  rsm_role <- review$model_roles[review$model_roles$Model == "RSM", , drop = FALSE]
  gpcm_role <- review$model_roles[review$model_roles$Model == "GPCM", , drop = FALSE]
  expect_identical(rsm_role$StepStructure, "shared_ladder")
  expect_identical(rsm_role$StepCoordinates, 3L)
  expect_identical(rsm_role$FreeStepParameters, 2L)
  expect_identical(rsm_role$SlopeStructure, "fixed_equal_at_1")
  expect_identical(rsm_role$SlopeCoordinates, 0L)
  expect_identical(rsm_role$FreeSlopeParameters, 0L)
  expect_identical(gpcm_role$StepStructure, "step_facet_specific_ladders")
  expect_identical(gpcm_role$StepCoordinates, 12L)
  expect_identical(gpcm_role$FreeStepParameters, 8L)
  expect_identical(
    gpcm_role$SlopeStructure,
    "relative_by_slope_facet_geometric_mean_1"
  )
  expect_identical(gpcm_role$SlopeCoordinates, 4L)
  expect_identical(gpcm_role$FreeSlopeParameters, 3L)
  expect_identical(rsm_role$FitReadiness, rsm_fit$summary$FitReadiness)
  expect_identical(gpcm_role$FitReadiness, gpcm_fit$summary$FitReadiness)
  expect_identical(
    gpcm_role$FormalInference,
    if (isTRUE(gpcm_fit$summary$InferenceReady)) "Yes" else "No"
  )
  capability <- gpcm_capability_matrix()
  gpcm_route <- review$downstream_routes[
    review$downstream_routes$Model == "GPCM", , drop = FALSE
  ]
  route_areas <- c(
    FullAPARoute = "APA writer and fit-based export bundles",
    ScoreSideExport = "Score-side scorefile export under bounded GPCM",
    LinkingSynthesis = "Operational linking synthesis",
    RecoveryChecks = "Direct simulation-spec generation and recovery",
    FairAverage = "Fair-average semantics under bounded GPCM (slope-aware)",
    BiasScreening = "Residual-bias screening under bounded GPCM",
    SummaryAppendix = "Checklist and summary-table appendix route"
  )
  for (field in names(route_areas)) {
    idx <- match(route_areas[[field]], capability$Area)
    expect_false(is.na(idx), info = paste("Missing capability area for", field))
    expect_identical(
      gpcm_route[[field]][1],
      capability$Status[idx],
      info = paste("Stale bounded-GPCM downstream status for", field)
    )
  }
  expect_identical(gpcm_route$FullAPARoute[1], "supported_with_caveat")
  expect_identical(gpcm_route$ScoreSideExport[1], "supported_with_caveat")
  expect_identical(gpcm_route$LinkingSynthesis[1], "supported_with_caveat")
  expect_true(any(review$downstream_routes$FairAverage == "supported_with_caveat" &
    review$downstream_routes$Model == "GPCM"))
  expect_false(any(grepl(
    "blocked for bounded GPCM|fit-based bundles remain RSM/PCM only",
    review$route_map$Interpretation,
    fixed = FALSE
  )))
  expect_true(any(grepl("automatic operational-scoring decision", review$key_warnings, fixed = TRUE)))
  expect_true(any(grepl(
    "requires optimizer or convergence review",
    review$comparison_warnings,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "screening/review-only",
    review$comparison_warnings,
    fixed = TRUE
  )))
  expect_true(nrow(review$support_status) > 0)
  expect_false(isTRUE(review$weighting_review_status$Requested[1]))
  expect_false("weighting_audit_status" %in% names(review))
  expect_false("weighting_audit" %in% names(review))
  expect_false("run_weighting_audit" %in% names(review$settings))

  sx <- summary(review)
  expect_s3_class(sx, "summary.mfrm_model_choice_review")
  expect_true(all(c(
    "overview", "comparison_table", "comparison_warnings", "model_roles",
    "downstream_routes",
    "report_templates", "route_map", "weighting_review_status"
  ) %in% names(sx)))
  expect_false("weighting_audit_status" %in% names(sx))
  printed <- capture.output(print(sx))
  expect_true(any(grepl("Model Contracts and Readiness", printed, fixed = TRUE)))
  expect_true(any(grepl("Step parameters (reported values versus independent parameters)",
                        printed, fixed = TRUE)))
  expect_true(any(grepl("Slope parameters (reported values versus independent parameters)",
                        printed, fixed = TRUE)))
  expect_true(any(grepl("Fit readiness", printed, fixed = TRUE)))
  expect_true(any(grepl("Full score contracts and readiness reasons", printed,
                        fixed = TRUE)))
  contract_start <- which(printed == "Model Contracts and Readiness")
  downstream_start <- which(printed == "Downstream Routes")
  expect_length(contract_start, 1L)
  expect_length(downstream_start, 1L)
  contract_lines <- printed[contract_start:(downstream_start - 1L)]
  expect_lte(max(nchar(contract_lines)), 80L)
  expect_true(any(grepl("Comparison Warnings", printed, fixed = TRUE)))
  expect_true(any(grepl("screening/review-only", printed, fixed = TRUE)))

  bundle <- build_summary_table_bundle(review)
  expect_s3_class(bundle, "mfrm_summary_table_bundle")
  expect_identical(bundle$source_class, "mfrm_model_choice_review")
  expect_identical(bundle$summary_class, "summary.mfrm_model_choice_review")
  expect_true(all(c(
    "overview",
    "comparison_table",
    "comparison_warnings",
    "model_roles",
    "downstream_routes",
    "report_templates",
    "route_map",
    "weighting_review_status",
    "support_status"
  ) %in% names(bundle$tables)))

  compact_bundle <- build_summary_table_bundle(sx, appendix_preset = "compact")
  expect_identical(compact_bundle$source_class, "summary.mfrm_model_choice_review")
  expect_true(all(c(
    "overview",
    "comparison_table",
    "comparison_warnings",
    "model_roles",
    "report_templates",
    "weighting_review_status"
  ) %in% names(compact_bundle$tables)))
  expect_false("downstream_routes" %in% names(compact_bundle$tables))

  out_dir <- tempfile("mfrmr-model-choice-appendix-")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  appendix <- export_summary_appendix(
    review,
    output_dir = out_dir,
    prefix = "model_choice",
    include_html = FALSE,
    overwrite = TRUE
  )
  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  warning_path <- appendix$written_files$Path[
    grepl("comparison_warnings", appendix$written_files$Component, fixed = TRUE)
  ]
  expect_length(warning_path, 1L)
  exported_warnings <- utils::read.csv(warning_path, stringsAsFactors = FALSE)
  expect_true(any(grepl(
    "screening/review-only",
    exported_warnings$Warning,
    fixed = TRUE
  )))
})

test_that("build_model_choice_review can attach the detailed weighting review", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:8]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  rsm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 18
  ))
  gpcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 5,
    maxit = 18
  ))

  review <- build_model_choice_review(
    RSM = rsm_fit,
    GPCM = gpcm_fit,
    run_weighting_review = TRUE,
    theta_points = 11,
    top_n = 3
  )

  expect_true(isTRUE(review$weighting_review_status$Requested[1]))
  expect_true(isTRUE(review$weighting_review_status$Available[1]))
  expect_s3_class(review$weighting_review, "mfrm_weighting_review")
})

test_that("build_model_choice_review handles RSM versus PCM without GPCM routes", {
  rsm_fit <- make_toy_fit(method = "JML", model = "RSM", maxit = 12)
  toy <- load_mfrmr_data("example_core")
  pcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    model = "PCM",
    step_facet = "Criterion",
    maxit = 12
  ))

  review <- build_model_choice_review(RSM = rsm_fit, PCM = pcm_fit)

  expect_s3_class(review, "mfrm_model_choice_review")
  expect_false(isTRUE(review$overview$HasBoundedGPCM[1]))
  expect_equal(nrow(review$support_status), 0)
  expect_true(all(review$downstream_routes$FullAPARoute == "supported"))
  pcm_role <- review$model_roles[review$model_roles$Model == "PCM", , drop = FALSE]
  expect_identical(pcm_role$StepStructure, "step_facet_specific_ladders")
  expect_identical(pcm_role$StepCoordinates, 12L)
  expect_identical(pcm_role$FreeStepParameters, 8L)
  expect_identical(pcm_role$SlopeStructure, "fixed_equal_at_1")
  expect_identical(pcm_role$SlopeCoordinates, 0L)
  expect_identical(pcm_role$FreeSlopeParameters, 0L)
  expect_true(any(grepl("score interpretation", review$next_actions, fixed = TRUE)))
  expect_true(any(grepl(
    "limited to converged MML fits",
    review$comparison_warnings,
    fixed = TRUE
  )))
  expect_false(any(review$comparison_table$ICComparable))
  expect_true(all(is.na(review$comparison_table$Delta_AIC)))
  expect_length(review$comparison$preferred, 0L)
})

test_that("PCM versus GPCM JML review remains descriptive and preserves its boundary", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:8]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  pcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    model = "PCM",
    step_facet = "Criterion",
    maxit = 18
  ))
  gpcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    maxit = 18
  ))

  review <- build_model_choice_review(PCM = pcm_fit, GPCM = gpcm_fit)

  expect_true(all(review$comparison_table$Method == "JML"))
  expect_false(any(review$comparison_table$ICComparable))
  expect_true(all(is.na(review$comparison_table$Delta_AIC)))
  expect_true(all(is.na(review$comparison_table$Delta_BIC)))
  expect_true(all(is.na(review$comparison_table$Delta_SABIC)))
  expect_length(review$comparison$preferred, 0L)
  expect_true(any(grepl(
    "limited to converged MML fits",
    review$comparison_warnings,
    fixed = TRUE
  )))
})
