fit_for_summary_profile <- make_toy_fit()
diagnostics_for_summary_profile <- make_toy_diagnostics(
  fit_for_summary_profile,
  residual_pca = "none"
)

test_that("fit profile stays lightweight and preserves the legacy fit contract", {
  testthat::local_mocked_bindings(
    mfrm_results = function(...) {
      stop("expanded results must not run for profile = 'fit'")
    },
    diagnose_mfrm = function(...) {
      stop("diagnostics must not run for profile = 'fit'")
    },
    .package = "mfrmr"
  )

  default_summary <- summary(fit_for_summary_profile)
  fit_summary <- summary(fit_for_summary_profile, profile = "fit")

  expect_s3_class(default_summary, "summary.mfrm_fit")
  expect_identical(default_summary$profile, "fit")
  expect_identical(default_summary$detail, "full")
  expect_null(default_summary$results)
  expect_equal(default_summary$overview, fit_summary$overview)
  expect_true(all(c(
    "overview", "facet_overview", "person_overview", "step_overview",
    "facet_extremes", "person_high", "person_low"
  ) %in% names(default_summary)))
  default_print <- capture.output(print(default_summary))
  expect_false(any(grepl("Highest person measures", default_print, fixed = TRUE)))
  expect_false(any(grepl("Lowest person measures", default_print, fixed = TRUE)))
  expect_identical(default_summary$provenance$ComputePolicy, "never_fit_profile")
  expect_true(any(
    default_summary$section_status$Section == "diagnostics" &
      default_summary$section_status$Status == "not_requested"
  ))
  expect_error(
    summary(fit_for_summary_profile, include_person = 1),
    "must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("step overview evaluates PCM and GPCM ladders within StepFacet", {
  multi_ladder_fit <- fit_for_summary_profile
  multi_ladder_fit$steps <- tibble::tibble(
    StepFacet = rep(c("Criterion_A", "Criterion_B"), each = 3L),
    Step = rep(c("Step_1", "Step_2", "Step_3"), times = 2L),
    Estimate = c(-1, 0, 1, -2, -0.5, 0.75)
  )

  for (model in c("PCM", "GPCM")) {
    model_fit <- multi_ladder_fit
    model_fit$config$model <- model
    model_fit$summary$Model <- model
    fit_summary <- summary(model_fit, profile = "fit")

    expect_identical(
      as.character(fit_summary$step_overview$StepFacet),
      c("Criterion_A", "Criterion_B"),
      info = model
    )
    expect_identical(
      fit_summary$step_overview$Steps,
      c(3L, 3L),
      info = model
    )
    expect_equal(
      fit_summary$step_overview$Span,
      c(2, 2.75),
      info = model
    )
    expect_true(all(fit_summary$step_overview$Monotonic), info = model)
    expect_false(any(grepl(
      "Step estimates are not monotonic",
      fit_summary$notes,
      fixed = TRUE
    )), info = model)
  }

  multi_ladder_fit$config$model <- "PCM"
  multi_ladder_fit$summary$Model <- "PCM"
  multi_ladder_fit$steps$Estimate[5:6] <- c(0.5, 0.25)
  nonmonotone_summary <- summary(multi_ladder_fit, profile = "fit")
  expect_identical(
    nonmonotone_summary$step_overview$Monotonic,
    c(TRUE, FALSE)
  )
  expect_true(any(grepl(
    "Step estimates are not monotonic",
    nonmonotone_summary$notes,
    fixed = TRUE
  )))

  rsm_summary <- summary(fit_for_summary_profile, profile = "fit")
  expect_equal(nrow(rsm_summary$step_overview), 1L)
  expect_equal(rsm_summary$step_overview$Steps, nrow(fit_for_summary_profile$steps))
})

test_that("brief fit summaries lead with Wright and do not print person IDs", {
  fit_summary <- summary(
    fit_for_summary_profile,
    profile = "fit",
    detail = "brief"
  )
  printed <- capture.output(print(fit_summary))
  person_ids <- as.character(fit_for_summary_profile$facets$person$Person)

  expect_identical(
    as.character(fit_summary$required_visual$Visual),
    c(
      "mfrmr Wright map with SE",
      "FACETS-style Wright map",
      "Infit pathway"
    )
  )
  expect_identical(
    fit_summary$required_visual$Required,
    c(TRUE, FALSE, FALSE)
  )
  expect_true(any(grepl("Visual workflow (in order)", printed, fixed = TRUE)))
  expect_true(any(grepl("Wright map", printed, fixed = TRUE)))
  expect_false(any(grepl("Highest person measures", printed, fixed = TRUE)))
  expect_false(any(grepl("Lowest person measures", printed, fixed = TRUE)))
  expect_false(any(vapply(
    person_ids,
    function(id) any(grepl(id, printed, fixed = TRUE)),
    logical(1)
  )))

  opted_in <- summary(
    fit_for_summary_profile,
    profile = "fit",
    detail = "full",
    include_person = TRUE
  )
  opted_in_print <- capture.output(print(opted_in))
  expect_true(any(grepl("Highest person measures", opted_in_print, fixed = TRUE)))
})

test_that("small terminal gradients remain visible in the printed status", {
  gradient_fit <- fit_for_summary_profile
  gradient_fit$summary$TerminalGradientSupNorm <- 0.000232887
  gradient_fit$summary$TerminalGradientRMS <- 0.000135901
  gradient_fit$summary$GradientReviewTolerance <- 0.0001
  gradient_fit$summary$ConvergenceStatus <- "converged_gradient_review"
  gradient_fit$summary$ConvergenceSeverity <- "review"
  gradient_fit$summary$InferenceReady <- FALSE

  printed <- capture.output(print(summary(
    gradient_fit,
    profile = "fit",
    detail = "brief"
  ), digits = 3))
  printed_text <- gsub("[[:space:]]+", " ", paste(printed, collapse = " "))

  expect_match(
    printed_text,
    "maximum absolute gradient: 0.000233",
    fixed = TRUE
  )
  expect_false(grepl(
    "maximum absolute gradient: 0)",
    printed_text,
    fixed = TRUE
  ))
})

test_that("facets profile reuses diagnostics and keeps explicit analysis boundaries", {
  testthat::local_mocked_bindings(
    diagnose_mfrm = function(...) {
      stop("supplied diagnostics must be reused without recomputation")
    },
    estimate_bias = function(...) stop("estimate_bias() must be explicit"),
    analyze_dif = function(...) stop("analyze_dif() must be explicit"),
    analyze_dff = function(...) stop("analyze_dff() must be explicit"),
    analyze_residual_pca = function(...) {
      stop("analyze_residual_pca() must be explicit")
    },
    detect_anchor_drift = function(...) {
      stop("detect_anchor_drift() must be explicit")
    },
    .package = "mfrmr"
  )

  facets_summary <- suppressWarnings(summary(
    fit_for_summary_profile,
    profile = "facets",
    diagnostics = diagnostics_for_summary_profile,
    include_person = TRUE
  ))

  expect_s3_class(facets_summary, "summary.mfrm_fit")
  expect_identical(facets_summary$profile, "facets")
  expect_identical(facets_summary$detail, "brief")
  expect_s3_class(facets_summary$results, "mfrm_results")
  expect_identical(
    facets_summary$results$diagnostics,
    diagnostics_for_summary_profile
  )
  nested_readiness <- facets_summary$results$readiness[
    facets_summary$results$readiness$Domain != "Plot",
    ,
    drop = FALSE
  ]
  rownames(nested_readiness) <- NULL
  outer_readiness <- facets_summary$readiness
  rownames(outer_readiness) <- NULL
  expect_equal(nested_readiness, outer_readiness)
  expect_identical(
    summary(facets_summary$results)$readiness$Status[
      summary(facets_summary$results)$readiness$Domain == "Reporting"
    ],
    facets_summary$readiness$Status[
      facets_summary$readiness$Domain == "Reporting"
    ]
  )
  expect_identical(facets_summary$provenance$DiagnosticsSource, "explicit")
  expect_identical(
    facets_summary$provenance$IdentityCheck,
    "structural_identity_match"
  )
  expect_match(
    facets_summary$provenance$OrganizationBoundary,
    "not a claim of numerical equivalence",
    fixed = TRUE
  )

  visual <- facets_summary$required_visual
  expect_identical(
    as.character(visual$Visual),
    c(
      "mfrmr Wright map with SE",
      "FACETS-style Wright map",
      "Infit pathway"
    )
  )
  expect_identical(visual$Required, c(TRUE, FALSE, FALSE))
  expect_true(all(c(
    "InterpretationStatus", "InterpretationReady"
  ) %in% names(visual)))
  expect_match(visual$Route[1], "show_ci = TRUE", fixed = TRUE)
  expect_match(visual$Route[1], "top_n = Inf", fixed = TRUE)
  expect_match(visual$Route[2], "renderer = \"facets\"", fixed = TRUE)
  expect_match(visual$Route[2], "show_ci = FALSE", fixed = TRUE)
  expect_false(any(grepl("summary(fit", visual$Route, fixed = TRUE), na.rm = TRUE))
  expect_match(
    paste(visual$Route[2], visual$Detail[2]),
    "category_labels",
    fixed = TRUE
  )
  expect_match(visual$Route[3], "fit_stat = \"Infit\"", fixed = TRUE)
  expect_match(visual$Route[3], "include_person = TRUE", fixed = TRUE)

  boundary_rows <- facets_summary$section_status[
    facets_summary$section_status$Section %in%
      c("bias_dif", "residual_pca", "linking_drift"),
    ,
    drop = FALSE
  ]
  expect_setequal(
    boundary_rows$Section,
    c("bias_dif", "residual_pca", "linking_drift")
  )
  expect_true(all(boundary_rows$Status == "not_computed_by_summary"))
  expect_false(any(c(
    "bias_screen", "linking_review", "residual_pca", "anchor_drift"
  ) %in% names(facets_summary$results$components)))

  threshold_table <- as.data.frame(
    facets_summary$results$components$rating_scale$threshold_table
  )
  expect_gt(nrow(threshold_table), 0L)
  expect_true(all(c(
    "LowerCategory", "UpperCategory", "ThresholdCaveat"
  ) %in% names(threshold_table)))

  printed <- capture.output(print(facets_summary))
  person_ids <- as.character(fit_for_summary_profile$facets$person$Person)
  expect_true(any(grepl(
    "Labeled step transitions",
    printed,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Analyses intentionally not run by summary",
    printed,
    fixed = TRUE
  )))
  expect_true(any(grepl("Bias / DIF", printed, fixed = TRUE)))
  expect_true(any(grepl("Residual PCA", printed, fixed = TRUE)))
  expect_true(any(grepl("Linking / anchor drift", printed, fixed = TRUE)))
  expect_false(any(grepl("not_computed_by_summary", printed, fixed = TRUE)))
  expect_false(any(vapply(
    person_ids,
    function(id) any(grepl(id, printed, fixed = TRUE)),
    logical(1)
  )))
})

test_that("expanded profiles reject diagnostics from a different fitted analysis", {
  mismatched_diagnostics <- diagnostics_for_summary_profile
  mismatched_diagnostics$obs$Person <- as.character(
    mismatched_diagnostics$obs$Person
  )
  mismatched_diagnostics$obs$Person[1] <- "not-the-fitted-person"

  expect_error(
    summary(
      fit_for_summary_profile,
      profile = "facets",
      diagnostics = mismatched_diagnostics
    ),
    "fit/diagnostics mismatch",
    fixed = TRUE
  )
})

test_that("automatic facets diagnostics are computed exactly once", {
  automatic_diagnostics <- suppressWarnings(diagnose_mfrm(
    fit_for_summary_profile,
    residual_pca = "none",
    diagnostic_mode = "both",
    fit_df_method = "both"
  ))
  diagnostic_calls <- 0L
  testthat::local_mocked_bindings(
    diagnose_mfrm = function(...) {
      diagnostic_calls <<- diagnostic_calls + 1L
      automatic_diagnostics
    },
    .package = "mfrmr"
  )

  facets_summary <- suppressWarnings(summary(
    fit_for_summary_profile,
    profile = "facets"
  ))

  expect_identical(diagnostic_calls, 1L)
  expect_identical(facets_summary$provenance$DiagnosticsSource, "computed")
  expect_identical(
    facets_summary$results$diagnostics,
    automatic_diagnostics
  )
})

test_that("compute never exposes unavailable diagnostics without hidden work", {
  testthat::local_mocked_bindings(
    diagnose_mfrm = function(...) {
      stop("diagnostics must not run with compute = 'never'")
    },
    .package = "mfrmr"
  )

  facets_summary <- suppressWarnings(summary(
    fit_for_summary_profile,
    profile = "facets",
    compute = "never"
  ))

  expect_s3_class(facets_summary$results, "mfrm_results")
  expect_null(facets_summary$results$diagnostics)
  expect_identical(facets_summary$provenance$ComputePolicy, "never")
  expect_identical(facets_summary$provenance$DiagnosticsSource, "not_computed")
  expect_true(any(
    facets_summary$section_status$Section == "diagnostics" &
      facets_summary$section_status$Status == "not_computed"
  ))
  expect_false(facets_summary$required_visual$Available[3])
  expect_true(is.na(facets_summary$required_visual$Route[3]))
  expect_false(any(
    grepl("summary(fit", facets_summary$required_visual$Route, fixed = TRUE),
    na.rm = TRUE
  ))
  expect_true(any(
    facets_summary$analysis$triage$Area == "Section availability" &
      facets_summary$analysis$triage$Severity == "review"
  ))
  dependent_sections <- c(
    "diagnostics", "diagnostics_summary", "fit_measures",
    "facet_statistics", "fair_average", "rating_scale", "unexpected",
    "facets_fit_review"
  )
  dependent_status <- facets_summary$section_status[
    facets_summary$section_status$Section %in% dependent_sections,
    c("Section", "Status"),
    drop = FALSE
  ]
  expect_setequal(dependent_status$Section, dependent_sections)
  expect_true(all(dependent_status$Status == "not_computed"))
  expect_true(any(
    facets_summary$analysis$triage$Area == "Diagnostics" &
      facets_summary$analysis$triage$Signal == "diagnostics_not_computed"
  ))
  results_overview <- summary(facets_summary$results)$overview
  expect_gt(results_overview$NotComputed, 0L)
  expect_match(
    facets_summary$results$input$reproducible_code,
    "compute = \"never\"",
    fixed = TRUE
  )
})

test_that("data-frame compute never fits without running diagnostic workflow", {
  fit_calls <- 0L
  testthat::local_mocked_bindings(
    fit_mfrm = function(...) {
      fit_calls <<- fit_calls + 1L
      fit_for_summary_profile
    },
    run_mfrm_facets = function(...) {
      stop("run_mfrm_facets() must not run with compute = 'never'")
    },
    diagnose_mfrm = function(...) {
      stop("diagnose_mfrm() must not run with compute = 'never'")
    },
    .package = "mfrmr"
  )

  data_result <- suppressWarnings(mfrm_results(
    load_mfrmr_data("example_core")[, c(
      "Person", "Rater", "Criterion", "Score"
    )],
    include = "fit",
    compute = "never"
  ))

  expect_identical(fit_calls, 1L)
  expect_s3_class(data_result$fit, "mfrm_fit")
  expect_null(data_result$diagnostics)
  expect_identical(data_result$input$mode, "data.frame")
  expect_match(
    data_result$input$reproducible_code,
    "compute = \"never\"",
    fixed = TRUE
  )
})

test_that("replay code preserves diagnostic reuse and Wright uses it", {
  result <- suppressWarnings(mfrm_results(
    fit_for_summary_profile,
    include = "facets",
    diagnostics = diagnostics_for_summary_profile,
    compute = "auto"
  ))

  expect_match(
    result$input$reproducible_code,
    "diagnostics = diagnostics",
    fixed = TRUE
  )
  expect_match(
    result$input$reproducible_code,
    "compute = \"auto\"",
    fixed = TRUE
  )

  expect_warning(
    wright <- plot(result, type = "wright", show_ci = TRUE, draw = FALSE),
    "Review-only display",
    fixed = TRUE
  )
  locations <- as.data.frame(wright$data$locations, stringsAsFactors = FALSE)
  expect_true(any(
    locations$Measure_Source %in% "diagnostics$measures",
    na.rm = TRUE
  ))
})

test_that("not-computed policy remains explicit in downstream reports", {
  result <- suppressWarnings(mfrm_results(
    fit_for_summary_profile,
    include = "standard",
    compute = "never"
  ))
  report <- mfrm_report(result, style = "qc")

  section_status <- setNames(report$sections$Status, report$sections$Section)
  expect_identical(
    unname(section_status[c(
      "First-screen diagnostics",
      "Fit, separation, and precision",
      "Category functioning"
    )]),
    rep("not_computed", 3L)
  )
  expect_true(any(report$first_screen$Status == "not_computed"))
  report_summary <- summary(report)
  expect_identical(report_summary$overview$OverallStatus, "not_computed")
  expect_gt(report_summary$overview$NotComputedAreas, 0L)
})

test_that("reporting profile returns the standard result organization", {
  reporting_summary <- suppressWarnings(summary(
    fit_for_summary_profile,
    profile = "reporting",
    diagnostics = diagnostics_for_summary_profile
  ))

  expect_s3_class(reporting_summary$results, "mfrm_results")
  expect_identical(reporting_summary$profile, "reporting")
  expect_identical(reporting_summary$detail, "brief")
  expect_true("reporting" %in% reporting_summary$results$include)
  expect_false("facets_fit" %in% reporting_summary$results$include)
  expect_identical(reporting_summary$provenance$DiagnosticsSource, "explicit")
})
