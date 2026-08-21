test_that("CRAN smoke covers the primary MML review and export route", {
  dat <- load_mfrmr_data("example_operational")
  roster_env <- new.env(parent = emptyenv())
  utils::data(
    list = "mfrmr_example_operational_design",
    package = "mfrmr",
    envir = roster_env
  )
  described <- describe_mfrm_data(
    dat,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    expected_design = roster_env$mfrmr_example_operational_design
  )
  expect_s3_class(described, "mfrm_data_description")
  expect_equal(described$structural_missingness$summary$MissingExpectedCells, 6L)
  expect_true(all(described$design_connectivity$Connected))

  fit <- fit_mfrm(
    dat,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM"
  )
  expect_s3_class(fit, "mfrm_fit")
  expect_equal(fit$summary$Method, "MML")
  expect_true(isTRUE(fit$summary$Converged))
  expect_true(isTRUE(fit$summary$InferenceReady))
  expect_identical(fit$summary$ConvergenceSeverity, "pass")
  # This is intentionally the exact README/default path (31 quadrature
  # points, maxit = 400, reltol = 1e-9). A reduced q = 7 smoke fit previously
  # hid the terminal-gradient review reproduced by the public example.
  expect_equal(fit$config$estimation_control$quad_points, 31L)
  expect_equal(fit$summary$RequestedReltol, 1e-9)
  expect_true(isTRUE(fit$summary$OptimizerPolished))
  expect_lt(fit$summary$EffectiveReltol, fit$summary$RequestedReltol)
  stages <- fit$opt$optimizer_polish$Stages
  expect_true(all(c("Reltol", "Factr", "Pgtol") %in% names(stages)))
  initial_stage <- stages[stages$StageLabel == "initial", , drop = FALSE]
  expect_identical(initial_stage$Method, "L-BFGS-B")
  expect_equal(initial_stage$Factr, 1e-9 / .Machine$double.eps)
  expect_equal(initial_stage$Pgtol, sqrt(.Machine$double.eps))
  selected_stage <- stages[stages$Selected, , drop = FALSE]
  expect_equal(fit$summary$OptimizerFactr, selected_stage$Factr)
  expect_equal(fit$summary$OptimizerPgtol, selected_stage$Pgtol)

  diagnostics <- diagnose_mfrm(
    fit,
    residual_pca = "none",
    diagnostic_mode = "both",
    fit_df_method = "both"
  )
  expect_s3_class(diagnostics, "mfrm_diagnostics")

  fit_summary <- summary(fit, profile = "fit", detail = "brief")
  facets_summary <- summary(
    fit,
    profile = "facets",
    detail = "brief",
    diagnostics = diagnostics,
    compute = "never"
  )
  expect_s3_class(fit_summary, "summary.mfrm_fit")
  expect_s3_class(facets_summary, "summary.mfrm_fit")
  expect_s3_class(facets_summary$results, "mfrm_results")
  expect_identical(
    fit_summary$decision$Interpretation,
    "Fit gates passed; formal precision review required"
  )
  expect_identical(fit_summary$decision$FormalInference, "No")
  expect_identical(fit_summary$decision$FitReadiness, "ready")
  expect_identical(facets_summary$decision$Interpretation,
                   "Ready for formal inference")
  expect_identical(facets_summary$decision$FormalInference, "Yes")
  expect_identical(
    fit_summary$decision$NextAction,
    fit_summary$next_actions[1L]
  )

  native <- plot(
    fit,
    type = "wright",
    renderer = "native",
    show_ci = TRUE,
    top_n = Inf,
    draw = FALSE
  )
  score_values <- fit$prep$score_map$OriginalScore
  rubric_labels <- stats::setNames(paste("Rubric", score_values), score_values)
  facets <- plot_wright_unified(
    fit,
    renderer = "facets",
    show_ci = FALSE,
    top_n = Inf,
    category_labels = rubric_labels,
    draw = FALSE
  )
  pathway <- plot(
    fit,
    type = "fit_pathway",
    diagnostics = diagnostics,
    fit_stat = "Infit",
    include_person = TRUE,
    draw = FALSE
  )
  expect_s3_class(native, "mfrm_plot_data")
  expect_type(facets, "list")
  expect_identical(facets$renderer, "facets")
  expect_s3_class(pathway, "mfrm_plot_data")
  expect_identical(
    native$data$interpretation_status,
    "ready_for_diagnostic_interpretation"
  )
  expect_identical(
    facets$interpretation_status,
    "ready_for_diagnostic_interpretation"
  )
  expect_identical(
    pathway$data$interpretation_status,
    "ready_for_diagnostic_interpretation"
  )

  output_dir <- tempfile("mfrmr-cran-smoke-")
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  expect_warning(
    written <- export_mfrm(
      fit,
      diagnostics = diagnostics,
      output_dir = output_dir,
      prefix = "smoke",
      tables = c("facets", "summary", "steps")
    ),
    "not deidentify|identif|sensitive"
  )
  expect_s3_class(written, "data.frame")
  expect_true(all(file.exists(written$Path)))
  expect_true(all(!written$Deidentified))
  expect_true(all(!written$ShareableWithoutReview))
})

test_that("the executable fit_mfrm example remains short and inference-ready", {
  dat <- load_mfrmr_data("example_operational")
  expect_silent(
    fit <- fit_mfrm(
      dat,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      method = "MML",
      model = "RSM",
      quad_points = 7,
      maxit = 30,
      reltol = 1e-11
    )
  )
  expect_true(isTRUE(fit$summary$InferenceReady))
  expect_identical(fit$summary$ConvergenceSeverity, "pass")
})

test_that("the canonical PCM first-use route separates coordinates from freedom", {
  dat <- load_mfrmr_data("example_operational")
  expect_silent(
    fit <- fit_mfrm(
      dat,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      method = "MML",
      model = "PCM",
      step_facet = "Criterion",
      quad_points = 7,
      maxit = 30,
      reltol = 1e-11
    )
  )
  fit_summary <- summary(fit, profile = "fit", detail = "brief")
  sizes <- mfrmr:::build_param_sizes(fit$config)

  expect_identical(fit_summary$decision$Interpretation,
                   "Fit gates passed; formal precision review required")
  expect_identical(fit_summary$decision$FormalInference, "No")
  expect_match(
    fit_summary$decision$NextAction,
    "diagnose_mfrm()",
    fixed = TRUE
  )
  expect_identical(nrow(fit$steps), 9L)
  expect_identical(as.integer(sizes$steps), 6L)
  expect_identical(nrow(fit$slopes), 0L)
  expect_null(sizes$log_slopes)
})
