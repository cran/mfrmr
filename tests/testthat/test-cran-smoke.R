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
    "Fit-readiness requirements satisfied; formal precision review required"
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
    written <- export_mfrm_tables(
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
  # Exercise the default archive route, including absence of optional predictions.
  bundle <- export_mfrm_bundle(
    fit, diagnostics = diagnostics, output_dir = file.path(output_dir, "archive"),
    acknowledge_sensitive = TRUE
  )
  expect_s3_class(bundle, "mfrm_export_bundle")
  expect_true(all(file.exists(bundle$written_files$Path)))
  expect_false(bundle$privacy_notice$Deidentified)
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
                   "Fit-readiness requirements satisfied; formal precision review required")
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

test_that("public adaptive fits use one objective for likelihood, gradient, covariance and EAP", {
  # GPCM exercises free slopes and category-specific steps in the live Hessian route.
  models <- if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) c("RSM", "PCM", "GPCM") else "GPCM"
  for (model in models) {
    fixture <- adaptive_api_fixture(model)
    fit <- fixture$fit
    config <- fit$config
    order <- config$estimation_control$quad_points
    expect_equal(order, if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) 31L else 15L)
    sizes <- mfrmr:::build_param_sizes(config)
    idx <- mfrmr:::build_indices(fit$prep, config$step_facet, config$slope_facet,
                                 config$interaction_specs)
    evaluate <- mfrmr:::mfrmr_make_adaptive_mml_evaluator(idx, config, sizes, order)
    expected <- evaluate(fit$opt$par)
    expect_identical(fit$summary$MMLIntegration, "adaptive")
    displayed <- paste(capture.output(print(summary(fit))), collapse = "\n")
    displayed <- gsub("[[:space:]]+", " ", displayed)
    expect_match(displayed, paste0("Integration: adaptive Gauss-Hermite, q=", order), fixed = TRUE)
    expect_match(fit$summary$IntegrationEvaluationId, paste0("aghq_normal_v1:q=", order), fixed = TRUE)
    expect_lt(abs(fit$summary$LogLik + expected$value), 1e-8)
    expect_lt(abs(fit$summary$TerminalGradientSupNorm - max(abs(expected$gradient))), 1e-8)
    # GPCM covers the general live Hessian route in ordinary checks.
    # The complete suite also verifies the RSM and PCM special cases.
    if (model == "GPCM" || isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
      covariance <- mfrmr:::compute_mml_parameter_covariance(fit)
      hessian <- stats::optimHess(fit$opt$par, function(x) evaluate(x)$value,
                                 function(x) evaluate(x)$gradient)
      expect_lt(max(abs(covariance$hessian - hessian)), 1e-7)
    }
    review <- mfrmr:::mfrmr_adaptive_quadrature_review(
      idx, config, mfrmr:::expand_params(fit$opt$par, sizes, config),
      mfrmr:::gauss_hermite_normal(order), fit$prep$levels$Person, order
    )
    person <- fit$facets$person
    expect_equal(nrow(person), nrow(review))
    index <- match(person$Person, review$Person)
    expect_false(anyNA(index))
    expect_lt(max(abs(person$Estimate - review$AdaptiveEAP[index])), 1e-9)
    expected_categories <- mfrmr:::compute_mml_expected_category_diagnostics(fit)
    expect_lt(abs(sum(expected_categories$posterior_bundle$person_bundle$log_marginal) -
                    fit$summary$LogLik), 1e-8)
    if (model == "GPCM") {
      boundary <- config$boundary_audit$gpcm_slope_boundary
      expect_identical(boundary$state, "not_evaluated_adaptive_quadrature")
      expect_false(boundary$fixed_quadrature_certificate)
      audit <- config$estimability_audit
      expect_identical(audit$mml_observed_pattern_score$status, "not_evaluated_adaptive_quadrature")
      expect_identical(audit$mml_all_pattern_information$status, "not_evaluated_adaptive_quadrature")
      expect_identical(audit$nonlinear_local_estimability$state, "not_evaluated")
      expect_identical(audit$nonlinear_local_estimability$reason_codes,
                       "adaptive_quadrature_probability_map_not_validated")
      expect_identical(audit$nonlinear_local_estimability$probability_model_scope,
                       "implemented_adaptive_quadrature_objective")
    }
  }
})

