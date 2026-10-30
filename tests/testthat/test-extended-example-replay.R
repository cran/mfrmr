test_that("ordinary saved inputs retain their data and replay matching diagnostics", {
  blocked <- function(...) stop("Unexpected statistical recomputation")
  local_mocked_bindings(fit_mfrm = blocked, diagnose_mfrm = blocked,
    .package = "mfrmr")
  for (name in c("core-jml", "bias-jml", "operational-mml")) {
    x <- readRDS(system.file("examples", paste0(name, ".rds"),
      package = "mfrmr", mustWork = TRUE))
    before <- unserialize(serialize(x, NULL))
    expected_dataset <- if (name == "core-jml") "example_core" else if (
      name == "bias-jml") "example_bias" else "example_operational"
    expect_identical(x$dataset, expected_dataset)
    ratings <- load_mfrmr_data(x$dataset)
    expect_equal(nrow(x$fit$prep$data), nrow(ratings))
    for (column in c("Person", "Rater", "Criterion", "Score"))
      expect_identical(as.character(x$fit$prep$data[[column]]),
        as.character(ratings[[column]]))
    expect_equal(x$fit$prep$n_person, length(unique(ratings$Person)))
    expect_identical(x$fit$config$model, "RSM")
    method <- if (grepl("jml", name)) "JML" else "MML"
    expect_identical(x$fit$config$method, method)
    expect_equal(x$fit$config$estimation_control$maxit,
      if (method == "JML") 300 else 400)
    expect_equal(x$fit$config$estimation_control$quad_points, 31)
    expect_identical(x$diagnostics$diagnostic_mode, "both")
    expect_identical(x$diagnostics$residual_pca_mode, "none")
    res <- mfrm_results(x$fit, diagnostics = x$diagnostics, compute = "never")
    expect_s3_class(res, "mfrm_results")
    expect_equal(plot(res, type = "wright", draw = FALSE)$data,
      plot_data(res, type = "wright"))
    expect_identical(x, before)
  }
  recipe <- new.env(parent = baseenv())
  sys.source(system.file("examples", "ordinary-fits.R", package = "mfrmr",
    mustWork = TRUE), envir = recipe)
  expect_true(is.function(recipe$ordinary_fit_examples))
})

test_that("ordinary saved analyses agree with fresh estimation and diagnostics", {
  skip_on_cran()
  recipe <- new.env(parent = baseenv())
  sys.source(system.file("examples", "ordinary-fits.R", package = "mfrmr",
    mustWork = TRUE), envir = recipe)
  fresh <- recipe$ordinary_fit_examples()
  expect_identical(names(fresh), c("core-jml", "bias-jml", "operational-mml"))
  controls <- c("maxit", "reltol", "quad_points", "mml_integration",
    "optimizer_requested", "mml_engine_requested")
  # Fair-score/displacement tables rank estimates. Floating-point ties can
  # change their row order across R versions; compare every value by identity.
  by_key <- function(table, keys) {
    expect_identical(anyDuplicated(table[keys]), 0L)
    table <- table[do.call(order, unname(table[keys])), , drop = FALSE]
    rownames(table) <- NULL
    table
  }
  align_tables <- function(diagnostics) {
    diagnostics$fair_average$raw_by_facet <- lapply(
      diagnostics$fair_average$raw_by_facet, by_key, keys = "Level")
    diagnostics$fair_average$by_facet <- lapply(
      diagnostics$fair_average$by_facet, by_key, keys = "Element")
    diagnostics$fair_average$stacked <- by_key(
      diagnostics$fair_average$stacked, c("Facet", "Element"))
    diagnostics$displacement$table <- by_key(
      diagnostics$displacement$table, c("Facet", "Level"))
    # Optimizer evaluation counts can differ across R/compiler builds even
    # when the fitted parameters, likelihood and convergence decisions agree.
    calls <- diagnostics$parameter_uncertainty$source$optimizer$FunctionEvaluations
    if (!is.null(calls)) {
      expect_true(is.numeric(calls) && all(is.finite(calls)) &&
                    all(calls >= 1 & calls == floor(calls)))
      diagnostics$parameter_uncertainty$source$optimizer$FunctionEvaluations <- NULL
    }
    diagnostics
  }
  max_numeric_error <- function(actual, expected) {
    if (is.numeric(actual) && is.numeric(expected)) {
      if (length(actual) != length(expected)) return(Inf)
      finite <- is.finite(actual) & is.finite(expected)
      return(max(c(0, abs(actual[finite] - expected[finite]))))
    }
    if (is.list(actual) && is.list(expected)) {
      if (length(actual) != length(expected)) return(Inf)
      return(max(c(0, vapply(seq_along(actual), function(i)
        max_numeric_error(actual[[i]], expected[[i]]), numeric(1)))))
    }
    0
  }
  for (name in names(fresh)) {
    saved <- readRDS(system.file("examples", paste0(name, ".rds"),
      package = "mfrmr", mustWork = TRUE))
    actual <- fresh[[name]]
    expect_identical(actual$fit$prep, saved$fit$prep)
    expect_identical(actual$fit$config$estimation_control[controls],
      saved$fit$config$estimation_control[controls])
    # Numerical tolerance permits floating-point variation, far below reporting
    # precision; data, controls and inference decisions must match exactly.
    expect_equal(actual$fit$opt$par, saved$fit$opt$par, tolerance = 1e-8)
    expect_equal(actual$fit$facets, saved$fit$facets, tolerance = 1e-8)
    expect_equal(actual$fit$summary$LogLik, saved$fit$summary$LogLik, tolerance = 1e-8)
    expect_identical(actual$fit$summary$InferenceReady, saved$fit$summary$InferenceReady)
    expect_identical(actual$fit$summary$ConvergenceSeverity,
      saved$fit$summary$ConvergenceSeverity)
    actual_diagnostics <- align_tables(actual$diagnostics)
    saved_diagnostics <- align_tables(saved$diagnostics)
    # Bias/displacement can be almost zero; a relative-only comparison divides
    # floating-point noise by that tiny value. Keep the intended 1e-8 absolute
    # precision, check every numeric entry, and still compare all other fields.
    comparison <- all.equal(actual_diagnostics, saved_diagnostics,
      tolerance = 1e-8, scale = 1)
    expect_true(isTRUE(comparison), info = paste(comparison, collapse = "\n"))
    expect_lte(max_numeric_error(actual_diagnostics, saved_diagnostics), 1e-8)
  }
})

test_that("saved measurement fits retain controls and support consumers without refitting", {
  x <- readRDS(system.file("examples", "measurement-fits.rds",
    package = "mfrmr", mustWork = TRUE))
  # Independent snapshot: R 4.5 may materialize deferred character names while
  # reading a fit. Byte serialization can change without any value/attribute change.
  before <- unserialize(serialize(x, NULL))
  expect_equal(vapply(x, function(f) f$config$estimation_control$quad_points,
    numeric(1)), c(gpcm = 61, conquest = 31, conquest_quick = 7))
  expect_equal(x$gpcm$config$slope_facet, "Rater")
  expect_equal(x$gpcm$config$step_facet, "Criterion")
  expect_identical(environment(x$gpcm$population$formula), baseenv())
  blocked <- function(...) stop("Unexpected refitting")
  local_mocked_bindings(fit_mfrm = blocked, .package = "mfrmr")
  recipe <- new.env(parent = baseenv())
  sys.source(system.file("examples", "measurement-fits.R", package = "mfrmr",
    mustWork = TRUE), envir = recipe)
  expect_true(is.function(recipe$measurement_fit_examples))
  ci <- confint(x$gpcm, parm = "slopes")
  expect_true(is.matrix(ci))
  expect_true(any(is.finite(ci)))
  for (name in c("conquest", "conquest_quick")) {
    b <- build_conquest_overlap_bundle(x[[name]], quad_points = 99)
    expect_equal(b$summary$MfrmrQuadraturePoints,
      x[[name]]$config$estimation_control$quad_points)
    expect_equal(b$summary$Persons, 60)
    expect_equal(b$summary$Items, 6)
  }
  legacy <- x$conquest
  legacy$config$quad_points <- legacy$config$estimation_control$quad_points
  legacy$config$estimation_control$quad_points <- NULL
  expect_equal(build_conquest_overlap_bundle(legacy)$summary$MfrmrQuadraturePoints, 31)
  legacy$config$quad_points <- NULL
  expect_error(build_conquest_overlap_bundle(legacy), "quadrature-point count is unavailable")
  expect_identical(x, before)
})

test_that("packaged extended examples are coherent and replay without recomputation", {
  path <- system.file("examples", "extended-models.rds", package = "mfrmr", mustWork = TRUE)
  x <- readRDS(path)
  before <- serialize(x, NULL)
  ratings <- load_mfrmr_data("example_core")
  for (model in x) {
    expect_true(model$fit$checks$NumericalReady)
    expect_true(model$fit$checks$InformationPositive)
    data <- model$fit$input$assigned_data
    if (is.null(data)) data <- model$fit$input$data
    expect_equal(data, ratings[names(data)], ignore_attr = TRUE)
    expect_equal(model$fit$settings$quad_points, 121)
    expect_true(all(is.na(model$fit$calibration_table$Lower)))
    expect_true(all(is.na(model$fit$calibration_table$Upper)))
  }
  expect_identical(x$random_rater$intervals$source, x$random_rater$fit)
  expect_true(all(is.na(x$random_rater$fit$raters$Lower)))
  expect_true(all(is.na(x$random_rater$fit$raters$Upper)))
  expect_equal(nrow(x$random_rater$intervals$trials), 19)
  expect_equal(nrow(x$random_rater$intervals$error), 19)
  expect_equal(x$testlet$scores$table$Person, sprintf("P%03d", 1:4))
  expect_equal(x$random_rater$scores$table$Person, sprintf("P%03d", 1:2))
  expect_equal(x$testlet$scores$scoring_data, x$testlet$fit$input$assigned_data)

  blocked <- function(...) stop("Unexpected statistical recomputation")
  local_mocked_bindings(fit_mfrm_testlet = blocked, fit_mfrm_random_rater = blocked,
    predict.mfrm_testlet = blocked, score_mfrm_random_rater = blocked,
    mfrm_random_rater_intervals = blocked, mfrm_response_diagnostics = blocked,
    .package = "mfrmr")
  testlet <- mfrm_results(x$testlet$fit, scores = x$testlet$scores,
    response_diagnostics = x$testlet$diagnostics, compute = "never")
  random <- mfrm_results(x$random_rater$fit, scores = x$random_rater$scores,
    intervals = x$random_rater$intervals, compute = "never")
  expect_s3_class(testlet, "mfrm_results")
  expect_s3_class(random, "mfrm_results")
  expect_equal(unname(summary(x$random_rater$intervals)$trials["Planned"]), 19)
  expect_equal(nrow(confint(x$random_rater$intervals, level = .90)), 4)
  for (result in list(testlet, random)) {
    p <- plot(result, type = "wright", draw = FALSE)
    expect_equal(p$data, plot_data(result, type = "wright"))
  }
  expect_equal(plot(x$testlet$diagnostics, draw = FALSE)$data,
    plot_data(x$testlet$diagnostics))
  expect_identical(serialize(x, NULL), before)

  # Reading the optional regeneration script must not start an estimator.
  recipe <- new.env(parent = baseenv())
  sys.source(system.file("examples", "extended-models.R", package = "mfrmr",
    mustWork = TRUE), envir = recipe)
  expect_true(is.function(recipe$extended_model_examples))
})

test_that("extended reports distinguish unsupplied response diagnostics from unavailable methods", {
  x <- readRDS(system.file("examples", "extended-models.rds", package = "mfrmr", mustWork = TRUE))
  blocked <- function(...) stop("Unexpected statistical recomputation")
  local_mocked_bindings(fit_mfrm_testlet = blocked, fit_mfrm_random_rater = blocked,
    diagnose_mfrm = blocked, mfrm_response_diagnostics = blocked, .package = "mfrmr")
  for (model in x) {
    res <- mfrm_results(model$fit, compute = "never")
    before <- serialize(res, NULL)
    report <- mfrm_report(res)
    row <- report$first_screen[report$first_screen$Area == "Response residuals", ]
    expect_equal(nrow(row), 1L)
    expect_identical(row$Status, "not_computed")
    expect_identical(row$Readiness, "Not supplied")
    expect_match(row$NextAction, "mfrm_response_diagnostics(fit)", fixed = TRUE)
    expect_match(row$MainIssue, "Ordinary residual diagnostics are not implemented", fixed = TRUE)
    expect_identical(res$status$Status[res$status$Section == "response_diagnostics"], "not_computed")
    expect_identical(res$status$Status[res$status$Section == "diagnostics"], "not_available")
    expect_equal(summary(report)$overview$NotComputedAreas, 1L)
    expect_identical(report$tables, res$tables)
    expect_identical(serialize(res, NULL), before)

    unresolved <- model$fit
    unresolved$checks$NumericalReady <- FALSE
    expect_identical(summary(mfrm_report(mfrm_results(unresolved, compute = "never")))$overview$OverallStatus,
      "review")
  }
  res <- mfrm_results(x$testlet$fit, response_diagnostics = x$testlet$diagnostics,
    compute = "never")
  report <- mfrm_report(res)
  row <- report$first_screen[report$first_screen$Area == "Response residuals", ]
  expect_identical(row$Status, "caveat")
  expect_identical(row$Readiness, "Descriptive only")
  expect_identical(report$tables$response_overview, res$tables$response_overview)
})
