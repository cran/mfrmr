# Deterministic execution failures, not estimator replications.
stage_path <- test_path("..", "..", "inst", "validation", "mfrm-wide-map-mml-stages-20261001.R")
skip_if_not(file.exists(stage_path), "Repository-only MML execution helper")
stage_env <- new.env(parent = globalenv())
sys.source(stage_path, envir = stage_env)
stage_mock <- function(...) {
  bindings <- list(...)
  old <- mget(names(bindings), envir = stage_env)
  list2env(bindings, envir = stage_env)
  withr::defer(list2env(old, envir = stage_env), envir = parent.frame())
}

stage_fixture <- function() {
  fit <- readRDS(test_path("fixtures", "native-location-fixed-gpcm.rds"))
  args <- list(data = data.frame(Person = "p", Rater = c("1", "2", "3"), Score = 0:2),
    person = "Person", facets = "Rater", score = "Score", model = "GPCM", method = "MML",
    slope_facet = "Rater", step_facet = "Rater", gpcm_mml_identification = "fixed_standard_normal",
    category_policy = "preserve", mml_engine = "direct", mml_integration = "fixed",
    optimizer = "BFGS", maxit = 400L, reltol = 1e-9)
  spec <- list(ConditionId = "execution-test", Arm = "fixed-GPCM", Replicate = 1L,
    InputId = "artificial-control-flow", args = args, facet = "Rater", level = .95,
    purpose = "workflow_witness_only")
  score <- function(ready = FALSE) list(ready = ready, parameter_ready = TRUE,
    reason_codes = if (ready) character() else "local_calibration_check_failed",
    local_calibration_review = list(basis = "local_mml_and_integration_v1", eligible = ready,
      review = "fixture", integration = data.frame(OriginalOrder = 31, ComparisonOrder = 61,
        NLLChange = if (ready) 0 else .01, GradientChange = 0, ComparisonGradient = 0)))
  ci <- function(fit, pass = TRUE) list(fit = fit,
    settings = list(procedure = "mml_native_location_model_v1"),
    checks = data.frame(Check = c("Joint information", "Quadrature comparison"), Passed = c(TRUE, pass)),
    numerical_checks = data.frame(InverseResidual = 0, ComparisonInverseResidual = 0,
      CurvatureScaledGradient = 0, QuadratureScoreShift = if (pass) 0 else .1,
      QuadratureCovarianceChange = 0),
    table = data.frame(Target = c("1", "2", "3"), Estimate = rep(fit$config$estimation_control$quad_points, 3),
      SE = if (pass) 1 else NA_real_, Lower = if (pass) 30 else NA_real_, Upper = if (pass) 32 else NA_real_,
      Status = if (pass) "available" else "numerical_review_failed", CIEligible = pass, InferenceReview = "fixture"))
  list(fit = fit, spec = spec, score = score, ci = ci)
}

test_that("refinement requires finite integration-only failures", {
  x <- stage_fixture(); s <- x$score(); ci <- x$ci(x$fit, FALSE)
  expect_true(stage_env$wide_mml_scoring_retry(s))
  expect_false(stage_env$wide_mml_scoring_retry(x$score(TRUE)))
  s$local_calibration_review$integration$ComparisonGradient <- NA_real_
  expect_false(stage_env$wide_mml_scoring_retry(s))
  s$local_calibration_review$integration <- NULL
  expect_false(stage_env$wide_mml_scoring_retry(s))
  expect_true(stage_env$wide_mml_interval_retry(ci))
  ci$numerical_checks$ComparisonInverseResidual <- .01
  expect_false(stage_env$wide_mml_interval_retry(ci))
  ci$numerical_checks$ComparisonInverseResidual <- 0
  ci$checks$Passed[1] <- FALSE
  expect_false(stage_env$wide_mml_interval_retry(ci))
  ci$checks$Passed[1] <- NA
  expect_false(stage_env$wide_mml_interval_retry(ci))
  expect_false(stage_env$wide_mml_interval_retry(NULL))
})

test_that("outputs keep their first qualifying stage and replay without computation", {
  x <- stage_fixture(); calls <- integer(); ci_calls <- integer()
  stage_mock(wide_mml_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) {
      calls <<- c(calls, args$quad_points); f <- x$fit
      f$config$estimation_control$quad_points <- args$quad_points; f
    }, wide_mml_score = function(fit) x$score(fit$config$estimation_control$quad_points == 61L),
    wide_mml_interval = function(fit, ...) {
      ci_calls <<- c(ci_calls, fit$config$estimation_control$quad_points); x$ci(fit)
    })
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  result <- stage_env$wide_mml_run(x$spec, out)
  expect_identical(calls, c(31L, 61L)); expect_identical(ci_calls, 31L)
  expect_match(result$point$stage, "^q61:"); expect_match(result$interval$stage, "^q31:")
  records <- stage_env$wide_mml_records(x$spec, result)
  expect_true(all(records$Available))
  expect_equal(records$Estimate[records$Output == "interval"], rep(31, 3))
  stage_env$wide_mml_fit <- stage_env$wide_mml_score <- stage_env$wide_mml_interval <- function(...) stop("Unexpected computation")
  expect_identical(stage_env$wide_mml_run(x$spec, out), result)
  bad <- x$spec; bad$args$data$Score[1] <- 1L
  expect_error(stage_env$wide_mml_run(bad, out), "data/settings/source changed")
  stage_env$wide_mml_runtime <- function() list(source = "changed")
  expect_error(stage_env$wide_mml_run(x$spec, out), "data/settings/source changed")
})

test_that("an interrupted review resumes the saved fit and completed review", {
  x <- stage_fixture(); calls <- 0L; score_calls <- 0L; interrupt <- TRUE
  stage_mock(wide_mml_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) { calls <<- calls + 1L; x$fit },
    wide_mml_score = function(fit) { score_calls <<- score_calls + 1L; x$score(TRUE) },
    wide_mml_interval = function(fit, ...) {
      if (interrupt) signalCondition(structure(list(message = "test interruption"), class = c("interrupt", "condition")))
      x$ci(fit)
    })
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  caught <- tryCatch(stage_env$wide_mml_run(x$spec, out), interrupt = function(e) "interrupted")
  expect_identical(caught, "interrupted")
  expect_identical(readRDS(file.path(out, "status.rds"))$status, "interrupted")
  expect_true(file.exists(file.path(out, "q31-fit.rds")))
  expect_false(file.exists(file.path(out, "q31-interval.rds")))
  interrupt <- FALSE
  result <- stage_env$wide_mml_run(x$spec, out)
  expect_identical(calls, 1L); expect_identical(score_calls, 1L)
  expect_identical(readRDS(file.path(out, "status.rds"))$status, "complete")
  # A corrupt saved output is refused even when selected.rds exists.
  path <- file.path(out, "q31-interval.rds"); old <- readRDS(path)
  old$payload$value$table$Estimate[1] <- 999; saveRDS(old, path)
  expect_error(stage_env$wide_mml_run(x$spec, out), "checksum mismatch")
})

test_that("refinement stops at q121 and retains earlier outputs on a fit error", {
  x <- stage_fixture(); calls <- integer(); fail <- FALSE
  stage_mock(wide_mml_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) {
      calls <<- c(calls, args$quad_points)
      if (fail && args$quad_points == 61L) stop("solver failure")
      f <- x$fit; f$config$estimation_control$quad_points <- args$quad_points; f
    }, wide_mml_score = function(fit) x$score(FALSE),
    wide_mml_interval = function(fit, ...) x$ci(fit))
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  result <- stage_env$wide_mml_run(x$spec, out)
  expect_identical(calls, c(31L, 61L, 121L))
  expect_false(result$point$value$ready); expect_match(result$point$stage, "^q121:")
  expect_match(result$interval$stage, "^q31:")
  failed_out <- tempfile(); on.exit(unlink(failed_out, recursive = TRUE), add = TRUE)
  fail <- TRUE; calls <- integer()
  result <- stage_env$wide_mml_run(x$spec, failed_out)
  expect_identical(calls, c(31L, 61L)); expect_identical(result$point$error, "solver failure")
  expect_match(result$interval$stage, "^q31:")
  records <- stage_env$wide_mml_records(x$spec, result)
  expect_true(all(records$Status[records$Output == "point"] == "error"))
  expect_true(all(records$Available[records$Output == "interval"]))
})

test_that("a later interval retains its own point without replacing the scoring source", {
  x <- stage_fixture(); calls <- integer(); score_calls <- 0L
  stage_mock(wide_mml_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) {
      calls <<- c(calls, args$quad_points)
      f <- x$fit; f$config$estimation_control$quad_points <- args$quad_points; f
    }, wide_mml_score = function(fit) { score_calls <<- score_calls + 1L; x$score(TRUE) },
    wide_mml_interval = function(fit, ...) x$ci(fit, fit$config$estimation_control$quad_points == 61L))
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  result <- stage_env$wide_mml_run(x$spec, out)
  expect_identical(calls, c(31L, 61L)); expect_identical(score_calls, 1L)
  expect_match(result$point$stage, "^q31:"); expect_match(result$interval$stage, "^q61:")
  records <- stage_env$wide_mml_records(x$spec, result)
  expect_equal(records$Estimate[records$Output == "interval"], rep(61, 3))
  expect_equal(records$Estimate[records$Output == "point"],
    unname(x$fit$facets$others$Estimate[x$fit$facets$others$Facet == "Rater"]))
})

test_that("a local-source refusal is settled without an integration retry", {
  x <- stage_fixture(); calls <- integer()
  stage_mock(wide_mml_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) { calls <<- c(calls, args$quad_points); x$fit },
    wide_mml_score = function(fit) {
      z <- x$score(FALSE); z$local_calibration_review$integration <- NULL; z
    }, wide_mml_interval = function(...) stop("unsupported source"))
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  result <- stage_env$wide_mml_run(x$spec, out)
  expect_identical(calls, 31L)
  records <- stage_env$wide_mml_records(x$spec, result)
  expect_false(any(records$Available))
  expect_true(all(records$Status[records$Output == "point"] == "returned"))
  expect_true(all(is.finite(records$Estimate[records$Output == "point"])))
  expect_true(all(records$Status[records$Output == "interval"] == "error"))
})

stage_product_spec <- function(n = 40L) {
  x <- stage_fixture()$spec
  x$args$data <- expand.grid(Person = paste0("p", seq_len(n)), Task = 1:3, Rater = 1:3)
  x$args$data$Score <- rep(0:2, length.out = nrow(x$args$data))
  x$args$facets <- x$args$slope_facet <- c("Task", "Rater")
  x$args$noncenter_facet <- x$args$step_facet <- "Rater"
  x$args$mml_integration <- "adaptive"; x$args$maxit <- 500L; x$args$reltol <- 1e-10
  x$args$gpcm_mml_start <- "neutral_em"; x$args$rating_min <- 0L; x$args$rating_max <- 2L
  x
}

test_that("two-family refinement retains its per-Person and interval thresholds", {
  s <- list(ready = FALSE, parameter_ready = TRUE,
    local_calibration_review = list(basis = "two_family_point_calibration_v1",
      integration = data.frame(OriginalOrder = 31, ComparisonOrder = 61,
        NLLChangePerPerson = 2e-6, ComparisonMeanGradient = 0)))
  expect_true(stage_env$wide_mml_scoring_retry(s))
  s$local_calibration_review$integration$NLLChangePerPerson <- 0
  s$local_calibration_review$integration$ComparisonMeanGradient <- 2e-6
  expect_true(stage_env$wide_mml_scoring_retry(s))
  s$local_calibration_review$integration$ComparisonMeanGradient <- Inf
  expect_false(stage_env$wide_mml_scoring_retry(s))
  s$local_calibration_review$integration <- NULL
  expect_false(stage_env$wide_mml_scoring_retry(s))
  x <- stage_fixture(); ci <- x$ci(x$fit, FALSE)
  ci$settings <- list(two_family = TRUE)
  ci$checks$Check[2] <- "Quadrature sensitivity"
  ci$numerical_checks$MaximumMeanScore <- 0
  expect_true(stage_env$wide_mml_interval_retry(ci))
  ci$checks$Check[2] <- "Local score rank"
  expect_false(stage_env$wide_mml_interval_retry(ci))
  ci$checks$Check[2] <- "Quadrature sensitivity"
  ci$numerical_checks$MaximumMeanScore <- 1e-4
  expect_false(stage_env$wide_mml_interval_retry(ci))
  expect_true(stage_env$wide_mml_validate(stage_product_spec()))
  bad <- stage_product_spec(); bad$args$noncenter_facet <- "Task"
  expect_error(stage_env$wide_mml_validate(bad))
  bad <- stage_product_spec(); bad$args$mml_integration <- "fixed"
  expect_error(stage_env$wide_mml_validate(bad))
})

test_that("two-family location execution uses its own consumer", {
  x <- stage_fixture(); spec <- stage_product_spec()
  stage_mock(wide_mml_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) x$fit,
    wide_mml_score = function(fit) x$score(TRUE),
    wide_mml_interval = function(fit, ...) {
      ci <- x$ci(fit); ci$settings <- list(two_family = TRUE); ci
    })
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  result <- stage_env$wide_mml_run(spec, out)
  expect_true(result$interval$value$settings$two_family)
  expect_true(all(stage_env$wide_mml_records(spec, result)$Available))
  bad <- spec; bad$InputId <- "another-input"
  expect_error(stage_env$wide_mml_records(bad, result), "exact job specification")
})

test_that("a proved score-rank exclusion does not exclude the point or launch intervals", {
  x <- stage_fixture(); spec <- stage_product_spec(10L)
  stage_mock(wide_mml_runtime = function() list(source = "test"),
    wide_mml_fit = function(args) x$fit,
    wide_mml_score = function(fit) x$score(TRUE),
    wide_mml_interval = function(...) stop("Excluded interval was attempted"))
  out <- tempfile(); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  result <- stage_env$wide_mml_run(spec, out)
  expect_null(result$interval)
  expect_equal(result$interval_plan$parameters, 13)
  expect_equal(result$interval_plan$persons, 10)
  plan <- stage_env$wide_mml_output_plan(spec)
  expect_true(all(plan$Eligibility[plan$Output == "point"] == "eligible"))
  expect_true(all(plan$Eligibility[plan$Output == "interval"] == "design_excluded"))
  records <- stage_env$wide_mml_records(spec, result)
  expect_identical(unique(records$Output), "point")
  expect_true(all(records$Available))
  # N >= p is necessary only; it must still run the actual numerical rank check.
  expect_identical(stage_env$wide_mml_interval_plan(stage_product_spec(13L))$eligibility, "eligible")
})
