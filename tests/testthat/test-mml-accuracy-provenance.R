accuracy_provenance_case <- function(fail_stage = Inf, check_error_stage = Inf,
                                     pass_stage = Inf, convergence = 0L) {
  stage <- 0L
  testthat::local_mocked_bindings(
    run_mfrm_optimization = function(start, ...) {
      stage <<- stage + 1L
      if (stage == fail_stage) stop("Controlled optimizer failure.")
      list(par = start + c(.25, -.5), value = 10 + stage, convergence = convergence)
    },
    mfrmr_check_mml_accuracy = function(...) {
      if (stage == check_error_stage) stop("Controlled accuracy-check failure.")
      list(passed = stage == pass_stage, per_person = data.frame(),
        metrics = c(MaxLogMarginalDifference = .1, MaxEAPDifference = .01,
          MaxPosteriorSDDifference = .01, MaxGradientDifference = .01, MaxGradient = .01))
    }, .package = "mfrmr")
  tryCatch(mfrmr:::mfrmr_optimize_accuracy_mml(start = c(0, 0), idx = NULL,
    config = list(model = "PCM", method = "MML", step_facet = "Item",
      estimation_control = list(quad_points = 7L, mml_integration = "auto")),
    sizes = list(facets = 2L), labels = "P1", maxit = 1L, reltol = 1e-8,
    optimizer = "auto", control = list(adaptive_points = c(3L, 5L))), error = identity)
}

test_that("failed later optimization cannot relabel retained parameters", {
  x <- accuracy_provenance_case(fail_stage = 2L)
  expect_s3_class(x, "mfrmr_mml_accuracy_error")
  expect_false(inherits(x, "mfrm_fit"))
  expect_identical(x$provenance_version, 2L)
  expect_identical(x$last_parameters, c(.25, -.5))
  expect_identical(x$last_integration, "fixed")
  expect_identical(x$last_quad_points, 7L)
  expect_identical(x$parameter_source, "optimizer_candidate")
  expect_identical(x$last_candidate$parameters, x$last_parameters)
  expect_identical(x$last_candidate$parameter_names, c("facets[1]", "facets[2]"))
  expect_identical(x$last_candidate$config$estimation_control$mml_integration, "fixed")
  expect_identical(x$last_candidate$diagnostics, x$history[1, , drop = FALSE])
  expect_identical(x$last_attempt$Stage, 2L)
  expect_identical(x$last_attempt$Integration, "adaptive")
  expect_identical(x$last_attempt$QuadPoints, 3L)
  expect_identical(x$last_attempt$FailurePhase, "optimization")
  expect_identical(x$last_attempt$Status, "unavailable")
  expect_match(conditionMessage(x), "Retained candidate: stage 1 (fixed, 7 points)", fixed = TRUE)
  expect_identical(unserialize(serialize(x, NULL)), x)
})

test_that("initial values are not presented as an optimizer candidate", {
  x <- accuracy_provenance_case(fail_stage = 1L)
  expect_s3_class(x, "mfrmr_mml_accuracy_error")
  expect_null(x$last_candidate)
  expect_identical(x$parameter_source, "initial_values")
  expect_identical(x$last_parameters, x$initial_parameters)
  expect_true(is.na(x$last_integration))
  expect_true(is.na(x$last_quad_points))
  expect_identical(x$last_attempt$Integration, "fixed")
  expect_identical(x$last_attempt$QuadPoints, 7L)
  expect_match(conditionMessage(x), "No optimizer candidate was returned", fixed = TRUE)
})

test_that("accuracy-check failure retains the current unqualified candidate", {
  x <- accuracy_provenance_case(fail_stage = 3L, check_error_stage = 2L, convergence = 1L)
  expect_identical(x$last_parameters, c(.5, -1))
  expect_identical(x$last_integration, "adaptive")
  expect_identical(x$last_quad_points, 3L)
  expect_identical(x$last_candidate$diagnostics$Stage, 2L)
  expect_identical(x$last_candidate$diagnostics$FailurePhase, "accuracy_check")
  expect_identical(x$last_candidate$diagnostics$Status, "unavailable")
  expect_identical(x$last_candidate$diagnostics$OptimizerConvergence, 1L)
  expect_identical(x$last_candidate$diagnostics$LogLik, -12)
  expect_identical(x$last_attempt$Stage, 3L)
  expect_identical(x$last_attempt$FailurePhase, "optimization")
  expect_match(x$last_candidate$diagnostics$Detail, "Controlled accuracy-check failure")
})

test_that("exhaustion and success keep their existing accuracy decisions", {
  x <- accuracy_provenance_case()
  expect_s3_class(x, "mfrmr_mml_accuracy_error")
  expect_identical(x$last_candidate$diagnostics$Status, "refinement_required")
  expect_identical(x$last_candidate$diagnostics, x$last_attempt)
  expect_identical(x$last_parameters, c(.75, -1.5))
  expect_identical(x$last_quad_points, 5L)
  ok <- accuracy_provenance_case(pass_stage = 2L)
  expect_false(inherits(ok, "error"))
  expect_identical(ok$opt$par, c(.5, -1))
  expect_identical(ok$config$estimation_control$quad_points, 3L)
  expect_identical(ok$config$mml_accuracy$status, "criteria_met")
})
