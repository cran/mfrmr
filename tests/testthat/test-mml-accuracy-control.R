accuracy_small_arguments <- function(mass = 1, nodes = 31L) {
  list(data = data.frame(Person = "P1", Rater = rep(c("A", "B"), each = 4),
    Score = c(1, 1, 1, 0, 1, 0, 0, 0), Weight = mass),
    person = "Person", facets = "Rater", score = "Score", weight = "Weight",
    repeated_ratings = "independent", anchor_policy = "silent",
    mml_integration = "auto", quad_points = nodes, reltol = 1e-12)
}

# Independent continuous integral and observed information for one joint
# response vector. No package probability, quadrature or derivative is used.
accuracy_binary_reference <- function(b, mass) {
  b <- unname(b)
  joint <- function(theta) mass * (
    3 * plogis(theta - b, log.p = TRUE) +
    plogis(theta - b, lower.tail = FALSE, log.p = TRUE) +
    plogis(theta + b, log.p = TRUE) +
    3 * plogis(theta + b, lower.tail = FALSE, log.p = TRUE)) +
    dnorm(theta, log = TRUE)
  offset <- joint(0)
  integral <- function(k) integrate(function(theta) {
    a <- plogis(theta - b); z <- plogis(theta + b)
    score <- mass * (-2 + 4*a - 4*z)
    multiplier <- switch(k, rep(1, length(theta)), theta, theta^2,
      score, score^2 - 4*mass*(a*(1-a) + z*(1-z)))
    exp(joint(theta) - offset) * multiplier
  }, -Inf, Inf, rel.tol = 1e-10, abs.tol = 1e-12)$value
  values <- vapply(1:5, integral, 0)
  c(nll = -offset - log(values[1]), eap = values[2]/values[1],
    sd = sqrt(values[3]/values[1] - (values[2]/values[1])^2),
    gradient = -values[4]/values[1],
    information = (values[4]/values[1])^2 - values[5]/values[1])
}

accuracy_operational_fixture <- local({
  result <- NULL
  function() {
    if (is.null(result)) {
      data <- load_mfrmr_data("example_operational")
      fit <- fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score",
        mml_integration = "auto", quad_points = 7L, reltol = 1e-12,
        mml_accuracy_control = list(adaptive_points = c(3L, 31L)))
      result <<- list(data = data, fit = fit)
    }
    result
  }
})

test_that("adequate fixed integration is retained and the default is unchanged", {
  args <- accuracy_small_arguments(nodes = 81L)
  auto <- suppressWarnings(do.call(fit_mfrm, args))
  args$mml_integration <- "fixed"
  fixed <- suppressWarnings(do.call(fit_mfrm, args))
  args$mml_integration <- NULL
  default <- suppressWarnings(do.call(fit_mfrm, args))
  expect_identical(default$opt$par, fixed$opt$par)
  expect_identical(auto$opt$par, fixed$opt$par)
  expect_identical(auto$summary$MMLIntegration, "fixed")
  expect_identical(fixed$summary$MMLAccuracyStatus, "not_requested")
  expect_identical(auto$summary$MMLAccuracyStatus, "criteria_met")
  expect_equal(nrow(auto$config$mml_accuracy$history), 1L)
  expect_identical(auto$config$estimation_control$quad_points_requested, 81L)
})

test_that("narrow posterior repair agrees with an independent continuous integral", {
  args <- accuracy_small_arguments(mass = 50)
  fit <- suppressWarnings(do.call(fit_mfrm, args))
  history <- fit$config$mml_accuracy$history
  expect_identical(history$Integration, c("fixed", "adaptive"))
  expect_identical(history$Status, c("refinement_required", "criteria_met"))
  expect_gt(history$MaxPosteriorSDDifference[1], 0.01)
  b <- fit$facets$others$Estimate[fit$facets$others$Level == "A"]
  expect_length(b, 1L)
  reference <- accuracy_binary_reference(b, 50)
  expect_equal(fit$opt$value, unname(reference["nll"]), tolerance = 1e-9)
  expect_equal(fit$facets$person$Estimate, unname(reference["eap"]), tolerance = 1e-9)
  expect_equal(fit$facets$person$SD, unname(reference["sd"]), tolerance = 1e-9)
  expect_lt(abs(reference["gradient"]), 1e-4)
  covariance <- mfrmr:::compute_mml_parameter_covariance(fit)
  expect_equal(as.numeric(covariance$hessian), unname(reference["information"]),
    tolerance = 1e-5)
  expect_equal(as.numeric(covariance$cov), 1/unname(reference["information"]),
    tolerance = 1e-7)
})

test_that("refinement retains a coherent final object and preserves requested replay", {
  fixture <- accuracy_operational_fixture()
  fit <- fixture$fit
  history <- fit$config$mml_accuracy$history
  expect_identical(history$QuadPoints, c(7L, 3L, 31L))
  expect_identical(history$Status,
    c("refinement_required", "refinement_required", "criteria_met"))
  expect_identical(fit$config$estimation_control$mml_integration, "adaptive")
  expect_identical(fit$config$estimation_control$quad_points, 31L)
  expect_match(fit$summary$IntegrationEvaluationId, "aghq_normal_v1:q=31", fixed = TRUE)
  expect_match(paste(capture.output(summary(fit)), collapse = " "),
    "Automatic integration checks: criteria met", fixed = TRUE)
  expect_equal(nrow(fit$config$mml_accuracy$per_person), fit$prep$n_person)
  diagnostics <- mfrmr:::compute_mml_expected_category_diagnostics(fit)
  expect_equal(sum(diagnostics$posterior_bundle$person_bundle$log_marginal),
    fit$summary$LogLik, tolerance = 1e-9)
  restored <- unserialize(serialize(fit, NULL))
  scored <- predict_mfrm_units(restored, fixture$data, readiness_policy = "review")
  position <- match(fit$facets$person$Person, scored$estimates$Person)
  expect_false(anyNA(position))
  expect_equal(unname(scored$estimates$Estimate[position]), fit$facets$person$Estimate, tolerance = 1e-8)
  expect_equal(unname(scored$estimates$SD[position]), fit$facets$person$SD, tolerance = 1e-8)
  lines <- mfrmr:::build_replay_fit_mfrm_lines(fit$config$replay_inputs, fit$population,
    NULL, fit$config$source_columns, fit$config)
  env <- new.env(parent = environment())
  env$data <- fixture$data
  env$anchors <- env$group_anchors <- NULL
  eval(parse(text = lines), env)
  expect_identical(env$fit$config$mml_accuracy$control, fit$config$mml_accuracy$control)
  expect_identical(env$fit$config$mml_accuracy$history, history)
  expect_equal(env$fit$opt$par, fit$opt$par, tolerance = 1e-10)
  refit <- mfrmr:::mfrmr_gqs_refit_arguments(fit, fixture$data, 41L)
  expect_identical(refit$mml_integration, "adaptive")
  expect_null(refit$mml_accuracy_control)
  expect_identical(refit$quad_points, 41L)
})

test_that("auto integration preserves PCM, GPCM and latent regression specifications", {
  data <- load_mfrmr_data("example_operational")
  for (kind in c("PCM", "GPCM", "regression")) {
    args <- list(data = data, person = "Person", facets = c("Rater", "Criterion"),
      score = "Score", quad_points = 7L, reltol = 1e-12, mml_integration = "auto")
    if (kind != "regression") { args$model <- kind; args$step_facet <- "Criterion" }
    if (kind == "GPCM") args$slope_facet <- "Criterion"
    if (kind == "regression") {
      args$population_formula <- ~ X
      args$person_data <- data.frame(Person = unique(data$Person),
        X = rep(c(-1, 1), length.out = length(unique(data$Person))))
    }
    fit <- do.call(fit_mfrm, args)
    expect_identical(fit$config$mml_accuracy$status, "criteria_met")
    expect_identical(fit$config$estimation_control$mml_integration, "adaptive")
    sizes <- mfrmr:::build_param_sizes(fit$config)
    idx <- mfrmr:::build_indices(fit$prep, fit$config$step_facet, fit$config$slope_facet,
      fit$config$interaction_specs)
    reference <- mfrmr:::mfrmr_make_adaptive_mml_evaluator(idx, fit$config, sizes, 63L)(fit$opt$par)
    expect_equal(fit$summary$LogLik, -reference$value, tolerance = 1e-7)
    expect_lt(max(abs(reference$gradient)), 1e-4)
    if (kind %in% c("regression", "GPCM")) {
      expect_true(fit$population$active)
      expect_identical(fit$config$estimability_audit$nonlinear_local_estimability$state,
        "not_evaluated")
      expect_false(isTRUE(fit$summary$InferenceReady))
    }
  }
})

test_that("exhausted precision and optimizer failures never return a qualified fit", {
  data <- load_mfrmr_data("example_operational")
  error <- tryCatch(fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score",
    mml_integration = "auto", quad_points = 7L,
    mml_accuracy_control = list(adaptive_points = 3L)), error = identity)
  expect_s3_class(error, "mfrmr_mml_accuracy_error")
  expect_false(inherits(error, "mfrm_fit"))
  expect_identical(error$history$Status, rep("refinement_required", 2L))
  expect_identical(error$control$adaptive_points, 3L)
  expect_true(all(is.finite(error$last_parameters)))
  error <- tryCatch(fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score",
    mml_integration = "auto", maxit = 1L,
    mml_accuracy_control = list(adaptive_points = 3L)), error = identity)
  expect_s3_class(error, "mfrmr_mml_accuracy_error")
  expect_true(any(error$history$OptimizerConvergence != 0L))
  testthat::local_mocked_bindings(run_mfrm_optimization = function(...) {
    stop("Numerical evaluation unavailable.")
  }, .package = "mfrmr")
  error <- suppressWarnings(tryCatch(do.call(fit_mfrm, accuracy_small_arguments()),
    error = identity))
  expect_s3_class(error, "mfrmr_mml_accuracy_error")
  expect_identical(error$history$Status, "unavailable")
  expect_match(error$history$Detail, "Numerical evaluation unavailable", fixed = TRUE)
})

test_that("automatic integration rejects unsupported routes and malformed controls", {
  args <- accuracy_small_arguments()
  for (setting in list(list(method = "JML"), list(method = "CML"),
                      list(mml_engine = "em"), list(checkpoint = "checkpoint.rds"),
                      list(fit_scope = "comparisons"), list(slope_facet = c("A", "B")))) {
    expect_error(do.call(fit_mfrm, utils::modifyList(args, setting)), "ordinary full MML")
  }
  for (control in list(list(31), list(unknown = 1), list(loglik_tol = NA_real_),
      list(posterior_tol = 0), list(gradient_tol = Inf), list(adaptive_points = NULL),
      list(adaptive_points = c(31, 3)), list(adaptive_points = c(3, 3)),
      list(adaptive_points = 3.1), list(adaptive_points = 1))) {
    expect_error(do.call(fit_mfrm, c(args, list(mml_accuracy_control = control))),
      "mml_accuracy_control|adaptive_points|loglik_tol|posterior_tol|gradient_tol")
  }
  args$mml_integration <- "fixed"
  expect_error(do.call(fit_mfrm, c(args, list(mml_accuracy_control = list()))), "requires")
})
