test_that("adaptive scoring survives artifact persistence and uses the same person grids", {
  for (model in c("RSM", "PCM")) {
    fixture <- adaptive_api_fixture(model)
    fit <- fixture$fit
    quadrature_review <- adaptive_sensitivity_fixture(model)
    fit <- quadrature_review$fits$q41
    rows <- fixture$data[rev(which(fixture$data$Person %in% unique(fixture$data$Person)[c(2, 7)])), ]
    rows$Weight <- rep(c(0.25, 1.75, 3), length.out = nrow(rows))
    predicted <- predict_mfrm_units(fit, rows, weight = "Weight", n_draws = 5L, seed = 421L)
    expect_identical(predicted$settings$scoring_algorithm, "adaptive_quadrature_eap_v2")
    expect_identical(predicted$draws,
      predict_mfrm_units(fit, rows, weight = "Weight", n_draws = 5L, seed = 421L)$draws)
    plausible <- sample_mfrm_plausible_values(fit, rows, weight = "Weight", n_draws = 5L, seed = 421L)
    expect_identical(plausible$values, predicted$draws)
    expect_false(any(grepl("fixed quadrature-grid", plausible$notes, fixed = TRUE)))
    draft <- extract_mfrm_calibration(fit, quadrature_review = quadrature_review)
    artifact <- freeze_mfrm_calibration(validate_mfrm_calibration(draft))
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    save_mfrm_calibration(artifact, path)
    restored <- load_mfrm_calibration(path)
    expect_identical(restored, artifact)
    scored <- score_mfrm_calibration(restored, rows, weight = "Weight")
    expect_identical(scored$settings$scoring_algorithm, "adaptive_quadrature_eap_v2")
    index <- match(scored$estimates$Person, predicted$estimates$Person)
    for (column in c("Estimate", "SD", "Lower", "Upper")) {
      expect_lt(max(abs(scored$estimates[[column]] - predicted$estimates[[column]][index])), 1e-10)
    }
    altered <- restored
    altered$scoring_basis$scoring_algorithm <- "quadrature_eap_v2"
    expect_error(score_mfrm_calibration(altered, rows), "SEMANTIC|INTEGRITY|IDENTITY")
  }
})

test_that("adaptive refits, replay and integration identities preserve the selected mode", {
  fixture <- adaptive_api_fixture()
  fit <- fixture$fit
  sensitivity <- adaptive_sensitivity_fixture()
  expect_identical(sensitivity$settings$mml_integration, "adaptive")
  expect_true(all(vapply(sensitivity$fits, function(x) {
    identical(x$config$estimation_control$mml_integration, "adaptive")
  }, TRUE)))
  lines <- mfrmr:::build_replay_fit_mfrm_lines(fit$config$replay_inputs, fit$population,
    NULL, fit$config$source_columns, fit$config)
  env <- new.env(parent = environment())
  env$data <- fixture$data
  env$anchors <- env$group_anchors <- NULL
  eval(parse(text = lines), env)
  expect_identical(env$fit$config$estimation_control$mml_integration, "adaptive")
  expect_equal(env$fit$opt$par, fit$opt$par, tolerance = 1e-10)
  expect_identical(mfrmr:::resolve_dff_refit_controls(fit)$mml_integration, "adaptive")
  expect_error(mfrmr:::validate_conquest_overlap_fit(fit), "requires a fixed-grid")
  unknown <- fit
  unknown$config$estimation_control$mml_integration <- "future_algorithm"
  expect_error(predict_mfrm_units(unknown, fixture$data), "Unrecognized MML integration mode")
  fixed_config <- fit$config
  fixed_config$estimation_control$mml_integration <- "fixed"
  expect_false(identical(mfrmr:::mfrm_ic_integration_evaluation_id("MML", fixed_config),
                         fit$summary$IntegrationEvaluationId))
})

test_that("adaptive integration refuses unsupported engines and preserves fixed defaults", {
  data <- load_mfrmr_data("example_core")
  args <- list(data = data, person = "Person", facets = c("Rater", "Criterion"), score = "Score",
               quad_points = 7L, maxit = 30L)
  for (extra in list(list(method = "JML"), list(mml_engine = "em"),
                    list(mml_engine = "hybrid"), list(checkpoint = list(path = tempfile())))) {
    expect_error(do.call(fit_mfrm, c(args, extra, list(mml_integration = "adaptive"))),
                 "requires MML with the direct engine and no checkpoint")
  }
  default <- suppressWarnings(do.call(fit_mfrm, args))
  explicit <- suppressWarnings(do.call(fit_mfrm, c(args, list(mml_integration = "fixed"))))
  expect_identical(default$opt$par, explicit$opt$par)
  expect_identical(default$summary, explicit$summary)
})
