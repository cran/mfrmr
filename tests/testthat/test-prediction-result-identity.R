prediction_identity_fixture <- local({
  cases <- readRDS(test_path("fixtures", "mml-comparison-intervals.rds"))
  fit <- cases$RSM$fit
  new_data <- data.frame(Person = c("NEW", "NEW"), Rater = fit$prep$levels$Rater[1],
    Criterion = fit$prep$levels$Criterion[1:2], Group = fit$prep$levels$Group[1], Score = c(1, 2))
  list(fit = fit, other = cases[["RSM-regression"]]$fit, new_data = new_data,
    diagnostics = diagnose_mfrm(fit, residual_pca = "none"),
    other_diagnostics = diagnose_mfrm(cases[["RSM-regression"]]$fit, residual_pca = "none"),
    unit_prediction = predict_mfrm_units(fit, new_data, readiness_policy = "review", n_draws = 2, seed = 1),
    plausible_values = sample_mfrm_plausible_values(fit, new_data, readiness_policy = "review", n_draws = 2, seed = 1))
})

test_that("fit-level prediction exports and replay reject another calibration", {
  x <- prediction_identity_fixture
  expect_identical(x$fit$facets$others[c("Facet", "Level")],
                   x$other$facets$others[c("Facet", "Level")])
  local_mocked_bindings(predict_mfrm_units = function(...) stop("Unexpected scoring"),
    prediction_source_scoring_readiness = function(...) stop("Unexpected numerical review"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"), .package = "mfrmr")
  for (arg in c("unit_prediction", "plausible_values")) {
    for (consume in list(build_mfrm_manifest, build_mfrm_replay_script)) {
      args <- c(list(fit = x$other, diagnostics = x$other_diagnostics), setNames(list(x[[arg]]), arg))
      expect_error(do.call(consume, args), paste0("fit/", arg, " mismatch"))
      args$fit <- x$fit
      args$diagnostics <- x$diagnostics
      expect_no_error(do.call(consume, args))
      args[[arg]]$fit_identity <- NULL
      expect_error(do.call(consume, args), "no model refit is needed")
      expect_s3_class(summary(args[[arg]]), paste0("summary.", class(x[[arg]])[1]))
    }
    directory <- tempfile("prediction-identity-")
    args <- c(list(fit = x$other, diagnostics = x$other_diagnostics,
      output_dir = directory, include = "predictions", acknowledge_sensitive = TRUE),
      setNames(list(x[[arg]]), arg))
    expect_error(do.call(export_mfrm_bundle, args), paste0("fit/", arg, " mismatch"))
    expect_false(dir.exists(directory))
    args$fit <- x$fit
    args$diagnostics <- x$diagnostics
    out <- do.call(export_mfrm_bundle, args)
    expect_true(all(file.exists(out$written_files$Path)))
    unlink(directory, recursive = TRUE)
  }
})

test_that("prediction identity covers future coding and scoring eligibility", {
  fit <- prediction_identity_fixture$other
  original <- prediction_fit_identity(fit)
  for (change in list(
    function(f) { f$steps$Estimate[1] <- f$steps$Estimate[1] + .1; f },
    function(f) { f$population$xlevels$Group <- rev(f$population$xlevels$Group); f },
    function(f) { f$population$contrasts$Group <- "contr.sum"; f },
    function(f) { attr(f$population$terms, "predvars") <- quote(list(relevel(Group, ref = "B"))); f },
    function(f) { f$readiness$fit$NumericalState <- "review"; f }
  )) {
    expect_false(identical(prediction_fit_identity(change(fit)), original))
  }
  # An unrelated formula workspace must not become part of saved-model identity.
  attr(fit$population$terms, ".Environment") <- new.env(parent = baseenv())
  environment(fit$population$terms)$unrelated <- 123
  expect_identical(prediction_fit_identity(fit), original)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  expect_identical(prediction_fit_identity(readRDS(path)), original)
  expect_identical(prediction_identity_fixture$unit_prediction$fit_identity,
                   prediction_identity_fixture$plausible_values$fit_identity)
})

test_that("fit-derived forecasts are bound while supplied scenarios remain separate", {
  fit <- readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits[[1]]
  pred <- suppressWarnings(predict_mfrm_population(fit, reps = 1, maxit = 10, seed = 812))
  expect_identical(pred$fit_identity, prediction_fit_identity(fit))
  check <- function(x, f = fit) export_validate_optional_object(x,
    "mfrm_population_prediction", "population_prediction", fit = f)
  expect_no_error(check(pred))
  expect_error(check(pred, prediction_identity_fixture$fit), "fit/population_prediction mismatch")
  old <- pred
  old$fit_identity <- NULL
  expect_error(check(old), "predict_mfrm_population")

  # Reuse the completed simulation to check the explicit-spec producer branch.
  local_mocked_bindings(evaluate_mfrm_design = function(...) pred$simulation, .package = "mfrmr")
  scenario <- predict_mfrm_population(sim_spec = pred$sim_spec, reps = 1, maxit = 10, seed = 812)
  expect_null(scenario$fit_identity)
  expect_identical(scenario$settings$source, "mfrm_sim_spec")
  expect_no_error(check(scenario, prediction_identity_fixture$fit))
  manifest <- build_mfrm_manifest(prediction_identity_fixture$fit,
    diagnostics = prediction_identity_fixture$diagnostics, population_prediction = scenario)
  expect_identical(manifest$settings$Value[manifest$settings$Setting == "population_prediction_source"], "mfrm_sim_spec")
  replay <- build_mfrm_replay_script(fit, population_prediction = pred)
  expect_match(replay$script, "not re-estimated from the replayed fit", fixed = TRUE)
  expect_match(replay$script, "sim_spec = population_prediction_sim_spec", fixed = TRUE)
})
