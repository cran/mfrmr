anchor_population_replay_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    data <- load_mfrmr_data("example_core")
    data <- data[data$Person %in% unique(data$Person)[seq_len(20L)], ]
    anchors <- data.frame(Facet = "Rater", Level = as.character(data$Rater[1]), Anchor = .2)
    fit <- suppressWarnings(fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score",
      model = "RSM", method = "MML", noncenter_facet = "Rater", anchors = anchors,
      population_formula = ~1, person_data = unique(data["Person"]), population_sd = .8,
      quad_points = 21L, maxit = 200L, reltol = 1e-10, anchor_policy = "silent"))
    cached <<- list(data = data, fit = fit, anchors = anchors)
    cached
  }
})

test_that("same-model refits and executable replay retain fixed population SD", {
  x <- anchor_population_replay_fixture()
  args <- mfrmr_gqs_refit_arguments(x$fit, x$data, 31L)
  expect_identical(args$population_sd, .8)
  expect_identical(args$noncenter_facet, "Rater")
  expect_equal(as.data.frame(args$anchors), x$anchors)
  expect_identical(args$quad_points, 31L)

  lines <- build_replay_fit_mfrm_lines(x$fit$config$replay_inputs, x$fit$population,
    "Person", x$fit$config$source_columns, x$fit$config)
  env <- new.env(parent = environment())
  env$data <- x$data
  env$anchors <- x$anchors
  env$group_anchors <- NULL
  env$fit_person_data <- x$fit$population$person_table_replay
  suppressWarnings(eval(parse(text = lines), env))
  expect_identical(env$fit$population$fixed_sd, .8)
  expect_equal(env$fit$population$sigma2, .64)
  expect_identical(build_param_sizes(env$fit$config)$log_sigma2, 0L)
  expect_equal(env$fit$opt$par, x$fit$opt$par, tolerance = 1e-8)
  expect_equal(env$fit$config$facet_specs, x$fit$config$facet_specs)

  free <- x$fit
  free$population$fixed_sd <- free$config$population_spec$fixed_sd <- NULL
  expect_false("population_sd" %in% names(mfrmr_gqs_population_arguments(free)))
})

test_that("new-person scoring and saved reuse retain the fixed population prior", {
  x <- anchor_population_replay_fixture()
  data <- x$data[x$data$Person == x$data$Person[1], ]
  data$Person <- "NEW"
  prepared <- prepare_mfrm_prediction_data(x$fit, data)
  pop <- prepare_mfrm_prediction_population(x$fit, prepared)
  expect_true(pop$auto_person_data)
  expect_identical(pop$scaffold$fixed_sd, .8)
  expect_identical(pop$spec$fixed_sd, .8)
  expect_equal(pop$spec$sigma2, .64)
  score <- suppressWarnings(predict_mfrm_units(x$fit, data, readiness_policy = "review",
    scoring_quad_points = 41L, n_draws = 0L))
  expect_identical(score$settings$population_fixed_sd, .8)
  expect_equal(score$estimates$PriorSD, .8)
  expect_equal(score$estimates$PriorMean, unname(x$fit$population$coefficients))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x$fit, path)
  restored <- suppressWarnings(predict_mfrm_units(readRDS(path), data,
    readiness_policy = "review", scoring_quad_points = 41L, n_draws = 0L))
  expect_identical(restored$estimates, score$estimates)
  expect_identical(restored$fit_identity, score$fit_identity)
})

test_that("export metadata states the fixed population constraint", {
  x <- anchor_population_replay_fixture()
  diagnostics <- suppressWarnings(diagnose_mfrm(x$fit, residual_pca = "none"))
  manifest <- build_mfrm_manifest(x$fit, diagnostics = diagnostics)
  expect_identical(manifest$model_settings$Value[
    manifest$model_settings$Setting == "population_fixed_sd"], "0.8")
  expect_identical(manifest$settings$Value[
    manifest$settings$Setting == "fit_population_fixed_sd"], "0.8")
  replay <- build_mfrm_replay_script(x$fit, diagnostics = diagnostics)
  expect_match(replay$script, "population_sd = 0.8", fixed = TRUE)
  expect_identical(replay$settings$Value[
    replay$settings$Setting == "fit_population_fixed_sd"], "0.8")
})

test_that("unsupported forecasting and imputation do not change fixed-SD estimation", {
  x <- anchor_population_replay_fixture()
  expect_error(predict_mfrm_population(x$fit, reps = 1L),
    "do not yet preserve fixed `population_sd`", fixed = TRUE)
  data <- x$data[seq_len(8L), ]
  data$Event <- paste0("E", seq_len(nrow(data)))
  completed <- data
  data$Score[1] <- NA
  imputed <- mfrm_response_imputations(data, list(completed, completed),
    "Person", c("Rater", "Criterion"), "Score", "Event", "E1", 1:4,
    imputation_model = list(method = "preservation test only"))
  expect_error(fit_mfrm_imputed(imputed, population_sd = .8), "supported shared fitting")
  expect_error(build_conquest_overlap_bundle(x$fit),
    "ConQuest overlap template estimates residual variance", fixed = TRUE)
})
