export_bundle_fixture <- local({
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

  dat <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    maxit = 20
  ))
  diagnostics <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "overall"))
  run <- suppressWarnings(run_mfrm_facets(
    dat,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 20
  ))
  bias_all <- suppressWarnings(estimate_bias(
    fit,
    diagnostics = diagnostics,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))

  list(
    fit = fit,
    diagnostics = diagnostics,
    run = run,
    bias_all = bias_all
  )
})

prediction_bundle_fixture <- local({
  dat <- load_mfrmr_data("example_core")
  keep_people <- unique(dat$Person)[1:18]
  dat <- dat[dat$Person %in% keep_people, , drop = FALSE]

  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    quad_points = 5,
    maxit = 15
  ))
  diagnostics <- diagnose_mfrm(fit, residual_pca = "none")

  spec <- build_mfrm_sim_spec(
    n_person = 20,
    n_rater = 4,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating"
  )

  population_prediction <- suppressWarnings(
    predict_mfrm_population(
      sim_spec = spec,
      reps = 2,
      maxit = 15,
      seed = 1
    )
  )

  new_units <- data.frame(
    Person = c("NEW01", "NEW01"),
    Rater = unique(dat$Rater)[1],
    Criterion = unique(dat$Criterion)[1:2],
    Score = c(2, 3)
  )

  unit_prediction <- predict_mfrm_units(
    fit,
    new_units,
    n_draws = 2,
    seed = 1
  )
  plausible_values <- sample_mfrm_plausible_values(
    fit,
    new_units,
    n_draws = 2,
    seed = 1
  )

  list(
    fit = fit,
    diagnostics = diagnostics,
    population_prediction = population_prediction,
    unit_prediction = unit_prediction,
    plausible_values = plausible_values
  )
})

test_that("build_mfrm_manifest captures reproducibility metadata", {
  manifest <- build_mfrm_manifest(
    fit = export_bundle_fixture$fit,
    diagnostics = export_bundle_fixture$diagnostics,
    bias_results = export_bundle_fixture$bias_all
  )

  expect_s3_class(manifest, "mfrm_manifest")
  expect_true(is.data.frame(manifest$summary))
  expect_true(is.data.frame(manifest$environment))
  expect_true(is.data.frame(manifest$available_outputs))
  expect_true(any(manifest$available_outputs$Component == "residual_pca"))
  expect_true(
    manifest$available_outputs$Available[manifest$available_outputs$Component == "bias_results"][1]
  )
  expect_equal(manifest$summary$Method[[1]], "JML")
  expect_equal(manifest$summary$MethodUsed[[1]], "JMLE")
  expect_equal(manifest$summary$Observations[[1]], nrow(export_bundle_fixture$fit$prep$data))
  expect_equal(manifest$summary$Persons[[1]], export_bundle_fixture$fit$config$n_person)
})

test_that("build_mfrm_manifest and replay script support FACETS-mode runs", {
  manifest <- build_mfrm_manifest(export_bundle_fixture$run)
  replay <- build_mfrm_replay_script(
    export_bundle_fixture$run,
    bias_results = export_bundle_fixture$bias_all,
    data_file = "analysis_data.csv"
  )

  expect_s3_class(manifest, "mfrm_manifest")
  expect_s3_class(replay, "mfrm_replay_script")
  expect_match(replay$script, "run_mfrm_facets\\(")
  expect_match(replay$script, "analysis_data\\.csv")
  expect_match(replay$script, "estimate_bias\\(")
})

test_that("build_mfrm_manifest records optional prediction artifacts", {
  manifest <- build_mfrm_manifest(
    fit = prediction_bundle_fixture$fit,
    diagnostics = prediction_bundle_fixture$diagnostics,
    population_prediction = prediction_bundle_fixture$population_prediction,
    unit_prediction = prediction_bundle_fixture$unit_prediction,
    plausible_values = prediction_bundle_fixture$plausible_values
  )

  expect_s3_class(manifest, "mfrm_manifest")
  expect_true(any(manifest$available_outputs$Component == "population_prediction"))
  expect_true(any(manifest$available_outputs$Component == "unit_prediction"))
  expect_true(any(manifest$available_outputs$Component == "plausible_values"))
  expect_true(
    manifest$available_outputs$Available[manifest$available_outputs$Component == "population_prediction"][1]
  )
  expect_true(
    manifest$available_outputs$Available[manifest$available_outputs$Component == "unit_prediction"][1]
  )
  expect_true(
    manifest$available_outputs$Available[manifest$available_outputs$Component == "plausible_values"][1]
  )
})

test_that("build_mfrm_replay_script reproduces optional prediction artifacts", {
  replay <- build_mfrm_replay_script(
    fit = prediction_bundle_fixture$fit,
    diagnostics = prediction_bundle_fixture$diagnostics,
    population_prediction = prediction_bundle_fixture$population_prediction,
    unit_prediction = prediction_bundle_fixture$unit_prediction,
    plausible_values = prediction_bundle_fixture$plausible_values,
    include_bundle = TRUE,
    bundle_prefix = "bundle_pred_test",
    data_file = "analysis_data.csv"
  )

  expect_s3_class(replay, "mfrm_replay_script")
  expect_match(replay$script, "predict_mfrm_population\\(")
  expect_match(replay$script, "predict_mfrm_units\\(")
  expect_match(replay$script, "sample_mfrm_plausible_values\\(")
  expect_match(
    replay$script,
    'include = c\\("core_tables", "checklist", "dashboard", "manifest", "html",\\s*"predictions"\\)'
  )
  expect_true(replay$summary$PopulationPrediction[[1]])
  expect_true(replay$summary$UnitPrediction[[1]])
  expect_true(replay$summary$PlausibleValues[[1]])
})

test_that("build_mfrm_replay_script preserves keep_original and rating range", {
  dat <- mfrmr:::sample_mfrm_data(seed = 42) |>
    dplyr::filter(.data$Score %in% c(1, 3, 5))

  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Task", "Criterion"),
    "Score",
    method = "JML",
    maxit = 25,
    keep_original = TRUE
  ))

  replay <- build_mfrm_replay_script(fit, data_file = "analysis_data.csv")

  expect_match(replay$script, "keep_original = TRUE", fixed = TRUE)
  expect_match(replay$script, "rating_min = 1", fixed = TRUE)
  expect_match(replay$script, "rating_max = 5", fixed = TRUE)
  expect_match(replay$script, "# Model: RSM | Method: JML | InternalMethod: JMLE", fixed = TRUE)
  expect_match(replay$script, 'method = "JML"', fixed = TRUE)
})

test_that("export_mfrm_bundle writes requested tables and html output", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-test")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expect_no_warning(
    bundle <- export_mfrm_bundle(
      fit = export_bundle_fixture$fit,
      diagnostics = export_bundle_fixture$diagnostics,
      bias_results = export_bundle_fixture$bias_all,
      output_dir = out_dir,
      prefix = "bundle_test",
      include = c("core_tables", "checklist", "dashboard", "apa", "anchors", "manifest", "visual_summaries", "script", "html"),
      overwrite = TRUE
    )
  )

  expect_s3_class(bundle, "mfrm_export_bundle")
  expect_true(is.data.frame(bundle$written_files))
  expect_true(any(bundle$written_files$Component == "bundle_html"))
  expect_true(any(grepl("bundle_test_manifest_summary.csv$", bundle$written_files$Path)))
  expect_true(any(grepl("bundle_test_checklist.csv$", bundle$written_files$Path)))
  expect_true(any(grepl("bundle_test_facet_dashboard_detail.csv$", bundle$written_files$Path)))
  expect_true(any(grepl("bundle_test_replay.R$", bundle$written_files$Path)))
  expect_true(any(grepl("bundle_test_visual_warning_counts.csv$", bundle$written_files$Path)))
  expect_true(file.exists(file.path(out_dir, "bundle_test_bundle.html")))
  expect_true(file.exists(file.path(out_dir, "bundle_test_manifest.txt")))
  expect_true(file.exists(file.path(out_dir, "bundle_test_replay.R")))
  expect_true(file.exists(file.path(out_dir, "bundle_test_visual_warning_map.txt")))
})

test_that("export_mfrm_bundle writes optional prediction artifacts", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-predictions")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  bundle <- export_mfrm_bundle(
    fit = prediction_bundle_fixture$fit,
    diagnostics = prediction_bundle_fixture$diagnostics,
    population_prediction = prediction_bundle_fixture$population_prediction,
    unit_prediction = prediction_bundle_fixture$unit_prediction,
    plausible_values = prediction_bundle_fixture$plausible_values,
    output_dir = out_dir,
    prefix = "bundle_pred_test",
    include = c("manifest", "predictions", "html"),
    overwrite = TRUE
  )

  expect_s3_class(bundle, "mfrm_export_bundle")
  expect_true(any(bundle$written_files$Component == "population_prediction_forecast"))
  expect_true(any(bundle$written_files$Component == "unit_prediction_estimates"))
  expect_true(any(bundle$written_files$Component == "plausible_values"))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_population_prediction_forecast.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_unit_prediction_estimates.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_plausible_values.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_population_prediction_ademp.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_population_prediction_sim_spec_settings.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_population_prediction_sim_spec_thresholds.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_unit_prediction_input.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_plausible_value_input.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_bundle.html")))

  html_lines <- readLines(file.path(out_dir, "bundle_pred_test_bundle.html"), warn = FALSE)
  html_text <- paste(html_lines, collapse = "\n")
  expect_match(html_text, "<h2>population_prediction_forecast</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>unit_prediction_estimates</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>plausible_value_summary</h2>", fixed = TRUE)

  unit_settings <- utils::read.csv(
    file.path(out_dir, "bundle_pred_test_unit_prediction_settings.csv"),
    stringsAsFactors = FALSE
  )
  expect_true(any(unit_settings$Setting == "source_columns.person"))
  expect_false(any(grepl("<list", unit_settings$Value, fixed = TRUE)))
})

test_that("export_mfrm_bundle requires explicit prediction objects for prediction export", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-predictions-missing")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expect_error(
    export_mfrm_bundle(
      fit = export_bundle_fixture$fit,
      diagnostics = export_bundle_fixture$diagnostics,
      output_dir = out_dir,
      prefix = "bundle_pred_missing",
      include = c("predictions"),
      overwrite = TRUE
    ),
    "`include = 'predictions'` requires at least one of `population_prediction`, `unit_prediction`, or `plausible_values`.",
    fixed = TRUE
  )
})

test_that("export_mfrm_bundle does not change the caller working directory", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-wd")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  original_wd <- getwd()
  bundle <- export_mfrm_bundle(
    fit = export_bundle_fixture$fit,
    diagnostics = export_bundle_fixture$diagnostics,
    output_dir = out_dir,
    prefix = "bundle_zip_test",
    include = c("manifest"),
    zip_bundle = TRUE,
    overwrite = TRUE
  )

  expect_identical(getwd(), original_wd)
  expect_s3_class(bundle, "mfrm_export_bundle")
  expect_true(any(bundle$written_files$Component == "bundle_zip"))
  expect_true(file.exists(file.path(out_dir, "bundle_zip_test_bundle.zip")))
})

test_that("export_mfrm_bundle respects overwrite for zip bundles", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-zip-overwrite")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expect_no_warning(
    export_mfrm_bundle(
      fit = export_bundle_fixture$fit,
      diagnostics = export_bundle_fixture$diagnostics,
      output_dir = out_dir,
      prefix = "bundle_zip_overwrite",
      include = c("manifest"),
      zip_bundle = TRUE,
      overwrite = TRUE
    )
  )

  expect_error(
    export_mfrm_bundle(
      fit = export_bundle_fixture$fit,
      diagnostics = export_bundle_fixture$diagnostics,
      output_dir = out_dir,
      prefix = "bundle_zip_overwrite",
      include = c("manifest"),
      zip_bundle = TRUE,
      overwrite = FALSE
    ),
    "File already exists:",
    fixed = TRUE
  )
})
