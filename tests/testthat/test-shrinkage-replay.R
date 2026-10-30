test_that("post-fit shrinkage replay preserves adjustment order and replacement settings", {
  data <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(suppressMessages(fit_mfrm(
    data, "Person", c("Rater", "Criterion"), "Score", method = "JML",
    maxit = 25, attach_diagnostics = TRUE)))
  original_inputs <- fit$config$replay_inputs
  adjusted <- apply_empirical_bayes_shrinkage(fit, facet_prior_sd = .7, shrink_person = TRUE)
  expect_identical(adjusted$config$replay_inputs, original_inputs)
  expect_identical(adjusted$config$shrinkage_settings,
                   list(facet_prior_sd = .7, shrink_person = TRUE, applied_after_fit = TRUE))
  expect_true(any(adjusted$facets$person$ShrinkageFactor > 0, na.rm = TRUE))

  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv(data, csv, row.names = FALSE)
  diagnostics <- suppressWarnings(suppressMessages(diagnose_mfrm(adjusted, residual_pca = "none")))
  script <- build_mfrm_replay_script(adjusted, diagnostics = diagnostics, data_file = csv)$script
  expect_match(script, 'facet_shrinkage = "none"', fixed = TRUE)
  expect_match(script, "fit <- apply_empirical_bayes_shrinkage(", fixed = TRUE)
  replay <- new.env(parent = globalenv())
  suppressWarnings(suppressMessages(eval(parse(text = script), envir = replay)))
  expect_equal(replay$fit$shrinkage_report, adjusted$shrinkage_report, tolerance = 1e-10)
  for (table in c("others", "person")) {
    cols <- c("ShrunkEstimate", "ShrunkSE", "ShrinkageFactor")
    expect_equal(replay$fit$facets[[table]][cols], adjusted$facets[[table]][cols], tolerance = 1e-10)
  }
  expect_equal(replay$fit$summary$LogLik, fit$summary$LogLik, tolerance = 1e-10)

  replaced <- apply_empirical_bayes_shrinkage(adjusted, facet_prior_sd = .4, shrink_person = FALSE)
  fresh <- apply_empirical_bayes_shrinkage(fit, facet_prior_sd = .4, shrink_person = FALSE)
  expect_identical(replaced$facets, fresh$facets)
  expect_identical(replaced$shrinkage_report, fresh$shrinkage_report)
  expect_false(any(c("ShrunkEstimate", "ShrunkSE", "ShrinkageFactor") %in% names(replaced$facets$person)))
  expect_false("Person" %in% replaced$shrinkage_report$Facet)
  expect_identical(replaced$config$replay_inputs, original_inputs)
  expect_identical(replaced$config$shrinkage_settings$shrink_person, FALSE)
  replaced_script <- build_mfrm_replay_script(replaced, diagnostics = diagnostics)$script
  expect_match(replaced_script, "facet_prior_sd = 0.4", fixed = TRUE)
  expect_match(replaced_script, "shrink_person = FALSE", fixed = TRUE)
  empirical <- apply_empirical_bayes_shrinkage(replaced, facet_prior_sd = NULL)
  expect_null(empirical$config$shrinkage_settings$facet_prior_sd)
  expect_match(build_mfrm_replay_script(empirical, diagnostics = diagnostics)$script,
               "facet_prior_sd = NULL", fixed = TRUE)

  old <- adjusted
  old$config$shrinkage_settings <- NULL
  expect_error(build_mfrm_replay_script(old, diagnostics = diagnostics), "Reapply.*original prior and Person")
  run <- structure(list(fit = adjusted, diagnostics = diagnostics), class = "mfrm_facets_run")
  expect_identical(build_mfrm_replay_script(run)$summary$ScriptMode, "fit")
  expect_error(build_mfrm_replay_script(run, script_mode = "facets"), "script_mode = 'fit'", fixed = TRUE)
  expect_match(build_mfrm_replay_script(run, script_mode = "fit")$script,
               "fit <- apply_empirical_bayes_shrinkage(", fixed = TRUE)
})
