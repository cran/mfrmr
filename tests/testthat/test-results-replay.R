test_that("ordinary result exports reopen without source-session objects", {
  fits <- readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits
  results <- lapply(fits, mfrm_results, include = "fit", compute = "never")
  fail <- function(...) stop("Replay must reload, not reconstruct or refit")
  local_mocked_bindings(mfrm_results = fail, fit_mfrm = fail,
    diagnose_mfrm = fail, compute_mml_parameter_covariance = fail,
    .package = "mfrmr")
  root <- tempfile("results-replay-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  modes <- list(default = "default", replay = "replay", both = c("rds", "replay"))
  for (model in names(results)) for (mode in names(modes)) {
    exported <- export_mfrm_results(results[[model]],
      output_dir = file.path(root, model, mode), prefix = "saved analysis",
      include = modes[[mode]], acknowledge_sensitive = TRUE)
    files <- exported$written_files
    expect_true("results_rds" %in% files$Component)
    script <- files$Path[files$Component == "replay_code"]
    expect_length(script, 1L)
    # The namespace parent provides functions, but no fit or diagnostics variables.
    fresh <- new.env(parent = asNamespace("mfrmr"))
    expect_no_error(capture.output(withr::with_dir(dirname(script),
      sys.source(basename(script), envir = fresh))))
    expect_identical(fresh$res, results[[model]])
    expect_identical(fresh$report, mfrm_report(results[[model]]))
  }
})
