secondary_result_fixture <- local({
  cases <- readRDS(test_path("fixtures", "mml-comparison-intervals.rds"))
  make <- function(fit) {
    diagnostics <- suppressWarnings(diagnose_mfrm(fit, diagnostic_mode = "legacy", residual_pca = "none"))
    list(fit = fit, diagnostics = diagnostics,
      bias = estimate_bias(fit, diagnostics, "Rater", "Criterion", max_iter = 1),
      unexpected = unexpected_response_table(fit, diagnostics),
      displacement = displacement_table(fit, diagnostics))
  }
  list(a = make(cases$RSM$fit), b = make(cases[["RSM-regression"]]$fit))
})

test_that("bias consumers reject another fit's results before fallbacks or writes", {
  x <- secondary_result_fixture
  expect_identical(x$a$fit$facets$others[c("Facet", "Level")],
                   x$b$fit$facets$others[c("Facet", "Level")])
  local_mocked_bindings(diagnose_mfrm = function(...) stop("Unexpected recomputation"),
                       fit_mfrm = function(...) stop("Unexpected refitting"), .package = "mfrmr")
  consumers <- list(unexpected_after_bias_table, facet_quality_dashboard,
    reporting_checklist, build_apa_outputs, facets_output_contract_review,
    run_qc_pipeline, build_mfrm_manifest, build_mfrm_replay_script)
  for (consume in consumers) {
    expect_error(consume(x$a$fit, diagnostics = x$a$diagnostics, bias_results = x$b$bias),
                 "fit/bias_results mismatch")
    expect_no_error(consume(x$b$fit, diagnostics = x$b$diagnostics, bias_results = x$b$bias))
  }
  directory <- tempfile("secondary-result-identity-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  expect_error(export_mfrm_bundle(x$a$fit, diagnostics = x$a$diagnostics,
    bias_results = x$b$bias, output_dir = directory, include = "manifest",
    acknowledge_sensitive = TRUE), "fit/bias_results mismatch")
  expect_false(dir.exists(directory))
  expect_no_error(export_mfrm_bundle(x$b$fit, diagnostics = x$b$diagnostics,
    bias_results = x$b$bias, output_dir = directory, include = "manifest",
    acknowledge_sensitive = TRUE))
})

test_that("casebooks bind independently supplied tables to their fitted analysis", {
  x <- secondary_result_fixture
  for (component in c("unexpected", "displacement")) {
    args <- list(fit = x$a$fit, diagnostics = x$a$diagnostics)
    args[[component]] <- x$b[[component]]
    expect_error(do.call(build_misfit_casebook, args), paste0("fit/", component, " mismatch"))
    args[[component]]$table <- args[[component]]$table[FALSE, ]
    expect_error(do.call(build_misfit_casebook, args), paste0("fit/", component, " mismatch"))
    args[[component]] <- x$a[[component]]
    expect_no_error(do.call(build_misfit_casebook, args))
    args[[component]]$fit_identity <- NULL
    expect_error(do.call(build_misfit_casebook, args), "no model refit is needed")
  }
})

test_that("all collection members and its primary result must match", {
  x <- secondary_result_fixture
  collection <- estimate_all_bias(x$a$fit, x$a$diagnostics,
    pairs = list(c("Rater", "Criterion")), max_iter = 1)
  check <- function(bias) mfrm_validate_bias_identity(x$a$fit, bias, "test")
  expect_no_error(check(collection))
  expect_no_error(check(list(matching = x$a$bias, absent = NULL)))
  expect_error(check(list(matching = x$a$bias, other = x$b$bias)), "fit/bias_results")
  wrong <- collection
  wrong$by_pair[[1]] <- x$b$bias
  expect_error(check(wrong), "fit/bias_results")
  wrong <- collection
  wrong$primary <- x$b$bias
  expect_error(check(wrong), "fit/bias_results\\$primary mismatch")
  collection$by_pair <- list()
  collection$primary <- NULL
  expect_no_error(check(collection))
  expect_error(mfrm_validate_bias_identity(x$b$fit, collection, "test"), "fit/bias_results mismatch")
  collection$fit_identity <- NULL
  expect_error(check(collection), "estimate_all_bias")
  old <- x$a$bias
  old$fit_identity <- NULL
  expect_error(check(old), "no model refit is needed")
  local_mocked_bindings(estimate_bias = function(...) list(), .package = "mfrmr")
  empty <- estimate_all_bias(x$a$fit, x$a$diagnostics,
    pairs = list(c("Rater", "Criterion")), drop_empty = FALSE)
  expect_length(empty$by_pair, 1L)
  expect_no_error(check(empty))
  expect_error(mfrm_validate_bias_identity(x$b$fit, empty, "test"), "fit/bias_results mismatch")
})

test_that("screening identity preserves saved reuse and detects calibration or scoring changes", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(secondary_result_fixture$b, path)
  x <- readRDS(path)
  expect_no_error(mfrm_validate_bias_identity(x$fit, x$bias, "saved reuse"))
  expect_identical(x$bias$fit_identity, x$unexpected$fit_identity)
  expect_identical(x$bias$fit_identity, x$displacement$fit_identity)
  variants <- list(
    function(f) { f$facets$person$Estimate[1] <- f$facets$person$Estimate[1] + .1; f },
    function(f) { f$steps$Estimate[1] <- f$steps$Estimate[1] + .1; f },
    function(f) { f$prep$data$Weight[1] <- 2; f },
    function(f) { f$config$population_spec$design_matrix[1, 2] <- 99; f }
  )
  for (change in variants) {
    expect_error(mfrm_validate_bias_identity(change(x$fit), x$bias, "changed fit"),
                 "fit/bias_results mismatch")
  }
  x$fit$call <- quote(changed_call())
  x$fit$opt$optimizer_diagnostics <- list(note = "annotation")
  expect_no_error(mfrm_validate_bias_identity(x$fit, x$bias, "annotations"))
  adjusted <- unexpected_after_bias_table(x$fit, x$bias, x$diagnostics)
  expect_identical(adjusted$fit_identity, x$bias$fit_identity)
})
