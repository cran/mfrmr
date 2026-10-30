one_family_scoring_results_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      x <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))
      x$new <- x$data[x$data$Person %in% c("P1", "P34"), ]
      x$scores <- predict_mfrm_units(x$fit, x$new, n_draws = 2, seed = 42)
      cached <<- x
    }
    cached
  }
})

test_that("one-family saved EAP, priors and checks reach reports without recalculation", {
  x <- one_family_scoring_results_fixture()
  local_mocked_bindings(fit_mfrm = function(...) stop("unexpected fit"),
    diagnose_mfrm = function(...) stop("unexpected diagnostics"),
    predict_mfrm_units = function(...) stop("unexpected scoring"),
    prediction_source_scoring_readiness = function(...) stop("unexpected source checks"),
    mfrmr_adaptive_person_kernel = function(...) stop("unexpected integration"), .package = "mfrmr")
  res <- mfrm_results(x$fit, scores = x$scores, include = "fit", compute = "never")
  expect_identical(res$scores, x$scores)
  expect_identical(res$components$scores, x$scores)
  expect_identical(res$tables$person_scores, x$scores$estimates)
  expect_identical(res$tables$scoring_priors, x$scores$settings$prior_comparison)
  expect_identical(res$tables$scoring_integration, x$scores$settings$score_integration_review)
  expect_identical(res$tables$scoring_population_review, x$scores$population_review)
  expect_equal(res$tables$scoring_settings$PriorMean, unname(x$fit$population$coefficients[1]))
  expect_equal(res$tables$scoring_settings$PriorSD, sqrt(x$fit$population$sigma2))
  expect_identical(res$tables$scoring_overview$Status, "conditional")
  expect_true(is.na(res$tables$scoring_overview$NotScored))
  expect_false(x$scores$settings$source_inference_ready)
  expect_match(paste(summary(res)$reproducible_code$Code, collapse = "\n"), "readRDS", fixed = TRUE)
  expect_false(grepl("fit_mfrm(", paste(summary(res)$reproducible_code$Code, collapse = "\n"), fixed = TRUE))
  for (style in c("qc", "apa", "technical", "validation", "reviewer")) {
    report <- mfrm_report(res, style = style)
    expect_identical(report$tables$person_scores, x$scores$estimates)
    expect_identical(report$tables$scoring_priors, x$scores$settings$prior_comparison)
    expect_identical(report$first_screen$Status[report$first_screen$Area == "New-Person scoring"], "conditional")
    expect_match(report$markdown, "Posterior SD", fixed = TRUE)
    expect_match(report$markdown, "diagnostic normal bands", fixed = TRUE)
    expect_match(report$markdown, "population transport", fixed = TRUE)
  }
  expect_identical(mfrm_results(x$fit, scores = x$scores, include = "fit", compute = "never",
    output = "tables"), res$tables)
})

test_that("one-family attachments reject different calibration, population and status", {
  x <- one_family_scoring_results_fixture()
  changes <- list(
    function(f) { f$opt$par[1] <- f$opt$par[1] + .01; f },
    function(f) { f$prep$data$score_k[1] <- (f$prep$data$score_k[1] + 1L) %% 3L; f },
    function(f) { f$prep$score_map$OriginalScore <- 1:3; f },
    function(f) { f$config$score_map$OriginalScore <- f$prep$score_map$OriginalScore <- 1:3; f },
    function(f) { f$config$slope_facet <- f$config$step_facet; f },
    function(f) { f$config$estimation_control$quad_points <- 63L; f },
    function(f) { f$population$coefficients[1] <- f$population$coefficients[1] + .1; f },
    function(f) { f$population$sigma2 <- f$population$sigma2 * 2; f },
    function(f) { f$readiness$fit$NumericalState <- "failed"; f })
  for (change in changes)
    expect_error(mfrm_gpcm_results_score_tables(change(x$fit), x$scores), "fit/scores mismatch")
  old <- x$scores; old$fit_identity <- NULL
  expect_s3_class(summary(old), "summary.mfrm_unit_prediction")
  expect_error(mfrm_gpcm_results_score_tables(x$fit, old), "lacks its fitted-analysis identity")
  expect_error(mfrm_gpcm_results_score_tables(x$fit, summary(x$scores)), "must be saved")
  portable <- readRDS(test_path("fixtures", "calibration-gpcm-separate.rds"))$calibration
  expect_error(mfrm_gpcm_results_score_tables(x$fit, portable), "must be saved")
})

test_that("one-family attachment cannot relabel saved intervals or numerical checks", {
  x <- one_family_scoring_results_fixture()
  changes <- list(
    function(s) { s$settings$interval_level <- .8; s },
    function(s) { s$estimates$IntervalLevel[1] <- .8; s },
    function(s) { s$settings$scoring_algorithm <- "quadrature_eap_v1"; s },
    function(s) { s$settings$method <- "JML"; s },
    function(s) { s$settings$scoring_quad_points <- 61L; s },
    function(s) { s$settings$adaptive_quad_points <- 121L; s },
    function(s) { s$settings$score_integration_review <- NULL; s },
    function(s) { s$settings$score_integration_review$EAPChange[1] <- 0; s },
    function(s) { s$settings$prior_comparison <- NULL; s },
    function(s) { s$estimates$PriorMean[1] <- 99; s })
  for (change in changes)
    expect_error(mfrm_gpcm_results_score_tables(x$fit, change(x$scores)), "records")
  res <- mfrm_results(x$fit, scores = x$scores, include = "fit", compute = "never")
  for (change in list(
      function(r) { r$tables$person_scores$Estimate[1] <- 99; r },
      function(r) { r$tables$scoring_priors <- NULL; r },
      function(r) { r$components$scores <- NULL; r },
      function(r) { r$status$Status[r$status$Section == "person_scores"] <- "ready"; r })) {
    bad <- change(res)
    expect_error(summary(bad), "inconsistent with the attached scores")
    expect_error(mfrm_report(bad), "inconsistent with the attached scores")
    expect_error(plot(bad, type = "tables", draw = FALSE), "inconsistent with the attached scores")
    path <- tempfile()
    expect_error(export_mfrm_results(bad, output_dir = path, include = "tables",
      acknowledge_sensitive = TRUE), "inconsistent with the attached scores")
    expect_false(dir.exists(path))
  }
})

test_that("selected priors, source/batch review and omissions stay explicit", {
  x <- one_family_scoring_results_fixture()
  prior <- list(mean = .4, sd = .9)
  selected <- predict_mfrm_units(x$fit, x$new, scoring_prior = prior, interval_level = .9)
  res <- mfrm_results(x$fit, scores = selected, include = "fit", compute = "never")
  expect_equal(res$tables$scoring_settings$PriorMean, prior$mean)
  expect_equal(res$tables$scoring_settings$PriorSD, prior$sd)
  expect_identical(res$tables$scoring_settings$PosteriorBasis, "user_supplied_normal")
  expect_identical(res$tables$person_scores, selected$estimates)
  expect_true(all(res$tables$person_scores$PriorSource == "user_supplied"))
  expect_match(mfrm_report(res)$markdown, "analyst assumption", fixed = TRUE)
  changed <- selected; changed$settings$scoring_prior$mean <- 1
  expect_error(mfrm_gpcm_results_score_tables(x$fit, changed), "prior records")
  batch <- predict_mfrm_units(x$fit, x$new, scoring_quad_points = 2L, readiness_policy = "review")
  source_fit <- x$fit; source_fit$readiness$fit$NumericalState <- "failed"
  source <- predict_mfrm_units(source_fit, x$new, readiness_policy = "review")
  for (pair in list(list(fit = x$fit, scores = batch), list(fit = source_fit, scores = source))) {
    res <- mfrm_results(pair$fit, scores = pair$scores, include = "fit", compute = "never")
    expect_identical(res$tables$scoring_overview$Status, "review")
    expect_identical(res$tables$person_scores$EstimateUse, pair$scores$estimates$EstimateUse)
    expect_identical(mfrm_report(res)$first_screen$Status[
      mfrm_report(res)$first_screen$Area == "New-Person scoring"], "review")
  }
  missing <- x$new; missing$Score[1] <- NA_real_
  expect_warning(omitted <- predict_mfrm_units(x$fit, missing), "Dropped 1 row")
  res <- mfrm_results(x$fit, scores = omitted, include = "fit", compute = "never")
  expect_identical(res$tables$scoring_row_review, omitted$row_review)
  expect_true(is.na(res$tables$scoring_overview$NotScored))
})

test_that("ordinary GPCM saved score exports replay in a fresh session", {
  skip_if_not_installed("callr")
  x <- one_family_scoring_results_fixture()
  res <- mfrm_results(x$fit, scores = x$scores, include = "fit", compute = "never")
  folder <- tempfile(); dir.create(folder); on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  file <- file.path(folder, "results.rds"); saveRDS(res, file)
  worker <- function(root, file, folder, locale) {
    Sys.setlocale("LC_COLLATE", locale)
    if (file.exists(file.path(root, "R", "api-results.R"))) pkgload::load_all(root, quiet = TRUE, compile = FALSE) else
      library("mfrmr", lib.loc = dirname(root), character.only = TRUE)
    testthat::local_mocked_bindings(fit_mfrm = function(...) stop("no refit"),
      diagnose_mfrm = function(...) stop("no diagnostics"),
      predict_mfrm_units = function(...) stop("no scoring"),
      prediction_source_scoring_readiness = function(...) stop("no source recalculation"),
      mfrmr_adaptive_person_kernel = function(...) stop("no integration"), .package = "mfrmr")
    res <- readRDS(file)
    exported <- mfrmr::export_mfrm_results(res, output_dir = folder,
      include = c("tables", "report", "replay"), acknowledge_sensitive = TRUE)
    csv <- exported$written_files$Path[exported$written_files$Component == "table_person_scores"]
    table <- read.csv(csv, colClasses = c(Person = "character"), na.strings = "", check.names = FALSE)
    code <- exported$written_files$Path[exported$written_files$Format == "R"]
    replay <- new.env()
    withr::with_dir(folder, sys.source(basename(code), envir = replay))
    list(result = replay$res, report = mfrmr::mfrm_report(res)$tables, csv = table)
  }
  environment(worker) <- baseenv()
  actual <- callr::r(worker, list(normalizePath(find.package("mfrmr")), file, folder, Sys.getlocale("LC_COLLATE")), libpath = .libPaths())
  expect_identical(actual$result, res)
  expect_identical(actual$report, mfrm_report(res)$tables)
  expect_identical(actual$csv$Person, res$scores$estimates$Person)
  fields <- c("Estimate", "SD", "Lower", "Upper", "IntervalLevel", "PriorMean", "PriorSD")
  expect_equal(unname(as.matrix(actual$csv[fields])), unname(as.matrix(res$scores$estimates[fields])), tolerance = 1e-13)
})
