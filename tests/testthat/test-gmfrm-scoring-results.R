gmfrm_scoring_results_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      x <- gmfrm_scoring_fixture()
      ids <- c(NEW1 = "001", NEW2 = "1", NEW3 = "NA")
      x$data$Person <- unname(ids[x$data$Person])
      x$native <- predict_mfrm_units(x$fit, x$data, scoring_quad_points = 121L)
      x$calibration <- freeze_mfrm_calibration(validate_mfrm_calibration(
        extract_mfrm_calibration(x$fit, scoring_quad_points = 121L)))
      x$portable <- score_mfrm_calibration(x$calibration, x$data)
      cached <<- x
    }
    cached
  }
})

test_that("two-family score identities exclude the writing R version", {
  a <- structure(list(), names = character(), row.names = integer(), class = "data.frame")
  b <- structure(list(), row.names = integer(), names = character(), class = "data.frame")
  expect_identical(a, b)
  expect_identical(mfrm_checkpoint_fingerprint(a, canonical = TRUE),
    mfrm_checkpoint_fingerprint(b, canonical = TRUE))
  x <- readRDS(test_path("fixtures", "gmfrm-explicit-start.rds"))
  expect_identical(prediction_gmfrm_source_identity(x$fit),
    x$scores$settings$local_calibration_review$source_identity)
  expect_match(prediction_gmfrm_source_identity(x$fit), "^gmfrm_scoring_source_v2:")
  # The same fixture is read on every supported R/OS CI cell. Also retain
  # same-runtime v1 compatibility without pretending old hashes are portable.
  legacy <- x$scores
  legacy$settings$local_calibration_review$source_identity <-
    prediction_gmfrm_source_identity(x$fit, legacy = TRUE)
  expect_identical(mfrm_results(x$fit, include = "fit", compute = "never", scores = legacy)$scores, legacy)
})

test_that("two-family results retain native and portable scores without recalculation", {
  x <- gmfrm_scoring_results_fixture()
  local_mocked_bindings(fit_mfrm = function(...) stop("unexpected fit"),
    predict_mfrm_units = function(...) stop("unexpected scoring"),
    score_mfrm_calibration = function(...) stop("unexpected portable scoring"),
    prediction_gmfrm_scoring_readiness = function(...) stop("unexpected source evaluation"),
    mfrm_gmfrm_problem = function(...) stop("unexpected likelihood construction"),
    mfrmr_adaptive_person_kernel = function(...) stop("unexpected integration"), .package = "mfrmr")
  for (scores in list(x$native, x$portable)) {
    before <- scores
    res <- mfrm_results(x$fit, scores = scores)
    expect_identical(scores, before)
    expect_identical(res$scores, scores)
    expect_identical(res$components$scores, scores)
    expect_identical(res$tables$person_scores, scores$estimates)
    expect_identical(res$tables$scoring_integration, scores$settings$score_integration_review)
    expect_identical(res$tables$scoring_row_review, scores$row_review)
    expect_identical(mfrm_results(x$fit, scores = scores, compute = "never")$tables, res$tables)
    expect_identical(mfrm_results(x$fit, scores = scores, output = "tables"), res$tables)
    expect_identical(mfrm_results(x$fit, scores = scores, output = "summary"), summary(res))
    for (style in c("qc", "apa", "technical", "validation", "reviewer")) {
      report <- mfrm_report(res, style = style)
      expect_identical(report$tables$person_scores, scores$estimates)
      expect_identical(report$tables$scoring_integration, scores$settings$score_integration_review)
      expect_match(report$markdown, "calibration uncertainty")
      expect_match(report$markdown, "population transport")
      expect_true("New-Person scoring" %in% report$first_screen$Area)
    }
    html <- mfrm_results(x$fit, scores = scores, output = "html")
    expect_match(html$html, "scoring_integration")
    unlink(html$path)
  }
})

test_that("score attachments reject different sources and inconsistent saved tables", {
  x <- gmfrm_scoring_results_fixture()
  changes <- list(
    function(f) { f$opt$par[1] <- f$opt$par[1] + .01; f },
    function(f) { f$prep$data$score_k[1] <- (f$prep$data$score_k[1] + 1L) %% 3L; f },
    function(f) { f$prep$score_map$OriginalScore <- 1:3; f },
    function(f) { f$config$slope_facet <- rev(f$config$slope_facet); f },
    function(f) { f$config$estimation_control$quad_points <- 63L; f },
    function(f) { f$facets$others$Estimate[1] <- 10; f },
    function(f) { f$summary$Converged <- FALSE; f })
  for (scores in list(x$native, x$portable)) {
    for (change in changes)
      expect_error(mfrm_gpcm_results_score_tables(change(x$fit), scores), "Saved scores must match")
    old <- scores
    if (inherits(old, "mfrm_unit_prediction")) {
      old$settings$local_calibration_review$source_identity <- NULL
    } else {
      old$settings$source_scoring_evidence$local_calibration_review$source_identity <- NULL
      old$settings$semantic_components$source_scoring_review$source_scoring_evidence <- old$settings$source_scoring_evidence
    }
    expect_s3_class(summary(old), class(summary(scores))[1])
    expect_error(mfrm_results(x$fit, scores = old), "Re-score older outputs")
    res <- mfrm_results(x$fit, scores = scores)
    for (mutate in list(
        function(r) { r$tables$person_scores$Estimate[1] <- 99; r },
        function(r) { r$tables$scoring_integration <- NULL; r },
        function(r) { r$components$scores <- NULL; r },
        function(r) { r$status$Status[r$status$Section == "person_scores"] <- "ready"; r })) {
      bad <- mutate(res)
      expect_error(summary(bad), "inconsistent with the attached scores")
      expect_error(mfrm_report(bad), "inconsistent with the attached scores")
      expect_error(plot(bad, type = "tables", draw = FALSE), "inconsistent with the attached scores")
      path <- tempfile()
      expect_error(export_mfrm_results(bad, output_dir = path, include = "tables",
        acknowledge_sensitive = TRUE), "inconsistent with the attached scores")
      expect_false(dir.exists(path))
    }
  }
  expect_error(mfrm_results(x$fit, scores = summary(x$native)), "must be saved")
})

test_that("review-only and not-scored outcomes survive result collection", {
  x <- gmfrm_scoring_results_fixture()
  d <- x$data[x$data$Person == "NA", ]
  source_review <- predict_mfrm_units(x$coarse, d, scoring_quad_points = 121L, readiness_policy = "review")
  batch_review <- predict_mfrm_units(x$fit, d, scoring_quad_points = 2L, readiness_policy = "review")
  pairs <- list(list(fit = x$coarse, scores = source_review), list(fit = x$fit, scores = batch_review))
  for (pair in pairs) {
    res <- mfrm_results(pair$fit, scores = pair$scores)
    expect_identical(res$tables$scoring_overview$Status, "review")
    expect_true(all(startsWith(res$tables$person_scores$EstimateUse, "review_only")))
    report <- mfrm_report(res)
    expect_identical(report$first_screen$Status[report$first_screen$Area == "New-Person scoring"], "review")
    expect_identical(report$tables$person_scores$EstimateUse, pair$scores$estimates$EstimateUse)
  }
  d$Score <- NA_real_
  omitted <- score_mfrm_calibration(x$calibration, d, missing_response = "omit")
  res <- mfrm_results(x$fit, scores = omitted)
  expect_equal(res$tables$scoring_overview$Scored, 0L)
  expect_equal(res$tables$scoring_overview$NotScored, 1L)
  expect_identical(res$tables$scoring_overview$Status, "review")
  expect_identical(mfrm_report(res)$tables$scoring_person_dispositions, omitted$person_dispositions)
})

test_that("saved two-family interval bounds cannot be relabelled", {
  x <- gmfrm_scoring_results_fixture()
  changes <- list(
    function(s) { s$settings$interval_level <- .8; s },
    function(s) { s$estimates$IntervalLevel[1] <- .8; s },
    function(s) { s$settings$interval_level <- NA_real_; s },
    function(s) { s$estimates$IntervalLevel <- NULL; s },
    function(s) { s$settings$scoring_algorithm <- "quadrature_eap_v1"; s },
    function(s) { s$estimates$ScoringAlgorithm[1] <- "adaptive_quadrature_eap_v2"; s })
  for (scores in list(x$native, x$portable)) {
    for (change in changes) {
      bad <- change(scores)
      expect_error(summary(bad), "inconsistent interval records")
      expect_error(mfrm_results(x$fit, scores = bad), "inconsistent interval records")
    }
    display <- summary(scores, digits = 0)
    expect_identical(display$estimates$IntervalLevel, scores$estimates$IntervalLevel)
    display$settings$interval_level <- .8
    expect_error(print(display), "inconsistent interval records")
  }
  # A different requested level is valid when the bounds are actually recomputed.
  native <- predict_mfrm_units(x$fit, x$data, scoring_quad_points = 121L, interval_level = c(requested = .8765))
  portable <- score_mfrm_calibration(x$calibration, x$data, interval_level = c(requested = .8765))
  for (scores in list(native, portable)) {
    expect_true(all(summary(scores, digits = 0)$estimates$IntervalLevel == .8765))
    expect_equal(mfrm_results(x$fit, scores = scores)$tables$scoring_settings$IntervalLevel, .8765)
    expect_true(all(scores$estimates$Lower > x$native$estimates$Lower))
    expect_true(all(scores$estimates$Upper < x$native$estimates$Upper))
  }
  bad <- x$portable
  bad$settings$quadrature_order <- 31L
  expect_error(summary(bad), "inconsistent prior or numerical-check records")
})

test_that("source identity cannot substitute for the saved scoring calibration", {
  x <- gmfrm_scoring_results_fixture()
  changes <- list(
    function(s) { s$settings$two_family_calibration$slope_components[[2]]$Slope[1] <- 2; s },
    function(s) { s$settings$two_family_calibration$locations$Estimate[1] <- 2; s },
    function(s) { s$settings$two_family_calibration$steps$Estimate[1] <- 2; s })
  for (change in changes)
    expect_error(mfrm_results(x$fit, scores = change(x$native)), "calibration records do not match")
  changes <- list(
    function(s) { s$settings$semantic_components$parameter_coordinates$Value[1] <- 2; s },
    function(s) { s$settings$semantic_components$facet_order_roles_levels_signs$facet_signs[1] <- 1; s },
    function(s) { s$settings$semantic_components$slope_specification$slope_owner <- c("Rater", "Task"); s },
    function(s) { s$settings$semantic_components$response_map$score_map$OriginalScore <- 1:3; s },
    function(s) { s$settings$semantic_components$product_structure$observed_contexts <- NULL; s })
  for (change in changes)
    expect_error(mfrm_results(x$fit, scores = change(x$portable)), "calibration records do not match")
  # An internally valid artifact with a changed free slope is also a different source.
  changed <- x$calibration
  at <- which(changed$parameters$coordinates$ParameterClass == "slope" &
    changed$parameters$coordinates$OwnerFacet == "Rater")[1]
  changed$parameters$coordinates$Value[at] <- changed$parameters$coordinates$Value[at] * 1.01
  changed$integrity$semantic_components <- mfrmr_calibration_semantic_components(changed)
  expect_equal(nrow(review_mfrm_calibration(changed)), 0L)
  scored <- score_mfrm_calibration(changed, x$data)
  expect_error(mfrm_results(x$fit, scores = scored), "calibration records do not match")
})

test_that("saved scoring attachments replay and export in a fresh session", {
  skip_if_not_installed("callr")
  x <- gmfrm_scoring_results_fixture()
  results <- lapply(list(x$native, x$portable), function(s) mfrm_results(x$fit, scores = s))
  folder <- tempfile(); dir.create(folder); on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  file <- file.path(folder, "results.rds"); saveRDS(results, file)
  worker <- function(root, file, folder) {
    if (file.exists(file.path(root, "R", "api-results.R"))) pkgload::load_all(root, quiet = TRUE, compile = FALSE) else
      library("mfrmr", lib.loc = dirname(root), character.only = TRUE)
    testthat::local_mocked_bindings(fit_mfrm = function(...) stop("no refit"),
      predict_mfrm_units = function(...) stop("no native scoring"),
      score_mfrm_calibration = function(...) stop("no portable scoring"),
      prediction_gmfrm_scoring_readiness = function(...) stop("no source recalculation"),
      mfrm_gmfrm_problem = function(...) stop("no likelihood"),
      mfrmr_adaptive_person_kernel = function(...) stop("no integration"), .package = "mfrmr")
    results <- readRDS(file)
    lapply(seq_along(results), function(i) {
      res <- results[[i]]
      path <- file.path(folder, paste0("export", i))
      exported <- mfrmr::export_mfrm_results(res, output_dir = path,
        include = c("tables", "report", "replay"), acknowledge_sensitive = TRUE)
      csv <- exported$written_files$Path[exported$written_files$Component == "table_person_scores"]
      table <- read.csv(csv, colClasses = c(Person = "character"), na.strings = "", check.names = FALSE)
      code <- exported$written_files$Path[exported$written_files$Format == "R"]
      replay <- new.env()
      withr::with_dir(path, sys.source(basename(code), envir = replay))
      list(result = replay$res, summary = summary(res), report = mfrmr::mfrm_report(res)$tables, csv = table)
    })
  }
  environment(worker) <- baseenv()
  actual <- callr::r(worker, list(normalizePath(find.package("mfrmr")), file, folder), libpath = .libPaths())
  for (i in seq_along(results)) {
    expect_identical(actual[[i]]$result, results[[i]])
    expect_identical(actual[[i]]$summary, summary(results[[i]]))
    expect_identical(actual[[i]]$report, mfrm_report(results[[i]])$tables)
    expect_identical(actual[[i]]$csv$Person, results[[i]]$scores$estimates$Person)
    fields <- c("Estimate", "SD", "Lower", "Upper", "IntervalLevel", "PriorMean", "PriorSD")
    expect_equal(unname(as.matrix(actual[[i]]$csv[fields])),
      unname(as.matrix(results[[i]]$scores$estimates[fields])), tolerance = 1e-13)
  }
})
