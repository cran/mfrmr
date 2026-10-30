portable_jml_fixture <- local({
  cache <- list()
  function(model) {
    if (is.null(cache[[model]])) {
      skip_if_not_installed("lpSolve")
      d <- load_mfrmr_data("example_core")
      d <- d[d$Person %in% unique(d$Person)[1:14], ]
      fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
        model = model, method = "JML", step_facet = if (model == "PCM") "Criterion" else NULL,
        quad_points = 1, maxit = 500, reltol = 1e-10))
      rows <- expand.grid(Person = c("001", "1", "NA", "sparse", "missing"),
        Rater = unique(d$Rater), Criterion = unique(d$Criterion), stringsAsFactors = FALSE)
      rows$Score <- rep(c(2, 3, 4), length.out = nrow(rows))
      rows$Score[rows$Person == "001"] <- min(d$Score)
      rows$Score[rows$Person == "1"] <- max(d$Score)
      rows$Score[rows$Person == "missing"] <- NA_real_
      sparse <- which(rows$Person == "sparse")
      rows <- rows[-sparse[-1], ]
      cache[[model]] <<- list(fit = fit, rows = rows)
    }
    cache[[model]]
  }
})

test_that("JML RSM and PCM artifacts replay reference-prior EAP without training Persons", {
  for (model in c("RSM", "PCM")) {
    f <- portable_jml_fixture(model)
    before <- serialize(f$fit, NULL)
    draft <- extract_mfrm_calibration(f$fit, scoring_quad_points = 141)
    expect_identical(draft$header$schema_version, 3L)
    expect_identical(draft$model$estimator, "JML")
    expect_identical(draft$scoring_basis$type, "post_hoc_standard_normal")
    expect_length(mfrmr:::mfrmr_calibration_find_prohibited(draft), 0)
    expect_identical(nrow(review_mfrm_calibration(draft)), 0L)
    frozen <- freeze_mfrm_calibration(validate_mfrm_calibration(draft))
    path <- tempfile(fileext = ".rds"); withr::defer(unlink(path))
    save_mfrm_calibration(frozen, path)
    loaded <- load_mfrm_calibration(path)
    expect_identical(loaded, frozen)
    for (prior in list(NULL, list(mean = .25, sd = 1.2))) {
      scored <- score_mfrm_calibration(loaded, f$rows, missing_response = "omit", scoring_prior = prior)
      native <- predict_mfrm_units(f$fit, f$rows[!is.na(f$rows$Score), ],
        scoring_quad_points = 141, scoring_prior = prior)
      fields <- c("Estimate", "SD", "Lower", "Upper")
      order <- match(scored$estimates$Person, native$estimates$Person)
      expect_equal(unname(as.matrix(scored$estimates[fields])),
        unname(as.matrix(native$estimates[order, fields])), tolerance = 1e-10)
      expect_true(all(scored$estimates$ScoreIntegrationReady))
      expect_setequal(scored$estimates$Person, c("001", "1", "NA", "sparse"))
      expect_identical(scored$person_dispositions$Disposition[
        scored$person_dispositions$Person == "missing"], "not_scored")
      expect_true(all(scored$person_dispositions$Disposition[
        scored$person_dispositions$Person %in% c("001", "1", "sparse")] == "scored_review"))
      expect_identical(scored$settings$source_scoring_evidence,
        frozen$eligibility$source_scoring_evidence)
      expect_s3_class(summary(scored), "summary.mfrm_calibration_score")
      expect_no_error(plot(scored, draw = FALSE, main = "New-cohort EAP scores"))
      expect_true(any(grepl("not estimated by JML", capture.output(print(scored)), fixed = TRUE)))
      saved <- tempfile(fileext = ".rds"); withr::defer(unlink(saved))
      saveRDS(scored, saved)
      expect_equal(summary(readRDS(saved)), summary(scored))
      broken <- scored; broken$settings$source_scoring_evidence$max_abs_gradient <- 1
      expect_error(summary(broken), "inconsistent")
    }
    expect_identical(serialize(f$fit, NULL), before)
    empty <- score_mfrm_calibration(frozen, f$rows[f$rows$Person == "missing", ], missing_response = "omit")
    expect_identical(nrow(empty$estimates), 0L)
    expect_s3_class(summary(empty), "summary.mfrm_calibration_score")
  }
})

test_that("JML portable source qualification cannot be supplied by flags alone", {
  f <- portable_jml_fixture("PCM")$fit
  original <- serialize(f, NULL)
  bad <- f; bad$opt$par[1] <- bad$opt$par[1] + 1
  expect_error(extract_mfrm_calibration(bad), "SOURCE_READINESS_INELIGIBLE")
  bad <- f; bad$opt$value <- bad$opt$value + 1
  expect_error(extract_mfrm_calibration(bad), "SOURCE_READINESS_INELIGIBLE")
  for (field in c("BoundaryState", "EstimabilityState", "NumericalState")) {
    bad <- f; bad$readiness$fit[[field]] <- "not_evaluated"
    expect_error(extract_mfrm_calibration(bad), "SOURCE_READINESS_INELIGIBLE")
  }
  bad <- f; bad$readiness$contract_version <- "old"
  expect_error(validate_mfrm_calibration(extract_mfrm_calibration(bad)), "SOURCE_READINESS")
  bad <- f; bad$prep$data$Weight[1] <- 2
  expect_error(extract_mfrm_calibration(bad), "SCORING_WEIGHT_UNSUPPORTED")
  bad <- f; bad$config$interaction_specs <- list(list())
  expect_error(extract_mfrm_calibration(bad), "MODEL_STRUCTURE_UNSUPPORTED")
  bad <- f; bad$config$method <- "unsupported"
  expect_error(extract_mfrm_calibration(bad), "MODEL_ESTIMATOR_UNSUPPORTED")
  expect_error(extract_mfrm_calibration(f, quadrature_review = list()), "does not use an MML")
  central <- portable_jml_fixture("PCM")$rows
  central <- central[central$Person == "NA", ]; central$Score <- 3
  low_grid <- freeze_mfrm_calibration(validate_mfrm_calibration(
    extract_mfrm_calibration(f, scoring_quad_points = 61)))
  expect_error(score_mfrm_calibration(low_grid, central), "SCORING_INTEGRATION_FAILED")
  draft <- extract_mfrm_calibration(f)
  for (change in c("evidence", "estimator", "prior", "schema")) {
    bad <- draft
    if (change == "evidence") bad$eligibility$source_scoring_evidence$max_abs_gradient <- 1
    if (change == "estimator") bad$model$estimator <- "MML"
    if (change == "prior") bad$scoring_basis$type <- "fixed_standard_normal"
    if (change == "schema") bad$header$schema_version <- 1L
    bad$integrity$semantic_components <- mfrmr:::mfrmr_calibration_semantic_components(bad)
    expect_gt(nrow(review_mfrm_calibration(bad)), 0)
    expect_error(validate_mfrm_calibration(bad), class = "mfrm_calibration_error")
  }
  expect_identical(serialize(f, NULL), original)
})

test_that("portable JML scoring runs in a new process with no source fit", {
  skip_if_not_installed("callr")
  root <- normalizePath(find.package("mfrmr"))
  for (model in c("RSM", "PCM")) {
    f <- portable_jml_fixture(model)
    frozen <- freeze_mfrm_calibration(validate_mfrm_calibration(
      extract_mfrm_calibration(f$fit, scoring_quad_points = 141)))
    path <- tempfile(fileext = ".rds"); withr::defer(unlink(path))
    save_mfrm_calibration(frozen, path)
    expected <- score_mfrm_calibration(frozen, f$rows, missing_response = "omit")
    worker <- function(root, path, rows) {
      if (file.exists(file.path(root, "R", "api-calibration.R"))) {
        pkgload::load_all(root, quiet = TRUE, compile = FALSE)
      } else {
        library("mfrmr", lib.loc = dirname(root), character.only = TRUE)
      }
      # Fitting and fitted-object scoring must not be used by this route.
      testthat::local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected fit"),
        predict_mfrm_units = function(...) stop("Unexpected native scoring"), .package = "mfrmr")
      scores <- mfrmr::score_mfrm_calibration(mfrmr::load_mfrm_calibration(path), rows, missing_response = "omit")
      list(estimates = scores$estimates, settings = scores$settings,
        review = summary(scores)$review)
    }
    environment(worker) <- baseenv()
    result <- callr::r(worker, args = list(root, path, f$rows), libpath = .libPaths())
    expect_identical(result$estimates, expected$estimates)
    expect_identical(result$settings, expected$settings)
    expect_identical(result$review, summary(expected)$review)
  }
})
