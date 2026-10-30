gmfrm_portable_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      x <- gmfrm_scoring_fixture()
      draft <- extract_mfrm_calibration(x$fit, scoring_quad_points = 121L,
        calibration_id = "gmfrm-portable", created_at_utc = "2026-10-01T00:00:00Z")
      validated <- validate_mfrm_calibration(draft, validated_at_utc = "2026-10-01T00:01:00Z")
      frozen <- freeze_mfrm_calibration(validated, frozen_at_utc = "2026-10-01T00:02:00Z")
      cached <<- c(x, list(draft = draft, validated = validated, frozen = frozen))
    }
    cached
  }
})

test_that("two-family portable scoring agrees with the fitted route and literal integrals", {
  x <- gmfrm_portable_fixture(); c <- x$frozen
  expect_identical(c$header$schema_version, 6L)
  expect_equal(nrow(review_mfrm_calibration(c)), 0L)
  expect_length(mfrmr_calibration_find_prohibited(c), 0L)
  expect_identical(c$model$slope_owner, c("Task", "Rater"))
  expect_identical(c$model$step_owner, "Rater")
  expect_identical(c$scoring_basis$type, "fixed_standard_normal")
  expect_identical(c$parameters$coordinates$Value[c$parameters$coordinates$ParameterClass == "slope"],
    as.numeric(x$fit$slopes$OptimizerEstimate))
  native <- predict_mfrm_units(x$fit, x$data, scoring_quad_points = 121L)
  portable <- score_mfrm_calibration(c, x$data)
  fields <- c("Estimate", "SD", "Lower", "Upper")
  expect_equal(unname(as.matrix(portable$estimates[fields])), unname(as.matrix(
    native$estimates[match(portable$estimates$Person, native$estimates$Person), fields])), tolerance = 1e-10)
  for (id in portable$estimates$Person) {
    reference <- gmfrm_scoring_reference(x$fit, x$data[x$data$Person == id, ])
    expect_equal(unname(unlist(portable$estimates[portable$estimates$Person == id, fields])),
      unname(reference), tolerance = 1e-7)
  }
  expect_true(all(portable$estimates$ScoreIntegrationReady))
  expect_identical(portable$settings$engine_identity, "artifact_coordinates_v6")
  expect_identical(portable$settings$semantic_components$product_structure$observed_contexts, c$model$observed_contexts)
  expect_s3_class(summary(portable), "summary.mfrm_calibration_score")
  expect_s3_class(plot(portable, draw = FALSE), "mfrm_plot_data")
  expect_match(plot_data(plot(portable, draw = FALSE))$caption, "Experimental two-family")
  expect_match(paste(capture.output(print(summary(portable))), collapse = " "), "Experimental two-family")
})

test_that("two-family artifacts reject structural corruption even after identity reconstruction", {
  x <- gmfrm_portable_fixture()
  changes <- list(
    function(c) { c$model$slope_owner <- rev(c$model$slope_owner); c },
    function(c) { c$model$step_owner <- c$model$slope_owner[1]; c },
    function(c) { c$model$slope_composition <- "single"; c },
    function(c) { c$model$facet_roles$Role[1] <- "facet_and_step_owner"; c },
    function(c) { c$scoring_basis$prior_sd <- 1.1; c },
    function(c) { c$constraints$identification$ConstraintType[2] <- "sum_to_zero"; c },
    function(c) { c$parameters$coordinates$Value[1] <- c$parameters$coordinates$Value[1] + .1; c },
    function(c) { at <- which(c$parameters$coordinates$ParameterClass == "owned_step")[1]; c$parameters$coordinates$Value[at] <- c$parameters$coordinates$Value[at] + .1; c },
    function(c) { at <- which(c$parameters$coordinates$ParameterClass == "slope")[1]; c$parameters$coordinates$Value[at] <- c$parameters$coordinates$Value[at] * 2; c },
    function(c) { c$parameters$coordinates <- head(c$parameters$coordinates, -1); c },
    function(c) { c$model$observed_contexts[[1]][1] <- "unknown"; c },
    function(c) { c$response$score_map$OriginalScore <- 1:3; c },
    function(c) { c$eligibility$source_scoring_evidence$local_calibration_review$integration$ComparisonMeanGradient <- 1e-4; c },
    function(c) { c$eligibility$source_scoring_evidence$inference_ready <- TRUE; c })
  for (change in changes) {
    bad <- change(x$frozen)
    bad$integrity$semantic_components <- mfrmr_calibration_semantic_components(bad)
    expect_gt(nrow(review_mfrm_calibration(bad)), 0L)
    expect_error(score_mfrm_calibration(bad, x$data), class = "mfrm_calibration_error")
  }
  bad <- x$frozen; bad$parameters$coordinates$Value[1] <- bad$parameters$coordinates$Value[1] + .1
  expect_true("IDENTITY_COMPONENT_MISMATCH" %in% review_mfrm_calibration(bad)$Code)
  bad <- x$frozen; bad$input_schema$source_columns$facets[] <- rev(bad$input_schema$source_columns$facets)
  expect_true("IDENTITY_COMPONENT_MISMATCH" %in% review_mfrm_calibration(bad)$Code)
  bad <- x$frozen; bad$parameters$person_estimates <- x$fit$facets$person
  expect_true(any(grepl("PROHIBITED", review_mfrm_calibration(bad)$Code)))
  bad <- x$frozen; bad$header$schema_version <- 2L
  expect_gt(nrow(review_mfrm_calibration(bad)), 0L)
  expect_error(score_mfrm_calibration(x$draft, x$data), "LIFECYCLE_NOT_FROZEN")
  expect_error(score_mfrm_calibration(x$validated, x$data), "LIFECYCLE_NOT_FROZEN")
  retired <- retire_mfrm_calibration(x$frozen, "retired-gmfrm")
  expect_error(score_mfrm_calibration(retired, x$data), "LIFECYCLE_NOT_FROZEN")
  expect_error(extract_mfrm_calibration(x$coarse), "SOURCE_READINESS_INELIGIBLE")
})

test_that("portable two-family input policies and integration refusals remain explicit", {
  x <- gmfrm_portable_fixture(); d <- x$data[x$data$Person == "NEW3", ]
  expect_error(score_mfrm_calibration(x$frozen, transform(d, Task = "unknown")), "SCORING_FACET_LEVEL_UNKNOWN")
  expect_error(score_mfrm_calibration(x$frozen, transform(d, Score = 3)), "SCORING_SCORE_UNKNOWN")
  expect_error(score_mfrm_calibration(x$frozen, transform(d, W = 2), weight = "W"), "SCORING_WEIGHT_UNSUPPORTED")
  expect_error(score_mfrm_calibration(x$frozen, rbind(d, d[1, ])), "SCORING_EVENT_DUPLICATE")
  expect_error(score_mfrm_calibration(x$frozen, d, event_id = "event"), "MODEL_STRUCTURE_UNSUPPORTED")
  expect_error(score_mfrm_calibration(x$frozen, d, scoring_prior = list(mean = 0, sd = 1)), "SCORING_BASIS_UNSUPPORTED")
  missing <- d[1, ]; missing$Person <- "missing"; missing$Score <- NA
  scored <- score_mfrm_calibration(x$frozen, rbind(d, missing), missing_response = "omit")
  expect_identical(scored$person_dispositions$Disposition[scored$person_dispositions$Person == "missing"], "not_scored")
  expect_false("missing" %in% scored$estimates$Person)
  empty <- score_mfrm_calibration(x$frozen, missing, missing_response = "omit")
  expect_equal(nrow(empty$estimates), 0L)
  expect_s3_class(summary(empty), "summary.mfrm_calibration_score")
  expect_error(score_mfrm_calibration(x$frozen, missing), "SCORING_SCORE_INVALID")
  low <- freeze_mfrm_calibration(validate_mfrm_calibration(extract_mfrm_calibration(x$fit, scoring_quad_points = 2L)))
  expect_error(score_mfrm_calibration(low, d), "SCORING_INTEGRATION_FAILED")
  stale <- scored; stale$settings$source_scoring_evidence <- NULL
  expect_error(summary(stale), "inconsistent prior or numerical-check records")
})

test_that("frozen two-family calibration scores in a fresh process without training state", {
  skip_if_not_installed("callr")
  x <- gmfrm_portable_fixture()
  path <- tempfile(); dir.create(path); on.exit(unlink(path, recursive = TRUE), add = TRUE)
  file <- file.path(path, "calibration.rds")
  save_mfrm_calibration(x$frozen, file)
  expect_identical(load_mfrm_calibration(file), x$frozen)
  expected <- score_mfrm_calibration(x$frozen, x$data)
  worker <- function(root, file, data) {
    if (file.exists(file.path(root, "R", "api-calibration.R"))) pkgload::load_all(root, quiet = TRUE, compile = FALSE) else
      library("mfrmr", lib.loc = dirname(root), character.only = TRUE)
    testthat::local_mocked_bindings(fit_mfrm = function(...) stop("no refitting"),
      predict_mfrm_units = function(...) stop("no native fit"),
      prediction_gmfrm_scoring_readiness = function(...) stop("no training data"),
      mfrm_gmfrm_problem = function(...) stop("no training data"),
      expand_params = function(...) stop("no optimizer parameters"), .package = "mfrmr")
    mfrmr::score_mfrm_calibration(mfrmr::load_mfrm_calibration(file), data)
  }
  environment(worker) <- baseenv()
  actual <- callr::r(worker, list(normalizePath(find.package("mfrmr")), file, x$data), libpath = .libPaths())
  expect_identical(actual, expected)
})
