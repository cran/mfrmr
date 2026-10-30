test_that("two-family new-Person moments and intervals match independent continuous integrals", {
  x <- gmfrm_scoring_fixture(); before <- x$fit
  pred <- predict_mfrm_units(x$fit, x$data, scoring_quad_points = 121L, n_draws = 5, seed = 42)
  expect_identical(x$fit, before)
  expect_identical(pred$settings$source_scoring_status, "conditional")
  expect_false(pred$settings$source_inference_ready)
  expect_true(all(pred$estimates$ScoreIntegrationReady))
  expect_true(all(pred$draws$ScoreIntegrationReady))
  expect_true(all(pred$estimates$PriorMean == 0 & pred$estimates$PriorSD == 1))
  columns <- c("Estimate", "SD", "Lower", "Upper")
  for (id in pred$estimates$Person) {
    reference <- gmfrm_scoring_reference(x$fit, x$data[x$data$Person == id, ])
    expect_equal(unname(unlist(pred$estimates[pred$estimates$Person == id, columns])),
      unname(reference), tolerance = 1e-7)
  }
  expect_identical(pred$settings$two_family_calibration$slope_facets, c("Task", "Rater"))
  expect_identical(pred$settings$two_family_calibration$score_map, x$fit$prep$score_map)
  expect_identical(pred$settings$two_family_calibration$locations, x$fit$facets$others)
  expect_match(paste(pred$notes, collapse = " "), "calibration.*excluded")
  expect_match(paste(pred$notes, collapse = " "), "sampling coverage")
  # Saved summaries/exports retain the scores and checks without integrating again.
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
  saveRDS(pred, file); restored <- readRDS(file)
  local_mocked_bindings(mfrmr_adaptive_person_kernel = function(...) stop("no integration"),
    mfrm_gmfrm_em = function(...) stop("no refitting"), .package = "mfrmr")
  expect_identical(summary(restored), summary(pred))
  expect_identical(restored$estimates, pred$estimates)
  expect_silent(export_validate_optional_object(restored, "mfrm_unit_prediction", "prediction", fit = x$fit))
  expect_s3_class(build_summary_table_bundle(restored), "mfrm_summary_table_bundle")
  missing <- restored; missing$settings$score_integration_review <- NULL
  expect_error(summary(missing), "inconsistent numerical review")
  missing <- restored; missing$settings$two_family_calibration <- NULL
  expect_error(summary(missing), "inconsistent calibration records")
  stale <- restored; stale$settings$local_calibration_review$integration$ComparisonMeanGradient <- 1
  expect_error(summary(stale), "inconsistent calibration records")
})

test_that("two-family source and scoring integration failures stay separate", {
  x <- gmfrm_scoring_fixture(); d <- x$data[x$data$Person == "NEW3", ]
  expect_error(predict_mfrm_units(x$coarse, d), "not ready for fitted-object scoring")
  review <- predict_mfrm_units(x$coarse, d, scoring_quad_points = 121L, readiness_policy = "review")
  expect_false(review$settings$source_scoring_ready)
  expect_true(all(review$estimates$ScoreIntegrationReady))
  expect_true(all(review$estimates$EstimateUse == "review_only_nonready_source"))
  expect_s3_class(summary(review), "summary.mfrm_unit_prediction")
  expect_error(predict_mfrm_units(x$fit, d, scoring_quad_points = 2L), "Posterior scoring integration did not pass")
  coarse <- predict_mfrm_units(x$fit, d, scoring_quad_points = 2L, readiness_policy = "review")
  expect_true(coarse$settings$source_scoring_ready)
  expect_false(all(coarse$estimates$ScoreIntegrationReady))
  expect_true(all(coarse$estimates$EstimateUse == "review_only_scoring_integration"))
  for (change in list(
      function(f) { f$opt$par[1] <- f$opt$par[1] + .1; f },
      function(f) { f$slopes$Estimate[1] <- NA_real_; f },
      function(f) { f$config$slope_facet <- rev(f$config$slope_facet); f },
      function(f) { f$prep$score_map$OriginalScore <- 1:3; f },
      function(f) { f$population$active <- TRUE; f },
      function(f) { f$opt$convergence <- 1L; f })) {
    for (policy in c("error", "review"))
      expect_error(predict_mfrm_units(change(x$fit), d, readiness_policy = policy))
  }
  unsupported <- x$fit; unsupported$readiness$fit$CategoryState <- "unsupported_coordinate"
  expect_error(predict_mfrm_units(unsupported, d), "not ready for fitted-object scoring")
})

test_that("two-family scoring enforces the retained response and prior contract", {
  x <- gmfrm_scoring_fixture(); d <- x$data[x$data$Person == "NEW3", ]
  d$Task <- as.character(d$Task); d$Rater <- as.character(d$Rater)
  expect_error(predict_mfrm_units(x$fit, transform(d, Task = "unknown")), "unseen levels")
  expect_error(predict_mfrm_units(x$fit, transform(d, Score = 3)), "score support")
  expect_error(predict_mfrm_units(x$fit, transform(d, w = 2), weight = "w"), "unit observation weights")
  expect_error(predict_mfrm_units(x$fit, rbind(d, d[1, ])), "Repeated person-facet")
  expect_error(predict_mfrm_units(x$fit, d, scoring_prior = list(mean = 1, sd = 2)), "retains the fitted")
  expect_error(predict_mfrm_units(x$fit, d, person_data = data.frame(Person = "NEW3")), "person-level")
  expected <- predict_mfrm_units(x$fit, d[-1, ], scoring_quad_points = 121L)
  missing <- d; missing$Score[1] <- NA
  expect_warning(actual <- predict_mfrm_units(x$fit, missing, scoring_quad_points = 121L), "Dropped 1 row")
  expect_equal(actual$estimates, expected$estimates)
  expect_equal(actual$row_review$DroppedRows, 1L)
  expect_warning(expect_error(predict_mfrm_units(x$fit, transform(d, Score = NA)),
    "No valid rows"), "Dropped")
  pv <- sample_mfrm_plausible_values(x$fit, d, scoring_quad_points = 121L, n_draws = 3, seed = 7)
  expect_identical(pv$settings$source_scoring_status, "conditional")
  expect_true(all(pv$values$ScoreIntegrationReady))
  expect_s3_class(summary(pv), "summary.mfrm_plausible_values")
  expect_identical(extract_mfrm_calibration(x$fit)$header$schema_version, 6L)
})

test_that("new crossings retain owner identities and arbitrary column mappings", {
  x <- gmfrm_scoring_fixture()
  d <- x$fit$gmfrm$specification$data
  d <- d[!(d$Task == "t3" & d$Rater == "r3"), ]
  d$Task <- sub("t", "level:", as.character(d$Task), fixed = TRUE)
  d$Rater <- sub("r", "level:", as.character(d$Rater), fixed = TRUE)
  names(d) <- c("Candidate", "観点 名", "Judge-ID", "Rating")
  p <- mfrm_gmfrm_problem(d, 2L, gauss_hermite_normal(61L),
    slope_facets = names(d)[2:3], person = "Candidate", score = "Rating")
  fit <- mfrm_gmfrm_fit_result(p, mfrm_gmfrm_em(p, start = x$fit$opt$par, maxit = 200L))
  # As in a public call whose facets argument has the reverse owner order.
  fit$config$source_columns <- list(person = "Candidate", facets = names(d)[3:2], score = "Rating")
  new <- data.frame(Candidate = c("NEWa", "NEWa", "NEWb"),
    first = c("level:3", "level:1", "level:3"),
    second = c("level:3", "level:2", "level:2"), Rating = c(2, 0, 1))
  names(new)[2:3] <- names(d)[2:3]
  pred <- predict_mfrm_units(fit, new, scoring_quad_points = 121L, readiness_policy = "review")
  expect_true(all(pred$estimates$ScoreIntegrationReady))
  expect_identical(pred$settings$two_family_calibration$slope_facets, names(d)[2:3])
  contexts <- pred$settings$two_family_calibration$observed_contexts
  expect_false(any(contexts[[1]] == "level:3" & contexts[[2]] == "level:3"))
  for (id in pred$estimates$Person) {
    reference_data <- new[new$Candidate == id, ]; reference_data$Score <- reference_data$Rating
    reference <- gmfrm_scoring_reference(fit, reference_data)
    expect_equal(unname(unlist(pred$estimates[pred$estimates$Person == id,
      c("Estimate", "SD", "Lower", "Upper")])), unname(reference), tolerance = 1e-7)
  }
  mapped <- new; names(mapped) <- c("ID", "first", "second", "Y")
  again <- predict_mfrm_units(fit, mapped, person = "ID", score = "Y",
    facets = setNames(c("second", "first"), names(d)[3:2]),
    scoring_quad_points = 121L, readiness_policy = "review")
  expect_identical(again$estimates, pred$estimates)
  expect_s3_class(summary(again), "summary.mfrm_unit_prediction")
  calibration <- freeze_mfrm_calibration(validate_mfrm_calibration(
    extract_mfrm_calibration(fit, scoring_quad_points = 121L)))
  portable <- score_mfrm_calibration(calibration, new)
  expect_identical(portable$row_dispositions$ObservedContext, c(FALSE, TRUE, TRUE))
  fields <- c("Estimate", "SD", "Lower", "Upper")
  expect_equal(unname(as.matrix(portable$estimates[fields])), unname(as.matrix(pred$estimates[fields])), tolerance = 1e-10)
  expect_identical(names(calibration$input_schema$source_columns$facets), fit$config$slope_facet)
  for (scores in list(pred, portable))
    expect_identical(mfrm_results(fit, scores = scores)$tables$person_scores, scores$estimates)
})

test_that("two-family saved fits and scores reopen in a fresh process", {
  skip_if_not_installed("callr")
  x <- gmfrm_scoring_fixture(); x$data <- x$data[x$data$Person == "NEW3", ]
  x$prediction <- predict_mfrm_units(x$fit, x$data, scoring_quad_points = 121L)
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE); saveRDS(x, file)
  worker <- function(root, file) {
    if (file.exists(file.path(root, "R", "api-prediction.R"))) pkgload::load_all(root, quiet = TRUE, compile = FALSE) else
      library("mfrmr", lib.loc = dirname(root), character.only = TRUE)
    testthat::local_mocked_bindings(fit_mfrm = function(...) stop("no refitting"),
      mfrm_gmfrm_em = function(...) stop("no refitting"),
      compute_mml_parameter_covariance = function(...) stop("no covariance"), .package = "mfrmr")
    x <- readRDS(file)
    scored <- mfrmr::predict_mfrm_units(x$fit, x$data, scoring_quad_points = 121L)
    testthat::local_mocked_bindings(mfrmr_adaptive_person_kernel = function(...) stop("no integration"),
      .package = "mfrmr")
    list(scored = scored, summary = summary(x$prediction))
  }
  environment(worker) <- baseenv()
  restored <- callr::r(worker, list(normalizePath(find.package("mfrmr")), file), libpath = .libPaths())
  expect_identical(restored$scored, x$prediction)
  expect_identical(restored$summary, summary(x$prediction))
})
