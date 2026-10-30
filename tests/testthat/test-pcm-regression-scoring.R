test_that("qualified PCM regression scores use the fitted conditional population", {
  fit <- readRDS(test_path("fixtures", "native-location-pcm-regression.rds"))$fit
  before <- fit$readiness
  persons <- fit$population$person_table[1:3, ]
  rows <- fit$prep$data[fit$prep$data$Person %in% persons$Person,
    c("Person", "Rater", "Criterion", "Score")]
  rows$Person <- paste0("new-", rows$Person)
  persons$Person <- paste0("new-", persons$Person)
  local_mocked_bindings(fit_mfrm=function(...) stop("No refitting"), .package="mfrmr")
  pred <- predict_mfrm_units(fit, rows, person_data=persons, scoring_quad_points=61)
  expect_true(pred$settings$source_scoring_ready)
  expect_identical(pred$settings$source_scoring_status,"conditional")
  expect_true(all(pred$settings$score_integration_review$Passed))
  expect_true(all(pred$estimates$EstimateUse=="fitted_object_scoring"))
  expect_equal(pred$estimates$PriorMean,
    unname(fit$population$coefficients[1]+fit$population$coefficients[2]*
      persons$X[match(pred$estimates$Person,persons$Person)]))
  expect_equal(pred$estimates$PriorSD,rep(sqrt(fit$population$sigma2),3))
  expect_identical(fit$readiness,before)
  expect_false(pred$settings$source_inference_ready)
  expect_identical(summary(unserialize(serialize(pred,NULL))),summary(pred))
  expect_silent(export_validate_optional_object(pred,"mfrm_unit_prediction","prediction", fit = fit))
})

test_that("regression scoring retains source and scope refusals", {
  fit <- readRDS(test_path("fixtures", "native-location-pcm-regression.rds"))$fit
  bad <- fit; bad$config$population_spec$design_matrix[,2] <- 1
  expect_false(prediction_source_scoring_readiness(bad)$ready)
  bad <- fit; bad$facets$others$Estimate[1] <- bad$facets$others$Estimate[1]+.1
  check <- prediction_source_scoring_readiness(bad)
  expect_false(check$ready)
  expect_match(check$local_calibration_review$review,"Saved locations")
  bad <- fit; bad$readiness$fit$BoundaryState <- "has_exclusions"
  expect_false(prediction_source_scoring_readiness(bad)$ready)
  bad <- fit; bad$config$model <- "GPCM"
  expect_false(prediction_source_scoring_readiness(bad)$ready)
})
