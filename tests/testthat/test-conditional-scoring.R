conditional_scoring_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      x <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))
      x$new <- x$data[x$data$Person %in% head(unique(x$data$Person), 2), ]
      # Retain the formula environment so an unchanged fit has stable identity.
      cached <<- x
    }
    local_cran_fixture_information(cached$fit, "conditional-scoring", .env = parent.frame())
    cached
  }
})

test_that("conditional population scoring carries checks through summaries, draws and exports", {
  x <- conditional_scoring_fixture()
  before <- x$fit$readiness
  pred <- predict_mfrm_units(x$fit,x$new,n_draws=3,seed=42)
  expect_true(pred$settings$source_scoring_ready)
  expect_identical(pred$settings$source_scoring_status,"conditional")
  expect_true(pred$settings$local_calibration_review$eligible)
  expect_true(all(pred$settings$score_integration_review$Passed))
  expect_true(all(pred$estimates$ScoreIntegrationReady))
  expect_true(all(pred$draws$ScoreIntegrationReady))
  expect_identical(x$fit$readiness,before)
  expect_false(pred$settings$source_inference_ready)
  expect_identical(pred$settings$source_audit_states[["boundary"]],"not_evaluated")
  expect_identical(summary(pred)$settings,pred$settings)
  stale <- pred
  stale$settings$score_integration_review <- NULL
  expect_error(summary(stale),"inconsistent numerical review records")
  expect_error(mfrmr:::export_validate_optional_object(stale,"mfrm_unit_prediction","prediction", fit = x$fit),
    "inconsistent numerical review records")
  expect_s3_class(build_summary_table_bundle(pred),"mfrm_summary_table_bundle")
  expect_silent(mfrmr:::export_validate_optional_object(pred,"mfrm_unit_prediction","prediction", fit = x$fit))
  expect_match(paste(capture.output(print(summary(pred))),collapse=" "),"local solution")
  pv <- sample_mfrm_plausible_values(x$fit,x$new,n_draws=3,seed=42)
  expect_identical(pv$settings$source_scoring_status,"conditional")
  expect_true(all(pv$values$ScoreIntegrationReady))
  expect_s3_class(summary(pv),"summary.mfrm_plausible_values")
})

test_that("source qualification cannot authorize an unresolved scoring grid", {
  x <- conditional_scoring_fixture()
  expect_error(predict_mfrm_units(x$fit,x$new,scoring_quad_points=2),
    "Posterior scoring integration did not pass")
  pred <- predict_mfrm_units(x$fit,x$new,scoring_quad_points=2,readiness_policy="review")
  expect_true(pred$settings$source_scoring_ready)
  expect_false(all(pred$estimates$ScoreIntegrationReady))
  expect_true(any(pred$estimates$EstimateUse=="review_only_scoring_integration"))
  expect_s3_class(summary(pred),"summary.mfrm_unit_prediction")
  expect_match(paste(capture.output(print(summary(pred))),collapse=" "),"Integration review required")
})

test_that("invalid calibrations and stale scoring priors cannot be enabled by review", {
  x <- conditional_scoring_fixture()
  for(change in c("parameter","variance","mean")) {
    fit <- x$fit
    if(change=="parameter")fit$opt$par[1]<-Inf
    if(change=="variance")fit$population$sigma2<-0
    if(change=="mean")fit$population$coefficients[1]<-fit$population$coefficients[1]+0.1
    for(policy in c("error","review")) expect_error(
      predict_mfrm_units(fit,x$new,readiness_policy=policy),"calibration or scoring prior is invalid")
  }
})

test_that("actual source failures stay ineligible and unevaluated audits are labelled accurately", {
  x <- conditional_scoring_fixture()
  for(component in c("NumericalState","EstimabilityState","CategoryState","BoundaryState")) {
    fit <- x$fit
    fit$readiness$fit[[component]] <- switch(component,
      NumericalState="failed",EstimabilityState="structurally_unidentified",
      CategoryState="unsupported_coordinate",BoundaryState="has_exclusions")
    expect_false(mfrmr:::prediction_source_scoring_readiness(fit)$ready)
    expect_error(predict_mfrm_units(fit,x$new),"not ready for fitted-object scoring")
  }
  reasons <- mfrmr:::prediction_scoring_reasons(c(
    "boundary_not_scoring_ready:not_evaluated","boundary_not_scoring_ready:has_exclusions"))
  expect_match(reasons[1],"has not been evaluated")
  expect_match(reasons[2],"boundary or be unbounded")
  failed <- x$fit
  failed$opt$par[1]<-failed$opt$par[1]+0.1
  expect_false(mfrmr:::prediction_source_scoring_readiness(failed)$ready)
})
