test_that("ordinary retained-prior scoring checks the batch independently of calibration readiness", {
  x <- readRDS(test_path("fixtures","ordinary-score-integration.rds"))
  for (fit in x$fits) {
    before <- serialize(fit,NULL)
    expect_true(mfrmr:::prediction_source_scoring_readiness(fit)$ready)
    expect_error(predict_mfrm_units(fit,x$new_data),"Posterior scoring integration did not pass")
    review <- predict_mfrm_units(fit,x$new_data,readiness_policy="review")
    expect_true(review$settings$score_integration_required)
    expect_false(all(review$estimates$ScoreIntegrationReady))
    expect_true(any(review$estimates$EstimateUse=="review_only_scoring_integration"))
    expect_s3_class(summary(review),"summary.mfrm_unit_prediction")
    retained <- predict_mfrm_units(fit,x$new_data,scoring_quad_points=121L,n_draws=2,seed=17)
    reference <- predict_mfrm_units(fit,x$new_data,scoring_quad_points=121L,
      scoring_prior=list(mean=0,sd=1),n_draws=2,seed=17)
    expect_true(all(retained$estimates$ScoreIntegrationReady))
    expect_true(all(retained$draws$ScoreIntegrationReady))
    expect_true(all(retained$estimates$PriorSource=="retained"))
    expect_true(all(reference$estimates$PriorSource=="user_supplied"))
    expect_equal(retained$estimates[c("Estimate","SD","Lower","Upper")],
      reference$estimates[c("Estimate","SD","Lower","Upper")],tolerance=1e-13)
    expect_identical(retained$draws$Value,reference$draws$Value)
    expect_identical(serialize(fit,NULL),before)
    expect_identical(summary(retained)$settings,retained$settings)
    stale <- retained; stale$settings$score_integration_review <- NULL
    expect_error(summary(stale),"inconsistent numerical review records")
    expect_error(mfrmr:::export_validate_optional_object(stale,"mfrm_unit_prediction","prediction", fit = fit),
      "inconsistent numerical review records")
    legacy <- retained
    legacy$settings$score_integration_required <- legacy$settings$score_integration_review <- NULL
    legacy$estimates$ScoreIntegrationReady <- NULL
    expect_s3_class(summary(legacy),"summary.mfrm_unit_prediction")
  }
})
