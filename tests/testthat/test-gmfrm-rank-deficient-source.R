test_that('a real N20 source with deficient scores has positive observed information', {
  saved <- readRDS(test_path('fixtures','gmfrm-rank-deficient-source.rds'))
  fit <- saved$fit
  before <- fit
  expect_length(fit$prep$levels$Person,20L)
  expect_length(fit$opt$par,22L)
  joint <- mfrm_gpcm_product_inference(fit)
  expect_true(joint$check$eligible)
  expect_true(joint$numerical_checks$ScoreRankDeficient)
  expect_lte(joint$numerical_checks$ScoreRank,20L)
  expect_true(all(joint$checks$Passed))
  expect_gt(min(eigen(saved$information,symmetric=TRUE,only.values=TRUE)$values),0)
  # The independent Louis calculation uses the refined integration order;
  # compare in information-scaled units under the recorded 1e-3 criterion.
  C <- chol(joint$information$hessian)
  expect_lt(norm(C %*% (solve(saved$information)-joint$covariance) %*% t(C),'2'),1e-3)
  local_mocked_bindings(mfrm_gpcm_product_inference=function(...)joint,.package='mfrmr')
  expect_warning(slopes <- confint(fit, parm = "slopes"),'Observed Person-score rank is deficient')
  expect_true(all(is.finite(slopes)))
  expect_warning(locations <- mfrm_facet_intervals(fit,'Rater'),
    'Observed Person-score rank is deficient')
  expect_true(all(is.finite(locations$table$Lower)))
  expect_match(plot_data(plot(locations,draw=FALSE))$subtitle,'Deficient Person-score rank')
  expect_match(paste(plot_data(as_ggplot(locations), "notes")$Text, collapse = "\n"),'Deficient Person-score rank')
  expect_null(plot_data(plot(locations,draw=FALSE,subtitle=NULL))$subtitle)
  expect_identical(fit,before)
})
