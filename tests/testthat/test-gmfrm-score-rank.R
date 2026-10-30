gmfrm_rank_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      saved <- readRDS(test_path('fixtures', 'gmfrm-joint-information.rds'))
      problem <- mfrm_gmfrm_problem(saved$data, 2L, gauss_hermite_normal(31L))
      fit <- mfrm_gmfrm_fit_result(problem,
        mfrm_gmfrm_em(problem, start=saved$parameters, maxit=1L))
      cached <<- list(fit=fit, joint=mfrm_gpcm_product_inference(fit))
    }
    cached
  }
})

test_that('empirical score rank is a diagnostic, not an observed-Hessian veto', {
  x <- gmfrm_rank_fixture()
  original_scores <- mfrm_mml_person_scores_numeric
  # Controlled unit counterexample: pool scores into two groups, preserving
  # their total (stationarity/integration) while making empirical rank <= 2.
  # This isolates the gate; actual saved N20 fits are checked independently.
  pooled_scores <- function(...) {
    scores <- original_scores(...)
    group <- seq_len(nrow(scores)) %% 2L
    for (g in 0:1) scores[group==g, ] <- rep(colMeans(scores[group==g, ,drop=FALSE]),
      each=sum(group==g))
    scores
  }
  local_mocked_bindings(mfrm_mml_person_scores_numeric=pooled_scores, .package='mfrmr')
  joint <- mfrm_gpcm_product_inference(x$fit)
  expect_true(joint$check$eligible)
  expect_true(joint$numerical_checks$ScoreRankDeficient)
  expect_lte(joint$numerical_checks$ScoreRank,2L)
  expect_lte(joint$numerical_checks$FineScoreRank,2L)
  expect_identical(joint$covariance,x$joint$covariance)
  expect_match(joint$check$caution,'Observed Person-score rank is deficient')
  expect_true(all(joint$checks$Passed))
  # Exercise both real consumers with the same verified result, without
  # recomputing information just to test reporting and persistence.
  local_mocked_bindings(mfrm_gpcm_product_inference=function(...)joint,.package='mfrmr')
  expect_warning(slopes <- confint(x$fit, parm = "slopes"),'Observed Person-score rank is deficient')
  expect_warning(locations <- mfrm_facet_intervals(x$fit,'Rater'),
    'Observed Person-score rank is deficient')
  expect_output(print(slopes),'Observed Person-score rank is deficient')
  expect_output(print(locations),'Observed Person-score rank is deficient')
  expect_match(paste(plot_data(as_ggplot(slopes), "notes")$Text, collapse = "\n"),'Observed Person-score rank is deficient')
  res <- mfrm_results(x$fit,include='fit',compute='never',
    intervals=list(slopes=slopes,locations=locations))
  expect_match(mfrm_report(res)$markdown,'Observed Person-score rank is deficient')
  expect_identical(res$tables$gpcm_slopes_numerical_checks,joint$numerical_checks)
  expect_identical(res$tables$gpcm_locations_numerical_checks,joint$numerical_checks)
  path <- tempfile(fileext='.rds'); withr::defer(unlink(path))
  saveRDS(res,path)
  expect_identical(readRDS(path),res)
  expect_match(mfrm_report(readRDS(path))$markdown,'Observed Person-score rank is deficient')
})

test_that('relaxing score rank leaves derivative and information safeguards intact', {
  x <- gmfrm_rank_fixture()
  original_scores <- mfrm_mml_person_scores_numeric
  local_mocked_bindings(mfrm_mml_person_scores_numeric=function(fit,relative_step=1e-5) {
    scores <- original_scores(fit,relative_step=relative_step)
    if (relative_step < 1e-5) scores[1,1] <- scores[1,1]+1
    scores
  },.package='mfrmr')
  bad <- mfrm_gpcm_product_inference(x$fit)
  expect_false(bad$check$eligible)
  expect_null(bad$covariance)
  expect_identical(tail(bad$checks$Check,1),'Score derivative stability')
  expect_gt(bad$numerical_checks$DerivativeDifference,1e-6)
  for (status in c('unavailable','ok')) {
    info <- x$joint$information; info$status <- status; info$regularized <- TRUE
    local_mocked_bindings(compute_mml_parameter_covariance=function(...)info,.package='mfrmr')
    bad <- mfrm_gpcm_product_inference(x$fit)
    expect_false(bad$check$eligible)
    expect_null(bad$covariance)
    expect_identical(tail(bad$checks$Check,1),'Joint information')
  }
})
