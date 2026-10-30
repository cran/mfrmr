test_that("slope intervals need fresh qualification and transport joint covariance", {
  fit <- list(config=list(model='GPCM',method='MML'),opt=list(par=0),steps=data.frame(),
    slopes=data.frame(Estimate=c(2,.5),LogEstimate=log(c(2,.5))))
  # The slope covariance is a block of the inverse JOINT information.
  # With nuisance correlation it differs from inverting the slope block alone.
  info <- matrix(c(25,10,10,20),2)
  covariance <- list(status='ok',cov=solve(info),param_slices=list(log_slopes=1L),
    detail='Known joint covariance',rank=2L,regularized=FALSE)
  qualified <- FALSE
  local_mocked_bindings(mfrm_gpcm_slope_inference_check=function(...) {
    list(eligible=qualified,review='Controlled qualification test')
  },.package='mfrmr')
  uncertainty <- function(object=fit, cv=covariance) {
    mfrmr:::compute_mml_structural_parameter_se(object,covariance=cv)$slopes
  }
  withheld <- function(out) {
    for (field in c('SE','LogSE','CI_Lower','CI_Upper','LogCI_Lower','LogCI_Upper')) {
      expect_true(all(is.na(out[[field]])))
    }
    expect_false(any(out$SEEligible)); expect_false(any(out$CIEligible))
  }
  legacy <- uncertainty(); withheld(legacy)
  se <- sqrt(solve(info)[1,1])
  expect_equal(legacy$OptimizerSE,c(2,.5)*se)
  expect_equal(legacy$OptimizerLogSE,rep(se,2))
  expect_false(isTRUE(all.equal(se,sqrt(1/info[1,1]))))
  for (version in c('legacy_or_unknown',mfrmr:::mfrmr_readiness_contract_version())) {
    fit$slopes$SEEligible <- fit$slopes$CIEligible <- TRUE
    fit$slopes$ReadinessContractVersion <- version
    withheld(uncertainty()) # Flags, even current ones, are not qualification.
  }
  qualified <- TRUE
  admitted <- uncertainty()
  expect_equal(admitted$SE,c(2,.5)*se)
  expect_equal(admitted$CI_Lower,c(2,.5)*exp(-qnorm(.975)*se))
  expect_equal(admitted$CI_Upper,c(2,.5)*exp(qnorm(.975)*se))
  expect_equal(admitted$PrimaryEstimate,c(2,.5))
  expect_true(all(admitted$CIEligible))
  expect_true(all(admitted$CIUse=='approximate_pointwise'))
  # Bounds that overflow cannot be reported, even with an otherwise valid SE.
  large <- covariance; large$cov <- diag(c(1e6,1))
  wide <- uncertainty(cv=large)
  expect_true(all(wide$SEEligible)); expect_false(any(wide$CIEligible))
  expect_true(all(is.na(wide$CI_Lower)))
  expect_true(all(grepl('not representable',wide$InferenceReview)))
  covariance$status <- 'regularized'
  withheld(uncertainty())
  covariance$status <- 'fallback'; covariance$cov <- NULL
  withheld(uncertainty())
})

test_that("slope qualification checks coordinates instead of stored interval flags", {
  spec <- mfrmr:::build_gpcm_slope_spec(c('a','b'),'Rater','Rater')
  fit <- list(config=list(model='GPCM',method='MML',gpcm_spec=spec,
    slope_facet='Rater',step_facet='Rater',facet_names='Rater',n_cat=2L,
    facet_levels=list(Rater=c('a','b'))),opt=list(par=log(2)),
    slopes=data.frame(SlopeFacet=c('a','b'),Estimate=c(2,.5),LogEstimate=log(c(2,.5))))
  cv <- list(status='ok',regularized=FALSE,param_slices=list(log_slopes=1L),
    solution_information=list(status='test'))
  calls <- 0L
  local_mocked_bindings(mfrm_extract_fit_ic_contract=function(...) data.frame(),
    mfrm_ic_fit_check=function(fit,contract,information,allow_singleton) {
      expect_false(allow_singleton)
      calls <<- calls+1L; expect_identical(information,cv$solution_information)
      list(eligible=TRUE,review='')
    },.package='mfrmr')
  check <- function(f=fit,v=cv) mfrmr:::mfrm_gpcm_slope_inference_check(f,v)
  expect_true(check()$eligible)
  old <- fit; old$slopes$Estimate[1] <- 3
  expect_false(check(old)$eligible)
  old <- fit; old$slopes <- old$slopes[2:1,]
  expect_false(check(old)$eligible)
  old <- fit; old$config$gpcm_spec$n_params <- 0L
  expect_false(check(old)$eligible)
  old <- fit; old$config$method <- 'JML'
  expect_false(check(old)$eligible)
  old <- fit; old$config$step_facet <- 'Criterion'
  expect_false(check(old)$eligible)
  old <- cv; old$solution_information <- NULL
  expect_false(check(v=old)$eligible)
  old <- cv; old$regularized <- TRUE
  expect_false(check(v=old)$eligible)
  expect_equal(calls,1L)
})

test_that("GPCM guidance separates global readiness from output-specific inference", {
  readiness <- data.frame(ReadinessContractVersion=mfrmr:::mfrmr_readiness_contract_version(),
    InputState='pass',EstimabilityState='not_evaluated',CategoryState='adequate',
    BoundaryState='not_evaluated',NumericalState='ready',FitReadiness='review',
    InferenceReady=FALSE,ReasonCodes='mml_gpcm_slope_boundary_not_evaluated')
  decision <- mfrmr:::mfrm_fit_decision_summary(readiness)
  expect_identical(decision$FormalInference,'No')
  expect_match(decision$Why,'confint(fit',fixed=TRUE)
  expect_match(decision$Why,'does not certify global boundary',fixed=TRUE)
  readiness$ReasonCodes <- 'boundary_audit_incomplete'
  expect_false(grepl('GPCM',mfrmr:::mfrm_fit_decision_summary(readiness)$Why))
})

test_that("confint uses the common calculation and validates the request", {
  fit <- structure(list(config=list(model='GPCM',method='MML',
    slope_facet='Rater',step_facet='Rater',facet_names='Rater',n_cat=2L,
    facet_levels=list(Rater=c('a','b')))),class='mfrm_fit')
  before <- fit
  local_mocked_bindings(compute_mml_structural_parameter_se=function(res,ci_level) {
    list(slopes=data.frame(SlopeFacet=c('a','b'),Estimate=c(1,1),CI_Lower=c(ci_level,NA),CI_Upper=c(2,NA),
      CIEligible=c(TRUE,FALSE),InferenceReview=c('Approximation','Unavailable')))
  },.package='mfrmr')
  ci <- confint(fit,parm='slopes',level=.90)
  expect_equal(unname(ci[1,]),c(.9,2))
  expect_true(all(is.na(ci[2,])))
  expect_identical(attr(ci,'level'),.9)
  expect_identical(rownames(ci),c('a','b'))
  expect_identical(fit,before)
  expect_s3_class(ci,'mfrm_slope_intervals')
  display <- capture.output(print(ci))
  expect_true(any(grepl('Approximate 90%',display,fixed=TRUE)))
  expect_true(any(grepl('Unavailable intervals: 1',display,fixed=TRUE)))
  expect_false(any(grepl('InferenceReview|ReadinessContractVersion',display)))
  expect_match(attr(ci,'target'),'geometric mean one')
  for (level in list(0,1,NA,Inf,c(.9,.95),'95%')) expect_error(confint(fit, parm = "slopes",level=level),'level')
  expect_error(confint(fit, parm = "slopes",unknown=TRUE),'must be empty')
  fit$config$method <- 'JML'; expect_error(confint(fit, parm = "slopes"),'MML')
  fit$config$model <- 'PCM'; expect_error(confint(fit, parm = "slopes"),'GPCM')
})
