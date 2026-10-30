anchored_gpcm_fixture <- local({
  cache <- NULL
  function(case) {
    if(is.null(cache)) cache <<- readRDS(test_path('fixtures','anchored-gpcm.rds'))
    local_cran_fixture_information(cache[[case]]$fit, paste0("anchor-", case), .env = parent.frame())
    cache[[case]]
  }
})

for (case in c('shared','separate_group','fixed_population')) {
  test_that(paste(case,'GPCM keeps anchored coordinates in joint uncertainty'), {
    # Grouped/separate ownership is the live joint-information representative.
    if (case != 'separate_group') skip_on_cran()
    x <- anchored_gpcm_fixture(case); fit <- x$fit; ref <- x$reference
    original_parameters <- fit$opt$par
    information <- compute_mml_parameter_covariance(fit)
    expect_identical(information$status,'ok')
    expect_lt(max(abs(information$cov-ref$covariance))/max(abs(ref$covariance)),1e-5)
    expect_equal(unname(fit$facets$others$Estimate),c(x$coefficients$rater,x$coefficients$criterion),tolerance=1e-12)
    expect_equal(sum(x$coefficients$rater),0,tolerance=1e-12)
    # Reuse fresh derivatives across targets; independent covariance checks stay live.
    scores <- mfrm_mml_person_scores_numeric(fit)
    source_fit <- fit
    local_mocked_bindings(compute_mml_parameter_covariance=function(res, ...) {
      stopifnot(identical(res, source_fit))
      information
    }, mfrm_mml_person_scores_numeric=function(fit, relative_step=1e-5) {
      stopifnot(identical(fit, source_fit), identical(relative_step, 1e-5))
      scores
    },
      fit_mfrm=function(...) stop('No refitting is required'),.package='mfrmr')
    default <- confint(fit, parm = "slopes")
    expect_true(all(attr(default,'diagnostics')$CIEligible))
    expect_match(attr(default,'anchor_note'),'Conditional on supplied anchors',fixed=TRUE)
    expect_true(all(grepl('source calibration',attr(default,'diagnostics')$InferenceReview,fixed=TRUE)))
    expect_output(print(default),'source calibration')
    expect_true(any(grepl('source calibration',plot_data(plot(default,draw=FALSE),'notes')$Text,fixed=TRUE)))
    for(method in c('model','sandwich')) {
      covariance <- if(method=='model') ref$covariance else ref$sandwich
      J <- ref$standardized_jacobian
      expected <- J %*% covariance %*% t(J)
      intervals <- confint(fit, parm = "slopes", scale='standardized',method=method)
      tab <- attr(intervals,'diagnostics')
      expect_true(all(tab$CIEligible))
      expect_lt(max(abs(attr(intervals,'covariance')-expected))/max(abs(expected)),2e-5)
      estimate <- x$coefficients$slopes*x$coefficients$sd
      expect_equal(tab$Estimate,estimate,tolerance=1e-12)
      expect_lt(max(abs(as.numeric(intervals[,1])-estimate*exp(-qnorm(.975)*sqrt(diag(expected))))),1e-5)
      expect_true(all(grepl('Conditional on supplied anchors',tab$InferenceReview,fixed=TRUE)))
      curves <- mfrm_curve_intervals(fit,ref$grid,method=method)
      expect_true(all(curves$table$CIEligible))
      expect_lt(max(abs(curves$table$Estimate-ref$probability)),1e-12)
      variance <- rowSums((ref$probability_jacobian %*% covariance)*ref$probability_jacobian)
      expect_lt(max(abs(curves$table$SE-sqrt(variance))),1e-5)
      expect_output(print(curves),'source calibration')
      p <- plot(curves,draw=FALSE)
      expect_null(p$labels$title); expect_null(p$labels$caption)
      expect_true(any(grepl('source calibration',plot_data(p,'notes')$Text,fixed=TRUE)))
      expect_true(all(grepl('source calibration',curves$table$InferenceReview,fixed=TRUE)))
      # Independent slope-squared category-variance information and its derivative.
      probs <- matrix(ref$probability,ncol=4,byrow=TRUE)
      means <- drop(probs %*% (0:3))
      v <- drop(probs %*% (0:3)^2)-means^2
      owner <- match(ref$grid$Criterion,fit$config$facet_levels$Criterion)
      a <- x$coefficients$slopes[owner]
      I <- a^2*v
      slope_J <- ref$standardized_jacobian
      pop <- information$param_slices$log_sigma2
      if(length(pop)) slope_J[,pop] <- 0
      derivative <- t(vapply(seq_len(nrow(probs)),function(i) {
        dp <- ref$probability_jacobian[(4*i-3):(4*i),,drop=FALSE]
        a[i]^2*drop(crossprod((0:3)^2-2*means[i]*(0:3),dp)) +
          2*I[i]*slope_J[owner[i],]
      },numeric(ncol(covariance))))
      curve_info <- mfrm_curve_intervals(fit,ref$grid,type='information',method=method)
      expect_lt(max(abs(curve_info$table$Estimate-I)),1e-12)
      expect_lt(max(abs(curve_info$table$SE-sqrt(rowSums((derivative %*% covariance)*derivative)))),1e-5)
    }
    diagnostics <- diagnose_mfrm(fit,residual_pca='none')
    rater <- subset(diagnostics$measures,Facet=='Rater')
    fixed <- if(case=='separate_group') c(TRUE,FALSE,FALSE,TRUE) else c(TRUE,FALSE,FALSE,FALSE)
    expect_identical(rater$Fixed,fixed)
    expect_true(all(is.na(rater$ModelSE[fixed])))
    expect_true(all(is.na(subset(diagnostics$reliability,Facet=='Rater')$Reliability)))
    expect_identical(fit$opt$par,original_parameters)
    expect_false(isTRUE(fit$summary$InferenceReady))
    expect_error(mfrm_facet_intervals(fit,facet='Rater'),'without anchors or groups')
    expect_error(confint(fit, parm = "slopes", method='profile',slope=fit$slopes$SlopeFacet[1]),'no anchors')
    expect_error(extract_mfrm_calibration(fit),if(case=='fixed_population')
      'SCORING_BASIS_UNSUPPORTED' else 'portable GPCM does not yet support anchors')
    saved <- unserialize(serialize(default,NULL))
    expect_identical(saved,default)
  })
}

test_that('GPCM refit controls preserve anchors, group totals and centering', {
  for(case in c('shared','separate_group','fixed_population')) {
    fit <- anchored_gpcm_fixture(case)$fit
    data <- load_mfrmr_data('example_core')
    calls <- list()
    local_mocked_bindings(fit_mfrm=function(...) {
      calls[[length(calls)+1L]] <<- list(...); fit
    },.package='mfrmr')
    replay <- mfrm_gpcm_bootstrap_refit(fit,data)
    expect_identical(replay$opt$par,fit$opt$par)
    call <- calls[[1]]
    expect_identical(call$centering_scope,'facet')
    expect_equal(call$anchors$Anchor[call$anchors$Facet=='Rater'],.7)
    expect_identical(call$step_facet,fit$config$step_facet)
    expect_identical(call$slope_facet,fit$config$slope_facet)
    if(case=='separate_group') expect_equal(call$group_anchors$GroupValue,rep(.1,3))
  }
})

test_that('unsupported extension anchors are rejected before optimization', {
  ratings <- load_mfrmr_data('example_core')
  anchor <- data.frame(Facet='Rater',Level='R01',Anchor=.7)
  expect_error(fit_mfrm(ratings,'Person',c('Rater','Criterion'),'Score',
    method='JML',model='GPCM',step_facet='Criterion',jml_correction_order=1,
    anchors=anchor),'not supported by corrected JML')
  expect_error(fit_mfrm(ratings,'Person',c('Rater','Criterion'),'Score',
    method='MML',model='GPCM',slope_facet=c('Criterion','Rater'),step_facet='Rater',
    noncenter_facet='Rater',gpcm_mml_identification='fixed_standard_normal',mml_engine='em',
    anchors=anchor),'not supported by the two-family route')
  expect_error(fit_mfrm_testlet(anchors=anchor),'unused argument')
  expect_error(fit_mfrm_random_rater(anchors=anchor),'unused argument')
  rsm <- readRDS(test_path('fixtures','anchored-facet-information.rds'))$RSM$fit
  for(format in c('stan','brms')) {
    destination <- tempfile()
    on.exit(unlink(destination, recursive=TRUE), add=TRUE)
    expect_silent(exported <- export_mfrm_model(rsm,destination,format=format,acknowledge_sensitive=TRUE))
    bundle <- readRDS(exported$files[['model-data.rds']])
    expect_identical(bundle$constraints$facets,rsm$config$facet_specs)
  }
})

test_that('anchored GPCM scoring preserves conditional values and source-review limits', {
  local_mocked_bindings(fit_mfrm=function(...) stop('No refitting'),.package='mfrmr')
  for(case in c('shared','separate_group','fixed_population')) {
    x <- anchored_gpcm_fixture(case); spec <- x$scoring
    if(!spec$ready) expect_error(predict_mfrm_units(x$fit,spec$rows,scoring_quad_points=161L),
      'not ready for fitted-object scoring')
    score <- predict_mfrm_units(x$fit,spec$rows,scoring_quad_points=161L,
      readiness_policy=if(spec$ready) 'error' else 'review')
    expect_lt(max(abs(c(score$estimates$Estimate,score$estimates$SD)-spec$reference)),1e-7)
    expect_identical(score$estimates$SourceScoringReady,spec$ready)
    expect_identical(score$estimates$EstimateUse,if(spec$ready)
      'fitted_object_scoring' else 'review_only_nonready_source')
    expect_false(isTRUE(x$fit$summary$InferenceReady))
  }
})
