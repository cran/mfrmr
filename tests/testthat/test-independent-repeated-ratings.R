independent_rating_fixture <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      withr::local_seed(76013)
      d <- expand.grid(Person=paste0('P',1:40),Rater=c('r1','r2'),
        Criterion=c('c1','c2'),Event=1:2,stringsAsFactors=FALSE)
      theta <- rnorm(40)
      eta <- theta[match(d$Person,paste0('P',1:40))] -
        c(-.4,.4)[match(d$Rater,c('r1','r2'))] -
        c(-.25,.25)[match(d$Criterion,c('c1','c2'))]
      d$Score <- rbinom(nrow(d),1,plogis(eta))
      value <<- d
    }
    value
  }
})

test_that('independent repeated ratings preserve events and other input reviews', {
  d <- independent_rating_fixture()
  args <- list(data=d,person_col='Person',facet_cols=c('Rater','Criterion'),
    score_col='Score',rating_min=0,rating_max=1,keep_original=TRUE)
  expect_warning(a <- do.call(prepare_mfrm_data,args),'duplicated Person')
  expect_warning(b <- do.call(prepare_mfrm_data,c(args,list(repeated_ratings='independent'))),NA)
  expect_identical(a$data,b$data)
  expect_identical(a$levels,b$levels)
  expect_identical(b$n_obs,nrow(d))
  expect_identical(b$repeated_ratings,'independent')
  expect_true('independent_repeated_ratings_declared' %in% b$preparation_notes$Condition)
  expect_false('duplicate_person_facet_cells' %in% b$preparation_notes$Condition)
  args$data$Score[1] <- NA_real_
  omitted <- do.call(prepare_mfrm_data,c(args,list(repeated_ratings='independent')))
  expect_true('missing_or_nonpositive_weight_rows_dropped' %in% omitted$preparation_notes$Condition)
  expect_identical(mfrmr_readiness_input_component(omitted,list(),NULL)$State,'review')
  expect_error(do.call(prepare_mfrm_data,c(args,list(repeated_ratings='dependent'))),'arg')
})

test_that('data description retains repeated-cell counts without a false dependence conclusion', {
  d <- independent_rating_fixture()
  z <- describe_mfrm_data(d,'Person',c('Rater','Criterion'),'Score',
    rating_min=0,rating_max=1,include_agreement=FALSE,repeated_ratings='independent')
  expect_true(z$duplicate_cell_summary$HasDuplicates)
  expect_equal(z$duplicate_cell_summary$ExtraRows,nrow(d)/2)
  expect_identical(z$repeated_ratings,'independent')
  caveats <- summary(z)$caveats
  expect_false('duplicate_person_facet_cells' %in% caveats$Condition)
  expect_true('independent_repeated_ratings_declared' %in% caveats$Condition)
})

test_that('the declaration changes ordinary admission without changing estimation', {
  d <- independent_rating_fixture()
  args <- list(data=d,person='Person',facets=c('Rater','Criterion'),score='Score',
    rating_min=0,rating_max=1,model='PCM',step_facet='Criterion',
    category_policy='preserve')
  for(method in c('MML','JML')) {
    args$method <- method
    a <- suppressWarnings(do.call(fit_mfrm,args))
    b <- suppressWarnings(do.call(fit_mfrm,c(args,list(repeated_ratings='independent'))))
    expect_identical(a$opt$par,b$opt$par)
    expect_identical(a$opt$value,b$opt$value)
    expect_identical(a$prep$data,b$prep$data)
    expect_identical(a$readiness$fit$InputState,'review')
    expect_identical(b$readiness$fit$InputState,'pass')
    expect_identical(b$config$replay_inputs$repeated_ratings,'independent')
    path <- tempfile(fileext='.rds');saveRDS(b,path);restored <- readRDS(path);unlink(path)
    expect_identical(restored,b)
    if(method=='MML') {
      expect_true(prediction_source_scoring_readiness(restored)$ready)
      rows <- d[d$Person %in% c('P1','P2','P3'),]
      scored <- predict_mfrm_units(restored,rows,scoring_quad_points=81L)
      expect_true(scored$settings$source_scoring_ready)
      expect_identical(scored$settings$source_repeated_ratings,'independent')
      expect_true(all(scored$settings$score_integration_review$Passed))
    }
  }
})

test_that('corrected JML accepts the same explicit declaration without changing its equation', {
  args <- list(data=independent_rating_fixture(),person='Person',facets=c('Rater','Criterion'),
    score='Score',rating_min=0,rating_max=1,model='GPCM',method='JML',
    step_facet='Criterion',slope_facet='Criterion',category_policy='preserve',jml_correction_order=2L)
  a <- suppressWarnings(do.call(fit_mfrm,args))
  b <- suppressWarnings(do.call(fit_mfrm,c(args,list(repeated_ratings='independent'))))
  expect_identical(a$jml_adjustment,b$jml_adjustment)
  expect_identical(a$prep$data,b$prep$data)
  expect_identical(b$config$repeated_ratings,'independent')
  expect_identical(b$config$replay_inputs$repeated_ratings,'independent')
})
