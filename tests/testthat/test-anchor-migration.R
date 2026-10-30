test_that('saved native constraints remain distinct in summaries and specifications', {
  old <- readRDS(test_path('fixtures','legacy-anchored-fit.rds'))
  current <- readRDS(test_path('fixtures','anchored-facet-information.rds'))
  fits <- c(list(legacy=old),lapply(current,`[[`,'fit'))
  local_mocked_bindings(fit_mfrm=function(...) stop('Do not refit saved summaries'),
    diagnose_mfrm=function(...) stop('No new diagnostics'),
    compute_mml_parameter_covariance=function(...) stop('No new information'),
    .package='mfrmr')
  for (name in names(fits)) {
    fit <- fits[[name]]
    before <- list(parameters=fit$opt$par,specs=fit$config$facet_specs,
      readiness=fit$readiness,signature=normalize_compare_signature(fit))
    expected <- if(name=='legacy') 'remainder' else 'facet'
    s <- summary(fit,compute='never',detail='brief')
    expect_identical(s$settings_overview$CenteringScope,expected)
    expect_identical(s$settings_overview$NoncenterFacet,fit$config$noncenter_facet)
    expect_output(print(s),if(name=='legacy') 'legacy ungrouped' else 'whole declared facet')
    spec <- specifications_report(fit,include_fixed=TRUE)
    expect_identical(spec$data_spec$Value[spec$data_spec$Setting=='CenteringScope'],expected)
    expect_match(paste(spec$fixed,collapse=' '),expected)
    expect_match(summarize_anchor_constraints(fit$config),if(name=='legacy')
      'Legacy centering' else 'including anchors and group totals',fixed=TRUE)
    bundle <- build_summary_table_bundle(s)
    expect_identical(bundle$tables$settings_overview$CenteringScope,expected)
    expect_identical(list(parameters=fit$opt$par,specs=fit$config$facet_specs,
      readiness=fit$readiness,signature=normalize_compare_signature(fit)),before)
  }
  expect_null(old$config$centering_scope)
  expect_null(old$config$replay_inputs$centering_scope)
  expect_true(all(vapply(old$config$facet_specs,function(x) is.null(x$centering_scope),TRUE)))
  r <- subset(old$facets$others,Facet=='Rater')
  expect_equal(sum(r$Estimate),.7,tolerance=1e-12)
  expect_lt(abs(sum(r$Estimate[r$Level!='R01'])),1e-12)
  expect_identical(mfrm_centering_scope(list()),'not_recorded')
  expect_match(summarize_anchor_constraints(list()),'not recorded',fixed=TRUE)
})

test_that('legacy anchored results export reopens without recomputing or relabeling', {
  fit <- readRDS(test_path('fixtures','legacy-anchored-fit.rds'))
  result <- mfrm_results(fit,include='fit',compute='never')
  root <- tempfile('anchor-migration-');dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  fail <- function(...) stop('Reopening saved results must not recompute')
  local_mocked_bindings(fit_mfrm=fail,mfrm_results=fail,diagnose_mfrm=fail,
    compute_mml_parameter_covariance=fail,.package='mfrmr')
  exported <- export_mfrm_results(result,output_dir=root,include='replay',
    acknowledge_sensitive=TRUE)
  files <- exported$written_files
  script <- files$Path[files$Component=='replay_code']
  expect_length(script,1L)
  fresh <- new.env(parent=asNamespace('mfrmr'))
  expect_no_error(capture.output(withr::with_dir(dirname(script),
    sys.source(basename(script),envir=fresh))))
  expect_identical(fresh$res,result)
  expect_identical(fresh$report,mfrm_report(result))
  expect_null(fit$config$centering_scope)
})

test_that('manifests and methods text retain the fitted centering convention', {
  # Detailed propagation remains in the extended local/CI suite.
  skip_on_cran()
  old <- readRDS(test_path('fixtures','legacy-anchored-fit.rds'))
  current <- readRDS(test_path('fixtures','anchored-facet-information.rds'))$RSM$fit
  local_mocked_bindings(fit_mfrm=function(...) stop('Do not refit reports'),.package='mfrmr')
  expect_error(extract_mfrm_calibration(old),'QUADRATURE_REVIEW_REQUIRED')
  for(fit in list(old,current)) {
    expected <- if(is.null(fit$config$centering_scope)) 'remainder' else 'facet'
    diagnostics <- diagnose_mfrm(fit,residual_pca='none')
    manifest <- build_mfrm_manifest(fit,diagnostics=diagnostics)
    expect_identical(manifest$model_settings$Value[manifest$model_settings$Setting=='centering_scope'],expected)
    text <- gsub('[[:space:]]+',' ',build_apa_reporting_contract(fit,diagnostics)$report_text)
    expect_match(text,if(expected=='remainder') 'Legacy centering' else
      'including anchors and group totals',fixed=TRUE)
    recipe <- build_mfrm_replay_script(fit,diagnostics=diagnostics)
    expect_match(recipe$script,paste0('centering_scope = "',expected,'"'),fixed=TRUE)
  }
})
