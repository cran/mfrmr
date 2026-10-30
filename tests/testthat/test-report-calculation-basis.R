report_basis_diagnostics <- function() {
  measures <- data.frame(Facet = c('Item','Person','Rater'), Level = c('I1','P1','R1'),
    Estimate = c(-.2,.4,.1), SE = c(.1,.3,.2), N = c(10,4,8),
    Infit = c(.8,.9,NA), Outfit = c(.7,1.1,NA),
    InfitZSTD = c(-1,0,NA), OutfitZSTD = c(-1,.5,NA),
    SE_Method = c('Observed information (MML)','Posterior SD (EAP)','Observed information (MML)'))
  structure(list(measures = measures, obs = data.frame(PersonMeasureBasis = 'saved EAP'),
    diagnostic_basis = data.frame(DiagnosticPath = 'legacy_residual_fit',
      Basis = 'plugin_residuals_and_eap_tables'),
    reliability = data.frame(Facet = c('Item','Person','Rater'), Levels = c(3L,300L,2L),
      Reliability = c(.99,.63,NA), Separation = c(10,1.3,NA), Strata = c(14,2,NA),
      RealReliability = c(.98,.60,NA), RealSeparation = c(9,1.2,NA), RealStrata = c(12,1.9,NA),
      SupportsFormalInference = FALSE, ReliabilityUse = 'review_before_reporting',
      ModelSummaryNote = c('available','available','Missing standard errors.'),
      RealSummaryNote = c('available','available','Missing standard errors.'))),
    class = 'mfrm_diagnostics')
}

test_that('fit tables and summaries retain values and the saved calculation target', {
  local_mocked_bindings(fit_mfrm = function(...) stop('No fitting'),
    diagnose_mfrm = function(...) stop('No diagnostics'))
  d <- report_basis_diagnostics()
  fm <- fit_measures_table(d, include_person = TRUE, sort_by = 'facet', fit_df_method = 'engine')
  expect_equal(fm$table$Infit, d$measures$Infit)
  expect_equal(fm$table$Outfit, d$measures$Outfit)
  expect_identical(fm$diagnostic_basis, d$diagnostic_basis)
  expect_true(all(fm$table$FitBasis == 'saved EAP'))
  expect_identical(fm$facets_table[['Fit basis']], fm$table$FitBasis)
  expect_identical(fm$facets_table[['Fit interpretation']], fm$table$FitInterpretation)
  expect_match(paste(summary(fm)$notes, collapse=' '), 'expectation-one')
  expect_identical(summary(fm)$diagnostic_basis, d$diagnostic_basis)
  expect_identical(summary(fm)$fit_basis, fm$fit_basis)
  expect_identical(unserialize(serialize(fm,NULL)), fm)
  p <- plot(fm, type='infit_outfit', draw=FALSE)
  expect_true('FitBasis' %in% names(p$data$table))
})

test_that('JML and imported row definitions are not relabelled as EAP residuals', {
  d <- report_basis_diagnostics()
  d$obs$PersonMeasureBasis <- 'JML optimizer'
  jml <- fit_measures_table(d,include_person=TRUE,fit_df_method='engine')
  expect_true(all(jml$table$FitBasis == 'JML optimizer'))
  expect_false(any(grepl('at saved EAP', jml$table$FitInterpretation)))
  d$measures$FitBasis <- c('Source posterior-averaged item fit',
    'Source WLE person fit; estimates remain EAP','Source posterior-averaged item fit')
  class(d) <- c('mfrm_imported_diagnostics','mfrm_diagnostics')
  imported <- fit_measures_table(d,include_person=TRUE,sort_by='facet',fit_df_method='engine')
  expect_identical(imported$table$FitBasis, d$measures$FitBasis)
  expect_equal(imported$table$Infit, d$measures$Infit)
  # Old table rows may be sorted differently from their saved diagnostics.
  old <- imported; old$fit_basis <- NULL
  old$table <- old$table[3:1, ]; old$table$FitBasis <- NULL
  x <- list(components=list(fit_measures=old), diagnostics=d)
  restored <- mfrmr:::mfrm_report_fit_calculation_basis(x)
  expect_identical(restored$FitBasis, rev(d$measures$FitBasis))
})

test_that('old saved tables preserve unrecorded and unavailable information', {
  d <- report_basis_diagnostics()
  fm <- fit_measures_table(d,include_person=TRUE,fit_df_method='engine')
  fm$fit_basis <- fm$diagnostic_basis <- NULL
  fm$table$FitBasis <- fm$table$FitInterpretation <- NULL
  x <- list(components=list(fit_measures=fm), diagnostics=d)
  expect_true(all(mfrmr:::mfrm_report_fit_calculation_basis(x)$FitBasis == 'saved EAP'))
  x$diagnostics <- NULL
  expect_true(all(mfrmr:::mfrm_report_fit_calculation_basis(x)$FitBasis == 'not recorded'))
  expect_equal(x$components$fit_measures$table$Infit, fm$table$Infit)
  expect_match(paste(summary(fm)$notes, collapse=' '), 'not recorded')
  missing <- mfrmr:::mfrm_report_fit_evidence_summary(list(diagnostics=d))
  templates <- mfrmr:::mfrm_report_fit_reporting_templates('technical', missing, data.frame(), data.frame())
  expect_match(templates$Template, 'output-collection state')
})

test_that('precision reporting preserves every facet value and its missing-value reasons', {
  d <- report_basis_diagnostics()
  d$reliability$UncertaintyBasis <- 'User supplied basis column'
  d$reliability$ReportingInterpretation <- 'User supplied interpretation column'
  x <- list(diagnostics=d)
  tab <- mfrmr:::mfrm_report_precision_by_facet(x)
  expect_identical(tab[names(d$reliability)], d$reliability)
  expect_identical(tab$InputUncertaintyBasis[2], 'Posterior SD (EAP)')
  expect_match(tab$ReportInterpretation[2], 'not posterior-variance EAP reliability')
  ev <- mfrmr:::mfrm_report_precision_evidence_summary(x)
  templates <- mfrmr:::mfrm_report_precision_reporting_templates('technical',ev,data.frame(),tab)
  rows <- templates[templates$Topic=='Reliability wording', ]
  expect_identical(rows$Facet, d$reliability$Facet)
  expect_match(rows$Template[2], 'model = 0.63; fit-adjusted = 0.60',fixed=TRUE)
  expect_match(rows$Template[3], 'model = NA',fixed=TRUE)
  expect_match(rows$Caveat[3], 'Missing standard errors',fixed=TRUE)
  expect_false(any(grepl('ranged from',templates$Template)))
  components_only <- list(components=list(diagnostics=d))
  expect_identical(mfrmr:::mfrm_report_precision_evidence_summary(components_only), ev)
  expect_identical(mfrmr:::mfrm_report_precision_by_facet(components_only), tab)
  expect_false(any(tab$SupportsFormalInference))
  expect_match(ev$Boundary, 'not an overall reliability coefficient',fixed=TRUE)
  expect_match(paste(templates$Caveat,collapse=' '),'covariance')
})

test_that('minimal old precision rows do not invent metadata or new reliability values', {
  x <- list(diagnostics=list(reliability=data.frame(Facet='Person',Reliability=.63)))
  tab <- mfrmr:::mfrm_report_precision_by_facet(x)
  expect_identical(tab$Reliability,.63)
  expect_identical(tab$InputUncertaintyBasis,'not recorded')
  expect_false('Levels' %in% names(tab))
  ev <- mfrmr:::mfrm_report_precision_evidence_summary(x)
  out <- mfrmr:::mfrm_report_precision_reporting_templates('technical',ev,data.frame(),tab)
  expect_true(any(grepl('levels: NA',out$Template,fixed=TRUE)))
  expect_true(any(grepl('model = 0.63',out$Template,fixed=TRUE)))
  only_review <- list(components=list(precision_review=structure(list(
    profile=data.frame(PrecisionTier='review', SupportsFormalInference=FALSE)),
    class='mfrm_precision_review')))
  empty <- mfrmr:::mfrm_report_precision_by_facet(only_review)
  expect_equal(nrow(empty),0)
  out <- mfrmr:::mfrm_report_precision_reporting_templates('technical',
    mfrmr:::mfrm_report_precision_evidence_summary(only_review),data.frame(),empty)
  expect_true(any(out$Topic=='Facet precision unavailable'))
  expect_false(any(out$Topic=='Reliability wording'))
})
