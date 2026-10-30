# FACETS 4.5.0, executed on Windows on 2026-10-23. Synthetic binary ratings,
# unit weights, Person positive, centered negative Rater/Task; ordinary JML.
# The reference reports print measures/SEs to six decimals, MnSq/PTMEA/ZSTD to
# two. Bias columns are narrower. Tolerances include those print resolutions
# and the retained FACETS stopping criterion (score/logit change 1e-6).
facets450_fit <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      x <- read.csv(test_path('fixtures','facets-450-ratings.csv'))
      value <<- fit_mfrm(x, 'Person', c('Rater','Task'), 'Score',
        method='JML', rating_min=0, rating_max=1, min_obs_per_element=1,
        min_obs_per_category=1, maxit=1500, reltol=1e-12)
    }
    value
  }
})

test_that('matched JML statistics reproduce executed FACETS YES and NO reports', {
  f <- facets450_fit()
  ref <- read.csv(test_path('fixtures','facets-450-measures.csv'))
  for (mode in c('yes','no')) {
    d <- diagnose_mfrm(f, whexact=mode=='no', diagnostic_mode='legacy')
    x <- ref[ref$Mode==mode,]
    y <- d$measures[match(paste(x$Facet,x$Level),paste(d$measures$Facet,d$measures$Level)),]
    expect_lt(max(abs(y$Estimate-x$Estimate)), 1e-5)
    expect_lt(max(abs(y$SE-x$SE)), 2e-6)
    for (nm in c('Infit','Outfit','PTMEA','InfitZSTD','OutfitZSTD'))
      expect_lt(max(abs(y[[nm]]-x[[nm]])), .0051)
    expect_identical(d$calculation_version, 2L)
  }
  expect_equal(mfrmr:::zstd_from_mnsq(2,10,TRUE),1.7436034,tolerance=1e-7)
  expect_equal(mfrmr:::zstd_from_mnsq(2,10,'legacy_linear'),sqrt(5))
  expect_error(mfrmr:::zstd_from_mnsq(2,10,'unknown'),'whexact')
})

test_that('bias reports retain baseline departure and test zero without centering', {
  f <- facets450_fit(); d <- diagnose_mfrm(f,diagnostic_mode='legacy')
  b <- estimate_bias(f,d,facet_a='Rater',facet_b='Task')
  ref <- read.csv(test_path('fixtures','facets-450-bias.csv'))
  y <- b$table[match(paste(ref$Rater,ref$Task),paste(b$table$FacetA_Level,b$table$FacetB_Level)),]
  expect_lt(max(abs(y$`Expctd Score`-ref$Expected)),.0051)
  expect_lt(max(abs(y$`Obs-Exp Average`-ref$ObsExpAverage)),.0051)
  expect_lt(max(abs(y$`Bias Size`-ref$Bias)),.0011)
  expect_lt(max(abs(y$`S.E.`-ref$SE)),.0011)
  expect_equal(y$`d.f.`,ref$df)
  expect_lt(max(abs(y$`Prob.`-ref$p)),.00011)
  expect_lt(max(abs(y$AfterBiasObsExpAverage)),1e-5)
  expect_equal(b$chi_sq$FixedChiSq,sum((y$`Bias Size`/y$`S.E.`)^2))
  expect_equal(b$chi_sq$FixedDF,12L)
  expect_equal(round(b$chi_sq$FixedChiSq,1),8.4)
  old <- d; old$calculation_version <- NULL
  expect_error(estimate_bias(f,old,facet_a='Rater',facet_b='Task'),'Recompute diagnose_mfrm')
  old <- b; old$calculation_version <- NULL
  expect_error(mfrmr:::mfrm_validate_bias_identity(f,old,'test'),'Recreate estimate_bias')
})

test_that('displacement follows the fitted coordinate and chi-square df follows eligible levels', {
  f <- facets450_fit(); d <- diagnose_mfrm(f,diagnostic_mode='legacy')
  obs <- d$obs
  obs$ObservedScoreDerivative <- 1
  obs$ScoreInformation <- 1
  x <- mfrmr:::calc_displacement_table(obs,f,d$measures)
  expect_equal(x$Displacement[x$Facet=='Person'],rep(1,8))
  expect_equal(x$Displacement[x$Facet=='Rater'],rep(-1,3))
  expect_equal(x$ReleasedEstimate,x$Estimate+x$Displacement)
  chi <- mfrmr:::calc_facets_chisq(data.frame(Facet='Rater',
    Estimate=c(-1,1,Inf),SE=c(1,1,NA_real_),Fixed=FALSE))
  expect_equal(chi$EligibleLevels,2)
  expect_equal(chi$ExcludedLevels,1)
  expect_equal(chi$FixedChiSq,2)
  expect_equal(chi$FixedDF,1)
  expect_equal(chi$FixedProb,pchisq(2,1,lower.tail=FALSE))
})

test_that('unexpected defaults and tied modes do not manufacture asymmetry', {
  obs <- data.frame(Person='p',score_k=0:4,Observed=0:4,Expected=2,
    Residual=(-2):2,StdResidual=((-2):2)/sqrt(2))
  probs <- matrix(.2,5,5)
  fun <- mfrmr:::calc_unexpected_response_table
  expect_equal(nrow(fun(obs,probs,character(),0)),0L)
  all <- fun(obs,probs,character(),0,abs_z_min=2,rule='either')
  expect_equal(nrow(all),5L)
  expect_true(all(is.na(all$MostLikely)))
  expect_equal(all$MostLikelyCount,rep(5,5))
  expect_equal(all$CategoryGap,rep(0,5))
  expect_equal(all$Severity[all$Observed==0],all$Severity[all$Observed==4])
})

test_that('malformed nonempty constraints cannot silently change the fitted model', {
  x <- read.csv(test_path('fixtures','facets-450-ratings.csv'))
  bad <- data.frame(facet='Rater',level=1,estimate=.7)
  for (policy in c('warn','silent','error')) {
    expect_error(fit_mfrm(x,'Person',c('Rater','Task'),'Score',method='JML',
      anchors=bad,anchor_policy=policy), 'Malformed anchor schema')
  }
})

test_that('PCM fair scores average administered scale expectations and preserve anchors', {
  x <- read.csv(test_path('fixtures','facets-450-ratings.csv'))
  withr::local_seed(241024)
  x$Score <- sample(0:2,nrow(x),replace=TRUE)
  f <- fit_mfrm(x,'Person',c('Rater','Task'),'Score',method='JML',model='PCM',
    step_facet='Task',rating_min=0,rating_max=2,min_obs_per_element=1,
    min_obs_per_category=1,maxit=1500,reltol=1e-11,
    anchors=data.frame(facet='Rater',level=1,anchor=.7))
  d <- diagnose_mfrm(f,diagnostic_mode='legacy')
  tables <- mfrmr:::calc_facets_report_tbls(f,d)
  pars <- mfrmr:::expand_params(f$opt$par,mfrmr:::build_param_sizes(f$config),f$config)
  person <- tables$Person
  eta <- person$PrimaryMeasure + sum(vapply(f$config$facet_names,function(nm)
    f$config$facet_signs[[nm]]*mean(pars$facets[[nm]]),numeric(1)))
  direct <- vapply(eta,function(e) mean(apply(pars$steps_mat,1,function(steps) {
    num <- exp((0:2)*e-c(0,cumsum(steps)))
    sum((0:2)*num)/sum(num)
  })), numeric(1))
  expect_equal(person$FairM,unname(direct),tolerance=1e-10)
  # Controlled display input: an extreme response pattern must not replace a fixed anchor.
  d$obs$Observed[d$obs$Rater=='1'] <- 0
  displayed <- mfrmr:::calc_facets_report_tbls(f,d,xtreme=.3)$Rater
  fixed <- displayed[displayed$Level=='1',]
  expect_equal(unname(fixed$Measure),.7)
  expect_equal(fixed$PrimaryMeasure,.7)
  expect_equal(fixed$ExtremeAdjustment,0)
  f$facets$person$Estimate[1] <- Inf
  finite <- mfrmr:::calc_facets_report_tbls(f,d)$Rater
  expect_true(all(is.finite(finite$FairM)))
  expect_equal(finite$FairMExcludedPersons,rep(1L,nrow(finite)))
})

# Removing four responses from one cell distinguishes the modified df from
# classical Welch; the FACETS tail areas also distinguish truncation from
# using continuous df. No external executable is needed to run this test.
test_that('unequal-count pairwise screens match executed FACETS modified integer df', {
  x <- read.csv(test_path('fixtures','facets-450-unequal-ratings.csv'))
  f <- fit_mfrm(x,'Person',c('Rater','Task'),'Score',method='JML',
    rating_min=0,rating_max=1,min_obs_per_element=1,min_obs_per_category=1,
    maxit=1500,reltol=1e-12)
  d <- diagnose_mfrm(f,diagnostic_mode='legacy',residual_pca='none')
  b <- estimate_bias(f,d,facet_a='Rater',facet_b='Task')
  report <- bias_pairwise_report(b,target_facet='Task',context_facet='Rater')
  p <- report$table
  ref <- read.csv(test_path('fixtures','facets-450-pairwise.csv'))
  key <- function(z) paste(z$Target,z$Context1,z$Context2)
  p <- p[match(key(ref),key(p)),]
  expect_lt(max(abs(p$Contrast-ref$Contrast)), 1.1e-4)
  expect_lt(max(abs(p$SE-ref$SE)), 1.1e-5)
  expect_equal(p$`d.f.`,ref$df)
  expect_lt(max(abs(p$`Prob.`-ref$p)), .000051)
  expect_true(all(p$CalculationVersion==3L))
  expect_true(all(p$DFMethod=='facets'))
  fixed <- build_fixed_reports(b,target_facet='Task')$pairwise_table
  expect_equal(unname(fixed$`Prob.`[match(key(p),key(fixed))]),p$`Prob.`)
  old <- bias_pairwise_report(b,target_facet='Task',context_facet='Rater',df_method='classical')$table
  row <- old$Target=='1' & old$Context1=='1' & old$Context2=='2'
  expect_equal(old$`d.f.`[row],6.878172,tolerance=1e-6)
  expect_gt(abs(old$`Prob.`[row]-p$`Prob.`[key(p)=='1 1 2']),.0005)
  legacy <- report; legacy$table$CalculationVersion <- 2L
  expect_error(summary(legacy),'Recreate build_fixed_reports')
})


test_that('table and plot defaults preserve FACETS df and supplied WHEXACT choices', {
  f <- facets450_fit()
  d <- diagnose_mfrm(f, diagnostic_mode='legacy')
  key <- function(x) paste(x$Facet,x$Level)
  tab <- fit_measures_table(f,include_person=TRUE)$table
  expected <- d$measures[match(key(tab),key(d$measures)),]
  expect_equal(tab$InfitZSTD,expected$InfitZSTD)
  expect_true(all(tab$FitDfMethod=='facets_wright_masters'))
  p <- plot(f,type='fit_pathway',fit_scale='zstd',include_person=TRUE,
    top_n_person=Inf,draw=FALSE)
  expected <- d$measures[match(key(p$data$table),key(d$measures)),]
  expect_identical(p$data$zstd_method,'facets')
  expect_equal(p$data$table$FitValue,expected$InfitZSTD)

  # Engine-only results lack companion columns. Adding them must retain the
  # chosen cube-root variant, rather than silently switching NO back to YES.
  old <- diagnose_mfrm(f,whexact=TRUE,fit_df_method='engine',diagnostic_mode='legacy')
  ref <- read.csv(test_path('fixtures','facets-450-measures.csv'))
  ref <- ref[ref$Mode=='no',]
  tab <- fit_measures_table(f,diagnostics=old,include_person=TRUE)$table
  expected <- ref[match(key(tab),key(ref)),]
  expect_lt(max(abs(tab$InfitZSTD-expected$InfitZSTD)),.0051)
  p <- plot(f,type='fit_pathway',fit_scale='zstd',diagnostics=old,
    include_person=TRUE,top_n_person=Inf,draw=FALSE)
  expected <- ref[match(key(p$data$table),key(ref)),]
  expect_lt(max(abs(p$data$table$FitValue-expected$InfitZSTD)),.0051)
})


test_that("unexpected plots retain the table rule and permit a completely evaluated empty screen", {
  f <- facets450_fit()
  d <- diagnose_mfrm(f, residual_pca="none")
  tab <- unexpected_response_table(f, d)
  p <- plot(f, type="unexpected", diagnostics=d, draw=FALSE)
  expect_identical(p$data$table, as.data.frame(tab$table))
  expect_identical(p$data$thresholds, tab$thresholds)
  expect_identical(p$data$thresholds$rule, "residual")
  expect_false(any(p$data$reference_lines$axis == "h"))
  saved <- unexpected_response_table(f, d, abs_z_min=.5, rule="both", prob_max=.8)
  p <- plot(saved, draw=FALSE)
  expect_identical(p$data$thresholds, saved$thresholds)
  expect_equal(p$data$reference_lines$value, c(-.5, .5, -log10(.8)))
  expect_match(p$data$subtitle, "both")
  empty <- unexpected_response_table(f, d, abs_z_min=100)
  expect_equal(empty$summary$UnexpectedN, 0L)
  file <- tempfile(fileext="pdf")
  grDevices::pdf(file)
  tryCatch(p <- plot(empty, draw=TRUE), finally=grDevices::dev.off())
  unlink(file)
  expect_equal(nrow(p$data$table), 0L)
  expect_equal(p$data$summary$UnavailableObservations, 0L)
  expect_match(p$data$notes, "No unexpected responses")
})
