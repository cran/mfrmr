test_that("zero fit and centered chi-square retain their mathematical boundaries", {
  expect_equal(zstd_from_mnsq_facets(0,10),-6.5591327,tolerance=1e-7)
  expect_true(is.na(zstd_from_mnsq_facets(-.1,10)))
  tab <- data.frame(Facet='Rater',Estimate=c(-1,0,1),SE=.2,Fixed=FALSE)
  a <- calc_facets_chisq(tab); tab$Estimate <- tab$Estimate+1e8
  b <- calc_facets_chisq(tab)
  expect_equal(a$FixedChiSq,50,tolerance=1e-12)
  expect_equal(a$RandomChiSq,2,tolerance=1e-12)
  expect_equal(a[c('FixedChiSq','FixedProb','RandomChiSq','RandomProb')],
    b[c('FixedChiSq','FixedProb','RandomChiSq','RandomProb')],tolerance=1e-12)
})

manual_boundary_fit <- local({
  f <- NULL
  function() {
    if (is.null(f)) f <<- readRDS(test_path('fixtures','anchored-facet-information.rds'))$RSM$fit
    f
  }
})

# These display tests use one immutable calibration and identical diagnostic
# settings. Build the fixture once per file; modified copies remain local to
# their test. The independent diagnostic/numerical tests are unchanged.
manual_boundary_diagnostics <- local({
  d <- NULL
  function() {
    if (is.null(d)) d <<- diagnose_mfrm(manual_boundary_fit(), residual_pca = "none")
    d
  }
})

test_that("chi-square unavailable flags and summary survive display truncation", {
  f <- manual_boundary_fit(); d <- manual_boundary_diagnostics()
  all <- facets_chisq_table(f,d); one <- facets_chisq_table(f,d,top_n=1)
  expect_equal(all$summary,one$summary)
  expect_equal(nrow(one$table),1L)
  anchor <- all$table[all$table$Facet=='Rater',]
  expect_true(is.na(anchor$FixedFlag))
  expect_true(is.na(anchor$RandomFlag))
  expect_gt(all$summary$FixedUnavailable,0)
})

test_that("weighted display preserves fractions and SE fallback scales once", {
  f <- manual_boundary_fit(); d <- manual_boundary_diagnostics()
  row <- which(d$measures$Facet=='Rater' & !d$measures$Fixed)[1]
  d$measures$ModelSE[row] <- .2; d$measures$RealSE[row] <- NA_real_
  d$measures$Infit[row] <- 1.8
  raw <- calc_facets_report_tbls(f,d,uscale=10)$Rater
  index <- match(d$measures$Level[row],raw$Level)
  expect_equal(unname(raw$RealSE[index]),10*.2*sqrt(1.8),tolerance=1e-12)
  raw$WeightdScore <- .501; raw$WeightdCount <- .192
  formatted <- format_facets_report_gt(raw,facet="Rater",decimals=4)
  expect_equal(formatted$`Weightd Count`,rep(.192,nrow(raw)))
  expect_equal(formatted$`Weightd Score`,rep(.501,nrow(raw)))
})

test_that("outcome-selected category residuals are not unit-baseline fit tests", {
  y <- c(rep(0,90),rep(1,10))
  obs <- data.frame(Observed=y,Expected=.1,PersonMeasure=2-2*y,
    MeasureSum=qlogis(.1),Residual=y-.1,Var=.09,StdSq=(y-.1)^2/.09)
  tab <- calc_category_stats(obs)
  expect_equal(tab$Outfit,c(1/9,9))
  expect_equal(tab$AvgMeasure,rep(qlogis(.1),2))
  expect_true(all(is.na(tab$InfitZSTD) & is.na(tab$OutfitZSTD)))
  expect_true(all(is.na(tab$InfitFlag) & is.na(tab$OutfitFlag) & is.na(tab$ZSTDFlag)))
})

test_that("median boundaries solve cumulative probabilities and missing crossings stay missing", {
  f <- manual_boundary_fit()
  out <- category_structure_report(f,theta_points=1001)
  steps <- f$steps$Estimate
  cumulative <- c(0,cumsum(steps)); k <- seq_along(cumulative)-1
  roots <- vapply(seq_along(steps),function(j) uniroot(function(eta) {
    z <- k*eta-cumulative; p <- exp(z-max(z));sum(p[seq_len(j)])/sum(p)-.5
  },c(-6,6),tol=1e-12)$root,numeric(1))
  expect_equal(out$median_thresholds$MedianThreshold,roots,tolerance=1e-4)
  expect_equal(out$median_thresholds$AndrichThreshold,steps)
  expect_gt(max(abs(roots-steps)),.1)
  old <- out; old$calculation_version <- NULL
  expect_error(summary(old),"Recreate category_structure_report")
  expect_error(plot(old,draw=FALSE),"Recreate category_structure_report")
  narrow <- category_structure_report(f,theta_range=c(5,6))
  expect_true(all(is.na(narrow$mean_halfscore_points$MeanBoundaryTheta)))
  expect_true(all(narrow$mean_halfscore_points$CrossingStatus=='outside_theta_range'))
  expect_true(all(is.na(narrow$median_thresholds$MedianThreshold)))
})

test_that("bias-count display never relabels logits as score points", {
  f <- manual_boundary_fit();d <- manual_boundary_diagnostics()
  b <- estimate_bias(f,d,facet_a='Rater',facet_b='Criterion')
  x <- bias_count_table(b,branch='facets')$table
  idx <- match(paste(x$Facet1,x$Facet2),paste(b$table$FacetA_Level,b$table$FacetB_Level))
  expect_equal(x$`Obs-Exp Average`,b$table$`Obs-Exp Average`[idx])
  expect_equal(x$BiasSize,b$table$`Bias Size`[idx])
  expect_equal(x$AfterBiasObsExpAverage,b$table$AfterBiasObsExpAverage[idx])
})

test_that("actual fixed-field suffixes parse and ambiguous layouts fail explicitly", {
  file <- test_path('fixtures','facets-450-score.2.txt')
  x <- read_facets_fit_table(file,facet='Rater')
  expect_equal(x$Level,as.character(1:3))
  ref <- read.csv(test_path('fixtures','facets-450-measures.csv'))
  ref <- ref[ref$Mode=='yes' & ref$Facet=='Rater',]
  expect_equal(x$Estimate,ref$Estimate,tolerance=1e-12)
  expect_equal(x$Infit,ref$Infit,tolerance=1e-12)
  lines <- readLines(file)
  path <- tempfile(fileext='.txt');on.exit(unlink(path),add=TRUE)
  writeLines(c(lines[1],paste0('   Element',lines[-1])),path)
  expect_error(read_facets_fit_table(path,facet='Rater',format='scorefile'),'Unsupported fixed-field')
  writeLines(lines[-2],path)
  expect_error(read_facets_fit_table(path,facet='Rater',format='scorefile'),'requires validated')
})
