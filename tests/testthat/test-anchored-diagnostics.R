anchored_diagnostic_fixture <- local({
  cache <- list()
  function(case) {
    if (is.null(cache[[case]])) {
      ref <- readRDS(test_path('fixtures','anchored-facet-information.rds'))[[case]]
      cache[[case]] <<- list(fit=ref$fit,reference=ref,
        diag=diagnose_mfrm(ref$fit,residual_pca='none'))
    }
    cache[[case]]
  }
})

# Independent adjacent-category probabilities from the public coefficient tables.
anchored_diagnostic_probabilities <- function(fit, rows, shift=0) {
  eta <- fit$facets$person$Estimate[match(rows$Person,fit$facets$person$Person)]
  for (facet in fit$config$facet_names) {
    tab <- fit$facets$others[fit$facets$others$Facet==facet,]
    eta <- eta + fit$config$facet_signs[[facet]] * tab$Estimate[match(rows[[facet]],tab$Level)]
  }
  K <- fit$config$n_cat-1L
  t(vapply(seq_len(nrow(rows)),function(i) {
    steps <- fit$steps
    if (fit$config$model=='PCM') steps <- steps[steps$StepFacet==rows[[fit$config$step_facet]][i],]
    logits <- (0:K)*(eta[i]+shift)-c(0,cumsum(steps$Estimate))
    weights <- exp(logits-max(logits)); weights/sum(weights)
  },numeric(K+1L)))
}

for (case in c('RSM','PCM','RSM_group')) test_that(paste(case,'anchored fit and bias use the stored coordinates'), {
  x <- anchored_diagnostic_fixture(case); fit <- x$fit; dx <- x$diag
  local_mocked_bindings(fit_mfrm=function(...) stop('No refitting'),.package='mfrmr')
  obs <- dx$obs; K <- fit$config$n_cat-1L
  p <- anchored_diagnostic_probabilities(fit,obs)
  expected <- drop(p %*% (0:K))
  centered <- sweep(matrix(rep(0:K,each=nrow(p)),nrow(p)),1,expected)
  variance <- rowSums(p*centered^2)
  fourth <- rowSums(p*centered^4)
  expect_equal(obs$Expected,expected+fit$prep$rating_min,tolerance=1e-12)
  expect_equal(obs$Var,variance,tolerance=1e-12)
  expect_equal(obs$FourthCentralMoment,fourth,tolerance=1e-12)
  residual <- obs$Observed-obs$Expected
  for (rater in unique(obs$Rater)) {
    index <- obs$Rater==rater
    row <- dx$measures[dx$measures$Facet=='Rater' & dx$measures$Level==rater,]
    expect_equal(row$Infit,sum(residual[index]^2)/sum(variance[index]),tolerance=1e-12)
    expect_equal(row$Outfit,mean(residual[index]^2/variance[index]),tolerance=1e-12)
  }
  rater <- dx$measures[dx$measures$Facet=='Rater',]
  fixed <- if (case=='RSM_group') c(TRUE,FALSE,FALSE,TRUE) else c(TRUE,FALSE,FALSE,FALSE)
  expect_identical(rater$Fixed,fixed)
  expect_true(all(is.na(rater$ModelSE[fixed])))
  J <- x$reference$rater_jacobian
  variance <- diag(J %*% x$reference$covariance %*% t(J))
  expect_equal(rater$ModelSE[!fixed],unname(sqrt(variance[!fixed])),tolerance=1e-6)
  for (facet in unique(dx$reliability$Facet)) {
    measure <- dx$measures[dx$measures$Facet==facet,]
    rel <- dx$reliability[dx$reliability$Facet==facet,]
    if (any(measure$Fixed)) {
      expect_true(all(is.na(unlist(rel[c('Reliability','Separation','Strata')]))))
      expect_match(rel$ModelSummaryNote,'fixed by constraints')
    } else {
      v <- var(measure$Estimate); e <- mean(measure$ModelSE^2)
      expect_equal(rel$Reliability,max(v-e,0)/v,tolerance=1e-12)
      expect_equal(rel$Strata,(4*sqrt(max(v-e,0)/e)+1)/3,tolerance=1e-12)
    }
  }
  bias <- estimate_bias(fit,dx,facet_a='Rater',facet_b='Criterion')
  for (i in seq_len(nrow(bias$table))) {
    b <- bias$table[i,]
    cell <- obs[obs$Rater==b$Facet1_Level & obs$Criterion==b$Facet2_Level,]
    y <- cell$Observed-fit$prep$rating_min
    score <- function(shift) sum(y-drop(anchored_diagnostic_probabilities(fit,cell,shift) %*% (0:K)))
    root <- uniroot(score,c(-10,10),tol=1e-10)$root
    # The bounded likelihood optimizer uses its default scalar tolerance.
    expect_lt(abs(b[['Bias Size']]-root),5e-5)
    probs <- anchored_diagnostic_probabilities(fit,cell,b[['Bias Size']])
    means <- drop(probs %*% (0:K))
    centered <- sweep(matrix(rep(0:K,each=nrow(probs)),nrow(probs)),1,means)
    expect_equal(b[['S.E.']],1/sqrt(sum(probs*centered^2)),tolerance=1e-12)
    baseline <- drop(anchored_diagnostic_probabilities(fit,cell) %*% (0:K))
    expect_lt(abs(b[['Obs-Exp Average']]-mean(y-baseline)),1e-12)
    expect_lt(abs(b[['AfterBiasObsExpAverage']]-mean(y-means)),1e-12)
  }
  expect_true(all(bias$table$ReportingUse=='screening_only'))
  expect_false(any(bias$table$SupportsFormalInference))
  expect_true(all(is.finite(bias$table[['S.E.']][bias$table$Facet1_Level=='R01'])))
  residual_dff <- analyze_dff(fit,facet='Rater',group='Group',data=load_mfrmr_data('example_core'))
  expect_true(all(is.na(residual_dff$dif_table$SE)))
  expect_true(all(residual_dff$dif_table$ReportingUse=='descriptive'))
  saved <- unserialize(serialize(dx,NULL))
  expect_identical(saved$measures,dx$measures)
  expect_identical(summary(saved)$reliability,summary(dx)$reliability)
})

test_that('refit DFF preserves target anchors and distinguishes fixed values from comparisons', {
  x <- anchored_diagnostic_fixture('RSM'); fit <- x$fit; calls <- list()
  # Test orchestration with a retained fit; actual subgroup estimation is
  # independently checked in the repository validation record.
  local_mocked_bindings(fit_mfrm=function(...) {calls[[length(calls)+1L]] <<- list(...); fit},
    diagnose_mfrm=function(...) x$diag,.package='mfrmr')
  result <- analyze_dff(fit,x$diag,facet='Rater',group='Group',
    data=load_mfrmr_data('example_core'),method='refit',min_obs=1)
  expect_length(calls,2L)
  for (call in calls) {
    a <- call$anchors
    expect_equal(a$Anchor[a$Facet=='Rater' & a$Level=='R01'],.7,tolerance=0)
    expect_equal(sum(a$Facet=='Criterion'),4L)
    expect_identical(call$centering_scope,'facet')
  }
  row <- result$dif_table[result$dif_table$Level=='R01',]
  expect_true(row$Fixed1 && row$Fixed2)
  expect_equal(row$Estimate1,.7); expect_equal(row$Estimate2,.7)
  expect_true(all(is.na(unlist(row[c('Contrast','SE','p_value','p_adjusted')]))))
  expect_identical(row$ReportingUse,'not_applicable_fixed_target')
  expect_identical(row$Classification,'Fixed target; no comparison')
  expect_false(row$ConditionalRefitScreenEligible)
  expect_true(all(is.finite(result$dif_table$Contrast[result$dif_table$Level!='R01'])))
  expect_equal(sum(result$summary$Count),nrow(result$dif_table))
  expect_identical(summary(result)$dif_table,result$dif_table)
  report <- dif_report(result)
  expect_equal(unname(report$counts['Fixed_target']),1L)
  expect_match(report$narrative,'equal anchors do not establish invariance',fixed=TRUE)
  graph <- plot(result,type='summary',draw=FALSE)
  expect_match(paste(graph$data$notes$Text,collapse=' '),'fixed-target row(s)',fixed=TRUE)
  expect_false(any(grepl('R01',graph$data$data$Pair,fixed=TRUE)))
  group <- anchored_diagnostic_fixture('RSM_group')
  expect_error(analyze_dff(group$fit,group$diag,facet='Rater',group='Group',
    data=load_mfrmr_data('example_core'),method='refit'),'group-anchor')
  expect_length(calls,2L)
})

test_that('linking retains baseline JML Person anchors as well as target anchors', {
  x <- anchored_diagnostic_fixture('RSM'); fit <- x$fit
  # Isolate anchor-table construction; this is not a fitted JML example.
  fit$config$method <- 'JML'
  fit$config$theta_spec <- build_facet_constraint(c('P1','P2'),anchors=c(P1=.2),centered=FALSE)
  link <- build_dff_linking_setup(fit,'Rater',fit$config$facet_names)
  expect_equal(link$anchor_tbl$Anchor[link$anchor_tbl$Facet=='Person'],.2)
  expect_equal(link$anchor_tbl$Anchor[link$anchor_tbl$Facet=='Rater'],.7)
  expect_identical(link$linking_facets,'Criterion')
})

test_that('fixed-zero SE conventions cannot produce a spurious separation reliability', {
  measures <- data.frame(Facet='Rater',Level=c('a','b','c'),Estimate=c(-1,0,1),
    SE=c(0,.1,.1),ModelSE=c(0,.1,.1),RealSE=c(0,.1,.1),Fixed=c(TRUE,FALSE,FALSE),
    PrecisionTier='model_based',Converged=TRUE)
  rel <- calc_reliability(measures)
  expect_true(all(is.na(unlist(rel[c('Reliability','Separation','Strata','TrueSD','RMSE')]))))
  precision <- build_facet_precision_summary(measures)
  expect_true(all(is.na(precision$Reliability)))
  expect_match(rel$ModelSummaryNote,'fixed by constraints')
  expect_equal(rel$ObservedVariance,var(measures$Estimate))
})
