test_that("agreement preserves matched context identity and expected coverage", {
  obs <- data.frame(Person = rep(c('a|b', 'a'), each = 2),
                    Task = rep(c('c', 'b|c'), each = 2),
                    Rater = rep(c('r1', 'r2'), 2), Observed = c(0, 0, 1, 0))
  probabilities <- rbind(c(.8,.2),c(.8,.2),c(.2,.8),c(.6,.4))
  res <- list(facets = list(person = data.frame(Person = unique(obs$Person), Estimate = 0)))
  local_mocked_bindings(compute_prob_matrix = function(...) probabilities, .package = 'mfrmr')
  a <- calc_interrater_agreement(obs, c('Person','Task','Rater'), 'Rater', res = res)
  expect_equal(a$pairs$N, 2)
  expect_equal(a$pairs$Exact, .5)
  expect_equal(a$pairs$ExpectedExact, mean(c(.68,.44)))
  single <- calc_interrater_agreement(obs[1:2,], c('Person','Task','Rater'), 'Rater')
  expect_equal(single$summary$Contexts, 1)
  probabilities[3,] <- NA_real_
  a <- calc_interrater_agreement(obs, c('Person','Task','Rater'), 'Rater', res = res)
  expect_equal(a$pairs$ExpectedAvailableContexts, 1)
  expect_equal(a$pairs$ExpectedUnavailableContexts, 1)
  expect_true(is.na(a$pairs$ExpectedExact))
  expect_true(is.na(a$summary$ExpectedAgreements))
  expect_true(is.na(a$summary$AgreementMinusExpected))
  wide <- rater_network_score_wide(obs, c('Person','Task','Rater'), 'Rater')
  expect_equal(nrow(wide$wide), 2)
  expect_equal(rater_network_direction_pairs(wide$wide, wide$raters)$DirectionN, 2)
  repeated <- obs[c(1,1,2,3,4),]
  probabilities <- rbind(c(.8,.2),c(.8,.2),c(.8,.2),c(.2,.8),c(.6,.4))
  a <- calc_interrater_agreement(repeated, c('Person','Task','Rater'), 'Rater', res = res)
  expect_equal(a$summary$RepeatedCells, 1)
  expect_true(is.na(a$pairs$ExpectedExact))
})

test_that("unavailable agreement rules cannot become passing flags", {
  pairs <- data.frame(Rater1 = 'r1', Rater2 = c('r2','r3','r4'), N = c(0,2,2),
    Exact = c(NA,.2,.9), Corr = c(NA,NA,NA), ExpectedExact = NA_real_,
    MeanDiff = NA_real_, MAD = NA_real_)
  local_mocked_bindings(
    mfrm_results_validate_diagnostics_identity = function(...) invisible(NULL),
    calc_interrater_agreement = function(...) list(pairs = pairs, summary = data.frame(Pairs = 3)),
    .package = 'mfrmr')
  fit <- structure(list(config = list(facet_names = 'Rater')), class='mfrm_fit')
  a <- interrater_agreement_table(fit, list(obs = data.frame(x=1)),
    rater_facet='Rater', include_precision=FALSE, top_n=1)
  expect_equal(a$summary$FlaggedPairs, 1)
  expect_equal(a$summary$ClassifiedPairs, 1)
  expect_equal(a$summary$UnclassifiedPairs, 2)
  expect_true(is.na(a$summary$FlaggedShare))
  expect_true(a$pairs$Flag)
  expect_true(is.na(a$pairs$LowCorrFlag))
  a <- interrater_agreement_table(fit, list(obs=data.frame(x=1)), rater_facet='Rater', include_precision=FALSE)
  expect_equal(sum(is.na(a$pairs$Flag)), 2)
  for (type in c('exact','corr','difference')) {
    p <- plot_interrater_agreement(a, plot_type=type, draw=FALSE)
    expect_equal(nrow(p$data$pairs), 3)
    expect_match(p$data$subtitle, 'unassessed', fixed=TRUE)
  }
  path <- tempfile(fileext='.pdf'); grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add=TRUE)
  expect_no_error(plot_interrater_agreement(a, plot_type='corr'))
  expect_no_error(plot_interrater_agreement(a, plot_type='difference'))
  old <- a; old$summary$ClassifiedPairs <- NULL
  expect_error(summary(old), 'Recreate')
  expect_error(plot(old, draw=FALSE), 'Recreate')
})

test_that("network review cannot approve unavailable design checks", {
  raw <- data.frame(Nodes=3, Edges=2, Connected=NA, Components=NA,
                    ArticulationPoints=NA, Bridges=NA)
  expect_identical(network_review_status(raw)$NetworkReviewStatus, 'review')
  raw$Connected <- TRUE; raw$Components <- 1; raw$ScopeComplete <- FALSE
  expect_identical(network_review_status(raw)$NetworkReviewStatus, 'review')
  old <- summary(structure(list(network_summary=raw),class='mfrm_network_review'))
  old$overview$NetworkReviewStatus <- 'ok'
  expect_identical(build_summary_table_bundle(old)$tables$overview$NetworkReviewStatus,'review')
  expect_identical(old$overview$NetworkReviewStatus,'ok')
  raw$Connected <- FALSE
  expect_identical(network_review_status(raw)$NetworkReviewStatus, 'warning')
  expect_false(grepl('design_diagnostic', paste(capture.output(
    print(mfrm_descriptive_display(network_review_overview(raw)))),collapse=' ')))
})

test_that("rater networks retain isolates without inventing directional balance", {
  skip_if_not_installed('igraph')
  obs <- data.frame(Person=rep(paste0('p',1:4),each=2), Rater=rep(c('r1','r2'),4),
                    Observed=rep(c(0,1),4))
  fit <- structure(list(config=list(facet_names='Rater')),class='mfrm_fit')
  local_mocked_bindings(mfrm_results_validate_diagnostics_identity=function(...) invisible(NULL),
    compute_prob_matrix=function(...) NULL,.package='mfrmr')
  a <- rater_network_analysis(fit,list(obs=obs),rater_facet='Rater',mode='agreement')
  expect_equal(a$summary$Edges,0)
  expect_equal(a$summary$Components,2)
  obs$Observed <- rep(c(0,0),4)
  a <- rater_network_analysis(fit,list(obs=obs),rater_facet='Rater',mode='severity_direction')
  expect_equal(a$summary$Edges,0)
  expect_true(all(is.na(a$node_metrics$SeverityIndex)))
  expect_true(all(a$node_metrics$RelativePattern=='insufficient_directional_edges'))
  expect_equal(a$summary$UnavailableDirectionalRaters,2)
  expect_true(is.na(a$summary$MeanSeverityIndex))
  expect_no_error(capture.output(print(summary(a))))
})

test_that("halo reviews keep missing comparisons and withhold dependent-edge tests", {
  skip_if_not_installed('igraph')
  obs <- expand.grid(Person=paste0('p',1:6),Rater=c('r1','r2'),Criterion=c('c1','c2'))
  obs$Observed <- as.numeric(factor(obs$Person))
  obs$Observed[obs$Rater=='r2'] <- 1
  fit <- structure(list(config=list(facet_names=c('Rater','Criterion'))),class='mfrm_fit')
  local_mocked_bindings(mfrm_results_validate_diagnostics_identity=function(...) invisible(NULL),.package='mfrmr')
  a <- rater_halo_network_analysis(fit,list(obs=obs),rater_facet='Rater',criterion_facet='Criterion')
  expect_true(all(is.na(a$summary[c('WelchT','WelchDF','WelchP')])))
  expect_equal(a$summary$AvailablePairs,1)
  expect_equal(a$summary$UnavailablePairs,5)
  expect_false(any(a$halo_summary_by_rater$ReviewStatus=='ok'))
  expect_true(any(a$halo_summary_by_rater$UnavailableIncidentNonHaloPairs > 0))
  expect_match(paste(a$caveats$Message,collapse=' '),'dependent',fixed=TRUE)
  old <- a; old$summary$AvailablePairs <- NULL
  expect_error(summary(old),'Recreate')
  collision <- data.frame(Person='p',Rater=c('a::b','a'),Criterion=c('c','b::c'),Observed=c(1,2))
  wide <- halo_network_wide_scores(collision,'Person','Rater','Criterion')
  expect_equal(nrow(wide$nodes),2)
  expect_equal(length(unique(wide$nodes$Node)),2)
})

test_that("agreement and design networks reject diagnostics from another fit", {
  fit <- make_toy_fit()
  dx <- make_toy_diagnostics(fit)
  other <- fit; other$facets$others$Estimate[1] <- other$facets$others$Estimate[1] + 1
  expect_error(interrater_agreement_table(other,dx),'mismatch')
  expect_error(subset_connectivity_report(other,dx),'mismatch')
  expect_error(plot_rater_agreement_heatmap(other,dx),'mismatch')
  if (requireNamespace('igraph',quietly=TRUE)) {
    expect_error(rater_network_analysis(other,dx),'mismatch')
    expect_error(rater_halo_network_analysis(other,dx),'mismatch')
  }
})

test_that("timing summaries preserve numeric factor values and excluded groups", {
  d <- data.frame(Person=c('p1','p1','p2'),Time=factor(c('10',NA,NA)),Rater=c('r1','r1','r2'))
  a <- response_time_review(d,'Person',facets='Rater',time='Time',rapid_threshold=NULL)
  expect_equal(a$observations$Time,10)
  expect_identical(a$thresholds$Basis,c('quantile_0.05','quantile_0.95'))
  expect_equal(a$person_summary$InputRows,c(2,1))
  expect_equal(a$person_summary$N,c(1,0))
  expect_equal(a$person_summary$ExcludedRows,c(1,1))
  expect_true(is.na(a$person_summary$RapidRate[2]))
  expect_true(is.na(a$person_summary$RapidResponses[2]))
  expect_equal(a$overview$UnassessedPersons,1)
  text <- paste(capture.output(print(summary(a))),collapse=' ')
  expect_false(grepl('quantile_|high_rapid_response_rate',text))
  expect_match(text,'valid timed rows',fixed=TRUE)
  expect_equal(summary(a)$availability$InputRows,c(2,1))
  path <- tempfile(fileext='.pdf'); grDevices::pdf(path)
  on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  expect_no_error(plot(a,type='distribution'))
  expect_no_error(plot(a,type='facet'))
  old <- a; old$overview$UnassessedPersons <- NULL
  expect_error(summary(old),'Recreate')
  expect_error(plot(old,draw=FALSE),'Recreate')
  expect_error(response_time_review(d,'Person',time='Time',rapid_threshold=20,slow_threshold=10),'must not exceed')
  expect_error(response_time_review(d,'Person',time='Time',min_n_flag=1.5),'positive integer')
  expect_error(response_time_review(d,'Person',time='Time',rapid_threshold=c(1,2)),'one finite number')
})
