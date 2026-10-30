gpcm_inference_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      data <- simulate_mfrm_data(n_person=80, n_rater=3, n_criterion=3,
        score_levels=3, model='GPCM', step_facet='Criterion', slope_facet='Criterion',
        slopes=c(C01=.8,C02=1,C03=1.25), seed=924090)
      persons <- data.frame(Person=unique(data$Person), X=rep(c(-1,1),40))
      fit <- fit_mfrm(data,'Person',c('Rater','Criterion'),'Score',model='GPCM',
        method='MML',step_facet='Criterion',slope_facet='Criterion',
        population_formula=~X,person_data=persons,person_id='Person',
        quad_points=31,maxit=400,reltol=1e-10)
      cached <<- list(fit=fit,data=data)
    }
    local_cran_fixture_information(cached$fit, "gpcm-inference", .env = parent.frame())
    cached
  }
})

test_that('verified weak-information cautions reach intervals, curves and reports', {
  # Extended local/CI regression; see CONTRIBUTING.md.
  skip_on_cran()
  f <- gpcm_inference_fixture()$fit
  info <- compute_mml_parameter_covariance(f)
  caution <- 'Ill-conditioned information: inspect interval width and boundary proximity.'
  info$solution_information$inverse_review <- data.frame(Verified=TRUE,
    RelativeChange=1e-5,InverseResidual=1e-8,CurvatureScaledGradient=1e-6,Detail=caution)
  local_mocked_bindings(compute_mml_parameter_covariance=function(...) info,.package='mfrmr')
  expect_warning(relative <- confint(f, parm = "slopes"),'Ill-conditioned')
  expect_warning(ci <- confint(f, parm = "slopes",scale='standardized'),'Ill-conditioned')
  expect_true(all(attr(ci,'diagnostics')$CIEligible))
  expect_identical(attr(ci,'cautions'),caution)
  expect_identical(attr(relative,'cautions'),caution)
  expect_match(attr(relative,'diagnostics')$InferenceReview[1],'Ill-conditioned')
  expect_output(print(ci),'Caution: Ill-conditioned')
  expect_match(apa_table(ci)$table$InferenceReview[1],'Ill-conditioned')
  expect_match(paste(plot_data(plot(ci,draw=FALSE), "notes")$Text, collapse = "\n"),'Ill-conditioned')
  expect_null(plot(ci,draw=FALSE,subtitle=NULL)$labels$subtitle)
  grid <- data.frame(Theta=c(-1,0,1),Rater='R01',Criterion='C01')
  expect_warning(curves <- mfrm_curve_intervals(f,grid),'Ill-conditioned')
  expect_identical(curves$cautions,caution)
  expect_match(curves$table$InferenceReview[1],'Ill-conditioned')
  expect_output(print(curves),'Caution: Ill-conditioned')
  expect_match(paste(plot_data(plot(curves,draw=FALSE), "notes")$Text, collapse = "\n"),'Ill-conditioned')
  expect_null(plot(curves,draw=FALSE,subtitle=NULL)$labels$subtitle)
  saved <- mfrm_results(f,intervals=ci,include='fit',compute='never')
  expect_match(mfrm_report(saved)$tables$gpcm_inference_intervals$InferenceReview[1],'Ill-conditioned')
  # The same explicit caution is retained in bootstrap diagnostics; logging
  # must not invoke another refinement or information calculation.
  inference <- list(information=info,check=mfrm_gpcm_slope_inference_check(f,info))
  check <- mfrm_gpcm_bootstrap_check(f,inference)
  expect_true(check$eligible)
  expect_identical(check$caution,caution)
  local_mocked_bindings(compute_mml_parameter_covariance=function(...) stop('recomputed'),.package='mfrmr')
  tab <- mfrm_gpcm_bootstrap_fit_checks(f,1L,'alternative',inference,check)
  expect_true(tab$InformationRefinementVerified)
  expect_equal(tab$InformationRelativeChange,1e-5)
  expect_equal(tab$InformationInverseResidual,1e-8)
  expect_equal(tab$InformationScaledGradient,1e-6)
  expect_identical(tab$BootstrapCaution,caution)
})

test_that('slope transformations retain full covariance, scale uncertainty and multiplicity', {
  f <- gpcm_inference_fixture()$fit
  cv <- mfrmr:::compute_mml_parameter_covariance(f)
  local_mocked_bindings(compute_mml_parameter_covariance = function(res, ...) {
    stopifnot(identical(res, f))
    cv
  }, .package = 'mfrmr')
  at <- cv$param_slices$log_slopes; sigma <- cv$param_slices$log_sigma2
  # Independent explicit coordinate matrix, including the scale cross terms.
  j <- matrix(0,3,length(f$opt$par)); j[,at] <- rbind(c(1,0),c(0,1),c(-1,-1))
  j[,sigma] <- .5
  expected <- j %*% cv$cov %*% t(j)
  ci <- confint(f, parm = "slopes",scale='standardized')
  expect_true(all(attr(ci,'diagnostics')$CIEligible))
  expect_equal(unname(attr(ci,'covariance')),unname(expected),tolerance=1e-10)
  estimate <- f$slopes$Estimate * sqrt(f$population$sigma2)
  expect_equal(attr(ci,'diagnostics')$Estimate,estimate)
  expect_equal(as.numeric(ci[,1]),estimate*exp(-qnorm(.975)*sqrt(diag(expected))),tolerance=1e-9)
  cmat <- rbind('A/B'=c(1,-1,0),'A/C'=c(1,0,-1))
  colnames(cmat) <- f$slopes$SlopeFacet
  ratio <- confint(f, parm = "slopes",contrasts=cmat)
  standardized <- confint(f, parm = "slopes",contrasts=cmat,scale='standardized')
  expect_equal(as.vector(ratio),as.vector(standardized),tolerance=1e-10)
  reordered <- confint(f, parm = "slopes",contrasts=cmat[,3:1,drop=FALSE])
  expect_equal(ratio,reordered)
  difference <- confint(f, parm = "slopes",contrasts=cmat,contrast_scale='difference')
  expect_equal(attr(difference,'diagnostics')$Estimate,unname(drop(cmat %*% f$slopes$Estimate)))
  simultaneous <- confint(f, parm = "slopes",contrasts=cmat,simultaneous='bonferroni')
  expect_true(all(simultaneous[,1] <= ratio[,1]))
  expect_true(all(simultaneous[,2] >= ratio[,2]))
  expect_equal(attr(simultaneous,'diagnostics')$PValue,pmin(1,2*attr(ratio,'diagnostics')$PValue))
  bad <- cmat; bad[1,1] <- 2
  expect_error(confint(f, parm = "slopes",contrasts=bad),'summing to zero')
  expect_error(confint(f, parm = "slopes",adjust=TRUE),'sandwich')
})

test_that('extreme valid confidence levels preserve slope and curve tail probabilities', {
  # Reuse the eligible live population-GPCM fixture in CRAN tests.
  # The complete suite keeps the separate conditional-scoring fixture.
  f <- if (isTRUE(as.logical(Sys.getenv('NOT_CRAN', 'false'))))
    readRDS(test_path('fixtures', 'mfrm-conditional-scoring-gpcm.rds'))$fit else
    gpcm_inference_fixture()$fit
  # Reuse fresh derivatives for this unchanged fit, retaining every target/tail check.
  information <- compute_mml_parameter_covariance(f)
  scores <- mfrm_mml_person_scores_numeric(f)
  local_mocked_bindings(compute_mml_parameter_covariance = function(res, ...) {
    stopifnot(identical(res, f))
    information
  }, mfrm_mml_person_scores_numeric = function(fit, relative_step = 1e-5) {
    stopifnot(identical(fit, f), identical(relative_step, 1e-5))
    scores
  }, .package = 'mfrmr')
  level <- 1 - .Machine$double.eps / 2
  grid <- unique(f$prep$data[f$config$facet_names])[1:2, , drop = FALSE]
  grid$Theta <- 0
  # Independent inversion of the normal survival probability, including a
  # Bonferroni tail too small to represent its complementary CDF below one.
  critical <- function(n) uniroot(function(z)
    pnorm(z, lower.tail = FALSE, log.p = TRUE) - log((1-level)/(2*n)),
    c(0, 20), tol = 1e-12)$root
  for (method in c('model', 'sandwich')) for (simultaneous in if (isTRUE(as.logical(Sys.getenv('NOT_CRAN', 'false')))) c('none', 'bonferroni') else if (method == 'model') 'none' else 'bonferroni') {
    for (scale in c('relative', 'standardized')) {
      ci <- confint(f, parm = "slopes", level = level, scale = scale, method = method,
        simultaneous = simultaneous)
      tab <- attr(ci, 'diagnostics')
      expect_true(all(tab$CIEligible))
      z <- critical(if (simultaneous == 'none') 1 else nrow(tab))
      expect_equal(log(as.numeric(ci[, 1])), log(tab$Estimate)-z*tab$LogSE, tolerance = 1e-10)
      expect_equal(log(as.numeric(ci[, 2])), log(tab$Estimate)+z*tab$LogSE, tolerance = 1e-10)
    }
    for (type in c('probability', 'information')) {
      ci <- mfrm_curve_intervals(f, grid, type = type, level = level,
        method = method, simultaneous = simultaneous)
      tab <- ci$table
      expect_true(all(tab$CIEligible))
      z <- critical(if (simultaneous == 'none') 1 else nrow(tab))
      link <- if (type == 'probability') qlogis else log
      link_se <- tab$SE / if (type == 'probability') tab$Estimate*(1-tab$Estimate) else tab$Estimate
      expect_equal(link(tab$Lower), link(tab$Estimate)-z*link_se, tolerance = 1e-8)
      expect_equal(link(tab$Upper), link(tab$Estimate)+z*link_se, tolerance = 1e-8)
    }
  }
})

test_that('GPCM sandwich aggregates Person vectors and preserves cluster failure reasons', {
  f <- gpcm_inference_fixture()$fit
  inf <- mfrmr:::mfrm_gpcm_inference(f,method='sandwich')
  expect_true(inf$check$eligible)
  expected <- inf$information$cov %*% crossprod(inf$person_scores) %*% inf$information$cov
  expect_equal(inf$covariance,expected,tolerance=1e-10,ignore_attr=TRUE)
  ids <- rownames(inf$person_scores)
  clusters <- data.frame(Person=ids,Cluster=rep(c('school1','school2'),length.out=length(ids)))
  limited <- confint(f, parm = "slopes",method='sandwich',clusters=clusters)
  expect_true(all(is.na(limited)))
  expect_match(attr(limited,'diagnostics')$InferenceReview[1],'span the free parameters')
  adjusted <- mfrmr:::mfrm_gpcm_inference(f,method='sandwich',adjust=TRUE)
  expect_equal(adjusted$covariance,inf$covariance*length(ids)/(length(ids)-1),tolerance=1e-10)
  expect_error(confint(f, parm = "slopes",method='sandwich',clusters=clusters[-1,]),'every fitted Person')
  # Check the derivative-sum identity away from the optimum as well.
  shifted <- f; shifted$opt$par[1] <- shifted$opt$par[1] + .02
  scores <- mfrmr:::mfrm_mml_person_scores_numeric(shifted)
  expect_true(all(is.finite(scores)))
  expect_identical(dim(scores),c(f$config$n_person,length(f$opt$par)))
})

test_that('one shared information budget governs interval and IC qualification', {
  f <- gpcm_inference_fixture()$fit
  withr::local_options(mfrmr.max_information_bytes=1)
  ci <- confint(f, parm = "slopes")
  expect_true(all(is.na(ci)))
  expect_match(attr(ci,'diagnostics')$InferenceReview[1],'workspace')
  check <- mfrmr:::mfrm_ic_fit_check(f,mfrmr:::mfrm_extract_fit_ic_contract(f))
  expect_false(check$eligible)
  expect_match(check$review,'workspace')
  expect_match(attr(confint(f, parm = "slopes",scale='standardized'),'diagnostics')$InferenceReview[1],'workspace')
  curve <- mfrm_curve_intervals(f,data.frame(Theta=0,Rater='R01',Criterion='C01'))
  expect_false(any(curve$table$CIEligible))
  expect_match(curve$table$InferenceReview[1],'workspace')
  withr::local_options(mfrmr.max_information_bytes=1+1i)
  expect_error(confint(f, parm = "slopes"),'finite positive byte count')
})

test_that('explicit population replay preserves covariates and rejects changed coding', {
  example <- gpcm_inference_fixture(); f <- example$fit
  args <- mfrmr:::mfrmr_gqs_refit_arguments(f,example$data,41L)
  expect_identical(args$population_formula,f$population$formula)
  expect_identical(args$person_data,f$population$person_table_replay)
  expect_equal(args$quad_points,41L)
  # Scope validation now accepts a user-supplied population before refitting.
  expect_no_error(mfrmr:::mfrmr_gqs_validate(f,example$data,c(31,41),c(-3,3),21))
  incomplete <- f; incomplete$population$person_table <- NULL
  incomplete$population$person_table_replay <- NULL
  expect_error(mfrmr:::mfrmr_gqs_refit_arguments(incomplete,example$data,41L),'person data')
})

test_that('curve intervals respect probability and information geometry', {
  f <- gpcm_inference_fixture()$fit
  grid <- expand.grid(Theta=c(-1,0,1),Rater='R01',Criterion='C01')
  prob <- mfrm_curve_intervals(f,grid)
  info <- mfrm_curve_intervals(f,grid,type='information')
  expect_true(all(prob$table$CIEligible)); expect_true(all(info$table$CIEligible))
  expect_true(all(prob$table$Lower>=0 & prob$table$Upper<=1))
  expect_true(all(info$table$Lower>0))
  expect_equal(as.vector(rowsum(prob$table$Estimate,prob$table$InputRow)),rep(1,3),tolerance=1e-12)
  categories <- 0:2
  expected <- vapply(split(prob$table,prob$table$InputRow),function(t) {
    p <- t$Estimate
    f$slopes$Estimate[match('C01',f$slopes$SlopeFacet)]^2 *
      (sum(p*categories^2)-sum(p*categories)^2)
  },numeric(1))
  expect_equal(info$table$Estimate,unname(expected),tolerance=1e-10)
  grid$Rater <- as.character(grid$Rater)
  grid$Rater[1] <- 'unknown'
  expect_error(mfrm_curve_intervals(f,grid),'Unknown')
  skip_if_not_installed('ggplot2')
  p <- plot(prob,draw=FALSE,title=NULL,subtitle=NULL)
  expect_null(p$labels$title); expect_null(p$labels$subtitle)
  expect_equal(p$data$Estimate,prob$table$Estimate)
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that('bootstrap uncertainty retains failed draws and null simulations cannot yield slope CIs', {
  f <- gpcm_inference_fixture()$fit
  draws <- matrix(rep(f$opt$par,each=4),4)
  draws[4,] <- NA_real_
  b <- structure(list(source=f,draws=draws,trials=data.frame(Available=c(TRUE,TRUE,TRUE,FALSE)),
    settings=list(purpose='slope_intervals')),class='mfrm_gpcm_bootstrap')
  ci <- confint(b)
  expect_equal(as.numeric(ci[,1]),rep(0,3))
  expect_true(all(is.infinite(ci[,2])))
  expect_equal(attr(ci,'availability')$Unresolved,rep(1L,3))
  expect_false(any(attr(ci,'diagnostics')$CIEligible))
  test <- mfrmr:::mfrm_gpcm_bootstrap_test(c(0,3,5,NA),4)
  expect_true(is.na(test$PValue)); expect_equal(test$PValueLower,2/5)
  expect_equal(test$PValueUpper,3/5)
  b$settings$purpose <- 'PCM_GPCM_test'
  expect_error(confint(b),'Null-model')
  # Generation keeps all analyzed assignment rows and restores RNG state.
  set.seed(999); before <- .Random.seed
  dat <- mfrmr:::mfrm_gpcm_bootstrap_generate(f,1234)
  expect_identical(.Random.seed,before)
  expect_equal(nrow(dat),nrow(f$prep$data))
  expect_true(all(dat$Score %in% f$prep$score_map$OriginalScore))
  expect_equal(as.integer(table(dat$Person)),as.integer(table(f$prep$data$Person)))
})

test_that('basic bootstrap reverses error quantiles and the runner retains failures', {
  f <- gpcm_inference_fixture()$fit
  at <- mfrmr:::build_param_slices(mfrmr:::build_param_sizes(f$config))$log_slopes[1]
  errors <- seq(-.1,.2,length.out=20)
  draws <- matrix(rep(f$opt$par,each=20),20)
  draws[,at] <- draws[,at]+errors
  b <- structure(list(source=f,draws=draws,trials=data.frame(Available=rep(TRUE,20)),
    settings=list(purpose='slope_intervals')),class='mfrm_gpcm_bootstrap')
  ci <- confint(b,level=.8)
  expected <- exp(f$opt$par[at]-quantile(errors,c(.9,.1),type=1,names=FALSE))
  expect_equal(as.numeric(ci[1,]),expected)
  attempts <- 0L
  local_mocked_bindings(mfrm_gpcm_bootstrap_refit=function(source,data) {
    attempts <<- attempts+1L
    if(attempts==2L) stop('deliberate failed refit')
    source
  },.package='mfrmr')
  set.seed(127); before <- .Random.seed
  result <- bootstrap_mfrm_gpcm(f,nsim=3,seed=131)
  expect_equal(attempts,3L)
  expect_identical(result$trials$Available,c(TRUE,FALSE,TRUE))
  expect_match(result$trials$Reason[2],'deliberate failed refit')
  expect_true(all(is.na(result$draws[2,])))
  expect_identical(result$trials$Stage,c('complete','alternative_refit','complete'))
  expect_identical(result$trials$AlternativeReturned,c(TRUE,FALSE,TRUE))
  expect_true(all(is.na(result$refit_draws[2,])))
  expect_equal(as.numeric(result$refit_draws[1,]),as.numeric(f$opt$par))
  expect_identical(result$checks$Replicate,c(1L,3L))
  expect_true(all(result$checks$WaldEligible))
  expect_equal(result$checks$PopulationSD,rep(sqrt(f$population$sigma2),2))
  standardized <- f$slopes$OptimizerEstimate * sqrt(f$population$sigma2)
  expect_equal(result$checks$MinimumStandardizedSlope,rep(min(standardized),2))
  expect_equal(result$checks$MaximumStandardizedSlope,rep(max(standardized),2))
  expect_equal(apa_table(result,which='checks')$table$Replicate,c(1L,3L))
  expect_identical(.Random.seed,before)
  expect_equal(result$settings$replicate_seeds,result$trials$Seed)
})

test_that('bootstrap diagnostic recording tolerates unavailable slope-scale estimates', {
  f <- gpcm_inference_fixture()$fit
  local_mocked_bindings(compute_mml_parameter_covariance=function(...) stop('Unexpected information calculation'),
    .package='mfrmr')
  f$population$sigma2 <- NULL
  check <- mfrm_gpcm_bootstrap_fit_checks(f,1L,'alternative')
  expect_true(is.na(check$PopulationSD))
  expect_true(is.na(check$MinimumStandardizedSlope))
  f$config$population_spec$active <- FALSE
  check <- mfrm_gpcm_bootstrap_fit_checks(f,1L,'alternative')
  expect_equal(check$PopulationSD,1)
  expect_equal(check$MinimumStandardizedSlope,min(f$slopes$OptimizerEstimate))
  f$slopes <- NULL
  check <- mfrm_gpcm_bootstrap_fit_checks(f,1L,'alternative')
  expect_true(is.na(check$MinimumStandardizedSlope))
  expect_true(is.na(check$MaximumStandardizedSlope))
})


test_that('returned but rejected bootstrap estimates remain diagnostic and never enter intervals', {
  f <- gpcm_inference_fixture()$fit
  information <- mfrm_gpcm_inference(f)
  rejected <- f
  rejected$readiness$fit$CategoryState <- 'weak_information'
  rejected$config$category_support_audit$readiness$ReasonCodes <- 'weak_category_information'
  calls <- 0L
  local_mocked_bindings(
    mfrm_gpcm_bootstrap_refit=function(source,data) {
      calls <<- calls+1L
      if(calls==2L) rejected else source
    },
    mfrm_gpcm_inference=function(fit,...) {
      out <- information
      if(identical(fit$readiness$fit$CategoryState,'weak_information'))
        out$check <- list(eligible=FALSE,review='Category support requires review.')
      out
    },.package='mfrmr')
  result <- bootstrap_mfrm_gpcm(f,nsim=3,seed=9153)
  expect_identical(result$trials$AlternativeReturned,rep(TRUE,3))
  expect_identical(result$trials$Stage,c('complete','slope_eligibility','complete'))
  expect_identical(result$checks$WaldEligible,c(TRUE,FALSE,TRUE))
  expect_identical(result$checks$InformationStatus,rep('ok',3))
  expect_identical(result$checks$CategoryState,c('adequate','weak_information','adequate'))
  expect_true(all(is.finite(result$refit_draws)))
  expect_true(all(is.na(result$draws[2,])))
  ci <- confint(result)
  expect_true(all(is.infinite(ci[,2])))
  expect_equal(attr(ci,'availability')$Unresolved,rep(1L,3))
  res <- mfrm_results(f,intervals=result,include='fit',compute='never')
  expect_identical(res$tables$gpcm_inference_checks,result$checks)
  expect_identical(mfrm_report(res)$tables$gpcm_inference_checks,result$checks)
})

test_that('bootstrap comparison stages retain a returned alternative when the null refit fails', {
  f <- gpcm_inference_fixture()$fit
  null <- f; null$config$model <- 'PCM'
  # The mock PCM generator retains population coordinates after removing slopes.
  slope_indices <- build_param_slices(build_param_sizes(f$config))[["log_slopes"]]
  null$opt$par <- f$opt$par[-slope_indices]
  comparison <- list(comparison_basis=list(lrt_status='computed'),lrt=data.frame(ChiSq=1))
  calls <- 0L
  local_mocked_bindings(compare_mfrm=function(...) comparison,
    mfrm_gpcm_bootstrap_refit=function(source,data) {
      calls <<- calls+1L
      if(identical(source$config$model,'PCM')) stop('deliberate null refit failure')
      source
    },.package='mfrmr')
  result <- bootstrap_mfrm_gpcm(f,nsim=2,seed=9154,null_fit=null)
  expect_identical(result$trials$Stage,rep('null_refit',2))
  expect_true(all(result$trials$AlternativeReturned))
  expect_false(any(result$trials$NullReturned))
  expect_equal(nrow(result$checks),2)
  expect_true(all(is.na(result$checks$WaldEligible)))
  expect_true(all(is.finite(result$refit_draws)))
  expect_true(all(is.na(result$draws)))
  expect_true(is.na(result$test$PValue))
  expect_equal(result$test$Unresolved,2)
})

test_that('singleton bootstrap warnings survive intervals and reports without granting Wald readiness', {
  f <- gpcm_inference_fixture()$fit
  inf <- mfrm_gpcm_inference(f)
  weak <- f
  weak$readiness$fit$CategoryState <- 'weak_information'
  # Isolate the decision/reporting branches; real singleton data are also
  # exercised by the retained-data validation replay, without a new simulation.
  audit <- f$config$category_support_audit
  audit$readiness$CategoryState <- 'weak_information'
  audit$readiness$ReasonCodes <- 'weak_category_information'
  audit$category_table$WithinScopeCount[1] <- 1L
  weak$config$category_support_audit <- audit
  calls <- 0L
  local_mocked_bindings(
    audit_mfrm_category_support=function(...) audit,
    mfrm_gpcm_inference=function(fit,...) {
      out <- inf
      out$check <- mfrm_gpcm_slope_inference_check(fit,out$information)
      out
    },
    mfrm_gpcm_bootstrap_refit=function(source,data) {
      calls <<- calls+1L
      if(calls==1L) weak else source
    },.package='mfrmr')
  expect_warning(result <- bootstrap_mfrm_gpcm(f,nsim=2,seed=9155),'1 bootstrap refit')
  expect_true(all(result$trials$Available))
  expect_true(all(result$checks$BootstrapEligible))
  expect_identical(result$checks$WaldEligible,c(FALSE,TRUE))
  expect_match(result$trials$Warnings[1],'Singleton')
  expect_identical(result$trials$Warnings[2],'')
  expect_true(all(is.finite(result$draws)))
  expect_identical(result$checks$CategoryState,c('weak_information','adequate'))
  ci <- confint(result)
  expect_true(all(is.finite(ci)))
  expect_match(attr(ci,'cautions'),'Singleton')
  expect_match(attr(ci,'diagnostics')$InferenceReview[1],'Singleton')
  expect_output(print(ci),'Caution: Singleton')
  expect_match(apa_table(ci)$table$InferenceReview[1],'Singleton')
  source_table <- apa_table(result,which='source_checks')$table
  expect_identical(source_table$BootstrapEligible,result$source_checks$BootstrapEligible)
  expect_identical(source_table$WaldEligible,result$source_checks$WaldEligible)
  expect_identical(source_table$BootstrapCaution,result$source_checks$BootstrapCaution)
  res <- mfrm_results(f,intervals=result,include='fit',compute='never')
  expect_match(mfrm_report(res)$tables$gpcm_inference_intervals$InferenceReview[1],'Singleton')
  expect_identical(res$tables$gpcm_inference_checks$BootstrapEligible,c(TRUE,TRUE))
  expect_identical(mfrm_gpcm_inference(weak)$check$eligible,FALSE)
  # The source fit follows the same policy. The LRT path keeps the strict gate.
  calls <- 1L
  messages <- character()
  from_weak <- withCallingHandlers(bootstrap_mfrm_gpcm(weak,nsim=2,seed=9156),
    warning=function(w) {messages <<- c(messages,conditionMessage(w)); invokeRestart('muffleWarning')})
  expect_match(messages[1],'Singleton')
  expect_true(from_weak$source_checks$BootstrapEligible)
  expect_false(from_weak$source_checks$WaldEligible)
  expect_identical(from_weak$source$readiness$fit$CategoryState,'weak_information')
  null <- f; null$config$model <- 'PCM'
  expect_error(bootstrap_mfrm_gpcm(weak,nsim=2,seed=9156,null_fit=null),'category support')
  skip_if_not_installed('ggplot2')
  p <- plot(ci,draw=FALSE)
  expect_match(paste(plot_data(p, "notes")$Text, collapse = "\n"),'Singleton')
  expect_null(plot(ci,subtitle=NULL,draw=FALSE)$labels$subtitle)
  expect_identical(plot_data(plot(ci,subtitle='Custom',draw=FALSE), "subtitle"),'Custom')
  expect_no_error(ggplot2::ggplot_build(p))
})
