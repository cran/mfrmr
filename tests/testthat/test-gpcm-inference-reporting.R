gpcm_reporting_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      d <- simulate_mfrm_data(n_person=80,n_rater=3,n_criterion=3,score_levels=3,
        model='GPCM',step_facet='Criterion',slope_facet='Criterion',seed=924090)
      f <- fit_mfrm(d,'Person',c('Rater','Criterion'),'Score',model='GPCM',method='MML',
        step_facet='Criterion',slope_facet='Criterion',quad_points=31,maxit=400,reltol=1e-10)
      ci <- confint(f, parm = "slopes",scale='standardized',level=.995,simultaneous='bonferroni')
      grid <- data.frame(Theta=seq(-2,2,length.out=7),Rater='R01',Criterion='C01')
      cached <<- list(fit=f,ci=ci,curves=mfrm_curve_intervals(f,grid,level=.9))
    }
    cached
  }
})

test_that('saved inference carries targets and exact sources through tables and reports', {
  x <- gpcm_reporting_fixture()
  res <- mfrm_results(x$fit,intervals=list(slopes=x$ci,curves=x$curves),
    include=c('fit','plots'),compute='never')
  expect_identical(res$gpcm_inference$slopes,x$ci)
  code <- paste(summary(res)$reproducible_code$Code,collapse='\n')
  expect_match(code,'saveRDS',fixed=TRUE)
  expect_match(code,'readRDS',fixed=TRUE)
  expect_identical(mfrm_results_viewer_payload(res)$replay_code,code)
  expect_identical(res$fit$slopes,x$fit$slopes)
  expect_identical(res$tables$gpcm_slopes_intervals$ConfidenceLevel,rep('99.5%',3))
  expect_true(all(res$tables$gpcm_slopes_intervals$Adjustment=='Bonferroni'))
  expect_true(all(c('gpcm_slopes','gpcm_curves') %in% res$plot_map$Type))
  expect_equal(apa_table(x$ci,digits=10)$table$Estimate,
    attr(x$ci,'diagnostics')$Estimate,tolerance=1e-9)
  expect_identical(apa_table(x$ci)$table$ConfidenceLevel,rep('99.5%',3))
  expect_identical(apa_table(x$curves)$table$ConfidenceLevel,rep('90%',21))
  report <- mfrm_report(res)
  expect_identical(report$tables$gpcm_slopes_intervals,res$tables$gpcm_slopes_intervals)
  expect_identical(report$tables$gpcm_curves_curves,res$tables$gpcm_curves_curves)
  expect_match(report$markdown,'GPCM uncertainty')
  contrasts <- matrix(c(1,-1,0),1,dimnames=list('first/second',x$fit$slopes$SlopeFacet))
  # Explicit names specify which fitted slopes the comparison combines.
  ratio <- confint(x$fit, parm = "slopes",contrasts=contrasts,method='sandwich',adjust=TRUE)
  tables <- mfrmr:::mfrm_gpcm_inference_tables(ratio)
  expect_equal(as.matrix(tables$contrasts[-1]),contrasts,ignore_attr=TRUE)
  expect_equal(tables$settings$IndependentClusters,80)
  expect_equal(tables$settings$SmallSampleFactor,80/79)
  expect_equal(nrow(tables$clusters),80)
  altered <- x$fit; altered$opt$par[1] <- altered$opt$par[1]+.01
  expect_error(mfrm_results(altered,intervals=x$ci,compute='never'),'must match')
  altered <- x$fit; altered$prep$data$Score[1] <- 99
  expect_error(mfrm_results(altered,intervals=x$curves,compute='never'),'must match')
  altered <- x$fit; altered$config$estimation_control$quad_points <- 81
  expect_error(mfrm_results(altered,intervals=x$ci,compute='never'),'must match')
  expect_error(mfrm_results(x$fit,intervals=list(A=x$ci,a=x$ci)),'distinct')
})

test_that('plot extraction and customization retain unavailable intervals without recalculation', {
  skip_if_not_installed('ggplot2')
  x <- gpcm_reporting_fixture()
  local_mocked_bindings(fit_mfrm=function(...) stop('unexpected refit'),
    bootstrap_mfrm_gpcm=function(...) stop('unexpected resampling'),
    compute_mml_parameter_covariance=function(...) stop('unexpected information calculation'),.package='mfrmr')
  p <- as_ggplot(x$ci,title=NULL,subtitle=NULL)
  expect_s3_class(p,'ggplot')
  expect_null(p$labels$title); expect_null(p$labels$subtitle)
  expect_equal(plot_data(p,'table')$Estimate,attr(x$ci,'diagnostics')$Estimate)
  expect_equal(plot_data(x$curves,'table')$Estimate,x$curves$table$Estimate)
  expect_identical(plot_data(x$curves,'settings')$level,.9)
  expect_match(paste(plot_data(as_ggplot(x$curves), "notes")$Text, collapse = "\n"),'90%')
  bad <- x$ci; tab <- attr(bad,'diagnostics')
  tab$CI_Lower[1:2] <- c(0,NA); tab$CI_Upper[1:2] <- c(Inf,NA)
  tab$CIEligible[1:2] <- FALSE; attr(bad,'diagnostics') <- tab
  shown <- plot_data(bad,'table')
  expect_identical(shown$Display,c('Unbounded','Unavailable','Available'))
  expect_identical(shown$CI_Upper[1],Inf)
  expect_no_error(ggplot2::ggplot_build(as_ggplot(bad)))
  expect_error(plot(x$ci,reference=Inf,draw=FALSE),'finite')
})

test_that('curve plots expose unavailable intervals and break ribbons at gaps', {
  skip_if_not_installed('ggplot2')
  x <- gpcm_reporting_fixture()$curves
  local_mocked_bindings(compute_mml_parameter_covariance=function(...) stop('unexpected recalculation'),
    fit_mfrm=function(...) stop('unexpected refit'),.package='mfrmr')
  middle <- which(x$table$Theta==0 & x$table$Category==1)
  expect_length(middle,1L)
  x$table$CIEligible[middle] <- FALSE
  x$table$Lower[middle] <- x$table$Upper[middle] <- NA_real_
  x$table$InferenceReview[middle] <- 'Interval unavailable at this point.'
  # An unsorted input must still break the ribbon on the ability axis.
  x$table <- x$table[nrow(x$table):1,]
  p <- plot(x,draw=FALSE)
  tab <- plot_data(p,'table')
  one <- tab[tab$Category==1 & tab$CIEligible,]
  expect_length(unique(one$IntervalGroup),2L)
  expect_identical(tab$InferenceReview,x$table$InferenceReview)
  expect_match(paste(plot_data(p, "notes")$Text, collapse = "\n"),'1 of 21')
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$data[[3]]),1L)
  expect_equal(built$data[[3]]$shape,4)
  hidden <- plot(x,title=NULL,subtitle=NULL,caption=NULL,draw=FALSE)
  expect_null(hidden$labels$title); expect_null(hidden$labels$subtitle)
  expect_null(hidden$labels$caption)
  expect_equal(nrow(ggplot2::ggplot_build(hidden)$data[[3]]),1L)
  x$table$CIEligible[] <- FALSE
  x$table$Lower[] <- x$table$Upper[] <- NA_real_
  all_missing <- plot(x,draw=FALSE)
  expect_match(paste(plot_data(all_missing, "notes")$Text, collapse = "\n"),'21 of 21')
  expect_equal(nrow(ggplot2::ggplot_build(all_missing)$data[[3]]),21L)
  expect_match(all_missing$labels$alt,'Crosses')
  x$cautions <- paste(rep('Weak information requires careful review.',8),collapse=' ')
  expect_true(max(nchar(strsplit(plot_data(plot(x,draw=FALSE), "subtitle"),'\n')[[1]]))<=75)
})

test_that('saved GPCM exports replay the exact selected inference', {
  x <- gpcm_reporting_fixture()
  res <- mfrm_results(x$fit,intervals=list(slopes=x$ci,curves=x$curves),
    include=c('fit','plots'),compute='never')
  folder <- tempfile('gpcm-export-'); on.exit(unlink(folder,recursive=TRUE))
  local_mocked_bindings(fit_mfrm=function(...) stop('unexpected refit'),
    bootstrap_mfrm_gpcm=function(...) stop('unexpected resampling'),.package='mfrmr')
  warnings <- character()
  exported <- withCallingHandlers(export_mfrm_results(res,output_dir=folder,
    include=c('tables','replay','report','plots'),preset='starter',acknowledge_sensitive=TRUE),
    warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')})
  expect_length(warnings,3)
  expect_true(all(startsWith(warnings,'Review-only display:')))
  expect_s3_class(as_ggplot(res,type='gpcm_slopes'),'ggplot')
  expect_s3_class(as_ggplot(res,type='gpcm_curves'),'ggplot')
  files <- list.files(folder,recursive=TRUE,full.names=TRUE)
  saved <- readRDS(files[grepl('\\.rds$',files)][1])
  expect_identical(saved$gpcm_inference,res$gpcm_inference)
  replay <- files[grepl('replay.*\\.R$',files)]
  expect_length(replay,1)
  index <- paste(readLines(file.path(folder,'index.html')),collapse='\n')
  expect_match(index,'Saved inference',fixed=TRUE)
  expect_match(index,'src="mfrmr_results_plot_gpcm_slopes.png"',fixed=TRUE)
  expect_match(index,'Saved GPCM inference',fixed=TRUE)
  expect_match(paste(readLines(replay),collapse='\n'),'readRDS')
  expect_false(any(grepl('gpcm_',exported$plot_errors$Type)))
  expect_true(any(grepl('gpcm_slopes.*\\.png$',files)))
  expect_equal(plot_data(saved,type='gpcm_slopes',component='table')$Estimate,
    attr(x$ci,'diagnostics')$Estimate)
})

test_that('report Markdown retains interval status and saved inference reasons', {
  x <- gpcm_reporting_fixture()
  ci <- x$ci; tab <- attr(ci,'diagnostics')
  ci[2:3,] <- rbind(c(0,Inf),c(NA_real_,NA_real_))
  tab$CI_Lower[2:3] <- c(0,NA_real_); tab$CI_Upper[2:3] <- c(Inf,NA_real_)
  tab$CIEligible[2:3] <- FALSE
  tab$InferenceReview <- c('Weak information: review the interval width.',
    'Unresolved bootstrap draws prevent finite bounds.','Joint information unavailable.')
  attr(ci,'diagnostics') <- tab
  local_mocked_bindings(compute_mml_parameter_covariance=function(...) stop('unexpected recalculation'),
    fit_mfrm=function(...) stop('unexpected refit'),.package='mfrmr')
  res <- mfrm_results(x$fit,intervals=ci,include='fit',compute='never')
  report <- mfrm_report(res)
  expect_match(report$markdown,'1 finite; 1 unbounded; 1 unavailable',fixed=TRUE)
  for(reason in tab$InferenceReview)expect_match(report$markdown,reason,fixed=TRUE)
  expect_identical(report$tables$gpcm_inference_intervals$CI_Upper,tab$CI_Upper)
})

test_that('bootstrap test plots retain all unresolved replicates', {
  skip_if_not_installed('ggplot2')
  x <- structure(list(test=data.frame(LR=1),
    trials=data.frame(Available=c(TRUE,FALSE),LR=c(.2,NA_real_)),
    settings=list(sampling='Independent Persons',seed=1)),class='mfrm_gpcm_bootstrap')
  expect_identical(plot_data(x,'table'),x$trials)
  expect_equal(nrow(apa_table(x,which='trials')$table),2)
  expect_match(paste(plot_data(as_ggplot(x), "notes")$Text, collapse = "\n"),'unresolved: 1')
  x$trials$Available[] <- FALSE; x$trials$LR[] <- NA_real_
  expect_no_error(ggplot2::ggplot_build(as_ggplot(x)))
})

test_that('bootstrap diagnostic tables preserve near-zero values through APA and kable', {
  tab <- data.frame(BootstrapEligible=FALSE,WaldEligible=FALSE,
    PopulationSD=.0042,MinimumStandardizedSlope=1.583352e-5,MaximumStandardizedSlope=3.161197,
    GradientMaxAbs=7.804191e-5,SmallestEigenvalue=8.544122e-7,InformationScale=97604.53,
    EmptyCategories=0L)
  before <- tab
  out <- apa_table(tab,digits=2)
  expect_identical(tab,before)
  expect_equal(as.numeric(out$table$MinimumStandardizedSlope),1.6e-5)
  expect_equal(as.numeric(out$table$SmallestEigenvalue),8.5e-7)
  expect_equal(as.numeric(out$table$GradientMaxAbs),7.8e-5)
  expect_equal(out$table$EmptyCategories,0)
  expect_true(all(as.numeric(apa_table(tab,digits=0)$table$MinimumStandardizedSlope)>0))
  skip_if_not_installed('knitr')
  html <- as.character(as_kable(out,format='html'))
  expect_match(html,'1.6e-05',fixed=TRUE)
  expect_match(html,'8.5e-07',fixed=TRUE)
})

test_that('selected bootstrap updates remain identified through saved output', {
  f <- gpcm_reporting_fixture()$fit
  local_mocked_bindings(fit_mfrm=function(...) stop('unexpected refit'),
    bootstrap_mfrm_gpcm=function(...) stop('unexpected resampling'),.package='mfrmr')
  b <- structure(list(source=f, draws=rbind(f$opt$par,f$opt$par,NA_real_),
    trials=data.frame(Available=c(TRUE,TRUE,FALSE)),
    settings=list(purpose='slope_intervals',sampling='Independent Persons',seed=17)),
    class='mfrm_gpcm_bootstrap')
  original <- confint(b)
  b$settings$recheck <- 'Original draws with selected refit updates.'
  b$settings$diagnostic_checks_scope <- 'Only selected refits were rechecked.'
  before <- serialize(b,NULL)
  ci <- confint(b)
  expect_identical(summary(b)$intervals, ci)
  expect_identical(summary(b)$settings, b$settings)
  expect_equal(summary(b)$overview$Unresolved, 1L)
  expect_identical(as.numeric(ci),as.numeric(original))
  expect_identical(attr(ci,'availability'),attr(original,'availability'))
  expect_match(paste(capture.output(print(b)),collapse=' '),'not a complete rerun')
  expect_true(all(grepl('not a complete rerun',attr(ci,'diagnostics')$InferenceReview)))
  expect_match(paste(plot_data(plot(ci,draw=FALSE), "notes")$Text, collapse = "\n"),'not a complete')
  expect_null(plot(ci,subtitle=NULL,draw=FALSE)$labels$subtitle)
  expect_identical(apa_table(b,which='sampling')$table$Reanalysis,b$settings$recheck)
  tables <- mfrmr:::mfrm_gpcm_inference_tables(ci)
  expect_identical(tables$settings$DiagnosticChecksScope,b$settings$diagnostic_checks_scope)
  res <- mfrm_results(f,intervals=list(bootstrap=b,intervals=ci),include='fit',compute='never')
  report <- mfrm_report(res)
  expect_match(report$markdown,'not a complete rerun')
  expect_identical(report$tables$gpcm_bootstrap_sampling$Reanalysis,b$settings$recheck)
  expect_identical(report$tables$gpcm_intervals_settings$Reanalysis,b$settings$recheck)
  expect_identical(serialize(b,NULL),before)
  # A saved null-test result uses its print/table path, not slope intervals.
  b$test <- data.frame(PValue=NA_real_)
  b$settings$purpose <- 'PCM_GPCM_test'
  expect_match(paste(capture.output(print(b)),collapse=' '),'not a complete rerun')
  expect_identical(mfrmr:::mfrm_gpcm_inference_tables(b)$sampling$Reanalysis,b$settings$recheck)
})


test_that("saved curve and bootstrap summaries retain unavailable targets and history", {
  curve <- structure(list(settings = list(level = .95, point_only = TRUE),
    table = data.frame(Estimate = .4, Lower = NA_real_, Upper = NA_real_,
      CIEligible = FALSE, InferenceReview = "Provisional curves only"),
    contexts = data.frame(ObservedContext = FALSE), cautions = "Retained caution",
    anchor_note = "Conditional on fixed anchors"), class = "mfrm_curve_intervals")
  out <- summary(curve)
  expect_identical(out$table, curve$table)
  expect_identical(out$settings, curve$settings)
  expect_identical(out$contexts, curve$contexts)
  expect_identical(out$cautions, curve$cautions)
  expect_identical(out$anchor_note, curve$anchor_note)
  boot <- structure(list(trials = data.frame(Available = c(TRUE, FALSE, NA)),
    settings = list(purpose = "PCM_GPCM_test", recheck = list(replicate = 2)),
    checks = data.frame(BootstrapEligible = c(TRUE, FALSE),
      BootstrapCaution = c("Weak information", "Rejected trial")),
    test = data.frame(PValue = NA_real_, Reason = "Unresolved trials")),
    class = "mfrm_gpcm_bootstrap")
  s <- summary(boot)
  expect_equal(s$overview, data.frame(Planned = 3L, Available = 1L, Unresolved = 2L))
  expect_identical(s$settings, boot$settings)
  expect_identical(s$test, boot$test)
  expect_null(s$intervals)
  expect_identical(s$cautions, "Weak information")
  expect_error(summary(curve, unused = TRUE), "must be empty")
  expect_error(summary(boot, unused = TRUE), "must be empty")
})
