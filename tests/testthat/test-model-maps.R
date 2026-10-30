model_map_fixture <- local({
  saved <- NULL
  function() {
    if(!is.null(saved)) return(saved)
    d <- expand.grid(Person=paste0('P',1:12),Rater=c('A','B'),Criterion=c('X','Y'))
    d$Score <- c(0,1,1,0,1,0,1,1,0,0,1,0, 1,0,0,1,0,1,1,0,1,0,0,1,
      1,0,1,0,1,0,1,0,1,1,0,1, 0,1,0,1,0,0,1,1,0,1,0,0)
    o <- fit_mfrm(d,'Person',c('Rater','Criterion'),'Score',quad_points=61)
    # Fixed nonzero-mean/nonunit population reference. No inferential fit claim.
    o$population <- list(active=TRUE,design_columns='(Intercept)',coefficients=.6,sigma2=.81)
    o$config$population_spec$design_matrix <- matrix(1,12,1)
    input <- mfrm_testlet_data(d,'Person','Score','Rater',c('Rater','Criterion'),0:1,'omit')
    a <- mfrm_ordinary_response_input(o);beta <- as.vector(qr.solve(input$X,a$offset))
    fixed <- as.data.frame(o$facets$others[c('Facet','Level','Estimate')]);fixed$Parameter<-'Fixed facet'
    fixed$SE<-fixed$Lower<-fixed$Upper<-NA_real_
    e <- structure(list(input=input,parameters=c(beta,a$steps-a$mean,0,a$sd^2),
      calibration=list(beta=beta,steps=a$steps-a$mean,variance=0,person_variance=a$sd^2,person_sd=a$sd),
      calibration_table=fixed,loglik=o$summary$LogLik,
      checks=list(NumericalReady=TRUE,InformationPositive=TRUE,EstimatedVarianceBoundary=FALSE),
      settings=list(quad_points=61L,missing='omit',method='Fixed numerical reference')),class='mfrm_testlet')
    saved <<- list(o=o,e=e,data=d)
    saved
  }
})

test_that('ordinary conditional scores agree with independent moments and continuous tails',{
  f<-model_map_fixture();a<-mfrm_ordinary_response_input(f$o)
  s<-score_mfrm_persons(f$o,persons='P1')
  ids<-which(a$person=='P1')
  density<-function(t) vapply(t,function(theta) {
    p<-plogis(theta-a$offset[ids]-a$steps)
    prod(ifelse(a$y[ids]==1,p,1-p))*dnorm(theta,a$mean,a$sd)
  },numeric(1))
  mass<-integrate(density,-Inf,Inf,rel.tol=1e-10,abs.tol=0)$value
  mu<-integrate(function(t) t*density(t),-Inf,Inf,rel.tol=1e-10,abs.tol=0)$value/mass
  variance<-integrate(function(t) (t-mu)^2*density(t),-Inf,Inf,rel.tol=1e-10,abs.tol=0)$value/mass
  expect_equal(s$table$Estimate,mu,tolerance=1e-7)
  expect_equal(s$table$ConditionalSD,sqrt(variance),tolerance=1e-7)
  expect_equal(integrate(density,-Inf,s$table$Lower,rel.tol=1e-10,abs.tol=0)$value/mass,.025,tolerance=1e-7)
  expect_equal(integrate(density,s$table$Upper,Inf,rel.tol=1e-10,abs.tol=0)$value/mass,.025,tolerance=1e-7)
  expect_s3_class(s,'mfrm_person_scores')
  expect_equal(nrow(s$scoring_data),48)
  expect_s3_class(as_ggplot(s),'ggplot')
  expect_error(score_mfrm_persons(f$o,persons='Absent'),'Person IDs')
  expect_error(score_mfrm_persons(f$o,quad_points=31),'continuous integration')
})

test_that('same-roster zero-local scores align by population origin and retain uncertainty',{
  f<-model_map_fixture()
  a<-score_mfrm_persons(f$o,persons=c('P2','P1'))
  b<-score_mfrm_persons(f$e,persons=c('P1','P2'))
  expect_equal(nrow(b$scoring_data),48)
  x<-compare_mfrm(f$o,f$e,person_scores=list(a,b))
  expect_equal(x$persons$table$Reference,x$persons$table$Comparison,tolerance=1e-7)
  expect_equal(x$persons$table$LowerReference,x$persons$table$LowerComparison,tolerance=1e-7)
  expect_equal(x$persons$table$ConditionalSDReference,x$persons$table$ConditionalSDComparison,tolerance=1e-7)
  expect_equal(x$persons$table$OriginReference,rep(.6,2))
  reverse<-compare_mfrm(f$e,f$o,person_scores=list(b,a))
  expect_equal(reverse$persons$table$Difference,-x$persons$table$Difference,tolerance=1e-12)
  expect_null(x$persons$DifferenceSE)
  expect_error(compare_mfrm(f$o,f$e,person_scores=list(a,predict(f$e,persons='P1'))),'same requested Persons')
  reduced<-predict(f$e,f$data[f$data$Person %in% c('P1','P2'),])
  expect_error(compare_mfrm(f$o,f$e,person_scores=list(a,reduced)),'complete source roster')
  bad<-b;bad$settings$level<-.9
  expect_error(compare_mfrm(f$o,f$e,person_scores=list(a,bad)),'interval level')
  old<-b;old$scoring_data<-NULL
  expect_error(compare_mfrm(f$o,f$e,person_scores=list(a,old)),'older scores')
  p<-plot(x,metric='person',draw=FALSE)
  expect_match(p$data$alt_text,'conditional Person')
  expect_equal(p$data$table$Y,x$persons$table$Comparison)
})

test_that('locations encode adjacent crossings without borrowing facet intervals',{
  f<-model_map_fixture();s<-score_mfrm_persons(f$e,persons=c('P1','P2'))
  # Multiple, non-centered steps exercise the declared reference transformation.
  e<-f$e;e$calibration$steps<-c(-.4,.6,1.3);e$input$score_levels<-1:4
  loc<-mfrm_model_locations(e)
  z<-loc[loc$Kind=='Fixed facet',]
  expect_equal(z$Estimate,z$SourceEstimate+mean(c(-.4,.6,1.3)))
  expect_true(all(is.na(z$Lower)))
  for(i in seq_len(nrow(z))) {
    thresholds<-z$SourceEstimate[i]+e$calibration$steps
    ratios<-exp(thresholds-z$SourceEstimate[i]-e$calibration$steps)
    expect_equal(ratios,rep(1,3))
    expect_equal(mean(thresholds),z$Estimate[i])
  }
  expect_equal(loc$Estimate[loc$Kind=='Category boundary'],c(-.4,.6,1.3))
  random<-e;class(random)<-'mfrm_random_rater';random$input$columns$rater<-'Rater'
  random$calibration_table<-random$calibration_table[random$calibration_table$Facet=='Criterion',]
  random$raters<-data.frame(Rater=c('A','B'),Estimate=c(-.2,.3),Lower=-1,Upper=1)
  rloc<-mfrm_model_locations(random)
  expect_equal(rloc$Estimate[rloc$Kind=='Observed rater'],c(-.2,.3)+mean(e$calibration$steps))
  expect_true(all(is.na(rloc$Lower[rloc$Kind=='Observed rater'])))
  bad<-f$e;bad$checks$NumericalReady<-FALSE
  expect_error(mfrm_model_locations(bad,s),'numerical checks')
})

test_that('maps use saved matching evidence and retain absent groups and prior-only Persons',{
  f<-model_map_fixture();d<-f$data;extra<-d[1,];extra$Person<-'Prior';extra$Score<-NA;d<-rbind(d,extra)
  e<-f$e;e$input<-mfrm_testlet_data(d,'Person','Score','Rater',c('Rater','Criterion'),0:1,'omit')
  s<-score_mfrm_persons(e,persons=c('P1','Prior'))
  diag<-mfrm_response_diagnostics(e,group_by=c('Person','Rater'))
  res<-mfrm_results(e,scores=s,response_diagnostics=diag)
  fail<-function(...) stop('Unexpected recomputation')
  local_mocked_bindings(score_mfrm_persons=fail,predict.mfrm_testlet=fail,
    mfrm_response_diagnostics=fail,mfrm_testlet_person_interval=fail,.package='mfrmr')
  w<-plot(res,type='wright',draw=FALSE)
  p<-plot(res,type='fit_pathway',facet='Person',draw=FALSE)
  expect_equal(w$data$table$Shape[w$data$table$Level=='Prior'],1)
  expect_false(p$data$table$Included[p$data$table$Level=='Prior'])
  expect_equal(p$data$table$Missing[p$data$table$Level=='Prior'],1)
  expect_true(all(w$data$table$IntervalShown == (w$data$table$Kind %in% c('Person EAP','Prior only'))))
  expect_equal(p$data$table$X[p$data$table$Level=='P1'],diag$measures$Infit[diag$measures$Facet=='Person' & diag$measures$Level=='P1'])
  expect_false(any(c('ZSTD','Flag','Cutoff') %in% names(p$data$table)))
  expect_match(p$data$alt_text,'no fit cutoffs')
  partial<-plot(res,type='fit_pathway',facet=c('Person','Rater'),persons='Prior',draw=FALSE)
  expect_equal(partial$data$person_selection$SelectedPersons,1)
  expect_equal(partial$data$person_selection$DisplayedPersons,0)
  built<-ggplot2::ggplot_build(as_ggplot(partial))
  expect_equal(nrow(built$layout$layout),2)
  expect_true(any(vapply(built$data,function(z) 'label' %in% names(z) &&
    any(z$label=='No available locations'),logical(1))))
  expect_error(plot(res,type='fit_pathway',facet='Criterion',draw=FALSE),'available location')
  expect_error(plot(res,type='wright',persons='Absent',draw=FALSE),'saved Person')
  no_labels<-plot(res,type='wright',draw=FALSE,show_intervals=FALSE,show_steps=FALSE,show_title=FALSE,show_notes=FALSE,palette='mono')
  expect_false(any(no_labels$data$table$IntervalShown))
  expect_false(any(no_labels$data$table$Kind=='Category boundary'))
  expect_true(all(no_labels$data$table$Colour=='#222222'))
  expect_null(as_ggplot(no_labels)$labels$title)
  for(type in c('wright','fit_pathway')) {
    raw<-plot(res,type=type,draw=FALSE)
    expect_silent(ggplot2::ggplotGrob(as_ggplot(raw)))
  }
  grDevices::pdf(tempfile(fileext='.pdf'));on.exit(grDevices::dev.off(),add=TRUE)
  before<-graphics::par(c('mar','oma','mfrow','cex'))
  expect_silent(plot(res,type='wright'));expect_silent(plot(res,type='fit_pathway',facet='Person'))
  expect_equal(graphics::par(c('mar','oma','mfrow','cex')),before)
  exported<-export_mfrm_results(res,output_dir=tempfile(),preset='starter',acknowledge_sensitive=TRUE)
  expect_equal(nrow(exported$plot_errors),0)
  expect_true(all(c('plot_wright','plot_fit_pathway') %in% exported$written_files$Component))
  expect_identical(mfrm_report(res)$tables$model_locations,res$tables$model_locations)
  expect_error(mfrm_results(e,scores=s,predictions=s),'once')
})

test_that('ordinary prior-only and zero-variance Persons never imply measured certainty',{
  f<-model_map_fixture();d<-f$data;extra<-d[1,];extra$Person<-'Prior';extra$Score<-NA;d<-rbind(d,extra)
  o<-f$o;o$prep<-prepare_mfrm_data(d,'Person',c('Rater','Criterion'),'Score')
  s<-score_mfrm_persons(o,persons='Prior')
  expect_equal(s$table$Status,'prior_only');expect_equal(s$table$Estimate,.6)
  expect_equal(s$table$ConditionalSD,.9)
  o$population$sigma2<-0
  zero<-score_mfrm_persons(o,persons=c('P1','Prior'))
  expect_true(all(zero$table$Status=='unavailable'))
  expect_true(all(is.na(zero$table$Lower)))
})
