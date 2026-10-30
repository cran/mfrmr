comparison_fixture <- local({
  saved <- NULL
  function() {
    if (!is.null(saved)) return(saved)
    d <- expand.grid(Person = paste0('P',1:12), Rater = c('A','B'), Criterion = c('X','Y'))
    d$Score <- c(0,1,0,0,1,1,0,1,1,0,1,0, 1,1,1,0,0,1,0,0,1,0,0,1,
                 1,0,1,0,1,0,1,0,1,1,0,1, 0,1,0,1,0,0,1,1,0,1,0,0)
    o <- fit_mfrm(d,'Person',c('Rater','Criterion'),'Score',method='MML',model='RSM',quad_points=31)
    expect_identical(o$summary$NumericalState, 'ready')
    input <- mfrm_random_rater_data(d,'Person','Rater','Criterion','Score',0:1,'fail')
    e <- structure(list(input=input, calibration=list(person_sd=1,person_variance=1,rater_sd=.7),
      calibration_table=data.frame(Parameter='Fixed facet',Facet='Criterion',Level=c('X','Y'),Estimate=c(-.4,.4)),
      raters=data.frame(Rater=c('A','B'),Estimate=c(-2,0)), parameters=numeric(),
      settings=list(fixed_person_sd=1,method='MML: shared-rater Laplace',fixed_rater_sd=.7),
      checks=list(NumericalReady=TRUE,InformationPositive=TRUE),loglik=-30),class='mfrm_random_rater')
    o$facets$others$Estimate <- c(3,5,-.7,.7)
    saved <<- list(ordinary=o,extension=e,data=d)
    saved
  }
})

test_that('comparison verifies matched events and preserves descriptive centering', {
  f <- comparison_fixture(); x <- compare_mfrm(f$ordinary,f$extension,labels=c('Fixed','Random'))
  expect_s3_class(x,'mfrm_extended_comparison')
  expect_false(inherits(x,'mfrm_comparison'))
  expect_null(x$preferred);expect_null(x$lrt)
  r <- subset(x$effects,Facet=='Rater');c <- subset(x$effects,Facet=='Criterion')
  expect_equal(r$Reference,c(-1,1));expect_equal(r$Comparison,c(-1,1));expect_equal(r$Difference,c(0,0))
  expect_equal(c$Difference,c(.3,-.3))
  expect_equal(r$CenterReference,c(4,4));expect_equal(r$CenterComparison,c(-1,-1))
  expect_identical(r$ComparisonKind,rep('Conditional rater mode',2))
  reverse <- compare_mfrm(f$extension,f$ordinary)
  expect_equal(reverse$effects$Difference,-x$effects$Difference)
  shifted <- f$extension;shifted$raters$Estimate <- shifted$raters$Estimate+7
  expect_equal(compare_mfrm(f$ordinary,shifted)$effects$Difference,x$effects$Difference)
  permuted <- f$ordinary;permuted$prep$data <- permuted$prep$data[nrow(permuted$prep$data):1,]
  expect_identical(compare_mfrm(permuted,f$extension)$effects,x$effects)
  expect_identical(summary(x)$effects,x$effects)
  file <- tempfile();saveRDS(x,file);expect_identical(readRDS(file),x)
})

test_that('equal counts do not bypass event, category, population or constraint checks', {
  f <- comparison_fixture(); o<-f$ordinary;e<-f$extension
  bad <- e;bad$input$data$Score[1] <- 1-bad$input$data$Score[1]
  expect_error(compare_mfrm(o,bad),'Observed rating events differ')
  # Replacing an event with another changes multiplicity, even at equal N.
  bad <- e;bad$input$data[1,] <- bad$input$data[2,]
  expect_error(compare_mfrm(o,bad),'multiplicities')
  bad <- e;bad$input$score_levels <- 1:2
  expect_error(compare_mfrm(o,bad),'categories')
  bad <- e;bad$settings$fixed_person_sd <- NULL
  expect_error(compare_mfrm(o,bad),'population assumption')
  bad <- e;bad$calibration$person_sd <- .7
  expect_error(compare_mfrm(o,bad),'population assumption')
  bad <- o;bad$config$method <- 'JML'
  expect_error(compare_mfrm(bad,e),'RSM MML')
  bad <- o;bad$prep$data$Weight[1] <- 2
  expect_error(compare_mfrm(bad,e),'unit weights')
  bad <- o;bad$config$facet_signs['Rater'] <- 1
  expect_error(compare_mfrm(bad,e),'severity facets')
  expect_error(compare_mfrm(o,e,nested=TRUE),'descriptive')
  expect_error(compare_mfrm(o,e,labels=c('a','a')),'distinct')
  expect_error(compare_mfrm(e,e),'exactly one')
  expect_error(compare_mfrm(o,e,o),'exactly one')
})

test_that('omission identities are recorded before filtering and compared', {
  f <- comparison_fixture(); d <- f$data;d$Score[c(1,5)] <- NA
  prep <- prepare_mfrm_data(d,'Person',c('Rater','Criterion'),'Score')
  expect_identical(prep$omitted_input_rows,c(1L,5L))
  expect_equal(as.character(prep$omitted_data$Person),c('P1','P5'))
  expect_true(all(is.na(prep$omitted_data$Score)))
  o<-f$ordinary;o$prep<-prep
  e<-f$extension;e$input<-mfrm_random_rater_data(d,'Person','Rater','Criterion','Score',0:1,'omit')
  x<-compare_mfrm(o,e)
  expect_equal(x$models$Omitted,c(2,2));expect_equal(nrow(x$omitted_events),2)
  bad<-o;bad$prep$omitted_data$Person[1]<-'Different'
  expect_error(compare_mfrm(bad,e),'Omitted rating-event identities differ')
  old<-o;old$prep$omitted_data<-NULL
  expect_error(compare_mfrm(old,e),'older fit')
  bad<-o;bad$prep$omitted_data$Weight[1]<-0
  expect_error(compare_mfrm(bad,e),'excluded weights')
  d$Weight<-1;d$Weight[8]<-0
  p<-prepare_mfrm_data(d,'Person',c('Rater','Criterion'),'Score',weight_col='Weight')
  expect_identical(p$omitted_input_rows,c(1L,5L,8L))
  expect_equal(nrow(p$data),nrow(d)-3)
})

test_that('unavailable outputs stay explicit and centering uses complete level sets', {
  f<-comparison_fixture();e<-f$extension;e$raters$Estimate[1]<-NA_real_
  x<-compare_mfrm(f$ordinary,e)
  expect_true(all(is.na(subset(x$effects,Facet=='Rater')$Difference)))
  expect_true(all(is.finite(subset(x$effects,Facet=='Criterion')$Difference)))
  e<-f$extension;e$checks$NumericalReady<-FALSE
  x<-compare_mfrm(f$ordinary,e)
  expect_true(all(is.na(x$effects$Difference)))
  expect_true(all(is.finite(x$effects$SourceComparison)))
  expect_true('review' %in% x$checks$Status)
  blocked<-f$ordinary;blocked$readiness$fit$FitReadiness<-'blocked'
  expect_true(all(is.na(compare_mfrm(blocked,f$extension)$effects$Difference)))
  flagged<-f$ordinary;flagged$facets$others$ParameterStatus<-'estimable'
  flagged$facets$others$ParameterStatus[1]<-'unbounded_high'
  expect_true(all(is.na(subset(compare_mfrm(flagged,f$extension)$effects,Facet=='Rater')$Difference)))
  before<-dev.cur();raw<-plot(x,draw=FALSE,show_title=FALSE,show_notes=FALSE)
  expect_identical(dev.cur(),before)
  if(requireNamespace('ggplot2',quietly=TRUE)) {
    p<-as_ggplot(raw);expect_null(p$labels$title);expect_null(p$labels$caption)
    expect_silent(ggplot2::ggplotGrob(p))
  }
})

test_that('comparison plots preserve values and controls in both renderers', {
  f<-comparison_fixture();x<-compare_mfrm(f$ordinary,f$extension)
  grDevices::pdf(tempfile(fileext='.pdf'),width=9,height=5);on.exit(grDevices::dev.off())
  before<-par(c('mar','oma','mfrow','cex'))
  for(style in c('paired','difference')) {
    expect_silent(plot(x,style=style,palette='mono',show_labels=FALSE,title='Custom',caption=''))
    expect_equal(par(c('mar','oma','mfrow','cex')),before)
    raw<-plot(x,style=style,facet='Rater',draw=FALSE)
    expect_equal(raw$data$table$Y,if(style=='paired') c(-1,1) else c(0,0))
    expect_match(raw$data$alt_text,'descriptive')
    if(requireNamespace('ggplot2',quietly=TRUE)) expect_silent(ggplot2::ggplotGrob(as_ggplot(raw)))
  }
  expect_error(plot(x,facet='Unknown',draw=FALSE),'facet names')
})
