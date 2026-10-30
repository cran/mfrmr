residual_sem_data <- function() {
  read.csv(system.file('extdata','mgenova-table12.csv',package='mfrmr'))
}

test_that("residual SEM retains failure of the common square-root image", {
  g <- mfrm_g_study(residual_sem_data(), "V", rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = c(1, 1)))
  local_mocked_bindings(.mfrm_mvgt_residual_intervals = function(...) {
    out <- data.frame(Source = "Residual", Estimate = 1, IntervalCenter = 1,
      Lower = 1, Upper = 1 + .Machine$double.eps, Available = TRUE,
      Status = "Exact normal-model interval", Validation = "Exact normal-model interval")
    attr(out, "residual_projection") <- list(df = 20)
    out
  })
  ci <- confint(d, parm = "RelativeSEM", method = "residual_contrast", assumption = "normal")
  expect_false(any(ci$Available))
  expect_true(all(ci$Status == "SEM interval below numeric resolution"))
  expect_true(attr(ci, "base_interval")$Available)
  identical_plan <- confint(d, parm = "RelativeSEM", method = "residual_contrast",
    reference = 1, assumption = "normal")
  expect_false(identical_plan$Available)
  expect_identical(identical_plan$Status, "SEM interval below numeric resolution")
})

test_that('one-facet residual SEM agrees with independent linear-model intervals', {
  for (incomplete in c(FALSE,TRUE)) {
    x <- residual_sem_data()
    if (incomplete) x <- x[-c(1,8,23),]
    x$Selected <- .4*x$V-.6*x$W
    g <- mfrm_g_study(x,c('V','W'),rater=NULL,
      method=if(incomplete)'minque0' else 'anova')
    weights <- cbind(Difference=c(V=.4,W=-.6),Mean=c(V=.5,W=.5))
    d <- mfrm_d_study(g,data.frame(Tasks=c(1,3,6,12)),weights)
    original <- serialize(d,NULL)
    ci <- confint(d,parm='RelativeSEM',composite='Difference',method='residual_contrast',assumption='normal')
    fit <- lm(Selected~factor(Person)+factor(Task),data=x)
    df <- df.residual(fit); ss <- sum(residuals(fit)^2)
    expected <- sqrt(ss/qchisq(c(.975,.025),df))
    expect_equal(ci$Lower,expected[1]/sqrt(ci$Tasks),tolerance=1e-10)
    expect_equal(ci$Upper,expected[2]/sqrt(ci$Tasks),tolerance=1e-10)
    expect_equal(ci$IntervalCenter,sqrt(ss/df)/sqrt(ci$Tasks),tolerance=1e-10)
    expect_identical(ci$Estimate,subset(d$coefficients,Kind=='Composite' & Score=='Difference')$RelativeSEM)
    expect_identical(attr(ci,'source_data'),g$data)
    expect_equal(attr(ci,'residual_projection')$df,df)
    expect_identical(serialize(d,NULL),original)
    expect_true(all(ci$Available))
    if (!incomplete) {
      previous <- confint(d,parm='RelativeSEM',composite='Difference',assumption='normal')
      expect_equal(ci[c('Estimate','Lower','Upper')],previous[c('Estimate','Lower','Upper')],ignore_attr=TRUE)
      comparison <- confint(d,parm='RelativeSEM',composite='Difference',
        method=c('mls','residual_contrast'),assumption='normal')
      expect_identical(attr(comparison,'method_results')$residual_contrast,ci)
      expect_identical(attr(comparison,'method_results')$mls,previous)
      expect_no_error(summary(d,intervals=comparison))
      expect_no_error(plot(comparison,draw=FALSE))
    }
  }
})

test_that('paired residual SEM limits map one common event with signs and equal counts', {
  x <- residual_sem_data(); x <- x[-c(1,8,23),]
  g <- mfrm_g_study(x,c('V','W'),rater=NULL,method='minque0')
  d <- mfrm_d_study(g,data.frame(Tasks=c(6,12,6,3)),c(V=.4,W=-.6))
  individual <- confint(d,parm='RelativeSEM',method='residual_contrast',assumption='normal')
  ci <- confint(d,parm='RelativeSEM',reference=1,method='residual_contrast',assumption='normal')
  multiplier <- 1/sqrt(ci$Tasks)-1/sqrt(6)
  fit <- lm(I(.4*V-.6*W)~factor(Person)+factor(Task),data=x)
  sd_bounds <- sqrt(sum(residuals(fit)^2)/qchisq(c(.975,.025),df.residual(fit)))
  expected <- t(vapply(multiplier,function(a)sort(a*sd_bounds),numeric(2)))
  expect_equal(ci$Lower,expected[,1]);expect_equal(ci$Upper,expected[,2])
  expect_equal(ci$Estimate,individual$Estimate[-1]-individual$Estimate[1])
  expect_equal(ci$IntervalCenter,individual$IntervalCenter[-1]-individual$IntervalCenter[1])
  expect_identical(ci$Status[2],'Identical counts')
  expect_equal(c(ci$Lower[2],ci$Upper[2]),c(0,0))
  reverse <- confint(d,parm='RelativeSEM',reference=2,method='residual_contrast',assumption='normal')
  expect_equal(c(reverse$Lower[1],reverse$Upper[1]),-c(ci$Upper[1],ci$Lower[1]))
  expect_no_error(summary(d,intervals=list(individual,ci)))
  expect_true(all(vapply(attr(ci,'diagnostics'),function(z)is.finite(z$multiplier),logical(1))))
  # These intervals are an image of one SD event, not subtracted marginal limits.
  expect_gt(abs(ci$Lower[1]-(individual$Lower[2]-individual$Upper[1])),1e-3)
  large <- mfrm_d_study(g,data.frame(Tasks=c(1e9,1e9+1)))
  delta <- confint(large,parm='RelativeSEM',reference=1,method='residual_contrast',assumption='normal')
  expect_lt(delta$Upper,0);expect_gt(delta$Upper,delta$Lower)
})

test_that('negative original estimates and unavailable projections retain distinct states', {
  x <- residual_sem_data();x <- x[-c(1,8,23),]
  x$V <- x$Person-x$Task+.01*sin(seq_len(nrow(x)))
  g <- mfrm_g_study(x,'V',rater=NULL,method='minque0')
  expect_lt(g$components$Residual[1],0)
  d <- mfrm_d_study(g,data.frame(Tasks=c(3,6)))
  ci <- confint(d,parm='RelativeSEM',method='residual_contrast',assumption='normal')
  expect_true(all(is.na(ci$Estimate)));expect_true(all(ci$Available))
  expect_true(all(is.finite(ci$IntervalCenter)));expect_true(all(ci$NegativeComponent))
  expect_no_error(summary(d,intervals=ci))
  expect_no_error(plot(d,type='sem',metric='RelativeSEM',intervals=ci,draw=FALSE))
  change <- confint(d,parm='RelativeSEM',reference=1,method='residual_contrast',assumption='normal')
  expect_true(is.na(change$Estimate));expect_true(change$Available)
  expect_match(change$PointStatus,'unavailable')
  # A connected but saturated incidence model has no residual information.
  x <- data.frame(Person=c(1,2,2,3,3,4),Task=c(1,1,2,2,3,3),V=c(2,1,4,3,6,5))
  d <- mfrm_d_study(mfrm_g_study(x,'V',rater=NULL,method='minque0'),
    data.frame(Tasks=c(3,6,3)))
  ci <- confint(d,parm='RelativeSEM',method='residual_contrast',assumption='normal')
  expect_false(any(ci$Available));expect_true(all(grepl('No residual contrasts',ci$Status)))
  change <- confint(d,parm='RelativeSEM',reference=1,method='residual_contrast',assumption='normal')
  expect_false(any(change$Available));expect_true(all(is.na(change$Lower)))
  x <- residual_sem_data();x$V <- x$Person+2*x$Task
  d <- mfrm_d_study(mfrm_g_study(x,'V',rater=NULL))
  ci <- confint(d,parm='RelativeSEM',method='residual_contrast',assumption='normal')
  expect_false(any(ci$Available));expect_true(all(grepl('numerical resolution',ci$Status)))
})

test_that('residual SEM supports named facet roles and scale changes without broadening scope', {
  x <- residual_sem_data();x <- x[-c(1,8,23),]
  run <- function(w) {
    g <- mfrm_g_study(x,c('V','W'),facets=c(Assessor='Task'),method='minque0')
    d <- mfrm_d_study(g,data.frame(Assessor=c(2,4)),weights=w)
    confint(d,parm='RelativeSEM',method='residual_contrast',assumption='normal')
  }
  ci <- run(c(V=.4,W=-.6))
  for (scale in c(-3,1e-80,1e80)) {
    other <- run(scale*c(W=-.6,V=.4))
    for (column in c('Estimate','IntervalCenter','Lower','Upper'))
      expect_equal(other[[column]]/abs(scale),ci[[column]],tolerance=1e-9)
  }
  d <- mfrm_d_study(mfrm_g_study(x,c('V','W'),rater=NULL,method='minque0'),
    data.frame(Tasks=c(3,6)))
  for (metric in c('G','Phi','AbsoluteSEM'))
    expect_error(confint(d,parm=metric,method='residual_contrast',assumption='normal'),'Unsupported')
  crossed <- mfrm_d_study(mfrm_g_study(mfrmr_example_gtheory_crossed,'Content'))
  expect_error(confint(crossed,parm='RelativeSEM',method='residual_contrast',assumption='normal'),'Unsupported')
  expect_error(confint(d,parm='RelativeSEM',method='residual_contrast'),'normal')
})

test_that('summaries ribbons and saved results preserve residual SEM source and center', {
  x <- residual_sem_data();x <- x[-c(1,8,23),]
  d <- mfrm_d_study(mfrm_g_study(x,c('V','W'),rater=NULL,method='minque0'),
    data.frame(Tasks=c(2,4,6)))
  s <- summary(d,intervals=TRUE,parm='RelativeSEM',method='residual_contrast',assumption='normal')
  ci <- s$intervals$RelativeSEM
  expect_identical(s$table$IntervalCenter,ci$IntervalCenter)
  p <- plot(d,type='sem',metric='RelativeSEM',intervals=TRUE,method='residual_contrast',assumption='normal',draw=FALSE)
  expect_equal(plot_data(p,'series')$Lower,ci$Lower)
  expect_equal(plot_data(p,'series')$Upper,ci$Upper)
  expect_identical(plot_data(p,'intervals')[[1]],ci)
  changed <- d;changed$gstudy$data$V[1] <- changed$gstudy$data$V[1]+1
  expect_error(summary(changed,intervals=ci),'source ratings')
  expect_error(plot(changed,type='sem',metric='RelativeSEM',intervals=ci,draw=FALSE),'source ratings')
  changed <- ci;attr(changed,'source_data') <- NULL
  expect_error(summary(d,intervals=changed),'source ratings')
  file <- tempfile(fileext='.rds');withr::defer(unlink(file));saveRDS(list(s=s,p=p,ci=ci),file)
  local_mocked_bindings(.mfrm_mvgt_residual_block=function(...)stop('Unexpected recalculation'),.package='mfrmr')
  saved <- readRDS(file)
  expect_identical(saved$s,s);expect_identical(saved$ci,ci)
  expect_no_error(summary(d,intervals=saved$ci));expect_output(print(saved$s),'IntervalCenter')
  if (requireNamespace('ggplot2',quietly=TRUE)) {
    gg <- as_ggplot(saved$p)
    expect_no_error(ggplot2::ggplot_build(gg))
    expect_null(gg$labels$title);expect_null(gg$labels$caption)
  }
})

test_that('method comparisons retain different residual and REML inference bases', {
  skip_if_not_installed('lme4');skip_if_not_installed('minqa')
  g <- mfrm_g_study(residual_sem_data(),'V',rater=NULL)
  d <- mfrm_d_study(g,data.frame(Tasks=c(3,6)))
  ci <- confint(d,parm='RelativeSEM',method=c('reml_profile','residual_contrast'),assumption='normal')
  parts <- attr(ci,'method_results')
  expect_identical(names(parts),c('reml_profile','residual_contrast'))
  expect_true(is.list(attr(parts$reml_profile,'reml_source')))
  expect_identical(attr(parts$residual_contrast,'base_interval')$Source,'Residual')
  expect_equal(ci$Estimate,rep(d$coefficients$RelativeSEM,2))
  expect_equal(ci$IntervalCenter,c(parts$reml_profile$IntervalCenter,parts$residual_contrast$IntervalCenter))
  expect_no_error(summary(d,intervals=ci))
  p <- plot(ci,draw=FALSE)
  expect_equal(plot_data(p,'table')$Estimate,ci$Estimate)
})
