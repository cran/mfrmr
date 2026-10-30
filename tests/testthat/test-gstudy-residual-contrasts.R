residual_tasks <- function() {
  read.csv(system.file('extdata', 'mgenova-table12.csv', package = 'mfrmr'))
}

test_that('residual contrasts reproduce balanced ANOVA and independent linear-model limits', {
  cases <- list(residual_tasks(), mfrmr_example_gtheory_crossed, mfrmr_example_gtheory_nested)
  for (i in seq_along(cases)) {
    x <- cases[[i]]; scores <- if (i==1) c('V','W') else c('Content','Organization')
    g <- mfrm_g_study(x,scores,rater=if(i==1)NULL else 'Rater',
      nesting=if(i==3)c(Rater='Task') else NULL)
    original <- serialize(g,NULL)
    ci <- confint(g,parm='Residual',score=scores[1],assumption='normal',method='residual_contrast')
    ordinary <- confint(g,parm='Residual',score=scores[1],assumption='normal')
    expect_equal(ci[c('Estimate','Lower','Upper')],ordinary[c('Estimate','Lower','Upper')],tolerance=1e-9)
    expect_equal(ci$IntervalCenter,ci$Estimate,tolerance=1e-9)
    cov <- confint(g,score=scores,assumption='normal',method='residual_contrast')
    old <- confint(g,parm='Residual',score=scores,assumption='normal',method='regression_region')
    expect_equal(cov[cov$Source=='Residual',c('Lower','Upper')],old[c('Lower','Upper')],ignore_attr=TRUE,tolerance=1e-9)
    expect_equal(attr(ci,'residual_projection')$df,tail(g$degrees_of_freedom$DF,1))
    expect_identical(cov$Available,c(FALSE,rep(FALSE,nrow(cov)-2),TRUE))
    expect_true(all(is.na(cov$SE)))
    expect_identical(serialize(g,NULL),original)
  }
})

test_that('incomplete sources use a new residual center and an independently fitted fixed design', {
  x <- residual_tasks(); x <- x[(x$Person+x$Task)%%3!=0,]
  g <- mfrm_g_study(x,c('V','W'),rater=NULL,method='minque0')
  ci <- confint(g,score='V',method='residual_contrast',assumption='normal')
  j <- match('Residual',ci$Source)
  fit <- lm(V~factor(Person)+factor(Task),data=x)
  df <- df.residual(fit); ss <- sum(residuals(fit)^2)
  expect_equal(ci$IntervalCenter[j],ss/df,tolerance=1e-10)
  expect_equal(c(ci$Lower[j],ci$Upper[j]),ss/qchisq(c(.975,.025),df),tolerance=1e-10)
  expect_equal(ci$Estimate,unname(vapply(g$components,function(a)a[1,1],numeric(1))))
  expect_gt(abs(ci$Estimate[j]-ci$IntervalCenter[j]),1e-4)
  expect_true(all(grepl('only to Residual',ci$Status[-j])))
  composite <- confint(g,parm='Residual',weights=c(V=.4,W=-.6),
    method='residual_contrast',assumption='normal')
  x$Composite <- .4*x$V-.6*x$W
  fit <- lm(Composite~factor(Person)+factor(Task),data=x)
  expect_equal(composite$IntervalCenter,sum(residuals(fit)^2)/df,tolerance=1e-10)
  bad <- g; bad$data$V[1] <- bad$data$V[1]+1
  expect_error(confint(bad,method='residual_contrast',assumption='normal'),'disagree')
  bad <- g; bad$design$counts[1] <- 1
  expect_error(confint(bad,method='residual_contrast',assumption='normal'),'counts')
  expect_error(confint(g,method='residual_contrast',assumption='normal',draws=2000),'only')
  expect_error(confint(g,method='residual_contrast'),'normal')
  expect_error(confint(g,score=c('V','V'),method='residual_contrast',assumption='normal'),'distinct')
})

test_that('rank-zero and rank-one residual spaces do not masquerade as the same information', {
  x <- data.frame(Person=c(1,2,2,3,3,4,1),Task=c(1,1,2,2,3,3,2),
    V=c(2,1,4,3,6,5,8),W=c(1,4,2,7,3,5,9))
  for (n in c(6,7)) {
    g <- mfrm_g_study(x[seq_len(n),],c('V','W'),rater=NULL,method='minque0')
    call <- function(score)confint(g,parm='Residual',score=score,
      method='residual_contrast',assumption='normal')
    scalar <- call('V'); pair <- call(c('V','W'))
    expect_equal(attr(scalar,'residual_projection')$df,n-6)
    if (n==6) {
      expect_false(scalar$Available);expect_false(pair$Available)
      expect_match(scalar$Status,'No residual contrasts')
      expect_true(is.na(scalar$IntervalCenter));expect_true(is.finite(scalar$Estimate))
    } else {
      expect_true(scalar$Available);expect_true(pair$Available)
      expect_identical(attr(pair,'residual_projection')$region$method,'Marginal variances and Cauchy-Schwarz')
    }
  }
  x <- residual_tasks(); x$V <- x$Person+2*x$Task
  g <- mfrm_g_study(x,c('V','W'),rater=NULL)
  ci <- confint(g,parm='Residual',score='V',method='residual_contrast',assumption='normal')
  expect_false(ci$Available);expect_match(ci$Status,'numerical resolution')
})

test_that('incomplete crossed and nested residual spaces remove every nuisance interaction', {
  for (nested in c(FALSE, TRUE)) {
    x <- if (nested) mfrmr_example_gtheory_nested else mfrmr_example_gtheory_crossed
    x <- x[-c(1, 8, 23), ]
    g <- mfrm_g_study(x, c('Content', 'Organization'), method='minque0',
      nesting=if(nested)c(Rater='Task') else NULL)
    ci <- confint(g, parm='Residual', score='Content',
      method='residual_contrast', assumption='normal')
    x[c('Person', 'Rater', 'Task')] <- lapply(x[c('Person', 'Rater', 'Task')], factor)
    fit <- if (nested) lm(Content ~ Person*Task + Task:Rater, data=x) else
      lm(Content ~ (Person+Task+Rater)^2, data=x)
    ss <- sum(residuals(fit)^2); df <- df.residual(fit)
    expect_equal(attr(ci, 'residual_projection')$df, df)
    expect_equal(ci$IntervalCenter, ss/df, tolerance=1e-9)
    expect_equal(c(ci$Lower,ci$Upper), ss/qchisq(c(.975,.025),df), tolerance=1e-9)
    if (nested) {
      reverse <- mfrm_g_study(x,c('Content','Organization'),
        facets=c(Rater='Rater',Task='Task'),nesting=c(Rater='Task'),method='minque0')
      reversed <- confint(reverse,parm='Residual',score='Content',
        method='residual_contrast',assumption='normal')
      expect_equal(ci[c('Estimate','IntervalCenter','Lower','Upper')],
        reversed[c('Estimate','IntervalCenter','Lower','Upper')], tolerance=1e-9)
    }
  }
})

test_that('residual covariance intervals preserve signs units and row order', {
  x <- residual_tasks(); x <- x[(x$Person+x$Task)%%3!=0,]
  run <- function(x,score=c('V','W'))confint(
    mfrm_g_study(x,c('V','W'),rater=NULL,method='minque0'),
    parm='Residual',score=score,method='residual_contrast',assumption='normal')
  ci <- run(x); reverse <- run(x,c('W','V'))
  columns <- c('Estimate','IntervalCenter','Lower','Upper')
  expect_equal(ci[columns],reverse[columns],tolerance=1e-10)
  expect_equal(ci[columns],run(x[nrow(x):1,])[columns],tolerance=1e-10)
  x$V <- -.01*x$V+10; x$W <- 100*x$W
  changed <- run(x)
  expect_equal(changed$IntervalCenter,-ci$IntervalCenter,tolerance=1e-10)
  expect_equal(c(changed$Lower,changed$Upper),-c(ci$Upper,ci$Lower),tolerance=1e-10)
})

test_that('residual covariance does not report a rounded point as an interval', {
  g <- mfrm_g_study(residual_tasks(),c('V','W'),rater=NULL)
  call <- function(level)confint(g,parm='Residual',score=c('V','W'),
    method='residual_contrast',assumption='normal',level=level)
  ci <- call(1e-100); base <- call(.95)
  expect_identical(ci$Estimate,base$Estimate)
  expect_identical(ci$IntervalCenter,base$IntervalCenter)
  expect_false(ci$Available)
  expect_true(is.na(ci$Lower) && is.na(ci$Upper))
  expect_match(ci$Status,'numeric resolution')
  expect_match(attr(ci,'residual_projection')$region$method,'numeric resolution')
})

test_that('residual covariance endpoints outside score units are unavailable', {
  x <- data.frame(Person=c(1,2,2,3,3,4,1),Task=c(1,1,2,2,3,3,2),
    V=c(2,1,4,3,6,5,8),W=c(1,4,2,7,3,5,9))
  x[c('V','W')] <- x[c('V','W')]*1e150
  g <- mfrm_g_study(x,c('V','W'),rater=NULL,method='minque0')
  ci <- confint(g,parm='Residual',score=c('V','W'),
    method='residual_contrast',assumption='normal',level=1-1e-12)
  expect_false(ci$Available)
  expect_true(is.finite(ci$Estimate) && is.finite(ci$IntervalCenter))
  expect_true(is.na(ci$Lower) && is.na(ci$Upper))
  expect_match(ci$Status,'Unrepresentable residual interval')
})

test_that('residual intervals retain their interpretation in saved plot data', {
  x <- residual_tasks();x <- x[-c(1,8,23),]
  g <- mfrm_g_study(x,c('V','W'),rater=NULL,method='minque0')
  for (type in c('components','covariance')) {
    p <- plot(g,type=type,score=if(type=='components')'V' else c('V','W'),
      intervals=TRUE,method='residual_contrast',assumption='normal',draw=FALSE)
    tab <- plot_data(p,'table')
    expect_equal(sum(tab$IntervalAvailable),1)
    expect_true('IntervalCenter' %in% names(tab))
    expect_equal(nrow(plot_data(p,'unavailable_intervals')),nrow(tab)-1)
    file <- tempfile(fileext='.rds');withr::defer(unlink(file))
    saveRDS(p,file);saved <- readRDS(file)
    expect_identical(plot_data(saved,'table'),tab)
    if (requireNamespace('ggplot2',quietly=TRUE)) {
      pp <- as_ggplot(saved)
      expect_no_error(ggplot2::ggplot_build(pp))
      expect_null(pp$labels$title);expect_null(pp$labels$caption)
    }
  }
})
