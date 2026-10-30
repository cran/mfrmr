sem_tasks <- function() {
  read.csv(system.file('extdata', 'mgenova-table12.csv', package = 'mfrmr', mustWork = TRUE))
}

test_that('one-facet relative SEM intervals equal independent residual chi-square limits', {
  x <- sem_tasks(); x$Total <- (x$V+x$W)/2
  g <- mfrm_g_study(x, c('V','W'), rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks=c(1,6,12,60)), weights=c(V=.5,W=.5))
  aov_table <- summary(aov(Total ~ factor(Person)+factor(Task), data=x))[[1]]
  ms <- aov_table[3,'Mean Sq']; df <- aov_table[3,'Df']
  ci <- confint(d, parm='RelativeSEM', assumption='normal')
  expect_equal(ci$Lower, sqrt(df*ms/qchisq(.975,df)/d$design_grid$Tasks))
  expect_equal(ci$Upper, sqrt(df*ms/qchisq(.025,df)/d$design_grid$Tasks))
  expect_identical(ci$Estimate, subset(d$coefficients,Kind=='Composite')$RelativeSEM)
  expect_equal(ci$ErrorVarianceSE^2, 2*ms^2/(df+2)/d$design_grid$Tasks^2)
  expect_true(all(ci$Method=='Chi-square SEM interval'))
  expect_true(all(grepl('Exact normal-model',ci$Validation)))
  expect_identical(ci,confint(d,parm='RelativeSEM',assumption='normal',method='mls'))
  wider<-confint(d,parm='RelativeSEM',level=.99,assumption='normal')
  expect_true(all(wider$Lower<=ci$Lower & wider$Upper>=ci$Upper))
  expect_equal(attr(ci,'mean_squares')$df[3],df)
  absolute<-confint(d,parm='AbsoluteSEM',assumption='normal')
  ems<-aov_table[c(2,3),'Mean Sq']; edf<-aov_table[c(2,3),'Df']
  # Positive sum: independent Graybill-Wang formula, with no component-CI addition.
  z<-c(1/10,9/10)*ems
  lo<-sum(z)-sqrt(sum(((1-edf/qchisq(.975,edf))*z)^2))
  hi<-sum(z)+sqrt(sum(((edf/qchisq(.025,edf)-1)*z)^2))
  expect_equal(absolute$Lower,sqrt(lo/d$design_grid$Tasks))
  expect_equal(absolute$Upper,sqrt(hi/d$design_grid$Tasks))
  expect_true(all(absolute$Method=='Graybill-Wang MLS SEM interval'))
  expect_true(all(grepl('approximation',absolute$Validation)))
})

test_that("chi-square intervals retain representable scale and identify collapsed widths", {
  # df=2 has an independent exponential-quantile expression.
  exact <- 3 / c(-log(.025), -log1p(-.025))
  result <- .mfrm_mvgt_linear_variance_interval(3, 2, 1, .95)
  expect_equal(result$bounds, exact, tolerance = 1e-14)
  base <- .mfrm_mvgt_linear_variance_interval(1, 10000, 1, .95)
  large <- .mfrm_mvgt_linear_variance_interval(1e308, 10000, 1, .95)
  expect_identical(large$status, "Available")
  expect_equal(large$bounds / 1e308, base$bounds, tolerance = 1e-14)
  narrow <- .mfrm_mvgt_linear_variance_interval(1, 20, 1, 1e-18)
  expect_identical(narrow$status, "Variance interval below numeric resolution")
  expect_true(all(is.na(narrow$bounds)))
  expect_true(narrow$se > 0)
  tiny <- .mfrm_mvgt_linear_variance_interval(.Machine$double.xmin * .Machine$double.eps, 2, 1, .95)
  expect_identical(tiny$status, "Unrepresentable chi-square endpoints")
  expect_true(all(is.na(tiny$bounds)))
})

test_that("SEM and its differences retain a failed variance-interval diagnosis", {
  g <- mfrm_g_study(sem_tasks(), "V", rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = c(6, 12, 6)))
  for (method in c("mls", "residual_contrast")) {
    ci <- confint(d, parm = "RelativeSEM", method = method, level = 1e-18, assumption = "normal")
    expect_true(all(is.finite(ci$Estimate)))
    expect_false(any(ci$Available))
    expect_true(all(grepl("below numeric resolution", ci$Status)))
    change <- confint(d, parm = "RelativeSEM", method = method, reference = 1,
      level = 1e-18, assumption = "normal")
    expect_false(any(change$Available))
    expect_true(all(grepl("below numeric resolution", change$Status)))
    expect_identical(change$Estimate[2], 0)
  }
})

test_that("square-root rounding does not turn a positive variance width into exact SEM", {
  g <- mfrm_g_study(sem_tasks(), "V", rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = 6))
  local_mocked_bindings(.mfrm_mvgt_linear_variance_interval = function(...) {
    list(estimate = 1, se = .1, raw = c(1, 1 + .Machine$double.eps),
      bounds = c(1, 1 + .Machine$double.eps), status = "Available", method = "Chi-square")
  })
  ci <- confint(d, parm = "RelativeSEM", assumption = "normal")
  expect_false(ci$Available)
  expect_identical(ci$Status, "SEM interval below numeric resolution")
  expect_true(is.finite(ci$Estimate))
  expect_true(diff(attr(ci, "diagnostics")[[1]]$normalized_raw_bounds) > 0)
})

test_that('linear-variance helper reproduces the published sum and edge cases', {
  # Hess & Iyer (2002), p.4, all four components summed. Printed MS values.
  ms<-c(.645160,.417296,.083924,.083860);df<-c(1,4,24,60)
  L<-rbind(c(45,15,3,1),c(0,15,3,1),c(0,0,3,1),c(0,0,0,1))
  a<-drop(rep(1,4)%*%solve(L))
  ci<-.mfrm_mvgt_linear_variance_interval(ms,df,a,.95)
  expect_lt(abs(ci$estimate-.1111699),5e-7) # Printed mean squares have six decimals.
  expect_equal(ci$raw,c(.0867974,14.69615),tolerance=1e-6)
  expect_identical(ci$method,'Graybill-Wang MLS')
  for (scale in c(1e-160,1e160)) {
    scaled<-.mfrm_mvgt_linear_variance_interval(ms*scale,df,a,.95)
    expect_equal(scaled$raw/scale,ci$raw)
    expect_equal(scaled$se/scale,ci$se)
  }
  empty<-.mfrm_mvgt_linear_variance_interval(c(0,1),c(9,45),c(1,-1),.95)
  expect_identical(empty$status,'Empty nonnegative interval')
  expect_true(all(is.na(empty$bounds)));expect_lt(empty$raw[2],0)
  zero<-.mfrm_mvgt_linear_variance_interval(c(0,0),c(9,45),c(1,-1),.95)
  expect_identical(zero$status,'Degenerate mean squares')
  expect_equal(zero$se,0);expect_true(all(is.na(zero$bounds)))
})

test_that('crossed and nested SEM intervals preserve their error maps and units', {
  for (nested in c(FALSE,TRUE)) {
    x <- if(nested) mfrmr_example_gtheory_nested else mfrmr_example_gtheory_crossed
    g <- mfrm_g_study(x,c('Content','Organization'),
      nesting=if(nested)c(Rater='Task') else NULL)
    if(nested) g <- mfrm_g_study(x,c('Content','Organization'),
      facets=c(Occasion='Task',Assessor='Rater'),nesting=c(Assessor='Occasion'))
    grid <- if(nested)data.frame(Occasion=c(3,6,12),Assessor=c(1,3,6)) else
      data.frame(Raters=c(1,2,4),Tasks=c(3,6,12))
    w<-c(Content=.3,Organization=-.7)
    d<-mfrm_d_study(g,grid,w)
    original<-serialize(d,NULL)
    for(metric in c('RelativeSEM','AbsoluteSEM')) {
      ci<-confint(d,parm=metric,assumption='normal')
      expect_identical(ci$Estimate,subset(d$coefficients,Kind=='Composite')[[metric]])
      for(i in seq_len(nrow(ci))) {
        diag<-attr(ci,'diagnostics')[[i]];ms<-attr(ci,'mean_squares')
        # Raw projected components and the separate MS map give the same error.
        theta<-vapply(g$components,function(a)drop(crossprod(w,a%*%w)),numeric(1))
        expect_equal(sum(diag$component_weights*theta),ci$ErrorVariance[i])
        expect_equal(unname(diag$mean_square_coefficients),
          drop(diag$component_weights %*% ms$component_map), ignore_attr=TRUE, tolerance=1e-14)
        expect_equal(sum(diag$mean_square_coefficients*ms$normalized)*ms$scale,ci$ErrorVariance[i])
        C<-ms$component_map%*%diag(2*(ms$normalized*ms$scale)^2/(ms$df+2))%*%t(ms$component_map)
        expect_equal(ci$ErrorVarianceSE[i]^2,drop(crossprod(diag$component_weights,C%*%diag$component_weights)))
      }
      scaled<-confint(mfrm_d_study(g,grid,-3*w),parm=metric,assumption='normal')
      expect_equal(scaled$Lower,3*ci$Lower);expect_equal(scaled$Upper,3*ci$Upper)
      expect_equal(scaled$ErrorVarianceSE,9*ci$ErrorVarianceSE)
      expect_identical(serialize(d,NULL),original)
    }
    if(nested) {
      ci<-confint(d,parm='RelativeSEM',assumption='normal')
      expect_identical(ci$Method[2],'Chi-square SEM interval')
      df<-g$degrees_of_freedom$DF[4]
      ms<-drop(crossprod(w,g$mean_products[[4]]%*%w))
      expect_equal(c(ci$Lower[2],ci$Upper[2]),sqrt(ms*df/qchisq(c(.975,.025),df)/(3*6)))
      expect_error(confint(d,assumption='normal'),'nested')
    }
    expect_error(confint(d,parm='RelativeSEM',reference=1,assumption='normal'),'reference')
  }
})

test_that('exact count cancellations do not accidentally choose signed MLS', {
  x<-expand.grid(Person=1:4,Rater=1:3,Task=1:6)
  x$Score<-sin(seq_len(nrow(x)))+cos(seq_len(nrow(x))/3)
  d<-mfrm_d_study(mfrm_g_study(x,'Score'),
    data.frame(Raters=c(1,1),Tasks=c(4,5)))
  for(metric in c('RelativeSEM','AbsoluteSEM')) {
    ci<-confint(d,parm=metric,assumption='normal')
    expect_identical(attr(ci,'diagnostics')[[1]]$mean_square_coefficients[['Residual']],0)
    expect_identical(ci$Method[1],'Graybill-Wang MLS SEM interval')
    expect_identical(ci$Method[2],'Ting MLS SEM interval')
  }
})

test_that('one-facet SEM differences map the same interval and reverse correctly', {
  g<-mfrm_g_study(sem_tasks(),c('V','W'),rater=NULL)
  d<-mfrm_d_study(g,data.frame(Tasks=c(6,12,6,3)),weights=c(V=.5,W=.5))
  for(metric in c('RelativeSEM','AbsoluteSEM')) {
    ci<-confint(d,parm=metric,assumption='normal')
    difference<-confint(d,parm=metric,reference=1,assumption='normal')
    k<-sqrt(6/d$design_grid$Tasks[-1])-1
    expected<-t(vapply(k,function(z)sort(z*c(ci$Lower[1],ci$Upper[1])),numeric(2)))
    expect_equal(difference$Lower,expected[,1]);expect_equal(difference$Upper,expected[,2])
    expect_equal(difference$Estimate,ci$Estimate[-1]-ci$Estimate[1])
    expect_identical(difference$Status[2],'Identical counts')
    expect_equal(c(difference$Lower[2],difference$Upper[2]),c(0,0))
    reverse<-confint(d,parm=metric,reference=2,assumption='normal')
    expect_equal(c(reverse$Lower[1],reverse$Upper[1]),-rev(c(difference$Lower[1],difference$Upper[1])))
    expect_equal(attr(difference,'reference')$sem_interval,c(ci$Lower[1],ci$Upper[1]))
    expect_false('ErrorVariance'%in%names(difference))
    figure<-plot(difference,draw=FALSE)
    expect_true(figure$data$xrange[1]<0);expect_true(figure$data$xrange[2]>0)
  }
})

test_that('SEM intervals keep degenerate states and save without raw ratings', {
  x<-sem_tasks();x$V<-5;x$W<-5
  d<-mfrm_d_study(mfrm_g_study(x,c('V','W'),rater=NULL),data.frame(Tasks=c(6,12)))
  for(metric in c('RelativeSEM','AbsoluteSEM')) {
    ci<-confint(d,parm=metric,assumption='normal')
    expect_equal(ci$Estimate,c(0,0));expect_false(any(ci$Available))
    expect_true(all(ci$Status=='Degenerate mean squares'))
    delta<-confint(d,parm=metric,reference=1,assumption='normal')
    expect_true(is.na(delta$Lower)&&is.na(delta$Upper))
  }
  g<-mfrm_g_study(sem_tasks(),c('V','W'),rater=NULL)
  d<-mfrm_d_study(g)
  ci<-confint(d,parm='AbsoluteSEM',assumption='normal')
  path<-tempfile(fileext='.rds');on.exit(unlink(path),add=TRUE)
  d$gstudy$data<-NULL;saveRDS(d,path)
  expect_identical(confint(readRDS(path),parm='AbsoluteSEM',assumption='normal'),ci)
  expect_error(confint(d,parm='AbsoluteSEM',assumption='normal',method='exact_f'),'Unsupported')
  expect_error(confint(d,parm='RelativeSEM'),'assumption')
  bad<-d;bad$gstudy$design$method<-'Multivariate MINQUE(0)'
  expect_error(confint(bad,parm='RelativeSEM',assumption='normal'),'MINQUE')
})

test_that('SEM ribbons and forests retain score units, endpoints and saved results', {
  skip_if_not_installed('ggplot2',minimum_version='3.4.0')
  g<-mfrm_g_study(sem_tasks(),c('V','W'),rater=NULL)
  d<-mfrm_d_study(g,data.frame(Tasks=c(1,3,6,12)))
  ci<-confint(d,parm='AbsoluteSEM',score='V',assumption='normal')
  grDevices::pdf(NULL,width=7,height=5);on.exit(grDevices::dev.off(),add=TRUE)
  expect_no_warning(p<-plot(d,type='sem',metric='AbsoluteSEM',score='V',intervals=TRUE,assumption='normal'))
  expect_equal(p$data$series$Lower,ci$Lower);expect_equal(p$data$series$Upper,ci$Upper)
  expect_identical(plot(d,type='sem',metric='AbsoluteSEM',score='V',intervals=ci,draw=FALSE),p)
  forest<-plot(ci,draw=FALSE)
  expect_gte(forest$data$xrange[2],max(ci$Upper));expect_match(forest$data$xlab,'score units')
  expect_no_warning(plot(ci))
  expect_no_warning(ggplot2::ggplotGrob(as_ggplot(forest)))
  saved<-tempfile(fileext='.rds');on.exit(unlink(saved),add=TRUE)
  saveRDS(p,saved)
  local_mocked_bindings(confint.mfrm_multivariate_d_study=function(...)stop('No interval recomputation'),
    .package='mfrmr')
  gg<-as_ggplot(readRDS(saved))
  expect_identical(plot_data(gg,'intervals'),p$data$intervals)
  expect_gte(gg$scales$get_scales('y')$limits[2],max(ci$Upper))
  expect_no_warning(ggplot2::ggplotGrob(gg))
})
