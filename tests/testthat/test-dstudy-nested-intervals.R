nested_interval_pair <- function(custom = FALSE) {
  x <- mfrmr_example_gtheory_nested
  scores <- c('Content','Organization')
  means <- aggregate(x[scores], x[c('Person','Task')], mean)
  if (custom) {
    g <- mfrm_g_study(x, scores,
      facets=c(Occasion='Task',Assessor='Rater'),nesting=c(Assessor='Occasion'))
    a <- mfrm_g_study(means, scores, facets=c(Occasion='Task'))
    grid <- data.frame(Occasion=c(1,6,12,6),Assessor=3)
    reduced <- grid['Occasion']
  } else {
    g <- mfrm_g_study(x,scores,nesting=c(Rater='Task'))
    a <- mfrm_g_study(means,scores,rater=NULL)
    grid <- data.frame(Raters=3,Tasks=c(1,6,12,6))
    reduced <- grid['Tasks']
  }
  w <- cbind(Equal=c(Content=.5,Organization=.5),Difference=c(Content=1,Organization=-1))
  list(g=g,averaged=a,d=mfrm_d_study(g,grid,w),
    da=mfrm_d_study(a,reduced,w))
}

test_that('nested targets match independently averaged raw ratings', {
  for (custom in c(FALSE,TRUE)) {
    z <- nested_interval_pair(custom)
    original <- serialize(z$d,NULL)
    # Compare the separate raw-data ANOVA, not a second call to the reduction.
    expect_equal(unname(z$averaged$mean_products),
      unname(lapply(z$g$mean_products[c(1,2,4)],function(x)x/3)))
    expect_equal(z$averaged$degrees_of_freedom$DF,z$g$degrees_of_freedom$DF[c(1,2,4)])
    for (choice in list(list(score='Content'),list(composite='Equal'),list(composite='Difference'))) {
      for (metric in c('G','Phi','RelativeSEM','AbsoluteSEM')) for (reference in list(NULL,2L)) {
        args <- c(choice,list(parm=metric,reference=reference,assumption='normal',draws=4000,seed=41))
        ci <- do.call(confint,c(list(object=z$d),args))
        average <- do.call(confint,c(list(object=z$da),args))
        fields <- c('Estimate','Lower','Upper','Status','PointStatus','Method','Validation')
        expect_equal(as.data.frame(ci)[fields],as.data.frame(average)[fields],tolerance=1e-11)
        expect_identical(attr(ci,'source_summary'),z$g[c('design','components','mean_products','degrees_of_freedom')])
        expect_identical(attr(ci,'source_design')$nesting,z$g$design$nesting)
        expect_equal(attr(ci,'nesting_context')$source_children_per_parent,3)
        if (!is.null(reference)) {
          expect_true(all(ci$ReferenceCount==6))
          expect_equal(attr(ci,'reference')$children_per_parent,3)
          expect_identical(attr(ci,'reference')$count_facet,unname(z$g$design$nesting))
          identical_plan <- ci$Scenario==4
          expect_equal(ci$Lower[identical_plan],0);expect_equal(ci$Upper[identical_plan],0)
        }
        if (metric=='Phi') {
          expect_equal(attr(ci,'monte_carlo'),attr(average,'monte_carlo'),tolerance=1e-11)
          expect_identical(attr(ci,'pivot')$df,attr(average,'pivot')$df)
        }
      }
    }
    expect_identical(serialize(z$d,NULL),original)
  }
})

test_that('nested reduction retains source degrees of freedom and formula-based F limits', {
  z <- nested_interval_pair()
  ci <- confint(z$d,score='Content',assumption='normal')
  x <- mfrmr_example_gtheory_nested
  means <- aggregate(x['Content'],x[c('Person','Task')],mean)
  a <- summary(aov(Content~factor(Person)+factor(Task),data=means))[[1]]
  df <- a[c(1,3),'Df']; ratio <- a[1,'Mean Sq']/a[3,'Mean Sq']
  bounds <- pmax(1,ratio/qf(c(.975,.025),df[1],df[2]))
  m <- z$d$design_grid$Tasks; k <- length(unique(x$Task))
  expect_equal(ci$Lower,m*(bounds[1]-1)/(m*(bounds[1]-1)+k))
  expect_equal(ci$Upper,m*(bounds[2]-1)/(m*(bounds[2]-1)+k))
  expect_equal(attr(ci,'pivot')$df,df)
  expect_identical(ci$MethodCode,rep('exact_f',4))
  expect_identical(ci,confint(z$d,score='Content',assumption='normal',method='exact_f'))
  expect_s3_class(confint(z$d,score='Content',assumption='normal',method='approximate_f'),'mfrm_dstudy_intervals')
  changed <- mfrm_d_study(z$g,data.frame(Raters=c(2,3),Tasks=6))
  expect_error(confint(changed,assumption='normal'),'child count per parent')
  expect_error(confint(changed,parm='Phi',reference=1,assumption='normal'),'child count per parent')
  expect_error(confint(changed,parm='RelativeSEM',reference=1,assumption='normal'),'same future child count')
})

test_that('nested SEM differences need common future child counts, not unchanged ones', {
  z <- nested_interval_pair()
  for (children in c(1,3,6)) for (metric in c('RelativeSEM','AbsoluteSEM')) {
    d <- mfrm_d_study(z$g,data.frame(Raters=children,Tasks=c(3,6,12)))
    individual <- confint(d,parm=metric,assumption='normal')
    ci <- confint(d,parm=metric,reference=2,assumption='normal')
    factor <- sqrt(6/c(3,12))-1
    expected <- t(vapply(factor,function(k)sort(k*c(individual$Lower[2],individual$Upper[2])),numeric(2)))
    expect_equal(ci$Lower,expected[,1]);expect_equal(ci$Upper,expected[,2])
    expect_equal(ci$Estimate,individual$Estimate[c(1,3)]-individual$Estimate[2])
    expect_equal(attr(ci,'reference')$children_per_parent,children)
    expect_equal(attr(ci,'reference')$count,6)
    expect_identical(ci$Validation,individual$Validation[c(1,3)])
    if (children==3 && metric=='RelativeSEM') expect_match(ci$Method[1],'Chi-square')
    else expect_match(ci$Method[1],'MLS')
  }
})

test_that('nested interval plots keep both counts and validate saved source identity', {
  z <- nested_interval_pair(TRUE)
  ci <- confint(z$d,parm='Phi',composite='Equal',assumption='normal',draws=4000)
  figure <- plot(z$d,x_var='Occasion',metric='Phi',composite='Equal',intervals=ci,draw=FALSE)
  expect_equal(figure$data$series$Lower,ci$Lower)
  expect_identical(figure$data$intervals[[1]],ci)
  change <- confint(z$d,composite='Equal',reference=2,assumption='normal')
  forest <- plot(change,draw=FALSE)
  expect_true(all(grepl('Assessor per Occasion',forest$data$table$Plan)))
  expect_true(any(grepl('Assessor per Occasion',forest$data$notes$Text)))
  path <- tempfile(fileext='.rds');on.exit(unlink(path))
  saveRDS(list(d=z$d,ci=ci,figure=figure,change=change),path)
  saved <- readRDS(path)
  expect_identical(saved$figure,figure)
  z$d$gstudy$data <- NULL
  expect_identical(ci,confint(z$d,parm='Phi',composite='Equal',assumption='normal',draws=4000))
  if (requireNamespace('ggplot2',quietly=TRUE)) {
    gg <- as_ggplot(saved$figure)
    expect_identical(plot_data(gg,'table'),figure$data$table)
    expect_s3_class(ggplot2::ggplotGrob(gg),'gtable')
    expect_s3_class(ggplot2::ggplotGrob(as_ggplot(forest)),'gtable')
  }
})

test_that('nested degenerate source intervals are not replaced by zero differences', {
  x <- mfrmr_example_gtheory_nested
  x$Content <- x$Organization <- 1
  g <- mfrm_g_study(x,c('Content','Organization'),nesting=c(Rater='Task'))
  d <- mfrm_d_study(g,data.frame(Raters=3,Tasks=c(6,6,12)))
  for (metric in c('G','Phi','RelativeSEM','AbsoluteSEM')) {
    ci <- confint(d,parm=metric,assumption='normal',draws=4000)
    delta <- confint(d,parm=metric,reference=1,assumption='normal',draws=4000)
    expect_false(any(ci$Available));expect_false(any(delta$Available))
    expect_true(all(is.na(delta$Lower) & is.na(delta$Upper)))
    if (metric=='Phi') expect_equal(attr(delta,'monte_carlo')$draws_used,0L)
  }
  # Incomplete-source estimators do not inherit the balanced pivot, even
  # when their observed count labels happen to match future counts.
  incomplete <- mfrmr_example_gtheory_nested[-1,]
  g <- mfrm_g_study(incomplete,c('Content','Organization'),
    nesting=c(Rater='Task'),method='minque0')
  d <- mfrm_d_study(g,data.frame(Raters=3,Tasks=6))
  expect_error(confint(d,assumption='normal'),'MINQUE')
})
