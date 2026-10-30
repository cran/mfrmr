test_that('general nested F uses the population target and independent EMS algebra', {
  n<-10;r<-2;t<-3
  L<-rbind(c(r*t,0,0,r,1),c(0,n*r,n,r,1),c(0,0,n,0,1),c(0,0,0,r,1),c(0,0,0,0,1))
  df<-c(n-1,t-1,t*(r-1),(n-1)*(t-1),t*(n-1)*(r-1))
  map<-list(A=solve(L),df=df)
  for (theta in list(c(.6,.2,.3,.1,.8),c(0,0,0,0,1))) {
    ms<-drop(L%*%theta)
    for (future in list(c(1,2),c(2,6),c(4,12))) for (absolute in c(FALSE,TRUE)) {
      R<-future[1];T<-future[2];mi<-r/R;mj<-t/T
      e<-c(0,0,0,1/T,1/(R*T))+as.integer(absolute)*c(0,1/T,1/(R*T),0,0)
      rho<-theta[1]/(theta[1]+sum(e*theta))
      weights<-function(rho) if(absolute)c(rho*mj/n,rho*mj*(mi-1)/n,
        1-rho+rho*mj*(1-1/n),rho*mj*(mi-1)*(1-1/n)) else
          c(0,0,1-rho+rho*mj,rho*mj*(mi-1))
      den<-sum(weights(rho)*ms[-1])
      expect_equal(den,(1-rho)*ms[1])
      ans<-.mfrm_mvgt_approximate_f(ms,map,e,.95)
      expect_equal(ans$unrestricted_estimate,rho,tolerance=1e-12)
      nu<-den^2/sum((weights(rho)*ms[-1])^2/df[-1])
      expect_equal(ans$denominator_df,nu)
      grid<-seq(0,1,length.out=2001)
      denominator<-vapply(grid,function(q)sum(weights(q)*ms[-1]),numeric(1))
      statistic<-(1-grid)*ms[1]/denominator
      accepted<-denominator>0 & statistic>=qf(.025,df[1],nu) & statistic<=qf(.975,df[1],nu)
      returned<-if(all(is.finite(ans$bounds)))grid>=ans$bounds[1] & grid<=ans$bounds[2] else rep(FALSE,length(grid))
      expect_identical(returned,accepted)
      if(!absolute && R==r) {
        exact<-.mfrm_mvgt_g_f_interval(ms[1],ms[4],df[c(1,4)],t,T,.95)
        expect_equal(ans$bounds,c(exact$lower,exact$upper))
      }
      if(absolute && rho>0) {
        # The printed p.278 c,d terms omit -rho*mj/n and +rho*mj/n.
        printed<-weights(rho)+c(0,0,rho*mj/n,-rho*mj/n)
        expect_gt(abs(sum(printed*ms[-1])-den),1e-4)
      }
    }
  }
})

test_that('region projections match independently enumerated nested polytope vertices', {
  skip_if_not_installed('lpSolve')
  n<-10;r<-2;t<-3
  L<-rbind(c(r*t,0,0,r,1),c(0,n*r,n,r,1),c(0,0,n,0,1),c(0,0,0,r,1),c(0,0,0,0,1))
  map<-list(A=solve(L),df=c(9,2,3,18,27))
  ms<-drop(L%*%c(.6,.2,.3,.1,.8));level<-.8
  tail<-(1-level^(1/5))/2
  lo<-map$df*ms/qchisq(1-tail,map$df);hi<-map$df*ms/qchisq(tail,map$df)
  C<-rbind(L,-L,-diag(5));rhs<-c(hi,-lo,rep(0,5))
  vertices<-lapply(combn(nrow(C),5,simplify=FALSE),function(i) {
    z<-tryCatch(solve(C[i,],rhs[i]),error=function(e)NULL)
    if(!is.null(z) && max(C%*%z-rhs)<1e-8)z
  })
  vertices<-do.call(rbind,Filter(Negate(is.null),vertices))
  expect_gt(nrow(vertices),1)
  setup<-.mfrm_mvgt_ms_region_setup(map,level)
  expect_equal((1-2*setup$tail)^5,level)
  for(e in list(c(0,0,0,.5,.5),c(0,.5,.5,.5,.5),c(0,0,0,1/6,1/24))) {
    variance<-drop(vertices%*%e)
    ratio<-vertices[,1]/(vertices[,1]+variance)
    expected<-range(ratio)
    got<-.mfrm_mvgt_ms_region(ms,setup,e)
    expect_equal(got$bounds,expected,tolerance=1e-8)
    got_variance<-.mfrm_mvgt_ms_region(ms,setup,e,ratio=FALSE)
    expect_equal(got_variance$bounds,range(variance),tolerance=1e-8)
    for(scale in c(1e-150,1e150)) {
      scaled<-.mfrm_mvgt_ms_region(ms*scale,setup,e,ratio=FALSE)
      expect_equal(scaled$bounds/scale,got_variance$bounds,tolerance=1e-8)
    }
  }
})

general_interval_study <- function(nested=TRUE,custom=FALSE) {
  x<-if(nested)mfrmr_example_gtheory_nested else mfrmr_example_gtheory_crossed
  if(custom) {
    g<-mfrm_g_study(x,c('Content','Organization'),
      facets=c(Occasion='Task',Assessor='Rater'),nesting=if(nested)c(Assessor='Occasion'))
    grid<-data.frame(Occasion=c(6,12,3,12,6),Assessor=c(3,1,6,6,3))
  } else {
    g<-mfrm_g_study(x,c('Content','Organization'),nesting=if(nested)c(Rater='Task'))
    grid<-data.frame(Raters=c(3,1,6,6,3),Tasks=c(6,12,3,12,6))
  }
  mfrm_d_study(g,grid,weights=c(Content=.5,Organization=.5))
}

test_that('general nested methods retain complete plans, diagnostics and method selection', {
  skip_if_not_installed('lpSolve')
  d<-general_interval_study();before<-serialize(d,NULL)
  custom<-general_interval_study(custom=TRUE)
  for(metric in c('G','Phi')) {
    ci<-confint(d,parm=metric,method=c('approximate_f','mean_square_region'),assumption='normal')
    other<-confint(custom,parm=metric,method=c('approximate_f','mean_square_region'),assumption='normal')
    expect_equal(ci$Lower,other$Lower,tolerance=1e-10)
    expect_equal(ci$Upper,other$Upper,tolerance=1e-10)
    expect_equal(ci$Estimate,rep(subset(d$coefficients,Kind=='Composite')[[metric]],2))
    expect_length(attr(ci,'diagnostics'),10)
    expect_identical(attr(ci,'source_design')$nesting,c(Rater='Task'))
    expect_identical(attr(ci,'source_summary'),d$gstudy[c('design','components','mean_products','degrees_of_freedom')])
    expect_error(confint(d,parm=metric,assumption='normal'),'explicitly select')
    expect_error(confint(d,parm=metric,method='exact_f',assumption='normal'),'Unsupported')
    expect_error(confint(d,parm=metric,reference=1,method='approximate_f',assumption='normal'),'not supported')
  }
  expect_identical(serialize(d,NULL),before)
  expect_error(confint(d,parm='RelativeSEM',method=c('mls','mean_square_region'),assumption='normal'),'Select one')
})

test_that('shared-region difference enclosures retain dependence and ordering for both designs', {
  skip_if_not_installed('lpSolve')
  for(nested in c(FALSE,TRUE)) for(metric in c('G','Phi','RelativeSEM','AbsoluteSEM')) {
    d<-general_interval_study(nested)
    ci<-confint(d,parm=metric,method='mean_square_region',assumption='normal')
    change<-confint(d,parm=metric,reference=1,method='mean_square_region',assumption='normal')
    sem<-grepl('SEM',metric)
    expect_equal(change$Estimate,ci$Estimate[-1]-ci$Estimate[1])
    expect_gte(min(change$Lower-(ci$Lower[-1]-ci$Upper[1])),-1e-10)
    expect_lte(max(change$Upper-(ci$Upper[-1]-ci$Lower[1])),1e-10)
    expect_equal(change$Lower[change$Scenario==5],0)
    expect_equal(change$Upper[change$Scenario==5],0)
    if(sem) expect_lte(change$Upper[change$Scenario==4],0) else
      expect_gte(change$Lower[change$Scenario==4],0)
    expect_identical(attr(change,'reference')$counts,d$design_grid[1,,drop=FALSE])
    expect_identical(attr(change,'individual_intervals')$Lower,ci$Lower)
    expect_length(attr(change,'diagnostics'),4)
    reverse<-confint(d,parm=metric,reference=2,method='mean_square_region',assumption='normal')
    expect_equal(reverse$Lower[reverse$Scenario==1],-change$Upper[change$Scenario==2])
    expect_equal(reverse$Upper[reverse$Scenario==1],-change$Lower[change$Scenario==2])
    if(!sem) {
      a<-.mfrm_mvgt_profile_error(d,1,metric)
      feasible<-t(vapply(seq_len(nrow(d$design_grid))[-1],function(i)
        .mfrm_mvgt_contrast_range(a,.mfrm_mvgt_profile_error(d,i,metric)),numeric(2)))
      expect_equal(change$FeasibleLower,feasible[,1])
      expect_equal(change$FeasibleUpper,feasible[,2])
      expect_true(all(change$Lower>=change$FeasibleLower))
      expect_true(all(change$Upper<=change$FeasibleUpper))
    }
    if(sem) {
      scaled<-mfrm_d_study(d$gstudy,d$design_grid,-7*d$weights)
      got<-confint(scaled,parm=metric,reference=1,method='mean_square_region',assumption='normal')
      expect_equal(got$Lower,7*change$Lower,tolerance=1e-7)
      expect_equal(got$Upper,7*change$Upper,tolerance=1e-7)
    }
  }
})

test_that('analytic coefficient limits remove impossible differences without a refit', {
  skip_if_not_installed('lpSolve')
  x<-mfrmr_example_gtheory_crossed
  x<-subset(x,Person %in% head(unique(Person),10) &
    Rater %in% head(unique(Rater),2) & Task %in% head(unique(Task),3))
  g<-mfrm_g_study(x,'Content')
  d<-mfrm_d_study(g,data.frame(Raters=c(2,3),Tasks=c(3,2)))
  ci<-confint(d,reference=1,method='mean_square_region',assumption='normal')
  limit<-(sqrt(1.5)-1)/(sqrt(1.5)+1)
  expect_equal(c(ci$Lower,ci$Upper),c(-limit,limit),tolerance=1e-12)
  # Near-extrema arise when one Person interaction dominates the other error.
  for(j in 4:5) {
    v<-c(0,0,0,0,0,0,1);v[j]<-1e10
    a<-.mfrm_mvgt_profile_error(d,1,'G');b<-.mfrm_mvgt_profile_error(d,2,'G')
    ea<-sum(a*v);eb<-sum(b*v);v[1]<-sqrt(ea*eb)
    delta<-v[1]/(v[1]+eb)-v[1]/(v[1]+ea)
    expect_equal(abs(delta),limit,tolerance=1e-10)
    expect_gte(delta,ci$Lower);expect_lte(delta,ci$Upper)
  }
})

test_that('unavailable shared regions remain unavailable, including duplicate plans', {
  skip_if_not_installed('lpSolve')
  x<-mfrmr_example_gtheory_nested;x$Content<-x$Organization<-1
  d<-mfrm_d_study(mfrm_g_study(x,c('Content','Organization'),nesting=c(Rater='Task')),
    data.frame(Raters=c(3,3,6),Tasks=c(6,6,12)))
  for(metric in c('G','Phi','RelativeSEM','AbsoluteSEM')) {
    ci<-confint(d,parm=metric,reference=1,method='mean_square_region',assumption='normal')
    expect_false(any(ci$Available));expect_true(all(is.na(ci$Lower)&is.na(ci$Upper)))
    expect_true(all(grepl('Reference interval unavailable',ci$Status)))
    expect_identical(attr(ci,'reference')$status,'Degenerate mean squares')
  }
})

test_that('general nested regions feed saved ribbons and difference forests without refitting', {
  skip_if_not_installed('lpSolve')
  skip_if_not_installed('ggplot2')
  d<-general_interval_study()
  ci<-confint(d,parm='AbsoluteSEM',method='mean_square_region',assumption='normal')
  figure<-plot(d,type='sem',metric='AbsoluteSEM',x_var='Tasks',intervals=ci,draw=FALSE)
  expect_equal(figure$data$series$Lower,ci$Lower)
  expect_equal(figure$data$series$Upper,ci$Upper)
  change<-confint(d,parm='Phi',method='mean_square_region',reference=1,assumption='normal')
  forest<-plot(change,draw=FALSE)
  expect_true(all(grepl('Raters per Task',forest$data$table$Plan)))
  expect_true(any(grepl('Reference plan',forest$data$notes$Text)))
  expect_identical(forest$data$intervals,change)
  path<-tempfile(fileext='.rds');on.exit(unlink(path),add=TRUE)
  saveRDS(list(figure=figure,forest=forest),path)
  saved<-readRDS(path)
  local_mocked_bindings(confint.mfrm_multivariate_d_study=function(...)stop('No interval recalculation'),
    mfrm_g_study=function(...)stop('No refit'),.package='mfrmr')
  for(p in saved) {
    gg<-as_ggplot(p)
    expect_identical(plot_data(gg,'table'),p$data$table)
    expect_s3_class(ggplot2::ggplotGrob(gg),'gtable')
  }
  expect_gte(as_ggplot(saved$figure)$scales$get_scales('y')$limits[2],max(ci$Upper))
})
