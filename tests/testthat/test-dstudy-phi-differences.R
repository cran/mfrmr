phi_difference_source <- function(counts=c(6,12,3,6), weights=c(V=.5,W=.5)) {
  path <- test_path("..","..","inst","extdata","mgenova-table12.csv")
  if(!file.exists(path)) path <- system.file("extdata","mgenova-table12.csv",package="mfrmr")
  g <- mfrm_g_study(read.csv(path),c("V","W"),rater=NULL)
  mfrm_d_study(g,data.frame(Tasks=counts),weights=weights)
}

test_that("Phi differences map the complete base set, including interior extrema", {
  d<-phi_difference_source();before<-serialize(d,NULL)
  base<-confint(d,parm="Phi",assumption="normal",draws=5000,seed=29)
  ci<-confint(d,parm="Phi",reference=1,assumption="normal",draws=5000,seed=29)
  bounds<-attr(base,"pivot")$single_rating_bounds
  counts<-d$design_grid$Tasks
  for(i in seq_len(nrow(ci))) {
    b<-counts[ci$Scenario[i]];a<-counts[1]
    f<-function(rho)rho/(rho+(1-rho)/b)-rho/(rho+(1-rho)/a)
    extrema<-c(f(bounds),optimize(f,bounds,tol=1e-12)$objective,
      optimize(f,bounds,maximum=TRUE,tol=1e-12)$objective)
    expect_equal(c(ci$Lower[i],ci$Upper[i]),range(extrema),tolerance=1e-10)
    grid<-seq(bounds[1],bounds[2],length.out=501)
    expect_true(all(f(grid)>=ci$Lower[i]-1e-12 & f(grid)<=ci$Upper[i]+1e-12))
  }
  expect_equal(ci$Upper[1],tanh(log(2)/4))
  expect_equal(ci$Lower[3],0);expect_equal(ci$Upper[3],0)
  expect_equal(ci$Status[3],"Identical counts")
  expect_equal(ci$Estimate,base$Estimate[-1]-base$Estimate[1])
  expect_identical(ci$Metric,rep("Phi difference",3))
  expect_identical(attr(ci,"pivot"),attr(base,"pivot"))
  expect_equal(attr(ci,"reference")$value,base$Estimate[1])
  expect_identical(serialize(d,NULL),before)
  expect_match(ci$Validation[1],"not exact")
  expect_match(paste(attr(ci,"notes"),collapse=" "),"interior extremum")
})

test_that("Phi reference reversal, scaling, selection, RNG and saved replay hold", {
  d<-phi_difference_source();set.seed(432);before<-.Random.seed
  first<-confint(d,parm="Phi",reference=1,assumption="normal",draws=5000,seed=29)
  expect_identical(.Random.seed,before)
  reverse<-confint(d,parm="Phi",reference=2,assumption="normal",draws=5000,seed=29)
  expect_equal(reverse$Lower[1],-first$Upper[1]);expect_equal(reverse$Upper[1],-first$Lower[1])
  expect_equal(reverse$Estimate[1],-first$Estimate[1])
  wider<-confint(d,parm="Phi",reference=1,assumption="normal",level=.99,draws=20000,seed=29)
  same_draws<-confint(d,parm="Phi",reference=1,assumption="normal",draws=20000,seed=29)
  expect_true(all(wider$Lower<=same_draws$Lower & wider$Upper>=same_draws$Upper))
  scaled<-phi_difference_source(weights=c(V=-5,W=-5))
  scaled_ci<-confint(scaled,parm="Phi",reference=1,assumption="normal",draws=5000,seed=29)
  expect_equal(scaled_ci$Lower,first$Lower);expect_equal(scaled_ci$Upper,first$Upper)
  path<-tempfile(fileext=".rds");on.exit(unlink(path),add=TRUE)
  saveRDS(list(d=d,ci=first),path);saved<-readRDS(path);saved$d$gstudy$data<-NULL
  expect_identical(confint(saved$d,parm="Phi",reference=1,assumption="normal",draws=5000,seed=29),saved$ci)
  for(ref in list(0,5,NA_real_,1.5,TRUE,matrix(1)))
    expect_error(confint(d,parm="Phi",reference=ref,assumption="normal"),"reference")
})

test_that("Phi difference MC bands protect both endpoints of the whole image", {
  d<-phi_difference_source();ci<-confint(d,parm="Phi",reference=1,
    assumption="normal",draws=5000,seed=29)
  mc<-attr(ci,"monte_carlo");base<-mc$base_endpoint_bands;bands<-mc$endpoint_bands
  expect_match(mc$band_scope,"jointly")
  expect_identical(bands$Scenario,ci$Scenario)
  # Independently rebuild order ranks from the corrected pivot's same draws.
  p<-.mfrm_mvgt_anova_projection(d$gstudy,weights=d$weights)
  z<-with_preserved_rng_seed(29,sort(.mfrm_mvgt_phi_pivots(p$ms,10,6,5000)))
  probs<-c(.025,.975)
  for(j in 1:2) {
    ranks<-qbinom(c(.0125,.9875),5000,probs[j])+c(0,1)
    expect_equal(unname(base[j,]),c(0,z,1)[ranks+1])
    coverage<-pbinom(ranks[2]-1,5000,probs[j])-pbinom(ranks[1]-1,5000,probs[j])
    expect_gte(coverage,.975)
  }
  for(lo in seq(base[1,1],base[1,2],length.out=7))
    for(hi in seq(base[2,1],base[2,2],length.out=7)) {
      mapped<-.mfrm_mvgt_count_difference_range(qlogis(c(lo,hi)),d$design_grid$Tasks,6,"Available")
      expect_true(all(mapped$lower[-1]>=bands$LowerMCLower-1e-12 &
        mapped$lower[-1]<=bands$LowerMCUpper+1e-12))
      expect_true(all(mapped$upper[-1]>=bands$UpperMCLower-1e-12 &
        mapped$upper[-1]<=bands$UpperMCUpper+1e-12))
    }
})

test_that("Phi differences preserve zero, unavailable and collapsed base intervals", {
  d<-phi_difference_source()
  # A deterministic pivot is used only to exercise atom/overlapping-band states.
  testthat::local_mocked_bindings(.mfrm_mvgt_phi_pivots=function(ms,n,k,draws)
    rep(0,draws),.package="mfrmr")
  zero<-confint(d,parm="Phi",reference=1,assumption="normal",draws=5000)
  expect_equal(zero$Lower,rep(0,3));expect_equal(zero$Upper,rep(0,3))
  expect_equal(zero$Status,c("Upper endpoint at zero","Upper endpoint at zero","Identical counts"))
  expect_equal(attr(zero,"monte_carlo")$zero_proportion,1)
  data<-expand.grid(Person=1:10,Task=1:6);data$Score<-data$Person
  g<-mfrm_g_study(data,"Score",rater=NULL)
  bad<-confint(mfrm_d_study(g,data.frame(Tasks=c(6,12,6))),
    parm="Phi",reference=1,assumption="normal",draws=5000)
  expect_true(all(is.na(bad$Lower)&is.na(bad$Upper)))
  expect_equal(bad$Status,rep("Degenerate mean squares",2))
  expect_equal(attr(bad,"monte_carlo")$draws_used,0L)
})

test_that("Phi difference forests retain negative changes and saved diagnostics", {
  ci<-confint(phi_difference_source(),parm="Phi",reference=1,assumption="normal",draws=5000)
  p<-plot(ci,draw=FALSE)
  expect_equal(p$data$xrange,c(-1,1));expect_identical(p$data$intervals,ci)
  expect_lt(p$data$table$Upper[p$data$table$Tasks==3],0)
  expect_equal(p$data$xlab,"Phi difference")
  expect_true(any(grepl("Reference plan: 1: Tasks = 6",p$data$notes$Text,fixed=TRUE)))
  if(requireNamespace("ggplot2",quietly=TRUE)) {
    gg<-as_ggplot(p)
    expect_s3_class(gg,"ggplot")
    expect_silent(ggplot2::ggplot_build(gg))
  }
})

test_that("Phi differences retain intervals when Person points are unavailable", {
  n<-10;k<-6;p<-qr.Q(qr(contr.helmert(n)));t<-qr.Q(qr(contr.helmert(k)))
  data<-expand.grid(Person=1:n,Task=1:k)
  data$Score<-as.vector(outer(p[,1]*sqrt(.5*(n-1)/k),rep(1,k)) +
    outer(rep(1,n),t[,1]*sqrt((k-1)/n)) +
    sqrt((n-1)*(k-1))*outer(p[,2],t[,2]))
  g<-mfrm_g_study(data,"Score",rater=NULL)
  ci<-confint(mfrm_d_study(g,data.frame(Tasks=c(6,12))),
    parm="Phi",reference=1,assumption="normal",draws=5000,seed=29)
  expect_true(is.na(ci$Estimate));expect_equal(ci$PointStatus,"Point difference unavailable")
  expect_equal(ci$Lower,0);expect_gt(ci$Upper,0)
  expect_true(ci$Available);expect_true(ci$NegativeComponent)
  expect_equal(attr(ci,"monte_carlo")$base_status,"Includes zero")
})

test_that("overlapping Monte Carlo bands return a conservative image envelope", {
  d<-phi_difference_source()
  testthat::local_mocked_bindings(.mfrm_mvgt_phi_pivots=function(ms,n,k,draws)
    seq(.1,.9,length.out=draws),.package="mfrmr")
  ci<-confint(d,parm="Phi",reference=1,assumption="normal",draws=1000,level=.001)
  mc<-attr(ci,"monte_carlo");base<-mc$base_endpoint_bands;bands<-mc$endpoint_bands
  expect_gt(base[1,2],base[2,1])
  expect_equal(bands$LowerMCLower,bands$UpperMCLower)
  expect_equal(bands$LowerMCUpper,bands$UpperMCUpper)
  expect_true(all(ci$Lower>=bands$LowerMCLower & ci$Lower<=bands$LowerMCUpper))
  expect_true(all(ci$Upper>=bands$UpperMCLower & ci$Upper<=bands$UpperMCUpper))
})
