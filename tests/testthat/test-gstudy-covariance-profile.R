joint_profile_fixture <- function() {
  tasks <- read.csv(system.file('extdata','mgenova-table12.csv',package='mfrmr',mustWork=TRUE))
  tasks <- tasks[(tasks$Person+tasks$Task)%%3!=0,]
  mfrm_g_study(tasks,c('V','W'),rater=NULL,method='minque0')
}

joint_dense_reference <- function(y,groups,factors) {
  n <- nrow(y); x <- kronecker(diag(2),rep(1,n)); yy <- as.vector(y)
  v <- Reduce('+',Map(function(f,g)kronecker(tcrossprod(f),outer(g,g,'==')*1),factors,groups))
  vi <- solve(v); vix <- vi%*%x; information <- crossprod(x,vix)
  p <- vi-vix%*%solve(information,t(vix)); py <- p%*%yy
  value <- as.numeric(determinant(v,logarithm=TRUE)$modulus)+
    as.numeric(determinant(information,logarithm=TRUE)$modulus)+
    drop(crossprod(yy,py))+(2*n-2)*log(2*pi)
  a <- p-tcrossprod(py)
  gradient <- lapply(groups,function(g) {
    kernel <- outer(g,g,'==')*1
    matrix(c(sum(a[1:n,1:n]*kernel),sum(a[(n+1):(2*n),1:n]*kernel),
      sum(a[1:n,(n+1):(2*n)]*kernel),sum(a[(n+1):(2*n),(n+1):(2*n)]*kernel)),2)
  })
  list(objective=value,gradient=gradient)
}

test_that('joint selection checks equivalent retained candidates without relaxing diagnostics', {
  # Selection logic only: numerical gradient/likelihood checks are independent below.
  candidate <- function(id,value=0,convergence=0L)
    list(par=c(id,0,.4,0,0,0),objective=value,convergence=convergence)
  first <- candidate(.4); good <- candidate(.5,5e-8)
  attempts <- list(first,good,candidate(.6,.01))
  source <- list(groups=list(1:3,1:3),initial=rep(list(diag(2)),2),objective=function(f)0)
  local_mocked_bindings(
    .mfrm_mvgt_profile_minimize=function(...)list(best=attempts[[1]],attempts=attempts),
    .mfrm_mvgt_profile_gradient=function(par,...)list(valid=par[1]>=.5),
    .mfrm_mvgt_joint_cone_gradient=function(factors,...)list(valid=TRUE,objective=0,
      gradient=rep(list(matrix(0,2,2)),2)))
  fit <- .mfrm_mvgt_joint_fit(source)
  expect_true(fit$valid);expect_identical(fit$fit,good)
  expect_identical(fit$attempts,attempts)
  expect_identical(fit$selection$selected,2L)
  expect_named(fit$selection$checks,c('1','2'))
  expect_false(fit$selection$checks[['1']]$valid)
  expect_equal(fit$selection$objective_tolerance,1e-7)
  # A more favorable-looking diagnostic cannot override a materially lower deviance.
  attempts[[2]]$objective <- 1e-5
  fit <- .mfrm_mvgt_joint_fit(source)
  expect_false(fit$valid);expect_identical(fit$fit,first)
  expect_named(fit$selection$checks,'1')
  # Failed optimizer convergence remains disqualifying, even at an equal objective.
  attempts[[2]] <- candidate(.5,5e-8,1L)
  fit <- .mfrm_mvgt_joint_fit(source)
  expect_false(fit$valid);expect_named(fit$selection$checks,'1')
  attempts <- list(good,candidate(.6,0))
  fit <- .mfrm_mvgt_joint_fit(source)
  expect_true(fit$valid);expect_identical(fit$fit,good)
  expect_named(fit$selection$checks,'1')
  # Preserve a failed finite search and an all-nonfinite optimizer result.
  attempts <- list(first,candidate(.3,0))
  fit <- .mfrm_mvgt_joint_fit(source)
  expect_false(fit$valid);expect_identical(fit$fit,first)
  expect_named(fit$selection$checks,c('1','2'))
  attempts <- list(candidate(.4,Inf,1L))
  fit <- .mfrm_mvgt_joint_fit(source)
  expect_false(fit$valid);expect_identical(fit$attempts,attempts)
})

test_that('joint profiles preserve original covariances and independently checked endpoints', {
  skip_if_not_installed('lme4');skip_if_not_installed('minqa')
  g <- joint_profile_fixture(); before <- serialize(g,NULL)
  ci <- confint(g,score=c('V','W'),method='reml_profile',assumption='normal')
  expect_identical(serialize(g,NULL),before)
  expect_equal(ci$Estimate,unname(vapply(g$components,function(m)m[1,2],numeric(1))))
  expect_true(all(ci$Available));expect_true(all(is.na(ci$SE)))
  expect_true(all(is.na(ci$RawLower)));expect_true(all(is.na(ci$RawUpper)))
  expect_true(any(ci$Lower<0));expect_true(all(ci$Lower<=ci$IntervalCenter & ci$IntervalCenter<=ci$Upper))
  expect_identical(attr(ci,'kind'),'Covariance');expect_identical(attr(ci,'score_pair'),c('V','W'))
  expect_match(paste(attr(ci,'notes'),collapse=' '),'more persons alone')
  source <- attr(ci,'reml_source');expect_null(source$objective)
  expected <- joint_dense_reference(source$y,source$groups,source$fit$factors)
  expect_equal(source$fit$fit$objective,expected$objective,tolerance=1e-9)
  diagnostics <- attr(ci,'diagnostics');expect_identical(names(diagnostics),ci$Source)
  for(z in diagnostics)for(endpoint in z$bounds) {
    f <- z$profiles[[sprintf('%.17g',endpoint)]]
    expect_true(f$valid);expect_equal(f$tau,endpoint)
    expect_equal(f$lr,z$critical,tolerance=1e-5)
    expect_equal(tcrossprod(f$factors[[z$target]])[1,2],endpoint,tolerance=1e-10)
    expect_equal(f$fit$objective,joint_dense_reference(source$y,source$groups,f$factors)$objective,tolerance=1e-9)
  }
  selected <- confint(g,parm=c('Residual','Task'),score=c('V','W'),method='reml_profile',assumption='normal')
  expect_identical(selected$Source,c('Residual','Task'))
  expect_identical(names(attr(selected,'diagnostics')),selected$Source)
  expect_equal(selected$Lower,ci$Lower[c(3,2)],tolerance=1e-6)
})

test_that('covariance-cone derivatives expose a factor-coordinate stationary nonoptimum', {
  skip_if_not_installed('lme4');skip_if_not_installed('minqa')
  source <- .mfrm_mvgt_joint_source(joint_profile_fixture(),c('V','W'))
  f <- rep(list(matrix(0,2,2)),length(source$groups));f[[length(f)]]<-diag(2)
  derivative <- .mfrm_mvgt_joint_cone_gradient(f,source$objective)
  independent <- joint_dense_reference(source$y,source$groups,f)
  expect_true(derivative$valid)
  expect_equal(unname(unlist(derivative$gradient)),unname(unlist(independent$gradient)),tolerance=1e-5)
  expect_lt(min(eigen(derivative$gradient[[1]],symmetric=TRUE)$values),-.002)
  coordinates <- .mfrm_mvgt_joint_coordinates(length(f));par <- coordinates$encode(f)
  factor_gradient <- .mfrm_mvgt_profile_gradient(par,function(p)source$objective(coordinates$decode(p)),
    lower=coordinates$lower)
  expect_equal(factor_gradient$gradient[seq_len(3*(length(f)-1))],rep(0,3*(length(f)-1)),tolerance=1e-8)
  cone <- .mfrm_mvgt_joint_cone_check(f,function(...)derivative)
  expect_false(all(cone$table$Valid))
})

test_that('joint profiles respect score order units and signs', {
  skip_if_not_installed('lme4');skip_if_not_installed('minqa')
  g <- joint_profile_fixture()
  call <- function(g,score=c('V','W'))confint(g,parm='Task',score=score,method='reml_profile',assumption='normal')
  original <- call(g); reversed <- call(g,c('W','V'))
  expect_true(original$Available && reversed$Available)
  expect_equal(reversed[c('Estimate','IntervalCenter','Lower','Upper')],
    original[c('Estimate','IntervalCenter','Lower','Upper')],tolerance=1e-5)
  d <- g$data; d$V <- d$V*-.01+10; d$W <- d$W*100
  changed <- call(mfrm_g_study(d,c('V','W'),rater=NULL,method='minque0'))
  expect_true(changed$Available)
  expect_equal(changed$Estimate,-original$Estimate,tolerance=1e-8)
  expect_equal(changed$IntervalCenter,-original$IntervalCenter,tolerance=1e-5)
  expect_equal(c(changed$Lower,changed$Upper),-c(original$Upper,original$Lower),tolerance=1e-5)
})

test_that('joint profiles reject mismatches and retain degenerate pairs', {
  skip_if_not_installed('lme4');skip_if_not_installed('minqa')
  g <- joint_profile_fixture()
  call <- function(x,...)confint(x,score=c('V','W'),method='reml_profile',assumption='normal',...)
  expect_error(call(g,draws=1000),'only')
  expect_error(call(g,parm='Absent'),'parm')
  bad <- g;bad$data <- NULL;expect_error(call(bad),'used ratings')
  bad <- g;bad$components$Person[1,2] <- bad$components$Person[2,1] <- 10
  expect_error(call(bad),'disagree')
  d <- g$data;d$W <- d$V
  singular <- mfrm_g_study(d,c('V','W'),rater=NULL,method='minque0')
  z <- call(singular)
  expect_false(any(z$Available));expect_true(all(is.na(z$Lower)))
  expect_true(all(grepl('Degenerate',z$Status)))
  source <- .mfrm_mvgt_joint_source(g,c('V','W'));source$valid<-FALSE;source$status<-'Unresolved source'
  rejected <- .mfrm_mvgt_joint_set(source,1L,.95)
  expect_true(all(is.na(rejected$bounds)));expect_match(rejected$status,'Unresolved')
})

test_that('covariance profile plots retain diagnostics and replay without fitting', {
  skip_if_not_installed('lme4');skip_if_not_installed('minqa');skip_if_not_installed('ggplot2')
  g <- joint_profile_fixture()
  p <- plot(g,type='covariance',score=c('V','W'),intervals=TRUE,
    method='reml_profile',assumption='normal',draw=FALSE)
  tab <- plot_data(p,'table');ci <- plot_data(p,'intervals')[[1]]
  expect_equal(tab$Estimate,ci$Estimate);expect_equal(tab$IntervalCenter,ci$IntervalCenter)
  expect_identical(tab$Lower,ci$Lower);expect_identical(tab$SEStatus,ci$SEStatus)
  expect_equal(nrow(p$data$panels),1)
  restored <- unserialize(serialize(p,NULL))
  local_mocked_bindings(.mfrm_mvgt_joint_source=function(...)stop('Unexpected refit'),
    confint.mfrm_multivariate_gstudy=function(...)stop('Unexpected interval calculation'))
  gg <- as_ggplot(restored)
  expect_identical(plot_data(gg,'table'),tab)
  expect_identical(plot_data(gg,'intervals')[[1]],ci)
  expect_null(gg$labels$title);expect_null(gg$labels$caption)
  expect_no_warning(ggplot2::ggplotGrob(gg))
})
