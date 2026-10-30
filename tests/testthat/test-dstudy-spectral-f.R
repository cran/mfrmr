spectral_data <- function() read.csv(system.file('extdata','mgenova-table12.csv',package='mfrmr'))

test_that('spectral G reduces to balanced F with original or composite scores', {
  x <- spectral_data()
  g <- mfrm_g_study(x,c('V','W'),rater=NULL)
  d <- mfrm_d_study(g,data.frame(Tasks=c(1,3,6,12)),c(V=.4,W=-.6))
  original <- serialize(d,NULL)
  for (score in list(NULL,'V')) {
    ci <- confint(d,score=score,method='spectral_f',assumption='normal')
    old <- confint(d,score=score,assumption='normal')
    expect_equal(ci[c('Estimate','Lower','Upper')],old[c('Estimate','Lower','Upper')],
      ignore_attr=TRUE,tolerance=1e-10)
    expect_equal(attr(ci,'pivot')$eigenvalues,rep(6,9),tolerance=1e-10)
    expect_equal(attr(ci,'pivot')$residual_df,45)
    expect_false('IntervalCenter' %in% names(ci))
    expect_identical(attr(ci,'source_data'),g$data)
    a <- confint(d,score=score,reference=2,method='spectral_f',assumption='normal')
    b <- confint(d,score=score,reference=2,assumption='normal')
    expect_equal(a[c('Estimate','Lower','Upper')],b[c('Estimate','Lower','Upper')],
      ignore_attr=TRUE,tolerance=1e-10)
  }
  both <- confint(d,method=c('exact_f','spectral_f'),assumption='normal')
  expect_equal(unique(both$MethodCode),c('exact_f','spectral_f'))
  expect_identical(attr(both,'method_results')$spectral_f,
    confint(d,method='spectral_f',assumption='normal'))
  expect_no_error(summary(d,intervals=both))
  expect_no_error(plot(both,draw=FALSE))
  expect_identical(serialize(d,NULL),original)
})

test_that('incomplete spectral endpoints agree with independent GLS quadratic forms', {
  x <- spectral_data();x <- x[-c(1,8,23),]
  g <- mfrm_g_study(x,c('V','W'),rater=NULL,method='minque0')
  d <- mfrm_d_study(g,data.frame(Tasks=c(2,6,12)),c(V=.4,W=-.6))
  ci <- confint(d,method='spectral_f',assumption='normal')
  y <- .4*x$V-.6*x$W
  z <- model.matrix(~factor(Person)-1,x);f <- model.matrix(~factor(Task)-1,x)
  fit <- lm(y~factor(Person)+factor(Task),x)
  nu <- df.residual(fit);sse <- sum(residuals(fit)^2)
  r <- nrow(x)-ncol(f)-nu
  pivot <- attr(ci,'pivot')
  for (j in 1:2) if(is.finite(pivot$log_bounds[j])) {
    tau <- exp(pivot$log_bounds[j]);v <- diag(nrow(x))+tau*tcrossprod(z)
    inverse <- solve(v)
    beta <- solve(crossprod(f,inverse%*%f),crossprod(f,inverse%*%y))
    e <- y-f%*%beta
    q <- drop(crossprod(e,inverse%*%e))
    expect_equal((q-sse)/r/(sse/nu),qf(c(.975,.025)[j],r,nu),tolerance=1e-9)
  }
  expected <- outer(ci$Tasks,exp(pivot$log_bounds),function(m,t)m*t/(1+m*t))
  expect_equal(ci$Lower,expected[,1]);expect_equal(ci$Upper,expected[,2])
  expect_identical(ci$Estimate,subset(d$coefficients,Kind=='Composite')$G)
  expect_true(all(ci$Available))
  expect_no_error(summary(d,intervals=ci))
  expect_no_error(plot(d,metric='G',intervals=ci,draw=FALSE))
  expect_no_error(plot(d,metric='G',intervals=TRUE,method='spectral_f',assumption='normal',draw=FALSE))
  expect_error(confint(d,assumption='normal'),'explicitly select')
})

test_that('spectral paired sets include interior extrema and preserve failed bases', {
  x <- spectral_data();x <- x[-c(1,8,23),]
  g <- mfrm_g_study(x,'V',rater=NULL,method='minque0')
  d <- mfrm_d_study(g,data.frame(Tasks=c(6,12,6,3)))
  ci <- confint(d,method='spectral_f',reference=1,assumption='normal')
  rho <- plogis(attr(ci,'pivot')$log_bounds)
  for(i in seq_len(nrow(ci))) {
    effect <- function(p) ci$Tasks[i]*p/(1+(ci$Tasks[i]-1)*p)-6*p/(1+5*p)
    values <- c(effect(rho),optimize(effect,rho,tol=1e-12)$objective,
      optimize(effect,rho,maximum=TRUE,tol=1e-12)$objective)
    expect_equal(c(ci$Lower[i],ci$Upper[i]),range(values),tolerance=1e-9)
  }
  expect_identical(ci$Status[2],'Identical counts')
  reverse <- confint(d,method='spectral_f',reference=2,assumption='normal')
  expect_equal(c(reverse$Lower[1],reverse$Upper[1]),-c(ci$Upper[1],ci$Lower[1]))
  # Residual-only observations give an empty nonnegative Person-ratio set.
  x$V <- residuals(lm(sin(seq_len(nrow(x)))~factor(Person)+factor(Task),x))
  g <- mfrm_g_study(x,'V',rater=NULL,method='minque0')
  d <- mfrm_d_study(g,data.frame(Tasks=c(3,6,3)))
  empty <- confint(d,method='spectral_f',reference=1,assumption='normal')
  expect_false(any(empty$Available));expect_true(all(is.na(empty$Lower)))
  expect_true(all(empty$Status=='Empty parameter-space interval'))
})

test_that('spectral G keeps numerical and original-point limitations visible', {
  x <- data.frame(Person=c(1,2,2,3,3,4),Task=c(1,1,2,2,3,3),V=c(2,1,4,3,6,5))
  g <- mfrm_g_study(x,'V',rater=NULL,method='minque0')
  d <- mfrm_d_study(g,data.frame(Tasks=c(3,6)))
  ci <- confint(d,method='spectral_f',assumption='normal')
  expect_false(any(ci$Available));expect_true(all(grepl('No residual contrasts',ci$Status)))
  x <- spectral_data();x$V <- x$Person-x$Task
  d <- mfrm_d_study(mfrm_g_study(x,'V',rater=NULL))
  ci <- confint(d,method='spectral_f',assumption='normal')
  expect_false(any(ci$Available));expect_true(all(grepl('numerical resolution',ci$Status)))
  x <- x[-c(1,8,23),];x$V <- x$Person-x$Task+.01*sin(seq_len(nrow(x)))
  g <- mfrm_g_study(x,'V',rater=NULL,method='minque0')
  d <- mfrm_d_study(g,data.frame(Tasks=c(3,6)))
  ci <- confint(d,method='spectral_f',assumption='normal')
  expect_true(all(is.na(ci$Estimate)));expect_true(all(ci$Available))
  expect_true(all(ci$NegativeComponent));expect_no_error(plot(d,metric='G',intervals=ci,draw=FALSE))
  x <- data.frame(Person=c(1,1,2,2,3),Task=c(1,2,1,2,1),V=c(1,3,2,3,5))
  g <- mfrm_g_study(x,'V',rater=NULL,method='minque0')
  d <- mfrm_d_study(g,data.frame(Tasks=c(1,2)))
  ci <- confint(d,method='spectral_f',assumption='normal')
  expect_equal(attr(ci,'pivot')$residual_df,1)
  expect_true(all(ci$Available));expect_true(all(ci$Lower<=ci$Upper))
  x <- spectral_data();x$V <- x$Person-x$Task+1e-9*sin(seq_len(nrow(x)))
  g <- mfrm_g_study(x,'V',rater=NULL)
  d <- mfrm_d_study(g,data.frame(Tasks=c(3,6)))
  ci <- confint(d,method='spectral_f',assumption='normal')
  expect_true(attr(ci,'pivot')$valid)
  expect_false(any(ci$Available))
  expect_true(all(grepl('below numeric resolution',ci$Status)))
})

test_that('spectral G is invariant to score units, labels and row order', {
  x <- spectral_data();x <- x[-c(1,8,23),]
  run <- function(x,w) {
    g <- mfrm_g_study(x,c('V','W'),facets=c(Assessor='Task'),method='minque0')
    d <- mfrm_d_study(g,data.frame(Assessor=c(2,4)),weights=w)
    confint(d,method='spectral_f',assumption='normal')
  }
  ci <- run(x,c(V=.4,W=-.6))
  for (scale in c(-3,1e-80,1e80)) {
    other <- run(x,scale*c(W=-.6,V=.4))
    expect_equal(other[c('Estimate','Lower','Upper')],ci[c('Estimate','Lower','Upper')],
      ignore_attr=TRUE,tolerance=1e-9)
    transformed <- x;transformed$V <- x$V*scale;transformed$W <- x$W*scale
    other <- run(transformed,c(V=.4,W=-.6))
    expect_equal(other[c('Lower','Upper')],ci[c('Lower','Upper')],ignore_attr=TRUE,tolerance=1e-9)
  }
  x$Person <- paste0('P',x$Person);x$Task <- paste0('T',x$Task)
  other <- run(x[nrow(x):1,],c(V=.4,W=-.6))
  expect_equal(other[c('Estimate','Lower','Upper')],ci[c('Estimate','Lower','Upper')],ignore_attr=TRUE,tolerance=1e-9)
})

test_that('spectral intervals bind to raw source data and replay without calculation', {
  x <- spectral_data();x <- x[-c(1,8,23),]
  g <- mfrm_g_study(x,'V',rater=NULL,method='minque0')
  d <- mfrm_d_study(g,data.frame(Tasks=c(3,6)))
  ci <- confint(d,method='spectral_f',assumption='normal')
  bad <- ci;attr(bad,'source_data')$V[1] <- 100
  expect_error(summary(d,intervals=bad),'ratings do not match')
  expect_error(plot(d,metric='G',intervals=bad,draw=FALSE),'ratings do not match')
  for(metric in c('Phi','RelativeSEM','AbsoluteSEM'))
    expect_error(confint(d,parm=metric,method='spectral_f',assumption='normal'),'Unsupported')
  crossed <- mfrm_d_study(mfrm_g_study(mfrmr_example_gtheory_crossed,'Content'))
  expect_error(confint(crossed,method='spectral_f',assumption='normal'),'Unsupported')
  path <- tempfile(fileext='.rds');on.exit(unlink(path))
  saveRDS(list(d=d,ci=ci),path);saved <- readRDS(path)
  local_mocked_bindings(.mfrm_mvgt_spectral_ratio=function(...)stop('Unexpected recalculation'),.package='mfrmr')
  expect_no_error(summary(saved$d,intervals=saved$ci))
  expect_no_error(plot(saved$d,metric='G',intervals=saved$ci,draw=FALSE))
  expect_no_error(plot(saved$ci,draw=FALSE))
})

test_that('spectral F quantile resolution cannot become a zero-width confidence interval', {
  for (sparse in c(FALSE, TRUE)) {
    x <- spectral_data()
    if (sparse) x <- x[-c(1,8,23),]
    g <- mfrm_g_study(x,'V',rater=NULL,method=if(sparse)'minque0' else 'anova')
    d <- mfrm_d_study(g,data.frame(Tasks=c(3,6,3)))
    for (level in c(1e-12,1e-16,1e-20)) {
      ci <- confint(d,method='spectral_f',assumption='normal',level=level)
      expect_false(any(ci$Available))
      expect_true(all(is.na(ci$Lower) & is.na(ci$Upper)))
      expect_identical(ci$Estimate,d$coefficients$G)
      expect_true(all(ci$Status=='F quantile interval below numeric resolution'))
      expect_true(all(ci$Validation=='No interval supplied'))
      pivot <- attr(ci,'pivot')
      expect_lt(pivot$endpoint_tolerance,pivot$quantile_resolution)
      difference <- confint(d,method='spectral_f',assumption='normal',level=level,reference=1)
      expect_false(any(difference$Available))
      expect_true(all(difference$Status==ci$Status[1]))
    }
    # A narrow but resolved interval still agrees with an independent GLS pivot.
    ci <- confint(d,method='spectral_f',assumption='normal',level=1e-6)
    expect_true(all(ci$Available & ci$Lower<ci$Upper))
    pivot <- attr(ci,'pivot')
    zp <- model.matrix(~factor(Person)-1,x)
    xf <- model.matrix(~factor(Task)-1,x)
    residual <- residuals(lm(V~factor(Person)+factor(Task),x))
    sse <- sum(residual^2)
    for (j in 1:2) {
      tau <- exp(pivot$log_bounds[j])
      inverse <- solve(diag(nrow(x))+tau*tcrossprod(zp))
      beta <- solve(crossprod(xf,inverse%*%xf),crossprod(xf,inverse%*%x$V))
      e <- x$V-xf%*%beta
      f <- (drop(crossprod(e,inverse%*%e))-sse)/pivot$person_df/(sse/pivot$residual_df)
      expect_lt(abs(log(f)-pivot$log_critical[j]),pivot$endpoint_tolerance)
    }
  }
})

test_that('spectral projection isolates independent normal quadratic forms including zero Person variance', {
  x <- spectral_data()[-c(1,8,23),]
  g <- mfrm_g_study(x,'V',rater=NULL,method='minque0')
  pivot <- mfrmr:::.mfrm_mvgt_spectral_ratio(g,c(V=1),.95)
  zp <- model.matrix(~factor(Person)-1,x)
  xf <- model.matrix(~factor(Task)-1,x)
  # Independent full facet-complement eigendecomposition, not the thin SVD.
  h <- qr.Q(qr(xf),complete=TRUE)[,-seq_len(ncol(xf)),drop=FALSE]
  decomposition <- eigen(crossprod(h,zp)%*%crossprod(zp,h),symmetric=TRUE)
  r <- pivot$person_df
  u <- h%*%decomposition$vectors[,seq_len(r),drop=FALSE]
  v <- h%*%decomposition$vectors[,-seq_len(r),drop=FALSE]
  expect_equal(decomposition$values[seq_len(r)],pivot$eigenvalues,tolerance=1e-12)
  expect_equal(ncol(v),pivot$residual_df)
  for (person_variance in c(0,.01,3)) {
    covariance <- person_variance*tcrossprod(zp)+7*tcrossprod(xf)+diag(2,nrow(x))
    expect_equal(unname(crossprod(u,covariance%*%u)),
      diag(2+person_variance*pivot$eigenvalues),tolerance=1e-12)
    expect_equal(unname(crossprod(v,covariance%*%v)),diag(2,ncol(v)),tolerance=1e-12)
    expect_lt(max(abs(crossprod(u,covariance%*%v))),1e-12)
  }
})

test_that('a genuine zero-only spectral set is distinct from collapsed finite bounds', {
  g <- mfrm_g_study(spectral_data(),'V',rater=NULL)
  d <- mfrm_d_study(g,data.frame(Tasks=c(3,6,3)))
  # F(0) equal to the lower critical value gives this legitimate intersection.
  local_mocked_bindings(.mfrm_mvgt_spectral_ratio=function(...)
    list(log_bounds=c(-Inf,-Inf),valid=TRUE,status='Upper endpoint at zero'),.package='mfrmr')
  ci <- confint(d,method='spectral_f',assumption='normal')
  expect_true(all(ci$Available))
  expect_equal(ci$Lower,rep(0,3)); expect_equal(ci$Upper,rep(0,3))
  expect_identical(ci$Status,rep('Upper endpoint at zero',3))
  change <- confint(d,method='spectral_f',assumption='normal',reference=1)
  expect_true(all(change$Available))
  expect_equal(change$Lower,c(0,0)); expect_equal(change$Upper,c(0,0))
  expect_identical(change$Status,c('Upper endpoint at zero','Identical counts'))
})
