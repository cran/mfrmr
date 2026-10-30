sem_profile_source <- function() {
  x <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
  x <- x[(x$Person+x$Task) %% 3 != 0, ]
  mfrm_g_study(x, c("V", "W"), rater = NULL, method = "minque0")
}

test_that("absolute REML constraints match a dense covariance calculation", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  g <- sem_profile_source()
  s <- .mfrm_mvgt_profile_source(g, c(V = .4, W = .6))
  expect_equal(s$score_scale, sd(.4*g$data$V+.6*g$data$W), tolerance = 1e-12)
  dense <- function(theta) {
    V <- Reduce(`+`, Map(function(group, v) v*(outer(group, group, `==`)*1), s$groups, theta))
    inv <- solve(V); one <- rowSums(inv); y <- s$y
    as.numeric(determinant(V, logarithm = TRUE)$modulus)+log(sum(one))+
      (length(y)-1)*log(2*pi)+drop(crossprod(y, inv %*% y))-sum(one*y)^2/sum(one)
  }
  expect_equal(dense(s$components), s$fit$objective, tolerance = 1e-8)
  for (a in list(c(0,0,1), c(0,1/6,1/6), c(1,0,0))) {
    for (target in c(if (a[3] == 0) 0, sum(a*s$components), 2*sum(a*s$components))) {
      p <- .mfrm_mvgt_variance_at(s, target, a)
      expect_true(p$valid)
      expect_equal(sum(a*p$components), target, tolerance = 1e-10)
      expect_equal(p$fit$objective, dense(p$components), tolerance = 1e-8)
    }
  }
  expect_false(.mfrm_mvgt_variance_at(s, 0, c(0,0,1))$valid)
})

test_that("incomplete-source SEM preserves points, score units and proportional plans", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  g <- sem_profile_source()
  d <- mfrm_d_study(g, data.frame(Tasks = c(6,12,6)), c(V = -1, W = 1))
  original <- serialize(d, NULL)
  for (metric in c("RelativeSEM", "AbsoluteSEM")) {
    ci <- confint(d, parm = metric, method = "reml_profile", assumption = "normal")
    expect_true(all(ci$Available))
    expect_identical(ci$Estimate, subset(summary(d), Kind == "Composite")[[metric]])
    expect_equal(unlist(ci[2,c("Lower","Upper","IntervalCenter")]),
      unlist(ci[1,c("Lower","Upper","IntervalCenter")])/sqrt(2), ignore_attr = TRUE)
    expect_identical(attr(ci,"diagnostics")[[1]], attr(ci,"diagnostics")[[3]])
    expect_identical(serialize(d, NULL), original)
    scaled <- mfrm_d_study(g,d$design_grid,c(V=3,W=-3))
    cj <- confint(scaled,parm=metric,method="reml_profile",assumption="normal")
    expect_equal(as.matrix(cj[c("Lower","Upper","IntervalCenter")]),
      3*as.matrix(ci[c("Lower","Upper","IntervalCenter")]), tolerance=1e-6)
    expect_match(paste(attr(ci,"notes"),collapse=" "),"disconnected remote")
    expect_match(paste(attr(ci,"notes"),collapse=" "),"SEM units")
    expect_null(attr(ci,"reml_source")$objective)
    for (i in seq_len(nrow(ci))) {
      z <- attr(ci,"diagnostics")[[i]]
      expect_true(all(z$scan$Valid))
      expect_equal(sqrt(z$bounds)*attr(ci,"reml_source")$score_scale,
        c(ci$Lower[i],ci$Upper[i]),ignore_attr=TRUE)
      for (p in z$profiles) expect_equal(sum(p$weights*p$components),p$target,tolerance=1e-9)
    }
  }
})

test_that("SEM profiles compare with MLS and replay without fitting", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  x <- read.csv(system.file("extdata","mgenova-table12.csv",package="mfrmr"))
  d <- mfrm_d_study(mfrm_g_study(x,"V",rater=NULL),data.frame(Tasks=c(6,12)))
  ci <- confint(d,parm="AbsoluteSEM",method=c("mls","reml_profile"),assumption="normal")
  expect_identical(ci$MethodCode,rep(c("mls","reml_profile"),each=2))
  expect_true(all(ci$Available))
  expect_equal(ci$Estimate[1:2],ci$Estimate[3:4])
  f <- tempfile(fileext=".rds"); on.exit(unlink(f),add=TRUE)
  saveRDS(ci,f); expect_identical(readRDS(f),ci)
  local_mocked_bindings(.mfrm_mvgt_profile_source=function(...)stop("Unexpected refit"))
  forest <- plot(ci,draw=FALSE)
  expect_identical(plot_data(forest,"table")$Estimate,ci$Estimate)
  ribbon <- plot(d,type="sem",metric="AbsoluteSEM",intervals=ci[3:4,],draw=FALSE)
  expect_s3_class(ribbon,"mfrm_plot_data")
  if (requireNamespace("ggplot2",quietly=TRUE)) expect_s3_class(as_ggplot(ribbon),"ggplot")
})

test_that("target-directed SEM inversion retains failures and its local scope", {
  s <- list(valid=TRUE,components=c(1,1,1),sigma2=1,status="test",fit=list(objective=0))
  mock <- function(source,target,weights) list(target=target,weights=weights,
    lr=(target-1)^2/.04,valid=TRUE,check=list(kkt=0,error=0))
  local_mocked_bindings(.mfrm_mvgt_variance_at=mock)
  z <- .mfrm_mvgt_variance_set(s,c(0,0,1),.95)
  expect_equal(z$bounds,1+c(-1,1)*sqrt(qchisq(.95,1)*.04),tolerance=1e-7)
  local_mocked_bindings(.mfrm_mvgt_variance_at=function(source,target,weights) {
    p <- mock(source,target,weights);p$valid <- target >= .8;p
  })
  bad <- .mfrm_mvgt_variance_set(s,c(0,0,1),.95)
  expect_true(all(is.na(bad$bounds)))
  expect_true(is.na(bad$candidate_bounds[1]) && is.finite(bad$candidate_bounds[2]))
  expect_match(bad$status,"Unresolved")
  local_mocked_bindings(.mfrm_mvgt_variance_at=function(source,target,weights) {
    p <- mock(source,target,weights);p$lr <- 0;p
  })
  flat <- .mfrm_mvgt_variance_set(s,c(0,0,1),.95)
  expect_true(all(is.na(flat$bounds)))
  expect_true(all(is.na(flat$candidate_bounds)))
  local_mocked_bindings(.mfrm_mvgt_variance_at=function(source,target,weights) {
    p <- mock(source,target,weights);p$lr <- 100*((target-1)*(target-3))^2;p
  })
  disconnected <- .mfrm_mvgt_variance_set(s,c(0,0,1),.95)
  expect_true(all(is.finite(disconnected$bounds)))
  expect_lt(disconnected$bounds[2],2)
})
