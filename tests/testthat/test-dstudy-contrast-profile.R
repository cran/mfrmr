contrast_test_study <- function(nested = FALSE) {
  x <- if (nested) mfrmr_example_gtheory_nested else mfrmr_example_gtheory_crossed
  x <- subset(x, Person %in% head(unique(Person), if (nested) 6 else 10) &
    Task %in% head(unique(Task), 3))
  if (!nested) x <- subset(x, Rater %in% head(unique(Rater), 2))
  x <- x[(seq_len(nrow(x))*17) %% 11 < 8, ]
  g <- mfrm_g_study(x, "Content", method = "minque0",
    nesting = if (nested) c(Rater = "Task") else NULL)
  mfrm_d_study(g, data.frame(Raters = c(2,3,2,3), Tasks = c(3,2,3,2)))
}

test_that("coefficient geometry retains both roots and constraint derivatives", {
  for (ea in c(1e-8, .2, 10, 1e8)) for (eb in c(1e-8, .7, 20, 1e8)) {
    if (ea == eb) next
    limit <- tanh(log(ea/eb)/4)
    roots <- .mfrm_mvgt_contrast_person_roots(ea, eb, limit/2)
    expect_length(roots, 2)
    expect_equal(roots/(roots+ea)*(ea-eb)/(roots+eb), rep(limit/2,2), tolerance=1e-12)
    expect_equal(prod(roots)/(ea*eb), 1, tolerance=1e-12)
  }
  a <- c(0,0,0,1/2,1/3,0,1/6); b <- c(0,0,0,1/3,1/2,0,1/6)
  v <- c(1,.1,.2,.3,.4,.5,1)
  numerical <- vapply(seq_along(v), function(i) {
    plus <- minus <- v; plus[i] <- plus[i]+1e-5; minus[i] <- minus[i]-1e-5
    (.mfrm_mvgt_contrast_value(plus,a,b)-.mfrm_mvgt_contrast_value(minus,a,b))/2e-5
  }, numeric(1))
  expect_equal(.mfrm_mvgt_contrast_gradient(v,a,b), numerical, tolerance=1e-9)
  expect_equal(.mfrm_mvgt_contrast_range(a,b), c(-1,1)*tanh(log(1.5)/4))
})

test_that("contrast inversion preserves disconnected sets and unresolved searches", {
  a <- c(0,0,0,1/2,1/3,0,1/6); b <- c(0,0,0,1/3,1/2,0,1/6)
  v <- c(1,.2,.2,.2,.2,.2,1)
  source <- list(valid=TRUE,relative_variances=v,sigma2=1,components=v,status="Checked",
    fit=list(objective=0))
  critical <- qchisq(.95,1)
  run <- function(lr, aa=a, bb=b) {
    local_mocked_bindings(.mfrm_mvgt_contrast_profile_at=function(source,a,b,target,...)
      list(target=target,valid=is.finite(lr(target)),lr=lr(target)))
    .mfrm_mvgt_contrast_profile_set(source,aa,bb)
  }
  one <- run(function(x)(x/.02)^2)
  expect_true(one$valid)
  expect_equal(one$bounds,c(-1,1)*sqrt(critical)*.02,tolerance=1e-8)
  mode <- .mfrm_mvgt_contrast_range(a,b)[2]*.735
  two <- run(function(x)1e7*(x*(x-mode))^2)
  s <- sqrt(critical/1e7)
  expected <- sort(c((mode+c(-1,1)*sqrt(mode^2+4*s))/2,
    (mode+c(-1,1)*sqrt(mode^2-4*s))/2))
  expect_true(two$valid); expect_equal(nrow(two$sets),2)
  expect_equal(two$roots,expected,tolerance=1e-8)
  failed <- run(function(x)if(abs(x-.5*.mfrm_mvgt_contrast_range(a,b)[2])<.002) NA_real_ else (x/.02)^2)
  expect_false(failed$valid); expect_true(all(is.na(failed$bounds)))
  expect_true(all(is.finite(failed$candidate_bounds)))
  flat <- run(function(x)0)
  expect_false(flat$valid); expect_false(any(flat$tail_resolved))
  expect_true(all(is.na(flat$bounds)))
  boundary <- run(function(x)0,bb=a/2)
  expect_true(boundary$valid)
  expect_equal(boundary$bounds,.mfrm_mvgt_contrast_range(a,a/2))
  jump <- run(function(x)if(abs(x)<.03)0 else 100)
  expect_false(jump$valid); expect_true(any(grepl("cutoff",jump$failures)))
  source$valid <- FALSE
  expect_false(.mfrm_mvgt_contrast_profile_set(source,a,b)$valid)
  expect_true(all(is.na(.mfrm_mvgt_contrast_profile_set(source,a,a)$bounds)))
})

test_that("direct paired profiles keep the original points, source and plan context", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa"); skip_if_not_installed("nloptr")
  for (nested in c(FALSE,TRUE)) {
    d <- contrast_test_study(nested)
    metric <- if (nested) "AbsoluteSEM" else "G"
    original <- serialize(d,NULL)
    ci <- confint(d,parm=metric,reference=1,method="reml_profile",assumption="normal")
    point <- subset(summary(d),Kind=="Score")[[metric]]
    expect_equal(ci$Estimate,point[-1]-point[1])
    expect_true(all(ci$Available))
    expect_equal(ci$IntervalCenter[2],0)
    expect_equal(ci$Lower[2],0); expect_equal(ci$Upper[2],0)
    expect_equal(ci[c(1,3),c("Lower","Upper","IntervalCenter")][1,],
      ci[c(1,3),c("Lower","Upper","IntervalCenter")][2,],ignore_attr=TRUE)
    expect_identical(serialize(d,NULL),original)
    expect_null(attr(ci,"reml_source")$objective)
    expect_null(attr(ci,"base_interval"))
    expect_equal(attr(ci,"reference")$counts,d$design_grid[1,,drop=FALSE])
    expect_true(any(grepl("pointwise, not simultaneous",attr(ci,"notes"))))
    expect_true(all(grepl("not established",ci$Validation)))
    diagnostics <- attr(ci,"diagnostics")
    expect_identical(diagnostics[[1]],diagnostics[[3]])
    for (root in diagnostics[[1]]$roots)
      expect_equal(diagnostics[[1]]$profiles[[sprintf("%.17g",root)]]$lr,qchisq(.95,1),tolerance=1e-5)
    p <- plot(ci,draw=FALSE)
    expect_equal(plot_data(p,"table")$Estimate,ci$Estimate)
    expect_equal(plot_data(p,"table")$Lower,ci$Lower)
    saved <- tempfile(fileext=".rds")
    saveRDS(p,saved); expect_identical(readRDS(saved),p); unlink(saved)
    if (nested) {
      expect_equal(attr(ci,"reference")$children_per_parent,2)
      expect_false(any(grepl("comparisons keep a common",attr(ci,"notes"))))
    }
  }
})

test_that("near-zero profiles retain the equal-error branch as a start", {
  a <- c(0,0,0,1/2,1/3,0,1/6); b <- c(0,0,0,1/3,1/2,0,1/6)
  v <- c(1,.2,.2,.2,.2,.2,1)
  equal <- c(2,.1,.3,.4,.4,.2,1)
  source <- list(valid=TRUE,relative_variances=v,sigma2=1,components=v,status="Checked",
    fit=list(objective=0))
  attempts <- list(list(valid=TRUE,branch="zero Person variance",objective=0),
    list(valid=TRUE,branch="equal errors",objective=1,relative_variances=equal))
  null <- list(target=0,valid=TRUE,lr=0,attempts=attempts)
  starts <- list(); previous_seen <- FALSE
  local_mocked_bindings(.mfrm_mvgt_contrast_profile_at=function(source,a,b,target,
      starts_extra,previous=NULL) {
    if (target==0) return(null)
    starts[[length(starts)+1L]] <<- starts_extra
    if (!is.null(previous)) previous_seen <<- TRUE
    list(target=target,valid=TRUE,lr=(target/.02)^2)
  })
  grid <- .mfrm_mvgt_contrast_range(a,b)[1]+diff(.mfrm_mvgt_contrast_range(a,b))*.01
  x <- .mfrm_mvgt_contrast_profile_set(source,a,b,
    initial=list(null,list(target=grid,valid=FALSE,lr=NA_real_,attempts=list())))
  expect_true(x$valid); expect_true(previous_seen)
  expect_true(all(vapply(starts,function(s)identical(s,list(v,equal)),logical(1))))
})

test_that("SEM constraints use the full residual scale and the equal-error null", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa"); skip_if_not_installed("nloptr")
  d <- contrast_test_study()
  source <- .mfrm_mvgt_profile_source(d$gstudy,c(Content=1))
  a <- .mfrm_mvgt_profile_error(d,1,"AbsoluteSEM")
  b <- .mfrm_mvgt_profile_error(d,2,"AbsoluteSEM")
  target <- .03
  p <- .mfrm_mvgt_sem_contrast_profile_at(source,a,b,target)
  expect_true(p$valid)
  v <- p$best$relative_variances; sigma2 <- p$best$sigma2
  V <- sigma2*Reduce(`+`,Map(function(g,s)s*(outer(g,g,`==`)*1),source$groups,v))
  inv <- solve(V); one <- rowSums(inv); y <- source$y
  q <- drop(crossprod(y,inv%*%y))-sum(one*y)^2/sum(one)
  dense <- as.numeric(determinant(V,logarithm=TRUE)$modulus)+log(sum(one))+q+
    (length(y)-1)*log(2*pi)
  expect_equal(p$best$objective,dense,tolerance=1e-8)
  expect_equal(sqrt(sigma2)*.mfrm_mvgt_sem_contrast_value(v,a,b),target,tolerance=1e-8)
  expect_equal(sigma2,(target/.mfrm_mvgt_sem_contrast_value(v,a,b))^2,tolerance=1e-7)
  zero <- .mfrm_mvgt_sem_contrast_profile_at(source,a,b,0)
  expect_true(zero$valid)
  expect_equal(sum((a-b)*zero$best$relative_variances),0,tolerance=1e-8)
  expect_false(any(vapply(zero$attempts,function(x)identical(x$branch,"zero Person variance"),logical(1))))
  impossible <- .mfrm_mvgt_sem_contrast_profile_at(source,a,a/2,0)
  expect_true(impossible$infeasible); expect_identical(impossible$lr,Inf)
})
