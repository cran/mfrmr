profile_one_source <- function(incomplete = FALSE) {
  x <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
  if (incomplete) x <- x[(x$Person+x$Task) %% 3 != 0, ]
  mfrm_g_study(x, c("V", "W"), rater = NULL,
    method = if (incomplete) "minque0" else "anova")
}

test_that("REML target objectives agree with independent dense restricted likelihood", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  g <- profile_one_source(TRUE)
  source <- .mfrm_mvgt_profile_source(g, c(V = .4, W = .6))
  expect_true(source$valid)
  dense <- function(v) {
    V <- Reduce(`+`, Map(function(group, variance)
      variance * (outer(group, group, `==`)*1), source$groups, v))
    inv <- solve(V); one <- rowSums(inv); y <- source$y
    q <- drop(crossprod(y, inv %*% y))-sum(one*y)^2/sum(one)
    as.numeric(determinant(V, logarithm = TRUE)$modulus)+log(sum(one))+
      (length(y)-1)*(1+log(2*pi*q/(length(y)-1)))
  }
  expect_equal(source$fit$objective, dense(source$relative_variances), tolerance = 1e-8)
  for (rho in c(0, .25, .8)) {
    error <- c(0, 1/6, 1/6)
    p <- .mfrm_mvgt_profile_at(source, rho, error)
    expect_true(p$valid)
    expect_equal(p$fit$objective, dense(p$relative_variances), tolerance = 1e-8)
    v <- p$relative_variances
    expect_equal(v[1]/(v[1]+sum(error*v)), rho, tolerance = 1e-12)
    canonical <- .mfrm_mvgt_profile_at(source, rho/(6-5*rho), 6*error)
    expect_true(canonical$valid)
    expect_equal(canonical$fit$objective,p$fit$objective,tolerance=1e-8)
  }
  theta <- vapply(g$components, function(m) as.numeric(crossprod(c(.4,.6), m %*% c(.4,.6))), numeric(1))
  expect_equal(source$estimate/tail(source$estimate,1), theta/tail(theta,1))
  # An unused score with very different units must not erase the selected one.
  x <- g$data; x$V <- x$V*1e-100; x$W <- x$W*1e100
  scaled <- mfrm_g_study(x,c("V","W"),rater=NULL,method="minque0")
  ordinary <- .mfrm_mvgt_profile_source(g,c(V=1,W=0))
  extreme <- .mfrm_mvgt_profile_source(scaled,c(V=1,W=0))
  expect_true(extreme$valid)
  expect_equal(extreme$y,ordinary$y,tolerance=1e-12)
  expect_equal(extreme$fit$objective,ordinary$fit$objective,tolerance=1e-8)
})

test_that("invalid variance trials cannot corrupt later REML evaluations", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  source <- .mfrm_mvgt_profile_source(profile_one_source(TRUE),c(V=1,W=0))
  objective <- source$objective; par <- source$fit$par
  before <- objective(par,details=TRUE)
  for (bad in c(Inf,NaN,NA_real_,-1)) {
    trial <- par;trial[1] <- bad
    expect_identical(objective(trial),Inf)
    expect_identical(objective(trial,details=TRUE),list(objective=Inf,sigma2=NA_real_))
    expect_identical(objective(par,details=TRUE),before)
  }
  expect_identical(objective(c(par,0)),Inf)
  expect_identical(objective(par,details=TRUE),before)
  for (scale in c(1e14,1e100,1e300)) {
    trial <- par;trial[1] <- scale
    value <- suppressWarnings(objective(trial))
    expect_true(is.finite(value) || identical(value,Inf))
    expect_equal(objective(par,details=TRUE),before,tolerance=1e-12)
  }
})

test_that("incomplete profile intervals retain original points and numerical evidence", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  g <- profile_one_source(TRUE)
  d <- mfrm_d_study(g, data.frame(Tasks = c(6,12,6)), c(V = -1, W = 1))
  original <- serialize(d, NULL)
  ci <- confint(d, parm = "Phi", method = "reml_profile", assumption = "normal")
  point <- subset(summary(d), Kind == "Composite")
  expect_identical(ci$Estimate, point$Phi)
  expect_identical(ci$PointStatus, point$PhiStatus)
  expect_true(all(ci$Available))
  expect_true(all(ci$Lower <= ci$IntervalCenter & ci$Upper >= ci$IntervalCenter))
  expect_identical(attr(ci,"diagnostics")[[1]], attr(ci,"diagnostics")[[3]])
  expect_null(attr(ci,"reml_source")$objective)
  expect_true(any(grepl("not.*automatic", attr(ci,"notes"))))
  expect_identical(serialize(d, NULL), original)
  for (z in attr(ci, "diagnostics")) {
    expect_true(all(z$scan$Valid))
    for (root in z$roots) {
      p <- z$profiles[[sprintf("%.17g",root)]]
      expect_equal(p$lr, qchisq(.95,1), tolerance = 1e-5)
    }
  }
  expect_error(confint(d, assumption="normal"), "balanced ANOVA")
  expect_error(confint(d, method="reml_profile", assumption="normal", reference=0), "reference")
  d$gstudy$data <- NULL
  expect_error(confint(d, method="reml_profile", assumption="normal"), "used ratings")
})

test_that("profile intervals respect crossed and nested facet definitions", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  x <- mfrmr_example_gtheory_crossed
  x <- subset(x, Person %in% head(unique(Person),10) &
    Rater %in% head(unique(Rater),2) & Task %in% head(unique(Task),3))
  x <- x[(seq_len(nrow(x))*17) %% 11 < 8, ]
  g <- mfrm_g_study(x, c("Content","Organization"), method="minque0")
  d <- mfrm_d_study(g, data.frame(Raters=2,Tasks=3), c(Content=.4,Organization=.6))
  ci <- confint(d, parm="Phi", method="reml_profile", assumption="normal")
  y <- x[nrow(x):1, ]; y$Content <- 100-3*y$Content; y$Organization <- 100-3*y$Organization
  gy <- mfrm_g_study(y,c("Content","Organization"),
    facets=c(Occasion="Task",Judge="Rater"), method="minque0")
  dy <- mfrm_d_study(gy,data.frame(Occasion=3,Judge=2),c(Content=.4,Organization=.6))
  other <- confint(dy,parm="Phi",method="reml_profile",assumption="normal")
  expect_true(all(ci$Available & other$Available))
  expect_equal(other[c("Lower","Upper","IntervalCenter")],
    ci[c("Lower","Upper","IntervalCenter")], ignore_attr=TRUE,tolerance=1e-6)

  z <- mfrmr_example_gtheory_nested
  z <- subset(z, Person %in% head(unique(Person),6) & Task %in% head(unique(Task),3))
  z <- z[(seq_len(nrow(z))*17) %% 11 < 8, ]
  gn <- mfrm_g_study(z,"Content",nesting=c(Rater="Task"),method="minque0")
  dn <- mfrm_d_study(gn,data.frame(Raters=2,Tasks=4))
  nci <- confint(dn,parm="Phi",method="reml_profile",assumption="normal")
  expect_equal(unname(attr(nci,"diagnostics")[[1]]$error_coefficients),c(0,1/4,1/8,1/4,1/8))
  expect_length(attr(nci,"reml_source")$groups,5)
  expect_true(nci$Available)
  expect_equal(attr(nci,"nesting_context")$future$ChildrenPerParent,2)
  expect_equal(attr(nci,"nesting_context")$future$Parents,4)
})

test_that("profile difference images include extrema and preserve failed base sets", {
  d <- mfrm_d_study(profile_one_source(),data.frame(Tasks=c(6,12,3,6)))
  calls <- 0L
  base <- NULL
  local_mocked_bindings(.mfrm_mvgt_profile_table=function(d,tab,weights,metric,level) {
    calls <<- calls+1L
    expect_equal(d$design_grid$Tasks,1)
    x <- data.frame(IntervalCenter=.2,Lower=.01,Upper=.8,
      CandidateLower=.01,CandidateUpper=.8,Available=TRUE,Status="Computed",
      Validation="Experimental")
    if (!is.null(base)) x <- base
    attr(x,"reml_source") <- list(estimate=c(1,1,1))
    x
  })
  ci <- confint(d,reference=1,method="reml_profile",assumption="normal")
  expect_identical(calls,1L)
  expect_equal(ci$Scenario,c(2,3,4))
  f <- function(t,b,a) b*t/(1+b*t)-a*t/(1+a*t)
  ts <- c(.01/.99,.8/.2,1/sqrt(6*12),1/sqrt(6*3))
  for (i in 1:2) {
    vals <- f(ts,ci$Tasks[i],6)
    expect_equal(c(ci$Lower[i],ci$Upper[i]),range(vals),tolerance=1e-12)
  }
  expect_equal(c(ci$Lower[3],ci$Upper[3]),c(0,0))
  expect_gt(ci$Upper[1],max(f(ts[1:2],12,6)))
  expect_equal(ci$IntervalCenter,f(.2/.8,ci$Tasks,6))
  base <- data.frame(IntervalCenter=.2,Lower=NA_real_,Upper=NA_real_,
    CandidateLower=0,CandidateUpper=1,Available=FALSE,Status="Unresolved upper tail",
    Validation="Experimental")
  failed <- confint(d,reference=1,method="reml_profile",assumption="normal")
  expect_false(any(failed$Available))
  expect_true(all(is.na(failed$Lower) & is.na(failed$Upper)))
  expect_true(all(is.finite(failed$CandidateLower) & is.finite(failed$CandidateUpper)))
  expect_match(failed$Status[3],"Unresolved")
  base$CandidateUpper <- NA_real_
  partial <- confint(d,reference=1,method="reml_profile",assumption="normal")
  expect_true(all(is.na(partial$CandidateLower) & is.na(partial$CandidateUpper)))
  expect_equal(attr(partial,"base_interval")$CandidateLower,0)
})

test_that("incomplete profile comparisons preserve one shared fit and original points", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  d <- mfrm_d_study(profile_one_source(TRUE),
    data.frame(Tasks=c(6,12,3,6)),c(V=-1,W=1))
  original <- subset(summary(d),Kind=="Composite")
  for (metric in c("Phi","AbsoluteSEM")) {
    ci <- confint(d,parm=metric,reference=1,method="reml_profile",assumption="normal")
    base <- attr(ci,"base_interval")
    expect_true(all(ci$Available))
    expect_equal(attr(base,"design_grid")$Tasks,1)
    expect_equal(attr(ci,"source_data"),d$gstudy$data)
    expect_identical(ci$Estimate,original[[metric]][-1]-original[[metric]][1])
    if (metric=="Phi") {
      rho <- seq(base$Lower,base$Upper,length.out=20001)
      for (i in 1:2) {
        delta <- ci$Tasks[i]*rho/(1-rho+ci$Tasks[i]*rho)-6*rho/(1-rho+6*rho)
        expect_equal(c(ci$Lower[i],ci$Upper[i]),range(delta),tolerance=1e-7)
      }
    } else {
      for (i in 1:2) expect_equal(c(ci$Lower[i],ci$Upper[i]),
        sort(c(base$Lower,base$Upper)*(1/sqrt(ci$Tasks[i])-1/sqrt(6))))
    }
    expect_equal(c(ci$Lower[3],ci$Upper[3]),c(0,0))
    expect_true(all(ci$Lower <= ci$IntervalCenter & ci$Upper >= ci$IntervalCenter))
    # Reversing the reference must reuse the same population base target.
    reverse <- (function() {
      local_mocked_bindings(.mfrm_mvgt_profile_table=function(...)base)
      confint(d,parm=metric,reference=2,method="reml_profile",assumption="normal")
    })()
    expect_equal(c(reverse$Lower[1],reverse$Upper[1]),-c(ci$Upper[1],ci$Lower[1]))
    expect_equal(reverse$IntervalCenter[1],-ci$IntervalCenter[1])
  }
})

test_that("method comparison and saved plots use existing point and interval semantics", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  d <- mfrm_d_study(profile_one_source(),data.frame(Tasks=c(6,12)))
  ci <- confint(d,method=c("exact_f","reml_profile"),assumption="normal")
  expect_identical(ci$MethodCode,rep(c("exact_f","reml_profile"),each=2))
  expect_length(attr(ci,"diagnostics"),4)
  expect_null(attr(ci,"diagnostics")[[1]])
  expect_true(all(is.na(ci$IntervalCenter[1:2])))
  expect_true(all(is.finite(ci$IntervalCenter[3:4])))
  expect_equal(ci$Estimate[1:2],ci$Estimate[3:4])
  tmp <- tempfile(fileext=".rds"); on.exit(unlink(tmp),add=TRUE)
  saveRDS(ci,tmp); restored <- readRDS(tmp)
  expect_identical(restored,ci)
  local_mocked_bindings(.mfrm_mvgt_profile_source=function(...)stop("Unexpected refit"))
  p <- plot(restored,draw=FALSE)
  expect_identical(plot_data(p,"table")$Estimate,ci$Estimate)
  single <- ci[3:4, ]
  curve <- plot(d,metric="G",score="V",intervals=single,draw=FALSE)
  expect_true(inherits(curve,"mfrm_plot_data"))
  if (requireNamespace("ggplot2",quietly=TRUE)) expect_s3_class(as_ggplot(p),"ggplot")
  expect_error(plot(d,metric="G",score="V",intervals=ci,draw=FALSE),"one interval method")
  altered <- d; altered$gstudy$data$V[1] <- altered$gstudy$data$V[1]+1
  expect_error(plot(altered,metric="G",score="V",intervals=single,draw=FALSE),"source ratings")
})

test_that("profile comparisons retain unequal nested assignments and count meanings", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  x <- subset(mfrmr_example_gtheory_nested,
    Person %in% head(unique(Person),6) & Task %in% head(unique(Task),3))
  first <- x$Task[1]
  x <- x[!(x$Task==first & x$Rater==tail(unique(x$Rater[x$Task==first]),1)), ]
  x <- x[(seq_len(nrow(x))*17) %% 11 < 8, ]
  g <- mfrm_g_study(x,"Content",facets=c(Occasion="Task",Judge="Rater"),
    nesting=c(Judge="Occasion"),method="minque0")
  d <- mfrm_d_study(g,data.frame(Occasion=c(3,6,3),Judge=3))
  ci <- confint(d,parm="Phi",reference=1,method="reml_profile",assumption="normal")
  base <- attr(ci,"base_interval")
  expect_true(all(ci$Available))
  expect_equal(attr(base,"design_grid"),data.frame(Occasion=1,Judge=3))
  expect_identical(attr(ci,"source_data"),g$data)
  expect_length(attr(base,"reml_source")$groups,5)
  expect_equal(unname(attr(base,"diagnostics")[[1]]$error_coefficients),c(0,1,1/3,1,1/3))
  expect_gt(length(unique(g$design$child_counts)),1)
  expect_identical(attr(ci,"nesting_context")$source_children_per_parent,g$design$child_counts)
  expect_equal(attr(ci,"reference")$count_facet,"Occasion")
  expect_equal(attr(ci,"reference")$children_per_parent,3)
  expect_equal(c(ci$Lower[2],ci$Upper[2]),c(0,0))
  restored <- unserialize(serialize(ci,NULL))
  expect_identical(restored,ci)
  local_mocked_bindings(.mfrm_mvgt_profile_source=function(...)stop("Unexpected refit"))
  p <- plot(restored,draw=FALSE)
  expect_identical(plot_data(p,"table")$Estimate,ci$Estimate)
  if (requireNamespace("ggplot2",quietly=TRUE)) expect_s3_class(as_ggplot(p),"ggplot")
})

test_that("method vectors preserve the reference in every profile comparison", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  d <- mfrm_d_study(profile_one_source(),data.frame(Tasks=c(6,12,3)))
  ci <- confint(d,reference=1,method=c("exact_f","reml_profile"),assumption="normal")
  expect_equal(ci$Scenario,rep(c(2,3),2))
  expect_identical(ci$Metric,rep("G difference",4))
  expect_equal(ci$Estimate[1:2],ci$Estimate[3:4])
  expect_length(attr(ci,"diagnostics"),4)
  expect_true(all(ci$Available))
  expect_equal(attr(ci,"base_interval")$Metric,"G")
  expect_equal(attr(ci,"reference")$scenario,1)
  expect_s3_class(plot(ci,draw=FALSE),"mfrm_plot_data")
})

test_that("profile failures retain candidates without inventing accepted bounds", {
  source <- list(valid=TRUE,relative_variances=c(1,1,1),status="test",fit=list(objective=0))
  error <- c(0,0,3)
  local_mocked_bindings(.mfrm_mvgt_profile_at=function(source,rho,error)
    list(target=rho,lr=0,valid=TRUE,check=list(kkt=0,error=0)))
  unlimited <- .mfrm_mvgt_profile_set(source,error,.95)
  expect_true(all(is.na(unlimited$bounds)))
  expect_equal(unlimited$candidate_bounds,c(0,1))
  expect_match(unlimited$status,"Unresolved")
  local_mocked_bindings(.mfrm_mvgt_profile_at=function(source,rho,error)
    list(target=rho,lr=1000*((rho-.25)*(rho-.75))^2,valid=TRUE,
      check=list(kkt=0,error=0)))
  separate <- .mfrm_mvgt_profile_set(source,error,.95)
  expect_equal(nrow(separate$sets),2)
  expect_equal(separate$bounds,c(min(separate$sets),max(separate$sets)))
  expect_match(separate$status,"disconnected")
})

test_that("boundary derivatives and degenerate data do not become zero uncertainty", {
  fn <- function(v) (v[1]-.4)^2+3*v[2]+v[2]^2
  check <- .mfrm_mvgt_profile_gradient(c(.4,0),fn)
  expect_equal(check$gradient,c(0,3),tolerance=1e-7)
  expect_true(check$valid)
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  g <- profile_one_source(); g$data$V[] <- 1
  g <- mfrm_g_study(g$data,"V",rater=NULL)
  ci <- confint(mfrm_d_study(g),method="reml_profile",assumption="normal")
  expect_false(ci$Available)
  expect_match(ci$Status,"Degenerate")
  expect_true(is.na(ci$Lower) && is.na(ci$Upper))
  f <- .mfrm_mvgt_profile_source
  environment(f) <- list2env(list(requireNamespace=function(...)FALSE),parent=environment(f))
  expect_error(f(g,c(V=1)),"optional lme4")
})

test_that("a rounded REML coefficient center explains why its interval is withheld", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  x <- profile_one_source()$data
  x$Score <- x$Person + 1e-4*x$V
  g <- mfrm_g_study(x, "Score", rater = NULL)
  d <- mfrm_d_study(g, data.frame(Tasks = 2e9))
  before <- serialize(d, NULL)
  ci <- confint(d, parm = "G", method = "reml_profile", assumption = "normal")
  expect_true(attr(ci, "reml_source")$valid)
  expect_identical(ci$Estimate, d$coefficients$G)
  expect_identical(ci$IntervalCenter, 1)
  expect_false(ci$Available)
  expect_true(all(is.na(ci[c("Lower", "Upper", "CandidateLower", "CandidateUpper")])))
  expect_match(ci$Status, "center rounds to one.*numeric resolution")
  expect_identical(attr(ci, "diagnostics")[[1]]$status, ci$Status)
  expect_match(paste(attr(ci, "notes"), collapse = " "), "does not establish perfect reliability")
  expect_identical(serialize(d, NULL), before)
  expect_identical(plot_data(plot(d, metric = "G", score = "Score", intervals = ci,
    draw = FALSE), "intervals")[[1L]], ci)
})
