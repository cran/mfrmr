test_that("profile nuisance optimization preserves a dependent sum-zero slope", {
  # The last expanded log slope is -(p1+p2), not a separate free coordinate.
  start <- c(.2, -.3, .5)
  contrast <- c(-1, -1, 0)
  ev <- list(value = function(p) sum((p - start)^2) / 2,
    gradient = function(p) p - start)
  v <- .8
  ans <- mfrmr:::mfrm_gpcm_profile_nuisance(ev, ev, start, contrast, v, 100)
  expected <- start + contrast * (v - sum(contrast * start)) / sum(contrast^2)
  expect_true(ans$check$Passed)
  expect_equal(ans$parameters, expected, tolerance = 1e-7)
  expect_equal(sum(ans$parameters * contrast), v, tolerance = 1e-12)
  map <- mfrmr:::mfrm_gpcm_profile_constraint(start, contrast, v)
  expect_equal(drop(crossprod(contrast, map$basis)), c(0, 0))
})

test_that("profile endpoints recover an independent quadratic likelihood interval", {
  center <- c(.2, -.1, .3)
  precision <- matrix(c(4,1,0,1,3,1,0,1,2), 3)
  ev <- list(value = function(p) drop(crossprod(p-center, precision %*% (p-center))) / 2,
    gradient = function(p) drop(precision %*% (p-center)))
  a <- c(-1,-1,0)
  ans <- mfrmr:::mfrm_gpcm_profile_search(ev, ev, center, a, .95, maxit=100)
  se <- sqrt(drop(a %*% solve(precision) %*% a))
  expect_equal(ans$endpoints$Bound, exp(sum(a*center) + c(-1,1)*qnorm(.975)*se), tolerance=1e-6)
  expect_true(all(ans$endpoints$Status == "computed"))
  expect_true(all(ans$profile$Passed))
  expect_length(ans$attempts, nrow(ans$profile))
  bad <- ev; bad$value <- function(p) ev$value(p) + .01
  failed <- mfrmr:::mfrm_gpcm_profile_search(ev,bad,center,a,.95,maxit=100)
  expect_true(all(is.na(failed$endpoints$Bound)))
  expect_false(any(failed$profile$Passed))
  limited <- mfrmr:::mfrm_gpcm_profile_search(ev,ev,center,a,.95,maxit=100,
    max_steps=1,initial_step=.001)
  expect_true(all(is.na(limited$endpoints$Bound)))
  expect_true(all(limited$endpoints$Status == "crossing_not_found"))
  expect_true(all(is.finite(limited$profile$NLL)))
})

test_that("profile points reject stationary nuisance maxima at either integration order", {
  # The source has positive information, but y=0 becomes a nuisance maximum
  # when cos(x)<0. Both retained and neutral starts stay at y=0. The true
  # profiled NLL there is x^2/2 - cos(x)^2/4, not x^2/2.
  ev <- list(value = function(p) p[1]^2/2 + cos(p[1])*p[2]^2/2 + p[2]^4/4,
    gradient = function(p) c(p[1]-sin(p[1])*p[2]^2/2, cos(p[1])*p[2]+p[2]^3))
  failed <- mfrm_gpcm_profile_nuisance(ev, ev, c(0,0), c(1,0), 2, 100)
  expect_identical(failed$check$OptimizerCode, 0L)
  expect_equal(failed$check$MaxGradient, 0)
  # Central differentiation of the cubic gradient adds h^2 (h = .001).
  expect_lte(abs(failed$check$NuisanceMinEigenvalue - cos(2)), 1.01e-6)
  expect_false(failed$check$Passed)
  expect_match(failed$check$Detail, "curvature.*not positive definite")
  # This is local qualification, not a proof that the interval does not exist.
  ans <- mfrm_gpcm_profile_search(ev, ev, c(0,0), c(1,0), .95)
  expect_true(all(is.na(ans$endpoints$Bound)))
  expect_true(any(grepl("curvature", ans$profile$Detail)))
  root <- uniroot(function(x) x^2 - cos(x)^2/2 - qchisq(.95,1), c(pi/2, 3))$root
  expect_gt(root, qnorm(.975))
  recovered <- mfrm_gpcm_profile_nuisance(ev, ev, c(0,1), c(1,0), 2, 100)
  expect_true(recovered$check$Passed)
  expect_equal(recovered$check$NLL, 2 - cos(2)^2/4, tolerance = 1e-8)
  expect_gt(recovered$check$NuisanceMinEigenvalue, 0)
  # Pointwise value/gradient agreement alone cannot qualify the reference grid.
  positive <- list(value = function(p) sum(p^2)/2, gradient = function(p) p)
  negative <- list(value = function(p) (p[1]^2-p[2]^2)/2,
    gradient = function(p) c(p[1], -p[2]))
  reference <- mfrm_gpcm_profile_nuisance(positive, negative, c(0,0), c(1,0), 0, 100)
  expect_equal(reference$check$NLLChange, 0)
  expect_equal(reference$check$GradientChange, 0)
  expect_gt(reference$check$NuisanceMinEigenvalue, 0)
  expect_lt(reference$check$ReferenceMinEigenvalue, 0)
  expect_false(reference$check$Passed)
  flat <- list(value = function(p) p[1]^2/2, gradient = function(p) c(p[1], 0))
  expect_false(mfrm_gpcm_profile_nuisance(flat, flat, c(0,0), c(1,0), 0, 100)$check$Passed)
  unavailable <- positive
  unavailable$gradient <- function(p) { if (p[2] != 0) stop("unavailable derivative"); p }
  broken <- mfrm_gpcm_profile_nuisance(positive, unavailable, c(0,0), c(1,0), 0, 100)
  expect_false(broken$check$Passed)
  expect_match(broken$check$Detail, "curvature could not be verified: unavailable derivative")
})

test_that("profile search scaling preserves total-likelihood cutoffs", {
  n <- 240
  center <- c(.2, -.3)
  ev <- list(value = function(p) n * sum((p-center)^2) / 2,
    gradient = function(p) n * (p-center), optimization_scale = n)
  ans <- mfrm_gpcm_profile_search(ev, ev, center, c(1,0), .95)
  expected <- exp(center[1] + c(-1,1) * qnorm(.975)/sqrt(n))
  expect_equal(ans$endpoints$Bound, expected, tolerance = 1e-6)
  for(bound in ans$endpoints$Bound) {
    at <- which.min(abs(ans$profile$LogSlope-log(bound)))
    expect_lte(abs(ans$profile$LR[at]-qchisq(.95,1)), 1e-4)
  }
})

test_that("a failed outward profile trial contracts without hiding failures", {
  ev <- list(value = function(p) 50 * sum(p^2), gradient = function(p) 100*p)
  reference <- ev
  reference$value <- function(p) ev$value(p) + if (abs(p[1]) > .3) .01 else 0
  ans <- mfrm_gpcm_profile_search(ev, reference, c(0,0), c(1,0), .95)
  expect_true(all(ans$endpoints$Status == "computed"))
  expect_equal(ans$endpoints$Bound, exp(c(-1,1)*qnorm(.975)/10), tolerance = 1e-6)
  failed <- ans$profile[!ans$profile$Passed,]
  expect_equal(sort(failed$LogSlope), c(-.5,.5))
  expect_length(ans$attempts,nrow(ans$profile))
  limited <- mfrm_gpcm_profile_search(ev, reference, c(0,0), c(1,0), .95,max_steps=1L)
  expect_true(all(is.na(limited$endpoints$Bound)))
  expect_true(all(limited$endpoints$Status == "numerical_or_start_failure"))
  expect_equal(nrow(limited$profile),3L) # Center plus one trial on each side.
  # Integration remains inadequate before the true endpoints: do not rescue.
  reference$value <- function(p) ev$value(p) + if (abs(p[1]) > .15) .01 else 0
  blocked <- mfrm_gpcm_profile_search(ev, reference, c(0,0), c(1,0), .95)
  expect_true(all(is.na(blocked$endpoints$Bound)))
  expect_true(all(blocked$endpoints$Status == "numerical_or_start_failure"))
  # A failed evaluation inside a root bracket is not skipped.
  reference$value <- function(p) ev$value(p) + if (abs(p[1]) > .01 && abs(p[1]) < .3) .01 else 0
  gap <- mfrm_gpcm_profile_search(ev, reference, c(0,0), c(1,0), .95)
  expect_true(all(is.na(gap$endpoints$Bound)))
  expect_true(all(gap$endpoints$Status == "endpoint_search_failed"))
})

test_that("profile curvature refinement also checks the higher-order nuisance score", {
  # The optimizer has stopped just inside the low-grid tolerance, but its
  # reference gradient is just outside it. A real Newton step resolves this.
  ev <- list(value = function(p) sum(p^2) / 2, gradient = function(p) p)
  hi <- list(value = function(p) ev$value(p) + 3e-6 * p[2],
    gradient = function(p) p + c(0, 3e-6))
  local_mocked_bindings(optim = function(...) list(par = 9.8e-5, convergence = 0L), .package = "mfrmr")
  ans <- mfrm_gpcm_profile_nuisance(ev, hi, c(0, 0), c(1, 0), 0, 100)
  expect_true(ans$check$Passed)
  expect_lt(abs(ans$parameters[2]), 1e-10)
  expect_identical(tail(ans$stages, 1)[[1]]$stage, "curvature_polish")
  # A substantive integration error cannot be repaired by optimizing the low grid.
  hi$gradient <- function(p) p + c(0, .001)
  failed <- mfrm_gpcm_profile_nuisance(ev, hi, c(0, 0), c(1, 0), 0, 100)
  expect_false(failed$check$Passed)
  expect_gt(failed$check$ReferenceGradient, 1e-4)
})

test_that("better source likelihoods and numerical errors are retained, not promoted", {
  ev <- list(value=function(p)sum(p^2)/2, gradient=function(p)p)
  bad_start <- mfrmr:::mfrm_gpcm_profile_search(ev,ev,c(0,2),c(1,0),.95,maxit=50)
  expect_true(all(bad_start$endpoints$Status == 'source_likelihood_improved'))
  expect_true(all(is.na(bad_start$endpoints$Bound)))
  expect_match(bad_start$profile$Detail[1], 'improved')
  # A better likelihood discovered outward invalidates an already found bound.
  escaped <- list(value=function(p) 50*sum(p^2)-if(p[1]>.4) 100 else 0,
    gradient=function(p) 100*p)
  outward <- mfrm_gpcm_profile_search(escaped,escaped,c(0,0),c(1,0),.95)
  expect_true(all(outward$endpoints$Status == "source_likelihood_improved"))
  expect_true(all(is.na(outward$endpoints$Bound)))
  expect_true(any(outward$profile$LR < 0))
  failed <- mfrmr:::mfrm_gpcm_profile_nuisance(
    list(value=function(p)stop('unavailable'),gradient=function(p)p),ev,
    c(0,0),c(1,0),1,10)
  expect_false(failed$check$Passed)
  expect_null(failed$parameters)
  expect_true(all(vapply(failed$stages,function(s)grepl('unavailable',s$error),TRUE)))
})

test_that("profile API validates its target and retains unresolved endpoint reasons", {
  fit <- readRDS(test_path('fixtures','mfrm-conditional-scoring-gpcm.rds'))$fit
  slope <- fit$config$gpcm_spec$levels[1]
  expect_error(confint(fit, parm = "slopes",method='profile'),'one exact slope-level')
  expect_error(confint(fit, parm = "slopes",method='profile',slope='unknown'),'one exact slope-level')
  expect_error(confint(fit, parm = "slopes",method='profile',slope=slope,scale='standardized'),'one relative slope')
  expect_error(confint(fit, parm = "slopes",slope=slope),'require method')
  expect_error(confint(fit, parm = "slopes",method='profile',slope=slope,profile_control=list(tolerance=1)),'accepts only')
  expect_error(confint(fit, parm = "slopes",method='profile',slope=slope,profile_control=list(maxit=NULL)),'positive finite')
  changed <- fit; changed$config$population_spec$design_columns <- 'X'
  expect_error(confint(changed, parm = "slopes",method='profile',slope=slope),'intercept-only')
  saved <- readRDS(test_path('fixtures','gpcm-profile-separate.rds'))
  result <- attr(saved,'profile')
  local_mocked_bindings(mfrm_gpcm_profile_search=function(...) result,.package='mfrmr')
  # Version 2 compares stored values without ALTREP materialization state.
  before <- serialize(fit, NULL, version = 2)
  out <- confint(fit, parm = "slopes",method='profile',slope=slope)
  expect_equal(as.vector(out),as.vector(saved))
  expect_match(attr(out,'target'),'Criterion C1')
  expect_identical(serialize(fit, NULL, version = 2), before)
  # The saved search is a test double: this checks only the Wald comparator,
  # not profile endpoints at a different likelihood-ratio cutoff.
  level <- 1 - .Machine$double.eps / 2
  high <- confint(fit, parm = "slopes", method='profile', slope=slope, level=level)
  wald <- attr(out, 'wald'); high_wald <- attr(high, 'wald')
  se <- (log(wald$Upper)-log(wald$Lower))/(2*sqrt(qchisq(.95, 1)))
  expect_equal(log(high_wald$Lower), log(wald$Estimate)-sqrt(qchisq(level, 1))*se)
  expect_equal(log(high_wald$Upper), log(wald$Estimate)+sqrt(qchisq(level, 1))*se)
  result$endpoints$Status[1] <- 'crossing_not_found'; result$endpoints$Bound[1] <- NA_real_
  expect_warning(unresolved <- confint(fit, parm = "slopes",method='profile',slope=slope),'no cutoff crossing')
  expect_true(is.na(unresolved[1,1])); expect_true(is.finite(unresolved[1,2]))
  expect_false(attr(unresolved,'diagnostics')$CIEligible)
  expect_false(grepl('crossing_not_found',attr(unresolved,'diagnostics')$InferenceReview,fixed=TRUE))
  expect_identical(attr(unresolved,'profile')$endpoints$Status[1],'crossing_not_found')
})

test_that("saved profile figures, reports and archives never restart the search", {
  fit <- readRDS(test_path('fixtures','mfrm-conditional-scoring-gpcm.rds'))$fit
  ci <- readRDS(test_path('fixtures','gpcm-profile-separate.rds'))
  local_mocked_bindings(mfrm_gpcm_profile_search=function(...) stop('unexpected profile'),
    mfrm_gpcm_inference=function(...) stop('unexpected covariance'),
    fit_mfrm=function(...) stop('unexpected fit'),.package='mfrmr')
  p <- as_ggplot(ci,title=NULL,subtitle=NULL,caption=NULL)
  expect_null(p$labels$title);expect_null(p$labels$subtitle);expect_null(p$labels$caption)
  expect_identical(plot_data(p,'endpoints'),attr(ci,'profile')$endpoints)
  payload <- attr(p,'mfrmr_plot_data')
  expect_identical(plot_data(as_ggplot(payload)),plot_data(p))
  expect_s3_class(as_ggplot(ci,type='interval'),'ggplot')
  expect_error(as_ggplot(payload,title='change'),'must be empty')
  expect_error(plot(ci,type='profile',reference=0,draw=FALSE),'positive')
  gap <- ci; pr <- attr(gap,'profile'); mid <- floor(nrow(pr$profile)/2)
  pr$profile$Passed[mid] <- FALSE; pr$profile$Detail[mid] <- 'Failed integration check'
  attr(gap,'profile') <- pr
  gp <- plot(gap,draw=FALSE); tab <- plot_data(gp,'table')
  expect_identical(tab$Display[mid],'Unresolved')
  expect_gt(length(unique(tab$Segment[tab$Passed])),1)
  expect_no_error(ggplot2::ggplot_build(gp))
  res <- mfrm_results(fit,intervals=list(profile=ci),include=c('fit','plots'),compute='never')
  expect_identical(res$tables$gpcm_profile_profile,attr(ci,'profile')$profile)
  expect_identical(mfrm_report(res)$tables$gpcm_profile_profile_endpoints,attr(ci,'profile')$endpoints)
  folder <- tempfile();on.exit(unlink(folder,recursive=TRUE),add=TRUE)
  suppressWarnings(export_mfrm_results(res,output_dir=folder,include=c('tables','replay','report','plots'),
    preset='starter',acknowledge_sensitive=TRUE))
  paths <- list.files(folder,recursive=TRUE,full.names=TRUE)
  archive <- readRDS(paths[grepl('[.]rds$',paths)][1])
  expect_identical(archive$gpcm_inference$profile,ci)
  expect_true(any(grepl('profile_checks.*[.]csv$',paths)))
  expect_true(any(grepl('plot_gpcm_profile[.]png$',paths)))
  expect_identical(plot_data(archive,type='gpcm_profile',component='endpoints'),attr(ci,'profile')$endpoints)
})
