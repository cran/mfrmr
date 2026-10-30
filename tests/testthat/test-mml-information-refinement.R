test_that("weak positive information is inverted without an eigenvalue floor", {
  h <- diag(c(1e-10, 2))
  fn <- function(p) sum(p * (h %*% p))/2
  gr <- function(p) drop(h %*% p)
  out <- mfrm_refine_mml_information(c(0,0),fn,gr,h)
  expect_true(out$review$Verified)
  expect_equal(out$cov,diag(c(1e10,.5)),tolerance=1e-10)
  expect_lt(out$review$RelativeChange,1e-8)
  expect_match(out$review$Detail,"Ill-conditioned")
  expect_false(isTRUE(all.equal(out$cov,invert_information_matrix(h)$cov)))
  # A change of units must give the corresponding covariance, not a ridge.
  scaling <- diag(c(1e4,.01))
  hs <- t(scaling)%*%h%*%scaling
  transformed <- mfrm_refine_mml_information(c(0,0),
    function(z) fn(drop(scaling%*%z)),function(z) drop(t(scaling)%*%gr(drop(scaling%*%z))),hs)
  expect_true(transformed$review$Verified)
  expect_equal(scaling%*%transformed$cov%*%t(scaling),out$cov,tolerance=1e-8)
})

test_that("unresolved information cannot become an accepted regularized inverse", {
  for(h in list(diag(c(0,1)),diag(c(-1e-8,1)))) {
    out <- mfrm_refine_mml_information(c(0,0),function(p) sum(p*(h%*%p))/2,
      function(p) drop(h%*%p),diag(c(1e-10,1)))
    expect_false(out$review$Verified)
    expect_null(out$cov)
    expect_match(out$review$Detail,"refinement failed")
  }
  # A raw gradient of 1e-6 is not small in the weak curvature direction.
  h <- diag(c(1e-10,1))
  out <- mfrm_refine_mml_information(c(1e4,0),function(p) sum(p*(h%*%p))/2,
    function(p) drop(h%*%p),h)
  expect_false(out$review$Verified)
  expect_gt(out$review$CurvatureScaledGradient,1e-4)
  expect_null(out$cov)
  expect_match(out$review$Detail,"standardized Newton displacement 0.1 (limit 1e-4)",fixed=TRUE)
  expect_match(out$review$Detail,"inverse residual",fixed=TRUE)
  bad <- mfrm_refine_mml_information(c(0,0),function(p) stop("domain"),
    function(p) stop("domain"),h)
  expect_false(bad$review$Verified)
  expect_null(bad$cov)
  expect_match(bad$review$Detail,"domain")
})

test_that("two refinement levels expose unstable weak-direction curvature", {
  h <- diag(c(1e-12,1))
  # Derivative cancellation may be small entrywise but large relative to the
  # weak direction. A deliberate oscillatory perturbation tests that failure.
  gr <- function(p) c(1e-12*p[1]+1e-15*sin(1e6*p[1]),p[2])
  out <- mfrm_refine_mml_information(c(0,0),function(p)
    5e-13*p[1]^2 - 1e-21*cos(1e6*p[1]) + p[2]^2/2,gr,h)
  expect_false(out$review$Verified)
  expect_null(out$cov)
})

test_that("the information budget is enforced before extra refinement allocations", {
  h <- diag(c(1e-10,2))
  local_mocked_bindings(
    build_param_sizes=function(...) list(x=2L),
    build_param_slices=function(...) list(x=1:2),
    build_indices=function(...) list(),
    make_param_cache=function(...) {
      cache <- new.env(parent=emptyenv())
      cache$ensure <- function(p) cache$par <- p
      cache
    },
    mfrm_loglik_mml_cached=function(cache,...) sum(cache$par*(h%*%cache$par))/2,
    mfrm_grad_mml_cached=function(cache,...) drop(h%*%cache$par),
    mfrm_refine_mml_information=function(...) stop('unexpected extra allocation'),
    .package='mfrmr')
  withr::local_options(mfrmr.max_information_bytes=8*8*4)
  for(model in c('RSM','PCM','GPCM')) {
    fit <- list(summary=data.frame(Method='MML'),opt=list(par=c(0,0)),prep=list(),
      config=list(model=model,estimation_control=list(quad_points=3)))
    out <- compute_mml_parameter_covariance(fit)
    expect_null(out$cov)
    expect_false(out$solution_information$inverse_review$Verified)
    expect_match(out$detail,'refinement workspace estimate')
    expect_match(out$detail,'mfrmr.max_information_bytes',fixed=TRUE)
  }
})
