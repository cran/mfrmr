test_that('testlet start selection prefers convergence only for numerical ties', {
  runs <- function(values, codes) Map(function(v, code)
    list(value=v, convergence=code), values, codes)
  eps <- .Machine$double.eps
  expect_identical(mfrm_testlet_select_run(runs(c(1,1+10*eps),c(52L,0L))),2L)
  expect_identical(mfrm_testlet_select_run(runs(c(1,1+1e-6),c(52L,0L))),1L)
  expect_identical(mfrm_testlet_select_run(runs(c(1,1),c(0L,52L))),1L)
  expect_identical(mfrm_testlet_select_run(runs(c(2,1),c(52L,52L))),2L)
  expect_identical(mfrm_testlet_select_run(runs(c(1000,1000+1e-12),c(52L,0L))),2L)
  expect_identical(mfrm_testlet_select_run(runs(c(1000,1000+1e-6),c(52L,0L))),1L)
  expect_identical(mfrm_testlet_select_run(runs(c(1,1+eps,1+2*eps),c(52L,0L,0L))),2L)
  # Failed/nonfinite alternatives must not displace an already converged minimum.
  expect_identical(mfrm_testlet_select_run(runs(c(Inf,1,NA_real_),c(52L,0L,0L))),2L)
})
