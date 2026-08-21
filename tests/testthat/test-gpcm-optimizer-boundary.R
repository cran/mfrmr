test_that("GPCM slope numeric-boundary conditions retain audit payloads", {
  condition <- mfrmr:::new_gpcm_slope_numeric_boundary_error(
    c(800, -800, 0)
  )

  expect_s3_class(condition, "mfrmr_gpcm_slope_numeric_boundary_error")
  expect_match(
    conditionMessage(condition),
    "finite and strictly positive"
  )
  expect_equal(condition$expanded_log_slopes, c(800, -800, 0))
  expect_true(any(!is.finite(condition$expanded_slopes)))
  expect_true(any(condition$expanded_slopes <= 0))
})

test_that("direct objectives reject only typed non-representable slope trials", {
  evaluator <- list(value = function(par) {
    if (abs(par[1]) > 0.1) {
      stop(mfrmr:::new_gpcm_slope_numeric_boundary_error(par[1]))
    }
    (par[1] - 1)^2
  })
  safe <- mfrmr:::make_mfrm_boundary_safe_objective(evaluator)

  expect_equal(safe$value(0), 1)
  expect_equal(safe$value(2), 1e100)
  expect_identical(safe$rejections(), 1L)
  expect_true(is.finite(safe$value(2)))
  expect_identical(safe$rejections(), 2L)

  unrelated <- mfrmr:::make_mfrm_boundary_safe_objective(list(
    value = function(par) stop("unrelated evaluator failure")
  ))
  expect_error(unrelated$value(0), "unrelated evaluator failure")
  expect_identical(unrelated$rejections(), 0L)
})

test_that("boundary-safe BFGS contracts an invalid trial instead of aborting", {
  evaluator <- list(value = function(par) {
    if (abs(par[1]) > 0.1) {
      stop(mfrmr:::new_gpcm_slope_numeric_boundary_error(par[1]))
    }
    (par[1] - 0.08)^2
  })
  safe <- mfrmr:::make_mfrm_boundary_safe_objective(evaluator)

  opt <- stats::optim(
    par = 0,
    fn = safe$value,
    gr = function(par) 2 * (par - 0.08),
    method = "BFGS",
    control = list(maxit = 50L, reltol = 1e-12)
  )

  expect_false(inherits(opt, "error"))
  expect_true(is.finite(opt$value))
  expect_lte(abs(opt$par[1]), 0.1)
  expect_gt(safe$rejections(), 0L)
})
