test_that("the declared per-Person review never relaxes the original gradient rule", {
  config <- list(n_person = 60L, estimation_control = list(gradient_review_per_person = 1e-6))
  expect_identical(mfrm_optimizer_gradient_tolerance(config, 1e-10), 60 * 1e-6)
  for (n in c(20L, 60L, 100L, 240L, 960L)) for (tol in c(1e-6, 1e-10, 1e-14)) {
    config$n_person <- n
    expect_equal(mfrm_optimizer_gradient_tolerance(config, tol),
      min(max(1e-4, 10 * tol), n * 1e-6))
  }
  config$estimation_control$gradient_review_per_person <- NULL
  expect_null(mfrm_optimizer_gradient_tolerance(config, 1e-10))
  for (bad in list(0, NA_real_, Inf, TRUE, 1i, c(1, 2))) {
    config$estimation_control$gradient_review_per_person <- bad
    expect_error(mfrm_optimizer_gradient_tolerance(config, 1e-10), "per-Person")
  }
  config$estimation_control$gradient_review_per_person <- 1e-6
  config$n_person <- 1.5
  expect_error(mfrm_optimizer_gradient_tolerance(config, 1e-10), "per-Person")
})

test_that("a stricter declared rule starts existing polishing without changing the objective", {
  fn <- function(p) 1000 + sum(p^2)/2
  gr <- function(p) p
  calls <- 0L; worse <- FALSE
  local_mocked_bindings(make_param_cache = function(...) list(),
    make_mfrm_direct_evaluator = function(...) list(value = fn, gradient = gr,
      diagnostics = function() list()),
    optim = function(par, ...) {
      calls <<- calls + 1L
      p <- if (calls == 1L) par else if (worse) rep(1, length(par)) else rep(0, length(par))
      list(par = p, value = fn(p), convergence = 0L,
        counts = c("function" = 1L, "gradient" = 1L))
    }, .package = "mfrmr")
  run <- function(n = 60L, declared = TRUE) {
    calls <<- 0L
    config <- list(model = "GPCM", n_person = n,
      estimation_control = list(mml_integration = "adaptive",
        gradient_review_per_person = if (declared) 1e-6 else NULL))
    run_mfrm_direct_optimization(c(8e-5, 0), "MML", list(), config, list(),
      quad_points = 3L, maxit = 20L, reltol = 1e-10, optimizer = "BFGS")
  }
  finished <- run()
  expect_gt(calls, 1L)
  expect_equal(finished$par, c(0, 0))
  expect_lte(finished$value, fn(c(8e-5, 0)))
  expect_equal(finished$optimizer_diagnostics$GradientReviewTolerance, 6e-5)
  expect_identical(finished$optimizer_diagnostics$ConvergenceSeverity, "pass")
  expect_true(finished$optimizer_polish$Triggered)
  for (n in c(240L, 960L)) {
    unchanged <- run(n)
    expect_identical(calls, 1L)
    expect_identical(unchanged$par, c(8e-5, 0))
    expect_identical(unchanged$optimizer_diagnostics$GradientReviewTolerance, 1e-4)
  }
  legacy <- run(declared = FALSE)
  expect_identical(calls, 1L)
  expect_identical(legacy$par, c(8e-5, 0))
  worse <- TRUE
  expect_warning(rejected <- run(), "convergence review")
  expect_identical(rejected$par, c(8e-5, 0))
  expect_identical(rejected$optimizer_diagnostics$ConvergenceSeverity, "review")
  expect_false(rejected$optimizer_polish$Succeeded)
})

test_that("gradient diagnostics retain the declared tolerance and reject invalid values", {
  opt <- list(par = c(0, 0), convergence = 0L)
  old <- build_optimizer_diagnostics(opt, c(8e-5, 0), reltol = 1e-10)
  new <- build_optimizer_diagnostics(opt, c(8e-5, 0), reltol = 1e-10,
    gradient_tolerance = 6e-5)
  expect_identical(old$ConvergenceSeverity, "pass")
  expect_identical(new$ConvergenceSeverity, "review")
  expect_identical(new$ConvergenceReason, "code_zero_large_gradient")
  expect_identical(new$GradientReviewTolerance, 6e-5)
  expect_identical(new, unserialize(serialize(new, NULL)))
  for (bad in list(0, NA_real_, Inf, TRUE, 1i, c(1, 2)))
    expect_error(build_optimizer_diagnostics(opt, c(0, 0), gradient_tolerance = bad),
      "finite positive")
})
