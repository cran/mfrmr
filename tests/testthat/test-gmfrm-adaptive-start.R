test_that("adaptive starts retain better unfinished solutions and failed alternatives", {
  problem <- list(start = 0, common = list(idx = list(), sizes = list()))
  config <- list(estimation_control = list(mml_integration = "adaptive"))
  mode <- "both"; calls <- numeric()
  local_mocked_bindings(
    mfrm_gmfrm_em = function(problem, maxit, score_tol) {
      expect_equal(maxit, 5L); expect_equal(score_tol, 1e-6)
      if (mode == "em_error") stop("seed unavailable")
      list(par = 1, converged = FALSE, reason = "iteration_limit")
    },
    mfrmr_make_adaptive_mml_evaluator = function(...) function(par) list(value = 10 + par),
    run_mfrm_optimization = function(start, method, idx, config, sizes, quad_points,
                                     maxit, reltol, optimizer) {
      calls <<- c(calls, start)
      if (mode == "all_error" || (mode == "neutral_error" && start == 0)) stop("optimizer unavailable")
      value <- if (mode == "worse_start") 12 + start else if (mode == "tie") 5 else 5 - start
      list(par = start, value = value, convergence = if (start == 1) 1L else 0L,
        optimizer_diagnostics = list(ConvergenceSeverity = if (start == 1) "review" else "pass"),
        optimizer_polish = list(Stages = data.frame(Stage = 1L)))
    }, .package = "mfrmr")
  run <- function(policy = "neutral_em") mfrm_gmfrm_adaptive_fit(problem, config,
    31L, 5L, 1e-10, "BFGS", policy)
  fit <- run()
  expect_identical(calls, c(0, 1))
  expect_identical(fit$par, 1)
  expect_identical(fit$convergence, 1L)
  expect_false(fit$mml_initialization$em$value$converged)
  expect_identical(fit$mml_initialization$selected, "em")
  expect_identical(fit$mml_initialization$table$Selected, c(FALSE, TRUE))
  expect_equal(fit$mml_initialization$table$InitialNLL, c(10, 11))
  expect_equal(fit$mml_initialization$table$FinalNLL, c(5, 4))
  mode <- "tie"
  expect_identical(run()$mml_initialization$selected, "neutral")
  mode <- "neutral_error"
  expect_warning(fit <- run(), "Some adaptive starting values failed")
  expect_match(fit$mml_initialization$table$Error[1], "optimizer unavailable")
  expect_identical(fit$par, 1)
  mode <- "em_error"
  expect_warning(fit <- run(), "Some adaptive starting values failed")
  expect_identical(fit$par, 0)
  expect_match(fit$mml_initialization$em$error, "seed unavailable")
  mode <- "all_error"
  error <- tryCatch(run(), error = identity)
  expect_s3_class(error, "mfrmr_gmfrm_start_error")
  expect_equal(nrow(error$initialization$table), 2L)
  expect_identical(error$initialization$controls$quad_points, 31L)
  mode <- "worse_start"
  expect_warning(fit <- run(), "starting value has a lower")
  expect_identical(fit$convergence, 2L)
  expect_identical(fit$optimizer_diagnostics$ConvergenceReason, "better_starting_objective")
  mode <- "em_error"; calls <- numeric()
  fit <- run("neutral")
  expect_identical(calls, 0)
  expect_null(fit$mml_initialization$em)
})

test_that("explicit native starts use the adaptive objective and retain every attempt", {
  problem <- list(start = c(native = 0), common = list(idx = list(), sizes = list()))
  config <- list(estimation_control = list(mml_integration = "adaptive"))
  starts <- list(external = c(native = 1), neutral = c(native = 0))
  seen <- list(); failed <- FALSE
  local_mocked_bindings(mfrm_gmfrm_em = function(...) stop("must not seed from EM"),
    mfrmr_make_adaptive_mml_evaluator = function(...) function(p) list(value = 10 + unname(p)),
    run_mfrm_optimization = function(start, ...) {
      seen[[length(seen) + 1L]] <<- start
      if (failed) stop("explicit-start failure")
      list(par = start, value = 5 - unname(start), convergence = if (start == 1) 1L else 0L,
        optimizer_diagnostics = list(ConvergenceSeverity = if (start == 1) "review" else "pass"))
    }, .package = "mfrmr")
  run <- function() mfrm_gmfrm_adaptive_fit(problem, config, 31L, 5L, 1e-9, "BFGS", starts)
  fit <- run()
  expect_identical(seen, unname(starts))
  expect_identical(fit$mml_initialization$policy, "explicit")
  expect_null(fit$mml_initialization$em)
  expect_identical(lapply(fit$mml_initialization$attempts, `[[`, "start"), starts)
  expect_equal(fit$mml_initialization$table$InitialNLL, c(11, 10))
  expect_identical(fit$mml_initialization$selected, "external")
  expect_identical(fit$convergence, 1L)
  failed <- TRUE
  error <- tryCatch(run(), error = identity)
  expect_s3_class(error, "mfrmr_gmfrm_start_error")
  expect_identical(error$initialization$table$Start, names(starts))
  expect_true(all(grepl("explicit-start failure", error$initialization$table$Error)))
})

test_that("explicit native starts reject malformed coordinates before numerical work", {
  problem <- list(start = c(a = 0, b = 0), common = list(idx = list(), sizes = list()))
  local_mocked_bindings(mfrm_gmfrm_em = function(...) stop("unexpected EM"),
    mfrmr_make_adaptive_mml_evaluator = function(...) stop("numerical stage reached"), .package = "mfrmr")
  run <- function(starts) mfrm_gmfrm_adaptive_fit(problem, list(), 31L, 5L, 1e-9, "BFGS", starts)
  for (bad in list(list(), list(only = c(0, 0)), list(c(0, 0), c(1, 1)),
      setNames(list(c(0, 0), c(1, 1)), c("same", "same")),
      setNames(list(c(0, 0), c(1, 1)), c("first", " ")),
      data.frame(a = c(0, 0), b = c(1, 1))))
    expect_error(run(bad), "plain named list")
  for (bad in list(0, c(0, NA), c(0, Inf), c(0, 1i), c("0", "1"),
      matrix(0, 1, 2), c(b = 0, a = 0), c(a = 0, wrong = 0)))
    expect_error(run(list(first = bad, second = c(0, 0))), "free-parameter length and order")
  expect_error(run(list(first = c(a = 0, b = 0), second = c(1, 1))), "numerical stage reached")
})

test_that("quadrature replay preserves explicit vectors and refuses altered start records", {
  fit <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))$fit
  starts <- list(retained = fit$opt$par, neutral = 0 * fit$opt$par)
  fit$config$estimation_control$gpcm_mml_start <- "explicit"
  fit$config$replay_inputs$gpcm_mml_start <- starts
  fit$opt$mml_initialization <- list(policy = "explicit",
    attempts = lapply(starts, function(x) list(start = x)), table = data.frame(Start = names(starts)))
  args <- mfrmr_gqs_refit_arguments(fit, fit$gmfrm$specification$data, 61L)
  expect_identical(args$gpcm_mml_start, starts)
  expect_identical(args$quad_points, 61L)
  bad <- fit; bad$config$replay_inputs$gpcm_mml_start[[1]][1] <- 999
  expect_error(mfrmr_gqs_refit_arguments(bad, args$data, 61L), "vectors and replay")
  bad <- fit; bad$opt$mml_initialization$attempts <- rev(bad$opt$mml_initialization$attempts)
  expect_error(mfrmr_gqs_refit_arguments(bad, args$data, 61L), "vectors and replay")
  bad <- fit; bad$config$replay_inputs$gpcm_mml_start <- NULL
  expect_error(mfrmr_gqs_refit_arguments(bad, args$data, 61L), "vectors and replay")
})

test_that("explicit public fitting replays real optimizer outcomes without importing source qualification", {
  x <- readRDS(test_path("fixtures", "gmfrm-explicit-start.rds"))
  seen <- character()
  local_mocked_bindings(mfrm_gmfrm_em = function(...) stop("no EM seed"),
    run_mfrm_optimization = function(start, method, idx, config, sizes, quad_points, maxit, reltol, optimizer) {
      expect_identical(method, "MML")
      expect_identical(config$estimation_control$mml_integration, "adaptive")
      at <- which(vapply(x$args$gpcm_mml_start, identical, TRUE, y = start))
      expect_length(at, 1L)
      name <- names(x$args$gpcm_mml_start)[at]
      seen <<- c(seen, name)
      x$fit$opt$mml_initialization$attempts[[name]]$value
    }, .package = "mfrmr")
  fit <- do.call(fit_mfrm, x$args)
  expect_identical(seen, names(x$args$gpcm_mml_start))
  expect_s3_class(fit, "mfrm_fit")
  expect_identical(fit$opt$par, x$fit$opt$par)
  expect_identical(fit$summary$Converged, TRUE)
  expect_identical(fit$summary$ConvergenceBasis, "optimizer_gradient")
  expect_identical(fit$summary$MMLIntegration, "adaptive")
  expect_identical(fit$summary$MMLEngineUsed, "direct")
  expect_equal(fit$opt$mml_initialization$table$InitialNLL,
    x$fit$opt$mml_initialization$table$InitialNLL, tolerance = 1e-12)
  expect_identical(fit$config$replay_inputs$gpcm_mml_start, x$args$gpcm_mml_start)
  expect_true(all(summary(fit)$settings_overview$GpcmMmlInitialization == "explicit"))
  expect_false(fit$summary$InferenceReady)
  expect_null(fit$opt$mml_initialization$em)
  args <- mfrmr_gqs_refit_arguments(fit, x$args$data, 61L)
  expect_identical(args$gpcm_mml_start, x$args$gpcm_mml_start)
  # The prior segmented result's local success has not become a public fit flag.
  source <- readRDS(test_path("fixtures", "gmfrm-result.rds"))
  expect_identical(source$run$admission$public_fit, FALSE)
  malformed <- x$args; names(malformed$gpcm_mml_start[[1]]) <- rep("wrong", length(fit$opt$par))
  expect_error(do.call(fit_mfrm, malformed), "names must match exactly")
  malformed$mml_integration <- "fixed"; malformed$mml_engine <- "em"
  expect_error(do.call(fit_mfrm, malformed), "only to two-family adaptive")
})

test_that("explicit starts persist through saved reports and retain ordinary scoring contracts", {
  x <- readRDS(test_path("fixtures", "gmfrm-explicit-start.rds"))
  forbidden <- function(...) stop("Saved output must not refit or recompute information")
  local_mocked_bindings(fit_mfrm = forbidden, run_mfrm_optimization = forbidden,
    compute_mml_parameter_covariance = forbidden, mfrm_gpcm_product_inference = forbidden,
    .package = "mfrmr")
  res <- mfrm_results(x$fit, include = "fit", compute = "never", intervals = x$intervals, scores = x$scores)
  report <- mfrm_report(res)
  expect_identical(report$tables$gpcm_initialization, x$fit$opt$mml_initialization$table)
  expect_identical(report$tables$model_settings$Value[report$tables$model_settings$Setting == "Initialization"], "explicit")
  expect_true(all(res$tables$gpcm_slopes_intervals$CIEligible))
  expect_false(any(res$tables$gpcm_slopes_intervals$CoverageQualified))
  expect_true(all(res$tables$gpcm_locations_intervals$CIEligible))
  expect_identical(x$scores$settings$source_scoring_policy_basis, "two_family_conditional_eap_v1")
  file <- tempfile(); on.exit(unlink(file), add = TRUE)
  saveRDS(res, file); restored <- readRDS(file)
  expect_identical(restored, res)
  expect_identical(mfrm_report(restored)$markdown, report$markdown)
  expect_identical(mfrmr_gqs_refit_arguments(restored$fit, x$args$data, 61L)$gpcm_mml_start,
    x$args$gpcm_mml_start)
})
