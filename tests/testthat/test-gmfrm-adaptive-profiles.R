test_that("adaptive component profiles retain the moving-node objective and constraints", {
  x <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))
  fit <- x$fit; p <- fit$opt$par
  before <- fit
  checked <- FALSE
  local_mocked_bindings(
    mfrm_gmfrm_profile_evaluator = function(...) stop("must not use fixed-grid EM"),
    mfrm_gpcm_profile_search = function(evaluator, reference, start, contrast, ...) {
      checked <<- TRUE
      expect_equal(evaluator$optimization_scale, length(fit$prep$levels$Person))
      expect_equal(reference$optimization_scale, evaluator$optimization_scale)
      expect_equal(evaluator$value(p), x$baseline$continuous_nll, tolerance = 1e-7)
      expect_lte(max(abs(evaluator$gradient(p) + colSums(x$baseline$scores))), 1e-4)
      expected <- numeric(length(p)); expected[15:16] <- -1
      expect_identical(unname(contrast), expected)
      # At a changed vector, finite differences include moving modes/scales;
      # no frozen posterior or fixed-grid evaluator supplies these gradients.
      trial <- p + seq(-.01, .015, length.out = length(p))
      fd <- vapply(seq_along(trial), function(j) {
        lo <- hi <- trial; lo[j] <- lo[j] - 1e-5; hi[j] <- hi[j] + 1e-5
        (evaluator$value(hi) - evaluator$value(lo)) / 2e-5
      }, 0)
      expect_lte(max(abs(fd - evaluator$gradient(trial))), 1e-5)
      expect_lte(abs(reference$value(p) - evaluator$value(p)), 1e-5)
      attr(x$profiles$Task, "profile")
    }, .package = "mfrmr")
  warnings <- character()
  ci <- withCallingHandlers(confint(fit, parm = "slopes", method = "profile", slope = c(Task = "t3")),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  expect_length(warnings, 2L)
  expect_match(warnings[1], "Experimental two-family|rank is deficient")
  expect_match(warnings[2], "Profile interval unresolved")
  expect_true(checked)
  expect_identical(fit, before)
  expect_identical(attr(ci, "settings")$integration, "adaptive")
  expect_identical(attr(ci, "settings")$engine, "direct")
  expect_identical(attr(ci, "profile")$endpoints, attr(x$profiles$Task, "profile")$endpoints)
  expect_false(attr(ci, "diagnostics")$CIEligible)
})

test_that("adaptive profile reports and replay preserve unavailable endpoints and integration", {
  x <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))
  local_mocked_bindings(mfrm_gpcm_profile_search = function(...) stop("must not reprofile"),
    mfrm_gpcm_product_inference = function(...) stop("must not recompute information"),
    fit_mfrm = function(...) stop("must not refit"), .package = "mfrmr")
  ci <- x$profiles$Task
  pr <- attr(ci, "profile")
  expect_true(any(!pr$profile$Passed))
  expect_true(all(pr$endpoints$Status == "numerical_or_start_failure"))
  expect_true(all(is.na(pr$endpoints$Bound)))
  figure <- plot(ci, draw = FALSE)
  expect_identical(plot_data(figure, "endpoints"), pr$endpoints)
  expect_true(any(plot_data(figure, "table")$Display == "Unresolved"))
  expect_no_error(ggplot2::ggplot_build(figure))
  res <- mfrm_results(x$fit, include = c("fit", "plots"), compute = "never",
    intervals = list(component = ci))
  report <- mfrm_report(res)
  expect_identical(report$tables$gpcm_component_settings$Integration, "adaptive")
  expect_identical(report$tables$gpcm_component_settings$Engine, "direct")
  expect_identical(report$tables$gpcm_component_profile_endpoints, pr$endpoints)
  path <- tempfile(); on.exit(unlink(path, recursive = TRUE), add = TRUE)
  exported <- export_mfrm_results(res, output_dir = path, preset = "starter",
    acknowledge_sensitive = TRUE)
  expect_equal(nrow(exported$plot_errors), 0L)
  saved <- readRDS(exported$written_files$Path[exported$written_files$Component == "results_rds"])
  expect_identical(saved$gpcm_inference, res$gpcm_inference)
  expect_identical(mfrm_report(saved)$tables$gpcm_component_settings,
    report$tables$gpcm_component_settings)
})
