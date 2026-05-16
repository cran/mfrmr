# Tests for the package-level MnSq misfit threshold pair
# (`mfrm_misfit_thresholds()`) and its propagation to summary,
# build_misfit_casebook, build_apa_outputs, and facet_quality_dashboard.

test_that("mfrm_misfit_thresholds() defaults to Linacre 0.5-1.5", {
  thr <- mfrm_misfit_thresholds()
  expect_named(thr, c("lower", "upper"))
  expect_equal(unname(thr), c(0.5, 1.5))
})

test_that("mfrm_misfit_thresholds() respects R options", {
  old <- options(
    mfrmr.misfit_lower = 0.7,
    mfrmr.misfit_upper = 1.3
  )
  on.exit(options(old), add = TRUE)
  thr <- mfrm_misfit_thresholds()
  expect_equal(unname(thr), c(0.7, 1.3))
})

test_that("mfrm_misfit_thresholds() respects per-call overrides", {
  thr <- mfrm_misfit_thresholds(lower = 0.6, upper = 1.4)
  expect_equal(unname(thr), c(0.6, 1.4))
})

test_that("mfrm_misfit_thresholds() rejects invalid bounds", {
  expect_error(mfrm_misfit_thresholds(lower = 1.5, upper = 0.5),
               "0 < lower < upper")
  expect_error(mfrm_misfit_thresholds(lower = -0.1, upper = 1.5),
               "0 < lower < upper")
})

test_that("summary(diag) inherits the option-driven band", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 25)
  ))
  diag <- suppressMessages(diagnose_mfrm(fit, residual_pca = "none",
                                          diagnostic_mode = "legacy"))
  old <- options(
    mfrmr.misfit_lower = 0.7,
    mfrmr.misfit_upper = 1.3
  )
  on.exit(options(old), add = TRUE)
  s <- summary(diag)
  expect_equal(unname(s$misfit_thresholds), c(0.7, 1.3))
})
