# Tests for empirical-Bayes facet shrinkage (added in 0.1.6).
# Cross-refs:
#  * `.compute_facet_shrinkage()` — core math (internal)
#  * `apply_empirical_bayes_shrinkage()` — standalone post-hoc wrapper
#  * `shrinkage_report()` — accessor
#  * `fit_mfrm(..., facet_shrinkage = ...)` — integrated path
#  * `build_mfrm_manifest()$shrinkage_review` — reproducibility trail
#  * `reporting_checklist()` — "Empirical-Bayes shrinkage" item

# --- Closed-form math --------------------------------------------------------

test_that(".compute_facet_shrinkage matches hand calculation", {
  fn <- getFromNamespace(".compute_facet_shrinkage", "mfrmr")
  # Known example: K = 4 estimates ~ N(0, 1), SEs all 0.1.
  # Expected: tau^2 = var(est) - mean(se^2); B = se^2/(tau^2 + se^2).
  est <- c(-1.0, -0.5, 0.5, 1.0)
  se  <- rep(0.1, 4)
  out <- fn(est, se, method = "empirical_bayes", min_levels = 3L)
  expect_equal(out$n_levels, 4L)
  expect_equal(out$n_levels_used, 4L)
  expected_tau2 <- mean(est^2) - mean(se^2)
  expect_equal(out$tau2, max(0, expected_tau2), tolerance = 1e-12)
  expected_B <- se^2 / (out$tau2 + se^2)
  expect_equal(out$shrinkage_factors, expected_B, tolerance = 1e-12)
  expect_equal(out$shrunk_estimates, (1 - expected_B) * est, tolerance = 1e-12)
  expect_equal(out$shrunk_ses, sqrt((1 - expected_B) * se^2), tolerance = 1e-12)
})

test_that(".compute_facet_shrinkage respects user prior_sd", {
  fn <- getFromNamespace(".compute_facet_shrinkage", "mfrmr")
  est <- c(-0.4, 0.1, 0.3, 0.0)
  se  <- rep(0.2, 4)
  out <- fn(est, se, method = "empirical_bayes", prior_sd = 1.5)
  expect_equal(out$tau2, 1.5^2)
  expect_identical(out$prior_sd_source, "user")
})

test_that(".compute_facet_shrinkage collapses when tau2 <= 0", {
  fn <- getFromNamespace(".compute_facet_shrinkage", "mfrmr")
  # Huge SE relative to spread -> MoM tau2 goes negative, clamped to 0.
  est <- c(-0.02, 0.01, 0.005, -0.015)
  se  <- rep(1.0, 4)
  out <- fn(est, se, method = "empirical_bayes")
  expect_equal(out$tau2, 0)
  expect_true(all(out$shrinkage_factors == 1))
  expect_true(all(out$shrunk_estimates == 0))
  expect_true(!is.na(out$note))
})

test_that(".compute_facet_shrinkage passes through when K < min_levels", {
  fn <- getFromNamespace(".compute_facet_shrinkage", "mfrmr")
  est <- c(-0.1, 0.2)
  se  <- c(0.1, 0.1)
  out <- fn(est, se, min_levels = 3L)
  expect_equal(out$shrunk_estimates, est)
  expect_equal(out$shrunk_ses, se)
  expect_equal(out$shrinkage_factors, c(0, 0))
  expect_true(!is.na(out$note))
})

# --- fit_mfrm integration ----------------------------------------------------

test_that("fit_mfrm(facet_shrinkage = 'none') is a no-op", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15)
  ))
  expect_identical(as.character(fit$config$facet_shrinkage), "none")
  expect_false("ShrunkEstimate" %in% names(fit$facets$others))
  expect_null(fit$shrinkage_report)
})

test_that("fit_mfrm(facet_shrinkage = 'empirical_bayes') populates schema", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15,
             facet_shrinkage = "empirical_bayes")
  ))
  expect_identical(as.character(fit$config$facet_shrinkage), "empirical_bayes")
  expect_true(all(c("ShrunkEstimate", "ShrunkSE", "ShrinkageFactor") %in%
                    names(fit$facets$others)))
  expect_s3_class(fit$shrinkage_report, "data.frame")
  expect_true(all(c("Facet", "NLevels", "Tau2", "MeanShrinkage",
                    "EffectiveDF") %in% names(fit$shrinkage_report)))
  # Shrinkage factors are in [0, 1].
  bf <- suppressWarnings(as.numeric(fit$facets$others$ShrinkageFactor))
  expect_true(all(bf >= 0 - 1e-10 & bf <= 1 + 1e-10, na.rm = TRUE))
})

test_that("apply_empirical_bayes_shrinkage gives identical output to integrated path", {
  toy <- load_mfrmr_data("example_core")
  fit_a <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15,
             facet_shrinkage = "empirical_bayes")
  ))
  fit_b0 <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15)
  ))
  fit_b <- apply_empirical_bayes_shrinkage(fit_b0)
  expect_equal(fit_a$shrinkage_report$Tau2, fit_b$shrinkage_report$Tau2)
  expect_equal(fit_a$shrinkage_report$MeanShrinkage,
               fit_b$shrinkage_report$MeanShrinkage)
})

test_that("shrinkage_report() accessor gives NULL when no shrinkage was applied", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15)
  ))
  expect_message(out <- shrinkage_report(fit), "No shrinkage applied")
  expect_null(out)
})

test_that("shrinkage_report() returns a data.frame after shrinkage", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15,
             facet_shrinkage = "empirical_bayes")
  ))
  r <- shrinkage_report(fit)
  expect_s3_class(r, "data.frame")
  expect_setequal(r$Facet, c("Rater", "Criterion"))
})

# --- Reporting integration ---------------------------------------------------

test_that("build_mfrm_manifest records shrinkage_review", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15,
             facet_shrinkage = "empirical_bayes")
  ))
  m <- build_mfrm_manifest(fit)
  expect_true("shrinkage_review" %in% names(m))
  expect_true(nrow(m$shrinkage_review) >= 1L)
  expect_true("Mode" %in% names(m$shrinkage_review))
  expect_identical(unique(m$shrinkage_review$Mode), "empirical_bayes")
})

test_that("reporting_checklist has a shrinkage item", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15,
             facet_shrinkage = "empirical_bayes")
  ))
  diag <- suppressMessages(suppressWarnings(
    diagnose_mfrm(fit, residual_pca = "none")
  ))
  chk <- reporting_checklist(fit, diagnostics = diag)
  hit <- grepl("shrinkage", chk$checklist$Item, ignore.case = TRUE)
  expect_true(any(hit))
})

# --- Polish: PosteriorSD + SE columns ---------------------------------------

test_that("MML fit exposes PosteriorSD and SE alongside legacy SD", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "MML", quad_points = 7)
  ))
  expect_true(all(c("SD", "PosteriorSD", "SE") %in% names(fit$facets$person)))
  # Aliases should carry the same numeric values.
  expect_equal(fit$facets$person$SD, fit$facets$person$PosteriorSD)
  expect_equal(fit$facets$person$SE, fit$facets$person$PosteriorSD)
})

test_that("JML fit exposes NA SE column for consistency with MML", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15)
  ))
  expect_true("SE" %in% names(fit$facets$person))
  expect_true(all(is.na(fit$facets$person$SE)))
})

# --- Shrinkage visualisation --------------------------------------------------

test_that("plot.mfrm_fit(type = 'shrinkage') returns a plot bundle", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15,
             facet_shrinkage = "empirical_bayes")
  ))
  out <- plot(fit, type = "shrinkage", draw = FALSE)
  expect_s3_class(out, "mfrm_plot_data")
  expect_true(is.data.frame(out$data$data))
  expect_true(all(c("Facet", "Level", "Estimate", "ShrunkEstimate",
                    "ShrinkageFactor") %in% names(out$data$data)))
  expect_identical(out$data$mode, "empirical_bayes")
})

test_that("plot.mfrm_fit(type = 'shrinkage') gracefully renders with no shrinkage", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15)
  ))
  # Without shrinkage, draw = FALSE should still return a structure
  # with an empty data frame and mode == "none".
  out <- plot(fit, type = "shrinkage", draw = FALSE)
  expect_s3_class(out, "mfrm_plot_data")
  expect_equal(nrow(out$data$data), 0L)
  expect_identical(out$data$mode, "none")
})

test_that("plot.mfrm_facet_sample_review runs and returns a data.frame", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15)
  ))
  audit <- facet_small_sample_review(fit)
  pdf(NULL)  # suppress device
  on.exit(dev.off(), add = TRUE)
  out <- plot(audit)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("Facet", "Level", "N", "SampleCategory") %in% names(out)))
})

test_that("plot.mfrm_facet_nesting runs and returns a matrix", {
  toy <- load_mfrmr_data("example_core")
  nest <- detect_facet_nesting(toy, c("Rater", "Criterion"))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  m <- plot(nest)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 2L))
})
