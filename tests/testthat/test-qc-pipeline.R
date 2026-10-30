make_qc_study1_fit <- function() {
  toy <- load_mfrmr_data("study1")
  expect_warning(
    fit <- .mfrmr_muffle_expected_warnings(
      fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
               method = "JML", maxit = 25),
      "^Optimization convergence review did not produce"),
    "^Category support is retained but requires review:")
  expect_identical(fit$summary$CategoryState, "weak_information")
  fit
}

test_that("run_qc_pipeline returns correct class and structure", {
  fit <- make_qc_study1_fit()
  qc <- run_qc_pipeline(fit)

  expect_s3_class(qc, "mfrm_qc_pipeline")
  expect_true(is.list(qc))
  expect_named(qc, c("verdicts", "overall", "details", "recommendations", "config", "gpcm_boundary"),
               ignore.order = TRUE)
})

test_that("verdicts table has 10 rows with valid verdict values", {
  fit <- make_qc_study1_fit()
  qc <- run_qc_pipeline(fit)


  expect_equal(nrow(qc$verdicts), 10)
  expect_true(all(qc$verdicts$Verdict %in% c("Pass", "Warn", "Fail", "Skip")))
  expected_checks <- c("Convergence", "Global Fit", "Reliability", "Separation",
                        "Element Misfit", "Unexpected Responses", "Category Structure",
                        "Connectivity", "Inter-rater Agreement", "Functioning/Bias Screen")
  expect_equal(qc$verdicts$Check, expected_checks)
  expect_true(all(c("Check", "Verdict", "Value", "Threshold", "Detail") %in%
                    names(qc$verdicts)))
})

test_that("overall verdict is one of Pass/Warn/Fail", {
  fit <- make_qc_study1_fit()
  qc <- run_qc_pipeline(fit)

  expect_true(qc$overall %in% c("Pass", "Warn", "Fail"))
})

test_that("recommendations is a character vector", {
  fit <- make_qc_study1_fit()
  qc <- run_qc_pipeline(fit)

  expect_type(qc$recommendations, "character")
})

test_that("print.mfrm_qc_pipeline produces output", {
  fit <- make_qc_study1_fit()
  qc <- run_qc_pipeline(fit)

  out <- capture.output(print(qc))
  expect_true(length(out) > 0)
  expect_true(any(grepl("QC Pipeline", out)))
  expect_true(any(grepl("Overall:", out)))
})

test_that("summary.mfrm_qc_pipeline produces correct output", {
  fit <- make_qc_study1_fit()
  qc <- run_qc_pipeline(fit)

  s <- summary(qc)
  expect_s3_class(s, "summary.mfrm_qc_pipeline")
  expect_true("pass_count" %in% names(s))
  expect_true("warn_count" %in% names(s))
  expect_true("fail_count" %in% names(s))
  expect_true("skip_count" %in% names(s))
  expect_equal(s$pass_count + s$warn_count + s$fail_count + s$skip_count, 10)

  out <- capture.output(print(s))
  expect_true(length(out) > 0)
  expect_true(any(grepl("QC Pipeline Summary", out)))
})

test_that("plot_qc_pipeline with draw=FALSE returns data", {
  fit <- make_qc_study1_fit()
  qc <- run_qc_pipeline(fit)

  vt <- plot_qc_pipeline(qc, draw = FALSE)
  expect_true(is.data.frame(vt) || tibble::is_tibble(vt))
  expect_equal(nrow(vt), 10)
})

test_that("threshold_profile strict and lenient work", {
  fit <- make_qc_study1_fit()

  qc_strict  <- run_qc_pipeline(fit, threshold_profile = "strict")
  qc_lenient <- run_qc_pipeline(fit, threshold_profile = "lenient")

  expect_s3_class(qc_strict, "mfrm_qc_pipeline")
  expect_s3_class(qc_lenient, "mfrm_qc_pipeline")
  expect_equal(qc_strict$config$threshold_profile, "strict")
  expect_equal(qc_lenient$config$threshold_profile, "lenient")

  # strict thresholds should be tighter
  expect_true(qc_strict$config$thresholds$reliability_pass >
                qc_lenient$config$thresholds$reliability_pass)
  expect_true(qc_strict$config$thresholds$separation_pass >
                qc_lenient$config$thresholds$separation_pass)
})

test_that("thresholds override works", {
  fit <- make_qc_study1_fit()

  qc <- run_qc_pipeline(fit, thresholds = list(reliability_pass = 0.99))
  expect_equal(qc$config$thresholds$reliability_pass, 0.99)
})

test_that("run_qc_pipeline works with pre-computed diagnostics", {
  fit <- make_qc_study1_fit()
  diag <- diagnose_mfrm(fit)
  qc <- run_qc_pipeline(fit, diagnostics = diag)

  expect_s3_class(qc, "mfrm_qc_pipeline")
  expect_equal(nrow(qc$verdicts), 10)
})

test_that("run_qc_pipeline records screening-tier bias metadata when bias results are provided", {
  toy <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    toy, "Person", c("Rater", "Criterion"), "Score",
    method = "JML", maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  bias <- suppressWarnings(estimate_bias(
    fit, diag,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))

  qc <- run_qc_pipeline(fit, diagnostics = diag, bias_results = bias)

  expect_match(qc$verdicts$Detail[10], "screened interactions crossed", fixed = TRUE)
  expect_identical(qc$details$bias$inference_tier, "screening")
  expect_true(qc$details$bias$available)
  expect_true(is.finite(qc$details$bias$total))
})

test_that("QC rejects invalid observation provenance before category assessment", {
  fit <- make_qc_study1_fit()
  diag <- diagnose_mfrm(fit, residual_pca = "none")
  diag$obs <- data.frame()
  expect_error(run_qc_pipeline(fit, diagnostics = diag), "fit/diagnostics mismatch", fixed = TRUE)
})

test_that("run_qc_pipeline surfaces incomplete bias collections as warn rather than skip/pass", {
  toy <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    toy, "Person", c("Rater", "Criterion"), "Score",
    method = "JML", maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  diag$interactions <- data.frame()
  bias_collection <- with_mocked_bindings(
    estimate_all_bias(fit, diagnostics = diag, pairs = list(c("Rater", "Criterion"))),
    estimate_bias = function(...) stop("forced pair failure"), .package = "mfrmr")

  qc <- run_qc_pipeline(fit, diagnostics = diag, bias_results = bias_collection)

  expect_identical(as.character(qc$verdicts$Verdict[10]), "Warn")
  expect_match(qc$verdicts$Detail[10], "incomplete", fixed = TRUE)
  expect_identical(qc$details$bias$error_count, 1L)
  expect_false(identical(qc$overall, "Pass"))
})

test_that("run_qc_pipeline rejects non-mfrm_fit input", {
  expect_error(run_qc_pipeline(list()), "mfrm_fit")
})

test_that("plot_qc_pipeline rejects non-qc-pipeline input", {
  expect_error(plot_qc_pipeline(list()), "mfrm_qc_pipeline")
})


test_that("QC cannot pass by ignoring unavailable facet reliability", {
  fit <- make_qc_study1_fit()
  diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  old <- diagnostics
  old$reliability$ModelSummaryNote <- NULL
  expect_error(mfrmr:::validate_diagnostics_precision(old),
               "levels used for reliability", fixed = TRUE)
  rel <- diagnostics$reliability
  idx <- which(rel$Facet != "Person")
  expect_gte(length(idx), 2L)
  rel$Reliability[idx] <- .99
  rel$Separation[idx] <- 10
  rel$Reliability[idx[1]] <- NA_real_
  rel$Separation[idx[1]] <- NA_real_
  diagnostics$reliability <- rel
  qc <- run_qc_pipeline(fit, diagnostics, include_bias = FALSE,
                        separation_facets = diagnostics$reliability$Facet[idx])
  rows <- qc$verdicts[qc$verdicts$Check %in% c("Reliability", "Separation"), ]
  expect_true(all(rows$Verdict == "Warn"))
  expect_true(all(rows$Value == "NA"))
  rel$Reliability <- rel$RealReliability <- NA_real_
  diagnostics$measures$ModelSE <- diagnostics$measures$RealSE <- NA_real_
  checks <- mfrmr:::audit_precision_outputs(
    fit, diagnostics$measures, rel, diagnostics$facet_precision, diagnostics$precision_profile
  )
  idx <- checks$Check %in% c("Fit-adjusted SE ordering", "Reliability ordering")
  expect_true(all(checks$Status[idx] == "review"))
  expect_true(all(grepl("No finite", checks$Detail[idx], fixed = TRUE)))
})

test_that("QC differentiation must be requested for named facets", {
  fit <- make_qc_study1_fit()
  diag <- diagnose_mfrm(fit, residual_pca = "none")
  diag$reliability$Reliability <- ifelse(diag$reliability$Facet == "Rater", 0, .99)
  diag$reliability$Separation <- ifelse(diag$reliability$Facet == "Rater", 0, 10)
  qc <- run_qc_pipeline(fit, diag, include_bias = FALSE)
  expect_identical(qc$verdicts$Verdict[3:4], c("Skip", "Skip"))
  expect_false(any(qc$verdicts$AffectsOverall[c(3, 4, 10)]))
  criterion <- run_qc_pipeline(fit, diag, separation_facets = "Criterion")
  expect_true(all(criterion$verdicts$Verdict[3:4] == "Pass"))
  rater <- run_qc_pipeline(fit, diag, separation_facets = "Rater")
  expect_true(all(rater$verdicts$Verdict[3:4] == "Fail"))
  expect_error(run_qc_pipeline(fit, diag, separation_facets = "Person"), "non-Person")
  expect_error(run_qc_pipeline(fit, diag, separation_facets = "missing"), "non-Person")
  stale <- qc
  stale$verdicts$AffectsOverall <- NULL
  expect_error(print(stale), "Recreate it with run_qc_pipeline", fixed = TRUE)
  expect_error(summary(stale), "Recreate it with run_qc_pipeline", fixed = TRUE)
  expect_error(plot_qc_pipeline(stale, draw = FALSE), "Recreate it with run_qc_pipeline", fixed = TRUE)
  sx <- summary(qc)
  sx$verdicts$AffectsOverall <- NULL
  expect_error(print(sx), "Recreate it with run_qc_pipeline", fixed = TRUE)
  # A non-rater facet must not be selected just because it appears first.
  local_mocked_bindings(
    infer_default_rater_facet = function(facet_names, fallback_first = TRUE) {
      if (fallback_first) facet_names[1] else NULL
    }, .package = "mfrmr"
  )
  no_rater <- run_qc_pipeline(fit, diag)
  expect_identical(no_rater$verdicts$Verdict[9], "Skip")
  expect_false(no_rater$verdicts$AffectsOverall[9])
  explicit_rater <- run_qc_pipeline(fit, diag, rater_facet = "Rater")
  expect_true(explicit_rater$verdicts$AffectsOverall[9])
})

test_that("QC cannot replace unavailable diagnostics with passing values", {
  fit <- make_qc_study1_fit()
  diag <- diagnose_mfrm(fit, residual_pca = "none")
  diag$overall_fit$Infit <- NA_real_
  diag$fit$Outfit[1] <- NA_real_
  diag$unexpected$summary <- data.frame()
  diag$subsets$summary <- data.frame()
  qc <- run_qc_pipeline(fit, diag)
  expect_true(all(qc$verdicts$Verdict[c(2, 5, 6, 8)] == "Warn"))
  expect_true(is.na(qc$details$global_fit$infit))
  expect_true(is.na(qc$details$element_misfit$misfit_pct))
  expect_true(is.na(qc$details$unexpected$unexpected_pct))
  expect_true(is.na(qc$details$connectivity$n_subsets))
  diag$fit <- data.frame()
  expect_identical(run_qc_pipeline(fit, diag)$verdicts$Verdict[5], "Warn")
  mismatched <- fit
  mismatched$facets$others$Estimate[1] <- mismatched$facets$others$Estimate[1] + 1
  expect_error(run_qc_pipeline(mismatched, diag), "fit/diagnostics mismatch", fixed = TRUE)
})

test_that("QC checks step order within ladders in numeric step order", {
  fit <- make_qc_study1_fit()
  diag <- diagnose_mfrm(fit, residual_pca = "none")
  # Isolate the category screen with two ordered ladders whose concatenation drops.
  fit$steps <- data.frame(
    StepFacet = rep(c("C1", "C2"), each = 3),
    Step = rep(c("Step_1", "Step_2", "Step_3"), 2),
    Estimate = c(-1, 0, 1, -2, 0, 2)
  )
  fit$steps <- fit$steps[c(3, 2, 1, 6, 5, 4), ]
  expect_error(run_qc_pipeline(fit, diag), "fit/diagnostics mismatch", fixed = TRUE)
  # This synthetic screen fixture records each mocked step table, not a refit.
  diag$parameter_uncertainty$source <- mfrm_native_location_identity(fit)
  qc <- run_qc_pipeline(fit, diag)
  expect_true(qc$details$category_structure$ordered)
  fit$steps$Estimate[1] <- -3
  diag$parameter_uncertainty$source <- mfrm_native_location_identity(fit)
  expect_false(run_qc_pipeline(fit, diag)$details$category_structure$ordered)
  fit$steps$Estimate[1] <- NA_real_
  diag$parameter_uncertainty$source <- mfrm_native_location_identity(fit)
  missing <- run_qc_pipeline(fit, diag)
  expect_true(is.na(missing$details$category_structure$ordered))
  expect_false(missing$verdicts$Verdict[7] == "Pass")
})
