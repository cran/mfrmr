test_that("run_qc_pipeline returns correct class and structure", {
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)
  qc <- run_qc_pipeline(fit)

  expect_s3_class(qc, "mfrm_qc_pipeline")
  expect_true(is.list(qc))
  expect_named(qc, c("verdicts", "overall", "details", "recommendations", "config"),
               ignore.order = TRUE)
})

test_that("verdicts table has 10 rows with valid verdict values", {
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)
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
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)
  qc <- run_qc_pipeline(fit)

  expect_true(qc$overall %in% c("Pass", "Warn", "Fail"))
})

test_that("recommendations is a character vector", {
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)
  qc <- run_qc_pipeline(fit)

  expect_type(qc$recommendations, "character")
})

test_that("print.mfrm_qc_pipeline produces output", {
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)
  qc <- run_qc_pipeline(fit)

  out <- capture.output(print(qc))
  expect_true(length(out) > 0)
  expect_true(any(grepl("QC Pipeline", out)))
  expect_true(any(grepl("Overall:", out)))
})

test_that("summary.mfrm_qc_pipeline produces correct output", {
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)
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
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)
  qc <- run_qc_pipeline(fit)

  vt <- plot_qc_pipeline(qc, draw = FALSE)
  expect_true(is.data.frame(vt) || tibble::is_tibble(vt))
  expect_equal(nrow(vt), 10)
})

test_that("threshold_profile strict and lenient work", {
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)

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
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)

  qc <- run_qc_pipeline(fit, thresholds = list(reliability_pass = 0.99))
  expect_equal(qc$config$thresholds$reliability_pass, 0.99)
})

test_that("run_qc_pipeline works with pre-computed diagnostics", {
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)
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

test_that("run_qc_pipeline does not convert category-count failures into pass verdicts", {
  toy <- load_mfrmr_data("study1")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 25)
  diag <- diagnose_mfrm(fit, residual_pca = "none")
  diag$obs <- structure(list(), class = "bogus_obs")

  qc <- run_qc_pipeline(fit, diagnostics = diag)

  expect_identical(as.character(qc$verdicts$Verdict[7]), "Skip")
  expect_match(qc$verdicts$Detail[7], "Category counts could not be computed", fixed = TRUE)
  expect_true(nzchar(qc$details$category_structure$error))
  expect_false(identical(qc$overall, "Pass"))
})

test_that("run_qc_pipeline surfaces incomplete bias collections as warn rather than skip/pass", {
  toy <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    toy, "Person", c("Rater", "Criterion"), "Score",
    method = "JML", maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  diag$interactions <- data.frame()
  bias_collection <- structure(
    list(
      by_pair = list(),
      errors = data.frame(
        Interaction = "Rater x Criterion",
        Facets = "Rater x Criterion",
        Error = "forced pair failure",
        stringsAsFactors = FALSE
      )
    ),
    class = c("mfrm_bias_collection", "mfrm_bundle", "list")
  )

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
