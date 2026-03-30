test_that("reference_case_benchmark returns package-native benchmark bundle", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = c("synthetic_truth", "synthetic_bias_contract", "study1_itercal_pair"),
    method = "MML",
    quad_points = 5,
    maxit = 30
  ))

  expect_s3_class(bench, "mfrm_reference_benchmark")
  expect_true(all(c(
    "overview", "summary", "table", "fit_runs", "case_summary",
    "design_checks", "recovery_checks", "bias_checks", "pair_checks",
    "linking_checks", "source_profile", "settings", "notes"
  ) %in% names(bench)))

  expect_true(is.data.frame(bench$overview))
  expect_true(is.data.frame(bench$fit_runs))
  expect_true(is.data.frame(bench$case_summary))
  expect_true(is.data.frame(bench$design_checks))
  expect_true(is.data.frame(bench$recovery_checks))
  expect_true(is.data.frame(bench$bias_checks))
  expect_true(is.data.frame(bench$pair_checks))
  expect_true(is.data.frame(bench$linking_checks))
  expect_true(is.data.frame(bench$source_profile))
  expect_identical(bench$settings$intended_use, "internal_benchmark")
  expect_false(isTRUE(bench$settings$external_validation))

  expect_equal(nrow(bench$case_summary), 3)
  expect_true(all(c("synthetic_truth", "synthetic_bias_contract", "study1_itercal_pair") %in% bench$case_summary$Case))
  expect_true(all(bench$fit_runs$PrecisionTier %in% c("model_based", "hybrid", "exploratory")))
})

test_that("reference_case_benchmark recovers synthetic truth under MML", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_truth",
    method = "MML",
    quad_points = 7,
    maxit = 40
  ))

  expect_true(all(bench$recovery_checks$Status %in% c("Pass", "Warn", "Fail")))
  expect_true(min(bench$recovery_checks$Correlation, na.rm = TRUE) > 0.95)
  expect_true(max(bench$recovery_checks$MeanAbsoluteDeviation, na.rm = TRUE) < 0.30)
  expect_true(all(bench$fit_runs$SupportsFormalInference))
})

test_that("reference_case_benchmark captures pair stability and summary output", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "study1_itercal_pair",
    method = "MML",
    quad_points = 5,
    maxit = 30
  ))

  expect_true(all(c("Criterion", "Rater", "OverallFit") %in% bench$pair_checks$Facet))
  expect_true(all(c("Criterion", "Rater") %in% bench$linking_checks$Facet))
  expect_identical(
    as.character(bench$linking_checks$Status[bench$linking_checks$Facet == "Rater"][1]),
    "Pass"
  )
  expect_identical(
    as.character(bench$linking_checks$Status[bench$linking_checks$Facet == "Criterion"][1]),
    "Warn"
  )
  criterion_row <- bench$pair_checks[bench$pair_checks$Facet == "Criterion", , drop = FALSE]
  expect_true(criterion_row$Pearson[1] > 0.95)
  expect_true(criterion_row$MeanAbsoluteDifference[1] < 0.10)

  s <- summary(bench)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_true(is.data.frame(s$overview))
  expect_identical(as.character(s$overview$Class[1]), "mfrm_reference_benchmark")

  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "mfrmr Internal Benchmark Summary", fixed = TRUE)
})

test_that("reference_case_benchmark verifies bias-contract identities", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_bias_contract",
    method = "MML",
    quad_points = 7,
    maxit = 40
  ))

  expect_true(nrow(bench$bias_checks) >= 4)
  stable_rows <- bench$bias_checks[
    bench$bias_checks$Metric %in% c("BiasDFIdentity", "LocalMeasureIdentity", "PairContrastIdentity"),
    ,
    drop = FALSE
  ]
  expect_true(nrow(stable_rows) == 3)
  expect_true(all(stable_rows$Status == "Pass"))
  expect_true(max(stable_rows$MaxError, na.rm = TRUE) < 1e-8)
  expect_true(any(bench$source_profile$RuleID == "bias_pairwise_welch"))
})

test_that("reference_case_benchmark reports exploratory precision under JML", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_truth",
    method = "JML",
    maxit = 30
  ))

  expect_true(all(bench$fit_runs$PrecisionTier == "exploratory"))
  expect_true(all(!bench$fit_runs$SupportsFormalInference))
})

test_that("reference_case_benchmark handles single-case runs without warnings", {
  expect_no_warning(
    reference_case_benchmark(
      cases = "synthetic_truth",
      method = "MML",
      quad_points = 5,
      maxit = 20
    )
  )
})

test_that("reference_case_benchmark does not pass bias-contract cases when checks are missing", {
  local_mocked_bindings(
    estimate_bias = function(...) {
      structure(list(table = data.frame()), class = "mfrm_bias")
    },
    .package = "mfrmr"
  )

  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_bias_contract",
    method = "MML",
    quad_points = 5,
    maxit = 20
  ))

  row <- bench$case_summary[bench$case_summary$Case == "synthetic_bias_contract", , drop = FALSE]
  expect_identical(as.character(row$Status[1]), "Warn")
  expect_match(row$KeySignal[1], "No bias-contract checks were produced", fixed = TRUE)
})
