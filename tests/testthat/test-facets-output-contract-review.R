test_that("facets_output_contract_review returns full-coverage review in facets branch", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- suppressWarnings(
    mfrmr::fit_mfrm(
      data = d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "JML",
      model = "RSM",
      maxit = 20
    )
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")
  bias <- mfrmr::estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Task", max_iter = 2)

  contract_review <- mfrmr::facets_output_contract_review(
    fit = fit,
    diagnostics = diag,
    bias_results = bias,
    branch = "facets"
  )

  expect_s3_class(contract_review, "mfrm_facets_contract_review")
  expect_true(all(c(
    "overall", "column_summary", "column_review", "missing_preview",
    "metric_summary", "metric_by_table", "metric_checks", "settings"
  ) %in% names(contract_review)))
  expect_true(is.data.frame(contract_review$overall))
  expect_true(is.data.frame(contract_review$column_review))
  expect_true(is.data.frame(contract_review$metric_checks))
  expect_identical(contract_review$settings$intended_use, "facets_output_contract_review")
  expect_false(isTRUE(contract_review$settings$external_validation))

  expect_equal(contract_review$overall$ColumnMismatches[1], 0)
  expect_equal(contract_review$overall$ColumnMismatchRate[1], 0)
  expect_equal(contract_review$overall$MeanColumnCoverage[1], 1)
  expect_equal(contract_review$overall$MinColumnCoverage[1], 1)
  expect_equal(contract_review$overall$MeanColumnCoverageAvailable[1], 1)
  expect_equal(contract_review$overall$MinColumnCoverageAvailable[1], 1)
  expect_true(contract_review$overall$MetricFailed[1] <= 0)
})

test_that("facets_output_contract_review integrates with summary() and plot()", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- suppressWarnings(
    mfrmr::fit_mfrm(
      data = d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "JML",
      model = "RSM",
      maxit = 20
    )
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")
  bias <- mfrmr::estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Task", max_iter = 2)

  contract_review <- mfrmr::facets_output_contract_review(
    fit = fit,
    diagnostics = diag,
    bias_results = bias,
    branch = "original"
  )

  s <- summary(contract_review)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_identical(s$summary_kind, "mfrm_facets_contract_review")

  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "mfrmr FACETS Output Contract Review Summary", fixed = TRUE)

  p1 <- plot(contract_review, draw = FALSE)
  p2 <- plot(contract_review, type = "table_coverage", draw = FALSE)
  p3 <- plot(contract_review, type = "metric_status", draw = FALSE)
  p4 <- plot(contract_review, type = "metric_by_table", draw = FALSE)

  expect_s3_class(p1, "mfrm_plot_data")
  expect_s3_class(p2, "mfrm_plot_data")
  expect_s3_class(p3, "mfrm_plot_data")
  expect_s3_class(p4, "mfrm_plot_data")
})

test_that("facets_output_contract_review contract coverage includes unavailable rows", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- suppressWarnings(
    mfrmr::fit_mfrm(
      data = d,
      person = "Person",
      facets = c("Rater"),
      score = "Score",
      method = "JML",
      model = "RSM",
      maxit = 20
    )
  )

  contract_review <- mfrmr::facets_output_contract_review(
    fit = fit,
    branch = "facets",
    include_metrics = FALSE
  )

  expect_gt(contract_review$overall$ColumnMismatches[1], 0)
  expect_lt(contract_review$overall$MeanColumnCoverage[1], 1)
  expect_lt(contract_review$overall$MinColumnCoverage[1], 1)
  expect_equal(contract_review$overall$MeanColumnCoverageAvailable[1], 1)
  expect_equal(contract_review$overall$MinColumnCoverageAvailable[1], 1)
})

test_that("reference_case_review exposes package-native review wording", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- suppressWarnings(
    mfrmr::fit_mfrm(
      data = d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "JML",
      model = "RSM",
      maxit = 20
    )
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")

  review <- mfrmr::reference_case_review(
    fit = fit,
    diagnostics = diag,
    reference_profile = "core"
  )

  expect_s3_class(review, "mfrm_reference_review")
  expect_true(all(c(
    "overall", "component_summary", "attention_items",
    "metric_summary", "metric_checks", "settings", "contract_review"
  ) %in% names(review)))
  expect_identical(as.character(review$overall$ReferenceProfile[1]), "core")
  expect_identical(as.character(review$overall$CompatibilityLayer[1]), "package-native")
  expect_identical(review$settings$intended_use, "reference_contract_review")
  expect_false(isTRUE(review$settings$external_validation))

  s <- summary(review)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_true(is.data.frame(s$overview))
  expect_identical(as.character(s$overview$Class[1]), "mfrm_reference_review")

  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "mfrmr Reference Review Summary", fixed = TRUE)
})
