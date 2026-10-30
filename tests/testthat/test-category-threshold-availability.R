test_that("threshold order uses adjacent indices within complete families", {
  steps <- data.frame(Step = paste0("Step_", 1:4), Estimate = c(-1, NA, 1, 2))
  ordered <- calc_step_order(steps, expected_steps = 4)
  expect_equal(ordered$Spacing, c(NA, NA, NA, 1))
  coverage <- summarize_threshold_order(ordered)
  expect_equal(coverage$Available, 1)
  expect_equal(coverage$Unavailable, 2)
  expect_true(is.na(coverage$ThresholdMonotonic))
  expect_match(summarize_step_estimates(steps, 4), "2 of 3 comparisons unavailable")
  removed <- calc_step_order(steps[-2, ], 4)
  expect_equal(removed$Spacing, c(NA, NA, 1))
  expect_equal(summarize_threshold_order(removed)$Unavailable, 2)
  expect_true(is.na(summarize_threshold_order(calc_step_order(steps[-4, ], 4))$ThresholdMonotonic))
  steps$Estimate <- c(0, 0, 1, 2)
  expect_true(summarize_threshold_order(calc_step_order(steps))$ThresholdMonotonic)
  many <- data.frame(StepFacet = rep(c("B", "A"), each = 10),
                     Step = rep(paste0("Step_", 1:10), 2), Estimate = rep(1:10, 2))
  expect_equal(summarize_threshold_order(calc_step_order(many[nrow(many):1, ]))$Available, 18)
  expect_true(summarize_threshold_order(calc_step_order(many))$ThresholdMonotonic)
  missing_family <- summarize_threshold_order(calc_step_order(many, 10, c("A", "B", "C")))
  expect_equal(missing_family$Unavailable, 9)
  expect_true(is.na(missing_family$ThresholdMonotonic))
  for (labels in list(c("Step_1", "Step_1"), c("first", "second"), c("Step_-1", "Step_2"))) {
    unknown <- summarize_threshold_order(calc_step_order(data.frame(Step = labels, Estimate = c(0, 1)), 2))
    expect_equal(unknown$Available, 0)
    expect_true(is.na(unknown$ThresholdMonotonic))
  }
  steps$Estimate <- c(1, 0, NA, 2)
  partial <- summarize_threshold_order(calc_step_order(steps))
  expect_false(partial$ThresholdMonotonic)
  expect_equal(partial$Decreasing, 1)
  binary <- summarize_threshold_order(calc_step_order(data.frame(Step = "Step_1", Estimate = 0), 1))
  expect_true(binary$NotApplicable)
  expect_true(is.na(binary$ThresholdMonotonic))
  expect_match(threshold_order_note(binary), "not applicable")
  expect_false(summarize_threshold_order(calc_step_order(data.frame(Step = "Step_1", Estimate = NA), 1))$NotApplicable)
  expect_true(is.na(summarize_threshold_order(calc_step_order(NULL))$ThresholdMonotonic))
})

test_that("category counts and fit flags preserve unavailable inputs", {
  obs <- data.frame(Observed = c(0, 0, 1, 1), Weight = c(1, 2, 1, 0.5),
                    PersonMeasure = 0, Expected = 0.5, StdSq = 1, Var = 0.25, Residual = 0)
  complete <- calc_category_stats(obs)
  expect_equal(complete$Count, c(3, 1.5))
  expect_equal(complete$Infit, c(1, 1))
  expect_true(all(is.na(complete$InfitFlag)))
  obs$StdSq[1] <- NA
  incomplete_fit <- calc_category_stats(obs)
  expect_true(is.na(incomplete_fit$Infit[1]))
  expect_true(is.na(incomplete_fit$InfitFlag[1]))
  expect_equal(incomplete_fit$Count, complete$Count)
  obs$Weight[1] <- NA
  unavailable <- calc_category_stats(obs)
  expect_true(all(is.na(unavailable$Count)))
  expect_true(all(is.na(unavailable$LowCount)))
  expect_true(all(is.na(unavailable$OutfitFlag)))
  expect_true(is.na(summarize_category_usage(unavailable)$UnusedCategories))
  expect_match(category_usage_note(summarize_category_usage(unavailable)), "0 of 2")
  expect_true(is.na(summarize_category_usage(NULL)$UsedCategories))
  expect_match(category_usage_note(summarize_category_usage(NULL)), "not available")
})

test_that("missing expected probability is not converted to a zero expected count", {
  fit <- make_toy_fit()
  # Use the native observation count rather than requiring a particular prep slot.
  probs <- matrix(1 / fit$config$n_cat, nrow = nrow(make_toy_diagnostics(fit)$obs), ncol = fit$config$n_cat)
  probs[1, 1] <- NA_real_
  testthat::local_mocked_bindings(compute_response_probability_bundle = function(...) list(probs = probs))
  expected <- calc_expected_category_counts(fit)
  expect_true(is.na(expected$ExpectedCount[1]))
  expect_true(all(is.na(expected$ExpectedPercent)))
  expect_true(all(is.finite(expected$ExpectedCount[-1])))
})

test_that("reports do not infer adequate usage or ordering from empty summaries", {
  fit <- make_toy_fit()
  diag <- make_toy_diagnostics(fit)
  testthat::local_mocked_bindings(calc_category_stats = function(...) tibble::tibble(),
                                 calc_step_order = function(...) tibble::tibble())
  contract <- build_apa_reporting_contract(fit, diag)
  text <- gsub("[[:space:]]+", " ", paste(contract$report_text, collapse = " "))
  expect_match(text, "Category counts were not available")
  expect_match(text, "no ordering conclusion")
  expect_false(grepl("category usage was adequate|thresholds were ordered|no disordered steps", text))
  expect_true(is.na(contract$summaries$unused_categories))
})

test_that("rating and category summaries expose coverage without internal mode codes", {
  fit <- make_toy_fit()
  diag <- make_toy_diagnostics(fit)
  rating <- rating_scale_table(fit, diag)
  expect_equal(rating$category_usage$Categories, fit$config$n_cat)
  expect_equal(rating$summary$ThresholdComparisons, fit$config$n_cat - 2)
  expect_equal(rating$summary$UnavailableThresholdComparisons, 0)
  s <- summary(rating)
  printed <- paste(capture.output(print(s)), collapse = " ")
  expect_false(grepl("legacy_plugin|DiagnosticMode|ExpectedCountBasis|CIEligible", printed))
  expect_match(printed, "Adjacent threshold comparisons")
  bad_fit <- fit
  bad_fit$facets$others$Estimate[1] <- bad_fit$facets$others$Estimate[1] + 1
  expect_error(rating_scale_table(bad_fit, diag), "fit/diagnostics mismatch", fixed = TRUE)
  expect_error(category_structure_report(bad_fit, diag), "fit/diagnostics mismatch", fixed = TRUE)
  rating$threshold_coverage <- NULL
  expect_error(summary(rating), "Recreate it")
  expect_error(plot(rating, draw = FALSE), "Recreate it")
  s$category_usage <- NULL
  expect_error(print(s), "Recreate it")
  cat_tbl <- calc_category_stats(diag$obs, fit)
  cat_tbl$Count[1] <- 0
  for (nm in c("LowCount", "InfitFlag", "OutfitFlag", "ZSTDFlag")) cat_tbl[[nm]][] <- NA
  testthat::local_mocked_bindings(calc_category_stats = function(...) cat_tbl)
  dropped <- rating_scale_table(fit, diag, drop_unused = TRUE)
  expect_equal(dropped$summary$Categories, nrow(cat_tbl))
  expect_equal(dropped$summary$DisplayedCategories, nrow(cat_tbl) - 1)
  structure <- category_structure_report(fit, diag, drop_unused = TRUE)
  expect_equal(summary(structure)$summary$Categories, nrow(cat_tbl))
  expect_true(is.na(summary(structure)$summary$FlaggedStats))
  expect_equal(summary(structure)$summary$UnavailableFlagDecisions, nrow(cat_tbl) - 1)
  printed <- paste(capture.output(print(summary(structure))), collapse = " ")
  expect_false(grepl("legacy|DiagnosticMode", printed))
})

test_that("QC distinguishes unavailable ordering from a binary scale", {
  fit <- make_toy_fit()
  diag <- make_toy_diagnostics(fit)
  changed <- fit
  changed$steps$Estimate[2] <- NA_real_
  expect_error(run_qc_pipeline(changed, diag), "fit/diagnostics mismatch")
  # Exercise the QC formatter with controlled threshold summaries, without
  # mixing a changed calibration and diagnostics from the original fit.
  step_order <- calc_step_order(changed$steps)
  order_steps <- calc_step_order
  testthat::local_mocked_bindings(calc_step_order = function(...) step_order)
  qc <- run_qc_pipeline(fit, diag)
  expect_identical(qc$verdicts$Verdict[7], "Warn")
  expect_gt(qc$details$category_structure$threshold_coverage$Unavailable, 0)
  step_order <- order_steps(data.frame(Step = "Step_1", Estimate = 0), expected_steps = 1L)
  # This branch tests the QC formatter using matching count data explicitly.
  testthat::local_mocked_bindings(calc_category_stats = function(...) data.frame(Category = 0:1, Count = c(100, 100)))
  qc <- run_qc_pipeline(fit, diag)
  expect_true(qc$details$category_structure$threshold_coverage$NotApplicable)
  expect_identical(qc$verdicts$Verdict[7], "Pass")
})

test_that("threshold plots connect only available adjacent pairs in one family", {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  segments <- NULL
  testthat::local_mocked_bindings(segments = function(x0, y0, x1, y1, ...) {
    segments <<- data.frame(from = x0, to = x1)
  }, .package = "mfrmr")
  steps <- data.frame(StepFacet = c("A", "A", "A", "B", "B"),
                     Step = c("Step_1", "Step_2", "Step_4", "Step_1", "Step_2"),
                     Estimate = c(0, 1, 2, -1, NA))
  expect_no_error(draw_step_plot(steps))
  expect_equal(segments, data.frame(from = 1L, to = 2L))
  steps$Estimate[] <- NA_real_
  expect_no_error(draw_step_plot(steps))
})
