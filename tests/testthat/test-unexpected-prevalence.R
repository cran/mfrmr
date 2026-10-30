test_that("unexpected-response prevalence is independent of the display limit", {
  toy <- load_mfrmr_data("example_operational")
  # Replicate the synthetic person profiles to exercise the diagnostic 100-row cap.
  toy <- do.call(rbind, lapply(1:3, function(i) {
    copy <- toy
    copy$Person <- paste0(copy$Person, "_", i)
    copy
  }))
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "MML", model = "RSM")
  diagnostics <- diagnose_mfrm(fit, residual_pca = "none")

  full <- unexpected_response_table(fit, diagnostics, top_n = nrow(toy),
                                    abs_z_min = 2, rule = "either")
  preview <- unexpected_response_table(fit, diagnostics, top_n = 1,
                                       abs_z_min = 2, rule = "either")
  expect_gt(full$summary$UnexpectedN, 100L)
  expect_equal(nrow(preview$table), 1L)
  expect_equal(preview$table, head(full$table, 1))
  expect_equal(preview$summary, full$summary)
  expect_equal(full$summary$UnexpectedN, nrow(full$table))
  expect_equal(full$summary$UnexpectedPercent, 100 * nrow(full$table) / nrow(toy))
  default_full <- unexpected_response_table(fit, diagnostics, top_n = nrow(toy))
  expect_equal(nrow(diagnostics$unexpected$table), min(100L, nrow(default_full$table)))
  expect_equal(diagnostics$unexpected$summary, default_full$summary)

  for (rule in c("either", "both")) {
    all_rows <- unexpected_response_table(
      fit, diagnostics, abs_z_min = 1.5, prob_max = 0.4,
      top_n = nrow(toy), rule = rule
    )
    one_row <- unexpected_response_table(
      fit, diagnostics, abs_z_min = 1.5, prob_max = 0.4,
      top_n = 1, rule = rule
    )
    expect_equal(one_row$summary, all_rows$summary)
    expect_equal(all_rows$summary$LowProbabilityN, sum(all_rows$table$FlagLowProbability))
    expect_equal(all_rows$summary$LargeResidualN, sum(all_rows$table$FlagLargeResidual))
  }

  bias <- estimate_bias(fit, diagnostics, facet_a = "Rater", facet_b = "Criterion",
                        max_iter = 2)
  # With zero adjustment the before/after count must agree, even for top_n = 1.
  bias$table[["Bias Size"]] <- 0
  adjusted <- unexpected_after_bias_table(fit, bias, diagnostics, top_n = 1,
                                          abs_z_min = 2, rule = "either")
  expect_equal(nrow(adjusted$table), 1L)
  expect_equal(adjusted$summary$BaselineUnexpectedN, full$summary$UnexpectedN)
  expect_equal(adjusted$summary$AfterBiasUnexpectedN, full$summary$UnexpectedN)
  expect_equal(adjusted$summary$UnexpectedPercent, full$summary$UnexpectedPercent)
  expect_equal(adjusted$summary$ReducedBy, 0)
  expect_equal(adjusted$summary$ReducedPercent, 0)

  none <- unexpected_response_table(fit, diagnostics, abs_z_min = Inf,
                                    prob_max = 0, top_n = 1)
  expect_equal(nrow(none$table), 0L)
  expect_equal(none$summary$UnexpectedN, 0L)
  expect_equal(none$summary$UnexpectedPercent, 0)

  diagnostics$obs$StdResidual[1] <- NA_real_
  partial <- unexpected_after_bias_table(fit, bias, diagnostics, prob_max = 0)
  expect_equal(partial$summary$BaselineUnavailableObservations, 1)
  expect_equal(partial$summary$UnavailableObservations, 0)
  expect_true(is.na(partial$summary$ReducedBy))
  expect_true(is.na(partial$summary$ReducedPercent))
  expect_error(plot(partial, type = "comparison", draw = FALSE), "some responses could not be classified")
  legacy <- preview
  legacy$summary$UnavailableObservations <- NULL
  expect_error(summary(legacy), "Recreate this unexpected-response result")
  expect_error(plot(legacy, draw = FALSE), "Recreate this unexpected-response result")
  legacy_summary <- summary(preview)
  legacy_summary$summary$UnavailableObservations <- NULL
  expect_error(print(legacy_summary), "Recreate this unexpected-response result")
})
