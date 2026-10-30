test_that("reliability wording survives report, APA and console routes", {
  fit <- make_toy_fit(); diagnostics <- make_toy_diagnostics(fit)
  precision <- precision_review_report(fit, diagnostics = diagnostics)
  results <- mfrm_results(fit, include = c("fit", "precision"), diagnostics = diagnostics)
  report <- mfrm_report(results)
  apa <- build_apa_reporting_contract(fit, diagnostics)$report_text
  console <- paste(capture.output(print(summary(diagnostics))), collapse = " ")
  # Check each independent reporting route, not just the shared sentence helper.
  for (text in list(precision$fit_separation_basis$Interpretation,
      report$precision_evidence_summary$Boundary,
      report$precision_reporting_templates$Caveat, unlist(apa), console)) {
    text <- paste(text, collapse = " ")
    text <- gsub("[[:space:]]+", " ", text)
    expect_match(text, "high rater values indicate distinguishable severities")
    expect_match(text, "Low rater values can also reflect imprecision")
    expect_match(text, "approximate compatibility indices")
    expect_match(text, "not full-covariance joint tests")
    expect_match(text, "not discovered groups")
  }
  expect_identical(results$diagnostics$reliability, diagnostics$reliability)
})
