test_that("agreement prose distinguishes zero, missing expectation and no overlap", {
  summary <- data.frame(RaterFacet = "Judge", Pairs = 3, AvailablePairs = 2,
    OpportunityCount = 8, ExactAgreement = 0, AdjacentAgreement = .75,
    ExpectedExactAgreement = NA_real_, RepeatedCells = 1,
    ExpectedUnavailableContexts = 2)
  diagnostics <- list(interrater = list(summary = summary))
  text <- reporting_interrater_text(diagnostics)
  expect_match(text, "agreement for Judge")
  expect_match(text, "8 matched-context opportunities")
  expect_match(text, "2 of 3 candidate rater pairs")
  expect_match(text, "Observed exact agreement = 0.00", fixed = TRUE)
  expect_match(text, "model-expected exact agreement = unavailable", fixed = TRUE)
  expect_match(text, "absolute score difference at most 1, including exact matches", fixed = TRUE)
  expect_match(text, "averaged using observation weights")
  expect_match(text, "unavailable for 2 matched opportunities")
  expect_match(text, "not chance agreement or an inferential test")
  expect_match(text, "lacks the primary-estimate availability audit")

  diagnostics$interrater$summary$OpportunityCount <- 0
  diagnostics$interrater$summary$AvailablePairs <- 0
  diagnostics$interrater$summary$ExactAgreement <- NA_real_
  diagnostics$interrater$summary$AdjacentAgreement <- NA_real_
  diagnostics$interrater$summary$RepeatedCells <- 0
  diagnostics$interrater$summary$ExpectedUnavailableContexts <- 0
  text <- reporting_interrater_text(diagnostics)
  expect_match(text, "agreement is unavailable, not zero")
  expect_match(text, "Observed exact agreement = unavailable", fixed = TRUE)
  expect_false(grepl("averaged|unavailable for 2", text))
})

test_that("absent and older agreement summaries do not acquire invented values", {
  expect_identical(reporting_interrater_text(NULL), "")
  expect_match(reporting_interrater_text(list()), "No inter-rater agreement summary was stored")
  old <- list(interrater = list(summary = data.frame(
    RaterFacet = "Judge", TotalPairs = 8, ExactAgreement = .5,
    ExpectedExactAgreement = .25, AdjacentAgreement = 1)))
  text <- reporting_interrater_text(old)
  expect_match(text, "8 matched-context opportunities")
  expect_match(text, "model-expected exact agreement = 0.25", fixed = TRUE)
  old$interrater$summary$TotalPairs <- NULL
  text <- reporting_interrater_text(old)
  expect_match(text, "denominator was not stored")
  expect_false(grepl("candidate rater pairs|no matched opportunities", text, ignore.case = TRUE))
})

test_that("agreement values and interpretation survive APA, reports and saved export", {
  fit <- readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits[[1L]]
  diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  local_mocked_bindings(fit_mfrm = function(...) stop("No refit"),
    diagnose_mfrm = function(...) stop("No new diagnostics"),
    calc_interrater_agreement = function(...) stop("No new agreement"),
    compute_prob_matrix = function(...) stop("No new probabilities"), .package = "mfrmr")
  before <- serialize(diagnostics, NULL)
  apa <- build_apa_outputs(fit, diagnostics, context = list(rater_facet = "Criterion"))
  text <- apa$contract$summaries$interrater_sentence
  expect_match(text, "agreement for Rater")
  expect_false(grepl("agreement for Criterion", text, fixed = TRUE))
  expect_match(paste(apa$report_text, collapse = " "), "model-expected exact agreement")
  expect_match(apa$table_figure_notes, "matched-context opportunities")
  result <- mfrm_results(fit, diagnostics = diagnostics,
    include = c("fit", "diagnostics", "apa"), compute = "never")
  saved <- unserialize(serialize(result, NULL))
  for (style in c("qc", "apa", "validation", "reviewer", "technical")) {
    report <- mfrm_report(saved, style = style)
    expect_identical(report$tables$interrater_summary, as.data.frame(diagnostics$interrater$summary))
    expect_identical(report$tables$interrater_interpretation$Text, text)
    expect_match(report$markdown, text, fixed = TRUE)
    expect_match(paste(report$precision_reporting_templates$Caveat, collapse = " "),
      "high rater values indicate distinguishable severities")
  }
  html <- mfrm_report(saved, output = "html")
  on.exit(unlink(html$path), add = TRUE)
  expect_match(html$html, text, fixed = TRUE)
  path <- tempfile("agreement-report-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  exported <- export_mfrm_results(saved, output_dir = path, include = "report",
    acknowledge_sensitive = TRUE)
  files <- exported$written_files
  note_file <- files$Path[files$Component == "report_interrater_interpretation"]
  expect_length(note_file, 1)
  expect_identical(read.csv(note_file)$Text, text)
  values <- read.csv(files$Path[files$Component == "report_interrater_summary"])
  expect_equal(values, as.data.frame(diagnostics$interrater$summary), tolerance = 1e-14)
  expect_identical(serialize(diagnostics, NULL), before)

  saved$diagnostics$interrater <- list(summary = data.frame())
  expect_match(mfrm_report(saved)$markdown, "No inter-rater agreement summary was stored")
  saved$diagnostics <- NULL
  expect_false("interrater_interpretation" %in% names(mfrm_report(saved)$tables))
  fit_only <- mfrm_results(fit, include = "fit", compute = "never")
  expect_false("interrater_interpretation" %in% names(mfrm_report(fit_only)$tables))
  expect_true(any(mfrm_report(fit_only)$sections$Status == "not_computed"))
})
