test_that("zero-reference summaries use FairZ values and retain diagnostic scope", {
  fit <- make_toy_fit(maxit = 10)
  for (style in c("native", "legacy", "both")) {
    b <- fair_average_table(fit, reference = "zero", label_style = style)
    s <- summary(b, digits = 10)
    expect_equal(s$summary$FairMetric, "FairZ")
    expect_true(all(is.finite(s$preview$StandardizedAdjustedAverage)))
    expect_false("AdjustedAverage" %in% names(s$preview))
    expect_false(s$summary$FairCIEligible)
    printed <- capture.output(print(s))
    expect_false(any(grepl("FairCIEligible|FairSEMethod|FairSEStatus|not_requested|qualified formal|full-refit", printed)))
  }
  # Deliberately supplied bundle: counts refer to FairZ, even without settings.
  b <- list(stacked = data.frame(Facet = "Rater", Level = "R1",
    ObservedAverage = 1, StandardizedAdjustedAverage = 1.2,
    StandardizedAdjustedAverageSE = .1,
    StandardizedAdjustedAverageCI_Lower = 1,
    StandardizedAdjustedAverageCI_Upper = 1.4,
    StandardizedAdjustedAverageCI_Level = .95,
    StandardizedAdjustedAverageSEStatus = "ok"), settings = list())
  class(b) <- class(fair_average_table(fit))
  s <- summary(b, digits = 10)
  expect_equal(s$summary$FairSEAvailableRows, 1)
  expect_equal(s$summary$MeanStandardizedAdjustedAverageSE, .1)
  expect_equal(s$summary$StandardizedAdjustedAverageCILevel, .95)
  expect_false(s$preview$FairCIEligible)
  expect_true(any(grepl("should not be used as confidence intervals for decisions", s$notes)))
})

test_that("extreme display measures retain fitted values without borrowing their SEs", {
  d <- simulate_mfrm_data(n_person = 24, n_rater = 4, n_criterion = 3,
    raters_per_person = 4, assignment = "crossed", seed = 803)
  d$Score[d$Person == "P001"] <- 4L
  d$Score[d$Person == "P002"] <- 1L
  for (method in c("JML", "MML")) {
    fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
      method = method, model = if (method == "JML") "PCM" else "RSM", maxit = 120,
      anchors = if (method == "JML") data.frame(Facet = "Person", Level = "P001", Anchor = .75) else NULL))
    dx <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "legacy")
    before <- serialize(fit, NULL)
    base <- fair_average_table(fit, dx)
    b <- fair_average_table(fit, dx, xtreme = .253, umean = 50, uscale = -10)
    p <- b$raw_by_facet$Person
    at <- match(p$Level, fit$facets$person$Person)
    expect_equal(p$PrimaryMeasure, 50 - 10 * fit$facets$person$Estimate[at])
    adjusted <- p$ExtremeAdjustment > 0
    expect_true(any(adjusted))
    expect_true(all(is.finite(p$Measure[adjusted])))
    expect_true(all(is.na(p$ModelSE[adjusted])))
    expect_true(all(is.na(p$RealSE[adjusted])))
    expect_true(all(p$MeasureBasis[adjusted] == "Extreme-score display only"))
    expect_equal(p$ExtremeAdjustment[adjusted], rep(.253, sum(adjusted)))
    for (metric in c("FairM", "FairZ")) {
      expect_equal(p[[metric]], base$raw_by_facet$Person[[metric]][match(p$Level, base$raw_by_facet$Person$Level)])
    }
    expect_identical(serialize(fit, NULL), before)
    if (method == "JML") {
      expect_equal(p$PrimaryMeasure[p$Level == "P001"], 42.5)
      expect_equal(p$PrimaryMeasure[p$Level == "P002"], Inf)
      expect_match(summary(fit)$estimation_note, "without an extreme-score adjustment", fixed = TRUE)
      printed <- gsub("[[:space:]]+", " ", paste(capture.output(print(fit)), collapse = " "))
      expect_match(printed, "finite-item bias correction", fixed = TRUE)
      expect_false(grepl("ready_with_exclusions|joint_person_coordinate_scale|descriptive_jml", printed))
      expect_true(all(is.finite(b$raw_by_facet$Rater$FairM)))
      expect_true(all(is.finite(b$raw_by_facet$Rater$FairZ)))
      expect_match(b$raw_by_facet$Rater$FairMReference[1], "finite fitted Person", fixed = TRUE)
      expect_no_error(plot_fair_average(b, facet = "Rater", draw = FALSE))
      expect_gt(b$raw_by_facet$Rater$FairMExcludedPersons[1], 0)
      expect_no_error(plot_fair_average(b, facet = "Rater", metric = "FairZ", draw = FALSE))
    } else {
      expect_true(all(is.finite(p$PrimaryMeasure)))
      expect_length(summary(fit)$estimation_note, 0L)
      expect_true(all(is.finite(base$raw_by_facet$Person$ModelSE)))
    }
    csv <- tempfile(fileext = ".csv")
    write.csv(b$stacked, csv, row.names = FALSE)
    exported <- read.csv(csv)
    expect_true(all(c("PrimaryMeasure", "MeasureBasis", "ExtremeAdjustment") %in% names(exported)))
    expect_true(all(exported$ExtremeAdjustment[exported$ExtremeAdjustment > 0] == .253))
    unlink(csv)
    expect_true(any(grepl("do not correct JML bias", summary(b)$notes, fixed = TRUE)))
    older_summary <- summary(b)
    older_summary$fair_measure_basis_recorded <- NULL
    expect_error(print(older_summary), "Recreate summary(fair_average_table", fixed = TRUE)
    older <- b
    older$raw_by_facet <- lapply(older$raw_by_facet, function(tbl) {
      tbl$PrimaryMeasure <- NULL
      tbl
    })
    expect_error(summary(older), "Recreate fair_average_table", fixed = TRUE)
    expect_error(plot_fair_average(older, draw = FALSE), "Recreate fair_average_table", fixed = TRUE)
    # Older unadjusted bundles also lack the original mean-reference basis.
    older$settings$xtreme <- 0
    expect_error(summary(older), "Recreate fair_average_table", fixed = TRUE)
    older_dx <- dx
    older_dx$fair_average$raw_by_facet$Person$FairMReference <- NULL
    expect_error(summary(older_dx), "Recompute diagnose_mfrm(fit)", fixed = TRUE)
  }
})

test_that("incomplete gradients cannot create a finite fair-score SE", {
  fit <- make_toy_fit(maxit = 10)
  fit$config$model <- "GPCM"
  raw <- list(Rater = data.frame(Level = "R1", FairM = 1, FairZ = 1))
  testthat::local_mocked_bindings(
    finite_difference_gradient = function(...) c(1, NA_real_), .package = "mfrmr")
  got <- mfrmr:::add_gpcm_fair_average_delta_se(raw, fit,
    covariance = list(status = "ok", cov = diag(2)))$Rater
  expect_true(is.na(got$FairMSE))
  expect_true(is.na(got$FairZSE))
  expect_match(got$FairZ_SE_Detail, "complete")
})

test_that("fair-score plot intervals cannot acquire formal eligibility", {
  fit <- make_toy_fit(maxit = 10)
  p <- plot_fair_average(fit, show_ci = TRUE, draw = FALSE,
    show_title = FALSE, show_notes = FALSE)
  printed <- capture.output(print(p))
  expect_false(any(grepl("mfrm_plot_data>|CI_Eligible|diagnostic_only|gpcm_boundary", printed)))
  d <- p$data$data
  expect_true(any(is.finite(d$CI_Lower)))
  expect_match(p$data$subtitle, "observed means held fixed", fixed = TRUE)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    visible <- plot_fair_average(fit, show_ci = TRUE, draw = FALSE)
    expect_match(gsub("[[:space:]]+", " ", paste(plot_data(as_ggplot(visible), "notes")$Text, collapse = "\n")), "observed means held fixed", fixed = TRUE)
  }
  expect_true(all(!d$CI_Eligible))
  expect_true(all(d$CI_ReportingUse[is.finite(d$CI_Lower)] == "diagnostic_only"))
  expect_true(any(grepl("These intervals omit some sources of estimation uncertainty", p$data$notes$Text)))
  b <- fair_average_table(fit, fair_se = TRUE, reference = "zero")
  expect_true(all(!b$stacked$FairCIEligible))
  expect_true(all(b$stacked$FairCIReportingUse == "unavailable"))
})
