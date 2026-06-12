test_that("evaluate_mfrm_diagnostic_screening separates well-specified and local-dependence scenarios", {
  for (model in c("RSM", "PCM")) {
    diag_eval <- suppressWarnings(mfrmr::evaluate_mfrm_diagnostic_screening(
      design = list(person = 18, rater = 3, criterion = 3, assignment = 3),
      reps = 2,
      model = model,
      maxit = 8,
      quad_points = 7,
      local_dependence_sd = 1.0,
      seed = 20260409
    ))

    expect_s3_class(diag_eval, "mfrm_diagnostic_screening")
    expect_true(is.data.frame(diag_eval$results))
    expect_true(is.data.frame(diag_eval$scenario_summary))
    expect_true(is.data.frame(diag_eval$performance_summary))
    expect_true(is.data.frame(diag_eval$scenario_contrast))
    expect_true(all(c("well_specified", "local_dependence") %in% diag_eval$scenario_summary$Scenario))
    expect_true(all(diag_eval$results$RunOK))
    expect_true(all(diag_eval$results$MarginalAvailable))
    expect_true(nrow(diag_eval$scenario_contrast) >= 1)
    expect_true(all(diag_eval$scenario_contrast$Scenario == "local_dependence"))
    expect_true("StrictSignalImproved" %in% names(diag_eval$scenario_contrast))
    expect_true("StrictSignalDominatesLegacy" %in% names(diag_eval$scenario_contrast))
    expect_true(all(diag_eval$scenario_contrast$DeltaLegacyMeanAbsZ >= 0))
    expect_true(all(c("type_I_proxy", "sensitivity_proxy") %in% diag_eval$performance_summary$EvaluationUse))
    expect_true(all(diag_eval$performance_summary$MeanElapsedSec >= 0))
    expect_true(all(diag_eval$scenario_summary$PairwiseAvailabilityRate > 0))
    expect_true(length(diag_eval$notes) >= 1)

    if (identical(model, "RSM")) {
      overview_plot <- plot(diag_eval, type = "overview", draw = FALSE)
      expect_s3_class(overview_plot, "mfrm_plot_data")
      expect_identical(overview_plot$name, "diagnostic_screening")
      expect_identical(overview_plot$data$type, "overview")
      expect_true(all(c("Signal", "Value", "SourceTable", "group") %in% names(overview_plot$data$plot_long)))
      expect_true(any(overview_plot$data$plot_long$Signal == "Strict combined any-flag rate"))
      expect_true(any(mfrmr::plot_data(diag_eval, type = "overview", component = "plot_long")$Signal == "Strict combined any-flag rate"))
      expect_true(all(c("overview", "reading_order", "next_actions", "reporting_notes", "figure_recipes") %in% names(overview_plot$data)))
      expect_true(any(mfrmr::plot_data(overview_plot, component = "next_actions")$Area == "Appendix and plot-data handoff"))
      expect_true(any(mfrmr::plot_data(overview_plot, component = "figure_recipes")$FigureID == "overview_rates"))
      expect_true(any(mfrmr::plot_data_components(overview_plot)$Component == "next_actions"))

      diag_summary <- summary(diag_eval)
      expect_s3_class(diag_summary, "summary.mfrm_diagnostic_screening")
      expect_true(all(c(
        "overview", "reading_order", "next_actions", "reporting_notes", "figure_recipes",
        "scenario_summary", "performance_summary",
        "plot_overview_rate", "plot_runtime"
      ) %in% names(diag_summary)))
      expect_true(any(diag_summary$reading_order$Table == "scenario_summary"))
      expect_true(any(diag_summary$reading_order$Table == "figure_recipes"))
      expect_true(any(diag_summary$next_actions$Area == "Appendix and plot-data handoff"))
      expect_true(any(diag_summary$figure_recipes$SummaryTable == "plot_overview_rate"))
      expect_true(any(grepl("export_summary_appendix", diag_summary$next_actions$Route, fixed = TRUE)))
      expect_true(any(grepl("release gates", diag_summary$reporting_notes$ReportingBoundary, fixed = TRUE)))
      expect_true(any(diag_summary$plot_overview_rate$Signal == "Strict combined any-flag rate"))
      expect_true(any(grepl("operating-characteristic", diag_summary$notes, fixed = TRUE)))

      contrast_plot <- plot(diag_eval, type = "contrast", metric = "count", draw = FALSE)
      expect_s3_class(contrast_plot, "mfrm_plot_data")
      expect_identical(contrast_plot$data$type, "contrast")
      expect_true(any(contrast_plot$data$plot_long$Metric == "DeltaPairwiseFlaggedLevelPairs"))
      expect_true(all(contrast_plot$data$plot_long$Scenario == "local_dependence"))
    }
  }
})

test_that("evaluate_mfrm_diagnostic_screening supports latent misspecification scenarios", {
  for (model in c("RSM", "PCM")) {
    diag_eval <- suppressWarnings(mfrmr::evaluate_mfrm_diagnostic_screening(
      design = list(person = 18, rater = 3, criterion = 3, assignment = 3),
      reps = 1,
      scenarios = c("well_specified", "latent_misspecification"),
      model = model,
      maxit = 8,
      quad_points = 7,
      seed = 20260409
    ))

    expect_s3_class(diag_eval, "mfrm_diagnostic_screening")
    expect_true(all(c("well_specified", "latent_misspecification") %in% diag_eval$scenario_summary$Scenario))
    expect_true(all(diag_eval$results$RunOK))
    expect_true(all(diag_eval$results$MarginalAvailable))
    expect_true(is.data.frame(diag_eval$performance_summary))
    expect_true(nrow(diag_eval$scenario_contrast) >= 1)
    expect_true(all(diag_eval$scenario_contrast$Scenario == "latent_misspecification"))
    expect_true(any(diag_eval$performance_summary$EvaluationUse == "type_I_proxy"))
    expect_true(any(diag_eval$performance_summary$EvaluationUse == "sensitivity_proxy"))
    expect_true(all(diag_eval$scenario_summary$PairwiseAvailabilityRate > 0))
    expect_true(length(diag_eval$notes) >= 1)
  }
})

test_that("evaluate_mfrm_diagnostic_screening supports step-structure misspecification scenarios", {
  for (model in c("RSM", "PCM")) {
    diag_eval <- suppressWarnings(mfrmr::evaluate_mfrm_diagnostic_screening(
      design = list(person = 30, rater = 3, criterion = 3, assignment = 3),
      reps = 1,
      scenarios = c("well_specified", "step_structure_misspecification"),
      model = model,
      maxit = 10,
      quad_points = 9,
      seed = 20260409
    ))

    expect_s3_class(diag_eval, "mfrm_diagnostic_screening")
    expect_true(all(c("well_specified", "step_structure_misspecification") %in% diag_eval$scenario_summary$Scenario))
    expect_true(all(diag_eval$results$RunOK))
    expect_true(all(diag_eval$results$MarginalAvailable))
    expect_true(is.data.frame(diag_eval$performance_summary))
    expect_true(nrow(diag_eval$scenario_contrast) >= 1)
    expect_true(all(diag_eval$scenario_contrast$Scenario == "step_structure_misspecification"))
    expect_true("StrictSignalImproved" %in% names(diag_eval$scenario_contrast))
    expect_true("StrictSignalDominatesLegacy" %in% names(diag_eval$scenario_contrast))
    expect_true(all(diag_eval$scenario_contrast$DeltaLegacyMeanAbsZ >= 0))
    expect_true(any(diag_eval$performance_summary$StrictSensitivityProxy >= diag_eval$performance_summary$LegacySensitivityProxy, na.rm = TRUE))
    expect_true(all(diag_eval$scenario_summary$PairwiseAvailabilityRate > 0))
    expect_true(length(diag_eval$notes) >= 1)
  }
})

test_that("evaluate_mfrm_diagnostic_screening can retain report-index review signals", {
  diag_eval <- suppressWarnings(mfrmr::evaluate_mfrm_diagnostic_screening(
    design = list(person = 14, rater = 3, criterion = 3, assignment = 3),
    reps = 1,
    scenarios = "well_specified",
    model = "RSM",
    maxit = 8,
    quad_points = 7,
    include_report = TRUE,
    report_include = c("fit", "diagnostics", "tables", "precision", "reporting"),
    seed = 20260409
  ))

  expect_s3_class(diag_eval, "mfrm_diagnostic_screening")
  expect_true(all(c(
    "ReportIndexAvailable", "ReportFitReadiness", "ReportFitReviewSignalCount",
    "ReportPrecisionReadiness", "ReportPrecisionReviewSignalCount"
  ) %in% names(diag_eval$results)))
  expect_true(all(diag_eval$results$ReportIndexAvailable))
  expect_true(is.data.frame(diag_eval$report_signal_summary))
  expect_true(all(c(
    "ReportIndexAvailabilityRate", "MeanReportReviewAreas",
    "FitReportReviewRate", "MeanFitReportSignals"
  ) %in% names(diag_eval$report_signal_summary)))
  expect_true(all(diag_eval$report_signal_summary$ReportIndexAvailabilityRate > 0))

  report_plot <- plot(diag_eval, type = "report", metric = "rate", draw = FALSE)
  expect_s3_class(report_plot, "mfrm_plot_data")
  expect_identical(report_plot$data$type, "report")
  expect_true(any(report_plot$data$plot_long$Signal == "Report-index availability rate"))
  expect_true(any(report_plot$data$plot_long$Signal == "Fit review rate"))
  expect_true(any(mfrmr::plot_data(report_plot, component = "next_actions")$Area == "Report-index signals" &
                    mfrmr::plot_data(report_plot, component = "next_actions")$Status == "ok"))

  count_plot <- plot(diag_eval, type = "overview", metric = "count", draw = FALSE)
  expect_s3_class(count_plot, "mfrm_plot_data")
  expect_true(any(count_plot$data$plot_long$Signal == "Mean report fit signals"))
  expect_true(any(mfrmr::plot_data_components(count_plot)$Component == "plot_long"))
  expect_true(any(mfrmr::plot_data_components(count_plot)$Component == "reporting_notes"))
  expect_true(any(mfrmr::plot_data_components(count_plot)$Component == "figure_recipes"))

  diag_summary <- summary(diag_eval)
  expect_s3_class(diag_summary, "summary.mfrm_diagnostic_screening")
  expect_true(nrow(diag_summary$report_signal_summary) > 0L)
  expect_true(any(diag_summary$reading_order$Table == "report_signal_summary"))
  expect_true(any(diag_summary$figure_recipes$FigureID == "report_review_rates" &
                    diag_summary$figure_recipes$Availability == "available_when_report_rows_exist"))
  expect_true(any(diag_summary$next_actions$Area == "Report-index signals" &
                    diag_summary$next_actions$Status == "ok"))
  expect_true(any(grepl("report-index availability", diag_summary$reporting_notes$Evidence, fixed = TRUE)))
  expect_true(any(diag_summary$plot_report_rate$Signal == "Report-index availability rate"))
})

test_that("evaluate_mfrm_diagnostic_screening supports caveated GPCM screening output", {
  diag_eval <- suppressWarnings(mfrmr::evaluate_mfrm_diagnostic_screening(
    design = list(person = 10, rater = 2, criterion = 2, assignment = 2),
    reps = 1,
    scenarios = c("well_specified", "local_dependence"),
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    slopes = c(C01 = 1, C02 = 1),
    maxit = 10,
    seed = 1
  ))

  expect_s3_class(diag_eval, "mfrm_diagnostic_screening")
  expect_identical(diag_eval$settings$model, "GPCM")
  expect_identical(diag_eval$settings$gpcm_screening_status, "supported_with_caveat")
  expect_true(nrow(diag_eval$gpcm_boundary) > 0L)
  expect_true(all(diag_eval$results$Model == "GPCM"))
  expect_true(any(grepl("slope-aware operating-characteristic", diag_eval$notes, fixed = TRUE)))

  diag_summary <- summary(diag_eval)
  expect_s3_class(diag_summary, "summary.mfrm_diagnostic_screening")
  expect_true(nrow(diag_summary$gpcm_boundary) > 0L)
  expect_true(any(grepl("slope-aware operating-characteristic", diag_summary$notes, fixed = TRUE)))
  printed <- capture.output(print(diag_summary))
  expect_true(any(grepl("Bounded GPCM boundary", printed, fixed = TRUE)))

  diag_plot <- plot(diag_eval, type = "overview", metric = "rate", draw = FALSE)
  expect_s3_class(diag_plot, "mfrm_plot_data")
  expect_true(nrow(diag_plot$data$gpcm_boundary) > 0L)
})
