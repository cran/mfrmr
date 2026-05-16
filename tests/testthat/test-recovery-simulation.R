test_that("evaluate_mfrm_recovery returns row-level and summary recovery tables", {
  rec <- suppressWarnings(
    evaluate_mfrm_recovery(
      n_person = 12,
      n_rater = 2,
      n_criterion = 2,
      raters_per_person = 2,
      reps = 1,
      maxit = 8,
      seed = 20260509
    )
  )

  expect_s3_class(rec, "mfrm_recovery_simulation")
  expect_s3_class(rec$recovery, "data.frame")
  expect_s3_class(rec$recovery_summary, "data.frame")
  expect_s3_class(rec$rep_overview, "data.frame")
  expect_true(all(c("Truth", "Estimate", "EstimateAligned", "ErrorAligned",
                    "RawTruth", "RawEstimate", "ComparisonScale") %in%
                    names(rec$recovery)))
  expect_true(all(c("ParameterType", "Facet", "RMSE", "Bias", "MAE",
                    "Correlation", "Coverage95", "McseBias", "McseRMSE") %in%
                    names(rec$recovery_summary)))
  expect_true(any(rec$recovery$ParameterType == "facet"))
  expect_true(any(rec$recovery$ParameterType == "step"))
  expect_true(all(is.finite(rec$recovery_summary$RMSE)))
  expect_true(is.list(rec$ademp))
  expect_match(rec$ademp$aims, "parameter recovery", ignore.case = TRUE)

  s <- summary(rec)
  expect_s3_class(s, "summary.mfrm_recovery_simulation")
  expect_true(all(c("overview", "recovery_summary", "rep_overview", "ademp") %in%
                    names(s)))
  expect_output(print(s), "MFRM Parameter Recovery Simulation Summary")
  expect_output(print(rec), "MFRM Parameter Recovery Simulation Summary")

  summary_plot <- plot(rec, type = "summary", metric = "rmse", draw = FALSE)
  expect_s3_class(summary_plot, "mfrm_plot_data")
  expect_identical(summary_plot$name, "recovery_simulation")
  expect_identical(summary_plot$data$type, "summary")
  expect_true(all(c("plot_table", "metric", "metric_label", "notes") %in%
                    names(summary_plot$data)))
  expect_true(nrow(summary_plot$data$plot_table) > 0)

  error_plot <- plot(rec, type = "errors", parameter_type = "facet", draw = FALSE)
  expect_s3_class(error_plot, "mfrm_plot_data")
  expect_identical(error_plot$data$type, "errors")
  expect_true(all(error_plot$data$plot_table$ParameterType == "facet"))

  scatter_plot <- plot(rec, type = "scatter", comparison = "unaligned", draw = FALSE)
  expect_s3_class(scatter_plot, "mfrm_plot_data")
  expect_identical(scatter_plot$data$comparison, "unaligned")

  status_plot <- plot(rec, type = "replications", draw = FALSE)
  expect_s3_class(status_plot, "mfrm_plot_data")
  expect_identical(status_plot$data$type, "replications")
  expect_true("rep_overview" %in% names(status_plot$data))

  assessment <- assess_mfrm_recovery(
    rec,
    min_reps = 1,
    max_rmse = c(default = 2),
    max_abs_bias = c(default = 1)
  )
  expect_s3_class(assessment, "mfrm_recovery_assessment")
  expect_s3_class(assessment$overview, "data.frame")
  expect_s3_class(assessment$checklist, "data.frame")
  expect_s3_class(assessment$metric_review, "data.frame")
  expect_true(all(c("Section", "Item", "Status", "Evidence", "NextAction") %in%
                    names(assessment$checklist)))
  expect_true(all(c("RMSEStatus", "BiasStatus", "CoverageStatus",
                    "OverallStatus", "NextAction") %in%
                    names(assessment$metric_review)))
  expect_true(length(assessment$next_actions) > 0)

  assessment_summary <- summary(assessment)
  expect_s3_class(assessment_summary, "summary.mfrm_recovery_assessment")
  expect_output(print(assessment_summary), "MFRM Recovery Adequacy Assessment")
  expect_output(print(assessment), "MFRM Recovery Adequacy Assessment")

  assessment_status_plot <- plot(assessment, type = "status", draw = FALSE)
  expect_s3_class(assessment_status_plot, "mfrm_plot_data")
  expect_identical(assessment_status_plot$name, "recovery_assessment")
  expect_identical(assessment_status_plot$data$type, "status")
  expect_true(all(c("plot_table", "section_status", "status_counts",
                    "checklist", "reading_order", "guidance") %in%
                    names(assessment_status_plot$data)))
  expect_true(nrow(assessment_status_plot$data$status_counts) > 0)
  expect_true("AttentionOrder" %in% names(assessment_status_plot$data$section_status))
  expect_true(any(grepl("Read this plot before metric-level plots",
                        assessment_status_plot$data$guidance, fixed = TRUE)))

  assessment_metric_plot <- plot(assessment, type = "metrics",
                                 metric = "rmse", draw = FALSE)
  expect_s3_class(assessment_metric_plot, "mfrm_plot_data")
  expect_identical(assessment_metric_plot$name, "recovery_assessment")
  expect_identical(assessment_metric_plot$data$type, "metrics")
  expect_identical(assessment_metric_plot$data$metric, "rmse")
  expect_true(all(c("Value", "Limit", "Status", "AttentionOrder") %in%
                    names(assessment_metric_plot$data$plot_table)))
  expect_true(all(c("metric_review", "reading_order", "guidance") %in%
                    names(assessment_metric_plot$data)))
  expect_true(nrow(assessment_metric_plot$data$plot_table) > 0)
  expect_true(any(grepl("sorted by status priority",
                        assessment_metric_plot$data$guidance, fixed = TRUE)))

  assessment_no_limits <- assess_mfrm_recovery(
    rec,
    min_reps = 1,
    min_se_available = NULL,
    max_mcse_rmse_ratio = NULL
  )
  rmse_row <- assessment_no_limits$checklist[
    assessment_no_limits$checklist$Item == "RMSE threshold",
    ,
    drop = FALSE
  ]
  expect_identical(rmse_row$Status[1], "not_assessed")
})

test_that("evaluate_mfrm_recovery refits on the declared generator score support", {
  sim <- simulate_mfrm_data(
    n_person = 10,
    n_rater = 2,
    n_criterion = 2,
    raters_per_person = 2,
    score_levels = 5,
    seed = 20260510
  )
  fit_args <- mfrmr:::simulation_add_fit_score_support(list(data = sim), sim)
  expect_identical(fit_args$rating_min, 1L)
  expect_identical(fit_args$rating_max, 5L)

  messages <- character()
  rec <- withCallingHandlers(
    suppressWarnings(
      evaluate_mfrm_recovery(
        n_person = 10,
        n_rater = 2,
        n_criterion = 2,
        raters_per_person = 2,
        score_levels = 5,
        reps = 1,
        maxit = 5,
        seed = 20260511
      )
    ),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_s3_class(rec, "mfrm_recovery_simulation")
  expect_false(any(grepl("Rating range inferred", messages, fixed = TRUE)))
})

test_that("evaluate_mfrm_recovery supports fitted bounded GPCM slope rows", {
  spec <- build_mfrm_sim_spec(
    n_person = 14,
    n_rater = 2,
    n_criterion = 2,
    raters_per_person = 2,
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    slopes = c(0.85, 1.15),
    assignment = "crossed"
  )

  rec <- suppressWarnings(
    evaluate_mfrm_recovery(
      sim_spec = spec,
      reps = 1,
      fit_method = "MML",
      quad_points = 5,
      maxit = 12,
      seed = 20260510,
      include_person = FALSE
    )
  )

  expect_s3_class(rec, "mfrm_recovery_simulation")
  expect_true(any(rec$recovery$ParameterType == "slope"))
  slope_rows <- rec$recovery[rec$recovery$ParameterType == "slope", , drop = FALSE]
  expect_true(all(slope_rows$ComparisonScale == "log_slope"))
  expect_true(all(!slope_rows$AlignWithinGroup))
  expect_equal(slope_rows$AlignmentShift, rep(0, nrow(slope_rows)), tolerance = 1e-12)
  expect_true(all(slope_rows$RecoveryBasis == "geometric_mean_one_log_slope"))
  expect_true(all(is.finite(slope_rows$RawTruth)))
  expect_true(all(is.finite(slope_rows$RawEstimate)))
  expect_equal(exp(mean(log(slope_rows$RawTruth))), 1, tolerance = 1e-12)
  expect_true(any(rec$recovery_summary$ParameterType == "slope"))
  expect_true(any(grepl("Bounded GPCM recovery", rec$notes, fixed = TRUE)))

  slope_plot <- plot(rec, type = "summary", metric = "bias",
                     parameter_type = "slope", draw = FALSE)
  expect_s3_class(slope_plot, "mfrm_plot_data")
  expect_true(all(slope_plot$data$plot_table$ParameterType == "slope"))

  assessment <- assess_mfrm_recovery(
    rec,
    min_reps = 1,
    max_rmse = c(slope = 2),
    max_abs_bias = c(slope = 1),
    min_se_available = NULL,
    max_mcse_rmse_ratio = NULL
  )
  expect_s3_class(assessment, "mfrm_recovery_assessment")
  slope_review <- assessment$metric_review[
    assessment$metric_review$ParameterType == "slope",
    ,
    drop = FALSE
  ]
  expect_true(nrow(slope_review) > 0)
  expect_true(all(is.finite(slope_review$RMSELimit)))
})
