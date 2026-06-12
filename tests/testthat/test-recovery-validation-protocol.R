recovery_validation_script_path <- function() {
  candidates <- c(
    file.path("inst", "validation", "recovery-validation.R"),
    file.path("..", "..", "inst", "validation", "recovery-validation.R"),
    system.file("validation", "recovery-validation.R", package = "mfrmr")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0L) return(NA_character_)
  candidates[1]
}

test_that("release recovery-validation script exposes review steps and case plan", {
  script <- recovery_validation_script_path()
  expect_true(file.exists(script))

  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)

  steps <- env$mfrmr_recovery_validation_prompt_steps()
  expect_s3_class(steps, "data.frame")
  expect_true(all(c("Step", "Label", "Prompt", "Output") %in% names(steps)))
  expect_true(any(grepl("threshold", steps$Prompt, ignore.case = TRUE)))
  expect_true(any(grepl("Markdown", steps$Output, ignore.case = TRUE)))

  plan <- env$mfrmr_recovery_validation_case_plan()
  expect_s3_class(plan, "data.frame")
  expect_true(all(c("CaseID", "Tier", "Model", "Purpose", "SummaryFocus") %in% names(plan)))
  expect_true(all(c("RSM", "PCM", "GPCM") %in% plan$Model))
  expect_true(any(plan$Tier == "extended"))
  expect_true("gpcm_high_dispersion_sparse" %in% plan$CaseID)
})

test_that("release recovery-validation GPCM case uses identified slope scale", {
  env <- new.env(parent = globalenv())
  sys.source(recovery_validation_script_path(), envir = env)

  spec <- env$mfrmr_recovery_validation_spec("gpcm_slope_profile")
  expect_s3_class(spec, "mfrm_sim_spec")
  expect_equal(spec$model, "GPCM")
  expect_true(is.data.frame(spec$slope_table))
  expect_equal(exp(mean(log(spec$slope_table$Estimate))), 1, tolerance = 1e-12)

  stress_spec <- env$mfrmr_recovery_validation_spec("gpcm_high_dispersion_sparse")
  expect_s3_class(stress_spec, "mfrm_sim_spec")
  expect_equal(stress_spec$model, "GPCM")
  expect_identical(stress_spec$slope_regime, "high_dispersion")
  expect_equal(stress_spec$score_levels, 5)
})

test_that("release recovery-validation summary handles completed case objects", {
  env <- new.env(parent = globalenv())
  sys.source(recovery_validation_script_path(), envir = env)

  plan <- env$mfrmr_recovery_validation_case_plan()[1, , drop = FALSE]
  assessment <- list(
    overview = data.frame(
      Reps = 2,
      SuccessfulRuns = 2,
      SuccessRate = 1,
      ConvergedRuns = 2,
      ConvergenceRate = 1,
      RecoveryRows = 6,
      RecoveryGroups = 1,
      OverallStatus = "ok"
    ),
    metric_review = data.frame(
      ParameterType = "facet",
      Facet = "Rater",
      ComparisonScale = "logit",
      RMSE = 0.2,
      Bias = 0.01,
      RMSEStatus = "ok",
      BiasStatus = "ok",
      CoverageStatus = "not_available",
      SEStatus = "not_assessed",
      MonteCarloStatus = "ok",
      OverallStatus = "ok",
      stringsAsFactors = FALSE
    ),
    condition_review = data.frame(
      Model = "GPCM",
      GPCMSlopeRegime = "high_dispersion",
      StressLevel = "high",
      SlopeLevels = 2L,
      MaxAbsCenteredLogSlope = 0.79,
      Replications = 2L,
      ScoreSupportReplications = 2L,
      MinScoreCount = 0L,
      MinScoreProportion = 0,
      MaxZeroScoreLevels = 1L,
      ScoreSupportStatus = "review",
      Status = "ok",
      Interpretation = "High-dispersion stress condition.",
      ScoreSupportInterpretation = "Sparse category stress condition.",
      ScoreSupportNextAction = "Report sparse support.",
      NextAction = "Report generator condition.",
      stringsAsFactors = FALSE
    ),
    condition_reporting_notes = data.frame(
      Model = c("GPCM", "GPCM"),
      GPCMSlopeRegime = c("high_dispersion", "high_dispersion"),
      StressLevel = c("high", "high"),
      ConditionArea = c("slope_regime", "score_support"),
      ReportingAttention = c("reporting_review", "reporting_review"),
      ConditionFinding = c("high_dispersion_slope_stress",
                           "omitted_generated_score_categories"),
      Evidence = c(
        "model=GPCM; slope_regime=high_dispersion; stress_level=high; slope_levels=2; max_abs_centered_log_slope=0.79",
        "score_support_replications=2; min_score_count=0; min_score_proportion=0; max_zero_score_levels=1"
      ),
      ReportingImplication = c(
        "Report high-dispersion as stress context, not an adequacy cut point.",
        "Report sparse generated score support as condition stress."
      ),
      NextAction = c(
        "Read slope recovery with uncertainty and score support.",
        "Inspect category-level recovery before generalizing."
      ),
      ValidationUse = "generator_condition_not_release_gate",
      stringsAsFactors = FALSE
    ),
    diagnostic_review = data.frame(
      Facet = "Rater",
      Replications = 2L,
      MeanSeparation = 2.1,
      MeanReliability = 0.82,
      MeanInfit = 1.00,
      MeanOutfit = 1.02,
      MeanMisfitRateAbsZ2 = 0,
      MeanDfSensitiveFlagRate = 0,
      DiagnosticAvailability = "available",
      Status = "not_assessed",
      ValidationUse = "diagnostic_only_not_release_gate",
      Interpretation = "Diagnostic context only.",
      NextAction = "Keep separate from release status.",
      stringsAsFactors = FALSE
    ),
    next_actions = "Keep this case in the validation bundle."
  )
  recovery <- list(notes = "No immediate warnings.")
  fake <- structure(
    list(
      prompt_steps = env$mfrmr_recovery_validation_prompt_steps(),
      plan = plan,
      results = list(list(
        case = plan,
        recovery = recovery,
        assessment = assessment,
        reps = 2,
        error = NA_character_
      )),
      started_at = as.POSIXct("2026-05-10 00:00:00", tz = "UTC"),
      completed_at = as.POSIXct("2026-05-10 00:01:00", tz = "UTC")
    ),
    class = "mfrmr_recovery_validation"
  )

  summary <- env$mfrmr_summarize_recovery_validation(fake)
  expect_s3_class(summary$reading_order, "data.frame")
  expect_true(all(c("Step", "Route", "WhatToRead", "Purpose") %in%
                    names(summary$reading_order)))
  expect_s3_class(summary$case_summary, "data.frame")
  expect_equal(summary$case_summary$OverallStatus, "ok")
  expect_equal(summary$case_summary$ReleaseRecoveryStatus, "ok")
  expect_equal(summary$case_summary$RecoveryMetricStatus, "ok")
  expect_equal(summary$case_summary$UncertaintyStatus, "review")
  expect_equal(summary$case_summary$MonteCarloStatus, "ok")
  expect_equal(summary$case_summary$DiagnosticStatus, "available")
  expect_equal(summary$case_summary$GPCMSlopeRegime, "high_dispersion")
  expect_equal(summary$case_summary$ScoreSupportStatus, "review")
  expect_equal(summary$case_summary$WorstMetric, "facet / Rater / logit")
  expect_s3_class(summary$topline_release_decision, "data.frame")
  expect_equal(summary$topline_release_decision$ReleaseRecoveryStatus, "ok")
  expect_equal(summary$topline_release_decision$UncertaintyStatus, "review")
  expect_s3_class(summary$release_decision_table, "data.frame")
  expect_true(all(c("ReleaseRecoveryStatus", "UncertaintyLimitation") %in%
                    names(summary$release_decision_table)))
  expect_s3_class(summary$metric_summary, "data.frame")
  expect_true(all(c("CaseID", "OverallStatus", "RMSE") %in% names(summary$metric_summary)))
  expect_s3_class(summary$condition_summary, "data.frame")
  expect_true(all(c("CaseID", "GPCMSlopeRegime", "ScoreSupportStatus") %in%
                    names(summary$condition_summary)))
  expect_s3_class(summary$condition_reporting_notes, "data.frame")
  expect_true(all(c("CaseID", "ConditionArea", "ReportingAttention",
                    "ConditionFinding", "ValidationUse") %in%
                    names(summary$condition_reporting_notes)))
  expect_equal(summary$condition_reporting_notes$ValidationUse,
               rep("generator_condition_not_release_gate", 2))
  expect_s3_class(summary$diagnostic_oc_summary, "data.frame")
  expect_true(all(c("CaseID", "Facet", "ValidationUse") %in%
                    names(summary$diagnostic_oc_summary)))
  expect_true(all(summary$diagnostic_oc_summary$ValidationUse ==
                    "diagnostic_only_not_release_gate"))
  expect_s3_class(summary$diagnostic_reporting_notes, "data.frame")
  expect_true(all(c("CaseID", "Facet", "ReportingAttention",
                    "DiagnosticFinding", "ValidationUse") %in%
                    names(summary$diagnostic_reporting_notes)))
  expect_equal(summary$diagnostic_reporting_notes$ReportingAttention, "context")
  expect_equal(summary$diagnostic_reporting_notes$DiagnosticFinding,
               "diagnostic_context_available")
  expect_s3_class(summary$domain_decision_table, "data.frame")
  expect_true(all(c("CaseID", "StatusDomain", "Status") %in% names(summary$domain_decision_table)))
  expect_true("score_support" %in% summary$domain_decision_table$StatusDomain)
  expect_true("diagnostic_operating_characteristics" %in%
                summary$domain_decision_table$StatusDomain)

  md <- env$mfrmr_recovery_validation_markdown(fake)
  expect_match(md, "Top-line release decision", fixed = TRUE)
  expect_match(md, "Recommended reading order", fixed = TRUE)
  expect_match(md, "Release decision by case", fixed = TRUE)
  expect_match(md, "Review steps", fixed = TRUE)
  expect_match(md, "Case summary", fixed = TRUE)
  expect_match(md, "Domain decision table", fixed = TRUE)
  expect_match(md, "Condition summary", fixed = TRUE)
  expect_match(md, "Condition reporting notes", fixed = TRUE)
  expect_match(md, "Diagnostic reporting notes", fixed = TRUE)
  expect_match(md, "Diagnostic operating-characteristic summary", fixed = TRUE)

  validation_summary <- env$summary.mfrmr_recovery_validation(fake)
  expect_s3_class(validation_summary, "summary.mfrmr_recovery_validation")
  expect_s3_class(validation_summary$reading_order, "data.frame")
  expect_s3_class(validation_summary$condition_summary, "data.frame")
  expect_s3_class(validation_summary$condition_reporting_notes, "data.frame")
  expect_s3_class(validation_summary$diagnostic_reporting_notes, "data.frame")
  expect_s3_class(validation_summary$diagnostic_oc_summary, "data.frame")
  expect_s3_class(validation_summary$topline_release_decision, "data.frame")
  expect_output(env$print.summary.mfrmr_recovery_validation(validation_summary), "mfrmr Recovery Validation Summary")
  expect_output(env$print.summary.mfrmr_recovery_validation(validation_summary), "Release recovery status: ok")
  expect_output(env$print.summary.mfrmr_recovery_validation(validation_summary), "Recommended reading order")
  expect_output(env$print.summary.mfrmr_recovery_validation(validation_summary), "Case decisions")
  expect_output(env$print.mfrmr_recovery_validation(fake), "mfrmr Recovery Validation Summary")

  out_dir <- tempfile("mfrmr-validation-output-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  written <- env$mfrmr_write_recovery_validation_outputs(
    fake,
    output_dir = out_dir,
    prefix = "fake_validation"
  )
  expect_true(dir.exists(written))
  expected_files <- file.path(out_dir, paste0("fake_validation_", c(
    "topline_release_decision.csv",
    "reading_order.csv",
    "release_decision_table.csv",
    "prompt_steps.csv",
    "case_plan.csv",
    "case_summary.csv",
    "condition_summary.csv",
    "condition_reporting_notes.csv",
    "diagnostic_reporting_notes.csv",
    "diagnostic_oc_summary.csv",
    "metric_summary.csv",
    "decision_table.csv",
    "domain_decision_table.csv",
    "run_notes.csv"
  )))
  expect_true(all(file.exists(expected_files)))
  expect_true(file.exists(file.path(out_dir, "fake_validation.md")))
  expect_true(file.exists(file.path(out_dir, "fake_validation.rds")))

  exported_condition_notes <- utils::read.csv(
    file.path(out_dir, "fake_validation_condition_reporting_notes.csv"),
    stringsAsFactors = FALSE
  )
  expect_true(all(c("ConditionArea", "ReportingAttention", "ConditionFinding",
                    "ValidationUse") %in% names(exported_condition_notes)))
  expect_true("omitted_generated_score_categories" %in%
                exported_condition_notes$ConditionFinding)

  exported_notes <- utils::read.csv(
    file.path(out_dir, "fake_validation_diagnostic_reporting_notes.csv"),
    stringsAsFactors = FALSE
  )
  expect_true(all(c("ReportingAttention", "DiagnosticFinding", "ValidationUse") %in%
                    names(exported_notes)))
  expect_equal(exported_notes$DiagnosticFinding, "diagnostic_context_available")
  exported_markdown <- paste(readLines(file.path(out_dir, "fake_validation.md"), warn = FALSE),
                             collapse = "\n")
  expect_match(exported_markdown, "## Diagnostic reporting notes", fixed = TRUE)
  expect_match(exported_markdown, "## Condition reporting notes", fixed = TRUE)
  expect_match(exported_markdown, "diagnostic_context_available", fixed = TRUE)
})

test_that("release recovery-validation top line separates core and extended cases", {
  env <- new.env(parent = globalenv())
  sys.source(recovery_validation_script_path(), envir = env)

  case_summary <- data.frame(
    CaseID = c("gpcm_slope_profile", "gpcm_high_dispersion_sparse"),
    Tier = c("core", "extended"),
    ReleaseRecoveryStatus = c("review", "concern"),
    RecoveryMetricStatus = c("ok", "concern"),
    MonteCarloStatus = c("not_available", "not_available"),
    UncertaintyStatus = c("review", "review"),
    DiagnosticStatus = c("available", "review"),
    stringsAsFactors = FALSE
  )
  decision <- env$mfrmr_recovery_validation_release_decision(case_summary)

  expect_equal(decision$Cases, 2)
  expect_equal(decision$CoreCases, 1)
  expect_equal(decision$ExtendedCases, 1)
  expect_equal(decision$ReleaseRecoveryStatus, "review")
  expect_equal(decision$IncludedCaseStatus, "concern")
  expect_equal(decision$ExtendedSensitivityStatus, "concern")
  expect_equal(decision$DiagnosticStatus, "review")
  expect_match(decision$Conclusion, "Extended sensitivity cases", fixed = TRUE)
})
