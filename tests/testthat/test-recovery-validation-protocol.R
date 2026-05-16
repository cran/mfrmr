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
})

test_that("release recovery-validation GPCM case uses identified slope scale", {
  env <- new.env(parent = globalenv())
  sys.source(recovery_validation_script_path(), envir = env)

  spec <- env$mfrmr_recovery_validation_spec("gpcm_slope_profile")
  expect_s3_class(spec, "mfrm_sim_spec")
  expect_equal(spec$model, "GPCM")
  expect_true(is.data.frame(spec$slope_table))
  expect_equal(exp(mean(log(spec$slope_table$Estimate))), 1, tolerance = 1e-12)
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
  expect_s3_class(summary$case_summary, "data.frame")
  expect_equal(summary$case_summary$OverallStatus, "ok")
  expect_equal(summary$case_summary$ReleaseRecoveryStatus, "ok")
  expect_equal(summary$case_summary$RecoveryMetricStatus, "ok")
  expect_equal(summary$case_summary$UncertaintyStatus, "review")
  expect_equal(summary$case_summary$MonteCarloStatus, "ok")
  expect_equal(summary$case_summary$WorstMetric, "facet / Rater / logit")
  expect_s3_class(summary$topline_release_decision, "data.frame")
  expect_equal(summary$topline_release_decision$ReleaseRecoveryStatus, "ok")
  expect_equal(summary$topline_release_decision$UncertaintyStatus, "review")
  expect_s3_class(summary$release_decision_table, "data.frame")
  expect_true(all(c("ReleaseRecoveryStatus", "UncertaintyLimitation") %in%
                    names(summary$release_decision_table)))
  expect_s3_class(summary$metric_summary, "data.frame")
  expect_true(all(c("CaseID", "OverallStatus", "RMSE") %in% names(summary$metric_summary)))
  expect_s3_class(summary$domain_decision_table, "data.frame")
  expect_true(all(c("CaseID", "StatusDomain", "Status") %in% names(summary$domain_decision_table)))

  md <- env$mfrmr_recovery_validation_markdown(fake)
  expect_match(md, "Top-line release decision", fixed = TRUE)
  expect_match(md, "Release decision by case", fixed = TRUE)
  expect_match(md, "Review steps", fixed = TRUE)
  expect_match(md, "Case summary", fixed = TRUE)
  expect_match(md, "Domain decision table", fixed = TRUE)

  validation_summary <- env$summary.mfrmr_recovery_validation(fake)
  expect_s3_class(validation_summary, "summary.mfrmr_recovery_validation")
  expect_s3_class(validation_summary$topline_release_decision, "data.frame")
  expect_output(env$print.summary.mfrmr_recovery_validation(validation_summary), "mfrmr Recovery Validation Summary")
  expect_output(env$print.summary.mfrmr_recovery_validation(validation_summary), "Release recovery status: ok")
  expect_output(env$print.summary.mfrmr_recovery_validation(validation_summary), "Case decisions")
  expect_output(env$print.mfrmr_recovery_validation(fake), "mfrmr Recovery Validation Summary")
})
