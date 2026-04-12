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

test_that("evaluate_mfrm_diagnostic_screening rejects unsupported GPCM paths", {
  expect_error(
    mfrmr::evaluate_mfrm_diagnostic_screening(
      design = list(person = 12, rater = 2, criterion = 2, assignment = 2),
      reps = 1,
      model = "GPCM",
      seed = 1
    ),
    "does not yet support bounded `GPCM`",
    fixed = TRUE
  )
})
