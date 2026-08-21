jml_profile_recovery_pilot_path <- function() {
  testthat::test_path(
    "..", "..", "inst", "validation",
    "jml-extreme-profile-recovery-pilot-0.2.3.R"
  )
}

load_jml_profile_recovery_pilot <- function() {
  runner <- jml_profile_recovery_pilot_path()
  skip_if_not(file.exists(runner),
              "repository-internal validation artifacts are excluded")
  env <- new.env(parent = globalenv())
  sys.source(runner, envir = env)
  env
}

test_that("profile-recovery manifests freeze the calibration design", {
  env <- load_jml_profile_recovery_pilot()
  smoke <- env$mfrmr_jml_profile_recovery_manifest("smoke")
  pilot <- env$mfrmr_run_jml_profile_recovery_pilot(
    "pilot", dry_run = TRUE
  )

  expect_equal(nrow(smoke), 3L)
  expect_equal(nrow(pilot), 90L)
  expect_setequal(smoke$Model, c("RSM", "PCM", "GPCM"))
  expect_setequal(pilot$Information, c("low", "high"))
  expect_setequal(pilot$ExtremeFraction, c(0, 0.10, 0.25))
  expect_equal(sort(unique(pilot$Replicate)), seq_len(5L))
  expect_false(anyDuplicated(pilot$Seed) > 0L)
  expect_true(all(pilot$ForcedExtremeN %% 2L == 0L))
  expect_true(all(!pilot$PersonRecoveryEligible))
  expect_true(all(pilot$StructuralRecoveryEligible))
  expect_true(all(pilot$ExtremeAdjustment == "none"))
  expect_true(all(pilot$BiasCorrection == "none"))
  expect_true(all(pilot$LimitCaps == "4,8,12,16,24,32,48,64"))
  expect_true(all(pilot$LimitTolerance == 1e-8))
})

test_that("the full profile-recovery pilot requires explicit authorization", {
  env <- load_jml_profile_recovery_pilot()
  expect_error(
    env$mfrmr_run_jml_profile_recovery_pilot("pilot"),
    "authorize_pilot = TRUE"
  )
})

test_that("paired profile-recovery smoke closes all three model families", {
  env <- load_jml_profile_recovery_pilot()
  result <- env$mfrmr_run_jml_profile_recovery_pilot(
    "smoke", progress = FALSE
  )

  expect_true(result$ContractPassed)
  expect_false(result$EvidenceReady)
  expect_identical(result$ReadinessEffect, "none_pilot_only")
  expect_equal(nrow(result$Runs), 3L)
  expect_true(all(result$Runs$RawFitOK))
  expect_true(all(result$Runs$ForcedTyped))
  expect_true(all(result$Runs$ForcedDirectionCorrect))
  expect_true(all(result$Runs$ProfileOK))
  expect_true(all(result$Runs$PairedRowsComplete))
  expect_true(all(!result$Runs$FalseReady))
  expect_true(all(result$Runs$ProfileState ==
                    "profile_limit_refit_verified"))
  expect_true(all(result$Runs$TerminalLimitGap <= 1e-8))
  expect_true(all(result$Runs$ForcedTypedN ==
                    result$Runs$ForcedExtremeN))
  expect_true(all(result$Runs$ActualFreeExtremeN >=
                    result$Runs$ForcedExtremeN))

  recovery <- result$Recovery
  expect_equal(nrow(recovery), 64L)
  expect_setequal(
    recovery$EstimatorVariant,
    c("raw_finite_jml", "extended_profile_limit")
  )
  expect_false(any(recovery$ParameterType == "person"))
  expect_true(all(recovery$PersonRecoveryEligible == FALSE))
  expect_true(all(recovery$RawConvergenceSeverity %in% c("pass", "review")))
  expect_true(all(recovery$RawFitReadiness %in%
                    c("ready", "ready_with_exclusions", "review")))
  expect_true(all(recovery$ProfileOK))
  expect_true(all(recovery$ActualFreeExtremeN >=
                    recovery$ForcedExtremeN))
  expect_true(all(is.na(
    recovery$SE[recovery$EstimatorVariant == "extended_profile_limit"]
  )))
  expect_true(all(is.na(
    recovery$Covered95[
      recovery$EstimatorVariant == "extended_profile_limit"
    ]
  )))

  gpcm_slopes <- recovery[
    recovery$Model == "GPCM" & recovery$ParameterType == "slope" &
      recovery$EstimatorVariant == "raw_finite_jml",
    , drop = FALSE
  ]
  expect_equal(nrow(gpcm_slopes), 2L)
  expect_true(all(is.finite(gpcm_slopes$Truth)))
  expect_gt(diff(range(gpcm_slopes$Truth)), 0.5)
})

test_that("profile recovery never inherits raw structural uncertainty", {
  env <- load_jml_profile_recovery_pilot()
  row <- env$mfrmr_jml_profile_recovery_manifest("smoke")
  row <- row[row$Model == "PCM", , drop = FALSE]
  result <- env$mfrmr_jml_profile_recovery_fit_one(row)

  profile <- result$recovery[
    result$recovery$EstimatorVariant == "extended_profile_limit",
    , drop = FALSE
  ]
  expect_gt(nrow(profile), 0L)
  expect_true(all(is.na(profile$SE)))
  expect_true(all(is.na(profile$Covered95)))
  expect_false(result$profile$FiniteOriginalJMLMaximum)
  expect_false(result$profile$OriginalLikelihoodMaximumAttained)
})
