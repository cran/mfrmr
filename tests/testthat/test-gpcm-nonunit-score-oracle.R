gpcm_nonunit_oracle_validation_dir <- function() {
  candidates <- c(
    file.path("inst", "validation"),
    testthat::test_path("..", "..", "inst", "validation")
  )
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0L) return(NA_character_)
  candidates[1]
}

load_gpcm_nonunit_score_oracle <- function() {
  validation_dir <- gpcm_nonunit_oracle_validation_dir()
  testthat::skip_if(
    is.na(validation_dir),
    "Repository-only GPCM non-unit score oracle is unavailable."
  )
  env <- new.env(parent = globalenv())
  script <- file.path(
    validation_dir,
    "gpcm-nonunit-score-oracle-0.2.3.R"
  )
  sys.source(script, envir = env)
  list(env = env, script = script)
}

test_that("the non-unit GPCM oracle contract is bounded and never auto-runs", {
  loaded <- load_gpcm_nonunit_score_oracle()
  env <- loaded$env

  expect_identical(
    env$mfrmr_gno_contract_version,
    "mfrmr_gpcm_nonunit_score_oracle_v1"
  )
  expect_identical(env$mfrmr_gno_primary_step, 1e-5)
  expect_identical(
    env$mfrmr_gno_expected_points,
    c(
      "retained_solution",
      "high_dispersion_probe",
      "finite_slope_stress_forward",
      "finite_slope_stress_reverse"
    )
  )
  expect_identical(
    unname(env$mfrmr_gno_limits),
    c(1e-12, 1e-12, 1e-10, 1e-6, 1e-12, 1e-12)
  )
  expect_false(exists("fit", envir = env, inherits = FALSE))
  expect_false(exists("audit", envir = env, inherits = FALSE))
  expect_false(exists("decision", envir = env, inherits = FALSE))
})

test_that("the independent slope expansion is identified and fails closed", {
  env <- load_gpcm_nonunit_score_oracle()$env
  out <- env$mfrmr_gno_expand_slopes(c(-3, -1, 1))

  expect_equal(out$log_slopes, c(-3, -1, 1, 3), tolerance = 0)
  expect_equal(out$slopes, exp(c(-3, -1, 1, 3)), tolerance = 0)
  expect_lt(out$geometric_mean_residual, 1e-15)
  expect_error(
    env$mfrmr_gno_expand_slopes(numeric(0)),
    "non-empty finite numeric vector",
    fixed = TRUE
  )
  expect_error(
    env$mfrmr_gno_expand_slopes(c(0, NA_real_)),
    "non-empty finite numeric vector",
    fixed = TRUE
  )
  expect_error(
    env$mfrmr_gno_expand_slopes(c(1000, 1000)),
    "not representable",
    fixed = TRUE
  )
})

test_that("four fixed non-unit points agree with the independent oracle", {
  env <- load_gpcm_nonunit_score_oracle()$env
  result <- env$mfrmr_run_gpcm_nonunit_score_oracle()
  audit <- result$audit
  decision <- result$decision

  expect_identical(nrow(audit), 4L)
  expect_identical(audit$Point, env$mfrmr_gno_expected_points)
  expect_true(all(audit$PrimaryRelativeStep == 1e-5))
  expect_true(all(audit$EvaluationComplete))
  expect_true(all(audit$LogProbabilityMaxAbsDifference <= 1e-12))
  expect_true(all(audit$ProbabilityMaxAbsDifference <= 1e-12))
  expect_true(all(audit$ObjectiveAbsDifference <= 1e-10))
  expect_true(all(audit$ScoreMaxAbsDifference <= 1e-6))
  expect_true(all(audit$SlopeScoreMaxAbsDifference <= 1e-6))
  expect_true(all(audit$ExpandedLogSlopeMaxAbsDifference <= 1e-12))
  expect_true(all(audit$ExpandedSlopeMaxAbsDifference <= 1e-12))
  expect_true(all(audit$GeometricMeanResidual <= 1e-12))

  stress <- audit[grepl("stress", audit$Point), , drop = FALSE]
  expect_equal(stress$MinSlope, rep(exp(-3), 2L), tolerance = 1e-15)
  expect_equal(stress$MaxSlope, rep(exp(3), 2L), tolerance = 1e-14)

  expect_true(decision$CompleteStructure)
  expect_true(decision$NumericComplete)
  expect_true(decision$NonunitOracleObserved)
  expect_identical(decision$Status, "review_oracle_agreement")
  expect_identical(decision$ScoreToleranceStatus, "pilot_required")
  expect_false(decision$BoundaryClaim)
  expect_false(decision$SelectionAuthorized)
  expect_false(decision$ConfirmationAuthorized)
  expect_identical(result$score_tolerance_status, "pilot_required")
  expect_false(result$boundary_claim)
  expect_false(result$selection_authorized)
  expect_false(result$confirmation_authorized)
})

test_that("the non-unit oracle decision rejects incomplete or mutated evidence", {
  env <- load_gpcm_nonunit_score_oracle()$env
  audit <- env$mfrmr_run_gpcm_nonunit_score_oracle()$audit

  expect_identical(env$mfrmr_gno_decision(audit[-1, ])$Status, "rejected")

  duplicate <- audit
  duplicate$Point[4] <- duplicate$Point[3]
  expect_identical(env$mfrmr_gno_decision(duplicate)$Status, "rejected")

  wrong_contract <- audit
  wrong_contract$ContractVersion <- "mutated"
  expect_identical(
    env$mfrmr_gno_decision(wrong_contract)$Status,
    "rejected"
  )

  wrong_step <- audit
  wrong_step$PrimaryRelativeStep <- 3e-5
  expect_identical(env$mfrmr_gno_decision(wrong_step)$Status, "rejected")

  incomplete <- audit
  incomplete$EvaluationComplete[1] <- FALSE
  expect_identical(env$mfrmr_gno_decision(incomplete)$Status, "rejected")

  score_drift <- audit
  score_drift$ScoreMaxAbsDifference[1] <- 2e-6
  expect_identical(env$mfrmr_gno_decision(score_drift)$Status, "rejected")

  transform_drift <- audit
  transform_drift$ExpandedSlopeMaxAbsDifference[1] <- 2e-12
  expect_identical(
    env$mfrmr_gno_decision(transform_drift)$Status,
    "rejected"
  )

  unauthorized <- audit
  unauthorized$SelectionAuthorized[1] <- TRUE
  expect_identical(env$mfrmr_gno_decision(unauthorized)$Status, "rejected")
})
