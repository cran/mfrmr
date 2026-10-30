

global_existence_source_fixture <- function() {
  list(
    config = list(
      model = "GPCM", method = "JML",
      gpcm_spec = list(n_params = 2L)
    ),
    person = list(
      contract_version = mfrmr:::mfrmr_readiness_contract_version(),
      method = "JML", model = "GPCM",
      parameter_status = data.frame(
        ParameterStatus = c("estimable", "estimable"),
        stringsAsFactors = FALSE
      )
    ),
    additive = list(
      contract_version = mfrmr:::mfrmr_readiness_contract_version(),
      method = "JML", model = "GPCM",
      prescreen = list(cone_certified = FALSE),
      cone_certificate = data.frame(
        Certified = FALSE,
        StrictContrastRows = 0L,
        DirectionL1 = 0,
        stringsAsFactors = FALSE
      )
    ),
    sequence = list(
      contract_version =
        "mfrmr-jml-gpcm-parameter-sequence-flag-0.2.3-v1",
      every_parameter_sequence_has_classifiable_further_subsequence_likelihood_limit =
        TRUE
    ),
    slope = list(
      contract_version = mfrmr:::mfrmr_readiness_contract_version(),
      method = "JML", model = "GPCM", certificates = data.frame()
    ),
    joint = list(
      contract_version = mfrmr:::mfrmr_readiness_contract_version(),
      method = "JML", model = "GPCM",
      certificates = data.frame(), rate_certificates = data.frame()
    )
  )
}

call_global_existence_audit <- function(source, ...) {
  mfrmr:::audit_mfrm_jml_gpcm_global_existence(
    config = source$config,
    person_boundary = source$person,
    joint_additive = source$additive,
    parameter_sequence_flag = source$sequence,
    slope_boundary = source$slope,
    joint_slope_boundary = source$joint,
    retained_log_likelihood = -10,
    ...
  )
}


test_that("production free extreme Persons certify finite JMLE nonexistence", {
  source <- global_existence_source_fixture()
  source$person$parameter_status$ParameterStatus[1] <- "unbounded_low"
  audit <- call_global_existence_audit(source)

  expect_identical(
    audit$state,
    "finite_jmle_nonexistence_certified_strict_additive_recession"
  )
  expect_true(audit$evaluated)
  expect_true(audit$source_parameter_sequence_flag_contract_matches)
  expect_true(
    audit$every_parameter_sequence_has_classifiable_further_subsequence_likelihood_limit
  )
  expect_identical(audit$free_extreme_person_recession_count, 1L)
  expect_false(audit$global_additive_recession_cone_certified)
  expect_true(audit$global_boundary_competitiveness_certified)
  expect_false(audit$global_supremum_value_classified)
  expect_true(is.na(audit$global_log_likelihood_supremum))
  expect_true(audit$finite_jmle_existence_classified)
  expect_false(audit$finite_jmle_existence_certified)
  expect_true(audit$finite_jmle_nonexistence_certified)
  expect_true(audit$global_boundary_classified)
  expect_true(audit$global_boundary_presence_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$complete_boundary_envelope_constructed)
  expect_false(audit$finite_over_complete_boundary_gap_certified)
  expect_false(audit$production_parameter_path_declared)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
})

test_that("a certified global additive cone proves nonexistence", {
  source <- global_existence_source_fixture()
  source$additive$cone_certificate <- data.frame(
    Certified = TRUE,
    StrictContrastRows = 3L,
    DirectionL1 = 2,
    stringsAsFactors = FALSE
  )
  source$additive$prescreen$cone_certified <- TRUE
  audit <- call_global_existence_audit(source)

  expect_identical(
    audit$state,
    "finite_jmle_nonexistence_certified_strict_additive_recession"
  )
  expect_identical(audit$free_extreme_person_recession_count, 0L)
  expect_true(audit$global_additive_recession_cone_certified)
  expect_true(audit$finite_jmle_nonexistence_certified)
  expect_match(audit$reason_codes, "global_additive_strict_recession_cone")

  source$additive$cone_certificate$StrictContrastRows <- 0L
  zero_margin <- call_global_existence_audit(source)
  expect_false(zero_margin$global_additive_recession_cone_certified)
  expect_false(zero_margin$finite_jmle_nonexistence_certified)
})

test_that("an analytic zero boundary classifies the global supremum", {
  source <- global_existence_source_fixture()
  source$slope$certificates <- data.frame(
    BoundaryLogLikelihood = c(-2, 0),
    Certified = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  audit <- call_global_existence_audit(source)

  expect_identical(
    audit$state,
    "global_supremum_zero_finite_jmle_nonexistence_certified"
  )
  expect_identical(
    audit$universal_upper_bound_zero_path_source,
    "slope_only_analytic_path"
  )
  expect_true(audit$global_supremum_value_classified)
  expect_identical(audit$global_log_likelihood_supremum, 0)
  expect_true(audit$finite_jmle_nonexistence_certified)

  source <- global_existence_source_fixture()
  source$joint$rate_certificates <- data.frame(
    BoundaryLogLikelihood = 0,
    Certified = TRUE,
    AnalyticTailCertified = TRUE,
    stringsAsFactors = FALSE
  )
  joint <- call_global_existence_audit(source)
  expect_identical(
    joint$universal_upper_bound_zero_path_source,
    "joint_slope_additive_analytic_path"
  )
  expect_true(joint$global_supremum_value_classified)

  source$joint$rate_certificates$AnalyticTailCertified <- FALSE
  unverified <- call_global_existence_audit(source)
  expect_false(unverified$global_supremum_value_classified)
})

test_that("absence of a positive certificate leaves finite existence open", {
  source <- global_existence_source_fixture()
  audit <- call_global_existence_audit(source)

  expect_identical(
    audit$state,
    "global_finite_jmle_existence_open_complete_boundary_envelope_unavailable"
  )
  expect_true(audit$evaluated)
  expect_true(audit$retained_point_finite)
  expect_true(audit$retained_point_is_optimizer_trace_only)
  expect_false(audit$global_boundary_competitiveness_certified)
  expect_false(audit$finite_jmle_existence_classified)
  expect_false(audit$finite_jmle_existence_certified)
  expect_false(audit$finite_jmle_nonexistence_certified)
  expect_false(audit$global_boundary_classified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$complete_boundary_envelope_constructed)
  expect_match(audit$reason_codes, "complete_boundary_envelope_unavailable")
})

test_that("current fits retain typed P2f states without readiness promotion", {
  skip_if_not_installed("lpSolve")
  crossed <- fit_jml_gpcm_global_existence_fixture("crossed")
  open <- crossed$config$boundary_audit$gpcm_global_existence
  expect_identical(
    open$contract_version,
    "mfrmr-jml-gpcm-global-existence-0.2.3-v1"
  )
  expect_true(open$evaluated)
  expect_true(open$source_parameter_sequence_flag_contract_matches)
  expect_false(open$finite_jmle_existence_certified)
  expect_false(open$global_boundary_absence_certified)
  expect_false(open$external_comparison_eligible)
  expect_identical(open$readiness_effect, "none_diagnostic_only")

  extreme <- fit_jml_gpcm_global_existence_fixture("all_extreme")
  absent <- extreme$config$boundary_audit$gpcm_global_existence
  expect_true(absent$evaluated)
  expect_true(absent$finite_jmle_nonexistence_certified)
  expect_true(absent$global_boundary_classified)
  expect_false(absent$finite_jmle_existence_certified)
  expect_false(absent$standard_error_eligible)
  expect_false(absent$confidence_interval_eligible)
  expect_false(absent$external_comparison_eligible)
})

test_that("model estimator unit-slope and source boundaries fail closed", {
  source <- global_existence_source_fixture()

  pcm <- source
  pcm$config$model <- "PCM"
  pcm_result <- call_global_existence_audit(pcm)
  expect_identical(pcm_result$state, "not_applicable_model")
  expect_false(pcm_result$evaluated)
  expect_identical(pcm_result$estimator_identity, "not_applicable")

  mml <- source
  mml$config$method <- "MML"
  mml_result <- call_global_existence_audit(mml)
  expect_identical(mml_result$state, "not_applicable_estimator")
  expect_false(mml_result$evaluated)
  expect_identical(
    mml_result$objective_identity, "not_applicable_marginal_estimator"
  )

  unit <- source
  unit$config$gpcm_spec$n_params <- 0L
  unit_result <- call_global_existence_audit(unit)
  expect_identical(unit_result$state, "not_required_unit_slope")
  expect_false(unit_result$evaluated)
  expect_false(unit_result$finite_jmle_nonexistence_certified)

  legacy <- source
  legacy$sequence$contract_version <- "legacy"
  legacy_result <- call_global_existence_audit(legacy)
  expect_identical(legacy_result$state, "not_evaluated_source_contract")
  expect_false(legacy_result$evaluated)
  expect_false(legacy_result$source_parameter_sequence_flag_contract_matches)
  expect_false(legacy_result$finite_jmle_nonexistence_certified)
  expect_false(legacy_result$global_boundary_classified)
})
