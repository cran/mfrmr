

test_that("current fits retain the further-subsequence theorem but no flag", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_parameter_sequence_flag_fixture()
  audit <- fit$config$boundary_audit$gpcm_parameter_sequence_flag

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-parameter-sequence-flag-0.2.3-v1"
  )
  expect_identical(
    audit$state,
    "parameter_sequence_further_subsequence_flag_theorem_available_no_sequence_declared"
  )
  expect_true(audit$evaluated)
  expect_gt(audit$free_parameter_dimension, 0L)
  expect_identical(audit$observations, nrow(fit$prep$data))
  expect_identical(audit$categories, as.integer(fit$config$n_cat))
  expect_identical(
    audit$response_contrast_dimension,
    nrow(fit$prep$data) * (fit$config$n_cat - 1)
  )
  expect_true(audit$finite_dimensional_contrast_image)
  expect_true(
    audit$exact_nonlinear_slope_utility_map_precedes_compactification
  )
  expect_true(audit$bounded_contrast_sequence_has_convergent_subsequence)
  expect_true(
    audit$unbounded_contrast_sequence_has_scale_separated_flag_subsequence
  )
  expect_true(audit$flag_directions_can_be_chosen_orthonormal)
  expect_true(audit$successive_flag_scale_ratios_tend_to_zero)
  expect_true(audit$terminal_reference_contrast_remainder_tends_to_zero)
  expect_identical(
    audit$maximum_divergent_flag_stages,
    audit$response_contrast_dimension
  )
  expect_true(
    audit$every_parameter_sequence_has_classifiable_further_subsequence_likelihood_limit
  )
  expect_true(audit$declared_contrast_flag_oracle_available)
  expect_true(audit$source_compactification_contract_matches)
  expect_true(audit$source_finite_sequence_contract_matches)
  expect_identical(audit$production_parameter_sequences_declared, 0L)
  expect_false(audit$production_parameter_sequence_mapped)
  expect_false(audit$production_subsequence_selected)
  expect_false(audit$production_flag_extracted)
  expect_false(
    audit$optimizer_stage_endpoints_used_as_asymptotic_sequence
  )
  expect_false(audit$original_sequence_limit_classified)
  expect_false(audit$common_subsequence_limit_certified)
  expect_false(audit$all_subsequence_limits_equal_certified)
  expect_false(audit$competitive_boundary_certified)
  expect_false(audit$global_boundary_classified)
  expect_false(audit$global_finite_maximum_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
})
