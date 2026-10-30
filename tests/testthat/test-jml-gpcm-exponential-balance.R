

test_that("production audit certifies the canonical nonproperness witness", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_exponential_balance_fixture()
  audit <- fit$config$boundary_audit$gpcm_exponential_balance

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-exponential-balance-0.2.3-v1"
  )
  expect_identical(
    audit$state, "canonical_zero_utility_bounded_image_escape_certified"
  )
  expect_true(audit$evaluated)
  expect_true(audit$declared_exponential_balance_oracle_available)
  expect_true(audit$exact_combined_exponent_aggregation_available)
  expect_true(audit$bounded_response_image_escape_witness_classifiable)
  expect_true(audit$source_parameter_sequence_flag_contract_matches)
  expect_true(audit$source_global_existence_contract_matches)
  expect_true(audit$canonical_zero_utility_center_certified)
  expect_true(audit$canonical_affine_additive_offset_exactly_zero)
  expect_true(audit$canonical_zero_utility_divergent_slope_path_certified)
  expect_true(audit$canonical_path_additive_coordinates_identically_zero)
  expect_true(audit$canonical_path_response_contrasts_identically_zero)
  expect_identical(audit$production_exponential_paths_declared, 0L)
  expect_identical(audit$production_canonical_paths_evaluated, 1L)
  expect_false(audit$production_bounded_image_escape_search_performed)
  expect_true(audit$production_bounded_image_escape_certified)
  expect_true(audit$parameter_to_response_contrast_map_properness_classified)
  expect_true(audit$parameter_to_response_contrast_map_properness_refuted)
  expect_false(audit$complete_bounded_image_escape_family_classified)
  expect_false(audit$complete_escaping_sequence_boundary_envelope_constructed)
  expect_false(audit$finite_jmle_existence_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
})
