

test_that("fit-level scope records an operator but declares no path", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_parameter_path_fixture()
  audit <- fit$config$boundary_audit$gpcm_parameter_path_reachability

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-parameter-path-reachability-0.2.3-v1"
  )
  expect_identical(
    audit$state,
    "retained_additive_parameter_operator_available_no_path_declared"
  )
  expect_true(audit$evaluated)
  expect_true(audit$declared_free_coordinate_path_wrapper_available)
  expect_true(audit$forward_parameter_reachability_operator_available)
  expect_gt(audit$additive_design_rows, 0L)
  expect_gt(audit$additive_free_coordinates, 0L)
  expect_gt(audit$additive_design_nonzeros, 0)
  expect_identical(
    nrow(audit$coordinate_map), audit$additive_free_coordinates
  )
  expect_identical(audit$declared_paths_evaluated, 0L)
  expect_false(audit$production_parameter_path_reachability_checked)
  expect_false(audit$production_path_limit_classified)
  expect_true(audit$source_lexicographic_limit_contract_matches)
  expect_false(audit$inverse_projection_used)
  expect_false(audit$path_inferred_from_optimizer_trace)
  expect_false(audit$scale_exponents_inferred)
  expect_false(audit$path_search_performed)
  expect_false(audit$arbitrary_path_limit_classified)
  expect_false(audit$common_subsequence_limit_certified)
  expect_false(audit$monotone_tail_certified)
  expect_false(audit$competitive_boundary_certified)
  expect_false(audit$global_boundary_classified)
  expect_false(audit$global_finite_maximum_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
})
