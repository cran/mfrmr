test_that("ordered owner-rate hierarchy counts and generators are exact", {
  expected <- c(1, 3, 13, 75, 541, 4683, 47293)
  observed <- vapply(seq_along(expected), function(elements) {
    mfrmr:::mfrmr_jml_gpcm_higher_order_lift_partition_count(elements)
  }, numeric(1))
  expect_identical(observed, expected)
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_higher_order_lift_partition_count(3, 1),
    1
  )
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_higher_order_lift_partition_count(3, 2),
    7
  )
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_higher_order_lift_partition_count(0),
    0
  )

  one_stage <- mfrmr:::mfrmr_jml_gpcm_higher_order_lift_surjections(3, 1)
  expect_length(one_stage, 1L)
  expect_identical(one_stage[[1]], rep(1L, 3))
  two_stage <- mfrmr:::mfrmr_jml_gpcm_higher_order_lift_surjections(3, 2)
  expect_length(two_stage, 6L)
  expect_true(all(vapply(two_stage, function(value) {
    identical(sort(unique(value)), 1:2)
  }, logical(1))))
  expect_identical(length(unique(vapply(
    two_stage,
    paste,
    character(1),
    collapse = ","
  ))), 6L)
  three_stage <- mfrmr:::mfrmr_jml_gpcm_higher_order_lift_surjections(3, 3)
  expect_length(three_stage, 6L)
  expect_true(all(vapply(three_stage, function(value) {
    identical(sort(value), 1:3)
  }, logical(1))))
  expect_identical(
    length(one_stage) + length(two_stage) + length(three_stage),
    13L
  )
})

test_that("finite and outer-excluded P2j sources short-circuit safely", {
  skip_if_not_installed("lpSolve")

  finite <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_p2i_fixture(
    c(2, 1, 4, 1)
  )
  expect_identical(
    finite$contract_version,
    "mfrmr-jml-gpcm-higher-order-face-lifts-0.2.3-v1"
  )
  expect_identical(
    finite$state,
    "finite_source_response_image_closure_certified"
  )
  expect_true(finite$evaluated)
  expect_true(finite$higher_order_face_lifts_classified)
  expect_true(finite$source_face_chart_contract_matches)
  expect_true(finite$source_face_chart_classified)
  expect_true(finite$zero_affine_offset_certified)
  expect_true(finite$semialgebraic_curve_selection_reduction_derived)
  expect_true(finite$puiseux_owner_leading_rate_reduction_derived)
  expect_true(
    finite$distinct_owner_rates_compressible_to_ordered_integer_stages
  )
  expect_true(finite$hierarchy_linear_conditions_necessary_for_closure)
  expect_true(
    finite$hierarchy_linear_conditions_sufficient_by_explicit_path
  )
  expect_false(finite$hierarchy_enumeration_required)
  expect_true(finite$all_required_ordered_owner_rate_hierarchies_enumerated)
  expect_true(
    finite$all_required_owner_rate_hierarchies_numerically_decided
  )
  expect_identical(finite$hierarchy_count, 0L)
  expect_identical(finite$expected_hierarchy_count, 0)
  expect_true(finite$finite_response_image_certified)
  expect_true(finite$target_in_finite_response_image_closure_certified)
  expect_false(finite$reachable_missing_finite_boundary_certified)
  expect_false(
    finite$target_outside_finite_response_image_closure_certified
  )
  expect_true(finite$targetwise_zero_offset_closure_membership_classified)
  expect_false(finite$symbolic_entire_operator_closure_stratified)
  expect_false(finite$complete_general_boundary_likelihood_envelope_constructed)
  expect_false(finite$finite_jmle_existence_certified)
  expect_false(finite$finite_jmle_nonexistence_certified)
  expect_false(finite$standard_error_eligible)
  expect_false(finite$confidence_interval_eligible)
  expect_false(finite$external_comparison_eligible)
  expect_identical(finite$readiness_effect, "none_theorem_only")
  expect_identical(
    finite$next_gate,
    "general_boundary_likelihood_envelope_over_rate_hierarchies"
  )

  outside <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_p2i_fixture(
    c(1, 2, 2, 1)
  )
  expect_identical(
    outside$state,
    "outside_closure_certified_by_source_simplex_outer_bound"
  )
  expect_true(outside$higher_order_face_lifts_classified)
  expect_false(outside$finite_response_image_certified)
  expect_false(outside$target_in_finite_response_image_closure_certified)
  expect_true(outside$target_outside_finite_response_image_closure_certified)
  expect_true(outside$targetwise_zero_offset_closure_membership_classified)
  expect_identical(outside$hierarchy_count, 0L)

  one_owner_design <- rbind(c(1, 0), c(0, 1), c(1, 1))
  one_owner_source <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    one_owner_design,
    rep(1, 3),
    c(1, 2, 4)
  )
  one_owner <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(
    one_owner_source
  )
  expect_identical(
    one_owner$state,
    "outside_closure_certified_by_source_simplex_outer_bound"
  )
  expect_true(one_owner$target_outside_finite_response_image_closure_certified)
})

test_that("higher-order chart recovers both minimal P2i missing axes", {
  skip_if_not_installed("lpSolve")

  owner_1 <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_p2i_fixture(
    c(log(2), log(4), 0, 0)
  )
  expect_identical(
    owner_1$state,
    "missing_boundary_hierarchy_lift_certified"
  )
  expect_true(owner_1$higher_order_face_lifts_classified)
  expect_true(owner_1$hierarchy_enumeration_required)
  expect_true(owner_1$all_required_ordered_owner_rate_hierarchies_enumerated)
  expect_true(
    owner_1$all_required_owner_rate_hierarchies_numerically_decided
  )
  expect_identical(owner_1$hierarchy_count, 1L)
  expect_identical(owner_1$expected_hierarchy_count, 1)
  expect_identical(owner_1$certified_hierarchy_count, 1L)
  expect_identical(owner_1$certified_first_order_hierarchy_count, 1L)
  expect_identical(owner_1$certified_higher_order_hierarchy_count, 0L)
  expect_false(owner_1$finite_response_image_certified)
  expect_true(owner_1$finite_response_image_excluded)
  expect_true(owner_1$target_in_finite_response_image_closure_certified)
  expect_true(owner_1$reachable_missing_finite_boundary_certified)
  expect_false(
    owner_1$target_outside_finite_response_image_closure_certified
  )
  expect_identical(owner_1$hierarchy_table$HierarchyId, "rates_1_0")
  expect_identical(owner_1$hierarchy_table$RateCode, "1,0")
  expect_true(owner_1$hierarchy_table$FirstOrderHierarchy)
  expect_true(owner_1$hierarchy_table$SourceFaceFirstOrderReachable)
  expect_true(owner_1$hierarchy_table$HierarchyLiftCertified)
  expect_gt(owner_1$hierarchy_table$MaximumPositiveMargin, 0)
  expect_lte(
    owner_1$hierarchy_table$MaximumAbsoluteEquationResidual,
    1e-12
  )

  owner_2 <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_p2i_fixture(
    c(0, 0, 1, 2)
  )
  expect_identical(
    owner_2$state,
    "missing_boundary_hierarchy_lift_certified"
  )
  expect_identical(owner_2$hierarchy_count, 1L)
  expect_identical(owner_2$certified_hierarchy_count, 1L)
  expect_identical(owner_2$hierarchy_table$HierarchyId, "rates_0_1")
  expect_identical(owner_2$hierarchy_table$RateCode, "0,1")
  expect_true(owner_2$hierarchy_table$HierarchyLiftCertified)
})

test_that("a hierarchy certificate gives a finite product-one path", {
  skip_if_not_installed("lpSolve")

  target <- c(log(2), log(4), 0, 0)
  result <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_p2i_fixture(target)
  indices <- c(0, 1, 3, 8, 20)
  paths <- lapply(indices, function(index) {
    mfrmr:::mfrmr_jml_gpcm_higher_order_face_path_at(
      result,
      "rates_1_0",
      index
    )
  })
  expect_true(all(vapply(paths, `[[`, logical(1), "valid")))
  expect_true(all(vapply(paths, `[[`, logical(1), "finite_parameter_point")))
  expect_true(all(vapply(
    paths,
    `[[`,
    logical(1),
    "parameter_path_declared_divergent"
  )))
  expect_true(all(vapply(
    paths,
    `[[`,
    logical(1),
    "response_contrast_limit_certified"
  )))
  expect_true(all(vapply(paths, function(path) {
    abs(sum(path$log_slopes)) <= 1e-12
  }, logical(1))))
  expect_true(all(vapply(paths, function(path) {
    abs(prod(path$inverse_slopes) - 1) <= 1e-12
  }, logical(1))))
  expect_true(all(vapply(paths, function(path) {
    path$maximum_absolute_raw_equation_residual <= 1e-14
  }, logical(1))))
  expect_true(all(vapply(paths, function(path) {
    mfrmr:::mfrmr_jml_gpcm_binary_closure_image_state(
      path$contrasts
    )$finite_response_image
  }, logical(1))))
  distances <- vapply(
    paths,
    `[[`,
    numeric(1),
    "contrast_distance_to_target"
  )
  expect_true(all(diff(distances) < 0))
  expect_lt(tail(distances, 1), 4e-09)
  expect_equal(tail(paths, 1)[[1]]$contrasts, target, tolerance = 4e-09)
  expect_gt(max(abs(tail(paths, 1)[[1]]$log_slopes)), 9)
})

test_that("higher stages resolve first-order-open rank-one faces", {
  skip_if_not_installed("lpSolve")

  design <- matrix(1, nrow = 3, ncol = 1)
  source <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    1:3,
    c(1, 0, 0)
  )
  expect_identical(
    source$state,
    "reachable_missing_boundary_face_certified"
  )
  result <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(source)
  expect_identical(
    result$state,
    "missing_boundary_hierarchy_lift_certified"
  )
  expect_identical(result$hierarchy_count, 7L)
  expect_identical(result$expected_hierarchy_count, 7)
  expect_identical(result$certified_hierarchy_count, 3L)
  expect_identical(result$certified_first_order_hierarchy_count, 1L)
  expect_identical(result$certified_higher_order_hierarchy_count, 2L)
  expect_true(result$source_first_order_open_face_resolved_by_higher_order)
  expect_true(result$target_in_finite_response_image_closure_certified)
  expect_false(result$target_outside_finite_response_image_closure_certified)
  expect_identical(
    result$hierarchy_table$HierarchyId,
    c(
      "rates_1_0_0",
      "rates_1_0_1",
      "rates_1_0_2",
      "rates_2_0_1",
      "rates_1_1_0",
      "rates_1_2_0",
      "rates_2_1_0"
    )
  )
  expect_identical(
    result$hierarchy_table$HierarchyId[
      result$hierarchy_table$HierarchyLiftCertified
    ],
    c("rates_1_0_0", "rates_2_0_1", "rates_2_1_0")
  )
  higher <- result$hierarchy_table$MaximumRate > 1L &
    result$hierarchy_table$HierarchyLiftCertified
  expect_true(all(!result$hierarchy_table$SourceFaceFirstOrderReachable[higher]))
  expect_true(all(result$hierarchy_table$StrictHierarchyExcluded[!result$hierarchy_table$HierarchyLiftCertified]))

  indices <- c(0, 1, 3, 8, 15)
  paths <- lapply(indices, function(index) {
    mfrmr:::mfrmr_jml_gpcm_higher_order_face_path_at(
      result,
      "rates_2_0_1",
      index
    )
  })
  expect_true(all(vapply(paths, `[[`, logical(1), "valid")))
  expect_true(all(vapply(paths, function(path) {
    identical(path$rates, c(2L, 0L, 1L))
  }, logical(1))))
  distances <- vapply(
    paths,
    `[[`,
    numeric(1),
    "contrast_distance_to_target"
  )
  expect_true(all(diff(distances) < 0))
  expect_lt(tail(distances, 1), 4e-07)
  expect_equal(
    tail(paths, 1)[[1]]$contrasts,
    c(1, 0, 0),
    tolerance = 4e-07
  )
  expect_gt(max(abs(tail(paths, 1)[[1]]$log_slopes)), 14)
})

test_that("all higher hierarchies exclude the P2j outer-only counterexample", {
  skip_if_not_installed("lpSolve")

  design <- matrix(1, nrow = 3, ncol = 1)
  source <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    1:3,
    c(1, -2, 0)
  )
  expect_identical(
    source$state,
    "outer_boundary_face_candidate_lift_open"
  )
  expect_true(source$nonnegative_simplex_face_witness_exists)
  expect_true(source$response_image_closure_membership_open)

  result <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(source)
  expect_identical(
    result$state,
    "outside_closure_all_rate_hierarchies_excluded"
  )
  expect_true(result$higher_order_face_lifts_classified)
  expect_true(result$all_required_ordered_owner_rate_hierarchies_enumerated)
  expect_identical(result$hierarchy_count, 3L)
  expect_identical(result$expected_hierarchy_count, 3)
  expect_identical(
    result$hierarchy_table$HierarchyId,
    c("rates_1_1_0", "rates_1_2_0", "rates_2_1_0")
  )
  expect_false(any(result$hierarchy_table$HierarchyLiftCertified))
  expect_true(all(result$hierarchy_table$StrictHierarchyExcluded))
  expect_false(any(result$hierarchy_table$NumericalMarginOpen))
  expect_false(result$target_in_finite_response_image_closure_certified)
  expect_false(result$reachable_missing_finite_boundary_certified)
  expect_true(result$target_outside_finite_response_image_closure_certified)
  expect_true(result$targetwise_zero_offset_closure_membership_classified)
  expect_false(result$nonnegative_face_condition_sufficient_for_general_closure)
  expect_false(result$finite_jmle_nonexistence_certified)
})

test_that("small positive hierarchy margins remain open until certified", {
  skip_if_not_installed("lpSolve")

  source <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    matrix(1, nrow = 3, ncol = 1),
    1:3,
    c(1, 1e6, 0)
  )
  expect_true(source$response_image_face_chart_classified)
  expect_true(source$reachable_missing_finite_boundary_certified)

  guarded <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(
    source,
    certificate_tolerance = 1e-05
  )
  expect_identical(
    guarded$state,
    "source_closure_certified_hierarchy_numerical_margin_open"
  )
  expect_true(guarded$evaluated)
  expect_false(guarded$higher_order_face_lifts_classified)
  expect_true(guarded$source_closure_membership_certified)
  expect_true(guarded$target_in_finite_response_image_closure_certified)
  expect_true(guarded$all_required_ordered_owner_rate_hierarchies_enumerated)
  expect_false(
    guarded$all_required_owner_rate_hierarchies_numerically_decided
  )
  expect_false(guarded$target_outside_finite_response_image_closure_certified)
  expect_true(any(guarded$hierarchy_table$NumericalMarginOpen))
  expect_gt(
    max(guarded$hierarchy_table$MaximumPositiveMargin, na.rm = TRUE),
    0
  )
  expect_lt(
    max(guarded$hierarchy_table$MaximumPositiveMargin, na.rm = TRUE),
    1e-05
  )

  resolved <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(
    source,
    certificate_tolerance = 1e-08
  )
  expect_identical(
    resolved$state,
    "missing_boundary_hierarchy_lift_certified"
  )
  expect_true(resolved$higher_order_face_lifts_classified)
  expect_true(resolved$target_in_finite_response_image_closure_certified)
  expect_identical(resolved$certified_hierarchy_count, 1L)
})

test_that("hierarchy classification is invariant to response-row order", {
  skip_if_not_installed("lpSolve")

  design <- matrix(1, nrow = 3, ncol = 1)
  owner <- 1:3
  target <- c(1, 0, 0)
  source <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    owner,
    target
  )
  result <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(source)

  permutation <- c(3, 1, 2)
  permuted_source <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design[permutation, , drop = FALSE],
    owner[permutation],
    target[permutation]
  )
  permuted <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(
    permuted_source
  )
  expect_identical(permuted$state, result$state)
  expect_identical(permuted$hierarchy_count, result$hierarchy_count)
  expect_identical(
    permuted$hierarchy_table$HierarchyId,
    result$hierarchy_table$HierarchyId
  )
  expect_identical(
    permuted$hierarchy_table$HierarchyLiftCertified,
    result$hierarchy_table$HierarchyLiftCertified
  )
  expect_identical(
    permuted$hierarchy_table$StrictHierarchyExcluded,
    result$hierarchy_table$StrictHierarchyExcluded
  )
})

test_that("higher-order source, workload, solver, and path controls fail closed", {
  skip_if_not_installed("lpSolve")

  source <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    matrix(1, nrow = 3, ncol = 1),
    1:3,
    c(1, 0, 0)
  )
  control_cases <- list(
    list(certificate_tolerance = 0),
    list(max_stages = 0),
    list(max_hierarchies = 0),
    list(max_equations = 0),
    list(max_lp_nonzeros = 0),
    list(lp_timeout = 0)
  )
  for (controls in control_cases) {
    result <- do.call(
      mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts,
      c(list(source), controls)
    )
    expect_identical(result$state, "not_evaluated_control")
    expect_false(result$evaluated)
    expect_false(result$higher_order_face_lifts_classified)
    expect_false(result$targetwise_zero_offset_closure_membership_classified)
  }

  invalid_sources <- list(
    1,
    within(source, contract_version <- "legacy"),
    within(source, response_image_face_chart_classified <- FALSE),
    within(source, zero_affine_offset_certified <- FALSE)
  )
  for (invalid in invalid_sources) {
    result <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(invalid)
    expect_identical(result$state, "not_evaluated_source")
    expect_false(result$evaluated)
    expect_false(result$higher_order_face_lifts_classified)
  }

  near_source <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    adjacent_design = rbind(
      c(1, 0, -1),
      c(0, 1, -1),
      c(1, 0, 1),
      c(0, 1, 1)
    ),
    slope_owner_by_row = c(1, 1, 2, 2),
    target_contrasts = c(1, 0, 1e-06, 0),
    certificate_tolerance = 1e-05
  )
  expect_false(near_source$response_image_face_chart_classified)
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(near_source)$state,
    "not_evaluated_source"
  )

  stage_cap <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(
    source,
    max_stages = 1
  )
  expect_identical(stage_cap$state, "not_evaluated_stage_cap")
  expect_false(stage_cap$evaluated)

  hierarchy_cap <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(
    source,
    max_hierarchies = 6
  )
  expect_identical(
    hierarchy_cap$state,
    "not_evaluated_hierarchy_workload"
  )
  expect_false(hierarchy_cap$evaluated)
  expect_identical(hierarchy_cap$expected_hierarchy_count, 7)

  equation_cap <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(
    source,
    max_equations = 1
  )
  expect_identical(equation_cap$state, "not_evaluated_solver")
  expect_true(equation_cap$evaluated)
  expect_false(equation_cap$higher_order_face_lifts_classified)

  nonzero_cap <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(
    source,
    max_lp_nonzeros = 1
  )
  expect_identical(nonzero_cap$state, "not_evaluated_solver")
  expect_true(nonzero_cap$evaluated)
  expect_false(nonzero_cap$higher_order_face_lifts_classified)

  valid <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(source)
  invalid_path_sources <- list(
    1,
    within(valid, contract_version <- "legacy"),
    valid
  )
  invalid_ids <- c("rates_2_0_1", "rates_2_0_1", "unknown")
  for (position in seq_along(invalid_path_sources)) {
    path <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_path_at(
      invalid_path_sources[[position]],
      invalid_ids[position],
      1
    )
    expect_false(path$valid)
    expect_identical(path$state, "not_evaluated_higher_order_path")
  }
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_higher_order_face_path_at(
      valid,
      "rates_2_0_1",
      -1
    )$valid
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_higher_order_face_path_at(
      valid,
      "rates_2_0_1",
      101
    )$valid
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_higher_order_face_path_at(
      valid,
      "rates_2_0_1",
      10,
      max_log_rate_span = 10
    )$valid
  )

  finite <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_p2i_fixture(
    c(2, 1, 4, 1)
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_higher_order_face_path_at(
      finite,
      "rates_0_0",
      1
    )$valid
  )
  excluded_source <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    matrix(1, nrow = 3, ncol = 1),
    1:3,
    c(1, -2, 0)
  )
  excluded <- mfrmr:::mfrmr_jml_gpcm_higher_order_face_lifts(
    excluded_source
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_higher_order_face_path_at(
      excluded,
      "rates_1_1_0",
      1
    )$valid
  )
})
