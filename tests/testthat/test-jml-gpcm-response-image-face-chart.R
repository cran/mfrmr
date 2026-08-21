test_that("general face chart reconstructs finite P2i response points", {
  skip_if_not_installed("lpSolve")

  result <- mfrmr:::mfrmr_jml_gpcm_response_image_face_p2i_fixture(
    c(2, 1, 4, 1)
  )
  expect_identical(
    result$contract_version,
    "mfrmr-jml-gpcm-response-image-face-chart-0.2.3-v1"
  )
  expect_identical(
    result$state,
    "finite_response_image_representative_certified"
  )
  expect_true(result$evaluated)
  expect_true(result$response_image_face_chart_classified)
  expect_true(result$zero_affine_offset_certified)
  expect_true(result$finite_image_positive_inverse_slope_equivalence_derived)
  expect_true(result$inverse_slope_simplex_compactification_derived)
  expect_true(result$every_finite_closure_point_has_nonnegative_face_witness)
  expect_true(result$proper_simplex_face_necessary_for_missing_boundary)
  expect_true(result$first_order_face_lift_sufficiency_derived)
  expect_false(result$nonnegative_face_condition_sufficient_for_general_closure)
  expect_true(result$all_nonempty_simplex_faces_enumerated)
  expect_identical(result$face_count, 3L)
  expect_identical(result$rows, 4L)
  expect_identical(result$additive_free_coordinates, 3L)
  expect_identical(result$slope_owners, 2L)
  expect_true(result$finite_response_image_certified)
  expect_false(result$finite_response_image_excluded_by_full_support_lp)
  expect_true(result$nonnegative_simplex_face_witness_exists)
  expect_true(result$target_in_finite_response_image_closure_certified)
  expect_false(result$reachable_missing_finite_boundary_certified)
  expect_false(result$target_outside_closure_outer_bound_certified)
  expect_false(result$response_image_closure_membership_open)
  expect_false(result$complete_general_response_image_closure_stratified)
  expect_false(result$complete_general_boundary_likelihood_envelope_constructed)
  expect_false(result$finite_jmle_existence_certified)
  expect_false(result$finite_jmle_nonexistence_certified)
  expect_false(result$standard_error_eligible)
  expect_false(result$confidence_interval_eligible)
  expect_false(result$external_comparison_eligible)
  expect_identical(result$readiness_effect, "none_theorem_only")
  expect_identical(
    result$next_gate,
    "higher_order_face_lifts_and_general_closure_completion"
  )

  representative <- result$finite_representative
  expect_equal(prod(representative$inverse_slopes), 1, tolerance = 1e-14)
  expect_equal(prod(representative$slopes), 1, tolerance = 1e-14)
  expect_equal(sum(representative$log_slopes), 0, tolerance = 1e-14)
  expect_equal(
    representative$reconstructed_contrasts,
    c(2, 1, 4, 1),
    tolerance = 1e-12
  )
  expect_lte(
    representative$maximum_absolute_reconstruction_residual,
    1e-12
  )
  expect_identical(
    result$face_table$FaceId,
    c("support_1_2", "support_1", "support_2")
  )
  expect_identical(
    result$face_table$RelativeInteriorCompatible,
    c(TRUE, FALSE, FALSE)
  )
  expect_equal(
    result$face_table$MaximumRelativeMargin[1],
    0.25,
    tolerance = 1e-12
  )

  negative <- mfrmr:::mfrmr_jml_gpcm_response_image_face_p2i_fixture(
    c(1, 2, 1, 4)
  )
  expect_true(negative$finite_response_image_certified)
  expect_equal(
    negative$finite_representative$reconstructed_contrasts,
    c(1, 2, 1, 4),
    tolerance = 1e-12
  )

  equal <- mfrmr:::mfrmr_jml_gpcm_response_image_face_p2i_fixture(
    c(3, 3, -2, -2)
  )
  expect_true(equal$finite_response_image_certified)
  expect_true(all(equal$face_table$RelativeInteriorCompatible))
  expect_identical(equal$compatible_proper_face_count, 2L)
  expect_identical(equal$first_order_reachable_face_count, 2L)
})

test_that("simplex chart reproduces every exact P2i image state", {
  skip_if_not_installed("lpSolve")

  targets <- list(
    c(2, 1, 4, 1),
    c(1, 2, 1, 4),
    c(3, 3, -2, -2),
    c(1, 2, 0, 0),
    c(0, 0, 1, 2),
    c(1, 2, 2, 1)
  )
  binary <- lapply(
    targets,
    mfrmr:::mfrmr_jml_gpcm_binary_closure_image_state
  )
  charts <- lapply(
    targets,
    mfrmr:::mfrmr_jml_gpcm_response_image_face_p2i_fixture
  )

  expect_identical(
    vapply(charts, `[[`, logical(1), "finite_response_image_certified"),
    vapply(binary, `[[`, logical(1), "finite_response_image")
  )
  expect_identical(
    vapply(
      charts,
      `[[`,
      logical(1),
      "target_in_finite_response_image_closure_certified"
    ),
    vapply(binary, `[[`, logical(1), "finite_response_image_closure")
  )
  expect_identical(
    vapply(
      charts,
      `[[`,
      logical(1),
      "reachable_missing_finite_boundary_certified"
    ),
    vapply(binary, `[[`, logical(1), "missing_finite_boundary")
  )
  expect_identical(
    vapply(
      charts,
      `[[`,
      logical(1),
      "target_outside_closure_outer_bound_certified"
    ),
    !vapply(binary, `[[`, logical(1), "finite_response_image_closure")
  )
  expect_identical(
    vapply(charts, `[[`, character(1), "state"),
    c(
      "finite_response_image_representative_certified",
      "finite_response_image_representative_certified",
      "finite_response_image_representative_certified",
      "reachable_missing_boundary_face_certified",
      "reachable_missing_boundary_face_certified",
      "outside_response_image_closure_outer_bound_certified"
    )
  )

  owner_1_axis <- charts[[4]]
  expect_identical(
    owner_1_axis$face_table$FaceId[
      owner_1_axis$face_table$RelativeInteriorCompatible
    ],
    "support_2"
  )
  expect_identical(
    owner_1_axis$face_table$FaceId[
      owner_1_axis$face_table$FirstOrderBoundaryReachable
    ],
    "support_2"
  )
  expect_true(owner_1_axis$finite_response_image_excluded_by_full_support_lp)

  owner_2_axis <- charts[[5]]
  expect_identical(
    owner_2_axis$face_table$FaceId[
      owner_2_axis$face_table$RelativeInteriorCompatible
    ],
    "support_1"
  )
  expect_identical(
    owner_2_axis$face_table$FaceId[
      owner_2_axis$face_table$FirstOrderBoundaryReachable
    ],
    "support_1"
  )
  expect_true(owner_2_axis$finite_response_image_excluded_by_full_support_lp)
})

test_that("first-order face lift gives an explicit finite parameter path", {
  skip_if_not_installed("lpSolve")

  target <- c(log(2), log(4), 0, 0)
  result <- mfrmr:::mfrmr_jml_gpcm_response_image_face_p2i_fixture(target)
  indices <- c(0, 1, 3, 8, 20)
  paths <- lapply(indices, function(index) {
    mfrmr:::mfrmr_jml_gpcm_response_image_face_path_at(
      result,
      "support_2",
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
    state <- mfrmr:::mfrmr_jml_gpcm_binary_closure_image_state(
      path$contrasts
    )
    isTRUE(state$finite_response_image)
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

test_that("three-owner rank-one operator separates outer faces from lifts", {
  skip_if_not_installed("lpSolve")

  design <- matrix(1, nrow = 3, ncol = 1)
  owner <- 1:3

  finite <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    owner,
    c(1, 2, 4)
  )
  expect_identical(
    finite$state,
    "finite_response_image_representative_certified"
  )
  expect_identical(finite$face_count, 7L)
  expect_identical(
    finite$face_table$SupportSize,
    c(3L, 2L, 2L, 2L, 1L, 1L, 1L)
  )
  expect_equal(
    finite$finite_representative$reconstructed_contrasts,
    c(1, 2, 4),
    tolerance = 1e-12
  )

  reachable <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    owner,
    c(1, 0, 0)
  )
  expect_identical(
    reachable$state,
    "reachable_missing_boundary_face_certified"
  )
  expect_false(reachable$finite_response_image_certified)
  expect_true(reachable$finite_response_image_excluded_by_full_support_lp)
  expect_true(reachable$nonnegative_simplex_face_witness_exists)
  expect_true(reachable$target_in_finite_response_image_closure_certified)
  expect_true(reachable$reachable_missing_finite_boundary_certified)
  expect_identical(
    reachable$face_table$FaceId[
      reachable$face_table$FirstOrderBoundaryReachable
    ],
    "support_2_3"
  )
  path <- mfrmr:::mfrmr_jml_gpcm_response_image_face_path_at(
    reachable,
    "support_2_3",
    20
  )
  expect_true(path$valid)
  expect_equal(path$contrasts, c(1, 0, 0), tolerance = 5e-09)

  outer_only <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    owner,
    c(1, -2, 0)
  )
  expect_identical(
    outer_only$state,
    "outer_boundary_face_candidate_lift_open"
  )
  expect_false(outer_only$finite_response_image_certified)
  expect_true(outer_only$finite_response_image_excluded_by_full_support_lp)
  expect_true(outer_only$nonnegative_simplex_face_witness_exists)
  expect_false(outer_only$target_in_finite_response_image_closure_certified)
  expect_false(outer_only$reachable_missing_finite_boundary_certified)
  expect_false(outer_only$target_outside_closure_outer_bound_certified)
  expect_true(outer_only$response_image_closure_membership_open)
  expect_false(outer_only$nonnegative_face_condition_sufficient_for_general_closure)
  expect_identical(
    outer_only$face_table$FaceId[
      outer_only$face_table$RelativeInteriorCompatible
    ],
    "support_3"
  )
  expect_identical(
    outer_only$face_table$LiftSolverStatus[
      outer_only$face_table$RelativeInteriorCompatible
    ],
    2L
  )

  no_face <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    owner,
    c(1, -2, 3)
  )
  expect_identical(
    no_face$state,
    "outside_response_image_closure_outer_bound_certified"
  )
  expect_false(no_face$nonnegative_simplex_face_witness_exists)
  expect_true(no_face$target_outside_closure_outer_bound_certified)
})

test_that("finite image reconstruction is invariant to rows and extra columns", {
  skip_if_not_installed("lpSolve")

  design <- rbind(
    c(1, 0, -1),
    c(0, 1, -1),
    c(1, 0, 1),
    c(0, 1, 1),
    c(1, -1, 0),
    c(-1, 1, 0)
  )
  owner <- c(1, 1, 2, 2, 3, 3)
  log_slopes <- c(log(2), -log(3), log(3 / 2))
  slopes <- exp(log_slopes)
  beta <- c(0.7, -0.2, 0.4)
  target <- slopes[owner] * as.numeric(design %*% beta)
  result <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    owner,
    target
  )
  expect_true(result$finite_response_image_certified)
  expect_equal(prod(slopes), 1, tolerance = 1e-14)
  expect_equal(
    result$finite_representative$reconstructed_contrasts,
    target,
    tolerance = 1e-11
  )

  permutation <- c(6, 2, 4, 1, 5, 3)
  permuted <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design[permutation, , drop = FALSE],
    owner[permutation],
    target[permutation]
  )
  expect_identical(permuted$state, result$state)
  expect_identical(
    permuted$face_table$RelativeInteriorCompatible,
    result$face_table$RelativeInteriorCompatible
  )
  expect_equal(
    permuted$finite_representative$reconstructed_contrasts,
    target[permutation],
    tolerance = 1e-11
  )

  redundant <- cbind(design, design[, 1])
  expanded <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    redundant,
    owner,
    target
  )
  expect_true(expanded$finite_response_image_certified)
  expect_equal(
    expanded$finite_representative$reconstructed_contrasts,
    target,
    tolerance = 1e-11
  )

  sparse <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    Matrix::Matrix(design, sparse = TRUE),
    owner,
    target
  )
  expect_identical(sparse$state, result$state)
  expect_equal(
    sparse$finite_representative$reconstructed_contrasts,
    target,
    tolerance = 1e-11
  )
})

test_that("one-owner chart reduces to the closed additive column space", {
  skip_if_not_installed("lpSolve")

  design <- rbind(c(1, 0), c(0, 1), c(1, 1))
  finite <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    rep(1, 3),
    c(1, 2, 3)
  )
  expect_identical(finite$face_count, 1L)
  expect_identical(finite$compatible_proper_face_count, 0L)
  expect_identical(finite$first_order_reachable_face_count, 0L)
  expect_true(finite$finite_response_image_certified)
  expect_equal(finite$finite_representative$slopes, 1, tolerance = 0)
  expect_equal(
    finite$finite_representative$reconstructed_contrasts,
    c(1, 2, 3),
    tolerance = 1e-12
  )

  outside <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    rep(1, 3),
    c(1, 2, 4)
  )
  expect_identical(
    outside$state,
    "outside_response_image_closure_outer_bound_certified"
  )
  expect_true(outside$finite_response_image_excluded_by_full_support_lp)
  expect_false(outside$nonnegative_simplex_face_witness_exists)
  expect_true(outside$target_outside_closure_outer_bound_certified)
})

test_that("near-face positive margins remain open until numerically certified", {
  skip_if_not_installed("lpSolve")

  design <- rbind(
    c(1, 0, -1),
    c(0, 1, -1),
    c(1, 0, 1),
    c(0, 1, 1)
  )
  target <- c(1, 0, 1e-06, 0)
  guarded <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    adjacent_design = design,
    slope_owner_by_row = c(1, 1, 2, 2),
    target_contrasts = target,
    certificate_tolerance = 1e-05
  )
  expect_identical(guarded$state, "not_certified_numerical_margin")
  expect_true(guarded$evaluated)
  expect_false(guarded$response_image_face_chart_classified)
  expect_false(guarded$finite_response_image_certified)
  expect_false(guarded$finite_response_image_excluded_by_full_support_lp)
  expect_false(guarded$reachable_missing_finite_boundary_certified)
  expect_false(guarded$target_outside_closure_outer_bound_certified)
  expect_gt(guarded$face_table$MaximumRelativeMargin[1], 0)
  expect_lt(guarded$face_table$MaximumRelativeMargin[1], 1e-05)

  resolved <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    adjacent_design = design,
    slope_owner_by_row = c(1, 1, 2, 2),
    target_contrasts = target,
    certificate_tolerance = 1e-08
  )
  expect_identical(
    resolved$state,
    "finite_response_image_representative_certified"
  )
  expect_true(resolved$finite_response_image_certified)
  expect_equal(
    resolved$finite_representative$reconstructed_contrasts,
    target,
    tolerance = 1e-12
  )
})

test_that("face chart and path controls fail closed", {
  skip_if_not_installed("lpSolve")

  design <- matrix(1, nrow = 3, ncol = 1)
  valid <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    1:3,
    c(1, 0, 0)
  )

  control_cases <- list(
    list(certificate_tolerance = 0),
    list(max_rows = 0),
    list(max_coordinates = 0),
    list(max_owners = 0),
    list(max_faces = 0),
    list(max_lp_nonzeros = 0),
    list(lp_timeout = 0)
  )
  for (controls in control_cases) {
    result <- do.call(
      mfrmr:::mfrmr_jml_gpcm_response_image_face_chart,
      c(list(design, 1:3, c(1, 0, 0)), controls)
    )
    expect_identical(result$state, "not_evaluated_control")
    expect_false(result$evaluated)
    expect_false(result$response_image_face_chart_classified)
    expect_false(result$finite_response_image_certified)
    expect_false(result$reachable_missing_finite_boundary_certified)
  }

  input_cases <- list(
    list(adjacent_design = matrix(numeric(0), 0, 0),
         slope_owner_by_row = 1:3, target_contrasts = c(1, 0, 0)),
    list(adjacent_design = design,
         slope_owner_by_row = c(1, 2), target_contrasts = c(1, 0, 0)),
    list(adjacent_design = design,
         slope_owner_by_row = c(1, 3, 3), target_contrasts = c(1, 0, 0)),
    list(adjacent_design = design,
         slope_owner_by_row = c(1, 2, 1e20), target_contrasts = c(1, 0, 0)),
    list(adjacent_design = design,
         slope_owner_by_row = 1:3, target_contrasts = c(1, 0)),
    list(adjacent_design = design,
         slope_owner_by_row = 1:3, target_contrasts = c(1, NA, 0)),
    list(adjacent_design = design,
         slope_owner_by_row = 1:3, target_contrasts = c(1, 0, 0),
         additive_offset = c(0, NA, 0))
  )
  for (arguments in input_cases) {
    result <- do.call(
      mfrmr:::mfrmr_jml_gpcm_response_image_face_chart,
      arguments
    )
    expect_identical(result$state, "not_evaluated_input")
    expect_false(result$evaluated)
    expect_false(result$response_image_face_chart_classified)
  }

  offset <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    1:3,
    c(1, 0, 0),
    additive_offset = c(0, 1, 0)
  )
  expect_identical(offset$state, "not_evaluated_nonzero_offset")
  expect_false(offset$evaluated)
  expect_false(offset$response_image_face_chart_classified)
  expect_false(offset$zero_affine_offset_certified)

  workload <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    1:3,
    c(1, 0, 0),
    max_faces = 6
  )
  expect_identical(workload$state, "not_evaluated_face_workload")
  expect_false(workload$evaluated)
  expect_false(workload$response_image_face_chart_classified)

  solver_cap <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    1:3,
    c(1, 0, 0),
    max_lp_nonzeros = 1
  )
  expect_identical(solver_cap$state, "not_evaluated_solver")
  expect_true(solver_cap$evaluated)
  expect_false(solver_cap$response_image_face_chart_classified)

  invalid_sources <- list(
    1,
    within(valid, contract_version <- "legacy"),
    valid
  )
  invalid_faces <- c("support_2_3", "support_2_3", "unknown")
  for (index in seq_along(invalid_sources)) {
    path <- mfrmr:::mfrmr_jml_gpcm_response_image_face_path_at(
      invalid_sources[[index]],
      invalid_faces[index],
      1
    )
    expect_false(path$valid)
    expect_identical(path$state, "not_evaluated_face_path")
  }
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_response_image_face_path_at(
      valid,
      "support_2_3",
      -1
    )$valid
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_response_image_face_path_at(
      valid,
      "support_2_3",
      301
    )$valid
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_response_image_face_path_at(
      valid,
      "support_2_3",
      1,
      max_index = 0
    )$valid
  )

  finite <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    1:3,
    c(1, 2, 4)
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_response_image_face_path_at(
      finite,
      "support_1_2_3",
      1
    )$valid
  )
  outer <- mfrmr:::mfrmr_jml_gpcm_response_image_face_chart(
    design,
    1:3,
    c(1, -2, 0)
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_response_image_face_path_at(
      outer,
      "support_3",
      1
    )$valid
  )
})
