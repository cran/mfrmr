p2l_binary_design <- function() {
  rbind(
    c(1, 0, -1),
    c(0, 1, -1),
    c(1, 0, 1),
    c(0, 1, 1)
  )
}

p2l_binary_owner <- function() {
  c(1, 1, 2, 2)
}

test_that("saturated ordered-category likelihood helpers are exact", {
  adjacent <- c(log(2), log(3), log(4), log(5))
  cumulative <- mfrmr:::mfrmr_jml_gpcm_saturated_envelope_cumulative(
    adjacent,
    observations = 2,
    categories = 3
  )
  expect_equal(
    cumulative,
    rbind(c(0, log(2), log(8)), c(0, log(3), log(15))),
    tolerance = 1e-15
  )

  mass <- rbind(c(1, 2, 8), c(1, 3, 15))
  expected <- sum(mass * log(mass / rowSums(mass)))
  expect_equal(
    mfrmr:::mfrmr_jml_gpcm_saturated_envelope_loglik(adjacent, mass),
    expected,
    tolerance = 1e-14
  )
  expect_equal(
    mfrmr:::mfrmr_jml_gpcm_saturated_envelope_logsumexp(c(1000, 999)),
    1000 + log1p(exp(-1)),
    tolerance = 0
  )
})

test_that("positive binary masses can have a nonattained saturated JML bound", {
  skip_if_not_installed("lpSolve")

  result <-
    mfrmr:::mfrmr_jml_gpcm_saturated_envelope_p2i_binary_fixture()
  expected <- c(log(2), log(4), 0, 0)
  mass <- cbind(rep(1, 4), c(2, 4, 1, 1))
  expected_value <- sum(mass * log(mass / rowSums(mass)))

  expect_identical(
    result$contract_version,
    "mfrmr-jml-gpcm-saturated-response-envelope-0.2.3-v1"
  )
  expect_identical(
    result$state,
    "global_supremum_nonattained_at_saturated_boundary_certified"
  )
  expect_true(result$evaluated)
  expect_true(result$saturated_response_envelope_classified)
  expect_identical(
    result$estimator_identity,
    "unpenalized_fixed_effects_jml_no_finite_box"
  )
  expect_identical(
    result$objective_identity,
    "identified_conditional_joint_log_likelihood"
  )
  expect_match(result$response_likelihood_identity, "ordered-category")
  expect_match(result$response_contrast_layout, "transition-major")
  expect_true(result$strictly_positive_category_mass)
  expect_true(result$ordinary_observation_category_separation_absent)
  expect_true(result$independent_saturated_optimum_unique)
  expect_true(
    result$positive_mass_response_log_likelihood_strictly_concave
  )
  expect_true(
    result$unbounded_response_contrast_likelihood_limit_negative_infinity
  )
  expect_true(
    result$every_finite_limsup_parameter_sequence_has_bounded_response_subsequence
  )
  expect_identical(result$observations, 4L)
  expect_identical(result$categories, 2L)
  expect_identical(result$adjacent_response_contrast_dimension, 4L)
  expect_identical(result$additive_free_coordinates, 3L)
  expect_identical(result$slope_owners, 2L)
  expect_equal(
    result$independent_saturated_optimum_adjacent_contrasts,
    expected,
    tolerance = 0
  )
  expect_equal(
    result$independent_saturated_optimum_cumulative_contrasts[, 2],
    expected,
    tolerance = 0
  )
  expect_equal(
    result$independent_saturated_optimum_log_likelihood,
    expected_value,
    tolerance = 1e-14
  )
  expect_identical(
    result$global_supremum_value,
    result$saturated_log_likelihood_upper_bound
  )
  expect_true(result$face_chart$response_image_face_chart_classified)
  expect_true(result$higher_order_closure$higher_order_face_lifts_classified)
  expect_false(result$saturated_target_finite_response_image_certified)
  expect_true(result$saturated_target_missing_finite_boundary_certified)
  expect_false(
    result$saturated_target_outside_response_image_closure_certified
  )
  expect_true(result$saturated_target_closure_membership_classified)
  expect_true(
    result$complete_instance_parameter_sequence_limsup_envelope_constructed
  )
  expect_true(
    result$global_supremum_exists_as_response_closure_maximum_certified
  )
  expect_true(result$global_supremum_equals_saturated_upper_bound_certified)
  expect_false(
    result$global_supremum_strictly_below_saturated_upper_bound_certified
  )
  expect_true(result$global_supremum_value_classified)
  expect_false(result$response_image_supremum_attained_by_finite_parameter)
  expect_false(result$finite_jmle_existence_certified)
  expect_true(result$finite_jmle_nonexistence_certified)
  expect_true(
    result$finite_jmle_nonexistence_without_observation_category_separation_certified
  )
  expect_true(result$saturated_boundary_competitive_certified)
  expect_false(result$global_constrained_optimum_location_and_value_open)
  expect_identical(result$supremum_path_source, "p2k_higher_order_hierarchy")
  expect_identical(result$supremum_path_id, "rates_1_0")
  expect_length(result$finite_global_representative, 0L)
  expect_false(result$optimizer_trace_used)
  expect_false(result$optimizer_trace_promoted_to_finite_estimate)
  expect_false(result$standard_error_eligible)
  expect_false(result$confidence_interval_eligible)
  expect_false(result$external_comparison_eligible)
  expect_identical(result$readiness_effect, "none_theorem_only")
  expect_identical(
    result$next_gate,
    "apply_saturated_envelope_to_exact_current_fit_operators"
  )
})

test_that("certified finite parameter paths attain the saturated supremum", {
  skip_if_not_installed("lpSolve")

  result <-
    mfrmr:::mfrmr_jml_gpcm_saturated_envelope_p2i_binary_fixture()
  indices <- c(0, 1, 3, 8, 20)
  paths <- lapply(indices, function(index) {
    mfrmr:::mfrmr_jml_gpcm_saturated_envelope_path_at(result, index)
  })
  expect_true(all(vapply(paths, `[[`, logical(1), "valid")))
  expect_true(all(vapply(
    paths, `[[`, logical(1), "finite_parameter_point"
  )))
  expect_true(all(vapply(
    paths, `[[`, logical(1), "parameter_path_declared_divergent"
  )))
  expect_true(all(vapply(
    paths, `[[`, logical(1), "response_contrast_limit_certified"
  )))
  expect_true(all(vapply(paths, function(path) {
    identical(path$source_path_type, "p2k_higher_order_hierarchy") &&
      identical(path$source_path_id, "rates_1_0")
  }, logical(1))))

  distances <- vapply(
    paths, `[[`, numeric(1), "contrast_distance_to_target"
  )
  gaps <- vapply(
    paths, `[[`, numeric(1), "target_minus_path_log_likelihood"
  )
  likelihoods <- vapply(paths, `[[`, numeric(1), "log_likelihood")
  expect_true(all(diff(distances) < 0))
  expect_true(all(diff(gaps) < 0))
  expect_true(all(diff(likelihoods) > 0))
  expect_true(all(gaps >= -1e-12))
  expect_lt(tail(distances, 1), 4e-09)
  expect_lt(tail(gaps, 1), 1e-12)
  expect_equal(
    tail(paths, 1)[[1]]$contrasts,
    result$independent_saturated_optimum_adjacent_contrasts,
    tolerance = 4e-09
  )
  expect_true(all(vapply(paths, function(path) {
    abs(sum(path$source_path$log_slopes)) <= 1e-12 &&
      abs(prod(path$source_path$inverse_slopes) - 1) <= 1e-12
  }, logical(1))))
})

test_that("finite saturated targets certify finite global JMLEs", {
  skip_if_not_installed("lpSolve")

  design <- p2l_binary_design()
  owner <- p2l_binary_owner()
  finite <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    design,
    owner,
    cbind(rep(1, 4), c(4, 2, 4, 2))
  )
  target <- log(c(4, 2, 4, 2))

  expect_identical(
    finite$state,
    "finite_global_jmle_at_saturated_response_optimum_certified"
  )
  expect_true(finite$saturated_response_envelope_classified)
  expect_true(finite$saturated_target_finite_response_image_certified)
  expect_false(finite$saturated_target_missing_finite_boundary_certified)
  expect_false(
    finite$saturated_target_outside_response_image_closure_certified
  )
  expect_true(finite$global_supremum_value_classified)
  expect_true(finite$response_image_supremum_attained_by_finite_parameter)
  expect_true(finite$finite_jmle_existence_certified)
  expect_false(finite$finite_jmle_nonexistence_certified)
  expect_false(finite$saturated_boundary_competitive_certified)
  expect_identical(finite$supremum_path_source, "none")
  expect_identical(finite$supremum_path_id, "")
  expect_true(is.list(finite$finite_global_representative))
  expect_equal(
    finite$finite_global_representative$reconstructed_contrasts,
    target,
    tolerance = 2e-12
  )
  expect_lte(
    finite$finite_global_representative$maximum_absolute_reconstruction_residual,
    2e-12
  )
  expect_equal(
    prod(finite$finite_global_representative$inverse_slopes),
    1,
    tolerance = 1e-14
  )

  intersection <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    design,
    owner,
    cbind(rep(1, 4), rep(2, 4))
  )
  expect_true(intersection$finite_jmle_existence_certified)
  expect_equal(
    intersection$independent_saturated_optimum_adjacent_contrasts,
    rep(log(2), 4),
    tolerance = 0
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_saturated_envelope_path_at(finite, 1)$valid
  )
})

test_that("outside saturated targets retain a strict but uncomputed bound", {
  skip_if_not_installed("lpSolve")

  target <- c(1, -2, 0, 1)
  result <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    p2l_binary_design(),
    p2l_binary_owner(),
    cbind(rep(1, 4), exp(target))
  )

  expect_identical(
    result$state,
    "saturated_upper_bound_outside_closure_constrained_envelope_open"
  )
  expect_true(result$saturated_response_envelope_classified)
  expect_false(result$saturated_target_finite_response_image_certified)
  expect_false(result$saturated_target_missing_finite_boundary_certified)
  expect_true(
    result$saturated_target_outside_response_image_closure_certified
  )
  expect_true(result$saturated_target_closure_membership_classified)
  expect_false(
    result$complete_instance_parameter_sequence_limsup_envelope_constructed
  )
  expect_true(
    result$global_supremum_exists_as_response_closure_maximum_certified
  )
  expect_false(result$global_supremum_equals_saturated_upper_bound_certified)
  expect_true(
    result$global_supremum_strictly_below_saturated_upper_bound_certified
  )
  expect_false(result$global_supremum_value_classified)
  expect_true(is.na(result$global_supremum_value))
  expect_true(is.finite(result$saturated_log_likelihood_upper_bound))
  expect_false(result$finite_jmle_existence_certified)
  expect_false(result$finite_jmle_nonexistence_certified)
  expect_true(result$global_constrained_optimum_location_and_value_open)
  expect_false(result$constrained_boundary_likelihood_envelope_constructed)
  expect_identical(result$supremum_path_source, "none")
  expect_identical(
    result$next_gate,
    "constrained_closure_likelihood_envelope_below_saturated_bound"
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_saturated_envelope_path_at(result, 1)$valid
  )
})

test_that("the theorem extends from binary to ordered polytomous mass", {
  skip_if_not_installed("lpSolve")

  binary_design <- p2l_binary_design()
  zero <- matrix(0, nrow(binary_design), ncol(binary_design))
  design <- rbind(
    cbind(binary_design, zero),
    cbind(zero, binary_design)
  )
  first_transition <- c(log(2), log(4), 0, 0)
  mass <- cbind(1, exp(first_transition), exp(first_transition))
  result <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    design,
    rep(p2l_binary_owner(), 2),
    mass
  )

  expect_identical(result$categories, 3L)
  expect_identical(result$adjacent_response_contrast_dimension, 8L)
  expect_equal(
    result$independent_saturated_optimum_adjacent_contrasts,
    c(first_transition, rep(0, 4)),
    tolerance = 0
  )
  expect_identical(
    result$state,
    "global_supremum_nonattained_at_saturated_boundary_certified"
  )
  expect_true(result$finite_jmle_nonexistence_certified)
  expect_true(result$global_supremum_value_classified)
  expect_equal(
    result$global_supremum_value,
    sum(mass * log(mass / rowSums(mass))),
    tolerance = 2e-14
  )

  path <- mfrmr:::mfrmr_jml_gpcm_saturated_envelope_path_at(result, 8)
  expect_true(path$valid)
  expect_lt(path$contrast_distance_to_target, 6e-04)
  expect_lt(path$target_minus_path_log_likelihood, 2e-07)

  finite_mass <- rbind(c(1, 2, 6), c(2, 3, 12))
  finite <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    diag(4),
    c(1, 2, 1, 2),
    finite_mass
  )
  expect_true(finite$finite_jmle_existence_certified)
  expect_equal(
    finite$independent_saturated_optimum_adjacent_contrasts,
    c(log(2), log(3 / 2), log(3), log(4)),
    tolerance = 1e-15
  )
})

test_that("positive mass supplies strict concavity and contrast coercivity", {
  mass <- cbind(rep(1, 4), c(2, 4, 1, 1))
  optimum <- log(mass[, 2] / mass[, 1])
  perturbation <- c(0.3, -0.4, 0.2, -0.1)
  loglik <- mfrmr:::mfrmr_jml_gpcm_saturated_envelope_loglik

  optimum_value <- loglik(optimum, mass)
  plus <- loglik(optimum + perturbation, mass)
  minus <- loglik(optimum - perturbation, mass)
  expect_gt(optimum_value, plus)
  expect_gt(optimum_value, minus)
  expect_gt(optimum_value, (plus + minus) / 2)

  direction <- c(1, -1, 1, -1)
  values <- vapply(c(1, 10, 100), function(scale) {
    loglik(optimum + scale * direction, mass)
  }, numeric(1))
  expect_true(all(diff(values) < 0))
  expect_lt(tail(values, 1), -300)
})

test_that("classification respects mass scaling and aligned row permutations", {
  skip_if_not_installed("lpSolve")

  design <- p2l_binary_design()
  owner <- p2l_binary_owner()
  mass <- cbind(rep(1, 4), c(2, 4, 1, 1))
  original <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    design, owner, mass
  )
  scaled <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    design, owner, 7 * mass
  )
  expect_identical(scaled$state, original$state)
  expect_equal(
    scaled$independent_saturated_optimum_adjacent_contrasts,
    original$independent_saturated_optimum_adjacent_contrasts,
    tolerance = 1e-15
  )
  expect_equal(
    scaled$global_supremum_value,
    7 * original$global_supremum_value,
    tolerance = 2e-14
  )

  permutation <- c(3, 1, 4, 2)
  permuted <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    design[permutation, , drop = FALSE],
    owner[permutation],
    mass[permutation, , drop = FALSE]
  )
  expect_identical(permuted$state, original$state)
  expect_equal(
    permuted$independent_saturated_optimum_adjacent_contrasts,
    original$independent_saturated_optimum_adjacent_contrasts[permutation],
    tolerance = 0
  )
  expect_equal(
    permuted$global_supremum_value,
    original$global_supremum_value,
    tolerance = 1e-14
  )

  sparse <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    Matrix::Matrix(design, sparse = TRUE), owner, mass
  )
  expanded <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    cbind(design, 0), owner, mass
  )
  expect_identical(sparse$state, original$state)
  expect_identical(expanded$state, original$state)
  expect_equal(sparse$global_supremum_value, original$global_supremum_value)
  expect_equal(expanded$global_supremum_value, original$global_supremum_value)
})

test_that("malformed saturated-envelope inputs fail closed", {
  skip_if_not_installed("lpSolve")

  design <- p2l_binary_design()
  owner <- p2l_binary_owner()
  mass <- cbind(rep(1, 4), c(2, 4, 1, 1))
  evaluate <- function(...) {
    mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
      design, owner, mass, ...
    )
  }

  invalid_controls <- list(
    list(certificate_tolerance = 0),
    list(certificate_tolerance = 1),
    list(max_observations = 1.5),
    list(max_categories = 0),
    list(max_mass_elements = Inf),
    list(max_coordinates = 0),
    list(max_owners = 1.5),
    list(max_faces = 0),
    list(max_stages = 0),
    list(max_hierarchies = 0),
    list(max_equations = 0),
    list(max_lp_nonzeros = 0),
    list(lp_timeout = 0)
  )
  for (control in invalid_controls) {
    result <- do.call(evaluate, control)
    expect_identical(result$state, "not_evaluated_control")
    expect_false(result$evaluated)
    expect_false(result$saturated_response_envelope_classified)
    expect_false(result$finite_jmle_existence_certified)
    expect_false(result$finite_jmle_nonexistence_certified)
  }

  invalid_inputs <- list(
    list(category_mass = cbind(c(1, 1, 0, 1), c(2, 4, 1, 1))),
    list(category_mass = cbind(c(1, 1, NA, 1), c(2, 4, 1, 1))),
    list(category_mass = matrix(letters[1:8], 4, 2)),
    list(category_mass = matrix(1, 4, 1)),
    list(category_mass = matrix(1, 3, 2)),
    list(slope_owner_by_row = c(1, 1, 3, 3)),
    list(slope_owner_by_row = c(1, 1, 2)),
    list(slope_owner_by_row = c(1, 1, 2, 0)),
    list(slope_owner_by_row = rep(1e20, 4)),
    list(additive_offset = c(0, 0, 0, 1e-06)),
    list(additive_offset = c(0, 0, 0)),
    list(adjacent_design = matrix(NA_real_, 4, 3)),
    list(adjacent_design = matrix(1, 3, 2)),
    list(max_observations = 3),
    list(max_categories = 1),
    list(max_mass_elements = 7),
    list(max_total_mass = 1),
    list(max_coordinates = 2),
    list(max_owners = 1)
  )
  for (replacement in invalid_inputs) {
    arguments <- list(
      adjacent_design = design,
      slope_owner_by_row = owner,
      category_mass = mass
    )
    arguments[names(replacement)] <- replacement
    result <- expect_silent(do.call(
      mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope,
      arguments
    ))
    expect_identical(result$state, "not_evaluated_input")
    expect_false(result$saturated_response_envelope_classified)
    expect_false(result$finite_jmle_existence_certified)
    expect_false(result$finite_jmle_nonexistence_certified)
  }
})

test_that("numerical and workload ambiguity never promotes an envelope", {
  skip_if_not_installed("lpSolve")

  design <- p2l_binary_design()
  owner <- p2l_binary_owner()
  near <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    design,
    owner,
    cbind(rep(1, 4), exp(c(1, 0, 1e-06, 0))),
    certificate_tolerance = 1e-05
  )
  expect_identical(near$state, "not_evaluated_face_chart")
  expect_true(near$evaluated)
  expect_false(near$saturated_response_envelope_classified)
  expect_true(near$strictly_positive_category_mass)
  expect_true(near$independent_saturated_optimum_unique)
  expect_false(near$finite_jmle_existence_certified)
  expect_false(near$finite_jmle_nonexistence_certified)

  face_cap <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    design,
    owner,
    cbind(rep(1, 4), c(2, 4, 1, 1)),
    max_faces = 1
  )
  expect_identical(face_cap$state, "not_evaluated_face_chart")
  expect_false(face_cap$saturated_response_envelope_classified)

  hierarchy_cap <- mfrmr:::mfrmr_jml_gpcm_saturated_response_envelope(
    matrix(1, nrow = 3, ncol = 1),
    1:3,
    cbind(rep(1, 3), exp(c(1, -2, 0))),
    max_hierarchies = 2
  )
  expect_identical(
    hierarchy_cap$state,
    "not_evaluated_closure_membership"
  )
  expect_true(hierarchy_cap$evaluated)
  expect_false(hierarchy_cap$saturated_response_envelope_classified)
  expect_false(hierarchy_cap$global_supremum_value_classified)
  expect_false(hierarchy_cap$finite_jmle_existence_certified)
  expect_false(hierarchy_cap$finite_jmle_nonexistence_certified)
  expect_false(hierarchy_cap$optimizer_trace_promoted_to_finite_estimate)
})

test_that("saturated-envelope path controls fail closed", {
  skip_if_not_installed("lpSolve")

  result <-
    mfrmr:::mfrmr_jml_gpcm_saturated_envelope_p2i_binary_fixture()
  invalid <- list(
    list(result, -1),
    list(result, 1.5),
    list(result, NA_real_),
    list(result, 101, max_index = 100),
    list(result, 1, max_index = -1),
    list(result, 1, max_log_rate_span = 0),
    list(list(), 1)
  )
  for (arguments in invalid) {
    path <- do.call(
      mfrmr:::mfrmr_jml_gpcm_saturated_envelope_path_at,
      arguments
    )
    expect_false(path$valid)
    expect_identical(path$state, "not_evaluated_saturated_envelope_path")
  }
})
