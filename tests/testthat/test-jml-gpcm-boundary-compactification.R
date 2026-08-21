fit_jml_gpcm_compactification_fixture <- function() {
  data <- expand.grid(
    Person = c("P1", "P2"),
    Criterion = c("C1", "C2", "C3"),
    Score = 0:1,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data <- data[rep(seq_len(nrow(data)), each = 2L), , drop = FALSE]
  suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = "Criterion",
    score = "Score",
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    rating_min = 0,
    rating_max = 1,
    maxit = 200,
    reltol = 1e-12
  ))
}

test_that("primary P/Z/L/D role counts close the finite combinatorics", {
  counts <- do.call(rbind, lapply(
    2:6, mfrmr:::mfrmr_jml_gpcm_primary_role_counts
  ))
  expect_equal(
    counts$CanonicalPrimaryRolePatterns,
    c(2, 18, 110, 570, 2702)
  )
  expect_equal(
    counts$AdditionalRolePatterns,
    c(0, 12, 98, 550, 2672)
  )
  expect_equal(
    counts$ZeroFreePrimaryRolePatterns,
    c(2, 12, 50, 180, 602)
  )
  expect_equal(
    counts$ZeroContainingPrimaryRolePatterns,
    c(0, 6, 60, 390, 2100)
  )
  expect_equal(
    counts$CanonicalPrimaryRolePatterns,
    vapply(
      2:6,
      mfrmr:::mfrmr_jml_gpcm_joint_rate_candidate_count,
      numeric(1)
    )
  )
  expect_equal(
    counts$OrderedPairRolePatterns,
    (2:6) * (1:5)
  )

  empty <- mfrmr:::mfrmr_jml_gpcm_primary_role_counts(1L)
  expect_identical(empty$CanonicalPrimaryRolePatterns, 0)
  expect_identical(empty$AdditionalRolePatterns, 0)
})

test_that("nonzero sum-zero rates receive canonical primary roles", {
  partition <- mfrmr:::mfrmr_jml_gpcm_rate_role_partition(c(3, -1, -2))
  expect_identical(
    partition$state, "nonzero_sum_zero_primary_rate_partitioned"
  )
  expect_true(partition$valid_direction)
  expect_true(partition$sum_zero)
  expect_identical(partition$roles, c("P", "L", "D"))
  expect_identical(partition$role_signature, "PLD")
  expect_identical(partition$positive_indices, 1L)
  expect_identical(partition$leading_negative_indices, 2L)
  expect_identical(partition$deeper_negative_indices, 3L)
  expect_length(partition$zero_indices, 0L)
  expect_equal(partition$canonical_equivalent_rates, c(3, -1, -2))
  expect_true(partition$canonical_equivalent_sum_zero)
  expect_false(partition$secondary_hierarchy_required)

  with_zero <-
    mfrmr:::mfrmr_jml_gpcm_rate_role_partition(c(1, -1, 0))
  expect_identical(with_zero$roles, c("P", "L", "Z"))
  expect_true(with_zero$secondary_hierarchy_required)
  expect_identical(with_zero$zero_indices, 3L)

  permutation <- c(3L, 1L, 2L)
  permuted <- mfrmr:::mfrmr_jml_gpcm_rate_role_partition(
    c(3, -1, -2)[permutation]
  )
  expect_identical(permuted$roles, partition$roles[permutation])
  expect_equal(
    permuted$canonical_equivalent_rates,
    partition$canonical_equivalent_rates[permutation]
  )
  expect_true(permuted$canonical_equivalent_sum_zero)

  zero <- mfrmr:::mfrmr_jml_gpcm_rate_role_partition(c(0, 0, 0))
  expect_identical(zero$state, "all_zero_primary_slope_rate")
  expect_false(zero$valid_direction)
  expect_true(zero$all_zero_primary_slope_rate)
  expect_identical(zero$roles, rep("Z", 3L))
  expect_true(zero$secondary_hierarchy_required)

  invalid <- mfrmr:::mfrmr_jml_gpcm_rate_role_partition(c(1, 0, -0.5))
  expect_identical(invalid$state, "invalid_rate_vector")
  expect_false(invalid$sum_zero)
  expect_false(invalid$valid_direction)
})

test_that("combined directions separate slope-primary and additive-primary sectors", {
  joint <- mfrmr:::mfrmr_jml_gpcm_compactify_direction(
    additive_displacement = c(2, -1),
    expanded_log_slope_displacement = c(3, -1, -2)
  )
  expect_identical(
    joint$state, "nonzero_slope_primary_direction_partitioned"
  )
  expect_true(joint$valid_direction)
  expect_equal(joint$normalized_norm, 1, tolerance = 1e-14)
  expect_false(joint$slope_primary_zero)
  expect_identical(joint$rate_partition$roles, c("P", "L", "D"))
  expect_false(joint$secondary_hierarchy_required)

  additive <- mfrmr:::mfrmr_jml_gpcm_compactify_direction(
    additive_displacement = c(1, -1),
    expanded_log_slope_displacement = c(0, 0, 0)
  )
  expect_identical(
    additive$state, "additive_primary_zero_slope_rate"
  )
  expect_true(additive$valid_direction)
  expect_true(additive$slope_primary_zero)
  expect_true(additive$secondary_hierarchy_required)
  expect_identical(additive$rate_partition$roles, rep("Z", 3L))

  malformed <- mfrmr:::mfrmr_jml_gpcm_compactify_direction(
    1, c(1, 0, 0)
  )
  expect_identical(malformed$state, "invalid_direction")
  expect_false(malformed$valid_direction)
})

test_that("nonconvergent normalized rates reduce to distinct subsequences", {
  rate_a <- c(1, -1, 0)
  rate_b <- c(0, 1, -1)
  indices <- seq_len(12L)
  directions <- t(vapply(indices, function(index) {
    rate <- if (index %% 2L == 1L) rate_a else rate_b
    compact <- mfrmr:::mfrmr_jml_gpcm_compactify_direction(
      numeric(0), index * rate
    )
    compact$normalized_expanded_log_slopes
  }, numeric(3)))

  expect_gt(min(sqrt(rowSums(diff(directions)^2))), 0.5)
  expect_equal(
    directions[indices %% 2L == 1L, , drop = FALSE],
    matrix(rate_a / sqrt(2), nrow = 6L, ncol = 3L, byrow = TRUE),
    tolerance = 1e-14
  )
  expect_equal(
    directions[indices %% 2L == 0L, , drop = FALSE],
    matrix(rate_b / sqrt(2), nrow = 6L, ncol = 3L, byrow = TRUE),
    tolerance = 1e-14
  )
  expect_false(isTRUE(all(diff(directions) == 0)))
})

test_that("a zero primary rate can carry a divergent secondary hierarchy", {
  primary <- c(1, -1, 0)
  secondary <- c(-0.5, -0.5, 1)
  distances <- c(1e2, 1e4, 1e6)
  paths <- t(vapply(distances, function(distance) {
    distance * primary + sqrt(distance) * secondary
  }, numeric(3)))
  primary_scaled <- paths / distances
  residuals <- paths - outer(distances, primary)
  residual_norms <- sqrt(rowSums(residuals^2))
  residual_directions <- residuals / residual_norms
  classified <-
    mfrmr:::mfrmr_jml_gpcm_rate_role_partition(primary)

  expect_equal(rowSums(paths), rep(0, length(distances)), tolerance = 1e-10)
  expect_lt(
    max(abs(primary_scaled[nrow(primary_scaled), ] - primary)),
    max(abs(primary_scaled[1, ] - primary))
  )
  expect_gt(min(diff(residual_norms)), 0)
  expect_equal(
    residual_directions,
    matrix(
      secondary / sqrt(sum(secondary^2)),
      nrow = length(distances), ncol = 3L, byrow = TRUE
    ),
    tolerance = 1e-12
  )
  expect_identical(classified$roles, c("P", "L", "Z"))
  expect_true(classified$secondary_hierarchy_required)
})

test_that("production fits retain diagnostic-only compactification", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_compactification_fixture()
  audit <- fit$config$boundary_audit$gpcm_boundary_compactification

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-boundary-compactification-0.2.3-v1"
  )
  expect_identical(
    audit$state,
    "primary_subsequence_partitioned_secondary_hierarchy_open"
  )
  expect_true(audit$evaluated)
  expect_true(audit$finite_parameter_dimension)
  expect_gt(audit$combined_free_dimension, 0L)
  expect_gt(audit$additive_free_coordinates, 0L)
  expect_identical(audit$slope_levels, 3L)
  expect_identical(audit$free_log_slope_coordinates, 2L)
  expect_true(audit$unit_sphere_compactness_applied)
  expect_true(audit$unbounded_sequence_has_convergent_normalized_subsequence)
  expect_false(
    audit$whole_sequence_normalized_direction_convergence_required
  )
  expect_true(audit$accumulation_direction_reduction_available)
  expect_true(audit$canonical_primary_role_partition_complete)
  expect_equal(audit$role_counts$CanonicalPrimaryRolePatterns, 18)
  expect_equal(audit$role_counts$ZeroContainingPrimaryRolePatterns, 6)
  expect_equal(audit$role_counts$AdditionalRolePatterns, 12)
  expect_true(audit$all_zero_slope_primary_sector_present)
  expect_true(audit$zero_primary_rate_secondary_hierarchy_required)
  expect_false(audit$secondary_hierarchy_classified)
  expect_false(audit$finite_secondary_hierarchy_depth_certified)
  expect_true(audit$canonical_role_boundary_search_complete)
  expect_true(audit$source_boundary_contract_matches)
  expect_true(audit$affine_transport_contract_matches)
  expect_false(audit$nonconvergent_normalized_rate_boundary_classified)
  expect_false(audit$arbitrary_nonvanishing_residual_boundary_classified)
  expect_false(audit$general_curved_path_classified)
  expect_false(audit$global_boundary_classified)
  expect_false(audit$global_finite_maximum_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
  expect_false(audit$criterion_frozen_for_scientific_claims)
  expect_match(audit$limitations, "slower divergent scales", fixed = TRUE)
})

test_that("scope controls fail closed without borrowing MML or FACETS claims", {
  fit <- fit_jml_gpcm_compactification_fixture()
  config <- fit$config
  sizes <- mfrmr:::build_param_sizes(config)
  classified <-
    config$boundary_audit$gpcm_fixed_objective_classification
  transport <-
    config$boundary_audit$gpcm_asymptotically_affine_transport

  mismatched <- classified
  mismatched$contract_version <- "mismatched"
  structural <-
    mfrmr:::audit_mfrm_jml_gpcm_boundary_compactification(
      config, sizes, mismatched, transport
    )
  expect_identical(
    structural$state,
    "primary_subsequence_partitioned_secondary_hierarchy_open"
  )
  expect_false(structural$source_boundary_contract_matches)
  expect_false(structural$global_boundary_classified)

  mml_config <- config
  mml_config$method <- "MML"
  mml <- mfrmr:::audit_mfrm_jml_gpcm_boundary_compactification(
    mml_config, sizes, classified, transport
  )
  expect_identical(mml$state, "not_applicable_estimator")
  expect_identical(mml$objective_identity, "not_applicable_marginal_estimator")
  expect_false(mml$external_comparison_eligible)

  pcm_config <- config
  pcm_config$model <- "PCM"
  pcm <- mfrmr:::audit_mfrm_jml_gpcm_boundary_compactification(
    pcm_config, sizes, classified, transport
  )
  expect_identical(pcm$state, "not_applicable_model")
  expect_identical(pcm$objective_identity, "not_applicable")

  unit_config <- config
  unit_config$gpcm_spec$levels <- "C1"
  unit_config$gpcm_spec$n_params <- 0L
  unit_sizes <- sizes
  unit_sizes$log_slopes <- 0L
  unit <- mfrmr:::audit_mfrm_jml_gpcm_boundary_compactification(
    unit_config, unit_sizes, classified, transport
  )
  expect_identical(unit$state, "not_required_unit_slope")
  expect_true(unit$unit_sphere_compactness_applied)
  expect_false(unit$canonical_primary_role_partition_complete)

  invalid_dimension <- sizes
  invalid_dimension$log_slopes <- 99L
  invalid <- mfrmr:::audit_mfrm_jml_gpcm_boundary_compactification(
    config, invalid_dimension, classified, transport
  )
  expect_identical(invalid$state, "not_evaluated_dimension")
  expect_false(invalid$evaluated)

  invalid_control <-
    mfrmr:::audit_mfrm_jml_gpcm_boundary_compactification(
      config, sizes, classified, transport, zero_tolerance = -1
    )
  expect_identical(invalid_control$state, "not_evaluated_control")
  expect_false(invalid_control$evaluated)
})
