fit_jml_gpcm_rate_hierarchy_fixture <- function() {
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

test_that("a two-stage hierarchy resolves a primary zero-rate coordinate", {
  rates <- rbind(
    c(1, -1, 0),
    c(-0.5, -0.5, 1)
  )
  hierarchy <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    rates,
    levels = c("C1", "C2", "C3"),
    scale_exponents = c(1, 0.5)
  )

  expect_identical(
    hierarchy$contract_version,
    "mfrmr-jml-gpcm-rate-hierarchy-0.2.3-v1"
  )
  expect_identical(hierarchy$state, "complete_finite_rate_hierarchy")
  expect_true(hierarchy$evaluated)
  expect_true(hierarchy$complete)
  expect_identical(hierarchy$slope_levels, 3L)
  expect_identical(hierarchy$declared_stages, 2L)
  expect_identical(hierarchy$used_stages, 2L)
  expect_identical(hierarchy$nonzero_stages, 2L)
  expect_identical(hierarchy$maximum_nonzero_stages, 2L)
  expect_true(hierarchy$finite_depth_bound_respected)
  expect_true(hierarchy$scale_order_declared)
  expect_true(hierarchy$scale_order_valid)
  expect_true(hierarchy$global_sum_zero_each_stage)
  expect_false(hierarchy$active_restriction_sum_zero_required)
  expect_identical(hierarchy$remaining_unresolved_count, 0L)
  expect_length(hierarchy$remaining_unresolved_levels, 0L)
  expect_true(hierarchy$structural_growth_hierarchy_classified)
  expect_false(hierarchy$likelihood_boundary_classified)
  expect_false(hierarchy$common_likelihood_limit_certified)
  expect_false(hierarchy$additive_hierarchy_classified)
  expect_false(hierarchy$global_boundary_classified)
  expect_false(hierarchy$criterion_frozen_for_scientific_claims)

  stages <- hierarchy$stage_table
  expect_identical(stages$ActiveCount, c(3L, 1L))
  expect_identical(stages$ResolvedCount, c(2L, 1L))
  expect_identical(stages$RemainingZeroCount, c(1L, 0L))
  expect_identical(stages$RoleSignature, c("PLZ", "..P"))
  expect_identical(stages$PositiveLevels, c("C1", "C3"))
  expect_identical(stages$LeadingNegativeLevels, c("C2", ""))
  expect_identical(stages$ZeroLevels, c("C3", ""))
  expect_equal(stages$GlobalSumZeroResidual, c(0, 0))
  expect_true(all(stages$GlobalSumZero))
  expect_identical(
    stages$StageState,
    c("resolved_subset", "resolved_all_remaining_levels")
  )

  resolution <- hierarchy$level_resolution
  expect_identical(resolution$Level, c("C1", "C2", "C3"))
  expect_identical(resolution$FirstNonzeroStage, c(1L, 1L, 2L))
  expect_identical(resolution$RoleAtResolution, c("P", "L", "P"))
  expect_equal(resolution$NormalizedLoadingAtResolution, c(1, -1, 1))
  expect_true(all(resolution$Resolved))
})

test_that("the sharp finite coordinate depth is J minus one", {
  rates <- rbind(
    c(1, -1, 0, 0, 0),
    c(-1, 0, 1, 0, 0),
    c(-1, 0, 0, 1, 0),
    c(-1, 0, 0, 0, 1)
  )
  levels <- paste0("C", 1:5)
  hierarchy <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    rates, levels, scale_exponents = c(1, 0.75, 0.5, 0.25)
  )

  expect_identical(hierarchy$state, "complete_finite_rate_hierarchy")
  expect_true(hierarchy$complete)
  expect_identical(hierarchy$nonzero_stages, 4L)
  expect_identical(hierarchy$maximum_nonzero_stages, 4L)
  expect_true(hierarchy$finite_depth_bound_respected)
  expect_identical(
    hierarchy$stage_table$ActiveCount, c(5L, 3L, 2L, 1L)
  )
  expect_identical(
    hierarchy$stage_table$ResolvedCount, c(2L, 1L, 1L, 1L)
  )
  expect_identical(
    hierarchy$stage_table$RoleSignature,
    c("PLZZZ", "..PZZ", "...PZ", "....P")
  )
  expect_identical(
    hierarchy$level_resolution$FirstNonzeroStage,
    c(1L, 1L, 2L, 3L, 4L)
  )

  permutation <- c(5L, 3L, 1L, 4L, 2L)
  permuted <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    rates[, permutation, drop = FALSE],
    levels[permutation],
    scale_exponents = c(1, 0.75, 0.5, 0.25)
  )
  original_resolution <- hierarchy$level_resolution[
    order(hierarchy$level_resolution$Level), , drop = FALSE
  ]
  permuted_resolution <- permuted$level_resolution[
    order(permuted$level_resolution$Level), , drop = FALSE
  ]
  expect_identical(permuted$state, hierarchy$state)
  expect_identical(
    permuted_resolution$FirstNonzeroStage,
    original_resolution$FirstNonzeroStage
  )
  expect_identical(
    permuted_resolution$RoleAtResolution,
    original_resolution$RoleAtResolution
  )
  expect_equal(
    permuted_resolution$NormalizedLoadingAtResolution,
    original_resolution$NormalizedLoadingAtResolution
  )
})

test_that("unresolved and redundant declared stages remain explicit", {
  exhausted <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    matrix(c(1, -1, 0), nrow = 1L),
    levels = c("C1", "C2", "C3")
  )
  expect_identical(exhausted$state, "declared_depth_exhausted")
  expect_true(exhausted$evaluated)
  expect_false(exhausted$complete)
  expect_identical(exhausted$remaining_unresolved_count, 1L)
  expect_identical(exhausted$remaining_unresolved_levels, "C3")
  expect_false(exhausted$scale_order_declared)
  expect_true(exhausted$scale_order_valid)
  expect_false(exhausted$structural_growth_hierarchy_classified)
  expect_match(exhausted$reason_codes, "additional_slower_scale_required")

  stalled <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    rbind(c(1, -1, 0, 0), c(1, -1, 0, 0)),
    levels = paste0("C", 1:4),
    scale_exponents = c(1, 0.5)
  )
  expect_identical(stalled$state, "incomplete_zero_active_loading")
  expect_false(stalled$complete)
  expect_identical(stalled$used_stages, 2L)
  expect_identical(stalled$nonzero_stages, 1L)
  expect_identical(stalled$remaining_unresolved_count, 2L)
  expect_identical(stalled$remaining_unresolved_levels, c("C3", "C4"))
  expect_identical(
    stalled$stage_table$StageState[2], "stalled_zero_active_loading"
  )
  expect_identical(stalled$stage_table$ResolvedCount[2], 0L)

  extra <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    rbind(
      c(1, -1, 0),
      c(-0.5, -0.5, 1),
      c(1, -1, 0)
    ),
    levels = c("C1", "C2", "C3"),
    scale_exponents = c(1, 0.5, 0.25)
  )
  expect_true(extra$complete)
  expect_identical(extra$used_stages, 2L)
  expect_identical(extra$unused_declared_stages, 1L)
  expect_identical(nrow(extra$stage_table), 2L)
})

test_that("hierarchy controls and sum-zero identification fail closed", {
  nonfinite <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    matrix(c(1, -1, NA_real_), nrow = 1L)
  )
  expect_identical(nonfinite$state, "not_evaluated_rate_matrix")
  expect_false(nonfinite$evaluated)

  one_level <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(1)
  expect_identical(one_level$state, "not_evaluated_rate_matrix")

  duplicate <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    c(1, -1, 0), levels = c("C1", "C1", "C3")
  )
  expect_identical(duplicate$state, "not_evaluated_levels")
  expect_false(duplicate$evaluated)

  scale <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    rbind(c(1, -1, 0), c(-0.5, -0.5, 1)),
    scale_exponents = c(1, 1)
  )
  expect_identical(scale$state, "not_evaluated_scale_order")
  expect_false(scale$evaluated)
  expect_true(scale$scale_order_declared)
  expect_false(scale$scale_order_valid)

  unidentified <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(c(1, -0.5, 0))
  expect_identical(unidentified$state, "not_evaluated_identification")
  expect_false(unidentified$evaluated)
  expect_false(unidentified$global_sum_zero_each_stage)

  control <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    c(1, -1, 0), zero_tolerance = -1
  )
  expect_identical(control$state, "not_evaluated_control")
  expect_false(control$evaluated)
})

test_that("secondary slope roles can change the limit at one primary direction", {
  primary <- c(1, -1, 0)
  secondary_high <- c(-0.5, -0.5, 1)
  secondary_low <- -secondary_high
  distances <- c(25, 100, 400)
  top_log_probability <- function(log_slope) {
    slope <- exp(log_slope)
    -log1p(exp(-slope))
  }
  high <- vapply(distances, function(distance) {
    ell <- distance * primary + sqrt(distance) * secondary_high
    top_log_probability(ell[3])
  }, numeric(1))
  low <- vapply(distances, function(distance) {
    ell <- distance * primary + sqrt(distance) * secondary_low
    top_log_probability(ell[3])
  }, numeric(1))
  high_primary <- t(vapply(distances, function(distance) {
    (distance * primary + sqrt(distance) * secondary_high) / distance
  }, numeric(3)))
  low_primary <- t(vapply(distances, function(distance) {
    (distance * primary + sqrt(distance) * secondary_low) / distance
  }, numeric(3)))
  high_hierarchy <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    rbind(primary, secondary_high), scale_exponents = c(1, 0.5)
  )
  low_hierarchy <- mfrmr:::mfrmr_jml_gpcm_rate_hierarchy(
    rbind(primary, secondary_low), scale_exponents = c(1, 0.5)
  )

  expect_equal(rowSums(high_primary), rep(0, 3L), tolerance = 1e-14)
  expect_equal(rowSums(low_primary), rep(0, 3L), tolerance = 1e-14)
  expect_lt(
    max(abs(high_primary[3, ] - primary)),
    max(abs(high_primary[1, ] - primary))
  )
  expect_lt(
    max(abs(low_primary[3, ] - primary)),
    max(abs(low_primary[1, ] - primary))
  )
  expect_lt(abs(high[3]), 1e-12)
  expect_lt(abs(low[3] + log(2)), 1e-8)
  expect_gt(abs(high[3] - low[3]), 0.6)
  expect_identical(
    high_hierarchy$stage_table$RoleSignature, c("PLZ", "..P")
  )
  expect_identical(
    low_hierarchy$stage_table$RoleSignature, c("PLZ", "..L")
  )
  expect_true(high_hierarchy$complete)
  expect_true(low_hierarchy$complete)
  expect_false(high_hierarchy$common_likelihood_limit_certified)
  expect_false(low_hierarchy$common_likelihood_limit_certified)
})

test_that("production scope records a finite depth but no likelihood claim", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_rate_hierarchy_fixture()
  audit <- fit$config$boundary_audit$gpcm_rate_hierarchy

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-rate-hierarchy-0.2.3-v1"
  )
  expect_identical(
    audit$state,
    "finite_coordinate_rate_hierarchy_available_likelihood_open"
  )
  expect_true(audit$evaluated)
  expect_identical(audit$slope_levels, 3L)
  expect_identical(audit$free_log_slope_coordinates, 2L)
  expect_true(audit$finite_coordinate_rate_hierarchy_theorem)
  expect_identical(audit$maximum_nonzero_slope_hierarchy_stages, 2L)
  expect_false(audit$later_active_restriction_sum_zero_required)
  expect_true(audit$global_coefficient_sum_zero_required)
  expect_true(audit$pathwise_hierarchy_extraction_available)
  expect_false(audit$arbitrary_fit_path_hierarchy_computed)
  expect_false(audit$arbitrary_sequence_likelihood_classified)
  expect_false(
    audit$common_accumulation_direction_likelihood_limit_certified
  )
  expect_false(audit$additive_hierarchy_classified)
  expect_true(audit$secondary_hierarchy_structurally_bounded)
  expect_false(audit$secondary_hierarchy_likelihood_classified)
  expect_false(audit$global_boundary_classified)
  expect_false(audit$global_finite_maximum_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
  expect_true(audit$source_compactification_contract_matches)
  expect_match(audit$limitations, "Additive directions", fixed = TRUE)
})

test_that("production hierarchy scope does not borrow MML or FACETS authority", {
  fit <- fit_jml_gpcm_rate_hierarchy_fixture()
  config <- fit$config
  sizes <- mfrmr:::build_param_sizes(config)
  compact <- config$boundary_audit$gpcm_boundary_compactification

  mismatch <- compact
  mismatch$contract_version <- "mismatched"
  structural <- mfrmr:::audit_mfrm_jml_gpcm_rate_hierarchy_scope(
    config, sizes, mismatch
  )
  expect_identical(
    structural$state,
    "finite_coordinate_rate_hierarchy_available_likelihood_open"
  )
  expect_false(structural$source_compactification_contract_matches)
  expect_false(structural$global_boundary_classified)

  mml_config <- config
  mml_config$method <- "MML"
  mml <- mfrmr:::audit_mfrm_jml_gpcm_rate_hierarchy_scope(
    mml_config, sizes, compact
  )
  expect_identical(mml$state, "not_applicable_estimator")
  expect_identical(mml$objective_identity, "not_applicable_marginal_estimator")
  expect_false(mml$external_comparison_eligible)

  pcm_config <- config
  pcm_config$model <- "PCM"
  pcm <- mfrmr:::audit_mfrm_jml_gpcm_rate_hierarchy_scope(
    pcm_config, sizes, compact
  )
  expect_identical(pcm$state, "not_applicable_model")
  expect_identical(pcm$objective_identity, "not_applicable")

  unit_config <- config
  unit_config$gpcm_spec$levels <- "C1"
  unit_config$gpcm_spec$n_params <- 0L
  unit_sizes <- sizes
  unit_sizes$log_slopes <- 0L
  unit <- mfrmr:::audit_mfrm_jml_gpcm_rate_hierarchy_scope(
    unit_config, unit_sizes, compact
  )
  expect_identical(unit$state, "not_required_unit_slope")
  expect_false(unit$finite_coordinate_rate_hierarchy_theorem)

  invalid_sizes <- sizes
  invalid_sizes$log_slopes <- 99L
  invalid <- mfrmr:::audit_mfrm_jml_gpcm_rate_hierarchy_scope(
    config, invalid_sizes, compact
  )
  expect_identical(invalid$state, "not_evaluated_dimension")
  expect_false(invalid$evaluated)
})
