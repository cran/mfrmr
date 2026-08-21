fit_jml_gpcm_sequence_remainder_fixture <- function() {
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

test_that("finite sequences receive a diagnostic but no asymptotic authority", {
  distance <- c(1, 4, 16, 64, 256)
  sequence <- cbind(
    primary = distance,
    secondary = sqrt(distance),
    unused = 0
  )
  result <- mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
    sequence,
    distances = distance,
    reference = c(0, 0, 0)
  )

  expect_identical(
    result$contract_version,
    "mfrmr-jml-gpcm-sequence-remainder-diagnostic-0.2.3-v1"
  )
  expect_identical(
    result$state, "finite_sequence_direction_scale_diagnostic"
  )
  expect_true(result$evaluated)
  expect_identical(result$sequence_points, 5L)
  expect_identical(result$free_coordinates, 3L)
  expect_identical(nrow(result$stage_table), 2L)
  expect_identical(nrow(result$direction_table), 6L)
  expect_equal(
    sum(result$stage_table$DisplacementEnergyFraction), 1,
    tolerance = 1e-12
  )
  expect_true(result$finite_low_rank_reconstruction_within_tolerance)
  expect_true(result$direction_estimates_are_coordinate_basis_dependent)
  expect_true(result$scale_exponents_estimated)
  expect_false(result$stage_endpoints_are_optimizer_iterations)
  expect_false(result$asymptotic_direction_certified)
  expect_false(result$scale_exponents_certified)
  expect_false(result$remainder_vanishing_certified)
  expect_false(result$exact_coefficient_comparisons_eligible)
  expect_false(result$p2b_handoff_eligible)
  expect_false(result$path_inferred_from_optimizer_trace)
  expect_false(result$path_search_performed)
  expect_false(result$global_boundary_classified)
  expect_false(result$external_comparison_eligible)
  expect_identical(result$readiness_effect, "none_diagnostic_only")
})

test_that("a rank-one power sequence reports only a finite scale estimate", {
  distance <- c(1, 2, 4, 8, 16)
  direction <- c(1, -2, 0.5)
  sequence <- outer(distance^1.5, direction)
  result <- mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
    sequence,
    distances = distance,
    reference = rep(0, 3),
    max_stages = 1L
  )

  expect_true(result$evaluated)
  expect_identical(nrow(result$stage_table), 1L)
  expect_equal(
    result$stage_table$EstimatedPowerExponent, 1.5,
    tolerance = 1e-12
  )
  expect_equal(result$stage_table$PowerFitRSquared, 1, tolerance = 1e-12)
  expect_lt(max(result$residual_ratio), 1e-12)
  expect_false(result$scale_exponents_certified)
  expect_false(result$p2b_handoff_eligible)
})

test_that("finite sequence diagnostics fail closed on invalid inputs", {
  short <- mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
    matrix(1:4, nrow = 2L)
  )
  expect_identical(short$state, "not_evaluated_sequence")

  sequence <- matrix(1:9, nrow = 3L)
  distances <- mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
    sequence, distances = c(1, 3, 2)
  )
  expect_identical(distances$state, "not_evaluated_distances")

  reference <- mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
    sequence, reference = c(0, 0)
  )
  expect_identical(reference$state, "not_evaluated_reference")

  stationary <- mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
    matrix(1, nrow = 3L, ncol = 2L)
  )
  expect_identical(stationary$state, "not_evaluated_no_displacement")

  bad_map <- mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
    sequence,
    coordinate_map = data.frame(
      OptimizerIndex = c(1L, 1L, 3L),
      Coordinate = letters[1:3]
    )
  )
  expect_identical(bad_map$state, "not_evaluated_coordinate_map")

  fractional_map <- mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
    sequence,
    coordinate_map = data.frame(
      OptimizerIndex = c(1, 2, 3.5),
      Coordinate = letters[1:3]
    )
  )
  expect_identical(fractional_map$state, "not_evaluated_coordinate_map")

  too_large <- mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
    sequence, max_elements = 8L
  )
  expect_identical(too_large$state, "not_evaluated_size_limit")

  bad_control <- mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
    sequence, max_stages = 0L
  )
  expect_identical(bad_control$state, "not_evaluated_control")
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_finite_sequence_diagnostic(
      sequence, max_stages = 3L
    )$state,
    "not_evaluated_control"
  )
})

test_that("fit audit consumes transient stage endpoints without retaining them", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_sequence_remainder_fixture()
  audit <- fit$config$boundary_audit$gpcm_optimizer_sequence_diagnostic

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-sequence-remainder-diagnostic-0.2.3-v1"
  )
  expect_gte(audit$stage_endpoint_rows, 1L)
  if (audit$stage_endpoint_rows < 3L) {
    expect_identical(
      audit$state, "not_evaluated_insufficient_stage_endpoints"
    )
    expect_false(audit$evaluated)
    expect_false(audit$finite_direction_scale_diagnostic_completed)
  } else {
    expect_identical(
      audit$state,
      "finite_stage_endpoint_direction_scale_diagnostic_only"
    )
    expect_true(audit$evaluated)
    expect_true(audit$finite_direction_scale_diagnostic_completed)
  }
  expect_gt(audit$free_coordinates, 0L)
  expect_true(audit$aggregate_stage_history_available)
  expect_true(audit$stage_endpoint_coordinates_available_during_fit)
  expect_false(audit$stage_endpoint_coordinates_retained_in_fit)
  expect_false(audit$stage_endpoints_are_optimizer_iterations)
  expect_true(audit$source_parameter_path_contract_matches)
  expect_false(audit$production_path_limit_classified)
  expect_false(audit$asymptotic_direction_certified)
  expect_false(audit$scale_exponents_certified)
  expect_false(audit$remainder_vanishing_certified)
  expect_false(audit$p2b_handoff_eligible)
  expect_false(audit$path_inferred_from_optimizer_trace)
  expect_false(audit$path_search_performed)
  expect_false(audit$global_boundary_classified)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
  expect_null(fit$opt$optimizer_stage_parameters)
  expect_null(fit$config$optimizer_stage_parameters)
})

test_that("three explicit stage endpoints stay diagnostic at fit scope", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_sequence_remainder_fixture()
  sizes <- mfrmr:::build_param_sizes(fit$config)
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  sequence <- matrix(0, nrow = 4L, ncol = free_dimension)
  sequence[, 1L] <- c(0, 1, 4, 9)
  stages <- data.frame(Stage = 1:4)

  result <- mfrmr:::audit_mfrm_jml_gpcm_optimizer_sequence_scope(
    config = fit$config,
    sizes = sizes,
    stages = stages,
    parameter_sequence = sequence,
    parameter_path =
      fit$config$boundary_audit$gpcm_parameter_path_reachability
  )

  expect_identical(
    result$state,
    "finite_stage_endpoint_direction_scale_diagnostic_only"
  )
  expect_true(result$evaluated)
  expect_true(result$finite_direction_scale_diagnostic_completed)
  expect_true(result$finite_sequence_diagnostic$evaluated)
  expect_false(result$asymptotic_direction_certified)
  expect_false(result$scale_exponents_certified)
  expect_false(result$remainder_vanishing_certified)
  expect_false(result$p2b_handoff_eligible)
  expect_false(result$production_path_limit_classified)
})

test_that("negative-power scaled-logit remainders retain the P2c limit", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_sequence_remainder_fixture()
  path <- mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(fit)
  expect_true(path$parameter_reachable_limit_contract_completed)

  base_utilities <- path$limit$path_definition$base_utilities
  n_observations <- nrow(base_utilities)
  n_categories <- ncol(base_utilities)
  primary <- matrix(0, n_observations, n_categories)
  primary[, 1L] <- 1
  primary[, 2L] <- -1
  secondary <- matrix(0, n_observations, n_categories)
  secondary[, n_categories] <- 0.5

  result <- mfrmr:::mfrmr_jml_gpcm_decaying_logit_remainder(
    path,
    residual_stage_directions = list(primary, secondary),
    residual_decay_exponents = c(0.5, 1.5)
  )

  expect_identical(
    result$state, "decaying_scaled_logit_remainder_same_limit"
  )
  expect_true(result$evaluated)
  expect_identical(result$remainder_stages, 2L)
  expect_equal(
    result$log_likelihood_limit, path$log_likelihood_limit,
    tolerance = 0
  )
  expect_true(result$within_row_logit_contrast_remainder_vanishes)
  expect_true(result$same_row_probability_limits_certified)
  expect_true(result$same_joint_log_likelihood_limit_certified)
  expect_true(result$p2b_limit_extended_only_within_declared_remainder_class)
  expect_false(result$finite_distance_ties_required)
  expect_false(result$arbitrary_utility_remainder_classified)
  expect_false(result$arbitrary_slope_remainder_classified)
  expect_false(result$arbitrary_logit_remainder_classified)
  expect_false(result$optimizer_sequence_remainder_certified)
  expect_false(result$global_boundary_classified)
  expect_false(result$external_comparison_eligible)
  expect_identical(result$readiness_effect, "none_diagnostic_only")

  distances <- c(1, 4, 16, 64, 256)
  differences <- vapply(distances, function(distance) {
    declared <- mfrmr:::mfrmr_jml_gpcm_parameter_hierarchy_loglik_at(
      path, distance
    )
    remainder <-
      mfrmr:::mfrmr_jml_gpcm_parameter_remainder_loglik_at(
        path, result, distance
      )
    expect_true(declared$valid)
    expect_true(remainder$valid)
    max(abs(
      declared$row_log_probability - remainder$row_log_probability
    ))
  }, numeric(1))
  expect_true(all(diff(differences) < 0))
  expect_lt(tail(differences, 1L), head(differences, 1L) / 10)
})

test_that("remainder theorem inputs and scopes fail closed", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_sequence_remainder_fixture()
  path <- mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(fit)
  dimensions <- dim(path$limit$path_definition$base_utilities)
  direction <- matrix(0, dimensions[1], dimensions[2])

  invalid_source <- path
  invalid_source$contract_version <- "legacy"
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_decaying_logit_remainder(
      invalid_source, list(direction), 1
    )$state,
    "not_evaluated_source_path"
  )
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_decaying_logit_remainder(
      path, list(direction), 0
    )$state,
    "not_evaluated_remainder_mapping"
  )
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_decaying_logit_remainder(
      path, list(direction, direction), c(1, 0.5)
    )$state,
    "not_evaluated_remainder_mapping"
  )
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_decaying_logit_remainder(
      path, list(matrix(0, 1L, 1L)), 1
    )$state,
    "not_evaluated_remainder_mapping"
  )
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_decaying_logit_remainder(
      path, list(direction), 1, max_stages = -1L
    )$state,
    "not_evaluated_control"
  )
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_decaying_logit_remainder(
      path, list(direction), 1, max_stages = 3L
    )$state,
    "not_evaluated_control"
  )

  invalid_evaluation <-
    mfrmr:::mfrmr_jml_gpcm_parameter_remainder_loglik_at(
      path,
      mfrmr:::mfrmr_jml_gpcm_decaying_logit_remainder(
        path, list(direction), 1
      ),
      0
    )
  expect_false(invalid_evaluation$valid)
})

test_that("optimizer sequence scope does not borrow other authorities", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_sequence_remainder_fixture()
  sizes <- mfrmr:::build_param_sizes(fit$config)
  dimension <- sum(vapply(sizes, as.integer, integer(1)))
  sequence <- matrix(0, nrow = 3L, ncol = dimension)
  stages <- data.frame(Stage = 1:3)
  path <- fit$config$boundary_audit$gpcm_parameter_path_reachability

  pcm <- fit$config
  pcm$model <- "PCM"
  expect_identical(
    mfrmr:::audit_mfrm_jml_gpcm_optimizer_sequence_scope(
      pcm, sizes, stages, sequence, path
    )$state,
    "not_applicable_model"
  )

  mml <- fit$config
  mml$method <- "MML"
  expect_identical(
    mfrmr:::audit_mfrm_jml_gpcm_optimizer_sequence_scope(
      mml, sizes, stages, sequence, path
    )$state,
    "not_applicable_estimator"
  )

  unit <- fit$config
  unit$gpcm_spec$n_params <- 0L
  expect_identical(
    mfrmr:::audit_mfrm_jml_gpcm_optimizer_sequence_scope(
      unit, sizes, stages, sequence, path
    )$state,
    "not_required_unit_slope"
  )

  legacy <- path
  legacy$contract_version <- "legacy"
  expect_identical(
    mfrmr:::audit_mfrm_jml_gpcm_optimizer_sequence_scope(
      fit$config, sizes, stages, sequence, legacy
    )$state,
    "not_evaluated_source_contract"
  )
})
