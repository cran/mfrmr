fit_jml_gpcm_exponential_balance_fixture <- function(anchors = NULL) {
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
    anchors = anchors,
    rating_min = 0,
    rating_max = 1,
    maxit = 180,
    reltol = 1e-11
  ))
}

exponential_balance_problem <- function() {
  list(
    adjacent_design = diag(2),
    additive_center = c(0, 1),
    additive_coordinate_directions = list(c(2, 0)),
    additive_exponential_rates = -1,
    score_k = c(1L, 0L),
    slope_index = c(1L, 2L),
    base_log_slopes = c(0, 0),
    slope_log_rates = c(1, -1),
    weight = c(1, 2)
  )
}

call_exponential_balance <- function(problem, ...) {
  do.call(
    mfrmr:::mfrmr_jml_gpcm_exponential_balance_limit,
    c(problem, list(...))
  )
}

test_that("inverse exponential utility decay gives a bounded response-image escape", {
  problem <- exponential_balance_problem()
  result <- call_exponential_balance(problem)
  expected <- 2 - mfrmr:::logsumexp(c(0, 2)) - 2 * log(2)

  expect_identical(
    result$contract_version,
    "mfrmr-jml-gpcm-exponential-balance-0.2.3-v1"
  )
  expect_identical(
    result$state, "bounded_response_image_parameter_escape_certified"
  )
  expect_true(result$evaluated)
  expect_true(result$exponential_balance_limit_classified)
  expect_true(result$parameter_space_reachability_checked)
  expect_true(result$parameter_space_reachability_certified)
  expect_true(result$exact_slope_utility_product_expanded)
  expect_true(result$equal_combined_exponents_aggregated_before_classification)
  expect_true(result$exact_combined_exponent_comparisons)
  expect_false(result$base_log_slopes_sum_zero_structurally_certified)
  expect_identical(result$base_log_slope_sum_zero_residual, 0)
  expect_identical(result$slope_log_rate_sum_zero_residual, 0)
  expect_identical(result$observations, 2L)
  expect_identical(result$categories, 2L)
  expect_identical(result$slope_levels, 2L)
  expect_identical(result$additive_free_coordinates, 2L)
  expect_identical(result$additive_exponential_components, 1L)
  expect_identical(result$combined_exponent_groups, 2L)
  expect_identical(result$positive_contrast_exponent_groups, 0L)
  expect_true(result$zero_contrast_exponent_present)
  expect_identical(result$negative_contrast_exponent_groups, 1L)
  expect_true(result$parameter_path_divergent)
  expect_true(result$divergent_log_slope_path)
  expect_false(result$divergent_additive_path)
  expect_true(result$path_escapes_every_bounded_parameter_set)
  expect_true(result$response_contrast_image_bounded)
  expect_true(result$bounded_response_image_escape_certified)
  expect_true(result$parameter_to_response_contrast_map_properness_refuted)
  expect_false(result$properness_globally_certified)
  expect_false(result$complete_bounded_image_escape_family_classified)
  expect_false(result$complete_escaping_sequence_boundary_envelope_constructed)
  expect_equal(result$log_likelihood_limit, expected, tolerance = 1e-14)
  expect_true(result$finite_log_likelihood_limit)
  expect_false(result$negative_infinite_log_likelihood_limit)
  expect_identical(
    result$rate_table$CombinedExponent, c(0, -1)
  )
  expect_identical(
    result$rate_table$AsymptoticRole,
    c("finite_contrast_base", "vanishing_contrast_remainder")
  )
  expect_equal(
    result$contrast_flag_limit$base_contrasts,
    matrix(c(2, 0), nrow = 2L), tolerance = 0
  )
  expect_identical(
    result$contrast_flag_limit$divergent_flag_stages, 0L
  )
  expect_false(result$path_inferred_from_optimizer_trace)
  expect_false(result$exponent_rates_inferred)
  expect_false(result$path_search_performed)
  expect_false(result$competitive_boundary_certified)
  expect_false(result$global_supremum_value_classified)
  expect_false(result$finite_jmle_existence_certified)
  expect_false(result$finite_jmle_nonexistence_certified)
  expect_false(result$global_boundary_absence_certified)
  expect_false(result$standard_error_eligible)
  expect_false(result$confidence_interval_eligible)
  expect_false(result$external_comparison_eligible)
  expect_identical(result$readiness_effect, "none_diagnostic_only")
})

test_that("finite-index evaluation follows the exact original parameter path", {
  result <- call_exponential_balance(exponential_balance_problem())
  for (index in c(0, 0.5, 2, 5)) {
    direct <- mfrmr:::mfrmr_jml_gpcm_exponential_balance_loglik_at(
      result, index
    )
    logits <- c(2, exp(-index))
    manual_rows <- c(
      logits[1] - mfrmr:::logsumexp(c(0, logits[1])),
      -mfrmr:::logsumexp(c(0, logits[2]))
    )
    expect_true(direct$valid)
    expect_equal(direct$row_log_probability, manual_rows, tolerance = 1e-14)
    expect_equal(
      direct$log_likelihood,
      manual_rows[1] + 2 * manual_rows[2], tolerance = 1e-14
    )
    expect_equal(direct$maximum_absolute_log_slope, index, tolerance = 0)
  }
  tail_value <- mfrmr:::mfrmr_jml_gpcm_exponential_balance_loglik_at(
    result, 20
  )
  expect_true(tail_value$valid)
  expect_equal(
    tail_value$log_likelihood,
    result$log_likelihood_limit, tolerance = 1e-8
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_exponential_balance_loglik_at(result, -1)$valid
  )
  expect_false(
    mfrmr:::mfrmr_jml_gpcm_exponential_balance_loglik_at(list(), 1)$valid
  )
})

test_that("equal combined exponents are aggregated before asymptotic roles", {
  problem <- list(
    adjacent_design = matrix(c(1, 1), nrow = 1L),
    additive_center = c(0, 0),
    additive_coordinate_directions = list(c(1, 0), c(0, -1)),
    additive_exponential_rates = c(-1, -1),
    score_k = 1L,
    slope_index = 1L,
    base_log_slopes = c(0, 0),
    slope_log_rates = c(1, -1),
    weight = 1
  )
  result <- call_exponential_balance(problem)

  expect_identical(
    result$state, "bounded_response_image_parameter_escape_certified"
  )
  expect_true(result$parameter_path_divergent)
  expect_true(result$response_contrast_image_bounded)
  expect_true(result$bounded_response_image_escape_certified)
  expect_identical(result$combined_exponent_groups, 0L)
  expect_identical(nrow(result$rate_table), 0L)
  expect_false(result$zero_contrast_exponent_present)
  expect_equal(
    result$contrast_flag_limit$base_contrasts,
    matrix(0, 1L, 1L), tolerance = 0
  )
  expect_equal(result$log_likelihood_limit, -log(2), tolerance = 1e-14)
  direct <- mfrmr:::mfrmr_jml_gpcm_exponential_balance_loglik_at(
    result, 4
  )
  expect_true(direct$valid)
  expect_equal(direct$log_likelihood, -log(2), tolerance = 1e-14)
})

test_that("positive combined rates become an exact P2e contrast flag", {
  problem <- exponential_balance_problem()
  problem$additive_center <- c(1, 1)
  problem$additive_coordinate_directions <- list()
  problem$additive_exponential_rates <- numeric(0)
  result <- call_exponential_balance(problem)

  expect_identical(
    result$state, "divergent_parameter_and_response_flag_limit_classified"
  )
  expect_true(result$parameter_path_divergent)
  expect_true(result$divergent_log_slope_path)
  expect_false(result$divergent_additive_path)
  expect_false(result$response_contrast_image_bounded)
  expect_false(result$bounded_response_image_escape_certified)
  expect_false(result$parameter_to_response_contrast_map_properness_refuted)
  expect_identical(result$positive_contrast_exponent_groups, 1L)
  expect_identical(
    result$rate_table$CombinedExponent, c(1, -1)
  )
  expect_identical(
    result$rate_table$AsymptoticRole,
    c("divergent_contrast_flag", "vanishing_contrast_remainder")
  )
  expect_identical(
    result$contrast_flag_limit$divergent_flag_stages, 1L
  )
  expect_true(result$finite_log_likelihood_limit)
  expect_equal(result$log_likelihood_limit, -2 * log(2), tolerance = 1e-14)
})

test_that("a nonzero affine offset participates in the exact exponent balance", {
  problem <- exponential_balance_problem()
  problem$additive_offset <- c(1, 0)
  result <- call_exponential_balance(problem)

  expect_identical(
    result$state, "divergent_parameter_and_response_flag_limit_classified"
  )
  expect_identical(result$rate_table$CombinedExponent, c(1, 0, -1))
  expect_identical(result$positive_contrast_exponent_groups, 1L)
  expect_false(result$response_contrast_image_bounded)
  expect_false(result$bounded_response_image_escape_certified)
  expect_false(result$parameter_to_response_contrast_map_properness_refuted)
  expect_equal(result$log_likelihood_limit, -2 * log(2), tolerance = 1e-14)

  direct <- mfrmr:::mfrmr_jml_gpcm_exponential_balance_loglik_at(result, 2)
  manual <- c(
    exp(2) + 2 - mfrmr:::logsumexp(c(0, exp(2) + 2)),
    -mfrmr:::logsumexp(c(0, exp(-2)))
  )
  expect_true(direct$valid)
  expect_equal(direct$row_log_probability, manual, tolerance = 1e-14)
  expect_equal(direct$log_likelihood, manual[1] + 2 * manual[2],
               tolerance = 1e-14)
})

test_that("divergent additive coordinates can balance a vanishing slope", {
  problem <- list(
    adjacent_design = matrix(1, 1L, 1L),
    additive_center = 0,
    additive_coordinate_directions = list(1),
    additive_exponential_rates = 1,
    score_k = 1L,
    slope_index = 1L,
    base_log_slopes = c(0, 0),
    slope_log_rates = c(-1, 1),
    weight = 1
  )
  result <- call_exponential_balance(problem)

  expect_identical(
    result$state, "bounded_response_image_parameter_escape_certified"
  )
  expect_true(result$parameter_path_divergent)
  expect_true(result$divergent_log_slope_path)
  expect_true(result$divergent_additive_path)
  expect_true(result$response_contrast_image_bounded)
  expect_true(result$bounded_response_image_escape_certified)
  expect_identical(result$rate_table$CombinedExponent, 0)
  expect_equal(result$log_likelihood_limit, 1 - log1p(exp(1)), tolerance = 1e-14)
  finite <- mfrmr:::mfrmr_jml_gpcm_exponential_balance_loglik_at(
    result, 5
  )
  expect_true(finite$valid)
  expect_equal(finite$log_likelihood, result$log_likelihood_limit,
               tolerance = 1e-14)
  expect_gt(finite$maximum_absolute_additive_coordinate, 100)
})

test_that("nondivergent exponential paths are classified without a boundary claim", {
  problem <- exponential_balance_problem()
  problem$slope_log_rates <- c(0, 0)
  problem$additive_exponential_rates <- -1
  result <- call_exponential_balance(problem)

  expect_identical(
    result$state, "bounded_parameter_exponential_path_limit_classified"
  )
  expect_false(result$parameter_path_divergent)
  expect_false(result$divergent_log_slope_path)
  expect_false(result$divergent_additive_path)
  expect_false(result$path_escapes_every_bounded_parameter_set)
  expect_true(result$response_contrast_image_bounded)
  expect_false(result$bounded_response_image_escape_certified)
  expect_false(result$parameter_to_response_contrast_map_properness_refuted)
})

test_that("malformed exponential balance contracts fail closed", {
  problem <- exponential_balance_problem()

  bad_design <- problem
  bad_design$adjacent_design <- matrix(1, 3L, 1L)
  expect_identical(
    call_exponential_balance(bad_design)$state, "not_evaluated_design"
  )

  bad_center <- problem
  bad_center$additive_center <- 0
  expect_identical(
    call_exponential_balance(bad_center)$state,
    "not_evaluated_additive_components"
  )

  bad_offset <- problem
  bad_offset$additive_offset <- 0
  expect_identical(
    call_exponential_balance(bad_offset)$state,
    "not_evaluated_additive_components"
  )

  nonfinite_offset <- problem
  nonfinite_offset$additive_offset <- c(0, NA_real_)
  expect_identical(
    call_exponential_balance(nonfinite_offset)$state,
    "not_evaluated_additive_components"
  )

  bad_direction <- problem
  bad_direction$additive_coordinate_directions <- list(c(1, NA))
  expect_identical(
    call_exponential_balance(bad_direction)$state,
    "not_evaluated_additive_components"
  )

  bad_additive_rate <- problem
  bad_additive_rate$additive_exponential_rates <- numeric(0)
  expect_identical(
    call_exponential_balance(bad_additive_rate)$state,
    "not_evaluated_additive_components"
  )

  zero_direction <- problem
  zero_direction$additive_coordinate_directions <- list(c(0, 0))
  expect_identical(
    call_exponential_balance(zero_direction)$state,
    "not_evaluated_parameter_reachability"
  )

  bad_base_slope <- problem
  bad_base_slope$base_log_slopes <- c(1, 0)
  expect_identical(
    call_exponential_balance(bad_base_slope)$state,
    "not_evaluated_slope_or_observation_map"
  )

  bad_slope_rate <- problem
  bad_slope_rate$slope_log_rates <- c(1, 0)
  expect_identical(
    call_exponential_balance(bad_slope_rate)$state,
    "not_evaluated_slope_or_observation_map"
  )

  bad_score <- problem
  bad_score$score_k <- c(2, 0)
  expect_identical(
    call_exponential_balance(bad_score)$state,
    "not_evaluated_slope_or_observation_map"
  )

  bad_weight <- problem
  bad_weight$weight <- c(0, 0)
  expect_identical(
    call_exponential_balance(bad_weight)$state, "not_evaluated_weight"
  )

  expect_identical(
    call_exponential_balance(problem, max_components = 0L)$state,
    "not_evaluated_additive_components"
  )
  expect_identical(
    call_exponential_balance(problem, max_design_nonzeros = 1L)$state,
    "not_evaluated_design"
  )
  expect_identical(
    call_exponential_balance(problem, max_exponent_groups = 1L)$state,
    "not_evaluated_exponent_group_limit"
  )
  expect_identical(
    call_exponential_balance(problem, max_contrast_elements = 1L)$state,
    "not_evaluated_contrast_element_limit"
  )
  expect_identical(
    call_exponential_balance(problem, max_categories = -1L)$state,
    "not_evaluated_control"
  )
})

test_that("current-fit wrapper reconstructs the canonical zero-utility escape", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_exponential_balance_fixture()
  scope <- fit$config$boundary_audit$gpcm_exponential_balance
  zero_center <- rep(0, scope$free_parameter_dimension -
    fit$config$gpcm_spec$n_params)
  rate <- c(1, -1, 0)
  wrapped <- mfrmr:::mfrmr_jml_gpcm_fit_exponential_balance_limit(
    fit,
    additive_coordinate_directions = list(),
    additive_exponential_rates = numeric(0),
    slope_log_rates = rate,
    additive_center = zero_center
  )

  expect_identical(
    wrapped$contract_version,
    "mfrmr-jml-gpcm-exponential-balance-0.2.3-v1"
  )
  expect_identical(
    wrapped$state,
    "current_fit_bounded_response_image_parameter_escape_certified"
  )
  expect_true(wrapped$evaluated)
  expect_true(wrapped$current_fit_reconstructed)
  expect_true(wrapped$exponential_balance_limit_classified)
  expect_true(wrapped$bounded_response_image_escape_certified)
  expect_true(wrapped$parameter_to_response_contrast_map_properness_refuted)
  expect_true(
    wrapped$result$base_log_slopes_sum_zero_structurally_certified
  )
  expect_equal(
    wrapped$log_likelihood_limit,
    -nrow(fit$prep$data) * log(fit$config$n_cat), tolerance = 1e-12
  )
  expect_true(is.finite(wrapped$retained_log_likelihood))
  expect_true(is.finite(wrapped$boundary_improvement_over_retained))
  expect_false(wrapped$path_inferred_from_optimizer_trace)
  expect_false(wrapped$path_search_performed)
  expect_false(wrapped$competitive_boundary_certified)
  expect_false(wrapped$finite_jmle_existence_certified)
  expect_false(wrapped$finite_jmle_nonexistence_certified)
  expect_false(wrapped$global_boundary_absence_certified)
  expect_false(wrapped$external_comparison_eligible)
  expect_identical(wrapped$readiness_effect, "none_diagnostic_only")
})

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

test_that("a nonzero anchor closes the canonical production certificate", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_exponential_balance_fixture(data.frame(
    Facet = "Criterion",
    Level = "C1",
    Anchor = 0.25,
    stringsAsFactors = FALSE
  ))
  audit <- fit$config$boundary_audit$gpcm_exponential_balance

  expect_identical(
    audit$state, "declared_exponential_balance_oracle_available_no_path_declared"
  )
  expect_true(audit$evaluated)
  expect_true(audit$declared_exponential_balance_oracle_available)
  expect_false(audit$canonical_zero_utility_center_certified)
  expect_false(audit$canonical_affine_additive_offset_exactly_zero)
  expect_false(audit$canonical_zero_utility_divergent_slope_path_certified)
  expect_false(audit$production_bounded_image_escape_certified)
  expect_false(audit$parameter_to_response_contrast_map_properness_classified)
  expect_false(audit$parameter_to_response_contrast_map_properness_refuted)

  zero_center <- rep(
    0,
    audit$free_parameter_dimension - fit$config$gpcm_spec$n_params
  )
  wrapped <- mfrmr:::mfrmr_jml_gpcm_fit_exponential_balance_limit(
    fit,
    slope_log_rates = c(1, -1, 0),
    additive_center = zero_center
  )
  expect_identical(
    wrapped$state,
    "current_fit_divergent_parameter_and_response_flag_limit_classified"
  )
  expect_false(wrapped$bounded_response_image_escape_certified)
  expect_false(wrapped$parameter_to_response_contrast_map_properness_refuted)
})

test_that("model estimator unit-slope source and workload boundaries are typed", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_exponential_balance_fixture()
  audit <- fit$config$boundary_audit$gpcm_exponential_balance
  sizes <- mfrmr:::build_param_sizes(fit$config)

  pcm_config <- fit$config
  pcm_config$model <- "PCM"
  pcm <- mfrmr:::audit_mfrm_jml_gpcm_exponential_balance_scope(
    pcm_config, sizes,
    fit$config$boundary_audit$gpcm_parameter_sequence_flag,
    fit$config$boundary_audit$gpcm_global_existence
  )
  expect_identical(pcm$state, "not_applicable_model")
  expect_false(pcm$evaluated)
  expect_identical(pcm$estimator_identity, "not_applicable")

  mml_config <- fit$config
  mml_config$method <- "MML"
  mml <- mfrmr:::audit_mfrm_jml_gpcm_exponential_balance_scope(
    mml_config, sizes,
    fit$config$boundary_audit$gpcm_parameter_sequence_flag,
    fit$config$boundary_audit$gpcm_global_existence
  )
  expect_identical(mml$state, "not_applicable_estimator")
  expect_false(mml$evaluated)
  expect_identical(mml$objective_identity, "not_applicable_marginal_estimator")

  unit_config <- fit$config
  unit_config$gpcm_spec$n_params <- 0L
  unit <- mfrmr:::audit_mfrm_jml_gpcm_exponential_balance_scope(
    unit_config, sizes,
    fit$config$boundary_audit$gpcm_parameter_sequence_flag,
    fit$config$boundary_audit$gpcm_global_existence
  )
  expect_identical(unit$state, "not_required_unit_slope")
  expect_false(unit$evaluated)
  expect_false(unit$parameter_to_response_contrast_map_properness_refuted)

  legacy <- fit$config$boundary_audit$gpcm_parameter_sequence_flag
  legacy$contract_version <- "legacy"
  source <- mfrmr:::audit_mfrm_jml_gpcm_exponential_balance_scope(
    fit$config, sizes, legacy,
    fit$config$boundary_audit$gpcm_global_existence
  )
  expect_identical(source$state, "not_evaluated_source_contract")
  expect_false(source$evaluated)
  expect_false(source$source_parameter_sequence_flag_contract_matches)

  workload <- mfrmr:::audit_mfrm_jml_gpcm_exponential_balance_scope(
    fit$config, sizes,
    fit$config$boundary_audit$gpcm_parameter_sequence_flag,
    fit$config$boundary_audit$gpcm_global_existence,
    max_oracle_free_coordinates = 1L
  )
  expect_identical(
    workload$state, "exponential_balance_oracle_workload_not_admitted"
  )
  expect_true(workload$evaluated)
  expect_false(workload$declared_exponential_balance_oracle_available)
  expect_false(workload$canonical_zero_utility_center_certified)
  expect_false(workload$canonical_zero_utility_divergent_slope_path_certified)
  expect_false(workload$parameter_to_response_contrast_map_properness_refuted)

  canonical_without_oracle <-
    mfrmr:::audit_mfrm_jml_gpcm_exponential_balance_scope(
      fit$config, sizes,
      fit$config$boundary_audit$gpcm_parameter_sequence_flag,
      fit$config$boundary_audit$gpcm_global_existence,
      canonical_zero_utility_center_certified = TRUE,
      max_oracle_free_coordinates = 1L
    )
  expect_identical(
    canonical_without_oracle$state,
    "canonical_zero_utility_bounded_image_escape_certified"
  )
  expect_false(
    canonical_without_oracle$declared_exponential_balance_oracle_available
  )
  expect_true(
    canonical_without_oracle$canonical_zero_utility_divergent_slope_path_certified
  )
  expect_true(
    canonical_without_oracle$parameter_to_response_contrast_map_properness_refuted
  )

  invalid <- mfrmr:::audit_mfrm_jml_gpcm_exponential_balance_scope(
    fit$config, sizes,
    fit$config$boundary_audit$gpcm_parameter_sequence_flag,
    fit$config$boundary_audit$gpcm_global_existence,
    max_oracle_free_coordinates = -1L
  )
  expect_identical(invalid$state, "not_evaluated_control")
  expect_false(invalid$evaluated)

  invalid_fit <- fit
  invalid_fit$config$boundary_audit$gpcm_exponential_balance$contract_version <-
    "legacy"
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_fit_exponential_balance_limit(
      invalid_fit, slope_log_rates = c(1, -1, 0)
    )$state,
    "not_evaluated_source_contract"
  )
  expect_true(audit$parameter_to_response_contrast_map_properness_refuted)
})
