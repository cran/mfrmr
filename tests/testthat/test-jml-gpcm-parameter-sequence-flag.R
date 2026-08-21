fit_jml_gpcm_parameter_sequence_flag_fixture <- function() {
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

parameter_sequence_flag_problem <- function() {
  list(
    base_contrasts = rbind(
      c(0.2, -0.1),
      c(-0.5, 0.3),
      c(0.2, 0.8),
      c(0.1, -0.2)
    ),
    score_k = c(2L, 1L, 2L, 0L),
    stage_directions = list(
      rbind(
        c(1, 1),
        c(2, -1),
        c(0, 0),
        c(0, 0)
      ),
      rbind(
        c(0, 1),
        c(0, 1),
        c(1, 1),
        c(0, 0)
      )
    ),
    weight = c(1, 1, 2, 0)
  )
}

call_parameter_sequence_flag <- function(problem, ...) {
  do.call(
    mfrmr:::mfrmr_jml_gpcm_declared_contrast_flag_limit,
    c(problem, list(...))
  )
}

test_that("reference contrasts remove exactly the rowwise logit gauge", {
  logits <- rbind(
    c(2, 4, -1),
    c(-3, 0.5, 7)
  )
  shifts <- c(100, -25)
  shifted <- logits + matrix(shifts, nrow = 2L, ncol = 3L)

  original <- mfrmr:::mfrmr_jml_gpcm_reference_contrasts(logits)
  transported <- mfrmr:::mfrmr_jml_gpcm_reference_contrasts(shifted)

  expect_equal(original, rbind(c(2, -3), c(3.5, 10)))
  expect_identical(transported, original)
  expect_null(
    mfrmr:::mfrmr_jml_gpcm_reference_contrasts(matrix(1, 2L, 1L))
  )
  expect_null(
    mfrmr:::mfrmr_jml_gpcm_reference_contrasts(
      matrix(c(0, Inf), 1L, 2L)
    )
  )
})

test_that("a two-stage contrast flag has an exact finite likelihood limit", {
  problem <- parameter_sequence_flag_problem()
  result <- call_parameter_sequence_flag(problem)
  expected_row_three <- 0.8 - mfrmr:::logsumexp(c(0.2, 0.8))

  expect_identical(
    result$contract_version,
    "mfrmr-jml-gpcm-parameter-sequence-flag-0.2.3-v1"
  )
  expect_identical(
    result$state,
    "finite_declared_contrast_flag_log_likelihood_limit"
  )
  expect_true(result$evaluated)
  expect_true(result$declared_contrast_flag_limit_classified)
  expect_identical(
    result$response_coordinate_system,
    "category_logit_contrasts_relative_to_category_zero"
  )
  expect_identical(result$observations, 4L)
  expect_identical(result$categories, 3L)
  expect_identical(result$contrast_dimension, 8)
  expect_identical(result$divergent_flag_stages, 2L)
  expect_identical(
    result$scale_contract,
    "s_g_to_infinity_and_s_(g+1)_over_s_g_to_zero"
  )
  expect_identical(
    result$contrast_remainder_contract,
    "uniform_reference_contrast_o(1)"
  )
  expect_true(result$exact_coefficient_comparisons)
  expect_true(result$exact_lexicographic_ties)
  expect_equal(
    result$row_limits$RowLogProbabilityLimit,
    c(0, 0, expected_row_three, -mfrmr:::logsumexp(c(0, 0.1, -0.2))),
    tolerance = 1e-14
  )
  expect_equal(
    result$log_likelihood_limit, 2 * expected_row_three,
    tolerance = 1e-14
  )
  expect_true(result$finite_log_likelihood_limit)
  expect_false(result$negative_infinite_log_likelihood_limit)
  expect_identical(result$stage_table$RowsResolved, c(2L, 2L))
  expect_identical(result$stage_table$CategoriesRemoved, c(3L, 2L))
  expect_true(result$direct_finite_scale_evaluation_available)
  expect_false(result$parameter_sequence_reachability_checked)
  expect_false(result$parameter_flag_extracted)
  expect_false(result$optimizer_sequence_used)
  expect_false(result$p2b_handoff_eligible)
  expect_false(result$original_sequence_limit_classified)
  expect_false(result$common_subsequence_limit_certified)
  expect_false(result$arbitrary_parameter_sequence_global_limit_classified)
  expect_false(result$competitive_boundary_certified)
  expect_false(result$global_boundary_classified)
  expect_false(result$global_finite_maximum_certified)
  expect_false(result$global_boundary_absence_certified)
  expect_false(result$standard_error_eligible)
  expect_false(result$confidence_interval_eligible)
  expect_false(result$external_comparison_eligible)
  expect_identical(result$readiness_effect, "none_diagnostic_only")
})

test_that("bounded contrast images reduce to their convergent base softmax", {
  base <- rbind(c(0.5, -0.2), c(-0.4, 0.7))
  result <- mfrmr:::mfrmr_jml_gpcm_declared_contrast_flag_limit(
    base_contrasts = base,
    score_k = c(1L, 0L),
    stage_directions = list(),
    weight = c(2, 1)
  )
  expected <- c(
    0.5 - mfrmr:::logsumexp(c(0, 0.5, -0.2)),
    -mfrmr:::logsumexp(c(0, -0.4, 0.7))
  )

  expect_true(result$declared_contrast_flag_limit_classified)
  expect_identical(result$divergent_flag_stages, 0L)
  expect_identical(nrow(result$stage_table), 0L)
  expect_equal(
    result$row_limits$RowLogProbabilityLimit, expected,
    tolerance = 1e-14
  )
  expect_equal(result$log_likelihood_limit, 2 * expected[1] + expected[2])

  direct <- mfrmr:::mfrmr_jml_gpcm_declared_contrast_flag_loglik_at(
    result, numeric(0)
  )
  expect_true(direct$valid)
  expect_equal(direct$row_log_probability, expected, tolerance = 1e-14)
  expect_equal(direct$log_likelihood, result$log_likelihood_limit)
})

test_that("observed-category exclusion yields a negative-infinite limit", {
  problem <- parameter_sequence_flag_problem()
  problem$score_k[1L] <- 0L
  result <- call_parameter_sequence_flag(problem)

  expect_identical(
    result$state,
    "negative_infinite_declared_contrast_flag_log_likelihood_limit"
  )
  expect_true(result$declared_contrast_flag_limit_classified)
  expect_true(result$negative_infinite_log_likelihood_limit)
  expect_false(result$finite_log_likelihood_limit)
  expect_identical(result$row_limits$ObservedExclusionStage[1L], 1L)
  expect_false(result$row_limits$ObservedInLimitSupport[1L])
  expect_identical(result$log_likelihood_limit, -Inf)

  problem$weight[1L] <- 0
  zero_weight <- call_parameter_sequence_flag(problem)
  expect_true(zero_weight$finite_log_likelihood_limit)
  expect_true(is.finite(zero_weight$log_likelihood_limit))
})

test_that("non-power scale flags converge with a vanishing contrast remainder", {
  result <- call_parameter_sequence_flag(parameter_sequence_flag_problem())
  indices <- c(9, 16, 25, 36, 49, 64)
  remainder_direction <- rbind(
    c(1, -1), c(-0.5, 0.5), c(2, -1), c(0, 0)
  )
  deviations <- vapply(indices, function(index) {
    scales <- c(exp(sqrt(index)), log(index))
    direct <- mfrmr:::mfrmr_jml_gpcm_declared_contrast_flag_loglik_at(
      result,
      scales = scales,
      contrast_remainder = remainder_direction / index
    )
    expect_true(direct$valid)
    max(abs(
      direct$row_log_probability[result$row_limits$Weight > 0] -
        result$row_limits$RowLogProbabilityLimit[
          result$row_limits$Weight > 0
        ]
    ))
  }, numeric(1))

  expect_true(all(diff(deviations) < 0))
  expect_lt(tail(deviations, 1L), head(deviations, 1L) / 3)

  invalid_scale <-
    mfrmr:::mfrmr_jml_gpcm_declared_contrast_flag_loglik_at(
      result, c(1, 2)
    )
  expect_false(invalid_scale$valid)
  invalid_remainder <-
    mfrmr:::mfrmr_jml_gpcm_declared_contrast_flag_loglik_at(
      result, c(4, 2), matrix(0, 1L, 1L)
    )
  expect_false(invalid_remainder$valid)
})

test_that("the contrast oracle agrees with the fixed-slope P2b special case", {
  base_utilities <- rbind(
    c(0, 0.2, -0.1),
    c(0, -0.5, 0.3),
    c(0, 0.2, 0.8)
  )
  directions <- list(
    rbind(c(0, 1, 1), c(0, 2, -1), c(0, 0, 0)),
    rbind(c(0, 0, 1), c(0, 0, 1), c(0, 1, 1))
  )
  scores <- c(2L, 1L, 2L)
  weights <- c(1, 1, 2)
  p2b <- mfrmr:::mfrmr_jml_gpcm_declared_hierarchy_limit(
    base_utilities = base_utilities,
    score_k = scores,
    slope_index = c(1L, 2L, 1L),
    base_log_slopes = c(0, 0),
    additive_stage_directions = directions,
    additive_scale_exponents = c(1, 0.5),
    weight = weights,
    slope_levels = c("A", "B")
  )
  p2e <- mfrmr:::mfrmr_jml_gpcm_declared_contrast_flag_limit(
    base_contrasts =
      mfrmr:::mfrmr_jml_gpcm_reference_contrasts(base_utilities),
    score_k = scores,
    stage_directions = lapply(
      directions, mfrmr:::mfrmr_jml_gpcm_reference_contrasts
    ),
    weight = weights
  )

  expect_true(p2b$declared_path_limit_classified)
  expect_true(p2e$declared_contrast_flag_limit_classified)
  expect_equal(
    p2e$row_limits$RowLogProbabilityLimit,
    p2b$row_limits$RowLogProbabilityLimit,
    tolerance = 1e-14
  )
  expect_equal(
    p2e$log_likelihood_limit, p2b$log_likelihood_limit,
    tolerance = 1e-14
  )
  expect_false(p2e$p2b_handoff_eligible)
})

test_that("declared contrast flags reject malformed contracts", {
  problem <- parameter_sequence_flag_problem()

  bad_base <- problem
  bad_base$base_contrasts <- matrix(NA_real_, 4L, 2L)
  expect_identical(
    call_parameter_sequence_flag(bad_base)$state,
    "not_evaluated_base"
  )

  bad_stage <- problem
  bad_stage$stage_directions <- list(matrix(0, 1L, 1L))
  expect_identical(
    call_parameter_sequence_flag(bad_stage)$state,
    "not_evaluated_stage_mapping"
  )

  zero_stage <- problem
  zero_stage$stage_directions <- list(matrix(0, 4L, 2L))
  expect_identical(
    call_parameter_sequence_flag(zero_stage)$state,
    "not_evaluated_zero_stage"
  )

  one_dimension <- list(
    base_contrasts = matrix(0, 1L, 1L),
    score_k = 0L,
    stage_directions = list(matrix(1, 1L, 1L), matrix(2, 1L, 1L)),
    weight = 1
  )
  expect_identical(
    call_parameter_sequence_flag(one_dimension)$state,
    "not_evaluated_stage_limit"
  )

  bad_score <- problem
  bad_score$score_k[1L] <- 3L
  expect_identical(
    call_parameter_sequence_flag(bad_score)$state,
    "not_evaluated_observation_map"
  )
  expect_identical(
    call_parameter_sequence_flag(problem, max_observations = 3L)$state,
    "not_evaluated_size_limit"
  )
  expect_identical(
    call_parameter_sequence_flag(
      problem, max_contrast_elements = 20L
    )$state,
    "not_evaluated_size_limit"
  )
  expect_identical(
    call_parameter_sequence_flag(problem, max_stages = -1L)$state,
    "not_evaluated_control"
  )
})

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

test_that("explicit finite parameter points use the exact nonlinear map", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_parameter_sequence_flag_fixture()
  sequence <- rbind(fit$opt$par, fit$opt$par)
  sequence[2L, 1L] <- sequence[2L, 1L] + 0.25
  result <- mfrmr:::mfrmr_jml_gpcm_parameter_sequence_contrast_image(
    fit, sequence
  )

  expect_identical(
    result$state, "finite_parameter_sequence_exact_contrast_image"
  )
  expect_true(result$evaluated)
  expect_true(result$source_sequence_declared_by_caller)
  expect_true(result$exact_nonlinear_parameter_to_contrast_map_evaluated)
  expect_true(result$slope_utility_product_evaluated_before_compactification)
  expect_identical(result$parameter_sequence_points, 2L)
  expect_identical(
    result$contrast_dimension,
    as.integer(nrow(fit$prep$data) * (fit$config$n_cat - 1L))
  )
  expect_identical(
    nrow(result$coordinate_map), result$contrast_dimension
  )

  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = fit$config$step_facet,
    slope_facet = fit$config$slope_facet,
    interaction_specs = fit$config$interaction_specs
  )
  sizes <- mfrmr:::build_param_sizes(fit$config)
  adjacent <- mfrmr:::mfrmr_gpcm_adjacent_logits(
    sequence[1L, ], idx, fit$config, sizes
  )
  adjacent <- matrix(
    adjacent,
    nrow = nrow(fit$prep$data),
    ncol = fit$config$n_cat - 1L
  )
  expected <- as.numeric(t(apply(adjacent, 1L, cumsum)))
  expect_equal(result$contrast_sequence[1L, ], expected, tolerance = 0)
  expect_true(result$finite_sequence_only)
  expect_false(result$asymptotic_flag_extracted)
  expect_false(result$asymptotic_scales_certified)
  expect_false(result$subsequence_selected)
  expect_false(result$p2b_handoff_eligible)
  expect_false(result$global_boundary_classified)
  expect_false(result$external_comparison_eligible)
  expect_identical(result$readiness_effect, "none_diagnostic_only")
})

test_that("finite parameter mapping and theorem scopes fail closed", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_parameter_sequence_flag_fixture()
  sizes <- mfrmr:::build_param_sizes(fit$config)
  compactification <-
    fit$config$boundary_audit$gpcm_boundary_compactification
  sequence_diagnostic <-
    fit$config$boundary_audit$gpcm_optimizer_sequence_diagnostic

  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_parameter_sequence_contrast_image(
      unclass(fit), matrix(fit$opt$par, nrow = 1L)
    )$state,
    "not_evaluated_fit"
  )
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_parameter_sequence_contrast_image(
      fit, matrix(0, 1L, 1L)
    )$state,
    "not_evaluated_sequence"
  )
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_parameter_sequence_contrast_image(
      fit, matrix(fit$opt$par, nrow = 1L), max_points = 0L
    )$state,
    "not_evaluated_control"
  )
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_parameter_sequence_contrast_image(
      fit, rbind(fit$opt$par, fit$opt$par), max_points = 1L
    )$state,
    "not_evaluated_size_limit"
  )

  pcm <- fit$config
  pcm$model <- "PCM"
  expect_identical(
    mfrmr:::audit_mfrm_jml_gpcm_parameter_sequence_flag_scope(
      pcm, sizes, nrow(fit$prep$data), fit$config$n_cat,
      compactification, sequence_diagnostic
    )$state,
    "not_applicable_model"
  )
  mml <- fit$config
  mml$method <- "MML"
  expect_identical(
    mfrmr:::audit_mfrm_jml_gpcm_parameter_sequence_flag_scope(
      mml, sizes, nrow(fit$prep$data), fit$config$n_cat,
      compactification, sequence_diagnostic
    )$state,
    "not_applicable_estimator"
  )
  unit <- fit$config
  unit$gpcm_spec$n_params <- 0L
  expect_identical(
    mfrmr:::audit_mfrm_jml_gpcm_parameter_sequence_flag_scope(
      unit, sizes, nrow(fit$prep$data), fit$config$n_cat,
      compactification, sequence_diagnostic
    )$state,
    "not_required_unit_slope"
  )
  expect_identical(
    mfrmr:::audit_mfrm_jml_gpcm_parameter_sequence_flag_scope(
      fit$config, sizes, 0L, fit$config$n_cat,
      compactification, sequence_diagnostic
    )$state,
    "not_evaluated_dimension"
  )
  legacy <- compactification
  legacy$contract_version <- "legacy"
  expect_identical(
    mfrmr:::audit_mfrm_jml_gpcm_parameter_sequence_flag_scope(
      fit$config, sizes, nrow(fit$prep$data), fit$config$n_cat,
      legacy, sequence_diagnostic
    )$state,
    "not_evaluated_source_contract"
  )
  workload <- mfrmr:::audit_mfrm_jml_gpcm_parameter_sequence_flag_scope(
    fit$config, sizes, nrow(fit$prep$data), fit$config$n_cat,
    compactification, sequence_diagnostic,
    max_oracle_contrast_dimension = 1L
  )
  expect_identical(
    workload$state,
    "parameter_sequence_flag_theorem_available_oracle_workload_not_admitted"
  )
  expect_true(workload$evaluated)
  expect_true(
    workload$every_parameter_sequence_has_classifiable_further_subsequence_likelihood_limit
  )
  expect_false(workload$declared_contrast_flag_oracle_available)
  expect_identical(
    mfrmr:::audit_mfrm_jml_gpcm_parameter_sequence_flag_scope(
      fit$config, sizes, nrow(fit$prep$data), fit$config$n_cat,
      compactification, sequence_diagnostic,
      max_oracle_contrast_dimension = 0L
    )$state,
    "not_evaluated_control"
  )
})
