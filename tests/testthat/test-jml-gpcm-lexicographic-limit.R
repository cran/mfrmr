fit_jml_gpcm_lexicographic_limit_fixture <- function() {
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

two_stage_lexicographic_problem <- function() {
  list(
    base_utilities = rbind(
      c(0, 0, 0),
      c(0, 3, -2),
      c(0, 1, 1)
    ),
    score_k = c(2L, 1L, 2L),
    slope_index = 1:3,
    base_log_slopes = c(0, 0, 0),
    slope_stage_rates = rbind(
      c(1, -1, 0),
      c(-0.5, -0.5, 1)
    ),
    slope_scale_exponents = c(1, 0.5),
    additive_stage_directions = list(
      rbind(c(0, 1, 1), c(0, -1, 2), c(0, 0, 0)),
      rbind(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0))
    ),
    additive_scale_exponents = c(1, 0.5),
    weight = c(1, 1, 1),
    slope_levels = c("C1", "C2", "C3")
  )
}

call_lexicographic_limit <- function(problem, ...) {
  do.call(
    mfrmr:::mfrmr_jml_gpcm_declared_hierarchy_limit,
    c(problem, list(...))
  )
}

call_lexicographic_at <- function(problem, distance, ...) {
  do.call(
    mfrmr:::mfrmr_jml_gpcm_declared_hierarchy_loglik_at,
    c(list(distance = distance), problem, list(...))
  )
}

test_that("two-stage joint hierarchy has an exact finite likelihood limit", {
  problem <- two_stage_lexicographic_problem()
  oracle <- call_lexicographic_limit(problem)

  expect_identical(
    oracle$contract_version,
    "mfrmr-jml-gpcm-lexicographic-limit-0.2.3-v1"
  )
  expect_identical(
    oracle$state, "finite_declared_hierarchy_log_likelihood_limit"
  )
  expect_true(oracle$evaluated)
  expect_true(oracle$declared_path_limit_classified)
  expect_true(oracle$exact_coefficient_comparisons)
  expect_identical(oracle$slope_stages, 2L)
  expect_identical(oracle$additive_stages, 2L)
  expect_identical(oracle$maximum_slope_stages, 2L)
  expect_identical(oracle$maximum_additive_stages, 2L)
  expect_identical(oracle$observations, 3L)
  expect_identical(oracle$categories, 3L)
  expect_identical(oracle$slope_levels, 3L)
  expect_identical(
    oracle$slope_hierarchy_state, "complete_finite_rate_hierarchy"
  )
  expect_equal(oracle$log_likelihood_limit, -log(6), tolerance = 1e-14)
  expect_true(oracle$finite_log_likelihood_limit)
  expect_false(oracle$negative_infinite_log_likelihood_limit)
  expect_false(oracle$positive_infinite_log_likelihood_limit)
  expect_true(oracle$direct_finite_distance_evaluation_available)
  expect_false(
    oracle$declared_utility_path_parameter_reachability_checked
  )
  expect_false(oracle$path_inferred_from_optimizer_trace)
  expect_false(oracle$path_search_performed)
  expect_false(oracle$monotone_tail_certified)
  expect_false(oracle$competitive_boundary_certified)
  expect_false(oracle$common_subsequence_limit_certified)
  expect_false(oracle$arbitrary_path_limit_classified)
  expect_false(oracle$global_boundary_classified)
  expect_false(oracle$global_finite_maximum_certified)
  expect_false(oracle$global_boundary_absence_certified)
  expect_false(oracle$standard_error_eligible)
  expect_false(oracle$confidence_interval_eligible)
  expect_false(oracle$external_comparison_eligible)
  expect_identical(oracle$readiness_effect, "none_diagnostic_only")

  slopes <- oracle$slope_table
  expect_identical(slopes$LimitState, c("infinite", "zero", "infinite"))
  expect_identical(slopes$ResolutionStage, c(1L, 1L, 2L))
  expect_equal(slopes$ResolutionCoefficient, c(1, -1, 1))

  rows <- oracle$row_limits
  expect_identical(
    rows$EvaluationState,
    c(
      "infinite_slope_lexicographic_concentration",
      "uniform_vanishing_slope",
      "infinite_slope_lexicographic_concentration"
    )
  )
  expect_identical(rows$AdditiveDivergentTieCategories, c("2", "2", "0;1;2"))
  expect_identical(rows$LimitSupportCategories, c("2", "0;1;2", "1;2"))
  expect_identical(rows$LimitSupportSize, c(1L, 3L, 2L))
  expect_true(all(rows$ObservedInLimitSupport))
  expect_equal(
    rows$RowLogProbabilityLimit,
    c(0, -log(3), -log(2)), tolerance = 1e-14
  )
  expect_true(all(rows$FiniteRowLimit))
})

test_that("direct finite-distance evaluation converges to the analytic oracle", {
  problem <- two_stage_lexicographic_problem()
  oracle <- call_lexicographic_limit(problem)
  distances <- c(4, 16, 64, 256)
  evaluated <- lapply(distances, function(distance) {
    call_lexicographic_at(problem, distance)
  })
  valid <- vapply(evaluated, `[[`, logical(1), "valid")
  values <- vapply(evaluated, `[[`, numeric(1), "log_likelihood")

  expect_true(all(valid))
  expect_true(all(is.finite(values)))
  expect_lt(abs(values[4] - oracle$log_likelihood_limit), 1e-8)
  expect_lt(
    abs(values[4] - oracle$log_likelihood_limit),
    abs(values[1] - oracle$log_likelihood_limit)
  )
  expect_equal(
    evaluated[[4]]$row_log_probability,
    oracle$row_limits$RowLogProbabilityLimit,
    tolerance = 1e-8
  )
  expect_identical(
    evaluated[[1]]$reason_codes,
    "declared_hierarchy_finite_distance_evaluated"
  )

  invalid <- call_lexicographic_at(problem, 0)
  expect_false(invalid$valid)
  expect_true(is.na(invalid$log_likelihood))
  expect_identical(
    invalid$reason_codes, "lexicographic_limit_distance_invalid"
  )
})

test_that("finite slopes use a softmax only on the additive tie set", {
  base <- matrix(c(0, 0, 1), nrow = 1L)
  direction <- matrix(c(0, 1, 1), nrow = 1L)
  oracle <- mfrmr:::mfrmr_jml_gpcm_declared_hierarchy_limit(
    base_utilities = base,
    score_k = 2L,
    slope_index = 3L,
    base_log_slopes = c(log(2), -log(2), 0),
    slope_stage_rates = c(1, -1, 0),
    slope_scale_exponents = 1,
    additive_stage_directions = list(direction),
    additive_scale_exponents = 1,
    slope_levels = c("C1", "C2", "C3")
  )
  expected <- -log1p(exp(-1))

  expect_identical(
    oracle$state, "finite_declared_hierarchy_log_likelihood_limit"
  )
  expect_identical(oracle$slope_table$LimitState, c("infinite", "zero", "finite"))
  expect_identical(
    oracle$row_limits$EvaluationState,
    "finite_slope_softmax_on_divergent_tie_set"
  )
  expect_identical(oracle$row_limits$AdditiveDivergentTieCategories, "1;2")
  expect_identical(oracle$row_limits$LimitSupportCategories, "1;2")
  expect_identical(oracle$row_limits$LimitSupportSize, 2L)
  expect_true(oracle$row_limits$ObservedInLimitSupport)
  expect_equal(oracle$row_limits$RowLogProbabilityLimit, expected)
  expect_equal(oracle$log_likelihood_limit, expected)

  direct <- vapply(c(4, 16, 64), function(distance) {
    mfrmr:::mfrmr_jml_gpcm_declared_hierarchy_loglik_at(
      distance = distance,
      base_utilities = base,
      score_k = 2L,
      slope_index = 3L,
      base_log_slopes = c(log(2), -log(2), 0),
      slope_stage_rates = c(1, -1, 0),
      slope_scale_exponents = 1,
      additive_stage_directions = list(direction),
      additive_scale_exponents = 1
    )$log_likelihood
  }, numeric(1))
  expect_gt(abs(direct[1] - expected), abs(direct[2] - expected))
  expect_gt(abs(direct[2] - expected), abs(direct[3] - expected))
  expect_equal(direct[3], expected, tolerance = 1e-12)
})

test_that("excluded observed categories yield a negative infinite limit", {
  oracle <- mfrmr:::mfrmr_jml_gpcm_declared_hierarchy_limit(
    base_utilities = matrix(c(0, 0, 0), nrow = 1L),
    score_k = 0L,
    slope_index = 1L,
    base_log_slopes = c(0, 0),
    slope_stage_rates = c(1, -1),
    slope_scale_exponents = 1,
    additive_stage_directions = list(matrix(c(0, 1, 2), nrow = 1L)),
    additive_scale_exponents = 1,
    slope_levels = c("C1", "C2")
  )

  expect_identical(
    oracle$state,
    "negative_infinite_declared_hierarchy_log_likelihood_limit"
  )
  expect_true(oracle$declared_path_limit_classified)
  expect_true(is.infinite(oracle$log_likelihood_limit))
  expect_lt(oracle$log_likelihood_limit, 0)
  expect_false(oracle$finite_log_likelihood_limit)
  expect_true(oracle$negative_infinite_log_likelihood_limit)
  expect_false(oracle$row_limits$ObservedInLimitSupport)
  expect_identical(oracle$row_limits$LimitSupportCategories, "2")
  expect_identical(
    oracle$row_limits$ReasonCodes,
    "observed_category_excluded_from_infinite_slope_limit_support"
  )

  zero_weight <- mfrmr:::mfrmr_jml_gpcm_declared_hierarchy_limit(
    base_utilities = rbind(c(0, 0, 0), c(0, 0, 0)),
    score_k = c(2L, 0L),
    slope_index = c(1L, 1L),
    base_log_slopes = c(0, 0),
    slope_stage_rates = c(1, -1),
    slope_scale_exponents = 1,
    additive_stage_directions = list(
      rbind(c(0, 1, 2), c(0, 1, 2))
    ),
    additive_scale_exponents = 1,
    weight = c(1, 0)
  )
  expect_true(zero_weight$finite_log_likelihood_limit)
  expect_equal(zero_weight$log_likelihood_limit, 0)
  expect_false(zero_weight$row_limits$FiniteRowLimit[2])
})

test_that("vanishing slopes dominate every admitted polynomial additive scale", {
  oracle <- mfrmr:::mfrmr_jml_gpcm_declared_hierarchy_limit(
    base_utilities = matrix(c(0, 5, -4), nrow = 1L),
    score_k = 1L,
    slope_index = 2L,
    base_log_slopes = c(0, 0, 0),
    slope_stage_rates = rbind(c(1, -1, 0), c(-0.5, -0.5, 1)),
    slope_scale_exponents = c(1, 0.5),
    additive_stage_directions = list(
      matrix(c(0, -100, 200), nrow = 1L),
      matrix(c(0, 500, -300), nrow = 1L)
    ),
    additive_scale_exponents = c(2, 0.25),
    slope_levels = c("C1", "C2", "C3")
  )

  expect_identical(oracle$row_limits$SlopeLimitState, "zero")
  expect_identical(
    oracle$row_limits$EvaluationState, "uniform_vanishing_slope"
  )
  expect_identical(oracle$row_limits$LimitSupportCategories, "0;1;2")
  expect_equal(oracle$row_limits$RowLogProbabilityLimit, -log(3))
  expect_equal(oracle$log_likelihood_limit, -log(3))

  direct <- vapply(c(16, 64, 256), function(distance) {
    mfrmr:::mfrmr_jml_gpcm_declared_hierarchy_loglik_at(
      distance = distance,
      base_utilities = matrix(c(0, 5, -4), nrow = 1L),
      score_k = 1L,
      slope_index = 2L,
      base_log_slopes = c(0, 0, 0),
      slope_stage_rates = rbind(c(1, -1, 0), c(-0.5, -0.5, 1)),
      slope_scale_exponents = c(1, 0.5),
      additive_stage_directions = list(
        matrix(c(0, -100, 200), nrow = 1L),
        matrix(c(0, 500, -300), nrow = 1L)
      ),
      additive_scale_exponents = c(2, 0.25)
    )$log_likelihood
  }, numeric(1))
  expect_lt(abs(direct[3] + log(3)), abs(direct[1] + log(3)))
  expect_lt(abs(direct[3] + log(3)), 1e-10)
})

test_that("declared-path contract failures remain typed and fail closed", {
  problem <- two_stage_lexicographic_problem()

  invalid_base <- problem
  invalid_base$base_utilities[1, 1] <- NA_real_
  expect_identical(
    call_lexicographic_limit(invalid_base)$state, "not_evaluated_base"
  )

  invalid_base_log <- problem
  invalid_base_log$base_log_slopes <- c(1, 0, 0)
  expect_identical(
    call_lexicographic_limit(invalid_base_log)$state,
    "not_evaluated_identification"
  )

  malformed_base_log <- problem
  malformed_base_log$base_log_slopes <- list("bad", 0, 0)
  expect_identical(
    call_lexicographic_limit(malformed_base_log)$state,
    "not_evaluated_base"
  )

  invalid_levels <- problem
  invalid_levels$slope_levels <- c("C1", "C1", "C3")
  expect_identical(
    call_lexicographic_limit(invalid_levels)$state, "not_evaluated_levels"
  )

  invalid_mapping <- problem
  invalid_mapping$additive_stage_directions[[1]] <- matrix(0, 1L, 3L)
  expect_identical(
    call_lexicographic_limit(invalid_mapping)$state,
    "not_evaluated_stage_mapping"
  )

  too_many <- problem
  too_many$slope_stage_rates <- rbind(
    problem$slope_stage_rates,
    c(1, 0, -1)
  )
  too_many$slope_scale_exponents <- c(1, 0.5, 0.25)
  expect_identical(
    call_lexicographic_limit(too_many)$state,
    "not_evaluated_stage_limit"
  )

  invalid_scale <- problem
  invalid_scale$additive_scale_exponents <- c(1, 1)
  expect_identical(
    call_lexicographic_limit(invalid_scale)$state,
    "not_evaluated_scale_order"
  )

  invalid_rate <- problem
  invalid_rate$slope_stage_rates[2, ] <- c(1, 0, 0)
  expect_identical(
    call_lexicographic_limit(invalid_rate)$state,
    "not_evaluated_identification"
  )

  redundant <- problem
  redundant$additive_stage_directions[[2]] <- matrix(
    rep(c(1, 1, 1), 3L), nrow = 3L, byrow = TRUE
  )
  expect_identical(
    call_lexicographic_limit(redundant)$state,
    "not_evaluated_additive_stage"
  )

  invalid_observation <- problem
  invalid_observation$score_k[1] <- 9L
  expect_identical(
    call_lexicographic_limit(invalid_observation)$state,
    "not_evaluated_observation_map"
  )

  noninteger_score <- problem
  noninteger_score$score_k[1] <- 0.5
  expect_identical(
    call_lexicographic_limit(noninteger_score)$state,
    "not_evaluated_observation_map"
  )

  noninteger_slope <- problem
  noninteger_slope$slope_index[1] <- 1.5
  expect_identical(
    call_lexicographic_limit(noninteger_slope)$state,
    "not_evaluated_observation_map"
  )

  invalid_control <- call_lexicographic_limit(
    problem, max_categories = 1L
  )
  expect_identical(invalid_control$state, "not_evaluated_control")

  size_limited <- call_lexicographic_limit(
    problem, max_observations = 2L
  )
  expect_identical(size_limited$state, "not_evaluated_size_limit")

  invalid_direct <- call_lexicographic_at(invalid_rate, 16)
  expect_false(invalid_direct$valid)
  expect_identical(
    invalid_direct$reason_codes,
    "lexicographic_limit_rate_sum_zero_invalid"
  )
})

test_that("production fits advertise the oracle without inferring a path", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_lexicographic_limit_fixture()
  audit <- fit$config$boundary_audit$gpcm_lexicographic_limit

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-lexicographic-limit-0.2.3-v1"
  )
  expect_identical(
    audit$state,
    "declared_two_stage_limit_oracle_available_no_path_inferred"
  )
  expect_true(audit$evaluated)
  expect_true(audit$declared_two_stage_limit_oracle_available)
  expect_true(audit$exact_coefficient_comparisons_required)
  expect_identical(audit$maximum_slope_stages, 2L)
  expect_identical(audit$maximum_additive_stages, 2L)
  expect_true(audit$positive_power_scales_only)
  expect_true(audit$exponential_slope_vs_polynomial_additive_dominance_applied)
  expect_false(
    audit$declared_utility_path_parameter_reachability_checked
  )
  expect_false(audit$path_inferred_from_optimizer_trace)
  expect_false(audit$path_search_performed)
  expect_false(audit$production_path_limit_classified)
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
  expect_true(audit$source_rate_hierarchy_contract_matches)
  expect_match(audit$limitations, "does not derive a path", fixed = TRUE)
})

test_that("fit-level oracle scope does not borrow MML or FACETS authority", {
  fit <- fit_jml_gpcm_lexicographic_limit_fixture()
  config <- fit$config
  sizes <- mfrmr:::build_param_sizes(config)
  hierarchy <- config$boundary_audit$gpcm_rate_hierarchy

  mismatch <- hierarchy
  mismatch$contract_version <- "mismatched"
  structural <- mfrmr:::audit_mfrm_jml_gpcm_lexicographic_limit_scope(
    config, sizes, mismatch
  )
  expect_identical(
    structural$state,
    "declared_two_stage_limit_oracle_available_no_path_inferred"
  )
  expect_false(structural$source_rate_hierarchy_contract_matches)
  expect_false(structural$production_path_limit_classified)

  mml_config <- config
  mml_config$method <- "MML"
  mml <- mfrmr:::audit_mfrm_jml_gpcm_lexicographic_limit_scope(
    mml_config, sizes, hierarchy
  )
  expect_identical(mml$state, "not_applicable_estimator")
  expect_identical(mml$objective_identity, "not_applicable_marginal_estimator")
  expect_false(mml$external_comparison_eligible)

  pcm_config <- config
  pcm_config$model <- "PCM"
  pcm <- mfrmr:::audit_mfrm_jml_gpcm_lexicographic_limit_scope(
    pcm_config, sizes, hierarchy
  )
  expect_identical(pcm$state, "not_applicable_model")
  expect_identical(pcm$objective_identity, "not_applicable")

  unit_config <- config
  unit_config$gpcm_spec$levels <- "C1"
  unit_config$gpcm_spec$n_params <- 0L
  unit_sizes <- sizes
  unit_sizes$log_slopes <- 0L
  unit <- mfrmr:::audit_mfrm_jml_gpcm_lexicographic_limit_scope(
    unit_config, unit_sizes, hierarchy
  )
  expect_identical(unit$state, "not_required_unit_slope")
  expect_false(unit$declared_two_stage_limit_oracle_available)

  invalid_sizes <- sizes
  invalid_sizes$log_slopes <- 99L
  invalid <- mfrmr:::audit_mfrm_jml_gpcm_lexicographic_limit_scope(
    config, invalid_sizes, hierarchy
  )
  expect_identical(invalid$state, "not_evaluated_dimension")
  expect_false(invalid$evaluated)
})
