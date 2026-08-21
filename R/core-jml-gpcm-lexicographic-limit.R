# JML GPCM declared lexicographic boundary-limit oracle
# ==============================================================================
#
# This oracle evaluates a deliberately narrow path class. Expanded log slopes
# have at most two positive power scales, with an exact sum-zero coefficient
# vector at each scale. Unscaled cumulative category utilities have at most two
# positive power scales. The path is supplied explicitly; no asymptotic path is
# inferred from a finite optimizer trace.
#
# If a slope tends to infinity, probabilities concentrate on the exact
# lexicographic maxima of the additive directions followed by the base
# utility. If a slope tends to zero, exponential decay dominates every
# admitted polynomial additive scale and the response becomes uniform. If a
# slope stays finite, divergent additive directions select a tie set and the
# finite base slope produces a softmax over that set. These row limits yield an
# analytic declared-path log-likelihood limit, possibly negative infinity.

mfrmr_jml_gpcm_lexicographic_limit_contract_version <- function() {
  "mfrmr-jml-gpcm-lexicographic-limit-0.2.3-v1"
}

mfrmr_jml_gpcm_lexicographic_limit_empty_row_table <- function() {
  data.frame(
    Observation = integer(0),
    SlopeIndex = integer(0),
    SlopeLevel = character(0),
    ObservedCategory = integer(0),
    Weight = numeric(0),
    SlopeLimitState = character(0),
    SlopeResolutionStage = integer(0),
    SlopeResolutionCoefficient = numeric(0),
    AdditiveDivergentTieCategories = character(0),
    LimitSupportCategories = character(0),
    LimitSupportSize = integer(0),
    ObservedInLimitSupport = logical(0),
    RowLogProbabilityLimit = numeric(0),
    FiniteRowLimit = logical(0),
    EvaluationState = character(0),
    ReasonCodes = character(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_lexicographic_limit_empty_slope_table <- function() {
  data.frame(
    SlopeIndex = integer(0),
    SlopeLevel = character(0),
    BaseLogSlope = numeric(0),
    LimitState = character(0),
    ResolutionStage = integer(0),
    ResolutionCoefficient = numeric(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_lexicographic_top <- function(values, candidates) {
  candidates <- as.integer(candidates)
  if (length(candidates) == 0L) return(integer(0))
  selected <- as.numeric(values)[candidates]
  candidates[selected == max(selected)]
}

mfrmr_jml_gpcm_lexicographic_stage_matrices <- function(
    directions,
    n_observations,
    n_categories) {
  if (is.null(directions) || length(directions) == 0L) return(list())
  if (is.matrix(directions) || is.data.frame(directions)) {
    directions <- list(directions)
  }
  if (!is.list(directions)) return(NULL)
  out <- lapply(directions, function(direction) {
    converted <- tryCatch(
      matrix(
        as.numeric(direction),
        nrow = nrow(direction), ncol = ncol(direction),
        dimnames = dimnames(direction)
      ),
      error = function(e) NULL
    )
    if (is.null(converted) ||
        !identical(dim(converted), c(n_observations, n_categories)) ||
        any(!is.finite(converted))) {
      return(NULL)
    }
    converted
  })
  if (any(vapply(out, is.null, logical(1)))) return(NULL)
  out
}

mfrmr_jml_gpcm_lexicographic_rate_matrix <- function(rates, n_levels) {
  if (is.null(rates) || length(rates) == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = n_levels))
  }
  converted <- tryCatch({
    if (is.null(dim(rates))) {
      matrix(as.numeric(rates), nrow = 1L)
    } else {
      matrix(
        as.numeric(rates),
        nrow = nrow(rates), ncol = ncol(rates),
        dimnames = dimnames(rates)
      )
    }
  }, error = function(e) NULL)
  if (is.null(converted) || ncol(converted) != n_levels ||
      any(!is.finite(converted))) {
    return(NULL)
  }
  converted
}

mfrmr_jml_gpcm_lexicographic_scale_order_valid <- function(
    exponents,
    stages) {
  exponents <- tryCatch(
    suppressWarnings(as.numeric(exponents)),
    error = function(e) NULL
  )
  length(exponents) == stages && all(is.finite(exponents)) &&
    all(exponents > 0) &&
    (stages <= 1L || all(diff(exponents) < 0))
}

mfrmr_jml_gpcm_declared_hierarchy_limit <- function(
    base_utilities,
    score_k,
    slope_index,
    base_log_slopes,
    slope_stage_rates = NULL,
    slope_scale_exponents = numeric(0),
    additive_stage_directions = list(),
    additive_scale_exponents = numeric(0),
    weight = NULL,
    slope_levels = NULL,
    max_slope_stages = 2L,
    max_additive_stages = 2L,
    max_observations = 20000L,
    max_categories = 20L) {
  empty_rows <- mfrmr_jml_gpcm_lexicographic_limit_empty_row_table()
  empty_slopes <- mfrmr_jml_gpcm_lexicographic_limit_empty_slope_table()
  base <- tryCatch(
    matrix(
      as.numeric(base_utilities),
      nrow = nrow(base_utilities), ncol = ncol(base_utilities),
      dimnames = dimnames(base_utilities)
    ),
    error = function(e) NULL
  )
  n_observations <- if (is.null(base)) 0L else nrow(base)
  n_categories <- if (is.null(base)) 0L else ncol(base)
  base_logs <- tryCatch(
    suppressWarnings(as.numeric(base_log_slopes)),
    error = function(e) NULL
  )
  n_levels <- length(base_logs)
  scalar_integer_or_na <- function(value) {
    converted <- tryCatch(
      suppressWarnings(as.integer(value)[1]),
      error = function(e) NA_integer_
    )
    if (length(converted) == 0L) NA_integer_ else converted
  }

  finish <- function(state, evaluated, detail, reason_codes,
                     classified = FALSE,
                     slope_table = empty_slopes,
                     row_limits = empty_rows,
                     slope_hierarchy_state = "not_evaluated",
                     log_likelihood_limit = NA_real_,
                     path_definition = list()) {
    effective <- if (nrow(row_limits) > 0L) row_limits$Weight > 0 else logical(0)
    finite_effective <- length(effective) > 0L &&
      all(row_limits$FiniteRowLimit[effective])
    negative_infinite <- length(effective) > 0L && any(
      is.infinite(row_limits$RowLogProbabilityLimit[effective]) &
        row_limits$RowLogProbabilityLimit[effective] < 0
    )
    list(
      contract_version =
        mfrmr_jml_gpcm_lexicographic_limit_contract_version(),
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      path_family = paste(
        "declared_positive-power additive and sum-zero expanded-log-slope",
        "hierarchy with at most two stages per block"
      ),
      state = state,
      evaluated = isTRUE(evaluated),
      declared_path_limit_classified = isTRUE(classified),
      exact_coefficient_comparisons = TRUE,
      slope_stages = as.integer(nrow(
        path_definition$slope_stage_rates %||%
          matrix(numeric(0), nrow = 0L, ncol = n_levels)
      )),
      additive_stages = as.integer(length(
        path_definition$additive_stage_directions %||% list()
      )),
      maximum_slope_stages = scalar_integer_or_na(max_slope_stages),
      maximum_additive_stages = scalar_integer_or_na(max_additive_stages),
      observations = as.integer(n_observations),
      categories = as.integer(n_categories),
      slope_levels = as.integer(n_levels),
      slope_hierarchy_state = slope_hierarchy_state,
      slope_table = slope_table,
      row_limits = row_limits,
      log_likelihood_limit = as.numeric(log_likelihood_limit),
      finite_log_likelihood_limit = isTRUE(classified) &&
        isTRUE(finite_effective) && is.finite(log_likelihood_limit),
      negative_infinite_log_likelihood_limit = isTRUE(classified) &&
        isTRUE(negative_infinite) &&
        is.infinite(log_likelihood_limit) && log_likelihood_limit < 0,
      positive_infinite_log_likelihood_limit = FALSE,
      direct_finite_distance_evaluation_available = isTRUE(classified),
      declared_utility_path_parameter_reachability_checked = FALSE,
      path_inferred_from_optimizer_trace = FALSE,
      path_search_performed = FALSE,
      monotone_tail_certified = FALSE,
      competitive_boundary_certified = FALSE,
      common_subsequence_limit_certified = FALSE,
      arbitrary_path_limit_classified = FALSE,
      global_boundary_classified = FALSE,
      global_finite_maximum_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      path_definition = path_definition,
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The oracle classifies only the exact declared path. It admits at",
        "most two positive power scales in each parameter block, uses exact",
        "coefficient zero and tie identities, and assumes no unreported",
        "remainder changes those identities. Parameter-space reachability of",
        "the caller-supplied cumulative-utility directions is not checked. It",
        "neither infers a path from a",
        "finite optimizer trace nor searches for one. Monotonicity,",
        "competitiveness, common limits across subsequences, arbitrary",
        "remainders, global boundary status, finite maxima, MML, uncertainty,",
        "readiness, and external comparison remain unclassified."
      )
    )
  }

  controls <- tryCatch(
    suppressWarnings(as.numeric(c(
      max_slope_stages, max_additive_stages,
      max_observations, max_categories
    ))),
    error = function(e) NULL
  )
  valid_controls <- length(controls) == 4L && all(is.finite(controls)) &&
    all(controls >= 0) && all(controls == floor(controls)) &&
    max_observations >= 1L && max_categories >= 2L
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE,
      "Stage and workload limits must be finite nonnegative integers.",
      "lexicographic_limit_control_invalid"
    ))
  }
  max_slope_stages <- as.integer(controls[1])
  max_additive_stages <- as.integer(controls[2])
  max_observations <- as.integer(controls[3])
  max_categories <- as.integer(controls[4])
  valid_base <- !is.null(base) && n_observations >= 1L &&
    n_categories >= 2L && all(is.finite(base)) && n_levels >= 2L &&
    all(is.finite(base_logs))
  if (!valid_base) {
    return(finish(
      "not_evaluated_base", FALSE,
      "Base utilities and expanded base log slopes were malformed or non-finite.",
      "lexicographic_limit_base_invalid"
    ))
  }
  if (n_observations > max_observations || n_categories > max_categories) {
    return(finish(
      "not_evaluated_size_limit", FALSE,
      "The declared path exceeded the observation or category workload limit.",
      "lexicographic_limit_size_limit"
    ))
  }
  if (!identical(sum(base_logs), 0)) {
    return(finish(
      "not_evaluated_identification", FALSE,
      "Expanded base log slopes must sum exactly to zero.",
      "lexicographic_limit_base_log_slope_sum_zero_invalid"
    ))
  }

  if (is.null(slope_levels)) {
    slope_levels <- names(base_log_slopes) %||% paste0("L", seq_len(n_levels))
  }
  slope_levels <- as.character(slope_levels)
  valid_levels <- length(slope_levels) == n_levels &&
    !anyNA(slope_levels) && all(nzchar(slope_levels)) &&
    !anyDuplicated(slope_levels)
  if (!valid_levels) {
    return(finish(
      "not_evaluated_levels", FALSE,
      "Slope level labels must be unique, nonempty, and dimension-matched.",
      "lexicographic_limit_slope_levels_invalid"
    ))
  }

  rates <- mfrmr_jml_gpcm_lexicographic_rate_matrix(
    slope_stage_rates, n_levels
  )
  directions <- mfrmr_jml_gpcm_lexicographic_stage_matrices(
    additive_stage_directions, n_observations, n_categories
  )
  if (is.null(rates) || is.null(directions)) {
    return(finish(
      "not_evaluated_stage_mapping", FALSE,
      "Slope-rate or additive-direction stage dimensions were malformed.",
      "lexicographic_limit_stage_mapping_invalid"
    ))
  }
  n_slope_stages <- nrow(rates)
  n_additive_stages <- length(directions)
  if (n_slope_stages > max_slope_stages ||
      n_additive_stages > max_additive_stages) {
    return(finish(
      "not_evaluated_stage_limit", FALSE,
      "The declared hierarchy exceeded the two-stage oracle scope.",
      "lexicographic_limit_stage_limit"
    ))
  }
  slope_scale_valid <- mfrmr_jml_gpcm_lexicographic_scale_order_valid(
    slope_scale_exponents, n_slope_stages
  )
  additive_scale_valid <- mfrmr_jml_gpcm_lexicographic_scale_order_valid(
    additive_scale_exponents, n_additive_stages
  )
  if (!slope_scale_valid || !additive_scale_valid) {
    return(finish(
      "not_evaluated_scale_order", FALSE,
      "Each declared exponent vector must be positive, finite, dimension-matched, and strictly decreasing.",
      "lexicographic_limit_scale_order_invalid"
    ))
  }
  if (n_slope_stages > 0L &&
      (!all(rowSums(rates) == 0) ||
       any(apply(abs(rates), 1L, max) == 0))) {
    return(finish(
      "not_evaluated_identification", FALSE,
      "Every slope-rate stage must be nonzero and sum exactly to zero.",
      "lexicographic_limit_rate_sum_zero_invalid"
    ))
  }
  if (n_additive_stages > 0L) {
    contrast_nonzero <- vapply(directions, function(direction) {
      any(apply(direction, 1L, function(row) max(row) > min(row)))
    }, logical(1))
    if (!all(contrast_nonzero)) {
      return(finish(
        "not_evaluated_additive_stage", FALSE,
        "Every additive stage must change at least one within-row category contrast.",
        "lexicographic_limit_additive_stage_redundant"
      ))
    }
  }

  score_values <- tryCatch(
    suppressWarnings(as.numeric(score_k)),
    error = function(e) NULL
  )
  slope_values <- tryCatch(
    suppressWarnings(as.numeric(slope_index)),
    error = function(e) NULL
  )
  if (is.null(weight)) weight <- rep(1, n_observations)
  weights <- tryCatch(
    suppressWarnings(as.numeric(weight)),
    error = function(e) NULL
  )
  valid_observation_map <- length(score_values) == n_observations &&
    length(slope_values) == n_observations &&
    length(weights) == n_observations &&
    !anyNA(score_values) && !anyNA(slope_values) && !anyNA(weights) &&
    all(is.finite(score_values)) && all(is.finite(slope_values)) &&
    all(score_values == floor(score_values)) &&
    all(slope_values == floor(slope_values)) &&
    all(score_values >= 0 & score_values < n_categories) &&
    all(slope_values >= 1 & slope_values <= n_levels) &&
    all(is.finite(weights)) && all(weights >= 0) && any(weights > 0)
  if (!valid_observation_map) {
    return(finish(
      "not_evaluated_observation_map", FALSE,
      "Scores, slope indices, or likelihood weights were malformed.",
      "lexicographic_limit_observation_map_invalid"
    ))
  }
  scores <- as.integer(score_values)
  slope_map <- as.integer(slope_values)

  path_definition <- list(
    base_utilities = base,
    score_k = scores,
    slope_index = slope_map,
    base_log_slopes = base_logs,
    slope_stage_rates = rates,
    slope_scale_exponents = as.numeric(slope_scale_exponents),
    additive_stage_directions = directions,
    additive_scale_exponents = as.numeric(additive_scale_exponents),
    weight = weights,
    slope_levels = slope_levels
  )

  rate_hierarchy <- if (n_slope_stages > 0L) {
    mfrmr_jml_gpcm_rate_hierarchy(
      rates,
      levels = slope_levels,
      scale_exponents = slope_scale_exponents,
      zero_tolerance = 0
    )
  } else {
    list(state = "no_divergent_log_slope_stage")
  }
  rate_normalized <- if (n_slope_stages > 0L) {
    rates / apply(abs(rates), 1L, max)
  } else {
    rates
  }
  slope_state <- rep("finite", n_levels)
  resolution_stage <- rep(NA_integer_, n_levels)
  resolution_coefficient <- rep(NA_real_, n_levels)
  for (level in seq_len(n_levels)) {
    if (n_slope_stages == 0L) next
    nonzero <- which(rate_normalized[, level] != 0)
    if (length(nonzero) == 0L) next
    first <- nonzero[1]
    coefficient <- rate_normalized[first, level]
    resolution_stage[level] <- first
    resolution_coefficient[level] <- coefficient
    slope_state[level] <- if (coefficient > 0) "infinite" else "zero"
  }
  slope_table <- data.frame(
    SlopeIndex = seq_len(n_levels),
    SlopeLevel = slope_levels,
    BaseLogSlope = base_logs,
    LimitState = slope_state,
    ResolutionStage = resolution_stage,
    ResolutionCoefficient = resolution_coefficient,
    stringsAsFactors = FALSE
  )

  category_label <- function(indices) {
    paste(as.integer(indices) - 1L, collapse = ";")
  }
  row_results <- lapply(seq_len(n_observations), function(observation) {
    candidates <- seq_len(n_categories)
    if (n_additive_stages > 0L) {
      for (stage in seq_len(n_additive_stages)) {
        candidates <- mfrmr_jml_gpcm_lexicographic_top(
          directions[[stage]][observation, ], candidates
        )
      }
    }
    divergent_ties <- candidates
    group <- slope_map[observation]
    observed <- scores[observation] + 1L
    state <- slope_state[group]
    if (identical(state, "zero")) {
      support <- seq_len(n_categories)
      row_limit <- -log(n_categories)
      evaluation <- "uniform_vanishing_slope"
      reason <- "negative_log_slope_hierarchy_dominates_polynomial_additive_path"
    } else if (identical(state, "infinite")) {
      support <- mfrmr_jml_gpcm_lexicographic_top(
        base[observation, ], divergent_ties
      )
      row_limit <- if (observed %in% support) {
        -log(length(support))
      } else {
        -Inf
      }
      evaluation <- "infinite_slope_lexicographic_concentration"
      reason <- if (observed %in% support) {
        "observed_category_in_infinite_slope_limit_support"
      } else {
        "observed_category_excluded_from_infinite_slope_limit_support"
      }
    } else {
      support <- divergent_ties
      if (!observed %in% support) {
        row_limit <- -Inf
      } else {
        scaled <- exp(base_logs[group]) * base[observation, support]
        maximum <- max(scaled)
        denominator <- maximum + log(sum(exp(scaled - maximum)))
        row_limit <-
          exp(base_logs[group]) * base[observation, observed] - denominator
      }
      evaluation <- "finite_slope_softmax_on_divergent_tie_set"
      reason <- if (observed %in% support) {
        "observed_category_in_finite_slope_tie_set"
      } else {
        "observed_category_excluded_by_additive_hierarchy"
      }
    }
    data.frame(
      Observation = as.integer(observation),
      SlopeIndex = as.integer(group),
      SlopeLevel = slope_levels[group],
      ObservedCategory = as.integer(scores[observation]),
      Weight = weights[observation],
      SlopeLimitState = state,
      SlopeResolutionStage = resolution_stage[group],
      SlopeResolutionCoefficient = resolution_coefficient[group],
      AdditiveDivergentTieCategories = category_label(divergent_ties),
      LimitSupportCategories = category_label(support),
      LimitSupportSize = as.integer(length(support)),
      ObservedInLimitSupport = observed %in% support,
      RowLogProbabilityLimit = as.numeric(row_limit),
      FiniteRowLimit = is.finite(row_limit),
      EvaluationState = evaluation,
      ReasonCodes = reason,
      stringsAsFactors = FALSE
    )
  })
  row_limits <- do.call(rbind, row_results)
  effective <- weights > 0
  log_likelihood_limit <- if (any(
    is.infinite(row_limits$RowLogProbabilityLimit[effective]) &
      row_limits$RowLogProbabilityLimit[effective] < 0
  )) {
    -Inf
  } else {
    sum(weights[effective] *
          row_limits$RowLogProbabilityLimit[effective])
  }
  state <- if (is.finite(log_likelihood_limit)) {
    "finite_declared_hierarchy_log_likelihood_limit"
  } else {
    "negative_infinite_declared_hierarchy_log_likelihood_limit"
  }
  finish(
    state, TRUE,
    paste(
      "The exact declared two-stage-or-smaller power hierarchy has an",
      "analytic response-level and joint log-likelihood limit."
    ),
    if (is.finite(log_likelihood_limit)) {
      "declared_hierarchy_finite_likelihood_limit"
    } else {
      "declared_hierarchy_negative_infinite_likelihood_limit"
    },
    classified = TRUE,
    slope_table = slope_table,
    row_limits = row_limits,
    slope_hierarchy_state = as.character(rate_hierarchy$state)[1],
    log_likelihood_limit = log_likelihood_limit,
    path_definition = path_definition
  )
}

mfrmr_jml_gpcm_declared_hierarchy_loglik_at <- function(
    distance,
    base_utilities,
    score_k,
    slope_index,
    base_log_slopes,
    slope_stage_rates = NULL,
    slope_scale_exponents = numeric(0),
    additive_stage_directions = list(),
    additive_scale_exponents = numeric(0),
    weight = NULL,
    slope_levels = NULL,
    max_slope_stages = 2L,
    max_additive_stages = 2L,
    max_observations = 20000L,
    max_categories = 20L) {
  distance <- tryCatch(
    suppressWarnings(as.numeric(distance)[1]),
    error = function(e) NA_real_
  )
  if (!is.finite(distance) || distance <= 0) {
    return(list(
      valid = FALSE,
      distance = distance,
      log_likelihood = NA_real_,
      row_log_probability = numeric(0),
      reason_codes = "lexicographic_limit_distance_invalid"
    ))
  }
  limit <- mfrmr_jml_gpcm_declared_hierarchy_limit(
    base_utilities = base_utilities,
    score_k = score_k,
    slope_index = slope_index,
    base_log_slopes = base_log_slopes,
    slope_stage_rates = slope_stage_rates,
    slope_scale_exponents = slope_scale_exponents,
    additive_stage_directions = additive_stage_directions,
    additive_scale_exponents = additive_scale_exponents,
    weight = weight,
    slope_levels = slope_levels,
    max_slope_stages = max_slope_stages,
    max_additive_stages = max_additive_stages,
    max_observations = max_observations,
    max_categories = max_categories
  )
  if (!isTRUE(limit$declared_path_limit_classified)) {
    return(list(
      valid = FALSE,
      distance = distance,
      log_likelihood = NA_real_,
      row_log_probability = numeric(0),
      reason_codes = limit$reason_codes
    ))
  }
  path <- limit$path_definition
  log_slopes <- path$base_log_slopes
  if (nrow(path$slope_stage_rates) > 0L) {
    for (stage in seq_len(nrow(path$slope_stage_rates))) {
      log_slopes <- log_slopes +
        distance^path$slope_scale_exponents[stage] *
        path$slope_stage_rates[stage, ]
    }
  }
  utilities <- path$base_utilities
  if (length(path$additive_stage_directions) > 0L) {
    for (stage in seq_along(path$additive_stage_directions)) {
      utilities <- utilities +
        distance^path$additive_scale_exponents[stage] *
        path$additive_stage_directions[[stage]]
    }
  }
  slopes <- exp(log_slopes)
  row_log_probability <- rep(NA_real_, nrow(utilities))
  for (observation in seq_len(nrow(utilities))) {
    group <- path$slope_index[observation]
    scaled <- slopes[group] * utilities[observation, ]
    if (any(!is.finite(scaled))) next
    maximum <- max(scaled)
    denominator <- maximum + log(sum(exp(scaled - maximum)))
    observed <- path$score_k[observation] + 1L
    row_log_probability[observation] <- scaled[observed] - denominator
  }
  effective <- path$weight > 0
  valid <- all(is.finite(row_log_probability[effective]))
  list(
    valid = isTRUE(valid),
    distance = distance,
    log_likelihood = if (valid) {
      sum(path$weight[effective] * row_log_probability[effective])
    } else {
      NA_real_
    },
    row_log_probability = row_log_probability,
    reason_codes = if (valid) {
      "declared_hierarchy_finite_distance_evaluated"
    } else {
      "declared_hierarchy_finite_distance_overflow_or_nonfinite"
    }
  )
}

audit_mfrm_jml_gpcm_lexicographic_limit_scope <- function(
    config,
    sizes,
    rate_hierarchy =
      config$boundary_audit$gpcm_rate_hierarchy %||% list()) {
  method <- as.character(config$method %||% "unknown")[1]
  model <- as.character(config$model %||% "unknown")[1]
  levels <- as.character(config$gpcm_spec$levels %||% character(0))
  n_levels <- length(levels)
  free_log_slopes <- as.integer(config$gpcm_spec$n_params %||% 0L)
  size_values <- suppressWarnings(vapply(sizes, as.integer, integer(1)))
  size_log_slopes <- if ("log_slopes" %in% names(size_values)) {
    as.integer(size_values[["log_slopes"]])
  } else {
    NA_integer_
  }
  dimension_valid <- length(size_values) > 0L &&
    all(is.finite(size_values)) && all(size_values >= 0L) &&
    is.finite(size_log_slopes) && is.finite(free_log_slopes) &&
    size_log_slopes == free_log_slopes &&
    free_log_slopes == max(n_levels - 1L, 0L)
  source_matches <- identical(
    as.character(rate_hierarchy$contract_version %||% "")[1],
    mfrmr_jml_gpcm_rate_hierarchy_contract_version()
  ) && identical(
    as.character(rate_hierarchy$objective_identity %||% "")[1],
    "identified_conditional_joint_log_likelihood"
  )

  finish <- function(state, evaluated, detail, reason_codes,
                     available = FALSE) {
    list(
      contract_version =
        mfrmr_jml_gpcm_lexicographic_limit_contract_version(),
      method = method,
      model = model,
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      state = state,
      evaluated = isTRUE(evaluated),
      declared_two_stage_limit_oracle_available = isTRUE(available),
      exact_coefficient_comparisons_required = TRUE,
      maximum_slope_stages = 2L,
      maximum_additive_stages = 2L,
      positive_power_scales_only = TRUE,
      exponential_slope_vs_polynomial_additive_dominance_applied =
        isTRUE(available),
      declared_utility_path_parameter_reachability_checked = FALSE,
      path_inferred_from_optimizer_trace = FALSE,
      path_search_performed = FALSE,
      production_path_limit_classified = FALSE,
      arbitrary_path_limit_classified = FALSE,
      common_subsequence_limit_certified = FALSE,
      monotone_tail_certified = FALSE,
      competitive_boundary_certified = FALSE,
      global_boundary_classified = FALSE,
      global_finite_maximum_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      source_rate_hierarchy_contract = as.character(
        rate_hierarchy$contract_version %||% "missing"
      )[1],
      source_rate_hierarchy_state = as.character(
        rate_hierarchy$state %||% "missing"
      )[1],
      source_rate_hierarchy_contract_matches = isTRUE(source_matches),
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The fit-level object advertises an internal declared-path oracle; it",
        "does not derive a path from the fit or classify the production",
        "likelihood boundary. The v1 oracle is limited to exact coefficients",
        "and two positive power scales per block. Parameter-space reachability",
        "of declared cumulative-utility directions is not checked. Remainders,",
        "broader scale families, path search, monotonicity, competitiveness, common",
        "subsequence limits, global status, MML, uncertainty, readiness, and",
        "external comparison remain open."
      )
    )
  }

  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "The lexicographic JML limit oracle applies only to GPCM.",
      "lexicographic_limit_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional fixed-effects limit oracle is not reused for MML.",
      "lexicographic_limit_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (!dimension_valid) {
    return(finish(
      "not_evaluated_dimension", FALSE,
      "The identified free log-slope dimension did not match the slope levels.",
      "lexicographic_limit_dimension_invalid"
    ))
  }
  if (free_log_slopes == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction uses the existing PCM boundary geometry.",
      "no_free_log_slope_coordinate"
    ))
  }
  finish(
    "declared_two_stage_limit_oracle_available_no_path_inferred", TRUE,
    paste(
      "An exact declared-path oracle is available for at most two positive",
      "power scales in the slope and additive blocks. No production path was",
      "inferred or classified."
    ),
    paste(
      "declared_two_stage_limit_oracle_available",
      "production_path_not_inferred",
      "global_boundary_not_classified",
      sep = ";"
    ),
    available = TRUE
  )
}
