# JML GPCM parameter-sequence contrast-flag compactification
# ==============================================================================
#
# A free-parameter sequence reaches response probabilities through a nonlinear
# map: cumulative utilities are linear in the additive coordinates, expanded
# log slopes are linear in their identified free coordinates, and GPCM logits
# multiply the utilities by exponentiated log slopes. Trying to transport an
# additive/log-slope Euclidean direction through that product can therefore
# miss cancellations and intermediate scales.
#
# P2e applies compactness after the exact nonlinear map. Probabilities depend
# only on category logits relative to a reference category, so a retained fit
# has a finite-dimensional contrast image with N(K-1) coordinates. Every
# sequence in that space has a further subsequence which is either convergent
# or admits a finite scale-separated flag, followed by a convergent base and
# an o(1) contrast remainder. Exact lexicographic maximization of the flag and
# a base softmax classify the response and joint likelihood limit along that
# further subsequence.
#
# This is not an optimizer-trace extractor and not a global-boundary theorem.
# Different subsequences may have different flags and different limits.

mfrmr_jml_gpcm_parameter_sequence_flag_contract_version <- function() {
  "mfrmr-jml-gpcm-parameter-sequence-flag-0.2.3-v1"
}

mfrmr_jml_gpcm_parameter_sequence_flag_empty_rows <- function() {
  data.frame(
    Observation = integer(0),
    ObservedCategory = integer(0),
    Weight = numeric(0),
    DivergentStagesUsed = integer(0),
    ObservedExclusionStage = integer(0),
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

mfrmr_jml_gpcm_parameter_sequence_flag_empty_stages <- function() {
  data.frame(
    Stage = integer(0),
    MaximumAbsoluteCoefficient = numeric(0),
    RowsResolved = integer(0),
    CategoriesRemoved = integer(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_reference_contrasts <- function(logits) {
  value <- tryCatch(
    matrix(
      as.numeric(logits),
      nrow = nrow(logits), ncol = ncol(logits),
      dimnames = dimnames(logits)
    ),
    error = function(e) NULL
  )
  if (is.null(value) || nrow(value) < 1L || ncol(value) < 2L ||
      any(!is.finite(value))) {
    return(NULL)
  }
  sweep(
    value[, -1L, drop = FALSE],
    1L, value[, 1L], FUN = "-"
  )
}

mfrmr_jml_gpcm_contrast_flag_matrices <- function(
    directions,
    n_observations,
    n_contrasts) {
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
    if (is.null(converted) || nrow(converted) != n_observations ||
        ncol(converted) != n_contrasts || any(!is.finite(converted))) {
      return(NULL)
    }
    converted
  })
  if (any(vapply(out, is.null, logical(1)))) return(NULL)
  out
}

mfrmr_jml_gpcm_declared_contrast_flag_limit <- function(
    base_contrasts,
    score_k,
    stage_directions = list(),
    weight = NULL,
    max_stages = 100L,
    max_observations = 20000L,
    max_categories = 20L,
    max_contrast_elements = 2e6) {
  empty_rows <- mfrmr_jml_gpcm_parameter_sequence_flag_empty_rows()
  empty_stages <- mfrmr_jml_gpcm_parameter_sequence_flag_empty_stages()
  base <- tryCatch(
    matrix(
      as.numeric(base_contrasts),
      nrow = nrow(base_contrasts), ncol = ncol(base_contrasts),
      dimnames = dimnames(base_contrasts)
    ),
    error = function(e) NULL
  )
  n_observations <- if (is.null(base)) 0L else nrow(base)
  n_contrasts <- if (is.null(base)) 0L else ncol(base)
  n_categories <- n_contrasts + 1L
  contrast_dimension <- as.double(n_observations) * as.double(n_contrasts)

  scalar_integer <- function(value) {
    converted <- tryCatch(
      suppressWarnings(as.numeric(value)[1]),
      error = function(e) NA_real_
    )
    if (length(converted) == 0L) NA_real_ else converted
  }
  controls <- c(
    scalar_integer(max_stages),
    scalar_integer(max_observations),
    scalar_integer(max_categories),
    scalar_integer(max_contrast_elements)
  )

  finish <- function(state, evaluated, classified, detail, reason_codes,
                     row_limits = empty_rows,
                     stage_table = empty_stages,
                     log_likelihood_limit = NA_real_,
                     directions = list()) {
    effective <- if (nrow(row_limits) > 0L) {
      row_limits$Weight > 0
    } else {
      logical(0)
    }
    negative_infinite <- length(effective) > 0L && any(
      is.infinite(row_limits$RowLogProbabilityLimit[effective]) &
        row_limits$RowLogProbabilityLimit[effective] < 0
    )
    list(
      contract_version =
        mfrmr_jml_gpcm_parameter_sequence_flag_contract_version(),
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      response_coordinate_system =
        "category_logit_contrasts_relative_to_category_zero",
      path_family = paste(
        "declared finite scale-separated contrast flag with divergent",
        "positive scales, successive scale ratios tending to zero, a finite",
        "base, and a uniform o(1) reference-contrast remainder"
      ),
      state = state,
      evaluated = isTRUE(evaluated),
      declared_contrast_flag_limit_classified = isTRUE(classified),
      observations = as.integer(n_observations),
      categories = as.integer(n_categories),
      contrast_dimension = as.double(contrast_dimension),
      divergent_flag_stages = as.integer(length(directions)),
      maximum_flag_stages = suppressWarnings(as.integer(controls[1])),
      theoretical_flag_depth_bound = as.double(contrast_dimension),
      scale_contract =
        "s_g_to_infinity_and_s_(g+1)_over_s_g_to_zero",
      contrast_remainder_contract = "uniform_reference_contrast_o(1)",
      exact_coefficient_comparisons = TRUE,
      exact_lexicographic_ties = TRUE,
      row_limits = row_limits,
      stage_table = stage_table,
      base_contrasts = base,
      stage_directions = directions,
      log_likelihood_limit = as.numeric(log_likelihood_limit),
      finite_log_likelihood_limit = isTRUE(classified) &&
        is.finite(log_likelihood_limit),
      negative_infinite_log_likelihood_limit = isTRUE(classified) &&
        isTRUE(negative_infinite) &&
        is.infinite(log_likelihood_limit) && log_likelihood_limit < 0,
      positive_infinite_log_likelihood_limit = FALSE,
      direct_finite_scale_evaluation_available = isTRUE(classified),
      parameter_sequence_reachability_checked = FALSE,
      parameter_flag_extracted = FALSE,
      optimizer_sequence_used = FALSE,
      p2b_handoff_eligible = FALSE,
      original_sequence_limit_classified = FALSE,
      common_subsequence_limit_certified = FALSE,
      arbitrary_parameter_sequence_global_limit_classified = FALSE,
      competitive_boundary_certified = FALSE,
      global_boundary_classified = FALSE,
      global_finite_maximum_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The oracle classifies only the declared contrast flag. The abstract",
        "subsequence theorem guarantees that an exact nonlinear contrast",
        "image has such a further-subsequence representation, but this call",
        "does not prove that caller-supplied directions came from a parameter",
        "sequence. It does not infer a flag from finite optimizer output,",
        "establish a common limit across subsequences, prove competitiveness,",
        "or classify a global boundary, finite maximum, MML, uncertainty,",
        "readiness, recovery, simulation, or external comparison."
      )
    )
  }

  valid_controls <- length(controls) == 4L && all(is.finite(controls)) &&
    all(controls == floor(controls)) && controls[1] >= 0 &&
    controls[2] >= 1 && controls[3] >= 2 && controls[4] >= 1
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "Flag and workload controls must be finite integers in their valid ranges.",
      "parameter_sequence_flag_control_invalid"
    ))
  }
  if (is.null(base) || n_observations < 1L || n_contrasts < 1L ||
      any(!is.finite(base))) {
    return(finish(
      "not_evaluated_base", FALSE, FALSE,
      "The finite base reference-contrast matrix was malformed.",
      "parameter_sequence_flag_base_invalid"
    ))
  }
  if (n_observations > controls[2] || n_categories > controls[3] ||
      contrast_dimension > controls[4]) {
    return(finish(
      "not_evaluated_size_limit", FALSE, FALSE,
      "The declared contrast flag exceeded a workload limit.",
      "parameter_sequence_flag_size_limit"
    ))
  }

  directions <- mfrmr_jml_gpcm_contrast_flag_matrices(
    stage_directions, n_observations, n_contrasts
  )
  if (is.null(directions)) {
    return(finish(
      "not_evaluated_stage_mapping", FALSE, FALSE,
      "A contrast-flag stage was non-finite or dimensionally inconsistent.",
      "parameter_sequence_flag_stage_mapping_invalid"
    ))
  }
  n_stages <- length(directions)
  if (n_stages > controls[1] || n_stages > contrast_dimension) {
    return(finish(
      "not_evaluated_stage_limit", FALSE, FALSE,
      "The declared flag exceeded its computational or finite-dimensional depth bound.",
      "parameter_sequence_flag_stage_limit"
    ))
  }
  if (contrast_dimension * (n_stages + 1) > controls[4]) {
    return(finish(
      "not_evaluated_size_limit", FALSE, FALSE,
      "The base plus declared flag stages exceeded the element workload limit.",
      "parameter_sequence_flag_size_limit"
    ))
  }
  if (n_stages > 0L && any(vapply(directions, function(direction) {
    max(abs(direction)) == 0
  }, logical(1)))) {
    return(finish(
      "not_evaluated_zero_stage", FALSE, FALSE,
      "Every declared divergent contrast stage must be globally nonzero.",
      "parameter_sequence_flag_zero_stage"
    ))
  }

  scores <- tryCatch(
    suppressWarnings(as.numeric(score_k)),
    error = function(e) NULL
  )
  if (is.null(weight)) weight <- rep(1, n_observations)
  weights <- tryCatch(
    suppressWarnings(as.numeric(weight)),
    error = function(e) NULL
  )
  observation_map_valid <- !is.null(scores) && !is.null(weights) &&
    length(scores) == n_observations && length(weights) == n_observations &&
    all(is.finite(scores)) && all(scores == floor(scores)) &&
    all(scores >= 0 & scores < n_categories) &&
    all(is.finite(weights)) && all(weights >= 0) && any(weights > 0)
  if (!observation_map_valid) {
    return(finish(
      "not_evaluated_observation_map", FALSE, FALSE,
      "Scores or fixed likelihood weights were malformed.",
      "parameter_sequence_flag_observation_map_invalid"
    ))
  }
  scores <- as.integer(scores)

  stage_rows_resolved <- integer(n_stages)
  stage_categories_removed <- integer(n_stages)
  category_label <- function(indices) {
    paste(as.integer(indices) - 1L, collapse = ";")
  }
  row_results <- lapply(seq_len(n_observations), function(observation) {
    candidates <- seq_len(n_categories)
    observed <- scores[observation] + 1L
    exclusion_stage <- NA_integer_
    stages_used <- 0L
    if (n_stages > 0L) {
      for (stage in seq_len(n_stages)) {
        full_direction <- c(0, directions[[stage]][observation, ])
        next_candidates <- mfrmr_jml_gpcm_lexicographic_top(
          full_direction, candidates
        )
        removed <- length(candidates) - length(next_candidates)
        if (removed > 0L) {
          stage_rows_resolved[stage] <<- stage_rows_resolved[stage] + 1L
          stage_categories_removed[stage] <<-
            stage_categories_removed[stage] + removed
        }
        candidates <- next_candidates
        stages_used <- stage
        if (is.na(exclusion_stage) && !observed %in% candidates) {
          exclusion_stage <- stage
        }
      }
    }
    base_logits <- c(0, base[observation, ])
    observed_in_support <- observed %in% candidates
    row_limit <- if (!observed_in_support) {
      -Inf
    } else {
      denominator <- logsumexp(base_logits[candidates])
      base_logits[observed] - denominator
    }
    data.frame(
      Observation = as.integer(observation),
      ObservedCategory = as.integer(scores[observation]),
      Weight = as.numeric(weights[observation]),
      DivergentStagesUsed = as.integer(stages_used),
      ObservedExclusionStage = as.integer(exclusion_stage),
      LimitSupportCategories = category_label(candidates),
      LimitSupportSize = as.integer(length(candidates)),
      ObservedInLimitSupport = observed_in_support,
      RowLogProbabilityLimit = as.numeric(row_limit),
      FiniteRowLimit = is.finite(row_limit),
      EvaluationState = if (observed_in_support) {
        "base_softmax_on_final_lexicographic_tie_set"
      } else {
        "observed_category_excluded_by_divergent_flag"
      },
      ReasonCodes = if (observed_in_support) {
        "observed_category_survives_all_divergent_contrast_stages"
      } else {
        "observed_category_removed_at_faster_contrast_scale"
      },
      stringsAsFactors = FALSE
    )
  })
  row_limits <- do.call(rbind, row_results)
  stage_table <- if (n_stages > 0L) {
    data.frame(
      Stage = seq_len(n_stages),
      MaximumAbsoluteCoefficient = vapply(
        directions, function(direction) max(abs(direction)), numeric(1)
      ),
      RowsResolved = stage_rows_resolved,
      CategoriesRemoved = stage_categories_removed,
      stringsAsFactors = FALSE
    )
  } else {
    empty_stages
  }
  effective <- weights > 0
  log_likelihood_limit <- if (any(
    is.infinite(row_limits$RowLogProbabilityLimit[effective]) &
      row_limits$RowLogProbabilityLimit[effective] < 0
  )) {
    -Inf
  } else {
    sum(weights[effective] * row_limits$RowLogProbabilityLimit[effective])
  }
  state <- if (is.finite(log_likelihood_limit)) {
    "finite_declared_contrast_flag_log_likelihood_limit"
  } else {
    "negative_infinite_declared_contrast_flag_log_likelihood_limit"
  }
  finish(
    state, TRUE, TRUE,
    paste(
      "Exact lexicographic maximization of the scale-separated reference-",
      "contrast flag followed by the finite base softmax classified every",
      "retained response-row limit."
    ),
    if (is.finite(log_likelihood_limit)) {
      "declared_contrast_flag_finite_likelihood_limit"
    } else {
      "declared_contrast_flag_negative_infinite_likelihood_limit"
    },
    row_limits = row_limits,
    stage_table = stage_table,
    log_likelihood_limit = log_likelihood_limit,
    directions = directions
  )
}

mfrmr_jml_gpcm_declared_contrast_flag_loglik_at <- function(
    result,
    scales,
    contrast_remainder = NULL) {
  valid_result <- is.list(result) && identical(
    as.character(result$contract_version %||% "")[1],
    mfrmr_jml_gpcm_parameter_sequence_flag_contract_version()
  ) && isTRUE(result$declared_contrast_flag_limit_classified)
  scale_values <- tryCatch(
    suppressWarnings(as.numeric(scales)),
    error = function(e) NULL
  )
  n_stages <- if (valid_result) length(result$stage_directions) else 0L
  valid_scales <- !is.null(scale_values) && length(scale_values) == n_stages &&
    all(is.finite(scale_values)) && all(scale_values > 0) &&
    (n_stages <= 1L || all(diff(scale_values) < 0))
  if (!valid_result || !valid_scales) {
    return(list(
      valid = FALSE,
      scales = scale_values,
      log_likelihood = NA_real_,
      row_log_probability = numeric(0),
      reason_codes = "declared_contrast_flag_evaluation_invalid"
    ))
  }
  contrasts <- result$base_contrasts
  if (n_stages > 0L) {
    for (stage in seq_len(n_stages)) {
      contrasts <- contrasts +
        scale_values[stage] * result$stage_directions[[stage]]
    }
  }
  if (is.null(contrast_remainder)) {
    remainder <- matrix(0, nrow(contrasts), ncol(contrasts))
  } else {
    remainder <- tryCatch(
      matrix(
        as.numeric(contrast_remainder),
        nrow = nrow(contrast_remainder), ncol = ncol(contrast_remainder)
      ),
      error = function(e) NULL
    )
  }
  if (is.null(remainder) || !identical(dim(remainder), dim(contrasts)) ||
      any(!is.finite(remainder))) {
    return(list(
      valid = FALSE,
      scales = scale_values,
      log_likelihood = NA_real_,
      row_log_probability = numeric(0),
      reason_codes = "declared_contrast_flag_remainder_invalid"
    ))
  }
  contrasts <- contrasts + remainder
  if (any(!is.finite(contrasts))) {
    return(list(
      valid = FALSE,
      scales = scale_values,
      log_likelihood = NA_real_,
      row_log_probability = numeric(0),
      reason_codes = "declared_contrast_flag_nonfinite_logits"
    ))
  }
  row_log_probability <- vapply(seq_len(nrow(contrasts)), function(observation) {
    logits <- c(0, contrasts[observation, ])
    observed <- result$row_limits$ObservedCategory[observation] + 1L
    logits[observed] - logsumexp(logits)
  }, numeric(1))
  effective <- result$row_limits$Weight > 0
  valid_numeric <- all(is.finite(row_log_probability[effective]))
  list(
    valid = isTRUE(valid_numeric),
    scales = scale_values,
    log_likelihood = if (valid_numeric) {
      sum(
        result$row_limits$Weight[effective] *
          row_log_probability[effective]
      )
    } else {
      NA_real_
    },
    row_log_probability = row_log_probability,
    maximum_absolute_contrast_remainder = max(abs(remainder)),
    reason_codes = if (valid_numeric) {
      "declared_contrast_flag_finite_scales_evaluated"
    } else {
      "declared_contrast_flag_finite_scale_probability_nonfinite"
    }
  )
}

mfrmr_jml_gpcm_parameter_sequence_contrast_image <- function(
    fit,
    parameter_sequence,
    max_points = 100L,
    max_contrast_elements = 2e6) {
  finish <- function(state, evaluated, detail, reason_codes,
                     contrast_sequence = matrix(numeric(0), 0L, 0L),
                     coordinate_map = data.frame()) {
    list(
      contract_version =
        mfrmr_jml_gpcm_parameter_sequence_flag_contract_version(),
      state = state,
      evaluated = isTRUE(evaluated),
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      source_sequence_declared_by_caller = TRUE,
      exact_nonlinear_parameter_to_contrast_map_evaluated = isTRUE(evaluated),
      slope_utility_product_evaluated_before_compactification = isTRUE(evaluated),
      parameter_sequence_points = as.integer(nrow(contrast_sequence)),
      contrast_dimension = as.integer(ncol(contrast_sequence)),
      contrast_sequence = contrast_sequence,
      coordinate_map = coordinate_map,
      finite_sequence_only = TRUE,
      asymptotic_flag_extracted = FALSE,
      asymptotic_scales_certified = FALSE,
      subsequence_selected = FALSE,
      p2b_handoff_eligible = FALSE,
      global_boundary_classified = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "This helper evaluates only the exact finite caller-supplied parameter",
        "points. It does not infer an asymptotic flag, choose a subsequence,",
        "certify scales or remainders, use optimizer history, or classify a",
        "global boundary, MML, readiness, uncertainty, or comparison."
      )
    )
  }

  if (!inherits(fit, "mfrm_fit") || is.null(fit$config) ||
      is.null(fit$prep)) {
    return(finish(
      "not_evaluated_fit", FALSE,
      "A current retained mfrm_fit is required.",
      "parameter_sequence_contrast_fit_invalid"
    ))
  }
  config <- fit$config
  if (!identical(config$model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "The nonlinear contrast image helper applies only to GPCM.",
      "parameter_sequence_contrast_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(config$method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional parameter contrast map is not reused for MML.",
      "parameter_sequence_contrast_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (as.integer(config$gpcm_spec$n_params %||% 0L) == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction uses the existing PCM geometry.",
      "no_free_log_slope_coordinate"
    ))
  }
  source_contract <- as.character(
    config$boundary_audit$gpcm_parameter_sequence_flag$contract_version %||%
      ""
  )[1]
  if (!identical(
    source_contract,
    mfrmr_jml_gpcm_parameter_sequence_flag_contract_version()
  )) {
    return(finish(
      "not_evaluated_source_contract", FALSE,
      "The current fit-level P2e contract was unavailable.",
      "parameter_sequence_contrast_source_mismatch"
    ))
  }

  sequence <- tryCatch(
    matrix(
      as.numeric(parameter_sequence),
      nrow = nrow(parameter_sequence), ncol = ncol(parameter_sequence),
      dimnames = dimnames(parameter_sequence)
    ),
    error = function(e) NULL
  )
  controls <- suppressWarnings(as.numeric(c(max_points, max_contrast_elements)))
  valid_controls <- length(controls) == 2L && all(is.finite(controls)) &&
    all(controls == floor(controls)) && all(controls >= 1)
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE,
      "Finite-sequence mapping controls were invalid.",
      "parameter_sequence_contrast_control_invalid"
    ))
  }
  built <- tryCatch({
    idx <- build_indices(
      fit$prep,
      step_facet = config$step_facet,
      slope_facet = config$slope_facet,
      interaction_specs = config$interaction_specs
    )
    sizes <- build_param_sizes(config)
    list(idx = idx, sizes = sizes)
  }, error = function(e) e)
  if (inherits(built, "error")) {
    return(finish(
      "not_evaluated_design", FALSE,
      paste0("The retained response map could not be rebuilt: ",
             conditionMessage(built)),
      "parameter_sequence_contrast_design_invalid"
    ))
  }
  free_dimension <- sum(vapply(built$sizes, as.integer, integer(1)))
  n_observations <- length(built$idx$score_k)
  n_contrasts <- as.integer(config$n_cat) - 1L
  contrast_dimension <- as.double(n_observations) * n_contrasts
  sequence_valid <- !is.null(sequence) && nrow(sequence) >= 1L &&
    ncol(sequence) == free_dimension && all(is.finite(sequence))
  if (!sequence_valid) {
    return(finish(
      "not_evaluated_sequence", FALSE,
      "The explicit parameter sequence was non-finite or dimensionally inconsistent.",
      "parameter_sequence_contrast_sequence_invalid"
    ))
  }
  if (nrow(sequence) > controls[1] ||
      as.double(nrow(sequence)) * contrast_dimension > controls[2]) {
    return(finish(
      "not_evaluated_size_limit", FALSE,
      "The explicit parameter sequence exceeded a mapping workload limit.",
      "parameter_sequence_contrast_size_limit"
    ))
  }

  mapped <- tryCatch(lapply(seq_len(nrow(sequence)), function(point) {
    adjacent <- mfrmr_gpcm_adjacent_logits(
      sequence[point, ], built$idx, config, built$sizes
    )
    adjacent <- matrix(
      adjacent, nrow = n_observations, ncol = n_contrasts
    )
    cumulative <- t(apply(adjacent, 1L, cumsum))
    if (n_observations == 1L) {
      cumulative <- matrix(cumulative, nrow = 1L)
    }
    as.numeric(cumulative)
  }), error = function(e) e)
  if (inherits(mapped, "error")) {
    return(finish(
      "not_evaluated_numeric_map", FALSE,
      paste0("The finite nonlinear contrast map failed: ",
             conditionMessage(mapped)),
      "parameter_sequence_contrast_numeric_map_failed"
    ))
  }
  contrast_sequence <- do.call(rbind, mapped)
  if (any(!is.finite(contrast_sequence))) {
    return(finish(
      "not_evaluated_numeric_map", FALSE,
      "The finite nonlinear contrast map produced a non-finite value.",
      "parameter_sequence_contrast_numeric_map_nonfinite"
    ))
  }
  coordinate_map <- expand.grid(
    Observation = seq_len(n_observations),
    Category = seq_len(n_contrasts),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  coordinate_map$ContrastCoordinate <- seq_len(nrow(coordinate_map))
  coordinate_map <- coordinate_map[
    , c("ContrastCoordinate", "Observation", "Category"), drop = FALSE
  ]
  finish(
    "finite_parameter_sequence_exact_contrast_image", TRUE,
    paste(
      "Every caller-supplied finite free-parameter vector was expanded and",
      "mapped through the exact GPCM slope-times-cumulative-utility response",
      "kernel before reference-category contrasts were retained."
    ),
    "finite_parameter_sequence_exact_nonlinear_contrast_map",
    contrast_sequence = contrast_sequence,
    coordinate_map = coordinate_map
  )
}

audit_mfrm_jml_gpcm_parameter_sequence_flag_scope <- function(
    config,
    sizes,
    n_observations,
    n_categories,
    compactification =
      config$boundary_audit$gpcm_boundary_compactification %||% list(),
    sequence_diagnostic =
      config$boundary_audit$gpcm_optimizer_sequence_diagnostic %||% list(),
    max_oracle_contrast_dimension = 2e6) {
  method <- as.character(config$method %||% "unknown")[1]
  model <- as.character(config$model %||% "unknown")[1]
  size_values <- suppressWarnings(vapply(sizes, as.integer, integer(1)))
  dimension_valid <- length(size_values) > 0L && all(is.finite(size_values)) &&
    all(size_values >= 0L)
  free_dimension <- if (dimension_valid) sum(size_values) else NA_integer_
  observations <- suppressWarnings(as.numeric(n_observations)[1])
  categories <- suppressWarnings(as.numeric(n_categories)[1])
  contrast_dimension <- observations * (categories - 1)
  control <- suppressWarnings(as.numeric(max_oracle_contrast_dimension)[1])
  control_valid <- is.finite(control) && control >= 1 && control == floor(control)
  response_dimension_valid <- is.finite(observations) &&
    observations >= 1 && observations == floor(observations) &&
    is.finite(categories) && categories >= 2 && categories == floor(categories)
  compactification_matches <- identical(
    as.character(compactification$contract_version %||% "")[1],
    mfrmr_jml_gpcm_boundary_compactification_contract_version()
  )
  sequence_diagnostic_matches <- identical(
    as.character(sequence_diagnostic$contract_version %||% "")[1],
    mfrmr_jml_gpcm_sequence_remainder_contract_version()
  )

  finish <- function(state, evaluated, detail, reason_codes,
                     theorem = FALSE, oracle_available = FALSE) {
    list(
      contract_version =
        mfrmr_jml_gpcm_parameter_sequence_flag_contract_version(),
      method = method,
      model = model,
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      state = state,
      evaluated = isTRUE(evaluated),
      free_parameter_dimension = as.integer(free_dimension),
      observations = suppressWarnings(as.integer(observations)),
      categories = suppressWarnings(as.integer(categories)),
      response_contrast_dimension = as.double(contrast_dimension),
      response_coordinate_system =
        "category_logit_contrasts_relative_to_category_zero",
      finite_dimensional_contrast_image = isTRUE(theorem),
      exact_nonlinear_slope_utility_map_precedes_compactification =
        isTRUE(theorem),
      bounded_contrast_sequence_has_convergent_subsequence = isTRUE(theorem),
      unbounded_contrast_sequence_has_scale_separated_flag_subsequence =
        isTRUE(theorem),
      flag_directions_can_be_chosen_orthonormal = isTRUE(theorem),
      successive_flag_scale_ratios_tend_to_zero = isTRUE(theorem),
      terminal_reference_contrast_remainder_tends_to_zero = isTRUE(theorem),
      maximum_divergent_flag_stages = if (isTRUE(theorem)) {
        as.double(contrast_dimension)
      } else {
        NA_real_
      },
      every_parameter_sequence_has_classifiable_further_subsequence_likelihood_limit =
        isTRUE(theorem),
      declared_contrast_flag_oracle_available = isTRUE(oracle_available),
      oracle_workload_admitted = isTRUE(oracle_available),
      source_compactification_contract = as.character(
        compactification$contract_version %||% "missing"
      )[1],
      source_compactification_contract_matches =
        isTRUE(compactification_matches),
      source_finite_sequence_contract = as.character(
        sequence_diagnostic$contract_version %||% "missing"
      )[1],
      source_finite_sequence_contract_matches =
        isTRUE(sequence_diagnostic_matches),
      production_parameter_sequences_declared = 0L,
      production_parameter_sequence_mapped = FALSE,
      production_subsequence_selected = FALSE,
      production_flag_extracted = FALSE,
      optimizer_stage_endpoints_used_as_asymptotic_sequence = FALSE,
      original_sequence_limit_classified = FALSE,
      common_subsequence_limit_certified = FALSE,
      all_subsequence_limits_equal_certified = FALSE,
      competitive_boundary_certified = FALSE,
      global_boundary_classified = FALSE,
      global_finite_maximum_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The theorem classifies a further subsequence, not the original",
        "parameter sequence. Different subsequences may induce different",
        "contrast flags and likelihood limits. The production fit declares",
        "no sequence and extracts no flag. The result does not prove a common",
        "limit, competitiveness, a finite maximizer, global boundary status,",
        "MML geometry, uncertainty, readiness, recovery, simulation, or",
        "external comparison."
      )
    )
  }

  if (!control_valid) {
    return(finish(
      "not_evaluated_control", FALSE,
      "The contrast-oracle workload control was invalid.",
      "parameter_sequence_flag_scope_control_invalid"
    ))
  }
  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "Parameter-sequence contrast flags are attached only to GPCM.",
      "parameter_sequence_flag_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional contrast compactification is not reused for MML.",
      "parameter_sequence_flag_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (as.integer(config$gpcm_spec$n_params %||% 0L) == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction uses the existing PCM geometry.",
      "no_free_log_slope_coordinate"
    ))
  }
  if (!dimension_valid || !response_dimension_valid ||
      !is.finite(contrast_dimension) || contrast_dimension < 1) {
    return(finish(
      "not_evaluated_dimension", FALSE,
      "The retained free-parameter or response-contrast dimension was invalid.",
      "parameter_sequence_flag_dimension_invalid"
    ))
  }
  if (!compactification_matches || !sequence_diagnostic_matches) {
    return(finish(
      "not_evaluated_source_contract", FALSE,
      "Current P1z compactification and P2d finite-sequence contracts are required.",
      "parameter_sequence_flag_source_mismatch"
    ))
  }
  oracle_available <- contrast_dimension <= control
  finish(
    if (oracle_available) {
      "parameter_sequence_further_subsequence_flag_theorem_available_no_sequence_declared"
    } else {
      "parameter_sequence_flag_theorem_available_oracle_workload_not_admitted"
    },
    TRUE,
    paste(
      "The exact nonlinear GPCM map sends every finite parameter point to a",
      "finite-dimensional reference-contrast vector. Recursive compactness",
      "and orthogonal projection give every parameter sequence a further",
      "subsequence with a finite scale-separated contrast flag or a convergent",
      "bounded contrast image, so its conditional likelihood limit is",
      "classifiable along that further subsequence."
    ),
    paste(
      "exact_nonlinear_parameter_to_contrast_map",
      "finite_contrast_flag_subsequence_compactification",
      "further_subsequence_likelihood_limit_classifiable",
      "original_sequence_and_global_boundary_open",
      sep = ";"
    ),
    theorem = TRUE,
    oracle_available = oracle_available
  )
}
