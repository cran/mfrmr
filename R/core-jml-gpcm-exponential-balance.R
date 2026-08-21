# JML GPCM exponential slope-utility balance paths
# ==============================================================================
#
# An affine expanded log-slope path and a finite exponential sum in the free
# additive coordinates give an exact finite exponential sum in every GPCM
# category-logit contrast. Terms with the same combined exponent are added
# before any comparison. Positive combined exponents form a scale-separated
# P2e flag, exponent zero forms the finite base, and negative exponents are a
# uniformly vanishing contrast remainder.
#
# This catches a boundary class which a parameter-space direction alone cannot:
# log slopes can diverge while utilities vanish at the inverse exponential rate,
# leaving a bounded response-contrast image. Such a path is an explicit witness
# that the parameter-to-response-contrast map is not proper on that design.

mfrmr_jml_gpcm_exponential_balance_contract_version <- function() {
  "mfrmr-jml-gpcm-exponential-balance-0.2.3-v1"
}

mfrmr_jml_gpcm_exponential_balance_empty_rates <- function() {
  data.frame(
    CombinedExponent = numeric(0),
    AsymptoticRole = character(0),
    MaximumAbsoluteCoefficient = numeric(0),
    NonzeroContrastCoordinates = integer(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_exponential_balance_empty_components <- function() {
  data.frame(
    Component = integer(0),
    ComponentType = character(0),
    AdditiveExponent = numeric(0),
    CoordinateL1 = numeric(0),
    CoordinateMaxAbs = numeric(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_exponential_balance_scalar <- function(value) {
  converted <- tryCatch(
    suppressWarnings(as.numeric(value)[1]),
    error = function(e) NA_real_
  )
  if (length(converted) == 0L) NA_real_ else converted
}

mfrmr_jml_gpcm_exponential_balance_cumulative <- function(
    adjacent,
    n_observations,
    n_categories) {
  adjacent_matrix <- matrix(
    as.numeric(adjacent),
    nrow = n_observations,
    ncol = n_categories - 1L
  )
  category <- matrix(0, nrow = n_observations, ncol = n_categories)
  for (transition in seq_len(n_categories - 1L)) {
    category[, transition + 1L] <-
      category[, transition] + adjacent_matrix[, transition]
  }
  category
}

mfrmr_jml_gpcm_exponential_balance_limit <- function(
    adjacent_design,
    additive_center,
    additive_offset = NULL,
    additive_coordinate_directions = list(),
    additive_exponential_rates = numeric(0),
    score_k,
    slope_index,
    base_log_slopes,
    slope_log_rates,
    weight = NULL,
    parameter_map = NULL,
    base_log_slopes_sum_zero_certified = FALSE,
    max_components = 100L,
    max_observations = 20000L,
    max_categories = 20L,
    max_design_nonzeros = 5e6,
    max_exponent_groups = 2000L,
    max_contrast_elements = 2e6) {
  empty_rates <- mfrmr_jml_gpcm_exponential_balance_empty_rates()
  empty_components <-
    mfrmr_jml_gpcm_exponential_balance_empty_components()
  design <- tryCatch(
    methods::as(
      methods::as(
        Matrix::Matrix(adjacent_design, sparse = TRUE), "generalMatrix"
      ),
      "CsparseMatrix"
    ),
    error = function(e) NULL
  )
  center <- tryCatch(
    suppressWarnings(as.numeric(additive_center)),
    error = function(e) NULL
  )
  offset <- if (is.null(additive_offset) && !is.null(design)) {
    rep(0, nrow(design))
  } else {
    tryCatch(
      suppressWarnings(as.numeric(additive_offset)),
      error = function(e) NULL
    )
  }
  directions <- mfrmr_jml_gpcm_parameter_path_directions(
    additive_coordinate_directions
  )
  additive_rates <- tryCatch(
    suppressWarnings(as.numeric(additive_exponential_rates)),
    error = function(e) NULL
  )
  base_logs <- tryCatch(
    suppressWarnings(as.numeric(base_log_slopes)),
    error = function(e) NULL
  )
  slope_rates <- tryCatch(
    suppressWarnings(as.numeric(slope_log_rates)),
    error = function(e) NULL
  )
  scores <- tryCatch(
    suppressWarnings(as.numeric(score_k)),
    error = function(e) NULL
  )
  slope_map <- tryCatch(
    suppressWarnings(as.numeric(slope_index)),
    error = function(e) NULL
  )
  controls <- suppressWarnings(as.numeric(c(
    max_components, max_observations, max_categories,
    max_design_nonzeros, max_exponent_groups, max_contrast_elements
  )))
  n_observations <- if (is.null(scores)) 0L else length(scores)
  n_categories <- if (is.null(design) || n_observations < 1L) {
    0L
  } else {
    as.integer(nrow(design) / n_observations) + 1L
  }
  n_levels <- if (is.null(base_logs)) 0L else length(base_logs)
  free_coordinates <- if (is.null(design)) 0L else ncol(design)
  n_directions <- if (is.null(directions)) 0L else length(directions)
  base_sum_zero_certified <- identical(
    base_log_slopes_sum_zero_certified, TRUE
  )

  finish <- function(state, evaluated, classified, detail, reason_codes,
                     reachability = NULL,
                     rate_table = empty_rates,
                     component_table = empty_components,
                     flag_limit = NULL,
                     path_definition = list(),
                     parameter_divergent = FALSE,
                     slope_divergent = FALSE,
                     additive_divergent = FALSE,
                     bounded_image_escape = FALSE) {
    limit <- mfrmr_jml_gpcm_exponential_balance_scalar(
      flag_limit$log_likelihood_limit %||% NA_real_
    )
    positive_stages <- if (is.null(flag_limit)) 0L else
      as.integer(flag_limit$divergent_flag_stages %||% 0L)
    list(
      contract_version =
        mfrmr_jml_gpcm_exponential_balance_contract_version(),
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      path_family = paste(
        "affine expanded log slopes plus a finite exponential sum in",
        "retained constrained free-additive coordinates"
      ),
      state = state,
      evaluated = isTRUE(evaluated),
      exponential_balance_limit_classified = isTRUE(classified),
      parameter_space_reachability_checked = isTRUE(
        reachability$parameter_space_reachability_checked
      ),
      parameter_space_reachability_certified = isTRUE(
        reachability$parameter_space_reachability_certified
      ),
      exact_slope_utility_product_expanded = isTRUE(classified),
      equal_combined_exponents_aggregated_before_classification =
        isTRUE(classified),
      exact_combined_exponent_comparisons = isTRUE(classified),
      base_log_slopes_sum_zero_structurally_certified =
        isTRUE(base_sum_zero_certified),
      base_log_slope_sum_zero_residual = if (
        length(base_logs) > 0L && all(is.finite(base_logs))
      ) sum(base_logs) else NA_real_,
      slope_log_rate_sum_zero_residual = if (
        length(slope_rates) > 0L && all(is.finite(slope_rates))
      ) sum(slope_rates) else NA_real_,
      observations = as.integer(n_observations),
      categories = as.integer(n_categories),
      slope_levels = as.integer(n_levels),
      additive_free_coordinates = as.integer(free_coordinates),
      additive_exponential_components = as.integer(n_directions),
      combined_exponent_groups = as.integer(nrow(rate_table)),
      positive_contrast_exponent_groups = as.integer(positive_stages),
      zero_contrast_exponent_present = any(
        rate_table$CombinedExponent == 0
      ),
      negative_contrast_exponent_groups = sum(
        rate_table$CombinedExponent < 0
      ),
      parameter_path_divergent = isTRUE(parameter_divergent),
      divergent_log_slope_path = isTRUE(slope_divergent),
      divergent_additive_path = isTRUE(additive_divergent),
      path_escapes_every_bounded_parameter_set = isTRUE(parameter_divergent),
      response_contrast_image_bounded = isTRUE(classified) &&
        positive_stages == 0L,
      bounded_response_image_escape_certified = isTRUE(bounded_image_escape),
      parameter_to_response_contrast_map_properness_refuted =
        isTRUE(bounded_image_escape),
      properness_globally_certified = FALSE,
      complete_bounded_image_escape_family_classified = FALSE,
      complete_escaping_sequence_boundary_envelope_constructed = FALSE,
      log_likelihood_limit = as.numeric(limit),
      finite_log_likelihood_limit = isTRUE(classified) && is.finite(limit),
      negative_infinite_log_likelihood_limit = isTRUE(classified) &&
        is.infinite(limit) && limit < 0,
      rate_table = rate_table,
      component_table = component_table,
      additive_reachability = reachability,
      contrast_flag_limit = flag_limit,
      path_definition = path_definition,
      direct_finite_index_evaluation_available = isTRUE(classified),
      path_inferred_from_optimizer_trace = FALSE,
      exponent_rates_inferred = FALSE,
      path_search_performed = FALSE,
      competitive_boundary_certified = FALSE,
      global_supremum_value_classified = FALSE,
      finite_jmle_existence_certified = FALSE,
      finite_jmle_nonexistence_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The oracle classifies only the declared finite exponential-sum",
        "parameter path. A bounded-image witness refutes properness, but",
        "failure to find one does not prove properness. The result neither",
        "searches all paths nor constructs the complete boundary envelope,",
        "and it does not promote finite attainment, non-attainment, MML,",
        "uncertainty, readiness, recovery, simulation, or comparison."
      )
    )
  }

  valid_controls <- length(controls) == 6L && all(is.finite(controls)) &&
    all(controls == floor(controls)) && controls[1] >= 0 &&
    controls[2] >= 1 && controls[3] >= 2 && all(controls[4:6] >= 1)
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "Path and workload controls must be finite integers in their valid ranges.",
      "exponential_balance_control_invalid"
    ))
  }
  valid_design <- !is.null(design) && n_observations >= 1L &&
    n_categories >= 2L &&
    nrow(design) == n_observations * (n_categories - 1L) &&
    all(is.finite(design@x)) && length(design@x) <= controls[4] &&
    n_observations <= controls[2] && n_categories <= controls[3]
  if (!valid_design) {
    return(finish(
      "not_evaluated_design", FALSE, FALSE,
      "The retained adjacent design was malformed or exceeded a workload limit.",
      "exponential_balance_design_invalid"
    ))
  }
  valid_components <- !is.null(center) && length(center) == free_coordinates &&
    all(is.finite(center)) && !is.null(offset) &&
    length(offset) == nrow(design) && all(is.finite(offset)) &&
    !is.null(directions) &&
    n_directions <= controls[1] && !is.null(additive_rates) &&
    length(additive_rates) == n_directions && all(is.finite(additive_rates)) &&
    all(vapply(directions, function(direction) {
      length(direction) == free_coordinates && all(is.finite(direction))
    }, logical(1)))
  if (!valid_components) {
    return(finish(
      "not_evaluated_additive_components", FALSE, FALSE,
      "The additive center, exponential directions, or rates were malformed.",
      "exponential_balance_additive_components_invalid"
    ))
  }
  valid_slopes <- n_levels >= 2L && !is.null(slope_rates) &&
    length(slope_rates) == n_levels && all(is.finite(base_logs)) &&
    all(is.finite(slope_rates)) &&
    (identical(sum(base_logs), 0) || base_sum_zero_certified) &&
    identical(sum(slope_rates), 0)
  valid_observations <- !is.null(scores) && !is.null(slope_map) &&
    length(slope_map) == n_observations && all(is.finite(scores)) &&
    all(scores == floor(scores)) && all(scores >= 0 & scores < n_categories) &&
    all(is.finite(slope_map)) && all(slope_map == floor(slope_map)) &&
    all(slope_map >= 1 & slope_map <= n_levels)
  if (!valid_slopes || !valid_observations) {
    return(finish(
      "not_evaluated_slope_or_observation_map", FALSE, FALSE,
      "Expanded log slopes, sum-zero rates, scores, or slope indices were invalid.",
      "exponential_balance_slope_observation_invalid"
    ))
  }
  if (is.null(weight)) weight <- rep(1, n_observations)
  weights <- tryCatch(
    suppressWarnings(as.numeric(weight)),
    error = function(e) NULL
  )
  if (is.null(weights) || length(weights) != n_observations ||
      any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
    return(finish(
      "not_evaluated_weight", FALSE, FALSE,
      "Likelihood weights were malformed.",
      "exponential_balance_weight_invalid"
    ))
  }
  scores <- as.integer(scores)
  slope_map <- as.integer(slope_map)

  reachability <- mfrmr_jml_gpcm_additive_utility_path(
    adjacent_design = design,
    coordinate_stage_directions = directions,
    n_observations = n_observations,
    n_categories = n_categories,
    parameter_map = parameter_map,
    max_stages = controls[1],
    max_design_nonzeros = controls[4]
  )
  if (!isTRUE(reachability$parameter_space_reachability_certified)) {
    return(finish(
      "not_evaluated_parameter_reachability", isTRUE(reachability$evaluated),
      FALSE,
      "At least one additive exponential direction was not forward-reachable.",
      reachability$reason_codes,
      reachability = reachability
    ))
  }

  center_category <- mfrmr_jml_gpcm_exponential_balance_cumulative(
    offset + as.numeric(design %*% center), n_observations, n_categories
  )
  category_components <- c(
    list(center_category),
    reachability$category_stage_directions
  )
  component_rates <- c(0, additive_rates)
  component_directions <- c(list(center), directions)
  component_table <- data.frame(
    Component = seq_along(category_components),
    ComponentType = c("finite_additive_center", rep(
      "exponential_additive_direction", n_directions
    )),
    AdditiveExponent = component_rates,
    CoordinateL1 = vapply(
      component_directions, function(value) sum(abs(value)), numeric(1)
    ),
    CoordinateMaxAbs = vapply(
      component_directions, function(value) {
        if (length(value) == 0L) 0 else max(abs(value))
      }, numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  exponent_grid <- outer(
    slope_rates[slope_map], component_rates, "+"
  )
  combined_exponents <- sort(unique(as.numeric(exponent_grid)), decreasing = TRUE)
  if (length(combined_exponents) > controls[5]) {
    return(finish(
      "not_evaluated_exponent_group_limit", FALSE, FALSE,
      "The exact combined exponent expansion exceeded its workload limit.",
      "exponential_balance_exponent_group_limit",
      reachability = reachability,
      component_table = component_table
    ))
  }
  slope_multiplier <- exp(base_logs[slope_map])
  coefficient_matrices <- lapply(combined_exponents, function(exponent) {
    coefficient <- matrix(
      0, nrow = n_observations, ncol = n_categories - 1L
    )
    for (component in seq_along(category_components)) {
      rows <- which(exponent_grid[, component] == exponent)
      if (length(rows) == 0L) next
      coefficient[rows, ] <- coefficient[rows, , drop = FALSE] +
        slope_multiplier[rows] *
        category_components[[component]][rows, -1L, drop = FALSE]
    }
    coefficient
  })
  nonzero <- vapply(coefficient_matrices, function(value) {
    max(abs(value)) > 0
  }, logical(1))
  combined_exponents <- combined_exponents[nonzero]
  coefficient_matrices <- coefficient_matrices[nonzero]
  if (length(combined_exponents) > 0L &&
      as.double(length(combined_exponents)) * n_observations *
        (n_categories - 1L) > controls[6]) {
    return(finish(
      "not_evaluated_contrast_element_limit", FALSE, FALSE,
      "The nonzero combined-exponent contrast expansion exceeded its workload limit.",
      "exponential_balance_contrast_element_limit",
      reachability = reachability,
      component_table = component_table
    ))
  }
  role <- ifelse(
    combined_exponents > 0, "divergent_contrast_flag",
    ifelse(combined_exponents == 0, "finite_contrast_base",
           "vanishing_contrast_remainder")
  )
  rate_table <- if (length(combined_exponents) == 0L) {
    empty_rates
  } else {
    data.frame(
      CombinedExponent = combined_exponents,
      AsymptoticRole = role,
      MaximumAbsoluteCoefficient = vapply(
        coefficient_matrices, function(value) max(abs(value)), numeric(1)
      ),
      NonzeroContrastCoordinates = vapply(
        coefficient_matrices, function(value) sum(value != 0), integer(1)
      ),
      stringsAsFactors = FALSE
    )
  }
  positive <- which(combined_exponents > 0)
  zero <- which(combined_exponents == 0)
  negative <- which(combined_exponents < 0)
  base_contrasts <- matrix(
    0, nrow = n_observations, ncol = n_categories - 1L
  )
  if (length(zero) > 0L) {
    for (position in zero) {
      base_contrasts <- base_contrasts + coefficient_matrices[[position]]
    }
  }
  stage_directions <- coefficient_matrices[positive]
  flag_limit <- mfrmr_jml_gpcm_declared_contrast_flag_limit(
    base_contrasts = base_contrasts,
    score_k = scores,
    stage_directions = stage_directions,
    weight = weights,
    max_stages = controls[5],
    max_observations = controls[2],
    max_categories = controls[3],
    max_contrast_elements = controls[6]
  )
  if (!isTRUE(flag_limit$declared_contrast_flag_limit_classified)) {
    return(finish(
      "not_evaluated_contrast_flag", isTRUE(flag_limit$evaluated), FALSE,
      "The exact combined-exponent response path did not satisfy the P2e oracle.",
      flag_limit$reason_codes,
      reachability = reachability,
      rate_table = rate_table,
      component_table = component_table,
      flag_limit = flag_limit
    ))
  }

  positive_additive_rates <- sort(unique(additive_rates[additive_rates > 0]),
                                  decreasing = TRUE)
  additive_divergent <- any(vapply(positive_additive_rates, function(rate) {
    members <- which(additive_rates == rate)
    aggregate <- Reduce(`+`, directions[members])
    max(abs(aggregate)) > 0
  }, logical(1)))
  slope_divergent <- any(slope_rates != 0)
  parameter_divergent <- slope_divergent || additive_divergent
  bounded_image_escape <- parameter_divergent && length(positive) == 0L
  path_definition <- list(
    adjacent_design = design,
    additive_center = center,
    additive_offset = offset,
    additive_coordinate_directions = directions,
    additive_exponential_rates = additive_rates,
    score_k = scores,
    slope_index = slope_map,
    base_log_slopes = base_logs,
    slope_log_rates = slope_rates,
    weight = weights,
    parameter_map = reachability$parameter_map,
    combined_exponents = combined_exponents,
    coefficient_matrices = coefficient_matrices,
    negative_exponent_positions = negative
  )
  finish(
    if (bounded_image_escape) {
      "bounded_response_image_parameter_escape_certified"
    } else if (parameter_divergent) {
      "divergent_parameter_and_response_flag_limit_classified"
    } else {
      "bounded_parameter_exponential_path_limit_classified"
    },
    TRUE, TRUE,
    paste(
      "The exact slope-times-utility product was grouped by combined",
      "exponent. Positive groups were classified as a P2e flag, the zero",
      "group as its finite base, and all negative groups as a uniformly",
      "vanishing finite exponential remainder."
    ),
    paste(
      "exact_exponential_slope_utility_balance",
      if (bounded_image_escape) {
        "parameter_to_response_contrast_map_nonproper_witness"
      } else {
        "declared_exponential_parameter_path_limit_classified"
      },
      sep = ";"
    ),
    reachability = reachability,
    rate_table = rate_table,
    component_table = component_table,
    flag_limit = flag_limit,
    path_definition = path_definition,
    parameter_divergent = parameter_divergent,
    slope_divergent = slope_divergent,
    additive_divergent = additive_divergent,
    bounded_image_escape = bounded_image_escape
  )
}

mfrmr_jml_gpcm_exponential_balance_loglik_at <- function(result, path_index) {
  valid_result <- is.list(result) && identical(
    as.character(result$contract_version %||% "")[1],
    mfrmr_jml_gpcm_exponential_balance_contract_version()
  ) && isTRUE(result$exponential_balance_limit_classified)
  index <- mfrmr_jml_gpcm_exponential_balance_scalar(path_index)
  if (!valid_result || !is.finite(index) || index < 0) {
    return(list(
      valid = FALSE,
      path_index = index,
      log_likelihood = NA_real_,
      row_log_probability = numeric(0),
      reason_codes = "exponential_balance_finite_index_invalid"
    ))
  }
  path <- result$path_definition
  additive <- path$additive_center
  if (length(path$additive_coordinate_directions) > 0L) {
    for (component in seq_along(path$additive_coordinate_directions)) {
      multiplier <- exp(path$additive_exponential_rates[component] * index)
      additive <- additive +
        multiplier * path$additive_coordinate_directions[[component]]
    }
  }
  log_slopes <- path$base_log_slopes + index * path$slope_log_rates
  adjacent <- tryCatch(
    path$additive_offset + as.numeric(path$adjacent_design %*% additive),
    error = function(e) numeric(0)
  )
  n_observations <- length(path$score_k)
  n_categories <- result$categories
  if (length(adjacent) != n_observations * (n_categories - 1L) ||
      any(!is.finite(adjacent)) || any(!is.finite(log_slopes))) {
    return(list(
      valid = FALSE,
      path_index = index,
      log_likelihood = NA_real_,
      row_log_probability = numeric(0),
      reason_codes = "exponential_balance_finite_index_nonfinite_parameter"
    ))
  }
  utilities <- mfrmr_jml_gpcm_exponential_balance_cumulative(
    adjacent, n_observations, n_categories
  )
  slopes <- exp(log_slopes)
  row_log_probability <- vapply(seq_len(n_observations), function(observation) {
    group <- path$slope_index[observation]
    logits <- slopes[group] * utilities[observation, ]
    if (any(!is.finite(logits))) return(NA_real_)
    observed <- path$score_k[observation] + 1L
    logits[observed] - logsumexp(logits)
  }, numeric(1))
  effective <- path$weight > 0
  numeric_valid <- all(is.finite(row_log_probability[effective]))
  list(
    valid = isTRUE(numeric_valid),
    path_index = index,
    log_likelihood = if (numeric_valid) {
      sum(path$weight[effective] * row_log_probability[effective])
    } else {
      NA_real_
    },
    row_log_probability = row_log_probability,
    maximum_absolute_additive_coordinate = if (length(additive) == 0L) {
      0
    } else {
      max(abs(additive))
    },
    maximum_absolute_log_slope = max(abs(log_slopes)),
    reason_codes = if (numeric_valid) {
      "exponential_balance_finite_index_evaluated"
    } else {
      "exponential_balance_finite_index_nonfinite_logit"
    }
  )
}

mfrmr_jml_gpcm_fit_exponential_balance_limit <- function(
    fit,
    additive_coordinate_directions = list(),
    additive_exponential_rates = numeric(0),
    slope_log_rates,
    additive_center = NULL,
    ...) {
  finish <- function(state, evaluated, detail, reason_codes, result = NULL) {
    list(
      contract_version =
        mfrmr_jml_gpcm_exponential_balance_contract_version(),
      state = state,
      evaluated = isTRUE(evaluated),
      current_fit_reconstructed = !is.null(result),
      exponential_balance_limit_classified = isTRUE(
        result$exponential_balance_limit_classified
      ),
      bounded_response_image_escape_certified = isTRUE(
        result$bounded_response_image_escape_certified
      ),
      parameter_to_response_contrast_map_properness_refuted = isTRUE(
        result$parameter_to_response_contrast_map_properness_refuted
      ),
      log_likelihood_limit = mfrmr_jml_gpcm_exponential_balance_scalar(
        result$log_likelihood_limit %||% NA_real_
      ),
      retained_log_likelihood = if (
        inherits(fit, "mfrm_fit") && is.finite(fit$opt$value %||% NA_real_)
      ) {
        -as.numeric(fit$opt$value)
      } else {
        NA_real_
      },
      boundary_improvement_over_retained = if (
        !is.null(result) && is.finite(result$log_likelihood_limit) &&
          inherits(fit, "mfrm_fit") && is.finite(fit$opt$value %||% NA_real_)
      ) {
        result$log_likelihood_limit + as.numeric(fit$opt$value)
      } else {
        NA_real_
      },
      result = result,
      path_declared_by_caller = TRUE,
      path_inferred_from_optimizer_trace = FALSE,
      path_search_performed = FALSE,
      competitive_boundary_certified = FALSE,
      finite_jmle_existence_certified = FALSE,
      finite_jmle_nonexistence_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "This wrapper reconstructs one caller-declared exponential path on",
        "the current retained fit. It does not search paths, infer rates,",
        "certify competitiveness, or construct a complete boundary envelope."
      )
    )
  }

  if (!inherits(fit, "mfrm_fit") || is.null(fit$config) ||
      is.null(fit$prep) || is.null(fit$opt$par)) {
    return(finish(
      "not_evaluated_fit", FALSE,
      "A current retained mfrm_fit is required.",
      "exponential_balance_fit_invalid"
    ))
  }
  config <- fit$config
  if (!identical(config$model, "GPCM")) {
    return(finish(
      "not_applicable_model", FALSE,
      "The exponential balance wrapper applies only to GPCM.",
      "exponential_balance_not_applicable_model"
    ))
  }
  if (!identical(config$method, "JML")) {
    return(finish(
      "not_applicable_estimator", FALSE,
      "The conditional exponential balance wrapper is not reused for MML.",
      "exponential_balance_not_applicable_estimator"
    ))
  }
  if (as.integer(config$gpcm_spec$n_params %||% 0L) == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction continues to use PCM geometry.",
      "no_free_log_slope_coordinate"
    ))
  }
  source_matches <- identical(
    as.character(
      config$boundary_audit$gpcm_exponential_balance$contract_version %||% ""
    )[1],
    mfrmr_jml_gpcm_exponential_balance_contract_version()
  )
  if (!source_matches) {
    return(finish(
      "not_evaluated_source_contract", FALSE,
      "The current fit-level P2g source contract was unavailable.",
      "exponential_balance_source_mismatch"
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
    adjacent <- mfrmr_estimability_adjacent_design(
      fit$prep, idx, config, sizes,
      include_person = TRUE, include_population_beta = FALSE
    )
    zero_kernel <- mfrmr_gpcm_response_kernel_design(
      fit$prep, idx, config, sizes, rep(0, length(fit$opt$par))
    )
    params <- expand_params(fit$opt$par, sizes, config)
    list(
      idx = idx, sizes = sizes, adjacent = adjacent,
      zero_kernel = zero_kernel, params = params
    )
  }, error = function(e) e)
  if (inherits(built, "error")) {
    return(finish(
      "not_evaluated_design", FALSE,
      paste0("The retained exponential path design failed: ",
             conditionMessage(built)),
      "exponential_balance_retained_design_invalid"
    ))
  }
  retained_additive <-
    fit$opt$par[as.integer(built$adjacent$map$OptimizerIndex)]
  offset <- as.numeric(built$zero_kernel$eta_minus_step)
  center <- if (is.null(additive_center)) {
    retained_additive
  } else {
    additive_center
  }
  result <- mfrmr_jml_gpcm_exponential_balance_limit(
    adjacent_design = built$adjacent$design,
    additive_center = center,
    additive_offset = offset,
    additive_coordinate_directions = additive_coordinate_directions,
    additive_exponential_rates = additive_exponential_rates,
    score_k = built$idx$score_k,
    slope_index = built$idx$slope_idx,
    base_log_slopes = built$params$log_slopes,
    slope_log_rates = slope_log_rates,
    weight = built$idx$weight,
    parameter_map = built$adjacent$map,
    base_log_slopes_sum_zero_certified = TRUE,
    ...
  )
  if (!isTRUE(result$exponential_balance_limit_classified)) {
    return(finish(
      "not_evaluated_exponential_balance", isTRUE(result$evaluated),
      "The caller-declared current-fit path was not classified.",
      result$reason_codes,
      result = result
    ))
  }
  finish(
    paste0("current_fit_", result$state), TRUE,
    paste(
      "The current fit's retained additive operator and base log slopes were",
      "reconstructed before the declared exponential balance path was",
      "classified."
    ),
    result$reason_codes,
    result = result
  )
}

audit_mfrm_jml_gpcm_exponential_balance_scope <- function(
    config,
    sizes,
    parameter_sequence_flag =
      config$boundary_audit$gpcm_parameter_sequence_flag %||% list(),
    global_existence =
      config$boundary_audit$gpcm_global_existence %||% list(),
    canonical_zero_utility_center_certified = FALSE,
    max_oracle_free_coordinates = 10000L) {
  method <- as.character(config$method %||% "unknown")[1]
  model <- as.character(config$model %||% "unknown")[1]
  free_dimension <- suppressWarnings(sum(vapply(sizes, as.integer, integer(1))))
  control <- mfrmr_jml_gpcm_exponential_balance_scalar(
    max_oracle_free_coordinates
  )
  p2e_matches <- identical(
    as.character(parameter_sequence_flag$contract_version %||% "")[1],
    mfrmr_jml_gpcm_parameter_sequence_flag_contract_version()
  )
  p2f_matches <- identical(
    as.character(global_existence$contract_version %||% "")[1],
    mfrmr_jml_gpcm_global_existence_contract_version()
  )
  canonical_zero <- identical(
    canonical_zero_utility_center_certified, TRUE
  )
  finish <- function(state, evaluated, available, detail, reason_codes) {
    list(
      contract_version =
        mfrmr_jml_gpcm_exponential_balance_contract_version(),
      method = method,
      model = model,
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      state = state,
      evaluated = isTRUE(evaluated),
      free_parameter_dimension = as.integer(free_dimension),
      declared_exponential_balance_oracle_available = isTRUE(available),
      exact_combined_exponent_aggregation_available = isTRUE(available),
      bounded_response_image_escape_witness_classifiable = isTRUE(available),
      source_parameter_sequence_flag_contract_matches = isTRUE(p2e_matches),
      source_global_existence_contract_matches = isTRUE(p2f_matches),
      canonical_zero_utility_center_certified =
        isTRUE(evaluated) && isTRUE(canonical_zero),
      canonical_affine_additive_offset_exactly_zero =
        isTRUE(evaluated) && isTRUE(canonical_zero),
      canonical_zero_utility_divergent_slope_path_certified =
        isTRUE(evaluated) && isTRUE(canonical_zero),
      canonical_path_additive_coordinates_identically_zero =
        isTRUE(evaluated) && isTRUE(canonical_zero),
      canonical_path_response_contrasts_identically_zero =
        isTRUE(evaluated) && isTRUE(canonical_zero),
      production_exponential_paths_declared = 0L,
      production_canonical_paths_evaluated = if (
        isTRUE(evaluated) && isTRUE(canonical_zero)
      ) 1L else 0L,
      production_bounded_image_escape_search_performed = FALSE,
      production_bounded_image_escape_certified =
        isTRUE(evaluated) && isTRUE(canonical_zero),
      parameter_to_response_contrast_map_properness_classified =
        isTRUE(evaluated) && isTRUE(canonical_zero),
      parameter_to_response_contrast_map_properness_refuted =
        isTRUE(evaluated) && isTRUE(canonical_zero),
      complete_bounded_image_escape_family_classified = FALSE,
      complete_escaping_sequence_boundary_envelope_constructed = FALSE,
      finite_jmle_existence_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "When the affine additive offset is exactly zero, the canonical all-",
        "zero free-additive path proves nonproperness without a search:",
        "sum-zero log slopes can diverge while every response contrast stays",
        "zero. Otherwise this canonical witness remains unclassified. One",
        "finite exponential family is not the complete escaping-sequence",
        "boundary envelope required for finite attainment."
      )
    )
  }
  valid_control <- is.finite(control) && control >= 1 &&
    control == floor(control) && is.finite(free_dimension) &&
    free_dimension >= 0
  if (!valid_control) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "The P2g scope dimension or workload control was invalid.",
      "exponential_balance_scope_control_invalid"
    ))
  }
  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE, FALSE,
      "The P2g exponential balance scope applies only to GPCM.",
      "exponential_balance_scope_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE, FALSE,
      "The conditional P2g scope is not reused for MML.",
      "exponential_balance_scope_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (as.integer(config$gpcm_spec$n_params %||% 0L) == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE, FALSE,
      "The exact unit-slope reduction uses PCM geometry.",
      "no_free_log_slope_coordinate"
    ))
  }
  if (!p2e_matches || !p2f_matches) {
    return(finish(
      "not_evaluated_source_contract", FALSE, FALSE,
      "Current P2e and P2f source contracts are required.",
      "exponential_balance_scope_source_mismatch"
    ))
  }
  available <- free_dimension <= control
  finish(
    if (canonical_zero) {
      "canonical_zero_utility_bounded_image_escape_certified"
    } else if (available) {
      "declared_exponential_balance_oracle_available_no_path_declared"
    } else {
      "exponential_balance_oracle_workload_not_admitted"
    },
    TRUE, available,
    paste(
      "Finite exponential sums in retained additive coordinates and affine",
      "sum-zero log-slope paths can be expanded exactly by combined exponent",
      "and handed to the P2e contrast-flag oracle. When the retained affine",
      "additive offset is exactly zero, the all-zero free-additive path leaves",
      "every response contrast identically zero while a nonzero sum-zero log-",
      "slope direction escapes to infinity, refuting properness."
    ),
    paste(
      "exact_slope_utility_combined_exponent_expansion",
      if (canonical_zero) {
        "canonical_zero_utility_bounded_image_escape_certified"
      } else {
        "canonical_zero_utility_center_not_certified"
      },
      if (canonical_zero) {
        "parameter_to_response_contrast_map_properness_refuted"
      } else {
        "parameter_to_response_contrast_map_properness_open"
      },
      "production_path_search_not_performed",
      sep = ";"
    )
  )
}
