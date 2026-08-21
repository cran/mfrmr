# JML GPCM finite optimizer-sequence and remainder diagnostics
# ==============================================================================
#
# The optimizer supplies only finitely many stage endpoints, not an asymptotic
# parameter path. This file permits an explicitly supplied finite coordinate
# sequence to receive a Euclidean low-rank direction/scale diagnostic, but no
# estimated direction or exponent is eligible for the exact P2b oracle.
#
# A separate analytic extension covers a narrow remainder class after P2c:
# finitely many scaled-logit residual directions multiplied by negative powers
# of distance. Their within-row category contrasts vanish uniformly, so the
# declared response probabilities and log likelihood retain the P2b limit.

mfrmr_jml_gpcm_sequence_remainder_contract_version <- function() {
  "mfrmr-jml-gpcm-sequence-remainder-diagnostic-0.2.3-v1"
}

mfrmr_jml_gpcm_sequence_empty_stage_table <- function() {
  data.frame(
    Stage = integer(0),
    SingularValue = numeric(0),
    DisplacementEnergyFraction = numeric(0),
    FinalScale = numeric(0),
    TailScaleMonotone = logical(0),
    EstimatedPowerExponent = numeric(0),
    PowerFitRSquared = numeric(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_sequence_empty_direction_table <- function() {
  data.frame(
    Stage = integer(0),
    OptimizerIndex = integer(0),
    Coordinate = character(0),
    Loading = numeric(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_finite_sequence_diagnostic <- function(
    parameter_sequence,
    distances = NULL,
    reference = NULL,
    coordinate_map = NULL,
    max_stages = 2L,
    tail_points = 4L,
    reconstruction_tolerance = 0.05,
    max_elements = 2e6) {
  empty_stages <- mfrmr_jml_gpcm_sequence_empty_stage_table()
  empty_directions <- mfrmr_jml_gpcm_sequence_empty_direction_table()
  sequence <- tryCatch(
    matrix(
      as.numeric(parameter_sequence),
      nrow = nrow(parameter_sequence), ncol = ncol(parameter_sequence),
      dimnames = dimnames(parameter_sequence)
    ),
    error = function(e) NULL
  )
  n_points <- if (is.null(sequence)) 0L else nrow(sequence)
  n_coordinates <- if (is.null(sequence)) 0L else ncol(sequence)
  max_stages_value <-
    mfrmr_jml_gpcm_parameter_path_scalar_numeric(max_stages)
  tail_points_value <-
    mfrmr_jml_gpcm_parameter_path_scalar_numeric(tail_points)
  tolerance_value <-
    mfrmr_jml_gpcm_parameter_path_scalar_numeric(reconstruction_tolerance)
  max_elements_value <-
    mfrmr_jml_gpcm_parameter_path_scalar_numeric(max_elements)

  finish <- function(state, evaluated, detail, reason_codes,
                     stage_table = empty_stages,
                     direction_table = empty_directions,
                     scales = matrix(numeric(0), nrow = 0L, ncol = 0L),
                     residual_norm = numeric(0),
                     residual_ratio = numeric(0),
                     low_rank_screen = FALSE) {
    list(
      contract_version =
        mfrmr_jml_gpcm_sequence_remainder_contract_version(),
      state = state,
      evaluated = isTRUE(evaluated),
      diagnostic_basis =
        "finite_euclidean_svd_of_explicit_free_parameter_sequence",
      sequence_points = as.integer(n_points),
      free_coordinates = as.integer(n_coordinates),
      reference_source = if (is.null(reference)) {
        "first_sequence_point"
      } else {
        "caller_supplied_reference"
      },
      distance_source = if (is.null(distances)) {
        "not_supplied"
      } else {
        "caller_supplied"
      },
      maximum_diagnostic_stages = suppressWarnings(as.integer(max_stages_value)),
      stage_table = stage_table,
      direction_table = direction_table,
      scale_trajectories = scales,
      residual_norm = residual_norm,
      residual_ratio = residual_ratio,
      finite_low_rank_reconstruction_within_tolerance = isTRUE(low_rank_screen),
      direction_estimates_are_coordinate_basis_dependent = TRUE,
      stage_endpoints_are_optimizer_iterations = FALSE,
      asymptotic_direction_certified = FALSE,
      scale_exponents_estimated = !is.null(distances) && nrow(stage_table) > 0L,
      scale_exponents_certified = FALSE,
      remainder_vanishing_certified = FALSE,
      exact_coefficient_comparisons_eligible = FALSE,
      p2b_handoff_eligible = FALSE,
      path_inferred_from_optimizer_trace = FALSE,
      path_search_performed = FALSE,
      global_boundary_classified = FALSE,
      readiness_effect = "none_diagnostic_only",
      external_comparison_eligible = FALSE,
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The SVD summarizes a finite coordinate sequence in its supplied",
        "Euclidean basis. It does not establish an asymptotic direction, an",
        "exact coefficient zero or tie, a power law, a vanishing remainder,",
        "optimizer-iteration history, a competitive path, or a global",
        "boundary. Estimated directions and exponents cannot enter P2b."
      )
    )
  }

  controls <- c(
    max_stages_value, tail_points_value,
    tolerance_value, max_elements_value
  )
  valid_controls <- all(is.finite(controls)) &&
    max_stages_value >= 1 && max_stages_value <= 2 &&
    max_stages_value == floor(max_stages_value) &&
    tail_points_value >= 2 && tail_points_value == floor(tail_points_value) &&
    tolerance_value >= 0 && max_elements_value >= 1 &&
    max_elements_value == floor(max_elements_value)
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE,
      "Sequence diagnostic controls were invalid.",
      "optimizer_sequence_control_invalid"
    ))
  }
  if (is.null(sequence) || n_points < 3L || n_coordinates < 1L ||
      any(!is.finite(sequence))) {
    return(finish(
      "not_evaluated_sequence", FALSE,
      "At least three finite free-parameter vectors are required.",
      "optimizer_sequence_invalid_or_too_short"
    ))
  }
  if (as.double(n_points) * as.double(n_coordinates) > max_elements_value) {
    return(finish(
      "not_evaluated_size_limit", FALSE,
      "The explicit coordinate sequence exceeded the diagnostic workload.",
      "optimizer_sequence_element_limit"
    ))
  }
  reference_value <- if (is.null(reference)) {
    as.numeric(sequence[1, ])
  } else {
    tryCatch(
      suppressWarnings(as.numeric(reference)),
      error = function(e) NULL
    )
  }
  if (is.null(reference_value) || length(reference_value) != n_coordinates ||
      any(!is.finite(reference_value))) {
    return(finish(
      "not_evaluated_reference", FALSE,
      "The sequence reference was dimensionally inconsistent or non-finite.",
      "optimizer_sequence_reference_invalid"
    ))
  }
  distance_value <- NULL
  if (!is.null(distances)) {
    distance_value <- tryCatch(
      suppressWarnings(as.numeric(distances)),
      error = function(e) NULL
    )
    if (is.null(distance_value) || length(distance_value) != n_points ||
        any(!is.finite(distance_value)) || any(distance_value <= 0) ||
        any(diff(distance_value) <= 0)) {
      return(finish(
        "not_evaluated_distances", FALSE,
        "Supplied distances must be positive, finite, and strictly increasing.",
        "optimizer_sequence_distances_invalid"
      ))
    }
  }

  optimizer_index <- seq_len(n_coordinates)
  coordinate <- colnames(sequence) %||%
    paste0("optimizer:", optimizer_index)
  if (!is.null(coordinate_map)) {
    map <- tryCatch(
      as.data.frame(coordinate_map, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    mapped_index <- if (!is.null(map) && "OptimizerIndex" %in% names(map)) {
      suppressWarnings(as.integer(as.character(map$OptimizerIndex)))
    } else {
      integer(0)
    }
    raw_mapped_index <- if (!is.null(map) &&
                            "OptimizerIndex" %in% names(map)) {
      suppressWarnings(as.numeric(as.character(map$OptimizerIndex)))
    } else {
      numeric(0)
    }
    mapped_coordinate <- if (!is.null(map) && "Coordinate" %in% names(map)) {
      as.character(map$Coordinate)
    } else {
      character(0)
    }
    map_valid <- !is.null(map) && nrow(map) == n_coordinates &&
      all(c("OptimizerIndex", "Coordinate") %in% names(map)) &&
      length(mapped_index) == n_coordinates &&
      length(raw_mapped_index) == n_coordinates &&
      all(is.finite(raw_mapped_index)) &&
      all(raw_mapped_index == mapped_index) &&
      identical(sort(mapped_index), seq_len(n_coordinates)) &&
      length(mapped_coordinate) == n_coordinates &&
      !anyNA(mapped_coordinate) && all(nzchar(mapped_coordinate))
    if (!map_valid) {
      return(finish(
        "not_evaluated_coordinate_map", FALSE,
        "The coordinate map was incomplete or duplicated.",
        "optimizer_sequence_coordinate_map_invalid"
      ))
    }
    optimizer_index <- mapped_index
    coordinate <- mapped_coordinate
  }

  displacement <- sweep(sequence, 2L, reference_value, FUN = "-")
  displacement_energy <- sum(displacement^2)
  if (!is.finite(displacement_energy) || displacement_energy == 0) {
    return(finish(
      "not_evaluated_no_displacement", FALSE,
      "The explicit sequence contains no displacement from its reference.",
      "optimizer_sequence_no_displacement"
    ))
  }
  decomposition <- tryCatch(
    base::svd(
      displacement,
      nu = 0L,
      nv = min(as.integer(max_stages_value), n_points, n_coordinates)
    ),
    error = function(e) NULL
  )
  if (is.null(decomposition) || length(decomposition$d) == 0L ||
      is.null(decomposition$v)) {
    return(finish(
      "not_evaluated_decomposition", FALSE,
      "The finite sequence SVD could not be computed.",
      "optimizer_sequence_svd_failed"
    ))
  }
  rank_threshold <- max(decomposition$d) * .Machine$double.eps *
    max(n_points, n_coordinates)
  available <- which(decomposition$d > rank_threshold)
  stages <- utils::head(available, as.integer(max_stages_value))
  if (length(stages) == 0L) {
    return(finish(
      "not_evaluated_no_displacement", FALSE,
      "The explicit sequence had no numerically nonzero displacement stage.",
      "optimizer_sequence_no_displacement"
    ))
  }
  directions <- decomposition$v[, stages, drop = FALSE]
  scales <- displacement %*% directions
  for (stage in seq_len(ncol(directions))) {
    if (scales[n_points, stage] < 0) {
      directions[, stage] <- -directions[, stage]
      scales[, stage] <- -scales[, stage]
    }
  }
  fitted <- scales %*% t(directions)
  residual <- displacement - fitted
  displacement_norm <- sqrt(rowSums(displacement^2))
  residual_norm <- sqrt(rowSums(residual^2))
  residual_ratio <- ifelse(
    displacement_norm > 0,
    residual_norm / displacement_norm,
    0
  )
  tail_n <- min(as.integer(tail_points_value), n_points)
  tail_index <- seq.int(n_points - tail_n + 1L, n_points)

  stage_rows <- lapply(seq_len(ncol(directions)), function(stage) {
    scale <- as.numeric(scales[, stage])
    exponent <- NA_real_
    r_squared <- NA_real_
    if (!is.null(distance_value)) {
      usable <- is.finite(scale) & abs(scale) > 0
      if (sum(usable) >= 3L) {
        fit <- suppressWarnings(stats::lm(
          log(abs(scale[usable])) ~ log(distance_value[usable])
        ))
        exponent <- as.numeric(stats::coef(fit)[2])
        r_squared <- as.numeric(suppressWarnings(summary(fit))$r.squared)
      }
    }
    data.frame(
      Stage = as.integer(stage),
      SingularValue = as.numeric(decomposition$d[stages[stage]]),
      DisplacementEnergyFraction =
        as.numeric(decomposition$d[stages[stage]]^2 / displacement_energy),
      FinalScale = as.numeric(scale[n_points]),
      TailScaleMonotone = all(diff(scale[tail_index]) >= 0),
      EstimatedPowerExponent = exponent,
      PowerFitRSquared = r_squared,
      stringsAsFactors = FALSE
    )
  })
  direction_rows <- lapply(seq_len(ncol(directions)), function(stage) {
    data.frame(
      Stage = rep(as.integer(stage), n_coordinates),
      OptimizerIndex = as.integer(optimizer_index),
      Coordinate = coordinate,
      Loading = as.numeric(directions[, stage]),
      stringsAsFactors = FALSE
    )
  })
  stage_table <- do.call(rbind, stage_rows)
  direction_table <- do.call(rbind, direction_rows)
  low_rank <- all(
    residual_ratio[tail_index] <= tolerance_value
  )
  finish(
    "finite_sequence_direction_scale_diagnostic", TRUE,
    paste(
      "A finite Euclidean low-rank decomposition was computed for the",
      "explicit coordinate sequence; it has no asymptotic evidential force."
    ),
    "finite_optimizer_sequence_diagnostic_only",
    stage_table = stage_table,
    direction_table = direction_table,
    scales = scales,
    residual_norm = residual_norm,
    residual_ratio = residual_ratio,
    low_rank_screen = low_rank
  )
}

mfrmr_jml_gpcm_remainder_stage_matrices <- function(
    directions, n_observations, n_categories) {
  mfrmr_jml_gpcm_lexicographic_stage_matrices(
    directions, n_observations, n_categories
  )
}

mfrmr_jml_gpcm_decaying_logit_remainder <- function(
    parameter_path_result,
    residual_stage_directions = list(),
    residual_decay_exponents = numeric(0),
    max_stages = 2L) {
  path_valid <- is.list(parameter_path_result) && identical(
    as.character(parameter_path_result$contract_version %||% "")[1],
    mfrmr_jml_gpcm_parameter_path_contract_version()
  ) && isTRUE(
    parameter_path_result$parameter_reachable_limit_contract_completed
  ) && isTRUE(parameter_path_result$limit$declared_path_limit_classified)
  limit <- parameter_path_result$limit %||% list()
  path <- limit$path_definition %||% list()
  n_observations <- nrow(path$base_utilities %||% matrix(numeric(0), 0L, 0L))
  n_categories <- ncol(path$base_utilities %||% matrix(numeric(0), 0L, 0L))
  directions <- mfrmr_jml_gpcm_remainder_stage_matrices(
    residual_stage_directions, n_observations, n_categories
  )
  exponents <- tryCatch(
    suppressWarnings(as.numeric(residual_decay_exponents)),
    error = function(e) NULL
  )
  max_stages_value <-
    mfrmr_jml_gpcm_parameter_path_scalar_numeric(max_stages)
  finish <- function(state, evaluated, certified, detail, reason_codes,
                     stage_table = data.frame(), residual_directions = list()) {
    list(
      contract_version =
        mfrmr_jml_gpcm_sequence_remainder_contract_version(),
      state = state,
      evaluated = isTRUE(evaluated),
      source_parameter_path_contract = as.character(
        parameter_path_result$contract_version %||% "missing"
      )[1],
      source_limit_contract = as.character(limit$contract_version %||% "missing")[1],
      scaled_logit_remainder_class =
        "finite_sum_of_fixed_directions_times_negative_power_scales",
      remainder_stages = as.integer(length(residual_directions)),
      residual_decay_exponents = exponents,
      stage_table = stage_table,
      residual_stage_directions = residual_directions,
      within_row_logit_contrast_remainder_vanishes = isTRUE(certified),
      same_row_probability_limits_certified = isTRUE(certified),
      same_joint_log_likelihood_limit_certified = isTRUE(certified),
      log_likelihood_limit = as.numeric(limit$log_likelihood_limit %||% NA_real_),
      exact_ties_preserved_in_limit = isTRUE(certified),
      finite_distance_ties_required = FALSE,
      arbitrary_utility_remainder_classified = FALSE,
      arbitrary_slope_remainder_classified = FALSE,
      arbitrary_logit_remainder_classified = FALSE,
      optimizer_sequence_remainder_certified = FALSE,
      p2b_limit_extended_only_within_declared_remainder_class = isTRUE(certified),
      global_boundary_classified = FALSE,
      readiness_effect = "none_diagnostic_only",
      external_comparison_eligible = FALSE,
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The theorem concerns additive perturbations on the already scaled",
        "category logits whose within-row contrasts decay as fixed negative",
        "powers. It does not separately classify utility or slope remainders,",
        "infer a remainder from optimizer points, cover bounded oscillations",
        "or slower divergent terms, prove competitiveness, or classify a",
        "global boundary, MML, uncertainty, readiness, or comparison."
      )
    )
  }
  if (!path_valid) {
    return(finish(
      "not_evaluated_source_path", FALSE, FALSE,
      "A completed P2c parameter-reachable limit is required.",
      "decaying_logit_remainder_source_invalid"
    ))
  }
  controls_valid <- is.finite(max_stages_value) && max_stages_value >= 0 &&
    max_stages_value <= 2 &&
    max_stages_value == floor(max_stages_value)
  if (!controls_valid) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "The remainder stage limit must be a finite nonnegative integer.",
      "decaying_logit_remainder_control_invalid"
    ))
  }
  if (is.null(directions) || is.null(exponents) ||
      length(directions) != length(exponents) ||
      length(directions) > max_stages_value ||
      any(!is.finite(exponents)) || any(exponents <= 0) ||
      (length(exponents) > 1L && any(diff(exponents) <= 0))) {
    return(finish(
      "not_evaluated_remainder_mapping", FALSE, FALSE,
      "Remainder directions and strictly increasing decay exponents were invalid.",
      "decaying_logit_remainder_mapping_invalid"
    ))
  }
  if (length(directions) == 0L) {
    return(finish(
      "no_logit_remainder_limit_unchanged", TRUE, TRUE,
      "No scaled-logit remainder was declared.",
      "decaying_logit_remainder_vacuous"
    ))
  }
  spans <- vapply(directions, function(direction) {
    max(apply(direction, 1L, function(row) max(row) - min(row)))
  }, numeric(1))
  stage_table <- data.frame(
    Stage = seq_along(directions),
    DecayExponent = exponents,
    MaximumWithinRowCoefficientSpan = spans,
    PairwiseContrastLimit = rep(0, length(directions)),
    stringsAsFactors = FALSE
  )
  finish(
    "decaying_scaled_logit_remainder_same_limit", TRUE, TRUE,
    paste(
      "Every within-row scaled-logit residual contrast is bounded by a finite",
      "sum of fixed coefficients times negative powers and converges to zero."
    ),
    "decaying_logit_remainder_same_probability_and_likelihood_limit",
    stage_table = stage_table,
    residual_directions = directions
  )
}

mfrmr_jml_gpcm_parameter_remainder_loglik_at <- function(
    parameter_path_result, remainder_result, distance) {
  valid <- is.list(parameter_path_result) && is.list(remainder_result) &&
    identical(
      as.character(parameter_path_result$contract_version %||% "")[1],
      mfrmr_jml_gpcm_parameter_path_contract_version()
    ) && identical(
      as.character(remainder_result$contract_version %||% "")[1],
      mfrmr_jml_gpcm_sequence_remainder_contract_version()
    ) && isTRUE(parameter_path_result$parameter_reachable_limit_contract_completed) &&
    isTRUE(remainder_result$same_joint_log_likelihood_limit_certified)
  distance_value <-
    mfrmr_jml_gpcm_parameter_path_scalar_numeric(distance)
  if (!valid || !is.finite(distance_value) || distance_value <= 0) {
    return(list(
      valid = FALSE,
      distance = distance_value,
      log_likelihood = NA_real_,
      row_log_probability = numeric(0),
      reason_codes = "parameter_remainder_evaluation_invalid"
    ))
  }
  path <- parameter_path_result$limit$path_definition
  log_slopes <- path$base_log_slopes
  if (nrow(path$slope_stage_rates) > 0L) {
    for (stage in seq_len(nrow(path$slope_stage_rates))) {
      log_slopes <- log_slopes +
        distance_value^path$slope_scale_exponents[stage] *
        path$slope_stage_rates[stage, ]
    }
  }
  utilities <- path$base_utilities
  if (length(path$additive_stage_directions) > 0L) {
    for (stage in seq_along(path$additive_stage_directions)) {
      utilities <- utilities +
        distance_value^path$additive_scale_exponents[stage] *
        path$additive_stage_directions[[stage]]
    }
  }
  logits <- utilities * matrix(
    exp(log_slopes[path$slope_index]),
    nrow = nrow(utilities), ncol = ncol(utilities)
  )
  if (length(remainder_result$residual_stage_directions) > 0L) {
    for (stage in seq_along(remainder_result$residual_stage_directions)) {
      logits <- logits +
        distance_value^(-remainder_result$residual_decay_exponents[stage]) *
        remainder_result$residual_stage_directions[[stage]]
    }
  }
  row_log_probability <- rep(NA_real_, nrow(logits))
  for (observation in seq_len(nrow(logits))) {
    row <- logits[observation, ]
    if (any(!is.finite(row))) next
    maximum <- max(row)
    observed <- path$score_k[observation] + 1L
    row_log_probability[observation] <-
      row[observed] - (maximum + log(sum(exp(row - maximum))))
  }
  effective <- path$weight > 0
  valid_numeric <- all(is.finite(row_log_probability[effective]))
  list(
    valid = isTRUE(valid_numeric),
    distance = distance_value,
    log_likelihood = if (valid_numeric) {
      sum(path$weight[effective] * row_log_probability[effective])
    } else {
      NA_real_
    },
    row_log_probability = row_log_probability,
    reason_codes = if (valid_numeric) {
      "parameter_path_with_decaying_logit_remainder_evaluated"
    } else {
      "parameter_path_with_remainder_nonfinite"
    }
  )
}

audit_mfrm_jml_gpcm_optimizer_sequence_scope <- function(
    config,
    sizes,
    stages,
    parameter_sequence,
    parameter_path =
      config$boundary_audit$gpcm_parameter_path_reachability %||% list()) {
  method <- as.character(config$method %||% "unknown")[1]
  model <- as.character(config$model %||% "unknown")[1]
  stage_table <- tryCatch(
    as.data.frame(stages, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  sequence <- tryCatch(
    matrix(
      as.numeric(parameter_sequence),
      nrow = nrow(parameter_sequence), ncol = ncol(parameter_sequence)
    ),
    error = function(e) NULL
  )
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  source_matches <- identical(
    as.character(parameter_path$contract_version %||% "")[1],
    mfrmr_jml_gpcm_parameter_path_contract_version()
  )
  finish <- function(state, evaluated, detail, reason_codes,
                     diagnostic = NULL) {
    list(
      contract_version =
        mfrmr_jml_gpcm_sequence_remainder_contract_version(),
      method = method,
      model = model,
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      state = state,
      evaluated = isTRUE(evaluated),
      stage_endpoint_rows = as.integer(nrow(stage_table)),
      free_coordinates = as.integer(free_dimension),
      aggregate_stage_history_available = nrow(stage_table) > 0L,
      stage_endpoint_coordinates_available_during_fit = !is.null(sequence),
      stage_endpoint_coordinates_retained_in_fit = FALSE,
      stage_endpoints_are_optimizer_iterations = FALSE,
      finite_sequence_diagnostic = diagnostic,
      finite_direction_scale_diagnostic_completed = isTRUE(diagnostic$evaluated),
      source_parameter_path_contract = as.character(
        parameter_path$contract_version %||% "missing"
      )[1],
      source_parameter_path_contract_matches = isTRUE(source_matches),
      production_path_limit_classified = FALSE,
      asymptotic_direction_certified = FALSE,
      scale_exponents_certified = FALSE,
      remainder_vanishing_certified = FALSE,
      p2b_handoff_eligible = FALSE,
      path_inferred_from_optimizer_trace = FALSE,
      path_search_performed = FALSE,
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
        "Only optimizer stage endpoints, not within-stage iterations, are",
        "available transiently. Fewer than three endpoints cannot receive a",
        "direction diagnostic. Any computed SVD is finite, Euclidean-basis",
        "dependent, discarded as an exact path input, and cannot establish",
        "scales, remainders, competitiveness, global status, readiness,",
        "uncertainty, recovery, or external comparison."
      )
    )
  }
  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "The optimizer sequence diagnostic applies only to GPCM.",
      "optimizer_sequence_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional optimizer sequence diagnostic is not reused for MML.",
      "optimizer_sequence_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (as.integer(config$gpcm_spec$n_params %||% 0L) == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction has no free log-slope sequence.",
      "no_free_log_slope_coordinate"
    ))
  }
  if (!source_matches) {
    return(finish(
      "not_evaluated_source_contract", FALSE,
      "The current P2c parameter-path contract was unavailable.",
      "optimizer_sequence_parameter_path_source_mismatch"
    ))
  }
  mapping_valid <- !is.null(sequence) && nrow(sequence) == nrow(stage_table) &&
    ncol(sequence) == free_dimension && all(is.finite(sequence))
  if (!mapping_valid) {
    return(finish(
      "not_evaluated_coordinate_history", FALSE,
      "Transient optimizer stage endpoints were unavailable or misaligned.",
      "optimizer_sequence_coordinate_history_invalid"
    ))
  }
  if (nrow(sequence) < 3L) {
    return(finish(
      "not_evaluated_insufficient_stage_endpoints", FALSE,
      "At least three optimizer stage endpoints are required for a finite sequence diagnostic.",
      "optimizer_sequence_too_few_stage_endpoints"
    ))
  }
  diagnostic <- mfrmr_jml_gpcm_finite_sequence_diagnostic(sequence)
  finish(
    if (isTRUE(diagnostic$evaluated)) {
      "finite_stage_endpoint_direction_scale_diagnostic_only"
    } else {
      "not_evaluated_stage_endpoint_diagnostic"
    },
    isTRUE(diagnostic$evaluated),
    paste(
      "Transient optimizer stage endpoints received a finite low-rank",
      "diagnostic and were not retained as an asymptotic path."
    ),
    diagnostic$reason_codes,
    diagnostic = diagnostic
  )
}
