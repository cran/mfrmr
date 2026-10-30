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
    mfrmr_jml_gpcm_exponential_balance_scalar(max_stages)
  tail_points_value <-
    mfrmr_jml_gpcm_exponential_balance_scalar(tail_points)
  tolerance_value <-
    mfrmr_jml_gpcm_exponential_balance_scalar(reconstruction_tolerance)
  max_elements_value <-
    mfrmr_jml_gpcm_exponential_balance_scalar(max_elements)

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
        "boundary. Estimated directions and exponents cannot enter the finite lexicographic-limit check."
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
      "The required parameter-path classification was unavailable.",
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
