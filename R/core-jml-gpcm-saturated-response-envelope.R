# Saturated-response likelihood envelope for zero-offset JML GPCM operators
# ==============================================================================
#
# Strictly positive category mass gives each observation a unique independent
# multinomial optimum.  If that saturated adjacent-logit target is in the
# finite GPCM response image, a finite global JMLE exists.  If it is in the
# response-image closure but not the finite image, its likelihood is the global
# supremum and no finite JMLE exists.  P2j/P2k classify that target and provide
# a finite path.  If the saturated target is outside the closure, the constrained
# optimum over the closure remains a separate nonlinear envelope problem.

mfrmr_jml_gpcm_saturated_envelope_contract_version <- function() {
  "mfrmr-jml-gpcm-saturated-response-envelope-0.2.3-v1"
}

mfrmr_jml_gpcm_saturated_envelope_logsumexp <- function(values) {
  maximum <- max(values)
  maximum + log(sum(exp(values - maximum)))
}

mfrmr_jml_gpcm_saturated_envelope_cumulative <- function(
    adjacent_contrasts,
    observations,
    categories) {
  adjacent <- matrix(
    as.numeric(adjacent_contrasts),
    nrow = observations,
    ncol = categories - 1L
  )
  cumulative <- matrix(
    0,
    nrow = observations,
    ncol = categories
  )
  for (transition in seq_len(categories - 1L)) {
    cumulative[, transition + 1L] <-
      cumulative[, transition] + adjacent[, transition]
  }
  cumulative
}

mfrmr_jml_gpcm_saturated_envelope_loglik <- function(
    adjacent_contrasts,
    category_mass) {
  mass <- as.matrix(category_mass)
  observations <- nrow(mass)
  categories <- ncol(mass)
  utilities <- mfrmr_jml_gpcm_saturated_envelope_cumulative(
    adjacent_contrasts,
    observations,
    categories
  )
  log_denominator <- apply(
    utilities,
    1L,
    mfrmr_jml_gpcm_saturated_envelope_logsumexp
  )
  sum(mass * (utilities - log_denominator))
}

mfrmr_jml_gpcm_saturated_response_envelope <- function(
    adjacent_design,
    slope_owner_by_row,
    category_mass,
    additive_offset = NULL,
    certificate_tolerance = 1e-8,
    max_observations = 20000L,
    max_categories = 20L,
    max_mass_elements = 2e6,
    max_total_mass = 1e12,
    max_coordinates = 2000L,
    max_owners = 12L,
    max_faces = 4095L,
    max_stages = 12L,
    max_hierarchies = 50000L,
    max_equations = 200000L,
    max_lp_nonzeros = 5e6,
    lp_timeout = 30L) {
  design <- mfrmr_jml_gpcm_response_image_sparse_design(adjacent_design)
  owner <- tryCatch(
    suppressWarnings(as.numeric(slope_owner_by_row)),
    error = function(e) NULL
  )
  mass <- tryCatch({
    value <- suppressWarnings(as.matrix(category_mass))
    suppressWarnings(storage.mode(value) <- "double")
    value
  }, error = function(e) NULL)
  offset <- if (is.null(additive_offset) && !is.null(design)) {
    rep(0, nrow(design))
  } else {
    tryCatch(
      suppressWarnings(as.numeric(additive_offset)),
      error = function(e) NULL
    )
  }
  controls <- tryCatch(
    suppressWarnings(as.numeric(c(
      certificate_tolerance, max_observations, max_categories,
      max_mass_elements, max_total_mass, max_coordinates, max_owners,
      max_faces, max_stages, max_hierarchies, max_equations,
      max_lp_nonzeros, lp_timeout
    ))),
    error = function(e) numeric(0)
  )
  observations <- if (is.null(mass)) 0L else nrow(mass)
  categories <- if (is.null(mass)) 0L else ncol(mass)
  rows <- if (is.null(design)) 0L else nrow(design)
  coordinates <- if (is.null(design)) 0L else ncol(design)
  owners <- if (is.null(owner) || length(owner) == 0L ||
                any(!is.finite(owner)) || max(owner) > .Machine$integer.max) {
    0L
  } else {
    as.integer(max(owner))
  }
  owner_integer <- tryCatch(
    suppressWarnings(as.integer(owner)),
    error = function(e) integer(0)
  )
  empty_face <- list()
  empty_higher <- list()
  problem_input_certified <- FALSE

  finish <- function(
      state,
      evaluated,
      classified,
      detail,
      reason_codes,
      independent_adjacent = numeric(0),
      independent_cumulative = matrix(numeric(0), 0L, 0L),
      independent_value = NA_real_,
      face_chart = empty_face,
      higher_order = empty_higher,
      finite_image = FALSE,
      missing_boundary = FALSE,
      outside_closure = FALSE,
      path_source = "none",
      path_id = "") {
    closure_member <- isTRUE(finite_image) || isTRUE(missing_boundary)
    global_classified <- closure_member
    list(
      contract_version =
        mfrmr_jml_gpcm_saturated_envelope_contract_version(),
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      response_likelihood_identity = paste(
        "grouped repeated-row ordered-category log likelihood in complete",
        "adjacent response contrasts"
      ),
      response_contrast_layout = paste(
        "transition-major column order: all observations for transition 1,",
        "then transition 2, and so on"
      ),
      state = state,
      evaluated = isTRUE(evaluated),
      saturated_response_envelope_classified = isTRUE(classified),
      strictly_positive_category_mass = isTRUE(problem_input_certified),
      ordinary_observation_category_separation_absent =
        isTRUE(problem_input_certified),
      independent_saturated_optimum_unique =
        isTRUE(problem_input_certified),
      independent_saturated_optimum_adjacent_contrasts =
        as.numeric(independent_adjacent),
      independent_saturated_optimum_cumulative_contrasts =
        independent_cumulative,
      independent_saturated_optimum_log_likelihood =
        as.numeric(independent_value),
      saturated_log_likelihood_upper_bound = as.numeric(independent_value),
      observations = as.integer(observations),
      categories = as.integer(categories),
      adjacent_response_contrast_dimension = as.integer(rows),
      additive_free_coordinates = as.integer(coordinates),
      slope_owners = as.integer(owners),
      category_mass = mass,
      adjacent_design = design,
      slope_owner_by_row = owner_integer,
      additive_offset = as.numeric(offset %||% numeric(0)),
      positive_mass_response_log_likelihood_strictly_concave =
        isTRUE(problem_input_certified),
      unbounded_response_contrast_likelihood_limit_negative_infinity =
        isTRUE(problem_input_certified),
      every_finite_limsup_parameter_sequence_has_bounded_response_subsequence =
        isTRUE(problem_input_certified),
      face_chart = face_chart,
      higher_order_closure = higher_order,
      saturated_target_finite_response_image_certified =
        isTRUE(finite_image),
      saturated_target_missing_finite_boundary_certified =
        isTRUE(missing_boundary),
      saturated_target_outside_response_image_closure_certified =
        isTRUE(outside_closure),
      saturated_target_closure_membership_classified =
        isTRUE(finite_image) || isTRUE(missing_boundary) ||
        isTRUE(outside_closure),
      complete_instance_parameter_sequence_limsup_envelope_constructed =
        isTRUE(global_classified),
      global_supremum_exists_as_response_closure_maximum_certified =
        isTRUE(global_classified) || isTRUE(outside_closure),
      global_supremum_equals_saturated_upper_bound_certified =
        isTRUE(global_classified),
      global_supremum_strictly_below_saturated_upper_bound_certified =
        isTRUE(outside_closure),
      global_supremum_value_classified = isTRUE(global_classified),
      global_supremum_value = if (global_classified) {
        as.numeric(independent_value)
      } else {
        NA_real_
      },
      response_image_supremum_attained_by_finite_parameter =
        isTRUE(finite_image),
      finite_jmle_existence_certified = isTRUE(finite_image),
      finite_jmle_nonexistence_certified = isTRUE(missing_boundary),
      finite_jmle_nonexistence_without_observation_category_separation_certified =
        isTRUE(missing_boundary),
      saturated_boundary_competitive_certified = isTRUE(missing_boundary),
      global_constrained_optimum_location_and_value_open =
        isTRUE(outside_closure),
      constrained_boundary_likelihood_envelope_constructed = FALSE,
      all_positive_mass_profiles_symbolically_classified = FALSE,
      finite_global_representative = if (isTRUE(finite_image)) {
        face_chart$finite_representative
      } else {
        list()
      },
      supremum_path_source = path_source,
      supremum_path_id = path_id,
      optimizer_trace_used = FALSE,
      optimizer_trace_promoted_to_finite_estimate = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_theorem_only",
      next_gate = if (isTRUE(outside_closure)) {
        "constrained_closure_likelihood_envelope_below_saturated_bound"
      } else if (isTRUE(global_classified)) {
        "apply_saturated_envelope_to_exact_current_fit_operators"
      } else {
        "classified_saturated_target_closure_membership"
      },
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The theorem gives the complete global parameter-sequence envelope",
        "only when the unique independent saturated response optimum is in",
        "the finite response image or its certified closure. If that target",
        "is outside the closure, its value is only an unattainable upper",
        "bound and the constrained closure maximum remains open. The result",
        "does not promote optimizer traces, uncertainty, readiness, MML,",
        "recovery, simulation, or external comparison."
      )
    )
  }

  valid_controls <- length(controls) == 13L && all(is.finite(controls)) &&
    controls[1] > 0 && controls[1] < 1 &&
    all(controls[c(2:4, 6:13)] == floor(controls[c(2:4, 6:13)])) &&
    all(controls[2:13] >= 1)
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "The response-envelope tolerance or workload controls were invalid.",
      "saturated_envelope_control_invalid"
    ))
  }
  valid_design <- !is.null(design) && rows >= 1L && coordinates >= 1L &&
    coordinates <= controls[6] && all(is.finite(design@x))
  valid_mass <- !is.null(mass) && observations >= 1L && categories >= 2L &&
    observations <= controls[2] && categories <= controls[3] &&
    length(mass) <= controls[4] && all(is.finite(mass)) && all(mass > 0) &&
    sum(mass) <= controls[5]
  valid_owner <- !is.null(owner) && length(owner) == rows &&
    all(is.finite(owner)) && all(owner == floor(owner)) && all(owner >= 1) &&
    owners >= 1L && owners <= controls[7] &&
    identical(sort(unique(as.integer(owner))), seq_len(owners))
  valid_offset <- !is.null(offset) && length(offset) == rows &&
    all(is.finite(offset)) && all(offset == 0)
  expected_rows <- observations * (categories - 1L)
  valid_layout <- is.finite(expected_rows) && rows == expected_rows
  if (!valid_design || !valid_mass || !valid_owner || !valid_offset ||
      !valid_layout) {
    return(finish(
      "not_evaluated_input", FALSE, FALSE,
      paste(
        "A finite zero-offset design, contiguous owner map, and strictly",
        "positive observation-by-category mass matrix with matching",
        "transition-major row dimension are required."
      ),
      "saturated_envelope_input_invalid"
    ))
  }
  problem_input_certified <- TRUE

  log_mass <- log(mass)
  independent_adjacent_matrix <-
    log_mass[, 2:categories, drop = FALSE] -
    log_mass[, 1:(categories - 1L), drop = FALSE]
  independent_adjacent <- as.numeric(independent_adjacent_matrix)
  independent_cumulative <-
    log_mass - log_mass[, 1L]
  row_total <- rowSums(mass)
  independent_value <- sum(
    mass * (log_mass - log(row_total))
  )
  direct_value <- mfrmr_jml_gpcm_saturated_envelope_loglik(
    independent_adjacent,
    mass
  )
  objective_scale <- max(1, abs(independent_value), abs(direct_value))
  if (!is.finite(independent_value) || !is.finite(direct_value) ||
      abs(independent_value - direct_value) > controls[1] * objective_scale) {
    return(finish(
      "not_evaluated_saturated_objective", TRUE, FALSE,
      "The independently derived saturated objective values did not agree.",
      "saturated_envelope_objective_reconstruction_mismatch",
      independent_adjacent,
      independent_cumulative,
      independent_value
    ))
  }

  face_chart <- mfrmr_jml_gpcm_response_image_face_chart(
    adjacent_design = design,
    slope_owner_by_row = as.integer(owner),
    target_contrasts = independent_adjacent,
    additive_offset = offset,
    certificate_tolerance = controls[1],
    max_rows = rows,
    max_coordinates = controls[6],
    max_owners = controls[7],
    max_faces = controls[8],
    max_lp_nonzeros = controls[12],
    lp_timeout = controls[13]
  )
  if (!isTRUE(face_chart$response_image_face_chart_classified)) {
    return(finish(
      "not_evaluated_face_chart", TRUE, FALSE,
      "The saturated target P2j face chart was not classified.",
      "saturated_envelope_face_chart_unclassified",
      independent_adjacent,
      independent_cumulative,
      independent_value,
      face_chart
    ))
  }

  higher_order <- mfrmr_jml_gpcm_higher_order_face_lifts(
    face_chart_result = face_chart,
    certificate_tolerance = controls[1],
    max_stages = controls[9],
    max_hierarchies = controls[10],
    max_equations = controls[11],
    max_lp_nonzeros = controls[12],
    lp_timeout = controls[13]
  )
  finite_image <- isTRUE(face_chart$finite_response_image_certified)
  closure_member <- isTRUE(
    higher_order$target_in_finite_response_image_closure_certified
  ) || isTRUE(face_chart$target_in_finite_response_image_closure_certified)
  outside_closure <- isTRUE(
    higher_order$target_outside_finite_response_image_closure_certified
  ) || isTRUE(face_chart$target_outside_closure_outer_bound_certified)
  missing_boundary <- !finite_image && closure_member
  closure_decided <- finite_image || missing_boundary || outside_closure
  if (!closure_decided) {
    return(finish(
      "not_evaluated_closure_membership", TRUE, FALSE,
      "The saturated target closure membership remained unclassified.",
      "saturated_envelope_closure_membership_open",
      independent_adjacent,
      independent_cumulative,
      independent_value,
      face_chart,
      higher_order
    ))
  }

  path_source <- "none"
  path_id <- ""
  if (missing_boundary) {
    higher_table <- higher_order$hierarchy_table
    if (is.data.frame(higher_table) && nrow(higher_table) > 0L &&
        any(higher_table$HierarchyLiftCertified)) {
      winner <- which(higher_table$HierarchyLiftCertified)[1]
      path_source <- "p2k_higher_order_hierarchy"
      path_id <- higher_table$HierarchyId[winner]
    } else {
      face_table <- face_chart$face_table
      winner <- which(face_table$FirstOrderBoundaryReachable)[1]
      if (length(winner) == 1L && is.finite(winner)) {
        path_source <- "p2j_first_order_face"
        path_id <- face_table$FaceId[winner]
      }
    }
    if (identical(path_source, "none")) {
      return(finish(
        "not_evaluated_supremum_path", TRUE, FALSE,
        "Closure membership was certified without an evaluable finite path.",
        "saturated_envelope_path_unavailable",
        independent_adjacent,
        independent_cumulative,
        independent_value,
        face_chart,
        higher_order
      ))
    }
  }

  state <- if (finite_image) {
    "finite_global_jmle_at_saturated_response_optimum_certified"
  } else if (missing_boundary) {
    "global_supremum_nonattained_at_saturated_boundary_certified"
  } else {
    "saturated_upper_bound_outside_closure_constrained_envelope_open"
  }
  detail <- if (finite_image) {
    paste(
      "The unique independent saturated response optimum has a finite",
      "product-one GPCM representative and is the global JML maximum."
    )
  } else if (missing_boundary) {
    paste(
      "The unique independent saturated response optimum is approached by",
      "a certified finite GPCM path but has no finite representative."
    )
  } else {
    paste(
      "The independent saturated response optimum is outside the response",
      "closure; its value is an upper bound, not the constrained supremum."
    )
  }
  finish(
    state, TRUE, TRUE, detail,
    paste0("saturated_envelope_", state),
    independent_adjacent,
    independent_cumulative,
    independent_value,
    face_chart,
    higher_order,
    finite_image,
    missing_boundary,
    outside_closure,
    path_source,
    path_id
  )
}

mfrmr_jml_gpcm_saturated_envelope_path_at <- function(
    result,
    index,
    max_index = 100L,
    max_log_rate_span = 600) {
  valid_source <- is.list(result) && identical(
    as.character(result$contract_version %||% "")[1],
    mfrmr_jml_gpcm_saturated_envelope_contract_version()
  ) && isTRUE(result$saturated_response_envelope_classified) &&
    isTRUE(result$saturated_target_missing_finite_boundary_certified)
  index_value <- tryCatch(
    suppressWarnings(as.numeric(index)[1]),
    error = function(e) NA_real_
  )
  max_value <- tryCatch(
    suppressWarnings(as.numeric(max_index)[1]),
    error = function(e) NA_real_
  )
  valid_index <- is.finite(index_value) && index_value == floor(index_value) &&
    index_value >= 0 && is.finite(max_value) && max_value == floor(max_value) &&
    max_value >= 0 && index_value <= max_value
  if (!valid_source || !valid_index) {
    return(list(
      contract_version =
        mfrmr_jml_gpcm_saturated_envelope_contract_version(),
      valid = FALSE,
      state = "not_evaluated_saturated_envelope_path",
      index = as.integer(if (is.finite(index_value)) index_value else NA),
      reason_codes = "saturated_envelope_path_invalid"
    ))
  }

  source_type <- result$supremum_path_source
  path <- if (identical(source_type, "p2k_higher_order_hierarchy")) {
    mfrmr_jml_gpcm_higher_order_face_path_at(
      result$higher_order_closure,
      result$supremum_path_id,
      index_value,
      max_index = max_value,
      max_log_rate_span = max_log_rate_span
    )
  } else if (identical(source_type, "p2j_first_order_face")) {
    mfrmr_jml_gpcm_response_image_face_path_at(
      result$face_chart,
      result$supremum_path_id,
      index_value,
      max_index = max_value
    )
  } else {
    list(valid = FALSE)
  }
  if (!isTRUE(path$valid)) {
    return(list(
      contract_version =
        mfrmr_jml_gpcm_saturated_envelope_contract_version(),
      valid = FALSE,
      state = "not_evaluated_saturated_envelope_path",
      index = as.integer(index_value),
      reason_codes = "saturated_envelope_source_path_invalid",
      source_path = path
    ))
  }

  log_likelihood <- mfrmr_jml_gpcm_saturated_envelope_loglik(
    path$contrasts,
    result$category_mass
  )
  gap <- result$global_supremum_value - log_likelihood
  list(
    contract_version =
      mfrmr_jml_gpcm_saturated_envelope_contract_version(),
    valid = TRUE,
    state = "finite_saturated_supremum_path_evaluated",
    index = as.integer(index_value),
    source_path_type = source_type,
    source_path_id = result$supremum_path_id,
    target_contrasts =
      result$independent_saturated_optimum_adjacent_contrasts,
    contrasts = path$contrasts,
    contrast_distance_to_target = path$contrast_distance_to_target,
    log_likelihood = log_likelihood,
    target_log_likelihood = result$global_supremum_value,
    target_minus_path_log_likelihood = gap,
    finite_parameter_point = isTRUE(path$finite_parameter_point),
    parameter_path_declared_divergent = isTRUE(
      path$parameter_path_declared_divergent
    ),
    response_contrast_limit_certified = isTRUE(
      path$response_contrast_limit_certified
    ),
    source_path = path,
    reason_codes = "saturated_envelope_finite_path_evaluated"
  )
}

mfrmr_jml_gpcm_saturated_envelope_p2i_binary_fixture <- function() {
  mfrmr_jml_gpcm_saturated_response_envelope(
    adjacent_design = rbind(
      c(1, 0, -1),
      c(0, 1, -1),
      c(1, 0, 1),
      c(0, 1, 1)
    ),
    slope_owner_by_row = c(1, 1, 2, 2),
    category_mass = cbind(
      c(1, 1, 1, 1),
      c(2, 4, 1, 1)
    )
  )
}
