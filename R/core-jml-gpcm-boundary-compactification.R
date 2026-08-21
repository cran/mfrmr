# JML GPCM boundary compactification and primary-rate hierarchy
# ==============================================================================
#
# An arbitrary unbounded sequence in the finite-dimensional identified JML
# parameter space need not have a globally convergent normalized direction.
# The unit sphere is compact, however, so it always has a subsequence whose
# normalized additive and expanded sum-zero log-slope displacement converges.
# A nonzero limiting log-slope direction has a finite asymptotic role pattern:
# positive (P), zero (Z), leading negative (L), or deeper negative (D).
#
# This is a scope theorem, not a boundary-likelihood certificate. Zero primary
# slope rates can carry a slower divergent scale, and the all-zero primary
# slope sector can do the same. Repeating the argument may therefore require a
# secondary hierarchy which this audit deliberately does not classify.

mfrmr_jml_gpcm_boundary_compactification_contract_version <- function() {
  "mfrmr-jml-gpcm-boundary-compactification-0.2.3-v1"
}

mfrmr_jml_gpcm_primary_role_counts <- function(n_levels) {
  n_levels <- suppressWarnings(as.integer(n_levels)[1])
  if (!is.finite(n_levels) || n_levels < 2L) {
    return(data.frame(
      SlopeLevels = as.integer(max(n_levels, 0L)),
      CanonicalPrimaryRolePatterns = 0,
      ZeroFreePrimaryRolePatterns = 0,
      ZeroContainingPrimaryRolePatterns = 0,
      OrderedPairRolePatterns = 0,
      AdditionalRolePatterns = 0,
      stringsAsFactors = FALSE
    ))
  }
  total <- 4^n_levels - 2 * 3^n_levels + 2^n_levels
  zero_free <- 3^n_levels - 2 * 2^n_levels + 1
  pairs <- as.double(n_levels) * as.double(n_levels - 1L)
  data.frame(
    SlopeLevels = n_levels,
    CanonicalPrimaryRolePatterns = as.double(total),
    ZeroFreePrimaryRolePatterns = as.double(zero_free),
    ZeroContainingPrimaryRolePatterns = as.double(total - zero_free),
    OrderedPairRolePatterns = as.double(pairs),
    AdditionalRolePatterns = as.double(total - pairs),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_rate_role_partition <- function(
    expanded_rates,
    zero_tolerance = 1e-10) {
  rates <- suppressWarnings(as.numeric(expanded_rates))
  tolerance <- suppressWarnings(as.numeric(zero_tolerance)[1])
  valid_tolerance <- length(tolerance) == 1L && is.finite(tolerance) &&
    tolerance >= 0
  finite_input <- length(rates) >= 2L && all(is.finite(rates))
  input_scale <- if (finite_input) max(1, max(abs(rates))) else NA_real_
  sum_zero_residual <- if (finite_input) sum(rates) else NA_real_
  sum_zero <- valid_tolerance && finite_input &&
    abs(sum_zero_residual) <= tolerance * input_scale * length(rates)
  magnitude <- if (finite_input) max(abs(rates)) else NA_real_
  all_zero <- valid_tolerance && finite_input &&
    is.finite(magnitude) && magnitude <= tolerance
  normalized <- if (finite_input && is.finite(magnitude) && magnitude > 0) {
    rates / magnitude
  } else {
    rep(NA_real_, length(rates))
  }
  roles <- rep("invalid", length(rates))
  canonical <- rep(NA_real_, length(rates))
  state <- "invalid_rate_vector"
  valid_direction <- FALSE

  if (valid_tolerance && finite_input && sum_zero && all_zero) {
    roles[] <- "Z"
    canonical[] <- 0
    state <- "all_zero_primary_slope_rate"
  } else if (valid_tolerance && finite_input && sum_zero &&
             is.finite(magnitude) && magnitude > tolerance) {
    positive <- normalized > tolerance
    negative <- normalized < -tolerance
    zero <- !(positive | negative)
    if (any(positive) && any(negative)) {
      leading_value <- max(normalized[negative])
      leading <- negative &
        abs(normalized - leading_value) <= tolerance
      deeper <- negative & !leading
      roles[positive] <- "P"
      roles[zero] <- "Z"
      roles[leading] <- "L"
      roles[deeper] <- "D"
      canonical[positive] <-
        (sum(leading) + 2 * sum(deeper)) / sum(positive)
      canonical[zero] <- 0
      canonical[leading] <- -1
      canonical[deeper] <- -2
      valid_direction <- TRUE
      state <- "nonzero_sum_zero_primary_rate_partitioned"
    }
  }

  list(
    state = state,
    valid_direction = isTRUE(valid_direction),
    finite_input = isTRUE(finite_input),
    sum_zero = isTRUE(sum_zero),
    sum_zero_residual = as.numeric(sum_zero_residual),
    all_zero_primary_slope_rate = isTRUE(all_zero) && isTRUE(sum_zero),
    expanded_rates = rates,
    normalized_rates = normalized,
    roles = roles,
    role_signature = paste(roles, collapse = ""),
    positive_indices = which(roles == "P"),
    zero_indices = which(roles == "Z"),
    leading_negative_indices = which(roles == "L"),
    deeper_negative_indices = which(roles == "D"),
    canonical_equivalent_rates = canonical,
    canonical_equivalent_sum_zero = all(is.finite(canonical)) &&
      abs(sum(canonical)) <= tolerance * max(1, length(canonical)),
    secondary_hierarchy_required = isTRUE(sum_zero) &&
      (isTRUE(all_zero) || any(roles == "Z")),
    zero_tolerance = tolerance,
    criterion_frozen_for_scientific_claims = FALSE
  )
}

mfrmr_jml_gpcm_compactify_direction <- function(
    additive_displacement,
    expanded_log_slope_displacement,
    zero_tolerance = 1e-10) {
  additive <- suppressWarnings(as.numeric(additive_displacement))
  slopes <- suppressWarnings(as.numeric(expanded_log_slope_displacement))
  tolerance <- suppressWarnings(as.numeric(zero_tolerance)[1])
  finite_input <- length(slopes) >= 2L &&
    all(is.finite(additive)) && all(is.finite(slopes)) &&
    is.finite(tolerance) && tolerance >= 0
  combined <- c(additive, slopes)
  scale <- if (finite_input && length(combined) > 0L) {
    max(abs(combined))
  } else {
    NA_real_
  }
  norm <- if (is.finite(scale) && scale > 0) {
    scale * sqrt(sum((combined / scale)^2))
  } else {
    NA_real_
  }
  valid <- finite_input && is.finite(norm) && norm > 0 &&
    abs(sum(slopes)) <= tolerance * max(1, max(abs(slopes))) *
      length(slopes)
  normalized_additive <- if (valid) additive / norm else rep(NA_real_, length(additive))
  normalized_slopes <- if (valid) slopes / norm else rep(NA_real_, length(slopes))
  rate_partition <- mfrmr_jml_gpcm_rate_role_partition(
    normalized_slopes, zero_tolerance = tolerance
  )
  slope_primary_zero <- valid &&
    max(abs(normalized_slopes)) <= tolerance
  state <- if (!valid) {
    "invalid_direction"
  } else if (slope_primary_zero) {
    "additive_primary_zero_slope_rate"
  } else if (isTRUE(rate_partition$valid_direction)) {
    "nonzero_slope_primary_direction_partitioned"
  } else {
    "indeterminate_primary_direction"
  }
  list(
    state = state,
    valid_direction = isTRUE(valid),
    coordinate_system =
      "free_additive_plus_expanded_sum_zero_log_slope_euclidean",
    combined_norm = as.numeric(norm),
    normalized_additive = normalized_additive,
    normalized_expanded_log_slopes = normalized_slopes,
    normalized_norm = if (valid) {
      sqrt(sum(normalized_additive^2) + sum(normalized_slopes^2))
    } else {
      NA_real_
    },
    slope_primary_zero = isTRUE(slope_primary_zero),
    rate_partition = rate_partition,
    secondary_hierarchy_required = isTRUE(valid) &&
      (isTRUE(slope_primary_zero) ||
         isTRUE(rate_partition$secondary_hierarchy_required)),
    zero_tolerance = tolerance,
    criterion_frozen_for_scientific_claims = FALSE
  )
}

audit_mfrm_jml_gpcm_boundary_compactification <- function(
    config,
    sizes,
    boundary_classification =
      config$boundary_audit$gpcm_fixed_objective_classification %||% list(),
    affine_transport =
      config$boundary_audit$gpcm_asymptotically_affine_transport %||% list(),
    zero_tolerance = 1e-10) {
  method <- as.character(config$method %||% "unknown")[1]
  model <- as.character(config$model %||% "unknown")[1]
  levels <- as.character(config$gpcm_spec$levels %||% character(0))
  n_levels <- length(levels)
  free_log_slopes <- as.integer(config$gpcm_spec$n_params %||% 0L)
  size_values <- suppressWarnings(vapply(sizes, as.integer, integer(1)))
  dimension_valid <- length(size_values) > 0L && all(is.finite(size_values)) &&
    all(size_values >= 0L)
  size_log_slopes <- if ("log_slopes" %in% names(size_values)) {
    as.integer(size_values[["log_slopes"]])
  } else {
    NA_integer_
  }
  total_dimension <- if (dimension_valid) sum(size_values) else NA_integer_
  additive_dimension <- if (dimension_valid) {
    total_dimension - size_log_slopes
  } else {
    NA_integer_
  }
  role_counts <- mfrmr_jml_gpcm_primary_role_counts(n_levels)
  tolerance <- suppressWarnings(as.numeric(zero_tolerance)[1])
  tolerance_valid <- is.finite(tolerance) && tolerance >= 0
  source_contract_matches <- identical(
    as.character(boundary_classification$contract_version %||% "")[1],
    mfrmr_jml_gpcm_boundary_classification_contract_version()
  ) && identical(
    as.character(boundary_classification$objective_identity %||% "")[1],
    "identified_conditional_joint_log_likelihood"
  )
  transport_contract_matches <- identical(
    as.character(affine_transport$contract_version %||% "")[1],
    mfrmr_jml_gpcm_asymptotically_affine_contract_version()
  )

  finish <- function(state, evaluated, detail, reason_codes,
                     compactness = FALSE, role_partition = FALSE,
                     hierarchy_required = FALSE) {
    list(
      contract_version =
        mfrmr_jml_gpcm_boundary_compactification_contract_version(),
      method = method,
      model = model,
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      coordinate_system =
        "free_additive_plus_expanded_sum_zero_log_slope_euclidean",
      state = state,
      evaluated = isTRUE(evaluated),
      finite_parameter_dimension = dimension_valid,
      combined_free_dimension = as.integer(total_dimension),
      additive_free_coordinates = as.integer(additive_dimension),
      slope_levels = as.integer(n_levels),
      free_log_slope_coordinates = as.integer(free_log_slopes),
      normalized_direction_space = "finite_dimensional_unit_sphere",
      unit_sphere_compactness_applied = isTRUE(compactness),
      unbounded_sequence_has_convergent_normalized_subsequence =
        isTRUE(compactness),
      whole_sequence_normalized_direction_convergence_required = FALSE,
      accumulation_direction_reduction_available = isTRUE(compactness),
      canonical_primary_role_partition_complete = isTRUE(role_partition),
      role_counts = role_counts,
      all_zero_slope_primary_sector_present = isTRUE(hierarchy_required),
      zero_primary_rate_secondary_hierarchy_required =
        isTRUE(hierarchy_required),
      secondary_hierarchy_classified = FALSE,
      finite_secondary_hierarchy_depth_certified = FALSE,
      canonical_role_boundary_search_complete = isTRUE(
        boundary_classification$canonical_constant_rate_scope_complete
      ),
      source_boundary_contract = as.character(
        boundary_classification$contract_version %||% "missing"
      )[1],
      source_boundary_state = as.character(
        boundary_classification$state %||% "missing"
      )[1],
      source_boundary_contract_matches = isTRUE(source_contract_matches),
      affine_transport_contract = as.character(
        affine_transport$contract_version %||% "missing"
      )[1],
      affine_transport_state = as.character(
        affine_transport$state %||% "missing"
      )[1],
      affine_transport_contract_matches = isTRUE(transport_contract_matches),
      nonconvergent_normalized_rate_boundary_classified = FALSE,
      arbitrary_nonvanishing_residual_boundary_classified = FALSE,
      general_curved_path_classified = FALSE,
      global_boundary_classified = FALSE,
      global_finite_maximum_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      zero_tolerance = tolerance,
      criterion_frozen_for_scientific_claims = FALSE,
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "Compactness supplies convergent normalized subsequences, not a",
        "common limit for the original sequence or a boundary-likelihood",
        "classification. The finite P/Z/L/D partition covers only the primary",
        "nonzero expanded log-slope direction. Zero primary-rate coordinates",
        "and the all-zero slope-primary sector can contain slower divergent",
        "scales. Their secondary hierarchy, arbitrary nonvanishing residuals,",
        "global curved paths, finite global maxima, boundary absence, MML",
        "geometry, uncertainty, readiness, and external comparison remain",
        "unclassified."
      )
    )
  }

  if (!tolerance_valid) {
    return(finish(
      "not_evaluated_control", FALSE,
      "The numerical zero-role tolerance must be finite and nonnegative.",
      "boundary_compactification_control_invalid"
    ))
  }
  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "Boundary compactification is attached only to GPCM fits.",
      "boundary_compactification_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional fixed-effects compactification is not reused for MML.",
      "boundary_compactification_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (!dimension_valid || !is.finite(size_log_slopes) ||
      !is.finite(free_log_slopes) || size_log_slopes != free_log_slopes ||
      free_log_slopes != max(n_levels - 1L, 0L)) {
    return(finish(
      "not_evaluated_dimension", FALSE,
      "The identified free dimension or expanded slope-level map was invalid.",
      "boundary_compactification_dimension_invalid"
    ))
  }
  if (free_log_slopes == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction has no relative log-slope hierarchy.",
      "no_free_log_slope_coordinate",
      compactness = TRUE
    ))
  }
  finish(
    "primary_subsequence_partitioned_secondary_hierarchy_open", TRUE,
    paste(
      "Every unbounded parameter sequence admits a convergent normalized",
      "subsequence. Every nonzero limiting expanded log-slope direction has",
      "one enumerated P/Z/L/D primary role pattern, while zero primary-rate",
      "coordinates retain an unresolved slower-scale hierarchy."
    ),
    paste(
      "normalized_subsequence_compactness",
      "canonical_primary_rate_roles_enumerated",
      "zero_primary_rate_secondary_hierarchy_open",
      "global_boundary_not_classified",
      sep = ";"
    ),
    compactness = TRUE,
    role_partition = TRUE,
    hierarchy_required = TRUE
  )
}
