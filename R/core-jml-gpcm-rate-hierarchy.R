# JML GPCM lexicographic log-slope rate hierarchy
# ==============================================================================
#
# P1z reduces an unbounded finite-dimensional parameter sequence to convergent
# normalized subsequences and a primary P/Z/L/D log-slope role vector. A zero
# primary-rate coordinate can still diverge at a slower scale. This file
# supplies a pathwise structural recursion for a declared finite hierarchy of
# sum-zero expanded log-slope coefficient vectors.
#
# At stage h only coordinates which were zero at every earlier stage remain
# active. A nonzero stage resolves at least one active coordinate. The first
# nonzero sum-zero stage resolves at least two, so a complete hierarchy has at
# most J - 1 nonzero stages for J slope levels. Later coefficient vectors stay
# sum-zero globally, but their restriction to the active coordinates need not
# be sum-zero: slower compensation may occur in coordinates already resolved
# at a faster scale.
#
# This is a growth-order decomposition, not a likelihood certificate. Additive
# directions and their own slower scales can still change category limits.

mfrmr_jml_gpcm_rate_hierarchy_contract_version <- function() {
  "mfrmr-jml-gpcm-rate-hierarchy-0.2.3-v1"
}


audit_mfrm_jml_gpcm_rate_hierarchy_scope <- function(
    config,
    sizes,
    compactification =
      config$boundary_audit$gpcm_boundary_compactification %||% list()) {
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
    as.character(compactification$contract_version %||% "")[1],
    mfrmr_jml_gpcm_boundary_compactification_contract_version()
  ) && identical(
    as.character(compactification$objective_identity %||% "")[1],
    "identified_conditional_joint_log_likelihood"
  )

  finish <- function(state, evaluated, detail, reason_codes,
                     theorem = FALSE) {
    list(
      contract_version = mfrmr_jml_gpcm_rate_hierarchy_contract_version(),
      method = method,
      model = model,
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      state = state,
      evaluated = isTRUE(evaluated),
      slope_levels = as.integer(n_levels),
      free_log_slope_coordinates = as.integer(free_log_slopes),
      finite_coordinate_rate_hierarchy_theorem = isTRUE(theorem),
      maximum_nonzero_slope_hierarchy_stages = as.integer(
        max(n_levels - 1L, 0L)
      ),
      later_active_restriction_sum_zero_required = FALSE,
      global_coefficient_sum_zero_required = TRUE,
      pathwise_hierarchy_extraction_available = isTRUE(theorem),
      arbitrary_fit_path_hierarchy_computed = FALSE,
      arbitrary_sequence_likelihood_classified = FALSE,
      common_accumulation_direction_likelihood_limit_certified = FALSE,
      additive_hierarchy_classified = FALSE,
      secondary_hierarchy_structurally_bounded = isTRUE(theorem),
      secondary_hierarchy_likelihood_classified = FALSE,
      global_boundary_classified = FALSE,
      global_finite_maximum_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      source_compactification_contract = as.character(
        compactification$contract_version %||% "missing"
      )[1],
      source_compactification_state = as.character(
        compactification$state %||% "missing"
      )[1],
      source_compactification_contract_matches = isTRUE(source_matches),
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "Finite slope-coordinate depth does not classify an arbitrary fit",
        "sequence from its terminal optimizer trace and does not establish a",
        "likelihood limit. Additive directions, ties, bounded remainders,",
        "different accumulation subsequences, and their interactions remain",
        "open. No finite global maximum, boundary absence, uncertainty,",
        "readiness, MML, or external-comparison conclusion follows."
      )
    )
  }

  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "The rate-hierarchy scope applies only to GPCM.",
      "rate_hierarchy_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional fixed-effects rate hierarchy is not reused for MML.",
      "rate_hierarchy_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (!dimension_valid) {
    return(finish(
      "not_evaluated_dimension", FALSE,
      "The identified free log-slope dimension did not match the slope levels.",
      "rate_hierarchy_dimension_invalid"
    ))
  }
  if (free_log_slopes == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction has no relative log-slope hierarchy.",
      "no_free_log_slope_coordinate"
    ))
  }
  finish(
    "finite_coordinate_rate_hierarchy_available_likelihood_open", TRUE,
    paste(
      "Along successively refined subsequences, unresolved zero-rate slope",
      "coordinates can be resolved at strictly slower scales in at most",
      max(n_levels - 1L, 0L), "nonzero stages. This is a structural depth",
      "bound, not an arbitrary-path likelihood classification."
    ),
    paste(
      "finite_slope_coordinate_hierarchy_depth",
      "pathwise_scale_extraction_required",
      "hierarchy_likelihood_open",
      "global_boundary_not_classified",
      sep = ";"
    ),
    theorem = TRUE
  )
}
