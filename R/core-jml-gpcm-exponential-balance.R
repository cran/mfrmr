# JML GPCM exponential slope-utility balance paths
# ==============================================================================
#
# An affine expanded log-slope path and a finite exponential sum in the free
# additive coordinates give an exact finite exponential sum in every GPCM
# category-logit contrast. Terms with the same combined exponent are added
# before any comparison. Positive combined exponents form a scale-separated
# divergent contrast group, exponent zero forms the finite base, and negative
# exponents form a uniformly vanishing contrast remainder.
#
# This catches a boundary class which a parameter-space direction alone cannot:
# log slopes can diverge while utilities vanish at the inverse exponential rate,
# leaving a bounded response-contrast image. Such a path is an explicit witness
# that the parameter-to-response-contrast map is not proper on that design.

mfrmr_jml_gpcm_exponential_balance_contract_version <- function() {
  "mfrmr-jml-gpcm-exponential-balance-0.2.3-v1"
}


mfrmr_jml_gpcm_exponential_balance_scalar <- function(value) {
  converted <- tryCatch(
    suppressWarnings(as.numeric(value)[1]),
    error = function(e) NA_real_
  )
  if (length(converted) == 0L) NA_real_ else converted
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
      "The response-image escape scope dimension or workload control was invalid.",
      "exponential_balance_scope_control_invalid"
    ))
  }
  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE, FALSE,
      "The exponential-balance response-image scope applies only to GPCM.",
      "exponential_balance_scope_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE, FALSE,
      "The conditional response-image escape scope is not reused for MML.",
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
      "Current finite-limit and global-existence results are required.",
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
      "and handed to the finite-limit contrast analysis. When the retained affine",
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
