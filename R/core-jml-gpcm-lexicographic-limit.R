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
        "The fit-level object advertises a non-public diagnostic path; it",
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
