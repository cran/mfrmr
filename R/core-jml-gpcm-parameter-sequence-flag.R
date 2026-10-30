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
      "Current rate compactification and finite-sequence classifications are required.",
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
