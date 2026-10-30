# JML GPCM parameter-path reachability for declared lexicographic limits
# ==============================================================================
#
# P2b accepts cumulative category-utility directions directly. An inverse
# column-space test would make reachability tolerance-dependent and could alter
# an exact category tie. This contract therefore takes the safer forward route:
# the caller declares directions in the retained constrained free-additive
# coordinates, the exact retained adjacent-utility operator maps them to
# adjacent directions, and cumulative utilities are built with category zero
# fixed at zero. The resulting utility path is reachable by construction in
# the implementation parameter space.
#
# No direction, scale, or asymptotic path is inferred from an optimizer trace.

mfrmr_jml_gpcm_parameter_path_contract_version <- function() {
  "mfrmr-jml-gpcm-parameter-path-reachability-0.2.3-v1"
}


audit_mfrm_jml_gpcm_parameter_path_scope <- function(
    prep,
    idx,
    config,
    sizes,
    lexicographic_limit =
      config$boundary_audit$gpcm_lexicographic_limit %||% list(),
    max_design_nonzeros = 5e6) {
  method <- as.character(config$method %||% "unknown")[1]
  model <- as.character(config$model %||% "unknown")[1]
  source_matches <- identical(
    as.character(lexicographic_limit$contract_version %||% "")[1],
    mfrmr_jml_gpcm_lexicographic_limit_contract_version()
  )
  empty_map <- data.frame()
  finish <- function(state, evaluated, detail, reason_codes,
                     available = FALSE, design = NULL, map = empty_map) {
    list(
      contract_version =
        mfrmr_jml_gpcm_parameter_path_contract_version(),
      method = method,
      model = model,
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      state = state,
      evaluated = isTRUE(evaluated),
      declared_free_coordinate_path_wrapper_available = isTRUE(available),
      forward_parameter_reachability_operator_available = isTRUE(available),
      additive_design_rows = if (is.null(design)) 0L else nrow(design),
      additive_free_coordinates = if (is.null(design)) 0L else ncol(design),
      additive_design_nonzeros = if (is.null(design)) 0 else length(design@x),
      coordinate_map = map,
      declared_paths_evaluated = 0L,
      production_parameter_path_reachability_checked = FALSE,
      production_path_limit_classified = FALSE,
      source_lexicographic_limit_contract = as.character(
        lexicographic_limit$contract_version %||% "missing"
      )[1],
      source_lexicographic_limit_contract_matches = isTRUE(source_matches),
      inverse_projection_used = FALSE,
      path_inferred_from_optimizer_trace = FALSE,
      scale_exponents_inferred = FALSE,
      path_search_performed = FALSE,
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
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The fit-level object records a forward operator and wrapper scope but",
        "does not declare, infer, or evaluate any path. Only caller-declared",
        "free-coordinate directions can receive a retained-design reachability",
        "certificate. Remainders, path search, competitiveness, common",
        "subsequence limits, global status, MML, uncertainty, readiness,",
        "recovery, and external comparison remain open."
      )
    )
  }

  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "The parameter-path reachability operator applies only to GPCM.",
      "parameter_path_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional free-coordinate operator is not reused for MML.",
      "parameter_path_not_applicable_estimator"
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
  if (!isTRUE(source_matches)) {
    return(finish(
      "not_evaluated_source_contract", FALSE,
      "The required lexicographic-limit result was unavailable.",
      "parameter_path_lexicographic_source_mismatch"
    ))
  }
  max_nonzeros <-
    mfrmr_jml_gpcm_exponential_balance_scalar(max_design_nonzeros)
  if (!is.finite(max_nonzeros) || max_nonzeros < 0 ||
      max_nonzeros != floor(max_nonzeros)) {
    return(finish(
      "not_evaluated_control", FALSE,
      "The sparse operator workload must be a finite nonnegative integer.",
      "parameter_path_control_invalid"
    ))
  }
  adjacent <- tryCatch(
    mfrmr_estimability_adjacent_design(
      prep, idx, config, sizes,
      include_person = TRUE, include_population_beta = FALSE
    ),
    error = function(e) e
  )
  if (inherits(adjacent, "error")) {
    return(finish(
      "not_evaluated_design", FALSE,
      paste0("The retained additive operator could not be constructed: ",
             conditionMessage(adjacent)),
      "parameter_path_retained_design_invalid"
    ))
  }
  design <- methods::as(adjacent$design, "dgCMatrix")
  if (length(design@x) > max_nonzeros) {
    return(finish(
      "not_evaluated_size_limit", FALSE,
      "The retained additive operator exceeded the sparse workload limit.",
      "parameter_path_design_nonzero_limit",
      design = design,
      map = adjacent$map
    ))
  }
  finish(
    "retained_additive_parameter_operator_available_no_path_declared", TRUE,
    paste(
      "The constrained retained additive operator and coordinate map are",
      "available for caller-declared forward path construction."
    ),
    paste(
      "parameter_path_forward_operator_available",
      "production_path_not_declared",
      "global_boundary_not_classified",
      sep = ";"
    ),
    available = TRUE,
    design = design,
    map = adjacent$map
  )
}
