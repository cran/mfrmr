# Fixed-objective JML GPCM boundary classification
# ==============================================================================
#
# The slope-only and joint additive/slope audits answer different sufficient
# questions.  This classifier combines their typed outcomes without turning a
# negative bounded search into a global finite-MLE claim.  The objective
# identity is part of the result because penalized JML, finite-box JML, and MML
# can select different finite points even on the same response model.

mfrmr_jml_gpcm_boundary_classification_contract_version <- function() {
  "mfrmr-jml-gpcm-fixed-objective-boundary-0.2.3-v2"
}

mfrmr_jml_gpcm_boundary_audit_objective_finite <- function(audit) {
  retained <- suppressWarnings(as.numeric(
    audit$retained_log_likelihood %||% NA_real_
  )[1])
  optimizer <- suppressWarnings(as.numeric(
    audit$optimizer_log_likelihood %||% NA_real_
  )[1])
  is.finite(retained) && is.finite(optimizer)
}

mfrmr_jml_gpcm_boundary_audit_numerically_indeterminate <- function(audit) {
  state <- as.character(audit$state %||% "")[1]
  if (state %in% c(
    "not_evaluated_numerical",
    "not_evaluated_likelihood_mismatch"
  )) {
    return(TRUE)
  }
  if (!state %in% c(
    "not_evaluated_pair_solver", "not_evaluated_rate_solver"
  )) return(FALSE)
  certificates <- if (identical(state, "not_evaluated_rate_solver")) {
    as.data.frame(
      audit$rate_certificates %||% data.frame(), stringsAsFactors = FALSE
    )
  } else {
    as.data.frame(
      audit$certificates %||% data.frame(), stringsAsFactors = FALSE
    )
  }
  if (nrow(certificates) == 0L ||
      !"ReasonCodes" %in% names(certificates)) {
    return(TRUE)
  }
  failure_reason <- if (identical(state, "not_evaluated_rate_solver")) {
    "joint_rate_linear_program_failed"
  } else {
    "joint_pair_linear_program_failed"
  }
  any(as.character(certificates$ReasonCodes) == failure_reason)
}

mfrmr_classify_jml_gpcm_fixed_objective_boundary <- function(
    config,
    slope_audit = config$boundary_audit$gpcm_slope_boundary %||% list(),
    joint_audit = config$boundary_audit$gpcm_joint_boundary %||% list()) {
  method <- as.character(config$method %||% "unknown")[1]
  model <- as.character(config$model %||% "unknown")[1]
  slope_state <- as.character(slope_audit$state %||% "missing")[1]
  joint_state <- as.character(joint_audit$state %||% "missing")[1]
  free_log_slopes <- as.integer(config$gpcm_spec$n_params %||% 0L) > 0L

  base <- list(
    contract_version =
      mfrmr_jml_gpcm_boundary_classification_contract_version(),
    method = method,
    model = model,
    estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
    objective_identity = "identified_conditional_joint_log_likelihood",
    person_treatment = "fixed_unknown_coordinates_jointly_estimated",
    slope_identification = "sum_zero_expanded_log_slopes",
    penalty = "none",
    finite_parameter_box = FALSE,
    state = "not_evaluated",
    slope_audit_state = slope_state,
    joint_audit_state = joint_state,
    retained_point_state = "not_evaluated",
    retained_point_finite = FALSE,
    slope_only_recession_certified = FALSE,
    competitive_joint_boundary_certified = FALSE,
    canonical_constant_rate_scope_complete = FALSE,
    canonical_constant_rate_boundary_certified = FALSE,
    positive_boundary_certificate_present = FALSE,
    path_family_classification_complete = FALSE,
    audited_path_family_result_classified = FALSE,
    fixed_objective_boundary_presence_classified = FALSE,
    no_path_certified_in_audited_families = FALSE,
    numerical_indeterminacy = FALSE,
    global_boundary_classified = FALSE,
    global_finite_maximum_certified = FALSE,
    global_boundary_absence_certified = FALSE,
    standard_error_eligible = FALSE,
    confidence_interval_eligible = FALSE,
    external_comparison_eligible = FALSE,
    readiness_effect = "none_diagnostic_only",
    reason_codes = "fixed_objective_boundary_not_evaluated",
    detail = paste(
      "The estimator-specific JML GPCM boundary classifier was not",
      "evaluated."
    ),
    limitations = paste(
      "The classifier combines the complete constant slope-only ray check",
      "with a workload-bounded joint linear-additive family whose canonical",
      "constant rates include positive, zero, leading-negative, and",
      "deeper-negative groups. Failure to certify these families does not",
      "establish a finite global maximum, boundary absence, nonlinear",
      "structural identification, uncertainty, or external-comparison",
      "eligibility. Curved and rate-nonconvergent paths remain open. The",
      "result does not apply to penalized JML, finite-box JML, or MML."
    )
  )

  if (!identical(model, "GPCM")) {
    base$state <- "not_applicable_model"
    base$estimator_identity <- "not_applicable"
    base$objective_identity <- "not_applicable"
    base$person_treatment <- "not_applicable"
    base$slope_identification <- "not_applicable"
    base$penalty <- "not_applicable"
    base$finite_parameter_box <- NA
    base$reason_codes <- "fixed_objective_boundary_not_applicable_model"
    base$detail <- "The classifier applies only to GPCM."
    return(base)
  }
  if (!identical(method, "JML")) {
    base$state <- "not_applicable_estimator"
    base$estimator_identity <- "not_applicable"
    base$objective_identity <- "not_applicable_marginal_estimator"
    base$person_treatment <- "integrated_or_estimator_specific"
    base$slope_identification <- as.character(
      config$gpcm_spec$identification %||% "unknown"
    )[1]
    base$penalty <- "not_applicable"
    base$finite_parameter_box <- NA
    base$reason_codes <- "fixed_objective_boundary_not_applicable_estimator"
    base$detail <- paste(
      "The conditional fixed-effects JML classifier is not reused for MML",
      "or another estimator."
    )
    return(base)
  }
  if (!free_log_slopes) {
    base$state <- "not_required_unit_slope"
    base$retained_point_state <- "identified_unit_slope_reduction"
    base$path_family_classification_complete <- TRUE
    base$audited_path_family_result_classified <- TRUE
    base$reason_codes <- "no_free_log_slope_coordinate"
    base$detail <- paste(
      "The unit-slope GPCM reduction has no free relative log-slope",
      "coordinate requiring this classifier."
    )
    return(base)
  }

  objective_consistent <- identical(
    as.character(slope_audit$method %||% "")[1], "JML"
  ) && identical(
    as.character(slope_audit$model %||% "")[1], "GPCM"
  ) && identical(
    as.character(joint_audit$method %||% "")[1], "JML"
  ) && identical(
    as.character(joint_audit$model %||% "")[1], "GPCM"
  )
  if (!objective_consistent) {
    base$state <- "not_evaluated_objective_identity"
    base$reason_codes <- "boundary_audit_objective_identity_missing_or_mismatched"
    base$detail <- paste(
      "Both component audits must identify the JML GPCM objective before",
      "their results can be combined."
    )
    return(base)
  }

  slope_scope_complete <- isTRUE(slope_audit$complete) &&
    isTRUE(slope_audit$scope_complete)
  joint_scope_complete <- isTRUE(joint_audit$complete) &&
    isTRUE(joint_audit$scope_complete)
  slope_state_classified <- slope_state %in% c(
    "certified_monotone_boundary_path", "none_certified"
  )
  joint_state_classified <- joint_state %in% c(
    "certified_competitive_joint_boundary_path",
    "none_certified",
    "no_free_additive_coordinates"
  )
  slope_certified <- slope_scope_complete && identical(
    slope_state, "certified_monotone_boundary_path"
  )
  joint_certified <- identical(
    joint_state, "certified_competitive_joint_boundary_path"
  ) && if (!is.null(joint_audit$positive_certificate_present)) {
    isTRUE(joint_audit$positive_certificate_present)
  } else {
    joint_scope_complete
  }
  path_complete <- slope_scope_complete && joint_scope_complete &&
    slope_state_classified && joint_state_classified
  retained_finite <-
    mfrmr_jml_gpcm_boundary_audit_objective_finite(slope_audit) ||
    mfrmr_jml_gpcm_boundary_audit_objective_finite(joint_audit)
  numerically_indeterminate <-
    mfrmr_jml_gpcm_boundary_audit_numerically_indeterminate(slope_audit) ||
    mfrmr_jml_gpcm_boundary_audit_numerically_indeterminate(joint_audit)

  base$retained_point_finite <- retained_finite
  base$retained_point_state <- if (retained_finite) {
    "finite_optimizer_trace_not_finite_maximum_claim"
  } else {
    "not_evaluated"
  }
  base$slope_only_recession_certified <- slope_certified
  base$competitive_joint_boundary_certified <- joint_certified
  base$canonical_constant_rate_scope_complete <-
    isTRUE(joint_audit$rate_scope_complete %||% joint_scope_complete)
  base$canonical_constant_rate_boundary_certified <-
    joint_certified && any(as.logical(
      as.data.frame(
        joint_audit$rate_certificates %||% data.frame(Certified = FALSE),
        stringsAsFactors = FALSE
      )$Certified %||% FALSE
    ), na.rm = TRUE)
  base$positive_boundary_certificate_present <-
    slope_certified || joint_certified
  base$path_family_classification_complete <- path_complete
  base$audited_path_family_result_classified <-
    slope_certified || joint_certified || path_complete
  base$fixed_objective_boundary_presence_classified <-
    slope_certified || joint_certified
  base$no_path_certified_in_audited_families <-
    path_complete && !slope_certified && !joint_certified
  base$numerical_indeterminacy <- numerically_indeterminate

  if (slope_certified && joint_certified) {
    base$state <- "certified_slope_and_joint_boundary"
    base$reason_codes <- paste(
      "certified_monotone_slope_only_recession",
      "certified_competitive_joint_boundary_candidate",
      "global_boundary_not_classified",
      sep = ";"
    )
    base$detail <- paste(
      "A monotone slope-only recession path and a competitive joint",
      "additive/constant-log-slope-rate boundary path were certified for the same",
      "unpenalized conditional JML objective."
    )
  } else if (slope_certified) {
    base$state <- "certified_slope_only_recession"
    base$reason_codes <- paste(
      "certified_monotone_slope_only_recession",
      "global_boundary_not_classified",
      sep = ";"
    )
    base$detail <- paste(
      "A monotone slope-only recession path was certified for the",
      "unpenalized conditional JML objective."
    )
  } else if (joint_certified) {
    base$state <- "certified_competitive_joint_boundary"
    base$reason_codes <- paste(
      "certified_competitive_joint_boundary_candidate",
      "global_status_open",
      sep = ";"
    )
    base$detail <- paste(
      "A competitive joint additive/constant-log-slope-rate asymptotic path was",
      "certified. It is a sufficient warning, not a global proof for the",
      "non-concave GPCM likelihood."
    )
  } else if (numerically_indeterminate) {
    base$state <- "indeterminate_numerical"
    base$audited_path_family_result_classified <- FALSE
    base$fixed_objective_boundary_presence_classified <- FALSE
    base$no_path_certified_in_audited_families <- FALSE
    base$reason_codes <- paste(
      "boundary_audit_numerically_indeterminate",
      "global_boundary_not_classified",
      sep = ";"
    )
    base$detail <- paste(
      "A likelihood reconstruction or numerical path computation was",
      "indeterminate, so absence of a boundary path is not classified."
    )
  } else if (path_complete) {
    base$state <- "finite_retained_point_no_path_certified"
    base$reason_codes <- paste(
      "no_path_certified_in_bounded_audited_families",
      "finite_optimizer_trace_not_finite_maximum_claim",
      "global_boundary_not_classified",
      sep = ";"
    )
    base$detail <- paste(
      "Both bounded path families were evaluated without a positive",
      "certificate. The retained optimizer point is finite, but no finite",
      "global maximum or global boundary absence is established."
    )
  } else {
    base$state <- "not_evaluated"
    base$audited_path_family_result_classified <- FALSE
    base$fixed_objective_boundary_presence_classified <- FALSE
    base$no_path_certified_in_audited_families <- FALSE
    base$reason_codes <- paste(
      "bounded_path_family_audit_incomplete",
      "global_boundary_not_classified",
      sep = ";"
    )
    base$detail <- paste(
      "At least one required bounded path-family audit was not completed,",
      "so boundary absence is not classified."
    )
  }
  base
}
