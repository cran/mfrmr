# JML GPCM asymptotically-affine boundary-certificate transport
# ==============================================================================
#
# A certified constant-rate path also certifies a narrow curved neighborhood.
# If every additive coordinate equals its certified affine path plus a residual
# converging to zero, and every expanded log slope equals its certified affine
# path plus a sum-zero residual converging to zero, strict positive/neutral/
# leading-negative inequalities persist and the likelihood has the same
# analytic boundary limit. This is a positive transport theorem only.
#
# Bounded residuals that do not vanish are deliberately excluded. In
# particular, a zero-rate slope can then oscillate among different finite
# values and change the limiting likelihood. Curved paths with nonconvergent
# normalized rates also remain outside this audit.

mfrmr_jml_gpcm_asymptotically_affine_contract_version <- function() {
  "mfrmr-jml-gpcm-asymptotically-affine-transport-0.2.3-v1"
}

mfrmr_jml_gpcm_asymptotically_affine_empty_registry <- function() {
  data.frame(
    PathId = character(0),
    SourceFamily = character(0),
    SourceState = character(0),
    BoundaryLogLikelihood = numeric(0),
    BoundaryImprovement = numeric(0),
    PositiveMinimumMargin = numeric(0),
    LeadingNegativeCoefficient = numeric(0),
    StrictCertificate = logical(0),
    SourceComplete = logical(0),
    TransportCertified = logical(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_asymptotically_affine_rows <- function(
    certificates,
    source_family,
    source_state,
    source_complete,
    id_column,
    margin_column = NULL,
    leading_column = NULL,
    strict_rows_column = NULL) {
  certificates <- as.data.frame(
    certificates %||% data.frame(), stringsAsFactors = FALSE
  )
  empty <- mfrmr_jml_gpcm_asymptotically_affine_empty_registry()
  if (nrow(certificates) == 0L ||
      !"Certified" %in% names(certificates)) {
    return(empty)
  }
  certificates <- certificates[
    !is.na(certificates$Certified) & certificates$Certified,
    , drop = FALSE
  ]
  if (nrow(certificates) == 0L) return(empty)

  get_numeric <- function(column, default = NA_real_) {
    if (is.null(column) || !column %in% names(certificates)) {
      return(rep(default, nrow(certificates)))
    }
    suppressWarnings(as.numeric(certificates[[column]]))
  }
  path_id <- if (id_column %in% names(certificates)) {
    as.character(certificates[[id_column]])
  } else {
    sprintf("PATH%05d", seq_len(nrow(certificates)))
  }
  boundary <- get_numeric("BoundaryLogLikelihood")
  improvement <- get_numeric("BoundaryImprovement")
  margin <- get_numeric(margin_column)
  leading <- get_numeric(leading_column)
  strict_rows <- get_numeric(strict_rows_column, default = 0)
  strict <- if (!is.null(strict_rows_column)) {
    is.finite(strict_rows) & strict_rows > 0
  } else {
    is.finite(margin) & margin > 0 &
      is.finite(leading) & leading > 0
  }
  transport <- strict & is.finite(boundary) & is.finite(improvement)
  data.frame(
    PathId = path_id,
    SourceFamily = rep(source_family, nrow(certificates)),
    SourceState = rep(source_state, nrow(certificates)),
    BoundaryLogLikelihood = boundary,
    BoundaryImprovement = improvement,
    PositiveMinimumMargin = margin,
    LeadingNegativeCoefficient = leading,
    StrictCertificate = strict,
    SourceComplete = rep(isTRUE(source_complete), nrow(certificates)),
    TransportCertified = transport,
    stringsAsFactors = FALSE
  )
}

audit_mfrm_jml_gpcm_asymptotically_affine_transport <- function(
    config,
    slope_audit = config$boundary_audit$gpcm_slope_boundary %||% list(),
    joint_audit = config$boundary_audit$gpcm_joint_boundary %||% list(),
    boundary_classification =
      config$boundary_audit$gpcm_fixed_objective_classification %||% list()) {
  method <- as.character(config$method %||% "unknown")[1]
  model <- as.character(config$model %||% "unknown")[1]
  empty <- mfrmr_jml_gpcm_asymptotically_affine_empty_registry()
  finish <- function(state, evaluated, detail, reason_codes,
                     registry = empty) {
    registry <- as.data.frame(registry, stringsAsFactors = FALSE)
    transported <- if (nrow(registry) > 0L) {
      sum(registry$TransportCertified, na.rm = TRUE)
    } else {
      0L
    }
    list(
      contract_version =
        mfrmr_jml_gpcm_asymptotically_affine_contract_version(),
      method = method,
      model = model,
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      source_boundary_contract = as.character(
        boundary_classification$contract_version %||% "missing"
      )[1],
      source_boundary_state = as.character(
        boundary_classification$state %||% "missing"
      )[1],
      state = state,
      evaluated = isTRUE(evaluated),
      source_path_count = as.integer(nrow(registry)),
      transported_path_count = as.integer(transported),
      path_registry = registry,
      additive_path_condition = paste(
        "certified additive affine path plus a residual whose",
        "free-coordinate sup norm converges to zero"
      ),
      log_slope_path_condition = paste(
        "certified expanded log-slope affine path plus a sum-zero residual",
        "whose sup norm converges to zero"
      ),
      same_boundary_limit_certified = transported > 0L,
      positive_transport_theorem_applied = transported > 0L,
      vanishing_residual_curved_neighborhood_classified =
        transported > 0L,
      bounded_nonvanishing_residual_classified = FALSE,
      rate_nonconvergent_path_classified = FALSE,
      general_curved_path_classified = FALSE,
      curved_path_absence_certified = FALSE,
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
        "This is a positive continuity transport for already certified",
        "constant-rate paths. A negative source search or failure to",
        "transport a path does not exclude curved boundaries. Bounded but",
        "nonvanishing residuals can change or destroy the boundary limit,",
        "especially for zero-rate slopes. Paths with nonconvergent normalized",
        "rates, unrestricted nonlinear additive motion, and MML marginal",
        "geometry remain unclassified."
      )
    )
  }

  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "The asymptotically-affine transport applies only to GPCM.",
      "asymptotically_affine_transport_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional JML path theorem is not reused for MML.",
      "asymptotically_affine_transport_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (as.integer(config$gpcm_spec$n_params %||% 0L) == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction has no free log-slope boundary path.",
      "no_free_log_slope_coordinate"
    ))
  }
  identity_matches <- identical(
    as.character(boundary_classification$contract_version %||% "")[1],
    mfrmr_jml_gpcm_boundary_classification_contract_version()
  ) && identical(
    as.character(boundary_classification$objective_identity %||% "")[1],
    "identified_conditional_joint_log_likelihood"
  )
  if (!identity_matches) {
    return(finish(
      "not_evaluated_objective_identity", FALSE,
      "The source fixed-objective boundary classification is missing or mismatched.",
      "asymptotically_affine_source_identity_mismatch"
    ))
  }

  slope_rows <- mfrmr_jml_gpcm_asymptotically_affine_rows(
    slope_audit$certificates,
    source_family = "slope_only_constant_rate",
    source_state = as.character(slope_audit$state %||% "missing")[1],
    source_complete = isTRUE(slope_audit$scope_complete),
    id_column = "PairId",
    strict_rows_column = "StrictRows"
  )
  pair_rows <- mfrmr_jml_gpcm_asymptotically_affine_rows(
    joint_audit$certificates,
    source_family = "joint_ordered_pair",
    source_state = as.character(joint_audit$state %||% "missing")[1],
    source_complete = isTRUE(joint_audit$scope_complete),
    id_column = "PairId",
    margin_column = "PositiveMinimumMargin",
    leading_column = "NegativeLeadingCoefficient"
  )
  rate_rows <- mfrmr_jml_gpcm_asymptotically_affine_rows(
    joint_audit$rate_certificates,
    source_family = "joint_canonical_constant_rate",
    source_state = as.character(joint_audit$rate_scope_state %||% "missing")[1],
    source_complete = isTRUE(joint_audit$rate_scope_complete),
    id_column = "RateId",
    margin_column = "PositiveMinimumMargin",
    leading_column = "LeadingNegativeCoefficient"
  )
  registry <- rbind(slope_rows, pair_rows, rate_rows)
  transported <- if (nrow(registry) > 0L) {
    any(registry$TransportCertified)
  } else {
    FALSE
  }
  source_positive <-
    isTRUE(boundary_classification$positive_boundary_certificate_present)
  if (transported) {
    return(finish(
      "certified_vanishing_residual_curved_neighborhood", TRUE,
      paste(
        "At least one fixed-objective constant-rate boundary certificate",
        "transports to every declared vanishing-residual asymptotically",
        "affine perturbation with the same source direction and boundary",
        "limit."
      ),
      paste(
        "positive_constant_rate_certificate_transported",
        "bounded_nonvanishing_residuals_not_classified",
        "general_curved_paths_not_classified",
        sep = ";"
      ),
      registry = registry
    ))
  }
  if (source_positive) {
    return(finish(
      "indeterminate_source_certificate", FALSE,
      paste(
        "The fixed-objective classifier reports a positive boundary, but",
        "no strict finite source certificate could be transported."
      ),
      "asymptotically_affine_source_certificate_incoherent",
      registry = registry
    ))
  }
  if (isTRUE(
    boundary_classification$no_path_certified_in_audited_families
  )) {
    return(finish(
      "no_positive_path_to_transport", TRUE,
      paste(
        "The admitted constant-rate source families completed without a",
        "positive certificate. This supplies no negative conclusion for",
        "curved paths."
      ),
      paste(
        "no_positive_constant_rate_certificate",
        "curved_path_absence_not_certified",
        sep = ";"
      ),
      registry = registry
    ))
  }
  finish(
    "not_evaluated_source", FALSE,
    paste(
      "The fixed-objective source boundary families were incomplete or",
      "indeterminate, so no curve transport conclusion is available."
    ),
    "asymptotically_affine_source_boundary_incomplete",
    registry = registry
  )
}
