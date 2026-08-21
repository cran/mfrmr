# Estimator-specific Person boundary audit
# ==============================================================================
#
# A finite optimizer iterate is not the primary JML estimate when an otherwise
# free Person has an extreme sufficient score.  This audit separates that
# numerical trace from the statistical result.  It deliberately covers Person
# coordinates only.  Separation of constrained non-Person elements and
# interactions requires a recession-direction audit and is not inferred from a
# response-constant level alone.

mfrmr_boundary_contract_version <- function() {
  mfrmr_readiness_contract_version()
}

mfrmr_person_constraint_rows <- function(config, persons) {
  spec <- config$theta_spec
  persons <- as.character(persons)
  if (is.null(spec)) {
    return(data.frame(
      Person = persons,
      DirectFixed = FALSE,
      ConstraintFixed = FALSE,
      ConstraintCoupled = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  spec_levels <- as.character(spec$levels %||% character(0))
  level_index <- match(persons, spec_levels)
  anchors <- as.numeric(spec$anchors %||% rep(NA_real_, length(spec_levels)))
  jacobian <- constraint_jacobian(spec)
  row_active <- rep(FALSE, length(persons))
  if (ncol(jacobian) > 0L) {
    present <- !is.na(level_index)
    row_active[present] <- apply(
      abs(jacobian[level_index[present], , drop = FALSE]) > 1e-12,
      1L,
      any
    )
  }
  direct_fixed <- rep(FALSE, length(persons))
  present <- !is.na(level_index)
  direct_fixed[present] <- is.finite(anchors[level_index[present]])
  constraint_fixed <- present & !row_active

  groups <- as.character(spec$groups %||% rep(NA_character_, length(spec_levels)))
  grouped <- rep(FALSE, length(persons))
  grouped[present] <- !is.na(groups[level_index[present]]) &
    nzchar(groups[level_index[present]])
  constraint_coupled <- row_active & (isTRUE(spec$centered) | grouped)

  data.frame(
    Person = persons,
    DirectFixed = direct_fixed,
    ConstraintFixed = constraint_fixed,
    ConstraintCoupled = constraint_coupled,
    stringsAsFactors = FALSE
  )
}

audit_mfrm_person_boundary <- function(person_table, config, prep) {
  tbl <- as.data.frame(person_table, stringsAsFactors = FALSE)
  if (!all(c("Person", "Estimate") %in% names(tbl))) {
    stop("Internal Person boundary audit received a malformed table.",
         call. = FALSE)
  }

  method <- as.character(config$method %||% NA_character_)
  response_flag <- .flag_extreme_persons(
    prep,
    rating_min = prep$rating_min,
    rating_max = prep$rating_max
  )
  response_extreme <- tolower(as.character(
    response_flag[match(as.character(tbl$Person), names(response_flag))]
  ))
  response_extreme[is.na(response_extreme) | !nzchar(response_extreme)] <- "none"
  response_extreme[!response_extreme %in% c("none", "low", "high")] <- "none"
  constraints <- mfrmr_person_constraint_rows(config, tbl$Person)

  n <- nrow(tbl)
  response_rows <- integer(n)
  weighted_response_total <- numeric(n)
  prep_data <- as.data.frame(prep$data %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(prep_data) > 0L && "Person" %in% names(prep_data)) {
    person_index <- match(as.character(prep_data$Person), as.character(tbl$Person))
    contributing <- !is.na(person_index)
    response_rows <- tabulate(person_index[contributing], nbins = n)
    weights <- if ("Weight" %in% names(prep_data)) {
      suppressWarnings(as.numeric(prep_data$Weight))
    } else {
      rep(1, nrow(prep_data))
    }
    weights[!is.finite(weights) | weights <= 0] <- 0
    weighted_by_person <- tapply(
      weights[contributing],
      person_index[contributing],
      sum
    )
    if (length(weighted_by_person) > 0L) {
      weighted_index <- suppressWarnings(as.integer(names(weighted_by_person)))
      keep_weight <- is.finite(weighted_index) &
        weighted_index >= 1L & weighted_index <= n
      weighted_response_total[weighted_index[keep_weight]] <-
        as.numeric(weighted_by_person[keep_weight])
    }
  }
  status <- rep("estimable", n)
  reason <- rep("", n)
  primary <- as.numeric(tbl$Estimate)
  optimizer_estimate <- rep(NA_real_, n)
  display_estimate <- primary
  boundary_direction <- rep("none", n)
  primary_basis <- if (identical(method, "MML")) {
    rep("posterior_eap", n)
  } else {
    rep("jml_finite_mle", n)
  }
  optimizer_use <- rep("not_applicable", n)

  extreme <- response_extreme %in% c("low", "high")
  if (identical(method, "MML")) {
    reason[extreme] <- "mml_extreme_response_prior_regularized"
    boundary_direction[extreme] <- response_extreme[extreme]
  } else if (identical(method, "JML")) {
    optimizer_estimate <- primary
    optimizer_use[] <- "same_as_primary"

    fixed <- constraints$DirectFixed | constraints$ConstraintFixed
    status[fixed] <- "fixed"
    primary_basis[fixed] <- "fixed_constraint"
    optimizer_estimate[fixed] <- NA_real_
    optimizer_use[fixed] <- "not_optimized"
    reason[fixed & extreme] <- "fixed_extreme_response"
    boundary_direction[fixed & extreme] <- response_extreme[fixed & extreme]

    coupled_extreme <- extreme & !fixed & constraints$ConstraintCoupled
    status[coupled_extreme] <- "weak_information"
    reason[coupled_extreme] <- "weak_design_information"
    primary_basis[coupled_extreme] <- "constraint_coupled_jml_review"
    optimizer_use[coupled_extreme] <- "review_only"
    boundary_direction[coupled_extreme] <- response_extreme[coupled_extreme]

    free_extreme <- extreme & !fixed & !constraints$ConstraintCoupled
    low <- free_extreme & response_extreme == "low"
    high <- free_extreme & response_extreme == "high"
    status[low] <- "unbounded_low"
    status[high] <- "unbounded_high"
    reason[low] <- "jml_extreme_low"
    reason[high] <- "jml_extreme_high"
    primary[low] <- -Inf
    primary[high] <- Inf
    display_estimate[free_extreme] <- NA_real_
    primary_basis[free_extreme] <- "jml_extreme_sufficient_score_boundary"
    optimizer_use[free_extreme] <- "numerical_trace_only"
    boundary_direction[free_extreme] <- response_extreme[free_extreme]
  } else {
    stop("Internal Person boundary audit requires method JML or MML.",
         call. = FALSE)
  }

  parameter_id <- if (identical(method, "MML")) {
    paste0("Person:", tbl$Person, ":EAP")
  } else {
    paste0("Person:", tbl$Person)
  }
  parameter_status <- data.frame(
    ReadinessContractVersion = mfrmr_boundary_contract_version(),
    ReadinessScope = "parameter",
    ParameterId = parameter_id,
    ParameterClass = "Person",
    Facet = "Person",
    Level = as.character(tbl$Person),
    PrimaryEstimate = primary,
    OptimizerEstimate = optimizer_estimate,
    DisplayEstimate = display_estimate,
    DisplayAdjustment = "none",
    ParameterStatus = status,
    BoundaryDirection = boundary_direction,
    ResponseExtreme = response_extreme,
    ResponseRows = as.integer(response_rows),
    WeightedResponseTotal = as.numeric(weighted_response_total),
    PrimaryEstimateBasis = primary_basis,
    OptimizerEstimateUse = optimizer_use,
    ReasonCodes = reason,
    AuditProvenance = "person_sufficient_score_and_constraint_jacobian_v1",
    stringsAsFactors = FALSE
  )

  unbounded <- status %in% c("unbounded_low", "unbounded_high")
  coupled_review <- status == "weak_information" & extreme
  boundary_state <- if (any(coupled_review)) {
    "not_evaluated"
  } else if (any(unbounded)) {
    "has_exclusions"
  } else {
    "finite"
  }
  fit_reason_rows <- unbounded | coupled_review
  fit_reasons <- unique(reason[fit_reason_rows & nzchar(reason)])
  readiness <- data.frame(
    ReadinessContractVersion = mfrmr_boundary_contract_version(),
    ReadinessScope = "fit",
    BoundaryState = boundary_state,
    ReasonCodes = paste(fit_reasons, collapse = ";"),
    Complete = !any(coupled_review),
    AuditedParameterClass = "Person",
    UnboundedLowN = sum(status == "unbounded_low"),
    UnboundedHighN = sum(status == "unbounded_high"),
    FixedExtremeN = sum(status == "fixed" & extreme),
    PriorRegularizedExtremeN = sum(
      reason == "mml_extreme_response_prior_regularized"
    ),
    ConstraintCoupledExtremeN = sum(coupled_review),
    stringsAsFactors = FALSE
  )

  list(
    contract_version = mfrmr_boundary_contract_version(),
    method = method,
    model = as.character(config$model %||% NA_character_),
    scope = "Person",
    complete = !any(coupled_review),
    readiness = readiness,
    parameter_status = parameter_status,
    limitations = paste(
      "This audit classifies Person response boundaries only.",
      "Non-Person facets and interactions require a separate constrained",
      "recession-direction audit."
    )
  )
}

apply_mfrm_person_boundary <- function(person_table, boundary_audit) {
  tbl <- as.data.frame(person_table, stringsAsFactors = FALSE)
  status <- as.data.frame(
    boundary_audit$parameter_status %||% data.frame(),
    stringsAsFactors = FALSE
  )
  if (nrow(tbl) == 0L || nrow(status) == 0L) return(tbl)
  index <- match(as.character(tbl$Person), status$Level)
  if (anyNA(index)) {
    stop("Internal Person boundary audit could not align Person levels.",
         call. = FALSE)
  }
  tbl$Estimate <- status$PrimaryEstimate[index]
  tbl$PrimaryEstimate <- status$PrimaryEstimate[index]
  tbl$OptimizerEstimate <- status$OptimizerEstimate[index]
  tbl$DisplayEstimate <- status$DisplayEstimate[index]
  tbl$DisplayAdjustment <- status$DisplayAdjustment[index]
  tbl$ParameterStatus <- status$ParameterStatus[index]
  tbl$BoundaryDirection <- status$BoundaryDirection[index]
  tbl$ResponseExtreme <- status$ResponseExtreme[index]
  tbl$Extreme <- status$ResponseExtreme[index]
  tbl$ResponseRows <- status$ResponseRows[index]
  tbl$WeightedResponseTotal <- status$WeightedResponseTotal[index]
  tbl$PrimaryEstimateBasis <- status$PrimaryEstimateBasis[index]
  tbl$OptimizerEstimateUse <- status$OptimizerEstimateUse[index]
  tbl$ReasonCodes <- status$ReasonCodes[index]
  tbl$ReadinessContractVersion <- status$ReadinessContractVersion[index]
  tbl
}
