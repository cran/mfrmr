# Higher-order owner-rate lifts for zero-offset JML GPCM response images
# ==============================================================================
#
# P2j gives every finite closure point a nonnegative inverse-slope simplex-face
# witness.  A first-order lift can be insufficient because different inactive
# slope owners may have to vanish at different rates.  For an ordered owner
# rate hierarchy r_g and positive coefficients c_g, set
#
#   w_g(epsilon) = c_g epsilon^r_g,
#   beta(epsilon) = sum_k epsilon^k beta_k.
#
# On rows owned by g, every A beta_k below r_g must vanish and A beta_r_g must
# equal c_g z.  These are linear conditions.  A feasible hierarchy gives an
# explicit finite path.  Conversely, semialgebraic curve selection and Puiseux
# leading terms reduce every zero-offset finite-image closure sequence to one
# of the finitely many ordered owner-rate hierarchies.  Thus exhaustive,
# numerically decided hierarchy enumeration classifies targetwise closure
# membership, while it does not yet construct a symbolic global stratification
# or a likelihood envelope.

mfrmr_jml_gpcm_higher_order_lift_contract_version <- function() {
  "mfrmr-jml-gpcm-higher-order-face-lifts-0.2.3-v1"
}

mfrmr_jml_gpcm_higher_order_lift_empty_table <- function() {
  data.frame(
    HierarchyId = character(0),
    ZeroStageSupport = character(0),
    RateCode = character(0),
    PositiveStages = integer(0),
    MaximumRate = integer(0),
    FirstOrderHierarchy = logical(0),
    SourceFaceFirstOrderReachable = logical(0),
    SolverStatus = integer(0),
    MaximumPositiveMargin = numeric(0),
    MaximumAbsoluteEquationResidual = numeric(0),
    SumCoefficientResidual = numeric(0),
    HierarchyLiftCertified = logical(0),
    StrictHierarchyExcluded = logical(0),
    NumericalMarginOpen = logical(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_higher_order_lift_partition_count <- function(
    elements,
    max_stages = elements) {
  n <- as.integer(elements)
  stage_cap <- min(as.integer(max_stages), n)
  if (n < 1L || stage_cap < 1L) {
    return(0)
  }
  stirling <- matrix(0, nrow = n + 1L, ncol = n + 1L)
  stirling[1, 1] <- 1
  for (row in seq_len(n)) {
    for (stage in seq_len(row)) {
      stirling[row + 1L, stage + 1L] <-
        stage * stirling[row, stage + 1L] +
        stirling[row, stage]
    }
  }
  sum(vapply(seq_len(stage_cap), function(stage) {
    factorial(stage) * stirling[n + 1L, stage + 1L]
  }, numeric(1)))
}

mfrmr_jml_gpcm_higher_order_lift_surjections <- function(
    elements,
    stages) {
  n <- as.integer(elements)
  stage_count <- as.integer(stages)
  output <- list()
  current <- integer(n)
  used <- logical(stage_count)
  output_index <- 0L

  recurse <- function(position) {
    remaining <- n - position + 1L
    if (sum(!used) > remaining) {
      return(invisible(NULL))
    }
    if (position > n) {
      if (all(used)) {
        output_index <<- output_index + 1L
        output[[output_index]] <<- current
      }
      return(invisible(NULL))
    }
    for (stage in seq_len(stage_count)) {
      was_used <- used[stage]
      current[position] <<- stage
      used[stage] <<- TRUE
      recurse(position + 1L)
      used[stage] <<- was_used
    }
    invisible(NULL)
  }
  recurse(1L)
  output
}

mfrmr_jml_gpcm_higher_order_lift_lp <- function(
    design,
    target,
    owner,
    rates,
    tolerance,
    max_equations,
    max_lp_nonzeros,
    lp_timeout) {
  rows <- nrow(design)
  coordinates <- ncol(design)
  owners <- max(owner)
  rates <- as.integer(rates)
  maximum_rate <- max(rates)
  row_rates <- rates[owner]
  equations <- sum(row_rates + 1L)
  coefficient_blocks <- maximum_rate + 1L
  beta_block <- coordinates * coefficient_blocks
  margin_index <- owners + 2L * beta_block + 1L
  variables <- margin_index
  if (equations > max_equations) {
    return(list(
      evaluated = FALSE,
      status = NA_integer_,
      certified = FALSE,
      reason = "higher_order_lift_equation_limit"
    ))
  }

  equation_index <- matrix(
    NA_integer_,
    nrow = rows,
    ncol = coefficient_blocks
  )
  counter <- 0L
  for (row in seq_len(rows)) {
    count <- row_rates[row] + 1L
    equation_index[row, seq_len(count)] <- counter + seq_len(count)
    counter <- counter + count
  }
  design_summary <- Matrix::summary(design)
  triplet_parts <- vector("list", 2L * coefficient_blocks + 3L)
  part <- 0L
  for (stage in 0:maximum_rate) {
    keep <- row_rates[design_summary$i] >= stage
    if (any(keep)) {
      selected <- design_summary[keep, , drop = FALSE]
      equation_rows <- equation_index[cbind(
        selected$i,
        rep(stage + 1L, nrow(selected))
      )]
      plus_columns <- owners + stage * coordinates + selected$j
      minus_columns <- owners + beta_block +
        stage * coordinates + selected$j
      part <- part + 1L
      triplet_parts[[part]] <- cbind(
        equation_rows,
        plus_columns,
        selected$x
      )
      part <- part + 1L
      triplet_parts[[part]] <- cbind(
        equation_rows,
        minus_columns,
        -selected$x
      )
    }
  }
  target_keep <- which(target != 0)
  if (length(target_keep) > 0L) {
    terminal_rows <- equation_index[cbind(
      target_keep,
      row_rates[target_keep] + 1L
    )]
    part <- part + 1L
    triplet_parts[[part]] <- cbind(
      terminal_rows,
      owner[target_keep],
      -target[target_keep]
    )
  }
  sum_row <- equations + 1L
  margin_rows <- sum_row + seq_len(owners)
  part <- part + 1L
  triplet_parts[[part]] <- cbind(sum_row, seq_len(owners), 1)
  part <- part + 1L
  triplet_parts[[part]] <- rbind(
    cbind(margin_rows, seq_len(owners), 1),
    cbind(margin_rows, margin_index, -1)
  )
  triplet <- do.call(rbind, triplet_parts[seq_len(part)])
  storage.mode(triplet) <- "double"
  if (nrow(triplet) > max_lp_nonzeros) {
    return(list(
      evaluated = FALSE,
      status = NA_integer_,
      certified = FALSE,
      reason = "higher_order_lift_lp_nonzero_limit"
    ))
  }

  fit <- tryCatch(
    lpSolve::lp(
      direction = "max",
      objective.in = c(rep(0, variables - 1L), 1),
      const.dir = c(
        rep("=", equations + 1L),
        rep(">=", owners)
      ),
      const.rhs = c(rep(0, equations), 1, rep(0, owners)),
      dense.const = triplet,
      timeout = as.integer(lp_timeout)
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(
      evaluated = FALSE,
      status = NA_integer_,
      certified = FALSE,
      reason = "higher_order_lift_lp_error",
      error = conditionMessage(fit)
    ))
  }
  status <- as.integer(fit$status %||% -1L)
  if (!identical(status, 0L)) {
    return(list(
      evaluated = status == 2L,
      status = status,
      certified = FALSE,
      strict_excluded = status == 2L,
      numerical_margin_open = FALSE,
      reason = if (status == 2L) {
        "higher_order_lift_hierarchy_infeasible"
      } else {
        "higher_order_lift_lp_failed"
      }
    ))
  }

  solution <- as.numeric(fit$solution)
  coefficients <- solution[seq_len(owners)]
  plus <- matrix(
    solution[owners + seq_len(beta_block)],
    nrow = coordinates,
    ncol = coefficient_blocks
  )
  minus <- matrix(
    solution[owners + beta_block + seq_len(beta_block)],
    nrow = coordinates,
    ncol = coefficient_blocks
  )
  beta_coefficients <- plus - minus
  equation_residuals <- numeric(equations)
  residual_position <- 0L
  residual_scale <- 1
  for (row in seq_len(rows)) {
    for (stage in 0:row_rates[row]) {
      residual_position <- residual_position + 1L
      value <- as.numeric(
        design[row, , drop = FALSE] %*%
          beta_coefficients[, stage + 1L]
      )
      expected <- if (stage == row_rates[row]) {
        target[row] * coefficients[owner[row]]
      } else {
        0
      }
      equation_residuals[residual_position] <- value - expected
      residual_scale <- max(residual_scale, abs(value), abs(expected))
    }
  }
  maximum_residual <- max(abs(equation_residuals))
  sum_residual <- abs(sum(coefficients) - 1)
  margin <- as.numeric(fit$objval)
  residual_limit <- tolerance * residual_scale
  certified <- is.finite(margin) && margin > tolerance &&
    min(coefficients) > tolerance &&
    maximum_residual <= residual_limit && sum_residual <= tolerance
  strict_excluded <- !certified && is.finite(margin) && isTRUE(margin == 0)
  numerical_open <- !certified && is.finite(margin) && margin > 0 &&
    margin <= tolerance

  list(
    evaluated = TRUE,
    status = status,
    certified = certified,
    strict_excluded = strict_excluded,
    numerical_margin_open = numerical_open,
    reason = if (certified) {
      "higher_order_lift_hierarchy_certified"
    } else if (strict_excluded) {
      "higher_order_lift_hierarchy_strictly_excluded"
    } else {
      "higher_order_lift_hierarchy_numerical_margin_open"
    },
    rates = rates,
    coefficients = coefficients,
    beta_coefficients = beta_coefficients,
    margin = margin,
    equation_residuals = equation_residuals,
    maximum_residual = maximum_residual,
    sum_coefficient_residual = sum_residual,
    residual_limit = residual_limit,
    maximum_rate = as.integer(maximum_rate),
    lp_nonzeros = as.integer(nrow(triplet)),
    lp_equations = as.integer(equations),
    lp_constraints = as.integer(equations + 1L + owners),
    lp_variables = as.integer(variables)
  )
}

mfrmr_jml_gpcm_higher_order_face_lifts <- function(
    face_chart_result,
    certificate_tolerance = 1e-8,
    max_stages = 12L,
    max_hierarchies = 50000L,
    max_equations = 200000L,
    max_lp_nonzeros = 5e6,
    lp_timeout = 30L) {
  empty_table <- mfrmr_jml_gpcm_higher_order_lift_empty_table()
  source <- face_chart_result
  source_matches <- is.list(source) && identical(
    as.character(source$contract_version %||% "")[1],
    mfrmr_jml_gpcm_response_image_face_contract_version()
  )
  source_record <- if (is.list(source)) source else list()
  design <- source_record$adjacent_design
  target <- tryCatch(
    suppressWarnings(as.numeric(source_record$target_contrasts)),
    error = function(e) NULL
  )
  owner <- tryCatch(
    suppressWarnings(as.integer(source_record$slope_owner_by_row)),
    error = function(e) NULL
  )
  rows <- if (inherits(design, "Matrix")) nrow(design) else 0L
  coordinates <- if (inherits(design, "Matrix")) ncol(design) else 0L
  owners <- if (is.null(owner) || length(owner) == 0L) {
    0L
  } else {
    max(owner)
  }
  controls <- suppressWarnings(as.numeric(c(
    certificate_tolerance, max_stages, max_hierarchies, max_equations,
    max_lp_nonzeros, lp_timeout
  )))

  finish <- function(
      state,
      evaluated,
      classified,
      detail,
      reason_codes,
      hierarchy_table = empty_table,
      certificates = list(),
      total_hierarchies = 0L,
      finite_image = FALSE,
      closure_member = FALSE,
      outside_closure = FALSE) {
    certified_rows <- if (nrow(hierarchy_table) > 0L) {
      hierarchy_table$HierarchyLiftCertified
    } else {
      logical(0)
    }
    higher_rows <- certified_rows &
      hierarchy_table$MaximumRate > 1L &
      !hierarchy_table$SourceFaceFirstOrderReachable
    list(
      contract_version =
        mfrmr_jml_gpcm_higher_order_lift_contract_version(),
      source_contract_version = as.character(
        source_record$contract_version %||% ""
      )[1],
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      state = state,
      evaluated = isTRUE(evaluated),
      higher_order_face_lifts_classified = isTRUE(classified),
      source_face_chart_contract_matches = isTRUE(source_matches),
      source_face_chart_classified = isTRUE(
        source_record$response_image_face_chart_classified
      ),
      zero_affine_offset_certified = isTRUE(
        source_record$zero_affine_offset_certified
      ),
      semialgebraic_curve_selection_reduction_derived = isTRUE(evaluated),
      puiseux_owner_leading_rate_reduction_derived = isTRUE(evaluated),
      distinct_owner_rates_compressible_to_ordered_integer_stages =
        isTRUE(evaluated),
      hierarchy_linear_conditions_necessary_for_closure = isTRUE(evaluated),
      hierarchy_linear_conditions_sufficient_by_explicit_path =
        isTRUE(evaluated),
      hierarchy_enumeration_required = total_hierarchies > 0,
      all_required_ordered_owner_rate_hierarchies_enumerated =
        isTRUE(evaluated) && nrow(hierarchy_table) == total_hierarchies,
      all_required_owner_rate_hierarchies_numerically_decided =
        isTRUE(classified) && nrow(hierarchy_table) == total_hierarchies,
      hierarchy_count = as.integer(nrow(hierarchy_table)),
      expected_hierarchy_count = as.numeric(total_hierarchies),
      certified_hierarchy_count = as.integer(sum(certified_rows)),
      certified_first_order_hierarchy_count = as.integer(sum(
        certified_rows & hierarchy_table$FirstOrderHierarchy
      )),
      certified_higher_order_hierarchy_count = as.integer(sum(higher_rows)),
      rows = as.integer(rows),
      additive_free_coordinates = as.integer(coordinates),
      slope_owners = as.integer(owners),
      target_contrasts = as.numeric(target %||% numeric(0)),
      adjacent_design = design,
      slope_owner_by_row = as.integer(owner %||% integer(0)),
      finite_response_image_certified = isTRUE(finite_image),
      finite_response_image_excluded = isTRUE(
        source_record$finite_response_image_excluded_by_full_support_lp
      ),
      source_closure_membership_certified = isTRUE(
        source_record$target_in_finite_response_image_closure_certified
      ),
      target_in_finite_response_image_closure_certified =
        isTRUE(closure_member),
      reachable_missing_finite_boundary_certified =
        !isTRUE(finite_image) && isTRUE(closure_member),
      target_outside_finite_response_image_closure_certified =
        isTRUE(outside_closure),
      targetwise_zero_offset_closure_membership_classified =
        isTRUE(closure_member) || isTRUE(outside_closure),
      source_first_order_open_face_resolved_by_higher_order =
        any(higher_rows),
      nonnegative_face_condition_sufficient_for_general_closure = FALSE,
      symbolic_entire_operator_closure_stratified = FALSE,
      complete_general_boundary_likelihood_envelope_constructed = FALSE,
      hierarchy_table = hierarchy_table,
      hierarchy_certificates = certificates,
      path_search_performed = FALSE,
      optimizer_trace_used = FALSE,
      competitive_boundary_certified = FALSE,
      global_supremum_value_classified = FALSE,
      finite_jmle_existence_certified = FALSE,
      finite_jmle_nonexistence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_theorem_only",
      next_gate = if (isTRUE(classified)) {
        "general_boundary_likelihood_envelope_over_rate_hierarchies"
      } else {
        "complete_numerically_decided_owner_rate_hierarchy_enumeration"
      },
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The contract classifies targetwise closure membership only for a",
        "fixed zero-offset operator when every ordered owner-rate hierarchy",
        "is workload-admitted and numerically decided. It does not produce a",
        "symbolic global stratum decomposition, construct the likelihood",
        "envelope over those strata, adjudicate a finite JMLE, or authorize",
        "MML, uncertainty, readiness, recovery, simulation, or comparison."
      )
    )
  }

  valid_controls <- length(controls) == 6L && all(is.finite(controls)) &&
    controls[1] > 0 && controls[1] < 1 &&
    all(controls[2:6] == floor(controls[2:6])) &&
    all(controls[2:6] >= 1)
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "The tolerance or hierarchy workload controls were invalid.",
      "higher_order_lift_control_invalid"
    ))
  }
  valid_source <- source_matches &&
    isTRUE(source_record$response_image_face_chart_classified) &&
    isTRUE(source_record$zero_affine_offset_certified) &&
    inherits(design, "Matrix") && rows >= 1L && coordinates >= 1L &&
    !is.null(target) && length(target) == rows && all(is.finite(target)) &&
    !is.null(owner) && length(owner) == rows && all(is.finite(owner)) &&
    owners >= 1L && all(owner >= 1L & owner <= owners)
  if (!valid_source) {
    return(finish(
      "not_evaluated_source", FALSE, FALSE,
      "A classified zero-offset P2j face chart is required.",
      "higher_order_lift_source_invalid"
    ))
  }
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    return(finish(
      "not_evaluated_dependency", FALSE, FALSE,
      "Package 'lpSolve' is required for higher-order hierarchy LPs.",
      "higher_order_lift_lp_dependency_unavailable"
    ))
  }

  if (isTRUE(source_record$finite_response_image_certified)) {
    return(finish(
      "finite_source_response_image_closure_certified", TRUE, TRUE,
      "The P2j source already supplies a finite full-support representative.",
      "higher_order_lift_finite_source",
      finite_image = TRUE,
      closure_member = TRUE
    ))
  }
  if (isTRUE(source_record$target_outside_closure_outer_bound_certified)) {
    return(finish(
      "outside_closure_certified_by_source_simplex_outer_bound", TRUE, TRUE,
      paste(
        "The complete necessary P2j simplex-face chart has no nonnegative",
        "witness, so no owner-rate hierarchy can exist."
      ),
      "higher_order_lift_source_outer_exclusion",
      outside_closure = TRUE
    ))
  }

  faces <- source_record$face_table
  compatible <- which(
    faces$ProperFace & faces$RelativeInteriorCompatible
  )
  if (length(compatible) < 1L) {
    return(finish(
      "not_evaluated_source_faces", TRUE, FALSE,
      "The source supplied no compatible proper face to refine.",
      "higher_order_lift_source_face_missing"
    ))
  }
  inactive_counts <- owners - faces$SupportSize[compatible]
  if (any(inactive_counts > controls[2])) {
    return(finish(
      "not_evaluated_stage_cap", FALSE, FALSE,
      paste(
        "The stage cap is smaller than an inactive-owner count, so the",
        "ordered hierarchy enumeration would be incomplete."
      ),
      "higher_order_lift_stage_cap_incomplete"
    ))
  }
  hierarchy_counts <- vapply(inactive_counts, function(count) {
    mfrmr_jml_gpcm_higher_order_lift_partition_count(
      count,
      max_stages = controls[2]
    )
  }, numeric(1))
  total_hierarchies <- sum(hierarchy_counts)
  if (!is.finite(total_hierarchies) || total_hierarchies > controls[3]) {
    return(finish(
      "not_evaluated_hierarchy_workload", FALSE, FALSE,
      "The complete ordered hierarchy enumeration exceeded its cap.",
      "higher_order_lift_hierarchy_workload_exceeded",
      total_hierarchies = total_hierarchies
    ))
  }

  hierarchy_rates <- list()
  hierarchy_source_faces <- character(0)
  hierarchy_source_first_order <- logical(0)
  hierarchy_position <- 0L
  for (face_position in compatible) {
    support <- as.integer(strsplit(
      faces$SupportCode[face_position],
      ",",
      fixed = TRUE
    )[[1]])
    inactive <- setdiff(seq_len(owners), support)
    for (stage_count in seq_len(length(inactive))) {
      assignments <- mfrmr_jml_gpcm_higher_order_lift_surjections(
        length(inactive),
        stage_count
      )
      for (assignment in assignments) {
        hierarchy_position <- hierarchy_position + 1L
        rates <- integer(owners)
        rates[inactive] <- assignment
        hierarchy_rates[[hierarchy_position]] <- rates
        hierarchy_source_faces[hierarchy_position] <-
          faces$FaceId[face_position]
        hierarchy_source_first_order[hierarchy_position] <-
          faces$FirstOrderBoundaryReachable[face_position]
      }
    }
  }
  if (hierarchy_position != as.integer(total_hierarchies)) {
    return(finish(
      "not_evaluated_hierarchy_generation", TRUE, FALSE,
      "The generated hierarchy count did not match the exact ordered count.",
      "higher_order_lift_hierarchy_count_mismatch",
      total_hierarchies = total_hierarchies
    ))
  }

  certificates <- vector("list", hierarchy_position)
  table_rows <- vector("list", hierarchy_position)
  hierarchy_ids <- vapply(hierarchy_rates, function(rates) {
    paste0("rates_", paste(rates, collapse = "_"))
  }, character(1))
  names(certificates) <- hierarchy_ids
  for (index in seq_len(hierarchy_position)) {
    rates <- hierarchy_rates[[index]]
    certificate <- mfrmr_jml_gpcm_higher_order_lift_lp(
      design = design,
      target = target,
      owner = owner,
      rates = rates,
      tolerance = controls[1],
      max_equations = controls[4],
      max_lp_nonzeros = controls[5],
      lp_timeout = controls[6]
    )
    certificate$hierarchy_id <- hierarchy_ids[index]
    certificate$source_face_id <- hierarchy_source_faces[index]
    certificate$source_face_first_order_reachable <-
      hierarchy_source_first_order[index]
    certificates[[index]] <- certificate
    table_rows[[index]] <- data.frame(
      HierarchyId = hierarchy_ids[index],
      ZeroStageSupport = paste(which(rates == 0L), collapse = ","),
      RateCode = paste(rates, collapse = ","),
      PositiveStages = as.integer(length(unique(rates[rates > 0L]))),
      MaximumRate = as.integer(max(rates)),
      FirstOrderHierarchy = max(rates) == 1L,
      SourceFaceFirstOrderReachable =
        hierarchy_source_first_order[index],
      SolverStatus = as.integer(certificate$status %||% NA_integer_),
      MaximumPositiveMargin = as.numeric(certificate$margin %||% NA_real_),
      MaximumAbsoluteEquationResidual =
        as.numeric(certificate$maximum_residual %||% NA_real_),
      SumCoefficientResidual =
        as.numeric(certificate$sum_coefficient_residual %||% NA_real_),
      HierarchyLiftCertified = isTRUE(certificate$certified),
      StrictHierarchyExcluded = isTRUE(certificate$strict_excluded),
      NumericalMarginOpen = isTRUE(certificate$numerical_margin_open),
      stringsAsFactors = FALSE
    )
  }
  hierarchy_table <- do.call(rbind, table_rows)
  rownames(hierarchy_table) <- NULL
  solver_complete <- all(vapply(certificates, function(certificate) {
    isTRUE(certificate$evaluated) &&
      as.integer(certificate$status %||% -1L) %in% c(0L, 2L)
  }, logical(1)))
  if (!solver_complete) {
    return(finish(
      "not_evaluated_solver", TRUE, FALSE,
      "At least one owner-rate hierarchy LP did not complete.",
      "higher_order_lift_solver_incomplete",
      hierarchy_table,
      certificates,
      total_hierarchies
    ))
  }

  certified <- hierarchy_table$HierarchyLiftCertified
  if (any(certified)) {
    return(finish(
      "missing_boundary_hierarchy_lift_certified", TRUE, TRUE,
      paste(
        "At least one exhaustive owner-rate hierarchy supplies an explicit",
        "finite product-one path to the missing response boundary."
      ),
      "higher_order_lift_missing_boundary_certified",
      hierarchy_table,
      certificates,
      total_hierarchies,
      closure_member = TRUE
    ))
  }
  if (any(hierarchy_table$NumericalMarginOpen)) {
    source_closure <- isTRUE(
      source_record$target_in_finite_response_image_closure_certified
    )
    return(finish(
      if (source_closure) {
        "source_closure_certified_hierarchy_numerical_margin_open"
      } else {
        "not_certified_numerical_margin"
      },
      TRUE, FALSE,
      if (source_closure) {
        paste(
          "The higher-order hierarchy margin is numerically open, but the",
          "source P2j first-order finite path remains a closure certificate."
        )
      } else {
        paste(
          "At least one hierarchy optimum lies strictly between zero and the",
          "certificate tolerance; closure membership remains open."
        )
      },
      "higher_order_lift_numerical_margin_open",
      hierarchy_table,
      certificates,
      total_hierarchies,
      closure_member = source_closure
    ))
  }
  all_excluded <- all(hierarchy_table$StrictHierarchyExcluded)
  if (!all_excluded) {
    return(finish(
      "not_certified_hierarchy_decision", TRUE, FALSE,
      "At least one hierarchy was neither certified nor strictly excluded.",
      "higher_order_lift_hierarchy_decision_open",
      hierarchy_table,
      certificates,
      total_hierarchies
    ))
  }
  if (isTRUE(
    source_record$target_in_finite_response_image_closure_certified
  )) {
    return(finish(
      "not_certified_source_hierarchy_conflict", TRUE, FALSE,
      paste(
        "Every hierarchy was numerically excluded although the source P2j",
        "record contains an explicit finite path; source evidence is retained."
      ),
      "higher_order_lift_source_hierarchy_conflict",
      hierarchy_table,
      certificates,
      total_hierarchies,
      closure_member = TRUE
    ))
  }
  finish(
    "outside_closure_all_rate_hierarchies_excluded", TRUE, TRUE,
    paste(
      "Every ordered positive owner-rate hierarchy on every compatible",
      "simplex face was strictly infeasible."
    ),
    "higher_order_lift_complete_outer_face_exclusion",
    hierarchy_table,
    certificates,
    total_hierarchies,
    outside_closure = TRUE
  )
}

mfrmr_jml_gpcm_higher_order_face_path_at <- function(
    result,
    hierarchy_id,
    index,
    max_index = 100L,
    max_log_rate_span = 600) {
  valid_source <- is.list(result) && identical(
    as.character(result$contract_version %||% "")[1],
    mfrmr_jml_gpcm_higher_order_lift_contract_version()
  ) && isTRUE(result$higher_order_face_lifts_classified)
  key <- as.character(hierarchy_id %||% "")[1]
  certificate <- if (valid_source && nzchar(key)) {
    result$hierarchy_certificates[[key]]
  } else {
    NULL
  }
  index_value <- tryCatch(
    suppressWarnings(as.numeric(index)[1]),
    error = function(e) NA_real_
  )
  max_value <- tryCatch(
    suppressWarnings(as.numeric(max_index)[1]),
    error = function(e) NA_real_
  )
  span_value <- tryCatch(
    suppressWarnings(as.numeric(max_log_rate_span)[1]),
    error = function(e) NA_real_
  )
  maximum_rate <- as.numeric(certificate$maximum_rate %||% NA_real_)[1]
  valid_index <- is.finite(index_value) && index_value == floor(index_value) &&
    index_value >= 0 && is.finite(max_value) && max_value == floor(max_value) &&
    max_value >= 0 && index_value <= max_value && is.finite(span_value) &&
    span_value >= 1 && is.finite(maximum_rate) &&
    index_value * maximum_rate <= span_value
  valid_hierarchy <- is.list(certificate) && isTRUE(certificate$certified)
  if (!valid_source || !valid_index || !valid_hierarchy) {
    return(list(
      contract_version =
        mfrmr_jml_gpcm_higher_order_lift_contract_version(),
      valid = FALSE,
      state = "not_evaluated_higher_order_path",
      hierarchy_id = key,
      index = as.integer(if (is.finite(index_value)) index_value else NA),
      reason_codes = "higher_order_lift_path_invalid"
    ))
  }

  rates <- certificate$rates
  coefficients <- certificate$coefficients
  beta_coefficients <- certificate$beta_coefficients
  design <- result$adjacent_design
  owner <- result$slope_owner_by_row
  target <- result$target_contrasts
  epsilon <- exp(-index_value)
  powers <- exp(-index_value * 0:certificate$maximum_rate)
  raw_beta <- as.numeric(beta_coefficients %*% powers)
  log_raw_inverse <- log(coefficients) - index_value * rates
  raw_inverse <- exp(log_raw_inverse)
  contrasts <- as.numeric(design %*% raw_beta) / raw_inverse[owner]
  raw_residual <- contrasts * raw_inverse[owner] -
    as.numeric(design %*% raw_beta)
  log_scale <- -mean(log_raw_inverse)
  log_inverse <- log_raw_inverse + log_scale
  inverse_slopes <- exp(log_inverse)
  slopes <- exp(-log_inverse)
  additive <- exp(log_scale) * raw_beta

  list(
    contract_version =
      mfrmr_jml_gpcm_higher_order_lift_contract_version(),
    valid = TRUE,
    state = "finite_higher_order_face_lift_path_evaluated",
    hierarchy_id = key,
    index = as.integer(index_value),
    epsilon = epsilon,
    rates = rates,
    target_contrasts = target,
    contrasts = contrasts,
    contrast_distance_to_target = sqrt(sum((contrasts - target)^2)),
    raw_inverse_slopes = raw_inverse,
    inverse_slopes = inverse_slopes,
    slopes = slopes,
    log_slopes = log(slopes),
    additive_coordinates = additive,
    maximum_absolute_raw_equation_residual = max(abs(raw_residual)),
    log_slope_sum_residual = sum(log(slopes)),
    inverse_slope_product_residual = prod(inverse_slopes) - 1,
    finite_parameter_point = all(is.finite(c(
      inverse_slopes, slopes, additive
    ))) && all(inverse_slopes > 0) && all(slopes > 0),
    parameter_path_declared_divergent = any(rates > 0L),
    response_contrast_limit_certified = TRUE,
    reason_codes = "higher_order_lift_finite_path_evaluated"
  )
}

mfrmr_jml_gpcm_higher_order_face_p2i_fixture <- function(contrasts) {
  source <- mfrmr_jml_gpcm_response_image_face_p2i_fixture(contrasts)
  mfrmr_jml_gpcm_higher_order_face_lifts(source)
}
