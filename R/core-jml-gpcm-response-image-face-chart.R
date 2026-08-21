# General zero-offset JML GPCM response-image simplex-face chart
# ==============================================================================
#
# For a fixed adjacent-category affine operator A, row-to-slope-owner incidence
# P, positive slopes alpha, and inverse slopes w = 1 / alpha, a response
# contrast vector z has a finite representation exactly when
#
#   diag(z) P w = A beta,  w > 0.
#
# With zero affine offset the product-one normalization of w is immaterial for
# membership: a common positive rescaling of (w, beta) restores product(w)=1.
# Normalizing instead by sum(w)=1 turns every bounded-contrast parameter
# sequence into a compact simplex sequence.  Consequently every finite closure
# point has a nonzero nonnegative face witness.  A proper face is reachable by
# an explicit first-order finite path when its inactive-owner target rows can
# be lifted by one additive direction.  The nonnegative face condition alone
# is necessary but is not sufficient for a general design.

mfrmr_jml_gpcm_response_image_face_contract_version <- function() {
  "mfrmr-jml-gpcm-response-image-face-chart-0.2.3-v1"
}

mfrmr_jml_gpcm_response_image_empty_faces <- function() {
  data.frame(
    FaceId = character(0),
    SupportCode = character(0),
    SupportSize = integer(0),
    InactiveCode = character(0),
    ProperFace = logical(0),
    SolverStatus = integer(0),
    MaximumRelativeMargin = numeric(0),
    MaximumAbsoluteEquationResidual = numeric(0),
    SumWeightResidual = numeric(0),
    RelativeInteriorCompatible = logical(0),
    StrictRelativeInteriorExcluded = logical(0),
    FirstOrderLiftEvaluated = logical(0),
    LiftSolverStatus = integer(0),
    LiftMargin = numeric(0),
    LiftMaximumAbsoluteResidual = numeric(0),
    FirstOrderBoundaryReachable = logical(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_response_image_sparse_design <- function(value) {
  tryCatch(
    methods::as(
      methods::as(Matrix::Matrix(value, sparse = TRUE), "generalMatrix"),
      "CsparseMatrix"
    ),
    error = function(e) NULL
  )
}

mfrmr_jml_gpcm_response_image_lp <- function(
    design,
    target,
    owner,
    support,
    tolerance,
    max_lp_nonzeros,
    lp_timeout) {
  rows <- nrow(design)
  coordinates <- ncol(design)
  owners <- max(owner)
  support <- as.integer(support)
  inactive <- setdiff(seq_len(owners), support)
  margin_index <- owners + 2L * coordinates + 1L
  variables <- margin_index
  sum_row <- rows + 1L
  margin_rows <- sum_row + seq_along(support)
  inactive_rows <- sum_row + length(support) + seq_along(inactive)
  constraints <- rows + 1L + length(support) + length(inactive)

  target_keep <- which(target != 0)
  triplet <- if (length(target_keep) > 0L) {
    cbind(target_keep, owner[target_keep], target[target_keep])
  } else {
    matrix(numeric(0), nrow = 0L, ncol = 3L)
  }
  design_summary <- Matrix::summary(design)
  if (nrow(design_summary) > 0L) {
    triplet <- rbind(
      triplet,
      cbind(
        design_summary$i,
        owners + design_summary$j,
        -design_summary$x
      ),
      cbind(
        design_summary$i,
        owners + coordinates + design_summary$j,
        design_summary$x
      )
    )
  }
  triplet <- rbind(
    triplet,
    cbind(sum_row, seq_len(owners), 1),
    cbind(margin_rows, support, 1),
    cbind(margin_rows, margin_index, -1)
  )
  if (length(inactive) > 0L) {
    triplet <- rbind(
      triplet,
      cbind(inactive_rows, inactive, 1)
    )
  }
  storage.mode(triplet) <- "double"
  if (nrow(triplet) > max_lp_nonzeros) {
    return(list(
      evaluated = FALSE,
      status = NA_integer_,
      reason = "response_image_face_lp_nonzero_limit"
    ))
  }

  fit <- tryCatch(
    lpSolve::lp(
      direction = "max",
      objective.in = c(rep(0, variables - 1L), 1),
      const.dir = c(
        rep("=", rows + 1L),
        rep(">=", length(support)),
        rep("=", length(inactive))
      ),
      const.rhs = c(rep(0, rows), 1, rep(0, owners)),
      dense.const = triplet,
      timeout = as.integer(lp_timeout)
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(
      evaluated = FALSE,
      status = NA_integer_,
      reason = "response_image_face_lp_error",
      error = conditionMessage(fit)
    ))
  }
  status <- as.integer(fit$status %||% -1L)
  if (!identical(status, 0L)) {
    return(list(
      evaluated = status == 2L,
      status = status,
      reason = if (status == 2L) {
        "response_image_face_lp_infeasible"
      } else {
        "response_image_face_lp_failed"
      }
    ))
  }

  solution <- as.numeric(fit$solution)
  weights <- solution[seq_len(owners)]
  beta <- solution[owners + seq_len(coordinates)] -
    solution[owners + coordinates + seq_len(coordinates)]
  equation_residual <- target * weights[owner] -
    as.numeric(design %*% beta)
  scale <- max(
    1,
    max(abs(target * weights[owner])),
    max(abs(as.numeric(design %*% beta)))
  )
  residual_limit <- tolerance * scale
  maximum_residual <- max(abs(equation_residual))
  sum_residual <- abs(sum(weights) - 1)
  inactive_maximum <- if (length(inactive) > 0L) {
    max(abs(weights[inactive]))
  } else {
    0
  }
  margin <- as.numeric(fit$objval)
  compatible <- is.finite(margin) && margin > tolerance &&
    maximum_residual <= residual_limit && sum_residual <= tolerance &&
    inactive_maximum <= tolerance &&
    min(weights[support]) > tolerance

  list(
    evaluated = TRUE,
    status = status,
    reason = if (compatible) {
      "response_image_face_relative_interior_certified"
    } else {
      "response_image_face_relative_interior_excluded"
    },
    support = support,
    inactive = inactive,
    weights = weights,
    beta = beta,
    margin = margin,
    equation_residual = equation_residual,
    maximum_residual = maximum_residual,
    sum_weight_residual = sum_residual,
    residual_limit = residual_limit,
    relative_interior_compatible = compatible,
    strict_relative_interior_excluded = !compatible &&
      is.finite(margin) && isTRUE(margin == 0),
    lp_nonzeros = as.integer(nrow(triplet)),
    lp_constraints = as.integer(constraints),
    lp_variables = as.integer(variables)
  )
}

mfrmr_jml_gpcm_response_image_lift_lp <- function(
    design,
    target,
    owner,
    face,
    tolerance,
    max_lp_nonzeros,
    lp_timeout) {
  inactive <- face$inactive
  inactive_rows <- which(owner %in% inactive)
  coordinates <- ncol(design)
  inactive_count <- length(inactive)
  owner_position <- match(owner[inactive_rows], inactive)
  margin_index <- inactive_count + 2L * coordinates + 1L
  variables <- margin_index
  equations <- length(inactive_rows)
  sum_row <- equations + 1L
  margin_rows <- sum_row + seq_len(inactive_count)

  target_keep <- which(target[inactive_rows] != 0)
  triplet <- if (length(target_keep) > 0L) {
    cbind(
      target_keep,
      owner_position[target_keep],
      -target[inactive_rows[target_keep]]
    )
  } else {
    matrix(numeric(0), nrow = 0L, ncol = 3L)
  }
  selected_design <- design[inactive_rows, , drop = FALSE]
  design_summary <- Matrix::summary(selected_design)
  if (nrow(design_summary) > 0L) {
    triplet <- rbind(
      triplet,
      cbind(
        design_summary$i,
        inactive_count + design_summary$j,
        design_summary$x
      ),
      cbind(
        design_summary$i,
        inactive_count + coordinates + design_summary$j,
        -design_summary$x
      )
    )
  }
  triplet <- rbind(
    triplet,
    cbind(sum_row, seq_len(inactive_count), 1),
    cbind(margin_rows, seq_len(inactive_count), 1),
    cbind(margin_rows, margin_index, -1)
  )
  storage.mode(triplet) <- "double"
  if (nrow(triplet) > max_lp_nonzeros) {
    return(list(
      evaluated = FALSE,
      status = NA_integer_,
      certified = FALSE,
      reason = "response_image_lift_lp_nonzero_limit"
    ))
  }

  fit <- tryCatch(
    lpSolve::lp(
      direction = "max",
      objective.in = c(rep(0, variables - 1L), 1),
      const.dir = c(
        rep("=", equations + 1L),
        rep(">=", inactive_count)
      ),
      const.rhs = c(rep(0, equations), 1, rep(0, inactive_count)),
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
      reason = "response_image_lift_lp_error",
      error = conditionMessage(fit)
    ))
  }
  status <- as.integer(fit$status %||% -1L)
  if (!identical(status, 0L)) {
    return(list(
      evaluated = status == 2L,
      status = status,
      certified = FALSE,
      reason = if (status == 2L) {
        "response_image_first_order_lift_infeasible"
      } else {
        "response_image_first_order_lift_failed"
      }
    ))
  }

  solution <- as.numeric(fit$solution)
  inactive_weights <- solution[seq_len(inactive_count)]
  gamma <- solution[inactive_count + seq_len(coordinates)] -
    solution[inactive_count + coordinates + seq_len(coordinates)]
  residual <- as.numeric(selected_design %*% gamma) -
    target[inactive_rows] * inactive_weights[owner_position]
  scale <- max(
    1,
    max(abs(as.numeric(selected_design %*% gamma))),
    max(abs(target[inactive_rows] * inactive_weights[owner_position]))
  )
  residual_limit <- tolerance * scale
  maximum_residual <- max(abs(residual))
  margin <- as.numeric(fit$objval)
  certified <- is.finite(margin) && margin > tolerance &&
    maximum_residual <= residual_limit &&
    abs(sum(inactive_weights) - 1) <= tolerance &&
    min(inactive_weights) > tolerance

  list(
    evaluated = TRUE,
    status = status,
    certified = certified,
    reason = if (certified) {
      "response_image_first_order_lift_certified"
    } else {
      "response_image_first_order_lift_not_certified"
    },
    inactive_owners = inactive,
    inactive_rows = inactive_rows,
    inactive_weights = inactive_weights,
    gamma = gamma,
    margin = margin,
    residual = residual,
    maximum_residual = maximum_residual,
    residual_limit = residual_limit,
    lp_nonzeros = as.integer(nrow(triplet)),
    lp_constraints = as.integer(equations + 1L + inactive_count),
    lp_variables = as.integer(variables)
  )
}

mfrmr_jml_gpcm_response_image_face_chart <- function(
    adjacent_design,
    slope_owner_by_row,
    target_contrasts,
    additive_offset = NULL,
    certificate_tolerance = 1e-8,
    max_rows = 20000L,
    max_coordinates = 2000L,
    max_owners = 12L,
    max_faces = 4095L,
    max_lp_nonzeros = 5e6,
    lp_timeout = 30L) {
  empty_faces <- mfrmr_jml_gpcm_response_image_empty_faces()
  design <- mfrmr_jml_gpcm_response_image_sparse_design(adjacent_design)
  target <- tryCatch(
    suppressWarnings(as.numeric(target_contrasts)),
    error = function(e) NULL
  )
  owner <- tryCatch(
    suppressWarnings(as.numeric(slope_owner_by_row)),
    error = function(e) NULL
  )
  offset <- if (is.null(additive_offset) && !is.null(design)) {
    rep(0, nrow(design))
  } else {
    tryCatch(
      suppressWarnings(as.numeric(additive_offset)),
      error = function(e) NULL
    )
  }
  controls <- suppressWarnings(as.numeric(c(
    certificate_tolerance, max_rows, max_coordinates, max_owners,
    max_faces, max_lp_nonzeros, lp_timeout
  )))
  rows <- if (is.null(design)) 0L else nrow(design)
  coordinates <- if (is.null(design)) 0L else ncol(design)
  owner_maximum <- if (is.null(owner) || length(owner) == 0L ||
                       any(!is.finite(owner))) {
    NA_real_
  } else {
    max(owner)
  }
  owner_count <- if (!is.finite(owner_maximum) || owner_maximum < 1 ||
                     owner_maximum > .Machine$integer.max) {
    0L
  } else {
    as.integer(owner_maximum)
  }
  owner_record <- if (owner_count >= 1L && !is.null(owner) &&
                      all(is.finite(owner)) && all(owner == floor(owner)) &&
                      all(owner >= 1) && all(owner <= .Machine$integer.max)) {
    as.integer(owner)
  } else {
    integer(0)
  }

  finish <- function(
      state,
      evaluated,
      classified,
      detail,
      reason_codes,
      faces = empty_faces,
      certificates = list(),
      finite_representative = list(),
      finite_image = FALSE,
      finite_excluded = FALSE,
      nonnegative_face = FALSE,
      reachable_boundary = FALSE,
      outside_outer = FALSE) {
    compatible_proper <- if (nrow(faces) > 0L) {
      sum(faces$ProperFace & faces$RelativeInteriorCompatible)
    } else {
      0L
    }
    reachable_faces <- if (nrow(faces) > 0L) {
      sum(faces$FirstOrderBoundaryReachable)
    } else {
      0L
    }
    list(
      contract_version =
        mfrmr_jml_gpcm_response_image_face_contract_version(),
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      response_contrast_identity =
        "complete retained adjacent-category logit contrast vector",
      state = state,
      evaluated = isTRUE(evaluated),
      response_image_face_chart_classified = isTRUE(classified),
      zero_affine_offset_required = TRUE,
      zero_affine_offset_certified = !is.null(offset) &&
        length(offset) == rows && all(offset == 0),
      finite_image_positive_inverse_slope_equivalence_derived =
        isTRUE(classified),
      inverse_slope_simplex_compactification_derived = isTRUE(classified),
      every_finite_closure_point_has_nonnegative_face_witness =
        isTRUE(classified),
      proper_simplex_face_necessary_for_missing_boundary =
        isTRUE(classified),
      first_order_face_lift_sufficiency_derived = isTRUE(classified),
      nonnegative_face_condition_sufficient_for_general_closure = FALSE,
      all_nonempty_simplex_faces_enumerated = isTRUE(classified),
      face_count = as.integer(nrow(faces)),
      compatible_proper_face_count = as.integer(compatible_proper),
      first_order_reachable_face_count = as.integer(reachable_faces),
      rows = as.integer(rows),
      additive_free_coordinates = as.integer(coordinates),
      slope_owners = as.integer(owner_count),
      target_contrasts = as.numeric(target %||% numeric(0)),
      adjacent_design = design,
      slope_owner_by_row = owner_record,
      additive_offset = as.numeric(offset %||% numeric(0)),
      finite_response_image_certified = isTRUE(finite_image),
      finite_response_image_excluded_by_full_support_lp =
        isTRUE(finite_excluded),
      nonnegative_simplex_face_witness_exists = isTRUE(nonnegative_face),
      target_in_finite_response_image_closure_certified =
        isTRUE(finite_image) || isTRUE(reachable_boundary),
      reachable_missing_finite_boundary_certified =
        isTRUE(reachable_boundary),
      target_outside_closure_outer_bound_certified = isTRUE(outside_outer),
      response_image_closure_membership_open = isTRUE(nonnegative_face) &&
        !isTRUE(finite_image) && !isTRUE(reachable_boundary),
      complete_general_response_image_closure_stratified = FALSE,
      complete_general_boundary_likelihood_envelope_constructed = FALSE,
      finite_representative = finite_representative,
      face_table = faces,
      face_certificates = certificates,
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
        "higher_order_face_lifts_and_general_closure_completion"
      } else {
        "valid_zero_offset_response_image_face_chart"
      },
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The chart gives an exact algebraic finite-image equivalence, a",
        "complete necessary simplex-face outer stratification, and explicit",
        "sufficient first-order boundary lifts only for a fixed zero-offset",
        "operator. A nonnegative face without a certified lift is not",
        "promoted to closure membership. Higher-order lifts, nonzero affine",
        "offsets, a complete general likelihood envelope, finite-JMLE",
        "existence, MML, uncertainty, readiness, recovery, simulation, and",
        "external comparison remain unclassified."
      )
    )
  }

  valid_controls <- length(controls) == 7L && all(is.finite(controls)) &&
    controls[1] > 0 && controls[1] < 1 &&
    all(controls[2:7] == floor(controls[2:7])) &&
    all(controls[2:7] >= 1)
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "The tolerance or workload controls were invalid.",
      "response_image_face_control_invalid"
    ))
  }
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    return(finish(
      "not_evaluated_dependency", FALSE, FALSE,
      "Package 'lpSolve' is required for the simplex-face chart.",
      "response_image_face_lp_dependency_unavailable"
    ))
  }
  valid_design <- !is.null(design) && rows >= 1L && coordinates >= 1L &&
    rows <= controls[2] && coordinates <= controls[3] &&
    all(is.finite(design@x))
  valid_target <- !is.null(target) && length(target) == rows &&
    all(is.finite(target))
  valid_owner <- !is.null(owner) && length(owner) == rows &&
    all(is.finite(owner)) && all(owner == floor(owner)) && all(owner >= 1) &&
    owner_count >= 1L && owner_count <= controls[4] &&
    identical(sort(unique(as.integer(owner))), seq_len(owner_count))
  valid_offset <- !is.null(offset) && length(offset) == rows &&
    all(is.finite(offset))
  if (!valid_design || !valid_target || !valid_owner || !valid_offset) {
    return(finish(
      "not_evaluated_input", FALSE, FALSE,
      "The fixed design, target, owner map, or affine offset was invalid.",
      "response_image_face_input_invalid"
    ))
  }
  if (any(offset != 0)) {
    return(finish(
      "not_evaluated_nonzero_offset", FALSE, FALSE,
      paste(
        "The product normalization cannot be removed from the affine",
        "membership system when the offset is nonzero."
      ),
      "response_image_face_nonzero_offset_unsupported"
    ))
  }
  face_count <- 2^owner_count - 1
  if (!is.finite(face_count) || face_count > controls[5]) {
    return(finish(
      "not_evaluated_face_workload", FALSE, FALSE,
      "The complete nonempty simplex-face enumeration exceeded its cap.",
      "response_image_face_workload_exceeded"
    ))
  }

  face_codes <- seq_len(as.integer(face_count))
  support_list <- lapply(face_codes, function(code) {
    which(as.logical(intToBits(code)[seq_len(owner_count)]))
  })
  ordering <- order(
    -vapply(support_list, length, integer(1)),
    face_codes
  )
  support_list <- support_list[ordering]
  certificates <- vector("list", length(support_list))
  face_rows <- vector("list", length(support_list))
  names(certificates) <- vapply(support_list, function(support) {
    paste0("support_", paste(support, collapse = "_"))
  }, character(1))

  for (index in seq_along(support_list)) {
    support <- support_list[[index]]
    face_id <- names(certificates)[index]
    face <- mfrmr_jml_gpcm_response_image_lp(
      design = design,
      target = target,
      owner = as.integer(owner),
      support = support,
      tolerance = controls[1],
      max_lp_nonzeros = controls[6],
      lp_timeout = controls[7]
    )
    face$face_id <- face_id
    face$lift <- list(
      evaluated = FALSE,
      status = NA_integer_,
      certified = FALSE,
      reason = "response_image_first_order_lift_not_applicable"
    )
    if (isTRUE(face$relative_interior_compatible) &&
        length(face$inactive) > 0L) {
      face$lift <- mfrmr_jml_gpcm_response_image_lift_lp(
        design = design,
        target = target,
        owner = as.integer(owner),
        face = face,
        tolerance = controls[1],
        max_lp_nonzeros = controls[6],
        lp_timeout = controls[7]
      )
    }
    certificates[[index]] <- face
    face_rows[[index]] <- data.frame(
      FaceId = face_id,
      SupportCode = paste(support, collapse = ","),
      SupportSize = as.integer(length(support)),
      InactiveCode = paste(setdiff(seq_len(owner_count), support), collapse = ","),
      ProperFace = length(support) < owner_count,
      SolverStatus = as.integer(face$status %||% NA_integer_),
      MaximumRelativeMargin = as.numeric(face$margin %||% NA_real_),
      MaximumAbsoluteEquationResidual =
        as.numeric(face$maximum_residual %||% NA_real_),
      SumWeightResidual = as.numeric(face$sum_weight_residual %||% NA_real_),
      RelativeInteriorCompatible = isTRUE(face$relative_interior_compatible),
      StrictRelativeInteriorExcluded = isTRUE(
        face$strict_relative_interior_excluded
      ) || identical(as.integer(face$status %||% -1L), 2L),
      FirstOrderLiftEvaluated = isTRUE(face$lift$evaluated),
      LiftSolverStatus = as.integer(face$lift$status %||% NA_integer_),
      LiftMargin = as.numeric(face$lift$margin %||% NA_real_),
      LiftMaximumAbsoluteResidual =
        as.numeric(face$lift$maximum_residual %||% NA_real_),
      FirstOrderBoundaryReachable = isTRUE(face$lift$certified),
      stringsAsFactors = FALSE
    )
  }
  faces <- do.call(rbind, face_rows)
  rownames(faces) <- NULL
  solver_complete <- all(vapply(certificates, function(face) {
    isTRUE(face$evaluated) &&
      as.integer(face$status %||% -1L) %in% c(0L, 2L)
  }, logical(1)))
  if (!solver_complete) {
    return(finish(
      "not_evaluated_solver", TRUE, FALSE,
      "At least one simplex-face linear program did not complete.",
      "response_image_face_solver_incomplete",
      faces = faces,
      certificates = certificates
    ))
  }
  face_decided <- all(
    faces$RelativeInteriorCompatible |
      faces$StrictRelativeInteriorExcluded
  )
  if (!face_decided) {
    return(finish(
      "not_certified_numerical_margin", TRUE, FALSE,
      paste(
        "At least one face optimum lay strictly between zero and the",
        "certificate tolerance; membership was left unclassified."
      ),
      "response_image_face_numerical_margin_open",
      faces = faces,
      certificates = certificates
    ))
  }

  full_id <- paste0("support_", paste(seq_len(owner_count), collapse = "_"))
  full_index <- match(full_id, faces$FaceId)
  full_face <- certificates[[full_index]]
  finite_image <- isTRUE(full_face$relative_interior_compatible)
  finite_excluded <- !finite_image && (
    identical(as.integer(full_face$status %||% -1L), 2L) ||
      isTRUE(full_face$strict_relative_interior_excluded)
  )
  nonnegative_face <- any(faces$RelativeInteriorCompatible)
  reachable_boundary <- finite_excluded &&
    any(faces$ProperFace & faces$FirstOrderBoundaryReachable)
  outside_outer <- finite_excluded && !nonnegative_face

  finite_representative <- list()
  if (finite_image) {
    raw_weights <- full_face$weights
    log_scale <- -mean(log(raw_weights))
    normalized_weights <- exp(log(raw_weights) + log_scale)
    normalized_beta <- exp(log_scale) * full_face$beta
    reconstructed <- as.numeric(design %*% normalized_beta) /
      normalized_weights[as.integer(owner)]
    finite_representative <- list(
      inverse_slopes = normalized_weights,
      slopes = 1 / normalized_weights,
      log_slopes = -log(normalized_weights),
      additive_coordinates = normalized_beta,
      reconstructed_contrasts = reconstructed,
      maximum_absolute_reconstruction_residual =
        max(abs(reconstructed - target)),
      log_slope_sum_residual = sum(-log(normalized_weights)),
      inverse_slope_product_residual = prod(normalized_weights) - 1
    )
  }

  state <- if (finite_image) {
    "finite_response_image_representative_certified"
  } else if (reachable_boundary) {
    "reachable_missing_boundary_face_certified"
  } else if (nonnegative_face) {
    "outer_boundary_face_candidate_lift_open"
  } else {
    "outside_response_image_closure_outer_bound_certified"
  }
  detail <- if (finite_image) {
    paste(
      "A strict-positive inverse-slope simplex solution and finite",
      "product-one representative were reconstructed."
    )
  } else if (reachable_boundary) {
    paste(
      "Full-support membership was excluded and at least one proper face",
      "has an explicit first-order finite-parameter lift."
    )
  } else if (nonnegative_face) {
    paste(
      "A necessary nonnegative face witness exists, but no first-order",
      "lift was certified; general closure membership remains open."
    )
  } else {
    paste(
      "The exhaustive simplex-face chart has no nonnegative witness, which",
      "excludes the target from the finite response-image closure."
    )
  }
  finish(
    state, TRUE, TRUE, detail,
    paste0("response_image_face_", state),
    faces = faces,
    certificates = certificates,
    finite_representative = finite_representative,
    finite_image = finite_image,
    finite_excluded = finite_excluded,
    nonnegative_face = nonnegative_face,
    reachable_boundary = reachable_boundary,
    outside_outer = outside_outer
  )
}

mfrmr_jml_gpcm_response_image_face_path_at <- function(
    result,
    face_id,
    index,
    max_index = 300L) {
  valid_source <- is.list(result) && identical(
    as.character(result$contract_version %||% "")[1],
    mfrmr_jml_gpcm_response_image_face_contract_version()
  ) && isTRUE(result$response_image_face_chart_classified)
  index_value <- tryCatch(
    suppressWarnings(as.numeric(index)[1]),
    error = function(e) NA_real_
  )
  max_value <- tryCatch(
    suppressWarnings(as.numeric(max_index)[1]),
    error = function(e) NA_real_
  )
  valid_index <- length(index_value) == 1L && is.finite(index_value) &&
    index_value == floor(index_value) && index_value >= 0 &&
    length(max_value) == 1L && is.finite(max_value) &&
    max_value == floor(max_value) && max_value >= 0 && index_value <= max_value
  key <- as.character(face_id %||% "")[1]
  face <- if (valid_source && nzchar(key)) {
    result$face_certificates[[key]]
  } else {
    NULL
  }
  valid_face <- is.list(face) && isTRUE(face$relative_interior_compatible) &&
    length(face$inactive) > 0L && isTRUE(face$lift$certified)
  if (!valid_source || !valid_index || !valid_face) {
    return(list(
      contract_version =
        mfrmr_jml_gpcm_response_image_face_contract_version(),
      valid = FALSE,
      state = "not_evaluated_face_path",
      index = as.integer(if (is.finite(index_value)) index_value else NA),
      face_id = key,
      reason_codes = "response_image_face_path_invalid"
    ))
  }

  design <- result$adjacent_design
  target <- result$target_contrasts
  owner <- result$slope_owner_by_row
  epsilon <- exp(-index_value)
  raw_weights <- face$weights
  raw_weights[face$inactive] <-
    epsilon * face$lift$inactive_weights
  raw_beta <- face$beta + epsilon * face$lift$gamma
  path_contrasts <- as.numeric(design %*% raw_beta) /
    raw_weights[owner]
  raw_equation_residual <- path_contrasts * raw_weights[owner] -
    as.numeric(design %*% raw_beta)
  log_scale <- -mean(log(raw_weights))
  normalized_log_inverse <- log(raw_weights) + log_scale
  normalized_inverse <- exp(normalized_log_inverse)
  slopes <- exp(-normalized_log_inverse)
  normalized_beta <- exp(log_scale) * raw_beta

  list(
    contract_version =
      mfrmr_jml_gpcm_response_image_face_contract_version(),
    valid = TRUE,
    state = "finite_response_image_face_lift_path_evaluated",
    face_id = key,
    index = as.integer(index_value),
    epsilon = epsilon,
    target_contrasts = target,
    contrasts = path_contrasts,
    contrast_distance_to_target = sqrt(sum((path_contrasts - target)^2)),
    raw_inverse_slopes = raw_weights,
    inverse_slopes = normalized_inverse,
    slopes = slopes,
    log_slopes = log(slopes),
    additive_coordinates = normalized_beta,
    maximum_absolute_raw_equation_residual =
      max(abs(raw_equation_residual)),
    log_slope_sum_residual = sum(log(slopes)),
    inverse_slope_product_residual = prod(normalized_inverse) - 1,
    finite_parameter_point = all(is.finite(c(
      normalized_inverse, slopes, normalized_beta
    ))) && all(normalized_inverse > 0) && all(slopes > 0),
    parameter_path_declared_divergent = TRUE,
    response_contrast_limit_certified = TRUE,
    reason_codes = "response_image_first_order_face_lift_path_evaluated"
  )
}

mfrmr_jml_gpcm_response_image_face_p2i_fixture <- function(contrasts) {
  mfrmr_jml_gpcm_response_image_face_chart(
    adjacent_design = rbind(
      c(1, 0, -1),
      c(0, 1, -1),
      c(1, 0, 1),
      c(0, 1, 1)
    ),
    slope_owner_by_row = c(1, 1, 2, 2),
    target_contrasts = contrasts
  )
}
