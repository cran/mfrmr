# Exact response-image closure and likelihood envelope for a minimal JML GPCM
# ==============================================================================
#
# This file completes the response geometry of the binary two-Person,
# two-slope-owner operator used by P2h.  With alpha_1 * alpha_2 = 1 and
# z = (alpha_1(theta_1-delta), alpha_1(theta_2-delta),
#      alpha_2(theta_1+delta), alpha_2(theta_2+delta)), write
# d_1 = z_1-z_2 and d_2 = z_3-z_4.  Every finite point satisfies
# d_1 = alpha_1^2 d_2.  Hence the finite image, its closure, and its missing
# axes all have exact sign descriptions.
#
# Strictly positive success and failure mass in every response cell makes the
# binary log likelihood coercive in z.  The complete boundary likelihood
# envelope then reduces to two exact pooled-logit axis candidates.

mfrmr_jml_gpcm_binary_closure_contract_version <- function() {
  "mfrmr-jml-gpcm-binary-closure-envelope-0.2.3-v1"
}

mfrmr_jml_gpcm_binary_closure_empty_strata <- function() {
  data.frame(
    Stratum = character(0),
    DifferenceCondition = character(0),
    Dimension = integer(0),
    FiniteResponseImage = logical(0),
    FiniteClosure = logical(0),
    BoundedResponseEscapeLimit = logical(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_binary_closure_strata <- function() {
  data.frame(
    Stratum = c(
      "finite_positive_difference_interior",
      "finite_negative_difference_interior",
      "finite_equal_difference_intersection",
      "missing_owner_1_dominant_axis",
      "missing_owner_2_dominant_axis"
    ),
    DifferenceCondition = c(
      "d1>0,d2>0", "d1<0,d2<0", "d1=0,d2=0",
      "d1!=0,d2=0", "d1=0,d2!=0"
    ),
    Dimension = c(4L, 4L, 2L, 3L, 3L),
    FiniteResponseImage = c(TRUE, TRUE, TRUE, FALSE, FALSE),
    FiniteClosure = rep(TRUE, 5L),
    BoundedResponseEscapeLimit = c(FALSE, FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_binary_closure_empty_axes <- function() {
  data.frame(
    Axis = character(0),
    DivergentSlopeOwner = integer(0),
    Constraint = character(0),
    Z1 = numeric(0),
    Z2 = numeric(0),
    Z3 = numeric(0),
    Z4 = numeric(0),
    DifferenceOwner1 = numeric(0),
    DifferenceOwner2 = numeric(0),
    LogLikelihood = numeric(0),
    FiniteResponseImage = logical(0),
    MissingFiniteBoundary = logical(0),
    BoundaryEnvelopeWinner = logical(0),
    GlobalClosureMaximizer = logical(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_binary_closure_image_state <- function(contrasts) {
  value <- tryCatch(
    suppressWarnings(as.numeric(contrasts)),
    error = function(e) NULL
  )
  valid <- !is.null(value) && length(value) == 4L &&
    all(is.finite(value))
  if (!valid) {
    return(list(
      valid = FALSE,
      state = "invalid_contrast_vector",
      differences = c(NA_real_, NA_real_),
      finite_response_image = FALSE,
      finite_response_image_closure = FALSE,
      missing_finite_boundary = FALSE,
      bounded_response_escape_limit = FALSE
    ))
  }
  differences <- c(value[1] - value[2], value[3] - value[4])
  both_zero <- all(differences == 0)
  same_positive <- all(differences > 0)
  same_negative <- all(differences < 0)
  opposite <- (differences[1] > 0 && differences[2] < 0) ||
    (differences[1] < 0 && differences[2] > 0)
  finite_image <- same_positive || same_negative || both_zero
  finite_closure <- !opposite
  missing_boundary <- finite_closure && !finite_image
  state <- if (same_positive) {
    "finite_positive_difference_interior"
  } else if (same_negative) {
    "finite_negative_difference_interior"
  } else if (both_zero) {
    "finite_equal_difference_intersection"
  } else if (opposite) {
    "outside_finite_response_image_closure"
  } else if (differences[1] != 0) {
    "missing_owner_1_dominant_axis"
  } else if (differences[2] != 0) {
    "missing_owner_2_dominant_axis"
  } else {
    "outside_finite_response_image_closure"
  }
  list(
    valid = TRUE,
    state = state,
    differences = differences,
    finite_response_image = finite_image,
    finite_response_image_closure = finite_closure,
    missing_finite_boundary = missing_boundary,
    bounded_response_escape_limit = both_zero || missing_boundary
  )
}

mfrmr_jml_gpcm_binary_closure_softplus <- function(value) {
  positive <- value > 0
  out <- numeric(length(value))
  out[positive] <-
    value[positive] + log1p(exp(-value[positive]))
  out[!positive] <- log1p(exp(value[!positive]))
  out
}

mfrmr_jml_gpcm_binary_closure_loglik <- function(
    contrasts,
    success_mass,
    failure_mass) {
  value <- as.numeric(contrasts)
  success <- as.numeric(success_mass)
  failure <- as.numeric(failure_mass)
  if (length(value) != 4L || length(success) != 4L ||
      length(failure) != 4L || any(!is.finite(value)) ||
      any(!is.finite(success)) || any(!is.finite(failure))) {
    return(NA_real_)
  }
  sum(
    success * value -
      (success + failure) *
        mfrmr_jml_gpcm_binary_closure_softplus(value)
  )
}

mfrmr_jml_gpcm_binary_closure_representative <- function(contrasts) {
  state <- mfrmr_jml_gpcm_binary_closure_image_state(contrasts)
  if (!isTRUE(state$valid) || !isTRUE(state$finite_response_image)) {
    return(list(
      available = FALSE,
      alpha = numeric(0),
      log_alpha = numeric(0),
      theta = numeric(0),
      delta = NA_real_,
      reconstructed_contrasts = numeric(0),
      maximum_absolute_reconstruction_residual = NA_real_
    ))
  }
  value <- as.numeric(contrasts)
  differences <- state$differences
  alpha_1 <- if (all(differences == 0)) {
    1
  } else {
    exp(0.5 * (log(abs(differences[1])) - log(abs(differences[2]))))
  }
  alpha <- c(alpha_1, 1 / alpha_1)
  theta <- c(
    (value[1] / alpha_1 + alpha_1 * value[3]) / 2,
    (value[2] / alpha_1 + alpha_1 * value[4]) / 2
  )
  delta <- (alpha_1 * value[3] - value[1] / alpha_1) / 2
  reconstructed <- c(
    alpha[1] * (theta[1] - delta),
    alpha[1] * (theta[2] - delta),
    alpha[2] * (theta[1] + delta),
    alpha[2] * (theta[2] + delta)
  )
  residual <- reconstructed - value
  list(
    available = all(is.finite(c(alpha, theta, delta, reconstructed))),
    alpha = alpha,
    log_alpha = log(alpha),
    theta = theta,
    delta = delta,
    reconstructed_contrasts = reconstructed,
    maximum_absolute_reconstruction_residual = max(abs(residual))
  )
}

mfrmr_jml_gpcm_binary_closure_envelope <- function(
    success_mass,
    failure_mass,
    source_response_quotient = NULL,
    max_total_mass = 1e9) {
  empty_axes <- mfrmr_jml_gpcm_binary_closure_empty_axes()
  empty_strata <- mfrmr_jml_gpcm_binary_closure_empty_strata()
  source <- if (is.null(source_response_quotient)) {
    mfrmr_jml_gpcm_response_quotient_binary_counterexample()
  } else {
    source_response_quotient
  }
  source_record <- if (is.list(source)) source else list()
  source_matches <- identical(
    as.character(source_record$contract_version %||% "")[1],
    mfrmr_jml_gpcm_response_quotient_contract_version()
  ) && isTRUE(
    source_record$response_quotient_closure_obstruction_certified
  )

  finish <- function(state, evaluated, classified, detail, reason_codes,
                     success = numeric(0), failure = numeric(0),
                     independent = numeric(0), independent_state = NULL,
                     strata = empty_strata, axes = empty_axes,
                     global_supremum = NA_real_,
                     boundary_supremum = NA_real_,
                     global_candidates = matrix(
                       numeric(0), nrow = 0L, ncol = 4L
                     ),
                     finite_representative = list()) {
    finite_exists <- isTRUE(classified) && isTRUE(
      independent_state$finite_response_image
    )
    finite_nonexists <- isTRUE(classified) && !finite_exists
    boundary_gap <- if (
      isTRUE(classified) && is.finite(global_supremum) &&
        is.finite(boundary_supremum)
    ) {
      global_supremum - boundary_supremum
    } else {
      NA_real_
    }
    list(
      contract_version =
        mfrmr_jml_gpcm_binary_closure_contract_version(),
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      fixture_identity =
        "binary_two_person_two_slope_owner_centered_location_gpcm",
      state = state,
      evaluated = isTRUE(evaluated),
      binary_closure_envelope_classified = isTRUE(classified),
      source_response_quotient_contract_matches = isTRUE(source_matches),
      strictly_positive_success_mass_per_cell = isTRUE(classified),
      strictly_positive_failure_mass_per_cell = isTRUE(classified),
      ordinary_single_cell_separation_absent = isTRUE(classified),
      integer_expanded_row_multiplicities = isTRUE(classified) &&
        all(c(success, failure) == floor(c(success, failure))),
      response_contrast_dimension = 4L,
      slope_owners = 2L,
      finite_response_image_condition = paste(
        "(d1>0 and d2>0) or (d1<0 and d2<0) or",
        "(d1=0 and d2=0)"
      ),
      finite_response_image_closure_condition =
        "not ((d1>0 and d2<0) or (d1<0 and d2>0))",
      finite_response_image_exactly_characterized = isTRUE(classified),
      complete_finite_response_image_closure_stratified = isTRUE(classified),
      complete_bounded_response_escape_limit_set_classified =
        isTRUE(classified),
      bounded_response_escape_limit_condition = "d1=0 or d2=0",
      unbounded_response_contrast_likelihood_limit_negative_infinity =
        isTRUE(classified),
      complete_fixture_parameter_sequence_limsup_envelope_constructed =
        isTRUE(classified),
      general_gpcm_response_image_closure_classified = FALSE,
      success_mass = success,
      failure_mass = failure,
      independent_contrast_optimum = independent,
      independent_difference_owner_1 = if (is.null(independent_state)) {
        NA_real_
      } else {
        independent_state$differences[1]
      },
      independent_difference_owner_2 = if (is.null(independent_state)) {
        NA_real_
      } else {
        independent_state$differences[2]
      },
      independent_optimum_image_state = if (is.null(independent_state)) {
        "not_evaluated"
      } else {
        independent_state$state
      },
      independent_optimum_in_finite_response_image = finite_exists,
      independent_optimum_in_finite_response_image_closure = isTRUE(
        independent_state$finite_response_image_closure
      ),
      independent_optimum_on_missing_finite_boundary = isTRUE(
        independent_state$missing_finite_boundary
      ),
      closure_strata = strata,
      boundary_axis_candidates = axes,
      boundary_axis_candidate_count = as.integer(nrow(axes)),
      boundary_envelope_winner_count = if (nrow(axes) == 0L) {
        0L
      } else {
        as.integer(sum(axes$BoundaryEnvelopeWinner))
      },
      global_closure_maximizers = global_candidates,
      global_closure_maximizer_count = as.integer(nrow(global_candidates)),
      global_supremum_value = as.numeric(global_supremum),
      boundary_limsup_supremum_value = as.numeric(boundary_supremum),
      finite_over_boundary_log_likelihood_gap = as.numeric(boundary_gap),
      strict_boundary_limsup_gap_certified = finite_exists &&
        !is.null(independent_state) &&
        independent_state$state %in% c(
          "finite_positive_difference_interior",
          "finite_negative_difference_interior"
        ) && is.finite(boundary_gap) && boundary_gap > 0,
      finite_jmle_existence_certified = finite_exists,
      finite_jmle_nonexistence_certified = finite_nonexists,
      finite_jmle_nonexistence_without_cell_separation_certified =
        finite_nonexists,
      response_image_supremum_attained_by_finite_parameter = finite_exists,
      finite_global_representative = finite_representative,
      competitive_boundary_certified = isTRUE(classified) &&
        is.finite(boundary_gap) && boundary_gap == 0,
      global_supremum_value_classified = isTRUE(classified),
      complete_escaping_sequence_boundary_envelope_constructed =
        isTRUE(classified),
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_theorem_only",
      next_gate = if (isTRUE(classified)) {
        "general_design_response_image_closure_stratification"
      } else {
        "unclassified_minimal_binary_closure_envelope"
      },
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "This is a complete closure and likelihood-envelope theorem only for",
        "the fixed binary two-Person, two-slope-owner affine operator with",
        "positive success and failure mass in all four cells. It does not",
        "classify a general retained GPCM design, MML, uncertainty, readiness,",
        "recovery, simulation, or FACETS comparison."
      )
    )
  }

  control <- mfrmr_jml_gpcm_exponential_balance_scalar(max_total_mass)
  success <- tryCatch(
    suppressWarnings(as.numeric(success_mass)),
    error = function(e) NULL
  )
  failure <- tryCatch(
    suppressWarnings(as.numeric(failure_mass)),
    error = function(e) NULL
  )
  valid_control <- is.finite(control) && control > 0
  if (!valid_control) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "The positive total response-mass workload control was invalid.",
      "binary_closure_envelope_control_invalid"
    ))
  }
  if (!source_matches) {
    return(finish(
      "not_evaluated_source", FALSE, FALSE,
      "The exact P2h binary response-quotient source was unavailable.",
      "binary_closure_envelope_source_invalid"
    ))
  }
  valid_mass <- !is.null(success) && !is.null(failure) &&
    length(success) == 4L && length(failure) == 4L &&
    all(is.finite(success)) && all(is.finite(failure)) &&
    all(success > 0) && all(failure > 0) &&
    sum(c(success, failure)) <= control
  if (!valid_mass) {
    return(finish(
      "not_evaluated_response_mass", FALSE, FALSE,
      paste(
        "Four finite strictly positive success and failure masses within the",
        "workload limit are required."
      ),
      "binary_closure_envelope_response_mass_invalid"
    ))
  }

  independent <- log(success / failure)
  independent_state <-
    mfrmr_jml_gpcm_binary_closure_image_state(independent)
  pooled_12 <- log(sum(success[1:2]) / sum(failure[1:2]))
  pooled_34 <- log(sum(success[3:4]) / sum(failure[3:4]))
  axis_contrasts <- rbind(
    owner_1_dominant = c(independent[1:2], pooled_34, pooled_34),
    owner_2_dominant = c(pooled_12, pooled_12, independent[3:4])
  )
  axis_states <- lapply(seq_len(nrow(axis_contrasts)), function(index) {
    mfrmr_jml_gpcm_binary_closure_image_state(axis_contrasts[index, ])
  })
  axis_loglik <- apply(axis_contrasts, 1L, function(value) {
    mfrmr_jml_gpcm_binary_closure_loglik(value, success, failure)
  })
  boundary_supremum <- max(axis_loglik)
  boundary_winner <- axis_loglik == boundary_supremum
  axes <- data.frame(
    Axis = c("owner_1_dominant", "owner_2_dominant"),
    DivergentSlopeOwner = c(1L, 2L),
    Constraint = c("d2=0", "d1=0"),
    Z1 = axis_contrasts[, 1],
    Z2 = axis_contrasts[, 2],
    Z3 = axis_contrasts[, 3],
    Z4 = axis_contrasts[, 4],
    DifferenceOwner1 = vapply(
      axis_states, function(value) value$differences[1], numeric(1)
    ),
    DifferenceOwner2 = vapply(
      axis_states, function(value) value$differences[2], numeric(1)
    ),
    LogLikelihood = as.numeric(axis_loglik),
    FiniteResponseImage = vapply(
      axis_states, function(value) value$finite_response_image, logical(1)
    ),
    MissingFiniteBoundary = vapply(
      axis_states, function(value) value$missing_finite_boundary, logical(1)
    ),
    BoundaryEnvelopeWinner = boundary_winner,
    GlobalClosureMaximizer = FALSE,
    stringsAsFactors = FALSE
  )

  if (isTRUE(independent_state$finite_response_image)) {
    global_supremum <- mfrmr_jml_gpcm_binary_closure_loglik(
      independent, success, failure
    )
    global_candidates <- matrix(independent, nrow = 1L, ncol = 4L)
    finite_representative <-
      mfrmr_jml_gpcm_binary_closure_representative(independent)
    state <- if (all(independent_state$differences == 0)) {
      "finite_jmle_exists_on_response_equivalence_intersection"
    } else {
      "finite_jmle_exists_with_strict_boundary_gap"
    }
  } else {
    global_supremum <- boundary_supremum
    global_candidates <- axis_contrasts[boundary_winner, , drop = FALSE]
    finite_representative <-
      mfrmr_jml_gpcm_binary_closure_representative(numeric(0))
    state <- "finite_jmle_nonexistence_missing_response_boundary_certified"
  }
  axes$GlobalClosureMaximizer <- axis_loglik == global_supremum

  finish(
    state, TRUE, TRUE,
    paste(
      "The exact finite image and its two-axis closure were classified.",
      "Positive success and failure mass makes every unbounded contrast path",
      "tend to minus infinity. Pooling each axis-constrained pair gives the",
      "complete bounded-response escape envelope and the global response-image",
      "supremum."
    ),
    paste(
      "minimal_binary_response_image_closure_complete",
      if (isTRUE(independent_state$finite_response_image)) {
        "finite_global_response_optimum_representable"
      } else {
        "global_supremum_only_on_missing_finite_boundary"
      },
      sep = ";"
    ),
    success, failure, independent, independent_state,
    mfrmr_jml_gpcm_binary_closure_strata(), axes,
    global_supremum, boundary_supremum, global_candidates,
    finite_representative
  )
}

mfrmr_jml_gpcm_binary_closure_path_at <- function(
    result,
    axis = c("owner_1_dominant", "owner_2_dominant"),
    path_index,
    max_path_index = 300) {
  axis <- tryCatch(match.arg(axis), error = function(e) NA_character_)
  index <- mfrmr_jml_gpcm_exponential_balance_scalar(path_index)
  control <- mfrmr_jml_gpcm_exponential_balance_scalar(max_path_index)
  valid_result <- is.list(result) && identical(
    as.character(result$contract_version %||% "")[1],
    mfrmr_jml_gpcm_binary_closure_contract_version()
  ) && isTRUE(result$binary_closure_envelope_classified)
  valid <- valid_result && !is.na(axis) && is.finite(index) && index >= 0 &&
    is.finite(control) && control >= 0 && index <= control
  if (!valid) {
    return(list(
      valid = FALSE,
      axis = axis,
      path_index = index,
      contrasts = numeric(0),
      log_likelihood = NA_real_,
      target_log_likelihood = NA_real_,
      target_minus_path_log_likelihood = NA_real_,
      contrast_distance_to_target = NA_real_,
      reason_codes = "binary_closure_path_invalid"
    ))
  }
  row <- result$boundary_axis_candidates[
    result$boundary_axis_candidates$Axis == axis, , drop = FALSE
  ]
  if (nrow(row) != 1L) {
    return(list(
      valid = FALSE,
      axis = axis,
      path_index = index,
      contrasts = numeric(0),
      log_likelihood = NA_real_,
      target_log_likelihood = NA_real_,
      target_minus_path_log_likelihood = NA_real_,
      contrast_distance_to_target = NA_real_,
      reason_codes = "binary_closure_path_axis_missing"
    ))
  }
  target <- as.numeric(row[1, c("Z1", "Z2", "Z3", "Z4")])
  if (identical(axis, "owner_1_dominant")) {
    log_alpha_1 <- index
    alpha_1 <- exp(index)
    difference <- target[1] - target[2]
    contrasts <- c(
      target[1:2],
      target[3] + difference * exp(-2 * index) / 2,
      target[4] - difference * exp(-2 * index) / 2
    )
  } else {
    log_alpha_1 <- -index
    alpha_1 <- exp(-index)
    difference <- target[3] - target[4]
    contrasts <- c(
      target[1] + difference * exp(-2 * index) / 2,
      target[2] - difference * exp(-2 * index) / 2,
      target[3:4]
    )
  }
  theta <- c(
    (contrasts[1] / alpha_1 + alpha_1 * contrasts[3]) / 2,
    (contrasts[2] / alpha_1 + alpha_1 * contrasts[4]) / 2
  )
  delta <-
    (alpha_1 * contrasts[3] - contrasts[1] / alpha_1) / 2
  reconstructed <- c(
    alpha_1 * (theta[1] - delta),
    alpha_1 * (theta[2] - delta),
    (theta[1] + delta) / alpha_1,
    (theta[2] + delta) / alpha_1
  )
  log_likelihood <- mfrmr_jml_gpcm_binary_closure_loglik(
    contrasts, result$success_mass, result$failure_mass
  )
  target_log_likelihood <- as.numeric(row$LogLikelihood)
  list(
    valid = all(is.finite(c(
      contrasts, theta, delta, reconstructed, log_likelihood
    ))),
    axis = axis,
    divergent_slope_owner = as.integer(row$DivergentSlopeOwner),
    path_index = index,
    parameter_path_declared_divergent = TRUE,
    log_alpha = c(log_alpha_1, -log_alpha_1),
    alpha = c(alpha_1, 1 / alpha_1),
    theta = theta,
    delta = delta,
    contrasts = contrasts,
    reconstructed_contrasts = reconstructed,
    maximum_absolute_reconstruction_residual = max(
      abs(reconstructed - contrasts)
    ),
    target_contrasts = target,
    contrast_distance_to_target = sqrt(sum((contrasts - target)^2)),
    log_likelihood = log_likelihood,
    target_log_likelihood = target_log_likelihood,
    target_minus_path_log_likelihood =
      target_log_likelihood - log_likelihood,
    reason_codes = "binary_closure_finite_boundary_path_evaluated"
  )
}

mfrmr_jml_gpcm_binary_closure_nonattainment_fixture <- function() {
  mfrmr_jml_gpcm_binary_closure_envelope(
    success_mass = c(2, 4, 1, 1),
    failure_mass = c(1, 1, 1, 1)
  )
}

mfrmr_jml_gpcm_fit_binary_closure_envelope <- function(
    fit,
    source_response_quotient = NULL,
    max_total_mass = 1e9) {
  finish <- function(state, evaluated, reconstructed, detail, reason_codes,
                     result = NULL, success = numeric(0),
                     failure = numeric(0), operator_rows = matrix(
                       numeric(0), nrow = 0L, ncol = 0L
                     ), retained_contrasts = numeric(0),
                     reconstructed_retained_loglik = NA_real_) {
    classified <- isTRUE(result$binary_closure_envelope_classified)
    retained_loglik <- if (
      inherits(fit, "mfrm_fit") && is.finite(fit$opt$value %||% NA_real_)
    ) {
      -as.numeric(fit$opt$value)
    } else {
      NA_real_
    }
    supremum <- if (classified) {
      as.numeric(result$global_supremum_value)
    } else {
      NA_real_
    }
    list(
      contract_version =
        mfrmr_jml_gpcm_binary_closure_contract_version(),
      estimator_identity = if (isTRUE(evaluated)) {
        "unpenalized_fixed_effects_jml_no_finite_box"
      } else {
        "not_classified"
      },
      objective_identity = if (isTRUE(evaluated)) {
        "identified_conditional_joint_log_likelihood"
      } else {
        "not_classified"
      },
      state = state,
      evaluated = isTRUE(evaluated),
      current_fit_reconstructed = isTRUE(reconstructed),
      exact_canonical_affine_operator_certified = classified,
      complete_two_person_by_two_owner_cell_support = classified,
      response_cell_order = c(
        "person_1_owner_1", "person_2_owner_1",
        "person_1_owner_2", "person_2_owner_2"
      ),
      response_cell_success_mass = success,
      response_cell_failure_mass = failure,
      canonical_adjacent_operator_rows = operator_rows,
      retained_response_cell_contrasts = retained_contrasts,
      retained_log_likelihood = retained_loglik,
      reconstructed_retained_log_likelihood = as.numeric(
        reconstructed_retained_loglik
      ),
      retained_objective_reconstruction_residual = if (
        is.finite(retained_loglik) &&
          is.finite(reconstructed_retained_loglik)
      ) {
        reconstructed_retained_loglik - retained_loglik
      } else {
        NA_real_
      },
      current_fit_objective_numerically_reconstructed =
        is.finite(retained_loglik) &&
        is.finite(reconstructed_retained_loglik) &&
        abs(reconstructed_retained_loglik - retained_loglik) <= 1e-10,
      classified_global_supremum = supremum,
      supremum_minus_retained_log_likelihood = if (
        is.finite(retained_loglik) && is.finite(supremum)
      ) {
        supremum - retained_loglik
      } else {
        NA_real_
      },
      finite_jmle_existence_certified = isTRUE(
        result$finite_jmle_existence_certified
      ),
      finite_jmle_nonexistence_certified = isTRUE(
        result$finite_jmle_nonexistence_certified
      ),
      finite_jmle_nonexistence_without_cell_separation_certified = isTRUE(
        result$finite_jmle_nonexistence_without_cell_separation_certified
      ),
      complete_escaping_sequence_boundary_envelope_constructed = isTRUE(
        result$complete_escaping_sequence_boundary_envelope_constructed
      ),
      optimizer_trace_promoted_to_finite_estimate = FALSE,
      result = result,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_theorem_only",
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The wrapper recognizes only the exact unanchored binary two-Person,",
        "two-owner retained operator with one common step/slope facet and no",
        "interaction. It performs no general design search and does not alter",
        "the fit's stored production readiness or primary parameter fields."
      )
    )
  }

  if (!inherits(fit, "mfrm_fit") || is.null(fit$config) ||
      is.null(fit$prep) || is.null(fit$opt)) {
    return(finish(
      "not_evaluated_fit", FALSE, FALSE,
      "A current retained mfrm_fit is required.",
      "binary_closure_fit_invalid"
    ))
  }
  config <- fit$config
  if (!identical(config$model, "GPCM")) {
    return(finish(
      "not_applicable_model", FALSE, FALSE,
      "The exact binary closure wrapper applies only to GPCM.",
      "binary_closure_not_applicable_model"
    ))
  }
  if (!identical(config$method, "JML")) {
    return(finish(
      "not_applicable_estimator", FALSE, FALSE,
      "The fixed-effects response-image closure is not reused for MML.",
      "binary_closure_not_applicable_estimator"
    ))
  }

  built <- tryCatch({
    idx <- build_indices(
      fit$prep,
      step_facet = config$step_facet,
      slope_facet = config$slope_facet,
      interaction_specs = config$interaction_specs
    )
    sizes <- build_param_sizes(config)
    adjacent <- mfrmr_estimability_adjacent_design(
      fit$prep, idx, config, sizes,
      include_person = TRUE, include_population_beta = FALSE
    )
    zero_kernel <- mfrmr_gpcm_response_kernel_design(
      fit$prep, idx, config, sizes, rep(0, length(fit$opt$par))
    )
    params <- expand_params(fit$opt$par, sizes, config)
    list(
      idx = idx, sizes = sizes, adjacent = adjacent,
      zero_kernel = zero_kernel, params = params
    )
  }, error = function(e) e)
  if (inherits(built, "error")) {
    return(finish(
      "not_evaluated_design", FALSE, FALSE,
      paste0(
        "The retained binary closure operator failed: ",
        conditionMessage(built)
      ),
      "binary_closure_retained_design_invalid"
    ))
  }

  design <- tryCatch(
    as.matrix(built$adjacent$design),
    error = function(e) NULL
  )
  offset <- tryCatch(
    suppressWarnings(as.numeric(built$zero_kernel$eta_minus_step)),
    error = function(e) NULL
  )
  person <- tryCatch(
    suppressWarnings(as.integer(built$idx$person)),
    error = function(e) NULL
  )
  slope <- tryCatch(
    suppressWarnings(as.integer(built$idx$slope_idx)),
    error = function(e) NULL
  )
  score <- tryCatch(
    suppressWarnings(as.integer(built$idx$score_k)),
    error = function(e) NULL
  )
  weight <- tryCatch(
    suppressWarnings(as.numeric(built$idx$weight)),
    error = function(e) NULL
  )
  n_observations <- length(score)
  structural_scope <- length(fit$prep$facet_names %||% character(0)) == 1L &&
    identical(config$step_facet, config$slope_facet) &&
    length(config$interaction_specs %||% list()) == 0L
  valid_vectors <- structural_scope && !is.null(design) &&
    !is.null(offset) && !is.null(person) && !is.null(slope) &&
    !is.null(score) && !is.null(weight) && n_observations >= 8L &&
    nrow(design) == n_observations && ncol(design) == 3L &&
    length(offset) == n_observations && all(offset == 0) &&
    length(person) == n_observations &&
    length(slope) == n_observations &&
    length(weight) == n_observations && all(is.finite(design)) &&
    all(is.finite(weight)) && all(weight > 0) &&
    identical(sort(unique(person)), 1:2) &&
    identical(sort(unique(slope)), 1:2) &&
    all(score %in% 0:1)
  if (!valid_vectors) {
    return(finish(
      "not_evaluated_exact_operator", TRUE, FALSE,
      paste(
        "The retained fit was outside the exact unanchored binary",
        "two-Person, two-owner affine-operator scope."
      ),
      "binary_closure_exact_operator_unavailable"
    ))
  }

  cell <- (slope - 1L) * 2L + person
  complete_cells <- identical(sort(unique(cell)), 1:4)
  operator_rows <- matrix(NA_real_, nrow = 4L, ncol = 3L)
  cell_constant <- complete_cells
  if (complete_cells) {
    for (index in 1:4) {
      rows <- which(cell == index)
      operator_rows[index, ] <- design[rows[1], ]
      cell_constant <- cell_constant && all(
        design[rows, , drop = FALSE] == matrix(
          operator_rows[index, ], nrow = length(rows), ncol = 3L,
          byrow = TRUE
        )
      )
    }
  }
  canonical <- rbind(
    c(1, 0, -1),
    c(0, 1, -1),
    c(1, 0, 1),
    c(0, 1, 1)
  )
  exact_operator <- complete_cells && cell_constant &&
    all(operator_rows == canonical)
  if (!exact_operator) {
    return(finish(
      "not_evaluated_exact_operator", TRUE, TRUE,
      "The grouped retained adjacent operator did not match the P2i fixture.",
      "binary_closure_exact_operator_mismatch",
      operator_rows = operator_rows
    ))
  }

  success <- vapply(1:4, function(index) {
    sum(weight[cell == index & score == 1L])
  }, numeric(1))
  failure <- vapply(1:4, function(index) {
    sum(weight[cell == index & score == 0L])
  }, numeric(1))
  optimizer_index <- tryCatch(
    suppressWarnings(as.integer(built$adjacent$map$OptimizerIndex)),
    error = function(e) integer(0)
  )
  retained_contrasts <- numeric(0)
  reconstructed_retained_loglik <- NA_real_
  valid_retained_map <- length(optimizer_index) == 3L &&
    all(optimizer_index >= 1L & optimizer_index <= length(fit$opt$par)) &&
    length(built$params$slopes) == 2L &&
    all(is.finite(built$params$slopes))
  if (valid_retained_map) {
    retained_additive <- fit$opt$par[optimizer_index]
    retained_adjacent <- offset + as.numeric(design %*% retained_additive)
    retained_row_contrasts <-
      built$params$slopes[slope] * retained_adjacent
    retained_contrasts <- vapply(1:4, function(index) {
      retained_row_contrasts[which(cell == index)[1]]
    }, numeric(1))
    reconstructed_retained_loglik <-
      mfrmr_jml_gpcm_binary_closure_loglik(
        retained_contrasts, success, failure
      )
  }
  result_record <- mfrmr_jml_gpcm_binary_closure_envelope(
    success_mass = success,
    failure_mass = failure,
    source_response_quotient = source_response_quotient,
    max_total_mass = max_total_mass
  )
  if (!isTRUE(result_record$binary_closure_envelope_classified)) {
    return(finish(
      "not_evaluated_response_envelope", TRUE, TRUE,
      "The exact retained response masses did not yield a P2i envelope.",
      result_record$reason_codes,
      result = result_record,
      success = success,
      failure = failure,
      operator_rows = operator_rows,
      retained_contrasts = retained_contrasts,
      reconstructed_retained_loglik = reconstructed_retained_loglik
    ))
  }
  finish(
    paste0("current_fit_", result_record$state), TRUE, TRUE,
    paste(
      "The exact retained binary operator and four expanded response-cell",
      "masses were reconstructed before applying the complete P2i envelope."
    ),
    result_record$reason_codes,
    result = result_record,
    success = success,
    failure = failure,
    operator_rows = operator_rows,
    retained_contrasts = retained_contrasts,
    reconstructed_retained_loglik = reconstructed_retained_loglik
  )
}
