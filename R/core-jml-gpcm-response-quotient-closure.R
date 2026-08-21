# JML GPCM response-equivalence quotient closure obstruction
# ==============================================================================
#
# Collapsing finite parameter points which induce the same response contrasts
# removes flat fibres, but it does not add limits which are absent from the
# finite response image.  This contract consumes one exact P2g bounded-image
# escape and an exact affine row relation.  The relation can prove that the
# finite contrast limit is in the closure of the response image but is not in
# the image itself.

mfrmr_jml_gpcm_response_quotient_contract_version <- function() {
  "mfrmr-jml-gpcm-response-quotient-closure-0.2.3-v1"
}

mfrmr_jml_gpcm_response_quotient_indices <- function(value) {
  converted <- tryCatch(
    suppressWarnings(as.numeric(value)),
    error = function(e) NULL
  )
  if (is.null(converted) || length(converted) < 1L ||
      any(!is.finite(converted)) || any(converted != floor(converted))) {
    return(NULL)
  }
  as.integer(converted)
}

mfrmr_jml_gpcm_response_quotient_closure_obstruction <- function(
    exponential_balance_result,
    zero_rows,
    witness_rows,
    witness_coefficients,
    zero_relation_coefficients,
    max_certificate_rows = 1000L,
    max_certificate_elements = 1e6) {
  source <- exponential_balance_result
  source_matches <- is.list(source) && identical(
    as.character(source$contract_version %||% "")[1],
    mfrmr_jml_gpcm_exponential_balance_contract_version()
  )
  source_record <- if (is.list(source)) source else list()
  empty_matrix <- matrix(numeric(0), nrow = 0L, ncol = 0L)

  finish <- function(state, evaluated, certified, detail, reason_codes,
                     target_cumulative = empty_matrix,
                     target_adjacent = numeric(0),
                     zero_index = integer(0),
                     witness_index = integer(0),
                     witness_owner = NA_integer_,
                     relation_residual = numeric(0),
                     witness_target_combination = NA_real_) {
    owner_certified <- isTRUE(certified) || state %in% c(
      "not_certified_nonzero_declared_zero_target",
      "not_certified_affine_relation",
      "not_certified_zero_target_witness"
    )
    zero_rows_certified <- isTRUE(certified) || state %in% c(
      "not_certified_affine_relation",
      "not_certified_zero_target_witness"
    )
    relation_checked <- zero_rows_certified
    relation_certified <- isTRUE(certified) || identical(
      state, "not_certified_zero_target_witness"
    )
    list(
      contract_version =
        mfrmr_jml_gpcm_response_quotient_contract_version(),
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      response_equivalence_identity = paste(
        "equality of the complete finite cumulative category-logit",
        "contrast vector"
      ),
      quotient_metric_identity =
        "Euclidean response-contrast metric on finite equivalence classes",
      state = state,
      evaluated = isTRUE(evaluated),
      response_quotient_closure_obstruction_certified = isTRUE(certified),
      source_exponential_balance_contract_matches = isTRUE(source_matches),
      source_exponential_balance_limit_classified = isTRUE(
        source_record$exponential_balance_limit_classified
      ),
      source_bounded_response_image_escape_certified = isTRUE(
        source_record$bounded_response_image_escape_certified
      ),
      exact_affine_row_relation_checked = isTRUE(relation_checked),
      exact_affine_row_relation_certified = isTRUE(relation_certified),
      zero_target_rows_certified = isTRUE(zero_rows_certified),
      common_slope_owner_witness_certified = isTRUE(owner_certified),
      nonzero_target_witness_certified = isTRUE(certified),
      finite_contrast_target_in_response_image_closure = isTRUE(evaluated),
      finite_contrast_target_representable_by_finite_parameter = if (
        isTRUE(certified)
      ) FALSE else NA,
      finite_parameter_response_image_closed = if (
        isTRUE(certified)
      ) FALSE else NA,
      response_metric_quotient_complete = if (
        isTRUE(certified)
      ) FALSE else NA,
      induced_quotient_to_contrast_map_proper = if (
        isTRUE(certified)
      ) FALSE else NA,
      collapsing_finite_response_equivalence_fibres_suffices = if (
        isTRUE(certified)
      ) FALSE else NA,
      quotient_completion_or_boundary_strata_required = isTRUE(certified),
      target_cumulative_contrasts = target_cumulative,
      target_adjacent_contrasts = as.numeric(target_adjacent),
      zero_rows = as.integer(zero_index),
      witness_rows = as.integer(witness_index),
      witness_slope_owner = as.integer(witness_owner),
      affine_relation_residual = as.numeric(relation_residual),
      witness_target_linear_combination = as.numeric(
        witness_target_combination
      ),
      boundary_log_likelihood_continuation = if (
        isTRUE(certified)
      ) {
        as.numeric(source_record$log_likelihood_limit)
      } else {
        NA_real_
      },
      complete_response_image_closure_stratified = FALSE,
      complete_bounded_image_escape_family_classified = FALSE,
      complete_escaping_sequence_boundary_envelope_constructed = FALSE,
      competitive_boundary_certified = FALSE,
      global_supremum_value_classified = FALSE,
      finite_jmle_existence_certified = FALSE,
      finite_jmle_nonexistence_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_theorem_only",
      next_gate = if (isTRUE(certified)) {
        paste(
          "complete_response_image_closure_stratification_and_boundary_",
          "limsup_envelope",
          sep = ""
        )
      } else {
        "unclassified_response_quotient_closure"
      },
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The certificate proves that one finite response-contrast target is",
        "missing from the finite image although it is approached by a P2g",
        "path. It refutes a fibre-collapse-only shortcut. It neither",
        "classifies every design nor proves that this target is competitive,",
        "that a finite JMLE fails to exist, or that the boundary envelope is",
        "complete. It is not reused for MML or FACETS comparison."
      )
    )
  }

  controls <- suppressWarnings(as.numeric(c(
    max_certificate_rows, max_certificate_elements
  )))
  valid_controls <- length(controls) == 2L && all(is.finite(controls)) &&
    all(controls == floor(controls)) && all(controls >= 1)
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "The exact quotient-certificate workload controls were invalid.",
      "response_quotient_control_invalid"
    ))
  }
  valid_source <- source_matches &&
    isTRUE(source_record$exponential_balance_limit_classified) &&
    isTRUE(source_record$bounded_response_image_escape_certified) &&
    isTRUE(source_record$response_contrast_image_bounded) &&
    isTRUE(source_record$parameter_path_divergent) &&
    identical(as.integer(source_record$positive_contrast_exponent_groups), 0L)
  if (!valid_source) {
    return(finish(
      "not_evaluated_source", FALSE, FALSE,
      paste(
        "A classified P2g bounded-response-image parameter escape with no",
        "positive contrast exponent is required."
      ),
      "response_quotient_source_invalid"
    ))
  }

  path <- source_record$path_definition
  design <- tryCatch(
    methods::as(
      methods::as(
        Matrix::Matrix(path$adjacent_design, sparse = TRUE), "generalMatrix"
      ),
      "CsparseMatrix"
    ),
    error = function(e) NULL
  )
  offset <- tryCatch(
    suppressWarnings(as.numeric(path$additive_offset)),
    error = function(e) NULL
  )
  slope_index <- tryCatch(
    suppressWarnings(as.numeric(path$slope_index)),
    error = function(e) NULL
  )
  target_cumulative <- tryCatch(
    as.matrix(source_record$contrast_flag_limit$base_contrasts),
    error = function(e) NULL
  )
  n_observations <- as.integer(
    mfrmr_jml_gpcm_exponential_balance_scalar(source_record$observations)
  )
  n_categories <- as.integer(
    mfrmr_jml_gpcm_exponential_balance_scalar(source_record$categories)
  )
  slope_levels <- as.integer(
    mfrmr_jml_gpcm_exponential_balance_scalar(source_record$slope_levels)
  )
  valid_dimensions <- length(n_observations) == 1L &&
    length(n_categories) == 1L && length(slope_levels) == 1L &&
    is.finite(n_observations) && is.finite(n_categories) &&
    is.finite(slope_levels) && n_observations >= 1L &&
    n_categories >= 2L && slope_levels >= 2L
  expected_rows <- if (valid_dimensions) {
    n_observations * (n_categories - 1L)
  } else {
    NA_integer_
  }
  valid_path <- valid_dimensions && !is.null(design) && !is.null(offset) &&
    !is.null(slope_index) && !is.null(target_cumulative) &&
    nrow(design) == expected_rows && length(offset) == expected_rows &&
    length(slope_index) == n_observations &&
    nrow(target_cumulative) == n_observations &&
    ncol(target_cumulative) == n_categories - 1L &&
    all(is.finite(design@x)) && all(is.finite(offset)) &&
    all(is.finite(slope_index)) &&
    all(slope_index == floor(slope_index)) && all(slope_index >= 1) &&
    all(slope_index <= slope_levels) &&
    all(is.finite(target_cumulative))
  if (!valid_path) {
    return(finish(
      "not_evaluated_path", FALSE, FALSE,
      "The P2g affine design, slope map, or finite contrast target was invalid.",
      "response_quotient_path_invalid"
    ))
  }

  target_adjacent_matrix <- target_cumulative
  if (n_categories > 2L) {
    target_adjacent_matrix[, 2:(n_categories - 1L)] <-
      target_cumulative[, 2:(n_categories - 1L), drop = FALSE] -
      target_cumulative[, 1:(n_categories - 2L), drop = FALSE]
  }
  target_adjacent <- as.numeric(target_adjacent_matrix)
  zero_index <- mfrmr_jml_gpcm_response_quotient_indices(zero_rows)
  witness_index <- mfrmr_jml_gpcm_response_quotient_indices(witness_rows)
  witness_coefficients <- tryCatch(
    suppressWarnings(as.numeric(witness_coefficients)),
    error = function(e) NULL
  )
  zero_coefficients <- tryCatch(
    suppressWarnings(as.numeric(zero_relation_coefficients)),
    error = function(e) NULL
  )
  all_indices <- c(witness_index, zero_index)
  requested_elements <- length(unique(all_indices)) * (ncol(design) + 1L)
  valid_declaration <- !is.null(zero_index) && !is.null(witness_index) &&
    length(zero_index) <= controls[1] &&
    length(witness_index) <= controls[1] &&
    length(unique(zero_index)) == length(zero_index) &&
    length(unique(witness_index)) == length(witness_index) &&
    length(intersect(zero_index, witness_index)) == 0L &&
    all(all_indices >= 1L & all_indices <= expected_rows) &&
    !is.null(witness_coefficients) && !is.null(zero_coefficients) &&
    length(witness_coefficients) == length(witness_index) &&
    length(zero_coefficients) == length(zero_index) &&
    all(is.finite(witness_coefficients)) &&
    all(is.finite(zero_coefficients)) &&
    any(witness_coefficients != 0) && any(zero_coefficients != 0) &&
    is.finite(requested_elements) && requested_elements <= controls[2]
  if (!valid_declaration) {
    return(finish(
      "not_evaluated_declaration", FALSE, FALSE,
      "The exact zero-row or affine-relation witness declaration was invalid.",
      "response_quotient_declaration_invalid",
      target_cumulative = target_cumulative,
      target_adjacent = target_adjacent
    ))
  }

  slope_owner_by_row <- rep(
    as.integer(slope_index), times = n_categories - 1L
  )
  witness_owners <- unique(slope_owner_by_row[witness_index])
  if (length(witness_owners) != 1L) {
    return(finish(
      "not_certified_mixed_witness_owner", TRUE, FALSE,
      "All nonzero witness rows must share one positive finite slope.",
      "response_quotient_witness_owner_mismatch",
      target_cumulative, target_adjacent, zero_index, witness_index
    ))
  }
  if (!all(target_adjacent[zero_index] == 0)) {
    return(finish(
      "not_certified_nonzero_declared_zero_target", TRUE, FALSE,
      "Every declared zero row must have an implementation-exact zero target.",
      "response_quotient_zero_target_mismatch",
      target_cumulative, target_adjacent, zero_index, witness_index,
      witness_owners
    ))
  }

  selected_design <- as.matrix(
    design[all_indices, , drop = FALSE]
  )
  selected_augmented <- cbind(offset[all_indices], selected_design)
  witness_positions <- seq_along(witness_index)
  zero_positions <- length(witness_index) + seq_along(zero_index)
  relation_residual <-
    colSums(sweep(
      selected_augmented[witness_positions, , drop = FALSE],
      1L, witness_coefficients, `*`
    )) -
    colSums(sweep(
      selected_augmented[zero_positions, , drop = FALSE],
      1L, zero_coefficients, `*`
    ))
  if (!all(relation_residual == 0)) {
    return(finish(
      "not_certified_affine_relation", TRUE, FALSE,
      "The declared affine utility-row relation did not hold exactly.",
      "response_quotient_affine_relation_mismatch",
      target_cumulative, target_adjacent, zero_index, witness_index,
      witness_owners, relation_residual
    ))
  }
  witness_target_combination <- sum(
    witness_coefficients * target_adjacent[witness_index]
  )
  if (witness_target_combination == 0 ||
      !is.finite(witness_target_combination)) {
    return(finish(
      "not_certified_zero_target_witness", TRUE, FALSE,
      "The common-owner target combination must be exactly nonzero.",
      "response_quotient_target_witness_zero",
      target_cumulative, target_adjacent, zero_index, witness_index,
      witness_owners, relation_residual, witness_target_combination
    ))
  }

  finish(
    "finite_response_image_not_closed_certified", TRUE, TRUE,
    paste(
      "The P2g path approaches the finite target. Exact zero target rows",
      "force their unscaled affine utilities to zero at every hypothetical",
      "finite representative. Their declared affine relation then forces a",
      "common-slope witness combination to zero, contradicting its nonzero",
      "target value. The target is in the response-image closure but not the",
      "finite response image."
    ),
    paste(
      "finite_response_image_closure_obstruction",
      "finite_response_equivalence_fibre_collapse_insufficient",
      sep = ";"
    ),
    target_cumulative, target_adjacent, zero_index, witness_index,
    witness_owners, relation_residual, witness_target_combination
  )
}

mfrmr_jml_gpcm_response_quotient_binary_counterexample <- function() {
  source <- mfrmr_jml_gpcm_exponential_balance_limit(
    adjacent_design = rbind(
      c(1, 0, -1),
      c(0, 1, -1),
      c(1, 0, 1),
      c(0, 1, 1)
    ),
    additive_center = c(0, 0, 0),
    additive_coordinate_directions = list(c(1, 2, 0)),
    additive_exponential_rates = -1,
    score_k = c(1L, 1L, 0L, 0L),
    slope_index = c(1L, 1L, 2L, 2L),
    base_log_slopes = c(0, 0),
    slope_log_rates = c(1, -1),
    weight = rep(1, 4)
  )
  certificate <- mfrmr_jml_gpcm_response_quotient_closure_obstruction(
    source,
    zero_rows = c(3L, 4L),
    witness_rows = c(1L, 2L),
    witness_coefficients = c(1, -1),
    zero_relation_coefficients = c(1, -1)
  )
  certificate$canonical_fixture_identity <-
    "binary_two_person_two_slope_owner_centered_location_gpcm"
  certificate$canonical_path_formula <- paste(
    "theta_1(t)=exp(-t), theta_2(t)=2*exp(-t), delta(t)=0;",
    "alpha_1(t)=exp(t), alpha_2(t)=exp(-t)"
  )
  certificate$source_exponential_balance <- source
  certificate
}
