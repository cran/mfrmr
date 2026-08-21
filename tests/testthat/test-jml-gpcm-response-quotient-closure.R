canonical_response_quotient_source <- function(
    first = 1,
    second = 2,
    additive_offset = rep(0, 4),
    additive_center = c(0, 0, 0)) {
  mfrmr:::mfrmr_jml_gpcm_exponential_balance_limit(
    adjacent_design = rbind(
      c(1, 0, -1),
      c(0, 1, -1),
      c(1, 0, 1),
      c(0, 1, 1)
    ),
    additive_center = additive_center,
    additive_offset = additive_offset,
    additive_coordinate_directions = list(c(first, second, 0)),
    additive_exponential_rates = -1,
    score_k = c(1L, 1L, 0L, 0L),
    slope_index = c(1L, 1L, 2L, 2L),
    base_log_slopes = c(0, 0),
    slope_log_rates = c(1, -1),
    weight = rep(1, 4)
  )
}

call_response_quotient_obstruction <- function(source, ...) {
  defaults <- list(
    exponential_balance_result = source,
    zero_rows = c(3L, 4L),
    witness_rows = c(1L, 2L),
    witness_coefficients = c(1, -1),
    zero_relation_coefficients = c(1, -1)
  )
  do.call(
    mfrmr:::mfrmr_jml_gpcm_response_quotient_closure_obstruction,
    utils::modifyList(defaults, list(...))
  )
}

test_that("a binary GPCM response-equivalence quotient image is not closed", {
  result <- mfrmr:::mfrmr_jml_gpcm_response_quotient_binary_counterexample()
  source <- result$source_exponential_balance

  expect_identical(
    result$contract_version,
    "mfrmr-jml-gpcm-response-quotient-closure-0.2.3-v1"
  )
  expect_identical(
    result$state, "finite_response_image_not_closed_certified"
  )
  expect_true(result$evaluated)
  expect_true(result$response_quotient_closure_obstruction_certified)
  expect_true(result$source_exponential_balance_contract_matches)
  expect_true(result$source_exponential_balance_limit_classified)
  expect_true(result$source_bounded_response_image_escape_certified)
  expect_true(result$exact_affine_row_relation_checked)
  expect_true(result$exact_affine_row_relation_certified)
  expect_true(result$zero_target_rows_certified)
  expect_true(result$common_slope_owner_witness_certified)
  expect_true(result$nonzero_target_witness_certified)
  expect_true(result$finite_contrast_target_in_response_image_closure)
  expect_false(result$finite_contrast_target_representable_by_finite_parameter)
  expect_false(result$finite_parameter_response_image_closed)
  expect_false(result$response_metric_quotient_complete)
  expect_false(result$induced_quotient_to_contrast_map_proper)
  expect_false(result$collapsing_finite_response_equivalence_fibres_suffices)
  expect_true(result$quotient_completion_or_boundary_strata_required)
  expect_identical(result$target_cumulative_contrasts, matrix(c(1, 2, 0, 0)))
  expect_identical(result$target_adjacent_contrasts, c(1, 2, 0, 0))
  expect_identical(result$zero_rows, c(3L, 4L))
  expect_identical(result$witness_rows, c(1L, 2L))
  expect_identical(result$witness_slope_owner, 1L)
  expect_identical(result$affine_relation_residual, c(0, 0, 0, 0))
  expect_identical(result$witness_target_linear_combination, -1)
  expect_equal(
    result$boundary_log_likelihood_continuation,
    source$log_likelihood_limit,
    tolerance = 0
  )
  expect_false(result$complete_response_image_closure_stratified)
  expect_false(result$complete_bounded_image_escape_family_classified)
  expect_false(result$complete_escaping_sequence_boundary_envelope_constructed)
  expect_false(result$competitive_boundary_certified)
  expect_false(result$global_supremum_value_classified)
  expect_false(result$finite_jmle_existence_certified)
  expect_false(result$finite_jmle_nonexistence_certified)
  expect_false(result$global_boundary_absence_certified)
  expect_false(result$standard_error_eligible)
  expect_false(result$confidence_interval_eligible)
  expect_false(result$external_comparison_eligible)
  expect_identical(result$readiness_effect, "none_theorem_only")
  expect_identical(
    result$next_gate,
    paste0(
      "complete_response_image_closure_stratification_and_boundary_",
      "limsup_envelope"
    )
  )
  expect_identical(
    result$canonical_fixture_identity,
    "binary_two_person_two_slope_owner_centered_location_gpcm"
  )
})

test_that("the canonical finite response sequence has the declared limit", {
  source <- canonical_response_quotient_source()
  design <- as.matrix(source$path_definition$adjacent_design)
  indices <- c(0, 1, 3, 8)
  contrasts <- vapply(indices, function(index) {
    additive <- c(exp(-index), 2 * exp(-index), 0)
    slopes <- exp(c(index, -index))
    adjacent <- as.numeric(design %*% additive)
    adjacent * slopes[source$path_definition$slope_index]
  }, numeric(4))
  expected <- rbind(
    rep(1, length(indices)),
    rep(2, length(indices)),
    exp(-2 * indices),
    2 * exp(-2 * indices)
  )

  expect_equal(contrasts, expected, tolerance = 1e-15)
  expect_equal(contrasts[, length(indices)], c(1, 2, 0, 0),
               tolerance = 3e-7)
  distances <- apply(
    contrasts - c(1, 2, 0, 0), 2L,
    function(column) sqrt(sum(column^2))
  )
  expect_true(all(diff(distances) < 0))
  expect_true(source$bounded_response_image_escape_certified)
  expect_true(source$parameter_path_divergent)

  finite <- lapply(indices, function(index) {
    mfrmr:::mfrmr_jml_gpcm_exponential_balance_loglik_at(source, index)
  })
  expect_true(all(vapply(finite, `[[`, logical(1), "valid")))
  expect_true(all(is.finite(vapply(
    finite, `[[`, numeric(1), "log_likelihood"
  ))))
})

test_that("the exact affine relation supplies the finite-image contradiction", {
  source <- canonical_response_quotient_source()
  result <- call_response_quotient_obstruction(source)
  design <- as.matrix(source$path_definition$adjacent_design)
  offset <- source$path_definition$additive_offset
  augmented <- cbind(offset, design)

  expect_identical(augmented[1, ] - augmented[2, ],
                   augmented[3, ] - augmented[4, ])
  expect_identical(result$affine_relation_residual, rep(0, 4))

  # If a finite parameter represented target rows 3 and 4 as zero, positivity
  # of their common finite slope would force u3 = u4 = 0.  The exact affine
  # relation then gives u1 - u2 = 0.  Rows 1 and 2 share one positive slope,
  # so their response contrasts must agree, contrary to target 1 versus 2.
  target <- result$target_adjacent_contrasts
  expect_identical(target[3:4], c(0, 0))
  expect_identical(target[1] - target[2], -1)
  expect_true(result$response_quotient_closure_obstruction_certified)

  offset_source <- canonical_response_quotient_source(
    additive_offset = c(3, 3, 4, 4),
    additive_center = c(0, 0, 3)
  )
  offset_result <- call_response_quotient_obstruction(offset_source)
  expect_true(any(offset_source$path_definition$additive_offset != 0))
  expect_true(offset_source$bounded_response_image_escape_certified)
  expect_identical(offset_result$affine_relation_residual, rep(0, 4))
  expect_true(offset_result$response_quotient_closure_obstruction_certified)
})

test_that("quotient closure certificates fail closed", {
  source <- canonical_response_quotient_source()

  bad_contract <- source
  bad_contract$contract_version <- "legacy"
  result <- call_response_quotient_obstruction(bad_contract)
  expect_identical(result$state, "not_evaluated_source")
  expect_false(result$evaluated)
  expect_false(result$response_quotient_closure_obstruction_certified)
  expect_true(is.na(result$finite_parameter_response_image_closed))
  expect_true(is.na(
    result$collapsing_finite_response_equivalence_fibres_suffices
  ))
  expect_identical(
    result$next_gate, "unclassified_response_quotient_closure"
  )

  result <- call_response_quotient_obstruction(1)
  expect_identical(result$state, "not_evaluated_source")
  expect_false(result$source_exponential_balance_contract_matches)

  not_escape <- source
  not_escape$bounded_response_image_escape_certified <- FALSE
  result <- call_response_quotient_obstruction(not_escape)
  expect_identical(result$state, "not_evaluated_source")

  missing_dimension <- source
  missing_dimension$observations <- NULL
  result <- call_response_quotient_obstruction(missing_dimension)
  expect_identical(result$state, "not_evaluated_path")
  expect_false(result$evaluated)

  result <- call_response_quotient_obstruction(source, max_certificate_rows = 0)
  expect_identical(result$state, "not_evaluated_control")
  expect_false(result$evaluated)

  result <- call_response_quotient_obstruction(
    source, zero_rows = c(3L, 3L)
  )
  expect_identical(result$state, "not_evaluated_declaration")

  result <- call_response_quotient_obstruction(
    source, witness_coefficients = 1
  )
  expect_identical(result$state, "not_evaluated_declaration")

  result <- call_response_quotient_obstruction(
    source, max_certificate_elements = 1
  )
  expect_identical(result$state, "not_evaluated_declaration")

  result <- call_response_quotient_obstruction(
    source,
    witness_rows = c(1L, 3L),
    zero_rows = c(2L, 4L)
  )
  expect_identical(result$state, "not_certified_mixed_witness_owner")
  expect_true(result$evaluated)
  expect_false(result$response_quotient_closure_obstruction_certified)

  result <- call_response_quotient_obstruction(
    source,
    zero_rows = c(1L, 2L),
    witness_rows = c(3L, 4L)
  )
  expect_identical(
    result$state, "not_certified_nonzero_declared_zero_target"
  )

  result <- call_response_quotient_obstruction(
    source, witness_coefficients = c(1, 1)
  )
  expect_identical(result$state, "not_certified_affine_relation")
  expect_true(result$exact_affine_row_relation_checked)
  expect_false(result$exact_affine_row_relation_certified)
  expect_true(result$zero_target_rows_certified)
  expect_true(any(result$affine_relation_residual != 0))

  equal_target <- canonical_response_quotient_source(1, 1)
  result <- call_response_quotient_obstruction(equal_target)
  expect_identical(result$state, "not_certified_zero_target_witness")
  expect_identical(result$witness_target_linear_combination, 0)
  expect_true(result$exact_affine_row_relation_certified)
  expect_false(result$nonzero_target_witness_certified)
  expect_false(result$response_quotient_closure_obstruction_certified)
  expect_false(result$finite_jmle_existence_certified)
  expect_false(result$finite_jmle_nonexistence_certified)
  expect_false(result$external_comparison_eligible)
})
