test_that("minimal binary GPCM finite image and closure have exact strata", {
  classify <- mfrmr:::mfrmr_jml_gpcm_binary_closure_image_state

  positive <- classify(c(2, 1, 4, 1))
  expect_true(positive$valid)
  expect_identical(
    positive$state, "finite_positive_difference_interior"
  )
  expect_identical(positive$differences, c(1, 3))
  expect_true(positive$finite_response_image)
  expect_true(positive$finite_response_image_closure)
  expect_false(positive$bounded_response_escape_limit)

  negative <- classify(c(1, 2, 1, 4))
  expect_identical(
    negative$state, "finite_negative_difference_interior"
  )
  expect_true(negative$finite_response_image)
  expect_true(negative$finite_response_image_closure)

  intersection <- classify(c(3, 3, -2, -2))
  expect_identical(
    intersection$state, "finite_equal_difference_intersection"
  )
  expect_true(intersection$finite_response_image)
  expect_true(intersection$bounded_response_escape_limit)

  owner_1_axis <- classify(c(1, 2, 0, 0))
  expect_identical(
    owner_1_axis$state, "missing_owner_1_dominant_axis"
  )
  expect_false(owner_1_axis$finite_response_image)
  expect_true(owner_1_axis$finite_response_image_closure)
  expect_true(owner_1_axis$missing_finite_boundary)
  expect_true(owner_1_axis$bounded_response_escape_limit)

  owner_2_axis <- classify(c(0, 0, 1, 2))
  expect_identical(
    owner_2_axis$state, "missing_owner_2_dominant_axis"
  )
  expect_false(owner_2_axis$finite_response_image)
  expect_true(owner_2_axis$finite_response_image_closure)
  expect_true(owner_2_axis$missing_finite_boundary)

  outside <- classify(c(1, 2, 2, 1))
  expect_identical(
    outside$state, "outside_finite_response_image_closure"
  )
  expect_false(outside$finite_response_image)
  expect_false(outside$finite_response_image_closure)
  expect_false(outside$missing_finite_boundary)
  expect_false(outside$bounded_response_escape_limit)

  invalid <- classify(c(1, 2, 3))
  expect_false(invalid$valid)
  expect_identical(invalid$state, "invalid_contrast_vector")

  strata <- mfrmr:::mfrmr_jml_gpcm_binary_closure_strata()
  expect_identical(nrow(strata), 5L)
  expect_identical(strata$Dimension, c(4L, 4L, 2L, 3L, 3L))
  expect_identical(
    strata$FiniteResponseImage, c(TRUE, TRUE, TRUE, FALSE, FALSE)
  )
  expect_true(all(strata$FiniteClosure))
  expect_identical(
    strata$BoundedResponseEscapeLimit,
    c(FALSE, FALSE, TRUE, TRUE, TRUE)
  )
})

test_that("finite image points have exact inverse parameter representatives", {
  points <- list(
    c(2, 1, 4, 1),
    c(1, 3, -2, -1),
    c(3, 3, -2, -2)
  )
  for (point in points) {
    representative <-
      mfrmr:::mfrmr_jml_gpcm_binary_closure_representative(point)
    expect_true(representative$available)
    expect_equal(prod(representative$alpha), 1, tolerance = 1e-15)
    expect_equal(
      sum(representative$log_alpha), 0, tolerance = 1e-15
    )
    expect_equal(
      representative$reconstructed_contrasts,
      point,
      tolerance = 1e-14
    )
    expect_lte(
      representative$maximum_absolute_reconstruction_residual, 1e-14
    )
  }

  missing <- mfrmr:::mfrmr_jml_gpcm_binary_closure_representative(
    c(1, 2, 0, 0)
  )
  expect_false(missing$available)
  expect_length(missing$alpha, 0L)
})

test_that("positive outcome masses give nonseparated boundary nonattainment", {
  result <-
    mfrmr:::mfrmr_jml_gpcm_binary_closure_nonattainment_fixture()
  expected <- c(log(2), log(4), 0, 0)

  expect_identical(
    result$contract_version,
    "mfrmr-jml-gpcm-binary-closure-envelope-0.2.3-v1"
  )
  expect_identical(
    result$state,
    "finite_jmle_nonexistence_missing_response_boundary_certified"
  )
  expect_true(result$evaluated)
  expect_true(result$binary_closure_envelope_classified)
  expect_true(result$source_response_quotient_contract_matches)
  expect_true(result$strictly_positive_success_mass_per_cell)
  expect_true(result$strictly_positive_failure_mass_per_cell)
  expect_true(result$ordinary_single_cell_separation_absent)
  expect_true(result$integer_expanded_row_multiplicities)
  expect_identical(result$response_contrast_dimension, 4L)
  expect_identical(result$slope_owners, 2L)
  expect_true(result$finite_response_image_exactly_characterized)
  expect_true(result$complete_finite_response_image_closure_stratified)
  expect_true(result$complete_bounded_response_escape_limit_set_classified)
  expect_true(
    result$unbounded_response_contrast_likelihood_limit_negative_infinity
  )
  expect_true(
    result$complete_fixture_parameter_sequence_limsup_envelope_constructed
  )
  expect_false(result$general_gpcm_response_image_closure_classified)
  expect_identical(result$success_mass, c(2, 4, 1, 1))
  expect_identical(result$failure_mass, rep(1, 4))
  expect_equal(
    result$independent_contrast_optimum, expected, tolerance = 0
  )
  expect_identical(
    result$independent_optimum_image_state,
    "missing_owner_1_dominant_axis"
  )
  expect_false(result$independent_optimum_in_finite_response_image)
  expect_true(result$independent_optimum_in_finite_response_image_closure)
  expect_true(result$independent_optimum_on_missing_finite_boundary)
  expect_identical(result$boundary_axis_candidate_count, 2L)
  expect_identical(result$boundary_envelope_winner_count, 1L)
  expect_identical(result$global_closure_maximizer_count, 1L)
  expect_equal(
    as.numeric(result$global_closure_maximizers[1, ]),
    expected,
    tolerance = 0
  )
  expect_equal(
    result$global_supremum_value, -7.184143, tolerance = 1e-6
  )
  expect_identical(
    result$global_supremum_value,
    result$boundary_limsup_supremum_value
  )
  expect_identical(result$finite_over_boundary_log_likelihood_gap, 0)
  expect_false(result$strict_boundary_limsup_gap_certified)
  expect_false(result$finite_jmle_existence_certified)
  expect_true(result$finite_jmle_nonexistence_certified)
  expect_true(
    result$finite_jmle_nonexistence_without_cell_separation_certified
  )
  expect_false(result$response_image_supremum_attained_by_finite_parameter)
  expect_true(result$competitive_boundary_certified)
  expect_true(result$global_supremum_value_classified)
  expect_true(result$complete_escaping_sequence_boundary_envelope_constructed)
  expect_false(result$global_boundary_absence_certified)
  expect_false(result$standard_error_eligible)
  expect_false(result$confidence_interval_eligible)
  expect_false(result$external_comparison_eligible)
  expect_identical(result$readiness_effect, "none_theorem_only")
  expect_identical(
    result$next_gate,
    "general_design_response_image_closure_stratification"
  )

  axes <- result$boundary_axis_candidates
  expect_identical(
    axes$Axis, c("owner_1_dominant", "owner_2_dominant")
  )
  expect_identical(axes$Constraint, c("d2=0", "d1=0"))
  expect_identical(axes$BoundaryEnvelopeWinner, c(TRUE, FALSE))
  expect_identical(axes$GlobalClosureMaximizer, c(TRUE, FALSE))
  expect_identical(axes$MissingFiniteBoundary, c(TRUE, FALSE))
  expect_equal(
    as.numeric(axes[1, c("Z1", "Z2", "Z3", "Z4")]),
    expected,
    tolerance = 0
  )
})

test_that("finite boundary paths converge to the nonattained supremum", {
  result <-
    mfrmr:::mfrmr_jml_gpcm_binary_closure_nonattainment_fixture()
  indices <- c(0, 1, 3, 8)
  paths <- lapply(indices, function(index) {
    mfrmr:::mfrmr_jml_gpcm_binary_closure_path_at(
      result, "owner_1_dominant", index
    )
  })

  expect_true(all(vapply(paths, `[[`, logical(1), "valid")))
  expect_true(all(vapply(
    paths, `[[`, logical(1), "parameter_path_declared_divergent"
  )))
  expect_true(all(vapply(paths, function(path) {
    abs(sum(path$log_alpha)) <= 1e-15 &&
      abs(prod(path$alpha) - 1) <= 1e-15
  }, logical(1))))
  expect_true(all(vapply(paths, function(path) {
    state <-
      mfrmr:::mfrmr_jml_gpcm_binary_closure_image_state(path$contrasts)
    isTRUE(state$finite_response_image)
  }, logical(1))))
  expect_true(all(vapply(paths, function(path) {
    path$maximum_absolute_reconstruction_residual <= 1e-14
  }, logical(1))))

  distances <- vapply(
    paths, `[[`, numeric(1), "contrast_distance_to_target"
  )
  likelihoods <- vapply(paths, `[[`, numeric(1), "log_likelihood")
  gaps <- vapply(
    paths, `[[`, numeric(1), "target_minus_path_log_likelihood"
  )
  expect_true(all(diff(distances) < 0))
  expect_true(all(diff(likelihoods) > 0))
  expect_true(all(gaps >= -1e-14))
  expect_lt(tail(distances, 1), 6e-8)
  expect_lt(tail(gaps, 1), 2e-14)
  expect_equal(
    tail(paths, 1)[[1]]$target_contrasts,
    c(log(2), log(4), 0, 0),
    tolerance = 0
  )
})

test_that("same-sign independent logits give finite global JMLEs", {
  positive <- mfrmr:::mfrmr_jml_gpcm_binary_closure_envelope(
    success_mass = c(4, 2, 4, 2),
    failure_mass = rep(1, 4)
  )
  expect_identical(
    positive$state, "finite_jmle_exists_with_strict_boundary_gap"
  )
  expect_identical(
    positive$independent_optimum_image_state,
    "finite_positive_difference_interior"
  )
  expect_true(positive$finite_jmle_existence_certified)
  expect_false(positive$finite_jmle_nonexistence_certified)
  expect_true(positive$response_image_supremum_attained_by_finite_parameter)
  expect_true(positive$strict_boundary_limsup_gap_certified)
  expect_gt(positive$finite_over_boundary_log_likelihood_gap, 0)
  expect_false(positive$competitive_boundary_certified)
  expect_true(positive$finite_global_representative$available)
  expect_equal(
    positive$finite_global_representative$reconstructed_contrasts,
    positive$independent_contrast_optimum,
    tolerance = 1e-14
  )

  negative <- mfrmr:::mfrmr_jml_gpcm_binary_closure_envelope(
    success_mass = c(2, 4, 2, 4),
    failure_mass = rep(1, 4)
  )
  expect_identical(
    negative$independent_optimum_image_state,
    "finite_negative_difference_interior"
  )
  expect_true(negative$finite_jmle_existence_certified)
  expect_true(negative$strict_boundary_limsup_gap_certified)
  expect_equal(
    negative$global_supremum_value,
    positive$global_supremum_value,
    tolerance = 1e-14
  )
})

test_that("equal logits have a finite maximum and an escaping equivalent fibre", {
  result <- mfrmr:::mfrmr_jml_gpcm_binary_closure_envelope(
    success_mass = rep(2, 4),
    failure_mass = rep(1, 4)
  )
  expect_identical(
    result$state,
    "finite_jmle_exists_on_response_equivalence_intersection"
  )
  expect_identical(
    result$independent_optimum_image_state,
    "finite_equal_difference_intersection"
  )
  expect_true(result$finite_jmle_existence_certified)
  expect_false(result$finite_jmle_nonexistence_certified)
  expect_identical(result$finite_over_boundary_log_likelihood_gap, 0)
  expect_false(result$strict_boundary_limsup_gap_certified)
  expect_true(result$competitive_boundary_certified)
  expect_true(all(result$boundary_axis_candidates$GlobalClosureMaximizer))
  expect_true(result$finite_global_representative$available)

  path <- mfrmr:::mfrmr_jml_gpcm_binary_closure_path_at(
    result, "owner_2_dominant", 8
  )
  expect_true(path$valid)
  expect_identical(path$contrasts, result$independent_contrast_optimum)
  expect_identical(path$contrast_distance_to_target, 0)
  expect_equal(path$log_likelihood, result$global_supremum_value,
               tolerance = 1e-14)
  expect_identical(path$log_alpha, c(-8, 8))
})

test_that("opposite independent differences select missing closure axes", {
  result <- mfrmr:::mfrmr_jml_gpcm_binary_closure_envelope(
    success_mass = c(2, 4, 4, 2),
    failure_mass = rep(1, 4)
  )
  expect_identical(
    result$independent_optimum_image_state,
    "outside_finite_response_image_closure"
  )
  expect_false(result$independent_optimum_in_finite_response_image)
  expect_false(result$independent_optimum_in_finite_response_image_closure)
  expect_false(result$independent_optimum_on_missing_finite_boundary)
  expect_true(result$finite_jmle_nonexistence_certified)
  expect_true(all(
    result$boundary_axis_candidates$MissingFiniteBoundary
  ))
  expect_true(all(
    result$boundary_axis_candidates$BoundaryEnvelopeWinner
  ))
  expect_identical(result$global_closure_maximizer_count, 2L)
  expect_identical(result$finite_over_boundary_log_likelihood_gap, 0)
})

binary_closure_fit_fixture <- function() {
  success <- c(2L, 4L, 1L, 1L)
  failure <- rep(1L, 4L)
  person <- c("P1", "P2", "P1", "P2")
  owner <- c("G1", "G1", "G2", "G2")
  data <- do.call(rbind, lapply(1:4, function(index) {
    data.frame(
      Person = person[index],
      Owner = owner[index],
      Score = c(rep(1L, success[index]), rep(0L, failure[index])),
      stringsAsFactors = FALSE
    )
  }))
  suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = "Owner",
    score = "Score",
    model = "GPCM",
    method = "JML",
    step_facet = "Owner",
    slope_facet = "Owner",
    rating_min = 0,
    rating_max = 1,
    maxit = 200,
    reltol = 1e-10
  ))
}

test_that("current exact GPCM fit reconstructs the P2i nonattainment envelope", {
  skip_if_not_installed("lpSolve")
  fit <- binary_closure_fit_fixture()
  wrapper <- mfrmr:::mfrmr_jml_gpcm_fit_binary_closure_envelope(fit)

  expect_identical(
    wrapper$state,
    paste0(
      "current_fit_",
      "finite_jmle_nonexistence_missing_response_boundary_certified"
    )
  )
  expect_true(wrapper$evaluated)
  expect_true(wrapper$current_fit_reconstructed)
  expect_true(wrapper$exact_canonical_affine_operator_certified)
  expect_true(wrapper$complete_two_person_by_two_owner_cell_support)
  expect_identical(
    wrapper$response_cell_success_mass, c(2, 4, 1, 1)
  )
  expect_identical(wrapper$response_cell_failure_mass, rep(1, 4))
  expect_identical(
    wrapper$canonical_adjacent_operator_rows,
    rbind(
      c(1, 0, -1), c(0, 1, -1),
      c(1, 0, 1), c(0, 1, 1)
    )
  )
  expect_true(is.finite(wrapper$retained_log_likelihood))
  expect_length(wrapper$retained_response_cell_contrasts, 4L)
  expect_true(is.finite(wrapper$reconstructed_retained_log_likelihood))
  expect_lte(abs(wrapper$retained_objective_reconstruction_residual), 1e-10)
  expect_true(wrapper$current_fit_objective_numerically_reconstructed)
  expect_equal(
    wrapper$classified_global_supremum, -7.184143, tolerance = 1e-6
  )
  expect_gt(wrapper$supremum_minus_retained_log_likelihood, 0)
  expect_lt(wrapper$supremum_minus_retained_log_likelihood, 1e-3)
  expect_false(wrapper$finite_jmle_existence_certified)
  expect_true(wrapper$finite_jmle_nonexistence_certified)
  expect_true(
    wrapper$finite_jmle_nonexistence_without_cell_separation_certified
  )
  expect_true(
    wrapper$complete_escaping_sequence_boundary_envelope_constructed
  )
  expect_false(wrapper$optimizer_trace_promoted_to_finite_estimate)
  expect_false(wrapper$standard_error_eligible)
  expect_false(wrapper$confidence_interval_eligible)
  expect_false(wrapper$external_comparison_eligible)
  expect_identical(wrapper$readiness_effect, "none_theorem_only")

  expect_identical(
    fit$config$boundary_audit$gpcm_global_existence$state,
    "global_finite_jmle_existence_open_complete_boundary_envelope_unavailable"
  )
  expect_false(
    fit$config$boundary_audit$gpcm_global_existence$
      finite_jmle_nonexistence_certified
  )
  expect_true(wrapper$result$finite_jmle_nonexistence_certified)
})

test_that("current-fit P2i wrapper keeps nonmatching routes fail closed", {
  skip_if_not_installed("lpSolve")
  fit <- binary_closure_fit_fixture()

  invalid <- mfrmr:::mfrmr_jml_gpcm_fit_binary_closure_envelope(1)
  expect_identical(invalid$state, "not_evaluated_fit")
  expect_false(invalid$evaluated)

  pcm <- fit
  pcm$config$model <- "PCM"
  result <- mfrmr:::mfrmr_jml_gpcm_fit_binary_closure_envelope(pcm)
  expect_identical(result$state, "not_applicable_model")
  expect_false(result$evaluated)

  mml <- fit
  mml$config$method <- "MML"
  result <- mfrmr:::mfrmr_jml_gpcm_fit_binary_closure_envelope(mml)
  expect_identical(result$state, "not_applicable_estimator")
  expect_false(result$evaluated)

  unsupported <- fit
  unsupported$prep$facet_names <- c("Owner", "SecondFacet")
  result <-
    mfrmr:::mfrmr_jml_gpcm_fit_binary_closure_envelope(unsupported)
  expect_false(result$finite_jmle_existence_certified)
  expect_false(result$finite_jmle_nonexistence_certified)

  legacy <-
    mfrmr:::mfrmr_jml_gpcm_response_quotient_binary_counterexample()
  legacy$contract_version <- "legacy"
  result <- mfrmr:::mfrmr_jml_gpcm_fit_binary_closure_envelope(
    fit, source_response_quotient = legacy
  )
  expect_identical(result$state, "not_evaluated_response_envelope")
  expect_true(result$evaluated)
  expect_true(result$current_fit_reconstructed)
  expect_false(result$finite_jmle_existence_certified)
  expect_false(result$finite_jmle_nonexistence_certified)
  expect_false(result$external_comparison_eligible)
})

test_that("binary closure envelope declarations fail closed", {
  valid <-
    mfrmr:::mfrmr_jml_gpcm_response_quotient_binary_counterexample()
  legacy <- valid
  legacy$contract_version <- "legacy"

  result <- mfrmr:::mfrmr_jml_gpcm_binary_closure_envelope(
    rep(1, 4), rep(1, 4), source_response_quotient = legacy
  )
  expect_identical(result$state, "not_evaluated_source")
  expect_false(result$evaluated)
  expect_false(result$binary_closure_envelope_classified)
  expect_false(result$finite_jmle_existence_certified)
  expect_false(result$finite_jmle_nonexistence_certified)
  expect_false(result$complete_escaping_sequence_boundary_envelope_constructed)
  expect_identical(
    result$next_gate, "unclassified_minimal_binary_closure_envelope"
  )

  cases <- list(
    list(success_mass = c(1, 1, 1), failure_mass = rep(1, 4)),
    list(success_mass = c(1, 1, 1, 0), failure_mass = rep(1, 4)),
    list(success_mass = rep(1, 4), failure_mass = c(1, 1, 1, NA)),
    list(success_mass = rep(1, 4), failure_mass = rep(1, 4),
         max_total_mass = 7)
  )
  for (arguments in cases) {
    result <- do.call(
      mfrmr:::mfrmr_jml_gpcm_binary_closure_envelope,
      arguments
    )
    expect_identical(result$state, "not_evaluated_response_mass")
    expect_false(result$evaluated)
    expect_false(result$global_supremum_value_classified)
    expect_false(result$finite_jmle_existence_certified)
    expect_false(result$finite_jmle_nonexistence_certified)
  }

  result <- mfrmr:::mfrmr_jml_gpcm_binary_closure_envelope(
    rep(1, 4), rep(1, 4), max_total_mass = 0
  )
  expect_identical(result$state, "not_evaluated_control")
  expect_false(result$evaluated)

  classified <-
    mfrmr:::mfrmr_jml_gpcm_binary_closure_nonattainment_fixture()
  invalid_contract <- classified
  invalid_contract$contract_version <- "legacy"
  expect_false(mfrmr:::mfrmr_jml_gpcm_binary_closure_path_at(
    invalid_contract, "owner_1_dominant", 1
  )$valid)
  expect_false(mfrmr:::mfrmr_jml_gpcm_binary_closure_path_at(
    classified, "unknown", 1
  )$valid)
  expect_false(mfrmr:::mfrmr_jml_gpcm_binary_closure_path_at(
    classified, "owner_1_dominant", -1
  )$valid)
  expect_false(mfrmr:::mfrmr_jml_gpcm_binary_closure_path_at(
    classified, "owner_1_dominant", 301
  )$valid)
})
