fit_jml_gpcm_parameter_path_fixture <- function() {
  data <- expand.grid(
    Person = c("P1", "P2"),
    Criterion = c("C1", "C2", "C3"),
    Score = 0:1,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data <- data[rep(seq_len(nrow(data)), each = 2L), , drop = FALSE]
  suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = "Criterion",
    score = "Score",
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    rating_min = 0,
    rating_max = 1,
    maxit = 200,
    reltol = 1e-12
  ))
}

parameter_path_design <- function() {
  matrix(
    c(
      1, 0,
      0, 1,
      1, 1,
      -1, 1
    ),
    nrow = 4L,
    byrow = TRUE
  )
}

parameter_path_map <- function() {
  data.frame(
    OptimizerIndex = c(3L, 7L),
    Block = c("theta", "steps"),
    Coordinate = c("theta:P1", "step:C1:1"),
    Facet = c("Person", "Criterion"),
    Level = c("P1", "C1"),
    stringsAsFactors = FALSE
  )
}

test_that("free-coordinate stages map exactly to cumulative utility paths", {
  mapped <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    adjacent_design = parameter_path_design(),
    coordinate_stage_directions = list(c(2, -1), c(-1, 0.5)),
    n_observations = 2L,
    n_categories = 3L,
    parameter_map = parameter_path_map()
  )

  expect_identical(
    mapped$contract_version,
    "mfrmr-jml-gpcm-parameter-path-reachability-0.2.3-v1"
  )
  expect_identical(mapped$state, "reachable_additive_utility_path")
  expect_true(mapped$evaluated)
  expect_true(mapped$parameter_space_reachability_checked)
  expect_true(mapped$parameter_space_reachability_certified)
  expect_identical(
    mapped$reachability_basis,
    "forward_map_from_retained_constrained_free_additive_coordinates"
  )
  expect_false(mapped$inverse_projection_used)
  expect_false(mapped$inverse_tolerance_reachability_claim)
  expect_true(mapped$exact_tie_contract_preserved_by_forward_construction)
  expect_identical(mapped$retained_observations, 2L)
  expect_identical(mapped$categories, 3L)
  expect_identical(mapped$adjacent_design_rows, 4L)
  expect_identical(mapped$additive_free_coordinates, 2L)
  expect_identical(mapped$stages, 2L)
  expect_identical(mapped$maximum_stages, 2L)
  expect_identical(
    mapped$category_zero_gauge,
    "first_category_utility_fixed_to_zero"
  )
  expect_false(mapped$probability_equivalent_row_shifts_required)
  expect_false(mapped$path_inferred_from_optimizer_trace)
  expect_false(mapped$scale_exponents_inferred)
  expect_false(mapped$path_search_performed)
  expect_false(mapped$arbitrary_utility_direction_reachability_classified)
  expect_false(mapped$external_comparison_eligible)
  expect_identical(mapped$readiness_effect, "none_diagnostic_only")

  expect_equal(
    mapped$adjacent_stage_directions[[1]],
    matrix(c(2, -1, 1, -3), nrow = 2L)
  )
  expect_equal(
    mapped$category_stage_directions[[1]],
    rbind(c(0, 2, 3), c(0, -1, -4))
  )
  expect_equal(
    mapped$category_stage_directions[[2]],
    rbind(c(0, -1, -1.5), c(0, 0.5, 2))
  )
  expect_identical(
    mapped$stage_table$ReachabilityState,
    rep("reachable_by_forward_construction", 2L)
  )
  expect_identical(mapped$stage_table$NonzeroCategoryContrasts, c(4L, 4L))
  expect_identical(mapped$loading_table$OptimizerIndex, rep(c(3L, 7L), 2L))
  expect_equal(mapped$loading_table$Loading, c(2, -1, -1, 0.5))
})

test_that("invertible free-coordinate changes preserve the utility path", {
  design <- parameter_path_design()
  transform <- matrix(c(2, 1, -1, 1), nrow = 2L)
  direction <- c(2, -1)
  transformed_direction <- solve(transform, direction)

  original <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    design, direction, 2L, 3L
  )
  transformed <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    design %*% transform, transformed_direction, 2L, 3L
  )

  expect_true(original$parameter_space_reachability_certified)
  expect_true(transformed$parameter_space_reachability_certified)
  expect_equal(
    transformed$adjacent_stage_directions[[1]],
    original$adjacent_stage_directions[[1]],
    tolerance = 1e-14
  )
  expect_equal(
    transformed$category_stage_directions[[1]],
    original$category_stage_directions[[1]],
    tolerance = 1e-14
  )
  expect_identical(
    transformed$category_stage_directions[[1]][, 1],
    c(0, 0)
  )
})

test_that("empty additive paths are reachable without manufacturing a stage", {
  mapped <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    parameter_path_design(), list(), 2L, 3L,
    parameter_map = parameter_path_map()
  )

  expect_identical(
    mapped$state, "no_additive_stage_reachability_vacuous"
  )
  expect_true(mapped$evaluated)
  expect_true(mapped$parameter_space_reachability_certified)
  expect_identical(mapped$stages, 0L)
  expect_length(mapped$coordinate_stage_directions, 0L)
  expect_length(mapped$adjacent_stage_directions, 0L)
  expect_length(mapped$category_stage_directions, 0L)
  expect_identical(mapped$parameter_map$OptimizerIndex, c(3L, 7L))
})

test_that("forward reachability failures remain typed and fail closed", {
  design <- parameter_path_design()

  invalid_design <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    design[-1, ], c(1, 0), 2L, 3L
  )
  expect_identical(invalid_design$state, "not_evaluated_design")
  expect_false(invalid_design$parameter_space_reachability_certified)

  invalid_direction <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    design, c(1, 0, 0), 2L, 3L
  )
  expect_identical(
    invalid_direction$state, "not_evaluated_coordinate_directions"
  )

  nonfinite_direction <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    design, c(1, Inf), 2L, 3L
  )
  expect_identical(
    nonfinite_direction$state, "not_evaluated_coordinate_directions"
  )

  too_many <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    design,
    rbind(c(1, 0), c(0, 1), c(1, -1)),
    2L, 3L
  )
  expect_identical(too_many$state, "not_evaluated_stage_limit")

  null_stage <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    design, c(0, 0), 2L, 3L
  )
  expect_identical(null_stage$state, "not_evaluated_null_utility_stage")
  expect_false(null_stage$evaluated)

  size_limited <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    design, c(1, 0), 2L, 3L, max_design_nonzeros = 0L
  )
  expect_identical(size_limited$state, "not_evaluated_size_limit")

  duplicate_map <- parameter_path_map()
  duplicate_map$OptimizerIndex <- 1L
  invalid_map <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    design, c(1, 0), 2L, 3L, parameter_map = duplicate_map
  )
  expect_identical(invalid_map$state, "not_evaluated_parameter_map")

  invalid_control <- mfrmr:::mfrmr_jml_gpcm_additive_utility_path(
    design, c(1, 0), list("bad"), 3L
  )
  expect_identical(invalid_control$state, "not_evaluated_control")
})

test_that("a current fit maps declared parameter paths into the P2b oracle", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_parameter_path_fixture()
  scope <- fit$config$boundary_audit$gpcm_parameter_path_reachability
  direction <- rep(0, scope$additive_free_coordinates)
  direction[1] <- 1

  result <- mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(
    fit = fit,
    additive_coordinate_directions = list(direction),
    additive_scale_exponents = 0.5,
    slope_stage_rates = c(1, -1, 0),
    slope_scale_exponents = 1
  )

  expect_identical(
    result$contract_version,
    "mfrmr-jml-gpcm-parameter-path-reachability-0.2.3-v1"
  )
  expect_true(result$evaluated)
  expect_true(result$retained_fit_base_reconstructed)
  expect_true(result$declared_parameter_path_reachability_checked)
  expect_true(result$declared_parameter_path_reachability_certified)
  expect_true(result$parameter_reachable_limit_contract_completed)
  expect_true(result$declared_path_limit_classified)
  expect_true(startsWith(result$state, "parameter_reachable_"))
  expect_true(result$additive_path$parameter_space_reachability_certified)
  expect_identical(result$additive_path$stages, 1L)
  expect_true(any(
    result$additive_path$stage_table$NonzeroCategoryContrasts > 0L
  ))
  expect_true(
    result$limit$declared_utility_path_parameter_reachability_checked
  )
  expect_true(
    result$limit$declared_utility_path_parameter_reachability_certified
  )
  expect_identical(
    result$limit$parameter_reachability_contract,
    result$contract_version
  )
  expect_false(result$path_inferred_from_optimizer_trace)
  expect_false(result$scale_exponents_inferred)
  expect_false(result$path_search_performed)
  expect_false(result$monotone_tail_certified)
  expect_false(result$competitive_boundary_certified)
  expect_false(result$common_subsequence_limit_certified)
  expect_false(result$arbitrary_path_limit_classified)
  expect_false(result$global_boundary_classified)
  expect_false(result$global_finite_maximum_certified)
  expect_false(result$global_boundary_absence_certified)
  expect_false(result$standard_error_eligible)
  expect_false(result$confidence_interval_eligible)
  expect_false(result$external_comparison_eligible)
  expect_identical(result$readiness_effect, "none_diagnostic_only")

  direct <- mfrmr:::mfrmr_jml_gpcm_parameter_hierarchy_loglik_at(
    result, 4
  )
  expect_true(direct$valid)
  expect_true(is.finite(direct$log_likelihood))
  expect_length(direct$row_log_probability, nrow(fit$prep$data))

  invalid_scale <- mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(
    fit = fit,
    additive_coordinate_directions = list(direction),
    additive_scale_exponents = numeric(0),
    slope_stage_rates = c(1, -1, 0),
    slope_scale_exponents = 1
  )
  expect_identical(
    invalid_scale$state, "not_evaluated_lexicographic_limit"
  )
  expect_true(invalid_scale$retained_fit_base_reconstructed)
  expect_true(invalid_scale$declared_parameter_path_reachability_checked)
  expect_true(invalid_scale$declared_parameter_path_reachability_certified)
  expect_false(invalid_scale$parameter_reachable_limit_contract_completed)
  expect_false(invalid_scale$declared_path_limit_classified)
})

test_that("fit-level scope records an operator but declares no path", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_parameter_path_fixture()
  audit <- fit$config$boundary_audit$gpcm_parameter_path_reachability

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-parameter-path-reachability-0.2.3-v1"
  )
  expect_identical(
    audit$state,
    "retained_additive_parameter_operator_available_no_path_declared"
  )
  expect_true(audit$evaluated)
  expect_true(audit$declared_free_coordinate_path_wrapper_available)
  expect_true(audit$forward_parameter_reachability_operator_available)
  expect_gt(audit$additive_design_rows, 0L)
  expect_gt(audit$additive_free_coordinates, 0L)
  expect_gt(audit$additive_design_nonzeros, 0)
  expect_identical(
    nrow(audit$coordinate_map), audit$additive_free_coordinates
  )
  expect_identical(audit$declared_paths_evaluated, 0L)
  expect_false(audit$production_parameter_path_reachability_checked)
  expect_false(audit$production_path_limit_classified)
  expect_true(audit$source_lexicographic_limit_contract_matches)
  expect_false(audit$inverse_projection_used)
  expect_false(audit$path_inferred_from_optimizer_trace)
  expect_false(audit$scale_exponents_inferred)
  expect_false(audit$path_search_performed)
  expect_false(audit$arbitrary_path_limit_classified)
  expect_false(audit$common_subsequence_limit_certified)
  expect_false(audit$monotone_tail_certified)
  expect_false(audit$competitive_boundary_certified)
  expect_false(audit$global_boundary_classified)
  expect_false(audit$global_finite_maximum_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
})

test_that("fit wrappers keep model, estimator, unit, and source scopes separate", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_parameter_path_fixture()

  invalid <- unclass(fit)
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(invalid)$state,
    "not_evaluated_fit"
  )

  pcm <- fit
  pcm$config$model <- "PCM"
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(pcm)$state,
    "not_applicable_model"
  )

  mml <- fit
  mml$config$method <- "MML"
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(mml)$state,
    "not_applicable_estimator"
  )

  unit <- fit
  unit$config$gpcm_spec$n_params <- 0L
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(unit)$state,
    "not_required_unit_slope"
  )

  legacy <- fit
  legacy$config$boundary_audit$gpcm_lexicographic_limit$contract_version <-
    "legacy"
  expect_identical(
    mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(legacy)$state,
    "not_evaluated_source_contract"
  )

  bad_result <- mfrmr:::mfrmr_jml_gpcm_parameter_hierarchy_loglik_at(
    list(), 2
  )
  expect_false(bad_result$valid)
  expect_identical(bad_result$reason_codes, "parameter_path_result_invalid")
})
