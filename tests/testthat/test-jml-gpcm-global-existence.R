fit_jml_gpcm_global_existence_fixture <- function(mode = c(
    "crossed", "all_extreme", "one_extreme")) {
  mode <- match.arg(mode)
  data <- expand.grid(
    Person = c("P1", "P2"),
    Criterion = c("C1", "C2", "C3"),
    Rep = 1:2,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$Score <- switch(
    mode,
    crossed = as.integer(
      (match(data$Criterion, c("C1", "C2", "C3")) + data$Rep) %% 2
    ),
    all_extreme = ifelse(data$Person == "P1", 0L, 1L),
    one_extreme = ifelse(
      data$Person == "P1", 0L,
      as.integer((match(data$Criterion, c("C1", "C2", "C3")) +
                    data$Rep) %% 2)
    )
  )
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
    maxit = 160,
    reltol = 1e-11
  ))
}

global_existence_person_direction <- function(fit, levels, signs) {
  scope <- fit$config$boundary_audit$gpcm_parameter_path_reachability
  direction <- rep(0, scope$additive_free_coordinates)
  for (i in seq_along(levels)) {
    keep <- scope$coordinate_map$Block == "Person" &
      scope$coordinate_map$Level == levels[i]
    direction[keep] <- signs[i]
  }
  direction
}

global_existence_source_fixture <- function() {
  list(
    config = list(
      model = "GPCM", method = "JML",
      gpcm_spec = list(n_params = 2L)
    ),
    person = list(
      contract_version = mfrmr:::mfrmr_boundary_contract_version(),
      method = "JML", model = "GPCM",
      parameter_status = data.frame(
        ParameterStatus = c("estimable", "estimable"),
        stringsAsFactors = FALSE
      )
    ),
    additive = list(
      contract_version = mfrmr:::mfrmr_boundary_contract_version(),
      method = "JML", model = "GPCM",
      prescreen = list(cone_certified = FALSE),
      cone_certificate = data.frame(
        Certified = FALSE,
        StrictContrastRows = 0L,
        DirectionL1 = 0,
        stringsAsFactors = FALSE
      )
    ),
    sequence = list(
      contract_version =
        "mfrmr-jml-gpcm-parameter-sequence-flag-0.2.3-v1",
      every_parameter_sequence_has_classifiable_further_subsequence_likelihood_limit =
        TRUE
    ),
    slope = list(
      contract_version = mfrmr:::mfrmr_boundary_contract_version(),
      method = "JML", model = "GPCM", certificates = data.frame()
    ),
    joint = list(
      contract_version = mfrmr:::mfrmr_boundary_contract_version(),
      method = "JML", model = "GPCM",
      certificates = data.frame(), rate_certificates = data.frame()
    )
  )
}

call_global_existence_audit <- function(source, ...) {
  mfrmr:::audit_mfrm_jml_gpcm_global_existence(
    config = source$config,
    person_boundary = source$person,
    joint_additive = source$additive,
    parameter_sequence_flag = source$sequence,
    slope_boundary = source$slope,
    joint_slope_boundary = source$joint,
    retained_log_likelihood = -10,
    ...
  )
}

test_that("a reachable divergent zero-limit path proves finite JMLE nonexistence", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_global_existence_fixture("all_extreme")
  direction <- global_existence_person_direction(
    fit, c("P1", "P2"), c(-1, 1)
  )
  path <- mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(
    fit,
    additive_coordinate_directions = list(direction),
    additive_scale_exponents = 1,
    slope_stage_rates = NULL,
    slope_scale_exponents = numeric(0)
  )
  adjudicated <- mfrmr:::mfrmr_jml_gpcm_reachable_supremum_zero(path)

  expect_true(path$declared_parameter_path_reachability_certified)
  expect_identical(path$log_likelihood_limit, 0)
  expect_identical(
    adjudicated$contract_version,
    "mfrmr-jml-gpcm-global-existence-0.2.3-v1"
  )
  expect_identical(
    adjudicated$state,
    "global_supremum_zero_finite_jmle_nonexistence_certified"
  )
  expect_true(adjudicated$evaluated)
  expect_true(adjudicated$source_parameter_path_contract_matches)
  expect_true(adjudicated$parameter_path_reachability_certified)
  expect_true(adjudicated$parameter_path_limit_classified)
  expect_identical(adjudicated$divergent_parameter_stages, 1L)
  expect_true(adjudicated$path_escapes_every_bounded_parameter_set)
  expect_identical(adjudicated$universal_finite_log_likelihood_upper_bound, 0)
  expect_true(adjudicated$finite_parameter_log_likelihood_strictly_below_zero)
  expect_identical(adjudicated$path_log_likelihood_limit, 0)
  expect_true(adjudicated$path_reaches_universal_upper_bound)
  expect_true(adjudicated$global_boundary_competitiveness_certified)
  expect_true(adjudicated$global_supremum_value_classified)
  expect_identical(adjudicated$global_log_likelihood_supremum, 0)
  expect_true(adjudicated$global_supremum_attained_on_declared_boundary)
  expect_false(adjudicated$global_supremum_attained_at_finite_parameter)
  expect_true(adjudicated$finite_jmle_existence_classified)
  expect_false(adjudicated$finite_jmle_existence_certified)
  expect_true(adjudicated$finite_jmle_nonexistence_certified)
  expect_true(adjudicated$global_boundary_classified)
  expect_true(adjudicated$global_boundary_presence_certified)
  expect_false(adjudicated$global_boundary_absence_certified)
  expect_false(adjudicated$standard_error_eligible)
  expect_false(adjudicated$confidence_interval_eligible)
  expect_false(adjudicated$external_comparison_eligible)
  expect_identical(adjudicated$readiness_effect, "none_diagnostic_only")

  finite <- vapply(c(1, 2, 4, 8), function(distance) {
    mfrmr:::mfrmr_jml_gpcm_parameter_hierarchy_loglik_at(
      path, distance
    )$log_likelihood
  }, numeric(1))
  expect_true(all(is.finite(finite)))
  expect_true(all(finite < 0))
  expect_gt(min(diff(finite)), 0)
  expect_lt(abs(tail(finite, 1)), abs(finite[1]))
})

test_that("a reachable boundary below zero does not classify the global supremum", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_global_existence_fixture("one_extreme")
  direction <- global_existence_person_direction(fit, "P1", -1)
  path <- mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(
    fit,
    additive_coordinate_directions = list(direction),
    additive_scale_exponents = 1,
    slope_stage_rates = NULL,
    slope_scale_exponents = numeric(0)
  )
  adjudicated <- mfrmr:::mfrmr_jml_gpcm_reachable_supremum_zero(path)

  expect_true(path$declared_parameter_path_reachability_certified)
  expect_true(is.finite(path$log_likelihood_limit))
  expect_lt(path$log_likelihood_limit, 0)
  expect_identical(
    adjudicated$state, "reachable_boundary_below_universal_supremum"
  )
  expect_true(adjudicated$evaluated)
  expect_false(adjudicated$path_reaches_universal_upper_bound)
  expect_false(adjudicated$global_supremum_value_classified)
  expect_true(is.na(adjudicated$global_log_likelihood_supremum))
  expect_false(adjudicated$finite_jmle_existence_classified)
  expect_false(adjudicated$finite_jmle_nonexistence_certified)
  expect_false(adjudicated$global_boundary_classified)
})

test_that("bounded and malformed P2c paths fail closed", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_global_existence_fixture("crossed")
  bounded <- mfrmr:::mfrmr_jml_gpcm_fit_parameter_hierarchy_limit(
    fit,
    additive_coordinate_directions = list(),
    additive_scale_exponents = numeric(0),
    slope_stage_rates = NULL,
    slope_scale_exponents = numeric(0)
  )
  bounded_result <- mfrmr:::mfrmr_jml_gpcm_reachable_supremum_zero(bounded)
  expect_identical(
    bounded_result$state, "bounded_parameter_path_not_global_boundary"
  )
  expect_true(bounded_result$evaluated)
  expect_identical(bounded_result$divergent_parameter_stages, 0L)
  expect_false(bounded_result$path_escapes_every_bounded_parameter_set)
  expect_false(bounded_result$finite_jmle_nonexistence_certified)

  malformed <- mfrmr:::mfrmr_jml_gpcm_reachable_supremum_zero(list())
  expect_identical(malformed$state, "not_evaluated_source_contract")
  expect_false(malformed$evaluated)
  expect_false(malformed$source_parameter_path_contract_matches)
  atomic <- mfrmr:::mfrmr_jml_gpcm_reachable_supremum_zero(1)
  expect_identical(atomic$state, "not_evaluated_source_contract")
  expect_false(atomic$evaluated)

  incomplete <- bounded
  incomplete$declared_parameter_path_reachability_certified <- FALSE
  incomplete_result <-
    mfrmr:::mfrmr_jml_gpcm_reachable_supremum_zero(incomplete)
  expect_identical(incomplete_result$state, "not_evaluated_parameter_path")
  expect_false(incomplete_result$evaluated)
  expect_false(incomplete_result$finite_jmle_nonexistence_certified)
})

test_that("the finite-attainment implication requires a strict complete-boundary gap", {
  conditional <- mfrmr:::mfrmr_jml_gpcm_finite_attainment_gap_theorem(
    finite_log_likelihood = -5,
    boundary_limsup_upper_bound = -6,
    complete_boundary_envelope_declared = TRUE
  )
  expect_identical(
    conditional$state,
    "conditional_finite_attainment_from_declared_complete_boundary_gap"
  )
  expect_true(conditional$evaluated)
  expect_true(conditional$strict_finite_over_boundary_gap)
  expect_true(conditional$complete_boundary_envelope_declared_by_caller)
  expect_false(conditional$complete_boundary_envelope_independently_certified)
  expect_true(conditional$objective_continuity_on_finite_parameter_space)
  expect_true(conditional$closed_finite_parameter_space)
  expect_true(conditional$compact_upper_level_set_follows_conditionally)
  expect_true(conditional$finite_attainment_follows_conditionally)
  expect_false(conditional$global_finite_jmle_existence_certified)
  expect_false(conditional$global_boundary_absence_certified)
  expect_false(conditional$production_decision_eligible)
  expect_false(conditional$external_comparison_eligible)
  expect_identical(conditional$readiness_effect, "none_diagnostic_only")

  incomplete <- mfrmr:::mfrmr_jml_gpcm_finite_attainment_gap_theorem(
    -5, -6, FALSE
  )
  expect_identical(
    incomplete$state, "boundary_envelope_incomplete_global_existence_open"
  )
  expect_true(incomplete$strict_finite_over_boundary_gap)
  expect_false(incomplete$finite_attainment_follows_conditionally)

  no_gap <- mfrmr:::mfrmr_jml_gpcm_finite_attainment_gap_theorem(
    -6, -6, TRUE
  )
  expect_identical(
    no_gap$state, "no_strict_boundary_gap_global_existence_open"
  )
  expect_false(no_gap$strict_finite_over_boundary_gap)
  expect_false(no_gap$finite_attainment_follows_conditionally)

  invalid <- mfrmr:::mfrmr_jml_gpcm_finite_attainment_gap_theorem(
    Inf, -6, TRUE
  )
  expect_identical(invalid$state, "not_evaluated_boundary_gap_inputs")
  expect_false(invalid$evaluated)
  expect_false(invalid$global_finite_jmle_existence_certified)
})

test_that("production free extreme Persons certify finite JMLE nonexistence", {
  source <- global_existence_source_fixture()
  source$person$parameter_status$ParameterStatus[1] <- "unbounded_low"
  audit <- call_global_existence_audit(source)

  expect_identical(
    audit$state,
    "finite_jmle_nonexistence_certified_strict_additive_recession"
  )
  expect_true(audit$evaluated)
  expect_true(audit$source_parameter_sequence_flag_contract_matches)
  expect_true(
    audit$every_parameter_sequence_has_classifiable_further_subsequence_likelihood_limit
  )
  expect_identical(audit$free_extreme_person_recession_count, 1L)
  expect_false(audit$global_additive_recession_cone_certified)
  expect_true(audit$global_boundary_competitiveness_certified)
  expect_false(audit$global_supremum_value_classified)
  expect_true(is.na(audit$global_log_likelihood_supremum))
  expect_true(audit$finite_jmle_existence_classified)
  expect_false(audit$finite_jmle_existence_certified)
  expect_true(audit$finite_jmle_nonexistence_certified)
  expect_true(audit$global_boundary_classified)
  expect_true(audit$global_boundary_presence_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$complete_boundary_envelope_constructed)
  expect_false(audit$finite_over_complete_boundary_gap_certified)
  expect_false(audit$production_parameter_path_declared)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
})

test_that("a certified global additive cone proves nonexistence", {
  source <- global_existence_source_fixture()
  source$additive$cone_certificate <- data.frame(
    Certified = TRUE,
    StrictContrastRows = 3L,
    DirectionL1 = 2,
    stringsAsFactors = FALSE
  )
  source$additive$prescreen$cone_certified <- TRUE
  audit <- call_global_existence_audit(source)

  expect_identical(
    audit$state,
    "finite_jmle_nonexistence_certified_strict_additive_recession"
  )
  expect_identical(audit$free_extreme_person_recession_count, 0L)
  expect_true(audit$global_additive_recession_cone_certified)
  expect_true(audit$finite_jmle_nonexistence_certified)
  expect_match(audit$reason_codes, "global_additive_strict_recession_cone")

  source$additive$cone_certificate$StrictContrastRows <- 0L
  zero_margin <- call_global_existence_audit(source)
  expect_false(zero_margin$global_additive_recession_cone_certified)
  expect_false(zero_margin$finite_jmle_nonexistence_certified)
})

test_that("an analytic zero boundary classifies the global supremum", {
  source <- global_existence_source_fixture()
  source$slope$certificates <- data.frame(
    BoundaryLogLikelihood = c(-2, 0),
    Certified = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  audit <- call_global_existence_audit(source)

  expect_identical(
    audit$state,
    "global_supremum_zero_finite_jmle_nonexistence_certified"
  )
  expect_identical(
    audit$universal_upper_bound_zero_path_source,
    "slope_only_analytic_path"
  )
  expect_true(audit$global_supremum_value_classified)
  expect_identical(audit$global_log_likelihood_supremum, 0)
  expect_true(audit$finite_jmle_nonexistence_certified)

  source <- global_existence_source_fixture()
  source$joint$rate_certificates <- data.frame(
    BoundaryLogLikelihood = 0,
    Certified = TRUE,
    AnalyticTailCertified = TRUE,
    stringsAsFactors = FALSE
  )
  joint <- call_global_existence_audit(source)
  expect_identical(
    joint$universal_upper_bound_zero_path_source,
    "joint_slope_additive_analytic_path"
  )
  expect_true(joint$global_supremum_value_classified)

  source$joint$rate_certificates$AnalyticTailCertified <- FALSE
  unverified <- call_global_existence_audit(source)
  expect_false(unverified$global_supremum_value_classified)
})

test_that("absence of a positive certificate leaves finite existence open", {
  source <- global_existence_source_fixture()
  audit <- call_global_existence_audit(source)

  expect_identical(
    audit$state,
    "global_finite_jmle_existence_open_complete_boundary_envelope_unavailable"
  )
  expect_true(audit$evaluated)
  expect_true(audit$retained_point_finite)
  expect_true(audit$retained_point_is_optimizer_trace_only)
  expect_false(audit$global_boundary_competitiveness_certified)
  expect_false(audit$finite_jmle_existence_classified)
  expect_false(audit$finite_jmle_existence_certified)
  expect_false(audit$finite_jmle_nonexistence_certified)
  expect_false(audit$global_boundary_classified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$complete_boundary_envelope_constructed)
  expect_match(audit$reason_codes, "complete_boundary_envelope_unavailable")
})

test_that("current fits retain typed P2f states without readiness promotion", {
  skip_if_not_installed("lpSolve")
  crossed <- fit_jml_gpcm_global_existence_fixture("crossed")
  open <- crossed$config$boundary_audit$gpcm_global_existence
  expect_identical(
    open$contract_version,
    "mfrmr-jml-gpcm-global-existence-0.2.3-v1"
  )
  expect_true(open$evaluated)
  expect_true(open$source_parameter_sequence_flag_contract_matches)
  expect_false(open$finite_jmle_existence_certified)
  expect_false(open$global_boundary_absence_certified)
  expect_false(open$external_comparison_eligible)
  expect_identical(open$readiness_effect, "none_diagnostic_only")

  extreme <- fit_jml_gpcm_global_existence_fixture("all_extreme")
  absent <- extreme$config$boundary_audit$gpcm_global_existence
  expect_true(absent$evaluated)
  expect_true(absent$finite_jmle_nonexistence_certified)
  expect_true(absent$global_boundary_classified)
  expect_false(absent$finite_jmle_existence_certified)
  expect_false(absent$standard_error_eligible)
  expect_false(absent$confidence_interval_eligible)
  expect_false(absent$external_comparison_eligible)
})

test_that("model estimator unit-slope and source boundaries fail closed", {
  source <- global_existence_source_fixture()

  pcm <- source
  pcm$config$model <- "PCM"
  pcm_result <- call_global_existence_audit(pcm)
  expect_identical(pcm_result$state, "not_applicable_model")
  expect_false(pcm_result$evaluated)
  expect_identical(pcm_result$estimator_identity, "not_applicable")

  mml <- source
  mml$config$method <- "MML"
  mml_result <- call_global_existence_audit(mml)
  expect_identical(mml_result$state, "not_applicable_estimator")
  expect_false(mml_result$evaluated)
  expect_identical(
    mml_result$objective_identity, "not_applicable_marginal_estimator"
  )

  unit <- source
  unit$config$gpcm_spec$n_params <- 0L
  unit_result <- call_global_existence_audit(unit)
  expect_identical(unit_result$state, "not_required_unit_slope")
  expect_false(unit_result$evaluated)
  expect_false(unit_result$finite_jmle_nonexistence_certified)

  legacy <- source
  legacy$sequence$contract_version <- "legacy"
  legacy_result <- call_global_existence_audit(legacy)
  expect_identical(legacy_result$state, "not_evaluated_source_contract")
  expect_false(legacy_result$evaluated)
  expect_false(legacy_result$source_parameter_sequence_flag_contract_matches)
  expect_false(legacy_result$finite_jmle_nonexistence_certified)
  expect_false(legacy_result$global_boundary_classified)
})
