fit_jml_gpcm_lexicographic_limit_fixture <- function() {
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


test_that("production fits advertise the oracle without inferring a path", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_lexicographic_limit_fixture()
  audit <- fit$config$boundary_audit$gpcm_lexicographic_limit

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-lexicographic-limit-0.2.3-v1"
  )
  expect_identical(
    audit$state,
    "declared_two_stage_limit_oracle_available_no_path_inferred"
  )
  expect_true(audit$evaluated)
  expect_true(audit$declared_two_stage_limit_oracle_available)
  expect_true(audit$exact_coefficient_comparisons_required)
  expect_identical(audit$maximum_slope_stages, 2L)
  expect_identical(audit$maximum_additive_stages, 2L)
  expect_true(audit$positive_power_scales_only)
  expect_true(audit$exponential_slope_vs_polynomial_additive_dominance_applied)
  expect_false(
    audit$declared_utility_path_parameter_reachability_checked
  )
  expect_false(audit$path_inferred_from_optimizer_trace)
  expect_false(audit$path_search_performed)
  expect_false(audit$production_path_limit_classified)
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
  expect_true(audit$source_rate_hierarchy_contract_matches)
  expect_match(audit$limitations, "does not derive a path", fixed = TRUE)
})

test_that("fit-level oracle scope does not borrow MML or FACETS authority", {
  fit <- fit_jml_gpcm_lexicographic_limit_fixture()
  config <- fit$config
  sizes <- mfrmr:::build_param_sizes(config)
  hierarchy <- config$boundary_audit$gpcm_rate_hierarchy

  mismatch <- hierarchy
  mismatch$contract_version <- "mismatched"
  structural <- mfrmr:::audit_mfrm_jml_gpcm_lexicographic_limit_scope(
    config, sizes, mismatch
  )
  expect_identical(
    structural$state,
    "declared_two_stage_limit_oracle_available_no_path_inferred"
  )
  expect_false(structural$source_rate_hierarchy_contract_matches)
  expect_false(structural$production_path_limit_classified)

  mml_config <- config
  mml_config$method <- "MML"
  mml <- mfrmr:::audit_mfrm_jml_gpcm_lexicographic_limit_scope(
    mml_config, sizes, hierarchy
  )
  expect_identical(mml$state, "not_applicable_estimator")
  expect_identical(mml$objective_identity, "not_applicable_marginal_estimator")
  expect_false(mml$external_comparison_eligible)

  pcm_config <- config
  pcm_config$model <- "PCM"
  pcm <- mfrmr:::audit_mfrm_jml_gpcm_lexicographic_limit_scope(
    pcm_config, sizes, hierarchy
  )
  expect_identical(pcm$state, "not_applicable_model")
  expect_identical(pcm$objective_identity, "not_applicable")

  unit_config <- config
  unit_config$gpcm_spec$levels <- "C1"
  unit_config$gpcm_spec$n_params <- 0L
  unit_sizes <- sizes
  unit_sizes$log_slopes <- 0L
  unit <- mfrmr:::audit_mfrm_jml_gpcm_lexicographic_limit_scope(
    unit_config, unit_sizes, hierarchy
  )
  expect_identical(unit$state, "not_required_unit_slope")
  expect_false(unit$declared_two_stage_limit_oracle_available)

  invalid_sizes <- sizes
  invalid_sizes$log_slopes <- 99L
  invalid <- mfrmr:::audit_mfrm_jml_gpcm_lexicographic_limit_scope(
    config, invalid_sizes, hierarchy
  )
  expect_identical(invalid$state, "not_evaluated_dimension")
  expect_false(invalid$evaluated)
})
