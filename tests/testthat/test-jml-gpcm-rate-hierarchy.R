fit_jml_gpcm_rate_hierarchy_fixture <- function() {
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


test_that("production scope records a finite depth but no likelihood claim", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_rate_hierarchy_fixture()
  audit <- fit$config$boundary_audit$gpcm_rate_hierarchy

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-rate-hierarchy-0.2.3-v1"
  )
  expect_identical(
    audit$state,
    "finite_coordinate_rate_hierarchy_available_likelihood_open"
  )
  expect_true(audit$evaluated)
  expect_identical(audit$slope_levels, 3L)
  expect_identical(audit$free_log_slope_coordinates, 2L)
  expect_true(audit$finite_coordinate_rate_hierarchy_theorem)
  expect_identical(audit$maximum_nonzero_slope_hierarchy_stages, 2L)
  expect_false(audit$later_active_restriction_sum_zero_required)
  expect_true(audit$global_coefficient_sum_zero_required)
  expect_true(audit$pathwise_hierarchy_extraction_available)
  expect_false(audit$arbitrary_fit_path_hierarchy_computed)
  expect_false(audit$arbitrary_sequence_likelihood_classified)
  expect_false(
    audit$common_accumulation_direction_likelihood_limit_certified
  )
  expect_false(audit$additive_hierarchy_classified)
  expect_true(audit$secondary_hierarchy_structurally_bounded)
  expect_false(audit$secondary_hierarchy_likelihood_classified)
  expect_false(audit$global_boundary_classified)
  expect_false(audit$global_finite_maximum_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
  expect_true(audit$source_compactification_contract_matches)
  expect_match(audit$limitations, "Additive directions", fixed = TRUE)
})

test_that("production hierarchy scope does not borrow MML or FACETS authority", {
  fit <- fit_jml_gpcm_rate_hierarchy_fixture()
  config <- fit$config
  sizes <- mfrmr:::build_param_sizes(config)
  compact <- config$boundary_audit$gpcm_boundary_compactification

  mismatch <- compact
  mismatch$contract_version <- "mismatched"
  structural <- mfrmr:::audit_mfrm_jml_gpcm_rate_hierarchy_scope(
    config, sizes, mismatch
  )
  expect_identical(
    structural$state,
    "finite_coordinate_rate_hierarchy_available_likelihood_open"
  )
  expect_false(structural$source_compactification_contract_matches)
  expect_false(structural$global_boundary_classified)

  mml_config <- config
  mml_config$method <- "MML"
  mml <- mfrmr:::audit_mfrm_jml_gpcm_rate_hierarchy_scope(
    mml_config, sizes, compact
  )
  expect_identical(mml$state, "not_applicable_estimator")
  expect_identical(mml$objective_identity, "not_applicable_marginal_estimator")
  expect_false(mml$external_comparison_eligible)

  pcm_config <- config
  pcm_config$model <- "PCM"
  pcm <- mfrmr:::audit_mfrm_jml_gpcm_rate_hierarchy_scope(
    pcm_config, sizes, compact
  )
  expect_identical(pcm$state, "not_applicable_model")
  expect_identical(pcm$objective_identity, "not_applicable")

  unit_config <- config
  unit_config$gpcm_spec$levels <- "C1"
  unit_config$gpcm_spec$n_params <- 0L
  unit_sizes <- sizes
  unit_sizes$log_slopes <- 0L
  unit <- mfrmr:::audit_mfrm_jml_gpcm_rate_hierarchy_scope(
    unit_config, unit_sizes, compact
  )
  expect_identical(unit$state, "not_required_unit_slope")
  expect_false(unit$finite_coordinate_rate_hierarchy_theorem)

  invalid_sizes <- sizes
  invalid_sizes$log_slopes <- 99L
  invalid <- mfrmr:::audit_mfrm_jml_gpcm_rate_hierarchy_scope(
    config, invalid_sizes, compact
  )
  expect_identical(invalid$state, "not_evaluated_dimension")
  expect_false(invalid$evaluated)
})
