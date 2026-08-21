gpcm_fixed_objective_classifier_config <- function(method = "JML",
                                                    model = "GPCM",
                                                    free = TRUE) {
  list(
    method = method,
    model = model,
    gpcm_spec = list(
      n_params = if (isTRUE(free)) 1L else 0L,
      identification = "sum_to_zero_log_slopes"
    )
  )
}

gpcm_fixed_objective_classifier_audit <- function(
    state = "none_certified",
    complete = TRUE,
    scope_complete = TRUE,
    retained_log_likelihood = -10,
    optimizer_log_likelihood = -10,
    certificates = data.frame()) {
  list(
    method = "JML",
    model = "GPCM",
    state = state,
    complete = complete,
    scope_complete = scope_complete,
    retained_log_likelihood = retained_log_likelihood,
    optimizer_log_likelihood = optimizer_log_likelihood,
    certificates = certificates
  )
}

test_that("fixed-objective classifier types a slope recession", {
  classified <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
    config = gpcm_fixed_objective_classifier_config(),
    slope_audit = gpcm_fixed_objective_classifier_audit(
      "certified_monotone_boundary_path"
    ),
    joint_audit = gpcm_fixed_objective_classifier_audit()
  )

  expect_identical(
    classified$contract_version,
    "mfrmr-jml-gpcm-fixed-objective-boundary-0.2.3-v2"
  )
  expect_identical(
    classified$state, "certified_slope_only_recession"
  )
  expect_identical(
    classified$estimator_identity,
    "unpenalized_fixed_effects_jml_no_finite_box"
  )
  expect_identical(
    classified$objective_identity,
    "identified_conditional_joint_log_likelihood"
  )
  expect_identical(classified$penalty, "none")
  expect_false(classified$finite_parameter_box)
  expect_true(classified$retained_point_finite)
  expect_identical(
    classified$retained_point_state,
    "finite_optimizer_trace_not_finite_maximum_claim"
  )
  expect_true(classified$slope_only_recession_certified)
  expect_false(classified$competitive_joint_boundary_certified)
  expect_true(classified$positive_boundary_certificate_present)
  expect_true(classified$audited_path_family_result_classified)
  expect_true(classified$fixed_objective_boundary_presence_classified)
  expect_false(classified$global_boundary_classified)
  expect_false(classified$global_finite_maximum_certified)
  expect_false(classified$global_boundary_absence_certified)
  expect_false(classified$standard_error_eligible)
  expect_false(classified$confidence_interval_eligible)
  expect_false(classified$external_comparison_eligible)
  expect_identical(classified$readiness_effect, "none_diagnostic_only")
})

test_that("fixed-objective classifier separates joint and combined positives", {
  config <- gpcm_fixed_objective_classifier_config()
  none <- gpcm_fixed_objective_classifier_audit()
  joint <- gpcm_fixed_objective_classifier_audit(
    "certified_competitive_joint_boundary_path"
  )
  slope <- gpcm_fixed_objective_classifier_audit(
    "certified_monotone_boundary_path"
  )

  joint_only <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
    config, none, joint
  )
  both <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
    config, slope, joint
  )

  expect_identical(
    joint_only$state, "certified_competitive_joint_boundary"
  )
  expect_false(joint_only$slope_only_recession_certified)
  expect_true(joint_only$competitive_joint_boundary_certified)
  expect_match(joint_only$reason_codes, "global_status_open", fixed = TRUE)
  expect_identical(both$state, "certified_slope_and_joint_boundary")
  expect_true(both$slope_only_recession_certified)
  expect_true(both$competitive_joint_boundary_certified)
  expect_false(both$global_boundary_classified)
})

test_that("two scoped negatives retain only a finite optimizer trace", {
  classified <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
    config = gpcm_fixed_objective_classifier_config(),
    slope_audit = gpcm_fixed_objective_classifier_audit(),
    joint_audit = gpcm_fixed_objective_classifier_audit()
  )

  expect_identical(
    classified$state, "finite_retained_point_no_path_certified"
  )
  expect_true(classified$retained_point_finite)
  expect_true(classified$path_family_classification_complete)
  expect_true(classified$audited_path_family_result_classified)
  expect_false(classified$fixed_objective_boundary_presence_classified)
  expect_true(classified$no_path_certified_in_audited_families)
  expect_false(classified$positive_boundary_certificate_present)
  expect_false(classified$global_boundary_classified)
  expect_false(classified$global_finite_maximum_certified)
  expect_false(classified$global_boundary_absence_certified)
  expect_match(
    classified$reason_codes,
    "finite_optimizer_trace_not_finite_maximum_claim",
    fixed = TRUE
  )
  expect_match(
    classified$limitations,
    "does not establish a finite global maximum",
    fixed = TRUE
  )
})

test_that("numerical indeterminacy and execution limits remain distinct", {
  config <- gpcm_fixed_objective_classifier_config()
  none <- gpcm_fixed_objective_classifier_audit()
  numerical <- gpcm_fixed_objective_classifier_audit(
    state = "not_evaluated_likelihood_mismatch",
    complete = FALSE,
    scope_complete = FALSE,
    retained_log_likelihood = NA_real_,
    optimizer_log_likelihood = NA_real_
  )
  limited <- gpcm_fixed_objective_classifier_audit(
    state = "not_evaluated_size_limit",
    complete = FALSE,
    scope_complete = FALSE,
    retained_log_likelihood = NA_real_,
    optimizer_log_likelihood = NA_real_
  )

  indeterminate <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
      config, none, numerical
    )
  not_evaluated <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
      config, none, limited
    )

  expect_identical(indeterminate$state, "indeterminate_numerical")
  expect_true(indeterminate$numerical_indeterminacy)
  expect_false(indeterminate$audited_path_family_result_classified)
  expect_false(indeterminate$fixed_objective_boundary_presence_classified)
  expect_false(indeterminate$no_path_certified_in_audited_families)
  expect_identical(not_evaluated$state, "not_evaluated")
  expect_false(not_evaluated$numerical_indeterminacy)
  expect_false(not_evaluated$path_family_classification_complete)
})

test_that("classifier fails closed across objective families and reductions", {
  none <- gpcm_fixed_objective_classifier_audit()
  mml <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
    gpcm_fixed_objective_classifier_config(method = "MML"), none, none
  )
  pcm <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
    gpcm_fixed_objective_classifier_config(model = "PCM"), none, none
  )
  unit <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
    gpcm_fixed_objective_classifier_config(free = FALSE), none, none
  )
  mismatched <- none
  mismatched$method <- "MML"
  mismatch <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
    gpcm_fixed_objective_classifier_config(), none, mismatched
  )
  incomplete_positive <- gpcm_fixed_objective_classifier_audit(
    state = "certified_monotone_boundary_path",
    complete = FALSE,
    scope_complete = FALSE
  )
  incomplete <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
    gpcm_fixed_objective_classifier_config(), incomplete_positive, none
  )

  expect_identical(mml$state, "not_applicable_estimator")
  expect_identical(mml$objective_identity, "not_applicable_marginal_estimator")
  expect_identical(pcm$state, "not_applicable_model")
  expect_identical(unit$state, "not_required_unit_slope")
  expect_true(unit$audited_path_family_result_classified)
  expect_false(unit$fixed_objective_boundary_presence_classified)
  expect_false(unit$global_finite_maximum_certified)
  expect_identical(mismatch$state, "not_evaluated_objective_identity")
  expect_false(mismatch$fixed_objective_boundary_presence_classified)
  expect_identical(incomplete$state, "not_evaluated")
  expect_false(incomplete$positive_boundary_certificate_present)
  expect_false(incomplete$audited_path_family_result_classified)
})
