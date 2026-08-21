make_jml_gpcm_terminal_gradient_fixture <- function(nonseparated = FALSE) {
  if (isTRUE(nonseparated)) {
    data <- expand.grid(
      Person = c("P1", "P2"),
      Criterion = c("C1", "C2"),
      Score = 0:1,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    return(data[rep(seq_len(nrow(data)), each = 2L), , drop = FALSE])
  }
  data <- expand.grid(
    Person = c("P1", "P2"),
    Criterion = c("C1", "C2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$Score <- c(0L, 1L, 1L, 0L)
  data[rep(seq_len(nrow(data)), each = 3L), , drop = FALSE]
}

fit_jml_gpcm_terminal_gradient_fixture <- function(
    boundary = c("joint", "slope", "negative")) {
  boundary <- match.arg(boundary)
  anchors <- if (identical(boundary, "slope")) {
    data.frame(
      Facet = "Person",
      Level = c("P1", "P2"),
      Anchor = c(-1, 1),
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
  suppressWarnings(fit_mfrm(
    make_jml_gpcm_terminal_gradient_fixture(
      nonseparated = identical(boundary, "negative")
    ),
    person = "Person",
    facets = "Criterion",
    score = "Score",
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    anchors = anchors,
    rating_min = 0,
    rating_max = 1,
    maxit = 300,
    reltol = 1e-12
  ))
}

jml_gpcm_terminal_gradient_problem <- function(fit) {
  config <- fit$config
  sizes <- mfrmr:::build_param_sizes(config)
  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = config$step_facet,
    slope_facet = config$slope_facet,
    interaction_specs = config$interaction_specs
  )
  list(config = config, sizes = sizes, idx = idx)
}

test_that("positive slope boundary takes precedence over a small gradient", {
  fit <- fit_jml_gpcm_terminal_gradient_fixture("slope")
  audit <- fit$config$boundary_audit$gpcm_terminal_gradient_stability

  expect_identical(
    audit$contract_version,
    "mfrmr-jml-gpcm-fixed-objective-terminal-gradient-0.2.3-v1"
  )
  expect_identical(audit$state, "boundary_path_gradient_nondecisive")
  expect_true(audit$evaluated)
  expect_identical(
    audit$estimator_identity,
    "unpenalized_fixed_effects_jml_no_finite_box"
  )
  expect_identical(
    audit$objective_identity,
    "identified_conditional_joint_log_likelihood"
  )
  expect_identical(
    audit$coordinate_system, "identified_optimizer_free_coordinates"
  )
  expect_true(audit$objective_reconstruction_agrees)
  expect_true(audit$stored_gradient_agrees)
  expect_true(audit$numeric_gradient_agrees)
  expect_true(audit$implementation_tolerance_met)
  expect_true(audit$boundary_precedence)
  expect_true(audit$optimizer_code_zero)
  expect_false(audit$retained_point_stationarity_supported)
  expect_false(audit$finite_interior_stationarity_certified)
  expect_false(audit$global_finite_maximum_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_false(audit$criterion_frozen_for_scientific_claims)
  expect_identical(audit$readiness_effect, "none_diagnostic_only")
  expect_identical(audit$polish_history_state, "consistent")
  expect_true(all(audit$numeric_probes$Evaluated))
  expect_true(all(
    audit$numeric_probes$ScaledDifference <= audit$numeric_probe_tolerance
  ))
  expect_true("log_slopes" %in% audit$block_summary$Block)
  expect_match(
    audit$limitations,
    "small gradient does not establish an interior finite maximum",
    fixed = TRUE
  )
})

test_that("competitive joint boundary also overrides a zero gradient", {
  fit <- fit_jml_gpcm_terminal_gradient_fixture("joint")
  boundary <-
    fit$config$boundary_audit$gpcm_fixed_objective_classification
  audit <- fit$config$boundary_audit$gpcm_terminal_gradient_stability

  expect_identical(
    boundary$state, "certified_competitive_joint_boundary"
  )
  expect_identical(audit$state, "boundary_path_gradient_nondecisive")
  expect_equal(audit$gradient_sup_norm, 0, tolerance = 1e-12)
  expect_true(audit$implementation_tolerance_met)
  expect_true(audit$boundary_precedence)
  expect_false(audit$retained_point_stationarity_supported)
  expect_false(audit$standard_error_eligible)
  expect_false(audit$confidence_interval_eligible)
  expect_false(audit$external_comparison_eligible)
})

test_that("scoped negatives allow retained-point first-order typing only", {
  fit <- fit_jml_gpcm_terminal_gradient_fixture("negative")
  boundary <-
    fit$config$boundary_audit$gpcm_fixed_objective_classification
  audit <- fit$config$boundary_audit$gpcm_terminal_gradient_stability

  expect_identical(
    boundary$state, "finite_retained_point_no_path_certified"
  )
  expect_identical(
    audit$state,
    "retained_point_gradient_within_implementation_tolerance"
  )
  expect_true(audit$evaluated)
  expect_true(audit$retained_point_stationarity_supported)
  expect_false(audit$boundary_precedence)
  expect_true(audit$implementation_tolerance_met)
  expect_true(all(audit$block_summary$WithinImplementationTolerance))
  expect_false(audit$finite_interior_stationarity_certified)
  expect_false(audit$global_finite_maximum_certified)
  expect_false(audit$global_boundary_absence_certified)
  expect_match(
    audit$reason_codes,
    "finite_global_maximum_not_certified",
    fixed = TRUE
  )
})

test_that("objective and stored-gradient drift fail closed", {
  fit <- fit_jml_gpcm_terminal_gradient_fixture("negative")
  problem <- jml_gpcm_terminal_gradient_problem(fit)
  boundary <-
    fit$config$boundary_audit$gpcm_fixed_objective_classification

  objective_drift <- fit$opt
  objective_drift$value <- objective_drift$value + 1
  objective_audit <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    par = fit$opt$par,
    opt = objective_drift,
    boundary_classification = boundary
  )
  expect_identical(objective_audit$state, "indeterminate_numerical")
  expect_false(objective_audit$evaluated)
  expect_false(objective_audit$objective_reconstruction_agrees)
  expect_false(objective_audit$retained_point_stationarity_supported)

  gradient_drift <- fit$opt
  gradient_drift$optimizer_diagnostics$TerminalGradientSupNorm <- 1
  gradient_audit <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    par = fit$opt$par,
    opt = gradient_drift,
    boundary_classification = boundary
  )
  expect_identical(gradient_audit$state, "indeterminate_numerical")
  expect_true(gradient_audit$objective_reconstruction_agrees)
  expect_false(gradient_audit$stored_gradient_agrees)
  expect_false(gradient_audit$retained_point_stationarity_supported)

  polish_drift <- fit$opt
  selected <- which(polish_drift$optimizer_polish$Stages$Selected)[1]
  polish_drift$optimizer_polish$Stages$Objective[selected] <-
    polish_drift$optimizer_polish$Stages$Objective[selected] + 1
  polish_audit <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    par = fit$opt$par,
    opt = polish_drift,
    boundary_classification = boundary
  )
  expect_identical(polish_audit$state, "indeterminate_numerical")
  expect_identical(polish_audit$polish_history_state, "inconsistent")
  expect_false(polish_audit$evaluated)
})

test_that("large gradients and optimizer warnings remain separate", {
  fit <- fit_jml_gpcm_terminal_gradient_fixture("negative")
  problem <- jml_gpcm_terminal_gradient_problem(fit)
  boundary <-
    fit$config$boundary_audit$gpcm_fixed_objective_classification

  moved_par <- fit$opt$par
  theta_index <- mfrmr:::build_param_slices(problem$sizes)$theta[1]
  moved_par[theta_index] <- moved_par[theta_index] + 0.5
  moved_gradient <- mfrmr:::mfrm_grad_jml(
    moved_par, problem$idx, problem$config, problem$sizes
  )
  moved_metrics <- mfrmr:::compute_gradient_metrics(moved_gradient)
  moved_opt <- fit$opt
  moved_opt$par <- moved_par
  moved_opt$value <- mfrmr:::mfrm_loglik_jml(
    moved_par, problem$idx, problem$config, problem$sizes
  )
  moved_opt$optimizer_diagnostics$TerminalGradientSupNorm <-
    moved_metrics$TerminalGradientSupNorm
  moved_opt$optimizer_diagnostics$TerminalGradientRMS <-
    moved_metrics$TerminalGradientRMS
  moved_opt$optimizer_diagnostics$GradientReviewTolerance <- 1e-8
  moved_opt$optimizer_polish <- NULL
  moved_audit <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    problem$idx, problem$config, problem$sizes, moved_par, moved_opt,
    boundary
  )
  expect_identical(
    moved_audit$state,
    "terminal_gradient_exceeds_implementation_tolerance"
  )
  expect_true(moved_audit$evaluated)
  expect_false(moved_audit$implementation_tolerance_met)
  expect_false(moved_audit$retained_point_stationarity_supported)
  expect_identical(moved_audit$reason_codes, "terminal_gradient_review")

  warning_opt <- fit$opt
  warning_opt$convergence <- 1L
  warning_opt$optimizer_diagnostics$ConvergenceCode <- 1L
  warning_audit <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    problem$idx, problem$config, problem$sizes, fit$opt$par, warning_opt,
    boundary
  )
  expect_identical(warning_audit$state, "small_gradient_optimizer_warning")
  expect_true(warning_audit$implementation_tolerance_met)
  expect_false(warning_audit$optimizer_code_zero)
  expect_false(warning_audit$retained_point_stationarity_supported)
})

test_that("open boundaries, limits, and estimator scope fail closed", {
  fit <- fit_jml_gpcm_terminal_gradient_fixture("negative")
  problem <- jml_gpcm_terminal_gradient_problem(fit)
  open_boundary <-
    fit$config$boundary_audit$gpcm_fixed_objective_classification
  open_boundary$no_path_certified_in_audited_families <- FALSE
  open_boundary$audited_path_family_result_classified <- FALSE
  open_audit <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    problem$idx, problem$config, problem$sizes, fit$opt$par, fit$opt,
    open_boundary
  )
  expect_identical(
    open_audit$state, "retained_point_gradient_small_boundary_open"
  )
  expect_true(open_audit$implementation_tolerance_met)
  expect_false(open_audit$retained_point_stationarity_supported)

  limited <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    problem$idx, problem$config, problem$sizes, fit$opt$par, fit$opt,
    open_boundary, max_free_dimension = 1L
  )
  expect_identical(limited$state, "not_evaluated_size_limit")
  expect_false(limited$evaluated)

  mml_config <- problem$config
  mml_config$method <- "MML"
  mml <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    problem$idx, mml_config, problem$sizes, fit$opt$par, fit$opt,
    open_boundary
  )
  expect_identical(mml$state, "not_applicable_estimator")
  expect_identical(mml$objective_identity, "not_applicable_marginal_estimator")

  pcm_config <- problem$config
  pcm_config$model <- "PCM"
  pcm <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    problem$idx, pcm_config, problem$sizes, fit$opt$par, fit$opt,
    open_boundary
  )
  expect_identical(pcm$state, "not_applicable_model")
  expect_identical(pcm$objective_identity, "not_applicable")

  unit_config <- problem$config
  unit_config$gpcm_spec$n_params <- 0L
  unit <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    problem$idx, unit_config, problem$sizes, fit$opt$par, fit$opt,
    open_boundary
  )
  expect_identical(unit$state, "not_required_unit_slope")
  expect_false(unit$finite_interior_stationarity_certified)

  malformed <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    problem$idx, problem$config, problem$sizes, fit$opt$par[-1], fit$opt,
    open_boundary
  )
  expect_identical(malformed$state, "not_evaluated_parameter_vector")
  expect_false(malformed$evaluated)

  invalid <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    problem$idx, problem$config, problem$sizes, fit$opt$par, fit$opt,
    open_boundary, derivative_tolerance = NA_real_
  )
  expect_identical(invalid$state, "not_evaluated_control")
  expect_false(invalid$evaluated)
})

test_that("over-strict numeric derivative agreement is indeterminate", {
  fit <- fit_jml_gpcm_terminal_gradient_fixture("slope")
  problem <- jml_gpcm_terminal_gradient_problem(fit)
  boundary <-
    fit$config$boundary_audit$gpcm_fixed_objective_classification
  audit <- mfrmr:::audit_mfrm_jml_gpcm_terminal_gradient(
    problem$idx, problem$config, problem$sizes, fit$opt$par, fit$opt,
    boundary, derivative_tolerance = 1e-15
  )

  expect_identical(audit$state, "boundary_path_gradient_indeterminate")
  expect_false(audit$numeric_gradient_agrees)
  expect_false(audit$evaluated)
  expect_true(audit$boundary_precedence)
  expect_false(audit$retained_point_stationarity_supported)
})
