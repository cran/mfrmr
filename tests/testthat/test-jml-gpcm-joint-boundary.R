make_gpcm_joint_boundary_fixture <- function(reverse = FALSE,
                                              nonseparated = FALSE) {
  if (isTRUE(nonseparated)) {
    data <- expand.grid(
      Person = c("P1", "P2"),
      Criterion = c("C1", "C2"),
      Score = 0:1,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    data <- data[rep(seq_len(nrow(data)), each = 2L), , drop = FALSE]
  } else {
    data <- expand.grid(
      Person = c("P1", "P2"),
      Criterion = c("C1", "C2"),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    data$Score <- c(0L, 1L, 1L, 0L)
    data <- data[rep(seq_len(nrow(data)), each = 3L), , drop = FALSE]
  }
  if (isTRUE(reverse)) data <- data[nrow(data):1L, , drop = FALSE]
  rownames(data) <- NULL
  data
}

fit_gpcm_joint_boundary_fixture <- function(reverse = FALSE,
                                             nonseparated = FALSE) {
  suppressWarnings(fit_mfrm(
    make_gpcm_joint_boundary_fixture(reverse, nonseparated),
    person = "Person",
    facets = "Criterion",
    score = "Score",
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    rating_min = 0,
    rating_max = 1,
    maxit = 300,
    reltol = 1e-12
  ))
}

gpcm_joint_boundary_problem <- function(fit) {
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

test_that("joint GPCM path catches the unanchored checkerboard", {
  skip_if_not_installed("lpSolve")
  fit <- fit_gpcm_joint_boundary_fixture()
  slope_only <- fit$config$boundary_audit$gpcm_slope_boundary
  audit <- fit$config$boundary_audit$gpcm_joint_boundary
  classified <-
    fit$config$boundary_audit$gpcm_fixed_objective_classification

  expect_identical(slope_only$state, "none_certified")
  expect_identical(
    audit$state, "certified_competitive_joint_boundary_path"
  )
  expect_true(isTRUE(audit$complete))
  expect_true(isTRUE(audit$scope_complete))
  expect_false(isTRUE(audit$structural_identification_complete))
  expect_lt(abs(audit$likelihood_difference), 1e-10)

  certified <- audit$certificates[audit$certificates$Certified, , drop = FALSE]
  expect_identical(nrow(certified), 2L)
  expect_true(all(certified$AnalyticTailCertified))
  expect_true(all(certified$CompetitiveBoundary))
  expect_true(all(certified$PositiveMinimumMargin > 0))
  expect_true(all(certified$NegativeLeadingCoefficient > 0))
  expect_true(all(certified$BoundaryImprovement > 0))
  expect_setequal(certified$PositiveLevel, c("C1", "C2"))
  expect_setequal(certified$NegativeLevel, c("C1", "C2"))
  expect_true(all(
    audit$target_status$CandidateStatus == "joint_boundary_path_both"
  ))

  for (pair_id in certified$PairId) {
    expanded <- audit$direction_loadings[
      audit$direction_loadings$PairId == pair_id &
        audit$direction_loadings$CoordinateType == "expanded_log_slope",
      , drop = FALSE
    ]
    expect_equal(sum(expanded$Loading), 0)
    expect_setequal(expanded$Loading, c(-1, 1))
    expect_true(any(
      audit$direction_loadings$PairId == pair_id &
        audit$direction_loadings$CoordinateType == "optimizer_additive"
    ))
  }

  # A competitive joint path is not promoted to a global unbounded estimate.
  expect_true(all(fit$slopes$ParameterStatus == "not_evaluated"))
  expect_true(all(is.na(fit$slopes$PrimaryLogEstimate)))
  expect_true(all(is.na(fit$slopes$PrimaryEstimate)))
  expect_true(all(!fit$slopes$SEEligible))
  expect_true(all(!fit$slopes$CIEligible))
  expect_true(all(fit$slopes$ComparisonEligibility == "ineligible"))
  expect_true(all(
    fit$slopes$ReasonCodes == "jml_gpcm_joint_boundary_candidate_both"
  ))
  expect_true(all(
    fit$slopes$PrimaryEstimateBasis ==
      "certified_joint_nonlinear_boundary_candidate_global_status_open"
  ))
  expect_match(
    fit$readiness$fit$ReasonCodes,
    "jml_gpcm_joint_boundary_candidate",
    fixed = TRUE
  )
  expect_identical(fit$readiness$fit$BoundaryState, "not_evaluated")
  expect_identical(fit$readiness$fit$FitReadiness, "review")
  expect_false(fit$readiness$fit$InferenceReady)
  expect_match(audit$limitations, "not a proof", fixed = TRUE)
  expect_identical(
    classified$state, "certified_competitive_joint_boundary"
  )
  expect_false(classified$slope_only_recession_certified)
  expect_true(classified$competitive_joint_boundary_certified)
  expect_false(classified$global_boundary_classified)
  expect_false(classified$global_finite_maximum_certified)
})

test_that("joint GPCM certificate agrees with a direct likelihood path", {
  skip_if_not_installed("lpSolve")
  fit <- fit_gpcm_joint_boundary_fixture()
  audit <- fit$config$boundary_audit$gpcm_joint_boundary
  certificate <- audit$certificates[
    which(audit$certificates$Certified)[1], , drop = FALSE
  ]
  loading <- audit$direction_loadings[
    audit$direction_loadings$PairId == certificate$PairId &
      audit$direction_loadings$CoordinateType %in%
        c("optimizer_additive", "optimizer_log_slope"),
    , drop = FALSE
  ]
  problem <- gpcm_joint_boundary_problem(fit)
  path_t <- c(0, 0.25, 0.5, 1, 2, 4, 8, 16)
  path_log_likelihood <- vapply(path_t, function(distance) {
    par <- fit$opt$par
    par[loading$OptimizerIndex] <-
      par[loading$OptimizerIndex] + distance * loading$Loading
    -mfrmr:::mfrm_loglik_jml(
      par, problem$idx, problem$config, problem$sizes
    )
  }, numeric(1))

  expect_gte(min(diff(path_log_likelihood)), -1e-10)
  expect_gt(path_log_likelihood[length(path_log_likelihood)],
            path_log_likelihood[1])
  expect_equal(
    path_log_likelihood[1], certificate$CurrentLogLikelihood,
    tolerance = 1e-10
  )
  expect_equal(
    path_log_likelihood[length(path_log_likelihood)],
    certificate$BoundaryLogLikelihood,
    tolerance = 1e-5
  )
})

test_that("joint GPCM result is invariant to retained row order", {
  skip_if_not_installed("lpSolve")
  forward <- fit_gpcm_joint_boundary_fixture()
  reversed <- fit_gpcm_joint_boundary_fixture(reverse = TRUE)
  a <- forward$config$boundary_audit$gpcm_joint_boundary
  b <- reversed$config$boundary_audit$gpcm_joint_boundary

  expect_identical(b$state, a$state)
  expect_identical(
    b$target_status$CandidateStatus,
    a$target_status$CandidateStatus
  )
  a_pairs <- a$certificates[
    order(a$certificates$PositiveLevel, a$certificates$NegativeLevel),
    , drop = FALSE
  ]
  b_pairs <- b$certificates[
    order(b$certificates$PositiveLevel, b$certificates$NegativeLevel),
    , drop = FALSE
  ]
  expect_identical(b_pairs$Certified, a_pairs$Certified)
  expect_equal(
    b_pairs$BoundaryImprovement, a_pairs$BoundaryImprovement,
    tolerance = 1e-10
  )
})

test_that("nonseparated data are negative only within the bounded path family", {
  skip_if_not_installed("lpSolve")
  fit <- fit_gpcm_joint_boundary_fixture(nonseparated = TRUE)
  audit <- fit$config$boundary_audit$gpcm_joint_boundary
  classified <-
    fit$config$boundary_audit$gpcm_fixed_objective_classification

  expect_identical(audit$state, "none_certified")
  expect_true(isTRUE(audit$complete))
  expect_true(isTRUE(audit$scope_complete))
  expect_false(isTRUE(audit$structural_identification_complete))
  expect_true(all(!audit$certificates$Certified))
  expect_true(all(
    audit$target_status$CandidateStatus ==
      "none_certified_in_joint_pair_paths"
  ))
  expect_true(all(fit$slopes$ParameterStatus == "not_evaluated"))
  expect_true(all(
    fit$slopes$ReasonCodes == "jml_gpcm_joint_boundary_none_certified"
  ))
  expect_match(audit$limitations, "negative result does not establish", fixed = TRUE)
  expect_identical(
    classified$state, "finite_retained_point_no_path_certified"
  )
  expect_true(classified$retained_point_finite)
  expect_true(classified$no_path_certified_in_audited_families)
  expect_false(classified$global_finite_maximum_certified)
  expect_false(classified$global_boundary_absence_certified)
})

test_that("joint GPCM audit fails closed on workload and estimator scope", {
  skip_if_not_installed("lpSolve")
  fit <- fit_gpcm_joint_boundary_fixture()
  problem <- gpcm_joint_boundary_problem(fit)

  limited <- mfrmr:::audit_mfrm_jml_gpcm_joint_boundary(
    fit$prep, problem$idx, problem$config, problem$sizes, fit$opt$par,
    max_direction_pairs = 0L
  )
  expect_identical(limited$state, "not_evaluated_size_limit")
  expect_false(isTRUE(limited$complete))
  expect_false(isTRUE(limited$scope_complete))
  expect_true(all(
    limited$target_status$CandidateStatus == "not_evaluated"
  ))

  sparse_limit <- mfrmr:::audit_mfrm_jml_gpcm_joint_boundary(
    fit$prep, problem$idx, problem$config, problem$sizes, fit$opt$par,
    max_lp_nonzeros = 0L
  )
  expect_identical(sparse_limit$state, "not_evaluated_pair_solver")
  expect_false(isTRUE(sparse_limit$complete))
  expect_true(any(
    sparse_limit$certificates$ReasonCodes == "joint_pair_lp_nonzero_limit"
  ))

  invalid <- mfrmr:::audit_mfrm_jml_gpcm_joint_boundary(
    fit$prep, problem$idx, problem$config, problem$sizes, fit$opt$par,
    certificate_tolerance = NA_real_
  )
  expect_identical(invalid$state, "not_evaluated_control")
  expect_false(isTRUE(invalid$complete))

  mml_config <- problem$config
  mml_config$method <- "MML"
  mml <- mfrmr:::audit_mfrm_jml_gpcm_joint_boundary(
    fit$prep, problem$idx, mml_config, problem$sizes, fit$opt$par
  )
  expect_identical(mml$state, "not_applicable_mml")
  expect_true(isTRUE(mml$complete))
  expect_false(isTRUE(mml$scope_complete))
  expect_false(isTRUE(mml$structural_identification_complete))
})
