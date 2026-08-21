make_mml_gpcm_slope_boundary_fixture <- function(reverse = FALSE) {
  data <- expand.grid(
    Person = sprintf("P%02d", seq_len(12L)),
    Rater = c("R1", "R2"),
    Criterion = c("C1", "C2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$Score <- ifelse(data$Criterion == "C1", 1L, 0L)
  if (isTRUE(reverse)) data <- data[nrow(data):1L, , drop = FALSE]
  rownames(data) <- NULL
  data
}

fit_mml_gpcm_slope_boundary_fixture <- function(reverse = FALSE) {
  anchors <- data.frame(
    Facet = "Rater",
    Level = c("R1", "R2"),
    Anchor = c(-10, -10),
    stringsAsFactors = FALSE
  )
  suppressWarnings(fit_mfrm(
    make_mml_gpcm_slope_boundary_fixture(reverse),
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    anchors = anchors,
    rating_min = 0,
    rating_max = 1,
    quad_points = 5L,
    maxit = 10L
  ))
}

mml_gpcm_slope_boundary_problem <- function(fit) {
  config <- fit$config
  sizes <- mfrmr:::build_param_sizes(config)
  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = config$step_facet,
    slope_facet = config$slope_facet,
    interaction_specs = config$interaction_specs
  )
  list(
    config = config,
    sizes = sizes,
    idx = idx,
    slices = mfrmr:::build_param_slices(sizes),
    quad = mfrmr:::gauss_hermite_normal(5L)
  )
}

test_that("MML fixed-quadrature marginal slope path is certified", {
  fit <- fit_mml_gpcm_slope_boundary_fixture()
  audit <- fit$config$boundary_audit$gpcm_slope_boundary

  expect_identical(
    audit$state,
    "certified_fixed_quadrature_marginal_boundary_path"
  )
  expect_identical(audit$method, "MML")
  expect_true(audit$complete)
  expect_true(audit$scope_complete)
  expect_true(audit$fixed_quadrature_certificate)
  expect_false(audit$continuous_integral_certificate)
  expect_false(audit$structural_identification_complete)
  expect_identical(audit$readiness_effect, "none_instrumentation_only")
  expect_lt(abs(audit$likelihood_difference), 1e-10)
  expect_identical(audit$dimensions$QuadratureNodes, 5L)

  certified <- audit$certificates[
    audit$certificates$Certified, , drop = FALSE
  ]
  expect_identical(nrow(certified), 1L)
  expect_identical(certified$PositiveLevel, "C1")
  expect_identical(certified$NegativeLevel, "C2")
  expect_equal(certified$DirectionSum, 0)
  expect_gt(certified$BoundaryImprovement, 0)
  expect_gt(certified$StrictRows, 0L)

  # Instrumentation cannot turn the still-unfinished global marginal geometry
  # into finite primary values or inference-ready output.
  expect_true(all(fit$slopes$ParameterStatus == "not_evaluated"))
  expect_true(all(is.na(fit$slopes$PrimaryEstimate)))
  expect_true(all(is.finite(fit$slopes$OptimizerEstimate)))
  expect_true(all(
    fit$slopes$ReasonCodes == "mml_gpcm_slope_boundary_not_evaluated"
  ))
  expect_identical(fit$readiness$fit$BoundaryState, "not_evaluated")
  expect_false(fit$readiness$fit$InferenceReady)
})

test_that("MML marginal certificate agrees with a direct path oracle", {
  fit <- fit_mml_gpcm_slope_boundary_fixture()
  audit <- fit$config$boundary_audit$gpcm_slope_boundary
  certificate <- audit$certificates[
    audit$certificates$Certified, , drop = FALSE
  ]
  problem <- mml_gpcm_slope_boundary_problem(fit)

  expanded_direction <- c(C1 = 1, C2 = -1)
  free_direction <- expanded_direction[
    seq_len(length(expanded_direction) - 1L)
  ]
  distance <- c(0, 0.25, 0.5, 1, 2, 4, 8)
  path_log_likelihood <- vapply(distance, function(value) {
    par <- fit$opt$par
    par[problem$slices$log_slopes] <-
      par[problem$slices$log_slopes] + value * free_direction
    -mfrmr:::mfrm_loglik_mml(
      par, problem$idx, problem$config, problem$sizes, problem$quad
    )
  }, numeric(1))

  expect_equal(sum(expanded_direction), 0)
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
    tolerance = 5e-5
  )
})

test_that("MML marginal slope audit is row-order invariant and fail-closed", {
  forward <- fit_mml_gpcm_slope_boundary_fixture()
  reversed <- fit_mml_gpcm_slope_boundary_fixture(reverse = TRUE)
  a <- forward$config$boundary_audit$gpcm_slope_boundary
  b <- reversed$config$boundary_audit$gpcm_slope_boundary
  expect_identical(b$state, a$state)
  expect_identical(
    b$certificates$Certified,
    a$certificates$Certified
  )
  expect_equal(
    b$certificates$BoundaryImprovement,
    a$certificates$BoundaryImprovement,
    tolerance = 1e-10
  )

  problem <- mml_gpcm_slope_boundary_problem(forward)
  limited <- mfrmr:::audit_mfrm_mml_gpcm_slope_boundary(
    prep = forward$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    par = forward$opt$par,
    quad_points = 5L,
    max_pair_node_elements = 0
  )
  expect_identical(limited$state, "not_evaluated_size_limit")
  expect_false(limited$complete)
  expect_false(limited$scope_complete)
  expect_true(all(
    limited$target_status$CandidateStatus == "not_evaluated_size_limit"
  ))

  invalid <- mfrmr:::audit_mfrm_mml_gpcm_slope_boundary(
    prep = forward$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    par = forward$opt$par,
    quad_points = NA_integer_
  )
  expect_identical(invalid$state, "not_evaluated_control")
  expect_false(invalid$complete)
})

test_that("ordinary MML GPCM negative result is not a finite-MLE claim", {
  data <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 5L,
    maxit = 25L
  )))
  audit <- fit$config$boundary_audit$gpcm_slope_boundary

  expect_identical(
    audit$state,
    "none_certified_fixed_quadrature_marginal"
  )
  expect_true(audit$scope_complete)
  expect_false(audit$structural_identification_complete)
  expect_match(audit$limitations, "does not establish finite GPCM slopes",
               fixed = TRUE)
  expect_true(all(fit$slopes$ParameterStatus == "not_evaluated"))
  expect_false(fit$readiness$fit$InferenceReady)
})
