make_structural_recession_rater_data <- function(separated = TRUE) {
  data <- expand.grid(
    Person = sprintf("P%02d", seq_len(12L)),
    Rater = c("R1", "R2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$Score <- if (isTRUE(separated)) {
    ifelse(data$Rater == "R1", 0L, 1L)
  } else {
    ifelse(data$Rater == "R1", 0L, rep(c(0L, 1L), 12L))
  }
  data
}

fit_structural_recession_rater <- function(data, ...) {
  suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = "Rater",
    score = "Score",
    rating_min = 0,
    rating_max = 1,
    method = "JML",
    model = "RSM",
    maxit = 80,
    ...
  ))
}

structural_recession_problem <- function(fit) {
  config <- fit$config
  sizes <- mfrmr:::build_param_sizes(config)
  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = config$step_facet,
    slope_facet = config$slope_facet,
    interaction_specs = config$interaction_specs
  )
  params <- mfrmr:::expand_params(fit$opt$par, sizes, config)
  adjacent <- mfrmr:::mfrmr_estimability_adjacent_design(
    prep = fit$prep,
    idx = idx,
    config = config,
    sizes = sizes,
    include_person = TRUE,
    include_population_beta = FALSE
  )
  structural_columns <- which(adjacent$map$Block != "Person")
  optimizer_index <- as.integer(
    adjacent$map$OptimizerIndex[structural_columns]
  )
  contrast <- mfrmr:::mfrmr_jml_observed_contrast_design(
    adjacent_design = adjacent$design[, structural_columns, drop = FALSE],
    score_k = idx$score_k,
    n_obs = adjacent$observation_rows,
    n_steps = adjacent$transitions
  )
  targets <- mfrmr:::mfrmr_jml_structural_target_system(
    config, sizes, params
  )
  list(
    config = config,
    sizes = sizes,
    idx = idx,
    params = params,
    contrast = contrast,
    target_metadata = targets$metadata,
    target_expansion = targets$expansion[, optimizer_index, drop = FALSE]
  )
}

finite_grid_recession_oracle <- function(contrast, target, tolerance = 1e-10) {
  contrast <- as.matrix(contrast)
  target <- as.numeric(target)
  n_parameters <- ncol(contrast)
  stopifnot(n_parameters >= 1L, n_parameters <= 8L)
  grid <- expand.grid(
    rep(list(c(-1, 0, 1)), n_parameters),
    KEEP.OUT.ATTRS = FALSE
  )
  direction <- as.matrix(grid)
  direction <- direction[rowSums(abs(direction)) > 0, , drop = FALSE]
  margins <- contrast %*% t(direction)
  feasible <- apply(margins >= -tolerance, 2L, all) &
    apply(margins > tolerance, 2L, any)
  target_change <- as.numeric(direction %*% target)
  c(
    positive = any(feasible & target_change > tolerance),
    negative = any(feasible & target_change < -tolerance)
  )
}

test_that("the additive recession fit policy is explicit and bounded", {
  policy <- mfrmr:::mfrmr_jml_recession_fit_policy()

  expect_identical(
    names(policy),
    c(
      "contract_version", "scope", "native_timeout_seconds",
      "attempts_per_stage", "retry_after_unaccepted",
      "maximum_native_seconds_per_positive_target"
    )
  )
  expect_identical(
    policy$contract_version,
    "mfrmr-jml-recession-fit-policy-v1"
  )
  expect_identical(
    policy$scope,
    "additive_structural_and_joint_recession"
  )
  expect_identical(policy$native_timeout_seconds, 10L)
  expect_identical(policy$attempts_per_stage, 1L)
  expect_false(policy$retry_after_unaccepted)
  expect_identical(
    policy$maximum_native_seconds_per_positive_target,
    20L
  )

  namespace <- asNamespace("mfrmr")
  expect_identical(
    eval(
      formals(mfrmr:::mfrmr_jml_recession_run_lp)$timeout,
      envir = namespace
    ),
    10L
  )
  expect_identical(
    eval(
      formals(mfrmr:::mfrmr_jml_recession_target_lp)$timeout,
      envir = namespace
    ),
    10L
  )
  expect_identical(
    eval(
      formals(mfrmr:::audit_mfrm_jml_additive_recession)$lp_timeout,
      envir = namespace
    ),
    10L
  )
  expect_identical(
    eval(
      formals(mfrmr:::audit_mfrm_jml_gpcm_joint_boundary)$lp_timeout,
      envir = namespace
    ),
    2L
  )
})

test_that("the target LP uses one policy-bound attempt per stage", {
  skip_if_not_installed("lpSolve")
  observed_timeout <- integer(0)
  testthat::local_mocked_bindings(
    mfrmr_jml_recession_run_lp = function(
        lp_base,
        objective,
        extra_constraint = NULL,
        extra_direction = NULL,
        extra_rhs = NULL,
        timeout) {
      observed_timeout <<- c(observed_timeout, as.integer(timeout))
      list(status = 0L, objval = 1, solution = c(1, 0))
    },
    .package = "mfrmr"
  )
  lp_base <- mfrmr:::mfrmr_jml_recession_lp_base(
    Matrix::sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(1L, 1L))
  )

  default <- mfrmr:::mfrmr_jml_recession_target_lp(
    lp_base = lp_base,
    target = 1
  )
  explicit <- mfrmr:::mfrmr_jml_recession_target_lp(
    lp_base = lp_base,
    target = 1,
    timeout = 3L
  )

  expect_true(default$certified)
  expect_true(explicit$certified)
  expect_identical(default$lp_calls, 2L)
  expect_identical(explicit$lp_calls, 2L)
  expect_identical(observed_timeout, c(10L, 10L, 3L, 3L))
})

test_that("the cone exclusion tolerance cannot hide a near-boundary target", {
  skip_if_not_installed("lpSolve")
  contrast <- Matrix::sparseMatrix(
    i = 1L, j = 1L, x = 5e-7, dims = c(1L, 1L)
  )
  lp_base <- mfrmr:::mfrmr_jml_recession_lp_base(contrast)

  target <- mfrmr:::mfrmr_jml_recession_target_lp(
    lp_base = lp_base, target = 1
  )
  unsafe_default_cone <- mfrmr:::mfrmr_jml_recession_target_lp(
    lp_base = lp_base, target = as.numeric(Matrix::colSums(contrast))
  )
  guarded_cone <- mfrmr:::mfrmr_jml_recession_target_lp(
    lp_base = lp_base,
    target = as.numeric(Matrix::colSums(contrast)),
    objective_tolerance = 1e-10,
    certificate_tolerance = 1e-7
  )

  expect_true(isTRUE(target$certified))
  expect_false(isTRUE(unsafe_default_cone$certified))
  expect_true(isTRUE(guarded_cone$certified))
  expect_gt(guarded_cone$positive_margin, 1e-7)
  expect_gt(guarded_cone$strict_rows, 0L)
})

test_that("sum-zero Rater separation receives an exact structural certificate", {
  skip_if_not_installed("lpSolve")
  fit <- fit_structural_recession_rater(
    make_structural_recession_rater_data(TRUE)
  )
  audit <- fit$config$boundary_audit$structural_additive
  rater <- audit$target_status[
    audit$target_status$ParameterClass == "facet", , drop = FALSE
  ]
  rater <- rater[match(c("R1", "R2"), rater$Level), , drop = FALSE]

  expect_identical(audit$state, "certified_recession")
  expect_true(isTRUE(audit$complete))
  expect_identical(
    rater$CandidateStatus,
    c("unbounded_high", "unbounded_low")
  )
  expect_identical(rater$PositiveRecession, c(TRUE, FALSE))
  expect_identical(rater$NegativeRecession, c(FALSE, TRUE))
  expect_true(all(is.finite(rater$OptimizerEstimate)))
  expect_true(all(audit$certificates$MinimumContrastMargin >= -1e-7))
  certified <- audit$certificates[audit$certificates$Certified, , drop = FALSE]
  expect_true(all(certified$StrictContrastRows > 0L))
  expect_true(all(certified$PositiveContrastMargin > 0))
  expect_true(nrow(audit$direction_loadings) > 0L)
  expect_identical(audit$dimensions$LPRepresentation, "sparse_triplet")
  expect_lte(
    audit$dimensions$SparseLPNonzeros,
    audit$dimensions$DenseReferenceElements
  )
  expect_identical(
    audit$prescreen$contract_version,
    "mfrmr-jml-global-cone-prescreen-v1"
  )
  expect_true(isTRUE(audit$prescreen$requested))
  expect_identical(
    audit$prescreen$state,
    "positive_cone_target_enumeration"
  )
  expect_true(isTRUE(audit$prescreen$evaluated))
  expect_true(isTRUE(audit$prescreen$cone_certified))
  expect_false(isTRUE(audit$prescreen$target_enumeration_skipped))
  expect_identical(
    audit$prescreen$target_directions_evaluated,
    audit$dimensions$TargetDirections
  )
  expect_gt(audit$prescreen$cone_lp_calls, 0L)
  expect_gte(
    audit$prescreen$target_lp_calls,
    audit$prescreen$target_directions_evaluated
  )
  expect_identical(
    audit$prescreen$total_lp_calls,
    audit$prescreen$cone_lp_calls + audit$prescreen$target_lp_calls
  )

  problem <- structural_recession_problem(fit)
  unscreened <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = fit$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params,
    screen_global_cone = FALSE
  )
  expect_identical(unscreened$state, audit$state)
  expect_identical(unscreened$complete, audit$complete)
  expect_identical(unscreened$target_status, audit$target_status)
  expect_identical(
    unscreened$certificates$Certified,
    audit$certificates$Certified
  )
  expect_identical(unscreened$prescreen$state, "not_requested")
  expect_identical(unscreened$prescreen$cone_lp_calls, 0L)
  expect_identical(
    unscreened$prescreen$target_directions_evaluated,
    unscreened$dimensions$TargetDirections
  )

  dense_reference <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = fit$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params,
    lp_representation = "dense_reference"
  )
  expect_identical(
    dense_reference$target_status$CandidateStatus,
    audit$target_status$CandidateStatus
  )
  expect_identical(
    dense_reference$certificates$Certified,
    audit$certificates$Certified
  )
  expect_equal(
    dense_reference$certificates$TargetCapacity,
    audit$certificates$TargetCapacity,
    tolerance = 1e-8
  )

  reordered_fit <- fit_structural_recession_rater(
    make_structural_recession_rater_data(TRUE)[24:1, , drop = FALSE]
  )
  reordered <- reordered_fit$config$boundary_audit$structural_additive
  target_order <- order(audit$target_status$ParameterId)
  reordered_target_order <- order(reordered$target_status$ParameterId)
  expect_identical(
    audit$target_status$CandidateStatus[target_order],
    reordered$target_status$CandidateStatus[reordered_target_order]
  )
  certificate_order <- order(
    audit$certificates$ParameterId,
    audit$certificates$RequestedDirection
  )
  reordered_certificate_order <- order(
    reordered$certificates$ParameterId,
    reordered$certificates$RequestedDirection
  )
  expect_equal(
    audit$certificates$TargetCapacity[certificate_order],
    reordered$certificates$TargetCapacity[reordered_certificate_order],
    tolerance = 1e-8
  )

  # This slice stores the certificate but deliberately does not yet mutate
  # public result tables. Cross-surface primary-value propagation is WP4.
  expect_true(all(is.finite(fit$facets$others$Estimate)))
})

test_that("a response-constant level alone does not prove structural recession", {
  skip_if_not_installed("lpSolve")
  fit <- fit_structural_recession_rater(
    make_structural_recession_rater_data(FALSE)
  )
  audit <- fit$config$boundary_audit$structural_additive
  rater <- audit$target_status[
    audit$target_status$ParameterClass == "facet", , drop = FALSE
  ]

  expect_identical(audit$state, "none_certified")
  expect_true(isTRUE(audit$complete))
  expect_true(all(rater$CandidateStatus == "finite_in_audited_subspace"))
  expect_false(any(rater$PositiveRecession))
  expect_false(any(rater$NegativeRecession))
  expect_true(any(
    as.character(fit$prep$data$Rater) == "R1" & fit$prep$data$Score == 0L
  ))
  expect_true(all(
    fit$prep$data$Score[as.character(fit$prep$data$Rater) == "R1"] == 0L
  ))

  problem <- structural_recession_problem(fit)
  unscreened <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = fit$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params,
    screen_global_cone = FALSE
  )
  expect_identical(unscreened$state, audit$state)
  expect_identical(unscreened$complete, audit$complete)
  expect_identical(unscreened$target_status, audit$target_status)
  expect_true(isTRUE(audit$prescreen$requested))
  expect_identical(
    audit$prescreen$state,
    "negative_cone_enumeration_skipped"
  )
  expect_true(isTRUE(audit$prescreen$evaluated))
  expect_false(isTRUE(audit$prescreen$cone_certified))
  expect_true(isTRUE(audit$prescreen$target_enumeration_skipped))
  expect_gt(nrow(unscreened$certificates), 0L)
  expect_identical(nrow(audit$certificates), 0L)
  expect_identical(audit$prescreen$target_directions_evaluated, 0L)
  expect_identical(audit$prescreen$target_lp_calls, 0L)
  expect_gt(audit$prescreen$cone_lp_calls, 0L)
  expect_lt(
    audit$prescreen$total_lp_calls,
    unscreened$prescreen$total_lp_calls
  )
})

test_that("structural contrasts use retained positive-weight observations", {
  skip_if_not_installed("lpSolve")
  data <- make_structural_recession_rater_data(TRUE)
  data$Weight <- 1
  zero_row <- which(data$Person == "P01" & data$Rater == "R1")
  missing_row <- which(data$Person == "P02" & data$Rater == "R2")
  data$Score[zero_row] <- 1L
  data$Weight[zero_row] <- 0
  data$Score[missing_row] <- NA_integer_

  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = "Rater",
    score = "Score",
    weight = "Weight",
    rating_min = 0,
    rating_max = 1,
    method = "JML",
    model = "RSM",
    maxit = 80
  ))
  audit <- fit$config$boundary_audit$structural_additive
  rater <- audit$target_status[
    audit$target_status$ParameterClass == "facet", , drop = FALSE
  ]
  rater <- rater[match(c("R1", "R2"), rater$Level), , drop = FALSE]

  expect_identical(nrow(fit$prep$data), 22L)
  expect_identical(audit$dimensions$Observations, 22L)
  expect_identical(
    rater$CandidateStatus,
    c("unbounded_high", "unbounded_low")
  )
})

test_that("facet signs and direct anchors enter the structural certificate", {
  skip_if_not_installed("lpSolve")
  data <- make_structural_recession_rater_data(TRUE)
  positive_fit <- fit_structural_recession_rater(
    data,
    positive_facets = "Rater"
  )
  positive <- positive_fit$config$boundary_audit$structural_additive$target_status
  positive <- positive[positive$ParameterClass == "facet", , drop = FALSE]
  positive <- positive[match(c("R1", "R2"), positive$Level), , drop = FALSE]
  expect_identical(
    positive$CandidateStatus,
    c("unbounded_low", "unbounded_high")
  )

  anchored_fit <- fit_structural_recession_rater(
    data,
    anchors = data.frame(
      Facet = "Rater",
      Level = "R1",
      Anchor = 0,
      stringsAsFactors = FALSE
    )
  )
  anchored <- anchored_fit$config$boundary_audit$structural_additive$target_status
  anchored <- anchored[anchored$ParameterClass == "facet", , drop = FALSE]
  anchored <- anchored[match(c("R1", "R2"), anchored$Level), , drop = FALSE]
  expect_identical(anchored$CandidateStatus, c("fixed", "fixed"))
  expect_false(any(anchored$PositiveRecession))
  expect_false(any(anchored$NegativeRecession))
  expect_identical(
    anchored_fit$config$boundary_audit$structural_additive$state,
    "no_free_structural_coordinates"
  )
})

test_that("checkerboard interaction cells receive sign-specific certificates", {
  skip_if_not_installed("lpSolve")
  data <- expand.grid(
    Person = sprintf("P%02d", seq_len(12L)),
    Rater = c("R1", "R2"),
    Criterion = c("C1", "C2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$Score <- as.integer(
    (data$Rater == "R1") == (data$Criterion == "C1")
  )
  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 0,
    rating_max = 1,
    method = "JML",
    model = "RSM",
    facet_interactions = "Rater:Criterion",
    min_obs_per_interaction = 1,
    maxit = 80
  ))
  audit <- fit$config$boundary_audit$structural_additive
  interaction <- audit$target_status[
    audit$target_status$ParameterClass == "interaction", , drop = FALSE
  ]
  interaction <- interaction[
    match(c("R1*C1", "R2*C1", "R1*C2", "R2*C2"), interaction$Level),
    , drop = FALSE
  ]

  expect_identical(audit$state, "certified_recession")
  expect_identical(
    interaction$CandidateStatus,
    c("unbounded_high", "unbounded_low", "unbounded_low", "unbounded_high")
  )
  expect_true(all(interaction$EvaluationState == "evaluated"))

  problem <- structural_recession_problem(fit)
  dense_reference <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = fit$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params,
    lp_representation = "dense_reference"
  )
  dense_interaction <- dense_reference$target_status[
    dense_reference$target_status$ParameterClass == "interaction",
    , drop = FALSE
  ]
  dense_interaction <- dense_interaction[
    match(interaction$ParameterId, dense_interaction$ParameterId),
    , drop = FALSE
  ]
  expect_identical(
    dense_interaction$CandidateStatus,
    interaction$CandidateStatus
  )

  interaction_index <- which(
    problem$target_metadata$ParameterClass == "interaction"
  )
  oracle <- t(vapply(interaction_index, function(target_index) {
    finite_grid_recession_oracle(
      problem$contrast,
      problem$target_expansion[target_index, ]
    )
  }, numeric(2)))
  expect_identical(
    as.logical(oracle[, "positive"]),
    interaction$PositiveRecession
  )
  expect_identical(
    as.logical(oracle[, "negative"]),
    interaction$NegativeRecession
  )
})

test_that("sparse triplet constraints avoid the dense reference allocation", {
  skip_if_not_installed("lpSolve")
  n_rows <- 20000L
  n_parameters <- 100L
  contrast <- Matrix::sparseMatrix(
    i = seq_len(n_rows),
    j = ((seq_len(n_rows) - 1L) %% n_parameters) + 1L,
    x = rep(c(-1, 1), length.out = n_rows),
    dims = c(n_rows, n_parameters)
  )
  base <- mfrmr:::mfrmr_jml_recession_lp_base(
    contrast,
    representation = "sparse_triplet"
  )

  expect_null(base$constraint_matrix)
  expect_identical(base$representation, "sparse_triplet")
  expect_identical(
    nrow(base$constraint_triplet),
    2L * length(contrast@x) + 2L * n_parameters
  )
  expect_lt(
    nrow(base$constraint_triplet),
    (n_rows + n_parameters) * (2L * n_parameters) / 50
  )
})

test_that("MML is not reduced to the JML structural recession contract", {
  config <- list(method = "MML", model = "RSM")
  audit <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = list(),
    idx = list(),
    config = config,
    sizes = list(),
    params = list()
  )

  expect_identical(audit$state, "not_applicable_mml")
  expect_true(isTRUE(audit$complete))
  expect_identical(nrow(audit$target_status), 0L)
})

test_that("bounded execution limits fail closed without inventing a certificate", {
  skip_if_not_installed("lpSolve")
  fit <- fit_structural_recession_rater(
    make_structural_recession_rater_data(TRUE)
  )
  config <- fit$config
  sizes <- mfrmr:::build_param_sizes(config)
  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = config$step_facet,
    slope_facet = config$slope_facet,
    interaction_specs = config$interaction_specs
  )
  params <- mfrmr:::expand_params(fit$opt$par, sizes, config)
  audit <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = fit$prep,
    idx = idx,
    config = config,
    sizes = sizes,
    params = params,
    max_lp_nonzeros = 1
  )

  expect_identical(audit$state, "not_evaluated_size_limit")
  expect_false(isTRUE(audit$complete))
  expect_false(any(audit$target_status$PositiveRecession, na.rm = TRUE))
  expect_false(any(audit$target_status$NegativeRecession, na.rm = TRUE))
  expect_true(any(
    audit$target_status$CandidateStatus == "not_evaluated_size_limit"
  ))
  expect_identical(nrow(audit$certificates), 0L)
  expect_identical(audit$prescreen$state, "not_evaluated_size_limit")
  expect_identical(audit$prescreen$total_lp_calls, 0L)
})

test_that("a failed structural cone solve cannot become a finite result", {
  skip_if_not_installed("lpSolve")
  fit <- fit_structural_recession_rater(
    make_structural_recession_rater_data(FALSE)
  )
  problem <- structural_recession_problem(fit)

  testthat::local_mocked_bindings(
    mfrmr_jml_recession_target_lp = function(lp_base, target, ...) {
      list(
        evaluated = FALSE,
        certified = FALSE,
        solver_status = 5L,
        target_capacity = NA_real_,
        target_change = NA_real_,
        minimum_margin = NA_real_,
        positive_margin = NA_real_,
        strict_rows = NA_integer_,
        direction = rep(NA_real_, length(target)),
        lp_calls = 1L,
        reason = "forced_solver_failure"
      )
    },
    .package = "mfrmr"
  )

  audit <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = fit$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params
  )

  expect_identical(audit$state, "not_evaluated_solver")
  expect_false(isTRUE(audit$complete))
  expect_identical(audit$prescreen$state, "not_evaluated_solver")
  expect_false(isTRUE(audit$prescreen$evaluated))
  expect_identical(audit$prescreen$cone_lp_calls, 1L)
  expect_identical(audit$prescreen$target_lp_calls, 0L)
  expect_true(any(is.na(audit$target_status$PositiveRecession)))
  expect_true(any(
    audit$target_status$CandidateStatus == "not_evaluated_solver"
  ))
  expect_false(any(
    audit$target_status$CandidateStatus == "finite_in_audited_subspace"
  ))
})
