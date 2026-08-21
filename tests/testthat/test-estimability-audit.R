.estimability_balanced_data <- function() {
  data <- expand.grid(
    Person = paste0("P", sprintf("%02d", 1:8)),
    Rater = c("R1", "R2"),
    Criterion = c("C1", "C2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$Score <- 1L + (
    as.integer(factor(data$Person)) +
      as.integer(factor(data$Rater)) +
      as.integer(factor(data$Criterion))
  ) %% 4L
  data
}

.estimability_zero_common_data <- function() {
  data <- expand.grid(
    Person = paste0("P", sprintf("%02d", 1:8)),
    Criterion = c("C1", "C2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  person_index <- as.integer(sub("P", "", data$Person))
  data$Rater <- ifelse(person_index <= 4L, "R1", "R2")
  data$Score <- 1L + (
    as.integer(factor(data$Person)) +
      as.integer(factor(data$Criterion))
  ) %% 4L
  data
}

.estimability_fit <- function(data, method = "JML", ...) {
  fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    model = "RSM",
    method = method,
    quad_points = 5,
    maxit = 5,
    ...
  )
}

test_that("nonlinear estimability block selection never recycles", {
  sizes <- c(
    theta = 8L, facets = 3L, interactions = 0L, steps = 6L,
    log_slopes = 2L, beta = 1L, log_sigma2 = 1L
  )
  expect_warning(
    blocks <- mfrmr:::mfrmr_estimability_nonlinear_blocks(sizes),
    NA
  )
  expect_identical(blocks, c("log_slopes", "log_sigma2"))

  sizes[c("log_slopes", "log_sigma2")] <- 0L
  expect_identical(
    mfrmr:::mfrmr_estimability_nonlinear_blocks(sizes),
    character(0)
  )
})

test_that("bounded response-pattern enumeration is complete and unique", {
  patterns <- mfrmr:::mfrmr_enumerate_response_patterns(
    n_observations = 2L,
    n_categories = 3L
  )
  expect_identical(dim(patterns), c(9L, 2L))
  expect_identical(nrow(unique(as.data.frame(patterns))), 9L)
  expect_identical(patterns[1, ], c(0L, 0L))
  expect_identical(patterns[nrow(patterns), ], c(2L, 2L))
  expect_identical(sort(unique(as.integer(patterns))), 0:2)
})

test_that("balanced JML constrained design is full rank and invariant", {
  data <- .estimability_balanced_data()
  fit <- suppressWarnings(.estimability_fit(data, "JML"))
  audit <- fit$data_review$estimability

  expect_identical(
    audit$contract_version,
    "mfrmr-readiness-0.2.3-v3"
  )
  expect_identical(audit$readiness$EstimabilityState, "identified")
  expect_true(audit$readiness$Complete)
  expect_identical(audit$design$rank, audit$design$free_dimension)
  expect_identical(audit$design$nullity, 0L)
  expect_identical(
    audit$readiness$AuditedFreeDimension,
    audit$readiness$OptimizerFreeDimension
  )
  expect_identical(
    audit$fitted_information$status,
    "pending_weak_information_calibration"
  )
  expect_identical(audit$nonlinear_transformation$status, "not_required")
  expect_identical(audit$gpcm_response_kernel$status,
                   "not_applicable_model")
  expect_identical(audit$mml_observed_pattern_score$status,
                   "not_applicable_estimator")
  expect_identical(audit$mml_all_pattern_information$status,
                   "not_applicable_estimator")
  expect_false(audit$fitted_information$weak_information_classified)
  expect_true(all(c("Person", "Rater", "Criterion", "steps") %in%
                    audit$parameter_blocks$Block))

  permuted <- data[c(seq(2L, nrow(data), by = 2L),
                     seq(1L, nrow(data), by = 2L)), , drop = FALSE]
  relabelled <- data
  relabelled$Person <- paste0("Candidate_", match(
    relabelled$Person, rev(unique(relabelled$Person))
  ))
  relabelled$Rater <- ifelse(relabelled$Rater == "R1", "Judge_B", "Judge_A")
  relabelled$Criterion <- ifelse(
    relabelled$Criterion == "C1", "Scale_Y", "Scale_X"
  )

  permuted_audit <- suppressWarnings(
    .estimability_fit(permuted, "JML")
  )$data_review$estimability
  relabelled_audit <- suppressWarnings(
    .estimability_fit(relabelled, "JML")
  )$data_review$estimability

  expect_identical(permuted_audit$readiness$EstimabilityState, "identified")
  expect_identical(relabelled_audit$readiness$EstimabilityState, "identified")
  expect_identical(permuted_audit$design$rank, audit$design$rank)
  expect_identical(relabelled_audit$design$rank, audit$design$rank)
  expect_identical(permuted_audit$design$nullity, audit$design$nullity)
  expect_identical(relabelled_audit$design$nullity, audit$design$nullity)
})

test_that("QR tolerance sensitivity is diagnostic until a weak-information rule is frozen", {
  design <- Matrix::sparseMatrix(
    i = c(1L, 1L, 2L),
    j = c(1L, 2L, 2L),
    x = c(1, 1, 1e-9),
    dims = c(2L, 2L)
  )
  map <- mfrmr:::mfrmr_estimability_map(
    Block = c("facet", "facet"),
    Coordinate = c("facet:a", "facet:b"),
    Facet = c("facet", "facet"),
    Level = c("a", "b"),
    ReferenceLevel = c("", ""),
    Constraint = c("test", "test"),
    OptimizerIndex = c(1L, 2L)
  )
  result <- mfrmr:::mfrmr_estimability_rank_audit(design, map)

  expect_identical(result$State, "identified")
  expect_true(result$ToleranceSensitive)
  expect_gt(length(unique(result$ToleranceRanks$Rank)), 1L)
})

test_that("sparse constraint Jacobian matches the optimizer expansion", {
  build <- mfrmr:::build_facet_constraint
  cases <- list(
    build(
      levels = c("A", "B", "C", "D"),
      centered = TRUE
    ),
    build(
      levels = c("A", "B", "C", "D"),
      anchors = c(A = 0),
      centered = TRUE
    ),
    build(
      levels = c("A", "B", "C", "D"),
      groups = c(A = "G1", B = "G1", C = "G2", D = "G2"),
      group_values = c(G1 = 0, G2 = 0.5),
      centered = TRUE
    ),
    build(
      levels = c("A", "B", "C", "D"),
      anchors = c(A = -0.2),
      centered = FALSE
    )
  )

  for (spec in cases) {
    sparse <- mfrmr:::mfrmr_constraint_jacobian_sparse(spec, "Facet")
    dense <- mfrmr:::constraint_jacobian(spec)
    expect_equal(unname(as.matrix(sparse$jacobian)), unname(dense),
                 tolerance = 0)
    expect_equal(ncol(sparse$jacobian), spec$n_params)
    expect_equal(nrow(sparse$map), spec$n_params)
  }
})

test_that("PCM step coordinates enter the constrained rank audit", {
  data <- .estimability_balanced_data()
  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    model = "PCM",
    method = "JML",
    step_facet = "Criterion",
    maxit = 5
  ))
  audit <- fit$data_review$estimability
  step_block <- audit$parameter_blocks[
    audit$parameter_blocks$Block == "steps", , drop = FALSE
  ]
  expect_identical(audit$readiness$EstimabilityState, "identified")
  expect_identical(step_block$FreeCoordinates, 4L)
  expect_true(all(grepl(
    "^steps:",
    audit$parameter_map$Coordinate[audit$parameter_map$Block == "steps"]
  )))
})

test_that("zero-common-Person panels distinguish JML alias from MML linkage", {
  data <- .estimability_zero_common_data()

  condition <- tryCatch(
    suppressWarnings(.estimability_fit(data, "JML")),
    error = identity
  )
  expect_s3_class(condition, "mfrmr_estimability_error")
  expect_s3_class(condition, "mfrmr_readiness_condition")
  expect_identical(
    condition$readiness$EstimabilityState,
    "structurally_unidentified"
  )
  expect_identical(condition$estimability$design$nullity, 1L)
  expect_match(condition$readiness$ReasonCodes, "design_rank_deficient",
               fixed = TRUE)
  expect_match(conditionMessage(condition), "Optimization was not run",
               fixed = TRUE)
  expect_false(grepl("P01|P05", conditionMessage(condition)))

  mml <- NULL
  expect_warning(
    mml <- .estimability_fit(data, "MML"),
    class = "mfrmr_estimability_warning"
  )
  audit <- mml$data_review$estimability
  expect_identical(
    audit$readiness$EstimabilityState,
    "population_assumption_linked"
  )
  expect_true(audit$readiness$PopulationAssumptionLinked)
  expect_identical(audit$counterfactual_jml$State,
                   "structurally_unidentified")
  expect_identical(
    mml$data_review$status$Status[
      mml$data_review$status$Domain == "Design"
    ],
    "review_population_assumption_linked"
  )
})

test_that("a defensible anchor removes the zero-common-Person JML alias", {
  data <- .estimability_zero_common_data()
  anchors <- data.frame(
    Facet = "Rater",
    Level = "R1",
    Anchor = 0,
    stringsAsFactors = FALSE
  )

  fit <- suppressWarnings(.estimability_fit(
    data,
    method = "JML",
    anchors = anchors
  ))
  audit <- fit$data_review$estimability
  expect_identical(audit$readiness$EstimabilityState, "identified")
  expect_identical(audit$design$nullity, 0L)
  expect_identical(audit$design$rank, audit$design$free_dimension)
  expect_equal(
    unname(fit$config$anchor_summary$AnchoredLevels[
      fit$config$anchor_summary$Facet == "Rater"
    ]),
    1L
  )

  alternative <- data.frame(
    Facet = "Rater",
    Level = "R2",
    Anchor = 0,
    stringsAsFactors = FALSE
  )
  alternative_fit <- suppressWarnings(.estimability_fit(
    data,
    method = "JML",
    anchors = alternative
  ))
  alternative_audit <- alternative_fit$data_review$estimability
  expect_identical(alternative_audit$readiness$EstimabilityState,
                   "identified")
  expect_identical(alternative_audit$design$rank, audit$design$rank)
  expect_identical(alternative_audit$design$nullity, audit$design$nullity)
})

test_that("an unsupported interaction direction is blocked before optimization", {
  data <- .estimability_balanced_data()
  data <- data[!(data$Rater == "R2" & data$Criterion == "C2"), , drop = FALSE]

  condition <- tryCatch(
    suppressWarnings(.estimability_fit(
      data,
      method = "JML",
      facet_interactions = "Rater:Criterion"
    )),
    error = identity
  )
  expect_s3_class(condition, "mfrmr_estimability_error")
  expect_identical(condition$estimability$design$nullity, 1L)
  expect_true("interactions" %in% condition$estimability$null_blocks$Block)
  expect_identical(condition$readiness$EstimabilityState,
                   "structurally_unidentified")
})

test_that("GPCM additive rank does not overclaim nonlinear completeness", {
  data <- .estimability_balanced_data()
  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    model = "GPCM",
    method = "MML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 5,
    maxit = 50
  ))
  audit <- fit$data_review$estimability
  expect_identical(audit$design$state, "identified")
  expect_false(audit$readiness$Complete)
  expect_true("log_slopes" %in% audit$nonlinear_blocks)
  expect_match(audit$readiness$ReasonCodes, "design_rank_not_evaluated",
               fixed = TRUE)
  expect_identical(
    audit$fitted_information$status,
    "evaluated_diagnostic_only"
  )
  expect_identical(
    audit$nonlinear_transformation$status,
    "evaluated_diagnostic_only"
  )
  expect_identical(
    sort(audit$nonlinear_transformation$block_summary$Block),
    c("log_sigma2", "log_slopes")
  )
  expect_true(
    all(audit$nonlinear_transformation$block_summary$LogFullColumnRank)
  )
  expect_true(
    all(audit$nonlinear_transformation$block_summary$NaturalFullColumnRank)
  )
  expect_lt(
    max(audit$nonlinear_transformation$block_summary$
      MaxAbsLogJacobianDifference),
    1e-8
  )
  expect_lt(
    max(audit$nonlinear_transformation$block_summary$
      MaxScaledNaturalJacobianDifference),
    1e-8
  )
  expect_true(audit$nonlinear_transformation$parameterization_only)
  expect_false(
    audit$nonlinear_transformation$likelihood_jacobian_evaluated
  )
  expect_false(
    audit$nonlinear_transformation$structural_identification_classified
  )
  expect_true(audit$fitted_information$attempted)
  expect_equal(nrow(audit$fitted_information$rank_ladder), 3L)
  expect_identical(
    audit$fitted_information$numerical_differentiation$free_coordinate_step,
    1e-3
  )
  expect_true(all(is.finite(
    unlist(audit$fitted_information$evaluation_summary[1, ], use.names = FALSE)
  )))
  expect_equal(
    audit$fitted_information$evaluation_summary$ObjectiveDifference,
    0,
    tolerance = 1e-10
  )
  expect_identical(
    sort(audit$fitted_information$block_summary$Block),
    c("log_sigma2", "log_slopes")
  )
  expect_false(audit$fitted_information$weak_information_classified)
  expect_identical(
    audit$fitted_information$readiness_effect,
    "none_pending_pilot_calibrated_rule"
  )
  expect_identical(
    audit$gpcm_response_kernel$status,
    "not_evaluated_marginal_person_pattern_required"
  )
  expect_false(
    audit$gpcm_response_kernel$
      conditional_response_kernel_jacobian_evaluated
  )
  expect_false(
    audit$gpcm_response_kernel$
      marginal_person_pattern_jacobian_evaluated
  )
  expect_match(
    audit$gpcm_response_kernel$detail,
    "cannot certify the marginal person-pattern",
    fixed = TRUE
  )

  pattern <- audit$mml_observed_pattern_score
  expect_identical(
    pattern$status,
    "evaluated_observed_pattern_diagnostic_only"
  )
  expect_true(pattern$attempted)
  expect_identical(pattern$person_rows, 8L)
  expect_identical(pattern$execution_limits$persons, 200L)
  expect_identical(pattern$execution_limits$numeric_persons, 100L)
  expect_identical(pattern$free_dimension,
                   audit$readiness$OptimizerFreeDimension)
  expect_true(pattern$observed_patterns_only)
  expect_false(pattern$all_possible_response_patterns_evaluated)
  expect_false(pattern$conditional_jml_kernel_reused)
  expect_false(pattern$structural_identification_classified)
  expect_false(pattern$weak_information_classified)
  expect_identical(
    pattern$readiness_effect,
    "none_observed_patterns_diagnostic_only"
  )
  expect_equal(
    pattern$reconstruction$NegativeLogLikelihoodDifference,
    0,
    tolerance = 1e-12
  )
  expect_lt(
    pattern$reconstruction$ScoreSumVsNegativeGradientMaxAbs,
    1e-12
  )
  expect_identical(pattern$numerical_differentiation$status, "evaluated")
  expect_lt(pattern$numerical_differentiation$max_abs_difference, 1e-7)
  expect_lt(pattern$numerical_differentiation$max_scaled_difference, 1e-7)
  expect_true("log_slopes" %in% pattern$parameter_blocks$Block)
  expect_false("person_ids" %in% names(pattern))
  expect_false("score" %in% names(pattern))

  all_pattern <- audit$mml_all_pattern_information
  expect_identical(
    all_pattern$status,
    "evaluated_all_patterns_local_diagnostic_only"
  )
  expect_true(all_pattern$attempted)
  expect_identical(all_pattern$person_designs, 8L)
  expect_identical(all_pattern$retained_observation_rows, nrow(data))
  expect_equal(all_pattern$total_response_patterns, 8 * 4^4,
               tolerance = 0)
  expect_equal(all_pattern$evaluated_response_patterns, 4^4,
               tolerance = 0)
  expect_identical(all_pattern$unique_person_designs, 1L)
  expect_identical(all_pattern$design_reuse$PersonDesigns, 8L)
  expect_identical(all_pattern$design_reuse$UniquePersonDesigns, 1L)
  expect_identical(all_pattern$design_reuse$ReusedPersonDesigns, 7L)
  expect_identical(all_pattern$design_reuse$LargestDesignGroup, 8L)
  expect_equal(
    all_pattern$design_reuse$ConceptualToEvaluatedPatternRatio,
    8,
    tolerance = 0
  )
  expect_true(all_pattern$unit_row_weights_observed)
  expect_true(all_pattern$retained_observation_designs)
  expect_false(all_pattern$missing_rows_imputed)
  expect_false(all_pattern$observed_patterns_only)
  expect_true(all_pattern$all_possible_response_patterns_evaluated)
  expect_true(all_pattern$expected_information_evaluated)
  expect_false(all_pattern$conditional_jml_kernel_reused)
  expect_identical(all_pattern$local_rank, all_pattern$free_dimension)
  expect_identical(all_pattern$local_nullity, 0L)
  expect_lt(
    all_pattern$probability_normalization$MaximumAbsoluteMassError,
    1e-12
  )
  expect_lt(
    all_pattern$expected_score_identity$MaximumPersonScoreAbs,
    1e-12
  )
  expect_lt(
    all_pattern$expected_score_identity$MaximumAggregateScoreAbs,
    1e-11
  )
  expect_gte(
    all_pattern$expected_information_summary$MinimumEigenvalue,
    -1e-10
  )
  expect_lt(
    all_pattern$expected_information_summary$MaximumAsymmetry,
    1e-12
  )
  expect_identical(all_pattern$numerical_differentiation$status,
                   "evaluated")
  expect_lt(all_pattern$numerical_differentiation$max_abs_difference,
            1e-7)
  expect_lt(all_pattern$numerical_differentiation$max_scaled_difference,
            1e-7)
  expect_false(all_pattern$structural_identification_classified)
  expect_false(all_pattern$weak_information_classified)
  expect_identical(
    all_pattern$readiness_effect,
    "none_all_patterns_local_diagnostic_only"
  )
  expect_false("person_ids" %in% names(all_pattern))
  expect_false("weighted_score" %in% names(all_pattern))
  expect_false("expected_information" %in% names(all_pattern))
  expect_false("design_signatures" %in% names(all_pattern))
  expect_false("population_design_rows" %in% names(all_pattern))

  local <- audit$nonlinear_local_estimability
  expect_identical(local$contract_version,
                   "mfrmr-nonlinear-local-estimability-0.2.3-v1")
  expect_identical(local$state, "locally_full_rank_sufficient")
  expect_identical(
    local$evidence_basis,
    "all_pattern_expected_score_information_fixed_quadrature"
  )
  expect_true(local$local_first_order_classified)
  expect_true(local$local_full_rank_sufficient)
  expect_false(local$global_identification_classified)
  expect_false(local$continuous_integral_identification_classified)
  expect_false(local$weak_information_classified)
  expect_false(local$boundary_classified)
  expect_identical(local$readiness_effect, "none_local_property_only")
  expect_false(audit$readiness$Complete)

  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = fit$config$step_facet,
    slope_facet = fit$config$slope_facet,
    interaction_specs = fit$config$interaction_specs
  )
  sizes <- mfrmr:::build_param_sizes(fit$config)
  quad <- mfrmr:::gauss_hermite_normal(5L)
  person_score <- mfrmr:::mfrmr_mml_observed_person_score_matrix(
    fit$opt$par, idx, fit$config, sizes, quad
  )
  full_gradient <- mfrmr:::mfrm_grad_mml(
    fit$opt$par, idx, fit$config, sizes, quad
  )
  expect_lt(
    max(abs(unname(colSums(person_score$score)) - unname(-full_gradient))),
    1e-12
  )

  expected_information <-
    mfrmr:::mfrmr_mml_all_pattern_expected_information(
      fit$opt$par, idx, fit$config, sizes, quad
    )
  expect_lt(
    max(abs(
      crossprod(expected_information$weighted_score) -
        expected_information$expected_information
    )),
    1e-12
  )
  expect_lt(max(abs(expected_information$probability_mass - 1)), 1e-12)
  expect_lt(max(abs(expected_information$expected_score)), 1e-12)
  expect_equal(expected_information$total_patterns, 8 * 4^4,
               tolerance = 0)
  expect_equal(expected_information$evaluated_patterns, 4^4,
               tolerance = 0)
  expect_identical(expected_information$unique_person_designs, 1L)
  expect_identical(expected_information$design_group_sizes, 8L)

  unreused_expected_information <-
    mfrmr:::mfrmr_mml_all_pattern_expected_information(
      fit$opt$par, idx, fit$config, sizes, quad,
      reuse_identical_designs = FALSE
    )
  expect_equal(unreused_expected_information$evaluated_patterns,
               expected_information$total_patterns, tolerance = 0)
  expect_identical(unreused_expected_information$unique_person_designs,
                   8L)
  expect_lt(
    max(abs(
      expected_information$expected_information -
        unreused_expected_information$expected_information
    )),
    1e-12
  )
  expect_lt(
    max(abs(
      expected_information$expected_score -
        unreused_expected_information$expected_score
    )),
    1e-12
  )
  expect_lt(
    max(abs(
      expected_information$probability_mass -
        unreused_expected_information$probability_mass
    )),
    1e-12
  )

  zero_information <-
    mfrmr:::mfrmr_mml_all_pattern_expected_information(
      rep(0, length(fit$opt$par)), idx, fit$config, sizes,
      mfrmr:::gauss_hermite_normal(1L)
    )
  zero_map <- mfrmr:::mfrmr_mml_optimizer_parameter_map(
    fit$prep, idx, fit$config, sizes
  )
  zero_rank <- suppressWarnings(
    mfrmr:::mfrmr_estimability_rank_audit(
      Matrix::Matrix(zero_information$weighted_score, sparse = TRUE),
      zero_map
    )
  )
  expect_lt(zero_rank$Rank, nrow(zero_map))
  expect_true(any(grepl("log_slopes", zero_rank$ZeroCoordinates,
                        fixed = TRUE)))
  expect_lt(max(abs(zero_information$expected_score)), 1e-12)

  permuted_idx <- mfrmr:::mfrmr_subset_observation_indices(
    idx, rev(seq_along(idx$score_k))
  )
  permuted_score <- mfrmr:::mfrmr_mml_observed_person_score_matrix(
    fit$opt$par, permuted_idx, fit$config, sizes, quad
  )
  params <- mfrmr:::expand_params(fit$opt$par, sizes, fit$config)
  eap <- mfrmr:::compute_person_eap(idx, fit$config, params, quad)
  permuted_eap <- mfrmr:::compute_person_eap(
    permuted_idx, fit$config, params, quad
  )
  expect_equal(permuted_eap$Estimate, eap$Estimate, tolerance = 1e-12)
  expect_equal(permuted_eap$SD, eap$SD, tolerance = 1e-12)
  person_match <- match(person_score$person_ids, permuted_score$person_ids)
  expect_false(anyNA(person_match))
  expect_lt(
    max(abs(
      person_score$score -
        permuted_score$score[person_match, , drop = FALSE]
    )),
    1e-10
  )
  expect_lt(
    max(abs(
      person_score$log_marginal -
        permuted_score$log_marginal[person_match]
    )),
    1e-12
  )
  permuted_expected_information <-
    mfrmr:::mfrmr_mml_all_pattern_expected_information(
      fit$opt$par, permuted_idx, fit$config, sizes, quad
    )
  expect_lt(
    max(abs(
      expected_information$expected_information -
        permuted_expected_information$expected_information
    )),
    1e-10
  )
  expect_identical(permuted_expected_information$unique_person_designs, 1L)
  expect_equal(permuted_expected_information$evaluated_patterns, 4^4,
               tolerance = 0)

  mixed_rows <- unlist(lapply(seq_len(fit$config$n_person), function(id) {
    rows <- which(as.integer(idx$person) == id)
    if (id %% 2L == 1L) rev(rows) else rows
  }), use.names = FALSE)
  mixed_idx <- mfrmr:::mfrmr_subset_observation_indices(idx, mixed_rows)
  mixed <- mfrmr:::audit_mfrm_mml_all_pattern_information(
    opt = fit$opt,
    prep = fit$prep,
    idx = mixed_idx,
    config = fit$config,
    sizes = sizes,
    quad_points = 5L,
    nonlinear_blocks = "log_slopes"
  )
  expect_identical(
    mixed$status,
    "evaluated_all_patterns_local_diagnostic_only"
  )
  expect_identical(mixed$unique_person_designs, 1L)
  expect_equal(mixed$evaluated_response_patterns, 4^4, tolerance = 0)
  expect_identical(mixed$numerical_differentiation$status, "evaluated")
  expect_lt(mixed$numerical_differentiation$max_scaled_difference, 1e-7)

  retained_rows <- seq_along(idx$score_k)[-1L]
  missing_idx <- mfrmr:::mfrmr_subset_observation_indices(idx, retained_rows)
  missing <- mfrmr:::audit_mfrm_mml_all_pattern_information(
    opt = fit$opt,
    prep = fit$prep,
    idx = missing_idx,
    config = fit$config,
    sizes = sizes,
    quad_points = 5L,
    nonlinear_blocks = "log_slopes"
  )
  expect_identical(
    missing$status,
    "evaluated_all_patterns_local_diagnostic_only"
  )
  expect_identical(
    missing$retained_observations_per_person$Minimum,
    3L
  )
  expect_identical(
    missing$retained_observations_per_person$Maximum,
    4L
  )
  expect_equal(missing$total_response_patterns, 7 * 4^4 + 4^3,
               tolerance = 0)
  expect_equal(missing$evaluated_response_patterns, 4^4 + 4^3,
               tolerance = 0)
  expect_identical(missing$unique_person_designs, 2L)
  expect_identical(missing$design_reuse$ReusedPersonDesigns, 6L)
  expect_identical(missing$design_reuse$LargestDesignGroup, 7L)
  expect_true(missing$retained_observation_designs)
  expect_false(missing$missing_rows_imputed)
  expect_lt(
    missing$probability_normalization$MaximumAbsoluteMassError,
    1e-12
  )

  weighted_idx <- idx
  weighted_idx$weight[1] <- 2
  weighted <- mfrmr:::audit_mfrm_mml_all_pattern_information(
    opt = fit$opt,
    prep = fit$prep,
    idx = weighted_idx,
    config = fit$config,
    sizes = sizes,
    quad_points = 5L,
    nonlinear_blocks = "log_slopes"
  )
  expect_identical(weighted$status,
                   "not_evaluated_nonunit_row_weights")
  expect_false(weighted$attempted)
  expect_false(weighted$unit_row_weights_observed)
  expect_false(weighted$structural_identification_classified)

  limited <- mfrmr:::audit_mfrm_mml_observed_pattern_score(
    opt = fit$opt,
    prep = fit$prep,
    idx = idx,
    config = fit$config,
    sizes = sizes,
    quad_points = 5L,
    nonlinear_blocks = "log_slopes",
    max_persons = 0L
  )
  expect_identical(limited$status, "not_evaluated_execution_limit")
  expect_false(limited$attempted)
  expect_false(limited$structural_identification_classified)
  expect_identical(limited$readiness_effect,
                   "none_observed_patterns_diagnostic_only")

  analytic_only <- mfrmr:::audit_mfrm_mml_observed_pattern_score(
    opt = fit$opt,
    prep = fit$prep,
    idx = idx,
    config = fit$config,
    sizes = sizes,
    quad_points = 5L,
    nonlinear_blocks = "log_slopes",
    max_numeric_persons = 0L
  )
  expect_identical(
    analytic_only$status,
    "evaluated_observed_pattern_diagnostic_only"
  )
  expect_true(analytic_only$attempted)
  expect_true(analytic_only$finite_parameter_vector_observed)
  expect_identical(
    analytic_only$numerical_differentiation$status,
    "not_evaluated_execution_limit"
  )

  all_pattern_limited <-
    mfrmr:::audit_mfrm_mml_all_pattern_information(
      opt = fit$opt,
      prep = fit$prep,
      idx = idx,
      config = fit$config,
      sizes = sizes,
      quad_points = 5L,
      nonlinear_blocks = "log_slopes",
      max_total_patterns = 0
    )
  expect_identical(all_pattern_limited$status,
                   "not_evaluated_execution_limit")
  expect_false(all_pattern_limited$attempted)
  expect_false(
    all_pattern_limited$all_possible_response_patterns_evaluated
  )
  expect_false(all_pattern_limited$structural_identification_classified)
  expect_identical(
    all_pattern_limited$readiness_effect,
    "none_all_patterns_local_diagnostic_only"
  )

  compressed_limit <-
    mfrmr:::audit_mfrm_mml_all_pattern_information(
      opt = fit$opt,
      prep = fit$prep,
      idx = idx,
      config = fit$config,
      sizes = sizes,
      quad_points = 5L,
      nonlinear_blocks = "log_slopes",
      max_total_patterns = 4^4
    )
  expect_identical(
    compressed_limit$status,
    "evaluated_all_patterns_local_diagnostic_only"
  )
  expect_equal(compressed_limit$total_response_patterns, 8 * 4^4,
               tolerance = 0)
  expect_equal(compressed_limit$evaluated_response_patterns, 4^4,
               tolerance = 0)
})

test_that("JML GPCM response-kernel Jacobian joins additive and slope coordinates", {
  data <- .estimability_balanced_data()
  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    maxit = 100
  ))
  estimability <- fit$data_review$estimability
  audit <- estimability$gpcm_response_kernel

  expect_identical(audit$status, "evaluated_local_diagnostic_only")
  expect_true(audit$attempted)
  expect_true(audit$conditional_response_kernel_jacobian_evaluated)
  expect_false(audit$marginal_person_pattern_jacobian_evaluated)
  expect_false(audit$structural_identification_classified)
  expect_identical(audit$local_rank_state, "locally_full_column_rank")
  expect_identical(audit$local_rank, audit$free_dimension)
  expect_identical(audit$local_nullity, 0L)
  expect_identical(audit$observation_rows, nrow(data))
  expect_identical(audit$transition_rows, 3L * nrow(data))
  expect_identical(
    audit$parameter_map$OptimizerIndex,
    seq_len(audit$free_dimension)
  )
  expect_true("log_slopes" %in% audit$parameter_blocks$Block)
  expect_identical(audit$numerical_differentiation$status, "evaluated")
  expect_lt(audit$numerical_differentiation$max_abs_difference, 1e-7)
  expect_lt(audit$numerical_differentiation$max_scaled_difference, 1e-7)
  expect_identical(
    audit$readiness_effect,
    "none_pending_property_and_marginal_model_audits"
  )
  local <- estimability$nonlinear_local_estimability
  expect_identical(local$state, "locally_full_rank_sufficient")
  expect_identical(
    local$evidence_basis,
    "conditional_adjacent_logit_jacobian_full_free_coordinates"
  )
  expect_true(local$local_first_order_classified)
  expect_true(local$local_full_rank_sufficient)
  expect_false(local$global_identification_classified)
  expect_false(local$weak_information_classified)
  expect_false(local$boundary_classified)
  expect_identical(local$readiness_effect, "none_local_property_only")
  expect_false(fit$readiness$fit$InferenceReady)

  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = fit$config$step_facet,
    slope_facet = fit$config$slope_facet,
    interaction_specs = fit$config$interaction_specs
  )
  sizes <- mfrmr:::build_param_sizes(fit$config)
  slices <- mfrmr:::build_param_slices(sizes)
  unit_par <- fit$opt$par
  unit_par[slices$log_slopes] <- 0
  combined <- mfrmr:::mfrmr_gpcm_response_kernel_design(
    fit$prep, idx, fit$config, sizes, unit_par
  )
  additive <- mfrmr:::mfrmr_estimability_adjacent_design(
    fit$prep, idx, fit$config, sizes,
    include_person = TRUE,
    include_population_beta = FALSE
  )
  additive_columns <- seq_len(ncol(additive$design))
  expect_equal(
    unname(as.matrix(combined$jacobian[, additive_columns, drop = FALSE])),
    unname(as.matrix(additive$design)),
    tolerance = 0
  )
  expect_equal(combined$adjacent_logits, combined$eta_minus_step,
               tolerance = 1e-12)

  limited <- mfrmr:::audit_mfrm_gpcm_response_kernel(
    prep = fit$prep,
    idx = idx,
    config = fit$config,
    sizes = sizes,
    par = fit$opt$par,
    max_numeric_free_dimension = 0L
  )
  expect_identical(limited$status, "evaluated_local_diagnostic_only")
  expect_identical(
    limited$numerical_differentiation$status,
    "not_evaluated_execution_limit"
  )
  expect_true(limited$conditional_response_kernel_jacobian_evaluated)
  expect_false(limited$structural_identification_classified)
})

test_that("GPCM selected log-slope expansion is sparse and exact", {
  slope_index <- c(1L, 4L, 2L, 4L, 3L)
  sparse <- mfrmr:::mfrmr_gpcm_selected_log_slope_jacobian_sparse(
    slope_index = slope_index,
    n_levels = 4L
  )
  dense <- mfrmr:::sum_zero_jacobian(4L)[slope_index, , drop = FALSE]

  expect_s4_class(sparse, "dgCMatrix")
  expect_equal(unname(as.matrix(sparse)), unname(dense), tolerance = 0)
  expect_identical(nrow(sparse), length(slope_index))
  expect_identical(ncol(sparse), 3L)
})

test_that("JML GPCM response-kernel rank is invariant to row permutation", {
  data <- .estimability_balanced_data()
  permuted <- data[c(seq(2L, nrow(data), by = 2L),
                     seq(1L, nrow(data), by = 2L)), , drop = FALSE]
  fit_one <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    maxit = 100
  ))
  fit_two <- suppressWarnings(fit_mfrm(
    permuted,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    maxit = 100
  ))
  audit_one <- fit_one$data_review$estimability$gpcm_response_kernel
  audit_two <- fit_two$data_review$estimability$gpcm_response_kernel

  expect_identical(audit_two$status, audit_one$status)
  expect_identical(audit_two$free_dimension, audit_one$free_dimension)
  expect_identical(audit_two$local_rank, audit_one$local_rank)
  expect_identical(audit_two$local_nullity, audit_one$local_nullity)
  expect_identical(
    audit_two$rank_ladder[, c("Rank", "Nullity")],
    audit_one$rank_ladder[, c("Rank", "Nullity")]
  )
  expect_lt(audit_two$numerical_differentiation$max_scaled_difference,
            1e-7)
})

test_that("latent residual variance enters fitted-information instrumentation", {
  set.seed(2718)
  n_person <- 50L
  n_item <- 5L
  persons <- paste0("P", sprintf("%03d", seq_len(n_person)))
  items <- paste0("I", seq_len(n_item))
  x <- stats::rnorm(n_person)
  theta <- 0.25 + 0.9 * x + stats::rnorm(n_person, sd = 0.6)
  item_beta <- seq(-1.2, 1.2, length.out = n_item)
  data <- expand.grid(
    Person = persons,
    Item = items,
    stringsAsFactors = FALSE
  )
  eta <- theta[match(data$Person, persons)] -
    item_beta[match(data$Item, items)]
  data$Score <- stats::rbinom(nrow(data), 1, stats::plogis(eta))
  person_data <- data.frame(Person = persons, X = x)

  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = "Item",
    score = "Score",
    model = "RSM",
    method = "MML",
    population_formula = ~ X,
    person_data = person_data,
    quad_points = 7,
    maxit = 150
  ))
  audit <- fit$data_review$estimability

  expect_identical(audit$nonlinear_blocks, "log_sigma2")
  expect_identical(
    audit$fitted_information$status,
    "evaluated_diagnostic_only"
  )
  expect_identical(
    audit$fitted_information$block_summary$Block,
    "log_sigma2"
  )
  expect_identical(
    audit$nonlinear_transformation$status,
    "evaluated_diagnostic_only"
  )
  expect_identical(
    audit$nonlinear_transformation$block_summary$Block,
    "log_sigma2"
  )
  expect_equal(
    audit$nonlinear_transformation$block_summary$NaturalMinimum,
    fit$population$sigma2,
    tolerance = 1e-12
  )
  expect_identical(
    audit$fitted_information$block_summary$FreeCoordinates,
    1L
  )
  expect_true(all(is.finite(
    unlist(audit$fitted_information$evaluation_summary[1, ], use.names = FALSE)
  )))
  expect_false(audit$fitted_information$weak_information_classified)
  expect_false(audit$readiness$Complete)
  pattern <- audit$mml_observed_pattern_score
  expect_identical(
    pattern$status,
    "evaluated_observed_pattern_diagnostic_only"
  )
  expect_true(all(c("beta", "log_sigma2") %in%
                    pattern$parameter_blocks$Block))
  expect_lt(
    pattern$reconstruction$ScoreSumVsNegativeGradientMaxAbs,
    1e-10
  )
  expect_lt(pattern$numerical_differentiation$max_scaled_difference, 1e-7)
  expect_false(pattern$structural_identification_classified)
  all_pattern <- audit$mml_all_pattern_information
  expect_identical(
    all_pattern$status,
    "evaluated_all_patterns_local_diagnostic_only"
  )
  expect_true(all(c("beta", "log_sigma2") %in%
                    all_pattern$parameter_blocks$Block))
  expect_equal(all_pattern$total_response_patterns, n_person * 2^n_item,
               tolerance = 0)
  expect_equal(all_pattern$evaluated_response_patterns,
               all_pattern$total_response_patterns, tolerance = 0)
  expect_identical(all_pattern$unique_person_designs, n_person)
  expect_identical(all_pattern$design_reuse$ReusedPersonDesigns, 0L)
  latent_idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = fit$config$step_facet,
    slope_facet = fit$config$slope_facet,
    interaction_specs = fit$config$interaction_specs
  )
  latent_groups <- mfrmr:::mfrmr_mml_person_design_groups(
    latent_idx, fit$config
  )
  expect_identical(latent_groups$unique_person_designs, n_person)
  repeated_config <- fit$config
  repeated_lookup <- repeated_config$population_spec$person_lookup
  repeated_config$population_spec$design_matrix[
    repeated_lookup[2L],
  ] <- repeated_config$population_spec$design_matrix[
    repeated_lookup[1L],
  ]
  repeated_groups <- mfrmr:::mfrmr_mml_person_design_groups(
    latent_idx, repeated_config
  )
  expect_identical(repeated_groups$unique_person_designs, n_person - 1L)
  expect_identical(max(repeated_groups$group_sizes), 2L)
  expect_lt(
    all_pattern$probability_normalization$MaximumAbsoluteMassError,
    1e-12
  )
  expect_lt(
    all_pattern$expected_score_identity$MaximumPersonScoreAbs,
    1e-10
  )
  expect_lt(all_pattern$numerical_differentiation$max_scaled_difference,
            1e-7)
  expect_false(all_pattern$structural_identification_classified)
})

test_that("nonlinear transformation Jacobians match known analytic maps", {
  free_slopes <- c(-0.4, 0.1, 0.2)
  free_log_sigma2 <- log(0.36)
  audit <- mfrmr:::mfrmr_nonlinear_transformation_audit(
    par = c(free_slopes, free_log_sigma2),
    sizes = list(log_slopes = 3L, log_sigma2 = 1L),
    config = list(
      method = "MML",
      model = "GPCM",
      gpcm_spec = list(
        levels = paste0("I", 1:4),
        n_params = 3L
      ),
      population_spec = list(active = TRUE)
    ),
    nonlinear_blocks = c("log_slopes", "log_sigma2")
  )

  expect_identical(audit$status, "evaluated_diagnostic_only")
  expect_identical(audit$block_summary$Block,
                   c("log_slopes", "log_sigma2"))
  expect_identical(audit$block_summary$ExpectedRank, c(3L, 1L))
  expect_true(all(audit$block_summary$LogFullColumnRank))
  expect_true(all(audit$block_summary$NaturalFullColumnRank))
  expect_lt(max(audit$block_summary$MaxAbsLogJacobianDifference), 1e-8)
  expect_lt(
    max(audit$block_summary$MaxScaledNaturalJacobianDifference),
    1e-8
  )
  expect_equal(
    audit$block_summary$InvariantResidual[audit$block_summary$Block ==
                                            "log_slopes"],
    0,
    tolerance = 1e-15
  )
  expect_equal(
    audit$block_summary$NaturalMinimum[audit$block_summary$Block ==
                                        "log_sigma2"],
    0.36,
    tolerance = 1e-12
  )
  expect_equal(nrow(audit$rank_ladder), 12L)
  expect_false(audit$structural_identification_classified)
  expect_identical(
    audit$readiness_effect,
    "none_parameterization_audit_only"
  )
})

test_that("nonlinear transformation audit fails closed on malformed vectors", {
  audit <- mfrmr:::mfrmr_nonlinear_transformation_audit(
    par = c(0, 0),
    sizes = list(log_slopes = 3L),
    config = list(method = "MML", model = "GPCM"),
    nonlinear_blocks = "log_slopes"
  )

  expect_identical(audit$status, "not_evaluated_parameter_vector")
  expect_false(audit$attempted)
  expect_false(audit$structural_identification_classified)
})

test_that("underflowed natural coordinates remain unavailable and nonbinding", {
  audit <- mfrmr:::mfrmr_nonlinear_transformation_audit(
    par = -1000,
    sizes = list(log_sigma2 = 1L),
    config = list(
      method = "MML",
      model = "RSM",
      population_spec = list(active = TRUE)
    ),
    nonlinear_blocks = "log_sigma2"
  )

  expect_identical(audit$status, "partially_unavailable")
  expect_identical(audit$block_summary$Status, "unavailable")
  expect_false(audit$block_summary$Finite)
  expect_false(audit$structural_identification_classified)
  expect_identical(audit$readiness_effect,
                   "none_parameterization_audit_only")
})

test_that("fitted-information eigenvalue ladder remains diagnostic", {
  result <- mfrmr:::mfrmr_information_rank_ladder(
    diag(c(4, 1, 1e-9, -1e-7)),
    tolerances = c(1e-10, 1e-8, 1e-6)
  )

  expect_true(result$valid)
  expect_equal(nrow(result$rank_ladder), 3L)
  expect_equal(result$rank_ladder$PositiveRank, c(3L, 2L, 2L))
  expect_equal(result$rank_ladder$NegativeCount, c(1L, 1L, 0L))
  expect_equal(result$rank_ladder$NearZeroCount, c(0L, 1L, 2L))
  expect_equal(result$smallest_eigenvalue, -1e-7)
  expect_equal(result$largest_eigenvalue, 4)
})

test_that("fitted-information execution limit is not a readiness classification", {
  opt <- list(
    par = rep(0, 81L),
    optimizer_diagnostics = list(ConvergenceSeverity = "pass")
  )
  audit <- mfrmr:::audit_mfrm_fitted_information(
    opt = opt,
    idx = NULL,
    config = list(method = "MML", model = "GPCM"),
    sizes = list(log_slopes = 81L),
    quad_points = 5L,
    nonlinear_blocks = "log_slopes",
    max_free_dimension = 80L
  )

  expect_identical(audit$status, "not_evaluated_dimension_limit")
  expect_false(audit$attempted)
  expect_false(audit$weak_information_classified)
  expect_identical(audit$readiness_effect,
                   "none_pending_pilot_calibrated_rule")
  expect_match(audit$detail, "execution limit, not an estimability classification",
               fixed = TRUE)
})

test_that("nonstationary retained states do not trigger fitted information", {
  opt <- list(
    par = rep(0, 2L),
    optimizer_diagnostics = list(ConvergenceSeverity = "review")
  )
  audit <- mfrmr:::audit_mfrm_fitted_information(
    opt = opt,
    idx = NULL,
    config = list(method = "MML", model = "GPCM"),
    sizes = list(log_slopes = 2L),
    quad_points = 5L,
    nonlinear_blocks = "log_slopes"
  )

  expect_identical(audit$status, "not_evaluated_nonstationary")
  expect_false(audit$attempted)
  expect_false(audit$weak_information_classified)
})
