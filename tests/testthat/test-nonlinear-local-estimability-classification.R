.nonlinear_local_parameter_map <- function(free_dimension) {
  data.frame(
    Block = rep("log_slopes", free_dimension),
    Coordinate = paste0("log_slopes:", seq_len(free_dimension)),
    Facet = rep("Criterion", free_dimension),
    Level = as.character(seq_len(free_dimension)),
    ReferenceLevel = rep("reference", free_dimension),
    Constraint = rep("sum_zero_log_slopes", free_dimension),
    OptimizerIndex = seq_len(free_dimension),
    stringsAsFactors = FALSE
  )
}

.nonlinear_local_mml_observed <- function(rank = 2L, nullity = 0L,
                                           tolerance_sensitive = FALSE) {
  list(
    status = "evaluated_observed_pattern_diagnostic_only",
    unit_row_weights_observed = TRUE,
    finite_parameter_vector_observed = TRUE,
    score_rank = rank,
    score_nullity = nullity,
    free_dimension = rank + nullity,
    tolerance_sensitive = tolerance_sensitive,
    parameter_map = .nonlinear_local_parameter_map(rank + nullity),
    quadrature_points = 31L,
    numerical_differentiation = list(status = "evaluated")
  )
}

.nonlinear_local_mml_all_pattern <- function(
    rank = 2L, nullity = 0L, tolerance_sensitive = FALSE) {
  list(
    status = "evaluated_all_patterns_local_diagnostic_only",
    all_possible_response_patterns_evaluated = TRUE,
    expected_information_evaluated = TRUE,
    unit_row_weights_observed = TRUE,
    local_rank = rank,
    local_nullity = nullity,
    free_dimension = rank + nullity,
    tolerance_sensitive = tolerance_sensitive,
    parameter_map = .nonlinear_local_parameter_map(rank + nullity),
    quadrature_points = 31L,
    numerical_differentiation = list(status = "evaluated")
  )
}

test_that("MML observed-pattern span is a sufficient local certificate", {
  audit <- list(
    nonlinear_blocks = c("log_slopes", "log_sigma2"),
    mml_observed_pattern_score = .nonlinear_local_mml_observed(),
    mml_all_pattern_information = list(status = "not_evaluated_execution_limit")
  )
  classified <- mfrmr:::mfrmr_classify_nonlinear_local_estimability(
    audit, list(method = "MML", model = "GPCM")
  )

  expect_identical(classified$state, "locally_full_rank_sufficient")
  expect_identical(
    classified$evidence_basis,
    "observed_pattern_subset_score_span_fixed_quadrature"
  )
  expect_true(classified$local_first_order_classified)
  expect_true(classified$local_full_rank_sufficient)
  expect_false(classified$global_identification_classified)
  expect_false(classified$continuous_integral_identification_classified)
  expect_false(classified$weak_information_classified)
  expect_false(classified$boundary_classified)
  expect_identical(classified$readiness_effect,
                   "none_local_property_only")
  expect_identical(classified$quadrature_points, 31L)
})

test_that("a deficient observed MML subset is not a nonidentification claim", {
  audit <- list(
    nonlinear_blocks = "log_slopes",
    mml_observed_pattern_score = .nonlinear_local_mml_observed(
      rank = 1L, nullity = 1L
    ),
    mml_all_pattern_information = list(status = "not_evaluated_execution_limit")
  )
  classified <- mfrmr:::mfrmr_classify_nonlinear_local_estimability(
    audit, list(method = "MML", model = "GPCM")
  )

  expect_identical(classified$state, "not_evaluated")
  expect_false(classified$local_first_order_classified)
  expect_false(classified$local_first_order_rank_deficient)
  expect_false(classified$local_nonidentifiability_established)
  expect_match(
    classified$reason_codes,
    "observed_pattern_subset_rank_deficient_is_inconclusive",
    fixed = TRUE
  )
})

test_that("complete MML information can classify first-order deficiency", {
  audit <- list(
    nonlinear_blocks = "log_slopes",
    mml_observed_pattern_score = .nonlinear_local_mml_observed(
      rank = 1L, nullity = 1L
    ),
    mml_all_pattern_information = .nonlinear_local_mml_all_pattern(
      rank = 1L, nullity = 1L
    )
  )
  classified <- mfrmr:::mfrmr_classify_nonlinear_local_estimability(
    audit, list(method = "MML", model = "GPCM")
  )

  expect_identical(
    classified$state, "locally_first_order_rank_deficient"
  )
  expect_identical(
    classified$evidence_basis,
    "all_pattern_expected_score_information_fixed_quadrature"
  )
  expect_true(classified$local_first_order_classified)
  expect_true(classified$local_first_order_rank_deficient)
  expect_false(classified$local_nonidentifiability_established)
  expect_false(classified$global_identification_classified)
})

test_that("rank tolerance sensitivity remains indeterminate", {
  audit <- list(
    nonlinear_blocks = "log_slopes",
    gpcm_response_kernel = list(
      status = "evaluated_local_diagnostic_only",
      conditional_response_kernel_jacobian_evaluated = TRUE,
      local_rank = 2L,
      local_nullity = 0L,
      free_dimension = 2L,
      tolerance_sensitive = TRUE,
      parameter_map = .nonlinear_local_parameter_map(2L),
      numerical_differentiation = list(status = "evaluated")
    )
  )
  classified <- mfrmr:::mfrmr_classify_nonlinear_local_estimability(
    audit, list(method = "JML", model = "GPCM")
  )

  expect_identical(classified$state,
                   "indeterminate_tolerance_sensitive")
  expect_false(classified$local_first_order_classified)
  expect_false(classified$local_full_rank_sufficient)
  expect_false(classified$global_identification_classified)
})

test_that("malformed local parameter maps fail closed", {
  observed <- .nonlinear_local_mml_observed()
  observed$parameter_map$OptimizerIndex <- c(1L, 1L)
  audit <- list(
    nonlinear_blocks = "log_slopes",
    mml_observed_pattern_score = observed,
    mml_all_pattern_information = list(status = "not_evaluated_execution_limit")
  )
  classified <- mfrmr:::mfrmr_classify_nonlinear_local_estimability(
    audit, list(method = "MML", model = "GPCM")
  )

  expect_identical(classified$state, "not_evaluated")
  expect_false(classified$parameter_map_complete)
  expect_false(classified$local_full_rank_sufficient)
})

test_that("a fitted MML GPCM can use the observed-pattern certificate", {
  set.seed(20260812)
  data <- expand.grid(
    Person = paste0("P", sprintf("%02d", 1:20)),
    Rater = c("R1", "R2"),
    Criterion = c("C1", "C2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  person_index <- as.integer(factor(data$Person))
  rater_index <- as.integer(factor(data$Rater))
  criterion_index <- as.integer(factor(data$Criterion))
  eta <- (person_index - 10) / 7 +
    0.3 * (rater_index - 1) - 0.4 * (criterion_index - 1)
  cumulative <- plogis(outer(eta, c(-1.2, -0.2, 0.9), "-"))
  uniform <- stats::runif(nrow(data))
  data$Score <- 1L + rowSums(uniform > cumulative)

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
    maxit = 100
  ))
  observed <- fit$data_review$estimability$mml_observed_pattern_score
  classified <-
    fit$data_review$estimability$nonlinear_local_estimability

  expect_identical(observed$score_rank, observed$free_dimension)
  expect_identical(observed$score_nullity, 0L)
  expect_false(observed$tolerance_sensitive)
  expect_identical(classified$state, "locally_full_rank_sufficient")
  expect_identical(
    classified$evidence_basis,
    "observed_pattern_subset_score_span_fixed_quadrature"
  )
  expect_true(classified$local_full_rank_sufficient)
  expect_false(classified$global_identification_classified)
  expect_identical(classified$readiness_effect,
                   "none_local_property_only")
  expect_false(fit$readiness$fit$InferenceReady)
})
