make_jml_gpcm_general_rate_negative_fixture <- function() {
  data <- expand.grid(
    Person = c("P1", "P2"),
    Criterion = c("C1", "C2", "C3"),
    Score = 0:1,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data[rep(seq_len(nrow(data)), each = 2L), , drop = FALSE]
}

fit_jml_gpcm_general_rate_negative_fixture <- function() {
  suppressWarnings(fit_mfrm(
    make_jml_gpcm_general_rate_negative_fixture(),
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

fit_jml_gpcm_general_rate_pair_positive_fixture <- function() {
  data <- expand.grid(
    Person = c("P1", "P2"),
    Criterion = c("C1", "C2", "C3"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$Score <- c(0L, 1L, 1L, 0L, 0L, 0L)
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
    maxit = 300,
    reltol = 1e-12
  ))
}

general_rate_problem <- function(fit) {
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

test_that("canonical rate partitions extend the ordered-pair family", {
  expect_equal(
    vapply(2:4, mfrmr:::mfrmr_jml_gpcm_joint_rate_candidate_count,
           numeric(1)),
    c(2, 18, 110)
  )
  expect_identical(
    vapply(2:4, function(levels) {
      length(mfrmr:::mfrmr_jml_gpcm_joint_rate_candidates(levels))
    }, integer(1)),
    c(0L, 12L, 98L)
  )

  candidates <- mfrmr:::mfrmr_jml_gpcm_joint_rate_candidates(4L)
  complete_partition <- vapply(candidates, function(candidate) {
    groups <- c(
      candidate$positive, candidate$zero,
      candidate$leading_negative, candidate$deeper_negative
    )
    identical(sort(groups), 1:4)
  }, logical(1))
  expect_true(all(complete_partition))
  expect_equal(
    vapply(candidates, function(candidate) sum(candidate$rates), numeric(1)),
    rep(0, length(candidates)), tolerance = 1e-14
  )
  expect_true(all(vapply(candidates, function(candidate) {
    all(candidate$rates[candidate$positive] > 0)
  }, logical(1))))
  expect_true(all(vapply(candidates, function(candidate) {
    all(candidate$rates[candidate$leading_negative] == -1)
  }, logical(1))))
  expect_true(all(vapply(candidates, function(candidate) {
    all(candidate$rates[candidate$deeper_negative] == -2)
  }, logical(1))))
  three_level_rates <- vapply(
    mfrmr:::mfrmr_jml_gpcm_joint_rate_candidates(3L),
    function(candidate) paste(candidate$rates, collapse = ","),
    character(1)
  )
  expect_true("3,-1,-2" %in% three_level_rates)
})

test_that("a deeper negative group yields a certificate missed by every pair", {
  skip_if_not_installed("lpSolve")
  contrast <- Matrix::Matrix(
    matrix(c(1, -1, 1, -1), ncol = 1L), sparse = TRUE
  )
  contrast_observation <- 1:4
  slope_index <- c(1L, 2L, 3L, 3L)
  weight <- rep(1, 4)

  general <- mfrmr:::mfrmr_jml_gpcm_joint_rate_lp(
    contrast = contrast,
    contrast_observation = contrast_observation,
    slope_index = slope_index,
    weight = weight,
    n_categories = 2L,
    positive_indices = 1L,
    zero_indices = integer(0),
    leading_negative_indices = 2L,
    max_lp_nonzeros = 1e6,
    lp_timeout = 2L,
    certificate_tolerance = 1e-7
  )
  expect_true(general$evaluated)
  expect_true(general$certified)
  expect_gt(general$positive_minimum_margin, 0)
  expect_gt(general$negative_leading_coefficient, 0)

  pair_results <- unlist(lapply(1:3, function(positive) {
    vapply(setdiff(1:3, positive), function(negative) {
      mfrmr:::mfrmr_jml_gpcm_joint_pair_lp(
        contrast, contrast_observation, slope_index, weight, 2L,
        positive, negative, 1e6, 2L, 1e-7
      )$certified
    }, logical(1))
  }))
  expect_true(all(!pair_results))

  rates <- c(3, -1, -2)
  expect_equal(sum(rates), 0)
  utilities <- matrix(0, nrow = 4L, ncol = 2L)
  category_direction <- cbind(0, c(1, -1, 1, -1))
  boundary <- mfrmr:::mfrmr_jml_gpcm_joint_rate_boundary_limit(
    utilities = utilities,
    category_direction = category_direction,
    score_k = rep(1L, 4L),
    slope_index = slope_index,
    slopes = rep(1, 3L),
    weight = weight,
    positive_indices = 1L,
    zero_indices = integer(0),
    negative_indices = c(2L, 3L),
    tolerance = 1e-7
  )
  expect_true(boundary$valid)
  expect_equal(boundary$log_likelihood, -3 * log(2), tolerance = 1e-12)
  expect_equal(
    boundary$log_likelihood - (-4 * log(2)), log(2),
    tolerance = 1e-12
  )

  tail_t <- c(2, 4, 8, 16)
  tail_log_likelihood <- vapply(tail_t, function(distance) {
    adjacent <- distance * c(1, -1, 1, -1) *
      exp(rates[slope_index] * distance)
    sum(-log1p(exp(-adjacent)))
  }, numeric(1))
  expect_gte(min(diff(tail_log_likelihood)), -1e-12)
  expect_equal(
    tail_log_likelihood[length(tail_log_likelihood)],
    boundary$log_likelihood,
    tolerance = 1e-6
  )
})

test_that("three-level production fits complete the canonical rate scope", {
  skip_if_not_installed("lpSolve")
  fit <- fit_jml_gpcm_general_rate_negative_fixture()
  audit <- fit$config$boundary_audit$gpcm_joint_boundary
  classified <-
    fit$config$boundary_audit$gpcm_fixed_objective_classification

  expect_identical(audit$state, "none_certified")
  expect_true(audit$complete)
  expect_true(audit$scope_complete)
  expect_identical(audit$rate_scope_state, "none_certified")
  expect_true(audit$rate_scope_complete)
  expect_false(audit$positive_certificate_present)
  expect_identical(nrow(audit$certificates), 6L)
  expect_identical(nrow(audit$rate_certificates), 12L)
  expect_true(all(audit$rate_certificates$EvaluationState == "evaluated"))
  expect_true(all(!audit$rate_certificates$Certified))
  expect_equal(audit$dimensions$CanonicalRateCandidates, 18)
  expect_equal(audit$dimensions$AdditionalRateCandidates, 12)
  expect_equal(audit$dimensions$EvaluatedRateCandidates, 12)
  expect_identical(
    classified$contract_version,
    "mfrmr-jml-gpcm-fixed-objective-boundary-0.2.3-v2"
  )
  expect_true(classified$canonical_constant_rate_scope_complete)
  expect_false(classified$canonical_constant_rate_boundary_certified)
  expect_identical(
    classified$state, "finite_retained_point_no_path_certified"
  )
  expect_false(classified$global_boundary_absence_certified)
})

test_that("rate workload limits stay open while positive pairs stay primary", {
  skip_if_not_installed("lpSolve")
  negative_fit <- fit_jml_gpcm_general_rate_negative_fixture()
  problem <- general_rate_problem(negative_fit)
  limited <- mfrmr:::audit_mfrm_jml_gpcm_joint_boundary(
    negative_fit$prep, problem$idx, problem$config, problem$sizes,
    negative_fit$opt$par, max_rate_candidates = 6L
  )
  expect_identical(limited$state, "not_evaluated_rate_size_limit")
  expect_false(limited$complete)
  expect_false(limited$scope_complete)
  expect_identical(limited$rate_scope_state, "not_evaluated_size_limit")
  expect_false(limited$rate_scope_complete)
  expect_false(limited$positive_certificate_present)
  limited_classification <-
    mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
      problem$config,
      negative_fit$config$boundary_audit$gpcm_slope_boundary,
      limited
    )
  expect_identical(limited_classification$state, "not_evaluated")
  expect_false(
    limited_classification$canonical_constant_rate_scope_complete
  )
  expect_false(limited_classification$no_path_certified_in_audited_families)

  positive_fit <- fit_jml_gpcm_general_rate_pair_positive_fixture()
  positive_problem <- general_rate_problem(positive_fit)
  positive_limited <- mfrmr:::audit_mfrm_jml_gpcm_joint_boundary(
    positive_fit$prep,
    positive_problem$idx,
    positive_problem$config,
    positive_problem$sizes,
    positive_fit$opt$par,
    max_rate_candidates = 6L
  )
  expect_identical(
    positive_limited$state,
    "certified_competitive_joint_boundary_path"
  )
  expect_false(positive_limited$scope_complete)
  expect_true(positive_limited$positive_certificate_present)
  expect_true(any(positive_limited$certificates$Certified))

  classified <- mfrmr:::mfrmr_classify_jml_gpcm_fixed_objective_boundary(
    positive_problem$config,
    positive_fit$config$boundary_audit$gpcm_slope_boundary,
    positive_limited
  )
  expect_identical(
    classified$state, "certified_competitive_joint_boundary"
  )
  expect_true(classified$positive_boundary_certificate_present)
  expect_false(classified$path_family_classification_complete)
  expect_false(classified$global_boundary_classified)
})
