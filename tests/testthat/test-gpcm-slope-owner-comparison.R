test_that("GPCM assigns level-specific slopes to exactly one facet", {
  facets <- c("Rater", "Task", "Criterion")

  criterion <- mfrmr:::resolve_step_and_slope_facets(
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    facet_names = facets
  )
  expect_identical(criterion$step_facet, "Criterion")
  expect_identical(criterion$slope_facet, "Criterion")

  rater <- mfrmr:::resolve_step_and_slope_facets(
    model = "GPCM",
    step_facet = "Rater",
    slope_facet = "Rater",
    facet_names = facets
  )
  expect_identical(rater$step_facet, "Rater")
  expect_identical(rater$slope_facet, "Rater")

  expect_error(
    mfrmr:::resolve_step_and_slope_facets(
      model = "GPCM",
      step_facet = "Criterion",
      slope_facet = "Rater",
      facet_names = facets
    ),
    "slope_facet == step_facet",
    fixed = TRUE
  )

  criterion_levels <- c("Content", "Organization", "Language")
  expect_identical(
    mfrmr:::sum_zero_param_count(length(criterion_levels)),
    2L
  )
  log_slopes <- mfrmr:::expand_sum_zero_vector(c(log(1.2), log(0.8)), 3L)
  slopes <- exp(log_slopes)
  expect_length(slopes, length(criterion_levels))
  expect_true(all(slopes > 0))
  expect_equal(exp(mean(log(slopes))), 1, tolerance = 1e-15)
})

test_that("model-choice contract names level-specific single-owner slopes", {
  contract <- mfrmr:::.model_choice_score_contract("GPCM")
  expect_match(contract, "every level of one designated slope facet")
  expect_match(contract, "slope_facet == step_facet", fixed = TRUE)
  expect_match(contract, "does not combine criterion and rater slope blocks")
  expect_false(grepl("one common slope for the whole model", contract, fixed = TRUE))
})

test_that("weighting comparison cannot enable GPCM ranking from supplied comparison flags", {
  make_fit <- function(model, method) {
    list(config = list(model = model, method = method))
  }
  make_comparison <- function(method,
                              ready,
                              ic_comparable,
                              ic_selectable = ic_comparable) {
    list(
      table = data.frame(
        Label = c(paste0("PCM/", method), paste0("GPCM/", method)),
        LogLik = c(-100, -97),
        stringsAsFactors = FALSE
      ),
      preferred = list(AIC = "GPCM/MML", BIC = "PCM/MML", SABIC = "GPCM/MML"),
      comparison_basis = list(
        same_data = TRUE,
        all_inference_ready = ready,
        ic_comparable = ic_comparable,
        all_ic_selectable = ic_selectable
      )
    )
  }

  mml <- mfrmr:::.weighting_review_comparison_contract(
    rasch_fit = make_fit("PCM", "MML"),
    gpcm_fit = make_fit("GPCM", "MML"),
    comparison = make_comparison("MML", ready = TRUE, ic_comparable = TRUE),
    reference_model = "PCM",
    aligned_pcm_owner = TRUE,
    pcm_gpcm_lrt = "withheld_current_scope"
  )
  expect_identical(mml$EvidenceTier, "mml_numerical_review_only")
  expect_false(mml$FormalModelSelectionAvailable)
  expect_equal(mml$ObservedLogLikDifference, 3)
  expect_true(is.na(mml$AICPreferred))
  expect_identical(
    mml$LogLikDifferenceStatus,
    "numerical_review_required"
  )

  jml <- mfrmr:::.weighting_review_comparison_contract(
    rasch_fit = make_fit("PCM", "JML"),
    gpcm_fit = make_fit("GPCM", "JML"),
    comparison = make_comparison("JML", ready = TRUE, ic_comparable = FALSE),
    reference_model = "PCM",
    aligned_pcm_owner = TRUE,
    pcm_gpcm_lrt = "withheld_current_scope"
  )
  expect_identical(jml$EvidenceTier, "jml_numerical_review_only")
  expect_false(jml$FormalModelSelectionAvailable)
  expect_true(is.na(jml$AICPreferred))
  expect_identical(
    jml$LogLikDifferenceStatus,
    "numerical_review_required"
  )
  expect_identical(
    jml$FACETSComparisonRole,
    "PCM_JML_side_only_no_FACETS_free_slope_GPCM_counterpart"
  )
})


test_that("weighting review distinguishes numerical convergence from inference readiness", {
  state <- "ready"
  testthat::local_mocked_bindings(
    mfrmr_get_readiness_record = function(...) list(fit = data.frame(NumericalState = state)),
    .package = "mfrmr"
  )
  contract <- function() mfrmr:::.weighting_review_comparison_contract(
    list(config = list(model = "PCM", method = "MML")),
    list(config = list(model = "GPCM", method = "MML")),
    list(table = data.frame(LogLik = c(-100, -97)),
         comparison_basis = list(same_data = TRUE, all_inference_ready = FALSE)),
    "PCM", TRUE, "withheld_current_scope"
  )
  ready <- contract()
  expect_true(ready$BothNumericallyReady)
  expect_false(ready$BothInferenceReady)
  expect_false(ready$FormalModelSelectionAvailable)
  expect_identical(ready$EvidenceTier, "mml_descriptive_not_inference_ready")
  expect_identical(ready$LogLikDifferenceStatus, "descriptive_noncomparable")
  expect_match(ready$NumericalReview, "checks passed for both fits", fixed = TRUE)
  state <- "review"
  pending <- contract()
  expect_false(pending$BothNumericallyReady)
  expect_identical(pending$LogLikDifferenceStatus, "numerical_review_required")
})
