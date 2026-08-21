test_that("build_weighting_review returns a structured review bundle", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  rasch_fit <- fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 7,
    maxit = 25
  )
  gpcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 7,
    maxit = 25
  ))

  audit <- build_weighting_review(rasch_fit, gpcm_fit, theta_points = 31, top_n = 6)

  expect_s3_class(audit, "mfrm_weighting_review")
  expect_s3_class(audit, "mfrm_bundle")
  expect_true(all(c(
    "overview", "status", "model_comparison", "facet_shift", "slope_profile",
    "comparison_contract",
    "information_redistribution", "top_reweighted_levels", "plot_map",
    "reporting_map", "support_status", "notes", "settings"
  ) %in% names(audit)))
  expect_true(all(c(
    "SlopeFacet", "StepFacet", "ReferenceStepFacet",
    "AlignedStepSlopeOwner", "SlopeLevelCount",
    "FreeRelativeSlopeContrasts", "SlopeComposition",
    "SimultaneousCriterionRaterSlopeBlocks", "UnitSlopePCMReduction",
    "PCMvsGPCMLRT", "EvidenceTier", "FormalModelSelectionAvailable",
    "ObservedLogLikDifference", "LogLikDifferenceStatus"
  ) %in% names(audit$overview)))
  expect_true(all(c(
    "EvidenceTier", "FormalModelSelectionAvailable", "SelectionRoute",
    "ObservedLogLikDifference", "LogLikDifferenceStatus",
    "FACETSComparisonRole", "RecommendedUse"
  ) %in% names(audit$comparison_contract)))
  expect_identical(
    audit$overview$EvidenceTier,
    audit$comparison_contract$EvidenceTier
  )
  expect_false(audit$overview$FormalModelSelectionAvailable)
  expect_identical(audit$overview$SlopeFacet, "Criterion")
  expect_identical(audit$overview$StepFacet, "Criterion")
  expect_true(is.na(audit$overview$ReferenceStepFacet))
  expect_false(audit$overview$AlignedStepSlopeOwner)
  expect_identical(
    audit$overview$SlopeLevelCount,
    as.integer(nrow(audit$slope_profile))
  )
  expect_identical(
    audit$overview$FreeRelativeSlopeContrasts,
    as.integer(nrow(audit$slope_profile) - 1L)
  )
  expect_false(audit$overview$SimultaneousCriterionRaterSlopeBlocks)
  expect_false(audit$overview$UnitSlopePCMReduction)
  expect_identical(
    audit$overview$PCMvsGPCMLRT,
    "not_applicable_reference_is_not_pcm"
  )
  expect_true(all(c(
    "Facet", "Level", "ReferenceEstimate", "ComparisonEstimate",
    "DeltaEstimate", "AbsDeltaEstimate", "ReferenceRank", "ComparisonRank",
    "RankShift"
  ) %in% names(audit$facet_shift)))
  expect_true(all(c(
    "SlopeFacet", "Estimate", "LogEstimate", "RelativeWeight",
    "WeightingDirection"
  ) %in% names(audit$slope_profile)))
  expect_true(all(c(
    "Facet", "Level", "ReferenceInfoShare", "ComparisonInfoShare",
    "InfoShareDelta"
  ) %in% names(audit$information_redistribution)))
})

test_that("summary methods for build_weighting_review expose front-door tables", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  rasch_fit <- fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "PCM",
    step_facet = "Criterion",
    quad_points = 7,
    maxit = 25
  )
  gpcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 7,
    maxit = 25
  ))

  audit <- build_weighting_review(rasch_fit, gpcm_fit, theta_points = 21, top_n = 5)
  sx <- summary(audit, top_n = 3)

  expect_true(audit$overview$UnitSlopePCMReduction)
  expect_identical(audit$overview$ReferenceStepFacet, "Criterion")
  expect_true(audit$overview$AlignedStepSlopeOwner)
  expect_identical(audit$overview$PCMvsGPCMLRT, "withheld_current_scope")
  expect_true(any(grepl(
    "PCM-versus-GPCM chi-square LRT",
    audit$key_warnings,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "other facets have no separate slope block",
    audit$notes,
    fixed = TRUE
  )))

  expect_s3_class(sx, "summary.mfrm_weighting_review")
  expect_true(all(c(
    "overview", "status", "comparison_contract", "key_warnings", "next_actions",
    "top_measure_shifts", "top_reweighted_levels",
    "plot_map", "reporting_map", "support_status"
  ) %in% names(sx)))
  expect_lte(nrow(sx$top_measure_shifts), 3)
  expect_lte(nrow(sx$top_reweighted_levels), 3)
})

test_that("JML weighting review is typed as descriptive or optimizer-trace evidence", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  pcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    model = "PCM",
    step_facet = "Criterion",
    maxit = 15
  ))
  gpcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    maxit = 15
  ))

  audit <- build_weighting_review(
    pcm_fit,
    gpcm_fit,
    theta_points = 21,
    top_n = 5
  )

  expect_match(audit$comparison_contract$EvidenceTier, "^jml_")
  expect_false(audit$comparison_contract$FormalModelSelectionAvailable)
  expect_identical(
    audit$comparison_contract$SelectionRoute,
    "withheld_JML_has_no_automatic_PCM_GPCM_selection"
  )
  expect_true(audit$comparison_contract$LogLikDifferenceStatus %in% c(
    "optimizer_trace_only_not_inference_ready",
    "descriptive_unpenalized_gain_not_selection",
    "unavailable"
  ))
  expect_identical(
    audit$comparison_contract$FACETSComparisonRole,
    "PCM_JML_side_only_no_FACETS_free_slope_GPCM_counterpart"
  )
  expect_true(any(grepl(
    "JML log-likelihood difference is unpenalized",
    audit$key_warnings,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "unit_slopes",
    audit$next_actions,
    fixed = TRUE
  )))
})

test_that("build_weighting_review requires shared prepared response data", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  rasch_fit <- fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 7,
    maxit = 25
  )

  toy_less <- toy[toy$Person != keep_people[[1]], , drop = FALSE]
  gpcm_fit <- suppressWarnings(fit_mfrm(
    toy_less,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 7,
    maxit = 25
  ))

  expect_error(
    build_weighting_review(rasch_fit, gpcm_fit, theta_points = 21),
    "same prepared response data"
  )
})
