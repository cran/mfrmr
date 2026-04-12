test_that("build_weighting_audit returns a structured review bundle", {
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

  audit <- build_weighting_audit(rasch_fit, gpcm_fit, theta_points = 31, top_n = 6)

  expect_s3_class(audit, "mfrm_weighting_audit")
  expect_s3_class(audit, "mfrm_bundle")
  expect_true(all(c(
    "overview", "status", "model_comparison", "facet_shift", "slope_profile",
    "information_redistribution", "top_reweighted_levels", "plot_map",
    "reporting_map", "support_status", "notes", "settings"
  ) %in% names(audit)))
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

test_that("summary methods for build_weighting_audit expose front-door tables", {
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

  audit <- build_weighting_audit(rasch_fit, gpcm_fit, theta_points = 21, top_n = 5)
  sx <- summary(audit, top_n = 3)

  expect_s3_class(sx, "summary.mfrm_weighting_audit")
  expect_true(all(c(
    "overview", "status", "key_warnings", "next_actions",
    "top_measure_shifts", "top_reweighted_levels",
    "plot_map", "reporting_map", "support_status"
  ) %in% names(sx)))
  expect_lte(nrow(sx$top_measure_shifts), 3)
  expect_lte(nrow(sx$top_reweighted_levels), 3)
})

test_that("build_weighting_audit requires shared prepared response data", {
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
    build_weighting_audit(rasch_fit, gpcm_fit, theta_points = 21),
    "same prepared response data"
  )
})
