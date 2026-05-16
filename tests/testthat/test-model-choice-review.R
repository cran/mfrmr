test_that("build_model_choice_review bundles comparison and user guidance", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:8]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  rsm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 18
  ))
  gpcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 5,
    maxit = 18
  ))

  review <- build_model_choice_review(RSM = rsm_fit, GPCM = gpcm_fit)

  expect_s3_class(review, "mfrm_model_choice_review")
  expect_s3_class(review, "mfrm_bundle")
  expect_true(all(c(
    "overview", "comparison", "comparison_table", "model_roles",
    "downstream_routes", "report_templates", "route_map", "support_status",
    "weighting_review_status", "weighting_review", "notes", "settings"
  ) %in% names(review)))
  expect_true(isTRUE(review$overview$HasBoundedGPCM[1]))
  expect_identical(review$overview$OperationalReference[1], "RSM")
  expect_identical(review$overview$SensitivityModel[1], "GPCM")
  expect_true(any(review$model_roles$RecommendedRole == "equal_weighting_reference"))
  expect_true(any(review$model_roles$RecommendedRole == "slope_aware_sensitivity"))
  expect_true(any(grepl("slope_facet == step_facet", review$model_roles$ScoreContract, fixed = TRUE)))
  expect_true(any(review$downstream_routes$FullAPARoute == "blocked" &
                    review$downstream_routes$Model == "GPCM"))
  expect_true(any(review$downstream_routes$FairAverage == "supported_with_caveat" &
                    review$downstream_routes$Model == "GPCM"))
  expect_true(any(grepl("automatic operational-scoring decision", review$key_warnings, fixed = TRUE)))
  expect_true(nrow(review$support_status) > 0)
  expect_false(isTRUE(review$weighting_review_status$Requested[1]))
  expect_false("weighting_audit_status" %in% names(review))
  expect_false("weighting_audit" %in% names(review))
  expect_false("run_weighting_audit" %in% names(review$settings))

  sx <- summary(review)
  expect_s3_class(sx, "summary.mfrm_model_choice_review")
  expect_true(all(c(
    "overview", "comparison_table", "model_roles", "downstream_routes",
    "report_templates", "route_map", "weighting_review_status"
  ) %in% names(sx)))
  expect_false("weighting_audit_status" %in% names(sx))
})

test_that("build_model_choice_review can attach the detailed weighting review", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:8]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  rsm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 18
  ))
  gpcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 5,
    maxit = 18
  ))

  review <- build_model_choice_review(
    RSM = rsm_fit,
    GPCM = gpcm_fit,
    run_weighting_review = TRUE,
    theta_points = 11,
    top_n = 3
  )

  expect_true(isTRUE(review$weighting_review_status$Requested[1]))
  expect_true(isTRUE(review$weighting_review_status$Available[1]))
  expect_s3_class(review$weighting_review, "mfrm_weighting_review")
})

test_that("build_model_choice_review handles RSM versus PCM without GPCM routes", {
  rsm_fit <- make_toy_fit(method = "JML", model = "RSM", maxit = 12)
  toy <- load_mfrmr_data("example_core")
  pcm_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    model = "PCM",
    step_facet = "Criterion",
    maxit = 12
  ))

  review <- build_model_choice_review(RSM = rsm_fit, PCM = pcm_fit)

  expect_s3_class(review, "mfrm_model_choice_review")
  expect_false(isTRUE(review$overview$HasBoundedGPCM[1]))
  expect_equal(nrow(review$support_status), 0)
  expect_true(all(review$downstream_routes$FullAPARoute == "supported"))
  expect_true(any(grepl("score interpretation", review$next_actions, fixed = TRUE)))
})
