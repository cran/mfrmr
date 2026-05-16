test_that("facets_fit_review compares engine and FACETS-style fit standardization", {
  d <- mfrmr:::sample_mfrm_data(seed = 321)
  fit <- suppressWarnings(
    mfrmr::fit_mfrm(
      data = d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "JML",
      model = "RSM",
      maxit = 20
    )
  )
  diag_engine <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")

  review <- mfrmr::facets_fit_review(fit, diagnostics = diag_engine)

  expect_s3_class(review, "mfrm_facets_fit_review")
  expect_s3_class(review, "mfrm_bundle")
  expect_true(all(c(
    "summary", "standardization", "df_sensitivity", "df_sensitive",
    "df_sensitivity_summary", "external_table_quality", "external_comparison",
    "df_conversion_guide", "guidance", "settings"
  ) %in% names(review)))
  expect_false("internal_comparison" %in% names(review))
  expect_s3_class(review$df_conversion_guide, "mfrm_facets_fit_df_guide")
  expect_true(any(grepl("MnSq first", review$df_conversion_guide$summary$PrimaryRule, fixed = TRUE)))
  expect_gt(nrow(review$df_sensitivity), 0)
  expect_equal(nrow(review$external_comparison), 0)
  expect_identical(review$settings$intended_use, "fit_standardization_review")
  expect_true(all(c(
    "DF_Infit_ENGINE", "DF_Infit_FACETS",
    "InfitZSTD_ENGINE", "InfitZSTD_FACETS",
    "MaxDFRelativeDifference_ENGINE_vs_FACETS",
    "DfSensitivityStatus", "Interpretation", "FlagChangedByDf"
  ) %in% names(review$df_sensitivity)))
  expect_true(all(c(
    "ComparedRows", "FlagChangedByDfRows", "DfConventionDifferenceRows"
  ) %in% names(review$df_sensitivity_summary)))
  expect_true(all(c(
    "DfComparedRows", "DfSensitiveRows", "DfSameOrRoundingRows"
  ) %in% names(review$summary)))
  expect_equal(review$summary$DfSensitiveRows, nrow(review$df_sensitive))
  expect_equal(review$settings$external_zstd_tolerance, 0.05)
  expect_equal(review$settings$df_zstd_tolerance, 0.05)
  expect_equal(review$settings$df_zstd_large_shift, 0.5)
  expect_error(
    mfrmr::facets_fit_review(fit, diagnostics = diag_engine, zstd_tolerance = 0.05),
    "unused argument"
  )

  s <- summary(review, top_n = 3)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_identical(s$summary_kind, "mfrm_facets_fit_review")
  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "mfrm_facets_fit_review", fixed = TRUE)
  p <- plot(review, type = "df_sensitivity", top_n = 4, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$data$plot, "df_sensitivity")
  expect_lte(nrow(p$data$table), 4L)
})

test_that("facets_fit_df_guide documents comparison routes", {
  guide <- mfrmr::facets_fit_df_guide()
  expect_s3_class(guide, "mfrm_facets_fit_df_guide")
  expect_s3_class(guide, "mfrm_bundle")
  expect_true(all(c(
    "summary", "formula_guide", "column_guide",
    "decision_guide", "interpretation_guide", "references"
  ) %in% names(guide)))
  expect_true(any(grepl("FACETS-style df", guide$formula_guide$Quantity, fixed = TRUE)))
  expect_true(any(grepl("diagnose_mfrm", guide$column_guide$Route, fixed = TRUE)))
  expect_true(any(grepl("MnSq", guide$decision_guide$Question, fixed = TRUE)))

  s <- summary(guide)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_identical(s$summary_kind, "mfrm_facets_fit_df_guide")
})

test_that("facets_fit_review matches an external FACETS-style fit table", {
  d <- mfrmr:::sample_mfrm_data(seed = 322)
  fit <- suppressWarnings(
    mfrmr::fit_mfrm(
      data = d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "JML",
      model = "RSM",
      maxit = 20
    )
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none", fit_df_method = "both")
  facets_like <- diag$fit |>
    dplyr::transmute(
      Facet = .data$Facet,
      Level = .data$Level,
      Infit = .data$Infit,
      Outfit = .data$Outfit,
      InfitZSTD = .data$InfitZSTD_FACETS,
      OutfitZSTD = .data$OutfitZSTD_FACETS,
      DF_Infit = .data$DF_Infit_FACETS,
      DF_Outfit = .data$DF_Outfit_FACETS,
      TCount = .data$N
    )

  review <- mfrmr::facets_fit_review(
    fit,
    diagnostics = diag,
    facets_fit = facets_like
  )

  expect_gt(nrow(review$external_comparison), 0)
  expect_true(all(review$external_comparison$ExternalMatched))
  expect_true(all(review$external_comparison$ExternalStatus == "same"))
  expect_equal(review$summary$ExternalRows, nrow(facets_like))
  expect_equal(review$external_table_quality$CompleteMnSqRows, nrow(facets_like))
  expect_equal(review$external_table_quality$CompleteZSTDRows, nrow(facets_like))
  expect_equal(review$external_table_quality$CompleteDFRows, nrow(facets_like))
  expect_equal(review$summary$ExternalDuplicateKeyRows, 0)
  expect_equal(review$summary$ExternalMatched, nrow(facets_like))
  expect_equal(review$summary$ExternalNeedsReview, 0)
})
