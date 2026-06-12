# --------------------------------------------------------------------------
# test-output-stability.R
# Regression guard: output structures must remain stable across versions.
# --------------------------------------------------------------------------

# === 5.1 fit_mfrm output structure ========================================

test_that("fit_mfrm returns all required components", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  expect_s3_class(fit, "mfrm_fit")
  # Required top-level components
  expect_true(all(c("summary", "facets", "steps", "population", "config", "prep", "opt") %in%
    names(fit)))
  # Person table structure
  expect_true(is.data.frame(fit$facets$person))
  expect_true(all(c("Person", "Estimate") %in% names(fit$facets$person)))
  # Others table structure
  expect_true(is.data.frame(fit$facets$others))
  expect_true(all(c("Facet", "Level", "Estimate") %in%
    names(fit$facets$others)))
  # Summary structure
  expect_true(is.data.frame(fit$summary))
  expect_true("LogLik" %in% names(fit$summary))
  # Steps structure
  expect_true(is.data.frame(fit$steps))
  expect_true(all(c("Step", "Estimate") %in% names(fit$steps)))
  # Population scaffold structure
  expect_true(is.list(fit$population))
  expect_false(isTRUE(fit$population$active))
  expect_identical(fit$population$posterior_basis, "legacy_mml")
  expect_false(isTRUE(fit$config$population_active))
  expect_identical(fit$config$posterior_basis, "legacy_mml")
  expect_true(is.list(fit$config$population_spec))
  expect_false(isTRUE(fit$config$population_spec$active))
  expect_identical(fit$config$population_spec$posterior_basis, "legacy_mml")
  # Column types
  expect_type(fit$facets$others$Estimate, "double")
  expect_type(fit$facets$others$Level, "character")
  expect_type(fit$facets$others$Facet, "character")
})

test_that("summary.mfrm_fit reports legacy population basis for ordinary fits", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", maxit = 30, quad_points = 7
  ))
  s <- summary(fit)

  expect_true("population_overview" %in% names(s))
  expect_true("population_design" %in% names(s))
  expect_true("population_coding" %in% names(s))
  expect_true("status" %in% names(s))
  expect_true("key_warnings" %in% names(s))
  expect_true("next_actions" %in% names(s))
  expect_true("settings_overview" %in% names(s))
  expect_true("reporting_map" %in% names(s))
  expect_true(is.data.frame(s$population_overview))
  expect_true(is.data.frame(s$population_design))
  expect_true(is.data.frame(s$population_coding))
  expect_true(is.data.frame(s$status))
  expect_true(is.character(s$key_warnings))
  expect_true(is.character(s$next_actions))
  expect_true(is.data.frame(s$settings_overview))
  expect_true(is.data.frame(s$reporting_map))
  expect_false(isTRUE(s$population_overview$PopulationModel[1]))
  expect_equal(nrow(s$population_design), 0L)
  expect_equal(nrow(s$population_coding), 0L)
  expect_identical(as.character(s$population_overview$PosteriorBasis[1]), "legacy_mml")
  expect_true(all(c("StepFacet", "SlopeFacet", "NoncenterFacet", "QuadPoints") %in% names(s$settings_overview)))
  expect_true(all(c("Area", "CoveredHere", "CompanionOutput") %in% names(s$reporting_map)))
  expect_true(any(grepl("legacy unconditional prior", s$notes, fixed = TRUE)))
})

test_that("GPCM summaries expose slope overview and diagnostics are now available", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    model = "GPCM", step_facet = "Criterion",
    method = "MML", quad_points = 5, maxit = 20
  ))
  s <- summary(fit)

  expect_true("slopes" %in% names(fit))
  expect_true(is.data.frame(fit$slopes))
  expect_true("slope_overview" %in% names(s))
  expect_true(is.data.frame(s$slope_overview))
  expect_true(isTRUE(s$slope_overview$Positive[1]))
  expect_equal(s$slope_overview$GeometricMean[1], 1, tolerance = 1e-6)
  dx <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "overall")
  pca <- analyze_residual_pca(dx, mode = "overall")
  unexp <- unexpected_response_table(fit, diagnostics = dx, top_n = 10)
  disp <- displacement_table(fit, diagnostics = dx, top_n = 10)
  meas <- measurable_summary_table(fit, diagnostics = dx)
  rate <- rating_scale_table(fit, diagnostics = dx)
  ir <- interrater_agreement_table(fit, diagnostics = dx, rater_facet = "Rater", top_n = 5)
  chk <- reporting_checklist(fit, diagnostics = dx)
  dash <- facet_quality_dashboard(fit, diagnostics = dx, facet = "Rater")
  qc <- plot_qc_dashboard(fit, diagnostics = dx, draw = FALSE)

  expect_s3_class(dx, "mfrm_diagnostics")
  expect_identical(dx$diagnostic_mode, "both")
  expect_true(isTRUE(dx$marginal_fit$available))
  expect_identical(as.character(dx$marginal_fit$summary$Model[1]), "GPCM")
  expect_s3_class(pca, "mfrm_residual_pca")
  expect_s3_class(unexp, "mfrm_unexpected")
  expect_s3_class(disp, "mfrm_displacement")
  expect_s3_class(meas, "mfrm_measurable")
  expect_s3_class(rate, "mfrm_rating_scale")
  expect_s3_class(ir, "mfrm_interrater")
  expect_s3_class(chk, "mfrm_reporting_checklist")
  expect_s3_class(dash, "mfrm_facet_dashboard")
  expect_s3_class(qc, "mfrm_plot_data")
  expect_true(isFALSE(dx$fair_average$available))
  expect_true(grepl("disabled for bounded `GPCM`", dx$fair_average$reason, fixed = TRUE))
  expect_true(grepl("slope-aware", dx$fair_average$reason, fixed = TRUE))
  expect_true(isFALSE(qc$data$fair_average$available))
  expect_true(grepl("slope-aware", qc$data$fair_average$reason, fixed = TRUE))

  qc_pipe <- run_qc_pipeline(fit)
  expect_s3_class(qc_pipe, "mfrm_qc_pipeline")
  expect_equal(nrow(qc_pipe$verdicts), 10)
  expect_true(any(qc_pipe$gpcm_boundary$Area == "APA writer and fit-based export bundles"))
  # `fair_average_table()` and `estimate_bias()` are both unblocked for
  # GPCM fits in 0.2.0 under the slope-aware element-conditional GPCM
  # construction. Confirm they return populated bundles with the GPCM
  # caveat instead of raising the legacy error.
  fa <- fair_average_table(fit)
  expect_s3_class(fa, "mfrm_fair_average")
  expect_true(nrow(fa$stacked) > 0)
  expect_identical(fa$settings$method, "GPCM-slope-aware")
  expect_true(!is.null(fa$caveat))
  expect_true(grepl("slope-aware element-conditional", fa$caveat, fixed = TRUE))
  expect_true(grepl("Standard errors", fa$caveat, fixed = TRUE))

  bias <- estimate_bias(fit, dx, facet_a = "Rater", facet_b = "Criterion")
  expect_s3_class(bias, "mfrm_bias")
  expect_true(nrow(bias$table) > 0)
  expect_identical(bias$method, "GPCM-slope-aware")
  expect_true(!is.null(bias$caveat))
  expect_true(grepl("slope-aware GPCM kernel", bias$caveat, fixed = TRUE))

  expect_error(
    facets_output_contract_review(fit),
    "does not support `GPCM` fits",
    fixed = TRUE
  )
})

test_that("curve/report GPCM workflows open where the probability kernel is already generalized", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    model = "GPCM", step_facet = "Criterion",
    method = "MML", quad_points = 5, maxit = 20
  ))

  p_bundle <- plot(fit, type = "bundle", draw = FALSE)
  p_ccc <- plot(fit, type = "ccc", draw = FALSE)
  p_path <- plot(fit, type = "pathway", draw = FALSE)
  cs <- category_structure_report(fit)
  cc <- category_curves_report(fit)
  out_graph <- facets_output_file_bundle(fit, include = "graph", write_files = FALSE)

  expect_s3_class(p_bundle, "mfrm_plot_bundle")
  expect_s3_class(p_ccc, "mfrm_plot_data")
  expect_s3_class(p_path, "mfrm_plot_data")
  expect_s3_class(cs, "mfrm_category_structure")
  expect_s3_class(cc, "mfrm_category_curves")
  expect_true(all(c("graphfile", "graphfile_syntactic", "settings") %in% names(out_graph)))

  expect_error(
    estimation_iteration_report(fit, max_iter = 3),
    "does not support `GPCM` fits",
    fixed = TRUE
  )
  out_score <- facets_output_file_bundle(fit, include = "score", write_files = FALSE)
  expect_s3_class(out_score, "mfrm_output_bundle")
  expect_true(all(c("scorefile", "gpcm_score_side_contract", "gpcm_boundary") %in% names(out_score)))
  expect_true(all(c(
    "ScoreSlope", "ObsProb", "ScoreSideStatus", "ScoreUncertaintyStatus",
    "ScoreUncertaintyMethod", "ExpectedScoreSE", "ExpectedScoreCI_Lower",
    "ExpectedScoreCI_Upper", "ResidualSE", "ScoreSideSE",
    "ScoreSideCI_Lower", "ScoreSideCI_Upper", "ScoreSideSE_Status",
    "ScoreSideSE_Method", "ScoreSideCaveat"
  ) %in% names(out_score$scorefile)))
  expect_identical(out_score$settings$score_se_method, "both")
  expect_true(any(out_score$scorefile$ScoreSideStatus == "supported_with_caveat"))
  expect_false(any(out_score$scorefile$ScoreUncertaintyStatus == "score_side_se_not_exported"))
  expect_true(any(is.finite(out_score$scorefile$ExpectedScoreSE)))
  expect_true(any(is.finite(out_score$scorefile$ScoreSideSE)))

  score_se_plot <- plot(out_score, type = "score_se", draw = FALSE)
  expect_s3_class(score_se_plot, "mfrm_plot_data")
  expect_identical(score_se_plot$data$se_column, "ScoreSideSE")

  out_summary <- summary(out_score)
  expect_true("ScoreSEMethod" %in% names(out_summary$summary))
  expect_identical(as.character(out_summary$summary$ScoreSEMethod[1]), "both")
})

# === 5.2 MML output has SD column =========================================

test_that("MML fit includes person SD column", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "MML", maxit = 30))
  expect_true("SD" %in% names(fit$facets$person))
  expect_true(all(fit$facets$person$SD > 0))
})

test_that("MML diagnostics expose model and real precision columns", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "MML", maxit = 30, quad_points = 7))
  dx <- diagnose_mfrm(fit, residual_pca = "none")
  se_mask <- is.finite(dx$measures$SE) & is.finite(dx$measures$ModelSE)
  real_mask <- is.finite(dx$measures$RealSE) & is.finite(dx$measures$ModelSE)
  rel_mask <- is.finite(dx$reliability$RealReliability) & is.finite(dx$reliability$Reliability)

  expect_true(all(c(
    "SE", "ModelSE", "RealSE", "SE_Method", "PrecisionTier",
    "Converged", "SupportsFormalInference", "SEUse", "CIBasis", "CIUse",
    "CI_Level", "CI_Method", "CIEligible", "CILabel"
  ) %in% names(dx$measures)))
  expect_true(all(dx$measures$SE[se_mask] == dx$measures$ModelSE[se_mask]))
  expect_true(all(dx$measures$RealSE[real_mask] >= dx$measures$ModelSE[real_mask]))
  expect_true(all(c(
    "ModelReliability", "RealReliability", "Converged", "PrecisionTier",
    "SupportsFormalInference", "ReliabilityUse"
  ) %in% names(dx$reliability)))
  expect_true(all(dx$reliability$RealReliability[rel_mask] <= dx$reliability$Reliability[rel_mask]))
  expect_true(all(dx$measures$CIEligible == dx$measures$SupportsFormalInference))
  expect_equal(unique(dx$measures$CI_Level), 0.95)
  expect_true(all(dx$measures$CI_Method == "Normal approximation"))
  expect_true(all(dx$measures$Converged == fit$summary$Converged[1]))
})

# === 5.3 diagnose_mfrm output structure ===================================

test_that("diagnose_mfrm returns all required components", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  dx <- diagnose_mfrm(fit, residual_pca = "none")
  expect_s3_class(dx, "mfrm_diagnostics")
  # Required components
  required <- c("obs", "measures", "overall_fit", "reliability", "precision_profile", "precision_review", "facet_precision")
  expect_true(all(required %in% names(dx)))
  expect_false("precision_audit" %in% names(dx))
  # measures columns
  expect_true(all(c(
    "Facet", "Level", "Estimate", "SE", "ModelSE", "RealSE",
    "Converged", "PrecisionTier", "SupportsFormalInference", "SEUse",
    "CIBasis", "CIUse", "CI_Level", "CI_Method", "CIEligible", "CILabel",
    "Infit", "Outfit"
  ) %in%
    names(dx$measures)))
  # reliability columns
  expect_true(all(c(
    "Facet", "Separation", "Reliability", "Converged", "PrecisionTier",
    "SupportsFormalInference", "ReliabilityUse"
  ) %in%
    names(dx$reliability)))
  expect_true(all(c("Method", "Converged", "PrecisionTier", "SupportsFormalInference", "HasFallbackSE", "RecommendedUse") %in%
    names(dx$precision_profile)))
  expect_true(all(c("Check", "Status", "Detail") %in%
    names(dx$precision_review)))
  expect_true(all(c("Facet", "DistributionBasis", "SEMode", "SEColumn", "Separation", "Reliability") %in%
    names(dx$facet_precision)))
  # obs should be data.frame
  expect_true(is.data.frame(dx$obs))
})

test_that("diagnostics and data-description summaries expose paper-reporting maps", {
  d_local <- mfrmr:::sample_mfrm_data(seed = 42)
  fit_local <- suppressWarnings(fit_mfrm(
    d_local, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30
  ))
  dx_local <- diagnose_mfrm(fit_local, residual_pca = "none")

  dx_sum <- summary(dx_local)
  expect_s3_class(dx_sum, "summary.mfrm_diagnostics")
  expect_true("reporting_map" %in% names(dx_sum))
  expect_true(is.data.frame(dx_sum$reporting_map))
  expect_true(all(c("Area", "CoveredHere", "CompanionOutput") %in% names(dx_sum$reporting_map)))

  ds <- describe_mfrm_data(d_local, "Person", c("Rater", "Task", "Criterion"), "Score")
  ds_sum <- summary(ds)
  expect_s3_class(ds_sum, "summary.mfrm_data_description")
  expect_true("reporting_map" %in% names(ds_sum))
  expect_true(is.data.frame(ds_sum$reporting_map))
  expect_true(all(c("Area", "CoveredHere", "CompanionOutput") %in% names(ds_sum$reporting_map)))
})

# === 5.4 Table builder output stability ====================================
#
# Shared fixture: fit + diagnostics computed once for table builder tests.

.stab_d   <- mfrmr:::sample_mfrm_data(seed = 42)
.stab_fit <- suppressWarnings(fit_mfrm(.stab_d, "Person",
  c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
.stab_dx  <- diagnose_mfrm(.stab_fit, residual_pca = "none")

test_that("unexpected_response_table returns table + summary", {
  res <- unexpected_response_table(.stab_fit, diagnostics = .stab_dx)
  expect_true(all(c("table", "summary") %in% names(res)))
  expect_true(is.data.frame(res$table))
  expect_true(is.data.frame(res$summary))
})

test_that("fair_average_table returns stacked + by_facet", {
  res <- fair_average_table(.stab_fit, diagnostics = .stab_dx)
  alias_mask <- is.finite(res$stacked$AdjustedAverage) & is.finite(res$stacked$`Fair(M) Average`)
  se_mask <- is.finite(res$stacked$ModelBasedSE) & is.finite(res$stacked$`Model S.E.`)
  fit_adj_mask <- is.finite(res$stacked$FitAdjustedSE) & is.finite(res$stacked$`Real S.E.`)
  expect_true(all(c("stacked", "by_facet") %in% names(res)))
  expect_true(is.data.frame(res$stacked))
  expect_true(is.list(res$by_facet))
  expect_true(all(c(
    "ObservedAverage", "AdjustedAverage", "StandardizedAdjustedAverage",
    "ModelBasedSE", "FitAdjustedSE"
  ) %in% names(res$stacked)))
  expect_true(all(res$stacked$AdjustedAverage[alias_mask] == res$stacked$`Fair(M) Average`[alias_mask]))
  expect_true(all(res$stacked$ModelBasedSE[se_mask] == res$stacked$`Model S.E.`[se_mask]))
  expect_true(all(res$stacked$FitAdjustedSE[fit_adj_mask] == res$stacked$`Real S.E.`[fit_adj_mask]))
})

test_that("fair_average_table native label style hides legacy aliases", {
  res <- fair_average_table(.stab_fit, diagnostics = .stab_dx, reference = "mean", label_style = "native")
  expect_true(all(c("ObservedAverage", "AdjustedAverage", "ModelBasedSE", "FitAdjustedSE") %in% names(res$stacked)))
  expect_false(any(c("Obsvd Average", "Fair(M) Average", "Fair(Z) Average", "Model S.E.", "Real S.E.") %in%
    names(res$stacked)))
  expect_false("StandardizedAdjustedAverage" %in% names(res$stacked))
})

test_that("displacement_table returns table + summary", {
  res <- displacement_table(.stab_fit, diagnostics = .stab_dx,
    anchored_only = FALSE)
  expect_true(all(c("table", "summary") %in% names(res)))
  expect_true(is.data.frame(res$table))
  expect_true(is.data.frame(res$summary))
})

test_that("rating_scale_table returns category_table + threshold_table", {
  res <- rating_scale_table(.stab_fit, diagnostics = .stab_dx)
  expect_true(all(c("category_table", "threshold_table") %in% names(res)))
  expect_true(is.data.frame(res$category_table))
  expect_true(is.data.frame(res$threshold_table))
})

test_that("measurable_summary_table returns summary + facet_coverage", {
  res <- measurable_summary_table(.stab_fit, diagnostics = .stab_dx)
  expect_true(all(c("summary", "facet_coverage") %in% names(res)))
  expect_true(is.data.frame(res$summary))
  expect_true(is.data.frame(res$facet_coverage))
})

test_that("interrater_agreement_table returns pairs + summary", {
  res <- interrater_agreement_table(.stab_fit, diagnostics = .stab_dx,
    rater_facet = "Rater")
  expect_true(all(c("pairs", "summary") %in% names(res)))
  expect_true(is.data.frame(res$pairs))
  expect_true(is.data.frame(res$summary))
  expect_true(all(c("OpportunityCount", "ExactCount", "ExpectedExactCount", "AdjacentCount") %in%
    names(res$pairs)))
  expect_true(all(c("AgreementMinusExpected", "RaterSeparation", "RaterReliability") %in%
    names(res$summary)))
})

test_that("facets_chisq_table returns table + summary", {
  res <- facets_chisq_table(.stab_fit, diagnostics = .stab_dx)
  expect_true(all(c("table", "summary") %in% names(res)))
  expect_true(is.data.frame(res$table))
  expect_true(is.data.frame(res$summary))
})

# === 5.5 Extreme data - all same score =====================================

test_that("all-same-score data produces informative error", {
  d <- data.frame(
    Person = rep(paste0("P", 1:6), each = 3),
    Rater  = rep(paste0("R", 1:3), 6),
    Score  = rep(3, 18),
    stringsAsFactors = FALSE
  )
  expect_error(
    fit_mfrm(d, "Person", "Rater", "Score", method = "JML"),
    "one score"
  )
})

# === 5.6 All max score =====================================================

test_that("all-max-score data produces informative error", {
  d <- data.frame(
    Person = rep(paste0("P", 1:6), each = 3),
    Rater  = rep(paste0("R", 1:3), 6),
    Score  = rep(5, 18),
    stringsAsFactors = FALSE
  )
  expect_error(
    fit_mfrm(d, "Person", "Rater", "Score", method = "JML"),
    "one score"
  )
})

# === 5.7 Output length matches facet levels ================================

test_that("output rows match number of facet levels", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  expect_equal(nrow(fit$facets$person), length(unique(d$Person)))
  n_rater <- nrow(fit$facets$others |>
    dplyr::filter(Facet == "Rater"))
  expect_equal(n_rater, length(unique(d$Rater)))
  n_task <- nrow(fit$facets$others |>
    dplyr::filter(Facet == "Task"))
  expect_equal(n_task, length(unique(d$Task)))
})

# === 5.8 Numeric types for estimates and SEs ===============================

test_that("all Estimate and SE columns are numeric", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  dx <- diagnose_mfrm(fit, residual_pca = "none")
  # Estimate columns are numeric/double
  expect_type(fit$facets$person$Estimate, "double")
  expect_type(fit$facets$others$Estimate, "double")
  expect_type(dx$measures$Estimate, "double")
  expect_type(fit$steps$Estimate, "double")
  # SE columns are numeric and positive where present
  se_vals <- dx$measures$SE
  se_finite <- se_vals[is.finite(se_vals)]
  expect_true(all(se_finite > 0))
  expect_true(all(dx$measures$ModelSE[is.finite(dx$measures$ModelSE)] > 0))
  expect_true(all(dx$measures$RealSE[is.finite(dx$measures$RealSE)] > 0))
})

# === 5.9 Fit statistics are in expected range ==============================

test_that("fit statistics are in expected range", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  dx <- diagnose_mfrm(fit, residual_pca = "none")
  # Infit and Outfit are positive (> 0) where finite
  infit_vals <- dx$measures$Infit[is.finite(dx$measures$Infit)]
  outfit_vals <- dx$measures$Outfit[is.finite(dx$measures$Outfit)]
  expect_true(all(infit_vals > 0))
  expect_true(all(outfit_vals > 0))
  # PTMEA is between -1 and 1 where present
  if ("PTMEA" %in% names(dx$measures)) {
    ptmea_vals <- dx$measures$PTMEA[is.finite(dx$measures$PTMEA)]
    expect_true(all(ptmea_vals >= -1 & ptmea_vals <= 1))
  }
})
