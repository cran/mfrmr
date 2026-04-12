# test-report-functions.R
# Tests for report-building and table functions in isolation.
# Uses a shared fixture fit + diagnostics for efficiency.

# ---- Shared fixture ----

local({
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  .fit <<- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  .diag <<- diagnose_mfrm(.fit, residual_pca = "both", pca_max_factors = 3)
  .bias <<- estimate_bias(.fit, .diag, facet_a = "Rater", facet_b = "Task")
})

# ---- specifications_report ----

test_that("specifications_report returns a bundle", {
  spec <- specifications_report(.fit)
  expect_s3_class(spec, "mfrm_bundle")
  s <- summary(spec)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- estimation_iteration_report ----

test_that("estimation_iteration_report returns a bundle", {
  expect_no_warning(
    iter <- estimation_iteration_report(.fit)
  )
  expect_s3_class(iter, "mfrm_bundle")
})

# ---- data_quality_report ----

test_that("data_quality_report returns a bundle", {
  dq <- data_quality_report(.fit)
  expect_s3_class(dq, "mfrm_bundle")
  s <- summary(dq)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- category_curves_report ----

test_that("category_curves_report returns a bundle", {
  cc <- category_curves_report(.fit)
  expect_s3_class(cc, "mfrm_bundle")
})

# ---- category_structure_report ----

test_that("category_structure_report returns a bundle", {
  cs <- category_structure_report(.fit, diagnostics = .diag)
  expect_s3_class(cs, "mfrm_bundle")
})

# ---- subset_connectivity_report ----

test_that("subset_connectivity_report returns a bundle", {
  sc <- subset_connectivity_report(.fit)
  expect_s3_class(sc, "mfrm_bundle")
  p_cov <- plot(sc, type = "linking_matrix", draw = FALSE)
  expect_s3_class(p_cov, "mfrm_plot_data")
  expect_identical(p_cov$data$plot, "coverage_matrix")
  expect_identical(p_cov$data$requested_type, "linking_matrix")
  expect_true(is.matrix(p_cov$data$matrix))
  expect_true(all(c("facet_summary", "subset_summary") %in% names(p_cov$data)))
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_cov$data)))
  expect_true(is.data.frame(p_cov$data$legend))
  expect_true(is.data.frame(p_cov$data$reference_lines))

  p_design <- plot(sc, type = "design_matrix", draw = FALSE)
  expect_s3_class(p_design, "mfrm_plot_data")
  expect_identical(p_design$data$plot, "coverage_matrix")
  expect_identical(p_design$data$requested_type, "design_matrix")
})

test_that("subset_connectivity linking_matrix draws without error", {
  sc <- subset_connectivity_report(.fit)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(sc, type = "linking_matrix", preset = "publication"))
  expect_no_error(plot(sc, type = "design_matrix", preset = "publication"))
})

# ---- facet_statistics_report ----

test_that("facet_statistics_report returns a bundle", {
  fs <- facet_statistics_report(.fit, diagnostics = .diag)
  expect_s3_class(fs, "mfrm_bundle")
  expect_true(all(c("precision_summary", "variability_tests", "se_modes") %in% names(fs)))
  expect_true(is.data.frame(fs$precision_summary))
  expect_true(is.data.frame(fs$variability_tests))
  expect_true(is.data.frame(fs$se_modes))
  s <- summary(fs)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("precision_audit_report returns a bundle", {
  pa <- precision_audit_report(.fit, diagnostics = .diag)
  expect_s3_class(pa, "mfrm_bundle")
  expect_true(all(c("profile", "checks", "approximation_notes", "settings") %in% names(pa)))
  expect_true(is.data.frame(pa$profile))
  expect_true(is.data.frame(pa$checks))
  expect_true(is.data.frame(pa$approximation_notes))
  s <- summary(pa)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("precision_audit_report marks JML runs as exploratory", {
  pa <- precision_audit_report(.fit, diagnostics = .diag)
  expect_identical(as.character(pa$profile$PrecisionTier[1]), "exploratory")
  expect_true(any(pa$checks$Status %in% c("review", "warn")))
})

test_that("build_precision_profile demotes MML runs with fallback SE to hybrid tier", {
  mock_fit <- .fit
  mock_fit$summary$Method[1] <- "MML"
  mock_fit$config$method <- "MML"

  measure_df <- data.frame(
    Facet = c("Person", "Rater"),
    Level = c("P1", "R1"),
    SE_Method = c("Posterior SD (EAP)", "Fallback observation-table information"),
    RealSE = c(0.4, 0.5),
    stringsAsFactors = FALSE
  )

  profile <- mfrmr:::build_precision_profile(
    res = mock_fit,
    measure_df = measure_df,
    reliability_tbl = data.frame(),
    facet_precision_tbl = data.frame()
  )

  expect_identical(as.character(profile$PrecisionTier[1]), "hybrid")
  expect_false(isTRUE(profile$SupportsFormalInference[1]))
  expect_true(isTRUE(profile$HasFallbackSE[1]))
})

test_that("facet_statistics_report filters precision summary by basis and SE mode", {
  fs <- facet_statistics_report(
    .fit,
    diagnostics = .diag,
    distribution_basis = "population",
    se_mode = "model"
  )
  expect_true(all(fs$precision_summary$DistributionBasis == "population"))
  expect_true(all(fs$precision_summary$SEMode == "model"))
})

# ---- unexpected_response_table ----

test_that("unexpected_response_table produces valid output", {
  ut <- unexpected_response_table(.fit, diagnostics = .diag)
  expect_s3_class(ut, "mfrm_unexpected")
  expect_true(all(c("table", "summary", "thresholds") %in% names(ut)))
  s <- summary(ut)
  expect_s3_class(s, "summary.mfrm_bundle")
  p <- plot(ut, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# ---- fair_average_table ----

test_that("fair_average_table produces valid output", {
  fa <- fair_average_table(.fit, diagnostics = .diag)
  expect_s3_class(fa, "mfrm_fair_average")
  expect_true("stacked" %in% names(fa))
  s <- summary(fa)
  expect_s3_class(s, "summary.mfrm_bundle")
  p <- plot(fa, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

test_that("fair_average_table supports native label and reference controls", {
  fa <- fair_average_table(.fit, diagnostics = .diag, reference = "mean", label_style = "native")
  expect_true(all(c("AdjustedAverage", "ObservedAverage", "ModelBasedSE", "FitAdjustedSE") %in%
    names(fa$stacked)))
  expect_false(any(c("Fair(M) Average", "Fair(Z) Average", "Obsvd Average", "Model S.E.", "Real S.E.") %in%
    names(fa$stacked)))
  expect_false("StandardizedAdjustedAverage" %in% names(fa$stacked))
})

# ---- displacement_table ----

test_that("displacement_table produces valid output", {
  dt <- displacement_table(.fit, diagnostics = .diag)
  expect_s3_class(dt, "mfrm_displacement")
  s <- summary(dt)
  expect_s3_class(s, "summary.mfrm_bundle")
  p <- plot(dt, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

test_that("estimate_bias exposes screening-tier metadata", {
  expect_true(all(c(
    "InferenceTier", "SupportsFormalInference", "FormalInferenceEligible",
    "PrimaryReportingEligible", "ReportingUse",
    "SEBasis", "ProbabilityMetric", "DFBasis", "StatisticLabel"
  ) %in% names(.bias$table)))
  expect_true(all(.bias$table$InferenceTier == "screening"))
  expect_true(all(!.bias$table$SupportsFormalInference))
  expect_true(all(!.bias$table$FormalInferenceEligible))
  expect_true(all(!.bias$table$PrimaryReportingEligible))
  expect_true(all(.bias$table$ReportingUse == "screening_only"))
  expect_true(all(.bias$table$StatisticLabel == "screening t"))
  expect_true(all(c(
    "InferenceTier", "SupportsFormalInference", "FormalInferenceEligible",
    "PrimaryReportingEligible", "ReportingUse", "TestBasis"
  ) %in% names(.bias$chi_sq)))
  expect_true(all(.bias$chi_sq$InferenceTier == "screening"))
  expect_true(all(.bias$chi_sq$ReportingUse == "screening_only"))
})

# ---- measurable_summary_table ----

test_that("measurable_summary_table produces valid output", {
  ms <- measurable_summary_table(.fit, diagnostics = .diag)
  expect_s3_class(ms, "mfrm_bundle")
  s <- summary(ms)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- rating_scale_table ----

test_that("rating_scale_table produces valid output", {
  rs <- rating_scale_table(.fit, diagnostics = .diag)
  expect_s3_class(rs, "mfrm_rating_scale")
  s <- summary(rs)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("rating_scale_table computes PCM threshold gaps within each step facet", {
  toy <- load_mfrmr_data("example_core")
  fit_pcm <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person",
      c("Rater", "Criterion"),
      "Score",
      method = "JML",
      model = "PCM",
      step_facet = "Rater",
      maxit = 20
    )
  )
  rs <- rating_scale_table(fit_pcm)

  expect_true("StepFacet" %in% names(rs$threshold_table))
  split_tbl <- split(rs$threshold_table, rs$threshold_table$StepFacet)
  expect_true(all(vapply(split_tbl, function(tbl) is.na(tbl$GapFromPrev[1]), logical(1))))
  expect_true(is.logical(rs$summary$ThresholdMonotonic) || is.na(rs$summary$ThresholdMonotonic))
})

# ---- bias_count_table ----

test_that("bias_count_table produces valid output", {
  bc <- bias_count_table(.bias)
  expect_s3_class(bc, "mfrm_bundle")
  s <- summary(bc)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- unexpected_after_bias_table ----

test_that("unexpected_after_bias_table produces valid output", {
  ub <- unexpected_after_bias_table(.fit, bias_results = .bias, diagnostics = .diag)
  expect_s3_class(ub, "mfrm_bundle")
})

# ---- interrater_agreement_table ----

test_that("interrater_agreement_table produces valid output", {
  ia <- interrater_agreement_table(.fit, diagnostics = .diag)
  expect_s3_class(ia, "mfrm_interrater")
  expect_true(all(c("OpportunityCount", "ExactCount", "ExpectedExactCount", "AdjacentCount") %in%
    names(ia$pairs)))
  expect_true(all(c("AgreementMinusExpected", "RaterSeparation", "RaterReliability") %in%
    names(ia$summary)))
  s <- summary(ia)
  expect_s3_class(s, "summary.mfrm_bundle")
  p <- plot(ia, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# ---- facets_chisq_table ----

test_that("facets_chisq_table produces valid output", {
  fc <- facets_chisq_table(.fit, diagnostics = .diag)
  expect_s3_class(fc, "mfrm_facets_chisq")
  s <- summary(fc)
  expect_s3_class(s, "summary.mfrm_bundle")
  p <- plot(fc, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# ---- bias_interaction_report ----

test_that("bias_interaction_report produces valid output", {
  bi <- bias_interaction_report(.fit, diagnostics = .diag,
                                facet_a = "Rater", facet_b = "Task")
  expect_s3_class(bi, "mfrm_bundle")
  s <- summary(bi)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("bias_iteration_report produces valid output", {
  bi <- bias_iteration_report(.bias)
  expect_s3_class(bi, "mfrm_bundle")
  expect_true(all(c("table", "summary", "orientation_audit", "settings") %in% names(bi)))
  expect_true(is.data.frame(bi$table))
  expect_true(is.data.frame(bi$summary))
  expect_true(is.data.frame(bi$orientation_audit))
  s <- summary(bi)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("bias_pairwise_report produces valid output", {
  bp <- bias_pairwise_report(.bias, top_n = 8)
  expect_s3_class(bp, "mfrm_bundle")
  expect_true(all(c("table", "summary", "orientation_audit", "settings") %in% names(bp)))
  expect_true(is.data.frame(bp$table))
  expect_true(is.data.frame(bp$summary))
  expect_true(is.data.frame(bp$orientation_audit))
  if (nrow(bp$table) > 0) {
    expect_true(all(c("ContrastBasis", "SEBasis", "StatisticLabel", "ProbabilityMetric", "DFBasis") %in% names(bp$table)))
    expect_true(all(bp$table$StatisticLabel == "Bias-contrast Welch screening t"))

    tgt_se_sq <- bp$table$`Target S.E.`^2
    bias1 <- bp$table$`Local Measure1` - bp$table$`Target Measure`
    bias2 <- bp$table$`Local Measure2` - bp$table$`Target Measure`
    bias_se1_sq <- pmax(bp$table$SE1^2 - tgt_se_sq, 0)
    bias_se2_sq <- pmax(bp$table$SE2^2 - tgt_se_sq, 0)
    expected_contrast <- bias1 - bias2
    expected_se <- sqrt(bias_se1_sq + bias_se2_sq)
    naive_se <- sqrt(bp$table$SE1^2 + bp$table$SE2^2)

    expect_equal(bp$table$Contrast, expected_contrast, tolerance = 1e-8)
    expect_equal(bp$table$SE, expected_se, tolerance = 1e-8)
    expect_true(all(bp$table$SE <= naive_se + 1e-10, na.rm = TRUE))
  }
  s <- summary(bp)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("estimate_bias surfaces optimization failures instead of using zero-bias fallback", {
  toy <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    toy, "Person", c("Rater", "Criterion"), "Score",
    method = "JML", maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))

  testthat::local_mocked_bindings(
    optimize = function(...) stop("forced optimize failure"),
    .package = "stats"
  )

  bias <- estimate_bias(
    fit, diag,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 1
  )

  expect_true(all(bias$table$OptimizationStatus == "failed"))
  expect_true(all(grepl("forced optimize failure", bias$table$OptimizationDetail, fixed = TRUE)))
  expect_true(all(is.na(bias$table$`Obs-Exp Average`)))
  expect_true(all(is.na(bias$table$`S.E.`)))
  expect_true(is.data.frame(bias$optimization_failures))
  expect_true(nrow(bias$optimization_failures) > 0)
})

test_that("plot_bias_interaction treats non-finite scatter and ranked inputs as no-data", {
  bundle <- list(
    ranked_table = data.frame(Pair = "A vs B", BiasSize = NA_real_, Flag = FALSE, stringsAsFactors = FALSE),
    scatter_data = data.frame(ObsExpAverage = NA_real_, BiasSize = NA_real_, Flag = FALSE, t = NA_real_, stringsAsFactors = FALSE),
    summary = data.frame(MeanAbsBias = NA_real_, PctFlagged = NA_real_, stringsAsFactors = FALSE),
    thresholds = list(abs_bias_warn = 0.2, abs_t_warn = 2)
  )

  dev_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(dev_path)
  on.exit({
    grDevices::dev.off()
    unlink(dev_path)
  }, add = TRUE)

  expect_silent(plot_bias_interaction(bundle, plot = "scatter", draw = TRUE))
  expect_silent(plot_bias_interaction(bundle, plot = "ranked", draw = TRUE))
})

test_that("bias reports flag mixed-sign orientation when facets mix score directions", {
  fit_pos <- suppressWarnings(
    fit_mfrm(
      mfrmr:::sample_mfrm_data(seed = 7),
      "Person",
      c("Rater", "Task", "Criterion"),
      "Score",
      method = "JML",
      maxit = 15,
      positive_facets = "Rater"
    )
  )
  diag_pos <- diagnose_mfrm(fit_pos, residual_pca = "none")
  bi <- bias_iteration_report(fit_pos, diagnostics = diag_pos, facet_a = "Rater", facet_b = "Task", max_iter = 2)
  expect_true(isTRUE(bi$summary$MixedSign[1]))
  expect_true(any(bi$orientation_audit$Orientation == "positive"))
  expect_true(any(bi$orientation_audit$Orientation == "negative"))
  expect_match(bi$direction_note, "higher-than-expected|lower-than-expected")
})

# ---- build_apa_outputs ----

test_that("build_apa_outputs produces structured APA text", {
  apa <- build_apa_outputs(.fit, diagnostics = .diag)
  expect_s3_class(apa, "mfrm_apa_outputs")
  expect_true("report_text" %in% names(apa))
  expect_true("section_map" %in% names(apa))
  expect_true(nchar(apa$report_text) > 50)
  s <- summary(apa)
  expect_s3_class(s, "summary.mfrm_apa_outputs")
  expect_true(is.data.frame(s$sections))
  expect_true("DraftContractPass" %in% names(s$overview))
  expect_true(any(grepl("contract completeness", s$notes, fixed = TRUE)))
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
})

test_that("build_apa_outputs with bias produces extended text", {
  apa <- build_apa_outputs(.fit, diagnostics = .diag, bias = .bias)
  expect_true(nchar(apa$report_text) > 100)
})

# ---- build_fixed_reports ----

test_that("build_fixed_reports produces text reports", {
  fr <- build_fixed_reports(.bias)
  expect_true(is.list(fr))
  expect_true(length(fr) > 0)
})

test_that("build_fixed_reports pvalue plot degrades gracefully when p-values are unavailable", {
  fr <- build_fixed_reports(.bias)
  fr$pairwise_table$`Prob.` <- NA_character_

  p <- plot(fr, type = "pvalue", draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$data$plot, "pvalue")
  expect_length(p$data$p_values, 0)
  expect_match(p$data$message, "No finite p-values available", fixed = TRUE)
})

# ---- build_visual_summaries ----

test_that("build_visual_summaries produces warning and summary maps", {
  vs <- build_visual_summaries(.fit, diagnostics = .diag)
  expect_true(is.list(vs))
  expect_true("warning_map" %in% names(vs) || "summary_map" %in% names(vs))
  expect_true("category_probability_surface" %in% names(vs$plot_payloads))
  expect_s3_class(vs$plot_payloads$category_probability_surface, "mfrm_plot_data")
  expect_true("category_probability_surface" %in% vs$public_plot_routes$Visual)
})

# ---- apa_table ----

test_that("apa_table produces structured output", {
  at <- apa_table(.fit, diagnostics = .diag)
  expect_s3_class(at, "apa_table")
  s <- summary(at)
  expect_s3_class(s, "summary.apa_table")
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
})

# ---- analyze_residual_pca ----

test_that("analyze_residual_pca produces eigenvalue and loading output", {
  pca <- analyze_residual_pca(.diag, mode = "both")
  expect_s3_class(pca, "mfrm_residual_pca")
  s <- summary(pca)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("analyze_residual_pca accepts fit object directly", {
  pca <- analyze_residual_pca(.fit, mode = "overall")
  expect_s3_class(pca, "mfrm_residual_pca")
})

test_that("analyze_residual_pca retains computation errors instead of dropping them", {
  local_mocked_bindings(
    compute_pca_overall = function(...) {
      list(
        pca = NULL,
        residual_matrix = NULL,
        cor_matrix = NULL,
        error = "forced PCA failure"
      )
    },
    .package = "mfrmr"
  )

  pca <- analyze_residual_pca(.diag, mode = "overall", facets = "Rater")

  expect_s3_class(pca, "mfrm_residual_pca")
  expect_equal(nrow(pca$overall_table), 0)
  expect_match(pca$errors$overall, "forced PCA failure", fixed = TRUE)
})

# ---- plot_residual_pca ----

test_that("plot_residual_pca produces plot bundles", {
  pca <- analyze_residual_pca(.diag, mode = "overall")
  p_scree <- plot_residual_pca(pca, plot_type = "scree", draw = FALSE)
  expect_s3_class(p_scree, "mfrm_plot_data")
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p_scree$data)))
  expect_true(is.data.frame(p_scree$data$legend))
  expect_true(is.data.frame(p_scree$data$reference_lines))
})

# ---- plot.mfrm_fit specific types ----

test_that("plot.mfrm_fit supports all named types", {
  p_wright <- plot(.fit, type = "wright", draw = FALSE)
  expect_s3_class(p_wright, "mfrm_plot_data")
  expect_true(all(c("person_hist", "person_stats", "label_points", "group_summary", "y_range") %in% names(p_wright$data)))

  p_pathway <- plot(.fit, type = "pathway", draw = FALSE)
  expect_s3_class(p_pathway, "mfrm_plot_data")
  expect_true(all(c("steps", "endpoint_labels", "dominance_regions") %in% names(p_pathway$data)))

  p_ccc <- plot(.fit, type = "ccc", draw = FALSE)
  expect_s3_class(p_ccc, "mfrm_plot_data")

  p_person <- plot(.fit, type = "person", draw = FALSE)
  expect_s3_class(p_person, "mfrm_plot_data")

  p_step <- plot(.fit, type = "step", draw = FALSE)
  expect_s3_class(p_step, "mfrm_plot_data")
})

test_that("plot_wright_unified returns enhanced Wright-map payload", {
  p_wright_unified <- plot_wright_unified(.fit, draw = FALSE, preset = "publication", show_thresholds = FALSE)
  expect_true(all(c("persons", "facets", "person_hist", "person_stats", "group_summary", "y_lim") %in%
    names(p_wright_unified)))
  expect_null(p_wright_unified$thresholds)
})

# ---- plot_qc_dashboard ----

test_that("plot_qc_dashboard returns a plot bundle", {
  p <- plot_qc_dashboard(.fit, diagnostics = .diag, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(as.character(p$data$preset), "standard")
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(p$data)))

  p_pub <- plot_qc_dashboard(.fit, diagnostics = .diag, draw = FALSE, preset = "publication")
  expect_identical(as.character(p_pub$data$preset), "publication")
})

# ---- make_anchor_table ----

test_that("make_anchor_table extracts anchors from fitted model", {
  at <- make_anchor_table(.fit)
  expect_true(is.data.frame(at))
  expect_true(all(c("Facet", "Level") %in% names(at)))
})

test_that("make_anchor_table includes persons when requested", {
  at <- make_anchor_table(.fit, include_person = TRUE)
  expect_true("Person" %in% at$Facet || nrow(at) > 0)
})

# ---- Formatting helpers (internal) ----

test_that("py_style_format converts Python-style format strings", {
  fmt <- mfrmr:::py_style_format
  expect_equal(fmt("{:.2f}", 3.14159), "3.14")
  expect_equal(fmt("{:.0f}", 42.7), "43")
})

test_that("fmt_num formats numbers correctly", {
  fn <- mfrmr:::fmt_num
  expect_equal(fn(3.14159, 2), "3.14")
  expect_equal(fn(NA, 2), "NA")
})

test_that("fmt_count formats integers correctly", {
  fc <- mfrmr:::fmt_count
  expect_equal(fc(42), "42")
  expect_equal(fc(NA), "NA")
})

test_that("fmt_pvalue formats p-values correctly", {
  fp <- mfrmr:::fmt_pvalue
  expect_true(grepl("< .001", fp(0.0001)))
  expect_true(grepl("= ", fp(0.05)))
})
