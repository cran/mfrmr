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
  expect_true("quality_overview" %in% names(dq))
  expect_true("quality_flags" %in% names(dq))
  expect_true("facet_response_patterns" %in% names(dq))
  expect_false("attention_items" %in% names(dq))
  expect_true(all(c("Area", "Status", "NextStep") %in% names(dq$quality_overview)))
  s <- summary(dq)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_error(
    data_quality_report(.fit, dominant_category_cutoff = 0),
    "`dominant_category_cutoff`"
  )
})

test_that("data_quality_report flags retained zero-frequency score categories", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d <- d[d$Score != 3, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = 8,
      rating_min = 1, rating_max = 5,
      keep_original = TRUE
    )
  )
  dq <- data_quality_report(
    fit,
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score"
  )
  zero_row <- dq$category_counts[dq$category_counts$Score == 3, , drop = FALSE]
  expect_equal(nrow(zero_row), 1L)
  expect_true(isTRUE(zero_row$ZeroCount[1]))
  expect_identical(zero_row$UnusedCategoryType[1], "internal")
  expect_true("zero_count_intermediate_score_category" %in% dq$caveats$Condition)
  expect_true(any(dq$quality_flags$Flag == "Intermediate score categories have zero observations"))
  expect_true(any(dq$quality_overview$Area == "Score support" & dq$quality_overview$Status == "high"))
  expect_identical(summary(dq)$preview_name, "quality_flags")

  p <- plot(dq, type = "score_support", preset = "monochrome", draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$data$plot, "score_support")
  expect_equal(p$data$preset, "monochrome")
  p_flags <- plot(dq, type = "quality_flags", draw = FALSE)
  expect_s3_class(p_flags, "mfrm_plot_data")
  expect_identical(p_flags$data$plot, "quality_flags")
})

test_that("data_quality_report flags original gaps hidden by score recoding", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d <- d[d$Score != 3, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = 8,
      rating_min = 1, rating_max = 5,
      keep_original = FALSE
    )
  )
  dq <- data_quality_report(
    fit,
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score"
  )
  expect_true("score_categories_recoded" %in% dq$caveats$Condition)
  expect_true("original_score_gap_before_recoding" %in% dq$caveats$Condition)
  expect_equal(
    dq$caveats$Categories[dq$caveats$Condition == "original_score_gap_before_recoding"],
    "3"
  )
  expect_true(any(dq$quality_flags$Flag == "Original score sequence had gaps before recoding"))
  expect_true(any(dq$quality_overview$Area == "Score support" & dq$quality_overview$Status == "high"))
  p_map <- plot(dq, type = "score_map", draw = FALSE)
  expect_s3_class(p_map, "mfrm_plot_data")
  expect_identical(p_map$data$plot, "score_map")
  expect_true(any(p_map$data$table$MappingStatus == "recoded"))
  expect_false("score_out_of_range" %in% dq$row_review$Status)
})

test_that("data_quality_report flags facet-level category usage gaps", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score[d$Rater == "R1" & d$Score == 3] <- 2L
  expect_true(any(d$Score == 3))
  expect_false(any(d$Rater == "R1" & d$Score == 3))

  fit <- suppressWarnings(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = 8,
      rating_min = 1, rating_max = 5,
      keep_original = TRUE
    )
  )
  dq <- data_quality_report(
    fit,
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    min_category_count = 10
  )

  r1_mid <- dq$category_usage_by_facet[
    dq$category_usage_by_facet$Facet == "Rater" &
      dq$category_usage_by_facet$Level == "R1" &
      dq$category_usage_by_facet$Score == 3,
    ,
    drop = FALSE
  ]
  expect_equal(nrow(r1_mid), 1L)
  expect_true(isTRUE(r1_mid$ZeroCount[1]))
  expect_identical(r1_mid$CategoryPosition[1], "intermediate")
  expect_identical(r1_mid$ReviewStatus[1], "warning")

  r1_summary <- dq$category_usage_summary[
    dq$category_usage_summary$Facet == "Rater" &
      dq$category_usage_summary$Level == "R1",
    ,
    drop = FALSE
  ]
  expect_equal(r1_summary$IntermediateZeroCategories[1], 1)
  expect_identical(r1_summary$ReviewStatus[1], "warning")
  expect_true(dq$summary$FacetLevelsWithIntermediateZeroCategories[1] >= 1)
  expect_true(any(dq$quality_flags$Flag == "Facet levels have intermediate zero-category use"))
  expect_true(any(dq$quality_overview$Area == "Facet category use" & dq$quality_overview$Status == "high"))

  p <- plot(dq, type = "facet_category_usage", top_n = 5, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$data$plot, "facet_category_usage")
  expect_equal(p$data$top_n, 5L)
})

test_that("data_quality_report flags facet levels with restricted response patterns", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score[d$Rater == "R1"] <- 1L

  fit <- suppressWarnings(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = 8,
      rating_min = 1, rating_max = 5,
      keep_original = TRUE
    )
  )
  dq <- data_quality_report(
    fit,
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score"
  )

  r1_pattern <- dq$facet_response_patterns[
    dq$facet_response_patterns$Facet == "Rater" &
      dq$facet_response_patterns$Level == "R1",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(r1_pattern), 1L)
  expect_true(isTRUE(r1_pattern$SingleCategoryUse[1]))
  expect_equal(r1_pattern$DominantScore[1], 1L)
  expect_identical(r1_pattern$PatternStatus[1], "high")
  expect_true(any(dq$quality_flags$Flag == "Facet levels use only one score category"))
  expect_true(any(dq$quality_overview$Area == "Facet response patterns" & dq$quality_overview$Status == "high"))

  p <- plot(dq, type = "facet_response_patterns", top_n = 5, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$data$plot, "facet_response_patterns")
})

test_that("fit_measures_table lists directional fit statuses", {
  fm <- fit_measures_table(.fit, diagnostics = .diag, facet = "Rater", ci_level = 0.90)
  expect_s3_class(fm, "mfrm_fit_measures")
  expect_true(all(c(
    "Facet", "Level", "Infit", "Outfit", "FitStatus",
    "Underfit", "Overfit", "ReviewReason", "CI_Lower", "CI_Upper", "CI_Level"
  ) %in% names(fm$table)))
  expect_true(all(c("Infit MnSq", "Outfit MnSq", "Lower CI", "Upper CI", "CI Level", "Fit Status") %in% names(fm$facets_table)))
  expect_true(all(fm$table$Facet == "Rater"))
  expect_true(all(fm$table$FitStatus %in% c("underfit", "overfit", "mixed", "within_band", "not_available")))
  expect_equal(unique(fm$table$CI_Level), 0.90)
  expect_true(all(c(
    "Profile", "ProfileLabel", "Lower", "Upper", "Facet",
    "UnderfitRate", "OverfitRate", "AnyFlagRate"
  ) %in% names(fm$profile_summary_by_facet)))
  expect_gt(length(unique(fm$profile_summary_by_facet$Profile)), 1)
  expect_true(all(fm$profile_summary_by_facet$Facet == "Rater"))

  s <- summary(fm)
  expect_s3_class(s, "summary.mfrm_bundle")
  p_status <- plot(fm, draw = FALSE)
  expect_s3_class(p_status, "mfrm_plot_data")
  expect_identical(p_status$data$plot, "status")
  p_scatter <- plot(fm, type = "infit_outfit", draw = FALSE)
  expect_s3_class(p_scatter, "mfrm_plot_data")
  expect_identical(p_scatter$data$plot, "infit_outfit")
  p_ci <- plot(fm, type = "measure_ci", ci_level = 0.80, preset = "monochrome", draw = FALSE)
  expect_s3_class(p_ci, "mfrm_plot_data")
  expect_identical(p_ci$data$plot, "measure_ci")
  expect_equal(p_ci$data$ci_level, 0.80)
  expect_equal(p_ci$data$preset, "monochrome")

  fm_df <- fit_measures_table(
    .fit,
    facet = "Rater",
    fit_df_method = "both",
    df_zstd_tolerance = 0.01,
    df_zstd_large_shift = 0.25,
    df_ratio_tolerance = 0.01,
    top_n = Inf
  )
  expect_true(all(c(
    "Infit df", "Outfit df", "Fit df method",
    "FACETS Infit df", "FACETS Outfit df",
    "FACETS Infit ZStd", "FACETS Outfit ZStd",
    "Max df rel shift"
  ) %in% names(fm_df$facets_table)))
  expect_true(all(c(
    "DF_Infit_ENGINE", "DF_Infit_FACETS",
    "InfitZSTD_FACETS", "FitDfMethod", "FitZSTDTransform"
  ) %in% names(fm_df$table)))
  expect_s3_class(fm_df$df_conversion_guide, "mfrm_facets_fit_df_guide")
  expect_equal(fm_df$settings$fit_df_method, "both")
  expect_true(all(c(
    "Facet", "Level", "InfitZSTD_ENGINE", "InfitZSTD_FACETS",
    "MaxAbsZSTDDiff_FACETS_vs_ENGINE", "FlagChangedByDf",
    "MaxDFRelativeDifference_ENGINE_vs_FACETS",
    "DfSensitivityStatus", "Interpretation"
  ) %in% names(fm_df$df_sensitivity)))
  expect_true(all(c(
    "ComparedRows", "FlagChangedByDfRows", "LargeZSTDShiftRows",
    "DfConventionDifferenceRows"
  ) %in% names(fm_df$df_sensitivity_summary)))
  expect_true(all(c(
    "DfComparedRows", "DfSensitiveRows", "FlagChangedByDfRows"
  ) %in% names(fm_df$summary)))
  expect_equal(fm_df$summary$DfSensitiveRows, nrow(fm_df$df_sensitive))
  expect_equal(fm_df$settings$df_zstd_tolerance, 0.01)
  expect_equal(fm_df$settings$df_zstd_large_shift, 0.25)
  expect_equal(fm_df$settings$df_ratio_tolerance, 0.01)
  s_df <- summary(fm_df)
  expect_s3_class(s_df, "summary.mfrm_bundle")
  expect_identical(s_df$summary_kind, "mfrm_fit_measures")
  p_df <- plot(fm_df, type = "df_sensitivity", top_n = 5, draw = FALSE)
  expect_s3_class(p_df, "mfrm_plot_data")
  expect_identical(p_df$data$plot, "df_sensitivity")
  expect_lte(nrow(p_df$data$table), 5L)

  fm_df_from_engine_diag <- fit_measures_table(
    .fit,
    diagnostics = .diag,
    facet = "Rater",
    fit_df_method = "both",
    top_n = Inf
  )
  expect_true("FACETS Infit df" %in% names(fm_df_from_engine_diag$facets_table))
  expect_error(
    fit_measures_table(.fit, facet = "Rater", fit_df_method = "both", df_zstd_large_shift = 0.01),
    "`df_zstd_large_shift`"
  )
})

test_that("fit-measure threshold profile rates use all selected rows", {
  fm_full <- fit_measures_table(.fit, diagnostics = .diag, facet = "Rater", top_n = Inf)
  fm_top <- fit_measures_table(.fit, diagnostics = .diag, facet = "Rater", top_n = 1)
  expect_equal(nrow(fm_top$table), 1L)
  expect_equal(fm_top$summary$Rows, fm_full$summary$Rows)
  expect_equal(fm_top$summary$DisplayedRows, 1L)
  expect_equal(fm_top$profile_summary_by_facet, fm_full$profile_summary_by_facet)
})

# ---- category_curves_report ----

test_that("category_curves_report returns a bundle", {
  cc <- category_curves_report(.fit)
  expect_s3_class(cc, "mfrm_bundle")
  cc_summary <- summary(cc)
  expect_s3_class(cc_summary, "summary.mfrm_bundle")
  expect_true(all(c("Metric", "Value") %in% names(cc_summary$summary)))
  expect_true(all(c(
    "cumulative_probability_rows",
    "cumulative_boundary_rows",
    "category_information_rows",
    "boundary_rows_needing_review"
  ) %in% cc_summary$summary$Metric))
  expect_true("category_information" %in% names(cc))
  expect_true(all(c("CategoryInformation", "CategoryInformationShare") %in%
                    names(cc$category_information)))
  expect_true(all(c("cumulative_probabilities", "cumulative_boundaries") %in% names(cc)))
  expect_true(all(c("Direction", "CumulativeProbability", "BoundaryCategory") %in%
                    names(cc$cumulative_probabilities)))
  expect_true(any(cc$cumulative_probabilities$Direction == "at_or_below"))
  expect_true(any(cc$cumulative_probabilities$Direction == "at_or_above"))
})

test_that("standalone residual and subset writers create files", {
  residual_path <- tempfile(fileext = ".csv")
  residual <- write_mfrm_residual_file(
    .fit,
    diagnostics = .diag,
    path = residual_path,
    overwrite = TRUE,
    include_probabilities = TRUE
  )
  expect_s3_class(residual, "mfrm_residual_file")
  expect_true(file.exists(residual_path))
  expect_true(all(c("Observed", "Expected", "Residual", "StdResidual") %in% names(residual$table)))
  expect_true(any(grepl("^PrCategory_", names(residual$table))))
  expect_true(all(c("Component", "Format", "Path") %in% names(residual$written_files)))

  subset_path <- tempfile(fileext = ".csv")
  subset <- write_mfrm_subset_file(
    .fit,
    diagnostics = .diag,
    path = subset_path,
    overwrite = TRUE
  )
  expect_s3_class(subset, "mfrm_subset_file")
  expect_true(file.exists(subset_path))
  node_path <- subset$written_files$Path[subset$written_files$Component == "subset_nodes"]
  expect_length(node_path, 1L)
  expect_true(file.exists(node_path))
  expect_true(all(c("Subset", "Observations") %in% names(subset$table)))
  expect_true(all(c("Node", "Subset", "Facet", "Level") %in% names(subset$nodes)))

  subset_tsv_path <- tempfile(fileext = ".tsv")
  subset_tsv <- write_mfrm_subset_file(
    .fit,
    diagnostics = .diag,
    path = subset_tsv_path,
    overwrite = TRUE,
    include_nodes = FALSE
  )
  expect_identical(subset_tsv$summary$Format[1], "tsv")
  expect_identical(subset_tsv$settings$include_nodes, FALSE)
  expect_error(
    write_mfrm_subset_file(
      .fit,
      diagnostics = .diag,
      path = subset_path,
      node_path = subset_path,
      overwrite = TRUE
    ),
    "`node_path` must differ"
  )
  expect_error(
    write_mfrm_residual_file(
      .fit,
      diagnostics = .diag,
      path = residual_path,
      overwrite = FALSE
    ),
    "already exists"
  )
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

  p_net <- plot(sc, type = "network", draw = FALSE)
  expect_s3_class(p_net, "mfrm_plot_data")
  expect_identical(p_net$data$plot, "network")
  expect_true(all(c("nodes", "edges") %in% names(p_net$data)))
  expect_s3_class(p_net$data$nodes, "data.frame")
  expect_s3_class(p_net$data$edges, "data.frame")
  expect_true(all(c("From", "To", "Weight") %in% names(p_net$data$edges)))
})

test_that("subset_connectivity linking_matrix draws without error", {
  sc <- subset_connectivity_report(.fit)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(sc, type = "linking_matrix", preset = "publication"))
  expect_no_error(plot(sc, type = "design_matrix", preset = "publication"))
})

test_that("mfrm_network_analysis returns graph metrics for the fitted design", {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    skip("igraph (Suggests) not installed.")
  }
  net <- mfrm_network_analysis(.fit, diagnostics = .diag)
  expect_s3_class(net, "mfrm_bundle")
  expect_s3_class(net, "mfrm_network_analysis")
  expect_true(all(c(
    "summary", "node_metrics", "edge_metrics", "facet_summary",
    "cut_nodes", "bridge_edges"
  ) %in% names(net)))
  expect_true(all(c(
    "Nodes", "Edges", "Components", "ArticulationPoints", "Bridges",
    "Connected"
  ) %in% names(net$summary)))
  expect_true(all(c(
    "Node", "Facet", "Degree", "Strength", "Betweenness",
    "IsArticulationPoint"
  ) %in% names(net$node_metrics)))
  expect_true(all(c("From", "To", "Weight", "EdgeBetweenness", "IsBridge")
                  %in% names(net$edge_metrics)))
  expect_true(all(c("Facet", "Levels", "ArticulationPoints", "BridgeIncidentEdges")
                  %in% names(net$facet_summary)))
  expect_true(all(net$node_metrics$Degree >= 0))
  p_cent <- plot(net, type = "centrality", draw = FALSE)
  expect_s3_class(p_cent, "mfrm_plot_data")
  expect_identical(p_cent$data$plot, "centrality")
  p_facet <- plot(net, type = "facet_summary", metric = "BridgeIncidentEdges", draw = FALSE)
  expect_identical(p_facet$data$plot, "facet_summary")
  p_network <- plot(net, type = "network", draw = FALSE)
  expect_identical(p_network$data$plot, "network")
  s <- summary(net)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("build_mfrm_network_review keeps design-network diagnostics separate from measurement results", {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    skip("igraph (Suggests) not installed.")
  }

  sparse_spec <- build_mfrm_sim_spec(
    n_person = 16,
    n_rater = 4,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "sparse_linked",
    sparse_controls = list(
      link_persons = 3,
      link_raters_per_person = 4,
      min_common_persons_per_rater_pair = 2
    )
  )
  sparse_sim <- simulate_mfrm_data(sim_spec = sparse_spec, seed = 20260526)
  peer_spec <- build_peer_review_sim_spec(
    n_submission = 8,
    n_criterion = 2,
    reviewers_per_submission = 2,
    anchor_submissions = 1
  )
  peer_sim <- simulate_mfrm_data(sim_spec = peer_spec, seed = 20260526)

  review <- build_mfrm_network_review(
    .fit,
    diagnostics = .diag,
    sparse_design = sparse_sim,
    peer_review_design = peer_sim,
    top_n = 3
  )

  expect_s3_class(review, "mfrm_bundle")
  expect_s3_class(review, "mfrm_network_review")
  expect_true(all(c(
    "overview", "network_summary", "facet_summary", "top_central_nodes",
    "top_cut_nodes", "top_bridge_edges", "sparse_review", "peer_review",
    "reporting_map", "source_network"
  ) %in% names(review)))
  expect_true(all(c(
    "NetworkReviewStatus", "NetworkReviewReason", "ReviewUse"
  ) %in% names(review$overview)))
  expect_identical(
    review$overview$ReviewUse[1],
    "design_diagnostic_not_measurement_gate"
  )
  expect_true(any(review$reporting_map$Area == "MFRM measurement model"))
  expect_true(any(review$reporting_map$Area == "Peer-review design"))
  expect_true(any(grepl("Design diagnostic", review$reporting_map$Boundary, fixed = TRUE)))
  expect_true(nrow(review$top_central_nodes) <= 3L)
  expect_true(nrow(review$sparse_review) == 1L)
  expect_identical(
    review$sparse_review$ReviewUse[1],
    "design_diagnostic_not_recovery_gate"
  )
  expect_true(nrow(review$peer_review) == 1L)
  expect_identical(
    review$peer_review$ReviewUse[1],
    "design_diagnostic_not_measurement_gate"
  )
  expect_equal(review$peer_review$SelfReviews[1], 0L)

  s <- summary(review, top_n = 2)
  expect_s3_class(s, "summary.mfrm_network_review")
  expect_true(nrow(s$top_central_nodes) <= 2L)
  expect_true(nrow(s$peer_review) == 1L)

  bundle <- build_summary_table_bundle(review)
  expect_identical(bundle$source_class, "mfrm_network_review")
  expect_identical(bundle$summary_class, "summary.mfrm_network_review")
  expect_true(all(c(
    "overview", "network_summary", "facet_summary", "top_central_nodes",
    "sparse_review", "peer_review", "reporting_map"
  ) %in% names(bundle$tables)))
  expect_identical(
    as.character(bundle$table_index$Role[bundle$table_index$Table == "overview"]),
    "network_review_overview"
  )
  expect_identical(
    as.character(bundle$table_index$Role[bundle$table_index$Table == "peer_review"]),
    "peer_review_design_diagnostics"
  )
})

test_that("rater_network_analysis returns rater relationship graph metrics", {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    skip("igraph (Suggests) not installed.")
  }
  rn <- rater_network_analysis(.fit, diagnostics = .diag, mode = "severity_direction")
  expect_s3_class(rn, "mfrm_bundle")
  expect_s3_class(rn, "mfrm_rater_network")
  expect_true(all(c(
    "summary", "node_metrics", "edge_metrics", "pair_metrics",
    "source_interrater", "caveats"
  ) %in% names(rn)))
  expect_true(all(c(
    "Rater", "InStrength", "OutStrength", "SeverityIndex",
    "RelativePattern"
  ) %in% names(rn$node_metrics)))
  expect_true(all(c("From", "To", "Weight", "Direction") %in% names(rn$edge_metrics)))
  expect_true(any(is.finite(rn$node_metrics$SeverityIndex)))
  p_sev <- plot(rn, type = "severity", draw = FALSE)
  expect_s3_class(p_sev, "mfrm_plot_data")
  expect_identical(p_sev$data$plot, "severity")
  p_net <- plot(rn, type = "network", draw = FALSE)
  expect_identical(p_net$data$plot, "network")
  p_mat <- plot(rn, type = "matrix", draw = FALSE)
  expect_identical(p_mat$data$plot, "matrix")
  rn_agree <- rater_network_analysis(.fit, diagnostics = .diag, mode = "agreement")
  expect_identical(rn_agree$summary$Mode[1], "agreement")
  expect_true(all(rn_agree$edge_metrics$Weight >= 0))
  s <- summary(rn)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("rater_halo_network_analysis returns rater-by-criterion halo diagnostics", {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    skip("igraph (Suggests) not installed.")
  }
  halo <- rater_halo_network_analysis(.fit, diagnostics = .diag)
  expect_s3_class(halo, "mfrm_bundle")
  expect_s3_class(halo, "mfrm_halo_network")
  expect_true(all(c(
    "summary", "node_metrics", "edge_metrics", "pair_metrics",
    "halo_summary_by_rater", "caveats"
  ) %in% names(halo)))
  expect_true(all(c(
    "RaterFacet", "CriterionFacet", "HaloEdges", "MeanHaloWeight",
    "MeanNonHaloWeight"
  ) %in% names(halo$summary)))
  expect_true(all(c("Node", "Rater", "Criterion", "HaloStrength") %in% names(halo$node_metrics)))
  expect_true(all(c("ReviewStatus", "ReviewReason", "HaloMinusIncidentNonHalo")
                  %in% names(halo$halo_summary_by_rater)))
  expect_true(all(c("From", "To", "EdgeType", "Estimate", "PAdjusted", "RetainedEdge")
                  %in% names(halo$pair_metrics)))
  expect_true(any(halo$pair_metrics$EdgeType == "halo"))
  expect_true(all(halo$halo_summary_by_rater$ReviewStatus %in%
                    c("warning", "review", "ok", "insufficient_data")))
  p_dist <- plot(halo, type = "edge_distribution", draw = FALSE)
  expect_s3_class(p_dist, "mfrm_plot_data")
  expect_identical(p_dist$data$plot, "edge_distribution")
  p_sum <- plot(halo, type = "halo_summary", draw = FALSE)
  expect_identical(p_sum$data$plot, "halo_summary")
  p_mat <- plot(halo, type = "matrix", draw = FALSE)
  expect_identical(p_mat$data$plot, "matrix")
  p_net <- plot(halo, type = "network", draw = FALSE)
  expect_identical(p_net$data$plot, "network")
  s <- summary(halo)
  expect_s3_class(s, "summary.mfrm_bundle")
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

test_that("precision_review_report returns a bundle", {
  pa <- precision_review_report(.fit, diagnostics = .diag)
  expect_s3_class(pa, "mfrm_bundle")
  expect_true(all(c(
    "profile", "checks", "fit_separation_basis", "approximation_notes", "settings"
  ) %in% names(pa)))
  expect_true(is.data.frame(pa$profile))
  expect_true(is.data.frame(pa$checks))
  expect_true(is.data.frame(pa$fit_separation_basis))
  expect_true(is.data.frame(pa$approximation_notes))
  expect_true(all(c(
    "Topic", "SourceBasis", "PackageSurface", "Interpretation",
    "ValidationUse", "Availability"
  ) %in% names(pa$fit_separation_basis)))
  expect_true(any(grepl("Wright & Linacre (1994)", pa$fit_separation_basis$SourceBasis, fixed = TRUE)))
  expect_true(any(grepl("Wright & Masters (1982)", pa$fit_separation_basis$SourceBasis, fixed = TRUE)))
  expect_true(any(grepl("ZSTD", pa$fit_separation_basis$Topic, fixed = TRUE)))
  expect_true(any(grepl("not a standalone validation success criterion",
                        pa$fit_separation_basis$ValidationUse, fixed = TRUE)))
  s <- summary(pa)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_equal(s$summary$FitSeparationRows[1], nrow(pa$fit_separation_basis))
})

test_that("precision_review_report marks JML runs as exploratory", {
  pa <- precision_review_report(.fit, diagnostics = .diag)
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
  expect_true(all(c(
    "FairSERequested", "FairSEAvailableRows", "FairSEMethod", "FairSEStatus"
  ) %in% names(s$summary)))
  expect_false(isTRUE(s$summary$FairSERequested[1]))
  expect_identical(s$summary$FairSEMethod[1], "not_requested")
  expect_identical(s$summary$FairSEStatus[1], "not_requested")
  expect_true(all(c("ObservedAverage", "AdjustedAverage") %in% names(s$preview)))
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
  expect_true(all(c("table", "summary", "orientation_review", "settings") %in% names(bi)))
  expect_true(is.data.frame(bi$table))
  expect_true(is.data.frame(bi$summary))
  expect_true(is.data.frame(bi$orientation_review))
  s <- summary(bi)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("bias_pairwise_report produces valid output", {
  bp <- bias_pairwise_report(.bias, top_n = 8)
  expect_s3_class(bp, "mfrm_bundle")
  expect_true(all(c("table", "summary", "orientation_review", "settings") %in% names(bp)))
  expect_true(is.data.frame(bp$table))
  expect_true(is.data.frame(bp$summary))
  expect_true(is.data.frame(bp$orientation_review))
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
  expect_true(any(bi$orientation_review$Orientation == "positive"))
  expect_true(any(bi$orientation_review$Orientation == "negative"))
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

test_that("analyze_residual_pca supports residual-permutation parallel analysis", {
  pca <- analyze_residual_pca(
    .diag,
    mode = "both",
    parallel = TRUE,
    parallel_reps = 5,
    seed = 101
  )

  expect_s3_class(pca, "mfrm_residual_pca")
  expect_true(all(c(
    "ParallelMean", "ParallelCutoff", "ExcessOverParallelCutoff",
    "ExceedsParallelCutoff", "SuccessfulParallelReps"
  ) %in% names(pca$overall_table)))
  expect_true(all(c("Component", "ParallelCutoff") %in% names(pca$parallel_overall_table)))
  expect_true("Facet" %in% names(pca$parallel_by_facet_table))
  expect_true(all(c("ParallelAvailable", "SuccessfulParallelReps", "Error", "Warning") %in%
    names(pca$parallel_status)))
  expect_true(any(pca$parallel_status$ParallelAvailable))
  expect_equal(pca$parallel_settings$Enabled, TRUE)
  expect_equal(pca$parallel_settings$Reps, 5L)

  p_scree <- plot_residual_pca(pca, plot_type = "parallel_scree", draw = FALSE)
  expect_s3_class(p_scree, "mfrm_plot_data")
  expect_equal(p_scree$data$plot, "parallel_scree")

  p_excess <- plot_residual_pca(pca, plot_type = "parallel_excess", draw = FALSE)
  expect_s3_class(p_excess, "mfrm_plot_data")
  expect_equal(p_excess$data$plot, "parallel_excess")
})

test_that("plot_residual_pca requires parallel results for parallel plots", {
  pca <- analyze_residual_pca(.diag, mode = "overall")
  expect_error(
    plot_residual_pca(pca, plot_type = "parallel_scree", draw = FALSE),
    "parallel = TRUE",
    fixed = TRUE
  )
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

  p_pathway <- plot(.fit, type = "pathway", diagnostics = .diag, draw = FALSE)
  expect_s3_class(p_pathway, "mfrm_plot_data")
  expect_true(all(c(
    "steps", "endpoint_labels", "dominance_regions",
    "pathway_long", "pathway_annotations", "fit_measures",
    "fit_status", "curve_fit_status", "fit_measure_status"
  ) %in% names(p_pathway$data)))
  expect_true(all(c(
    "Layer", "CurveGroup", "Theta", "Value", "ValueName"
  ) %in% names(p_pathway$data$pathway_long)))
  expect_true(all(c(
    "AnnotationType", "CurveGroup", "X", "Y", "Label"
  ) %in% names(p_pathway$data$pathway_annotations)))
  expect_true(all(c(
    "Facet", "Level", "Measure", "SE", "FitStatus"
  ) %in% names(p_pathway$data$fit_measures)))
  expect_true(all(c(
    "CurveGroup", "FitStatus", "MatchedFitRow"
  ) %in% names(p_pathway$data$curve_fit_status)))
  expect_true(isTRUE(p_pathway$data$fit_measure_status$Available[1]))

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
