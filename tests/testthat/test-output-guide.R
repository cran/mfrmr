test_that("mfrmr_output_guide returns a stable purpose-to-helper map", {
  guide <- mfrmr_output_guide()

  expect_s3_class(guide, "data.frame")
  expect_true(all(c(
    "Scope",
    "Question",
    "OutputFamily",
    "MainFunction",
    "UseWhen",
    "TypicalInput",
    "NextStep",
    "GPCMStatus",
    "Notes"
  ) %in% names(guide)))
  expect_true(nrow(guide) >= 10L)
  expect_true(any(grepl("precision_review_report", guide$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("build_summary_table_bundle", guide$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("facets_output_file_bundle", guide$MainFunction, fixed = TRUE)))
})

test_that("mfrmr_output_guide supports focused scopes", {
  reviews <- mfrmr_output_guide("reviews")
  expect_true(nrow(reviews) > 0L)
  expect_true(all(reviews$Scope == "reviews"))
  expect_true(all(reviews$OutputFamily == "review"))

  exports <- mfrmr_output_guide("exports")
  expect_true(nrow(exports) > 0L)
  expect_true(all(exports$Scope == "exports"))
  expect_true(any(grepl("export_summary_appendix", exports$MainFunction, fixed = TRUE)))

  gpcm <- mfrmr_output_guide("gpcm")
  expect_true(nrow(gpcm) > 0L)
  expect_false(any(gpcm$GPCMStatus == "supported"))
  expect_true(any(grepl("GPCM", gpcm$NextStep, ignore.case = TRUE)))
})

test_that("facets_feature_coverage separates implemented and unsupported FACETS surfaces", {
  coverage <- facets_feature_coverage()

  expect_s3_class(coverage, "data.frame")
  expect_true(all(c(
    "FACETSArea",
    "FACETSFeature",
    "FACETSReference",
    "mfrmrRoute",
    "Status",
    "Scope",
    "GapOrBoundary",
    "Priority"
  ) %in% names(coverage)))
  expect_true(nrow(coverage) >= 40L)
  expect_true(any(grepl("Table 14", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented"))
  expect_true(any(grepl("Wright map", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented"))
  expect_true(any(grepl("Generalizability Theory", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented" &
                    grepl("mfrm_d_study", coverage$mfrmrRoute, fixed = TRUE)))
  expect_true(any(grepl("Connectivity network graph", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented" &
                    grepl("mfrm_network_analysis", coverage$mfrmrRoute, fixed = TRUE) &
                    grepl("type = \"network\"", coverage$mfrmrRoute, fixed = TRUE)))
  expect_true(any(grepl("Residuals output file", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented" &
                    grepl("write_mfrm_residual_file", coverage$mfrmrRoute, fixed = TRUE)))
  expect_true(any(grepl("Category information function", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented" &
                    grepl("category_information", coverage$mfrmrRoute, fixed = TRUE)))
  expect_true(any(grepl("Cumulative probability curves", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "implemented" &
                    grepl("type = \"cumulative\"", coverage$mfrmrRoute, fixed = TRUE)))
  expect_true(any(grepl("Winsteps", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "not_implemented"))
  expect_true(any(grepl("command-file parser", coverage$FACETSFeature, fixed = TRUE) &
                    coverage$Status == "not_targeted"))

  partial <- facets_feature_coverage("partial")
  expect_true(nrow(partial) > 0L)
  expect_true(all(partial$Status == "partial"))

  missing <- facets_feature_coverage("not_implemented")
  expect_true(nrow(missing) > 0L)
  expect_true(all(missing$Status == "not_implemented"))
})

test_that("facets_positioning_guide prevents FACETS numerical-clone wording", {
  guide <- facets_positioning_guide()

  expect_s3_class(guide, "data.frame")
  expect_true(all(c(
    "Topic",
    "Position",
    "RecommendedWording",
    "PrimaryRoute"
  ) %in% names(guide)))
  expect_true(any(guide$Topic == "Estimation authority"))
  expect_true(any(guide$Topic == "External FACETS comparison"))
  expect_true(any(grepl("package-native", guide$Position, fixed = TRUE)))
  expect_true(any(grepl("not evidence of FACETS numerical equivalence",
                        guide$RecommendedWording,
                        fixed = TRUE)))
  expect_true(any(grepl("read_facets_fit_table", guide$PrimaryRoute, fixed = TRUE)))
  expect_false(any(grepl("\\baudit\\b", unlist(guide), ignore.case = TRUE)))
})

test_that("mfrmr_output_guide gives FACETS, ConQuest, and R user pathways", {
  facets <- mfrmr_output_guide("facets")
  expect_true(nrow(facets) >= 10L)
  expect_true(all(facets$Scope == "facets"))
  expect_true(any(grepl("facets_positioning_guide", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("not a FACETS numerical clone", facets$Notes, fixed = TRUE)))
  expect_true(any(grepl("facets_feature_coverage", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("review_mfrm_anchors", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("make_anchor_table", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("group_anchors", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("detect_anchor_drift", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_anchor_drift", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("fit_measures_table", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("facets_fit_df_guide", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("df_sensitivity", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("df_sensitive", facets$NextStep, fixed = TRUE)))
  expect_true(any(grepl("rating_scale_table", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("fair_average_table", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("estimate_bias", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("bias_interaction_report", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("bias_pairwise_report", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_bias_interaction", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_wright_unified", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("type = \"wright\"", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("data_quality_report", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("write_mfrm_residual_file", facets$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("write_mfrm_subset_file", facets$MainFunction, fixed = TRUE)))

  conquest <- mfrmr_output_guide("conquest")
  expect_true(nrow(conquest) >= 3L)
  expect_true(all(conquest$Scope == "conquest"))
  expect_true(any(grepl("build_conquest_overlap_bundle", conquest$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("review_conquest_overlap", conquest$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("less free than ConQuest", conquest$Question, fixed = TRUE)))

  r_path <- mfrmr_output_guide("r")
  expect_true(nrow(r_path) >= 3L)
  expect_true(all(r_path$Scope == "r"))
  expect_true(any(grepl("plot_data", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_data_components", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_long", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("type = \"pathway\"", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_bias_interaction", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_information", r_path$MainFunction, fixed = TRUE)))
  expect_true(any(grepl("plot_data_components", r_path$NextStep, fixed = TRUE)))
  expect_true(any(grepl("pathway_long", r_path$NextStep, fixed = TRUE)))
  expect_true(any(grepl("annotations/settings", r_path$NextStep, fixed = TRUE)))
  expect_true(any(grepl("ggplot2", r_path$UseWhen, fixed = TRUE)))
})

test_that("plot_data extracts reusable payloads and selected components", {
  fit <- make_toy_fit(maxit = 8)
  dq <- data_quality_report(fit)

  dashboard_payload <- plot_data(dq, type = "dashboard")
  expect_true(is.list(dashboard_payload))
  expect_true(all(c(
    "quality_flags",
    "score_support",
    "facet_response_patterns"
  ) %in% names(dashboard_payload)))

  flags <- plot_data(dq, type = "dashboard", component = "quality_flags")
  expect_s3_class(flags, "data.frame")
  dashboard_components <- plot_data_components(dq, type = "dashboard")
  expect_s3_class(dashboard_components, "data.frame")
  expect_true(all(c(
    "PlotName", "Component", "Role", "ObjectType", "Rows",
    "Columns", "Accessor", "Notes"
  ) %in% names(dashboard_components)))
  expect_true(any(dashboard_components$Component == "quality_flags"))

  plotted <- plot(dq, type = "dashboard", draw = FALSE)
  expect_s3_class(plotted, "mfrm_plot_data")
  score_support <- plot_data(plotted, component = "score_support")
  expect_s3_class(score_support, "data.frame")
  plotted_components <- plot_data_components(plotted)
  expect_true(any(plotted_components$Component == "score_support"))
  expect_error(plot_data(plotted, component = "not_a_component"), "component")

  curves <- category_curves_report(fit, theta_points = 21)
  cumulative_payload <- plot_data(curves, type = "cumulative")
  expect_true(all(c(
    "cumulative_probabilities",
    "cumulative_boundaries",
    "cumulative_direction"
  ) %in% names(cumulative_payload)))
  cumulative_table <- plot_data(curves, type = "cumulative", component = "cumulative_probabilities")
  expect_s3_class(cumulative_table, "data.frame")
  category_probability_payload <- plot_data(curves, type = "category_probability")
  expect_identical(category_probability_payload$plot, "ccc")
  expect_identical(category_probability_payload$plot_settings$RequestedType[1], "category_probability")
  expect_true(any(category_probability_payload$plot_long$PlotType == "ccc" &
                    category_probability_payload$plot_long$DisplayedByDefault))
  expect_true(all(c("plot_annotations", "curve_summary") %in%
                    names(category_probability_payload)))
  category_probability_components <- plot_data_components(curves, type = "category_probability")
  expect_true(any(category_probability_components$Component == "plot_long" &
                    category_probability_components$Role == "primary_data"))
  expect_true(any(category_probability_components$Component == "plot_annotations" &
                    category_probability_components$Role == "annotation"))
  category_information <- plot_data(
    curves,
    type = "category_information",
    component = "category_information"
  )
  expect_s3_class(category_information, "data.frame")
  expect_true(all(c("CategoryInformation", "CategoryInformationShare") %in% names(category_information)))

  curve_long <- plot_data(curves, component = "plot_long")
  expect_s3_class(curve_long, "data.frame")
  expect_true(all(c(
    "PlotType", "Panel", "CurveGroup", "Theta", "Series",
    "ValueName", "Value", "DisplayedByDefault"
  ) %in% names(curve_long)))
  expect_true(all(c(
    "ogive", "ccc", "cumulative", "information", "category_information"
  ) %in% unique(curve_long$PlotType)))
  expect_true(any(curve_long$PlotType == "cumulative" &
                    curve_long$Direction == "at_or_below" &
                    curve_long$DisplayedByDefault))
  curve_style <- plot_data(curves, preset = "monochrome", component = "curve_style")
  expect_s3_class(curve_style, "data.frame")
  expect_true(all(c("Series", "Colour", "LineType", "Preset") %in% names(curve_style)))
  expect_true(all(curve_style$Preset == "monochrome"))
  expect_gt(length(unique(curve_style$LineType)), 1L)

  diag <- diagnose_mfrm(fit, residual_pca = "none")
  pathway_long <- plot_data(
    fit,
    type = "pathway",
    diagnostics = diag,
    component = "pathway_long"
  )
  expect_s3_class(pathway_long, "data.frame")
  expect_true(all(c("Layer", "CurveGroup", "Theta", "Value") %in% names(pathway_long)))
  expect_true(any(pathway_long$Layer == "expected_score"))

  pathway_fit <- plot_data(
    fit,
    type = "pathway",
    diagnostics = diag,
    component = "fit_measures"
  )
  expect_s3_class(pathway_fit, "data.frame")
  expect_true(all(c("Facet", "Level", "Infit", "Outfit", "FitStatus") %in% names(pathway_fit)))
  pathway_components <- plot_data_components(fit, type = "pathway", diagnostics = diag)
  expect_true(any(pathway_components$Component == "pathway_long" &
                    pathway_components$Role == "primary_data"))
  expect_true(any(pathway_components$Component == "fit_measures" &
                    pathway_components$Role == "fit_review"))

  pathway_without_fit <- plot_data(
    fit,
    type = "pathway",
    include_fit_measures = FALSE,
    component = "fit_measure_status"
  )
  expect_s3_class(pathway_without_fit, "data.frame")
  expect_false(isTRUE(pathway_without_fit$Available[1]))
  expect_identical(pathway_without_fit$Status[1], "not_requested")

  info <- compute_information(fit, theta_points = 21)
  sem_long <- plot_data(
    plot_information(info, type = "sem", draw = FALSE),
    component = "plot_long"
  )
  expect_s3_class(sem_long, "data.frame")
  expect_true(any(sem_long$ValueName == "ConditionalSEM"))
  sem_components <- plot_data_components(plot_information(info, type = "sem", draw = FALSE))
  expect_true(any(sem_components$Component == "conditional_sem"))

  bias <- suppressWarnings(suppressMessages(
    estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 1)
  ))
  bias_long <- plot_data(
    plot_bias_interaction(bias, plot = "heatmap", draw = FALSE),
    component = "plot_long"
  )
  expect_s3_class(bias_long, "data.frame")
  expect_true(any(bias_long$Layer == "heatmap_cell"))
  bias_components <- plot_data_components(
    plot_bias_interaction(bias, plot = "heatmap", draw = FALSE)
  )
  expect_true(any(bias_components$Component == "flag_summary" &
                    bias_components$Role == "summary_or_guidance"))
})
