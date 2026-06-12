test_that("visual_reporting_template returns beginner-oriented figure guidance", {
  tbl <- visual_reporting_template()

  expect_s3_class(tbl, "data.frame")
  expect_true(all(c(
    "FigureFamily", "Scope", "PrimaryHelper", "DefaultPlacement",
    "WhatToReport", "CaptionSkeleton", "ResultsWording", "WhatNotToClaim",
    "BeginnerCheck", "ThreeDPolicy"
  ) %in% names(tbl)))
  expect_true(all(c(
    "manuscript", "appendix", "diagnostic", "surface"
  ) %in% tbl$Scope))
  expect_true("Wright map" %in% tbl$FigureFamily)
  expect_true("Category probability surface" %in% tbl$FigureFamily)
  expect_match(
    tbl$ThreeDPolicy[tbl$FigureFamily == "Category probability surface"][1],
    "advanced surface data only",
    fixed = TRUE
  )
  expect_match(
    tbl$BeginnerCheck[tbl$FigureFamily == "Category probability surface"][1],
    "category_support",
    fixed = TRUE
  )
  expect_match(
    tbl$WhatNotToClaim[tbl$FigureFamily == "Residual PCA"][1],
    "standalone dimensionality test",
    fixed = TRUE
  )
  expect_match(
    tbl$CaptionSkeleton[tbl$FigureFamily == "Wright map"][1],
    "Figure X. Wright map",
    fixed = TRUE
  )
  expect_match(
    tbl$ResultsWording[tbl$FigureFamily == "Category probability surface"][1],
    "exploratory support",
    fixed = TRUE
  )
  expect_true("Empirical-Bayes shrinkage funnel" %in% tbl$FigureFamily)
  expect_match(
    tbl$PrimaryHelper[tbl$FigureFamily == "Empirical-Bayes shrinkage funnel"][1],
    "show_ci = TRUE",
    fixed = TRUE
  )
  expect_match(
    tbl$WhatNotToClaim[tbl$FigureFamily == "Empirical-Bayes shrinkage funnel"][1],
    "automatic evidence of rater quality",
    fixed = TRUE
  )
  expect_true("Response-time review" %in% tbl$FigureFamily)
  expect_match(
    tbl$PrimaryHelper[tbl$FigureFamily == "Response-time review"][1],
    "plot_response_time_review",
    fixed = TRUE
  )
  expect_match(
    tbl$WhatNotToClaim[tbl$FigureFamily == "Response-time review"][1],
    "speed-accuracy parameters",
    fixed = TRUE
  )
})

test_that("visual_reporting_template covers second-wave plots", {
  tbl <- visual_reporting_template()
  for (fam in c("Guttman scalogram", "Residual Q-Q",
                "Rater trajectory (linked waves)",
                "Rater agreement heatmap",
                "Response-time review",
                "Empirical-Bayes shrinkage funnel")) {
    expect_true(fam %in% tbl$FigureFamily,
                info = paste("missing row for", fam))
  }
  expect_match(
    tbl$WhatNotToClaim[tbl$FigureFamily == "Rater trajectory (linked waves)"][1],
    "anchor-linking",
    fixed = TRUE
  )
})

test_that("visual_reporting_template filters by reporting scope", {
  full_tbl <- visual_reporting_template()
  manuscript_tbl <- visual_reporting_template("manuscript")
  surface_tbl <- visual_reporting_template("surface")

  expect_gt(nrow(manuscript_tbl), 0)
  expect_gt(nrow(surface_tbl), 0)
  expect_true(all(manuscript_tbl$Scope == "manuscript"))
  expect_true(all(surface_tbl$Scope == "surface"))
  expect_equal(nrow(manuscript_tbl), sum(full_tbl$Scope == "manuscript"))
  expect_identical(surface_tbl$FigureFamily, "Category probability surface")
})

test_that("mfrmr_interval_guide maps CI-capable routes and filters scopes", {
  guide <- mfrmr_interval_guide()

  expect_s3_class(guide, "data.frame")
  expect_true(all(c(
    "Route", "Scope", "PrimaryHelper", "DisplayRoute", "DefaultLevel",
    "IntervalColumns", "Basis", "UseFor", "InterpretationBoundary",
    "GPCMStatus", "Notes"
  ) %in% names(guide)))
  expect_true(nrow(guide) >= 15L)
  expect_true(all(guide$DefaultLevel == 0.95))
  expect_true(any(grepl("show_ci = TRUE", guide$PrimaryHelper, fixed = TRUE)))
  expect_true(any(grepl("ci_level = 0.95", guide$PrimaryHelper, fixed = TRUE)))
  expect_true(any(grepl("delta-method", guide$Basis, fixed = TRUE)))
  expect_true(any(grepl("profile", guide$Basis, ignore.case = TRUE)))
  expect_true(any(grepl("not global model-fit proof",
                        guide$InterpretationBoundary,
                        fixed = TRUE)))
  guide_text <- paste(unlist(guide, use.names = FALSE), collapse = "\n")
  expect_true(grepl("plot_shrinkage_funnel(..., show_ci = TRUE, ci_level = 0.95)",
                    guide_text, fixed = TRUE))
  expect_true(grepl("plot(fit, type = \"shrinkage\", show_ci = TRUE, ci_level = 0.95)",
                    guide_text, fixed = TRUE))
  ci_public_helpers <- c(
    "analyze_facet_equivalence",
    "analyze_hierarchical_structure",
    "compute_facet_icc",
    "fair_average_table",
    "fit_measures_table",
    "plot_anchor_drift",
    "plot_apa_figure_one",
    "plot_bias_interaction",
    "plot_dif_summary",
    "plot_displacement",
    "plot_fair_average",
    "plot_shrinkage_funnel",
    "plot_rater_severity_profile",
    "plot_rater_trajectory",
    "plot_wright_unified"
  )
  for (helper in ci_public_helpers) {
    expect_true(grepl(helper, guide_text, fixed = TRUE), info = helper)
  }

  visual <- mfrmr_interval_guide("visual")
  expect_true(nrow(visual) > 0L)
  expect_true(all(grepl("(^|,)visual(,|$)", visual$Scope)))
  expect_true(any(visual$Route == "Wright map uncertainty overlay"))
  expect_true(any(visual$Route == "Bias-interaction interval overlay"))
  expect_true(any(visual$Route == "Rater severity profile"))
  expect_true(any(visual$Route == "Manuscript Figure 1 composite"))

  gpcm <- mfrmr_interval_guide("gpcm")
  expect_true(nrow(gpcm) > 0L)
  expect_true(all(grepl("(^|,)gpcm(,|$)", gpcm$Scope)))
  expect_true(any(grepl("bounded GPCM", gpcm$Notes, fixed = TRUE)))
  expect_true(any(gpcm$Route == "DFF / DIF contrast summary" &
                    gpcm$GPCMStatus == "supported_with_caveat"))
  expect_true(any(gpcm$Route == "Anchor drift forest plot" &
                    grepl("linking synthesis supported_with_caveat",
                          gpcm$GPCMStatus, fixed = TRUE)))
  expect_true(any(grepl("screening evidence", gpcm$InterpretationBoundary,
                        fixed = TRUE)))
})
