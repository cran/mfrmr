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
    "plot(timing",
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
  expect_true(all(guide$DefaultLevel[!is.na(guide$DefaultLevel)] == 0.95))
  expect_true(any(grepl("show_ci = TRUE", guide$PrimaryHelper, fixed = TRUE)))
  expect_true(any(grepl("ci_level = 0.95", guide$PrimaryHelper, fixed = TRUE)))
  expect_true(any(grepl("delta-method", guide$Basis, fixed = TRUE)))
  expect_true(any(grepl("profile", guide$Basis, ignore.case = TRUE)))
  expect_true(any(grepl("not global model-fit proof",
                        guide$InterpretationBoundary,
                        fixed = TRUE)))
  guide_text <- paste(unlist(guide, use.names = FALSE), collapse = "\n")
  expect_true(grepl("plot(..., type = \"shrinkage\", view = \"funnel\", show_ci = TRUE, ci_level = 0.95)",
                    guide_text, fixed = TRUE))
  expect_true(grepl("plot(fit, type = \"shrinkage\", show_ci = TRUE, ci_level = 0.95)",
                    guide_text, fixed = TRUE))
  ci_public_helpers <- c(
    "analyze_facet_equivalence",
    "analyze_hierarchical_structure",
    "compute_facet_icc",
    "fair_average_table",
    "fit_measures_table",
    "plot(drift, type = \"forest\"",
    "type = \"apa_summary\"",
    "type = \"bias\"",
    "type = \"summary\"",
    "type = \"displacement\"",
    "type = \"fair_average\"",
    "view = \"funnel\"",
    "type = \"rater_severity\"",
    "type = \"trajectory\"",
    "view = \"unified\""
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
  expect_true(any(grepl("GPCM", gpcm$Notes, fixed = TRUE)))
  expect_true(any(gpcm$Route == "Group contrast summary" &
                    gpcm$GPCMStatus == "supported_with_caveat"))
  expect_true(any(gpcm$Route == "Anchor drift forest plot" &
                    grepl("linking synthesis supported_with_caveat",
                          gpcm$GPCMStatus, fixed = TRUE)))
  expect_true(any(grepl("screening evidence", gpcm$InterpretationBoundary,
                        fixed = TRUE)))
})

test_that("the interval guide separates GPCM targets and uncertainty-only output", {
  guide <- mfrmr_interval_guide("gpcm")
  routes <- c("Fixed-facet model/sandwich intervals", "GPCM relative-slope intervals",
    "GMFRM component-slope intervals", "GPCM slope profile intervals",
    "Corrected-JML local root uncertainty")
  expect_true(all(routes %in% guide$Route))
  expect_false(anyDuplicated(guide$Route) > 0L)
  selected <- guide[match(routes, guide$Route), ]
  expect_match(selected$GPCMStatus[1], "supported_with_caveat one-family MML locations")
  expect_match(selected$GPCMStatus[1], "experimental two-family")
  expect_match(selected$GPCMStatus[2], "one-family GPCM MML only")
  expect_match(selected$IntervalColumns[3], "CoverageQualified")
  expect_match(selected$PrimaryHelper[4], 'method = "profile"', fixed = TRUE)
  expect_match(selected$InterpretationBoundary[4], "unresolved side")
  expect_true(is.na(selected$DefaultLevel[5]))
  expect_match(selected$IntervalColumns[5], "no structural Lower or Upper")
  expect_false(routes[5] %in% mfrmr_interval_guide("visual")$Route)
})
