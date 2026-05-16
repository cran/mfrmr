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
})

test_that("visual_reporting_template covers second-wave plots", {
  tbl <- visual_reporting_template()
  for (fam in c("Guttman scalogram", "Residual Q-Q",
                "Rater trajectory (linked waves)",
                "Rater agreement heatmap")) {
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
