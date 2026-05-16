# Tests for hierarchical structure / small-sample review helpers (0.1.6).
# Cross-refs:
#  * `detect_facet_nesting()` - entropy-based classification
#  * `facet_small_sample_review()` - Linacre (1994) bands
#  * `compute_facet_icc()` - lme4 variance components (Suggests)
#  * `compute_facet_design_effect()` - Kish (1965) design effect
#  * `analyze_hierarchical_structure()` - bundle report
#  * `fit$summary$FacetSampleSizeFlag` - post-fit summary flag

skip_on_cran_check <- function() {
  if (!is.null(getOption("mfrmr.skip_nested_tests"))) skip("disabled via option")
}

test_that("detect_facet_nesting runs on example_core", {
  toy <- load_mfrmr_data("example_core")
  nested <- detect_facet_nesting(toy, c("Rater", "Criterion"))
  expect_s3_class(nested, "mfrm_facet_nesting")
  expect_true(nrow(nested$pairwise_table) >= 1L)
  expect_true(all(c("NestingIndex_AinB", "NestingIndex_BinA", "Direction") %in%
                    names(nested$pairwise_table)))
  expect_true(nested$summary$NFacets >= 2L)
})

test_that("detect_facet_nesting detects perfect nesting", {
  # Construct Rater nested in Region: each rater appears in exactly one region.
  d <- data.frame(
    Person = rep(sprintf("P%02d", 1:20), each = 6),
    Rater = rep(sprintf("R%d", 1:6), 20),
    Region = rep(rep(c("A", "A", "B", "B", "C", "C"), 20)),
    Score = sample(0:4, 120, replace = TRUE),
    stringsAsFactors = FALSE
  )
  nested <- detect_facet_nesting(d, c("Rater", "Region"))
  rater_in_region <- nested$pairwise_table[
    nested$pairwise_table$FacetA == "Rater" &
      nested$pairwise_table$FacetB == "Region", ]
  expect_gte(rater_in_region$NestingIndex_AinB, 0.99)
  expect_equal(rater_in_region$ClassificationAinB, "Fully nested")
})

test_that("detect_facet_nesting flags crossed designs as Crossed", {
  toy <- load_mfrmr_data("example_core")
  nested <- detect_facet_nesting(toy, c("Rater", "Criterion"))
  pairs <- nested$pairwise_table
  # example_core is a fully-crossed rater x criterion design.
  expect_true(any(pairs$Direction == "crossed"))
})

test_that("facet_small_sample_review classifies Linacre bands correctly", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 20)
  ))
  audit <- facet_small_sample_review(fit)
  expect_s3_class(audit, "mfrm_facet_sample_review")
  expect_true(all(c("Facet", "Level", "N", "SampleCategory") %in% names(audit$table)))
  expect_true(all(audit$table$SampleCategory %in%
                    c("sparse", "marginal", "standard", "strong", NA_character_)))
})

test_that("fit$summary carries FacetSampleSizeFlag", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 20)
  ))
  expect_true("FacetSampleSizeFlag" %in% names(fit$summary))
  expect_true(as.character(fit$summary$FacetSampleSizeFlag[1]) %in%
                c("sparse", "marginal", "standard", "strong"))
  expect_true(is.finite(suppressWarnings(
    as.integer(fit$summary$FacetMinLevelN[1])
  )))
})

test_that("analyze_hierarchical_structure returns a bundle", {
  skip_on_cran_check()
  toy <- load_mfrmr_data("example_core")
  hs <- suppressMessages(analyze_hierarchical_structure(
    toy, facets = c("Rater", "Criterion"),
    person = "Person", score = "Score",
    compute_icc = FALSE, igraph_layout = FALSE
  ))
  expect_s3_class(hs, "mfrm_hierarchical_structure")
  expect_true("nesting" %in% names(hs))
  expect_true("crosstabs" %in% names(hs))
  expect_true(length(hs$crosstabs) >= 1L)
})

test_that("compute_facet_icc works when lme4 is available", {
  skip_if_not_installed("lme4")
  toy <- load_mfrmr_data("example_core")
  icc <- suppressMessages(suppressWarnings(
    compute_facet_icc(toy, facets = c("Rater", "Criterion"),
                      score = "Score", person = "Person")
  ))
  expect_s3_class(icc, "mfrm_facet_icc")
  expect_true(all(c("Facet", "Variance", "ICC", "Interpretation") %in% names(icc)))
  expect_true(all(is.na(icc$ICC) | (icc$ICC >= 0 & icc$ICC <= 1)))
})

test_that("compute_facet_design_effect computes Kish deff", {
  skip_if_not_installed("lme4")
  toy <- load_mfrmr_data("example_core")
  icc <- suppressMessages(suppressWarnings(
    compute_facet_icc(toy, facets = c("Rater", "Criterion"),
                      score = "Score", person = "Person")
  ))
  deff <- compute_facet_design_effect(toy, facets = c("Rater", "Criterion"),
                                      icc_table = icc)
  expect_s3_class(deff, "mfrm_facet_design_effect")
  expect_true(all(c("Facet", "AvgClusterSize", "ICC", "DesignEffect", "EffectiveN") %in%
                    names(deff)))
  # Deff should be >= 1 whenever ICC and avg cluster are finite and non-negative.
  finite_rows <- is.finite(deff$DesignEffect)
  if (any(finite_rows)) {
    expect_true(all(deff$DesignEffect[finite_rows] >= 1 - 1e-8 |
                      !is.finite(deff$ICC[finite_rows]) |
                      deff$ICC[finite_rows] <= 0))
  }
})

test_that("reporting_checklist surfaces hierarchical audit rows", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 20)
  ))
  diag <- suppressMessages(suppressWarnings(
    diagnose_mfrm(fit, residual_pca = "none")
  ))
  chk <- reporting_checklist(fit, diagnostics = diag)
  items <- chk$checklist$Item
  expect_true("Facet sample-size adequacy" %in% items)
  expect_true("Hierarchical structure review" %in% items)
})

test_that("build_mfrm_manifest carries hierarchical_review table", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 20)
  ))
  manifest <- build_mfrm_manifest(fit)
  expect_true("hierarchical_review" %in% names(manifest))
  expect_true("FacetSampleSizeFlag" %in% names(manifest$hierarchical_review))
})
