test_that("analyze_facet_equivalence returns the expected bundle structure", {
  toy <- load_mfrmr_data("example_core")
  fit <- fit_mfrm(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  )

  eq <- analyze_facet_equivalence(fit, facet = "Rater")

  expect_s3_class(eq, "mfrm_facet_equivalence")
  expect_true(all(c("summary", "chi_square", "pairwise", "rope", "forest", "settings") %in% names(eq)))
  expect_equal(nrow(eq$summary), 1)
  expect_true(all(c("FixedChiSq", "FixedProb", "Separation", "Reliability") %in% names(eq$chi_square)))
  expect_true(all(c("P_TOST", "Equivalent") %in% names(eq$pairwise)))
  expect_true(all(c("ROPEPct", "ROPEStatus") %in% names(eq$rope)))
  expect_true(eq$summary$Facet[1] == "Rater")
  expect_true(all(c("AllPairsEquivalent", "AnyPairEquivalent", "PairwiseDecisionBasis") %in% names(eq$summary)))
  expect_true(eq$summary$Decision[1] %in% c(
    "all_pairs_equivalent",
    "partial_pairwise_equivalence",
    "no_pairwise_equivalence_established"
  ))
  expect_true(eq$summary$BF01[1] >= 0)
})

test_that("plot_facet_equivalence prepares forest and rope views without drawing", {
  toy <- load_mfrmr_data("example_core")
  fit <- fit_mfrm(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  )
  eq <- analyze_facet_equivalence(fit, facet = "Rater")

  forest_obj <- plot_facet_equivalence(eq, draw = FALSE)
  rope_obj <- plot(eq, type = "rope", draw = FALSE)

  expect_type(forest_obj, "list")
  expect_equal(forest_obj$type, "forest")
  expect_equal(rope_obj$type, "rope")
  expect_true(is.data.frame(forest_obj$data))
  expect_true(is.data.frame(rope_obj$data))
})
