test_that("supplied diagnostics cannot bypass reliability model restrictions", {
  diagnostics <- list(reliability = data.frame(Facet = "Rater",
    Reliability = .8, Separation = 2, Strata = 3))
  corrected <- structure(list(config = list(jml_adjustment = TRUE)), class = "mfrm_fit")
  two_family <- structure(list(config = list(gpcm_spec =
    list(component_levels = list(Rater = c("a", "b"), Task = c("x", "y"))))),
    class = "mfrm_fit")
  local_mocked_bindings(diagnose_mfrm = function(...) stop("Unexpected diagnostics"), .package = "mfrmr")
  for (metric in c("reliability", "separation", "strata")) {
    for (supplied in list(NULL, diagnostics)) {
      expect_error(plot_reliability_snapshot(corrected, supplied, metric, draw = FALSE), "corrected JML")
      expect_error(plot_reliability_snapshot(two_family, supplied, metric, draw = FALSE), "two slope families")
    }
  }
})

test_that("reliability snapshots retain row-aligned precision interpretation", {
  fit <- structure(list(config = list(method = "JML")), class = "mfrm_fit")
  rel <- data.frame(Facet = c("Rater", "Person", "Criterion"),
    Reliability = c(.4, .8, NA), Separation = c(1, 2, NA), Strata = c(2, 3, NA),
    PrecisionTier = "exploratory", ReliabilityUse = "screening_only",
    ReliabilityBasis = "Approximate marginal-SE variance subtraction; covariance not propagated",
    SupportsFormalInference = FALSE,
    ModelSummaryNote = c("Rater scope", "Person scope", "Unavailable"))
  precision <- c("PrecisionTier", "ReliabilityUse", "ReliabilityBasis",
    "SupportsFormalInference", "ModelSummaryNote")
  local_mocked_bindings(diagnose_mfrm = function(...) stop("Unexpected diagnostics"), .package = "mfrmr")
  for (metric in c("reliability", "separation", "strata")) {
    p <- plot_reliability_snapshot(fit, list(reliability = rel), metric, draw = FALSE)
    rows <- match(p$data$table$Facet, rel$Facet)
    expect_identical(p$data$table$Facet, c("Person", "Rater"))
    expect_equal(p$data$table[precision], rel[rows, precision], ignore_attr = TRUE)
    expect_match(p$data$notes$Text, "Ordinary JML.*screening summaries")
    expect_identical(unserialize(serialize(p, NULL))$data$table, p$data$table)
  }
  # Previously saved minimal diagnostic tables remain usable without inventing flags.
  p <- plot_reliability_snapshot(fit, list(reliability = rel[c("Facet", "Reliability")]), draw = FALSE)
  expect_identical(names(p$data$table), c("Facet", "Metric", "Value"))
  expect_identical(p$data$table$Value, c(.8, .4))
})
