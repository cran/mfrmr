# Smoke tests for the secondary plot helpers in
# R/api-plotting-secondary.R: plot_local_dependence_heatmap,
# plot_reliability_snapshot, plot_residual_matrix, plot_shrinkage_funnel.
# Assertions focus on the mfrm_plot_data payload contract; pixel-exact
# output is not part of the contract.

# Shared toy fixtures. `make_toy_fit()` and `make_toy_diagnostics()`
# live in helper-fixtures.R and cache their result so that re-running
# the file does not pay the fitting cost more than once.
.toy <- load_mfrmr_data("example_core")
.fit <- make_toy_fit()
.diag <- make_toy_diagnostics(.fit)

# --- plot_local_dependence_heatmap ----------------------------------------

test_that("plot_local_dependence_heatmap builds a symmetric matrix", {
  p <- plot_local_dependence_heatmap(.fit, diagnostics = .diag, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(is.matrix(p$data$matrix))
  expect_equal(dim(p$data$matrix)[1], dim(p$data$matrix)[2])
  expect_true(all(diag(p$data$matrix) == 1))
})

test_that("plot_local_dependence_heatmap rejects unknown facet", {
  expect_error(
    plot_local_dependence_heatmap(.fit, diagnostics = .diag,
                                   facet = "NotARealFacet", draw = FALSE),
    "must be one of"
  )
})

test_that("plot_local_dependence_heatmap draws without error", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(
    plot_local_dependence_heatmap(.fit, diagnostics = .diag, draw = TRUE)
  ))
})

# --- plot_reliability_snapshot --------------------------------------------

test_that("plot_reliability_snapshot returns a tidy table", {
  p <- plot_reliability_snapshot(.fit, diagnostics = .diag, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(is.data.frame(p$data$table))
  expect_true(all(c("Facet", "Metric", "Value") %in% names(p$data$table)))
})

test_that("plot_reliability_snapshot supports separation and strata metrics", {
  for (m in c("separation", "strata")) {
    p <- plot_reliability_snapshot(.fit, diagnostics = .diag,
                                     metric = m, draw = FALSE)
    expect_identical(p$data$metric, m)
  }
})

# --- plot_residual_matrix -------------------------------------------------

test_that("plot_residual_matrix caps rows at top_n_persons", {
  p <- plot_residual_matrix(.fit, diagnostics = .diag,
                              top_n_persons = 6L, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_lte(nrow(p$data$matrix), 6L)
})

test_that("plot_residual_matrix rejects unknown facet", {
  expect_error(
    plot_residual_matrix(.fit, diagnostics = .diag,
                          facet = "NotARealFacet", draw = FALSE),
    "must be one of"
  )
})

# --- plot_shrinkage_funnel -----------------------------------------------

test_that("plot_shrinkage_funnel works on an EB-augmented fit", {
  fit_eb <- suppressMessages(apply_empirical_bayes_shrinkage(.fit))
  p <- plot_shrinkage_funnel(fit_eb, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(is.data.frame(p$data$table))
  expect_true(all(c("RawEstimate", "ShrunkEstimate", "ShrinkageFactor")
                  %in% names(p$data$table)))
})

test_that("plot_shrinkage_funnel rejects fits without EB columns", {
  expect_error(
    plot_shrinkage_funnel(.fit, draw = FALSE),
    "empirical-Bayes shrinkage"
  )
})
