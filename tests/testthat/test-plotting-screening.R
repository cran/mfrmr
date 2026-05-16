# Smoke tests for the 4 screening / case-level plot helpers in
# R/api-plotting-screening.R (plot_guttman_scalogram, plot_residual_qq,
# plot_rater_trajectory, plot_rater_agreement_heatmap). Assertions
# focus on the mfrm_plot_data payload contract; pixel-exact output
# is not part of the contract.

# Shared toy fixtures. `make_toy_fit()` and `make_toy_diagnostics()`
# live in helper-fixtures.R and cache their result so that re-running
# the file does not pay the fitting cost more than once.
.toy <- load_mfrmr_data("example_core")
.fit <- make_toy_fit(maxit = 15)
.diag <- make_toy_diagnostics(.fit)

# --- plot_guttman_scalogram ------------------------------------------------

test_that("plot_guttman_scalogram returns matrix + unexpected overlay", {
  p <- plot_guttman_scalogram(.fit, diagnostics = .diag, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(is.matrix(p$data$matrix))
  expect_gte(nrow(p$data$matrix), 2L)
  expect_gte(ncol(p$data$matrix), 2L)
  expect_s3_class(p$data$unexpected, "data.frame")
})

test_that("plot_guttman_scalogram accepts a custom column_facet", {
  p <- plot_guttman_scalogram(.fit, diagnostics = .diag,
                               column_facet = "Rater", draw = FALSE)
  expect_identical(p$data$column_facet, "Rater")
})

test_that("plot_guttman_scalogram caps rows at top_n_persons", {
  p <- plot_guttman_scalogram(.fit, diagnostics = .diag,
                               top_n_persons = 10L, draw = FALSE)
  expect_lte(nrow(p$data$matrix), 10L)
})

test_that("plot_guttman_scalogram draws without error", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(
    plot_guttman_scalogram(.fit, diagnostics = .diag, draw = TRUE)
  ))
})

# --- plot_residual_qq ------------------------------------------------------

test_that("plot_residual_qq returns theoretical vs sample quantiles", {
  p <- plot_residual_qq(.fit, diagnostics = .diag, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(all(c("Person", "Theoretical", "Sample") %in%
                    names(p$data$data)))
  expect_gte(nrow(p$data$data), 2L)
  # Theoretical quantiles are monotone.
  expect_true(!is.unsorted(p$data$data$Theoretical))
})

test_that("plot_residual_qq draws without error", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(
    plot_residual_qq(.fit, diagnostics = .diag, draw = TRUE)
  ))
})

# --- plot_rater_trajectory -------------------------------------------------

test_that("plot_rater_trajectory assembles per-wave rows", {
  p <- plot_rater_trajectory(list(T1 = .fit, T2 = .fit),
                              ci_level = 0.95, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(all(c("Wave", "Level", "Estimate", "SE",
                    "CI_Lower", "CI_Upper") %in%
                    names(p$data$data)))
  expect_gte(nrow(p$data$data), 2L)
  expect_setequal(unique(p$data$data$Wave), c("T1", "T2"))
})

test_that("plot_rater_trajectory requires at least two fits", {
  expect_error(plot_rater_trajectory(list(T1 = .fit), draw = FALSE),
               "at least two")
})

test_that("plot_rater_trajectory requires named list", {
  expect_error(plot_rater_trajectory(list(.fit, .fit), draw = FALSE),
               "distinct non-empty names")
})

test_that("plot_rater_trajectory draws without error", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(
    plot_rater_trajectory(list(T1 = .fit, T2 = .fit), draw = TRUE)
  ))
})

# --- plot_rater_agreement_heatmap ------------------------------------------

test_that("plot_rater_agreement_heatmap returns symmetric matrix", {
  p <- plot_rater_agreement_heatmap(.fit, diagnostics = .diag,
                                     draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(is.matrix(p$data$matrix))
  expect_equal(dim(p$data$matrix)[1], dim(p$data$matrix)[2])
  # Symmetric off-diagonal.
  m <- p$data$matrix
  valid <- is.finite(m)
  expect_true(all(m[valid] == t(m)[valid]))
  # Diagonal is 1 (identity agreement).
  expect_true(all(diag(m) == 1))
})

test_that("plot_rater_agreement_heatmap draws without error", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(
    plot_rater_agreement_heatmap(.fit, diagnostics = .diag, draw = TRUE)
  ))
})

test_that("plot_rater_agreement_heatmap metric='correlation' maps to Corr", {
  p <- plot_rater_agreement_heatmap(.fit, diagnostics = .diag,
                                     metric = "correlation", draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$data$metric, "correlation")
  expect_true(grepl("Corr", p$data$subtitle))
  m <- p$data$matrix
  finite_off <- is.finite(m) & row(m) != col(m)
  # correlations stay within [-1, 1]
  expect_true(all(m[finite_off] >= -1 & m[finite_off] <= 1))
  # symmetric off-diagonal
  expect_true(all(m[finite_off] == t(m)[finite_off]))
})

test_that("plot_rater_agreement_heatmap rejects unknown metrics", {
  expect_error(
    plot_rater_agreement_heatmap(.fit, diagnostics = .diag,
                                  metric = "kappa", draw = FALSE),
    "'arg' should be one of"
  )
})
