local({
  d <- mfrmr:::sample_mfrm_data(seed = 2026)

  .marg_fit <<- suppressWarnings(
    fit_mfrm(
      d,
      "Person",
      c("Rater", "Task", "Criterion"),
      "Score",
      method = "MML",
      model = "RSM",
      quad_points = 7,
      maxit = 30
    )
  )
  .marg_diag <<- diagnose_mfrm(.marg_fit, residual_pca = "none", diagnostic_mode = "both")
  .marg_diag_legacy <<- diagnose_mfrm(.marg_fit, residual_pca = "none", diagnostic_mode = "legacy")
})

with_null_device_local <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  testthat::expect_gt(grDevices::dev.cur(), 1)
  value <- force(expr)
  invisible(value)
}

test_that("plot_marginal_fit returns reusable plot payload", {
  p <- plot_marginal_fit(.marg_diag, draw = FALSE, preset = "publication")
  p_prop <- plot_marginal_fit(.marg_fit, diagnostics = .marg_diag, plot_type = "prop_diff", draw = FALSE)

  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$name, "marginal_fit")
  expect_identical(as.character(p$data$preset), "publication")
  expect_true(all(c(
    "plot", "table", "full_table", "summary", "facet_summary",
    "step_summary", "guidance", "thresholds",
    "title", "subtitle", "legend", "reference_lines"
  ) %in% names(p$data)))
  expect_gt(nrow(p$data$table), 0)
  expect_s3_class(p_prop, "mfrm_plot_data")
  expect_identical(as.character(p_prop$data$plot), "prop_diff")
})

test_that("plot_marginal_pairwise returns reusable plot payload", {
  p <- plot_marginal_pairwise(.marg_diag, draw = FALSE, preset = "publication")
  p_adj <- plot_marginal_pairwise(
    .marg_fit,
    diagnostics = .marg_diag,
    metric = "adjacent",
    facet = "Rater",
    draw = FALSE
  )

  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$name, "marginal_pairwise")
  expect_identical(as.character(p$data$preset), "publication")
  expect_true(all(c(
    "plot", "table", "full_table", "summary", "pair_stats",
    "guidance", "thresholds", "title", "subtitle",
    "legend", "reference_lines"
  ) %in% names(p$data)))
  expect_gt(nrow(p$data$table), 0)
  expect_s3_class(p_adj, "mfrm_plot_data")
  expect_identical(as.character(p_adj$data$plot), "adjacent")
  expect_true(all(as.character(p_adj$data$table$Facet) == "Rater"))
})

test_that("plot_marginal helpers require strict marginal diagnostics", {
  expect_error(
    plot_marginal_fit(.marg_diag_legacy, draw = FALSE),
    "Strict marginal diagnostics are not available",
    fixed = TRUE
  )
  expect_error(
    plot_marginal_pairwise(.marg_diag_legacy, draw = FALSE),
    "Strict marginal diagnostics are not available",
    fixed = TRUE
  )
})

test_that("plot_marginal helpers draw on a null graphics device", {
  p_fit <- with_null_device_local(
    plot_marginal_fit(.marg_diag, draw = TRUE)
  )
  p_pair <- with_null_device_local(
    plot_marginal_pairwise(.marg_diag, metric = "adjacent", draw = TRUE)
  )

  expect_s3_class(p_fit, "mfrm_plot_data")
  expect_s3_class(p_pair, "mfrm_plot_data")
})
