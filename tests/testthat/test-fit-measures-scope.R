test_that("supplied diagnostics cannot bypass ordinary fit-measure restrictions", {
  diagnostics <- list(measures = data.frame(Facet = "Rater", Level = "r1",
    Estimate = .2, SE = .1, Infit = 1, Outfit = 1, N = 12L))
  corrected <- structure(list(config = list(jml_adjustment = TRUE)), class = "mfrm_fit")
  two_family <- structure(list(config = list(gpcm_spec =
    list(component_levels = list(Task = c("a", "b"), Rater = c("x", "y"))))),
    class = "mfrm_fit")
  testlet <- structure(list(), class = "mfrm_testlet")
  shared <- structure(list(), class = "mfrm_random_rater")
  local_mocked_bindings(diagnose_mfrm = function(...) stop("Unexpected diagnostics"),
    .package = "mfrmr")
  for (supplied in list(NULL, diagnostics)) for (df in c("engine", "facets", "both")) {
    expect_error(fit_measures_table(two_family, supplied, fit_df_method = df), "two slope families")
    expect_error(fit_measures_table(corrected, supplied, fit_df_method = df), "corrected JML")
    expect_error(fit_measures_table(testlet, supplied, fit_df_method = df), "testlet/shared-rater")
    expect_error(fit_measures_table(shared, supplied, fit_df_method = df), "testlet/shared-rater")
  }
  # Detached legacy tables remain a supported descriptive input; they do not
  # acquire an uncertainty qualification or require a fitted object.
  result <- fit_measures_table(diagnostics)
  expect_equal(result$table$Measure, .2)
  expect_equal(result$table$Infit, 1)
  expect_false(any(result$table$CIEligible | result$table$SupportsFormalInference))
  expect_identical(result$table$PrecisionTier, "unknown")
})

test_that("ggplot conversion keeps saved measure bands and individual levels", {
  skip_if_not_installed("ggplot2")
  table <- data.frame(Facet = c("Task", rep("Rater", 4)),
    Level = c("a", "a", "b", "c", "d"), Measure = c(.1, .2, .3, .4, NA),
    SE = 10, CI_Lower = c(-.1, -.3, .3, NA, -1),
    CI_Upper = c(.25, .7, .3, NA, 1),
    Fixed = c(FALSE, FALSE, TRUE, FALSE, FALSE),
    FitStatus = c("within_band", "underfit", "not_available", "mixed", "overfit"),
    CIEligible = FALSE, CIUse = "screening_only")
  saved <- mfrmr:::new_mfrm_plot_data("fit_measures", list(plot = "measure_ci",
    table = table, title = "Exploratory bands", caption = "Screening only.",
    preset = "standard", display = list(show_notes = FALSE)))
  p <- as_ggplot(saved)
  layers <- ggplot2::ggplot_build(p)$data
  expect_equal(layers[[2]]$x, c(-.1, -.3))
  expect_equal(layers[[2]]$xend, c(.25, .7))
  expect_equal(layers[[2]]$y, c(-1, -2))
  expect_equal(layers[[3]]$x, c(.1, .2, .3, .4))
  expect_equal(layers[[3]]$y, -(1:4))
  expect_equal(layers[[3]]$shape, c(21, 21, 5, 4))
  expect_identical(p$scales$get_scales("y")$labels,
    c("Task: a", "Rater: a", "Rater: b", "Rater: c"))
  expect_identical(plot_data(p, "table"), table)
  expect_identical(attr(p, "mfrmr_plot_data"), saved)
  expect_null(p$labels$title)
  expect_null(p$labels$caption)
  expect_equal(ggplot2::ggplot_build(as_ggplot(saved, component = "table"))$data, layers)
  saved$data$display$show_notes <- TRUE
  expect_identical(as_ggplot(saved)$labels$caption, "Screening only.")
  expect_error(as_ggplot(saved, component = "caption"), "plot_data")
  expect_error(as_ggplot(saved, ci_level = .9), "must be empty")
  # Older detached tables may have no Fixed column, or no available intervals.
  saved$data$table$Fixed <- NULL
  saved$data$table$CI_Lower <- saved$data$table$CI_Upper <- NA_real_
  empty <- ggplot2::ggplot_build(as_ggplot(saved))$data
  expect_equal(nrow(empty[[2]]), 0L)
  expect_equal(empty[[3]]$shape, rep(4, 4))
})
