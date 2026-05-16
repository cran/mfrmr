# Tests covering the 0.1.6 CI API consistency surface:
#   * lifecycle::deprecate_warn() path for analyze_facet_equivalence's
#     conf_level alias,
#   * show_ci / ci_level on plot_fair_average / plot_displacement /
#     plot_bias_interaction,
#   * CI column contract (CI_Lower / CI_Upper / CI_Level) populated
#     only when show_ci = TRUE.

local({
  .toy <<- load_mfrmr_data("example_core")
  .fit <<- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15)
  ))
  .diag <<- suppressMessages(
    diagnose_mfrm(.fit, residual_pca = "none",
                  diagnostic_mode = "legacy")
  )
})

# --- analyze_facet_equivalence deprecation ---------------------------------

test_that("conf_level emits a lifecycle deprecation warning", {
  old_opt <- options(lifecycle_verbosity = "warning")
  on.exit(options(old_opt), add = TRUE)
  w <- testthat::capture_warnings(
    analyze_facet_equivalence(.fit, diagnostics = .diag,
                              facet = "Rater", conf_level = 0.9)
  )
  expect_true(any(grepl("conf_level", paste(w, collapse = " "))))
})

test_that("ci_level is silent and becomes the active bound", {
  old_opt <- options(lifecycle_verbosity = "warning")
  on.exit(options(old_opt), add = TRUE)
  expect_no_warning(
    eq <- analyze_facet_equivalence(.fit, diagnostics = .diag,
                                    facet = "Rater", ci_level = 0.9)
  )
  expect_equal(unname(eq$settings$ci_level), 0.9)
})

test_that("supplying both routes honors conf_level and warns", {
  old_opt <- options(lifecycle_verbosity = "warning")
  on.exit(options(old_opt), add = TRUE)
  w <- testthat::capture_warnings(
    eq <- analyze_facet_equivalence(.fit, diagnostics = .diag,
                                    facet = "Rater",
                                    ci_level = 0.95, conf_level = 0.8)
  )
  expect_true(any(grepl("conf_level", paste(w, collapse = " "))))
  expect_equal(unname(eq$settings$ci_level), 0.8)
})

# --- plot_fair_average show_ci ---------------------------------------------

test_that("plot_fair_average(show_ci = FALSE) omits CI columns (default)", {
  p <- plot_fair_average(.fit, plot_type = "scatter", draw = FALSE)
  expect_false("CI_Lower" %in% names(p$data$data))
})

test_that("plot_fair_average(show_ci = TRUE) adds clipped delta-method CI", {
  p <- plot_fair_average(.fit, plot_type = "scatter",
                         show_ci = TRUE, ci_level = 0.95, draw = FALSE)
  expect_true(all(c("CI_Lower", "CI_Upper", "CI_Level") %in%
                    names(p$data$data)))
  expect_true(all(unique(p$data$data$CI_Level) == 0.95))
  # Every finite CI should respect the rating range.
  r_min <- .fit$prep$rating_min
  r_max <- .fit$prep$rating_max
  valid <- is.finite(p$data$data$CI_Lower) &
    is.finite(p$data$data$CI_Upper)
  if (any(valid)) {
    expect_true(all(p$data$data$CI_Lower[valid] >= r_min - 1e-8))
    expect_true(all(p$data$data$CI_Upper[valid] <= r_max + 1e-8))
  }
})

test_that("plot_fair_average validates ci_level", {
  expect_error(
    plot_fair_average(.fit, show_ci = TRUE, ci_level = 1.1, draw = FALSE),
    "ci_level"
  )
})

# --- plot_displacement show_ci ---------------------------------------------

test_that("plot_displacement(show_ci = TRUE) adds CI columns", {
  p <- plot_displacement(.fit, diagnostics = .diag,
                         show_ci = TRUE, ci_level = 0.9, draw = FALSE)
  expect_true(all(c("CI_Lower", "CI_Upper", "CI_Level") %in%
                    names(p$data$table)))
  valid <- is.finite(p$data$table$CI_Lower) &
    is.finite(p$data$table$CI_Upper)
  expect_true(any(valid))
  if (any(valid)) {
    # Each finite CI should be a proper interval (upper >= lower) and
    # be centred on the point estimate (|Disp - mid| < 1e-6).
    widths <- p$data$table$CI_Upper[valid] - p$data$table$CI_Lower[valid]
    expect_true(all(widths >= -1e-8))
    mid <- (p$data$table$CI_Lower[valid] +
              p$data$table$CI_Upper[valid]) / 2
    expect_true(all(abs(p$data$table$Displacement[valid] - mid) <= 1e-6))
  }
})

test_that("plot_displacement default omits CI columns", {
  p <- plot_displacement(.fit, diagnostics = .diag, draw = FALSE)
  expect_false("CI_Lower" %in% names(p$data$table))
})

# --- plot_bias_interaction show_ci ----------------------------------------

test_that("plot_bias_interaction(show_ci) populates CI on ranked + scatter", {
  bias <- suppressMessages(suppressWarnings(
    estimate_bias(.fit, .diag, facet_a = "Rater",
                  facet_b = "Criterion", max_iter = 1)
  ))
  p_r <- plot_bias_interaction(bias, plot = "ranked",
                               show_ci = TRUE, ci_level = 0.95,
                               draw = FALSE)
  p_s <- plot_bias_interaction(bias, plot = "scatter",
                               show_ci = TRUE, ci_level = 0.95,
                               draw = FALSE)
  expect_true(all(c("CI_Lower", "CI_Upper", "CI_Level") %in%
                    names(p_r$data$ranked_table)))
  expect_true(all(c("CI_Lower", "CI_Upper", "CI_Level") %in%
                    names(p_s$data$scatter_data)))
})

test_that("plot_bias_interaction CI omitted when show_ci = FALSE", {
  bias <- suppressMessages(suppressWarnings(
    estimate_bias(.fit, .diag, facet_a = "Rater",
                  facet_b = "Criterion", max_iter = 1)
  ))
  p <- plot_bias_interaction(bias, plot = "ranked", draw = FALSE)
  expect_false("CI_Lower" %in% names(p$data$ranked_table))
})

# --- plot(fit, type = "wright", show_ci = TRUE) ---------------------------

test_that("plot(fit, type = 'wright', show_ci = TRUE) builds SE table", {
  p <- plot(.fit, type = "wright", show_ci = TRUE, ci_level = 0.95,
            draw = FALSE)
  # Wright map overlays SE via build_wright_map_data; the returned
  # payload surfaces `CI_Lower` / `CI_Upper` on the locations table.
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(any(grepl("CI", names(p$data$locations))))
})
