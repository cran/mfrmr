# Smoke tests for the plot helpers in R/api-plotting-extras.R
# (plot_threshold_ladder, plot_person_fit, plot_rater_severity_profile,
# plot_apa_figure_one, plot_dif_summary). These pin the contract of the
# 0.1.6 visualization surface; the assertions focus on contract
# (class, slot names, reference-free draw) rather than exact pixel
# output so they remain stable across graphics-device variants.

local({
  toy <- load_mfrmr_data("example_core")
  .fit <<- suppressWarnings(suppressMessages(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 20)
  ))
  .diag <<- suppressMessages(
    diagnose_mfrm(.fit, residual_pca = "none",
                  diagnostic_mode = "legacy")
  )
  .toy <<- toy
})

# --- plot_threshold_ladder --------------------------------------------------

test_that("plot_threshold_ladder returns structured payload", {
  p <- plot_threshold_ladder(.fit, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(all(c("data", "n_disorder_groups", "title", "subtitle",
                    "legend", "reference_lines") %in% names(p$data)))
  expect_s3_class(p$data$data, "data.frame")
  expect_true(all(c("Group", "Step", "Threshold", "Disordered") %in%
                    names(p$data$data)))
  expect_gte(nrow(p$data$data), 1L)
})

test_that("plot_threshold_ladder draw = TRUE runs without error", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_silent(plot_threshold_ladder(.fit, draw = TRUE))
})

test_that("plot_threshold_ladder honours highlight_disorder = FALSE", {
  p <- plot_threshold_ladder(.fit, highlight_disorder = FALSE, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# --- plot_person_fit --------------------------------------------------------

test_that("plot_person_fit returns one row per Person", {
  p <- plot_person_fit(.fit, diagnostics = .diag, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(all(c("Person", "Infit", "Outfit", "N", "Status") %in%
                    names(p$data$data)))
  expect_true(all(c("plot_long", "person_fit_indices", "flag_summary",
                    "plot_settings") %in% names(p$data)))
  expect_true(all(c("ReportIndex", "ReportValue", "ReviewStatus",
                    "ReportCaveat") %in% names(p$data$data)))
  person_n <- length(unique(.fit$facets$person$Person))
  expect_lte(nrow(p$data$data), person_n)
  expect_true(all(p$data$data$Status %in%
                    c("in_band", "one_outside", "both_outside")))
  expect_true(any(p$data$plot_long$Metric == "ReportValue"))
})

test_that("plot_person_fit honours custom fit envelope", {
  p <- plot_person_fit(.fit, diagnostics = .diag,
                       lower = 0.7, upper = 1.3, draw = FALSE)
  expect_equal(p$data$lower, 0.7)
  expect_equal(p$data$upper, 1.3)
})

test_that("plot_person_fit draws without error", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_silent(plot_person_fit(.fit, diagnostics = .diag, draw = TRUE,
                                 top_n_label = 3))
})

test_that("plot_person_fit exposes log-likelihood report view and monochrome preset", {
  p <- plot_person_fit(.fit, diagnostics = .diag,
                       fit_index = "loglik", preset = "monochrome",
                       draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$data$fit_index, "loglik")
  expect_identical(p$data$preset, "monochrome")
  expect_true(any(p$data$plot_long$Metric == "ReportValue" &
                    p$data$plot_long$DisplayedByDefault))
  expect_true(all(p$data$data$ReviewStatus %in% c(
    "review_1pct", "review_5pct", "not_flagged", "not_available", NA
  )))
  expect_true(any(p$data$reference_lines$role == "person-fit z threshold"))
})

# --- plot_rater_severity_profile --------------------------------------------

test_that("plot_rater_severity_profile produces CI whiskers", {
  p <- plot_rater_severity_profile(.fit, diagnostics = .diag,
                                    facet = "Rater", ci_level = 0.95,
                                    draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(all(c("Level", "Estimate", "SE", "CI_Lower",
                    "CI_Upper", "Band") %in% names(p$data$data)))
  expect_equal(p$data$ci_level, 0.95)
  # CI_Lower <= Estimate <= CI_Upper when SE is finite.
  valid <- is.finite(p$data$data$SE)
  if (any(valid)) {
    expect_true(all(p$data$data$CI_Lower[valid] -
                      p$data$data$Estimate[valid] <= 1e-6))
    expect_true(all(p$data$data$Estimate[valid] -
                      p$data$data$CI_Upper[valid] <= 1e-6))
  }
})

test_that("plot_rater_severity_profile validates ci_level", {
  expect_error(
    plot_rater_severity_profile(.fit, ci_level = -0.5, draw = FALSE),
    "ci_level"
  )
  expect_error(
    plot_rater_severity_profile(.fit, ci_level = c(0.9, 0.95), draw = FALSE),
    "ci_level"
  )
})

test_that("plot_rater_severity_profile errors on unknown facet", {
  expect_error(
    plot_rater_severity_profile(.fit, facet = "NoSuchFacet", draw = FALSE),
    "NoSuchFacet"
  )
})

# --- plot_dif_summary -------------------------------------------------------

test_that("plot_dif_summary accepts mfrm_dff output", {
  dff <- suppressWarnings(suppressMessages(
    analyze_dff(.fit, diagnostics = .diag,
                facet = "Rater", group = "Group",
                data = .toy, method = "residual")
  ))
  p <- plot_dif_summary(dff, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(all(c("Pair", "Effect", "SE", "Classification", "Color") %in%
                    names(p$data$data)))
  expect_gte(nrow(p$data$data), 1L)
})

test_that("plot_dif_summary rejects non-DIF inputs", {
  expect_error(plot_dif_summary(list(dif_table = data.frame()), draw = FALSE),
               "analyze_dff|analyze_dif")
})

test_that("plot_dif_summary validates top_n", {
  dff <- suppressWarnings(suppressMessages(
    analyze_dff(.fit, diagnostics = .diag,
                facet = "Rater", group = "Group",
                data = .toy, method = "residual")
  ))
  expect_error(plot_dif_summary(dff, top_n = 0, draw = FALSE), "`top_n`")
  expect_error(plot_dif_summary(dff, top_n = 1.5, draw = FALSE), "`top_n`")
  expect_error(plot_dif_summary(dff, top_n = Inf, draw = FALSE), "`top_n`")
})

test_that("plot_dif_summary supports CI and threshold guides", {
  dff <- suppressWarnings(suppressMessages(
    analyze_dff(.fit, diagnostics = .diag,
                facet = "Rater", group = "Group",
                data = .toy, method = "residual")
  ))
  p <- plot_dif_summary(dff, draw = FALSE, ci_level = 0.90,
                        effect_thresholds = c(screen = 0.5),
                        effect_axis_label = "Screening contrast")

  expect_true(all(c("CI_Lower", "CI_Upper", "ClassificationSystem") %in%
                    names(p$data$data)))
  expect_equal(p$data$settings$ci_level, 0.90)
  expect_equal(unname(p$data$settings$effect_thresholds), 0.5)
  expect_true(any(p$data$reference_lines$role == "threshold"))
  expect_true(is.data.frame(p$data$interpretation_guide))

  expect_error(plot_dif_summary(dff, draw = FALSE, ci_level = 1), "`ci_level`")
  expect_error(plot_dif_summary(dff, draw = FALSE,
                                effect_thresholds = c(-0.1)),
               "`effect_thresholds`")
  expect_error(plot_dif_summary(dff, draw = FALSE, effect_axis_label = ""),
               "`effect_axis_label`")
})

test_that("plot_dif_summary sort_by = 'effect' orders by signed contrast", {
  dff <- suppressWarnings(suppressMessages(
    analyze_dff(.fit, diagnostics = .diag,
                facet = "Rater", group = "Group",
                data = .toy, method = "residual")
  ))
  p <- plot_dif_summary(dff, sort_by = "effect", draw = FALSE)
  effects <- p$data$data$Effect
  expect_true(!is.unsorted(effects))
})

# --- plot_apa_figure_one ----------------------------------------------------

test_that("plot_apa_figure_one bundles the four panels", {
  expect_warning(
    p <- plot_apa_figure_one(.fit, diagnostics = .diag, draw = FALSE),
    NA
  )
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(all(c("wright", "severity", "threshold", "summary") %in%
                    names(p$data$data)))
  expect_s3_class(p$data$data$wright, "mfrm_plot_data")
  expect_s3_class(p$data$data$severity, "mfrm_plot_data")
  expect_s3_class(p$data$data$threshold, "mfrm_plot_data")
  expect_true(length(p$data$data$summary) >= 1L)
  expect_true(any(grepl("N obs = 768", p$data$data$summary, fixed = TRUE)))
  expect_true(any(grepl("Persons = 48", p$data$data$summary, fixed = TRUE)))
})

test_that("plot_apa_figure_one draws the 2x2 composite", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  # Allow warnings on composite draw (layout() may warn on null
  # devices); we only require the call to succeed without error.
  expect_no_error(suppressWarnings(
    plot_apa_figure_one(.fit, diagnostics = .diag, draw = TRUE)
  ))
})

# --- plot(fit, type = "ccc_overlay") ---------------------------------------

test_that("plot(fit, type = 'ccc_overlay') carries overlay rows", {
  p <- plot(.fit, type = "ccc_overlay", draw = FALSE)
  expect_identical(p$name, "category_characteristic_curves_overlay")
  expect_true("overlay" %in% names(p$data))
  expect_s3_class(p$data$overlay, "data.frame")
  expect_true(all(c("Bin", "Theta", "Category", "Proportion", "N") %in%
                    names(p$data$overlay)))
  expect_gte(nrow(p$data$overlay), 1L)
})

# --- plot(fit, type = "wright", group = ...) -------------------------------

test_that("plot(fit, type = 'wright') subgroup overlay via group_data", {
  p <- plot(.fit, type = "wright", group = "Group",
            group_data = .toy, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_s3_class(p$data$group, "data.frame")
  expect_true(all(c("Group", "Theta", "Density") %in%
                    names(p$data$group)))
  expect_gte(length(unique(p$data$group$Group)), 2L)
})

test_that("plot(fit, type = 'wright') without group has no group payload", {
  p <- plot(.fit, type = "wright", draw = FALSE)
  expect_null(p$data$group)
})

# --- plot_bias_interaction(plot = 'heatmap') -------------------------------

test_that("plot_bias_interaction heatmap mode returns bundle", {
  bias <- suppressWarnings(suppressMessages(
    estimate_bias(.fit, .diag, facet_a = "Rater",
                  facet_b = "Criterion", max_iter = 1)
  ))
  p <- plot_bias_interaction(bias, plot = "heatmap", draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$name, "table13_bias")
  expect_identical(p$data$plot, "heatmap")
  expect_s3_class(p$data$heatmap_cells, "data.frame")
  expect_true(is.matrix(p$data$heatmap_matrix))
  expect_true(is.matrix(p$data$heatmap_flag_matrix))
  expect_equal(dim(p$data$heatmap_matrix), dim(p$data$heatmap_flag_matrix))
  expect_true(all(c("Facet1_Level", "Facet2_Level", "BiasSize",
                    "CollapsedRows", "RepresentativeRule") %in%
                    names(p$data$heatmap_cells)))
  expect_true(all(c(
    "plot_long", "plot_annotations", "flag_summary", "plot_settings"
  ) %in% names(p$data)))
  expect_true(all(c(
    "Layer", "PlotType", "X", "Y", "Value", "ValueName", "Flag"
  ) %in% names(p$data$plot_long)))
  expect_true(any(p$data$plot_long$Layer == "heatmap_cell"))
  expect_true(all(c("AnnotationType", "Axis", "Value", "Label") %in%
                    names(p$data$plot_annotations)))
  expect_true(all(c("Trigger", "Rows", "Percent") %in% names(p$data$flag_summary)))
  expect_true(all(c("Plot", "AbsTWarn", "AbsBiasWarn", "CIRequested", "Preset") %in%
                    names(p$data$plot_settings)))
  expect_true(any(grepl("Table 13", p$data$reference_notes$Reference,
                         fixed = TRUE)))

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_bias_interaction(bias, plot = "heatmap", draw = TRUE))
})

# --- plot_anchor_drift(type = 'forest') -----------------------------------

test_that("plot_anchor_drift forest draws from mfrm_anchor_drift", {
  fit_a <- suppressWarnings(suppressMessages(
    fit_mfrm(load_mfrmr_data("example_bias"), "Person",
             c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5)
  ))
  fit_b <- suppressWarnings(suppressMessages(
    fit_mfrm(load_mfrmr_data("example_bias"), "Person",
             c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5)
  ))
  drift <- suppressWarnings(suppressMessages(
    detect_anchor_drift(list(W1 = fit_a, W2 = fit_b))
  ))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  p <- plot_anchor_drift(drift, type = "forest",
                          ci_level = 0.95, draw = TRUE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(all(c("data", "ci_level") %in% names(p$data)))
  expect_gte(nrow(p$data$data), 1L)
})

# --- plot.mfrm_equating_chain -----------------------------------------------

test_that("plot.mfrm_equating_chain common_anchors falls back cleanly", {
  fit_a <- suppressWarnings(suppressMessages(
    fit_mfrm(load_mfrmr_data("example_bias"), "Person",
             c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5)
  ))
  fit_b <- suppressWarnings(suppressMessages(
    fit_mfrm(load_mfrmr_data("example_bias"), "Person",
             c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5)
  ))
  chain <- suppressWarnings(suppressMessages(
    build_equating_chain(list(WaveA = fit_a, WaveB = fit_b))
  ))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  p <- plot(chain, type = "common_anchors", draw = TRUE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$name, "equating_chain_common_anchors")
})
