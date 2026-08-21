test_that("as_ggplot converts the core fitted-result plots", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  fit <- make_toy_fit(maxit = 20)

  payloads <- .mfrmr_muffle_expected_warnings(
    list(
      plot(fit, type = "wright", draw = FALSE),
      plot(fit, type = "pathway", draw = FALSE),
      plot(fit, type = "fit_pathway", include_person = TRUE,
           top_n_person = 3, draw = FALSE),
      plot(fit, type = "ccc", draw = FALSE)
    ),
    "^Review-only display:"
  )
  plots <- lapply(payloads, as_ggplot)
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  for (p in plots) expect_no_error(ggplot2::ggplot_build(p))
  expect_gte(length(ggplot2::ggplot_build(plots[[1]])$data), 8L)
})

test_that("as_ggplot accepts fitted objects and CCC slope styling", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  fit <- make_toy_fit(maxit = 20)

  .mfrmr_muffle_expected_warnings({
    direct <- as_ggplot(fit, type = "fit_pathway", fit_stat = "Outfit")
    ccc <- as_ggplot(fit, type = "ccc", slope_aes = "linewidth")
  }, "^Review-only display:")
  expect_s3_class(direct, "ggplot")
  expect_s3_class(ccc, "ggplot")
  expect_no_error(ggplot2::ggplot_build(direct))
  expect_no_error(ccc_build <- ggplot2::ggplot_build(ccc))
  linewidths <- suppressWarnings(as.numeric(
    ccc_build$data[[1L]]$linewidth %||% numeric(0)
  ))
  linewidths <- linewidths[is.finite(linewidths)]
  expect_gt(length(linewidths), 0L)
  expect_gte(min(linewidths), 0.45)
  expect_lte(max(linewidths), 1.35)
})

test_that("as_ggplot converts DIF summary and heatmap payload shapes", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  summary_payload <- mfrmr:::new_mfrm_plot_data(
    "dif_summary",
    list(
      data = data.frame(
        Pair = c("A", "B"), Effect = c(-0.4, 0.7),
        CI_Lower = c(-0.6, 0.4), CI_Upper = c(-0.2, 1),
        Color = c("#2B6CB0", "#C53030")
      ),
      settings = list(effect_axis_label = "Contrast")
    )
  )
  heat_payload <- mfrmr:::new_mfrm_plot_data(
    "dif_heatmap",
    list(
      matrix = matrix(c(-1, 0.2, 0.4, 1.1), 2, 2,
                      dimnames = list(c("L1", "L2"), c("G1", "G2"))),
      flag_matrix = matrix(c(TRUE, FALSE, FALSE, TRUE), 2, 2),
      settings = list(show_values = TRUE, value_digits = 1L)
    )
  )
  expect_no_error(ggplot2::ggplot_build(as_ggplot(summary_payload)))
  expect_no_error(ggplot2::ggplot_build(as_ggplot(heat_payload)))
})

test_that("as_ggplot converts simulation handoff lists", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  payload <- list(
    plot = "design_evaluation",
    metric = "separation",
    metric_col = "MeanSeparation",
    x_var = "n_person",
    x_label = "Persons",
    data = data.frame(
      n_person = c(20, 40, 20, 40),
      y = c(1.1, 1.8, 1.3, 2.1),
      group = c("2 raters", "2 raters", "4 raters", "4 raters")
    )
  )
  class(payload) <- "mfrm_design_evaluation_plot_data"
  expect_no_error(ggplot2::ggplot_build(as_ggplot(payload)))

  signal <- payload
  signal$plot <- "signal_detection"
  signal$display_metric <- "Bias screening hit rate"
  class(signal) <- "mfrm_signal_detection_plot_data"
  expect_no_error(ggplot2::ggplot_build(as_ggplot(signal)))
})
