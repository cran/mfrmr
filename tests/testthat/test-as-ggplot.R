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
  for (p in plots) expect_null(p$labels$caption)
  for (i in c(2L, 4L)) {
    subtitle <- plot_data(plots[[i]], "subtitle")
    expect_null(plots[[i]]$labels$subtitle)
    expect_match(gsub("\n", " ", subtitle),
      "reference profile fixes additive facet effects and fitted interactions at zero", fixed = TRUE)
  }
  expect_gte(length(ggplot2::ggplot_build(plots[[1]])$data), 8L)
})

test_that("fit curves share colour and line encodings across renderers", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  fit <- make_toy_fit(model = "PCM", maxit = 20)
  grDevices::pdf(NULL, width = 9, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  drawn <- list()
  original_lines <- graphics::lines
  testthat::local_mocked_bindings(lines = function(...) {
    args <- list(...)
    if (!is.null(args$lty)) drawn[[length(drawn) + 1L]] <<- args
    original_lines(...)
  }, .package = "mfrmr")
  for (type in c("ccc", "pathway")) {
    original <- .mfrmr_muffle_expected_warnings(
      plot(fit, type = type, draw = FALSE), "^Review-only display:")
    table_name <- if (type == "ccc") "probabilities" else "expected"
    key <- if (type == "ccc") "Category" else "CurveGroup"
    series <- unique(as.character(original$data[[table_name]][[key]]))
    for (preset in c("standard", "publication", "monochrome")) {
      drawn <- list()
      p <- .mfrmr_muffle_expected_warnings(plot(fit, type = type, preset = preset),
                                         "^Review-only display:")
      expected <- mfrmr:::.plot_series_colors(series, preset)
      lty <- mfrmr:::.plot_series_linetypes(series)
      expect_identical(p$data[[table_name]], original$data[[table_name]])
      expect_identical(p$data$fit_readiness, original$data$fit_readiness)
      colors <- vapply(drawn, function(a) unname(as.character(a$col)), character(1))
      types <- vapply(drawn, function(a) unname(as.character(a$lty)), character(1))
      expect_identical(unique(colors), unname(expected))
      expect_identical(unique(types), unique(unname(lty)))
      g <- ggplot2::ggplot_build(as_ggplot(p))$plot
      expect_identical(unname(g$scales$get_scales("colour")$map(series)), unname(expected))
      expect_identical(unname(g$scales$get_scales("linetype")$map(series)), unname(lty))
    }
    custom <- stats::setNames(rep("#333333", length(series)), series)
    p <- .mfrmr_muffle_expected_warnings(
      plot(fit, type = type, palette = custom, draw = FALSE), "^Review-only display:")
    expect_identical(p$data$palette, custom)
    g <- ggplot2::ggplot_build(as_ggplot(p))$plot
    expect_identical(unname(g$scales$get_scales("colour")$map(series)), unname(custom))
    expect_identical(p$data[[table_name]], original$data[[table_name]])
  }
  # Check opaque line ink against the two built-in light backgrounds.
  luminance <- function(col) {
    rgb <- grDevices::col2rgb(col) / 255
    linear <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
    colSums(linear * c(0.2126, 0.7152, 0.0722))
  }
  for (n in c(8, 10)) for (background in c("white", "#fcfdff")) {
    contrast <- (luminance(background) + 0.05) /
      (luminance(mfrmr:::.plot_series_colors(seq_len(n))) + 0.05)
    expect_true(all(contrast >= 3))
  }
  overlay <- .mfrmr_muffle_expected_warnings(
    plot(fit, type = "ccc_overlay", draw = FALSE), "^Review-only display:")
  for (slope in c("none", "linewidth", "alpha", "colour")) {
    expect_no_error(built <- ggplot2::ggplot_build(as_ggplot(overlay, slope_aes = slope)))
    points <- which(vapply(built$plot$layers, function(layer)
      inherits(layer$geom, "GeomPoint"), logical(1)))
    expect_gt(length(unique(built$data[[points[1L]]]$shape)), 1L)
  }
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


test_that("residual Q-Q payloads cannot become generic bar charts", {
  p <- new_mfrm_plot_data("residual_qq", list(
    table = data.frame(Person = c("A", "B"), Theoretical = c(-1, 1), Observed = c(-.5, 2))))
  expect_error(as_ggplot(p), "both quantile axes")
  expect_error(as_ggplot(p, component = "table"), "both quantile axes")
  expect_identical(plot_data(p, "table")$Observed, c(-.5, 2))
})
