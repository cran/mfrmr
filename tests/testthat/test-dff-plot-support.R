dff_plot_support_fixture <- function(method = "residual") {
  pairs <- data.frame(Level = paste0("R", 1:4), Group1 = "A", Group2 = "B",
    Contrast = c(.1, -.4, NA_real_, 0), SE = NA_real_, Method = method,
    N_Group1 = c(20L, 10L, 5L, 20L), N_Group2 = c(30L, 20L, 0L, 20L),
    NPersons_Group1 = c(10L, 5L, 3L, 10L), NPersons_Group2 = c(15L, 10L, 0L, 10L),
    sparse = c(FALSE, FALSE, TRUE, FALSE), Classification = "Residual contrast",
    ClassificationSystem = "descriptive", FormalInferenceEligible = FALSE,
    ScaleLinkStatus = "not_applicable")
  cells <- data.frame(Level = rep(paste0("R", 1:4), each = 2),
    GroupValue = rep(c("A", "B"), 4), N = c(20L, 30L, 10L, 20L, 5L, 0L, 20L, 20L),
    NPersons = c(10L, 15L, 5L, 10L, 3L, 0L, 10L, 10L),
    ObsExpAvg = c(.2, .1, -.3, .1, .2, NA_real_, .1, .1),
    StdResidual = c(1, .5, -1.5, .5, NA_real_, NA_real_, .5, .5),
    sparse = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))
  cells <- cells[cells$N > 0, ]
  cases <- data.frame(InputRow = 1:3, FitRow = c(1L, 2L, NA_integer_),
    GroupValue = c("A", "B", NA), Included = c(TRUE, TRUE, FALSE),
    Reason = c("included", "included", "missing_group"))
  if (method == "refit") {
    pairs$ScaleLinkStatus <- c("linked", "weak_link", "linked", "linked")
    pairs$Fixed1 <- c(FALSE, FALSE, TRUE, FALSE)
    pairs$Fixed2 <- pairs$Fixed1
    pairs$Classification[3] <- "Fixed target; no comparison"
    pairs$SE[1] <- .2
    cells <- cases <- NULL
  }
  structure(list(dif_table = pairs, cell_table = cells, cases = cases,
    config = list(facet = "Rater", group = "Group", method = method,
      min_obs = 10L, residual_output_version = 3L,
      refit_output_version = 1L)), class = c("mfrm_dff", "mfrm_dif"))
}

test_that("DFF summaries retain full comparisons and distinguish display omissions", {
  x <- dff_plot_support_fixture()
  original <- serialize(x, NULL)
  local_mocked_bindings(fit_mfrm = function(...) stop("No fitting"),
    analyze_dff = function(...) stop("No new analysis"), .package = "mfrmr")
  p <- plot(x, top_n = 2, draw = FALSE)
  pairs <- plot_data(p, "pairs")
  expect_identical(pairs[names(x$dif_table)], x$dif_table)
  expect_identical(pairs$PlotStatus,
    c("Displayed", "Displayed", "Contrast unavailable", "Not selected by top_n"))
  expect_equal(p$data$data$Effect, c(-.4, .1))
  expect_identical(plot_data(p, "cases"), x$cases)
  expect_identical(plot_data(p, "analysis_config"), x$config)
  expect_match(paste(plot_data(p, "notes")$Text, collapse = " "),
    "2 of 4 source comparisons displayed; 1 unavailable and 1 not selected")
  expect_match(paste(capture.output(print(p)), collapse = " "), "Group1 minus Group2")
  all <- plot(x, draw = FALSE)
  expect_identical(all$data$pairs$PlotStatus[4], "Displayed")
  expect_equal(tail(all$data$data$Effect, 1), 0)
  expect_true(all(is.na(all$data$data$CI_Lower)))
  expect_error(plot(x, ci_level = .95, draw = FALSE), "Confidence intervals are unavailable")
  expect_identical(serialize(x, NULL), original)
  x$dif_table$Contrast[] <- NA_real_
  expect_error(plot(x, draw = FALSE), "No finite contrast values")
})

test_that("DFF heatmaps retain sparse cells and missing contrasts on their original scales", {
  x <- dff_plot_support_fixture()
  for (metric in c("obs_exp", "t", "contrast")) {
    p <- plot(x, type = "heatmap", metric = metric, draw = FALSE)
    expect_identical(plot_data(p, "cases"), x$cases)
    expect_identical(plot_data(p, "analysis_config"), x$config)
    if (metric == "contrast") {
      expect_identical(p$data$pairs[names(x$dif_table)], x$dif_table)
      expect_true(is.na(p$data$matrix["R3", "A vs B"]))
      expect_equal(p$data$matrix["R4", "A vs B"], 0)
    } else {
      expect_identical(p$data$pairs, x$cell_table)
      expect_true(is.na(p$data$matrix["R3", "B"]))
      if (metric == "obs_exp") expect_equal(p$data$matrix["R3", "A"], .2)
      else expect_true(is.na(p$data$matrix["R3", "A"]))
    }
  }
})

test_that("refit plot support preserves weak links and fixed targets without qualifying inference", {
  x <- dff_plot_support_fixture("refit")
  for (type in c("summary", "heatmap")) {
    p <- if (type == "summary") plot(x, ci_level = .95, draw = FALSE) else
      plot(x, type = type, metric = "contrast", draw = FALSE)
    expect_identical(p$data$pairs[names(x$dif_table)], x$dif_table)
    expect_null(plot_data(p, "cases"))
    notes <- paste(plot_data(p, "notes")$Text, collapse = " ")
    expect_match(notes, "Group1 minus Group2 in facet locations")
    expect_match(notes, "estimated anchors")
    expect_match(notes, "larger Group1 location")
    if (type == "summary") {
      expect_identical(p$data$settings$effect_axis_label, "Subgroup location difference (logits)")
      expect_identical(p$data$pairs$PlotStatus[3], "Contrast unavailable")
      expect_true(is.na(p$data$data$CI_Lower[1])) # Weak link has a point, no interval.
      expect_equal(p$data$data$CI_Lower[2], .1 - qnorm(.975) * .2)
    } else {
      expect_match(p$data$subtitle, "formal DFF inference remains unavailable")
    }
  }
})

test_that("saved DFF support survives ggplot conversion without recalculation", {
  skip_if_not_installed("ggplot2")
  local_mocked_bindings(fit_mfrm = function(...) stop("No fitting"),
    analyze_dff = function(...) stop("No new analysis"), .package = "mfrmr")
  for (type in c("summary", "heatmap")) {
    x <- dff_plot_support_fixture()
    p <- if (type == "summary") plot(x, top_n = 2, draw = FALSE) else
      plot(x, type = type, metric = "contrast", draw = FALSE)
    saved <- unserialize(serialize(p, NULL))
    g <- as_ggplot(saved)
    ggplot2::ggplot_build(g)
    expect_identical(plot_data(g), plot_data(p))
    expect_identical(plot_data(g, "notes"), plot_data(p, "notes"))
    expect_null(g$labels$title)
    expect_null(g$labels$caption)
  }
})

test_that("an entirely unavailable contrast heatmap stays missing without range warnings", {
  x <- dff_plot_support_fixture()
  x$dif_table$Contrast[] <- NA_real_
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_warning(p <- plot(x, type = "heatmap", metric = "contrast"), NA)
  expect_true(all(is.na(p$data$matrix)))
  expect_identical(p$data$pairs[names(x$dif_table)], x$dif_table)
  skip_if_not_installed("ggplot2")
  expect_warning(ggplot2::ggplot_build(as_ggplot(p)), NA)
})

test_that("DFF comparison labels fit the device without changing the plotted results", {
  x <- dff_plot_support_fixture()
  x$dif_table$Level <- c("Organization", "Language", "Accuracy", "Content")
  grDevices::pdf(file = NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  axis <- graphics::axis
  left_margin <- label_width <- NULL
  local_mocked_bindings(axis = function(side, ...) {
    if (side == 2) {
      args <- list(...)
      left_margin <<- graphics::par("mai")[2]
      label_width <<- max(graphics::strwidth(args$labels,
        units = "inches", cex = args$cex.axis))
    }
    axis(side, ...)
  }, .package = "mfrmr")
  p <- plot(x)
  expect_gt(left_margin, label_width)
  expect_identical(plot_data(p), plot_data(plot(x, draw = FALSE)))
  # Labels that cannot leave a usable plotting region must not silently clip.
  x$dif_table$Level[1] <- paste(rep("long comparison label", 20), collapse = " ")
  expect_error(plot(x), "Comparison labels need more horizontal space")
  expect_s3_class(plot(x, draw = FALSE), "mfrm_plot_data")
})
