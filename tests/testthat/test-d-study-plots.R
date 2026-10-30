test_that("D-study curves keep every condition fixed and surfaces refuse mixed conditions", {
  # Stored main-effects components isolate plot semantics without fitting.
  gt <- structure(list(
    variance_components = data.frame(Source = c("Person", "Rater", "Task", "Occasion", "Residual"),
      Variance = c(1, .3, .2, .1, .8)),
    design = list(object_facet = "Person", random_facets = c("Rater", "Task", "Occasion"),
      calculation_version = 2L, identification_status = "identified")
  ), class = c("mfrm_generalizability", "list"))
  ds <- mfrm_d_study(gt, expand.grid(Rater = c(2, 4), Task = c(3, 6), Occasion = c(1, 2)))
  for (extra_group in list(NULL, "n_Task", "Metric")) {
    p <- plot_data(plot(ds, group_var = extra_group, draw = FALSE))
    curves <- split(p$series, p$series$Series)
    expect_length(curves, 8L)
    expect_true(all(vapply(curves, function(z) nrow(unique(z[c("n_Task", "n_Occasion", "Metric")])) == 1L, logical(1))))
    expect_true(all(grepl("Occasion =", names(curves), fixed = TRUE)))
  }
  expect_error(plot(ds, type = "heatmap", draw = FALSE), "Varying: n_Occasion")
  expect_error(plot(ds, type = "contour", draw = FALSE), "panel_by")
  expect_error(plot(ds, x_var = "G", draw = FALSE), "must be one of")
  expect_error(plot(ds, type = "heatmap", y_var = "Scenario", draw = FALSE), "must be one of")
  surface <- plot_data(plot(ds, type = "heatmap", panel_by = "n_Occasion", draw = FALSE))
  expect_equal(anyDuplicated(surface$surface[c("Panel", "X", "Y")]), 0L)
  expect_equal(nrow(surface$surface), nrow(ds))
  expect_identical(unique(surface$surface$Panel), c("Occasion = 1", "Occasion = 2"))
  expect_equal(nrow(surface$legend), 5L)
  expect_true(all(surface$legend$aesthetic == "fill"))
  sensitive <- mfrm_d_study(gt, expand.grid(Rater = c(2, 4), Task = c(3, 6), Occasion = c(1, 2)),
    residual_scaling = "sensitivity")
  expect_error(plot(sensitive, type = "heatmap", panel_by = "n_Occasion", draw = FALSE), "Varying: ResidualScaling")
  s <- plot_data(plot(sensitive, type = "heatmap", panel_grid = c("n_Occasion", "ResidualScaling"), draw = FALSE))
  expect_equal(anyDuplicated(s$surface[c("Panel", "X", "Y")]), 0L)
  expect_equal(length(unique(s$surface$Panel)), 6L)
  one_panel <- plot_data(plot(ds, panel_grid = "n_Occasion", draw = FALSE))
  expect_identical(one_panel$panel_by, "n_Occasion")
  expect_identical(unique(one_panel$series$Panel), c("Occasion = 1", "Occasion = 2"))
  expect_no_error(plot(ds[ds$n_Occasion == 1, ], type = "heatmap", draw = FALSE))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_error(as_ggplot(ds), "Use plot(x) to preserve", fixed = TRUE)
  }

  # Keep unavailable points in their original series; do not join across them.
  gap <- mfrm_d_study(gt, data.frame(Rater = c(2, 3, 4), Task = 3, Occasion = 1))
  gap$G[2] <- NA_real_
  p <- plot_data(plot(gap, draw = FALSE))
  expect_equal(nrow(p$series), 6L)
  expect_true(is.na(p$series$Value[p$series$Metric == "G" & p$series$X == 3]))
  empty <- gap
  empty$G <- empty$Phi <- NA_real_
  expect_error(plot(empty, draw = FALSE), "Inspect the D-study table")
  grDevices::pdf(NULL, width = 8, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(c("mfrow", "mar", "oma", "cex", "mex"))
  expect_no_warning(plot(gap, panel_by = "Metric"))
  expect_no_warning(plot(ds, type = "heatmap", panel_by = "n_Occasion"))
  expect_equal(graphics::par(names(old)), old)
  partial <- ds
  partial$Phi[partial$n_Occasion == 1] <- NA_real_
  labels <- character()
  original_text <- graphics::text
  local_mocked_bindings(text = function(...) {
    args <- list(...)
    if (length(args) >= 3L) labels <<- c(labels, as.character(args[[3L]]))
    original_text(...)
  }, .package = "mfrmr")
  expect_no_warning(plot(partial, type = "contour", panel_by = "n_Occasion"))
  expect_true(any(grepl("Estimates unavailable", labels, fixed = TRUE)))
  expect_false(any(grepl("constant surface", labels, fixed = TRUE)))

  titles <- character()
  original_title <- graphics::title
  local_mocked_bindings(title = function(main = NULL, ...) {
    if (!is.null(main)) titles <<- c(titles, main)
    original_title(main = main, ...)
  }, .package = "mfrmr")
  local_graphics_headings(function(main) titles <<- c(titles, main))
  expect_no_warning(plot(sensitive, panel_grid = c("Metric", "ResidualScaling")))
  expect_true(all(startsWith(titles[1:3], "G / ")))
  expect_true(all(startsWith(titles[4:6], "Phi / ")))
})


test_that("main-effects D-study guides are optional and never facet-count lines", {
  gt <- structure(list(
    variance_components = data.frame(Source = c("Person", "Rater", "Task", "Residual"),
      Variance = c(1, .3, .2, .8)),
    design = list(object_facet = "Person", random_facets = c("Rater", "Task"),
      calculation_version = 2L, identification_status = "identified")
  ), class = c("mfrm_generalizability", "list"))
  d <- mfrm_d_study(gt, expand.grid(Rater = c(2, 4), Task = c(3, 6)))
  expect_equal(nrow(plot_data(plot(d, draw = FALSE))$reference_lines), 0L)
  p <- plot_data(plot(d, metric = "Phi", reference = c(.7, .8), draw = FALSE))
  expect_equal(p$reference_lines$value, c(.7, .8))
  expect_match(paste(p$notes$Text, collapse = " "), "not universal acceptance rules")
  expect_error(plot(d, reference = Inf, draw = FALSE), "distinct finite numeric")
  expect_error(plot(d, reference = 1.1, draw = FALSE), "between 0 and 1")
  for (type in c("heatmap", "contour", "surface3d")) {
    expect_error(plot(d, type = type, reference = .8, draw = FALSE), "not surface plots")
  }
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(d, metric = "Phi", reference = c(.7, .8)))
  expect_no_warning(plot(d, type = "error_variance", reference = 2))
})
