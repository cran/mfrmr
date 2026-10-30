extended_view_fixture <- function() {
  structure(list(table = data.frame(Person = c("P1", "P2", "Prior", "Missing", "P3"),
    Estimate = c(.5, -.5, 0, NA, .5), Lower = c(-.2, -1, -1.96, NA, -Inf),
    Upper = c(1.2, .2, 1.96, NA, Inf),
    Status = c("available_conditional", "available_conditional", "prior_only", "unavailable", "available_conditional")),
    settings = list(level = .95)), class = "mfrm_random_rater_scores")
}

test_that("views retain source rows and explicitly describe exclusions", {
  x <- extended_view_fixture(); before <- dev.cur()
  interval <- plot(x, draw = FALSE, sort = "estimate")$data
  expect_identical(interval$table, x$table)
  expect_identical(interval$display$rows, c(2L, 3L, 1L, 5L, 4L))
  precision <- plot(x, draw = FALSE, style = "precision")$data
  expect_equal(precision$display_data$Width[1:3], c(1.4, 1.2, 3.92))
  expect_identical(precision$display_data$Included, c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_identical(precision$display_data$Shape[3], 1)
  expect_identical(precision$display_data$Reason[5], "Finite interval unavailable")
  distribution <- plot(x, draw = FALSE, style = "distribution")$data
  expect_equal(distribution$distribution, data.frame(Estimate = c(-.5, .5), Cumulative = c(1/3, 1)))
  expect_identical(distribution$display_data$Included, c(TRUE, TRUE, FALSE, FALSE, TRUE))
  expect_match(distribution$alt_text, "3 of 5")
  expect_match(distribution$notes$Text, "not a latent population")
  expect_identical(dev.cur(), before)
  path <- tempfile(); saveRDS(distribution, path)
  expect_identical(readRDS(path), distribution)
  for (arg in list(list(show_title = NA), list(reference = Inf), list(text_scale = 0),
      list(point_size = -1), list(title = FALSE), list(caption = NA), list(show_labels = 1))) {
    expect_error(do.call(plot, c(list(x, draw = FALSE), arg)), "must be")
  }
  expect_error(plot(x, draw = FALSE, style = "violin"), "arg")
})

test_that("both renderers honor controls and retain interpretation metadata", {
  x <- extended_view_fixture()
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(c("mar", "cex"))
  for (style in c("interval", "precision", "distribution")) {
    expect_silent(plot(x, style = style, palette = "mono", show_title = FALSE,
      show_notes = FALSE, show_labels = FALSE, reference = NULL, text_scale = 1.2))
    expect_equal(graphics::par(c("mar", "cex")), old)
    if (!requireNamespace("ggplot2", quietly = TRUE)) next
    raw <- plot(x, draw = FALSE, style = style, palette = "mono", title = "Custom title",
      caption = "Custom note", show_title = FALSE, show_notes = FALSE,
      show_labels = FALSE, reference = NULL, text_scale = 1.2, point_size = 4)
    p <- as_ggplot(raw)
    expect_null(p$labels$title); expect_null(p$labels$caption)
    expect_match(p$labels$alt, style)
    expect_match(paste(plot_data(p, "notes")$Text, collapse = "\n"), "calibration fixed")
    expect_equal(unique(attr(p, "mfrmr_display_data")$Colour), "#222222")
    expect_equal(p$theme$text$size, 13.2)
    expect_silent(ggplot2::ggplotGrob(p))
    expect_null(as_ggplot(x, style = style, title = "", caption = "")$labels$title)
  }
  x$table$Estimate <- x$table$Lower <- x$table$Upper <- NA_real_
  x$table$Status <- "unavailable"
  for (style in c("interval", "precision", "distribution")) {
    expect_silent(plot(x, style = style))
    if (requireNamespace("ggplot2", quietly = TRUE)) expect_silent(ggplot2::ggplotGrob(as_ggplot(x, style = style)))
  }
})

test_that("calibration and rater plots use the same view contract", {
  scores <- extended_view_fixture()
  fit <- structure(list(raters = transform(scores$table, Rater = Person),
    settings = list(), checks = list(NumericalReady = TRUE, InformationPositive = TRUE)), class = "mfrm_random_rater")
  fit$raters$PredictionSE <- .2
  expected <- fit$raters; expected$Lower <- expected$Upper <- NA_real_
  for (style in c("interval", "precision", "distribution")) {
    p <- plot(fit, style = style, draw = FALSE,
      intervals = if (style == "precision") "normal" else "none")
    if (style != "precision") {
      expect_identical(p$data$table[names(expected)], expected)
      expect_identical(p$data$table$IntervalStatus,
        c(rep("not_requested", 3), "invalid_estimate_or_se", "not_requested"))
    }
    expect_identical(p$data$checks, fit$checks)
    if (requireNamespace("ggplot2", quietly = TRUE)) expect_silent(ggplot2::ggplotGrob(as_ggplot(p)))
  }
  class(scores) <- "mfrm_testlet_scores"
  testlet <- structure(list(calibration_table = transform(fit$raters,
    Parameter = "Fixed facet", Facet = "Rater", Level = Rater, SE = PredictionSE),
    input = list(columns = list(facets = "Rater")), checks = fit$checks), class = "mfrm_testlet")
  for (object in list(testlet, scores)) for (style in c("precision", "distribution")) {
    args <- list(object, style = style, draw = FALSE)
    if (inherits(object, "mfrm_testlet")) args$intervals <- "normal"
    p <- do.call(plot, args)
    expect_identical(p$data$display$style, style)
    if (requireNamespace("ggplot2", quietly = TRUE)) expect_silent(ggplot2::ggplotGrob(as_ggplot(p)))
  }
  old <- new_mfrm_plot_data("random_rater_severity", list(table = fit$raters,
    settings = fit$settings, checks = fit$checks))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_match(paste(plot_data(as_ggplot(old), "notes")$Text, collapse = "\n"), "first[[:space:]]+order")
    old$data$checks$EstimatedVarianceBoundary <- TRUE
    expect_match(paste(plot_data(as_ggplot(old), "notes")$Text, collapse = "\n"), "zero")
  }
})

test_that("bootstrap conversion retains infinities and distinguishable comparison lines", {
  roots <- matrix(seq(-1, 1, length.out = 20), ncol = 1, dimnames = list(NULL, "R1"))
  roots[1:4, ] <- NA
  x <- structure(list(source = list(raters = data.frame(Rater = "R1", Estimate = .3, PredictionSE = .2)),
    studentized = roots, error = roots, settings = list(level = .8, nsim = 20),
    trials = data.frame(FitReady = rep(TRUE, 20))), class = "mfrm_random_rater_intervals")
  raw <- plot(x, draw = FALSE, palette = "mono", show_notes = FALSE)
  expect_identical(raw$data$table$Lower, -Inf)
  expect_identical(raw$data$table$Upper, Inf)
  expect_match(raw$data$notes$Text, "Solid: bootstrap; dashed: ordinary")
  expect_error(plot(x, style = "precision", draw = FALSE), "empty")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- as_ggplot(raw)
    expect_identical(p$data$Lower, -Inf)
    expect_equal(sum(vapply(p$layers, function(z) !is.null(z$geom_params$arrow), logical(1))), 2)
    expect_silent(ggplot2::ggplotGrob(p))
    legacy <- raw; legacy$data <- legacy$data[c("table", "settings", "availability")]
    expect_silent(ggplot2::ggplotGrob(as_ggplot(legacy)))
  }
})
