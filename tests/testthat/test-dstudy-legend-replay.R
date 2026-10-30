test_that("D-study legends recalculate text spacing on a resized recorded plot", {
  g <- mfrm_g_study(mfrmr_example_gtheory_nested, "Content", nesting = c(Rater = "Task"))
  d <- mfrm_d_study(g, expand.grid(Raters = c(2, 3), Tasks = c(4, 6)))
  paths <- c(tempfile(fileext = ".pdf"), tempfile(fileext = ".pdf"))
  on.exit(unlink(paths), add = TRUE)
  calls <- list()
  original <- graphics::legend
  testthat::local_mocked_bindings(legend = function(...) {
    value <- original(...)
    y <- sort(unique(value$text$y))
    calls[[length(calls) + 1L]] <<- c(spacing = min(diff(y)),
      height = graphics::strheight("M", cex = .8), device = graphics::par("din")[2])
    value
  }, .package = "mfrmr")
  grDevices::pdf(paths[1], width = 7, height = 7)
  device <- grDevices::dev.cur()
  withr::defer(if (device %in% grDevices::dev.list()) grDevices::dev.off(device))
  grDevices::dev.control("enable")
  expected <- plot(d, x_var = "Raters")
  saved <- unserialize(serialize(grDevices::recordPlot(), NULL))
  grDevices::dev.off()
  count <- length(calls)
  grDevices::pdf(paths[2], width = 7, height = 4.33)
  grDevices::replayPlot(saved)
  grDevices::dev.off()
  expect_equal(length(calls) - count, 2L)
  replay <- do.call(rbind, tail(calls, 2L))
  expect_true(all(replay[, "spacing"] > replay[, "height"]))
  expect_equal(replay[, "device"], c(4.33, 4.33))
  expect_identical(plot_data(plot(d, x_var = "Raters", draw = FALSE)), plot_data(expected))
})

test_that("D-study interval forests leave room for numbered and nested plan labels", {
  paths <- c(tempfile(fileext = ".pdf"), tempfile(fileext = ".pdf"))
  on.exit(unlink(paths), add = TRUE)
  calls <- list()
  original <- graphics::axis
  testthat::local_mocked_bindings(axis = function(side, at, labels, ...) {
    if (side == 2 && is.character(labels)) calls[[length(calls)+1L]] <<- c(
      margin = graphics::par("mai")[2],
      required = max(graphics::strwidth(labels, units = "inches", cex = .8))+
        graphics::par("mgp")[2]*graphics::par("csi")*graphics::par("mex"))
    original(side = side, at = at, labels = labels, ...)
  }, .package = "mfrmr")
  for (nested in c(FALSE, TRUE)) {
    g <- if (nested) mfrm_g_study(mfrmr_example_gtheory_nested, "Content",
      nesting = c(Rater = "Task")) else
      mfrm_g_study(mfrmr_example_gtheory_crossed, "Content")
    d <- mfrm_d_study(g, data.frame(Raters = c(2, 2, 4), Tasks = c(3, 3, 6)))
    ci <- confint(d, method = "mean_square_region", assumption = "normal")
    grDevices::pdf(paths[1L+nested], width = 7, height = 4.5)
    device <- grDevices::dev.cur()
    withr::defer(if (device %in% grDevices::dev.list()) grDevices::dev.off(device))
    expect_identical(plot(ci), plot(ci, draw = FALSE))
    grDevices::dev.off()
  }
  expect_length(calls, 2L)
  measured <- do.call(rbind, calls)
  expect_true(all(measured[, "margin"] > measured[, "required"]))
})
