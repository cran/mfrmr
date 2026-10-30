test_that("portable score conversion preserves its selected view with explicit data", {
  skip_if_not_installed("ggplot2")
  data <- data.frame(Person = c("001", "NA"), Estimate = c(-.5, .2),
    SD = c(.2, .1), Lower = c(-1.1, -.1), Upper = c(.25, .35),
    Observations = c(4, 12), QuadratureEdgeMass = c(.01, .03),
    ReviewFlag = c(FALSE, TRUE), DisplayOrder = c(2, 1))
  local_mocked_bindings(score_mfrm_calibration = function(...) stop("No scoring"),
    fit_mfrm = function(...) stop("No fitting"))
  for (type in c("interval", "precision", "edge_mass")) {
    p <- new_mfrm_plot_data(paste0("calibration_score_", type), list(data = data,
      settings = list(type = type, label_review = FALSE, interval_level = .9),
      unplotted_dispositions = data.frame(Person = "missing", Disposition = "not_scored"),
      caption = "Conditional on frozen calibration; calibration uncertainty excluded.",
      reference_lines = if(type == "edge_mass") new_reference_lines("h", .02, "Review", "dashed", "review") else new_reference_lines()))
    g <- as_ggplot(p, component = "data")
    b <- ggplot2::ggplot_build(g)$data
    expect_equal(b, ggplot2::ggplot_build(as_ggplot(p))$data)
    expect_identical(plot_data(g), p$data)
    expect_identical(plot_data(g, "unplotted_dispositions")$Person, "missing")
    expect_null(g$labels$title); expect_null(g$labels$caption)
    if(type == "interval") {
      expect_equal(b[[1]]$x, data$Lower); expect_equal(b[[1]]$xend, data$Upper)
      expect_equal(b[[2]]$x, data$Estimate)
      expect_equal(b[[2]]$shape, c(16, 17))
    } else {
      expect_equal(b[[1]]$x, if(type == "precision") data$Observations else data$Estimate)
      expect_equal(b[[1]]$y, if(type == "precision") data$SD else data$QuadratureEdgeMass)
      expect_equal(b[[1]]$shape, c(16, 17))
      if(type == "edge_mass") expect_equal(b[[2]]$yintercept, .02)
    }
    expect_error(as_ggplot(p, component = "unplotted_dispositions"), "plot_data")
    expect_error(as_ggplot(p, interval_level = .8), "Unsupported arguments")
  }
})
