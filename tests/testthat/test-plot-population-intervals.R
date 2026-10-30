population_plot_fixture <- function() {
  tab <- data.frame(Target = c("beta:(Intercept)", "beta:GroupB", "contrast:B minus A", "log_residual_variance"),
    Type = c("coefficient", "coefficient", "coefficient_contrast", "residual_variance"),
    Term = c("(Intercept)", "GroupB", "B minus A", "Residual variance"),
    Estimate = c(.2, -.4, -.4, 1.2), SE = c(.1, .2, .2, .4),
    Lower = c(.03, -.73, -.73, .6), Upper = c(.37, -.07, -.07, 2.4),
    CIEligible = TRUE, Status = "available", InferenceReview = "Experimental test fixture")
  structure(list(table = tab, settings = list(level = .9, model = "RSM", method = "model"),
    contrasts = matrix(c(0, 1), 1, dimnames = list("B minus A", c("(Intercept)", "GroupB"))),
    fit = list(population = list(formula = ~ Group, terms = terms(~ Group),
      xlevels = list(Group = c("A", "B")), contrasts = list(Group = "contr.treatment"),
      design_columns = c("(Intercept)", "GroupB"))),
    checks = data.frame(Check = "Fixture numerical check", Passed = TRUE),
    numerical_checks = data.frame(FullGradient = 0),
    cautions = "Experimental source qualification.", source = list(id = "plot-fixture")),
    class = "mfrm_population_intervals")
}

test_that("population plots select saved targets and retain interpretation metadata", {
  x <- population_plot_fixture(); original <- serialize(x, NULL)
  local_mocked_bindings(fit_mfrm = function(...) stop("No fitting"),
    mfrm_population_intervals = function(...) stop("No new intervals"),
    compute_mml_parameter_covariance = function(...) stop("No information calculation"), .package = "mfrmr")
  for (type in c("coefficients", "contrasts")) {
    p <- plot(x, type = type, draw = FALSE)
    selected <- if (type == "coefficients") 1:2 else 3L
    expect_identical(p$data$table[names(x$table)], x$table[selected, ])
    expect_identical(p$data$intervals, x$table)
    expect_identical(p$data$contrast_weights, x$contrasts)
    expect_identical(p$data$population_spec, x$fit$population)
    expect_identical(p$data$checks, x$checks)
    expect_identical(p$data$numerical_checks, x$numerical_checks)
    expect_identical(p$data$source, x$source)
    expect_match(p$data$population_coding$Levels, "A")
    expect_match(paste(p$data$notes$Text, collapse = " "), "Experimental")
    expect_false(any(p$data$table$Type == "residual_variance"))
    expect_true(all(p$data$table$PlotStatus == "Interval available"))
  }
  expect_identical(serialize(x, NULL), original)
  expect_identical(summary(x), x$table)
  expect_output(print(x), "Conditional-normal population")
})

test_that("population plot conversion preserves unavailable rows and stored bounds", {
  skip_if_not_installed("ggplot2")
  x <- population_plot_fixture()
  x$table$CIEligible[1] <- FALSE
  x$table$Status[1] <- "numerical_review_failed"
  x$table$InferenceReview[1] <- "Check failed; retain the point"
  x$table$Upper[2] <- NA_real_
  p <- plot(x, draw = FALSE)
  expect_identical(p$data$table[names(x$table)], x$table[1:2, ])
  expect_true(all(p$data$table$PlotStatus == "CI unavailable"))
  expect_true(all(grepl("CI unavailable", p$data$table$Label)))
  saved <- unserialize(serialize(p, NULL))
  g <- as_ggplot(saved)
  layers <- ggplot2::ggplot_build(g)$data
  expect_equal(nrow(layers[[1]]), 0L)
  expect_equal(layers[[2]]$x, x$table$Estimate[1:2])
  expect_equal(layers[[2]]$shape, c(4, 4))
  expect_identical(plot_data(g), plot_data(saved))
  expect_identical(saved$data$table, p$data$table)
  # Serialization can recreate the fixture's test environment; compare its coding.
  design <- function(spec) model.matrix(spec$terms, data.frame(Group = c("A", "B")),
    contrasts.arg = spec$contrasts)
  expect_identical(design(saved$data$population_spec), design(p$data$population_spec))
  expect_null(g$labels$title)
  expect_null(g$labels$caption)
  expect_identical(plot_data(as_ggplot(x, type = "contrasts"))$table$Target, x$table$Target[3])
  expect_error(as_ggplot(p, component = "table"), "plot_data")
})

test_that("population plots keep missing points without inventing zeros", {
  x <- population_plot_fixture()
  x$table$Estimate[1:2] <- NA_real_
  x$table$Lower[1:2] <- x$table$Upper[1:2] <- NA_real_
  p <- plot(x, reference = NULL, draw = FALSE)
  expect_identical(p$data$table$PlotStatus, rep("Estimate unavailable", 2))
  expect_identical(p$data$x_range, c(-1, 1))
  expect_true(all(is.na(p$data$table$Estimate)))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    layers <- ggplot2::ggplot_build(as_ggplot(p))$data
    expect_true(all(vapply(layers, nrow, integer(1)) == 0))
  }
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(x, reference = NULL))
})

test_that("population plot controls and guidance cannot silently change the target", {
  x <- population_plot_fixture()
  x$contrasts <- NULL; x$table <- x$table[x$table$Type != "coefficient_contrast", ]
  expect_error(plot(x, type = "contrasts", draw = FALSE), "No saved contrasts")
  expect_error(plot(x, type = "variance", draw = FALSE), "arg")
  expect_error(plot(x, draw = NA), "TRUE or FALSE")
  expect_error(plot(x, reference = Inf, draw = FALSE), "finite number")
  expect_error(plot(x, ci_level = .8, draw = FALSE), "must be empty")
  bad <- x; bad$table$CIEligible <- NULL
  expect_error(plot(bad, draw = FALSE), "qualification columns")
  withr::local_options(mfrmr.plot_preset = "monochrome")
  expect_identical(plot(x, draw = FALSE)$data$preset, "monochrome")
  guide <- subset(mfrmr_interval_guide("visual"), Route == "Conditional-normal population intervals")
  expect_equal(nrow(guide), 1L)
  expect_match(guide$DisplayRoute, "plot(ci)", fixed = TRUE)
  routes <- subset(mfrmr_output_guide("plots"), InputClass == "mfrm_population_intervals")
  expect_equal(nrow(routes), 2L)
  expect_true(all(routes$GGPlot == "dedicated"))
})
