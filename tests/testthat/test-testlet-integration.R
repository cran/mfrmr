testlet_display_fixture <- function() {
  structure(list(table = data.frame(Person = c("Measured", "Missing", "Unresolved"),
    Estimate = c(.4, 0, NA), Lower = c(-.5, -1.96, NA), Upper = c(1.2, 1.96, NA),
    Status = c("available_conditional", "prior_only", "unavailable"),
    Reason = c("", "", "Integration unresolved")),
    settings = list(level = .95, calibration_uncertainty = FALSE)), class = "mfrm_testlet_scores")
}

test_that("testlet display data retain interpretation and work with shared accessors", {
  scores <- testlet_display_fixture(); before <- dev.cur()
  table <- plot_data(scores, component = "table")
  expect_identical(table, scores$table)
  expect_identical(dev.cur(), before)
  components <- plot_data_components(scores)
  expect_true(all(c("table", "settings", "notes", "caption") %in% components$Component))
  payload <- plot(scores, draw = FALSE)
  expect_match(plot_data(payload)$caption, "calibration fixed")
  expect_match(plot_data(payload)$notes$Text, "prior only")
  path <- tempfile(); saveRDS(payload, path)
  expect_identical(plot_data(readRDS(path)), plot_data(payload))
})

test_that("testlet ggplots preserve intervals, unavailable rows and prior-only symbols", {
  skip_if_not_installed("ggplot2")
  scores <- testlet_display_fixture(); before <- dev.cur()
  p <- as_ggplot(scores)
  expect_identical(dev.cur(), before)
  expect_s3_class(p, "ggplot")
  expect_identical(p$data[names(scores$table)], scores$table)
  expect_equal(p$layers[[2]]$data$Lower, scores$table$Lower[1:2])
  expect_equal(p$layers[[2]]$data$Upper, scores$table$Upper[1:2])
  expect_equal(p$layers[[3]]$data$.Shape, c(16, 1))
  expect_match(paste(plot_data(p, "notes")$Text, collapse = "\n"), "calibration fixed")
  expect_match(attr(p, "mfrmr_notes")$Text, "prior only")
  b <- ggplot2::ggplot_build(p)
  expect_identical(b$layout$panel_params[[1]]$y$get_labels(), scores$table$Person)
  expect_silent(ggplot2::ggplotGrob(p))
  draw_free <- plot(scores, draw = FALSE)
  expect_identical(as_ggplot(draw_free, component = "table")$labels, p$labels)
  expect_error(as_ggplot(draw_free, component = "settings"), "plot_data")
  expect_error(as_ggplot(draw_free, level = .9), "empty")
  old <- new_mfrm_plot_data("testlet_scores", list(table = scores$table, settings = scores$settings))
  expect_identical(as_ggplot(old)$labels, p$labels)
  expect_equal(as_ggplot(old)$data, p$data)
  scores$table[, c("Estimate", "Lower", "Upper")] <- NA_real_
  scores$table$Status <- "unavailable"
  empty <- as_ggplot(scores)
  expect_equal(nrow(empty$data), 3)
  expect_equal(nrow(empty$layers[[3]]$data), 0)
  expect_silent(ggplot2::ggplotGrob(empty))
})

test_that("fixed-facet testlet plots and unsupported ordinary-model routes stay distinct", {
  tab <- data.frame(Parameter = "Fixed facet", Facet = "Rater", Level = c("A", "B"),
    Estimate = c(-.3,.3), Lower = NA_real_, Upper = NA_real_)
  fit <- structure(list(calibration_table = tab,
    input = list(columns = list(facets = "Rater")), checks = list(EstimatedVarianceBoundary = TRUE)),
    class = "mfrm_testlet")
  displayed <- plot_data(fit, "table", facet = "Rater")
  expect_identical(displayed[names(tab)], tab)
  expect_true(all(displayed$IntervalStatus == "numerical_review_failed"))
  expect_false(any(displayed$CIEligible))
  expect_true(all(nzchar(displayed$InferenceReview)))
  for (object in list(fit, structure(list(), class = "mfrm_random_rater"))) {
    expect_error(diagnose_mfrm(object), "numerical checks are not model-fit diagnostics", fixed = TRUE)
  }
  expect_error(mfrm_results(testlet_display_fixture()), "testlet")
  expect_error(mfrm_report(fit), "mfrm_results")
  expect_error(export_mfrm_results(fit, output_dir = tempfile()), "mfrm_results")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- as_ggplot(fit, facet = "Rater")
    expect_equal(nrow(p$layers[[2]]$data), 0)
    expect_equal(nrow(p$layers[[3]]$data), 2)
    expect_match(paste(plot_data(p, "notes")$Text, collapse = "\n"), "bounds not requested")
    expect_silent(ggplot2::ggplotGrob(p))
  }
})

test_that("the model guide distinguishes prediction targets and preserves the beginner route", {
  guide <- mfrmr_output_guide("models")
  expect_identical(guide$MainFunction, c("fit_mfrm()", "fit_mfrm_random_rater()", "fit_mfrm_testlet()"))
  expect_match(guide$NextStep[2], "probabilities at supplied abilities")
  expect_match(guide$NextStep[3], "score_mfrm_persons()", fixed = TRUE)
  expect_match(guide$DecisionBoundary[3], "exclude calibration uncertainty")
  expect_identical(guide$Lifecycle, c("stable", "advanced", "advanced"))
  expect_equal(nrow(mfrmr_output_guide("beginner")), 6)
  expect_true(all(guide$Question %in% mfrmr_output_guide()$Question))
})
