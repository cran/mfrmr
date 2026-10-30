cml_interval_display_fixture <- function() readRDS(test_path("fixtures", "cml-comparison.rds"))

test_that("CML comparison intervals reuse saved estimates, joint covariance and qualifications", {
  fit <- cml_interval_display_fixture()
  local_mocked_bindings(fit_mfrm = function(...) stop("unexpected refit"),
    mfrmr_cml_evaluate = function(...) stop("unexpected information"),
    mfrmr_cml_optimize = function(...) stop("unexpected optimization"))
  ci <- mfrm_facet_intervals(fit)
  expect_error(mfrm_results(fit), "does not support CML fits.*saved comparisons")
  expect_error(mfrm_results(fit, intervals = ci), "plot\\(ci\\).*saveRDS")
  for (name in c("Estimate", "SE", "Lower", "Upper", "ConditionalState", "EstimateReason"))
    expect_identical(ci$table[[name]], fit$comparisons[[name]])
  expect_identical(ci$table$CIEligible, c(FALSE, TRUE))
  expect_identical(ci$table$InferenceReview, ci$table$EstimateReason)
  expect_identical(ci$table$Status, c("unavailable", "available"))
  expect_identical(ci$table$IdentificationBasis, c("unidentified", "conditional_support"))
  expect_false(any(ci$table$CoverageQualified))
  expect_identical(ci$covariance, fit$covariance)
  expect_identical(ci$parameter_covariance, fit$optimization$covariance)
  expect_identical(ci$comparison_review, fit$comparison_review)
  expect_identical(ci$contrasts, fit$comparison_review$weights)
  expect_identical(ci$fit_checks, fit$checks)
  expect_identical(ci$fit, fit)
  for (level in c(.9, .995, 1 - .Machine$double.eps)) {
    other <- mfrm_facet_intervals(fit, level = level)
    expect_identical(other$table$Estimate, ci$table$Estimate)
    expect_identical(other$table$SE, ci$table$SE)
    expect_identical(other$covariance, ci$covariance)
    expect_equal(other$table$Upper[2], ci$table$Estimate[2] + qnorm((1-level)/2, lower.tail=FALSE)*ci$table$SE[2])
    expect_true(all(is.na(unlist(other$table[1, c("Estimate", "SE", "Lower", "Upper")]))))
    expect_match(other$cautions[3], paste0(format(100*level, trim=TRUE), "%"), fixed=TRUE)
  }
  expect_identical(summary(ci), ci$table)
  saved <- tempfile(fileext = ".rds")
  withr::defer(unlink(saved))
  saveRDS(ci, saved)
  expect_identical(readRDS(saved), ci)
  expect_equal(apa_table(ci, digits = 10)$table$SE, ci$table$SE, tolerance=1e-9)
  expect_identical(apa_table(ci)$table$CIEligible, c(FALSE, TRUE))
  expect_equal(apa_table(ci, which="comparison_cells")$table, ci$comparison_review$cells)
  expect_identical(apa_table(ci, which="comparison_review")$table$ConditionalState,
    fit$comparisons$ConditionalState)
})

test_that("CML display cannot change the fitted targets or supply an absent population SD", {
  fit <- cml_interval_display_fixture()
  expect_error(mfrm_facet_intervals(fit, scale="standardized"), "no.*population SD")
  expect_error(mfrm_facet_intervals(fit, method="sandwich"), "method = 'model'")
  expect_error(mfrm_facet_intervals(fit, clusters=data.frame()), "clusters = NULL")
  expect_error(mfrm_facet_intervals(fit, adjust=TRUE), "adjust = FALSE")
  expect_error(mfrm_facet_intervals(fit, contrasts=fit$comparison_review$weights), "saved targets")
  expect_error(mfrm_facet_intervals(fit, comparison_review=fit$comparison_review), "saved targets")
  expect_error(mfrm_facet_intervals(fit, facet="Criterion"), "saved targets")
  expect_identical(mfrm_facet_intervals(fit, facet="Rater"), mfrm_facet_intervals(fit))
  for (level in list(0, 1, NA, Inf, c(.9,.95), "95", .95+0i, matrix(.95)))
    expect_error(mfrm_facet_intervals(fit, level=level), "0 < level < 1")
  x <- readRDS(test_path("fixtures", "partial-comparisons.rds"))[["RSM-MML"]]
  expect_error(mfrm_facet_intervals(x$fit, comparison_review=x$review), "full native")
})

test_that("missing CML comparisons remain labeled without invented points in both plot routes", {
  skip_if_not_installed("ggplot2")
  fit <- cml_interval_display_fixture()
  # All-extreme totals are a real no-information fit, not an optimizer failure.
  data <- load_mfrmr_data("example_bias"); data$Score <- ifelse(data$Group == "A", 1L, 4L)
  review <- describe_mfrm_data(data, "Person", c("Rater", "Criterion", "Group"), "Score",
    rating_min=1, rating_max=4, category_policy="preserve", include_agreement=FALSE,
    comparisons=fit$comparison_review$specification)
  unavailable <- do.call(fit_mfrm, c(list(data=data, fit_scope="comparisons",
    comparison_review=review, interaction_policy="silent"), review$comparison_review$fit_arguments))
  expect_null(unavailable$optimization)
  local_mocked_bindings(fit_mfrm = function(...) stop("unexpected refit"),
    mfrmr_cml_evaluate = function(...) stop("unexpected information"))
  path <- tempfile(fileext=".pdf")
  grDevices::pdf(path, width=7, height=4.2)
  withr::defer({grDevices::dev.off(); unlink(path)})
  for (source in list(fit, unavailable)) {
    ci <- mfrm_facet_intervals(source)
    for (reference in list(0, NULL)) {
      expect_no_warning(p <- plot(ci, reference=reference))
      expect_identical(plot_data(p, "table"), ci$table)
      expect_identical(plot_data(p, "comparison_review"), source$comparison_review)
      expect_identical(plot_data(p, "fit_checks"), source$checks)
      labels <- plot_data(p, "target_labels")
      expect_true(all(grepl("not estimated", labels[is.na(ci$table$Estimate)])))
      expect_match(plot_data(p, "subtitle"), "Conditional information")
      expect_match(plot_data(p, "xlab"), "Interaction contrast")
      gg <- as_ggplot(p)
      expect_no_warning(print(gg))
      expect_identical(plot_data(gg, "table"), ci$table)
      built <- ggplot2::ggplot_build(gg)
      expect_equal(nrow(built$data[[2L]]), sum(is.finite(ci$table$Estimate)))
      axis <- built$layout$panel_params[[1]]$y
      expect_identical(axis$get_labels(), labels)
      expect_true(all(is.finite(axis$get_breaks())))
    }
  }
})
