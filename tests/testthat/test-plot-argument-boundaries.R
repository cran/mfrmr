test_that("saved plot types cannot be silently changed in dedicated converters", {
  local_mocked_bindings(.require_mfrmr_ggplot2 = function() stop("unexpected conversion"))
  for (name in c("extended_model_map", "response_diagnostics", "screening_sensitivity",
    "extended_model_comparison", "testlet_calibration", "testlet_scores",
    "random_rater_scores", "person_scores", "random_rater_severity",
    "random_rater_interval_bootstrap", "multivariate_d_study")) {
    payload <- new_mfrm_plot_data(name, list())
    expect_error(as_ggplot(payload, type = "another_view"), "already selects a view")
  }
})

local({
  fit <- make_toy_fit(method = "MML", maxit = 500)
  diagnostics <- make_toy_diagnostics(fit)
  res <- mfrm_results(fit, diagnostics = diagnostics, compute = "never")

  test_that("named level changes fitted-plot intervals through each entry", {
    old <- plot(fit, diagnostics = diagnostics, ci_level = .8, draw = FALSE)
    new <- plot(fit, diagnostics = diagnostics, level = .8, draw = FALSE)
    expect_identical(new, old)
    expect_identical(plot(res, level = .8, draw = FALSE), old)
    expect_identical(plot_data(res, type = "wright", component = "locations", level = .8),
      old$data$locations)
    default <- plot(fit, diagnostics = diagnostics, draw = FALSE)
    expect_identical(default,
      plot(fit, diagnostics = diagnostics, level = .95, draw = FALSE))
    # Verify actual geometry, not just stored confidence-level labels.
    width80 <- old$data$locations$CI_Upper - old$data$locations$CI_Lower
    width95 <- default$data$locations$CI_Upper - default$data$locations$CI_Lower
    usable <- is.finite(width80) & is.finite(width95) & width95 > 0
    expect_true(any(usable))
    expect_equal(width80[usable] / width95[usable],
      rep(qnorm(.9) / qnorm(.975), sum(usable)))
    skip_if_not_installed("ggplot2")
    expect_equal(ggplot2::ggplot_build(as_ggplot(res, level = .8))$data,
      ggplot2::ggplot_build(as_ggplot(old))$data)
  })

  test_that("ambiguous or invalid level requests fail explicitly", {
    expect_error(plot(fit, ci_level = .8, level = .8, draw = FALSE), "only one")
    for (bad in list(NULL, NA_real_, .8 + 0i, c(.8, .9), 0, 1, Inf, "0.8")) {
      expect_error(plot(fit, level = bad, draw = FALSE), "single number")
    }
    ci <- mfrm_facet_intervals(fit, "Rater", level = .8)
    saved <- mfrm_results(fit, intervals = list(raters = ci), compute = "never")
    expect_error(plot(saved, type = "facet_raters", level = .9, draw = FALSE),
      "must be empty")
    expect_identical(plot(saved, type = "facet_raters", draw = FALSE)$data$table, ci$table)
  })

  test_that("saved conversion rejects ignored settings and retains supported CCC styling", {
    skip_if_not_installed("ggplot2")
    payload <- plot(fit, diagnostics = diagnostics, draw = FALSE)
    for (args in list(list(level = .8), list(ci_level = .8), list(title = NULL),
      list(preset = "monochrome"), list(slope_aes = "linewidth"), list(.8))) {
      expect_error(do.call(as_ggplot, c(list(payload), args)), "saved plot")
    }
    expect_error(as_ggplot(payload, type = "ccc"), "already selects")
    expect_identical(as_ggplot(unserialize(serialize(payload, NULL)))$data,
      as_ggplot(payload)$data)
    ccc <- plot(fit, type = "ccc", draw = FALSE)
    expect_s3_class(as_ggplot(ccc, slope_aes = "linewidth"), "ggplot")
    expect_error(as_ggplot(ccc, slope_aes = "linewidth", title = "Ignored"), "saved plot")
    expect_error(as_ggplot(ccc, component = "probabilities", facet_by = "none"), "saved plot")
    expect_error(do.call(as_ggplot,
      c(list(ccc), setNames(list("none", "none"), c("slope_aes", "slope_aes")))), "saved plot")
  })
})
