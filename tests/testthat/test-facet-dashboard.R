facet_dashboard_fixture <- local({
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

  dat <- mfrmr:::sample_mfrm_data(seed = 321)
  fit <- suppressWarnings(mfrmr::fit_mfrm(
    data = dat,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 20
  ))
  diagnostics <- suppressWarnings(mfrmr::diagnose_mfrm(fit, residual_pca = "none"))
  bias_rater <- suppressWarnings(mfrmr::estimate_bias(
    fit,
    diagnostics,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))
  bias_task <- suppressWarnings(mfrmr::estimate_bias(
    fit,
    diagnostics,
    facet_a = "Task",
    facet_b = "Criterion",
    max_iter = 2
  ))
  dashboard_single <- mfrmr::facet_quality_dashboard(
    fit,
    diagnostics = diagnostics,
    bias_results = bias_rater
  )
  dashboard_list <- mfrmr::facet_quality_dashboard(
    fit,
    diagnostics = diagnostics,
    bias_results = list(rater_criterion = bias_rater, task_criterion = bias_task)
  )

  list(
    fit = fit,
    diagnostics = diagnostics,
    bias_rater = bias_rater,
    bias_task = bias_task,
    dashboard_single = dashboard_single,
    dashboard_list = dashboard_list
  )
})

test_that("facet_quality_dashboard constructs a dashboard bundle with inferred facet", {
  expect_s3_class(facet_dashboard_fixture$dashboard_single, "mfrm_facet_dashboard")
  expect_identical(facet_dashboard_fixture$dashboard_single$facet, "Rater")
  expect_true(all(c("summary", "detail", "flagged", "settings") %in% names(facet_dashboard_fixture$dashboard_single)))
  expect_true(is.data.frame(facet_dashboard_fixture$dashboard_single$overview))
  expect_true(is.data.frame(facet_dashboard_fixture$dashboard_single$summary))
})

test_that("facet_quality_dashboard handles single and named-list bias bundles", {
  dash_single <- facet_dashboard_fixture$dashboard_single
  dash_list <- facet_dashboard_fixture$dashboard_list

  expect_identical(dash_single$detail$BiasCount, dash_list$detail$BiasCount)
  expect_identical(dash_single$detail$BiasSources, dash_list$detail$BiasSources)
  expect_true(any(!dash_list$bias_sources$Used))
  expect_true(any(grepl("target facet not involved", dash_list$bias_sources$Reason, fixed = TRUE)))
  expect_identical(sum(dash_single$detail$AnyFlag, na.rm = TRUE), nrow(dash_single$flagged))
})

test_that("facet_quality_dashboard surfaces failed bias-collection pairs in notes", {
  bias_collection <- structure(
    list(
      by_pair = list(rater_criterion = facet_dashboard_fixture$bias_rater),
      errors = data.frame(
        Interaction = "Task x Criterion",
        Facets = "Task x Criterion",
        Error = "forced pair failure",
        stringsAsFactors = FALSE
      )
    ),
    class = c("mfrm_bias_collection", "mfrm_bundle", "list")
  )

  dash <- mfrmr::facet_quality_dashboard(
    facet_dashboard_fixture$fit,
    diagnostics = facet_dashboard_fixture$diagnostics,
    bias_results = bias_collection
  )

  expect_true(any(grepl("^pair error:", dash$bias_sources$Reason)))
  expect_true(any(grepl("failed", dash$notes, fixed = TRUE)))
})

test_that("summary() returns a compact facet dashboard summary", {
  sum_dash <- summary(facet_dashboard_fixture$dashboard_single, top_n = 5)

  expect_s3_class(sum_dash, "summary.mfrm_facet_dashboard")
  expect_identical(sum_dash$summary_kind, "facet_dashboard")
  expect_true(is.data.frame(sum_dash$overview))
  expect_true(nrow(sum_dash$overview) == 1)
  expect_true(nrow(sum_dash$preview) <= 5)
})

test_that("plot_facet_quality_dashboard returns mfrm_plot_data for severity and flags", {
  severity_plot <- plot_facet_quality_dashboard(
    facet_dashboard_fixture$dashboard_single,
    plot_type = "severity",
    draw = FALSE
  )
  flags_plot <- plot_facet_quality_dashboard(
    facet_dashboard_fixture$fit,
    diagnostics = facet_dashboard_fixture$diagnostics,
    bias_results = facet_dashboard_fixture$bias_rater,
    plot_type = "flags",
    draw = FALSE
  )

  expect_s3_class(severity_plot, "mfrm_plot_data")
  expect_s3_class(flags_plot, "mfrm_plot_data")
  expect_identical(severity_plot$name, "facet_quality_dashboard")
  expect_identical(severity_plot$data$plot, "severity")
  expect_identical(flags_plot$data$plot, "flags")
  expect_true(is.data.frame(severity_plot$data$table))
  expect_true(is.data.frame(flags_plot$data$table))
})
