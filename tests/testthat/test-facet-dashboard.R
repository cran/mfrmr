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
  expect_false(any(facet_dashboard_fixture$dashboard_single$detail$CentralTendencyFlag))
  expect_true(any(grepl("disabled by default", facet_dashboard_fixture$dashboard_single$notes, fixed = TRUE)))
})

test_that("legacy central-tendency cutoff is labelled as origin proximity", {
  dash <- mfrmr::facet_quality_dashboard(
    facet_dashboard_fixture$fit,
    diagnostics = facet_dashboard_fixture$diagnostics,
    central_tendency_max = 0.25
  )

  expect_identical(
    dash$detail$CentralTendencyFlag,
    is.finite(abs(dash$detail$Estimate)) & abs(dash$detail$Estimate) <= 0.25
  )
  expect_true(any(grepl("origin-proximity marker only", dash$notes, fixed = TRUE)))
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
      fit_identity = facet_dashboard_fixture$bias_rater$fit_identity,
      calculation_version = facet_dashboard_fixture$bias_rater$calculation_version,
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


test_that("unavailable diagnostics cannot silently pass screening", {
  fit <- facet_dashboard_fixture$fit
  dx <- facet_dashboard_fixture$diagnostics
  dx$measures <- as.data.frame(dx$measures)
  rows <- which(dx$measures$Facet == "Rater")
  # Keep the missing estimate consistent with the fitted source being displayed.
  level <- dx$measures$Level[rows[1]]
  fit$facets$others$Estimate[fit$facets$others$Facet == "Rater" &
                            fit$facets$others$Level == level] <- NA_real_
  dx$measures$Estimate[rows[1]] <- NA_real_
  dx$parameter_uncertainty$source <- mfrm_native_location_identity(fit)
  dx$measures$Infit[rows] <- 1
  dx$measures$Outfit[rows] <- NA_real_
  dx$measures$Infit[rows[1]] <- 2
  dash <- facet_quality_dashboard(fit, diagnostics = dx)
  expect_true(is.na(dash$detail$SeverityFlag[1]))
  expect_true(dash$detail$MisfitFlag[1])
  expect_true(all(is.na(dash$detail$MisfitFlag[-1])))
  expect_true(all(grepl("Outfit", dash$detail$MissingMetrics)))
  expect_equal(dash$overview$IncompleteLevels, length(rows))
  p <- plot_facet_quality_dashboard(dash, plot_type = "flags", draw = FALSE)
  expect_true(dash$detail$Level[1] %in% p$data$table$Level)
  expect_identical(p$data$notes$Text, dash$notes)
})

test_that("nonfinite fit indices cannot hide or create observed misfit flags", {
  dx <- facet_dashboard_fixture$diagnostics
  dx$measures <- as.data.frame(dx$measures)
  rows <- which(dx$measures$Facet == "Rater")
  expected <- rep(c(TRUE, TRUE, NA), length.out = length(rows))
  for (unavailable in c(NA_real_, NaN, Inf, -Inf)) {
    for (missing_metric in c("Infit", "Outfit")) {
      observed_metric <- setdiff(c("Infit", "Outfit"), missing_metric)
      dx$measures[[missing_metric]][rows] <- unavailable
      dx$measures[[observed_metric]][rows] <- rep(c(2, 0.2, 1), length.out = length(rows))
      dash <- facet_quality_dashboard(facet_dashboard_fixture$fit, diagnostics = dx)
      expect_identical(dash$detail$MisfitFlag, expected)
      expect_true(all(grepl(missing_metric, dash$detail$MissingMetrics, fixed = TRUE)))
    }
  }
})

test_that("dashboard plots and saved outputs retain the screening basis", {
  dx <- facet_dashboard_fixture$diagnostics
  dx$measures <- as.data.frame(dx$measures)
  rows <- which(dx$measures$Facet == "Rater")
  dx$measures$Infit[rows] <- 0.6
  dx$measures$Outfit[rows] <- 1
  dash <- facet_quality_dashboard(facet_dashboard_fixture$fit, diagnostics = dx,
                                  severity_warn = 2.25)
  p <- plot_facet_quality_dashboard(dash, severity_warn = 99, draw = FALSE)
  expect_equal(p$data$thresholds$severity_warn, 2.25)
  expect_equal(p$data$thresholds$misfit_lower, 0.5)
  direct <- plot_facet_quality_dashboard(facet_dashboard_fixture$fit,
                                         diagnostics = dx, draw = FALSE)
  expect_false(any(direct$data$table$MisfitFlag))
  expect_identical(p$data$fit_readiness, dash$fit_readiness)
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  saveRDS(p, saved)
  expect_identical(readRDS(saved), p)
  expect_output(print(p), "screening prompts")

  out_dir <- tempfile("dashboard-export-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  bundle <- export_mfrm_bundle(facet_dashboard_fixture$fit, diagnostics = dx,
    output_dir = out_dir, prefix = "feedback", include = c("dashboard", "html"),
    acknowledge_sensitive = TRUE)
  settings <- read.csv(file.path(out_dir, "feedback_facet_dashboard_settings.csv"))
  expect_equal(as.numeric(settings$Value[settings$Setting == "misfit_lower"]), 0.5)
  notes <- readLines(file.path(out_dir, "feedback_facet_dashboard_notes.txt"))
  expect_true(any(grepl("zero flags does not mean", notes, fixed = TRUE)))
  html <- paste(readLines(file.path(out_dir, "feedback_bundle.html")), collapse = "\n")
  expect_match(html, "zero flags does not mean", fixed = TRUE)
  expect_match(html, "misfit_lower", fixed = TRUE)
  expect_true(all(c("dashboard_settings", "dashboard_notes") %in%
                    bundle$written_files$Component))
})


test_that("dashboard HTML retains unflagged missing diagnostics without internal columns", {
  dx <- facet_dashboard_fixture$diagnostics
  dx$measures <- as.data.frame(dx$measures)
  rows <- which(dx$measures$Facet == "Rater")
  dx$measures$Infit[rows] <- 1
  dx$measures$Outfit[rows] <- NA_real_
  dx$measures$ReadinessContractVersion <- "dashboard-internal-test"
  dash <- facet_quality_dashboard(facet_dashboard_fixture$fit, diagnostics = dx,
                                  facet = "Rater")
  expect_true(any(!dash$detail$AnyFlag & nzchar(dash$detail$MissingMetrics)))
  out_dir <- tempfile("dashboard-html-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  export_mfrm_bundle(facet_dashboard_fixture$fit, diagnostics = dx, facet = "Rater",
    output_dir = out_dir, prefix = "review", include = c("dashboard", "html"),
    acknowledge_sensitive = TRUE)
  html <- paste(readLines(file.path(out_dir, "review_bundle.html")), collapse = "\n")
  detail <- strsplit(html, "<h2>facet_dashboard_detail</h2>", fixed = TRUE)[[1]][2]
  detail <- strsplit(detail, "</table>", fixed = TRUE)[[1]][1]
  expect_true(all(vapply(dash$detail$Level, grepl, logical(1), x = detail, fixed = TRUE)))
  expect_match(detail, "MissingMetrics", fixed = TRUE)
  expect_match(detail, "Outfit", fixed = TRUE)
  expect_match(html, "IncompleteLevels", fixed = TRUE)
  expect_false(grepl("dashboard-internal-test|ReadinessContractVersion|N[.]x", html))
  csv <- read.csv(file.path(out_dir, "review_facet_dashboard_detail.csv"))
  expect_true(all(csv$ReadinessContractVersion == "dashboard-internal-test"))
  expect_identical(csv$MissingMetrics, dash$detail$MissingMetrics)
})

test_that("rater feedback preserves a restricted design after saving", {
  fit <- facet_dashboard_fixture$fit
  fit$data_review$status$Status[fit$data_review$status$Domain == "Design"] <-
    "review_population_assumption_linked"
  dx <- facet_dashboard_fixture$diagnostics
  dash <- facet_quality_dashboard(fit, diagnostics = dx)
  severity <- plot_rater_severity_profile(fit, diagnostics = dx, draw = FALSE)
  expect_identical(dash$interpretation_status, "review_only")
  expect_identical(severity$data$interpretation_status, "review_only")
  expect_match(severity$data$title, "REVIEW ONLY", fixed = TRUE)
  expect_true(any(grepl("review_population_assumption_linked", dash$notes, fixed = TRUE)))
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  saveRDS(dash, saved)
  restored <- readRDS(saved)
  expect_output(print(summary(restored)), "review_population_assumption_linked")
  p <- plot_facet_quality_dashboard(restored, draw = FALSE)
  expect_match(p$data$title, "REVIEW ONLY", fixed = TRUE)
  expect_identical(p$data$fit_readiness, severity$data$fit_readiness)
})
