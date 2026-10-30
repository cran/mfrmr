gmfrm_response_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      saved <- readRDS(test_path("fixtures", "gmfrm-joint-information.rds"))
      problem <- mfrm_gmfrm_problem(saved$data, 2L, gauss_hermite_normal(31L))
      fit <- mfrm_gmfrm_fit_result(problem, mfrm_gmfrm_em(problem, start = saved$parameters, maxit = 1L))
      persons <- unique(saved$data$Person)[c(1, 20, 160)]
      rows <- which(saved$data$Person %in% persons)
      diagnostics <- mfrm_response_diagnostics(fit, rows = rows, group_by = c("Rater", "Task"), quad_points = 121L)
      cached <<- list(fit = fit, diagnostics = diagnostics, rows = rows, saved = saved)
    }
    cached
  }
})

test_that("two-family response moments integrate the complete Person posterior", {
  x <- gmfrm_response_fixture(); f <- x$fit; d <- x$diagnostics
  expect_true(all(d$rows$Status == "available_conditional"))
  expect_identical(d$rows$InputRow, x$rows)
  expect_identical(lapply(d$source_data, as.character), lapply(f$gmfrm$specification$data, as.character))
  expect_false(d$settings$calibration_uncertainty)
  expect_match(d$settings$target, "Same-data posterior")
  expect_match(d$settings$limitation, "no expectation-one")
  for (at in 1:3) {
    reference <- gmfrm_response_reference(f, x$rows[at])
    expect_equal(unname(d$probabilities[at, ]), unname(reference$probabilities[1, ]), tolerance = 1e-8)
    expect_equal(d$rows$ExpectedScore[at], reference$mean, tolerance = 1e-8)
    expect_equal(d$rows$PredictiveVariance[at], reference$variance, tolerance = 1e-8)
    expect_gt(reference$variance - reference$conditional_variance, .001)
    expect_gt(max(abs(reference$plugin - reference$probabilities)), .001)
  }
  for (i in seq_len(nrow(d$measures))) {
    tab <- d$measures[i, ]
    at <- d$source_data[[tab$Facet]][x$rows] == tab$Level
    expect_equal(tab$Infit, sum(d$rows$SquaredResidual[at]) / sum(d$rows$PredictiveVariance[at]))
    expect_equal(tab$Outfit, mean(d$rows$StandardizedResidual[at]^2))
  }
  selected <- rev(x$rows[c(2, 8, 14)])
  subset <- mfrm_response_diagnostics(f, rows = selected, quad_points = 121L)
  expect_equal(unname(subset$probabilities), unname(d$probabilities[match(selected, x$rows), ]))
  expect_identical(subset$rows$InputRow, selected)
  expect_identical(f, x$fit)
})

test_that("integration failures remain in rows and group denominators", {
  x <- gmfrm_response_fixture()
  bad <- mfrm_response_diagnostics(x$fit, rows = x$rows, quad_points = 7L, group_by = "Rater")
  expect_true(any(bad$rows$Status == "unavailable"))
  expect_true(all(bad$rows$Reason[bad$rows$Status == "unavailable"] != ""))
  expect_true(all(is.na(bad$probabilities[bad$rows$Status == "unavailable", ])))
  expect_equal(sum(bad$measures$Observed), length(x$rows))
  unresolved <- bad$measures$Available < bad$measures$Observed
  expect_true(all(is.na(bad$measures$Infit[unresolved])))
  expect_true(all(is.na(bad$measures$Outfit[unresolved])))
  result <- mfrm_results(x$fit, response_diagnostics = bad)
  overview <- result$tables$response_overview
  expect_equal(overview$Available + overview$Unresolved, nrow(bad$rows))
  expect_equal(overview$Unresolved, sum(bad$rows$Status == "unavailable"))
  expect_equal(overview$Selected + overview$NotIncluded, nrow(bad$source_data))
  expect_identical(result$status$Status[result$status$Section == "response_diagnostics"], overview$Status)
  report <- mfrm_report(result)
  expect_identical(report$first_screen$MainIssue[report$first_screen$Area == "Response residuals"], overview$Detail)
  expect_match(report$markdown, overview$Detail, fixed = TRUE)
  expect_error(mfrm_response_diagnostics(x$fit, rows = c(1, 1)), "distinct original")
  expect_error(mfrm_response_diagnostics(x$fit, group_by = "Score"), "identifier")
  expect_error(mfrm_response_diagnostics(x$fit, quad_points = 6), "Invalid quad_points")
})

test_that("all-unavailable two-family residuals remain unavailable in report overviews", {
  x <- gmfrm_response_fixture()
  local_mocked_bindings(mfrm_gmfrm_response_probabilities = function(...) stop("unresolved integration"),
    .package = "mfrmr")
  d <- mfrm_response_diagnostics(x$fit, rows = x$rows, group_by = "Rater")
  result <- mfrm_results(x$fit, response_diagnostics = d)
  expect_identical(result$status$Status[result$status$Section == "response_diagnostics"], "not_available")
  report <- mfrm_report(result)
  expect_identical(report$first_screen$Status[report$first_screen$Area == "Response residuals"], "unavailable")
  expect_equal(report$tables$response_overview$Unresolved, length(x$rows))
  expect_equal(report$tables$response_overview$Available, 0L)
  expect_identical(report$tables$response_residuals, d$rows)
  # Older saved results obtain the same overview from their retained rows.
  result$tables$response_overview <- NULL
  result$status$Status[result$status$Section == "response_diagnostics"] <- "available"
  before <- serialize(result, NULL)
  expect_identical(summary(result)$status$Status[summary(result)$status$Section == "response_diagnostics"], "not_available")
  expect_identical(mfrm_report(result)$tables$response_overview, report$tables$response_overview)
  viewer <- mfrm_results_viewer_payload(result)
  expect_identical(viewer$summary$status, summary(result)$status)
  expect_identical(viewer$tables$response_overview, report$tables$response_overview)
  expect_identical(serialize(result, NULL), before)
})

test_that("two-family diagnostics retain arbitrary owner names and incomplete assignments", {
  x <- gmfrm_response_fixture(); data <- x$saved$data
  names(data) <- c("Candidate", "観点 名", "Judge-ID", "Rating")
  data[["観点 名"]] <- sub("t", "level:", data[["観点 名"]], fixed = TRUE)
  data[["Judge-ID"]] <- sub("r", "level:", data[["Judge-ID"]], fixed = TRUE)
  p <- mfrm_gmfrm_problem(data, 2L, gauss_hermite_normal(31L),
    slope_facets = c("観点 名", "Judge-ID"), person = "Candidate", score = "Rating")
  f <- mfrm_gmfrm_fit_result(p, mfrm_gmfrm_em(p, start = x$saved$parameters, maxit = 1L))
  d <- mfrm_response_diagnostics(f, rows = x$rows, group_by = "Judge-ID", quad_points = 121L)
  expect_equal(d$probabilities, x$diagnostics$probabilities)
  expect_identical(d$settings$group_by, "Judge-ID")
  expect_true(all(d$measures$Facet == "Judge-ID"))
  expect_identical(names(d$source_data), names(data))
  # A saved genuinely sparse fit: no refitting or fabricated unassigned cells.
  sparse <- readRDS(test_path("fixtures", "gmfrm-small-optimization-residual.rds"))
  p <- mfrm_gmfrm_problem(sparse$data, 2L, gauss_hermite_normal(31L))
  f <- mfrm_gmfrm_fit_result(p, mfrm_gmfrm_em(p, start = sparse$parameters, maxit = 1L))
  expect_lt(nrow(sparse$data), prod(lengths(p$levels)))
  selected <- which(sparse$data$Person == sparse$data$Person[1])
  d <- mfrm_response_diagnostics(f, rows = selected, quad_points = 121L)
  expect_identical(lapply(d$source_data, as.character), lapply(sparse$data, as.character))
  expect_equal(nrow(d$rows), length(selected))
  expect_true(all(d$rows$Status == "available_conditional"))
})

test_that("two-family diagnostics reject stale calibration and unsupported fitting scope", {
  x <- gmfrm_response_fixture(); f <- x$fit
  changes <- list(
    function(z) {z$opt$par[1] <- z$opt$par[1] + .1; z},
    function(z) {z$slopes$Estimate[1] <- z$slopes$Estimate[1] + .1; z},
    function(z) {z$steps$Estimate[1] <- z$steps$Estimate[1] + .1; z},
    function(z) {z$facets$others$Estimate[1] <- NA_real_; z},
    function(z) {z$config$slope_facet <- rev(z$config$slope_facet); z},
    function(z) {z$prep$data$Weight[1] <- 2; z},
    function(z) {z$prep$data$score_k[1] <- 2; z},
    function(z) {z$config$estimation_control$quad_points <- 15; z},
    function(z) {z$population$active <- TRUE; z},
    function(z) {z$opt$convergence <- 1L; z})
  for (change in changes) expect_error(mfrm_response_diagnostics(change(f), rows = 1))
  # Interval eligibility does not decide whether point predictions exist.
  f$slopes$CIEligible <- FALSE
  d <- mfrm_response_diagnostics(f, rows = x$rows, quad_points = 121L)
  expect_equal(d$probabilities, x$diagnostics$probabilities)
  f <- x$fit; f$gmfrm$controls$maxit <- f$gmfrm$controls$maxit + 1L
  expect_error(mfrm_results(f, include = "fit", compute = "never",
    response_diagnostics = x$diagnostics), "matching calibration")
})

test_that("saved GMFRM residuals connect to plots reports and exports without recalculation", {
  x <- gmfrm_response_fixture(); f <- x$fit; d <- x$diagnostics
  path <- tempfile("gmfrm-response-export-"); on.exit(unlink(path, recursive = TRUE), add = TRUE)
  local_mocked_bindings(mfrm_gmfrm_em = function(...) stop("must not refit"),
    mfrm_gpcm_response_evaluator = function(...) stop("must not recompute"),
    mfrm_jml_probability_bundle = function(...) stop("must not integrate"), .package = "mfrmr")
  res <- mfrm_results(f, include = c("fit", "plots"), compute = "never", response_diagnostics = d)
  expect_identical(res$tables$response_measures, d$measures)
  expect_identical(mfrm_results(f, include = "fit", compute = "never", diagnostics = d)$response_diagnostics, d)
  expect_error(mfrm_results(f, include = "fit", compute = "never", diagnostics = d, response_diagnostics = d), "only once")
  for (style in c("paired", "scatter")) {
    p <- plot(res, type = "response_diagnostics", style = style, draw = FALSE)
    expect_identical(plot_data(plot(res, style=style, draw=FALSE)), plot_data(p))
    expect_s3_class(p, "mfrm_plot_data")
    expect_s3_class(as_ggplot(p), "ggplot")
    expect_identical(ggplot2::ggplot_build(as_ggplot(res, style=style))$data,
      ggplot2::ggplot_build(as_ggplot(p))$data)
    expect_match(plot_data(p)$caption, "no reference cutoffs")
  }
  report <- mfrm_report(res)
  expect_identical(report$tables$response_residuals, d$rows)
  expect_match(report$markdown, "Same-data posterior predictive")
  expect_match(report$markdown, "no expectation-one")
  exported <- export_mfrm_results(res, output_dir = path, preset = "starter", acknowledge_sensitive = TRUE)
  expect_equal(nrow(exported$plot_errors), 0L)
  expect_true("plot_response_diagnostics" %in% exported$written_files$Component)
  restored <- readRDS(exported$written_files$Path[exported$written_files$Component == "results_rds"])
  expect_identical(restored$response_diagnostics, d)
  expect_identical(summary(restored$response_diagnostics), summary(d))
  replay <- exported$written_files$Path[exported$written_files$Component == "replay_code"]
  withr::with_dir(path, sys.source(replay, envir = new.env(parent = globalenv())))
})

test_that("two-family saved response reports reopen in a fresh process", {
  skip_if_not_installed("callr")
  x <- gmfrm_response_fixture()
  res <- mfrm_results(x$fit, include = c("fit", "plots"), compute = "never", response_diagnostics = x$diagnostics)
  path <- tempfile(fileext = ".rds"); withr::defer(unlink(path)); saveRDS(res, path)
  worker <- function(root, path) {
    if (file.exists(file.path(root, "R", "api-results.R"))) pkgload::load_all(root, quiet = TRUE, compile = FALSE) else
      library("mfrmr", lib.loc = dirname(root), character.only = TRUE)
    testthat::local_mocked_bindings(fit_mfrm = function(...) stop("unexpected refit"),
      mfrm_jml_probability_bundle = function(...) stop("unexpected integration"), .package = "mfrmr")
    res <- readRDS(path)
    list(summary = summary(res$response_diagnostics),
      report = mfrmr::mfrm_report(res)$tables$response_measures,
      caption = mfrmr::plot_data(plot(res, draw = FALSE))$caption)
  }
  environment(worker) <- baseenv()
  out <- callr::r(worker, args = list(normalizePath(find.package("mfrmr")), path), libpath = .libPaths())
  expect_identical(out$summary, summary(x$diagnostics))
  expect_identical(out$report, x$diagnostics$measures)
  expect_match(out$caption, "no reference cutoffs")
})
