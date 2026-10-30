test_that("plot guide calls resolve without changing the route catalogue", {
  guide <- mfrmr_output_guide("plots")
  expect_identical(names(guide), c("Question", "InputClass", "PlotCall", "PlotName",
    "DataComponent", "GGPlot", "GGPlotCall", "Notes", "ResultFunction", "NextStep"))
  expect_false(anyDuplicated(guide$Question) > 0)
  expect_false(anyNA(guide[setdiff(names(guide), "GGPlotCall")]))
  expect_identical(is.na(guide$GGPlotCall), guide$GGPlot == "unavailable")
  expect_setequal(guide$GGPlot, c("dedicated", "native", "generic", "unavailable"))
  # These result classes inherit the registered plot.mfrm_bundle() method.
  bundle_classes <- c("mfrm_rating_scale", "mfrm_category_curves",
    "mfrm_bias_interaction", "mfrm_halo_network", "mfrm_rater_network")
  for (i in seq_len(nrow(guide))) {
    call <- str2lang(guide$PlotCall[i])
    fun <- if (identical(call[[1]], as.name("plot"))) {
      getS3method("plot", if (guide$InputClass[i] %in% bundle_classes)
        "mfrm_bundle" else guide$InputClass[i])
    } else get(as.character(call[[1]]), mode = "function")
    expect_true(is.function(fun))
    expect_true(is.function(get(guide$ResultFunction[i], envir = asNamespace("mfrmr"))))
    expect_identical(call$draw, FALSE)
    expect_true("..." %in% names(formals(fun)) ||
      all(names(as.list(call)[-1L])[-1L] %in% names(formals(fun))))
  }
  expect_identical(mfrmr_output_guide("beginner"), mfrmr_output_guide("public"))
  expect_identical(names(mfrmr_output_guide()), names(mfrmr_output_guide("feedback")))
  locations <- subset(mfrmr_interval_guide("gpcm"), Route == "Fixed-facet model/sandwich intervals")
  expect_match(locations$GPCMStatus, "supported_with_caveat one-family")
  expect_match(locations$GPCMStatus, "experimental two-family")
  expect_match(locations$InterpretationBoundary, "regular MML assumptions")
})

# Exercise the calls printed in the guide, including extraction and conversion.
check_plot_guide_route <- function(row, x) {
  device <- grDevices::dev.cur()
  p <- eval(str2lang(row$PlotCall))
  saved <- mfrmr:::as_mfrm_plot_data_object(p)
  expect_identical(saved$name, row$PlotName)
  expect_true(row$DataComponent %in% names(plot_data(p)))
  expect_identical(plot_data(p, component = row$DataComponent), saved$data[[row$DataComponent]])
  if (row$GGPlot != "unavailable") {
    g <- eval(str2lang(row$GGPlotCall))
    expect_s3_class(g, "ggplot")
    expect_no_error(ggplot2::ggplot_build(g))
  } else expect_error(as_ggplot(p), "conversion.*(unavailable|not available)")
  expect_identical(grDevices::dev.cur(), device)
  invisible(p)
}

test_that("interval guide calls reuse CML and full MML comparisons on their declared scales", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  guide <- mfrmr_interval_guide("fit")
  routes <- c("CML conditional-information intervals", "MML interaction-contrast intervals")
  rows <- guide[match(routes, guide$Route), ]
  expect_false(anyNA(rows))
  expect_false(any(routes %in% mfrmr_interval_guide("gpcm")$Route))
  plot_row <- subset(mfrmr_output_guide("plots"), InputClass == "mfrm_facet_intervals")
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected refit"),
    mfrmr_cml_evaluate = function(...) stop("Unexpected conditional information"))

  fit <- readRDS(test_path("fixtures", "cml-comparison.rds"))
  x <- eval(str2lang(rows$PrimaryHelper[1]))
  expect_identical(x$covariance, fit$covariance)
  for (column in c("Estimate", "SE", "Lower", "Upper"))
    expect_identical(x$table[[column]], fit$comparisons[[column]])
  expect_identical(x$table$CIEligible, c(FALSE, TRUE))
  p <- check_plot_guide_route(plot_row, x)
  expect_identical(plot_data(p, "table"), summary(x))
  expect_identical(apa_table(x)$table$CIEligible, x$table$CIEligible)
  saved <- tempfile(fileext = ".rds")
  withr::defer(unlink(saved))
  saveRDS(x, saved)
  expect_identical(readRDS(saved), x)

  source <- readRDS(test_path("fixtures", "mml-comparison-intervals.rds"))$RSM
  fit <- source$fit
  review <- source$review
  info <- source$intervals
  local_mocked_bindings(mfrm_native_location_inference = function(...) list(
    covariance = info$parameter_covariance, information = list(status = "ok"),
    checks = info$checks, numerical_checks = info$numerical_checks,
    check = list(eligible = TRUE, review = "Saved joint information qualification")))
  x <- suppressWarnings(eval(str2lang(rows$PrimaryHelper[2])))
  expect_identical(x$settings$scale, "standardized_population_sd")
  expect_equal(x$settings$population_sd, sqrt(fit$population$sigma2))
  # These intervals are recomputed from a saved cross-platform fixture;
  # standardized-scale multiplication can differ by a few floating-point ulps.
  for (column in c("Estimate", "SE", "Lower", "Upper"))
    expect_equal(x$table[[column]], info$table[[column]], tolerance = 1e-12)
  p <- check_plot_guide_route(plot_row, x)
  expect_identical(plot_data(p, "table")$IdentificationBasis, x$table$IdentificationBasis)
  result <- x
  # The reporting call and named plot type are the ones printed in DisplayRoute.
  bundle <- mfrm_results(fit, intervals = list(comparisons = result), compute = "never")
  expect_identical(plot_data(bundle, type = "facet_comparisons", component = "table"), x$table)
  reported <- mfrm_report(bundle)$tables$facet_comparisons_intervals
  expect_identical(reported[names(x$table)], x$table)
  expect_identical(reported$ConfidenceLevel, rep("95%", nrow(x$table)))
})

test_that("core plot guide routes produce the advertised payloads", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  x <- make_toy_fit(maxit = 20)
  guide <- subset(mfrmr_output_guide("plots"), InputClass == "mfrm_fit")
  for (i in seq_len(nrow(guide))) {
    .mfrmr_muffle_expected_warnings(check_plot_guide_route(guide[i, ], x),
      "^Review-only display:")
  }
})

test_that("native and generic guide routes describe different rendering contracts", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  guide <- mfrmr_output_guide("plots")
  x <- structure(matrix(0, 2, 2), class = c("mfrm_slope_intervals", "matrix", "array"),
    target = "slope", method = "model", level = .95, simultaneous = "none",
    diagnostics = data.frame(SlopeFacet = c("C1", "C2"), Estimate = c(1, 1.2),
      CI_Lower = c(.8, .9), CI_Upper = c(1.2, 1.5), CIEligible = TRUE))
  p <- check_plot_guide_route(subset(guide, InputClass == "mfrm_slope_intervals"), x)
  expect_equal(plot_data(p, "table")$Estimate, attr(x, "diagnostics")$Estimate)
  x <- structure(list(table = data.frame(Theta = c(-1, 0, 1), Rater = "R1",
    Category = 1, Estimate = c(.2, .5, .8), Lower = c(.1, .4, .7),
    Upper = c(.3, .6, .9), CIEligible = TRUE, InputRow = 1:3),
    settings = list(facets = "Rater", type = "probability", method = "model",
      level = .95, simultaneous = "none")), class = "mfrm_curve_intervals")
  p <- check_plot_guide_route(subset(guide, InputClass == "mfrm_curve_intervals"), x)
  expect_equal(plot_data(p, "table")$Estimate, x$table$Estimate)
  x <- structure(list(tif = data.frame(Theta = c(-1, 0, 1),
    Information = c(2, 4, 3), SE = 1 / sqrt(c(2, 4, 3)))), class = "mfrm_information")
  row <- subset(guide, InputClass == "mfrm_information")
  p <- check_plot_guide_route(row, x)
  g <- eval(str2lang(row$GGPlotCall))
  points <- ggplot2::ggplot_build(g)$data[[1]]
  expect_equal(points$x, x$tif$Theta)
  expect_equal(points$y, x$tif$Information)
  expect_s3_class(g$layers[[1]]$geom, "GeomPoint")
})

test_that("selected extended-model and D-study routes use saved results", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  guide <- mfrmr_output_guide("plots")
  # Fixed saved-score shape: no expensive estimation is needed to test display.
  x <- structure(list(table = data.frame(Person = c("P1", "P2"),
    Estimate = c(.5, -.5), Lower = c(-.2, -1), Upper = c(1.2, .2),
    Status = "available_conditional"), settings = list(level = .95)),
    class = "mfrm_random_rater_scores")
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected refit"))
  check_plot_guide_route(subset(guide, InputClass == "mfrm_random_rater_scores"), x)
  tasks <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
  g <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
  x <- mfrm_d_study(g, data.frame(Tasks = c(3, 6, 12)), c(V = -1, W = 1))
  rows <- subset(guide, InputClass == "mfrm_multivariate_d_study")
  for (i in seq_len(nrow(rows))) {
    p <- check_plot_guide_route(rows[i, ], x)
    expect_equal(plot_data(p)$weights, c(V = -1, W = 1))
  }
})

test_that("every unavailable guide route refuses even a convertible component", {
  rows <- subset(mfrmr_output_guide("plots"), GGPlot == "unavailable")
  # Valid generic-converter shapes must not bypass a plot-specific refusal.
  for (i in seq_len(nrow(rows))) {
    call <- str2lang(rows$PlotCall[i])
    p <- mfrmr:::new_mfrm_plot_data(rows$PlotName[i], list(
      plot = call$type %||% NULL, metric = call$metric %||% NULL,
      table = data.frame(Label = c("A", "B"), Value = c(1, 2)),
      matrix = matrix(1:4, 2)))
    message <- if (rows$InputClass[i] == "mfrm_residual_network")
      "supports its type = 'matrix' view only" else "conversion.*(unavailable|not available)"
    expect_error(as_ggplot(p), message)
    expect_error(as_ggplot(p, component = "table"), message)
    expect_error(as_ggplot(p, component = "matrix"), message)
  }
})


test_that("subset coverage guidance distinguishes a matrix from the full base figure", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  fit <- make_toy_fit(maxit = 20)
  x <- subset_connectivity_report(fit, diagnostics = make_toy_diagnostics(fit))
  row <- subset(mfrmr_output_guide("plots"), PlotName == "subset_connectivity")
  p <- check_plot_guide_route(row, x)
  expect_identical(row$GGPlot, "generic")
  expect_match(row$Notes, "omits the observation-share panel", fixed = TRUE)
  matrix <- plot_data(p, "matrix")
  g <- eval(str2lang(row$GGPlotCall))
  expect_equal(sort(g$data$Value), sort(as.numeric(matrix)))
})

test_that("all three PCA guide routes use dedicated axis-preserving conversion", {
  skip_if_not_installed("ggplot2")
  x <- mfrm_pca(mfrm_features(data.frame(ID = letters[1:6],
    Years = c(1, 2, 5, 7, 10, 12), Hours = c(4, 8, 6, 10, 16, 14)),
    "ID", c("Years", "Hours")))
  rows <- subset(mfrmr_output_guide("plots"), InputClass == "mfrm_pca")
  expect_equal(nrow(rows), 3)
  expect_true(all(rows$GGPlot == "dedicated"))
  for (i in seq_len(nrow(rows))) check_plot_guide_route(rows[i, ], x)
})

test_that("the hierarchy guide converts the complete tree and extracts its source", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cluster")
  x <- mfrm_cluster_hierarchical(mfrm_features(data.frame(ID = letters[1:5],
    Years = c(1, 2, 5, 10, 12)), "ID", "Years"), 2)
  row <- subset(mfrmr_output_guide("plots"), PlotName == "cluster_dendrogram")
  expect_identical(row$GGPlot, "dedicated")
  check_plot_guide_route(row, x)
})

test_that("the co-membership guide preserves the complete fixed-scale view", {
  skip_if_not_installed("ggplot2")
  ids <- c("b", "a", "excluded")
  mat <- matrix(c(1, .5, NA, .5, 1, NA, NA, NA, NA), 3, dimnames = list(ids, ids))
  x <- structure(list(co_membership = mat, settings = list(imputations = 2L)),
    class = "mfrm_imputed_clusters")
  row <- subset(mfrmr_output_guide("plots"), PlotName == "cluster_co_membership")
  expect_identical(row$GGPlot, "dedicated")
  check_plot_guide_route(row, x)
})

test_that("partition guide routes preserve silhouettes and both feature profile types", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cluster")
  x <- mfrm_cluster(mfrm_features(data.frame(ID = letters[1:6],
    Years = c(1, 2, 3, 10, 12, 18), Specialty = rep(c("A", "B"), each = 3)),
    "ID", c("Years", "Specialty")), 2)
  rows <- subset(mfrmr_output_guide("plots"), InputClass == "mfrm_clusters")
  for (i in seq_len(nrow(rows))) check_plot_guide_route(rows[i, ], x)
  row <- subset(rows, PlotName == "cluster_profile")
  row$PlotCall <- sub("Years", "Specialty", row$PlotCall)
  check_plot_guide_route(row, x)
})

test_that("pooled interval guide preserves source intervals through dedicated conversion", {
  skip_if_not_installed("ggplot2")
  x <- structure(list(table = data.frame(Target = "R02 minus R01", Estimate = .4,
    Lower = -.1, Upper = .9, DF = 8, Status = "pooled"),
    settings = list(facet = "Rater", imputations = 5L, ci_level = .95)), class = "mfrm_pooled")
  row <- subset(mfrmr_output_guide("plots"), PlotName == "pooled_facet_intervals")
  expect_identical(row$GGPlot, "dedicated")
  check_plot_guide_route(row, x)
})


test_that("network intervals enter the interval guide without changing visual templates", {
  guide <- mfrmr_interval_guide()
  row <- guide[guide$Route == "Binary rater-network coefficient intervals", ]
  expect_equal(nrow(row), 1L)
  expect_match(row$PrimaryHelper, "confint(model", fixed = TRUE)
  expect_match(row$InterpretationBoundary, "estimated-edge uncertainty", fixed = TRUE)
  expect_no_error(template <- visual_reporting_template())
  expect_false("Binary rater-network coefficient intervals" %in% unlist(template))
})
