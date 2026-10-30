test_that("partial slope bounds survive print, tables, plots and saved reports", {
  source <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))
  ci <- source$profiles$Task
  ci[1, 1] <- .4
  attr(ci, "diagnostics")$CI_Lower <- .4
  attr(ci, "profile")$endpoints$Bound[1] <- .4
  attr(ci, "profile")$endpoints$Status[1] <- "computed"
  before <- serialize(ci, NULL)
  local_mocked_bindings(fit_mfrm = function(...) stop("No refit"),
    compute_mml_parameter_covariance = function(...) stop("No information calculation"),
    mfrm_gpcm_profile_search = function(...) stop("No profile search"), .package = "mfrmr")
  text <- paste(capture.output(print(ci)), collapse = " ")
  expect_match(text, "Incomplete")
  expect_match(text, "NA is not infinity")
  tables <- mfrm_gpcm_inference_tables(ci)
  expect_identical(tables$intervals$IntervalStatus, "Incomplete")
  expect_identical(tables$intervals$CIUse, "unavailable")
  expect_equal(tables$intervals$CI_Lower, .4)
  expect_true(is.na(tables$intervals$CI_Upper))
  expect_identical(apa_table(ci)$table$IntervalStatus, "Incomplete")
  p <- plot(ci, type = "interval", draw = FALSE)
  expect_null(p$labels$title); expect_null(p$labels$caption)
  expect_identical(plot_data(p, "table")$Display, "Incomplete")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$data[[1]]), 0)
  expect_equal(built$data[[2]]$x, .4)
  expect_match(paste(plot_data(p, "notes")$Text, collapse = " "), "neither an infinite")
  result <- mfrm_results(source$fit, intervals = list(partial = ci), include = "fit", compute = "never")
  report <- mfrm_report(result)
  expect_identical(report$tables$gpcm_partial_intervals$IntervalStatus, "Incomplete")
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
  saveRDS(list(result = result, plot = p), file)
  saved <- readRDS(file)
  expect_identical(mfrm_report(saved$result)$tables, report$tables)
  expect_identical(plot_data(saved$plot, "table"), plot_data(p, "table"))
  expect_identical(serialize(ci, NULL), before)
})

test_that("computed bounds, eligibility and unboundedness remain distinct", {
  tab <- data.frame(CI_Lower = c(.4, .4, NA, NA, 0, .4),
    CI_Upper = c(1, NA, 1, NA, Inf, 1), CIEligible = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE))
  expect_identical(mfrm_slope_interval_display(tab),
    c("Available", "Incomplete", "Incomplete", "Unavailable", "Unbounded", "Review"))
})

test_that("bias plots retain conditional search limits and do not replace failed profiles", {
  tab <- data.frame(Pair = c("r1 | t1", "r2 | t1"),
    FacetA = "Rater", FacetB = "Task", FacetA_Level = c("r1", "r2"), FacetB_Level = "t1",
    ObsExpAverage = c(.1, -.1), BiasSize = c(.2, -.2), SE = .1,
    t = c(2, -2), Prob = .05, AbsT = 2, Flag = TRUE,
    TFlag = TRUE, BiasFlag = FALSE, PFlag = TRUE, ObservedCount = 10,
    ProfileCILower = c(0, -5), ProfileCIUpper = c(.4, .1), ProfileCILevel = .95,
    ProfileCIStatus = c("ok", "limited by search range"))
  bundle <- list(ranked_table = tab, scatter_data = tab,
    summary = data.frame(InteractionFacets = "Rater x Task"),
    thresholds = list(abs_t_warn = 2, abs_bias_warn = .5, p_max = .05))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  lines <- list(); original <- graphics::segments
  local_mocked_bindings(segments = function(..., lty = 1) {
    lines[[length(lines) + 1L]] <<- lty
    original(..., lty = lty)
  }, .package = "mfrmr")
  p <- plot_bias_interaction(bundle, plot = "scatter", show_ci = TRUE)
  expect_equal(p$data$scatter_data$CI_Lower, c(0,-5))
  expect_identical(p$data$scatter_data$CI_Status, tab$ProfileCIStatus)
  expect_identical(p$data$scatter_data$CI_Complete, c(TRUE,FALSE))
  expect_true(any(vapply(lines, identical, TRUE, c(1,2))))
  expect_match(paste(plot_data(p, "notes")$Text, collapse = " "), "calibration uncertainty is not propagated")
  expect_match(paste(plot_data(as_ggplot(p), "notes")$Text, collapse = " "), "search limits are not")
  bundle$ranked_table$ProfileCILower <- bundle$ranked_table$ProfileCIUpper <- NA_real_
  bundle$scatter_data <- bundle$ranked_table
  q <- plot_bias_interaction(bundle, plot = "ranked", show_ci = TRUE, draw = FALSE)
  expect_true(all(is.na(q$data$ranked_table$CI_Lower)))
  expect_true(all(q$data$ranked_table$CI_Method == "conditional profile likelihood"))
  expect_false(q$data$plot_settings$CIAvailable)
  expect_false(any(q$data$ranked_table$CI_Complete))
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
  saveRDS(p, file)
  expect_identical(plot_data(readRDS(file), "notes"), plot_data(p, "notes"))
})

test_that("diagonal heterogeneity metadata follows facet precision tables", {
  measures <- data.frame(Facet = "Rater", Estimate = c(-1, 0, 1), SE = .2,
    ModelSE = .2, RealSE = .3)
  q <- calc_facets_chisq(measures)
  expect_equal(q$FixedChiSq, 50)
  expect_false(q$SupportsFormalInference)
  expect_identical(q$InferenceTier, "screening")
  precision <- build_facet_precision_summary(measures, q)
  expect_true(all(precision$StatisticBasis == q$StatisticBasis))
  expect_equal(precision$FixedChiSq, rep(50, nrow(precision)))
})
