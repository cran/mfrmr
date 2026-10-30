category_owner_fixture <- function() {
  readRDS(test_path("fixtures", "native-location-pcm-regression.rds"))$fit
}

category_owner_diagnostics <- function(fit) {
  steps <- as.data.frame(fit$steps)
  steps$SE <- seq_len(nrow(steps)) / 100
  steps$CI_Lower <- steps$Estimate - 1.96 * steps$SE
  steps$CI_Upper <- steps$Estimate + 1.96 * steps$SE
  steps$SE_Status <- "regularized"
  steps$CIEligible <- FALSE
  steps$CIUse <- "review_before_reporting"
  structure(list(calculation_version = 2L, obs = compute_obs_table(fit), facet_names = fit$config$facet_names,
    fit_readiness = as.data.frame(mfrmr_get_readiness_record(fit)$fit),
    parameter_uncertainty = list(steps = steps,
      source = mfrm_native_location_identity(fit))), class = "mfrm_diagnostics")
}

test_that("ordinary threshold ownership is scalar and retains documented defaults", {
  resolve <- function(step, model = "PCM", method = "MML", slope = NULL)
    resolve_step_and_slope_facets(model, step, slope, c("Rater", "Task"), method)
  for (bad in list(c("Task", "Rater"), character(), NA_character_, "", list("Task"))) {
    expect_error(resolve(bad), "exactly one")
  }
  expect_error(resolve("Person"), "not among")
  expect_identical(resolve(NULL)$step_facet, "Task")
  expect_identical(resolve(NULL)$step_facet_source, "inferred_from_name")
  expect_identical(resolve("Rater", "GPCM", slope = "Task")$slope_facet, "Task")
  expect_error(resolve("Rater", "GPCM", "JML", "Task"), "require MML")
  expect_warning(fallback <- resolve_step_and_slope_facets("PCM", NULL, NULL, c("Judge", "Domain")), "first")
  expect_identical(fallback$step_facet, "Judge")
})

test_that("owner counts distinguish exposure, weight, unavailable and unobserved categories", {
  fit <- list(config = list(model = "PCM", method = "JML", step_facet = "Judge",
    facet_levels = list(Judge = c("A", "B", "C"))),
    prep = list(rating_min = 1, rating_max = 3,
      score_map = data.frame(InternalScore = 1:3, OriginalScore = c(2, 4, 8))))
  obs <- data.frame(Judge = c("A", "A", "A", "B", "B"),
    Observed = c(1, 1, 3, 2, 3), Weight = c(1, 2, 1, 0, 2))
  out <- category_owner_usage(fit, obs)
  expect_equal(out$Count, c(3, 0, 1, 0, 0, 2, 0, 0, 0))
  expect_equal(out$RawCount, c(2, 0, 1, 0, 1, 1, 0, 0, 0))
  expect_equal(out$Percent[1:6], c(75, 0, 25, 0, 0, 100))
  expect_true(all(is.na(out$Percent[7:9])))
  expect_equal(out$UsageStatus, c("observed", "not_observed", "observed",
    "not_observed", "zero_weight", "observed", rep("no_exposure", 3)))
  expect_identical(out$CategoryLabel, rep(c("2", "4", "8"), 3))
  obs$Weight[1] <- NA_real_
  missing <- category_owner_usage(fit, obs)
  expect_true(all(is.na(missing$Count[1:3])))
  expect_equal(missing$Count[4:9], out$Count[4:9])
  obs$Judge[1] <- NA_character_
  expect_true(all(is.na(category_owner_usage(fit, obs)$Count)))
})

test_that("stored step uncertainty is joined by keys and retains review status", {
  f <- category_owner_fixture()
  d <- category_owner_diagnostics(f)
  original <- d$parameter_uncertainty$steps
  d$parameter_uncertainty$steps <- original[nrow(original):1, ]
  out <- rating_scale_table(f, d)
  expect_equal(out$threshold_table$SE, original$SE)
  expect_true(all(out$threshold_table$SE_Status == "regularized"))
  expect_true(all(out$threshold_table$CIUse == "review_before_reporting"))
  expect_true(all(out$threshold_table$UncertaintySourceStatus == "source_verified"))
  d$parameter_uncertainty$steps <- original[-2, ]
  expect_true(is.na(rating_scale_table(f, d)$threshold_table$SE[2]))
  d$parameter_uncertainty$steps <- rbind(original, original[1, ])
  expect_error(rating_scale_table(f, d), "unique matching")
  d$parameter_uncertainty$steps <- original
  d$parameter_uncertainty$steps$Estimate[1] <- original$Estimate[1] + .1
  expect_error(rating_scale_table(f, d), "stored step estimates")
  d <- category_owner_diagnostics(f)
  d$parameter_uncertainty$source$source$parameters[1] <- 99
  expect_error(category_structure_report(f, d), "fit/diagnostics mismatch")
  d$parameter_uncertainty$source <- NULL
  d$parameter_uncertainty$steps$CIEligible <- TRUE
  d$parameter_uncertainty$steps$CIUse <- "primary_reporting"
  expect_error(rating_scale_table(f, d), "do not record their fitted calibration")
})

test_that("saved PCM and GPCM category workflows do not request fitting or covariance", {
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
    simulate_mfrm_data = function(...) stop("No generation"),
    compute_mml_parameter_covariance = function(...) stop("No covariance calculation"),
    diagnose_mfrm = function(...) stop("No implicit diagnostics"))
  for (f in list(category_owner_fixture(),
    readRDS(test_path("fixtures", "native-location-rsm-zero-boundary.rds"))$fit, readRDS(test_path("fixtures", "native-location-fixed-gpcm.rds")),
    readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))$fit)) {
    rating <- rating_scale_table(f)
    owner <- rating$owner_category_table
    expect_equal(aggregate(Count ~ Category, owner, sum)$Count, rating$category_table$Count)
    expect_true(all(is.na(rating$threshold_table$SE)))
    structure <- category_structure_report(f)
    expect_equal(structure$threshold_table, rating$threshold_table)
    expect_identical(structure$owner_category_table, owner)
    curves <- category_curves_report(f)
    expect_true(all(curves$curve_profiles$StepOwner == rating$model_context$StepOwner))
    for (x in list(rating, structure, curves)) {
      payload <- plot(x, draw = FALSE)$data
      expect_identical(payload$model_context, x$model_context)
      path <- tempfile(); on.exit(unlink(path), add = TRUE)
      saveRDS(x, path)
      expect_identical(readRDS(path), x)
    }
  }
  two <- readRDS(test_path("fixtures", "gmfrm-explicit-start.rds"))$fit
  expect_error(rating_scale_table(two), "two|Two")
  expect_error(category_structure_report(two), "two|Two")
  expect_error(category_curves_report(two), "two|Two")
})

test_that("a binary centered step has no free sampling contrast", {
  f <- category_owner_fixture()
  f$config$n_cat <- 2L
  f$prep$rating_max <- 1
  f$steps <- f$steps[f$steps$Step == "Step_1", ]
  f$steps$Estimate <- 0
  usage <- category_owner_usage(f, data.frame(Criterion = "1", Observed = 0, Weight = 1))
  out <- category_threshold_table(f, NULL, usage)
  expect_true(out$coverage$NotApplicable)
  expect_true(all(out$table$SE_Status == "fixed_by_step_centering"))
  expect_true(all(out$table$CIUse == "not_applicable"))
})

test_that("threshold ladders follow step owners rather than slope profiles", {
  local_mocked_bindings(fit_mfrm = function(...) stop("No fitting"),
    build_step_curve_spec = function(...) stop("No curve calculation"))
  fit <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))$fit
  p <- plot(fit, type = "thresholds", draw = FALSE)
  expect_equal(nrow(p$data$data), nrow(fit$steps))
  expect_equal(length(unique(p$data$data$Group)), length(unique(fit$steps$StepFacet)))
  expect_true(all(p$data$data$StepOwner == fit$config$step_facet))
  key <- match(paste(p$data$data$Group, p$data$data$Step), paste(fit$steps$StepFacet, fit$steps$Step))
  expect_equal(p$data$data$Threshold, fit$steps$Estimate[key])
  expect_match(paste(plot_data(p, "notes")$Text, collapse = " "), "not an automatic category-merging rule")
  expect_error(as_ggplot(p, component = "data"), "unavailable")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  fit$steps$Estimate[2] <- NA_real_
  expect_silent(p <- plot(fit, type = "thresholds"))
  expect_true(anyNA(p$data$data$Threshold))
  fit$steps$Estimate[] <- NA_real_
  expect_silent(plot(fit, type = "thresholds"))
  fit$steps$Estimate[] <- 0
  expect_silent(plot(fit, type = "thresholds"))
})

test_that("category conversion preserves counts, step values and missing gaps", {
  skip_if_not_installed("ggplot2")
  x <- as_mfrm_bundle(list(category_table = data.frame(Category = 0:2,
    Count = c(12, 0, NA), ExpectedCount = c(10, 2, NA)),
    threshold_table = data.frame(StepFacet = rep(c("A", "B"), each = 3),
      Step = rep(paste0("Step_", 1:3), 2), Estimate = c(-1, NA, 1, 1, 0, -1)),
    notes = "Stored values; step intervals are not curve intervals."), "mfrm_rating_scale")
  x$category_usage <- summarize_category_usage(x$category_table)
  x$threshold_coverage <- summarize_threshold_order(calc_step_order(x$threshold_table))
  p <- plot(x, draw = FALSE, palette = c(counts = "orange", expected = "purple"))
  g <- as_ggplot(p)
  built <- ggplot2::ggplot_build(g)
  expect_equal(built$data[[1]]$y, c(12, 0))
  expect_equal(built$data[[2]]$y[1:2], c(10, 2))
  expect_true(all(built$data[[1]]$fill == "orange"))
  expect_true(all(built$data[[2]]$colour == "purple"))
  expect_equal(built$data[[4]]$label, "Unavailable")
  expect_null(g$labels$title)
  expect_null(g$labels$caption)
  expect_identical(plot_data(g), p$data)
  expect_error(as_ggplot(p, component = "threshold_table"), "selected category view")
  p <- plot(x, type = "thresholds", draw = FALSE, palette = c(step_line = "purple"))
  g <- as_ggplot(p, component = "threshold_table")
  built <- ggplot2::ggplot_build(g)
  expect_equal(built$data[[1]]$x, c(4, 5))
  expect_equal(built$data[[1]]$xend, c(5, 6))
  expect_equal(built$data[[2]]$y, x$threshold_table$Estimate)
  expect_true(all(built$data[[2]]$colour == "purple"))
  expect_identical(plot_data(g), p$data)
  expect_error(as_ggplot(p, component = "category_table"), "selected category view")
  old <- p; old$data$palette <- old$data$label_angle <- NULL
  expect_s3_class(as_ggplot(old), "ggplot")
  x$category_table$Count[] <- NA_real_
  x$category_table$ExpectedCount[] <- NA_real_
  expect_s3_class(as_ggplot(plot(x, draw = FALSE)), "ggplot")

  seen <- character(); title <- graphics::title
  local_mocked_bindings(title = function(main = NULL, ...) {
    seen <<- c(seen, as.character(main)); title(main = main, ...)
  }, .package = "mfrmr")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  local_graphics_headings(function(main) seen <<- c(seen, as.character(main)))
  plot(x, type = "thresholds")
  expect_false(any(nzchar(seen)))
})

test_that("category curve and boundary conversion cannot mix targets through a component", {
  for (name in c("category_curves", "category_structure")) {
    for (view in c("overview", "ccc", "cumulative", "mode_boundaries", "mean_halfscore")) {
      p <- new_mfrm_plot_data(name, list(plot = view, plot_long = data.frame(Theta = 0, Value = 1)))
      expect_error(as_ggplot(p), "unavailable")
      expect_error(as_ggplot(p, component = "plot_long"), "unavailable")
    }
  }
})
