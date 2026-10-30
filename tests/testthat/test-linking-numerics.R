linking_numeric_fits <- function(diffs, se = rep(1, length(diffs))) {
  first <- data.frame(Facet = "Rater", Level = paste0("R", seq_along(diffs)),
    Anchor = rep(0, length(diffs)))
  second <- first
  second$Anchor <- diffs
  list(A = structure(list(estimates = first, se = transform(first[c("Facet", "Level")],
      SE = se), config = list(model = "RSM")), class = "mfrm_fit"),
    B = structure(list(estimates = second, se = transform(second[c("Facet", "Level")],
      SE = se), config = list(model = "RSM")), class = "mfrm_fit"))
}

mock_linking_inputs <- function(env = parent.frame()) {
  local_mocked_bindings(
    make_anchor_table = function(fit, ...) fit$estimates,
    measure_se_table = function(fit, ...) fit$se,
    .package = "mfrmr", .env = env)
}

test_that("link offsets are invariant to a common change of SE units", {
  f <- mfrmr:::compute_equating_offset
  d <- c(-1, .5, 2)
  se <- c(1, 2, 4)
  expected <- weighted.mean(d, 1 / (2 * se^2))
  for (scale in c(1e-200, 1, 1e200)) {
    x <- f(d, se * scale, se * scale, Inf)
    expect_equal(x$offset, expected, tolerance = 1e-12)
    expect_true(all(x$contributing))
    expect_identical(x$status, "computed")
  }
  # Avoid overflow of weighted sums and loss of tiny weights before
  # multiplication by a large source difference.
  x <- f(c(1e308, -1e308), c(1, 2), c(1, 2), Inf)
  expect_equal(x$offset / 1e308, .6, tolerance = 1e-12)
  x <- f(c(0, 1e300), c(1e-100, 1e100), c(1e-100, 1e100), Inf)
  expect_equal(x$offset / 1e-100, 1, tolerance = 1e-12)
  expect_true(all(x$contributing))
  x <- f(rep(1e308, 3), c(1, 2, 3), c(1, 2, 3), 0)
  expect_identical(x$offset, 1e308)
  expect_identical(x$residual, rep(0, 3))
})

test_that("unavailable link spread differs from zero spread", {
  mock_linking_inputs()
  x <- suppressWarnings(build_equating_chain(linking_numeric_fits(1)))
  expect_equal(x$links$Offset, 1)
  expect_true(is.na(x$links$Offset_SD))
  expect_identical(x$links$OffsetSDStatus, "insufficient_elements")
  x <- suppressWarnings(build_equating_chain(linking_numeric_fits(c(NA_real_, Inf))))
  expect_true(is.na(x$links$Offset))
  expect_true(is.na(x$links$Offset_SD))
  expect_true(is.na(x$links$Max_Residual))
  expect_identical(x$links$OffsetStatus, "no_finite_differences")
  expect_true(all(is.na(x$element_detail$Residual)))
  expect_true(is.na(summary(x)$n_flagged))
  expect_equal(summary(x)$n_unavailable, 2L)
  expect_no_error(capture.output(print(summary(x))))
  x <- suppressWarnings(build_equating_chain(linking_numeric_fits(rep(.3, 2))))
  expect_equal(x$links$Offset_SD, 0)
  expect_identical(x$links$OffsetSDStatus, "available")
})

test_that("link builders retain screening fallback and actual contributors", {
  mock_linking_inputs()
  fits <- linking_numeric_fits(rep(c(-2, 2), each = 3))
  chain <- build_equating_chain(fits, drift_threshold = .1)
  drift <- detect_anchor_drift(fits, drift_threshold = .1)
  expect_true(chain$links$ScreeningFallback)
  expect_true(all(drift$drift_table$ScreeningFallback))
  expect_equal(chain$links$N_Contributing, 6L)
  expect_true(all(chain$element_detail$Contributing))
  expect_true(all(chain$element_detail$Flag))
  # This remains a count convention, not a successful invariance test.
  expect_true(chain$links$LinkSupportAdequate)
  fits <- linking_numeric_fits(c(0, 1, 2), c(1, 2, NA_real_))
  chain <- suppressWarnings(build_equating_chain(fits, drift_threshold = Inf))
  drift <- suppressWarnings(detect_anchor_drift(fits, drift_threshold = Inf))
  expect_equal(chain$links$Offset, .2)
  expect_equal(chain$links$N_Retained, 3L)
  expect_equal(chain$links$N_Contributing, 2L)
  expect_identical(chain$element_detail$Contributing, c(TRUE, TRUE, FALSE))
  expect_true(all(drift$drift_table$N_Contributing == 2L))
  expect_identical(drift$drift_table$FlagBasis[drift$drift_table$Level == "R3"],
    "magnitude_only")
})

test_that("drift preserves representable SEs and exposes unavailable offsets", {
  mock_linking_inputs()
  for (scale in c(1e-200, 1e200)) {
    x <- suppressWarnings(detect_anchor_drift(linking_numeric_fits(c(0, 1), rep(scale, 2)),
      drift_threshold = Inf))
    expect_equal(x$drift_table$SE / scale, rep(sqrt(2), 2), tolerance = 1e-12)
    expect_true(all(is.finite(x$drift_table$Drift_SE_Ratio)))
  }
  x <- suppressWarnings(detect_anchor_drift(linking_numeric_fits(c(NA_real_, Inf))))
  expect_true(all(is.na(x$drift_table$Drift)))
  expect_true(all(x$drift_table$OffsetStatus == "no_finite_differences"))
  expect_true(all(x$drift_table$FlagBasis == "unavailable"))
  expect_true(all(is.na(x$summary$N_Flagged)))
  expect_equal(summary(x)$n_unavailable, 2L)
  expect_true(is.na(summary(x)$n_flagged))
  expect_no_error(capture.output(print(summary(x))))
  review <- build_linking_review(drift = x)
  expect_gt(nrow(review$top_linking_risks), 0L)
})

test_that("linking reviews retain support warnings when no drift is flagged", {
  mock_linking_inputs()
  for (n_common in c(2L, 6L)) {
    drift <- suppressWarnings(detect_anchor_drift(linking_numeric_fits(rep(0, n_common))))
    expect_equal(drift$summary$N_Flagged, 0L)
    expect_equal(drift$summary$N_Unavailable, 0L)
    review <- build_linking_review(drift = drift)
    expect_type(review$drift_risks$Signal, "character")
    if (n_common == 2L) {
      expect_identical(review$drift_risks$SourceFamily, "thin_link_support")
      expect_match(review$drift_risks$Signal, "Retained common-element support")
      expect_equal(review$drift_risks$Magnitude, 3)
    } else {
      expect_equal(nrow(review$drift_risks), 0L)
    }
    expect_s3_class(summary(review), "summary.mfrm_linking_review")
  }
})

test_that("unrepresentable spread and cumulative offsets are not zeros or infinities", {
  mock_linking_inputs()
  x <- suppressWarnings(build_equating_chain(
    linking_numeric_fits(c(-1.5e308, 1.5e308)), drift_threshold = Inf))
  expect_equal(x$links$Offset, 0)
  expect_true(is.na(x$links$Offset_SD))
  expect_identical(x$links$OffsetSDStatus, "numerical_failure")
  expect_equal(x$links$Max_Residual, 1.5e308)
  x <- mfrmr:::.linking_cumulative_offsets(c(1e308, 1e308, -1e308))
  expect_true(is.na(x[3]))
  expect_false(any(is.infinite(x)))
  expect_identical(mfrmr:::.linking_cumulative_offsets(c(.2, NA_real_, .3)),
    c(0, .2, NA_real_, NA_real_))
  expect_equal(mfrmr:::.linking_difference_se(c(0, -1, 1e200), c(0, 1, 1e200)),
    c(0, NA_real_, sqrt(2) * 1e200))
})

test_that("drift drawings retain unavailable rows and handle all-zero heatmaps", {
  mock_linking_inputs()
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (values in list(rep(0, 6), rep(NA_real_, 6), c(rep(0, 5), NA_real_))) {
    x <- suppressWarnings(detect_anchor_drift(linking_numeric_fits(values)))
    for (type in c("drift", "heatmap", "forest")) {
      expect_no_warning(p <- plot(x, type = type))
      q <- plot(x, type = type, draw = FALSE)
      expect_identical(p, q)
      tab <- plot_data(p, if (type == "forest") "data" else "table")
      expect_identical(tab$OffsetStatus,
        x$drift_table$OffsetStatus[match(tab$Level, x$drift_table$Level)])
    }
  }
})

test_that("linking inputs fail before extraction instead of recycling controls", {
  mock_linking_inputs()
  fits <- linking_numeric_fits(c(0, 1))
  for (bad in list(-1, NA_real_, c(.1, .2), "small")) {
    expect_error(build_equating_chain(fits, drift_threshold = bad), "drift_threshold")
    expect_error(detect_anchor_drift(fits, drift_threshold = bad), "drift_threshold")
    expect_error(detect_anchor_drift(fits, flag_se_ratio = bad), "flag_se_ratio")
  }
  expect_error(detect_anchor_drift(fits, reference = 1.5), "reference")
  expect_error(detect_anchor_drift(fits, reference = c(1, 2)), "reference")
  expect_error(mfrmr:::compute_equating_offset(c(0, 1), 1, c(1, 1)), "same length")
})

test_that("ggplot linking views use their named statistic and preserve saved results", {
  skip_if_not_installed("ggplot2")
  mock_linking_inputs()
  drift <- suppressWarnings(detect_anchor_drift(linking_numeric_fits(c(0, .2, .7))))
  original <- serialize(drift, NULL)
  for (view in c("drift", "heatmap", "forest")) {
    p <- plot(drift, type = view, draw = FALSE)
    g <- as_ggplot(p)
    expect_identical(attr(g, "mfrmr_plot_data"), p)
    if (view == "heatmap") {
      expect_equal(g$data$Drift, as.vector(p$data$matrix))
      expect_equal(toupper(g$scales$get_scales("fill")$palette(c(0, 1))),
        toupper(p$data$legend$value))
    } else if (view == "forest") {
      expect_equal(g$data$Value, p$data$data$Estimate)
      expect_equal(g$data$CI_Lower, p$data$data$CI_Lower)
      expect_equal(g$data$CI_Upper, p$data$data$CI_Upper)
    } else {
      expect_equal(sort(g$data$Value), sort(drift$drift_table$Drift))
      expect_false(identical(g$data$Value, rep(0, 3)))
    }
    expect_no_warning(ggplot2::ggplot_build(g))
    expect_error(as_ggplot(p, component = "links"), "selected linking view")
  }
  expect_identical(serialize(drift, NULL), original)
  chain <- suppressWarnings(build_equating_chain(linking_numeric_fits(c(.3, .3))))
  p <- plot(chain, type = "chain", draw = FALSE)
  expect_equal(as_ggplot(p)$data$Cumulative_Offset, c(0, .3))
  p <- plot(chain, type = "common_anchors", draw = FALSE)
  expect_equal(as_ggplot(p)$data$Value, 2)
  # Do not bridge an unavailable intermediate wave.
  p <- plot(chain, type = "chain", draw = FALSE)
  p$data$table <- data.frame(Wave = c("A", "B", "C", "D"),
    Cumulative_Offset = c(0, NA_real_, .5, .8))
  g <- as_ggplot(p)
  expect_equal(nrow(g$layers[[1]]$data), 1L)
  expect_equal(g$layers[[1]]$data$Position, 3L)
  expect_true(is.na(g$data$Cumulative_Offset[2]))
  expect_no_warning(ggplot2::ggplot_build(g))
})
