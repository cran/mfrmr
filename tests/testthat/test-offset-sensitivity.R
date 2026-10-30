offset_chain <- function(diffs = c(0, 0.4, 0.8, 2), se = rep(NA_real_, length(diffs)), threshold = 0.5) {
  info <- mfrmr:::compute_equating_offset(diffs, se, se, threshold)
  structure(list(
    config = list(waves = c("A", "B"), method = "screened_common_element_alignment",
      drift_threshold = threshold, min_common_per_facet = 5L),
    links = data.frame(Link = 1L, FromID = 1L, ToID = 2L, From = "A", To = "B",
      Offset = info$offset, N_Common = length(diffs), N_Retained = info$n_retained),
    element_detail = data.frame(LinkID = rep(1L, length(diffs)), Facet = rep("Item", length(diffs)),
      Level = sprintf("I%d", seq_along(diffs)), Est_From = rep(0, length(diffs)), Est_To = diffs,
      SE_From = se, SE_To = se, Retained = info$retained, Flag = abs(info$residual) > threshold)
  ), class = c("mfrm_equating_chain", "list"))
}

test_that("removing an excluded element can change offsets through rescreening", {
  x <- offset_chain()
  before <- x
  d <- plot(x, type = "offset_sensitivity", draw = FALSE)$data$data
  expect_identical(x, before)
  expect_equal(d$baseline_links$Offset, 0.6)
  expect_equal(d$baseline_elements$Retained, c(FALSE, TRUE, TRUE, FALSE))
  id <- d$removal$AnchorId[d$removal$Level == "I4"]
  after <- d$links[d$links$RemovedAnchorId == id, ]
  expect_equal(after$Offset, 0.4)
  expect_equal(after$Change, -0.2)
  expect_equal(after$N_Retained, 3L)
  change <- d$retention_changes[d$retention_changes$RemovedAnchorId == id, ]
  expect_equal(change$Level, "I1")
  expect_false(change$Retained_Before)
  expect_true(change$Retained_After)
  expect_equal(d$removal$RetentionChanges[d$removal$AnchorId == id], 1L)
  expect_equal(d$cumulative$Change[d$cumulative$RemovedAnchorId == id], c(0, -0.2))
})

test_that("weighted offsets preserve contributing rows and expose unweighted fallback", {
  x <- offset_chain(c(0, 1, 2), c(1, 2, NA_real_), Inf)
  d <- plot(x, type = "offset_sensitivity", draw = FALSE)$data$data
  expect_equal(d$baseline_links$Offset, 0.2) # weights 1/2, 1/8; missing SE contributes zero
  expect_equal(d$baseline_links$N_Retained, 3L)
  expect_equal(d$baseline_links$N_Contributing, 2L)
  expect_equal(d$baseline_elements$Contributing, c(TRUE, TRUE, FALSE))
  expect_equal(d$links$Offset, c(1, 0, 0.2))
  x <- offset_chain(c(0, 2), c(1, NA_real_), Inf)
  d <- plot(x, type = "offset_sensitivity", draw = FALSE)$data$data
  expect_equal(d$links$Offset, c(2, 0))
  expect_equal(d$links$Offset_Method, c("unweighted", "inverse_variance"))
  x <- offset_chain(c(0, 2), threshold = 0.1)
  d <- plot(x, type = "offset_sensitivity", draw = FALSE)$data$data
  expect_equal(d$baseline_links$Offset, 1)
  expect_true(d$baseline_links$ScreeningFallback)
  expect_false(any(d$links$ScreeningFallback))
})

test_that("missing links propagate and are not reported as zero changes", {
  x <- offset_chain()
  x$config$waves <- c("A", "B", "C")
  x$links <- rbind(x$links, data.frame(Link = 2L, FromID = 2L, ToID = 3L,
    From = "B", To = "C", Offset = 0.25, N_Common = 1L, N_Retained = 1L))
  x$element_detail <- rbind(x$element_detail, data.frame(LinkID = 2L, Facet = "Item",
    Level = "I1", Est_From = 0, Est_To = 0.25, SE_From = NA_real_, SE_To = NA_real_, Retained = TRUE, Flag = FALSE))
  d <- plot(x, type = "offset_sensitivity", draw = FALSE)$data$data
  expect_equal(d$baseline_cumulative$Offset, c(0, 0.6, 0.85))
  id <- d$removal$AnchorId[d$removal$Level == "I1"]
  expect_equal(d$cumulative$Offset[d$cumulative$RemovedAnchorId == id], c(0, 0.8, NA))
  expect_equal(d$links$Status[d$links$RemovedAnchorId == id], c("computed", "no_common"))
  expect_false(d$removal$CompleteComparison[d$removal$AnchorId == id])
  expect_equal(d$removal$NewlyUnavailableLinks[d$removal$AnchorId == id], 1L)
  support <- d$common_by_facet[d$common_by_facet$RemovedAnchorId == id & d$common_by_facet$LinkID == "2", ]
  expect_equal(support$Facet, "Item")
  expect_equal(support$N_Common, 0L)
  expect_false(support$LinkSupportAdequate)
  d <- plot(offset_chain(0.5), type = "offset_sensitivity", draw = FALSE)$data$data
  expect_true(is.na(d$removal$MaxAbsCumulativeChange))
  expect_true(is.na(d$removal$MaxAbsLinkChange))
  expect_equal(d$cumulative$Change, c(0, NA))
})

test_that("nonfinite inputs and empty candidate sets retain failure status", {
  d <- plot(offset_chain(c(NA_real_, Inf)), type = "offset_sensitivity", draw = FALSE)$data$data
  expect_equal(d$baseline_links$Status, "no_finite_differences")
  expect_true(all(is.na(d$removal$MaxAbsCumulativeChange)))
  d <- plot(offset_chain(c(0, 1), c(1e200, 1e200)), type = "offset_sensitivity", draw = FALSE)$data$data
  expect_equal(d$baseline_links$Status, "computed")
  expect_equal(d$baseline_links$Offset, .5)
  expect_true(all(d$links$Status == "computed"))
  expect_equal(d$links$Offset, c(1, 0))
  d <- plot(offset_chain(numeric()), type = "offset_sensitivity", draw = FALSE)$data$data
  expect_equal(d$baseline_links$Status, "no_common")
  expect_equal(nrow(d$removal), 0L)
})

test_that("missing or inconsistent source metadata fails explicitly", {
  x <- offset_chain()
  x$element_detail$SE_From <- NULL
  expect_error(plot(x, type = "offset_sensitivity", draw = FALSE), "numeric Est_From")
  x <- offset_chain()
  x$config$drift_threshold <- NULL
  expect_error(plot(x, type = "offset_sensitivity", draw = FALSE), "recorded screened-linking")
  x <- offset_chain()
  x$links$Offset <- 10
  expect_error(plot(x, type = "offset_sensitivity", draw = FALSE), "baseline does not match")
  x <- offset_chain()
  x$config$waves <- c("A", "B", "C")
  expect_error(plot(x, type = "offset_sensitivity", draw = FALSE), "ordered adjacent chain")
})

test_that("sensitivity drawing retains notes, missing values and caller state", {
  grDevices::pdf(NULL, width = 5, height = 4)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (x in list(offset_chain(), offset_chain(0.5), offset_chain(numeric()))) {
    p <- plot(x, type = "offset_sensitivity", draw = FALSE,
      preset = "monochrome", show_title = FALSE, show_notes = FALSE)
    old <- graphics::par(c("mar", "xpd", "bg"))
    expect_no_warning(q <- plot(x, type = "offset_sensitivity",
      preset = "monochrome", show_title = FALSE, show_notes = FALSE))
    expect_identical(q, p)
    expect_identical(graphics::par(c("mar", "xpd", "bg")), old)
    expect_true(any(grepl("No source model is refitted", p$data$notes$Text)))
    expect_identical(mfrmr:::.mfrm_plot_display_title(p$data), "")
    rgb <- grDevices::col2rgb(p$data$styles$Color)
    expect_equal(rgb[1, ], rgb[2, ])
    expect_equal(rgb[2, ], rgb[3, ])
  }
})
