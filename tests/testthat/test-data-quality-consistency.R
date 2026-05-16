make_qc_score_gap_case <- function(keep_original = TRUE, maxit = 8) {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d <- d[d$Score != 3, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = maxit,
      rating_min = 1, rating_max = 5,
      keep_original = keep_original
    )
  )
  dq <- data_quality_report(
    fit,
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    min_category_count = 10
  )
  list(data = d, fit = fit, dq = dq)
}

make_qc_response_pattern_case <- function(maxit = 8) {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score[d$Rater == "R1"] <- 1L
  fit <- suppressWarnings(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = maxit,
      rating_min = 1, rating_max = 5,
      keep_original = TRUE
    )
  )
  dq <- data_quality_report(
    fit,
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score"
  )
  list(data = d, fit = fit, dq = dq)
}

count_flag_area <- function(flags, area, severity = NULL) {
  if (nrow(flags) == 0L) return(0L)
  keep <- flags$Area %in% area
  if (!is.null(severity)) {
    keep <- keep & tolower(as.character(flags$Severity)) %in% severity
  }
  sum(keep, na.rm = TRUE)
}

test_that("data-quality summary counts are recomputable from returned detail tables", {
  case <- make_qc_score_gap_case(keep_original = TRUE)
  dq <- case$dq
  summary_tbl <- as.data.frame(dq$summary, stringsAsFactors = FALSE)
  category_counts <- as.data.frame(dq$category_counts, stringsAsFactors = FALSE)
  usage_summary <- as.data.frame(dq$category_usage_summary, stringsAsFactors = FALSE)
  patterns <- as.data.frame(dq$facet_response_patterns, stringsAsFactors = FALSE)
  caveats <- as.data.frame(dq$caveats, stringsAsFactors = FALSE)

  expect_equal(
    summary_tbl$ZeroCountScoreCategories,
    sum(category_counts$ZeroCount %in% TRUE, na.rm = TRUE)
  )
  expect_equal(
    summary_tbl$IntermediateZeroCountScoreCategories,
    sum(category_counts$UnusedCategoryType == "internal", na.rm = TRUE)
  )
  expect_equal(
    summary_tbl$FacetLevelsWithZeroCategories,
    sum(usage_summary$ZeroCategories > 0, na.rm = TRUE)
  )
  expect_equal(
    summary_tbl$FacetLevelsWithIntermediateZeroCategories,
    sum(usage_summary$IntermediateZeroCategories > 0, na.rm = TRUE)
  )
  expect_equal(
    summary_tbl$FacetLevelsWithSparseCategories,
    sum(usage_summary$SparseCategories > 0, na.rm = TRUE)
  )
  expect_equal(
    summary_tbl$FacetLevelsWithSingleCategoryUse,
    sum(patterns$SingleCategoryUse %in% TRUE, na.rm = TRUE)
  )
  expect_equal(
    summary_tbl$FacetLevelsWithDominantCategoryUse,
    sum(patterns$DominantCategoryUse %in% TRUE, na.rm = TRUE)
  )
  expect_equal(
    summary_tbl$FacetLevelsWithBoundaryOnlyUse,
    sum(patterns$BoundaryOnlyUse %in% TRUE, na.rm = TRUE)
  )
  expect_equal(summary_tbl$ScoreSupportCaveats, nrow(caveats))

  support_expected <- category_counts[
    order(!as.logical(category_counts$ZeroCount), category_counts$Score),
    names(dq$score_support_review),
    drop = FALSE
  ]
  row.names(support_expected) <- NULL
  expect_equal(dq$score_support_review, support_expected)
})

test_that("quality flags and overview summarize the same QC evidence", {
  case <- make_qc_score_gap_case(keep_original = TRUE)
  dq <- case$dq
  flags <- as.data.frame(dq$quality_flags, stringsAsFactors = FALSE)
  overview <- as.data.frame(dq$quality_overview, stringsAsFactors = FALSE)
  summary_tbl <- as.data.frame(dq$summary, stringsAsFactors = FALSE)

  score_support <- overview[overview$Area == "Score support", , drop = FALSE]
  facet_use <- overview[overview$Area == "Facet category use", , drop = FALSE]
  patterns <- overview[overview$Area == "Facet response patterns", , drop = FALSE]

  expect_equal(score_support$Status, "high")
  expect_equal(score_support$Count, summary_tbl$ZeroCountScoreCategories)
  expect_equal(facet_use$Status, "high")
  expect_equal(
    facet_use$Count,
    max(
      summary_tbl$FacetLevelsWithZeroCategories,
      summary_tbl$FacetLevelsWithIntermediateZeroCategories,
      summary_tbl$FacetLevelsWithSparseCategories
    )
  )
  expect_equal(
    patterns$Count,
    max(
      summary_tbl$FacetLevelsWithSingleCategoryUse,
      summary_tbl$FacetLevelsWithDominantCategoryUse,
      summary_tbl$FacetLevelsWithBoundaryOnlyUse
    )
  )

  for (area in overview$Area) {
    row <- overview[overview$Area == area, , drop = FALSE]
    expect_equal(row$QualityFlags, count_flag_area(flags, area))
    expect_equal(row$HighSeverityFlags, count_flag_area(flags, area, "high"))
  }

  p_flags <- plot(dq, type = "quality_flags", draw = FALSE)
  plot_tbl <- as.data.frame(p_flags$data$table, stringsAsFactors = FALSE)
  manual_plot <- flags |>
    dplyr::group_by(.data$Area) |>
    dplyr::summarise(
      Flags = dplyr::n(),
      HighSeverityFlags = sum(tolower(as.character(.data$Severity)) %in% "high", na.rm = TRUE),
      ReviewFlags = sum(tolower(as.character(.data$Severity)) %in% "review", na.rm = TRUE),
      TotalReferencedCount = sum(suppressWarnings(as.numeric(.data$Count)), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$HighSeverityFlags), dplyr::desc(.data$Flags), .data$Area) |>
    as.data.frame(stringsAsFactors = FALSE)
  row.names(manual_plot) <- NULL
  expect_equal(plot_tbl, manual_plot)
  expect_equal(p_flags$data$quality_flags, dq$quality_flags)
  expect_equal(p_flags$data$quality_overview, dq$quality_overview)
})

test_that("facet category usage plot data follows category-usage summary ordering", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score[d$Rater == "R1" & d$Score == 3] <- 2L
  fit <- suppressWarnings(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = 8,
      rating_min = 1, rating_max = 5,
      keep_original = TRUE
    )
  )
  dq <- data_quality_report(
    fit,
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    min_category_count = 10
  )
  usage <- as.data.frame(dq$category_usage_summary, stringsAsFactors = FALSE)
  expected <- usage |>
    dplyr::arrange(
      dplyr::desc(.data$IntermediateZeroCategories),
      dplyr::desc(.data$ZeroCategories),
      dplyr::desc(.data$SparseCategories),
      .data$Facet,
      .data$Level
    ) |>
    dplyr::slice_head(n = 4) |>
    as.data.frame(stringsAsFactors = FALSE)
  plot_obj <- plot(dq, type = "facet_category_usage", top_n = 4, draw = FALSE)

  expect_equal(plot_obj$data$table, expected)
  expect_equal(plot_obj$data$top_n, 4L)
  expect_true(any(plot_obj$data$table$Facet == "Rater" &
                    plot_obj$data$table$Level == "R1" &
                    plot_obj$data$table$IntermediateZeroCategories > 0))
})

test_that("facet response-pattern flags and plot data preserve restricted scoring evidence", {
  case <- make_qc_response_pattern_case()
  dq <- case$dq
  flags <- as.data.frame(dq$quality_flags, stringsAsFactors = FALSE)
  patterns <- as.data.frame(dq$facet_response_patterns, stringsAsFactors = FALSE)
  r1 <- patterns[patterns$Facet == "Rater" & patterns$Level == "R1", , drop = FALSE]

  expect_equal(nrow(r1), 1L)
  expect_true(isTRUE(r1$SingleCategoryUse[1]))
  expect_false(isTRUE(r1$DominantCategoryUse[1]))
  expect_equal(r1$DominantPercent[1], 1, tolerance = 1e-12)
  expect_equal(r1$PatternStatus[1], "high")
  expect_true(any(flags$Flag == "Facet levels use only one score category"))
  expect_equal(
    flags$Count[flags$Flag == "Facet levels use only one score category"],
    sum(patterns$SingleCategoryUse %in% TRUE, na.rm = TRUE)
  )

  expected <- patterns |>
    dplyr::arrange(
      match(.data$PatternStatus, c("high", "review", "ok")),
      dplyr::desc(.data$DominantPercent),
      .data$Facet,
      .data$Level
    ) |>
    dplyr::slice_head(n = 5) |>
    as.data.frame(stringsAsFactors = FALSE)
  plot_obj <- plot(dq, type = "facet_response_patterns", top_n = 5, draw = FALSE)
  expect_equal(plot_obj$data$table, expected)
  expect_true(plot_obj$data$table$PatternStatus[1] %in% "high")
})

test_that("score-map plot data exposes original-label gaps hidden by recoding", {
  case <- make_qc_score_gap_case(keep_original = FALSE)
  dq <- case$dq
  score_map <- as.data.frame(dq$score_map, stringsAsFactors = FALSE)
  caveats <- as.data.frame(dq$caveats, stringsAsFactors = FALSE)
  flags <- as.data.frame(dq$quality_flags, stringsAsFactors = FALSE)

  expect_true(any(caveats$Condition == "original_score_gap_before_recoding"))
  expect_equal(
    caveats$Categories[caveats$Condition == "original_score_gap_before_recoding"],
    "3"
  )
  expect_true(any(flags$Flag == "Original score sequence had gaps before recoding"))

  plot_obj <- plot(dq, type = "score_map", draw = FALSE)
  plot_tbl <- as.data.frame(plot_obj$data$table, stringsAsFactors = FALSE)
  expected <- score_map[, c("OriginalScore", "InternalScore"), drop = FALSE]
  expected$OriginalNumeric <- suppressWarnings(as.numeric(expected$OriginalScore))
  expected$InternalNumeric <- suppressWarnings(as.numeric(expected$InternalScore))
  expected$MappingStatus <- ifelse(
    as.character(expected$OriginalScore) == as.character(expected$InternalScore),
    "identity",
    "recoded"
  )
  expected <- expected[order(expected$OriginalNumeric, as.character(expected$OriginalScore), na.last = TRUE), , drop = FALSE]
  row.names(expected) <- NULL

  expect_equal(plot_tbl, expected)
  expect_true(any(plot_tbl$MappingStatus == "recoded"))
  expect_false(3 %in% plot_tbl$OriginalNumeric)
})
