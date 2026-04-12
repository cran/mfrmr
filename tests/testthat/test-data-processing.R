# test-data-processing.R
# Tests for the data processing layer -- the most important bug hotspot.
# Covers: facet ID conversion, factor vs character equivalence, label
# permutation invariance, missing data, score recoding, anchors, weights,
# and edge cases.

# ---------------------------------------------------------------------------
# 1.1  Facet ID conversion accuracy
# ---------------------------------------------------------------------------

test_that("character facet labels are preserved in fitted output", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
    "Score", method = "JML", maxit = 30))

  # Rater labels

  rater_levels <- sort(unique(d$Rater))
  fitted_rater <- fit$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::pull(Level) |>
    sort()
  expect_equal(fitted_rater, rater_levels)

  # Task labels
  task_levels <- sort(unique(d$Task))
  fitted_task <- fit$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::pull(Level) |>
    sort()
  expect_equal(fitted_task, task_levels)

  # Criterion labels
  crit_levels <- sort(unique(d$Criterion))
  fitted_crit <- fit$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::pull(Level) |>
    sort()
  expect_equal(fitted_crit, crit_levels)

  # Person labels preserved
  expect_setequal(fit$facets$person$Person, unique(d$Person))
})

test_that("facet level ordering follows alphabetical (factor default)", {

  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
    "Score", method = "JML", maxit = 30))

  rater_tbl <- fit$facets$others |>
    dplyr::filter(Facet == "Rater")
  expect_equal(rater_tbl$Level, sort(unique(d$Rater)))
})

# ---------------------------------------------------------------------------
# 1.2  Factor vs character equivalence
# ---------------------------------------------------------------------------

test_that("factor vs character facet columns produce identical estimates", {
  d <- mfrmr:::sample_mfrm_data(seed = 100)
  d_factor <- d |>
    dplyr::mutate(Rater = factor(Rater), Task = factor(Task))

  fit_char <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))
  fit_factor <- suppressWarnings(fit_mfrm(d_factor, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))

  expect_equal(fit_char$summary$LogLik, fit_factor$summary$LogLik, tolerance = 1e-6)

  for (facet in c("Rater", "Task", "Criterion")) {
    est_c <- fit_char$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::arrange(Level) |>
      dplyr::pull(Estimate)
    est_f <- fit_factor$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::arrange(Level) |>
      dplyr::pull(Estimate)
    expect_equal(est_c, est_f, tolerance = 1e-8,
      label = paste("factor vs char for", facet))
  }
})

test_that("factor with explicit level ordering matches character", {
  d <- mfrmr:::sample_mfrm_data(seed = 55)
  d_fac <- d |>
    dplyr::mutate(Rater = factor(Rater, levels = c("R3", "R1", "R2")))

  fit_chr <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))
  fit_fac <- suppressWarnings(fit_mfrm(d_fac, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))

  expect_equal(fit_chr$summary$LogLik, fit_fac$summary$LogLik, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 1.3  Integer ID vs string ID equivalence
# ---------------------------------------------------------------------------

test_that("integer facet IDs and string IDs produce equivalent estimates", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  rater_map <- c(R1 = "1", R2 = "2", R3 = "3")
  d_int <- d |> dplyr::mutate(Rater = rater_map[Rater])

  fit_str <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))
  fit_int <- suppressWarnings(fit_mfrm(d_int, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))

  expect_equal(fit_str$summary$LogLik, fit_int$summary$LogLik, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 1.4  Label permutation invariance (CRITICAL)
# ---------------------------------------------------------------------------

test_that("shuffling rater labels preserves estimates after remapping", {
  d <- mfrmr:::sample_mfrm_data(seed = 200)
  fit_orig <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  shuffle <- c(R1 = "RaterC", R2 = "RaterA", R3 = "RaterB")
  d_shuf <- d |> dplyr::mutate(Rater = shuffle[Rater])
  fit_shuf <- suppressWarnings(fit_mfrm(d_shuf, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  rev_map <- c(RaterC = "R1", RaterA = "R2", RaterB = "R3")
  est_orig <- fit_orig$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  est_remapped <- fit_shuf$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::mutate(OrigLevel = rev_map[Level]) |>
    dplyr::arrange(OrigLevel) |>
    dplyr::pull(Estimate)

  expect_equal(unname(est_orig), unname(est_remapped), tolerance = 0.01)
})

test_that("shuffling task labels preserves estimates after remapping", {
  d <- mfrmr:::sample_mfrm_data(seed = 200)
  fit_orig <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  shuffle <- c(T1 = "TaskD", T2 = "TaskB", T3 = "TaskA", T4 = "TaskC")
  d_shuf <- d |> dplyr::mutate(Task = shuffle[Task])
  fit_shuf <- suppressWarnings(fit_mfrm(d_shuf, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  rev_map <- c(TaskD = "T1", TaskB = "T2", TaskA = "T3", TaskC = "T4")
  est_orig <- fit_orig$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  est_remapped <- fit_shuf$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::mutate(OrigLevel = rev_map[Level]) |>
    dplyr::arrange(OrigLevel) |>
    dplyr::pull(Estimate)

  expect_equal(unname(est_orig), unname(est_remapped), tolerance = 0.01)
})

test_that("shuffling criterion labels preserves estimates after remapping", {
  d <- mfrmr:::sample_mfrm_data(seed = 200)
  fit_orig <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  shuffle <- c(C1 = "CritZ", C2 = "CritX", C3 = "CritY")
  d_shuf <- d |> dplyr::mutate(Criterion = shuffle[Criterion])
  fit_shuf <- suppressWarnings(fit_mfrm(d_shuf, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  rev_map <- c(CritZ = "C1", CritX = "C2", CritY = "C3")
  est_orig <- fit_orig$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  est_remapped <- fit_shuf$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::mutate(OrigLevel = rev_map[Level]) |>
    dplyr::arrange(OrigLevel) |>
    dplyr::pull(Estimate)

  expect_equal(unname(est_orig), unname(est_remapped), tolerance = 0.01)
})

# ---------------------------------------------------------------------------
# 1.5  Missing data handling
# ---------------------------------------------------------------------------

test_that("NA in person column: rows dropped, fit succeeds", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Person[1:5] <- NA
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
    "Score", method = "JML", maxit = 30))

  expect_s3_class(fit, "mfrm_fit")
  expect_lt(nrow(fit$prep$data), nrow(d))
})

test_that("NA in score column: rows dropped, fit succeeds", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score[10:20] <- NA
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
    "Score", method = "JML", maxit = 30))

  expect_s3_class(fit, "mfrm_fit")
  expect_lt(nrow(fit$prep$data), nrow(d))
})

test_that("5% random NA insertion: fit still succeeds, prep data rows < original", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  n_orig <- nrow(d)
  set.seed(999)
  na_idx <- sample(seq_len(n_orig), size = ceiling(0.05 * n_orig))
  d$Score[na_idx] <- NA
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
    "Score", method = "JML", maxit = 30))

  expect_s3_class(fit, "mfrm_fit")
  expect_lt(nrow(fit$prep$data), n_orig)
})

test_that("all NA after filtering raises explicit error", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score <- NA
  expect_error(
    suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
      "Score", method = "JML", maxit = 30)),
    "No valid observations"
  )
})

# ---------------------------------------------------------------------------
# 1.6  Score category handling
# ---------------------------------------------------------------------------

test_that("unused intermediate category recoded to contiguous (keep_original=FALSE)", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  # Only keep scores 1, 3, 5 (drop 2 and 4)
  d_gap <- d |> dplyr::filter(Score %in% c(1, 3, 5))
  fit <- suppressWarnings(fit_mfrm(d_gap, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30, keep_original = FALSE))

  expect_s3_class(fit, "mfrm_fit")
  # After recoding, the data should have contiguous scores
  observed_scores <- sort(unique(fit$prep$data$Score))
  expected_contiguous <- seq(min(observed_scores), max(observed_scores))
  expect_equal(observed_scores, expected_contiguous)
})

test_that("non-consecutive score recoding is surfaced before estimation", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d_gap <- d |> dplyr::filter(Score %in% c(1, 3, 5))

  warning_messages <- character(0)
  prep <- withCallingHandlers(
    mfrmr:::prepare_mfrm_data(
      d_gap,
      "Person",
      c("Rater", "Task", "Criterion"),
      "Score",
      keep_original = FALSE
    ),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("were recoded internally to a contiguous scale", warning_messages, fixed = TRUE)))
  expect_equal(prep$score_map$OriginalScore, c(1, 3, 5))
  expect_equal(prep$score_map$InternalScore, c(1, 2, 3))
})

test_that("keep_original=TRUE preserves original score codes", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  # Only keep scores 1, 3, 5
  d_gap <- d |> dplyr::filter(Score %in% c(1, 3, 5))
  fit <- suppressWarnings(fit_mfrm(d_gap, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30, keep_original = TRUE))

  expect_s3_class(fit, "mfrm_fit")
  observed_scores <- sort(unique(fit$prep$data$Score))
  expect_true(all(observed_scores %in% c(1, 3, 5)))
})

test_that("explicit rating_min/rating_max override auto-detection", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit_auto <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30))
  fit_explicit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30,
    rating_min = 1, rating_max = 5))

  # With the same data and explicit range matching the data, results should match
  expect_equal(fit_auto$summary$LogLik, fit_explicit$summary$LogLik, tolerance = 1e-10)
})

test_that("explicit rating range preserves unused boundary categories", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d_gap <- d |> dplyr::filter(Score %in% 2:5)

  prep <- mfrmr:::prepare_mfrm_data(
    d_gap,
    "Person",
    c("Rater", "Task", "Criterion"),
    "Score",
    rating_min = 1,
    rating_max = 5
  )
  fit <- suppressWarnings(fit_mfrm(
    d_gap, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30,
    rating_min = 1,
    rating_max = 5
  ))
  ds <- describe_mfrm_data(
    d_gap, "Person", c("Rater", "Task", "Criterion"), "Score",
    rating_min = 1,
    rating_max = 5
  )
  ds_summary <- summary(ds)
  rs <- rating_scale_table(fit, drop_unused = FALSE)
  rs_drop <- rating_scale_table(fit, drop_unused = TRUE)
  fit_summary <- summary(fit)
  surface <- plot(fit, type = "ccc_surface", draw = FALSE, theta_points = 55)

  expect_equal(prep$rating_min, 1L)
  expect_equal(prep$rating_max, 5L)
  expect_equal(prep$score_map$OriginalScore, 1:5)
  expect_equal(prep$score_map$InternalScore, 1:5)
  expect_equal(prep$unused_score_categories, 1L)
  expect_equal(fit$config$n_cat, 5L)
  expect_equal(fit$prep$unused_score_categories, 1L)
  expect_true(1 %in% ds$score_distribution$Score)
  expect_equal(ds$score_distribution$RawN[ds$score_distribution$Score == 1], 0L)
  expect_match(paste(ds_summary$notes, collapse = " "), "zero-count")
  expect_equal(fit_summary$settings_overview$UnusedScoreCategories[1], "1")
  expect_equal(fit_summary$settings_overview$UnusedScoreCategoryType[1], "boundary")
  expect_match(paste(fit_summary$key_warnings, collapse = " "), "Unused boundary score")
  expect_true(1 %in% rs$category_table$Category)
  expect_equal(rs$category_table$Count[rs$category_table$Category == 1], 0)
  expect_false(1 %in% rs_drop$category_table$Category)
  expect_equal(rs_drop$summary$UnusedScoreCategories[1], "1")
  expect_true(any(grepl("Unused boundary score", rs_drop$caveats$Message, fixed = TRUE)))
  expect_true(all(c("ZeroCount", "UnusedCategoryType", "WeaklyIdentified", "CategoryCaveat") %in% names(rs$category_table)))
  expect_true(all(c("LowerCategory", "UpperCategory", "WeaklyIdentified", "ThresholdCaveat") %in% names(rs$threshold_table)))
  expect_true(isTRUE(rs$category_table$ZeroCount[rs$category_table$Category == 1]))
  expect_equal(rs$category_table$UnusedCategoryType[rs$category_table$Category == 1], "boundary")
  expect_true(any(rs$threshold_table$WeaklyIdentified, na.rm = TRUE))
  expect_true(any(grepl("zero-count boundary", rs$threshold_table$ThresholdCaveat, fixed = TRUE)))
  expect_true(any(grepl("Unused boundary score", rs$caveats$Message, fixed = TRUE)))
  expect_true(rs$category_table$Percent[rs$category_table$Category == 2] > 0)
  expect_gt(sum(rs$category_table$ExpectedPercent, na.rm = TRUE), 99)
  expect_lt(sum(rs$category_table$ExpectedPercent, na.rm = TRUE), 101)
  expect_s3_class(surface, "mfrm_plot_data")
  expect_true("1" %in% surface$data$surface$Category)
  expect_false("6" %in% surface$data$surface$Category)
  expect_true(all(sort(unique(surface$data$surface$Category)) == as.character(1:5)))
  expect_true("category_support" %in% names(surface$data))
  expect_identical(
    surface$data$category_support$ZeroObserved[surface$data$category_support$Category == "1"][1],
    TRUE
  )
  expect_match(
    surface$data$interpretation_guide$Guidance[
      surface$data$interpretation_guide$Topic == "Zero-frequency categories"
    ][1],
    "Retained categories with zero observed responses: 1",
    fixed = TRUE
  )
})

test_that("explicit rating range still collapses internal gaps unless original codes are requested", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d_gap <- d |> dplyr::filter(Score %in% c(1, 3, 5))

  prep <- suppressWarnings(mfrmr:::prepare_mfrm_data(
    d_gap,
    "Person",
    c("Rater", "Task", "Criterion"),
    "Score",
    rating_min = 1,
    rating_max = 5
  ))
  prep_keep <- mfrmr:::prepare_mfrm_data(
    d_gap,
    "Person",
    c("Rater", "Task", "Criterion"),
    "Score",
    rating_min = 1,
    rating_max = 5,
    keep_original = TRUE
  )

  expect_equal(prep$rating_min, 1L)
  expect_equal(prep$rating_max, 3L)
  expect_equal(prep$score_map$OriginalScore, c(1L, 3L, 5L))
  expect_equal(prep$score_map$InternalScore, 1:3)
  expect_equal(prep_keep$rating_min, 1L)
  expect_equal(prep_keep$rating_max, 5L)
  expect_equal(prep_keep$unused_score_categories, c(2L, 4L))

  fit_keep <- suppressWarnings(fit_mfrm(
    d_gap, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30,
    rating_min = 1,
    rating_max = 5,
    keep_original = TRUE
  ))
  fit_keep_summary <- summary(fit_keep)
  rs_keep <- rating_scale_table(fit_keep, drop_unused = FALSE)
  cs_keep <- category_structure_report(fit_keep)
  rs_keep_summary <- summary(rs_keep)
  cs_keep_summary <- summary(cs_keep)
  surface_keep <- plot(fit_keep, type = "ccc_surface", draw = FALSE, theta_points = 55)
  expect_equal(fit_keep_summary$settings_overview$UnusedScoreCategories[1], "2, 4")
  expect_equal(fit_keep_summary$settings_overview$UnusedScoreCategoryType[1], "internal")
  expect_match(paste(fit_keep_summary$key_warnings, collapse = " "), "Unused intermediate score")
  expect_true(all(rs_keep$category_table$WeaklyIdentified[rs_keep$category_table$Category %in% c(2, 4)]))
  expect_true(any(grepl("zero-count intermediate", rs_keep$threshold_table$ThresholdCaveat, fixed = TRUE)))
  expect_true(any(grepl("Unused intermediate score", rs_keep$caveats$Message, fixed = TRUE)))
  expect_true("caveats" %in% names(rs_keep_summary))
  expect_true(any(grepl("Unused intermediate score", rs_keep_summary$caveats$Message, fixed = TRUE)))
  expect_true("caveats" %in% names(cs_keep_summary))
  expect_true(any(grepl("Unused intermediate score", cs_keep_summary$caveats$Message, fixed = TRUE)))
  expect_true("WeaklyIdentified" %in% names(cs_keep$median_thresholds))
  expect_true(any(cs_keep$median_thresholds$WeaklyIdentified, na.rm = TRUE))
  expect_s3_class(surface_keep, "mfrm_plot_data")
  expect_true(all(as.character(1:5) %in% surface_keep$data$surface$Category))
  expect_true(all(c("2", "4") %in% surface_keep$data$surface$Category))
  expect_identical(
    surface_keep$data$category_support$ZeroObserved[surface_keep$data$category_support$Category == "2"][1],
    TRUE
  )
  expect_match(
    surface_keep$data$interpretation_guide$Guidance[
      surface_keep$data$interpretation_guide$Topic == "Zero-frequency categories"
    ][1],
    "2, 4",
    fixed = TRUE
  )
})

test_that("single intermediate gap is recoded by default and retained when requested", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d_gap <- d |> dplyr::filter(Score %in% c(1, 2, 4, 5))

  prep <- suppressWarnings(mfrmr:::prepare_mfrm_data(
    d_gap,
    "Person",
    c("Rater", "Task", "Criterion"),
    "Score",
    rating_min = 1,
    rating_max = 5
  ))
  prep_keep <- mfrmr:::prepare_mfrm_data(
    d_gap,
    "Person",
    c("Rater", "Task", "Criterion"),
    "Score",
    rating_min = 1,
    rating_max = 5,
    keep_original = TRUE
  )

  expect_equal(prep$rating_min, 1L)
  expect_equal(prep$rating_max, 4L)
  expect_equal(prep$score_map$OriginalScore, c(1L, 2L, 4L, 5L))
  expect_equal(prep$score_map$InternalScore, 1:4)
  expect_length(prep$unused_score_categories, 0L)
  expect_equal(prep_keep$rating_min, 1L)
  expect_equal(prep_keep$rating_max, 5L)
  expect_equal(prep_keep$unused_score_categories, 3L)

  ds <- suppressWarnings(describe_mfrm_data(
    d_gap, "Person", c("Rater", "Task", "Criterion"), "Score",
    rating_min = 1,
    rating_max = 5
  ))
  ds_keep <- describe_mfrm_data(
    d_gap, "Person", c("Rater", "Task", "Criterion"), "Score",
    rating_min = 1,
    rating_max = 5,
    keep_original = TRUE
  )
  ds_summary <- summary(ds)
  ds_keep_summary <- summary(ds_keep)

  expect_true(any(ds_summary$caveats$Condition == "score_categories_recoded"))
  expect_false(any(ds_summary$caveats$Condition == "zero_count_intermediate_score_category"))
  expect_true(3 %in% ds_keep$score_distribution$Score)
  expect_equal(ds_keep$score_distribution$RawN[ds_keep$score_distribution$Score == 3], 0L)
  expect_true(any(ds_keep_summary$caveats$Condition == "zero_count_intermediate_score_category"))
  expect_true(any(ds_keep_summary$caveats$CategoryType == "internal"))
  expect_true(any(grepl("Unused intermediate score", ds_keep_summary$caveats$Message, fixed = TRUE)))
})

test_that("score-support caveats cover upper-boundary and mixed gaps", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d_upper <- d |> dplyr::filter(Score %in% 1:4)
  d_mixed <- d |> dplyr::filter(Score %in% c(2, 4))

  prep_upper <- mfrmr:::prepare_mfrm_data(
    d_upper,
    "Person",
    c("Rater", "Task", "Criterion"),
    "Score",
    rating_min = 1,
    rating_max = 5
  )
  ds_upper <- describe_mfrm_data(
    d_upper, "Person", c("Rater", "Task", "Criterion"), "Score",
    rating_min = 1,
    rating_max = 5
  )
  ds_upper_summary <- summary(ds_upper)

  expect_equal(prep_upper$unused_score_categories, 5L)
  expect_true(5 %in% ds_upper$score_distribution$Score)
  expect_equal(ds_upper$score_distribution$RawN[ds_upper$score_distribution$Score == 5], 0L)
  expect_true(any(ds_upper_summary$caveats$Condition == "zero_count_boundary_score_category"))
  expect_true(any(grepl("5", ds_upper_summary$caveats$Categories, fixed = TRUE)))

  prep_mixed <- mfrmr:::prepare_mfrm_data(
    d_mixed,
    "Person",
    c("Rater", "Task", "Criterion"),
    "Score",
    rating_min = 1,
    rating_max = 5,
    keep_original = TRUE
  )
  ds_mixed <- describe_mfrm_data(
    d_mixed, "Person", c("Rater", "Task", "Criterion"), "Score",
    rating_min = 1,
    rating_max = 5,
    keep_original = TRUE
  )
  ds_mixed_summary <- summary(ds_mixed)

  expect_equal(prep_mixed$unused_score_categories, c(1L, 3L, 5L))
  expect_true(any(ds_mixed_summary$caveats$Condition == "zero_count_boundary_score_category"))
  expect_true(any(ds_mixed_summary$caveats$Condition == "zero_count_intermediate_score_category"))
  expect_true(any(ds_mixed_summary$caveats$CategoryType == "boundary"))
  expect_true(any(ds_mixed_summary$caveats$CategoryType == "internal"))
})

test_that("explicit rating range rejects observed categories outside support", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  expect_error(
    mfrmr:::prepare_mfrm_data(
      d,
      "Person",
      c("Rater", "Task", "Criterion"),
      "Score",
      rating_min = 1,
      rating_max = 4
    ),
    "outside the supplied rating range"
  )
})

# ---------------------------------------------------------------------------
# 1.7  Constraint (anchor) settings
# ---------------------------------------------------------------------------

test_that("anchor R2=0 forces Rater R2 estimate to zero", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  anchor_tbl <- data.frame(
    Facet = "Rater", Level = "R2", Anchor = 0,
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 60,
    anchors = anchor_tbl))

  r2_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R2") |>
    dplyr::pull(Estimate)

  expect_equal(unname(r2_est), 0, tolerance = 1e-8)
})

test_that("invalid anchor (non-existent level) with warn policy triggers warning", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  anchor_tbl <- data.frame(
    Facet = "Rater", Level = "R99_NONEXISTENT", Anchor = 0,
    stringsAsFactors = FALSE
  )
  expect_warning(
    fit_mfrm(d, "Person",
      c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = 60,
      anchors = anchor_tbl,
      anchor_policy = "warn"),
    regex = "Anchor audit"
  )
})

test_that("group anchors constrain group mean to target", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  # Group R1 and R3 together with target mean = 0
  group_tbl <- data.frame(
    Facet = c("Rater", "Rater"),
    Level = c("R1", "R3"),
    Group = c("G1", "G1"),
    GroupValue = c(0, 0),
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 60,
    group_anchors = group_tbl))

  r1_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R1") |>
    dplyr::pull(Estimate)
  r3_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R3") |>
    dplyr::pull(Estimate)

  group_mean <- mean(c(r1_est, r3_est))
  expect_equal(group_mean, 0, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# 1.8  Weight handling
# ---------------------------------------------------------------------------

test_that("zero-weight rows are excluded from estimation", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Weight <- 1
  d$Weight[1:50] <- 0

  fit_wt <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30, weight = "Weight"))

  expect_lt(nrow(fit_wt$prep$data), nrow(d))
})

test_that("explicit weight column produces different results from unit weights", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  set.seed(777)
  d$Weight <- sample(c(1, 2, 3), nrow(d), replace = TRUE)

  fit_unit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 40))
  fit_wt <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 40, weight = "Weight"))

  # LogLik should differ when weights are non-uniform
  expect_false(isTRUE(all.equal(fit_unit$summary$LogLik, fit_wt$summary$LogLik,
    tolerance = 1e-4)))
})

# ---------------------------------------------------------------------------
# 1.9  Edge cases
# ---------------------------------------------------------------------------

test_that("duplicate column names in data raise error with 'duplicate'", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d2 <- cbind(d, Rater = d$Rater)
  # This data.frame now has two columns named "Rater"
  expect_error(
    suppressWarnings(fit_mfrm(d2, "Person",
      c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = 20)),
    regex = "[Dd]uplic"
  )
})

test_that("large number of facet levels (20+ persons) does not crash", {
  set.seed(1234)
  n_person <- 25
  d_large <- expand.grid(
    Person = paste0("P", sprintf("%02d", seq_len(n_person))),
    Rater = paste0("R", 1:2),
    Task = paste0("T", 1:2),
    stringsAsFactors = FALSE
  )
  d_large$Score <- sample(1:3, nrow(d_large), replace = TRUE)

  fit <- suppressWarnings(fit_mfrm(d_large, "Person",
    c("Rater", "Task"), "Score", method = "JML", maxit = 30))

  expect_s3_class(fit, "mfrm_fit")
  expect_equal(nrow(fit$facets$person), n_person)
})

test_that("non-numeric score column (character) is coerced correctly", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score <- as.character(d$Score)

  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30))

  expect_s3_class(fit, "mfrm_fit")
  expect_true(is.integer(fit$prep$data$Score) || is.numeric(fit$prep$data$Score))
})

test_that("non-numeric score or weight entries are surfaced before rows are dropped", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score[1] <- "bad-score"
  d$Weight <- 1
  d$Weight[2] <- "bad-weight"

  warning_messages <- character(0)
  prep <- withCallingHandlers(
    mfrmr:::prepare_mfrm_data(d, "Person", c("Rater", "Task", "Criterion"), "Score", weight_col = "Weight"),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("`Score` contained", warning_messages, fixed = TRUE)))
  expect_true(any(grepl("`Weight` contained", warning_messages, fixed = TRUE)))
  expect_equal(nrow(prep$data), nrow(d) - 2L)
})

test_that("fractional score inputs are rejected before integer coercion", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score[1] <- 1.5
  d$Score[2] <- 3.25

  expect_error(
    fit_mfrm(
      d,
      "Person",
      c("Rater", "Task", "Criterion"),
      "Score",
      method = "JML",
      maxit = 10
    ),
    "`Score` must contain ordered integer category codes",
    fixed = TRUE
  )
})

test_that("prepare_mfrm_data returns expected structure", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  prep <- mfrmr:::prepare_mfrm_data(d, "Person",
    c("Rater", "Task", "Criterion"), "Score")

  expect_true(is.list(prep))
  expect_true("data" %in% names(prep))
  expect_true("levels" %in% names(prep))
  expect_true("facet_names" %in% names(prep))
  expect_true("rating_min" %in% names(prep))
  expect_true("rating_max" %in% names(prep))
  expect_equal(prep$facet_names, c("Rater", "Task", "Criterion"))
  expect_true("Person" %in% names(prep$levels))
  expect_true("Rater" %in% names(prep$levels))
})

test_that("duplicate person/facets/score arguments raise error", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  expect_error(
    mfrmr:::prepare_mfrm_data(d, "Person", c("Person", "Task"), "Score"),
    regex = "[Dd]uplic"
  )
})
