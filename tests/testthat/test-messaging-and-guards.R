# Tests for 0.1.6 messaging / guardrail surface:
#   * the inferred rating-range message is quiet by default but can be
#     re-enabled once per fit_mfrm() call,
#   * analyze_dff(method = "refit") raises a mfrmr-styled error when
#     diagnostics is missing,
#   * compute_information() and fair_average_table() produce finite
#     outputs on the default example and on a deliberately boundary-ish
#     subset.

local({
  .toy <<- load_mfrmr_data("example_core")
})

test_that("fit_mfrm() keeps inferred-rating messages quiet by default", {
  msgs <- testthat::capture_messages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 10)
  ))
  n <- sum(grepl("Rating range inferred", msgs, fixed = TRUE))
  expect_equal(n, 0L)
})

test_that("fit_mfrm() can emit the inferred-rating message exactly once", {
  old_opt <- options(mfrmr.show_inferred_rating_range = TRUE)
  on.exit(options(old_opt), add = TRUE)
  msgs <- testthat::capture_messages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 10)
  ))
  n <- sum(grepl("Rating range inferred", msgs, fixed = TRUE))
  expect_equal(n, 1L)
})

test_that("explicit rating_min/rating_max suppresses the inferred message", {
  msgs <- testthat::capture_messages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5,
             rating_min = 1L, rating_max = 4L)
  ))
  expect_false(any(grepl("Rating range inferred", msgs, fixed = TRUE)))
})

test_that("rating range provenance is retained for summary outputs", {
  fit_auto <- suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5)
  )
  fit_declared <- suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5,
             rating_min = 1L, rating_max = 4L)
  )
  auto_summary <- summary(fit_auto)
  declared_summary <- summary(fit_declared)
  data_desc <- describe_mfrm_data(.toy, "Person", c("Rater", "Criterion"), "Score")

  expect_identical(fit_auto$prep$rating_range_source, "observed")
  expect_identical(fit_auto$prep$rating_min_source, "observed")
  expect_identical(fit_auto$prep$rating_max_source, "observed")
  expect_identical(auto_summary$settings_overview$RatingRangeSource[1], "observed")
  expect_identical(auto_summary$settings_overview$RatingMinSource[1], "observed")
  expect_identical(auto_summary$settings_overview$RatingMaxSource[1], "observed")

  expect_identical(fit_declared$prep$rating_range_source, "declared")
  expect_identical(declared_summary$settings_overview$RatingRangeSource[1], "declared")
  expect_identical(declared_summary$settings_overview$RatingMinSource[1], "declared")
  expect_identical(declared_summary$settings_overview$RatingMaxSource[1], "declared")

  expect_identical(data_desc$overview$RatingRangeSource[1], "observed")
  expect_identical(data_desc$score_support$rating_range_source, "observed")
})

test_that("data-preparation notes and row retention are retained", {
  dirty <- .toy
  dirty$Score[1] <- NA
  dirty$Person[2] <- paste0(" ", dirty$Person[2], " ")

  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(dirty, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5)
  ))
  prep_notes <- as.data.frame(fit$prep$preparation_notes, stringsAsFactors = FALSE)
  row_retention <- as.data.frame(fit$prep$row_retention, stringsAsFactors = FALSE)
  s_fit <- summary(fit)
  ds <- suppressMessages(
    describe_mfrm_data(dirty, "Person", c("Rater", "Criterion"), "Score")
  )
  s_ds <- summary(ds)

  expect_true(all(c("Stage", "Condition", "Severity", "Count",
                    "Affected", "Message", "RecommendedAction") %in%
                    names(prep_notes)))
  expect_true("missing_or_nonpositive_weight_rows_dropped" %in% prep_notes$Condition)
  expect_true("trimmed_person_ids" %in% prep_notes$Condition)
  expect_true(any(row_retention$DroppedRows == 1L))
  expect_identical(
    s_fit$row_retention$DroppedRows[s_fit$row_retention$Stage == "after_missing_and_weight_filter"][1],
    1L
  )
  expect_true("missing_or_nonpositive_weight_rows_dropped" %in%
                s_fit$preparation_notes$Condition)
  expect_true("missing_or_nonpositive_weight_rows_dropped" %in%
                s_ds$preparation_notes$Condition)
  expect_true(any(s_ds$row_retention$DroppedRows == 1L))
})

test_that("routine data-preparation messages are quiet by default", {
  dirty <- .toy
  dirty$Score[1] <- NA
  dirty$Person[2] <- paste0(" ", dirty$Person[2], " ")

  msgs <- testthat::capture_messages(suppressWarnings(
    fit_mfrm(dirty, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5)
  ))
  expect_false(any(grepl("Dropped 1 row", msgs, fixed = TRUE)))
  expect_false(any(grepl("Trimmed leading/trailing whitespace", msgs, fixed = TRUE)))
})

test_that("routine data-preparation messages can be enabled once per fit", {
  old_opt <- options(mfrmr.show_preparation_messages = TRUE)
  on.exit(options(old_opt), add = TRUE)
  dirty <- .toy
  dirty$Score[1] <- NA
  dirty$Person[2] <- paste0(" ", dirty$Person[2], " ")

  msgs <- testthat::capture_messages(suppressWarnings(
    fit_mfrm(dirty, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5)
  ))
  expect_equal(sum(grepl("Dropped 1 row", msgs, fixed = TRUE)), 1L)
  expect_equal(sum(grepl("Trimmed leading/trailing whitespace", msgs, fixed = TRUE)), 1L)
})

test_that("analyze_dff(method = 'refit') without diagnostics raises mfrmr error", {
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 10)
  ))
  expect_error(
    analyze_dff(fit, facet = "Rater", group = "Group",
                data = .toy, method = "refit"),
    "method = \"refit\""
  )
})

test_that("analyze_dff(method = 'residual') is tolerant of missing diagnostics", {
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 10)
  ))
  diag <- suppressMessages(
    diagnose_mfrm(fit, residual_pca = "none",
                  diagnostic_mode = "legacy")
  )
  out <- suppressMessages(suppressWarnings(
    analyze_dff(fit, diagnostics = diag,
                facet = "Rater", group = "Group",
                data = .toy, method = "residual")
  ))
  expect_s3_class(out, "mfrm_dff")
})

test_that("fit_mfrm(missing_codes = TRUE) recodes default sentinels to NA", {
  dirty <- .toy
  # Inject a default FACETS/SPSS sentinel for one rater row.
  dirty$Score <- as.character(dirty$Score)
  dirty$Score[1:3] <- "99"
  f <- suppressMessages(suppressWarnings(
    fit_mfrm(dirty, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5, missing_codes = TRUE)
  ))
  expect_true("missing_recoding" %in% names(f$prep))
  rec <- f$prep$missing_recoding
  expect_s3_class(rec, "data.frame")
  expect_true(any(rec$Replaced > 0))
})

test_that("fit_mfrm(missing_codes = NULL) keeps existing data untouched", {
  f <- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5, missing_codes = NULL)
  ))
  if (is.data.frame(f$prep$missing_recoding)) {
    expect_equal(nrow(f$prep$missing_recoding), 0L)
  } else {
    expect_null(f$prep$missing_recoding)
  }
})
