# Tests for the `missing_codes` integration (added in 0.1.6).
# Cross-refs:
#  * `fit_mfrm(..., missing_codes = ...)` — opt-in pre-processing
#  * `prepare_mfrm_data()` — stores `prep$missing_recoding`
#  * `build_mfrm_manifest()$missing_recoding` — reads prep audit
#  * `review_mfrm_anchors(..., missing_codes = ...)` — pass-through
#  * `describe_mfrm_data(..., missing_codes = ...)` — pass-through
#  * `recode_missing_codes()` — standalone helper (existing)

make_dirty_data <- function(seed = 42L, n_persons = 20L) {
  set.seed(seed)
  d <- expand.grid(
    Person = paste0("P", seq_len(n_persons)),
    Rater = c("R1", "R2"),
    Task = c("T1", "T2", "T3", "T4"),
    stringsAsFactors = FALSE
  )
  d$Score <- sample(c(0, 1, 2, 3, 99),
                    size = nrow(d),
                    replace = TRUE,
                    prob = c(0.20, 0.30, 0.25, 0.20, 0.05))
  d
}

test_that("missing_codes = TRUE drops FACETS sentinels from the fit", {
  d <- make_dirty_data(seed = 1L)
  fit_raw <- suppressMessages(suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task"), "Score",
             method = "JML", maxit = 20)
  ))
  fit_clean <- suppressMessages(suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task"), "Score",
             method = "JML", maxit = 20, missing_codes = TRUE)
  ))
  # 99 values should have become a top category in the raw fit but been
  # dropped when missing_codes is active.
  expect_gt(as.integer(fit_raw$summary$Categories),
            as.integer(fit_clean$summary$Categories))
  expect_gt(as.integer(fit_raw$summary$N),
            as.integer(fit_clean$summary$N))
})

test_that("prep$missing_recoding stores per-column replacement counts", {
  d <- make_dirty_data(seed = 2L)
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task"), "Score",
             method = "JML", maxit = 15, missing_codes = TRUE)
  ))
  audit <- fit$prep$missing_recoding
  expect_s3_class(audit, "data.frame")
  expect_true(all(c("Column", "Replaced") %in% names(audit)))
  expect_setequal(audit$Column, c("Person", "Rater", "Task", "Score"))
  expect_true(is.numeric(audit$Replaced) || is.integer(audit$Replaced))
  # 99 is a valid sentinel; expect at least one Score replacement.
  expect_gt(audit$Replaced[audit$Column == "Score"], 0L)
})

test_that("missing_codes accepts custom character vectors", {
  d <- data.frame(
    Person = rep(paste0("P", 1:10), each = 4),
    Rater  = rep(c("R1", "R2"), 20),
    Task   = rep(c("T1", "T2"), 20),
    Score  = rep(c(0, 1, 2, "MISSING"), 10),
    stringsAsFactors = FALSE
  )
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task"), "Score",
             method = "JML", maxit = 15,
             missing_codes = c("MISSING"))
  ))
  audit <- fit$prep$missing_recoding
  expect_equal(audit$Replaced[audit$Column == "Score"], 10L)
})

test_that("missing_codes = NULL is a strict no-op", {
  d <- make_dirty_data(seed = 3L)
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task"), "Score",
             method = "JML", maxit = 15)
  ))
  expect_null(fit$prep$missing_recoding)
})

test_that("build_mfrm_manifest exposes missing_recoding from the fit", {
  d <- make_dirty_data(seed = 4L)
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task"), "Score",
             method = "JML", maxit = 15, missing_codes = TRUE)
  ))
  m <- build_mfrm_manifest(fit)
  expect_true("missing_recoding" %in% names(m))
  expect_s3_class(m$missing_recoding, "data.frame")
  expect_true("Score" %in% m$missing_recoding$Column)
  expect_true(m$missing_recoding$Replaced[
    m$missing_recoding$Column == "Score"
  ] > 0L)
})

test_that("describe_mfrm_data respects missing_codes", {
  d <- make_dirty_data(seed = 5L)
  desc_raw <- suppressMessages(suppressWarnings(
    describe_mfrm_data(d, "Person", c("Rater", "Task"), "Score")
  ))
  desc_clean <- suppressMessages(suppressWarnings(
    describe_mfrm_data(d, "Person", c("Rater", "Task"), "Score",
                       missing_codes = TRUE)
  ))
  # Observation count must drop when 99 is coerced to NA.
  expect_gt(as.integer(desc_raw$overview$Observations),
            as.integer(desc_clean$overview$Observations))
  # Declared support should no longer include 99.
  raw_max <- as.integer(desc_raw$score_support$rating_max)
  clean_max <- as.integer(desc_clean$score_support$rating_max)
  expect_gt(raw_max, clean_max)
})

test_that("review_mfrm_anchors respects missing_codes", {
  d <- make_dirty_data(seed = 6L)
  aud_raw <- suppressMessages(suppressWarnings(
    review_mfrm_anchors(d, "Person", c("Rater", "Task"), "Score")
  ))
  aud_clean <- suppressMessages(suppressWarnings(
    review_mfrm_anchors(d, "Person", c("Rater", "Task"), "Score",
                       missing_codes = TRUE)
  ))
  # Facet summaries exist in both cases.
  expect_true(is.data.frame(aud_raw$facet_summary))
  expect_true(is.data.frame(aud_clean$facet_summary))
})
