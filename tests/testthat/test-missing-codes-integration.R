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
  expect_true(all(c("Column", "Replaced", "Scope") %in% names(audit)))
  expect_setequal(audit$Column, c("Person", "Rater", "Task", "Score"))
  expect_true(is.numeric(audit$Replaced) || is.integer(audit$Replaced))
  # 99 is a valid sentinel; expect at least one Score replacement.
  expect_gt(audit$Replaced[audit$Column == "Score"], 0L)
  expect_identical(
    audit$Scope[audit$Column == "Score"],
    "default_score_only"
  )
})

test_that("default missing codes preserve legitimate person and facet IDs", {
  d <- make_dirty_data(seed = 22L)
  d$Person[d$Person == "P1"] <- "N"
  d$Rater[d$Rater == "R1"] <- "N"

  prep <- suppressWarnings(mfrmr:::prepare_mfrm_data(
    d,
    person_col = "Person",
    facet_cols = c("Rater", "Task"),
    score_col = "Score",
    missing_codes = TRUE
  ))

  expect_true("N" %in% prep$data$Person)
  expect_true("N" %in% prep$data$Rater)
  expect_equal(
    prep$missing_recoding$Replaced[
      prep$missing_recoding$Column %in% c("Person", "Rater", "Task")
    ],
    rep(0L, 3L)
  )
  expect_true(all(
    prep$missing_recoding$Scope[
      prep$missing_recoding$Column %in% c("Person", "Rater", "Task")
    ] == "identifier_preserved_by_default"
  ))
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

test_that("unhandled FACETS sentinels are warned once and remain auditable", {
  d <- make_dirty_data(seed = 31L)
  warnings <- character(0)
  fit <- withCallingHandlers(
    suppressMessages(
      fit_mfrm(
        d, "Person", c("Rater", "Task"), "Score",
        method = "JML", maxit = 15, missing_codes = NULL
      )
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  score_warnings <- warnings[grepl("non-consecutive", warnings, fixed = TRUE)]
  expect_length(score_warnings, 1L)
  expect_match(score_warnings, "commonly used as missing-value codes", fixed = TRUE)

  notes <- as.data.frame(fit$prep$preparation_notes, stringsAsFactors = FALSE)
  expect_true(any(notes$Condition == "score_categories_recoded"))
  expect_match(
    notes$Message[notes$Condition == "score_categories_recoded"][1],
    "missing-value codes",
    fixed = TRUE
  )
  expect_true(any(fit$prep$score_map$OriginalScore == 99L))
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
  expect_s3_class(desc_clean$missing_recoding, "data.frame")
  expect_gt(
    desc_clean$missing_recoding$Replaced[
      desc_clean$missing_recoding$Column == "Score"
    ],
    0L
  )
})

test_that("data description does not silently call an item facet a rater", {
  d <- expand.grid(
    Person = paste0("P", seq_len(12)),
    Item = paste0("I", seq_len(4)),
    stringsAsFactors = FALSE
  )
  d$Score <- rep(c(0L, 1L), length.out = nrow(d))

  auto <- suppressWarnings(
    describe_mfrm_data(d, "Person", "Item", "Score")
  )
  expect_false(auto$agreement$settings$included)
  expect_identical(auto$agreement$settings$status, "not_computed")
  expect_equal(nrow(auto$agreement$summary), 0L)
  expect_identical(
    summary(auto)$reporting_map$CoveredHere[
      summary(auto)$reporting_map$Area == "Rater-facet agreement"
    ],
    "no"
  )

  explicit <- suppressWarnings(
    describe_mfrm_data(
      d, "Person", "Item", "Score",
      rater_facet = "Item"
    )
  )
  expect_true(explicit$agreement$settings$included)
  expect_identical(explicit$agreement$settings$status, "computed")
  expect_identical(explicit$agreement$settings$rater_facet, "Item")
  explicit_summary <- summary(explicit)
  expect_identical(explicit_summary$agreement_settings$rater_facet, "Item")
  expect_match(
    paste(capture.output(print(explicit_summary)), collapse = "\n"),
    "Observed agreement by Item",
    fixed = TRUE
  )
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
