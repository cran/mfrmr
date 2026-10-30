test_that("PAM entry names preserve partitions and saved-result methods", {
  skip_if_not_installed("cluster")
  x <- mfrm_features(data.frame(ID = letters[1:6],
    Experience = c(1, 2, 3, 12, 13, 14), Specialty = rep(c("A", "B"), each = 3)),
    "ID", c("Experience", "Specialty"))
  old <- mfrm_cluster(x, 2, c(Experience = 2, Specialty = 1))
  new <- mfrm_cluster_pam(x, 2, c(Experience = 2, Specialty = 1))
  expect_identical(new, old)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(old, path)
  expect_identical(summary(readRDS(path)), summary(new))
  expect_identical(plot_data(readRDS(path)), plot_data(new))
})

test_that("imputation-review names preserve selected events and reject ambiguous inputs", {
  d <- data.frame(Event = c("001", "1", "3", "4"), Person = c("A", "A", "B", "B"),
    Rater = c("X", "Y", "X", "Y"), Score = c(0, NA, 1, 2))
  a <- b <- d
  a$Score[2] <- 1; b$Score[2] <- 2
  provenance <- list(method = "supplied test completions")
  old <- mfrm_response_imputations(d, list(a, b), "Person", "Rater", "Score",
    "Event", "1", 0:2, imputation_model = provenance)
  new <- review_mfrm_imputations(d, list(a, b), "Person", "Rater", "Score",
    "Event", impute_ids = "1", categories = 0:2, imputation_model = provenance)
  expect_identical(new, old)
  expect_identical(new$events$Imputed, c(FALSE, TRUE, FALSE, FALSE))
  expect_error(review_mfrm_imputations(d, list(a, b), "Person", "Rater", "Score",
    "Event", impute_ids = TRUE, categories = 0:2, imputation_model = provenance), "event IDs")
  expect_error(review_mfrm_imputations(d, list(a, b), "Person", "Rater", "Score",
    "Event", impute_ids = "1", impute = "001", categories = 0:2,
    imputation_model = provenance), "unused argument")
  expect_error(fit_mfrm_imputed(new, category_policy = "collapse"), "supported shared fitting")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(old, path)
  expect_identical(summary(readRDS(path)), summary(new))
})

test_that("explicit category policies agree across review and fitting without changing defaults", {
  d <- load_mfrmr_data("example_core")
  args <- list(data = d, person = "Person", facets = c("Rater", "Criterion"),
    score = "Score", rating_min = 1, rating_max = 4)
  gapped <- args
  gapped$data$Score[gapped$data$Score == 2] <- 3
  collapsed_review <- function(fun, arguments) {
    expect_warning(result <- do.call(fun, arguments), "recoded internally")
    result
  }
  for (fun in list(describe_mfrm_data, review_mfrm_anchors)) {
    for (preserve in c(FALSE, TRUE)) {
      run <- if (preserve) do.call else collapsed_review
      expect_identical(run(fun, c(gapped, list(keep_original = preserve))),
        run(fun, c(gapped, list(category_policy = if (preserve) "preserve" else "collapse"))))
    }
    expect_identical(collapsed_review(fun, gapped),
      collapsed_review(fun, c(gapped, list(category_policy = "collapse"))))
  }
  for (fun in list(fit_mfrm, describe_mfrm_data, review_mfrm_anchors)) {
    expect_error(do.call(fun, c(args, list(keep_original = FALSE, category_policy = "preserve"))),
      "conflict")
    expect_error(do.call(fun, c(args, list(keep_original = TRUE, category_policy = "collapse"))),
      "conflict")
    for (invalid in list(NA_character_, TRUE, character(), c("collapse", "preserve"))) {
      expect_error(do.call(fun, c(args, list(category_policy = invalid))), "category_policy")
    }
  }
  old <- do.call(fit_mfrm, c(args, list(keep_original = TRUE, method = "JML", maxit = 400)))
  new <- do.call(fit_mfrm, c(args, list(category_policy = "preserve", method = "JML", maxit = 400)))
  expect_equal(new$opt$par, old$opt$par, tolerance = 0)
  expect_identical(new$prep$score_map, old$prep$score_map)
  expect_identical(new$config$replay_inputs, old$config$replay_inputs)
  expect_true(new$config$replay_inputs$keep_original)
  settings <- summary(new)$settings_overview
  expect_identical(settings$CategoryPolicy, "preserve")
  expect_false(settings$ScoreRecoded)
  expect_identical(build_summary_table_bundle(summary(new))$tables$settings_overview$CategoryPolicy,
    "preserve")
  expect_error(do.call(fit_mfrm, c(gapped, list(category_policy = "preserve"))),
    class = "mfrmr_category_readiness_error")
  # A matched old/new choice is valid and retains the old data-review result.
  expect_identical(do.call(describe_mfrm_data, c(gapped, list(keep_original = TRUE))),
    do.call(describe_mfrm_data, c(gapped, list(keep_original = TRUE, category_policy = "preserve"))))
})

test_that("category summaries distinguish chosen policy, recoding and missing legacy records", {
  d <- load_mfrmr_data("example_core")
  review <- function(data, ...) describe_mfrm_data(data, "Person",
    c("Rater", "Criterion"), "Score", rating_min = 1, rating_max = 4, ...)
  contiguous <- review(d)
  expect_identical(contiguous$overview$CategoryPolicy, "collapse")
  expect_false(contiguous$overview$ScoreRecoded)
  d$Score[d$Score == 2] <- 3
  expect_warning(collapsed <- review(d), "recoded internally")
  preserved <- review(d, category_policy = "preserve")
  expect_identical(collapsed$overview$CategoryPolicy, "collapse")
  expect_true(collapsed$overview$ScoreRecoded)
  expect_identical(preserved$overview$CategoryPolicy, "preserve")
  expect_false(preserved$overview$ScoreRecoded)
  expect_true(build_summary_table_bundle(summary(collapsed))$tables$overview$ScoreRecoded)

  # Identical category labels cannot identify the policy in an older result.
  legacy_review <- contiguous
  legacy_review$overview[c("CategoryPolicy", "ScoreRecoded")] <- NULL
  legacy_review$score_support$keep_original <- NULL
  expect_identical(summary(legacy_review)$overview$CategoryPolicy, "not_recorded")
  expect_false(summary(legacy_review)$overview$ScoreRecoded)
  legacy_review$score_support$score_map <- NULL
  expect_true(is.na(summary(legacy_review)$overview$ScoreRecoded))

  fit <- make_toy_fit()
  expect_identical(summary(fit)$settings_overview$CategoryPolicy, "collapse")
  fit$prep$keep_original <- fit$config$keep_original <- NULL
  expect_identical(summary(fit)$settings_overview$CategoryPolicy, "collapse")
  fit$config$replay_inputs$keep_original <- NULL
  expect_identical(summary(fit)$settings_overview$CategoryPolicy, "not_recorded")
  expect_false(summary(fit)$settings_overview$ScoreRecoded)
  fit$prep$score_map <- fit$config$score_map <- NULL
  expect_true(is.na(summary(fit)$settings_overview$ScoreRecoded))
})

test_that("the migration table exposes canonical functions and category meaning", {
  aliases <- compatibility_alias_table()
  expect_identical(aliases$PreferredName[match(c("mfrm_cluster", "mfrm_response_imputations"), aliases$Alias)],
    c("mfrm_cluster_pam", "review_mfrm_imputations"))
  expect_identical(aliases$PreferredName[aliases$Alias == "keep_original"], "category_policy")
})
