make_prediction_fixture <- function() {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Rater", "Criterion"), "Score",
      method = "MML",
      quad_points = 5,
      maxit = 15
    )
  )
  raters <- unique(toy$Rater)[1:2]
  criteria <- unique(toy$Criterion)[1:2]
  new_units <- data.frame(
    Person = c("NEW01", "NEW01", "NEW02", "NEW02"),
    Rater = c(raters[1], raters[2], raters[1], raters[2]),
    Criterion = c(criteria[1], criteria[2], criteria[1], criteria[2]),
    Score = c(2, 3, 2, 4)
  )
  list(fit = fit, new_units = new_units)
}

test_that("predict_mfrm_units returns posterior summaries and optional draws", {
  fixture <- make_prediction_fixture()
  pred <- predict_mfrm_units(
    fixture$fit,
    fixture$new_units,
    interval_level = 0.8,
    n_draws = 3,
    seed = 42
  )

  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_true(is.data.frame(pred$estimates))
  expect_true(all(c("Person", "Estimate", "SD", "Lower", "Upper",
                    "Observations", "WeightedN") %in% names(pred$estimates)))
  expect_equal(sort(unique(pred$estimates$Person)), c("NEW01", "NEW02"))
  expect_true(all(pred$estimates$Lower <= pred$estimates$Estimate))
  expect_true(all(pred$estimates$Estimate <= pred$estimates$Upper))
  expect_true(is.data.frame(pred$draws))
  expect_true(is.data.frame(pred$input_data))
  expect_true(all(c("Person", "Rater", "Criterion", "Score", "Weight") %in% names(pred$input_data)))
  expect_equal(nrow(pred$draws), 6)

  s <- summary(pred)
  expect_s3_class(s, "summary.mfrm_unit_prediction")
  expect_true(is.data.frame(s$estimates))
  expect_true(is.data.frame(s$audit))
})

test_that("predict_mfrm_units supports explicit column remapping and weights", {
  fixture <- make_prediction_fixture()
  remapped <- data.frame(
    Candidate = fixture$new_units$Person,
    Judge = fixture$new_units$Rater,
    Dimension = fixture$new_units$Criterion,
    Rating = fixture$new_units$Score,
    Wt = c(1, 2, 1, 2)
  )

  pred <- predict_mfrm_units(
    fixture$fit,
    remapped,
    person = "Candidate",
    facets = c(Rater = "Judge", Criterion = "Dimension"),
    score = "Rating",
    weight = "Wt",
    n_draws = 2,
    seed = 7
  )

  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_true(all(c("Person", "Rater", "Criterion", "Score", "Weight") %in% names(pred$input_data)))
  expect_equal(pred$input_data$Weight, remapped$Wt)
  expect_equal(pred$settings$source_columns$person, "Candidate")
  expect_equal(unname(pred$settings$source_columns$facets), c("Judge", "Dimension"))
  expect_equal(pred$settings$source_columns$score, "Rating")
  expect_equal(pred$settings$source_columns$weight, "Wt")
})

test_that("prediction draws are reproducible for a fixed seed", {
  fixture <- make_prediction_fixture()

  pred_a <- predict_mfrm_units(fixture$fit, fixture$new_units, n_draws = 3, seed = 91)
  pred_b <- predict_mfrm_units(fixture$fit, fixture$new_units, n_draws = 3, seed = 91)
  pv_a <- sample_mfrm_plausible_values(fixture$fit, fixture$new_units, n_draws = 3, seed = 92)
  pv_b <- sample_mfrm_plausible_values(fixture$fit, fixture$new_units, n_draws = 3, seed = 92)

  expect_identical(pred_a$draws, pred_b$draws)
  expect_identical(pv_a$values, pv_b$values)
})

test_that("predict_mfrm_units rejects unseen facet levels", {
  fixture <- make_prediction_fixture()
  bad_new <- fixture$new_units
  bad_new$Rater[1] <- "UNKNOWN_RATER"

  expect_error(
    predict_mfrm_units(fixture$fit, bad_new),
    "unseen levels for facet `Rater`",
    fixed = TRUE
  )
})

test_that("predict_mfrm_units requires MML calibration", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:14]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit_jml <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Rater", "Criterion"), "Score",
      method = "JML",
      maxit = 15
    )
  )

  expect_error(
    predict_mfrm_units(
      fit_jml,
      data.frame(
        Person = c("NEW01", "NEW01"),
        Rater = unique(toy$Rater)[1],
        Criterion = unique(toy$Criterion)[1:2],
        Score = c(2, 3)
      )
    ),
    "supports only fits estimated with method = 'MML'",
    fixed = TRUE
  )
})

test_that("predict_mfrm_units respects stored score mapping from compressed fits", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:16]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  score_map <- c(`1` = 1L, `2` = 3L, `3` = 5L, `4` = 7L)
  toy$Score <- unname(score_map[as.character(toy$Score)])

  fit <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Rater", "Criterion"), "Score",
      method = "MML",
      keep_original = FALSE,
      quad_points = 5,
      maxit = 15
    )
  )

  new_units <- data.frame(
    Person = c("NEW01", "NEW01"),
    Rater = unique(toy$Rater)[1],
    Criterion = unique(toy$Criterion)[1:2],
    Score = c(1, 7)
  )

  pred <- predict_mfrm_units(fit, new_units)
  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_equal(nrow(pred$estimates), 1)
  expect_true(isTRUE(all.equal(sort(fit$prep$score_map$OriginalScore), c(1, 3, 5, 7))))
})

test_that("sample_mfrm_plausible_values exposes fixed-calibration posterior draws", {
  fixture <- make_prediction_fixture()

  pv <- sample_mfrm_plausible_values(
    fixture$fit,
    fixture$new_units,
    n_draws = 4,
    seed = 99
  )

  expect_s3_class(pv, "mfrm_plausible_values")
  expect_true(is.data.frame(pv$values))
  expect_true(is.data.frame(pv$estimates))
  expect_true(is.data.frame(pv$input_data))
  expect_equal(sort(unique(pv$values$Person)), c("NEW01", "NEW02"))
  expect_equal(nrow(pv$values), 8)

  s <- summary(pv)
  expect_s3_class(s, "summary.mfrm_plausible_values")
  expect_true(is.data.frame(s$draw_summary))
  expect_true(all(c("Person", "Draws", "MeanValue", "SDValue",
                    "LowerValue", "UpperValue") %in% names(s$draw_summary)))
  expect_true(all(s$draw_summary$Draws == 4))
})

test_that("sample_mfrm_plausible_values requires positive draw count", {
  fixture <- make_prediction_fixture()

  expect_error(
    sample_mfrm_plausible_values(
      fixture$fit,
      fixture$new_units,
      n_draws = 0
    ),
    "`n_draws` must be a positive integer.",
    fixed = TRUE
  )
})

test_that("prediction preprocessing warns and audits dropped invalid rows", {
  fixture <- make_prediction_fixture()
  new_units <- data.frame(
    Candidate = c("NEW01", "NEW01", "NEW02", "NEW03"),
    Judge = c(fixture$new_units$Rater[1], fixture$new_units$Rater[2], fixture$new_units$Rater[1], fixture$new_units$Rater[2]),
    Dimension = c(fixture$new_units$Criterion[1], fixture$new_units$Criterion[2], fixture$new_units$Criterion[1], fixture$new_units$Criterion[2]),
    Rating = c(2, 3, NA, 4),
    Wt = c(1, 1, 1, 0)
  )

  expect_warning(
    pred <- predict_mfrm_units(
      fixture$fit,
      new_units,
      person = "Candidate",
      facets = c(Rater = "Judge", Criterion = "Dimension"),
      score = "Rating",
      weight = "Wt"
    ),
    "Dropped 2 row\\(s\\) from `new_data` before posterior scoring"
  )

  expect_equal(pred$audit$InputRows, 4)
  expect_equal(pred$audit$KeptRows, 2)
  expect_equal(pred$audit$DroppedRows, 2)
  expect_equal(pred$audit$DroppedMissing, 1)
  expect_equal(pred$audit$DroppedNonpositiveWeight, 1)
})

test_that("prediction integer validation does not leak coercion warnings", {
  fixture <- make_prediction_fixture()
  pred <- predict_mfrm_units(fixture$fit, fixture$new_units, n_draws = 0)

  expect_no_warning(
    expect_error(
      predict_mfrm_units(fixture$fit, fixture$new_units, n_draws = "foo"),
      "`n_draws` must be a non-negative integer.",
      fixed = TRUE
    )
  )

  expect_no_warning(
    expect_error(
      sample_mfrm_plausible_values(fixture$fit, fixture$new_units, n_draws = "foo"),
      "`n_draws` must be a positive integer.",
      fixed = TRUE
    )
  )

  expect_no_warning(
    expect_error(
      summary(pred, digits = "foo"),
      "`digits` must be a non-negative integer.",
      fixed = TRUE
    )
  )
})
