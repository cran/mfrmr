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

make_population_model_prediction_fixture <- function() {
  mfrmr:::with_preserved_rng_seed(20260403, {
    persons <- paste0("P", sprintf("%02d", 1:60))
    items <- paste0("I", 1:6)
    x <- stats::rnorm(length(persons))
    theta <- 0.25 + 0.9 * x + stats::rnorm(length(persons), sd = 0.6)
    item_beta <- seq(-1.0, 1.0, length.out = length(items))

    dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
    eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
    dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))

    person_tbl <- data.frame(
      Person = persons,
      X = x,
      stringsAsFactors = FALSE
    )

    fit <- suppressWarnings(
      fit_mfrm(
        dat,
        "Person", "Item", "Score",
        method = "MML",
        model = "RSM",
        population_formula = ~ X,
        person_data = person_tbl,
        quad_points = 7,
        maxit = 80
      )
    )

    new_units <- data.frame(
      Person = c("NEW_LOW", "NEW_LOW", "NEW_HIGH", "NEW_HIGH"),
      Item = c(items[1], items[2], items[1], items[2]),
      Score = c(1, 0, 1, 0),
      stringsAsFactors = FALSE
    )
    new_person_data <- data.frame(
      Person = c("NEW_LOW", "NEW_HIGH"),
      X = c(-1.5, 1.5),
      stringsAsFactors = FALSE
    )

    list(
      fit = fit,
      new_units = new_units,
      person_data = new_person_data
    )
  })
}

make_mean_only_population_model_prediction_fixture <- function() {
  mfrmr:::with_preserved_rng_seed(20260410, {
    persons <- paste0("P", sprintf("%02d", 1:60))
    items <- paste0("I", 1:6)
    theta <- 0.35 + stats::rnorm(length(persons), sd = 0.7)
    item_beta <- seq(-1.0, 1.0, length.out = length(items))

    dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
    eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
    dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))

    fit <- suppressWarnings(
      fit_mfrm(
        dat,
        "Person", "Item", "Score",
        method = "MML",
        model = "RSM",
        population_formula = ~ 1,
        person_data = data.frame(Person = persons, stringsAsFactors = FALSE),
        quad_points = 7,
        maxit = 80
      )
    )

    new_units <- data.frame(
      Person = c("NEW01", "NEW01", "NEW02", "NEW02"),
      Item = c(items[1], items[2], items[1], items[2]),
      Score = c(1, 0, 1, 0),
      stringsAsFactors = FALSE
    )

    list(
      fit = fit,
      new_units = new_units
    )
  })
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
  printed <- capture.output(print(summary(pred)))
  expect_true(any(grepl("mfrmr Unit Prediction Summary", printed, fixed = TRUE)))
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

test_that("predict_mfrm_units requires person_data when latent-regression scoring uses covariates", {
  fixture <- make_population_model_prediction_fixture()

  expect_error(
    predict_mfrm_units(fixture$fit, fixture$new_units),
    "`person_data` must be supplied when scoring a latent-regression fit with background covariates.",
    fixed = TRUE
  )
})

test_that("predict_mfrm_units scores latent-regression fits under the fitted population model", {
  fixture <- make_population_model_prediction_fixture()
  swapped_person_data <- fixture$person_data
  swapped_person_data$X <- rev(swapped_person_data$X)
  pred <- predict_mfrm_units(
    fixture$fit,
    fixture$new_units,
    person_data = fixture$person_data,
    n_draws = 3,
    seed = 17
  )
  pred_swapped <- predict_mfrm_units(
    fixture$fit,
    fixture$new_units,
    person_data = swapped_person_data,
    n_draws = 0
  )
  pv <- sample_mfrm_plausible_values(
    fixture$fit,
    fixture$new_units,
    person_data = fixture$person_data,
    n_draws = 3,
    seed = 18
  )

  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_identical(pred$settings$posterior_basis, "population_model")
  expect_true(is.data.frame(pred$population_audit))
  expect_equal(pred$population_audit$RetainedPersons, 2)
  expect_equal(pred$population_audit$OmittedPersons, 0)
  expect_true(is.data.frame(pred$person_data))
  expect_equal(sort(unique(pred$person_data$Person)), c("NEW_HIGH", "NEW_LOW"))
  expect_identical(as.character(pred$estimates$Person), c("NEW_LOW", "NEW_HIGH"))
  est <- pred$estimates$Estimate[match(c("NEW_LOW", "NEW_HIGH"), pred$estimates$Person)]
  est_swapped <- pred_swapped$estimates$Estimate[match(c("NEW_LOW", "NEW_HIGH"), pred_swapped$estimates$Person)]
  expect_true(is.finite(sum(abs(est - est_swapped))))
  expect_gt(est[2], est[1])
  expect_gt(est_swapped[1], est[1])
  expect_lt(est_swapped[2], est[2])
  expect_true(any(grepl("conditional normal population model", pred$notes, fixed = TRUE)))

  expect_s3_class(pv, "mfrm_plausible_values")
  expect_identical(pv$settings$posterior_basis, "population_model")
  expect_true(is.data.frame(pv$population_audit))
  expect_true(is.data.frame(pv$person_data))
  expect_equal(nrow(pv$values), 6)
})

test_that("latent-regression prediction reuses categorical model-matrix levels", {
  mfrmr:::with_preserved_rng_seed(20260411, {
    persons <- paste0("P", sprintf("%02d", 1:50))
    items <- paste0("I", 1:4)
    group <- rep(c("low", "high"), length.out = length(persons))
    theta <- ifelse(group == "high", 0.8, -0.4) + stats::rnorm(length(persons), sd = 0.6)
    item_beta <- seq(-0.8, 0.8, length.out = length(items))

    dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
    eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
    dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))
    person_tbl <- data.frame(
      Person = persons,
      Group = group,
      stringsAsFactors = FALSE
    )

    fit <- suppressWarnings(
      fit_mfrm(
        dat,
        "Person", "Item", "Score",
        method = "MML",
        model = "RSM",
        population_formula = ~ Group,
        person_data = person_tbl,
        quad_points = 5,
        maxit = 60
      )
    )

    new_units <- data.frame(
      Person = c("NEW_HIGH", "NEW_HIGH"),
      Item = items[1:2],
      Score = c(1, 0),
      stringsAsFactors = FALSE
    )
    single_level_person_data <- data.frame(
      Person = "NEW_HIGH",
      Group = "high",
      stringsAsFactors = FALSE
    )

    pred <- predict_mfrm_units(
      fit,
      new_units,
      person_data = single_level_person_data,
      n_draws = 0
    )

    expect_s3_class(pred, "mfrm_unit_prediction")
    expect_identical(pred$settings$posterior_basis, "population_model")
    expect_equal(pred$population_audit$RetainedPersons, 1)
    expect_equal(pred$population_audit$OmittedPersons, 0)
    expect_equal(pred$estimates$Person, "NEW_HIGH")
    expect_true("Group" %in% names(fit$population$xlevels))
    expect_true(any(startsWith(fit$population$design_columns, "Group")))

    bad_person_data <- single_level_person_data
    bad_person_data$Group <- "unseen"
    expect_error(
      predict_mfrm_units(
        fit,
        new_units,
        person_data = bad_person_data,
        n_draws = 0
      ),
      "Could not build the latent-regression model frame",
      fixed = TRUE
    )
  })
})

test_that("latent-regression prediction can omit scored persons with incomplete background data", {
  fixture <- make_population_model_prediction_fixture()
  person_data <- fixture$person_data
  person_data$X[1] <- NA_real_

  pred <- predict_mfrm_units(
    fixture$fit,
    fixture$new_units,
    person_data = person_data,
    population_policy = "omit",
    n_draws = 0
  )
  pv <- sample_mfrm_plausible_values(
    fixture$fit,
    fixture$new_units,
    person_data = person_data,
    population_policy = "omit",
    n_draws = 3,
    seed = 19
  )

  expect_equal(pred$population_audit$RetainedPersons, 1)
  expect_equal(pred$population_audit$OmittedPersons, 1)
  expect_equal(pred$population_audit$OmittedRows, 2)
  expect_equal(pred$estimates$Person, "NEW_HIGH")
  expect_equal(sort(unique(pred$input_data$Person)), "NEW_HIGH")
  expect_equal(pv$population_audit$RetainedPersons, 1)
  expect_equal(pv$population_audit$OmittedPersons, 1)
  expect_equal(pv$population_audit$OmittedRows, 2)
  expect_equal(pv$estimates$Person, "NEW_HIGH")
  expect_equal(sort(unique(pv$values$Person)), "NEW_HIGH")
  expect_equal(nrow(pv$values), 3)
})

test_that("latent-regression prediction errors on incomplete covariates under the default policy", {
  fixture <- make_population_model_prediction_fixture()
  person_data <- fixture$person_data
  person_data$X[1] <- NA_real_

  expect_error(
    predict_mfrm_units(
      fixture$fit,
      fixture$new_units,
      person_data = person_data
    ),
    "contains missing covariate values",
    fixed = TRUE
  )
})

test_that("latent-regression scoring changes scored contrasts relative to legacy MML", {
  fixture <- make_population_model_prediction_fixture()

  legacy_fit <- suppressWarnings(
    fit_mfrm(
      fixture$fit$prep$data,
      "Person", "Item", "Score",
      method = "MML",
      model = "RSM",
      quad_points = 7,
      maxit = 80
    )
  )

  pred_population <- predict_mfrm_units(
    fixture$fit,
    fixture$new_units,
    person_data = fixture$person_data,
    n_draws = 0
  )
  pred_legacy <- predict_mfrm_units(
    legacy_fit,
    fixture$new_units,
    n_draws = 0
  )

  est_population <- pred_population$estimates$Estimate[
    match(c("NEW_LOW", "NEW_HIGH"), pred_population$estimates$Person)
  ]
  est_legacy <- pred_legacy$estimates$Estimate[
    match(c("NEW_LOW", "NEW_HIGH"), pred_legacy$estimates$Person)
  ]

  expect_identical(pred_population$settings$posterior_basis, "population_model")
  expect_identical(pred_legacy$settings$posterior_basis, "legacy_mml")
  expect_lt(abs(diff(est_legacy)), 0.01)
  expect_gt(diff(est_population), 0.50)
  expect_gt(diff(est_population), diff(est_legacy) + 0.50)
})

test_that("predict_mfrm_units supports intercept-only latent-regression scoring without person_data", {
  fixture <- make_mean_only_population_model_prediction_fixture()

  pred <- predict_mfrm_units(
    fixture$fit,
    fixture$new_units,
    n_draws = 2,
    seed = 27
  )
  pv <- sample_mfrm_plausible_values(
    fixture$fit,
    fixture$new_units,
    n_draws = 2,
    seed = 28
  )

  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_identical(pred$settings$posterior_basis, "population_model")
  expect_equal(pred$population_audit$RetainedPersons, 2)
  expect_equal(pred$population_audit$OmittedPersons, 0)
  expect_true(is.data.frame(pred$person_data))
  expect_identical(names(pred$person_data), "Person")
  expect_equal(sort(unique(pred$person_data$Person)), c("NEW01", "NEW02"))
  expect_true(any(grepl("intercept-only population model", pred$notes, fixed = TRUE)))

  expect_s3_class(pv, "mfrm_plausible_values")
  expect_identical(pv$settings$posterior_basis, "population_model")
  expect_true(is.data.frame(pv$person_data))
  expect_identical(names(pv$person_data), "Person")
  expect_equal(nrow(pv$values), 4)
})

test_that("predict_mfrm_units supports JML PCM calibrations with custom facet names", {
  toy <- load_mfrmr_data("example_core")
  names(toy)[names(toy) == "Rater"] <- "Judge"
  names(toy)[names(toy) == "Criterion"] <- "Task"
  keep_people <- unique(toy$Person)[1:14]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit_jml <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Judge", "Task"), "Score",
      method = "JML",
      model = "PCM",
      step_facet = "Task",
      maxit = 15
    )
  )

  new_units <- data.frame(
    Person = c("NEW01", "NEW01"),
    Judge = unique(toy$Judge)[1],
    Task = unique(toy$Task)[1:2],
    Score = c(2, 3)
  )

  pred <- predict_mfrm_units(fit_jml, new_units, n_draws = 2, seed = 11)
  pv <- sample_mfrm_plausible_values(fit_jml, new_units, n_draws = 2, seed = 12)

  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_equal(pred$settings$method, "JML")
  expect_true(any(grepl("standard normal reference prior", pred$notes, fixed = TRUE)))
  expect_true(all(c("Person", "Judge", "Task", "Score", "Weight") %in% names(pred$input_data)))
  expect_s3_class(pv, "mfrm_plausible_values")
  expect_equal(pv$settings$method, "JML")
  expect_true(any(grepl("standard normal reference prior", pv$notes, fixed = TRUE)))
  expect_equal(nrow(pv$values), 2)
})

test_that("predict_mfrm_units supports bounded GPCM fixed-calibration scoring", {
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:16]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]
  fit_gpcm <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Rater", "Criterion"), "Score",
      method = "MML",
      model = "GPCM",
      step_facet = "Criterion",
      quad_points = 5,
      maxit = 20
    )
  )

  new_units <- data.frame(
    Person = c("NEW01", "NEW01", "NEW02", "NEW02"),
    Rater = c(unique(toy$Rater)[1], unique(toy$Rater)[2], unique(toy$Rater)[1], unique(toy$Rater)[2]),
    Criterion = c(unique(toy$Criterion)[1], unique(toy$Criterion)[2], unique(toy$Criterion)[1], unique(toy$Criterion)[2]),
    Score = c(2, 3, 2, 4)
  )

  pred <- predict_mfrm_units(fit_gpcm, new_units, n_draws = 2, seed = 31)
  pv <- sample_mfrm_plausible_values(fit_gpcm, new_units, n_draws = 2, seed = 32)

  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_equal(pred$settings$method, "MML")
  expect_identical(pred$settings$posterior_basis, "legacy_mml")
  expect_equal(sort(unique(pred$estimates$Person)), c("NEW01", "NEW02"))
  expect_true(all(is.finite(pred$estimates$Estimate)))
  expect_true(any(grepl("fixed fitted MML calibration", pred$notes, fixed = TRUE)))
  expect_s3_class(pv, "mfrm_plausible_values")
  expect_equal(nrow(pv$values), 4)
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
  printed <- capture.output(print(summary(pv)))
  expect_true(any(grepl("mfrmr Plausible Values Summary", printed, fixed = TRUE)))
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
