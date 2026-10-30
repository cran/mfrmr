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
  expect_true(is.data.frame(s$row_review))
  printed <- capture.output(print(summary(pred)))
  expect_true(any(grepl("mfrmr Unit Prediction Summary", printed, fixed = TRUE)))
  alternate <- predict_mfrm_units(fixture$fit, fixture$new_units,
    scoring_prior = list(mean = 0.25, sd = 1.2))
  expect_true(all(alternate$estimates$RetainedPriorMean == 0))
  expect_true(all(alternate$estimates$RetainedPriorSD == 1))
  expect_true(all(alternate$estimates$ScoreIntegrationReady))
  expect_s3_class(summary(alternate), "summary.mfrm_unit_prediction")
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

  expect_error(predict_mfrm_units(fixture$fit, fixture$new_units,
    scoring_prior = list(mean = 0, sd = 1), readiness_policy = "review"),
    "cannot override a population model with background covariates", fixed = TRUE)

  expect_error(
    predict_mfrm_units(fixture$fit, fixture$new_units, readiness_policy = "review"),
    "`person_data` must be supplied when scoring a latent-regression fit with background covariates.",
    fixed = TRUE
  )
})

test_that("predict_mfrm_units scores latent-regression fits under the fitted population model", {
  fixture <- make_population_model_prediction_fixture()
  swapped_person_data <- fixture$person_data
  swapped_person_data$X <- rev(swapped_person_data$X)
  pred <- predict_mfrm_units(
    fixture$fit, readiness_policy = "review",
    fixture$new_units,
    person_data = fixture$person_data,
    n_draws = 3,
    seed = 17
  )
  pred_swapped <- predict_mfrm_units(
    fixture$fit, readiness_policy = "review",
    fixture$new_units,
    person_data = swapped_person_data,
    n_draws = 0
  )
  reviewed <- predict_mfrm_units(
    fixture$fit, fixture$new_units, readiness_policy = "review",
    person_data = fixture$person_data, n_draws = 3, seed = 17,
    adaptive_quad_points = c(15L, 31L)
  )
  expect_identical(reviewed$estimates, pred$estimates)
  expect_identical(reviewed$draws, pred$draws)
  expect_true(all(reviewed$quadrature_review$Status == "computed"))
  expect_equal(reviewed$quadrature_review$FixedEAP[
                 match(pred$estimates$Person, reviewed$quadrature_review$Person)],
               unname(pred$estimates$Estimate), tolerance = 1e-12)
  expect_identical(reviewed$settings$source_scoring_ready, FALSE)
  pv <- sample_mfrm_plausible_values(
    fixture$fit, readiness_policy = "review",
    fixture$new_units,
    person_data = fixture$person_data,
    n_draws = 3,
    seed = 18
  )

  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_identical(pred$settings$posterior_basis, "population_model")
  expected_prior <- as.numeric(cbind(1, fixture$person_data$X) %*% fixture$fit$population$coefficients)
  expect_equal(pred$estimates$PriorMean[match(fixture$person_data$Person, pred$estimates$Person)], expected_prior)
  expect_equal(pred$estimates$PriorSD, rep(sqrt(fixture$fit$population$sigma2), nrow(pred$estimates)))
  expect_true(is.data.frame(pred$population_review))
  expect_equal(pred$population_review$RetainedPersons, 2)
  expect_equal(pred$population_review$OmittedPersons, 0)
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
  expect_true(is.data.frame(pv$population_review))
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
      fit, readiness_policy = "review",
      new_units,
      person_data = single_level_person_data,
      n_draws = 0
    )

    expect_s3_class(pred, "mfrm_unit_prediction")
    expect_identical(pred$settings$posterior_basis, "population_model")
    expect_equal(pred$population_review$RetainedPersons, 1)
    expect_equal(pred$population_review$OmittedPersons, 0)
    expect_equal(pred$estimates$Person, "NEW_HIGH")
    expect_true("Group" %in% names(fit$population$xlevels))
    expect_true(any(startsWith(fit$population$design_columns, "Group")))

    bad_person_data <- single_level_person_data
    bad_person_data$Group <- "unseen"
    expect_error(
      predict_mfrm_units(
        fit, readiness_policy = "review",
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
    fixture$fit, readiness_policy = "review",
    fixture$new_units,
    person_data = person_data,
    population_policy = "omit",
    n_draws = 0
  )
  pv <- sample_mfrm_plausible_values(
    fixture$fit, readiness_policy = "review",
    fixture$new_units,
    person_data = person_data,
    population_policy = "omit",
    n_draws = 3,
    seed = 19
  )

  expect_equal(pred$population_review$RetainedPersons, 1)
  expect_equal(pred$population_review$OmittedPersons, 1)
  expect_equal(pred$population_review$OmittedRows, 2)
  expect_equal(pred$estimates$Person, "NEW_HIGH")
  expect_equal(sort(unique(pred$input_data$Person)), "NEW_HIGH")
  expect_equal(pv$population_review$RetainedPersons, 1)
  expect_equal(pv$population_review$OmittedPersons, 1)
  expect_equal(pv$population_review$OmittedRows, 2)
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
      fixture$fit, readiness_policy = "review",
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
    fixture$fit, readiness_policy = "review",
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
    fixture$fit, readiness_policy = "review",
    fixture$new_units,
    n_draws = 2,
    seed = 27
  )
  pv <- sample_mfrm_plausible_values(
    fixture$fit, readiness_policy = "review",
    fixture$new_units,
    n_draws = 2,
    seed = 28
  )

  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_identical(pred$settings$posterior_basis, "population_model")
  expect_equal(pred$population_review$RetainedPersons, 2)
  expect_equal(pred$population_review$OmittedPersons, 0)
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
  skip_if_not_installed("lpSolve")
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
      quad_points = 1,
      maxit = 50
    )
  )

  new_units <- data.frame(
    Person = c("NEW01", "NEW01"),
    Judge = unique(toy$Judge)[1],
    Task = unique(toy$Task)[1:2],
    Score = c(2, 3)
  )

  original_fit <- serialize(fit_jml, NULL)
  pred <- predict_mfrm_units(fit_jml, new_units, n_draws = 2, seed = 11)
  pv <- sample_mfrm_plausible_values(fit_jml, new_units, n_draws = 2, seed = 12)

  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_equal(pred$settings$method, "JML")
  expect_identical(pred$settings$scoring_quad_points, 31L)
  expect_true(all(pred$estimates$SD > 0))
  expect_true(any(grepl("standard normal reference prior", pred$notes, fixed = TRUE)))
  expect_true(all(c("Person", "Judge", "Task", "Score", "Weight") %in% names(pred$input_data)))
  expect_s3_class(pv, "mfrm_plausible_values")
  expect_equal(pv$settings$method, "JML")
  expect_true(any(grepl("standard normal reference prior", pv$notes, fixed = TRUE)))
  expect_equal(nrow(pv$values), 2)
  alternate <- sample_mfrm_plausible_values(fit_jml, new_units,
    scoring_prior = list(mean = 0.25, sd = 1.2), seed = 12)
  expect_true(all(alternate$estimates$PriorSource == "user_supplied"))
  expect_true(all(alternate$estimates$RetainedPriorSD == 1))
  expect_true(all(alternate$estimates$ScoreIntegrationReady))
  expect_false(any(grepl("standard normal reference prior", alternate$notes, fixed = TRUE)))
  expect_s3_class(summary(alternate), "summary.mfrm_plausible_values")
  expect_gt(max(abs(alternate$estimates$Estimate - pred$estimates$Estimate)), 1e-6)
  expect_true(any(grepl("post hoc fitted-object EAP", pred$notes, fixed = TRUE)))
  expect_identical(serialize(fit_jml, NULL), original_fit)
})

test_that("predict_mfrm_units supports GPCM fixed-calibration scoring", {
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

  pred <- predict_mfrm_units(
    fit_gpcm, new_units, readiness_policy = "review",
    n_draws = 2, seed = 31
  )
  pv <- sample_mfrm_plausible_values(
    fit_gpcm, new_units, readiness_policy = "review",
    n_draws = 2, seed = 32
  )

  expect_s3_class(pred, "mfrm_unit_prediction")
  expect_equal(pred$settings$method, "MML")
  expect_identical(pred$settings$posterior_basis, "population_model")
  expect_equal(sort(unique(pred$estimates$Person)), c("NEW01", "NEW02"))
  expect_true(all(is.finite(pred$estimates$Estimate)))
  expect_true(all(!pred$estimates$SourceScoringReady))
  expect_identical(unique(pred$estimates$EstimateUse),
                   "review_only_nonready_source")
  expect_true(any(grepl("conditional normal population model", pred$notes,
                        fixed = TRUE)))
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

test_that("prediction preprocessing audits missing responses and refuses invalid weights", {
  fixture <- make_prediction_fixture()
  new_units <- data.frame(
    Candidate = c("NEW01", "NEW01", "NEW02", "NEW03"),
    Judge = c(fixture$new_units$Rater[1], fixture$new_units$Rater[2], fixture$new_units$Rater[1], fixture$new_units$Rater[2]),
    Dimension = c(fixture$new_units$Criterion[1], fixture$new_units$Criterion[2], fixture$new_units$Criterion[1], fixture$new_units$Criterion[2]),
    Rating = c(2, 3, NA, 4),
    Wt = c(1, 1, 1, 1)
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
    "Dropped 1 row\\(s\\) from `new_data` before posterior scoring"
  )

  expect_equal(pred$row_review$InputRows, 4)
  expect_equal(pred$row_review$KeptRows, 3)
  expect_equal(pred$row_review$DroppedRows, 1)
  expect_equal(pred$row_review$DroppedMissing, 1)
  expect_equal(pred$row_review$DroppedNonpositiveWeight, 0)

  for (invalid in list(0, Inf, -Inf, NaN, NA_real_)) {
    bad <- new_units
    bad$Wt[1] <- invalid
    expect_error(
      predict_mfrm_units(
        fixture$fit,
        bad,
        person = "Candidate",
        facets = c(Rater = "Judge", Criterion = "Dimension"),
        score = "Rating",
        weight = "Wt"
      ),
      "Scoring weights must be finite and strictly positive",
      fixed = TRUE
    )
  }
})

test_that("unit scoring fails closed for non-ready fits unless review is explicit", {
  fixture <- make_prediction_fixture()
  nonready <- fixture$fit
  nonready$readiness$fit$FitReadiness[1] <- "blocked"
  nonready$readiness$fit$InferenceReady[1] <- FALSE
  nonready$readiness$fit$NumericalState[1] <- "failed"
  nonready$readiness$fit$ReasonCodes[1] <- "test_nonready_source"
  nonready$summary$FitReadiness[1] <- "blocked"
  nonready$summary$InferenceReady[1] <- FALSE
  nonready$summary$NumericalState[1] <- "failed"
  nonready$summary$ReadinessReasonCodes[1] <- "test_nonready_source"

  expect_error(
    predict_mfrm_units(nonready, fixture$new_units),
    "not ready for fitted-object scoring",
    fixed = TRUE
  )
  reviewed <- predict_mfrm_units(
    nonready, fixture$new_units, readiness_policy = "review"
  )
  expect_true(all(!reviewed$estimates$SourceScoringReady))
  expect_identical(reviewed$settings$source_scoring_status, "review_only")
  expect_true(any(grepl("review-only", reviewed$notes, fixed = TRUE)))
})

test_that("estimated-population scoring cannot bypass unresolved source validity", {
  dat <- data.frame(Person = paste0("P", 1:100), Rater = rep(c("R1", "R2"), each = 50),
    Score = c(rep(0, 35), rep(1, 15), rep(0, 15), rep(1, 35)))
  fit <- suppressWarnings(fit_mfrm(dat, "Person", "Rater", "Score", method = "MML",
    population_formula = ~ 1, person_data = dat["Person"], rating_min = 0,
    rating_max = 1, quad_points = 61, maxit = 300, reltol = 1e-10))
  new <- data.frame(Person = "NEW", Rater = "R1", Score = 1)
  expect_true(fit$data_review$estimability$nonlinear_local_estimability$local_first_order_rank_deficient)
  expect_error(predict_mfrm_units(fit, new), "estimated population distribution")
  expect_error(sample_mfrm_plausible_values(fit, new, n_draws = 3),
               "estimated population distribution")
  reviewed <- predict_mfrm_units(fit, new, readiness_policy = "review")
  pv <- sample_mfrm_plausible_values(fit, new, n_draws = 3, seed = 2,
                                    readiness_policy = "review")
  expect_false(reviewed$settings$source_scoring_ready)
  expect_false(any(reviewed$estimates$SourceScoringReady))
  expect_true(any(grepl("review-only", pv$notes, fixed = TRUE)))
  local_only <- fit
  local_only$readiness$fit$EstimabilityState <- "identified"
  expect_false(mfrmr:::prediction_source_scoring_readiness(local_only)$ready)

  for (x in list(reviewed, pv, summary(reviewed), summary(pv))) {
    stale <- x
    stale$settings$source_scoring_ready <- TRUE
    stale$settings$source_scoring_status <- "ready"
    stale$settings$readiness_policy <- "error"
    stale$estimates$SourceScoringReady <- TRUE
    if (inherits(stale, c("mfrm_unit_prediction", "mfrm_plausible_values"))) {
      expect_error(summary(stale), "requires explicit review-only eligibility")
      expect_error(mfrmr:::export_validate_optional_object(stale, class(stale)[1], "prediction", fit = fit),
                   "requires explicit review-only eligibility")
    }
    expect_error(print(stale), "requires explicit review-only eligibility")
    expect_error(build_summary_table_bundle(stale), "requires explicit review-only eligibility")
  }
  expect_s3_class(build_mfrm_manifest(fit, unit_prediction = reviewed, plausible_values = pv),
                  "mfrm_manifest")
  # A passed batch check cannot promote an unresolved population calibration.
  promoted <- reviewed
  promoted$settings$source_scoring_ready <- TRUE
  promoted$settings$source_scoring_status <- "ready"
  promoted$settings$readiness_policy <- "error"
  promoted$estimates$SourceScoringReady <- TRUE
  promoted$estimates$EstimateUse <- "fitted_object_scoring"
  expect_error(summary(promoted), "requires explicit review-only eligibility")
})

test_that("scoring quadrature is explicit and cannot degenerate to one point", {
  fixture <- make_prediction_fixture()
  expect_error(
    predict_mfrm_units(
      fixture$fit, fixture$new_units, scoring_quad_points = 1
    ),
    "greater than or equal to 2",
    fixed = TRUE
  )
  scored <- predict_mfrm_units(
    fixture$fit, fixture$new_units, scoring_quad_points = 17
  )
  expect_identical(scored$settings$quad_points, 17L)
  expect_identical(scored$settings$scoring_algorithm, "quadrature_eap_v2")
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

test_that("scoring tables retain interval, prior and draw meanings without internal console codes", {
  fixture <- make_prediction_fixture()
  level <- .81234567
  pred <- predict_mfrm_units(fixture$fit, fixture$new_units, interval_level = level,
                             n_draws = 3, seed = 47)
  table <- summary(pred, digits = 2)$estimates
  expect_identical(unique(table$IntervalLevel), level)
  expect_true(all(table$PriorMean == 0 & table$PriorSD == 1))
  expect_true(all(!table$WeightedLikelihood))
  expect_true(all(grepl("Continuous equal-tail", table$IntervalBasis)))
  expect_true(all(grepl("uncertainty excluded", table$UncertaintyBasis)))
  expect_true(all(grepl("Discrete quadrature", pred$draws$DrawBasis)))
  expect_identical(pred$draws$Person, rep(pred$estimates$Person, each = 3L))
  bundled <- build_summary_table_bundle(pred, digits = 2)$tables$estimates
  expect_identical(unique(bundled$IntervalLevel), level)
  text <- paste(capture.output(print(pred)), collapse = " ")
  expect_false(grepl("legacy_mml|quadrature_eap_v|SourceScoringReady|EstimateUse|readiness_contract", text))
  expect_match(text, "N\\(0,1\\)")

  weighted <- fixture$new_units
  weighted$Weight <- c(2, 1, 1, 1)
  scored <- predict_mfrm_units(fixture$fit, weighted, weight = "Weight")
  expect_identical(scored$estimates$WeightedLikelihood, c(TRUE, FALSE))
  expect_true(any(grepl("exponentiate likelihood", scored$notes)))

  # A known ordered draw sample distinguishes the requested central range
  # (9 to 90) from the former hard-coded 95% range (2 to 97).
  pv <- structure(list(values = data.frame(Person = pred$estimates$Person[1],
    Draw = 1:100, Value = 0:99), estimates = pred$estimates[1, ],
    settings = pred$settings, notes = pred$notes, row_review = pred$row_review),
    class = "mfrm_plausible_values")
  # Keep the artificial one-person draw sample's prior record aligned.
  pv$settings$prior_comparison <- pv$settings$prior_comparison[1, , drop = FALSE]
  if (!is.null(pv$settings$score_integration_review))
    pv$settings$score_integration_review <- pv$settings$score_integration_review[1, , drop = FALSE]
  draw_summary <- summary(pv, digits = 2)
  expect_equal(draw_summary$draw_summary$LowerValue, 9)
  expect_equal(draw_summary$draw_summary$UpperValue, 90)
  expect_identical(draw_summary$draw_summary$IntervalLevel, level)
  expect_match(draw_summary$draw_summary$DrawSummaryBasis, "Empirical quantiles")
  expect_identical(draw_summary$draw_summary$PriorMean, pv$estimates$PriorMean)
  expect_match(draw_summary$draw_summary$UncertaintyBasis, "estimation uncertainty excluded")
  expect_match(paste(capture.output(print(draw_summary)), collapse = " "), "With few draws")
  stale <- draw_summary
  stale$draw_summary$DrawSummaryBasis <- NULL
  expect_error(print(stale), "no refitting or resampling")
  expect_error(build_summary_table_bundle(stale), "no refitting or resampling")

  older <- pred
  # Reconstruct the old schema, which predates separate retained-prior records.
  older$settings$prior_comparison <- NULL
  older$settings$retained_posterior_basis <- NULL
  older$settings$score_integration_required <- NULL
  older$settings$score_integration_review <- NULL
  older$estimates <- older$estimates[c("Person", "Estimate", "SD", "Lower", "Upper", "Observations", "WeightedN", "SourceScoringReady", "EstimateUse")]
  expect_no_warning(summary(older))
  expect_identical(unique(summary(older, digits = 0)$estimates$IntervalLevel), level)
  expect_true(all(summary(older)$estimates$PriorSD == 1))
  missing_method <- older
  missing_method$settings$scoring_algorithm <- NULL
  expect_true(all(grepl("not recorded", summary(missing_method)$estimates$IntervalBasis)))
  expect_error(predict_mfrm_units(fixture$fit, fixture$new_units, interval_level = c(.8, .95)), "single number")
  imported <- structure(list(), class = c("mfrm_imported_fit", "mfrm_fit"))
  expect_error(predict_mfrm_units(imported, fixture$new_units, readiness_policy = "review"), "native fit_mfrm")
})
