population_terms_fixture <- function(formula) {
  persons <- data.frame(Person = paste0("P", 1:12), X = seq(-2, 3, length.out = 12),
    Group = rep(c("a", "b"), 6))
  pop <- prepare_mfrm_population_scaffold(persons, "Person", formula,
    persons, "Person")
  pop$coefficients <- setNames(seq_len(ncol(pop$design_matrix)) / 5, pop$design_columns)
  pop$sigma2 <- .7
  list(population = pop, config = list(posterior_basis = "population_model",
    population_spec = compact_population_spec(pop, persons$Person)))
}

population_terms_prediction <- function(fit, persons, policy = "error") {
  prepared <- list(input_data = persons,
    prep = list(levels = list(Person = persons$Person)))
  prepare_mfrm_prediction_population(fit, prepared, persons, "Person", policy)$scaffold
}

test_that("latent-regression transforms retain their training basis in new cohorts", {
  future <- data.frame(Person = paste0("NEW", 1:3), X = c(4, 6, 9), Group = c("b", "a", "b"))
  for (formula in list(~ scale(X), ~ poly(X, 2), ~ splines::ns(X, df = 3), ~ scale(X) * Group)) {
    fit <- population_terms_fixture(formula)
    train <- fit$population$person_table
    mf <- stats::model.frame(formula, train)
    training_terms <- attr(mf, "terms")
    expected <- stats::model.matrix(training_terms,
      stats::model.frame(training_terms, future, xlev = fit$population$xlevels),
      contrasts.arg = fit$population$contrasts)
    batch <- population_terms_prediction(fit, future)
    expect_equal(unname(batch$design_matrix), unname(expected), tolerance = 1e-13,
      ignore_attr = TRUE)
    one <- population_terms_prediction(fit, future[1, ])
    expect_equal(unname(one$design_matrix), unname(expected[1, , drop = FALSE]), tolerance = 1e-13)
    replay <- unserialize(serialize(fit, NULL))
    expect_equal(population_terms_prediction(replay, future)$design_matrix,
      batch$design_matrix, tolerance = 0)
    expect_equal(drop(batch$design_matrix %*% fit$population$coefficients),
      drop(expected %*% fit$population$coefficients), tolerance = 1e-13)
  }
})

test_that("legacy population transforms are recovered only from matching training data", {
  fit <- population_terms_fixture(~ scale(X))
  future <- data.frame(Person = c("NEW1", "NEW2"), X = c(4, 7))
  expected <- (future$X - mean(fit$population$person_table$X)) / sd(fit$population$person_table$X)
  fit$population$terms <- fit$config$population_spec$terms <- NULL
  expect_equal(unname(population_terms_prediction(fit, future)$design_matrix[, 2]), expected)
  changed <- fit
  changed$population$person_table_replay$X[1] <- 100
  expect_error(population_terms_prediction(changed, future), "training design")
  missing <- fit
  missing$population$person_table <- missing$population$person_table_replay <- NULL
  expect_error(population_terms_prediction(missing, future), "training.*terms|terms.*training")
  plain <- population_terms_fixture(~ X + Group)
  plain$population$terms <- plain$config$population_spec$terms <- NULL
  plain$population$person_table <- plain$population$person_table_replay <- NULL
  expect_equal(unname(population_terms_prediction(plain,
    data.frame(Person = "NEW", X = 4, Group = "b"))$design_matrix), matrix(c(1, 4, 1), 1))
})

test_that("population terms preserve fitted scaling across complete-case omission", {
  persons <- data.frame(Person = paste0("P", 1:8), X = 1:8, Z = c(NA, rep(c(0, 1), 3), 0))
  pop <- prepare_mfrm_population_scaffold(persons, "Person", ~ scale(X) + Z,
    persons, "Person", population_policy = "omit")
  expect_s3_class(pop$terms, "terms")
  spec <- compact_population_spec(pop, persons$Person)
  expect_identical(spec$terms, pop$terms)
  fit <- list(population = pop, config = list(posterior_basis = "population_model", population_spec = spec))
  future <- data.frame(Person = "NEW", X = 9, Z = 0)
  expected <- (9 - mean(persons$X)) / sd(persons$X)
  expect_equal(unname(population_terms_prediction(fit, future)$design_matrix[, 2]), expected)
  fit$population$terms <- fit$config$population_spec$terms <- NULL
  expect_equal(unname(population_terms_prediction(fit, future)$design_matrix[, 2]), expected)
})

test_that("public latent-regression scoring is unchanged by the other scored Persons", {
  data <- simulate_mfrm_data(n_person = 24, n_rater = 3, n_criterion = 2,
    score_levels = 3, model = "RSM", seed = 100226)
  persons <- data.frame(Person = unique(data$Person), X = seq(-2, 3, length.out = 24))
  fit <- fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score",
    model = "RSM", method = "MML", population_formula = ~ scale(X),
    person_data = persons, quad_points = 61)
  expect_s3_class(fit$population$terms, "terms")
  future <- persons[1:3, ]
  future$X <- c(4, 6, 9)
  responses <- data[data$Person %in% future$Person, ]
  batch <- predict_mfrm_units(fit, responses, person_data = future,
    scoring_quad_points = 121, readiness_policy = "review")
  single <- predict_mfrm_units(fit, responses[responses$Person == future$Person[1], ],
    person_data = future[1, ], scoring_quad_points = 121, readiness_policy = "review")
  fields <- c("Estimate", "SD", "Lower", "Upper")
  expect_true(all(is.finite(as.matrix(single$estimates[fields]))))
  expect_equal(single$estimates[fields], batch$estimates[batch$estimates$Person == future$Person[1], fields],
    tolerance = 1e-12, ignore_attr = TRUE)
  legacy <- unserialize(serialize(fit, NULL))
  legacy$population$terms <- legacy$config$population_spec$terms <- NULL
  replay <- predict_mfrm_units(legacy, responses, person_data = future,
    scoring_quad_points = 121, readiness_policy = "review")
  expect_equal(replay$estimates[fields], batch$estimates[fields], tolerance = 1e-12)
  expect_identical(replay$settings$source_scoring_ready, batch$settings$source_scoring_ready)
})
