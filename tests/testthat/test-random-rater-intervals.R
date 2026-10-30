bootstrap_interval_fixture <- function() {
  source <- list(raters = data.frame(Rater = c("B", "A"), Estimate = c(1, -2), PredictionSE = c(2, .5)))
  roots <- cbind(B = seq(-2, 2, length.out = 20), A = seq(-1, 3, length.out = 20))
  structure(list(source = source, studentized = roots, error = roots * 2,
    settings = list(level = .8, nsim = 20), trials = data.frame(FitReady = rep(TRUE, 20),
      EstimatedBoundary = rep(FALSE, 20))), class = "mfrm_random_rater_intervals")
}

test_that("bootstrap endpoints use truth-minus-estimate direction and planned roots", {
  x <- bootstrap_interval_fixture()
  ci <- confint(x, level = .8)
  q <- apply(x$studentized, 2, quantile, probs = c(.1, .9), type = 1, names = FALSE)
  expected <- cbind(Lower = c(1, -2) + c(2, .5) * q[1, ], Upper = c(1, -2) + c(2, .5) * q[2, ])
  expect_equal(unname(ci[, , drop = FALSE]), unname(expected), ignore_attr = TRUE)
  expect_equal(attr(ci, "expected_tail_draws"), 2)
  expect_identical(rownames(confint(x, parm = c("A", "B"))), c("A", "B"))
  expect_error(confint(x, parm = "missing"), "IDs")
  expect_error(confint(x, level = NA), "level")
  x$studentized[c(1, 6), 1] <- NA
  bound <- confint(x, level = .8)
  expect_identical(bound[1, 1], -Inf)
  expect_true(attr(bound, "availability")$Unresolved[1] == 2)
  # Enumerate extreme finite completions: the retained envelope contains them all.
  for (a in c(-100, 0, 100)) for (b in c(-100, 0, 100)) {
    completed <- x; completed$studentized[c(1, 6), 1] <- c(a, b)
    inner <- confint(completed, level = .8)
    expect_lte(bound[1, 1], inner[1, 1]); expect_gte(bound[1, 2], inner[1, 2])
  }
  raw <- confint(x, level = .8, method = "error")
  expect_equal(attr(raw, "availability")$Known, c(20L, 20L))
  x$studentized[,] <- NA
  unavailable <- confint(x)
  expect_true(all(unavailable[, 1] == -Inf & unavailable[, 2] == Inf))
})

test_that("bootstrap plots keep unbounded endpoints and source labels", {
  x <- bootstrap_interval_fixture(); x$studentized[1:4, 1] <- NA
  x$intervals <- confint(x, level = .8)
  before <- grDevices::dev.cur(); p <- plot(x, draw = FALSE)
  expect_identical(grDevices::dev.cur(), before)
  tab <- plot_data(p)$table
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(as_ggplot(p), "ggplot")
  expect_identical(tab$Rater, c("B", "A"))
  expect_identical(tab$Lower[1], -Inf); expect_identical(tab$Upper[1], Inf)
  level <- 1 - .Machine$double.eps / 2
  high <- plot_data(plot(x, level = level, draw = FALSE))$table
  z <- sqrt(qchisq(level, 1))
  expect_equal(high$OrdinaryLower, x$source$raters$Estimate-z*x$source$raters$PredictionSE)
  expect_equal(high$OrdinaryUpper, x$source$raters$Estimate+z*x$source$raters$PredictionSE)
  expect_identical(high$Lower[1], -Inf)
  expect_identical(high$Upper[1], Inf)
  grDevices::pdf(tempfile(fileext = ".pdf")); old <- graphics::par("mar")
  expect_silent(plot(x)); expect_identical(graphics::par("mar"), old); grDevices::dev.off()
  path <- tempfile(); saveRDS(x, path)
  expect_identical(confint(readRDS(path)), confint(x))
})

test_that("boundary and failed refits retain their prediction-error availability", {
  object <- bootstrap_interval_fixture()$source
  object$input <- list(columns = list(person = "Person", rater = "Rater", score = "Score",
    facets = NULL), score_levels = 0:1)
  object$settings <- list(fixed_rater_sd = NULL, quad_points = 31L, maxit = 300L)
  local_mocked_bindings(mfrm_random_rater_generate = function(...) list(
    data = data.frame(Score = 0), truth = c(A = -.4, B = .7)))
  boundary <- list(checks = list(NumericalReady = TRUE, InformationPositive = TRUE,
    EstimatedVarianceBoundary = TRUE, MaxGradient = 0), calibration = list(rater_sd = 0),
    raters = data.frame(Rater = c("A", "B"), Estimate = c(0, 0), PredictionSE = c(NA_real_, NA_real_)))
  local_mocked_bindings(fit_mfrm_random_rater = function(...) boundary)
  z <- mfrm_random_rater_bootstrap_one(object, 7L)
  expect_true(z$trial$FitReady); expect_true(z$trial$EstimatedBoundary)
  expect_equal(z$error, c(.7, -.4)); expect_true(all(is.na(z$studentized)))
  local_mocked_bindings(fit_mfrm_random_rater = function(...) {
    warning("unresolved curvature"); stop("refit failed")
  })
  failed <- mfrm_random_rater_bootstrap_one(object, 8L)
  expect_false(failed$trial$FitReady)
  expect_true(all(is.na(failed$error))); expect_true(all(is.na(failed$studentized)))
  expect_equal(failed$truth, c(.7, -.4))
  expect_identical(failed$trial$Error, "refit failed")
  expect_identical(failed$trial$Warnings, "unresolved curvature")
  expect_true(is.na(failed$trial$NumericalReady))
  expect_true(is.na(failed$trial$InformationPositive))
  expect_true(is.na(failed$trial$OptimizerCode))
})

test_that("bootstrap records distinguish numerical and information failures without changing roots", {
  object <- bootstrap_interval_fixture()$source
  object$input <- list(columns = list(person = "Person", rater = "Rater", score = "Score",
    facets = NULL), score_levels = 0:1)
  object$settings <- list(fixed_rater_sd = NULL, quad_points = 61L, maxit = 300L)
  local_mocked_bindings(mfrm_random_rater_generate = function(...) list(
    data = data.frame(Score = 0), truth = c(A = -.4, B = .7)))
  candidate <- list(checks = data.frame(OptimizerCode = 0L, NumericalReady = FALSE,
    InformationPositive = TRUE, EstimatedVarianceBoundary = FALSE,
    EstimatedPersonVarianceBoundary = FALSE, PersonVarianceUpperBoundary = FALSE,
    PersonQuadratureStable = FALSE,
    MaxGradient = 1e-8, QuadraturePoints = 61L, CheckPoints = 123L,
    LogLikDifference = 2e-5, GradientDifference = 2e-4),
    calibration = list(rater_sd = .4, person_sd = 1.3),
    raters = data.frame(Rater = c("B", "A"), Estimate = c(.5, -.3), PredictionSE = c(.2, .1)))
  local_mocked_bindings(fit_mfrm_random_rater = function(...) candidate)
  numeric_failure <- mfrm_random_rater_bootstrap_one(object, 7L)
  expect_false(numeric_failure$trial$FitReady)
  expect_identical(numeric_failure$trial[names(candidate$checks)], candidate$checks)
  expect_true(all(is.na(numeric_failure$error)))
  expect_true(all(is.na(numeric_failure$studentized)))
  candidate$checks$NumericalReady <- TRUE
  candidate$checks$PersonQuadratureStable <- TRUE
  candidate$checks$InformationPositive <- FALSE
  candidate$checks$LogLikDifference <- 1e-8
  candidate$checks$GradientDifference <- 1e-7
  information_failure <- mfrm_random_rater_bootstrap_one(object, 7L)
  expect_false(information_failure$trial$FitReady)
  expect_identical(information_failure$trial[names(candidate$checks)], candidate$checks)
  expect_true(all(is.na(information_failure$studentized)))
  candidate$checks$InformationPositive <- TRUE
  ready <- mfrm_random_rater_bootstrap_one(object, 7L)
  expect_true(ready$trial$FitReady)
  expect_equal(ready$error, c(.2, -.1))
  expect_equal(ready$studentized, c(1, -1))
  path <- tempfile(fileext = ".rds"); withr::defer(unlink(path))
  saveRDS(list(numeric_failure, information_failure, ready), path)
  expect_identical(readRDS(path), list(numeric_failure, information_failure, ready))
})

test_that("ordinal generation preserves fixed facets, repeated identities and category labels", {
  person <- c(2L, 1L, 2L, 3L, 1L, 2L, 3L, 3L)
  rater <- c(1L, 2L, 2L, 1L, 1L, 1L, 2L, 1L)
  X <- cbind(c(-1, 0, 1, -1, 1, 0, 0, 1), c(0, 1, -1, 1, 0, -1, 1, 0))
  object <- list(input = list(y = rep(0, 8), person = person, rater = rater, X = X,
    data = data.frame(Participant = person, Judge = rater, Rating = NA_integer_),
    columns = list(rater = "Judge", score = "Rating"),
    levels = list(Judge = c("z", "a")), score_levels = c(1, 3, 7)),
    calibration = list(rater_sd = .7, beta = c(.2, -.4), steps = c(-.5, .8)))
  before <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
  withr::defer(if (is.null(before)) rm(".Random.seed", envir = .GlobalEnv) else
    assign(".Random.seed", before, envir = .GlobalEnv))
  set.seed(907); theta <- rnorm(3); u <- rnorm(2, sd = .7); p <- runif(8)
  expected <- vapply(seq_len(8), function(i) {
    eta <- theta[person[i]] - u[rater[i]] - sum(X[i, ] * c(.2, -.4))
    adjacent <- exp(eta - c(-.5, .8))
    weights <- c(1, adjacent[1], prod(adjacent))
    c(1, 3, 7)[which(p[i] <= cumsum(weights / sum(weights)))[1]]
  }, numeric(1))
  generated <- mfrm_random_rater_generate(object, 907)
  expect_identical(generated$data$Rating, expected)
  expect_identical(generated$truth, setNames(u, c("z", "a")))
  expect_identical(generated$data[c("Participant", "Judge")], object$input$data[c("Participant", "Judge")])
})

test_that("shared generation, refits and saved roots preserve the experiment", {
  skip_if_not_installed("RTMB", "2.0")
  # No fixtures from ignored validation artifacts; a reproducible binary design.
  d <- expand.grid(Person = 1:24, Rater = c("R:1", "R:2", "R:3"))
  d$Score <- (d$Person + match(d$Rater, unique(d$Rater))) %% 2
  fit <- fit_mfrm_random_rater(d, "Person", "Rater", "Score", score_levels = 0:1, rater_sd = .4, person_sd = 1)
  set.seed(371); before <- .Random.seed
  a <- mfrm_random_rater_generate(fit, 178)
  expect_identical(.Random.seed, before)
  set.seed(178); theta <- rnorm(24); u <- rnorm(3, sd = .4); uniforms <- runif(nrow(d))
  eta <- theta[fit$input$person] - u[fit$input$rater] - fit$calibration$steps
  expected <- as.integer(uniforms > 1 / (1 + exp(eta)))
  expect_identical(a$data$Score, expected)
  expect_equal(unname(a$truth), u)
  expect_equal(a$data[c("Person", "Rater")], fit$input$data[c("Person", "Rater")])
  before <- .Random.seed
  b <- mfrm_random_rater_intervals(fit, nsim = 3, seed = 924013)
  expect_identical(.Random.seed, before)
  expect_identical(dim(b$studentized), c(3L, 3L))
  expect_identical(colnames(b$studentized), fit$raters$Rater)
  expect_true(all(b$trials$FitReady))
  expect_true(all(b$trials$RaterSD == .4))
  expect_equal(b$error, b$truth - b$estimates)
  expect_equal(b$studentized, b$error / b$prediction_se)
  reproduced <- mfrm_random_rater_bootstrap_one(fit, b$trials$Seed[2])
  expect_equal(unname(b$error[2, ]), reproduced$error)
  # RNGkind() and native code may initialize a previously absent RNG state.
  old <- .Random.seed
  withr::defer(assign(".Random.seed", old, envir = .GlobalEnv))
  rm(".Random.seed", envir = .GlobalEnv)
  no_state <- mfrm_random_rater_intervals(fit, nsim = 2, seed = 924013)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_equal(no_state$error, b$error[1:2, ], ignore_attr = TRUE)
  expect_error(mfrm_random_rater_intervals(fit, nsim = 1, seed = 1), "nsim")
  expect_error(mfrm_random_rater_intervals(fit, nsim = 2), "seed")
  fit$calibration$rater_sd <- 0
  expect_error(mfrm_random_rater_intervals(fit, nsim = 2, seed = 1), "boundary")
})
