# Small display fixtures: no fitted statistical conclusions are tested here.
extended_report_fixture <- function(testlet = TRUE) {
  data <- data.frame(Person = c("P1", "P1", "P2", "P2"), Rater = c("A", "B", "A", "B"), Score = c(0, 1, 1, 0))
  input <- list(data = data, assigned_data = data, input_rows = 5L, omitted_rows = 5L,
    columns = list(person = "Person", score = "Score", facets = if (testlet) "Rater" else character()),
    levels = list(Rater = c("A", "B")), basis = list(Rater = matrix(c(-1, 1), 2, dimnames = list(c("A", "B"), "B1"))),
    score_levels = 0:1, blocks = data.frame(Person = c("P1", "P2"), Testlet = "T", Observed = 2L))
  input$columns[[if (testlet) "testlet" else "rater"]] <- "Rater"
  structure(list(input = input, parameters = c(0, 0, 0),
    settings = list(method = "Fixture", model = "RSM", person_distribution = "N(0,1)",
      missing = "omit", quad_points = 7L),
    checks = list(NumericalReady = TRUE, InformationPositive = TRUE, EstimatedVarianceBoundary = FALSE),
    calibration = list(beta = if (testlet) 0 else numeric(), steps = 0, variance = .2, rater_sd = .7),
    calibration_table = data.frame(Parameter = "Fixed facet", Facet = "Rater", Level = c("A", "B"),
      Estimate = c(-.3, .3), SE = .2, Lower = c(-.7, -.1), Upper = c(.1, .7)),
    raters = data.frame(Rater = c("A", "B"), Persons = 2, Estimate = c(-.3, .3),
      ConditionalSD = .2, PredictionSE = .25, Lower = c(-.8, -.2), Upper = c(.2, .8)),
    rater_mode = c(-.3, .3), conditional_rater_covariance = diag(.04, 2)),
    class = if (testlet) "mfrm_testlet" else "mfrm_random_rater")
}

extended_scores_fixture <- function(fit) {
  structure(list(table = data.frame(Person = c("P1", "Missing", "Unresolved"),
    Estimate = c(.4, 0, NA_real_), Lower = c(-.5, -1.96, NA_real_), Upper = c(1.2, 1.96, NA_real_),
    Status = c("available_conditional", "prior_only", "unavailable"), Reason = c("", "", "Integration unresolved")),
    source = mfrm_extended_prediction_source(fit), blocks = fit$input$blocks,
    data_usage = c(Input = 6L, Observed = 4L, Omitted = 2L), omitted_rows = 5:6,
    settings = list(level = .95, calibration_uncertainty = FALSE)), class = "mfrm_testlet_scores")
}

test_that("extended reporting collects saved values without model computations", {
  fit <- extended_report_fixture(); scores <- extended_scores_fixture(fit)
  fail <- function(...) stop("Unexpected computation")
  local_mocked_bindings(fit_mfrm_testlet = fail, fit_mfrm_random_rater = fail,
    predict.mfrm_testlet = fail, predict.mfrm_random_rater = fail, diagnose_mfrm = fail,
    mfrm_random_rater_intervals = fail)
  res <- mfrm_results(fit, predictions = scores, include = "all")
  expect_identical(res, mfrm_results(fit, predictions = scores, include = "all", compute = "never"))
  expect_identical(res$tables$person_scores, scores$table)
  expect_identical(res$tables$scoring_omitted_rows$InputRow, 5:6)
  expected <- fit$calibration_table; expected$Lower <- expected$Upper <- NA_real_
  expect_identical(res$tables$calibration[names(expected)], expected)
  expect_true(all(res$tables$calibration$IntervalStatus == "not_requested"))
  expect_false(any(res$tables$calibration$CIEligible))
  expected <- res$tables$calibration
  expect_identical(res$tables$data_usage$Omitted, 1L)
  expect_true(all(res$status$Status[res$status$Section %in% c("diagnostics", "bias", "linking", "response_time", "apa")] == "not_available"))
  expect_match(summary(res)$decision$Interpretation, "model adequacy not assessed")
  for (style in c("qc", "apa", "validation", "reviewer", "technical")) {
    report <- mfrm_report(res, style = style)
    expect_identical(report$tables$person_scores, scores$table)
    expect_identical(summary(report)$overview$OverallStatus, "caveat")
    expect_match(report$markdown, "prior_only")
    expect_match(report$markdown, "Integration unresolved")
    expect_match(report$markdown, "not implemented")
    expect_identical(mfrm_report(res, style, output = "tables"), report$tables)
  }
  expect_identical(mfrm_results(fit, predictions = scores, output = "tables"), mfrm_results(fit, predictions = scores)$tables)
  expect_s3_class(mfrm_results(fit, output = "summary"), "summary.mfrm_results")
  expect_s3_class(mfrm_results(fit, output = "html"), "mfrm_results_html")
  expect_s3_class(mfrm_report(res, output = "html"), "mfrm_report_html")
  expect_output(print(res), "calibration")
  expect_identical(plot_data(plot(res, type = "scores", draw = FALSE))$table, scores$table)
  expect_identical(plot_data(plot(res, type = "calibration", draw = FALSE))$table, expected)
  expect_error(plot(res, type = "wright"), "available type")
  expect_error(plot(mfrm_results(fit, include = "fit")), "available type")
  expect_error(launch_mfrmr_viewer(res, return_app = TRUE), "Shiny viewer does not support")
})

test_that("source matching rejects different calibration semantics and unsupported attachments", {
  fit <- extended_report_fixture(); scores <- extended_scores_fixture(fit)
  old <- scores; old$source <- NULL
  expect_error(mfrm_results(fit, predictions = old), "Regenerate older predictions")
  for (field in c("calibration", "columns", "levels", "basis", "settings", "checks")) {
    changed <- fit
    if (field %in% c("columns", "levels", "basis")) changed$input[[field]] <- list(changed = TRUE) else changed[[field]] <- list(changed = TRUE)
    expect_error(mfrm_results(changed, predictions = scores), "matching source metadata")
  }
  expect_error(mfrm_results(fit, predictions = list()), "matching source metadata")
  expect_error(mfrm_results(fit, intervals = list()), "exact random-rater fit")
  expect_error(mfrm_results(fit, diagnostics = list()), "mfrm_response_diagnostics")
  expect_error(mfrm_results(fit, response_time_score = "Score"), "unavailable")
  expect_error(mfrm_results(structure(list(), class = "mfrm_fit"), predictions = scores), "only for testlet")
  expect_error(mfrm_results(structure(list(), class = "mfrm_fit"), intervals = list()), "native GPCM fit")
  expect_error(mfrm_results(scores), "first argument")
  # Numerical failures remain reportable; no successful-inference label is introduced.
  fit$checks$NumericalReady <- FALSE
  report <- mfrm_report(mfrm_results(fit))
  expect_identical(summary(report)$overview$OverallStatus, "review")
  expect_match(report$decision$Interpretation, "Numerical review")
  fit$checks$NumericalReady <- TRUE; fit$checks$EstimatedVarianceBoundary <- TRUE
  fit$calibration_table$Lower <- fit$calibration_table$Upper <- NA_real_
  res <- mfrm_results(fit)
  expect_match(res$decision$Interpretation, "boundary")
  expect_true(all(is.na(res$tables$calibration$Lower)))
  expect_true(subset(res$tables$variance, Effect == "Person-specific testlet")$EstimatedBoundary)
  expect_false(subset(res$tables$variance, Effect == "Person ability")$EstimatedBoundary)
})

test_that("prediction calls attach calibration semantics without retaining training rows", {
  fit <- extended_report_fixture(FALSE)
  p <- predict(fit, data.frame(Rater = c("A", "B")), ability = c(-1, 1))
  expect_identical(p$source, mfrm_extended_prediction_source(fit))
  expect_null(p$source$input)
  res <- mfrm_results(fit, predictions = p)
  expect_identical(res$tables$category_probabilities$Probability, as.vector(p$probabilities))
  expect_identical(res$tables$expected_scores, p$expected_scores)
  expect_identical(res$tables$prediction_inputs, p$newdata)
  changed <- fit; changed$rater_mode[1] <- 99
  expect_error(mfrm_results(changed, predictions = p), "matching source metadata")
  # With no scored rows the testlet scorer returns the prior and adds provenance.
  fit <- extended_report_fixture()
  new <- data.frame(Person = "Prior", Rater = "A", Score = NA_real_)
  p <- predict(fit, new, missing = "omit")
  expect_identical(p$source, mfrm_extended_prediction_source(fit))
  expect_identical(p$table$Status, "prior_only")
  expect_identical(mfrm_results(fit, predictions = p)$tables$person_scores, p$table)
})

test_that("random-rater report retains every bootstrap trial and infinite endpoints", {
  fit <- extended_report_fixture(FALSE)
  roots <- matrix(seq(-2, 2, length.out = 40), 20, 2)
  roots[1:4, 1] <- NA_real_
  intervals <- structure(list(source = fit, studentized = roots, error = roots,
    settings = list(level = .8, nsim = 20),
    trials = data.frame(Replicate = 1:20, FitReady = c(rep(FALSE, 4), rep(TRUE, 16)))), class = "mfrm_random_rater_intervals")
  intervals$intervals <- confint(intervals, level = .8)
  res <- mfrm_results(fit, intervals = intervals)
  expect_identical(res$tables$bootstrap_trials, intervals$trials)
  expect_identical(res$tables$bootstrap_intervals$Lower[1], -Inf)
  expect_identical(res$tables$bootstrap_intervals$Upper[1], Inf)
  expect_identical(res$tables$bootstrap_availability, attr(intervals$intervals, "availability"))
  expect_match(tail(res$tables$interval_basis$Interval, 1), "80%")
  expect_identical(plot_data(plot(res, type = "intervals", draw = FALSE))$table$Lower,
    res$tables$bootstrap_intervals$Lower)
  expected <- fit$raters; expected$Lower <- expected$Upper <- NA_real_
  expect_identical(res$tables$raters[names(expected)], expected)
  expect_true(all(res$tables$raters$IntervalStatus == "not_requested"))
  expect_false(any(res$tables$raters$CIEligible))
  expected <- res$tables$raters
  expect_identical(plot_data(plot(res, draw = FALSE))$table, expected)
  expect_match(res$tables$interval_basis$Interval[3], "No automatic")
  bad <- intervals; bad$source$raters$Estimate[1] <- 77
  expect_error(mfrm_results(fit, intervals = bad), "exact random-rater fit")
  out <- export_mfrm_results(res, tempfile(), include = c("tables", "report"), acknowledge_sensitive = TRUE)
  path <- out$written_files$Path[out$written_files$Component == "table_bootstrap_intervals"]
  expect_identical(read.csv(path)$Lower[1], -Inf)
  expect_identical(mfrm_report(res)$tables$bootstrap_trials$FitReady, intervals$trials$FitReady)
})

test_that("extended archives replay without fitting and retain missing and unavailable rows", {
  fit <- extended_report_fixture(); scores <- extended_scores_fixture(fit)
  res <- mfrm_results(fit, predictions = scores)
  dir <- tempfile(); out <- export_mfrm_results(res, dir, prefix = "testlet", preset = "starter", acknowledge_sensitive = TRUE)
  expect_equal(nrow(out$plot_errors), 0L)
  expect_identical(readRDS(file.path(dir, "testlet_results.rds")), res)
  index <- paste(readLines(file.path(dir, "index.html")), collapse = "\n")
  expect_false(grepl("required Wright|rerun the export|Infit-versus-measure", index))
  expect_match(index, "testlet_plot_scores.png", fixed = TRUE)
  expect_match(index, "descriptive posterior residuals without ordinary-model cutoffs")
  expect_match(index, "Only saved Person conditional intervals")
  csv <- out$written_files$Path[out$written_files$Component == "table_person_scores"]
  expect_identical(read.csv(csv)$Person, scores$table$Person)
  expect_true(is.na(read.csv(csv)$Estimate[3]))
  fail <- function(...) stop("Unexpected recomputation")
  local_mocked_bindings(fit_mfrm_testlet = fail, predict.mfrm_testlet = fail, diagnose_mfrm = fail)
  withr::local_dir(dir)
  env <- new.env()
  capture.output(sys.source("testlet_replay.R", envir = env))
  expect_identical(env$res, res)
  expect_identical(env$report$tables, res$tables)
  out2 <- export_mfrm_results(res, tempfile(), include = "replay", acknowledge_sensitive = TRUE)
  expect_true(all(c("rds", "R") %in% out2$written_files$Format))
})


test_that("explicit calibration intervals agree across saved-fit output routes", {
  for (testlet in c(TRUE, FALSE)) {
    fit <- extended_report_fixture(testlet)
    fit$calibration_table <- rbind(fit$calibration_table,
      data.frame(Parameter = c("Step", "Population SD"), Facet = c("Step", "Person"),
        Level = c("1", "SD"), Estimate = c(.4, 1), SE = c(.1, .2), Lower = -999, Upper = 999))
    original <- fit
    expect_true(all(is.na(summary(fit)$calibration$Lower)))
    expect_output(print(fit), "Calibration bounds omitted")
    for (level in c(.8, .925, .95, 1 - .Machine$double.eps / 2)) {
      expected <- fit$calibration_table
      expected$Lower <- expected$Upper <- NA_real_
      z <- qnorm((1 - level) / 2, lower.tail = FALSE)
      expected$Lower[1:3] <- expected$Estimate[1:3] - z * expected$SE[1:3]
      expected$Upper[1:3] <- expected$Estimate[1:3] + z * expected$SE[1:3]
      ci <- confint(fit, parm = "calibration", level = level)
      expect_equal(as.numeric(ci), as.numeric(as.matrix(expected[1:3, c("Lower", "Upper")])))
      expect_identical(attr(ci, "level"), level)
      expect_match(attr(ci, "note"), "nominal coverage is not established")
      actual <- summary(fit, calibration_intervals = "normal", level = level)$calibration
      expect_identical(actual[names(expected)], expected)
      expect_identical(actual$IntervalStatus, c(rep("available", 3), "unsupported_target"))
      expect_identical(actual$CIEligible, c(TRUE, TRUE, TRUE, FALSE))
      expect_identical(attr(ci, "diagnostics"), actual[1:3, ])
      expected <- actual
      res <- mfrm_results(fit, calibration_intervals = "normal", calibration_level = level)
      expect_identical(res$tables$calibration, expected)
      expect_identical(mfrm_report(res)$tables$calibration, expected)
      if (testlet) {
        p <- plot(res, type = "calibration", draw = FALSE)
        expect_identical(p$data$table, expected[1:2, ])
        expect_identical(p$data$settings$level, level)
        expect_identical(p$data$settings$calibration_intervals, "normal")
      }
    }
    expect_identical(fit, original)
    expect_error(summary(fit, intervals = "normal"), "empty")
    for (level in list(0, 1, NA_real_, c(.8, .9), "95%")) {
      expect_error(confint(fit, parm = "calibration", level = level), "0 < level < 1")
      expect_error(mfrm_results(fit, calibration_level = level), "0 < level < 1")
    }
    for (flag in c("NumericalReady", "InformationPositive", "EstimatedVarianceBoundary", "EstimatedPersonVarianceBoundary")) {
      invalid <- fit
      invalid$checks[[flag]] <- !flag %in% c("NumericalReady", "InformationPositive")
      expect_true(all(is.na(confint(invalid, parm = "calibration"))))
    }
    invalid <- fit; invalid$calibration_table$SE <- c(NA, Inf, -1, .2)
    expect_true(all(is.na(confint(invalid, parm = "calibration"))))
  }
  fit <- extended_report_fixture()
  expect_error(confint(fit, parm = "variance"), "not supplied")
  expect_error(plot(fit, style = "precision", draw = FALSE), "explicit intervals")
  expect_error(plot(fit, sort = "uncertainty", draw = FALSE), "explicit intervals")
  expect_true(all(is.na(plot(fit, draw = FALSE)$data$table$Lower)))
  expect_s3_class(plot(fit, style = "precision", intervals = "normal", draw = FALSE), "mfrm_plot_data")
  p <- plot(fit, intervals = "normal", level = .925, draw = FALSE)
  expect_match(p$data$notes$Text, "92.5%", fixed = TRUE)
  ordinary <- structure(list(), class = "mfrm_fit")
  expect_error(mfrm_results(ordinary, calibration_intervals = "normal"), "only for testlet")
  expect_error(mfrm_results(ordinary, calibration_level = .9), "only for testlet")
})

test_that("explicit calibration selection survives archive and replay without fitting", {
  fit <- extended_report_fixture()
  res <- mfrm_results(fit, calibration_intervals = "normal", calibration_level = .8)
  fail <- function(...) stop("Unexpected recomputation")
  local_mocked_bindings(fit_mfrm_testlet = fail, predict.mfrm_testlet = fail,
    mfrm_testlet_fit = fail, diagnose_mfrm = fail)
  dir <- tempfile()
  out <- export_mfrm_results(res, dir, prefix = "normal", preset = "starter", acknowledge_sensitive = TRUE)
  expect_equal(nrow(out$plot_errors), 0L)
  csv <- out$written_files$Path[out$written_files$Component == "table_calibration"]
  expect_equal(read.csv(csv)$Lower, res$tables$calibration$Lower)
  withr::local_dir(dir)
  env <- new.env()
  capture.output(sys.source("normal_replay.R", envir = env))
  expect_identical(env$res, res)
  expect_identical(env$report$tables$calibration, res$tables$calibration)
  expect_match(env$report$markdown, "Explicit 80% observed-information normal approximation", fixed = TRUE)
  expect_identical(plot(env$res, type = "calibration", draw = FALSE)$data$table, res$tables$calibration)
  # Earlier bundles keep their recorded tables; rebuilding from the source fit
  # applies current defaults without mutating the original object.
  legacy <- res; legacy$calibration_intervals <- NULL
  expect_identical(mfrm_report(legacy)$tables$calibration, legacy$tables$calibration)
  expect_true(all(is.na(mfrm_results(legacy$fit)$tables$calibration$Lower)))
})

test_that("testlet reports preserve the scoring roster and Person-local sharing", {
  fit <- extended_report_fixture()
  fit$input$columns$person <- "Candidate ID"
  roster <- data.frame(Rater = c("A", "A", "A", "B", "A"),
    Score = c(0, 1, 1, 0, NA_real_), check.names = FALSE)
  roster[["Candidate ID"]] <- c("p:1", "p:1", "p2", "p2", "Missing")
  scored <- predict(fit, roster, missing = "omit", persons = c("p:1", "Missing"))
  expect_identical(scored$table$Status, c("available_conditional", "prior_only"))
  expect_match(scored$settings$effect_sharing, "Person/testlet pair", fixed = TRUE)
  expect_match(scored$settings$local_effect_integration, "no stored local-effect mode")
  # Another Person's scores cannot update this Person-local testlet effect.
  isolated <- predict(fit, roster[roster[["Candidate ID"]] == "p:1", ], persons = "p:1")
  expect_equal(scored$table[1, ], isolated$table, tolerance = 1e-10)
  fail <- function(...) stop("Unexpected scoring or fitting")
  local_mocked_bindings(predict.mfrm_testlet = fail, fit_mfrm_testlet = fail)
  result <- mfrm_results(fit, predictions = scored, compute = "never")
  expect_identical(result$tables$scoring_roster, scored$scoring_data)
  expect_identical(result$tables$scoring_roster[["Candidate ID"]], roster[["Candidate ID"]])
  expect_true(is.na(result$tables$scoring_roster$Score[5]))
  expect_identical(result$tables$scoring_blocks, scored$blocks)
  expect_identical(mfrm_report(result)$tables$scoring_roster, scored$scoring_data)
  expect_identical(result$tables$scoring_settings,
    mfrm_extended_settings_table(scored$settings))
  directory <- tempfile(); withr::defer(unlink(directory, recursive = TRUE))
  archive <- export_mfrm_results(result, directory, include = c("tables", "replay"),
    acknowledge_sensitive = TRUE)
  file <- archive$written_files$Path[archive$written_files$Component == "table_scoring_roster"]
  expect_length(file, 1L)
  csv <- read.csv(file, check.names = FALSE)
  expect_identical(names(csv), names(scored$scoring_data))
  expect_equal(csv, scored$scoring_data, ignore_attr = TRUE)
  file <- tempfile(); withr::defer(unlink(file))
  saveRDS(result, file)
  expect_identical(readRDS(file)$tables$scoring_roster, scored$scoring_data)
  legacy <- scored; legacy$scoring_data <- NULL
  legacy$settings[c("effect_sharing", "local_effect_integration", "roster")] <- NULL
  older <- mfrm_results(fit, predictions = legacy, compute = "never")
  expect_null(older$tables$scoring_roster)
  expect_identical(older$tables$person_scores, legacy$table)
})

test_that("random-rater probability settings preserve the integration target", {
  fit <- extended_report_fixture(FALSE)
  ability <- c(-.8, .6)
  known <- predict(fit, data.frame(Rater = c("A", "A")), ability)
  replacement <- predict(fit, data.frame(Rater = c("new", "new")), ability, rater = "new")
  expect_match(known$settings$rater_integration, "Conditional Laplace")
  expect_match(replacement$settings$rater_integration, "replacement-rater population")
  expect_identical(known$settings$effect_sharing, replacement$settings$effect_sharing)
  # Independent integration checks that the labelled distributions are used.
  reference <- function(mu, sd) vapply(ability, function(theta)
    stats::integrate(function(z) stats::plogis(theta - mu - sd * z) * stats::dnorm(z),
      -Inf, Inf, rel.tol = 1e-10)$value, numeric(1))
  expect_equal(known$probabilities[, 2], reference(-.3, .2), tolerance = 1e-9)
  expect_equal(replacement$probabilities[, 2], reference(0, .7), tolerance = 1e-9)
  fail <- function(...) stop("Unexpected prediction or fitting")
  local_mocked_bindings(predict.mfrm_random_rater = fail, fit_mfrm_random_rater = fail)
  for (prediction in list(known, replacement)) {
    res <- mfrm_results(fit, predictions = prediction, compute = "never")
    expect_identical(res$tables$prediction_settings,
      mfrm_extended_settings_table(prediction$settings))
    expect_identical(mfrm_report(res)$tables$prediction_settings, res$tables$prediction_settings)
    file <- tempfile(); withr::defer(unlink(file))
    saveRDS(res, file)
    expect_identical(readRDS(file)$tables$prediction_settings, res$tables$prediction_settings)
  }
})
