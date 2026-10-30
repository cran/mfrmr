test_that("extended normal bounds retain specific causes without changing saved points", {
  fit <- structure(list(calibration_table = data.frame(
    Parameter = c("Fixed facet", "Step", "Population SD"), Facet = c("Task", "Score", "Person"),
    Level = c("A", "1", "Population"), Estimate = c(.4, -.2, 1), SE = c(.3, .1, .2)),
    raters = data.frame(Rater = c("A", "B"), Estimate = c(.4, -.2), PredictionSE = c(.3, .1)),
    checks = list(NumericalReady = TRUE, InformationPositive = TRUE)), class = "mfrm_random_rater")
  before <- serialize(fit, NULL)
  blocked <- function(...) stop("Unexpected recomputation")
  local_mocked_bindings(mfrm_random_rater_fit = blocked, mfrm_random_rater_objective = blocked)
  default <- mfrm_extended_calibration_table(fit)
  expect_identical(default$IntervalStatus, c("not_requested", "not_requested", "unsupported_target"))
  expect_false(any(default$CIEligible))
  expect_true(all(is.na(default$Lower)))
  expect_match(default$InferenceReview[3], "population SDs")
  for (parm in c("calibration", "raters")) {
    good <- confint(fit, parm = parm, level = .8)
    expect_equal(unname(good[, "Lower"]), c(.4, -.2) - qnorm(.9) * c(.3, .1))
    expect_true(all(attr(good, "diagnostics")$CIEligible))
    expect_identical(attr(good, "checks"), fit$checks)
    for (key in c("NumericalReady", "InformationPositive", "EstimatedVarianceBoundary", "EstimatedPersonVarianceBoundary")) {
      changed <- fit; boundary <- grepl("Boundary", key)
      changed$checks[[key]] <- boundary
      ci <- confint(changed, parm = parm)
      tab <- attr(ci, "diagnostics")
      expect_true(all(is.na(ci)))
      expect_equal(tab$Estimate, c(.4, -.2))
      expect_false(any(tab$CIEligible))
      expect_true(all(tab$IntervalStatus == if (boundary) "estimated_variance_boundary" else "numerical_review_failed"))
      expect_match(tab$InferenceReview[1], switch(key,
        NumericalReady = "Numerical readiness", InformationPositive = "Positive information",
        EstimatedVarianceBoundary = "random-effect variance is at zero",
        EstimatedPersonVarianceBoundary = "ability variance is at zero"))
      expect_match(paste(capture.output(print(ci)), collapse = "\n"), tab$IntervalStatus[1])
    }
  }
  combined <- fit
  combined$checks$NumericalReady <- FALSE
  combined$checks$EstimatedVarianceBoundary <- TRUE
  tab <- attr(confint(combined, parm = "calibration"), "diagnostics")
  expect_match(tab$InferenceReview[1], "Numerical readiness.*variance is at zero.*More iterations")
  # A specified zero random-effect variance is not an estimated boundary.
  fit$settings <- list(fixed_rater_sd = 0)
  expect_true(all(attr(confint(fit, parm = "calibration"), "diagnostics")$CIEligible))
  fit$settings <- NULL
  expect_identical(serialize(fit, NULL), before)
  path <- tempfile(); saveRDS(combined, path)
  expect_identical(confint(readRDS(path), parm = "calibration"), confint(combined, parm = "calibration"))
  csv <- tempfile(fileext = ".csv"); write.csv(tab, csv, row.names = FALSE)
  expect_identical(read.csv(csv)$InferenceReview, tab$InferenceReview)
  unlink(c(path, csv))
})

test_that("row-specific invalid normal calculations are not reported as available", {
  fit <- structure(list(checks = list(NumericalReady = TRUE, InformationPositive = TRUE),
    raters = data.frame(Rater = LETTERS[1:5], Estimate = c(.4, Inf, .2, 1e308, .1),
      PredictionSE = c(.3, .1, NA, 1e308, -.1))), class = "mfrm_random_rater")
  ci <- confint(fit, parm = "raters")
  tab <- attr(ci, "diagnostics")
  expect_identical(tab$IntervalStatus, c("available", "invalid_estimate_or_se",
    "invalid_estimate_or_se", "nonfinite_bounds", "invalid_estimate_or_se"))
  expect_true(all(is.na(ci[-1, ])))
  expect_identical(tab$Estimate, fit$raters$Estimate)
  expect_identical(mfrm_random_rater_table(fit, "normal"), tab)
})

test_that("CML uncertainty distinguishes a constrained value from an estimated target", {
  result <- mfrmr_cml_wald(c(.4, 0), matrix(c(1, 0), ncol = 1), matrix(.09), .95)
  expect_identical(result$table$Status, c("available", "fixed"))
  expect_identical(result$table$CIEligible, c(TRUE, FALSE))
  expect_match(result$table$InferenceReview[1], "conditional-information")
  expect_match(result$table$InferenceReview[2], "Constraint-fixed")
  expect_equal(result$table$SE, c(.3, 0))
  expect_equal(result$table$Lower, c(.4 - qnorm(.975) * .3, NA_real_))
  expect_equal(result$table$Upper, c(.4 + qnorm(.975) * .3, NA_real_))
})

test_that("replaying real extended fits preserves interval reasons in summaries and plots", {
  saved <- readRDS(system.file("examples", "extended-models.rds", package = "mfrmr", mustWork = TRUE))
  for (model in saved) {
    fit <- model$fit
    fit$checks$EstimatedPersonVarianceBoundary <- TRUE
    before <- serialize(fit, NULL)
    ci <- confint(fit, parm = "calibration")
    tab <- summary(fit, calibration_intervals = "normal")$calibration
    selected <- tab$Parameter %in% c("Fixed facet", "Step")
    expect_identical(attr(ci, "diagnostics"), tab[selected, , drop = FALSE])
    expect_true(all(tab$IntervalStatus[selected] == "estimated_variance_boundary"))
    expect_match(tab$InferenceReview[which(selected)[1]], "ability variance is at zero")
    p <- if (inherits(fit, "mfrm_testlet")) plot(fit, facet = "Criterion",
      intervals = "normal", draw = FALSE) else plot(fit, intervals = "normal", draw = FALSE)
    expect_true(all(is.na(p$data$table$Lower)))
    expect_true(all(p$data$table$IntervalStatus == "estimated_variance_boundary"))
    expect_identical(serialize(fit, NULL), before)
  }
})

test_that("prior-only scoring explains the retained population values in every route", {
  input <- list(assigned_data = data.frame(Person = "unobserved", Score = NA_real_),
    columns = list(person = "Person"), person = character(), mean = .5, sd = 1.2,
    data = data.frame(), input_rows = 1L, omitted_rows = 1L)
  local_mocked_bindings(mfrm_ordinary_response_input = function(...) input,
    mfrm_response_source = function(...) list())
  ordinary <- score_mfrm_persons(list(), cut_scores = c(Pass = 0))
  expect_identical(ordinary$table$Status, "prior_only")
  expect_match(ordinary$table$Reason, "No observed responses.*assumed population")
  expect_true(all(grepl("No observed responses", ordinary$decisions$Reason)))
  expect_equal(ordinary$table$Estimate, .5)
  saved <- readRDS(system.file("examples", "extended-models.rds", package = "mfrmr", mustWork = TRUE))
  for (model in saved) {
    fit <- model$fit; roster <- fit$input$assigned_data %||% fit$input$data
    extra <- roster[1, ]; extra$Person <- "unobserved"; extra$Score <- NA_real_
    roster <- rbind(roster, extra)
    scores <- if (inherits(fit, "mfrm_testlet")) predict(fit, roster,
      persons = "unobserved", missing = "omit") else score_mfrm_random_rater(fit, roster,
      persons = "unobserved", missing = "omit")
    expect_identical(scores$table$Status, "prior_only")
    expect_match(scores$table$Reason, "No observed responses.*assumed population")
    expect_identical(summary(scores)$scores$Reason, scores$table$Reason)
    expect_equal(scores$table$ConditionalSD, sqrt(scores$settings$person_variance))
  }
})
