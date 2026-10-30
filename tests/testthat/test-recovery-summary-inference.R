recovery_summary_fixture <- function() {
  x <- data.frame(rep = rep(1:3, each = 2), ParameterType = "facet", Facet = "Item",
    ComparisonScale = "logit", Truth = rep(c(-1, 1), 3),
    ErrorAligned = c(-1, 1, -3, 3, -5, 5), SE = 1,
    RecoveryComparable = TRUE, RecoveryBasis = "test", CIEligible = TRUE,
    CoverageBasis = "centered_facet_model_interval", CoverageTruth = 0,
    Lower95 = -1, Upper95 = 1, CIReason = "available")
  x$EstimateAligned <- x$Truth + x$ErrorAligned
  x$ErrorRaw <- x$ErrorAligned
  x
}

test_that("MCSE treats each trial as the independent unit", {
  x <- recovery_summary_fixture()
  runs <- data.frame(rep = 1:3, RunOK = TRUE)
  targets <- unique(x[c("rep", "ParameterType", "Facet", "ComparisonScale")])
  targets$ExpectedRows <- 2L
  s <- mfrmr:::recovery_summarize_rows(x, runs, targets)
  expect_equal(s$Bias, 0)
  expect_equal(s$McseBias, 0)
  expect_equal(s$RMSE, sqrt(mean(c(1, 9, 25))))
  expect_equal(s$McseRMSE, sd(c(1, 9, 25)) / sqrt(3) / (2 * s$RMSE))
  doubled <- rbind(x, x)
  targets$ExpectedRows <- 4L
  s2 <- mfrmr:::recovery_summarize_rows(doubled, runs, targets)
  expect_equal(s2$McseRMSE, s$McseRMSE)
  expect_equal(s2$McseBias, s$McseBias)
  expect_true(is.na(mfrmr:::recovery_summarize_rows(x[x$rep == 1, ])$McseRMSE))
  expect_true(is.na(mfrmr:::recovery_summarize_rows(x[x$rep == 1, ])$McseBias))
  big <- mfrmr:::recovery_trial_metric(c(1e200, -1e200), 1:2, 1:2, TRUE)
  expect_equal(big[["value"]], 1e200)
  expect_equal(big[["mcse"]], 0)
})

test_that("unequal available counts use the joint numerator-denominator delta method", {
  x <- c(1, 2, 4, 9, 10, NA)
  rep <- c(1, 1, 2, 3, 3, 3)
  totals <- matrix(c(3, 4, 19, 0, 2, 1, 2, 0), ncol = 2)
  mu <- sum(totals[, 1]) / sum(totals[, 2])
  gradient <- c(1, -mu) / mean(totals[, 2])
  expected <- sqrt(drop(t(gradient) %*% cov(totals) %*% gradient) / 4)
  value <- mfrmr:::recovery_trial_metric(x, rep, 1:4)
  expect_equal(value[["value"]], mu)
  expect_equal(value[["mcse"]], expected)
  expect_equal(value[["reps"]], 3)
  expect_error(mfrmr:::recovery_trial_metric(x, rep, 1:2), "unrecorded")
})

test_that("failed trials and unavailable intervals retain separate denominators", {
  x <- recovery_summary_fixture()[1:4, ]
  x$Lower95[1] <- .1
  x$CIEligible[4] <- FALSE
  x$SE[4] <- NA
  runs <- data.frame(rep = 1:4, RunOK = c(TRUE, TRUE, FALSE, FALSE))
  targets <- data.frame(rep = 1:4, ParameterType = "facet", Facet = "Item",
    ComparisonScale = "logit", ExpectedRows = 2L)
  s <- mfrmr:::recovery_summarize_rows(x, runs, targets)
  expect_equal(s$AttemptedReps, 4)
  expect_equal(s$FailedReps, 2)
  expect_equal(s$ExpectedRows, 8)
  expect_equal(s$EstimateAvailableRate, .5)
  expect_equal(s$SEAvailableRate, 3/8)
  expect_equal(s$CoverageAvailabilityRate, 3/8)
  expect_equal(s$CoverageRows, 3)
  expect_equal(s$Coverage95, 2/3)
  expect_equal(s$UnavailableCIRows, 5)
  expect_equal(s$CoverageReps, 2)
  totals <- cbind(c(1, 1, 0, 0), c(2, 1, 0, 0))
  gradient <- c(1, -2/3) / mean(totals[, 2])
  expect_equal(s$McseCoverage95, sqrt(drop(t(gradient) %*% cov(totals) %*% gradient) / 4))
  unknown <- mfrmr:::recovery_summarize_rows(x, runs, targets[-4, ])
  expect_true(is.na(unknown$ExpectedRows))
  expect_true(is.na(unknown$CoverageAvailabilityRate))
  expect_equal(unknown$KnownTargetRows, 6)
  failed <- mfrmr:::recovery_summarize_rows(NULL, transform(runs, RunOK = FALSE), targets)
  expect_equal(failed$ExpectedRows, 8)
  expect_equal(failed$FailedReps, 4)
  expect_equal(failed$CoverageAvailabilityRate, 0)
  expect_equal(failed$UnavailableCIRows, 8)
  expect_true(is.na(failed$Coverage95))
  expect_true(is.na(failed$RMSE))
})

test_that("centered intervals use the full covariance and keep source uncertainty", {
  v <- matrix(c(4, 1, 1, 9), 2)
  estimates <- c(2, 5)
  truth <- list(facets = list(Item = c(A = 1, B = 6)))
  fit <- list(config = list(method = "MML", model = "RSM"),
    facets = list(others = data.frame(Facet = "Item", Level = c("A", "B"),
      Estimate = estimates, SE = sqrt(diag(v)))))
  local_mocked_bindings(mfrm_facet_intervals = function(fit, facet, contrasts, ...) {
    expect_equal(rowSums(contrasts), c(A = 0, B = 0))
    centered <- drop(contrasts %*% estimates)
    se <- sqrt(diag(contrasts %*% v %*% t(contrasts)))
    list(table = data.frame(Target = rownames(contrasts), Estimate = centered, SE = se,
      Lower = centered - qnorm(.975) * se, Upper = centered + qnorm(.975) * se,
      CIEligible = TRUE, Status = "available", InferenceReview = "test"))
  }, .package = "mfrmr")
  rows <- mfrmr:::recovery_rows_from_fit(fit, truth, 1)
  expect_equal(rows$RawSE, c(2, 3))
  expect_equal(rows$SE, rep(sqrt(11/4), 2))
  expect_equal(rows$CoverageEstimate, c(-1.5, 1.5), ignore_attr = TRUE)
  expect_equal(rows$CoverageTruth, c(-2.5, 2.5))
  expect_true(all(rows$CIEligible))
  expect_true(all(rows$Covered95))
  incomplete_truth <- truth
  incomplete_truth$facets$Item <- c(truth$facets$Item, C = 8)
  incomplete <- mfrmr:::recovery_rows_from_fit(fit, incomplete_truth, 1)
  expect_true(all(incomplete$CIReason == "incomplete_alignment_group"))
  expect_true(all(is.na(incomplete$Covered95)))
  fit$config$method <- "JML"
  jml <- mfrmr:::recovery_rows_from_fit(fit, truth, 1)
  expect_equal(jml$EstimateAligned, rows$EstimateAligned)
  expect_equal(jml$RawSE, rows$RawSE)
  expect_true(all(is.na(jml$SE)))
  expect_true(all(is.na(jml$Covered95)))
  expect_true(all(grepl("jml_point", jml$CIReason)))
})

test_that("unrepresentable intervals retain valid target SEs", {
  fit <- list(config = list(method = "MML", model = "RSM"),
    facets = list(others = data.frame(Facet = "Item", Level = c("A", "B"),
      Estimate = c(-1, 1), SE = 2)))
  truth <- list(facets = list(Item = c(A = -1, B = 1)))
  local_mocked_bindings(mfrm_facet_intervals = function(...) list(table = data.frame(
    Target = c("A", "B"), Estimate = c(-1, 1), SE = 2,
    Lower = NA_real_, Upper = NA_real_, CIEligible = FALSE,
    Status = "interval_unavailable", InferenceReview = "numeric_resolution")), .package = "mfrmr")
  rows <- mfrmr:::recovery_rows_from_fit(fit, truth, 1)
  expect_equal(rows$SE, c(2, 2))
  expect_true(all(is.na(rows$Covered95)))
  expect_true(all(grepl("numeric_resolution", rows$CIReason)))
})

test_that("truth mapping preserves PCM/GPCM adjacent logits and population coordinates", {
  truth <- list(facets = list(Rater = c(R1 = .2, R2 = .6),
    Item = c(A = -1, B = 1)),
    step_table = data.frame(StepFacet = rep(c("A", "B"), each = 2),
      Step = rep(c("Step_1", "Step_2"), 2), Estimate = c(-1, 3, -4, 0)),
    population = list(coefficients = c("(Intercept)" = .8, X = 1.5)))
  for (model in c("PCM", "GPCM")) {
    mapped <- mfrmr:::recovery_identified_truth(truth, model, "Item", TRUE)
    expect_equal(mapped$facets$Item, c(A = 0, B = -1))
    expect_equal(mapped$population$coefficients, c("(Intercept)" = .9, X = 1.5))
    for (item in c("A", "B")) {
      tau <- truth$step_table$Estimate[truth$step_table$StepFacet == item]
      original <- .8 - truth$facets$Rater[1] - truth$facets$Item[item] - tau
      identified <- mapped$population$coefficients[1] -
        (mapped$facets$Rater[1] - mean(mapped$facets$Rater)) -
        (mapped$facets$Item[item] - mean(mapped$facets$Item)) - (tau - mean(tau))
      expect_equal(unname(identified), unname(original))
    }
  }
  common <- truth
  common$step_table$StepFacet <- "Common"
  common$step_table <- common$step_table[1:2, ]
  mapped <- mfrmr:::recovery_identified_truth(common, "RSM", NULL, TRUE)
  expect_equal(mapped$population$coefficients[[1]], -.6)
  observed <- mfrmr:::recovery_identified_truth(truth, "PCM", "Item", TRUE,
    levels = list(Rater = "R1", Item = c("A", "B")))
  expect_equal(observed$population$coefficients[[1]], 1.1)
})

test_that("legacy summaries cannot retain invalid coverage or stale MCSE", {
  rows <- recovery_summary_fixture()
  rows$CIEligible <- NULL
  rows$Covered95 <- TRUE
  object <- structure(list(recovery = rows, rep_overview = data.frame(rep = 1:3,
    RunOK = TRUE, Converged = TRUE, ElapsedSec = 0, RecoveryRows = 2),
    recovery_summary = data.frame(RMSE = 999, McseRMSE = 0, Coverage95 = 1)),
    class = "mfrm_recovery_simulation")
  s <- summary(object)
  expect_equal(s$recovery_summary$McseBias, 0)
  expect_true(is.na(s$recovery_summary$Coverage95))
  expect_true(is.na(s$recovery_summary$ExpectedRows))
  expect_match(s$recovery_summary$CIReasons, "legacy_interval_basis_unverified")
  expect_true(any(grepl("Legacy point comparisons", s$notes)))
  expect_equal(assess_mfrm_recovery(object)$metric_review$RMSE, s$recovery_summary$RMSE)
  expect_true(any(grepl("Legacy point comparisons", assess_mfrm_recovery(object)$notes)))
  expect_equal(plot(object, draw = FALSE)$data$plot_table$Value, s$recovery_summary$RMSE)
  expect_true(any(grepl("Legacy point comparisons", plot(object, draw = FALSE)$data$notes)))
  expect_error(plot(object, type = "coverage", draw = FALSE), "No finite coverage")
})

test_that("identified zero bias is not an adequacy pass and reviews retain denominators", {
  rows <- recovery_summary_fixture()
  rows$AlignWithinGroup <- TRUE
  targets <- unique(rows[c("rep", "ParameterType", "Facet", "ComparisonScale")])
  targets$ExpectedRows <- 2L
  object <- structure(list(recovery = rows, recovery_targets = targets,
    rep_overview = data.frame(rep = 1:3, RunOK = TRUE, Converged = TRUE,
      ElapsedSec = 0, RecoveryRows = 2)), class = "mfrm_recovery_simulation")
  review <- assess_mfrm_recovery(object, max_abs_bias = .1)
  expect_equal(review$metric_review$BiasBasis, "zero_by_group_identification")
  expect_equal(review$metric_review$Bias, 0)
  expect_equal(review$metric_review$BiasStatus, "not_assessed")
  expect_equal(review$uncertainty_review$ExpectedRows, 6)
  expect_equal(review$uncertainty_review$CoverageRows, 6)
  expect_equal(review$uncertainty_review$McseCoverage95, 0)
  expect_match(review$checklist$NextAction[review$checklist$Item == "Bias threshold"],
    "fixed by group identification")
  printed <- paste(capture.output(print(review)), collapse = "\n")
  expect_match(printed, "CoverageRows")
  expect_match(printed, "zero by construction")
  expect_equal(summary(object)$recovery_targets, tibble::as_tibble(targets))
  p <- plot(object, draw = FALSE)
  expect_equal(p$data$availability_table$CoverageAvailabilityRate, 1)
})

test_that("named recovery thresholds keep their parameter-specific meaning", {
  f <- mfrmr:::recovery_assessment_threshold
  thresholds <- c("facet:Item:logit" = .1, slope = .7, default = 2)
  expect_equal(unname(f(thresholds, "facet", "Item", "logit")), .1)
  expect_equal(unname(f(thresholds, "slope", "SlopeFacet", "log_slope")), .7)
  expect_equal(unname(f(thresholds, "person", "Person", "logit")), 2)
})
