test_that("person likelihoods preserve small probabilities and response coverage", {
  obs <- data.frame(Person = c("p", "p"), PrObserved = c(1e-50, 0.2),
                    ItemEntropy = c(-1, -1), ItemVarLogP = c(1, 1))
  result <- compute_person_fit_indices(list(obs = obs))
  expect_equal(result$LogLik, log(1e-50) + log(0.2))
  expect_equal(result$lz, (log(1e-50) + log(0.2) + 2) / sqrt(2))
  for (invalid in c(NA_real_, 0, -0.2, 1.2)) {
    obs$PrObserved[2] <- invalid
    unavailable <- compute_person_fit_indices(list(obs = obs))
    expect_equal(unavailable$TotalN, 2)
    expect_equal(unavailable$N, 1)
    expect_equal(unavailable$UnavailableN, 1)
    expect_true(is.na(unavailable$LogLik))
    expect_true(is.na(unavailable$lz_flag_5pct))
    expect_true(is.na(unavailable$lz_star_flag_5pct))
    expect_true(is.na(unavailable$ReportFlag))
    expect_true(is.na(summary(unavailable)$overview$FlagRate))
  }
})

test_that("native person-fit inputs retain identity and current probability moments", {
  fit <- make_toy_fit()
  diag <- make_toy_diagnostics(fit)
  other <- diag
  other$obs$Score[1] <- other$obs$Score[1] + 1
  expect_error(compute_person_fit_indices(other, fit), "fit/diagnostics mismatch")
  expect_error(plot_person_fit(fit, other, draw = FALSE), "fit/diagnostics mismatch")
  current <- compute_person_fit_indices(diag, fit)
  old <- diag
  attr(old$obs, "person_fit_moments_version") <- NULL
  expect_equal(compute_person_fit_indices(old, fit), current)
  expect_error(compute_person_fit_indices(old), "Recreate these diagnostics")

  probability_bundle <- compute_response_probability_bundle
  local_mocked_bindings(compute_response_probability_bundle = function(...) {
    bundle <- probability_bundle(...)
    bundle$probs[1, ] <- c(1, 1e-20, rep(0, ncol(bundle$probs) - 2L))
    bundle
  })
  obs <- compute_obs_table(fit)
  expected_mean <- 1e-20 * log(1e-20)
  expected_variance <- 1e-20 * log(1e-20)^2 - expected_mean^2
  expect_equal(obs$ItemEntropy[1] / expected_mean, 1, tolerance = 1e-12)
  expect_equal(obs$ItemVarLogP[1] / expected_variance, 1, tolerance = 1e-12)
})

test_that("ability correction does not silently discard derivative terms or ignore weights", {
  fit <- make_toy_fit()
  obs <- make_toy_diagnostics(fit)$obs
  person <- as.character(obs$Person[1])
  obs$ItemLogPScoreCov[1] <- NA_real_
  correction <- compute_snijders_lz_star(obs, person, fit)
  expect_true(is.na(correction$lz_star))
  expect_equal(correction$lz_star_status, "incomplete_observation_terms")
  obs$Weight <- 2
  expect_equal(compute_snijders_lz_star(obs, person, fit)$lz_star_status, "nonunit_weights")
  obs <- make_toy_diagnostics(fit)$obs
  obs$PersonMeasure[1] <- Inf
  expect_equal(compute_snijders_lz_star(obs, person, fit)$lz_star_status, "nonfinite_person_estimate")
})

test_that("ability correction requires established numerical convergence", {
  fit <- make_toy_fit(maxit = 1)
  diag <- make_toy_diagnostics(fit)
  expect_false(identical(mfrmr_get_readiness_record(fit)$fit$NumericalState[1], "ready"))
  result <- compute_person_fit_indices(diag, fit)
  expect_true(all(is.na(result$lz_star)))
  expect_true(all(result$lz_star_status == "fit_not_converged"))
  expect_true(any(is.finite(result$lz)))
  expect_true(all(grepl("Numerical convergence", result$ReportCaveat)))
})

test_that("person-fit plots keep unavailable persons and separate their indices", {
  fit <- make_toy_fit()
  diag <- make_toy_diagnostics(fit)
  persons <- diag$measures$Facet == "Person"
  diag$measures$Infit[persons] <- NA_real_
  diag$measures$Outfit[persons] <- NA_real_
  p <- plot_person_fit(fit, diag, fit_index = "loglik", draw = FALSE)
  expect_equal(nrow(p$data$data), sum(persons))
  expect_true(any(is.finite(p$data$data$ReportValue)))
  expect_equal(p$data$flag_summary$AvailableRows[1], 0)
  expect_true(is.na(p$data$flag_summary$FlaggedRows[1]))
  expect_equal(p$data$flag_summary$UnavailableRows[1], sum(persons))
  expect_false(any(grepl("1% threshold|5% threshold|Review 5%", p$data$legend$label)))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_person_fit(fit, diag, fit_index = "meansquare"))
  expect_no_error(plot_person_fit(fit, diag, fit_index = "loglik"))
})

test_that("person-fit console output uses plain explanations and guards saved results", {
  fit <- make_toy_fit()
  result <- compute_person_fit_indices(make_toy_diagnostics(fit), fit)
  printed <- paste(capture.output(print(result)), capture.output(print(summary(result, include_person = TRUE))), collapse = "\n")
  expect_match(printed, "Person-Fit Summary")
  expect_false(grepl("computed_jml_conditional_calibration|lz_star_status|review_1pct|ReportIndex", printed))
  old <- result
  attr(old, "person_fit_calculation_version") <- NULL
  expect_error(summary(old), "Recreate this person-fit result")
  expect_error(print(old), "Recreate this person-fit result")
  saved <- summary(result)
  saved$calculation_version <- NULL
  expect_error(print(saved), "Recreate this summary")
})

test_that("unexpected-response counts preserve three-valued screening and QC coverage", {
  obs <- data.frame(Person = c("a", "b", "c", "d"), score_k = 0,
                    StdResidual = c(NA, NA, 3, 0), Observed = 0,
                    Residual = c(NA, NA, 1, 0))
  probs <- rbind(c(.1, .9), c(.9, .1), c(NA, NA), c(.9, .1))
  either <- calc_unexpected_response_table(obs, probs, character(), 0, rule = "either")
  es <- summarize_unexpected_response_table(either, nrow(obs))
  expect_equal(es$EvaluatedObservations, 3)
  expect_equal(es$UnavailableObservations, 1)
  expect_equal(es$UnexpectedN, 2)
  expect_true(is.na(es$UnexpectedPercent))
  both <- calc_unexpected_response_table(obs, probs, character(), 0, rule = "both")
  bs <- summarize_unexpected_response_table(both, nrow(obs), rule = "both")
  expect_equal(bs$EvaluatedObservations, 2)
  expect_equal(bs$UnavailableObservations, 2)
  expect_true(is.na(bs$UnexpectedPercent))

  fit <- make_toy_fit()
  diag <- make_toy_diagnostics(fit)
  diag$obs$StdResidual[] <- NA_real_
  u <- unexpected_response_table(fit, diag, prob_max = 0)
  expect_equal(u$summary$EvaluatedObservations, 0)
  expect_true(is.na(u$summary$UnexpectedN))
  diag$unexpected <- u
  qc <- run_qc_pipeline(fit, diag, include_bias = FALSE)
  row <- qc$verdicts[qc$verdicts$Check == "Unexpected Responses", ]
  expect_equal(row$Verdict, "Warn")
  expect_match(row$Detail, "unavailable")
  stale <- qc
  stale$config$screening_coverage_version <- NULL
  expect_error(print(stale), "Recreate it with run_qc_pipeline", fixed = TRUE)
  expect_error(plot_qc_pipeline(stale, draw = FALSE), "Recreate it with run_qc_pipeline", fixed = TRUE)
})
