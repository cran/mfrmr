scoring_prior_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      x <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))
      x$new <- x$data[x$data$Person %in% c("P1", "P34"), ]
      x$prior <- list(mean = unname(x$fit$population$coefficients[1]),
        sd = sqrt(x$fit$population$sigma2))
      x$review <- mfrmr:::prediction_local_calibration_review(x$fit)
      cached <<- x
    }
    cached
  }
})

test_that("explicit scoring priors match independently checked retained patterns", {
  x <- scoring_prior_fixture()
  expect_true(x$review$eligible)
  # This matrix checks scoring priors, not repeated Hessians of an unchanged
  # calibration. Reuse one real review only while its exact input is preserved.
  local_mocked_bindings(prediction_local_calibration_review = function(fit, covariance = NULL) {
    expect_identical(fit, x$fit)
    x$review
  }, .package = "mfrmr")
  before <- serialize(x$fit, NULL)
  reference <- read.csv(test_path("fixtures", "mfrm-scoring-prior-reference.csv"))
  people <- unique(as.character(x$data$Person))
  for (n in c(1L, 3L, 6L)) {
    rows <- do.call(rbind, lapply(c("P1", "P34"), function(person) {
      i <- match(person, people)
      d <- x$data[x$data$Person == person, ]
      d <- d[order(d$Rater, d$Criterion), ]
      head(d[(seq_len(6L) + i - 2L) %% 6L + 1L, ], n)
    }))
    for (scenario in c("retained", "higher_mean", "wider")) {
      expected <- reference[reference$Ratings == n & reference$Scenario == scenario, ]
      prior <- list(mean = expected$PriorMean[1], sd = expected$PriorSD[1])
      pred <- predict_mfrm_units(x$fit, rows, scoring_prior = prior)
      expected <- expected[match(pred$estimates$Person, expected$Person), ]
      cols <- c("Estimate", "SD", "Lower", "Upper")
      expect_equal(unname(as.matrix(pred$estimates[cols])),
        unname(as.matrix(expected[cols])), tolerance = 1e-10)
      expect_true(all(pred$estimates$ScoreIntegrationReady))
      expect_equal(pred$estimates$RetainedPriorSD, rep(x$prior$sd, 2))
    }
  }
  ordinary <- predict_mfrm_units(x$fit, x$new, n_draws = 3, seed = 42)
  omitted <- predict_mfrm_units(x$fit, x$new, n_draws = 3, seed = 42, scoring_prior = NULL)
  identical_prior <- predict_mfrm_units(x$fit, x$new, n_draws = 3, seed = 42, scoring_prior = x$prior)
  expect_identical(ordinary, omitted)
  expect_identical(ordinary$draws$Value, identical_prior$draws$Value)
  expect_equal(ordinary$estimates$Estimate, identical_prior$estimates$Estimate, tolerance = 1e-14)
  expect_identical(serialize(x$fit, NULL), before)
})

test_that("prior identity accompanies summaries, draws, figure data and exported estimates", {
  x <- scoring_prior_fixture()
  local_mocked_bindings(prediction_local_calibration_review = function(fit, covariance = NULL) {
    expect_identical(fit, x$fit)
    x$review
  }, .package = "mfrmr")
  prior <- x$prior; prior$sd <- 1.5 * prior$sd
  p <- predict_mfrm_units(x$fit, x$new, scoring_prior = prior, n_draws = 3, seed = 42)
  pv <- sample_mfrm_plausible_values(x$fit, x$new, scoring_prior = prior, n_draws = 3, seed = 42)
  expect_identical(p$draws, pv$values)
  expect_identical(p$settings, pv$settings)
  for (table in list(summary(p, digits = 0)$estimates, summary(pv, digits = 0)$estimates,
                     summary(pv, digits = 0)$draw_summary, pv$values)) {
    expect_true(all(table$PriorSource == "user_supplied"))
    expect_true(all(table$PriorSD == prior$sd))
    expect_true(all(table$RetainedPriorSD == x$prior$sd))
  }
  expect_match(paste(capture.output(print(summary(p))), collapse = " "), "User-supplied scoring prior")
  expect_match(paste(pv$notes, collapse = " "), "analyst assumption")
  plot <- ggplot2::ggplot(p$estimates, ggplot2::aes(x = Person, y = Estimate)) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = Lower, ymax = Upper))
  expect_identical(plot$data$PriorSource, p$estimates$PriorSource)
  expect_silent(ggplot2::ggplot_build(plot))
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  export_summary_appendix(p, output_dir = path, include_html = FALSE)
  files <- list.files(path, "estimates.*[.]csv$", recursive = TRUE, full.names = TRUE)
  expect_length(files, 1)
  exported <- read.csv(files[1])
  expect_true(all(exported$PriorSource == "user_supplied"))
  expect_equal(exported$RetainedPriorSD, rep(x$prior$sd, 2), tolerance = 1e-14)
  expect_equal(exported$PriorSD, rep(prior$sd, 2), tolerance = 1e-14)
  for (field in c("prior_comparison", "retained_posterior_basis")) {
    bad <- p; bad$settings[[field]] <- NULL
    expect_error(summary(bad), "prior records")
    expect_error(export_summary_appendix(bad, output_dir = path), "prior records")
  }
  bad <- p; bad$settings$scoring_prior$mean <- prior$mean + 1
  expect_error(summary(bad), "prior records")
  bad <- p; bad$estimates$RetainedPriorMean[1] <- 99
  expect_error(summary(bad), "prior records")
})

test_that("explicit priors cannot override source or integration failures", {
  # No mocked review here: invalid sources and scoring grids must be refused
  # using fresh likelihood, curvature and integration checks.
  x <- scoring_prior_fixture()
  invalid <- list(list(mean = 0), c(mean = 0, sd = 1), list(mean = 0, sd = 1, extra = 2),
    list(mean = Inf, sd = 1), list(mean = NA_real_, sd = 1),
    list(mean = 0, sd = 0), list(mean = 0, sd = -1), list(mean = 0, sd = Inf),
    list(mean = 0, sd = 1e300), list(mean = 0, sd = 1e-300),
    list(mean = 0, sd = "1"), list(mean = c(0, 1), sd = 1),
    list(mean = matrix(0), sd = 1))
  for (prior in invalid) for (policy in c("error", "review")) {
    expect_error(predict_mfrm_units(x$fit, x$new, scoring_prior = prior,
      readiness_policy = policy), "scoring_prior")
  }
  failed <- x$fit; failed$readiness$fit$NumericalState <- "failed"
  expect_error(predict_mfrm_units(failed, x$new, scoring_prior = x$prior), "not ready")
  review <- predict_mfrm_units(failed, x$new, scoring_prior = x$prior, readiness_policy = "review")
  expect_true(all(!review$estimates$SourceScoringReady))
  expect_true(all(review$estimates$EstimateUse == "review_only_nonready_source"))
  expect_s3_class(summary(review), "summary.mfrm_unit_prediction")
  failed$opt$par[1] <- Inf
  expect_error(predict_mfrm_units(failed, x$new, scoring_prior = x$prior,
    readiness_policy = "review"), "invalid")
  expect_error(predict_mfrm_units(x$fit, x$new, scoring_prior = x$prior,
    scoring_quad_points = 2), "integration did not pass")
  coarse <- predict_mfrm_units(x$fit, x$new, scoring_prior = x$prior,
    scoring_quad_points = 2, readiness_policy = "review")
  expect_true(any(!coarse$estimates$ScoreIntegrationReady))
  expect_s3_class(summary(coarse), "summary.mfrm_unit_prediction")
})
