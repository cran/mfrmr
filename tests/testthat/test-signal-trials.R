signal_trial_bindings <- function(.env = parent.frame()) {
  local_mocked_bindings(
    simulate_mfrm_data = function(seed = NULL, ...) with_preserved_rng_seed(seed, {
      runif(1)
      data.frame(Person = rep(c("P1", "P2"), each = 4),
        Rater = rep(c("R01", "R02"), 4), Criterion = rep(c("C01", "C02"), each = 2),
        Group = rep(c("A", "B"), each = 4), Score = rep(1:4, 2))
    }),
    fit_mfrm = function(...) list(draw = rnorm(1)),
    diagnose_mfrm = function(...) list(draw = rnorm(1)),
    mfrm_inference_ready = function(...) TRUE,
    analyze_dff = function(fit, diagnostics, ...) list(dif_table = data.frame(
      Level = c("C01", "C02"), Group1 = "A", Group2 = "B",
      Contrast = fit$draw + diagnostics$draw + rnorm(1),
      ContrastComparable = TRUE, ScaleLinkStatus = "linked",
      p_adjusted = .01, p_value = .001, ETS = NA_character_,
      Classification = "Screen positive", ClassificationSystem = "screening")),
    estimate_bias = function(...) list(table = data.frame(
      FacetA_Level = c("R01", "R02"), FacetB_Level = "C02",
      `Bias Size` = rnorm(1), `Prob.` = .01, t = 3, check.names = FALSE)),
    .env = .env)
}

signal_trial_run <- function(...) evaluate_mfrm_signal_detection(
  n_person = 8, n_rater = 2, n_criterion = 2, raters_per_person = 2,
  dif_method = "refit", dif_abs_cut = 0, seed = 124, ...)

test_that("subgroup refits preserve declared independent repeated ratings without inferring them", {
  fit <- list(config = list(model = "RSM", method = "MML", rating_min = 1,
    rating_max = 4, repeated_ratings = "review"))
  expect_identical(resolve_dff_refit_controls(fit)$repeated_ratings, "review")
  fit$config$replay_inputs <- list(repeated_ratings = "independent")
  expect_identical(resolve_dff_refit_controls(fit)$repeated_ratings, "independent")
  fit$config$replay_inputs <- NULL
  fit$config$repeated_ratings <- "independent"
  expect_identical(resolve_dff_refit_controls(fit)$repeated_ratings, "independent")
  fit$config$repeated_ratings <- NULL
  expect_identical(resolve_dff_refit_controls(fit)$repeated_ratings, "review")
})

test_that("signal contrast direction follows requested groups, not table order", {
  signal_trial_bindings()
  x <- signal_trial_run(reps = 2)
  reversed <- signal_trial_run(reps = 2, reference_group = "B", focal_group = "A")
  expect_equal(reversed$results$DIFContrast, -x$results$DIFContrast)
  expect_identical(reversed$results$DIFReferenceGroup, rep("B", 2))
  expect_identical(reversed$results$DIFFocalGroup, rep("A", 2))
  expect_identical(reversed$results$DIFContrastDirection, rep("reference_minus_focal", 2))
  expect_identical(reversed$results$DIFDetected, x$results$DIFDetected)
  expect_equal(summary(reversed)$detection_summary$MeanTargetContrast,
    -summary(x)$detection_summary$MeanTargetContrast)
  old <- x; old$results$DIFContrastDirection <- NULL
  expect_true(is.na(summary(old)$detection_summary$MeanTargetContrast))
  expect_equal(summary(old)$detection_summary$MeanTargetContrastAbs,
    summary(x)$detection_summary$MeanTargetContrastAbs)
  expect_equal(summary(old)$detection_summary$DIFPower, summary(x)$detection_summary$DIFPower)
  expect_true(all(is.na(plot(old, metric = "estimate", draw = FALSE)$data$y)))
  expect_true(any(grepl("direction", plot(old, metric = "estimate", draw = FALSE)$notes)))
})

test_that("signal seeds cover all trial stages and earlier failures cannot shift later trials", {
  signal_trial_bindings()
  set.seed(18); before <- .Random.seed
  x <- signal_trial_run(reps = 3)
  expect_identical(.Random.seed, before)
  set.seed(99); before <- .Random.seed
  y <- signal_trial_run(reps = 3)
  expect_identical(.Random.seed, before)
  expect_identical(x$results$DIFContrast, y$results$DIFContrast)
  seeds <- with_preserved_rng_seed(124, sample.int(.Machine$integer.max, 3))
  expect_identical(x$rep_overview$Seed, seeds)
  expect_identical(x$results$Seed, seeds)
  expected <- vapply(seeds, function(seed) with_preserved_rng_seed(seed, {
    runif(1); sum(rnorm(3))
  }), numeric(1))
  expect_equal(x$results$DIFContrast, expected)
  generate <- simulate_mfrm_data
  calls <- 0L
  local_mocked_bindings(simulate_mfrm_data = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) { runif(7); stop("Injected generation failure") }
    generate(...)
  })
  z <- signal_trial_run(reps = 3)
  expect_identical(.Random.seed, before)
  expect_equal(nrow(z$results), 3)
  expect_identical(z$results$DIFContrast[-1], x$results$DIFContrast[-1])
  expect_identical(z$results$BiasSize[-1], x$results$BiasSize[-1])
  expect_identical(z$rep_overview$RunOK, c(FALSE, TRUE, TRUE))
  expect_identical(z$rep_overview$ErrorStage, c("generation", NA_character_, NA_character_))
  expect_match(z$rep_overview$Error[1], "Injected generation failure")
  expect_match(z$results$DIFScreenReason[1], "Injected generation failure")
  expect_match(z$results$BiasScreenReason[1], "Injected generation failure")
  expect_true(is.na(z$rep_overview$Observations[1]))
  expect_equal(z$results$DIFNonTargetAvailable, c(0, 1, 1))
  s <- summary(z)$detection_summary
  expect_equal(s$DIFPlanned, 3)
  expect_equal(s$DIFAvailable, 2)
  expect_equal(s$DIFAllTrialsLower, 2/3)
  expect_equal(s$DIFAllTrialsUpper, 1)
})

test_that("all failed generation retains unavailable summaries and failure stages", {
  signal_trial_bindings()
  local_mocked_bindings(simulate_mfrm_data = function(...) stop("No cohort"))
  x <- signal_trial_run(reps = 2)
  s <- summary(x)$detection_summary
  expect_equal(s$Reps, 2)
  expect_equal(s$DIFPlanned, 2)
  expect_equal(s$DIFAvailable, 0)
  expect_equal(s$DIFNonTargetPlanned, 2)
  expect_equal(s$DIFNonTargetAvailable, 0)
  expect_equal(s$DIFNonTargetAvailableReps, 0)
  expect_true(is.na(s$MeanTargetContrast) && !is.nan(s$MeanTargetContrast))
  expect_true(is.na(s$MeanTargetBias) && !is.nan(s$MeanTargetBias))
  expect_true(all(is.na(plot(x, metric = "estimate", draw = FALSE)$data$y)))
  expect_equal(summary(x)$overview$SuccessfulRuns, 0)
})

test_that("non-target averages use available trials and report unequal cell denominators", {
  signal_trial_bindings()
  x <- signal_trial_run(reps = 3)
  x$results$DIFContrast <- c(1, 2, NA)
  x$results$DIFDetected <- c(TRUE, FALSE, NA)
  x$results$DIFClassDetected <- c(TRUE, NA, NA)
  x$results$DIFNonTargetPlanned <- 2L
  x$results$DIFNonTargetAvailable <- c(1L, 2L, 0L)
  x$results$DIFFalsePositiveRate <- c(1, 0, NA)
  x$results$BiasNonTargetAvailable <- c(1L, 2L, 0L)
  x$results$BiasScreenFalsePositiveRate <- c(1, 0, NA)
  s <- summary(x)$detection_summary
  expect_equal(s$DIFFalsePositiveRate, .5) # Not the pooled-cell rate, 1/3.
  expect_equal(s$McseDIFFalsePositiveRate, .5)
  expect_equal(s$McseBiasScreenFalsePositiveRate, .5)
  expect_equal(s$BiasScreenFalsePositiveRate, .5)
  expect_equal(s$DIFNonTargetPlanned, 6)
  expect_equal(s$DIFNonTargetAvailable, 3)
  expect_equal(s$DIFNonTargetAvailableReps, 2)
  expect_equal(s$DIFClassAvailable, 1)
  expect_equal(s$TargetContrastAvailable, 2)
  expect_equal(s$DIFPower, .5)
  expect_equal(s$DIFAllTrialsLower, 1/3)
  expect_equal(s$DIFAllTrialsUpper, 2/3)
  p <- plot(x, metric = "estimate", draw = FALSE)
  expect_equal(p$summary_table, s)
  expect_match(p$display_metric, "reference - focal", fixed = TRUE)
  save <- tempfile(fileext = ".rds"); on.exit(unlink(save), add = TRUE)
  saveRDS(x, save)
  expect_identical(summary(readRDS(save)), summary(x))
})

test_that("signal stage failures preserve other usable screens and absent caller RNG", {
  signal_trial_bindings()
  before <- if (exists(".Random.seed", .GlobalEnv)) .Random.seed else NULL
  on.exit(if (is.null(before)) {
    if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", before, envir = .GlobalEnv), add = TRUE)
  if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  local_mocked_bindings(analyze_dff = function(...) { runif(5); stop("DIF unavailable") })
  x <- signal_trial_run(reps = 1)
  expect_false(exists(".Random.seed", .GlobalEnv))
  expect_identical(x$rep_overview$ErrorStage, "dif")
  expect_false(x$rep_overview$RunOK)
  expect_true(x$results$BiasDetected)
  expect_true(is.na(x$results$DIFDetected))
  local_mocked_bindings(diagnose_mfrm = function(...) stop("Diagnostics unavailable"))
  y <- signal_trial_run(reps = 1)
  expect_identical(y$rep_overview$ErrorStage, "diagnostics")
  expect_true(y$rep_overview$Converged)
  expect_true(is.na(y$results$BiasDetected))
  local_mocked_bindings(fit_mfrm = function(...) stop("Fit unavailable"))
  z <- signal_trial_run(reps = 1)
  expect_identical(z$rep_overview$ErrorStage, "fit")
  expect_false(z$rep_overview$Converged)
  expect_equal(z$rep_overview$Observations, 8)
})

test_that("signal groups inherit a supplied specification unless explicitly overridden", {
  signal_trial_bindings()
  spec <- build_mfrm_sim_spec(n_person = 8, n_rater = 2, n_criterion = 2,
    group_levels = c("Reference", "Focal"))
  x <- signal_trial_run(reps = 1, sim_spec = spec)
  expect_identical(x$settings$group_levels, c("Reference", "Focal"))
  expect_identical(x$settings$reference_group, "Reference")
  expect_identical(x$settings$focal_group, "Focal")
  y <- signal_trial_run(reps = 1, sim_spec = spec, group_levels = c("A", "B"))
  expect_identical(y$settings$group_levels, c("A", "B"))
})

test_that("missing adjusted screens never fall back silently to raw tail areas", {
  signal_trial_bindings()
  local_mocked_bindings(analyze_dff = function(...) list(dif_table = data.frame(
    Level = c("C01", "C02"), Group1 = "A", Group2 = "B", Contrast = 1,
    p_adjusted = NA_real_, p_value = .001, ETS = NA_character_,
    Classification = "Screen positive", ClassificationSystem = "screening")))
  x <- signal_trial_run(reps = 1)
  expect_true(is.na(x$results$DIFP))
  expect_true(is.na(x$results$DIFDetected))
  expect_equal(x$results$DIFNonTargetAvailable, 0)
  expect_true(is.na(x$results$DIFFalsePositiveRate))
  expect_error(signal_trial_run(dif_effect = NA_real_), "dif_effect")
  expect_error(signal_trial_run(bias_effect = c(0, 1)), "bias_effect")
  expect_error(signal_trial_run(reference_group = c("A", "B")), "reference_group")
})

test_that("unlinked refit contrasts remain recorded but cannot be averaged as common-scale effects", {
  signal_trial_bindings()
  local_mocked_bindings(analyze_dff = function(...) list(dif_table = data.frame(
    Level = c("C01", "C02"), Group1 = "A", Group2 = "B", Contrast = 1,
    ContrastComparable = FALSE, ScaleLinkStatus = "unlinked",
    ConditionalScreenReason = "No shared anchors", p_adjusted = NA_real_,
    p_value = NA_real_, ETS = NA_character_, Classification = "Descriptive",
    ClassificationSystem = "descriptive")))
  x <- signal_trial_run(reps = 1)
  expect_equal(x$results$DIFContrast, 1)
  expect_false(x$results$DIFContrastComparable)
  expect_identical(x$results$DIFScaleLinkStatus, "unlinked")
  expect_identical(x$results$DIFScreenReason, "No shared anchors")
  s <- summary(x)$detection_summary
  expect_true(is.na(s$MeanTargetContrast))
  expect_true(is.na(s$MeanTargetContrastAbs))
  expect_equal(s$TargetContrastAvailable, 0)
  old <- x; old$results$DIFContrastComparable <- NULL
  expect_true(is.na(summary(old)$detection_summary$MeanTargetContrastAbs))
  expect_true(any(grepl("common-scale comparability", summary(old)$notes)))
})
