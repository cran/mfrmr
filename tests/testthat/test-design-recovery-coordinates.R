test_that("design recovery uses identified step locations through the public runner", {
  dat <- expand.grid(Person = c("P1", "P2", "P3"), Rater = c("R1", "R2"),
    Criterion = c("A", "B"), stringsAsFactors = FALSE)
  dat$Score <- rep(1:3, length.out = nrow(dat))
  attr(dat, "mfrm_truth") <- list(person = c(P1 = -1, P2 = 0, P3 = 1),
    facets = list(Rater = c(R1 = .2, R2 = .6), Criterion = c(A = -1, B = 1)),
    step_table = data.frame(StepFacet = rep(c("A", "B"), each = 2),
      Step = rep(c("Step_1", "Step_2"), 2), Estimate = c(-1, 3, -4, 0)))
  diagnostics <- list(
    reliability = data.frame(Facet = c("Rater", "Criterion"), Levels = 2,
      Separation = 3, Strata = 4, Reliability = .9, MeanInfit = 1, MeanOutfit = 1),
    fit = data.frame(Facet = c("Rater", "Criterion"), InfitZSTD = 0, OutfitZSTD = 0),
    measures = data.frame(Facet = rep(c("Rater", "Criterion"), each = 2),
      Level = c("R1", "R2", "A", "B"), Estimate = c(-.2, .2, .5, -.5)))
  local_mocked_bindings(simulate_mfrm_data = function(...) dat,
    fit_mfrm = function(..., model) list(config = list(model = model,
      facet_names = c("Rater", "Criterion"), step_facet = "Criterion")),
    diagnose_mfrm = function(...) diagnostics, mfrm_inference_ready = function(...) TRUE,
    .package = "mfrmr")
  for (model in c("PCM", "GPCM")) {
    x <- suppressWarnings(evaluate_mfrm_design(n_person = 3, n_rater = 2, n_criterion = 2,
      raters_per_person = 2, reps = 2, score_levels = 3, model = model,
      step_facet = "Criterion", fit_method = "MML", seed = 3, progress = FALSE))
    expect_equal(x$results$SeverityRMSE, rep(0, 4), tolerance = 1e-14)
    expect_equal(x$results$SeverityBiasRaw[x$results$Facet == "Criterion"], c(.5, .5))
    expect_true(all(grepl("identified_step_locations", x$results$RecoveryBasis)))
    s <- summary(x)
    expect_equal(s$design_summary$RecoveryAvailableReps, c(2L, 2L))
    expect_true(any(grepl("zero by construction", s$notes)))
    expect_equal(nrow(recommend_mfrm_design(s)$recommended), 1)
    old_summary <- s
    old_summary$recovery_coordinates <- NULL
    expect_error(recommend_mfrm_design(old_summary), "recovery-coordinate review")
    legacy <- x
    legacy$results$RecoveryBasis <- "generator_fit_contract_aligned"
    legacy$results$SeverityRMSE <- .01
    reviewed <- summary(legacy)
    expect_true(all(is.na(reviewed$design_summary$MeanSeverityRMSE)))
    expect_false(any(is.nan(reviewed$design_summary$MeanSeverityRMSE)))
    expect_equal(reviewed$design_summary$RecoveryAvailableReps, c(0L, 0L))
    expect_true(any(grepl("Legacy PCM/GPCM", reviewed$notes)))
    expect_equal(nrow(recommend_mfrm_design(legacy)$recommended), 0)
    expect_true(is.na(recommend_mfrm_design(legacy)$design_table$MaxSeverityRMSE))
    expect_equal(legacy$results$SeverityRMSE, rep(.01, 4))
  }
})

test_that("design RMSE is finite for representable large errors and excludes unknown truth", {
  metric <- mfrmr:::design_eval_recovery_metrics(c("A", "B", "C"),
    c(1e200, -1e200, 10), c(A = 0, B = 0, C = NA))
  expect_equal(metric$raw_rmse, 1e200)
  expect_equal(metric$aligned_rmse, 1e200)
  expect_equal(metric$raw_bias, 0)
  expect_equal(metric$aligned_bias, 0)
})
