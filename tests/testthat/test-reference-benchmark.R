test_that("reference_case_benchmark returns package-native benchmark bundle", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = c("synthetic_truth", "synthetic_bias_contract", "study1_itercal_pair"),
    method = "MML",
    quad_points = 5,
    maxit = 30
  ))

  expect_s3_class(bench, "mfrm_reference_benchmark")
  expect_true(all(c(
    "overview", "summary", "table", "fit_runs", "case_summary",
    "design_checks", "recovery_checks", "bias_checks", "pair_checks",
    "linking_checks", "conquest_overlap_checks", "population_policy_checks",
    "source_profile", "settings", "notes"
  ) %in% names(bench)))

  expect_true(is.data.frame(bench$overview))
  expect_true(is.data.frame(bench$fit_runs))
  expect_true(is.data.frame(bench$case_summary))
  expect_true(is.data.frame(bench$design_checks))
  expect_true(is.data.frame(bench$recovery_checks))
  expect_true(is.data.frame(bench$bias_checks))
  expect_true(is.data.frame(bench$pair_checks))
  expect_true(is.data.frame(bench$linking_checks))
  expect_true(is.data.frame(bench$conquest_overlap_checks))
  expect_true(is.data.frame(bench$population_policy_checks))
  expect_true(is.data.frame(bench$source_profile))
  expect_identical(bench$settings$intended_use, "internal_benchmark")
  expect_false(isTRUE(bench$settings$external_validation))
  expect_true(all(c(
    "MMLEngineRequested", "MMLEngineUsed", "EMIterations",
    "PosteriorBasis", "PopulationModelActive", "PopulationDesignColumns"
  ) %in% names(bench$fit_runs)))

  expect_equal(nrow(bench$case_summary), 3)
  expect_true(all(c("synthetic_truth", "synthetic_bias_contract", "study1_itercal_pair") %in% bench$case_summary$Case))
  expect_true(all(bench$fit_runs$PrecisionTier %in% c("model_based", "hybrid", "exploratory")))
})

test_that("reference_case_benchmark includes latent-regression omission contract case", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_latent_regression_omit",
    method = "MML",
    quad_points = 7,
    maxit = 40
  ))

  expect_s3_class(bench, "mfrm_reference_benchmark")
  expect_true("synthetic_latent_regression_omit" %in% bench$case_summary$Case)
  expect_true(is.data.frame(bench$population_policy_checks))
  expect_true(all(bench$population_policy_checks$Status == "Pass"))
  expect_true(all(c(
    "PopulationPolicy",
    "PopulationOmittedPersons",
    "PopulationResponseRowsOmitted",
    "PopulationResponseRowsRetained",
    "OmittedPersonExcludedFromEstimates",
    "OmittedPersonPreservedForReplay"
  ) %in% bench$population_policy_checks$Metric))

  run <- bench$fit_runs[
    bench$fit_runs$Case == "synthetic_latent_regression_omit",
    ,
    drop = FALSE
  ]
  expect_identical(as.character(run$PopulationPolicy[1]), "omit")
  expect_identical(run$PopulationOmittedPersons[1], 1L)
  expect_identical(run$PopulationResponseRowsOmitted[1], 6L)
  expect_identical(run$PopulationIncludedPersons[1], 59L)

  row <- bench$case_summary[
    bench$case_summary$Case == "synthetic_latent_regression_omit",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(row), 1)
  expect_identical(as.character(row$Status[1]), "Pass")
  expect_identical(row$PopulationPolicyChecks[1], nrow(bench$population_policy_checks))
  expect_match(row$KeySignal[1], "Population response rows omitted = 6", fixed = TRUE)
})

test_that("reference_case_benchmark includes ConQuest-overlap package-side check", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_conquest_overlap_dry_run",
    method = "MML",
    model = "RSM",
    quad_points = 7,
    maxit = 40
  ))

  expect_s3_class(bench, "mfrm_reference_benchmark")
  expect_true("synthetic_conquest_overlap_dry_run" %in% bench$case_summary$Case)
  expect_true(is.data.frame(bench$conquest_overlap_checks))
  expect_true(all(bench$conquest_overlap_checks$Status == "Pass"))
  expect_true(all(c(
    "AttentionItems",
    "PopulationMaxAbsDifference",
    "ItemCenteredMaxAbsDifference",
    "CaseMaxAbsDifference"
  ) %in% bench$conquest_overlap_checks$Metric))
  expect_identical(bench$settings$external_validation, FALSE)

  row <- bench$case_summary[
    bench$case_summary$Case == "synthetic_conquest_overlap_dry_run",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(row), 1)
  expect_identical(as.character(row$Status[1]), "Pass")
  expect_identical(row$ConQuestOverlapChecks[1], nrow(bench$conquest_overlap_checks))
  expect_match(row$KeySignal[1], "package-side attention items = 0", fixed = TRUE)
})

test_that("reference_case_benchmark enforces the ConQuest-overlap package-side route", {
  expect_error(
    reference_case_benchmark(
      cases = "synthetic_conquest_overlap_dry_run",
      method = "JML",
      model = "RSM"
    ),
    "requires `method = \"MML\"`",
    fixed = TRUE
  )
  expect_error(
    reference_case_benchmark(
      cases = "synthetic_conquest_overlap_dry_run",
      method = "MML",
      model = "GPCM"
    ),
    "requires `model = \"RSM\"` or `model = \"PCM\"`",
    fixed = TRUE
  )
})

test_that("reference_case_benchmark summary surfaces specialized contract checks", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = c("synthetic_latent_regression_omit", "synthetic_conquest_overlap_dry_run"),
    method = "MML",
    model = "RSM",
    quad_points = 5,
    maxit = 30
  ))

  s <- summary(bench)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_equal(nrow(s$summary), 2)
  expect_true(all(c(
    "synthetic_latent_regression_omit",
    "synthetic_conquest_overlap_dry_run"
  ) %in% s$summary$Case))
  expect_true(is.data.frame(s$validation_scope))
  expect_true(is.data.frame(s$conquest_overlap_checks))
  expect_true(is.data.frame(s$population_policy_checks))
  expect_true(all(c(
    "Package reference check",
    "Latent-regression omission policy",
    "ConQuest-overlap package-side check",
    "External ConQuest validation"
  ) %in% s$validation_scope$Area))
  external_row <- s$validation_scope[
    s$validation_scope$Area == "External ConQuest validation",
    ,
    drop = FALSE
  ]
  expect_identical(as.character(external_row$Status[1]), "not performed")
  expect_match(
    external_row$Interpretation[1],
    "Actual external ConQuest output tables are required",
    fixed = TRUE
  )
  expect_true(nrow(s$conquest_overlap_checks) > 0)
  expect_true(nrow(s$population_policy_checks) > 0)
  expect_true(any(grepl("actual external ConQuest output is still required", s$notes, fixed = TRUE)))

  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "synthetic_conquest_overlap_dry_run", fixed = TRUE)
  expect_match(printed, "Validation scope", fixed = TRUE)
  expect_match(printed, "ConQuest-overlap checks", fixed = TRUE)
  expect_match(printed, "Population-policy checks", fixed = TRUE)
  expect_match(printed, "actual external ConQuest output is still required", fixed = TRUE)
  expect_match(printed, "external validation.\n - Bias checks", fixed = TRUE)
})

test_that("reference_case_benchmark includes latent-regression benchmark case", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_latent_regression",
    method = "MML",
    quad_points = 7,
    maxit = 40
  ))

  expect_s3_class(bench, "mfrm_reference_benchmark")
  expect_true("synthetic_latent_regression" %in% bench$case_summary$Case)
  expect_true(any(grepl("Population:", bench$recovery_checks$Facet, fixed = TRUE)))
  expect_true(any(bench$recovery_checks$Facet == "Population:posterior_shift"))
  expect_true(all(bench$fit_runs$SupportsFormalInference))
  expect_true(all(bench$fit_runs$PopulationModelActive))
  expect_identical(as.character(bench$fit_runs$PosteriorBasis[1]), "population_model")
  expect_match(bench$fit_runs$PopulationFormula[1], "X", fixed = TRUE)
  expect_identical(as.character(bench$fit_runs$PopulationDesignColumns[1]), "(Intercept), X")
  expect_identical(as.character(bench$fit_runs$PopulationXlevelVariables[1]), "")
  expect_identical(as.character(bench$fit_runs$PopulationContrastVariables[1]), "")
  expect_identical(bench$fit_runs$PopulationCoefficientCount[1], 2L)
  expect_true(is.finite(bench$fit_runs$PopulationResidualVariance[1]))
  expect_identical(bench$fit_runs$PopulationIncludedPersons[1], 60L)
  expect_identical(bench$fit_runs$PopulationOmittedPersons[1], 0L)
  expect_identical(bench$fit_runs$PopulationResponseRowsOmitted[1], 0L)
})

test_that("reference_case_benchmark recovers synthetic truth under MML", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_truth",
    method = "MML",
    quad_points = 7,
    maxit = 40
  ))

  expect_true(all(bench$recovery_checks$Status %in% c("Pass", "Warn", "Fail")))
  expect_true(min(bench$recovery_checks$Correlation, na.rm = TRUE) > 0.95)
  expect_true(max(bench$recovery_checks$MeanAbsoluteDeviation, na.rm = TRUE) < 0.30)
  expect_true(all(bench$fit_runs$SupportsFormalInference))
  expect_true(all(!bench$fit_runs$PopulationModelActive))
  expect_identical(as.character(bench$fit_runs$PosteriorBasis[1]), "legacy_mml")
  expect_identical(bench$fit_runs$PopulationOmittedPersons[1], 0L)
  expect_identical(bench$fit_runs$PopulationResponseRowsOmitted[1], 0L)
})

test_that("reference_case_benchmark recovers latent-regression synthetic case under MML", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_latent_regression",
    method = "MML",
    quad_points = 7,
    maxit = 40
  ))

  intercept_row <- bench$recovery_checks[bench$recovery_checks$Facet == "Population:(Intercept)", , drop = FALSE]
  slope_row <- bench$recovery_checks[bench$recovery_checks$Facet == "Population:X", , drop = FALSE]
  sigma_row <- bench$recovery_checks[bench$recovery_checks$Facet == "Population:sigma2", , drop = FALSE]
  crit_row <- bench$recovery_checks[bench$recovery_checks$Facet == "Criterion", , drop = FALSE]
  shift_row <- bench$recovery_checks[bench$recovery_checks$Facet == "Population:posterior_shift", , drop = FALSE]

  expect_equal(nrow(intercept_row), 1)
  expect_equal(nrow(slope_row), 1)
  expect_equal(nrow(sigma_row), 1)
  expect_equal(nrow(crit_row), 1)
  expect_equal(nrow(shift_row), 1)
  expect_true(intercept_row$MeanAbsoluteDeviation[1] < 0.20)
  expect_true(slope_row$MeanAbsoluteDeviation[1] < 0.35)
  expect_true(sigma_row$MeanAbsoluteDeviation[1] < 0.35)
  expect_true(crit_row$Correlation[1] > 0.95)
  expect_identical(as.character(shift_row$Status[1]), "Pass")
})

test_that("reference_case_benchmark recovers synthetic GPCM case under MML", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_gpcm",
    method = "MML",
    model = "GPCM",
    quad_points = 11,
    maxit = 120
  ))

  slope_row <- bench$recovery_checks[bench$recovery_checks$Facet == "GPCM:slopes", , drop = FALSE]
  step_row <- bench$recovery_checks[bench$recovery_checks$Facet == "GPCM:steps", , drop = FALSE]
  crit_row <- bench$recovery_checks[bench$recovery_checks$Facet == "Criterion", , drop = FALSE]

  expect_equal(nrow(slope_row), 1)
  expect_equal(nrow(step_row), 1)
  expect_equal(nrow(crit_row), 1)
  expect_true(slope_row$Correlation[1] > 0.95)
  expect_true(slope_row$MeanAbsoluteDeviation[1] < 0.15)
  expect_true(step_row$Correlation[1] > 0.98)
  expect_true(step_row$MeanAbsoluteDeviation[1] < 0.20)
  expect_true(crit_row$Correlation[1] > 0.98)
})

test_that("reference_case_benchmark enforces the GPCM benchmark contract", {
  expect_error(
    reference_case_benchmark(cases = "synthetic_gpcm", method = "MML"),
    "requires `model = \"GPCM\"`",
    fixed = TRUE
  )
  expect_error(
    reference_case_benchmark(cases = "synthetic_gpcm", method = "JML", model = "GPCM"),
    "validated only for `method = \"MML\"`",
    fixed = TRUE
  )
})

test_that("reference_case_benchmark captures pair stability and summary output", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "study1_itercal_pair",
    method = "MML",
    quad_points = 5,
    maxit = 30
  ))

  expect_true(all(c("Criterion", "Rater", "OverallFit") %in% bench$pair_checks$Facet))
  expect_true(all(c("Criterion", "Rater") %in% bench$linking_checks$Facet))
  expect_identical(
    as.character(bench$linking_checks$Status[bench$linking_checks$Facet == "Rater"][1]),
    "Pass"
  )
  expect_identical(
    as.character(bench$linking_checks$Status[bench$linking_checks$Facet == "Criterion"][1]),
    "Warn"
  )
  criterion_row <- bench$pair_checks[bench$pair_checks$Facet == "Criterion", , drop = FALSE]
  expect_true(criterion_row$Pearson[1] > 0.95)
  expect_true(criterion_row$MeanAbsoluteDifference[1] < 0.10)

  s <- summary(bench)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_true(is.data.frame(s$overview))
  expect_identical(as.character(s$overview$Class[1]), "mfrm_reference_benchmark")

  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "mfrmr Reference Case Check Summary", fixed = TRUE)
})

test_that("reference_case_benchmark verifies bias-contract identities", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_bias_contract",
    method = "MML",
    quad_points = 7,
    maxit = 40
  ))

  expect_true(nrow(bench$bias_checks) >= 4)
  stable_rows <- bench$bias_checks[
    bench$bias_checks$Metric %in% c("BiasDFIdentity", "LocalMeasureIdentity", "PairContrastIdentity"),
    ,
    drop = FALSE
  ]
  expect_true(nrow(stable_rows) == 3)
  expect_true(all(stable_rows$Status == "Pass"))
  expect_true(max(stable_rows$MaxError, na.rm = TRUE) < 1e-8)
  expect_true(any(bench$source_profile$RuleID == "bias_pairwise_welch"))
})

test_that("reference_case_benchmark reports exploratory precision under JML", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_truth",
    method = "JML",
    maxit = 30
  ))

  expect_true(all(bench$fit_runs$PrecisionTier == "exploratory"))
  expect_true(all(!bench$fit_runs$SupportsFormalInference))
})

test_that("reference_case_benchmark handles single-case runs without warnings", {
  expect_no_warning(
    reference_case_benchmark(
      cases = "synthetic_truth",
      method = "MML",
      quad_points = 5,
      maxit = 20
    )
  )
})

test_that("reference_case_benchmark carries MML engine metadata", {
  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_truth",
    method = "MML",
    quad_points = 5,
    maxit = 20,
    mml_engine = "em"
  ))

  expect_identical(bench$settings$mml_engine, "em")
  expect_true(all(bench$fit_runs$MMLEngineRequested == "em"))
  expect_true(all(bench$fit_runs$MMLEngineUsed == "em"))
  expect_true(all(is.finite(bench$fit_runs$EMIterations)))
})

test_that("reference_case_benchmark does not pass bias-contract cases when checks are missing", {
  local_mocked_bindings(
    estimate_bias = function(...) {
      structure(list(table = data.frame()), class = "mfrm_bias")
    },
    .package = "mfrmr"
  )

  bench <- suppressWarnings(reference_case_benchmark(
    cases = "synthetic_bias_contract",
    method = "MML",
    quad_points = 5,
    maxit = 20
  ))

  row <- bench$case_summary[bench$case_summary$Case == "synthetic_bias_contract", , drop = FALSE]
  expect_identical(as.character(row$Status[1]), "Warn")
  expect_match(row$KeySignal[1], "No bias checks were produced", fixed = TRUE)
})
