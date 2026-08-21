gpcm_estimator_asymptotics_runner <- function() {
  path <- testthat::test_path(
    "..", "..", "inst", "validation",
    "gpcm-estimator-asymptotics-0.2.3.R"
  )
  skip_if_not(file.exists(path),
              "repository-internal validation artifacts are excluded")
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  env
}

test_that("GPCM asymptotic sequences separate Persons from exposure", {
  env <- gpcm_estimator_asymptotics_runner()
  smoke <- env$mfrmr_gas_manifest("smoke")
  pilot <- env$mfrmr_gas_manifest("pilot")
  membership <- env$mfrmr_gas_sequence_membership("smoke")

  expect_equal(nrow(smoke), 10L)
  expect_equal(nrow(pilot), 200L)
  expect_equal(length(unique(smoke$CellId)), 5L)
  expect_setequal(smoke$Method, c("JML", "MML"))
  expect_true(all(table(smoke$CellId) == 2L))
  expect_true(all(!smoke$EstimatorSelectionAuthorized))
  expect_true(all(!smoke$BiasCorrectionAuthorized))
  expect_true(all(is.na(smoke$BayesianEstimatorRequired)))
  expect_true(all(!smoke$ConfirmationAuthorized))
  expect_true(all(
    pilot$IncidentalBiasDecision ==
      "pilot_planned_no_prespecified_decision_rule"
  ))

  fixed_exposure <- membership[
    membership$Sequence == "persons_increase_exposure_fixed", , drop = FALSE
  ]
  fixed_persons <- membership[
    membership$Sequence == "exposure_increases_persons_fixed", , drop = FALSE
  ]
  cells <- env$mfrmr_gas_cells("smoke")
  first <- merge(fixed_exposure, cells, by = "CellId")
  second <- merge(fixed_persons, cells, by = "CellId")
  expect_identical(sort(first$NPersons), c(24L, 48L, 96L))
  expect_true(all(first$ObservationsPerPerson == 8L))
  expect_true(all(second$NPersons == 48L))
  expect_identical(sort(second$ObservationsPerPerson), c(8L, 16L, 24L))
  expect_equal(length(intersect(fixed_exposure$CellId,
                                fixed_persons$CellId)), 1L)
})

test_that("nested GPCM cells preserve truth and realized responses", {
  env <- gpcm_estimator_asymptotics_runner()
  manifest <- env$mfrmr_gas_manifest("smoke")
  base <- env$mfrmr_gas_generate_base(manifest, replicate = 1L)
  cells <- env$mfrmr_gas_cells("smoke")
  generated <- lapply(seq_len(nrow(cells)), function(index) {
    env$mfrmr_gas_subset_cell(base, cells[index, , drop = FALSE])
  })
  audits <- lapply(seq_along(generated), function(index) {
    env$mfrmr_gas_design_audit(
      generated[[index]], cells[index, , drop = FALSE]
    )
  })
  expect_true(all(vapply(audits, function(x) {
    isTRUE(x$DesignContractPassed)
  }, logical(1))))
  expect_true(all(vapply(generated, function(x) {
    identical(attr(x, "mfrm_truth")$facets, base$truth$facets) &&
      identical(attr(x, "mfrm_truth")$step_table, base$truth$step_table) &&
      identical(attr(x, "mfrm_truth")$slope_table, base$truth$slope_table) &&
      identical(attr(x, "mfrm_truth")$population, base$truth$population)
  }, logical(1))))
  expect_equal(base$truth$population$coefficients[["(Intercept)"]], 0)
  expect_equal(base$truth$population$sigma2, 1)
  expect_match(
    base$truth$population$basis,
    "not_realized_sample_moments",
    fixed = TRUE
  )

  row_key <- function(data) {
    paste(data$Person, data$Rater, data$Criterion, data$Score, sep = "\r")
  }
  low_n <- generated[[1L]]
  middle_n <- generated[[2L]]
  high_n <- generated[[3L]]
  medium_exposure <- generated[[4L]]
  high_exposure <- generated[[5L]]
  expect_true(all(row_key(low_n) %in% row_key(middle_n)))
  expect_true(all(row_key(middle_n) %in% row_key(high_n)))
  expect_true(all(row_key(middle_n) %in% row_key(medium_exposure)))
  expect_true(all(row_key(medium_exposure) %in% row_key(high_exposure)))
})

test_that("warnings and errors remain explicit in asymptotic runs", {
  env <- gpcm_estimator_asymptotics_runner()
  warning_result <- env$mfrmr_gas_capture({
    warning("retained asymptotic warning")
    1L
  })
  error_result <- env$mfrmr_gas_capture(stop("retained asymptotic error"))
  expect_identical(warning_result$value, 1L)
  expect_identical(warning_result$warnings, "retained asymptotic warning")
  expect_s3_class(error_result$value, "error")
  expect_identical(conditionMessage(error_result$value),
                   "retained asymptotic error")
})

test_that("component aggregation cannot present cancelling pooled bias", {
  env <- gpcm_estimator_asymptotics_runner()
  fixture <- data.frame(
    Sequence = "fixed_exposure",
    XName = "NPersons",
    XValue = 48L,
    Method = "JML",
    Replicate = c(1L, 1L, 2L, 2L),
    ErrorAligned = c(-0.5, 0.5, -0.25, 0.25),
    InferenceEligible = FALSE,
    stringsAsFactors = FALSE
  )
  summary <- env$mfrmr_gas_aggregate(
    fixture, c("Sequence", "XName", "XValue", "Method"),
    error_column = "ErrorAligned"
  )
  expect_false("Bias" %in% names(summary))
  expect_equal(summary$Replicates, 2L)
  expect_equal(summary$MAE, 0.375)
  expect_equal(summary$RMSE, sqrt(mean(fixture$ErrorAligned^2)))
  expect_equal(summary$InferenceEligibleRate, 0)
})

test_that("pilot summaries retain Monte Carlo and paired uncertainty", {
  env <- gpcm_estimator_asymptotics_runner()
  replicate_metrics <- data.frame(
    Replicate = rep(1:2, 2),
    CellId = rep(c("low", "high"), each = 2),
    NPersons = rep(c(60L, 240L), each = 2),
    ObservationsPerPerson = 8L,
    Method = "JML",
    Component = "step",
    Coordinates = 3L,
    RMSE = c(0.4, 0.6, 0.2, 0.3),
    MAE = c(0.3, 0.5, 0.1, 0.2),
    stringsAsFactors = FALSE
  )
  component <- env$mfrmr_gas_component_mc(replicate_metrics)
  low <- component[component$CellId == "low", , drop = FALSE]
  expect_equal(low$Replicates, 2L)
  expect_equal(low$MeanRMSE, 0.5)
  expect_equal(low$MCSERMSE, 0.1)

  membership <- data.frame(
    CellId = c("low", "high"),
    Sequence = "persons_increase_exposure_fixed",
    XName = "NPersons",
    XValue = c(60L, 240L),
    stringsAsFactors = FALSE
  )
  contrast <- env$mfrmr_gas_endpoint_contrasts(
    replicate_metrics, membership
  )
  expect_equal(contrast$MeanRMSEDifference, -0.25)
  expect_equal(contrast$MCSERMSEDifference, 0.05)
  expect_equal(contrast$RMSEImprovementRate, 1)
  expect_true(contrast$Lower95RMSEDifference < -0.25)
  expect_true(contrast$Upper95RMSEDifference > -0.25)

  method_metrics <- rbind(
    replicate_metrics[replicate_metrics$CellId == "low", , drop = FALSE],
    transform(
      replicate_metrics[replicate_metrics$CellId == "low", , drop = FALSE],
      Method = "MML", RMSE = c(0.3, 0.4), MAE = c(0.2, 0.3)
    )
  )
  method_contrast <- env$mfrmr_gas_method_contrasts(method_metrics)
  expect_equal(method_contrast$MeanRMSEDifference, 0.15)
  expect_equal(method_contrast$MMLLowerRMSERate, 1)

  mml_population <- transform(
    method_metrics[method_metrics$Method == "MML", , drop = FALSE],
    Component = "population_variance", RMSE = c(0.2, 0.1)
  )
  shared_only <- env$mfrmr_gas_method_contrasts(
    rbind(method_metrics, mml_population)
  )
  expect_identical(unique(shared_only$Component), "step")

  coordinates <- data.frame(
    Replicate = rep(1:2, 2),
    CellId = "low",
    NPersons = 60L,
    ObservationsPerPerson = 8L,
    Method = "JML",
    ParameterType = "facet",
    Facet = "Rater",
    Level = rep(c("R01", "R02"), each = 2),
    Subparameter = "measure",
    ErrorAligned = c(-0.5, -0.25, 0.5, 0.25),
    stringsAsFactors = FALSE
  )
  coordinate <- env$mfrmr_gas_coordinate_mc(
    coordinates,
    keys = c(
      "CellId", "NPersons", "ObservationsPerPerson", "Method",
      "ParameterType", "Facet", "Level", "Subparameter"
    ),
    error_column = "ErrorAligned"
  )
  expect_equal(nrow(coordinate), 2L)
  expect_equal(sort(coordinate$MeanError), c(-0.375, 0.375))
})

test_that("the asymptotic contract has no machine-identity gate", {
  env <- gpcm_estimator_asymptotics_runner()
  plan <- env$mfrmr_run_gpcm_estimator_asymptotics("smoke")
  expect_s3_class(plan, "mfrmr_gpcm_estimator_asymptotics_plan")
  expect_equal(nrow(plan$manifest), 10L)

  path <- testthat::test_path(
    "..", "..", "inst", "validation",
    "gpcm-estimator-asymptotics-0.2.3.R"
  )
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_false(grepl("sha-?256|md5|digest::|serialize\\(", source,
                     ignore.case = TRUE))
})

test_that("the smoke record preserves its non-decision", {
  record <- testthat::test_path(
    "..", "..", "inst", "validation",
    "gpcm-estimator-asymptotics-smoke-record-0.2.3.md"
  )
  skip_if_not(file.exists(record),
              "repository-internal validation artifacts are excluded")
  text <- paste(readLines(record, warn = FALSE), collapse = "\n")
  expect_match(text, "All ten fits returned objects", fixed = TRUE)
  expect_match(text, "Every fit remained `InferenceReady=FALSE`", fixed = TRUE)
  expect_match(
    text,
    "IncidentalBiasDecision = not_assigned_replicated_pilot_required",
    fixed = TRUE
  )
  expect_match(text, "BayesianEstimatorRequired = NA", fixed = TRUE)
  expect_match(text, "No Mac library", fixed = TRUE)
})

test_that("the pilot record distinguishes evidence from authorization", {
  record <- testthat::test_path(
    "..", "..", "inst", "validation",
    "gpcm-estimator-asymptotics-pilot-record-0.2.3.md"
  )
  skip_if_not(file.exists(record),
              "repository-internal validation artifacts are excluded")
  text <- paste(readLines(record, warn = FALSE), collapse = "\n")
  expect_match(text, "MML optimizer convergence was 100/100", fixed = TRUE)
  expect_match(text, "JML optimizer convergence was 98/100", fixed = TRUE)
  expect_match(
    text,
    "generating population truth was mean zero and variance one",
    fixed = TRUE
  )
  expect_match(
    text,
    "pilot_completed_no_prespecified_decision_rule",
    fixed = TRUE
  )
  expect_match(text, "EstimatorSelectionAuthorized = FALSE", fixed = TRUE)
  expect_match(text, "BayesianEstimatorRequired = NA", fixed = TRUE)
  expect_false(grepl("SHA[^\n]*=|MD5[^\n]*=", text))
})
