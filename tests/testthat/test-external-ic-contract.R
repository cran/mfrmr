test_that("the shared common IC panel has one exact Person-basis formula", {
  panel <- mfrmr:::mfrm_ic_common_panel(
    deviance = 800,
    npar = 12L,
    persons = 100L
  )
  expect_equal(panel$AIC, 824)
  expect_equal(panel$BIC, 800 + log(100) * 12, tolerance = 1e-12)
  expect_equal(
    panel$SABIC,
    800 + log((100 + 2) / 24) * 12,
    tolerance = 1e-12
  )
  expect_identical(panel$BICFormula, "bic_person_count")
  expect_identical(panel$SABICFormula, "sclove_n_plus_2_over_24")
  expect_true(panel$SABICSelectable)
  expect_error(
    mfrmr:::mfrm_ic_common_panel(800, 12.5, 100),
    "non-negative integer"
  )
  expect_error(
    mfrmr:::mfrm_ic_common_panel(800, 12, 0),
    "positive integer"
  )
})

test_that("the public TAM importer fails closed outside unidimensional MML", {
  tam_jml <- structure(list(), class = "tam.jml")
  expect_error(
    mfrmr:::.tam_validate_import_scope(tam_jml),
    "does not import `tam.jml` objects as MML",
    fixed = TRUE
  )

  tam_two_dimensional <- structure(
    list(ndim = 2L, beta = matrix(0, nrow = 1L, ncol = 2L)),
    class = "tam.mml"
  )
  expect_error(
    mfrmr:::.tam_validate_import_scope(tam_two_dimensional),
    "supports unidimensional TAM MML only",
    fixed = TRUE
  )

  tam_one_dimensional <- structure(
    list(ndim = 1L, beta = matrix(0, nrow = 1L, ncol = 1L)),
    class = "tam.mml"
  )
  expect_identical(
    mfrmr:::.tam_validate_import_scope(tam_one_dimensional),
    list(Method = "MML", Dimensions = 1L)
  )
})

test_that("repository external IC fixtures preserve native and common panels", {
  candidates <- c(
    file.path("inst", "validation"),
    testthat::test_path("..", "..", "inst", "validation")
  )
  validation_dir <- candidates[dir.exists(candidates)][1]
  testthat::skip_if(
    length(validation_dir) == 0L || is.na(validation_dir),
    "Repository-only validation files are not installed with the package."
  )
  env <- new.env(parent = globalenv())
  sys.source(
    file.path(validation_dir, "external-ic-normalizer-0.2.3.R"),
    envir = env
  )
  sys.source(
    file.path(validation_dir, "external-ic-audit-0.2.3.R"),
    envir = env
  )
  audit <- env$mfrmr_run_external_ic_fixture_audit(
    pkg_dir = normalizePath(
      file.path(validation_dir, "..", ".."),
      winslash = "/",
      mustWork = TRUE
    )
  )
  expect_identical(audit$status, "ok")
  expect_identical(audit$contract_version, "mfrmr_external_ic_v1")
  expect_true(audit$ready_comparison$comparable)
  expect_false(audit$mismatch_comparison$comparable)
  tam_record <- audit$records[["EXTIC-001"]]$record
  expect_true(tam_record$NativeABICFormulaVerified)
  expect_false(isTRUE(all.equal(
    tam_record$NativeABIC,
    tam_record$CommonSABIC,
    tolerance = 1e-12
  )))
})

test_that("the repository ConQuest adapter audits objective, free dimension, and Persons", {
  candidates <- c(
    file.path("inst", "validation"),
    testthat::test_path("..", "..", "inst", "validation")
  )
  validation_dir <- candidates[dir.exists(candidates)][1]
  testthat::skip_if(
    length(validation_dir) == 0L || is.na(validation_dir),
    "Repository-only validation files are not installed with the package."
  )
  env <- new.env(parent = globalenv())
  sys.source(
    file.path(validation_dir, "external-ic-normalizer-0.2.3.R"),
    envir = env
  )

  out_dir <- file.path(tempdir(), "mfrmr-conquest-ic-adapter")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  final_vector <- c(
    0.379656, 0.881274, 0.422549,
    -0.847056, -0.653147, -0.039260, -0.122987, 0.534787
  )
  history <- data.frame(
    RowLabels = c("row_1", "row_42"),
    `Run Number` = c(1L, 1L),
    Iteration = c(1L, 42L),
    LogLikelihood = c(483.141322, 424.739512),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  history_columns <- c(
    "Dim 1 Var 1", "Dim 2 Var 2", "wvar 1 1",
    "xsi 1", "xsi 2", "xsi 3", "xsi 4", "xsi 5"
  )
  for (index in seq_along(history_columns)) {
    history[[history_columns[index]]] <- c(0, final_vector[index])
  }
  parameter <- data.frame(
    P = 1:5,
    Estimate = final_vector[4:8],
    Label = paste("item", sprintf("i%03d", 1:5)),
    stringsAsFactors = FALSE
  )
  regression <- data.frame(
    Dimension = c(1L, 1L),
    Regressor = 1:2,
    Estimate = final_vector[1:2]
  )
  covariance <- data.frame(
    Dim1 = 1L,
    Dim2 = 1L,
    Covariance = final_vector[3]
  )
  cases <- data.frame(
    SeqNum = 1:60,
    PID = sprintf("P%03d", 1:60),
    EAP_1 = 0,
    weight_raw = 1,
    weight_scaled = 1,
    stringsAsFactors = FALSE
  )
  paths <- file.path(
    out_dir,
    c("history.csv", "parameters.csv", "regression.csv", "covariance.csv",
      "cases.csv")
  )
  utils::write.csv(history, paths[1], row.names = FALSE)
  utils::write.csv(parameter, paths[2], row.names = FALSE)
  utils::write.csv(regression, paths[3], row.names = FALSE)
  utils::write.csv(covariance, paths[4], row.names = FALSE)
  utils::write.csv(cases, paths[5], row.names = FALSE)

  adapter_args <- list(
    history_file = paths[1],
    parameter_file = paths[2],
    regression_file = paths[3],
    covariance_file = paths[4],
    case_file = paths[5],
    engine_version = "5.47.5 Demonstration Version",
    run_date = as.Date("2026-07-27"),
    run_id = "cq-synthetic-31-node",
    model_id = "EXT-CQ-BINARY-PILOT",
    quadrature_nodes = 31L,
    expected_person_ids = cases$PID,
    observation_set_id = "cq-synthetic-60-persons",
    likelihood_basis_id = "conquest-binary-mml-observed-v1",
    constraint_basis_id = "conquest-default-item-sum-zero-v1",
    integration_comparison_id = "cq-31-node-review-v1",
    convergence_status = "pass",
    convergence_evidence_id = "table1-deviance-change-termination-review",
    integration_stability_status = "pass",
    candidate_id = "working-tree-pilot"
  )
  record <- do.call(env$mfrmr_external_ic_from_conquest, adapter_args)
  expect_s3_class(record, "mfrmr_external_ic_record")
  expect_true(record$record$ComparisonReady)
  expect_equal(record$record$Deviance, 424.739512, tolerance = 1e-12)
  expect_identical(record$record$Npar, 8L)
  expect_identical(record$record$Persons, 60L)
  expect_equal(record$record$CommonAIC, 440.739512, tolerance = 1e-12)
  expect_identical(record$audit$NparHistory, 8L)
  expect_identical(record$audit$NparExports, 8L)
  expect_true(record$audit$FreeVectorMatched)
  expect_true(record$audit$PersonIdsMatched)
  expect_identical(record$audit$ObjectiveNativeHeader, "LogLikelihood")
  expect_match(record$audit$ObjectiveInterpretation, "deviance_by_ConQuest")
  expect_equal(record$audit$InternalHandoffTolerance, 1e-6)
  expect_true(is.na(record$audit$ObjectiveExportResolution))
  expect_false(record$audit$ObjectiveExportResolutionEstablished)
  expect_identical(record$audit$ObjectiveExportRoundingRule, "unknown")
  expect_false(record$audit$HandoffToleranceIsCrossEngineTolerance)
  expect_true(all(nzchar(record$audit$ArtifactFingerprints$Digest)))

  legacy_tolerance <- adapter_args
  legacy_tolerance$export_tolerance <- 1e-7
  legacy_record <- NULL
  expect_warning(
    legacy_record <- do.call(
      env$mfrmr_external_ic_from_conquest,
      legacy_tolerance
    ),
    "deprecated"
  )
  expect_equal(legacy_record$audit$InternalHandoffTolerance, 1e-7)
  expect_true(is.na(legacy_record$audit$ObjectiveExportResolution))

  conflicting_tolerances <- adapter_args
  conflicting_tolerances$handoff_tolerance <- 1e-5
  conflicting_tolerances$export_tolerance <- 1e-7
  expect_error(
    do.call(
      env$mfrmr_external_ic_from_conquest,
      conflicting_tolerances
    ),
    "Supply only `handoff_tolerance`",
    fixed = TRUE
  )

  unaudited_version <- adapter_args
  unaudited_version$engine_version <- "5.48.0"
  expect_error(
    do.call(env$mfrmr_external_ic_from_conquest, unaudited_version),
    "audited for version 5.47.5 only",
    fixed = TRUE
  )

  missing_convergence_evidence <- adapter_args
  missing_convergence_evidence$convergence_evidence_id <- NA_character_
  expect_error(
    do.call(
      env$mfrmr_external_ic_from_conquest,
      missing_convergence_evidence
    ),
    "convergence_evidence_id",
    fixed = TRUE
  )

  bad_history <- history
  bad_history$unexpected_free_coordinate <- 0
  bad_history_file <- file.path(out_dir, "bad-history.csv")
  utils::write.csv(bad_history, bad_history_file, row.names = FALSE)
  mismatched_dimension <- adapter_args
  mismatched_dimension$history_file <- bad_history_file
  expect_error(
    do.call(env$mfrmr_external_ic_from_conquest, mismatched_dimension),
    "free dimension disagrees",
    fixed = TRUE
  )

  weighted_cases <- cases
  weighted_cases$weight_raw[1] <- 2
  weighted_case_file <- file.path(out_dir, "weighted-cases.csv")
  utils::write.csv(weighted_cases, weighted_case_file, row.names = FALSE)
  nonunit_weight <- adapter_args
  nonunit_weight$case_file <- weighted_case_file
  expect_error(
    do.call(env$mfrmr_external_ic_from_conquest, nonunit_weight),
    "requires unit case weights",
    fixed = TRUE
  )

  mismatched_persons <- adapter_args
  mismatched_persons$expected_person_ids[1] <- "MISSING_PERSON"
  expect_error(
    do.call(env$mfrmr_external_ic_from_conquest, mismatched_persons),
    "do not exactly match",
    fixed = TRUE
  )
})

test_that("the TAM adapter preserves native aBIC without calling it SABIC", {
  testthat::skip_if_not_installed("TAM")
  candidates <- c(
    file.path("inst", "validation"),
    testthat::test_path("..", "..", "inst", "validation")
  )
  validation_dir <- candidates[dir.exists(candidates)][1]
  testthat::skip_if(
    length(validation_dir) == 0L || is.na(validation_dir),
    "Repository-only validation files are not installed with the package."
  )
  env <- new.env(parent = globalenv())
  sys.source(
    file.path(validation_dir, "external-ic-normalizer-0.2.3.R"),
    envir = env
  )

  set.seed(2701L)
  response <- matrix(sample(0:2, 240L, replace = TRUE), nrow = 60L)
  colnames(response) <- paste0("I", seq_len(ncol(response)))
  tam_fit <- suppressMessages(suppressWarnings(TAM::tam.mml(
    resp = response,
    irtmodel = "PCM",
    control = list(progress = FALSE)
  )))
  record <- env$mfrmr_external_ic_from_tam(
    tam_fit,
    run_id = "tam-adapter-test",
    model_id = "TAM-1D-PCM",
    observation_set_id = "tam-adapter-data",
    likelihood_basis_id = "tam-pcm-observed-data-v1",
    constraint_basis_id = "tam-default-test",
    integration_comparison_id = "tam-test-qmc-review",
    convergence_status = "pass",
    integration_stability_status = "pass"
  )
  expect_true(record$record$ComparisonReady)
  expect_true(record$record$NativeAICFormulaVerified)
  expect_true(record$record$NativeBICFormulaVerified)
  expect_true(record$record$NativeABICFormulaVerified)
  expect_identical(
    record$record$NativeABICFormula,
    "tam_deviance_plus_log_n_minus_2_over_24_k"
  )
  expect_false(isTRUE(all.equal(
    record$record$NativeABIC,
    record$record$CommonSABIC,
    tolerance = 1e-12
  )))

  stochastic_fit <- tam_fit
  stochastic_fit$control$QMC <- FALSE
  stochastic_fit$control$seed <- 1977L
  stochastic_record <- env$mfrmr_external_ic_from_tam(
    stochastic_fit,
    run_id = "tam-stochastic-id-test",
    model_id = "TAM-1D-PCM",
    observation_set_id = "tam-adapter-data",
    likelihood_basis_id = "tam-pcm-observed-data-v1",
    constraint_basis_id = "tam-default-test",
    integration_comparison_id = "tam-stochastic-review"
  )
  expect_match(
    stochastic_record$record$IntegrationEvaluationId,
    "qmc=false:nnodes=[0-9]+:seed=1977"
  )
})
