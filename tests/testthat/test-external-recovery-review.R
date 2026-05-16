source_external_recovery_review_helpers <- function() {
  script <- testthat::test_path("..", "..", "inst", "validation", "external-recovery-audit.R")
  if (!file.exists(script)) {
    script <- system.file("validation", "external-recovery-audit.R", package = "mfrmr")
  }
  env <- new.env(parent = globalenv())
  source(script, local = env)
  env
}

make_external_recovery_review_fixture <- function(include_review_groups = TRUE) {
  base_dir <- tempfile("mfrmr-external-recovery-")
  dir.create(base_dir, recursive = TRUE)

  write_fixture_csv <- function(relative_path, data) {
    path <- file.path(base_dir, relative_path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(data, path, row.names = FALSE)
  }

  write_fixture_csv(
    "analysis/dataset_manifest.csv",
    data.frame(
      DesignPattern = c("baseline", "weak_bridge"),
      Rows = c(120, 96),
      ObservedDensity = c(0.95, 0.87)
    )
  )
  write_fixture_csv(
    "analysis/smoke_check_summary.csv",
    data.frame(Check = c("manifest", "engine_status"), Passed = c(TRUE, TRUE))
  )
  write_fixture_csv(
    "analysis/engine_status_summary.csv",
    data.frame(
      Engine = c("R", "Python"),
      Runs = c(2, 2),
      ErrorRuns = c(0, 0),
      ConvergenceRate = c(1, 1)
    )
  )

  agreement <- data.frame(
    Source = c("level_recovery", "step_recovery", "fit_statistics"),
    Metric = c("Measure", "Step", "Infit"),
    EnginePair = c("R_vs_Python", "R_vs_Python", "R_vs_Python"),
    Groups = c(4, 3, 5),
    ComparedRows = c(40, 12, 60),
    MaxRMSE = c(0.001, 0.002, 0.02),
    PracticalOrBetterGroups = c(4, 3, 4)
  )
  if (include_review_groups) {
    agreement$ReviewGroups <- c(0, 0, 1)
  }
  write_fixture_csv(
    file.path("analysis", paste0("engine_", "par", "ity", "_overview.csv")),
    agreement
  )

  write_fixture_csv(
    "analysis/key_findings.csv",
    data.frame(
      Priority = c("high", "medium"),
      Rank = c(1, 1),
      Domain = c("Recovery risk", "Fit risk"),
      Metric = c("MeanRMSE", "MisfitRate"),
      Value = c(0.41, 0.12),
      Message = c("Recovery risk | MeanRMSE=0.41", "Fit risk | MisfitRate=0.12")
    )
  )
  write_fixture_csv(
    "analysis/key_findings_counts.csv",
    data.frame(Priority = c("high", "medium"), Domain = c("Recovery risk", "Fit risk"), Rows = c(1, 1))
  )
  write_fixture_csv(
    "sample_size_dstudy/analysis/sample_size_decision_summary.csv",
    data.frame(SampleSizeFacet = c("Person", "Rater"), CandidateN = c(60, 8))
  )
  write_fixture_csv(
    "sample_size_dstudy/analysis/sample_size_classification_summary.csv",
    data.frame(ClassificationRisk = c("low", "moderate"), Cells = c(2, 1))
  )

  base_dir
}

test_that("external recovery helper summarizes complete local outputs", {
  helpers <- source_external_recovery_review_helpers()
  fixture <- make_external_recovery_review_fixture()

  review <- helpers$mfrmr_review_external_recovery_simulation(fixture)
  review_summary <- helpers$summary.mfrmr_external_recovery_review(review)

  expect_s3_class(review, "mfrmr_external_recovery_review")
  expect_equal(review$decision$EvidenceStatus, "review")
  expect_equal(review$overview$SmokeChecksPassed, 2)
  expect_equal(review$overview$ErrorRuns, 0)
  expect_equal(review$overview$MinConvergenceRate, 1)
  expect_true(all(review$schema_status$SchemaOK))
  expect_true(all(nzchar(review$file_status$MD5[review$file_status$Exists])))
  expect_equal(
    review$agreement_summary$ReviewGroups[review$agreement_summary$Source == "level_recovery"],
    0
  )
  expect_equal(
    review$agreement_summary$ReviewGroups[review$agreement_summary$Source == "fit_statistics"],
    1
  )
  expect_true(review$sample_size_summary$Available)
  expect_equal(nrow(review$top_findings), 1)
  expect_named(
    review_summary,
    c("decision", "schema_status", "overview", "agreement_summary", "sample_size_summary", "top_findings")
  )
  expect_output(helpers$print.summary.mfrmr_external_recovery_review(review_summary), "Schema status")
})

test_that("external recovery helper flags missing required outputs", {
  helpers <- source_external_recovery_review_helpers()
  fixture <- tempfile("mfrmr-external-recovery-missing-")
  dir.create(fixture, recursive = TRUE)

  review <- helpers$mfrmr_review_external_recovery_simulation(fixture)

  expect_equal(review$decision$EvidenceStatus, "concern")
  expect_true(any(review$file_status$Required & !review$file_status$Exists))
  expect_false(any(review$schema_status$SchemaOK[review$schema_status$Required]))
})

test_that("external recovery helper does not treat incomplete agreement tables as OK", {
  helpers <- source_external_recovery_review_helpers()
  fixture <- make_external_recovery_review_fixture(include_review_groups = FALSE)

  review <- helpers$mfrmr_review_external_recovery_simulation(fixture)

  expect_equal(review$decision$EvidenceStatus, "concern")
  expect_false(all(review$schema_status$SchemaOK[review$schema_status$Required]))
  expect_match(review$decision$Interpretation, "ReviewGroups", fixed = TRUE)
  expect_true(all(is.na(review$agreement_summary$ReviewGroups)))
})
