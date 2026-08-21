test_that("Draft.63 registers owner-specific GPCM evidence identities", {
  pkg_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  identity_path <- file.path(
    pkg_root, "inst", "validation", "gpcm-model-identity-contract-0.2.3.csv"
  )
  gate_path <- file.path(
    pkg_root, "inst", "validation", "release-gate-spec-0.2.3.md"
  )
  checklist_path <- file.path(
    pkg_root, "inst", "validation", "release-evidence-checklist-0.2.3.csv"
  )
  execution_path <- file.path(
    pkg_root, "inst", "validation",
    "gpcm-owner-specific-execution-contract-0.2.3.md"
  )
  record_path <- file.path(
    pkg_root, "inst", "validation",
    "gpcm-owner-specific-pilot-record-0.2.3.md"
  )
  skip_if_not(all(file.exists(c(
    identity_path, gate_path, checklist_path, execution_path, record_path
  ))), "repository-internal validation artifacts are excluded")
  expect_true(all(file.exists(c(
    identity_path, gate_path, checklist_path, execution_path, record_path
  ))))

  identity <- utils::read.csv(
    identity_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  required <- c(
    "ModelStratum", "ScenarioId", "SlopeOwner", "StepOwner",
    "SlopeComposition", "LatentDimensionCount", "Estimator",
    "AbilityScaleContract", "ImplementationStatus", "EvidenceStatus",
    "ClaimUse"
  )
  expect_true(all(required %in% names(identity)))
  expect_identical(anyDuplicated(identity$ModelStratum), 0L)

  criterion <- identity[
    identity$ScenarioId == "NUM-GPCM-ALIGN-CRITERION", , drop = FALSE
  ]
  rater <- identity[
    identity$ScenarioId == "NUM-GPCM-ALIGN-RATER", , drop = FALSE
  ]
  expect_identical(sort(criterion$Estimator), c("JML", "MML"))
  expect_identical(sort(rater$Estimator), c("JML", "MML"))
  expect_true(all(criterion$SlopeOwner == "Criterion"))
  expect_true(all(criterion$StepOwner == "Criterion"))
  expect_true(all(rater$SlopeOwner == "Rater"))
  expect_true(all(rater$StepOwner == "Rater"))
  expect_true(all(c(criterion$LatentDimensionCount,
                    rater$LatentDimensionCount) == "1"))
  expect_true(all(c(criterion$SlopeComposition,
                    rater$SlopeComposition) ==
                  "single_owner_relative_gm1"))
  expect_true(all(criterion$EvidenceStatus == "pilot_required"))
  expect_true(all(rater$EvidenceStatus == "not_run"))

  gate <- paste(readLines(gate_path, warn = FALSE, encoding = "UTF-8"),
                collapse = "\n")
  registered <- unique(identity$ScenarioId)
  expect_true(all(vapply(registered, function(id) {
    grepl(paste0("`", id, "`"), gate, fixed = TRUE)
  }, logical(1))))
  expect_match(gate, "NUM-GPCM-BOUND` is retained only", fixed = TRUE)
  expect_match(gate, "RuntimeIdentity", fixed = TRUE)

  execution <- paste(
    readLines(execution_path, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  expect_match(execution, "mfrmr-gpcm-owner-checkpoint-v1", fixed = TRUE)
  expect_match(execution, "mfrmr-gpcm-owner-completion-v1", fixed = TRUE)
  expect_match(execution, "120 manifest rows", fixed = TRUE)
  expect_match(execution, "There is no outcome-adaptive early stopping",
               fixed = TRUE)

  record <- paste(
    readLines(record_path, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  expect_match(record, "completed `0.2.3-draft.66`", fixed = TRUE)
  expect_match(record, "120/120 rows executed", fixed = TRUE)
  expect_match(record, "Draft.65 defect and correction", fixed = TRUE)

  checklist <- utils::read.csv(
    checklist_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  new_items <- c(
    "gpcm_owner_evidence_partition",
    "gpcm_rater_owner_category_support",
    "gpcm_fit_operating_characteristics",
    "gpcm_dff_estimand_specificity",
    "gpcm_generalized_family_scope_guards"
  )
  expect_true(all(new_items %in% checklist$Item))
  owner_row <- checklist[
    checklist$Item == "gpcm_owner_evidence_partition", , drop = FALSE
  ]
  expect_identical(nrow(owner_row), 1L)
  expect_match(owner_row$ScenarioId, "NUM-GPCM-ALIGN-CRITERION",
               fixed = TRUE)
  expect_match(owner_row$ScenarioId, "NUM-GPCM-ALIGN-RATER", fixed = TRUE)
  expect_false(grepl("NUM-GPCM-BOUND", owner_row$ScenarioId, fixed = TRUE))
  expect_match(
    owner_row$PackageSurface,
    "gpcm-owner-current-default-smoke-p1s-record-0.2.3.md",
    fixed = TRUE
  )
  expect_identical(owner_row$EvidenceStatus, "review")
})
