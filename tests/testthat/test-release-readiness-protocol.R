release_readiness_protocol_path <- function() {
  source_path <- testthat::test_path("..", "..", "inst", "validation", "release-readiness.R")
  if (file.exists(source_path)) {
    source_path
  } else {
    system.file("validation", "release-readiness.R", package = "mfrmr")
  }
}

test_that("release-readiness protocol exposes review steps and parses check logs", {
  protocol <- release_readiness_protocol_path()
  expect_true(nzchar(protocol))

  env <- new.env(parent = globalenv())
  source(protocol, local = env)

  steps <- env$mfrmr_release_readiness_prompt_steps()
  expect_s3_class(steps, "data.frame")
  expect_equal(nrow(steps), 8L)
  expect_true(all(c("Step", "Label", "Prompt", "Evidence", "Gate") %in% names(steps)))
  expect_true(all(c("blocker", "caveat") %in% steps$Gate))

  log_file <- tempfile(fileext = ".log")
  writeLines(c(
    "* this is package ‘mfrmr’ version ‘0.2.1’",
    "* checking package namespace information ... OK",
    "* checking tests ... OK",
    "Status: 1 NOTE"
  ), log_file)
  parsed <- env$mfrmr_release_readiness_parse_check_log(log_file, target_version = "0.2.1")
  expect_identical(parsed$PackageVersion, "0.2.1")
  expect_true(parsed$VersionMatchesTarget)
  expect_true(parsed$CheckPassed)
  expect_true(parsed$NeedsExplanation)
  expect_equal(parsed$Errors, 0L)
  expect_equal(parsed$Warnings, 0L)
  expect_equal(parsed$Notes, 1L)
})

test_that("release-readiness protocol rejects stale check logs by version", {
  protocol <- release_readiness_protocol_path()
  env <- new.env(parent = globalenv())
  source(protocol, local = env)

  log_file <- tempfile(fileext = ".log")
  writeLines(c(
    "* this is package ‘mfrmr’ version ‘0.2.0’",
    "* checking tests ... OK",
    "Status: OK"
  ), log_file)
  parsed <- env$mfrmr_release_readiness_parse_check_log(log_file, target_version = "0.2.1")
  expect_true(parsed$CheckPassed)
  expect_false(parsed$VersionMatchesTarget)

  version_status <- data.frame(
    TargetVersion = "0.2.1",
    DescriptionVersion = "0.2.1",
    NewsHeading = "# mfrmr 0.2.1",
    DevelopmentLabelPresent = FALSE,
    VersionOK = TRUE
  )
  term_status <- data.frame(
    FilesScanned = 1L,
    DisallowedRemovedTerms = 0L,
    TerminologyOK = TRUE,
    Examples = ""
  )
  checklist_status <- data.frame(
    Checklist = tempfile(),
    Rows = 1L,
    BlockerRows = 1L,
    CaveatRows = 0L,
    RoadmapRows = 0L,
    ChecklistAvailable = TRUE
  )
  ci_status <- data.frame(
    WorkflowAvailable = TRUE,
    PackageCheckStepPresent = TRUE,
    WarningsAreFailures = TRUE,
    CheckArtifactsUploaded = TRUE,
    ReadinessGatePresent = TRUE,
    CIWorkflowOK = TRUE
  )
  paths <- list(
    evidence_map = log_file,
    gpcm_roadmap = log_file,
    external_recovery_evidence = log_file,
    external_recovery_helper = log_file
  )
  gate <- env$mfrmr_release_readiness_gate_summary(
    version_status = version_status,
    check_status = parsed,
    term_status = term_status,
    checklist_status = checklist_status,
    ci_workflow_status = ci_status,
    paths = paths
  )
  expect_identical(gate$Status[gate$Gate == "package_check"], "review")
})

test_that("release-readiness protocol finds common check-log locations", {
  protocol <- release_readiness_protocol_path()
  env <- new.env(parent = globalenv())
  source(protocol, local = env)

  root <- tempfile("pkg")
  dir.create(file.path(root, "check", "mfrmr.Rcheck"), recursive = TRUE)
  log_file <- file.path(root, "check", "mfrmr.Rcheck", "00check.log")
  writeLines("Status: OK", log_file)

  found <- env$mfrmr_release_readiness_find_check_log(root)
  expect_identical(normalizePath(found, winslash = "/", mustWork = TRUE),
                   normalizePath(log_file, winslash = "/", mustWork = TRUE))

  stale_root_log <- file.path(root, "mfrmr.Rcheck", "00check.log")
  dir.create(dirname(stale_root_log), recursive = TRUE)
  writeLines(c(
    "* this is package ‘mfrmr’ version ‘0.2.0’",
    "Status: OK"
  ), stale_root_log)
  writeLines(c(
    "* this is package ‘mfrmr’ version ‘0.2.1’",
    "Status: OK"
  ), log_file)
  found_current <- env$mfrmr_release_readiness_find_check_log(
    root,
    target_version = "0.2.1"
  )
  expect_identical(normalizePath(found_current, winslash = "/", mustWork = TRUE),
                   normalizePath(log_file, winslash = "/", mustWork = TRUE))
})

test_that("release-readiness protocol checks CI workflow contract", {
  protocol <- release_readiness_protocol_path()
  env <- new.env(parent = globalenv())
  source(protocol, local = env)

  root <- tempfile("pkg")
  dir.create(file.path(root, ".github", "workflows"), recursive = TRUE)
  workflow <- file.path(root, ".github", "workflows", "R-CMD-check.yaml")
  writeLines(c(
    "name: R-CMD-check",
    "matrix:",
    "  config:",
    "    - {os: macos-latest, r: 'release'}",
    "    - {os: windows-latest, r: 'release'}",
    "    - {os: ubuntu-latest, r: 'devel'}",
    "    - {os: ubuntu-latest, r: 'oldrel-1'}",
    "- uses: r-lib/actions/check-r-package@v2",
    "  with:",
    "    error-on: '\"warning\"'",
    "- name: Upload check results",
    "  uses: actions/upload-artifact@v4",
    "  with:",
    "    path: check",
    "- name: Release-readiness gate",
    "  run: mfrmr_release_readiness_review(pkg_dir = \".\")"
  ), workflow)

  status <- env$mfrmr_release_readiness_ci_workflow_status(workflow)
  expect_true(status$WorkflowAvailable)
  expect_true(status$MatrixIncludesMainOS)
  expect_true(status$MatrixIncludesRDevelOldrelRelease)
  expect_true(status$PackageCheckStepPresent)
  expect_true(status$WarningsAreFailures)
  expect_true(status$CheckArtifactsUploaded)
  expect_true(status$ReadinessGatePresent)
  expect_true(status$CIWorkflowOK)
})

test_that("release-readiness protocol checks GPCM scope alignment", {
  protocol <- release_readiness_protocol_path()
  env <- new.env(parent = globalenv())
  source(protocol, local = env)

  pkg_root <- normalizePath(test_path("..", ".."), winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    pkg_root <- system.file(package = "mfrmr")
  }
  paths <- env$mfrmr_release_readiness_paths(pkg_root, target_version = "0.2.1")
  checklist_status <- env$mfrmr_release_readiness_checklist_status(paths$evidence_checklist)
  status <- env$mfrmr_release_readiness_gpcm_scope_status(
    paths = paths,
    checklist_status = checklist_status
  )

  expect_s3_class(status, "data.frame")
  expect_equal(status$GPCMScopeStatus[1], "ok")
  expect_gt(status$OutstandingRows[1], 0L)
  expect_true(status$GuidanceComplete[1])
  expect_true(status$RoadmapCoversOutstanding[1])
  expect_true(status$RuntimeGuardCoverageOK[1])
  expect_true(status$RuntimeGuardStatusOK[1])
  expect_gt(status$RuntimeGuardRows[1], 0L)
  expect_true(status$RuntimeGuardAreas[1] >= status$OutstandingRows[1])
  expect_identical(status$MissingRuntimeGuardAreas[1], "")
  expect_true(status$ChecklistRoadmapRows[1] >= status$OutstandingRows[1])
})

test_that("release-readiness protocol reviews the source tree shape", {
  protocol <- release_readiness_protocol_path()
  env <- new.env(parent = globalenv())
  source(protocol, local = env)

  pkg_root <- normalizePath(test_path("..", ".."), winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    pkg_root <- system.file(package = "mfrmr")
  }
  expect_true(nzchar(pkg_root))
  review <- env$mfrmr_release_readiness_review(pkg_dir = pkg_root)
  expect_s3_class(review, "mfrmr_release_readiness_review")
  expect_true(all(c(
    "prompt_steps", "gate_summary", "release_decision",
    "version_status", "check_status", "ci_workflow_status", "terminology_status",
    "checklist_status", "gpcm_scope_status", "external_recovery_status"
  ) %in% names(review)))
  expect_false(review$external_recovery_status$ExternalRecoveryRequested[1])
  expect_true(isTRUE(review$version_status$VersionOK[1]))
  expect_true(file.exists(review$paths$gpcm_roadmap))
  expect_equal(review$gpcm_scope_status$GPCMScopeStatus[1], "ok")
  if (file.exists(file.path(pkg_root, ".github", "workflows", "R-CMD-check.yaml"))) {
    expect_true(isTRUE(review$ci_workflow_status$CIWorkflowOK[1]))
  }
  expect_true(isTRUE(review$terminology_status$TerminologyOK[1]))
  expect_true(isTRUE(review$checklist_status$ChecklistAvailable[1]))
})
