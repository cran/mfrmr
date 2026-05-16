# mfrmr 0.2.0 release-readiness protocol
#
# Source this file in a development or release-check session:
#
#   source(system.file("validation", "release-readiness.R", package = "mfrmr"))
#   readiness <- mfrmr_release_readiness_review(pkg_dir = ".")
#   summary(readiness)
#
# The functions are intentionally not exported. They provide a reproducible
# release-readiness review of release evidence without adding work to routine tests.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

mfrmr_release_readiness_prompt_steps <- function() {
  data.frame(
    Step = seq_len(8L),
    Label = c(
      "Version contract",
      "Mathematical blockers",
      "GPCM boundary",
      "FACETS relationship",
      "UX and plot-data access",
      "Terminology",
      "Package check",
      "Submission handoff"
    ),
    Prompt = c(
      "Does DESCRIPTION, NEWS, and generated help all describe the same 0.2.0 release rather than a development snapshot?",
      "Do blocker rows for identification, GPCM slope/information kernels, fair-average uncertainty, person fit, and recovery validation have explicit evidence?",
      "Are supported, caveated, blocked, and deferred bounded-GPCM routes visible before users reach unsupported score-side workflows?",
      "Does the package describe FACETS-style output as comparison and handoff support rather than numerical FACETS reproduction?",
      "Can users start from summaries, status tables, and reusable draw-free plot data before reading row-level internals?",
      "Do public-facing docs use review/check/traceability wording and avoid exposing removed helper names as current API?",
      "Does R CMD build/check complete with zero errors and zero warnings, and does CI preserve cross-platform check evidence?",
      "Do cran-comments, NEWS, and validation artifacts tell the same story about release scope, caveats, and deferred work?"
    ),
    Evidence = c(
      "DESCRIPTION Version; first NEWS heading; absence of 0.2.0.9000 in current release files",
      "release-evidence-checklist blocker rows; targeted mathematical tests; recovery-validation summary",
      "gpcm_capability_matrix(); README; vignettes; NEWS deferred-work section",
      "facets_positioning_guide(); facets_fit_review(); read_facets_fit_table(); output guide",
      "summary methods; plot(..., draw = FALSE); plot_data(); summary-table bundles",
      "README/vignettes/man/cheatsheet terminology scan",
      "mfrmr.Rcheck/00check.log or attached check log; GitHub Actions warning policy and check artifacts",
      "cran-comments.md; NEWS.md; release-evidence map; external parameter-recovery summary and local review helper"
    ),
    Gate = c(
      "blocker",
      "blocker",
      "blocker",
      "caveat",
      "caveat",
      "caveat",
      "blocker",
      "caveat"
    ),
    stringsAsFactors = FALSE
  )
}

mfrmr_release_readiness_paths <- function(pkg_dir = ".") {
  pkg_dir <- normalizePath(pkg_dir, winslash = "/", mustWork = FALSE)
  if (!file.exists(file.path(pkg_dir, "DESCRIPTION")) &&
      identical(basename(pkg_dir), "inst") &&
      file.exists(file.path(dirname(pkg_dir), "DESCRIPTION"))) {
    pkg_dir <- dirname(pkg_dir)
  }
  validation_dir <- file.path(pkg_dir, "inst", "validation")
  if (!dir.exists(validation_dir)) {
    validation_dir <- file.path(pkg_dir, "validation")
  }
  list(
    pkg_dir = pkg_dir,
    description = file.path(pkg_dir, "DESCRIPTION"),
    news = file.path(pkg_dir, "NEWS.md"),
    cran_comments = file.path(pkg_dir, "cran-comments.md"),
    ci_workflow = file.path(pkg_dir, ".github", "workflows", "R-CMD-check.yaml"),
    evidence_map = file.path(validation_dir, "release-evidence-map-0.2.0.md"),
    evidence_checklist = file.path(validation_dir, "release-evidence-checklist-0.2.0.csv"),
    external_recovery_evidence = file.path(validation_dir, "external-parameter-recovery-simulation-0.2.0.md"),
    external_recovery_helper = file.path(validation_dir, "external-recovery-audit.R"),
    check_log = file.path(pkg_dir, "mfrmr.Rcheck", "00check.log")
  )
}

mfrmr_release_readiness_find_check_log <- function(pkg_dir) {
  candidates <- c(
    file.path(pkg_dir, "mfrmr.Rcheck", "00check.log"),
    file.path(pkg_dir, "check", "mfrmr.Rcheck", "00check.log")
  )
  recursive <- if (dir.exists(pkg_dir)) {
    list.files(pkg_dir, pattern = "^00check[.]log$", recursive = TRUE, full.names = TRUE)
  } else {
    character(0)
  }
  candidates <- unique(c(candidates, recursive))
  candidates[file.exists(candidates)][1] %||% file.path(pkg_dir, "mfrmr.Rcheck", "00check.log")
}

mfrmr_release_readiness_read_lines <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

mfrmr_release_readiness_count_status <- function(status_line, label) {
  if (length(status_line) == 0L || !nzchar(status_line[1])) {
    return(NA_integer_)
  }
  pattern <- paste0("([0-9]+) ", label, "S?")
  hit <- regexec(pattern, status_line[1], ignore.case = TRUE)
  match <- regmatches(status_line[1], hit)[[1]]
  if (length(match) >= 2L) {
    return(as.integer(match[2]))
  }
  0L
}

mfrmr_release_readiness_parse_check_log <- function(path) {
  lines <- mfrmr_release_readiness_read_lines(path)
  if (length(lines) == 0L) {
    return(data.frame(
      CheckLog = path,
      StatusLine = NA_character_,
      Errors = NA_integer_,
      Warnings = NA_integer_,
      Notes = NA_integer_,
      CheckPassed = FALSE,
      NeedsExplanation = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  status <- grep("^Status:", lines, value = TRUE)
  status <- if (length(status) > 0L) tail(status, 1L) else "Status: OK"
  errors <- mfrmr_release_readiness_count_status(status, "ERROR")
  warnings <- mfrmr_release_readiness_count_status(status, "WARNING")
  notes <- mfrmr_release_readiness_count_status(status, "NOTE")
  if (identical(status, "Status: OK")) {
    errors <- warnings <- notes <- 0L
  }
  data.frame(
    CheckLog = path,
    StatusLine = status,
    Errors = errors,
    Warnings = warnings,
    Notes = notes,
    CheckPassed = isTRUE(errors == 0L && warnings == 0L),
    NeedsExplanation = isTRUE(notes > 0L),
    stringsAsFactors = FALSE
  )
}

mfrmr_release_readiness_version_status <- function(paths, target_version = "0.2.0") {
  desc_version <- NA_character_
  if (file.exists(paths$description)) {
    dcf <- read.dcf(paths$description)
    if ("Version" %in% colnames(dcf)) {
      desc_version <- as.character(dcf[1, "Version"])
    }
  }
  news_lines <- mfrmr_release_readiness_read_lines(paths$news)
  first_heading <- news_lines[grep("^# ", news_lines)][1] %||% NA_character_
  current_files <- c(paths$description, paths$news, paths$cran_comments, paths$evidence_map)
  current_lines <- unlist(lapply(current_files[file.exists(current_files)], mfrmr_release_readiness_read_lines), use.names = FALSE)
  dev_label_present <- any(grepl(paste0("\\b", gsub(".", "\\\\.", target_version, fixed = TRUE), "\\.9000\\b"), current_lines))
  data.frame(
    TargetVersion = target_version,
    DescriptionVersion = desc_version,
    NewsHeading = first_heading,
    DevelopmentLabelPresent = dev_label_present,
    VersionOK = identical(desc_version, target_version) &&
      identical(first_heading, paste("# mfrmr", target_version)) &&
      !isTRUE(dev_label_present),
    stringsAsFactors = FALSE
  )
}

mfrmr_release_readiness_public_doc_files <- function(pkg_dir) {
  list_rmd <- function(path) {
    if (!dir.exists(path)) {
      return(character(0))
    }
    list.files(path, pattern = "\\.Rmd$", recursive = TRUE, full.names = TRUE)
  }
  list_rd <- function(path) {
    if (!dir.exists(path)) {
      return(character(0))
    }
    list.files(path, pattern = "\\.Rd$", recursive = TRUE, full.names = TRUE)
  }
  candidates <- c(
    file.path(pkg_dir, "README.md"),
    list_rmd(file.path(pkg_dir, "vignettes")),
    list_rd(file.path(pkg_dir, "man")),
    list_rmd(file.path(pkg_dir, "inst", "cheatsheet"))
  )
  candidates[file.exists(candidates)]
}

mfrmr_release_readiness_term_status <- function(pkg_dir) {
  files <- mfrmr_release_readiness_public_doc_files(pkg_dir)
  allow_source_header <- "^% Please edit documentation in R/.*audit.*\\.R$"
  hits <- character(0)
  for (path in files) {
    lines <- mfrmr_release_readiness_read_lines(path)
    idx <- grep("\\baudit\\b|\\bAudit\\b|_audit|audit_", lines, perl = TRUE)
    if (length(idx) == 0L) {
      next
    }
    rel <- sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", pkg_dir), "/?"), "", path)
    for (i in idx) {
      line <- lines[[i]]
      if (!grepl(allow_source_header, line, perl = TRUE)) {
        hits <- c(hits, paste0(rel, ":", i, ": ", line))
      }
    }
  }
  data.frame(
    FilesScanned = length(files),
    DisallowedRemovedTerms = length(hits),
    TerminologyOK = length(hits) == 0L,
    Examples = paste(utils::head(hits, 5L), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

mfrmr_release_readiness_checklist_status <- function(path) {
  if (!file.exists(path)) {
    return(data.frame(
      Checklist = path,
      Rows = 0L,
      BlockerRows = 0L,
      CaveatRows = 0L,
      RoadmapRows = 0L,
      ChecklistAvailable = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  checklist <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  decision <- as.character(checklist$ReleaseDecision %||% character(0))
  data.frame(
    Checklist = path,
    Rows = nrow(checklist),
    BlockerRows = sum(decision == "blocker_if_failed", na.rm = TRUE),
    CaveatRows = sum(decision == "caveat_if_incomplete", na.rm = TRUE),
    RoadmapRows = sum(decision == "roadmap_if_missing", na.rm = TRUE),
    ChecklistAvailable = nrow(checklist) > 0L,
    stringsAsFactors = FALSE
  )
}

mfrmr_release_readiness_ci_workflow_status <- function(path) {
  if (!file.exists(path)) {
    return(data.frame(
      Workflow = path,
      WorkflowAvailable = FALSE,
      MatrixIncludesMainOS = FALSE,
      MatrixIncludesRDevelOldrelRelease = FALSE,
      PackageCheckStepPresent = FALSE,
      WarningsAreFailures = FALSE,
      CheckArtifactsUploaded = FALSE,
      ReadinessGatePresent = FALSE,
      CIWorkflowOK = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  lines <- mfrmr_release_readiness_read_lines(path)
  contains <- function(pattern) {
    any(grepl(pattern, lines, fixed = TRUE))
  }
  matrix_os <- all(vapply(
    c("ubuntu-latest", "macos-latest", "windows-latest"),
    contains,
    logical(1)
  ))
  matrix_r <- all(vapply(
    c("devel", "release", "oldrel-1"),
    contains,
    logical(1)
  ))
  package_check <- contains("r-lib/actions/check-r-package@v2")
  warning_policy <- any(grepl("error-on:", lines, fixed = TRUE) &
    grepl("warning", lines, fixed = TRUE))
  artifact_upload <- contains("actions/upload-artifact@v4") &&
    any(grepl("check", lines, fixed = TRUE) | grepl("Rcheck", lines, fixed = TRUE))
  readiness_gate <- contains("Release-readiness gate") &&
    contains("mfrmr_release_readiness_review")
  data.frame(
    Workflow = path,
    WorkflowAvailable = TRUE,
    MatrixIncludesMainOS = matrix_os,
    MatrixIncludesRDevelOldrelRelease = matrix_r,
    PackageCheckStepPresent = package_check,
    WarningsAreFailures = warning_policy,
    CheckArtifactsUploaded = artifact_upload,
    ReadinessGatePresent = readiness_gate,
    CIWorkflowOK = isTRUE(matrix_os && matrix_r && package_check && warning_policy &&
      artifact_upload && readiness_gate),
    stringsAsFactors = FALSE
  )
}

mfrmr_release_readiness_gate_summary <- function(version_status,
                                                 check_status,
                                                 term_status,
                                                 checklist_status,
                                                 ci_workflow_status,
                                                 paths) {
  evidence_available <- file.exists(paths$evidence_map) &&
    file.exists(paths$external_recovery_evidence) &&
    file.exists(paths$external_recovery_helper) &&
    isTRUE(checklist_status$ChecklistAvailable[1])
  rows <- data.frame(
    Gate = c("version_contract", "package_check", "ci_workflow", "terminology", "evidence_artifacts"),
    Status = c(
      if (isTRUE(version_status$VersionOK[1])) "ok" else "concern",
      if (isTRUE(check_status$CheckPassed[1]) && !isTRUE(check_status$NeedsExplanation[1])) {
        "ok"
      } else if (isTRUE(check_status$CheckPassed[1])) {
        "review"
      } else {
        "concern"
      },
      if (isTRUE(ci_workflow_status$CIWorkflowOK[1])) "ok" else "review",
      if (isTRUE(term_status$TerminologyOK[1])) "ok" else "concern",
      if (evidence_available) "ok" else "concern"
    ),
    Detail = c(
      paste0("DESCRIPTION=", version_status$DescriptionVersion[1], "; NEWS=", version_status$NewsHeading[1]),
      check_status$StatusLine[1],
      paste0(
        "workflow=", ci_workflow_status$WorkflowAvailable[1],
        "; check_step=", ci_workflow_status$PackageCheckStepPresent[1],
        "; warnings_fail=", ci_workflow_status$WarningsAreFailures[1],
        "; artifacts=", ci_workflow_status$CheckArtifactsUploaded[1],
        "; gate=", ci_workflow_status$ReadinessGatePresent[1]
      ),
      paste0(term_status$DisallowedRemovedTerms[1], " disallowed removed-name hit(s)"),
      paste0(
        "evidence_map=", file.exists(paths$evidence_map),
        "; external_recovery=", file.exists(paths$external_recovery_evidence),
        "; external_helper=", file.exists(paths$external_recovery_helper),
        "; checklist_rows=", checklist_status$Rows[1]
      )
    ),
    stringsAsFactors = FALSE
  )
  rows
}

mfrmr_release_readiness_decision <- function(gate_summary) {
  status <- as.character(gate_summary$Status)
  if (any(status == "concern", na.rm = TRUE)) {
    "concern"
  } else if (any(status == "review", na.rm = TRUE)) {
    "review"
  } else {
    "ok"
  }
}

mfrmr_release_readiness_external_recovery_status <- function(paths, external_recovery_dir = NULL) {
  if (is.null(external_recovery_dir) || !nzchar(external_recovery_dir)) {
    return(data.frame(
      ExternalRecoveryRequested = FALSE,
      ExternalRecoveryDir = NA_character_,
      EvidenceStatus = NA_character_,
      RequiredSchemaOK = NA,
      Detail = "not requested",
      stringsAsFactors = FALSE
    ))
  }
  external_recovery_dir <- normalizePath(external_recovery_dir, winslash = "/", mustWork = FALSE)
  if (!file.exists(paths$external_recovery_helper)) {
    return(data.frame(
      ExternalRecoveryRequested = TRUE,
      ExternalRecoveryDir = external_recovery_dir,
      EvidenceStatus = "concern",
      RequiredSchemaOK = FALSE,
      Detail = "external recovery helper is missing",
      stringsAsFactors = FALSE
    ))
  }
  env <- new.env(parent = globalenv())
  source(paths$external_recovery_helper, local = env)
  if (!exists("mfrmr_review_external_recovery_simulation", envir = env, inherits = FALSE)) {
    return(data.frame(
      ExternalRecoveryRequested = TRUE,
      ExternalRecoveryDir = external_recovery_dir,
      EvidenceStatus = "concern",
      RequiredSchemaOK = FALSE,
      Detail = "external recovery helper did not define the review function",
      stringsAsFactors = FALSE
    ))
  }
  review <- env$mfrmr_review_external_recovery_simulation(external_recovery_dir)
  required_schema <- review$schema_status$Required
  required_schema_ok <- if (length(required_schema) == 0L) {
    FALSE
  } else {
    all(review$schema_status$SchemaOK[required_schema], na.rm = FALSE)
  }
  data.frame(
    ExternalRecoveryRequested = TRUE,
    ExternalRecoveryDir = external_recovery_dir,
    EvidenceStatus = as.character(review$decision$EvidenceStatus[1]),
    RequiredSchemaOK = isTRUE(required_schema_ok),
    Detail = as.character(review$decision$Interpretation[1]),
    stringsAsFactors = FALSE
  )
}

mfrmr_release_readiness_review <- function(pkg_dir = ".",
                                           check_log = NULL,
                                           checklist = NULL,
                                           target_version = "0.2.0",
                                           external_recovery_dir = NULL) {
  paths <- mfrmr_release_readiness_paths(pkg_dir)
  if (!is.null(check_log)) {
    paths$check_log <- check_log
  } else {
    paths$check_log <- mfrmr_release_readiness_find_check_log(paths$pkg_dir)
  }
  if (!is.null(checklist)) {
    paths$evidence_checklist <- checklist
  }
  version_status <- mfrmr_release_readiness_version_status(paths, target_version = target_version)
  check_status <- mfrmr_release_readiness_parse_check_log(paths$check_log)
  term_status <- mfrmr_release_readiness_term_status(paths$pkg_dir)
  checklist_status <- mfrmr_release_readiness_checklist_status(paths$evidence_checklist)
  ci_workflow_status <- mfrmr_release_readiness_ci_workflow_status(paths$ci_workflow)
  gate_summary <- mfrmr_release_readiness_gate_summary(
    version_status = version_status,
    check_status = check_status,
    term_status = term_status,
    checklist_status = checklist_status,
    ci_workflow_status = ci_workflow_status,
    paths = paths
  )
  external_recovery_status <- mfrmr_release_readiness_external_recovery_status(
    paths = paths,
    external_recovery_dir = external_recovery_dir
  )
  out <- list(
    prompt_steps = mfrmr_release_readiness_prompt_steps(),
    gate_summary = gate_summary,
    release_decision = data.frame(
      ReleaseReadinessStatus = mfrmr_release_readiness_decision(gate_summary),
      Explanation = paste(gate_summary$Gate, gate_summary$Status, sep = "=", collapse = "; "),
      stringsAsFactors = FALSE
    ),
    version_status = version_status,
    check_status = check_status,
    ci_workflow_status = ci_workflow_status,
    terminology_status = term_status,
    checklist_status = checklist_status,
    external_recovery_status = external_recovery_status,
    paths = paths
  )
  class(out) <- "mfrmr_release_readiness_review"
  out
}

summary.mfrmr_release_readiness_review <- function(object, ...) {
  if (!inherits(object, "mfrmr_release_readiness_review")) {
    stop("`object` must be output from mfrmr_release_readiness_review().", call. = FALSE)
  }
  out <- list(
    release_decision = object$release_decision,
    gate_summary = object$gate_summary,
    prompt_steps = object$prompt_steps,
    check_status = object$check_status,
    ci_workflow_status = object$ci_workflow_status,
    external_recovery_status = object$external_recovery_status
  )
  class(out) <- "summary.mfrmr_release_readiness_review"
  out
}

print.mfrmr_release_readiness_review <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

print.summary.mfrmr_release_readiness_review <- function(x, ...) {
  cat("mfrmr release-readiness review\n\n")
  print(x$release_decision, row.names = FALSE)
  cat("\nGate summary:\n")
  print(x$gate_summary, row.names = FALSE)
  cat("\nReview steps:\n")
  print(x$prompt_steps[, c("Step", "Label", "Gate")], row.names = FALSE)
  if (isTRUE(x$external_recovery_status$ExternalRecoveryRequested[1])) {
    cat("\nExternal recovery status:\n")
    print(x$external_recovery_status, row.names = FALSE)
  }
  invisible(x)
}
