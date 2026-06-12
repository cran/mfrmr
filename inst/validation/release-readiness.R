# mfrmr release-readiness protocol
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

mfrmr_release_readiness_has_value <- function(x) {
  !is.null(x) && length(x) > 0L && !is.na(x[1]) && nzchar(x[1])
}

mfrmr_release_readiness_prompt_steps <- function(target_version = NULL) {
  target_label <- if (!mfrmr_release_readiness_has_value(target_version)) {
    "the target release"
  } else {
    target_version[1]
  }
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
      paste0(
        "Does DESCRIPTION, NEWS, generated help, and the selected check log ",
        "all describe ", target_label, " rather than a development snapshot ",
        "or stale release artifact?"
      ),
      "Do blocker rows for identification, GPCM slope/information kernels, fair-average uncertainty, person fit, and recovery validation have explicit evidence?",
      "Are supported, caveated, blocked, and deferred bounded-GPCM routes visible before users reach unsupported score-side workflows?",
      "Does the package describe FACETS-style output as comparison and handoff support rather than numerical FACETS reproduction?",
      "Can users start from summaries, status tables, and reusable draw-free plot data before reading row-level internals?",
      "Do public-facing docs use review/check/traceability wording and avoid exposing removed helper names as current API?",
      "Does R CMD build/check complete with zero errors and zero warnings, and does CI preserve cross-platform check evidence?",
      "Do cran-comments, NEWS, and validation artifacts tell the same story about release scope, caveats, and deferred work?"
    ),
    Evidence = c(
      "DESCRIPTION Version; first NEWS heading; check-log package version; absence of development labels in current release files",
      "release-evidence-checklist blocker rows; targeted mathematical tests; recovery-validation summary",
      "gpcm_capability_matrix(); README; vignettes; NEWS deferred-work section; post-0.2.1 GPCM roadmap",
      "facets_positioning_guide(); facets_fit_review(); read_facets_fit_table(); output guide",
      "summary methods; plot(..., draw = FALSE); plot_data(); summary-table bundles",
      "README/vignettes/man/cheatsheet terminology scan",
      "mfrmr.Rcheck/00check.log or attached check log; GitHub Actions warning policy and check artifacts",
      "cran-comments.md; NEWS.md; release-evidence map; GPCM roadmap; external parameter-recovery summary and local review helper"
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

mfrmr_release_readiness_description_version <- function(pkg_dir) {
  description <- file.path(pkg_dir, "DESCRIPTION")
  if (!file.exists(description)) {
    return(NA_character_)
  }
  dcf <- read.dcf(description)
  if ("Version" %in% colnames(dcf)) {
    as.character(dcf[1, "Version"])
  } else {
    NA_character_
  }
}

mfrmr_release_readiness_resolve_target_version <- function(pkg_dir,
                                                           target_version = NULL) {
  if (mfrmr_release_readiness_has_value(target_version)) {
    return(target_version[1])
  }
  mfrmr_release_readiness_description_version(pkg_dir)
}

mfrmr_release_readiness_versioned_file <- function(validation_dir,
                                                   prefix,
                                                   target_version,
                                                   fallback_version = "0.2.0",
                                                   ext) {
  candidates <- character(0)
  if (mfrmr_release_readiness_has_value(target_version)) {
    candidates <- c(candidates, file.path(validation_dir, paste0(prefix, target_version, ext)))
  }
  candidates <- c(
    candidates,
    file.path(validation_dir, paste0(prefix, fallback_version, ext))
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits) > 0L) {
    hits[1]
  } else {
    candidates[1]
  }
}

mfrmr_release_readiness_paths <- function(pkg_dir = ".",
                                          target_version = NULL) {
  pkg_dir <- normalizePath(pkg_dir, winslash = "/", mustWork = FALSE)
  if (!file.exists(file.path(pkg_dir, "DESCRIPTION")) &&
      identical(basename(pkg_dir), "inst") &&
      file.exists(file.path(dirname(pkg_dir), "DESCRIPTION"))) {
    pkg_dir <- dirname(pkg_dir)
  }
  target_version <- mfrmr_release_readiness_resolve_target_version(
    pkg_dir,
    target_version = target_version
  )
  validation_dir <- file.path(pkg_dir, "inst", "validation")
  if (!dir.exists(validation_dir)) {
    validation_dir <- file.path(pkg_dir, "validation")
  }
  list(
    target_version = target_version,
    pkg_dir = pkg_dir,
    description = file.path(pkg_dir, "DESCRIPTION"),
    news = file.path(pkg_dir, "NEWS.md"),
    cran_comments = file.path(pkg_dir, "cran-comments.md"),
    ci_workflow = file.path(pkg_dir, ".github", "workflows", "R-CMD-check.yaml"),
    evidence_map = mfrmr_release_readiness_versioned_file(
      validation_dir,
      prefix = "release-evidence-map-",
      target_version = target_version,
      ext = ".md"
    ),
    evidence_checklist = mfrmr_release_readiness_versioned_file(
      validation_dir,
      prefix = "release-evidence-checklist-",
      target_version = target_version,
      ext = ".csv"
    ),
    gpcm_roadmap = file.path(validation_dir, "gpcm-post-0.2.1-roadmap.md"),
    gpcm_capability_source = file.path(pkg_dir, "R", "help_gpcm_scope.R"),
    external_recovery_evidence = mfrmr_release_readiness_versioned_file(
      validation_dir,
      prefix = "external-parameter-recovery-simulation-",
      target_version = target_version,
      ext = ".md"
    ),
    external_recovery_helper = file.path(validation_dir, "external-recovery-audit.R"),
    check_log = file.path(pkg_dir, "mfrmr.Rcheck", "00check.log")
  )
}

mfrmr_release_readiness_check_log_package_version <- function(lines) {
  version_line <- grep("\\* this is package .* version ", lines, value = TRUE)
  if (length(version_line) == 0L) {
    return(NA_character_)
  }
  version_line <- tail(version_line, 1L)
  match <- regexec("version [‘'`](.*)[’'`]", version_line)
  parsed <- regmatches(version_line, match)[[1]]
  if (length(parsed) >= 2L) {
    parsed[2]
  } else {
    NA_character_
  }
}

mfrmr_release_readiness_find_check_log <- function(pkg_dir,
                                                   target_version = NULL) {
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
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0L &&
      mfrmr_release_readiness_has_value(target_version)) {
    versions <- vapply(existing, function(path) {
      mfrmr_release_readiness_check_log_package_version(
        mfrmr_release_readiness_read_lines(path)
      )
    }, character(1))
    matching <- existing[!is.na(versions) & versions == target_version[1]]
    if (length(matching) > 0L) {
      return(matching[1])
    }
  }
  existing[1] %||% file.path(pkg_dir, "mfrmr.Rcheck", "00check.log")
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

mfrmr_release_readiness_parse_check_log <- function(path,
                                                    target_version = NULL) {
  lines <- mfrmr_release_readiness_read_lines(path)
  if (length(lines) == 0L) {
    return(data.frame(
      CheckLog = path,
      PackageVersion = NA_character_,
      TargetVersion = target_version %||% NA_character_,
      VersionMatchesTarget = if (is.null(target_version)) NA else FALSE,
      StatusLine = NA_character_,
      Errors = NA_integer_,
      Warnings = NA_integer_,
      Notes = NA_integer_,
      CheckPassed = FALSE,
      NeedsExplanation = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  package_version <- mfrmr_release_readiness_check_log_package_version(lines)
  version_matches_target <- if (!mfrmr_release_readiness_has_value(target_version)) {
    NA
  } else {
    identical(package_version, target_version[1])
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
    PackageVersion = package_version,
    TargetVersion = target_version %||% NA_character_,
    VersionMatchesTarget = version_matches_target,
    StatusLine = status,
    Errors = errors,
    Warnings = warnings,
    Notes = notes,
    CheckPassed = isTRUE(errors == 0L && warnings == 0L),
    NeedsExplanation = isTRUE(notes > 0L),
    stringsAsFactors = FALSE
  )
}

mfrmr_release_readiness_version_status <- function(paths, target_version = NULL) {
  target_version <- target_version %||% paths$target_version
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

mfrmr_release_readiness_gpcm_scope_status <- function(paths,
                                                      checklist_status = NULL) {
  empty_status <- function(status, detail) {
    data.frame(
      GPCMScopeStatus = status,
      MatrixRows = 0L,
      OutstandingRows = 0L,
      RoadmapRows = 0L,
      ChecklistRoadmapRows = if (!is.null(checklist_status)) {
        checklist_status$RoadmapRows[1] %||% NA_integer_
      } else {
        NA_integer_
      },
      RuntimeGuardRows = 0L,
      RuntimeGuardAreas = 0L,
      GuidanceComplete = FALSE,
      RoadmapCoversOutstanding = FALSE,
      RuntimeGuardCoverageOK = FALSE,
      RuntimeGuardStatusOK = FALSE,
      MissingRoadmapAreas = "",
      MissingRuntimeGuardAreas = "",
      Detail = detail,
      stringsAsFactors = FALSE
    )
  }

  if (!file.exists(paths$gpcm_roadmap)) {
    return(empty_status("concern", "GPCM roadmap is missing"))
  }

  env <- new.env(parent = globalenv())
  if (file.exists(paths$gpcm_capability_source)) {
    source(paths$gpcm_capability_source, local = env)
  } else if (isNamespaceLoaded("mfrmr") ||
             requireNamespace("mfrmr", quietly = TRUE)) {
    # Installed-package review context (for example R CMD check runs or CI
    # artifact reviews): installed packages do not retain `R/` source files,
    # so read the capability matrix and guard coverage from the installed
    # namespace instead of the source tree.
    env$gpcm_capability_matrix <-
      getExportedValue("mfrmr", "gpcm_capability_matrix")
    env$gpcm_runtime_guard_coverage <-
      getExportedValue("mfrmr", "gpcm_runtime_guard_coverage")
  } else {
    return(empty_status("concern", "GPCM capability source is missing"))
  }
  if (!exists("gpcm_capability_matrix", envir = env, inherits = FALSE)) {
    return(empty_status("concern", "GPCM capability matrix function is missing"))
  }
  if (!exists("gpcm_runtime_guard_coverage", envir = env, inherits = FALSE)) {
    return(empty_status("concern", "GPCM runtime guard coverage function is missing"))
  }

  matrix <- env$gpcm_capability_matrix()
  required_columns <- c(
    "Area", "Status", "RecommendedRoute", "NextValidationStep"
  )
  missing_columns <- setdiff(required_columns, names(matrix))
  if (length(missing_columns) > 0L) {
    return(empty_status(
      "concern",
      paste("GPCM capability matrix missing columns:", paste(missing_columns, collapse = ", "))
    ))
  }
  guard_coverage <- env$gpcm_runtime_guard_coverage()
  required_guard_columns <- c(
    "Area", "Helper", "Status", "GuardMode", "ExpectedConditionClass",
    "RecommendedRoute", "NextValidationStep"
  )
  missing_guard_columns <- setdiff(required_guard_columns, names(guard_coverage))
  if (length(missing_guard_columns) > 0L) {
    return(empty_status(
      "concern",
      paste("GPCM runtime guard coverage missing columns:",
            paste(missing_guard_columns, collapse = ", "))
    ))
  }

  outstanding <- matrix[matrix$Status %in% c("blocked", "deferred"), , drop = FALSE]
  guard_idx <- match(guard_coverage$Area, matrix$Area)
  guard_status_ok <- all(!is.na(guard_idx)) &&
    all(guard_coverage$Status == matrix$Status[guard_idx]) &&
    all(guard_coverage$RecommendedRoute == matrix$RecommendedRoute[guard_idx]) &&
    all(guard_coverage$NextValidationStep == matrix$NextValidationStep[guard_idx])
  runtime_rows <- guard_coverage[guard_coverage$GuardMode == "runtime_error", , drop = FALSE]
  runtime_condition_ok <- nrow(runtime_rows) > 0L &&
    all(runtime_rows$ExpectedConditionClass == "mfrmr_gpcm_scope_error")
  covered_guard_areas <- unique(guard_coverage$Area[
    guard_coverage$GuardMode %in% c("runtime_error", "roadmap_only")
  ])
  missing_guard_areas <- setdiff(outstanding$Area, covered_guard_areas)
  runtime_guard_coverage_ok <- length(missing_guard_areas) == 0L &&
    isTRUE(guard_status_ok) &&
    isTRUE(runtime_condition_ok)
  roadmap <- paste(mfrmr_release_readiness_read_lines(paths$gpcm_roadmap), collapse = "\n")
  area_present <- vapply(outstanding$Area, function(area) {
    grepl(area, roadmap, fixed = TRUE)
  }, logical(1))
  missing_areas <- outstanding$Area[!area_present]
  guidance_complete <- all(nzchar(outstanding$RecommendedRoute)) &&
    all(nzchar(outstanding$NextValidationStep))
  checklist_rows <- if (!is.null(checklist_status)) {
    checklist_status$RoadmapRows[1] %||% NA_integer_
  } else {
    NA_integer_
  }
  checklist_covers <- if (is.na(checklist_rows)) {
    NA
  } else {
    checklist_rows >= nrow(outstanding)
  }
  ok <- guidance_complete &&
    length(missing_areas) == 0L &&
    (is.na(checklist_covers) || isTRUE(checklist_covers)) &&
    isTRUE(runtime_guard_coverage_ok)

  data.frame(
    GPCMScopeStatus = if (ok) "ok" else "concern",
    MatrixRows = nrow(matrix),
    OutstandingRows = nrow(outstanding),
    RoadmapRows = length(grep("^### ", mfrmr_release_readiness_read_lines(paths$gpcm_roadmap))),
    ChecklistRoadmapRows = checklist_rows,
    RuntimeGuardRows = nrow(runtime_rows),
    RuntimeGuardAreas = length(unique(guard_coverage$Area)),
    GuidanceComplete = guidance_complete,
    RoadmapCoversOutstanding = length(missing_areas) == 0L,
    RuntimeGuardCoverageOK = runtime_guard_coverage_ok,
    RuntimeGuardStatusOK = isTRUE(guard_status_ok) && isTRUE(runtime_condition_ok),
    MissingRoadmapAreas = paste(missing_areas, collapse = " | "),
    MissingRuntimeGuardAreas = paste(missing_guard_areas, collapse = " | "),
    Detail = paste0(
      "outstanding_rows=", nrow(outstanding),
      "; guidance_complete=", guidance_complete,
      "; roadmap_covers_outstanding=", length(missing_areas) == 0L,
      "; checklist_roadmap_rows=", checklist_rows,
      "; runtime_guard_coverage=", runtime_guard_coverage_ok,
      "; runtime_guard_rows=", nrow(runtime_rows)
    ),
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
                                                 paths,
                                                 gpcm_scope_status = NULL) {
  gpcm_scope_ok <- if (is.null(gpcm_scope_status)) {
    TRUE
  } else {
    identical(gpcm_scope_status$GPCMScopeStatus[1], "ok")
  }
  evidence_available <- file.exists(paths$evidence_map) &&
    file.exists(paths$gpcm_roadmap) &&
    file.exists(paths$external_recovery_evidence) &&
    file.exists(paths$external_recovery_helper) &&
    isTRUE(checklist_status$ChecklistAvailable[1]) &&
    isTRUE(gpcm_scope_ok)
  rows <- data.frame(
    Gate = c("version_contract", "package_check", "ci_workflow", "terminology", "evidence_artifacts"),
    Status = c(
      if (isTRUE(version_status$VersionOK[1])) "ok" else "concern",
      if (isTRUE(check_status$CheckPassed[1]) &&
          !isTRUE(check_status$NeedsExplanation[1]) &&
          !identical(check_status$VersionMatchesTarget[1], FALSE)) {
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
      paste0(
        check_status$StatusLine[1],
        "; check_version=", check_status$PackageVersion[1],
        "; target=", check_status$TargetVersion[1],
        "; version_match=", check_status$VersionMatchesTarget[1]
      ),
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
        "; gpcm_roadmap=", file.exists(paths$gpcm_roadmap),
        "; external_recovery=", file.exists(paths$external_recovery_evidence),
        "; external_helper=", file.exists(paths$external_recovery_helper),
        "; checklist_rows=", checklist_status$Rows[1],
        "; gpcm_scope=", if (is.null(gpcm_scope_status)) {
          "not_checked"
        } else {
          gpcm_scope_status$GPCMScopeStatus[1]
        },
        "; gpcm_runtime_guard=", if (is.null(gpcm_scope_status)) {
          "not_checked"
        } else {
          gpcm_scope_status$RuntimeGuardCoverageOK[1]
        }
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
                                           target_version = NULL,
                                           external_recovery_dir = NULL) {
  paths <- mfrmr_release_readiness_paths(pkg_dir, target_version = target_version)
  target_version <- target_version %||% paths$target_version
  if (!is.null(check_log)) {
    paths$check_log <- check_log
  } else {
    paths$check_log <- mfrmr_release_readiness_find_check_log(
      paths$pkg_dir,
      target_version = target_version
    )
  }
  if (!is.null(checklist)) {
    paths$evidence_checklist <- checklist
  }
  version_status <- mfrmr_release_readiness_version_status(paths, target_version = target_version)
  check_status <- mfrmr_release_readiness_parse_check_log(
    paths$check_log,
    target_version = target_version
  )
  term_status <- mfrmr_release_readiness_term_status(paths$pkg_dir)
  checklist_status <- mfrmr_release_readiness_checklist_status(paths$evidence_checklist)
  ci_workflow_status <- mfrmr_release_readiness_ci_workflow_status(paths$ci_workflow)
  gpcm_scope_status <- mfrmr_release_readiness_gpcm_scope_status(
    paths = paths,
    checklist_status = checklist_status
  )
  gate_summary <- mfrmr_release_readiness_gate_summary(
    version_status = version_status,
    check_status = check_status,
    term_status = term_status,
    checklist_status = checklist_status,
    ci_workflow_status = ci_workflow_status,
    paths = paths,
    gpcm_scope_status = gpcm_scope_status
  )
  external_recovery_status <- mfrmr_release_readiness_external_recovery_status(
    paths = paths,
    external_recovery_dir = external_recovery_dir
  )
  out <- list(
    prompt_steps = mfrmr_release_readiness_prompt_steps(target_version = target_version),
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
    gpcm_scope_status = gpcm_scope_status,
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
    gpcm_scope_status = object$gpcm_scope_status,
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
  if (!is.null(x$gpcm_scope_status)) {
    cat("\nGPCM scope status:\n")
    print(x$gpcm_scope_status, row.names = FALSE)
  }
  if (isTRUE(x$external_recovery_status$ExternalRecoveryRequested[1])) {
    cat("\nExternal recovery status:\n")
    print(x$external_recovery_status, row.names = FALSE)
  }
  invisible(x)
}
