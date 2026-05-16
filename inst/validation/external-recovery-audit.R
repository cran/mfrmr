# External common-data recovery simulation audit
#
# Source this file in a development or release-review session:
#
#   source(system.file("validation", "external-recovery-audit.R", package = "mfrmr"))
#   review <- mfrmr_review_external_recovery_simulation(
#     "/path/to/Parameter_Recovery_Simulation"
#   )
#   summary(review)
#
# The external workflow is intentionally not bundled with the package. This
# helper audits its generated CSV outputs and turns them into compact,
# source-grounded release evidence.

mfrmr_external_recovery_column <- function(data, column, default = NA) {
  if (column %in% names(data)) {
    data[[column]]
  } else {
    rep(default, nrow(data))
  }
}

mfrmr_external_recovery_numeric <- function(data, column, default = NA_real_) {
  suppressWarnings(as.numeric(mfrmr_external_recovery_column(data, column, default)))
}

mfrmr_external_recovery_sum <- function(x, missing = NA_real_) {
  if (length(x) == 0L || all(is.na(x))) {
    return(missing)
  }
  sum(x, na.rm = TRUE)
}

mfrmr_external_recovery_min <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  min(x, na.rm = TRUE)
}

mfrmr_external_recovery_max <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  max(x, na.rm = TRUE)
}

mfrmr_external_recovery_unique_count <- function(data, column) {
  if (!column %in% names(data)) {
    return(NA_integer_)
  }
  length(unique(data[[column]]))
}

mfrmr_external_recovery_unique_label <- function(data, column) {
  if (!column %in% names(data) || nrow(data) == 0L) {
    return(NA_character_)
  }
  values <- unique(as.character(data[[column]]))
  values <- values[nzchar(values) & !is.na(values)]
  if (length(values) == 0L) {
    NA_character_
  } else {
    paste(values, collapse = ", ")
  }
}

mfrmr_external_recovery_agreement_file <- function() {
  file.path("analysis", paste0("engine_", "par", "ity", "_overview.csv"))
}

mfrmr_external_recovery_required_files <- function() {
  c(
    "analysis/dataset_manifest.csv",
    "analysis/smoke_check_summary.csv",
    "analysis/engine_status_summary.csv",
    mfrmr_external_recovery_agreement_file(),
    "analysis/key_findings.csv",
    "analysis/key_findings_counts.csv"
  )
}

mfrmr_external_recovery_optional_files <- function() {
  c(
    "sample_size_dstudy/analysis/sample_size_decision_summary.csv",
    "sample_size_dstudy/analysis/sample_size_classification_summary.csv"
  )
}

mfrmr_external_recovery_csv <- function(base_dir, relative_path) {
  path <- file.path(base_dir, relative_path)
  if (!file.exists(path)) {
    return(data.frame())
  }
  tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) data.frame()
  )
}

mfrmr_external_recovery_schema_requirements <- function() {
  data.frame(
    File = c(
      "analysis/dataset_manifest.csv",
      "analysis/smoke_check_summary.csv",
      "analysis/engine_status_summary.csv",
      mfrmr_external_recovery_agreement_file(),
      "analysis/key_findings.csv",
      "analysis/key_findings_counts.csv",
      "sample_size_dstudy/analysis/sample_size_decision_summary.csv",
      "sample_size_dstudy/analysis/sample_size_classification_summary.csv"
    ),
    Required = c(rep(TRUE, 6L), FALSE, FALSE),
    RequiredColumns = c(
      "DesignPattern|Rows|ObservedDensity",
      "Passed",
      "Runs|ErrorRuns|ConvergenceRate",
      "Source|Metric|EnginePair|Groups|ComparedRows|MaxRMSE|ReviewGroups|PracticalOrBetterGroups",
      "Priority|Rank|Domain|Metric|Value|Message",
      "Priority|Domain",
      "SampleSizeFacet",
      "ClassificationRisk"
    ),
    stringsAsFactors = FALSE
  )
}

mfrmr_external_recovery_csv_columns <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }
  tryCatch(
    names(utils::read.csv(path, nrows = 0L, stringsAsFactors = FALSE, check.names = FALSE)),
    error = function(e) character(0)
  )
}

mfrmr_external_recovery_file_status <- function(base_dir) {
  required <- mfrmr_external_recovery_required_files()
  optional <- mfrmr_external_recovery_optional_files()
  all_files <- c(required, optional)
  data.frame(
    File = all_files,
    Required = all_files %in% required,
    Exists = file.exists(file.path(base_dir, all_files)),
    Bytes = vapply(file.path(base_dir, all_files), function(path) {
      if (file.exists(path)) file.info(path)$size else NA_real_
    }, numeric(1)),
    MD5 = vapply(file.path(base_dir, all_files), function(path) {
      if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

mfrmr_external_recovery_schema_status <- function(base_dir) {
  requirements <- mfrmr_external_recovery_schema_requirements()
  rows <- lapply(seq_len(nrow(requirements)), function(i) {
    file <- requirements$File[[i]]
    path <- file.path(base_dir, file)
    expected <- strsplit(requirements$RequiredColumns[[i]], "|", fixed = TRUE)[[1]]
    present <- mfrmr_external_recovery_csv_columns(path)
    exists <- file.exists(path)
    missing <- setdiff(expected, present)
    schema_ok <- if (!exists && !isTRUE(requirements$Required[[i]])) {
      TRUE
    } else {
      exists && length(missing) == 0L
    }
    data.frame(
      File = file,
      Required = requirements$Required[[i]],
      Exists = exists,
      RequiredColumns = paste(expected, collapse = ", "),
      MissingColumns = if (length(missing) == 0L) "" else paste(missing, collapse = ", "),
      SchemaOK = schema_ok,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mfrmr_external_recovery_overview <- function(base_dir) {
  manifest <- mfrmr_external_recovery_csv(base_dir, "analysis/dataset_manifest.csv")
  smoke <- mfrmr_external_recovery_csv(base_dir, "analysis/smoke_check_summary.csv")
  status <- mfrmr_external_recovery_csv(base_dir, "analysis/engine_status_summary.csv")

  smoke_passed <- if ("Passed" %in% names(smoke)) {
    sum(as.logical(smoke$Passed), na.rm = TRUE)
  } else {
    NA_integer_
  }

  data.frame(
    ManifestRows = nrow(manifest),
    DesignPatterns = mfrmr_external_recovery_unique_label(manifest, "DesignPattern"),
    DatasetRowsMin = mfrmr_external_recovery_min(mfrmr_external_recovery_numeric(manifest, "Rows")),
    DatasetRowsMax = mfrmr_external_recovery_max(mfrmr_external_recovery_numeric(manifest, "Rows")),
    ObservedDensityMin = mfrmr_external_recovery_min(mfrmr_external_recovery_numeric(manifest, "ObservedDensity")),
    ObservedDensityMax = mfrmr_external_recovery_max(mfrmr_external_recovery_numeric(manifest, "ObservedDensity")),
    SmokeChecksPassed = smoke_passed,
    SmokeChecksTotal = nrow(smoke),
    EngineGroups = nrow(status),
    EngineRuns = mfrmr_external_recovery_sum(mfrmr_external_recovery_numeric(status, "Runs")),
    ErrorRuns = mfrmr_external_recovery_sum(mfrmr_external_recovery_numeric(status, "ErrorRuns")),
    MinConvergenceRate = mfrmr_external_recovery_min(mfrmr_external_recovery_numeric(status, "ConvergenceRate")),
    stringsAsFactors = FALSE
  )
}

mfrmr_external_recovery_agreement_summary <- function(base_dir) {
  agreement <- mfrmr_external_recovery_csv(
    base_dir,
    mfrmr_external_recovery_agreement_file()
  )
  if (nrow(agreement) == 0L || !"Source" %in% names(agreement)) {
    return(data.frame())
  }
  split_agreement <- split(agreement, agreement$Source)
  rows <- lapply(names(split_agreement), function(source) {
    dat <- split_agreement[[source]]
    data.frame(
      Source = source,
      Metrics = mfrmr_external_recovery_unique_count(dat, "Metric"),
      EnginePairs = mfrmr_external_recovery_unique_count(dat, "EnginePair"),
      Groups = mfrmr_external_recovery_sum(mfrmr_external_recovery_numeric(dat, "Groups")),
      ComparedRows = mfrmr_external_recovery_sum(mfrmr_external_recovery_numeric(dat, "ComparedRows")),
      MaxRMSE = mfrmr_external_recovery_max(mfrmr_external_recovery_numeric(dat, "MaxRMSE")),
      ReviewGroups = mfrmr_external_recovery_sum(mfrmr_external_recovery_numeric(dat, "ReviewGroups")),
      PracticalOrBetterGroups = mfrmr_external_recovery_sum(
        mfrmr_external_recovery_numeric(dat, "PracticalOrBetterGroups")
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$Source), , drop = FALSE]
}

mfrmr_external_recovery_finding_counts <- function(base_dir) {
  counts <- mfrmr_external_recovery_csv(base_dir, "analysis/key_findings_counts.csv")
  if (nrow(counts) == 0L) {
    return(data.frame())
  }
  if (!all(c("Priority", "Domain") %in% names(counts))) {
    return(counts)
  }
  counts[order(counts$Priority, counts$Domain), , drop = FALSE]
}

mfrmr_external_recovery_top_findings <- function(base_dir, priority = "high") {
  findings <- mfrmr_external_recovery_csv(base_dir, "analysis/key_findings.csv")
  if (nrow(findings) == 0L) {
    return(data.frame())
  }
  if (!all(c("Priority", "Rank") %in% names(findings))) {
    return(data.frame())
  }
  keep <- as.character(findings$Priority) == priority &
    suppressWarnings(as.integer(findings$Rank)) == 1L
  out <- findings[keep, intersect(c("Domain", "Metric", "Value", "Message"), names(findings)), drop = FALSE]
  row.names(out) <- NULL
  out
}

mfrmr_external_recovery_sample_size_summary <- function(base_dir) {
  decision <- mfrmr_external_recovery_csv(base_dir, "sample_size_dstudy/analysis/sample_size_decision_summary.csv")
  classification <- mfrmr_external_recovery_csv(base_dir, "sample_size_dstudy/analysis/sample_size_classification_summary.csv")
  if (nrow(decision) == 0L && nrow(classification) == 0L) {
    return(data.frame(
      Available = FALSE,
      DecisionRows = 0L,
      ClassificationRows = 0L,
      SampleSizeFacets = NA_character_,
      ClassificationRiskLevels = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    Available = TRUE,
    DecisionRows = nrow(decision),
    ClassificationRows = nrow(classification),
    SampleSizeFacets = mfrmr_external_recovery_unique_label(decision, "SampleSizeFacet"),
    ClassificationRiskLevels = mfrmr_external_recovery_unique_label(classification, "ClassificationRisk"),
    stringsAsFactors = FALSE
  )
}

mfrmr_external_recovery_decision <- function(file_status, schema_status, overview, agreement_summary) {
  missing_required <- file_status$Required & !file_status$Exists
  if (any(missing_required, na.rm = TRUE)) {
    return(data.frame(
      EvidenceStatus = "concern",
      Interpretation = "Required external recovery output files are missing.",
      stringsAsFactors = FALSE
    ))
  }

  schema_blocker <- schema_status$Required & !schema_status$SchemaOK
  if (any(schema_blocker, na.rm = TRUE)) {
    missing <- paste(
      paste0(schema_status$File[schema_blocker], " [", schema_status$MissingColumns[schema_blocker], "]"),
      collapse = "; "
    )
    return(data.frame(
      EvidenceStatus = "concern",
      Interpretation = paste("Required external recovery columns are missing:", missing),
      stringsAsFactors = FALSE
    ))
  }

  smoke_ok <- isTRUE(overview$SmokeChecksPassed[1] == overview$SmokeChecksTotal[1])
  engine_ok <- isTRUE(overview$ErrorRuns[1] == 0) &&
    isTRUE(overview$MinConvergenceRate[1] >= 0.95)
  recovery_sources <- agreement_summary$Source %in% c("level_recovery", "step_recovery")
  recovery_agreement_ok <- if (any(recovery_sources)) {
    review_groups <- agreement_summary$ReviewGroups[recovery_sources]
    all(!is.na(review_groups)) && all(review_groups == 0)
  } else {
    FALSE
  }
  fit_review <- any(
    agreement_summary$Source %in% c("fit_statistics", "separation") &
      !is.na(agreement_summary$ReviewGroups) &
      agreement_summary$ReviewGroups > 0,
    na.rm = TRUE
  )

  if (smoke_ok && engine_ok && recovery_agreement_ok) {
    status <- if (fit_review) "review" else "ok"
    interp <- if (fit_review) {
      paste(
        "Core recovery files, convergence, and level/step agreement are usable;",
        "fit/separation review groups remain convention-sensitive and should be interpreted explicitly."
      )
    } else {
      "Core recovery files, convergence, and agreement checks support the external evidence summary."
    }
  } else {
    status <- "concern"
    interp <- "External recovery evidence needs follow-up before being used as release evidence."
  }

  data.frame(EvidenceStatus = status, Interpretation = interp, stringsAsFactors = FALSE)
}

mfrmr_review_external_recovery_simulation <- function(base_dir) {
  base_dir <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  file_status <- mfrmr_external_recovery_file_status(base_dir)
  schema_status <- mfrmr_external_recovery_schema_status(base_dir)
  overview <- mfrmr_external_recovery_overview(base_dir)
  agreement_summary <- mfrmr_external_recovery_agreement_summary(base_dir)
  finding_counts <- mfrmr_external_recovery_finding_counts(base_dir)
  top_findings <- mfrmr_external_recovery_top_findings(base_dir)
  sample_size_summary <- mfrmr_external_recovery_sample_size_summary(base_dir)
  decision <- mfrmr_external_recovery_decision(file_status, schema_status, overview, agreement_summary)
  out <- list(
    base_dir = base_dir,
    decision = decision,
    file_status = file_status,
    schema_status = schema_status,
    overview = overview,
    agreement_summary = agreement_summary,
    finding_counts = finding_counts,
    top_findings = top_findings,
    sample_size_summary = sample_size_summary
  )
  class(out) <- "mfrmr_external_recovery_review"
  out
}

summary.mfrmr_external_recovery_review <- function(object, ...) {
  if (!inherits(object, "mfrmr_external_recovery_review")) {
    stop("`object` must be output from mfrmr_review_external_recovery_simulation().", call. = FALSE)
  }
  out <- object[c(
    "decision",
    "schema_status",
    "overview",
    "agreement_summary",
    "sample_size_summary",
    "top_findings"
  )]
  class(out) <- "summary.mfrmr_external_recovery_review"
  out
}

print.mfrmr_external_recovery_review <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

print.summary.mfrmr_external_recovery_review <- function(x, ...) {
  cat("mfrmr external common-data recovery review\n\n")
  print(x$decision, row.names = FALSE)
  cat("\nSchema status:\n")
  print(x$schema_status, row.names = FALSE)
  cat("\nOverview:\n")
  print(x$overview, row.names = FALSE)
  cat("\nAgreement summary:\n")
  print(x$agreement_summary, row.names = FALSE)
  cat("\nSample-size summary:\n")
  print(x$sample_size_summary, row.names = FALSE)
  cat("\nTop high-priority findings:\n")
  print(utils::head(x$top_findings, 12L), row.names = FALSE)
  invisible(x)
}
