
bundle_settings_table <- function(settings) {
  if (is.null(settings) || !is.list(settings) || length(settings) == 0) return(data.frame())
  keys <- names(settings)
  if (is.null(keys) || any(!nzchar(keys))) {
    keys <- paste0("Setting", seq_along(settings))
  }
  vals <- vapply(settings, function(v) {
    if (is.null(v)) return("NULL")
    if (is.data.frame(v)) return(paste0("<table ", nrow(v), "x", ncol(v), ">"))
    if (is.list(v)) return(paste0("<list ", length(v), ">"))
    paste(as.character(v), collapse = ", ")
  }, character(1))
  data.frame(Setting = keys, Value = vals, stringsAsFactors = FALSE)
}

bundle_preview_table <- function(object, top_n = 10L) {
  keys <- c(
    "table", "pairs", "stacked", "ranked_table", "facet_profile", "graphfile",
    "category_table", "facet_coverage", "listing", "overall_table", "by_facet_table",
    "missing_preview", "column_audit", "metric_audit", "column_summary", "metric_summary"
  )
  nm <- names(object)
  if (is.null(nm) || length(nm) == 0) {
    return(list(name = NA_character_, table = data.frame()))
  }
  key <- keys[keys %in% nm][1]
  if (is.na(key) || length(key) == 0) {
    return(list(name = NA_character_, table = data.frame()))
  }
  tbl <- object[[key]]
  if (!is.data.frame(tbl) || nrow(tbl) == 0) {
    return(list(name = key, table = data.frame()))
  }
  top_n <- max(1L, as.integer(top_n))
  list(name = key, table = utils::head(as.data.frame(tbl, stringsAsFactors = FALSE), n = top_n))
}

summarize_bias_count_bundle <- function(object, digits = 3, top_n = 10) {
  tbl <- as.data.frame(object$table %||% data.frame(), stringsAsFactors = FALSE)
  if ("Observd Count" %in% names(tbl) && !"Count" %in% names(tbl)) {
    tbl$Count <- suppressWarnings(as.numeric(tbl$`Observd Count`))
  }
  if (!"LowCountFlag" %in% names(tbl)) {
    tbl$LowCountFlag <- FALSE
  }
  if (!is.logical(tbl$LowCountFlag)) {
    tbl$LowCountFlag <- as.logical(tbl$LowCountFlag)
  }
  if (!"Count" %in% names(tbl)) {
    tbl$Count <- suppressWarnings(as.numeric(tbl$Count))
  }

  cnt <- suppressWarnings(as.numeric(tbl$Count))
  cnt <- cnt[is.finite(cnt)]
  count_distribution <- if (length(cnt) == 0) {
    data.frame()
  } else {
    data.frame(
      Min = min(cnt),
      Q1 = stats::quantile(cnt, 0.25, names = FALSE),
      Median = stats::median(cnt),
      Mean = mean(cnt),
      Q3 = stats::quantile(cnt, 0.75, names = FALSE),
      Max = max(cnt),
      stringsAsFactors = FALSE
    )
  }

  low_tbl <- tbl[tbl$LowCountFlag %in% TRUE, , drop = FALSE]
  if (nrow(low_tbl) > 0 && "Count" %in% names(low_tbl)) {
    low_tbl <- low_tbl |>
      dplyr::arrange(.data$Count) |>
      dplyr::slice_head(n = top_n)
  }

  summary_tbl <- as.data.frame(object$summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(summary_tbl) == 0) {
    summary_tbl <- data.frame(
      Branch = as.character(object$branch %||% "original"),
      Cells = nrow(tbl),
      LowCountCells = sum(tbl$LowCountFlag %in% TRUE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    summary_kind = "bias_count",
    overview = summary_tbl,
    count_distribution = count_distribution,
    low_count_cells = low_tbl,
    thresholds = bundle_settings_table(object$thresholds),
    notes = if (identical(object$branch, "facets")) {
      "Legacy-compatible branch: table columns mirror the compatibility contract naming."
    } else {
      "Original branch: compact count/bias columns for QC screening."
    },
    digits = digits
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

summarize_fixed_reports_bundle <- function(object, digits = 3, top_n = 10) {
  pair_tbl <- as.data.frame(object$pairwise_table %||% data.frame(), stringsAsFactors = FALSE)
  n_bias_lines <- length(strsplit(as.character(object$bias_fixed %||% ""), "\n", fixed = TRUE)[[1]])
  n_pair_lines <- length(strsplit(as.character(object$pairwise_fixed %||% ""), "\n", fixed = TRUE)[[1]])

  overview <- data.frame(
    Branch = as.character(object$branch %||% "facets"),
    Style = as.character(object$style %||% "facets_manual"),
    Interaction = as.character(object$interaction_label %||% ""),
    PairwiseRows = nrow(pair_tbl),
    BiasTextLines = n_bias_lines,
    PairwiseTextLines = n_pair_lines,
    stringsAsFactors = FALSE
  )

  out <- list(
    summary_kind = "fixed_reports",
    overview = overview,
    summary = data.frame(),
    preview_name = if (nrow(pair_tbl) > 0) "pairwise_table" else "",
    preview = utils::head(pair_tbl, n = top_n),
    settings = data.frame(),
    notes = if (nrow(pair_tbl) == 0) {
      "No pairwise contrasts available in this interaction mode."
    } else if (identical(object$branch, "facets")) {
      "Legacy-compatible branch: fixed-width text follows the compatibility layout."
    } else {
      "Original branch: sectioned fixed-width text optimized for quick review."
    },
    digits = digits
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

summarize_visual_summaries_bundle <- function(object, digits = 3, top_n = 10) {
  warning_counts <- as.data.frame(object$warning_counts %||% data.frame(), stringsAsFactors = FALSE)
  summary_counts <- as.data.frame(object$summary_counts %||% data.frame(), stringsAsFactors = FALSE)
  crosswalk <- as.data.frame(object$crosswalk %||% data.frame(), stringsAsFactors = FALSE)

  overview <- data.frame(
    Branch = as.character(object$branch %||% "original"),
    Style = as.character(object$style %||% "original"),
    ThresholdProfile = as.character(object$threshold_profile %||% ""),
    WarningVisuals = nrow(warning_counts),
    SummaryVisuals = nrow(summary_counts),
    stringsAsFactors = FALSE
  )

  preview_tbl <- warning_counts
  if (nrow(preview_tbl) == 0) preview_tbl <- summary_counts
  preview_tbl <- utils::head(preview_tbl, n = top_n)

  notes <- if (identical(object$branch, "facets")) {
    "Legacy-compatible branch includes crosswalk metadata to compatibility-oriented output names."
  } else {
    "Original branch keeps package-native warning/summary map organization."
  }

  out <- list(
    summary_kind = "visual_summaries",
    overview = overview,
    summary = warning_counts,
    preview_name = if (nrow(preview_tbl) > 0) "warning_counts" else "",
    preview = preview_tbl,
    settings = crosswalk,
    notes = notes,
    digits = digits,
    summary_counts = summary_counts
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

bundle_component_table <- function(object, name) {
  if (!is.list(object) || is.null(name) || !nzchar(name) || !name %in% names(object)) {
    return(data.frame())
  }
  value <- object[[name]]
  if (!is.data.frame(value)) return(data.frame())
  as.data.frame(value, stringsAsFactors = FALSE)
}

bundle_first_table <- function(object, candidates, top_n = 10L) {
  candidates <- as.character(candidates %||% character(0))
  if (length(candidates) == 0) {
    return(list(name = NA_character_, table = data.frame()))
  }
  for (nm in candidates) {
    tbl <- bundle_component_table(object, nm)
    if (nrow(tbl) > 0) {
      return(list(name = nm, table = utils::head(tbl, n = top_n)))
    }
  }
  for (nm in candidates) {
    tbl <- bundle_component_table(object, nm)
    if (ncol(tbl) > 0) {
      return(list(name = nm, table = tbl))
    }
  }
  list(name = NA_character_, table = data.frame())
}

bundle_known_overview <- function(object, obj_class, preview_name, preview_rows) {
  comp_names <- names(object)
  if (is.null(comp_names)) comp_names <- character(0)
  data.frame(
    Class = obj_class,
    Components = length(comp_names),
    ComponentNames = if (length(comp_names) == 0) "" else paste(comp_names, collapse = ", "),
    PreviewComponent = ifelse(is.na(preview_name), "", preview_name),
    PreviewRows = as.integer(preview_rows),
    stringsAsFactors = FALSE
  )
}

summarize_known_bundle <- function(object,
                                   obj_class,
                                   summary_candidates = "summary",
                                   preview_candidates = NULL,
                                   settings_candidates = "settings",
                                   notes = NULL,
                                   digits = 3,
                                   top_n = 10,
                                   summary_override = NULL) {
  top_n <- max(1L, as.integer(top_n))

  summary_tbl <- if (!is.null(summary_override)) {
    as.data.frame(summary_override, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  if (nrow(summary_tbl) == 0 && ncol(summary_tbl) == 0) {
    summary_pick <- bundle_first_table(object, summary_candidates, top_n = 1L)
    summary_tbl <- summary_pick$table
  }

  preview_pick <- bundle_first_table(object, preview_candidates, top_n = top_n)
  if (is.na(preview_pick$name) || nrow(preview_pick$table) == 0) {
    preview_pick <- bundle_preview_table(object, top_n = top_n)
  }

  settings_tbl <- data.frame()
  for (nm in as.character(settings_candidates %||% character(0))) {
    if (!nm %in% names(object)) next
    value <- object[[nm]]
    if (is.data.frame(value)) {
      settings_tbl <- as.data.frame(value, stringsAsFactors = FALSE)
      break
    }
    if (is.list(value)) {
      settings_tbl <- bundle_settings_table(value)
      break
    }
  }

  notes <- as.character(notes %||% "")
  notes <- notes[nzchar(notes)]
  if (length(notes) == 0) {
    if (nrow(summary_tbl) > 0 && nrow(preview_pick$table) > 0) {
      notes <- "Summary and preview tables were extracted for this bundle."
    } else if (nrow(preview_pick$table) > 0) {
      notes <- "Preview rows were extracted from the main table component."
    } else {
      notes <- "No populated table components were found in this bundle."
    }
  }

  out <- list(
    summary_kind = obj_class,
    overview = bundle_known_overview(
      object = object,
      obj_class = obj_class,
      preview_name = preview_pick$name,
      preview_rows = nrow(preview_pick$table)
    ),
    summary = summary_tbl,
    preview_name = preview_pick$name,
    preview = preview_pick$table,
    settings = settings_tbl,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

summarize_measurable_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_measurable",
    summary_candidates = "summary",
    preview_candidates = c("facet_coverage", "category_stats", "subsets"),
    settings_candidates = character(0),
    notes = "Measurable-data summary with facet coverage, category diagnostics, and subset/connectivity checks.",
    digits = digits,
    top_n = top_n
  )
}

summarize_unexpected_after_bias_bundle <- function(object, digits = 3, top_n = 10) {
  facet_note <- if (!is.null(object$facets) && length(object$facets) > 0) {
    paste("Bias interaction:", paste(as.character(object$facets), collapse = " x "))
  } else {
    "Bias interaction facets are not attached in this object."
  }
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_unexpected_after_bias",
    summary_candidates = "summary",
    preview_candidates = "table",
    settings_candidates = "thresholds",
    notes = c(
      "Unexpected-response summary after interaction adjustment.",
      facet_note
    ),
    digits = digits,
    top_n = top_n
  )
}

summarize_output_bundle <- function(object, digits = 3, top_n = 10) {
  settings <- object$settings %||% list()
  summary_tbl <- data.frame(
    GraphRows = nrow(bundle_component_table(object, "graphfile")),
    ScoreRows = nrow(bundle_component_table(object, "scorefile")),
    WrittenFiles = nrow(bundle_component_table(object, "written_files")),
    IncludeFixed = as.logical(settings$include_fixed %||% FALSE),
    WriteFiles = as.logical(settings$write_files %||% FALSE),
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_output_bundle",
    summary_candidates = character(0),
    preview_candidates = c("scorefile", "graphfile", "graphfile_syntactic", "written_files"),
    settings_candidates = "settings",
    notes = "Graphfile/SCORE-style export bundle (table output and optional file-write metadata).",
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_residual_pca_bundle <- function(object, digits = 3, top_n = 10) {
  mode <- as.character(object$mode %||% "unknown")
  facet_names <- as.character(object$facet_names %||% character(0))
  summary_tbl <- data.frame(
    Mode = mode,
    Facets = length(facet_names),
    OverallComponents = nrow(bundle_component_table(object, "overall_table")),
    FacetComponentRows = nrow(bundle_component_table(object, "by_facet_table")),
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_residual_pca",
    summary_candidates = character(0),
    preview_candidates = c("overall_table", "by_facet_table"),
    settings_candidates = character(0),
    notes = "Residual PCA summary for unidimensionality checks (overall and/or by facet).",
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_specifications_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_specifications",
    summary_candidates = "header",
    preview_candidates = c("data_spec", "facet_labels", "output_spec", "convergence_control", "anchor_summary"),
    settings_candidates = character(0),
    notes = "Model specification summary for method and run documentation.",
    digits = digits,
    top_n = top_n
  )
}

summarize_data_quality_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_data_quality",
    summary_candidates = "summary",
    preview_candidates = c("row_audit", "category_counts", "model_match", "unknown_elements"),
    settings_candidates = character(0),
    notes = "Legacy-compatible Table 2 data quality summary and row-level audit.",
    digits = digits,
    top_n = top_n
  )
}

summarize_iteration_report_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_iteration_report",
    summary_candidates = "summary",
    preview_candidates = "table",
    settings_candidates = "settings",
    notes = "Legacy-compatible Table 3 replay of estimation iterations.",
    digits = digits,
    top_n = top_n
  )
}

summarize_subset_connectivity_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_subset_connectivity",
    summary_candidates = "summary",
    preview_candidates = c("listing", "nodes"),
    settings_candidates = "settings",
    notes = "Legacy-compatible Table 6 subset/connectivity report with subset and node listings.",
    digits = digits,
    top_n = top_n
  )
}

summarize_facet_statistics_bundle <- function(object, digits = 3, top_n = 10) {
  table_tbl <- bundle_component_table(object, "table")
  range_tbl <- bundle_component_table(object, "ranges")
  precision_tbl <- bundle_component_table(object, "precision_summary")
  variability_tbl <- bundle_component_table(object, "variability_tests")
  summary_tbl <- data.frame(
    Facets = if ("Facet" %in% names(table_tbl)) length(unique(table_tbl$Facet)) else NA_integer_,
    Rows = nrow(table_tbl),
    Metrics = if ("Metric" %in% names(table_tbl)) length(unique(table_tbl$Metric)) else NA_integer_,
    PrecisionRows = nrow(precision_tbl),
    VariabilityRows = nrow(variability_tbl),
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_facet_statistics",
    summary_candidates = character(0),
    preview_candidates = c("precision_summary", "variability_tests", "table", "ranges"),
    settings_candidates = "settings",
    notes = if (nrow(precision_tbl) > 0) {
      "Facet profile summary including distribution basis, SE mode, and variability tests."
    } else if (nrow(range_tbl) > 0) {
      "Facet profile summary including range rulers."
    } else {
      "Facet profile summary."
    },
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_precision_audit_bundle <- function(object, digits = 3, top_n = 10) {
  profile_tbl <- bundle_component_table(object, "profile")
  checks_tbl <- bundle_component_table(object, "checks")
  notes_tbl <- bundle_component_table(object, "approximation_notes")

  flagged_n <- if ("Status" %in% names(checks_tbl)) {
    sum(as.character(checks_tbl$Status) %in% c("review", "warn"), na.rm = TRUE)
  } else {
    NA_integer_
  }

  summary_tbl <- data.frame(
    Method = as.character(profile_tbl$Method[1] %||% NA_character_),
    PrecisionTier = as.character(profile_tbl$PrecisionTier[1] %||% NA_character_),
    SupportsFormalInference = isTRUE(profile_tbl$SupportsFormalInference[1] %||% FALSE),
    Checks = nrow(checks_tbl),
    ReviewOrWarn = flagged_n,
    NoteRows = nrow(notes_tbl),
    stringsAsFactors = FALSE
  )

  notes <- if (nrow(profile_tbl) > 0 && identical(as.character(profile_tbl$PrecisionTier[1]), "exploratory")) {
    "Exploratory precision path detected; use this run for screening and calibration triage, not as the package's primary inferential summary."
  } else if (nrow(profile_tbl) > 0 && identical(as.character(profile_tbl$PrecisionTier[1]), "hybrid")) {
    "Hybrid precision path detected; at least one level fell back to observation-table information, so formal inference should be limited to the model-based rows."
  } else {
    "Model-based precision path detected for the current run."
  }

  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_precision_audit",
    summary_candidates = character(0),
    preview_candidates = c("checks", "approximation_notes", "profile"),
    settings_candidates = "settings",
    notes = notes,
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_parity_bundle <- function(object, digits = 3, top_n = 10) {
  overall_tbl <- as.data.frame(object$overall %||% data.frame(), stringsAsFactors = FALSE)
  missing_tbl <- as.data.frame(object$missing_preview %||% data.frame(), stringsAsFactors = FALSE)
  metric_summary <- as.data.frame(object$metric_summary %||% data.frame(), stringsAsFactors = FALSE)

  notes <- character(0)
  if (nrow(overall_tbl) > 0) {
    mismatch <- suppressWarnings(as.integer(overall_tbl$ColumnMismatches[1]))
    if (is.finite(mismatch) && mismatch == 0) {
      notes <- c(notes, "All contract rows reached full column coverage.")
    } else if (is.finite(mismatch)) {
      notes <- c(notes, paste0("Column mismatches detected: ", mismatch, "."))
    }
  }
  if (nrow(metric_summary) > 0) {
    failed <- suppressWarnings(as.integer(metric_summary$Failed[1]))
    if (is.finite(failed) && failed == 0) {
      notes <- c(notes, "All evaluated metric checks passed.")
    } else if (is.finite(failed)) {
      notes <- c(notes, paste0("Metric checks failed: ", failed, "."))
    }
  }
  if (length(notes) == 0) {
    notes <- "Parity checks completed."
  }

  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_parity_report",
    summary_candidates = character(0),
    preview_candidates = c("missing_preview", "column_audit", "metric_audit"),
    settings_candidates = "settings",
    notes = notes,
    digits = digits,
    top_n = top_n,
    summary_override = overall_tbl
  )
}

summarize_unexpected_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_unexpected",
    summary_candidates = "summary",
    preview_candidates = "table",
    settings_candidates = "thresholds",
    notes = "Unexpected-response summary for quick residual screening.",
    digits = digits,
    top_n = top_n
  )
}

summarize_fair_average_bundle <- function(object, digits = 3, top_n = 10) {
  stacked <- bundle_component_table(object, "stacked")
  obs_avg <- if ("Obsvd Average" %in% names(stacked)) suppressWarnings(as.numeric(stacked[["Obsvd Average"]])) else numeric(0)
  fair_m <- if ("Fair(M) Average" %in% names(stacked)) suppressWarnings(as.numeric(stacked[["Fair(M) Average"]])) else numeric(0)
  mean_abs_gap <- NA_real_
  if (length(obs_avg) == length(fair_m) && length(obs_avg) > 0) {
    dif <- abs(obs_avg - fair_m)
    dif <- dif[is.finite(dif)]
    if (length(dif) > 0) mean_abs_gap <- mean(dif)
  }
  summary_tbl <- data.frame(
    Facets = if ("Facet" %in% names(stacked)) length(unique(as.character(stacked$Facet))) else length(object$by_facet %||% list()),
    Levels = nrow(stacked),
    MeanAbsObservedFairM = mean_abs_gap,
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_fair_average",
    summary_candidates = character(0),
    preview_candidates = c("stacked", "raw_by_facet"),
    settings_candidates = "settings",
    notes = "Adjusted-score reference summary by facet level.",
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_displacement_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_displacement",
    summary_candidates = "summary",
    preview_candidates = "table",
    settings_candidates = "thresholds",
    notes = "Displacement summary for anchor drift and baseline drift checks.",
    digits = digits,
    top_n = top_n
  )
}

summarize_interrater_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_interrater",
    summary_candidates = "summary",
    preview_candidates = "pairs",
    settings_candidates = "settings",
    notes = "Inter-rater agreement summary across matched scoring contexts; severity spread is reported separately from agreement when available.",
    digits = digits,
    top_n = top_n
  )
}

summarize_facets_chisq_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_facets_chisq",
    summary_candidates = "summary",
    preview_candidates = "table",
    settings_candidates = "thresholds",
    notes = "Facet variability summary with fixed/random reference tests.",
    digits = digits,
    top_n = top_n
  )
}

summarize_bias_interaction_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_bias_interaction",
    summary_candidates = "summary",
    preview_candidates = c("ranked_table", "facet_profile"),
    settings_candidates = "thresholds",
    notes = "Bias interaction report with ranked cells and facet-level profiles.",
    digits = digits,
    top_n = top_n
  )
}

summarize_rating_scale_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_rating_scale",
    summary_candidates = "summary",
    preview_candidates = c("category_table", "threshold_table"),
    settings_candidates = character(0),
    notes = "Rating-scale diagnostics with category usage, fit, and threshold ordering.",
    digits = digits,
    top_n = top_n
  )
}

summarize_category_structure_bundle <- function(object, digits = 3, top_n = 10) {
  cat_tbl <- bundle_component_table(object, "category_table")
  flags <- integer(0)
  for (nm in c("LowCount", "InfitFlag", "OutfitFlag", "ZSTDFlag")) {
    if (nm %in% names(cat_tbl)) {
      v <- as.logical(cat_tbl[[nm]])
      flags <- c(flags, sum(v, na.rm = TRUE))
    }
  }
  summary_tbl <- data.frame(
    Categories = nrow(cat_tbl),
    UsedCategories = if ("Count" %in% names(cat_tbl)) sum(suppressWarnings(as.numeric(cat_tbl$Count)) > 0, na.rm = TRUE) else NA_integer_,
    FlaggedStats = if (length(flags) > 0) sum(flags, na.rm = TRUE) else NA_integer_,
    ModeBoundaries = nrow(bundle_component_table(object, "mode_boundaries")),
    MeanHalfscorePoints = nrow(bundle_component_table(object, "mean_halfscore_points")),
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_category_structure",
    summary_candidates = character(0),
    preview_candidates = c("category_table", "mode_boundaries", "mean_halfscore_points"),
    settings_candidates = "settings",
    notes = "Category-structure diagnostics with mode boundaries and half-score reference points.",
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_category_curves_bundle <- function(object, digits = 3, top_n = 10) {
  graph_tbl <- bundle_component_table(object, "graphfile")
  prob_cols <- grep("^Prob:", names(graph_tbl), value = TRUE)
  summary_tbl <- data.frame(
    Rows = nrow(graph_tbl),
    CurveGroups = if ("CurveGroup" %in% names(graph_tbl)) length(unique(as.character(graph_tbl$CurveGroup))) else NA_integer_,
    ThetaPoints = if ("Scale" %in% names(graph_tbl)) length(unique(suppressWarnings(as.numeric(graph_tbl$Scale)))) else NA_integer_,
    ProbabilityColumns = length(prob_cols),
    stringsAsFactors = FALSE
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_category_curves",
    summary_candidates = character(0),
    preview_candidates = c("expected_ogive", "graphfile", "probabilities"),
    settings_candidates = "settings",
    notes = "Expected-score and category-probability curve bundle for scale-structure review.",
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

#' Summarize report/table bundles in a user-friendly format
#'
#' @param object Any report bundle produced by `mfrmr` table/report helpers.
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of preview rows shown from the main table component.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method provides a compact summary for bundle-like outputs
#' (for example: unexpected-response, fair-average, chi-square, and
#' category report objects). It extracts:
#' - object class and available components
#' - one-row summary table when available
#' - preview rows from the main data component
#' - resolved settings/options
#'
#' Branch-aware summaries are provided for:
#' - `mfrm_bias_count` (`branch = "original"` / `"facets"`)
#' - `mfrm_fixed_reports` (`branch = "original"` / `"facets"`)
#' - `mfrm_visual_summaries` (`branch = "original"` / `"facets"`)
#'
#' Additional class-aware summaries are provided for:
#' - `mfrm_unexpected`, `mfrm_fair_average`, `mfrm_displacement`
#' - `mfrm_interrater`, `mfrm_facets_chisq`, `mfrm_bias_interaction`
#' - `mfrm_rating_scale`, `mfrm_category_structure`, `mfrm_category_curves`
#' - `mfrm_measurable`, `mfrm_unexpected_after_bias`, `mfrm_output_bundle`
#' - `mfrm_residual_pca`, `mfrm_specifications`, `mfrm_data_quality`
#' - `mfrm_iteration_report`, `mfrm_subset_connectivity`, `mfrm_facet_statistics`
#' - `mfrm_parity_report`
#'
#' @section Interpreting output:
#' - `overview`: class, component count, and selected preview component.
#' - `summary`: one-row aggregate block when supplied by the bundle.
#' - `preview`: first `top_n` rows from the main table-like component.
#' - `settings`: resolved option values if available.
#'
#' @section Typical workflow:
#' 1. Generate a bundle table/report helper output.
#' 2. Run `summary(bundle)` for compact QA.
#' 3. Drill into specific components via `$` and visualize with `plot(bundle, ...)`.
#'
#' @return An object of class `summary.mfrm_bundle`.
#' @seealso [unexpected_response_table()], [fair_average_table()], `plot()`
#' @examples
#' \donttest{
#' toy_full <- load_mfrmr_data("example_core")
#' toy_people <- unique(toy_full$Person)[1:12]
#' toy <- toy_full[toy_full$Person %in% toy_people, , drop = FALSE]
#' fit <- suppressWarnings(
#'   fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 10)
#' )
#' t4 <- unexpected_response_table(fit, abs_z_min = 1.5, prob_max = 0.4, top_n = 5)
#' summary(t4)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' t11 <- bias_count_table(bias, branch = "facets")
#' summary(t11)
#' }
#' @export
summary.mfrm_bundle <- function(object, digits = 3, top_n = 10, ...) {
  if (!is.list(object)) {
    stop("`object` must be a bundle-like list output.")
  }
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  if (inherits(object, "mfrm_bias_count")) {
    return(summarize_bias_count_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_fixed_reports")) {
    return(summarize_fixed_reports_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_visual_summaries")) {
    return(summarize_visual_summaries_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_unexpected")) {
    return(summarize_unexpected_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_fair_average")) {
    return(summarize_fair_average_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_displacement")) {
    return(summarize_displacement_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_interrater")) {
    return(summarize_interrater_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_facets_chisq")) {
    return(summarize_facets_chisq_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_bias_interaction")) {
    return(summarize_bias_interaction_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_rating_scale")) {
    return(summarize_rating_scale_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_category_structure")) {
    return(summarize_category_structure_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_category_curves")) {
    return(summarize_category_curves_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_measurable")) {
    return(summarize_measurable_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_unexpected_after_bias")) {
    return(summarize_unexpected_after_bias_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_output_bundle")) {
    return(summarize_output_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_residual_pca")) {
    return(summarize_residual_pca_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_specifications")) {
    return(summarize_specifications_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_data_quality")) {
    return(summarize_data_quality_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_iteration_report")) {
    return(summarize_iteration_report_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_subset_connectivity")) {
    return(summarize_subset_connectivity_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_facet_statistics")) {
    return(summarize_facet_statistics_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_precision_audit")) {
    return(summarize_precision_audit_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_parity_report")) {
    return(summarize_parity_bundle(object, digits = digits, top_n = top_n))
  }

  cls <- class(object)
  cls <- cls[!cls %in% c("list", "mfrm_bundle")]
  obj_class <- if (length(cls) == 0) "mfrm_bundle" else cls[1]

  comp_names <- names(object)
  if (is.null(comp_names)) comp_names <- character(0)

  summary_tbl <- if ("summary" %in% comp_names && is.data.frame(object$summary)) {
    as.data.frame(object$summary, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  preview <- bundle_preview_table(object, top_n = top_n)
  settings_tbl <- if ("settings" %in% comp_names) bundle_settings_table(object$settings) else data.frame()

  overview <- data.frame(
    Class = obj_class,
    Components = length(comp_names),
    ComponentNames = if (length(comp_names) == 0) "" else paste(comp_names, collapse = ", "),
    PreviewComponent = ifelse(is.na(preview$name), "", preview$name),
    PreviewRows = nrow(preview$table),
    stringsAsFactors = FALSE
  )

  notes <- if (nrow(summary_tbl) > 0) {
    "Summary table and preview rows were extracted."
  } else if (nrow(preview$table) > 0) {
    "No `summary` component found; showing preview rows from the main table."
  } else {
    "No tabular components available for preview."
  }

  out <- list(
    overview = overview,
    summary = summary_tbl,
    preview_name = preview$name,
    preview = preview$table,
    settings = settings_tbl,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

bundle_summary_labels <- function(summary_kind, overview = NULL) {
  class_name <- NA_character_
  if (!is.null(overview) && is.data.frame(overview) && nrow(overview) > 0 && "Class" %in% names(overview)) {
    class_name <- as.character(overview$Class[1])
  }
  key <- as.character(summary_kind %||% class_name %||% "")
  if (!nzchar(key) || identical(key, "NA")) key <- "mfrm_bundle"

  defaults <- list(
    title = "mfrmr Bundle Summary",
    summary = "Summary table",
    preview = "Preview",
    settings = "Settings"
  )

  maps <- list(
    mfrm_unexpected = list(title = "mfrmr Unexpected Response Summary", summary = "Threshold summary", preview = "Flagged responses"),
    mfrm_fair_average = list(title = "mfrmr Adjusted Score Summary", summary = "Overview", preview = "Facet-level adjusted-score rows"),
    mfrm_displacement = list(title = "mfrmr Displacement Summary", summary = "Displacement summary", preview = "Displacement rows"),
    mfrm_interrater = list(title = "mfrmr Agreement Summary", summary = "Agreement summary", preview = "Rater-pair rows"),
    mfrm_facets_chisq = list(title = "mfrmr Facet Variability Summary", summary = "Facet variability summary", preview = "Facet rows"),
    mfrm_bias_interaction = list(title = "mfrmr Bias Interaction Summary", summary = "Interaction summary", preview = "Ranked interaction rows"),
    mfrm_bias_iteration = list(title = "mfrmr Bias Iteration Summary", summary = "Iteration summary", preview = "Iteration rows"),
    mfrm_bias_pairwise = list(title = "mfrmr Bias Pairwise Summary", summary = "Pairwise summary", preview = "Contrast rows"),
    mfrm_rating_scale = list(title = "mfrmr Rating Scale Summary", summary = "Category/threshold summary", preview = "Category rows"),
    mfrm_category_structure = list(title = "mfrmr Category Structure Summary", summary = "Category structure overview", preview = "Category structure rows"),
    mfrm_category_curves = list(title = "mfrmr Category Curves Summary", summary = "Curve grid summary", preview = "Expected-score / curve rows"),
    mfrm_measurable = list(title = "mfrmr Measurable Summary", summary = "Run overview", preview = "Facet/category rows"),
    mfrm_unexpected_after_bias = list(title = "mfrmr Unexpected-after-Bias Summary", summary = "After-bias threshold summary", preview = "After-bias flagged rows"),
    mfrm_output_bundle = list(title = "mfrmr Output File Bundle Summary", summary = "Output overview", preview = "Output preview rows"),
    mfrm_residual_pca = list(title = "mfrmr Residual PCA Summary", summary = "PCA overview", preview = "Eigenvalue / loading rows"),
    mfrm_specifications = list(title = "mfrmr Specifications Summary", summary = "Specification header", preview = "Specification rows"),
    mfrm_data_quality = list(title = "mfrmr Data Quality Summary", summary = "Data quality overview", preview = "Audit rows"),
    mfrm_iteration_report = list(title = "mfrmr Iteration Report Summary", summary = "Iteration overview", preview = "Iteration rows"),
    mfrm_subset_connectivity = list(title = "mfrmr Subset Connectivity Summary", summary = "Subset overview", preview = "Subset/node rows"),
    mfrm_facet_statistics = list(title = "mfrmr Facet Profile Summary", summary = "Facet-profile overview", preview = "Facet-profile rows"),
    mfrm_precision_audit = list(title = "mfrmr Precision Audit Summary", summary = "Precision overview", preview = "Audit checks"),
    mfrm_parity_report = list(title = "mfrmr Compatibility Contract Audit Summary", summary = "Compatibility audit overview", preview = "Lowest-coverage contract items"),
    mfrm_reference_audit = list(title = "mfrmr Internal Reference Audit Summary", summary = "Internal audit overview", preview = "Attention items"),
    mfrm_reference_benchmark = list(title = "mfrmr Internal Benchmark Summary", summary = "Case audit summary", preview = "Internal benchmark fit runs"),
    mfrm_reporting_checklist = list(title = "mfrmr Reporting Checklist Summary", summary = "Checklist coverage", preview = "Checklist items"),
    mfrm_bias_collection = list(title = "mfrmr Bias Collection Summary", summary = "Interaction summary", preview = "Per-pair results"),
    mfrm_manifest = list(title = "mfrmr Manifest Summary", summary = "Analysis overview", preview = "Manifest tables"),
    mfrm_replay_script = list(title = "mfrmr Replay Script Summary", summary = "Replay settings", preview = "Script text"),
    mfrm_export_bundle = list(title = "mfrmr Export Bundle Summary", summary = "Export overview", preview = "Written files"),
    mfrm_facet_equivalence = list(title = "mfrmr Facet Equivalence Summary", summary = "Equivalence overview", preview = "Pairwise / ROPE rows")
  )

  if (key %in% names(maps)) {
    out <- utils::modifyList(defaults, maps[[key]])
  } else {
    out <- defaults
  }
  out
}

print_bundle_section <- function(title, table, digits = 3, round_numeric = TRUE) {
  if (is.null(table) || !is.data.frame(table) || nrow(table) == 0) return(invisible(NULL))
  cat("\n", title, "\n", sep = "")
  if (isTRUE(round_numeric)) {
    print(round_numeric_df(as.data.frame(table), digits = digits), row.names = FALSE)
  } else {
    print(as.data.frame(table), row.names = FALSE)
  }
  invisible(NULL)
}

#' @export
print.summary.mfrm_bundle <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  if (identical(x$summary_kind, "bias_count")) {
    cat("mfrmr Bias Count Summary\n")
    if (!is.null(x$overview) && nrow(x$overview) > 0) {
      cat("\nOverview\n")
      print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$count_distribution) && nrow(x$count_distribution) > 0) {
      cat("\nCount distribution\n")
      print(round_numeric_df(as.data.frame(x$count_distribution), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$low_count_cells) && nrow(x$low_count_cells) > 0) {
      cat("\nLow-count cells (preview)\n")
      print(round_numeric_df(as.data.frame(x$low_count_cells), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$thresholds) && nrow(x$thresholds) > 0) {
      cat("\nThresholds\n")
      print(as.data.frame(x$thresholds), row.names = FALSE)
    }
    if (length(x$notes) > 0) {
      cat("\nNotes\n")
      cat(" - ", x$notes, "\n", sep = "")
    }
    return(invisible(x))
  }

  if (identical(x$summary_kind, "visual_summaries")) {
    cat("mfrmr Visual Summary Bundle\n")
    if (!is.null(x$overview) && nrow(x$overview) > 0) {
      cat("\nOverview\n")
      print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$summary) && nrow(x$summary) > 0) {
      cat("\nWarning counts\n")
      print(round_numeric_df(as.data.frame(x$summary), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$summary_counts) && nrow(x$summary_counts) > 0) {
      cat("\nSummary counts\n")
      print(round_numeric_df(as.data.frame(x$summary_counts), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$settings) && nrow(x$settings) > 0) {
      cat("\nFACETS crosswalk\n")
      print(as.data.frame(x$settings), row.names = FALSE)
    }
    if (length(x$notes) > 0) {
      cat("\nNotes\n")
      cat(" - ", x$notes, "\n", sep = "")
    }
    return(invisible(x))
  }

  if (identical(x$summary_kind, "fixed_reports")) {
    cat("mfrmr Fixed-Report Bundle\n")
    if (!is.null(x$overview) && nrow(x$overview) > 0) {
      cat("\nOverview\n")
      print(round_numeric_df(as.data.frame(x$overview), digits = digits), row.names = FALSE)
    }
    if (!is.null(x$preview) && nrow(x$preview) > 0) {
      cat("\nPairwise preview\n")
      print(round_numeric_df(as.data.frame(x$preview), digits = digits), row.names = FALSE)
    }
    if (length(x$notes) > 0) {
      cat("\nNotes\n")
      cat(" - ", x$notes, "\n", sep = "")
    }
    return(invisible(x))
  }

  labels <- bundle_summary_labels(summary_kind = x$summary_kind, overview = x$overview)
  cat(labels$title, "\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    ov <- x$overview[1, , drop = FALSE]
    cat(sprintf("  Class: %s\n", ov$Class))
    cat(sprintf("  Components (%s): %s\n", ov$Components, ov$ComponentNames))
  }
  print_bundle_section(labels$summary, x$summary, digits = digits, round_numeric = TRUE)
  if (!is.null(x$preview) && nrow(x$preview) > 0) {
    preview_title <- labels$preview
    if (!is.null(x$preview_name) && !is.na(x$preview_name) && nzchar(x$preview_name)) {
      preview_title <- paste0(preview_title, ": ", x$preview_name)
    }
    print_bundle_section(preview_title, x$preview, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$settings) && nrow(x$settings) > 0) {
    print_bundle_section(labels$settings, x$settings, digits = digits, round_numeric = FALSE)
  }
  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    cat(" - ", x$notes, "\n", sep = "")
  }
  invisible(x)
}

draw_category_structure_bundle <- function(x,
                                           type = c("counts", "mode_boundaries", "mean_halfscore"),
                                           draw = TRUE,
                                           main = NULL,
                                           palette = NULL,
                                           label_angle = 45) {
  type <- match.arg(tolower(type), c("counts", "mode_boundaries", "mean_halfscore"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      counts = "#9ecae1",
      expected = "#08519c",
      mode = "#2b8cbe",
      mean = "#238b45"
    )
  )
  cat_tbl <- as.data.frame(x$category_table %||% data.frame(), stringsAsFactors = FALSE)
  mode_tbl <- as.data.frame(x$mode_boundaries %||% data.frame(), stringsAsFactors = FALSE)
  half_tbl <- as.data.frame(x$mean_halfscore_points %||% data.frame(), stringsAsFactors = FALSE)

  if (isTRUE(draw)) {
    if (type == "counts") {
      if (nrow(cat_tbl) == 0 || !all(c("Category", "Count") %in% names(cat_tbl))) stop("No category count data available.")
      bp <- barplot_rot45(
        height = suppressWarnings(as.numeric(cat_tbl$Count)),
        labels = as.character(cat_tbl$Category),
        col = pal["counts"],
        main = if (is.null(main)) "Category counts" else as.character(main[1]),
        ylab = "Count",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
      if ("ExpectedCount" %in% names(cat_tbl)) {
        exp_ct <- suppressWarnings(as.numeric(cat_tbl$ExpectedCount))
        if (any(is.finite(exp_ct))) {
          graphics::points(bp, exp_ct, pch = 21, bg = "white", col = pal["expected"])
          graphics::lines(bp, exp_ct, col = pal["expected"], lwd = 1.4)
        }
      }
    } else if (type == "mode_boundaries") {
      if (nrow(mode_tbl) == 0 || !all(c("CurveGroup", "ModeBoundaryTheta") %in% names(mode_tbl))) {
        stop("No mode-boundary data available.")
      }
      grp <- as.factor(mode_tbl$CurveGroup)
      y <- as.numeric(grp)
      graphics::plot(
        x = suppressWarnings(as.numeric(mode_tbl$ModeBoundaryTheta)),
        y = y,
        pch = 16,
        col = pal["mode"],
        xlab = "Theta / Logit",
        ylab = "",
        yaxt = "n",
        main = if (is.null(main)) "Mode boundaries" else as.character(main[1])
      )
      graphics::axis(side = 2, at = seq_along(levels(grp)), labels = levels(grp), las = 2)
    } else {
      if (nrow(half_tbl) == 0 || !all(c("CurveGroup", "MeanBoundaryTheta") %in% names(half_tbl))) {
        stop("No mean half-score data available.")
      }
      grp <- as.factor(half_tbl$CurveGroup)
      y <- as.numeric(grp)
      graphics::plot(
        x = suppressWarnings(as.numeric(half_tbl$MeanBoundaryTheta)),
        y = y,
        pch = 16,
        col = pal["mean"],
        xlab = "Theta / Logit",
        ylab = "",
        yaxt = "n",
        main = if (is.null(main)) "Mean half-score boundaries" else as.character(main[1])
      )
      graphics::axis(side = 2, at = seq_along(levels(grp)), labels = levels(grp), las = 2)
    }
  }

  new_mfrm_plot_data(
    "category_structure",
    list(
      plot = type,
      category_table = cat_tbl,
      mode_boundaries = mode_tbl,
      mean_halfscore_points = half_tbl
    )
  )
}

draw_category_curves_bundle <- function(x,
                                        type = c("ogive", "ccc"),
                                        draw = TRUE,
                                        main = NULL,
                                        palette = NULL) {
  type <- match.arg(tolower(type), c("ogive", "ccc"))
  ogive <- as.data.frame(x$expected_ogive %||% data.frame(), stringsAsFactors = FALSE)
  probs <- as.data.frame(x$probabilities %||% data.frame(), stringsAsFactors = FALSE)

  if (isTRUE(draw)) {
    if (type == "ogive") {
      if (nrow(ogive) == 0 || !all(c("Theta", "ExpectedScore", "CurveGroup") %in% names(ogive))) {
        stop("No expected-ogive data available.")
      }
      groups <- unique(as.character(ogive$CurveGroup))
      defaults <- stats::setNames(grDevices::hcl.colors(max(3L, length(groups)), "Dark 3")[seq_along(groups)], groups)
      cols <- resolve_palette(palette = palette, defaults = defaults)
      graphics::plot(
        x = range(ogive$Theta, finite = TRUE),
        y = range(ogive$ExpectedScore, finite = TRUE),
        type = "n",
        xlab = "Theta / Logit",
        ylab = "Expected score",
        main = if (is.null(main)) "Expected-score ogive" else as.character(main[1])
      )
      for (i in seq_along(groups)) {
        sub <- ogive[ogive$CurveGroup == groups[i], , drop = FALSE]
        graphics::lines(sub$Theta, sub$ExpectedScore, col = cols[groups[i]], lwd = 2)
      }
      graphics::legend("topleft", legend = groups, col = cols[groups], lty = 1, lwd = 2, bty = "n")
    } else {
      if (nrow(probs) == 0 || !all(c("Theta", "Probability", "Category", "CurveGroup") %in% names(probs))) {
        stop("No category-curve data available.")
      }
      traces <- unique(paste(probs$CurveGroup, probs$Category, sep = " | Cat "))
      defaults <- stats::setNames(grDevices::hcl.colors(max(3L, length(traces)), "Dark 3")[seq_along(traces)], traces)
      cols <- resolve_palette(palette = palette, defaults = defaults)
      graphics::plot(
        x = range(probs$Theta, finite = TRUE),
        y = c(0, 1),
        type = "n",
        xlab = "Theta / Logit",
        ylab = "Probability",
        main = if (is.null(main)) "Category characteristic curves" else as.character(main[1])
      )
      for (i in seq_along(traces)) {
        parts <- strsplit(traces[i], " \\| Cat ")[[1]]
        sub <- probs[probs$CurveGroup == parts[1] & probs$Category == parts[2], , drop = FALSE]
        graphics::lines(sub$Theta, sub$Probability, col = cols[traces[i]], lwd = 1.4)
      }
    }
  }

  new_mfrm_plot_data(
    "category_curves",
    list(plot = type, expected_ogive = ogive, probabilities = probs)
  )
}

draw_rating_scale_bundle <- function(x,
                                     type = c("counts", "thresholds"),
                                     draw = TRUE,
                                     main = NULL,
                                     palette = NULL,
                                     label_angle = 45) {
  type <- match.arg(tolower(type), c("counts", "thresholds"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      counts = "#c7e9c0",
      expected = "#08519c",
      step_line = "#1B9E77"
    )
  )
  cat_tbl <- as.data.frame(x$category_table %||% data.frame(), stringsAsFactors = FALSE)
  thr_tbl <- as.data.frame(x$threshold_table %||% data.frame(), stringsAsFactors = FALSE)

  if (isTRUE(draw)) {
    if (type == "counts") {
      if (nrow(cat_tbl) == 0 || !all(c("Category", "Count") %in% names(cat_tbl))) {
        stop("No category count data available.")
      }
      bp <- barplot_rot45(
        height = suppressWarnings(as.numeric(cat_tbl$Count)),
        labels = as.character(cat_tbl$Category),
        col = pal["counts"],
        main = if (is.null(main)) "Rating-scale category counts" else as.character(main[1]),
        ylab = "Count",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
      if ("ExpectedCount" %in% names(cat_tbl)) {
        exp_ct <- suppressWarnings(as.numeric(cat_tbl$ExpectedCount))
        if (any(is.finite(exp_ct))) {
          graphics::points(bp, exp_ct, pch = 21, bg = "white", col = pal["expected"])
          graphics::lines(bp, exp_ct, col = pal["expected"], lwd = 1.3)
        }
      }
    } else {
      if (nrow(thr_tbl) == 0 || !all(c("Step", "Estimate") %in% names(thr_tbl))) {
        stop("No threshold data available.")
      }
      draw_step_plot(
        thr_tbl,
        title = if (is.null(main)) "Rating-scale thresholds" else as.character(main[1]),
        palette = c(step_line = pal["step_line"]),
        label_angle = label_angle
      )
    }
  }

  new_mfrm_plot_data(
    "rating_scale",
    list(plot = type, category_table = cat_tbl, threshold_table = thr_tbl)
  )
}

draw_measurable_bundle <- function(x,
                                   type = c("facet_coverage", "category_counts", "subset_observations"),
                                   draw = TRUE,
                                   main = NULL,
                                   palette = NULL,
                                   label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("facet_coverage", "category_counts", "subset_observations"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      facet = "#2b8cbe",
      category = "#31a354",
      subset = "#756bb1"
    )
  )
  facet_tbl <- as.data.frame(x$facet_coverage %||% data.frame(), stringsAsFactors = FALSE)
  cat_tbl <- as.data.frame(x$category_stats %||% data.frame(), stringsAsFactors = FALSE)
  sub_tbl <- as.data.frame(x$subsets %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "facet_coverage") {
    if (nrow(facet_tbl) == 0 || !all(c("Facet", "Levels") %in% names(facet_tbl))) {
      stop("No facet-coverage table available.")
    }
    vals <- suppressWarnings(as.numeric(facet_tbl$Levels))
    labels <- as.character(facet_tbl$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["facet"],
        main = if (is.null(main)) "Facet coverage (levels per facet)" else as.character(main[1]),
        ylab = "Levels",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "measurable",
      list(plot = "facet_coverage", table = facet_tbl)
    )))
  }

  if (type == "category_counts") {
    if (nrow(cat_tbl) == 0 || !all(c("Category", "Count") %in% names(cat_tbl))) {
      stop("No category-statistics table available.")
    }
    vals <- suppressWarnings(as.numeric(cat_tbl$Count))
    labels <- as.character(cat_tbl$Category)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["category"],
        main = if (is.null(main)) "Category counts (measurable data)" else as.character(main[1]),
        ylab = "Count",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "measurable",
      list(plot = "category_counts", table = cat_tbl)
    )))
  }

  if (nrow(sub_tbl) == 0 || !all(c("Subset", "Observations") %in% names(sub_tbl))) {
    stop("No subset summary available.")
  }
  vals <- suppressWarnings(as.numeric(sub_tbl$Observations))
  labels <- paste0("Subset ", as.character(sub_tbl$Subset))
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["subset"],
      main = if (is.null(main)) "Observations by subset" else as.character(main[1]),
      ylab = "Observations",
      label_angle = label_angle,
      mar_bottom = 7.8
    )
  }
  invisible(new_mfrm_plot_data(
    "measurable",
    list(plot = "subset_observations", table = sub_tbl)
  ))
}

draw_unexpected_after_bias_bundle <- function(x,
                                              type = c("scatter", "severity", "comparison"),
                                              top_n = 40,
                                              draw = TRUE,
                                              main = NULL,
                                              palette = NULL,
                                              label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("scatter", "severity", "comparison"))
  top_n <- max(1L, as.integer(top_n))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      higher = "#d95f02",
      lower = "#1b9e77",
      severity = "#2b8cbe",
      baseline = "#9ecae1",
      after = "#3182bd"
    )
  )
  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  summary_tbl <- as.data.frame(x$summary %||% data.frame(), stringsAsFactors = FALSE)
  thr <- x$thresholds %||% list(abs_z_min = 2, prob_max = 0.30)

  if (type == "comparison") {
    if (nrow(summary_tbl) == 0) stop("No summary table available.")
    base_n <- suppressWarnings(as.numeric(summary_tbl$BaselineUnexpectedN[1] %||% NA_real_))
    after_n <- suppressWarnings(as.numeric(summary_tbl$AfterBiasUnexpectedN[1] %||% NA_real_))
    vals <- c(base_n, after_n)
    if (!all(is.finite(vals))) stop("Baseline/after-bias counts are not available.")
    labels <- c("Baseline", "After bias")
    if (isTRUE(draw)) {
      mids <- graphics::barplot(
        height = vals,
        col = c(pal["baseline"], pal["after"]),
        names.arg = labels,
        ylab = "Unexpected responses",
        main = if (is.null(main)) "Unexpected responses: baseline vs after bias" else as.character(main[1]),
        border = "white"
      )
      graphics::text(mids, vals, labels = as.integer(vals), pos = 3, cex = 0.85)
    }
    return(invisible(new_mfrm_plot_data(
      "unexpected_after_bias",
      list(plot = "comparison", baseline = base_n, after = after_n)
    )))
  }

  if (nrow(tbl) == 0) stop("No unexpected-after-bias rows available.")

  if (type == "scatter") {
    x_vals <- suppressWarnings(as.numeric(tbl$StdResidual))
    y_vals <- -log10(pmax(suppressWarnings(as.numeric(tbl$ObsProb)), .Machine$double.xmin))
    dirs <- as.character(tbl$Direction %||% rep(NA_character_, nrow(tbl)))
    cols <- ifelse(dirs == "Higher than expected", pal["higher"], pal["lower"])
    cols[!is.finite(x_vals) | !is.finite(y_vals)] <- "gray60"
    if (isTRUE(draw)) {
      graphics::plot(
        x = x_vals,
        y = y_vals,
        xlab = "Standardized residual",
        ylab = expression(-log[10](P[obs])),
        main = if (is.null(main)) "Unexpected responses after bias adjustment" else as.character(main[1]),
        pch = 16,
        col = cols
      )
      z_thr <- as.numeric(thr$abs_z_min %||% 2)
      p_thr <- as.numeric(thr$prob_max %||% 0.30)
      graphics::abline(v = c(-z_thr, z_thr), lty = 2, col = "gray45")
      graphics::abline(h = -log10(p_thr), lty = 2, col = "gray45")
      graphics::legend(
        "topleft",
        legend = c("Higher than expected", "Lower than expected"),
        col = c(pal["higher"], pal["lower"]),
        pch = 16,
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "unexpected_after_bias",
      list(plot = "scatter", table = tbl, thresholds = thr)
    )))
  }

  sev <- suppressWarnings(as.numeric(tbl$Severity))
  sev <- sev[is.finite(sev)]
  if (length(sev) == 0) stop("No finite severity values available.")
  ord <- order(suppressWarnings(as.numeric(tbl$Severity)), decreasing = TRUE, na.last = NA)
  use <- ord[seq_len(min(length(ord), top_n))]
  sub <- tbl[use, , drop = FALSE]
  labels <- if ("Row" %in% names(sub)) paste0("Row ", sub$Row) else paste0("Case ", seq_len(nrow(sub)))
  vals <- suppressWarnings(as.numeric(sub$Severity))
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["severity"],
      main = if (is.null(main)) "Unexpected-response severity after bias" else as.character(main[1]),
      ylab = "Severity",
      label_angle = label_angle,
      mar_bottom = 8.2
    )
  }
  invisible(new_mfrm_plot_data(
    "unexpected_after_bias",
    list(plot = "severity", table = sub)
  ))
}

draw_output_bundle <- function(x,
                               type = c("graph_expected", "score_residuals", "obs_probability"),
                               draw = TRUE,
                               main = NULL,
                               palette = NULL) {
  type <- match.arg(tolower(as.character(type[1])), c("graph_expected", "score_residuals", "obs_probability"))
  graph_tbl <- as.data.frame(x$graphfile %||% data.frame(), stringsAsFactors = FALSE)
  score_tbl <- as.data.frame(x$scorefile %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "graph_expected") {
    if (nrow(graph_tbl) == 0 || !all(c("Measure", "Expected") %in% names(graph_tbl))) {
      stop("Graphfile table with `Measure` and `Expected` is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    groups <- if ("CurveGroup" %in% names(graph_tbl)) unique(as.character(graph_tbl$CurveGroup)) else "All"
    if (!"CurveGroup" %in% names(graph_tbl)) graph_tbl$CurveGroup <- "All"
    defaults <- stats::setNames(grDevices::hcl.colors(max(3L, length(groups)), "Dark 3")[seq_along(groups)], groups)
    cols <- resolve_palette(palette = palette, defaults = defaults)
    if (isTRUE(draw)) {
      graphics::plot(
        x = range(graph_tbl$Measure, finite = TRUE),
        y = range(graph_tbl$Expected, finite = TRUE),
        type = "n",
        xlab = "Theta / Logit",
        ylab = "Expected score",
        main = if (is.null(main)) "Graphfile expected-score curves" else as.character(main[1])
      )
      for (g in groups) {
        sub <- graph_tbl[as.character(graph_tbl$CurveGroup) == g, , drop = FALSE]
        sub <- sub[order(sub$Measure), , drop = FALSE]
        graphics::lines(sub$Measure, sub$Expected, col = cols[g], lwd = 1.8)
      }
      if (length(groups) > 1) {
        graphics::legend("topleft", legend = groups, col = cols[groups], lty = 1, lwd = 2, bty = "n", cex = 0.85)
      }
    }
    return(invisible(new_mfrm_plot_data(
      "output_bundle",
      list(plot = "graph_expected", table = graph_tbl)
    )))
  }

  if (nrow(score_tbl) == 0) stop("Scorefile table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)

  if (type == "score_residuals") {
    if (!"Residual" %in% names(score_tbl)) stop("`Residual` column is not available in scorefile. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    vals <- suppressWarnings(as.numeric(score_tbl$Residual))
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) stop("No finite residual values available.")
    if (isTRUE(draw)) {
      graphics::hist(
        x = vals,
        breaks = "FD",
        col = "#9ecae1",
        border = "white",
        main = if (is.null(main)) "Scorefile residual distribution" else as.character(main[1]),
        xlab = "Residual",
        ylab = "Count"
      )
      graphics::abline(v = 0, lty = 2, col = "gray45")
    }
    return(invisible(new_mfrm_plot_data(
      "output_bundle",
      list(plot = "score_residuals", values = vals)
    )))
  }

  if (!"ObsProb" %in% names(score_tbl)) stop("`ObsProb` column is not available in scorefile. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  vals <- suppressWarnings(as.numeric(score_tbl$ObsProb))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) stop("No finite observed-probability values available.")
  if (isTRUE(draw)) {
    graphics::hist(
      x = vals,
      breaks = "FD",
      col = "#c7e9c0",
      border = "white",
      main = if (is.null(main)) "Observed probability distribution" else as.character(main[1]),
      xlab = "Observed probability",
      ylab = "Count"
    )
  }
  invisible(new_mfrm_plot_data(
    "output_bundle",
    list(plot = "obs_probability", values = vals)
  ))
}

draw_specifications_bundle <- function(x,
                                       type = c("facet_elements", "anchor_constraints", "convergence"),
                                       draw = TRUE,
                                       main = NULL,
                                       palette = NULL,
                                       label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("facet_elements", "anchor_constraints", "convergence"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      facet = "#2b8cbe",
      anchor = "#756bb1",
      group = "#9ecae1",
      free = "#d9d9d9",
      convergence = "#31a354"
    )
  )
  facet_tbl <- as.data.frame(x$facet_labels %||% data.frame(), stringsAsFactors = FALSE)
  anchor_tbl <- as.data.frame(x$anchor_summary %||% data.frame(), stringsAsFactors = FALSE)
  conv_tbl <- as.data.frame(x$convergence_control %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "facet_elements") {
    if (nrow(facet_tbl) == 0 || !all(c("Facet", "Elements") %in% names(facet_tbl))) {
      stop("Facet-label table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    vals <- suppressWarnings(as.numeric(facet_tbl$Elements))
    labels <- as.character(facet_tbl$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["facet"],
        main = if (is.null(main)) "Facet elements in model specification" else as.character(main[1]),
        ylab = "Elements",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "specifications",
      list(plot = "facet_elements", table = facet_tbl)
    )))
  }

  if (type == "anchor_constraints") {
    if (nrow(anchor_tbl) == 0 || !all(c("Facet", "AnchoredLevels", "GroupAnchors") %in% names(anchor_tbl))) {
      stop("Anchor summary table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    base_tbl <- anchor_tbl |>
      dplyr::transmute(
        Facet = as.character(.data$Facet),
        Anchored = suppressWarnings(as.numeric(.data$AnchoredLevels)),
        Grouped = suppressWarnings(as.numeric(.data$GroupAnchors))
      )
    if (nrow(facet_tbl) > 0 && all(c("Facet", "Elements") %in% names(facet_tbl))) {
      base_tbl <- base_tbl |>
        dplyr::left_join(
          facet_tbl |>
            dplyr::transmute(Facet = as.character(.data$Facet), Elements = suppressWarnings(as.numeric(.data$Elements))),
          by = "Facet"
        )
      base_tbl$Free <- pmax(0, base_tbl$Elements - base_tbl$Anchored - base_tbl$Grouped)
    } else {
      base_tbl$Elements <- NA_real_
      base_tbl$Free <- NA_real_
    }
    base_tbl <- base_tbl[order(base_tbl$Facet), , drop = FALSE]
    if (isTRUE(draw)) {
      old_mar <- graphics::par("mar")
      on.exit(graphics::par(mar = old_mar), add = TRUE)
      mar <- old_mar
      mar[1] <- max(mar[1], 8.8)
      graphics::par(mar = mar)
      mat <- rbind(
        Anchored = base_tbl$Anchored,
        Grouped = base_tbl$Grouped,
        Free = ifelse(is.finite(base_tbl$Free), base_tbl$Free, 0)
      )
      mids <- graphics::barplot(
        height = mat,
        beside = FALSE,
        names.arg = FALSE,
        col = c(pal["anchor"], pal["group"], pal["free"]),
        border = "white",
        ylab = "Levels",
        main = if (is.null(main)) "Anchor constraints by facet" else as.character(main[1])
      )
      draw_rotated_x_labels(
        at = mids,
        labels = base_tbl$Facet,
        srt = label_angle,
        cex = 0.82,
        line_offset = 0.085
      )
      graphics::legend(
        "topright",
        legend = c("Anchored", "Grouped", "Free"),
        fill = c(pal["anchor"], pal["group"], pal["free"]),
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "specifications",
      list(plot = "anchor_constraints", table = base_tbl)
    )))
  }

  if (nrow(conv_tbl) == 0 || !all(c("Setting", "Value") %in% names(conv_tbl))) {
    stop("Convergence-control table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  }
  keep <- c("MaxIterations", "RelativeTolerance", "QuadPoints", "FunctionEvaluations")
  sub <- conv_tbl[as.character(conv_tbl$Setting) %in% keep, , drop = FALSE]
  if (nrow(sub) == 0) stop("No numeric convergence settings found.")
  vals <- suppressWarnings(as.numeric(sub$Value))
  ok <- is.finite(vals)
  if (!any(ok)) stop("No finite numeric values in convergence settings.")
  sub <- sub[ok, , drop = FALSE]
  vals <- vals[ok]
  labels <- as.character(sub$Setting)
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["convergence"],
      main = if (is.null(main)) "Convergence controls and counts" else as.character(main[1]),
      ylab = "Value",
      label_angle = label_angle,
      mar_bottom = 8.2
    )
  }
  invisible(new_mfrm_plot_data(
    "specifications",
    list(plot = "convergence", table = sub)
  ))
}

draw_data_quality_bundle <- function(x,
                                     type = c("row_audit", "category_counts", "missing_rows"),
                                     draw = TRUE,
                                     main = NULL,
                                     palette = NULL,
                                     label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("row_audit", "category_counts", "missing_rows"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      row_audit = "#2b8cbe",
      category = "#31a354",
      missing = "#756bb1"
    )
  )
  row_tbl <- as.data.frame(x$row_audit %||% data.frame(), stringsAsFactors = FALSE)
  cat_tbl <- as.data.frame(x$category_counts %||% data.frame(), stringsAsFactors = FALSE)
  sum_tbl <- as.data.frame(x$summary %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "row_audit") {
    if (nrow(row_tbl) == 0 || !all(c("Status", "N") %in% names(row_tbl))) {
      stop("Row-audit table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    vals <- suppressWarnings(as.numeric(row_tbl$N))
    labels <- as.character(row_tbl$Status)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["row_audit"],
        main = if (is.null(main)) "Row-audit status counts" else as.character(main[1]),
        ylab = "Rows",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(plot = "row_audit", table = row_tbl)
    )))
  }

  if (type == "category_counts") {
    if (nrow(cat_tbl) == 0 || !all(c("Score", "Count") %in% names(cat_tbl))) {
      stop("Category-count table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    vals <- suppressWarnings(as.numeric(cat_tbl$Count))
    labels <- as.character(cat_tbl$Score)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["category"],
        main = if (is.null(main)) "Observed category counts" else as.character(main[1]),
        ylab = "Count",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(plot = "category_counts", table = cat_tbl)
    )))
  }

  if (nrow(sum_tbl) == 0) stop("Summary table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  row_cols <- grep("Rows$", names(sum_tbl), value = TRUE)
  if (length(row_cols) == 0) stop("No row-count columns found in summary table.")
  vals <- suppressWarnings(as.numeric(sum_tbl[1, row_cols, drop = TRUE]))
  labels <- row_cols
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["missing"],
      main = if (is.null(main)) "Missing/invalid row counts" else as.character(main[1]),
      ylab = "Rows",
      label_angle = label_angle,
      mar_bottom = 9.0
    )
  }
  invisible(new_mfrm_plot_data(
    "data_quality",
    list(plot = "missing_rows", table = data.frame(Field = labels, Count = vals, stringsAsFactors = FALSE))
  ))
}

draw_iteration_report_bundle <- function(x,
                                         type = c("residual", "logit_change", "objective"),
                                         draw = TRUE,
                                         main = NULL,
                                         palette = NULL) {
  type <- match.arg(tolower(as.character(type[1])), c("residual", "logit_change", "objective"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      residual_element = "#2b8cbe",
      residual_category = "#31a354",
      change_element = "#756bb1",
      change_step = "#d95f02",
      objective = "#1b9e77"
    )
  )
  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !"Iteration" %in% names(tbl)) {
    stop("Iteration table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  }
  it <- suppressWarnings(as.numeric(tbl$Iteration))
  if (!all(is.finite(it))) it <- seq_len(nrow(tbl))

  if (type == "residual") {
    y1 <- suppressWarnings(as.numeric(tbl$MaxScoreResidualElements))
    y2 <- suppressWarnings(as.numeric(tbl$MaxScoreResidualCategories))
    if (!any(is.finite(y1)) && !any(is.finite(y2))) {
      stop("No residual metrics available.")
    }
    if (isTRUE(draw)) {
      yr <- range(c(y1, y2), finite = TRUE)
      graphics::plot(
        x = it,
        y = y1,
        type = "b",
        pch = 16,
        col = pal["residual_element"],
        ylim = yr,
        xlab = "Iteration",
        ylab = "Residual metric",
        main = if (is.null(main)) "Iteration residual trajectory" else as.character(main[1])
      )
      graphics::lines(it, y2, type = "b", pch = 17, col = pal["residual_category"])
      graphics::legend(
        "topright",
        legend = c("Elements", "Categories"),
        col = c(pal["residual_element"], pal["residual_category"]),
        pch = c(16, 17),
        lty = 1,
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "iteration_report",
      list(plot = "residual", iteration = it, element = y1, category = y2)
    )))
  }

  if (type == "logit_change") {
    y1 <- suppressWarnings(as.numeric(tbl$MaxLogitChangeElements))
    y2 <- suppressWarnings(as.numeric(tbl$MaxLogitChangeSteps))
    if (!any(is.finite(y1)) && !any(is.finite(y2))) {
      stop("No logit-change metrics available.")
    }
    if (isTRUE(draw)) {
      yr <- range(c(y1, y2), finite = TRUE)
      graphics::plot(
        x = it,
        y = y1,
        type = "b",
        pch = 16,
        col = pal["change_element"],
        ylim = yr,
        xlab = "Iteration",
        ylab = "Max absolute change",
        main = if (is.null(main)) "Iteration logit-change trajectory" else as.character(main[1])
      )
      graphics::lines(it, y2, type = "b", pch = 17, col = pal["change_step"])
      graphics::legend(
        "topright",
        legend = c("Elements", "Steps"),
        col = c(pal["change_element"], pal["change_step"]),
        pch = c(16, 17),
        lty = 1,
        bty = "n",
        cex = 0.85
      )
    }
    return(invisible(new_mfrm_plot_data(
      "iteration_report",
      list(plot = "logit_change", iteration = it, element = y1, step = y2)
    )))
  }

  vals <- suppressWarnings(as.numeric(tbl$Objective))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) stop("No objective values available.")
  it2 <- it[is.finite(suppressWarnings(as.numeric(tbl$Objective)))]
  if (isTRUE(draw)) {
    graphics::plot(
      x = it2,
      y = vals,
      type = "b",
      pch = 16,
      col = pal["objective"],
      xlab = "Iteration",
      ylab = "Objective (log-likelihood proxy)",
      main = if (is.null(main)) "Iteration objective trajectory" else as.character(main[1])
    )
  }
  invisible(new_mfrm_plot_data(
    "iteration_report",
    list(plot = "objective", iteration = it2, objective = vals)
  ))
}

draw_subset_connectivity_bundle <- function(x,
                                            type = c("subset_observations", "facet_levels", "coverage_matrix", "linking_matrix", "design_matrix"),
                                            draw = TRUE,
                                            main = NULL,
                                            palette = NULL,
                                            label_angle = 45,
                                            preset = c("standard", "publication", "compact")) {
  requested_type <- match.arg(tolower(as.character(type[1])), c("subset_observations", "facet_levels", "coverage_matrix", "linking_matrix", "design_matrix"))
  type <- requested_type
  if (type %in% c("linking_matrix", "design_matrix")) type <- "coverage_matrix"
  style <- resolve_plot_preset(preset)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      subset = style$accent_secondary,
      facet = style$accent_primary,
      low = "#f1f5f9",
      high = style$accent_tertiary
    )
  )
  summary_tbl <- as.data.frame(x$summary %||% data.frame(), stringsAsFactors = FALSE)
  listing_tbl <- as.data.frame(x$listing %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "subset_observations") {
    if (nrow(summary_tbl) == 0 || !all(c("Subset", "Observations") %in% names(summary_tbl))) {
      stop("Subset summary table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    vals <- suppressWarnings(as.numeric(summary_tbl$Observations))
    labels <- paste0("Subset ", as.character(summary_tbl$Subset))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["subset"],
        main = if (is.null(main)) "Observations by subset" else as.character(main[1]),
        ylab = "Observations",
        label_angle = label_angle,
        mar_bottom = 8.0
      )
    }
    return(invisible(new_mfrm_plot_data(
      "subset_connectivity",
      list(
        plot = "subset_observations",
        table = summary_tbl,
        title = main %||% "Observations by subset",
        subtitle = "Observation counts for connected subsets",
        legend = new_plot_legend("Observation count", "subset", "bar", pal["subset"]),
        reference_lines = new_reference_lines(),
        preset = style$name
      )
    )))
  }

  if (nrow(listing_tbl) == 0 || !all(c("Subset", "Facet", "LevelsN") %in% names(listing_tbl))) {
    stop("Subset facet-listing table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  }
  if (identical(type, "coverage_matrix")) {
    listing_tbl$Subset <- as.character(listing_tbl$Subset)
    cov_tbl <- listing_tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::mutate(
        MaxLevels = max(.data$LevelsN, na.rm = TRUE),
        CoverageRatio = dplyr::if_else(.data$MaxLevels > 0, .data$LevelsN / .data$MaxLevels, NA_real_)
      ) |>
      dplyr::ungroup()
    subset_summary <- summary_tbl |>
      dplyr::mutate(
        Subset = as.character(.data$Subset),
        ObservationPercent = dplyr::coalesce(.data$ObservationPercent, 0)
      ) |>
      dplyr::arrange(dplyr::desc(.data$Observations), .data$Subset)
    facet_summary <- cov_tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(
        MeanCoverage = mean(.data$CoverageRatio, na.rm = TRUE),
        CoveredSubsets = sum(.data$CoverageRatio > 0, na.rm = TRUE),
        CompleteSubsets = sum(.data$CoverageRatio >= 0.999, na.rm = TRUE),
        TotalSubsets = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(.data$MeanCoverage), dplyr::desc(.data$CoveredSubsets), .data$Facet)
    subset_order <- subset_summary$Subset
    facet_order <- facet_summary$Facet
    cov_tbl$Subset <- factor(cov_tbl$Subset, levels = subset_order)
    cov_tbl$Facet <- factor(cov_tbl$Facet, levels = facet_order)
    cov_wide <- tryCatch({
      tidyr::pivot_wider(
        cov_tbl[, c("Facet", "Subset", "CoverageRatio")],
        names_from = "Subset",
        values_from = "CoverageRatio",
        values_fill = list(CoverageRatio = 0)
      ) |>
        tibble::column_to_rownames("Facet") |>
        as.matrix()
    }, error = function(e) NULL)
    label_wide <- tryCatch({
      tidyr::pivot_wider(
        cov_tbl |>
          dplyr::mutate(CellLabel = paste0(.data$LevelsN, "/", .data$MaxLevels)),
        names_from = "Subset",
        values_from = "CellLabel",
        values_fill = list(CellLabel = "")
      ) |>
        tibble::column_to_rownames("Facet") |>
        as.matrix()
    }, error = function(e) NULL)
    if (is.null(cov_wide) || nrow(cov_wide) == 0) {
      stop("Coverage matrix could not be constructed from the subset listing table.", call. = FALSE)
    }
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(0.8, 2.7))
      graphics::par(mar = c(1.2, 7.5, 2.5, 2.2))
      obs_vals <- subset_summary$ObservationPercent
      mids <- graphics::barplot(
        height = obs_vals,
        col = grDevices::adjustcolor(pal["subset"], alpha.f = 0.72),
        border = NA,
        axes = FALSE,
        space = 0.2,
        ylim = c(0, max(c(obs_vals, 1), na.rm = TRUE) * 1.15),
        main = if (is.null(main)) "Linking design matrix" else as.character(main[1])
      )
      graphics::axis(2, las = 1)
      graphics::mtext("Observation share (%)", side = 2, line = 5.2)
      graphics::abline(h = pretty(obs_vals, n = 4), col = grDevices::adjustcolor(style$grid, alpha.f = 0.8), lty = 1)
      graphics::text(
        x = mids,
        y = obs_vals,
        labels = paste0(round(obs_vals, 1), "%"),
        pos = 3,
        cex = 0.7,
        col = style$foreground
      )

      graphics::par(mar = c(6.5, 7.5, 1.2, 4.4))
      cols <- grDevices::colorRampPalette(c(pal["low"], pal["high"]))(21)
      graphics::image(
        x = seq_len(ncol(cov_wide)),
        y = seq_len(nrow(cov_wide)),
        z = t(cov_wide[nrow(cov_wide):1, , drop = FALSE]),
        xaxt = "n",
        yaxt = "n",
        xlab = "Subset",
        ylab = "",
        col = cols,
        zlim = c(0, 1),
        main = ""
      )
      graphics::axis(1, at = seq_len(ncol(cov_wide)), labels = colnames(cov_wide), las = 2)
      graphics::axis(2, at = seq_len(nrow(cov_wide)), labels = rev(rownames(cov_wide)), las = 1)
      graphics::abline(v = seq(0.5, ncol(cov_wide) + 0.5, by = 1), col = grDevices::adjustcolor("white", alpha.f = 0.7))
      graphics::abline(h = seq(0.5, nrow(cov_wide) + 0.5, by = 1), col = grDevices::adjustcolor("white", alpha.f = 0.7))
      if (!is.null(label_wide)) {
        for (i in seq_len(nrow(cov_wide))) {
          for (j in seq_len(ncol(cov_wide))) {
            graphics::text(
              x = j,
              y = nrow(cov_wide) - i + 1,
              labels = label_wide[i, j],
              cex = 0.72,
              col = if (is.finite(cov_wide[i, j]) && cov_wide[i, j] >= 0.65) "white" else style$foreground
            )
          }
        }
      }
      graphics::mtext("Facet", side = 2, line = 6)
      graphics::axis(
        4,
        at = seq_len(nrow(cov_wide)),
        labels = rev(sprintf("%.0f%%", 100 * facet_summary$MeanCoverage)),
        las = 1,
        cex.axis = 0.82
      )
      graphics::mtext("Mean coverage", side = 4, line = 3)
    }
    return(invisible(new_mfrm_plot_data(
      "subset_connectivity",
      list(
        plot = "coverage_matrix",
        requested_type = requested_type,
        matrix = cov_wide,
        labels = label_wide,
        table = cov_tbl,
        facet_summary = as.data.frame(facet_summary, stringsAsFactors = FALSE),
        subset_summary = as.data.frame(subset_summary, stringsAsFactors = FALSE),
        title = main %||% "Linking design matrix",
        subtitle = "Subset observation share and facet-by-subset coverage",
        legend = new_plot_legend(
          label = c("Low facet coverage", "High facet coverage"),
          role = c("coverage", "coverage"),
          aesthetic = c("heatmap", "heatmap"),
          value = c(pal["low"], pal["high"])
        ),
        reference_lines = new_reference_lines(),
        preset = style$name
      )
    )))
  }
  vals <- suppressWarnings(as.numeric(listing_tbl$LevelsN))
  labels <- paste0("S", listing_tbl$Subset, ":", listing_tbl$Facet)
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["facet"],
      main = if (is.null(main)) "Facet levels by subset" else as.character(main[1]),
      ylab = "Levels",
      label_angle = label_angle,
      mar_bottom = 8.8
    )
  }
  invisible(new_mfrm_plot_data(
    "subset_connectivity",
    list(
      plot = "facet_levels",
      table = listing_tbl,
      title = main %||% "Facet levels by subset",
      subtitle = "Observed levels per facet within each connected subset",
      legend = new_plot_legend("Observed levels", "facet", "bar", pal["facet"]),
      reference_lines = new_reference_lines(),
      preset = style$name
    )
  ))
}

draw_facet_statistics_bundle <- function(x,
                                         type = c("means", "sds", "ranges"),
                                         metric = NULL,
                                         draw = TRUE,
                                         main = NULL,
                                         palette = NULL,
                                         label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("means", "sds", "ranges"))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      mean = "#2b8cbe",
      sd = "#756bb1",
      range = "#9ecae1"
    )
  )
  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0 || !all(c("Metric", "Facet", "Mean", "SD", "Min", "Max") %in% names(tbl))) {
    stop("Facet-statistics table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  }
  metrics <- unique(as.character(tbl$Metric))
  if (is.null(metric)) metric <- metrics[1]
  metric <- as.character(metric[1])
  if (!metric %in% metrics) {
    stop("Requested `metric` not found. Available: ", paste(metrics, collapse = ", "))
  }
  sub <- tbl[as.character(tbl$Metric) == metric, , drop = FALSE]
  sub <- sub[order(as.character(sub$Facet)), , drop = FALSE]

  if (type == "means") {
    vals <- suppressWarnings(as.numeric(sub$Mean))
    labels <- as.character(sub$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["mean"],
        main = if (is.null(main)) paste0("Facet means (", metric, ")") else as.character(main[1]),
        ylab = "Mean",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "facet_statistics",
      list(plot = "means", metric = metric, table = sub)
    )))
  }

  if (type == "sds") {
    vals <- suppressWarnings(as.numeric(sub$SD))
    labels <- as.character(sub$Facet)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["sd"],
        main = if (is.null(main)) paste0("Facet SDs (", metric, ")") else as.character(main[1]),
        ylab = "SD",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "facet_statistics",
      list(plot = "sds", metric = metric, table = sub)
    )))
  }

  y <- seq_len(nrow(sub))
  mn <- suppressWarnings(as.numeric(sub$Min))
  mx <- suppressWarnings(as.numeric(sub$Max))
  md <- suppressWarnings(as.numeric(sub$Mean))
  if (isTRUE(draw)) {
    xr <- range(c(mn, mx), finite = TRUE)
    graphics::plot(
      x = xr,
      y = c(1, nrow(sub)),
      type = "n",
      yaxt = "n",
      xlab = metric,
      ylab = "",
      main = if (is.null(main)) paste0("Facet ranges (", metric, ")") else as.character(main[1])
    )
    graphics::segments(x0 = mn, y0 = y, x1 = mx, y1 = y, col = pal["range"], lwd = 2)
    graphics::points(md, y, pch = 16, col = pal["mean"])
    graphics::axis(side = 2, at = y, labels = as.character(sub$Facet), las = 2, cex.axis = 0.8)
  }
  invisible(new_mfrm_plot_data(
    "facet_statistics",
    list(plot = "ranges", metric = metric, table = sub)
  ))
}

draw_residual_pca_bundle <- function(x,
                                     type = c("overall_scree", "facet_scree", "overall_loadings", "facet_loadings"),
                                     facet = NULL,
                                     component = 1L,
                                     top_n = 20L,
                                     draw = TRUE) {
  type <- match.arg(tolower(as.character(type[1])), c("overall_scree", "facet_scree", "overall_loadings", "facet_loadings"))
  if (type == "overall_scree") {
    return(invisible(plot_residual_pca(
      x,
      mode = "overall",
      plot_type = "scree",
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  if (type == "facet_scree") {
    return(invisible(plot_residual_pca(
      x,
      mode = "facet",
      facet = facet,
      plot_type = "scree",
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  if (type == "overall_loadings") {
    return(invisible(plot_residual_pca(
      x,
      mode = "overall",
      plot_type = "loadings",
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  invisible(plot_residual_pca(
    x,
    mode = "facet",
    facet = facet,
    plot_type = "loadings",
    component = component,
    top_n = top_n,
    draw = draw
  ))
}

draw_parity_bundle <- function(x,
                               type = c("column_coverage", "table_coverage", "metric_status", "metric_by_table"),
                               top_n = 40,
                               draw = TRUE,
                               main = NULL,
                               palette = NULL,
                               label_angle = 45) {
  type <- match.arg(tolower(as.character(type[1])), c("column_coverage", "table_coverage", "metric_status", "metric_by_table"))
  top_n <- max(1L, as.integer(top_n))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      pass = "#31a354",
      fail = "#cb181d",
      missing = "#969696",
      coverage = "#3182bd",
      metric = "#756bb1"
    )
  )

  column_audit <- as.data.frame(x$column_audit %||% data.frame(), stringsAsFactors = FALSE)
  column_summary <- as.data.frame(x$column_summary %||% data.frame(), stringsAsFactors = FALSE)
  metric_audit <- as.data.frame(x$metric_audit %||% data.frame(), stringsAsFactors = FALSE)
  metric_by_table <- as.data.frame(x$metric_by_table %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "column_coverage") {
    if (nrow(column_audit) == 0 || !all(c("table_id", "component", "coverage", "available", "full_match") %in% names(column_audit))) {
      stop("Column-audit table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    tbl <- column_audit
    tbl$coverage <- suppressWarnings(as.numeric(tbl$coverage))
    tbl <- tbl |>
      dplyr::arrange(.data$coverage, .data$table_id, .data$component)
    if (nrow(tbl) > top_n) tbl <- tbl |> dplyr::slice_head(n = top_n)
    vals <- ifelse(is.finite(tbl$coverage), tbl$coverage, 0)
    labels <- paste0(tbl$table_id, ":", tbl$component)
    cols <- ifelse(!tbl$available, pal["missing"], ifelse(tbl$full_match, pal["pass"], pal["fail"]))
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = cols,
        main = if (is.null(main)) "Column contract coverage (lowest first)" else as.character(main[1]),
        ylab = "Coverage",
        label_angle = label_angle,
        mar_bottom = 9.2
      )
      graphics::abline(h = 1, lty = 3, col = "#999999")
    }
    return(invisible(new_mfrm_plot_data(
      "parity_report",
      list(plot = "column_coverage", table = tbl, labels = labels)
    )))
  }

  if (type == "table_coverage") {
    if (nrow(column_summary) == 0 || !all(c("table_id", "MeanCoverage") %in% names(column_summary))) {
      stop("Column-summary table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    tbl <- column_summary |>
      dplyr::arrange(.data$table_id)
    vals <- suppressWarnings(as.numeric(tbl$MeanCoverage))
    vals[!is.finite(vals)] <- 0
    labels <- as.character(tbl$table_id)
    if (isTRUE(draw)) {
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["coverage"],
        main = if (is.null(main)) "Mean column coverage by table" else as.character(main[1]),
        ylab = "Mean coverage",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
      graphics::abline(h = 1, lty = 3, col = "#999999")
    }
    return(invisible(new_mfrm_plot_data(
      "parity_report",
      list(plot = "table_coverage", table = tbl)
    )))
  }

  if (type == "metric_status") {
    if (nrow(metric_audit) == 0 || !"Pass" %in% names(metric_audit)) {
      stop("Metric-audit table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    status <- ifelse(is.na(metric_audit$Pass), "Not evaluated", ifelse(metric_audit$Pass %in% TRUE, "Pass", "Fail"))
    cnt <- table(factor(status, levels = c("Pass", "Fail", "Not evaluated")))
    vals <- as.numeric(cnt)
    labels <- names(cnt)
    cols <- c(pal["pass"], pal["fail"], pal["missing"])
    if (isTRUE(draw)) {
      graphics::barplot(
        height = vals,
        names.arg = labels,
        col = cols,
        las = 2,
        ylab = "Checks",
        main = if (is.null(main)) "Metric-check status counts" else as.character(main[1])
      )
    }
    return(invisible(new_mfrm_plot_data(
      "parity_report",
      list(plot = "metric_status", table = data.frame(Status = labels, Checks = vals, stringsAsFactors = FALSE))
    )))
  }

  if (nrow(metric_by_table) == 0 || !all(c("Table", "PassRate") %in% names(metric_by_table))) {
    stop("Metric-by-table summary is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  }
  tbl <- metric_by_table |>
    dplyr::arrange(.data$Table)
  vals <- suppressWarnings(as.numeric(tbl$PassRate))
  vals[!is.finite(vals)] <- 0
  labels <- as.character(tbl$Table)
  if (isTRUE(draw)) {
    barplot_rot45(
      height = vals,
      labels = labels,
      col = pal["metric"],
      main = if (is.null(main)) "Metric pass rate by table" else as.character(main[1]),
      ylab = "Pass rate",
      label_angle = label_angle,
      mar_bottom = 7.8
    )
    graphics::abline(h = 1, lty = 3, col = "#999999")
  }
  invisible(new_mfrm_plot_data(
    "parity_report",
    list(plot = "metric_by_table", table = tbl)
  ))
}

plot_bias_count_bundle <- function(x,
                                   plot_type = c("cell_counts", "lowcount_by_facet"),
                                   top_n = 40,
                                   draw = TRUE,
                                   main = NULL,
                                   palette = NULL,
                                   label_angle = 45) {
  plot_type <- match.arg(tolower(as.character(plot_type[1])), c("cell_counts", "lowcount_by_facet"))
  top_n <- max(1L, as.integer(top_n))
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      count = "#2b8cbe",
      low = "#cb181d",
      rate = "#756bb1",
      grid = "#ececec"
    )
  )

  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  if ("Observd Count" %in% names(tbl) && !"Count" %in% names(tbl)) {
    tbl$Count <- suppressWarnings(as.numeric(tbl$`Observd Count`))
  }
  if (!"Count" %in% names(tbl)) {
    stop("Bias-count table does not include a count column.")
  }
  if (!"LowCountFlag" %in% names(tbl)) {
    tbl$LowCountFlag <- FALSE
  }
  tbl$LowCountFlag <- as.logical(tbl$LowCountFlag)

  if (plot_type == "cell_counts") {
    tbl <- tbl[is.finite(suppressWarnings(as.numeric(tbl$Count))), , drop = FALSE]
    if (nrow(tbl) == 0) stop("No finite count rows available.")
    tbl$Count <- suppressWarnings(as.numeric(tbl$Count))
    ord <- order(tbl$Count, decreasing = TRUE, na.last = NA)
    use <- ord[seq_len(min(length(ord), top_n))]
    tbl <- tbl[use, , drop = FALSE]

    facet_cols <- names(x$by_facet %||% list())
    facet_cols <- facet_cols[facet_cols %in% names(tbl)]
    if (length(facet_cols) == 0) {
      facet_cols <- names(tbl)[vapply(tbl, is.character, logical(1))]
      facet_cols <- setdiff(facet_cols, c("Count", "LowCountFlag"))
      facet_cols <- facet_cols[seq_len(min(2L, length(facet_cols)))]
    }
    labels <- if (length(facet_cols) > 0) {
      apply(tbl[, facet_cols, drop = FALSE], 1, paste, collapse = " | ")
    } else {
      paste0("Cell ", seq_len(nrow(tbl)))
    }

    if (isTRUE(draw)) {
      barplot_rot45(
        height = tbl$Count,
        labels = labels,
        col = ifelse(tbl$LowCountFlag %in% TRUE, pal["low"], pal["count"]),
        main = if (is.null(main)) "Bias cell counts" else as.character(main[1]),
        ylab = "Observed count",
        label_angle = label_angle,
        mar_bottom = 9.0
      )
    }
    return(invisible(new_mfrm_plot_data(
      "bias_count",
      list(plot = "cell_counts", table = tbl, labels = labels)
    )))
  }

  by_facet <- x$by_facet %||% list()
  rate_tbl <- lapply(names(by_facet), function(facet_nm) {
    df <- as.data.frame(by_facet[[facet_nm]], stringsAsFactors = FALSE)
    if (!all(c("Level", "Cells", "LowCountCells") %in% names(df))) return(NULL)
    data.frame(
      Facet = facet_nm,
      Level = as.character(df$Level),
      Cells = suppressWarnings(as.numeric(df$Cells)),
      LowCountCells = suppressWarnings(as.numeric(df$LowCountCells)),
      LowCountRate = ifelse(
        suppressWarnings(as.numeric(df$Cells)) > 0,
        suppressWarnings(as.numeric(df$LowCountCells)) / suppressWarnings(as.numeric(df$Cells)),
        NA_real_
      ),
      stringsAsFactors = FALSE
    )
  })
  rate_tbl <- rate_tbl[!vapply(rate_tbl, is.null, logical(1))]
  if (length(rate_tbl) == 0) {
    stop("No by-facet low-count summary available.")
  }
  rate_tbl <- dplyr::bind_rows(rate_tbl)
  rate_tbl <- rate_tbl[is.finite(rate_tbl$LowCountRate), , drop = FALSE]
  if (nrow(rate_tbl) == 0) {
    stop("No finite low-count rates available.")
  }
  rate_tbl <- rate_tbl |>
    dplyr::arrange(dplyr::desc(.data$LowCountRate), dplyr::desc(.data$LowCountCells), .data$Facet, .data$Level) |>
    dplyr::slice_head(n = top_n)
  labels <- paste0(rate_tbl$Facet, ":", rate_tbl$Level)

  if (isTRUE(draw)) {
    barplot_rot45(
      height = rate_tbl$LowCountRate,
      labels = labels,
      col = pal["rate"],
      main = if (is.null(main)) "Low-count rate by facet level" else as.character(main[1]),
      ylab = "Low-count rate",
      label_angle = label_angle,
      mar_bottom = 9.0
    )
    graphics::abline(h = 0, col = pal["grid"], lty = 1)
  }
  invisible(new_mfrm_plot_data(
    "bias_count",
    list(plot = "lowcount_by_facet", table = rate_tbl, labels = labels)
  ))
}

plot_fixed_reports_bundle <- function(x,
                                      plot_type = c("contrast", "pvalue"),
                                      top_n = 30,
                                      draw = TRUE,
                                      main = NULL,
                                      palette = NULL,
                                      label_angle = 45) {
  plot_type <- match.arg(tolower(as.character(plot_type[1])), c("contrast", "pvalue"))
  top_n <- max(1L, as.integer(top_n))
  pair_tbl <- as.data.frame(x$pairwise_table %||% data.frame(), stringsAsFactors = FALSE)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      pos = "#1b9e77",
      neg = "#d95f02",
      hist = "#756bb1"
    )
  )
  no_data_plot <- function(message) {
    if (isTRUE(draw)) {
      graphics::plot.new()
      graphics::title(main = if (is.null(main)) {
        if (identical(plot_type, "contrast")) "Pairwise contrasts" else "Pairwise p-value distribution"
      } else {
        as.character(main[1])
      })
      graphics::text(0.5, 0.5, "No data")
    }
    payload <- list(plot = plot_type, message = message)
    if (identical(plot_type, "contrast")) {
      payload$table <- pair_tbl[0, , drop = FALSE]
      payload$labels <- character(0)
    } else {
      payload$p_values <- numeric(0)
    }
    invisible(new_mfrm_plot_data("fixed_reports", payload))
  }
  if (nrow(pair_tbl) == 0) {
    return(no_data_plot("Pairwise table is empty; no plot is available."))
  }

  if (plot_type == "contrast") {
    if (!"Contrast" %in% names(pair_tbl)) {
      stop("Pairwise table does not include `Contrast`.")
    }
    pair_tbl$Contrast <- suppressWarnings(as.numeric(pair_tbl$Contrast))
    pair_tbl <- pair_tbl[is.finite(pair_tbl$Contrast), , drop = FALSE]
    if (nrow(pair_tbl) == 0) {
      return(no_data_plot("No finite contrast values available."))
    }
    pair_tbl <- pair_tbl |>
      dplyr::mutate(.abs = abs(.data$Contrast)) |>
      dplyr::arrange(dplyr::desc(.data$.abs)) |>
      dplyr::slice_head(n = top_n)
    labels <- if (all(c("Target", "Context1", "Context2") %in% names(pair_tbl))) {
      paste0(pair_tbl$Target, ": ", pair_tbl$Context1, " vs ", pair_tbl$Context2)
    } else {
      paste0("Pair ", seq_len(nrow(pair_tbl)))
    }
    if (isTRUE(draw)) {
      barplot_rot45(
        height = pair_tbl$Contrast,
        labels = labels,
        col = ifelse(pair_tbl$Contrast >= 0, pal["pos"], pal["neg"]),
        main = if (is.null(main)) "Pairwise contrasts" else as.character(main[1]),
        ylab = "Contrast (logit)",
        label_angle = label_angle,
        mar_bottom = 9.2
      )
      graphics::abline(h = 0, lty = 2, col = "gray50")
    }
    return(invisible(new_mfrm_plot_data(
      "fixed_reports",
      list(plot = "contrast", table = pair_tbl, labels = labels)
    )))
  }

  p_col <- if ("Prob." %in% names(pair_tbl)) "Prob." else if ("p.value" %in% names(pair_tbl)) "p.value" else NA_character_
  if (is.na(p_col)) {
    stop("Pairwise table does not include p-value column (`Prob.` or `p.value`).")
  }
  p_vals <- suppressWarnings(as.numeric(pair_tbl[[p_col]]))
  p_vals <- p_vals[is.finite(p_vals)]
  if (length(p_vals) == 0) {
    return(no_data_plot("No finite p-values available."))
  }
  if (isTRUE(draw)) {
    graphics::hist(
      x = p_vals,
      breaks = "FD",
      col = pal["hist"],
      border = "white",
      main = if (is.null(main)) "Pairwise p-value distribution" else as.character(main[1]),
      xlab = "p-value",
      ylab = "Count"
    )
    graphics::abline(v = 0.05, lty = 2, col = "gray45")
  }
  invisible(new_mfrm_plot_data(
    "fixed_reports",
    list(plot = "pvalue", p_values = p_vals)
  ))
}

plot_visual_summaries_bundle <- function(x,
                                         plot_type = c("comparison", "warning_counts", "summary_counts"),
                                         draw = TRUE,
                                         main = NULL,
                                         palette = NULL,
                                         label_angle = 45) {
  plot_type <- match.arg(tolower(as.character(plot_type[1])), c("comparison", "warning_counts", "summary_counts"))
  warning_counts <- as.data.frame(x$warning_counts %||% data.frame(), stringsAsFactors = FALSE)
  summary_counts <- as.data.frame(x$summary_counts %||% data.frame(), stringsAsFactors = FALSE)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      warning = "#cb181d",
      summary = "#2b8cbe",
      single = "#756bb1"
    )
  )

  if (plot_type == "warning_counts" || plot_type == "summary_counts") {
    tbl <- if (plot_type == "warning_counts") warning_counts else summary_counts
    if (nrow(tbl) == 0 || !all(c("Visual", "Messages") %in% names(tbl))) {
      stop("Requested count table is empty.")
    }
    if (isTRUE(draw)) {
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl$Messages)),
        labels = as.character(tbl$Visual),
        col = pal["single"],
        main = if (is.null(main)) {
          if (plot_type == "warning_counts") "Warning message counts by visual" else "Summary message counts by visual"
        } else {
          as.character(main[1])
        },
        ylab = "Messages",
        label_angle = label_angle,
        mar_bottom = 8.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "visual_summaries",
      list(plot = plot_type, table = tbl)
    )))
  }

  vis <- sort(unique(c(as.character(warning_counts$Visual), as.character(summary_counts$Visual))))
  if (length(vis) == 0) {
    stop("No warning/summary counts available.")
  }
  warn <- stats::setNames(rep(0, length(vis)), vis)
  summ <- stats::setNames(rep(0, length(vis)), vis)
  if (nrow(warning_counts) > 0) {
    warn[as.character(warning_counts$Visual)] <- suppressWarnings(as.numeric(warning_counts$Messages))
  }
  if (nrow(summary_counts) > 0) {
    summ[as.character(summary_counts$Visual)] <- suppressWarnings(as.numeric(summary_counts$Messages))
  }
  mat <- rbind(warn, summ)
  rownames(mat) <- c("Warning", "Summary")
  if (isTRUE(draw)) {
    old_mar <- graphics::par("mar")
    on.exit(graphics::par(mar = old_mar), add = TRUE)
    mar <- old_mar
    mar[1] <- max(mar[1], 9.0)
    graphics::par(mar = mar)
    mids <- graphics::barplot(
      height = mat,
      beside = TRUE,
      names.arg = FALSE,
      col = c(pal["warning"], pal["summary"]),
      ylab = "Messages",
      main = if (is.null(main)) "Warning vs summary counts by visual" else as.character(main[1]),
      border = "white"
    )
    centers <- vapply(split(as.numeric(mids), rep(seq_along(vis), each = 2L)), mean, numeric(1))
    draw_rotated_x_labels(
      at = centers,
      labels = vis,
      srt = label_angle,
      cex = 0.82,
      line_offset = 0.085
    )
    graphics::legend(
      "topright",
      legend = c("Warning", "Summary"),
      fill = c(pal["warning"], pal["summary"]),
      bty = "n",
      cex = 0.85
    )
  }
  invisible(new_mfrm_plot_data(
    "visual_summaries",
    list(plot = "comparison", matrix = mat, visuals = vis)
  ))
}

#' Plot report/table bundles with base R defaults
#'
#' @param x A bundle object returned by mfrmr table/report helpers.
#' @param y Reserved for generic compatibility.
#' @param type Optional plot type. Available values depend on bundle class.
#' @param ... Additional arguments forwarded to class-specific plotters.
#'
#' @details
#' `plot()` dispatches by bundle class:
#' - `mfrm_unexpected` -> [plot_unexpected()]
#' - `mfrm_fair_average` -> [plot_fair_average()]
#' - `mfrm_displacement` -> [plot_displacement()]
#' - `mfrm_interrater` -> [plot_interrater_agreement()]
#' - `mfrm_facets_chisq` -> [plot_facets_chisq()]
#' - `mfrm_bias_interaction` -> [plot_bias_interaction()]
#' - `mfrm_bias_count` -> bias-count plots (cell counts / low-count rates)
#' - `mfrm_fixed_reports` -> pairwise-contrast diagnostics
#' - `mfrm_visual_summaries` -> warning/summary message count plots
#' - `mfrm_category_structure` -> default base-R category plots
#' - `mfrm_category_curves` -> default ogive/CCC plots
#' - `mfrm_rating_scale` -> category-counts/threshold plots
#' - `mfrm_measurable` -> measurable-data coverage/count plots
#' - `mfrm_unexpected_after_bias` -> post-bias unexpected-response plots
#' - `mfrm_output_bundle` -> graph/score output-file diagnostics
#' - `mfrm_residual_pca` -> residual PCA scree/loadings via [plot_residual_pca()]
#' - `mfrm_specifications` -> facet/anchor/convergence plots
#' - `mfrm_data_quality` -> row-audit/category/missing-row plots
#' - `mfrm_iteration_report` -> replayed-iteration trajectories
#' - `mfrm_subset_connectivity` -> subset-observation/connectivity plots
#' - `mfrm_facet_statistics` -> facet statistic profile plots
#'
#' If a class is outside these families, use dedicated plotting helpers
#' or custom base R graphics on component tables.
#'
#' @section Interpreting output:
#' The returned object is plotting data (`mfrm_plot_data`) that captures
#' the selected route and payload; set `draw = TRUE` for immediate base graphics.
#'
#' @section Typical workflow:
#' 1. Create bundle output (e.g., `unexpected_response_table()`).
#' 2. Inspect routing with `summary(bundle)` if needed.
#' 3. Call `plot(bundle, type = ..., draw = FALSE)` to obtain reusable plot data.
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso `summary()`, [plot_unexpected()], [plot_fair_average()], [plot_displacement()]
#' @examples
#' \donttest{
#' toy_full <- load_mfrmr_data("example_core")
#' toy_people <- unique(toy_full$Person)[1:12]
#' toy <- toy_full[toy_full$Person %in% toy_people, , drop = FALSE]
#' fit <- suppressWarnings(
#'   fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 10)
#' )
#' t4 <- unexpected_response_table(fit, abs_z_min = 1.5, prob_max = 0.4, top_n = 5)
#' p <- plot(t4, draw = FALSE)
#' vis <- build_visual_summaries(fit, diagnose_mfrm(fit, residual_pca = "none"))
#' p_vis <- plot(vis, type = "comparison", draw = FALSE)
#' spec <- specifications_report(fit)
#' p_spec <- plot(spec, type = "facet_elements", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     t4,
#'     type = "severity",
#'     draw = TRUE,
#'     main = "Unexpected Response Severity (Customized)",
#'     palette = c(higher = "#d95f02", lower = "#1b9e77", bar = "#2b8cbe"),
#'     label_angle = 45
#'   )
#'   plot(
#'     vis,
#'     type = "comparison",
#'     draw = TRUE,
#'     main = "Warning vs Summary Counts (Customized)",
#'     palette = c(warning = "#cb181d", summary = "#3182bd"),
#'     label_angle = 45
#'   )
#' }
#' }
#' @export
plot.mfrm_bundle <- function(x, y = NULL, type = NULL, ...) {
  dots <- list(...)

  if (inherits(x, "mfrm_unexpected")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_unexpected, args))
  }
  if (inherits(x, "mfrm_fair_average")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_fair_average, args))
  }
  if (inherits(x, "mfrm_displacement")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_displacement, args))
  }
  if (inherits(x, "mfrm_interrater")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_interrater_agreement, args))
  }
  if (inherits(x, "mfrm_facets_chisq")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_facets_chisq, args))
  }
  if (inherits(x, "mfrm_bias_interaction")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot <- type
    return(do.call(plot_bias_interaction, args))
  }
  if (inherits(x, "mfrm_bias_count")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_bias_count_bundle, args))
  }
  if (inherits(x, "mfrm_fixed_reports")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_fixed_reports_bundle, args))
  }
  if (inherits(x, "mfrm_visual_summaries")) {
    args <- c(list(x = x), dots)
    if (!is.null(type)) args$plot_type <- type
    return(do.call(plot_visual_summaries_bundle, args))
  }
  if (inherits(x, "mfrm_parity_report")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "column_coverage" else as.character(type[1])
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 40L
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_parity_bundle(
      x,
      type = ptype,
      top_n = top_n,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_category_structure")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "counts" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_category_structure_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_category_curves")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "ogive" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    return(invisible(draw_category_curves_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette
    )))
  }
  if (inherits(x, "mfrm_rating_scale")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "counts" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_rating_scale_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_measurable")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "facet_coverage" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_measurable_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_unexpected_after_bias")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "scatter" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 40L
    return(invisible(draw_unexpected_after_bias_bundle(
      x,
      type = ptype,
      top_n = top_n,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_output_bundle")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "graph_expected" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    return(invisible(draw_output_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette
    )))
  }
  if (inherits(x, "mfrm_residual_pca")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "overall_scree" else as.character(type[1])
    facet <- dots$facet %||% NULL
    component <- if ("component" %in% names(dots)) dots$component else 1L
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 20L
    return(invisible(draw_residual_pca_bundle(
      x,
      type = ptype,
      facet = facet,
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  if (inherits(x, "mfrm_specifications")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "facet_elements" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_specifications_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_data_quality")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "row_audit" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_data_quality_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }
  if (inherits(x, "mfrm_iteration_report")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "residual" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    return(invisible(draw_iteration_report_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette
    )))
  }
  if (inherits(x, "mfrm_subset_connectivity")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "subset_observations" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    preset <- dots$preset %||% "standard"
    return(invisible(draw_subset_connectivity_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle,
      preset = preset
    )))
  }
  if (inherits(x, "mfrm_facet_statistics")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "means" else as.character(type[1])
    metric <- dots$metric %||% NULL
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_facet_statistics_bundle(
      x,
      type = ptype,
      metric = metric,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle
    )))
  }

  stop(
    "No default plot method for class `", class(x)[1], "`.\n",
    "Use a dedicated plot helper (for example, `plot_unexpected()`, `plot_fair_average()`, or `plot_bias_interaction()`)."
  )
}

#' Summarize an `mfrm_diagnostics` object in a user-friendly format
#'
#' @param object Output from [diagnose_mfrm()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of highest-absolute-Z fit rows to keep.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method returns a compact diagnostics summary designed for quick review:
#' - design overview (observations, persons, facets, categories, subsets)
#' - global fit statistics
#' - approximate reliability/separation by facet
#' - top facet/person fit rows by absolute ZSTD
#' - counts of flagged diagnostics (unexpected, displacement, interactions)
#'
#' @section Interpreting output:
#' - `overview`: analysis scale, subset count, and residual-PCA mode.
#' - `overall_fit`: global fit indices.
#' - `reliability`: facet separation/reliability block, including model and
#'   real bounds when available.
#' - `top_fit`: highest `|ZSTD|` elements for immediate inspection.
#' - `flags`: compact counts for key warning domains.
#'
#' @section Typical workflow:
#' 1. Run diagnostics with [diagnose_mfrm()].
#' 2. Review `summary(diag)` for major warnings.
#' 3. Follow up with dedicated tables/plots for flagged domains.
#'
#' @return An object of class `summary.mfrm_diagnostics` with:
#' - `overview`: design-level counts and residual-PCA mode
#' - `overall_fit`: global fit block
#' - `reliability`: facet-level separation/reliability summary
#' - `top_fit`: top `|ZSTD|` rows
#' - `flags`: compact flag counts for major diagnostics
#' - `notes`: short interpretation notes
#' @seealso [diagnose_mfrm()], [summary.mfrm_fit()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' summary(diag)
#' @export
summary.mfrm_diagnostics <- function(object, digits = 3, top_n = 10, ...) {
  if (!is.list(object) || is.null(object$obs)) {
    stop("`object` must be output from diagnose_mfrm().")
  }

  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  obs_tbl <- tibble::as_tibble(object$obs)
  fit_tbl <- tibble::as_tibble(object$fit %||% tibble::tibble())
  reliability_tbl <- tibble::as_tibble(object$reliability %||% tibble::tibble())
  precision_profile_tbl <- tibble::as_tibble(object$precision_profile %||% tibble::tibble())
  precision_audit_tbl <- tibble::as_tibble(object$precision_audit %||% tibble::tibble())
  overall_fit <- tibble::as_tibble(object$overall_fit %||% tibble::tibble())
  subset_summary <- tibble::as_tibble(object$subsets$summary %||% tibble::tibble())
  approximation_tbl <- tibble::as_tibble(object$approximation_notes %||% tibble::tibble())

  n_obs <- nrow(obs_tbl)
  n_person <- if ("Person" %in% names(obs_tbl)) dplyr::n_distinct(obs_tbl$Person) else NA_integer_
  n_cat <- if ("Observed" %in% names(obs_tbl)) dplyr::n_distinct(obs_tbl$Observed) else NA_integer_
  n_subsets <- if ("Subset" %in% names(subset_summary)) dplyr::n_distinct(subset_summary$Subset) else 0L

  overview <- tibble::tibble(
    Observations = n_obs,
    Persons = n_person,
    Facets = length(object$facet_names %||% character(0)),
    Categories = n_cat,
    Subsets = n_subsets,
    ResidualPCA = as.character(object$residual_pca_mode %||% "none"),
    Method = as.character(precision_profile_tbl$Method[1] %||% NA_character_),
    PrecisionTier = as.character(precision_profile_tbl$PrecisionTier[1] %||% NA_character_)
  )

  reliability_overview <- tibble::tibble()
  keep_rel <- c(
    "Facet", "Levels",
    "Separation", "Strata", "Reliability",
    "RealSeparation", "RealStrata", "RealReliability",
    "MeanInfit", "MeanOutfit"
  )
  if (nrow(reliability_tbl) > 0) {
    keep <- intersect(keep_rel, names(reliability_tbl))
    reliability_overview <- reliability_tbl |>
      dplyr::select(dplyr::all_of(keep)) |>
      dplyr::arrange(.data$Facet)
  }

  top_fit <- tibble::tibble()
  fit_need <- c("Facet", "Level", "Infit", "Outfit", "InfitZSTD", "OutfitZSTD")
  if (nrow(fit_tbl) > 0 && all(fit_need %in% names(fit_tbl))) {
    top_fit <- fit_tbl |>
      dplyr::mutate(
        AbsZ = pmax(abs(.data$InfitZSTD), abs(.data$OutfitZSTD), na.rm = TRUE)
      ) |>
      dplyr::arrange(dplyr::desc(.data$AbsZ)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select("Facet", "Level", "Infit", "Outfit", "InfitZSTD", "OutfitZSTD", "AbsZ")
  }

  unexpected_n <- suppressWarnings(as.integer(object$unexpected$summary$UnexpectedN[1] %||% NA_integer_))
  displacement_flagged <- suppressWarnings(as.integer(object$displacement$summary$FlaggedLevels[1] %||% NA_integer_))
  interaction_n <- if (!is.null(object$interactions)) nrow(object$interactions) else NA_integer_
  interrater_pairs <- suppressWarnings(as.integer(object$interrater$summary$Pairs[1] %||% NA_integer_))

  flags <- tibble::tibble(
    Metric = c(
      "Unexpected responses",
      "Flagged displacement levels",
      "Interaction rows",
      "Inter-rater pairs"
    ),
    Count = c(unexpected_n, displacement_flagged, interaction_n, interrater_pairs)
  )

  notes <- character(0)
  if (nrow(precision_profile_tbl) > 0 && identical(as.character(precision_profile_tbl$PrecisionTier[1]), "exploratory")) {
    notes <- c(notes, "Precision outputs are exploratory for this run; prefer MML for formal SE, CI, and reliability reporting.")
  }
  if (nrow(precision_profile_tbl) > 0 && identical(as.character(precision_profile_tbl$PrecisionTier[1]), "hybrid")) {
    notes <- c(notes, "Precision outputs are hybrid for this run; inspect levels that fell back to observation-table information before treating the run as fully inferential.")
  }
  if (isTRUE(n_subsets > 1L)) {
    notes <- c(notes, "Multiple disconnected subsets were detected.")
  }
  if (isTRUE(!is.na(unexpected_n) && unexpected_n > 0L)) {
    notes <- c(notes, "Unexpected responses were flagged under current thresholds.")
  }
  if (nrow(approximation_tbl) > 0) {
    notes <- c(
      notes,
      "SE/ModelSE, CI, and reliability conventions depend on the estimation path; see diagnostics$approximation_notes for MML-vs-JML details."
    )
  }
  if (nrow(precision_audit_tbl) > 0 && "Status" %in% names(precision_audit_tbl)) {
    flagged_checks <- precision_audit_tbl |>
      dplyr::filter(.data$Status %in% c("review", "warn"))
    if (nrow(flagged_checks) > 0) {
      notes <- c(notes, paste0("Precision audit flagged ", nrow(flagged_checks), " review/warn checks."))
    }
  }
  if (length(notes) == 0) {
    notes <- "No immediate warnings from diagnostics summary."
  }

  out <- list(
    overview = overview,
    overall_fit = overall_fit,
    precision_profile = precision_profile_tbl,
    precision_audit = precision_audit_tbl,
    reliability = reliability_overview,
    top_fit = top_fit,
    flags = flags,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_diagnostics"
  out
}

#' @export
print.summary.mfrm_diagnostics <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("Many-Facet Rasch Diagnostics Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    ov <- round_numeric_df(as.data.frame(x$overview), digits = digits)[1, , drop = FALSE]
    cat(sprintf(
      "  Observations: %s | Persons: %s | Facets: %s | Categories: %s | Subsets: %s\n",
      ov$Observations, ov$Persons, ov$Facets, ov$Categories, ov$Subsets
    ))
    cat(sprintf("  Residual PCA mode: %s\n", ov$ResidualPCA))
    if ("Method" %in% names(ov) && "PrecisionTier" %in% names(ov)) {
      cat(sprintf("  Method: %s | Precision tier: %s\n", ov$Method, ov$PrecisionTier))
    }
  }

  if (!is.null(x$overall_fit) && nrow(x$overall_fit) > 0) {
    cat("\nOverall fit\n")
    print(round_numeric_df(as.data.frame(x$overall_fit), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$precision_profile) && nrow(x$precision_profile) > 0) {
    cat("\nPrecision basis\n")
    print(round_numeric_df(as.data.frame(x$precision_profile), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$precision_audit) && nrow(x$precision_audit) > 0) {
    flagged_checks <- as.data.frame(x$precision_audit)
    if ("Status" %in% names(flagged_checks)) {
      flagged_checks <- flagged_checks[flagged_checks$Status %in% c("review", "warn"), , drop = FALSE]
    }
    if (nrow(flagged_checks) > 0) {
      cat("\nPrecision audit checks to review\n")
      print(flagged_checks, row.names = FALSE)
    }
  }
  if (!is.null(x$reliability) && nrow(x$reliability) > 0) {
    cat("\nFacet precision and spread\n")
    print(round_numeric_df(as.data.frame(x$reliability), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$top_fit) && nrow(x$top_fit) > 0) {
    cat("\nLargest |ZSTD| rows\n")
    print(round_numeric_df(as.data.frame(x$top_fit), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$flags) && nrow(x$flags) > 0) {
    cat("\nFlag counts\n")
    print(as.data.frame(x$flags), row.names = FALSE)
  }

  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize an `mfrm_bias` object in a user-friendly format
#'
#' @param object Output from [estimate_bias()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of strongest bias rows to keep.
#' @param p_cut Significance cutoff used for counting flagged rows.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method returns a compact interaction-bias summary:
#' - interaction facets/order and analyzed cell counts
#' - effect-size profile (`|bias|` mean/max, significant cell count)
#' - fixed-effect chi-square block
#' - iteration-end convergence indicators
#' - top rows ranked by absolute t
#'
#' @section Interpreting output:
#' - `overview`: interaction order, analyzed cells, and effect-size profile.
#' - `chi_sq`: fixed-effect test block.
#' - `final_iteration`: end-of-loop status from the bias routine.
#' - `top_rows`: strongest bias contrasts by `|t|`.
#'
#' @section Typical workflow:
#' 1. Estimate interactions with [estimate_bias()].
#' 2. Check `summary(bias)` for screen-positive and unstable cells.
#' 3. Use [bias_interaction_report()] or [plot_bias_interaction()] for details.
#'
#' @return An object of class `summary.mfrm_bias` with:
#' - `overview`: interaction facets/order, cell counts, and effect-size profile
#' - `chi_sq`: fixed-effect chi-square block
#' - `final_iteration`: end-of-iteration status row
#' - `top_rows`: highest-`|t|` interaction rows
#' - `notes`: short interpretation notes
#' @seealso [estimate_bias()], [bias_interaction_report()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' summary(bias)
#' @export
summary.mfrm_bias <- function(object, digits = 3, top_n = 10, p_cut = 0.05, ...) {
  if (!is.list(object) || is.null(object$table) || nrow(object$table) == 0) {
    stop("`object` must be non-empty output from estimate_bias().")
  }

  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))
  p_cut <- max(0, min(1, as.numeric(p_cut[1])))

  bias_tbl <- tibble::as_tibble(object$table)
  chi_tbl <- tibble::as_tibble(object$chi_sq %||% tibble::tibble())
  iter_tbl <- tibble::as_tibble(object$iteration %||% tibble::tibble())
  spec <- extract_bias_facet_spec(object)
  interaction_facets <- if (!is.null(spec)) spec$facets else unique(c(
    as.character(object$facet_a[1] %||% NA_character_),
    as.character(object$facet_b[1] %||% NA_character_)
  ))
  interaction_facets <- interaction_facets[!is.na(interaction_facets) & nzchar(interaction_facets)]
  if (length(interaction_facets) == 0) interaction_facets <- c("Unknown")
  interaction_label <- paste(interaction_facets, collapse = " x ")
  interaction_order <- length(interaction_facets)
  interaction_mode <- ifelse(interaction_order > 2, "higher_order", "pairwise")

  abs_bias <- abs(suppressWarnings(as.numeric(bias_tbl$`Bias Size`)))
  p_vals <- suppressWarnings(as.numeric(bias_tbl$`Prob.`))
  sig_n <- sum(is.finite(p_vals) & p_vals <= p_cut, na.rm = TRUE)

  overview <- tibble::tibble(
    FacetPair = interaction_label,
    InteractionOrder = interaction_order,
    InteractionMode = interaction_mode,
    Cells = nrow(bias_tbl),
    MeanAbsBias = mean(abs_bias, na.rm = TRUE),
    MaxAbsBias = max(abs_bias, na.rm = TRUE),
    Significant = sig_n,
    SignificantCut = p_cut,
    ScreenPositive = sig_n,
    ScreeningCut = p_cut
  )

  final_iteration <- tibble::tibble()
  if (nrow(iter_tbl) > 0) {
    final_iteration <- iter_tbl |>
      dplyr::slice_tail(n = 1)
  }

  top_rows <- tibble::tibble()
  level_cols <- if (!is.null(spec)) {
    spec$level_cols
  } else if (all(c("FacetA_Level", "FacetB_Level") %in% names(bias_tbl))) {
    c("FacetA_Level", "FacetB_Level")
  } else {
    character(0)
  }
  keep <- c(level_cols, "Bias Size", "S.E.", "t", "Prob.", "Obs-Exp Average")
  if (all(keep %in% names(bias_tbl))) {
    top_rows <- bias_tbl |>
      dplyr::mutate(AbsT = abs(.data$t)) |>
      dplyr::arrange(dplyr::desc(.data$AbsT)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select(dplyr::all_of(c(keep, "AbsT")))
    if (length(level_cols) == length(interaction_facets)) {
      names(top_rows)[seq_along(level_cols)] <- interaction_facets
      top_rows <- dplyr::mutate(
        top_rows,
        Pair = do.call(paste, c(top_rows[interaction_facets], sep = " | ")),
        .before = 1
      )
    }
  }

  notes <- character(0)
  if (nrow(iter_tbl) > 0) {
    tail_row <- iter_tbl[nrow(iter_tbl), , drop = FALSE]
    tail_cells <- suppressWarnings(as.numeric(tail_row$BiasCells[1]))
    if (is.finite(tail_cells) && tail_cells > 0) {
      notes <- c(notes, "Bias iteration may not have fully stabilized (BiasCells > 0 at final step).")
    }
  }
  if (interaction_order > 2) {
    notes <- c(notes, "Higher-order interaction mode is active; pairwise contrasts should be interpreted from dedicated 2-way runs.")
  }
  if (isTRUE(object$mixed_sign)) {
    notes <- c(notes, as.character(object$direction_note %||% "Selected interaction facets mix score orientations; use neutral higher/lower-than-expected wording."))
  }
  if (length(notes) == 0) {
    notes <- "No immediate warnings from bias summary."
  }

  out <- list(
    overview = overview,
    chi_sq = chi_tbl,
    final_iteration = final_iteration,
    top_rows = top_rows,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_bias"
  out
}

#' @export
print.summary.mfrm_bias <- function(x, ...) {
  digits <- as.integer(x$digits %||% 3L)
  if (!is.finite(digits)) digits <- 3L

  cat("Many-Facet Rasch Bias Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0) {
    ov <- round_numeric_df(as.data.frame(x$overview), digits = digits)[1, , drop = FALSE]
    cat(sprintf("  Interaction facets: %s | Cells: %s\n", ov$FacetPair, ov$Cells))
    if ("InteractionOrder" %in% names(ov) && "InteractionMode" %in% names(ov)) {
      cat(sprintf("  Order: %s | Mode: %s\n", ov$InteractionOrder, ov$InteractionMode))
    }
    cat(sprintf(
      "  Mean |Bias|: %s | Max |Bias|: %s | Screen-positive (p <= %.3f): %s\n",
      ov$MeanAbsBias, ov$MaxAbsBias, as.numeric(ov$ScreeningCut), ov$ScreenPositive
    ))
  }

  if (!is.null(x$chi_sq) && nrow(x$chi_sq) > 0) {
    cat("\nFixed-effect chi-square\n")
    print(round_numeric_df(as.data.frame(x$chi_sq), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$final_iteration) && nrow(x$final_iteration) > 0) {
    cat("\nFinal iteration status\n")
    print(round_numeric_df(as.data.frame(x$final_iteration), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$top_rows) && nrow(x$top_rows) > 0) {
    cat("\nTop |t| bias rows\n")
    print(round_numeric_df(as.data.frame(x$top_rows), digits = digits), row.names = FALSE)
  }

  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize an `mfrm_fit` object in a user-friendly format
#'
#' @param object Output from [fit_mfrm()].
#' @param digits Number of digits for printed numeric values.
#' @param top_n Number of extreme facet/person rows shown in summaries.
#' @param ... Reserved for generic compatibility.
#'
#' @details
#' This method provides a compact, human-readable summary oriented to reporting.
#' It returns a structured object and prints:
#' - model fit overview (N, LogLik, AIC/BIC, convergence)
#' - facet-level estimate distribution (mean/SD/range)
#' - person measure distribution
#' - step/threshold checks
#' - high/low person measures and extreme facet levels
#'
#' @section Interpreting output:
#' - `overview`: convergence and information criteria.
#' - `facet_overview`: per-facet spread and range of estimates.
#' - `person_overview`: distribution of person measures.
#' - `step_overview`: threshold spread and monotonicity checks.
#' - `top_person` / `top_facet`: extreme estimates for quick triage.
#'
#' @section Typical workflow:
#' 1. Fit model with [fit_mfrm()].
#' 2. Run `summary(fit)` for first-pass diagnostics.
#' 3. Continue with [diagnose_mfrm()] for element-level fit checks.
#'
#' @return An object of class `summary.mfrm_fit` with:
#' - `overview`: global model/fit indicators
#' - `facet_overview`: per-facet estimate distribution summary
#' - `person_overview`: person-measure distribution summary
#' - `step_overview`: threshold/step diagnostics
#' - `top_person`: highest/lowest person measures
#' - `top_facet`: extreme facet-level estimates
#' - `notes`: short interpretation notes
#' @seealso [fit_mfrm()], [diagnose_mfrm()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' summary(fit)
#' @export
summary.mfrm_fit <- function(object, digits = 3, top_n = 5, ...) {
  if (is.null(object$summary) || nrow(object$summary) == 0) {
    stop("`object` does not contain fit summary information.")
  }

  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  overview <- tibble::as_tibble(object$summary)
  person_raw <- object$facets$person
  if (is.null(person_raw)) person_raw <- tibble::tibble()
  other_raw <- object$facets$others
  if (is.null(other_raw)) other_raw <- tibble::tibble()
  step_raw <- object$steps
  if (is.null(step_raw)) step_raw <- tibble::tibble()

  person_tbl <- tibble::as_tibble(person_raw)
  other_tbl <- tibble::as_tibble(other_raw)
  step_tbl <- tibble::as_tibble(step_raw)

  facet_overview <- tibble::tibble()
  if (nrow(other_tbl) > 0 && all(c("Facet", "Estimate") %in% names(other_tbl))) {
    facet_overview <- other_tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(
        Levels = dplyr::n(),
        MeanEstimate = mean(.data$Estimate, na.rm = TRUE),
        SDEstimate = stats::sd(.data$Estimate, na.rm = TRUE),
        MinEstimate = min(.data$Estimate, na.rm = TRUE),
        MaxEstimate = max(.data$Estimate, na.rm = TRUE),
        Span = .data$MaxEstimate - .data$MinEstimate,
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$Facet)
  }

  person_overview <- tibble::tibble()
  if (nrow(person_tbl) > 0 && "Estimate" %in% names(person_tbl)) {
    person_overview <- tibble::tibble(
      Persons = nrow(person_tbl),
      Mean = mean(person_tbl$Estimate, na.rm = TRUE),
      SD = stats::sd(person_tbl$Estimate, na.rm = TRUE),
      Median = stats::median(person_tbl$Estimate, na.rm = TRUE),
      Min = min(person_tbl$Estimate, na.rm = TRUE),
      Max = max(person_tbl$Estimate, na.rm = TRUE),
      Span = max(person_tbl$Estimate, na.rm = TRUE) - min(person_tbl$Estimate, na.rm = TRUE)
    )

    if ("SD" %in% names(person_tbl)) {
      person_overview$MeanPosteriorSD <- mean(person_tbl$SD, na.rm = TRUE)
    }
  }

  step_overview <- tibble::tibble()
  if (nrow(step_tbl) > 0 && all(c("Step", "Estimate") %in% names(step_tbl))) {
    ord <- order(step_tbl$Step)
    step_vals <- as.numeric(step_tbl$Estimate[ord])
    monotonic <- if (length(step_vals) <= 1) TRUE else all(diff(step_vals) >= -sqrt(.Machine$double.eps))
    step_overview <- tibble::tibble(
      Steps = nrow(step_tbl),
      Min = min(step_tbl$Estimate, na.rm = TRUE),
      Max = max(step_tbl$Estimate, na.rm = TRUE),
      Span = max(step_tbl$Estimate, na.rm = TRUE) - min(step_tbl$Estimate, na.rm = TRUE),
      Monotonic = monotonic
    )
  }

  facet_extremes <- tibble::tibble()
  if (nrow(other_tbl) > 0 && all(c("Facet", "Level", "Estimate") %in% names(other_tbl))) {
    facet_extremes <- other_tbl |>
      dplyr::mutate(AbsEstimate = abs(.data$Estimate)) |>
      dplyr::arrange(dplyr::desc(.data$AbsEstimate)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select("Facet", "Level", "Estimate")
  }

  person_high <- tibble::tibble()
  person_low <- tibble::tibble()
  if (nrow(person_tbl) > 0 && all(c("Person", "Estimate") %in% names(person_tbl))) {
    person_high <- person_tbl |>
      dplyr::arrange(dplyr::desc(.data$Estimate)) |>
      dplyr::slice_head(n = top_n)
    person_low <- person_tbl |>
      dplyr::arrange(.data$Estimate) |>
      dplyr::slice_head(n = top_n)
  }

  notes <- character(0)
  if ("Converged" %in% names(overview) && !isTRUE(overview$Converged[1])) {
    notes <- c(notes, "Optimization did not converge; interpret parameter estimates cautiously.")
  }
  if (nrow(step_overview) > 0 && !isTRUE(step_overview$Monotonic[1])) {
    notes <- c(notes, "Step estimates are not monotonic; verify category functioning.")
  }
  if (length(notes) == 0) {
    notes <- "No immediate warnings from fit-level summary checks."
  }

  out <- list(
    overview = overview,
    facet_overview = facet_overview,
    person_overview = person_overview,
    step_overview = step_overview,
    facet_extremes = facet_extremes,
    person_high = person_high,
    person_low = person_low,
    notes = notes,
    digits = digits
  )
  class(out) <- "summary.mfrm_fit"
  out
}

round_numeric_df <- function(df, digits = 3L) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  out <- df
  numeric_cols <- vapply(out, is.numeric, logical(1))
  out[numeric_cols] <- lapply(out[numeric_cols], round, digits = digits)
  out
}

#' @export
print.summary.mfrm_fit <- function(x, ...) {
  digits <- x$digits
  if (is.null(digits) || !is.finite(digits)) digits <- 3L
  overview <- round_numeric_df(as.data.frame(x$overview), digits = digits)

  cat("Many-Facet Rasch Model Summary\n")
  if (nrow(overview) > 0) {
    ov <- overview[1, , drop = FALSE]
    cat(sprintf("  Model: %s | Method: %s\n", ov$Model, ov$Method))
    cat(sprintf("  N: %s | Persons: %s | Facets: %s | Categories: %s\n", ov$N, ov$Persons, ov$Facets, ov$Categories))
    cat(sprintf("  LogLik: %s | AIC: %s | BIC: %s\n", ov$LogLik, ov$AIC, ov$BIC))
    cat(sprintf("  Converged: %s | Iterations: %s\n", ifelse(isTRUE(ov$Converged), "Yes", "No"), ov$Iterations))
  }

  if (nrow(x$facet_overview) > 0) {
    cat("\nFacet overview\n")
    print(round_numeric_df(as.data.frame(x$facet_overview), digits = digits), row.names = FALSE)
  }

  if (nrow(x$person_overview) > 0) {
    cat("\nPerson measure distribution\n")
    print(round_numeric_df(as.data.frame(x$person_overview), digits = digits), row.names = FALSE)
  }

  if (nrow(x$step_overview) > 0) {
    cat("\nStep parameter summary\n")
    print(round_numeric_df(as.data.frame(x$step_overview), digits = digits), row.names = FALSE)
  }

  if (nrow(x$facet_extremes) > 0) {
    cat("\nMost extreme facet levels (|estimate|)\n")
    print(round_numeric_df(as.data.frame(x$facet_extremes), digits = digits), row.names = FALSE)
  }

  if (nrow(x$person_high) > 0) {
    cat("\nHighest person measures\n")
    print(round_numeric_df(as.data.frame(x$person_high), digits = digits), row.names = FALSE)
  }

  if (nrow(x$person_low) > 0) {
    cat("\nLowest person measures\n")
    print(round_numeric_df(as.data.frame(x$person_low), digits = digits), row.names = FALSE)
  }

  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }

  invisible(x)
}
