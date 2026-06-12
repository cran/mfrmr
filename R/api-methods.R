
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
    "missing_preview", "column_review", "metric_checks", "column_summary", "metric_summary",
    "conquest_population", "conquest_item_estimates", "conquest_case_eap"
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
      "FACETS-style branch: table columns mirror the output-contract naming."
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
  plot_routes <- as.data.frame(object$public_plot_routes %||% data.frame(), stringsAsFactors = FALSE)
  payload_names <- names(object$plot_payloads %||% list())

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
  if (length(payload_names) > 0) {
    notes <- c(
      notes,
      paste0(
        "Reusable draw-free plot data are available in `plot_payloads`: ",
        paste(payload_names, collapse = ", "),
        "."
      )
    )
  }

  out <- list(
    summary_kind = "visual_summaries",
    overview = overview,
    summary = warning_counts,
    preview_name = if (nrow(preview_tbl) > 0) "warning_counts" else "",
    preview = preview_tbl,
    settings = crosswalk,
    plot_routes = utils::head(plot_routes, n = top_n),
    notes = notes,
    digits = digits,
    summary_counts = summary_counts
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

summarize_export_bundle <- function(object, digits = 3, top_n = 10) {
  written <- bundle_component_table(object, "written_files")
  summary_tbl <- bundle_component_table(object, "summary")
  settings <- bundle_settings_table(object$settings)
  format_summary <- bundle_export_format_summary(written)
  artifact_catalog <- bundle_export_artifact_catalog(written)
  artifact_preview <- utils::head(artifact_catalog, n = top_n)
  reporting_map <- bundle_export_reporting_map("bundle")

  out <- list(
    summary_kind = "export_bundle",
    overview = bundle_known_overview(
      object,
      obj_class = "mfrm_export_bundle",
      preview_name = if (nrow(written) > 0) "written_files" else NA_character_,
      preview_rows = min(nrow(written), top_n)
    ),
    summary = summary_tbl,
    preview_name = if (nrow(written) > 0) "written_files" else "",
    preview = utils::head(written, n = top_n),
    settings = settings,
    format_summary = format_summary,
    artifact_catalog = artifact_catalog,
    artifact_preview = artifact_preview,
    reporting_map = reporting_map,
    notes = as.character(object$notes %||% "Bundle export completed successfully."),
    digits = digits
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

summarize_summary_appendix_export <- function(object, digits = 3, top_n = 10) {
  written <- bundle_component_table(object, "written_files")
  summary_tbl <- bundle_component_table(object, "summary")
  selection_summary <- bundle_component_table(object, "selection_summary")
  selection_table_summary <- bundle_component_table(object, "selection_table_summary")
  selection_handoff_table_summary <- bundle_component_table(object, "selection_handoff_table_summary")
  selection_handoff_preset_summary <- bundle_component_table(object, "selection_handoff_preset_summary")
  selection_handoff_summary <- bundle_component_table(object, "selection_handoff_summary")
  selection_handoff_bundle_summary <- bundle_component_table(object, "selection_handoff_bundle_summary")
  selection_handoff_role_summary <- bundle_component_table(object, "selection_handoff_role_summary")
  selection_handoff_role_section_summary <- bundle_component_table(object, "selection_handoff_role_section_summary")
  selection_role_summary <- bundle_component_table(object, "selection_role_summary")
  selection_section_summary <- bundle_component_table(object, "selection_section_summary")
  selection_catalog <- bundle_component_table(object, "selection_catalog")
  settings <- bundle_settings_table(object$settings)
  format_summary <- bundle_export_format_summary(written)
  artifact_catalog <- bundle_export_artifact_catalog(written)
  artifact_preview <- utils::head(artifact_catalog, n = top_n)
  reporting_map <- bundle_export_reporting_map("appendix")

  out <- list(
    summary_kind = "summary_appendix_export",
    overview = bundle_known_overview(
      object,
      obj_class = "mfrm_summary_appendix_export",
      preview_name = if (nrow(written) > 0) "written_files" else NA_character_,
      preview_rows = min(nrow(written), top_n)
    ),
    summary = summary_tbl,
    preview_name = if (nrow(written) > 0) "written_files" else "",
    preview = utils::head(written, n = top_n),
    settings = settings,
    format_summary = format_summary,
    artifact_catalog = artifact_catalog,
    artifact_preview = artifact_preview,
    selection_summary = selection_summary,
    selection_table_summary = selection_table_summary,
    selection_handoff_table_summary = selection_handoff_table_summary,
    selection_handoff_preset_summary = selection_handoff_preset_summary,
    selection_handoff_summary = selection_handoff_summary,
    selection_handoff_bundle_summary = selection_handoff_bundle_summary,
    selection_handoff_role_summary = selection_handoff_role_summary,
    selection_handoff_role_section_summary = selection_handoff_role_section_summary,
    selection_role_summary = selection_role_summary,
    selection_section_summary = selection_section_summary,
    selection_catalog = utils::head(selection_catalog, n = top_n),
    reporting_map = reporting_map,
    notes = as.character(object$notes %||% "Appendix export completed successfully."),
    digits = digits
  )
  class(out) <- "summary.mfrm_bundle"
  out
}

bundle_export_format_summary <- function(written_files) {
  written_files <- as.data.frame(written_files %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(written_files) == 0L || !"Format" %in% names(written_files)) {
    return(data.frame())
  }

  formats <- sort(table(as.character(written_files$Format)), decreasing = TRUE)
  data.frame(
    Format = names(formats),
    Files = as.integer(formats),
    stringsAsFactors = FALSE
  )
}

bundle_export_artifact_catalog <- function(written_files) {
  written_files <- as.data.frame(written_files %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(written_files) == 0L) {
    return(data.frame())
  }

  component <- as.character(written_files$Component %||% "")
  format <- as.character(written_files$Format %||% "")
  path_base <- basename(as.character(written_files$Path %||% ""))
  artifact_group <- ifelse(
    startsWith(component, "summary_"),
    "summary_surface",
    ifelse(
      grepl("_html$", component),
      "html_review",
      ifelse(
        grepl("_zip$", component),
        "archive",
        "analysis_bundle"
      )
    )
  )
  recommended_use <- ifelse(
    artifact_group == "summary_surface",
    "manuscript / appendix handoff",
    ifelse(
      artifact_group == "html_review",
      "human-readable review",
      ifelse(
        artifact_group == "archive",
        "single-file transfer",
        "analysis export"
      )
    )
  )

  data.frame(
    Component = component,
    Format = format,
    ArtifactGroup = artifact_group,
    RecommendedUse = recommended_use,
    File = path_base,
    stringsAsFactors = FALSE
  )
}

bundle_export_reporting_map <- function(kind = c("bundle", "appendix")) {
  kind <- match.arg(kind)
  if (identical(kind, "appendix")) {
    return(data.frame(
      Area = c(
        "Export counts",
        "Artifact catalog / handoff",
        "Human-readable review",
        "Downstream manuscript bridge"
      ),
      CoveredHere = c("yes", "yes", "yes", "yes"),
      CompanionOutput = c(
        "summary(appendix)$summary / format_summary",
        "summary(appendix)$artifact_catalog",
        "appendix HTML artifact when include_html = TRUE",
        "CSV appendix tables + apa_table() / manuscript QA"
      ),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    Area = c(
      "Export counts",
      "Artifact catalog / file inventory",
      "Human-readable review",
      "Replay / archival bridge"
    ),
    CoveredHere = c("yes", "yes", "yes", "yes"),
    CompanionOutput = c(
      "summary(bundle)$summary / format_summary",
      "summary(bundle)$artifact_catalog",
      "bundle HTML artifact when requested",
      "manifest / replay script / zip bundle"
    ),
    stringsAsFactors = FALSE
  )
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
  caveats_tbl <- as.data.frame(object$caveats %||% data.frame(), stringsAsFactors = FALSE)

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
    caveats = caveats_tbl,
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
    ScoreSEMethod = as.character(settings$score_se_method %||% NA_character_),
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
  quality_overview <- as.data.frame(object$quality_overview %||% data.frame(), stringsAsFactors = FALSE)
  quality_flags <- as.data.frame(object$quality_flags %||% data.frame(), stringsAsFactors = FALSE)
  caveats <- as.data.frame(object$caveats %||% data.frame(), stringsAsFactors = FALSE)
  usage_summary <- as.data.frame(object$category_usage_summary %||% data.frame(), stringsAsFactors = FALSE)
  notes <- "Data quality summary for missingness, row status, score support, and category usage."
  if (nrow(quality_overview) > 0 && "Status" %in% names(quality_overview)) {
    high_areas <- sum(tolower(as.character(quality_overview$Status)) %in% "high", na.rm = TRUE)
    review_areas <- sum(tolower(as.character(quality_overview$Status)) %in% "review", na.rm = TRUE)
    notes <- c(
      notes,
      paste0(
        "QC overview: ", high_areas, " high-priority area(s), ",
        review_areas, " review area(s)."
      )
    )
  }
  if (nrow(quality_flags) > 0) {
    high_flags <- if ("Severity" %in% names(quality_flags)) {
      sum(tolower(as.character(quality_flags$Severity)) %in% "high", na.rm = TRUE)
    } else {
      0L
    }
    notes <- c(
      notes,
      paste0(
        "Priority QC flags: ", nrow(quality_flags),
        " flag(s), including ", high_flags, " high-severity flag(s)."
      )
    )
  } else {
    notes <- c(notes, "No priority QC flags were found in the supplied data-quality checks.")
  }
  if (nrow(usage_summary) > 0) {
    zero_levels <- sum(usage_summary$ZeroCategories > 0, na.rm = TRUE)
    internal_zero_levels <- sum(usage_summary$IntermediateZeroCategories > 0, na.rm = TRUE)
    sparse_levels <- sum(usage_summary$SparseCategories > 0, na.rm = TRUE)
    if (zero_levels > 0L || sparse_levels > 0L) {
      notes <- c(
        notes,
        paste0(
          "Facet-level category use: ", zero_levels,
          " level(s) have zero-count categories; ", internal_zero_levels,
          " have zero-count intermediate categories; ", sparse_levels,
          " have sparse non-zero categories."
        )
      )
    }
  }
  if (nrow(caveats) > 0 && "Message" %in% names(caveats)) {
    notes <- c(
      notes,
      paste(utils::head(as.character(caveats$Message), 3L), collapse = " ")
    )
  }
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_data_quality",
    summary_candidates = "summary",
    preview_candidates = c("quality_flags", "quality_overview", "facet_response_patterns", "caveats", "category_usage_summary", "score_support_review", "row_review", "category_counts", "model_match", "unknown_elements"),
    settings_candidates = "settings",
    notes = notes,
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

summarize_network_analysis_bundle <- function(object, digits = 3, top_n = 10) {
  summary_tbl <- as.data.frame(object$summary %||% data.frame(), stringsAsFactors = FALSE)
  notes <- "Network metrics summarize design connectedness and linking vulnerability, not person ability or rater quality."
  if (nrow(summary_tbl) > 0L) {
    components <- suppressWarnings(as.integer(summary_tbl$Components[1]))
    cut_n <- suppressWarnings(as.integer(summary_tbl$ArticulationPoints[1]))
    bridge_n <- suppressWarnings(as.integer(summary_tbl$Bridges[1]))
    notes <- c(
      notes,
      paste0(
        "Graph: ", summary_tbl$Nodes[1], " node(s), ",
        summary_tbl$Edges[1], " edge(s), ", components, " component(s), ",
        cut_n, " articulation point(s), ", bridge_n, " bridge edge(s)."
      )
    )
  }
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_network_analysis",
    summary_candidates = "summary",
    preview_candidates = c("cut_nodes", "bridge_edges", "facet_summary", "node_metrics"),
    settings_candidates = "settings",
    notes = notes,
    digits = digits,
    top_n = top_n
  )
}

summarize_rater_network_bundle <- function(object, digits = 3, top_n = 10) {
  summary_tbl <- as.data.frame(object$summary %||% data.frame(), stringsAsFactors = FALSE)
  notes <- "Rater-network metrics summarize observed pairwise rater relationships; they are not Rasch logit estimates or formal fit statistics."
  if (nrow(summary_tbl) > 0L) {
    notes <- c(
      notes,
      paste0(
        "Graph: ", summary_tbl$Raters[1], " rater node(s), ",
        summary_tbl$Edges[1], " edge(s), mode = ",
        summary_tbl$Mode[1], "."
      )
    )
  }
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_rater_network",
    summary_candidates = "summary",
    preview_candidates = c("node_metrics", "edge_metrics", "pair_metrics", "caveats"),
    settings_candidates = "settings",
    notes = notes,
    digits = digits,
    top_n = top_n
  )
}

summarize_halo_network_bundle <- function(object, digits = 3, top_n = 10) {
  summary_tbl <- as.data.frame(object$summary %||% data.frame(), stringsAsFactors = FALSE)
  notes <- "Halo-network metrics summarize rater-by-criterion score-profile similarity; they are screening diagnostics, not causal halo evidence by themselves."
  if (nrow(summary_tbl) > 0L) {
    notes <- c(
      notes,
      paste0(
        "Graph: ", summary_tbl$Nodes[1], " rater-by-criterion node(s), ",
        summary_tbl$Edges[1], " retained edge(s), halo minus non-halo mean = ",
        format(round(suppressWarnings(as.numeric(summary_tbl$HaloMinusNonHalo[1])), digits), nsmall = digits),
        "."
      )
    )
  }
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_halo_network",
    summary_candidates = "summary",
    preview_candidates = c("halo_summary_by_rater", "edge_metrics", "node_metrics", "caveats"),
    settings_candidates = "settings",
    notes = notes,
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

summarize_precision_review_bundle <- function(object, digits = 3, top_n = 10) {
  profile_tbl <- bundle_component_table(object, "profile")
  checks_tbl <- bundle_component_table(object, "checks")
  fit_sep_tbl <- bundle_component_table(object, "fit_separation_basis")
  notes_tbl <- bundle_component_table(object, "approximation_notes")

  flagged_n <- if ("Status" %in% names(checks_tbl)) {
    sum(as.character(checks_tbl$Status) %in% c("review", "warn"), na.rm = TRUE)
  } else {
    NA_integer_
  }

  summary_tbl <- data.frame(
    Method = resolve_public_mfrm_method(
      summary_method = profile_tbl$Method[1] %||% NA_character_
    ),
    PrecisionTier = as.character(profile_tbl$PrecisionTier[1] %||% NA_character_),
    SupportsFormalInference = isTRUE(profile_tbl$SupportsFormalInference[1] %||% FALSE),
    Checks = nrow(checks_tbl),
    ReviewOrWarn = flagged_n,
    FitSeparationRows = nrow(fit_sep_tbl),
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

  out <- summarize_known_bundle(
    object = object,
    obj_class = "mfrm_precision_review",
    summary_candidates = character(0),
    preview_candidates = c("checks", "fit_separation_basis", "approximation_notes", "profile"),
    settings_candidates = "settings",
    notes = c(notes, "Fit/separation basis rows state source grounding and validation-use boundaries."),
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
  out$profile <- profile_tbl
  out$checks <- checks_tbl
  out$fit_separation_basis <- fit_sep_tbl
  out$approximation_notes <- notes_tbl
  class(out) <- unique(c("summary.mfrm_precision_review", class(out)))
  out
}

summarize_facets_contract_bundle <- function(object, digits = 3, top_n = 10) {
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
    notes <- "Compatibility checks completed."
  }

  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_facets_contract_review",
    summary_candidates = character(0),
    preview_candidates = c("missing_preview", "column_review", "metric_checks"),
    settings_candidates = "settings",
    notes = notes,
    digits = digits,
    top_n = top_n,
    summary_override = overall_tbl
  )
}

summarize_facets_fit_review_bundle <- function(object, digits = 3, top_n = 10) {
  summary_tbl <- as.data.frame(object$summary %||% data.frame(), stringsAsFactors = FALSE)
  external_tbl <- as.data.frame(object$external_comparison %||% data.frame(), stringsAsFactors = FALSE)

  notes <- character(0)
  if (nrow(summary_tbl) > 0) {
    changed <- suppressWarnings(as.integer(summary_tbl$FlagChangedByDf[1]))
    df_sensitive <- suppressWarnings(as.integer(summary_tbl$DfSensitiveRows[1]))
    duplicate_external <- suppressWarnings(as.integer(summary_tbl$ExternalDuplicateKeyRows[1] %||% 0L))
    external_review <- suppressWarnings(as.integer(summary_tbl$ExternalNeedsReview[1]))
    if (is.finite(changed) && changed > 0) {
      notes <- c(notes, paste0("FACETS-style df changed |ZSTD|>2 flags for ", changed, " element(s)."))
    }
    if (is.finite(df_sensitive) && df_sensitive > 0) {
      notes <- c(notes, paste0("Engine-vs-FACETS-style df/ZSTD differences need review for ", df_sensitive, " element(s)."))
    }
    if (is.finite(external_review) && external_review > 0) {
      notes <- c(notes, paste0("External FACETS rows need review for ", external_review, " matched element(s)."))
    }
    if (is.finite(duplicate_external) && duplicate_external > 0) {
      notes <- c(notes, paste0("External FACETS table has duplicate Facet x Level rows: ", duplicate_external, " row(s)."))
    }
  }
  if (length(notes) == 0) {
    notes <- "Fit-standardization review completed; inspect df_sensitivity before interpreting ZSTD flags."
  }

  preview_candidates <- if (nrow(external_tbl) > 0) {
    c("external_table_quality", "external_comparison", "df_sensitive", "df_sensitivity", "guidance")
  } else {
    c("df_sensitive", "df_sensitivity", "guidance", "standardization")
  }

  out <- summarize_known_bundle(
    object = object,
    obj_class = "mfrm_facets_fit_review",
    summary_candidates = character(0),
    preview_candidates = preview_candidates,
    settings_candidates = "settings",
    notes = notes,
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
  out$standardization <- bundle_component_table(object, "standardization")
  out$df_sensitivity <- bundle_component_table(object, "df_sensitivity")
  out$df_sensitive <- bundle_component_table(object, "df_sensitive")
  out$df_sensitivity_summary <- bundle_component_table(object, "df_sensitivity_summary")
  out$external_table_quality <- bundle_component_table(object, "external_table_quality")
  out$external_comparison <- external_tbl
  out$guidance <- bundle_component_table(object, "guidance")
  class(out) <- unique(c("summary.mfrm_facets_fit_review", class(out)))
  out
}

summarize_facets_fit_df_guide_bundle <- function(object, digits = 3, top_n = 10) {
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_facets_fit_df_guide",
    summary_candidates = "summary",
    preview_candidates = c("decision_guide", "column_guide", "formula_guide"),
    settings_candidates = "settings",
    notes = "FACETS-style fit df and ZSTD comparison guide. Compare MnSq first, then df and ZSTD.",
    digits = digits,
    top_n = top_n
  )
}

summarize_fit_measures_bundle <- function(object, digits = 3, top_n = 10) {
  summary_tbl <- bundle_component_table(object, "summary")
  df_summary <- bundle_component_table(object, "df_sensitivity_summary")
  df_sensitive <- bundle_component_table(object, "df_sensitive")

  notes <- "Fit-measures table with directional underfit/overfit screening."
  if (nrow(df_summary) > 0L) {
    changed <- suppressWarnings(as.integer(df_summary$FlagChangedByDfRows[1] %||% 0L))
    sensitive <- nrow(df_sensitive)
    if (is.finite(changed) && changed > 0L) {
      notes <- c(notes, paste0("FACETS-style df changed |ZSTD| flag status for ", changed, " row(s)."))
    } else if (sensitive > 0L) {
      notes <- c(notes, paste0("FACETS-style df changed ZSTD interpretation for ", sensitive, " row(s); inspect df_sensitive."))
    }
  }

  out <- summarize_known_bundle(
    object = object,
    obj_class = "mfrm_fit_measures",
    summary_candidates = character(0),
    preview_candidates = c("df_sensitive", "table", "facets_table"),
    settings_candidates = "settings",
    notes = notes,
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
  out$table <- bundle_component_table(object, "table")
  out$facets_table <- bundle_component_table(object, "facets_table")
  out$status_summary <- bundle_component_table(object, "status_summary")
  out$threshold_profiles <- bundle_component_table(object, "threshold_profiles")
  out$profile_summary <- bundle_component_table(object, "profile_summary")
  out$profile_summary_by_facet <- bundle_component_table(object, "profile_summary_by_facet")
  out$profile_summary_overall <- bundle_component_table(object, "profile_summary_overall")
  out$df_sensitivity <- bundle_component_table(object, "df_sensitivity")
  out$df_sensitive <- bundle_component_table(object, "df_sensitive")
  out$df_sensitivity_summary <- bundle_component_table(object, "df_sensitivity_summary")
  out$underfit <- bundle_component_table(object, "underfit")
  out$overfit <- bundle_component_table(object, "overfit")
  out$mixed <- bundle_component_table(object, "mixed")
  class(out) <- unique(c("summary.mfrm_fit_measures", class(out)))
  out
}

reference_benchmark_validation_scope <- function(object) {
  settings <- object$settings %||% list()
  cases <- as.character(settings$cases %||% character(0))
  has_conquest_dry_run <- "synthetic_conquest_overlap_dry_run" %in% cases ||
    nrow(bundle_component_table(object, "conquest_overlap_checks")) > 0L
  has_population_policy <- "synthetic_latent_regression_omit" %in% cases ||
    nrow(bundle_component_table(object, "population_policy_checks")) > 0L
  external_validation <- isTRUE(settings$external_validation %||% FALSE)

  data.frame(
    Area = c(
      "Package reference check",
      "Latent-regression omission policy",
      "ConQuest-overlap package-side check",
      "External ConQuest validation"
    ),
    Status = c(
      "active",
      if (has_population_policy) "checked" else "not requested",
      if (has_conquest_dry_run) "package-side check only" else "not requested",
      if (external_validation) "claimed" else "not performed"
    ),
    Evidence = c(
      "case_summary and component check tables",
      if (has_population_policy) "`population_policy_checks`" else "request `synthetic_latent_regression_omit`",
      if (has_conquest_dry_run) "`conquest_overlap_checks`" else "request `synthetic_conquest_overlap_dry_run`",
      if (external_validation) "settings$external_validation is TRUE" else "settings$external_validation is FALSE"
    ),
    Interpretation = c(
      "Use as a package reference check.",
      "Complete-case omission behavior is reviewed only when the omission case is requested.",
      "The check covers package-side export, normalization, and review plumbing only.",
      "Actual external ConQuest output tables are required before making an external validation claim."
    ),
    stringsAsFactors = FALSE
  )
}

summarize_reference_benchmark_bundle <- function(object, digits = 3, top_n = 10) {
  case_summary <- bundle_component_table(object, "case_summary")
  if (nrow(case_summary) == 0L && ncol(case_summary) == 0L) {
    case_summary <- bundle_component_table(object, "summary")
  }
  out <- summarize_known_bundle(
    object = object,
    obj_class = "mfrm_reference_benchmark",
    summary_candidates = character(0),
    preview_candidates = c("fit_runs", "table"),
    settings_candidates = "settings",
    notes = object$notes %||% "Reference-case check; not an external validation study.",
    digits = digits,
    top_n = top_n,
    summary_override = utils::head(case_summary, n = top_n)
  )
  out$validation_scope <- reference_benchmark_validation_scope(object)
  out$conquest_overlap_checks <- utils::head(
    bundle_component_table(object, "conquest_overlap_checks"),
    n = top_n
  )
  out$population_policy_checks <- utils::head(
    bundle_component_table(object, "population_policy_checks"),
    n = top_n
  )
  out
}

conquest_overlap_command_scope <- function(object) {
  command <- paste(as.character(object$conquest_command %||% ""), collapse = "\n")
  summary_tbl <- bundle_component_table(object, "summary")
  facet <- as.character(summary_tbl$Facet[1] %||% "item")
  covariate <- as.character(summary_tbl$Covariate[1] %||% "")
  has_widths <- grepl("pidwidth=", command, fixed = TRUE) &&
    grepl("keepswidth=", command, fixed = TRUE)
  has_block_comment <- grepl("^\\s*/\\*", command) &&
    grepl("\\*/", command)
  has_exports <- all(vapply(
    c(
      "export parameters ! filetype=csv",
      "export reg_coefficients ! filetype=csv",
      "export covariance ! filetype=csv",
      "show cases ! estimates=eap, filetype=csv"
    ),
    function(pattern) grepl(pattern, command, fixed = TRUE),
    logical(1)
  ))

  data.frame(
    Area = c(
      "ConQuest command template",
      "Command-comment syntax",
      "Official command-reference alignment",
      "Overlap model scope",
      "External output requirements",
      "External comparison scope"
    ),
    Status = c(
      "template only",
      if (has_block_comment) "block comments" else "review required",
      if (has_widths) "explicit CSV widths" else "review required",
      "narrow overlap only",
      if (has_exports) "requested" else "review required",
      "not claimed"
    ),
    Evidence = c(
      "bundle$conquest_command",
      if (has_block_comment) "command text starts with /* and closes with */" else "ConQuest block comment markers not detected",
      if (has_widths) "pidwidth and keepswidth are present" else "pidwidth or keepswidth is missing",
      paste0("binary ", facet, " facet with numeric covariate `", covariate, "`"),
      "parameters, reg_coefficients, covariance, and cases EAP CSV outputs",
      "requires external ConQuest execution and extracted output-table review"
    ),
    Interpretation = c(
      "Use the command text as a starting point for a local ConQuest run, not as an executed benchmark.",
      "Generated comments follow the documented ConQuest block-comment style rather than FACETS-style leading asterisks.",
      "CSV input with PID/keeps variables needs explicit widths in the command template.",
      "The bundle does not generalize to full many-facet or polytomous ConQuest workflows.",
      "Review and combine external parameter, beta, sigma, and case outputs before review normalization.",
      "External comparison remains scoped until external outputs are reviewed and tolerances are justified."
    ),
    stringsAsFactors = FALSE
  )
}

summarize_conquest_overlap_bundle <- function(object, digits = 3, top_n = 10) {
  out <- summarize_known_bundle(
    object = object,
    obj_class = "mfrm_conquest_overlap_bundle",
    summary_candidates = "summary",
    preview_candidates = c("comparison_targets", "item_map", "mfrmr_population", "mfrmr_case_eap"),
    settings_candidates = "settings",
    notes = object$notes %||% "ConQuest-overlap bundle prepared.",
    digits = digits,
    top_n = top_n
  )
  out$conquest_command_scope <- conquest_overlap_command_scope(object)
  out$conquest_output_contract <- utils::head(
    bundle_component_table(object, "conquest_output_contract"),
    n = top_n
  )
  out
}

conquest_overlap_normalization_scope <- function(object) {
  summary_tbl <- bundle_component_table(object, "summary")
  first_numeric_sum <- function(cols) {
    vals <- numeric(0)
    for (nm in cols) {
      if (!nm %in% names(summary_tbl) || nrow(summary_tbl) == 0L) next
      vals <- c(vals, suppressWarnings(as.numeric(summary_tbl[[nm]][1])))
    }
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) 0L else as.integer(sum(vals))
  }

  row_n <- first_numeric_sum(c("PopulationRows", "ItemRows", "CaseRows"))
  duplicate_n <- first_numeric_sum(c("PopulationDuplicateIDs", "ItemDuplicateIDs", "CaseDuplicateIDs"))
  non_numeric_n <- first_numeric_sum(c("PopulationNonNumeric", "ItemNonNumeric", "CaseNonNumeric"))
  review_n <- duplicate_n + non_numeric_n

  data.frame(
    Area = c(
      "Extracted table normalization",
      "Raw ConQuest text parsing",
      "Bundle matching",
      "Pre-review table check"
    ),
    Status = c(
      "active",
      "not performed",
      "deferred to review",
      if (review_n > 0L) "review required" else "none detected"
    ),
    Evidence = c(
      paste0(row_n, " standardized row(s)"),
      "already extracted CSV/TSV/TXT or data.frame inputs only",
      "review_conquest_overlap() matches rows against the exported bundle",
      paste0(duplicate_n, " duplicate ID(s); ", non_numeric_n, " non-numeric estimate cell(s)")
    ),
    Interpretation = c(
      "Population, item, and case tables have been converted to the mfrmr review contract.",
      "This object does not prove that raw ConQuest report text was parsed correctly.",
      "Identifier matching and numerical comparison are intentionally handled by the review step.",
      "Resolve duplicate IDs or non-numeric estimates before treating the review as clean."
    ),
    stringsAsFactors = FALSE
  )
}

summarize_conquest_overlap_tables_bundle <- function(object, digits = 3, top_n = 10) {
  out <- summarize_known_bundle(
    object = object,
    obj_class = "mfrm_conquest_overlap_tables",
    summary_candidates = "summary",
    preview_candidates = c("conquest_population", "conquest_item_estimates", "conquest_case_eap"),
    settings_candidates = "settings",
    notes = object$notes %||% "ConQuest-overlap extracted tables normalized.",
    digits = digits,
    top_n = top_n
  )
  out$normalization_scope <- conquest_overlap_normalization_scope(object)
  out
}

conquest_overlap_review_scope <- function(object) {
  overall <- bundle_component_table(object, "overall")
  attention <- bundle_component_table(object, "attention_items")
  attention_n <- if ("AttentionItems" %in% names(overall)) {
    suppressWarnings(as.integer(overall$AttentionItems[1]))
  } else {
    nrow(attention)
  }
  if (!is.finite(attention_n)) attention_n <- nrow(attention)

  data.frame(
    Area = c(
      "User-supplied table review",
      "Raw ConQuest text parsing",
      "External comparison scope",
      "Attention items"
    ),
    Status = c(
      "active",
      "not performed",
      "not claimed",
      if (attention_n > 0L) "review required" else "none detected"
    ),
    Evidence = c(
      "population, item, and case comparison tables",
      "use normalize_conquest_overlap_files() / normalize_conquest_overlap_tables() first",
      "inspect differences, attention items, extraction steps, and tolerances",
      paste0(attention_n, " attention item(s)")
    ),
    Interpretation = c(
      "The review compares supplied normalized tables against the mfrmr overlap bundle.",
      "This helper does not parse raw ConQuest report text.",
      "Numerical agreement is limited to the documented overlap and supplied tables.",
      "Nonzero attention items indicate missing, duplicate, non-numeric, or unmatched rows to resolve."
    ),
    stringsAsFactors = FALSE
  )
}

summarize_conquest_overlap_review_bundle <- function(object, digits = 3, top_n = 10) {
  out <- summarize_known_bundle(
    object = object,
    obj_class = "mfrm_conquest_overlap_review",
    summary_candidates = "overall",
    preview_candidates = "attention_items",
    settings_candidates = "settings",
    notes = object$notes %||% "ConQuest-overlap review completed.",
    digits = digits,
    top_n = top_n
  )
  out$review_scope <- conquest_overlap_review_scope(object)
  out
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
  top_n <- max(1L, as.integer(top_n))
  stacked <- bundle_component_table(object, "stacked")

  first_present <- function(candidates) {
    candidates <- as.character(candidates %||% character(0))
    hit <- candidates[candidates %in% names(stacked)]
    if (length(hit) == 0L) NA_character_ else hit[1]
  }
  numeric_column <- function(candidates) {
    nm <- first_present(candidates)
    if (is.na(nm)) return(rep(NA_real_, nrow(stacked)))
    suppressWarnings(as.numeric(stacked[[nm]]))
  }
  character_column <- function(candidates) {
    nm <- first_present(candidates)
    if (is.na(nm)) return(rep(NA_character_, nrow(stacked)))
    as.character(stacked[[nm]])
  }
  first_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else x[1]
  }
  mean_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
  }
  max_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else max(x)
  }
  compact_counts <- function(x, preferred = character(0)) {
    x <- as.character(x %||% character(0))
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) return(character(0))
    counts <- table(x)
    keys <- names(counts)
    key_order <- c(intersect(preferred, keys), setdiff(keys, preferred))
    paste(paste0(key_order, "=", as.integer(counts[key_order])), collapse = ", ")
  }

  obs_avg <- numeric_column(c("ObservedAverage", "Obsvd Average"))
  fair_m <- numeric_column(c("AdjustedAverage", "Fair(M) Average"))
  fair_m_se <- numeric_column(c("AdjustedAverageSE", "Fair(M) S.E."))
  fair_m_ci_lower <- numeric_column(c("AdjustedAverageCI_Lower", "Fair(M) CI Lower"))
  fair_m_ci_upper <- numeric_column(c("AdjustedAverageCI_Upper", "Fair(M) CI Upper"))
  fair_m_ci_level <- numeric_column(c("AdjustedAverageCI_Level", "Fair(M) CI Level"))
  fair_m_method <- character_column(c("AdjustedAverageSEMethod", "Fair(M) S.E. Method"))
  fair_m_status <- character_column(c("AdjustedAverageSEStatus", "Fair(M) S.E. Status"))

  mean_abs_gap <- NA_real_
  if (length(obs_avg) == length(fair_m) && length(obs_avg) > 0) {
    dif <- abs(obs_avg - fair_m)
    dif <- dif[is.finite(dif)]
    if (length(dif) > 0) mean_abs_gap <- mean(dif)
  }

  has_fair_se_columns <- any(c(
    "AdjustedAverageSE", "Fair(M) S.E.",
    "AdjustedAverageCI_Lower", "Fair(M) CI Lower",
    "AdjustedAverageCI_Upper", "Fair(M) CI Upper",
    "AdjustedAverageSEStatus", "Fair(M) S.E. Status"
  ) %in% names(stacked))
  fair_se_requested <- isTRUE(object$settings$fair_se %||% FALSE) || has_fair_se_columns
  status_values <- fair_m_status[!is.na(fair_m_status) & nzchar(fair_m_status)]
  method_values <- fair_m_method[!is.na(fair_m_method) & nzchar(fair_m_method)]
  fair_se_status <- if (!fair_se_requested) {
    "not_requested"
  } else if (length(status_values) > 0L) {
    compact_counts(status_values, preferred = c("ok", "regularized", "not available", "not_available"))
  } else if (any(is.finite(fair_m_se))) {
    "available"
  } else {
    "not_available"
  }
  fair_se_method <- if (!fair_se_requested) {
    "not_requested"
  } else if (length(method_values) > 0L) {
    paste(unique(method_values), collapse = "; ")
  } else {
    "not_available"
  }

  summary_tbl <- data.frame(
    Facets = if ("Facet" %in% names(stacked)) length(unique(as.character(stacked$Facet))) else length(object$by_facet %||% list()),
    Levels = nrow(stacked),
    MeanAbsObservedFairM = mean_abs_gap,
    FairSERequested = fair_se_requested,
    FairSEAvailableRows = sum(is.finite(fair_m_se)),
    FairSEUnavailableRows = if (fair_se_requested) sum(!is.finite(fair_m_se)) else 0L,
    FairSEMethod = fair_se_method,
    FairSEStatus = fair_se_status,
    MeanAdjustedAverageSE = mean_finite(fair_m_se),
    MaxAdjustedAverageSE = max_finite(fair_m_se),
    AdjustedAverageCILevel = first_finite(fair_m_ci_level),
    stringsAsFactors = FALSE
  )

  preview_tbl <- data.frame()
  if (nrow(stacked) > 0L) {
    preview_tbl <- data.frame(
      Facet = character_column("Facet"),
      Level = character_column(c("Level", "Element")),
      ObservedAverage = obs_avg,
      AdjustedAverage = fair_m,
      stringsAsFactors = FALSE
    )
    if (fair_se_requested) {
      preview_tbl$AdjustedAverageSE <- fair_m_se
      preview_tbl$AdjustedAverageCI_Lower <- fair_m_ci_lower
      preview_tbl$AdjustedAverageCI_Upper <- fair_m_ci_upper
      preview_tbl$AdjustedAverageSEStatus <- fair_m_status
    }
    abs_gap <- abs(obs_avg - fair_m)
    if (fair_se_requested && any(is.finite(fair_m_se))) {
      ord <- order(!is.finite(fair_m_se), -abs_gap, na.last = TRUE)
    } else {
      ord <- order(-abs_gap, na.last = TRUE)
    }
    preview_tbl <- preview_tbl[ord, , drop = FALSE]
    preview_tbl <- utils::head(preview_tbl, n = top_n)
  }

  notes <- "Adjusted-score reference summary by facet level."
  if (fair_se_requested) {
    notes <- c(
      notes,
      "Fair-average SE columns summarize structural delta-method uncertainty when available; unavailable rows are reported explicitly."
    )
  } else {
    notes <- c(
      notes,
      "Fair-average structural SE columns are omitted unless requested by `fair_se = TRUE`."
    )
  }

  out <- summarize_known_bundle(
    object = object,
    obj_class = "mfrm_fair_average",
    summary_candidates = character(0),
    preview_candidates = c("stacked", "raw_by_facet"),
    settings_candidates = "settings",
    notes = notes,
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
  if (nrow(preview_tbl) > 0L) {
    out$preview_name <- "stacked"
    out$preview <- preview_tbl
    if (!is.null(out$overview) && nrow(out$overview) > 0L) {
      out$overview$PreviewComponent[1] <- "stacked"
      out$overview$PreviewRows[1] <- nrow(preview_tbl)
    }
  }
  out
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
  summary_tbl <- bundle_component_table(object, "summary")
  base_note <- if ("MarginalFitAvailable" %in% names(summary_tbl) &&
    isTRUE(summary_tbl$MarginalFitAvailable[1] %||% FALSE)) {
    "Rating-scale diagnostics with category usage, fit, threshold ordering, and strict marginal-fit companions."
  } else {
    "Rating-scale diagnostics with category usage, fit, and threshold ordering."
  }
  caveats <- as.data.frame(object$caveats %||% data.frame(), stringsAsFactors = FALSE)
  notes <- c(
    base_note,
    if (nrow(caveats) > 0 && "Message" %in% names(caveats)) as.character(caveats$Message) else character(0)
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_rating_scale",
    summary_candidates = "summary",
    preview_candidates = c("category_table", "threshold_table"),
    settings_candidates = character(0),
    notes = notes,
    digits = digits,
    top_n = top_n
  )
}

summarize_category_structure_bundle <- function(object, digits = 3, top_n = 10) {
  cat_tbl <- bundle_component_table(object, "category_table")
  marginal_available <- is.list(object$marginal_fit) && isTRUE(object$marginal_fit$available)
  marginal_summary <- as.data.frame(object$marginal_fit$summary %||% data.frame(), stringsAsFactors = FALSE)
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
    DiagnosticMode = as.character(object$diagnostic_mode %||% "legacy"),
    MarginalFitAvailable = marginal_available,
    MarginalFlaggedCategories = if ("MarginalFitFlag" %in% names(cat_tbl)) {
      sum(as.logical(cat_tbl$MarginalFitFlag), na.rm = TRUE)
    } else {
      NA_integer_
    },
    MarginalOverallRMSD = if (marginal_available) marginal_summary$OverallRMSD[1] %||% NA_real_ else NA_real_,
    MarginalMaxAbsStdResidual = if (marginal_available) marginal_summary$OverallMaxAbsStdResidual[1] %||% NA_real_ else NA_real_,
    stringsAsFactors = FALSE
  )
  base_note <- if (marginal_available) {
    "Category-structure diagnostics with mode boundaries, half-score reference points, and strict marginal-fit companions."
  } else {
    "Category-structure diagnostics with mode boundaries and half-score reference points."
  }
  caveats <- as.data.frame(object$caveats %||% data.frame(), stringsAsFactors = FALSE)
  notes <- c(
    base_note,
    if (nrow(caveats) > 0 && "Message" %in% names(caveats)) as.character(caveats$Message) else character(0)
  )
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_category_structure",
    summary_candidates = character(0),
    preview_candidates = c("category_table", "mode_boundaries", "mean_halfscore_points"),
    settings_candidates = "settings",
    notes = notes,
    digits = digits,
    top_n = top_n,
    summary_override = summary_tbl
  )
}

summarize_category_curves_bundle <- function(object, digits = 3, top_n = 10) {
  graph_tbl <- bundle_component_table(object, "graphfile")
  ogive <- bundle_component_table(object, "expected_ogive")
  probs <- bundle_component_table(object, "probabilities")
  cumulative <- bundle_component_table(object, "cumulative_probabilities")
  cumulative_boundaries <- bundle_component_table(object, "cumulative_boundaries")
  cat_info <- bundle_component_table(object, "category_information")
  prob_cols <- grep("^Prob:", names(graph_tbl), value = TRUE)
  count_unique <- function(tbl, col) {
    if (!is.data.frame(tbl) || nrow(tbl) == 0L || !col %in% names(tbl)) return(NA_integer_)
    length(unique(as.character(tbl[[col]][!is.na(tbl[[col]])])))
  }
  numeric_max <- function(tbl, col) {
    if (!is.data.frame(tbl) || nrow(tbl) == 0L || !col %in% names(tbl)) return(NA_real_)
    x <- suppressWarnings(as.numeric(tbl[[col]]))
    if (!any(is.finite(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  }
  boundary_status <- if (nrow(cumulative_boundaries) > 0L && "BoundaryStatus" %in% names(cumulative_boundaries)) {
    as.character(cumulative_boundaries$BoundaryStatus)
  } else {
    character(0)
  }
  boundary_review <- sum(!is.na(boundary_status) & nzchar(boundary_status) & boundary_status != "in_range")
  boundary_outside <- sum(boundary_status == "outside_theta_range", na.rm = TRUE)
  boundary_multiple <- sum(boundary_status == "multiple_crossings", na.rm = TRUE)
  summary_tbl <- data.frame(
    Metric = c(
      "curve_groups",
      "theta_points",
      "categories",
      "legacy_graph_rows",
      "expected_ogive_rows",
      "probability_rows",
      "probability_columns",
      "cumulative_probability_rows",
      "cumulative_boundary_rows",
      "category_information_rows",
      "boundary_rows_needing_review",
      "boundary_rows_outside_theta_range",
      "boundary_rows_with_multiple_crossings",
      "max_total_information",
      "max_category_information"
    ),
    Value = c(
      count_unique(probs, "CurveGroup"),
      count_unique(probs, "Theta"),
      count_unique(probs, "Category"),
      nrow(graph_tbl),
      nrow(ogive),
      nrow(probs),
      length(prob_cols),
      nrow(cumulative),
      nrow(cumulative_boundaries),
      nrow(cat_info),
      boundary_review,
      boundary_outside,
      boundary_multiple,
      numeric_max(ogive, "Information"),
      numeric_max(cat_info, "CategoryInformation")
    ),
    stringsAsFactors = FALSE
  )
  boundary_note <- if (nrow(cumulative_boundaries) == 0L) {
    "No cumulative .5 boundary rows were returned."
  } else if (boundary_review > 0L) {
    paste0(
      "Review ", boundary_review,
      " cumulative .5 boundary row(s): boundaries outside the theta grid or with multiple crossings should not be read as single stable thresholds."
    )
  } else {
    "Cumulative .5 boundary rows are in range with a single crossing where reported."
  }
  summarize_known_bundle(
    object = object,
    obj_class = "mfrm_category_curves",
    summary_candidates = character(0),
    preview_candidates = c("cumulative_boundaries", "category_information", "expected_ogive", "probabilities", "graphfile"),
    settings_candidates = "settings",
    notes = c(
      "Category-curve bundle with probabilities, cumulative probabilities, total information, and category-specific information.",
      "Category-specific information contributions sum to total information at the same curve and theta point.",
      boundary_note
    ),
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
#' - `mfrm_residual_pca`, `mfrm_specifications`, `mfrm_data_quality`,
#'   `mfrm_fit_measures`
#' - `mfrm_iteration_report`, `mfrm_subset_connectivity`, `mfrm_facet_statistics`
#' - `mfrm_facets_contract_review`, `mfrm_facets_fit_review`,
#'   `mfrm_facets_fit_df_guide`, `mfrm_reference_benchmark`
#'
#' @section Interpreting output:
#' - `overview`: class, component count, and selected preview component.
#' - `summary`: one-row aggregate block when supplied by the bundle.
#' - `preview`: first `top_n` rows from the main table-like component.
#' - `settings`: resolved option values if available.
#' - `validation_scope`: internal-versus-external validation scope when
#'   summarizing `mfrm_reference_benchmark`.
#' - `conquest_command_scope`: ConQuest command-template scope when summarizing
#'   `mfrm_conquest_overlap_bundle`.
#' - `conquest_output_contract`: requested ConQuest outputs and review handoff
#'   when summarizing `mfrm_conquest_overlap_bundle`.
#' - `normalization_scope`: extracted-table normalization scope when summarizing
#'   `mfrm_conquest_overlap_tables`.
#' - `review_scope`: supplied-table review scope when summarizing
#'   `mfrm_conquest_overlap_review`.
#' - `conquest_overlap_checks` / `population_policy_checks`: specialized
#'   benchmark check previews when summarizing `mfrm_reference_benchmark`.
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
#'   fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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
  if (inherits(object, "mfrm_fit_measures")) {
    return(summarize_fit_measures_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_iteration_report")) {
    return(summarize_iteration_report_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_subset_connectivity")) {
    return(summarize_subset_connectivity_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_network_analysis")) {
    return(summarize_network_analysis_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_rater_network")) {
    return(summarize_rater_network_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_halo_network")) {
    return(summarize_halo_network_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_facet_statistics")) {
    return(summarize_facet_statistics_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_precision_review")) {
    return(summarize_precision_review_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_facets_contract_review")) {
    return(summarize_facets_contract_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_facets_fit_review")) {
    return(summarize_facets_fit_review_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_facets_fit_df_guide")) {
    return(summarize_facets_fit_df_guide_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_reference_benchmark")) {
    return(summarize_reference_benchmark_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_conquest_overlap_bundle")) {
    return(summarize_conquest_overlap_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_conquest_overlap_tables")) {
    return(summarize_conquest_overlap_tables_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_conquest_overlap_review")) {
    return(summarize_conquest_overlap_review_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_export_bundle")) {
    return(summarize_export_bundle(object, digits = digits, top_n = top_n))
  }
  if (inherits(object, "mfrm_summary_appendix_export")) {
    return(summarize_summary_appendix_export(object, digits = digits, top_n = top_n))
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
    mfrm_category_curves = list(title = "mfrmr Category Curves Summary", summary = "Curve grid summary", preview = "Boundary / curve rows"),
    mfrm_measurable = list(title = "mfrmr Measurable Summary", summary = "Run overview", preview = "Facet/category rows"),
    mfrm_unexpected_after_bias = list(title = "mfrmr Unexpected-after-Bias Summary", summary = "After-bias threshold summary", preview = "After-bias flagged rows"),
    mfrm_output_bundle = list(title = "mfrmr Output File Bundle Summary", summary = "Output overview", preview = "Output preview rows"),
    mfrm_residual_pca = list(title = "mfrmr Residual PCA Summary", summary = "PCA overview", preview = "Eigenvalue / loading rows"),
    mfrm_specifications = list(title = "mfrmr Specifications Summary", summary = "Specification header", preview = "Specification rows"),
    mfrm_data_quality = list(title = "mfrmr Data Quality Summary", summary = "Data quality overview", preview = "Review rows"),
    mfrm_fit_measures = list(title = "mfrmr Fit Measures Summary", summary = "Fit-status overview", preview = "Fit-measure rows"),
    mfrm_facets_fit_df_guide = list(title = "mfrmr FACETS Fit df Guide", summary = "Guide overview", preview = "Comparison steps"),
    mfrm_iteration_report = list(title = "mfrmr Iteration Report Summary", summary = "Iteration overview", preview = "Iteration rows"),
    mfrm_subset_connectivity = list(title = "mfrmr Subset Connectivity Summary", summary = "Subset overview", preview = "Subset/node rows"),
    mfrm_facet_statistics = list(title = "mfrmr Facet Profile Summary", summary = "Facet-profile overview", preview = "Facet-profile rows"),
    mfrm_precision_review = list(title = "mfrmr Precision Review Summary", summary = "Precision overview", preview = "Review checks"),
    mfrm_facets_contract_review = list(title = "mfrmr FACETS Output Contract Review Summary", summary = "Contract review overview", preview = "Lowest-coverage items"),
    mfrm_reference_review = list(title = "mfrmr Reference Review Summary", summary = "Reference review overview", preview = "Attention items"),
    mfrm_reference_benchmark = list(title = "mfrmr Reference Case Check Summary", summary = "Case check summary", preview = "Reference-case fit runs"),
    mfrm_reporting_checklist = list(title = "mfrmr Reporting Checklist Summary", summary = "Checklist coverage", preview = "Checklist items"),
    mfrm_bias_collection = list(title = "mfrmr Bias Collection Summary", summary = "Interaction summary", preview = "Per-pair results"),
    mfrm_manifest = list(title = "mfrmr Manifest Summary", summary = "Analysis overview", preview = "Manifest tables"),
    mfrm_replay_script = list(title = "mfrmr Replay Script Summary", summary = "Replay settings", preview = "Script text"),
    mfrm_conquest_overlap_bundle = list(title = "mfrmr ConQuest Overlap Bundle Summary", summary = "Overlap scope summary", preview = "Comparison targets"),
    mfrm_conquest_overlap_tables = list(title = "mfrmr ConQuest Overlap Table Normalization Summary", summary = "Normalization overview", preview = "Standardized tables"),
    mfrm_conquest_overlap_review = list(title = "mfrmr ConQuest Overlap Review Summary", summary = "Comparison overview", preview = "Attention items"),
    export_bundle = list(title = "mfrmr Export Bundle Summary", summary = "Export overview", preview = "Written files"),
    mfrm_export_bundle = list(title = "mfrmr Export Bundle Summary", summary = "Export overview", preview = "Written files"),
    summary_appendix_export = list(title = "mfrmr Summary Appendix Export Summary", summary = "Appendix export overview", preview = "Written appendix files"),
    mfrm_summary_appendix_export = list(title = "mfrmr Summary Appendix Export Summary", summary = "Appendix export overview", preview = "Written appendix files"),
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

collect_mfrm_population_caveats <- function(population_overview = NULL,
                                            population_design = NULL,
                                            population_coefficients = NULL) {
  out <- empty_mfrm_caveats()
  population_overview <- as.data.frame(population_overview %||% data.frame(), stringsAsFactors = FALSE)
  population_design <- as.data.frame(population_design %||% data.frame(), stringsAsFactors = FALSE)
  population_coefficients <- as.data.frame(population_coefficients %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(population_overview) == 0L || !"PopulationModel" %in% names(population_overview) ||
      !isTRUE(population_overview$PopulationModel[1])) {
    return(out)
  }

  add_caveat <- function(severity, condition, details, message, action) {
    out <<- rbind(
      out,
      data.frame(
        Area = "population_model",
        Severity = severity,
        Condition = condition,
        Categories = "",
        CategoryType = "",
        Message = message,
        RecommendedAction = action,
        Details = details,
        stringsAsFactors = FALSE
      )
    )
    invisible(NULL)
  }

  omitted_persons <- suppressWarnings(as.integer(population_overview$OmittedPersons[1] %||% 0L))
  omitted_rows <- suppressWarnings(as.integer(population_overview$OmittedRows[1] %||% NA_integer_))
  if (is.finite(omitted_persons) && omitted_persons > 0L) {
    add_caveat(
      severity = "warning",
      condition = "population_complete_case_omission",
      details = paste0("OmittedPersons=", omitted_persons, "; OmittedRows=", omitted_rows),
      message = paste0(
        "Latent-regression fit omitted ",
        omitted_persons,
        " person(s) and ",
        if (is.finite(omitted_rows)) omitted_rows else NA_integer_,
        " response row(s) under the population-data policy."
      ),
      action = "Document the complete-case policy and verify that omitted persons do not change the intended estimation population."
    )
  }

  if (nrow(population_design) > 0L && all(c("Column", "IsIntercept", "ZeroVariance") %in% names(population_design))) {
    is_intercept <- as.logical(population_design$IsIntercept %||% FALSE)
    zero_variance <- as.logical(population_design$ZeroVariance %||% FALSE)
    flagged <- population_design$Column[zero_variance & !is_intercept]
    flagged <- flagged[!is.na(flagged) & nzchar(as.character(flagged))]
    if (length(flagged) > 0L) {
      add_caveat(
        severity = "warning",
        condition = "population_design_zero_variance",
        details = paste(as.character(flagged), collapse = ", "),
        message = paste0(
          "Latent-regression design matrix contains non-intercept zero-variance column(s): ",
          paste(as.character(flagged), collapse = ", "),
          "."
        ),
        action = "Remove or recode zero-variance background variables before interpreting population-model coefficients."
      )
    }
  }

  if (nrow(population_design) > 0L && all(c("Column", "Complete") %in% names(population_design))) {
    incomplete <- population_design$Column[!as.logical(population_design$Complete %||% TRUE)]
    incomplete <- incomplete[!is.na(incomplete) & nzchar(as.character(incomplete))]
    if (length(incomplete) > 0L) {
      add_caveat(
        severity = "warning",
        condition = "population_design_incomplete",
        details = paste(as.character(incomplete), collapse = ", "),
        message = paste0(
          "Latent-regression design matrix still has incomplete column(s): ",
          paste(as.character(incomplete), collapse = ", "),
          "."
        ),
        action = "Review `population_policy`, person-level background data, and complete-case filtering before reporting the population model."
      )
    }
  }

  residual_variance <- suppressWarnings(as.numeric(population_overview$ResidualVariance[1] %||% NA_real_))
  if (!is.finite(residual_variance) || residual_variance <= 0) {
    add_caveat(
      severity = "warning",
      condition = "population_residual_variance_unstable",
      details = paste0("ResidualVariance=", residual_variance),
      message = "Latent-regression residual variance is non-finite or non-positive.",
      action = "Review population-model convergence and covariate design before using population-model posterior scoring."
    )
  }

  if (nrow(population_coefficients) == 0L) {
    add_caveat(
      severity = "warning",
      condition = "population_coefficients_missing",
      details = "",
      message = "Population model is marked active but no population coefficients are available in the summary.",
      action = "Inspect `fit$population$coefficients` before reporting latent-regression effects."
    )
  }

  out
}

format_population_coding_value <- function(x) {
  if (is.null(x) || length(x) == 0L) return("")
  if (is.function(x)) return("<function>")
  if (is.matrix(x) || is.data.frame(x)) {
    return(paste0(class(x)[1], "[", nrow(x), "x", ncol(x), "]"))
  }
  vals <- as.character(x)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0L) return("")
  paste(vals, collapse = ", ")
}

population_coding_summary_table <- function(population) {
  population <- population %||% list()
  xlevels <- population$xlevels %||% list()
  contrasts <- population$contrasts %||% list()
  if (is.null(xlevels)) xlevels <- list()
  if (!is.list(xlevels)) {
    xlevel_names <- names(xlevels)
    xlevels <- as.list(xlevels)
    names(xlevels) <- xlevel_names
  }
  if (is.null(contrasts)) contrasts <- list()
  if (!is.list(contrasts)) {
    contrast_names <- names(contrasts)
    contrasts <- as.list(contrasts)
    names(contrasts) <- contrast_names
  }
  xlevel_names <- names(xlevels) %||% character(0)
  contrast_names <- names(contrasts) %||% character(0)
  vars <- unique(c(xlevel_names[nzchar(xlevel_names)], contrast_names[nzchar(contrast_names)]))
  if (length(vars) == 0L) {
    return(tibble::tibble(
      Variable = character(0),
      LevelCount = integer(0),
      Levels = character(0),
      Contrast = character(0),
      EncodedColumns = character(0),
      CodingNote = character(0)
    ))
  }

  design_columns <- as.character(population$design_columns %||% character(0))
  if (length(design_columns) == 0L && !is.null(population$design_matrix)) {
    design_columns <- as.character(colnames(population$design_matrix) %||% character(0))
  }
  design_columns <- design_columns[!is.na(design_columns) & nzchar(design_columns)]

  encoded_columns <- vapply(vars, function(var) {
    cols <- design_columns[design_columns == var | startsWith(design_columns, var)]
    paste(cols, collapse = ", ")
  }, character(1))
  level_text <- vapply(vars, function(var) format_population_coding_value(xlevels[[var]]), character(1))
  level_count <- vapply(vars, function(var) length(xlevels[[var]] %||% character(0)), integer(1))
  contrast_text <- vapply(vars, function(var) format_population_coding_value(contrasts[[var]]), character(1))
  has_xlevels <- nzchar(level_text)
  has_contrasts <- nzchar(contrast_text)
  tibble::tibble(
    Variable = as.character(vars),
    LevelCount = as.integer(level_count),
    Levels = level_text,
    Contrast = contrast_text,
    EncodedColumns = encoded_columns,
    CodingNote = dplyr::case_when(
      has_xlevels & has_contrasts ~ "stored levels and contrast",
      has_xlevels ~ "stored levels",
      has_contrasts ~ "stored contrast",
      TRUE ~ "not recorded"
    )
  )
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
      for (line in x$notes) cat(" - ", line, "\n", sep = "")
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
    if (!is.null(x$plot_routes) && nrow(x$plot_routes) > 0) {
      cat("\nPublic plot routes\n")
      print(as.data.frame(x$plot_routes), row.names = FALSE)
    }
    if (length(x$notes) > 0) {
      cat("\nNotes\n")
      for (line in x$notes) cat(" - ", line, "\n", sep = "")
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
      for (line in x$notes) cat(" - ", line, "\n", sep = "")
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
  if (!is.null(x$format_summary) && nrow(x$format_summary) > 0) {
    print_bundle_section("Format summary", x$format_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$artifact_preview) && nrow(x$artifact_preview) > 0) {
    print_bundle_section("Artifact catalog", x$artifact_preview, digits = digits, round_numeric = FALSE)
  } else if (!is.null(x$artifact_catalog) && nrow(x$artifact_catalog) > 0) {
    print_bundle_section("Artifact catalog", x$artifact_catalog, digits = digits, round_numeric = FALSE)
  }
  if (!is.null(x$preview) && nrow(x$preview) > 0) {
    preview_title <- labels$preview
    if (!is.null(x$preview_name) && !is.na(x$preview_name) && nzchar(x$preview_name)) {
      preview_title <- paste0(preview_title, ": ", x$preview_name)
    }
    print_bundle_section(preview_title, x$preview, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$validation_scope) && nrow(x$validation_scope) > 0) {
    print_bundle_section("Validation scope", x$validation_scope, digits = digits, round_numeric = FALSE)
  }
  if (!is.null(x$conquest_command_scope) && nrow(x$conquest_command_scope) > 0) {
    print_bundle_section("ConQuest command scope", x$conquest_command_scope, digits = digits, round_numeric = FALSE)
  }
  if (!is.null(x$conquest_output_contract) && nrow(x$conquest_output_contract) > 0) {
    print_bundle_section("ConQuest output contract", x$conquest_output_contract, digits = digits, round_numeric = FALSE)
  }
  if (!is.null(x$normalization_scope) && nrow(x$normalization_scope) > 0) {
    print_bundle_section("Normalization scope", x$normalization_scope, digits = digits, round_numeric = FALSE)
  }
  if (!is.null(x$review_scope) && nrow(x$review_scope) > 0) {
    print_bundle_section("Review scope", x$review_scope, digits = digits, round_numeric = FALSE)
  }
  if (!is.null(x$conquest_overlap_checks) && nrow(x$conquest_overlap_checks) > 0) {
    print_bundle_section("ConQuest-overlap checks", x$conquest_overlap_checks, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$population_policy_checks) && nrow(x$population_policy_checks) > 0) {
    print_bundle_section("Population-policy checks", x$population_policy_checks, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$settings) && nrow(x$settings) > 0) {
    print_bundle_section(labels$settings, x$settings, digits = digits, round_numeric = FALSE)
  }
  if (!is.null(x$selection_summary) && nrow(x$selection_summary) > 0) {
    print_bundle_section("Selection summary", x$selection_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$selection_table_summary) && nrow(x$selection_table_summary) > 0) {
    print_bundle_section("Selection tables", x$selection_table_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$selection_handoff_table_summary) && nrow(x$selection_handoff_table_summary) > 0) {
    print_bundle_section("Selection handoff tables", x$selection_handoff_table_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$selection_handoff_preset_summary) && nrow(x$selection_handoff_preset_summary) > 0) {
    print_bundle_section("Selection handoff presets", x$selection_handoff_preset_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$selection_handoff_summary) && nrow(x$selection_handoff_summary) > 0) {
    print_bundle_section("Selection handoff", x$selection_handoff_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$selection_handoff_bundle_summary) && nrow(x$selection_handoff_bundle_summary) > 0) {
    print_bundle_section("Selection handoff bundles", x$selection_handoff_bundle_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$selection_handoff_role_summary) && nrow(x$selection_handoff_role_summary) > 0) {
    print_bundle_section("Selection handoff roles", x$selection_handoff_role_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$selection_handoff_role_section_summary) && nrow(x$selection_handoff_role_section_summary) > 0) {
    print_bundle_section("Selection handoff role-sections", x$selection_handoff_role_section_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$selection_role_summary) && nrow(x$selection_role_summary) > 0) {
    print_bundle_section("Selection roles", x$selection_role_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$selection_section_summary) && nrow(x$selection_section_summary) > 0) {
    print_bundle_section("Selection sections", x$selection_section_summary, digits = digits, round_numeric = TRUE)
  }
  if (!is.null(x$selection_catalog) && nrow(x$selection_catalog) > 0) {
    print_bundle_section("Selection catalog", x$selection_catalog, digits = digits, round_numeric = FALSE)
  }
  if (!is.null(x$reporting_map) && nrow(x$reporting_map) > 0) {
    print_bundle_section("Reporting map", x$reporting_map, digits = digits, round_numeric = FALSE)
  }
  if (!is.null(x$caveats) && nrow(x$caveats) > 0) {
    print_bundle_section("Caveats", x$caveats, digits = digits, round_numeric = FALSE)
  }
  if (length(x$notes) > 0) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

draw_export_handoff_bundle <- function(x,
                                       type = c("formats", "artifact_groups", "selection_handoff_presets", "selection_tables", "selection_handoff", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections", "selection_bundles", "selection_roles", "selection_sections"),
                                       selection_value = c("count", "fraction"),
                                       draw = TRUE,
                                       main = NULL,
                                       palette = NULL,
                                       label_angle = 45) {
  type <- match.arg(tolower(type), c("formats", "artifact_groups", "selection_handoff_presets", "selection_tables", "selection_handoff", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections", "selection_bundles", "selection_roles", "selection_sections"))
  selection_value <- match.arg(selection_value)
  measure <- NULL
  sx <- summary(x, digits = 6)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      formats = "#4c78a8",
      artifact_groups = "#f58518",
      selection_handoff_presets = "#f4a261",
      selection_tables = "#e76f51",
      selection_handoff = "#ff9f1c",
      selection_handoff_bundles = "#5c677d",
      selection_handoff_roles = "#9c6644",
      selection_handoff_role_sections = "#7f5539",
      selection_bundles = "#54a24b",
      selection_roles = "#b279a2",
      selection_sections = "#2a9d8f"
    )
  )

  if (type == "formats") {
    tbl <- as.data.frame(sx$format_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L || !all(c("Format", "Files") %in% names(tbl))) {
      stop("No export format summary is available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$Format)
    values <- suppressWarnings(as.numeric(tbl$Files))
    cols <- pal["formats"]
    plot_title <- if (is.null(main)) "Export formats" else as.character(main[1])
  } else if (type == "artifact_groups") {
    catalog <- as.data.frame(sx$artifact_catalog %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(catalog) == 0L || !"ArtifactGroup" %in% names(catalog)) {
      stop("No export artifact catalog is available for plotting.", call. = FALSE)
    }
    counts <- sort(table(as.character(catalog$ArtifactGroup)), decreasing = TRUE)
    labels <- names(counts)
    values <- as.numeric(counts)
    cols <- pal["artifact_groups"]
    tbl <- data.frame(ArtifactGroup = labels, Files = values, stringsAsFactors = FALSE)
    plot_title <- if (is.null(main)) "Export artifact groups" else as.character(main[1])
  } else if (type == "selection_handoff_presets") {
    tbl <- as.data.frame(sx$selection_handoff_preset_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L || !all(c("Preset", "PlotReadyTables") %in% names(tbl))) {
      stop("No appendix handoff-preset summary is available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$Preset)
    measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
    values <- measure$values
    cols <- pal["selection_handoff_presets"]
    plot_title <- if (is.null(main)) "Appendix handoff by preset" else as.character(main[1])
  } else if (type == "selection_tables") {
    tbl <- as.data.frame(sx$selection_table_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L || !all(c("Table", "Rows") %in% names(tbl))) {
      stop("No appendix table selection summary is available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$Table)
    measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
    values <- measure$values
    cols <- pal["selection_tables"]
    plot_title <- if (is.null(main)) "Selected appendix tables" else as.character(main[1])
  } else if (type == "selection_bundles") {
    tbl <- as.data.frame(sx$selection_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L || !all(c("Bundle", "TablesSelected") %in% names(tbl))) {
      stop("No appendix selection summary is available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$Bundle)
    measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
    values <- measure$values
    cols <- pal["selection_bundles"]
    plot_title <- if (is.null(main)) "Appendix tables by bundle" else as.character(main[1])
  } else if (type == "selection_roles") {
    tbl <- as.data.frame(sx$selection_role_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L || !all(c("Role", "Tables") %in% names(tbl))) {
      stop("No appendix role-selection summary is available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$Role)
    measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
    values <- measure$values
    cols <- pal["selection_roles"]
    plot_title <- if (is.null(main)) "Selected appendix roles" else as.character(main[1])
  } else if (type == "selection_handoff") {
    tbl <- as.data.frame(sx$selection_handoff_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L || !all(c("AppendixSection", "PlotReadyTables") %in% names(tbl))) {
      stop("No appendix handoff summary is available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$AppendixSection)
    measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
    values <- measure$values
    cols <- pal["selection_handoff"]
    plot_title <- if (is.null(main)) "Appendix handoff by section" else as.character(main[1])
  } else if (type == "selection_handoff_bundles") {
    tbl <- as.data.frame(sx$selection_handoff_bundle_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L || !all(c("AppendixSection", "Bundle", "PlotReadyTables") %in% names(tbl))) {
      stop("No appendix handoff-bundle summary is available for plotting.", call. = FALSE)
    }
    labels <- paste0(as.character(tbl$AppendixSection), " :: ", as.character(tbl$Bundle))
    measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
    values <- measure$values
    cols <- pal["selection_handoff_bundles"]
    plot_title <- if (is.null(main)) "Appendix handoff by section and bundle" else as.character(main[1])
  } else if (type == "selection_handoff_roles") {
    tbl <- as.data.frame(sx$selection_handoff_role_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L || !all(c("Role", "PlotReadyTables") %in% names(tbl))) {
      stop("No appendix handoff-role summary is available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$Role)
    measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
    values <- measure$values
    cols <- pal["selection_handoff_roles"]
    plot_title <- if (is.null(main)) "Appendix handoff by role" else as.character(main[1])
  } else if (type == "selection_handoff_role_sections") {
    tbl <- as.data.frame(sx$selection_handoff_role_section_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L || !all(c("AppendixSection", "Role", "PlotReadyTables") %in% names(tbl))) {
      stop("No appendix handoff role-section summary is available for plotting.", call. = FALSE)
    }
    labels <- paste0(as.character(tbl$AppendixSection), " :: ", as.character(tbl$Role))
    measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
    values <- measure$values
    cols <- pal["selection_handoff_role_sections"]
    plot_title <- if (is.null(main)) "Appendix handoff by role and section" else as.character(main[1])
  } else {
    tbl <- as.data.frame(sx$selection_section_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L || !all(c("AppendixSection", "Tables") %in% names(tbl))) {
      stop("No appendix section summary is available for plotting.", call. = FALSE)
    }
    labels <- as.character(tbl$AppendixSection)
    measure <- resolve_selection_plot_measure(tbl, type, selection_value = selection_value)
    values <- measure$values
    cols <- pal["selection_sections"]
    plot_title <- if (is.null(main)) "Selected appendix sections" else as.character(main[1])
  }

  if (isTRUE(draw)) {
      barplot_rot45(
        height = values,
        labels = labels,
        col = cols,
        main = plot_title,
        ylab = if (type == "selection_tables") {
          measure$ylab
        } else if (type %in% c("selection_handoff", "selection_handoff_presets", "selection_handoff_bundles", "selection_handoff_roles", "selection_handoff_role_sections")) {
          measure$ylab
        } else if (type %in% c("selection_bundles", "selection_roles", "selection_sections")) {
          measure$ylab
        } else {
          "Files"
        },
        label_angle = label_angle,
        mar_bottom = 8.2
      )
  }

  new_mfrm_plot_data(
    "export_bundle",
    list(
      plot = type,
      selection_value = if (exists("measure")) measure$selection_value else NULL,
      table = tbl,
      title = plot_title,
      legend = new_plot_legend(
        if (type == "formats") {
          "Export format counts"
        } else if (type == "artifact_groups") {
          "Export artifact-group counts"
        } else if (!is.null(measure)) {
          measure$legend_label
        } else {
          "Appendix selection surface"
        },
        if (type == "formats") {
          "format"
        } else if (type == "artifact_groups") {
          "artifact_group"
        } else if (type == "selection_tables") {
          "table"
        } else if (type == "selection_handoff") {
          "appendix_section"
        } else if (type == "selection_handoff_bundles") {
          "appendix_bundle_section"
        } else if (type == "selection_handoff_roles") {
          "role"
        } else if (type == "selection_handoff_role_sections") {
          "appendix_role_section"
        } else if (type == "selection_bundles") {
          "bundle"
        } else if (type == "selection_roles") {
          "role"
        } else {
          "section"
        },
        "bar",
        cols
      ),
      reference_lines = new_reference_lines()
    )
  )
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
                                        type = c(
                                          "overview", "ogive", "ccc",
                                          "category_probability",
                                          "conditional_probability",
                                          "cumulative", "information",
                                          "category_information"
                                        ),
                                        draw = TRUE,
                                        main = NULL,
                                        palette = NULL,
                                        cumulative_direction = c("at_or_below", "at_or_above"),
                                        preset = c("standard", "publication", "compact", "monochrome"),
                                        show_cumulative_boundaries = TRUE,
                                        boundary_status = c("in_range", "all", "none")) {
  requested_type <- tolower(as.character(type[1] %||% "overview"))
  type_aliases <- c(
    category_probability = "ccc",
    category_probabilities = "ccc",
    conditional_probability = "ccc",
    conditional_probabilities = "ccc",
    probability = "ccc",
    probabilities = "ccc"
  )
  type <- if (requested_type %in% names(type_aliases)) type_aliases[[requested_type]] else requested_type
  type <- match.arg(type, c("overview", "ogive", "ccc", "cumulative", "information", "category_information"))
  cumulative_direction <- match.arg(tolower(as.character(cumulative_direction[1] %||% "at_or_below")),
                                    c("at_or_below", "at_or_above"))
  style <- resolve_plot_preset(preset)
  boundary_status <- match.arg(tolower(as.character(boundary_status[1] %||% "in_range")),
                               c("in_range", "all", "none"))
  show_cumulative_boundaries <- isTRUE(show_cumulative_boundaries) && !identical(boundary_status, "none")
  ogive <- as.data.frame(x$expected_ogive %||% data.frame(), stringsAsFactors = FALSE)
  probs <- as.data.frame(x$probabilities %||% data.frame(), stringsAsFactors = FALSE)
  cumulative <- as.data.frame(x$cumulative_probabilities %||% data.frame(), stringsAsFactors = FALSE)
  cumulative_boundaries <- as.data.frame(x$cumulative_boundaries %||% data.frame(), stringsAsFactors = FALSE)
  cat_info <- as.data.frame(x$category_information %||% probs, stringsAsFactors = FALSE)
  overview_panels <- data.frame(
    Panel = c(
      "Category probability",
      "Cumulative probability",
      "Total information",
      "Category-specific information"
    ),
    PlotType = c("ccc", "cumulative", "information", "category_information"),
    DataComponent = c(
      "probabilities",
      "cumulative_probabilities",
      "expected_ogive",
      "category_information"
    ),
    stringsAsFactors = FALSE
  )

  finite_range <- function(x, fallback = c(0, 1), pad = 0.05) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(fallback)
    out <- range(x)
    if (diff(out) <= sqrt(.Machine$double.eps)) {
      center <- mean(out)
      spread <- max(0.5, abs(center) * pad)
      return(center + c(-spread, spread))
    }
    out
  }
  line_palette <- function(keys) {
    keys <- unique(as.character(keys))
    defaults <- if (identical(style$name, "monochrome")) {
      stats::setNames(grDevices::gray.colors(max(3L, length(keys)), start = 0.15, end = 0.62)[seq_along(keys)], keys)
    } else {
      stats::setNames(grDevices::hcl.colors(max(3L, length(keys)), "Dark 3")[seq_along(keys)], keys)
    }
    resolve_palette(palette = palette, defaults = defaults)
  }
  line_types <- function(keys) {
    keys <- unique(as.character(keys))
    lty <- if (identical(style$name, "monochrome") && is.null(palette)) {
      rep(c(1, 2, 3, 4, 5, 6), length.out = length(keys))
    } else {
      rep(1L, length(keys))
    }
    stats::setNames(lty, keys)
  }
  boundary_lines <- function() {
    empty <- if (nrow(cumulative_boundaries) > 0L) cumulative_boundaries[0, , drop = FALSE] else data.frame()
    if (!isTRUE(show_cumulative_boundaries) || nrow(cumulative_boundaries) == 0L ||
        !all(c("CumulativeDirection", "ThurstonianThreshold") %in% names(cumulative_boundaries))) {
      return(empty)
    }
    out <- cumulative_boundaries[
      cumulative_boundaries$CumulativeDirection == cumulative_direction,
      ,
      drop = FALSE
    ]
    if (identical(boundary_status, "in_range")) {
      if ("InThetaRange" %in% names(out)) {
        out <- out[out$InThetaRange %in% TRUE, , drop = FALSE]
      }
      if ("BoundaryStatus" %in% names(out)) {
        out <- out[out$BoundaryStatus %in% "in_range", , drop = FALSE]
      }
    }
    threshold <- suppressWarnings(as.numeric(out$ThurstonianThreshold))
    out <- out[is.finite(threshold), , drop = FALSE]
    if (nrow(out) == 0L) return(empty)
    out
  }
  boundary_line_tbl <- boundary_lines()
  plot_settings <- data.frame(
    RequestedType = requested_type,
    PlotType = type,
    Preset = style$name,
    CumulativeDirection = cumulative_direction,
    ShowCumulativeBoundaries = isTRUE(show_cumulative_boundaries),
    BoundaryStatus = boundary_status,
    stringsAsFactors = FALSE
  )
  reference_line_tbl <- new_reference_lines("v", 0, "Centered theta reference", "dashed", "reference")
  reference_line_tbl <- rbind(
    reference_line_tbl,
    new_reference_lines("h", 0.5, "Cumulative probability .5 target", "dotted", "cumulative_target")
  )
  if (nrow(boundary_line_tbl) > 0L) {
    boundary_values <- suppressWarnings(as.numeric(boundary_line_tbl$ThurstonianThreshold))
    reference_line_tbl <- rbind(
      reference_line_tbl,
      new_reference_lines(
        axis = rep("v", length(boundary_values)),
        value = boundary_values,
        label = rep("Cumulative .5 boundary", length(boundary_values)),
        linetype = rep("dotted", length(boundary_values)),
        role = rep("cumulative_boundary", length(boundary_values))
      )
    )
  }
  plot_annotations <- data.frame(
    AnnotationType = as.character(reference_line_tbl$role %||% character(0)),
    Axis = as.character(reference_line_tbl$axis %||% character(0)),
    Value = suppressWarnings(as.numeric(reference_line_tbl$value %||% numeric(0))),
    Label = as.character(reference_line_tbl$label %||% character(0)),
    LineType = as.character(reference_line_tbl$linetype %||% character(0)),
    stringsAsFactors = FALSE
  )
  character_col <- function(tbl, col, default = NA_character_) {
    n <- nrow(tbl)
    if (n == 0L) return(character(0))
    if (col %in% names(tbl)) return(as.character(tbl[[col]]))
    rep(default, n)
  }
  numeric_col <- function(tbl, col, default = NA_real_) {
    n <- nrow(tbl)
    if (n == 0L) return(numeric(0))
    if (col %in% names(tbl)) return(suppressWarnings(as.numeric(tbl[[col]])))
    rep(default, n)
  }
  make_curve_long <- function() {
    rows <- list()
    if (nrow(ogive) > 0L && all(c("Theta", "ExpectedScore", "CurveGroup") %in% names(ogive))) {
      rows[[length(rows) + 1L]] <- data.frame(
        PlotType = "ogive",
        Panel = "Expected score",
        CurveGroup = character_col(ogive, "CurveGroup"),
        Theta = numeric_col(ogive, "Theta"),
        Series = character_col(ogive, "CurveGroup"),
        Category = NA_character_,
        BoundaryCategory = NA_character_,
        BoundaryOrder = NA_real_,
        CategorySet = NA_character_,
        Direction = NA_character_,
        ValueName = "ExpectedScore",
        Value = numeric_col(ogive, "ExpectedScore"),
        DisplayedByDefault = FALSE,
        Model = character_col(ogive, "Model"),
        Slope = numeric_col(ogive, "Slope"),
        stringsAsFactors = FALSE
      )
    }
    if (nrow(probs) > 0L && all(c("Theta", "Probability", "Category", "CurveGroup") %in% names(probs))) {
      rows[[length(rows) + 1L]] <- data.frame(
        PlotType = "ccc",
        Panel = "Category probability",
        CurveGroup = character_col(probs, "CurveGroup"),
        Theta = numeric_col(probs, "Theta"),
        Series = paste(character_col(probs, "CurveGroup"), character_col(probs, "Category"), sep = " | Cat "),
        Category = character_col(probs, "Category"),
        BoundaryCategory = NA_character_,
        BoundaryOrder = NA_real_,
        CategorySet = character_col(probs, "Category"),
        Direction = NA_character_,
        ValueName = "Probability",
        Value = numeric_col(probs, "Probability"),
        DisplayedByDefault = TRUE,
        Model = character_col(probs, "Model"),
        Slope = numeric_col(probs, "Slope"),
        stringsAsFactors = FALSE
      )
    }
    if (nrow(cumulative) > 0L &&
        all(c("Theta", "CumulativeProbability", "BoundaryCategory", "CurveGroup", "Direction") %in% names(cumulative))) {
      category_set <- character_col(cumulative, "CategorySet")
      missing_set <- is.na(category_set) | !nzchar(category_set)
      if (any(missing_set)) category_set[missing_set] <- character_col(cumulative, "BoundaryCategory")[missing_set]
      direction <- character_col(cumulative, "Direction")
      rows[[length(rows) + 1L]] <- data.frame(
        PlotType = "cumulative",
        Panel = "Cumulative probability",
        CurveGroup = character_col(cumulative, "CurveGroup"),
        Theta = numeric_col(cumulative, "Theta"),
        Series = paste(character_col(cumulative, "CurveGroup"), category_set, direction, sep = " | "),
        Category = NA_character_,
        BoundaryCategory = character_col(cumulative, "BoundaryCategory"),
        BoundaryOrder = numeric_col(cumulative, "BoundaryOrder"),
        CategorySet = category_set,
        Direction = direction,
        ValueName = "CumulativeProbability",
        Value = numeric_col(cumulative, "CumulativeProbability"),
        DisplayedByDefault = direction == cumulative_direction,
        Model = character_col(cumulative, "Model"),
        Slope = numeric_col(cumulative, "Slope"),
        stringsAsFactors = FALSE
      )
    }
    if (nrow(ogive) > 0L && all(c("Theta", "Information", "CurveGroup") %in% names(ogive))) {
      rows[[length(rows) + 1L]] <- data.frame(
        PlotType = "information",
        Panel = "Total information",
        CurveGroup = character_col(ogive, "CurveGroup"),
        Theta = numeric_col(ogive, "Theta"),
        Series = character_col(ogive, "CurveGroup"),
        Category = NA_character_,
        BoundaryCategory = NA_character_,
        BoundaryOrder = NA_real_,
        CategorySet = NA_character_,
        Direction = NA_character_,
        ValueName = "Information",
        Value = numeric_col(ogive, "Information"),
        DisplayedByDefault = TRUE,
        Model = character_col(ogive, "Model"),
        Slope = numeric_col(ogive, "Slope"),
        stringsAsFactors = FALSE
      )
    }
    if (nrow(cat_info) > 0L &&
        all(c("Theta", "CategoryInformation", "Category", "CurveGroup") %in% names(cat_info))) {
      rows[[length(rows) + 1L]] <- data.frame(
        PlotType = "category_information",
        Panel = "Category-specific information",
        CurveGroup = character_col(cat_info, "CurveGroup"),
        Theta = numeric_col(cat_info, "Theta"),
        Series = paste(character_col(cat_info, "CurveGroup"), character_col(cat_info, "Category"), sep = " | Cat "),
        Category = character_col(cat_info, "Category"),
        BoundaryCategory = NA_character_,
        BoundaryOrder = NA_real_,
        CategorySet = character_col(cat_info, "Category"),
        Direction = NA_character_,
        ValueName = "CategoryInformation",
        Value = numeric_col(cat_info, "CategoryInformation"),
        DisplayedByDefault = TRUE,
        Model = character_col(cat_info, "Model"),
        Slope = numeric_col(cat_info, "Slope"),
        stringsAsFactors = FALSE
      )
    }
    if (length(rows) == 0L) return(data.frame())
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }
  curve_long <- make_curve_long()
  curve_summary <- if (is.data.frame(curve_long) && nrow(curve_long) > 0L) {
    curve_long |>
      dplyr::group_by(.data$PlotType, .data$Panel, .data$ValueName) |>
      dplyr::summarise(
        Rows = dplyr::n(),
        Series = dplyr::n_distinct(.data$Series),
        CurveGroups = dplyr::n_distinct(.data$CurveGroup),
        DisplayedRows = sum(.data$DisplayedByDefault %in% TRUE, na.rm = TRUE),
        .groups = "drop"
      ) |>
      as.data.frame(stringsAsFactors = FALSE)
  } else {
    data.frame(
      PlotType = character(),
      Panel = character(),
      ValueName = character(),
      Rows = integer(),
      Series = integer(),
      CurveGroups = integer(),
      DisplayedRows = integer(),
      stringsAsFactors = FALSE
    )
  }
  make_curve_style <- function(curve_long) {
    if (!is.data.frame(curve_long) || nrow(curve_long) == 0L || !"Series" %in% names(curve_long)) {
      return(data.frame())
    }
    keys <- unique(as.character(curve_long$Series))
    keys <- keys[!is.na(keys) & nzchar(keys)]
    if (length(keys) == 0L) return(data.frame())
    cols <- line_palette(keys)
    ltys <- line_types(keys)
    data.frame(
      Series = keys,
      Colour = unname(cols[keys]),
      LineType = unname(as.integer(ltys[keys])),
      Preset = style$name,
      stringsAsFactors = FALSE
    )
  }
  curve_style <- make_curve_style(curve_long)
  draw_ogive_panel <- function(main_title) {
    if (nrow(ogive) == 0 || !all(c("Theta", "ExpectedScore", "CurveGroup") %in% names(ogive))) {
      stop("No expected-ogive data available.")
    }
    groups <- unique(as.character(ogive$CurveGroup))
    cols <- line_palette(groups)
    ltys <- line_types(groups)
    graphics::plot(
      x = finite_range(ogive$Theta),
      y = finite_range(ogive$ExpectedScore),
      type = "n",
      xlab = "Theta / Logit",
      ylab = "Expected score",
      main = main_title
    )
    graphics::grid(col = style$grid)
    for (i in seq_along(groups)) {
      sub <- ogive[ogive$CurveGroup == groups[i], , drop = FALSE]
      graphics::lines(sub$Theta, sub$ExpectedScore, col = cols[groups[i]], lwd = 2, lty = ltys[groups[i]])
    }
    if (length(groups) <= 8L) {
      graphics::legend("topleft", legend = groups, col = cols[groups], lty = ltys[groups], lwd = 2, bty = "n", cex = 0.8)
    }
  }
  draw_ccc_panel <- function(main_title) {
    if (nrow(probs) == 0 || !all(c("Theta", "Probability", "Category", "CurveGroup") %in% names(probs))) {
      stop("No category-curve data available.")
    }
    plot_tbl <- probs
    plot_tbl$Trace <- paste(plot_tbl$CurveGroup, plot_tbl$Category, sep = " | Cat ")
    traces <- unique(plot_tbl$Trace)
    cols <- line_palette(traces)
    ltys <- line_types(traces)
    graphics::plot(
      x = finite_range(plot_tbl$Theta),
      y = c(0, 1),
      type = "n",
      xlab = "Theta / Logit",
      ylab = "Probability",
      main = main_title
    )
    graphics::grid(col = style$grid)
    for (i in seq_along(traces)) {
      sub <- plot_tbl[plot_tbl$Trace == traces[i], , drop = FALSE]
      graphics::lines(sub$Theta, sub$Probability, col = cols[traces[i]], lwd = 1.4, lty = ltys[traces[i]])
    }
  }
  draw_cumulative_panel <- function(main_title) {
    if (nrow(cumulative) == 0 ||
        !all(c("Theta", "CumulativeProbability", "BoundaryCategory", "CurveGroup", "Direction") %in% names(cumulative))) {
      stop("No cumulative probability data available.")
    }
    cumulative_plot <- cumulative[cumulative$Direction == cumulative_direction, , drop = FALSE]
    if (nrow(cumulative_plot) == 0) {
      stop("No cumulative probability data available for `cumulative_direction = ", cumulative_direction, "`.")
    }
    if (!"CategorySet" %in% names(cumulative_plot)) {
      cumulative_plot$CategorySet <- as.character(cumulative_plot$BoundaryCategory)
    }
    cumulative_plot$Trace <- paste(cumulative_plot$CurveGroup, cumulative_plot$CategorySet, sep = " | ")
    traces <- unique(cumulative_plot$Trace)
    cols <- line_palette(traces)
    ltys <- line_types(traces)
    graphics::plot(
      x = finite_range(cumulative_plot$Theta),
      y = c(0, 1),
      type = "n",
      xlab = "Theta / Logit",
      ylab = "Cumulative probability",
      main = main_title
    )
    graphics::grid(col = style$grid)
    graphics::abline(h = 0.5, lty = 3, col = style$neutral)
    if (nrow(boundary_line_tbl) > 0L) {
      boundary_x <- unique(suppressWarnings(as.numeric(boundary_line_tbl$ThurstonianThreshold)))
      boundary_x <- boundary_x[is.finite(boundary_x)]
      if (length(boundary_x) > 0L) {
        graphics::abline(v = boundary_x, lty = 3, col = if (identical(style$name, "monochrome")) "gray55" else "gray75")
      }
    }
    for (i in seq_along(traces)) {
      sub <- cumulative_plot[cumulative_plot$Trace == traces[i], , drop = FALSE]
      graphics::lines(sub$Theta, sub$CumulativeProbability, col = cols[traces[i]], lwd = 1.4, lty = ltys[traces[i]])
    }
  }
  draw_information_panel <- function(main_title) {
    if (nrow(ogive) == 0 || !all(c("Theta", "Information", "CurveGroup") %in% names(ogive))) {
      stop("No total information data available.")
    }
    groups <- unique(as.character(ogive$CurveGroup))
    cols <- line_palette(groups)
    ltys <- line_types(groups)
    graphics::plot(
      x = finite_range(ogive$Theta),
      y = finite_range(ogive$Information),
      type = "n",
      xlab = "Theta / Logit",
      ylab = "Information",
      main = main_title
    )
    graphics::grid(col = style$grid)
    for (i in seq_along(groups)) {
      sub <- ogive[ogive$CurveGroup == groups[i], , drop = FALSE]
      graphics::lines(sub$Theta, sub$Information, col = cols[groups[i]], lwd = 2, lty = ltys[groups[i]])
    }
    if (length(groups) <= 8L) {
      graphics::legend("topleft", legend = groups, col = cols[groups], lty = ltys[groups], lwd = 2, bty = "n", cex = 0.8)
    }
  }
  draw_category_information_panel <- function(main_title) {
    if (nrow(cat_info) == 0 ||
        !all(c("Theta", "CategoryInformation", "Category", "CurveGroup") %in% names(cat_info))) {
      stop("No category-specific information data available.")
    }
    plot_tbl <- cat_info
    plot_tbl$Trace <- paste(plot_tbl$CurveGroup, plot_tbl$Category, sep = " | Cat ")
    traces <- unique(plot_tbl$Trace)
    cols <- line_palette(traces)
    ltys <- line_types(traces)
    graphics::plot(
      x = finite_range(plot_tbl$Theta),
      y = finite_range(plot_tbl$CategoryInformation),
      type = "n",
      xlab = "Theta / Logit",
      ylab = "Category information",
      main = main_title
    )
    graphics::grid(col = style$grid)
    for (i in seq_along(traces)) {
      sub <- plot_tbl[plot_tbl$Trace == traces[i], , drop = FALSE]
      graphics::lines(sub$Theta, sub$CategoryInformation, col = cols[traces[i]], lwd = 1.4, lty = ltys[traces[i]])
    }
  }

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    if (type == "overview") {
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.2, 1.2), oma = c(0, 0, 2.2, 0))
      draw_ccc_panel("Category probability")
      draw_cumulative_panel("Cumulative probability")
      draw_information_panel("Total information")
      draw_category_information_panel("Category information")
      graphics::mtext(
        if (is.null(main)) "Category curve overview" else as.character(main[1]),
        outer = TRUE,
        cex = 1.1,
        font = 2
      )
    } else if (type == "ogive") {
      draw_ogive_panel(if (is.null(main)) "Expected-score ogive" else as.character(main[1]))
    } else if (type == "ccc") {
      draw_ccc_panel(if (is.null(main)) "Category characteristic curves" else as.character(main[1]))
    } else if (type == "cumulative") {
      draw_cumulative_panel(if (is.null(main)) "Cumulative probability curves" else as.character(main[1]))
    } else if (type == "information") {
      draw_information_panel(if (is.null(main)) "Total information curves" else as.character(main[1]))
    } else {
      draw_category_information_panel(if (is.null(main)) "Category-specific information curves" else as.character(main[1]))
    }
  }

  new_mfrm_plot_data(
    "category_curves",
    list(
      plot = type,
      expected_ogive = ogive,
      probabilities = probs,
      cumulative_probabilities = cumulative,
      cumulative_boundaries = cumulative_boundaries,
      cumulative_direction = cumulative_direction,
      category_information = cat_info,
      overview_panels = overview_panels,
      plot_long = curve_long,
      plot_annotations = plot_annotations,
      curve_summary = curve_summary,
      curve_style = curve_style,
      boundary_lines = boundary_line_tbl,
      plot_settings = plot_settings,
      preset = style$name,
      title = switch(
        type,
        overview = "Category curve overview",
        ogive = "Expected-score ogive",
        ccc = if (requested_type %in% c("category_probability", "category_probabilities", "conditional_probability", "conditional_probabilities", "probability", "probabilities")) {
          "Category probability curves"
        } else {
          "Category characteristic curves"
        },
        cumulative = "Cumulative probability curves",
        information = "Total information curves",
        category_information = "Category-specific information curves"
      ),
      subtitle = switch(
        type,
        overview = "Category probabilities, cumulative probabilities, total information, and category-specific information",
        ogive = "Expected score across theta",
        ccc = "Category response probabilities conditional on theta",
        cumulative = "Modeled probability accumulated across ordered categories",
        information = "Total per-curve information; GPCM uses a^2 times score variance",
        category_information = "Category contributions sum to the total information at each theta"
      ),
      legend = new_plot_legend(
        label = switch(
          type,
          overview = "Category curve overview",
          ogive = "Expected score",
          ccc = "Category probability",
          cumulative = "Cumulative probability",
          information = "Information",
          category_information = "Category information contribution"
        ),
        role = switch(
          type,
          overview = "overview",
          ogive = "expected_score",
          ccc = "probability",
          cumulative = "cumulative_probability",
          information = "information",
          category_information = "category_information"
        ),
        aesthetic = "line",
        value = "curve_group_palette"
      ),
      reference_lines = reference_line_tbl
    )
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
                               type = c("graph_expected", "score_residuals", "obs_probability", "score_se"),
                               draw = TRUE,
                               main = NULL,
                               palette = NULL) {
  type <- match.arg(
    tolower(as.character(type[1])),
    c("graph_expected", "score_residuals", "obs_probability", "score_se")
  )
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

  if (type == "score_se") {
    se_col <- if ("ScoreSideSE" %in% names(score_tbl) &&
                  any(is.finite(suppressWarnings(as.numeric(score_tbl$ScoreSideSE))))) {
      "ScoreSideSE"
    } else if ("ExpectedScoreSE" %in% names(score_tbl) &&
               any(is.finite(suppressWarnings(as.numeric(score_tbl$ExpectedScoreSE))))) {
      "ExpectedScoreSE"
    } else {
      NA_character_
    }
    if (is.na(se_col) || !nzchar(se_col)) {
      stop("No finite scorefile SE column is available. Use `score_se_method` to request SE output.", call. = FALSE)
    }
    vals <- suppressWarnings(as.numeric(score_tbl[[se_col]]))
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) stop("No finite scorefile SE values available.")
    if (isTRUE(draw)) {
      graphics::hist(
        x = vals,
        breaks = "FD",
        col = "#fdd49e",
        border = "white",
        main = if (is.null(main)) paste0(se_col, " distribution") else as.character(main[1]),
        xlab = se_col,
        ylab = "Count"
      )
    }
    return(invisible(new_mfrm_plot_data(
      "output_bundle",
      list(plot = "score_se", se_column = se_col, values = vals)
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
                                     type = c("dashboard", "quality_flags", "row_review", "category_counts", "score_support", "facet_category_usage", "facet_response_patterns", "score_map", "missing_rows"),
                                     draw = TRUE,
                                     main = NULL,
                                     palette = NULL,
                                     label_angle = 45,
                                     preset = c("standard", "publication", "compact", "monochrome"),
                                     top_n = 30L) {
  type <- match.arg(tolower(as.character(type[1])), c("dashboard", "quality_flags", "row_review", "category_counts", "score_support", "facet_category_usage", "facet_response_patterns", "score_map", "missing_rows"))
  style <- resolve_plot_preset(preset)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      row_review = style$accent_primary,
      category = style$accent_tertiary,
      missing = style$accent_secondary,
      zero_internal = style$fail,
      zero_boundary = style$warn,
      sparse = style$accent_secondary,
      ok = style$accent_tertiary,
      high = style$fail,
      review = style$warn,
      recoded = style$accent_secondary
    )
  )
  row_tbl <- as.data.frame(x$row_review %||% data.frame(), stringsAsFactors = FALSE)
  quality_tbl <- as.data.frame(x$quality_overview %||% data.frame(), stringsAsFactors = FALSE)
  cat_tbl <- as.data.frame(x$category_counts %||% data.frame(), stringsAsFactors = FALSE)
  support_tbl <- as.data.frame(x$score_support_review %||% cat_tbl, stringsAsFactors = FALSE)
  usage_summary <- as.data.frame(x$category_usage_summary %||% data.frame(), stringsAsFactors = FALSE)
  pattern_tbl <- as.data.frame(x$facet_response_patterns %||% data.frame(), stringsAsFactors = FALSE)
  flags_tbl <- as.data.frame(x$quality_flags %||% data.frame(), stringsAsFactors = FALSE)
  score_map_tbl <- as.data.frame(x$score_map %||% data.frame(), stringsAsFactors = FALSE)
  sum_tbl <- as.data.frame(x$summary %||% data.frame(), stringsAsFactors = FALSE)
  top_n <- suppressWarnings(as.integer(top_n[1]))
  if (!is.finite(top_n) || top_n <= 0L) top_n <- 30L

  build_missing_rows_table <- function(sum_tbl) {
    if (nrow(sum_tbl) == 0) {
      return(data.frame(Field = character(0), Count = numeric(0), stringsAsFactors = FALSE))
    }
    row_cols <- grep("Rows$", names(sum_tbl), value = TRUE)
    if (length(row_cols) == 0) {
      return(data.frame(Field = character(0), Count = numeric(0), stringsAsFactors = FALSE))
    }
    data.frame(
      Field = row_cols,
      Count = suppressWarnings(as.numeric(sum_tbl[1, row_cols, drop = TRUE])),
      stringsAsFactors = FALSE
    )
  }

  prepare_score_support_plot <- function(cat_tbl, support_tbl) {
    plot_tbl <- if (nrow(support_tbl) > 0) support_tbl else cat_tbl
    if (nrow(plot_tbl) == 0 || !"Score" %in% names(plot_tbl)) return(plot_tbl)
    score_order <- order(suppressWarnings(as.numeric(plot_tbl$Score)), na.last = TRUE)
    plot_tbl[score_order, , drop = FALSE]
  }

  score_support_colors <- function(plot_tbl) {
    cols <- rep(pal["category"], nrow(plot_tbl))
    if ("ZeroCount" %in% names(plot_tbl)) {
      zero <- as.logical(plot_tbl$ZeroCount)
      zero[is.na(zero)] <- FALSE
      unused_type <- as.character(plot_tbl$UnusedCategoryType %||% rep("none", nrow(plot_tbl)))
      unused_type[is.na(unused_type)] <- "none"
      cols[zero & unused_type == "internal"] <- pal["zero_internal"]
      cols[zero & unused_type != "internal"] <- pal["zero_boundary"]
    }
    cols
  }

  prepare_facet_usage_plot <- function(usage_summary, top_n) {
    if (nrow(usage_summary) == 0) return(usage_summary)
    usage_summary |>
      dplyr::arrange(
        dplyr::desc(.data$IntermediateZeroCategories),
        dplyr::desc(.data$ZeroCategories),
        dplyr::desc(.data$SparseCategories),
        .data$Facet,
        .data$Level
      ) |>
      dplyr::slice_head(n = top_n) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  facet_usage_colors <- function(plot_tbl) {
    if (nrow(plot_tbl) == 0) return(character(0))
    ifelse(
      plot_tbl$IntermediateZeroCategories > 0,
      pal["zero_internal"],
      ifelse(plot_tbl$ZeroCategories > 0, pal["zero_boundary"],
             ifelse(plot_tbl$SparseCategories > 0, pal["sparse"], pal["ok"]))
    )
  }

  prepare_quality_flags_plot <- function(flags_tbl) {
    if (nrow(flags_tbl) == 0L || !"Area" %in% names(flags_tbl)) {
      return(data.frame(
        Area = character(0),
        Flags = integer(0),
        HighSeverityFlags = integer(0),
        ReviewFlags = integer(0),
        TotalReferencedCount = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    if (!"Severity" %in% names(flags_tbl)) flags_tbl$Severity <- "review"
    if (!"Count" %in% names(flags_tbl)) flags_tbl$Count <- NA_real_
    out <- flags_tbl |>
      dplyr::group_by(.data$Area) |>
      dplyr::summarise(
        Flags = dplyr::n(),
        HighSeverityFlags = sum(tolower(as.character(.data$Severity)) %in% "high", na.rm = TRUE),
        ReviewFlags = sum(tolower(as.character(.data$Severity)) %in% "review", na.rm = TRUE),
        TotalReferencedCount = sum(suppressWarnings(as.numeric(.data$Count)), na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(.data$HighSeverityFlags), dplyr::desc(.data$Flags), .data$Area) |>
      as.data.frame(stringsAsFactors = FALSE)
    row.names(out) <- NULL
    out
  }

  quality_flag_colors <- function(plot_tbl) {
    if (nrow(plot_tbl) == 0L) return(character(0))
    ifelse(plot_tbl$HighSeverityFlags > 0, pal["high"], pal["review"])
  }

  prepare_pattern_plot <- function(pattern_tbl, top_n) {
    if (nrow(pattern_tbl) == 0L) return(pattern_tbl)
    pattern_tbl |>
      dplyr::arrange(
        match(.data$PatternStatus, c("high", "review", "ok")),
        dplyr::desc(.data$DominantPercent),
        .data$Facet,
        .data$Level
      ) |>
      dplyr::slice_head(n = top_n) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  pattern_colors <- function(plot_tbl) {
    if (nrow(plot_tbl) == 0L) return(character(0))
    status <- as.character(plot_tbl$PatternStatus %||% rep("ok", nrow(plot_tbl)))
    ifelse(status == "high", pal["high"],
           ifelse(status == "review", pal["review"], pal["ok"]))
  }

  prepare_score_map_plot <- function(score_map_tbl, cat_tbl) {
    if (nrow(score_map_tbl) > 0L &&
        all(c("OriginalScore", "InternalScore") %in% names(score_map_tbl))) {
      out <- score_map_tbl[, c("OriginalScore", "InternalScore"), drop = FALSE]
    } else if (nrow(cat_tbl) > 0L && "Score" %in% names(cat_tbl)) {
      out <- data.frame(
        OriginalScore = cat_tbl$Score,
        InternalScore = cat_tbl$Score,
        stringsAsFactors = FALSE
      )
    } else {
      return(data.frame(
        OriginalScore = character(0),
        InternalScore = character(0),
        OriginalNumeric = numeric(0),
        InternalNumeric = numeric(0),
        MappingStatus = character(0),
        stringsAsFactors = FALSE
      ))
    }
    out$OriginalNumeric <- suppressWarnings(as.numeric(out$OriginalScore))
    out$InternalNumeric <- suppressWarnings(as.numeric(out$InternalScore))
    out$MappingStatus <- ifelse(
      as.character(out$OriginalScore) == as.character(out$InternalScore),
      "identity",
      "recoded"
    )
    order_idx <- order(out$OriginalNumeric, as.character(out$OriginalScore), na.last = TRUE)
    out <- out[order_idx, , drop = FALSE]
    row.names(out) <- NULL
    out
  }

  if (type == "quality_flags") {
    plot_tbl <- prepare_quality_flags_plot(flags_tbl)
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      if (nrow(plot_tbl) > 0L) {
        barplot_rot45(
          height = suppressWarnings(as.numeric(plot_tbl$Flags)),
          labels = as.character(plot_tbl$Area),
          col = quality_flag_colors(plot_tbl),
          main = if (is.null(main)) "Data quality flags" else as.character(main[1]),
          ylab = "Flags",
          label_angle = label_angle,
          mar_bottom = 7.8
        )
      } else {
        graphics::plot.new()
        graphics::title(if (is.null(main)) "Data quality flags" else as.character(main[1]))
        graphics::text(0.5, 0.5, "No priority QC flags")
      }
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(
        plot = "quality_flags",
        table = plot_tbl,
        quality_flags = flags_tbl,
        quality_overview = quality_tbl,
        preset = style$name
      )
    )))
  }

  if (type == "row_review") {
    if (nrow(row_tbl) == 0 || !all(c("Status", "N") %in% names(row_tbl))) {
      stop("Row-review table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    vals <- suppressWarnings(as.numeric(row_tbl$N))
    labels <- as.character(row_tbl$Status)
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pal["row_review"],
        main = if (is.null(main)) "Row-review status counts" else as.character(main[1]),
        ylab = "Rows",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(plot = "row_review", table = row_tbl, preset = style$name)
    )))
  }

  if (type == "facet_category_usage") {
    if (nrow(usage_summary) == 0 ||
        !all(c("Facet", "Level", "IssueCategories") %in% names(usage_summary))) {
      stop("Facet-category usage summary is not available. Run data_quality_report() first.", call. = FALSE)
    }
    plot_tbl <- prepare_facet_usage_plot(usage_summary, top_n = top_n)
    vals <- suppressWarnings(as.numeric(plot_tbl$IssueCategories))
    labels <- paste(plot_tbl$Facet, plot_tbl$Level, sep = ": ")
    cols <- facet_usage_colors(plot_tbl)
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      barplot_rot45(
        height = vals,
        labels = labels,
        col = cols,
        main = if (is.null(main)) "Facet-level category usage issues" else as.character(main[1]),
        ylab = "Zero or sparse categories",
        label_angle = label_angle,
        mar_bottom = 9.0
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(plot = "facet_category_usage", table = plot_tbl, preset = style$name, top_n = top_n)
    )))
  }

  if (type == "score_map") {
    plot_tbl <- prepare_score_map_plot(score_map_tbl, cat_tbl)
    if (nrow(plot_tbl) == 0L) {
      stop("Score-map table is not available. Run data_quality_report() from a fitted object with score information.", call. = FALSE)
    }
    cols <- ifelse(plot_tbl$MappingStatus == "recoded", pal["recoded"], pal["ok"])
    y <- plot_tbl$InternalNumeric
    if (any(!is.finite(y))) y <- seq_len(nrow(plot_tbl))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      graphics::plot(
        seq_len(nrow(plot_tbl)),
        y,
        type = "b",
        pch = 16,
        col = cols,
        xaxt = "n",
        xlab = "Original score",
        ylab = "Internal score",
        main = if (is.null(main)) "Score-map review" else as.character(main[1])
      )
      graphics::axis(1, at = seq_len(nrow(plot_tbl)), labels = as.character(plot_tbl$OriginalScore), las = 2)
      if (any(plot_tbl$MappingStatus == "recoded", na.rm = TRUE)) {
        graphics::legend(
          "topleft",
          legend = c("Identity", "Recoded"),
          col = c(pal["ok"], pal["recoded"]),
          pch = 16,
          bty = "n"
        )
      }
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(plot = "score_map", table = plot_tbl, preset = style$name)
    )))
  }

  if (type == "facet_response_patterns") {
    if (nrow(pattern_tbl) == 0 ||
        !all(c("Facet", "Level", "DominantPercent") %in% names(pattern_tbl))) {
      stop("Facet-response pattern table is not available. Run data_quality_report() first.", call. = FALSE)
    }
    plot_tbl <- prepare_pattern_plot(pattern_tbl, top_n = top_n)
    vals <- 100 * suppressWarnings(as.numeric(plot_tbl$DominantPercent))
    labels <- paste(plot_tbl$Facet, plot_tbl$Level, sep = ": ")
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      barplot_rot45(
        height = vals,
        labels = labels,
        col = pattern_colors(plot_tbl),
        main = if (is.null(main)) "Facet response-pattern dominance" else as.character(main[1]),
        ylab = "Dominant score (%)",
        label_angle = label_angle,
        mar_bottom = 9.0
      )
      graphics::abline(h = 95, lty = 2, col = "gray50")
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(plot = "facet_response_patterns", table = plot_tbl, preset = style$name, top_n = top_n)
    )))
  }

  if (type == "dashboard") {
    support_plot_tbl <- prepare_score_support_plot(cat_tbl, support_tbl)
    facet_plot_tbl <- prepare_facet_usage_plot(usage_summary, top_n = top_n)
    missing_tbl <- build_missing_rows_table(sum_tbl)

    if (isTRUE(draw)) {
      apply_plot_preset(style)
      graphics::par(mfrow = c(2, 2))
      if (nrow(row_tbl) > 0 && all(c("Status", "N") %in% names(row_tbl))) {
        barplot_rot45(
          height = suppressWarnings(as.numeric(row_tbl$N)),
          labels = as.character(row_tbl$Status),
          col = pal["row_review"],
          main = if (is.null(main)) "Data quality dashboard: row review" else as.character(main[1]),
          ylab = "Rows",
          label_angle = label_angle,
          mar_bottom = 7.2
        )
      } else {
        graphics::plot.new()
        graphics::title("Row review")
        graphics::text(0.5, 0.5, "No row-review table")
      }

      if (nrow(support_plot_tbl) > 0 && all(c("Score", "Count") %in% names(support_plot_tbl))) {
        barplot_rot45(
          height = suppressWarnings(as.numeric(support_plot_tbl$Count)),
          labels = as.character(support_plot_tbl$Score),
          col = score_support_colors(support_plot_tbl),
          main = "Score support",
          ylab = "Count",
          label_angle = label_angle,
          mar_bottom = 6.8
        )
      } else {
        graphics::plot.new()
        graphics::title("Score support")
        graphics::text(0.5, 0.5, "No score-support table")
      }

      if (nrow(facet_plot_tbl) > 0 && all(c("Facet", "Level", "IssueCategories") %in% names(facet_plot_tbl))) {
        barplot_rot45(
          height = suppressWarnings(as.numeric(facet_plot_tbl$IssueCategories)),
          labels = paste(facet_plot_tbl$Facet, facet_plot_tbl$Level, sep = ": "),
          col = facet_usage_colors(facet_plot_tbl),
          main = "Facet category use",
          ylab = "Issues",
          label_angle = label_angle,
          mar_bottom = 8.0
        )
      } else {
        graphics::plot.new()
        graphics::title("Facet category use")
        graphics::text(0.5, 0.5, "No facet category summary")
      }

      if (nrow(missing_tbl) > 0) {
        barplot_rot45(
          height = suppressWarnings(as.numeric(missing_tbl$Count)),
          labels = as.character(missing_tbl$Field),
          col = pal["missing"],
          main = "Missing/invalid rows",
          ylab = "Rows",
          label_angle = label_angle,
          mar_bottom = 8.0
        )
      } else {
        graphics::plot.new()
        graphics::title("Missing/invalid rows")
        graphics::text(0.5, 0.5, "No row-count fields")
      }
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(
        plot = "dashboard",
        row_review = row_tbl,
        quality_overview = quality_tbl,
        quality_flags = flags_tbl,
        score_map = prepare_score_map_plot(score_map_tbl, cat_tbl),
        score_support = support_plot_tbl,
        facet_category_usage = facet_plot_tbl,
        facet_response_patterns = prepare_pattern_plot(pattern_tbl, top_n = top_n),
        missing_rows = missing_tbl,
        preset = style$name,
        top_n = top_n
      )
    )))
  }

  if (type %in% c("category_counts", "score_support")) {
    if (nrow(cat_tbl) == 0 || !all(c("Score", "Count") %in% names(cat_tbl))) {
      stop("Category-count table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    plot_tbl <- if (type == "score_support") prepare_score_support_plot(cat_tbl, support_tbl) else prepare_score_support_plot(cat_tbl, cat_tbl)
    vals <- suppressWarnings(as.numeric(plot_tbl$Count))
    labels <- as.character(plot_tbl$Score)
    cols <- score_support_colors(plot_tbl)
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      barplot_rot45(
        height = vals,
        labels = labels,
        col = cols,
        main = if (is.null(main)) {
          if (type == "score_support") "Score-support category review" else "Observed category counts"
        } else {
          as.character(main[1])
        },
        ylab = "Count",
        label_angle = label_angle,
        mar_bottom = 7.8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "data_quality",
      list(plot = type, table = plot_tbl, preset = style$name)
    )))
  }

  if (nrow(sum_tbl) == 0) stop("Summary table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
  missing_tbl <- build_missing_rows_table(sum_tbl)
  if (nrow(missing_tbl) == 0) stop("No row-count columns found in summary table.")
  vals <- suppressWarnings(as.numeric(missing_tbl$Count))
  labels <- missing_tbl$Field
  if (isTRUE(draw)) {
    apply_plot_preset(style)
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
    list(
      plot = "missing_rows",
      table = missing_tbl,
      preset = style$name
    )
  ))
}

draw_fit_measures_bundle <- function(x,
                                     type = c("status", "infit_outfit", "measure_ci", "df_sensitivity"),
                                     draw = TRUE,
                                     main = NULL,
                                     palette = NULL,
                                     label_angle = 45,
                                     preset = c("standard", "publication", "compact", "monochrome"),
                                     ci_level = NULL,
                                     top_n = 30L) {
  type <- match.arg(tolower(as.character(type[1])), c("status", "infit_outfit", "measure_ci", "df_sensitivity"))
  style <- resolve_plot_preset(preset)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      underfit = style$warn,
      overfit = style$accent_secondary,
      mixed = style$fail,
      within_band = style$success,
      not_available = style$neutral,
      flag_changed_by_df = style$fail,
      large_zstd_shift = style$warn,
      df_convention_difference = style$accent_secondary,
      small_zstd_shift = style$accent_primary,
      same_or_rounding = style$neutral,
      reference = style$neutral,
      ci = style$neutral,
      point = style$accent_primary
    )
  )
  tbl <- as.data.frame(x$table %||% data.frame(), stringsAsFactors = FALSE)
  if (!identical(type, "df_sensitivity") &&
      (nrow(tbl) == 0 || !"FitStatus" %in% names(tbl))) {
    stop("Fit-measures table is not available. Run fit_measures_table() first.", call. = FALSE)
  }

  if (type == "status") {
    counts <- as.data.frame(table(tbl$FitStatus), stringsAsFactors = FALSE)
    names(counts) <- c("FitStatus", "Rows")
    preferred <- c("underfit", "overfit", "mixed", "within_band", "not_available")
    counts <- counts[order(match(counts$FitStatus, preferred)), , drop = FALSE]
    cols <- unname(pal[match(counts$FitStatus, names(pal))])
    cols[is.na(cols)] <- pal["not_available"]
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      barplot_rot45(
        height = suppressWarnings(as.numeric(counts$Rows)),
        labels = as.character(counts$FitStatus),
        col = cols,
        main = if (is.null(main)) "Fit-measure status counts" else as.character(main[1]),
        ylab = "Rows",
        label_angle = label_angle,
        mar_bottom = 8.2
      )
    }
    return(invisible(new_mfrm_plot_data(
      "fit_measures",
      list(plot = "status", table = counts, preset = style$name)
    )))
  }

  if (type == "measure_ci") {
    if (!all(c("Measure", "SE") %in% names(tbl))) {
      stop("Measure and SE columns are not available in the fit-measures table.", call. = FALSE)
    }
    active_ci <- ci_level
    if (is.null(active_ci)) {
      active_ci <- x$settings$ci_level %||% if ("CI_Level" %in% names(tbl)) tbl$CI_Level[1] else 0.95
    }
    active_ci <- fit_measure_validate_ci_level(active_ci)
    measure <- suppressWarnings(as.numeric(tbl$Measure))
    se <- suppressWarnings(as.numeric(tbl$SE))
    z_ci <- stats::qnorm(1 - (1 - active_ci) / 2)
    tbl$CI_Lower <- ifelse(is.finite(measure) & is.finite(se) & se >= 0, measure - z_ci * se, NA_real_)
    tbl$CI_Upper <- ifelse(is.finite(measure) & is.finite(se) & se >= 0, measure + z_ci * se, NA_real_)
    tbl$CI_Level <- active_ci
    ok <- is.finite(measure) & is.finite(tbl$CI_Lower) & is.finite(tbl$CI_Upper)
    if (!any(ok)) stop("No finite measure confidence intervals are available for plotting.", call. = FALSE)
    labels <- paste(tbl$Facet, tbl$Level, sep = ": ")
    y <- seq_len(nrow(tbl))
    cols <- unname(pal[match(tbl$FitStatus, names(pal))])
    cols[is.na(cols)] <- pal["point"]
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      xr <- range(c(tbl$CI_Lower[ok], tbl$CI_Upper[ok], measure[ok]), finite = TRUE)
      graphics::plot(
        x = measure[ok],
        y = y[ok],
        pch = 21,
        bg = cols[ok],
        col = style$background,
        xlim = xr,
        ylim = rev(range(y[ok])),
        yaxt = "n",
        xlab = "Measure (logits)",
        ylab = "",
        main = if (is.null(main)) {
          paste0("Fit-measure estimates with ", round(100 * active_ci), "% CI")
        } else {
          as.character(main[1])
        }
      )
      graphics::segments(
        x0 = tbl$CI_Lower[ok],
        y0 = y[ok],
        x1 = tbl$CI_Upper[ok],
        y1 = y[ok],
        col = pal["ci"],
        lwd = 1.6
      )
      graphics::points(
        x = measure[ok],
        y = y[ok],
        pch = 21,
        bg = cols[ok],
        col = style$background
      )
      graphics::axis(2, at = y[ok], labels = labels[ok], las = 2, cex.axis = 0.76)
      graphics::abline(v = 0, lty = 3, col = pal["reference"])
    }
    return(invisible(new_mfrm_plot_data(
      "fit_measures",
      list(
        plot = "measure_ci",
        table = tbl,
        ci_level = active_ci,
        preset = style$name
      )
    )))
  }

  if (type == "df_sensitivity") {
    sens <- as.data.frame(x$df_sensitivity %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(sens) == 0 || !"MaxAbsZSTDDiff_FACETS_vs_ENGINE" %in% names(sens)) {
      stop("df_sensitivity rows are not available. Run fit_measures_table(..., fit_df_method = \"both\").", call. = FALSE)
    }
    z_shift <- suppressWarnings(as.numeric(sens$MaxAbsZSTDDiff_FACETS_vs_ENGINE))
    ok <- is.finite(z_shift)
    if (!any(ok)) {
      stop("No finite engine-vs-FACETS-style ZSTD differences are available for plotting.", call. = FALSE)
    }
    sens <- sens[ok, , drop = FALSE]
    z_shift <- z_shift[ok]
    flag <- sens$FlagChangedByDf %in% TRUE
    status <- as.character(sens$DfSensitivityStatus %||% "not_available")
    ord <- order(-as.integer(flag), -z_shift, sens$Facet, sens$Level)
    sens <- sens[ord, , drop = FALSE]
    z_shift <- z_shift[ord]
    top_n_num <- suppressWarnings(as.numeric(top_n[1]))
    if (is.finite(top_n_num)) {
      keep <- seq_len(min(nrow(sens), max(1L, as.integer(top_n_num))))
      sens <- sens[keep, , drop = FALSE]
      z_shift <- z_shift[keep]
    }
    labels <- paste(sens$Facet, sens$Level, sep = ": ")
    status <- as.character(sens$DfSensitivityStatus %||% "not_available")
    cols <- unname(pal[match(status, names(pal))])
    cols[is.na(cols)] <- pal["not_available"]
    thresholds <- c(
      tolerance = suppressWarnings(as.numeric(x$settings$df_zstd_tolerance %||% NA_real_)),
      large_shift = suppressWarnings(as.numeric(x$settings$df_zstd_large_shift %||% NA_real_))
    )
    thresholds <- thresholds[is.finite(thresholds)]
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      ylim <- range(c(0, z_shift, thresholds), finite = TRUE)
      bars <- graphics::barplot(
        height = z_shift,
        names.arg = FALSE,
        col = cols,
        border = "white",
        ylim = ylim,
        ylab = "Max |ZSTD difference|",
        main = if (is.null(main)) "FACETS-style df sensitivity" else as.character(main[1])
      )
      if (length(thresholds) > 0L) {
        graphics::abline(h = thresholds, lty = c(3, 2)[seq_along(thresholds)], col = pal["reference"])
      }
      draw_rotated_x_labels(
        at = bars,
        labels = labels,
        srt = label_angle,
        cex = 0.78,
        line_offset = 0.085
      )
    }
    return(invisible(new_mfrm_plot_data(
      "fit_measures",
      list(
        plot = "df_sensitivity",
        table = sens,
        thresholds = thresholds,
        preset = style$name
      )
    )))
  }

  if (!all(c("Infit", "Outfit") %in% names(tbl))) {
    stop("Infit/Outfit columns are not available in the fit-measures table.", call. = FALSE)
  }
  xv <- suppressWarnings(as.numeric(tbl$Infit))
  yv <- suppressWarnings(as.numeric(tbl$Outfit))
  ok <- is.finite(xv) & is.finite(yv)
  if (!any(ok)) stop("No finite Infit/Outfit rows available for plotting.", call. = FALSE)
  band <- x$settings %||% list()
  lower <- suppressWarnings(as.numeric(band$lower %||% 0.5))
  upper <- suppressWarnings(as.numeric(band$upper %||% 1.5))
  cols <- unname(pal[match(tbl$FitStatus, names(pal))])
  cols[is.na(cols)] <- pal["not_available"]
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    xr <- range(c(xv[ok], lower, 1, upper), finite = TRUE)
    yr <- range(c(yv[ok], lower, 1, upper), finite = TRUE)
    graphics::plot(
      x = xv[ok],
      y = yv[ok],
      pch = 21,
      bg = cols[ok],
      col = "white",
      xlim = xr,
      ylim = yr,
      xlab = "Infit MnSq",
      ylab = "Outfit MnSq",
      main = if (is.null(main)) "Infit vs Outfit by fit status" else as.character(main[1])
    )
    graphics::abline(v = c(lower, 1, upper), lty = c(3, 2, 3), col = pal["reference"])
    graphics::abline(h = c(lower, 1, upper), lty = c(3, 2, 3), col = pal["reference"])
  }
  invisible(new_mfrm_plot_data(
    "fit_measures",
    list(
      plot = "infit_outfit",
      table = tbl,
      fit_range = c(lower = lower, upper = upper),
      preset = style$name
    )
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

draw_network_analysis_bundle <- function(x,
                                         type = c("centrality", "facet_summary", "network"),
                                         metric = NULL,
                                         top_n = 20,
                                         draw = TRUE,
                                         main = NULL,
                                         palette = NULL,
                                         label_angle = 45,
                                         preset = c("standard", "publication", "compact", "monochrome")) {
  type <- match.arg(tolower(as.character(type[1])), c("centrality", "facet_summary", "network"))
  if (identical(type, "network")) {
    sc <- x$source_connectivity %||% NULL
    if (is.null(sc)) {
      stop("Network source connectivity data are not available.", call. = FALSE)
    }
    return(draw_subset_connectivity_bundle(
      sc,
      type = "network",
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle,
      preset = preset
    ))
  }

  style <- resolve_plot_preset(preset)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      centrality = style$accent_primary,
      vulnerability = style$accent_secondary
    )
  )
  top_n <- max(1L, as.integer(top_n))

  if (identical(type, "centrality")) {
    tbl <- as.data.frame(x$node_metrics %||% data.frame(), stringsAsFactors = FALSE)
    metric <- as.character(metric %||% "Betweenness")[1L]
    numeric_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
    if (!metric %in% numeric_cols) {
      stop("`metric` must be a numeric node_metrics column: ",
           paste(numeric_cols, collapse = ", "), call. = FALSE)
    }
    tbl <- tbl |>
      dplyr::arrange(dplyr::desc(.data[[metric]]), .data$Facet, .data$Level) |>
      dplyr::slice_head(n = top_n)
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      labels <- paste(tbl$Facet, tbl$Level, sep = ":")
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl[[metric]])),
        labels = labels,
        col = pal["centrality"],
        main = main %||% paste("Design-network", metric),
        ylab = metric,
        label_angle = label_angle,
        mar_bottom = 9
      )
    }
    return(invisible(new_mfrm_plot_data(
      "network_analysis",
      list(
        plot = "centrality",
        table = tbl,
        metric = metric,
        title = main %||% paste("Design-network", metric),
        subtitle = "Node-level graph metric; interpret as design-link dependence",
        legend = new_plot_legend(metric, "node_metric", "bar", pal["centrality"]),
        reference_lines = new_reference_lines(),
        preset = style$name
      )
    )))
  }

  tbl <- as.data.frame(x$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
  metric <- as.character(metric %||% "ArticulationPoints")[1L]
  numeric_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
  if (!metric %in% numeric_cols) {
    stop("`metric` must be a numeric facet_summary column: ",
         paste(numeric_cols, collapse = ", "), call. = FALSE)
  }
  tbl <- tbl |>
    dplyr::arrange(dplyr::desc(.data[[metric]]), .data$Facet)
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    barplot_rot45(
      height = suppressWarnings(as.numeric(tbl[[metric]])),
      labels = tbl$Facet,
      col = pal["vulnerability"],
      main = main %||% paste("Facet network", metric),
      ylab = metric,
      label_angle = label_angle,
      mar_bottom = 7
    )
  }
  invisible(new_mfrm_plot_data(
    "network_analysis",
    list(
      plot = "facet_summary",
      table = tbl,
      metric = metric,
      title = main %||% paste("Facet network", metric),
      subtitle = "Facet-level aggregation of design-network vulnerability indicators",
      legend = new_plot_legend(metric, "facet_metric", "bar", pal["vulnerability"]),
      reference_lines = new_reference_lines(),
      preset = style$name
    )
  ))
}

draw_rater_network_bundle <- function(x,
                                      type = c("network", "centrality", "severity", "matrix"),
                                      metric = NULL,
                                      top_n = 20,
                                      draw = TRUE,
                                      main = NULL,
                                      palette = NULL,
                                      label_angle = 45,
                                      preset = c("standard", "publication", "compact", "monochrome")) {
  type <- match.arg(tolower(as.character(type[1])),
                    c("network", "centrality", "severity", "matrix"))
  style <- resolve_plot_preset(preset)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      agreement = style$accent_primary,
      disagreement = style$accent_secondary,
      lenient = style$accent_tertiary,
      severe = style$fail,
      neutral = style$neutral
    )
  )
  top_n <- max(1L, as.integer(top_n))
  nodes <- as.data.frame(x$node_metrics %||% data.frame(), stringsAsFactors = FALSE)
  edges <- as.data.frame(x$edge_metrics %||% data.frame(), stringsAsFactors = FALSE)
  mode <- as.character(x$settings$mode %||% x$summary$Mode %||% "network")[1L]
  directed <- isTRUE(x$settings$directed %||% identical(mode, "severity_direction"))

  if (identical(type, "centrality")) {
    if (nrow(nodes) == 0L) stop("No rater node metrics are available.", call. = FALSE)
    metric <- as.character(metric %||% if (directed) "Strength" else "Strength")[1L]
    numeric_cols <- names(nodes)[vapply(nodes, is.numeric, logical(1))]
    if (!metric %in% numeric_cols) {
      stop("`metric` must be a numeric node_metrics column: ",
           paste(numeric_cols, collapse = ", "), call. = FALSE)
    }
    tbl <- nodes |>
      dplyr::arrange(dplyr::desc(.data[[metric]]), .data$Rater) |>
      dplyr::slice_head(n = top_n)
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl[[metric]])),
        labels = truncate_axis_label(as.character(tbl$Rater), width = 24L),
        col = pal["agreement"],
        main = main %||% paste("Rater-network", metric),
        ylab = metric,
        label_angle = label_angle,
        mar_bottom = 8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "rater_network",
      list(
        plot = "centrality",
        table = tbl,
        metric = metric,
        title = main %||% paste("Rater-network", metric),
        subtitle = "Node-level network metric based on observed rater-pair relationships",
        legend = new_plot_legend(metric, "node_metric", "bar", pal["agreement"]),
        reference_lines = new_reference_lines(),
        preset = style$name
      )
    )))
  }

  if (identical(type, "severity")) {
    if (nrow(nodes) == 0L || !"SeverityIndex" %in% names(nodes)) {
      stop("No rater severity-network metrics are available.", call. = FALSE)
    }
    tbl <- nodes |>
      dplyr::filter(is.finite(.data$SeverityIndex)) |>
      dplyr::arrange(dplyr::desc(.data$SeverityIndex), .data$Rater) |>
      dplyr::slice_head(n = top_n)
    if (nrow(tbl) == 0L) {
      stop("No finite SeverityIndex values are available.", call. = FALSE)
    }
    cols <- ifelse(tbl$SeverityIndex > 0, pal["severe"],
                   ifelse(tbl$SeverityIndex < 0, pal["lenient"], pal["neutral"]))
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl$SeverityIndex)),
        labels = truncate_axis_label(as.character(tbl$Rater), width = 24L),
        col = cols,
        main = main %||% "Rater-network severity index",
        ylab = "-log(out-strength / in-strength)",
        label_angle = label_angle,
        mar_bottom = 8
      )
      graphics::abline(h = 0, lty = 2,
                       col = grDevices::adjustcolor(style$foreground, alpha.f = 0.65))
    }
    return(invisible(new_mfrm_plot_data(
      "rater_network",
      list(
        plot = "severity",
        table = tbl,
        metric = "SeverityIndex",
        title = main %||% "Rater-network severity index",
        subtitle = "Positive values indicate relatively severe; negative values indicate relatively lenient",
        legend = new_plot_legend(
          label = c("More severe", "More lenient", "Balanced"),
          role = c("status", "status", "status"),
          aesthetic = c("bar", "bar", "bar"),
          value = c(pal["severe"], pal["lenient"], pal["neutral"])
        ),
        reference_lines = new_reference_lines("h", 0, "Balanced", "dashed", "reference"),
        preset = style$name
      )
    )))
  }

  if (identical(type, "matrix")) {
    raters <- sort(unique(c(as.character(nodes$Rater), as.character(edges$From), as.character(edges$To))))
    mat <- matrix(NA_real_, nrow = length(raters), ncol = length(raters),
                  dimnames = list(raters, raters))
    if (nrow(edges) > 0L) {
      for (i in seq_len(nrow(edges))) {
        from <- as.character(edges$From[i])
        to <- as.character(edges$To[i])
        w <- suppressWarnings(as.numeric(edges$Weight[i]))
        if (!from %in% raters || !to %in% raters || !is.finite(w)) next
        mat[from, to] <- w
        if (!directed) mat[to, from] <- w
      }
    }
    if (isTRUE(draw) && length(raters) > 0L) {
      apply_plot_preset(style)
      graphics::image(
        x = seq_along(raters),
        y = seq_along(raters),
        z = t(mat[nrow(mat):1, , drop = FALSE]),
        axes = FALSE,
        col = grDevices::colorRampPalette(c(style$fill_soft, pal["agreement"], pal["disagreement"]))(64),
        xlab = "Rater",
        ylab = "Rater",
        main = main %||% paste("Rater-network", mode, "matrix")
      )
      graphics::axis(1, at = seq_along(raters),
                     labels = truncate_axis_label(raters, width = 14L),
                     las = 2, cex.axis = style$axis_cex)
      graphics::axis(2, at = seq_along(raters),
                     labels = rev(truncate_axis_label(raters, width = 14L)),
                     las = 2, cex.axis = style$axis_cex)
      graphics::box()
    }
    return(invisible(new_mfrm_plot_data(
      "rater_network",
      list(
        plot = "matrix",
        matrix = mat,
        edges = edges,
        title = main %||% paste("Rater-network", mode, "matrix"),
        subtitle = "Edge-weight matrix for custom heatmap or graph visualization",
        legend = new_plot_legend("Edge weight", "edge_weight", "fill", pal["agreement"]),
        reference_lines = new_reference_lines(),
        preset = style$name
      )
    )))
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("`plot(..., type = 'network')` requires the `igraph` package.",
         call. = FALSE)
  }
  vertices <- data.frame(
    name = sort(unique(c(as.character(nodes$Rater), as.character(edges$From), as.character(edges$To)))),
    stringsAsFactors = FALSE
  )
  graph_edges <- edges |>
    dplyr::select("From", "To", dplyr::everything())
  graph <- igraph::graph_from_data_frame(
    d = graph_edges,
    directed = directed,
    vertices = vertices
  )
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    weights <- suppressWarnings(as.numeric(igraph::E(graph)$Weight))
    edge_width <- if (length(weights) > 0L && any(is.finite(weights) & weights > 0)) {
      1 + 5 * weights / max(weights, na.rm = TRUE)
    } else {
      1
    }
    sev <- nodes$SeverityIndex[match(igraph::V(graph)$name, nodes$Rater)]
    vertex_col <- if (directed && length(sev) == igraph::vcount(graph)) {
      ifelse(is.finite(sev) & sev > 0, pal["severe"],
             ifelse(is.finite(sev) & sev < 0, pal["lenient"], pal["neutral"]))
    } else {
      pal["agreement"]
    }
    igraph::plot.igraph(
      graph,
      layout = igraph::layout_nicely(graph),
      vertex.label = truncate_axis_label(igraph::V(graph)$name, width = 16L),
      vertex.color = vertex_col,
      vertex.frame.color = style$foreground,
      vertex.label.color = style$foreground,
      vertex.size = 18,
      edge.width = edge_width,
      edge.color = grDevices::adjustcolor(style$axis, alpha.f = 0.55),
      edge.arrow.size = if (directed) 0.35 else 0,
      main = main %||% paste("Rater-network", mode)
    )
  }
  invisible(new_mfrm_plot_data(
    "rater_network",
    list(
      plot = "network",
      nodes = nodes,
      edges = edges,
      mode = mode,
      directed = directed,
      title = main %||% paste("Rater-network", mode),
      subtitle = "Graph of observed pairwise rater relationships",
      legend = new_plot_legend("Edge weight", "edge_weight", "edge_width", pal["agreement"]),
      reference_lines = new_reference_lines(),
      preset = style$name
    )
  ))
}

draw_halo_network_bundle <- function(x,
                                     type = c("edge_distribution", "halo_summary", "network", "matrix"),
                                     metric = NULL,
                                     top_n = 20,
                                     draw = TRUE,
                                     main = NULL,
                                     palette = NULL,
                                     label_angle = 45,
                                     preset = c("standard", "publication", "compact", "monochrome")) {
  type <- match.arg(tolower(as.character(type[1])),
                    c("edge_distribution", "halo_summary", "network", "matrix"))
  style <- resolve_plot_preset(preset)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      halo = style$accent_secondary,
      non_halo = style$accent_primary,
      node = style$accent_tertiary,
      warning = style$fail,
      review = style$warn,
      ok = style$accent_primary,
      neutral = style$neutral
    )
  )
  top_n <- max(1L, as.integer(top_n))
  nodes <- as.data.frame(x$node_metrics %||% data.frame(), stringsAsFactors = FALSE)
  edges <- as.data.frame(x$edge_metrics %||% data.frame(), stringsAsFactors = FALSE)
  pairs <- as.data.frame(x$pair_metrics %||% data.frame(), stringsAsFactors = FALSE)

  if (identical(type, "edge_distribution")) {
    if (nrow(pairs) == 0L || !"AbsEstimate" %in% names(pairs)) {
      stop("No halo pair metrics are available.", call. = FALSE)
    }
    tbl <- pairs |>
      dplyr::filter(.data$RetainedByN, is.finite(.data$AbsEstimate)) |>
      dplyr::mutate(EdgeType = factor(.data$EdgeType, levels = c("halo", "non_halo")))
    if (nrow(tbl) == 0L) {
      stop("No finite halo/non-halo weights are available.", call. = FALSE)
    }
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      graphics::boxplot(
        AbsEstimate ~ EdgeType,
        data = tbl,
        col = c(pal["halo"], pal["non_halo"]),
        border = style$foreground,
        ylab = "Absolute correlation",
        xlab = "",
        main = main %||% "Halo vs non-halo edge weights"
      )
      graphics::stripchart(
        AbsEstimate ~ EdgeType,
        data = tbl,
        vertical = TRUE,
        method = "jitter",
        pch = 16,
        col = grDevices::adjustcolor(style$foreground, alpha.f = 0.35),
        add = TRUE
      )
    }
    return(invisible(new_mfrm_plot_data(
      "halo_network",
      list(
        plot = "edge_distribution",
        table = tbl,
        metric = "AbsEstimate",
        title = main %||% "Halo vs non-halo edge weights",
        subtitle = "Same-rater cross-criterion edges are compared with other rater-by-criterion edges",
        legend = new_plot_legend(
          label = c("Halo edge", "Non-halo edge"),
          role = c("edge_type", "edge_type"),
          aesthetic = c("box", "box"),
          value = c(pal["halo"], pal["non_halo"])
        ),
        reference_lines = new_reference_lines(),
        preset = style$name
      )
    )))
  }

  if (identical(type, "halo_summary")) {
    tbl <- as.data.frame(x$halo_summary_by_rater %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L) stop("No per-rater halo summary is available.", call. = FALSE)
    metric <- as.character(metric %||% "MeanHaloWeight")[1L]
    numeric_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
    if (!metric %in% numeric_cols) {
      stop("`metric` must be a numeric halo_summary_by_rater column: ",
           paste(numeric_cols, collapse = ", "), call. = FALSE)
    }
    tbl <- tbl |>
      dplyr::arrange(dplyr::desc(.data[[metric]]), .data$Rater) |>
      dplyr::slice_head(n = top_n)
    status_cols <- if ("ReviewStatus" %in% names(tbl)) {
      ifelse(
        tbl$ReviewStatus == "warning",
        pal["warning"],
        ifelse(tbl$ReviewStatus == "review", pal["review"], pal["ok"])
      )
    } else {
      pal["halo"]
    }
    if (isTRUE(draw)) {
      apply_plot_preset(style)
      barplot_rot45(
        height = suppressWarnings(as.numeric(tbl[[metric]])),
        labels = truncate_axis_label(as.character(tbl$Rater), width = 24L),
        col = status_cols,
        main = main %||% paste("Rater halo", metric),
        ylab = metric,
        label_angle = label_angle,
        mar_bottom = 8
      )
    }
    return(invisible(new_mfrm_plot_data(
      "halo_network",
      list(
        plot = "halo_summary",
        table = tbl,
        metric = metric,
        title = main %||% paste("Rater halo", metric),
        subtitle = "Same-rater cross-criterion edge summary by rater",
        legend = new_plot_legend(
          label = c("Warning", "Review", "OK"),
          role = c("status", "status", "status"),
          aesthetic = c("bar", "bar", "bar"),
          value = c(pal["warning"], pal["review"], pal["ok"])
        ),
        reference_lines = new_reference_lines(),
        preset = style$name
      )
    )))
  }

  if (identical(type, "matrix")) {
    node_names <- sort(unique(c(as.character(nodes$Node), as.character(edges$From), as.character(edges$To))))
    mat <- matrix(NA_real_, nrow = length(node_names), ncol = length(node_names),
                  dimnames = list(node_names, node_names))
    if (nrow(edges) > 0L) {
      for (i in seq_len(nrow(edges))) {
        from <- as.character(edges$From[i])
        to <- as.character(edges$To[i])
        w <- suppressWarnings(as.numeric(edges$SignedWeight[i]))
        if (!from %in% node_names || !to %in% node_names || !is.finite(w)) next
        mat[from, to] <- w
        mat[to, from] <- w
      }
    }
    if (isTRUE(draw) && length(node_names) > 0L) {
      apply_plot_preset(style)
      graphics::image(
        x = seq_along(node_names),
        y = seq_along(node_names),
        z = t(mat[nrow(mat):1, , drop = FALSE]),
        axes = FALSE,
        col = grDevices::colorRampPalette(c(style$fill_soft, pal["non_halo"], pal["halo"]))(64),
        xlab = "Rater x criterion node",
        ylab = "Rater x criterion node",
        main = main %||% "Halo-network edge matrix"
      )
      graphics::axis(1, at = seq_along(node_names),
                     labels = truncate_axis_label(node_names, width = 14L),
                     las = 2, cex.axis = style$axis_cex)
      graphics::axis(2, at = seq_along(node_names),
                     labels = rev(truncate_axis_label(node_names, width = 14L)),
                     las = 2, cex.axis = style$axis_cex)
      graphics::box()
    }
    return(invisible(new_mfrm_plot_data(
      "halo_network",
      list(
        plot = "matrix",
        matrix = mat,
        edges = edges,
        title = main %||% "Halo-network edge matrix",
        subtitle = "Signed retained edge weights for custom heatmaps",
        legend = new_plot_legend("Signed edge weight", "edge_weight", "fill", pal["halo"]),
        reference_lines = new_reference_lines(),
        preset = style$name
      )
    )))
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("`plot(..., type = 'network')` requires the `igraph` package.",
         call. = FALSE)
  }
  vertices <- if (nrow(nodes) > 0L) {
    data.frame(
      name = as.character(nodes$Node),
      Rater = as.character(nodes$Rater),
      Criterion = as.character(nodes$Criterion),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(name = character(), stringsAsFactors = FALSE)
  }
  graph_edges <- edges |>
    dplyr::select("From", "To", dplyr::everything())
  graph <- igraph::graph_from_data_frame(
    d = graph_edges,
    directed = FALSE,
    vertices = vertices
  )
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    weights <- suppressWarnings(as.numeric(igraph::E(graph)$Weight))
    edge_width <- if (length(weights) > 0L && any(is.finite(weights) & weights > 0)) {
      1 + 5 * weights / max(weights, na.rm = TRUE)
    } else {
      1
    }
    edge_type <- as.character(igraph::E(graph)$EdgeType)
    edge_col <- ifelse(edge_type == "halo",
                       grDevices::adjustcolor(pal["halo"], alpha.f = 0.75),
                       grDevices::adjustcolor(pal["non_halo"], alpha.f = 0.45))
    igraph::plot.igraph(
      graph,
      layout = igraph::layout_nicely(graph),
      vertex.label = truncate_axis_label(igraph::V(graph)$name, width = 16L),
      vertex.color = pal["node"],
      vertex.frame.color = style$foreground,
      vertex.label.color = style$foreground,
      vertex.size = 15,
      edge.width = edge_width,
      edge.color = edge_col,
      main = main %||% "Rater-by-criterion halo network"
    )
  }
  invisible(new_mfrm_plot_data(
    "halo_network",
    list(
      plot = "network",
      nodes = nodes,
      edges = edges,
      title = main %||% "Rater-by-criterion halo network",
      subtitle = "Halo edges connect criteria scored by the same rater",
      legend = new_plot_legend(
        label = c("Halo edge", "Non-halo edge"),
        role = c("edge_type", "edge_type"),
        aesthetic = c("edge", "edge"),
        value = c(pal["halo"], pal["non_halo"])
      ),
      reference_lines = new_reference_lines(),
      preset = style$name
    )
  ))
}

draw_subset_connectivity_bundle <- function(x,
                                            type = c("subset_observations", "facet_levels", "coverage_matrix", "linking_matrix", "design_matrix", "network"),
                                            draw = TRUE,
                                            main = NULL,
                                            palette = NULL,
                                            label_angle = 45,
                                            preset = c("standard", "publication", "compact", "monochrome")) {
  requested_type <- match.arg(tolower(as.character(type[1])), c("subset_observations", "facet_levels", "coverage_matrix", "linking_matrix", "design_matrix", "network"))
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
  nodes_tbl <- as.data.frame(x$nodes %||% data.frame(), stringsAsFactors = FALSE)
  edges_tbl <- as.data.frame(x$edges %||% data.frame(), stringsAsFactors = FALSE)

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

  if (identical(type, "network")) {
    if (nrow(nodes_tbl) == 0 || nrow(edges_tbl) == 0 ||
        !all(c("Node", "Facet", "Level", "Subset") %in% names(nodes_tbl)) ||
        !all(c("From", "To", "Weight", "Subset") %in% names(edges_tbl))) {
      stop("Subset connectivity node/edge tables are not available. Rebuild with subset_connectivity_report() from a fitted object.", call. = FALSE)
    }
    node_order <- nodes_tbl |>
      dplyr::mutate(
        Subset = suppressWarnings(as.integer(.data$Subset)),
        Facet = as.character(.data$Facet),
        Level = as.character(.data$Level)
      ) |>
      dplyr::arrange(.data$Subset, .data$Facet, .data$Level)
    edge_order <- edges_tbl |>
      dplyr::mutate(
        Subset = suppressWarnings(as.integer(.data$Subset)),
        Weight = suppressWarnings(as.numeric(.data$Weight))
      ) |>
      dplyr::arrange(.data$Subset, dplyr::desc(.data$Weight), .data$From, .data$To)
    if (isTRUE(draw)) {
      if (!requireNamespace("igraph", quietly = TRUE)) {
        message("`plot(..., type = \"network\")` requires the `igraph` package for drawing. Returning plot data only.")
      } else {
        apply_plot_preset(style)
        vertices <- node_order
        row.names(vertices) <- vertices$Node
        g <- igraph::graph_from_data_frame(
          d = edge_order[, c("From", "To", "Weight"), drop = FALSE],
          directed = FALSE,
          vertices = vertices
        )
        facet_levels <- sort(unique(as.character(igraph::V(g)$Facet)))
        facet_cols <- stats::setNames(
          grDevices::hcl.colors(max(3L, length(facet_levels)), palette = "Dark 3")[seq_along(facet_levels)],
          facet_levels
        )
        igraph::V(g)$color <- facet_cols[as.character(igraph::V(g)$Facet)]
        igraph::V(g)$size <- ifelse(as.character(igraph::V(g)$Facet) == "Person", 4.5, 7)
        igraph::V(g)$label <- if (igraph::vcount(g) <= 80L) as.character(igraph::V(g)$Level) else NA_character_
        igraph::V(g)$label.cex <- 0.65
        igraph::V(g)$frame.color <- grDevices::adjustcolor(style$foreground, alpha.f = 0.35)
        igraph::E(g)$width <- pmax(0.4, log1p(as.numeric(igraph::E(g)$Weight)))
        igraph::E(g)$color <- grDevices::adjustcolor(style$grid, alpha.f = 0.55)
        layout <- igraph::layout_with_fr(g, weights = sqrt(pmax(1, as.numeric(igraph::E(g)$Weight))))
        graphics::plot(
          g,
          layout = layout,
          main = if (is.null(main)) "Connectivity network" else as.character(main[1]),
          vertex.label.color = style$foreground,
          vertex.label.family = "sans"
        )
        graphics::legend(
          "topleft",
          legend = names(facet_cols),
          pch = 21,
          pt.bg = unname(facet_cols),
          col = grDevices::adjustcolor(style$foreground, alpha.f = 0.35),
          bty = "n",
          cex = 0.8
        )
      }
    }
    return(invisible(new_mfrm_plot_data(
      "subset_connectivity",
      list(
        plot = "network",
        nodes = as.data.frame(node_order, stringsAsFactors = FALSE),
        edges = as.data.frame(edge_order, stringsAsFactors = FALSE),
        title = main %||% "Connectivity network",
        subtitle = "Facet-level co-observation network; edge width is observation count",
        legend = new_plot_legend(
          label = sort(unique(as.character(node_order$Facet))),
          role = "facet",
          aesthetic = "node_color",
          value = sort(unique(as.character(node_order$Facet)))
        ),
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
                                     type = c("overall_scree", "facet_scree",
                                              "overall_parallel_scree", "facet_parallel_scree",
                                              "overall_parallel_excess", "facet_parallel_excess",
                                              "overall_loadings", "facet_loadings"),
                                     facet = NULL,
                                     component = 1L,
                                     top_n = 20L,
                                     draw = TRUE) {
  type <- match.arg(
    tolower(as.character(type[1])),
    c("overall_scree", "facet_scree",
      "overall_parallel_scree", "facet_parallel_scree",
      "overall_parallel_excess", "facet_parallel_excess",
      "overall_loadings", "facet_loadings")
  )
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
  if (type == "overall_parallel_scree") {
    return(invisible(plot_residual_pca(
      x,
      mode = "overall",
      plot_type = "parallel_scree",
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  if (type == "facet_parallel_scree") {
    return(invisible(plot_residual_pca(
      x,
      mode = "facet",
      facet = facet,
      plot_type = "parallel_scree",
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  if (type == "overall_parallel_excess") {
    return(invisible(plot_residual_pca(
      x,
      mode = "overall",
      plot_type = "parallel_excess",
      component = component,
      top_n = top_n,
      draw = draw
    )))
  }
  if (type == "facet_parallel_excess") {
    return(invisible(plot_residual_pca(
      x,
      mode = "facet",
      facet = facet,
      plot_type = "parallel_excess",
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

draw_facets_contract_bundle <- function(x,
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

  column_review <- as.data.frame(x$column_review %||% data.frame(), stringsAsFactors = FALSE)
  column_summary <- as.data.frame(x$column_summary %||% data.frame(), stringsAsFactors = FALSE)
  metric_checks <- as.data.frame(x$metric_checks %||% data.frame(), stringsAsFactors = FALSE)
  metric_by_table <- as.data.frame(x$metric_by_table %||% data.frame(), stringsAsFactors = FALSE)

  if (type == "column_coverage") {
    if (nrow(column_review) == 0 || !all(c("table_id", "component", "coverage", "available", "full_match") %in% names(column_review))) {
      stop("Column-review table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    tbl <- column_review
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
      "facets_contract_review",
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
      "facets_contract_review",
      list(plot = "table_coverage", table = tbl)
    )))
  }

  if (type == "metric_status") {
    if (nrow(metric_checks) == 0 || !"Pass" %in% names(metric_checks)) {
      stop("Metric-review table is not available. Run the full workflow (fit_mfrm -> diagnose_mfrm) first.", call. = FALSE)
    }
    status <- ifelse(is.na(metric_checks$Pass), "Not evaluated", ifelse(metric_checks$Pass %in% TRUE, "Pass", "Fail"))
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
      "facets_contract_review",
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
    "facets_contract_review",
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
#' - `mfrm_category_curves` -> overview (default), ogive, CCC / category
#'   probability / conditional probability, cumulative, total-information, and
#'   category-specific-information plots
#' - `mfrm_rating_scale` -> category-counts/threshold plots
#' - `mfrm_measurable` -> measurable-data coverage/count plots
#' - `mfrm_unexpected_after_bias` -> post-bias unexpected-response plots
#' - `mfrm_output_bundle` -> graph/score output-file diagnostics,
#'   including `type = "score_se"` when scorefile SE columns are available
#' - `mfrm_residual_pca` -> residual PCA scree, parallel-analysis, or
#'   loadings views via [plot_residual_pca()]
#' - `mfrm_specifications` -> facet/anchor/convergence plots
#' - `mfrm_data_quality` -> dashboard, quality-flag, score-map, facet-pattern,
#'   and row/category/missing-row plots
#' - `mfrm_facets_fit_review` -> FACETS-style df-sensitivity plot
#' - `mfrm_fit_measures` -> fit-status counts, Infit/Outfit scatter, measure
#'   intervals, and FACETS-style df-sensitivity plots
#' - `mfrm_iteration_report` -> replayed-iteration trajectories
#' - `mfrm_subset_connectivity` -> subset-observation/connectivity plots
#' - `mfrm_facet_statistics` -> facet statistic profile plots
#' - `mfrm_export_bundle` / `mfrm_summary_appendix_export` -> export handoff
#'   plots (`formats`, `artifact_groups`, `selection_tables`,
#'   `selection_handoff`, `selection_handoff_bundles`,
#'   `selection_handoff_roles`,
#'   `selection_handoff_role_sections`, `selection_bundles`,
#'   `selection_roles`, `selection_sections`)
#'
#' If a class is outside these families, use dedicated plotting helpers
#' or custom base R graphics on component tables.
#'
#' For `mfrm_category_curves`, pass `preset = "monochrome"` for
#' grayscale/line-type output. Cumulative `.5` boundary lines are shown
#' only for interpretable in-range boundaries by default; use
#' `boundary_status = "all"` to show every finite boundary estimate or
#' `boundary_status = "none"` / `show_cumulative_boundaries = FALSE` to
#' suppress those vertical boundary lines. Use
#' `plot_data(x, component = "plot_long")` on a category-curve bundle when
#' you want one ggplot2/plotly-friendly table across all curve families.
#'
#' @section Interpreting output:
#' The returned object is plotting data (`mfrm_plot_data`) that captures
#' the selected route and reusable data; set `draw = TRUE` for immediate base graphics.
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
#'   fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
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
  if (inherits(x, "mfrm_facets_contract_review")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "column_coverage" else as.character(type[1])
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 40L
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    return(invisible(draw_facets_contract_bundle(
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
    ptype <- if (is.null(type)) "overview" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    cumulative_direction <- dots$cumulative_direction %||% "at_or_below"
    preset <- dots$preset %||% "standard"
    show_cumulative_boundaries <- dots$show_cumulative_boundaries %||% TRUE
    boundary_status <- dots$boundary_status %||% "in_range"
    return(invisible(draw_category_curves_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      cumulative_direction = cumulative_direction,
      preset = preset,
      show_cumulative_boundaries = show_cumulative_boundaries,
      boundary_status = boundary_status
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
    ptype <- if (is.null(type)) "dashboard" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    preset <- dots$preset %||% "standard"
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 30L
    return(invisible(draw_data_quality_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle,
      preset = preset,
      top_n = top_n
    )))
  }
  if (inherits(x, "mfrm_facets_fit_review")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "df_sensitivity" else as.character(type[1])
    if (!identical(tolower(ptype), "df_sensitivity")) {
      stop("`type` must be \"df_sensitivity\" for mfrm_facets_fit_review objects.", call. = FALSE)
    }
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    preset <- dots$preset %||% "standard"
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 30L
    return(invisible(draw_fit_measures_bundle(
      x,
      type = "df_sensitivity",
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle,
      preset = preset,
      top_n = top_n
    )))
  }
  if (inherits(x, "mfrm_fit_measures")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "status" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    preset <- dots$preset %||% "standard"
    ci_level <- dots$ci_level %||% NULL
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 30L
    return(invisible(draw_fit_measures_bundle(
      x,
      type = ptype,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle,
      preset = preset,
      ci_level = ci_level,
      top_n = top_n
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
  if (inherits(x, "mfrm_network_analysis")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "centrality" else as.character(type[1])
    metric <- dots$metric %||% NULL
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 20L
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    preset <- dots$preset %||% "standard"
    return(invisible(draw_network_analysis_bundle(
      x,
      type = ptype,
      metric = metric,
      top_n = top_n,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle,
      preset = preset
    )))
  }
  if (inherits(x, "mfrm_rater_network")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "network" else as.character(type[1])
    metric <- dots$metric %||% NULL
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 20L
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    preset <- dots$preset %||% "standard"
    return(invisible(draw_rater_network_bundle(
      x,
      type = ptype,
      metric = metric,
      top_n = top_n,
      draw = draw,
      main = main,
      palette = palette,
      label_angle = label_angle,
      preset = preset
    )))
  }
  if (inherits(x, "mfrm_halo_network")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "edge_distribution" else as.character(type[1])
    metric <- dots$metric %||% NULL
    top_n <- if ("top_n" %in% names(dots)) dots$top_n else 20L
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    preset <- dots$preset %||% "standard"
    return(invisible(draw_halo_network_bundle(
      x,
      type = ptype,
      metric = metric,
      top_n = top_n,
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
  if (inherits(x, "mfrm_export_bundle") || inherits(x, "mfrm_summary_appendix_export")) {
    draw <- if ("draw" %in% names(dots)) isTRUE(dots$draw) else TRUE
    ptype <- if (is.null(type)) "formats" else as.character(type[1])
    main <- dots$main %||% NULL
    palette <- dots$palette %||% NULL
    label_angle <- as.numeric(dots$label_angle %||% 45)
    selection_value <- as.character(dots$selection_value %||% "count")
    return(invisible(draw_export_handoff_bundle(
      x,
      type = ptype,
      selection_value = selection_value,
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
#' - diagnostic-basis guide for legacy versus strict fit paths
#' - global fit statistics
#' - approximate reliability/separation by facet
#' - top facet/person fit rows by absolute ZSTD
#' - counts of flagged diagnostics (unexpected, displacement, interactions)
#'
#' @section Interpreting output:
#' - `overview`: analysis scale, subset count, and residual-PCA mode.
#' - `diagnostic_basis`: plain-language map of which fit path was computed and
#'   what each path means statistically.
#' - `overall_fit`: global fit indices.
#' - `reliability`: facet separation/reliability block, including model and
#'   real bounds when available.
#' - `top_fit`: highest `|ZSTD|` elements for immediate inspection.
#' - `flags`: compact counts for key warning domains.
#'
#' @section Typical workflow:
#' 1. Run diagnostics with [diagnose_mfrm()], using `diagnostic_mode = "both"`
#'    for `RSM` / `PCM` when you want legacy continuity plus strict marginal screening.
#' 2. Review `summary(diag)` for major warnings and inspect `diagnostic_basis`
#'    before comparing legacy and strict outputs.
#' 3. Follow up with dedicated tables/plots for flagged domains.
#'
#' @return An object of class `summary.mfrm_diagnostics` with:
#' - `overview`: design-level counts and residual-PCA mode
#' - `status`: concise front-door status block for quick review
#' - `key_warnings`: highest-priority warnings to review first
#' - `next_actions`: recommended follow-up helpers
#' - `diagnostic_basis`: guide to legacy versus strict diagnostic targets
#' - `fit_standardization`: guide to the df convention used for fit ZSTD
#' - `overall_fit`: global fit block
#' - `precision_profile`: design-weighted precision summary across the
#'   information curve at decile theta points
#' - `precision_review`: separation / reliability / strata review for the
#'   sample- and population-basis modes (paired with `precision_profile`)
#' - `reliability`: facet-level separation/reliability summary
#' - `facets_chisq`: facets-style fixed-effect chi-square heterogeneity
#'   screen across non-person facets
#' - `interrater`: inter-rater agreement / pairwise correlation / rater
#'   separation overview when a Rater facet is present
#' - `misfit_flagged`: rows flagged by the Infit / Outfit / ZSTD
#'   misfit thresholds active for this fit
#' - `misfit_thresholds`: named numeric vector with the misfit
#'   `lower` / `upper` thresholds used to populate `misfit_flagged`
#' - `category_usage`: per-category response-frequency summary used
#'   to flag empty / collapsed categories
#' - `top_fit`: top `|ZSTD|` rows
#' - `marginal_fit`: optional strict marginal-fit overview when requested
#' - `top_marginal_cells`: largest strict marginal residual cells when requested
#' - `marginal_pairwise`: optional strict pairwise local-dependence overview
#' - `top_marginal_pairs`: largest strict pairwise residual summaries
#' - `marginal_guidance`: interpretation labels for strict marginal diagnostics
#' - `reporting_map`: manuscript-oriented guide to what is covered here versus
#'   which companion outputs should be consulted
#' - `flags`: compact flag counts for major diagnostics
#' - `notes`: short interpretation notes
#' - `digits`: numeric-print precision threaded through to
#'   `print.summary.mfrm_diagnostics()`
#' @seealso [diagnose_mfrm()], [summary.mfrm_fit()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' toy <- toy[toy$Person %in% unique(toy$Person)[1:4], ]
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' s <- summary(diag, top_n = 3)
#' s$key_warnings
#' # Look for: lines beginning with "MnSq misfit:" name the worst
#' #   element + Infit / Outfit values; "Unexpected responses flagged"
#' #   counts how many cell-level surprises the screen returned.
#' s$top_fit
#' # Look for: rows with |InfitZSTD| or |OutfitZSTD| > 2 are misfitting
#' #   at the 5% level; > 3 is misfitting at the 1% level. Investigate
#' #   in order of the AbsZ column.
#' s$facets_chisq
#' # Look for: FixedProb < 0.05 in each non-Person facet means the
#' #   facet contributes meaningful spread; FixedProb >= 0.05 means
#' #   that facet is statistically indistinguishable.
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
  precision_review_tbl <- tibble::as_tibble(precision_review(object, required = FALSE) %||% tibble::tibble())
  overall_fit <- tibble::as_tibble(object$overall_fit %||% tibble::tibble())
  subset_summary <- tibble::as_tibble(object$subsets$summary %||% tibble::tibble())
  approximation_tbl <- tibble::as_tibble(object$approximation_notes %||% tibble::tibble())
  diagnostic_basis_tbl <- tibble::as_tibble(object$diagnostic_basis %||% tibble::tibble())
  fit_standardization_tbl <- tibble::as_tibble(object$fit_standardization %||% tibble::tibble())
  marginal_fit_tbl <- tibble::as_tibble(object$marginal_fit$summary %||% tibble::tibble())
  marginal_step_summary <- tibble::as_tibble(object$marginal_fit$step_or_scale$summary_stats %||% tibble::tibble())
  marginal_facet_summary <- tibble::as_tibble(object$marginal_fit$facet_level$summary_stats %||% tibble::tibble())
  marginal_top_cells <- tibble::as_tibble(object$marginal_fit$top_cells %||% tibble::tibble())
  marginal_pairwise_tbl <- tibble::as_tibble(object$marginal_fit$pairwise$facet_summary %||% tibble::tibble())
  marginal_top_pairs_src <- tibble::as_tibble(object$marginal_fit$pairwise$top_pairs %||% tibble::tibble())
  marginal_guidance_tbl <- tibble::as_tibble(object$marginal_fit$guidance %||% tibble::tibble())
  interrater_summary_tbl <- tibble::as_tibble(object$interrater$summary %||% tibble::tibble())
  facets_chisq_tbl <- tibble::as_tibble(object$facets_chisq %||% tibble::tibble())
  marginal_available <- isTRUE(object$marginal_fit$available)
  marginal_pairwise_available <- isTRUE(object$marginal_fit$pairwise$available)
  diagnostic_mode <- as.character(object$diagnostic_mode %||% "legacy")

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
    DiagnosticMode = diagnostic_mode,
    Method = resolve_public_mfrm_method(
      summary_method = precision_profile_tbl$Method[1] %||% NA_character_
    ),
    PrecisionTier = as.character(precision_profile_tbl$PrecisionTier[1] %||% NA_character_),
    MarginalFit = if (marginal_available) "available" else "not_available"
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
    top_keep <- intersect(
      c(
        "Facet", "Level", "Infit", "Outfit", "InfitZSTD", "OutfitZSTD",
        "DF_Infit", "DF_Outfit",
        "InfitZSTD_FACETS", "OutfitZSTD_FACETS",
        "DF_Infit_FACETS", "DF_Outfit_FACETS"
      ),
      names(fit_tbl)
    )
    top_fit <- fit_tbl |>
      dplyr::mutate(
        AbsZ = pmax(abs(.data$InfitZSTD), abs(.data$OutfitZSTD), na.rm = TRUE)
      ) |>
      dplyr::arrange(dplyr::desc(.data$AbsZ)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select(dplyr::all_of(top_keep), "AbsZ")
  }

  # MnSq misfit threshold (Linacre, 0.5-1.5 acceptance band).
  # Tracked separately from `top_fit` so the auto-flag below can name the
  # offending element in `key_warnings` without depending on |ZSTD| ranking.
  misfit_thresholds <- mfrm_misfit_thresholds()
  misfit_lower <- as.numeric(misfit_thresholds["lower"])
  misfit_upper <- as.numeric(misfit_thresholds["upper"])
  misfit_flagged <- tibble::tibble()
  if (nrow(fit_tbl) > 0 && all(c("Facet", "Level", "Infit", "Outfit") %in% names(fit_tbl))) {
    misfit_flagged <- fit_tbl |>
      dplyr::filter(
        (is.finite(.data$Infit) & (.data$Infit < misfit_lower | .data$Infit > misfit_upper)) |
          (is.finite(.data$Outfit) & (.data$Outfit < misfit_lower | .data$Outfit > misfit_upper))
      )
  }

  # Category usage table: count per observed score, average person measure
  # within each category, and a disordering flag derived from any step
  # estimates that the diagnostics object exposes via the steps slot.
  category_usage <- tibble::tibble()
  if (nrow(obs_tbl) > 0 &&
      all(c("Person", "Observed") %in% names(obs_tbl))) {
    person_measure <- if (!is.null(object$measures) &&
                          all(c("Facet", "Level", "Estimate") %in%
                                names(object$measures))) {
      m <- as.data.frame(object$measures, stringsAsFactors = FALSE)
      m <- m[m$Facet == "Person", c("Level", "Estimate"), drop = FALSE]
      stats::setNames(suppressWarnings(as.numeric(m$Estimate)), m$Level)
    } else {
      numeric(0)
    }
    obs_local <- obs_tbl
    if (length(person_measure) > 0L) {
      obs_local$.PersonMeasure <- as.numeric(
        person_measure[as.character(obs_local$Person)]
      )
    } else {
      obs_local$.PersonMeasure <- NA_real_
    }
    category_usage <- obs_local |>
      dplyr::group_by(Category = .data$Observed) |>
      dplyr::summarize(
        Count = dplyr::n(),
        AvgMeasure = mean(.data$.PersonMeasure, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$Category)
    if (nrow(category_usage) > 1) {
      avg <- category_usage$AvgMeasure
      category_usage$Disordering <- ifelse(
        is.finite(avg) & c(FALSE, diff(avg) < 0),
        "AvgMeasure decreases vs previous category",
        ""
      )
    } else if (nrow(category_usage) == 1L) {
      category_usage$Disordering <- ""
    }
  }

  top_marginal_cells <- tibble::tibble()
  if (nrow(marginal_top_cells) > 0) {
    keep_marginal <- intersect(
      c(
        "CellType", "StepFacet", "Facet", "Level", "Category",
        "ObservedCount", "ExpectedCount", "PropDiff", "StdResidual", "AbsStdResidual"
      ),
      names(marginal_top_cells)
    )
    top_marginal_cells <- marginal_top_cells |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select(dplyr::all_of(keep_marginal))
  }

  top_marginal_pairs <- tibble::tibble()
  if (nrow(marginal_top_pairs_src) > 0) {
    keep_pairwise <- intersect(
      c(
        "Facet", "Level1", "Level2", "LevelPairCount", "OpportunityWeight",
        "ExactAgreement", "ExpectedExactAgreement", "ExactGap", "ExactStdResidual",
        "AdjacentAgreement", "ExpectedAdjacentAgreement", "AdjacentGap", "AdjacentStdResidual",
        "Flagged"
      ),
      names(marginal_top_pairs_src)
    )
    top_marginal_pairs <- marginal_top_pairs_src |>
      dplyr::slice_head(n = top_n) |>
      dplyr::select(dplyr::all_of(keep_pairwise))
  }

  unexpected_n <- suppressWarnings(as.integer(object$unexpected$summary$UnexpectedN[1] %||% NA_integer_))
  displacement_flagged <- suppressWarnings(as.integer(object$displacement$summary$FlaggedLevels[1] %||% NA_integer_))
  interaction_n <- if (!is.null(object$interactions)) nrow(object$interactions) else NA_integer_
  interrater_pairs <- suppressWarnings(as.integer(object$interrater$summary$Pairs[1] %||% NA_integer_))
  marginal_flagged_groups <- if (marginal_available) {
    sum(as.logical(marginal_step_summary$Flagged %||% logical(0)), na.rm = TRUE) +
      sum(as.logical(marginal_facet_summary$Flagged %||% logical(0)), na.rm = TRUE)
  } else {
    NA_integer_
  }
  marginal_flagged_pairs <- if (marginal_pairwise_available) {
    sum(as.logical(object$marginal_fit$pairwise$pair_stats$Flagged %||% logical(0)), na.rm = TRUE)
  } else {
    NA_integer_
  }

  flags <- tibble::tibble(
    Metric = c(
      "Unexpected responses",
      "Flagged displacement levels",
      "Interaction rows",
      "Inter-rater pairs",
      "Marginal fit flagged groups",
      "Marginal pairwise flagged level pairs"
    ),
    Count = c(unexpected_n, displacement_flagged, interaction_n, interrater_pairs, marginal_flagged_groups, marginal_flagged_pairs)
  )

  precision_tier <- as.character(precision_profile_tbl$PrecisionTier[1] %||% NA_character_)
  strict_path_status <- if (marginal_available) {
    "available"
  } else if (identical(diagnostic_mode, "legacy")) {
    "not_requested"
  } else {
    "requested_not_available"
  }
  primary_screen <- dplyr::case_when(
    identical(diagnostic_mode, "both") && marginal_available ~
      "Read strict marginal fit first; use legacy residuals for continuity and follow-up.",
    identical(diagnostic_mode, "marginal_fit") && marginal_available ~
      "Use strict marginal fit as the primary screen for first-order and pairwise follow-up.",
    identical(diagnostic_mode, "legacy") ~
      "Legacy residual diagnostics only; no strict marginal screen was requested.",
    TRUE ~
      "Requested strict marginal fit is not available for this run; fall back to the legacy path with caution."
  )

  key_warnings <- character(0)
  if (nrow(precision_review_tbl) > 0 && "Status" %in% names(precision_review_tbl)) {
    flagged_checks <- precision_review_tbl |>
      dplyr::filter(.data$Status %in% c("review", "warn"))
    if (nrow(flagged_checks) > 0) {
      key_warnings <- c(
        key_warnings,
        paste0("Precision review flagged ", nrow(flagged_checks), " review/warn checks.")
      )
    }
  }
  if (isTRUE(n_subsets > 1L)) {
    key_warnings <- c(key_warnings, "Multiple disconnected subsets were detected.")
  }
  if (isTRUE(!is.na(unexpected_n) && unexpected_n > 0L)) {
    key_warnings <- c(key_warnings, paste0("Unexpected responses flagged: ", unexpected_n, "."))
  }
  if (isTRUE(!is.na(displacement_flagged) && displacement_flagged > 0L)) {
    key_warnings <- c(key_warnings, paste0("Flagged displacement levels: ", displacement_flagged, "."))
  }
  # Name the worst MnSq offenders explicitly so the user does not have to
  # mentally apply the 0.5/1.5 acceptance band against the sorted top_fit
  # table.
  if (nrow(misfit_flagged) > 0) {
    worst <- misfit_flagged |>
      dplyr::mutate(
        WorstMnSq = pmax(
          ifelse(is.finite(.data$Outfit), abs(log(pmax(.data$Outfit, 1e-6))), 0),
          ifelse(is.finite(.data$Infit), abs(log(pmax(.data$Infit, 1e-6))), 0),
          na.rm = TRUE
        )
      ) |>
      dplyr::arrange(dplyr::desc(.data$WorstMnSq)) |>
      dplyr::slice_head(n = 3L)
    msgs <- vapply(seq_len(nrow(worst)), function(i) {
      sprintf(
        "MnSq misfit: %s:%s (Infit=%.2f, Outfit=%.2f; outside %.1f-%.1f).",
        worst$Facet[i], worst$Level[i],
        as.numeric(worst$Infit[i]), as.numeric(worst$Outfit[i]),
        misfit_lower, misfit_upper
      )
    }, character(1))
    key_warnings <- c(
      key_warnings,
      sprintf(
        "MnSq misfit flagged %d element(s) outside %.1f-%.1f (Linacre threshold).",
        nrow(misfit_flagged), misfit_lower, misfit_upper
      ),
      msgs
    )
  }
  if (isTRUE(marginal_available) && isTRUE(!is.na(marginal_flagged_groups) && marginal_flagged_groups > 0L)) {
    key_warnings <- c(
      key_warnings,
      paste0("Strict marginal fit flagged ", marginal_flagged_groups, " group-level summaries.")
    )
  }
  if (isTRUE(marginal_pairwise_available) && isTRUE(!is.na(marginal_flagged_pairs) && marginal_flagged_pairs > 0L)) {
    key_warnings <- c(
      key_warnings,
      paste0("Strict pairwise local dependence flagged ", marginal_flagged_pairs, " level pairs.")
    )
  }
  if (!marginal_available && !identical(diagnostic_mode, "legacy") && nrow(marginal_fit_tbl) > 0 && "Reason" %in% names(marginal_fit_tbl)) {
    key_warnings <- c(key_warnings, as.character(marginal_fit_tbl$Reason[1]))
  }
  key_warnings <- clean_summary_lines(key_warnings, max_n = 5L)
  if (length(key_warnings) == 0) {
    key_warnings <- "No immediate warnings from diagnostics summary."
  }

  next_actions <- c(
    "Inspect `diagnostic_basis` before comparing legacy residual evidence with strict marginal evidence."
  )
  if (!identical(as.character(precision_profile_tbl$Method[1] %||% NA_character_), "MML")) {
    next_actions <- c(
      next_actions,
      "Re-fit with `method = \"MML\"` if strict marginal diagnostics or formal SE/CI are required."
    )
  }
  if (marginal_available) {
    next_actions <- c(
      next_actions,
      "Review `top_marginal_cells` and `rating_scale_table(..., diagnostics = diag)` for first-order strict marginal follow-up."
    )
  }
  if (marginal_pairwise_available) {
    next_actions <- c(
      next_actions,
      "Review `top_marginal_pairs` for pairwise local-dependence follow-up."
    )
  }
  if (isTRUE(!is.na(unexpected_n) && unexpected_n > 0L) || isTRUE(!is.na(displacement_flagged) && displacement_flagged > 0L)) {
    next_actions <- c(
      next_actions,
      "Use `unexpected_response_table()` / `plot_unexpected()` and `displacement_table()` / `plot_displacement()` for case-level follow-up."
    )
  }
  if (!identical(as.character(object$residual_pca_mode %||% "none"), "none")) {
    next_actions <- c(
      next_actions,
      "Use `analyze_residual_pca()` if residual structure needs deeper follow-up."
    )
  }
  next_actions <- clean_summary_lines(next_actions, max_n = 4L)

  overall_status <- dplyr::case_when(
    any(key_warnings != "No immediate warnings from diagnostics summary.") ~ "follow_up_needed",
    identical(precision_tier, "exploratory") ~ "exploratory_screen_only",
    TRUE ~ "no_major_screening_flags"
  )
  status <- make_summary_block(
    "Overall status" = overall_status,
    "Diagnostic path" = diagnostic_mode,
    "Strict marginal fit" = strict_path_status,
    "Precision tier" = precision_tier,
    "Primary screen" = primary_screen
  )

  reporting_map <- tibble::tibble(
    Area = c(
      "Overall fit / reliability",
      "Precision basis / inferential caveats",
      "Strict marginal fit",
      "Residual PCA / local structure",
      "Unexpected responses / displacement",
      "Connectivity / subsets",
      "Manuscript checklist / export"
    ),
    CoveredHere = c(
      "yes",
      "yes",
      if (marginal_available) "yes" else if (identical(diagnostic_mode, "legacy")) "no" else "requested_not_available",
      "partial",
      "partial",
      "partial",
      "no"
    ),
    CompanionOutput = c(
      "summary(diagnostics)",
      "summary(diagnostics)",
      if (marginal_available) {
        "diagnostics$marginal_fit / rating_scale_table(..., diagnostics = diagnostics)"
      } else {
        "diagnose_mfrm(..., diagnostic_mode = \"both\")"
      },
      "analyze_residual_pca() / diagnostics$pca details",
      "unexpected_response_table() / displacement_table() / interaction tables",
      "subset_connectivity_report() / measurable_summary_table()",
      "reporting_checklist() / summary(build_apa_outputs(...))"
    )
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
  if (nrow(reliability_tbl) > 0) {
    notes <- c(
      notes,
      "Use `diagnostics$reliability` for facet-level separation/reliability. Use `diagnostics$interrater` only for observed agreement across matched rater contexts."
    )
  }
  if (nrow(precision_review_tbl) > 0 && "Status" %in% names(precision_review_tbl)) {
    flagged_checks <- precision_review_tbl |>
      dplyr::filter(.data$Status %in% c("review", "warn"))
    if (nrow(flagged_checks) > 0) {
      notes <- c(notes, paste0("Precision review flagged ", nrow(flagged_checks), " review/warn checks."))
    }
  }
  if (nrow(interrater_summary_tbl) == 0) {
    notes <- c(
      notes,
      "No inter-rater agreement block was produced for this run. Facet-level reliability summarizes separation/precision, not observed rater agreement."
    )
  }
  if (marginal_available) {
    notes <- c(
      notes,
      "Strict marginal fit was computed from latent-integrated first-order category counts."
    )
    if (marginal_pairwise_available) {
      notes <- c(
        notes,
        "Strict pairwise local-dependence checks were computed from posterior-integrated expected exact and adjacent agreement and should be read as exploratory screening summaries."
      )
    }
    notes <- c(
      notes,
      "Posterior predictive checking remains a planned corroborating follow-up for strict marginal flags and practical-significance review."
    )
    if (identical(diagnostic_mode, "both")) {
      notes <- c(
        notes,
        "Legacy residual diagnostics and strict marginal diagnostics target different quantities; do not compare their residual magnitudes directly."
      )
    }
  } else if (!identical(diagnostic_mode, "legacy") && nrow(marginal_fit_tbl) > 0 && "Reason" %in% names(marginal_fit_tbl)) {
    notes <- c(notes, as.character(marginal_fit_tbl$Reason[1]))
  }
  if (length(notes) == 0) {
    notes <- "No immediate warnings from diagnostics summary."
  }

  # Surface the fixed/random chi-square block (Linacre / Eckes
  # "are all elements equal?" headline) instead of leaving it in
  # diag$facets_chisq only.
  facets_chisq_overview <- tibble::tibble()
  if (nrow(facets_chisq_tbl) > 0) {
    keep_chi <- intersect(
      c("Facet", "Levels", "MeanMeasure", "SD",
        "FixedChiSq", "FixedDF", "FixedProb",
        "RandomChiSq", "RandomDF", "RandomProb"),
      names(facets_chisq_tbl)
    )
    facets_chisq_overview <- facets_chisq_tbl |>
      dplyr::select(dplyr::all_of(keep_chi))
  }

  # Lift inter-rater agreement summary into the printed output
  # rather than only counting flagged pairs.
  interrater_overview <- tibble::tibble()
  if (nrow(interrater_summary_tbl) > 0) {
    keep_ir <- intersect(
      c("RaterFacet", "Raters", "Pairs", "OpportunityCount",
        "ExactAgreement", "ExpectedExactAgreement",
        "AgreementMinusExpected", "AdjacentAgreement",
        "MeanAbsDiff", "MeanCorr",
        "RaterSeparation", "RaterReliability"),
      names(interrater_summary_tbl)
    )
    interrater_overview <- interrater_summary_tbl |>
      dplyr::select(dplyr::all_of(keep_ir))
  }

  out <- list(
    overview = overview,
    status = status,
    key_warnings = key_warnings,
    next_actions = next_actions,
    diagnostic_basis = diagnostic_basis_tbl,
    fit_standardization = fit_standardization_tbl,
    overall_fit = overall_fit,
    precision_profile = precision_profile_tbl,
    precision_review = precision_review_tbl,
    reliability = reliability_overview,
    facets_chisq = facets_chisq_overview,
    interrater = interrater_overview,
    misfit_flagged = misfit_flagged,
    misfit_thresholds = c(lower = misfit_lower, upper = misfit_upper),
    category_usage = category_usage,
    top_fit = top_fit,
    marginal_fit = marginal_fit_tbl,
    top_marginal_cells = top_marginal_cells,
    marginal_pairwise = marginal_pairwise_tbl,
    top_marginal_pairs = top_marginal_pairs,
    marginal_guidance = marginal_guidance_tbl,
    reporting_map = reporting_map,
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
  key_warning_lines <- if (summary_lines_are_default(
    x$key_warnings,
    "No immediate warnings from diagnostics summary."
  )) {
    "None."
  } else {
    x$key_warnings
  }

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
    if ("DiagnosticMode" %in% names(ov) && "MarginalFit" %in% names(ov)) {
      cat(sprintf("  Diagnostic mode: %s | Strict marginal fit: %s\n", ov$DiagnosticMode, ov$MarginalFit))
    }
  }
  if (!is.null(x$status) && nrow(x$status) > 0) {
    cat("\nStatus\n")
    for (i in seq_len(nrow(x$status))) {
      cat(" - ", x$status$Item[i], ": ", x$status$Value[i], "\n", sep = "")
    }
  }
  print_bullet_section("Key warnings", key_warning_lines)
  print_bullet_section("Next actions", x$next_actions)

  if (!is.null(x$overall_fit) && nrow(x$overall_fit) > 0) {
    cat("\nOverall fit\n")
    print(round_numeric_df(as.data.frame(x$overall_fit), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$fit_standardization) && nrow(x$fit_standardization) > 0) {
    cat("\nFit ZSTD standardization\n")
    print(round_numeric_df(as.data.frame(x$fit_standardization), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$diagnostic_basis) && nrow(x$diagnostic_basis) > 0) {
    cat("\nDiagnostic basis guide\n")
    print(as.data.frame(x$diagnostic_basis), row.names = FALSE)
  }
  if (!is.null(x$precision_profile) && nrow(x$precision_profile) > 0) {
    cat("\nPrecision basis\n")
    print(round_numeric_df(as.data.frame(x$precision_profile), digits = digits), row.names = FALSE)
  }
  precision_review_tbl <- precision_review(x, required = FALSE)
  if (!is.null(precision_review_tbl) && nrow(precision_review_tbl) > 0) {
    flagged_checks <- as.data.frame(precision_review_tbl)
    if ("Status" %in% names(flagged_checks)) {
      flagged_checks <- flagged_checks[flagged_checks$Status %in% c("review", "warn"), , drop = FALSE]
    }
    if (nrow(flagged_checks) > 0) {
      cat("\nPrecision review checks\n")
      print(flagged_checks, row.names = FALSE)
    }
  }
  if (!is.null(x$reliability) && nrow(x$reliability) > 0) {
    cat("\nFacet precision and spread\n")
    print(round_numeric_df(as.data.frame(x$reliability), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$facets_chisq) && nrow(x$facets_chisq) > 0) {
    cat("\nFacet variability (fixed-effect chi-square; null = all elements equal)\n")
    print(round_numeric_df(as.data.frame(x$facets_chisq), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$interrater) && nrow(x$interrater) > 0) {
    cat("\nInter-rater agreement summary\n")
    print(round_numeric_df(as.data.frame(x$interrater), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$top_fit) && nrow(x$top_fit) > 0) {
    cat("\nLargest |ZSTD| rows\n")
    print(round_numeric_df(as.data.frame(x$top_fit), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$misfit_flagged) && nrow(x$misfit_flagged) > 0) {
    thr <- x$misfit_thresholds %||% c(lower = 0.5, upper = 1.5)
    cat(sprintf(
      "\nMnSq misfit (outside %.1f-%.1f Linacre band; %d element(s))\n",
      as.numeric(thr["lower"]), as.numeric(thr["upper"]),
      nrow(x$misfit_flagged)
    ))
    keep_show <- intersect(
      c("Facet", "Level", "Infit", "InfitZSTD", "Outfit", "OutfitZSTD"),
      names(x$misfit_flagged)
    )
    print(round_numeric_df(
      as.data.frame(x$misfit_flagged)[, keep_show, drop = FALSE],
      digits = digits
    ), row.names = FALSE)
  }
  if (!is.null(x$category_usage) && nrow(x$category_usage) > 0) {
    cat("\nCategory usage (Count + AvgMeasure within each observed score)\n")
    print(round_numeric_df(as.data.frame(x$category_usage), digits = digits),
          row.names = FALSE)
  }
  if (!is.null(x$marginal_fit) && nrow(x$marginal_fit) > 0) {
    cat("\nStrict marginal fit\n")
    print(round_numeric_df(as.data.frame(x$marginal_fit), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$top_marginal_cells) && nrow(x$top_marginal_cells) > 0) {
    cat("\nLargest marginal residual cells\n")
    print(round_numeric_df(as.data.frame(x$top_marginal_cells), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$marginal_pairwise) && nrow(x$marginal_pairwise) > 0) {
    cat("\nStrict pairwise local dependence\n")
    print(round_numeric_df(as.data.frame(x$marginal_pairwise), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$top_marginal_pairs) && nrow(x$top_marginal_pairs) > 0) {
    cat("\nLargest marginal pairwise residuals\n")
    print(round_numeric_df(as.data.frame(x$top_marginal_pairs), digits = digits), row.names = FALSE)
  }
  if (!is.null(x$marginal_guidance) && nrow(x$marginal_guidance) > 0) {
    cat("\nStrict marginal guidance\n")
    print(as.data.frame(x$marginal_guidance), row.names = FALSE)
  }
  if (!is.null(x$reporting_map) && nrow(x$reporting_map) > 0) {
    cat("\nPaper reporting map\n")
    print(as.data.frame(x$reporting_map), row.names = FALSE)
  }
  if (!is.null(x$flags) && nrow(x$flags) > 0) {
    cat("\nFlag counts\n")
    print(as.data.frame(x$flags), row.names = FALSE)
  }

  if (length(x$notes) > 0 &&
      !summary_lines_are_default(x$notes, "No immediate warnings from diagnostics summary.")) {
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
#' - `top_rows`: strongest bias contrasts by `|t|`; bounded `GPCM`
#'   summaries also retain the profile-likelihood review columns when present.
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
#' toy <- toy[toy$Person %in% unique(toy$Person)[1:8], ]
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 1)
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
  lr_p_vals <- if ("LR Prob." %in% names(bias_tbl)) {
    suppressWarnings(as.numeric(bias_tbl$`LR Prob.`))
  } else {
    numeric(0)
  }
  lr_sig_n <- sum(is.finite(lr_p_vals) & lr_p_vals <= p_cut, na.rm = TRUE)

  # Multiple-testing corrected significant counts. We report Bonferroni and
  # Holm because they are the two corrections most commonly reported in
  # MFRM-bias screens; FDR is sometimes used in larger DIF-style
  # applications but is rarely reported in bias-by-rater contexts.
  finite_p <- p_vals[is.finite(p_vals)]
  bonferroni_n <- if (length(finite_p) > 0L) {
    sum(p.adjust(finite_p, method = "bonferroni") <= p_cut, na.rm = TRUE)
  } else 0L
  holm_n <- if (length(finite_p) > 0L) {
    sum(p.adjust(finite_p, method = "holm") <= p_cut, na.rm = TRUE)
  } else 0L

  overview <- tibble::tibble(
    FacetPair = interaction_label,
    InteractionOrder = interaction_order,
    InteractionMode = interaction_mode,
    Cells = nrow(bias_tbl),
    MeanAbsBias = mean(abs_bias, na.rm = TRUE),
    MaxAbsBias = max(abs_bias, na.rm = TRUE),
    Significant = sig_n,
    SignificantCut = p_cut,
    BonferroniSignificant = bonferroni_n,
    HolmSignificant = holm_n,
    ScreenPositive = sig_n,
    ScreeningCut = p_cut
  )
  if (any(is.finite(lr_p_vals))) {
    overview$LRScreenPositive <- lr_sig_n
  }

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
  base_keep <- c(level_cols, "Bias Size", "S.E.", "t", "Prob.", "Obs-Exp Average")
  likelihood_keep <- intersect(
    c("LR ChiSq", "LR Prob.", "Profile CI Lower", "Profile CI Upper",
      "Profile CI Status"),
    names(bias_tbl)
  )
  keep <- c(base_keep, likelihood_keep)
  if (all(base_keep %in% names(bias_tbl))) {
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
    if ("LRScreenPositive" %in% names(ov)) {
      cat(sprintf(
        "  GPCM profile-LR screen-positive (p <= %.3f): %s\n",
        as.numeric(ov$ScreeningCut), ov$LRScreenPositive
      ))
    }
    if (all(c("BonferroniSignificant", "HolmSignificant") %in% names(ov))) {
      cat(sprintf(
        "  Bonferroni significant: %s | Holm significant: %s (alpha = %.3f, m = %s)\n",
        ov$BonferroniSignificant, ov$HolmSignificant,
        as.numeric(ov$ScreeningCut), ov$Cells
      ))
    }
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
#' - estimation settings that affect identification/scoring interpretation
#' - facet-level estimate distribution (mean/SD/range)
#' - person measure distribution
#' - step/threshold checks
#' - a reporting map showing which companion summaries/tables should be used for
#'   manuscript-oriented data description, diagnostics, category checks, and draft
#'   reporting
#' - high/low person measures and extreme facet levels
#'
#' @section Interpreting output:
#' - `overview`: convergence and information criteria.
#' - `facet_overview`: per-facet spread and range of estimates.
#' - `person_overview`: distribution of person measures.
#' - `step_overview`: threshold spread and monotonicity checks.
#' - `settings_overview`: estimation settings that affect interpretation.
#' - `population_coding`: fitted categorical levels and contrasts that must be
#'   reused when scoring new persons under the population-model posterior.
#' - `key_warnings` / `notes`: short triage subset of retained zero-count score
#'   categories and latent-regression population-model caveats such as
#'   complete-case omissions, zero-variance design columns, missing
#'   coefficients, or unstable residual variance when present. Incomplete or
#'   non-finite covariates are normally handled before fitting as input errors
#'   or complete-case omissions; they appear here only if retained in a
#'   population-design check row.
#' - `caveats`: structured rows behind those warnings for appendix/export use;
#'   `print(summary(fit))` shows a compact `Caveats` block when rows are present.
#' - `reporting_map`: where to get companion outputs for manuscript reporting.
#' - `top_person` / `top_facet`: extreme estimates for quick triage.
#'
#' @section Typical workflow:
#' 1. Fit model with [fit_mfrm()].
#' 2. Run `summary(fit)` for first-pass diagnostics.
#' 3. For `RSM` / `PCM`, continue with [diagnose_mfrm()] for element-level fit
#'    checks. For bounded `GPCM`, continue with [compute_information()] /
#'    [plot_information()] or the fixed-calibration posterior scoring helpers.
#'
#' @return An object of class `summary.mfrm_fit` with:
#' - `overview`: global model/fit indicators
#' - `status`: concise front-door status block for quick review
#' - `key_warnings`: highest-priority warnings to review first
#' - `next_actions`: recommended follow-up helpers
#' - `population_overview`: current population-model basis, residual variance,
#'   and omission review
#' - `population_coefficients`: fitted latent-regression coefficients when a
#'   population model is active
#' - `population_design`: latent-regression design-matrix column check when a
#'   population model is active
#' - `population_coding`: categorical covariate levels and contrast provenance
#'   when a population model uses model-matrix coding
#' - `facet_overview`: per-facet estimate distribution summary
#' - `person_overview`: person-measure distribution summary
#' - `targeting`: person-versus-non-person facet targeting overview
#'   (Wright-map-style mean/SD comparison)
#' - `step_overview`: threshold/step diagnostics
#' - `slope_overview`: discrimination summary for `GPCM` fits
#' - `interaction_overview`: model-estimated facet-interaction summary
#'   when the fit was specified with `facet_interactions`
#' - `settings_overview`: estimation-settings overview that pins the
#'   configuration that affects identification/scoring
#' - `attached_diagnostics`: logical flag indicating whether the
#'   `mfrm_fit` was returned with diagnostics already attached
#' - `attached_diagnostics_cols`: character vector of diagnostic
#'   columns attached to `fit$facets$person` when
#'   `attached_diagnostics = TRUE`
#' - `row_retention`: row counts before and after preparation filters
#' - `preparation_notes`: structured preparation notes retained from
#'   `fit$prep`
#' - `reporting_map`: routing map showing which companion summaries
#'   and tables should be used for the four manuscript-oriented
#'   reporting sections (data description, diagnostics, category
#'   checks, draft reporting)
#' - `person_high` / `person_low`: highest and lowest person measures
#' - `facet_extremes`: extreme facet-level estimates
#' - `caveats`: structured warning/review rows for score-support and
#'   latent-regression population-model issues
#' - `notes`: short interpretation notes
#' - `digits`: numeric-print precision threaded through to
#'   `print.summary.mfrm_fit()`
#' @seealso [fit_mfrm()], [diagnose_mfrm()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(
#'   toy, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", quad_points = 5
#' )
#' s <- summary(fit)
#' s$overview[, c("Model", "Method", "Converged")]
#' # Look for: Converged = TRUE. If FALSE the fit is not safe to report;
#' #   raise `maxit`, relax `reltol`, or rerun with `quad_points = 31`.
#' s$person_overview
#' # Look for: Mean ~ 0 (logits) and SD ~ 1 are typical when the sample
#' #   is centred on the test difficulty. Min < -3 or Max > 3 with
#' #   `Extreme = "min"/"max"` rows indicates ceiling / floor cases.
#' s$targeting
#' # Look for: |Targeting| < ~0.5 logits across non-person facets is
#' #   comfortable. Larger absolute values mean the test is systematically
#' #   easier or harder than the person sample. SpreadRatio > 2 means
#' #   persons dominate facet variability; < 0.5 means facets dominate.
#' @export
summary.mfrm_fit <- function(object, digits = 3, top_n = 5, ...) {
  if (is.null(object$summary) || nrow(object$summary) == 0) {
    stop("`object` does not contain fit summary information.")
  }

  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  config <- object$config %||% list()
  overview <- tibble::as_tibble(object$summary)
  if ("Method" %in% names(overview)) {
    overview$Method <- public_mfrm_method_label(overview$Method)
  }
  if (!"MethodUsed" %in% names(overview)) {
    overview$MethodUsed <- as.character(config$method %||% NA_character_)
  } else {
    missing_used <- is.na(overview$MethodUsed) | !nzchar(trimws(as.character(overview$MethodUsed)))
    overview$MethodUsed[missing_used] <- as.character(config$method %||% NA_character_)
  }
  prep <- object$prep %||% list()
  population_raw <- object$population %||% list()
  person_raw <- object$facets$person
  if (is.null(person_raw)) person_raw <- tibble::tibble()
  other_raw <- object$facets$others
  if (is.null(other_raw)) other_raw <- tibble::tibble()
  step_raw <- object$steps
  if (is.null(step_raw)) step_raw <- tibble::tibble()
  slope_raw <- object$slopes
  if (is.null(slope_raw)) slope_raw <- tibble::tibble()
  interaction_raw <- object$interactions$effects
  if (is.null(interaction_raw)) interaction_raw <- tibble::tibble()

  person_tbl <- tibble::as_tibble(person_raw)
  other_tbl <- tibble::as_tibble(other_raw)
  step_tbl <- tibble::as_tibble(step_raw)
  slope_tbl <- tibble::as_tibble(slope_raw)
  interaction_tbl <- tibble::as_tibble(interaction_raw)
  population_coding <- population_coding_summary_table(population_raw)
  population_coding_variables <- if (nrow(population_coding) > 0L) {
    paste(population_coding$Variable, collapse = ", ")
  } else {
    ""
  }
  population_contrast_variables <- if (nrow(population_coding) > 0L) {
    paste(population_coding$Variable[nzchar(population_coding$Contrast)], collapse = ", ")
  } else {
    ""
  }
  population_overview <- tibble::tibble(
    PopulationModel = isTRUE(population_raw$active),
    PosteriorBasis = as.character(population_raw$posterior_basis %||% "legacy_mml"),
    Formula = if (!is.null(population_raw$formula)) paste(deparse(population_raw$formula), collapse = " ") else NA_character_,
    PersonRows = if (!is.null(population_raw$person_table)) nrow(as.data.frame(population_raw$person_table, stringsAsFactors = FALSE)) else NA_integer_,
    DesignColumns = if (!is.null(population_raw$design_matrix)) ncol(population_raw$design_matrix) else NA_integer_,
    CodingVariables = population_coding_variables,
    ContrastVariables = population_contrast_variables,
    Policy = as.character(population_raw$policy %||% NA_character_),
    ResidualVariance = as.numeric(population_raw$sigma2 %||% NA_real_),
    OmittedPersons = length(population_raw$omitted_persons %||% character(0)),
    OmittedRows = as.integer(population_raw$response_rows_omitted %||% NA_integer_)
  )
  population_coefficients <- tibble::tibble()
  if (isTRUE(population_raw$active) && length(population_raw$coefficients %||% numeric(0)) > 0) {
    coeff <- as.numeric(population_raw$coefficients)
    coeff_names <- names(population_raw$coefficients)
    if (is.null(coeff_names) || length(coeff_names) != length(coeff)) {
      coeff_names <- paste0("beta_", seq_along(coeff))
    }
    population_coefficients <- tibble::tibble(
      Term = as.character(coeff_names),
      Estimate = coeff
    )
  }
  population_design <- tibble::tibble()
  if (isTRUE(population_raw$active) && !is.null(population_raw$design_matrix)) {
    design_matrix <- as.matrix(population_raw$design_matrix)
    storage.mode(design_matrix) <- "double"
    if (ncol(design_matrix) > 0L) {
      design_columns <- colnames(design_matrix)
      if (is.null(design_columns) || length(design_columns) != ncol(design_matrix)) {
        design_columns <- paste0("X", seq_len(ncol(design_matrix)))
      }
      col_finite <- function(j) {
        vals <- as.numeric(design_matrix[, j])
        vals[is.finite(vals)]
      }
      col_stat <- function(j, fn) {
        vals <- col_finite(j)
        if (length(vals) == 0L) return(NA_real_)
        fn(vals)
      }
      col_sd <- function(j) {
        vals <- col_finite(j)
        if (length(vals) <= 1L) return(0)
        stats::sd(vals)
      }
      col_seq <- seq_len(ncol(design_matrix))
      non_missing <- vapply(col_seq, function(j) length(col_finite(j)), integer(1))
      sd_values <- vapply(col_seq, col_sd, numeric(1))
      population_design <- tibble::tibble(
        Column = as.character(design_columns),
        IsIntercept = as.character(design_columns) == "(Intercept)",
        PersonRows = nrow(design_matrix),
        NonMissing = non_missing,
        Complete = non_missing == nrow(design_matrix),
        Mean = vapply(col_seq, col_stat, numeric(1), fn = mean),
        SD = sd_values,
        Min = vapply(col_seq, col_stat, numeric(1), fn = min),
        Max = vapply(col_seq, col_stat, numeric(1), fn = max),
        ZeroVariance = non_missing <= 1L | sd_values <= sqrt(.Machine$double.eps)
      )
    }
  }

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

  # Targeting block.
  # Under the package's sum-to-zero identification on every non-person facet
  # the per-facet mean is constrained to 0, so a single "Person mean - Facet
  # mean" number collapses to the person mean. We still print the comparison
  # facet-by-facet because (i) the spread (SD/Span) tells the user whether
  # persons or facets dominate the test, and (ii) the row labels make the
  # implicit identification explicit.
  targeting_overview <- tibble::tibble()
  if (nrow(person_overview) > 0 && nrow(facet_overview) > 0) {
    person_mean <- as.numeric(person_overview$Mean[1])
    person_sd <- as.numeric(person_overview$SD[1])
    targeting_overview <- facet_overview |>
      dplyr::transmute(
        Facet = .data$Facet,
        PersonMean = person_mean,
        FacetMean = .data$MeanEstimate,
        Targeting = person_mean - .data$MeanEstimate,
        PersonSD = person_sd,
        FacetSD = .data$SDEstimate,
        SpreadRatio = person_sd / dplyr::na_if(.data$SDEstimate, 0)
      )
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

  slope_overview <- tibble::tibble()
  if (nrow(slope_tbl) > 0 && "Estimate" %in% names(slope_tbl)) {
    slope_overview <- tibble::tibble(
      Slopes = nrow(slope_tbl),
      Min = min(slope_tbl$Estimate, na.rm = TRUE),
      Max = max(slope_tbl$Estimate, na.rm = TRUE),
      GeometricMean = exp(mean(log(slope_tbl$Estimate), na.rm = TRUE)),
      Positive = all(is.finite(slope_tbl$Estimate) & slope_tbl$Estimate > 0)
    )
  }

  interaction_overview <- tibble::tibble()
  if (nrow(interaction_tbl) > 0 &&
      all(c("Interaction", "Estimate", "N", "Sparse") %in% names(interaction_tbl))) {
    interaction_overview <- interaction_tbl |>
      dplyr::group_by(.data$Interaction) |>
      dplyr::summarise(
        Cells = dplyr::n(),
        SparseCells = sum(.data$Sparse, na.rm = TRUE),
        MinN = suppressWarnings(min(.data$N, na.rm = TRUE)),
        MaxAbsEstimate = max(abs(.data$Estimate), na.rm = TRUE),
        .groups = "drop"
      )
  }

  settings_overview <- tibble::tibble(
    StepFacet = as.character(config$step_facet %||% NA_character_),
    SlopeFacet = as.character(config$slope_facet %||% NA_character_),
    NoncenterFacet = as.character(config$noncenter_facet %||% "Person"),
    WeightColumn = as.character(config$weight_col %||% NA_character_),
    QuadPoints = as.integer(config$estimation_control$quad_points %||% NA_integer_),
    RatingMin = as.numeric(config$rating_min %||% prep$rating_min %||% NA_real_),
    RatingMax = as.numeric(config$rating_max %||% prep$rating_max %||% NA_real_),
    RatingRangeSource = as.character(
      config$rating_range_source %||% prep$rating_range_source %||% "unknown"
    ),
    RatingMinSource = as.character(
      config$rating_min_source %||% prep$rating_min_source %||% "unknown"
    ),
    RatingMaxSource = as.character(
      config$rating_max_source %||% prep$rating_max_source %||% "unknown"
    ),
    DummyFacets = if (length(config$dummy_facets %||% character(0)) > 0L) {
      paste(sort(as.character(config$dummy_facets)), collapse = ", ")
    } else {
      ""
    },
    PositiveFacets = if (length(config$positive_facets %||% character(0)) > 0L) {
      paste(sort(as.character(config$positive_facets)), collapse = ", ")
    } else {
      ""
    },
    FacetInteractions = if (length(config$facet_interactions %||% character(0)) > 0L) {
      paste(as.character(config$facet_interactions), collapse = ", ")
    } else {
      ""
    },
    UnusedScoreCategories = "",
    UnusedScoreCategoryCount = 0L,
    UnusedScoreCategoryType = "none"
  )

  score_category_profile <- score_category_support_profile(prep = prep)
  score_category_caveats <- collect_mfrm_caveats(fit = object)
  population_caveats <- collect_mfrm_population_caveats(
    population_overview = population_overview,
    population_design = population_design,
    population_coefficients = population_coefficients
  )
  fit_caveats <- dplyr::bind_rows(score_category_caveats, population_caveats)
  unused_rows <- if (nrow(score_category_profile) > 0 && "ZeroCount" %in% names(score_category_profile)) {
    score_category_profile[as.logical(score_category_profile$ZeroCount), , drop = FALSE]
  } else {
    data.frame()
  }
  unused_score_categories <- sort(unique(as.integer(unused_rows$Category %||% integer(0))))
  unused_score_categories <- unused_score_categories[is.finite(unused_score_categories)]
  unused_score_category_type <- "none"
  if (nrow(unused_rows) > 0L && "UnusedCategoryType" %in% names(unused_rows)) {
    internal_unused <- unused_rows$UnusedCategoryType == "internal"
    unused_score_category_type <- dplyr::case_when(
      any(internal_unused) && any(!internal_unused) ~ "mixed",
      any(internal_unused) ~ "internal",
      TRUE ~ "boundary"
    )
    settings_overview$UnusedScoreCategories[1] <- paste(unused_score_categories, collapse = ", ")
    settings_overview$UnusedScoreCategoryCount[1] <- length(unused_score_categories)
    settings_overview$UnusedScoreCategoryType[1] <- unused_score_category_type
  }

  reporting_map <- tibble::tibble(
    Area = c(
      "Model identification / convergence",
      "Data structure / missingness",
      "Reliability / fit / residual PCA",
      "Category functioning",
      "Bias / DIF / interaction checks",
      "Draft reporting / checklist"
    ),
    CoveredHere = c("yes", "no", "no", "partial", "no", "no"),
    CompanionOutput = c(
      "summary(fit)",
      "summary(describe_mfrm_data(...))",
      "summary(diagnose_mfrm(fit))",
      "rating_scale_table() / category_structure_report() / category_curves_report()",
      "summary(estimate_bias(...)) / analyze_dff() / related bundle summaries",
      "reporting_checklist() / summary(build_apa_outputs(...))"
    )
  )

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
    status <- as.character(overview$ConvergenceStatus[1] %||% NA_character_)
    if (identical(status, "reviewable_warning")) {
      notes <- c(
        notes,
        "Optimizer returned a nonzero code, but the terminal gradient was already small; treat the fit as reviewable rather than an immediate hard failure."
      )
    } else {
      notes <- c(notes, "Optimization did not converge; interpret parameter estimates cautiously.")
    }
  }
  if (identical(as.character(overview$Method[1] %||% NA_character_), "MML")) {
    engine_requested <- as.character(overview$MMLEngineRequested[1] %||% NA_character_)
    engine_used <- as.character(overview$MMLEngineUsed[1] %||% NA_character_)
    if (!is.na(engine_requested) && !is.na(engine_used) && !identical(engine_requested, engine_used)) {
      notes <- c(
        notes,
        paste0("Requested mml_engine = '", engine_requested,
               "' was not available for this fit; the run used '", engine_used, "' instead.")
      )
    }
  }
  if (nrow(population_overview) > 0 && !isTRUE(population_overview$PopulationModel[1])) {
    fit_method <- as.character(overview$Method[1] %||% NA_character_)
    if (identical(fit_method, "MML")) {
      notes <- c(notes, "No population model was requested; current MML output uses the package's legacy unconditional prior.")
    }
  } else if (nrow(population_overview) > 0 && isTRUE(population_overview$PopulationModel[1])) {
    notes <- c(
      notes,
      "Population-model coefficients and residual variance describe the conditional normal prior used in the latent-regression MML branch."
    )
    if (nrow(population_coding) > 0L) {
      notes <- c(
        notes,
        "Categorical population-model coding is stored in `population_coding`; reuse those fitted levels and contrasts when scoring new persons."
      )
    }
  }
  if (nrow(step_overview) > 0 && !isTRUE(step_overview$Monotonic[1])) {
    notes <- c(notes, "Step estimates are not monotonic; verify category functioning.")
  }
  fit_caveat_messages <- if (nrow(fit_caveats) > 0 && "Message" %in% names(fit_caveats)) {
    as.character(fit_caveats$Message)
  } else {
    character(0)
  }
  if (length(fit_caveat_messages) > 0L) {
    notes <- c(notes, fit_caveat_messages)
  }
  if (nrow(slope_overview) > 0) {
    notes <- c(
      notes,
      "GPCM discriminations are reported under the package's positive log-slope identification with geometric-mean-one scaling."
    )
  }
  if (nrow(interaction_overview) > 0) {
    notes <- c(
      notes,
      "Facet interactions are model-estimated fixed effects with zero marginal sums; interpret them as deviations from the additive main-effects MFRM."
    )
  }
  if (length(notes) == 0) {
    notes <- "No immediate warnings from fit-level summary checks."
  }

  method_label <- as.character(overview$Method[1] %||% NA_character_)
  model_label <- as.character(overview$Model[1] %||% NA_character_)
  convergence_status <- as.character(overview$ConvergenceStatus[1] %||% NA_character_)
  convergence_severity <- as.character(overview$ConvergenceSeverity[1] %||% NA_character_)
  converged <- isTRUE(overview$Converged[1])
  engine_requested <- as.character(overview$MMLEngineRequested[1] %||% NA_character_)
  engine_used <- as.character(overview$MMLEngineUsed[1] %||% NA_character_)

  reporting_readiness <- dplyr::case_when(
    !converged && !identical(convergence_status, "reviewable_warning") ~ "review_before_reporting",
    identical(method_label, "MML") ~ "ready_for_diagnostics_and_reporting_follow_up",
    TRUE ~ "exploratory_fit_ready_for_diagnostics"
  )
  overall_status <- dplyr::case_when(
    converged ~ "usable_fit",
    identical(convergence_status, "reviewable_warning") ~ "reviewable_fit",
    TRUE ~ "fit_needs_review"
  )
  convergence_line <- if (isTRUE(is.finite(overview$TerminalGradientSupNorm[1] %||% NA_real_))) {
    paste0(
      convergence_status,
      " (severity: ", convergence_severity,
      ", sup-norm: ", round(as.numeric(overview$TerminalGradientSupNorm[1]), digits = digits), ")"
    )
  } else {
    paste0(convergence_status, " (severity: ", convergence_severity, ")")
  }
  engine_line <- if (identical(method_label, "MML") && !is.na(engine_used)) {
    if (!is.na(engine_requested) && !identical(engine_requested, engine_used)) {
      paste0(engine_used, " (requested ", engine_requested, ")")
    } else {
      engine_used
    }
  } else if (!is.na(method_label)) {
    method_label
  } else {
    NA_character_
  }
  status <- make_summary_block(
    "Overall status" = overall_status,
    "Convergence" = convergence_line,
    "Estimation path" = paste(model_label, "/", engine_line),
    "Reporting readiness" = reporting_readiness
  )

  row_retention <- as.data.frame(prep$row_retention %||% data.frame(), stringsAsFactors = FALSE)
  preparation_notes <- as.data.frame(prep$preparation_notes %||% data.frame(), stringsAsFactors = FALSE)
  preparation_review_messages <- if (nrow(preparation_notes) > 0L &&
                                     all(c("Severity", "Message") %in% names(preparation_notes))) {
    preparation_notes$Message[
      tolower(as.character(preparation_notes$Severity)) %in% c("review", "warning", "error")
    ]
  } else {
    character(0)
  }

  key_warnings <- clean_summary_lines(c(fit_caveat_messages, preparation_review_messages, notes), max_n = 4L)
  next_actions <- character(0)
  if (!identical(method_label, "MML")) {
    next_actions <- c(
      next_actions,
      "If formal SE/CI or strict marginal diagnostics are needed, re-fit with `method = \"MML\"`."
    )
  }
  if (nrow(population_overview) > 0 && isTRUE(population_overview$PopulationModel[1])) {
    next_actions <- c(
      next_actions,
      "Inspect `summary(fit)$population_coefficients` and `summary(fit)$population_coding` before reporting latent-regression effects.",
      "For new persons, use `predict_mfrm_units(..., person_data = ...)` or `sample_mfrm_plausible_values(..., person_data = ...)` so scoring reuses the fitted population-model coding.",
      "Report population coefficients as conditional-normal latent-regression parameters, not as a post hoc regression on EAP/MLE scores."
    )
  }
  if (identical(model_label, "GPCM")) {
    next_actions <- c(
      next_actions,
      "Run `diagnose_mfrm(fit, diagnostic_mode = \"both\")` for exploratory residual and strict-marginal screening.",
      "Use `compute_information()` / `plot_information()` for reporting-oriented precision follow-up."
    )
  } else {
    next_actions <- c(
      next_actions,
      "Run `diagnose_mfrm(fit, diagnostic_mode = \"both\")` for element-level fit review."
    )
  }
  if (nrow(interaction_overview) > 0) {
    next_actions <- c(
      next_actions,
      "Inspect `interaction_effect_table(fit)` and compare the additive and interaction fits on the same likelihood basis before reporting interaction claims."
    )
  }
  next_actions <- c(
    next_actions,
    "Use `plot(fit, type = \"wright\", preset = \"publication\")` for targeting and scale review.",
    "After diagnostics, use `reporting_checklist(fit, diagnostics = diagnostics)` for reporting readiness."
  )
  next_actions <- clean_summary_lines(next_actions, max_n = 6L)

  attached_diagnostics_flag <- isTRUE(config$attached_diagnostics)
  attached_diagnostics_cols <- as.character(config$attached_diagnostics_cols %||% character(0))

  out <- list(
    overview = overview,
    status = status,
    key_warnings = key_warnings,
    next_actions = next_actions,
    population_overview = population_overview,
    population_design = population_design,
    population_coefficients = population_coefficients,
    population_coding = population_coding,
    facet_overview = facet_overview,
    person_overview = person_overview,
    targeting = targeting_overview,
    step_overview = step_overview,
    slope_overview = slope_overview,
    interaction_overview = interaction_overview,
    settings_overview = settings_overview,
    attached_diagnostics = attached_diagnostics_flag,
    attached_diagnostics_cols = attached_diagnostics_cols,
    row_retention = row_retention,
    preparation_notes = preparation_notes,
    reporting_map = reporting_map,
    caveats = fit_caveats,
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

clean_summary_lines <- function(lines, max_n = Inf) {
  if (length(lines) == 0) return(character(0))
  lines <- trimws(as.character(lines))
  lines <- lines[!is.na(lines) & nzchar(lines)]
  lines <- unique(lines)
  if (is.finite(max_n) && length(lines) > max_n) {
    lines <- lines[seq_len(max_n)]
  }
  lines
}

make_summary_block <- function(...) {
  vals <- list(...)
  keys <- names(vals)
  if (is.null(keys)) {
    keys <- rep.int("", length(vals))
  }
  keep <- vapply(vals, function(x) {
    if (length(x) == 0) return(FALSE)
    if (all(is.na(x))) return(FALSE)
    if (is.character(x)) return(any(nzchar(trimws(x))))
    TRUE
  }, logical(1))
  if (!any(keep)) {
    return(tibble::tibble(Item = character(0), Value = character(0)))
  }
  tibble::tibble(
    Item = keys[keep],
    Value = vapply(vals[keep], function(x) paste(x, collapse = " "), character(1))
  )
}

summary_lines_are_default <- function(lines, default_line) {
  lines <- clean_summary_lines(lines)
  length(lines) == 1L && identical(lines, default_line)
}

print_bullet_section <- function(title, lines, prefix = " - ") {
  lines <- clean_summary_lines(lines)
  if (length(lines) == 0) return(invisible(NULL))
  cat("\n", title, "\n", sep = "")
  for (line in lines) cat(prefix, line, "\n", sep = "")
  invisible(NULL)
}

summary_caveat_display_table <- function(caveats) {
  caveats <- as.data.frame(caveats %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(caveats) == 0L) return(caveats)
  keep <- intersect(
    c("Area", "Severity", "Condition", "Categories", "Message", "RecommendedAction"),
    names(caveats)
  )
  if (length(keep) == 0L) return(caveats)
  caveats[, keep, drop = FALSE]
}

print_caveat_section <- function(caveats, title = "Caveats") {
  caveats <- summary_caveat_display_table(caveats)
  if (nrow(caveats) == 0L) return(invisible(NULL))
  cat("\n", title, "\n", sep = "")
  print(caveats, row.names = FALSE)
  invisible(NULL)
}

summary_preparation_display_table <- function(notes) {
  notes <- as.data.frame(notes %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(notes) == 0L) return(notes)
  keep <- intersect(
    c("Stage", "Severity", "Condition", "Count", "Affected", "RecommendedAction"),
    names(notes)
  )
  if (length(keep) == 0L) return(notes)
  notes[, keep, drop = FALSE]
}

print_preparation_section <- function(notes, title = "Data preparation notes") {
  notes <- summary_preparation_display_table(notes)
  if (nrow(notes) == 0L) return(invisible(NULL))
  cat("\n", title, "\n", sep = "")
  print(notes, row.names = FALSE)
  invisible(NULL)
}

#' @export
print.summary.mfrm_fit <- function(x, ...) {
  digits <- x$digits
  if (is.null(digits) || !is.finite(digits)) digits <- 3L
  overview <- round_numeric_df(as.data.frame(x$overview), digits = digits)
  key_warning_lines <- if (summary_lines_are_default(
    x$key_warnings,
    "No immediate warnings from fit-level summary checks."
  )) {
    "None."
  } else {
    x$key_warnings
  }

  cat("Many-Facet Rasch Model Summary\n")
  if (nrow(overview) > 0) {
    ov <- overview[1, , drop = FALSE]
    cat(sprintf(
      "  Model: %s | Method: %s | N: %s | Persons: %s | Facets: %s | Categories: %s\n",
      ov$Model, ov$Method, ov$N, ov$Persons, ov$Facets, ov$Categories
    ))
    if (isTRUE(x$attached_diagnostics)) {
      attached_cols <- as.character(x$attached_diagnostics_cols %||% character(0))
      if (length(attached_cols) > 0L) {
        cat(sprintf(
          "  Attached diagnostics: %s\n",
          paste(attached_cols, collapse = ", ")
        ))
      } else {
        cat("  Attached diagnostics: yes\n")
      }
    }
    used_public <- public_mfrm_method_label(ov$MethodUsed %||% NA_character_)
    if (!is.na(ov$MethodUsed %||% NA_character_) &&
        nzchar(as.character(ov$MethodUsed %||% "")) &&
        !identical(as.character(used_public), as.character(ov$Method))) {
      cat(sprintf("  Internal method label: %s\n", ov$MethodUsed))
    }
    if (identical(as.character(ov$Method %||% NA_character_), "MML") &&
        !is.na(ov$MMLEngineUsed %||% NA_character_)) {
      cat(sprintf(
        "  MML engine: %s (requested: %s)\n",
        ov$MMLEngineUsed %||% NA_character_,
        ov$MMLEngineRequested %||% NA_character_
      ))
      if (is.finite(ov$EMIterations %||% NA_real_)) {
        cat(sprintf(
          "  EM iterations: %s | EM converged: %s | Last relative change: %s\n",
          ov$EMIterations %||% NA,
          ifelse(isTRUE(ov$EMConverged), "Yes", "No"),
          ov$EMRelativeChange %||% NA_real_
        ))
      }
    }
  }
  if (nrow(x$status %||% data.frame()) > 0) {
    cat("\nStatus\n")
    for (i in seq_len(nrow(x$status))) {
      cat(" - ", x$status$Item[i], ": ", x$status$Value[i], "\n", sep = "")
    }
  }
  print_bullet_section("Key warnings", key_warning_lines)
  print_bullet_section("Next actions", x$next_actions)
  print_caveat_section(x$caveats)
  if (nrow(x$row_retention %||% data.frame()) > 0L &&
      "DroppedRows" %in% names(x$row_retention) &&
      any(suppressWarnings(as.numeric(x$row_retention$DroppedRows)) > 0, na.rm = TRUE)) {
    cat("\nRow retention\n")
    print(round_numeric_df(as.data.frame(x$row_retention), digits = digits), row.names = FALSE)
  }
  print_preparation_section(x$preparation_notes)

  if (nrow(overview) > 0) {
    ov <- overview[1, , drop = FALSE]
    cat("\nFit overview\n")
    cat(sprintf("  LogLik: %s | AIC: %s | BIC: %s\n", ov$LogLik, ov$AIC, ov$BIC))
    cat(sprintf(
      "  Converged: %s | Status: %s | Basis: %s | Fn evals: %s | Gr evals: %s\n",
      ifelse(isTRUE(ov$Converged), "Yes", "No"),
      ov$ConvergenceStatus %||% NA_character_,
      ov$ConvergenceBasis %||% NA_character_,
      ov$FunctionEvaluations %||% ov$Iterations %||% NA,
      ov$GradientEvaluations %||% NA
    ))
    if (is.finite(ov$TerminalGradientSupNorm %||% NA_real_)) {
      cat(sprintf(
        "  Terminal gradient: sup-norm = %s | RMS = %s | Review tol = %s\n",
        ov$TerminalGradientSupNorm,
        ov$TerminalGradientRMS %||% NA_real_,
        ov$GradientReviewTolerance %||% NA_real_
      ))
    }
    if (!is.na(ov$ConvergenceDetail %||% NA_character_) && nzchar(ov$ConvergenceDetail %||% "")) {
      cat(sprintf("  Optimization note: %s\n", ov$ConvergenceDetail))
    }
  }

  if (nrow(x$population_overview) > 0) {
    cat("\nPopulation basis\n")
    print(round_numeric_df(as.data.frame(x$population_overview), digits = digits), row.names = FALSE)
  }
  if (nrow(x$population_design %||% data.frame()) > 0) {
    cat("\nPopulation design matrix\n")
    print(round_numeric_df(as.data.frame(x$population_design), digits = digits), row.names = FALSE)
  }
  if (nrow(x$population_coding %||% data.frame()) > 0) {
    cat("\nPopulation covariate coding\n")
    print(as.data.frame(x$population_coding), row.names = FALSE)
  }
  if (nrow(x$population_coefficients) > 0) {
    cat("\nPopulation coefficients\n")
    print(round_numeric_df(as.data.frame(x$population_coefficients), digits = digits), row.names = FALSE)
  }

  if (nrow(x$facet_overview) > 0) {
    cat("\nFacet overview\n")
    print(round_numeric_df(as.data.frame(x$facet_overview), digits = digits), row.names = FALSE)
  }

  if (nrow(x$person_overview) > 0) {
    cat("\nPerson measure distribution\n")
    print(round_numeric_df(as.data.frame(x$person_overview), digits = digits), row.names = FALSE)
  }

  if (nrow(x$targeting %||% data.frame()) > 0) {
    cat("\nTargeting (Person vs facet means; sum-to-zero ID makes Targeting = Person mean)\n")
    print(round_numeric_df(as.data.frame(x$targeting), digits = digits), row.names = FALSE)
  }

  if (nrow(x$step_overview) > 0) {
    cat("\nStep parameter summary\n")
    print(round_numeric_df(as.data.frame(x$step_overview), digits = digits), row.names = FALSE)
  }
  if (nrow(x$slope_overview %||% data.frame()) > 0) {
    cat("\nSlope summary\n")
    print(round_numeric_df(as.data.frame(x$slope_overview), digits = digits), row.names = FALSE)
  }
  if (nrow(x$interaction_overview %||% data.frame()) > 0) {
    cat("\nFacet interaction summary\n")
    print(round_numeric_df(as.data.frame(x$interaction_overview), digits = digits), row.names = FALSE)
  }
  if (nrow(x$settings_overview %||% data.frame()) > 0) {
    cat("\nEstimation settings\n")
    print(round_numeric_df(as.data.frame(x$settings_overview), digits = digits), row.names = FALSE)
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

  if (nrow(x$reporting_map %||% data.frame()) > 0) {
    cat("\nPaper reporting map\n")
    print(as.data.frame(x$reporting_map), row.names = FALSE)
  }

  if (length(x$notes) > 0 &&
      !summary_lines_are_default(x$notes, "No immediate warnings from fit-level summary checks.")) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }

  invisible(x)
}
