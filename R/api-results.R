# High-level results entry point.

mfrm_results_status_row <- function(section, status, detail = "") {
  data.frame(
    Section = as.character(section),
    Status = as.character(status),
    Detail = as.character(detail %||% ""),
    stringsAsFactors = FALSE
  )
}

mfrm_results_component_class <- function(x) {
  if (is.null(x)) return("")
  paste(class(x), collapse = ", ")
}

mfrm_results_safe <- function(expr) {
  tryCatch(
    list(ok = TRUE, value = expr, message = ""),
    error = function(e) {
      list(ok = FALSE, value = NULL, message = conditionMessage(e))
    }
  )
}

mfrm_results_include_presets <- function() {
  standard <- c("fit", "diagnostics", "tables", "precision", "reporting", "categories", "plots")
  list(
    standard = standard,
    publication = unique(c(standard, "apa")),
    validation = unique(c(standard, "facets_fit")),
    facets = unique(c("fit", "diagnostics", "tables", "categories", "plots", "facets_fit")),
    bias = unique(c(standard, "bias")),
    bias_review = unique(c(standard, "bias")),
    misfit = unique(c(standard, "misfit")),
    misfit_review = unique(c(standard, "misfit")),
    linking = unique(c(standard, "linking")),
    anchors = unique(c(standard, "linking")),
    network = unique(c(standard, "network")),
    gpcm_review = standard,
    all = unique(c(standard, "facets_fit", "bias", "misfit", "linking", "network", "apa", "response_time"))
  )
}

mfrm_results_include_preset_table <- function() {
  presets <- mfrm_results_include_presets()
  data.frame(
    Preset = names(presets),
    Sections = vapply(presets, paste, collapse = ", ", character(1)),
    stringsAsFactors = FALSE
  )
}

mfrm_results_resolve_include <- function(include) {
  if (is.null(include) || length(include) == 0L) include <- "standard"
  requested <- unique(tolower(trimws(as.character(include))))
  requested <- requested[nzchar(requested)]
  include <- unique(tolower(trimws(as.character(include))))
  include <- include[nzchar(include)]
  if (length(include) == 0L) include <- "standard"

  aliases <- c(
    core = "standard",
    report = "reporting",
    reports = "reporting",
    publish = "publication",
    manuscript = "publication",
    validate = "validation",
    validation_review = "validation",
    facet = "facets",
    facet_review = "facets",
    fairness = "bias",
    fairness_screen = "bias",
    bias_screen = "bias",
    anchor = "linking",
    anchors = "linking",
    equating = "linking",
    drift = "linking",
    link = "linking",
    unexpected = "misfit",
    pathway = "misfit",
    misfit_screen = "misfit",
    rt = "response_time",
    timing = "response_time",
    response_times = "response_time",
    response_time_review = "response_time",
    gpcm = "gpcm_review",
    gpcm_validation = "gpcm_review",
    review = "precision",
    reviews = "precision",
    category = "categories",
    table = "tables",
    diagnostics = "diagnostics",
    diagnostic = "diagnostics",
    separation = "precision",
    fit = "fit"
  )
  include <- vapply(
    include,
    function(x) if (x %in% names(aliases)) aliases[[x]] else x,
    character(1)
  )

  presets <- mfrm_results_include_presets()
  preset_hits <- intersect(include, names(presets))
  include <- unique(unlist(
    lapply(include, function(x) {
      if (x %in% names(presets)) presets[[x]] else x
    }),
    use.names = FALSE
  ))

  allowed <- unique(unlist(presets, use.names = FALSE))
  bad <- setdiff(include, allowed)
  if (length(bad) > 0L) {
    stop(
      "Unsupported `include` value(s): ", paste(bad, collapse = ", "),
      ". Use a preset such as 'standard', 'publication', 'validation', ",
      "'facets', 'bias', 'misfit_review', 'linking', 'network', 'gpcm_review', ",
      "or 'all', or any section from: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  attr(include, "requested") <- if (length(requested) > 0L) requested else "standard"
  attr(include, "presets") <- preset_hits
  include
}

mfrm_results_response_time_aliases <- function() {
  c(
    "responsetime", "responsetimeseconds", "responsetimemilliseconds",
    "responsetimems", "responsemilliseconds", "responsems",
    "rt", "rtseconds", "rtmilliseconds", "rtms",
    "duration", "durationseconds", "durationmilliseconds",
    "elapsedtime", "elapsedseconds", "elapsedmilliseconds",
    "latency", "latencyseconds", "latencymilliseconds"
  )
}

mfrm_results_normalized_column_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

mfrm_results_response_time_candidates <- function(cols) {
  cols <- as.character(cols %||% character(0))
  norm <- mfrm_results_normalized_column_name(cols)
  cols[norm %in% mfrm_results_response_time_aliases()]
}

mfrm_results_maybe_response_time_column <- function(data) {
  if (!is.data.frame(data)) return(NULL)
  hits <- mfrm_results_response_time_candidates(names(data))
  if (length(hits) == 1L) hits[[1L]] else NULL
}

mfrm_results_resolve_response_time_column <- function(data, response_time = NULL) {
  if (!is.data.frame(data)) {
    stop(
      "Response-time review requires `response_time_data` or data-frame input ",
      "with the timing column.",
      call. = FALSE
    )
  }
  if (!is.null(response_time) && length(response_time) > 0L &&
      !is.na(response_time[1]) && nzchar(response_time[1])) {
    response_time <- as.character(response_time[1])
    if (!response_time %in% names(data)) {
      stop("Response-time column not found: ", response_time, call. = FALSE)
    }
    return(response_time)
  }
  hits <- mfrm_results_response_time_candidates(names(data))
  if (length(hits) == 1L) return(hits[[1L]])
  if (length(hits) > 1L) {
    stop(
      "Ambiguous response-time columns: ", paste(hits, collapse = ", "),
      ". Supply `response_time` explicitly.",
      call. = FALSE
    )
  }
  stop(
    "No response-time column was supplied or safely detected. ",
    "Use `response_time = \"ResponseTime\"` and, for fitted objects, ",
    "`response_time_data = original_data`.",
    call. = FALSE
  )
}

mfrm_results_find_standard_column <- function(cols, aliases, role) {
  lower <- tolower(cols)
  hits <- which(lower %in% aliases)
  if (length(hits) == 1L) return(cols[hits])
  if (length(hits) == 0L) {
    stop(
      "Could not safely infer the ", role, " column from a data.frame. ",
      "Use fit_mfrm() with explicit `person`, `facets`, and `score`, or use ",
      "mfrm_results_interactive() in an interactive session.",
      call. = FALSE
    )
  }
  stop(
    "Ambiguous ", role, " columns in data.frame: ",
    paste(cols[hits], collapse = ", "),
    ". Use fit_mfrm() explicitly or mfrm_results_interactive().",
    call. = FALSE
  )
}

mfrm_results_infer_standard_mapping <- function(data) {
  dat <- normalize_facets_mode_data(data)
  cols <- names(dat)
  person <- mfrm_results_find_standard_column(
    cols,
    aliases = c("person", "persons", "participant", "participants",
                "student", "students", "subject", "subjects",
                "examinee", "examinees"),
    role = "person"
  )
  score <- mfrm_results_find_standard_column(
    cols,
    aliases = c("score", "scores", "rating", "ratings", "mark", "marks"),
    role = "score"
  )

  weight_hits <- which(tolower(cols) %in% c("weight", "weights"))
  if (length(weight_hits) > 1L) {
    stop(
      "Ambiguous weight columns in data.frame: ",
      paste(cols[weight_hits], collapse = ", "),
      ". Use fit_mfrm() explicitly.",
      call. = FALSE
    )
  }
  weight <- if (length(weight_hits) == 1L) cols[weight_hits] else NULL
  metadata_cols <- mfrm_results_response_time_candidates(cols)
  facets <- setdiff(cols, c(person, score, weight, metadata_cols))
  if (length(facets) == 0L) {
    stop(
      "Could not infer facet columns after removing person, score, and weight. ",
      "Use fit_mfrm() explicitly or mfrm_results_interactive().",
      call. = FALSE
    )
  }

  list(person = person, facets = facets, score = score, weight = weight)
}

mfrm_results_mapping_table <- function(mapping) {
  if (is.null(mapping) || !is.list(mapping)) return(data.frame())
  data.frame(
    Key = c("Person", "Score", "Facets", "Weight"),
    Value = c(
      as.character(mapping$person %||% ""),
      as.character(mapping$score %||% ""),
      paste(as.character(mapping$facets %||% character(0)), collapse = ", "),
      as.character(mapping$weight %||% "")
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_results_resolve_input <- function(x) {
  if (inherits(x, "mfrm_facets_run")) {
    return(list(
      fit = x$fit,
      diagnostics = x$diagnostics,
      input_mode = "mfrm_facets_run",
      run = x,
      mapping = x$mapping %||% NULL,
      source_columns = x$fit$config$source_columns %||% x$mapping %||% NULL,
      source_data = NULL,
      notes = "Diagnostics were reused from the supplied run_mfrm_facets() object."
    ))
  }
  if (inherits(x, "mfrm_fit")) {
    return(list(
      fit = x,
      diagnostics = NULL,
      input_mode = "mfrm_fit",
      run = NULL,
      mapping = NULL,
      source_columns = x$config$source_columns %||% NULL,
      source_data = NULL,
      notes = character(0)
    ))
  }
  if (is.data.frame(x)) {
    mapping <- mfrm_results_infer_standard_mapping(x)
    run <- run_mfrm_facets(
      data = x,
      person = mapping$person,
      facets = mapping$facets,
      score = mapping$score,
      weight = mapping$weight
    )
    return(list(
      fit = run$fit,
      diagnostics = run$diagnostics,
      input_mode = "data.frame",
      run = run,
      mapping = mapping,
      source_columns = mapping,
      source_data = x,
      notes = paste0(
        "Input data.frame was estimated with run_mfrm_facets() using inferred columns: ",
        "person = ", mapping$person, ", score = ", mapping$score,
        ", facets = ", paste(mapping$facets, collapse = ", "), "."
      )
    ))
  }
  stop(
    "`fit` must be an mfrm_fit, mfrm_facets_run, or a standard long-format data.frame.",
    call. = FALSE
  )
}

mfrm_results_diagnose <- function(fit, diagnostics = NULL) {
  if (inherits(diagnostics, "mfrm_diagnostics")) {
    return(list(
      diagnostics = diagnostics,
      status = mfrm_results_status_row("diagnostics", "ok", "Reused supplied diagnostics.")
    ))
  }
  first <- mfrm_results_safe(
    diagnose_mfrm(
      fit,
      residual_pca = "none",
      diagnostic_mode = "both",
      fit_df_method = "both"
    )
  )
  if (isTRUE(first$ok)) {
    return(list(
      diagnostics = first$value,
      status = mfrm_results_status_row(
        "diagnostics", "ok",
        paste0(
          "Computed automatically with residual_pca = 'none', ",
          "diagnostic_mode = 'both', and fit_df_method = 'both'."
        )
      )
    ))
  }
  second <- mfrm_results_safe(diagnose_mfrm(fit, residual_pca = "none"))
  if (isTRUE(second$ok)) {
    return(list(
      diagnostics = second$value,
      status = mfrm_results_status_row(
        "diagnostics", "ok",
        paste0(
          "Computed automatically with residual_pca = 'none'. ",
          "The broader diagnostic_mode route was unavailable: ", first$message
        )
      )
    ))
  }
  list(
    diagnostics = NULL,
    status = mfrm_results_status_row(
      "diagnostics", "not_available",
      paste0("Automatic diagnostics failed: ", second$message)
    )
  )
}

mfrm_results_flatten_data_frames <- function(x, prefix, max_depth = 2L, depth = 0L) {
  if (is.data.frame(x)) {
    out <- list(as.data.frame(x, stringsAsFactors = FALSE))
    names(out) <- prefix
    return(out)
  }
  if (!is.list(x) || depth >= max_depth) return(list())
  nm <- names(x)
  if (is.null(nm)) nm <- paste0("component", seq_along(x))
  pieces <- list()
  for (i in seq_along(x)) {
    child <- x[[i]]
    child_name <- paste(prefix, nm[[i]], sep = "_")
    child_pieces <- mfrm_results_flatten_data_frames(
      child,
      prefix = child_name,
      max_depth = max_depth,
      depth = depth + 1L
    )
    if (length(child_pieces) > 0L) {
      pieces <- c(pieces, child_pieces)
    }
  }
  pieces
}

mfrm_results_add_component <- function(name, result, components, tables, status,
                                       table_prefix = name) {
  if (isTRUE(result$ok)) {
    components[[name]] <- result$value
    status <- rbind(
      status,
      mfrm_results_status_row(name, "ok", "Available.")
    )
    tbls <- mfrm_results_flatten_data_frames(result$value, table_prefix)
    if (length(tbls) > 0L) {
      tables <- c(tables, tbls)
    }
  } else {
    status <- rbind(
      status,
      mfrm_results_status_row(name, "not_available", result$message)
    )
  }
  list(components = components, tables = tables, status = status)
}

mfrm_results_bias_screen_bundle <- function(fit, diagnostics) {
  bias_tbl <- as.data.frame(diagnostics$bias %||% data.frame(), stringsAsFactors = FALSE)
  facet_names <- as.character(fit$config$facet_names %||% character(0))
  guidance <- data.frame(
    Area = c("Facet-level bias screen", "Interaction bias screen", "Reporting boundary"),
    Status = c(
      if (nrow(bias_tbl) > 0L) "available" else "not_available",
      if (length(facet_names) >= 2L) "requires_explicit_facet_pair" else "not_available",
      "screening_not_final_fairness_decision"
    ),
    Route = c(
      "diagnose_mfrm(fit)$bias",
      "estimate_bias(fit, diagnostics, facet_a = ..., facet_b = ...) -> bias_interaction_report()",
      "Use reporting_checklist() and substantive review before fairness claims."
    ),
    Detail = c(
      if (nrow(bias_tbl) > 0L) "Facet-level observed-minus-expected bias screen is available." else "No facet-level bias table is available.",
      if (length(facet_names) >= 2L) paste0("Available non-person facets: ", paste(facet_names, collapse = ", "), ". The results wrapper does not choose the contrast automatically.") else "At least two non-person facets are needed for an interaction-bias screen.",
      "Bias outputs are conditional screening layers, not standalone validity or fairness decisions."
    ),
    stringsAsFactors = FALSE
  )
  out <- list(
    table = bias_tbl,
    guidance = guidance,
    available_facets = data.frame(Facet = facet_names, stringsAsFactors = FALSE)
  )
  class(out) <- c("mfrm_bias_screen", "list")
  out
}

mfrm_results_misfit_review_bundle <- function(fit, diagnostics, top_n = 50L) {
  unexpected <- unexpected_response_table(fit, diagnostics = diagnostics, top_n = top_n)
  displacement <- displacement_table(fit, diagnostics = diagnostics, top_n = top_n)
  pathway <- plot(fit, type = "pathway", draw = FALSE)
  pathway_data <- pathway$data %||% list()
  guidance <- data.frame(
    Area = c("Unexpected responses", "Displacement", "Pathway map", "Interpretation boundary"),
    Status = c(
      if (nrow(as.data.frame(unexpected$table %||% data.frame())) > 0L) "available" else "empty",
      if (nrow(as.data.frame(displacement$table %||% data.frame())) > 0L) "available" else "empty",
      "available",
      "case_review_not_validity_decision"
    ),
    Route = c(
      "unexpected_response_table(fit, diagnostics)",
      "displacement_table(fit, diagnostics)",
      "plot(fit, type = \"pathway\")",
      "Use build_misfit_casebook() or substantive review before person/rater/item conclusions."
    ),
    Detail = c(
      "Observation-level surprising rows with residual, probability, and direction metadata.",
      "Facet/person displacement review for large measure changes or anchor-related shifts.",
      "Expected-score pathway with fit annotations for custom inspection.",
      "Misfit rows are follow-up prompts; they are not by themselves exclusion or fairness decisions."
    ),
    stringsAsFactors = FALSE
  )
  out <- list(
    unexpected = unexpected,
    displacement = displacement,
    pathway_fit_measures = as.data.frame(pathway_data$fit_measures %||% data.frame(), stringsAsFactors = FALSE),
    pathway_fit_status = as.data.frame(pathway_data$fit_status %||% data.frame(), stringsAsFactors = FALSE),
    pathway_curve_fit_status = as.data.frame(pathway_data$curve_fit_status %||% data.frame(), stringsAsFactors = FALSE),
    guidance = guidance
  )
  class(out) <- c("mfrm_misfit_review", "list")
  out
}

mfrm_results_linking_review_bundle <- function(fit) {
  anchor_review <- fit$config$anchor_review %||% NULL
  if (!inherits(anchor_review, "mfrm_anchor_review")) {
    stop(
      "Anchor-review metadata is not available in this fit. ",
      "Run review_mfrm_anchors() directly with the source data, or refit with fit_mfrm().",
      call. = FALSE
    )
  }

  review <- build_linking_review(anchor_review = anchor_review)
  model <- toupper(as.character(fit$config$model %||% ""))
  support_status <- as.data.frame(review$support_status %||% data.frame(), stringsAsFactors = FALSE)
  gpcm_status <- if (nrow(support_status) > 0L &&
      all(c("Scope", "Status") %in% names(support_status))) {
    as.character(support_status$Status[support_status$Scope %in% "bounded GPCM"][1] %||% "")
  } else {
    ""
  }
  review$first_screen_guidance <- data.frame(
    Area = c(
      "Anchor readiness",
      "Drift review",
      "Equating chain",
      "GPCM boundary",
      "Reporting boundary"
    ),
    Status = c(
      "available",
      "requires_multiple_fits",
      "requires_ordered_fit_list",
      if (identical(model, "GPCM")) gpcm_status else "rsm_pcm_route",
      "operational_review_not_validation_decision"
    ),
    Route = c(
      "summary(res$components$linking_review); plot(res, type = \"anchors\")",
      "detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2)); plot_anchor_drift(drift)",
      "build_equating_chain(list(Form1 = fit1, Form2 = fit2)); plot_anchor_drift(chain, type = \"chain\")",
      "gpcm_capability_matrix()",
      "build_summary_table_bundle(res$components$linking_review)"
    ),
    Detail = c(
      "The fitted object's stored anchor review is available for first-screen linking readiness.",
      "Drift checks compare at least two separately fitted calibrations and are not inferred from one fit.",
      "Screened chain offsets require an ordered list of fitted forms or waves.",
      if (identical(model, "GPCM")) {
        "Bounded GPCM linking synthesis is available as caveated exploratory review; do not treat it as an operational linking decision or anchor-drift absence claim."
      } else {
        "RSM/PCM operational linking review is the supported synthesis route."
      },
      "Linking outputs are operational support and scale-maintenance prompts, not standalone validity decisions."
    ),
    stringsAsFactors = FALSE
  )
  review
}

mfrm_results_response_time_bundle <- function(ctx,
                                              response_time = NULL,
                                              response_time_data = NULL,
                                              response_time_facets = NULL,
                                              response_time_score = NULL) {
  fit <- ctx$fit
  source_columns <- ctx$source_columns %||% fit$config$source_columns %||% ctx$mapping %||% list()
  data <- response_time_data %||% ctx$source_data %||% NULL
  if (!is.data.frame(data)) {
    prep_data <- fit$prep$data %||% NULL
    if (is.data.frame(prep_data) && !is.null(response_time) &&
        as.character(response_time[1]) %in% names(prep_data)) {
      data <- prep_data
    }
  }
  time_col <- mfrm_results_resolve_response_time_column(data, response_time = response_time)
  person_col <- as.character(source_columns$person %||% "Person")
  if (!person_col %in% names(data) && "Person" %in% names(data)) {
    person_col <- "Person"
  }
  if (!person_col %in% names(data)) {
    stop("Person column for response-time review not found: ", person_col, call. = FALSE)
  }
  facet_cols <- response_time_facets %||%
    source_columns$facets %||%
    fit$config$facet_names %||%
    character(0)
  facet_cols <- as.character(facet_cols)
  facet_cols <- facet_cols[nzchar(facet_cols) & facet_cols %in% names(data)]
  facet_cols <- setdiff(facet_cols, time_col)
  if (length(facet_cols) == 0L) facet_cols <- NULL
  score_col <- response_time_score %||% source_columns$score %||% "Score"
  score_col <- as.character(score_col[1] %||% "")
  if (!nzchar(score_col) || !score_col %in% names(data) || identical(score_col, time_col)) {
    score_col <- NULL
  }
  review <- response_time_review(
    data = data,
    person = person_col,
    facets = facet_cols,
    time = time_col,
    score = score_col
  )
  review$guidance <- data.frame(
    Area = c("Timing data", "Interpretation boundary", "Reusable plot data"),
    Status = c("available", "descriptive_qc_only", "available"),
    Route = c(
      "summary(res$components$response_time_review)",
      "response_time_review(); fit_mfrm() remains unchanged",
      "plot(res, type = \"response_time\", draw = FALSE); plot_data_components()"
    ),
    Detail = c(
      paste0("Response-time column `", time_col, "` was reviewed outside the fitted MFRM likelihood."),
      "Response-time summaries do not estimate speed parameters, modify logits, or create automatic exclusion rules.",
      "Distribution, person, facet, and score views can be extracted as mfrm_plot_data objects."
    ),
    stringsAsFactors = FALSE
  )
  review$config$results_route <- data.frame(
    TimeColumn = time_col,
    PersonColumn = person_col,
    FacetColumns = paste(facet_cols %||% character(0), collapse = ", "),
    ScoreColumn = as.character(score_col %||% ""),
    Source = if (!is.null(response_time_data)) "response_time_data" else if (!is.null(ctx$source_data)) "data_frame_input" else "fit_prep_data",
    stringsAsFactors = FALSE
  )
  review
}

mfrm_results_table_index <- function(tables) {
  if (!is.list(tables) || length(tables) == 0L) return(data.frame())
  rows <- lapply(names(tables), function(nm) {
    tbl <- as.data.frame(tables[[nm]], stringsAsFactors = FALSE)
    numeric_cols <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
    data.frame(
      Table = nm,
      Rows = nrow(tbl),
      Cols = ncol(tbl),
      NumericColumns = length(numeric_cols),
      PlotReady = nrow(tbl) > 0L && length(numeric_cols) > 0L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mfrm_results_plot_map <- function(has_fit, has_diagnostics, tables, components) {
  table_names <- names(tables %||% list())
  data.frame(
    Type = c("fit", "wright", "pathway", "qc", "category", "anchors", "response_time", "tables"),
    Available = c(
      isTRUE(has_fit),
      isTRUE(has_fit),
      isTRUE(has_fit),
      isTRUE(has_fit) && isTRUE(has_diagnostics),
      "rating_scale" %in% names(components) || any(grepl("^rating_scale", table_names)),
      "linking_review" %in% names(components),
      "response_time_review" %in% names(components),
      length(table_names) > 0L
    ),
    Route = c(
      "plot(res, type = 'fit')",
      "plot(res, type = 'wright')",
      "plot(res, type = 'pathway')",
      "plot(res, type = 'qc')",
      "plot(res, type = 'category')",
      "plot(res, type = 'anchors')",
      "plot(res, type = 'response_time')",
      "plot(res, type = 'tables')"
    ),
    Detail = c(
      "Model-level visual bundle from plot.mfrm_fit().",
      "Wright map from plot.mfrm_fit().",
      "Expected-score pathway map from plot.mfrm_fit().",
      "Quality-control dashboard from plot_qc_dashboard().",
      "Rating-scale/category plot when rating_scale_table() is available.",
      "Anchor-review plot from the stored fit_mfrm() anchor review.",
      "Descriptive response-time QC plot when response_time_review() is available.",
      "Numeric table-profile plot from the summary-table bundle."
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_results_deparse_one <- function(x) {
  paste(deparse(x), collapse = "")
}

mfrm_results_include_expr <- function(include) {
  requested <- attr(include, "requested", exact = TRUE)
  if (!is.null(requested) && length(requested) > 0L) {
    include <- requested
  }
  include <- as.character(include %||% "standard")
  if (length(include) == 1L) return(mfrm_results_deparse_one(include))
  paste0("c(", paste(vapply(include, mfrm_results_deparse_one, character(1)), collapse = ", "), ")")
}

mfrm_results_reproducible_response_time_args <- function(ctx) {
  rt <- ctx$response_time %||% NULL
  if (!is.list(rt) || !isTRUE(rt$requested)) return(character(0))
  lines <- character(0)
  if (!is.null(rt$time) && length(rt$time) > 0L &&
      !is.na(rt$time[1]) && nzchar(rt$time[1])) {
    lines <- c(lines, paste0("  response_time = ", mfrm_results_deparse_one(as.character(rt$time[1])), ","))
  }
  if (isTRUE(rt$data_supplied)) {
    lines <- c(lines, "  response_time_data = response_time_data,")
  } else if (identical(ctx$input_mode, "data.frame") && !is.null(ctx$source_data)) {
    lines <- c(lines, "  response_time_data = data,")
  }
  if (!is.null(rt$facets) && length(rt$facets) > 0L) {
    facet_expr <- paste(deparse(as.character(rt$facets)), collapse = "")
    lines <- c(lines, paste0("  response_time_facets = ", facet_expr, ","))
  }
  if (!is.null(rt$score) && length(rt$score) > 0L &&
      !is.na(rt$score[1]) && nzchar(rt$score[1])) {
    lines <- c(lines, paste0("  response_time_score = ", mfrm_results_deparse_one(as.character(rt$score[1])), ","))
  }
  lines
}

mfrm_results_reproducible_code <- function(ctx, include, output = "object") {
  include_expr <- mfrm_results_include_expr(include)
  output_expr <- mfrm_results_deparse_one(output)
  mapping <- ctx$mapping %||% NULL
  rt_lines <- mfrm_results_reproducible_response_time_args(ctx)
  if (is.list(mapping) && !is.null(mapping$person) && !is.null(mapping$score) &&
      length(mapping$facets %||% character(0)) > 0L) {
    fit_code <- mfrm_results_render_code(
      person = mapping$person,
      facets = as.character(mapping$facets),
      score = mapping$score,
      weight = mapping$weight %||% NULL,
      include = include,
      output = output,
      response_time_lines = rt_lines
    )
    return(fit_code)
  }
  lines <- c(
    "res <- mfrm_results(",
    "  fit = fit,",
    paste0("  include = ", include_expr, ","),
    rt_lines,
    paste0("  output = ", output_expr),
    ")"
  )
  paste(lines, collapse = "\n")
}

mfrm_results_code_table <- function(code) {
  code <- as.character(code %||% "")
  if (!nzchar(code)) return(data.frame())
  lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  data.frame(
    Line = seq_along(lines),
    Code = lines,
    stringsAsFactors = FALSE
  )
}

mfrm_results_status_detail <- function(status, section) {
  status <- as.data.frame(status %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(status) == 0L || !all(c("Section", "Status") %in% names(status))) return(NULL)
  row <- status[status$Section %in% section, , drop = FALSE]
  if (nrow(row) == 0L) return(NULL)
  row[1, , drop = FALSE]
}

mfrm_results_summary_key_warnings <- function(diagnostics, summaries) {
  diag_sum <- summaries$diagnostics %||% NULL
  if (is.null(diag_sum) && inherits(diagnostics, "mfrm_diagnostics")) {
    maybe <- mfrm_results_safe(summary(diagnostics))
    if (isTRUE(maybe$ok)) diag_sum <- maybe$value
  }
  key <- as.character(diag_sum$key_warnings %||% character(0))
  key <- key[nzchar(key)]
  key
}

mfrm_results_triage <- function(status, plot_map, components, table_index,
                                fit, diagnostics, summaries) {
  status <- as.data.frame(status %||% data.frame(), stringsAsFactors = FALSE)
  plot_map <- as.data.frame(plot_map %||% data.frame(), stringsAsFactors = FALSE)
  table_index <- as.data.frame(table_index %||% data.frame(), stringsAsFactors = FALSE)
  components <- components %||% list()
  summaries <- summaries %||% list()
  rows <- list()
  add <- function(area, severity, signal, route, detail) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Area = as.character(area),
      Severity = as.character(severity),
      Signal = as.character(signal),
      Route = as.character(route),
      Detail = as.character(detail %||% ""),
      stringsAsFactors = FALSE
    )
  }

  key_warnings <- mfrm_results_summary_key_warnings(diagnostics, summaries)
  no_diag_warning <- length(key_warnings) == 0L ||
    identical(key_warnings[1], "No immediate warnings from diagnostics summary.")
  if (inherits(diagnostics, "mfrm_diagnostics")) {
    add(
      "Diagnostics",
      if (isTRUE(no_diag_warning)) "ok" else "review",
      if (isTRUE(no_diag_warning)) "diagnostics_available" else "diagnostic_warnings_present",
      "summary(res$diagnostics)$key_warnings",
      if (isTRUE(no_diag_warning)) "Diagnostics are available and no immediate summary warning was reported." else paste(utils::head(key_warnings, 2L), collapse = " | ")
    )
  } else {
    add(
      "Diagnostics",
      "not_available",
      "diagnostics_missing",
      "diagnose_mfrm(res$fit, residual_pca = \"none\")",
      "Automatic diagnostics were not available; dependent result sections should be read as unavailable."
    )
  }

  not_available <- if (nrow(status) > 0L && "Status" %in% names(status)) {
    sum(status$Status %in% "not_available", na.rm = TRUE)
  } else {
    0L
  }
  add(
    "Section availability",
    if (not_available > 0L) "review" else "ok",
    if (not_available > 0L) "some_sections_unavailable" else "requested_sections_available",
    "summary(res)$status",
    if (not_available > 0L) paste0(not_available, " requested section(s) were not available; review status before treating omissions as evidence.") else "Requested sections that could be computed were available."
  )

  qc_available <- nrow(plot_map) > 0L &&
    any(plot_map$Type %in% "qc" & plot_map$Available %in% TRUE)
  add(
    "Visual screen",
    if (qc_available) "ok" else "not_available",
    if (qc_available) "qc_plot_available" else "qc_plot_missing",
    if (qc_available) "plot(res, type = \"qc\", preset = \"publication\")" else "summary(res)$plot_map",
    if (qc_available) "The QC dashboard route is available for first visual review." else "QC plotting is unavailable for this result object."
  )

  add(
    "Tables",
    if (nrow(table_index) > 0L) "ok" else "review",
    if (nrow(table_index) > 0L) "tables_collected" else "no_tables_collected",
    "build_summary_table_bundle(res)",
    if (nrow(table_index) > 0L) paste0(nrow(table_index), " data-frame table(s) were collected for appendix or handoff use.") else "No table-like outputs were collected."
  )

  precision_status <- mfrm_results_status_detail(status, "precision_review")
  if ("precision_review" %in% names(components)) {
    add(
      "Precision / separation",
      "review",
      "precision_review_available",
      "summary(res$components$precision_review)",
      "Precision review is available; inspect fit, separation, reliability, and ZSTD wording boundaries before reporting."
    )
  } else if (!is.null(precision_status) && precision_status$Status[1] %in% "not_available") {
    add(
      "Precision / separation",
      "not_available",
      "precision_review_missing",
      "summary(res)$status",
      as.character(precision_status$Detail[1] %||% "Precision review was not available.")
    )
  }

  if ("reporting_checklist" %in% names(components)) {
    add(
      "Reporting",
      "ok",
      "reporting_checklist_available",
      "summary(res$components$reporting_checklist)",
      "Reporting checklist is available as the manuscript-routing surface."
    )
  }

  if ("bias_screen" %in% names(components)) {
    add(
      "Bias screening",
      "review",
      "bias_screen_available",
      "res$components$bias_screen$guidance",
      "Bias screening is available; interaction-bias review still requires an explicit facet-pair choice."
    )
  }

  if ("misfit_review" %in% names(components)) {
    add(
      "Pathway / misfit",
      "review",
      "misfit_review_available",
      "res$components$misfit_review$guidance",
      "Unexpected-response, displacement, and pathway evidence are case-review prompts, not automatic exclusion decisions."
    )
  }

  if ("linking_review" %in% names(components)) {
    add(
      "Linking / anchors",
      "review",
      "linking_review_available",
      "summary(res$components$linking_review); plot(res, type = \"anchors\")",
      "Anchor readiness is available from the stored fit review; drift and screened-chain checks require multiple fitted calibrations."
    )
  }

  if ("response_time_review" %in% names(components)) {
    add(
      "Response-time QC",
      "info",
      "response_time_review_available",
      "summary(res$components$response_time_review); plot(res, type = \"response_time\")",
      "Response-time metadata were reviewed descriptively outside the fitted MFRM likelihood."
    )
  }

  if ("network_review" %in% names(components)) {
    add(
      "Network / connectivity",
      "info",
      "network_review_available",
      "summary(res$components$network_review)",
      "Network review is design/connectivity evidence and does not replace fit, separation, or fairness diagnostics."
    )
  }

  model <- as.character(fit$config$model %||% "")
  if (identical(toupper(model), "GPCM")) {
    add(
      "Model scope",
      "info",
      "bounded_gpcm_scope",
      "gpcm_capability_matrix()",
      "Bounded GPCM results should be read through the documented capability matrix and caveats."
    )
  }

  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  severity_rank <- c(not_available = 1L, review = 2L, info = 3L, ok = 4L)
  ord <- severity_rank[out$Severity]
  ord[is.na(ord)] <- 5L
  out[order(ord, out$Area), , drop = FALSE]
}

mfrm_results_next_actions <- function(status, plot_map, components, table_index,
                                      triage = NULL) {
  status <- as.data.frame(status %||% data.frame(), stringsAsFactors = FALSE)
  plot_map <- as.data.frame(plot_map %||% data.frame(), stringsAsFactors = FALSE)
  table_index <- as.data.frame(table_index %||% data.frame(), stringsAsFactors = FALSE)
  triage <- as.data.frame(triage %||% data.frame(), stringsAsFactors = FALSE)
  components <- components %||% list()
  rows <- list()
  add <- function(priority, area, action, route, reason) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Priority = as.integer(priority),
      Area = as.character(area),
      Action = as.character(action),
      Route = as.character(route),
      Reason = as.character(reason),
      stringsAsFactors = FALSE
    )
  }

  add(
    1L,
    "Overview",
    "Read the compact results summary.",
    "summary(res)",
    "Confirms input mode, model, method, section status, table coverage, and plot routes."
  )
  if (nrow(triage) > 0L &&
      "Severity" %in% names(triage) &&
      any(triage$Severity %in% c("not_available", "review"), na.rm = TRUE)) {
    add(
      2L,
      "Triage",
      "Read the first-screen triage before branching.",
      "summary(res)$triage",
      "Triage orders unavailable, review, information, and OK signals across diagnostics, tables, plots, and reporting surfaces."
    )
  }
  if (nrow(triage) > 0L &&
      all(c("Area", "Severity") %in% names(triage)) &&
      any(triage$Area %in% "Diagnostics" & triage$Severity %in% "review", na.rm = TRUE)) {
    add(
      3L,
      "Diagnostics",
      "Review diagnostic key warnings before report drafting.",
      "summary(res$diagnostics)$key_warnings",
      "Diagnostic warnings identify the highest-priority fit, precision, residual, or category follow-up surfaces."
    )
  }
  if (nrow(plot_map) > 0L && any(plot_map$Type %in% "qc" & plot_map$Available %in% TRUE)) {
    add(
      4L,
      "Visual diagnostics",
      "Open the QC dashboard before drilling into individual tables.",
      "plot(res, type = \"qc\", preset = \"publication\")",
      "The QC route gives a first visual check of fit, residual, and category surfaces."
    )
  }
  if ("precision_review" %in% names(components)) {
    add(
      5L,
      "Precision",
      "Inspect fit, separation, reliability, and ZSTD wording boundaries.",
      "summary(res$components$precision_review)",
      "Precision review keeps fit-size, standardized fit, and separation evidence in separate reporting lanes."
    )
  }
  if ("reporting_checklist" %in% names(components)) {
    add(
      6L,
      "Reporting",
      "Use the reporting checklist as the manuscript-routing surface.",
      "summary(res$components$reporting_checklist)",
      "Checklist rows identify report-ready, missing, and caveated sections."
    )
  }
  if ("bias_screen" %in% names(components)) {
    add(
      7L,
      "Bias screening",
      "Review facet-level bias screens and choose any interaction contrast explicitly.",
      "res$components$bias_screen$guidance",
      "The wrapper surfaces bias evidence but does not choose a facet pair for interaction-bias claims."
    )
  }
  if ("misfit_review" %in% names(components)) {
    add(
      8L,
      "Pathway / misfit",
      "Inspect unexpected responses alongside pathway and displacement evidence.",
      "res$components$misfit_review",
      "Misfit rows are follow-up prompts and should be read with pathway and substantive context."
    )
  }
  if ("linking_review" %in% names(components)) {
    add(
      9L,
      "Linking / anchors",
      "Inspect anchor readiness before making linking, drift, or DFF claims.",
      "summary(res$components$linking_review); plot(res, type = \"anchors\")",
      "The first-screen linking review uses the stored anchor review; drift and equating checks still need an explicit list of fitted waves or forms."
    )
  }
  if ("response_time_review" %in% names(components)) {
    add(
      10L,
      "Response-time QC",
      "Inspect timing summaries separately from model estimates.",
      "summary(res$components$response_time_review); plot(res, type = \"response_time\", draw = FALSE)",
      "Timing summaries are descriptive QC context and do not change fit_mfrm estimates or define exclusion rules."
    )
  }
  if (nrow(table_index) > 0L) {
    add(
      11L,
      "Tables",
      "Create an appendix-ready summary-table bundle.",
      "build_summary_table_bundle(res)",
      "The bundle exposes table roles, plot readiness, and conservative appendix presets."
    )
  }
  if (nrow(status) > 0L && any(status$Status %in% "not_available")) {
    add(
      12L,
      "Availability",
      "Review unavailable sections before interpreting missing output as evidence.",
      "summary(res)$status",
      "Unavailable sections usually reflect model scope, missing dependencies, or insufficient data rather than a psychometric result."
    )
  }
  if (nrow(triage) > 0L &&
      all(c("Area", "Signal") %in% names(triage)) &&
      any(triage$Area %in% "Model scope" & triage$Signal %in% "bounded_gpcm_scope", na.rm = TRUE)) {
    add(
      13L,
      "GPCM scope",
      "Check bounded-GPCM support boundaries before interpreting advanced outputs.",
      "gpcm_capability_matrix()",
      "The GPCM route is supported with documented helper-specific caveats and should not be treated as a universal Rasch-family replacement."
    )
  }
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(out$Priority, out$Area), , drop = FALSE]
}

mfrm_results_build <- function(ctx, include) {
  fit <- ctx$fit
  status <- mfrm_results_status_row("input", "ok", paste0("Input mode: ", ctx$input_mode, "."))
  components <- list(fit = fit)
  summaries <- list()
  tables <- list()
  notes <- as.character(ctx$notes %||% character(0))

  diag_info <- mfrm_results_diagnose(fit, diagnostics = ctx$diagnostics)
  diagnostics <- diag_info$diagnostics
  status <- rbind(status, diag_info$status)
  if (inherits(diagnostics, "mfrm_diagnostics")) {
    components$diagnostics <- diagnostics
  } else {
    notes <- c(notes, "Diagnostics were not available; dependent sections are omitted.")
  }

  if ("fit" %in% include) {
    fit_sum <- mfrm_results_safe(summary(fit))
    if (isTRUE(fit_sum$ok)) {
      summaries$fit <- fit_sum$value
      tables <- c(tables, mfrm_results_flatten_data_frames(fit_sum$value, "fit_summary"))
      status <- rbind(status, mfrm_results_status_row("fit_summary", "ok", "Available."))
    } else {
      status <- rbind(status, mfrm_results_status_row("fit_summary", "not_available", fit_sum$message))
    }
  }

  if ("diagnostics" %in% include && inherits(diagnostics, "mfrm_diagnostics")) {
    diag_sum <- mfrm_results_safe(summary(diagnostics))
    if (isTRUE(diag_sum$ok)) {
      summaries$diagnostics <- diag_sum$value
      tables <- c(tables, mfrm_results_flatten_data_frames(diag_sum$value, "diagnostics_summary"))
      status <- rbind(status, mfrm_results_status_row("diagnostics_summary", "ok", "Available."))
    } else {
      status <- rbind(status, mfrm_results_status_row("diagnostics_summary", "not_available", diag_sum$message))
    }
  }

  if ("tables" %in% include || "categories" %in% include) {
    table_calls <- list(
      iteration = quote(estimation_iteration_report(fit)),
      fit_measures = quote(fit_measures_table(
        fit,
        diagnostics = diagnostics,
        threshold_profiles = "all",
        fit_df_method = "both"
      )),
      facet_statistics = quote(facet_statistics_report(fit, diagnostics = diagnostics)),
      fair_average = quote(fair_average_table(fit, diagnostics = diagnostics)),
      rating_scale = quote(rating_scale_table(fit, diagnostics = diagnostics)),
      unexpected = quote(unexpected_response_table(fit, diagnostics = diagnostics, top_n = 20))
    )
    if (!inherits(diagnostics, "mfrm_diagnostics")) {
      table_calls <- table_calls["iteration"]
    }
    for (nm in names(table_calls)) {
      result <- mfrm_results_safe(eval(table_calls[[nm]]))
      added <- mfrm_results_add_component(nm, result, components, tables, status, table_prefix = nm)
      components <- added$components
      tables <- added$tables
      status <- added$status
    }
  }

  if ("precision" %in% include && inherits(diagnostics, "mfrm_diagnostics")) {
    result <- mfrm_results_safe(precision_review_report(fit, diagnostics = diagnostics))
    added <- mfrm_results_add_component("precision_review", result, components, tables, status)
    components <- added$components
    tables <- added$tables
    status <- added$status
  }

  if ("facets_fit" %in% include && inherits(diagnostics, "mfrm_diagnostics")) {
    result <- mfrm_results_safe(facets_fit_review(fit, diagnostics = diagnostics))
    added <- mfrm_results_add_component("facets_fit_review", result, components, tables, status)
    components <- added$components
    tables <- added$tables
    status <- added$status
  }

  if ("bias" %in% include && inherits(diagnostics, "mfrm_diagnostics")) {
    result <- mfrm_results_safe(mfrm_results_bias_screen_bundle(fit, diagnostics))
    added <- mfrm_results_add_component("bias_screen", result, components, tables, status)
    components <- added$components
    tables <- added$tables
    status <- added$status
  }

  if ("misfit" %in% include && inherits(diagnostics, "mfrm_diagnostics")) {
    result <- mfrm_results_safe(mfrm_results_misfit_review_bundle(fit, diagnostics, top_n = 50L))
    added <- mfrm_results_add_component("misfit_review", result, components, tables, status)
    components <- added$components
    tables <- added$tables
    status <- added$status
  }

  if ("linking" %in% include) {
    result <- mfrm_results_safe(mfrm_results_linking_review_bundle(fit))
    added <- mfrm_results_add_component("linking_review", result, components, tables, status)
    components <- added$components
    tables <- added$tables
    status <- added$status
  }

  if ("reporting" %in% include && inherits(diagnostics, "mfrm_diagnostics")) {
    result <- mfrm_results_safe(reporting_checklist(fit, diagnostics = diagnostics))
    added <- mfrm_results_add_component("reporting_checklist", result, components, tables, status)
    components <- added$components
    tables <- added$tables
    status <- added$status
  }

  if ("apa" %in% include && inherits(diagnostics, "mfrm_diagnostics")) {
    result <- mfrm_results_safe(build_apa_outputs(fit, diagnostics = diagnostics))
    added <- mfrm_results_add_component("apa_outputs", result, components, tables, status)
    components <- added$components
    tables <- added$tables
    status <- added$status
  }

  if ("network" %in% include && inherits(diagnostics, "mfrm_diagnostics")) {
    result <- mfrm_results_safe(build_mfrm_network_review(fit, diagnostics = diagnostics, include_graph = FALSE))
    added <- mfrm_results_add_component("network_review", result, components, tables, status)
    components <- added$components
    tables <- added$tables
    status <- added$status
  }

  if ("response_time" %in% include) {
    rt <- ctx$response_time %||% list()
    result <- mfrm_results_safe(mfrm_results_response_time_bundle(
      ctx = ctx,
      response_time = rt$time %||% NULL,
      response_time_data = rt$data %||% NULL,
      response_time_facets = rt$facets %||% NULL,
      response_time_score = rt$score %||% NULL
    ))
    added <- mfrm_results_add_component("response_time_review", result, components, tables, status)
    components <- added$components
    tables <- added$tables
    status <- added$status
  }

  table_index <- mfrm_results_table_index(tables)
  plot_map <- mfrm_results_plot_map(
    has_fit = inherits(fit, "mfrm_fit"),
    has_diagnostics = inherits(diagnostics, "mfrm_diagnostics"),
    tables = tables,
    components = components
  )
  triage <- mfrm_results_triage(
    status = status,
    plot_map = plot_map,
    components = components,
    table_index = table_index,
    fit = fit,
    diagnostics = diagnostics,
    summaries = summaries
  )
  next_actions <- mfrm_results_next_actions(
    status = status,
    plot_map = plot_map,
    components = components,
    table_index = table_index,
    triage = triage
  )
  reproducible_code <- mfrm_results_reproducible_code(
    ctx = ctx,
    include = include,
    output = "object"
  )

  out <- list(
    fit = fit,
    diagnostics = diagnostics,
    components = components,
    summaries = summaries,
    tables = tables,
    table_index = table_index,
    plot_map = plot_map,
    triage = triage,
    next_actions = next_actions,
    status = status,
    include = include,
    input = list(
      mode = ctx$input_mode,
      mapping = ctx$mapping %||% NULL,
      reproducible_code = reproducible_code
    ),
    notes = unique(notes[nzchar(notes)])
  )
  class(out) <- "mfrm_results"
  out
}

mfrm_results_html <- function(x) {
  summary_obj <- summary(x)
  html_tables <- c(
    list(
      overview = summary_obj$overview,
      status = summary_obj$status,
      component_index = summary_obj$component_index,
      table_index = summary_obj$table_index,
      plot_map = summary_obj$plot_map,
      triage = summary_obj$triage,
      next_actions = summary_obj$next_actions,
      reproducible_code = summary_obj$reproducible_code
    ),
    x$tables
  )
  text_sections <- list()
  if (length(summary_obj$notes %||% character(0)) > 0L) {
    text_sections$notes <- paste(summary_obj$notes, collapse = "\n")
  }
  html <- build_mfrm_bundle_html(
    title = "mfrmr Results",
    tables = html_tables,
    text_sections = text_sections
  )
  path <- tempfile("mfrmr_results_", fileext = ".html")
  writeLines(enc2utf8(html), con = path, useBytes = TRUE)
  out <- list(
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    summary = summary_obj,
    html = html
  )
  class(out) <- "mfrm_results_html"
  out
}

mfrm_report_title <- function(style) {
  switch(
    style,
    apa = "mfrmr APA Reporting Template",
    qc = "mfrmr QC Report",
    validation = "mfrmr Validation Report",
    reviewer = "mfrmr Reviewer Report",
    technical = "mfrmr Technical Report",
    "mfrmr Report"
  )
}

mfrm_report_style_focus <- function(style) {
  switch(
    style,
    apa = "Manuscript-ready wording scaffold with explicit evidence boundaries.",
    qc = "Quality-control triage before manuscript, appendix, or reviewer handoff.",
    validation = "Validity-argument evidence map with limits on what each output can support.",
    reviewer = "Reviewer-facing response map that separates checked evidence from caveats.",
    technical = "Technical appendix map for tables, routes, reproducibility, and diagnostics.",
    "Report-ready synthesis from mfrm_results()."
  )
}

mfrm_report_component_status <- function(x, component, available = "available",
                                         absent = "not_requested") {
  if (component %in% names(x$components %||% list())) return(available)
  status <- as.data.frame(x$status %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(status) > 0L && all(c("Section", "Status") %in% names(status))) {
    hit <- status[status$Section %in% component, , drop = FALSE]
    if (nrow(hit) > 0L && hit$Status[1] %in% "not_available") {
      return("not_available")
    }
  }
  absent
}

mfrm_report_has_category <- function(x) {
  comps <- names(x$components %||% list())
  tbls <- names(x$tables %||% list())
  any(comps %in% c("rating_scale", "category_structure", "category_curves")) ||
    any(grepl("rating_scale|category", tbls))
}

mfrm_report_plot_route <- function(x, type) {
  plot_map <- as.data.frame(x$plot_map %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(plot_map) == 0L || !all(c("Type", "Available", "Route") %in% names(plot_map))) {
    return("")
  }
  row <- plot_map[plot_map$Type %in% type & plot_map$Available %in% TRUE, , drop = FALSE]
  if (nrow(row) == 0L) "" else as.character(row$Route[1])
}

mfrm_report_section_plan <- function(x, sx, style) {
  overview <- as.data.frame(sx$overview %||% data.frame(), stringsAsFactors = FALSE)
  model <- as.character(overview$Model[1] %||% x$fit$config$model %||% "")
  method <- as.character(overview$Method[1] %||% x$fit$config$method %||% "")
  n_obs <- as.character(overview$N[1] %||% "")
  categories <- as.character(overview$Categories[1] %||% "")
  table_count <- length(x$tables %||% list())
  plot_count <- sum(as.data.frame(x$plot_map %||% data.frame())$Available %in% TRUE, na.rm = TRUE)
  rows <- list()
  add <- function(section, status, evidence, route, report_use, boundary) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Section = as.character(section),
      Status = as.character(status),
      Evidence = as.character(evidence %||% ""),
      Route = as.character(route %||% ""),
      ReportUse = as.character(report_use %||% ""),
      Boundary = as.character(boundary %||% ""),
      stringsAsFactors = FALSE
    )
  }

  add(
    "Model and data setup",
    "available",
    paste0("Model = ", model, "; method = ", method, "; N = ", n_obs,
           "; categories = ", categories, "."),
    "summary(res)$overview; specifications_report(fit)",
    "Use for method and analysis-setup wording.",
    "Confirm scoring, column roles, missing-data handling, anchoring, and estimation settings in the analysis script."
  )
  diag_status <- if (inherits(x$diagnostics, "mfrm_diagnostics")) "available" else "not_available"
  add(
    "First-screen diagnostics",
    diag_status,
    if (identical(diag_status, "available")) "Diagnostics object and triage rows are available." else "Diagnostics were not available in this result object.",
    "summary(res)$triage; summary(res$diagnostics)",
    "Use for QC ordering and report-readiness checks.",
    "Diagnostics are evidence for follow-up and wording strength, not a standalone validity decision."
  )
  add(
    "Fit, separation, and precision",
    mfrm_report_component_status(x, "precision_review", available = "review"),
    "Precision review keeps fit size, standardized fit, separation, and reliability in separate lanes.",
    "summary(res$components$precision_review); precision_review_report(fit, diagnostics)",
    "Use for cautious wording about fit, precision, and separation.",
    "Do not collapse fit, separation, reliability, and ZSTD into one pass/fail rule."
  )
  add(
    "Category functioning",
    if (mfrm_report_has_category(x)) "available" else "not_requested",
    "Rating-scale/category tables or curves are available when collected by mfrm_results().",
    "rating_scale_table(fit, diagnostics); category_structure_report(fit)",
    "Use for score-scale interpretation and category-functioning prose.",
    "Category evidence supports score-scale review; it does not by itself establish validity."
  )
  add(
    "Bias screening",
    mfrm_report_component_status(x, "bias_screen", available = "review"),
    "Facet-level bias screening is available only when requested in mfrm_results().",
    "mfrm_results(fit, include = \"bias\"); estimate_bias(); bias_interaction_report()",
    "Use for screening language and targeted follow-up contrasts.",
    "Treat positive screens as prompts for substantive review, not final fairness conclusions."
  )
  add(
    "Misfit and pathway review",
    mfrm_report_component_status(x, "misfit_review", available = "review"),
    "Unexpected-response, displacement, and pathway-map surfaces can localize observations for review.",
    "mfrm_results(fit, include = \"misfit_review\"); plot(res, type = \"pathway\")",
    "Use for case-review notes and reviewer-facing diagnostic follow-up.",
    "Observation-level misfit is not an automatic exclusion or bias decision."
  )
  add(
    "Anchors and linking",
    mfrm_report_component_status(x, "linking_review", available = "review"),
    "Anchor-readiness is available from stored fit metadata when the linking preset is requested.",
    "mfrm_results(fit, include = \"linking\"); plot(res, type = \"anchors\")",
    "Use for operational scale-maintenance checks.",
    "Drift and equating require multiple fitted forms or waves; they are not inferred from one fit."
  )
  add(
    "Response-time QC",
    mfrm_report_component_status(x, "response_time_review", available = "available"),
    "Response-time summaries are available only when timing metadata are explicitly supplied to mfrm_results().",
    "mfrm_results(fit, include = \"response_time\", response_time = ..., response_time_data = ...); plot(res, type = \"response_time\")",
    "Use for descriptive timing context, rapid/slow-response screening, and QC appendices.",
    "Response-time review does not alter MFRM estimates, fit speed parameters, or define automatic exclusion rules."
  )
  add(
    "Network and connectivity",
    mfrm_report_component_status(x, "network_review", available = "available"),
    "Network review describes design connectivity and overlap structure.",
    "mfrm_results(fit, include = \"network\"); build_mfrm_network_review()",
    "Use for design and sparseness documentation.",
    "Connectivity evidence does not replace model fit, precision, or bias diagnostics."
  )
  add(
    "APA and manuscript wording",
    mfrm_report_component_status(x, "apa_outputs", available = "available"),
    "APA output assembly is available for supported RSM/PCM manuscript routes.",
    "mfrm_results(fit, include = \"publication\"); build_apa_outputs()",
    "Use as a draft wording template and table/caption scaffold.",
    "APA text must be edited against the actual study design, model choice, and validation argument."
  )
  add(
    "Tables, plots, and handoff",
    if (table_count > 0L || plot_count > 0L) "available" else "review",
    paste0(table_count, " table(s) and ", plot_count, " plot route(s) were indexed."),
    "build_summary_table_bundle(res); export_mfrm_results(res)",
    "Use for appendix, reviewer supplement, or reproducible handoff.",
    "Exported tables preserve evidence surfaces; they do not add new analyses."
  )
  if (identical(toupper(model), "GPCM")) {
    add(
      "GPCM scope",
      "caveat",
      "The fitted model is GPCM.",
      "gpcm_capability_matrix(); mfrmr_output_guide(\"gpcm\")",
      "Use direct outputs with documented helper-specific caveats.",
      "Do not present bounded GPCM helper coverage as universal equivalence with all RSM/PCM report routes."
    )
  }
  out <- do.call(rbind, rows)
  out$Focus <- mfrm_report_style_focus(style)
  out
}

mfrm_report_evidence_boundary <- function() {
  data.frame(
    EvidenceSource = c(
      "Model setup and convergence",
      "Fit, separation, and precision",
      "Category functioning",
      "Bias and DFF screening",
      "Misfit and pathway maps",
      "Anchors, linking, and drift",
      "Network and connectivity",
      "Response-time QC",
      "GPCM helper coverage",
      "APA-style wording"
    ),
    Use = c(
      "Document estimation settings, data roles, and run stability.",
      "Separate fit-size, standardized fit, separation, reliability, and uncertainty evidence.",
      "Describe how score categories function and where thresholds or curves need review.",
      "Identify candidate contrasts for follow-up fairness or interaction review.",
      "Localize unexpected observations for case review and substantive interpretation.",
      "Assess anchor readiness and operational scale-maintenance workflow.",
      "Describe design overlap, sparseness, and connectedness.",
      "Describe rapid/slow-response timing patterns as separate QC context.",
      "State which GPCM summaries are supported and which are caveated.",
      "Draft report prose, captions, and section maps."
    ),
    DoNotUseAs = c(
      "Proof that the construct interpretation is valid.",
      "A single pass/fail psychometric rule.",
      "Standalone validity evidence.",
      "Final fairness, bias, or invariance conclusion.",
      "Automatic exclusion rule for persons, raters, items, or observations.",
      "Single-fit evidence of drift or completed equating.",
      "A replacement for fit, precision, or bias diagnostics.",
      "A fitted speed parameter, speed-accuracy model, or automatic exclusion rule.",
      "A claim that every RSM/PCM report pathway has an equivalent GPCM route.",
      "Final manuscript text without study-specific editing."
    ),
    RecommendedRoute = c(
      "summary(res)$overview; specifications_report(fit)",
      "summary(res$components$precision_review); precision_review_report()",
      "rating_scale_table(); category_structure_report(); category_curves_report()",
      "mfrm_results(fit, include = \"bias\"); estimate_bias(); bias_interaction_report()",
      "mfrm_results(fit, include = \"misfit_review\"); unexpected_response_table(); plot(res, type = \"pathway\")",
      "mfrm_results(fit, include = \"linking\"); review_mfrm_anchors(); detect_anchor_drift(); build_equating_chain()",
      "mfrm_results(fit, include = \"network\"); build_mfrm_network_review()",
      "mfrm_results(fit, include = \"response_time\", response_time = ..., response_time_data = ...); response_time_review()",
      "gpcm_capability_matrix(); mfrmr_output_guide(\"gpcm\")",
      "mfrm_report(res, style = \"apa\"); build_apa_outputs()"
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_report_fit_criteria <- function(x) {
  band <- mfrm_misfit_thresholds()
  lower <- suppressWarnings(as.numeric(band["lower"]))
  upper <- suppressWarnings(as.numeric(band["upper"]))
  profiles <- fit_measure_threshold_profile_table(
    lower = lower,
    upper = upper,
    threshold_profiles = "all"
  )
  if (nrow(profiles) == 0L) {
    profiles <- data.frame(
      Profile = "active",
      ProfileLabel = "Active review band",
      Lower = lower,
      Upper = upper,
      Source = "Current package defaults",
      SuggestedUse = "Main fit-measure screen",
      stringsAsFactors = FALSE
    )
  }
  profiles$Metric <- "Infit/Outfit MnSq"
  profiles$ZSTDCut <- 2
  profiles$DecisionRole <- ifelse(
    profiles$Profile %in% "active",
    "main_report_screen",
    "sensitivity_or_context"
  )
  profiles$ReportBoundary <- paste0(
    "Use as a screening band. Report the selected profile and do not treat ",
    "a different published band as a contradiction unless it changes the substantive conclusion."
  )
  profiles$Route <- "fit_measures_table(fit, threshold_profiles = \"all\", fit_df_method = \"both\")"
  profiles[, c(
    "Profile", "ProfileLabel", "Metric", "Lower", "Upper", "ZSTDCut",
    "Source", "SuggestedUse", "DecisionRole", "ReportBoundary", "Route"
  ), drop = FALSE]
}

mfrm_report_zstd_conventions <- function() {
  data.frame(
    Convention = c(
      "engine df",
      "FACETS-style df",
      "Wilson-Hilferty ZSTD",
      "WHEXACT / linear approximation",
      "Report comparison route"
    ),
    FormulaOrRule = c(
      "Infit df = sum(Var * Weight); Outfit df = sum(Weight).",
      "Fourth-moment Wright-Masters-style df: 2 * numerator^2 / denominator; package columns DF_*_FACETS.",
      "(MnSq^(1/3) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df)).",
      "(MnSq - 1) * sqrt(df / 2) when whexact = TRUE.",
      "Keep engine and FACETS-style columns side by side with fit_df_method = \"both\"."
    ),
    PackageConstraint = c(
      "zstd_from_mnsq() returns NA when df < 1 to avoid unstable Wilson-Hilferty signs.",
      "zstd_from_mnsq_facets() allows positive df below 1 and caps reported ZSTD at +/-9.",
      "Requires finite positive MnSq and usable df; small df can dominate the transformation.",
      "Still requires usable df; use only when the analysis intentionally follows that convention.",
      "Compare MnSq first, then df, then ZSTD; classify flag changes as convention-sensitive."
    ),
    ReportingImplication = c(
      "Routine mfrmr diagnostics are conservative for very small df cells.",
      "FACETS-style review can reproduce the small-df behavior users expect from FACETS-like output, but it must be caveated.",
      "ZSTD is a standardization of MnSq, not a separate residual-fit statistic.",
      "State the transformation setting before interpreting ZSTD.",
      "Do not explain fit decisions from ZSTD alone when MnSq, df, or threshold profile differs."
    ),
    SourceBasis = c(
      "Package-native numerical guard for Wilson-Hilferty stability.",
      "FACETS/Winsteps fit-standardization documentation and Wright-Masters fourth-moment df convention.",
      "Wilson-Hilferty cube-root approximation used in Rasch fit standardization.",
      "Winsteps/FACETS WHEXACT documentation.",
      "facets_fit_df_guide(); facets_fit_review(); fit_measures_table()."
    ),
    Route = c(
      "diagnose_mfrm(fit_df_method = \"engine\")",
      "diagnose_mfrm(fit_df_method = \"facets\")",
      "fit_measures_table(..., fit_df_method = \"engine\" or \"facets\")",
      "diagnose_mfrm(..., whexact = TRUE)",
      "fit_measures_table(..., fit_df_method = \"both\")"
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_report_fit_decision_policy <- function() {
  data.frame(
    Step = seq_len(6L),
    Rule = c(
      "Choose and state the MnSq band",
      "Read MnSq before ZSTD",
      "Keep ZSTD convention visible",
      "Separate fit from precision",
      "Treat profile disagreement as sensitivity evidence",
      "Use local review before action"
    ),
    Rationale = c(
      "Published bands differ by setting; the active band and sensitivity profiles should both be visible.",
      "MnSq is the size of the fit signal; ZSTD is a df-dependent standardization of that signal.",
      "FACETS-style df can change |ZSTD| flags even when MnSq is unchanged.",
      "Fit, separation, reliability, and strata answer different questions.",
      "A row can be flagged under one defensible band and not another; report this as a review sensitivity.",
      "Element fit flags are prompts for inspecting responses, raters, items, categories, or design links."
    ),
    RecommendedRoute = c(
      "fit_measures_table(threshold_profiles = \"all\")",
      "fit_measures_table()$table[, c(\"Infit\", \"Outfit\", \"FitStatus\")]",
      "fit_measures_table(fit_df_method = \"both\")$df_sensitivity",
      "precision_review_report()$fit_separation_basis",
      "fit_measures_table()$profile_summary_by_facet",
      "unexpected_response_table(); displacement_table(); plot(res, type = \"pathway\")"
    ),
    ReportBoundary = c(
      "Do not silently mix bands across reports.",
      "Do not treat ZSTD as independent evidence from MnSq.",
      "Do not call a df-sensitive ZSTD change a substantive fit change without MnSq/context evidence.",
      "Do not use high reliability to excuse misfit or good fit to imply high precision.",
      "Do not present one threshold profile as universal.",
      "Do not remove levels or observations from fit flags alone."
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_report_fit_bundle <- function(x) {
  comps <- x$components %||% list()
  if ("fit_measures" %in% names(comps) &&
      inherits(comps$fit_measures, "mfrm_fit_measures")) {
    return(list(
      available = TRUE,
      bundle = comps$fit_measures,
      source = "res$components$fit_measures"
    ))
  }
  list(
    available = FALSE,
    bundle = NULL,
    source = "not_available"
  )
}

mfrm_report_column_or <- function(df, column, default) {
  df <- as.data.frame(df %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(df) == 0L || !(column %in% names(df))) return(default)
  value <- df[[column]][1]
  if (length(value) == 0L || is.null(value)) default else value
}

mfrm_report_has_facets_companion_fit <- function(fit_bundle) {
  tbl <- as.data.frame(fit_bundle$table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L) return(FALSE)
  companion <- intersect(
    c(
      "DF_Infit_FACETS", "DF_Outfit_FACETS",
      "InfitZSTD_FACETS", "OutfitZSTD_FACETS"
    ),
    names(tbl)
  )
  if (length(companion) == 0L) return(FALSE)
  any(vapply(tbl[companion], function(col) {
    any(is.finite(suppressWarnings(as.numeric(col))))
  }, logical(1)))
}

mfrm_report_fit_unavailable_row <- function(kind) {
  boundary <- paste0(
    "This report does not recompute diagnostics. Rebuild the result with ",
    "mfrm_results(fit, include = c(\"diagnostics\", \"tables\")) before ",
    "writing result-specific fit claims."
  )
  switch(
    kind,
    evidence = data.frame(
      Status = "not_available",
      Rows = NA_integer_,
      DisplayedRows = NA_integer_,
      UnderfitRows = NA_integer_,
      OverfitRows = NA_integer_,
      MixedRows = NA_integer_,
      WithinBandRows = NA_integer_,
      NotAvailableRows = NA_integer_,
      DfComparedRows = NA_integer_,
      DfSensitiveRows = NA_integer_,
      FlagChangedByDfRows = NA_integer_,
      LargeZSTDShiftRows = NA_integer_,
      DfConventionDifferenceRows = NA_integer_,
      FitDfMethod = NA_character_,
      ThresholdProfiles = NA_character_,
      FacetsCompanionAvailable = FALSE,
      Source = "not_available",
      Route = "mfrm_results(fit, include = c(\"diagnostics\", \"tables\"))",
      Boundary = boundary,
      stringsAsFactors = FALSE
    ),
    threshold = data.frame(
      Status = "not_available",
      Profile = NA_character_,
      ProfileLabel = NA_character_,
      Lower = NA_real_,
      Upper = NA_real_,
      Facet = NA_character_,
      Rows = NA_integer_,
      AvailableRows = NA_integer_,
      UnderfitRate = NA_real_,
      OverfitRate = NA_real_,
      MixedRate = NA_real_,
      AnyFlagRate = NA_real_,
      Source = "not_available",
      Route = "mfrm_results(fit, include = c(\"diagnostics\", \"tables\"))",
      ReportBoundary = boundary,
      stringsAsFactors = FALSE
    ),
    df_summary = data.frame(
      Status = "not_available",
      ComparedRows = NA_integer_,
      SameOrRoundingRows = NA_integer_,
      FlagChangedByDfRows = NA_integer_,
      LargeZSTDShiftRows = NA_integer_,
      DfConventionDifferenceRows = NA_integer_,
      DfSensitiveRows = NA_integer_,
      FitDfMethod = NA_character_,
      Source = "not_available",
      Route = "mfrm_results(fit, include = c(\"diagnostics\", \"tables\"))",
      ReportBoundary = boundary,
      stringsAsFactors = FALSE
    ),
    df_rows = data.frame(
      Status = "not_available",
      Facet = NA_character_,
      Level = NA_character_,
      DfSensitivityStatus = NA_character_,
      FlagChangedByDf = NA,
      MaxAbsZSTDDiff_FACETS_vs_ENGINE = NA_real_,
      MaxDFRelativeDifference_ENGINE_vs_FACETS = NA_real_,
      InfitZSTD_ENGINE = NA_real_,
      InfitZSTD_FACETS = NA_real_,
      OutfitZSTD_ENGINE = NA_real_,
      OutfitZSTD_FACETS = NA_real_,
      Interpretation = boundary,
      Source = "not_available",
      Route = "mfrm_results(fit, include = c(\"diagnostics\", \"tables\"))",
      ReportBoundary = boundary,
      stringsAsFactors = FALSE
    )
  )
}

mfrm_report_fit_evidence_summary <- function(x) {
  fit_info <- mfrm_report_fit_bundle(x)
  if (!isTRUE(fit_info$available)) {
    return(mfrm_report_fit_unavailable_row("evidence"))
  }
  fit_bundle <- fit_info$bundle
  sum_tbl <- as.data.frame(fit_bundle$summary %||% data.frame(), stringsAsFactors = FALSE)
  settings <- fit_bundle$settings %||% list()
  threshold_profiles <- as.character(settings$threshold_profiles %||% NA_character_)
  fit_df_method <- as.character(settings$fit_df_method %||% NA_character_)
  data.frame(
    Status = "available",
    Rows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "Rows", NA_integer_))),
    DisplayedRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "DisplayedRows", NA_integer_))),
    UnderfitRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "UnderfitRows", NA_integer_))),
    OverfitRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "OverfitRows", NA_integer_))),
    MixedRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "MixedRows", NA_integer_))),
    WithinBandRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "WithinBandRows", NA_integer_))),
    NotAvailableRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "NotAvailableRows", NA_integer_))),
    DfComparedRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "DfComparedRows", NA_integer_))),
    DfSensitiveRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "DfSensitiveRows", NA_integer_))),
    FlagChangedByDfRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "FlagChangedByDfRows", NA_integer_))),
    LargeZSTDShiftRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "LargeZSTDShiftRows", NA_integer_))),
    DfConventionDifferenceRows = suppressWarnings(as.integer(mfrm_report_column_or(sum_tbl, "DfConventionDifferenceRows", NA_integer_))),
    FitDfMethod = paste(fit_df_method[nzchar(fit_df_method)], collapse = ", "),
    ThresholdProfiles = paste(threshold_profiles[nzchar(threshold_profiles)], collapse = ", "),
    FacetsCompanionAvailable = mfrm_report_has_facets_companion_fit(fit_bundle),
    Source = fit_info$source,
    Route = "res$components$fit_measures$summary",
    Boundary = paste0(
      "Counts summarize the stored fit-measures component. Interpret MnSq ",
      "status, df-sensitive ZSTD shifts, separation, and reliability as ",
      "separate evidence streams."
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_report_fit_threshold_sensitivity <- function(x) {
  fit_info <- mfrm_report_fit_bundle(x)
  if (!isTRUE(fit_info$available)) {
    return(mfrm_report_fit_unavailable_row("threshold"))
  }
  fit_bundle <- fit_info$bundle
  profile <- as.data.frame(
    fit_bundle$profile_summary_overall %||% data.frame(),
    stringsAsFactors = FALSE
  )
  if (nrow(profile) == 0L) {
    profile <- as.data.frame(fit_bundle$profile_summary %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(profile) > 0L && "Facet" %in% names(profile)) {
      profile <- profile[as.character(profile$Facet) %in% "All facets", , drop = FALSE]
    }
  }
  if (nrow(profile) == 0L) {
    out <- mfrm_report_fit_unavailable_row("threshold")
    out$Status <- "not_available"
    out$Source <- fit_info$source
    out$ReportBoundary <- "No threshold-profile summary was stored in the fit-measures component."
    return(out)
  }
  needed <- c(
    "Profile", "ProfileLabel", "Lower", "Upper", "Facet", "Rows",
    "AvailableRows", "UnderfitRate", "OverfitRate", "MixedRate", "AnyFlagRate"
  )
  for (nm in setdiff(needed, names(profile))) {
    profile[[nm]] <- if (nm %in% c("Lower", "Upper", "UnderfitRate", "OverfitRate", "MixedRate", "AnyFlagRate")) {
      NA_real_
    } else if (nm %in% c("Rows", "AvailableRows")) {
      NA_integer_
    } else {
      NA_character_
    }
  }
  profile <- profile[, needed, drop = FALSE]
  profile$Status <- "available"
  profile$Source <- fit_info$source
  profile$Route <- "res$components$fit_measures$profile_summary_overall"
  profile$ReportBoundary <- paste0(
    "Use profile disagreement as sensitivity evidence. Do not present one ",
    "published MnSq band as universal."
  )
  profile[, c("Status", needed, "Source", "Route", "ReportBoundary"), drop = FALSE]
}

mfrm_report_fit_df_sensitivity_summary <- function(x) {
  fit_info <- mfrm_report_fit_bundle(x)
  if (!isTRUE(fit_info$available)) {
    return(mfrm_report_fit_unavailable_row("df_summary"))
  }
  fit_bundle <- fit_info$bundle
  df_sum <- as.data.frame(fit_bundle$df_sensitivity_summary %||% data.frame(), stringsAsFactors = FALSE)
  df_sensitive <- as.data.frame(fit_bundle$df_sensitive %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(df_sum) == 0L) {
    out <- mfrm_report_fit_unavailable_row("df_summary")
    out$Status <- "not_available"
    out$Source <- fit_info$source
    out$ReportBoundary <- "No df-sensitivity summary was stored in the fit-measures component."
    return(out)
  }
  needed <- c(
    "ComparedRows", "SameOrRoundingRows", "FlagChangedByDfRows",
    "LargeZSTDShiftRows", "DfConventionDifferenceRows"
  )
  for (nm in setdiff(needed, names(df_sum))) df_sum[[nm]] <- NA_integer_
  settings <- fit_bundle$settings %||% list()
  out <- df_sum[, needed, drop = FALSE]
  out$DfSensitiveRows <- nrow(df_sensitive)
  out$FitDfMethod <- as.character(settings$fit_df_method %||% NA_character_)
  out$Status <- if (nrow(df_sensitive) > 0L) "review" else "available"
  out$Source <- fit_info$source
  out$Route <- "res$components$fit_measures$df_sensitivity_summary"
  out$ReportBoundary <- paste0(
    "A df-sensitive row means the ZSTD interpretation changed or moved ",
    "materially under engine-vs-FACETS-style standardization; it is not a ",
    "different MnSq fit statistic."
  )
  out[, c(
    "Status", needed, "DfSensitiveRows", "FitDfMethod",
    "Source", "Route", "ReportBoundary"
  ), drop = FALSE]
}

mfrm_report_fit_df_sensitive_rows <- function(x, top_n = 10L) {
  fit_info <- mfrm_report_fit_bundle(x)
  if (!isTRUE(fit_info$available)) {
    return(mfrm_report_fit_unavailable_row("df_rows"))
  }
  fit_bundle <- fit_info$bundle
  rows <- as.data.frame(fit_bundle$df_sensitive %||% data.frame(), stringsAsFactors = FALSE)
  needed <- c(
    "Facet", "Level", "DfSensitivityStatus", "FlagChangedByDf",
    "MaxAbsZSTDDiff_FACETS_vs_ENGINE",
    "MaxDFRelativeDifference_ENGINE_vs_FACETS",
    "InfitZSTD_ENGINE", "InfitZSTD_FACETS",
    "OutfitZSTD_ENGINE", "OutfitZSTD_FACETS",
    "Interpretation"
  )
  if (nrow(rows) == 0L) {
    out <- mfrm_report_fit_unavailable_row("df_rows")
    out$Status <- "none"
    out$Source <- fit_info$source
    out$Interpretation <- "No df-sensitive rows were stored in the fit-measures component."
    out$ReportBoundary <- "No row-level df/ZSTD sensitivity prompt is available for follow-up."
    return(out)
  }
  for (nm in setdiff(needed, names(rows))) rows[[nm]] <- NA
  rows <- rows[, needed, drop = FALSE]
  z_shift <- suppressWarnings(as.numeric(rows$MaxAbsZSTDDiff_FACETS_vs_ENGINE))
  flag <- rows$FlagChangedByDf %in% TRUE
  ord <- order(-as.integer(flag), -ifelse(is.finite(z_shift), z_shift, -Inf), rows$Facet, rows$Level)
  rows <- rows[ord, , drop = FALSE]
  top_n <- suppressWarnings(as.integer(top_n[1]))
  if (is.finite(top_n) && top_n > 0L) {
    rows <- utils::head(rows, top_n)
  }
  rows$Status <- "review"
  rows$Source <- fit_info$source
  rows$Route <- "res$components$fit_measures$df_sensitive"
  rows$ReportBoundary <- paste0(
    "Inspect the row context before using a ZSTD-only flag in report text; ",
    "MnSq size and substantive role remain primary."
  )
  rows[, c("Status", needed, "Source", "Route", "ReportBoundary"), drop = FALSE]
}

mfrm_report_fmt_int <- function(x) {
  x <- suppressWarnings(as.integer(x[1]))
  if (!is.finite(x)) return("NA")
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

mfrm_report_fmt_pct <- function(x) {
  x <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x)) return("NA")
  paste0(format(round(100 * x, 1), nsmall = 1, trim = TRUE), "%")
}

mfrm_report_fit_template_audience <- function(style) {
  switch(
    style,
    apa = "APA manuscript",
    validation = "Validity argument",
    reviewer = "Reviewer response",
    technical = "Technical appendix",
    qc = "QC report",
    "Report"
  )
}

mfrm_report_fit_threshold_sentence <- function(threshold_tbl) {
  threshold_tbl <- as.data.frame(threshold_tbl %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(threshold_tbl) == 0L || !"AnyFlagRate" %in% names(threshold_tbl)) {
    return("Threshold-profile sensitivity was not available in the stored result.")
  }
  rates <- suppressWarnings(as.numeric(threshold_tbl$AnyFlagRate))
  rates <- rates[is.finite(rates)]
  active <- threshold_tbl[as.character(threshold_tbl$Profile %||% "") %in% "active", , drop = FALSE]
  active_rate <- if (nrow(active) > 0L) {
    suppressWarnings(as.numeric(active$AnyFlagRate[1]))
  } else {
    NA_real_
  }
  if (length(rates) == 0L) {
    return("Threshold-profile sensitivity was stored, but no finite flag rates were available.")
  }
  paste0(
    "Across the stored mean-square threshold profiles, any-flag rates ranged from ",
    mfrm_report_fmt_pct(min(rates, na.rm = TRUE)), " to ",
    mfrm_report_fmt_pct(max(rates, na.rm = TRUE)), ". The active profile rate was ",
    mfrm_report_fmt_pct(active_rate), "."
  )
}

mfrm_report_template_evidence_table <- function(evidence_used) {
  evidence_used <- as.character(evidence_used %||% "")
  out <- vapply(evidence_used, function(x) {
    x <- trimws(strsplit(x, ";", fixed = TRUE)[[1]][1] %||% "")
    x <- sub("^report\\$", "", x)
    x <- sub("\\$.*$", "", x)
    x <- sub("\\(.*$", "", x)
    x <- trimws(x)
    if (!nzchar(x)) "evidence_not_recorded" else x
  }, character(1))
  unname(out)
}

mfrm_report_template_claim_strength <- function(topic, default = "write_with_caveat") {
  topic <- tolower(as.character(topic %||% ""))
  out <- rep(default, length(topic))
  out[grepl("unavailable|not requested", topic)] <- "not_supported_without_followup"
  out[grepl("boundary", topic)] <- "descriptive_only"
  out[grepl("dff follow-up|fairness|interaction|drift|equating|gpcm", topic)] <- "not_supported_without_followup"
  out[grepl("fit-status|precision-tier|separation|reliability|strata|facet-level|anchor-readiness", topic)] <- "descriptive_only"
  out[grepl("threshold|zstd|df/zstd|screen-positive|unexpected|displacement|pathway|linking-risk", topic)] <- "write_with_caveat"
  unname(out)
}

mfrm_report_template_recommended_use <- function(topic) {
  topic <- tolower(as.character(topic %||% ""))
  out <- rep("report_sentence_scaffold", length(topic))
  out[grepl("unavailable|not requested", topic)] <- "request_evidence_before_writing"
  out[grepl("boundary", topic)] <- "reporting_guardrail"
  out[grepl("threshold|zstd|df/zstd|drift|equating|gpcm", topic)] <- "methods_or_appendix_caveat"
  out[grepl("dff|fairness|interaction", topic)] <- "targeted_followup_before_claim"
  unname(out)
}

mfrm_report_template_enrich <- function(df,
                                        boundary_type,
                                        default_claim_strength = "write_with_caveat") {
  df <- as.data.frame(df %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(df) == 0L) return(df)
  evidence_used <- as.character(df$EvidenceUsed %||% "")
  topic <- as.character(df$Topic %||% "")
  if (!"EvidenceTable" %in% names(df)) {
    df$EvidenceTable <- mfrm_report_template_evidence_table(evidence_used)
  }
  if (!"EvidenceRoute" %in% names(df)) {
    report_table <- grepl(
      paste(c(
        "_evidence_summary",
        "fit_threshold_sensitivity",
        "fit_df_sensitivity_summary",
        "zstd_conventions",
        "fit_decision_policy",
        "precision_basis",
        "evidence_boundary"
      ), collapse = "|"),
      df$EvidenceTable
    )
    df$EvidenceRoute <- ifelse(
      grepl("^report\\$", evidence_used),
      evidence_used,
      ifelse(report_table, paste0("report$", df$EvidenceTable), as.character(df$Route %||% ""))
    )
  }
  if (!"BoundaryType" %in% names(df)) {
    df$BoundaryType <- rep(as.character(boundary_type), nrow(df))
  }
  if (!"ClaimStrength" %in% names(df)) {
    df$ClaimStrength <- mfrm_report_template_claim_strength(
      topic,
      default = default_claim_strength
    )
  }
  if (!"RecommendedUse" %in% names(df)) {
    df$RecommendedUse <- mfrm_report_template_recommended_use(topic)
  }
  df
}

mfrm_report_fit_reporting_templates <- function(style,
                                                fit_evidence_summary,
                                                fit_threshold_sensitivity,
                                                fit_df_sensitivity_summary) {
  ev <- as.data.frame(fit_evidence_summary %||% data.frame(), stringsAsFactors = FALSE)
  th <- as.data.frame(fit_threshold_sensitivity %||% data.frame(), stringsAsFactors = FALSE)
  df <- as.data.frame(fit_df_sensitivity_summary %||% data.frame(), stringsAsFactors = FALSE)
  available <- nrow(ev) > 0L && as.character(ev$Status[1] %||% "") %in% "available"
  audience <- mfrm_report_fit_template_audience(style)
  style_lead <- switch(
    style,
    apa = "Report in the manuscript as a descriptive diagnostic, then cite the exact table route.",
    validation = "Use as fit evidence inside a validity argument, not as a standalone validity proof.",
    reviewer = "Use as a direct response map: state what was checked, what changed, and what remains caveated.",
    technical = "Use as appendix wording tied to reproducible table routes.",
    qc = "Use as first-screen QC wording before manuscript or reviewer handoff.",
    "Use as report wording tied to the stored result object."
  )

  if (!available) {
    return(mfrm_report_template_enrich(data.frame(
      Audience = audience,
      Topic = "Fit evidence unavailable",
      Template = paste0(
        "Result-specific fit wording was not generated because the stored ",
        "mfrm_results object does not contain a fit-measures component."
      ),
      EvidenceUsed = "fit_evidence_summary",
      Caveat = "Rebuild the result with mfrm_results(fit, include = c(\"diagnostics\", \"tables\")) before reporting fit counts.",
      Route = "mfrm_results(fit, include = c(\"diagnostics\", \"tables\"))",
      Style = style,
      stringsAsFactors = FALSE
    ), boundary_type = "fit_not_validity", default_claim_strength = "not_supported_without_followup"))
  }

  rows <- list()
  add <- function(topic, template, evidence, caveat, route) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Audience = audience,
      Topic = as.character(topic),
      Template = as.character(template),
      EvidenceUsed = as.character(evidence),
      Caveat = as.character(caveat),
      Route = as.character(route),
      Style = style,
      stringsAsFactors = FALSE
    )
  }

  rows_n <- mfrm_report_fmt_int(ev$Rows[1])
  underfit <- mfrm_report_fmt_int(ev$UnderfitRows[1])
  overfit <- mfrm_report_fmt_int(ev$OverfitRows[1])
  mixed <- mfrm_report_fmt_int(ev$MixedRows[1])
  within <- mfrm_report_fmt_int(ev$WithinBandRows[1])
  fit_df_method <- as.character(ev$FitDfMethod[1] %||% "")
  threshold_profiles <- as.character(ev$ThresholdProfiles[1] %||% "")
  facets_companion <- isTRUE(ev$FacetsCompanionAvailable[1])

  add(
    "Fit-status wording",
    paste0(
      style_lead, " Element fit was screened for ", rows_n,
      " facet-element row(s) using the stored mean-square profile set",
      if (nzchar(threshold_profiles)) paste0(" (", threshold_profiles, ")") else "",
      "; ", underfit, " row(s) were flagged for underfit, ", overfit,
      " for overfit, ", mixed, " as mixed, and ", within,
      " remained within the selected band."
    ),
    "fit_evidence_summary",
    "This sentence reports a screening table, not a global model-validity decision.",
    "report$fit_evidence_summary"
  )

  add(
    "Threshold-profile wording",
    mfrm_report_fit_threshold_sentence(th),
    "fit_threshold_sensitivity",
    "Use profile disagreement as sensitivity evidence; do not silently mix fit bands across reports.",
    "report$fit_threshold_sensitivity"
  )

  add(
    "ZSTD-convention wording",
    paste0(
      "ZSTD values were treated as df-dependent standardizations of the same MnSq values. ",
      if (facets_companion) {
        paste0(
          "Engine and FACETS-style companion df/ZSTD columns were available",
          if (nzchar(fit_df_method)) paste0(" under fit_df_method = \"", fit_df_method, "\"") else "",
          "."
        )
      } else {
        "FACETS-style companion df/ZSTD columns were not available in this stored result."
      }
    ),
    "zstd_conventions; fit_evidence_summary",
    "Read MnSq size first; use ZSTD to explain standardization, not as independent residual evidence.",
    "report$zstd_conventions"
  )

  compared <- mfrm_report_fmt_int(df$ComparedRows[1])
  sensitive <- mfrm_report_fmt_int(df$DfSensitiveRows[1])
  changed <- mfrm_report_fmt_int(df$FlagChangedByDfRows[1])
  large <- mfrm_report_fmt_int(df$LargeZSTDShiftRows[1])
  convention <- mfrm_report_fmt_int(df$DfConventionDifferenceRows[1])
  add(
    "DF/ZSTD sensitivity wording",
    paste0(
      "Engine-vs-FACETS-style df comparison covered ", compared,
      " row(s): ", sensitive, " row(s) were df-sensitive, ", changed,
      " changed the |ZSTD| flag status, ", large,
      " had a large ZSTD shift without necessarily changing flag status, and ",
      convention, " showed a df-convention difference."
    ),
    "fit_df_sensitivity_summary",
    "A df-sensitive ZSTD result is a convention-sensitive review prompt, not a different MnSq fit signal.",
    "report$fit_df_sensitivity_summary; report$fit_df_sensitive_rows"
  )

  add(
    "Boundary wording",
    paste0(
      "Report fit, ZSTD standardization, separation/reliability, and local ",
      "case review in separate sentences. Avoid wording such as 'the model ",
      "passed fit' unless the stated threshold profile, df convention, and ",
      "follow-up review all support that narrower claim."
    ),
    "fit_decision_policy",
    "This boundary is intentionally conservative because published MnSq bands and ZSTD conventions differ.",
    "report$fit_decision_policy"
  )

  mfrm_report_template_enrich(
    do.call(rbind, rows),
    boundary_type = "fit_not_validity",
    default_claim_strength = "write_with_caveat"
  )
}

mfrm_report_precision_bundle <- function(x) {
  comps <- x$components %||% list()
  if ("precision_review" %in% names(comps) &&
      inherits(comps$precision_review, "mfrm_precision_review")) {
    return(list(
      available = TRUE,
      bundle = comps$precision_review,
      source = "res$components$precision_review"
    ))
  }
  list(
    available = FALSE,
    bundle = NULL,
    source = "not_available"
  )
}

mfrm_report_finite_range <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(c(Min = NA_real_, Max = NA_real_))
  c(Min = min(x), Max = max(x))
}

mfrm_report_fmt_num <- function(x, digits = 2L) {
  x <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

mfrm_report_precision_evidence_summary <- function(x) {
  precision_info <- mfrm_report_precision_bundle(x)
  reliability_tbl <- as.data.frame(
    x$diagnostics$reliability %||% data.frame(),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(precision_info$available) && nrow(reliability_tbl) == 0L) {
    return(data.frame(
      Status = "not_available",
      PrecisionTier = NA_character_,
      SupportsFormalInference = NA,
      ReliabilityRows = NA_integer_,
      ReviewOrWarnChecks = NA_integer_,
      MinSeparation = NA_real_,
      MaxSeparation = NA_real_,
      MinReliability = NA_real_,
      MaxReliability = NA_real_,
      MinStrata = NA_real_,
      MaxStrata = NA_real_,
      ZeroSeparationRows = NA_integer_,
      ZeroReliabilityRows = NA_integer_,
      ReliabilityUse = NA_character_,
      Source = "not_available",
      Route = "mfrm_results(fit, include = \"precision\")",
      Boundary = paste0(
        "Precision/separation wording requires the stored precision review ",
        "or diagnostics$reliability. Rebuild the result with precision included."
      ),
      stringsAsFactors = FALSE
    ))
  }

  precision <- precision_info$bundle %||% list()
  profile <- as.data.frame(precision$profile %||% data.frame(), stringsAsFactors = FALSE)
  checks <- as.data.frame(precision$checks %||% data.frame(), stringsAsFactors = FALSE)
  sep_range <- mfrm_report_finite_range(reliability_tbl$Separation %||% numeric(0))
  rel_range <- mfrm_report_finite_range(reliability_tbl$Reliability %||% numeric(0))
  strata_range <- mfrm_report_finite_range(reliability_tbl$Strata %||% numeric(0))
  sep <- suppressWarnings(as.numeric(reliability_tbl$Separation %||% numeric(0)))
  rel <- suppressWarnings(as.numeric(reliability_tbl$Reliability %||% numeric(0)))
  reliability_use <- if ("ReliabilityUse" %in% names(reliability_tbl)) {
    paste(sort(unique(as.character(reliability_tbl$ReliabilityUse))), collapse = ", ")
  } else {
    NA_character_
  }
  review_checks <- if (nrow(checks) > 0L && "Status" %in% names(checks)) {
    sum(as.character(checks$Status) %in% c("review", "warn"), na.rm = TRUE)
  } else {
    NA_integer_
  }
  source <- if (isTRUE(precision_info$available)) {
    precision_info$source
  } else {
    "res$diagnostics$reliability"
  }
  data.frame(
    Status = if (isTRUE(precision_info$available)) "available" else "available_without_precision_review",
    PrecisionTier = as.character(profile$PrecisionTier[1] %||% NA_character_),
    SupportsFormalInference = isTRUE(profile$SupportsFormalInference[1] %||% FALSE),
    ReliabilityRows = nrow(reliability_tbl),
    ReviewOrWarnChecks = review_checks,
    MinSeparation = sep_range[["Min"]],
    MaxSeparation = sep_range[["Max"]],
    MinReliability = rel_range[["Min"]],
    MaxReliability = rel_range[["Max"]],
    MinStrata = strata_range[["Min"]],
    MaxStrata = strata_range[["Max"]],
    ZeroSeparationRows = sum(is.finite(sep) & sep <= 0, na.rm = TRUE),
    ZeroReliabilityRows = sum(is.finite(rel) & rel <= 0, na.rm = TRUE),
    ReliabilityUse = reliability_use,
    Source = source,
    Route = "report$precision_evidence_summary; res$components$precision_review; res$diagnostics$reliability",
    Boundary = paste0(
      "Separation, reliability, and strata summarize spread relative to ",
      "measurement error. They are not inter-rater agreement, model fit, or ",
      "standalone validity evidence."
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_report_precision_basis <- function(x) {
  precision_info <- mfrm_report_precision_bundle(x)
  if (!isTRUE(precision_info$available)) {
    return(data.frame(
      Topic = "Separation reliability and strata",
      SourceBasis = "Wright & Masters G/R/H convention",
      PackageSurface = "diagnostics$reliability; precision_review_report()",
      Interpretation = "Precision review was not stored in this mfrm_results object.",
      ValidationUse = "Rebuild with include = \"precision\" before using source-grounded precision wording.",
      Availability = "not_requested",
      Source = "not_available",
      stringsAsFactors = FALSE
    ))
  }
  basis <- as.data.frame(
    precision_info$bundle$fit_separation_basis %||% data.frame(),
    stringsAsFactors = FALSE
  )
  if (nrow(basis) == 0L) {
    return(data.frame(
      Topic = "Separation reliability and strata",
      SourceBasis = NA_character_,
      PackageSurface = "precision_review_report()$fit_separation_basis",
      Interpretation = "No precision-basis rows were stored.",
      ValidationUse = "Avoid source-grounded precision wording until the basis table is available.",
      Availability = "not_available",
      Source = precision_info$source,
      stringsAsFactors = FALSE
    ))
  }
  keep <- grepl("Separation|reliability|strata|QC thresholds", basis$Topic, ignore.case = TRUE)
  basis <- basis[keep, , drop = FALSE]
  if (nrow(basis) == 0L) {
    basis <- as.data.frame(
      precision_info$bundle$fit_separation_basis %||% data.frame(),
      stringsAsFactors = FALSE
    )
  }
  basis$Source <- precision_info$source
  basis
}

mfrm_report_precision_reporting_templates <- function(style,
                                                      precision_evidence_summary,
                                                      precision_basis) {
  ev <- as.data.frame(precision_evidence_summary %||% data.frame(), stringsAsFactors = FALSE)
  basis <- as.data.frame(precision_basis %||% data.frame(), stringsAsFactors = FALSE)
  audience <- mfrm_report_fit_template_audience(style)
  if (nrow(ev) == 0L || as.character(ev$Status[1] %||% "") %in% "not_available") {
    return(mfrm_report_template_enrich(data.frame(
      Audience = audience,
      Topic = "Precision evidence unavailable",
      Template = paste0(
        "Precision, separation, reliability, and strata wording was not ",
        "generated because no stored precision evidence was available."
      ),
      EvidenceUsed = "precision_evidence_summary",
      Caveat = "Rebuild the result with mfrm_results(fit, include = \"precision\") before writing precision claims.",
      Route = "mfrm_results(fit, include = \"precision\")",
      Style = style,
      stringsAsFactors = FALSE
    ), boundary_type = "precision_not_agreement", default_claim_strength = "not_supported_without_followup"))
  }

  rows <- list()
  add <- function(topic, template, evidence, caveat, route) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Audience = audience,
      Topic = as.character(topic),
      Template = as.character(template),
      EvidenceUsed = as.character(evidence),
      Caveat = as.character(caveat),
      Route = as.character(route),
      Style = style,
      stringsAsFactors = FALSE
    )
  }

  tier <- as.character(ev$PrecisionTier[1] %||% "")
  supports <- isTRUE(ev$SupportsFormalInference[1])
  rows_n <- mfrm_report_fmt_int(ev$ReliabilityRows[1])
  checks <- mfrm_report_fmt_int(ev$ReviewOrWarnChecks[1])
  reliability_use <- as.character(ev$ReliabilityUse[1] %||% "")
  style_lead <- switch(
    style,
    apa = "Report as precision evidence, not as inter-rater agreement.",
    validation = "Use as precision evidence inside the validity argument, with explicit limits.",
    reviewer = "Use as a response-ready explanation of what the precision indices do and do not support.",
    technical = "Use as appendix wording tied to diagnostics$reliability and precision_review_report().",
    qc = "Use as first-screen precision wording before stronger report claims.",
    "Use as precision wording tied to the stored result object."
  )

  add(
    "Precision-tier wording",
    paste0(
      style_lead, " The precision review classified the run as ",
      if (nzchar(tier)) tier else "unspecified",
      "; formal inference support was ", if (supports) "available" else "not supported",
      ", and ", checks, " precision check(s) were marked review/warn."
    ),
    "precision_evidence_summary; precision_review_report()$profile",
    "A favorable precision tier does not override misfit, convergence, linking, or design problems.",
    "report$precision_evidence_summary"
  )

  add(
    "Separation wording",
    paste0(
      "Facet separation was available for ", rows_n,
      " facet row(s), ranging from ",
      mfrm_report_fmt_num(ev$MinSeparation[1]), " to ",
      mfrm_report_fmt_num(ev$MaxSeparation[1]),
      ". Interpret separation as spread relative to average measurement error."
    ),
    "diagnostics$reliability$Separation",
    "Do not describe separation as observed rater agreement or as proof of construct validity.",
    "res$diagnostics$reliability"
  )

  add(
    "Reliability wording",
    paste0(
      "Facet separation reliability ranged from ",
      mfrm_report_fmt_num(ev$MinReliability[1]), " to ",
      mfrm_report_fmt_num(ev$MaxReliability[1]),
      if (nzchar(reliability_use)) paste0("; reliability-use labels were: ", reliability_use, ".") else "."
    ),
    "diagnostics$reliability$Reliability",
    "This is Rasch/FACETS-style separation reliability, not classical inter-rater agreement.",
    "res$diagnostics$reliability"
  )

  add(
    "Strata wording",
    paste0(
      "Facet strata ranged from ",
      mfrm_report_fmt_num(ev$MinStrata[1]), " to ",
      mfrm_report_fmt_num(ev$MaxStrata[1]),
      " under the Wright/Masters G/R/H convention."
    ),
    "diagnostics$reliability$Strata; precision_basis",
    "Use strata as a precision-spread summary; do not turn it into an independent quality gate.",
    "report$precision_basis"
  )

  source_basis <- if (nrow(basis) > 0L && "SourceBasis" %in% names(basis)) {
    paste(utils::head(unique(as.character(basis$SourceBasis)), 2L), collapse = " | ")
  } else {
    "Wright/Masters G/R/H convention"
  }
  add(
    "Boundary wording",
    paste0(
      "State the precision tier and source convention before interpreting ",
      "separation, reliability, or strata. Source basis: ", source_basis, "."
    ),
    "precision_basis",
    "Do not use high reliability to excuse misfit, and do not use good fit to imply high precision.",
    "report$precision_basis; report$fit_decision_policy"
  )

  mfrm_report_template_enrich(
    do.call(rbind, rows),
    boundary_type = "precision_not_agreement",
    default_claim_strength = "write_with_caveat"
  )
}

mfrm_report_bias_bundle <- function(x) {
  comps <- x$components %||% list()
  if ("bias_screen" %in% names(comps) &&
      inherits(comps$bias_screen, "mfrm_bias_screen")) {
    return(list(
      available = TRUE,
      bundle = comps$bias_screen,
      source = "res$components$bias_screen"
    ))
  }
  list(
    available = FALSE,
    bundle = NULL,
    source = "not_requested"
  )
}

mfrm_report_bias_evidence_summary <- function(x) {
  bias_info <- mfrm_report_bias_bundle(x)
  if (!isTRUE(bias_info$available)) {
    return(data.frame(
      Status = "not_requested",
      Rows = NA_integer_,
      Facets = NA_integer_,
      NonPersonFacets = NA_integer_,
      MaxAbsBias = NA_real_,
      MaxAbsStdResidual = NA_real_,
      ResidualTScreenPositiveRows = NA_integer_,
      ChiSqScreenPositiveRows = NA_integer_,
      ExplicitInteractionSelected = FALSE,
      InteractionStatus = "not_requested",
      Source = "not_requested",
      Route = "mfrm_results(fit, include = \"bias\")",
      Boundary = paste0(
        "Bias/DFF wording requires the bias preset or an explicit bias/DFF ",
        "helper call. Do not infer fairness conclusions from omitted sections."
      ),
      stringsAsFactors = FALSE
    ))
  }

  bias_bundle <- bias_info$bundle
  tbl <- as.data.frame(bias_bundle$table %||% data.frame(), stringsAsFactors = FALSE)
  guidance <- as.data.frame(bias_bundle$guidance %||% data.frame(), stringsAsFactors = FALSE)
  available_facets <- as.data.frame(
    bias_bundle$available_facets %||% data.frame(),
    stringsAsFactors = FALSE
  )
  bias <- suppressWarnings(as.numeric(tbl$Bias %||% numeric(0)))
  std_resid <- suppressWarnings(as.numeric(tbl$MeanStdResidual %||% numeric(0)))
  t_resid <- suppressWarnings(as.numeric(tbl$t_Residual %||% numeric(0)))
  p_resid <- suppressWarnings(as.numeric(tbl$p_Residual %||% numeric(0)))
  chi_p <- suppressWarnings(as.numeric(tbl$ChiP %||% numeric(0)))
  residual_screen <- (is.finite(t_resid) & abs(t_resid) >= 2) |
    (is.finite(p_resid) & p_resid <= 0.05)
  chi_screen <- is.finite(chi_p) & chi_p <= 0.05
  interaction_status <- if (nrow(guidance) > 0L &&
      all(c("Area", "Status") %in% names(guidance))) {
    as.character(guidance$Status[guidance$Area %in% "Interaction bias screen"][1] %||% "not_available")
  } else {
    "not_available"
  }
  facets <- if ("Facet" %in% names(tbl)) unique(as.character(tbl$Facet)) else character(0)
  data.frame(
    Status = "available",
    Rows = nrow(tbl),
    Facets = length(facets),
    NonPersonFacets = nrow(available_facets),
    MaxAbsBias = if (any(is.finite(bias))) max(abs(bias), na.rm = TRUE) else NA_real_,
    MaxAbsStdResidual = if (any(is.finite(std_resid))) max(abs(std_resid), na.rm = TRUE) else NA_real_,
    ResidualTScreenPositiveRows = sum(residual_screen, na.rm = TRUE),
    ChiSqScreenPositiveRows = sum(chi_screen, na.rm = TRUE),
    ExplicitInteractionSelected = FALSE,
    InteractionStatus = interaction_status,
    Source = bias_info$source,
    Route = "res$components$bias_screen; estimate_bias(fit, diagnostics, facet_a = ..., facet_b = ...)",
    Boundary = paste0(
      "Facet-level bias rows are screening prompts. Interaction bias, DFF, ",
      "and fairness claims require explicit facet or group contrasts and ",
      "substantive review."
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_report_bias_reporting_templates <- function(style, bias_evidence_summary) {
  ev <- as.data.frame(bias_evidence_summary %||% data.frame(), stringsAsFactors = FALSE)
  audience <- mfrm_report_fit_template_audience(style)
  if (nrow(ev) == 0L || !as.character(ev$Status[1] %||% "") %in% "available") {
    return(mfrm_report_template_enrich(data.frame(
      Audience = audience,
      Topic = "Bias/DFF evidence not requested",
      Template = paste0(
        "Bias, DFF, or fairness-screen wording was not generated because ",
        "the mfrm_results object was not built with include = \"bias\"."
      ),
      EvidenceUsed = "bias_evidence_summary",
      Caveat = "Request the bias preset or run an explicit bias/DFF helper before writing fairness language.",
      Route = "mfrm_results(fit, include = \"bias\"); estimate_bias(); analyze_dff()",
      Style = style,
      stringsAsFactors = FALSE
    ), boundary_type = "screen_not_fairness_decision", default_claim_strength = "not_supported_without_followup"))
  }

  rows <- list()
  add <- function(topic, template, evidence, caveat, route) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Audience = audience,
      Topic = as.character(topic),
      Template = as.character(template),
      EvidenceUsed = as.character(evidence),
      Caveat = as.character(caveat),
      Route = as.character(route),
      Style = style,
      stringsAsFactors = FALSE
    )
  }
  style_lead <- switch(
    style,
    apa = "Report as a screening diagnostic only when bias or DFF was a study question.",
    validation = "Use as a fairness-related screening layer inside the validity argument.",
    reviewer = "Use as reviewer-facing evidence about which bias screens were requested and what remains untested.",
    technical = "Use as appendix wording tied to bias-screen table routes.",
    qc = "Use as first-screen fairness/bias triage wording.",
    "Use as bias-screen wording tied to the stored result object."
  )
  add(
    "Facet-level bias-screen wording",
    paste0(
      style_lead, " The stored facet-level bias screen contained ",
      mfrm_report_fmt_int(ev$Rows[1]), " row(s) across ",
      mfrm_report_fmt_int(ev$Facets[1]), " facet(s), with maximum absolute ",
      "observed-minus-expected bias of ", mfrm_report_fmt_num(ev$MaxAbsBias[1]),
      " score units and maximum absolute mean standardized residual of ",
      mfrm_report_fmt_num(ev$MaxAbsStdResidual[1]), "."
    ),
    "bias_evidence_summary",
    "Facet-level screens are prompts for follow-up, not final fairness or invariance conclusions.",
    "report$bias_evidence_summary; res$components$bias_screen$table"
  )
  add(
    "Screen-positive wording",
    paste0(
      "Using the package's first-screen residual and chi-square summaries, ",
      mfrm_report_fmt_int(ev$ResidualTScreenPositiveRows[1]),
      " row(s) were residual-screen-positive and ",
      mfrm_report_fmt_int(ev$ChiSqScreenPositiveRows[1]),
      " row(s) were chi-square-screen-positive."
    ),
    "bias_evidence_summary",
    "Screen-positive counts depend on the screening statistic and should not be reported as calibrated Type-I error or power.",
    "res$components$bias_screen$table"
  )
  add(
    "Interaction contrast wording",
    paste0(
      "No interaction-bias facet pair was selected by mfrm_results(); the ",
      "interaction-bias route status was ", as.character(ev$InteractionStatus[1] %||% "not_available"),
      ". Choose facet_a and facet_b explicitly before writing interaction-bias claims."
    ),
    "bias_evidence_summary; bias_screen guidance",
    "Do not let the wrapper choose the fairness contrast; the analyst must specify the substantive facet pair.",
    "estimate_bias(fit, diagnostics, facet_a = ..., facet_b = ...) -> bias_interaction_report()"
  )
  add(
    "DFF follow-up wording",
    paste0(
      "For group-by-facet functioning claims, run the DFF route explicitly ",
      "and report the method, grouping variable, linking/anchor support, and ",
      "screening thresholds."
    ),
    "bias_evidence_summary",
    "DFF/DIF labels require a documented grouping variable and method; facet-level residual bias alone is not a DFF conclusion.",
    "analyze_dff(); dif_report(); plot_dif_summary()"
  )
  add(
    "Fairness boundary wording",
    paste0(
      "Bias outputs should be described as conditional screening evidence. ",
      "Final fairness, bias, or invariance conclusions require targeted ",
      "contrasts, design context, low-count review, and substantive judgment."
    ),
    "bias_evidence_summary; evidence_boundary",
    "Do not present screen positives or null screens as standalone fairness decisions.",
    "report$evidence_boundary; reporting_checklist()"
  )
  mfrm_report_template_enrich(
    do.call(rbind, rows),
    boundary_type = "screen_not_fairness_decision",
    default_claim_strength = "write_with_caveat"
  )
}

mfrm_report_misfit_bundle <- function(x) {
  comps <- x$components %||% list()
  if ("misfit_review" %in% names(comps) &&
      inherits(comps$misfit_review, "mfrm_misfit_review")) {
    return(list(
      available = TRUE,
      bundle = comps$misfit_review,
      source = "res$components$misfit_review"
    ))
  }
  list(
    available = FALSE,
    bundle = NULL,
    source = "not_requested"
  )
}

mfrm_report_first_column <- function(df, candidates) {
  df <- as.data.frame(df %||% data.frame(), stringsAsFactors = FALSE)
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) return(NULL)
  df[[hit[1]]]
}

mfrm_report_first_numeric_column <- function(df, candidates) {
  col <- mfrm_report_first_column(df, candidates)
  if (is.null(col)) return(numeric(0))
  suppressWarnings(as.numeric(col))
}

mfrm_report_flag_rows <- function(df, candidates) {
  df <- as.data.frame(df %||% data.frame(), stringsAsFactors = FALSE)
  hit <- intersect(candidates, names(df))
  if (nrow(df) == 0L || length(hit) == 0L) return(NA_integer_)
  flags <- lapply(df[hit], function(col) {
    if (is.logical(col)) return(col %in% TRUE)
    value <- tolower(trimws(as.character(col)))
    value %in% c("true", "t", "1", "yes", "y", "flag", "flagged", "review", "warn")
  })
  flag_mat <- do.call(cbind, flags)
  if (is.null(dim(flag_mat))) flag_mat <- matrix(flag_mat, ncol = 1L)
  sum(rowSums(flag_mat, na.rm = TRUE) > 0L, na.rm = TRUE)
}

mfrm_report_misfit_evidence_summary <- function(x) {
  misfit_info <- mfrm_report_misfit_bundle(x)
  if (!isTRUE(misfit_info$available)) {
    return(data.frame(
      Status = "not_requested",
      UnexpectedRows = NA_integer_,
      DisplacementRows = NA_integer_,
      PathwayFitRows = NA_integer_,
      PathwayStatusRows = NA_integer_,
      CurveFitStatusRows = NA_integer_,
      UnexpectedScreenPositiveRows = NA_integer_,
      DisplacementFlaggedRows = NA_integer_,
      MaxAbsStdResidual = NA_real_,
      MaxAbsDisplacement = NA_real_,
      MaxAbsDisplacementT = NA_real_,
      PathwayAvailable = FALSE,
      Source = "not_requested",
      Route = "mfrm_results(fit, include = \"misfit_review\")",
      Boundary = paste0(
        "Misfit/pathway wording requires the misfit_review preset or explicit ",
        "unexpected-response, displacement, and pathway helper calls."
      ),
      stringsAsFactors = FALSE
    ))
  }

  bundle <- misfit_info$bundle
  unexpected <- as.data.frame(bundle$unexpected$table %||% data.frame(), stringsAsFactors = FALSE)
  displacement <- as.data.frame(bundle$displacement$table %||% data.frame(), stringsAsFactors = FALSE)
  pathway_fit <- as.data.frame(bundle$pathway_fit_measures %||% data.frame(), stringsAsFactors = FALSE)
  pathway_status <- as.data.frame(bundle$pathway_fit_status %||% data.frame(), stringsAsFactors = FALSE)
  pathway_curve_status <- as.data.frame(bundle$pathway_curve_fit_status %||% data.frame(), stringsAsFactors = FALSE)

  std_resid <- mfrm_report_first_numeric_column(
    unexpected,
    c("StdResidual", "StandardizedResidual", "ZResidual", "ResidualZ",
      "MeanStdResidual", "AbsStdResidual")
  )
  displacement_value <- mfrm_report_first_numeric_column(
    displacement,
    c("Displacement", "MeasureDisplacement", "DisplacementLogit", "AnchorGap")
  )
  displacement_t <- mfrm_report_first_numeric_column(
    displacement,
    c("DisplacementT", "t_Displacement", "DisplacementZ", "ZDisplacement")
  )
  unexpected_flags <- mfrm_report_flag_rows(
    unexpected,
    c("Flag", "ReviewFlag", "FlagLowProbability", "FlagLargeResidual")
  )
  if (!is.finite(unexpected_flags)) {
    unexpected_flags <- sum(is.finite(std_resid) & abs(std_resid) >= 2, na.rm = TRUE)
  }
  displacement_flags <- mfrm_report_flag_rows(
    displacement,
    c("Flag", "ReviewFlag", "FlagDisplacement", "FlagT")
  )
  if (!is.finite(displacement_flags)) {
    displacement_flags <- sum(is.finite(displacement_t) & abs(displacement_t) >= 2, na.rm = TRUE)
  }
  data.frame(
    Status = "available",
    UnexpectedRows = nrow(unexpected),
    DisplacementRows = nrow(displacement),
    PathwayFitRows = nrow(pathway_fit),
    PathwayStatusRows = nrow(pathway_status),
    CurveFitStatusRows = nrow(pathway_curve_status),
    UnexpectedScreenPositiveRows = as.integer(unexpected_flags),
    DisplacementFlaggedRows = as.integer(displacement_flags),
    MaxAbsStdResidual = if (any(is.finite(std_resid))) max(abs(std_resid), na.rm = TRUE) else NA_real_,
    MaxAbsDisplacement = if (any(is.finite(displacement_value))) max(abs(displacement_value), na.rm = TRUE) else NA_real_,
    MaxAbsDisplacementT = if (any(is.finite(displacement_t))) max(abs(displacement_t), na.rm = TRUE) else NA_real_,
    PathwayAvailable = nrow(pathway_fit) > 0L || nrow(pathway_status) > 0L ||
      nrow(pathway_curve_status) > 0L,
    Source = misfit_info$source,
    Route = paste0(
      "res$components$misfit_review; unexpected_response_table(); ",
      "displacement_table(); plot(fit, type = \"pathway\")"
    ),
    Boundary = paste0(
      "Unexpected responses, displacement rows, and pathway-map fit annotations ",
      "are case-review prompts. They are not automatic exclusion, fairness, or ",
      "validity decisions."
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_report_misfit_reporting_templates <- function(style, misfit_evidence_summary) {
  ev <- as.data.frame(misfit_evidence_summary %||% data.frame(), stringsAsFactors = FALSE)
  audience <- mfrm_report_fit_template_audience(style)
  if (nrow(ev) == 0L || !as.character(ev$Status[1] %||% "") %in% "available") {
    return(mfrm_report_template_enrich(data.frame(
      Audience = audience,
      Topic = "Misfit/pathway evidence not requested",
      Template = paste0(
        "Misfit/pathway wording was not generated because the mfrm_results ",
        "object was not built with include = \"misfit_review\"."
      ),
      EvidenceUsed = "misfit_evidence_summary",
      Caveat = "Request the misfit_review preset or run the local misfit helpers before writing case-review language.",
      Route = "mfrm_results(fit, include = \"misfit_review\"); build_misfit_casebook()",
      Style = style,
      stringsAsFactors = FALSE
    ), boundary_type = "misfit_not_exclusion_rule", default_claim_strength = "not_supported_without_followup"))
  }

  rows <- list()
  add <- function(topic, template, evidence, caveat, route) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Audience = audience,
      Topic = as.character(topic),
      Template = as.character(template),
      EvidenceUsed = as.character(evidence),
      Caveat = as.character(caveat),
      Route = as.character(route),
      Style = style,
      stringsAsFactors = FALSE
    )
  }
  style_lead <- switch(
    style,
    apa = "Report local misfit as follow-up case review only when it was part of the analysis plan.",
    validation = "Use local misfit evidence as interpretive support inside the validity argument, with explicit limits.",
    reviewer = "Use as reviewer-facing evidence about which local misfit traces were inspected.",
    technical = "Use as appendix wording tied to unexpected-response, displacement, and pathway routes.",
    qc = "Use as first-screen case-review wording before manuscript or reviewer handoff.",
    "Use as local misfit wording tied to the stored result object."
  )

  add(
    "Unexpected-response wording",
    paste0(
      style_lead, " The unexpected-response review retained ",
      mfrm_report_fmt_int(ev$UnexpectedRows[1]), " observation-level row(s); ",
      mfrm_report_fmt_int(ev$UnexpectedScreenPositiveRows[1]),
      " row(s) were screen-positive by stored low-probability or large-residual flags, ",
      "and the largest absolute standardized residual was ",
      mfrm_report_fmt_num(ev$MaxAbsStdResidual[1]), "."
    ),
    "misfit_evidence_summary; unexpected_response_table()",
    "Unexpected rows identify observations for review; they are not automatic person, rater, item, or response exclusions.",
    "res$components$misfit_review$unexpected$table; build_misfit_casebook()"
  )
  add(
    "Displacement wording",
    paste0(
      "The displacement review retained ",
      mfrm_report_fmt_int(ev$DisplacementRows[1]), " facet-level row(s); ",
      mfrm_report_fmt_int(ev$DisplacementFlaggedRows[1]),
      " row(s) were flagged by stored displacement criteria. The maximum absolute ",
      "displacement was ", mfrm_report_fmt_num(ev$MaxAbsDisplacement[1]),
      " logits and the maximum absolute displacement t value was ",
      mfrm_report_fmt_num(ev$MaxAbsDisplacementT[1]), "."
    ),
    "misfit_evidence_summary; displacement_table()",
    "Displacement is a measure-sensitivity prompt; inspect anchors, sparseness, and response context before writing substantive conclusions.",
    "res$components$misfit_review$displacement$table"
  )
  add(
    "Pathway-map wording",
    paste0(
      "The pathway-map surface was ",
      if (isTRUE(ev$PathwayAvailable[1])) "available" else "not available",
      " with ", mfrm_report_fmt_int(ev$PathwayFitRows[1]),
      " fit-measure row(s), ", mfrm_report_fmt_int(ev$PathwayStatusRows[1]),
      " facet-status row(s), and ",
      mfrm_report_fmt_int(ev$CurveFitStatusRows[1]), " curve-status row(s)."
    ),
    "misfit_evidence_summary; plot(fit, type = \"pathway\")",
    "A pathway map is visual context for expected-score trajectories and fit annotations, not proof of category functioning by itself.",
    "plot(res, type = \"pathway\", draw = FALSE)"
  )
  add(
    "Case-review wording",
    paste0(
      "For reportable local misfit claims, link unexpected rows and displacement ",
      "signals to a documented casebook review before interpreting individual ",
      "persons, raters, items, criteria, or observations."
    ),
    "misfit_evidence_summary",
    "Case-review prompts should be reconciled with study design, substantive records, and low-count checks before action.",
    "build_misfit_casebook(); unexpected_after_bias_table()"
  )
  add(
    "Boundary wording",
    paste0(
      "Keep local misfit, bias/fairness, fit status, and precision claims in ",
      "separate sentences. A screen-positive observation can motivate review ",
      "without justifying exclusion, fairness, or validity conclusions by itself."
    ),
    "misfit_evidence_summary; evidence_boundary",
    "Do not present unexpected-response, displacement, or pathway flags as automatic exclusion or acceptance rules.",
    "report$evidence_boundary; reporting_checklist()"
  )
  mfrm_report_template_enrich(
    do.call(rbind, rows),
    boundary_type = "misfit_not_exclusion_rule",
    default_claim_strength = "write_with_caveat"
  )
}

mfrm_report_linking_bundle <- function(x) {
  comps <- x$components %||% list()
  if ("linking_review" %in% names(comps) &&
      inherits(comps$linking_review, "mfrm_linking_review")) {
    return(list(
      available = TRUE,
      bundle = comps$linking_review,
      source = "res$components$linking_review"
    ))
  }
  list(
    available = FALSE,
    bundle = NULL,
    source = "not_requested"
  )
}

mfrm_report_linking_guidance_status <- function(review, area, default = "not_available") {
  guidance <- as.data.frame(review$first_screen_guidance %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(guidance) == 0L || !all(c("Area", "Status") %in% names(guidance))) {
    return(default)
  }
  as.character(guidance$Status[as.character(guidance$Area) %in% area][1] %||% default)
}

mfrm_report_sum_numeric_column <- function(df, column) {
  df <- as.data.frame(df %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(df) == 0L || !(column %in% names(df))) return(NA_integer_)
  value <- suppressWarnings(as.numeric(df[[column]]))
  if (!any(is.finite(value))) return(NA_integer_)
  as.integer(sum(value[is.finite(value)], na.rm = TRUE))
}

mfrm_report_linking_evidence_summary <- function(x) {
  linking_info <- mfrm_report_linking_bundle(x)
  if (!isTRUE(linking_info$available)) {
    return(data.frame(
      Status = "not_requested",
      AnchorReviewAvailable = FALSE,
      DriftAvailable = FALSE,
      ChainAvailable = FALSE,
      ReviewStatus = "not_requested",
      TopRiskRows = NA_integer_,
      AnchorRiskRows = NA_integer_,
      DriftRiskRows = NA_integer_,
      ChainRiskRows = NA_integer_,
      GroupViews = NA_integer_,
      AnchorFacetRows = NA_integer_,
      AnchoredLevels = NA_integer_,
      GroupAnchoredLevels = NA_integer_,
      OverlapLevels = NA_integer_,
      AnchorIssueTypes = NA_integer_,
      AnchorIssueRows = NA_integer_,
      LowObservationLevels = NA_integer_,
      LowCategoryRows = NA_integer_,
      DriftReviewStatus = "not_requested",
      EquatingChainStatus = "not_requested",
      GPCMSupport = NA_character_,
      SourceModels = NA_character_,
      Source = "not_requested",
      Route = "mfrm_results(fit, include = \"linking\")",
      Boundary = paste0(
        "Linking/anchor wording requires the linking preset or explicit ",
        "anchor, drift, or equating-chain helper calls."
      ),
      stringsAsFactors = FALSE
    ))
  }

  review <- linking_info$bundle
  overview <- as.data.frame(review$overview %||% data.frame(), stringsAsFactors = FALSE)
  anchor_review <- x$fit$config$anchor_review %||% NULL
  facet_summary <- if (inherits(anchor_review, "mfrm_anchor_review")) {
    as.data.frame(anchor_review$facet_summary %||% data.frame(), stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  design_checks <- if (inherits(anchor_review, "mfrm_anchor_review")) {
    anchor_review$design_checks %||% list()
  } else {
    list()
  }
  issue_counts <- if (inherits(anchor_review, "mfrm_anchor_review")) {
    as.data.frame(anchor_review$issue_counts %||% data.frame(), stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  issue_n <- suppressWarnings(as.numeric(issue_counts$N %||% numeric(0)))
  issue_positive <- is.finite(issue_n) & issue_n > 0
  source_profile <- as.data.frame(review$settings$source_profile %||% data.frame(), stringsAsFactors = FALSE)
  overview_value <- function(column, default) {
    mfrm_report_column_or(overview, column, default)
  }
  data.frame(
    Status = "available",
    AnchorReviewAvailable = isTRUE(overview_value("AnchorReviewAvailable", FALSE)),
    DriftAvailable = isTRUE(overview_value("DriftAvailable", FALSE)),
    ChainAvailable = isTRUE(overview_value("ChainAvailable", FALSE)),
    ReviewStatus = as.character(overview_value("ReviewStatus", "not_available")),
    TopRiskRows = suppressWarnings(as.integer(overview_value("TopRiskRows", NA_integer_))),
    AnchorRiskRows = nrow(as.data.frame(review$prefit_anchor_risks %||% data.frame(), stringsAsFactors = FALSE)),
    DriftRiskRows = nrow(as.data.frame(review$drift_risks %||% data.frame(), stringsAsFactors = FALSE)),
    ChainRiskRows = nrow(as.data.frame(review$chain_risks %||% data.frame(), stringsAsFactors = FALSE)),
    GroupViews = suppressWarnings(as.integer(overview_value("GroupViews", NA_integer_))),
    AnchorFacetRows = nrow(facet_summary),
    AnchoredLevels = mfrm_report_sum_numeric_column(facet_summary, "AnchoredLevels"),
    GroupAnchoredLevels = mfrm_report_sum_numeric_column(facet_summary, "GroupedLevels"),
    OverlapLevels = mfrm_report_sum_numeric_column(facet_summary, "OverlapLevels"),
    AnchorIssueTypes = sum(issue_positive, na.rm = TRUE),
    AnchorIssueRows = if (any(is.finite(issue_n))) as.integer(sum(issue_n[is.finite(issue_n)], na.rm = TRUE)) else NA_integer_,
    LowObservationLevels = nrow(as.data.frame(design_checks$low_observation_levels %||% data.frame(), stringsAsFactors = FALSE)),
    LowCategoryRows = nrow(as.data.frame(design_checks$low_categories %||% data.frame(), stringsAsFactors = FALSE)),
    DriftReviewStatus = mfrm_report_linking_guidance_status(review, "Drift review", "not_available"),
    EquatingChainStatus = mfrm_report_linking_guidance_status(review, "Equating chain", "not_available"),
    GPCMSupport = as.character(overview_value("GPCMSupport", NA_character_)),
    SourceModels = as.character(overview_value("SourceModels", NA_character_)),
    Source = linking_info$source,
    Route = paste0(
      "res$components$linking_review; plot(res, type = \"anchors\"); ",
      "detect_anchor_drift(); build_equating_chain()"
    ),
    Boundary = paste0(
      "Anchor readiness is a first-screen linking support check. Drift and ",
      "equating claims require explicit multi-fit wave/form comparisons; they ",
      "are not inferred from one fitted object."
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_report_linking_reporting_templates <- function(style, linking_evidence_summary) {
  ev <- as.data.frame(linking_evidence_summary %||% data.frame(), stringsAsFactors = FALSE)
  audience <- mfrm_report_fit_template_audience(style)
  if (nrow(ev) == 0L || !as.character(ev$Status[1] %||% "") %in% "available") {
    return(mfrm_report_template_enrich(data.frame(
      Audience = audience,
      Topic = "Linking evidence not requested",
      Template = paste0(
        "Anchor/linking wording was not generated because the mfrm_results ",
        "object was not built with include = \"linking\"."
      ),
      EvidenceUsed = "linking_evidence_summary",
      Caveat = "Request the linking preset or run explicit anchor, drift, or equating-chain helpers before writing linking claims.",
      Route = "mfrm_results(fit, include = \"linking\"); review_mfrm_anchors(); detect_anchor_drift(); build_equating_chain()",
      Style = style,
      stringsAsFactors = FALSE
    ), boundary_type = "anchor_not_drift_absence", default_claim_strength = "not_supported_without_followup"))
  }

  rows <- list()
  add <- function(topic, template, evidence, caveat, route) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Audience = audience,
      Topic = as.character(topic),
      Template = as.character(template),
      EvidenceUsed = as.character(evidence),
      Caveat = as.character(caveat),
      Route = as.character(route),
      Style = style,
      stringsAsFactors = FALSE
    )
  }
  style_lead <- switch(
    style,
    apa = "Report anchor/linking evidence only for the scale-maintenance question addressed by the design.",
    validation = "Use anchor/linking evidence as operational scale-maintenance support inside the validity argument.",
    reviewer = "Use as reviewer-facing evidence about which anchor, drift, and equating checks were available.",
    technical = "Use as appendix wording tied to anchor-review, drift, and chain routes.",
    qc = "Use as first-screen anchor-readiness wording before drift or equating claims.",
    "Use as anchor/linking wording tied to the stored result object."
  )

  add(
    "Anchor-readiness wording",
    paste0(
      style_lead, " The stored linking review classified anchor readiness as ",
      as.character(ev$ReviewStatus[1] %||% "not_available"),
      "; anchor-review metadata was ",
      if (isTRUE(ev$AnchorReviewAvailable[1])) "available" else "not available",
      ", with ", mfrm_report_fmt_int(ev$AnchorFacetRows[1]),
      " facet-summary row(s), ", mfrm_report_fmt_int(ev$AnchoredLevels[1]),
      " directly anchored level(s), ", mfrm_report_fmt_int(ev$GroupAnchoredLevels[1]),
      " group-anchored level(s), and ", mfrm_report_fmt_int(ev$OverlapLevels[1]),
      " overlap level(s)."
    ),
    "linking_evidence_summary; review_mfrm_anchors()",
    "Anchor readiness is a support check for scale maintenance, not proof that forms or waves are already equated.",
    "report$linking_evidence_summary; plot(res, type = \"anchors\")"
  )
  add(
    "Linking-risk wording",
    paste0(
      "The operational linking review retained ",
      mfrm_report_fmt_int(ev$TopRiskRows[1]), " top-risk row(s): ",
      mfrm_report_fmt_int(ev$AnchorRiskRows[1]), " anchor risk row(s), ",
      mfrm_report_fmt_int(ev$DriftRiskRows[1]), " drift risk row(s), and ",
      mfrm_report_fmt_int(ev$ChainRiskRows[1]), " screened-chain risk row(s). ",
      "Anchor schema/value issue rows totaled ",
      mfrm_report_fmt_int(ev$AnchorIssueRows[1]), " across ",
      mfrm_report_fmt_int(ev$AnchorIssueTypes[1]), " issue type(s)."
    ),
    "linking_evidence_summary; res$components$linking_review$top_linking_risks",
    "Risk rows are operational triage prompts; they should not be collapsed into a hidden composite linking score.",
    "summary(res$components$linking_review); build_summary_table_bundle(res$components$linking_review)"
  )
  add(
    "Drift wording",
    paste0(
      "Wave/form drift review status was ",
      as.character(ev$DriftReviewStatus[1] %||% "not_available"),
      "; drift evidence was ",
      if (isTRUE(ev$DriftAvailable[1])) "available" else "not available in this single-result report",
      ". Use separately fitted waves or forms before writing drift claims."
    ),
    "linking_evidence_summary; first_screen_guidance",
    "Do not infer rater, item, criterion, or anchor drift from a single fitted object.",
    "detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2)); plot_anchor_drift(drift)"
  )
  add(
    "Equating-chain wording",
    paste0(
      "Screened equating-chain status was ",
      as.character(ev$EquatingChainStatus[1] %||% "not_available"),
      "; chain evidence was ",
      if (isTRUE(ev$ChainAvailable[1])) "available" else "not available in this single-result report",
      ". Use an ordered list of fitted forms or waves before reporting cumulative offsets."
    ),
    "linking_evidence_summary; first_screen_guidance",
    "A first-screen anchor review is not a completed equating study.",
    "build_equating_chain(list(Form1 = fit1, Form2 = fit2)); plot_anchor_drift(chain, type = \"chain\")"
  )
  add(
    "GPCM boundary wording",
    paste0(
      "Bounded GPCM linking support was recorded as ",
      as.character(ev$GPCMSupport[1] %||% "not_available"),
      ". For bounded GPCM, read direct anchor/drift helper outputs with the capability-matrix caveat."
    ),
    "linking_evidence_summary; support_status",
    "Do not extend RSM/PCM linking synthesis wording to bounded GPCM without an explicit validation route.",
    "gpcm_capability_matrix(); mfrmr_output_guide(\"gpcm\")"
  )
  add(
    "Boundary wording",
    paste0(
      "Keep anchor readiness, drift, equating-chain evidence, and DFF/fairness ",
      "claims in separate sentences. Linking support can justify follow-up ",
      "workflow choices, but it does not by itself establish scale invariance."
    ),
    "linking_evidence_summary; evidence_boundary",
    "Do not report anchor evidence as automatic drift absence, completed equating, DFF support, or validity proof.",
    "report$evidence_boundary; reporting_checklist()"
  )
  mfrm_report_template_enrich(
    do.call(rbind, rows),
    boundary_type = "anchor_not_drift_absence",
    default_claim_strength = "write_with_caveat"
  )
}

mfrm_report_section_row <- function(sections, section) {
  sections <- as.data.frame(sections %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(sections) == 0L || !"Section" %in% names(sections)) return(data.frame())
  hit <- sections[sections$Section %in% section, , drop = FALSE]
  if (nrow(hit) == 0L) data.frame() else hit[1, , drop = FALSE]
}

mfrm_report_section_status <- function(sections, section, default = "not_requested") {
  row <- mfrm_report_section_row(sections, section)
  if (nrow(row) == 0L || !"Status" %in% names(row)) default else as.character(row$Status[1] %||% default)
}

mfrm_report_readiness_label <- function(status) {
  status <- as.character(status %||% "")
  if (status %in% "available") return("ready")
  if (status %in% "review") return("write_with_caveat")
  if (status %in% "caveat") return("caveated")
  if (status %in% "not_requested") return("needs_requested_section")
  if (status %in% "not_available") return("unavailable")
  "review"
}

mfrm_report_claim_readiness <- function(x, sections, style) {
  overview <- as.data.frame(summary(x)$overview %||% data.frame(), stringsAsFactors = FALSE)
  model <- as.character(overview$Model[1] %||% x$fit$config$model %||% "")
  rows <- list()
  add <- function(claim, section, evidence_needed, suggested_wording,
                  follow_up, boundary, status = NULL) {
    if (is.null(status)) status <- mfrm_report_section_status(sections, section)
    rows[[length(rows) + 1L]] <<- data.frame(
      Claim = as.character(claim),
      Section = as.character(section),
      CurrentStatus = as.character(status),
      Readiness = mfrm_report_readiness_label(status),
      EvidenceNeeded = as.character(evidence_needed %||% ""),
      SuggestedWording = as.character(suggested_wording %||% ""),
      FollowUp = as.character(follow_up %||% ""),
      Boundary = as.character(boundary %||% ""),
      Style = as.character(style),
      stringsAsFactors = FALSE
    )
  }

  add(
    "Model specification",
    "Model and data setup",
    "Model, method, facets, score coding, categories, sample size, and missing-data handling.",
    "Report the fitted MFRM specification, estimation method, scoring scale, and facet roles explicitly.",
    "Use specifications_report(fit) and the analysis script for final methods wording.",
    "This documents the analysis setup; it is not validity evidence by itself.",
    status = "available"
  )
  add(
    "Diagnostic review completed",
    "First-screen diagnostics",
    "A diagnostics object, triage rows, and any key warning text.",
    "State that diagnostics were inspected, then report only the specific supported findings.",
    "Inspect summary(res)$triage and summary(res$diagnostics)$key_warnings.",
    "Diagnostic availability is a starting point, not a global quality guarantee."
  )
  add(
    "Fit and precision evidence",
    "Fit, separation, and precision",
    "Selected MnSq threshold profile, observed fit-status counts, ZSTD convention, fit df, df-sensitivity rows, separation, reliability, strata, and uncertainty/context notes.",
    "Report MnSq fit, ZSTD standardization, separation, and reliability as separate evidence streams with the selected threshold profile stated.",
    "Use report$fit_evidence_summary, report$fit_threshold_sensitivity, report$fit_df_sensitivity_summary, precision_review_report(), and facets_fit_df_guide().",
    "Do not reduce these indices to one pass/fail claim, and do not interpret a df-sensitive ZSTD flag without MnSq and context."
  )
  add(
    "Category functioning",
    "Category functioning",
    "Rating-scale, category-structure, or category-curve evidence.",
    "Describe whether score categories behaved as intended and identify any category-level caveats.",
    "Use rating_scale_table(), category_structure_report(), and category_curves_report().",
    "Category evidence supports score-scale review but not a standalone validity claim."
  )
  add(
    "Bias or DFF screening",
    "Bias screening",
    "Facet-level bias table and any explicitly selected interaction or DFF contrast.",
    "Use screening language unless a targeted contrast and substantive review have been completed.",
    "Use mfrm_results(fit, include = \"bias\") and then estimate_bias() for explicit facet pairs.",
    "Do not present screen positives as final fairness conclusions."
  )
  add(
    "Misfit case review",
    "Misfit and pathway review",
    "Unexpected-response rows, displacement evidence, and pathway-map context.",
    "Frame misfit rows as case-review prompts and report the follow-up basis.",
    "Use mfrm_results(fit, include = \"misfit_review\") and build_misfit_casebook() when needed.",
    "Do not use observation-level misfit as an automatic exclusion rule."
  )
  add(
    "Anchor, linking, or drift claim",
    "Anchors and linking",
    "Anchor-readiness output for one fit; multiple fitted waves/forms for drift or equating.",
    "Report anchor readiness separately from drift or equating claims.",
    "Use mfrm_results(fit, include = \"linking\"); for drift/equating use detect_anchor_drift() or build_equating_chain().",
    "Do not infer drift or equating from a single fitted object."
  )
  add(
    "Design connectivity",
    "Network and connectivity",
    "Network/connectivity review and design overlap evidence.",
    "Use connectivity language to describe design support and sparseness, not model fit.",
    "Use mfrm_results(fit, include = \"network\") and build_mfrm_network_review().",
    "Connectivity evidence does not replace fit, precision, or bias diagnostics."
  )
  add(
    "APA-style manuscript text",
    "APA and manuscript wording",
    "Supported APA output object or a manually edited report template.",
    "Treat generated APA-style text as a draft and edit against the study design.",
    "Use mfrm_results(fit, include = \"publication\") or build_apa_outputs().",
    "Generated prose is not a substitute for study-specific reporting judgment."
  )
  add(
    "Appendix or reviewer supplement",
    "Tables, plots, and handoff",
    "Collected result tables, plot routes, replay code, and a written-files manifest if exported.",
    "Provide tables and replay routes so readers can inspect the evidence surface.",
    "Use build_summary_table_bundle(res), export_mfrm_results(res), or mfrm_report(res, output = \"html\").",
    "Appendix files preserve evidence; they do not add new analyses."
  )
  if (identical(toupper(model), "GPCM")) {
    add(
      "Bounded GPCM interpretation",
      "GPCM scope",
      "Capability-matrix row for each helper used with GPCM outputs.",
      "Report GPCM outputs through documented helper-specific caveats.",
      "Use gpcm_capability_matrix() and mfrmr_output_guide(\"gpcm\").",
      "Do not imply full equivalence with every RSM/PCM report route.",
      status = "caveat"
    )
  }
  out <- do.call(rbind, rows)
  rank <- c(unavailable = 1L, needs_requested_section = 2L, write_with_caveat = 3L,
            caveated = 4L, review = 5L, ready = 6L)
  ord <- rank[out$Readiness]
  ord[is.na(ord)] <- 7L
  out[order(ord, out$Claim), , drop = FALSE]
}

mfrm_report_gap_action <- function(section, status, route) {
  if (status %in% "not_requested") {
    if (grepl("Bias", section, fixed = TRUE)) return("Rebuild the result with mfrm_results(fit, include = \"bias\") before writing bias or fairness-screen text.")
    if (grepl("Misfit", section, fixed = TRUE)) return("Rebuild the result with mfrm_results(fit, include = \"misfit_review\") before writing observation-level misfit text.")
    if (grepl("Anchors", section, fixed = TRUE)) return("Rebuild the result with mfrm_results(fit, include = \"linking\") before writing anchor-readiness text.")
    if (grepl("Network", section, fixed = TRUE)) return("Rebuild the result with mfrm_results(fit, include = \"network\") before writing connectivity text.")
    if (grepl("APA", section, fixed = TRUE)) return("Rebuild the result with mfrm_results(fit, include = \"publication\") before using APA-style output.")
    if (grepl("Category", section, fixed = TRUE)) return("Rebuild the result with table/category sections or call the category helper directly.")
    return("Request the relevant mfrm_results() section or call the route-specific helper before reporting this claim.")
  }
  if (status %in% "not_available") {
    return("Document why this section was unavailable and avoid treating the omission as evidence.")
  }
  if (status %in% "review") {
    return("Write only a caveated claim and inspect the route-specific table before manuscript use.")
  }
  if (status %in% "caveat") {
    return("Keep the caveat in the report and cite the helper-specific support boundary.")
  }
  if (nzchar(route)) paste0("Inspect ", route, " before final wording.") else "Inspect the report section before final wording."
}

mfrm_report_gap_type <- function(status) {
  status <- as.character(status %||% "")
  if (status %in% "not_available") return("unavailable")
  if (status %in% "not_requested") return("not_requested")
  if (status %in% "review") return("caveated_evidence")
  if (status %in% "caveat") return("scope_caveat")
  "none"
}

mfrm_report_gaps <- function(sections) {
  sections <- as.data.frame(sections %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(sections) == 0L || !"Status" %in% names(sections)) {
    return(data.frame(
      Priority = 1L,
      GapType = "unavailable",
      Section = "Report",
      CurrentStatus = "not_available",
      RecommendedAction = "Rebuild mfrm_results() before creating a report.",
      Route = "mfrm_results(fit)",
      Reason = "The report section plan was empty.",
      stringsAsFactors = FALSE
    ))
  }
  keep <- sections$Status %in% c("not_available", "not_requested", "review", "caveat")
  gaps <- sections[keep, , drop = FALSE]
  if (nrow(gaps) == 0L) {
    return(data.frame(
      Priority = 1L,
      GapType = "none",
      Section = "Report",
      CurrentStatus = "ready",
      RecommendedAction = "No immediate report gaps were found in the requested sections.",
      Route = "mfrm_report(res)",
      Reason = "All section-plan rows were available without review or caveat status.",
      stringsAsFactors = FALSE
    ))
  }
  priority_rank <- c(not_available = 1L, not_requested = 2L, review = 3L, caveat = 4L)
  priority <- priority_rank[gaps$Status]
  priority[is.na(priority)] <- 5L
  out <- data.frame(
    Priority = as.integer(priority),
    GapType = vapply(gaps$Status, mfrm_report_gap_type, character(1)),
    Section = as.character(gaps$Section),
    CurrentStatus = as.character(gaps$Status),
    RecommendedAction = mapply(
      mfrm_report_gap_action,
      section = as.character(gaps$Section),
      status = as.character(gaps$Status),
      route = as.character(gaps$Route %||% ""),
      USE.NAMES = FALSE
    ),
    Route = as.character(gaps$Route %||% ""),
    Reason = as.character(gaps$Boundary %||% gaps$Evidence %||% ""),
    stringsAsFactors = FALSE
  )
  out[order(out$Priority, out$Section), , drop = FALSE]
}

mfrm_report_narrative <- function(x, sx, style, sections) {
  overview <- as.data.frame(sx$overview %||% data.frame(), stringsAsFactors = FALSE)
  triage <- as.data.frame(sx$triage %||% data.frame(), stringsAsFactors = FALSE)
  top_signal <- if (nrow(triage) > 0L && all(c("Area", "Severity", "Signal") %in% names(triage))) {
    paste(utils::head(paste0(triage$Area, "=", triage$Severity, " (", triage$Signal, ")"), 3L), collapse = "; ")
  } else {
    "No triage rows were available."
  }
  overview_text <- if (nrow(overview) > 0L) {
    paste0(
      "The source result uses model ", as.character(overview$Model[1] %||% ""),
      " with method ", as.character(overview$Method[1] %||% ""),
      " and ", as.character(overview$N[1] %||% ""), " observations."
    )
  } else {
    "The source result did not expose a compact overview table."
  }
  data.frame(
    Paragraph = c("Scope", "Overview", "FirstScreen", "Boundary"),
    Text = c(
      paste0(
        mfrm_report_style_focus(style),
        " The report is generated from an existing mfrm_results object and does not refit the model."
      ),
      overview_text,
      paste0("Highest-priority first-screen signals: ", top_signal, "."),
      paste0(
        "Use the section plan and evidence-boundary table to decide what can be written now, ",
        "what needs a targeted follow-up helper, and what should remain caveated."
      )
    ),
    stringsAsFactors = FALSE
  )
}

mfrm_report_action_items <- function(sx, sections, style) {
  actions <- as.data.frame(sx$next_actions %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(actions) == 0L) {
    actions <- data.frame(
      Priority = 1L,
      Area = "Report",
      Action = "Read the section plan and evidence-boundary table.",
      Route = "mfrm_report(res)",
      Reason = "No next-action table was available from summary(res).",
      stringsAsFactors = FALSE
    )
  }
  actions$ReportDecision <- switch(
    style,
    apa = "Draft manuscript wording only after checking the matching evidence boundary.",
    qc = "Clear QC blockers before report export or reviewer handoff.",
    validation = "Map this action to the validity argument rather than treating it as standalone proof.",
    reviewer = "Use this action to prepare a direct reviewer response or caveat.",
    technical = "Use this action to decide which appendix or reproducibility artifact to include.",
    "Use this action as report follow-up."
  )
  actions
}

mfrm_report_index_status <- function(df, default = "not_available") {
  df <- as.data.frame(df %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(df) == 0L || !"Status" %in% names(df)) return(default)
  as.character(df$Status[1] %||% default)
}

mfrm_report_signal_count <- function(...) {
  values <- unlist(list(...), use.names = FALSE)
  values <- suppressWarnings(as.numeric(values))
  if (!any(is.finite(values))) return(NA_integer_)
  as.integer(sum(values[is.finite(values)], na.rm = TRUE))
}

mfrm_report_index_readiness <- function(status, review_signals) {
  status <- as.character(status %||% "")
  review_signals <- suppressWarnings(as.integer(review_signals[1]))
  if (status %in% "not_requested") return("request_if_needed")
  if (status %in% "not_available") return("unavailable")
  if (status %in% c("available_without_precision_review", "caveat")) return("write_with_caveat")
  if (status %in% c("review", "warn")) return("review")
  if (status %in% "available") {
    if (is.finite(review_signals) && review_signals > 0L) return("review")
    return("ready")
  }
  "review"
}

mfrm_report_plot_route_or <- function(x, type, fallback = "") {
  route <- if (inherits(x, "mfrm_results")) mfrm_report_plot_route(x, type) else ""
  if (length(route) > 0L && nzchar(route[1])) return(as.character(route[1]))
  as.character(fallback %||% "")
}

mfrm_report_index <- function(sections,
                              fit_evidence_summary,
                              precision_evidence_summary,
                              bias_evidence_summary,
                              misfit_evidence_summary,
                              linking_evidence_summary,
                              x = NULL) {
  sections <- as.data.frame(sections %||% data.frame(), stringsAsFactors = FALSE)
  fit <- as.data.frame(fit_evidence_summary %||% data.frame(), stringsAsFactors = FALSE)
  precision <- as.data.frame(precision_evidence_summary %||% data.frame(), stringsAsFactors = FALSE)
  bias <- as.data.frame(bias_evidence_summary %||% data.frame(), stringsAsFactors = FALSE)
  misfit <- as.data.frame(misfit_evidence_summary %||% data.frame(), stringsAsFactors = FALSE)
  linking <- as.data.frame(linking_evidence_summary %||% data.frame(), stringsAsFactors = FALSE)
  rows <- list()
  add <- function(area, section, evidence_status, review_signals, evidence_detail,
                  primary_table, template_table, route, boundary,
                  include_preset, plot_route = "", export_route = "export_mfrm_results(res, include = \"report\")") {
    primary_table <- as.character(primary_table)
    template_table <- as.character(template_table)
    rows[[length(rows) + 1L]] <<- data.frame(
      Area = as.character(area),
      Section = as.character(section),
      SectionStatus = mfrm_report_section_status(sections, section, default = "not_requested"),
      EvidenceStatus = as.character(evidence_status),
      Readiness = mfrm_report_index_readiness(evidence_status, review_signals),
      ReviewSignalCount = suppressWarnings(as.integer(review_signals[1])),
      EvidenceDetail = as.character(evidence_detail %||% ""),
      PrimaryTable = primary_table,
      TemplateTable = template_table,
      Route = as.character(route),
      EvidenceRoute = paste0("report$", primary_table),
      TemplateRoute = paste0("report$", template_table),
      PlotRoute = as.character(plot_route %||% ""),
      ExportRoute = as.character(export_route %||% ""),
      IncludePreset = as.character(include_preset %||% ""),
      Boundary = as.character(boundary %||% ""),
      stringsAsFactors = FALSE
    )
  }

  fit_status <- mfrm_report_index_status(fit)
  fit_signals <- mfrm_report_signal_count(
    mfrm_report_column_or(fit, "UnderfitRows", NA_integer_),
    mfrm_report_column_or(fit, "OverfitRows", NA_integer_),
    mfrm_report_column_or(fit, "MixedRows", NA_integer_),
    mfrm_report_column_or(fit, "DfSensitiveRows", NA_integer_)
  )
  add(
    "Fit",
    "Fit, separation, and precision",
    fit_status,
    fit_signals,
    paste0(
      "underfit=", mfrm_report_fmt_int(mfrm_report_column_or(fit, "UnderfitRows", NA_integer_)),
      "; overfit=", mfrm_report_fmt_int(mfrm_report_column_or(fit, "OverfitRows", NA_integer_)),
      "; df_sensitive=", mfrm_report_fmt_int(mfrm_report_column_or(fit, "DfSensitiveRows", NA_integer_))
    ),
    "fit_evidence_summary",
    "fit_reporting_templates",
    "report$fit_evidence_summary; report$fit_reporting_templates",
    mfrm_report_column_or(fit, "Boundary", "Keep MnSq, ZSTD, df sensitivity, and precision separate."),
    "mfrm_results(fit, include = c(\"fit\", \"diagnostics\", \"precision\", \"reporting\"))",
    mfrm_report_plot_route_or(x, "qc", "plot(res, type = \"qc\")")
  )

  precision_status <- mfrm_report_index_status(precision)
  precision_signals <- mfrm_report_signal_count(
    mfrm_report_column_or(precision, "ReviewOrWarnChecks", NA_integer_),
    mfrm_report_column_or(precision, "ZeroSeparationRows", NA_integer_),
    mfrm_report_column_or(precision, "ZeroReliabilityRows", NA_integer_)
  )
  add(
    "Precision",
    "Fit, separation, and precision",
    precision_status,
    precision_signals,
    paste0(
      "tier=", as.character(mfrm_report_column_or(precision, "PrecisionTier", NA_character_)),
      "; review_warn=", mfrm_report_fmt_int(mfrm_report_column_or(precision, "ReviewOrWarnChecks", NA_integer_)),
      "; reliability_rows=", mfrm_report_fmt_int(mfrm_report_column_or(precision, "ReliabilityRows", NA_integer_))
    ),
    "precision_evidence_summary",
    "precision_reporting_templates",
    "report$precision_evidence_summary; report$precision_reporting_templates",
    mfrm_report_column_or(precision, "Boundary", "Keep separation reliability distinct from agreement, fit, and validity."),
    "mfrm_results(fit, include = c(\"fit\", \"diagnostics\", \"precision\", \"reporting\"))",
    mfrm_report_plot_route_or(x, "qc", "plot(res, type = \"qc\")")
  )

  bias_status <- mfrm_report_index_status(bias, default = "not_requested")
  bias_signals <- mfrm_report_signal_count(
    mfrm_report_column_or(bias, "ResidualTScreenPositiveRows", NA_integer_),
    mfrm_report_column_or(bias, "ChiSqScreenPositiveRows", NA_integer_)
  )
  add(
    "Bias / DFF",
    "Bias screening",
    bias_status,
    bias_signals,
    paste0(
      "screen_rows=", mfrm_report_fmt_int(mfrm_report_column_or(bias, "Rows", NA_integer_)),
      "; residual_screen=", mfrm_report_fmt_int(mfrm_report_column_or(bias, "ResidualTScreenPositiveRows", NA_integer_)),
      "; chi_sq_screen=", mfrm_report_fmt_int(mfrm_report_column_or(bias, "ChiSqScreenPositiveRows", NA_integer_))
    ),
    "bias_evidence_summary",
    "bias_reporting_templates",
    "report$bias_evidence_summary; report$bias_reporting_templates",
    mfrm_report_column_or(bias, "Boundary", "Request targeted contrasts before fairness or DFF language."),
    "mfrm_results(fit, include = \"bias\")",
    mfrm_report_plot_route_or(x, "tables", "plot(res, type = \"tables\")")
  )

  misfit_status <- mfrm_report_index_status(misfit, default = "not_requested")
  misfit_signals <- mfrm_report_signal_count(
    mfrm_report_column_or(misfit, "UnexpectedScreenPositiveRows", NA_integer_),
    mfrm_report_column_or(misfit, "DisplacementFlaggedRows", NA_integer_)
  )
  add(
    "Misfit / pathway",
    "Misfit and pathway review",
    misfit_status,
    misfit_signals,
    paste0(
      "unexpected=", mfrm_report_fmt_int(mfrm_report_column_or(misfit, "UnexpectedRows", NA_integer_)),
      "; displacement=", mfrm_report_fmt_int(mfrm_report_column_or(misfit, "DisplacementRows", NA_integer_)),
      "; pathway=", as.character(mfrm_report_column_or(misfit, "PathwayAvailable", FALSE))
    ),
    "misfit_evidence_summary",
    "misfit_reporting_templates",
    "report$misfit_evidence_summary; report$misfit_reporting_templates",
    mfrm_report_column_or(misfit, "Boundary", "Treat local misfit as case-review prompts."),
    "mfrm_results(fit, include = \"misfit_review\")",
    mfrm_report_plot_route_or(x, "pathway", "plot(res, type = \"pathway\")")
  )

  linking_status <- mfrm_report_index_status(linking, default = "not_requested")
  linking_signals <- mfrm_report_signal_count(
    mfrm_report_column_or(linking, "TopRiskRows", NA_integer_),
    mfrm_report_column_or(linking, "AnchorIssueRows", NA_integer_)
  )
  add(
    "Linking / anchors",
    "Anchors and linking",
    linking_status,
    linking_signals,
    paste0(
      "review_status=", as.character(mfrm_report_column_or(linking, "ReviewStatus", "not_requested")),
      "; drift=", as.character(mfrm_report_column_or(linking, "DriftReviewStatus", "not_requested")),
      "; chain=", as.character(mfrm_report_column_or(linking, "EquatingChainStatus", "not_requested"))
    ),
    "linking_evidence_summary",
    "linking_reporting_templates",
    "report$linking_evidence_summary; report$linking_reporting_templates",
    mfrm_report_column_or(linking, "Boundary", "Use multiple fitted waves or forms for drift/equating claims."),
    "mfrm_results(fit, include = \"linking\")",
    mfrm_report_plot_route_or(x, "anchors", "plot(res, type = \"anchors\")")
  )

  do.call(rbind, rows)
}

mfrm_report_template_index <- function(fit_reporting_templates,
                                       precision_reporting_templates,
                                       bias_reporting_templates,
                                       misfit_reporting_templates,
                                       linking_reporting_templates) {
  specs <- list(
    list(area = "Fit", table = "fit_reporting_templates", data = fit_reporting_templates),
    list(area = "Precision", table = "precision_reporting_templates", data = precision_reporting_templates),
    list(area = "Bias / DFF", table = "bias_reporting_templates", data = bias_reporting_templates),
    list(area = "Misfit / pathway", table = "misfit_reporting_templates", data = misfit_reporting_templates),
    list(area = "Linking / anchors", table = "linking_reporting_templates", data = linking_reporting_templates)
  )
  rows <- list()
  keep <- c(
    "Topic", "BoundaryType", "ClaimStrength", "RecommendedUse",
    "EvidenceTable", "EvidenceRoute", "Route", "Caveat"
  )
  for (spec in specs) {
    tbl <- as.data.frame(spec$data %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(tbl) == 0L) next
    for (nm in keep) {
      if (!nm %in% names(tbl)) tbl[[nm]] <- ""
    }
    out <- tbl[, keep, drop = FALSE]
    out$Area <- spec$area
    out$TemplateTable <- spec$table
    out$TemplateRow <- seq_len(nrow(out))
    rows[[length(rows) + 1L]] <- out[, c(
      "Area", "TemplateTable", "TemplateRow", keep
    ), drop = FALSE]
  }
  if (length(rows) == 0L) {
    return(data.frame(
      Area = character(),
      TemplateTable = character(),
      TemplateRow = integer(),
      Topic = character(),
      BoundaryType = character(),
      ClaimStrength = character(),
      RecommendedUse = character(),
      EvidenceTable = character(),
      EvidenceRoute = character(),
      Route = character(),
      Caveat = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  priority <- c(
    not_supported_without_followup = 1L,
    write_with_caveat = 2L,
    descriptive_only = 3L,
    ready_when_supported = 4L
  )
  ord <- priority[as.character(out$ClaimStrength)]
  ord[is.na(ord)] <- 99L
  out <- out[order(ord, out$Area, out$TemplateRow), , drop = FALSE]
  rownames(out) <- NULL
  out
}

mfrm_report_first_screen_status <- function(readiness) {
  readiness <- as.character(readiness %||% "")
  if (readiness %in% "ready") return("ok")
  if (readiness %in% "review") return("review")
  if (readiness %in% "request_if_needed") return("request_if_needed")
  if (readiness %in% "unavailable") return("unavailable")
  if (readiness %in% "write_with_caveat") return("caveat")
  "review"
}

mfrm_report_first_screen_action <- function(status) {
  switch(
    as.character(status %||% ""),
    ok = "Use the listed template route if this area is reported.",
    review = "Inspect the primary evidence table and template boundary before writing.",
    request_if_needed = "Request this evidence only if the claim is needed.",
    unavailable = "Rebuild the result object or run the listed route before reporting.",
    caveat = "Use caveated wording and inspect the boundary before writing.",
    "Inspect the listed evidence route before writing."
  )
}

mfrm_report_first_screen_issue <- function(row) {
  readiness <- as.character(row$Readiness[1] %||% "")
  signals <- suppressWarnings(as.integer(row$ReviewSignalCount[1] %||% NA_integer_))
  detail <- as.character(row$EvidenceDetail[1] %||% "")
  if (readiness %in% "ready") return("No report-index review signals.")
  if (readiness %in% "request_if_needed") return("Evidence was not requested.")
  if (readiness %in% "unavailable") return("Evidence is unavailable in this report object.")
  if (is.finite(signals) && signals > 0L) {
    return(paste0("ReviewSignalCount = ", signals, "; ", detail))
  }
  if (nzchar(detail)) return(detail)
  as.character(row$Boundary[1] %||% "Review the evidence boundary.")
}

mfrm_report_first_screen <- function(report_index, template_index) {
  report_index <- as.data.frame(report_index %||% data.frame(), stringsAsFactors = FALSE)
  template_index <- as.data.frame(template_index %||% data.frame(), stringsAsFactors = FALSE)
  empty <- data.frame(
    Area = character(),
    Status = character(),
    Readiness = character(),
    MainIssue = character(),
    NextAction = character(),
    PrimaryRoute = character(),
    TemplateRoute = character(),
    PlotRoute = character(),
    BoundaryType = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(report_index) == 0L) return(empty)

  rows <- lapply(seq_len(nrow(report_index)), function(i) {
    row <- report_index[i, , drop = FALSE]
    area <- as.character(row$Area[1] %||% "")
    readiness <- as.character(row$Readiness[1] %||% "")
    status <- mfrm_report_first_screen_status(readiness)
    template_rows <- if (nrow(template_index) > 0L && "Area" %in% names(template_index)) {
      template_index[as.character(template_index$Area) %in% area, , drop = FALSE]
    } else {
      data.frame()
    }
    boundary_type <- if (nrow(template_rows) > 0L && "BoundaryType" %in% names(template_rows)) {
      paste(unique(as.character(template_rows$BoundaryType)), collapse = " | ")
    } else {
      ""
    }
    primary_route <- if (status %in% "request_if_needed" &&
                           nzchar(as.character(row$IncludePreset[1] %||% ""))) {
      as.character(row$IncludePreset[1])
    } else if (nzchar(as.character(row$EvidenceRoute[1] %||% ""))) {
      as.character(row$EvidenceRoute[1])
    } else {
      as.character(row$Route[1] %||% "")
    }
    data.frame(
      Area = area,
      Status = status,
      Readiness = readiness,
      MainIssue = mfrm_report_first_screen_issue(row),
      NextAction = mfrm_report_first_screen_action(status),
      PrimaryRoute = primary_route,
      TemplateRoute = as.character(row$TemplateRoute[1] %||% ""),
      PlotRoute = as.character(row$PlotRoute[1] %||% ""),
      BoundaryType = boundary_type,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  status_rank <- c(unavailable = 1L, review = 2L, caveat = 3L, request_if_needed = 4L, ok = 5L)
  rank <- status_rank[out$Status]
  rank[is.na(rank)] <- 99L
  out <- out[order(rank, out$Area), , drop = FALSE]
  rownames(out) <- NULL

  counts <- table(factor(out$Status, levels = names(status_rank)))
  overall_status <- if (any(out$Status %in% "unavailable")) {
    "unavailable"
  } else if (any(out$Status %in% "review")) {
    "review"
  } else if (any(out$Status %in% "caveat")) {
    "caveat"
  } else if (any(out$Status %in% "request_if_needed")) {
    "request_if_needed"
  } else {
    "ok"
  }
  first_area <- out$Area[1] %||% "Report"
  overall <- data.frame(
    Area = "Overall",
    Status = overall_status,
    Readiness = overall_status,
    MainIssue = paste0(
      "ok=", counts[["ok"]], "; review=", counts[["review"]],
      "; caveat=", counts[["caveat"]], "; request_if_needed=", counts[["request_if_needed"]],
      "; unavailable=", counts[["unavailable"]], "."
    ),
    NextAction = paste0("Start with ", first_area, "."),
    PrimaryRoute = "report$report_index; report$template_index",
    TemplateRoute = "report$template_index",
    PlotRoute = "",
    BoundaryType = "first_screen_summary",
    stringsAsFactors = FALSE
  )
  rbind(overall, out)
}

mfrm_report_markdown_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\n", " ", x, fixed = TRUE)
  x <- gsub("|", "\\|", x, fixed = TRUE)
  x
}

mfrm_report_markdown_table <- function(df, max_rows = 20L) {
  df <- as.data.frame(df %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(df) == 0L || ncol(df) == 0L) return("_No rows._")
  df <- utils::head(df, max_rows)
  df[] <- lapply(df, mfrm_report_markdown_escape)
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  rule <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  body <- apply(df, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  paste(c(header, rule, body), collapse = "\n")
}

mfrm_report_markdown <- function(report) {
  narrative <- as.data.frame(report$narrative %||% data.frame(), stringsAsFactors = FALSE)
  narrative_lines <- if (nrow(narrative) > 0L && "Text" %in% names(narrative)) {
    paste0("- ", narrative$Text)
  } else {
    "- No narrative rows were available."
  }
  paste(
    c(
      paste0("# ", report$title),
      "",
      "## Narrative",
      narrative_lines,
      "",
      "## First Screen",
      mfrm_report_markdown_table(report$first_screen, max_rows = 20L),
      "",
      "## Report Index",
      mfrm_report_markdown_table(report$report_index, max_rows = 20L),
      "",
      "## Template Index",
      mfrm_report_markdown_table(report$template_index, max_rows = 30L),
      "",
      "## Section Plan",
      mfrm_report_markdown_table(report$sections, max_rows = 30L),
      "",
      "## Claim Readiness",
      mfrm_report_markdown_table(report$claim_readiness, max_rows = 30L),
      "",
      "## Report Gaps",
      mfrm_report_markdown_table(report$report_gaps, max_rows = 30L),
      "",
      "## Fit Criteria",
      mfrm_report_markdown_table(report$fit_criteria, max_rows = 30L),
      "",
      "## Fit Evidence Summary",
      mfrm_report_markdown_table(report$fit_evidence_summary, max_rows = 10L),
      "",
      "## Fit Threshold Sensitivity",
      mfrm_report_markdown_table(report$fit_threshold_sensitivity, max_rows = 30L),
      "",
      "## Fit Reporting Templates",
      mfrm_report_markdown_table(report$fit_reporting_templates, max_rows = 20L),
      "",
      "## Precision Evidence Summary",
      mfrm_report_markdown_table(report$precision_evidence_summary, max_rows = 10L),
      "",
      "## Precision Basis",
      mfrm_report_markdown_table(report$precision_basis, max_rows = 10L),
      "",
      "## Precision Reporting Templates",
      mfrm_report_markdown_table(report$precision_reporting_templates, max_rows = 20L),
      "",
      "## Bias Evidence Summary",
      mfrm_report_markdown_table(report$bias_evidence_summary, max_rows = 10L),
      "",
      "## Bias Reporting Templates",
      mfrm_report_markdown_table(report$bias_reporting_templates, max_rows = 20L),
      "",
      "## Misfit Evidence Summary",
      mfrm_report_markdown_table(report$misfit_evidence_summary, max_rows = 10L),
      "",
      "## Misfit Reporting Templates",
      mfrm_report_markdown_table(report$misfit_reporting_templates, max_rows = 20L),
      "",
      "## Linking Evidence Summary",
      mfrm_report_markdown_table(report$linking_evidence_summary, max_rows = 10L),
      "",
      "## Linking Reporting Templates",
      mfrm_report_markdown_table(report$linking_reporting_templates, max_rows = 20L),
      "",
      "## ZSTD Conventions",
      mfrm_report_markdown_table(report$zstd_conventions, max_rows = 30L),
      "",
      "## Fit Decision Policy",
      mfrm_report_markdown_table(report$fit_decision_policy, max_rows = 30L),
      "",
      "## Fit DF Sensitivity",
      mfrm_report_markdown_table(report$fit_df_sensitivity_summary, max_rows = 10L),
      "",
      "## Fit DF Sensitive Rows",
      mfrm_report_markdown_table(report$fit_df_sensitive_rows, max_rows = 10L),
      "",
      "## Evidence Boundary",
      mfrm_report_markdown_table(report$evidence_boundary, max_rows = 30L),
      "",
      "## Next Actions",
      mfrm_report_markdown_table(report$action_items, max_rows = 20L)
    ),
    collapse = "\n"
  )
}

mfrm_report_build <- function(x, style) {
  sx <- summary(
    x,
    top_n = max(10L, nrow(as.data.frame(x$table_index %||% data.frame())) + 1L)
  )
  sections <- mfrm_report_section_plan(x, sx, style)
  evidence <- mfrm_report_evidence_boundary()
  claim_readiness <- mfrm_report_claim_readiness(x, sections, style)
  report_gaps <- mfrm_report_gaps(sections)
  fit_criteria <- mfrm_report_fit_criteria(x)
  fit_evidence_summary <- mfrm_report_fit_evidence_summary(x)
  fit_threshold_sensitivity <- mfrm_report_fit_threshold_sensitivity(x)
  zstd_conventions <- mfrm_report_zstd_conventions()
  fit_decision_policy <- mfrm_report_fit_decision_policy()
  fit_df_sensitivity_summary <- mfrm_report_fit_df_sensitivity_summary(x)
  fit_df_sensitive_rows <- mfrm_report_fit_df_sensitive_rows(x)
  fit_reporting_templates <- mfrm_report_fit_reporting_templates(
    style = style,
    fit_evidence_summary = fit_evidence_summary,
    fit_threshold_sensitivity = fit_threshold_sensitivity,
    fit_df_sensitivity_summary = fit_df_sensitivity_summary
  )
  precision_evidence_summary <- mfrm_report_precision_evidence_summary(x)
  precision_basis <- mfrm_report_precision_basis(x)
  precision_reporting_templates <- mfrm_report_precision_reporting_templates(
    style = style,
    precision_evidence_summary = precision_evidence_summary,
    precision_basis = precision_basis
  )
  bias_evidence_summary <- mfrm_report_bias_evidence_summary(x)
  bias_reporting_templates <- mfrm_report_bias_reporting_templates(
    style = style,
    bias_evidence_summary = bias_evidence_summary
  )
  misfit_evidence_summary <- mfrm_report_misfit_evidence_summary(x)
  misfit_reporting_templates <- mfrm_report_misfit_reporting_templates(
    style = style,
    misfit_evidence_summary = misfit_evidence_summary
  )
  linking_evidence_summary <- mfrm_report_linking_evidence_summary(x)
  linking_reporting_templates <- mfrm_report_linking_reporting_templates(
    style = style,
    linking_evidence_summary = linking_evidence_summary
  )
  report_index <- mfrm_report_index(
    sections = sections,
    fit_evidence_summary = fit_evidence_summary,
    precision_evidence_summary = precision_evidence_summary,
    bias_evidence_summary = bias_evidence_summary,
    misfit_evidence_summary = misfit_evidence_summary,
    linking_evidence_summary = linking_evidence_summary,
    x = x
  )
  template_index <- mfrm_report_template_index(
    fit_reporting_templates = fit_reporting_templates,
    precision_reporting_templates = precision_reporting_templates,
    bias_reporting_templates = bias_reporting_templates,
    misfit_reporting_templates = misfit_reporting_templates,
    linking_reporting_templates = linking_reporting_templates
  )
  first_screen <- mfrm_report_first_screen(
    report_index = report_index,
    template_index = template_index
  )
  narrative <- mfrm_report_narrative(x, sx, style, sections)
  actions <- mfrm_report_action_items(sx, sections, style)
  tables <- list(
    overview = sx$overview,
    first_screen = first_screen,
    report_index = report_index,
    template_index = template_index,
    section_plan = sections,
    claim_readiness = claim_readiness,
    report_gaps = report_gaps,
    fit_criteria = fit_criteria,
    fit_evidence_summary = fit_evidence_summary,
    fit_threshold_sensitivity = fit_threshold_sensitivity,
    fit_reporting_templates = fit_reporting_templates,
    precision_evidence_summary = precision_evidence_summary,
    precision_basis = precision_basis,
    precision_reporting_templates = precision_reporting_templates,
    bias_evidence_summary = bias_evidence_summary,
    bias_reporting_templates = bias_reporting_templates,
    misfit_evidence_summary = misfit_evidence_summary,
    misfit_reporting_templates = misfit_reporting_templates,
    linking_evidence_summary = linking_evidence_summary,
    linking_reporting_templates = linking_reporting_templates,
    zstd_conventions = zstd_conventions,
    fit_decision_policy = fit_decision_policy,
    fit_df_sensitivity_summary = fit_df_sensitivity_summary,
    fit_df_sensitive_rows = fit_df_sensitive_rows,
    evidence_boundary = evidence,
    action_items = actions,
    triage = sx$triage,
    status = sx$status,
    plot_map = sx$plot_map,
    table_index = as.data.frame(x$table_index %||% sx$table_index, stringsAsFactors = FALSE),
    reproducible_code = sx$reproducible_code
  )
  if (nrow(as.data.frame(sx$mapping %||% data.frame())) > 0L) {
    tables$mapping <- sx$mapping
  }
  out <- list(
    style = style,
    title = mfrm_report_title(style),
    source_include = as.character(x$include %||% character(0)),
    summary = sx,
    first_screen = first_screen,
    report_index = report_index,
    template_index = template_index,
    sections = sections,
    claim_readiness = claim_readiness,
    report_gaps = report_gaps,
    fit_criteria = fit_criteria,
    fit_evidence_summary = fit_evidence_summary,
    fit_threshold_sensitivity = fit_threshold_sensitivity,
    fit_reporting_templates = fit_reporting_templates,
    precision_evidence_summary = precision_evidence_summary,
    precision_basis = precision_basis,
    precision_reporting_templates = precision_reporting_templates,
    bias_evidence_summary = bias_evidence_summary,
    bias_reporting_templates = bias_reporting_templates,
    misfit_evidence_summary = misfit_evidence_summary,
    misfit_reporting_templates = misfit_reporting_templates,
    linking_evidence_summary = linking_evidence_summary,
    linking_reporting_templates = linking_reporting_templates,
    zstd_conventions = zstd_conventions,
    fit_decision_policy = fit_decision_policy,
    fit_df_sensitivity_summary = fit_df_sensitivity_summary,
    fit_df_sensitive_rows = fit_df_sensitive_rows,
    evidence_boundary = evidence,
    action_items = actions,
    narrative = narrative,
    tables = tables,
    source = x
  )
  out$markdown <- mfrm_report_markdown(out)
  class(out) <- "mfrm_report"
  out
}

mfrm_report_status_count_table <- function(first_screen) {
  first_screen <- as.data.frame(first_screen %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(first_screen) == 0L || !"Status" %in% names(first_screen)) {
    return(data.frame(Status = character(), Areas = integer(), stringsAsFactors = FALSE))
  }
  rows <- first_screen
  if ("Area" %in% names(rows)) {
    rows <- rows[as.character(rows$Area) != "Overall", , drop = FALSE]
  }
  levels <- c("unavailable", "review", "caveat", "request_if_needed", "ok")
  counts <- table(factor(as.character(rows$Status), levels = levels))
  data.frame(
    Status = names(counts),
    Areas = as.integer(counts),
    stringsAsFactors = FALSE
  )
}

mfrm_report_claim_count_table <- function(claim_readiness) {
  claims <- as.data.frame(claim_readiness %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(claims) == 0L || !"Readiness" %in% names(claims)) {
    return(data.frame(Readiness = character(), Claims = integer(), ExampleClaim = character(), stringsAsFactors = FALSE))
  }
  readiness <- unique(as.character(claims$Readiness))
  out <- lapply(readiness, function(level) {
    rows <- claims[as.character(claims$Readiness) %in% level, , drop = FALSE]
    data.frame(
      Readiness = level,
      Claims = nrow(rows),
      ExampleClaim = as.character(rows$Claim[1] %||% ""),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

mfrm_report_route_table <- function() {
  data.frame(
    Surface = c(
      "First screen",
      "Detailed evidence routes",
      "Template boundaries",
      "Markdown report",
      "Download folder"
    ),
    Route = c(
      "report$first_screen",
      "report$report_index",
      "report$template_index",
      "mfrm_report(res, output = \"markdown\")",
      "export_mfrm_results(res, include = \"report\")"
    ),
    Use = c(
      "Start here; read the overall row and next-action route.",
      "Open when the first screen points to a specific evidence area.",
      "Open before using APA/QC/validation wording.",
      "Use after the evidence status and caveats have been reviewed.",
      "Use when report tables, Markdown, and HTML need to be handed off."
    ),
    stringsAsFactors = FALSE
  )
}

#' @export
summary.mfrm_report <- function(object, top_n = 8, ...) {
  if (!inherits(object, "mfrm_report")) {
    stop("`object` must be an mfrm_report object.", call. = FALSE)
  }
  top_n <- max(1L, as.integer(top_n))
  first_screen <- as.data.frame(object$first_screen %||% data.frame(), stringsAsFactors = FALSE)
  first_rows <- first_screen
  if (nrow(first_rows) > 0L) {
    first_rows <- utils::head(first_rows[, intersect(c(
      "Area", "Status", "Readiness", "MainIssue", "NextAction", "PrimaryRoute"
    ), names(first_rows)), drop = FALSE], top_n)
  }
  detail_rows <- first_screen
  if (nrow(detail_rows) > 0L && "Area" %in% names(detail_rows)) {
    detail_rows <- detail_rows[as.character(detail_rows$Area) != "Overall", , drop = FALSE]
  }
  status_counts <- mfrm_report_status_count_table(first_screen)
  status_value <- function(level) {
    idx <- match(level, status_counts$Status)
    if (is.na(idx)) 0L else as.integer(status_counts$Areas[idx])
  }
  overall <- if (nrow(first_screen) > 0L && "Area" %in% names(first_screen)) {
    first_screen[as.character(first_screen$Area) == "Overall", , drop = FALSE]
  } else {
    data.frame()
  }
  if (nrow(overall) == 0L && nrow(first_screen) > 0L) {
    overall <- first_screen[1, , drop = FALSE]
  }
  overview <- data.frame(
    Style = as.character(object$style %||% ""),
    OverallStatus = as.character(overall$Status[1] %||% ""),
    FirstAction = as.character(overall$NextAction[1] %||% ""),
    ReviewAreas = status_value("review"),
    CaveatAreas = status_value("caveat"),
    OptionalAreas = status_value("request_if_needed"),
    UnavailableAreas = status_value("unavailable"),
    OkAreas = status_value("ok"),
    SourceInclude = paste(as.character(object$source_include %||% character(0)), collapse = ", "),
    stringsAsFactors = FALSE
  )

  immediate_actions <- if (nrow(detail_rows) > 0L && "Status" %in% names(detail_rows)) {
    rows <- detail_rows[as.character(detail_rows$Status) %in% c("unavailable", "review", "caveat"), , drop = FALSE]
    utils::head(rows[, intersect(c(
      "Area", "Status", "MainIssue", "NextAction", "PrimaryRoute", "TemplateRoute"
    ), names(rows)), drop = FALSE], top_n)
  } else {
    data.frame()
  }
  optional_sections <- if (nrow(detail_rows) > 0L && "Status" %in% names(detail_rows)) {
    rows <- detail_rows[as.character(detail_rows$Status) %in% "request_if_needed", , drop = FALSE]
    utils::head(rows[, intersect(c(
      "Area", "Status", "MainIssue", "NextAction", "PrimaryRoute"
    ), names(rows)), drop = FALSE], top_n)
  } else {
    data.frame()
  }
  claim_readiness <- mfrm_report_claim_count_table(object$claim_readiness)
  report_gaps <- as.data.frame(object$report_gaps %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(report_gaps) > 0L) {
    report_gaps <- utils::head(report_gaps[, intersect(c(
      "Priority", "GapType", "Section", "RecommendedAction", "Route"
    ), names(report_gaps)), drop = FALSE], top_n)
  }
  template_index <- as.data.frame(object$template_index %||% data.frame(), stringsAsFactors = FALSE)
  boundary_index <- if (nrow(template_index) > 0L) {
    rows <- template_index
    if ("ClaimStrength" %in% names(rows)) {
      priority <- c(
        not_supported_without_followup = 1L,
        caveated_interpretation = 2L,
        descriptive_only = 3L,
        ready_when_supported = 4L
      )
      ord <- priority[as.character(rows$ClaimStrength)]
      ord[is.na(ord)] <- 99L
      rows <- rows[order(ord), , drop = FALSE]
    }
    utils::head(rows[, intersect(c(
      "Area", "Topic", "BoundaryType", "ClaimStrength", "RecommendedUse", "EvidenceRoute"
    ), names(rows)), drop = FALSE], top_n)
  } else {
    data.frame()
  }

  out <- list(
    overview = overview,
    first_screen = first_rows,
    status_counts = status_counts,
    immediate_actions = immediate_actions,
    optional_sections = optional_sections,
    claim_readiness = claim_readiness,
    report_gaps = report_gaps,
    boundary_index = boundary_index,
    routes = mfrm_report_route_table(),
    top_n = top_n
  )
  class(out) <- "summary.mfrm_report"
  out
}

mfrm_report_html_guidance <- function(summary_obj) {
  overview <- as.data.frame(summary_obj$overview %||% data.frame(), stringsAsFactors = FALSE)
  status <- as.character(overview$OverallStatus[1] %||% "")
  action <- as.character(overview$FirstAction[1] %||% "")
  c(
    paste0("Overall status: ", if (nzchar(status)) status else "not available"),
    paste0("First action: ", if (nzchar(action)) action else "Inspect report$first_screen."),
    "Read order: Report Summary Overview -> Report Summary First Screen -> Immediate Actions -> Optional Sections -> Detailed report tables.",
    "Interpretation boundary: this HTML report summarizes existing mfrm_results() evidence; it does not refit the model or create a new pass/fail decision."
  )
}

mfrm_report_html_tables <- function(report) {
  sx <- summary(report)
  summary_tables <- list(
    report_summary_overview = sx$overview,
    report_summary_first_screen = sx$first_screen,
    report_summary_status_counts = sx$status_counts,
    report_summary_immediate_actions = sx$immediate_actions,
    report_summary_optional_sections = sx$optional_sections,
    report_summary_claim_readiness = sx$claim_readiness,
    report_summary_report_gaps = sx$report_gaps,
    report_summary_boundary_index = sx$boundary_index,
    report_summary_routes = sx$routes
  )
  c(summary_tables, report$tables %||% list())
}

mfrm_report_html <- function(report) {
  if (!inherits(report, "mfrm_report")) {
    stop("`report` must be an mfrm_report object.", call. = FALSE)
  }
  sx <- summary(report)
  html <- build_mfrm_bundle_html(
    title = report$title,
    tables = mfrm_report_html_tables(report),
    text_sections = list("Reader guidance" = mfrm_report_html_guidance(sx)),
    text_sections_after = list("Report Markdown" = report$markdown)
  )
  path <- tempfile("mfrmr_report_", fileext = ".html")
  writeLines(enc2utf8(html), con = path, useBytes = TRUE)
  out <- list(
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    report = report,
    html = html
  )
  class(out) <- "mfrm_report_html"
  out
}

#' Build report-ready output from `mfrm_results()`
#'
#' @description
#' `mfrm_report()` is a report-synthesis layer for an existing
#' [mfrm_results()] object. It does not refit the model, recompute diagnostics,
#' or add new validity rules. Instead, it turns the comprehensive first-screen
#' result into a first-screen table, section plan, claim-readiness table,
#' report-gap table, report-index table, template-index table,
#' fit-criteria table, result-specific fit evidence summaries,
#' fit-reporting wording templates, precision/separation reporting templates,
#' bias/DFF reporting templates, misfit/pathway reporting templates,
#' linking/anchor reporting templates, ZSTD-convention table,
#' evidence-boundary table, next-action table, and optional Markdown or HTML
#' report.
#'
#' @param x An [mfrm_results()] object.
#' @param style Report emphasis. `"qc"` is the default first-screen report.
#'   `"apa"` emphasizes manuscript wording, `"validation"` emphasizes the
#'   validity-argument boundary, `"reviewer"` emphasizes reviewer response
#'   preparation, and `"technical"` emphasizes appendix/reproducibility routes.
#' @param output Return format: `"object"` for an `mfrm_report` object,
#'   `"markdown"` for a character scalar, `"html"` for a temporary HTML file,
#'   or `"tables"` for the report's named data-frame list.
#'
#' @details
#' The intended workflow is:
#' 1. Create `res <- mfrm_results(fit, include = ...)`.
#' 2. Inspect `summary(res)$triage` and `summary(res)$next_actions`.
#' 3. Create `report <- mfrm_report(res, style = "qc")`.
#' 4. Read `summary(report)` and `report$first_screen` before opening detailed
#'    report tables.
#' 5. Use `report$report_index` to choose the next `PrimaryTable`,
#'    `TemplateTable`, plot route, or export route.
#' 6. Use `report$template_index` before copying APA/QC/validation wording.
#' 7. Use `style = "apa"`, `"validation"`, `"reviewer"`, or `"technical"` only
#'    when that reporting question is needed.
#'
#' Report rows deliberately distinguish evidence from claims. The
#' `first_screen` table is the compact entry point: it gives an overall row and
#' one row per major evidence area with status, readiness, main issue, next
#' action, and primary route. The
#' `summary.mfrm_report` method summarizes that first screen into immediate
#' actions, optional not-requested sections, claim-readiness counts, report
#' gaps, and template-boundary rows without introducing a new pass/fail
#' decision. The default print method follows the same short reading order and
#' does not print every detailed evidence table. HTML output places the same
#' reader guidance and report-summary tables before the full Markdown text so
#' the browser view starts from the first-screen route. The
#' `report_index` table is the detailed evidence-route index: it lists the
#' major report areas, evidence status, readiness label, review-signal count,
#' and the primary/template tables, evidence routes, template routes, plot
#' routes, export route, and `mfrm_results(include = ...)` preset to inspect
#' next. In ordinary use, open detailed tables through the `PrimaryTable` and
#' `TemplateTable` columns rather than scanning every element of `report$tables`.
#' The
#' `template_index` table then stacks all fit, precision, bias, misfit, and
#' linking wording templates into a single boundary/claim-strength index before
#' users drill into the area-specific template tables. The
#' `claim_readiness` table marks which report claims are ready, caveated,
#' unavailable, or require additional requested sections. The `report_gaps`
#' table turns those statuses into follow-up actions. The fit-specific tables
#' keep multiple MnSq threshold profiles, observed fit-status counts, and
#' engine-vs-FACETS-style ZSTD conventions visible, including the
#' small-df/capping boundary used for FACETS-style ZSTD review. They summarize
#' the stored `fit_measures` component from `mfrm_results()`; `mfrm_report()`
#' itself does not recompute diagnostics. The `fit_reporting_templates` table
#' turns those counts into cautious APA/QC/validation/reviewer wording
#' scaffolds while keeping MnSq, ZSTD standardization, df sensitivity, and
#' separation/reliability in separate sentences. All reporting-template tables
#' share `EvidenceTable`, `EvidenceRoute`, `BoundaryType`, `ClaimStrength`, and
#' `RecommendedUse` columns so each template can be traced back to its evidence
#' and claim boundary. `template_index` stacks those columns across all
#' template areas so report authors can review unsupported or caveated wording
#' before opening the full template text. The `precision_reporting_templates`
#' table does the same for separation, reliability, and strata using the
#' stored precision review and `diagnostics$reliability`. The
#' `bias_reporting_templates` table is
#' available when the source result was built with `include = "bias"` and keeps
#' facet-level screens, interaction-bias contrasts, DFF follow-up, and fairness
#' conclusions in separate lanes. The `misfit_reporting_templates` table is
#' available when the source result was built with `include = "misfit_review"`
#' and keeps unexpected responses, displacement, pathway-map evidence, and
#' case-review actions separate. The `linking_reporting_templates` table is
#' available when the source result was built with `include = "linking"` and
#' keeps anchor readiness, drift review, equating-chain review, and GPCM
#' support boundaries separate. For example, fit and separation are not
#' collapsed into a single pass/fail statement; bias screens are not treated
#' as final fairness conclusions; pathway/misfit rows are case-review prompts;
#' and drift/equating claims require multiple fitted forms or waves.
#'
#' @return Depending on `output`, an `mfrm_report` object, a Markdown character
#'   scalar, an `mfrm_report_html` object, or a named list of data frames.
#' @seealso [mfrm_results()], [export_mfrm_results()],
#'   [build_apa_outputs()], [reporting_checklist()],
#'   [mfrmr_output_guide()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' toy_small <- toy[toy$Person %in% unique(toy$Person)[1:6], , drop = FALSE]
#' fit <- fit_mfrm(toy_small, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' res <- mfrm_results(fit, include = c("fit", "diagnostics", "tables"))
#'
#' report <- mfrm_report(res, style = "qc")
#' summary(report)
#' report$first_screen
#' report$report_index[, c("Area", "Readiness", "PrimaryTable",
#'                         "TemplateTable", "PlotRoute")]
#' report$template_index[, c("Area", "Topic", "BoundaryType",
#'                           "ClaimStrength", "EvidenceRoute")]
#'
#' # Open detailed evidence only after the index points to it.
#' fit_primary <- report$report_index$PrimaryTable[
#'   report$report_index$Area == "Fit"
#' ][1]
#' report$tables[[fit_primary]]
#'
#' mfrm_report(res, output = "markdown")
#' mfrm_report(res, output = "html")
#' }
#' @export
mfrm_report <- function(x,
                        style = c("qc", "apa", "validation", "reviewer", "technical"),
                        output = c("object", "markdown", "html", "tables")) {
  if (!inherits(x, "mfrm_results")) {
    stop("`x` must be an mfrm_results object. Call `mfrm_results()` first.", call. = FALSE)
  }
  style <- match.arg(tolower(as.character(style[1])), c("qc", "apa", "validation", "reviewer", "technical"))
  output <- match.arg(tolower(as.character(output[1])), c("object", "markdown", "html", "tables"))
  report <- mfrm_report_build(x, style = style)
  switch(
    output,
    object = report,
    markdown = report$markdown,
    html = mfrm_report_html(report),
    tables = report$tables
  )
}

mfrm_results_export_include <- function(include) {
  include <- unique(tolower(as.character(include %||% "default")))
  default <- c("summary", "tables", "html", "rds", "replay", "manifest")
  allowed <- c(default, "report", "plots", "all", "default")
  bad <- setdiff(include, allowed)
  if (length(bad) > 0L) {
    stop(
      "Unsupported `include` values: ", paste(bad, collapse = ", "),
      ". Allowed: ", paste(setdiff(allowed, c("default", "all")), collapse = ", "),
      ", default, all.",
      call. = FALSE
    )
  }
  if ("all" %in% include) {
    return(c(default, "report", "plots"))
  }
  include <- unique(c(if ("default" %in% include) default else character(0), setdiff(include, "default")))
  if (length(include) == 0L) {
    stop("`include` must contain at least one export component.", call. = FALSE)
  }
  include
}

mfrm_results_export_table <- function(x) {
  out <- as.data.frame(x %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(out) == 0L && ncol(out) == 0L) {
    return(out)
  }
  for (nm in names(out)) {
    if (is.list(out[[nm]])) {
      out[[nm]] <- vapply(out[[nm]], function(value) {
        paste(deparse(value, width.cutoff = 80L), collapse = " ")
      }, character(1))
    }
  }
  out
}

mfrm_results_export_summary_tables <- function(x) {
  sx <- summary(
    x,
    top_n = max(10L, nrow(as.data.frame(x$table_index %||% data.frame())) + 1L)
  )
  tables <- list(
    overview = sx$overview,
    status = sx$status,
    component_index = sx$component_index,
    table_index = as.data.frame(x$table_index %||% sx$table_index, stringsAsFactors = FALSE),
    plot_map = sx$plot_map,
    triage = sx$triage,
    next_actions = sx$next_actions,
    mapping = sx$mapping,
    reproducible_code = sx$reproducible_code
  )
  notes <- as.character(sx$notes %||% character(0))
  if (length(notes) > 0L) {
    tables$notes <- data.frame(
      Line = seq_along(notes),
      Note = notes,
      stringsAsFactors = FALSE
    )
  }
  tables
}

mfrm_results_export_write_csv <- function(df, path) {
  utils::write.csv(
    mfrm_results_export_table(df),
    file = path,
    row.names = FALSE,
    na = ""
  )
}

mfrm_results_export_add_written <- function(written_files, component, format, path, note = "") {
  rbind(
    written_files,
    data.frame(
      Component = as.character(component),
      Format = as.character(format),
      Path = normalizePath(path, winslash = "/", mustWork = FALSE),
      Note = as.character(note %||% ""),
      stringsAsFactors = FALSE
    )
  )
}

#' Export a lightweight mfrm_results archive
#'
#' @description
#' `export_mfrm_results()` writes the contents of an existing [mfrm_results()]
#' object to a small shareable folder. It is a results-download helper for the
#' comprehensive first-screen workflow, not a new estimation, diagnostics, or
#' validation step.
#'
#' @param x An [mfrm_results()] object.
#' @param output_dir Directory where files should be written.
#' @param prefix File-name prefix. Non-alphanumeric characters are converted to
#'   underscores.
#' @param include Export components. `"default"` expands to `"summary"`,
#'   `"tables"`, `"html"`, `"rds"`, `"replay"`, and `"manifest"`. Add
#'   `"report"` to write [mfrm_report()] tables plus Markdown and HTML; add
#'   `"plots"` to write available plot routes as PNG files, or use `"all"`.
#' @param overwrite Logical; if `FALSE`, existing files stop the export.
#' @param zip_bundle Logical; if `TRUE`, create a best-effort zip archive of
#'   the written files.
#' @param zip_name Optional zip file name. When omitted,
#'   `{prefix}_mfrm_results.zip` is used.
#' @param plot_width,plot_height,plot_res PNG device settings used when
#'   `include` contains `"plots"`.
#'
#' @details
#' The helper writes:
#' - summary CSVs from `summary(x)` such as overview, status, triage, plot
#'   routes, next actions, mapping, and replay-code lines;
#' - collected `x$tables` as CSV files;
#' - optional report artifacts from `mfrm_report(x)`, including report-index,
#'   evidence-summary, and reporting-template CSVs plus Markdown and HTML;
#' - a lightweight HTML report equivalent to `mfrm_results(x, output = "html")`
#'   for the already-created object;
#' - an `.rds` copy of the `mfrm_results` object;
#' - a replay `.R` scaffold from `x$input$reproducible_code`;
#' - a written-files manifest and compact export summary.
#'
#' Plot export is intentionally optional because some plot routes can be
#' comparatively slow or require richer graphics devices. Plot failures are
#' recorded in the returned `plot_errors` table rather than stopping the export.
#'
#' @return An `mfrm_results_export` object with `summary`, `written_files`,
#'   `plot_errors`, and zip status fields.
#' @seealso [mfrm_results()], [launch_mfrmr_viewer()],
#'   [export_mfrm_bundle()], [export_summary_appendix()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' toy_small <- toy[toy$Person %in% unique(toy$Person)[1:6], , drop = FALSE]
#' fit <- fit_mfrm(toy_small, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' res <- mfrm_results(fit, include = c("fit", "diagnostics", "tables"))
#'
#' exported <- export_mfrm_results(
#'   res,
#'   output_dir = tempdir(),
#'   prefix = "mfrmr_results_example",
#'   overwrite = TRUE
#' )
#' exported$summary[, c("FilesWritten", "CsvWritten", "HtmlWritten")]
#' }
#' @export
export_mfrm_results <- function(x,
                                output_dir = ".",
                                prefix = "mfrmr_results",
                                include = "default",
                                overwrite = FALSE,
                                zip_bundle = FALSE,
                                zip_name = NULL,
                                plot_width = 1200,
                                plot_height = 900,
                                plot_res = 144) {
  if (!inherits(x, "mfrm_results")) {
    stop("`x` must be an mfrm_results object. Call `mfrm_results()` first.", call. = FALSE)
  }
  include <- mfrm_results_export_include(include)
  overwrite <- isTRUE(overwrite)
  zip_bundle <- isTRUE(zip_bundle)
  prefix <- export_sanitize_component_tag(prefix, fallback = "mfrmr_results")
  output_dir <- as.character(output_dir[1] %||% ".")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir, call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  written_files <- data.frame(
    Component = character(0),
    Format = character(0),
    Path = character(0),
    Note = character(0),
    stringsAsFactors = FALSE
  )
  plot_errors <- data.frame(
    Plot = character(0),
    Error = character(0),
    stringsAsFactors = FALSE
  )

  ensure_path <- function(filename) {
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("File already exists: ", path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    path
  }
  add_written <- function(component, format, path, note = "") {
    written_files <<- mfrm_results_export_add_written(written_files, component, format, path, note)
    invisible(path)
  }
  write_csv <- function(df, filename, component, note = "") {
    path <- ensure_path(filename)
    mfrm_results_export_write_csv(df, path)
    add_written(component, "csv", path, note)
  }
  write_text <- function(text, filename, component, format = "txt", note = "") {
    path <- ensure_path(filename)
    writeLines(enc2utf8(as.character(text %||% "")), con = path, useBytes = TRUE)
    add_written(component, format, path, note)
  }

  if ("summary" %in% include) {
    summary_tables <- mfrm_results_export_summary_tables(x)
    tags <- make.unique(
      vapply(names(summary_tables), export_sanitize_component_tag, character(1), fallback = "summary"),
      sep = "_"
    )
    for (i in seq_along(summary_tables)) {
      write_csv(
        summary_tables[[i]],
        paste0(prefix, "_summary_", tags[[i]], ".csv"),
        paste0("summary_", tags[[i]])
      )
    }
  }

  if ("tables" %in% include) {
    table_list <- x$tables %||% list()
    tags <- make.unique(
      vapply(names(table_list), export_sanitize_component_tag, character(1), fallback = "table"),
      sep = "_"
    )
    for (i in seq_along(table_list)) {
      write_csv(
        table_list[[i]],
        paste0(prefix, "_table_", tags[[i]], ".csv"),
        paste0("table_", tags[[i]])
      )
    }
  }

  if ("report" %in% include) {
    report <- mfrm_report(x, style = "qc")
    report_tables <- report$tables %||% list()
    tags <- make.unique(
      vapply(names(report_tables), export_sanitize_component_tag, character(1), fallback = "report_table"),
      sep = "_"
    )
    for (i in seq_along(report_tables)) {
      write_csv(
        report_tables[[i]],
        paste0(prefix, "_report_", tags[[i]], ".csv"),
        paste0("report_", tags[[i]]),
        note = "Table from mfrm_report(x, style = \"qc\")."
      )
    }
    write_text(
      report$markdown,
      paste0(prefix, "_report.md"),
      "report_markdown",
      format = "md",
      note = "Markdown from mfrm_report(x, style = \"qc\", output = \"markdown\")."
    )
    report_html <- mfrm_report_html(report)
    html_path <- ensure_path(paste0(prefix, "_report.html"))
    writeLines(enc2utf8(as.character(report_html$html %||% "")), con = html_path, useBytes = TRUE)
    add_written(
      "report_html",
      "html",
      html_path,
      "HTML from mfrm_report(x, style = \"qc\", output = \"html\")."
    )
  }

  if ("html" %in% include) {
    html_obj <- mfrm_results_html(x)
    html_path <- ensure_path(paste0(prefix, "_results.html"))
    ok <- file.copy(html_obj$path, html_path, overwrite = overwrite)
    if (!isTRUE(ok)) {
      stop("Could not write HTML export: ", html_path, call. = FALSE)
    }
    add_written("results_html", "html", html_path)
  }

  if ("rds" %in% include) {
    rds_path <- ensure_path(paste0(prefix, "_results.rds"))
    saveRDS(x, rds_path)
    add_written("results_rds", "rds", rds_path)
  }

  if ("replay" %in% include) {
    replay_code <- as.character(x$input$reproducible_code %||% "")
    if (!nzchar(replay_code)) {
      replay_code <- paste(as.character(summary(x)$reproducible_code$Code %||% ""), collapse = "\n")
    }
    write_text(
      c(
        "# Replay scaffold generated by export_mfrm_results().",
        "# Review data paths, model settings, and reporting choices before rerunning.",
        "",
        replay_code
      ),
      paste0(prefix, "_replay.R"),
      "replay_code",
      format = "R"
    )
  }

  if ("plots" %in% include) {
    plot_map <- as.data.frame(x$plot_map %||% data.frame(), stringsAsFactors = FALSE)
    plot_types <- if (nrow(plot_map) > 0L && all(c("Type", "Available") %in% names(plot_map))) {
      unique(as.character(plot_map$Type[plot_map$Available %in% TRUE]))
    } else {
      character(0)
    }
    for (type in plot_types) {
      tag <- export_sanitize_component_tag(type, fallback = "plot")
      plot_path <- ensure_path(paste0(prefix, "_plot_", tag, ".png"))
      dev_before <- grDevices::dev.cur()
      result <- tryCatch(
        {
          grDevices::png(filename = plot_path, width = plot_width, height = plot_height, res = plot_res)
          plot(x, type = type)
          TRUE
        },
        error = function(e) e,
        finally = {
          while (grDevices::dev.cur() != dev_before && grDevices::dev.cur() > 1L) {
            grDevices::dev.off()
          }
        }
      )
      if (isTRUE(result) && file.exists(plot_path)) {
        add_written(paste0("plot_", tag), "png", plot_path)
      } else {
        if (file.exists(plot_path)) unlink(plot_path)
        plot_errors <- rbind(
          plot_errors,
          data.frame(
            Plot = type,
            Error = if (inherits(result, "error")) conditionMessage(result) else "Plot file was not created.",
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  if ("manifest" %in% include) {
    export_summary_path <- ensure_path(paste0(prefix, "_export_summary.csv"))
    export_summary <- data.frame(
      Created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      PackageVersion = as.character(utils::packageVersion("mfrmr")),
      Prefix = prefix,
      Include = paste(include, collapse = ","),
      FilesWritten = nrow(written_files) + 2L,
      CsvWritten = sum(written_files$Format %in% "csv") + 2L,
      HtmlWritten = sum(written_files$Format %in% "html"),
      RdsWritten = sum(written_files$Format %in% "rds"),
      ReplayWritten = sum(written_files$Format %in% "R"),
      PlotWritten = sum(written_files$Format %in% "png"),
      PlotErrors = nrow(plot_errors),
      stringsAsFactors = FALSE
    )
    mfrm_results_export_write_csv(export_summary, export_summary_path)
    add_written("export_summary", "csv", export_summary_path)

    manifest_path <- ensure_path(paste0(prefix, "_written_files.csv"))
    manifest_rows <- mfrm_results_export_add_written(
      written_files,
      "written_files",
      "csv",
      manifest_path,
      "Manifest of files written by export_mfrm_results()."
    )
    mfrm_results_export_write_csv(manifest_rows, manifest_path)
    written_files <- manifest_rows

    if (nrow(plot_errors) > 0L) {
      plot_error_path <- ensure_path(paste0(prefix, "_plot_errors.csv"))
      mfrm_results_export_write_csv(plot_errors, plot_error_path)
      add_written("plot_errors", "csv", plot_error_path)
    }
  }

  zip_written <- FALSE
  zip_path <- NULL
  zip_note <- NULL
  if (isTRUE(zip_bundle) && nrow(written_files) > 0L) {
    zip_file <- if (is.null(zip_name) || !nzchar(as.character(zip_name[1] %||% ""))) {
      paste0(prefix, "_mfrm_results.zip")
    } else {
      as.character(zip_name[1])
    }
    zip_path <- ensure_path(zip_file)
    zip_inputs <- unique(normalizePath(written_files$Path, winslash = "/", mustWork = TRUE))
    zip_result <- tryCatch(
      {
        utils::zip(zipfile = zip_path, files = zip_inputs, extras = "-j")
        TRUE
      },
      error = function(e) e
    )
    if (isTRUE(zip_result) && file.exists(zip_path)) {
      add_written("results_zip", "zip", zip_path)
      zip_written <- TRUE
    } else if (inherits(zip_result, "error")) {
      zip_note <- conditionMessage(zip_result)
    } else {
      zip_note <- "Zip archive was requested but was not created."
    }
  }

  summary_out <- data.frame(
    OutputDir = output_dir,
    Prefix = prefix,
    FilesWritten = nrow(written_files),
    CsvWritten = sum(written_files$Format %in% "csv"),
    HtmlWritten = sum(written_files$Format %in% "html"),
    RdsWritten = sum(written_files$Format %in% "rds"),
    ReplayWritten = sum(written_files$Format %in% "R"),
    PlotWritten = sum(written_files$Format %in% "png"),
    PlotErrors = nrow(plot_errors),
    ZipWritten = isTRUE(zip_written),
    stringsAsFactors = FALSE
  )

  out <- list(
    output_dir = output_dir,
    prefix = prefix,
    include = include,
    summary = summary_out,
    written_files = written_files,
    plot_errors = plot_errors,
    zip_path = if (isTRUE(zip_written)) normalizePath(zip_path, winslash = "/", mustWork = FALSE) else NULL,
    zip_note = zip_note
  )
  class(out) <- "mfrm_results_export"
  out
}

#' Build comprehensive first-screen MFRM results
#'
#' @param fit Output from [fit_mfrm()] or [run_mfrm_facets()]. A standard
#'   long-format `data.frame` is also accepted when person and score columns can
#'   be inferred unambiguously from common names such as `Person` and `Score`;
#'   all remaining columns are treated as facets.
#' @param include Result sections or purpose presets to include. Purpose
#'   presets are `"standard"`, `"publication"`, `"validation"`, `"facets"`,
#'   `"bias"`, `"misfit_review"`, `"linking"`, `"network"`,
#'   `"gpcm_review"`, and `"all"`.
#'   Section names include `"fit"`, `"diagnostics"`, `"tables"`,
#'   `"precision"`, `"reporting"`, `"categories"`, `"plots"`,
#'   `"facets_fit"`, `"bias"`, `"misfit"`, `"linking"`, `"network"`,
#'   and `"apa"`.
#' @param output Return format: `"object"` for an `mfrm_results` object,
#'   `"summary"` for its compact summary, `"tables"` for a named list of
#'   available data frames, or `"html"` for a temporary HTML report.
#'
#' @details
#' `mfrm_results()` is a high-level result object. It does not introduce a new
#' estimator or a new validity rule. It fits only when `fit` is a data frame,
#' computes diagnostics automatically when needed, and collects output from
#' existing helpers such as [diagnose_mfrm()],
#' [fit_measures_table()], [precision_review_report()], and
#' [reporting_checklist()]. Sections that are unsupported for a particular fit
#' are retained in the `status` table as `not_available` rather than stopping
#' the whole results workflow. The returned object also carries
#' `next_actions` and `input$reproducible_code` so users can move from the
#' comprehensive first screen to explicit reporting or replay code.
#'
#' @section Include presets:
#' - `"standard"`: fit, diagnostics, tables, precision, reporting, categories,
#'   and plot routes
#' - `"publication"`: standard sections plus APA output assembly
#' - `"validation"`: standard sections plus FACETS-fit/df-sensitivity review
#' - `"facets"`: fit, diagnostics, tables, categories, plots, and FACETS-fit
#'   review for FACETS-facing migration work
#' - `"bias"` / `"bias_review"`: standard sections plus facet-level bias-screen
#'   guidance; interaction bias still requires explicit facet-pair selection
#' - `"misfit"` / `"misfit_review"`: standard sections plus unexpected-response,
#'   displacement, and pathway-map case-review surfaces
#' - `"linking"` / `"anchors"`: standard sections plus anchor-readiness and
#'   operational linking-review surfaces from the fitted object's stored
#'   anchor review; drift and screened-chain review still require multiple
#'   fitted forms or waves
#' - `"network"`: standard sections plus network/connectivity review
#' - `"response_time"`: descriptive response-time QC review when timing
#'   metadata are supplied through `response_time` / `response_time_data`
#' - `"gpcm_review"`: standard sections with bounded-`GPCM` caveats retained
#'   in the collected summaries and reports
#' - `"all"`: standard sections plus FACETS-fit, network, APA, and
#'   response-time sections
#'
#' @section Response-time metadata:
#' Response-time review is opt-in and descriptive. It does not change fitted
#' MFRM estimates, fit a joint speed-accuracy model, or create automatic
#' exclusion rules. Use `include = "response_time"` together with
#' `response_time = "ResponseTime"`. When `fit` is an already fitted object,
#' also supply `response_time_data = original_data` because fitted objects keep
#' only the measurement columns needed for estimation.
#'
#' @section What to inspect first:
#' Start with `summary(res)`. The most useful fields are:
#' - `overview`: input mode, model, method, table count, and plot-route count
#' - `triage`: first-screen signals ordered by unavailable/review/info/ok
#' - `status`: which sections were available, skipped, or unsupported
#' - `plot_map`: the supported `plot(res, type = ...)` routes for this object
#' - `next_actions`: recommended follow-up calls
#' - `reproducible_code`: replay scaffold for the first-screen route
#'
#' @section Data-frame input:
#' Direct data-frame input is intentionally conservative. It is intended for
#' standard columns such as `Person`, `Score`, `Rater`, and `Criterion`. For
#' research scripts, use [fit_mfrm()] or [run_mfrm_facets()] explicitly when
#' column roles, model, method, anchors, or missing-data rules need to be
#' documented. Use [mfrm_results_interactive()] only when you want an opt-in
#' column-selection wizard in an interactive session.
#'
#' @section Visualization and HTML:
#' `plot(res)` routes to a FACETS-style model-level visual bundle by default.
#' Other routes include `plot(res, type = "wright")`, `"pathway"`, `"qc"`,
#' `"category"`, `"anchors"`, and `"tables"`. `output = "html"` writes a
#' lightweight temporary HTML file;
#' use [launch_mfrmr_viewer()] when you want an optional local Shiny reader
#' for an already-created `mfrm_results` object. Use
#' [export_mfrm_results()] for a lightweight download of the comprehensive
#' results object, or [export_mfrm_bundle()] when a fit-centered durable
#' analysis archive is needed.
#'
#' @section Typical workflow:
#' 1. Fit explicitly with [fit_mfrm()] in scripts and manuscripts.
#' 2. Call `res <- mfrm_results(fit)`.
#' 3. Read `summary(res)$triage`, `summary(res)$status`,
#'    `summary(res)$plot_map`, and `summary(res)$next_actions`.
#' 4. Use `plot(res, type = "qc")` for the first visual screen.
#' 5. Optionally inspect the same result with [launch_mfrmr_viewer()] in an
#'    interactive session.
#' 6. Use [build_summary_table_bundle()] or the helper named in
#'    `summary(res)$next_actions` for report-specific follow-up.
#'
#' @param response_time Optional response-time column name. When `NULL` and
#'   `include` contains `"response_time"`, conservative column names such as
#'   `ResponseTime`, `response_time`, or `RT` are detected when available.
#' @param response_time_data Optional original long-format data containing the
#'   timing column. Required for already fitted objects unless the timing
#'   column is still present in `fit$prep$data`.
#' @param response_time_facets Optional facet columns for response-time
#'   summaries. Defaults to the fitted model's source facet columns when
#'   available.
#' @param response_time_score Optional score column for response-time
#'   summaries. Defaults to the fitted model's source score column when
#'   available.
#'
#' @return Depending on `output`, an `mfrm_results` object, a
#'   `summary.mfrm_results` object, a named table list, or an
#'   `mfrm_results_html` object.
#' @seealso [fit_mfrm()], [run_mfrm_facets()], [diagnose_mfrm()],
#'   [reporting_checklist()], [build_summary_table_bundle()],
#'   [export_mfrm_results()],
#'   [launch_mfrmr_viewer()], [mfrmr_output_guide()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' toy_small <- toy[toy$Person %in% unique(toy$Person)[1:8], , drop = FALSE]
#'
#' # JML keeps the help example fast; use the recommended workflow settings
#' # for final analyses.
#' fit <- fit_mfrm(toy_small, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' res <- mfrm_results(fit)
#'
#' sx <- summary(res)
#' sx$overview
#' sx$triage
#' sx$plot_map
#' sx$next_actions
#' mfrm_results(fit, include = "validation", output = "summary")$status
#'
#' plot(res, type = "qc", draw = FALSE)
#'
#' # Direct data-frame input is available for conservative exploratory use
#' # when Person and Score columns are unambiguous.
#' mfrm_results(
#'   toy_small,
#'   include = c("fit", "diagnostics"),
#'   output = "summary"
#' )$mapping
#' }
#' @export
mfrm_results <- function(fit,
                         include = "standard",
                         response_time = NULL,
                         response_time_data = NULL,
                         response_time_facets = NULL,
                         response_time_score = NULL,
                         output = c("object", "summary", "tables", "html")) {
  output <- match.arg(tolower(as.character(output[1])), c("object", "summary", "tables", "html"))
  include <- mfrm_results_resolve_include(include)
  ctx <- mfrm_results_resolve_input(fit)
  rt_requested <- "response_time" %in% include ||
    !is.null(response_time) ||
    !is.null(response_time_data)
  if (isTRUE(rt_requested) && !"response_time" %in% include) {
    requested <- attr(include, "requested", exact = TRUE)
    presets <- attr(include, "presets", exact = TRUE)
    include <- unique(c(include, "response_time"))
    attr(include, "requested") <- unique(c(requested %||% character(0), "response_time"))
    attr(include, "presets") <- presets %||% character(0)
  }
  if (isTRUE(rt_requested) && is.null(response_time)) {
    response_time <- mfrm_results_maybe_response_time_column(response_time_data %||% ctx$source_data)
  }
  ctx$response_time <- list(
    requested = isTRUE(rt_requested),
    time = response_time,
    data = response_time_data,
    data_supplied = !is.null(response_time_data),
    facets = response_time_facets,
    score = response_time_score
  )
  out <- mfrm_results_build(ctx, include = include)

  switch(
    output,
    object = out,
    summary = summary(out),
    tables = out$tables,
    html = mfrm_results_html(out)
  )
}

mfrm_results_menu_one <- function(cols, role) {
  idx <- utils::menu(cols, title = paste0("Choose the ", role, " column"))
  if (!is.finite(idx) || idx < 1L) {
    stop("No ", role, " column was selected.", call. = FALSE)
  }
  cols[[idx]]
}

mfrm_results_menu_many <- function(cols, role) {
  if (length(cols) == 0L) {
    stop("No candidate ", role, " columns are available.", call. = FALSE)
  }
  cat("\nChoose ", role, " columns by number, separated by commas:\n", sep = "")
  for (i in seq_along(cols)) {
    cat(sprintf("  %d: %s\n", i, cols[[i]]))
  }
  ans <- readline("Columns: ")
  idx <- suppressWarnings(as.integer(strsplit(ans, ",", fixed = TRUE)[[1]]))
  idx <- idx[is.finite(idx) & idx >= 1L & idx <= length(cols)]
  idx <- unique(idx)
  if (length(idx) == 0L) {
    stop("No ", role, " columns were selected.", call. = FALSE)
  }
  cols[idx]
}

mfrm_results_render_code <- function(person, facets, score, weight, include, output,
                                     response_time_lines = character(0)) {
  facet_expr <- paste(sprintf("%s", deparse(facets)), collapse = "")
  if (!grepl("^c\\(", facet_expr)) {
    facet_expr <- paste0("c(", paste(vapply(facets, function(x) deparse(x), character(1)), collapse = ", "), ")")
  }
  include_expr <- mfrm_results_include_expr(include)
  output_expr <- mfrm_results_deparse_one(output)
  lines <- c(
    "fit <- fit_mfrm(",
    "  data = data,",
    paste0("  person = ", deparse(person), ","),
    paste0("  facets = ", facet_expr, ","),
    paste0("  score = ", deparse(score), ","),
    if (!is.null(weight)) paste0("  weight = ", deparse(weight), ",") else NULL,
    "  model = \"RSM\",",
    "  method = \"JML\"",
    ")",
    if (length(response_time_lines) > 0L) {
      c(
        "res <- mfrm_results(",
        "  fit,",
        paste0("  include = ", include_expr, ","),
        response_time_lines,
        paste0("  output = ", output_expr),
        ")"
      )
    } else {
      paste0(
        "res <- mfrm_results(fit, include = ", include_expr,
        ", output = ", output_expr, ")"
      )
    }
  )
  paste(lines, collapse = "\n")
}

#' Interactively choose data-frame columns before calling `mfrm_results()`
#'
#' @param data A long-format data frame.
#' @param include Passed to [mfrm_results()].
#' @param output Passed to [mfrm_results()].
#'
#' @details
#' This helper is deliberately opt-in and stops in non-interactive sessions.
#' It asks the user to choose the person, score, optional weight, and facet
#' columns, prints reproducible code for the selected roles, then fits the
#' default legacy-compatible `RSM`/`JML` route before calling
#' [mfrm_results()]. Use explicit [fit_mfrm()] calls in scripts, Quarto
#' documents, tests, and reproducible analyses.
#'
#' @section Why this helper is opt-in:
#' Interactive prompts are useful at the console but are unsafe defaults for
#' reproducible analysis, package checks, batch scripts, and manuscripts. The
#' helper therefore prints replay code and leaves the scripted route explicit.
#'
#' @return The selected [mfrm_results()] output.
#' @seealso [mfrm_results()], [fit_mfrm()], [run_mfrm_facets()]
#' @examples
#' if (interactive()) {
#'   toy <- load_mfrmr_data("example_core")
#'   res <- mfrm_results_interactive(toy)
#' }
#' @export
mfrm_results_interactive <- function(data,
                                     include = "standard",
                                     output = c("object", "summary", "tables", "html")) {
  if (!interactive()) {
    stop("`mfrm_results_interactive()` can only be used in an interactive session.", call. = FALSE)
  }
  output <- match.arg(tolower(as.character(output[1])), c("object", "summary", "tables", "html"))
  include <- mfrm_results_resolve_include(include)
  dat <- normalize_facets_mode_data(data)
  cols <- names(dat)
  person <- mfrm_results_menu_one(cols, "person")
  score <- mfrm_results_menu_one(setdiff(cols, person), "score")
  remaining <- setdiff(cols, c(person, score))
  weight <- NULL
  if (length(remaining) > 0L) {
    use_weight <- utils::menu(c("No weight column", remaining), title = "Choose a weight column, if any")
    if (is.finite(use_weight) && use_weight > 1L) {
      weight <- remaining[[use_weight - 1L]]
    }
  }
  facets <- mfrm_results_menu_many(setdiff(cols, c(person, score, weight)), "facet")
  code <- mfrm_results_render_code(person, facets, score, weight, include, output)
  message("Reproducible code for this selection:\n", code)
  run <- run_mfrm_facets(
    data = dat,
    person = person,
    facets = facets,
    score = score,
    weight = weight
  )
  out <- mfrm_results(run, include = include, output = output)
  attr(out, "mfrm_results_code") <- code
  out
}

#' @export
summary.mfrm_results <- function(object, digits = 3, top_n = 10, ...) {
  if (!inherits(object, "mfrm_results")) {
    stop("`object` must be an mfrm_results object.", call. = FALSE)
  }
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))
  fit <- object$fit
  fit_summary <- as.data.frame(fit$summary %||% data.frame(), stringsAsFactors = FALSE)
  ov <- if (nrow(fit_summary) > 0L) fit_summary[1, , drop = FALSE] else data.frame()
  overview <- data.frame(
    InputMode = as.character(object$input$mode %||% ""),
    Model = as.character(ov$Model[1] %||% fit$config$model %||% ""),
    Method = as.character(ov$Method[1] %||% fit$config$method %||% ""),
    N = suppressWarnings(as.integer(ov$N[1] %||% fit$prep$n_obs %||% NA_integer_)),
    Persons = suppressWarnings(as.integer(ov$Persons[1] %||% NA_integer_)),
    Facets = suppressWarnings(as.integer(ov$Facets[1] %||% length(fit$config$facet_names %||% character(0)))),
    Categories = suppressWarnings(as.integer(ov$Categories[1] %||% NA_integer_)),
    Components = length(object$components %||% list()),
    Tables = length(object$tables %||% list()),
    PlotRoutes = sum(object$plot_map$Available %in% TRUE, na.rm = TRUE),
    NotAvailable = sum(object$status$Status %in% "not_available", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  status <- as.data.frame(object$status %||% data.frame(), stringsAsFactors = FALSE)
  component_names <- names(object$components %||% list())
  component_index <- if (length(component_names) == 0L) {
    data.frame()
  } else {
    data.frame(
      Component = component_names,
      Class = vapply(object$components, mfrm_results_component_class, character(1)),
      stringsAsFactors = FALSE
    )
  }
  table_index <- as.data.frame(object$table_index %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(table_index) > top_n) {
    table_index <- utils::head(table_index, n = top_n)
  }
  plot_map <- as.data.frame(object$plot_map %||% data.frame(), stringsAsFactors = FALSE)
  triage <- as.data.frame(object$triage %||% data.frame(), stringsAsFactors = FALSE)
  next_actions <- as.data.frame(object$next_actions %||% data.frame(), stringsAsFactors = FALSE)
  mapping <- mfrm_results_mapping_table(object$input$mapping %||% NULL)
  reproducible_code <- mfrm_results_code_table(object$input$reproducible_code %||% "")

  out <- list(
    overview = round_numeric_df(overview, digits = digits),
    status = status,
    component_index = component_index,
    table_index = table_index,
    plot_map = plot_map,
    triage = triage,
    next_actions = next_actions,
    mapping = mapping,
    reproducible_code = reproducible_code,
    notes = object$notes %||% character(0),
    digits = digits
  )
  class(out) <- "summary.mfrm_results"
  out
}

#' @export
print.summary.mfrm_results <- function(x, ...) {
  cat("mfrmr Results Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    print(as.data.frame(x$overview), row.names = FALSE)
  }
  if (!is.null(x$status) && nrow(x$status) > 0L) {
    cat("\nSection status\n")
    print(as.data.frame(x$status), row.names = FALSE)
  }
  if (!is.null(x$triage) && nrow(x$triage) > 0L) {
    cat("\nTriage\n")
    print(as.data.frame(x$triage), row.names = FALSE)
  }
  if (!is.null(x$plot_map) && nrow(x$plot_map) > 0L) {
    cat("\nPlot routes\n")
    print(as.data.frame(x$plot_map), row.names = FALSE)
  }
  if (!is.null(x$next_actions) && nrow(x$next_actions) > 0L) {
    cat("\nNext actions\n")
    print(as.data.frame(x$next_actions), row.names = FALSE)
  }
  if (!is.null(x$notes) && length(x$notes) > 0L) {
    cat("\nNotes\n")
    for (line in x$notes) cat(" - ", line, "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.mfrm_results <- function(x, ...) {
  print(summary(x, ...), ...)
  invisible(x)
}

#' @export
print.mfrm_results_export <- function(x, ...) {
  cat("mfrmr Results Export\n")
  if (!is.null(x$summary)) {
    print(as.data.frame(x$summary), row.names = FALSE)
  }
  if (nrow(as.data.frame(x$written_files %||% data.frame())) > 0L) {
    cat("\nWritten files\n")
    print(utils::head(as.data.frame(x$written_files), 10L), row.names = FALSE)
    if (nrow(as.data.frame(x$written_files)) > 10L) {
      cat("... ", nrow(as.data.frame(x$written_files)) - 10L, " more file(s)\n", sep = "")
    }
  }
  if (length(x$zip_note %||% character(0)) > 0L) {
    cat("\nZip note: ", x$zip_note, "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.mfrm_report <- function(x, ...) {
  cat(x$title, "\n", sep = "")
  cat("  Style: ", x$style, "\n", sep = "")
  if (length(x$source_include %||% character(0)) > 0L) {
    cat("  Source include: ", paste(x$source_include, collapse = ", "), "\n", sep = "")
  }
  cat("  Read order: summary(report) -> report$first_screen -> report$report_index -> report$template_index\n", sep = "")
  first_screen <- as.data.frame(x$first_screen %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(first_screen) > 0L) {
    cat("\nFirst screen\n")
    print(utils::head(first_screen[, intersect(c(
      "Area", "Status", "Readiness", "MainIssue", "NextAction", "PrimaryRoute"
    ), names(first_screen)), drop = FALSE], 8L), row.names = FALSE)
  }
  report_index <- as.data.frame(x$report_index %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(report_index) > 0L) {
    cat("\nReport index\n")
    print(utils::head(report_index[, intersect(c(
      "Area", "EvidenceStatus", "Readiness", "ReviewSignalCount",
      "PrimaryTable", "TemplateTable", "PlotRoute"
    ), names(report_index)), drop = FALSE], 8L), row.names = FALSE)
  }
  sx <- summary(x)
  actions <- as.data.frame(sx$immediate_actions %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(actions) > 0L) {
    cat("\nImmediate actions\n")
    print(utils::head(actions[, intersect(c(
      "Area", "Status", "MainIssue", "NextAction", "PrimaryRoute", "TemplateRoute"
    ), names(actions)), drop = FALSE], 8L), row.names = FALSE)
  }
  optional <- as.data.frame(sx$optional_sections %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(optional) > 0L) {
    cat("\nOptional sections not requested\n")
    print(utils::head(optional[, intersect(c(
      "Area", "Status", "MainIssue", "NextAction", "PrimaryRoute"
    ), names(optional)), drop = FALSE], 8L), row.names = FALSE)
  }
  cat("\nDetailed tables are available in report$tables.\n", sep = "")
  cat("Use report$report_index$PrimaryTable and report$report_index$TemplateTable to choose the next table.\n", sep = "")
  invisible(x)
}

#' @export
print.summary.mfrm_report <- function(x, ...) {
  cat("mfrmr Report Summary\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    print(as.data.frame(x$overview), row.names = FALSE)
  }
  if (!is.null(x$first_screen) && nrow(x$first_screen) > 0L) {
    cat("\nFirst screen\n")
    print(as.data.frame(x$first_screen), row.names = FALSE)
  }
  if (!is.null(x$immediate_actions) && nrow(x$immediate_actions) > 0L) {
    cat("\nImmediate actions\n")
    print(as.data.frame(x$immediate_actions), row.names = FALSE)
  }
  if (!is.null(x$optional_sections) && nrow(x$optional_sections) > 0L) {
    cat("\nOptional sections not requested\n")
    print(as.data.frame(x$optional_sections), row.names = FALSE)
  }
  if (!is.null(x$claim_readiness) && nrow(x$claim_readiness) > 0L) {
    cat("\nClaim readiness\n")
    print(as.data.frame(x$claim_readiness), row.names = FALSE)
  }
  if (!is.null(x$report_gaps) && nrow(x$report_gaps) > 0L) {
    cat("\nReport gaps\n")
    print(as.data.frame(x$report_gaps), row.names = FALSE)
  }
  if (!is.null(x$boundary_index) && nrow(x$boundary_index) > 0L) {
    cat("\nBoundary index\n")
    print(as.data.frame(x$boundary_index), row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.mfrm_report_html <- function(x, ...) {
  cat("mfrmr Report HTML\n")
  cat("  Path: ", x$path, "\n", sep = "")
  invisible(x)
}

#' @export
plot.mfrm_results <- function(x,
                              y = NULL,
                              type = c("fit", "wright", "pathway", "qc", "category", "anchors", "response_time", "tables"),
                              ...) {
  if (!inherits(x, "mfrm_results")) {
    stop("`x` must be an mfrm_results object.", call. = FALSE)
  }
  type <- match.arg(type)
  available <- as.data.frame(x$plot_map %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(available) > 0L && "Type" %in% names(available)) {
    row <- available[available$Type %in% type, , drop = FALSE]
    if (nrow(row) > 0L && !isTRUE(row$Available[1])) {
      stop("Plot route `", type, "` is not available for this mfrm_results object.", call. = FALSE)
    }
  }
  if (identical(type, "fit")) {
    dots <- list(...)
    if (is.null(dots$type)) dots$type <- "bundle"
    return(do.call(plot, c(list(x$fit), dots)))
  }
  if (identical(type, "wright")) {
    dots <- list(...)
    dots$type <- "wright"
    return(do.call(plot, c(list(x$fit), dots)))
  }
  if (identical(type, "pathway")) {
    dots <- list(...)
    dots$type <- "pathway"
    return(do.call(plot, c(list(x$fit), dots)))
  }
  if (identical(type, "qc")) {
    if (!inherits(x$diagnostics, "mfrm_diagnostics")) {
      stop("QC plotting requires available diagnostics.", call. = FALSE)
    }
    return(plot_qc_dashboard(fit = x$fit, diagnostics = x$diagnostics, ...))
  }
  if (identical(type, "category")) {
    rating <- x$components$rating_scale %||%
      rating_scale_table(x$fit, diagnostics = x$diagnostics)
    return(plot(rating, ...))
  }
  if (identical(type, "anchors")) {
    anchor_review <- x$fit$config$anchor_review %||% NULL
    if (!inherits(anchor_review, "mfrm_anchor_review")) {
      stop("Anchor plotting requires stored anchor-review metadata.", call. = FALSE)
    }
    return(plot(anchor_review, ...))
  }
  if (identical(type, "response_time")) {
    rt <- x$components$response_time_review %||% NULL
    if (!inherits(rt, "mfrm_response_time_review")) {
      stop("Response-time plotting requires include = \"response_time\" with timing metadata.", call. = FALSE)
    }
    return(plot_response_time_review(rt, ...))
  }
  bundle <- build_summary_table_bundle(summary(x), include_empty = FALSE)
  plot(bundle, ...)
}

#' @export
print.mfrm_results_html <- function(x, ...) {
  cat("mfrmr Results HTML\n")
  cat("  Path: ", x$path, "\n", sep = "")
  invisible(x)
}
