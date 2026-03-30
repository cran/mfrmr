resolve_mfrm_export_context <- function(x,
                                        diagnostics = NULL,
                                        residual_pca = c("none", "overall", "facet", "both")) {
  residual_pca <- match.arg(tolower(as.character(residual_pca[1] %||% "none")),
                            c("none", "overall", "facet", "both"))
  diagnostics_supplied <- !is.null(diagnostics)

  run_obj <- NULL
  mapping <- NULL
  run_info <- data.frame()

  if (inherits(x, "mfrm_facets_run")) {
    run_obj <- x
    fit <- x$fit
    mapping <- x$mapping %||% NULL
    run_info <- as.data.frame(x$run_info %||% data.frame(), stringsAsFactors = FALSE)
  } else if (inherits(x, "mfrm_fit")) {
    fit <- x
  } else {
    stop("`fit` must be an mfrm_fit or mfrm_facets_run object.", call. = FALSE)
  }

  if (!inherits(fit, "mfrm_fit")) {
    stop("Resolved `fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }

  if (is.null(diagnostics) && !is.null(run_obj$diagnostics)) {
    diagnostics <- run_obj$diagnostics
    if (residual_pca %in% c("overall", "both") &&
        is.null(diagnostics$residual_pca_overall)) {
      diagnostics <- NULL
    }
    if (!is.null(diagnostics) &&
        residual_pca %in% c("facet", "both") &&
        (is.null(diagnostics$residual_pca_by_facet) ||
         length(diagnostics$residual_pca_by_facet) == 0)) {
      diagnostics <- NULL
    }
  }

  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = residual_pca)
  }
  if (!is.list(diagnostics) || is.null(diagnostics$obs)) {
    stop("`diagnostics` must be output from diagnose_mfrm().", call. = FALSE)
  }

  list(
    fit = fit,
    diagnostics = diagnostics,
    diagnostics_supplied = diagnostics_supplied,
    input_mode = if (is.null(run_obj)) "fit" else "facets_run",
    run = run_obj,
    mapping = mapping,
    run_info = run_info
  )
}

infer_export_residual_pca_mode <- function(diagnostics) {
  has_overall <- !is.null(diagnostics$residual_pca_overall)
  has_facet <- !is.null(diagnostics$residual_pca_by_facet) &&
    length(diagnostics$residual_pca_by_facet) > 0
  if (has_overall && has_facet) return("both")
  if (has_overall) return("overall")
  if (has_facet) return("facet")
  "none"
}

render_r_object_literal <- function(x) {
  paste(utils::capture.output(dput(x)), collapse = "\n")
}

render_named_text_map <- function(x, title = NULL) {
  parts <- character(0)
  if (!is.null(title) && nzchar(as.character(title[1]))) {
    parts <- c(parts, as.character(title[1]))
  }
  if (length(x) == 0) return(paste(c(parts, "No sections available."), collapse = "\n"))
  for (nm in names(x)) {
    vals <- as.character(x[[nm]] %||% character(0))
    vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
    if (length(vals) == 0) next
    parts <- c(parts, "", paste0("[", nm, "]"), paste0("- ", vals))
  }
  paste(parts, collapse = "\n")
}

export_extract_bias_pairs <- function(bias_results) {
  add_pair <- function(out, pair) {
    pair <- unique(as.character(pair))
    pair <- pair[!is.na(pair) & nzchar(pair)]
    if (length(pair) < 2L) return(out)
    key <- paste(pair, collapse = " x ")
    if (!key %in% names(out)) out[[key]] <- pair
    out
  }

  out <- list()
  if (inherits(bias_results, "mfrm_bias")) {
    return(unname(add_pair(out, bias_results$interaction_facets %||% c(bias_results$facet_a, bias_results$facet_b))))
  }

  if (inherits(bias_results, "mfrm_bias_collection")) {
    bias_results <- bias_results$by_pair %||% list()
  }

  if (is.list(bias_results) && !is.data.frame(bias_results)) {
    for (nm in names(bias_results)) {
      obj <- bias_results[[nm]]
      if (inherits(obj, "mfrm_bias")) {
        out <- add_pair(out, obj$interaction_facets %||% c(obj$facet_a, obj$facet_b))
      } else {
        pair <- strsplit(as.character(nm[1] %||% ""), "\\s+x\\s+")[[1]]
        out <- add_pair(out, pair)
      }
    }
  }

  unname(out)
}

#' Build a reproducibility manifest for an MFRM analysis
#'
#' @param fit Output from [fit_mfrm()] or [run_mfrm_facets()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. When `NULL`,
#'   diagnostics are computed with `residual_pca = "none"`.
#' @param bias_results Optional output from [estimate_bias()] or a named list of
#'   bias bundles.
#' @param population_prediction Optional output from
#'   [predict_mfrm_population()].
#' @param unit_prediction Optional output from [predict_mfrm_units()].
#' @param plausible_values Optional output from [sample_mfrm_plausible_values()].
#' @param include_person_anchors If `TRUE`, include person measures in the
#'   exported anchor table.
#'
#' @details
#' This helper captures the package-native equivalent of the Streamlit app's
#' configuration export. It summarizes analysis settings, source columns,
#' anchoring information, and which downstream outputs are currently available.
#'
#' @section When to use this:
#' Use `build_mfrm_manifest()` when you want a compact, machine-readable record
#' of how an analysis was run. Compared with related helpers:
#' - [export_mfrm()] writes analysis tables only.
#' - `build_mfrm_manifest()` records settings and available outputs.
#' - [build_mfrm_replay_script()] creates an executable R script.
#' - [export_mfrm_bundle()] writes a shareable folder of files.
#'
#' @section Output:
#' The returned bundle has class `mfrm_manifest` and includes:
#' - `summary`: one-row analysis overview
#' - `environment`: package/R/platform metadata
#' - `model_settings`: key-value model settings table
#' - `source_columns`: key-value data-column table
#' - `estimation_control`: key-value optimizer settings table
#' - `anchor_summary`: facet-level anchor summary
#' - `anchors`: machine-readable anchor table
#' - `available_outputs`: availability table for diagnostics/bias/PCA/prediction
#'   outputs
#' - `settings`: manifest build settings
#'
#' @section Interpreting output:
#' The `summary` table is the quickest place to confirm that you are looking at
#' the intended analysis. The `model_settings`, `source_columns`, and
#' `estimation_control` tables are designed for audit trails and method write-up.
#' The `available_outputs` table is especially useful before building bundles,
#' because it tells you whether residual PCA, anchors, bias results, or
#' prediction-side artifacts are already available. A practical reading order is
#' `summary` first, `available_outputs` second, and `anchors` last when
#' reproducibility depends on fixed constraints.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()] or [run_mfrm_facets()].
#' 2. Compute diagnostics once with [diagnose_mfrm()] if you want explicit
#'    control over residual PCA.
#' 3. Build a manifest and inspect `summary` plus `available_outputs`.
#' 4. If you need files on disk, pass the same objects to
#'    [export_mfrm_bundle()].
#'
#' @return A named list with class `mfrm_manifest`.
#' @seealso [export_mfrm_bundle()], [build_mfrm_replay_script()],
#'   [make_anchor_table()], [reporting_checklist()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' manifest <- build_mfrm_manifest(fit, diagnostics = diag)
#' manifest$summary[, c("Model", "Method", "Observations", "Facets")]
#' manifest$available_outputs[, c("Component", "Available")]
#' @export
build_mfrm_manifest <- function(fit,
                                diagnostics = NULL,
                                bias_results = NULL,
                                population_prediction = NULL,
                                unit_prediction = NULL,
                                plausible_values = NULL,
                                include_person_anchors = FALSE) {
  ctx <- resolve_mfrm_export_context(
    x = fit,
    diagnostics = diagnostics,
    residual_pca = "none"
  )
  fit <- ctx$fit
  diagnostics <- ctx$diagnostics
  diagnostics_supplied <- ctx$diagnostics_supplied

  bias_inputs <- export_normalize_bias_inputs(bias_results)
  population_prediction <- export_validate_optional_object(
    population_prediction,
    "mfrm_population_prediction",
    "population_prediction"
  )
  unit_prediction <- export_validate_optional_object(
    unit_prediction,
    "mfrm_unit_prediction",
    "unit_prediction"
  )
  plausible_values <- export_validate_optional_object(
    plausible_values,
    "mfrm_plausible_values",
    "plausible_values"
  )
  anchor_tbl <- make_anchor_table(
    fit = fit,
    include_person = isTRUE(include_person_anchors)
  )

  cfg <- fit$config %||% list()
  prep <- fit$prep %||% list()
  est_ctl <- cfg$estimation_control %||% list()

  summary_tbl <- data.frame(
    Model = as.character(cfg$model %||% NA_character_),
    Method = as.character(cfg$method_input %||% cfg$method %||% NA_character_),
    MethodUsed = as.character(cfg$method %||% NA_character_),
    Observations = as.integer(prep$n_obs %||% fit$summary$N[1] %||% NA_integer_),
    Persons = as.integer(prep$n_person %||% cfg$n_person %||% fit$summary$Persons[1] %||% NA_integer_),
    Facets = length(as.character(cfg$facet_names %||% character(0))),
    Categories = as.integer(cfg$n_cat %||% NA_integer_),
    BiasBundles = length(bias_inputs),
    HasResidualPCA = export_has_residual_pca(diagnostics),
    Converged = isTRUE(fit$summary$Converged %||% FALSE),
    stringsAsFactors = FALSE
  )

  environment_tbl <- data.frame(
    Package = "mfrmr",
    PackageVersion = as.character(utils::packageVersion("mfrmr")),
    RVersion = as.character(getRversion()),
    Platform = R.version$platform,
    Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    stringsAsFactors = FALSE
  )

  model_settings <- dashboard_settings_table(list(
    model = as.character(cfg$model %||% NA_character_),
    method = as.character(cfg$method_input %||% cfg$method %||% NA_character_),
    method_used = as.character(cfg$method %||% NA_character_),
    facet_names = paste(as.character(cfg$facet_names %||% character(0)), collapse = ", "),
    noncenter_facet = as.character(cfg$noncenter_facet %||% ""),
    step_facet = as.character(cfg$step_facet %||% ""),
    dummy_facets = paste(as.character(cfg$dummy_facets %||% character(0)), collapse = ", "),
    n_categories = as.character(cfg$n_cat %||% NA_character_)
  ))

  source_columns <- dashboard_settings_table(list(
    person = as.character(cfg$source_columns$person %||% "Person"),
    facets = paste(as.character(cfg$source_columns$facets %||% cfg$facet_names %||% character(0)), collapse = ", "),
    score = as.character(cfg$source_columns$score %||% "Score"),
    weight = as.character(cfg$source_columns$weight %||% "")
  ))

  estimation_control <- dashboard_settings_table(list(
    maxit = as.integer(est_ctl$maxit %||% NA_integer_),
    reltol = as.numeric(est_ctl$reltol %||% NA_real_),
    quad_points = as.integer(est_ctl$quad_points %||% NA_integer_)
  ))

  anchor_summary <- as.data.frame(cfg$anchor_summary %||% data.frame(), stringsAsFactors = FALSE)
  available_outputs <- export_available_outputs_table(
    diagnostics = diagnostics,
    bias_inputs = bias_inputs,
    anchor_tbl = anchor_tbl,
    population_prediction = population_prediction,
    unit_prediction = unit_prediction,
    plausible_values = plausible_values
  )

  settings <- dashboard_settings_table(list(
    diagnostics_supplied = diagnostics_supplied,
    include_person_anchors = isTRUE(include_person_anchors),
    bias_collection = inherits(bias_results, "mfrm_bias_collection"),
    population_prediction = !is.null(population_prediction),
    unit_prediction = !is.null(unit_prediction),
    plausible_values = !is.null(plausible_values),
    input_mode = ctx$input_mode
  ))

  out <- list(
    summary = summary_tbl,
    environment = environment_tbl,
    model_settings = model_settings,
    source_columns = source_columns,
    estimation_control = estimation_control,
    anchor_summary = anchor_summary,
    anchors = as.data.frame(anchor_tbl, stringsAsFactors = FALSE),
    available_outputs = available_outputs,
    settings = settings
  )
  as_mfrm_bundle(out, "mfrm_manifest")
}

#' Build a package-native replay script for an MFRM analysis
#'
#' @param fit Output from [fit_mfrm()] or [run_mfrm_facets()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. When `NULL`,
#'   diagnostics are reused from `run_mfrm_facets()` when available, otherwise
#'   recomputed.
#' @param bias_results Optional output from [estimate_bias()] or a named list of
#'   bias bundles. When supplied, the generated script includes package-native
#'   bias estimation calls.
#' @param population_prediction Optional output from
#'   [predict_mfrm_population()] to recreate in the generated script.
#' @param unit_prediction Optional output from [predict_mfrm_units()] to
#'   recreate in the generated script.
#' @param plausible_values Optional output from
#'   [sample_mfrm_plausible_values()] to recreate in the generated script.
#' @param data_file Path to the analysis data file used in the generated script.
#' @param script_mode One of `"auto"`, `"fit"`, or `"facets"`. `"auto"` uses
#'   `run_mfrm_facets()` when the input object came from that workflow.
#' @param include_bundle If `TRUE`, append an [export_mfrm_bundle()] call to the
#'   generated script.
#' @param bundle_dir Output directory used when `include_bundle = TRUE`.
#' @param bundle_prefix Prefix used by the generated bundle exporter call.
#'
#' @details
#' This helper mirrors the Streamlit app's reproducible-download idea, but uses
#' `mfrmr`'s installed API rather than embedding a separate estimation engine.
#' The generated script assumes the user has the package installed and provides
#' a data file at `data_file`.
#'
#' Anchor and group-anchor constraints are embedded directly from the fitted
#' object's stored configuration, so the script can replay anchored analyses
#' without manual table reconstruction.
#'
#' @section When to use this:
#' Use `build_mfrm_replay_script()` when you want a package-native recipe that
#' another analyst can rerun later. Compared with related helpers:
#' - [build_mfrm_manifest()] records settings but does not run anything.
#' - `build_mfrm_replay_script()` produces executable R code.
#' - [export_mfrm_bundle()] can optionally write the replay script to disk.
#'
#' @section Interpreting output:
#' The returned object contains:
#' - `summary`: a one-row overview of the chosen replay mode and whether bundle
#'   export was included
#' - `script`: the generated R code as a single string
#' - `anchors` and `group_anchors`: the exact stored constraints that were
#'   embedded into the script
#'
#' If `ScriptMode` is `"facets"`, the script replays the higher-level
#' [run_mfrm_facets()] workflow. If it is `"fit"`, the script uses
#' [fit_mfrm()] directly.
#'
#' @section Mode guide:
#' - `"auto"` is the safest default and follows the structure of the supplied
#'   object.
#' - `"fit"` is useful when you want a minimal script centered on
#'   [fit_mfrm()].
#' - `"facets"` is useful when you want to preserve the higher-level
#'   [run_mfrm_facets()] workflow, including stored column mapping.
#'
#' @section Typical workflow:
#' 1. Finalize a fit and diagnostics object.
#' 2. Generate the replay script with the path you want users to read from.
#' 3. Write `replay$script` to disk, or let [export_mfrm_bundle()] do it for
#'    you.
#' 4. Rerun the script in a fresh R session to confirm reproducibility.
#'
#' @return A named list with class `mfrm_replay_script`.
#' @seealso [build_mfrm_manifest()], [export_mfrm_bundle()], [run_mfrm_facets()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 25)
#' replay <- build_mfrm_replay_script(fit, data_file = "your_data.csv")
#' replay$summary[, c("ScriptMode", "ResidualPCA", "BiasPairs")]
#' cat(substr(replay$script, 1, 120))
#' @export
build_mfrm_replay_script <- function(fit,
                                     diagnostics = NULL,
                                     bias_results = NULL,
                                     population_prediction = NULL,
                                     unit_prediction = NULL,
                                     plausible_values = NULL,
                                     data_file = "your_data.csv",
                                     script_mode = c("auto", "fit", "facets"),
                                     include_bundle = FALSE,
                                     bundle_dir = "analysis_bundle",
                                     bundle_prefix = "mfrmr_replay") {
  script_mode <- match.arg(tolower(as.character(script_mode[1] %||% "auto")),
                           c("auto", "fit", "facets"))
  ctx <- resolve_mfrm_export_context(
    x = fit,
    diagnostics = diagnostics,
    residual_pca = "none"
  )
  fit <- ctx$fit
  diagnostics <- ctx$diagnostics

  resolved_mode <- if (script_mode == "auto") {
    if (ctx$input_mode == "facets_run") "facets" else "fit"
  } else {
    script_mode
  }
  if (resolved_mode == "facets" && is.null(ctx$run)) {
    resolved_mode <- "fit"
  }

  cfg <- fit$config %||% list()
  src <- cfg$source_columns %||% list(
    person = "Person",
    facets = as.character(cfg$facet_names %||% character(0)),
    score = "Score",
    weight = NULL
  )
  est_ctl <- cfg$estimation_control %||% list()
  anchor_tables <- extract_anchor_tables(cfg)
  anchor_df <- as.data.frame(anchor_tables$anchors %||% data.frame(), stringsAsFactors = FALSE)
  group_anchor_df <- as.data.frame(anchor_tables$groups %||% data.frame(), stringsAsFactors = FALSE)
  residual_pca_mode <- infer_export_residual_pca_mode(diagnostics)
  bias_pairs <- export_extract_bias_pairs(bias_results)
  include_diagnostics <- TRUE
  population_prediction <- export_validate_optional_object(
    population_prediction,
    "mfrm_population_prediction",
    "population_prediction"
  )
  unit_prediction <- export_validate_optional_object(
    unit_prediction,
    "mfrm_unit_prediction",
    "unit_prediction"
  )
  plausible_values <- export_validate_optional_object(
    plausible_values,
    "mfrm_plausible_values",
    "plausible_values"
  )

  render_classed_literal <- function(x, class_name) {
    paste0(
      "structure(",
      render_r_object_literal(unclass(x)),
      ", class = ",
      render_r_object_literal(as.character(class_name[1])),
      ")"
    )
  }

  lines <- c(
    "#!/usr/bin/env Rscript",
    "# Generated by mfrmr::build_mfrm_replay_script()",
    paste0("# Model: ", as.character(cfg$model %||% NA_character_),
           " | Method: ", as.character(cfg$method_input %||% cfg$method %||% NA_character_),
           if (!identical(cfg$method_input %||% NULL, cfg$method %||% NULL)) {
             paste0(" | InternalMethod: ", as.character(cfg$method %||% NA_character_))
           } else {
             ""
           }),
    "",
    "library(mfrmr)",
    "",
    paste0("data <- utils::read.csv(", render_r_object_literal(as.character(data_file[1])), ", stringsAsFactors = FALSE)")
  )

  lines <- c(
    lines,
    "",
    "# Stored constraints from the fitted analysis",
    paste0("anchors <- ", if (nrow(anchor_df) > 0) render_r_object_literal(anchor_df) else "NULL"),
    paste0("group_anchors <- ", if (nrow(group_anchor_df) > 0) render_r_object_literal(group_anchor_df) else "NULL")
  )

  if (resolved_mode == "facets") {
    mapping <- ctx$mapping %||% src
    top_n_interactions <- 20L
    if (nrow(ctx$run_info) > 0 && all(c("key", "value") %in% names(ctx$run_info))) {
      idx <- which(ctx$run_info$key == "top_n_interactions")[1]
      if (is.finite(idx)) {
        parsed <- suppressWarnings(as.integer(ctx$run_info$value[idx]))
        if (is.finite(parsed)) top_n_interactions <- parsed
      }
    }

	  lines <- c(
	    lines,
	    "",
	    "# Legacy-compatible workflow",
	    "run <- run_mfrm_facets(",
	    "  data = data,",
	    paste0("  person = ", render_r_object_literal(as.character(mapping$person %||% src$person)), ","),
	    paste0("  facets = ", render_r_object_literal(as.character(mapping$facets %||% src$facets)), ","),
	    paste0("  score = ", render_r_object_literal(as.character(mapping$score %||% src$score)), ","),
	    paste0("  weight = ", if (!is.null(mapping$weight %||% src$weight)) render_r_object_literal(as.character(mapping$weight %||% src$weight)) else "NULL", ","),
      paste0("  rating_min = ", render_r_object_literal(as.integer(cfg$rating_min %||% NA_integer_)), ","),
      paste0("  rating_max = ", render_r_object_literal(as.integer(cfg$rating_max %||% NA_integer_)), ","),
	    paste0("  keep_original = ", render_r_object_literal(isTRUE(cfg$keep_original)), ","),
	    paste0("  model = ", render_r_object_literal(as.character(cfg$model %||% "RSM")), ","),
	    paste0("  method = ", render_r_object_literal(as.character(cfg$method_input %||% cfg$method %||% "JML")), ","),
	    paste0("  step_facet = ", if (!is.null(cfg$step_facet) && nzchar(as.character(cfg$step_facet))) render_r_object_literal(as.character(cfg$step_facet)) else "NULL", ","),
	    "  anchors = anchors,",
	    "  group_anchors = group_anchors,",
	    paste0("  noncenter_facet = ", render_r_object_literal(as.character(cfg$noncenter_facet %||% "Person")), ","),
      paste0("  dummy_facets = ", if (length(cfg$dummy_facets %||% character(0)) > 0) render_r_object_literal(as.character(cfg$dummy_facets)) else "NULL", ","),
      paste0("  positive_facets = ", if (length(cfg$positive_facets %||% character(0)) > 0) render_r_object_literal(as.character(cfg$positive_facets)) else "NULL", ","),
      paste0("  quad_points = ", render_r_object_literal(as.integer(est_ctl$quad_points %||% 15L)), ","),
      paste0("  maxit = ", render_r_object_literal(as.integer(est_ctl$maxit %||% 400L)), ","),
      paste0("  reltol = ", render_r_object_literal(as.numeric(est_ctl$reltol %||% 1e-6)), ","),
      paste0("  top_n_interactions = ", render_r_object_literal(as.integer(top_n_interactions))),
      ")",
      "fit <- run$fit",
      "diagnostics <- run$diagnostics"
    )
    if (residual_pca_mode != "none") {
      lines <- c(
        lines,
        paste0("diagnostics <- diagnose_mfrm(fit, residual_pca = ", render_r_object_literal(residual_pca_mode), ")")
      )
    }
  } else {
	    lines <- c(
	      lines,
	      "",
	      "# Fit the model",
	      "fit <- fit_mfrm(",
	      "  data = data,",
	      paste0("  person = ", render_r_object_literal(as.character(src$person %||% "Person")), ","),
	      paste0("  facets = ", render_r_object_literal(as.character(src$facets %||% character(0))), ","),
	      paste0("  score = ", render_r_object_literal(as.character(src$score %||% "Score")), ","),
	      paste0("  weight = ", if (!is.null(src$weight)) render_r_object_literal(as.character(src$weight)) else "NULL", ","),
        paste0("  rating_min = ", render_r_object_literal(as.integer(cfg$rating_min %||% NA_integer_)), ","),
        paste0("  rating_max = ", render_r_object_literal(as.integer(cfg$rating_max %||% NA_integer_)), ","),
	      paste0("  keep_original = ", render_r_object_literal(isTRUE(cfg$keep_original)), ","),
	      paste0("  model = ", render_r_object_literal(as.character(cfg$model %||% "RSM")), ","),
	      paste0("  method = ", render_r_object_literal(as.character(cfg$method_input %||% cfg$method %||% "JML")), ","),
	      paste0("  step_facet = ", if (!is.null(cfg$step_facet) && nzchar(as.character(cfg$step_facet))) render_r_object_literal(as.character(cfg$step_facet)) else "NULL", ","),
	      "  anchors = anchors,",
	      "  group_anchors = group_anchors,",
	      paste0("  noncenter_facet = ", render_r_object_literal(as.character(cfg$noncenter_facet %||% "Person")), ","),
      paste0("  dummy_facets = ", if (length(cfg$dummy_facets %||% character(0)) > 0) render_r_object_literal(as.character(cfg$dummy_facets)) else "NULL", ","),
      paste0("  positive_facets = ", if (length(cfg$positive_facets %||% character(0)) > 0) render_r_object_literal(as.character(cfg$positive_facets)) else "NULL", ","),
      paste0("  quad_points = ", render_r_object_literal(as.integer(est_ctl$quad_points %||% 15L)), ","),
      paste0("  maxit = ", render_r_object_literal(as.integer(est_ctl$maxit %||% 400L)), ","),
      paste0("  reltol = ", render_r_object_literal(as.numeric(est_ctl$reltol %||% 1e-6))),
      ")"
    )

    if (include_diagnostics) {
      lines <- c(
        lines,
        "",
        "# Diagnostics",
        paste0("diagnostics <- diagnose_mfrm(fit, residual_pca = ", render_r_object_literal(residual_pca_mode), ")")
      )
    }
  }

  if (length(bias_pairs) > 0) {
    bias_lines <- vapply(seq_along(bias_pairs), function(i) {
      pair <- as.character(bias_pairs[[i]])
      pair <- pair[!is.na(pair) & nzchar(pair)][seq_len(min(2L, length(pair)))]
      paste0(
        "  bias_", i, " = estimate_bias(",
        "fit, diagnostics = diagnostics, facet_a = ",
        render_r_object_literal(pair[1]), ", facet_b = ",
        render_r_object_literal(pair[2]), ")",
        if (i < length(bias_pairs)) "," else ""
      )
    }, character(1))
    lines <- c(
      lines,
      "",
      "# Bias / interaction analysis",
      "bias_results <- list(",
      bias_lines,
      ")"
    )
  } else {
    lines <- c(lines, "", "bias_results <- NULL")
  }

  if (!is.null(population_prediction)) {
    pop_settings <- population_prediction$settings %||% list()
    lines <- c(
      lines,
      "",
      "# Scenario-level population forecast",
      paste0(
        "population_prediction_sim_spec <- ",
        render_classed_literal(population_prediction$sim_spec, "mfrm_sim_spec")
      ),
      "population_prediction <- predict_mfrm_population(",
      "  sim_spec = population_prediction_sim_spec,",
      paste0("  reps = ", render_r_object_literal(as.integer(pop_settings$reps %||% 50L)), ","),
      paste0("  fit_method = ", render_r_object_literal(as.character(pop_settings$fit_method %||% "MML")), ","),
      paste0("  model = ", render_r_object_literal(as.character(pop_settings$model %||% cfg$model %||% "RSM")), ","),
      paste0("  maxit = ", render_r_object_literal(as.integer(pop_settings$maxit %||% est_ctl$maxit %||% 25L)), ","),
      paste0("  quad_points = ", render_r_object_literal(as.integer(pop_settings$quad_points %||% est_ctl$quad_points %||% 7L)), ","),
      paste0("  residual_pca = ", render_r_object_literal(as.character(pop_settings$residual_pca %||% "none")), ","),
      paste0("  seed = ", if (!is.null(pop_settings$seed)) render_r_object_literal(pop_settings$seed) else "NULL"),
      ")"
    )
  } else {
    lines <- c(lines, "", "population_prediction <- NULL")
  }

  if (!is.null(unit_prediction)) {
    unit_settings <- unit_prediction$settings %||% list()
    unit_cols <- unit_settings$source_columns %||% list(
      person = "Person",
      facets = character(0),
      score = "Score",
      weight = NULL
    )
    lines <- c(
      lines,
      "",
      "# Fixed-calibration scoring for future or partially observed units",
      paste0(
        "unit_prediction_input <- ",
        render_r_object_literal(as.data.frame(unit_prediction$input_data %||% data.frame(), stringsAsFactors = FALSE))
      ),
      "unit_prediction <- predict_mfrm_units(",
      "  fit = fit,",
      "  new_data = unit_prediction_input,",
      paste0("  person = ", render_r_object_literal(as.character(unit_cols$person %||% "Person")), ","),
      paste0("  facets = ", render_r_object_literal(as.character(unit_cols$facets %||% character(0))), ","),
      paste0("  score = ", render_r_object_literal(as.character(unit_cols$score %||% "Score")), ","),
      paste0("  weight = ", if (!is.null(unit_cols$weight) && nzchar(as.character(unit_cols$weight))) render_r_object_literal(as.character(unit_cols$weight)) else "NULL", ","),
      paste0("  interval_level = ", render_r_object_literal(as.numeric(unit_settings$interval_level %||% 0.95)), ","),
      paste0("  n_draws = ", render_r_object_literal(as.integer(unit_settings$n_draws %||% 0L)), ","),
      paste0("  seed = ", if (!is.null(unit_settings$seed)) render_r_object_literal(unit_settings$seed) else "NULL"),
      ")"
    )
  } else {
    lines <- c(lines, "", "unit_prediction <- NULL")
  }

  if (!is.null(plausible_values)) {
    pv_settings <- plausible_values$settings %||% list()
    pv_cols <- pv_settings$source_columns %||% list(
      person = "Person",
      facets = character(0),
      score = "Score",
      weight = NULL
    )
    lines <- c(
      lines,
      "",
      "# Approximate plausible-value summaries under the fixed calibration",
      paste0(
        "plausible_value_input <- ",
        render_r_object_literal(as.data.frame(plausible_values$input_data %||% data.frame(), stringsAsFactors = FALSE))
      ),
      "plausible_values <- sample_mfrm_plausible_values(",
      "  fit = fit,",
      "  new_data = plausible_value_input,",
      paste0("  person = ", render_r_object_literal(as.character(pv_cols$person %||% "Person")), ","),
      paste0("  facets = ", render_r_object_literal(as.character(pv_cols$facets %||% character(0))), ","),
      paste0("  score = ", render_r_object_literal(as.character(pv_cols$score %||% "Score")), ","),
      paste0("  weight = ", if (!is.null(pv_cols$weight) && nzchar(as.character(pv_cols$weight))) render_r_object_literal(as.character(pv_cols$weight)) else "NULL", ","),
      paste0("  n_draws = ", render_r_object_literal(as.integer(pv_settings$n_draws %||% 5L)), ","),
      paste0("  interval_level = ", render_r_object_literal(as.numeric(pv_settings$interval_level %||% 0.95)), ","),
      paste0("  seed = ", if (!is.null(pv_settings$seed)) render_r_object_literal(pv_settings$seed) else "NULL"),
      ")"
    )
  } else {
    lines <- c(lines, "", "plausible_values <- NULL")
  }

  if (isTRUE(include_bundle)) {
    bundle_include <- c("core_tables", "checklist", "dashboard", "manifest", "html")
    if (!is.null(population_prediction) || !is.null(unit_prediction) || !is.null(plausible_values)) {
      bundle_include <- c(bundle_include, "predictions")
    }
    lines <- c(
      lines,
      "",
      "# Export a package-native bundle",
      "bundle <- export_mfrm_bundle(",
      "  fit = fit,",
      "  diagnostics = diagnostics,",
      "  bias_results = bias_results,",
      "  population_prediction = population_prediction,",
      "  unit_prediction = unit_prediction,",
      "  plausible_values = plausible_values,",
      paste0("  output_dir = ", render_r_object_literal(as.character(bundle_dir[1])) , ","),
      paste0("  prefix = ", render_r_object_literal(as.character(bundle_prefix[1])) , ","),
      paste0("  include = ", render_r_object_literal(bundle_include), ","),
      "  overwrite = TRUE",
      ")"
    )
  }

  script_text <- paste(lines, collapse = "\n")
  summary_tbl <- data.frame(
    ScriptMode = resolved_mode,
    ResidualPCA = residual_pca_mode,
    BiasPairs = length(bias_pairs),
    PopulationPrediction = !is.null(population_prediction),
    UnitPrediction = !is.null(unit_prediction),
    PlausibleValues = !is.null(plausible_values),
    Anchors = nrow(anchor_df),
    GroupAnchors = nrow(group_anchor_df),
    IncludeBundle = isTRUE(include_bundle),
    stringsAsFactors = FALSE
  )
  settings <- dashboard_settings_table(list(
    data_file = as.character(data_file[1]),
    script_mode = resolved_mode,
    input_mode = ctx$input_mode,
    include_bundle = isTRUE(include_bundle),
    bundle_dir = as.character(bundle_dir[1]),
    bundle_prefix = as.character(bundle_prefix[1])
  ))

  out <- list(
    summary = summary_tbl,
    script = script_text,
    settings = settings,
    anchors = anchor_df,
    group_anchors = group_anchor_df
  )
  as_mfrm_bundle(out, "mfrm_replay_script")
}

#' Export an analysis bundle for sharing or archiving
#'
#' @param fit Output from [fit_mfrm()] or [run_mfrm_facets()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. When `NULL`,
#'   diagnostics are reused from `run_mfrm_facets()` when available, otherwise
#'   computed with `residual_pca = "none"` (or `"both"` when visual summaries
#'   are requested).
#' @param bias_results Optional output from [estimate_bias()] or a named list of
#'   bias bundles.
#' @param population_prediction Optional output from
#'   [predict_mfrm_population()].
#' @param unit_prediction Optional output from [predict_mfrm_units()].
#' @param plausible_values Optional output from [sample_mfrm_plausible_values()].
#' @param output_dir Directory where files will be written.
#' @param prefix File-name prefix.
#' @param include Components to export. Supported values are
#'   `"core_tables"`, `"checklist"`, `"dashboard"`, `"apa"`, `"anchors"`,
#'   `"manifest"`, `"visual_summaries"`, `"predictions"`, `"script"`, and
#'   `"html"`.
#' @param facet Optional facet for [facet_quality_dashboard()].
#' @param include_person_anchors If `TRUE`, include person measures in the
#'   exported anchor table.
#' @param overwrite If `FALSE`, refuse to overwrite existing files.
#' @param zip_bundle If `TRUE`, attempt to zip the written files into a single
#'   archive using [utils::zip()]. This is best-effort and may depend on the
#'   local R installation.
#' @param zip_name Optional zip-file name. Defaults to `"{prefix}_bundle.zip"`.
#'
#' @details
#' This function is the package-native counterpart to the app's download bundle.
#' It reuses existing `mfrmr` helpers instead of reimplementing estimation or
#' diagnostics.
#'
#' @section Choosing exports:
#' The `include` argument lets you assemble a bundle for different audiences:
#' - `"core_tables"` for analysts who mainly want CSV output.
#' - `"manifest"` for a compact analysis record.
#' - `"script"` for reproducibility and reruns.
#' - `"html"` for a light, shareable summary page.
#' - `"visual_summaries"` when you want warning maps or residual PCA summaries
#'   to travel with the bundle.
#'
#' @section Recommended presets:
#' Common starting points are:
#' - minimal tables: `include = c("core_tables", "manifest")`
#' - reporting bundle: `include = c("core_tables", "checklist", "dashboard",
#'   "html")`
#' - archival bundle: `include = c("core_tables", "manifest", "script",
#'   "visual_summaries", "html")`
#'
#' @section Written outputs:
#' Depending on `include`, the exporter can write:
#' - core CSV tables via [export_mfrm()]
#' - checklist CSVs via [reporting_checklist()]
#' - facet-dashboard CSVs via [facet_quality_dashboard()]
#' - APA text files via [build_apa_outputs()]
#' - anchor CSV via [make_anchor_table()]
#' - manifest CSV/TXT via [build_mfrm_manifest()]
#' - visual warning/summary artifacts via [build_visual_summaries()]
#' - prediction/forecast CSVs via [predict_mfrm_population()],
#'   [predict_mfrm_units()], and [sample_mfrm_plausible_values()]
#' - a package-native replay script via [build_mfrm_replay_script()]
#' - a lightweight HTML report that bundles the exported tables/text
#'
#' @section Interpreting output:
#' The returned object reports both high-level bundle status and the exact files
#' written. In practice, `bundle$summary` is the quickest sanity check, while
#' `bundle$written_files` is the file inventory to inspect or hand off to other
#' tools.
#'
#' @section Typical workflow:
#' 1. Fit a model and compute diagnostics once.
#' 2. Decide whether the audience needs tables only, or also a manifest,
#'    replay script, and HTML summary.
#' 3. Call `export_mfrm_bundle()` with a dedicated output directory.
#' 4. Inspect `bundle$written_files` or open the generated HTML file.
#'
#' @return A named list with class `mfrm_export_bundle`.
#' @seealso [build_mfrm_manifest()], [build_mfrm_replay_script()],
#'   [export_mfrm()], [reporting_checklist()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bundle <- export_mfrm_bundle(
#'   fit,
#'   diagnostics = diag,
#'   output_dir = tempdir(),
#'   prefix = "mfrmr_bundle_example",
#'   include = c("core_tables", "manifest", "script", "html"),
#'   overwrite = TRUE
#' )
#' bundle$summary[, c("FilesWritten", "HtmlWritten", "ScriptWritten")]
#' head(bundle$written_files)
#' @export
export_mfrm_bundle <- function(fit,
                               diagnostics = NULL,
                               bias_results = NULL,
                               population_prediction = NULL,
                               unit_prediction = NULL,
                               plausible_values = NULL,
                               output_dir = ".",
                               prefix = "mfrmr_bundle",
                               include = c("core_tables", "checklist", "dashboard", "apa", "anchors", "manifest", "visual_summaries", "predictions", "script", "html"),
                               facet = NULL,
                               include_person_anchors = FALSE,
                               overwrite = FALSE,
                               zip_bundle = FALSE,
                               zip_name = NULL) {
  include <- unique(tolower(as.character(include)))
  allowed <- c("core_tables", "checklist", "dashboard", "apa", "anchors", "manifest", "visual_summaries", "predictions", "script", "html")
  bad <- setdiff(include, allowed)
  if (length(bad) > 0) {
    stop("Unsupported `include` values: ", paste(bad, collapse = ", "), ". Allowed: ", paste(allowed, collapse = ", "), call. = FALSE)
  }
  if (length(include) == 0) {
    stop("`include` must contain at least one component.", call. = FALSE)
  }

  prefix <- as.character(prefix[1])
  if (!nzchar(prefix)) prefix <- "mfrmr_bundle"
  output_dir <- as.character(output_dir[1])
  overwrite <- isTRUE(overwrite)
  zip_bundle <- isTRUE(zip_bundle)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir, call. = FALSE)
  }

  ctx <- resolve_mfrm_export_context(
    x = fit,
    diagnostics = diagnostics,
    residual_pca = if ("visual_summaries" %in% include) "both" else "none"
  )
  fit <- ctx$fit
  diagnostics <- ctx$diagnostics
  diagnostics_supplied <- ctx$diagnostics_supplied

  bias_inputs <- export_normalize_bias_inputs(bias_results)
  population_prediction <- export_validate_optional_object(
    population_prediction,
    "mfrm_population_prediction",
    "population_prediction"
  )
  unit_prediction <- export_validate_optional_object(
    unit_prediction,
    "mfrm_unit_prediction",
    "unit_prediction"
  )
  plausible_values <- export_validate_optional_object(
    plausible_values,
    "mfrm_plausible_values",
    "plausible_values"
  )
  if ("predictions" %in% include &&
      all(vapply(list(population_prediction, unit_prediction, plausible_values), is.null, logical(1)))) {
    stop(
      "`include = 'predictions'` requires at least one of `population_prediction`, `unit_prediction`, or `plausible_values`.",
      call. = FALSE
    )
  }
  written_files <- data.frame(
    Component = character(0),
    Format = character(0),
    Path = character(0),
    stringsAsFactors = FALSE
  )
  html_tables <- list()
  html_text <- list()
  visual <- NULL
  replay <- NULL

  add_written <- function(component, format, path) {
    written_files <<- rbind(
      written_files,
      data.frame(
        Component = as.character(component),
        Format = as.character(format),
        Path = normalizePath(path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    invisible(NULL)
  }

  write_csv <- function(df, filename, component) {
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("File already exists: ", path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(df, file = path, row.names = FALSE, na = "")
    add_written(component, "csv", path)
    invisible(path)
  }

  write_text <- function(text, filename, component) {
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("File already exists: ", path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    writeLines(enc2utf8(as.character(text)), con = path, useBytes = TRUE)
    add_written(component, "txt", path)
    invisible(path)
  }

  write_script <- function(text, filename, component) {
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("File already exists: ", path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    writeLines(enc2utf8(as.character(text)), con = path, useBytes = TRUE)
    add_written(component, "r", path)
    invisible(path)
  }

  write_settings_table <- function(settings, filename, component) {
    tbl <- export_flatten_named_values(settings %||% list())
    names(tbl) <- c("Setting", "Value")
    write_csv(tbl, filename, component)
  }

  write_sim_spec_bundle <- function(sim_spec, prefix_base) {
    if (!inherits(sim_spec, "mfrm_sim_spec")) return(invisible(NULL))

    scalar_settings <- list(
      source = sim_spec$source %||% "unknown",
      model = sim_spec$model %||% NA_character_,
      step_facet = sim_spec$step_facet %||% NA_character_,
      assignment = sim_spec$assignment %||% NA_character_,
      latent_distribution = sim_spec$latent_distribution %||% NA_character_,
      n_person = sim_spec$n_person %||% NA_integer_,
      n_rater = sim_spec$n_rater %||% NA_integer_,
      n_criterion = sim_spec$n_criterion %||% NA_integer_,
      raters_per_person = sim_spec$raters_per_person %||% NA_integer_,
      score_levels = sim_spec$score_levels %||% NA_integer_,
      theta_sd = sim_spec$theta_sd %||% NA_real_,
      rater_sd = sim_spec$rater_sd %||% NA_real_,
      criterion_sd = sim_spec$criterion_sd %||% NA_real_,
      noise_sd = sim_spec$noise_sd %||% NA_real_,
      step_span = sim_spec$step_span %||% NA_real_,
      group_levels = paste(as.character(sim_spec$group_levels %||% character(0)), collapse = ", ")
    )
    write_settings_table(
      scalar_settings,
      paste0(prefix, "_", prefix_base, "_settings.csv"),
      paste0(prefix_base, "_settings")
    )

    threshold_tbl <- as.data.frame(sim_spec$threshold_table %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(threshold_tbl) > 0) {
      write_csv(
        threshold_tbl,
        paste0(prefix, "_", prefix_base, "_thresholds.csv"),
        paste0(prefix_base, "_thresholds")
      )
    }

    if (is.list(sim_spec$empirical_support) && length(sim_spec$empirical_support) > 0) {
      empirical_tbl <- do.call(
        rbind,
        lapply(names(sim_spec$empirical_support), function(nm) {
          vals <- suppressWarnings(as.numeric(sim_spec$empirical_support[[nm]]))
          vals <- vals[is.finite(vals)]
          if (length(vals) == 0) return(NULL)
          data.frame(Facet = nm, Value = vals, stringsAsFactors = FALSE)
        })
      )
      if (is.data.frame(empirical_tbl) && nrow(empirical_tbl) > 0) {
        write_csv(
          empirical_tbl,
          paste0(prefix, "_", prefix_base, "_empirical_support.csv"),
          paste0(prefix_base, "_empirical_support")
        )
      }
    }

    assignment_profiles_tbl <- as.data.frame(sim_spec$assignment_profiles %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(assignment_profiles_tbl) > 0) {
      write_csv(
        assignment_profiles_tbl,
        paste0(prefix, "_", prefix_base, "_assignment_profiles.csv"),
        paste0(prefix_base, "_assignment_profiles")
      )
    }

    design_skeleton_tbl <- as.data.frame(sim_spec$design_skeleton %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(design_skeleton_tbl) > 0) {
      write_csv(
        design_skeleton_tbl,
        paste0(prefix, "_", prefix_base, "_design_skeleton.csv"),
        paste0(prefix_base, "_design_skeleton")
      )
    }

    dif_tbl <- as.data.frame(sim_spec$dif_effects %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(dif_tbl) > 0) {
      write_csv(
        dif_tbl,
        paste0(prefix, "_", prefix_base, "_dif_effects.csv"),
        paste0(prefix_base, "_dif_effects")
      )
    }

    interaction_tbl <- as.data.frame(sim_spec$interaction_effects %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(interaction_tbl) > 0) {
      write_csv(
        interaction_tbl,
        paste0(prefix, "_", prefix_base, "_interaction_effects.csv"),
        paste0(prefix_base, "_interaction_effects")
      )
    }

    if (length(sim_spec$source_summary %||% list()) > 0) {
      write_csv(
        export_flatten_named_values(sim_spec$source_summary),
        paste0(prefix, "_", prefix_base, "_source_summary.csv"),
        paste0(prefix_base, "_source_summary")
      )
    }

    invisible(NULL)
  }

  if ("core_tables" %in% include) {
    export_mfrm(
      fit = fit,
      diagnostics = diagnostics,
      output_dir = output_dir,
      prefix = prefix,
      overwrite = overwrite
    )
    add_core <- list(
      person = as.data.frame(fit$facets$person, stringsAsFactors = FALSE),
      facets = as.data.frame(fit$facets$others, stringsAsFactors = FALSE),
      summary = as.data.frame(fit$summary, stringsAsFactors = FALSE),
      measures = as.data.frame(diagnostics$measures %||% data.frame(), stringsAsFactors = FALSE)
    )
    if (nrow(as.data.frame(fit$steps, stringsAsFactors = FALSE)) > 0) {
      add_core$steps <- as.data.frame(fit$steps, stringsAsFactors = FALSE)
    }
    html_tables <- utils::modifyList(html_tables, add_core)
    core_paths <- list(
      person = file.path(output_dir, paste0(prefix, "_person_estimates.csv")),
      facets = file.path(output_dir, paste0(prefix, "_facet_estimates.csv")),
      summary = file.path(output_dir, paste0(prefix, "_fit_summary.csv")),
      measures = file.path(output_dir, paste0(prefix, "_measures.csv")),
      steps = file.path(output_dir, paste0(prefix, "_step_parameters.csv"))
    )
    for (nm in names(core_paths)) {
      if (file.exists(core_paths[[nm]])) add_written(paste0("core_", nm), "csv", core_paths[[nm]])
    }
  }

  if ("checklist" %in% include) {
    checklist <- reporting_checklist(
      fit = fit,
      diagnostics = diagnostics,
      bias_results = if (inherits(bias_results, "mfrm_bias_collection")) bias_results else bias_inputs
    )
    write_csv(checklist$checklist, paste0(prefix, "_checklist.csv"), "checklist")
    write_csv(checklist$summary, paste0(prefix, "_checklist_summary.csv"), "checklist_summary")
    if (nrow(as.data.frame(checklist$references, stringsAsFactors = FALSE)) > 0) {
      write_csv(checklist$references, paste0(prefix, "_checklist_references.csv"), "checklist_references")
    }
    html_tables$checklist <- checklist$checklist
    html_tables$checklist_summary <- checklist$summary
  }

  if ("dashboard" %in% include) {
    dash <- facet_quality_dashboard(
      fit = fit,
      diagnostics = diagnostics,
      facet = facet,
      bias_results = if (inherits(bias_results, "mfrm_bias_collection")) bias_results else bias_inputs
    )
    write_csv(dash$overview, paste0(prefix, "_facet_dashboard_overview.csv"), "dashboard_overview")
    write_csv(dash$summary, paste0(prefix, "_facet_dashboard_summary.csv"), "dashboard_summary")
    write_csv(dash$detail, paste0(prefix, "_facet_dashboard_detail.csv"), "dashboard_detail")
    if (nrow(as.data.frame(dash$flagged, stringsAsFactors = FALSE)) > 0) {
      write_csv(dash$flagged, paste0(prefix, "_facet_dashboard_flagged.csv"), "dashboard_flagged")
    }
    if (nrow(as.data.frame(dash$bias_sources, stringsAsFactors = FALSE)) > 0) {
      write_csv(dash$bias_sources, paste0(prefix, "_facet_dashboard_bias_sources.csv"), "dashboard_bias_sources")
    }
    html_tables$facet_dashboard_summary <- dash$summary
    html_tables$facet_dashboard_flagged <- dash$flagged
  }

  if ("apa" %in% include) {
    apa <- build_apa_outputs(
      fit = fit,
      diagnostics = diagnostics,
      bias_results = export_primary_bias_result(bias_results, bias_inputs)
    )
    apa_summary <- summary(apa)
    note_map_tbl <- data.frame(
      Component = names(apa$contract$note_map %||% list()),
      Text = unname(vapply(apa$contract$note_map %||% list(), as.character, character(1))),
      stringsAsFactors = FALSE
    )
    caption_map_tbl <- data.frame(
      Component = names(apa$contract$caption_map %||% list()),
      Text = unname(vapply(apa$contract$caption_map %||% list(), as.character, character(1))),
      stringsAsFactors = FALSE
    )
    section_map_tbl <- as.data.frame(apa$section_map %||% data.frame(), stringsAsFactors = FALSE)
    write_text(apa$report_text, paste0(prefix, "_apa_report.txt"), "apa_report")
    write_text(apa$table_figure_notes, paste0(prefix, "_apa_notes.txt"), "apa_notes")
    write_text(apa$table_figure_captions, paste0(prefix, "_apa_captions.txt"), "apa_captions")
    if (nrow(note_map_tbl) > 0) {
      write_csv(note_map_tbl, paste0(prefix, "_apa_note_map.csv"), "apa_note_map")
    }
    if (nrow(caption_map_tbl) > 0) {
      write_csv(caption_map_tbl, paste0(prefix, "_apa_caption_map.csv"), "apa_caption_map")
    }
    if (nrow(section_map_tbl) > 0) {
      write_csv(section_map_tbl, paste0(prefix, "_apa_sections.csv"), "apa_sections")
      html_tables$apa_sections <- section_map_tbl
    }
    if (nrow(as.data.frame(apa_summary$content_checks, stringsAsFactors = FALSE)) > 0) {
      write_csv(apa_summary$content_checks, paste0(prefix, "_apa_content_checks.csv"), "apa_content_checks")
      html_tables$apa_content_checks <- apa_summary$content_checks
    }
    html_text$apa_report <- as.character(apa$report_text)
    html_text$apa_notes <- paste(as.character(apa$table_figure_notes), collapse = "\n")
    html_text$apa_captions <- paste(as.character(apa$table_figure_captions), collapse = "\n")
  }

  if ("anchors" %in% include) {
    anchor_tbl <- make_anchor_table(
      fit = fit,
      include_person = isTRUE(include_person_anchors)
    )
    write_csv(anchor_tbl, paste0(prefix, "_anchors.csv"), "anchors")
    html_tables$anchors <- anchor_tbl
  }

  if ("predictions" %in% include) {
    if (!is.null(population_prediction)) {
      pop_sum <- summary(population_prediction, digits = 6)
      write_csv(population_prediction$design, paste0(prefix, "_population_prediction_design.csv"), "population_prediction_design")
      write_csv(population_prediction$forecast, paste0(prefix, "_population_prediction_forecast.csv"), "population_prediction_forecast")
      write_csv(population_prediction$overview, paste0(prefix, "_population_prediction_overview.csv"), "population_prediction_overview")
      write_settings_table(population_prediction$settings, paste0(prefix, "_population_prediction_settings.csv"), "population_prediction_settings")
      write_sim_spec_bundle(population_prediction$sim_spec, "population_prediction_sim_spec")
      if (!is.null(population_prediction$ademp) && length(population_prediction$ademp) > 0) {
        write_csv(
          export_flatten_named_values(population_prediction$ademp),
          paste0(prefix, "_population_prediction_ademp.csv"),
          "population_prediction_ademp"
        )
      }
      if (length(population_prediction$notes %||% character(0)) > 0) {
        write_text(
          paste(population_prediction$notes, collapse = "\n"),
          paste0(prefix, "_population_prediction_notes.txt"),
          "population_prediction_notes"
        )
      }
      html_tables$population_prediction_forecast <- pop_sum$forecast
      html_tables$population_prediction_overview <- pop_sum$overview
    }

    if (!is.null(unit_prediction)) {
      unit_sum <- summary(unit_prediction, digits = 6)
      write_csv(unit_prediction$estimates, paste0(prefix, "_unit_prediction_estimates.csv"), "unit_prediction_estimates")
      write_csv(unit_prediction$audit, paste0(prefix, "_unit_prediction_audit.csv"), "unit_prediction_audit")
      write_settings_table(unit_prediction$settings, paste0(prefix, "_unit_prediction_settings.csv"), "unit_prediction_settings")
      if (nrow(as.data.frame(unit_prediction$input_data %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
        write_csv(unit_prediction$input_data, paste0(prefix, "_unit_prediction_input.csv"), "unit_prediction_input")
      }
      if (nrow(as.data.frame(unit_prediction$draws %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
        write_csv(unit_prediction$draws, paste0(prefix, "_unit_prediction_draws.csv"), "unit_prediction_draws")
      }
      if (length(unit_prediction$notes %||% character(0)) > 0) {
        write_text(
          paste(unit_prediction$notes, collapse = "\n"),
          paste0(prefix, "_unit_prediction_notes.txt"),
          "unit_prediction_notes"
        )
      }
      html_tables$unit_prediction_estimates <- unit_sum$estimates
      html_tables$unit_prediction_audit <- unit_sum$audit
    }

    if (!is.null(plausible_values)) {
      pv_sum <- summary(plausible_values, digits = 6)
      write_csv(plausible_values$values, paste0(prefix, "_plausible_values.csv"), "plausible_values")
      write_csv(plausible_values$estimates, paste0(prefix, "_plausible_value_estimates.csv"), "plausible_value_estimates")
      write_csv(plausible_values$audit, paste0(prefix, "_plausible_value_audit.csv"), "plausible_value_audit")
      write_settings_table(plausible_values$settings, paste0(prefix, "_plausible_value_settings.csv"), "plausible_value_settings")
      if (nrow(as.data.frame(plausible_values$input_data %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
        write_csv(plausible_values$input_data, paste0(prefix, "_plausible_value_input.csv"), "plausible_value_input")
      }
      if (length(plausible_values$notes %||% character(0)) > 0) {
        write_text(
          paste(plausible_values$notes, collapse = "\n"),
          paste0(prefix, "_plausible_value_notes.txt"),
          "plausible_value_notes"
        )
      }
      if (nrow(as.data.frame(pv_sum$draw_summary %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
        html_tables$plausible_value_summary <- pv_sum$draw_summary
      }
    }
  }

  manifest <- NULL
  if ("manifest" %in% include || "html" %in% include) {
    manifest <- build_mfrm_manifest(
      fit = fit,
      diagnostics = diagnostics,
      bias_results = if (inherits(bias_results, "mfrm_bias_collection")) bias_results else bias_inputs,
      population_prediction = population_prediction,
      unit_prediction = unit_prediction,
      plausible_values = plausible_values,
      include_person_anchors = include_person_anchors
    )
  }

  if ("manifest" %in% include && !is.null(manifest)) {
    write_csv(manifest$summary, paste0(prefix, "_manifest_summary.csv"), "manifest_summary")
    write_csv(manifest$environment, paste0(prefix, "_manifest_environment.csv"), "manifest_environment")
    write_csv(manifest$model_settings, paste0(prefix, "_manifest_model_settings.csv"), "manifest_model_settings")
    write_csv(manifest$source_columns, paste0(prefix, "_manifest_source_columns.csv"), "manifest_source_columns")
    write_csv(manifest$estimation_control, paste0(prefix, "_manifest_estimation_control.csv"), "manifest_estimation_control")
    if (nrow(as.data.frame(manifest$anchor_summary, stringsAsFactors = FALSE)) > 0) {
      write_csv(manifest$anchor_summary, paste0(prefix, "_manifest_anchor_summary.csv"), "manifest_anchor_summary")
    }
    if (nrow(as.data.frame(manifest$anchors, stringsAsFactors = FALSE)) > 0) {
      write_csv(manifest$anchors, paste0(prefix, "_manifest_anchors.csv"), "manifest_anchors")
    }
    write_csv(manifest$available_outputs, paste0(prefix, "_manifest_available_outputs.csv"), "manifest_available_outputs")
    write_text(render_mfrm_manifest_text(manifest), paste0(prefix, "_manifest.txt"), "manifest_text")
    html_tables$manifest_summary <- manifest$summary
    html_tables$manifest_available_outputs <- manifest$available_outputs
  }

  if ("visual_summaries" %in% include) {
    visual <- build_visual_summaries(
      fit = fit,
      diagnostics = diagnostics,
      branch = if (ctx$input_mode == "facets_run") "facets" else "original"
    )
    write_csv(visual$warning_counts, paste0(prefix, "_visual_warning_counts.csv"), "visual_warning_counts")
    write_csv(visual$summary_counts, paste0(prefix, "_visual_summary_counts.csv"), "visual_summary_counts")
    if (nrow(as.data.frame(visual$crosswalk, stringsAsFactors = FALSE)) > 0) {
      write_csv(visual$crosswalk, paste0(prefix, "_visual_crosswalk.csv"), "visual_crosswalk")
    }
    write_text(
      render_named_text_map(visual$warning_map, title = "mfrmr Visual Warning Map"),
      paste0(prefix, "_visual_warning_map.txt"),
      "visual_warning_map"
    )
    write_text(
      render_named_text_map(visual$summary_map, title = "mfrmr Visual Summary Map"),
      paste0(prefix, "_visual_summary_map.txt"),
      "visual_summary_map"
    )
    html_text$visual_warning_map <- render_named_text_map(visual$warning_map)
    html_text$visual_summary_map <- render_named_text_map(visual$summary_map)
    html_tables$visual_warning_counts <- visual$warning_counts
    html_tables$visual_summary_counts <- visual$summary_counts
  }

  if ("script" %in% include) {
    replay <- build_mfrm_replay_script(
      fit = if (ctx$input_mode == "facets_run") ctx$run else fit,
      diagnostics = diagnostics,
      bias_results = if (inherits(bias_results, "mfrm_bias_collection")) bias_results else bias_inputs,
      population_prediction = population_prediction,
      unit_prediction = unit_prediction,
      plausible_values = plausible_values,
      include_bundle = TRUE,
      bundle_dir = "replayed_bundle",
      bundle_prefix = prefix
    )
    write_script(replay$script, paste0(prefix, "_replay.R"), "replay_script")
    html_text$replay_script <- replay$script
  }

  if ("html" %in% include) {
    html_text$manifest <- if (!is.null(manifest)) render_mfrm_manifest_text(manifest) else NULL
    html_path <- file.path(output_dir, paste0(prefix, "_bundle.html"))
    if (file.exists(html_path) && !overwrite) {
      stop("File already exists: ", html_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    html_doc <- build_mfrm_bundle_html(
      title = paste0("mfrmr Analysis Bundle: ", prefix),
      tables = html_tables,
      text_sections = html_text
    )
    writeLines(enc2utf8(html_doc), con = html_path, useBytes = TRUE)
    add_written("bundle_html", "html", html_path)
  }

  zip_written <- FALSE
  zip_path <- NULL
  zip_note <- NULL
  if (isTRUE(zip_bundle) && nrow(written_files) > 0) {
    zip_file <- if (is.null(zip_name) || !nzchar(as.character(zip_name[1]))) {
      paste0(prefix, "_bundle.zip")
    } else {
      as.character(zip_name[1])
    }
    zip_path <- file.path(output_dir, zip_file)
    if (file.exists(zip_path) && !overwrite) {
      stop("File already exists: ", zip_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    zip_inputs <- unique(normalizePath(written_files$Path, winslash = "/", mustWork = TRUE))
    zip_result <- tryCatch(
      {
        utils::zip(zipfile = zip_path, files = zip_inputs, extras = "-j")
        TRUE
      },
      error = function(e) e
    )
    if (isTRUE(zip_result) && file.exists(zip_path)) {
      add_written("bundle_zip", "zip", zip_path)
      zip_written <- TRUE
    } else if (inherits(zip_result, "error")) {
      zip_note <- conditionMessage(zip_result)
    }
  }

  summary_tbl <- data.frame(
    FilesWritten = nrow(written_files),
    CsvWritten = sum(written_files$Format == "csv"),
    TextWritten = sum(written_files$Format == "txt"),
    ScriptWritten = sum(written_files$Format == "r"),
    HtmlWritten = sum(written_files$Format == "html"),
    ZipWritten = sum(written_files$Format == "zip"),
    stringsAsFactors = FALSE
  )

  settings <- dashboard_settings_table(list(
    diagnostics_supplied = diagnostics_supplied,
    include = paste(include, collapse = ", "),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    prefix = prefix,
    overwrite = overwrite,
    zip_bundle = zip_bundle,
    zip_written = zip_written
  ))

  notes <- character(0)
  if (!is.null(zip_note) && nzchar(zip_note)) {
    notes <- c(notes, paste0("ZIP bundle was not created: ", zip_note))
  }
  if (length(notes) == 0) {
    notes <- "Bundle export completed successfully."
  }

  out <- list(
    summary = summary_tbl,
    written_files = written_files,
    manifest = manifest,
    visual_summaries = visual,
    replay_script = replay,
    settings = settings,
    notes = notes
  )
  as_mfrm_bundle(out, "mfrm_export_bundle")
}

export_normalize_bias_inputs <- function(bias_results) {
  if (is.null(bias_results)) return(list())
  if (inherits(bias_results, "mfrm_bias_collection")) {
    return(bias_results$by_pair %||% list())
  }
  if (inherits(bias_results, "mfrm_bias")) {
    return(list(bias_1 = bias_results))
  }
  if (is.list(bias_results) && !is.data.frame(bias_results)) {
    return(bias_results)
  }
  list()
}

export_primary_bias_result <- function(bias_results, bias_inputs = NULL) {
  if (inherits(bias_results, "mfrm_bias")) return(bias_results)
  if (inherits(bias_results, "mfrm_bias_collection")) {
    vals <- bias_results$by_pair %||% list()
    if (length(vals) > 0) return(vals[[1]])
    return(NULL)
  }
  vals <- bias_inputs %||% export_normalize_bias_inputs(bias_results)
  if (length(vals) > 0) vals[[1]] else NULL
}

export_validate_optional_object <- function(x, class_name, arg_name) {
  if (is.null(x)) return(NULL)
  if (!inherits(x, class_name)) {
    stop("`", arg_name, "` must be output from ", class_name, " helpers.", call. = FALSE)
  }
  x
}

export_flatten_named_values <- function(x, parent = NULL) {
  key <- if (is.null(parent) || !nzchar(parent)) "value" else parent
  if (is.null(x)) {
    return(data.frame(Key = character(0), Value = character(0), stringsAsFactors = FALSE))
  }
  if (is.data.frame(x)) {
    value <- paste(utils::capture.output(print(x, row.names = FALSE)), collapse = "\n")
    return(data.frame(Key = key, Value = value, stringsAsFactors = FALSE))
  }
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- rep("", length(x))
    parts <- lapply(seq_along(x), function(i) {
      child_nm <- nms[i]
      if (!nzchar(child_nm)) child_nm <- paste0("item", i)
      child_key <- if (identical(key, "value")) child_nm else paste(key, child_nm, sep = ".")
      export_flatten_named_values(x[[i]], child_key)
    })
    if (length(parts) == 0) {
      return(data.frame(Key = character(0), Value = character(0), stringsAsFactors = FALSE))
    }
    return(do.call(rbind, parts))
  }
  value <- if (length(x) == 0) "" else paste(as.character(x), collapse = " | ")
  data.frame(Key = key, Value = value, stringsAsFactors = FALSE)
}

export_has_residual_pca <- function(diagnostics) {
  overall <- diagnostics$residual_pca_overall %||% NULL
  by_facet <- diagnostics$residual_pca_by_facet %||% NULL
  (!is.null(overall) && length(overall) > 0) || (!is.null(by_facet) && length(by_facet) > 0)
}

export_available_outputs_table <- function(diagnostics,
                                           bias_inputs,
                                           anchor_tbl,
                                           population_prediction = NULL,
                                           unit_prediction = NULL,
                                           plausible_values = NULL) {
  pop_forecast <- as.data.frame(population_prediction$forecast %||% data.frame(), stringsAsFactors = FALSE)
  unit_estimates <- as.data.frame(unit_prediction$estimates %||% data.frame(), stringsAsFactors = FALSE)
  pv_values <- as.data.frame(plausible_values$values %||% data.frame(), stringsAsFactors = FALSE)
  data.frame(
    Component = c(
      "observed_residuals", "measures", "reliability", "residual_pca",
      "bias_results", "anchors", "population_prediction",
      "unit_prediction", "plausible_values"
    ),
    Available = c(
      !is.null(diagnostics$obs) && nrow(as.data.frame(diagnostics$obs, stringsAsFactors = FALSE)) > 0,
      !is.null(diagnostics$measures) && nrow(as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)) > 0,
      !is.null(diagnostics$reliability) && nrow(as.data.frame(diagnostics$reliability, stringsAsFactors = FALSE)) > 0,
      export_has_residual_pca(diagnostics),
      length(bias_inputs) > 0,
      nrow(as.data.frame(anchor_tbl, stringsAsFactors = FALSE)) > 0,
      !is.null(population_prediction) && nrow(pop_forecast) > 0,
      !is.null(unit_prediction) && nrow(unit_estimates) > 0,
      !is.null(plausible_values) && nrow(pv_values) > 0
    ),
    Detail = c(
      "diagnostics$obs",
      "diagnostics$measures",
      "diagnostics$reliability",
      "diagnostics$residual_pca_overall / residual_pca_by_facet",
      if (length(bias_inputs) > 0) paste0(length(bias_inputs), " bundle(s)") else "none",
      paste0(nrow(as.data.frame(anchor_tbl, stringsAsFactors = FALSE)), " row(s)"),
      if (!is.null(population_prediction)) paste0(nrow(pop_forecast), " forecast row(s)") else "none",
      if (!is.null(unit_prediction)) paste0(nrow(unit_estimates), " estimate row(s)") else "none",
      if (!is.null(plausible_values)) paste0(nrow(pv_values), " draw row(s)") else "none"
    ),
    stringsAsFactors = FALSE
  )
}

render_mfrm_manifest_text <- function(manifest) {
  sections <- list(
    Summary = manifest$summary,
    Environment = manifest$environment,
    ModelSettings = manifest$model_settings,
    SourceColumns = manifest$source_columns,
    EstimationControl = manifest$estimation_control,
    AvailableOutputs = manifest$available_outputs
  )
  if (!is.null(manifest$anchor_summary) && nrow(as.data.frame(manifest$anchor_summary, stringsAsFactors = FALSE)) > 0) {
    sections$AnchorSummary <- manifest$anchor_summary
  }
  parts <- c("mfrmr Analysis Manifest")
  for (nm in names(sections)) {
    tbl <- as.data.frame(sections[[nm]], stringsAsFactors = FALSE)
    parts <- c(parts, "", paste0("[", nm, "]"), utils::capture.output(print(tbl, row.names = FALSE)))
  }
  paste(parts, collapse = "\n")
}

build_mfrm_bundle_html <- function(title, tables = list(), text_sections = list()) {
  parts <- c(
    "<!DOCTYPE html>",
    "<html><head>",
    "<meta charset='utf-8'>",
    paste0("<title>", html_escape(title), "</title>"),
    "<style>",
    "body{font-family:system-ui,sans-serif;margin:2em;color:#222}",
    "h1{border-bottom:2px solid #333}",
    "h2{margin-top:1.5em;color:#444}",
    "table{border-collapse:collapse;margin:1em 0;width:100%}",
    "th,td{border:1px solid #ccc;padding:4px 8px;text-align:left;font-size:0.85em;vertical-align:top}",
    "th{background:#f5f5f5}",
    "tr:nth-child(even){background:#fafafa}",
    "pre{background:#f7f7f7;border:1px solid #ddd;padding:1em;white-space:pre-wrap}",
    "</style></head><body>",
    paste0("<h1>", html_escape(title), "</h1>")
  )

  text_sections <- text_sections[!vapply(text_sections, is.null, logical(1))]
  for (nm in names(text_sections)) {
    txt <- paste(as.character(text_sections[[nm]]), collapse = "\n")
    if (!nzchar(trimws(txt))) next
    parts <- c(
      parts,
      paste0("<h2>", html_escape(nm), "</h2>"),
      paste0("<pre>", html_escape(txt), "</pre>")
    )
  }

  for (nm in names(tables)) {
    tbl <- as.data.frame(tables[[nm]], stringsAsFactors = FALSE)
    if (nrow(tbl) == 0 || ncol(tbl) == 0) next
    parts <- c(
      parts,
      paste0("<h2>", html_escape(nm), "</h2>"),
      dataframe_to_html_table(tbl)
    )
  }

  parts <- c(parts, "</body></html>")
  paste(parts, collapse = "\n")
}

dataframe_to_html_table <- function(df) {
  head_html <- paste0(
    "<tr>",
    paste0("<th>", html_escape(names(df)), "</th>", collapse = ""),
    "</tr>"
  )
  body_html <- if (nrow(df) == 0) {
    "<tr><td><em>No rows</em></td></tr>"
  } else {
    paste(
      apply(df, 1, function(row) {
        paste0(
          "<tr>",
          paste0("<td>", html_escape(as.character(row)), "</td>", collapse = ""),
          "</tr>"
        )
      }),
      collapse = "\n"
    )
  }
  paste0("<table><thead>", head_html, "</thead><tbody>", body_html, "</tbody></table>")
}

html_escape <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}
