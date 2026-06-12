resolve_mfrm_export_context <- function(x,
                                        diagnostics = NULL,
                                        residual_pca = c("none", "overall", "facet", "both"),
                                        gpcm_helper = "export bundle helpers") {
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
  stop_if_gpcm_out_of_scope(fit, gpcm_helper)

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

# Build the multi-line `fit_mfrm(...)` call used in `replay.R`, from
# the `replay_inputs` list that `fit_mfrm()` stores on `fit$config`.
# The helper emits every argument that materially affects the fit
# (including the ones that were silently dropped before 0.1.6:
# `missing_codes`, `mml_engine`, `slope_facet`, `anchor_policy`,
# `min_obs_per_*`, `min_common_anchors`, `facet_shrinkage`,
# `facet_prior_sd`, `shrink_person`, `attach_diagnostics`).
#
# `anchors` and `group_anchors` are passed by name through the
# script-local variables defined earlier; latent-regression args
# (`population_formula`, `person_data`, `person_id`,
# `population_policy`) are emitted only when the population model is
# active, so the call signature mirrors the user's original.
build_replay_fit_mfrm_lines <- function(replay_inputs,
                                        fit_population,
                                        fit_population_person_id,
                                        src,
                                        cfg = list()) {
  ri <- replay_inputs %||% list()
  # If the user did not pass `rating_min` / `rating_max`, fall back
  # to the resolved value the original `fit_mfrm()` recorded on the
  # config so the replay reproduces the same range rather than
  # re-inferring from data that may have changed.
  if (is.null(ri$rating_min) && !is.null(cfg$rating_min) &&
      is.finite(suppressWarnings(as.numeric(cfg$rating_min)))) {
    ri$rating_min <- as.integer(cfg$rating_min)
  }
  if (is.null(ri$rating_max) && !is.null(cfg$rating_max) &&
      is.finite(suppressWarnings(as.numeric(cfg$rating_max)))) {
    ri$rating_max <- as.integer(cfg$rating_max)
  }

  # If the fit predates 0.1.6 (no `replay_inputs`), fall back to the
  # subset that older releases recorded. This keeps cross-version
  # bundles loadable; it does not magically recover dropped args.
  if (length(ri) == 0L) {
    ri <- list(
      person = as.character(src$person %||% "Person"),
      facets = as.character(src$facets %||% character(0)),
      score = as.character(src$score %||% "Score")
    )
  }

  emit <- function(name, value, suffix = ",", literal = NULL) {
    if (!is.null(literal)) {
      paste0("  ", name, " = ", literal, suffix)
    } else {
      paste0("  ", name, " = ", render_r_object_literal(value), suffix)
    }
  }

  # Required (non-NULL) args.
  out <- list(
    "fit <- fit_mfrm(",
    "  data = data,",
    emit("person", as.character(ri$person)),
    emit("facets", as.character(ri$facets)),
    emit("score", as.character(ri$score))
  )

  # Optional, NULL-when-absent args.
  add_optional <- function(name, value, render_fn = identity) {
    if (is.null(value)) return(NULL)
    if (is.character(value) && length(value) == 1L && !nzchar(value)) return(NULL)
    emit(name, render_fn(value))
  }

  out <- c(
    out,
    add_optional("weight", ri$weight),
    add_optional("rating_min", ri$rating_min, as.integer),
    add_optional("rating_max", ri$rating_max, as.integer),
    emit("keep_original", isTRUE(ri$keep_original)),
    add_optional("missing_codes", ri$missing_codes),
    emit("model", as.character(ri$model %||% "RSM")),
    emit("method", as.character(ri$method %||% "JML")),
    add_optional("step_facet", ri$step_facet),
    add_optional("slope_facet", ri$slope_facet),
    "  anchors = anchors,",
    "  group_anchors = group_anchors,",
    emit("noncenter_facet", as.character(ri$noncenter_facet %||% "Person")),
    add_optional("dummy_facets", ri$dummy_facets),
    add_optional("positive_facets", ri$positive_facets),
    emit("anchor_policy", as.character(ri$anchor_policy %||% "warn")),
    emit("min_common_anchors", as.integer(ri$min_common_anchors %||% 5L)),
    emit("min_obs_per_element", as.numeric(ri$min_obs_per_element %||% 30)),
    emit("min_obs_per_category", as.numeric(ri$min_obs_per_category %||% 10)),
    emit("quad_points", as.integer(ri$quad_points %||% 31L)),
    emit("maxit", as.integer(ri$maxit %||% 400L)),
    emit("reltol", as.numeric(ri$reltol %||% 1e-6)),
    emit("mml_engine", as.character(ri$mml_engine %||% "direct")),
    emit("facet_shrinkage", as.character(ri$facet_shrinkage %||% "none"))
  )

  if (!is.null(ri$facet_prior_sd)) {
    out <- c(out, emit("facet_prior_sd", as.numeric(ri$facet_prior_sd)))
  }
  out <- c(out,
    emit("shrink_person", isTRUE(ri$shrink_person)),
    emit("attach_diagnostics", isTRUE(ri$attach_diagnostics))
  )

  # Latent-regression block, conditional on the original fit having
  # used a population model.
  if (isTRUE(fit_population$active)) {
    out <- c(out,
      paste0("  population_formula = ",
             render_r_object_literal(
               fit_population$formula %||% stats::as.formula(~ 1)
             ), ","),
      "  person_data = fit_person_data,",
      paste0("  person_id = ",
             render_r_object_literal(fit_population_person_id), ","),
      paste0("  population_policy = ",
             render_r_object_literal(
               as.character(fit_population$policy %||% "error")
             ))
    )
  } else {
    # Drop the trailing comma from the last emitted line to keep the
    # call syntactically valid.
    last <- length(out)
    out[[last]] <- sub(",[[:space:]]*$", "", out[[last]])
  }

  c(unlist(out, use.names = FALSE), ")")
}

compact_population_coding_variables <- function(x) {
  x <- x %||% list()
  if (!is.list(x) || length(x) == 0L) {
    return("")
  }
  nm <- names(x)
  if (is.null(nm)) {
    nm <- character(0)
  }
  nm <- nm[nzchar(nm)]
  paste(unique(nm), collapse = ", ")
}

render_sim_spec_rebuild_call <- function(sim_spec) {
  compact_df <- function(x, require_rows = FALSE) {
    if (!is.data.frame(x)) {
      return(NULL)
    }
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (require_rows && nrow(x) == 0L) {
      return(NULL)
    }
    x
  }
  compact_vec <- function(x) {
    if (is.null(x) || length(x) == 0L) {
      return(NULL)
    }
    x
  }

  empirical_support <- sim_spec$empirical_support %||% list()
  population <- sim_spec$population %||% list(active = FALSE)
  args <- list(
    n_person = as.integer(sim_spec$n_person %||% 50L),
    n_rater = as.integer(sim_spec$n_rater %||% 4L),
    n_criterion = as.integer(sim_spec$n_criterion %||% 4L),
    raters_per_person = as.integer(sim_spec$raters_per_person %||% sim_spec$n_rater %||% 4L),
    score_levels = as.integer(sim_spec$score_levels %||% 4L),
    theta_sd = as.numeric(sim_spec$theta_sd %||% 1),
    rater_sd = as.numeric(sim_spec$rater_sd %||% 0.35),
    criterion_sd = as.numeric(sim_spec$criterion_sd %||% 0.25),
    noise_sd = as.numeric(sim_spec$noise_sd %||% 0),
    step_span = as.numeric(sim_spec$step_span %||% 1.4),
    thresholds = compact_df(sim_spec$threshold_table, require_rows = TRUE),
    model = as.character(sim_spec$model %||% "RSM"),
    step_facet = as.character(sim_spec$step_facet %||% NA_character_),
    slope_facet = compact_vec(as.character(sim_spec$slope_facet %||% NULL)),
    slopes = compact_df(sim_spec$slope_table, require_rows = TRUE),
    facet_names = compact_vec(as.character(sim_spec$facet_names %||% NULL)),
    assignment = as.character(sim_spec$assignment %||% "crossed"),
    latent_distribution = as.character(sim_spec$latent_distribution %||% "normal"),
    empirical_person = compact_vec(empirical_support$person %||% NULL),
    empirical_rater = compact_vec(empirical_support$rater %||% NULL),
    empirical_criterion = compact_vec(empirical_support$criterion %||% NULL),
    assignment_profiles = compact_df(sim_spec$assignment_profiles, require_rows = TRUE),
    design_skeleton = compact_df(sim_spec$design_skeleton, require_rows = TRUE),
    group_levels = compact_vec(as.character(sim_spec$group_levels %||% NULL)),
    dif_effects = compact_df(sim_spec$dif_effects, require_rows = TRUE),
    interaction_effects = compact_df(sim_spec$interaction_effects, require_rows = TRUE),
    population_formula = if (isTRUE(population$active)) population$formula else NULL,
    population_coefficients = if (isTRUE(population$active)) population$coefficients else NULL,
    population_sigma2 = if (isTRUE(population$active)) population$sigma2 else NULL,
    population_covariates = if (isTRUE(population$active)) compact_df(population$covariate_template) else NULL
  )
  args <- args[!vapply(args, is.null, logical(1))]
  arg_lines <- vapply(seq_along(args), function(i) {
    paste0(
      "  ",
      names(args)[i],
      " = ",
      render_r_object_literal(args[[i]]),
      if (i < length(args)) "," else ""
    )
  }, character(1))
  paste(c("build_mfrm_sim_spec(", arg_lines, ")"), collapse = "\n")
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

validate_bias_results_input <- function(bias_results,
                                        helper = "helper()") {
  if (is.null(bias_results)) {
    return(NULL)
  }
  if (inherits(bias_results, "mfrm_bias") ||
      inherits(bias_results, "mfrm_bias_collection")) {
    return(bias_results)
  }
  if (is.list(bias_results) && !is.data.frame(bias_results)) {
    if (length(bias_results) == 0L) {
      return(bias_results)
    }
    valid <- vapply(
      bias_results,
      function(obj) is.null(obj) || inherits(obj, "mfrm_bias"),
      logical(1)
    )
    if (all(valid)) {
      return(bias_results)
    }
    bad_names <- names(bias_results)
    if (is.null(bad_names)) {
      bad_names <- paste0("[[", seq_along(bias_results), "]]")
    } else {
      bad_names[is.na(bad_names) | !nzchar(bad_names)] <- paste0("[[", which(is.na(bad_names) | !nzchar(bad_names)), "]]")
    }
    bad_labels <- bad_names[!valid]
    stop(
      "`bias_results` in ", helper, " must be NULL, output from estimate_bias(), an `mfrm_bias_collection`, or a list of `mfrm_bias` objects. Invalid entries: ",
      paste(bad_labels, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  stop(
    "`bias_results` in ", helper, " must be NULL, output from estimate_bias(), an `mfrm_bias_collection`, or a list of `mfrm_bias` objects.",
    call. = FALSE
  )
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
#' @param data Optional original analysis data frame. When supplied,
#'   the manifest's `input_hash` row for `data` is computed against
#'   the user's untouched input rather than the package's internal
#'   `prep$data` (which carries synthesised `Weight` / `score_k`
#'   columns) so the recorded fingerprint matches what
#'   `read.csv()` will produce in a replay session.
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
#' - `hierarchical_review`: retained traceability table for hierarchical /
#'   small-sample design flags
#' - `missing_recoding`: retained traceability table for missing-code recoding
#' - `shrinkage_review`: retained traceability table for shrinkage settings
#' - `available_outputs`: availability table for diagnostics/bias/PCA/prediction
#'   outputs
#' - `dependencies`, `input_hash`, and `session_info`: reproducibility
#'   metadata tables
#' - `settings`: manifest build settings
#'
#' @section Interpreting output:
#' The `summary` table is the direct place to confirm that you are looking at
#' the intended analysis. The `model_settings`, `source_columns`, and
#' `estimation_control` tables are designed for reproducibility records and
#' method write-up.
#' Active latent-regression fits also record their population-model provenance
#' there, including the fitted scoring basis, stored `population_formula`, and
#' person-level contract used by the fitted population model. When categorical
#' background variables are expanded through `stats::model.matrix()`,
#' `population_xlevel_variables` and `population_contrast_variables` identify
#' the variables whose fitted coding must be preserved for replay/scoring.
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
#' For bounded `GPCM` fits, the manifest is available with an explicit
#' `gpcm_boundary` table. It records supported direct diagnostics/reporting
#' surfaces while keeping full FACETS score-side contract review blocked and
#' routing design forecasting through its separate caveated capability row.
#'
#' @return A named list with class `mfrm_manifest`.
#' @seealso [export_mfrm_bundle()], [build_mfrm_replay_script()],
#'   [make_anchor_table()], [reporting_checklist()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
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
                                include_person_anchors = FALSE,
                                data = NULL) {
  ctx <- resolve_mfrm_export_context(
    x = fit,
    diagnostics = diagnostics,
    residual_pca = "none",
    gpcm_helper = "build_mfrm_manifest()"
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
  unit_settings <- unit_prediction$settings %||% list()
  pv_settings <- plausible_values$settings %||% list()
  anchor_tbl <- make_anchor_table(
    fit = fit,
    include_person = isTRUE(include_person_anchors)
  )

  cfg <- fit$config %||% list()
  prep <- fit$prep %||% list()
  est_ctl <- cfg$estimation_control %||% list()
  fit_population <- fit$population %||% list(
    active = isTRUE(cfg$population_active),
    posterior_basis = cfg$posterior_basis %||% "legacy_mml",
    formula = cfg$population_formula %||% NULL,
    policy = cfg$population_policy %||% NULL,
    person_id = cfg$source_columns$person %||% "Person",
    design_columns = NULL,
    person_table = NULL,
    person_table_replay = NULL,
    person_table_replay_scope = NULL
  )

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
    FitPopulationActive = isTRUE(fit_population$active),
    FitPosteriorBasis = as.character(fit_population$posterior_basis %||% "legacy_mml"),
    Converged = isTRUE(fit$summary$Converged %||% FALSE),
    stringsAsFactors = FALSE
  )

  # RNG seed snapshot. .Random.seed may not exist yet (e.g. if no RNG
  # has been touched); `set.seed(NULL)` to prime it would change state,
  # so we just record the absence instead.
  rng_state <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    seed_vec <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    paste(utils::head(seed_vec, 8L), collapse = ",")
  } else {
    NA_character_
  }
  rng_kind <- paste(RNGkind(), collapse = ", ")

  environment_tbl <- data.frame(
    Package = "mfrmr",
    PackageVersion = as.character(utils::packageVersion("mfrmr")),
    RVersion = as.character(getRversion()),
    Platform = R.version$platform,
    Timestamp = format(as.POSIXct(Sys.time()), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    RNGKind = rng_kind,
    RNGSeedDigest = rng_state,
    Locale = Sys.getlocale("LC_CTYPE"),
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
    n_categories = as.character(cfg$n_cat %||% NA_character_),
    population_active = isTRUE(fit_population$active),
    posterior_basis = as.character(fit_population$posterior_basis %||% "legacy_mml"),
    population_formula = if (!is.null(fit_population$formula)) {
      paste(deparse(fit_population$formula), collapse = " ")
    } else {
      NULL
    },
    population_person_id = as.character(fit_population$person_id %||% NULL),
    population_policy = as.character(fit_population$policy %||% NULL),
    population_design_columns = paste(as.character(fit_population$design_columns %||% character(0)), collapse = ", "),
    population_xlevel_variables = compact_population_coding_variables(fit_population$xlevels),
    population_contrast_variables = compact_population_coding_variables(fit_population$contrasts),
    population_person_rows = if (is.data.frame(fit_population$person_table)) {
      as.integer(nrow(fit_population$person_table))
    } else {
      NA_integer_
    },
    population_person_rows_replay = if (is.data.frame(fit_population$person_table_replay)) {
      as.integer(nrow(fit_population$person_table_replay))
    } else {
      NA_integer_
    },
    population_person_replay_scope = as.character(fit_population$person_table_replay_scope %||% NULL),
    population_omitted_persons = length(as.character(fit_population$omitted_persons %||% character(0))),
    population_response_rows_omitted = as.integer(fit_population$response_rows_omitted %||% 0L)
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
    fit_population_active = isTRUE(fit_population$active),
    fit_posterior_basis = as.character(fit_population$posterior_basis %||% "legacy_mml"),
    fit_population_formula = if (!is.null(fit_population$formula)) {
      paste(deparse(fit_population$formula), collapse = " ")
    } else {
      NULL
    },
    fit_population_person_id = as.character(fit_population$person_id %||% NULL),
    fit_population_policy = as.character(fit_population$policy %||% NULL),
    fit_population_xlevel_variables = compact_population_coding_variables(fit_population$xlevels),
    fit_population_contrast_variables = compact_population_coding_variables(fit_population$contrasts),
    fit_population_person_rows_replay = if (is.data.frame(fit_population$person_table_replay)) {
      as.integer(nrow(fit_population$person_table_replay))
    } else {
      NA_integer_
    },
    fit_population_person_replay_scope = as.character(fit_population$person_table_replay_scope %||% NULL),
    fit_population_omitted_persons = length(as.character(fit_population$omitted_persons %||% character(0))),
    fit_population_response_rows_omitted = as.integer(fit_population$response_rows_omitted %||% 0L),
    unit_prediction_posterior_basis = as.character(unit_settings$posterior_basis %||% NULL),
    unit_prediction_person_id = as.character(unit_settings$person_id %||% NULL),
    unit_prediction_population_policy = as.character(unit_settings$population_policy %||% NULL),
    unit_prediction_population_formula = as.character(unit_settings$population_formula %||% NULL),
    plausible_value_posterior_basis = as.character(pv_settings$posterior_basis %||% NULL),
    plausible_value_person_id = as.character(pv_settings$person_id %||% NULL),
    plausible_value_population_policy = as.character(pv_settings$population_policy %||% NULL),
    plausible_value_population_formula = as.character(pv_settings$population_formula %||% NULL),
    input_mode = ctx$input_mode
  ))

  # Hierarchical structure review (0.1.6). Lightweight descriptive
  # summary so manifests document the observed nesting / sample-adequacy
  # state even when the full review table is not stored alongside the fit.
  hierarchical_review <- tryCatch({
    flag <- as.character(fit$summary$FacetSampleSizeFlag %||% NA_character_)
    min_n <- suppressWarnings(as.integer(fit$summary$FacetMinLevelN %||% NA_integer_))
    sparse_n <- suppressWarnings(as.integer(fit$summary$FacetSparseCount %||% NA_integer_))
    ext_hi <- suppressWarnings(as.integer(fit$summary$ExtremeHighN %||% NA_integer_))
    ext_lo <- suppressWarnings(as.integer(fit$summary$ExtremeLowN %||% NA_integer_))
    data.frame(
      FacetSampleSizeFlag = flag,
      FacetMinLevelN = min_n,
      FacetSparseCount = sparse_n,
      # Extreme-person counts (added in 0.1.6). These were computed
      # in fit$summary but not surfaced in any bundle; manifest consumers
      # need them to review JML extreme-score behaviour.
      ExtremeHighN = ext_hi,
      ExtremeLowN = ext_lo,
      Note = paste0(
        "FacetSampleSizeFlag aggregates non-Person facets only (see ",
        "?fit_mfrm 'Fixed effects assumption'). For Person-level and ",
        "per-level detail, run `facet_small_sample_review(fit)` and ",
        "`analyze_hierarchical_structure(data, facets)`."
      ),
      stringsAsFactors = FALSE
    )
  }, error = function(e) data.frame(
    FacetSampleSizeFlag = NA_character_,
    FacetMinLevelN = NA_integer_,
    FacetSparseCount = NA_integer_,
    ExtremeHighN = NA_integer_,
    ExtremeLowN = NA_integer_,
    Note = "Hierarchical review summary unavailable.",
    stringsAsFactors = FALSE
  ))

  # Surface any missing-code recoding that `recode_missing_codes()` or
  # `fit_mfrm(..., missing_codes = ...)` performed before the fit. The
  # integrated path stores the review in `prep$missing_recoding`; the
  # older helper-only path leaves it on an attr() of `prep$data`.
  missing_recoding_review <- tryCatch({
    review <- prep$missing_recoding
    if (is.null(review)) {
      review <- attr(prep$data, "mfrm_missing_recoding")
    }
    if (is.null(review) || !is.data.frame(review) || nrow(review) == 0L) {
      data.frame(
        Column = character(0),
        Replaced = integer(0),
        stringsAsFactors = FALSE
      )
    } else {
      review
    }
  }, error = function(e) data.frame(
    Column = character(0),
    Replaced = integer(0),
    stringsAsFactors = FALSE
  ))

  # Reproducibility additions (0.1.6 polish). These tables are always
  # populated (never NULL) so downstream writers can rely on the column
  # contract, but individual fields may be NA when the underlying
  # helper is unavailable (e.g. `digest` in Suggests but not installed).
  dependencies_tbl <- build_mfrm_dependency_table(
    c("dplyr", "tidyr", "tibble", "purrr", "stringr", "psych",
      "lifecycle", "rlang", "cpp11",
      "igraph", "lme4", "digest", "kableExtra", "flextable")
  )

  # Hash the user's *original* data when supplied, so the recorded
  # fingerprint matches what `read.csv()` will produce in the replay
  # session. `prep$data` carries synthesized `Weight` / `score_k`
  # columns that the user could not reconstruct from their CSV.
  hash_data <- if (is.data.frame(data) && nrow(data) > 0L) {
    data
  } else {
    prep$data %||% NULL
  }
  input_hash_tbl <- build_mfrm_input_hash_table(
    data = hash_data,
    anchors = anchor_tbl,
    group_anchors = fit$config$group_anchors %||% NULL,
    score_map = prep$score_map %||% NULL
  )

  session_info_tbl <- build_mfrm_session_info_table()

  # Empirical-Bayes shrinkage trail (0.1.6). Records the mode and per-facet
  # tau^2 / mean shrinkage values so reproducibility bundles mirror
  # what `build_apa_outputs()` puts in the Methods narrative. When no
  # shrinkage is applied, return a single-row sentinel so downstream
  # writers always see a schema with Mode + NA values rather than a
  # zero-row frame that confuses column-count checks.
  empty_shrinkage_review <- function(mode_val) {
    data.frame(
      Mode = mode_val,
      Facet = NA_character_,
      Tau2 = NA_real_,
      MeanShrinkage = NA_real_,
      EffectiveDF = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  shrinkage_review <- tryCatch({
    mode <- as.character(fit$config$facet_shrinkage %||% "none")
    report <- fit$shrinkage_report
    if (identical(mode, "none") || is.null(report) || !is.data.frame(report) ||
        nrow(report) == 0L) {
      empty_shrinkage_review(mode)
    } else {
      data.frame(
        Mode = rep(mode, nrow(report)),
        Facet = as.character(report$Facet),
        Tau2 = suppressWarnings(as.numeric(report$Tau2)),
        MeanShrinkage = suppressWarnings(as.numeric(report$MeanShrinkage)),
        EffectiveDF = suppressWarnings(as.numeric(report$EffectiveDF)),
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) empty_shrinkage_review(NA_character_))

  out <- list(
    summary = summary_tbl,
    environment = environment_tbl,
    model_settings = model_settings,
    source_columns = source_columns,
    estimation_control = estimation_control,
    anchor_summary = anchor_summary,
    anchors = as.data.frame(anchor_tbl, stringsAsFactors = FALSE),
    hierarchical_review = hierarchical_review,
    missing_recoding = missing_recoding_review,
    shrinkage_review = shrinkage_review,
    dependencies = dependencies_tbl,
    input_hash = input_hash_tbl,
    session_info = session_info_tbl,
    available_outputs = available_outputs,
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "build_mfrm_manifest()",
      extra_areas = c(
        "Score-side scorefile export under bounded GPCM",
        "FACETS output-contract score-side review",
        "Design planning and forecasting"
      )
    ),
    settings = settings
  )
  as_mfrm_bundle(out, "mfrm_manifest")
}

# ---- internal helpers for manifest reproducibility (0.1.6 polish) --------

#' @keywords internal
#' @noRd
build_mfrm_dependency_table <- function(pkgs) {
  pkgs <- as.character(pkgs %||% character(0))
  if (length(pkgs) == 0L) {
    return(data.frame(Package = character(0),
                      Version = character(0),
                      Role = character(0),
                      stringsAsFactors = FALSE))
  }
  desc <- utils::packageDescription("mfrmr")
  imports <- trimws(strsplit(as.character(desc$Imports %||% ""), ",")[[1]])
  suggests <- trimws(strsplit(as.character(desc$Suggests %||% ""), ",")[[1]])
  # Strip version tags like "pkg (>= 1.0)" to bare names.
  imports <- sub("\\s*\\(.*\\)$", "", imports)
  suggests <- sub("\\s*\\(.*\\)$", "", suggests)
  imports <- imports[nzchar(imports)]
  suggests <- suggests[nzchar(suggests)]

  rows <- lapply(pkgs, function(p) {
    role <- if (p %in% imports) "imports"
            else if (p %in% suggests) "suggests"
            else "other"
    ver <- tryCatch(as.character(utils::packageVersion(p)),
                    error = function(e) NA_character_)
    data.frame(Package = p,
               Version = ver,
               Role = role,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' @keywords internal
#' @noRd
build_mfrm_input_hash_table <- function(data = NULL, anchors = NULL,
                                        group_anchors = NULL,
                                        score_map = NULL) {
  hash_one <- function(x) {
    if (is.null(x)) return(NA_character_)
    if (requireNamespace("digest", quietly = TRUE)) {
      return(digest::digest(x, algo = "sha256"))
    }
    # Fallback: serialize + tools::md5sum via a temp file so we still
    # record a stable fingerprint even when digest is absent.
    tf <- tempfile()
    on.exit(unlink(tf), add = TRUE)
    saveRDS(x, tf)
    unname(tools::md5sum(tf))
  }
  n_rows <- function(x) if (is.null(x)) NA_integer_ else {
    if (is.data.frame(x)) nrow(x) else length(x)
  }

  data.frame(
    Object = c("data", "anchors", "group_anchors", "score_map"),
    Hash = c(hash_one(data),
             hash_one(anchors),
             hash_one(group_anchors),
             hash_one(score_map)),
    NRows = c(n_rows(data),
              n_rows(anchors),
              n_rows(group_anchors),
              n_rows(score_map)),
    Algorithm = c(rep(
      if (requireNamespace("digest", quietly = TRUE)) "sha256"
      else "md5-of-rds",
      4L)),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
build_mfrm_session_info_table <- function() {
  info <- utils::sessionInfo()
  attached <- if (!is.null(info$otherPkgs)) {
    vapply(info$otherPkgs, function(p) as.character(p$Version %||% NA_character_),
           character(1))
  } else character(0)
  loaded <- if (!is.null(info$loadedOnly)) {
    vapply(info$loadedOnly, function(p) as.character(p$Version %||% NA_character_),
           character(1))
  } else character(0)

  pkg_rows <- function(tbl, scope) {
    if (length(tbl) == 0L) return(NULL)
    data.frame(
      Scope = scope,
      Package = names(tbl),
      Version = as.character(unname(tbl)),
      stringsAsFactors = FALSE
    )
  }
  pkg_df <- rbind(pkg_rows(attached, "attached"),
                  pkg_rows(loaded, "loaded"))

  meta_df <- data.frame(
    Scope = c("platform", "running", "base"),
    Package = c("platform", "os", "base_R"),
    Version = c(as.character(info$platform %||% NA_character_),
                as.character(info$running %||% NA_character_),
                paste0(info$R.version$major, ".", info$R.version$minor)),
    stringsAsFactors = FALSE
  )
  res <- rbind(meta_df, pkg_df %||%
                 data.frame(Scope = character(0),
                            Package = character(0),
                            Version = character(0),
                            stringsAsFactors = FALSE))
  rownames(res) <- NULL
  res
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
#' @param fit_person_data_file Optional CSV filename to read for the fit-level
#'   latent-regression replay person table. When `NULL`, the replay script
#'   embeds that table inline. [export_mfrm_bundle()] uses this to keep replay
#'   scripts portable while avoiding large inline literals.
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
#' When the supplied fit uses the latent-regression `MML` branch, the generated
#' fit-mode script also carries the stored replay-ready person table together
#' with the corresponding `population_formula` / `person_id` /
#' `population_policy` arguments needed to recreate the population model.
#' By default that replay-ready table is embedded inline; when
#' `fit_person_data_file` is supplied, the generated script reads it from that
#' sidecar CSV relative to the replay script location.
#'
#' For bounded `GPCM`, replay scripts are available with an explicit
#' `gpcm_boundary` table. The generated script records `step_facet` and
#' `slope_facet` settings, but full FACETS score-side contract review remains
#' outside this replay contract. Role-based design forecasting is available
#' through [predict_mfrm_population()] as a separate caveated helper route.
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
#'                 method = "JML", maxit = 30)
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
                                     fit_person_data_file = NULL,
                                     script_mode = c("auto", "fit", "facets"),
                                     include_bundle = FALSE,
                                     bundle_dir = "analysis_bundle",
                                     bundle_prefix = "mfrmr_replay") {
  script_mode <- match.arg(tolower(as.character(script_mode[1] %||% "auto")),
                           c("auto", "fit", "facets"))
  ctx <- resolve_mfrm_export_context(
    x = fit,
    diagnostics = diagnostics,
    residual_pca = "none",
    gpcm_helper = "build_mfrm_replay_script()"
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
  unit_settings <- unit_prediction$settings %||% list()
  pv_settings <- plausible_values$settings %||% list()
  fit_population <- fit$population %||% list(
    active = isTRUE(cfg$population_active),
    posterior_basis = cfg$posterior_basis %||% "legacy_mml",
    formula = cfg$population_formula %||% NULL,
    policy = cfg$population_policy %||% NULL,
    person_id = cfg$source_columns$person %||% "Person",
    person_table = NULL,
    person_table_replay = NULL,
    person_table_replay_scope = NULL
  )
  fit_population_person_id <- as.character(fit_population$person_id %||% src$person %||% "Person")
  fit_population_person_table <- fit_population$person_table
  fit_population_person_table_replay <- fit_population$person_table_replay %||% fit_population_person_table
  fit_person_data_file <- as.character(fit_person_data_file[1] %||% "")
  if (!nzchar(fit_person_data_file)) fit_person_data_file <- NULL
  if (isTRUE(fit_population$active)) {
    if (!is.data.frame(fit_population_person_table_replay)) {
      stop(
        "Replay generation for latent-regression fits requires the stored replay-ready person-level background-data table.",
        call. = FALSE
      )
    }
    if (!(fit_population_person_id %in% names(fit_population_person_table_replay))) {
      stop(
        "Replay generation for latent-regression fits requires the stored replay-ready person table to retain the fitted `person_id` column.",
        call. = FALSE
      )
    }
    replay_person_ids <- unique(as.character(fit_population_person_table_replay[[fit_population_person_id]]))
    required_person_ids <- unique(c(
      as.character(fit_population$included_persons %||% character(0)),
      as.character(fit_population$omitted_persons %||% character(0))
    ))
    missing_replay_ids <- setdiff(required_person_ids, replay_person_ids)
    if (length(missing_replay_ids) > 0L ||
        (!is.data.frame(fit_population$person_table_replay) &&
          (length(as.character(fit_population$omitted_persons %||% character(0))) > 0L ||
            isTRUE(as.integer(fit_population$response_rows_omitted %||% 0L) > 0L)))) {
      stop(
        "Replay generation cannot reproduce this latent-regression fit because the stored replay-ready person table no longer covers every observed person in `data`. ",
        "This commonly occurs when the fit used `population_policy = 'omit'` before replay-ready background-data provenance was stored.",
        call. = FALSE
      )
    }
  }

  recorded_pkg_version <- as.character(
    cfg$replay_inputs$package_version %||% utils::packageVersion("mfrmr")
  )
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
    paste0("# Recorded mfrmr version: ", recorded_pkg_version),
    "",
    "library(mfrmr)",
    "",
    "# Version-mismatch warning. The replay script was generated against",
    "# the recorded mfrmr version above. If the installed version differs",
    "# the fit will still run, but estimates may diverge.",
    paste0(
      "if (utils::packageVersion(\"mfrmr\") != \"", recorded_pkg_version,
      "\") warning(\"replay.R was generated under mfrmr ",
      recorded_pkg_version,
      "; you are running \", as.character(utils::packageVersion(\"mfrmr\")), \". \", ",
      "\"Estimates may differ.\")"
    ),
    "",
    paste0("data <- utils::read.csv(", render_r_object_literal(as.character(data_file[1])), ", stringsAsFactors = FALSE)")
  )

  lines <- c(
    lines,
    "",
    "# Stored constraints from the fitted analysis",
    paste0("anchors <- ", if (nrow(anchor_df) > 0) render_r_object_literal(anchor_df) else "NULL"),
    paste0("group_anchors <- ", if (nrow(group_anchor_df) > 0) render_r_object_literal(group_anchor_df) else "NULL"),
    "",
    "# Population-model basis recorded in the fitted analysis",
    paste0("# population_active = ", ifelse(isTRUE(fit_population$active), "TRUE", "FALSE")),
    paste0("# posterior_basis = ", as.character(fit_population$posterior_basis %||% "legacy_mml"))
  )
  if (isTRUE(fit_population$active)) {
    lines <- c(
      lines,
      paste0("# population_formula = ", as.character(fit_population$formula %||% NA_character_)),
      paste0("# population_person_id = ", fit_population_person_id),
      paste0("# population_policy = ", as.character(fit_population$policy %||% NA_character_)),
      paste0("# population_person_replay_scope = ", as.character(fit_population$person_table_replay_scope %||% NA_character_))
    )
  }

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
      if (isTRUE(fit_population$active) &&
          !is.null(fit_population_person_table_replay) &&
          is.data.frame(fit_population_person_table_replay)) {
        if (!is.null(fit_person_data_file)) {
          lines <- c(
            lines,
            "",
            "# Person-level background data used by the fitted latent-regression branch",
            "replay_script_args <- commandArgs(trailingOnly = FALSE)",
            "replay_script_file <- replay_script_args[grepl('^--file=', replay_script_args)]",
            "replay_script_path <- if (length(replay_script_file) > 0) {",
            "  sub('^--file=', '', replay_script_file[1])",
            "} else {",
            "  tryCatch(sys.frames()[[1]]$ofile, error = function(e) '')",
            "}",
            "replay_script_dir <- if (is.character(replay_script_path) && nzchar(replay_script_path)) {",
            "  dirname(normalizePath(replay_script_path, winslash = '/', mustWork = FALSE))",
            "} else {",
            "  getwd()",
            "}",
            paste0(
              "fit_person_data <- utils::read.csv(file.path(replay_script_dir, ",
              render_r_object_literal(basename(fit_person_data_file)),
              "), stringsAsFactors = FALSE)"
            )
          )
        } else {
          lines <- c(
            lines,
            "",
            "# Person-level background data used by the fitted latent-regression branch",
            paste0(
              "fit_person_data <- ",
              render_r_object_literal(as.data.frame(fit_population_person_table_replay, stringsAsFactors = FALSE))
            )
          )
        }
      }
      # Build the `fit_mfrm()` call from the captured `replay_inputs`
      # so every argument that affected the original fit appears in the
      # script, not only the most common subset. `replay_inputs` is
      # written by `fit_mfrm()` itself (see api-estimation.R) and is
      # the single source of truth for replay.
      # `cfg$rating_min` / `cfg$rating_max` carry the resolved range
      # used by the original fit (whether the user supplied them or
      # `prepare_mfrm_data()` inferred them). Pass them as fallback
      # so the replayed fit pins the same range and does not depend
      # on re-inference from a possibly-changed input file.
      fit_call_lines <- build_replay_fit_mfrm_lines(
        replay_inputs = cfg$replay_inputs,
        fit_population = fit_population,
        fit_population_person_id = fit_population_person_id,
        src = src,
        cfg = cfg
      )
      lines <- c(
        lines,
        "",
        "# Fit the model",
        fit_call_lines
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
        render_sim_spec_rebuild_call(population_prediction$sim_spec)
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
    unit_person_id <- unit_settings$person_id %||% NULL
    unit_population_policy <- unit_settings$population_policy %||% NULL
    lines <- c(
      lines,
      "",
      "# Scoring for future or partially observed units",
      paste0(
        "unit_prediction_input <- ",
        render_r_object_literal(as.data.frame(unit_prediction$input_data %||% data.frame(), stringsAsFactors = FALSE))
      )
    )
    if (!is.null(unit_prediction$person_data)) {
      lines <- c(
        lines,
        paste0(
          "unit_prediction_person_data <- ",
          render_r_object_literal(as.data.frame(unit_prediction$person_data, stringsAsFactors = FALSE))
        )
      )
    }
    lines <- c(
      lines,
      "unit_prediction <- predict_mfrm_units(",
      "  fit = fit,",
      "  new_data = unit_prediction_input,",
      paste0("  person = ", render_r_object_literal(as.character(unit_cols$person %||% "Person")), ","),
      paste0("  facets = ", render_r_object_literal(as.character(unit_cols$facets %||% character(0))), ","),
      paste0("  score = ", render_r_object_literal(as.character(unit_cols$score %||% "Score")), ","),
      paste0("  weight = ", if (!is.null(unit_cols$weight) && nzchar(as.character(unit_cols$weight))) render_r_object_literal(as.character(unit_cols$weight)) else "NULL", ","),
      paste0("  person_data = ", if (!is.null(unit_prediction$person_data)) "unit_prediction_person_data" else "NULL", ","),
      paste0("  person_id = ", if (!is.null(unit_person_id) && nzchar(as.character(unit_person_id))) render_r_object_literal(as.character(unit_person_id)) else "NULL", ","),
      paste0("  population_policy = ", if (!is.null(unit_population_policy) && nzchar(as.character(unit_population_policy))) render_r_object_literal(as.character(unit_population_policy)) else "NULL", ","),
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
    pv_person_id <- pv_settings$person_id %||% NULL
    pv_population_policy <- pv_settings$population_policy %||% NULL
    lines <- c(
      lines,
      "",
      "# Approximate plausible-value summaries under posterior scoring",
      paste0(
        "plausible_value_input <- ",
        render_r_object_literal(as.data.frame(plausible_values$input_data %||% data.frame(), stringsAsFactors = FALSE))
      )
    )
    if (!is.null(plausible_values$person_data)) {
      lines <- c(
        lines,
        paste0(
          "plausible_value_person_data <- ",
          render_r_object_literal(as.data.frame(plausible_values$person_data, stringsAsFactors = FALSE))
        )
      )
    }
    lines <- c(
      lines,
      "plausible_values <- sample_mfrm_plausible_values(",
      "  fit = fit,",
      "  new_data = plausible_value_input,",
      paste0("  person = ", render_r_object_literal(as.character(pv_cols$person %||% "Person")), ","),
      paste0("  facets = ", render_r_object_literal(as.character(pv_cols$facets %||% character(0))), ","),
      paste0("  score = ", render_r_object_literal(as.character(pv_cols$score %||% "Score")), ","),
      paste0("  weight = ", if (!is.null(pv_cols$weight) && nzchar(as.character(pv_cols$weight))) render_r_object_literal(as.character(pv_cols$weight)) else "NULL", ","),
      paste0("  person_data = ", if (!is.null(plausible_values$person_data)) "plausible_value_person_data" else "NULL", ","),
      paste0("  person_id = ", if (!is.null(pv_person_id) && nzchar(as.character(pv_person_id))) render_r_object_literal(as.character(pv_person_id)) else "NULL", ","),
      paste0("  population_policy = ", if (!is.null(pv_population_policy) && nzchar(as.character(pv_population_policy))) render_r_object_literal(as.character(pv_population_policy)) else "NULL", ","),
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
    FitPopulationActive = isTRUE(fit_population$active),
    FitPosteriorBasis = as.character(fit_population$posterior_basis %||% "legacy_mml"),
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
    fit_population_active = isTRUE(fit_population$active),
    fit_posterior_basis = as.character(fit_population$posterior_basis %||% "legacy_mml"),
    fit_population_formula = if (!is.null(fit_population$formula)) {
      paste(deparse(fit_population$formula), collapse = " ")
    } else {
      NULL
    },
    fit_population_person_id = fit_population_person_id,
    fit_population_policy = as.character(fit_population$policy %||% NULL),
    fit_population_person_rows_replay = if (is.data.frame(fit_population$person_table_replay)) {
      as.integer(nrow(fit_population$person_table_replay))
    } else {
      NA_integer_
    },
    fit_population_person_data_mode = if (!is.null(fit_person_data_file)) "sidecar_csv" else if (isTRUE(fit_population$active)) "inline_literal" else NULL,
    fit_population_person_data_file = as.character(fit_person_data_file %||% NULL),
    fit_population_person_replay_scope = as.character(fit_population$person_table_replay_scope %||% NULL),
    fit_population_omitted_persons = length(as.character(fit_population$omitted_persons %||% character(0))),
    fit_population_response_rows_omitted = as.integer(fit_population$response_rows_omitted %||% 0L),
    unit_prediction_posterior_basis = as.character(unit_settings$posterior_basis %||% NULL),
    unit_prediction_person_id = as.character(unit_settings$person_id %||% NULL),
    unit_prediction_population_policy = as.character(unit_settings$population_policy %||% NULL),
    unit_prediction_population_formula = as.character(unit_settings$population_formula %||% NULL),
    plausible_value_posterior_basis = as.character(pv_settings$posterior_basis %||% NULL),
    plausible_value_person_id = as.character(pv_settings$person_id %||% NULL),
    plausible_value_population_policy = as.character(pv_settings$population_policy %||% NULL),
    plausible_value_population_formula = as.character(pv_settings$population_formula %||% NULL),
    include_bundle = isTRUE(include_bundle),
    bundle_dir = as.character(bundle_dir[1]),
    bundle_prefix = as.character(bundle_prefix[1])
  ))

  out <- list(
    summary = summary_tbl,
    script = script_text,
    settings = settings,
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "build_mfrm_replay_script()",
      extra_areas = c(
        "Score-side scorefile export under bounded GPCM",
        "FACETS output-contract score-side review",
        "Design planning and forecasting"
      )
    ),
    anchors = anchor_df,
    group_anchors = group_anchor_df
  )
  as_mfrm_bundle(out, "mfrm_replay_script")
}

resolve_conquest_overlap_input <- function(fit = NULL,
                                           case = c("synthetic_latent_regression"),
                                           quad_points = 7L,
                                           maxit = 40L,
                                           reltol = 1e-6) {
  if (!is.null(fit)) {
    if (inherits(fit, "mfrm_facets_run")) {
      fit <- fit$fit
    }
    if (!inherits(fit, "mfrm_fit")) {
      stop("`fit` must be an mfrm_fit or mfrm_facets_run object.", call. = FALSE)
    }
    return(list(
      fit = fit,
      source = "fit_input",
      case = "fit_input",
      truth = NULL
    ))
  }

  case <- match.arg(case)
  fit_obj <- fit_reference_benchmark_dataset(
    dataset = case,
    method = "MML",
    model = "RSM",
    quad_points = as.integer(quad_points[1] %||% 7L),
    maxit = as.integer(maxit[1] %||% 40L),
    reltol = as.numeric(reltol[1] %||% 1e-6)
  )

  list(
    fit = fit_obj$fit,
    source = "internal_case",
    case = case,
    truth = fit_obj$truth %||% NULL
  )
}

validate_conquest_overlap_fit <- function(fit) {
  cfg <- fit$config %||% list()
  pop <- fit$population %||% list()
  method <- toupper(as.character(cfg$method_input %||% cfg$method %||% NA_character_))
  model <- toupper(as.character(cfg$model %||% NA_character_))
  posterior_basis <- as.character(cfg$posterior_basis %||% pop$posterior_basis %||% "legacy_mml")
  facet_names <- as.character(cfg$facet_names %||% character(0))

  if (!identical(method, "MML")) {
    stop("ConQuest overlap bundle currently supports only `MML` fits.", call. = FALSE)
  }
  if (!(model %in% c("RSM", "PCM"))) {
    stop(
      "ConQuest overlap bundle is currently restricted to the ordered-response `RSM` / `PCM` model scope.",
      call. = FALSE
    )
  }
  if (!isTRUE(pop$active) || !identical(posterior_basis, "population_model")) {
    stop(
      "ConQuest overlap bundle currently requires an active latent-regression `MML` fit (`posterior_basis = 'population_model'`).",
      call. = FALSE
    )
  }
  if (length(facet_names) != 1L) {
    stop(
      "ConQuest overlap bundle is currently restricted to exact-overlap item-only fits with exactly one non-person facet.",
      call. = FALSE
    )
  }

  dat <- as.data.frame(fit$prep$data %||% data.frame(), stringsAsFactors = FALSE)
  facet_name <- facet_names[1]
  required_cols <- c("Person", facet_name, "Score")
  if (!all(required_cols %in% names(dat))) {
    stop("The fitted object does not retain the canonical response columns needed for the overlap bundle.",
         call. = FALSE)
  }
  if ("Weight" %in% names(dat)) {
    wt <- suppressWarnings(as.numeric(dat$Weight))
    if (any(is.finite(wt) & abs(wt - 1) > 1e-12, na.rm = TRUE)) {
      stop("ConQuest overlap bundle currently assumes unit weights only.", call. = FALSE)
    }
  }

  score_vals <- sort(unique(suppressWarnings(as.numeric(dat$Score))))
  score_vals <- score_vals[is.finite(score_vals)]
  if (length(score_vals) != 2L) {
    stop("ConQuest overlap bundle currently requires binary responses.", call. = FALSE)
  }

  person_tbl <- as.data.frame(pop$person_table %||% data.frame(), stringsAsFactors = FALSE)
  person_id <- as.character(pop$person_id %||% "Person")
  design_cols <- as.character(pop$design_columns %||% character(0))
  covariates <- setdiff(design_cols, "(Intercept)")
  if (length(covariates) != 1L) {
    stop(
      "ConQuest overlap bundle currently requires exactly one numeric person covariate beyond the intercept.",
      call. = FALSE
    )
  }
  covariate <- covariates[1]
  if (!all(c(person_id, covariate) %in% names(person_tbl))) {
    if (length(fit$population$xlevels %||% list()) > 0L ||
        length(fit$population$contrasts %||% list()) > 0L) {
      stop(
        "ConQuest overlap bundle currently requires one raw numeric person covariate column. ",
        "Categorical/model-matrix-expanded covariates are stored for mfrmr scoring, but are not yet supported by this exact-overlap ConQuest export.",
        call. = FALSE
      )
    }
    stop("Stored person-level background data do not contain the required ConQuest-overlap covariate columns.",
         call. = FALSE)
  }
  if (!is.numeric(person_tbl[[covariate]])) {
    stop("ConQuest overlap bundle currently requires a numeric person covariate.", call. = FALSE)
  }

  person_levels <- as.character(fit$prep$levels$Person %||% unique(as.character(dat$Person)))
  item_levels <- as.character(cfg$facet_levels[[facet_name]] %||% unique(as.character(dat[[facet_name]])))

  long <- data.frame(
    Person = as.character(dat$Person),
    ItemLevel = as.character(dat[[facet_name]]),
    ScoreOriginal = suppressWarnings(as.numeric(dat$Score)),
    stringsAsFactors = FALSE
  )
  key <- paste(long$Person, long$ItemLevel, sep = "\r")
  if (anyDuplicated(key)) {
    stop("ConQuest overlap bundle currently requires exactly one response per person-by-item cell.",
         call. = FALSE)
  }

  expected_n <- length(person_levels) * length(item_levels)
  if (nrow(long) != expected_n) {
    stop(
      "ConQuest overlap bundle currently requires a complete rectangular person-by-item response matrix with no missing cells.",
      call. = FALSE
    )
  }

  score_map <- stats::setNames(c(0L, 1L), as.character(score_vals))
  long$Score <- unname(score_map[as.character(long$ScoreOriginal)])
  if (anyNA(long$Score)) {
    stop("Binary responses could not be standardized to {0, 1} for the ConQuest overlap export.",
         call. = FALSE)
  }

  person_tbl <- person_tbl[, c(person_id, covariate), drop = FALSE]
  names(person_tbl) <- c("Person", covariate)
  person_tbl <- person_tbl[match(person_levels, as.character(person_tbl$Person)), , drop = FALSE]
  if (anyNA(person_tbl$Person)) {
    stop("Stored person-level background data do not align with all fitted persons.", call. = FALSE)
  }

  list(
    fit = fit,
    config = cfg,
    population = pop,
    facet_name = facet_name,
    covariate = covariate,
    person_levels = person_levels,
    item_levels = item_levels,
    long = long,
    person_data = person_tbl
  )
}

build_conquest_overlap_response_wide <- function(long, person_levels, item_levels, covariate_data) {
  long <- long[order(match(long$Person, person_levels), match(long$ItemLevel, item_levels)), , drop = FALSE]
  wide <- stats::reshape(
    long[, c("Person", "ItemLevel", "Score"), drop = FALSE],
    idvar = "Person",
    timevar = "ItemLevel",
    direction = "wide"
  )
  wide <- wide[match(person_levels, as.character(wide$Person)), , drop = FALSE]
  response_cols <- paste0("Score.", item_levels)
  response_vars <- sprintf("I%03d", seq_along(item_levels))
  scores <- wide[, response_cols, drop = FALSE]
  names(scores) <- response_vars
  out <- cbind(covariate_data, scores, stringsAsFactors = FALSE)
  list(
    response_wide = as.data.frame(out, stringsAsFactors = FALSE),
    item_map = data.frame(
      ResponseVar = response_vars,
      OriginalLevel = item_levels,
      stringsAsFactors = FALSE
    )
  )
}

conquest_overlap_response_spec <- function(response_vars) {
  response_vars <- as.character(response_vars %||% character(0))
  if (length(response_vars) == 0L) return("")
  if (length(response_vars) == 1L) return(response_vars[1])
  paste0(response_vars[1], " to ", utils::tail(response_vars, 1))
}

build_conquest_overlap_command_template <- function(data_file,
                                                    response_vars,
                                                    covariate,
                                                    prefix = "conquest_overlap") {
  response_spec <- conquest_overlap_response_spec(response_vars)
  c(
    "/*",
    "Generated by mfrmr::build_conquest_overlap_bundle()",
    "Scope: ordered-response RSM/PCM, presently operationalized as binary item-only latent regression with one numeric person covariate.",
    "Official-manual alignment: block comments, CSV PID/keeps widths, and machine-readable export/show CSV outputs are requested.",
    "Combine exported reg_coefficients and covariance outputs into the population table before review normalization.",
    "Confirm local ConQuest syntax/options before treating this as an external benchmark run.",
    "*/",
    paste0(
      "datafile ",
      as.character(data_file[1]),
      " ! filetype=csv, columnlabels=yes, pid=Person, pidwidth=32, responses=",
      response_spec,
      ", keeps=",
      as.character(covariate[1]),
      ", keepswidth=32",
      ";"
    ),
    "score (0,1);",
    paste0("regression ", as.character(covariate[1]), ";"),
    "model item;",
    "estimate;",
    paste0("export parameters ! filetype=csv >> ", as.character(prefix[1]), "_conquest_parameters.csv;"),
    paste0("export reg_coefficients ! filetype=csv >> ", as.character(prefix[1]), "_conquest_reg_coefficients.csv;"),
    paste0("export covariance ! filetype=csv >> ", as.character(prefix[1]), "_conquest_covariance.csv;"),
    paste0("show cases ! estimates=eap, filetype=csv, regressors=yes >> ", as.character(prefix[1]), "_conquest_cases_eap.csv;"),
    paste0("show parameters ! tables=1:4, estimates=eap >> ", as.character(prefix[1]), "_conquest_parameters_review.txt;")
  )
}

build_conquest_overlap_output_contract <- function(prefix = "conquest_overlap") {
  prefix <- as.character(prefix[1] %||% "conquest_overlap")
  data.frame(
    ExternalFile = c(
      paste0(prefix, "_conquest_parameters.csv"),
      paste0(prefix, "_conquest_reg_coefficients.csv"),
      paste0(prefix, "_conquest_covariance.csv"),
      paste0(prefix, "_conquest_cases_eap.csv"),
      paste0(prefix, "_conquest_parameters_review.txt")
    ),
    ConQuestCommand = c(
      "export parameters ! filetype=csv",
      "export reg_coefficients ! filetype=csv",
      "export covariance ! filetype=csv",
      "show cases ! estimates=eap, filetype=csv, regressors=yes",
      "show parameters ! tables=1:4, estimates=eap"
    ),
    ReviewHandoff = c(
      "Use as the extracted item-estimate table after confirming item labels/parameter rows.",
      "Combine into the normalized population table as regression-coefficient rows.",
      "Combine into the normalized population table as the residual variance/covariance row.",
      "Use as the normalized case-level EAP table after selecting the EAP column.",
      "Human-readable review only; do not treat this text file as a parsed review table."
    ),
    RequiredForReview = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}

build_conquest_overlap_readme <- function(summary_tbl,
                                          comparison_targets,
                                          output_contract,
                                          command_file,
                                          data_file) {
  c(
    "mfrmr ConQuest overlap bundle",
    "",
    "Purpose:",
    "This bundle packages one exact-overlap comparison case for ConQuest and mfrmr.",
    "The supported scope is intentionally narrow: ordered-response RSM/PCM,",
    "currently operationalized as binary item-only latent regression,",
    "unidimensional MML, and one numeric person covariate.",
    "",
    "Key files:",
    paste0("- Wide CSV for ConQuest template: ", as.character(data_file[1])),
    paste0("- Command template: ", as.character(command_file[1])),
    "- mfrmr_population.csv",
    "- mfrmr_item_estimates.csv",
    "- mfrmr_case_eap.csv",
    "- item_map.csv",
    "- *_conquest_output_contract.csv",
    "- *_conquest_parameters.csv / *_conquest_reg_coefficients.csv / *_conquest_covariance.csv are requested by the generated command template",
    "- *_conquest_cases_eap.csv is requested by the generated command template",
    "",
    "Requested external ConQuest outputs:",
    paste0(
      "- ",
      output_contract$ExternalFile,
      ": ",
      output_contract$ReviewHandoff
    ),
    "",
    "Comparison rules:",
    paste0(
      "- ",
      comparison_targets$Target,
      ": ",
      comparison_targets$ComparisonRule
    ),
    "",
    "Caution:",
    "This bundle is not a claim of ConQuest numerical equivalence.",
    "Use it only where model family, dimensionality, response coding, and covariate coding match exactly.",
    "",
    "Summary:",
    paste(utils::capture.output(print(summary_tbl, row.names = FALSE)), collapse = "\n")
  )
}

#' Build a scoped ConQuest-overlap bundle
#'
#' @param fit Optional output from [fit_mfrm()] or [run_mfrm_facets()]. When
#'   omitted, the helper builds the package's
#'   `"synthetic_latent_regression"` overlap case.
#' @param case Overlap case used when `fit = NULL`. Currently only
#'   `"synthetic_latent_regression"` is supported.
#' @param output_dir Optional directory where the bundle files should be
#'   written. When `NULL`, the helper returns the in-memory bundle only.
#' @param prefix File-name prefix used when writing the bundle to disk.
#' @param overwrite If `FALSE`, refuse to overwrite existing files.
#' @param quad_points Quadrature points used when `fit = NULL` and the
#'   overlap case is fit on the fly.
#' @param maxit Maximum optimizer iterations used when `fit = NULL`.
#' @param reltol Relative convergence tolerance used when `fit = NULL`.
#'
#' @details
#' This helper prepares a narrow ConQuest comparison bundle for an `RSM` / `PCM`
#' latent-regression `MML` fit and records the `mfrmr`-side tables to compare
#' after an external ConQuest run. The supported overlap is intentionally
#' narrow:
#'
#' - ordered-response `RSM` / `PCM` only;
#' - binary responses only;
#' - exactly one non-person facet, treated as the item facet;
#' - active latent-regression `MML`;
#' - exactly one numeric person covariate beyond the intercept;
#' - complete person-by-item rectangular data.
#'
#' The returned bundle standardizes the responses to `{0, 1}`, pivots them to a
#' one-row-per-person wide CSV, stores the corresponding person covariates, and
#' records the `mfrmr` estimates that should be compared externally.
#'
#' The `conquest_command` component is a conservative starting template, not a
#' guaranteed version-invariant automation. The `conquest_output_contract`
#' component records which requested external output should feed each
#' normalized review table.
#' Use [normalize_conquest_overlap_files()] or
#' [normalize_conquest_overlap_tables()] and then [review_conquest_overlap()] only
#' after the matching ConQuest run has been executed externally and the relevant
#' output tables have been extracted. The bundle and command template alone are
#' not external validation evidence.
#'
#' @section Comparison targets:
#' - regression slope: compare directly;
#' - residual variance `sigma2`: compare directly;
#' - item estimates: compare after centering because the Rasch location origin
#'   remains constraint-dependent;
#' - case EAP estimates: compare as posterior summaries under the fitted
#'   population model.
#'
#' @section Output:
#' The returned object has class `mfrm_conquest_overlap_bundle` and includes:
#' - `summary`: one-row scope summary with posterior-basis and
#'   population-model review fields
#' - `comparison_targets`: comparison rules for the exported tables
#' - `conquest_output_contract`: requested ConQuest outputs and review handoff
#' - `response_long`: long-format binary response data used by the bundle
#' - `response_wide`: wide CSV-ready response matrix for the ConQuest template
#' - `person_data`: one-row-per-person covariate table
#' - `item_map`: mapping from exported response columns to original item levels
#' - `mfrmr_population`: fitted population-model coefficients plus `sigma2`
#' - `mfrmr_item_estimates`: fitted item estimates with centered values
#' - `mfrmr_case_eap`: posterior EAP summaries for the fitted persons
#' - `conquest_command`: conservative ConQuest command template
#' - `written_files`: file inventory when `output_dir` is supplied
#' - `settings`: bundle settings
#' - `notes`: interpretation notes
#'
#' @return A named list with class `mfrm_conquest_overlap_bundle`.
#' @seealso [normalize_conquest_overlap_files()],
#'   [normalize_conquest_overlap_tables()], [review_conquest_overlap()],
#'   [reference_case_benchmark()], [build_mfrm_replay_script()],
#'   [export_mfrm_bundle()]
#' @examples
#' \donttest{
#' bundle <- build_conquest_overlap_bundle(quad_points = 3, maxit = 30)
#' bundle$summary[, c("Case", "Facet", "Covariate", "Persons", "Items")]
#' summary(bundle)$conquest_command_scope
#' summary(bundle)$conquest_output_contract
#' cat(substr(bundle$conquest_command, 1, 120))
#' }
#' @export
build_conquest_overlap_bundle <- function(fit = NULL,
                                          case = c("synthetic_latent_regression"),
                                          output_dir = NULL,
                                          prefix = "conquest_overlap",
                                          overwrite = FALSE,
                                          quad_points = 7L,
                                          maxit = 40L,
                                          reltol = 1e-6) {
  resolved <- resolve_conquest_overlap_input(
    fit = fit,
    case = case,
    quad_points = quad_points,
    maxit = maxit,
    reltol = reltol
  )
  prepared <- validate_conquest_overlap_fit(resolved$fit)
  wide_parts <- build_conquest_overlap_response_wide(
    long = prepared$long,
    person_levels = prepared$person_levels,
    item_levels = prepared$item_levels,
    covariate_data = prepared$person_data
  )

  fit <- prepared$fit
  pop <- fit$population
  population_formula <- pop$formula %||% fit$config$population_formula %||% NULL
  population_formula_label <- if (is.null(population_formula)) {
    NA_character_
  } else {
    paste(deparse(population_formula), collapse = " ")
  }
  population_design_columns <- as.character(pop$design_columns %||% character(0))
  population_coefficients <- pop$coefficients %||% numeric(0)
  item_tbl <- as.data.frame(fit$facets$others, stringsAsFactors = FALSE)
  item_tbl <- item_tbl[item_tbl$Facet == prepared$facet_name, c("Level", "Estimate"), drop = FALSE]
  item_tbl <- item_tbl[match(prepared$item_levels, as.character(item_tbl$Level)), , drop = FALSE]
  item_tbl$ResponseVar <- wide_parts$item_map$ResponseVar
  item_tbl$CenteredEstimate <- suppressWarnings(as.numeric(item_tbl$Estimate))
  item_tbl$CenteredEstimate <- item_tbl$CenteredEstimate - mean(item_tbl$CenteredEstimate, na.rm = TRUE)

  population_tbl <- data.frame(
    Parameter = c(names(pop$coefficients %||% numeric(0)), "sigma2"),
    Estimate = c(as.numeric(pop$coefficients %||% numeric(0)), as.numeric(pop$sigma2)),
    ComparisonRule = c(
      rep("Direct comparison for numeric regression coefficients.", length(pop$coefficients %||% numeric(0))),
      "Direct comparison for residual latent variance."
    ),
    stringsAsFactors = FALSE
  )
  if ("(Intercept)" %in% population_tbl$Parameter) {
    population_tbl$ComparisonRule[population_tbl$Parameter == "(Intercept)"] <-
      "Compare only if item estimates are centered under the same location constraint."
  }

  case_eap <- predict_mfrm_units(
    fit = fit,
    new_data = fit$prep$data[, c("Person", prepared$facet_name, "Score"), drop = FALSE],
    person = "Person",
    facets = prepared$facet_name,
    score = "Score",
    person_data = prepared$person_data,
    person_id = "Person",
    n_draws = 0
  )$estimates
  case_eap <- as.data.frame(case_eap, stringsAsFactors = FALSE)
  case_eap <- case_eap[match(prepared$person_levels, as.character(case_eap$Person)), , drop = FALSE]
  case_eap[[prepared$covariate]] <- prepared$person_data[[prepared$covariate]][match(case_eap$Person, prepared$person_data$Person)]

  comparison_targets <- data.frame(
    Target = c(
      paste0("Population coefficient: ", prepared$covariate),
      "Population coefficient: (Intercept)",
      "Population variance: sigma2",
      paste0("Item estimates for facet `", prepared$facet_name, "`"),
      "Case EAP estimates"
    ),
    MfrmrComponent = c(
      "mfrmr_population",
      "mfrmr_population",
      "mfrmr_population",
      "mfrmr_item_estimates",
      "mfrmr_case_eap"
    ),
    ComparisonRule = c(
      "Compare directly.",
      "Compare only if the external calibration uses the same item-centering/location convention.",
      "Compare directly.",
      "Compare after centering the item estimates.",
      "Compare as posterior EAP summaries under the fitted population model."
    ),
    stringsAsFactors = FALSE
  )
  output_contract <- build_conquest_overlap_output_contract(prefix = prefix)

  wide_file <- paste0(prefix, "_wide.csv")
  command_file <- paste0(prefix, ".cqc")
  command_text <- paste(
    build_conquest_overlap_command_template(
      data_file = wide_file,
      response_vars = wide_parts$item_map$ResponseVar,
      covariate = prepared$covariate,
      prefix = prefix
    ),
    collapse = "\n"
  )

  written_files <- data.frame(
    Component = character(0),
    Format = character(0),
    Path = character(0),
    stringsAsFactors = FALSE
  )
  add_written <- function(component, format, path) {
    written_files <<- rbind(
      written_files,
      data.frame(Component = component, Format = format, Path = path, stringsAsFactors = FALSE)
    )
  }

  if (!is.null(output_dir)) {
    output_dir <- normalizePath(as.character(output_dir[1]), winslash = "/", mustWork = FALSE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    overwrite <- isTRUE(overwrite)

    write_csv <- function(df, file, component) {
      path <- file.path(output_dir, file)
      if (file.exists(path) && !overwrite) {
        stop("File already exists: ", path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
      }
      utils::write.csv(df, path, row.names = FALSE, na = "")
      add_written(component, "csv", path)
    }
    write_text <- function(text, file, component) {
      path <- file.path(output_dir, file)
      if (file.exists(path) && !overwrite) {
        stop("File already exists: ", path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
      }
      writeLines(enc2utf8(as.character(text)), con = path, useBytes = TRUE)
      add_written(component, tools::file_ext(path), path)
    }

    write_csv(prepared$long[, c("Person", "ItemLevel", "Score"), drop = FALSE], paste0(prefix, "_long.csv"), "response_long")
    write_csv(wide_parts$response_wide, wide_file, "response_wide")
    write_csv(prepared$person_data, paste0(prefix, "_person_data.csv"), "person_data")
    write_csv(wide_parts$item_map, paste0(prefix, "_item_map.csv"), "item_map")
    write_csv(population_tbl, paste0(prefix, "_mfrmr_population.csv"), "mfrmr_population")
    write_csv(item_tbl, paste0(prefix, "_mfrmr_item_estimates.csv"), "mfrmr_item_estimates")
    write_csv(case_eap, paste0(prefix, "_mfrmr_case_eap.csv"), "mfrmr_case_eap")
    write_csv(comparison_targets, paste0(prefix, "_comparison_targets.csv"), "comparison_targets")
    write_csv(output_contract, paste0(prefix, "_conquest_output_contract.csv"), "conquest_output_contract")
    write_text(command_text, command_file, "conquest_command")
    write_text(
      build_conquest_overlap_readme(
        summary_tbl = data.frame(
          Case = resolved$case,
          Facet = prepared$facet_name,
          Covariate = prepared$covariate,
          Persons = length(prepared$person_levels),
          Items = length(prepared$item_levels),
          PosteriorBasis = as.character(fit$config$posterior_basis %||% pop$posterior_basis %||% "population_model"),
          PopulationFormula = population_formula_label,
          PopulationDesignColumns = paste(population_design_columns, collapse = ", "),
          PopulationResidualVariance = suppressWarnings(as.numeric(pop$sigma2 %||% NA_real_)),
          PopulationOmittedPersons = as.integer(length(pop$omitted_persons %||% character(0))),
          PopulationResponseRowsOmitted = suppressWarnings(as.integer(pop$response_rows_omitted %||% NA_integer_)),
          stringsAsFactors = FALSE
        ),
        comparison_targets = comparison_targets,
        output_contract = output_contract,
        command_file = command_file,
        data_file = wide_file
      ),
      paste0(prefix, "_README.txt"),
      "readme"
    )
  }

  summary_tbl <- data.frame(
    Case = resolved$case,
    InputSource = resolved$source,
    Facet = prepared$facet_name,
    Covariate = prepared$covariate,
    Persons = length(prepared$person_levels),
    Items = length(prepared$item_levels),
    RowsLong = nrow(prepared$long),
    Method = as.character(fit$config$method_input %||% fit$config$method %||% NA_character_),
    Model = as.character(fit$config$model %||% NA_character_),
    PosteriorBasis = as.character(fit$config$posterior_basis %||% pop$posterior_basis %||% "population_model"),
    PopulationFormula = population_formula_label,
    PopulationDesignColumns = paste(population_design_columns, collapse = ", "),
    PopulationCoefficientCount = as.integer(length(population_coefficients)),
    PopulationResidualVariance = suppressWarnings(as.numeric(pop$sigma2 %||% NA_real_)),
    PopulationIncludedPersons = as.integer(length(pop$included_persons %||% character(0))),
    PopulationOmittedPersons = as.integer(length(pop$omitted_persons %||% character(0))),
    PopulationResponseRowsRetained = suppressWarnings(as.integer(pop$response_rows_retained %||% NA_integer_)),
    PopulationResponseRowsOmitted = suppressWarnings(as.integer(pop$response_rows_omitted %||% NA_integer_)),
    FilesWritten = nrow(written_files),
    stringsAsFactors = FALSE
  )

  notes <- c(
    "This bundle is intentionally restricted to a binary item-only latent-regression overlap case.",
    "Use the ConQuest command text as a conservative template, not as a claim of guaranteed version-invariant automation.",
    "Item estimates should be centered before external comparison."
  )

  settings <- dashboard_settings_table(list(
    case = resolved$case,
    input_source = resolved$source,
    facet = prepared$facet_name,
    covariate = prepared$covariate,
    output_dir = if (!is.null(output_dir)) output_dir else NA_character_,
    prefix = as.character(prefix[1]),
    overwrite = isTRUE(overwrite)
  ))

  out <- list(
    summary = summary_tbl,
    comparison_targets = comparison_targets,
    conquest_output_contract = output_contract,
    response_long = prepared$long[, c("Person", "ItemLevel", "Score"), drop = FALSE],
    response_wide = wide_parts$response_wide,
    person_data = prepared$person_data,
    item_map = wide_parts$item_map,
    mfrmr_population = population_tbl,
    mfrmr_item_estimates = item_tbl,
    mfrmr_case_eap = case_eap,
    conquest_command = command_text,
    written_files = written_files,
    settings = settings,
    notes = notes
  )
  as_mfrm_bundle(out, "mfrm_conquest_overlap_bundle")
}

validate_conquest_overlap_bundle_object <- function(bundle) {
  if (!inherits(bundle, "mfrm_conquest_overlap_bundle")) {
    stop("`bundle` must be output from build_conquest_overlap_bundle().", call. = FALSE)
  }
  bundle
}

validate_conquest_overlap_input_df <- function(x, arg) {
  if (!is.data.frame(x)) {
    stop("`", arg, "` must be a data.frame containing normalized ConQuest output columns.",
         call. = FALSE)
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}

normalize_conquest_overlap_column_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

conquest_overlap_column_aliases <- function(arg) {
  switch(
    as.character(arg[1]),
    conquest_population_term = c("Parameter", "Term", "Label"),
    conquest_population_estimate = c("Estimate", "Est", "Value"),
    conquest_item_id = c("ResponseVar", "ItemID", "Item", "Label", "OriginalLevel"),
    conquest_item_estimate = c("Estimate", "Est", "Facility"),
    conquest_case_person = c("Person", "PID", "Sequence ID", "SequenceID"),
    conquest_case_estimate = c("Estimate", "EAP_1", "EAP"),
    character(0)
  )
}

resolve_conquest_overlap_column <- function(df, column, arg, role) {
  column <- as.character(column[1] %||% NA_character_)
  if (!nzchar(column) || is.na(column)) {
    stop("`", arg, "` must name an existing column in `", role, "` or use \"auto\".", call. = FALSE)
  }
  if (!identical(column, "auto")) {
    if (!column %in% names(df)) {
      stop("`", arg, "` must name an existing column in `", role, "`.", call. = FALSE)
    }
    return(column)
  }

  aliases <- conquest_overlap_column_aliases(arg)
  if (length(aliases) == 0L) {
    stop("`", arg, "` does not support automatic alias resolution.", call. = FALSE)
  }

  norm_names <- normalize_conquest_overlap_column_key(names(df))
  norm_aliases <- normalize_conquest_overlap_column_key(aliases)
  matched_idx <- which(norm_names %in% norm_aliases)

  if (length(matched_idx) == 0L) {
    stop(
      "`", arg, "` could not be resolved automatically in `", role, "`. ",
      "Supply the column name explicitly.",
      call. = FALSE
    )
  }
  if (length(matched_idx) > 1L) {
    stop(
      "`", arg, "` matched multiple columns in `", role, "` under automatic alias resolution: ",
      paste(names(df)[matched_idx], collapse = ", "),
      ". Supply the column name explicitly.",
      call. = FALSE
    )
  }

  names(df)[matched_idx]
}

normalize_conquest_overlap_table_component <- function(df,
                                                       id_col,
                                                       estimate_col,
                                                       id_name,
                                                       keep_extra_columns = TRUE) {
  estimate_numeric <- suppressWarnings(as.numeric(df[[estimate_col]]))
  out <- data.frame(
    Identifier = as.character(df[[id_col]]),
    Estimate = estimate_numeric,
    EstimateNonNumeric = flag_conquest_overlap_non_numeric(df[[estimate_col]], estimate_numeric),
    stringsAsFactors = FALSE
  )
  names(out)[1] <- id_name

  if (isTRUE(keep_extra_columns)) {
    extra_cols <- setdiff(names(df), c(id_col, estimate_col))
    extras <- df[, extra_cols, drop = FALSE]
    extras <- extras[, setdiff(names(extras), c(id_name, "Estimate", "EstimateNonNumeric")), drop = FALSE]
    if (ncol(extras) > 0L) {
      out <- cbind(out, extras, stringsAsFactors = FALSE)
    }
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}

count_conquest_overlap_non_numeric <- function(raw_values, numeric_values) {
  sum(flag_conquest_overlap_non_numeric(raw_values, numeric_values))
}

flag_conquest_overlap_non_numeric <- function(raw_values, numeric_values) {
  raw_chr <- trimws(as.character(raw_values))
  nzchar(raw_chr) & !is.na(raw_chr) & is.na(numeric_values)
}

resolve_conquest_overlap_non_numeric_flag <- function(df, raw_values, numeric_values) {
  if ("EstimateNonNumeric" %in% names(df)) {
    out <- suppressWarnings(as.logical(df$EstimateNonNumeric))
    out[is.na(out)] <- FALSE
    return(out)
  }
  flag_conquest_overlap_non_numeric(raw_values, numeric_values)
}

resolve_conquest_overlap_delimiter <- function(path, delimiter = c("auto", "comma", "tab", "semicolon", ",", "\t", ";")) {
  delimiter <- match.arg(delimiter)
  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }
  if (!identical(delimiter, "auto")) {
    return(
      switch(
        delimiter,
        comma = ",",
        tab = "\t",
        semicolon = ";",
        delimiter
      )
    )
  }

  ext <- tolower(tools::file_ext(path))
  first_line <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  first_line <- first_line[1] %||% ""

  if (grepl("\t", first_line, fixed = TRUE) || ext %in% c("tsv")) return("\t")
  if (grepl(";", first_line, fixed = TRUE) && !grepl(",", first_line, fixed = TRUE)) return(";")
  if (ext %in% c("txt") && grepl(";", first_line, fixed = TRUE) && !grepl(",", first_line, fixed = TRUE)) return(";")
  ","
}

read_conquest_overlap_file <- function(path,
                                       role,
                                       delimiter = c("auto", "comma", "tab", "semicolon", ",", "\t", ";")) {
  path <- normalizePath(as.character(path[1]), winslash = "/", mustWork = FALSE)
  if (!file.exists(path)) {
    stop("`", role, "` does not exist: ", path, call. = FALSE)
  }
  sep <- resolve_conquest_overlap_delimiter(path, delimiter = delimiter)
  utils::read.csv(
    file = path,
    sep = sep,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

coerce_conquest_overlap_tables_for_audit <- function(conquest_population,
                                                     conquest_item_estimates,
                                                     conquest_case_eap) {
  if (inherits(conquest_population, "mfrm_conquest_overlap_tables")) {
    if (!is.null(conquest_item_estimates) || !is.null(conquest_case_eap)) {
      stop(
        "When `conquest_population` is output from normalize_conquest_overlap_tables(), do not also supply `conquest_item_estimates` or `conquest_case_eap` separately.",
        call. = FALSE
      )
    }
    return(list(
      conquest_population = as.data.frame(conquest_population$conquest_population, stringsAsFactors = FALSE),
      conquest_item_estimates = as.data.frame(conquest_population$conquest_item_estimates, stringsAsFactors = FALSE),
      conquest_case_eap = as.data.frame(conquest_population$conquest_case_eap, stringsAsFactors = FALSE),
      standardized = TRUE
    ))
  }

  if (is.null(conquest_population) || is.null(conquest_item_estimates) || is.null(conquest_case_eap)) {
    stop(
      "Supply either three extracted ConQuest tables or a single object from normalize_conquest_overlap_tables().",
      call. = FALSE
    )
  }

  list(
    conquest_population = conquest_population,
    conquest_item_estimates = conquest_item_estimates,
    conquest_case_eap = conquest_case_eap,
    standardized = FALSE
  )
}

choose_conquest_item_match <- function(bundle, conquest_ids, item_id_source = c("auto", "response_var", "level")) {
  item_id_source <- match.arg(item_id_source)
  conquest_ids <- as.character(conquest_ids)
  item_map <- as.data.frame(bundle$item_map, stringsAsFactors = FALSE)
  response_overlap <- sum(as.character(item_map$ResponseVar) %in% conquest_ids)
  level_overlap <- sum(as.character(item_map$OriginalLevel) %in% conquest_ids)

  chosen <- item_id_source
  if (identical(chosen, "auto")) {
    chosen <- if (response_overlap >= level_overlap) "response_var" else "level"
  }
  list(
    source = chosen,
    match_ids = if (identical(chosen, "response_var")) {
      as.character(item_map$ResponseVar)
    } else {
      as.character(item_map$OriginalLevel)
    }
  )
}

build_conquest_overlap_attention <- function(section, id_values, issue, detail) {
  if (length(id_values) == 0L) {
    return(data.frame(
      Section = character(0),
      ID = character(0),
      Issue = character(0),
      Detail = character(0),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    Section = rep(as.character(section[1]), length(id_values)),
    ID = as.character(id_values),
    Issue = rep(as.character(issue[1]), length(id_values)),
    Detail = rep(as.character(detail[1]), length(id_values)),
    stringsAsFactors = FALSE
  )
}

count_conquest_overlap_attention <- function(attention, prefix) {
  if (!is.data.frame(attention) || nrow(attention) == 0L || !"Issue" %in% names(attention)) {
    return(0L)
  }
  as.integer(sum(startsWith(as.character(attention$Issue), as.character(prefix[1])), na.rm = TRUE))
}

max_abs_conquest_overlap_difference <- function(x) {
  x <- abs(suppressWarnings(as.numeric(x)))
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

max_abs_conquest_overlap_label <- function(x, label) {
  x <- abs(suppressWarnings(as.numeric(x)))
  keep <- is.finite(x)
  if (!any(keep)) {
    return(NA_character_)
  }
  idx <- which(keep)[which.max(x[keep])]
  out <- as.character(label)[idx]
  if (length(out) == 0L || is.na(out) || !nzchar(out)) NA_character_ else out
}

#' Normalize extracted ConQuest overlap tables to the `mfrmr` review contract
#'
#' @param conquest_population Extracted ConQuest population-parameter table as a
#'   data.frame.
#' @param conquest_item_estimates Extracted ConQuest item-estimate table as a
#'   data.frame.
#' @param conquest_case_eap Extracted ConQuest case-level EAP table as a
#'   data.frame.
#' @param conquest_population_term Column in `conquest_population` that stores
#'   parameter names. `"auto"` tries conservative aliases such as `Parameter`
#'   and `Term`.
#' @param conquest_population_estimate Column in `conquest_population` that
#'   stores parameter estimates. `"auto"` tries aliases such as `Estimate` and
#'   `Est`.
#' @param conquest_item_id Column in `conquest_item_estimates` that stores the
#'   item identifier as exported or extracted by the user. `"auto"` tries
#'   aliases such as `ResponseVar`, `ItemID`, `Item`, and `Label`.
#' @param conquest_item_estimate Column in `conquest_item_estimates` that stores
#'   item estimates. `"auto"` tries aliases such as `Estimate`, `Est`, and
#'   `Facility`.
#' @param conquest_case_person Column in `conquest_case_eap` that stores person
#'   IDs. `"auto"` tries conservative aliases such as `Person`, `PID`, and
#'   `Sequence ID`.
#' @param conquest_case_estimate Column in `conquest_case_eap` that stores case
#'   EAP estimates. `"auto"` tries conservative aliases such as `Estimate`,
#'   `EAP_1`, and `EAP`.
#' @param keep_extra_columns If `TRUE`, keep all remaining columns after the
#'   standardized identifier and estimate columns.
#'
#' @details
#' This helper does not parse raw ConQuest text output. It standardizes already
#' extracted tables to the contract used by [review_conquest_overlap()]:
#'
#' - population parameters become columns `Parameter`, `Estimate`, and
#'   `EstimateNonNumeric`;
#' - item estimates become columns `ItemID`, `Estimate`, and
#'   `EstimateNonNumeric`;
#' - case summaries become columns `Person`, `Estimate`, and
#'   `EstimateNonNumeric`.
#'
#' The resulting object is intentionally conservative. It does not infer
#' whether item IDs correspond to exported response variables or original item
#' levels; that matching step remains part of [review_conquest_overlap()], where
#' the standardized ConQuest tables are compared against a concrete overlap
#' bundle.
#'
#' @section Output:
#' The returned object has class `mfrm_conquest_overlap_tables` and includes:
#' - `summary`: one-row normalization summary
#' - `conquest_population`: standardized population table
#' - `conquest_item_estimates`: standardized item table
#' - `conquest_case_eap`: standardized case table
#' - `settings`: source-column metadata
#' - `notes`: interpretation notes
#'
#' Read `summary(normalized)$normalization_scope` before review to confirm
#' that the object contains extracted tabular inputs, not parsed raw ConQuest
#' report text, and to check duplicate-ID / non-numeric-estimate pre-review
#' flags.
#'
#' @return A named list with class `mfrm_conquest_overlap_tables`.
#' @seealso [build_conquest_overlap_bundle()], [review_conquest_overlap()]
#' @examples
#' normalized <- normalize_conquest_overlap_tables(
#'   conquest_population = data.frame(
#'     Term = c("(Intercept)", "GroupB", "sigma2"),
#'     Est = c(0, 0.2, 1)
#'   ),
#'   conquest_item_estimates = data.frame(
#'     Item = c("I1", "I2"),
#'     Est = c(-0.2, 0.2)
#'   ),
#'   conquest_case_eap = data.frame(
#'     PID = c("P001", "P002"),
#'     EAP = c(-0.1, 0.1)
#'   ),
#'   conquest_population_term = "Term",
#'   conquest_population_estimate = "Est",
#'   conquest_item_id = "Item",
#'   conquest_item_estimate = "Est",
#'   conquest_case_person = "PID",
#'   conquest_case_estimate = "EAP"
#' )
#' summary(normalized)$normalization_scope
#' @export
normalize_conquest_overlap_tables <- function(conquest_population,
                                              conquest_item_estimates,
                                              conquest_case_eap,
                                              conquest_population_term = "auto",
                                              conquest_population_estimate = "auto",
                                              conquest_item_id = "auto",
                                              conquest_item_estimate = "auto",
                                              conquest_case_person = "auto",
                                              conquest_case_estimate = "auto",
                                              keep_extra_columns = TRUE) {
  conquest_population <- validate_conquest_overlap_input_df(conquest_population, "conquest_population")
  conquest_item_estimates <- validate_conquest_overlap_input_df(conquest_item_estimates, "conquest_item_estimates")
  conquest_case_eap <- validate_conquest_overlap_input_df(conquest_case_eap, "conquest_case_eap")

  pop_term_col <- resolve_conquest_overlap_column(conquest_population, conquest_population_term, "conquest_population_term", "conquest_population")
  pop_est_col <- resolve_conquest_overlap_column(conquest_population, conquest_population_estimate, "conquest_population_estimate", "conquest_population")
  item_id_col <- resolve_conquest_overlap_column(conquest_item_estimates, conquest_item_id, "conquest_item_id", "conquest_item_estimates")
  item_est_col <- resolve_conquest_overlap_column(conquest_item_estimates, conquest_item_estimate, "conquest_item_estimate", "conquest_item_estimates")
  case_person_col <- resolve_conquest_overlap_column(conquest_case_eap, conquest_case_person, "conquest_case_person", "conquest_case_eap")
  case_est_col <- resolve_conquest_overlap_column(conquest_case_eap, conquest_case_estimate, "conquest_case_estimate", "conquest_case_eap")

  pop_tbl <- normalize_conquest_overlap_table_component(
    conquest_population,
    id_col = pop_term_col,
    estimate_col = pop_est_col,
    id_name = "Parameter",
    keep_extra_columns = keep_extra_columns
  )
  item_tbl <- normalize_conquest_overlap_table_component(
    conquest_item_estimates,
    id_col = item_id_col,
    estimate_col = item_est_col,
    id_name = "ItemID",
    keep_extra_columns = keep_extra_columns
  )
  case_tbl <- normalize_conquest_overlap_table_component(
    conquest_case_eap,
    id_col = case_person_col,
    estimate_col = case_est_col,
    id_name = "Person",
    keep_extra_columns = keep_extra_columns
  )

  summary_tbl <- data.frame(
    PopulationRows = nrow(pop_tbl),
    PopulationDuplicateIDs = sum(duplicated(pop_tbl$Parameter)),
    PopulationNonNumeric = count_conquest_overlap_non_numeric(conquest_population[[pop_est_col]], pop_tbl$Estimate),
    ItemRows = nrow(item_tbl),
    ItemDuplicateIDs = sum(duplicated(item_tbl$ItemID)),
    ItemNonNumeric = count_conquest_overlap_non_numeric(conquest_item_estimates[[item_est_col]], item_tbl$Estimate),
    CaseRows = nrow(case_tbl),
    CaseDuplicateIDs = sum(duplicated(case_tbl$Person)),
    CaseNonNumeric = count_conquest_overlap_non_numeric(conquest_case_eap[[case_est_col]], case_tbl$Estimate),
    KeepExtraColumns = isTRUE(keep_extra_columns),
    stringsAsFactors = FALSE
  )

  settings <- dashboard_settings_table(list(
    conquest_population_term = pop_term_col,
    conquest_population_estimate = pop_est_col,
    conquest_item_id = item_id_col,
    conquest_item_estimate = item_est_col,
    conquest_case_person = case_person_col,
    conquest_case_estimate = case_est_col,
    keep_extra_columns = isTRUE(keep_extra_columns)
  ))

  notes <- c(
    "This helper standardizes extracted ConQuest tables but does not parse raw ConQuest text output.",
    "Item identifiers remain user-supplied labels until review_conquest_overlap() matches them against the exported overlap bundle.",
    "Non-numeric estimate cells are converted to NA and counted in the summary table."
  )

  out <- list(
    summary = summary_tbl,
    conquest_population = pop_tbl,
    conquest_item_estimates = item_tbl,
    conquest_case_eap = case_tbl,
    settings = settings,
    notes = notes
  )
  as_mfrm_bundle(out, "mfrm_conquest_overlap_tables")
}

#' Normalize extracted ConQuest overlap files to the `mfrmr` review contract
#'
#' @param population_file Path to an extracted ConQuest population-parameter
#'   table in CSV/TSV/TXT form.
#' @param item_file Path to an extracted ConQuest item-estimate table in
#'   CSV/TSV/TXT form.
#' @param case_file Path to an extracted ConQuest case-level EAP table in
#'   CSV/TSV/TXT form.
#' @param population_delimiter Delimiter for `population_file`. `"auto"`
#'   chooses comma, tab, or semicolon from the file extension/header line.
#' @param item_delimiter Delimiter for `item_file`. `"auto"` chooses from the
#'   file extension/header line.
#' @param case_delimiter Delimiter for `case_file`. `"auto"` chooses from the
#'   file extension/header line.
#' @param conquest_population_term Column in `population_file` that stores
#'   parameter names. `"auto"` tries conservative aliases such as `Parameter`
#'   and `Term`.
#' @param conquest_population_estimate Column in `population_file` that stores
#'   parameter estimates. `"auto"` tries aliases such as `Estimate` and `Est`.
#' @param conquest_item_id Column in `item_file` that stores the item
#'   identifier as extracted by the user. `"auto"` tries aliases such as
#'   `ResponseVar`, `ItemID`, `Item`, and `Label`.
#' @param conquest_item_estimate Column in `item_file` that stores item
#'   estimates. `"auto"` tries aliases such as `Estimate`, `Est`, and
#'   `Facility`.
#' @param conquest_case_person Column in `case_file` that stores person IDs.
#'   `"auto"` tries conservative aliases such as `Person`, `PID`, and
#'   `Sequence ID`.
#' @param conquest_case_estimate Column in `case_file` that stores case EAP
#'   estimates. `"auto"` tries conservative aliases such as `Estimate`,
#'   `EAP_1`, and `EAP`.
#' @param keep_extra_columns If `TRUE`, keep all remaining columns after the
#'   standardized identifier and estimate columns.
#'
#' @details
#' This helper is a thin file-wrapper around [normalize_conquest_overlap_tables()].
#' It is intentionally limited to already extracted tabular files and does not
#' parse raw ConQuest report text.
#'
#' The recommended workflow is:
#'
#' 1. export an exact-overlap bundle with [build_conquest_overlap_bundle()];
#' 2. extract the relevant ConQuest tables to CSV/TSV/TXT files;
#' 3. call `normalize_conquest_overlap_files()` on those files;
#' 4. pass the result to [review_conquest_overlap()].
#'
#' Read `summary(normalized)$normalization_scope` before review to confirm
#' that the files were treated as extracted tables, not raw ConQuest report
#' text, and to check duplicate-ID / non-numeric-estimate pre-review flags.
#'
#' @return A named list with class `mfrm_conquest_overlap_tables`.
#' @seealso [normalize_conquest_overlap_tables()], [review_conquest_overlap()]
#' @examples
#' \donttest{
#' bundle <- build_conquest_overlap_bundle()
#' tmp_dir <- tempdir()
#' pop_path <- file.path(tmp_dir, "cq_pop.csv")
#' item_path <- file.path(tmp_dir, "cq_item.tsv")
#' case_path <- file.path(tmp_dir, "cq_case.csv")
#' utils::write.csv(
#'   data.frame(
#'     Term = bundle$mfrmr_population$Parameter,
#'     Est = bundle$mfrmr_population$Estimate
#'   ),
#'   pop_path,
#'   row.names = FALSE
#' )
#' utils::write.table(
#'   data.frame(
#'     Item = bundle$mfrmr_item_estimates$ResponseVar,
#'     Est = bundle$mfrmr_item_estimates$Estimate
#'   ),
#'   item_path,
#'   sep = "\t",
#'   row.names = FALSE
#' )
#' utils::write.csv(
#'   data.frame(
#'     PID = bundle$mfrmr_case_eap$Person,
#'     EAP = bundle$mfrmr_case_eap$Estimate
#'   ),
#'   case_path,
#'   row.names = FALSE
#' )
#' normalized <- normalize_conquest_overlap_files(
#'   population_file = pop_path,
#'   item_file = item_path,
#'   case_file = case_path,
#'   conquest_population_term = "Term",
#'   conquest_population_estimate = "Est",
#'   conquest_item_id = "Item",
#'   conquest_item_estimate = "Est",
#'   conquest_case_person = "PID",
#'   conquest_case_estimate = "EAP"
#' )
#' summary(normalized)$normalization_scope
#' review <- review_conquest_overlap(bundle, normalized)
#' summary(review)$summary
#' }
#' @export
normalize_conquest_overlap_files <- function(population_file,
                                             item_file,
                                             case_file,
                                             population_delimiter = c("auto", "comma", "tab", "semicolon", ",", "\t", ";"),
                                             item_delimiter = c("auto", "comma", "tab", "semicolon", ",", "\t", ";"),
                                             case_delimiter = c("auto", "comma", "tab", "semicolon", ",", "\t", ";"),
                                             conquest_population_term = "auto",
                                             conquest_population_estimate = "auto",
                                             conquest_item_id = "auto",
                                             conquest_item_estimate = "auto",
                                             conquest_case_person = "auto",
                                             conquest_case_estimate = "auto",
                                             keep_extra_columns = TRUE) {
  pop_sep <- resolve_conquest_overlap_delimiter(population_file, delimiter = population_delimiter)
  item_sep <- resolve_conquest_overlap_delimiter(item_file, delimiter = item_delimiter)
  case_sep <- resolve_conquest_overlap_delimiter(case_file, delimiter = case_delimiter)

  normalized <- normalize_conquest_overlap_tables(
    conquest_population = read_conquest_overlap_file(population_file, "population_file", delimiter = population_delimiter),
    conquest_item_estimates = read_conquest_overlap_file(item_file, "item_file", delimiter = item_delimiter),
    conquest_case_eap = read_conquest_overlap_file(case_file, "case_file", delimiter = case_delimiter),
    conquest_population_term = conquest_population_term,
    conquest_population_estimate = conquest_population_estimate,
    conquest_item_id = conquest_item_id,
    conquest_item_estimate = conquest_item_estimate,
    conquest_case_person = conquest_case_person,
    conquest_case_estimate = conquest_case_estimate,
    keep_extra_columns = keep_extra_columns
  )

  normalized$source_files <- data.frame(
    Role = c("population", "items", "cases"),
    Path = c(
      normalizePath(as.character(population_file[1]), winslash = "/", mustWork = FALSE),
      normalizePath(as.character(item_file[1]), winslash = "/", mustWork = FALSE),
      normalizePath(as.character(case_file[1]), winslash = "/", mustWork = FALSE)
    ),
    Delimiter = c(pop_sep, item_sep, case_sep),
    stringsAsFactors = FALSE
  )
  normalized$settings <- rbind(
    normalized$settings,
    dashboard_settings_table(list(
      population_file = normalized$source_files$Path[1],
      item_file = normalized$source_files$Path[2],
      case_file = normalized$source_files$Path[3],
      population_delimiter = pop_sep,
      item_delimiter = item_sep,
      case_delimiter = case_sep
    ))
  )
  normalized$notes <- c(
    normalized$notes,
    "normalize_conquest_overlap_files() reads extracted CSV/TSV/TXT tables only; it does not parse raw ConQuest report text."
  )
  normalized
}

#' Review an exact-overlap ConQuest comparison against an `mfrmr` overlap bundle
#'
#' @param bundle Output from [build_conquest_overlap_bundle()].
#' @param conquest_population Normalized ConQuest population-parameter table as a
#'   data.frame, or output from [normalize_conquest_overlap_tables()].
#' @param conquest_item_estimates Normalized ConQuest item-estimate table as a
#'   data.frame. Leave `NULL` when `conquest_population` is an object from
#'   [normalize_conquest_overlap_tables()].
#' @param conquest_case_eap Normalized ConQuest case-level EAP table as a
#'   data.frame. Leave `NULL` when `conquest_population` is an object from
#'   [normalize_conquest_overlap_tables()].
#' @param conquest_population_term Column in `conquest_population` that stores
#'   parameter names. `"auto"` tries conservative aliases such as `Parameter`
#'   and `Term`.
#' @param conquest_population_estimate Column in `conquest_population` that
#'   stores parameter estimates. `"auto"` tries aliases such as `Estimate` and
#'   `Est`.
#' @param conquest_item_id Column in `conquest_item_estimates` that stores the
#'   item identifier. This may be the exported response variable (for example
#'   `I001`) or the original item/facet level. `"auto"` tries aliases such as
#'   `ResponseVar`, `ItemID`, `Item`, and `Label`.
#' @param conquest_item_estimate Column in `conquest_item_estimates` that stores
#'   the item estimate. `"auto"` tries aliases such as `Estimate`, `Est`, and
#'   `Facility`.
#' @param item_id_source How `conquest_item_id` should be matched. `"auto"`
#'   chooses the larger overlap between exported response variables and original
#'   item levels, with ties resolved toward exported response variables.
#' @param conquest_case_person Column in `conquest_case_eap` that stores person
#'   IDs. `"auto"` tries conservative aliases such as `Person`, `PID`, and
#'   `Sequence ID`.
#' @param conquest_case_estimate Column in `conquest_case_eap` that stores case
#'   EAP estimates. `"auto"` tries conservative aliases such as `Estimate`,
#'   `EAP_1`, and `EAP`.
#'
#' @details
#' This helper compares normalized ConQuest output tables against the exact-
#' overlap bundle produced by [build_conquest_overlap_bundle()]. It is
#' intentionally conservative:
#'
#' - it does **not** parse raw ConQuest text output automatically;
#' - it expects already normalized data frames or output from
#'   [normalize_conquest_overlap_tables()];
#' - and it reports numerical differences and missing elements without claiming
#'   that any fixed tolerance implies software equivalence.
#'
#' This is the package's external-table review path. It is distinct from
#' `reference_case_benchmark(cases = "synthetic_conquest_overlap_dry_run")`,
#' which only round-trips package-native tables through the same normalization
#' and review contract without executing ConQuest.
#'
#' The intended workflow is:
#'
#' 1. export an exact-overlap bundle with [build_conquest_overlap_bundle()];
#' 2. run the narrow matching case in ConQuest;
#' 3. normalize the resulting ConQuest outputs into data frames;
#' 4. pass those tables here to inspect direct differences, centered item
#'    agreement, and case-level EAP agreement.
#'
#' @section Output:
#' The returned object has class `mfrm_conquest_overlap_review` and includes:
#' - `overall`: one-row comparison summary with missing/duplicate/non-numeric
#'   attention-item counts and worst-row labels
#' - `population_comparison`: parameter-by-parameter comparison table
#' - `item_comparison`: centered item-estimate comparison table
#' - `case_comparison`: case-level EAP comparison table
#' - `attention_items`: missing, malformed, or unmatched elements
#' - `settings`: review settings
#' - `notes`: interpretation notes
#'
#' @section Interpretation:
#' - Read `summary(review)$review_scope` first to confirm that the result is a
#'   supplied-table review, not raw ConQuest text parsing or a software-
#'   equivalence claim.
#' - Population slopes and `sigma2` are intended for direct comparison.
#' - Item estimates should be interpreted after centering.
#' - Case estimates should be interpreted as posterior EAP summaries under the
#'   fitted population model.
#' - The `overall` table reports both mean and maximum absolute differences for
#'   compared population, centered item, and case rows. The
#'   `PopulationMaxAbsParameter`, `ItemCenteredMaxAbsItem`, and
#'   `CaseMaxAbsPerson` columns identify the row where each maximum absolute
#'   difference occurs.
#' - Missing or non-numeric rows in `attention_items` indicate that the external
#'   tables do not yet align cleanly with the exported overlap bundle.
#'
#' @return A named list with class `mfrm_conquest_overlap_review`.
#' @seealso [build_conquest_overlap_bundle()],
#'   [normalize_conquest_overlap_files()], [normalize_conquest_overlap_tables()],
#'   [reference_case_benchmark()]
#' @examples
#' \donttest{
#' bundle <- build_conquest_overlap_bundle()
#' raw_pop <- data.frame(
#'   Term = bundle$mfrmr_population$Parameter,
#'   Est = bundle$mfrmr_population$Estimate
#' )
#' raw_item <- data.frame(
#'   Item = bundle$mfrmr_item_estimates$ResponseVar,
#'   Est = bundle$mfrmr_item_estimates$Estimate
#' )
#' raw_case <- data.frame(
#'   PID = bundle$mfrmr_case_eap$Person,
#'   EAP = bundle$mfrmr_case_eap$Estimate
#' )
#' normalized <- normalize_conquest_overlap_tables(
#'   conquest_population = raw_pop,
#'   conquest_item_estimates = raw_item,
#'   conquest_case_eap = raw_case,
#'   conquest_population_term = "Term",
#'   conquest_population_estimate = "Est",
#'   conquest_item_id = "Item",
#'   conquest_item_estimate = "Est",
#'   conquest_case_person = "PID",
#'   conquest_case_estimate = "EAP"
#' )
#' review <- review_conquest_overlap(bundle, normalized)
#' summary(review)$summary
#' }
#' @name review_conquest_overlap
#' @export
review_conquest_overlap <- function(bundle,
                                    conquest_population = NULL,
                                    conquest_item_estimates = NULL,
                                    conquest_case_eap = NULL,
                                    conquest_population_term = "auto",
                                    conquest_population_estimate = "auto",
                                    conquest_item_id = "auto",
                                    conquest_item_estimate = "auto",
                                    item_id_source = c("auto", "response_var", "level"),
                                    conquest_case_person = "auto",
                                    conquest_case_estimate = "auto") {
  item_id_source <- match.arg(item_id_source)
  bundle <- validate_conquest_overlap_bundle_object(bundle)
  inputs <- coerce_conquest_overlap_tables_for_audit(
    conquest_population = conquest_population,
    conquest_item_estimates = conquest_item_estimates,
    conquest_case_eap = conquest_case_eap
  )
  conquest_population <- validate_conquest_overlap_input_df(inputs$conquest_population, "conquest_population")
  conquest_item_estimates <- validate_conquest_overlap_input_df(inputs$conquest_item_estimates, "conquest_item_estimates")
  conquest_case_eap <- validate_conquest_overlap_input_df(inputs$conquest_case_eap, "conquest_case_eap")

  if (isTRUE(inputs$standardized)) {
    conquest_population_term <- "Parameter"
    conquest_population_estimate <- "Estimate"
    conquest_item_id <- "ItemID"
    conquest_item_estimate <- "Estimate"
    conquest_case_person <- "Person"
    conquest_case_estimate <- "Estimate"
  }

  pop_term_col <- resolve_conquest_overlap_column(conquest_population, conquest_population_term, "conquest_population_term", "conquest_population")
  pop_est_col <- resolve_conquest_overlap_column(conquest_population, conquest_population_estimate, "conquest_population_estimate", "conquest_population")
  item_id_col <- resolve_conquest_overlap_column(conquest_item_estimates, conquest_item_id, "conquest_item_id", "conquest_item_estimates")
  item_est_col <- resolve_conquest_overlap_column(conquest_item_estimates, conquest_item_estimate, "conquest_item_estimate", "conquest_item_estimates")
  case_person_col <- resolve_conquest_overlap_column(conquest_case_eap, conquest_case_person, "conquest_case_person", "conquest_case_eap")
  case_est_col <- resolve_conquest_overlap_column(conquest_case_eap, conquest_case_estimate, "conquest_case_estimate", "conquest_case_eap")

  attention <- data.frame(
    Section = character(0),
    ID = character(0),
    Issue = character(0),
    Detail = character(0),
    stringsAsFactors = FALSE
  )

  pop_tbl <- as.data.frame(bundle$mfrmr_population, stringsAsFactors = FALSE)
  pop_raw_estimate <- conquest_population[[pop_est_col]]
  pop_numeric_estimate <- suppressWarnings(as.numeric(pop_raw_estimate))
  cq_pop <- data.frame(
    Parameter = as.character(conquest_population[[pop_term_col]]),
    ConQuestEstimate = pop_numeric_estimate,
    ConQuestEstimateNonNumeric = resolve_conquest_overlap_non_numeric_flag(
      conquest_population,
      pop_raw_estimate,
      pop_numeric_estimate
    ),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(cq_pop$Parameter)) {
    dup <- unique(cq_pop$Parameter[duplicated(cq_pop$Parameter)])
    attention <- rbind(
      attention,
      build_conquest_overlap_attention(
        "population",
        dup,
        "duplicate_conquest_parameter",
        "ConQuest population table should have one row per parameter."
      )
    )
    cq_pop <- cq_pop[!duplicated(cq_pop$Parameter), , drop = FALSE]
  }
  population_comparison <- merge(
    pop_tbl,
    cq_pop,
    by.x = "Parameter",
    by.y = "Parameter",
    all.x = TRUE,
    sort = FALSE
  )
  population_comparison$Difference <- population_comparison$ConQuestEstimate - population_comparison$Estimate
  non_numeric_pop_flag <- !is.na(population_comparison$ConQuestEstimateNonNumeric) &
    population_comparison$ConQuestEstimateNonNumeric
  population_comparison$Status <- ifelse(
    non_numeric_pop_flag,
    "NonNumericInConQuest",
    ifelse(
      is.na(population_comparison$ConQuestEstimate),
      "MissingInConQuest",
      ifelse(population_comparison$Parameter == "(Intercept)", "ConstraintDependent", "Compared")
    )
  )
  non_numeric_pop <- population_comparison$Parameter[population_comparison$Status == "NonNumericInConQuest"]
  if (length(non_numeric_pop) > 0L) {
    attention <- rbind(
      attention,
      build_conquest_overlap_attention(
        "population",
        non_numeric_pop,
        "non_numeric_conquest_parameter",
        "ConQuest population table has a non-numeric estimate for a parameter exported by the overlap bundle."
      )
    )
  }
  missing_pop <- population_comparison$Parameter[population_comparison$Status == "MissingInConQuest"]
  if (length(missing_pop) > 0L) {
    attention <- rbind(
      attention,
      build_conquest_overlap_attention(
        "population",
        missing_pop,
        "missing_conquest_parameter",
        "ConQuest population table is missing a parameter exported by the overlap bundle."
      )
    )
  }

  item_match <- choose_conquest_item_match(bundle, conquest_item_estimates[[item_id_col]], item_id_source = item_id_source)
  item_tbl <- as.data.frame(bundle$mfrmr_item_estimates, stringsAsFactors = FALSE)
  item_tbl$MatchID <- item_match$match_ids
  item_raw_estimate <- conquest_item_estimates[[item_est_col]]
  item_numeric_estimate <- suppressWarnings(as.numeric(item_raw_estimate))
  cq_item <- data.frame(
    MatchID = as.character(conquest_item_estimates[[item_id_col]]),
    ConQuestEstimate = item_numeric_estimate,
    ConQuestEstimateNonNumeric = resolve_conquest_overlap_non_numeric_flag(
      conquest_item_estimates,
      item_raw_estimate,
      item_numeric_estimate
    ),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(cq_item$MatchID)) {
    dup <- unique(cq_item$MatchID[duplicated(cq_item$MatchID)])
    attention <- rbind(
      attention,
      build_conquest_overlap_attention(
        "items",
        dup,
        "duplicate_conquest_item",
        "ConQuest item table should have one row per exported item identifier."
      )
    )
    cq_item <- cq_item[!duplicated(cq_item$MatchID), , drop = FALSE]
  }
  item_comparison <- merge(
    item_tbl,
    cq_item,
    by = "MatchID",
    all.x = TRUE,
    sort = FALSE
  )
  item_comparison$ConQuestCentered <- item_comparison$ConQuestEstimate - mean(item_comparison$ConQuestEstimate, na.rm = TRUE)
  item_comparison$CenteredDifference <- item_comparison$ConQuestCentered - item_comparison$CenteredEstimate
  non_numeric_item_flag <- !is.na(item_comparison$ConQuestEstimateNonNumeric) &
    item_comparison$ConQuestEstimateNonNumeric
  item_comparison$Status <- ifelse(
    non_numeric_item_flag,
    "NonNumericInConQuest",
    ifelse(is.na(item_comparison$ConQuestEstimate), "MissingInConQuest", "Compared")
  )
  non_numeric_items <- item_comparison$MatchID[item_comparison$Status == "NonNumericInConQuest"]
  if (length(non_numeric_items) > 0L) {
    attention <- rbind(
      attention,
      build_conquest_overlap_attention(
        "items",
        non_numeric_items,
        "non_numeric_conquest_item",
        "ConQuest item table has a non-numeric estimate for an item exported by the overlap bundle."
      )
    )
  }
  missing_items <- item_comparison$MatchID[item_comparison$Status == "MissingInConQuest"]
  if (length(missing_items) > 0L) {
    attention <- rbind(
      attention,
      build_conquest_overlap_attention(
        "items",
        missing_items,
        "missing_conquest_item",
        "ConQuest item table is missing an exported item identifier."
      )
    )
  }

  case_tbl <- as.data.frame(bundle$mfrmr_case_eap, stringsAsFactors = FALSE)
  case_raw_estimate <- conquest_case_eap[[case_est_col]]
  case_numeric_estimate <- suppressWarnings(as.numeric(case_raw_estimate))
  cq_case <- data.frame(
    Person = as.character(conquest_case_eap[[case_person_col]]),
    ConQuestEstimate = case_numeric_estimate,
    ConQuestEstimateNonNumeric = resolve_conquest_overlap_non_numeric_flag(
      conquest_case_eap,
      case_raw_estimate,
      case_numeric_estimate
    ),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(cq_case$Person)) {
    dup <- unique(cq_case$Person[duplicated(cq_case$Person)])
    attention <- rbind(
      attention,
      build_conquest_overlap_attention(
        "cases",
        dup,
        "duplicate_conquest_case",
        "ConQuest case table should have one row per person."
      )
    )
    cq_case <- cq_case[!duplicated(cq_case$Person), , drop = FALSE]
  }
  case_comparison <- merge(
    case_tbl[, c("Person", "Estimate", "SD", "Lower", "Upper", "Observations", "WeightedN", setdiff(names(case_tbl), c("Person", "Estimate", "SD", "Lower", "Upper", "Observations", "WeightedN"))), drop = FALSE],
    cq_case,
    by = "Person",
    all.x = TRUE,
    sort = FALSE
  )
  case_comparison$Difference <- case_comparison$ConQuestEstimate - case_comparison$Estimate
  non_numeric_case_flag <- !is.na(case_comparison$ConQuestEstimateNonNumeric) &
    case_comparison$ConQuestEstimateNonNumeric
  case_comparison$Status <- ifelse(
    non_numeric_case_flag,
    "NonNumericInConQuest",
    ifelse(is.na(case_comparison$ConQuestEstimate), "MissingInConQuest", "Compared")
  )
  non_numeric_cases <- case_comparison$Person[case_comparison$Status == "NonNumericInConQuest"]
  if (length(non_numeric_cases) > 0L) {
    attention <- rbind(
      attention,
      build_conquest_overlap_attention(
        "cases",
        non_numeric_cases,
        "non_numeric_conquest_case",
        "ConQuest case table has a non-numeric estimate for a person exported by the overlap bundle."
      )
    )
  }
  missing_cases <- case_comparison$Person[case_comparison$Status == "MissingInConQuest"]
  if (length(missing_cases) > 0L) {
    attention <- rbind(
      attention,
      build_conquest_overlap_attention(
        "cases",
        missing_cases,
        "missing_conquest_case",
        "ConQuest case table is missing a person exported by the overlap bundle."
      )
    )
  }

  direct_pop <- population_comparison[population_comparison$Status == "Compared", , drop = FALSE]
  direct_pop <- direct_pop[direct_pop$Parameter != "(Intercept)", , drop = FALSE]
  compared_items <- item_comparison[item_comparison$Status == "Compared", , drop = FALSE]
  compared_cases <- case_comparison[case_comparison$Status == "Compared", , drop = FALSE]

  overall <- data.frame(
    PopulationParametersExpected = nrow(population_comparison),
    PopulationParametersCompared = nrow(direct_pop),
    PopulationMae = if (nrow(direct_pop) > 0) mean(abs(direct_pop$Difference), na.rm = TRUE) else NA_real_,
    PopulationMaxAbsDifference = max_abs_conquest_overlap_difference(direct_pop$Difference),
    PopulationMaxAbsParameter = max_abs_conquest_overlap_label(direct_pop$Difference, direct_pop$Parameter),
    ItemRowsExpected = nrow(item_comparison),
    ItemRowsCompared = nrow(compared_items),
    ItemCenteredCorrelation = if (nrow(compared_items) > 1) suppressWarnings(stats::cor(compared_items$CenteredEstimate, compared_items$ConQuestCentered)) else NA_real_,
    ItemCenteredMae = if (nrow(compared_items) > 0) mean(abs(compared_items$CenteredDifference), na.rm = TRUE) else NA_real_,
    ItemCenteredMaxAbsDifference = max_abs_conquest_overlap_difference(compared_items$CenteredDifference),
    ItemCenteredMaxAbsItem = max_abs_conquest_overlap_label(compared_items$CenteredDifference, compared_items$MatchID),
    CaseRowsExpected = nrow(case_comparison),
    CaseRowsCompared = nrow(compared_cases),
    CaseCorrelation = if (nrow(compared_cases) > 1) suppressWarnings(stats::cor(compared_cases$Estimate, compared_cases$ConQuestEstimate)) else NA_real_,
    CaseMae = if (nrow(compared_cases) > 0) mean(abs(compared_cases$Difference), na.rm = TRUE) else NA_real_,
    CaseMaxAbsDifference = max_abs_conquest_overlap_difference(compared_cases$Difference),
    CaseMaxAbsPerson = max_abs_conquest_overlap_label(compared_cases$Difference, compared_cases$Person),
    AttentionItems = nrow(attention),
    AttentionMissing = count_conquest_overlap_attention(attention, "missing_"),
    AttentionDuplicate = count_conquest_overlap_attention(attention, "duplicate_"),
    AttentionNonNumeric = count_conquest_overlap_attention(attention, "non_numeric_"),
    stringsAsFactors = FALSE
  )

  settings <- dashboard_settings_table(list(
    item_id_source = item_match$source,
    conquest_population_term = pop_term_col,
    conquest_population_estimate = pop_est_col,
    conquest_item_id = item_id_col,
    conquest_item_estimate = item_est_col,
    conquest_case_person = case_person_col,
    conquest_case_estimate = case_est_col
  ))

  notes <- c(
    "This review compares normalized ConQuest tables against the exact-overlap mfrmr bundle.",
    "No raw ConQuest text parsing is assumed here; normalize external tables before review.",
    "Population slopes and sigma2 are intended for direct comparison, whereas item estimates are compared after centering.",
    "Non-numeric external estimate cells are treated as attention items rather than silently as ordinary missing rows."
  )

  out <- list(
    overall = overall,
    summary = overall,
    population_comparison = population_comparison,
    item_comparison = item_comparison,
    case_comparison = case_comparison,
    attention_items = attention,
    settings = settings,
    notes = notes
  )
  as_mfrm_bundle(out, "mfrm_conquest_overlap_review")
}

#' Export manuscript appendix tables from validated summary surfaces
#'
#' @param x A supported `summary()` source, a prebuilt
#'   [build_summary_table_bundle()] result, or a named list of such objects.
#' @param output_dir Directory where files will be written.
#' @param prefix File-name prefix for written artifacts.
#' @param include_html If `TRUE`, also write a lightweight HTML appendix page.
#' @param preset Appendix table-selection preset:
#'   `"all"` keeps every returned summary table,
#'   `"recommended"` keeps manuscript-facing summary tables while dropping
#'   bridge-only or preview-only surfaces, and
#'   `"compact"` keeps a smaller reviewer-facing subset.
#'   Section-aware presets `"methods"`, `"results"`, `"diagnostics"`, and
#'   `"reporting"` keep only the returned tables classified to those appendix
#'   sections in the summary-table catalog.
#' @param overwrite If `FALSE`, refuse to overwrite existing files.
#' @param zip_bundle If `TRUE`, attempt to zip the written appendix artifacts.
#' @param zip_name Optional zip-file name. Defaults to `"{prefix}_appendix.zip"`.
#' @param digits Digits forwarded when raw objects must be normalized through
#'   [build_summary_table_bundle()].
#' @param top_n Row cap forwarded when raw objects must be normalized through
#'   [build_summary_table_bundle()].
#' @param preview_chars Character cap forwarded when APA-output summaries must
#'   be normalized through [build_summary_table_bundle()].
#'
#' @details
#' This helper is the narrow public bridge from validated `summary()` surfaces
#' to manuscript appendix artifacts. It accepts the same reporting objects that
#' [build_summary_table_bundle()] supports, exports their table bundles as CSV,
#' and optionally assembles a lightweight HTML appendix page.
#'
#' Fit-level caveats are exported through the `analysis_caveats` role, and
#' pre-fit score-support caveats are exported through the
#' `score_category_caveats` role. Both roles are classified as diagnostics, so
#' they remain available under `"recommended"` and `"diagnostics"` presets when
#' the source summary contains caveat rows.
#'
#' Precision-review summaries keep `fit_separation_basis` in the exported
#' precision-review role so fit, ZSTD, separation/reliability/strata, and
#' package QC thresholds can be reported without turning them into release or
#' recovery success gates.
#' Fit-measure and FACETS fit-review summaries keep df/ZSTD sensitivity and
#' optional external FACETS matching tables in the same precision-review lane.
#'
#' Parameter-recovery studies can be exported by passing
#' [evaluate_mfrm_recovery()] or [assess_mfrm_recovery()] output directly. The
#' exported bundle keeps the ADEMP-style simulation basis, recovery metrics,
#' replication status, adequacy checklist, thresholds, and next actions in
#' separate appendix-ready tables.
#'
#' Recovery-validation summaries from the packaged validation protocol can be
#' exported by passing `summary(validation)`, including top-line release
#' decisions, condition notes, diagnostic notes, and domain decisions.
#'
#' Unlike [export_mfrm_bundle()], this helper does not require a fitted model.
#' It is intended for the stage where compact reporting summaries already exist
#' and the task is to hand off appendix-ready tables, catalogs, and reporting
#' maps.
#'
#' @section Typical workflow:
#' 1. Build `summary(...)` objects from fit, diagnostics, data description,
#'    reporting checklist, or APA outputs.
#' 2. Call `export_summary_appendix(...)` on one object or a named list.
#' 3. Hand off the written CSV/HTML appendix artifacts to manuscript or QA
#'    workflows.
#'
#' @return A named list of class `mfrm_summary_appendix_export` with:
#' - `summary`
#' - `written_files`
#' - `selection_summary`
#' - `selection_table_summary`
#' - `selection_section_table_summary`
#' - `selection_handoff_table_summary`
#' - `selection_handoff_preset_summary`
#' - `selection_handoff_summary`
#' - `selection_handoff_bundle_summary`
#' - `selection_handoff_role_summary`
#' - `selection_handoff_role_section_summary`
#' - `selection_role_summary`
#' - `selection_section_summary`
#' - `selection_catalog`
#' - `settings`
#' - `notes`
#'
#' @seealso [build_summary_table_bundle()], [export_mfrm_bundle()],
#'   [apa_table()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' appendix <- export_summary_appendix(
#'   list(fit = fit, diagnostics = diag),
#'   output_dir = tempdir(),
#'   prefix = "mfrmr_appendix_example",
#'   include_html = TRUE,
#'   overwrite = TRUE
#' )
#' appendix$summary
#' }
#' @export
export_summary_appendix <- function(x,
                                    output_dir = ".",
                                    prefix = "mfrmr_appendix",
                                    include_html = TRUE,
                                    preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting"),
                                    overwrite = FALSE,
                                    zip_bundle = FALSE,
                                    zip_name = NULL,
                                    digits = 3,
                                    top_n = 10,
                                    preview_chars = 160) {
  output_dir <- as.character(output_dir[1])
  prefix <- as.character(prefix[1] %||% "mfrmr_appendix")
  if (!nzchar(prefix)) prefix <- "mfrmr_appendix"
  include_html <- isTRUE(include_html)
  preset <- match.arg(preset)
  overwrite <- isTRUE(overwrite)
  zip_bundle <- isTRUE(zip_bundle)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir, call. = FALSE)
  }

  summary_table_bundles <- export_normalize_summary_table_inputs(
    x,
    digits = digits,
    top_n = top_n,
    preview_chars = preview_chars
  )
  if (length(summary_table_bundles) == 0L) {
    stop("`x` did not resolve to any summary-table bundle inputs.", call. = FALSE)
  }
  original_summary_table_bundles <- summary_table_bundles
  summary_table_bundles <- export_select_summary_table_bundles_for_appendix(
    summary_table_bundles,
    preset = preset
  )
  selection_catalog <- export_summary_table_selection_catalog(
    original_bundles = original_summary_table_bundles,
    selected_bundles = summary_table_bundles,
    preset = preset
  )
  selection_summary <- export_summary_table_selection_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )
  selection_table_summary <- export_summary_table_selection_table_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )
  selection_section_table_summary <- export_summary_table_selection_section_table_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )
  selection_handoff_table_summary <- export_summary_table_selection_handoff_table_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )
  selection_handoff_preset_summary <- export_summary_table_selection_handoff_preset_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )
  selection_handoff_summary <- export_summary_table_selection_handoff_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )
  selection_handoff_bundle_summary <- export_summary_table_selection_handoff_bundle_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )
  selection_handoff_role_summary <- export_summary_table_selection_handoff_role_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )
  selection_handoff_role_section_summary <- export_summary_table_selection_handoff_role_section_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )
  selection_role_summary <- export_summary_table_selection_role_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )
  selection_section_summary <- export_summary_table_selection_section_summary(
    selection_catalog = selection_catalog,
    preset = preset
  )

  emitted <- export_write_summary_table_bundles(
    summary_table_bundles = summary_table_bundles,
    output_dir = output_dir,
    prefix = prefix,
    overwrite = overwrite,
    html_tables = list(),
    html_text = list()
  )

  written_files <- emitted$written_files
  html_tables <- emitted$html_tables
  html_text <- emitted$html_text
  if (nrow(selection_summary) > 0L) {
    selection_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_summary.csv"))
    if (file.exists(selection_summary_path) && !overwrite) {
      stop("File already exists: ", selection_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_summary, file = selection_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_summary",
        Format = "csv",
        Path = normalizePath(selection_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_summary"]] <- selection_summary
  }
  if (nrow(selection_table_summary) > 0L) {
    selection_table_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_table_summary.csv"))
    if (file.exists(selection_table_summary_path) && !overwrite) {
      stop("File already exists: ", selection_table_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_table_summary, file = selection_table_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_table_summary",
        Format = "csv",
        Path = normalizePath(selection_table_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_table_summary"]] <- selection_table_summary
  }
  if (nrow(selection_section_table_summary) > 0L) {
    selection_section_table_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_section_table_summary.csv"))
    if (file.exists(selection_section_table_summary_path) && !overwrite) {
      stop("File already exists: ", selection_section_table_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_section_table_summary, file = selection_section_table_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_section_table_summary",
        Format = "csv",
        Path = normalizePath(selection_section_table_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_section_table_summary"]] <- selection_section_table_summary
  }
  if (nrow(selection_handoff_table_summary) > 0L) {
    selection_handoff_table_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_handoff_table_summary.csv"))
    if (file.exists(selection_handoff_table_summary_path) && !overwrite) {
      stop("File already exists: ", selection_handoff_table_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_handoff_table_summary, file = selection_handoff_table_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_handoff_table_summary",
        Format = "csv",
        Path = normalizePath(selection_handoff_table_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_handoff_table_summary"]] <- selection_handoff_table_summary
  }
  if (nrow(selection_handoff_preset_summary) > 0L) {
    selection_handoff_preset_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_handoff_preset_summary.csv"))
    if (file.exists(selection_handoff_preset_summary_path) && !overwrite) {
      stop("File already exists: ", selection_handoff_preset_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_handoff_preset_summary, file = selection_handoff_preset_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_handoff_preset_summary",
        Format = "csv",
        Path = normalizePath(selection_handoff_preset_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_handoff_preset_summary"]] <- selection_handoff_preset_summary
  }
  if (nrow(selection_handoff_summary) > 0L) {
    selection_handoff_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_handoff_summary.csv"))
    if (file.exists(selection_handoff_summary_path) && !overwrite) {
      stop("File already exists: ", selection_handoff_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_handoff_summary, file = selection_handoff_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_handoff_summary",
        Format = "csv",
        Path = normalizePath(selection_handoff_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_handoff_summary"]] <- selection_handoff_summary
  }
  if (nrow(selection_handoff_bundle_summary) > 0L) {
    selection_handoff_bundle_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_handoff_bundle_summary.csv"))
    if (file.exists(selection_handoff_bundle_summary_path) && !overwrite) {
      stop("File already exists: ", selection_handoff_bundle_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_handoff_bundle_summary, file = selection_handoff_bundle_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_handoff_bundle_summary",
        Format = "csv",
        Path = normalizePath(selection_handoff_bundle_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_handoff_bundle_summary"]] <- selection_handoff_bundle_summary
  }
  if (nrow(selection_handoff_role_summary) > 0L) {
    selection_handoff_role_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_handoff_role_summary.csv"))
    if (file.exists(selection_handoff_role_summary_path) && !overwrite) {
      stop("File already exists: ", selection_handoff_role_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_handoff_role_summary, file = selection_handoff_role_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_handoff_role_summary",
        Format = "csv",
        Path = normalizePath(selection_handoff_role_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_handoff_role_summary"]] <- selection_handoff_role_summary
  }
  if (nrow(selection_handoff_role_section_summary) > 0L) {
    selection_handoff_role_section_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_handoff_role_section_summary.csv"))
    if (file.exists(selection_handoff_role_section_summary_path) && !overwrite) {
      stop("File already exists: ", selection_handoff_role_section_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_handoff_role_section_summary, file = selection_handoff_role_section_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_handoff_role_section_summary",
        Format = "csv",
        Path = normalizePath(selection_handoff_role_section_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_handoff_role_section_summary"]] <- selection_handoff_role_section_summary
  }
  if (nrow(selection_catalog) > 0L) {
    selection_catalog_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_catalog.csv"))
    if (file.exists(selection_catalog_path) && !overwrite) {
      stop("File already exists: ", selection_catalog_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_catalog, file = selection_catalog_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_catalog",
        Format = "csv",
        Path = normalizePath(selection_catalog_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_catalog"]] <- selection_catalog
  }
  if (nrow(selection_role_summary) > 0L) {
    selection_role_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_role_summary.csv"))
    if (file.exists(selection_role_summary_path) && !overwrite) {
      stop("File already exists: ", selection_role_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_role_summary, file = selection_role_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_role_summary",
        Format = "csv",
        Path = normalizePath(selection_role_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_role_summary"]] <- selection_role_summary
  }
  if (nrow(selection_section_summary) > 0L) {
    selection_section_summary_path <- file.path(output_dir, paste0(prefix, "_appendix_selection_section_summary.csv"))
    if (file.exists(selection_section_summary_path) && !overwrite) {
      stop("File already exists: ", selection_section_summary_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(selection_section_summary, file = selection_section_summary_path, row.names = FALSE, na = "")
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_selection_section_summary",
        Format = "csv",
        Path = normalizePath(selection_section_summary_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
    html_tables[["appendix_selection_section_summary"]] <- selection_section_summary
  }

  if (isTRUE(include_html)) {
    html_path <- file.path(output_dir, paste0(prefix, "_appendix.html"))
    if (file.exists(html_path) && !overwrite) {
      stop("File already exists: ", html_path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    html_doc <- build_mfrm_bundle_html(
      title = paste0("mfrmr Manuscript Appendix: ", prefix),
      tables = html_tables,
      text_sections = html_text
    )
    writeLines(enc2utf8(html_doc), con = html_path, useBytes = TRUE)
    written_files <- rbind(
      written_files,
      data.frame(
        Component = "appendix_html",
        Format = "html",
        Path = normalizePath(html_path, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    )
  }

  zip_written <- FALSE
  zip_path <- NULL
  zip_note <- NULL
  if (isTRUE(zip_bundle) && nrow(written_files) > 0) {
    zip_file <- if (is.null(zip_name) || !nzchar(as.character(zip_name[1]))) {
      paste0(prefix, "_appendix.zip")
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
      written_files <- rbind(
        written_files,
        data.frame(
          Component = "appendix_zip",
          Format = "zip",
          Path = normalizePath(zip_path, winslash = "/", mustWork = FALSE),
          stringsAsFactors = FALSE
        )
      )
      zip_written <- TRUE
    } else if (inherits(zip_result, "error")) {
      zip_note <- conditionMessage(zip_result)
    }
  }

  summary_tbl <- data.frame(
    BundlesWritten = length(summary_table_bundles),
    FilesWritten = nrow(written_files),
    CsvWritten = sum(written_files$Format == "csv"),
    TextWritten = sum(written_files$Format == "txt"),
    HtmlWritten = sum(written_files$Format == "html"),
    ZipWritten = sum(written_files$Format == "zip"),
    stringsAsFactors = FALSE
  )

  settings <- dashboard_settings_table(list(
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    prefix = prefix,
    include_html = include_html,
    preset = preset,
    overwrite = overwrite,
    zip_bundle = zip_bundle,
    zip_written = zip_written,
    bundles_written = length(summary_table_bundles)
  ))

  notes <- paste0(
    "Exported summary-table appendix bundle(s) with preset `", preset, "`: ",
    paste(names(summary_table_bundles), collapse = ", "),
    "."
  )
  if (!is.null(zip_note) && nzchar(zip_note)) {
    notes <- c(notes, paste0("ZIP bundle was not created: ", zip_note))
  }

  out <- list(
    summary = summary_tbl,
    written_files = written_files,
    selection_summary = selection_summary,
    selection_table_summary = selection_table_summary,
    selection_section_table_summary = selection_section_table_summary,
    selection_handoff_table_summary = selection_handoff_table_summary,
    selection_handoff_preset_summary = selection_handoff_preset_summary,
    selection_handoff_summary = selection_handoff_summary,
    selection_handoff_bundle_summary = selection_handoff_bundle_summary,
    selection_handoff_role_summary = selection_handoff_role_summary,
    selection_handoff_role_section_summary = selection_handoff_role_section_summary,
    selection_role_summary = selection_role_summary,
    selection_section_summary = selection_section_summary,
    selection_catalog = selection_catalog,
    settings = settings,
    notes = notes
  )
  as_mfrm_bundle(out, "mfrm_summary_appendix_export")
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
#' @param summary_tables Optional manuscript-summary bundle input. Can be
#'   [build_summary_table_bundle()] output, any object supported by
#'   `build_summary_table_bundle()`, or a named list of such objects. When
#'   `NULL` and `"summary_tables"` is requested in `include`, a default set is
#'   built from `fit`, `diagnostics`, [reporting_checklist()], and
#'   [build_apa_outputs()]. Recovery-validation summaries can be supplied here
#'   to co-locate release-review appendix tables with a fit-based export bundle.
#' @param output_dir Directory where files will be written.
#' @param prefix File-name prefix.
#' @param include Components to export. Supported values are
#'   `"core_tables"`, `"checklist"`, `"dashboard"`, `"apa"`, `"anchors"`,
#'   `"manifest"`, `"visual_summaries"`, `"predictions"`, `"summary_tables"`,
#'   `"script"`, and `"html"`.
#' @param facet Optional facet for [facet_quality_dashboard()].
#' @param include_person_anchors If `TRUE`, include person measures in the
#'   exported anchor table.
#' @param overwrite If `FALSE`, refuse to overwrite existing files.
#' @param zip_bundle If `TRUE`, attempt to zip the written files into a single
#'   archive using [utils::zip()]. This is best-effort and may depend on the
#'   local R installation.
#' @param zip_name Optional zip-file name. Defaults to `"{prefix}_bundle.zip"`.
#' @param data Optional original analysis data frame. When supplied,
#'   `export_mfrm_bundle()` co-locates a CSV copy of the data
#'   alongside the replay script and updates the script's
#'   `read.csv()` path to point at it. The manifest's `input_hash`
#'   row for `data` is also computed against the user's untouched
#'   input so the recorded fingerprint matches what the replay
#'   script will load. Default `NULL` falls back to the legacy
#'   `your_data.csv` placeholder path.
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
#' - `"script"` for reproducibility and reruns. For latent-regression fits,
#'   this also writes the fit-level replay person-data sidecar when available.
#' - `"html"` for a light, shareable summary page. When replay sidecars are
#'   present, the HTML shows an artifact index for them rather than embedding
#'   the raw person-level replay table.
#' - `"summary_tables"` for manuscript-facing CSV exports of validated
#'   `summary()` surfaces and their compact indexes.
#' - `"visual_summaries"` when you want warning maps or residual PCA summaries
#'   to travel with the bundle.
#'
#' @section Recommended presets:
#' Common starting points are:
#' - minimal tables: `include = c("core_tables", "manifest")`
#' - reporting bundle: `include = c("core_tables", "checklist", "dashboard",
#'   "summary_tables", "html")`
#' - archival bundle: `include = c("core_tables", "manifest", "script",
#'   "visual_summaries", "html")`
#'
#' @section Written outputs:
#' Depending on `include`, the exporter can write:
#' - core CSV tables via [export_mfrm()]
#' - checklist CSVs via [reporting_checklist()]
#' - facet-dashboard CSVs via [facet_quality_dashboard()]
#' - APA text files via [build_apa_outputs()]
#' - manuscript-summary CSVs via [build_summary_table_bundle()]
#' - anchor CSV via [make_anchor_table()]
#' - manifest CSV/TXT via [build_mfrm_manifest()]
#' - visual warning/summary artifacts via [build_visual_summaries()]
#' - prediction/forecast CSVs via [predict_mfrm_population()],
#'   [predict_mfrm_units()], and [sample_mfrm_plausible_values()]
#' - a package-native replay script via [build_mfrm_replay_script()]
#' - for latent-regression fits, a replay-side person-data CSV paired with the
#'   replay script
#' - a lightweight HTML report that bundles the exported tables/text and, for
#'   replay sidecars, an artifact summary instead of raw person-level rows
#'
#' For latent-regression fits, prediction-side artifacts can carry the fitted
#' population-model scoring basis when you explicitly supply the corresponding
#' prediction objects. [predict_mfrm_population()] remains the scenario-level
#' forecast helper, whereas [predict_mfrm_units()] and
#' [sample_mfrm_plausible_values()] are the scoring layer.
#' To keep exports and replay scripts practical, large future-planning schemas
#' from scenario-level population predictions are not flattened into
#' `*_population_prediction_settings.csv` or ADeMP CSVs; the compact simulation
#' specification files carry the replay-relevant settings instead.
#'
#' For bounded `GPCM`, this exporter is available as a caveated partial bundle
#' over supported diagnostics, report text, visual summaries, manifests, and
#' replay scripts. The returned object and manifest include `gpcm_boundary`.
#' Package-native bounded-`GPCM` scorefile export is available with caveats,
#' while full FACETS-style score-side contract review and design forecasting
#' remain outside this bundle contract.
#'
#' @section Interpreting output:
#' The returned object reports both high-level bundle status and the exact files
#' written. In practice, `bundle$summary` is the direct status check, while
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
#'   [export_mfrm()], [reporting_checklist()], [export_summary_appendix()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
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
                               summary_tables = NULL,
                               output_dir = ".",
                               prefix = "mfrmr_bundle",
                               include = c("core_tables", "checklist", "dashboard", "apa", "anchors", "manifest", "visual_summaries", "predictions", "summary_tables", "script", "html"),
                               facet = NULL,
                               include_person_anchors = FALSE,
                               overwrite = FALSE,
                               zip_bundle = FALSE,
                               zip_name = NULL,
                               data = NULL) {
  include <- unique(tolower(as.character(include)))
  allowed <- c("core_tables", "checklist", "dashboard", "apa", "anchors", "manifest", "visual_summaries", "predictions", "summary_tables", "script", "html")
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
    residual_pca = if ("visual_summaries" %in% include) "both" else "none",
    gpcm_helper = "export_mfrm_bundle()"
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
  fit_population <- fit$population %||% list(
    active = FALSE,
    person_table_replay = NULL
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
  checklist_obj <- NULL
  apa_obj <- NULL
  summary_table_bundles <- list()

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

  round_export_df <- function(df, digits = 6) {
    if (!is.data.frame(df) || nrow(df) == 0L) {
      return(df)
    }
    out <- as.data.frame(df, stringsAsFactors = FALSE)
    numeric_cols <- vapply(out, is.numeric, logical(1))
    out[numeric_cols] <- lapply(out[numeric_cols], round, digits = digits)
    out
  }

  compact_population_prediction_settings <- function(settings) {
    settings <- settings %||% list()
    if (is.list(settings) && "planning_schema" %in% names(settings)) {
      settings$planning_schema <- "omitted_from_export_flattening; see compact population_prediction_sim_spec files"
    }
    settings
  }

  compact_population_prediction_ademp <- function(ademp) {
    ademp <- ademp %||% list()
    if (is.list(ademp) && is.list(ademp$data_generating_mechanism)) {
      if ("planning_schema" %in% names(ademp$data_generating_mechanism)) {
        ademp$data_generating_mechanism$planning_schema <-
          "omitted_from_export_flattening; see compact population_prediction_sim_spec files"
      }
    }
    ademp
  }

  get_checklist_obj <- function() {
    if (is.null(checklist_obj)) {
      checklist_obj <<- reporting_checklist(
        fit = fit,
        diagnostics = diagnostics,
        bias_results = if (inherits(bias_results, "mfrm_bias_collection")) bias_results else bias_inputs
      )
    }
    checklist_obj
  }

  get_apa_obj <- function() {
    if (is.null(apa_obj)) {
      apa_obj <<- build_apa_outputs(
        fit = fit,
        diagnostics = diagnostics,
        bias_results = export_primary_bias_result(bias_results, bias_inputs)
      )
    }
    apa_obj
  }

  write_sim_spec_bundle <- function(sim_spec, prefix_base) {
    if (!inherits(sim_spec, "mfrm_sim_spec")) return(invisible(NULL))

    scalar_settings <- list(
      source = sim_spec$source %||% "unknown",
      model = sim_spec$model %||% NA_character_,
      step_facet = sim_spec$step_facet %||% NA_character_,
      slope_facet = sim_spec$slope_facet %||% NA_character_,
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

    slope_tbl <- as.data.frame(sim_spec$slope_table %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(slope_tbl) > 0) {
      write_csv(
        slope_tbl,
        paste0(prefix, "_", prefix_base, "_slopes.csv"),
        paste0(prefix_base, "_slopes")
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
    checklist <- get_checklist_obj()
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
    apa <- get_apa_obj()
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

  if ("summary_tables" %in% include) {
    if (is.null(summary_tables)) {
      summary_table_bundles <- export_normalize_summary_table_inputs(
        list(
          fit = fit,
          diagnostics = diagnostics,
          checklist = get_checklist_obj(),
          apa = get_apa_obj()
        ),
        digits = 3,
        top_n = 10,
        preview_chars = 160
      )
    } else {
      summary_table_bundles <- export_normalize_summary_table_inputs(
        summary_tables,
        digits = 3,
        top_n = 10,
        preview_chars = 160
      )
    }
    emitted <- export_write_summary_table_bundles(
      summary_table_bundles = summary_table_bundles,
      output_dir = output_dir,
      prefix = prefix,
      overwrite = overwrite,
      html_tables = html_tables,
      html_text = html_text
    )
    written_files <- rbind(written_files, emitted$written_files)
    html_tables <- emitted$html_tables
    html_text <- emitted$html_text
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
      pop_forecast <- round_export_df(population_prediction$forecast, digits = 6)
      pop_overview <- round_export_df(population_prediction$overview, digits = 6)
      write_csv(population_prediction$design, paste0(prefix, "_population_prediction_design.csv"), "population_prediction_design")
      write_csv(population_prediction$forecast, paste0(prefix, "_population_prediction_forecast.csv"), "population_prediction_forecast")
      write_csv(population_prediction$overview, paste0(prefix, "_population_prediction_overview.csv"), "population_prediction_overview")
      write_settings_table(
        compact_population_prediction_settings(population_prediction$settings),
        paste0(prefix, "_population_prediction_settings.csv"),
        "population_prediction_settings"
      )
      write_sim_spec_bundle(population_prediction$sim_spec, "population_prediction_sim_spec")
      if (!is.null(population_prediction$ademp) && length(population_prediction$ademp) > 0) {
        write_csv(
          export_flatten_named_values(compact_population_prediction_ademp(population_prediction$ademp)),
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
      html_tables$population_prediction_forecast <- pop_forecast
      html_tables$population_prediction_overview <- pop_overview
    }

    if (!is.null(unit_prediction)) {
      unit_sum <- summary(unit_prediction, digits = 6)
      write_csv(unit_prediction$estimates, paste0(prefix, "_unit_prediction_estimates.csv"), "unit_prediction_estimates")
      write_csv(unit_prediction$row_review, paste0(prefix, "_unit_prediction_row_review.csv"), "unit_prediction_row_review")
      write_settings_table(unit_prediction$settings, paste0(prefix, "_unit_prediction_settings.csv"), "unit_prediction_settings")
      if (nrow(as.data.frame(unit_prediction$population_review %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
        write_csv(
          unit_prediction$population_review,
          paste0(prefix, "_unit_prediction_population_review.csv"),
          "unit_prediction_population_review"
        )
      }
      if (nrow(as.data.frame(unit_prediction$person_data %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
        write_csv(
          unit_prediction$person_data,
          paste0(prefix, "_unit_prediction_person_data.csv"),
          "unit_prediction_person_data"
        )
      }
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
      html_tables$unit_prediction_row_review <- unit_sum$row_review
      if (nrow(as.data.frame(unit_sum$population_review %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
        html_tables$unit_prediction_population_review <- unit_sum$population_review
      }
    }

    if (!is.null(plausible_values)) {
      pv_sum <- summary(plausible_values, digits = 6)
      write_csv(plausible_values$values, paste0(prefix, "_plausible_values.csv"), "plausible_values")
      write_csv(plausible_values$estimates, paste0(prefix, "_plausible_value_estimates.csv"), "plausible_value_estimates")
      write_csv(plausible_values$row_review, paste0(prefix, "_plausible_value_row_review.csv"), "plausible_value_row_review")
      write_settings_table(plausible_values$settings, paste0(prefix, "_plausible_value_settings.csv"), "plausible_value_settings")
      if (nrow(as.data.frame(plausible_values$population_review %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
        write_csv(
          plausible_values$population_review,
          paste0(prefix, "_plausible_value_population_review.csv"),
          "plausible_value_population_review"
        )
      }
      if (nrow(as.data.frame(plausible_values$person_data %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
        write_csv(
          plausible_values$person_data,
          paste0(prefix, "_plausible_value_person_data.csv"),
          "plausible_value_person_data"
        )
      }
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
      if (nrow(as.data.frame(pv_sum$population_review %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
        html_tables$plausible_value_population_review <- pv_sum$population_review
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
      include_person_anchors = include_person_anchors,
      data = data
    )
  }

  if ("manifest" %in% include && !is.null(manifest)) {
    write_csv(manifest$summary, paste0(prefix, "_manifest_summary.csv"), "manifest_summary")
    write_csv(manifest$environment, paste0(prefix, "_manifest_environment.csv"), "manifest_environment")
    write_csv(manifest$model_settings, paste0(prefix, "_manifest_model_settings.csv"), "manifest_model_settings")
    write_csv(manifest$source_columns, paste0(prefix, "_manifest_source_columns.csv"), "manifest_source_columns")
    write_csv(manifest$estimation_control, paste0(prefix, "_manifest_estimation_control.csv"), "manifest_estimation_control")
    write_csv(manifest$settings, paste0(prefix, "_manifest_settings.csv"), "manifest_settings")
    if (nrow(as.data.frame(manifest$anchor_summary, stringsAsFactors = FALSE)) > 0) {
      write_csv(manifest$anchor_summary, paste0(prefix, "_manifest_anchor_summary.csv"), "manifest_anchor_summary")
    }
    if (nrow(as.data.frame(manifest$anchors, stringsAsFactors = FALSE)) > 0) {
      write_csv(manifest$anchors, paste0(prefix, "_manifest_anchors.csv"), "manifest_anchors")
    }
    write_csv(manifest$available_outputs, paste0(prefix, "_manifest_available_outputs.csv"), "manifest_available_outputs")
    if (nrow(as.data.frame(manifest$gpcm_boundary %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
      write_csv(manifest$gpcm_boundary, paste0(prefix, "_manifest_gpcm_boundary.csv"), "manifest_gpcm_boundary")
      html_tables$manifest_gpcm_boundary <- manifest$gpcm_boundary
    }
    write_text(render_mfrm_manifest_text(manifest), paste0(prefix, "_manifest.txt"), "manifest_text")
    html_tables$manifest_summary <- manifest$summary
    html_tables$manifest_available_outputs <- manifest$available_outputs
    html_tables$manifest_settings <- manifest$settings
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
    replay_fit_person_data_file <- if (nrow(as.data.frame(fit_population$person_table_replay %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
      paste0(prefix, "_replay_fit_person_data.csv")
    } else {
      NULL
    }
    if (!is.null(replay_fit_person_data_file)) {
      write_csv(
        fit_population$person_table_replay,
        replay_fit_person_data_file,
        "replay_fit_person_data"
      )
    }
    # Bundle the user's original input data when supplied so the
    # replay script and its hash record point to a co-located file.
    # Falls back to the legacy `your_data.csv` placeholder when the
    # user did not pass `data` to `export_mfrm_bundle()`.
    replay_data_file <- "your_data.csv"
    if (is.data.frame(data) && nrow(data) > 0L) {
      replay_data_file <- paste0(prefix, "_replay_data.csv")
      write_csv(data, replay_data_file, "replay_data")
    }
    replay <- build_mfrm_replay_script(
      fit = if (ctx$input_mode == "facets_run") ctx$run else fit,
      diagnostics = diagnostics,
      bias_results = if (inherits(bias_results, "mfrm_bias_collection")) bias_results else bias_inputs,
      population_prediction = population_prediction,
      unit_prediction = unit_prediction,
      plausible_values = plausible_values,
      data_file = replay_data_file,
      fit_person_data_file = replay_fit_person_data_file,
      include_bundle = TRUE,
      bundle_dir = "replayed_bundle",
      bundle_prefix = prefix
    )
    write_script(replay$script, paste0(prefix, "_replay.R"), "replay_script")
    html_text$replay_script <- replay$script
    replay_artifacts <- data.frame(
      Artifact = "replay_script",
      Format = "r",
      File = paste0(prefix, "_replay.R"),
      Detail = "Executable package-native replay script.",
      stringsAsFactors = FALSE
    )
    if (!identical(replay_data_file, "your_data.csv")) {
      replay_artifacts <- rbind(
        replay_artifacts,
        data.frame(
          Artifact = "replay_data",
          Format = "csv",
          File = replay_data_file,
          Detail = paste0(
            "Original analysis data co-located with the replay script (",
            nrow(data), " row(s)."
          ),
          stringsAsFactors = FALSE
        )
      )
    }
    if (!is.null(replay_fit_person_data_file)) {
      replay_artifacts <- rbind(
        replay_artifacts,
        data.frame(
          Artifact = "replay_fit_person_data",
          Format = "csv",
          File = replay_fit_person_data_file,
          Detail = paste0(
            "Fit-level latent-regression replay person-data sidecar (",
            nrow(as.data.frame(fit_population$person_table_replay, stringsAsFactors = FALSE)),
            " row(s))."
          ),
          stringsAsFactors = FALSE
        )
      )
    }
    html_tables$replay_artifacts <- replay_artifacts
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
    gpcm_boundary = gpcm_capability_boundary_table(
      fit,
      helper = "export_mfrm_bundle()",
      extra_areas = c(
        "Score-side scorefile export under bounded GPCM",
        "FACETS output-contract score-side review",
        "Design planning and forecasting"
      )
    ),
    settings = settings,
    notes = notes
  )
  as_mfrm_bundle(out, "mfrm_export_bundle")
}

export_normalize_bias_inputs <- function(bias_results) {
  bias_results <- validate_bias_results_input(
    bias_results,
    helper = "export helpers"
  )
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

export_sanitize_component_tag <- function(x, fallback = "bundle") {
  x <- as.character(x[1] %||% "")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- tolower(x)
  if (!nzchar(x)) fallback else x
}

export_as_summary_table_bundle <- function(x,
                                           label = NULL,
                                           digits = 3,
                                           top_n = 10,
                                           preview_chars = 160) {
  bundle <- if (inherits(x, "mfrm_summary_table_bundle")) {
    x
  } else {
    build_summary_table_bundle(
      x,
      include_empty = TRUE,
      digits = digits,
      top_n = top_n,
      preview_chars = preview_chars
    )
  }

  if (is.null(label) || !nzchar(as.character(label[1] %||% ""))) {
    label <- bundle$source_class %||% class(x)[1] %||% "bundle"
  }
  list(
    label = export_sanitize_component_tag(label, fallback = "bundle"),
    bundle = bundle
  )
}

export_normalize_summary_table_inputs <- function(summary_tables,
                                                  digits = 3,
                                                  top_n = 10,
                                                  preview_chars = 160) {
  if (is.null(summary_tables)) return(list())

  summary_table_input_classes <- c(
    "mfrm_fit",
    "mfrm_diagnostics",
    "mfrm_precision_review",
    "mfrm_fit_measures",
    "mfrm_facets_fit_review",
    "mfrm_person_fit_indices",
    "mfrm_data_description",
    "mfrm_reporting_checklist",
    "mfrm_apa_outputs",
    "mfrm_design_evaluation",
    "mfrm_signal_detection",
    "mfrm_recovery_simulation",
    "mfrm_recovery_assessment",
    "mfrm_population_prediction",
    "mfrm_future_branch_active_branch",
    "mfrm_facets_run",
    "mfrm_bias",
    "mfrm_anchor_review",
    "mfrm_linking_review",
    "mfrm_misfit_casebook",
    "mfrm_weighting_review",
    "mfrm_unit_prediction",
    "mfrm_plausible_values",
    "summary.mfrm_fit",
    "summary.mfrm_diagnostics",
    "summary.mfrm_precision_review",
    "summary.mfrm_fit_measures",
    "summary.mfrm_facets_fit_review",
    "summary.mfrm_person_fit_indices",
    "summary.mfrm_data_description",
    "summary.mfrm_reporting_checklist",
    "summary.mfrm_apa_outputs",
    "summary.mfrm_design_evaluation",
    "summary.mfrm_signal_detection",
    "summary.mfrm_recovery_simulation",
    "summary.mfrm_recovery_assessment",
    "summary.mfrmr_recovery_validation",
    "summary.mfrm_population_prediction",
    "summary.mfrm_future_branch_active_branch",
    "summary.mfrm_facets_run",
    "summary.mfrm_bias",
    "summary.mfrm_anchor_review",
    "summary.mfrm_linking_review",
    "summary.mfrm_misfit_casebook",
    "summary.mfrm_weighting_review",
    "summary.mfrm_unit_prediction",
    "summary.mfrm_plausible_values"
  )

  if (inherits(summary_tables, "mfrm_summary_table_bundle")) {
    one <- export_as_summary_table_bundle(
      summary_tables,
      digits = digits,
      top_n = top_n,
      preview_chars = preview_chars
    )
    out <- list(one$bundle)
    names(out) <- one$label
    return(out)
  }

  if (inherits(summary_tables, summary_table_input_classes)) {
    one <- export_as_summary_table_bundle(
      summary_tables,
      digits = digits,
      top_n = top_n,
      preview_chars = preview_chars
    )
    out <- list(one$bundle)
    names(out) <- one$label
    return(out)
  }

  if (is.list(summary_tables) && !is.data.frame(summary_tables)) {
    if (!is.null(summary_tables$overview) &&
        !is.null(summary_tables$table_index) &&
        !is.null(summary_tables$tables)) {
      one <- export_as_summary_table_bundle(
        summary_tables,
        digits = digits,
        top_n = top_n,
        preview_chars = preview_chars
      )
      out <- list(one$bundle)
      names(out) <- one$label
      return(out)
    }

    nms <- names(summary_tables)
    if (is.null(nms)) nms <- rep("", length(summary_tables))
    out <- vector("list", length(summary_tables))
    out_names <- character(length(summary_tables))
    for (i in seq_along(summary_tables)) {
      one <- export_as_summary_table_bundle(
        summary_tables[[i]],
        label = nms[i],
        digits = digits,
        top_n = top_n,
        preview_chars = preview_chars
      )
      out[[i]] <- one$bundle
      out_names[i] <- one$label
    }
    if (length(out_names) > 0L) {
      dup <- duplicated(out_names)
      if (any(dup)) {
        counts <- stats::ave(seq_along(out_names), out_names, FUN = seq_along)
        out_names[dup] <- paste0(out_names[dup], "_", counts[dup])
      }
    }
    names(out) <- out_names
    return(out)
  }

  one <- export_as_summary_table_bundle(
    summary_tables,
    digits = digits,
    top_n = top_n,
    preview_chars = preview_chars
  )
  out <- list(one$bundle)
  names(out) <- one$label
  out
}

export_select_summary_table_bundles_for_appendix <- function(summary_table_bundles,
                                                             preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  if (identical(preset, "all")) {
    return(summary_table_bundles)
  }
  lapply(summary_table_bundles, function(bundle) {
    summary_table_bundle_select_for_appendix(bundle, preset = preset)
  })
}

export_summary_table_selection_catalog <- function(original_bundles,
                                                   selected_bundles,
                                                   preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  if (length(original_bundles) == 0L) {
    return(data.frame())
  }

  out <- do.call(
    rbind,
    lapply(names(original_bundles), function(label) {
      catalog <- summary_table_bundle_catalog(original_bundles[[label]])
      if (nrow(catalog) == 0L) {
        return(data.frame())
      }
      selected_names <- names(selected_bundles[[label]]$tables %||% list())
      catalog$Bundle <- as.character(label)
      catalog$Preset <- preset
      catalog$Selected <- as.character(catalog$Table) %in% selected_names
      catalog
    })
  )
  rownames(out) <- NULL
  out
}

export_exact_fraction <- function(num, den) {
  num <- suppressWarnings(as.numeric(num))
  den <- suppressWarnings(as.numeric(den))
  if (!is.finite(num) || !is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  num / den
}

export_summary_table_selection_summary <- function(selection_catalog,
                                                  preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(selection_catalog) == 0L) {
    return(data.frame())
  }

  out <- do.call(
    rbind,
    lapply(split(selection_catalog, selection_catalog$Bundle), function(part) {
      selected <- part[part$Selected %in% TRUE, , drop = FALSE]
      tables_selected <- nrow(selected)
      tables_available <- nrow(part)
      plot_ready_selected <- sum(selected$PlotReady %in% TRUE, na.rm = TRUE)
      numeric_selected <- sum(suppressWarnings(as.numeric(selected$NumericColumns)) > 0, na.rm = TRUE)
      data.frame(
        Bundle = as.character(part$Bundle[1]),
        Preset = preset,
        TablesAvailable = tables_available,
        TablesSelected = tables_selected,
        SelectionFraction = export_exact_fraction(tables_selected, tables_available),
        PlotReadySelected = plot_ready_selected,
        PlotReadyFraction = export_exact_fraction(plot_ready_selected, tables_selected),
        NumericSelected = numeric_selected,
        NumericFraction = export_exact_fraction(numeric_selected, tables_selected),
        RolesCovered = length(unique(as.character(selected$Role))),
        SectionsCovered = summary_table_bundle_compact_labels(unique(as.character(selected$AppendixSection)), max_n = 4L),
        KeyTables = summary_table_bundle_compact_labels(selected$Table, max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out
}

export_summary_table_selection_table_summary <- function(selection_catalog,
                                                         preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  selection_catalog <- selection_catalog[selection_catalog$Selected %in% TRUE, , drop = FALSE]
  if (nrow(selection_catalog) == 0L || !"Table" %in% names(selection_catalog)) {
    return(data.frame())
  }

  split_tbl <- split(selection_catalog, as.character(selection_catalog$Table))
  out <- do.call(
    rbind,
    lapply(names(split_tbl), function(table_nm) {
      part <- split_tbl[[table_nm]]
      data.frame(
        Preset = preset,
        Table = as.character(table_nm),
        PresetsSelected = length(unique(as.character(part$Preset))),
        Rows = suppressWarnings(as.integer(part$Rows[[1]] %||% 0L)),
        Role = as.character(part$Role[[1]] %||% ""),
        AppendixSection = as.character(part$AppendixSection[[1]] %||% ""),
        PreferredAppendixOrder = suppressWarnings(as.integer(part$PreferredAppendixOrder[[1]] %||% NA_integer_)),
        PlotReady = any(part$PlotReady %in% TRUE, na.rm = TRUE),
        ExportReady = any(part$ExportReady %in% TRUE, na.rm = TRUE),
        ApaTableReady = any(part$ApaTableReady %in% TRUE, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$PreferredAppendixOrder, out$Table, na.last = TRUE), , drop = FALSE]
}

export_summary_table_selection_handoff_table_summary <- function(selection_catalog,
                                                                 preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  if ("Preset" %in% names(selection_catalog)) {
    selection_catalog <- selection_catalog[
      as.character(selection_catalog$Preset %||% "") %in% preset,
      ,
      drop = FALSE
    ]
  }
  selection_catalog <- selection_catalog[selection_catalog$Selected %in% TRUE, , drop = FALSE]
  if (nrow(selection_catalog) == 0L || !"Table" %in% names(selection_catalog)) {
    return(data.frame())
  }

  out <- data.frame(
    Preset = as.character(selection_catalog$Preset %||% preset),
    AppendixSection = as.character(selection_catalog$AppendixSection %||% ""),
    Role = as.character(selection_catalog$Role %||% ""),
    Bundle = as.character(selection_catalog$Bundle %||% ""),
    Table = as.character(selection_catalog$Table %||% ""),
    Rows = suppressWarnings(as.integer(selection_catalog$Rows %||% NA_integer_)),
    NumericColumns = suppressWarnings(as.integer(selection_catalog$NumericColumns %||% NA_integer_)),
    PlotReady = selection_catalog$PlotReady %in% TRUE,
    ExportReady = selection_catalog$ExportReady %in% TRUE,
    ApaTableReady = selection_catalog$ApaTableReady %in% TRUE,
    RecommendedAppendix = selection_catalog$RecommendedAppendix %in% TRUE,
    CompactAppendix = selection_catalog$CompactAppendix %in% TRUE,
    PreferredAppendixOrder = suppressWarnings(as.integer(selection_catalog$PreferredAppendixOrder %||% NA_integer_)),
    AppendixRationale = as.character(selection_catalog$AppendixRationale %||% ""),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out <- unique(out)
  out[order(
    out$Preset,
    out$AppendixSection,
    out$Role,
    out$PreferredAppendixOrder,
    out$Bundle,
    out$Table,
    na.last = TRUE
  ), , drop = FALSE]
}

export_summary_table_selection_section_table_summary <- function(selection_catalog,
                                                                preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  if ("Preset" %in% names(selection_catalog)) {
    selection_catalog <- selection_catalog[
      as.character(selection_catalog$Preset %||% "") %in% preset,
      ,
      drop = FALSE
    ]
  }
  selection_catalog <- selection_catalog[selection_catalog$Selected %in% TRUE, , drop = FALSE]
  if (nrow(selection_catalog) == 0L ||
      !"AppendixSection" %in% names(selection_catalog) ||
      !"Table" %in% names(selection_catalog)) {
    return(data.frame())
  }

  split_tbl <- split(
    selection_catalog,
    paste(
      as.character(selection_catalog$AppendixSection %||% ""),
      as.character(selection_catalog$Table %||% ""),
      sep = "\r"
    )
  )
  out <- do.call(
    rbind,
    lapply(split_tbl, function(part) {
      rows_num <- suppressWarnings(as.numeric(part$Rows))
      rows_num <- rows_num[is.finite(rows_num)]
      num_cols <- suppressWarnings(as.numeric(part$NumericColumns))
      num_cols <- num_cols[is.finite(num_cols)]
      data.frame(
        Preset = preset,
        AppendixSection = as.character(part$AppendixSection[[1]] %||% ""),
        Table = as.character(part$Table[[1]] %||% ""),
        Rows = if (length(rows_num) == 0L) NA_integer_ else as.integer(max(rows_num, na.rm = TRUE)),
        NumericColumns = if (length(num_cols) == 0L) NA_integer_ else as.integer(max(num_cols, na.rm = TRUE)),
        PlotReady = any(part$PlotReady %in% TRUE, na.rm = TRUE),
        BundlesCovered = summary_table_bundle_compact_labels(unique(as.character(part$Bundle)), max_n = 4L),
        RolesCovered = summary_table_bundle_compact_labels(unique(as.character(part$Role)), max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Preset, out$AppendixSection, out$Table), , drop = FALSE]
}

export_summary_table_selection_handoff_summary <- function(selection_catalog,
                                                           preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  selection_catalog <- selection_catalog[selection_catalog$Selected %in% TRUE, , drop = FALSE]
  if (nrow(selection_catalog) == 0L || !"AppendixSection" %in% names(selection_catalog)) {
    return(data.frame())
  }

  split_tbl <- split(selection_catalog, as.character(selection_catalog$AppendixSection %||% ""))
  out <- do.call(
    rbind,
    lapply(names(split_tbl), function(section_nm) {
      part <- split_tbl[[section_nm]]
      tables <- nrow(part)
      plot_ready <- sum(part$PlotReady %in% TRUE, na.rm = TRUE)
      numeric_tables <- sum(suppressWarnings(as.numeric(part$NumericColumns)) > 0, na.rm = TRUE)
      data.frame(
        Preset = preset,
        AppendixSection = as.character(section_nm),
        Tables = tables,
        PlotReadyTables = plot_ready,
        PlotReadyFraction = export_exact_fraction(plot_ready, tables),
        NumericTables = numeric_tables,
        NumericFraction = export_exact_fraction(numeric_tables, tables),
        RolesCovered = summary_table_bundle_compact_labels(unique(as.character(part$Role)), max_n = 4L),
        KeyTables = summary_table_bundle_compact_labels(as.character(part$Table), max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Preset, out$AppendixSection), , drop = FALSE]
}

export_summary_table_selection_handoff_bundle_summary <- function(selection_catalog,
                                                                  preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  selection_catalog <- selection_catalog[selection_catalog$Selected %in% TRUE, , drop = FALSE]
  if (nrow(selection_catalog) == 0L ||
      !"AppendixSection" %in% names(selection_catalog) ||
      !"Bundle" %in% names(selection_catalog)) {
    return(data.frame())
  }

  split_tbl <- split(
    selection_catalog,
    paste(
      as.character(selection_catalog$AppendixSection %||% ""),
      as.character(selection_catalog$Bundle %||% ""),
      sep = "\r"
    )
  )
  out <- do.call(
    rbind,
    lapply(split_tbl, function(part) {
      tables <- nrow(part)
      plot_ready <- sum(part$PlotReady %in% TRUE, na.rm = TRUE)
      numeric_tables <- sum(suppressWarnings(as.numeric(part$NumericColumns)) > 0, na.rm = TRUE)
      data.frame(
        Preset = preset,
        AppendixSection = as.character(part$AppendixSection[[1]] %||% ""),
        Bundle = as.character(part$Bundle[[1]] %||% ""),
        Tables = tables,
        PlotReadyTables = plot_ready,
        PlotReadyFraction = export_exact_fraction(plot_ready, tables),
        NumericTables = numeric_tables,
        NumericFraction = export_exact_fraction(numeric_tables, tables),
        RolesCovered = summary_table_bundle_compact_labels(unique(as.character(part$Role)), max_n = 4L),
        KeyTables = summary_table_bundle_compact_labels(as.character(part$Table), max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Preset, out$AppendixSection, out$Bundle), , drop = FALSE]
}

export_summary_table_selection_handoff_preset_summary <- function(selection_catalog,
                                                                  preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  selection_catalog <- selection_catalog[selection_catalog$Selected %in% TRUE, , drop = FALSE]
  if (nrow(selection_catalog) == 0L) {
    return(data.frame())
  }

  tables <- nrow(selection_catalog)
  plot_ready <- sum(selection_catalog$PlotReady %in% TRUE, na.rm = TRUE)
  numeric_tables <- sum(suppressWarnings(as.numeric(selection_catalog$NumericColumns)) > 0, na.rm = TRUE)
  data.frame(
    Preset = preset,
    Tables = tables,
    PlotReadyTables = plot_ready,
    PlotReadyFraction = export_exact_fraction(plot_ready, tables),
    NumericTables = numeric_tables,
    NumericFraction = export_exact_fraction(numeric_tables, tables),
    BundlesCovered = length(unique(as.character(selection_catalog$Bundle %||% ""))),
    SectionsCovered = length(unique(as.character(selection_catalog$AppendixSection %||% ""))),
    RolesCovered = summary_table_bundle_compact_labels(unique(as.character(selection_catalog$Role)), max_n = 4L),
    KeySections = summary_table_bundle_compact_labels(unique(as.character(selection_catalog$AppendixSection)), max_n = 4L),
    KeyTables = summary_table_bundle_compact_labels(as.character(selection_catalog$Table), max_n = 4L),
    stringsAsFactors = FALSE
  )
}

export_summary_table_selection_handoff_role_summary <- function(selection_catalog,
                                                                preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  selection_catalog <- selection_catalog[selection_catalog$Selected %in% TRUE, , drop = FALSE]
  if (nrow(selection_catalog) == 0L || !"Role" %in% names(selection_catalog)) {
    return(data.frame())
  }

  split_tbl <- split(selection_catalog, as.character(selection_catalog$Role %||% ""))
  out <- do.call(
    rbind,
    lapply(names(split_tbl), function(role_nm) {
      part <- split_tbl[[role_nm]]
      tables <- nrow(part)
      plot_ready <- sum(part$PlotReady %in% TRUE, na.rm = TRUE)
      numeric_tables <- sum(suppressWarnings(as.numeric(part$NumericColumns)) > 0, na.rm = TRUE)
      data.frame(
        Preset = preset,
        Role = as.character(role_nm),
        Tables = tables,
        PlotReadyTables = plot_ready,
        PlotReadyFraction = export_exact_fraction(plot_ready, tables),
        NumericTables = numeric_tables,
        NumericFraction = export_exact_fraction(numeric_tables, tables),
        SectionsCovered = summary_table_bundle_compact_labels(unique(as.character(part$AppendixSection)), max_n = 4L),
        KeyTables = summary_table_bundle_compact_labels(as.character(part$Table), max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Preset, out$Tables, out$Role, decreasing = TRUE), , drop = FALSE]
}

export_summary_table_selection_handoff_role_section_summary <- function(selection_catalog,
                                                                        preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  selection_catalog <- selection_catalog[selection_catalog$Selected %in% TRUE, , drop = FALSE]
  if (nrow(selection_catalog) == 0L ||
      !"Role" %in% names(selection_catalog) ||
      !"AppendixSection" %in% names(selection_catalog)) {
    return(data.frame())
  }

  split_tbl <- split(
    selection_catalog,
    paste(as.character(selection_catalog$AppendixSection %||% ""),
          as.character(selection_catalog$Role %||% ""),
          sep = "\r")
  )
  out <- do.call(
    rbind,
    lapply(split_tbl, function(part) {
      tables <- nrow(part)
      plot_ready <- sum(part$PlotReady %in% TRUE, na.rm = TRUE)
      numeric_tables <- sum(suppressWarnings(as.numeric(part$NumericColumns)) > 0, na.rm = TRUE)
      data.frame(
        Preset = preset,
        AppendixSection = as.character(part$AppendixSection[[1]] %||% ""),
        Role = as.character(part$Role[[1]] %||% ""),
        Tables = tables,
        PlotReadyTables = plot_ready,
        PlotReadyFraction = export_exact_fraction(plot_ready, tables),
        NumericTables = numeric_tables,
        NumericFraction = export_exact_fraction(numeric_tables, tables),
        KeyTables = summary_table_bundle_compact_labels(as.character(part$Table), max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Preset, out$AppendixSection, out$Role), , drop = FALSE]
}

export_summary_table_selection_role_summary <- function(selection_catalog,
                                                        preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  selection_catalog <- selection_catalog[selection_catalog$Selected %in% TRUE, , drop = FALSE]
  if (nrow(selection_catalog) == 0L || !"Role" %in% names(selection_catalog)) {
    return(data.frame())
  }

  split_tbl <- split(selection_catalog, as.character(selection_catalog$Role %||% ""))
  out <- do.call(
    rbind,
    lapply(names(split_tbl), function(role_nm) {
      part <- split_tbl[[role_nm]]
      tables <- nrow(part)
      plot_ready <- sum(part$PlotReady %in% TRUE, na.rm = TRUE)
      numeric_tables <- sum(suppressWarnings(as.numeric(part$NumericColumns)) > 0, na.rm = TRUE)
      data.frame(
        Preset = preset,
        Role = as.character(role_nm),
        Tables = tables,
        PlotReadyTables = plot_ready,
        PlotReadyFraction = export_exact_fraction(plot_ready, tables),
        NumericTables = numeric_tables,
        NumericFraction = export_exact_fraction(numeric_tables, tables),
        BundlesCovered = length(unique(as.character(part$Bundle))),
        SectionsCovered = summary_table_bundle_compact_labels(unique(as.character(part$AppendixSection)), max_n = 4L),
        KeyTables = summary_table_bundle_compact_labels(as.character(part$Table), max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Tables, out$Role, decreasing = TRUE), , drop = FALSE]
}

export_summary_table_selection_section_summary <- function(selection_catalog,
                                                           preset = c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")) {
  preset <- match.arg(preset)
  selection_catalog <- as.data.frame(selection_catalog %||% data.frame(), stringsAsFactors = FALSE)
  selection_catalog <- selection_catalog[selection_catalog$Selected %in% TRUE, , drop = FALSE]
  if (nrow(selection_catalog) == 0L || !"AppendixSection" %in% names(selection_catalog)) {
    return(data.frame())
  }

  split_tbl <- split(selection_catalog, as.character(selection_catalog$AppendixSection %||% ""))
  out <- do.call(
    rbind,
    lapply(names(split_tbl), function(section_nm) {
      part <- split_tbl[[section_nm]]
      tables <- nrow(part)
      plot_ready <- sum(part$PlotReady %in% TRUE, na.rm = TRUE)
      numeric_tables <- sum(suppressWarnings(as.numeric(part$NumericColumns)) > 0, na.rm = TRUE)
      data.frame(
        Preset = preset,
        AppendixSection = as.character(section_nm),
        Tables = tables,
        PlotReadyTables = plot_ready,
        PlotReadyFraction = export_exact_fraction(plot_ready, tables),
        NumericTables = numeric_tables,
        NumericFraction = export_exact_fraction(numeric_tables, tables),
        BundlesCovered = length(unique(as.character(part$Bundle))),
        RolesCovered = summary_table_bundle_compact_labels(unique(as.character(part$Role)), max_n = 4L),
        KeyTables = summary_table_bundle_compact_labels(as.character(part$Table), max_n = 4L),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$Tables, out$AppendixSection, decreasing = TRUE), , drop = FALSE]
}

export_write_summary_table_bundles <- function(summary_table_bundles,
                                               output_dir,
                                               prefix,
                                               overwrite = FALSE,
                                               html_tables = list(),
                                               html_text = list()) {
  written_files <- data.frame(
    Component = character(0),
    Format = character(0),
    Path = character(0),
    stringsAsFactors = FALSE
  )

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

  write_csv_local <- function(df, filename, component) {
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("File already exists: ", path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    utils::write.csv(df, file = path, row.names = FALSE, na = "")
    add_written(component, "csv", path)
    invisible(path)
  }

  write_text_local <- function(text, filename, component) {
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("File already exists: ", path, ". Set `overwrite = TRUE` to replace.", call. = FALSE)
    }
    writeLines(enc2utf8(as.character(text)), con = path, useBytes = TRUE)
    add_written(component, "txt", path)
    invisible(path)
  }

  for (nm in names(summary_table_bundles)) {
    bundle_obj <- summary_table_bundles[[nm]]
    bundle_sum <- summary(bundle_obj)
    label <- export_sanitize_component_tag(nm, fallback = "summary")
    reserved_table_tags <- c("overview", "reporting_map")

    write_csv_local(bundle_sum$overview, paste0(prefix, "_summary_", label, "_overview.csv"), paste0("summary_", label, "_overview"))
    if (nrow(as.data.frame(bundle_sum$role_summary %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
      write_csv_local(bundle_sum$role_summary, paste0(prefix, "_summary_", label, "_role_summary.csv"), paste0("summary_", label, "_role_summary"))
    }
    if (nrow(as.data.frame(bundle_sum$table_catalog %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
      write_csv_local(bundle_sum$table_catalog, paste0(prefix, "_summary_", label, "_table_catalog.csv"), paste0("summary_", label, "_table_catalog"))
    }
    if (nrow(as.data.frame(bundle_sum$table_profile %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
      write_csv_local(bundle_sum$table_profile, paste0(prefix, "_summary_", label, "_table_profile.csv"), paste0("summary_", label, "_table_profile"))
    }
    if (nrow(as.data.frame(bundle_sum$plot_index %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
      write_csv_local(bundle_sum$plot_index, paste0(prefix, "_summary_", label, "_plot_index.csv"), paste0("summary_", label, "_plot_index"))
    }
    if (nrow(as.data.frame(bundle_sum$reporting_map %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
      write_csv_local(bundle_sum$reporting_map, paste0(prefix, "_summary_", label, "_reporting_map.csv"), paste0("summary_", label, "_reporting_map"))
    }
    write_csv_local(bundle_obj$table_index, paste0(prefix, "_summary_", label, "_table_index.csv"), paste0("summary_", label, "_table_index"))

    for (tbl_name in names(bundle_obj$tables %||% list())) {
      if (tbl_name %in% reserved_table_tags) next
      tbl <- as.data.frame(bundle_obj$tables[[tbl_name]], stringsAsFactors = FALSE)
      tag <- export_sanitize_component_tag(tbl_name, fallback = "table")
      write_csv_local(
        tbl,
        paste0(prefix, "_summary_", label, "_", tag, ".csv"),
        paste0("summary_", label, "_", tag)
      )
      html_tables[[paste0("summary_", label, "_", tbl_name)]] <- tbl
    }
    html_tables[[paste0("summary_", label, "_table_index")]] <- bundle_obj$table_index
    if (nrow(as.data.frame(bundle_sum$table_catalog %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
      html_tables[[paste0("summary_", label, "_table_catalog")]] <- bundle_sum$table_catalog
    }
    if (nrow(as.data.frame(bundle_sum$plot_index %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
      html_tables[[paste0("summary_", label, "_plot_index")]] <- bundle_sum$plot_index
    }
    if (nrow(as.data.frame(bundle_sum$reporting_map %||% data.frame(), stringsAsFactors = FALSE)) > 0) {
      html_tables[[paste0("summary_", label, "_reporting_map")]] <- bundle_sum$reporting_map
    }
    if (length(bundle_obj$notes %||% character(0)) > 0L) {
      note_text <- paste(as.character(bundle_obj$notes), collapse = "\n")
      write_text_local(
        note_text,
        paste0(prefix, "_summary_", label, "_notes.txt"),
        paste0("summary_", label, "_notes")
      )
      html_text[[paste0("summary_", label, "_notes")]] <- note_text
    }
  }

  list(
    written_files = written_files,
    html_tables = html_tables,
    html_text = html_text
  )
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

build_mfrm_bundle_html <- function(title,
                                   tables = list(),
                                   text_sections = list(),
                                   text_sections_after = list()) {
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

  text_sections_after <- text_sections_after[!vapply(text_sections_after, is.null, logical(1))]
  for (nm in names(text_sections_after)) {
    txt <- paste(as.character(text_sections_after[[nm]]), collapse = "\n")
    if (!nzchar(trimws(txt))) next
    parts <- c(
      parts,
      paste0("<h2>", html_escape(nm), "</h2>"),
      paste0("<pre>", html_escape(txt), "</pre>")
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
