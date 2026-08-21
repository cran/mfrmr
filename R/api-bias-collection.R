#' Estimate bias across multiple facet pairs
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. When `NULL`,
#'   diagnostics are computed with `residual_pca = "none"`.
#' @param pairs Optional list of facet specifications. Each element should be a
#'   character vector of length 2 or more, for example
#'   `list(c("Rater", "Criterion"), c("Task", "Criterion"))`.
#'   Explicit specifications may include canonical `"Person"` regardless of
#'   `include_person`. When `NULL`, all 2-way combinations of modeled facets
#'   are used.
#' @param include_person If `TRUE` and `pairs = NULL`, include `"Person"` in the
#'   automatically generated pair set. `"Person"` denotes the source column
#'   supplied as `person` to [fit_mfrm()] and remains opt-in because these
#'   cells condition on fitted JML person measures or MML EAP scores.
#' @param drop_empty If `TRUE`, omit empty bias tables from `by_pair` while still
#'   recording them in the summary table.
#' @param keep_errors If `TRUE`, retain per-pair error rows in the returned
#'   `errors` table instead of failing the whole batch.
#' @param max_abs Passed to [estimate_bias()].
#' @param omit_extreme Passed to [estimate_bias()].
#' @param max_iter Passed to [estimate_bias()].
#' @param tol Passed to [estimate_bias()].
#'
#' @details
#' This function orchestrates repeated calls to [estimate_bias()] across
#' multiple facet pairs and returns a consolidated bundle.
#'
#' **Bias/interaction** in MFRM refers to a systematic departure from
#' the additive model for a specific combination of facet elements
#' (e.g., a particular rater is unexpectedly harsh on a particular
#' criterion).  See [estimate_bias()] for the mathematical formulation.
#'
#' When `pairs = NULL`, the function builds all 2-way combinations of
#' modelled facets automatically.  For a model with facets Rater,
#' Criterion, and Task, this yields Rater\eqn{\times}Criterion,
#' Rater\eqn{\times}Task, and Criterion\eqn{\times}Task.
#' With `include_person = TRUE`, Person\eqn{\times}facet pairs are prepended in
#' the fit's facet order. They are conditional plug-in likelihood screens, not
#' jointly fitted Person\eqn{\times}facet terms. A single warning records this
#' boundary for the batch.
#'
#' The `summary` table aggregates results across pairs:
#' - `Rows`: number of interaction cells estimated
#' - `ScreenPositive`: count of cells with \eqn{|t| \ge 2}
#' - `Significant`: backward-compatible alias for `ScreenPositive`
#' - `MeanAbsBias`: average absolute bias magnitude (logits)
#' - `FormalInferenceEligible` / `PrimaryReportingEligible`: always `FALSE`;
#'   these conditional plug-in interaction statistics are screening evidence
#'
#' Per-pair failures (e.g., insufficient data for a sparse pair) are
#' captured in `errors` rather than stopping the entire batch.
#'
#' @section Output:
#' The returned object is a bundle-like list with class
#' `mfrm_bias_collection` and components such as:
#' - `summary`: one row per requested interaction
#' - `by_pair`: named list of successful [estimate_bias()] outputs
#' - `errors`: per-pair error log
#' - `settings`: resolved execution settings
#' - `primary`: first successful bias bundle, useful for downstream helpers
#'
#' @section Typical workflow:
#' 1. Fit with [fit_mfrm()] and diagnose with [diagnose_mfrm()]. For
#'    `RSM` / `PCM` reporting runs, prefer `method = "MML"` plus
#'    `diagnostic_mode = "both"` in the diagnostics call.
#' 2. Run `estimate_all_bias()` to compute the requested multi-pair interactions.
#' 3. Pass the resulting `by_pair` list into [reporting_checklist()] or
#'    [facet_quality_dashboard()].
#'
#' @return A named list with class `mfrm_bias_collection`.
#' @seealso [estimate_bias()], [reporting_checklist()], [facet_quality_dashboard()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", quad_points = 7, maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "both")
#' bias_all <- estimate_all_bias(fit, diagnostics = diag)
#' bias_all$summary[, c("Interaction", "Rows", "ScreenPositive")]
#' }
#' @export
estimate_all_bias <- function(fit,
                              diagnostics = NULL,
                              pairs = NULL,
                              include_person = FALSE,
                              drop_empty = TRUE,
                              keep_errors = TRUE,
                              max_abs = 10,
                              omit_extreme = TRUE,
                              max_iter = 4,
                              tol = 1e-3) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }

  diagnostics_supplied <- !is.null(diagnostics)
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (!is.list(diagnostics) || is.null(diagnostics$obs)) {
    stop("`diagnostics` must be output from diagnose_mfrm().", call. = FALSE)
  }

  pair_specs <- resolve_bias_collection_pairs(
    fit = fit,
    pairs = pairs,
    include_person = include_person
  )
  person_requested <- any(vapply(
    pair_specs,
    function(spec) "Person" %in% as.character(spec$facets),
    logical(1)
  ))
  if (person_requested) {
    warn_bias_person_screen(fit, helper = "estimate_all_bias()")
  }

  by_pair <- list()
  summary_rows <- vector("list", length(pair_specs))
  error_rows <- list()

  # Show an interactive progress bar when more than one pair is being
  # estimated. cli::cli_progress_bar is silent when not needed
  # (e.g. unattended R CMD check) and respects options(cli.progress_show_after).
  progress_id <- NULL
  if (length(pair_specs) > 1L) {
    progress_id <- cli::cli_progress_bar(
      name = "estimate_all_bias",
      total = length(pair_specs),
      format = "{cli::pb_spin} bias interactions: {cli::pb_current}/{cli::pb_total} [{cli::pb_elapsed}]",
      clear = TRUE,
      .envir = parent.frame()
    )
    on.exit(cli::cli_progress_done(id = progress_id), add = TRUE)
  }

  for (i in seq_along(pair_specs)) {
    spec <- pair_specs[[i]]
    label <- spec$label
    facets <- spec$facets

    if (!is.null(progress_id)) {
      cli::cli_progress_update(id = progress_id, set = i - 1L)
    }

    bias_obj <- tryCatch(
      withCallingHandlers(
        estimate_bias(
          fit = fit,
          diagnostics = diagnostics,
          interaction_facets = facets,
          max_abs = max_abs,
          omit_extreme = omit_extreme,
          max_iter = max_iter,
          tol = tol
        ),
        mfrmr_person_bias_screen_warning = function(cnd) {
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )

    if (inherits(bias_obj, "error")) {
      msg <- conditionMessage(bias_obj)
      summary_rows[[i]] <- data.frame(
        Interaction = label,
        Order = length(facets),
        Facets = paste(facets, collapse = " x "),
        Rows = 0L,
        ScreenPositive = 0L,
        Significant = 0L,
        MaxAbsT = NA_real_,
        MeanAbsBias = NA_real_,
        SupportsFormalInference = FALSE,
        FormalInferenceEligible = FALSE,
        PrimaryReportingEligible = FALSE,
        ClassificationSystem = "screening",
        ReportingUse = "screening_only",
        Kept = FALSE,
        Error = msg,
        stringsAsFactors = FALSE
      )
      if (isTRUE(keep_errors)) {
        error_rows[[length(error_rows) + 1L]] <- data.frame(
          Interaction = label,
          Facets = paste(facets, collapse = " x "),
          Error = msg,
          stringsAsFactors = FALSE
        )
      }
      next
    }

    tbl <- as.data.frame(bias_obj$table %||% data.frame(), stringsAsFactors = FALSE)
    t_vals <- if ("t" %in% names(tbl)) suppressWarnings(as.numeric(tbl$t)) else rep(NA_real_, nrow(tbl))
    bias_vals <- if ("Bias Size" %in% names(tbl)) {
      suppressWarnings(as.numeric(tbl[["Bias Size"]]))
    } else if ("BiasSize" %in% names(tbl)) {
      suppressWarnings(as.numeric(tbl[["BiasSize"]]))
    } else {
      rep(NA_real_, nrow(tbl))
    }

    kept <- nrow(tbl) > 0 || !isTRUE(drop_empty)
    if (kept) {
      by_pair[[label]] <- bias_obj
    }

    screen_positive <- sum(is.finite(t_vals) & abs(t_vals) >= 2, na.rm = TRUE)
    summary_rows[[i]] <- data.frame(
      Interaction = label,
      Order = length(facets),
      Facets = paste(facets, collapse = " x "),
      Rows = nrow(tbl),
      ScreenPositive = screen_positive,
      Significant = screen_positive,
      MaxAbsT = if (any(is.finite(t_vals))) max(abs(t_vals), na.rm = TRUE) else NA_real_,
      MeanAbsBias = if (any(is.finite(bias_vals))) mean(abs(bias_vals), na.rm = TRUE) else NA_real_,
      SupportsFormalInference = FALSE,
      FormalInferenceEligible = FALSE,
      PrimaryReportingEligible = FALSE,
      ClassificationSystem = "screening",
      ReportingUse = "screening_only",
      Kept = kept,
      Error = "",
      stringsAsFactors = FALSE
    )
  }

  summary_tbl <- if (length(summary_rows) == 0) {
    data.frame(
      Interaction = character(0),
      Order = integer(0),
      Facets = character(0),
      Rows = integer(0),
      ScreenPositive = integer(0),
      Significant = integer(0),
      MaxAbsT = numeric(0),
      MeanAbsBias = numeric(0),
      SupportsFormalInference = logical(0),
      FormalInferenceEligible = logical(0),
      PrimaryReportingEligible = logical(0),
      ClassificationSystem = character(0),
      ReportingUse = character(0),
      Kept = logical(0),
      Error = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    dplyr::bind_rows(summary_rows)
  }

  errors_tbl <- if (length(error_rows) == 0) {
    data.frame(
      Interaction = character(0),
      Facets = character(0),
      Error = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    dplyr::bind_rows(error_rows)
  }

  settings <- dashboard_settings_table(list(
    requested_pairs = length(pair_specs),
    successful_pairs = length(by_pair),
    diagnostics_supplied = diagnostics_supplied,
    include_person = isTRUE(include_person),
    person_pairs_requested = person_requested,
    drop_empty = isTRUE(drop_empty),
    keep_errors = isTRUE(keep_errors),
    max_abs = as.numeric(max_abs[1]),
    omit_extreme = isTRUE(omit_extreme),
    max_iter = as.integer(max_iter[1]),
    tol = as.numeric(tol[1])
  ))

  out <- list(
    summary = as.data.frame(summary_tbl, stringsAsFactors = FALSE),
    by_pair = by_pair,
    errors = as.data.frame(errors_tbl, stringsAsFactors = FALSE),
    settings = settings,
    primary = if (length(by_pair) > 0) by_pair[[1]] else NULL
  )
  as_mfrm_bundle(out, "mfrm_bias_collection")
}

resolve_bias_collection_pairs <- function(fit, pairs = NULL, include_person = FALSE) {
  modeled_facets <- as.character(fit$config$facet_names %||% character(0))
  modeled_facets <- modeled_facets[!is.na(modeled_facets) & nzchar(modeled_facets)]
  modeled_facets <- unique(modeled_facets)

  # `include_person` controls automatic enumeration only. Explicit `pairs`
  # may name the canonical Person role directly.
  known_facets <- unique(c("Person", modeled_facets))
  available <- modeled_facets
  if (isTRUE(include_person)) {
    available <- c("Person", available)
  }
  available <- unique(available)

  if (is.null(pairs)) {
    if (length(available) < 2L) {
      stop("At least two available facets are required for multi-pair bias estimation.", call. = FALSE)
    }
    auto_pairs <- utils::combn(available, 2, simplify = FALSE)
    return(lapply(auto_pairs, function(x) list(facets = x, label = paste(x, collapse = " x "))))
  }

  if (!is.list(pairs) || length(pairs) == 0) {
    stop("`pairs` must be NULL or a non-empty list of character vectors.", call. = FALSE)
  }

  out <- vector("list", length(pairs))
  for (i in seq_along(pairs)) {
    pair <- unique(as.character(pairs[[i]]))
    pair <- pair[!is.na(pair) & nzchar(pair)]
    if (length(pair) < 2L) {
      stop("Each element of `pairs` must contain at least two facet names.", call. = FALSE)
    }
    bad <- setdiff(pair, known_facets)
    if (length(bad) > 0) {
      stop(
        "Unknown facet(s) in `pairs[[", i, "]]`: ",
        paste(bad, collapse = ", "),
        ". Available: ",
        paste(known_facets, collapse = ", "),
        call. = FALSE
      )
    }
    out[[i]] <- list(
      facets = pair,
      label = paste(pair, collapse = " x ")
    )
  }
  out
}
