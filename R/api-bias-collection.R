#' Estimate bias across multiple facet pairs
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. When `NULL`,
#'   diagnostics are computed with `residual_pca = "none"`.
#'   Must match this fitted analysis; see [mfrmr_compatibility_layer] for
#'   migration of saved diagnostics.
#' @param pairs Optional list of facet specifications. Each element should be a
#'   character vector of length 2 or more, for example
#'   `list(c("Rater", "Criterion"), c("Task", "Criterion"))`.
#'   Explicit specifications may include canonical `"Person"` regardless of
#'   `include_person`. When `NULL`, all 2-way combinations of modeled facets
#'   are used. Each combination may appear only once, regardless of facet order.
#' @param include_person If `TRUE` and `pairs = NULL`, include `"Person"` in the
#'   automatically generated pair set. `"Person"` denotes the source column
#'   supplied as `person` to [fit_mfrm()] and remains opt-in because these
#'   cells condition on fitted JML person measures or MML EAP scores.
#' @param drop_empty If `TRUE`, omit empty bias tables from `by_pair` while still
#'   recording them in the summary table.
#' @param keep_errors If `TRUE`, also retain per-pair error rows in the returned
#'   `errors` table. Either way, errors remain in `summary$Error` and do not
#'   stop the batch.
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
#' - `Points`: cells with an available bias estimate
#'   (including bounded optima; inspect `BiasStatus` in the combined cells)
#' - `ScreenAvailable`, `ScreenUnavailable`: cells with/without a usable
#'   conditional screening statistic
#' - `ScreenPositive`: count of available cells with \eqn{|t| \ge 2};
#'   `NA` when no screening statistic is available, including failed pairs
#' - `Significant`: backward-compatible alias for `ScreenPositive`
#' - `MeanAbsBias`: average absolute available bias magnitude (native ability units)
#' - `Status`: `available`, `partial`, `unavailable`, `empty`, or `failed`;
#'   `available` refers to conditional screening, not formal inference
#' - `FormalInferenceEligible` / `PrimaryReportingEligible`: always `FALSE`;
#'   these conditional plug-in interaction statistics are screening evidence
#'
#' Per-pair failures (e.g., insufficient data for a sparse pair) are
#' captured in `errors` rather than stopping the entire batch.
#' There is no collection-wide multiplicity correction. Cells and pairs reuse
#' the same observations and fitted calibration. Adjusting their screening tail
#' areas would not establish formal family-wise error or false-discovery control.
#' A screen-positive result is a prompt for review, not proof of unfairness;
#' no screen-positive cells does not establish absence of interaction.
#' Search-boundary points remain visible, with missing screening statistics
#' and the `ScreenStatus`/`ScreenReason` returned by [estimate_bias()].
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
#' 3. Inspect `summary(bias_all)$overview` for availability and magnitude, then
#'    `summary(bias_all)$cells` for all retained cells and their observed counts.
#'    Use [plot.mfrm_bias_collection()] for the overview or a named pair.
#' 4. Use `report <- mfrm_report(bias_all)` to retain the full search, failed
#'    pairs and all cells in one report; [export_summary_appendix()] saves its
#'    tables. Supply `study_context` to state rubric endpoints and score direction.
#'    The `by_pair` list also feeds [reporting_checklist()] or [facet_quality_dashboard()].
#' 5. For a substantively justified joint interaction model, use
#'    `facet_interactions` in [fit_mfrm()] and [compare_mfrm()] on compatible
#'    fits. Selecting terms from these screens makes same-data comparisons
#'    exploratory; it does not create a prespecified confirmatory test.
#'
#' @return A named list with class `mfrm_bias_collection`.
#' @seealso [estimate_bias()], [mfrm_report()], [compare_mfrm()],
#'   [reporting_checklist()], [facet_quality_dashboard()]
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("estimate_all_bias", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML")
#' diag <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "both")
#' bias_all <- estimate_all_bias(fit, diagnostics = diag)
#' summary(bias_all)
#' head(summary(bias_all)$cells)
#' report <- mfrm_report(bias_all, study_context = c(
#'   score_direction = "Higher coded scores indicate better performance"))
#' report$tables$overview
#' head(report$tables$bias_cells)
#' }
#' @section Saved-result reuse:
#' The returned `fit_identity` records the fitted analysis for reuse with
#' reports and exports. Recreate older results from the saved fit and matching
#' diagnostics using the original screening options; no model refit is needed.
#' See [mfrmr_compatibility_layer] for migration instructions.
#' @export
#' @order 1
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
  stop_if_product_slopes(fit, "estimate_all_bias()")
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
  mfrm_results_validate_diagnostics_identity(fit, diagnostics, helper = "estimate_all_bias()")

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
    kept <- nrow(tbl) > 0 || !isTRUE(drop_empty)
    if (kept) {
      by_pair[[label]] <- bias_obj
    }

    summary_rows[[i]] <- data.frame(
      Interaction = label,
      Order = length(facets),
      Facets = paste(facets, collapse = " x "),
      Rows = nrow(tbl),
      ScreenPositive = NA_integer_,
      Significant = NA_integer_,
      MaxAbsT = NA_real_,
      MeanAbsBias = NA_real_,
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

  settings <- bundle_settings_table(list(
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
  out$calculation_version <- 2L
  out$fit_identity <- mfrm_result_fit_identity(fit)
  out <- as_mfrm_bundle(out, "mfrm_bias_collection")
  out$summary <- summary(out)$overview
  out
}

#' @rdname estimate_all_bias
#' @details `summary()` method: Summarize all retained bias screens
#' @param object A result from [estimate_all_bias()].
#' @param digits Number of significant digits used for printing; calculations
#'   and returned tables retain their precision.
#' @return For `summary()`: `summary()` returns a `summary.mfrm_bias_collection` list:
#'   `overview` has one row per requested interaction, `cells` combines retained
#'   cell tables with an `Interaction` identifier and explicit `PointAvailable`,
#'   `ScreenAvailable` and `ScreenPositive` fields, `errors` lists failed pairs,
#'   and `optimization_failures` retains cell-level optimizer records.
#'   Original cell columns, including facet levels, observed counts, signed
#'   `Bias Size`, conditional `S.E.`, `t` and any GPCM profile status, are
#'   preserved. `settings`, `fit_identity` and `notes` retain provenance and
#'   interpretation. Print methods return their input invisibly.
#' @details For `summary()`: No model, diagnostics or cell calculation is repeated. The overview
#'   counts finite bias estimates from successful cell optimizations separately
#'   from cells with a finite screening statistic and positive finite conditional
#'   SE. `ScreenPositive` is missing when no statistic can be screened, rather
#'   than reporting zero. Failed pairs remain visible even when the original
#'   batch used `keep_errors = FALSE`. Older collections without retained tables
#'   are labeled `not_retained`; they are not recalculated.
#'
#'   A positive bias shift increases expected coded scores conditional on the
#'   fitted model. It does not necessarily mean better performance or unfair
#'   leniency. Baseline calibration and Person uncertainty are held fixed.
#'   Different pairs use overlapping data, so do not add their observed counts
#'   to obtain a study sample size. The overview's mean absolute bias describes
#'   magnitude only; use the signed cell estimates to study direction.
#'   Counts cover retained [estimate_bias()] cells: unobserved combinations and
#'   cells omitted by the chosen screening options are not reconstructed.
#'   No multiplicity-adjusted confirmatory decisions are supplied.
#' @seealso [estimate_all_bias()], [plot.mfrm_bias_collection()], [mfrm_report()], [estimate_bias()]
#' @examples # Long-running workflow; run explicitly when needed.
#' # example("summary.mfrm_bias_collection", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML")
#' screens <- estimate_all_bias(fit)
#' review <- summary(screens)
#' review$overview
#' head(review$cells[, c("Interaction", "Bias Size", "ObsN", "ScreenAvailable")])
#' }
#' @export
summary.mfrm_bias_collection <- function(object, ...) {
  rlang::check_dots_empty()
  if (!is.null(object$fit_identity)) for (label in names(object$by_pair)) {
    if (length(object$by_pair[[label]])) mfrm_validate_result_identity(
      object$by_pair[[label]], object$fit_identity, paste0("by_pair[[", label, "]]"),
      "estimate_bias()", "summary.mfrm_bias_collection()")
  }
  overview <- as.data.frame(object$summary)
  cells <- dplyr::bind_rows(lapply(object$by_pair, function(x) x$table),
    .id = "Interaction")
  n <- nrow(cells)
  if (!"Interaction" %in% names(cells)) cells$Interaction <- character(n)
  numeric_column <- function(name) if (name %in% names(cells)) cells[[name]] else rep(NA_real_, n)
  bias <- numeric_column("Bias Size")
  se <- numeric_column("S.E.")
  statistic <- numeric_column("t")
  ok <- if ("OptimizationStatus" %in% names(cells)) cells$OptimizationStatus %in% "ok" else rep(TRUE, n)
  cells$PointAvailable <- is.finite(bias) & ok
  cells$ScreenAvailable <- cells$PointAvailable & is.finite(se) & se > 0 & is.finite(statistic)
  cells$ScreenPositive <- ifelse(cells$ScreenAvailable, abs(statistic) >= 2, NA)
  overview$Points <- overview$ScreenAvailable <- overview$ScreenUnavailable <- NA_integer_
  overview$Status <- rep("not_retained", nrow(overview))
  for (i in seq_len(nrow(overview))) {
    if (!is.na(overview$Error[i]) && nzchar(overview$Error[i])) {
      overview$Status[i] <- "failed"
      overview$ScreenPositive[i] <- overview$Significant[i] <- NA_integer_
      next
    }
    label <- overview$Interaction[i]
    retained <- label %in% names(object$by_pair)
    if (!retained && overview$Rows[i] > 0L) next
    at <- which(cells$Interaction == label)
    points <- cells$PointAvailable[at]
    screens <- cells$ScreenAvailable[at]
    overview$Points[i] <- sum(points)
    overview$ScreenAvailable[i] <- sum(screens)
    overview$ScreenUnavailable[i] <- length(at) - sum(screens)
    overview$ScreenPositive[i] <- overview$Significant[i] <-
      if (any(screens)) sum(cells$ScreenPositive[at], na.rm = TRUE) else NA_integer_
    overview$MeanAbsBias[i] <- if (any(points)) mean(abs(bias[at][points])) else NA_real_
    overview$MaxAbsT[i] <- if (any(screens)) max(abs(statistic[at][screens])) else NA_real_
    overview$Status[i] <- if (!length(at)) "empty" else if (all(screens)) "available" else
      if (any(points)) "partial" else "unavailable"
  }
  errors <- overview[!is.na(overview$Error) & nzchar(overview$Error),
    c("Interaction", "Facets", "Error"), drop = FALSE]
  failures <- dplyr::bind_rows(lapply(object$by_pair, function(x) x$optimization_failures),
    .id = "Interaction")
  structure(list(overview = overview, cells = as.data.frame(cells),
    errors = errors, optimization_failures = as.data.frame(failures),
    settings = object$settings, fit_identity = object$fit_identity,
    notes = c(
      if (n > 0L && !all(c("BiasStatus", "ScreenStatus", "FitStatus") %in% names(cells)))
        "Saved cells predate the boundary and complete-cell information/fit checks. Recompute estimate_all_bias() from the saved fit and matching diagnostics; these values have not been recalculated.",
      "Conditional plug-in screens: calibration and Person uncertainty are held fixed; no formal inference or collection-wide multiplicity control.",
      "Positive bias increases expected coded scores. Review signed magnitude, observed counts and unavailable cells before interpreting a flag.",
      "Pairs reuse observations. Do not sum cell counts across pairs as a study sample size; no flags does not establish no interaction.")),
    class = "summary.mfrm_bias_collection")
}

#' @rdname estimate_all_bias
#' @export
print.mfrm_bias_collection <- function(x, digits = 3, ...) {
  print(summary(x), digits = digits, ...)
  invisible(x)
}

#' @rdname estimate_all_bias
#' @export
print.summary.mfrm_bias_collection <- function(x, digits = 3, ...) {
  rlang::check_dots_empty()
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 1 || digits > 22 || digits != floor(digits))
    stop("`digits` must be a positive integer from 1 to 22.", call. = FALSE)
  cat("Exploratory bias screens\n")
  print(x$overview[, c("Interaction", "Rows", "Points", "ScreenAvailable",
    "ScreenUnavailable", "ScreenPositive", "MeanAbsBias", "Status")],
    digits = digits, row.names = FALSE)
  if (nrow(x$errors)) {
    cat("\nPair errors:\n")
    print(x$errors[, c("Interaction", "Error")], row.names = FALSE)
  }
  cat("\n", paste(x$notes, collapse = "\n"), "\n", sep = "")
  invisible(x)
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
    if (any(vapply(out[seq_len(i - 1L)],
        function(spec) setequal(spec$facets, pair), logical(1)))) {
      stop("Each facet combination may appear only once in `pairs`; reordered duplicates are also repeated screens.",
        call. = FALSE)
    }
    out[[i]] <- list(
      facets = pair,
      label = paste(pair, collapse = " x ")
    )
  }
  out
}
