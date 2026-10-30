# Reader-facing methods for portable fixed-calibration score batches.

mfrmr_validate_calibration_score_assumptions <- function(x) {
  s <- x$settings
  if (identical(s$schema_version, 6L)) prediction_validate_interval_record(x)
  if (!identical(s$estimator, "JML") && !identical(s$family, "GPCM") && is.null(s$scoring_prior) &&
      is.null(s$score_integration_review)) return(invisible(x))
  valid <- tryCatch({
    prior <- prediction_validate_scoring_prior(s$scoring_prior)
    retained <- s$retained_prior_identity
    actual <- s$prior_identity
    original <- s$semantic_components$scoring_prior
    r <- s$score_integration_review
    e <- x$estimates
    fields <- c("Person", "PriorSource", "PriorMean", "PriorSD", "RetainedPriorMean", "RetainedPriorSD", "ScoreIntegrationReady")
    all(fields %in% names(e)) &&
      identical(s$estimator, s$semantic_components$model_family_estimator$estimator) &&
      is.list(retained) && is.list(actual) && is.list(original) &&
      identical(retained, list(type = original$type, mean = original$prior_mean, sd = original$prior_sd)) &&
      identical(actual, if (is.null(prior)) retained else list(type = "user_supplied_normal", mean = prior$mean, sd = prior$sd)) &&
      all(e$PriorSource == if (is.null(prior)) "retained" else "user_supplied") &&
      all(e$PriorMean == actual$mean & e$PriorSD == actual$sd &
        e$RetainedPriorMean == retained$mean & e$RetainedPriorSD == retained$sd) &&
      (!identical(s$family, "GPCM") ||
        ((if (identical(s$schema_version,6L)) mfrmr_calibration_gmfrm_evidence_valid else if (identical(s$schema_version,5L)) mfrm_jml_scoring_evidence_valid else if (identical(s$estimator, "JML")) mfrmr_calibration_gpcm_jml_evidence_valid else mfrmr_calibration_source_evidence_valid)(s$source_scoring_evidence) &&
          identical(s$source_scoring_evidence, s$semantic_components$source_scoring_review$source_scoring_evidence))) &&
      (!identical(s$estimator, "JML") || identical(s$family, "GPCM") ||
        (mfrmr_calibration_jml_evidence_valid(s$source_scoring_evidence) &&
          identical(s$source_scoring_evidence, s$semantic_components$source_scoring_review$source_scoring_evidence))) &&
      is.data.frame(r) && identical(names(r), c("Person", "EAPChange", "SDChange", "ReferenceChange", "Passed")) &&
      identical(as.character(r$Person), as.character(e$Person)) &&
      identical(e$ScoreIntegrationReady, r$Passed) && all(r$Passed) &&
      all(vapply(r[c("EAPChange", "SDChange", "ReferenceChange")], function(v)
        is.numeric(v) && all(is.finite(v) & v >= 0 & v <= 1e-5), TRUE))
  }, error = function(e) FALSE)
  if (identical(s$schema_version,5L)) valid <- valid &&
    all(c("CalibrationMethod","CorrectionOrder") %in% names(x$estimates)) &&
    isTRUE(all(x$estimates$CalibrationMethod=="Corrected JML" & x$estimates$CorrectionOrder==s$source_scoring_evidence$local_calibration_review$correction_order))
  if (identical(s$schema_version,6L)) valid <- valid &&
    identical(s$semantic_components$schema$schema_version, 6L) &&
    identical(s$family, "GPCM") && identical(s$estimator, "MML") &&
    is.null(s$scoring_prior) && is.null(s$event_id_column) &&
    identical(s$prior_identity, list(type = "fixed_standard_normal", mean = 0, sd = 1)) &&
    identical(s$semantic_components$product_structure$slope_composition, "product") &&
    identical(s$scoring_algorithm, s$semantic_components$scoring_algorithm) &&
    identical(s$quadrature_order, s$semantic_components$quadrature$quadrature_order) &&
    identical(s$engine_identity, "artifact_coordinates_v6")
  if (!isTRUE(valid)) stop("Portable scoring output has missing or inconsistent prior or numerical-check records; score again from the frozen calibration.", call. = FALSE)
  invisible(x)
}

mfrmr_validate_calibration_score <- function(x, arg = "x") {
  mfrm_validate_decisions(x)
  if (!inherits(x, "mfrm_calibration_score")) {
    stop(
      "`", arg, "` must be output from score_mfrm_calibration().",
      call. = FALSE
    )
  }
  required <- c("estimates", "row_review", "person_dispositions", "settings")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L ||
      !is.data.frame(x$estimates) ||
      !is.data.frame(x$row_review) ||
      !is.data.frame(x$person_dispositions) ||
      !is.list(x$settings)) {
    stop(
      "`", arg, "` is not a complete portable calibration score result.",
      call. = FALSE
    )
  }
  mfrmr_validate_calibration_score_assumptions(x)
}

mfrmr_calibration_score_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L ||
      !is.finite(digits) || digits < 0 || digits != floor(digits)) {
    stop("`digits` must be one non-negative integer.", call. = FALSE)
  }
  as.integer(digits)
}

mfrmr_validate_calibration_score_summary <- function(x) {
  if (!inherits(x, "summary.mfrm_calibration_score")) {
    stop(
      "`x` must be output from summary() on a portable calibration score.",
      call. = FALSE
    )
  }
  required <- c(
    "overview", "estimates", "review", "row_review", "settings", "notes",
    "digits"
  )
  required_columns <- list(
    overview = c(
      "CalibrationId", "Model", "Estimator", "Scored", "Review",
      "NotScored"
    ),
    estimates = c(
      "Person", "Estimate", "SD", "Lower", "Upper", "Disposition"
    ),
    review = c("Person", "Disposition", "ReasonCodes"),
    row_review = c(
      "InputRows", "ScoredRows", "OmittedRows", "RefusedRows",
      "MissingResponsePolicy"
    )
  )
  missing <- setdiff(required, names(x))
  missing_columns <- if (length(missing) == 0L) {
    any(vapply(names(required_columns), function(name) {
      !all(required_columns[[name]] %in% names(x[[name]]))
    }, logical(1)))
  } else {
    TRUE
  }
  if (length(missing) > 0L ||
      !is.data.frame(x$overview) || nrow(x$overview) != 1L ||
      !is.data.frame(x$estimates) || !is.data.frame(x$review) ||
      !is.data.frame(x$row_review) || !is.list(x$settings) ||
      missing_columns) {
    stop(
      "`x` is not a complete portable calibration score summary.",
      call. = FALSE
    )
  }
  mfrmr_validate_calibration_score_assumptions(x)
}

mfrmr_calibration_score_print_value <- function(value) {
  value <- value[1L]
  if (length(value) == 0L || is.na(value)) return("NA")
  format(value, trim = TRUE, scientific = FALSE)
}

mfrmr_calibration_score_print_prior <- function(settings) {
  if (is.null(settings$retained_prior_identity)) return(invisible(NULL))
  prior <- settings$prior_identity
  retained <- settings$retained_prior_identity
  print_wrapped_line(paste0("Scoring prior: ", if (is.null(settings$scoring_prior)) "retained normal" else "user-supplied normal",
    " (mean ", prior$mean, ", SD ", prior$sd, "). Retained prior: mean ", retained$mean, ", SD ", retained$sd, "."))
  if (identical(settings$schema_version,5L)) print_wrapped_line(mfrm_jml_scoring_note(settings$source_scoring_evidence))
  if (identical(settings$schema_version,6L)) print_wrapped_line(
    "Experimental two-family EAP on the fixed N(0,1) scale; sampling coverage and population transport are not established.")
  if (identical(settings$estimator, "JML")) print_wrapped_line(
    "JML-calibrated post-hoc EAP: the reference prior was not estimated by JML; these are not ML/WLE scores.")
  if (identical(settings$family, "GPCM")) print_wrapped_line(
    "Source calibration: conditional checks passed at extraction; the source fit is not re-estimated on replay.")
  if (any(settings$source_scoring_evidence$source_audit_states %in% "not_evaluated"))
    print_wrapped_line("Global identification or boundary checks remain incomplete; these are conditional scores.")
  invisible(NULL)
}

mfrmr_calibration_score_interval_note <- function(settings) {
  algorithm <- settings$scoring_algorithm %||% ""
  level <- settings$interval_level
  label <- if (length(level) == 1L && is.finite(level)) {
    paste0(format(100 * level, trim = TRUE), "% intervals: ")
  } else "Intervals: "
  paste0(label, if (endsWith(algorithm, "_v1")) {
    "saved grid endpoints; posterior mass may differ from the requested level."
  } else if (endsWith(algorithm, "_v2")) {
    "continuous posterior quantiles."
  } else {
    "calculation method not recorded in this result."
  })
}

mfrmr_calibration_score_print_preview <- function(data, kind, limit = 10L) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (nrow(data) == 0L) return(invisible(NULL))
  shown <- utils::head(data, limit)
  cat("\n", kind, " (", nrow(shown), " of ", nrow(data), ")\n", sep = "")
  for (i in seq_len(nrow(shown))) {
    person <- as.character(shown$Person[i] %||% "<unknown>")
    disposition <- switch(as.character(shown$Disposition[i] %||% "unknown"),
      scored = "scored", scored_review = "review required",
      not_scored = "not scored", "status not recorded")
    if (identical(kind, "Posterior estimates")) {
      line <- paste0(
        person,
        ": estimate ", mfrmr_calibration_score_print_value(shown$Estimate[i]),
        ", SD ", mfrmr_calibration_score_print_value(shown$SD[i]),
        ", interval [", mfrmr_calibration_score_print_value(shown$Lower[i]),
        ", ", mfrmr_calibration_score_print_value(shown$Upper[i]), "], ",
        disposition
      )
    } else {
      reasons <- as.character(shown$ReasonCodes[i] %||% "")
      if (is.na(reasons) || !nzchar(reasons)) reasons <- "none recorded"
      labels <- c(
        RESPONSES_MISSING_OMITTED = "missing responses omitted",
        ALL_RESPONSES_LOWER_ENDPOINT = "all responses at the minimum score",
        ALL_RESPONSES_UPPER_ENDPOINT = "all responses at the maximum score",
        VERY_SPARSE_RESPONSE_PATTERN = "very few valid responses",
        QUADRATURE_EDGE_MASS_REVIEW = "posterior mass near the integration limits",
        ZERO_VALID_RESPONSES = "no valid responses"
      )
      codes <- strsplit(reasons, ";", fixed = TRUE)[[1L]]
      descriptions <- unname(labels[codes])
      descriptions[is.na(descriptions)] <- "see ReasonCodes in the review table"
      reasons <- paste(descriptions, collapse = ", ")
      line <- paste0(
        person, ": ", disposition,
        "; reasons: ", reasons
      )
    }
    print_wrapped_line(line)
  }
  if (nrow(data) > nrow(shown)) {
    print_wrapped_line(
      paste0(
        nrow(data) - nrow(shown),
        " additional row(s) remain in the summary object."
      )
    )
  }
  invisible(NULL)
}

mfrmr_calibration_score_overview <- function(x) {
  dispositions <- as.data.frame(
    x$person_dispositions, stringsAsFactors = FALSE
  )
  row_review <- as.data.frame(x$row_review, stringsAsFactors = FALSE)
  row_value <- function(name, default = 0L) {
    if (!name %in% names(row_review) || nrow(row_review) == 0L) return(default)
    value <- suppressWarnings(as.numeric(row_review[[name]][1L]))
    if (!is.finite(value)) default else value
  }
  disposition <- as.character(dispositions$Disposition %||% character(0))
  data.frame(
    CalibrationId = as.character(
      x$settings$calibration_id %||% NA_character_
    ),
    Model = as.character(x$settings$family %||% NA_character_),
    Estimator = if (identical(x$settings$schema_version,5L)) "Corrected JML" else as.character(x$settings$estimator %||% NA_character_),
    ScoringBasis = as.character(
      x$settings$scoring_basis %||% NA_character_
    ),
    IntervalLevel = suppressWarnings(as.numeric(
      x$settings$interval_level %||% NA_real_
    )),
    Persons = as.integer(nrow(dispositions)),
    Scored = as.integer(sum(disposition %in% c("scored", "scored_review"))),
    Review = as.integer(sum(disposition == "scored_review")),
    NotScored = as.integer(sum(disposition == "not_scored")),
    InputRows = as.integer(row_value("InputRows")),
    ScoredRows = as.integer(row_value("ScoredRows")),
    OmittedRows = as.integer(row_value("OmittedRows")),
    stringsAsFactors = FALSE
  )
}

#' Review and plot portable fixed-calibration scores
#'
#' These methods provide the first review surface for the result returned by
#' [score_mfrm_calibration()]. `print()` gives a compact batch disposition,
#' `summary()` exposes readable score and review tables, and `plot()` shows
#' conditional score uncertainty or numerical review quantities without
#' refitting a model.
#'
#' The default interval plot shows posterior EAP estimates and central
#' intervals. `type = "precision"` plots the number of valid response rows
#' against posterior SD. `type = "edge_mass"` compares posterior mass on the
#' two outer quadrature nodes with the recorded review threshold. Review
#' dispositions are highlighted in every view.
#'
#' These are score-batch review displays, not calibration-fit diagnostics.
#' Posterior SDs and intervals are conditional on the frozen point calibration
#' and the actual scoring prior. They exclude calibration-parameter uncertainty.
#' For GPCM or an explicitly supplied scoring prior, summaries and plot data
#' retain both original and actual prior values without rounding, the per-Person
#' integration checks, and (for GPCM) the source decision recorded at extraction.
#' These records describe conditional scoring, not a new evaluation of the
#' training fit or the suitability of its population for another cohort.
#' `summary(x)$settings$score_integration_review` compares the reported EAP/SD
#' with higher-order adaptive references. The separate `quadrature_review`
#' compares adaptive and fixed grids; for an adaptive scorer, that fixed grid
#' did not produce the reported scores. Printing distinguishes the two checks.
#' Persons with no valid responses have no score coordinate and are retained
#' in `summary(x)$review` and in the plot payload's
#' `unplotted_dispositions` component.
#'
#' The summary object retains every returned score and review disposition.
#' Its estimates table preserves the available estimate and uncertainty basis
#' columns, calibration identifiers, scoring algorithm and requested interval
#' level when extracted or written to CSV. The interval level is not rounded.
#' Older score results recover these fields from their recorded settings when
#' re-summarized. Saved grid-based intervals retain their original endpoints;
#' printed results and interval plots explain that their posterior mass can
#' differ from the requested level.
#' Its print method shows at most ten rows from each table so routine console
#' output stays compact. Base and ggplot2 renderers distinguish scored and
#' review states by shape as well as colour.
#' `as_ggplot(p)` and `as_ggplot(p, component = "data")` preserve the same
#' selected view, saved bounds, review symbols and unplotted dispositions.
#' Use [plot_data()] to extract other components. Changing the interval level
#' requires an explicit scoring call; conversion does not recompute intervals.
#' Notes remain in the returned object and are not automatically drawn.
#'
#' @param object,x An `mfrm_calibration_score` returned by
#'   [score_mfrm_calibration()].
#' @param digits Number of decimal places in the summary.
#' @param type One of `"interval"`, `"precision"`, or `"edge_mass"`.
#' @param top_n Maximum number of scored Persons shown. Review rows are
#'   selected first when truncation is necessary. Use `Inf` to show all.
#' @param sort_by Selection priority after review rows: absolute estimate,
#'   posterior SD, or Person identifier.
#' @param label_review Whether review points are labelled in the precision and
#'   edge-mass views. Interval plots always label the displayed Persons.
#' @param main Optional plot title.
#' @param draw If `TRUE`, draw with base R graphics. The returned
#'   `mfrm_plot_data` is available invisibly in either case.
#' @param preset Visual preset: `"standard"`, `"publication"`, `"compact"`,
#'   or `"monochrome"`.
#' @param ... Reserved for generic compatibility.
#'
#' @return `print()` returns its input invisibly. `summary()` returns a
#'   `summary.mfrm_calibration_score` containing
#'   `overview`, `estimates`, `review`, `row_review`, `settings`, and `notes`.
#'   When adaptive integration was requested and scored rows exist, it also
#'   retains unrounded `quadrature_review` and a compact `quadrature_overview`.
#'   When cut scores were requested, `decisions` preserves every probability,
#'   cut boundary and actual scoring prior without rounding; unavailable
#'   Persons remain in that table with `NA` probabilities.
#'   `plot()` returns an `mfrm_plot_data` with the selected plotting table,
#'   selection accounting, unplotted Person dispositions, interpretation
#'   guidance, and plotting settings.
#'
#' @seealso [mfrm_calibration_workflow], [mfrm_calibration_capabilities],
#'   [plot_data()], [plot_data_components()], [as_ggplot()], [apa_table()],
#'   [build_summary_table_bundle()], [export_summary_appendix()],
#'   [predict_mfrm_units()] for conditional-decision interpretation
#' @name mfrm_calibration_score_methods
#'
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("mfrm_calibration_score_methods", run.dontrun = TRUE)
#' \dontrun{
#' dat <- load_mfrmr_data("example_core")
#' ids <- unique(dat$Person)
#' training <- dat[dat$Person %in% ids[1:18], , drop = FALSE]
#' fit <- fit_mfrm(
#'   training, "Person", c("Rater", "Criterion"), "Score",
#'   model = "RSM", method = "MML", quad_points = 31
#' )
#' q_review <- mml_quadrature_sensitivity(
#'   fit, training, quad_points = c(31, 61), theta_points = 41
#' )
#' summary(q_review) # Inspect movement; these grids are not universal cutoffs.
#' fit <- q_review$fits$q61
#' calibration <- freeze_mfrm_calibration(
#'   validate_mfrm_calibration(extract_mfrm_calibration(
#'     fit, quadrature_review = q_review, scoring_quad_points = 121
#'   ))
#' )
#' new_rows <- dat[dat$Person %in% ids[19:20], , drop = FALSE]
#' scores <- score_mfrm_calibration(calibration, new_rows,
#'   cut_scores = c(Basic = -0.5, Advanced = 0.8))
#' summary(scores)
#' apa_table(scores, which = "decisions")
#' build_summary_table_bundle(scores)$table_index
#' plot(scores, type = "interval", draw = FALSE)
#' plot(scores, type = "edge_mass", draw = FALSE)
#' }
NULL

#' @rdname mfrm_calibration_score_methods
#' @method summary mfrm_calibration_score
#' @export
summary.mfrm_calibration_score <- function(object, digits = 3L, ...) {
  mfrmr_validate_calibration_score(object, "object")
  digits <- mfrmr_calibration_score_digits(digits)

  estimates <- as.data.frame(object$estimates, stringsAsFactors = FALSE)
  if (!"ScoringAlgorithm" %in% names(estimates)) {
    estimates$ScoringAlgorithm <- rep(
      object$settings$scoring_algorithm %||% NA_character_, nrow(estimates)
    )
  }
  if (!"IntervalLevel" %in% names(estimates)) {
    estimates$IntervalLevel <- rep(
      object$settings$interval_level %||% NA_real_, nrow(estimates)
    )
  }
  estimate_columns <- intersect(
    c(
      "Person", "Estimate", "SD", "Lower", "Upper", "Observations",
      "WeightedN", "Disposition", "ReasonCodes", "ReadinessStatus",
      "EstimateBasis", "UncertaintyBasis", "CalibrationId", "SchemaVersion",
      "ScoringBasis", "ScoringAlgorithm", "IntervalLevel", "PriorSource", "PriorMean", "PriorSD",
      "RetainedPriorMean", "RetainedPriorSD", "ScoreIntegrationReady",
      "CalibrationMethod", "CorrectionOrder"
    ),
    names(estimates)
  )
  estimates <- estimates[, estimate_columns, drop = FALSE]

  dispositions <- as.data.frame(
    object$person_dispositions, stringsAsFactors = FALSE
  )
  review <- dispositions[
    as.character(dispositions$Disposition) != "scored",
    ,
    drop = FALSE
  ]
  review_columns <- intersect(
    c(
      "Person", "Disposition", "ReasonCodes", "ValidResponses",
      "OmittedResponses", "EndpointStatus", "VerySparsePattern",
      "QuadratureEdgeMass", "QuadratureEdgeThreshold", "ReadinessStatus",
      "AdministrationCompleteness"
    ),
    names(review)
  )
  review <- review[, review_columns, drop = FALSE]

  round_numeric <- function(data) {
    numeric_columns <- vapply(data, is.double, logical(1))
    numeric_columns[names(data) %in% c("IntervalLevel", "PriorMean", "PriorSD", "RetainedPriorMean", "RetainedPriorSD")] <- FALSE
    data[numeric_columns] <- lapply(
      data[numeric_columns], round, digits = digits
    )
    data
  }

  out <- list(
    overview = round_numeric(mfrmr_calibration_score_overview(object)),
    estimates = round_numeric(estimates),
    review = round_numeric(review),
    row_review = round_numeric(as.data.frame(
      object$row_review, stringsAsFactors = FALSE
    )),
    settings = list(
      calibration_id = object$settings$calibration_id,
      family = object$settings$family,
      estimator = object$settings$estimator,
      scoring_basis = object$settings$scoring_basis,
      scoring_algorithm = object$settings$scoring_algorithm,
      quadrature_order = object$settings$quadrature_order,
      interval_level = object$settings$interval_level,
      uncertainty_basis = if (nrow(estimates) > 0L &&
          "UncertaintyBasis" %in% names(object$estimates)) {
        unique(as.character(object$estimates$UncertaintyBasis))
      } else {
        "conditional_on_frozen_point_calibration"
      }
    ),
    notes = object$notes %||% character(0),
    digits = digits
  )
  if (!is.null(object$quadrature_review)) {
    out$quadrature_review <- object$quadrature_review
    out$quadrature_overview <- mfrmr_adaptive_quadrature_overview(object$quadrature_review)
  }
  if (!is.null(object$settings$score_integration_review)) {
    for (name in c("schema_version", "semantic_components", "source_scoring_evidence", "score_integration_review",
                    "scoring_prior", "prior_identity", "retained_prior_identity",
                    "engine_identity", "event_id_column")) {
      out$settings[name] <- object$settings[name]
    }
  }
  class(out) <- c("summary.mfrm_calibration_score", "list")
  if (!is.null(object$decisions)) {
    out$decisions <- object$decisions
    out$settings$cut_scores <- object$settings$cut_scores
  }
  out
}

#' @rdname mfrm_calibration_score_methods
#' @method print mfrm_calibration_score
#' @export
print.mfrm_calibration_score <- function(x, ...) {
  summary_object <- summary.mfrm_calibration_score(x)
  overview <- summary_object$overview[1L, , drop = FALSE]
  cat("<mfrm_calibration_score>\n")
  print_wrapped_line(paste0("Calibration: ", overview$CalibrationId))
  print_wrapped_line(
    paste0("Model / estimator: ", overview$Model, " / ", overview$Estimator)
  )
  cat(
    "  Persons: ", overview$Scored, " scored (", overview$Review,
    " requiring review); ", overview$NotScored, " not scored\n", sep = ""
  )
  cat(
    "  Rows: ", overview$ScoredRows, " scored; ", overview$OmittedRows,
    " omitted\n", sep = ""
  )
  print_wrapped_line("Intervals: conditional on frozen point calibration;")
  print_wrapped_line("calibration-parameter uncertainty excluded")
  print_wrapped_line(mfrmr_calibration_score_interval_note(x$settings))
  mfrmr_calibration_score_print_prior(x$settings)
  if (overview$Scored > 0L) {
    print_wrapped_line(
      "Use `summary(x)` for review tables and `plot(x)` for score intervals.",
      prefix = ""
    )
  } else {
    print_wrapped_line(
      "No score coordinate is available; inspect `summary(x)$review`.",
      prefix = ""
    )
  }
  invisible(x)
}

#' @rdname mfrm_calibration_score_methods
#' @method print summary.mfrm_calibration_score
#' @export
print.summary.mfrm_calibration_score <- function(x, ...) {
  mfrm_validate_decisions(x)
  mfrmr_validate_calibration_score_summary(x)
  overview <- x$overview[1L, , drop = FALSE]
  cat("mfrmr Portable Calibration Score Summary\n")
  mfrmr_calibration_score_print_prior(x$settings)
  checks <- x$settings$score_integration_review
  if (!is.null(checks)) {
    cat("\nReported-score integration check\n")
    if (nrow(checks)) {
      print_wrapped_line(sprintf(paste0(
        "%d of %d scored Persons passed. Against higher-order adaptive integration: ",
        "largest absolute EAP change %.3g; posterior SD change %.3g; reference change %.3g."),
        sum(checks$Passed), nrow(checks), max(checks$EAPChange),
        max(checks$SDChange), max(checks$ReferenceChange)))
    } else print_wrapped_line("No scored Persons; no score integration comparison is available.")
  }
  if (!is.null(x$quadrature_overview)) {
    cat("\nFixed-parameter integration review (adaptive minus fixed)\n")
    if (startsWith(x$settings$scoring_algorithm %||% "", "adaptive_"))
      print_wrapped_line(paste(
        "This compares alternative grids. The fixed grid in this table did not produce",
        "the reported scores, which use adaptive integration."))
    print(x$quadrature_overview, row.names = FALSE)
  }
  print_wrapped_line(paste0("Calibration: ", overview$CalibrationId))
  print_wrapped_line(
    paste0("Model / estimator: ", overview$Model, " / ", overview$Estimator)
  )
  print_wrapped_line(
    paste0(
      "Persons: ", overview$Scored, " scored (", overview$Review,
      " requiring review); ", overview$NotScored, " not scored"
    )
  )
  mfrmr_calibration_score_print_preview(
    x$estimates, "Posterior estimates"
  )
  mfrmr_calibration_score_print_preview(
    x$review, "Persons requiring review or not scored"
  )
  if (nrow(x$row_review) > 0L) {
    cat("\nResponse-row disposition\n")
    row <- x$row_review[1L, , drop = FALSE]
    print_wrapped_line(paste0(
      mfrmr_calibration_score_print_value(row$InputRows), " input; ",
      mfrmr_calibration_score_print_value(row$ScoredRows), " scored; ",
      mfrmr_calibration_score_print_value(row$OmittedRows), " omitted; ",
      mfrmr_calibration_score_print_value(row$RefusedRows), " refused; ",
      "missing-response policy: ",
      as.character(row$MissingResponsePolicy %||% "not recorded")
    ))
  }
  if (!is.null(x$decisions)) {
    cat("\nConditional decision probabilities (first 10; full table in $decisions)\n")
    print(utils::head(x$decisions, 10L), row.names=FALSE)
    print_wrapped_line(mfrm_decision_note())
  }
  cat("\nInterval interpretation\n")
  print_wrapped_line(mfrmr_calibration_score_interval_note(x$settings))
  print_wrapped_line(
    paste(
      "Posterior SDs and intervals are conditional on the frozen point",
      "calibration and exclude calibration-parameter uncertainty."
    )
  )
  invisible(x)
}

mfrmr_calibration_score_top_n <- function(top_n) {
  if (is.numeric(top_n) && length(top_n) == 1L && is.infinite(top_n) &&
      top_n > 0) {
    return(Inf)
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
      top_n < 1 || top_n != floor(top_n)) {
    stop("`top_n` must be one positive integer or `Inf`.", call. = FALSE)
  }
  as.integer(top_n)
}

mfrmr_calibration_score_plot_table <- function(x, top_n, sort_by) {
  estimates <- as.data.frame(x$estimates, stringsAsFactors = FALSE)
  required <- c(
    "Person", "Estimate", "SD", "Lower", "Upper", "Observations",
    "Disposition", "ReasonCodes", "ReadinessStatus"
  )
  missing <- setdiff(required, names(estimates))
  if (length(missing) > 0L) {
    stop(
      "The score result is missing plot field(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (nrow(estimates) == 0L) {
    stop(
      "No scored Person coordinate is available to plot; inspect ",
      "`summary(x)$review`.",
      call. = FALSE
    )
  }
  numeric_fields <- c("Estimate", "SD", "Lower", "Upper", "Observations")
  estimates[numeric_fields] <- lapply(
    estimates[numeric_fields], function(value) suppressWarnings(as.numeric(value))
  )
  if (any(!is.finite(as.matrix(estimates[numeric_fields])))) {
    stop("The score result contains non-finite plot coordinates.", call. = FALSE)
  }

  dispositions <- as.data.frame(
    x$person_dispositions, stringsAsFactors = FALSE
  )
  disposition_index <- match(estimates$Person, dispositions$Person)
  for (field in c(
    "EndpointStatus", "VerySparsePattern", "QuadratureEdgeMass",
    "QuadratureEdgeThreshold", "AdministrationCompleteness"
  )) {
    estimates[[field]] <- if (field %in% names(dispositions)) {
      dispositions[[field]][disposition_index]
    } else {
      NA
    }
  }
  estimates$QuadratureEdgeMass <- suppressWarnings(as.numeric(
    estimates$QuadratureEdgeMass
  ))
  estimates$QuadratureEdgeThreshold <- suppressWarnings(as.numeric(
    estimates$QuadratureEdgeThreshold
  ))
  estimates$ReviewFlag <- as.character(estimates$Disposition) != "scored" |
    as.character(estimates$ReadinessStatus) != "conditional_score_ready"
  estimates$ReviewLabel <- ifelse(
    estimates$ReviewFlag, as.character(estimates$Person), ""
  )

  priority <- switch(
    sort_by,
    estimate = -abs(estimates$Estimate),
    sd = -estimates$SD,
    person = xtfrm(as.character(estimates$Person))
  )
  selection_order <- order(!estimates$ReviewFlag, priority, method = "radix")
  selected <- if (is.infinite(top_n) || nrow(estimates) <= top_n) {
    seq_len(nrow(estimates))
  } else {
    selection_order[seq_len(top_n)]
  }
  estimates$Selected <- seq_len(nrow(estimates)) %in% selected
  selected_table <- estimates[selected, , drop = FALSE]
  display_order <- order(
    selected_table$Estimate, as.character(selected_table$Person),
    method = "radix"
  )
  selected_table <- selected_table[display_order, , drop = FALSE]
  selected_table$DisplayOrder <- seq_len(nrow(selected_table))

  list(full = estimates, selected = selected_table)
}

#' @rdname mfrm_calibration_score_methods
#' @method plot mfrm_calibration_score
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_calibration_score <- function(
    x,
    type = c("interval", "precision", "edge_mass"),
    top_n = 40L,
    sort_by = c("estimate", "sd", "person"),
    label_review = TRUE,
    main = NULL,
    draw = TRUE,
    preset = c("standard", "publication", "compact", "monochrome"),
    ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  mfrmr_validate_calibration_score(x)
  type <- match.arg(type)
  sort_by <- match.arg(sort_by)
  top_n <- mfrmr_calibration_score_top_n(top_n)
  if (!is.logical(label_review) || length(label_review) != 1L ||
      is.na(label_review)) {
    stop("`label_review` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw)) {
    stop("`draw` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.null(main) &&
      (!is.character(main) || length(main) != 1L || is.na(main) ||
       !nzchar(main))) {
    stop("`main` must be `NULL` or one non-empty string.", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  tables <- mfrmr_calibration_score_plot_table(x, top_n, sort_by)
  plotted <- tables$selected
  full <- tables$full
  plotted$LabelPosition <- 3L
  plotted$LabelX <- NA_real_
  plotted$LabelY <- NA_real_
  review_index <- which(plotted$ReviewFlag)
  if (length(review_index) > 0L && identical(type, "precision")) {
    plotted$LabelX[review_index] <- plotted$Observations[review_index]
    sd_span <- diff(range(plotted$SD, finite = TRUE))
    if (!is.finite(sd_span) || sd_span <= 0) sd_span <- 0.1
    label_spacing <- max(sd_span * 0.07, 0.01)
    label_order <- review_index[order(
      plotted$SD[review_index], as.character(plotted$Person[review_index]),
      method = "radix"
    )]
    previous <- -Inf
    for (row in label_order) {
      label_y <- max(plotted$SD[row] + label_spacing,
                     previous + label_spacing)
      plotted$LabelY[row] <- label_y
      previous <- label_y
    }
  }
  if (length(review_index) > 0L && identical(type, "edge_mass")) {
    midpoint <- mean(range(plotted$Estimate, finite = TRUE))
    plotted$LabelPosition[review_index] <- ifelse(
      plotted$Estimate[review_index] <= midpoint, 4L, 2L
    )
    estimate_span <- diff(range(plotted$Estimate, finite = TRUE))
    if (!is.finite(estimate_span) || estimate_span <= 0) estimate_span <- 1
    label_nudge <- estimate_span * 0.015
    plotted$LabelX[review_index] <- plotted$Estimate[review_index] +
      ifelse(plotted$LabelPosition[review_index] == 4L,
             label_nudge, -label_nudge)
    plotted$LabelY[review_index] <- plotted$QuadratureEdgeMass[review_index]
  }
  dispositions <- as.data.frame(
    x$person_dispositions, stringsAsFactors = FALSE
  )
  not_scored <- dispositions[
    as.character(dispositions$Disposition) == "not_scored",
    ,
    drop = FALSE
  ]
  review_available <- sum(full$ReviewFlag)
  review_plotted <- sum(plotted$ReviewFlag)
  title <- as.character(main %||% switch(
    type,
    interval = "Portable calibration score intervals",
    precision = "Conditional score precision",
    edge_mass = "Quadrature edge-mass review"
  ))
  subtitle <- paste0(
    nrow(plotted), " of ", nrow(full), " scored Persons shown; ",
    review_plotted, " review disposition",
    if (review_plotted == 1L) "" else "s", " shown"
  )
  if (nrow(not_scored)) subtitle <- paste0(subtitle, "; ", nrow(not_scored), " not scored")
  uncertainty_note <- paste(
    "Conditional on the frozen point calibration;",
    "calibration-parameter uncertainty excluded."
  )
  if (identical(type, "interval")) {
    uncertainty_note <- paste(
      mfrmr_calibration_score_interval_note(x$settings), uncertainty_note,
      sep = "\n"
    )
  }
  if (identical(x$settings$schema_version,5L)) uncertainty_note <- paste(
    uncertainty_note,paste0("Corrected JML (order ",x$settings$source_scoring_evidence$local_calibration_review$correction_order,
      "); residual calibration bias may remain."),sep="\n")
  if (identical(x$settings$schema_version,6L)) uncertainty_note <- paste(
    uncertainty_note, "Experimental two-family model; sampling coverage and population transport are not established.", sep = "\n")
  selection_summary <- data.frame(
    ScoredPersons = as.integer(nrow(full)),
    PlottedPersons = as.integer(nrow(plotted)),
    OmittedFromPlot = as.integer(nrow(full) - nrow(plotted)),
    ReviewPersonsAvailable = as.integer(review_available),
    ReviewPersonsPlotted = as.integer(review_plotted),
    NotScoredPersons = as.integer(nrow(not_scored)),
    ReviewRowsPrioritized = TRUE,
    SelectionPriority = sort_by,
    stringsAsFactors = FALSE
  )
  display_encoding <- data.frame(
    Status = c("Scored", "Review"),
    ReviewFlag = c(FALSE, TRUE),
    Colour = c(style$accent_primary, style$warn),
    Shape = c(16L, 17L),
    stringsAsFactors = FALSE
  )
  legend <- new_plot_legend(
    label = c("Scored", "Review"),
    role = c("conditional score ready", "inspect disposition and reason codes"),
    aesthetic = c("colour and shape", "colour and shape"),
    value = paste0(
      display_encoding$Colour, "; pch=", display_encoding$Shape
    )
  )
  edge_threshold <- unique(plotted$QuadratureEdgeThreshold[
    is.finite(plotted$QuadratureEdgeThreshold)
  ])
  reference_lines <- switch(
    type,
    interval = new_reference_lines(
      axis = "v", value = 0, label = "Scale origin",
      linetype = "dashed", role = "reference"
    ),
    precision = new_reference_lines(),
    edge_mass = if (length(edge_threshold) > 0L) {
      new_reference_lines(
        axis = "h", value = edge_threshold[1L],
        label = "Recorded edge-mass review threshold",
        linetype = "dashed", role = "review threshold"
      )
    } else {
      new_reference_lines()
    }
  )

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    point_color <- ifelse(
      plotted$ReviewFlag, style$warn, style$accent_primary
    )
    point_shape <- ifelse(plotted$ReviewFlag, 17L, 16L)
    if (identical(type, "interval")) {
      x_range <- range(c(plotted$Lower, plotted$Upper, 0), finite = TRUE)
      if (diff(x_range) <= 0) x_range <- x_range + c(-0.5, 0.5)
      old_mar <- graphics::par("mar")
      on.exit(graphics::par(mar = old_mar), add = TRUE)
      graphics::par(mar = c(6.4, 8.2, 4.2, 1.5))
      graphics::plot(
        plotted$Estimate, plotted$DisplayOrder,
        xlim = x_range, ylim = c(0.5, nrow(plotted) + 0.5),
        yaxt = "n", ylab = "", xlab = "Posterior EAP (logits)",
        main = "", pch = point_shape, col = point_color
      )
      graphics::abline(v = 0, lty = 2, col = style$neutral)
      graphics::segments(
        plotted$Lower, plotted$DisplayOrder,
        plotted$Upper, plotted$DisplayOrder,
        col = point_color, lwd = 1.5
      )
      graphics::points(
        plotted$Estimate, plotted$DisplayOrder,
        pch = point_shape, col = point_color
      )
      graphics::axis(
        2, at = plotted$DisplayOrder, labels = plotted$Person,
        las = 1, cex.axis = style$axis_cex
      )
    } else if (identical(type, "precision")) {
      x_range <- range(plotted$Observations, finite = TRUE)
      if (diff(x_range) <= 0) {
        padding <- max(1, abs(x_range[1L]) * 0.06)
        x_range <- x_range + c(-padding, padding)
      }
      y_range <- range(c(plotted$SD, plotted$LabelY), finite = TRUE)
      if (diff(y_range) <= 0) y_range <- y_range + c(-0.05, 0.05)
      y_padding <- max(diff(y_range) * 0.14, 0.01)
      y_range <- y_range + c(-y_padding * 0.2, y_padding)
      graphics::plot(
        plotted$Observations, plotted$SD,
        xlim = x_range, ylim = y_range,
        xlab = "Valid response rows", ylab = "Posterior SD (logits)",
        main = "", pch = point_shape, col = point_color
      )
      if (isTRUE(label_review) && any(plotted$ReviewFlag)) {
        graphics::segments(
          plotted$Observations[plotted$ReviewFlag],
          plotted$SD[plotted$ReviewFlag],
          plotted$LabelX[plotted$ReviewFlag],
          plotted$LabelY[plotted$ReviewFlag],
          col = style$warn, lwd = 0.7
        )
        graphics::text(
          plotted$LabelX[plotted$ReviewFlag],
          plotted$LabelY[plotted$ReviewFlag],
          labels = plotted$Person[plotted$ReviewFlag],
          pos = 3, offset = 0.35, cex = 0.75, col = style$warn
        )
      }
    } else {
      edge_mass <- plotted$QuadratureEdgeMass
      edge_mass[!is.finite(edge_mass)] <- NA_real_
      maximum <- max(c(edge_mass, edge_threshold, 0), na.rm = TRUE)
      if (!is.finite(maximum) || maximum <= 0) maximum <- 1
      graphics::plot(
        plotted$Estimate, edge_mass,
        ylim = c(0, maximum * 1.08),
        xlab = "Posterior EAP (logits)",
        ylab = "Posterior mass at outer quadrature nodes",
        main = "", pch = point_shape, col = point_color
      )
      if (length(edge_threshold) > 0L) {
        graphics::abline(h = edge_threshold[1L], lty = 2, col = style$warn)
      }
      if (isTRUE(label_review) && any(plotted$ReviewFlag)) {
        graphics::text(
          plotted$Estimate[plotted$ReviewFlag],
          edge_mass[plotted$ReviewFlag],
          labels = plotted$Person[plotted$ReviewFlag],
          pos = plotted$LabelPosition[plotted$ReviewFlag],
          offset = 0.6, cex = 0.75, col = style$warn
        )
      }
    }
    graphics::legend(
      "topleft", legend = c("Scored", "Review"),
      pch = c(16L, 17L), col = c(style$accent_primary, style$warn),
      bty = "o", bg = style$background, box.col = style$background,
      cex = 0.78
    )
  }

  interpretation <- data.frame(
    Topic = c(
      "Review color", "Not scored", "Interval basis", "Use boundary"
    ),
    Guidance = c(
      "Inspect ReasonCodes before using a highlighted score in a decision.",
      "Persons without a score coordinate remain in unplotted_dispositions.",
      uncertainty_note,
      "This display reviews a score batch; it does not establish calibration fit or validity."
    ),
    stringsAsFactors = FALSE
  )
  out <- new_mfrm_plot_data(
    paste0("calibration_score_", type),
    list(
      data = plotted,
      selection_summary = selection_summary,
      unplotted_dispositions = not_scored,
      interpretation = interpretation,
      display_encoding = display_encoding,
      settings = list(
        type = type,
        top_n = top_n,
        sort_by = sort_by,
        label_review = label_review,
        interval_level = x$settings$interval_level,
        scoring_algorithm = x$settings$scoring_algorithm,
        uncertainty_basis = "conditional_on_frozen_point_calibration",
        calibration_id = x$settings$calibration_id,
        preset = style$name
      ),
      title = title,
      subtitle = subtitle,
      caption = uncertainty_note,
      legend = legend,
      reference_lines = reference_lines
    )
  )
  for (name in c("prior_identity", "retained_prior_identity", "scoring_prior",
                 "source_scoring_evidence", "score_integration_review")) {
    if (!is.null(x$settings[[name]])) out$data$settings[[name]] <- x$settings[[name]]
  }
  invisible(out)
}
