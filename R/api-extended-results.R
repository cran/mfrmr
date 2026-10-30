# Stored-result reporting for the two separate latent-effect RSM classes.
# No fitting, scoring, diagnostics or resampling takes place in this adapter.

mfrm_extended_fit <- function(x) inherits(x, c("mfrm_testlet", "mfrm_random_rater"))

mfrm_extended_prediction_source <- function(fit) {
  list(class = class(fit), columns = fit$input$columns, levels = fit$input$levels,
    basis = fit$input$basis, score_levels = fit$input$score_levels,
    calibration = fit$calibration, parameters = fit$parameters,
    settings = fit$settings, checks = fit$checks,
    rater_mode = fit$rater_mode, conditional_rater_covariance = fit$conditional_rater_covariance)
}

mfrm_extended_settings_table <- function(settings) {
  data.frame(Setting = names(settings), Value = vapply(settings, function(z) {
    if (is.null(z)) "NULL" else paste(as.character(z), collapse = ", ")
  }, character(1)), row.names = NULL)
}

mfrm_extended_results <- function(fit, include, predictions, intervals, scores = NULL, comparison = NULL, diagnostics = NULL,
    calibration_intervals = "none", calibration_level = .95) {
  if (!is.data.frame(fit$calibration_table) || is.null(fit$input$columns) ||
      is.null(fit$settings) || is.null(fit$checks)) {
    stop("Supply a complete fit_mfrm_testlet() or fit_mfrm_random_rater() result.", call. = FALSE)
  }
  testlet <- inherits(fit, "mfrm_testlet")
  if (!is.null(diagnostics) && (!inherits(diagnostics, "mfrm_response_diagnostics") ||
      !identical(diagnostics$source, mfrm_extended_prediction_source(fit)) ||
      !identical(diagnostics$source_data, fit$input$assigned_data %||% fit$input$data) ||
      !identical(diagnostics$source_observed, fit$input$data))) {
    stop("`diagnostics` must be saved mfrm_response_diagnostics() output with matching calibration and the exact source roster.", call. = FALSE)
  }
  if (!is.null(predictions)) {
    correct_type <- if (testlet) inherits(predictions, "mfrm_testlet_scores") else
      is.list(predictions) && is.matrix(predictions$probabilities) && is.data.frame(predictions$expected_scores)
    if (!correct_type || !identical(predictions$source, mfrm_extended_prediction_source(fit))) {
      stop("`predictions` must come from predict() with this calibration and matching source metadata. Regenerate older predictions from the saved fit, or omit predictions.", call. = FALSE)
    }
  }
  if (!is.null(scores) && (testlet || !inherits(scores, "mfrm_random_rater_scores") ||
      !identical(scores$source, mfrm_extended_prediction_source(fit)))) {
    stop("`scores` must come from score_mfrm_random_rater() with matching calibration and source metadata.", call. = FALSE)
  }
  if (!is.null(intervals) && (testlet || !inherits(intervals, "mfrm_random_rater_intervals") ||
      !identical(intervals$source, fit) || !is.matrix(intervals$intervals))) {
    stop("`intervals` must be a saved mfrm_random_rater_intervals() result from this exact random-rater fit.", call. = FALSE)
  }
  if (!is.null(comparison)) {
    if (!inherits(comparison, "mfrm_extended_comparison") ||
        !identical(comparison$source, mfrm_extended_prediction_source(fit)) ||
        !identical(comparison$source_events, mfrm_compare_events(fit$input$data, comparison$event_columns))) {
      stop("`comparison` must come from compare_mfrm() with this calibration, settings and observed rating events.", call. = FALSE)
    }
    assigned <- fit$input$assigned_data %||% fit$input$data
    omitted <- assigned[is.na(assigned[[fit$input$columns$score]]), , drop = FALSE]
    if (!identical(comparison$omitted_events, mfrm_compare_events(omitted, comparison$event_columns))) stop(
      "`comparison` has different omitted rating events.", call. = FALSE)
  }
  numerical <- isTRUE(fit$checks$NumericalReady) && isTRUE(fit$checks$InformationPositive)
  boundary <- isTRUE(fit$checks$EstimatedVarianceBoundary) ||
    isTRUE(fit$checks$EstimatedPersonVarianceBoundary)
  model <- if (testlet) "Testlet RSM" else "Shared-rater RSM"
  interpretation <- if (!numerical) "Numerical review required" else if (boundary)
    "Estimated variance boundary; regular intervals restricted" else
      "Numerical checks satisfied; model adequacy not assessed"
  notes <- c("Numerical checks are not model-fit diagnostics or evidence of rater quality.",
    "No fitting, scoring, diagnostics or resampling is performed while collecting or exporting these results.",
    "Missing assigned scores and absent assignments are different: omitted rows are recorded; absent rows are not imputed.",
    if (testlet) "Person intervals condition on calibration and the fitted or specified normal ability population; their estimation uncertainty and general coverage guarantees are not included. Prior-only rows are not measured abilities." else
      "Individual-rater intervals are not supplied automatically; nominal coverage remains unresolved. PredictionSE is a first-order approximation, not a rater-quality classification.",
    if (!testlet) "Probability predictions use supplied abilities; each row is marginal, not a joint rating distribution. Saved Person scores use joint conditional rater Laplace integration, holding calibration fixed; numerical agreement does not certify approximation accuracy or coverage.",
    "Ordinary residual diagnostics, response-MI pooling, portable calibration and the Shiny viewer are unavailable for these model classes. Extended Wright/pathway maps use explicitly matched saved scores and posterior diagnostics.")
  decision <- data.frame(Interpretation = interpretation,
    FormalInference = "Target-specific limits; no general clearance", FitReadiness = "Not assessed",
    Why = notes[1], NextAction = "Review numerical checks, interval meanings and the rating design.")
  tables <- list(
    interpretation = decision,
    calibration = mfrm_extended_calibration_table(fit, calibration_intervals, calibration_level),
    numerical_checks = as.data.frame(fit$checks),
    data_usage = data.frame(Input = fit$input$input_rows, Observed = nrow(fit$input$data),
      Omitted = length(fit$input$omitted_rows)),
    omitted_rows = data.frame(InputRow = fit$input$omitted_rows),
    model_settings = mfrm_extended_settings_table(fit$settings),
    column_roles = data.frame(Role = rep(names(fit$input$columns), lengths(fit$input$columns)),
      Column = unlist(fit$input$columns, use.names = FALSE)),
    interpretation_notes = data.frame(Note = notes))
  fixed_variance <- if (testlet) fit$settings$fixed_variance else fit$settings$fixed_rater_sd
  tables$variance <- data.frame(Effect = if (testlet) "Person-specific testlet" else "Shared rater",
    Variance = if (testlet) fit$calibration$variance else fit$calibration$rater_sd^2,
    Estimated = is.null(fixed_variance), EstimatedBoundary = isTRUE(fit$checks$EstimatedVarianceBoundary))
  basis <- data.frame(Target = c("Fixed-facet and step calibration", "Variance component"),
    Interval = c(if (calibration_intervals == "none") "No automatic calibration interval" else
      paste0("Explicit ", format(100 * calibration_level, trim = TRUE, digits = 15), "% observed-information normal approximation"), "No regular variance interval supplied"),
    Limitation = c("SEs use observed information. Explicit normal bounds require resolved numerical/information checks and interior estimated variances; nominal coverage is not established. Missing endpoints remain unavailable.",
      "An estimated zero variance is a boundary, not proof that dependence is absent."))
  tables$variance <- rbind(tables$variance, data.frame(Effect = "Person ability",
      Variance = fit$calibration$person_variance %||% 1,
      Estimated = !is.null(fit$calibration$person_variance) && is.null(fit$settings$fixed_person_sd),
      EstimatedBoundary = isTRUE(fit$checks$EstimatedPersonVarianceBoundary)))
  if (testlet) {
    tables$blocks <- fit$input$blocks
  } else {
    tables$raters <- mfrm_random_rater_table(fit)
    tables$rater_workload <- fit$input$workload
    basis <- rbind(basis, data.frame(Target = "Observed-rater severity",
      Interval = "No automatic individual-rater interval",
      Limitation = "PredictionSE is a first-order approximation; nominal coverage unresolved. Explicit normal or saved bootstrap calculations require their own interpretation."))
  }
  if (!is.null(predictions)) {
    tables$prediction_settings <- mfrm_extended_settings_table(predictions$settings)
    if (!testlet) {
      tables$prediction_inputs <- predictions$newdata
      tables$expected_scores <- predictions$expected_scores
      probabilities <- predictions$probabilities
      tables$category_probabilities <- data.frame(Row = rep(seq_len(nrow(probabilities)), ncol(probabilities)),
        Score = rep(fit$input$score_levels, each = nrow(probabilities)), Probability = as.vector(probabilities))
    }
  }
  person_scores <- if (testlet) predictions else scores
  if (!is.null(person_scores)) {
    tables$person_scores <- person_scores$table
    tables$scoring_settings <- mfrm_extended_settings_table(person_scores$settings)
    if (testlet) tables$scoring_blocks <- person_scores$blocks
    # Keep the supplied event roster as well as aggregated block counts.
    # Older saved scores without a roster remain readable; do not reconstruct
    # their scoring responses from the calibration data.
    if (is.data.frame(person_scores$scoring_data)) tables$scoring_roster <- person_scores$scoring_data
    tables$scoring_data_usage <- as.data.frame(as.list(person_scores$data_usage))
    tables$scoring_omitted_rows <- data.frame(InputRow = person_scores$omitted_rows)
    tables$scoring_status <- as.data.frame(table(factor(person_scores$table$Status,
      levels = c("available_conditional", "prior_only", "unavailable"))))
    names(tables$scoring_status) <- c("Status", "Persons")
    basis <- rbind(basis, data.frame(Target = "Person ability",
      Interval = paste0(100 * person_scores$settings$level, "% conditional equal-tail posterior intervals"),
      Limitation = "Calibration held fixed; prior_only is prior information; unavailable rows retain their reason."))
  }
  if (!is.null(intervals)) {
    ci <- intervals$intervals
    tables$bootstrap_intervals <- data.frame(Rater = rownames(ci), Lower = ci[, 1], Upper = ci[, 2], row.names = NULL)
    tables$bootstrap_availability <- attr(ci, "availability")
    tables$bootstrap_trials <- intervals$trials
    tables$bootstrap_settings <- mfrm_extended_settings_table(c(intervals$settings,
      list(method = attr(ci, "method"), level_used = attr(ci, "level"),
        expected_tail_draws = attr(ci, "expected_tail_draws"), note = attr(ci, "note"))))
    basis <- rbind(basis, data.frame(Target = "Observed-rater severity (bootstrap)",
      Interval = paste0(100 * attr(ci, "level"), "% pointwise ", attr(ci, "method"), " prediction intervals"),
      Limitation = "Model-based bootstrap; unresolved refits and infinite endpoints are retained; no general coverage guarantee."))
  }
  tables$interval_basis <- basis
  core <- c("fit", "tables", "reporting", "precision", "plots")
  status <- do.call(rbind, lapply(include, function(section) mfrm_results_status_row(section,
    if (section %in% core) "available" else "not_available",
    if (section %in% core) "Stored model-specific results only; see interval_basis and numerical_checks." else
      "This diagnostic or ordinary-model section is not implemented for this model class; nothing was computed.")))
  status <- rbind(status,
    mfrm_results_status_row("predictions", if (is.null(predictions)) "not_computed" else "available",
      if (is.null(predictions)) "No saved predictions supplied; predict() is never called automatically." else "Matching saved predictions retained, including unavailable rows."),
    mfrm_results_status_row("bootstrap_intervals", if (testlet) "not_available" else if (is.null(intervals)) "not_computed" else "available",
      if (testlet) "This interval API is for random-rater fits only." else if (is.null(intervals))
        "No saved random-rater bootstrap supplied; no resampling performed." else "Matching saved intervals and every planned trial retained."))
  if (!testlet) status <- rbind(status, mfrm_results_status_row("person_scores",
    if (is.null(scores)) "not_computed" else "available",
    if (is.null(scores)) "No saved Person scores supplied; no scoring performed." else
      "Saved conditional scores retain the complete roster, all requested Persons and unavailable rows."))
  if (!is.null(comparison)) {
    tables$comparison_models <- comparison$models
    tables$comparison_checks <- comparison$checks
    tables$comparison_effects <- comparison$effects
    tables$comparison_notes <- comparison$notes
    tables$comparison_omitted_events <- comparison$omitted_events
    if (!is.null(comparison$persons)) {
      tables$comparison_person_scores <- comparison$persons$table
      tables$comparison_person_interpretation <- data.frame(Note=comparison$persons$note,Level=comparison$persons$level)
    }
    if (!is.null(comparison$responses)) for (name in c("events", "rows", "probabilities", "measures", "settings")) {
      tables[[paste0("comparison_response_", name)]] <- comparison$responses[[name]]
    }
    status <- rbind(status, mfrm_results_status_row("comparison", "available",
      "Matched-event descriptive facet comparison; no automatic model ranking or difference intervals."))
  }
  if (!is.null(diagnostics)) {
    tables <- c(tables, mfrm_response_diagnostic_tables(diagnostics))
  }
  status <- rbind(status, mfrm_results_status_row("response_diagnostics",
    if (is.null(diagnostics)) "not_computed" else tables$response_overview$Status,
    if (is.null(diagnostics)) "No saved response diagnostics supplied; no integration performed." else
      tables$response_overview$Detail))
  tables$section_status <- status
  matching_scores <- !is.null(person_scores) && isTRUE(tryCatch(mfrm_validate_person_scores(fit,person_scores),error=function(e) FALSE))
  locations <- tryCatch(mfrm_model_locations(fit,if(matching_scores) person_scores else NULL),error=function(e) NULL)
  if (!is.null(locations)) tables$model_locations <- locations
  plot_types <- if (testlet) c("calibration", "scores") else c("raters", "intervals", "scores")
  available <- if (testlet) c(length(fit$input$columns$facets) > 0, !is.null(predictions)) else c(TRUE, !is.null(intervals), !is.null(scores))
  plot_map <- data.frame(Type = plot_types, Available = available & "plots" %in% include,
    Route = paste0('plot(res, type = "', plot_types, '")'),
    RequiredArtifact = if (testlet) c("Fitted fixed facet", "Matching saved Person scores") else c("Fitted raters", "Matching saved bootstrap intervals", "Matching saved Person scores"),
    Detail = if (testlet) c("First fixed facet by default; select others with facet=.", "Conditional intervals; prior-only and unavailable rows retained.") else
      c("Rater point estimates; normal intervals require explicit selection and do not classify quality.", "Pointwise bootstrap intervals; arrows denote infinite endpoints.",
        "Conditional Person intervals with joint rater Laplace integration; unavailable rows retained."))
  plot_map <- rbind(plot_map, data.frame(Type = "comparison", Available = !is.null(comparison) && "plots" %in% include,
    Route = 'plot(res, type = "comparison")', RequiredArtifact = "Matching saved extended-model comparison",
    Detail = "Centered facet effects; paired or mean-versus-difference display, without difference intervals."))
  plot_map <- rbind(plot_map, data.frame(Type = "response_comparison", Available = !is.null(comparison$responses) && "plots" %in% include,
    Route = 'plot(res, type = "response_comparison")', RequiredArtifact = "Matching saved predictive comparison",
    Detail = "Expected scores by default; select variance, probability, infit or outfit with metric=."))
  plot_map <- rbind(plot_map, data.frame(Type = "response_diagnostics", Available = !is.null(diagnostics) && "plots" %in% include,
    Route = 'plot(res, type = "response_diagnostics")', RequiredArtifact = "Matching saved response diagnostics",
    Detail = "Descriptive posterior predictive Infit/Outfit; paired or scatter view; no reference cutoffs."))
  pathway_facets <- if(!is.null(locations) && !is.null(diagnostics)) intersect(stats::na.omit(locations$Facet),diagnostics$settings$group_by) else character()
  plot_map <- rbind(plot_map,data.frame(Type=c("wright","fit_pathway"),
    Available=c(matching_scores && !is.null(locations),length(pathway_facets)>0L) & "plots" %in% include,
    Route=c('plot(res, type = "wright")','plot(res, type = "fit_pathway")'),
    RequiredArtifact=c("Matching source-roster Person scores","Matching posterior diagnostics and located groups"),
    Detail=c("Conditional reference locations; Person intervals only. Older scores need complete scoring_data.",
      "Selected-row descriptive Infit/Outfit versus reference locations; no classic cutoffs.")))
  out <- list(fit = fit, model_family = model, predictions = predictions, intervals = intervals, scores = scores, comparison = comparison, diagnostics = diagnostics,
    calibration_intervals = list(method = calibration_intervals, level = calibration_level),
    include = include, tables = tables, table_index = mfrm_results_table_index(tables), plot_map = plot_map,
    status = status, decision = decision, notes = notes,
    input = list(mode = class(fit)[1], reproducible_code = "# Export with include = c('rds', 'replay') to replay stored results without refitting."))
  class(out) <- "mfrm_results"
  out
}

mfrm_extended_results_summary <- function(object, digits, top_n, view) {
  fit <- object$fit
  actions <- data.frame(Priority = 1L, Area = "Interpretation", Action = object$decision$NextAction,
    Route = "res$tables$numerical_checks; res$tables$interval_basis; res$tables$section_status")
  structure(list(model_family = object$model_family,
    overview = data.frame(Model = object$model_family, Method = fit$settings$method,
      N = nrow(fit$input$data), Persons = length(unique(fit$input$data[[fit$input$columns$person]])),
      Tables = length(object$tables), PlotRoutes = sum(object$plot_map$Available)),
    decision = object$decision, status = object$status,
    component_index = data.frame(Component = c("fit", if (!is.null(object$predictions)) "predictions",
      if (!is.null(object$intervals)) "intervals", if (!is.null(object$scores)) "scores",
      if (!is.null(object$comparison)) "comparison", if (!is.null(object$diagnostics)) "response_diagnostics")),
    table_index = utils::head(object$table_index, top_n), plot_map = object$plot_map,
    readiness = data.frame(Domain = "Model adequacy", Status = "Not assessed"),
    fit_readiness = data.frame(), fit_readiness_components = data.frame(), fit_readiness_parameters = data.frame(),
    triage = data.frame(Area = "Interpretation", Severity = "Review", Detail = object$decision$Interpretation),
    next_actions = actions, mapping = object$tables$column_roles,
    reproducible_code = mfrm_results_code_table(object$input$reproducible_code),
    notes = object$notes, digits = digits, view = view), class = "summary.mfrm_results")
}

mfrm_extended_results_plot <- function(x, type, ...) {
  if (is.null(type)) type <- x$plot_map$Type[which(x$plot_map$Available)[1L]]
  if (length(type) != 1L || is.na(type) || !type %in% x$plot_map$Type ||
      !isTRUE(x$plot_map$Available[match(type, x$plot_map$Type)])) {
    stop("Choose an available type from res$plot_map; supply saved predictions, Person scores or intervals when required.", call. = FALSE)
  }
  if (type == "intervals") {
    args <- list(...)
    if (is.null(args$level)) args$level <- attr(x$intervals$intervals, "level")
    if (is.null(args$method)) args$method <- attr(x$intervals$intervals, "method")
    return(do.call(plot, c(list(x = x$intervals), args)))
  }
  if (type == "comparison") return(plot(x$comparison, ...))
  if (type == "calibration") {
    args <- list(...)
    # Older saved result bundles retain their original automatic 95% display.
    # Rebuilding from the saved fit applies the current default policy.
    selection <- x$calibration_intervals %||% list(method = "normal", level = .95)
    if (is.null(args$intervals)) args$intervals <- selection$method
    if (is.null(args$level)) args$level <- selection$level
    return(do.call(plot, c(list(x$fit), args)))
  }
  if (type %in% c("wright","fit_pathway")) return(mfrm_model_map(x,type,...))
  if (type == "response_diagnostics") return(plot(x$diagnostics, ...))
  if (type == "response_comparison") {
    args <- list(...)
    if (is.null(args$metric)) args$metric <- "expected_score"
    return(do.call(plot, c(list(x$comparison), args)))
  }
  switch(type, raters = plot(x$fit, ...), scores = plot(x$scores %||% x$predictions, ...))
}

mfrm_extended_report <- function(x, style) {
  tables <- x$tables
  numerical <- isTRUE(x$fit$checks$NumericalReady) && isTRUE(x$fit$checks$InformationPositive)
  first_screen <- data.frame(Area = c("Overall", "Numerical checks", "Interval interpretation", "Response residuals"),
    Status = c(if (numerical) "caveat" else "review", if (numerical) "ok" else "review", "caveat", "not_computed"),
    Readiness = c("Target-specific review", "Numerical only", "Conditional on stated assumptions", "Not supplied"),
    MainIssue = c(x$decision$Interpretation, "Numerical convergence does not establish model adequacy.",
      "Read the target and uncertainty basis for each reported quantity.", "No saved response diagnostics were supplied. Ordinary residual diagnostics are not implemented for this model class."),
    NextAction = c(x$decision$NextAction, "Inspect all stored checks and boundary status.",
      "Retain interval limitations, unavailable rows and omitted-score counts.", "If needed, review numerical readiness, then compute mfrm_response_diagnostics(fit) and attach it with mfrm_results(fit, response_diagnostics = diagnostics). These summaries are descriptive; ordinary-model cutoffs do not apply."),
    PrimaryRoute = c("report$tables$interpretation", "report$tables$numerical_checks",
      "report$tables$interval_basis", "report$tables$section_status"))
  if (!is.null(x$diagnostics)) {
    overview <- mfrm_response_diagnostic_overview(x$diagnostics)
    tables$response_overview <- overview
    first_screen$Status[4L] <- switch(overview$Status,
      available = "caveat", not_available = "unavailable", "review")
    first_screen$Readiness[4L] <- "Descriptive only"
    first_screen$MainIssue[4L] <- overview$Detail
    first_screen$NextAction[4L] <- "Review selected-row summaries and unavailable rows; do not apply ordinary-model cutoffs."
    first_screen$PrimaryRoute[4L] <- "report$tables$response_overview"
  }
  report_index <- data.frame(Area = names(tables), PrimaryTable = paste0("report$tables$", names(tables)))
  out <- structure(list(title = paste("mfrmr", x$model_family, "Report"), style = style,
    source_include = x$include, decision = x$decision, first_screen = first_screen,
    report_index = report_index, template_index = data.frame(),
    claim_readiness = data.frame(Claim = "Model adequacy and general interval coverage", Readiness = "Not established"),
    tables = tables, source = x), class = "mfrm_report")
  out$markdown <- paste(c(paste0("# ", out$title),
    "This report collects stored results. Report style does not change estimation or establish model adequacy.",
    mfrm_report_markdown_table(first_screen),
    unlist(lapply(names(tables), function(nm) c(paste0("## ", gsub("_", " ", nm)),
      mfrm_report_markdown_table(tables[[nm]]),
      if (nrow(tables[[nm]]) > 20) "First 20 rows shown; all rows are retained in report$tables and CSV exports.")))), collapse = "\n\n")
  out
}
