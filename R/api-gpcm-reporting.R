# Saved GPCM inference is displayed and exported without refitting or changing
# the ordinary fit diagnostics. Tables always retain the requested target/method.
mfrm_gpcm_product_note <- function() paste(
  "Method qualification is incomplete (experimental two-family GPCM).",
  "Numerical convergence and available bounds do not establish general interval coverage or rater-quality decisions.",
  "Inspect the saved numerical checks and target-specific restrictions; curve/step intervals and ordinary fit/bias tests remain unavailable.")

mfrm_slope_interval_display <- function(tab) {
  lower <- tab$CI_Lower; upper <- tab$CI_Upper
  ifelse(is.infinite(lower) | is.infinite(upper), "Unbounded",
    ifelse(is.finite(lower) & is.finite(upper),
      ifelse(tab$CIEligible %in% TRUE, "Available", "Review"),
      ifelse(is.finite(lower) | is.finite(upper), "Incomplete", "Unavailable")))
}

mfrm_gpcm_inference_tables <- function(x) {
  if (inherits(x, "mfrm_facet_intervals") && (isTRUE(x$settings$two_family) ||
      identical(x$settings$procedure, "mml_native_location_model_v1")))
    return(mfrm_facet_interval_tables(x))
  if (inherits(x, "mfrm_slope_intervals")) {
    tab <- attr(x, "diagnostics")
    tab <- tab[intersect(c("SlopeFacet", "SlopeOwner", "SlopeLevel", "ScaleReference", "Estimate", "SE", "LogSE", "CI_Lower", "CI_Upper",
      "NullValue", "PValue", "CIEligible", "CIUse", "CoverageQualified", "InferenceReview"), names(tab))]
    tab$IntervalStatus <- mfrm_slope_interval_display(tab)
    tab$Target <- attr(x, "target"); tab$Method <- attr(x, "method")
    tab$ConfidenceLevel <- paste0(format(100 * attr(x, "level"), trim = TRUE), "%")
    tab$Adjustment <- if (attr(x, "simultaneous") == "none") "Pointwise" else "Bonferroni"
    tables <- list(intervals = tab)
    for (name in c("checks", "numerical_checks", "score_rank")) {
      if (!is.null(attr(x,name))) tables[[name]] <- attr(x,name)
    }
    profile <- attr(x, "profile")
    if (!is.null(profile)) {
      tables$profile <- profile$profile
      tables$profile_endpoints <- profile$endpoints
      tables$wald <- attr(x, "wald")
      tables$profile_checks <- do.call(rbind, lapply(names(profile$attempts), function(key)
        do.call(rbind, lapply(names(profile$attempts[[key]]), function(start)
          cbind(LogSlope = as.numeric(key), Start = start,
            profile$attempts[[key]][[start]]$check)))))
      rownames(tables$profile_checks) <- NULL
    }
    if (!is.null(attr(x, "availability"))) tables$availability <- attr(x, "availability")
    settings <- attr(x, "settings")
    if (!is.null(settings)) {
      tables$settings <- data.frame(Scale = settings$scale,
        SmallSampleFactor = settings$adjustment_factor %||% NA_real_,
        IndependentClusters = if (is.null(settings$clusters)) NA_integer_ else length(unique(settings$clusters$Cluster)),
        BootstrapSeed = settings$seed %||% NA_integer_)
      if (!is.null(settings$recheck)) tables$settings$Reanalysis <- settings$recheck
      if (!is.null(settings$integration)) tables$settings$Integration <- settings$integration
      if (!is.null(settings$engine)) tables$settings$Engine <- settings$engine
      if (!is.null(settings$diagnostic_checks_scope))
        tables$settings$DiagnosticChecksScope <- settings$diagnostic_checks_scope
      if (!is.null(settings$clusters)) tables$clusters <- settings$clusters
      if (!is.null(settings$contrasts)) tables$contrasts <- data.frame(
        Comparison = rownames(settings$contrasts), settings$contrasts, check.names = FALSE)
      if (!is.null(profile)) {
        tables$settings$Slope <- unname(settings$slope)
        if (length(names(settings$slope))) tables$settings$SlopeOwner <- names(settings$slope)
        tables$settings$QuadratureOrder <- settings$quad_points
        tables$settings$ComparisonOrder <- settings$comparison_quad_points
        tables$settings$MaxIterations <- settings$maxit
        tables$settings$MaxSearchSteps <- settings$max_steps
        tables$settings$InitialLogSlopeStep <- settings$initial_step
        tables$settings$ElapsedSeconds <- profile$elapsed
      }
    }
    return(tables)
  }
  if (inherits(x, "mfrm_curve_intervals")) {
    tab <- x$table
    tab$Target <- x$settings$target
    tab$Method <- if (isTRUE(x$settings$point_only)) "Fitted values only" else x$settings$method
    tab$ConfidenceLevel <- if (isTRUE(x$settings$point_only)) "Not available" else
      paste0(format(100 * x$settings$level, trim = TRUE), "%")
    tab$Adjustment <- if (isTRUE(x$settings$point_only)) "Not applicable" else
      if (x$settings$simultaneous == "none") "Pointwise" else "Bonferroni"
    tables <- list(curves = tab, settings = data.frame(
      Covariance = if (isTRUE(x$settings$point_only)) "Not available" else x$settings$method,
      SmallSampleAdjustment = x$settings$adjust,
      IndependentClusters = if (is.null(x$settings$clusters)) NA_integer_ else length(unique(x$settings$clusters$Cluster))))
    if (!is.null(x$settings$clusters)) tables$clusters <- x$settings$clusters
    if (!is.null(x$contexts)) tables$contexts <- x$contexts
    return(tables)
  }
  if (inherits(x, "mfrm_gpcm_bootstrap")) {
    tables <- if (is.null(x$test)) mfrm_gpcm_inference_tables(confint(x)) else list(test = x$test)
    tables$trials <- x$trials
    if (!is.null(x$checks)) tables$checks <- x$checks
    if (!is.null(x$source_checks)) tables$source_checks <- x$source_checks
    tables$sampling <- data.frame(Sampling = x$settings$sampling, Seed = x$settings$seed,
      Planned = nrow(x$trials), Available = sum(x$trials$Available))
    if (!is.null(x$settings$recheck)) tables$sampling$Reanalysis <- x$settings$recheck
    if (!is.null(x$settings$diagnostic_checks_scope))
      tables$sampling$DiagnosticChecksScope <- x$settings$diagnostic_checks_scope
    return(tables)
  }
  stop("Supply saved GPCM slope intervals, curve intervals or a GPCM bootstrap result.", call. = FALSE)
}

mfrm_gpcm_results_inputs <- function(fit, intervals) {
  if (is.null(intervals)) return(NULL)
  if (!inherits(fit, "mfrm_fit") || !identical(fit$config$model, "GPCM")) {
    stop("Saved GPCM inference requires its native GPCM fit.", call. = FALSE)
  }
  supported <- c("mfrm_slope_intervals", "mfrm_curve_intervals", "mfrm_gpcm_bootstrap", "mfrm_facet_intervals")
  if (inherits(intervals, supported)) intervals <- list(inference = intervals)
  if (!is.list(intervals) || !length(intervals) || is.null(names(intervals)) ||
      anyNA(names(intervals)) || anyDuplicated(tolower(names(intervals))) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", names(intervals))) ||
      !all(vapply(intervals, inherits, logical(1), what = supported))) {
    stop("`intervals` must be a saved GPCM result or a named list of saved results; use distinct letter/number/underscore names.", call. = FALSE)
  }
  names(intervals) <- tolower(names(intervals))
  source <- mfrm_gpcm_inference_source(fit)
  for (x in intervals) {
    native <- inherits(x, "mfrm_facet_intervals") &&
      identical(x$settings$procedure, "mml_native_location_model_v1")
    if (inherits(x, "mfrm_facet_intervals") && !isTRUE(x$settings$two_family) && !native)
      stop("Attach native location intervals with a declared procedure and source identity.", call. = FALSE)
    saved <- if (inherits(x, "mfrm_slope_intervals")) attr(x, "source") else
      if (inherits(x, "mfrm_gpcm_bootstrap")) mfrm_gpcm_inference_source(x$source) else x$source
    expected <- if (native) mfrm_native_location_identity(fit) else source
    if (!identical(expected, saved)) stop(
      "Saved GPCM inference must match the fitted parameters, data, population and integration settings. Recompute older intervals without source metadata.", call. = FALSE)
  }
  intervals
}

mfrm_gpcm_results_score_tables <- function(fit, scores) {
  product <- mfrm_has_product_slopes(fit)
  portable <- inherits(scores, "mfrm_calibration_score")
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      !identical(fit$config$model, "GPCM") ||
      !inherits(scores, c("mfrm_unit_prediction", "mfrm_calibration_score")) ||
      (portable && !product))
    stop("GPCM `scores` must be saved predict_mfrm_units() output; two-family fits also accept a format-6 score_mfrm_calibration() result.", call. = FALSE)
  if (!product) {
    if (!same_signature_component(fit$config$score_map %||% fit$prep$score_map, fit$prep$score_map) ||
        !identical(as.numeric(c(fit$config$rating_min, fit$config$rating_max)),
          as.numeric(c(fit$prep$rating_min, fit$prep$rating_max))))
      stop("fit/scores mismatch: the saved calibration category coding is internally inconsistent.", call. = FALSE)
    mfrm_validate_result_identity(scores, prediction_fit_identity(fit), "scores",
      "predict_mfrm_units()", "mfrm_results()")
    prediction_validate_interval_record(scores)
    if (!identical(scores$settings$method, prediction_resolve_fit_method(fit)) ||
        !isTRUE(scores$settings$score_integration_required) ||
        is.null(scores$settings$prior_comparison) ||
        !is.data.frame(scores$quadrature_review) ||
        !identical(unique(scores$quadrature_review$FixedNodes), scores$settings$scoring_quad_points) ||
        !identical(scores$settings$quad_points, scores$settings$scoring_quad_points) ||
        !identical(sort(unique(scores$quadrature_review$AdaptiveNodes)),
          sort(scores$settings$adaptive_quad_points)) ||
        !identical(scores$settings$score_integration_review,
          prediction_score_integration_review(scores$estimates, scores$quadrature_review)) ||
        !is.null(scores$settings$two_family_calibration))
      stop("Saved scores require matching estimator, prior and scoring-integration records. Recreate them with predict_mfrm_units() from the saved fit.", call. = FALSE)
  }
  if (portable) {
    mfrmr_validate_calibration_score(scores)
    if (!identical(scores$settings$schema_version, 6L))
      stop("Two-family results require portable scoring format 6.", call. = FALSE)
    source <- scores$settings$source_scoring_evidence
    estimates <- scores$estimates
    review <- sum(scores$person_dispositions$Disposition == "scored_review")
    not_scored <- sum(scores$person_dispositions$Disposition == "not_scored")
    notes <- scores$notes
  } else {
    if (product && !identical(scores$settings$source_scoring_policy_basis, "two_family_conditional_eap_v1"))
      stop("Two-family results require scores from the two-family calibration.", call. = FALSE)
    prediction_validate_population_output(scores)
    source <- list(status = scores$settings$source_scoring_status,
      local_calibration_review = scores$settings$local_calibration_review,
      source_audit_states = scores$settings$source_audit_states)
    estimates <- prediction_estimate_table(scores)
    review <- sum(estimates$EstimateUse != "fitted_object_scoring")
    # Native scoring records row omissions; it does not return a not-scored roster.
    not_scored <- NA_integer_
    notes <- prediction_output_notes(scores)
  }
  identity <- source$local_calibration_review$source_identity
  if (product) {
    legacy <- is.character(identity) && length(identity) == 1L &&
      isTRUE(startsWith(identity, "gmfrm_scoring_source_v1:"))
    if (is.null(identity) || !identical(identity, prediction_gmfrm_source_identity(fit, legacy)))
      stop("Saved scores must match this fit's calibration, data, owner/category coding, integration and source status. Re-score older outputs if their source identity is absent or incompatible; the fit need not be re-estimated.", call. = FALSE)
    settings <- scores$settings
    params <- expand_params(fit$opt$par, build_param_sizes(fit$config), fit$config)
    if (portable) {
      # Compare the saved operational coordinates as well as the source fingerprint.
      # A valid edited artifact can still describe a different calibration.
      semantic <- settings$semantic_components
      model <- semantic$facet_order_roles_levels_signs
      facets <- as.character(fit$config$facet_names)
      levels <- do.call(rbind, lapply(facets, function(facet) data.frame(
        Facet = facet, Level = as.character(fit$config$facet_levels[[facet]]),
        LevelIndex = seq_along(fit$config$facet_levels[[facet]]))))
      roles <- data.frame(Facet = facets,
        Role = ifelse(facets == fit$config$step_facet, "facet_and_step_owner", "facet"),
        OrderIndex = seq_along(facets))
      score_map <- as.data.frame(fit$prep$score_map)
      score_map$OrderIndex <- seq_len(nrow(score_map))
      score_map <- score_map[c("OriginalScore", "InternalScore", "OrderIndex")]
      contexts <- unique(as.data.frame(unname(lapply(fit$prep$data[fit$config$slope_facet],
        as.character)), check.names = FALSE))
      names(contexts) <- fit$config$slope_facet
      rownames(contexts) <- NULL
      coordinates <- mfrmr_calibration_extract_coordinates(fit, params)
      matching <- identical(semantic$parameter_coordinates, coordinates) &&
        identical(model, list(facet_names = facets, facet_roles = roles, facet_levels = levels,
          facet_signs = fit$config$facet_signs[facets], step_owner = fit$config$step_facet)) &&
        identical(semantic$slope_specification, list(slope_owner = fit$config$slope_facet,
          slope_action = "full_adjacent_predictor")) &&
        identical(semantic$response_map, list(score_map = score_map,
          rating_min = as.integer(fit$prep$rating_min), rating_max = as.integer(fit$prep$rating_max),
          n_categories = as.integer(fit$config$n_cat))) &&
        identical(semantic$product_structure, list(slope_composition = "product",
          observed_contexts = contexts, source_columns = mfrmr_calibration_normalize_source_columns(fit))) &&
        identical(semantic$identification_constraints, mfrmr_calibration_product_identification(
          fit$config$slope_facet, fit$config$facet_levels, fit$config$n_cat - 1L)) &&
        identical(semantic$typed_anchors, mfrmr_calibration_extract_anchors(fit)) &&
        identical(semantic$interaction_map, mfrmr_calibration_extract_interactions(fit, coordinates))
    } else matching <- identical(settings$two_family_calibration,
      prediction_gmfrm_calibration(fit, params))
    if (!matching)
      stop("Saved scoring calibration records do not match this fit. Attach scores from the matching calibration; editing a source identity does not change the calibration used for scoring.", call. = FALSE)
  } else {
    settings <- scores$settings
    identity <- mfrm_checkpoint_fingerprint(scores$fit_identity, canonical = TRUE)
  }
  tables <- list(
    person_scores = estimates,
    scoring_overview = data.frame(Scored = nrow(estimates), Review = review,
      NotScored = not_scored, Status = if (review > 0 || isTRUE(not_scored > 0)) "review" else "conditional",
      SourceStatus = source$status),
    scoring_settings = data.frame(Route = if (portable) "portable_format_6" else "fitted_object",
      SourceIdentity = identity, CalibrationId = settings$calibration_id %||% NA_character_,
      PriorMean = if (product) 0 else if (length(unique(estimates$PriorMean)) == 1L) estimates$PriorMean[1L] else NA_real_,
      PriorSD = if (product) 1 else if (length(unique(estimates$PriorSD)) == 1L) estimates$PriorSD[1L] else NA_real_,
      ScoringAlgorithm = settings$scoring_algorithm,
      QuadratureOrder = settings$scoring_quad_points %||% settings$quadrature_order,
      IntervalLevel = settings$interval_level,
      UncertaintyBasis = "conditional_on_fixed_calibration"),
    scoring_source_checks = source$local_calibration_review$integration %||% data.frame(),
    scoring_source_audits = data.frame(Domain = names(source$source_audit_states),
      Status = unname(source$source_audit_states)),
    scoring_integration = settings$score_integration_review,
    scoring_row_review = scores$row_review,
    scoring_notes = data.frame(Note = unique(c(notes,
      if (!portable) "Native row omissions are recorded in scoring_row_review; NotScored is not enumerated by this route.")))
  )
  if (!product) {
    tables$scoring_settings$CalibrationMethod <- settings$method
    tables$scoring_settings$PosteriorBasis <- settings$posterior_basis
    tables$scoring_settings$ReadinessPolicy <- settings$readiness_policy
    tables$scoring_priors <- settings$prior_comparison
    tables$scoring_population_review <- scores$population_review
  }
  if (portable) {
    tables$scoring_person_dispositions <- scores$person_dispositions
    tables$scoring_row_dispositions <- scores$row_dispositions
  }
  tables
}

mfrm_gpcm_results_attach_scores <- function(out, scores) {
  if (is.null(scores)) return(out)
  tables <- mfrm_gpcm_results_score_tables(out$fit, scores)
  out$scores <- out$components$scores <- scores
  out$tables <- c(out$tables, tables)
  out$table_index <- mfrm_results_table_index(out$tables)
  out$status <- rbind(out$status, mfrm_results_status_row("person_scores",
    tables$scoring_overview$Status,
    if (mfrm_has_product_slopes(out$fit))
      "Saved two-family conditional EAP; source/batch checks and review labels retained. Calibration uncertainty excluded; no rescoring." else
      "Saved conditional EAP; source/batch checks and review labels retained. Calibration and population estimation uncertainty excluded; no rescoring."))
  out$notes <- unique(c(out$notes, if (mfrm_has_product_slopes(out$fit))
    "Saved new-Person scores are conditional on fixed calibration and the N(0,1) prior. Review-only results remain review-only; coverage and population transport are not established." else
    "Saved Person scores are conditional on fixed calibration and their recorded scoring priors. Review-only results remain review-only; calibration/population estimation uncertainty and frequentist coverage are not supplied by these intervals."))
  out
}

mfrm_gpcm_results_validate_scores <- function(x) {
  if (!identical(x$fit$config$model, "GPCM")) return(invisible(x))
  if (is.null(x$scores) && is.null(x$components$scores) &&
      !any(names(x$tables) == "person_scores" | startsWith(names(x$tables), "scoring_"))) return(invisible(x))
  tables <- mfrm_gpcm_results_score_tables(x$fit, x$scores)
  if (!identical(x$components$scores, x$scores) ||
      !identical(x$tables[names(tables)], tables) ||
      !identical(x$status$Status[x$status$Section == "person_scores"], tables$scoring_overview$Status))
    stop("Saved scoring tables are inconsistent with the attached scores. Rebuild mfrm_results() with the matching saved scores; no rescoring is needed.", call. = FALSE)
  invisible(x)
}

mfrm_gpcm_product_point_tables <- function(fit) {
  list(
    fitted_slopes = data.frame(Facet = fit$slopes$SlopeOwner, Level = fit$slopes$SlopeFacet,
      FittedSlope = fit$slopes$OptimizerEstimate,
      Interval = "See attached slope-interval tables, if requested"),
    fitted_locations = data.frame(Facet = fit$facets$others$Facet,
      Level = fit$facets$others$Level, FittedLocation = fit$facets$others$Estimate,
      Constraint = ifelse(fit$facets$others$Facet == fit$config$slope_facet[1],
        "Sum to zero across levels", "Uncentered; fixed N(0,1) ability"),
      Interval = "See attached location-interval tables, if requested"),
    fitted_steps = data.frame(StepOwner = fit$config$step_facet,
      StepLevel = fit$steps$StepFacet, Step = fit$steps$Step,
      FittedStep = fit$steps$Estimate, Interval = "Not available"))
}

mfrm_gpcm_results_attach <- function(out, inputs) {
  if (mfrm_has_product_slopes(out$fit)) {
    out$tables$gpcm_interpretation_notes <- data.frame(Note = mfrm_gpcm_product_note())
    out$notes <- unique(c(out$notes, mfrm_gpcm_product_note()))
    if ("fit" %in% out$include) {
      points <- mfrm_gpcm_product_point_tables(out$fit)
      names(points) <- paste0("gpcm_", names(points))
      out$tables <- c(out$tables, points)
    }
    if (!is.null(out$fit$opt$mml_initialization))
      out$tables$gpcm_initialization <- out$fit$opt$mml_initialization$table
    out$table_index <- mfrm_results_table_index(out$tables)
  }
  if (is.null(inputs)) return(out)
  out$gpcm_inference <- inputs
  for (name in names(inputs)) {
    point_only <- inherits(inputs[[name]], "mfrm_curve_intervals") && isTRUE(inputs[[name]]$settings$point_only)
    key <- paste0("gpcm_", name)
    tables <- mfrm_gpcm_inference_tables(inputs[[name]])
    names(tables) <- paste(key, names(tables), sep = "_")
    out$tables <- c(out$tables, tables)
    out$components[[key]] <- inputs[[name]]
    out$status <- rbind(out$status, mfrm_results_status_row(key, "review",
      if (point_only) "Saved provisional fitted curves; calibration intervals are unavailable." else
        "Saved GPCM inference; target, approximation, multiplicity and unavailable outcomes retained. No automatic rater-quality decision."))
    out$plot_map <- dplyr::bind_rows(out$plot_map, data.frame(Type = key,
      Available = "plots" %in% out$include, RequiredArtifact = FALSE,
      Route = paste0('plot(res, type = "', key, '")'),
      Detail = if (point_only) "Saved provisional fitted curves without intervals; no recalculation." else
        "Saved GPCM uncertainty; no recalculation or replacement of ordinary fit diagnostics.",
      InterpretationStatus = if (point_only) "provisional_fitted_values" else "approximate_inference", InterpretationReady = FALSE,
      ReadinessRoute = paste0("res$tables$", names(tables)[1])))
  }
  out$table_index <- mfrm_results_table_index(out$tables)
  out$notes <- unique(c(out$notes,
    unlist(lapply(Filter(function(x) inherits(x, "mfrm_facet_intervals"), inputs),
      `[[`, "cautions"), use.names = FALSE),
    if (mfrm_has_product_slopes(out$fit)) "Two-family curves remain provisional without intervals; separately requested slope or location intervals use experimental approximations with saved numerical checks. No automatic rater-quality decision is supported." else
    "GPCM slope/curve uncertainty is separate from Wright/Pathway location and fit displays. Approximate intervals do not classify raters or choose scoring weights."))
  out
}

# Keep the ordinary fit-report templates from recommending unsupported outputs.
mfrm_gpcm_product_report <- function(x, style) {
  mfrm_gpcm_results_validate_scores(x)
  points <- mfrm_gpcm_product_point_tables(x$fit)
  adaptive <- mfrmr_adaptive_integration(x$fit$config)
  numerical <- isTRUE(x$fit$summary$Converged) &&
    identical(x$fit$opt$optimizer_diagnostics$ConvergenceSeverity, "pass")
  intervals <- Filter(function(z) inherits(z,"mfrm_slope_intervals"),x$gpcm_inference %||% list())
  n_intervals <- sum(vapply(intervals,function(z) sum(attr(z,"diagnostics")$CIEligible),integer(1)))
  interval_note <- if (length(intervals)) paste(n_intervals,
    "experimental slope intervals are available; inspect the saved checks and unavailable rows. Curve intervals remain unavailable.") else
    "No slope intervals were attached. Curve intervals are unavailable. Use confint(fit, parm = 'slopes') for the separately checked experimental slope approximation."
  locations <- Filter(function(z) inherits(z, "mfrm_facet_intervals"), x$gpcm_inference %||% list())
  n_locations <- sum(vapply(locations, function(z) sum(z$table$CIEligible), integer(1)))
  interval_note <- paste(interval_note, if (length(locations)) paste(n_locations,
    "experimental location/contrast intervals are available; inspect their fixed-scale references and saved checks.") else
    "No location intervals were attached; use mfrm_facet_intervals() for the separately checked experimental approximation.")
  first_screen <- data.frame(
    Area = c("Overall", "Numerical fit", "Calibration intervals", "Other outputs"),
    Status = c("review", if (numerical) "ok" else "review", if (n_intervals + n_locations) "caveat" else "unavailable", "unavailable"),
    Readiness = c("Provisional analysis", "Numerical only", "Output-specific checks", "Not supported"),
    MainIssue = c(mfrm_gpcm_product_note(),
      if (adaptive) {
        if (numerical) "Direct adaptive MML met its optimizer and gradient checks." else
          "Direct adaptive MML did not pass its numerical convergence checks."
      } else if (numerical) "The per-Person marginal-score stopping rule was met." else
        "EM stopped without meeting the per-Person marginal-score tolerance.",
      interval_note,
      "Ordinary diagnostics, Wright/Pathway maps and model ranking are unavailable."),
    NextAction = c("Read the model settings and numerical status before the saved curves.",
      if (numerical) "Convergence does not establish identification, model adequacy or inferential reliability." else
        "Review the stopping result and retain this fit as numerically unresolved.",
      "Retain unavailable bounds and labels for contexts not observed together.",
      "Conditional new-Person EAP is separately available through predict_mfrm_units() or the portable calibration workflow with source/batch checks. Attach matching saved output through scores in mfrm_results()."),
    PrimaryRoute = c("report$tables$fit_summary_settings_overview", "report$tables$fit_summary_readiness",
      "report$tables", "gpcm_capability_matrix()"))
  selected <- c("fit_summary_overview", "fit_summary_settings_overview", "fit_summary_readiness",
    "fit_summary_decision", "fit_summary_slope_overview", "fit_summary_caveats")
  tables <- c(x$tables[intersect(selected, names(x$tables))],
    list(slope_estimates = x$fit$slopes),
    x$tables[(startsWith(names(x$tables), "gpcm_") | startsWith(names(x$tables), "response_") |
      startsWith(names(x$tables), "scoring_") | names(x$tables) == "person_scores") &
      !names(x$tables) %in% paste0("gpcm_", names(points))])
  if (!is.null(x$scores)) {
    first_screen <- rbind(first_screen, data.frame(Area = "New-Person scoring",
      Status = x$tables$scoring_overview$Status, Readiness = "Conditional on fixed calibration",
      MainIssue = "Saved EAP, posterior SD and intervals retain their source/batch checks and any review-only or not-scored outcomes.",
      NextAction = "Review scoring_source_checks, scoring_integration and person_scores. Intervals exclude calibration uncertainty; coverage and population transport are not established.",
      PrimaryRoute = "report$tables$person_scores"))
  }
  if (!is.null(x$response_diagnostics)) {
    overview <- mfrm_response_diagnostic_overview(x$response_diagnostics)
    tables$response_overview <- overview
    first_screen <- rbind(first_screen, data.frame(Area = "Response residuals",
      Status = switch(overview$Status, available = "caveat", not_available = "unavailable", "review"),
      Readiness = "Descriptive only",
      MainIssue = overview$Detail,
      NextAction = "Review unavailable rows and integration differences before inspecting grouped summaries.",
      PrimaryRoute = "report$tables$response_overview"))
  }
  fit <- x$fit
  settings <- data.frame(Setting = c("Model", "Estimation", "Ability population",
    "First slope facet (geometric mean one)", "Second slope facet (free slopes)",
    "Step facet", "Quadrature points", if (adaptive) "Optimizer stages (selected start)" else "EM iterations", "Numerically converged",
    if (adaptive) "Maximum total gradient" else "Maximum per-Person marginal score", "Stopping tolerance"),
    Value = c("Two-family GPCM", if (adaptive) "MML with direct adaptive quadrature" else "MML with generalized EM", "Fixed standard normal",
      fit$config$slope_facet, fit$config$step_facet,
      fit$config$estimation_control$quad_points,
      if (adaptive) nrow(fit$opt$optimizer_polish$Stages) %||% NA_integer_ else fit$summary$EMIterations,
      fit$summary$Converged, fit$summary$TerminalGradientSupNorm,
      fit$summary$GradientReviewTolerance))
  if (adaptive) settings <- rbind(settings, data.frame(
    Setting = c("Initialization", "Selected start", "Initialization and optimization seconds"),
    Value = c(fit$config$estimation_control$gpcm_mml_start %||% "neutral (legacy)",
      fit$opt$mml_initialization$selected %||% "neutral",
      fit$opt$mml_initialization$total_seconds %||% NA_real_)))
  tables <- c(list(model_settings = settings), points, tables)
  tables$gpcm_interpretation_notes <- data.frame(Note = mfrm_gpcm_product_note())
  out <- structure(list(title = "mfrmr Two-family GPCM Report", style = style,
    source_include = x$include, decision = x$tables$fit_summary_decision,
    fit_readiness = x$fit_readiness, fit_readiness_components = x$fit_readiness_components,
    first_screen = first_screen,
    report_index = data.frame(Area = names(tables), PrimaryTable = paste0("report$tables$", names(tables))),
    template_index = data.frame(),
    claim_readiness = data.frame(Claim = "Formal inference, model adequacy and general coverage", Readiness = "Not established"),
    tables = tables, source = x), class = "mfrm_report")
  out$markdown <- paste(c(paste0("# ", out$title),
    "This report retains provisional numerical results. Report style does not qualify estimates or supply unavailable analyses.",
    mfrm_report_markdown_table(first_screen),
    "## Provisional GPCM analysis",
    interval_note,
    "Calibration intervals are unavailable for the response curves. Experimental slope and location intervals do not establish general coverage. Numerical convergence and curve shape do not qualify rater-quality decisions.",
    "The two fitted slopes multiply each other and the complete adjacent-category predictor. A large slope is not evidence of assessor competence or scoring accuracy.",
    "Fitted locations and category steps are provisional numerical values on the fixed standard-normal ability scale. Separately requested location intervals retain their own checks; step intervals remain unavailable. Location differences do not imply uniform differences in expected ratings when slopes or steps differ. Step values are offsets centered to sum to zero within each step-owner level, not standalone category thresholds. An adjacent-category threshold also includes both facet locations.",
    unlist(lapply(c("model_settings", names(points), intersect("gpcm_initialization", names(tables)),
      names(tables)[grepl("^gpcm_.*_(curves|contexts|intervals|contrasts|profile_endpoints)$", names(tables))],
      intersect(c("response_overview", "response_diagnostic_settings", "response_measures", "scoring_overview",
        "scoring_settings", "scoring_source_checks", "scoring_source_audits", "scoring_integration",
        "person_scores", "scoring_row_review", "scoring_person_dispositions", "scoring_notes"), names(tables))),
      function(name) {
        tab <- tables[[name]]
        if (endsWith(name, "profile_endpoints")) tab$Status <- NULL
        c(paste0("## ", gsub("_", " ", name)), mfrm_report_markdown_table(tab),
          if (nrow(tab) > 20) "First 20 rows shown; all rows remain in report$tables and CSV exports.")
      })),
    "Complete numerical settings, source checks and unavailable primary estimates remain in report$tables and the saved source result."), collapse = "\n\n")
  out
}

#' Display saved GPCM intervals and bootstrap results
#'
#' Plot the selected inference result without fitting or recomputing covariance.
#' These plots describe discrimination or specified contrasts, not rater severity,
#' agreement, quality or recommended scoring weights.
#' @param x Saved [confint.mfrm_fit()] output from `parm = "slopes"`, or
#'   [bootstrap_mfrm_gpcm()] output.
#' @param title,subtitle Optional returned metadata; NULL removes it. The interval
#'   description defaults to the saved target, method, confidence level and any
#'   saved numerical or bootstrap cautions. Custom text overrides this default, while the
#'   diagnostics remain in the plot data and interval tables. Titles and notes
#'   are not drawn. Use `plot_data(p, "notes")` for descriptions and cautions,
#'   and `plot_data(p, "title")` for the suggested figure title.
#' @param reference Optional finite numeric vertical reference. Default NULL;
#'   no quality threshold or scoring-policy cutoff is selected.
#' @param draw TRUE draws the ggplot; FALSE returns it only.
#' @param ... For bootstrap slope results, passed to `confint()` to select its
#'   target and level. Otherwise unused.
#' @details Crosses retain unavailable intervals. Open circles retain intervals
#'   with infinite endpoints; the plot does not replace them with finite bounds.
#'   Incomplete intervals retain the computed endpoint as a tick, without a
#'   connecting interval line. A missing endpoint is not an infinite endpoint
#'   or a calibrated one-sided interval. Finite review-only bounds use a dashed
#'   line. The table retains `IntervalStatus`, `CIEligible`, `CIUse` (when saved),
#'   all bounds and failure reasons; returned notes explain the display.
#'   For a bootstrap LRT, a histogram shows resolved null statistics and a line
#'   the observed statistic. The returned notes count every unresolved replicate.
#'
#'   Use [apa_table()] for target-aware tables, [plot_data()] for plotted values,
#'   and [as_ggplot()] for customization. Attach one result or a named list with
#'   `mfrm_results(fit, intervals = list(slopes = ci, curves = curves),
#'   compute = "never")`. Named plot routes become `gpcm_slopes`, `gpcm_curves`,
#'   etc. Exact source matching is required; default fit diagnostics are preserved.
#'   [export_mfrm_results()] saves these tables and plots; replay reloads the
#'   saved RDS, preserving the selected methods without rerunning estimation.
#'   Markdown report summaries include the finite, unbounded and unavailable
#'   interval counts and the saved inference cautions/reasons.
#'
#'   For positive slope or ratio intervals spanning many orders of magnitude,
#'   use `as_ggplot(ci) + ggplot2::scale_x_log10()`. Finite bounds and estimates
#'   must be strictly positive. Keep a linear axis for signed differences.
#'   Changing the axis does not narrow the interval or make it more reliable.
#' @param type For saved slope intervals, `"interval"` displays bounds;
#'   `"profile"` displays saved likelihood evaluations and the same-target Wald
#'   limits. The default is `"profile"` for saved profile calculations and
#'   `"interval"` otherwise. This never initiates profiling.
#' @param caption For saved slope intervals, optional returned caption. NULL removes it.
#'   Two-family intervals state the distinct family references by default.
#'   Profile plots supply a
#'   default explanation of line styles and missing coordinates when omitted.
#' @return Invisibly a ggplot, retaining an `mfrm_plot_data` attribute for
#'   [plot_data()]. No fitting or resampling is performed.
#' @examples
#' # ci <- confint(fit, parm = "slopes", scale = "standardized", method = "sandwich")
#' # plot(ci, title = NULL)
#' # apa_table(ci)
#' # plot_data(ci, "table")
#' @export
plot.mfrm_slope_intervals <- function(x, title = "GPCM slope uncertainty",
    subtitle = paste(attr(x, "target"), attr(x, "method"),
      paste0(100*attr(x, "level"), "%"),
      if (attr(x, "simultaneous") == "none") "pointwise" else "Bonferroni", sep = " | "),
    reference = NULL, draw = TRUE, ..., type = c("interval", "profile"), caption = NULL) {
  rlang::check_dots_empty(); .require_mfrmr_ggplot2()
  if (missing(type) && !is.null(attr(x, "profile"))) type <- "profile"
  type <- match.arg(type)
  if (!is.null(reference) && (!is.numeric(reference) || length(reference) != 1L ||
      is.complex(reference) || !is.finite(reference))) stop("`reference` must be NULL or one finite number.", call. = FALSE)
  tab <- mfrm_gpcm_inference_tables(x)$intervals
  tab$Display <- tab$IntervalStatus
  tab$Level <- factor(tab$SlopeFacet, levels = rev(tab$SlopeFacet))
  finite <- tab[tab$Display %in% c("Available", "Review"), , drop = FALSE]
  partial <- tab[tab$Display == "Incomplete", , drop = FALSE]
  partial$Bound <- ifelse(is.finite(partial$CI_Lower), partial$CI_Lower, partial$CI_Upper)
  interval_note <- paste(
    "Incomplete: only the computed endpoint is shown as a tick; the other endpoint remains unknown.",
    "This is neither an infinite endpoint nor a calibrated one-sided interval.",
    "Review: finite bounds retained with a dashed line; consult CIEligible, CIUse and InferenceReview.")
  if (missing(subtitle) && length(attr(x, "cautions"))) {
    subtitle <- paste(subtitle, paste(attr(x, "cautions"), collapse = " "), sep = "\n")
  }
  if (!is.null(subtitle)) subtitle <- paste(strwrap(subtitle, width = 75), collapse = "\n")
  if (type != "profile" && missing(caption) && !is.null(attr(x,"scale_note"))) {
    caption <- paste(strwrap(attr(x,"scale_note"),width=80),collapse="\n")
  }
  if (type == "profile") {
    saved <- attr(x, "profile")
    if (is.null(saved)) stop("type = 'profile' requires saved method = 'profile' intervals.", call. = FALSE)
    if (!is.null(reference) && reference <= 0) stop("A profile reference must be positive on the log axis.", call. = FALSE)
    tab <- saved$profile
    tab$Display <- ifelse(tab$Passed, "Passed", "Unresolved")
    tab$Segment <- cumsum(!tab$Passed)
    visible <- is.finite(tab$Slope) & tab$Slope > 0 & is.finite(tab$LR)
    if (missing(caption)) caption <- paste(
      "Dashed: likelihood-ratio cutoff. Dotted: Wald limits. Crosses: unresolved points.",
      sum(!visible), "evaluations without a plotting coordinate remain in the saved tables.",
      attr(x, "scale_note") %||% "")
    if (!is.null(caption)) caption <- paste(strwrap(caption, width = 80), collapse = "\n")
    wald <- attr(x, "wald")
    payload <- new_mfrm_plot_data("gpcm_slope_profile",
      list(table = tab, endpoints = saved$endpoints, wald = wald, title = title,
        subtitle = subtitle, caption = caption, cautions = attr(x, "cautions"), reference = reference,
        target = attr(x, "target"), slope_metadata = attr(x, "diagnostics"),
        x_label = if ("SlopeOwner" %in% names(attr(x, "diagnostics"))) "Component discrimination (log scale)" else "Relative discrimination (log scale)",
        level = attr(x, "level"), cutoff = saved$cutoff))
    p <- mfrm_gg_gpcm_profile(payload)
    if (isTRUE(draw)) print(p)
    return(invisible(p))
  }
  p <- ggplot2::ggplot(tab, ggplot2::aes(y = .data$Level, x = .data$Estimate)) +
    ggplot2::geom_segment(data = finite, ggplot2::aes(x = .data$CI_Lower,
      xend = .data$CI_Upper, yend = .data$Level), color = "#0072B2", linewidth = .8,
      linetype = ifelse(finite$Display == "Review", "dashed", "solid")) +
    ggplot2::geom_point(data = partial, ggplot2::aes(x = .data$Bound),
      shape = 124, size = 4, color = "#0072B2") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$Display), size = 2.8, color = "#0072B2") +
    ggplot2::scale_shape_manual(values = c(Available = 16, Incomplete = 4, Review = 2, Unavailable = 4, Unbounded = 1)) +
    ggplot2::labs(x = attr(x, "target"), y = NULL,
      shape = "Interval", alt = paste("Saved GPCM interval plot;", sum(tab$Display == "Available"),
        "finite intervals,", sum(tab$Display == "Unbounded"), "unbounded,",
        sum(tab$Display == "Incomplete"), "incomplete,", sum(tab$Display == "Review"), "review,",
        sum(tab$Display == "Unavailable"), "unavailable. Ticks retain computed partial bounds.")) +
    .mfrmr_gg_theme() + ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9))
  if (!is.null(reference)) p <- p + ggplot2::geom_vline(xintercept = reference, linetype = "dashed", color = "grey40")
  attr(p, "mfrmr_plot_data") <- new_mfrm_plot_data("gpcm_slope_intervals",
    list(table = tab, title = title, subtitle = subtitle, caption = caption,
      cautions = c(attr(x, "cautions"), interval_note), target = attr(x, "target"),
      uncertainty_note = attr(x, "anchor_note"),
      method = attr(x, "method"), level = attr(x, "level"), adjustment = attr(x, "simultaneous"),
      checks = attr(x, "checks"), numerical_checks = attr(x, "numerical_checks"),
      score_rank = attr(x, "score_rank")))
  p <- .mfrmr_gg_labs(p, attr(p, "mfrmr_plot_data")$data, x = attr(x, "target"))
  if (isTRUE(draw)) print(p)
  invisible(p)
}

#' @rdname plot.mfrm_slope_intervals
#' @export
plot.mfrm_gpcm_bootstrap <- function(x, title = "GPCM bootstrap inference",
    subtitle = "Fitted-model parametric bootstrap; unresolved refits retained", draw = TRUE, ...) {
  if (is.null(x$test)) {
    ci <- confint(x, ...)
    if (missing(subtitle)) return(plot(ci, title = title, draw = draw))
    return(plot(ci, title = title, subtitle = subtitle, draw = draw))
  }
  rlang::check_dots_empty(); .require_mfrmr_ggplot2()
  tab <- x$trials
  caption <- paste(sum(tab$Available), "resolved of", nrow(tab),
    "planned replicates; unresolved:", sum(!tab$Available))
  p <- ggplot2::ggplot(tab[tab$Available, , drop = FALSE], ggplot2::aes(x = .data$LR)) +
    ggplot2::geom_histogram(bins = 25, fill = "#0072B2", color = "white") +
    ggplot2::geom_vline(xintercept = x$test$LR, linetype = "dashed", linewidth = .8) +
    ggplot2::labs(x = "Null likelihood-ratio statistic", y = "Replicates",
      alt = "Histogram of resolved PCM-null bootstrap likelihood-ratio statistics with the observed statistic marked. Unresolved replicate counts are retained.") +
    .mfrmr_gg_theme()
  attr(p, "mfrmr_plot_data") <- new_mfrm_plot_data("gpcm_bootstrap_test",
    list(table = tab, test = x$test, title = title, subtitle = subtitle, caption = caption))
  p <- .mfrmr_gg_labs(p, attr(p, "mfrmr_plot_data")$data,
    x = "Null likelihood-ratio statistic", y = "Replicates")
  if (isTRUE(draw)) print(p)
  invisible(p)
}

#' @rdname as_ggplot
#' @export
as_ggplot.mfrm_slope_intervals <- function(x, type = NULL, component = NULL, ...) {
  if (inherits(x, "mfrm_slope_intervals") && length(type) == 1L && type %in% c("interval", "profile") && is.null(component))
    return(plot(x, type = type, draw = FALSE, ...))
  if (!is.null(type) || !is.null(component)) stop("Use plot_data() to select a GPCM inference table.", call. = FALSE)
  plot(x, draw = FALSE, ...)
}

#' @rdname as_ggplot
#' @name as_ggplot.mfrm_curve_intervals
#' @aliases as_ggplot.mfrm_curve_intervals
#' @usage \method{as_ggplot}{mfrm_curve_intervals}(x, type = NULL, component = NULL, ...)
#' @rawNamespace S3method(as_ggplot,mfrm_curve_intervals,as_ggplot.mfrm_slope_intervals)
NULL

#' @rdname as_ggplot
#' @name as_ggplot.mfrm_gpcm_bootstrap
#' @aliases as_ggplot.mfrm_gpcm_bootstrap
#' @usage \method{as_ggplot}{mfrm_gpcm_bootstrap}(x, type = NULL, component = NULL, ...)
#' @rawNamespace S3method(as_ggplot,mfrm_gpcm_bootstrap,as_ggplot.mfrm_slope_intervals)
NULL

#' @rdname as_ggplot
#' @export
as_ggplot.mfrm_results <- function(x, type = NULL, component = NULL, ...) {
  if ((mfrm_has_product_slopes(x$fit) && is.null(type)) ||
      (!is.null(type) && length(type) == 1L && !is.na(type) &&
      tolower(type) %in% paste0("gpcm_", names(x$gpcm_inference)))) {
    if (!is.null(component)) stop("Use plot_data() to select a GPCM inference table.", call. = FALSE)
    p <- plot(x, type = type, draw = FALSE, ...)
    return(if (inherits(p, "ggplot")) p else as_ggplot(p))
  }
  NextMethod()
}


mfrm_gg_gpcm_profile <- function(x) {
  .require_mfrmr_ggplot2()
  d <- x$data; tab <- d$table
  visible <- is.finite(tab$Slope) & tab$Slope > 0 & is.finite(tab$LR)
  p <- ggplot2::ggplot(tab[visible, ], ggplot2::aes(x = .data$Slope, y = .data$LR)) +
    ggplot2::geom_line(data = tab[visible & tab$Passed, ],
      ggplot2::aes(group = .data$Segment), color = "#0072B2") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$Display), size = 2, color = "#0072B2") +
    ggplot2::scale_shape_manual(values = c(Passed = 16, Unresolved = 4)) +
    ggplot2::geom_hline(yintercept = d$cutoff, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = unlist(d$wald[c("Lower", "Upper")]), linetype = "dotted") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(shape = "Numerical check",
      alt = "Saved local slope profile. Dashed horizontal cutoff and dotted Wald limits; unresolved points use crosses. Gaps are not interpolated.") +
    .mfrmr_gg_theme() + ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9))
  if (!is.null(d$reference)) p <- p + ggplot2::geom_vline(xintercept = d$reference, linetype = "dotdash")
  p <- .mfrmr_gg_labs(p, d, x = d$x_label %||% "Relative discrimination (log scale)",
    y = "Twice the log-likelihood loss")
  attr(p, "mfrmr_plot_data") <- x
  p
}
