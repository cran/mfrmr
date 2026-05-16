# mfrmr release recovery-validation protocol
#
# Source this file in a development or release-check session:
#
#   source(system.file("validation", "recovery-validation.R", package = "mfrmr"))
#   result <- mfrmr_run_recovery_validation(quick = TRUE)
#   summary <- mfrmr_summarize_recovery_validation(result)
#
# The functions are intentionally not exported. They provide a reproducible
# long-run validation protocol without adding heavy Monte Carlo work to CRAN
# or routine package tests.

mfrmr_recovery_validation_prompt_steps <- function() {
  data.frame(
    Step = seq_len(7L),
    Label = c(
      "Aim",
      "Case selection",
      "Generator contract",
      "Estimator contract",
      "Recovery metrics",
      "Decision thresholds",
      "Report handoff"
    ),
    Prompt = c(
      "What mathematical claim is being checked: equal-slope locations, step profiles, slope profiles, or population parameters?",
      "Which validation tier is needed now: quick smoke, core release evidence, or extended sensitivity evidence?",
      "Are the true locations, thresholds, slopes, assignment design, and population model explicit and identified?",
      "Does the fitted model use the same model family, step facet, slope facet, quadrature tier, and person-population contract as the generator?",
      "Which row-level recovery quantities are estimands, and which are only diagnostics or uncertainty checks?",
      "What RMSE, bias, convergence, coverage, and Monte Carlo precision limits define ok, review, or concern for this use case?",
      "Which tables should be kept for release review: review steps, case summary, metric summary, decision table, and run notes?"
    ),
    Output = c(
      "validation aim",
      "case plan",
      "identified simulation specification",
      "fit settings",
      "recovery and assessment objects",
      "status-coded metric review",
      "CSV/RDS/Markdown validation bundle"
    ),
    stringsAsFactors = FALSE
  )
}

mfrmr_recovery_validation_case_plan <- function() {
  data.frame(
    CaseID = c(
      "rsm_equal_slope_location",
      "pcm_step_profiles",
      "gpcm_slope_profile",
      "rsm_latent_regression"
    ),
    Tier = c("core", "core", "core", "extended"),
    Model = c("RSM", "PCM", "GPCM", "RSM"),
    Estimation = c("JML", "JML", "MML", "MML"),
    StepFacet = c(NA, "Criterion", "Criterion", NA),
    SlopeFacet = c(NA, NA, "Criterion", NA),
    RecommendedReps = c(50L, 50L, 40L, 30L),
    SmokeReps = c(2L, 2L, 1L, 1L),
    Maxit = c(60L, 70L, 90L, 80L),
    QuadPoints = c(NA_integer_, NA_integer_, 11L, 11L),
    IncludePerson = c(TRUE, TRUE, FALSE, FALSE),
    Purpose = c(
      "Equal-slope location and common-step recovery under the baseline ordered MFRM route.",
      "Step-facet-specific threshold recovery under the partial-credit route.",
      "Bounded GPCM relative-discrimination and threshold recovery on the identified log-slope scale.",
      "Latent-regression coefficient and variance recovery under the MML population branch."
    ),
    PrimaryRisk = c(
      "Location indeterminacy or step-centering mistakes.",
      "Mixing common-step and step-facet threshold contracts.",
      "Comparing unnormalized generator slopes with identified fitted slopes.",
      "Treating population recovery as a design-forecasting claim."
    ),
    SummaryFocus = c(
      "facet, person, and step RMSE/Bias after location alignment",
      "step-profile RMSE/Bias and facet recovery",
      "identified log-slope RMSE/Bias plus GPCM step recovery",
      "population coefficient and variance recovery"
    ),
    stringsAsFactors = FALSE
  )
}

mfrmr_recovery_validation_threshold_table <- function(levels,
                                                      score_levels = 4L,
                                                      base = NULL,
                                                      offsets = NULL) {
  if (is.null(base)) {
    base <- seq(-1.2, 1.2, length.out = score_levels - 1L)
  }
  if (is.null(offsets)) {
    offsets <- seq(-0.25, 0.25, length.out = length(levels))
  }
  if (length(offsets) != length(levels)) {
    stop("`offsets` must have one value per step-facet level.", call. = FALSE)
  }
  rows <- lapply(seq_along(levels), function(i) {
    data.frame(
      StepFacet = levels[i],
      StepIndex = seq_along(base),
      Estimate = as.numeric(base + offsets[i]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mfrmr_recovery_validation_population_covariates <- function(n = 80L) {
  data.frame(
    TemplatePerson = sprintf("T%03d", seq_len(n)),
    X = seq(-1.5, 1.5, length.out = n),
    G = rep(c("A", "B"), length.out = n),
    stringsAsFactors = FALSE
  )
}

mfrmr_recovery_validation_fun <- function(name) {
  ns <- asNamespace("mfrmr")
  if (exists(name, envir = ns, inherits = FALSE)) {
    return(get(name, envir = ns, inherits = FALSE))
  }
  getExportedValue("mfrmr", name)
}

mfrmr_recovery_validation_spec <- function(case_id) {
  if (!requireNamespace("mfrmr", quietly = TRUE)) {
    stop("The `mfrmr` package must be available before running validation.", call. = FALSE)
  }
  case_id <- as.character(case_id[1])
  build_spec <- mfrmr_recovery_validation_fun("build_mfrm_sim_spec")

  switch(
    case_id,
    rsm_equal_slope_location = build_spec(
      n_person = 60,
      n_rater = 4,
      n_criterion = 4,
      raters_per_person = 2,
      score_levels = 4,
      model = "RSM",
      assignment = "rotating",
      theta_sd = 1.0,
      rater_sd = 0.35,
      criterion_sd = 0.25,
      step_span = 1.2
    ),
    pcm_step_profiles = build_spec(
      n_person = 60,
      n_rater = 4,
      n_criterion = 4,
      raters_per_person = 2,
      score_levels = 4,
      model = "PCM",
      step_facet = "Criterion",
      assignment = "rotating",
      thresholds = mfrmr_recovery_validation_threshold_table(
        levels = sprintf("C%02d", 1:4),
        offsets = c(-0.30, -0.10, 0.10, 0.30)
      )
    ),
    gpcm_slope_profile = build_spec(
      n_person = 72,
      n_rater = 3,
      n_criterion = 3,
      raters_per_person = 3,
      score_levels = 4,
      model = "GPCM",
      step_facet = "Criterion",
      slope_facet = "Criterion",
      assignment = "crossed",
      thresholds = mfrmr_recovery_validation_threshold_table(
        levels = sprintf("C%02d", 1:3),
        offsets = c(-0.25, 0.00, 0.25)
      ),
      slopes = c(C01 = 0.70, C02 = 1.00, C03 = 1.45)
    ),
    rsm_latent_regression = build_spec(
      n_person = 70,
      n_rater = 4,
      n_criterion = 4,
      raters_per_person = 2,
      score_levels = 4,
      model = "RSM",
      assignment = "rotating",
      population_formula = ~ X + G,
      population_coefficients = c("(Intercept)" = 0.10, X = 0.45, GB = 0.35),
      population_sigma2 = 0.55,
      population_covariates = mfrmr_recovery_validation_population_covariates()
    ),
    stop("Unknown recovery-validation case: ", case_id, call. = FALSE)
  )
}

mfrmr_recovery_validation_case_thresholds <- function(case_id) {
  switch(
    as.character(case_id[1]),
    rsm_equal_slope_location = list(
      max_rmse = c(person = 1.25, facet = 0.75, step = 0.95, default = 1.25),
      max_abs_bias = c(default = 0.55),
      min_se_available = NULL
    ),
    pcm_step_profiles = list(
      max_rmse = c(person = 1.35, facet = 0.85, step = 1.00, default = 1.35),
      max_abs_bias = c(default = 0.65),
      min_se_available = NULL
    ),
    gpcm_slope_profile = list(
      max_rmse = c(slope = 0.40, facet = 0.90, step = 1.05, default = 1.35),
      max_abs_bias = c(slope = 0.25, default = 0.70),
      min_se_available = NULL
    ),
    rsm_latent_regression = list(
      max_rmse = c(population = 0.65, facet = 0.90, step = 1.10, default = 1.35),
      max_abs_bias = c(population = 0.45, default = 0.70),
      min_se_available = NULL
    ),
    list(max_rmse = NULL, max_abs_bias = NULL, min_se_available = NULL)
  )
}

mfrmr_recovery_validation_resolve_cases <- function(case_ids = NULL,
                                                   tier = c("core", "extended", "all")) {
  plan <- mfrmr_recovery_validation_case_plan()
  if (!is.null(case_ids)) {
    case_ids <- as.character(case_ids)
    missing <- setdiff(case_ids, plan$CaseID)
    if (length(missing) > 0L) {
      stop("Unknown recovery-validation case(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    return(plan[match(case_ids, plan$CaseID), , drop = FALSE])
  }
  tier <- match.arg(tier)
  if (identical(tier, "all")) {
    plan
  } else {
    plan[plan$Tier == tier, , drop = FALSE]
  }
}

mfrmr_run_recovery_validation <- function(case_ids = NULL,
                                          tier = c("core", "extended", "all"),
                                          reps = NULL,
                                          quick = FALSE,
                                          seed = 20260510L,
                                          output_dir = NULL,
                                          write_files = !is.null(output_dir),
                                          verbose = TRUE) {
  if (!requireNamespace("mfrmr", quietly = TRUE)) {
    stop("The `mfrmr` package must be available before running validation.", call. = FALSE)
  }
  tier <- match.arg(tier)
  plan <- mfrmr_recovery_validation_resolve_cases(case_ids = case_ids, tier = tier)
  eval_recovery <- mfrmr_recovery_validation_fun("evaluate_mfrm_recovery")
  assess_recovery <- mfrmr_recovery_validation_fun("assess_mfrm_recovery")
  prompt_steps <- mfrmr_recovery_validation_prompt_steps()
  started_at <- Sys.time()

  results <- vector("list", nrow(plan))
  names(results) <- plan$CaseID
  for (i in seq_len(nrow(plan))) {
    case <- plan[i, , drop = FALSE]
    case_reps <- if (isTRUE(quick)) case$SmokeReps else case$RecommendedReps
    if (!is.null(reps)) case_reps <- as.integer(reps[1])
    case_seed <- as.integer(seed) + i - 1L
    if (isTRUE(verbose)) {
      message("[", case$CaseID, "] running ", case_reps, " replication(s)")
    }

    spec <- mfrmr_recovery_validation_spec(case$CaseID)
    thresholds <- mfrmr_recovery_validation_case_thresholds(case$CaseID)
    recovery <- tryCatch(
      eval_recovery(
        sim_spec = spec,
        reps = case_reps,
        fit_method = case$Estimation,
        maxit = case$Maxit,
        quad_points = if (is.na(case$QuadPoints)) 7L else case$QuadPoints,
        include_person = isTRUE(case$IncludePerson),
        seed = case_seed
      ),
      error = function(e) e
    )

    assessment <- NULL
    if (!inherits(recovery, "error")) {
      assessment <- tryCatch(
        assess_recovery(
          recovery,
          min_reps = case_reps,
          min_se_available = thresholds$min_se_available,
          max_mcse_rmse_ratio = if (case_reps <= 2L) NULL else 0.35,
          max_rmse = thresholds$max_rmse,
          max_abs_bias = thresholds$max_abs_bias
        ),
        error = function(e) e
      )
    }

    results[[i]] <- list(
      case = case,
      sim_spec = spec,
      recovery = recovery,
      assessment = assessment,
      seed = case_seed,
      reps = case_reps,
      thresholds = thresholds,
      error = if (inherits(recovery, "error")) conditionMessage(recovery) else if (inherits(assessment, "error")) conditionMessage(assessment) else NA_character_
    )
  }

  out <- structure(
    list(
      prompt_steps = prompt_steps,
      plan = plan,
      results = results,
      started_at = started_at,
      completed_at = Sys.time(),
      quick = isTRUE(quick),
      seed = seed
    ),
    class = "mfrmr_recovery_validation"
  )
  if (isTRUE(write_files)) {
    mfrmr_write_recovery_validation_outputs(out, output_dir = output_dir)
  }
  out
}

mfrmr_recovery_validation_status_rank <- function(x) {
  ranks <- c(ok = 1L, not_assessed = 2L, not_available = 2L,
             review = 3L, concern = 4L, fail = 4L)
  out <- unname(ranks[as.character(x)])
  out[is.na(out)] <- 0L
  out
}

mfrmr_recovery_validation_compact_status <- function(status,
                                                    unavailable = c("ignore", "review")) {
  unavailable <- match.arg(unavailable)
  status <- as.character(status)
  status <- status[!is.na(status) & nzchar(status)]
  if (length(status) == 0L) return("not_available")
  if (any(status %in% c("concern", "fail"))) return("concern")
  if (any(status == "review")) return("review")
  unavailable_values <- c("not_available", "not_assessed")
  if (any(status %in% unavailable_values)) {
    if (identical(unavailable, "review")) return("review")
    if (all(status %in% unavailable_values)) return("not_available")
  }
  "ok"
}

mfrmr_recovery_validation_domain_status <- function(metric) {
  metric <- as.data.frame(metric, stringsAsFactors = FALSE)
  if (nrow(metric) == 0L) {
    return(list(
      recovery_metric_status = "not_available",
      uncertainty_status = "not_available",
      monte_carlo_status = "not_available"
    ))
  }
  recovery_status <- mfrmr_recovery_validation_compact_status(c(
    metric$RMSEStatus %||% character(0),
    metric$BiasStatus %||% character(0)
  ))
  uncertainty_status <- mfrmr_recovery_validation_compact_status(c(
    metric$CoverageStatus %||% character(0),
    metric$SEStatus %||% character(0)
  ), unavailable = "review")
  monte_carlo_status <- mfrmr_recovery_validation_compact_status(
    metric$MonteCarloStatus %||% character(0)
  )
  list(
    recovery_metric_status = recovery_status,
    uncertainty_status = uncertainty_status,
    monte_carlo_status = monte_carlo_status
  )
}

mfrmr_recovery_validation_release_status <- function(recovery_metric_status,
                                                     monte_carlo_status,
                                                     successful_runs,
                                                     converged_runs,
                                                     reps) {
  recovery_metric_status <- as.character(recovery_metric_status[1] %||% "not_available")
  monte_carlo_status <- as.character(monte_carlo_status[1] %||% "not_available")
  successful_runs <- suppressWarnings(as.numeric(successful_runs[1] %||% NA_real_))
  converged_runs <- suppressWarnings(as.numeric(converged_runs[1] %||% NA_real_))
  reps <- suppressWarnings(as.numeric(reps[1] %||% NA_real_))

  run_ok <- is.finite(reps) && reps > 0 &&
    is.finite(successful_runs) && successful_runs >= reps &&
    is.finite(converged_runs) && converged_runs >= reps
  if (!run_ok) return("concern")
  if (recovery_metric_status %in% c("concern", "fail", "not_available", "not_assessed")) {
    return("concern")
  }
  if (monte_carlo_status %in% c("concern", "fail")) {
    return("concern")
  }
  if (recovery_metric_status == "review" ||
      monte_carlo_status %in% c("review", "not_available", "not_assessed")) {
    return("review")
  }
  "ok"
}

mfrmr_recovery_validation_release_interpretation <- function(status) {
  status <- as.character(status[1] %||% "not_available")
  switch(
    status,
    ok = "Recovery metrics, convergence, and Monte Carlo precision support release-level recovery evidence.",
    review = "Recovery evidence is usable but needs review before it is treated as release-level evidence.",
    concern = "Recovery validation did not meet release-level recovery, convergence, or Monte Carlo criteria.",
    "Release-recovery status is not available."
  )
}

mfrmr_recovery_validation_uncertainty_limitation <- function(uncertainty_status) {
  uncertainty_status <- as.character(uncertainty_status[1] %||% "not_available")
  if (uncertainty_status %in% c("ok")) {
    return("No separate uncertainty limitation flagged by this summary.")
  }
  "SE/coverage evidence is reported separately and should not be read as a recovery-metric failure."
}

mfrmr_recovery_validation_add_release_status <- function(case_summary) {
  case_summary <- as.data.frame(case_summary, stringsAsFactors = FALSE)
  if (nrow(case_summary) == 0L) return(case_summary)
  release_status <- vapply(seq_len(nrow(case_summary)), function(i) {
    mfrmr_recovery_validation_release_status(
      recovery_metric_status = case_summary$RecoveryMetricStatus[i],
      monte_carlo_status = case_summary$MonteCarloStatus[i],
      successful_runs = case_summary$SuccessfulRuns[i],
      converged_runs = case_summary$ConvergedRuns[i],
      reps = case_summary$Reps[i]
    )
  }, character(1))
  release_interpretation <- vapply(
    release_status,
    mfrmr_recovery_validation_release_interpretation,
    character(1)
  )
  uncertainty_limitation <- vapply(
    case_summary$UncertaintyStatus,
    mfrmr_recovery_validation_uncertainty_limitation,
    character(1)
  )
  out <- cbind(
    case_summary[, c("CaseID", "Tier", "Model", "Estimation", "Reps",
                     "SuccessfulRuns", "ConvergedRuns"), drop = FALSE],
    ReleaseRecoveryStatus = release_status,
    RecoveryMetricStatus = case_summary$RecoveryMetricStatus,
    MonteCarloStatus = case_summary$MonteCarloStatus,
    UncertaintyStatus = case_summary$UncertaintyStatus,
    OverallStatus = case_summary$OverallStatus,
    ConcernMetrics = case_summary$ConcernMetrics,
    ReviewMetrics = case_summary$ReviewMetrics,
    MaxRMSE = case_summary$MaxRMSE,
    MaxAbsBias = case_summary$MaxAbsBias,
    WorstMetric = case_summary$WorstMetric,
    ReleaseInterpretation = release_interpretation,
    UncertaintyLimitation = uncertainty_limitation,
    NextAction = case_summary$NextAction
  )
  row.names(out) <- NULL
  out
}

mfrmr_recovery_validation_release_decision <- function(case_summary) {
  case_summary <- as.data.frame(case_summary, stringsAsFactors = FALSE)
  if (nrow(case_summary) == 0L) {
    return(data.frame(
      Cases = 0L,
      ReleaseRecoveryStatus = "not_available",
      CasesOK = 0L,
      CasesReview = 0L,
      CasesConcern = 0L,
      RecoveryMetricStatus = "not_available",
      MonteCarloStatus = "not_available",
      UncertaintyStatus = "not_available",
      PrimaryDecisionBasis = "recovery metrics, convergence, and Monte Carlo precision",
      Conclusion = "No validation cases were summarized.",
      stringsAsFactors = FALSE
    ))
  }
  release_status <- as.character(case_summary$ReleaseRecoveryStatus %||% "not_available")
  overall_release_status <- mfrmr_recovery_validation_compact_status(
    release_status,
    unavailable = "review"
  )
  conclusion <- switch(
    overall_release_status,
    ok = "Release recovery evidence is supported for the summarized cases.",
    review = "Release recovery evidence is partly supported, but at least one case needs review.",
    concern = "Release recovery evidence is not yet sufficient for at least one summarized case.",
    "Release recovery evidence could not be summarized."
  )
  data.frame(
    Cases = nrow(case_summary),
    ReleaseRecoveryStatus = overall_release_status,
    CasesOK = sum(release_status == "ok", na.rm = TRUE),
    CasesReview = sum(release_status == "review", na.rm = TRUE),
    CasesConcern = sum(release_status == "concern", na.rm = TRUE),
    RecoveryMetricStatus = mfrmr_recovery_validation_compact_status(
      case_summary$RecoveryMetricStatus,
      unavailable = "review"
    ),
    MonteCarloStatus = mfrmr_recovery_validation_compact_status(
      case_summary$MonteCarloStatus,
      unavailable = "review"
    ),
    UncertaintyStatus = mfrmr_recovery_validation_compact_status(
      case_summary$UncertaintyStatus,
      unavailable = "review"
    ),
    PrimaryDecisionBasis = "recovery metrics, convergence, and Monte Carlo precision",
    Conclusion = conclusion,
    stringsAsFactors = FALSE
  )
}

mfrmr_summarize_recovery_validation <- function(x) {
  if (!inherits(x, "mfrmr_recovery_validation")) {
    stop("`x` must be output from mfrmr_run_recovery_validation().", call. = FALSE)
  }

  case_rows <- lapply(x$results, function(res) {
    case <- res$case
    err <- res$error
    assessment <- res$assessment
    if (inherits(assessment, "error") || is.null(assessment)) {
      return(data.frame(
        CaseID = case$CaseID,
        Tier = case$Tier,
        Model = case$Model,
        Estimation = case$Estimation,
        Reps = res$reps,
        SuccessfulRuns = NA_integer_,
        ConvergedRuns = NA_integer_,
        OverallStatus = "concern",
        RecoveryMetricStatus = "concern",
        UncertaintyStatus = "not_available",
        MonteCarloStatus = "not_available",
        ConcernMetrics = NA_integer_,
        ReviewMetrics = NA_integer_,
        MaxRMSE = NA_real_,
        MaxAbsBias = NA_real_,
        WorstMetric = NA_character_,
        NextAction = ifelse(is.na(err), "Review failed validation output.", err),
        stringsAsFactors = FALSE
      ))
    }
    overview <- as.data.frame(assessment$overview, stringsAsFactors = FALSE)
    metric <- as.data.frame(assessment$metric_review, stringsAsFactors = FALSE)
    domain_status <- mfrmr_recovery_validation_domain_status(metric)
    worst_metric <- NA_character_
    if (nrow(metric) > 0L && "OverallStatus" %in% names(metric)) {
      worst_idx <- order(
        mfrmr_recovery_validation_status_rank(metric$OverallStatus),
        suppressWarnings(as.numeric(metric$RMSE)),
        decreasing = TRUE
      )[1]
      worst_metric <- paste(
        metric$ParameterType[worst_idx],
        metric$Facet[worst_idx],
        metric$ComparisonScale[worst_idx],
        sep = " / "
      )
    }
    data.frame(
      CaseID = case$CaseID,
      Tier = case$Tier,
      Model = case$Model,
      Estimation = case$Estimation,
      Reps = res$reps,
      SuccessfulRuns = overview$SuccessfulRuns[1] %||% NA_integer_,
      ConvergedRuns = overview$ConvergedRuns[1] %||% NA_integer_,
      OverallStatus = overview$OverallStatus[1] %||% "not_available",
      RecoveryMetricStatus = domain_status$recovery_metric_status,
      UncertaintyStatus = domain_status$uncertainty_status,
      MonteCarloStatus = domain_status$monte_carlo_status,
      ConcernMetrics = if (nrow(metric) > 0L) sum(metric$OverallStatus == "concern", na.rm = TRUE) else 0L,
      ReviewMetrics = if (nrow(metric) > 0L) sum(metric$OverallStatus == "review", na.rm = TRUE) else 0L,
      MaxRMSE = if (nrow(metric) > 0L) max(suppressWarnings(as.numeric(metric$RMSE)), na.rm = TRUE) else NA_real_,
      MaxAbsBias = if (nrow(metric) > 0L) max(abs(suppressWarnings(as.numeric(metric$Bias))), na.rm = TRUE) else NA_real_,
      WorstMetric = worst_metric,
      NextAction = paste(assessment$next_actions %||% character(0), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  case_summary <- mfrmr_recovery_validation_add_release_status(do.call(rbind, case_rows))
  release_decision_table <- case_summary[, c(
    "CaseID",
    "Tier",
    "Model",
    "ReleaseRecoveryStatus",
    "RecoveryMetricStatus",
    "MonteCarloStatus",
    "UncertaintyStatus",
    "OverallStatus",
    "ReleaseInterpretation",
    "UncertaintyLimitation"
  ), drop = FALSE]
  topline_release_decision <- mfrmr_recovery_validation_release_decision(case_summary)

  metric_rows <- lapply(x$results, function(res) {
    assessment <- res$assessment
    if (inherits(assessment, "error") || is.null(assessment)) return(NULL)
    metric <- as.data.frame(assessment$metric_review, stringsAsFactors = FALSE)
    if (nrow(metric) == 0L) return(NULL)
    case <- res$case
    out <- cbind(
      data.frame(
        CaseID = case$CaseID,
        Tier = case$Tier,
        Model = case$Model,
        Estimation = case$Estimation,
        stringsAsFactors = FALSE
      ),
      metric
    )
    row.names(out) <- NULL
    out
  })
  metric_summary <- do.call(rbind, metric_rows)
  if (is.null(metric_summary)) metric_summary <- data.frame()

  decision_table <- if (nrow(metric_summary) > 0L) {
    as.data.frame(
      stats::xtabs(~ CaseID + OverallStatus, data = metric_summary),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(CaseID = character(0), OverallStatus = character(0), Freq = integer(0))
  }

  domain_decision_table <- if (nrow(case_summary) > 0L) {
    domain_rows <- list(
      data.frame(CaseID = case_summary$CaseID,
                 StatusDomain = "recovery_metrics",
                 Status = case_summary$RecoveryMetricStatus,
                 stringsAsFactors = FALSE),
      data.frame(CaseID = case_summary$CaseID,
                 StatusDomain = "uncertainty",
                 Status = case_summary$UncertaintyStatus,
                 stringsAsFactors = FALSE),
      data.frame(CaseID = case_summary$CaseID,
                 StatusDomain = "monte_carlo",
                 Status = case_summary$MonteCarloStatus,
                 stringsAsFactors = FALSE),
      data.frame(CaseID = case_summary$CaseID,
                 StatusDomain = "overall",
                 Status = case_summary$OverallStatus,
                 stringsAsFactors = FALSE)
    )
    do.call(rbind, domain_rows)
  } else {
    data.frame(CaseID = character(0), StatusDomain = character(0), Status = character(0))
  }

  run_notes <- do.call(rbind, lapply(x$results, function(res) {
    rec <- res$recovery
    notes <- if (inherits(rec, "error")) conditionMessage(rec) else paste(rec$notes %||% character(0), collapse = " | ")
    data.frame(
      CaseID = res$case$CaseID,
      Notes = notes,
      Error = res$error,
      stringsAsFactors = FALSE
    )
  }))

  list(
    prompt_steps = x$prompt_steps,
    case_plan = x$plan,
    topline_release_decision = topline_release_decision,
    release_decision_table = release_decision_table,
    case_summary = case_summary,
    metric_summary = metric_summary,
    decision_table = decision_table,
    domain_decision_table = domain_decision_table,
    run_notes = run_notes,
    started_at = x$started_at,
    completed_at = x$completed_at
  )
}

mfrmr_validation_markdown_table <- function(df, max_rows = 20L) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) == 0L || ncol(df) == 0L) return("_No rows._")
  df <- df[seq_len(min(nrow(df), max_rows)), , drop = FALSE]
  df[] <- lapply(df, function(col) {
    col <- as.character(col)
    col[is.na(col)] <- ""
    gsub("\\|", "/", col)
  })
  header <- paste(names(df), collapse = " | ")
  rule <- paste(rep("---", ncol(df)), collapse = " | ")
  rows <- apply(df, 1L, paste, collapse = " | ")
  paste(c(paste0("| ", header, " |"), paste0("| ", rule, " |"), paste0("| ", rows, " |")), collapse = "\n")
}

mfrmr_recovery_validation_markdown <- function(x) {
  s <- mfrmr_summarize_recovery_validation(x)
  lines <- c(
    "# mfrmr recovery validation summary",
    "",
    paste0("- Started: ", format(s$started_at, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- Completed: ", format(s$completed_at, "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## Top-line release decision",
    "",
    mfrmr_validation_markdown_table(s$topline_release_decision, max_rows = 5L),
    "",
    "## Release decision by case",
    "",
    mfrmr_validation_markdown_table(s$release_decision_table, max_rows = 20L),
    "",
    "## Review steps",
    "",
    mfrmr_validation_markdown_table(s$prompt_steps, max_rows = 20L),
    "",
    "## Case summary",
    "",
    mfrmr_validation_markdown_table(s$case_summary, max_rows = 20L),
    "",
    "## Decision table",
    "",
    mfrmr_validation_markdown_table(s$decision_table, max_rows = 50L),
    "",
    "## Domain decision table",
    "",
    mfrmr_validation_markdown_table(s$domain_decision_table, max_rows = 80L),
    "",
    "## Metric summary",
    "",
    mfrmr_validation_markdown_table(s$metric_summary, max_rows = 50L),
    "",
    "## Run notes",
    "",
    mfrmr_validation_markdown_table(s$run_notes, max_rows = 20L)
  )
  paste(lines, collapse = "\n")
}

summary.mfrmr_recovery_validation <- function(object, ...) {
  s <- mfrmr_summarize_recovery_validation(object)
  structure(
    list(
      topline_release_decision = s$topline_release_decision,
      release_decision_table = s$release_decision_table,
      case_summary = s$case_summary,
      domain_decision_table = s$domain_decision_table,
      started_at = s$started_at,
      completed_at = s$completed_at
    ),
    class = "summary.mfrmr_recovery_validation"
  )
}

print.mfrmr_recovery_validation <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

print.summary.mfrmr_recovery_validation <- function(x, ...) {
  cat("mfrmr Recovery Validation Summary\n")
  cat("Started: ", format(x$started_at, "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
  cat("Completed: ", format(x$completed_at, "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
  top <- as.data.frame(x$topline_release_decision, stringsAsFactors = FALSE)
  if (nrow(top) > 0L) {
    cat("Release recovery status: ", as.character(top$ReleaseRecoveryStatus[1]), "\n", sep = "")
    cat("Cases: ", as.character(top$Cases[1]),
        " (ok=", as.character(top$CasesOK[1]),
        ", review=", as.character(top$CasesReview[1]),
        ", concern=", as.character(top$CasesConcern[1]), ")\n", sep = "")
    cat("Primary basis: ", as.character(top$PrimaryDecisionBasis[1]), "\n", sep = "")
    cat(as.character(top$Conclusion[1]), "\n", sep = "")
  }
  tbl <- as.data.frame(x$release_decision_table, stringsAsFactors = FALSE)
  keep <- intersect(
    c("CaseID", "ReleaseRecoveryStatus", "RecoveryMetricStatus",
      "MonteCarloStatus", "UncertaintyStatus", "OverallStatus"),
    names(tbl)
  )
  if (nrow(tbl) > 0L && length(keep) > 0L) {
    cat("\nCase decisions:\n")
    display <- tbl[, keep, drop = FALSE]
    names(display) <- c("Case", "Release", "Recovery", "MC",
                        "Uncertainty", "Overall")[seq_along(keep)]
    print(display, row.names = FALSE)
  }
  invisible(x)
}

mfrmr_write_recovery_validation_outputs <- function(x,
                                                   output_dir,
                                                   prefix = "mfrmr_recovery_validation") {
  if (is.null(output_dir) || !nzchar(as.character(output_dir[1]))) {
    stop("`output_dir` must be a non-empty directory path.", call. = FALSE)
  }
  output_dir <- as.character(output_dir[1])
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  s <- mfrmr_summarize_recovery_validation(x)
  utils::write.csv(s$topline_release_decision, file.path(output_dir, paste0(prefix, "_topline_release_decision.csv")), row.names = FALSE)
  utils::write.csv(s$release_decision_table, file.path(output_dir, paste0(prefix, "_release_decision_table.csv")), row.names = FALSE)
  utils::write.csv(s$prompt_steps, file.path(output_dir, paste0(prefix, "_prompt_steps.csv")), row.names = FALSE)
  utils::write.csv(s$case_plan, file.path(output_dir, paste0(prefix, "_case_plan.csv")), row.names = FALSE)
  utils::write.csv(s$case_summary, file.path(output_dir, paste0(prefix, "_case_summary.csv")), row.names = FALSE)
  utils::write.csv(s$metric_summary, file.path(output_dir, paste0(prefix, "_metric_summary.csv")), row.names = FALSE)
  utils::write.csv(s$decision_table, file.path(output_dir, paste0(prefix, "_decision_table.csv")), row.names = FALSE)
  utils::write.csv(s$domain_decision_table, file.path(output_dir, paste0(prefix, "_domain_decision_table.csv")), row.names = FALSE)
  utils::write.csv(s$run_notes, file.path(output_dir, paste0(prefix, "_run_notes.csv")), row.names = FALSE)
  saveRDS(x, file.path(output_dir, paste0(prefix, ".rds")))
  writeLines(mfrmr_recovery_validation_markdown(x), file.path(output_dir, paste0(prefix, ".md")))
  invisible(normalizePath(output_dir, mustWork = FALSE))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
