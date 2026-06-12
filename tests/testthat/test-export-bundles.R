delayed_export_fixture <- function(builder) {
  structure(
    list(
      .builder = builder,
      .cache = new.env(parent = emptyenv())
    ),
    class = "mfrmr_test_delayed_fixture"
  )
}

force_export_fixture <- function(x) {
  if (!inherits(x, "mfrmr_test_delayed_fixture")) {
    return(x)
  }
  parts <- unclass(x)
  cache <- parts[[".cache"]]
  if (!exists("value", envir = cache, inherits = FALSE)) {
    assign("value", parts[[".builder"]](), envir = cache)
  }
  get("value", envir = cache, inherits = FALSE)
}

`$.mfrmr_test_delayed_fixture` <- function(x, name) {
  force_export_fixture(x)[[name]]
}

`[[.mfrmr_test_delayed_fixture` <- function(x, i, ...) {
  force_export_fixture(x)[[i, ...]]
}

export_core_fixture <- delayed_export_fixture(function() {
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

  dat <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    maxit = 20
  ))
  diagnostics <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "overall"))

  list(
    fit = fit,
    diagnostics = diagnostics
  )
})

export_facets_fixture <- delayed_export_fixture(function() {
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

  dat <- mfrmr:::sample_mfrm_data(seed = 123)
  run <- suppressWarnings(run_mfrm_facets(
    dat,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 20
  ))

  list(run = run)
})

export_bias_fixture <- delayed_export_fixture(function() {
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

  bias_all <- suppressWarnings(estimate_bias(
    export_core_fixture$fit,
    diagnostics = export_core_fixture$diagnostics,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))

  list(bias_all = bias_all)
})

weighting_appendix_fixture <- delayed_export_fixture(function() {
  dat <- load_mfrmr_data("example_core")
  keep_people <- unique(dat$Person)[1:12]
  dat <- dat[dat$Person %in% keep_people, , drop = FALSE]

  rasch_fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 7,
    maxit = 25
  ))
  gpcm_fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    quad_points = 7,
    maxit = 25
  ))

  build_weighting_review(rasch_fit, gpcm_fit, theta_points = 21, top_n = 5)
})

misfit_appendix_fixture <- delayed_export_fixture(function() {
  dat <- load_mfrmr_data("example_core")
  keep_people <- unique(dat$Person)[1:12]
  dat <- dat[dat$Person %in% keep_people, , drop = FALSE]

  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 7,
    maxit = 25
  ))
  diagnostics <- suppressWarnings(diagnose_mfrm(
    fit,
    diagnostic_mode = "both",
    residual_pca = "none"
  ))

  build_misfit_casebook(fit, diagnostics = diagnostics, top_n = 5)
})

recovery_appendix_fixture <- delayed_export_fixture(function() {
  recovery <- structure(
    list(
      recovery = data.frame(
        rep = 1L,
        ParameterType = "facet",
        Facet = "Rater",
        Level = "R1",
        Truth = 0.1,
        Estimate = 0.12,
        EstimateAligned = 0.11,
        ErrorAligned = 0.01,
        stringsAsFactors = FALSE
      ),
      recovery_summary = data.frame(
        ParameterType = "facet",
        Facet = "Rater",
        ComparisonScale = "logit",
        Rows = 1L,
        Reps = 1L,
        Bias = 0.01,
        RMSE = 0.01,
        MAE = 0.01,
        Coverage95 = 1,
        McseBias = NA_real_,
        McseRMSE = NA_real_,
        stringsAsFactors = FALSE
      ),
      rep_overview = data.frame(
        rep = 1L,
        RunOK = TRUE,
        Converged = TRUE,
        RecoveryRows = 1L,
        ElapsedSec = 0.2,
        stringsAsFactors = FALSE
      ),
      settings = list(reps = 1L, fit_method = "MML", model = "RSM"),
      notes = "Use more replications before treating this as a final recovery study.",
      ademp = list(
        aims = "Assess parameter recovery under a fixed MFRM design.",
        data_generating_mechanism = list(
          model = "RSM",
          assignment = "rotating",
          step_facet = "Criterion"
        ),
        methods = list(fit_method = "MML", fitted_model = "RSM"),
        estimands = "Facet location recovery",
        performance_measures = c("Bias", "RMSE")
      )
    ),
    class = "mfrm_recovery_simulation"
  )

  assessment <- structure(
    list(
      overview = data.frame(
        Reps = 1L,
        SuccessfulRuns = 1L,
        SuccessRate = 1,
        ConvergedRuns = 1L,
        ConvergenceRate = 1,
        RecoveryRows = 1L,
        RecoveryGroups = 1L,
        OverallStatus = "review",
        stringsAsFactors = FALSE
      ),
      checklist = data.frame(
        Section = "Run completion",
        Item = "Replication count",
        Status = "review",
        Evidence = "1 replication; requested minimum is 30.",
        NextAction = "Increase the replication count.",
        stringsAsFactors = FALSE
      ),
      metric_review = data.frame(
        ParameterType = "facet",
        Facet = "Rater",
        ComparisonScale = "logit",
        RMSE = 0.01,
        Bias = 0.01,
        RMSEStatus = "ok",
        BiasStatus = "ok",
        CoverageStatus = "ok",
        OverallStatus = "ok",
        NextAction = "No immediate action.",
        stringsAsFactors = FALSE
      ),
      next_actions = "Increase the replication count.",
      thresholds = list(min_reps = 30L, max_rmse = c(default = 0.5)),
      notes = "RMSE and bias statuses depend on supplied practical thresholds.",
      source = recovery,
      digits = 3L
    ),
    class = "mfrm_recovery_assessment"
  )

  list(recovery = recovery, assessment = assessment)
})

export_recovery_validation_summary_fixture <- function() {
  validation_summary <- list(
    topline_release_decision = data.frame(
      Cases = 1L,
      ReleaseRecoveryStatus = "review",
      PrimaryDecisionBasis = "recovery metrics, convergence, and Monte Carlo precision",
      stringsAsFactors = FALSE
    ),
    reading_order = data.frame(
      Step = 1L,
      Route = "summary(validation)$topline_release_decision",
      WhatToRead = "Release status.",
      Purpose = "Start with the top-line decision.",
      stringsAsFactors = FALSE
    ),
    release_decision_table = data.frame(
      CaseID = "gpcm_high_dispersion_sparse",
      ReleaseRecoveryStatus = "review",
      RecoveryMetricStatus = "ok",
      UncertaintyStatus = "review",
      ScoreSupportStatus = "review",
      stringsAsFactors = FALSE
    ),
    case_summary = data.frame(
      CaseID = "gpcm_high_dispersion_sparse",
      GPCMSlopeRegime = "high_dispersion",
      ScoreSupportStatus = "review",
      stringsAsFactors = FALSE
    ),
    condition_reporting_notes = data.frame(
      CaseID = "gpcm_high_dispersion_sparse",
      ConditionArea = "score_support",
      ReportingAttention = "reporting_review",
      ConditionFinding = "omitted_generated_score_categories",
      Evidence = "max_zero_score_levels=1",
      ReportingImplication = "Report sparse generated score support as condition stress.",
      NextAction = "Inspect category-level recovery before generalizing.",
      ValidationUse = "generator_condition_not_release_gate",
      stringsAsFactors = FALSE
    ),
    condition_summary = data.frame(
      CaseID = "gpcm_high_dispersion_sparse",
      GPCMSlopeRegime = "high_dispersion",
      ScoreSupportStatus = "review",
      stringsAsFactors = FALSE
    ),
    diagnostic_reporting_notes = data.frame(
      CaseID = "gpcm_high_dispersion_sparse",
      Facet = "Rater",
      ReportingAttention = "reporting_review",
      DiagnosticFinding = "zero_separation_or_reliability",
      Evidence = "mean_separation=0; mean_reliability=0",
      ReportingImplication = "Report zero separation as diagnostic context.",
      NextAction = "Inspect the generated condition before using reliability language.",
      ValidationUse = "diagnostic_only_not_release_gate",
      stringsAsFactors = FALSE
    ),
    diagnostic_oc_summary = data.frame(
      CaseID = "gpcm_high_dispersion_sparse",
      Facet = "Rater",
      MeanSeparation = 0,
      MeanReliability = 0,
      DiagnosticAvailability = "available",
      ValidationUse = "diagnostic_only_not_release_gate",
      stringsAsFactors = FALSE
    ),
    domain_decision_table = data.frame(
      CaseID = "gpcm_high_dispersion_sparse",
      StatusDomain = "score_support",
      Status = "review",
      stringsAsFactors = FALSE
    ),
    started_at = Sys.time(),
    completed_at = Sys.time()
  )
  class(validation_summary) <- "summary.mfrmr_recovery_validation"
  validation_summary
}

diagnostic_screening_appendix_summary_fixture <- function() {
  out <- list(
    overview = data.frame(
      Designs = 1L,
      Reps = 1L,
      Scenarios = "well_specified, local_dependence",
      Models = "RSM",
      ReplicateRows = 2L,
      ScenarioRows = 2L,
      PerformanceRows = 2L,
      ReportSignalRows = 2L,
      ContrastRows = 1L,
      RunOKRate = 1,
      ConvergenceRate = 1,
      IncludeReport = TRUE,
      PlotDataContract = "mfrm_plot_data",
      stringsAsFactors = FALSE
    ),
    reading_order = data.frame(
      Step = 1L,
      Table = "scenario_summary",
      WhatToRead = "Scenario-by-design screening summaries.",
      Purpose = "Compare screening surfaces across scenarios.",
      InterpretationBoundary = "Scenario means and flag rates are operating-characteristic readouts.",
      stringsAsFactors = FALSE
    ),
    next_actions = data.frame(
      Priority = 1L,
      Area = "Appendix and plot-data handoff",
      Status = "ok",
      Evidence = "Summary tables and draw-free plot-data tables are available.",
      Action = "Use build_summary_table_bundle() or export_summary_appendix().",
      Route = "build_summary_table_bundle(diag_eval); export_summary_appendix(diag_eval)",
      ReportingBoundary = "Exports and plot data are presentation handoffs over the same summary evidence.",
      stringsAsFactors = FALSE
    ),
    reporting_notes = data.frame(
      Area = "Diagnostic-screening scope",
      Evidence = "Repeated simulation, fitting, diagnosis, and aggregation.",
      ReportingBoundary = "Describe as an operating-characteristic study, not a calibrated hypothesis test.",
      RecommendedAction = "Report scenarios, design grid, replication count, and fitting method.",
      stringsAsFactors = FALSE
    ),
    figure_recipes = data.frame(
      FigureID = "overview_rates",
      RecommendedUse = "main_text_or_primary_supplement",
      PrimaryQuestion = "How often do screening signals fire across scenarios?",
      PlotCall = "plot(diag_eval, type = \"overview\", metric = \"rate\", draw = FALSE)",
      PlotDataCall = "plot_data(diag_eval, type = \"overview\", metric = \"rate\", component = \"plot_long\")",
      SummaryTable = "plot_overview_rate",
      DisplaySuggestion = "Line or point plot by design variable.",
      CaptionFocus = "Describe operating-characteristic signal rates.",
      InterpretationBoundary = "Rates are simulation summaries, not calibrated inferential tests.",
      Availability = "available_when_plot_rows_exist",
      stringsAsFactors = FALSE
    ),
    scenario_summary = data.frame(
      design_id = "D1",
      Scenario = c("well_specified", "local_dependence"),
      n_person = 12L,
      LegacyAnyFlagRate = c(0, 1),
      MarginalAnyFlagRate = c(0, 1),
      PairwiseAnyFlagRate = c(0, 1),
      stringsAsFactors = FALSE
    ),
    performance_summary = data.frame(
      design_id = "D1",
      Scenario = c("well_specified", "local_dependence"),
      n_person = 12L,
      EvaluationUse = c("type_I_proxy", "sensitivity_proxy"),
      StrictAnyFlagRate = c(0, 1),
      MeanElapsedSec = c(0.10, 0.12),
      stringsAsFactors = FALSE
    ),
    report_signal_summary = data.frame(
      design_id = "D1",
      Scenario = c("well_specified", "local_dependence"),
      n_person = 12L,
      ReportIndexAvailabilityRate = c(1, 1),
      FitReportReviewRate = c(0, 1),
      stringsAsFactors = FALSE
    ),
    scenario_contrast = data.frame(
      design_id = "D1",
      Scenario = "local_dependence",
      n_person = 12L,
      DeltaLegacyFlaggedLevels = 2,
      DeltaMarginalFlaggedGroups = 1,
      stringsAsFactors = FALSE
    ),
    plot_overview_rate = data.frame(
      design_id = "D1",
      Scenario = c("well_specified", "local_dependence"),
      n_person = 12L,
      Signal = "Strict combined any-flag rate",
      Value = c(0, 1),
      SourceTable = "performance_summary",
      Metric = "StrictAnyFlagRate",
      stringsAsFactors = FALSE
    ),
    plot_overview_count = data.frame(),
    plot_report_rate = data.frame(
      design_id = "D1",
      Scenario = "local_dependence",
      n_person = 12L,
      Signal = "Fit review rate",
      Value = 1,
      SourceTable = "report_signal_summary",
      Metric = "FitReportReviewRate",
      stringsAsFactors = FALSE
    ),
    plot_contrast_count = data.frame(),
    plot_runtime = data.frame(
      design_id = "D1",
      Scenario = c("well_specified", "local_dependence"),
      n_person = 12L,
      Signal = "Mean elapsed seconds",
      Value = c(0.10, 0.12),
      SourceTable = "performance_summary",
      Metric = "MeanElapsedSec",
      stringsAsFactors = FALSE
    ),
    ademp = list(
      aims = "Screen diagnostic operating characteristics",
      data_generating_mechanism = list(model = "RSM", assignment = "complete", step_facet = NA_character_),
      methods = list(fit_method = "MML", fitted_model = "RSM"),
      estimands = "Screening signal rates",
      performance_measures = "Type I and sensitivity proxies"
    ),
    settings = list(reps = 1L, model = "RSM", include_report = TRUE),
    notes = "Draw-free diagnostic-screening plot tables are exported as operating-characteristic readouts."
  )
  class(out) <- "summary.mfrm_diagnostic_screening"
  out
}

prediction_bundle_fixture <- delayed_export_fixture(function() {
  dat <- load_mfrmr_data("example_core")
  keep_people <- unique(dat$Person)[1:18]
  dat <- dat[dat$Person %in% keep_people, , drop = FALSE]

  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    quad_points = 5,
    maxit = 15
  ))
  diagnostics <- diagnose_mfrm(fit, residual_pca = "none")

  spec <- build_mfrm_sim_spec(
    n_person = 20,
    n_rater = 4,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating"
  )

  population_prediction <- suppressWarnings(
    predict_mfrm_population(
      sim_spec = spec,
      reps = 2,
      maxit = 15,
      seed = 1
    )
  )

  new_units <- data.frame(
    Person = c("NEW01", "NEW01"),
    Rater = unique(dat$Rater)[1],
    Criterion = unique(dat$Criterion)[1:2],
    Score = c(2, 3)
  )

  unit_prediction <- predict_mfrm_units(
    fit,
    new_units,
    n_draws = 2,
    seed = 1
  )
  plausible_values <- sample_mfrm_plausible_values(
    fit,
    new_units,
    n_draws = 2,
    seed = 1
  )

  list(
    fit = fit,
    diagnostics = diagnostics,
    population_prediction = population_prediction,
    unit_prediction = unit_prediction,
    plausible_values = plausible_values
  )
})

latent_prediction_bundle_fixture <- delayed_export_fixture(function() {
  fixture <- mfrmr:::with_preserved_rng_seed(20260403, {
    persons <- paste0("P", sprintf("%02d", 1:60))
    items <- paste0("I", 1:6)
    x <- stats::rnorm(length(persons))
    theta <- 0.25 + 0.9 * x + stats::rnorm(length(persons), sd = 0.6)
    item_beta <- seq(-1.0, 1.0, length.out = length(items))

    dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
    eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
    dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))

    person_tbl <- data.frame(
      Person = persons,
      X = x,
      stringsAsFactors = FALSE
    )

    fit <- suppressWarnings(fit_mfrm(
      dat,
      "Person",
      "Item",
      "Score",
      method = "MML",
      model = "RSM",
      population_formula = ~ X,
      person_data = person_tbl,
      quad_points = 7,
      maxit = 80
    ))

    new_units <- data.frame(
      Person = c("NEW_LOW", "NEW_LOW", "NEW_HIGH", "NEW_HIGH"),
      Item = c(items[1], items[2], items[1], items[2]),
      Score = c(1, 0, 1, 0),
      stringsAsFactors = FALSE
    )
    new_person_data <- data.frame(
      Person = c("NEW_LOW", "NEW_HIGH"),
      X = c(-1.5, 1.5),
      stringsAsFactors = FALSE
    )

    unit_prediction <- predict_mfrm_units(
      fit,
      new_units,
      person_data = new_person_data,
      n_draws = 2,
      seed = 1
    )
    plausible_values <- sample_mfrm_plausible_values(
      fit,
      new_units,
      person_data = new_person_data,
      n_draws = 2,
      seed = 1
    )

    list(
      fit = fit,
      unit_prediction = unit_prediction,
      plausible_values = plausible_values
    )
  })

  fixture
})

latent_prediction_omit_bundle_fixture <- delayed_export_fixture(function() {
  fixture <- mfrmr:::with_preserved_rng_seed(20260411, {
    persons <- paste0("P", sprintf("%02d", 1:60))
    items <- paste0("I", 1:6)
    x <- stats::rnorm(length(persons))
    theta <- 0.15 + 0.8 * x + stats::rnorm(length(persons), sd = 0.65)
    item_beta <- seq(-1.0, 1.0, length.out = length(items))

    dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
    eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
    dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))

    fit_person_tbl <- data.frame(
      Person = persons,
      X = x,
      stringsAsFactors = FALSE
    )
    fit_person_tbl$X[1] <- NA_real_

    fit <- suppressWarnings(fit_mfrm(
      dat,
      "Person",
      "Item",
      "Score",
      method = "MML",
      model = "RSM",
      population_formula = ~ X,
      person_data = fit_person_tbl,
      population_policy = "omit",
      quad_points = 7,
      maxit = 80
    ))

    new_units <- data.frame(
      Person = c("NEW_LOW", "NEW_LOW", "NEW_HIGH", "NEW_HIGH"),
      Item = c(items[1], items[2], items[1], items[2]),
      Score = c(1, 0, 1, 0),
      stringsAsFactors = FALSE
    )
    scoring_person_tbl <- data.frame(
      Person = c("NEW_LOW", "NEW_HIGH"),
      X = c(NA_real_, 1.5),
      stringsAsFactors = FALSE
    )

    unit_prediction <- predict_mfrm_units(
      fit,
      new_units,
      person_data = scoring_person_tbl,
      population_policy = "omit",
      n_draws = 2,
      seed = 2
    )
    plausible_values <- sample_mfrm_plausible_values(
      fit,
      new_units,
      person_data = scoring_person_tbl,
      population_policy = "omit",
      n_draws = 2,
      seed = 2
    )

    list(
      fit = fit,
      unit_prediction = unit_prediction,
      plausible_values = plausible_values
    )
  })

  fixture
})

test_that("build_mfrm_manifest captures reproducibility metadata", {
  manifest <- build_mfrm_manifest(
    fit = export_core_fixture$fit,
    diagnostics = export_core_fixture$diagnostics,
    bias_results = export_bias_fixture$bias_all
  )

  expect_s3_class(manifest, "mfrm_manifest")
  expect_true(is.data.frame(manifest$summary))
  expect_true(is.data.frame(manifest$environment))
  expect_true(is.data.frame(manifest$available_outputs))
  expect_true(any(manifest$available_outputs$Component == "residual_pca"))
  expect_true(
    manifest$available_outputs$Available[manifest$available_outputs$Component == "bias_results"][1]
  )
  expect_equal(manifest$summary$Method[[1]], "JML")
  expect_equal(manifest$summary$MethodUsed[[1]], "JMLE")
  expect_equal(manifest$summary$Observations[[1]], nrow(export_core_fixture$fit$prep$data))
  expect_equal(manifest$summary$Persons[[1]], export_core_fixture$fit$config$n_person)
})

test_that("build_mfrm_manifest and replay script support FACETS-mode runs", {
  manifest <- build_mfrm_manifest(export_facets_fixture$run)
  replay <- build_mfrm_replay_script(
    export_facets_fixture$run,
    bias_results = export_bias_fixture$bias_all,
    data_file = "analysis_data.csv"
  )

  expect_s3_class(manifest, "mfrm_manifest")
  expect_s3_class(replay, "mfrm_replay_script")
  expect_match(replay$script, "run_mfrm_facets\\(")
  expect_match(replay$script, "analysis_data\\.csv")
  expect_match(replay$script, "estimate_bias\\(")
  expect_match(replay$script, "# posterior_basis = legacy_mml", fixed = TRUE)
})

test_that("build_mfrm_manifest records optional prediction artifacts", {
  manifest <- build_mfrm_manifest(
    fit = prediction_bundle_fixture$fit,
    diagnostics = prediction_bundle_fixture$diagnostics,
    population_prediction = prediction_bundle_fixture$population_prediction,
    unit_prediction = prediction_bundle_fixture$unit_prediction,
    plausible_values = prediction_bundle_fixture$plausible_values
  )

  expect_s3_class(manifest, "mfrm_manifest")
  expect_true(any(manifest$available_outputs$Component == "population_prediction"))
  expect_true(any(manifest$available_outputs$Component == "unit_prediction"))
  expect_true(any(manifest$available_outputs$Component == "plausible_values"))
  expect_true(
    manifest$available_outputs$Available[manifest$available_outputs$Component == "population_prediction"][1]
  )
  expect_true(
    manifest$available_outputs$Available[manifest$available_outputs$Component == "unit_prediction"][1]
  )
  expect_true(
    manifest$available_outputs$Available[manifest$available_outputs$Component == "plausible_values"][1]
  )
})

test_that("build_mfrm_manifest records latent-regression provenance", {
  manifest <- build_mfrm_manifest(
    fit = latent_prediction_bundle_fixture$fit,
    unit_prediction = latent_prediction_bundle_fixture$unit_prediction,
    plausible_values = latent_prediction_bundle_fixture$plausible_values
  )

  value_of <- function(tbl, key) {
    as.character(tbl$Value[tbl$Setting == key][1])
  }

  expect_s3_class(manifest, "mfrm_manifest")
  expect_true(isTRUE(manifest$summary$FitPopulationActive[[1]]))
  expect_identical(as.character(manifest$summary$FitPosteriorBasis[[1]]), "population_model")
  expect_identical(value_of(manifest$model_settings, "population_active"), "TRUE")
  expect_identical(value_of(manifest$model_settings, "posterior_basis"), "population_model")
  expect_match(value_of(manifest$model_settings, "population_formula"), "~\\s*X")
  expect_identical(value_of(manifest$model_settings, "population_person_id"), "Person")
  expect_identical(value_of(manifest$model_settings, "population_policy"), "error")
  expect_match(value_of(manifest$model_settings, "population_design_columns"), "\\(Intercept\\)")
  expect_match(value_of(manifest$model_settings, "population_design_columns"), "X")
  expect_identical(value_of(manifest$model_settings, "population_xlevel_variables"), "")
  expect_identical(value_of(manifest$model_settings, "population_contrast_variables"), "")
  expect_equal(as.integer(value_of(manifest$model_settings, "population_person_rows")), 60L)
  expect_equal(as.integer(value_of(manifest$model_settings, "population_person_rows_replay")), 60L)
  expect_identical(value_of(manifest$model_settings, "population_person_replay_scope"), "observed_person_subset_pre_omit")
  expect_equal(as.integer(value_of(manifest$model_settings, "population_omitted_persons")), 0L)
  expect_equal(as.integer(value_of(manifest$model_settings, "population_response_rows_omitted")), 0L)
  expect_identical(value_of(manifest$settings, "fit_population_active"), "TRUE")
  expect_identical(value_of(manifest$settings, "fit_posterior_basis"), "population_model")
  expect_match(value_of(manifest$settings, "fit_population_formula"), "~\\s*X")
  expect_identical(value_of(manifest$settings, "fit_population_person_id"), "Person")
  expect_identical(value_of(manifest$settings, "fit_population_policy"), "error")
  expect_identical(value_of(manifest$settings, "fit_population_xlevel_variables"), "")
  expect_identical(value_of(manifest$settings, "fit_population_contrast_variables"), "")
  expect_equal(as.integer(value_of(manifest$settings, "fit_population_person_rows_replay")), 60L)
  expect_identical(value_of(manifest$settings, "fit_population_person_replay_scope"), "observed_person_subset_pre_omit")
  expect_equal(as.integer(value_of(manifest$settings, "fit_population_omitted_persons")), 0L)
  expect_equal(as.integer(value_of(manifest$settings, "fit_population_response_rows_omitted")), 0L)
  expect_identical(value_of(manifest$settings, "unit_prediction_posterior_basis"), "population_model")
  expect_identical(value_of(manifest$settings, "unit_prediction_person_id"), "Person")
  expect_identical(value_of(manifest$settings, "unit_prediction_population_policy"), "error")
  expect_match(value_of(manifest$settings, "unit_prediction_population_formula"), "~\\s*X")
  expect_identical(value_of(manifest$settings, "plausible_value_posterior_basis"), "population_model")
  expect_identical(value_of(manifest$settings, "plausible_value_person_id"), "Person")
  expect_identical(value_of(manifest$settings, "plausible_value_population_policy"), "error")
  expect_match(value_of(manifest$settings, "plausible_value_population_formula"), "~\\s*X")
})

test_that("build_mfrm_manifest records latent-regression omit provenance", {
  manifest <- build_mfrm_manifest(
    fit = latent_prediction_omit_bundle_fixture$fit,
    unit_prediction = latent_prediction_omit_bundle_fixture$unit_prediction,
    plausible_values = latent_prediction_omit_bundle_fixture$plausible_values
  )

  value_of <- function(tbl, key) {
    as.character(tbl$Value[tbl$Setting == key][1])
  }

  expect_s3_class(manifest, "mfrm_manifest")
  expect_true(isTRUE(manifest$summary$FitPopulationActive[[1]]))
  expect_identical(as.character(manifest$summary$FitPosteriorBasis[[1]]), "population_model")
  expect_identical(value_of(manifest$model_settings, "population_person_id"), "Person")
  expect_identical(value_of(manifest$model_settings, "population_policy"), "omit")
  expect_equal(as.integer(value_of(manifest$model_settings, "population_person_rows")), 59L)
  expect_equal(as.integer(value_of(manifest$model_settings, "population_person_rows_replay")), 60L)
  expect_identical(value_of(manifest$model_settings, "population_person_replay_scope"), "observed_person_subset_pre_omit")
  expect_equal(as.integer(value_of(manifest$model_settings, "population_omitted_persons")), 1L)
  expect_equal(as.integer(value_of(manifest$model_settings, "population_response_rows_omitted")), 6L)
  expect_identical(value_of(manifest$settings, "fit_population_person_id"), "Person")
  expect_identical(value_of(manifest$settings, "fit_population_policy"), "omit")
  expect_equal(as.integer(value_of(manifest$settings, "fit_population_person_rows_replay")), 60L)
  expect_identical(value_of(manifest$settings, "fit_population_person_replay_scope"), "observed_person_subset_pre_omit")
  expect_equal(as.integer(value_of(manifest$settings, "fit_population_omitted_persons")), 1L)
  expect_equal(as.integer(value_of(manifest$settings, "fit_population_response_rows_omitted")), 6L)
  expect_identical(value_of(manifest$settings, "unit_prediction_posterior_basis"), "population_model")
  expect_identical(value_of(manifest$settings, "unit_prediction_person_id"), "Person")
  expect_identical(value_of(manifest$settings, "unit_prediction_population_policy"), "omit")
  expect_match(value_of(manifest$settings, "unit_prediction_population_formula"), "~\\s*X")
  expect_identical(value_of(manifest$settings, "plausible_value_posterior_basis"), "population_model")
  expect_identical(value_of(manifest$settings, "plausible_value_person_id"), "Person")
  expect_identical(value_of(manifest$settings, "plausible_value_population_policy"), "omit")
  expect_match(value_of(manifest$settings, "plausible_value_population_formula"), "~\\s*X")
})

test_that("build_mfrm_manifest returns bounded GPCM caveat boundary", {
  dat <- load_mfrmr_data("example_core")
  keep_people <- unique(dat$Person)[1:14]
  dat <- dat[dat$Person %in% keep_people, , drop = FALSE]
  fit_gpcm <- suppressWarnings(
    fit_mfrm(
      dat,
      "Person",
      c("Rater", "Criterion"),
      "Score",
      method = "MML",
      model = "GPCM",
      step_facet = "Criterion",
      quad_points = 5,
      maxit = 20
    )
  )

  manifest <- build_mfrm_manifest(fit_gpcm)

  expect_s3_class(manifest, "mfrm_manifest")
  expect_true(is.data.frame(manifest$gpcm_boundary))
  expect_true(any(manifest$gpcm_boundary$Area == "APA writer and fit-based export bundles"))
  expect_true(any(manifest$gpcm_boundary$Status == "supported_with_caveat"))
})

test_that("build_mfrm_replay_script reproduces optional prediction artifacts", {
  replay <- build_mfrm_replay_script(
    fit = prediction_bundle_fixture$fit,
    diagnostics = prediction_bundle_fixture$diagnostics,
    population_prediction = prediction_bundle_fixture$population_prediction,
    unit_prediction = prediction_bundle_fixture$unit_prediction,
    plausible_values = prediction_bundle_fixture$plausible_values,
    include_bundle = TRUE,
    bundle_prefix = "bundle_pred_test",
    data_file = "analysis_data.csv"
  )

  expect_s3_class(replay, "mfrm_replay_script")
  expect_match(replay$script, "predict_mfrm_population\\(")
  expect_match(replay$script, "population_prediction_sim_spec <- build_mfrm_sim_spec", fixed = TRUE)
  expect_false(grepl("planning_schema", replay$script, fixed = TRUE))
  expect_match(replay$script, "predict_mfrm_units\\(")
  expect_match(replay$script, "sample_mfrm_plausible_values\\(")
  expect_match(
    replay$script,
    'include = c\\("core_tables", "checklist", "dashboard", "manifest", "html",\\s*"predictions"\\)'
  )
  expect_true(replay$summary$PopulationPrediction[[1]])
  expect_true(replay$summary$UnitPrediction[[1]])
  expect_true(replay$summary$PlausibleValues[[1]])
  expect_identical(as.character(replay$summary$FitPosteriorBasis[[1]]), "legacy_mml")
})

test_that("build_mfrm_replay_script preserves latent-regression scoring inputs", {
  replay <- build_mfrm_replay_script(
    fit = latent_prediction_bundle_fixture$fit,
    unit_prediction = latent_prediction_bundle_fixture$unit_prediction,
    plausible_values = latent_prediction_bundle_fixture$plausible_values,
    data_file = "analysis_data.csv"
  )

  expect_s3_class(replay, "mfrm_replay_script")
  expect_match(replay$script, "fit_person_data <-")
  expect_match(replay$script, "population_formula = ~X", fixed = TRUE)
  expect_match(replay$script, "person_data = fit_person_data", fixed = TRUE)
  expect_match(replay$script, 'person_id = "Person"', fixed = TRUE)
  expect_match(replay$script, "unit_prediction_person_data <-")
  expect_match(replay$script, "plausible_value_person_data <-")
  expect_match(replay$script, "person_data = unit_prediction_person_data", fixed = TRUE)
  expect_match(replay$script, "person_data = plausible_value_person_data", fixed = TRUE)
  expect_match(replay$script, "# population_person_id = Person", fixed = TRUE)
  expect_equal(length(gregexpr('person_id = "Person"', replay$script, fixed = TRUE)[[1]]), 3L)
  expect_equal(length(gregexpr('population_policy = "error"', replay$script, fixed = TRUE)[[1]]), 3L)
  expect_true(isTRUE(replay$summary$FitPopulationActive[[1]]))
  expect_identical(as.character(replay$summary$FitPosteriorBasis[[1]]), "population_model")
  expect_identical(as.character(replay$summary$ScriptMode[[1]]), "fit")
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_population_active"][1]),
    "TRUE"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_posterior_basis"][1]),
    "population_model"
  )
  expect_match(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_population_formula"][1]),
    "~\\s*X"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_population_person_id"][1]),
    "Person"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_population_policy"][1]),
    "error"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_population_person_data_mode"][1]),
    "inline_literal"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "unit_prediction_posterior_basis"][1]),
    "population_model"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "unit_prediction_person_id"][1]),
    "Person"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "unit_prediction_population_policy"][1]),
    "error"
  )
  expect_match(
    as.character(replay$settings$Value[replay$settings$Setting == "unit_prediction_population_formula"][1]),
    "~\\s*X"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "plausible_value_posterior_basis"][1]),
    "population_model"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "plausible_value_person_id"][1]),
    "Person"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "plausible_value_population_policy"][1]),
    "error"
  )
  expect_match(
    as.character(replay$settings$Value[replay$settings$Setting == "plausible_value_population_formula"][1]),
    "~\\s*X"
  )
  expect_false(grepl("predict_mfrm_population\\(", replay$script))
})

test_that("build_mfrm_replay_script preserves replay-ready omit latent-regression inputs", {
  replay <- build_mfrm_replay_script(
    fit = latent_prediction_omit_bundle_fixture$fit,
    unit_prediction = latent_prediction_omit_bundle_fixture$unit_prediction,
    plausible_values = latent_prediction_omit_bundle_fixture$plausible_values,
    data_file = "analysis_data.csv"
  )

  expect_s3_class(replay, "mfrm_replay_script")
  expect_match(replay$script, "fit_person_data <-")
  expect_match(replay$script, 'population_policy = "omit"', fixed = TRUE)
  expect_match(replay$script, "# population_person_replay_scope = observed_person_subset_pre_omit", fixed = TRUE)
  expect_true(isTRUE(replay$summary$FitPopulationActive[[1]]))
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_population_person_rows_replay"][1]),
    "60"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_population_person_data_mode"][1]),
    "inline_literal"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_population_person_replay_scope"][1]),
    "observed_person_subset_pre_omit"
  )
})

test_that("build_mfrm_replay_script can externalize fit-level latent replay person data", {
  replay <- build_mfrm_replay_script(
    fit = latent_prediction_omit_bundle_fixture$fit,
    unit_prediction = latent_prediction_omit_bundle_fixture$unit_prediction,
    plausible_values = latent_prediction_omit_bundle_fixture$plausible_values,
    data_file = "analysis_data.csv",
    fit_person_data_file = "latent_fit_person_data.csv"
  )

  expect_s3_class(replay, "mfrm_replay_script")
  expect_match(replay$script, "replay_script_args <- commandArgs", fixed = TRUE)
  expect_match(replay$script, "sys.frames()[[1]]$ofile", fixed = TRUE)
  expect_match(replay$script, "latent_fit_person_data.csv", fixed = TRUE)
  expect_match(replay$script, "fit_person_data <- utils::read.csv", fixed = TRUE)
  expect_false(grepl("fit_person_data <- structure\\(", replay$script))
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_population_person_data_mode"][1]),
    "sidecar_csv"
  )
  expect_identical(
    as.character(replay$settings$Value[replay$settings$Setting == "fit_population_person_data_file"][1]),
    "latent_fit_person_data.csv"
  )
})

test_that("build_mfrm_replay_script rejects omit latent-regression fits without replay-ready person data", {
  fit <- latent_prediction_omit_bundle_fixture$fit
  fit$population$person_table_replay <- NULL
  fit$population$person_table_replay_scope <- NULL

  expect_error(
    build_mfrm_replay_script(
      fit = fit,
      unit_prediction = latent_prediction_omit_bundle_fixture$unit_prediction,
      plausible_values = latent_prediction_omit_bundle_fixture$plausible_values,
      data_file = "analysis_data.csv"
    ),
    "stored replay-ready person table no longer covers every observed person",
    fixed = TRUE
  )
})

test_that("build_conquest_overlap_bundle returns a minimal exact-overlap bundle", {
  bundle <- build_conquest_overlap_bundle()

  expect_s3_class(bundle, "mfrm_conquest_overlap_bundle")
  expect_true(is.data.frame(bundle$summary))
  expect_true(is.data.frame(bundle$comparison_targets))
  expect_true(is.data.frame(bundle$conquest_output_contract))
  expect_true(is.data.frame(bundle$response_wide))
  expect_true(is.data.frame(bundle$person_data))
  expect_true(is.data.frame(bundle$item_map))
  expect_true(is.data.frame(bundle$mfrmr_population))
  expect_true(is.data.frame(bundle$mfrmr_item_estimates))
  expect_true(is.data.frame(bundle$mfrmr_case_eap))
  expect_equal(bundle$summary$Case[[1]], "synthetic_latent_regression")
  expect_equal(bundle$summary$Persons[[1]], 60)
  expect_equal(bundle$summary$Items[[1]], 6)
  expect_identical(as.character(bundle$summary$PosteriorBasis[[1]]), "population_model")
  expect_identical(as.character(bundle$summary$PopulationDesignColumns[[1]]), "(Intercept), X")
  expect_identical(bundle$summary$PopulationCoefficientCount[[1]], 2L)
  expect_true(is.finite(bundle$summary$PopulationResidualVariance[[1]]))
  expect_identical(bundle$summary$PopulationIncludedPersons[[1]], 60L)
  expect_identical(bundle$summary$PopulationOmittedPersons[[1]], 0L)
  expect_identical(bundle$summary$PopulationResponseRowsOmitted[[1]], 0L)
  expect_true(all(c(
    "ExternalFile",
    "ConQuestCommand",
    "ReviewHandoff",
    "RequiredForReview"
  ) %in% names(bundle$conquest_output_contract)))
  expect_equal(sum(bundle$conquest_output_contract$RequiredForReview %in% TRUE), 4)
  expect_true(any(grepl("_conquest_parameters_review.txt", bundle$conquest_output_contract$ExternalFile, fixed = TRUE)))
  expect_true(all(sort(unique(unlist(bundle$response_wide[sprintf("I%03d", 1:6)]))) %in% c(0, 1)))
  expect_identical(names(bundle$person_data), c("Person", "X"))
  expect_true(is.numeric(bundle$person_data$X))
  expect_true(all(c("Person", "X", sprintf("I%03d", 1:6)) %in% names(bundle$response_wide)))
  expect_match(bundle$conquest_command, "filetype=csv", fixed = TRUE)
  expect_match(bundle$conquest_command, "/*", fixed = TRUE)
  expect_match(bundle$conquest_command, "*/", fixed = TRUE)
  expect_false(any(grepl("^\\*\\s", strsplit(bundle$conquest_command, "\n", fixed = TRUE)[[1]])))
  expect_match(bundle$conquest_command, "pidwidth=32", fixed = TRUE)
  expect_match(bundle$conquest_command, "keepswidth=32", fixed = TRUE)
  expect_match(bundle$conquest_command, "regression X;", fixed = TRUE)
  expect_match(bundle$conquest_command, "model item;", fixed = TRUE)
  expect_match(bundle$conquest_command, "export parameters ! filetype=csv", fixed = TRUE)
  expect_match(bundle$conquest_command, "export reg_coefficients ! filetype=csv", fixed = TRUE)
  expect_match(bundle$conquest_command, "export covariance ! filetype=csv", fixed = TRUE)
  expect_match(bundle$conquest_command, "show cases ! estimates=eap, filetype=csv, regressors=yes", fixed = TRUE)

  s <- summary(bundle)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_identical(as.character(s$overview$Class[1]), "mfrm_conquest_overlap_bundle")
  expect_true(is.data.frame(s$conquest_command_scope))
  expect_true(is.data.frame(s$conquest_output_contract))
  expect_true(all(c(
    "ConQuest command template",
    "Command-comment syntax",
    "Official command-reference alignment",
    "Overlap model scope",
    "External output requirements",
    "External comparison scope"
  ) %in% s$conquest_command_scope$Area))
  expect_identical(
    as.character(s$conquest_command_scope$Status[s$conquest_command_scope$Area == "ConQuest command template"][1]),
    "template only"
  )
  expect_identical(
    as.character(s$conquest_command_scope$Status[s$conquest_command_scope$Area == "Command-comment syntax"][1]),
    "block comments"
  )
  expect_identical(
    as.character(s$conquest_command_scope$Status[s$conquest_command_scope$Area == "External comparison scope"][1]),
    "not claimed"
  )
  expect_true(any(grepl("reg_coefficients", s$conquest_output_contract$ConQuestCommand, fixed = TRUE)))
  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "ConQuest command scope", fixed = TRUE)
  expect_match(printed, "ConQuest output contract", fixed = TRUE)
  expect_match(printed, "block comments", fixed = TRUE)
  expect_match(printed, "explicit CSV widths", fixed = TRUE)
  expect_match(printed, "not claimed", fixed = TRUE)
})

test_that("build_conquest_overlap_bundle rejects intercept-only latent-regression fits", {
  fixture <- mfrmr:::with_preserved_rng_seed(20260410, {
    persons <- paste0("P", sprintf("%02d", 1:40))
    items <- paste0("I", 1:5)
    theta <- 0.25 + stats::rnorm(length(persons), sd = 0.7)
    item_beta <- seq(-0.8, 0.8, length.out = length(items))

    dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
    eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
    dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))

    fit_mfrm(
      dat,
      "Person",
      "Item",
      "Score",
      method = "MML",
      model = "RSM",
      population_formula = ~ 1,
      person_data = data.frame(Person = persons, stringsAsFactors = FALSE),
      quad_points = 7,
      maxit = 60
    )
  })

  expect_error(
    build_conquest_overlap_bundle(fixture),
    "requires exactly one numeric person covariate beyond the intercept",
    fixed = TRUE
  )
})

test_that("build_conquest_overlap_bundle rejects unsupported overlap contracts", {
  base_fit <- latent_prediction_bundle_fixture$fit

  fit_jml <- export_core_fixture$fit
  expect_error(
    build_conquest_overlap_bundle(fit_jml),
    "supports only `MML` fits",
    fixed = TRUE
  )

  fit_gpcm <- base_fit
  fit_gpcm$config$model <- "GPCM"
  expect_error(
    build_conquest_overlap_bundle(fit_gpcm),
    "restricted to the ordered-response `RSM` / `PCM` model scope",
    fixed = TRUE
  )

  fit_inactive <- base_fit
  fit_inactive$population$active <- FALSE
  fit_inactive$config$posterior_basis <- "legacy_mml"
  expect_error(
    build_conquest_overlap_bundle(fit_inactive),
    "requires an active latent-regression `MML` fit",
    fixed = TRUE
  )

  fit_legacy_basis <- base_fit
  fit_legacy_basis$config$posterior_basis <- "legacy_mml"
  fit_legacy_basis$population$posterior_basis <- "legacy_mml"
  expect_error(
    build_conquest_overlap_bundle(fit_legacy_basis),
    "requires an active latent-regression `MML` fit",
    fixed = TRUE
  )

  fit_multifacet <- base_fit
  fit_multifacet$config$facet_names <- c("Item", "Task")
  expect_error(
    build_conquest_overlap_bundle(fit_multifacet),
    "exactly one non-person facet",
    fixed = TRUE
  )

  fit_weighted <- base_fit
  fit_weighted$prep$data$Weight <- 1
  fit_weighted$prep$data$Weight[1] <- 2
  expect_error(
    build_conquest_overlap_bundle(fit_weighted),
    "unit weights only",
    fixed = TRUE
  )

  fit_nonbinary <- base_fit
  fit_nonbinary$prep$data$Score[1] <- 2
  expect_error(
    build_conquest_overlap_bundle(fit_nonbinary),
    "requires binary responses",
    fixed = TRUE
  )

  fit_incomplete <- base_fit
  fit_incomplete$prep$data <- fit_incomplete$prep$data[-1, , drop = FALSE]
  expect_error(
    build_conquest_overlap_bundle(fit_incomplete),
    "complete rectangular person-by-item response matrix",
    fixed = TRUE
  )

  fit_nonnumeric_cov <- base_fit
  fit_nonnumeric_cov$population$person_table$X <- as.character(fit_nonnumeric_cov$population$person_table$X)
  expect_error(
    build_conquest_overlap_bundle(fit_nonnumeric_cov),
    "requires a numeric person covariate",
    fixed = TRUE
  )

  fit_two_covariates <- base_fit
  fit_two_covariates$population$person_table$Z <- seq_len(nrow(fit_two_covariates$population$person_table))
  fit_two_covariates$population$design_columns <- c("(Intercept)", "X", "Z")
  expect_error(
    build_conquest_overlap_bundle(fit_two_covariates),
    "requires exactly one numeric person covariate beyond the intercept",
    fixed = TRUE
  )

  fit_missing_score <- base_fit
  fit_missing_score$prep$data$Score <- NULL
  expect_error(
    build_conquest_overlap_bundle(fit_missing_score),
    "canonical response columns",
    fixed = TRUE
  )

  fit_missing_cov <- base_fit
  fit_missing_cov$population$person_table$X <- NULL
  expect_error(
    build_conquest_overlap_bundle(fit_missing_cov),
    "required ConQuest-overlap covariate columns",
    fixed = TRUE
  )

  fit_categorical_cov <- base_fit
  fit_categorical_cov$population$person_table$Group <- rep(
    c("A", "B"),
    length.out = nrow(fit_categorical_cov$population$person_table)
  )
  fit_categorical_cov$population$person_table$X <- NULL
  fit_categorical_cov$population$design_columns <- c("(Intercept)", "GroupB")
  fit_categorical_cov$population$xlevels <- list(Group = c("A", "B"))
  fit_categorical_cov$population$contrasts <- list(Group = "contr.treatment")
  expect_error(
    build_conquest_overlap_bundle(fit_categorical_cov),
    "Categorical/model-matrix-expanded covariates",
    fixed = TRUE
  )

  fit_duplicate_cell <- base_fit
  fit_duplicate_cell$prep$data <- rbind(
    fit_duplicate_cell$prep$data,
    fit_duplicate_cell$prep$data[1, , drop = FALSE]
  )
  expect_error(
    build_conquest_overlap_bundle(fit_duplicate_cell),
    "exactly one response per person-by-item cell",
    fixed = TRUE
  )

  fit_bad_standardize <- base_fit
  fit_bad_standardize$prep$data$Score[1] <- Inf
  expect_error(
    build_conquest_overlap_bundle(fit_bad_standardize),
    "could not be standardized to {0, 1}",
    fixed = TRUE
  )

  fit_bad_alignment <- base_fit
  fit_bad_alignment$population$person_table$Person[1] <- "MISSING_PERSON"
  expect_error(
    build_conquest_overlap_bundle(fit_bad_alignment),
    "do not align with all fitted persons",
    fixed = TRUE
  )
})

test_that("build_conquest_overlap_bundle writes expected external-comparison files", {
  out_dir <- file.path(tempdir(), "mfrmr-conquest-overlap")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  bundle <- build_conquest_overlap_bundle(
    output_dir = out_dir,
    prefix = "cq_overlap_test",
    overwrite = TRUE
  )

  expect_s3_class(bundle, "mfrm_conquest_overlap_bundle")
  expect_true(is.data.frame(bundle$written_files))
  expect_true(any(bundle$written_files$Component == "response_wide"))
  expect_true(any(bundle$written_files$Component == "conquest_command"))
  expect_true(file.exists(file.path(out_dir, "cq_overlap_test_wide.csv")))
  expect_true(file.exists(file.path(out_dir, "cq_overlap_test.cqc")))
  expect_true(file.exists(file.path(out_dir, "cq_overlap_test_mfrmr_population.csv")))
  expect_true(file.exists(file.path(out_dir, "cq_overlap_test_mfrmr_item_estimates.csv")))
  expect_true(file.exists(file.path(out_dir, "cq_overlap_test_mfrmr_case_eap.csv")))
  expect_true(file.exists(file.path(out_dir, "cq_overlap_test_conquest_output_contract.csv")))
  expect_true(file.exists(file.path(out_dir, "cq_overlap_test_README.txt")))
  readme <- paste(readLines(file.path(out_dir, "cq_overlap_test_README.txt"), warn = FALSE), collapse = "\n")
  expect_match(readme, "PopulationDesignColumns", fixed = TRUE)
  expect_match(readme, "PopulationResponseRowsOmitted", fixed = TRUE)
  expect_match(readme, "Requested external ConQuest outputs", fixed = TRUE)
  expect_match(readme, "cq_overlap_test_conquest_reg_coefficients.csv", fixed = TRUE)
})

test_that("build_conquest_overlap_bundle accepts the documented PCM overlap surface", {
  fixture <- mfrmr:::with_preserved_rng_seed(20260412, {
    persons <- paste0("P", sprintf("%02d", 1:24))
    items <- paste0("I", 1:4)
    x <- stats::rnorm(length(persons))
    theta <- 0.2 + 0.7 * x + stats::rnorm(length(persons), sd = 0.65)
    item_beta <- seq(-0.9, 0.9, length.out = length(items))
    dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
    eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
    dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))
    person_tbl <- data.frame(Person = persons, X = x, stringsAsFactors = FALSE)

    fit <- suppressWarnings(fit_mfrm(
      dat,
      "Person",
      "Item",
      "Score",
      method = "MML",
      model = "PCM",
      step_facet = "Item",
      population_formula = ~ X,
      person_data = person_tbl,
      quad_points = 5,
      maxit = 40
    ))

    list(fit = fit, n_person = length(persons))
  })

  bundle <- build_conquest_overlap_bundle(fixture$fit)

  expect_s3_class(bundle, "mfrm_conquest_overlap_bundle")
  expect_identical(as.character(bundle$summary$Model[1]), "PCM")
  expect_identical(as.character(bundle$summary$PopulationDesignColumns[1]), "(Intercept), X")
  expect_identical(bundle$summary$PopulationOmittedPersons[1], 0L)
  expect_equal(nrow(bundle$response_wide), fixture$n_person)
  expect_equal(nrow(bundle$person_data), fixture$n_person)

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      ItemID = bundle$mfrmr_item_estimates$ResponseVar,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemID",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 0)
  expect_true(all(abs(audit$item_comparison$CenteredDifference[audit$item_comparison$Status == "Compared"]) < 1e-10))
})

test_that("review_conquest_overlap compares normalized ConQuest tables", {
  bundle <- build_conquest_overlap_bundle()

  cq_pop <- data.frame(
    Term = bundle$mfrmr_population$Parameter,
    Est = bundle$mfrmr_population$Estimate,
    stringsAsFactors = FALSE
  )
  cq_item <- data.frame(
    ItemID = bundle$mfrmr_item_estimates$ResponseVar,
    Est = bundle$mfrmr_item_estimates$Estimate,
    stringsAsFactors = FALSE
  )
  cq_case <- data.frame(
    PID = bundle$mfrmr_case_eap$Person,
    EAP = bundle$mfrmr_case_eap$Estimate,
    stringsAsFactors = FALSE
  )

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = cq_pop,
    conquest_item_estimates = cq_item,
    conquest_case_eap = cq_case,
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemID",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_true(is.data.frame(audit$overall))
  expect_true(is.data.frame(audit$population_comparison))
  expect_true(is.data.frame(audit$item_comparison))
  expect_true(is.data.frame(audit$case_comparison))
  expect_equal(audit$overall$AttentionItems[[1]], 0)
  expect_equal(audit$overall$AttentionMissing[[1]], 0)
  expect_equal(audit$overall$AttentionDuplicate[[1]], 0)
  expect_equal(audit$overall$AttentionNonNumeric[[1]], 0)
  expect_true(all(abs(audit$population_comparison$Difference[audit$population_comparison$Status == "Compared"]) < 1e-10))
  expect_true(all(abs(audit$item_comparison$CenteredDifference[audit$item_comparison$Status == "Compared"]) < 1e-10))
  expect_true(all(abs(audit$case_comparison$Difference[audit$case_comparison$Status == "Compared"]) < 1e-10))
  expect_equal(audit$overall$PopulationMaxAbsDifference[[1]], 0, tolerance = 1e-10)
  expect_equal(audit$overall$ItemCenteredMaxAbsDifference[[1]], 0, tolerance = 1e-10)
  expect_equal(audit$overall$CaseMaxAbsDifference[[1]], 0, tolerance = 1e-10)
  expect_true(all(c(
    "PopulationMaxAbsParameter",
    "ItemCenteredMaxAbsItem",
    "CaseMaxAbsPerson"
  ) %in% names(audit$overall)))

  s <- summary(audit)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_identical(as.character(s$overview$Class[1]), "mfrm_conquest_overlap_review")
  expect_true(is.data.frame(s$review_scope))
  expect_true(is.data.frame(s$review_scope))
  expect_true(all(c(
    "User-supplied table review",
    "Raw ConQuest text parsing",
    "External comparison scope",
    "Attention items"
  ) %in% s$review_scope$Area))
  expect_true(any(grepl("normalized ConQuest tables", s$notes, fixed = TRUE)))
  attention_row <- s$review_scope[
    s$review_scope$Area == "Attention items",
    ,
    drop = FALSE
  ]
  expect_identical(as.character(attention_row$Status[1]), "none detected")
  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "Review scope", fixed = TRUE)
  expect_match(printed, "External comparison scope", fixed = TRUE)
  expect_match(printed, "not claimed", fixed = TRUE)
})

test_that("review_conquest_overlap separates constraint shifts from direct differences", {
  bundle <- build_conquest_overlap_bundle()

  pop_delta <- ifelse(
    bundle$mfrmr_population$Parameter == "(Intercept)",
    2.50,
    ifelse(bundle$mfrmr_population$Parameter == "sigma2", -0.03, 0.02)
  )
  cq_pop <- data.frame(
    Term = bundle$mfrmr_population$Parameter,
    Est = bundle$mfrmr_population$Estimate + pop_delta,
    stringsAsFactors = FALSE
  )
  cq_item <- data.frame(
    ItemID = bundle$mfrmr_item_estimates$ResponseVar,
    Est = bundle$mfrmr_item_estimates$Estimate + 1.75,
    stringsAsFactors = FALSE
  )
  cq_case <- data.frame(
    PID = bundle$mfrmr_case_eap$Person,
    EAP = bundle$mfrmr_case_eap$Estimate + 0.125,
    stringsAsFactors = FALSE
  )

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = cq_pop,
    conquest_item_estimates = cq_item,
    conquest_case_eap = cq_case,
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemID",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 0)
  expect_equal(audit$overall$AttentionMissing[[1]], 0)
  expect_equal(audit$overall$AttentionDuplicate[[1]], 0)
  expect_equal(audit$overall$AttentionNonNumeric[[1]], 0)
  expect_identical(
    audit$population_comparison$Status[audit$population_comparison$Parameter == "(Intercept)"][1],
    "ConstraintDependent"
  )
  direct_pop <- audit$population_comparison[audit$population_comparison$Status == "Compared", , drop = FALSE]
  expect_equal(sort(abs(direct_pop$Difference)), c(0.02, 0.03), tolerance = 1e-10)
  expect_equal(audit$overall$PopulationMae[[1]], mean(c(0.02, 0.03)), tolerance = 1e-10)
  expect_equal(audit$overall$PopulationMaxAbsDifference[[1]], 0.03, tolerance = 1e-10)
  expect_identical(audit$overall$PopulationMaxAbsParameter[[1]], "sigma2")
  expect_true(all(abs(audit$item_comparison$CenteredDifference) < 1e-10))
  expect_equal(audit$overall$ItemCenteredMae[[1]], 0, tolerance = 1e-10)
  expect_equal(audit$overall$ItemCenteredMaxAbsDifference[[1]], 0, tolerance = 1e-10)
  expect_true(all(abs(audit$case_comparison$Difference - 0.125) < 1e-10))
  expect_equal(audit$overall$CaseMae[[1]], 0.125, tolerance = 1e-10)
  expect_equal(audit$overall$CaseMaxAbsDifference[[1]], 0.125, tolerance = 1e-10)
})

test_that("review_conquest_overlap records worst compared rows in overall", {
  bundle <- build_conquest_overlap_bundle()
  skip_if(!("sigma2" %in% bundle$mfrmr_population$Parameter))
  skip_if(nrow(bundle$mfrmr_item_estimates) <= 2L)
  skip_if(nrow(bundle$mfrmr_case_eap) <= 3L)

  worst_item <- as.character(bundle$mfrmr_item_estimates$ResponseVar[2])
  worst_person <- as.character(bundle$mfrmr_case_eap$Person[3])

  cq_pop <- data.frame(
    Term = bundle$mfrmr_population$Parameter,
    Est = bundle$mfrmr_population$Estimate,
    stringsAsFactors = FALSE
  )
  cq_pop$Est[cq_pop$Term == "sigma2"] <- cq_pop$Est[cq_pop$Term == "sigma2"] - 0.4
  cq_item <- data.frame(
    ItemID = bundle$mfrmr_item_estimates$ResponseVar,
    Est = bundle$mfrmr_item_estimates$Estimate,
    stringsAsFactors = FALSE
  )
  cq_item$Est[cq_item$ItemID == worst_item] <- cq_item$Est[cq_item$ItemID == worst_item] + 0.4
  cq_case <- data.frame(
    PID = bundle$mfrmr_case_eap$Person,
    EAP = bundle$mfrmr_case_eap$Estimate,
    stringsAsFactors = FALSE
  )
  cq_case$EAP[cq_case$PID == worst_person] <- cq_case$EAP[cq_case$PID == worst_person] - 0.5

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = cq_pop,
    conquest_item_estimates = cq_item,
    conquest_case_eap = cq_case,
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemID",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_identical(audit$overall$PopulationMaxAbsParameter[[1]], "sigma2")
  expect_identical(audit$overall$ItemCenteredMaxAbsItem[[1]], worst_item)
  expect_identical(audit$overall$CaseMaxAbsPerson[[1]], worst_person)
  expect_gt(audit$overall$PopulationMaxAbsDifference[[1]], audit$overall$PopulationMae[[1]])
  expect_gt(audit$overall$ItemCenteredMaxAbsDifference[[1]], audit$overall$ItemCenteredMae[[1]])
  expect_gt(audit$overall$CaseMaxAbsDifference[[1]], audit$overall$CaseMae[[1]])
})

test_that("review_conquest_overlap records non-numeric extracted estimates", {
  bundle <- build_conquest_overlap_bundle()

  nonnumeric_parameter <- bundle$mfrmr_population$Parameter[bundle$mfrmr_population$Parameter != "(Intercept)"][1]
  nonnumeric_item <- bundle$mfrmr_item_estimates$ResponseVar[2]
  nonnumeric_person <- bundle$mfrmr_case_eap$Person[3]

  cq_pop <- data.frame(
    Term = bundle$mfrmr_population$Parameter,
    Est = as.character(bundle$mfrmr_population$Estimate),
    stringsAsFactors = FALSE
  )
  cq_pop$Est[cq_pop$Term == nonnumeric_parameter] <- "not_reported"
  cq_item <- data.frame(
    ItemID = bundle$mfrmr_item_estimates$ResponseVar,
    Est = as.character(bundle$mfrmr_item_estimates$Estimate),
    stringsAsFactors = FALSE
  )
  cq_item$Est[cq_item$ItemID == nonnumeric_item] <- "****"
  cq_case <- data.frame(
    PID = bundle$mfrmr_case_eap$Person,
    EAP = as.character(bundle$mfrmr_case_eap$Estimate),
    stringsAsFactors = FALSE
  )
  cq_case$EAP[cq_case$PID == nonnumeric_person] <- "omitted"

  raw_audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = cq_pop,
    conquest_item_estimates = cq_item,
    conquest_case_eap = cq_case,
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemID",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_equal(raw_audit$overall$AttentionItems[[1]], 3)
  expect_equal(raw_audit$overall$AttentionMissing[[1]], 0)
  expect_equal(raw_audit$overall$AttentionDuplicate[[1]], 0)
  expect_equal(raw_audit$overall$AttentionNonNumeric[[1]], 3)
  expect_true(all(c(
    "non_numeric_conquest_parameter",
    "non_numeric_conquest_item",
    "non_numeric_conquest_case"
  ) %in% raw_audit$attention_items$Issue))

  normalized <- normalize_conquest_overlap_tables(
    conquest_population = cq_pop,
    conquest_item_estimates = cq_item,
    conquest_case_eap = cq_case,
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemID",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_equal(normalized$summary$PopulationNonNumeric[[1]], 1)
  expect_equal(normalized$summary$ItemNonNumeric[[1]], 1)
  expect_equal(normalized$summary$CaseNonNumeric[[1]], 1)
  expect_true(all(c("EstimateNonNumeric") %in% names(normalized$conquest_population)))
  expect_true(all(c("EstimateNonNumeric") %in% names(normalized$conquest_item_estimates)))
  expect_true(all(c("EstimateNonNumeric") %in% names(normalized$conquest_case_eap)))
  ns <- summary(normalized)
  expect_s3_class(ns, "summary.mfrm_bundle")
  expect_true(is.data.frame(ns$normalization_scope))
  review_row <- ns$normalization_scope[
    ns$normalization_scope$Area == "Pre-review table check",
    ,
    drop = FALSE
  ]
  expect_identical(as.character(review_row$Status[1]), "review required")
  expect_match(as.character(review_row$Evidence[1]), "3 non-numeric estimate cell", fixed = TRUE)

  audit <- review_conquest_overlap(bundle, normalized)

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 3)
  expect_equal(audit$overall$AttentionMissing[[1]], 0)
  expect_equal(audit$overall$AttentionDuplicate[[1]], 0)
  expect_equal(audit$overall$AttentionNonNumeric[[1]], 3)
  expect_true(all(c(
    "non_numeric_conquest_parameter",
    "non_numeric_conquest_item",
    "non_numeric_conquest_case"
  ) %in% audit$attention_items$Issue))
  expect_identical(
    audit$population_comparison$Status[audit$population_comparison$Parameter == nonnumeric_parameter][1],
    "NonNumericInConQuest"
  )
  expect_identical(
    audit$item_comparison$Status[audit$item_comparison$MatchID == nonnumeric_item][1],
    "NonNumericInConQuest"
  )
  expect_identical(
    audit$case_comparison$Status[audit$case_comparison$Person == nonnumeric_person][1],
    "NonNumericInConQuest"
  )
  s <- summary(audit)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_true(all(c("AttentionMissing", "AttentionDuplicate", "AttentionNonNumeric") %in% names(s$summary)))
  expect_true(all(c(
    "PopulationMaxAbsDifference",
    "ItemCenteredMaxAbsDifference",
    "CaseMaxAbsDifference",
    "PopulationMaxAbsParameter",
    "ItemCenteredMaxAbsItem",
    "CaseMaxAbsPerson"
  ) %in% names(s$summary)))
  expect_equal(s$summary$AttentionNonNumeric[[1]], 3)
  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "AttentionNonNumeric", fixed = TRUE)
  expect_match(printed, "CaseMaxAbsDifference", fixed = TRUE)
})

test_that("normalize_conquest_overlap_tables standardizes extracted tables", {
  bundle <- build_conquest_overlap_bundle()

  normalized <- normalize_conquest_overlap_tables(
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      Group = "population",
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      ItemCode = bundle$mfrmr_item_estimates$ResponseVar,
      Est = bundle$mfrmr_item_estimates$Estimate,
      Source = "items",
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      Batch = "cases",
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemCode",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(normalized, "mfrm_conquest_overlap_tables")
  expect_true(all(c("Parameter", "Estimate", "Group") %in% names(normalized$conquest_population)))
  expect_true(all(c("ItemID", "Estimate", "Source") %in% names(normalized$conquest_item_estimates)))
  expect_true(all(c("Person", "Estimate", "Batch") %in% names(normalized$conquest_case_eap)))
  expect_true(is.numeric(normalized$conquest_population$Estimate))
  expect_true(is.numeric(normalized$conquest_item_estimates$Estimate))
  expect_true(is.numeric(normalized$conquest_case_eap$Estimate))
  expect_equal(normalized$summary$PopulationRows[[1]], nrow(bundle$mfrmr_population))
  expect_equal(normalized$summary$ItemRows[[1]], nrow(bundle$mfrmr_item_estimates))
  expect_equal(normalized$summary$CaseRows[[1]], nrow(bundle$mfrmr_case_eap))

  s <- summary(normalized)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_identical(as.character(s$overview$Class[1]), "mfrm_conquest_overlap_tables")
  expect_true(is.data.frame(s$normalization_scope))
  expect_true(all(c(
    "Extracted table normalization",
    "Raw ConQuest text parsing",
    "Bundle matching",
    "Pre-review table check"
  ) %in% s$normalization_scope$Area))
  review_row <- s$normalization_scope[
    s$normalization_scope$Area == "Pre-review table check",
    ,
    drop = FALSE
  ]
  expect_identical(as.character(review_row$Status[1]), "none detected")
  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "Normalization scope", fixed = TRUE)
  expect_match(printed, "Raw ConQuest text parsing", fixed = TRUE)
  expect_match(printed, "deferred to review", fixed = TRUE)
})

test_that("normalize_conquest_overlap_tables rejects unresolved automatic aliases", {
  expect_error(
    normalize_conquest_overlap_tables(
      conquest_population = data.frame(
        UnknownTerm = c("(Intercept)", "X"),
        UnknownEstimate = c(0.1, 0.2),
        stringsAsFactors = FALSE
      ),
      conquest_item_estimates = data.frame(
        ItemID = c("I001", "I002"),
        Estimate = c(-0.1, 0.1),
        stringsAsFactors = FALSE
      ),
      conquest_case_eap = data.frame(
        Person = c("P01", "P02"),
        Estimate = c(-0.2, 0.2),
        stringsAsFactors = FALSE
      )
    ),
    "`conquest_population_term` could not be resolved automatically",
    fixed = TRUE
  )
})

test_that("normalize_conquest_overlap_tables rejects ambiguous automatic aliases", {
  expect_error(
    normalize_conquest_overlap_tables(
      conquest_population = data.frame(
        Parameter = c("(Intercept)", "X"),
        Term = c("(Intercept)", "X"),
        Estimate = c(0.1, 0.2),
        stringsAsFactors = FALSE
      ),
      conquest_item_estimates = data.frame(
        ItemID = c("I001", "I002"),
        Estimate = c(-0.1, 0.1),
        stringsAsFactors = FALSE
      ),
      conquest_case_eap = data.frame(
        Person = c("P01", "P02"),
        Estimate = c(-0.2, 0.2),
        stringsAsFactors = FALSE
      )
    ),
    "`conquest_population_term` matched multiple columns",
    fixed = TRUE
  )
})

test_that("normalize_conquest_overlap_tables rejects missing explicit columns", {
  expect_error(
    normalize_conquest_overlap_tables(
      conquest_population = data.frame(
        Term = c("(Intercept)", "X"),
        Estimate = c(0.1, 0.2),
        stringsAsFactors = FALSE
      ),
      conquest_item_estimates = data.frame(
        ItemID = c("I001", "I002"),
        Estimate = c(-0.1, 0.1),
        stringsAsFactors = FALSE
      ),
      conquest_case_eap = data.frame(
        Person = c("P01", "P02"),
        Estimate = c(-0.2, 0.2),
        stringsAsFactors = FALSE
      ),
      conquest_item_id = "MissingItemID"
    ),
    "`conquest_item_id` must name an existing column in `conquest_item_estimates`.",
    fixed = TRUE
  )
})

test_that("review_conquest_overlap accepts normalized contract objects", {
  bundle <- build_conquest_overlap_bundle()

  normalized <- normalize_conquest_overlap_tables(
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      ItemCode = bundle$mfrmr_item_estimates$ResponseVar,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemCode",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  audit <- review_conquest_overlap(bundle, normalized)

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 0)
  expect_true(all(abs(audit$population_comparison$Difference[audit$population_comparison$Status == "Compared"]) < 1e-10))
  expect_true(all(abs(audit$item_comparison$CenteredDifference[audit$item_comparison$Status == "Compared"]) < 1e-10))
  expect_true(all(abs(audit$case_comparison$Difference[audit$case_comparison$Status == "Compared"]) < 1e-10))
  expect_identical(as.character(audit$settings$Value[audit$settings$Setting == "conquest_item_id"][1]), "ItemID")
})

test_that("review_conquest_overlap rejects mixed normalized and raw external tables", {
  bundle <- build_conquest_overlap_bundle()

  normalized <- normalize_conquest_overlap_tables(
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      ItemCode = bundle$mfrmr_item_estimates$ResponseVar,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemCode",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_error(
    review_conquest_overlap(
      bundle,
      normalized,
      conquest_item_estimates = data.frame(ItemID = character(0), Estimate = numeric(0))
    ),
    "do not also supply `conquest_item_estimates` or `conquest_case_eap` separately",
    fixed = TRUE
  )
})

test_that("review_conquest_overlap requires either a normalized object or all three raw tables", {
  bundle <- build_conquest_overlap_bundle()

  expect_error(
    review_conquest_overlap(
      bundle,
      conquest_population = data.frame(
        Term = bundle$mfrmr_population$Parameter,
        Est = bundle$mfrmr_population$Estimate,
        stringsAsFactors = FALSE
      ),
      conquest_item_estimates = data.frame(
        ItemID = bundle$mfrmr_item_estimates$ResponseVar,
        Est = bundle$mfrmr_item_estimates$Estimate,
        stringsAsFactors = FALSE
      )
    ),
    "Supply either three extracted ConQuest tables or a single object from normalize_conquest_overlap_tables().",
    fixed = TRUE
  )
})

test_that("review_conquest_overlap records duplicate external rows as attention items", {
  bundle <- build_conquest_overlap_bundle()

  dup_pop <- rbind(
    data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Term = bundle$mfrmr_population$Parameter[1],
      Est = bundle$mfrmr_population$Estimate[1],
      stringsAsFactors = FALSE
    )
  )
  dup_item <- rbind(
    data.frame(
      ItemID = bundle$mfrmr_item_estimates$ResponseVar,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    data.frame(
      ItemID = bundle$mfrmr_item_estimates$ResponseVar[1],
      Est = bundle$mfrmr_item_estimates$Estimate[1],
      stringsAsFactors = FALSE
    )
  )
  dup_case <- rbind(
    data.frame(
      PersonID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    ),
    data.frame(
      PersonID = bundle$mfrmr_case_eap$Person[1],
      EAP = bundle$mfrmr_case_eap$Estimate[1],
      stringsAsFactors = FALSE
    )
  )

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = dup_pop,
    conquest_item_estimates = dup_item,
    conquest_case_eap = dup_case,
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemID",
    conquest_item_estimate = "Est",
    conquest_case_person = "PersonID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 3)
  expect_equal(audit$overall$AttentionMissing[[1]], 0)
  expect_equal(audit$overall$AttentionDuplicate[[1]], 3)
  expect_equal(audit$overall$AttentionNonNumeric[[1]], 0)
  expect_true(all(c(
    "duplicate_conquest_parameter",
    "duplicate_conquest_item",
    "duplicate_conquest_case"
  ) %in% audit$attention_items$Issue))
  expect_true(all(abs(audit$population_comparison$Difference[audit$population_comparison$Status == "Compared"]) < 1e-10))
  expect_true(all(abs(audit$item_comparison$CenteredDifference[audit$item_comparison$Status == "Compared"]) < 1e-10))
  expect_true(all(abs(audit$case_comparison$Difference[audit$case_comparison$Status == "Compared"]) < 1e-10))
})

test_that("review_conquest_overlap records missing external rows as attention items", {
  bundle <- build_conquest_overlap_bundle()

  missing_parameter <- bundle$mfrmr_population$Parameter[bundle$mfrmr_population$Parameter != "(Intercept)"][1]
  missing_item <- bundle$mfrmr_item_estimates$ResponseVar[1]
  missing_person <- bundle$mfrmr_case_eap$Person[1]

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter[bundle$mfrmr_population$Parameter != missing_parameter],
      Est = bundle$mfrmr_population$Estimate[bundle$mfrmr_population$Parameter != missing_parameter],
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      ItemID = bundle$mfrmr_item_estimates$ResponseVar[bundle$mfrmr_item_estimates$ResponseVar != missing_item],
      Est = bundle$mfrmr_item_estimates$Estimate[bundle$mfrmr_item_estimates$ResponseVar != missing_item],
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PersonID = bundle$mfrmr_case_eap$Person[bundle$mfrmr_case_eap$Person != missing_person],
      EAP = bundle$mfrmr_case_eap$Estimate[bundle$mfrmr_case_eap$Person != missing_person],
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemID",
    conquest_item_estimate = "Est",
    conquest_case_person = "PersonID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 3)
  expect_equal(audit$overall$AttentionMissing[[1]], 3)
  expect_equal(audit$overall$AttentionDuplicate[[1]], 0)
  expect_equal(audit$overall$AttentionNonNumeric[[1]], 0)
  expect_true(all(c(
    "missing_conquest_parameter",
    "missing_conquest_item",
    "missing_conquest_case"
  ) %in% audit$attention_items$Issue))
  expect_identical(
    audit$population_comparison$Status[audit$population_comparison$Parameter == missing_parameter][1],
    "MissingInConQuest"
  )
  expect_identical(
    audit$item_comparison$Status[audit$item_comparison$MatchID == missing_item][1],
    "MissingInConQuest"
  )
  expect_identical(
    audit$case_comparison$Status[audit$case_comparison$Person == missing_person][1],
    "MissingInConQuest"
  )

  s <- summary(audit)
  attention_row <- s$review_scope[
    s$review_scope$Area == "Attention items",
    ,
    drop = FALSE
  ]
  expect_identical(as.character(attention_row$Status[1]), "review required")
  expect_identical(as.character(attention_row$Evidence[1]), "3 attention item(s)")
  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "review required", fixed = TRUE)
  expect_match(printed, "missing_conquest_parameter", fixed = TRUE)
})

test_that("review_conquest_overlap matches item IDs by response variable when requested", {
  bundle <- build_conquest_overlap_bundle()

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      Label = bundle$mfrmr_item_estimates$ResponseVar,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "Label",
    conquest_item_estimate = "Est",
    item_id_source = "response_var",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 0)
  expect_true(all(abs(audit$item_comparison$CenteredDifference[audit$item_comparison$Status == "Compared"]) < 1e-10))
  expect_identical(as.character(audit$settings$Value[audit$settings$Setting == "item_id_source"][1]), "response_var")
})

test_that("review_conquest_overlap auto-detects response-variable item IDs", {
  bundle <- build_conquest_overlap_bundle()

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      Label = bundle$mfrmr_item_estimates$ResponseVar,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "Label",
    conquest_item_estimate = "Est",
    item_id_source = "auto",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 0)
  expect_true(all(abs(audit$item_comparison$CenteredDifference[audit$item_comparison$Status == "Compared"]) < 1e-10))
  expect_identical(as.character(audit$settings$Value[audit$settings$Setting == "item_id_source"][1]), "response_var")
})

test_that("review_conquest_overlap matches item IDs by original level when requested", {
  bundle <- build_conquest_overlap_bundle()

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      OriginalLevel = bundle$mfrmr_item_estimates$Level,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "OriginalLevel",
    conquest_item_estimate = "Est",
    item_id_source = "level",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 0)
  expect_true(all(abs(audit$item_comparison$CenteredDifference[audit$item_comparison$Status == "Compared"]) < 1e-10))
  expect_identical(as.character(audit$settings$Value[audit$settings$Setting == "item_id_source"][1]), "level")
})

test_that("review_conquest_overlap auto-detects original-level item IDs", {
  bundle <- build_conquest_overlap_bundle()

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      OriginalLevel = bundle$mfrmr_item_estimates$Level,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "OriginalLevel",
    conquest_item_estimate = "Est",
    item_id_source = "auto",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 0)
  expect_true(all(abs(audit$item_comparison$CenteredDifference[audit$item_comparison$Status == "Compared"]) < 1e-10))
  expect_identical(as.character(audit$settings$Value[audit$settings$Setting == "item_id_source"][1]), "level")
})

test_that("review_conquest_overlap auto item matching breaks ties toward response variables", {
  bundle <- build_conquest_overlap_bundle()
  item_ids <- as.character(bundle$mfrmr_item_estimates$ResponseVar)
  item_ids[(floor(length(item_ids) / 2) + 1):length(item_ids)] <-
    as.character(bundle$mfrmr_item_estimates$Level[(floor(length(item_ids) / 2) + 1):length(item_ids)])

  audit <- review_conquest_overlap(
    bundle = bundle,
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      Label = item_ids,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "Label",
    conquest_item_estimate = "Est",
    item_id_source = "auto",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_identical(as.character(audit$settings$Value[audit$settings$Setting == "item_id_source"][1]), "response_var")
  expect_true(any(audit$item_comparison$Status == "MissingInConQuest"))
  expect_true(any(audit$attention_items$Issue == "missing_conquest_item"))
})

test_that("review_conquest_overlap treats explicit item-source mismatches as review attention", {
  bundle <- build_conquest_overlap_bundle()
  audit_with_items <- function(item_ids, item_id_source) {
    review_conquest_overlap(
      bundle = bundle,
      conquest_population = data.frame(
        Term = bundle$mfrmr_population$Parameter,
        Est = bundle$mfrmr_population$Estimate,
        stringsAsFactors = FALSE
      ),
      conquest_item_estimates = data.frame(
        Label = item_ids,
        Est = bundle$mfrmr_item_estimates$Estimate,
        stringsAsFactors = FALSE
      ),
      conquest_case_eap = data.frame(
        PID = bundle$mfrmr_case_eap$Person,
        EAP = bundle$mfrmr_case_eap$Estimate,
        stringsAsFactors = FALSE
      ),
      conquest_population_term = "Term",
      conquest_population_estimate = "Est",
      conquest_item_id = "Label",
      conquest_item_estimate = "Est",
      item_id_source = item_id_source,
      conquest_case_person = "PID",
      conquest_case_estimate = "EAP"
    )
  }

  response_mismatch <- audit_with_items(bundle$mfrmr_item_estimates$Level, "response_var")
  level_mismatch <- audit_with_items(bundle$mfrmr_item_estimates$ResponseVar, "level")

  expect_s3_class(response_mismatch, "mfrm_conquest_overlap_review")
  expect_identical(
    as.character(response_mismatch$settings$Value[response_mismatch$settings$Setting == "item_id_source"][1]),
    "response_var"
  )
  expect_true(any(response_mismatch$item_comparison$Status == "MissingInConQuest"))
  expect_true(any(response_mismatch$attention_items$Issue == "missing_conquest_item"))
  expect_s3_class(level_mismatch, "mfrm_conquest_overlap_review")
  expect_identical(
    as.character(level_mismatch$settings$Value[level_mismatch$settings$Setting == "item_id_source"][1]),
    "level"
  )
  expect_true(any(level_mismatch$item_comparison$Status == "MissingInConQuest"))
  expect_true(any(level_mismatch$attention_items$Issue == "missing_conquest_item"))
})

test_that("normalize_conquest_overlap_files reads extracted csv/tsv tables", {
  bundle <- build_conquest_overlap_bundle()
  out_dir <- file.path(tempdir(), "mfrmr-conquest-normalize-files")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pop_path <- file.path(out_dir, "cq_pop.csv")
  item_path <- file.path(out_dir, "cq_item.tsv")
  case_path <- file.path(out_dir, "cq_case.txt")
  nonnumeric_case <- bundle$mfrmr_case_eap$Person[2]
  cq_case <- data.frame(
    PID = bundle$mfrmr_case_eap$Person,
    EAP = as.character(bundle$mfrmr_case_eap$Estimate),
    stringsAsFactors = FALSE
  )
  cq_case$EAP[cq_case$PID == nonnumeric_case] <- "not_a_number"

  utils::write.csv(
    data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    pop_path,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(
      Item = bundle$mfrmr_item_estimates$ResponseVar,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    item_path,
    sep = "\t",
    row.names = FALSE
  )
  utils::write.table(
    cq_case,
    case_path,
    sep = ";",
    row.names = FALSE
  )

  normalized <- normalize_conquest_overlap_files(
    population_file = pop_path,
    item_file = item_path,
    case_file = case_path,
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "Item",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )

  expect_s3_class(normalized, "mfrm_conquest_overlap_tables")
  expect_true(is.data.frame(normalized$source_files))
  expect_equal(normalized$source_files$Delimiter[[1]], ",")
  expect_equal(normalized$source_files$Delimiter[[2]], "\t")
  expect_equal(normalized$source_files$Delimiter[[3]], ";")
  expect_equal(normalized$summary$CaseNonNumeric[[1]], 1)
  expect_true(any(normalized$conquest_case_eap$EstimateNonNumeric))

  audit <- review_conquest_overlap(bundle, normalized)
  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 1)
  expect_equal(audit$overall$AttentionMissing[[1]], 0)
  expect_equal(audit$overall$AttentionDuplicate[[1]], 0)
  expect_equal(audit$overall$AttentionNonNumeric[[1]], 1)
  expect_identical(
    audit$case_comparison$Status[audit$case_comparison$Person == nonnumeric_case][1],
    "NonNumericInConQuest"
  )
})

test_that("ConQuest overlap helpers auto-resolve conservative alias columns", {
  bundle <- build_conquest_overlap_bundle()

  normalized <- normalize_conquest_overlap_tables(
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      Label = bundle$mfrmr_item_estimates$ResponseVar,
      Facility = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP_1 = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    )
  )

  expect_s3_class(normalized, "mfrm_conquest_overlap_tables")
  expect_equal(names(normalized$conquest_population)[1:2], c("Parameter", "Estimate"))
  expect_equal(names(normalized$conquest_item_estimates)[1:2], c("ItemID", "Estimate"))
  expect_equal(names(normalized$conquest_case_eap)[1:2], c("Person", "Estimate"))

  audit <- review_conquest_overlap(
    bundle,
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      Label = bundle$mfrmr_item_estimates$ResponseVar,
      Facility = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP_1 = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    )
  )

  expect_s3_class(audit, "mfrm_conquest_overlap_review")
  expect_equal(audit$overall$AttentionItems[[1]], 0)
})

test_that("build_mfrm_replay_script preserves keep_original and rating range", {
  dat <- mfrmr:::sample_mfrm_data(seed = 42) |>
    dplyr::filter(.data$Score %in% c(1, 3, 5))

  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Task", "Criterion"),
    "Score",
    method = "JML",
    maxit = 25,
    keep_original = TRUE
  ))

  replay <- build_mfrm_replay_script(fit, data_file = "analysis_data.csv")

  expect_match(replay$script, "keep_original = TRUE", fixed = TRUE)
  expect_match(replay$script, "rating_min = 1", fixed = TRUE)
  expect_match(replay$script, "rating_max = 5", fixed = TRUE)
  expect_match(replay$script, "# Model: RSM | Method: JML | InternalMethod: JMLE", fixed = TRUE)
  expect_match(replay$script, "# population_active = FALSE", fixed = TRUE)
  expect_match(replay$script, "# posterior_basis = legacy_mml", fixed = TRUE)
  expect_match(replay$script, 'method = "JML"', fixed = TRUE)
})

test_that("export_mfrm_bundle writes requested tables and html output", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-test")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expect_no_warning(
    bundle <- export_mfrm_bundle(
      fit = export_core_fixture$fit,
      diagnostics = export_core_fixture$diagnostics,
      bias_results = export_bias_fixture$bias_all,
      output_dir = out_dir,
      prefix = "bundle_test",
      include = c("core_tables", "checklist", "dashboard", "apa", "anchors", "manifest", "visual_summaries", "script", "html"),
      overwrite = TRUE
    )
  )

  expect_s3_class(bundle, "mfrm_export_bundle")
  expect_true(is.data.frame(bundle$written_files))
  expect_true(any(bundle$written_files$Component == "bundle_html"))
  expect_true(any(grepl("bundle_test_manifest_summary.csv$", bundle$written_files$Path)))
  expect_true(any(grepl("bundle_test_checklist.csv$", bundle$written_files$Path)))
  expect_true(any(grepl("bundle_test_facet_dashboard_detail.csv$", bundle$written_files$Path)))
  expect_true(any(grepl("bundle_test_replay.R$", bundle$written_files$Path)))
  expect_true(any(grepl("bundle_test_visual_warning_counts.csv$", bundle$written_files$Path)))
  expect_true(file.exists(file.path(out_dir, "bundle_test_bundle.html")))
  expect_true(file.exists(file.path(out_dir, "bundle_test_manifest.txt")))
  expect_true(file.exists(file.path(out_dir, "bundle_test_replay.R")))
  expect_true(file.exists(file.path(out_dir, "bundle_test_visual_warning_map.txt")))
})

test_that("export_mfrm_bundle writes optional prediction artifacts", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-predictions")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  bundle <- export_mfrm_bundle(
    fit = prediction_bundle_fixture$fit,
    diagnostics = prediction_bundle_fixture$diagnostics,
    population_prediction = prediction_bundle_fixture$population_prediction,
    unit_prediction = prediction_bundle_fixture$unit_prediction,
    plausible_values = prediction_bundle_fixture$plausible_values,
    output_dir = out_dir,
    prefix = "bundle_pred_test",
    include = c("manifest", "predictions", "html"),
    overwrite = TRUE
  )

  expect_s3_class(bundle, "mfrm_export_bundle")
  expect_true(any(bundle$written_files$Component == "population_prediction_forecast"))
  expect_true(any(bundle$written_files$Component == "unit_prediction_estimates"))
  expect_true(any(bundle$written_files$Component == "plausible_values"))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_population_prediction_forecast.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_unit_prediction_estimates.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_plausible_values.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_population_prediction_ademp.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_population_prediction_sim_spec_settings.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_population_prediction_sim_spec_thresholds.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_unit_prediction_input.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_plausible_value_input.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_pred_test_bundle.html")))
  population_settings <- utils::read.csv(
    file.path(out_dir, "bundle_pred_test_population_prediction_settings.csv"),
    stringsAsFactors = FALSE
  )
  population_ademp <- utils::read.csv(
    file.path(out_dir, "bundle_pred_test_population_prediction_ademp.csv"),
    stringsAsFactors = FALSE
  )
  expect_identical(
    population_settings$Value[population_settings$Setting == "planning_schema"][1],
    "omitted_from_export_flattening; see compact population_prediction_sim_spec files"
  )
  expect_identical(
    population_ademp$Value[population_ademp$Key == "data_generating_mechanism.planning_schema"][1],
    "omitted_from_export_flattening; see compact population_prediction_sim_spec files"
  )

  html_lines <- readLines(file.path(out_dir, "bundle_pred_test_bundle.html"), warn = FALSE)
  html_text <- paste(html_lines, collapse = "\n")
  expect_match(html_text, "<h2>population_prediction_forecast</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>unit_prediction_estimates</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>plausible_value_summary</h2>", fixed = TRUE)

  unit_settings <- utils::read.csv(
    file.path(out_dir, "bundle_pred_test_unit_prediction_settings.csv"),
    stringsAsFactors = FALSE
  )
  expect_true(any(unit_settings$Setting == "source_columns.person"))
  expect_false(any(grepl("<list", unit_settings$Value, fixed = TRUE)))
})

test_that("export_mfrm_bundle writes latent-regression scoring provenance artifacts", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-latent-predictions")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  bundle <- export_mfrm_bundle(
    fit = latent_prediction_omit_bundle_fixture$fit,
    unit_prediction = latent_prediction_omit_bundle_fixture$unit_prediction,
    plausible_values = latent_prediction_omit_bundle_fixture$plausible_values,
    output_dir = out_dir,
    prefix = "bundle_latent_pred_test",
    include = c("manifest", "predictions", "script", "html"),
    overwrite = TRUE
  )

  expect_s3_class(bundle, "mfrm_export_bundle")
  expect_true(any(bundle$written_files$Component == "manifest_settings"))
  expect_true(any(bundle$written_files$Component == "replay_fit_person_data"))
  expect_true(any(bundle$written_files$Component == "unit_prediction_population_review"))
  expect_true(any(bundle$written_files$Component == "unit_prediction_person_data"))
  expect_true(any(bundle$written_files$Component == "plausible_value_population_review"))
  expect_true(any(bundle$written_files$Component == "plausible_value_person_data"))
  expect_true(file.exists(file.path(out_dir, "bundle_latent_pred_test_replay_fit_person_data.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_latent_pred_test_manifest_settings.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_latent_pred_test_replay.R")))
  expect_true(file.exists(file.path(out_dir, "bundle_latent_pred_test_bundle.html")))
  expect_true(file.exists(file.path(out_dir, "bundle_latent_pred_test_unit_prediction_population_review.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_latent_pred_test_unit_prediction_person_data.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_latent_pred_test_plausible_value_population_review.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_latent_pred_test_plausible_value_person_data.csv")))

  fit_person_data <- utils::read.csv(
    file.path(out_dir, "bundle_latent_pred_test_replay_fit_person_data.csv"),
    stringsAsFactors = FALSE
  )
  replay_lines <- readLines(
    file.path(out_dir, "bundle_latent_pred_test_replay.R"),
    warn = FALSE
  )
  replay_text <- paste(replay_lines, collapse = "\n")
  html_lines <- readLines(
    file.path(out_dir, "bundle_latent_pred_test_bundle.html"),
    warn = FALSE
  )
  html_text <- paste(html_lines, collapse = "\n")
  pop_audit <- utils::read.csv(
    file.path(out_dir, "bundle_latent_pred_test_unit_prediction_population_review.csv"),
    stringsAsFactors = FALSE
  )
  person_data <- utils::read.csv(
    file.path(out_dir, "bundle_latent_pred_test_unit_prediction_person_data.csv"),
    stringsAsFactors = FALSE
  )
  manifest_settings <- utils::read.csv(
    file.path(out_dir, "bundle_latent_pred_test_manifest_settings.csv"),
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(fit_person_data), 60L)
  expect_true(all(c("Person", "X") %in% names(fit_person_data)))
  expect_equal(sum(is.na(fit_person_data$X)), 1L)
  expect_identical(pop_audit$Policy[[1]], "omit")
  expect_equal(pop_audit$OmittedPersons[[1]], 1)
  expect_equal(nrow(person_data), 2L)
  expect_true(all(c("Person", "X") %in% names(person_data)))
  expect_true(any(is.na(person_data$X)))
  expect_true(nrow(fit_person_data) > nrow(person_data))
  expect_match(replay_text, "replay_script_args <- commandArgs", fixed = TRUE)
  expect_match(replay_text, "sys.frames()[[1]]$ofile", fixed = TRUE)
  expect_match(replay_text, "bundle_latent_pred_test_replay_fit_person_data.csv", fixed = TRUE)
  expect_match(replay_text, "fit_person_data <- utils::read.csv", fixed = TRUE)
  expect_match(html_text, "<h2>replay_artifacts</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>replay_script</h2>", fixed = TRUE)
  expect_match(html_text, "bundle_latent_pred_test_replay_fit_person_data.csv", fixed = TRUE)
  expect_match(html_text, "fit_person_data &lt;- utils::read.csv", fixed = TRUE)
  expect_false(grepl("fit_person_data &lt;- structure\\(", html_text))
  expect_false(grepl("<h2>replay_fit_person_data</h2>", html_text, fixed = TRUE))
  expect_identical(
    as.character(manifest_settings$Value[manifest_settings$Setting == "unit_prediction_population_policy"][1]),
    "omit"
  )
  expect_identical(
    as.character(manifest_settings$Value[manifest_settings$Setting == "plausible_value_population_policy"][1]),
    "omit"
  )
})

test_that("export_mfrm_bundle writes default summary-table bundles for manuscript handoff", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-summary-tables")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  bundle <- export_mfrm_bundle(
    fit = export_core_fixture$fit,
    diagnostics = export_core_fixture$diagnostics,
    bias_results = export_bias_fixture$bias_all,
    output_dir = out_dir,
    prefix = "bundle_summary_test",
    include = c("summary_tables", "html"),
    overwrite = TRUE
  )
  bundle_summary <- summary(bundle)

  expect_s3_class(bundle, "mfrm_export_bundle")
  expect_s3_class(bundle_summary, "summary.mfrm_bundle")
  expect_identical(bundle_summary$preview_name, "written_files")
  expect_equal(bundle_summary$summary$HtmlWritten[[1]], 1)
  expect_true(is.data.frame(bundle_summary$format_summary))
  expect_true(is.data.frame(bundle_summary$artifact_catalog))
  expect_true(is.data.frame(bundle_summary$reporting_map))
  expect_true(any(bundle_summary$artifact_catalog$ArtifactGroup == "summary_surface"))
  expect_true(any(bundle_summary$artifact_catalog$ArtifactGroup == "html_review"))
  bundle_plot <- plot(bundle, type = "formats", draw = FALSE)
  bundle_artifact_plot <- plot(bundle, type = "artifact_groups", draw = FALSE)
  expect_s3_class(bundle_plot, "mfrm_plot_data")
  expect_identical(bundle_plot$name, "export_bundle")
  expect_identical(bundle_plot$data$plot, "formats")
  expect_s3_class(bundle_artifact_plot, "mfrm_plot_data")
  expect_identical(bundle_artifact_plot$data$plot, "artifact_groups")
  expect_true(any(bundle$written_files$Component == "summary_fit_table_index"))
  expect_true(any(bundle$written_files$Component == "summary_fit_table_catalog"))
  expect_true(any(bundle$written_files$Component == "summary_fit_reporting_map"))
  expect_true(any(bundle$written_files$Component == "summary_fit_facet_overview"))
  expect_true(any(bundle$written_files$Component == "summary_diagnostics_flags"))
  expect_true(any(bundle$written_files$Component == "summary_checklist_action_items"))
  expect_true(any(bundle$written_files$Component == "summary_apa_components"))
  expect_true(file.exists(file.path(out_dir, "bundle_summary_test_summary_fit_table_index.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_summary_test_summary_fit_table_catalog.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_summary_test_summary_fit_reporting_map.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_summary_test_summary_fit_facet_overview.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_summary_test_summary_diagnostics_flags.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_summary_test_summary_checklist_action_items.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_summary_test_summary_apa_components.csv")))
  expect_true(file.exists(file.path(out_dir, "bundle_summary_test_bundle.html")))

  html_lines <- readLines(file.path(out_dir, "bundle_summary_test_bundle.html"), warn = FALSE)
  html_text <- paste(html_lines, collapse = "\n")
  expect_match(html_text, "<h2>summary_fit_table_index</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>summary_fit_table_catalog</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>summary_fit_reporting_map</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>summary_checklist_action_items</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>summary_apa_components</h2>", fixed = TRUE)
})

test_that("export_summary_appendix writes appendix-ready summary artifacts without requiring fit export inputs", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-export")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  chk <- reporting_checklist(export_core_fixture$fit, diagnostics = export_core_fixture$diagnostics)
  appendix <- export_summary_appendix(
    list(
      fit = summary(export_core_fixture$fit),
      diagnostics = export_core_fixture$diagnostics,
      checklist = chk
    ),
    output_dir = out_dir,
    prefix = "appendix_test",
    include_html = TRUE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  appendix_summary <- summary(appendix)
  expect_s3_class(appendix_summary, "summary.mfrm_bundle")
  expect_identical(appendix_summary$preview_name, "written_files")
  expect_true(is.data.frame(appendix_summary$format_summary))
  expect_true(is.data.frame(appendix_summary$artifact_catalog))
  expect_true(is.data.frame(appendix_summary$selection_summary))
  expect_true(is.data.frame(appendix_summary$selection_table_summary))
  expect_true(is.data.frame(appendix_summary$selection_handoff_table_summary))
  expect_true(is.data.frame(appendix_summary$selection_handoff_preset_summary))
  expect_true(is.data.frame(appendix_summary$selection_handoff_summary))
  expect_true(is.data.frame(appendix_summary$selection_handoff_role_summary))
  expect_true(is.data.frame(appendix_summary$selection_handoff_role_section_summary))
  expect_true(is.data.frame(appendix_summary$selection_role_summary))
  expect_true(is.data.frame(appendix_summary$selection_section_summary))
  expect_true(is.data.frame(appendix_summary$selection_catalog))
  expect_true(is.data.frame(appendix_summary$reporting_map))
  expect_true(all(c("TablesAvailable", "SelectionFraction", "PlotReadyFraction", "NumericFraction") %in% names(appendix_summary$selection_summary)))
  expect_true(all(c("Preset", "AppendixSection", "Role", "Bundle", "Table", "Rows", "NumericColumns", "PlotReady", "ExportReady", "ApaTableReady") %in%
                    names(appendix_summary$selection_handoff_table_summary)))
  expect_true(all(c("PlotReadyFraction", "NumericFraction") %in% names(appendix_summary$selection_handoff_summary)))
  expect_true(all(c("PlotReadyFraction", "NumericFraction") %in% names(appendix_summary$selection_handoff_role_summary)))
  expect_true(all(c("PlotReadyFraction", "NumericFraction") %in% names(appendix_summary$selection_handoff_role_section_summary)))
  expect_true(all(c("PlotReadyFraction", "NumericFraction") %in% names(appendix_summary$selection_role_summary)))
  expect_true(all(c("PlotReadyFraction", "NumericFraction") %in% names(appendix_summary$selection_section_summary)))
  expect_true(any(appendix_summary$artifact_catalog$ArtifactGroup == "summary_surface"))
  expect_true(any(appendix_summary$artifact_catalog$ArtifactGroup == "html_review"))
  appendix_plot <- plot(appendix, type = "formats", draw = FALSE)
  expect_s3_class(appendix_plot, "mfrm_plot_data")
  expect_identical(appendix_plot$name, "export_bundle")
  expect_identical(appendix_plot$data$plot, "formats")
  appendix_selection_plot <- plot(appendix, type = "selection_bundles", draw = FALSE)
  expect_s3_class(appendix_selection_plot, "mfrm_plot_data")
  expect_identical(appendix_selection_plot$data$plot, "selection_bundles")
  appendix_table_plot <- plot(appendix, type = "selection_tables", draw = FALSE)
  expect_s3_class(appendix_table_plot, "mfrm_plot_data")
  expect_identical(appendix_table_plot$data$plot, "selection_tables")
  appendix_handoff_preset_plot <- plot(appendix, type = "selection_handoff_presets", draw = FALSE)
  expect_s3_class(appendix_handoff_preset_plot, "mfrm_plot_data")
  expect_identical(appendix_handoff_preset_plot$data$plot, "selection_handoff_presets")
  appendix_handoff_plot <- plot(appendix, type = "selection_handoff", draw = FALSE)
  expect_s3_class(appendix_handoff_plot, "mfrm_plot_data")
  expect_identical(appendix_handoff_plot$data$plot, "selection_handoff")
  appendix_handoff_fraction_plot <- plot(appendix, type = "selection_handoff", selection_value = "fraction", draw = FALSE)
  expect_s3_class(appendix_handoff_fraction_plot, "mfrm_plot_data")
  expect_identical(appendix_handoff_fraction_plot$data$plot, "selection_handoff")
  expect_identical(appendix_handoff_fraction_plot$data$selection_value, "fraction")
  appendix_handoff_bundle_plot <- plot(appendix, type = "selection_handoff_bundles", draw = FALSE)
  expect_s3_class(appendix_handoff_bundle_plot, "mfrm_plot_data")
  expect_identical(appendix_handoff_bundle_plot$data$plot, "selection_handoff_bundles")
  appendix_handoff_role_plot <- plot(appendix, type = "selection_handoff_roles", draw = FALSE)
  expect_s3_class(appendix_handoff_role_plot, "mfrm_plot_data")
  expect_identical(appendix_handoff_role_plot$data$plot, "selection_handoff_roles")
  appendix_handoff_role_section_plot <- plot(appendix, type = "selection_handoff_role_sections", draw = FALSE)
  expect_s3_class(appendix_handoff_role_section_plot, "mfrm_plot_data")
  expect_identical(appendix_handoff_role_section_plot$data$plot, "selection_handoff_role_sections")
  appendix_section_plot <- plot(appendix, type = "selection_sections", draw = FALSE)
  expect_s3_class(appendix_section_plot, "mfrm_plot_data")
  expect_identical(appendix_section_plot$data$plot, "selection_sections")
  expect_error(
    plot(appendix, type = "selection_tables", selection_value = "fraction", draw = FALSE),
    "not available for `type = \"selection_tables\"`",
    fixed = TRUE
  )
  expect_identical(
    unique(appendix$written_files$Component[appendix$written_files$Component == "summary_fit_reporting_map"]),
    "summary_fit_reporting_map"
  )
  expect_equal(
    sum(appendix$written_files$Component == "summary_fit_reporting_map"),
    1L
  )
  expect_true(any(appendix$written_files$Component == "appendix_selection_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_table_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_handoff_table_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_handoff_preset_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_handoff_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_handoff_role_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_handoff_role_section_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_role_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_section_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_catalog"))
  expect_true(any(appendix$written_files$Component == "summary_fit_table_catalog"))
  expect_true(any(appendix$written_files$Component == "summary_fit_reporting_map"))
  expect_true(any(appendix$written_files$Component == "appendix_html"))
  expect_true(file.exists(file.path(out_dir, "appendix_test_summary_fit_table_catalog.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_summary_fit_reporting_map.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix_selection_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix_selection_table_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix_selection_handoff_table_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix_selection_handoff_preset_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix_selection_handoff_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix_selection_handoff_role_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix_selection_handoff_role_section_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix_selection_role_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix_selection_section_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix_selection_catalog.csv")))
  expect_true(file.exists(file.path(out_dir, "appendix_test_appendix.html")))
  expect_match(as.character(appendix_summary$overview$Class[1]), "mfrm_summary_appendix_export", fixed = TRUE)

  html_lines <- readLines(file.path(out_dir, "appendix_test_appendix.html"), warn = FALSE)
  html_text <- paste(html_lines, collapse = "\n")
  expect_match(html_text, "<h2>summary_fit_table_catalog</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>summary_fit_reporting_map</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>summary_checklist_action_items</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>appendix_selection_summary</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>appendix_selection_table_summary</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>appendix_selection_handoff_table_summary</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>appendix_selection_handoff_preset_summary</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>appendix_selection_handoff_summary</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>appendix_selection_handoff_bundle_summary</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>appendix_selection_handoff_role_summary</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>appendix_selection_handoff_role_section_summary</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>appendix_selection_role_summary</h2>", fixed = TRUE)
  expect_match(html_text, "<h2>appendix_selection_section_summary</h2>", fixed = TRUE)
})

test_that("export_summary_appendix preset trims bridge-only and preview-only tables", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-recommended")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  chk <- reporting_checklist(export_core_fixture$fit, diagnostics = export_core_fixture$diagnostics)
  apa <- build_apa_outputs(export_core_fixture$fit, diagnostics = export_core_fixture$diagnostics)

  appendix <- export_summary_appendix(
    list(
      fit = summary(export_core_fixture$fit),
      checklist = chk,
      apa = apa
    ),
    output_dir = out_dir,
    prefix = "appendix_recommended",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(all(appendix$selection_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_table_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_handoff_preset_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_handoff_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_handoff_role_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_handoff_role_section_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_role_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_section_summary$Preset == "recommended"))
  expect_true(any(appendix$written_files$Component == "summary_fit_reporting_map"))
  expect_false(any(appendix$written_files$Component == "summary_checklist_action_items"))
  expect_false(any(appendix$written_files$Component == "summary_apa_preview"))
  expect_true(any(appendix$written_files$Component == "summary_fit_overview"))
  expect_true(any(appendix$written_files$Component == "summary_checklist_section_summary"))
  expect_true(any(appendix$written_files$Component == "summary_apa_components"))
  expect_false(any(
    appendix$selection_catalog$Bundle == "fit" &
      appendix$selection_catalog$Table == "reporting_map" &
      appendix$selection_catalog$Selected %in% TRUE
  ))
  expect_true(any(appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(appendix$selection_catalog$Selected %in% FALSE))
  expect_true(all(appendix$selection_catalog$Preset == "recommended"))
  expect_true(any(appendix$selection_section_summary$AppendixSection %in% c("methods", "results", "diagnostics", "reporting")))
})

test_that("export_summary_appendix supports weighting-review inputs", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-weighting")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  appendix <- export_summary_appendix(
    force_export_fixture(weighting_appendix_fixture),
    output_dir = out_dir,
    prefix = "appendix_weighting",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(all(appendix$selection_summary$Preset == "recommended"))
  expect_true(any(appendix$selection_catalog$Table == "overview" & appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(appendix$selection_catalog$Table == "top_reweighted_levels" & appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(grepl("weighting", appendix$written_files$Component, fixed = TRUE)))
  expect_true(any(grepl("top_reweighted_levels", appendix$written_files$Component, fixed = TRUE)))
})

test_that("export_summary_appendix supports misfit-casebook review inputs", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-misfit-casebook")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  appendix <- export_summary_appendix(
    force_export_fixture(misfit_appendix_fixture),
    output_dir = out_dir,
    prefix = "appendix_casebook",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(all(appendix$selection_summary$Preset == "recommended"))
  expect_true(any(appendix$selection_catalog$Table == "overview" & appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(appendix$selection_catalog$Table == "top_cases" & appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(appendix$selection_catalog$Table == "case_rollup" & appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(appendix$selection_catalog$Table == "group_view_index" & appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(grepl("casebook", appendix$written_files$Component, fixed = TRUE)))
  expect_true(any(grepl("top_cases", appendix$written_files$Component, fixed = TRUE)))
})

test_that("export_summary_appendix supports recovery simulation and assessment inputs", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-recovery")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  recovery_inputs <- force_export_fixture(recovery_appendix_fixture)

  direct_appendix <- export_summary_appendix(
    recovery_inputs$recovery,
    output_dir = out_dir,
    prefix = "appendix_recovery_direct",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )
  expect_s3_class(direct_appendix, "mfrm_summary_appendix_export")
  expect_true(any(direct_appendix$written_files$Component ==
                    "summary_mfrm_recovery_simulation_recovery_summary"))
  expect_true(any(direct_appendix$selection_catalog$Role == "recovery_performance" &
                    direct_appendix$selection_catalog$Selected %in% TRUE))

  appendix <- export_summary_appendix(
    list(
      recovery = recovery_inputs$recovery,
      assessment = recovery_inputs$assessment
    ),
    output_dir = out_dir,
    prefix = "appendix_recovery",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(all(appendix$selection_summary$Preset == "recommended"))
  expect_true(any(appendix$written_files$Component == "summary_recovery_recovery_summary"))
  expect_true(any(appendix$written_files$Component == "summary_recovery_rep_overview"))
  expect_true(any(appendix$written_files$Component == "summary_recovery_ademp"))
  expect_true(any(appendix$written_files$Component == "summary_assessment_checklist"))
  expect_true(any(appendix$written_files$Component == "summary_assessment_metric_review"))
  expect_true(any(appendix$written_files$Component == "summary_assessment_thresholds"))
  expect_true(any(appendix$selection_catalog$Role == "recovery_design_basis" &
                    appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(appendix$selection_catalog$Role == "recovery_assessment_checklist" &
                    appendix$selection_catalog$Selected %in% TRUE))
  expect_true(file.exists(file.path(
    out_dir,
    "appendix_recovery_summary_assessment_metric_review.csv"
  )))
})

test_that("export_summary_appendix supports recovery-validation summaries", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-recovery-validation")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  validation_summary <- export_recovery_validation_summary_fixture()

  appendix <- export_summary_appendix(
    list(validation = validation_summary),
    output_dir = out_dir,
    prefix = "appendix_recovery_validation",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(any(appendix$written_files$Component ==
                    "summary_validation_topline_release_decision"))
  expect_true(any(appendix$written_files$Component ==
                    "summary_validation_condition_reporting_notes"))
  expect_true(any(appendix$written_files$Component ==
                    "summary_validation_diagnostic_reporting_notes"))
  expect_true(any(appendix$selection_catalog$Role ==
                    "recovery_validation_diagnostic_reporting_notes" &
                    appendix$selection_catalog$Selected %in% TRUE))
  expect_true(file.exists(file.path(
    out_dir,
    "appendix_recovery_validation_summary_validation_diagnostic_reporting_notes.csv"
  )))
})

test_that("export_summary_appendix supports diagnostic-screening summaries", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-diagnostic-screening")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  diag_summary <- diagnostic_screening_appendix_summary_fixture()

  appendix <- export_summary_appendix(
    list(diag = diag_summary),
    output_dir = out_dir,
    prefix = "appendix_diagnostic",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(any(appendix$written_files$Component ==
                    "summary_diag_scenario_summary"))
  expect_true(any(appendix$written_files$Component ==
                    "summary_diag_reading_order"))
  expect_true(any(appendix$written_files$Component ==
                    "summary_diag_next_actions"))
  expect_true(any(appendix$written_files$Component ==
                    "summary_diag_reporting_notes"))
  expect_true(any(appendix$written_files$Component ==
                    "summary_diag_figure_recipes"))
  expect_true(any(appendix$selection_catalog$Role ==
                    "diagnostic_screening_next_actions" &
                    appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(appendix$written_files$Component ==
                    "summary_diag_plot_overview_rate"))
  expect_true(any(appendix$selection_catalog$Role ==
                    "diagnostic_screening_reporting_notes" &
                    appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(appendix$selection_catalog$Role ==
                    "diagnostic_screening_figure_recipes" &
                    appendix$selection_catalog$Selected %in% TRUE))
  expect_true(any(appendix$selection_catalog$Role ==
                    "diagnostic_screening_plot_data" &
                    appendix$selection_catalog$Selected %in% TRUE))
  expect_true(file.exists(file.path(
    out_dir,
    "appendix_diagnostic_summary_diag_plot_overview_rate.csv"
  )))
})

test_that("export_summary_appendix supports person-fit summary inputs", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-person-fit")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  person_fit_summary <- list(
    overview = data.frame(
      Persons = 2L,
      ReportableRows = 1L,
      FlaggedRows = 1L,
      stringsAsFactors = FALSE
    ),
    status_summary = data.frame(
      Status = "review",
      Rows = 1L,
      stringsAsFactors = FALSE
    ),
    report_index_summary = data.frame(
      ReportIndex = "lz_star",
      Rows = 1L,
      stringsAsFactors = FALSE
    ),
    lz_star_status_summary = data.frame(
      Status = "available",
      Rows = 1L,
      stringsAsFactors = FALSE
    ),
    top_review = data.frame(
      Person = "P1",
      ReportIndex = "lz_star",
      ReportValue = -2.5,
      Status = "review",
      stringsAsFactors = FALSE
    ),
    caveats = data.frame(),
    thresholds = data.frame(
      Threshold = "z",
      Value = 2,
      stringsAsFactors = FALSE
    ),
    reporting_map = data.frame(
      Output = "top_review",
      Use = "response-level follow-up",
      stringsAsFactors = FALSE
    ),
    notes = "Person-fit rows are diagnostic follow-up, not automatic exclusions."
  )
  class(person_fit_summary) <- "summary.mfrm_person_fit_indices"

  appendix <- export_summary_appendix(
    person_fit_summary,
    output_dir = out_dir,
    prefix = "appendix_person_fit",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(any(appendix$written_files$Component ==
                    "summary_summary_mfrm_person_fit_indices_status_summary"))
  expect_true(any(appendix$selection_catalog$Role == "review_status" &
                    appendix$selection_catalog$Selected %in% TRUE))
  expect_true(file.exists(file.path(
    out_dir,
    "appendix_person_fit_summary_summary_mfrm_person_fit_indices_status_summary.csv"
  )))
})

test_that("export_summary_appendix supports precision-review summaries", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-precision-review")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  precision <- precision_review_report(
    export_core_fixture$fit,
    diagnostics = export_core_fixture$diagnostics
  )

  appendix <- export_summary_appendix(
    list(precision = summary(precision)),
    output_dir = out_dir,
    prefix = "appendix_precision",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(any(appendix$written_files$Component ==
                    "summary_precision_fit_separation_basis"))
  expect_true(any(appendix$selection_catalog$Role == "precision_review" &
                    appendix$selection_catalog$Selected %in% TRUE))
  expect_true(file.exists(file.path(
    out_dir,
    "appendix_precision_summary_precision_fit_separation_basis.csv"
  )))
})

test_that("export_summary_appendix supports fit-measure and FACETS fit-review summaries", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-fit-review")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fm <- fit_measures_table(
    export_core_fixture$fit,
    diagnostics = export_core_fixture$diagnostics,
    fit_df_method = "both",
    top_n = 5
  )
  facets_review <- facets_fit_review(
    export_core_fixture$fit,
    diagnostics = export_core_fixture$diagnostics
  )

  appendix <- export_summary_appendix(
    list(
      fit_measures = summary(fm),
      facets_fit = facets_review
    ),
    output_dir = out_dir,
    prefix = "appendix_fit_review",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(any(appendix$written_files$Component ==
                    "summary_fit_measures_df_sensitivity"))
  expect_true(any(appendix$written_files$Component ==
                    "summary_facets_fit_df_sensitivity"))
  expect_true(any(appendix$selection_catalog$Role == "precision_review" &
                    appendix$selection_catalog$Selected %in% TRUE))
  expect_true(file.exists(file.path(
    out_dir,
    "appendix_fit_review_summary_fit_measures_df_sensitivity.csv"
  )))
  expect_true(file.exists(file.path(
    out_dir,
    "appendix_fit_review_summary_facets_fit_df_sensitivity.csv"
  )))
})

test_that("export_mfrm_bundle supports recovery-validation summary tables", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-recovery-validation")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  validation_summary <- export_recovery_validation_summary_fixture()
  bundle <- export_mfrm_bundle(
    fit = export_core_fixture$fit,
    diagnostics = export_core_fixture$diagnostics,
    summary_tables = list(validation = validation_summary),
    output_dir = out_dir,
    prefix = "bundle_recovery_validation",
    include = "summary_tables",
    overwrite = TRUE
  )

  expect_s3_class(bundle, "mfrm_export_bundle")
  expect_true(any(bundle$written_files$Component ==
                    "summary_validation_topline_release_decision"))
  expect_true(any(bundle$written_files$Component ==
                    "summary_validation_condition_reporting_notes"))
  expect_true(any(bundle$written_files$Component ==
                    "summary_validation_diagnostic_reporting_notes"))
  expect_true(file.exists(file.path(
    out_dir,
    "bundle_recovery_validation_summary_validation_diagnostic_reporting_notes.csv"
  )))
})

test_that("export_summary_appendix supports section-aware appendix presets", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-methods")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  appendix <- export_summary_appendix(
    list(
      fit = summary(export_core_fixture$fit),
      diagnostics = export_core_fixture$diagnostics,
      checklist = reporting_checklist(export_core_fixture$fit, diagnostics = export_core_fixture$diagnostics)
    ),
    output_dir = out_dir,
    prefix = "appendix_methods",
    preset = "methods",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(all(appendix$selection_summary$Preset == "methods"))
  expect_true(all(appendix$selection_table_summary$Preset == "methods"))
  expect_true(all(appendix$selection_handoff_preset_summary$Preset == "methods"))
  expect_true(all(appendix$selection_handoff_summary$Preset == "methods"))
  expect_true(all(appendix$selection_handoff_bundle_summary$Preset == "methods"))
  expect_true(all(appendix$selection_handoff_role_summary$Preset == "methods"))
  expect_true(all(appendix$selection_handoff_role_section_summary$Preset == "methods"))
  expect_true(all(appendix$selection_role_summary$Preset == "methods"))
  expect_true(all(appendix$selection_section_summary$Preset == "methods"))
  expect_true(any(appendix$selection_catalog$Selected %in% TRUE))
  expect_true(all(
    appendix$selection_catalog$AppendixSection[appendix$selection_catalog$Selected %in% TRUE] %in% "methods"
  ))
  expect_true(all(appendix$selection_section_summary$AppendixSection %in% "methods"))
  expect_true(any(appendix$written_files$Component == "summary_fit_overview"))
  expect_true(any(appendix$written_files$Component == "summary_diagnostics_overview"))
})

test_that("export_summary_appendix supports future arbitrary-facet active-branch inputs", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-future-branch")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  active <- spec$planning_schema$future_branch_active_branch
  appendix <- export_summary_appendix(
    active,
    output_dir = out_dir,
    prefix = "appendix_future_branch",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  appendix_summary <- summary(appendix)
  expect_true(is.data.frame(appendix_summary$selection_handoff_preset_summary))
  expect_true(is.data.frame(appendix_summary$selection_handoff_bundle_summary))
  expect_true(is.data.frame(appendix_summary$selection_handoff_role_summary))
  expect_true(is.data.frame(appendix_summary$selection_handoff_role_section_summary))
  expect_true(any(appendix$written_files$Component == "summary_mfrm_future_branch_active_branch_future_branch_overview"))
  expect_true(any(appendix$written_files$Component == "summary_mfrm_future_branch_active_branch_future_branch_recommendation"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_handoff_preset_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_handoff_bundle_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_handoff_role_summary"))
  expect_true(any(appendix$written_files$Component == "appendix_selection_handoff_role_section_summary"))
  expect_true(file.exists(file.path(
    out_dir,
    "appendix_future_branch_summary_mfrm_future_branch_active_branch_future_branch_recommendation.csv"
  )))
})

test_that("export_summary_appendix applies appendix presets to future arbitrary-facet active-branch inputs", {
  out_dir <- file.path(tempdir(), "mfrmr-summary-appendix-future-branch-recommended")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  active <- spec$planning_schema$future_branch_active_branch
  appendix <- export_summary_appendix(
    active,
    output_dir = out_dir,
    prefix = "appendix_future_branch_recommended",
    preset = "recommended",
    include_html = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(appendix, "mfrm_summary_appendix_export")
  expect_true(all(appendix$selection_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_table_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_handoff_preset_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_handoff_bundle_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_handoff_role_summary$Preset == "recommended"))
  expect_true(all(appendix$selection_handoff_role_section_summary$Preset == "recommended"))
  expect_true(any(appendix$written_files$Component == "summary_mfrm_future_branch_active_branch_future_branch_overview"))
  expect_true(any(appendix$written_files$Component == "summary_mfrm_future_branch_active_branch_future_branch_profile"))
  expect_true(any(appendix$written_files$Component == "summary_mfrm_future_branch_active_branch_future_branch_readiness"))
  expect_true(any(appendix$written_files$Component == "summary_mfrm_future_branch_active_branch_future_branch_recommendation"))
  expect_false(any(appendix$written_files$Component == "summary_mfrm_future_branch_active_branch_future_branch_load_balance"))
  expect_false(any(appendix$written_files$Component == "summary_mfrm_future_branch_active_branch_future_branch_coverage"))
  expect_false(any(appendix$written_files$Component == "summary_mfrm_future_branch_active_branch_future_branch_guardrails"))
  expect_true(all(
    appendix$selection_catalog$AppendixSection[appendix$selection_catalog$Selected %in% TRUE] %in%
      c("methods", "diagnostics")
  ))
})

test_that("export_mfrm_bundle requires explicit prediction objects for prediction export", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-predictions-missing")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expect_error(
    export_mfrm_bundle(
      fit = export_core_fixture$fit,
      diagnostics = export_core_fixture$diagnostics,
      output_dir = out_dir,
      prefix = "bundle_pred_missing",
      include = c("predictions"),
      overwrite = TRUE
    ),
    "`include = 'predictions'` requires at least one of `population_prediction`, `unit_prediction`, or `plausible_values`.",
    fixed = TRUE
  )
})

test_that("export_mfrm_bundle rejects malformed bias_results inputs early", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-bad-bias")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expect_error(
    export_mfrm_bundle(
      fit = export_core_fixture$fit,
      diagnostics = export_core_fixture$diagnostics,
      bias_results = list(bad = data.frame(x = 1)),
      output_dir = out_dir,
      prefix = "bundle_bad_bias",
      include = c("manifest"),
      overwrite = TRUE
    ),
    "`bias_results` in export helpers must be NULL, output from estimate_bias\\(\\), an `mfrm_bias_collection`, or a list of `mfrm_bias` objects."
  )
})

test_that("export_mfrm_bundle does not change the caller working directory", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-wd")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  original_wd <- getwd()
  bundle <- export_mfrm_bundle(
    fit = export_core_fixture$fit,
    diagnostics = export_core_fixture$diagnostics,
    output_dir = out_dir,
    prefix = "bundle_zip_test",
    include = c("manifest"),
    zip_bundle = TRUE,
    overwrite = TRUE
  )

  expect_identical(getwd(), original_wd)
  expect_s3_class(bundle, "mfrm_export_bundle")
  expect_true(any(bundle$written_files$Component == "bundle_zip"))
  expect_true(file.exists(file.path(out_dir, "bundle_zip_test_bundle.zip")))
})

test_that("export_mfrm_bundle respects overwrite for zip bundles", {
  out_dir <- file.path(tempdir(), "mfrmr-export-bundle-zip-overwrite")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expect_no_warning(
    export_mfrm_bundle(
      fit = export_core_fixture$fit,
      diagnostics = export_core_fixture$diagnostics,
      output_dir = out_dir,
      prefix = "bundle_zip_overwrite",
      include = c("manifest"),
      zip_bundle = TRUE,
      overwrite = TRUE
    )
  )

  expect_error(
    export_mfrm_bundle(
      fit = export_core_fixture$fit,
      diagnostics = export_core_fixture$diagnostics,
      output_dir = out_dir,
      prefix = "bundle_zip_overwrite",
      include = c("manifest"),
      zip_bundle = TRUE,
      overwrite = FALSE
    ),
    "File already exists:",
    fixed = TRUE
  )
})
