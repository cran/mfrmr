summary_table_bundle_workflow_fixture <- local({
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:14]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  fit <- suppressWarnings(fit_mfrm(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  run <- suppressWarnings(run_mfrm_facets(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 20
  ))
  bias <- suppressWarnings(estimate_bias(
    fit,
    diagnostics = diag,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))
  audit <- audit_mfrm_anchors(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score"
  )
  toy_wave_a <- toy[toy$Person %in% keep_people[1:7], , drop = FALSE]
  toy_wave_b <- toy[toy$Person %in% keep_people[8:14], , drop = FALSE]
  fit_wave_a <- suppressWarnings(fit_mfrm(
    toy_wave_a,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 20
  ))
  fit_wave_b <- suppressWarnings(fit_mfrm(
    toy_wave_b,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 20
  ))
  drift <- suppressWarnings(detect_anchor_drift(list(W1 = fit_wave_a, W2 = fit_wave_b)))
  chain <- suppressWarnings(build_equating_chain(list(W1 = fit_wave_a, W2 = fit_wave_b)))
  linking_review <- build_linking_review(anchor_audit = audit, drift = drift, chain = chain)

  list(
    run = run,
    bias = bias,
    audit = audit,
    linking_review = linking_review
  )
})

summary_table_bundle_prediction_fixture <- local({
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:18]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    quad_points = 5,
    maxit = 15
  ))

  new_units <- data.frame(
    Person = c("NEW01", "NEW01"),
    Rater = unique(toy$Rater)[1],
    Criterion = unique(toy$Criterion)[1:2],
    Score = c(2, 3)
  )

  list(
    unit_prediction = predict_mfrm_units(
      fit,
      new_units,
      n_draws = 2,
      seed = 1
    ),
    plausible_values = sample_mfrm_plausible_values(
      fit,
      new_units,
      n_draws = 2,
      seed = 1
    )
  )
})

summary_table_bundle_weighting_fixture <- local({
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  rasch_fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 7,
    maxit = 25
  ))
  gpcm_fit <- suppressWarnings(fit_mfrm(
    toy,
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

  build_weighting_audit(rasch_fit, gpcm_fit, theta_points = 21, top_n = 5)
})

summary_table_bundle_casebook_fixture <- local({
  toy <- load_mfrmr_data("example_core")
  keep_people <- unique(toy$Person)[1:12]
  toy <- toy[toy$Person %in% keep_people, , drop = FALSE]

  fit <- suppressWarnings(fit_mfrm(
    toy,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "MML",
    model = "RSM",
    quad_points = 7,
    maxit = 25
  ))
  diag <- suppressWarnings(diagnose_mfrm(
    fit,
    diagnostic_mode = "both",
    residual_pca = "none"
  ))

  build_misfit_casebook(fit, diagnostics = diag, top_n = 5)
})

expect_summary_bundle_roles_registered <- function(...) {
  bundles <- list(...)
  registry <- mfrmr:::summary_table_bundle_appendix_role_registry()
  expect_equal(nrow(registry), length(unique(registry$Role)))
  expect_equal(nrow(registry), length(unique(registry$PreferredAppendixOrder)))
  roles <- unique(unlist(lapply(bundles, function(bundle) {
    if (!inherits(bundle, "mfrm_summary_table_bundle")) return(character(0))
    idx <- as.data.frame(bundle$table_index, stringsAsFactors = FALSE)
    if (!"Role" %in% names(idx)) return(character(0))
    as.character(idx$Role)
  }), use.names = FALSE))
  roles <- roles[nzchar(roles)]
  missing_registry_roles <- setdiff(roles, registry$Role)
  expect_length(missing_registry_roles, 0L)
}

test_that("summary table bundle role registry covers every supported spec role", {
  df <- data.frame(stringsAsFactors = FALSE)
  minimal_summaries <- list(
    summary.mfrm_fit = list(
      overview = df,
      population_overview = df,
      population_design = df,
      population_coefficients = df,
      population_coding = df,
      facet_overview = df,
      person_overview = df,
      step_overview = df,
      slope_overview = df,
      settings_overview = df,
      reporting_map = df,
      caveats = df,
      facet_extremes = df,
      person_high = df,
      person_low = df
    ),
    summary.mfrm_diagnostics = list(
      overview = df,
      overall_fit = df,
      precision_profile = df,
      precision_audit = df,
      flags = df,
      reliability = df,
      top_fit = df,
      reporting_map = df
    ),
    summary.mfrm_data_description = list(
      overview = df,
      missing = df,
      score_distribution = df,
      facet_overview = df,
      agreement = df,
      reporting_map = df,
      caveats = df
    ),
    summary.mfrm_reporting_checklist = list(
      overview = df,
      section_summary = df,
      action_items = df,
      settings = list()
    ),
    summary.mfrm_apa_outputs = list(
      overview = df,
      components = df,
      sections = df,
      content_checks = df,
      preview = df
    ),
    summary.mfrm_design_evaluation = list(
      overview = df,
      design_summary = df
    ),
    summary.mfrm_signal_detection = list(
      overview = df,
      detection_summary = df
    ),
    summary.mfrm_population_prediction = list(
      design = df,
      overview = df,
      forecast = df
    ),
    summary.mfrm_future_branch_active_branch = list(
      overview = df,
      profile_summary = df,
      load_balance_summary = df,
      coverage_summary = df,
      guardrail_summary = df,
      readiness = df,
      recommendation_table = df,
      appendix_presets = df,
      appendix_role_summary = df,
      appendix_section_summary = df,
      selection_table_preset_summary = df,
      selection_handoff_table_summary = df,
      selection_handoff_preset_summary = df,
      selection_handoff_summary = df,
      selection_handoff_bundle_summary = df,
      selection_handoff_role_summary = df,
      selection_handoff_role_section_summary = df,
      selection_table_summary = df,
      selection_summary = df,
      selection_role_summary = df,
      selection_section_summary = df,
      selection_catalog = df,
      reporting_map = df
    ),
    summary.mfrm_facets_run = list(
      overview = df,
      mapping = df,
      run_info = df,
      fit = list(overview = df, reporting_map = df),
      diagnostics = list(overview = df, flags = df, reporting_map = df)
    ),
    summary.mfrm_bias = list(
      overview = df,
      chi_sq = df,
      final_iteration = df,
      top_rows = df,
      notes = character()
    ),
    summary.mfrm_anchor_audit = list(
      issue_counts = df,
      facet_summary = df,
      level_observation_summary = df,
      category_counts = df,
      recommendations = character(),
      notes = character()
    ),
    summary.mfrm_linking_review = list(
      overview = df,
      status = df,
      top_linking_risks = df,
      group_view_index = df,
      prefit_anchor_risks = df,
      drift_risks = df,
      chain_risks = df,
      plot_map = df,
      reporting_map = df,
      support_status = df,
      next_actions = character(),
      notes = character(),
      settings = list()
    ),
    summary.mfrm_misfit_casebook = list(
      overview = df,
      status = df,
      top_cases = df,
      case_rollup = df,
      group_view_index = df,
      source_summary = df,
      plot_map = df,
      reporting_map = df,
      support_status = df,
      key_warnings = character(),
      next_actions = character(),
      notes = character(),
      settings = list()
    ),
    summary.mfrm_weighting_audit = list(
      overview = df,
      status = df,
      top_measure_shifts = df,
      top_reweighted_levels = df,
      plot_map = df,
      reporting_map = df,
      support_status = df,
      key_warnings = character(),
      next_actions = character(),
      notes = character(),
      settings = list()
    ),
    summary.mfrm_unit_prediction = list(
      estimates = df,
      audit = df,
      population_audit = df,
      settings = list(),
      notes = character()
    ),
    summary.mfrm_plausible_values = list(
      draw_summary = df,
      estimates = df,
      audit = df,
      population_audit = df,
      settings = list(),
      notes = character()
    )
  )

  emitted_roles <- character()
  for (cls in names(minimal_summaries)) {
    obj <- minimal_summaries[[cls]]
    class(obj) <- cls
    emitted_roles <- c(
      emitted_roles,
      unname(mfrmr:::summary_table_bundle_spec(obj)$roles)
    )
  }
  emitted_roles <- unique(emitted_roles)
  registry <- mfrmr:::summary_table_bundle_appendix_role_registry()

  expect_setequal(registry$Role, emitted_roles)
  expect_equal(nrow(registry), length(unique(registry$Role)))
  expect_equal(nrow(registry), length(unique(registry$PreferredAppendixOrder)))
})

test_that("build_summary_table_bundle converts supported reporting summaries into named tables", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(fit_mfrm(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  ds <- describe_mfrm_data(toy, "Person", c("Rater", "Criterion"), "Score")
  chk <- reporting_checklist(fit, diagnostics = diag)
  apa <- build_apa_outputs(fit, diagnostics = diag)

  fit_bundle <- build_summary_table_bundle(fit)
  expect_s3_class(fit_bundle, "mfrm_summary_table_bundle")
  expect_identical(fit_bundle$source_class, "mfrm_fit")
  expect_identical(fit_bundle$summary_class, "summary.mfrm_fit")
  expect_true(all(c("overview", "facet_overview", "reporting_map") %in% names(fit_bundle$tables)))
  expect_true(all(c("Table", "Rows", "Cols", "Role", "Description") %in% names(fit_bundle$table_index)))
  expect_true(all(c("Table", "PlotReady", "NumericColumns", "DefaultPlotTypes") %in% names(fit_bundle$plot_index)))
  printed <- capture.output(print(fit_bundle))
  expect_true(any(grepl("mfrmr Summary Table Bundle", printed, fixed = TRUE)))

  fit_bundle_summary <- summary(fit_bundle)
  expect_s3_class(fit_bundle_summary, "summary.mfrm_summary_table_bundle")
  expect_true(is.data.frame(fit_bundle_summary$overview))
  expect_true(is.data.frame(fit_bundle_summary$table_catalog))
  expect_true(is.data.frame(fit_bundle_summary$table_profile))
  expect_true(is.data.frame(fit_bundle_summary$plot_index))
  expect_true(is.data.frame(fit_bundle_summary$appendix_presets))
  expect_true(is.data.frame(fit_bundle_summary$appendix_role_summary))
  expect_true(is.data.frame(fit_bundle_summary$appendix_section_summary))
  expect_true(is.data.frame(fit_bundle_summary$reporting_map))
  expect_true("AnyNumericTable" %in% names(fit_bundle_summary$overview))
  expect_true(all(c("RecommendedAppendixTables", "CompactAppendixTables") %in%
                    names(fit_bundle_summary$overview)))
  expect_true(all(c("Table", "ExportReady", "ApaTableReady", "RecommendedBridge") %in%
                    names(fit_bundle_summary$table_catalog)))
  expect_true(all(c("AppendixSection", "RecommendedAppendix", "CompactAppendix", "PreferredAppendixOrder", "AppendixRationale") %in%
                    names(fit_bundle_summary$table_catalog)))
  expect_identical(
    as.character(fit_bundle_summary$appendix_presets$Preset),
    c("all", "recommended", "compact", "methods", "results", "diagnostics", "reporting")
  )
  expect_true(all(c("Role", "Tables", "RecommendedTables", "CompactTables") %in%
                    names(fit_bundle_summary$appendix_role_summary)))
  expect_true(all(c("AppendixSection", "Tables", "RolesCovered") %in%
                    names(fit_bundle_summary$appendix_section_summary)))
  expect_true(all(c("Area", "CoveredHere", "CompanionOutput") %in%
                    names(fit_bundle_summary$reporting_map)))

  diag_bundle <- build_summary_table_bundle(diag, which = c("overview", "flags"))
  expect_identical(names(diag_bundle$tables), c("overview", "flags"))
  expect_identical(diag_bundle$source_class, "mfrm_diagnostics")

  ds_bundle <- build_summary_table_bundle(summary(ds))
  expect_identical(ds_bundle$source_class, "summary.mfrm_data_description")
  expect_true(all(c("overview", "missing", "score_distribution") %in% names(ds_bundle$tables)))

  chk_bundle <- build_summary_table_bundle(chk)
  expect_true(all(c("overview", "section_summary", "action_items") %in% names(chk_bundle$tables)))

  apa_bundle <- build_summary_table_bundle(apa, which = c("overview", "components", "preview"))
  expect_identical(names(apa_bundle$tables), c("overview", "components", "preview"))
  expect_summary_bundle_roles_registered(fit_bundle, diag_bundle, ds_bundle, chk_bundle, apa_bundle)
})

test_that("build_summary_table_bundle carries fit and score-support caveats", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d_gap <- d |> dplyr::filter(Score %in% 2:5)
  fit <- suppressWarnings(fit_mfrm(
    d_gap, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30,
    rating_min = 1,
    rating_max = 5
  ))
  ds <- describe_mfrm_data(
    d_gap, "Person", c("Rater", "Task", "Criterion"), "Score",
    rating_min = 1,
    rating_max = 5
  )
  fit_summary <- summary(fit)
  ds_summary <- summary(ds)

  fit_bundle <- build_summary_table_bundle(fit_summary)
  ds_bundle <- build_summary_table_bundle(ds_summary)
  fit_rec_bundle <- build_summary_table_bundle(fit_summary, appendix_preset = "recommended")
  fit_diag_bundle <- build_summary_table_bundle(fit_summary, appendix_preset = "diagnostics")
  ds_diag_bundle <- build_summary_table_bundle(ds_summary, appendix_preset = "diagnostics")
  fit_diag_catalog <- mfrmr:::summary_table_bundle_catalog(fit_diag_bundle)
  ds_diag_catalog <- mfrmr:::summary_table_bundle_catalog(ds_diag_bundle)

  expect_true("caveats" %in% names(fit_bundle$tables))
  expect_true("caveats" %in% names(ds_bundle$tables))
  expect_true("caveats" %in% names(fit_rec_bundle$tables))
  expect_true("caveats" %in% names(fit_diag_bundle$tables))
  expect_true("caveats" %in% names(ds_diag_bundle$tables))
  expect_identical(as.character(fit_bundle$table_index$Role[fit_bundle$table_index$Table == "caveats"]), "analysis_caveats")
  expect_identical(as.character(ds_bundle$table_index$Role[ds_bundle$table_index$Table == "caveats"]), "score_category_caveats")
  expect_identical(as.character(fit_diag_catalog$AppendixSection[fit_diag_catalog$Table == "caveats"]), "diagnostics")
  expect_identical(as.character(ds_diag_catalog$AppendixSection[ds_diag_catalog$Table == "caveats"]), "diagnostics")
  expect_true(any(grepl("Unused boundary score", fit_bundle$tables$caveats$Message, fixed = TRUE)))
  expect_true(any(grepl("prepared score support", ds_bundle$tables$caveats$Message, fixed = TRUE)))
  expect_true(any(grepl("Caveats", capture.output(print(fit_summary)), fixed = TRUE)))
  expect_true(any(grepl("Caveats", capture.output(print(ds_summary)), fixed = TRUE)))
  expect_summary_bundle_roles_registered(fit_bundle, ds_bundle, fit_rec_bundle, fit_diag_bundle, ds_diag_bundle)
})

test_that("build_summary_table_bundle exports population-model caveats as analysis caveats", {
  population_caveats <- mfrmr:::collect_mfrm_population_caveats(
    population_overview = data.frame(
      PopulationModel = TRUE,
      ResidualVariance = 0,
      OmittedPersons = 1L,
      OmittedRows = 4L,
      stringsAsFactors = FALSE
    ),
    population_design = data.frame(
      Column = c("(Intercept)", "X_zero", "X_incomplete"),
      IsIntercept = c(TRUE, FALSE, FALSE),
      ZeroVariance = c(TRUE, TRUE, FALSE),
      Complete = c(TRUE, TRUE, FALSE),
      stringsAsFactors = FALSE
    ),
    population_coefficients = data.frame()
  )
  fit_summary <- structure(
    list(
      overview = data.frame(),
      population_overview = data.frame(),
      population_design = data.frame(),
      population_coefficients = data.frame(),
      population_coding = data.frame(),
      facet_overview = data.frame(),
      person_overview = data.frame(),
      step_overview = data.frame(),
      slope_overview = data.frame(),
      settings_overview = data.frame(),
      reporting_map = data.frame(),
      caveats = population_caveats,
      facet_extremes = data.frame(),
      person_high = data.frame(),
      person_low = data.frame(),
      notes = population_caveats$Message
    ),
    class = "summary.mfrm_fit"
  )

  bundle <- build_summary_table_bundle(fit_summary)
  expect_true("caveats" %in% names(bundle$tables))
  expect_identical(as.character(bundle$table_index$Role[bundle$table_index$Table == "caveats"]), "analysis_caveats")
  expect_setequal(bundle$tables$caveats$Condition, population_caveats$Condition)
  expect_summary_bundle_roles_registered(bundle)
})

test_that("build_summary_table_bundle supports workflow, bias, anchor, linking, and prediction summaries", {
  run_bundle <- build_summary_table_bundle(summary_table_bundle_workflow_fixture$run)
  expect_identical(run_bundle$source_class, "mfrm_facets_run")
  expect_identical(run_bundle$summary_class, "summary.mfrm_facets_run")
  expect_true(all(c(
    "overview",
    "mapping",
    "run_info",
    "fit_overview",
    "diagnostic_flags"
  ) %in% names(run_bundle$tables)))

  bias_bundle <- build_summary_table_bundle(summary(summary_table_bundle_workflow_fixture$bias))
  expect_identical(bias_bundle$source_class, "summary.mfrm_bias")
  expect_identical(bias_bundle$summary_class, "summary.mfrm_bias")
  expect_true(all(c("overview", "chi_sq", "top_rows", "notes") %in% names(bias_bundle$tables)))

  audit_bundle <- build_summary_table_bundle(summary_table_bundle_workflow_fixture$audit)
  expect_identical(audit_bundle$source_class, "mfrm_anchor_audit")
  expect_identical(audit_bundle$summary_class, "summary.mfrm_anchor_audit")
  expect_true(all(c(
    "overview",
    "issue_counts",
    "facet_summary",
    "recommendations"
  ) %in% names(audit_bundle$tables)))

  linking_bundle <- build_summary_table_bundle(summary_table_bundle_workflow_fixture$linking_review)
  expect_identical(linking_bundle$source_class, "mfrm_linking_review")
  expect_identical(linking_bundle$summary_class, "summary.mfrm_linking_review")
  expect_true(all(c(
    "overview",
    "status",
    "top_linking_risks",
    "group_view_index",
    "plot_map",
    "reporting_map"
  ) %in% names(linking_bundle$tables)))

  unit_bundle <- build_summary_table_bundle(summary_table_bundle_prediction_fixture$unit_prediction)
  expect_identical(unit_bundle$source_class, "mfrm_unit_prediction")
  expect_identical(unit_bundle$summary_class, "summary.mfrm_unit_prediction")
  expect_true(all(c("overview", "estimates", "settings", "notes") %in% names(unit_bundle$tables)))

  pv_bundle <- build_summary_table_bundle(summary(summary_table_bundle_prediction_fixture$plausible_values))
  expect_identical(pv_bundle$source_class, "summary.mfrm_plausible_values")
  expect_identical(pv_bundle$summary_class, "summary.mfrm_plausible_values")
  expect_true(all(c(
    "overview",
    "draw_summary",
    "estimates",
    "settings",
    "notes"
  ) %in% names(pv_bundle$tables)))

  weighting_bundle <- build_summary_table_bundle(summary_table_bundle_weighting_fixture)
  expect_identical(weighting_bundle$source_class, "mfrm_weighting_audit")
  expect_identical(weighting_bundle$summary_class, "summary.mfrm_weighting_audit")
  expect_true(all(c(
    "overview",
    "status",
    "top_measure_shifts",
    "top_reweighted_levels",
    "plot_map",
    "reporting_map",
    "support_status"
  ) %in% names(weighting_bundle$tables)))

  casebook_bundle <- build_summary_table_bundle(summary_table_bundle_casebook_fixture)
  expect_identical(casebook_bundle$source_class, "mfrm_misfit_casebook")
  expect_identical(casebook_bundle$summary_class, "summary.mfrm_misfit_casebook")
  expect_true(all(c(
    "overview",
    "status",
    "top_cases",
    "case_rollup",
    "group_view_index",
    "source_summary",
    "plot_map",
    "reporting_map",
    "support_status"
  ) %in% names(casebook_bundle$tables)))
  expect_summary_bundle_roles_registered(
    run_bundle,
    bias_bundle,
    audit_bundle,
    linking_bundle,
    unit_bundle,
    pv_bundle,
    weighting_bundle,
    casebook_bundle
  )
})

test_that("build_summary_table_bundle keeps explicitly requested empty tables and rejects unknown names", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(fit_mfrm(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  ))

  empty_bundle <- build_summary_table_bundle(fit, which = "population_coefficients")
  expect_identical(names(empty_bundle$tables), "population_coefficients")
  expect_true(is.data.frame(empty_bundle$tables$population_coefficients))
  expect_equal(nrow(empty_bundle$table_index), 1L)

  empty_design_bundle <- build_summary_table_bundle(fit, which = "population_design")
  expect_identical(names(empty_design_bundle$tables), "population_design")
  expect_true(is.data.frame(empty_design_bundle$tables$population_design))
  expect_equal(nrow(empty_design_bundle$table_index), 1L)
  expect_identical(as.character(empty_design_bundle$table_index$Role[1]), "population_design")
  empty_design_summary <- summary(empty_design_bundle)
  expect_identical(as.character(empty_design_summary$table_catalog$AppendixSection[1]), "methods")
  expect_true(isTRUE(empty_design_summary$table_catalog$RecommendedAppendix[1]))
  expect_false(isTRUE(empty_design_summary$table_catalog$CompactAppendix[1]))

  empty_coding_bundle <- build_summary_table_bundle(fit, which = "population_coding")
  expect_identical(names(empty_coding_bundle$tables), "population_coding")
  expect_true(is.data.frame(empty_coding_bundle$tables$population_coding))
  expect_equal(nrow(empty_coding_bundle$table_index), 1L)
  expect_identical(as.character(empty_coding_bundle$table_index$Role[1]), "population_coding")
  empty_coding_summary <- summary(empty_coding_bundle)
  expect_identical(as.character(empty_coding_summary$table_catalog$AppendixSection[1]), "methods")
  expect_true(isTRUE(empty_coding_summary$table_catalog$RecommendedAppendix[1]))
  expect_false(isTRUE(empty_coding_summary$table_catalog$CompactAppendix[1]))

  expect_error(
    build_summary_table_bundle(fit, which = "not_a_table"),
    "received unknown `which` table name"
  )

  expect_error(
    build_summary_table_bundle(fit, which = "overview", appendix_preset = "recommended"),
    "requires `appendix_preset` and `which` to be used separately"
  )
})

test_that("build_summary_table_bundle validates front-door inputs before bundle conversion", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(fit_mfrm(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  ))

  expect_error(
    build_summary_table_bundle(NULL),
    "requires `x` to be a supported package object"
  )

  expect_error(
    build_summary_table_bundle(fit, which = 1),
    "requires `which` to be `NULL` or a non-empty character vector"
  )

  expect_error(
    build_summary_table_bundle(fit, include_empty = NA),
    "requires `include_empty` to be either `TRUE` or `FALSE`"
  )

  expect_error(
    build_summary_table_bundle(fit, digits = -1),
    "requires `digits` to be a single non-negative number"
  )

  broken_fit_summary <- structure(
    list(
      overview = data.frame(
        Model = "RSM",
        Method = "JML",
        stringsAsFactors = FALSE
      )
    ),
    class = "summary.mfrm_fit"
  )
  expect_error(
    build_summary_table_bundle(broken_fit_summary),
    "Missing required component\\(s\\): reporting_map"
  )

  broken_apa_summary <- structure(
    list(
      overview = data.frame(Components = 1, stringsAsFactors = FALSE),
      components = data.frame(Component = "report_text", stringsAsFactors = FALSE)
    ),
    class = "summary.mfrm_apa_outputs"
  )
  expect_error(
    build_summary_table_bundle(broken_apa_summary),
    "Missing required component\\(s\\): preview"
  )
})

test_that("build_summary_table_bundle supports planning and forecast summaries with future-branch tables", {
  spec <- build_mfrm_sim_spec(
    n_person = 10,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  sim_eval <- suppressWarnings(evaluate_mfrm_design(
    sim_spec = spec,
    n_person = c(10, 12),
    reps = 1,
    maxit = 5,
    seed = 901
  ))
  sig_eval <- suppressWarnings(evaluate_mfrm_signal_detection(
    sim_spec = spec,
    n_person = 10,
    reps = 1,
    maxit = 5,
    bias_max_iter = 1,
    seed = 902
  ))
  pred <- suppressWarnings(predict_mfrm_population(
    sim_spec = spec,
    design = list(person = c(10, 12)),
    reps = 1,
    maxit = 5,
    seed = 903
  ))

  design_bundle <- build_summary_table_bundle(sim_eval)
  expect_identical(design_bundle$source_class, "mfrm_design_evaluation")
  expect_true(all(c("overview", "design_summary", "future_branch_overview",
                    "future_branch_recommendation") %in% names(design_bundle$tables)))

  signal_bundle <- build_summary_table_bundle(summary(sig_eval))
  expect_identical(signal_bundle$source_class, "summary.mfrm_signal_detection")
  expect_true(all(c("overview", "detection_summary", "future_branch_readiness") %in%
                    names(signal_bundle$tables)))

  pred_bundle <- build_summary_table_bundle(pred)
  expect_identical(pred_bundle$source_class, "mfrm_population_prediction")
  expect_true(all(c("design", "forecast", "future_branch_profile",
                    "future_branch_load_balance", "future_branch_coverage") %in%
                    names(pred_bundle$tables)))
  expect_summary_bundle_roles_registered(design_bundle, signal_bundle, pred_bundle)
})

test_that("build_summary_table_bundle supports future arbitrary-facet active-branch inputs", {
  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  active <- spec$planning_schema$future_branch_active_branch
  active_bundle <- build_summary_table_bundle(active)
  expect_s3_class(active_bundle, "mfrm_summary_table_bundle")
  expect_identical(active_bundle$source_class, "mfrm_future_branch_active_branch")
  expect_identical(active_bundle$summary_class, "summary.mfrm_future_branch_active_branch")
  expect_true(all(c(
    "future_branch_overview",
    "future_branch_profile",
    "future_branch_load_balance",
    "future_branch_coverage",
    "future_branch_guardrails",
    "future_branch_readiness",
    "future_branch_recommendation",
    "future_branch_appendix_presets",
    "future_branch_appendix_roles",
    "future_branch_appendix_sections",
    "future_branch_selection_table_presets",
    "future_branch_selection_handoff_tables",
    "future_branch_selection_handoff_presets",
    "future_branch_selection_handoff",
    "future_branch_selection_handoff_bundles",
    "future_branch_selection_handoff_roles",
    "future_branch_selection_handoff_role_sections",
    "future_branch_selection_tables",
    "future_branch_selection_summary",
    "future_branch_selection_roles",
    "future_branch_selection_sections",
    "future_branch_selection_catalog",
    "future_branch_reporting_map"
  ) %in% names(active_bundle$tables)))

  summary_bundle <- build_summary_table_bundle(summary(active))
  expect_identical(summary_bundle$source_class, "summary.mfrm_future_branch_active_branch")
  expect_identical(summary_bundle$summary_class, "summary.mfrm_future_branch_active_branch")
  expect_true(all(c(
    "future_branch_overview",
    "future_branch_profile",
    "future_branch_selection_table_presets",
    "future_branch_selection_handoff_tables",
    "future_branch_selection_handoff_presets",
    "future_branch_selection_handoff",
    "future_branch_selection_handoff_bundles",
    "future_branch_selection_handoff_roles",
    "future_branch_selection_handoff_role_sections",
    "future_branch_selection_tables",
    "future_branch_recommendation",
    "future_branch_selection_summary",
    "future_branch_reporting_map"
  ) %in% names(summary_bundle$tables)))

  active_bundle_summary <- summary(active_bundle)
  expect_s3_class(active_bundle_summary, "summary.mfrm_summary_table_bundle")
  expect_true(is.data.frame(active_bundle_summary$selection_handoff_table_summary))
  expect_true(is.data.frame(active_bundle_summary$selection_handoff_preset_summary))
  expect_true(is.data.frame(active_bundle_summary$selection_handoff_summary))
  expect_true(is.data.frame(active_bundle_summary$selection_handoff_bundle_summary))
  expect_true(is.data.frame(active_bundle_summary$selection_handoff_role_summary))
  expect_true(is.data.frame(active_bundle_summary$selection_handoff_role_section_summary))
  expect_true(is.data.frame(active_bundle_summary$selection_table_summary))
  expect_true(is.data.frame(active_bundle_summary$selection_role_summary))
  expect_true(is.data.frame(active_bundle_summary$selection_section_summary))
  expect_true(all(c("Preset", "SectionsCovered", "PlotReadyTables", "PlotReadyFraction", "NumericFraction") %in%
                    names(active_bundle_summary$selection_handoff_preset_summary)))
  expect_true(all(c("Preset", "AppendixSection", "PlotReadyTables", "PlotReadyFraction", "NumericFraction") %in%
                    names(active_bundle_summary$selection_handoff_summary)))
  expect_true(all(c("Preset", "AppendixSection", "Bundle", "PlotReadyTables", "PlotReadyFraction", "NumericFraction") %in%
                    names(active_bundle_summary$selection_handoff_bundle_summary)))
  expect_true(all(c("Preset", "Role", "PlotReadyTables", "PlotReadyFraction", "NumericFraction") %in%
                    names(active_bundle_summary$selection_handoff_role_summary)))
  expect_true(all(c("Preset", "AppendixSection", "Role", "PlotReadyTables", "PlotReadyFraction", "NumericFraction") %in%
                    names(active_bundle_summary$selection_handoff_role_section_summary)))
  expect_true(all(c("Preset", "AppendixSection", "Role", "Bundle", "Table", "Rows", "NumericColumns", "PlotReady", "ExportReady", "ApaTableReady") %in%
                    names(active_bundle_summary$selection_handoff_table_summary)))
  expect_true(all(c("Preset", "Bundle", "TablesAvailable", "SelectionFraction", "PlotReadyFraction", "NumericFraction") %in%
                    names(active_bundle_summary$selection_summary)))
  expect_true(all(c("Preset", "Role", "PlotReadyFraction", "NumericFraction") %in%
                    names(active_bundle_summary$selection_role_summary)))
  expect_true(all(c("Preset", "AppendixSection", "PlotReadyFraction", "NumericFraction") %in%
                    names(active_bundle_summary$selection_section_summary)))
  expect_summary_bundle_roles_registered(active_bundle, summary_bundle)
})

test_that("future arbitrary-facet active-branch bundles support appendix presets", {
  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  active <- spec$planning_schema$future_branch_active_branch
  full_bundle <- build_summary_table_bundle(active)
  rec_bundle <- build_summary_table_bundle(active, appendix_preset = "recommended")
  compact_bundle <- build_summary_table_bundle(active, appendix_preset = "compact")
  methods_bundle <- build_summary_table_bundle(active, appendix_preset = "methods")
  diagnostics_bundle <- build_summary_table_bundle(active, appendix_preset = "diagnostics")

  expect_true(all(c(
    "future_branch_overview",
    "future_branch_profile",
    "future_branch_readiness",
    "future_branch_recommendation"
  ) %in% names(rec_bundle$tables)))
  expect_false("future_branch_selection_table_presets" %in% names(rec_bundle$tables))
  expect_false("future_branch_selection_handoff" %in% names(rec_bundle$tables))
  expect_false("future_branch_selection_summary" %in% names(rec_bundle$tables))
  expect_false("future_branch_load_balance" %in% names(rec_bundle$tables))
  expect_false("future_branch_coverage" %in% names(rec_bundle$tables))
  expect_false("future_branch_guardrails" %in% names(rec_bundle$tables))

  expect_true(all(c(
    "future_branch_overview",
    "future_branch_readiness",
    "future_branch_recommendation"
  ) %in% names(compact_bundle$tables)))
  expect_false("future_branch_selection_table_presets" %in% names(compact_bundle$tables))
  expect_false("future_branch_selection_handoff" %in% names(compact_bundle$tables))
  expect_false("future_branch_appendix_presets" %in% names(compact_bundle$tables))
  expect_false("future_branch_profile" %in% names(compact_bundle$tables))

  expect_true(all(methods_bundle$table_index$AppendixSection %in% "methods"))
  expect_true(all(diagnostics_bundle$table_index$AppendixSection %in% "diagnostics"))
  expect_true(nrow(compact_bundle$table_index) <= nrow(rec_bundle$table_index))
  expect_true(nrow(rec_bundle$table_index) <= nrow(full_bundle$table_index))
})

test_that("build_summary_table_bundle applies appendix presets at bundle-construction time", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(fit_mfrm(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  ))

  full_bundle <- build_summary_table_bundle(fit)
  rec_bundle <- build_summary_table_bundle(fit, appendix_preset = "recommended")
  compact_bundle <- build_summary_table_bundle(fit, appendix_preset = "compact")
  methods_bundle <- build_summary_table_bundle(fit, appendix_preset = "methods")
  results_bundle <- build_summary_table_bundle(fit, appendix_preset = "results")

  expect_identical(rec_bundle$appendix_preset, "recommended")
  expect_identical(compact_bundle$appendix_preset, "compact")
  expect_identical(methods_bundle$appendix_preset, "methods")
  expect_identical(results_bundle$appendix_preset, "results")
  expect_true("AppendixPreset" %in% names(rec_bundle$overview))
  expect_identical(as.character(rec_bundle$overview$AppendixPreset[1]), "recommended")
  expect_identical(as.character(summary(rec_bundle)$overview$AppendixPreset[1]), "recommended")
  expect_true(nrow(rec_bundle$table_index) <= nrow(full_bundle$table_index))
  expect_true(nrow(compact_bundle$table_index) <= nrow(rec_bundle$table_index))
  expect_true(all(rec_bundle$table_index$Table %in% full_bundle$table_index$Table))
  expect_true(all(compact_bundle$table_index$Table %in% rec_bundle$table_index$Table))
  expect_true(all(methods_bundle$table_index$AppendixSection %in% "methods"))
  expect_true(all(results_bundle$table_index$AppendixSection %in% "results"))
  expect_false("reporting_map" %in% names(rec_bundle$tables))
})

test_that("apa_table consumes summary outputs and summary table bundles directly", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(fit_mfrm(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  chk <- reporting_checklist(fit, diagnostics = diag)

  tbl_from_summary <- apa_table(summary(fit), which = "facet_overview")
  expect_s3_class(tbl_from_summary, "apa_table")
  expect_identical(tbl_from_summary$which, "facet_overview")
  expect_true(nrow(tbl_from_summary$table) > 0)
  expect_true(nzchar(tbl_from_summary$caption))

  bundle <- build_summary_table_bundle(chk)
  tbl_from_bundle <- apa_table(bundle, which = "section_summary")
  expect_s3_class(tbl_from_bundle, "apa_table")
  expect_identical(tbl_from_bundle$which, "section_summary")
  expect_true(nrow(tbl_from_bundle$table) > 0)
  expect_true(nzchar(tbl_from_bundle$note))

  expect_error(
    apa_table(bundle, which = "not_present"),
    "Requested `which` not found in summary table bundle"
  )
})

test_that("apa_table consumes workflow, bias, anchor, and prediction summaries directly", {
  workflow_tbl <- apa_table(summary(summary_table_bundle_workflow_fixture$run), which = "mapping")
  expect_s3_class(workflow_tbl, "apa_table")
  expect_identical(workflow_tbl$which, "mapping")
  expect_true(nrow(workflow_tbl$table) > 0L)

  bias_tbl <- apa_table(summary(summary_table_bundle_workflow_fixture$bias), which = "top_rows")
  expect_s3_class(bias_tbl, "apa_table")
  expect_identical(bias_tbl$which, "top_rows")
  expect_true(nrow(bias_tbl$table) > 0L)

  anchor_tbl <- apa_table(summary(summary_table_bundle_workflow_fixture$audit), which = "facet_summary")
  expect_s3_class(anchor_tbl, "apa_table")
  expect_identical(anchor_tbl$which, "facet_summary")
  expect_true(nrow(anchor_tbl$table) > 0L)

  unit_tbl <- apa_table(summary(summary_table_bundle_prediction_fixture$unit_prediction), which = "estimates")
  expect_s3_class(unit_tbl, "apa_table")
  expect_identical(unit_tbl$which, "estimates")
  expect_true(nrow(unit_tbl$table) > 0L)

  pv_tbl <- apa_table(summary(summary_table_bundle_prediction_fixture$plausible_values), which = "draw_summary")
  expect_s3_class(pv_tbl, "apa_table")
  expect_identical(pv_tbl$which, "draw_summary")
  expect_true(nrow(pv_tbl$table) > 0L)

  weighting_tbl <- apa_table(summary(summary_table_bundle_weighting_fixture), which = "top_reweighted_levels")
  expect_s3_class(weighting_tbl, "apa_table")
  expect_identical(weighting_tbl$which, "top_reweighted_levels")
  expect_true(nrow(weighting_tbl$table) > 0L)
})

test_that("apa_table consumes future arbitrary-facet active-branch summaries directly", {
  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  active <- spec$planning_schema$future_branch_active_branch
  tbl <- apa_table(summary(active), which = "future_branch_readiness")
  expect_s3_class(tbl, "apa_table")
  expect_identical(tbl$which, "future_branch_readiness")
  expect_true(nrow(tbl$table) > 0L)
  expect_true(nzchar(tbl$caption))
})

test_that("plot methods consume summary table bundles directly", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(fit_mfrm(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  ))

  bundle <- build_summary_table_bundle(fit)
  rows_plot <- plot(bundle, type = "table_rows", draw = FALSE)
  expect_s3_class(rows_plot, "mfrm_plot_data")
  expect_identical(rows_plot$name, "summary_table_bundle")
  expect_identical(rows_plot$data$plot, "table_rows")

  roles_plot <- plot(bundle, type = "role_tables", draw = FALSE)
  expect_s3_class(roles_plot, "mfrm_plot_data")
  expect_identical(roles_plot$data$plot, "role_tables")

  appendix_roles_plot <- plot(bundle, type = "appendix_roles", draw = FALSE)
  expect_s3_class(appendix_roles_plot, "mfrm_plot_data")
  expect_identical(appendix_roles_plot$data$plot, "appendix_roles")

  sections_plot <- plot(bundle, type = "appendix_sections", draw = FALSE)
  expect_s3_class(sections_plot, "mfrm_plot_data")
  expect_identical(sections_plot$data$plot, "appendix_sections")

  presets_plot <- plot(bundle, type = "appendix_presets", draw = FALSE)
  expect_s3_class(presets_plot, "mfrm_plot_data")
  expect_identical(presets_plot$data$plot, "appendix_presets")

  numeric_plot <- plot(bundle, type = "numeric_profile", which = "facet_overview", draw = FALSE)
  expect_s3_class(numeric_plot, "mfrm_plot_data")
  expect_identical(numeric_plot$name, "summary_table_bundle")
  expect_identical(numeric_plot$data$source_table, "facet_overview")

  first_numeric_plot <- plot(bundle, type = "first_numeric", which = "facet_overview", draw = FALSE)
  expect_s3_class(first_numeric_plot, "mfrm_plot_data")
  expect_identical(first_numeric_plot$data$source_table, "facet_overview")
})

test_that("future-branch summary table bundles expose selection plot surfaces", {
  spec <- build_mfrm_sim_spec(
    n_person = 12,
    n_rater = 3,
    n_criterion = 4,
    raters_per_person = 2,
    assignment = "rotating",
    facet_names = c("Judge", "Task")
  )

  active_bundle <- build_summary_table_bundle(spec$planning_schema$future_branch_active_branch)

  handoff_preset_plot <- plot(active_bundle, type = "selection_handoff_presets", appendix_preset = "all", draw = FALSE)
  expect_s3_class(handoff_preset_plot, "mfrm_plot_data")
  expect_identical(handoff_preset_plot$name, "summary_table_bundle")
  expect_identical(handoff_preset_plot$data$plot, "selection_handoff_presets")
  expect_identical(handoff_preset_plot$data$appendix_preset, "all")

  handoff_plot <- plot(active_bundle, type = "selection_handoff", appendix_preset = "recommended", draw = FALSE)
  expect_s3_class(handoff_plot, "mfrm_plot_data")
  expect_identical(handoff_plot$name, "summary_table_bundle")
  expect_identical(handoff_plot$data$plot, "selection_handoff")
  expect_identical(handoff_plot$data$appendix_preset, "recommended")

  handoff_fraction_plot <- plot(active_bundle, type = "selection_handoff", appendix_preset = "recommended", selection_value = "fraction", draw = FALSE)
  expect_s3_class(handoff_fraction_plot, "mfrm_plot_data")
  expect_identical(handoff_fraction_plot$data$plot, "selection_handoff")
  expect_identical(handoff_fraction_plot$data$selection_value, "fraction")

  handoff_bundle_plot <- plot(active_bundle, type = "selection_handoff_bundles", appendix_preset = "recommended", draw = FALSE)
  expect_s3_class(handoff_bundle_plot, "mfrm_plot_data")
  expect_identical(handoff_bundle_plot$name, "summary_table_bundle")
  expect_identical(handoff_bundle_plot$data$plot, "selection_handoff_bundles")
  expect_identical(handoff_bundle_plot$data$appendix_preset, "recommended")

  handoff_role_plot <- plot(active_bundle, type = "selection_handoff_roles", appendix_preset = "recommended", draw = FALSE)
  expect_s3_class(handoff_role_plot, "mfrm_plot_data")
  expect_identical(handoff_role_plot$name, "summary_table_bundle")
  expect_identical(handoff_role_plot$data$plot, "selection_handoff_roles")
  expect_identical(handoff_role_plot$data$appendix_preset, "recommended")

  handoff_role_section_plot <- plot(active_bundle, type = "selection_handoff_role_sections", appendix_preset = "recommended", draw = FALSE)
  expect_s3_class(handoff_role_section_plot, "mfrm_plot_data")
  expect_identical(handoff_role_section_plot$name, "summary_table_bundle")
  expect_identical(handoff_role_section_plot$data$plot, "selection_handoff_role_sections")
  expect_identical(handoff_role_section_plot$data$appendix_preset, "recommended")

  section_plot <- plot(active_bundle, type = "selection_sections", appendix_preset = "compact", draw = FALSE)
  expect_s3_class(section_plot, "mfrm_plot_data")
  expect_identical(section_plot$data$plot, "selection_sections")
  expect_identical(section_plot$data$appendix_preset, "compact")

  expect_error(
    plot(active_bundle, type = "selection_tables", appendix_preset = "recommended", selection_value = "fraction", draw = FALSE),
    "not available for `type = \"selection_tables\"`",
    fixed = TRUE
  )
})
