.results_readiness_base_data <- function(seed = 20260724L) {
  simulate_mfrm_data(
    n_person = 20,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 4,
    assignment = "crossed",
    seed = seed
  )
}

.results_readiness_person_index <- function(data) {
  suppressWarnings(as.integer(sub("^P", "", as.character(data$Person))))
}

.results_readiness_disconnected_data <- function() {
  data <- .results_readiness_base_data()
  person_index <- .results_readiness_person_index(data)
  component_a <- person_index <= 10L &
    data$Rater %in% c("R01", "R02") &
    data$Criterion %in% "C01"
  component_b <- person_index > 10L &
    data$Rater %in% c("R03", "R04") &
    data$Criterion %in% c("C02", "C03")
  data[component_a | component_b, , drop = FALSE]
}

.results_readiness_shared_criterion_data <- function() {
  data <- .results_readiness_base_data()
  person_index <- .results_readiness_person_index(data)
  component_a <- person_index <= 10L & data$Rater %in% c("R01", "R02")
  component_b <- person_index > 10L & data$Rater %in% c("R03", "R04")
  data[component_a | component_b, , drop = FALSE]
}

.results_readiness_fit <- function(data) {
  suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 1,
    rating_max = 4,
    method = "MML"
  ))
}

.results_readiness_value <- function(readiness, domain) {
  as.character(readiness$Status[match(domain, readiness$Domain)])[1L]
}

.results_readiness_status_value <- function(status, section) {
  as.character(status$Status[match(section, status$Section)])[1L]
}

test_that("mfrm_results propagates a disconnected-design hold", {
  fit <- .results_readiness_fit(.results_readiness_disconnected_data())
  expect_true(isTRUE(fit$summary$InferenceReady))
  expect_equal(fit$data_review$overall_connectivity$components, 2L)

  diagnostics <- suppressWarnings(diagnose_mfrm(
    fit,
    residual_pca = "none",
    diagnostic_mode = "legacy"
  ))
  res <- suppressWarnings(mfrm_results(
    fit,
    diagnostics = diagnostics,
    include = c("fit", "diagnostics", "reporting")
  ))

  expect_identical(.results_readiness_value(res$readiness, "Data"), "pass")
  expect_identical(
    .results_readiness_value(res$readiness, "Design"),
    "hold_disconnected"
  )
  expect_identical(.results_readiness_value(res$readiness, "Stability"), "pass")
  expect_identical(.results_readiness_value(res$readiness, "Plot"), "review_only")
  expect_identical(
    .results_readiness_value(res$readiness, "Reporting"),
    "hold_for_design_review"
  )

  # Existing computation-availability rows retain their original meaning.
  expect_identical(.results_readiness_status_value(res$status, "input"), "ok")
  expect_identical(.results_readiness_status_value(res$status, "fit_summary"), "ok")
  expect_identical(.results_readiness_status_value(res$status, "data_readiness"), "ok")
  expect_identical(.results_readiness_status_value(res$status, "design_readiness"), "review")
  expect_identical(.results_readiness_status_value(res$status, "stability_readiness"), "ok")
  expect_identical(.results_readiness_status_value(res$status, "plot_interpretation"), "review")
  expect_identical(.results_readiness_status_value(res$status, "reporting_readiness"), "review")

  wright <- res$plot_map[res$plot_map$Type %in% "wright", , drop = FALSE]
  expect_true(wright$Available)
  expect_true(wright$RequiredArtifact)
  expect_identical(wright$InterpretationStatus, "review_only")
  expect_false(wright$InterpretationReady)

  design_triage <- res$triage[
    res$triage$Area %in% "Design / connectivity",
    ,
    drop = FALSE
  ]
  data_triage <- res$triage[res$triage$Area %in% "Data review", , drop = FALSE]
  stability_triage <- res$triage[res$triage$Area %in% "Stability", , drop = FALSE]
  wright_triage <- res$triage[res$triage$Area %in% "Wright map", , drop = FALSE]
  reporting_triage <- res$triage[res$triage$Area %in% "Reporting", , drop = FALSE]
  expect_identical(data_triage$Severity, "ok")
  expect_identical(design_triage$Severity, "review")
  expect_match(design_triage$Detail, "2 disconnected subsets", fixed = TRUE)
  expect_identical(stability_triage$Severity, "ok")
  expect_identical(wright_triage$Severity, "review")
  expect_identical(wright_triage$Signal, "required_wright_map_review_only")
  expect_identical(reporting_triage$Severity, "review")
  expect_identical(
    reporting_triage$Signal,
    "reporting_checklist_available_but_held"
  )

  sx <- summary(res)
  expect_identical(.results_readiness_value(sx$readiness, "Design"), "hold_disconnected")
})

test_that("a shared Criterion links otherwise disjoint Rater groups", {
  data <- .results_readiness_shared_criterion_data()
  person_index <- .results_readiness_person_index(data)
  raters_a <- unique(data$Rater[person_index <= 10L])
  raters_b <- unique(data$Rater[person_index > 10L])
  criteria_a <- unique(data$Criterion[person_index <= 10L])
  criteria_b <- unique(data$Criterion[person_index > 10L])
  expect_length(intersect(raters_a, raters_b), 0L)
  expect_gt(length(intersect(criteria_a, criteria_b)), 0L)

  fit <- .results_readiness_fit(data)
  expect_true(isTRUE(fit$summary$InferenceReady))
  expect_equal(fit$data_review$overall_connectivity$components, 1L)
  res <- suppressWarnings(mfrm_results(
    fit,
    include = "fit",
    compute = "never"
  ))

  expect_identical(.results_readiness_value(res$readiness, "Design"), "pass_linked")
  expect_identical(
    .results_readiness_value(res$readiness, "Plot"),
    "ready_for_diagnostic_interpretation"
  )
  expect_identical(.results_readiness_status_value(res$status, "design_readiness"), "ok")
  expect_identical(.results_readiness_status_value(res$status, "plot_interpretation"), "ok")

  wright <- res$plot_map[res$plot_map$Type %in% "wright", , drop = FALSE]
  expect_true(wright$Available)
  expect_identical(
    wright$InterpretationStatus,
    "ready_for_diagnostic_interpretation"
  )
  expect_true(wright$InterpretationReady)
  expect_true(any(
    res$triage$Area %in% "Design / connectivity" &
      res$triage$Severity %in% "ok"
  ))
  expect_true(any(
    res$triage$Area %in% "Wright map" & res$triage$Severity %in% "ok"
  ))
})

test_that("fully connected results remain interpretation-ready", {
  fit <- .results_readiness_fit(.results_readiness_base_data())
  expect_true(isTRUE(fit$summary$InferenceReady))
  expect_equal(fit$data_review$overall_connectivity$components, 1L)
  res <- suppressWarnings(mfrm_results(
    fit,
    include = "fit",
    compute = "never"
  ))

  expect_identical(.results_readiness_value(res$readiness, "Data"), "pass")
  expect_identical(.results_readiness_value(res$readiness, "Design"), "pass_linked")
  expect_identical(.results_readiness_value(res$readiness, "Stability"), "pass")
  expect_identical(
    .results_readiness_value(res$readiness, "Plot"),
    "ready_for_diagnostic_interpretation"
  )
  expect_true(all(c(
    "InterpretationStatus", "InterpretationReady", "ReadinessRoute"
  ) %in% names(res$plot_map)))
  expect_true(res$plot_map$Available[res$plot_map$Type %in% "wright"])
  expect_true(res$plot_map$InterpretationReady[res$plot_map$Type %in% "wright"])

  # Summaries of older serialized result objects can reconstruct readiness.
  legacy <- res
  legacy$readiness <- NULL
  legacy_summary <- summary(legacy)
  expect_identical(
    .results_readiness_value(legacy_summary$readiness, "Design"),
    "pass_linked"
  )
})
