.make_disconnected_first_use_data <- function(seed = 20260724L) {
  full <- simulate_mfrm_data(
    n_person = 20,
    n_rater = 4,
    n_criterion = 4,
    raters_per_person = 4,
    assignment = "crossed",
    seed = seed
  )
  person_index <- suppressWarnings(as.integer(sub("^P", "", full$Person)))
  component_a <- person_index <= 10L &
    full$Rater %in% c("R01", "R02") &
    full$Criterion %in% c("C01", "C02")
  component_b <- person_index > 10L &
    full$Rater %in% c("R03", "R04") &
    full$Criterion %in% c("C03", "C04")
  full[component_a | component_b, , drop = FALSE]
}

.summary_status_value <- function(x, item) {
  unname(x$status$Value[match(item, x$status$Item)][1L])
}

test_that("disconnected designs remain reporting holds after numerical convergence", {
  dat <- .make_disconnected_first_use_data()
  expect_warning(
    fit <- fit_mfrm(
      dat,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      rating_min = 1,
      rating_max = 4,
      method = "MML"
    ),
    "disconnected subset"
  )

  expect_true(isTRUE(fit$summary$InferenceReady))
  expect_equal(fit$data_review$overall_connectivity$components, 2L)
  expect_identical(
    fit$data_review$status$Status[fit$data_review$status$Domain == "Design"],
    "hold_disconnected"
  )

  fit_summary <- summary(fit, profile = "fit", detail = "brief")
  expect_identical(
    .summary_status_value(fit_summary, "Overall status"),
    "review_required"
  )
  expect_identical(
    .summary_status_value(fit_summary, "Reporting readiness"),
    "hold_for_design_review"
  )
  expect_identical(
    fit_summary$required_visual$InterpretationStatus[1L],
    "review_only"
  )
  expect_false(fit_summary$required_visual$InterpretationReady[1L])

  expect_warning(
    wright <- plot(fit, type = "wright", top_n = Inf, draw = FALSE),
    "Review-only display"
  )
  expect_identical(wright$data$interpretation_status, "review_only")
  expect_match(wright$data$subtitle, "REVIEW ONLY", fixed = TRUE)
  expect_true(any(grepl("CI", wright$data$legend$label, fixed = TRUE)))
  expect_identical(
    mfrmr:::.mfrm_plot_display_title(
      wright$data,
      title = NULL,
      fallback = "Wright map"
    ),
    "REVIEW ONLY - Wright map"
  )

  expect_warning(
    unified <- plot_wright_unified(fit, draw = FALSE),
    "Review-only display"
  )
  expect_identical(unified$interpretation_status, "review_only")
  expect_identical(
    mfrmr:::.mfrm_plot_display_title(
      unified,
      title = unified$title,
      fallback = "Unified Wright Map"
    ),
    paste0("REVIEW ONLY - ", unified$title)
  )

  diagnostics <- diagnose_mfrm(
    fit,
    residual_pca = "none",
    diagnostic_mode = "legacy"
  )
  checklist <- reporting_checklist(fit, diagnostics = diagnostics)
  connectivity <- checklist$checklist[
    checklist$checklist$Item == "Connectivity assessed",
    ,
    drop = FALSE
  ]
  expect_true(connectivity$Available)
  expect_false(connectivity$DraftReady)
  expect_identical(connectivity$Severity, "required")
  expect_match(connectivity$Detail, "2 disconnected subsets", fixed = TRUE)
})

test_that("boundary separation is a stability hold, not a numerical pass", {
  dat <- simulate_mfrm_data(
    n_person = 18,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 4,
    assignment = "crossed",
    seed = 77
  )
  dat$Score[dat$Rater == "R01"] <- 1L
  dat$Score[dat$Rater == "R02"] <- 4L

  expect_warning(
    fit <- fit_mfrm(
      dat,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      rating_min = 1,
      rating_max = 4,
      method = "MML"
    ),
    "Boundary-constant"
  )

  expect_true(isTRUE(fit$summary$InferenceReady))
  expect_identical(
    fit$data_review$status$Status[
      fit$data_review$status$Domain == "Stability"
    ],
    "hold_boundary_separation"
  )
  expect_setequal(
    fit$data_review$boundary_levels$Level,
    c("R01", "R02")
  )
  expect_gt(
    max(abs(fit$facets$others$Estimate[fit$facets$others$Facet == "Rater"])),
    10
  )
  expect_identical(
    .summary_status_value(summary(fit), "Reporting readiness"),
    "hold_for_stability_review"
  )
  expect_warning(
    facets_map <- plot(
      fit,
      type = "wright",
      renderer = "facets",
      show_ci = FALSE,
      draw = FALSE
    ),
    "Review-only display"
  )
  expect_identical(
    facets_map$data$facets_style$settings$AutoRangePolicy,
    "boundary_levels_at_ends"
  )
  expect_lt(
    diff(range(facets_map$data$facets_style$ruler_rows$Logit)),
    10
  )
})

test_that("first-use controls and non-finite inputs fail before fitting", {
  dat <- load_mfrmr_data("example_core")
  common <- list(
    data = dat,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML"
  )

  expect_error(
    do.call(fit_mfrm, c(common, list(maxit = 5.5))),
    "finite positive integer"
  )
  expect_error(
    do.call(fit_mfrm, c(common, list(quad_points = 7.5))),
    "finite positive integer"
  )
  expect_error(
    do.call(fit_mfrm, c(common, list(attach_diagnostics = NA))),
    "single logical value"
  )

  score_inf <- dat
  score_inf$Score[1L] <- Inf
  expect_error(
    do.call(fit_mfrm, utils::modifyList(common, list(data = score_inf))),
    "finite ordered integer category codes"
  )
  score_nan <- dat
  score_nan$Score[1L] <- NaN
  expect_error(
    do.call(fit_mfrm, utils::modifyList(common, list(data = score_nan))),
    "finite ordered integer category codes"
  )

  weight_inf <- dat
  weight_inf$Weight <- 1
  weight_inf$Weight[1L] <- Inf
  expect_error(
    do.call(
      fit_mfrm,
      c(utils::modifyList(common, list(data = weight_inf)), list(weight = "Weight"))
    ),
    "finite numeric values"
  )
  weight_nan <- dat
  weight_nan$Weight <- 1
  weight_nan$Weight[1L] <- NaN
  expect_error(
    do.call(
      fit_mfrm,
      c(utils::modifyList(common, list(data = weight_nan)), list(weight = "Weight"))
    ),
    "finite numeric values"
  )

  blank_id <- dat
  blank_id$Rater[1L] <- "   "
  expect_error(
    do.call(fit_mfrm, utils::modifyList(common, list(data = blank_id))),
    "Blank person or facet identifier"
  )
})

test_that("duplicate cells warn once and propagate a data-review state", {
  dat <- simulate_mfrm_data(
    n_person = 12,
    n_rater = 3,
    n_criterion = 2,
    raters_per_person = 2,
    assignment = "rotating",
    seed = 191
  )
  dat <- rbind(dat, dat[1:2, , drop = FALSE])
  warnings <- character(0)
  fit <- withCallingHandlers(
    fit_mfrm(
      dat,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      method = "JML",
      maxit = 20
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  duplicate_warnings <- warnings[grepl("duplicated Person", warnings, fixed = TRUE)]
  expect_length(duplicate_warnings, 1L)
  expect_match(duplicate_warnings, "4 row(s) across 2", fixed = TRUE)
  expect_identical(
    fit$data_review$status$Status[fit$data_review$status$Domain == "Data"],
    "review"
  )
  notes <- fit$prep$preparation_notes
  expect_equal(sum(notes$Condition == "duplicate_person_facet_cells"), 1L)
})
