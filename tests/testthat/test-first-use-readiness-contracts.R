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

test_that("exactly aliased disconnected designs stop before optimization", {
  dat <- .make_disconnected_first_use_data()
  condition <- tryCatch(
    suppressWarnings(fit_mfrm(
      dat,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      rating_min = 1,
      rating_max = 4,
      method = "MML"
    )),
    error = identity
  )
  expect_s3_class(condition, "mfrmr_estimability_error")
  expect_s3_class(condition, "mfrmr_readiness_condition")
  expect_identical(
    condition$readiness$EstimabilityState,
    "structurally_unidentified"
  )
  expect_identical(condition$estimability$observed_components, 2L)
  expect_gt(condition$estimability$design$nullity, 0L)
  expect_match(condition$readiness$ReasonCodes,
               "disconnected_without_link", fixed = TRUE)
  expect_match(conditionMessage(condition), "Optimization was not run",
               fixed = TRUE)
  expect_false(grepl("P01|P11", conditionMessage(condition)))
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

  expect_false(isTRUE(fit$summary$InferenceReady))
  expect_identical(fit$summary$FitReadiness, "review")
  expect_identical(fit$summary$CategoryState, "weak_information")
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
