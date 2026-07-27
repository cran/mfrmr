make_wright_ready_fit <- function(...) {
  fit <- make_toy_fit(...)
  fit$summary$Converged <- TRUE
  fit$summary$InferenceReady <- TRUE
  fit$summary$ConvergenceSeverity <- "pass"
  fit$summary$ConvergenceStatus <- "converged"
  fit
}

test_that("plot readiness preserves a failed numerical state", {
  fit <- make_wright_ready_fit()
  fit$summary$Converged <- FALSE
  fit$summary$InferenceReady <- FALSE
  fit$summary$ConvergenceSeverity <- "fail"
  fit$summary$ConvergenceStatus <- "iteration_limit"

  expect_warning(
    out <- plot(fit, type = "wright", draw = FALSE),
    "Numerical=fail",
    fixed = TRUE
  )
  numerical <- out$data$fit_readiness[
    out$data$fit_readiness$Domain == "Numerical",
    ,
    drop = FALSE
  ]
  expect_identical(numerical$Status, "fail")
  expect_identical(out$data$interpretation_status, "review_only")
  expect_match(out$data$subtitle, "REVIEW ONLY", fixed = TRUE)
})

test_that("native top_n compacts facets but retains every step transition", {
  fit <- make_wright_ready_fit()
  full <- plot(fit, type = "wright", top_n = 100L, draw = FALSE)
  compact <- plot(fit, type = "wright", top_n = 2L, draw = FALSE)
  complete <- plot(fit, type = "wright", top_n = Inf, draw = FALSE)
  facets <- plot(fit, type = "wright", renderer = "facets", top_n = 2L, draw = FALSE)

  expect_identical(compact$data$wright_style, "native")
  expect_true(compact$data$show_ci)
  expect_identical(compact$data$uncertainty_display, "native_mfrmr_ci")
  expect_identical(
    full$data$display_settings$AutoRangePolicy,
    "all_fitted_locations"
  )
  expect_equal(full$data$locations$DisplayEstimate, full$data$locations$Estimate)
  expect_equal(
    full$data$y_range,
    range(c(full$data$locations$Estimate, full$data$person$Estimate))
  )
  full_steps <- full$data$locations[
    full$data$locations$PlotType == "Step threshold",
    ,
    drop = FALSE
  ]
  compact_steps <- compact$data$locations[
    compact$data$locations$PlotType == "Step threshold",
    ,
    drop = FALSE
  ]
  compact_step_labels <- compact$data$label_points[
    compact$data$label_points$PlotType == "Step threshold",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(compact_steps), nrow(full_steps))
  expect_setequal(as.character(compact_steps$Label), as.character(full_steps$Label))
  expect_setequal(
    as.character(compact_step_labels$Label),
    as.character(full_steps$Label)
  )
  expect_true(all(grepl(" -> ", compact_step_labels$Label, fixed = TRUE)))
  expect_identical(compact$data$label_limit, 2L)
  expect_true(all(c("Component", "Shown", "Total", "Omitted", "Complete", "RequestedTopN") %in%
                    names(compact$data$retention)))
  compact_all <- compact$data$retention[
    compact$data$retention$Component == "All locations",
    ,
    drop = FALSE
  ]
  expect_equal(compact_all$Shown + compact_all$Omitted, compact_all$Total)
  expect_gt(compact_all$Omitted, 0L)
  expect_false(compact_all$Complete)
  expect_match(compact$data$retention_note, "Use `top_n = Inf`", fixed = TRUE)
  expect_match(compact$data$subtitle, "facet locations", fixed = TRUE)
  complete_all <- complete$data$retention[
    complete$data$retention$Component == "All locations",
    ,
    drop = FALSE
  ]
  expect_equal(complete_all$Shown, complete_all$Total)
  expect_equal(complete_all$Omitted, 0L)
  expect_true(complete_all$Complete)
  expect_identical(complete$data$retention_note, "")
  expect_equal(nrow(complete$data$locations), complete_all$Total)
  expect_identical(facets$data$renderer, "facets")
  expect_false(facets$data$show_ci)
  expect_identical(facets$data$uncertainty_display, "facets_style_no_ci")
  expect_equal(nrow(facets$data$locations), nrow(full$data$locations))
})

test_that("Wright CI payload stores absolute endpoints", {
  fit <- make_wright_ready_fit()
  out <- suppressWarnings(plot(
    fit,
    type = "wright",
    show_ci = TRUE,
    ci_level = 0.90,
    draw = FALSE
  ))
  loc <- out$data$locations
  loc <- loc[loc$PlotType == "Facet level" & is.finite(loc$SE), , drop = FALSE]
  expect_gt(nrow(loc), 0L)
  z <- stats::qnorm(0.95)
  expect_equal(loc$CI_Lower, loc$Estimate - z * loc$SE)
  expect_equal(loc$CI_Upper, loc$Estimate + z * loc$SE)
  expect_true(all(loc$CI_Lower <= loc$Estimate & loc$Estimate <= loc$CI_Upper))
})

test_that("supplied diagnostics drive Wright SE metadata but not coordinates", {
  fit <- make_wright_ready_fit()
  diagnostics <- make_toy_diagnostics(fit)
  row <- which(as.character(diagnostics$measures$Facet) != "Person")[1]
  facet <- as.character(diagnostics$measures$Facet[row])
  level <- as.character(diagnostics$measures$Level[row])
  fitted_estimate <- fit$facets$others$Estimate[
    as.character(fit$facets$others$Facet) == facet &
      as.character(fit$facets$others$Level) == level
  ]
  fitted_estimate <- unname(fitted_estimate)
  diagnostics$measures$Estimate[row] <- 1.2345
  diagnostics$measures$SE[row] <- 0.123

  out <- plot(
    fit,
    type = "wright",
    diagnostics = diagnostics,
    show_ci = TRUE,
    ci_level = 0.95,
    draw = FALSE
  )
  loc <- out$data$locations[
    as.character(out$data$locations$Group) == facet &
      as.character(out$data$locations$Label) == level,
    , drop = FALSE
  ]
  expect_equal(nrow(loc), 1L)
  expect_equal(loc$Estimate, fitted_estimate)
  expect_equal(loc$SE, 0.123)
  expect_identical(loc$Measure_Source, "diagnostics$measures")
  expect_equal(loc$CI_Lower, fitted_estimate - stats::qnorm(0.975) * 0.123)

  unified <- plot_wright_unified(
    fit,
    diagnostics = diagnostics,
    show_ci = TRUE,
    draw = FALSE
  )
  unified_loc <- unified$locations[
    as.character(unified$locations$Group) == facet &
      as.character(unified$locations$Label) == level,
    , drop = FALSE
  ]
  expect_equal(unified_loc$Estimate, fitted_estimate)
  expect_equal(unified_loc$SE, 0.123)
})

test_that("FACETS-style Wright payload is complete and explicitly visual", {
  fit <- make_wright_ready_fit()
  diagnostics <- make_toy_diagnostics(fit)
  score_values <- as.character(fit$prep$score_map$OriginalScore)
  rubric <- stats::setNames(paste("Rubric", score_values), score_values)
  out <- plot(
    fit,
    type = "wright",
    diagnostics = diagnostics,
    show_ci = TRUE,
    renderer = "facets",
    category_labels = rubric,
    rows_per_logit = 4L,
    persons_per_star = 2,
    draw = FALSE
  )

  expect_identical(out$data$wright_style, "facets_style")
  expect_identical(out$data$renderer, "facets")
  expect_true(out$data$show_ci)
  expect_identical(out$data$uncertainty_display, "hybrid_mfrmr_ci")
  expect_match(out$data$subtitle, "Hybrid", fixed = TRUE)
  expect_true(any(grepl("CI", out$data$legend$label, fixed = TRUE)))
  expect_true(any(out$data$legend$role == "uncertainty"))
  expect_match(out$data$visual_contract, "not FACETS numerical equivalence", fixed = TRUE)
  fs <- out$data$facets_style
  expect_true(all(c(
    "ruler_rows", "person_frequency", "person_placements", "facet_ruler",
    "facet_cells", "step_ruler", "score_transitions", "category_labels",
    "headers", "settings"
  ) %in% names(fs)))
  expect_equal(sum(fs$person_frequency$Count), nrow(out$data$person))
  expect_equal(nrow(fs$facet_ruler), sum(out$data$locations$PlotType == "Facet level"))
  expect_setequal(
    paste(fs$facet_ruler$Facet, fs$facet_ruler$Level, sep = "::"),
    paste(
      out$data$locations$Group[out$data$locations$PlotType == "Facet level"],
      out$data$locations$Label[out$data$locations$PlotType == "Facet level"],
      sep = "::"
    )
  )
  expect_true(all(fs$headers$Header[fs$headers$ColumnType == "facet"] %in%
    paste0(ifelse(fs$facet_ruler$Sign[match(
      fs$headers$Group[fs$headers$ColumnType == "facet"],
      fs$facet_ruler$Facet
    )] >= 0, "+", "-"), fs$headers$Group[fs$headers$ColumnType == "facet"])))
  expect_gt(nrow(fs$step_ruler), 0L)
  expect_true(all(grepl(" -> ", fs$step_ruler$TransitionLabel, fixed = TRUE)))
  expect_true(all(grepl("Rubric", fs$step_ruler$TransitionLabel, fixed = TRUE)))
  expect_equal(fs$step_ruler$DrawValue, fs$step_ruler$Estimate)
  expect_true(all(grepl("[", fs$step_ruler$DrawLabel, fixed = TRUE)))
  expect_gt(nrow(fs$score_transitions), 0L)
  expect_true(all(is.finite(fs$score_transitions$Estimate)))
  expect_equal(
    fs$score_transitions$DrawValue,
    fs$score_transitions$Estimate
  )
  expect_true(all(fs$score_transitions$CrossingStatus == "in_range"))
  curve_spec <- mfrmr:::build_step_curve_spec(fit)
  for (i in seq_len(nrow(fs$score_transitions))) {
    transition <- fs$score_transitions[i, , drop = FALSE]
    curve <- mfrmr:::build_curve_tables(curve_spec, transition$Estimate)$expected
    expected <- curve$ExpectedScore[curve$CurveGroup == transition$CurveGroup]
    expect_equal(expected, transition$MeanHalfScore, tolerance = 1e-5)
  }
  expect_identical(fs$settings$RowsPerLogit, 4L)
  expect_equal(fs$settings$PersonsPerStar, 2)
  expect_equal(fs$settings$StarsPerPerson, 0.5)
  expect_match(fs$settings$VisualCorrespondence, "not a claim", fixed = TRUE)
})

test_that("FACETS-style headers and transitions respect fit metadata", {
  fit <- make_wright_ready_fit()
  positive_facet <- as.character(fit$config$facet_names[1])
  fit$config$facet_signs[[positive_facet]] <- 1
  mapped_scores <- seq_len(nrow(fit$prep$score_map)) * 10
  fit$prep$score_map$OriginalScore <- mapped_scores
  rubric <- stats::setNames(paste("Level", mapped_scores), as.character(mapped_scores))

  out <- plot(
    fit,
    type = "wright",
    renderer = "facets",
    category_labels = rubric,
    draw = FALSE
  )
  fs <- out$data$facets_style
  expect_true(paste0("+", positive_facet) %in% fs$headers$Header)
  expect_equal(fs$category_labels$OriginalScore, as.character(mapped_scores))
  expect_true(all(fs$step_ruler$LowerScore %in% as.character(mapped_scores)))
  expect_true(all(fs$step_ruler$UpperScore %in% as.character(mapped_scores)))
  expect_true(all(grepl("Level", fs$step_ruler$TransitionLabel, fixed = TRUE)))
})

test_that("named Wright category labels must cover original score keys after recoding", {
  fit <- make_wright_ready_fit()
  n_cat <- nrow(fit$prep$score_map)
  internal_scores <- 0:(n_cat - 1L)
  original_scores <- seq_len(n_cat)
  fit$prep$rating_min <- 0L
  fit$prep$rating_max <- n_cat - 1L
  fit$config$n_cat <- n_cat
  fit$prep$score_map <- data.frame(
    OriginalScore = original_scores,
    InternalScore = internal_scores
  )

  valid <- stats::setNames(paste("Level", original_scores), original_scores)
  out <- plot(
    fit,
    type = "wright",
    renderer = "facets",
    category_labels = valid,
    draw = FALSE
  )
  expect_equal(
    out$data$facets_style$category_labels$OriginalScore,
    as.character(original_scores)
  )
  expect_equal(
    out$data$facets_style$category_labels$RubricLabel,
    unname(valid)
  )

  internal_keyed <- stats::setNames(paste("Internal", internal_scores), internal_scores)
  expect_error(
    plot(
      fit,
      type = "wright",
      renderer = "facets",
      category_labels = internal_keyed,
      draw = FALSE
    ),
    "appear to use internal score keys"
  )
  expect_error(
    plot(
      fit,
      type = "wright",
      renderer = "facets",
      category_labels = valid[-n_cat],
      draw = FALSE
    ),
    paste0("original score key\\(s\\): ", n_cat)
  )
  expect_error(
    plot(
      fit,
      type = "wright",
      renderer = "facets",
      category_labels = data.frame(
        Score = original_scores[-n_cat],
        Label = unname(valid[-n_cat])
      ),
      draw = FALSE
    ),
    paste0("original score key\\(s\\): ", n_cat)
  )
})

test_that("FACETS-style range, ruler resolution, and extremes are controllable", {
  fit <- make_wright_ready_fit()
  fit$facets$person$Extreme[1:2] <- c("low", "high")
  out <- plot(
    fit,
    type = "wright",
    renderer = "facets",
    wright_range = c(-2, 2),
    rows_per_logit = 2L,
    extreme_placement = "ends",
    draw = FALSE
  )
  fs <- out$data$facets_style
  expect_equal(range(fs$ruler_rows$Logit), c(-2, 2))
  expect_equal(unique(diff(fs$ruler_rows$Logit)), 0.5)
  expect_equal(fs$person_placements$RulerValue[1:2], c(-2, 2))

  estimates <- plot(
    fit,
    type = "wright",
    renderer = "facets",
    wright_range = c(-2, 2),
    extreme_placement = "estimate",
    draw = FALSE
  )
  expect_false(identical(
    estimates$data$facets_style$person_placements$DisplayEstimate[1:2],
    c(-2, 2)
  ))
})

test_that("FACETS-style thresholds can be omitted and base renderer draws", {
  fit <- make_wright_ready_fit()
  no_steps <- plot_wright_unified(
    fit,
    renderer = "facets",
    show_thresholds = FALSE,
    draw = FALSE
  )
  expect_equal(nrow(no_steps$facets_style$step_ruler), 0L)
  expect_equal(nrow(no_steps$facets_style$score_transitions), 0L)

  grDevices::pdf(nullfile())
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(
    fit,
    type = "wright",
    renderer = "native",
    top_n = 2L
  ))
  expect_silent(plot(
    fit,
    type = "wright",
    renderer = "facets",
    show_ci = FALSE
  ))
})

test_that("plot_wright_unified uses the same renderer-specific uncertainty contract", {
  fit <- make_wright_ready_fit()

  native <- plot_wright_unified(fit, top_n = Inf, draw = FALSE)
  expect_true(native$show_ci)
  expect_identical(native$uncertainty_display, "native_mfrmr_ci")
  expect_true(all(native$retention$Omitted == 0L))

  facets <- plot_wright_unified(fit, renderer = "facets", draw = FALSE)
  expect_false(facets$show_ci)
  expect_identical(facets$uncertainty_display, "facets_style_no_ci")

  hybrid <- plot_wright_unified(
    fit,
    renderer = "facets",
    show_ci = TRUE,
    draw = FALSE
  )
  expect_true(hybrid$show_ci)
  expect_identical(hybrid$uncertainty_display, "hybrid_mfrmr_ci")
  expect_match(hybrid$title, "hybrid", ignore.case = TRUE)
})

test_that("boundary separation uses one robust range and explicit CI display metadata", {
  fit <- make_wright_ready_fit()
  boundary_rows <- which(
    as.character(fit$facets$others$Facet) == "Rater"
  )[1:2]
  boundary_levels <- as.character(fit$facets$others$Level[boundary_rows])
  fit$facets$others$Estimate[boundary_rows] <- c(-25, 25)
  fit$data_review$boundary_levels <- data.frame(
    Facet = "Rater",
    Level = boundary_levels,
    stringsAsFactors = FALSE
  )

  native <- suppressWarnings(plot(
    fit,
    type = "wright",
    top_n = Inf,
    show_ci = TRUE,
    draw = FALSE
  ))
  settings <- native$data$display_settings
  expect_identical(settings$AutoRangePolicy, "boundary_levels_at_ends")
  expect_equal(native$data$y_range, c(settings$LowerLogit, settings$UpperLogit))
  expect_lt(diff(native$data$y_range), 10)

  boundary_native <- native$data$locations[
    native$data$locations$BoundarySeparated,
    ,
    drop = FALSE
  ]
  expect_equal(boundary_native$OriginalEstimate, boundary_native$Estimate)
  expect_equal(boundary_native$CI_Lower, boundary_native$OriginalCI_Lower)
  expect_equal(boundary_native$CI_Upper, boundary_native$OriginalCI_Upper)
  expect_true(all(boundary_native$DisplayEstimate >= native$data$y_range[1]))
  expect_true(all(boundary_native$DisplayEstimate <= native$data$y_range[2]))
  expect_setequal(boundary_native$BoundaryEnd, c("lower", "upper"))
  expect_true(all(boundary_native$CIDisplayStatus == "boundary_endpoint_marker"))
  expect_true(all(grepl("(", boundary_native$DisplayLabel, fixed = TRUE)))

  drawn_key <- mfrmr:::.wright_native_legend_spec(
    native$data,
    show_ci = TRUE,
    ci_level = 0.95
  )
  expect_identical(as.character(native$data$legend$label), drawn_key$label)
  expect_true(all(c(
    "Person density", "Facet level", "Step threshold", "IQR / range",
    "95% mfrmr CI", "Boundary level / CI beyond ruler"
  ) %in% native$data$legend$label))

  facets <- suppressWarnings(plot(
    fit,
    type = "wright",
    renderer = "facets",
    show_ci = TRUE,
    draw = FALSE
  ))
  fs <- facets$data$facets_style
  expect_equal(
    c(fs$settings$LowerLogit, fs$settings$UpperLogit),
    native$data$y_range
  )
  expect_identical(fs$settings$AutoRangePolicy, "boundary_levels_at_ends")
  boundary_facets <- fs$facet_ruler[fs$facet_ruler$BoundarySeparated, , drop = FALSE]
  expect_equal(boundary_facets$OriginalEstimate, boundary_facets$Estimate)
  expect_equal(boundary_facets$OriginalCI_Lower, boundary_facets$CI_Lower)
  expect_true(all(boundary_facets$CISuppressed))

  grDevices::pdf(nullfile(), width = 12, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(mfrmr:::draw_wright_map(
    native$data,
    show_ci = TRUE,
    ci_level = 0.95
  ))
  expect_silent(mfrmr:::draw_wright_map(
    facets$data,
    show_ci = TRUE,
    ci_level = 0.95
  ))
})

test_that("FACETS-style Wright arguments are validated", {
  fit <- make_wright_ready_fit()
  expect_error(
    plot(fit, type = "wright", top_n = 0, draw = FALSE),
    "top_n"
  )
  expect_error(
    plot(fit, type = "wright", wright_style = "facets", draw = FALSE),
    "not recognized"
  )
  expect_error(
    plot(fit, type = "wright", wright_style = "facets_style", rows_per_logit = 0, draw = FALSE),
    "rows_per_logit"
  )
  expect_error(
    plot(fit, type = "wright", wright_style = "facets_style", wright_range = c(2, -2), draw = FALSE),
    "wright_range"
  )
  expect_error(
    plot(fit, type = "wright", wright_style = "facets_style", persons_per_star = 0, draw = FALSE),
    "persons_per_star"
  )
  expect_warning(
    plot(
      fit,
      type = "wright",
      renderer = "facets",
      group = rep(c("A", "B"), length.out = nrow(fit$facets$person)),
      draw = FALSE
    ),
    "available only"
  )
  expect_error(
    plot(fit, type = "wright", renderer = "facets_style", draw = FALSE),
    "renderer"
  )
})
