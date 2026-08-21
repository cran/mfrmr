make_person_boundary_fixture <- function(method = "JML", anchors = NULL,
                                         groups = NULL,
                                         group_values = NULL,
                                         centered = FALSE) {
  persons <- c("P1", "P2", "P3")
  config <- list(
    method = method,
    model = "RSM",
    theta_spec = mfrmr:::build_facet_constraint(
      levels = persons,
      anchors = anchors,
      groups = groups,
      group_values = group_values,
      centered = centered
    )
  )
  prep <- list(
    data = data.frame(
      Person = rep(persons, each = 2L),
      Score = c(1L, 1L, 1L, 4L, 4L, 4L)
    ),
    rating_min = 1L,
    rating_max = 4L
  )
  person_table <- data.frame(
    Person = persons,
    Estimate = c(-12, 0, 13),
    SE = NA_real_,
    Extreme = c("low", "none", "high"),
    stringsAsFactors = FALSE
  )
  list(config = config, prep = prep, person_table = person_table)
}

test_that("free JML extreme Persons have typed unbounded primary estimates", {
  fixture <- make_person_boundary_fixture()
  fixture$person_table$Extreme <- "none"
  audit <- mfrmr:::audit_mfrm_person_boundary(
    fixture$person_table, fixture$config, fixture$prep
  )
  result <- mfrmr:::apply_mfrm_person_boundary(
    fixture$person_table, audit
  )

  expect_identical(audit$readiness$BoundaryState, "has_exclusions")
  expect_true(isTRUE(audit$readiness$Complete))
  expect_identical(
    audit$parameter_status$ParameterStatus,
    c("unbounded_low", "estimable", "unbounded_high")
  )
  expect_identical(result$Estimate, c(-Inf, 0, Inf))
  expect_equal(result$OptimizerEstimate, c(-12, 0, 13))
  expect_true(all(is.na(result$DisplayEstimate[c(1L, 3L)])))
  expect_identical(
    result$OptimizerEstimateUse[c(1L, 3L)],
    c("numerical_trace_only", "numerical_trace_only")
  )
  expect_identical(
    result$ReasonCodes[c(1L, 3L)],
    c("jml_extreme_low", "jml_extreme_high")
  )
  expect_identical(result$Extreme, c("low", "none", "high"))
  expect_identical(result$ResponseRows, rep(2L, 3L))
  expect_equal(result$WeightedResponseTotal, rep(2, 3L))
})

test_that("fixed extreme JML Persons are not relabelled as unbounded", {
  fixture <- make_person_boundary_fixture(anchors = c(P3 = 0.75))
  fixture$prep$data$Score[fixture$prep$data$Person == "P1"] <- c(1L, 4L)
  fixture$person_table$Extreme[fixture$person_table$Person == "P1"] <- "none"
  fixture$person_table$Estimate[fixture$person_table$Person == "P3"] <- 0.75
  audit <- mfrmr:::audit_mfrm_person_boundary(
    fixture$person_table, fixture$config, fixture$prep
  )
  p3 <- audit$parameter_status[
    audit$parameter_status$Level == "P3", , drop = FALSE
  ]

  expect_identical(p3$ParameterStatus, "fixed")
  expect_equal(p3$PrimaryEstimate, 0.75)
  expect_true(is.na(p3$OptimizerEstimate))
  expect_identical(p3$ReasonCodes, "fixed_extreme_response")
  expect_identical(audit$readiness$BoundaryState, "finite")
  expect_identical(audit$readiness$UnboundedLowN, 0L)
  expect_identical(audit$readiness$UnboundedHighN, 0L)
})

test_that("MML extreme responses retain finite prior-regularized EAP status", {
  fixture <- make_person_boundary_fixture(method = "MML")
  fixture$person_table$Estimate <- c(-1.4, 0.1, 1.7)
  audit <- mfrmr:::audit_mfrm_person_boundary(
    fixture$person_table, fixture$config, fixture$prep
  )
  result <- mfrmr:::apply_mfrm_person_boundary(
    fixture$person_table, audit
  )

  expect_identical(audit$readiness$BoundaryState, "finite")
  expect_identical(audit$readiness$ReasonCodes, "")
  expect_true(all(result$ParameterStatus == "estimable"))
  expect_true(all(is.finite(result$Estimate)))
  expect_true(all(is.na(result$OptimizerEstimate)))
  expect_identical(
    result$ReasonCodes[c(1L, 3L)],
    rep("mml_extreme_response_prior_regularized", 2L)
  )
})

test_that("constraint-coupled JML extremes fail closed as weak information", {
  group <- setNames(rep("G1", 3L), c("P1", "P2", "P3"))
  fixture <- make_person_boundary_fixture(
    groups = group,
    group_values = c(G1 = 0)
  )
  audit <- mfrmr:::audit_mfrm_person_boundary(
    fixture$person_table, fixture$config, fixture$prep
  )
  status <- audit$parameter_status

  expect_identical(audit$readiness$BoundaryState, "not_evaluated")
  expect_false(isTRUE(audit$readiness$Complete))
  expect_identical(audit$readiness$ConstraintCoupledExtremeN, 2L)
  expect_true(all(
    status$ParameterStatus[c(1L, 3L)] == "weak_information"
  ))
  expect_true(all(
    status$ReasonCodes[c(1L, 3L)] == "weak_design_information"
  ))
  expect_true(all(is.finite(status$PrimaryEstimate)))
})

test_that("fitted JML objects separate boundary estimates from optimizer traces", {
  data <- simulate_mfrm_data(
    n_person = 24,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 4,
    assignment = "crossed",
    seed = 803
  )
  data$Score[data$Person == "P001"] <- 4L
  data$Score[data$Person == "P002"] <- 1L
  # Missing rows are removed before the sufficient-score audit. The retained
  # positive-weight contributing observations for P001 remain all-maximum.
  p001_rows <- which(data$Person == "P001")
  data$Weight <- 1
  data$Score[p001_rows[1L]] <- NA_integer_
  data$Score[p001_rows[2L]] <- 1L
  data$Weight[p001_rows[2L]] <- 0

  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    weight = "Weight",
    rating_min = 1,
    rating_max = 4,
    method = "JML",
    model = "RSM",
    maxit = 120
  ))
  persons <- fit$facets$person[
    fit$facets$person$Person %in% c("P001", "P002"), , drop = FALSE
  ]
  persons <- persons[match(c("P001", "P002"), persons$Person), , drop = FALSE]

  expect_identical(fit$summary$BoundaryState, "has_exclusions")
  expect_identical(fit$summary$InputState, "review")
  expect_identical(fit$summary$FitReadiness, "review")
  expect_false(isTRUE(fit$summary$InferenceReady))
  expect_identical(
    fit$readiness$fit$FitReadiness,
    as.character(fit$summary$FitReadiness)
  )
  expect_identical(persons$Estimate, c(Inf, -Inf))
  expect_true(all(is.finite(persons$OptimizerEstimate)))
  expect_identical(
    persons$ParameterStatus,
    c("unbounded_high", "unbounded_low")
  )
  expect_identical(
    persons$ResponseRows,
    as.integer(vapply(c("P001", "P002"), function(person) {
      sum(as.character(fit$prep$data$Person) == person)
    }, integer(1)))
  )
  expect_equal(
    persons$WeightedResponseTotal,
    unname(vapply(c("P001", "P002"), function(person) {
      sum(fit$prep$data$Weight[as.character(fit$prep$data$Person) == person])
    }, numeric(1)))
  )
  expect_identical(
    fit$data_review$boundary$readiness$AuditedParameterClass,
    "Person"
  )

  diagnostics <- suppressWarnings(diagnose_mfrm(fit))
  diagnostic_persons <- diagnostics$measures[
    diagnostics$measures$Facet == "Person" &
      diagnostics$measures$Level %in% c("P001", "P002"), , drop = FALSE
  ]
  expect_true(all(diagnostic_persons$BoundaryExcluded))
  expect_true(all(is.na(diagnostic_persons$SE)))
  expect_true(all(is.na(diagnostic_persons$CI_Lower)))
  expect_true(all(is.na(diagnostic_persons$CI_Upper)))
  expect_false(any(diagnostic_persons$SupportsFormalInference))

  ends <- suppressWarnings(plot(
    fit,
    type = "wright",
    renderer = "facets",
    show_ci = FALSE,
    draw = FALSE,
    extreme_placement = "ends"
  ))
  end_persons <- ends$data$facets_style$person[
    ends$data$facets_style$person$Person %in% c("P001", "P002"),
    , drop = FALSE
  ]
  expect_true(all(is.finite(end_persons$DisplayEstimate)))
  expect_true(all(is.infinite(end_persons$OriginalEstimate)))

  estimate <- suppressWarnings(plot(
    fit,
    type = "wright",
    renderer = "facets",
    show_ci = FALSE,
    draw = FALSE,
    extreme_placement = "estimate"
  ))
  estimate_persons <- estimate$data$facets_style$person[
    estimate$data$facets_style$person$Person %in% c("P001", "P002"),
    , drop = FALSE
  ]
  expect_true(all(is.na(estimate_persons$DisplayEstimate)))
  expect_true(all(is.infinite(estimate_persons$OriginalEstimate)))

  all_boundary <- fit
  all_boundary$facets$person$Estimate <- rep(
    c(-Inf, Inf),
    length.out = nrow(all_boundary$facets$person)
  )
  all_boundary$facets$person$ParameterStatus <- rep(
    c("unbounded_low", "unbounded_high"),
    length.out = nrow(all_boundary$facets$person)
  )
  all_boundary$facets$person$Extreme <- rep(
    c("low", "high"),
    length.out = nrow(all_boundary$facets$person)
  )
  all_boundary_plot <- suppressWarnings(plot(
    all_boundary,
    type = "wright",
    renderer = "facets",
    show_ci = FALSE,
    draw = FALSE,
    extreme_placement = "ends"
  ))
  expect_identical(all_boundary_plot$data$person_stats$FiniteN, 0L)
  expect_identical(
    all_boundary_plot$data$person_stats$BoundaryExcludedN,
    nrow(all_boundary$facets$person)
  )
  expect_true(is.na(all_boundary_plot$data$person_stats$Mean))
})
