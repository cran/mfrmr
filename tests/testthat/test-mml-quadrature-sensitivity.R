ordered_quadrature_fixture <- local({
  cache <- new.env(parent = emptyenv())
  function(model = c("RSM", "PCM")) {
    model <- match.arg(model)
    if (exists(model, envir = cache, inherits = FALSE)) {
      return(get(model, envir = cache, inherits = FALSE))
    }
    data <- load_mfrmr_data("example_core")
    persons <- unique(as.character(data$Person))[seq_len(12L)]
    data <- data[as.character(data$Person) %in% persons, , drop = FALSE]
    fit <- suppressWarnings(fit_mfrm(
      data,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      method = "MML",
      model = model,
      step_facet = if (identical(model, "PCM")) "Criterion" else NULL,
      quad_points = 5L,
      maxit = 80L,
      anchor_policy = "silent"
    ))
    sensitivity <- suppressWarnings(mml_quadrature_sensitivity(
      fit, data, quad_points = c(5L, 7L), theta_points = 41L
    ))
    out <- list(data = data, fit = fit, sensitivity = sensitivity)
    assign(model, out, envir = cache)
    out
  }
})

test_that("RSM and PCM expose continuous same-data quadrature movement", {
  for (model in c("RSM", "PCM")) {
    fixture <- ordered_quadrature_fixture(model)
    sensitivity <- fixture$sensitivity
    comparison <- sensitivity$summary

    expect_s3_class(sensitivity, "mfrm_quadrature_sensitivity")
    expect_identical(sensitivity$settings$model, model)
    expect_identical(sensitivity$settings$quad_points, c(5L, 7L))
    expect_identical(names(sensitivity$fits), c("q5", "q7"))
    expect_identical(nrow(sensitivity$slopes), 0L)
    expect_true(all(is.na(comparison$SlopeMaxAbsChange)))
    expect_true(all(is.finite(comparison$MeasurementParameterMaxAbsChange)))
    expect_true(all(is.finite(comparison$ProbabilityMaxAbsChange)))
    expect_true(all(is.finite(comparison$EAPMaxAbsChange)))
    expect_true(all(is.finite(comparison$PosteriorSDMaxAbsChange)))
    expect_equal(
      as.numeric(comparison[comparison$IsReference, c(
        "NLLAbsChangePerPerson", "MeasurementParameterMaxAbsChange",
        "ProbabilityMaxAbsChange", "EAPMaxAbsChange",
        "PosteriorSDMaxAbsChange"
      )]),
      rep(0, 5L),
      tolerance = 1e-12
    )
    expect_identical(
      summary(sensitivity)$overview$StabilityClassification,
      "not_assigned_continuous_evidence_only"
    )
  }
})

test_that("RSM quadrature probabilities retain fitted interactions", {
  data <- simulate_mfrm_data(
    n_person = 24L,
    n_rater = 3L,
    n_criterion = 3L,
    raters_per_person = 3L,
    score_levels = 4L,
    model = "RSM",
    seed = 90603L
  )
  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    facet_interactions = "Rater:Criterion",
    min_obs_per_interaction = 0,
    quad_points = 5L,
    maxit = 80L,
    mml_engine = "direct"
  ))
  sensitivity <- suppressWarnings(mml_quadrature_sensitivity(
    fit, data, quad_points = c(5L, 7L), theta_points = 41L
  ))

  expect_true(nrow(fit$interactions$effects) > 0L)
  expect_true(all(is.finite(
    sensitivity$summary$MeasurementParameterMaxAbsChange
  )))
  expect_true(all(is.finite(sensitivity$summary$ProbabilityMaxAbsChange)))
  expect_identical(
    sensitivity$settings$probability_contract,
    "all_observed_facet_cells_and_fitted_interactions_on_common_theta_grid"
  )
})

test_that("the general entry point keeps the GPCM-specific scope guard", {
  fixture <- ordered_quadrature_fixture("RSM")
  expect_error(
    gpcm_mml_quadrature_sensitivity(
      fixture$fit, fixture$data, quad_points = c(5L, 7L),
      theta_points = 41L
    ),
    "requires a GPCM MML fit",
    fixed = TRUE
  )
})
