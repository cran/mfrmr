test_that("adaptive review resolves narrow and tail posteriors without changing the fixed grid", {
  for (case in c("narrow", "tail")) {
    scores <- if (case == "narrow") rep(0:3, 250L) else rep(3L, 100L)
    difficulty <- if (case == "narrow") 0 else 16
    reference <- if (case == "narrow") {
      c(log_marginal = -1389.8600406660846, eap = 0, sd = 0.028282562846515886)
    } else {
      c(log_marginal = -177.87927093394117, eap = 17.903577120328180, sd = 0.21878638089256411)
    }
    q <- mfrmr:::gauss_hermite_normal(301L)
    before <- q
    out <- mfrmr:::mfrmr_adaptive_quadrature_review(
      list(person = rep(1L, length(scores)), score_k = scores),
      list(model = "RSM", n_cat = 4L), list(steps = c(0, 0, 0)),
      q, "P", c(31L, 61L), base_eta = rep(-difficulty, length(scores))
    )
    expect_identical(q, before)
    expect_true(all(out$Status == "computed"))
    expect_lt(abs(out$AdaptiveLogMarginal[2L] - reference["log_marginal"]), 1e-8)
    expect_lt(abs(out$AdaptiveEAP[2L] - reference["eap"]), 1e-8)
    expect_lt(abs(out$AdaptivePosteriorSD[2L] - reference["sd"]), 1e-8)
    expect_lt(abs(out$AdaptiveEAPChangeFromPrevious[2L]), 1e-6)
    expect_lt(abs(out$AdaptiveSDChangeFromPrevious[2L]), 1e-6)
    if (case == "narrow") {
      expect_gt(abs(out$LogMarginalChange[2L]), 0.9)
      expect_gt(out$PosteriorSDChange[2L], 0.028)
    }
  }
})

test_that("adaptive review validates orders and retains failed Persons", {
  for (points in list(31, c(31, 31), c(1, 31), c(3, NA), c(3, 4.5), TRUE)) {
    expect_error(mfrmr:::mfrmr_validate_adaptive_quad_points(points), "adaptive_quad_points")
  }
  expect_error(mfrmr:::mfrmr_validate_adaptive_quad_points(c(31, 400)), "cannot all be represented")
  out <- mfrmr:::mfrmr_adaptive_quadrature_review(
    list(person = c(1L, 2L), score_k = c(1L, 1L), slope_idx = c(1L, 2L), step_idx = c(1L, 2L)),
    list(model = "GPCM", n_cat = 3L),
    list(slopes = c(1, Inf), steps_mat = matrix(0, 2, 2)),
    mfrmr:::gauss_hermite_normal(31L), c("valid", "invalid"), c(15L, 31L),
    base_eta = c(0, 0)
  )
  expect_identical(out$Status, c("computed", "computed", "unavailable", "unavailable"))
  expect_true(all(is.na(out$AdaptiveEAP[out$Person == "invalid"])))
  expect_true(all(nzchar(out$Detail[out$Person == "invalid"])))
  expect_equal(mfrmr:::mfrmr_adaptive_quadrature_overview(out)$Unavailable, c(1L, 1L))
})

test_that("adaptive review honors transformed priors and zero observation weights", {
  out <- mfrmr:::mfrmr_adaptive_quadrature_review(
    list(person = 1L, score_k = 3L, weight = 0),
    list(model = "RSM", n_cat = 4L), list(steps = c(-1, 0, 1)),
    mfrmr:::gauss_hermite_normal(31L), "P", c(15L, 31L), base_eta = -5,
    population_spec = list(active = TRUE, design_matrix = matrix(1, 1, 1),
                           coefficients = 2, sigma2 = 9, person_lookup = 1L)
  )
  expect_true(all(out$Status == "computed"))
  expect_equal(out$AdaptiveLogMarginal, c(0, 0), tolerance = 1e-12)
  expect_equal(out$AdaptiveEAP, c(2, 2), tolerance = 1e-12)
  expect_equal(out$AdaptivePosteriorSD, c(3, 3), tolerance = 1e-12)
})

test_that("public fit and artifact scoring share an optional unrounded review", {
  # Extended local/CI regression; see CONTRIBUTING.md.
  skip_on_cran()
  data <- load_mfrmr_data("example_operational")
  scoring_data <- data
  for (model in c("RSM", "PCM")) {
    fit <- fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score", model = model,
                    step_facet = if (model == "PCM") "Criterion" else NULL)
    # This comparison needs a scoring grid that passes the batch accuracy check.
    # Calibration stays fixed; coarse-grid refusals have separate regressions.
    ordinary <- predict_mfrm_units(fit, scoring_data, scoring_quad_points = 121L)
    reviewed <- predict_mfrm_units(fit, scoring_data, scoring_quad_points = 121L,
      adaptive_quad_points = c(121L, 241L))
    expect_identical(ordinary$estimates, reviewed$estimates)
    expect_identical(ordinary$settings$source_scoring_ready, reviewed$settings$source_scoring_ready)
    expect_identical(summary(reviewed)$quadrature_review, reviewed$quadrature_review)
    expect_equal(reviewed$quadrature_review$FixedEAP[
                   match(ordinary$estimates$Person, reviewed$quadrature_review$Person)],
                 unname(ordinary$estimates$Estimate), tolerance = 1e-12)
    draft <- mfrmr:::mfrmr_extract_calibration_draft(fit, scoring_quad_points = 121L)
    calibration <- freeze_mfrm_calibration(validate_mfrm_calibration(draft))
    before <- calibration
    scored <- score_mfrm_calibration(calibration, scoring_data, adaptive_quad_points = c(121L, 241L))
    expect_identical(calibration, before)
    expect_identical(scored$estimates, score_mfrm_calibration(calibration, scoring_data)$estimates)
    numerical <- vapply(reviewed$quadrature_review, is.numeric, TRUE)
    expect_identical(scored$quadrature_review[!numerical], reviewed$quadrature_review[!numerical])
    expect_identical(is.na(scored$quadrature_review), is.na(reviewed$quadrature_review))
    # Artifact grids are sorted; accumulation order can change fixed-grid
    # results by a few ulps. Check absolute error, including near-zero deltas.
    expect_lt(max(abs(as.matrix(scored$quadrature_review[numerical]) -
                      as.matrix(reviewed$quadrature_review[numerical])), na.rm = TRUE), 1e-12)
    expect_identical(summary(scored)$quadrature_review, scored$quadrature_review)
    if (model == "RSM") {
      selected <- data[rev(which(data$Person %in% unique(data$Person)[c(2L, 7L)])), ]
      selected$Wt <- rep(c(0.25, 1.75, 3), length.out = nrow(selected))
      fitted_weighted <- predict_mfrm_units(
        fit, selected, weight = "Wt", scoring_quad_points = 121L,
        adaptive_quad_points = c(121L, 241L)
      )
      artifact_weighted <- score_mfrm_calibration(
        calibration, selected, weight = "Wt", adaptive_quad_points = c(121L, 241L)
      )
      a <- artifact_weighted$quadrature_review
      b <- fitted_weighted$quadrature_review
      expect_equal(a$AdaptiveEAP,
                   b$AdaptiveEAP[match(paste(a$Person, a$AdaptiveNodes),
                                       paste(b$Person, b$AdaptiveNodes))], tolerance = 1e-12)
      selected$Score[selected$Person == selected$Person[1L]] <- NA
      omitted <- score_mfrm_calibration(
        calibration, selected, weight = "Wt", missing_response = "omit",
        adaptive_quad_points = c(121L, 241L)
      )
      expect_equal(nrow(omitted$quadrature_review), 2L)
      expect_true(all(omitted$quadrature_review$Person != selected$Person[1L]))
      expect_true(all(omitted$quadrature_review$Status == "computed"))
    }
    sensitivity <- mml_quadrature_sensitivity(
      fit, data, quad_points = c(31L, 41L), adaptive_quad_points = c(15L, 31L)
    )
    expect_identical(sensitivity$fits$q31, fit)
    expect_equal(nrow(sensitivity$quadrature_review), 48L * 2L * 2L)
    expect_true(all(sensitivity$quadrature_review$Status == "computed"))
    expect_identical(summary(sensitivity)$quadrature_review, sensitivity$quadrature_review)
    expect_identical(summary(sensitivity)$overview$ReadinessEffect, "none_diagnostic_only")
  }
})

test_that("moving-node gradients differentiate the whole adaptive objective", {
  data <- load_mfrmr_data("example_operational")
  data <- data[data$Person %in% unique(data$Person)[1:6], ]
  data$Weight <- rep(c(0.25, 1, 1.75), length.out = nrow(data))
  for (model in c("RSM", "PCM", "GPCM")) {
    # Only the configuration is needed; differentiation is checked away from a
    # solution, with low orders that expose omitted mode/scale derivatives.
    fit <- suppressWarnings(fit_mfrm(
      data, "Person", c("Rater", "Criterion"), "Score", weight = "Weight",
      model = model, step_facet = if (model != "RSM") "Criterion" else NULL,
      slope_facet = if (model == "GPCM") "Criterion" else NULL,
      facet_interactions = if (model == "RSM") "Rater:Criterion" else NULL,
      min_obs_per_interaction = 0, quad_points = 7L, maxit = 30L,
      anchors = data.frame(Facet = "Rater", Level = as.character(data$Rater[1L]), Anchor = -0.4),
      anchor_policy = "silent",
      population_formula = if (model == "RSM") ~ X else NULL,
      person_data = if (model == "RSM") data.frame(
        Person = unique(data$Person), X = rep(c(-1, 1), 3)) else NULL
    ))
    config <- fit$config
    # Exercise the actual typed step Jacobian, including a fixed transition.
    specs <- mfrmr:::mfrmr_step_specs(config)
    specs[[1L]] <- mfrmr:::build_step_constraint(
      config$n_cat - 1L, anchors = c(`1` = -0.5), scope = specs[[1L]]$scope
    )
    config$step_specs <- specs
    sizes <- mfrmr:::build_param_sizes(config)
    idx <- mfrmr:::build_indices(fit$prep, config$step_facet, config$slope_facet,
                                 config$interaction_specs)
    par <- mfrmr:::build_initial_param_vector(config, sizes)
    par <- par + 0.15 * sin(seq_along(par))
    orders <- if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) c(1L, 3L, 15L) else 3L
    for (order in orders) {
      evaluate <- mfrmr:::mfrmr_make_adaptive_mml_evaluator(idx, config, sizes, order)
      objective <- function(x) {
        review <- mfrmr:::mfrmr_adaptive_quadrature_review(
          idx, config, mfrmr:::expand_params(x, sizes, config),
          list(nodes = 0, weights = 1), fit$prep$levels$Person, order
        )
        stopifnot(all(review$Status == "computed"))
        -sum(review$AdaptiveLogMarginal)
      }
      difference <- function(step) vapply(seq_along(par), function(j) {
        plus <- minus <- par
        plus[j] <- plus[j] + step
        minus[j] <- minus[j] - step
        (objective(plus) - objective(minus)) / (2 * step)
      }, 0)
      result <- evaluate(par)
      g1 <- difference(1e-4)
      g2 <- difference(5e-5)
      expect_lt(abs(result$value - objective(par)), 1e-10)
      expect_lt(max(abs(result$gradient - (4 * g2 - g1) / 3)), 1e-6)
      expect_error(evaluate(c(par, 0)), "expected length")
    }
    idx$weight[] <- 0
    zero <- mfrmr:::mfrmr_make_adaptive_mml_evaluator(idx, config, sizes, 3L)(par)
    expect_lt(abs(zero$value), 1e-12)
    expect_lt(max(abs(zero$gradient)), 1e-12)
    idx$weight[1L] <- -1
    expect_error(mfrmr:::mfrmr_make_adaptive_mml_evaluator(idx, config, sizes, 3L),
                 "non-negative")
  }
})
