gpcm_quadrature_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      data <- load_mfrmr_data("example_core")
      persons <- unique(as.character(data$Person))[seq_len(12L)]
      data <- data[as.character(data$Person) %in% persons, , drop = FALSE]
      fit <- fit_mfrm(
        data,
        person = "Person",
        facets = c("Rater", "Criterion"),
        score = "Score",
        method = "MML",
        model = "GPCM",
        step_facet = "Criterion",
        slope_facet = "Criterion",
        quad_points = 5L,
        maxit = 100L,
        anchor_policy = "silent"
      )
      sensitivity <- gpcm_mml_quadrature_sensitivity(
        fit,
        data,
        quad_points = c(5L, 7L),
        theta_points = 41L
      )
      cached <<- list(data = data, fit = fit, sensitivity = sensitivity)
    }
    cached
  }
})

test_that("GPCM quadrature sensitivity separates numerical evidence", {
  fixture <- gpcm_quadrature_fixture()
  sensitivity <- fixture$sensitivity

  expect_s3_class(sensitivity, "mfrm_quadrature_sensitivity")
  expect_identical(sensitivity$settings$reference_nodes, 5L)
  expect_identical(sensitivity$settings$quad_points, c(5L, 7L))
  expect_identical(nrow(sensitivity$summary), 2L)
  expect_identical(nrow(sensitivity$runs), 2L)
  expect_identical(nrow(sensitivity$slopes), 8L)
  expect_identical(names(sensitivity$fits), c("q5", "q7"))
  expect_identical(sensitivity$fits$q5, fixture$fit)

  reference <- sensitivity$summary$IsReference
  expect_identical(reference, c(TRUE, FALSE))
  change_columns <- c(
    "NLLAbsChangePerPerson", "SlopeMaxAbsChange",
    "RawSlopeSEMaxAbsChange", "PopulationSDAbsChange",
    "RawPopulationSDSEAbsChange", "ProbabilityMaxAbsChange"
  )
  expect_equal(
    as.numeric(sensitivity$summary[reference, change_columns]),
    rep(0, length(change_columns)),
    tolerance = 1e-12
  )
  expect_true(all(is.finite(
    as.numeric(sensitivity$summary[!reference, change_columns])
  )))
  expect_true(all(
    as.numeric(sensitivity$summary[!reference, change_columns]) >= 0
  ))

  expect_true(all(sensitivity$runs$EstimationConverged))
  expect_true(all(
    sensitivity$runs$HessianRank == sensitivity$runs$HessianDimension
  ))
  expect_true(all(sensitivity$runs$MinimumEigenvalue > 0))
  expect_true(all(is.finite(sensitivity$runs$HessianConditionNumber)))
  expect_true(all(!sensitivity$runs$InferenceReady))
  expect_true(all(grepl(
    "mml_gpcm_slope_boundary_not_evaluated",
    sensitivity$runs$ReadinessReasons,
    fixed = TRUE
  )))
  expect_false(sensitivity$settings$stability_threshold_frozen)
  expect_identical(
    sensitivity$settings$readiness_effect,
    "none_diagnostic_only"
  )
})

test_that("quadrature summaries and APA tables remain non-decisional", {
  sensitivity <- gpcm_quadrature_fixture()$sensitivity
  summary_object <- summary(sensitivity)

  expect_s3_class(summary_object, "summary.mfrm_quadrature_sensitivity")
  expect_true(summary_object$overview$AllEstimationConverged)
  expect_true(summary_object$overview$AllHessiansFullRank)
  expect_identical(
    summary_object$overview$StabilityClassification,
    "not_assigned_continuous_evidence_only"
  )
  expect_identical(
    summary_object$overview$ReadinessEffect,
    "none_diagnostic_only"
  )
  expect_match(
    paste(capture.output(print(sensitivity)), collapse = "\n"),
    "No automatic stability classification",
    fixed = TRUE
  )
  expect_match(
    paste(capture.output(print(summary_object)), collapse = "\n"),
    "Grid comparisons",
    fixed = TRUE
  )

  table <- apa_table(sensitivity, digits = 5)
  expect_s3_class(table, "apa_table")
  expect_identical(table$which, "summary")
  expect_identical(nrow(table$table), 2L)
  expect_identical(
    as.data.frame(sensitivity),
    as.data.frame(sensitivity$summary, stringsAsFactors = FALSE)
  )
})

test_that("same-data and scope contracts fail closed", {
  fixture <- gpcm_quadrature_fixture()
  fit <- fixture$fit
  candidate <- fixture$sensitivity$fits$q7
  expect_true(mfrmr_gqs_same_prepared_data(fit, candidate))

  changed <- candidate
  changed$prep$data$score_k[1L] <- changed$prep$data$score_k[1L] + 1L
  expect_false(mfrmr_gqs_same_prepared_data(fit, changed))

  reordered <- candidate
  reordered$prep$data <- reordered$prep$data[
    nrow(reordered$prep$data):1L, , drop = FALSE
  ]
  expect_true(mfrmr_gqs_same_prepared_data(fit, reordered))

  expect_error(
    mfrmr_gqs_validate(
      fit, fixture$data, c(7L, 9L), c(-4, 4), 41L
    ),
    "must include the reference fit's quadrature count",
    fixed = TRUE
  )
  interaction_fit <- fit
  interaction_fit$config$interaction_specs <- list("Rater:Criterion" = list())
  expect_error(
    mfrmr_gqs_validate(
      interaction_fit, fixture$data, c(5L, 7L), c(-4, 4), 41L
    ),
    "facet interactions is not yet implemented",
    fixed = TRUE
  )
  regression_fit <- fit
  regression_fit$population$source <- "user_supplied"
  expect_error(
    mfrmr_gqs_validate(
      regression_fit, fixture$data, c(5L, 7L), c(-4, 4), 41L
    ),
    "user-supplied latent-regression",
    fixed = TRUE
  )
})

test_that("refit conditions are retained without error suppression", {
  expect_warning(
    captured <- mfrmr_gqs_capture_conditions({
      warning("quadrature warning remains visible")
      1L
    }),
    "quadrature warning remains visible",
    fixed = TRUE
  )
  expect_identical(captured$value, 1L)
  expect_identical(captured$warnings, "quadrature warning remains visible")
  expect_error(
    mfrmr_gqs_capture_conditions(stop("quadrature error remains visible")),
    "quadrature error remains visible",
    fixed = TRUE
  )
})

test_that("the public diagnostic has no machine-identity dependency", {
  namespace <- asNamespace("mfrmr")
  targets <- c(
    "gpcm_mml_quadrature_sensitivity",
    grep("^mfrmr_gqs_", ls(namespace, all.names = TRUE), value = TRUE),
    "as.data.frame.mfrm_quadrature_sensitivity",
    "summary.mfrm_quadrature_sensitivity",
    "print.mfrm_quadrature_sensitivity",
    "print.summary.mfrm_quadrature_sensitivity"
  )
  targets <- unique(targets)
  expect_true(all(vapply(
    targets,
    exists,
    logical(1),
    envir = namespace,
    inherits = FALSE
  )))
  runtime_text <- paste(vapply(
    targets,
    function(name) paste(deparse(get(name, envir = namespace)), collapse = "\n"),
    character(1)
  ), collapse = "\n")
  expect_false(grepl(
    "sha-?256|md5|digest::",
    runtime_text,
    ignore.case = TRUE
  ))
})
