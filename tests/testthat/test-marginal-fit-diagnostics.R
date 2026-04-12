fit_for_marginal_fit_tests <- function(model = c("RSM", "PCM", "GPCM")) {
  model <- match.arg(model)
  dat <- mfrmr:::sample_mfrm_data(seed = 20260409)
  keep_persons <- unique(dat$Person)[1:12]
  dat <- dat[dat$Person %in% keep_persons, , drop = FALSE]

  fit_args <- list(
    data = dat,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "MML",
    model = model,
    maxit = 8,
    quad_points = 7
  )
  if (identical(model, "GPCM")) {
    fit_args$step_facet <- "Criterion"
  }

  fit <- suppressWarnings(do.call(mfrmr::fit_mfrm, fit_args))

  list(data = dat, fit = fit)
}

test_that("strict marginal fit diagnostics are available for MML RSM, PCM, and GPCM", {
  for (model in c("RSM", "PCM", "GPCM")) {
    obj <- fit_for_marginal_fit_tests(model)
    diag <- suppressWarnings(mfrmr::diagnose_mfrm(
      obj$fit,
      diagnostic_mode = "both",
      residual_pca = "none"
    ))

    expect_identical(diag$diagnostic_mode, "both")
    expect_true(is.data.frame(diag$diagnostic_basis))
    expect_true(isTRUE(diag$marginal_fit$available))
    expect_true(is.data.frame(diag$marginal_fit$overall$cell_stats))
    expect_true(is.data.frame(diag$marginal_fit$summary))
    expect_identical(as.character(diag$marginal_fit$summary$Model[1]), model)
    expect_true(is.data.frame(diag$marginal_fit$guidance))
    expect_true(isTRUE(diag$marginal_fit$pairwise$available))
    expect_true(is.data.frame(diag$marginal_fit$pairwise$pair_stats))
    expect_true(is.data.frame(diag$marginal_fit$pairwise$facet_summary))
    expect_true(all(diag$marginal_fit$guidance$InferenceTier == "exploratory"))
    expect_true(all(diag$marginal_fit$guidance$ReportingUse == "screening_only"))
    expect_true(all(diag$marginal_fit$summary$ReportingUse == "screening_only"))
    expect_true(
      "limited_information_inspired_marginal_screen" %in%
        diag$marginal_fit$guidance$LiteratureSeries
    )
    expect_true(
      any(grepl(
        "generalized_residual_logic",
        diag$marginal_fit$guidance$LiteratureSeries,
        fixed = TRUE
      ))
    )
    expect_true("legacy_residual_fit" %in% diag$diagnostic_basis$DiagnosticPath)
    expect_true("strict_marginal_fit" %in% diag$diagnostic_basis$DiagnosticPath)
    expect_identical(
      diag$diagnostic_basis$Status[diag$diagnostic_basis$DiagnosticPath == "legacy_residual_fit"],
      "computed"
    )
    expect_identical(
      diag$diagnostic_basis$Status[diag$diagnostic_basis$DiagnosticPath == "strict_marginal_fit"],
      "computed"
    )
    expect_true("posterior_predictive_follow_up" %in% diag$marginal_fit$guidance$Component)
    expect_identical(
      diag$marginal_fit$summary$PosteriorPredictiveFollowUp[1],
      "planned_not_implemented"
    )

    overall_cells <- diag$marginal_fit$overall$cell_stats
    expect_equal(
      sum(overall_cells$ObservedCount, na.rm = TRUE),
      sum(overall_cells$ExpectedCount, na.rm = TRUE),
      tolerance = 1e-6
    )

    diag_summary <- summary(diag)
    expect_s3_class(diag_summary, "summary.mfrm_diagnostics")
    expect_true("marginal_fit" %in% names(diag_summary))
    expect_true("diagnostic_basis" %in% names(diag_summary))
    expect_true("top_marginal_cells" %in% names(diag_summary))
    expect_true("marginal_pairwise" %in% names(diag_summary))
    expect_true("top_marginal_pairs" %in% names(diag_summary))
    expect_true("marginal_guidance" %in% names(diag_summary))
    expect_true(nrow(diag_summary$marginal_fit) == 1)
    expect_true(nrow(diag_summary$marginal_pairwise) >= 1)
    expect_true(nrow(diag_summary$marginal_guidance) >= 3)
    expect_true(nrow(diag_summary$diagnostic_basis) >= 4)
    expect_true("posterior_predictive_follow_up" %in% diag_summary$marginal_guidance$Component)

    printed <- paste(capture.output(print(diag_summary)), collapse = "\n")
    expect_match(printed, "Diagnostic basis guide", fixed = TRUE)
    expect_match(printed, "Strict marginal fit", fixed = TRUE)
    expect_match(printed, "Strict pairwise local dependence", fixed = TRUE)
    expect_match(printed, "Strict marginal guidance", fixed = TRUE)
    expect_match(printed, "planned_not_implemented", fixed = TRUE)

    rs <- mfrmr::rating_scale_table(obj$fit, diagnostics = diag)
    expect_true(all(c(
      "MarginalObservedCount",
      "MarginalExpectedCount",
      "MarginalResidualCount",
      "MarginalStdResidual",
      "MarginalFitFlag"
    ) %in% names(rs$category_table)))
    expect_true(isTRUE(rs$summary$MarginalFitAvailable[1]))

    cs <- mfrmr::category_structure_report(obj$fit, diagnostics = diag)
    expect_true(all(c(
      "MarginalObservedCount",
      "MarginalExpectedCount",
      "MarginalResidualCount",
      "MarginalStdResidual",
      "MarginalFitFlag"
    ) %in% names(cs$category_table)))
  }
})

test_that("legacy-only diagnostics keep strict path explicit but not computed", {
  obj <- fit_for_marginal_fit_tests("RSM")
  diag <- suppressWarnings(mfrmr::diagnose_mfrm(
    obj$fit,
    diagnostic_mode = "legacy",
    residual_pca = "none"
  ))

  expect_identical(diag$diagnostic_mode, "legacy")
  expect_true(is.data.frame(diag$diagnostic_basis))
  expect_false(isTRUE(diag$marginal_fit$available))
  expect_identical(
    diag$diagnostic_basis$Status[diag$diagnostic_basis$DiagnosticPath == "legacy_residual_fit"],
    "computed"
  )
  expect_identical(
    diag$diagnostic_basis$Status[diag$diagnostic_basis$DiagnosticPath == "strict_marginal_fit"],
    "not_requested"
  )

  diag_summary <- summary(diag)
  expect_true("diagnostic_basis" %in% names(diag_summary))
  expect_true("not_requested" %in% diag_summary$diagnostic_basis$Status)
})
