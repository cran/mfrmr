console_output_internal_terms <- c(
  "internal_case",
  "structural_design",
  "planned_not_implemented",
  "source-backed",
  "scaffold",
  "exact-overlap",
  "legacy-compatible",
  "InternalMethod",
  "current release",
  "later release",
  "follow_up_needed",
  "model_based",
  "screening_only",
  "reviewable_fit",
  "df_convention_difference",
  "not_supplied",
  "synthetic_latent_regression",
  "bundled_synthetic_example",
  "fit_standardization_review",
  "intended_use",
  "input_source"
)

expect_public_console_output <- function(lines, person_ids = character(0)) {
  person_ids <- as.character(person_ids)
  person_ids <- person_ids[!is.na(person_ids) & nzchar(person_ids)]
  expect_lte(max(nchar(lines), 0L), 80L)
  expect_false(any(vapply(
    console_output_internal_terms,
    function(term) any(grepl(tolower(term), tolower(lines), fixed = TRUE)),
    logical(1)
  )))
  expect_false(any(vapply(
    person_ids,
    function(id) any(grepl(id, lines, fixed = TRUE)),
    logical(1)
  )))
}

test_that("fit output leads with a plain-language decision", {
  fit <- make_toy_fit(method = "JML", maxit = 25)
  fit_lines <- capture.output(print(fit))
  brief <- summary(fit, profile = "fit", detail = "brief")
  brief_lines <- capture.output(print(brief))

  expect_identical(
    names(brief$decision),
    c(
      "Interpretation", "FormalInference", "FitReadiness", "Why",
      "NextAction"
    )
  )
  expect_identical(brief$decision$Interpretation, "Do not interpret this fit")
  expect_identical(brief$decision$FormalInference, "No")
  expect_match(brief$decision$Why, "Numerical convergence failed", fixed = TRUE)
  expect_lt(
    which(grepl("^Decision$", fit_lines)),
    which(grepl("^  Scale:", fit_lines))
  )
  expect_lt(
    which(grepl("^Decision$", brief_lines)),
    which(grepl("^Visual workflow", brief_lines))
  )
  expect_false(any(grepl("^  Summary status:", fit_lines)))
  expect_false(any(grepl("^  Key warning:", fit_lines)))
  expect_public_console_output(
    fit_lines,
    person_ids = unique(as.character(fit$facets$person$Person))
  )
})

test_that("major beginner-facing summaries respect the 80-column privacy contract", {
  old_options <- options(width = 80L)
  on.exit(options(old_options), add = TRUE)

  fit <- make_toy_fit(method = "JML", maxit = 25)
  diagnostics <- make_toy_diagnostics(
    fit,
    diagnostic_mode = "legacy",
    residual_pca = "none"
  )
  results <- suppressWarnings(mfrm_results(
    fit,
    include = c("fit", "diagnostics", "tables", "precision", "reporting")
  ))
  facets_review <- suppressWarnings(facets_fit_review(
    fit,
    diagnostics = diagnostics
  ))
  conquest_bundle <- suppressWarnings(build_conquest_overlap_bundle(
    quad_points = 3,
    maxit = 5
  ))

  person_ids <- unique(as.character(fit$facets$person$Person))
  conquest_person_ids <- unique(as.character(conquest_bundle$mfrmr_case_eap$Person))
  outputs <- list(
    diagnostics = list(
      lines = capture.output(print(summary(diagnostics))),
      ids = person_ids
    ),
    results_brief = list(
      lines = capture.output(print(summary(results, view = "brief"))),
      ids = person_ids
    ),
    facets_review = list(
      lines = capture.output(print(summary(facets_review))),
      ids = person_ids
    ),
    conquest_bundle = list(
      lines = capture.output(print(summary(conquest_bundle))),
      ids = conquest_person_ids
    )
  )

  diagnostics_lines <- outputs$diagnostics$lines
  results_lines <- outputs$results_brief$lines
  expect_true(any(grepl("^Decision$", diagnostics_lines)))
  expect_true(any(grepl("Formal inference:", diagnostics_lines, fixed = TRUE)))
  expect_true(any(grepl("^Decision$", results_lines)))
  expect_true(any(grepl("Formal inference:", results_lines, fixed = TRUE)))

  conquest_lines <- outputs$conquest_bundle$lines
  expect_false(conquest_bundle$summary$MfrmrInferenceReady[[1]])
  expect_true(any(grepl("mfrmr fit status", conquest_lines, fixed = TRUE)))
  expect_true(any(grepl("Inference ready", conquest_lines, fixed = TRUE)))
  expect_true(any(grepl("not inference-ready", conquest_lines, fixed = TRUE)))
  expect_true(any(grepl("Resolve the mfrmr convergence issue", conquest_lines, fixed = TRUE)))
  expect_false(any(grepl("Run the generated command locally", conquest_lines, fixed = TRUE)))

  for (output in outputs) {
    expect_public_console_output(output$lines, output$ids)
  }
})
