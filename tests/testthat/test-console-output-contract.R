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
