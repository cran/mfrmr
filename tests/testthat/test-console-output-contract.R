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

test_that("fit summaries keep one inference decision and readable scoring restrictions", {
  fit <- make_toy_fit(method = "MML", maxit = 25)
  s <- summary(fit, detail = "full", compute = "never")
  # Reproduce the former contradiction: fit readiness alone printed Ready
  # even when the precision-aware decision explicitly withheld inference.
  s$overview$InferenceReady <- TRUE
  s$decision$FormalInference <- "No"
  s$decision$Why <- "Precision assumptions require review (tier: model_based)."
  s$overview$Converged <- TRUE
  s$overview$ConvergenceCode <- 1L
  original <- serialize(s, NULL)
  lines <- capture.output(returned <- print(s))
  text <- paste(lines, collapse = " ")
  expect_length(grep("Formal inference:", lines, fixed = TRUE), 1L)
  expect_match(text, "Formal inference: No", fixed = TRUE)
  expect_match(text, "Optimizer code 0: No", fixed = TRUE)
  expect_false(grepl("legacy_mml|model_based|optimizer_gradient|ready_for_diagnostic_interpretation|source_fit_ready|not_applicable|LegacyConvergedBasis", text))
  expect_match(text, "fixed standard normal Person distribution", fixed = TRUE)
  expect_identical(serialize(s, NULL), original)
  expect_identical(returned, s)

  s$readiness$Status[1] <- "NEW_INTERNAL_STATE"
  s$required_visual$InterpretationStatus[1] <- "NEW_INTERNAL_STATE"
  lines <- capture.output(print(s))
  expect_false(any(grepl("NEW_INTERNAL_STATE", lines, fixed = TRUE)))
  expect_true(any(grepl("Review the stored result", lines, fixed = TRUE)))

  jml <- summary(make_toy_fit(method = "JML", maxit = 25), detail = "full", compute = "never")
  text <- gsub("[[:space:]]+", " ", paste(capture.output(print(jml)), collapse = " "))
  expect_match(text, "does not estimate a normal population distribution", fixed = TRUE)
  expect_false(grepl("legacy_mml|not_applicable_jml|source_fit_not_ready", text))
})

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

test_that("MML console output states its integration and population contracts", {
  fit <- make_toy_fit(method = "MML", maxit = 25)

  for (lines in list(
    capture.output(print(fit)),
    capture.output(print(summary(fit, detail = "brief")))
  )) {
    text <- gsub("[[:space:]]+", " ", paste(lines, collapse = " "))
    expect_match(text, "fixed non-adaptive Gauss-Hermite", fixed = TRUE)
    expect_match(text, "q=31", fixed = TRUE)
    expect_match(text, "latent dimensions=1", fixed = TRUE)
    expect_match(text, "Population identification: fixed N(0,1)", fixed = TRUE)
    expect_match(text, "discrimination=1", fixed = TRUE)
    contract_start <- grep("MML engine:", lines, fixed = TRUE)
    contract_end <- grep("Population identification:", lines, fixed = TRUE)
    expect_length(contract_start, 1L)
    expect_length(contract_end, 1L)
    expect_public_console_output(lines[contract_start:contract_end])
  }
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
  diagnostics_summary <- summary(diagnostics)
  diagnostics_original <- serialize(diagnostics_summary, NULL)
  outputs <- list(
    diagnostics = list(
      lines = capture.output(printed_diagnostics <- print(diagnostics_summary)),
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
  expect_identical(printed_diagnostics, diagnostics_summary)
  expect_identical(serialize(diagnostics_summary, NULL), diagnostics_original)
  diagnostics_text <- gsub("[[:space:]]+", " ", paste(diagnostics_lines, collapse = " "))
  for (reason in unique(diagnostics_summary$overall_fit$FitReason)) {
    expect_match(diagnostics_text, reason, fixed = TRUE)
  }
  expect_true(any(grepl("^Decision$", diagnostics_lines)))
  expect_true(any(grepl("Formal inference:", diagnostics_lines, fixed = TRUE)))
  expect_true(any(grepl("^Decision$", results_lines)))
  expect_true(any(grepl("Formal inference:", results_lines, fixed = TRUE)))

  conquest_lines <- outputs$conquest_bundle$lines
  expect_false(conquest_bundle$summary$MfrmrInferenceReady[[1]])
  expect_true(any(grepl("mfrmr fit status", conquest_lines, fixed = TRUE)))
  expect_true(any(grepl("Inference ready", conquest_lines, fixed = TRUE)))
  expect_true(any(grepl("not inference-ready", conquest_lines, fixed = TRUE)))
  expect_true(any(grepl("Inspect why the mfrmr fit is not inference-ready", conquest_lines, fixed = TRUE)))
  expect_false(any(grepl("Run the generated command locally", conquest_lines, fixed = TRUE)))

  for (output in outputs) {
    expect_public_console_output(output$lines, output$ids)
  }
  display <- diagnostics_summary
  uses <- c("approximate_summary", "review_before_reporting", "screening_only")
  display$reliability$ReliabilityUse <- rep(uses, length.out = nrow(display$reliability))
  for (detail in c("brief", "full")) {
    display$detail <- detail
    before <- serialize(display, NULL)
    text <- paste(capture.output(print(display)), collapse = " ")
    expect_false(any(vapply(uses, grepl, logical(1), x = text, fixed = TRUE)))
    for (label in c("Approximate summary", "Review before reporting", "Screening only"))
      expect_match(text, label, fixed = TRUE)
    expect_identical(serialize(display, NULL), before)
  }
})
