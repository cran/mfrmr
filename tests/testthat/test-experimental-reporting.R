test_that("two-family method notes survive reporting without interval computation", {
  fit <- readRDS(test_path("fixtures", "gmfrm-explicit-start.rds"))$fit
  before <- serialize(fit, NULL)
  local_mocked_bindings(
    fit_mfrm = function(...) stop("No refitting"),
    confint.mfrm_fit = function(...) stop("No intervals"),
    compute_mml_parameter_covariance = function(...) stop("No information"),
    mfrm_gpcm_response_evaluator = function(...) stop("No response calculation"),
    .package = "mfrmr")
  res <- mfrm_results(fit, include = "fit", compute = "never")
  note <- res$tables$gpcm_interpretation_notes$Note
  expect_match(note, "qualification is incomplete")
  expect_match(note, "curve/step intervals.*unavailable")
  expect_true(note %in% summary(res)$notes)
  expect_output(print(summary(res, view = "brief")), "qualification is incomplete")
  expect_null(res$gpcm_inference)
  reopened <- unserialize(serialize(res, NULL))
  report <- mfrm_report(reopened)
  expect_identical(report$tables$gpcm_interpretation_notes$Note, note)
  expect_match(report$markdown, "qualification is incomplete")
  expect_identical(report$tables$fitted_slopes$FittedSlope, fit$slopes$OptimizerEstimate)

  # Old result objects remain readable and receive the current report explanation.
  legacy <- reopened
  legacy$tables$gpcm_interpretation_notes <- NULL
  expect_identical(mfrm_report(legacy)$tables$gpcm_interpretation_notes$Note, note)

  path <- tempfile("experimental-notes-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  export_mfrm_results(reopened, output_dir = path, include = "tables",
    acknowledge_sensitive = TRUE)
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  notes <- files[grepl("gpcm_interpretation_notes[.]csv$", files)]
  expect_length(notes, 1L)
  expect_identical(read.csv(notes)$Note, note)
  expect_identical(serialize(fit, NULL), before)
})

test_that("rater table handoffs keep recipient-safe interpretation for both model routes", {
  fits <- list(
    pcm = readRDS(test_path("fixtures", "native-location-pcm-regression.rds"))$fit,
    gpcm = readRDS(test_path("fixtures", "gmfrm-explicit-start.rds"))$fit)
  local_mocked_bindings(
    fit_mfrm = function(...) stop("No refitting"),
    diagnose_mfrm = function(...) stop("No diagnostics"),
    confint.mfrm_fit = function(...) stop("No intervals"),
    compute_mml_parameter_covariance = function(...) stop("No information"),
    mfrm_gpcm_response_evaluator = function(...) stop("No response calculation"),
    .package = "mfrmr")
  for (fit in fits) {
    res <- mfrm_results(fit, include = "fit", compute = "never")
    res$notes <- "PRIVATE_SOURCE_NOTE"
    res$fit$notes <- "PRIVATE_FIT_NOTE"
    facet <- if (mfrm_has_product_slopes(fit)) fit$config$slope_facet[2] else "Rater"
    rater <- as.character(fit$facets$others$Level[fit$facets$others$Facet == facet])[1]
    make_sheet <- function(audience) mfrm_report(res, style = "rater", facet = facet,
      rater = rater, audience = audience, max_cases = 0)
    for (audience in c("rater", "researcher")) {
      sheet <- make_sheet(audience)
      tables <- mfrm_report(res, "rater", "tables", facet, rater,
        audience = audience, max_cases = 0)
      expect_identical(tables, sheet$tables)
      notes <- tables$interpretation_notes
      expect_identical(notes$Note[notes$Section == "review"], sheet$review)
      expect_identical(notes$Note[notes$Section == "guidance"], sheet$guidance)
      expect_identical(notes$Note[match(names(sheet$notes), notes$Section)], unname(sheet$notes))
      expect_match(notes$Note[notes$Section == "uncertainty"], "No saved")
      expect_match(notes$Note[notes$Section == "cases"], "omitted")
      expect_false(grepl("PRIVATE_", paste(capture.output(dput(sheet)), collapse = "\n")))

      # The added export table must not duplicate or change the rendered sheet.
      legacy <- sheet; legacy$tables$interpretation_notes <- NULL
      expect_identical(mfrm_rater_feedback_render(sheet), mfrm_rater_feedback_render(legacy))
      expect_identical(mfrm_rater_feedback_render(sheet, html = FALSE), sheet$markdown)
      file <- tempfile(fileext = ".csv")
      write.csv(notes, file, row.names = FALSE)
      expect_identical(read.csv(file), notes)
      unlink(file)
    }
    rater_sheet <- make_sheet("rater"); researcher_sheet <- make_sheet("researcher")
    numeric_tables <- setdiff(names(rater_sheet$tables), "interpretation_notes")
    expect_identical(rater_sheet$tables[numeric_tables], researcher_sheet$tables[numeric_tables])
    expect_false(identical(rater_sheet$tables$interpretation_notes,
      researcher_sheet$tables$interpretation_notes))
  }
})
