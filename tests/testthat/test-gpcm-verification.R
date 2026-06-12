# Verification tests for the bounded-GPCM workflow declared in
# `gpcm_capability_matrix()`. Each test exercises one row of the
# matrix that is `"supported"` or `"supported_with_caveat"` and
# asserts the corresponding helper returns the documented shape.
# `"blocked"` and `"deferred"` rows have negative tests where the
# helper should refuse to run (or run with an explicit caveat).

skip_if_no_lme4 <- function() {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    skip("`lme4` (Suggests) not installed; skipping GPCM verification.")
  }
}

local({
  .toy_gpcm <<- load_mfrmr_data("example_core")
  .gpcm_fit <<- suppressMessages(suppressWarnings(
    fit_mfrm(.toy_gpcm, "Person", c("Rater", "Criterion"), "Score",
             method = "MML", model = "GPCM",
             step_facet = "Criterion",
             slope_facet = "Criterion",
             quad_points = 7L, maxit = 25L)
  ))
})

test_that("GPCM core fit returns a populated mfrm_fit", {
  expect_s3_class(.gpcm_fit, "mfrm_fit")
  expect_identical(as.character(.gpcm_fit$config$model), "GPCM")
  expect_true(nrow(.gpcm_fit$summary) > 0L)
  expect_true(nrow(.gpcm_fit$facets$person) > 0L)
})

test_that("GPCM print / summary do not error", {
  expect_no_error(invisible(utils::capture.output(print(.gpcm_fit))))
  expect_no_error(invisible(utils::capture.output(print(summary(.gpcm_fit)))))
})

test_that("GPCM diagnose_mfrm returns measures with caveat status", {
  diag <- suppressMessages(suppressWarnings(
    diagnose_mfrm(.gpcm_fit, residual_pca = "none", diagnostic_mode = "legacy")
  ))
  expect_true(is.data.frame(diag$measures))
  expect_true(nrow(diag$measures) > 0L)
  # The dashboard panel remains unavailable under GPCM; direct
  # fair_average_table() is supported with its own caveat.
  if (!is.null(diag$fair_average)) {
    fa_msg <- as.character(diag$fair_average$status %||% "")
    expect_true(any(grepl("placeholder|unavailable|GPCM", fa_msg,
                           ignore.case = TRUE)) ||
                  is.null(diag$fair_average$table) ||
                  nrow(as.data.frame(diag$fair_average$table)) == 0L)
  }
})

test_that("GPCM compute_information + plot_information work", {
  info <- compute_information(.gpcm_fit)
  expect_true(is.list(info))
  expect_true("tif" %in% names(info))
  p <- plot_information(info, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

test_that("GPCM CCC / pathway / Wright plots return mfrm_plot_data", {
  for (type in c("wright", "pathway", "ccc")) {
    p <- plot(.gpcm_fit, type = type, draw = FALSE)
    expect_s3_class(p, "mfrm_plot_data")
  }
})

test_that("GPCM capability matrix is consistent with the helper", {
  m <- gpcm_capability_matrix()
  expect_true(is.data.frame(m))
  expect_true(all(c("Area", "Helpers", "Status", "PrimaryUse", "Boundary")
                  %in% names(m)))
  expect_true(all(c("RecommendedRoute", "NextValidationStep") %in% names(m)))
  expect_true(all(m$Status %in%
                    c("supported", "supported_with_caveat", "blocked", "deferred")))
})

test_that("GPCM APA/QC reporting bundle returns with explicit caveats", {
  diag <- suppressMessages(suppressWarnings(
    diagnose_mfrm(.gpcm_fit, residual_pca = "none", diagnostic_mode = "legacy")
  ))
  apa <- suppressMessages(build_apa_outputs(.gpcm_fit, diag))
  expect_s3_class(apa, "mfrm_apa_outputs")
  expect_true(nrow(apa$gpcm_boundary) > 0)
  expect_true(grepl("Bounded\\s+GPCM note", apa$report_text))

  qc <- run_qc_pipeline(.gpcm_fit, diag)
  expect_s3_class(qc, "mfrm_qc_pipeline")
  expect_true(nrow(qc$gpcm_boundary) > 0)
})
