test_that("reporting_checklist returns a bundle with checklist coverage tables", {
  dat <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "overall"))
  bias <- suppressWarnings(estimate_bias(
    fit,
    diag,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))

  chk <- reporting_checklist(fit, diagnostics = diag, bias_results = list(rater_criterion = bias))

  expect_s3_class(chk, "mfrm_reporting_checklist")
  expect_true(is.data.frame(chk$checklist))
  expect_true(is.data.frame(chk$summary))
  expect_true(is.data.frame(chk$section_summary))
  expect_true(all(c(
    "Section", "Item", "Available", "DraftReady", "ReadyForAPA", "Severity",
    "Priority", "SourceComponent", "Detail", "NextAction"
  ) %in% names(chk$checklist)))
  expect_identical(chk$checklist$DraftReady, chk$checklist$ReadyForAPA)
  expect_true(any(chk$checklist$Item == "PCA of residuals"))
  expect_true(any(chk$checklist$Item == "Facet pairs tested"))
  expect_true(chk$checklist$Available[chk$checklist$Item == "PCA of residuals"][1])
  expect_false(chk$checklist$ReadyForAPA[chk$checklist$Item == "95% confidence intervals"][1])
  expect_false(chk$checklist$ReadyForAPA[chk$checklist$Item == "Separation / strata / reliability"][1])
  expect_true(any(nzchar(chk$checklist$NextAction)))
  expect_true(all(c("DraftReady", "ReadyForAPA", "NeedsDraftWork", "NeedsAction") %in% names(chk$summary)))

  s_chk <- summary(chk)
  expect_s3_class(s_chk, "summary.mfrm_bundle")
  expect_true(is.data.frame(s_chk$summary))
  expect_true(nrow(s_chk$summary) > 0)
})

test_that("reporting_checklist surfaces non-numeric bias screening statistics", {
  dat <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "overall"))
  bias <- suppressWarnings(estimate_bias(
    fit,
    diag,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))
  bias$table$t <- "not-a-number"

  chk <- reporting_checklist(fit, diagnostics = diag, bias_results = list(rater_criterion = bias))
  row <- chk$checklist[chk$checklist$Item == "Screen-positive interactions", , drop = FALSE]

  expect_match(row$Detail[1], "non-numeric screening statistics", fixed = TRUE)
  expect_false(row$ReadyForAPA[1])
})

test_that("reporting_checklist surfaces failed bias-collection pairs", {
  dat <- load_mfrmr_data("example_bias")
  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score",
    method = "JML",
    maxit = 20
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "overall"))
  bias <- suppressWarnings(estimate_bias(
    fit,
    diag,
    facet_a = "Rater",
    facet_b = "Criterion",
    max_iter = 2
  ))
  bias_collection <- structure(
    list(
      by_pair = list(rater_criterion = bias),
      errors = data.frame(
        Interaction = "Task x Criterion",
        Facets = "Task x Criterion",
        Error = "forced pair failure",
        stringsAsFactors = FALSE
      )
    ),
    class = c("mfrm_bias_collection", "mfrm_bundle", "list")
  )

  chk <- reporting_checklist(fit, diagnostics = diag, bias_results = bias_collection)
  row <- chk$checklist[chk$checklist$Item == "Screen-positive interactions", , drop = FALSE]

  expect_match(row$Detail[1], "failed", fixed = TRUE)
  expect_false(row$ReadyForAPA[1])
  expect_identical(chk$settings$bias_error_count, 1L)
})
