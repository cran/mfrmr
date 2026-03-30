bias_collection_fixture <- local({
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

  dat <- mfrmr:::sample_mfrm_data(seed = 654)
  fit <- suppressWarnings(mfrmr::fit_mfrm(
    data = dat,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 20
  ))
  diagnostics <- suppressWarnings(mfrmr::diagnose_mfrm(fit, residual_pca = "none"))
  bias_all <- suppressWarnings(mfrmr::estimate_all_bias(
    fit,
    diagnostics = diagnostics,
    max_iter = 2
  ))

  list(
    fit = fit,
    diagnostics = diagnostics,
    bias_all = bias_all
  )
})

test_that("estimate_all_bias batches all modeled facet pairs", {
  bias_all <- bias_collection_fixture$bias_all

  expect_s3_class(bias_all, "mfrm_bias_collection")
  expect_true(is.data.frame(bias_all$summary))
  expect_true(is.list(bias_all$by_pair))
  expect_true(all(c("Interaction", "Rows", "Significant", "Kept") %in% names(bias_all$summary)))
  expect_identical(
    sort(bias_all$summary$Interaction),
    sort(c("Rater x Task", "Rater x Criterion", "Task x Criterion"))
  )
  expect_true(length(bias_all$by_pair) >= 1)
})

test_that("estimate_all_bias accepts explicit pair specifications", {
  fit <- bias_collection_fixture$fit
  diagnostics <- bias_collection_fixture$diagnostics

  bias_subset <- suppressWarnings(mfrmr::estimate_all_bias(
    fit,
    diagnostics = diagnostics,
    pairs = list(c("Rater", "Criterion")),
    max_iter = 2
  ))

  expect_s3_class(bias_subset, "mfrm_bias_collection")
  expect_identical(bias_subset$summary$Interaction, "Rater x Criterion")
  expect_true("Rater x Criterion" %in% names(bias_subset$by_pair) || !bias_subset$summary$Kept[1])
})
