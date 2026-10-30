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
  expect_true(all(c(
    "Interaction", "Rows", "ScreenPositive", "Significant",
    "FormalInferenceEligible", "PrimaryReportingEligible", "Kept"
  ) %in% names(bias_all$summary)))
  expect_identical(bias_all$summary$ScreenPositive, bias_all$summary$Significant)
  expect_true(all(!bias_all$summary$FormalInferenceEligible))
  expect_true(all(!bias_all$summary$PrimaryReportingEligible))
  expect_true(all(bias_all$summary$ReportingUse == "screening_only"))
  expect_identical(
    sort(bias_all$summary$Interaction),
    sort(c("Rater x Task", "Rater x Criterion", "Task x Criterion"))
  )
  expect_true(length(bias_all$by_pair) >= 1)

  one_summary <- summary(bias_all$by_pair[[1]])$overview
  expect_true(all(c(
    "ScreenPositive", "BonferroniScreenPositive", "HolmScreenPositive",
    "Significant", "BonferroniSignificant", "HolmSignificant",
    "FormalInferenceEligible", "PrimaryReportingEligible"
  ) %in% names(one_summary)))
  expect_identical(one_summary$ScreenPositive, one_summary$Significant)
  expect_identical(
    one_summary$BonferroniScreenPositive,
    one_summary$BonferroniSignificant
  )
  expect_identical(one_summary$HolmScreenPositive, one_summary$HolmSignificant)
  expect_false(one_summary$FormalInferenceEligible)
  expect_false(one_summary$PrimaryReportingEligible)
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

test_that("collection summaries preserve every cell and source without refitting", {
  x <- bias_collection_fixture$bias_all
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected refit"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"),
    estimate_bias = function(...) stop("Unexpected cell estimation"), .package = "mfrmr")
  review <- summary(x)
  expect_s3_class(review, "summary.mfrm_bias_collection")
  expect_identical(review$fit_identity, x$fit_identity)
  expect_identical(review$settings, x$settings)
  expect_identical(review$overview, x$summary)
  expect_equal(nrow(review$cells), sum(x$summary$Rows))
  for (label in names(x$by_pair)) {
    source <- as.data.frame(x$by_pair[[label]]$table)
    combined <- review$cells[review$cells$Interaction == label, names(source)]
    rownames(combined) <- NULL
    expect_identical(combined, source)
  }
  expect_output(print(x), "Conditional plug-in screens")
  expect_output(print(review), "ScreenUnavailable")
  expect_error(print(x, digits = NA_real_), "positive integer")
  path <- tempfile(fileext = ".rds")
  withr::defer(unlink(path))
  saveRDS(x, path)
  expect_identical(summary(readRDS(path)), review)
})

test_that("unavailable cell statistics are not counted as negative screens", {
  x <- bias_collection_fixture$bias_all
  label <- names(x$by_pair)[1]
  tab <- x$by_pair[[label]]$table
  tab[["S.E."]][] <- NA_real_
  tab$t[] <- NA_real_
  tab[["Bias Size"]][1] <- NA_real_
  tab$OptimizationStatus[1] <- "failed"
  x$by_pair[[label]]$table <- tab
  review <- summary(x)
  row <- review$overview[review$overview$Interaction == label, ]
  expect_equal(row$Points, nrow(tab) - 1L)
  expect_equal(row$ScreenAvailable, 0)
  expect_equal(row$ScreenUnavailable, nrow(tab))
  expect_true(is.na(row$ScreenPositive))
  expect_true(is.na(row$Significant))
  expect_identical(row$Status, "partial")
  cells <- subset(review$cells, Interaction == label)
  expect_true(all(is.na(cells$ScreenPositive)))
  expect_false(any(cells$ScreenAvailable))
  expect_identical(cells[["Bias Size"]], tab[["Bias Size"]])
})

test_that("failed pairs remain visible even without the duplicate error table", {
  local_mocked_bindings(estimate_bias = function(...) stop("Injected sparse-pair failure"),
    .package = "mfrmr")
  x <- estimate_all_bias(bias_collection_fixture$fit,
    diagnostics = bias_collection_fixture$diagnostics, keep_errors = FALSE)
  review <- summary(x)
  expect_equal(nrow(x$errors), 0L)
  expect_equal(nrow(review$errors), 3L)
  expect_true(all(review$overview$Status == "failed"))
  expect_true(all(is.na(review$overview$ScreenPositive)))
  expect_true(all(is.na(review$overview$ScreenAvailable)))
  expect_equal(nrow(review$cells), 0L)
  p <- plot(x, draw = FALSE)
  expect_identical(plot_data(p, "overview"), review$overview)
  expect_equal(nrow(plot_data(p, "table")), 0L)
  expect_true(all(grepl("failed", plot_data(p, "pair_labels"))))
  expect_error(plot(x, pair = 1, draw = FALSE), "Injected sparse-pair failure")
  skip_if_not_installed("ggplot2")
  expect_no_error(ggplot2::ggplot_build(as_ggplot(p)))
})

test_that("empty and missing retained tables remain distinct from successful screens", {
  x <- bias_collection_fixture$bias_all
  label <- names(x$by_pair)[1]
  x$by_pair[[label]]$table <- x$by_pair[[label]]$table[FALSE, ]
  x$summary$Rows[1] <- 0L
  expect_identical(summary(x)$overview$Status[1], "empty")
  expect_true(is.na(summary(x)$overview$ScreenPositive[1]))
  expect_error(plot(x, pair = 1, draw = FALSE), "No retained cells")
  x$by_pair[[label]] <- list()
  expect_identical(summary(x)$overview$Status[1], "empty")
  expect_error(plot(x, pair = 1, draw = FALSE), "No retained cells")
  x$summary$Rows[1] <- 9L
  x$by_pair[[label]] <- NULL
  expect_identical(summary(x)$overview$Status[1], "not_retained")
})

test_that("bias collection plots retain signed magnitudes, counts and selected-pair provenance", {
  x <- bias_collection_fixture$bias_all
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected refit"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"),
    estimate_bias = function(...) stop("Unexpected cell estimation"), .package = "mfrmr")
  device <- grDevices::dev.cur()
  p <- plot(x, draw = FALSE)
  expect_identical(grDevices::dev.cur(), device)
  expect_identical(p$name, "bias_collection")
  expect_identical(plot_data(p, "table")[["Bias Size"]], summary(x)$cells[["Bias Size"]])
  expect_identical(plot_data(p, "table")$ObsN, summary(x)$cells$ObsN)
  expect_identical(plot_data(p, "fit_identity"), x$fit_identity)
  expect_identical(plot_data(p, "preset"), "standard")
  expect_equal(length(plot_data(p, "pair_labels")), nrow(x$summary))
  one <- plot(x, pair = 1, type = "ranked", draw = FALSE)
  direct <- plot(x$by_pair[[1]], type = "ranked", draw = FALSE)
  expect_identical(plot_data(one, "ranked_table"), plot_data(direct, "ranked_table"))
  expect_identical(plot_data(one, "collection_pair"), names(x$by_pair)[1])
  expect_identical(plot_data(one, "fit_identity"), x$fit_identity)
  expect_identical(plot(x, pair = names(x$by_pair)[1], draw = FALSE)$name, "table13_bias")
  expect_error(plot(x, pair = "missing", draw = FALSE), "exact Interaction")
  expect_error(plot(x, pair = 0, draw = FALSE), "exact Interaction")
  expect_error(plot(x, type = "ranked", draw = FALSE), "Select `pair`")
  expect_error(plot(x, top_n = 2, draw = FALSE), "must be empty")
  skip_if_not_installed("ggplot2")
  gg <- ggplot2::ggplot_build(as_ggplot(p))
  table <- plot_data(p, "table")
  expect_equal(gg$data[[2]]$x, table[["Bias Size"]][table$PointAvailable])
  expect_equal(gg$data[[2]]$y, table$PlotY[table$PointAvailable])
  expect_error(as_ggplot(p, component = "overview"), "complete overview")
})

test_that("duplicate pair specifications cannot overwrite retained results", {
  fit <- bias_collection_fixture$fit
  expect_error(resolve_bias_collection_pairs(fit,
    pairs = list(c("Rater", "Task"), c("Rater", "Task"))), "only once")
  expect_error(resolve_bias_collection_pairs(fit,
    pairs = list(c("Rater", "Task"), c("Task", "Rater"))), "only once")
})

test_that("collection presentation refuses mixed fitted-analysis identities", {
  x <- bias_collection_fixture$bias_all
  x$by_pair[[1]]$fit_identity$fingerprint <- "different calibration"
  expect_error(summary(x), "mismatch")
  expect_error(plot(x, draw = FALSE), "mismatch")
  expect_error(plot(x, pair = 1, draw = FALSE), "mismatch")
  expect_error(mfrm_report(x), "mismatch")
})

test_that("computed collections feed reports without repeating estimation", {
  x <- bias_collection_fixture$bias_all
  s <- summary(x)
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected refit"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"),
    estimate_bias = function(...) stop("Unexpected cell estimation"))
  r <- mfrm_report(x)
  expect_identical(r$tables$overview, s$overview)
  expect_identical(r$tables$bias_cells, s$cells)
  expect_identical(r$source$fit_identity, x$fit_identity)
  expect_identical(plot_data(plot(r$source, draw = FALSE), "table")[["Bias Size"]],
    r$tables$bias_cells[["Bias Size"]])
})

test_that("explicit Person pairs do not require automatic Person enumeration", {
  one_facet_fit <- list(config = list(facet_names = "Rater"))

  resolved <- mfrmr:::resolve_bias_collection_pairs(
    one_facet_fit,
    pairs = list(c("Person", "Rater")),
    include_person = FALSE
  )

  expect_identical(resolved[[1]]$facets, c("Person", "Rater"))
  expect_identical(resolved[[1]]$label, "Person x Rater")
})

test_that("estimate_bias supports explicit Person screens with provenance", {
  fit <- bias_collection_fixture$fit
  diagnostics <- bias_collection_fixture$diagnostics

  expect_warning(
    person_bias <- mfrmr::estimate_bias(
      fit,
      diagnostics = diagnostics,
      facet_a = "Person",
      facet_b = "Rater",
      max_iter = 1
    ),
    regexp = "conditional plug-in likelihood screen",
    class = "mfrmr_person_bias_screen_warning"
  )

  expect_s3_class(person_bias, "mfrm_bias")
  expect_identical(person_bias$interaction_facets, c("Person", "Rater"))
  expect_true(isTRUE(person_bias$person_interaction))
  expect_identical(person_bias$person_source_column, "Person")
  expect_match(person_bias$person_interaction_note, "not a jointly fitted Person x facet")
  expect_gt(nrow(person_bias$table), 0L)
  expect_true(all(person_bias$table$FacetA == "Person"))
  expect_true(all(!person_bias$table$SupportsFormalInference))
  expect_true(all(!person_bias$table$FormalInferenceEligible))

  orientation <- as.data.frame(person_bias$orientation_review)
  expect_identical(orientation$Sign[orientation$Facet == "Person"], 1)
  expect_identical(
    orientation$Orientation[orientation$Facet == "Person"],
    "positive"
  )
  expect_true(isTRUE(person_bias$mixed_sign))

  aliased_note <- mfrmr:::bias_person_screen_note(list(
    config = list(source_columns = list(person = "Candidate"))
  ))
  expect_match(aliased_note, "source person column `Candidate`", fixed = TRUE)
})

test_that("estimate_all_bias warns once and runs Person pairs in fit order", {
  fit <- bias_collection_fixture$fit
  diagnostics <- bias_collection_fixture$diagnostics

  person_warnings <- list()
  bias_with_person <- withCallingHandlers(
    mfrmr::estimate_all_bias(
      fit,
      diagnostics = diagnostics,
      include_person = TRUE,
      max_iter = 1
    ),
    mfrmr_person_bias_screen_warning = function(cnd) {
      person_warnings[[length(person_warnings) + 1L]] <<- cnd
      invokeRestart("muffleWarning")
    }
  )
  expect_length(person_warnings, 1L)
  expect_match(
    conditionMessage(person_warnings[[1]]),
    "estimate_all_bias.*included.*Person"
  )

  expected <- c(
    "Person x Rater", "Person x Task", "Person x Criterion",
    "Rater x Task", "Rater x Criterion", "Task x Criterion"
  )
  expect_identical(bias_with_person$summary$Interaction, expected)
  expect_true(all(bias_with_person$summary$Rows > 0L))
  expect_true(all(bias_with_person$summary$Kept))
  expect_true(all(bias_with_person$summary$Error == ""))
  expect_identical(names(bias_with_person$by_pair), expected)
  expect_equal(nrow(bias_with_person$errors), 0L)
  expect_true(any(
    bias_with_person$settings$Setting == "person_pairs_requested" &
      bias_with_person$settings$Value == "TRUE"
  ))
})

test_that("MML Person screens use fitted EAP measures and posterior SEs", {
  dat <- mfrmr:::sample_mfrm_data(seed = 655)
  keep <- unique(dat$Person)[seq_len(12)]
  dat <- dat[dat$Person %in% keep, , drop = FALSE]
  fit <- suppressWarnings(mfrmr::fit_mfrm(
    dat,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "MML",
    quad_points = 7,
    maxit = 20
  ))
  diagnostics <- suppressWarnings(mfrmr::diagnose_mfrm(
    fit,
    residual_pca = "none"
  ))

  expect_warning(
    person_bias <- mfrmr::estimate_bias(
      fit,
      diagnostics,
      interaction_facets = c("Person", "Rater"),
      max_iter = 1
    ),
    class = "mfrmr_person_bias_screen_warning"
  )

  person_measures <- diagnostics$measures[
    diagnostics$measures$Facet == "Person",
    c("Level", "Estimate", "SE", "SE_Method"),
    drop = FALSE
  ]
  matched <- match(person_bias$table$FacetA_Level, person_measures$Level)
  expect_false(anyNA(matched))
  expect_equal(
    person_bias$table$FacetA_Measure,
    person_measures$Estimate[matched],
    tolerance = 1e-10
  )
  expect_equal(
    person_bias$table$FacetA_SE,
    person_measures$SE[matched],
    tolerance = 1e-10
  )
  expect_true(all(grepl(
    "Posterior SD",
    person_measures$SE_Method[matched],
    fixed = TRUE
  )))
})
