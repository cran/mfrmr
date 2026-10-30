test_that("anchor, linking and bias reports keep complete evidence without fitting", {
  dat <- load_mfrmr_data("example_core")
  anchors <- data.frame(Facet = "Criterion", Level = unique(dat$Criterion)[1], Value = .123456789)
  a <- review_mfrm_anchors(dat, "Person", c("Rater", "Criterion"), "Score", anchors = anchors)
  linking <- build_linking_review(anchor_review = a)
  bias <- structure(list(facet_a = "Rater", facet_b = "Criterion",
    table = data.frame(FacetA_Level = paste0("r", 1:12), FacetB_Level = "c1",
      `Bias Size` = seq(-.5, .6, .1), `S.E.` = .2, t = c(NA, 1:11),
      `Prob.` = c(NA, rep(.04, 11)), `Obs-Exp Average` = .1,
      BiasStatus = "available", ScreenStatus = "screening_only", FitStatus = "available",
      check.names = FALSE)), class = "mfrm_bias")
  before <- serialize(list(a, linking, bias), NULL)
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected fitting"),
    estimate_bias = function(...) stop("Unexpected bias estimation"))
  ar <- mfrm_report(a)
  expect_identical(ar$tables$anchors, as.data.frame(a$anchors))
  expect_identical(ar$tables$anchor_issues, as.data.frame(a$issues))
  expect_match(paste(ar$notes, collapse = " "), "does not establish invariance")
  lr <- mfrm_report(linking)
  expect_identical(lr$tables$status, as.data.frame(linking$status))
  expect_false(linking$overview$DriftAvailable)
  expect_match(paste(lr$notes, collapse = " "), "Missing drift")
  br <- mfrm_report(bias)
  expect_identical(br$tables$bias_cells, bias$table)
  expect_equal(nrow(br$tables$top_rows), 12)
  expect_true(is.na(br$tables$bias_cells$`Prob.`[1]))
  expect_false(br$tables$overview$SupportsFormalInference)
  for (r in list(ar, lr, br)) {
    expect_identical(build_summary_table_bundle(summary(r, top_n = 1), include_empty = TRUE)$tables, r$tables)
    expect_true(all(build_summary_table_bundle(r)$table_index$Role %in%
      mfrmr:::summary_table_bundle_appendix_role_registry()$Role))
    expect_match(r$markdown, r$title, fixed = TRUE)
  }
  expect_identical(serialize(list(a, linking, bias), NULL), before)
  expect_error(mfrm_report(a, intervals = TRUE), "only to G/D")
  expect_error(mfrm_report(a, style = "rater"), "requires mfrm_results")
})

test_that("declarations survive G/D reports, appendix export and saving", {
  data("mfrmr_example_gtheory_crossed", package = "mfrmr", envir = environment())
  g <- mfrm_g_study(mfrmr_example_gtheory_crossed, c("Content", "Organization"))
  d <- mfrm_d_study(g, data.frame(Raters = 2, Tasks = 3))
  context <- c(question = "Compare scoring plans", score_direction = "Higher means harder to understand",
    transformation = "No reversal; do not infer favorable performance from the variable name")
  plain <- mfrm_report(d)
  r <- mfrm_report(d, study_context = context)
  expect_identical(r$tables[names(plain$tables)], plain$tables)
  expect_identical(r$study_context$Value, unname(context))
  expect_identical(mfrm_report(d, output = "tables", study_context = context), r$tables)
  expect_match(mfrm_report(d, output = "markdown", study_context = context), context[2], fixed = TRUE)
  h <- mfrm_report(d, output = "html", study_context = context)
  expect_match(paste(readLines(h$path, warn = FALSE), collapse = ""), context[2], fixed = TRUE)
  expect_identical(build_summary_table_bundle(summary(r, top_n = 1))$tables$study_context, r$study_context)
  out <- tempfile(); dir.create(out)
  export_summary_appendix(r, output_dir = out)
  files <- list.files(out, "study_context.*csv$", full.names = TRUE)
  expect_length(files, 1)
  expect_identical(read.csv(files)$Value, unname(context))
  saveRDS(r, file.path(out, "report.rds"))
  expect_identical(readRDS(file.path(out, "report.rds")), r)
  for (bad in list(character(), "unnamed", c(a = NA_character_), c(a = " "),
      c(a = "one", a = "two"), stats::setNames("value", ""), list(a = "text")))
    expect_error(mfrm_report(d, study_context = bad), "study_context")
})

test_that("collection reports retain the search, exact cells and unavailable outcomes", {
  cell <- data.frame(FacetA_Level = paste0("r", 1:5), FacetB_Level = "c1",
    `Bias Size` = c(-.123456789, .4, -.6, 0, NA_real_),
    `S.E.` = c(rep(.2, 4), NA_real_), t = c(-.617283945, 2, -3, 0, NA_real_),
    ObsN = c(5, 6, 7, 8, 9), OptimizationStatus = c(rep("ok", 4), "failed"),
    BiasStatus = c(rep("available", 4), "unavailable"),
    ScreenStatus = c(rep("screening_only", 4), "unavailable"),
    FitStatus = "available", check.names = FALSE)
  other <- cell[1:2, ]
  other[["Bias Size"]] <- c(.1, NA_real_)
  other[["S.E."]] <- other$t <- NA_real_
  other$OptimizationStatus[2] <- "failed"
  other$ScreenStatus <- "unavailable"
  pair <- list(table = cell, iteration = data.frame(Iteration = 1:2, MaxChange = c(.6, 0)),
    orientation_review = data.frame(Facet = c("Rater", "Criterion"), Sign = -1),
    chi_sq = data.frame(Statistic = NA_real_, Status = "unavailable"),
    optimization_failures = data.frame(FacetA_Level = "r5", Reason = "No finite optimum"))
  x <- structure(list(
    summary = data.frame(Interaction = c("Rater x Criterion", "Task x Criterion", "Rater x Task"),
      Facets = c("Rater, Criterion", "Task, Criterion", "Rater, Task"), Rows = c(5L, 2L, 0L),
      ScreenPositive = NA_integer_, Significant = NA_integer_,
      MeanAbsBias = NA_real_, MaxAbsT = NA_real_,
      FormalInferenceEligible = FALSE, PrimaryReportingEligible = FALSE,
      Error = c(NA_character_, NA_character_, "Sparse-pair failure")),
    by_pair = list(`Rater x Criterion` = pair, `Task x Criterion` = list(table = other)),
    errors = data.frame(), settings = data.frame(Setting = c("requested_pairs", "keep_errors"),
      Value = c("3", "FALSE"))), class = "mfrm_bias_collection")
  before <- serialize(x, NULL)
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected fitting"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"),
    estimate_bias = function(...) stop("Unexpected bias estimation"))
  r <- mfrm_report(x, study_context = c(score_direction = "Higher means harder to understand"))
  expect_identical(r$tables$overview$Points, c(4L, 1L, NA_integer_))
  expect_identical(r$tables$overview$ScreenAvailable, c(4L, 0L, NA_integer_))
  expect_equal(r$tables$overview$ScreenPositive, c(2, NA, NA))
  expect_identical(r$tables$overview$Status, c("partial", "partial", "failed"))
  expect_equal(r$tables$overview$MeanAbsBias[1], (.123456789 + .4 + .6) / 4)
  expect_identical(r$tables$bias_cells[["Bias Size"]], c(cell[["Bias Size"]], other[["Bias Size"]]))
  expect_identical(r$tables$bias_cells$ObsN, c(5, 6, 7, 8, 9, 5, 6))
  expect_identical(r$tables$bias_cells$Interaction, c(rep("Rater x Criterion", 5), rep("Task x Criterion", 2)))
  expect_identical(r$tables$errors$Error, "Sparse-pair failure")
  expect_identical(r$tables$settings, x$settings)
  expect_identical(r$tables$iteration$MaxChange, c(.6, 0))
  expect_identical(r$tables$orientation_review$Sign, c(-1, -1))
  expect_identical(r$tables$optimization_failures$Reason, "No finite optimum")
  expect_true(is.na(r$tables$chi_sq$Statistic))
  expect_match(paste(r$notes, collapse = " "), "not prespecified confirmatory")
  expect_match(r$markdown, "Sparse-pair failure", fixed = TRUE)
  expect_identical(build_summary_table_bundle(summary(r, top_n = 1), include_empty = TRUE)$tables, r$tables)
  expect_identical(build_summary_table_bundle(x, top_n = 1, include_empty = TRUE)$tables$bias_cells,
    r$tables$bias_cells)
  expect_s3_class(apa_table(x, which = "bias_cells"), "apa_table")
  expect_true(all(build_summary_table_bundle(r)$table_index$Role %in%
    mfrmr:::summary_table_bundle_appendix_role_registry()$Role))
  for (style in c("apa", "technical", "reviewer", "validation"))
    expect_identical(mfrm_report(x, style = style, output = "tables"),
      r$tables[setdiff(names(r$tables), "study_context")])
  out <- tempfile(); dir.create(out)
  withr::defer(unlink(out, recursive = TRUE))
  export_summary_appendix(r, output_dir = out)
  csv <- list.files(out, "bias_cells.*csv$", full.names = TRUE)
  expect_length(csv, 1)
  expect_equal(read.csv(csv, check.names = FALSE)[["Bias Size"]], r$tables$bias_cells[["Bias Size"]], tolerance = 1e-14)
  h <- mfrm_report(x, output = "html")
  withr::defer(unlink(h$path))
  expect_match(paste(readLines(h$path, warn = FALSE), collapse = ""), "Sparse-pair failure", fixed = TRUE)
  path <- file.path(out, "report.rds")
  saveRDS(r, path)
  expect_identical(readRDS(path), r)
  expect_identical(serialize(x, NULL), before)

  x$by_pair <- list()
  x$summary$Error <- "No retained result"
  failed <- mfrm_report(x)
  expect_equal(nrow(failed$tables$errors), 3)
  expect_equal(nrow(failed$tables$bias_cells), 0)
  expect_true(all(is.na(failed$tables$overview$ScreenPositive)))
  expect_no_error(mfrm_report(x, output = "markdown"))
  expect_error(mfrm_report(x, intervals = TRUE), "only to G/D")
  expect_error(mfrm_report(x, style = "rater"), "requires mfrm_results")
})
