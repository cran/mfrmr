dstudy_summary_example <- function() {
  ratings <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
  mfrm_d_study(mfrm_g_study(ratings, c("V", "W"), rater = NULL),
    data.frame(Tasks = c(6, 12)), c(V = .5, W = .5))
}

test_that("ordinary summaries stay compatible and interval summaries preserve all targets", {
  d <- dstudy_summary_example(); original <- d
  expect_identical(summary(d), d$coefficients)
  expect_error(summary(d, intervals = TRUE), "normal")
  expect_error(summary(d, parm = "G"), "only when")
  expect_error(summary(d, intervals = TRUE, parm = c("G", "G")), "distinct")
  expect_error(summary(d, intervals = NA), "must be FALSE")
  expect_error(summary(d, intervals = TRUE, typo = TRUE, assumption = "normal"), "must be empty")
  s <- summary(d, intervals = TRUE, assumption = "normal")
  expect_s3_class(s, "summary.mfrm_multivariate_d_study")
  expect_identical(d, original)
  expect_identical(s$source_design, d$gstudy$design)
  expect_identical(s$source_rows, d$gstudy$data_usage$counts)
  expect_identical(s$weights, d$weights)
  expect_identical(s$design_grid, d$design_grid)
  expect_identical(s$component_diagnostics, d$component_diagnostics)
  expect_identical(names(s$intervals), c("G", "Phi", "RelativeSEM", "AbsoluteSEM"))
  for (metric in names(s$intervals)) {
    ci <- confint(d, parm = metric, assumption = "normal")
    expect_identical(s$intervals[[metric]], ci)
    row <- s$table[s$table$Metric == metric, ]
    for (field in names(ci)) expect_identical(row[[field]], ci[[field]])
    expect_true(all(attr(ci, "notes") %in% s$notes))
  }
  expect_output(print(s), "Method assumptions and limitations")
  expect_error(summary(d, intervals = s$intervals, parm = "Phi"), "own metrics")
  expect_error(summary(d, intervals = s$intervals, assumption = "normal"), "must be empty")
})

test_that("paired results retain dependence and can be collected with individual plans", {
  d <- dstudy_summary_example()
  s <- summary(d, intervals = TRUE, reference = 1, assumption = "normal")
  expect_identical(s$table$Metric, paste(names(s$intervals), "difference"))
  expect_true(all(s$table$Reference == 1L))
  for (metric in names(s$intervals)) {
    expect_identical(s$intervals[[metric]], confint(d, parm = metric,
      reference = 1, assumption = "normal"))
  }
  ci <- confint(d, assumption = "normal")
  both <- summary(d, intervals = c(list(G = ci), s$intervals))
  expect_false(anyDuplicated(names(both$intervals)) > 0L)
  expect_true(all(c("G", "G difference") %in% both$table$Metric))
  expect_identical(both$intervals[[2]], s$intervals$G)
  p <- plot(both$intervals[[2]], draw = FALSE)
  expect_identical(plot_data(p, "table")$Estimate, s$intervals$G$Estimate)
})

test_that("saved summaries reject unrelated sources, targets and altered values", {
  d <- dstudy_summary_example(); ci <- confint(d, assumption = "normal")
  changed <- d; changed$gstudy$components$Person <- changed$gstudy$components$Person * 2
  expect_error(summary(changed, intervals = ci), "source summaries")
  scaled <- mfrm_d_study(d$gstudy, d$design_grid, 2*d$weights)
  expect_error(summary(scaled, intervals = ci), "weights")
  changed <- ci; changed$Tasks[1] <- 9
  expect_error(summary(d, intervals = changed), "facet counts")
  changed <- ci; changed$Estimate[1] <- .9
  expect_error(summary(d, intervals = changed), "point estimates")
  changed <- ci; changed$Available[1] <- FALSE
  expect_error(summary(d, intervals = changed), "availability")
  expect_error(summary(d, intervals = ci[1, ]), "one interval row")
  changed <- ci; attr(changed, "source_summary") <- NULL
  expect_error(summary(d, intervals = changed), "source summaries")
  changed <- confint(d, reference = 1, assumption = "normal")
  changed$ReferenceValue <- 0
  expect_error(summary(d, intervals = changed), "comparison values")
})

test_that("saved results keep unavailable bounds, methods and replay without calculation", {
  d <- dstudy_summary_example()
  ci <- confint(d, assumption = "normal")
  # Preserve a calculation-failure state without constructing new numerical evidence.
  ci$Lower[1] <- ci$Upper[1] <- NA_real_; ci$Available[1] <- FALSE
  ci$Status[1] <- "Unavailable for replay test"
  s <- summary(d, intervals = ci)
  expect_false(s$table$Available[1]); expect_true(is.na(s$table$Lower[1]))
  expect_identical(s$table$Status, ci$Status)
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(s, path)
  local_mocked_bindings(confint.mfrm_multivariate_d_study = function(...) stop("Unexpected calculation"))
  reopened <- readRDS(path)
  expect_identical(reopened, s)
  expect_identical(summary(d, intervals = reopened$intervals), s)
  expect_output(print(reopened), "Unavailable for replay test")
  p <- plot(d, metric = "G", intervals = reopened$intervals[[1]], draw = FALSE)
  expect_false(plot_data(p, "series")$IntervalAvailable[1])
})

test_that("crossed method comparisons and nested facet roles remain explicit", {
  g <- mfrm_g_study(mfrmr_example_gtheory_crossed, c("Content", "Organization"))
  d <- mfrm_d_study(g, data.frame(Raters = c(2,4), Tasks = 6))
  s <- summary(d, intervals = TRUE, parm = "G", score = "Content",
    method = c("mls", "approximate_f"), assumption = "normal")
  expect_setequal(s$table$MethodCode, c("mls", "approximate_f"))
  expect_true(all(s$table$Score == "Content"))
  expect_true("FormulaComplete" %in% names(s$table))
  expect_true(all(s$table$Validation == s$intervals$G$Validation))
  expect_true(all(s$table$Available))
  expect_false(any(s$table$FormulaComplete == FALSE, na.rm = TRUE))
  expect_false(any(grepl("Rows needing review", capture.output(print(s)), fixed = TRUE)))
  incomplete <- s$intervals$G
  incomplete$FormulaComplete[1] <- FALSE
  incomplete$Status[1] <- "Conservative enclosure includes undefined MLS widths"
  reviewed <- summary(d, intervals = incomplete)
  expect_output(print(reviewed), "Conservative enclosure includes undefined MLS widths")
  expect_identical(reviewed$table$Status, incomplete$Status)
  gn <- mfrm_g_study(mfrmr_example_gtheory_nested,
    c("Content", "Organization"), nesting = c(Rater = "Task"))
  dn <- mfrm_d_study(gn, data.frame(Raters = 3, Tasks = c(6,12)))
  sn <- summary(dn, intervals = TRUE, parm = c("G", "RelativeSEM"), assumption = "normal")
  expect_identical(sn$source_design$nesting, c(Rater = "Task"))
  expect_identical(attr(sn$intervals$G, "nesting_context")$future$ChildrenPerParent, c(3,3))
})

test_that("experimental REML summaries retain separate centers and source ratings", {
  skip_if_not_installed("lme4"); skip_if_not_installed("minqa")
  ratings <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
  g <- mfrm_g_study(ratings[-1, ], c("V", "W"), rater = NULL, method = "minque0")
  d <- mfrm_d_study(g, data.frame(Tasks = c(6,12)))
  s <- summary(d, intervals = TRUE, parm = "G", method = "reml_profile", assumption = "normal")
  expect_identical(s$table$Estimate, subset(d$coefficients, Score == "V")$G)
  expect_identical(s$table$IntervalCenter, s$intervals$G$IntervalCenter)
  expect_true(all(grepl("Experimental", s$table$Validation)))
  expect_false(s$source_design$complete)
  local_mocked_bindings(.mfrm_mvgt_profile_source = function(...) stop("Unexpected refit"))
  expect_identical(summary(d, intervals = s$intervals), s)
  changed <- d; changed$gstudy$data[1, "V"] <- changed$gstudy$data[1, "V"] + 1
  expect_error(summary(changed, intervals = s$intervals), "source ratings")
})

test_that("G/D print methods honor display options without calculation or mutation", {
  d <- dstudy_summary_example()
  s <- summary(d, intervals = TRUE, parm = "G", assumption = "normal")
  g_crossed <- mfrm_g_study(mfrmr_example_gtheory_crossed, c("Content", "Organization"))
  d_crossed <- mfrm_d_study(g_crossed, data.frame(Raters = c(2,3), Tasks = c(6,4)))
  comparison <- mfrm_multivariate_d_compare(d_crossed, reference = 1, assumption = "normal")
  local_mocked_bindings(confint.mfrm_multivariate_d_study = function(...) stop("Unexpected interval calculation"),
    mfrm_g_study = function(...) stop("Unexpected G-study"),
    mfrm_d_study = function(...) stop("Unexpected D-study"))
  for (object in list(d$gstudy, d, s, comparison)) {
    before <- serialize(object, NULL)
    lo <- capture.output(value <- withVisible(print(object, digits = 2, row.names = FALSE)))
    hi <- capture.output(print(object, digits = 8, row.names = FALSE))
    expect_false(identical(lo, hi))
    expect_false(value$visible)
    expect_identical(value$value, object)
    expect_identical(serialize(object, NULL), before)
    expect_no_error(capture.output(print(object, row.names = TRUE, quote = FALSE)))
  }
  expect_identical(summary(d$gstudy), d$gstudy$component_diagnostics)
  expect_identical(summary(comparison), comparison$comparisons)
  expect_error(summary(d$gstudy, intervals = TRUE), "must be empty")
  expect_error(summary(comparison, intervals = TRUE), "must be empty")
})

test_that("printed summaries identify confidence levels and nested future count units", {
  d <- dstudy_summary_example()
  s <- summary(d, intervals = list(
    confint(d, level = .8, assumption = "normal"),
    confint(d, level = .95, assumption = "normal")))
  output <- capture.output(print(s, digits = 3))
  expect_true(any(grepl("MethodCode", output, fixed = TRUE)))
  expect_true(any(grepl("\\bLevel\\b", output)))
  expect_true(any(grepl("0.80", output, fixed = TRUE)))
  expect_true(any(grepl("0.95", output, fixed = TRUE)))
  g <- mfrm_g_study(mfrmr_example_gtheory_nested,
    c("Content", "Organization"), nesting = c(Rater = "Task"))
  d <- mfrm_d_study(g, data.frame(Raters = 3, Tasks = c(6,12)))
  s <- summary(d, intervals = TRUE, parm = "G", assumption = "normal")
  for (object in list(d, s)) {
    expect_output(print(object), "Future Raters counts levels per Task")
    expect_output(print(object), "total across parents")
  }
})
