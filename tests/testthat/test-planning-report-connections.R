planning_report_fixture <- function() {
  saved <- readRDS(system.file("examples", "simulation-design.rds", package = "mfrmr"))
  s <- summary(saved$result)
  row <- s$design_summary
  row <- row[!duplicated(row$Facet), ]
  s$design_summary <- dplyr::bind_rows(lapply(1:3, function(i) {
    z <- row
    z$design_id <- paste0("D", i)
    z$n_person <- 10L * i
    z$MeanSeparation <- 3
    z$MeanReliability <- .9
    z$MeanSeverityRMSE <- .1
    z$MeanMisfitRate <- .01
    z$ConvergenceRate <- 1
    z$MaxRatings <- c(30, 40, 50)[i]
    z$MaxRatingsPerRater <- c(10, 20, NA)[i]
    z$MaxRaterComponents <- z$MaxCriterionComponents <- if (i == 3) NA_real_ else 1
    z$DisconnectedReps <- 0
    if (i == 1) z$MeanReliability[z$Facet == "Person"] <- .7
    if (i == 3) z$MeanReliability[z$Facet == "Person"] <- NA_real_
    z
  }))
  s
}

test_that("planning target and criterion outcomes are explicit without changing selection", {
  s <- planning_report_fixture()
  x <- recommend_mfrm_design(s, facets = "Person", max_ratings = 40)
  z <- summary(x)
  expect_identical(x$recommended$design_id, "D2")
  expect_identical(z$recommended, x$recommended)
  expect_identical(z$target$Target, "Person differentiation")
  expect_identical(z$thresholds$Source[z$thresholds$Setting == "min_reliability"], "default")
  expect_identical(z$thresholds$Source[z$thresholds$Setting == "max_ratings"], "user_supplied")
  row_ids <- z$design_table$DesignRow[match(c("D1", "D2", "D3"), z$design_table$design_id)]
  rel <- z$checks[z$checks$Criterion == "Reliability", ]
  expect_identical(rel$Status[match(row_ids, rel$DesignRow)], c("fail", "pass", "unavailable"))
  expect_equal(rel$Value[match(row_ids, rel$DesignRow)], c(.7, .9, NA))
  expect_equal(sum(z$design_table$Selected), 1L)
  expect_identical(z$design_table$design_id[z$design_table$Selected], "D2")
  expect_true(all(z$checks$Status[z$checks$Criterion == "Ratings per rater"] == "not_requested"))
  expect_true(all(z$checks$Status[z$checks$Criterion == "Facet coverage"] == "pass"))
  calibration <- recommend_mfrm_design(s)
  expect_identical(calibration$recommended$design_id, "D1")
  expect_true(all(summary(calibration)$target$Target == "Facet calibration"))
  expect_match(paste(capture.output(print(x)), collapse = " "), "Person differentiation")
  empty <- recommend_mfrm_design(s, facets = "Person", max_ratings = 0)
  expect_equal(nrow(empty$recommended), 0L)
  expect_false(any(summary(empty)$design_table$Selected))
  expect_match(paste(capture.output(print(empty)), collapse = " "), "No evaluated design")
  older <- x; older$threshold_sources <- NULL
  expect_true(all(summary(older)$thresholds$Source == "not_recorded"))
  expect_error(recommend_mfrm_design(s, min_reliability = c(.7, .8)), "one numeric threshold")
  expect_error(recommend_mfrm_design(s, min_reliability = NA_real_), "one numeric threshold")
  expect_error(summary(x, unused = TRUE), "must be empty")
  # Facet names are labels, not comma-delimited lists to parse for counts.
  named <- s
  named$design_summary$Facet[named$design_summary$Facet == "Rater"] <- "Judge, examiner"
  named$design_summary <- named$design_summary[!(named$design_summary$design_id == "D1" &
    named$design_summary$Facet == "Judge, examiner"), ]
  ns <- summary(recommend_mfrm_design(named, facets = c("Judge, examiner", "Criterion")))
  row <- ns$design_table$DesignRow[ns$design_table$design_id == "D1"]
  coverage <- ns$checks[ns$checks$DesignRow == row & ns$checks$Criterion == "Facet coverage", ]
  expect_equal(coverage$Value, 1)
  expect_equal(coverage$Threshold, 2)
  expect_identical(coverage$Status, "fail")
})

test_that("recommendation figures, reports and exports retain all candidates and unavailable values", {
  x <- recommend_mfrm_design(planning_report_fixture(), facets = "Person", max_ratings = 40)
  s <- summary(x)
  r <- mfrm_report(x, study_context = c(score_direction = "Higher is more comprehensible"))
  expect_identical(r$tables$checks, s$checks)
  expect_identical(r$tables$designs, s$design_table)
  expect_identical(r$tables$facets, s$facet_table)
  expect_identical(r$tables$thresholds, s$thresholds)
  for (input in list(x, r, summary(r, top_n = 1))) {
    bundle <- build_summary_table_bundle(input, top_n = 1)
    expect_identical(bundle$tables$checks, s$checks)
    expect_true(all(bundle$table_index$Role %in% summary_table_bundle_appendix_role_registry()$Role))
  }
  expect_s3_class(apa_table(x, which = "checks"), "apa_table")
  expect_error(mfrm_report(x, intervals = TRUE), "only to G/D")
  p <- plot(x, draw = FALSE)
  expect_identical(plot_data(p, "table")[names(s$checks)], s$checks)
  long <- plot_data(x, "designs", format = "long", value_cols = c("MinReliability", "MaxRatings"))
  expect_equal(nrow(long), 6L)
  expect_equal(long$Value[long$Variable == "MinReliability"], s$design_table$MinReliability)
  out <- tempfile(); dir.create(out)
  export_summary_appendix(r, output_dir = out, prefix = "planning")
  csv <- list.files(out, "checks.*csv$", full.names = TRUE)
  expect_length(csv, 1L)
  tab <- read.csv(csv)
  expect_equal(tab$Value, s$checks$Value, tolerance = 1e-12)
  expect_identical(tab$Status, s$checks$Status)
  html <- mfrm_report(x, output = "html")
  expect_true(any(grepl("Facet|Person", readLines(html$path, warn = FALSE))))
  saveRDS(x, file.path(out, "recommendation.rds"))
  expect_identical(summary(readRDS(file.path(out, "recommendation.rds"))), s)
  skip_if_not_installed("ggplot2")
  g <- as_ggplot(x)
  expect_equal(nrow(ggplot2::ggplot_build(g)$data[[1]]), nrow(s$checks))
  expect_identical(plot_data(g, "table"), plot_data(p, "table"))
  expect_null(g$labels$title)
})

test_that("main-effects D reports preserve source, assumptions and unavailable projections", {
  skip_if_not_installed("lme4")
  f <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))$fit
  g <- mfrm_generalizability(f)
  d <- mfrm_d_study(g, data.frame(Rater = c(2, 4), Criterion = 3), residual_scaling = "sensitivity")
  r <- mfrm_report(d)
  expect_equal(nrow(r$tables$coefficients), 6L)
  expect_identical(r$tables$coefficients, summary(d))
  expect_identical(r$tables$components, as.data.frame(g$variance_components))
  expect_identical(r$tables$source_counts$ObservedLevels, unname(g$design$observed_levels))
  expect_identical(r$tables$uncertainty$Status, "not_available")
  expect_true(all(c("highest_order", "single_condition", "none") %in% r$tables$coefficients$ResidualScaling))
  expect_s3_class(apa_table(d, which = "coefficients"), "apa_table")
  expect_equal(build_summary_table_bundle(d)$tables$coefficients$G, d$G)
  expect_error(mfrm_report(d, intervals = TRUE), "only to G/D")
  old <- d; attr(old, "source_design") <- attr(old, "source_components") <- NULL
  ro <- mfrm_report(old)
  expect_identical(ro$tables$source_design$SourceDetails, "not_recorded")
  expect_false("components" %in% names(ro$tables))
  expect_identical(ro$tables$coefficients$Phi, d$Phi)
  expect_identical(attr(d[1, ], "source_components"), g$variance_components)
  # Stored missing components remain missing; the report must not repair them.
  bad <- g; bad$variance_components$Variance[1] <- NA_real_
  bd <- mfrm_d_study(bad, data.frame(Rater = 2, Criterion = 3))
  br <- mfrm_report(bd)
  expect_identical(br$tables$coefficients$GReason, bd$GReason)
  expect_true(anyNA(br$tables$components$Variance))
  out <- tempfile(); dir.create(out)
  export_summary_appendix(r, output_dir = out, prefix = "approximation")
  csv <- read.csv(list.files(out, "coefficients.*csv$", full.names = TRUE))
  expect_equal(csv$G, d$G, tolerance = 1e-12)
  saveRDS(r, file.path(out, "report.rds"))
  expect_identical(readRDS(file.path(out, "report.rds"))$tables, r$tables)
})

test_that("numerical next actions separate evidence without changing readiness", {
  saved <- readRDS(system.file("examples", "operational-mml.rds", package = "mfrmr"))
  f <- saved$fit
  expected <- .mfrm_numerical_review(f)
  expect_identical(summary(f)$numerical_review, expected)
  res <- mfrm_results(f, diagnostics = saved$diagnostics)
  expect_identical(mfrm_report(res)$tables$numerical_review, expected)
  changed <- f
  changed$opt$optimizer_diagnostics$ConvergenceSeverity <- "review"
  changed$opt$optimizer_diagnostics$ConvergenceStatus <- "iteration_limit"
  q <- .mfrm_numerical_review(changed)
  expect_match(q$NextAction[q$Area == "Optimization"], "maxit ceiling")
  changed$opt$mml_engine$Used <- "em"
  changed$opt$optimizer_diagnostics$ConvergenceStatus <- "gradient_review"
  q <- .mfrm_numerical_review(changed)
  expect_match(q$NextAction[q$Area == "Optimization"], "direct or hybrid")
  changed$opt$start_review$table <- data.frame(NumericalStatus = c("pass", "error"))
  expect_identical(.mfrm_numerical_review(changed)$Status[3], "review")
  changed$opt$optimizer_polish <- NULL
  expect_identical(.mfrm_numerical_review(changed)$Status[2], "not_assessed")
  expect_match(q$NextAction[q$Area == "Integration"], "does not promote readiness")
  expect_identical(f, saved$fit)
})
