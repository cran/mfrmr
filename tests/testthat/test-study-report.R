study_report_fixture <- function() {
  data("mfrmr_example_gtheory_crossed", package = "mfrmr", envir = environment())
  dat <- mfrmr_example_gtheory_crossed
  g <- mfrm_g_study(dat, c("Content", "Organization"))
  d <- mfrm_d_study(g, expand.grid(Raters = c(2, 4), Tasks = c(3, 6)),
    weights = c(Content = .5, Organization = .5))
  list(g = g, d = d)
}

test_that("G/D report handoff retains estimands, errors and full table values", {
  z <- study_report_fixture()
  r <- mfrm_report(z$d)
  expect_s3_class(r, "mfrm_report")
  expect_identical(r$tables$coefficients, z$d$coefficients)
  expect_identical(r$tables$uncertainty$Status, "not_requested")
  expect_identical(r$tables$weights$Estimate, c(.5, .5))
  expect_equal(r$tables$components$Estimate, unlist(z$g$components, use.names = FALSE))
  expect_equal(r$tables$source_counts$ObservedLevels, c(40, 8, 6))
  expect_identical(r$tables$design_grid$Tasks, z$d$design_grid$Tasks)
  for (input in list(z$d, r, summary(r, top_n = 1))) {
    b <- build_summary_table_bundle(input, top_n = 1)
    expect_identical(b$tables$coefficients, z$d$coefficients)
    expect_true(all(b$table_index$Role %in%
      mfrmr:::summary_table_bundle_appendix_role_registry()$Role))
  }
  expect_s3_class(apa_table(r, which = "coefficients"), "apa_table")
  expect_s3_class(apa_table(z$d, which = "coefficients"), "apa_table")
  expect_true(grepl("error", mfrm_report(z$g, output = "markdown")))
  html <- mfrm_report(z$d, output = "html")
  expect_true(file.exists(html$path))
  expect_true(any(grepl("RelativeSEM", readLines(html$path, warn = FALSE))))
  out <- tempfile(); dir.create(out)
  export_summary_appendix(r, output_dir = out, prefix = "study")
  csv <- list.files(out, "coefficients.*csv$", full.names = TRUE)
  expect_length(csv, 1)
  tab <- read.csv(csv)
  expect_equal(tab$G, z$d$coefficients$G, tolerance = 1e-12)
  expect_identical(tab$Status, z$d$coefficients$Status)
  saveRDS(r, file.path(out, "report.rds"))
  expect_identical(readRDS(file.path(out, "report.rds"))$tables, r$tables)
  expect_error(mfrm_report(z$d, style = "rater"), "rater")
  expect_error(mfrm_report(z$d, facet = "Rater"), "Rater-specific")
  expect_error(mfrm_report(z$d, intervals = TRUE), "saved")
})

test_that("saved intervals are source checked and unavailable rows survive", {
  z <- study_report_fixture()
  ci <- confint(z$d, parm = "G", assumption = "normal", method = "delta")
  r <- mfrm_report(z$d, intervals = ci)
  expect_equal(r$tables$intervals$Estimate, ci$Estimate)
  expect_equal(r$tables$intervals$Lower, ci$Lower)
  expect_identical(r$intervals[[1]], ci)
  other <- z$d; other$gstudy$components$Person[1, 1] <- other$gstudy$components$Person[1, 1] + .1
  expect_error(mfrm_report(other, intervals = ci), "source summaries")
  gc <- confint(z$g, score = "Content", assumption = "normal")
  gr <- mfrm_report(z$g, intervals = gc)
  expect_identical(gr$intervals[[1]], gc)
  expect_identical(gr$tables$interval_targets$Weight, c(1, 0))
  altered <- gc; attr(altered, "weights") <- c(Content = 0, Organization = 1)
  expect_error(mfrm_report(z$g, intervals = altered), "point estimates")
  covariance <- confint(z$g, score = c("Content", "Organization"),
    method = "regression_region", assumption = "normal")
  cr <- mfrm_report(z$g, intervals = covariance)
  expect_identical(cr$tables$interval_targets$Score, c("Content", "Organization"))
  expect_equal(cr$tables$intervals$Estimate, covariance$Estimate)
  otherg <- z$g; otherg$mean_products$Person[1, 1] <- otherg$mean_products$Person[1, 1] + 1
  expect_error(mfrm_report(otherg, intervals = gc), "source does not match")
  residual_ci <- confint(z$g, score = "Content", method = "residual_contrast", assumption = "normal")
  otherg <- z$g; otherg$data$Content[1] <- otherg$data$Content[1] + 1
  expect_error(mfrm_report(otherg, intervals = residual_ci), "source does not match")
  old <- gc; attr(old, "source_summary") <- NULL
  expect_error(mfrm_report(z$g, intervals = old), "recreate")
  # Constant scores create non-estimable coefficient/interval rows, not zeros.
  data("mfrmr_example_gtheory_crossed", package = "mfrmr", envir = environment())
  dat <- mfrmr_example_gtheory_crossed; dat$Content <- 1
  flat <- mfrm_d_study(mfrm_g_study(dat, "Content"), data.frame(Raters = 2, Tasks = 3))
  unavailable <- confint(flat, assumption = "normal", method = "delta")
  fr <- mfrm_report(flat, intervals = unavailable)
  expect_identical(fr$tables$intervals$Available, unavailable$Available)
  expect_true(any(!fr$tables$intervals$Available))
  expect_equal(fr$tables$coefficients$G, flat$coefficients$G)
})

test_that("signed components, incomplete source and paired comparisons retain their meaning", {
  z <- study_report_fixture()
  # Report a saved signed component without silently repairing it.
  z$g$components$Rater[1, 2] <- z$g$components$Rater[2, 1] <- -.25
  r <- mfrm_report(z$g)
  expect_equal(r$tables$components$Estimate[r$tables$components$Source == "Rater"],
    as.vector(z$g$components$Rater))
  dat <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
  dat <- dat[-1, ]
  g <- mfrm_g_study(dat, c("V", "W"), rater = NULL, method = "minque0")
  expect_false(mfrm_report(g)$tables$source_design$Complete)
  cmp <- mfrm_multivariate_d_compare(z$d, assumption = "normal")
  r <- mfrm_report(cmp)
  expect_identical(r$tables$comparisons, cmp$comparisons)
  expect_equal(r$tables$comparison_covariance$Estimate, as.vector(cmp$covariance))
})
