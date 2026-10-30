gmfrm_feedback_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      fit <- gmfrm_scoring_fixture()$fit
      data <- fit$gmfrm$specification$data
      rows <- which(data$Person %in% unique(data$Person)[c(1, 20, 160)])
      diagnostics <- mfrm_response_diagnostics(fit, rows = rows, quad_points = 121L)
      intervals <- suppressWarnings(confint(fit, parm = "slopes"))
      cached <<- mfrm_results(fit, response_diagnostics = diagnostics, intervals = list(slopes = intervals))
    }
    cached
  }
})

test_that("two-family recipient sheets preserve component targets and saved residuals", {
  x <- gmfrm_feedback_fixture(); before <- x
  local_mocked_bindings(fit_mfrm = function(...) stop("no fitting"),
    confint.mfrm_fit = function(...) stop("no intervals"),
    mfrm_response_diagnostics = function(...) stop("no diagnostics"),
    mfrm_gpcm_response_evaluator = function(...) stop("no probabilities"), .package = "mfrmr")
  for (facet in c("Task", "Rater")) {
    id <- x$fit$config$facet_levels[[facet]][1]
    sheet <- mfrm_report(x, style = "rater", facet = facet, rater = id)
    expected <- x$fit$facets$others
    expect_equal(sheet$tables$severity$Location, unname(expected$Estimate[expected$Facet == facet & expected$Level == id]))
    expected <- x$fit$slopes
    expect_equal(sheet$tables$slope$ComponentSlope, expected$OptimizerEstimate[expected$SlopeOwner == facet & expected$SlopeFacet == id])
    tab <- attr(x$gpcm_inference$slopes, "diagnostics")
    at <- tab$SlopeOwner == facet & tab$SlopeLevel == id
    expect_equal(sheet$tables$uncertainty$Lower, tab$CI_Lower[at])
    expect_equal(sheet$tables$uncertainty$Upper, tab$CI_Upper[at])
    expect_identical(sheet$tables$uncertainty$Available, tab$CIEligible[at])
    expect_equal(nrow(sheet$tables$steps), if (facet == "Rater") 2L else 0L)
    expect_match(sheet$notes["slope"], if (facet == "Task") "geometric mean one" else "need not equal one")
    d <- x$response_diagnostics
    rows <- d$rows[as.character(d$source_data[[facet]][d$rows$InputRow]) == id, ]
    take <- x$fit$prep$data[[facet]] == id
    expect_equal(sheet$tables$exposure$RatingRows, sum(take))
    expect_equal(sheet$tables$residuals$NotIncluded, sum(take) - nrow(rows))
    expect_equal(sheet$tables$residuals$MeanResidual, mean(rows$Score - rows$ExpectedScore))
    expect_equal(sheet$tables$residuals$MeanSquaredStandardizedResidual, mean(rows$StandardizedResidual^2))
    ord <- head(order(-abs(rows$StandardizedResidual), rows$InputRow), 5L)
    expect_equal(sheet$tables$cases$Expected, rows$ExpectedScore[ord])
    expect_equal(sheet$tables$cases$StandardizedResidual, rows$StandardizedResidual[ord])
    expect_equal(sum(sheet$tables$categories$Ratings), sum(take))
    expect_equal(sum(sheet$tables$categories$Percent), 100)
    expect_match(sheet$notes["residuals"], "no expectation-one")
    expect_identical(mfrm_report(x, "rater", "tables", facet, id), sheet$tables)
    expect_identical(mfrm_report(x, "rater", "markdown", facet, id), sheet$markdown)
    numeric_tables <- setdiff(names(sheet$tables), "interpretation_notes")
    expect_identical(mfrm_report(x, "rater", facet = facet, rater = id, audience = "researcher")$tables[numeric_tables],
      sheet$tables[numeric_tables])
    expect_output(print(sheet), "Two-family GPCM")
  }
  expect_identical(x, before)
})

test_that("recipient sheets retain missing sections, failed rows and interval choices", {
  x <- gmfrm_feedback_fixture()
  sheet <- function(x, ...) mfrm_report(x, "rater", facet = "Rater", rater = "r1", ...)
  missing <- x; missing$response_diagnostics <- missing$gpcm_inference <- NULL
  s <- sheet(missing)
  expect_equal(nrow(s$tables$uncertainty), 0L)
  expect_equal(nrow(s$tables$residuals), 0L)
  expect_match(s$notes["uncertainty"], "No saved")
  expect_match(s$notes["residuals"], "No saved")
  missing$fit$summary$Converged <- FALSE
  expect_match(sheet(missing)$review, "convergence is unresolved")
  bad <- x; bad$fit$slopes$OptimizerEstimate[1] <- 999
  expect_error(sheet(bad), "must match")
  bad <- x; bad$response_diagnostics$source_data$Score[1] <- -1
  expect_error(sheet(bad), "matching calibration")
  bad <- x
  at <- which(bad$response_diagnostics$source_data$Rater[bad$response_diagnostics$rows$InputRow] == "r1")[1]
  bad$response_diagnostics$rows$Status[at] <- "unavailable"
  bad$response_diagnostics$rows$ExpectedScore[at] <- NA_real_
  s <- sheet(bad)
  expect_equal(s$tables$residuals$Unavailable, 1L)
  expect_true(is.na(s$tables$residuals$MeanResidual))
  expect_true(is.na(s$tables$residuals$MeanSquaredStandardizedResidual))
  expect_match(s$review, "unresolved rows")
  expect_equal(s$tables$residuals$Available + s$tables$residuals$Unavailable, s$tables$residuals$SavedRows)
  expect_equal(nrow(sheet(x, max_cases = 0)$tables$cases), 0L)
  expect_error(sheet(x, interval = "missing"), "individual component slope")
  bad <- x; bad$gpcm_inference$second <- bad$gpcm_inference$slopes
  expect_error(sheet(bad), "Several saved intervals")
  expect_identical(sheet(bad, interval = "second")$tables, sheet(x)$tables)
  bad <- x; attr(bad$gpcm_inference$slopes, "source")$parameters[1] <- 999
  expect_error(sheet(bad), "must match")
  bad <- x; attr(bad$gpcm_inference$slopes, "diagnostics")$CIEligible[] <- FALSE
  expect_false(sheet(bad)$tables$uncertainty$Available)
  expect_match(sheet(bad)$notes["uncertainty"], "unavailable")
  bad <- x; attr(bad$gpcm_inference$slopes, "level") <- .8
  expect_error(sheet(bad), "estimate, scale and method")
  bad <- x; attr(bad$gpcm_inference$slopes, "simultaneous") <- "bonferroni"
  expect_identical(sheet(bad)$tables$uncertainty$Adjustment, "Bonferroni")
  expect_error(mfrm_report(x, "rater", facet = "Person", rater = "p1"), "non-Person")
})

test_that("arbitrary owner names and private fields stay out of distributed sheets", {
  x <- gmfrm_feedback_fixture(); fit <- x$fit
  data <- fit$gmfrm$specification$data
  data$Person <- paste0("PRIVATE_PERSON_", data$Person)
  data$Task <- sub("t", "PRIVATE_LEVEL_", as.character(data$Task), fixed = TRUE)
  data$Rater <- sub("r", "PRIVATE_LEVEL_", as.character(data$Rater), fixed = TRUE)
  names(data) <- c("Candidate", "観点 名", "Judge-ID", "Rating")
  p <- mfrm_gmfrm_problem(data, 2L, gauss_hermite_normal(61L),
    slope_facets = names(data)[2:3], person = "Candidate", score = "Rating")
  renamed <- mfrm_gmfrm_fit_result(p, mfrm_gmfrm_em(p, start = fit$opt$par, maxit = 1L))
  rows <- which(data$Candidate %in% unique(data$Candidate)[1:2])
  d <- mfrm_response_diagnostics(renamed, rows = rows, quad_points = 121L)
  res <- mfrm_results(renamed, response_diagnostics = d)
  res$notes <- "PRIVATE_SOURCE_NOTE"
  res$response_diagnostics$rows$Reason <- "PRIVATE_REASON"
  attr(res$response_diagnostics$rows$ExpectedScore, "private") <- "PRIVATE_ATTRIBUTE"
  for (facet in names(data)[2:3]) {
    s <- mfrm_report(res, "rater", facet = facet, rater = "PRIVATE_LEVEL_1")
    expect_false(grepl("PRIVATE_|観点 名|Judge-ID|Candidate", paste(capture.output(dput(s)), collapse = "\n")))
    expect_equal(s$tables$categories$Ratings, as.numeric(table(factor(
      data$Rating[data[[facet]] == "PRIVATE_LEVEL_1"], levels = 0:2))))
  }
  h <- mfrm_report(res, "rater", "html", "Judge-ID", "PRIVATE_LEVEL_1", label = "<script>example</script>")
  on.exit(unlink(h$path), add = TRUE)
  expect_false(grepl("PRIVATE_", paste(capture.output(dput(h)), collapse = "\n")))
  expect_false(grepl("<script>", h$html, fixed = TRUE))
  expect_match(h$html, "&lt;script&gt;", fixed = TRUE)
  expect_match(h$html, 'scope="col"', fixed = TRUE)
  expect_match(h$html, 'lang="en"', fixed = TRUE)
  expect_match(h$html, "@media print", fixed = TRUE)
  expect_false(grepl("src=|href=", h$html))
  expect_identical(paste(readLines(h$path, warn = FALSE), collapse = "\n"), h$html)
})

test_that("two-family sheets explain failed checks and withhold partial intervals", {
  x <- gmfrm_feedback_fixture()
  sheet <- function(x) mfrm_report(x, "rater", facet = "Rater", rater = "r1")
  # Saved-state presentation tests; no claim of a newly fitted boundary case.
  ci <- x$gpcm_inference$slopes
  attr(ci, "diagnostics")$CIEligible[] <- FALSE
  attr(ci, "diagnostics")$InferenceReview[] <- "PRIVATE_REVIEW"
  attr(ci, "checks") <- data.frame(Check = "Joint information", Passed = FALSE, Detail = "PRIVATE_DETAIL")
  x$gpcm_inference$slopes <- ci
  s <- sheet(x)
  expect_match(s$tables$uncertainty$Reason, "joint information")
  expect_true(all(is.na(s$tables$uncertainty[c("Lower", "Upper")])))
  expect_false(grepl("PRIVATE_", paste(capture.output(dput(s)), collapse = "\n")))
  ci <- x$gpcm_inference$slopes
  attr(ci, "settings")$method <- "profile"
  attr(ci, "profile") <- list(endpoints = data.frame(Side = c("PRIVATE_LOWER", "PRIVATE_UPPER"),
    Status = c("computed", "crossing_not_found"), Detail = "PRIVATE_DETAIL"))
  attr(ci, "diagnostics")$CI_Upper[] <- NA_real_
  x$gpcm_inference$slopes <- ci
  original <- serialize(x, NULL)
  s <- sheet(x)
  expect_false(s$tables$uncertainty$Available)
  expect_match(s$tables$uncertainty$Reason, "cutoff crossing")
  expect_match(s$tables$uncertainty$Reason, "does not establish an unbounded interval")
  expect_true(all(is.na(s$tables$uncertainty[c("Lower", "Upper")])))
  expect_match(s$markdown, "cutoff crossing")
  expect_output(print(s), "cutoff crossing")
  expect_false(grepl("PRIVATE_", paste(capture.output(dput(s)), collapse = "\n")))
  expect_identical(serialize(x, NULL), original)
  h <- mfrm_report(x, "rater", "html", "Rater", "r1")
  on.exit(unlink(h$path))
  expect_match(h$html, "cutoff crossing")
  expect_false(grepl("PRIVATE_", h$html))
})

test_that("saved two-family recipient sheets reopen without fitting or integration", {
  skip_if_not_installed("callr")
  x <- gmfrm_feedback_fixture()
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE); saveRDS(x, file)
  worker <- function(root, file) {
    if (file.exists(file.path(root, "R", "api-results.R"))) pkgload::load_all(root, quiet = TRUE, compile = FALSE) else
      library("mfrmr", lib.loc = dirname(root), character.only = TRUE)
    testthat::local_mocked_bindings(fit_mfrm = function(...) stop("no fitting"),
      mfrm_gmfrm_em = function(...) stop("no EM"),
      confint.mfrm_fit = function(...) stop("no intervals"),
      mfrm_gpcm_response_evaluator = function(...) stop("no probabilities"),
      mfrm_response_diagnostics = function(...) stop("no diagnostics"), .package = "mfrmr")
    res <- readRDS(file)
    mfrmr::mfrm_report(res, "rater", facet = "Rater", rater = "r1")
  }
  environment(worker) <- baseenv()
  reopened <- callr::r(worker, list(normalizePath(find.package("mfrmr")), file), libpath = .libPaths())
  expect_identical(reopened, mfrm_report(x, "rater", facet = "Rater", rater = "r1"))
})
