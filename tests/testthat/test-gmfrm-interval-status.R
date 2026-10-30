# These output-contract tests use saved outcomes and explicit test doubles, not
# a new source qualification. Fitting, information and profile search are blocked.
test_that("both two-family owners retain their fixed-unit-SD contrasts", {
  x <- readRDS(test_path("fixtures", "gmfrm-explicit-start.rds"))
  base <- x$intervals$locations
  local_mocked_bindings(mfrm_gpcm_product_inference = function(...) list(
    covariance = base$parameter_covariance, checks = base$checks,
    numerical_checks = base$numerical_checks, score_rank = base$score_rank,
    check = list(eligible = TRUE, review = "Retained qualification")),
    compute_mml_parameter_covariance = function(...) stop("No new information"), .package = "mfrmr")
  for (facet in x$fit$config$slope_facet) {
    C <- matrix(c(1, -1), 1, dimnames = list("Difference", x$fit$config$facet_levels[[facet]]))
    native <- suppressWarnings(mfrm_facet_intervals(x$fit, facet, C))
    std <- suppressWarnings(mfrm_facet_intervals(x$fit, facet, C, scale = "standardized"))
    expect_equal(std$table$Estimate, native$table$Estimate)
    expect_equal(std$table$SE, native$table$SE)
    expect_equal(std$table$Lower, native$table$Lower)
    expect_false(std$settings$population_sd_estimated)
    expect_true(std$table$CIEligible)
    expect_false(std$table$CoverageQualified)
    p <- plot_data(plot(std, draw = FALSE))
    expect_identical(p$xlab, "Contrast (population SD units)")
    expect_false(grepl("Uncentered locations|Locations sum to zero", p$subtitle))
  }
})

gmfrm_status_outputs <- function(rejected = FALSE, level = .9) {
  x <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))
  fit <- x$fit; profile <- x$profiles$Task
  V <- solve(x$baseline$information)
  inference <- list(covariance = if (rejected) NULL else V,
    checks = attr(profile, "checks"), numerical_checks = attr(profile, "numerical_checks"),
    score_rank = attr(profile, "score_rank"),
    check = list(eligible = !rejected, review = if (rejected) "Retained numerical refusal" else
      "Retained covariance for output regression only; coverage is unqualified."))
  forbidden <- function(...) stop("Output regression must not repeat numerical estimation")
  local_mocked_bindings(mfrm_gpcm_product_inference = function(...) inference,
    fit_mfrm = forbidden, compute_mml_parameter_covariance = forbidden,
    mfrm_gpcm_profile_search = forbidden, .package = "mfrmr")
  lev <- fit$config$facet_levels$Rater
  C <- rbind(difference = c(1, -1, rep(0, length(lev) - 2)),
    overflow = c(1e200, -1e200, rep(0, length(lev) - 2)))
  colnames(C) <- lev
  slopes <- suppressWarnings(confint(fit, parm = "slopes", level = level))
  locations <- suppressWarnings(mfrm_facet_intervals(fit, "Rater", C, level = .9))
  ordinary <- suppressWarnings(mfrm_facet_intervals(fit, "Rater", C[1, , drop = FALSE], level = .9))
  fixed <- matrix(1, 1, length(fit$config$facet_levels$Task),
    dimnames = list("sum", fit$config$facet_levels$Task))
  fixed <- suppressWarnings(mfrm_facet_intervals(fit, "Task", fixed))
  # Reuse the saved profile outcomes; do not search or re-evaluate its likelihood.
  if (!rejected) {
    local_mocked_bindings(make_mfrm_direct_evaluator = function(...) list(
      value = function(p) fit$opt$value, gradient = function(p) numeric(length(p))),
      mfrm_gpcm_profile_search = function(...) attr(profile, "profile"), .package = "mfrmr")
    profile <- suppressWarnings(confint(fit, parm = "slopes", method = "profile", slope = c(Task = "t3")))
  }
  list(fit = fit, slopes = slopes, locations = locations, fixed = fixed,
    ordinary = ordinary, profile = profile, covariance = V)
}

test_that("numerical availability is distinct from coverage in two-family public outputs", {
  x <- gmfrm_status_outputs()
  slopes <- attr(x$slopes, "diagnostics")
  expect_true(all(slopes$CIEligible))
  expect_identical(slopes$CoverageQualified, rep(FALSE, nrow(slopes)))
  expect_identical(x$locations$table$CIEligible, c(TRUE, FALSE))
  expect_identical(x$locations$table$CoverageQualified, c(FALSE, FALSE))
  expect_identical(x$fixed$table$Status, "fixed")
  expect_identical(x$fixed$table$CIEligible, FALSE)
  expect_identical(x$fixed$table$CoverageQualified, FALSE)
  expect_identical(attr(x$profile, "diagnostics")$CoverageQualified, FALSE)
  expect_false(attr(x$profile, "diagnostics")$CIEligible)
  # Independent coordinate map for this retained three-task/six-rater model.
  J <- matrix(0, 9, 22)
  J[1:3, 15:16] <- rbind(c(1, 0), c(0, 1), c(-1, -1))
  J[4:9, 17:22] <- diag(6)
  expected <- drop(J %*% x$fit$opt$par)
  se <- sqrt(diag(J %*% x$covariance %*% t(J)))
  expect_equal(slopes$CI_Lower, exp(expected - qnorm(.95) * se), tolerance = 1e-12)
  expect_equal(slopes$CI_Upper, exp(expected + qnorm(.95) * se), tolerance = 1e-12)
  expect_equal(x$locations$table$SE[1]^2,
    x$covariance[3,3] + x$covariance[4,4] - 2 * x$covariance[3,4], tolerance = 1e-12)
  y <- gmfrm_status_outputs(TRUE)
  expect_false(any(attr(y$slopes, "diagnostics")$CIEligible))
  expect_false(any(attr(y$slopes, "diagnostics")$CoverageQualified))
  expect_true(all(is.na(y$locations$table$Lower)))
  expect_true(all(is.finite(y$locations$table$Estimate)))
  expect_match(y$locations$table$InferenceReview[1], "numerical refusal")
})

test_that("two-family log-Wald intervals retain finite extreme-level bounds and qualification", {
  level <- 1 - .Machine$double.eps / 2
  x <- gmfrm_status_outputs(level = level)
  tab <- attr(x$slopes, "diagnostics")
  expect_true(all(tab$CIEligible))
  expect_false(any(tab$CoverageQualified))
  z <- sqrt(qchisq(level, 1))
  expect_equal(log(tab$CI_Lower), log(tab$Estimate)-z*tab$LogSE)
  expect_equal(log(tab$CI_Upper), log(tab$Estimate)+z*tab$LogSE)
})

test_that("nonrepresentable two-family contrast uncertainty cannot leak into model columns", {
  x <- gmfrm_status_outputs(); tab <- x$locations$table
  expect_true(all(is.finite(tab$Estimate)))
  expect_identical(tab$Status[2], "nonpositive_variance")
  fields <- c("SE", "Lower", "Upper", "ModelSE", "ModelLower", "ModelUpper")
  expect_true(all(is.na(tab[2, fields])))
  expect_true(all(is.finite(unlist(tab[1, fields]))))
  expect_match(tab$InferenceReview[2], "finite positive variance")
  expect_match(x$fixed$table$InferenceReview, "Constraint-fixed")
  # Every coefficient is finite, but the resulting point is not representable.
  location <- subset(x$fit$facets$others, Facet == "Rater")
  C <- matrix(.Machine$double.xmax * sign(location$Estimate), 1,
    dimnames = list("overflow", location$Level))
  expect_gt(sum(abs(location$Estimate)), 1)
  expect_error(mfrm_facet_intervals(x$fit, "Rater", C), "estimates must be finite")
})

test_that("computed two-family profiles retain their endpoints without acquiring coverage support", {
  x <- readRDS(test_path("fixtures", "gmfrm-profile-intervals.rds"))
  fit <- x$fit; saved <- x$profiles$Task
  forbidden <- function(...) stop("Output regression must not repeat numerical estimation")
  # The comparison-Wald covariance is a test double; only the actual retained
  # profile endpoints and their availability are asserted in this display test.
  local_mocked_bindings(mfrm_gpcm_product_inference = function(...) list(
    covariance = diag(length(fit$opt$par)), check = list(eligible = TRUE)),
    make_mfrm_direct_evaluator = function(...) list(value = function(p) fit$opt$value,
      gradient = function(p) numeric(length(p))),
    mfrm_gmfrm_profile_evaluator = function(problem, direct) direct,
    mfrm_gpcm_profile_search = function(...) attr(saved, "profile"),
    fit_mfrm = forbidden, compute_mml_parameter_covariance = forbidden, .package = "mfrmr")
  ci <- confint(fit, parm = "slopes", method = "profile", slope = c(Task = "t3"))
  tab <- attr(ci, "diagnostics")
  expect_true(tab$CIEligible)
  expect_identical(tab$CoverageQualified, FALSE)
  expect_identical(tab$CI_Lower, attr(saved, "diagnostics")$CI_Lower)
  expect_identical(tab$CI_Upper, attr(saved, "diagnostics")$CI_Upper)
})

test_that("coverage metadata survives tables reports CSV and saved replay without estimation", {
  x <- gmfrm_status_outputs()
  forbidden <- function(...) stop("Saved output must not repeat numerical estimation")
  local_mocked_bindings(fit_mfrm = forbidden, compute_mml_parameter_covariance = forbidden,
    mfrm_gpcm_product_inference = forbidden, mfrm_gpcm_profile_search = forbidden, .package = "mfrmr")
  # The ordinary row is used for figures; the extreme row remains covered above.
  ci <- x$ordinary
  res <- mfrm_results(x$fit, include = c("fit", "plots"), compute = "never",
    intervals = list(slopes = x$slopes, locations = ci, profile = x$profile))
  report <- mfrm_report(res)
  for (name in c("slopes", "locations", "profile")) {
    key <- paste0("gpcm_", name, "_intervals")
    tab <- res$tables[[key]]
    expect_identical(tab$CoverageQualified, rep(FALSE, nrow(tab)))
    expect_identical(report$tables[[key]]$CoverageQualified, tab$CoverageQualified)
  }
  expect_true("CoverageQualified" %in% names(apa_table(x$slopes)$table))
  expect_true("CoverageQualified" %in% names(apa_table(ci)$table))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    for (z in list(x$slopes, ci)) {
      figure <- as_ggplot(z)
      expect_false(any(plot_data(figure, "table")$CoverageQualified))
      expect_null(figure$labels$title)
      expect_null(figure$labels$caption)
      expect_no_warning(ggplot2::ggplotGrob(figure))
    }
  }
  folder <- withr::local_tempdir()
  exported <- export_mfrm_results(res, output_dir = folder,
    include = c("tables", "report", "replay"), acknowledge_sensitive = TRUE)
  for (name in c("slopes", "locations", "profile")) {
    key <- paste0("table_gpcm_", name, "_intervals")
    csv <- exported$written_files$Path[exported$written_files$Component == key]
    expect_length(csv, 1L)
    tab <- read.csv(csv, check.names = FALSE)
    expect_identical(tab$CoverageQualified, rep(FALSE, nrow(tab)))
  }
  replay <- list.files(folder, pattern = "_replay[.]R$", full.names = FALSE)
  e <- new.env(parent = globalenv())
  withr::with_dir(folder, sys.source(replay, envir = e))
  expect_identical(e$res, res)
  expect_identical(e$report$markdown, report$markdown)
  # Existing saved results remain readable; absence is not filled with TRUE.
  old <- x$slopes; attr(old, "diagnostics")$CoverageQualified <- NULL
  expect_false("CoverageQualified" %in% names(mfrm_gpcm_inference_tables(old)$intervals))
})
