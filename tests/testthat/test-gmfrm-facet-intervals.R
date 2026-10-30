gmfrm_location_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      saved <- readRDS(test_path("fixtures", "gmfrm-joint-information.rds"))
      problem <- mfrm_gmfrm_problem(saved$data, 2L, gauss_hermite_normal(31L))
      fit <- mfrm_gmfrm_fit_result(problem,
        mfrm_gmfrm_em(problem, start = saved$parameters, maxit = 1L))
      cached <<- list(fit = fit, saved = saved,
        inference = mfrm_gpcm_product_inference(fit))
    }
    cached
  }
})

test_that("two-family locations use constrained blocks of the inverse full information", {
  x <- gmfrm_location_fixture(); fit <- x$fit; before <- fit
  # These targets share one freshly checked source, not separate calibrations.
  local_mocked_bindings(mfrm_gpcm_product_inference = function(source) {
    expect_identical(source, fit)
    x$inference
  }, .package = "mfrmr")
  V <- solve(x$saved$information)
  for (facet in c("Task", "Rater")) {
    ci <- suppressWarnings(mfrm_facet_intervals(fit, facet))
    J <- matrix(0, 3, 13)
    if (facet == "Task") J[, 1:2] <- rbind(c(1, 0), c(0, 1), c(-1, -1)) else J[, 3:5] <- diag(3)
    expected <- J %*% V %*% t(J)
    expect_equal(ci$parameter_covariance, V, ignore_attr = TRUE, tolerance = 1e-8)
    expect_equal(ci$covariance, expected, ignore_attr = TRUE, tolerance = 1e-8)
    expect_equal(ci$table$Estimate, drop(J %*% fit$opt$par))
    expect_equal(ci$table$SE, sqrt(diag(expected)), tolerance = 1e-8)
    expect_equal(ci$table$Lower, ci$table$Estimate - qnorm(.975) * ci$table$SE)
    expect_true(all(ci$table$CIEligible & ci$table$Status == "available"))
    expect_true(all(ci$checks$Passed))
    expect_match(ci$cautions[length(ci$cautions)], "Location differences do not imply uniform")
    expect_identical(ci$settings$location_reference, if (facet == "Task")
      "sum_to_zero_across_levels" else "uncentered_on_fixed_N01")
    if (facet == "Task") expect_equal(unname(colSums(ci$covariance)), rep(0, 3), tolerance = 1e-12)
    C <- rbind(Difference = c(1, -1, 0), Reverse = c(-1, 1, 0), Sum = c(1, 1, 1))
    colnames(C) <- ci$table$Target
    pair <- suppressWarnings(mfrm_facet_intervals(fit, facet, C[, 3:1], level = .9))
    expect_equal(pair$table$Estimate, as.numeric(C %*% ci$table$Estimate))
    expect_equal(pair$covariance, C %*% expected %*% t(C), ignore_attr = TRUE, tolerance = 1e-8)
    expect_equal(pair$table$Lower[1], -pair$table$Upper[2])
    expect_equal(pair$table$SE[1]^2, sum(expected[1:2, 1:2] * outer(c(1, -1), c(1, -1))), tolerance = 1e-8)
    expect_gt(abs(pair$table$SE[1]^2 - sum(diag(expected)[1:2])), 1e-5)
    if (facet == "Task") {
      expect_identical(pair$table$Status[3], "fixed")
      expect_identical(pair$table$SE[3], 0)
      expect_true(is.na(pair$table$Lower[3]))
    }
  }
  expect_identical(fit, before)
  # Inverting a location-only Hessian would incorrectly treat all nuisance terms as fixed.
  expect_gt(max(abs(V[1:5, 1:5] - solve(x$saved$information[1:5, 1:5]))), 1e-4)
})

test_that("two-family location failures keep points and explain missing intervals", {
  fit <- gmfrm_location_fixture()$fit
  failed <- fit; failed$summary$Converged <- FALSE
  ci <- suppressWarnings(mfrm_facet_intervals(failed, "Rater"))
  expect_true(all(ci$table$Status == "numerical_review_failed"))
  expect_true(all(is.na(ci$table$SE) & is.na(ci$table$ModelSE) & is.na(ci$table$Lower)))
  expect_true(all(is.finite(ci$table$Estimate)))
  expect_match(ci$table$InferenceReview[1], "EM must meet")
  expect_false(all(ci$checks$Passed))
  expect_error(mfrm_facet_intervals(fit, "Rater", method = "sandwich"), "sandwich intervals are unavailable")
  expect_error(mfrm_facet_intervals(fit, "Rater", adjust = TRUE), "only to method")
  expect_error(mfrm_facet_intervals(fit, "Person"), "non-person")
  for (level in list(0, 1, NA_real_, matrix(.9), .9 + 1i))
    expect_error(mfrm_facet_intervals(fit, "Rater", level = level), "level")
  bad <- fit; bad$facets$others$Estimate[1] <- 999
  expect_error(mfrm_facet_intervals(bad, "Rater"), "must match")
  # Refinement failure must not leak the original, otherwise invertible covariance.
  rejected <- gmfrm_location_fixture()$inference
  rejected$covariance <- NULL; rejected$check$eligible <- FALSE
  rejected$check$review <- "Quadrature comparison failed"
  rejected$checks$Passed[nrow(rejected$checks)] <- FALSE
  rejected$checks$Detail[nrow(rejected$checks)] <- rejected$check$review
  local_mocked_bindings(mfrm_gpcm_product_inference = function(...) rejected, .package = "mfrmr")
  ci <- suppressWarnings(mfrm_facet_intervals(fit, "Rater"))
  expect_true(all(is.na(ci$table$ModelLower) & is.na(ci$table$Lower)))
  expect_match(ci$table$InferenceReview[1], "Quadrature comparison failed")
  res <- mfrm_results(fit, intervals = list(locations = ci))
  expect_true(all(is.na(res$tables$gpcm_locations_intervals$Lower)))
  expect_true(all(is.na(plot_data(plot(res, draw = FALSE))$table$Lower)))
  expect_match(mfrm_report(res)$markdown, "0 experimental location/contrast intervals are available", fixed = TRUE)
})

test_that("location targets keep arbitrary owner names and shared level labels separate", {
  x <- gmfrm_location_fixture(); d <- x$saved$data
  d$Task <- paste0("shared", match(as.character(d$Task), c("t1", "t2", "t3")))
  d$Rater <- paste0("shared", match(as.character(d$Rater), c("r1", "r2", "r3")))
  names(d)[match(c("Person", "Task", "Rater", "Score"), names(d))] <- c("Candidate", "観点 名", "Judge-ID", "Rating")
  problem <- mfrm_gmfrm_problem(d, 2L, gauss_hermite_normal(31L),
    slope_facets = c("観点 名", "Judge-ID"), person = "Candidate", score = "Rating")
  fit <- mfrm_gmfrm_fit_result(problem,
    mfrm_gmfrm_em(problem, start = x$saved$parameters, maxit = 1L))
  C <- matrix(c(1, -1, 0), 1, dimnames = list("Difference", paste0("shared", 1:3)))
  ci <- suppressWarnings(mfrm_facet_intervals(fit, "Judge-ID", C, level = c(nominal = .9)))
  J <- c(0, 0, 1, -1, rep(0, 9))
  expect_equal(ci$table$Estimate, sum(J * fit$opt$par))
  expect_equal(ci$table$SE^2, drop(J %*% solve(x$saved$information) %*% J), tolerance = 1e-8)
  expect_identical(ci$settings$facet, "Judge-ID")
  expect_identical(ci$settings$level, .9)
  expect_identical(ci$settings$location_reference, "uncentered_on_fixed_N01")
  expect_true(ci$table$CIEligible)
})

test_that("two-family location outputs preserve scale and checks through reports and replay", {
  x <- gmfrm_location_fixture(); fit <- x$fit
  local_mocked_bindings(mfrm_gpcm_product_inference = function(source) {
    expect_identical(source, fit)
    x$inference
  }, .package = "mfrmr")
  ci <- suppressWarnings(mfrm_facet_intervals(fit, "Rater", level = .9))
  slopes <- suppressWarnings(confint(fit, parm = "slopes"))
  local_mocked_bindings(fit_mfrm = function(...) stop("no refitting"),
    compute_mml_parameter_covariance = function(...) stop("no covariance"),
    mfrm_gpcm_product_inference = function(...) stop("no inference"), .package = "mfrmr")
  res <- mfrm_results(fit, intervals = list(locations = ci, slopes = slopes))
  expect_identical(res$tables$gpcm_locations_intervals[names(ci$table)], ci$table)
  expect_identical(res$tables$gpcm_locations_checks, ci$checks)
  expect_equal(as.matrix(res$tables$gpcm_locations_contrasts[-1]), ci$contrasts, ignore_attr = TRUE)
  expect_identical(summary(ci), ci$table)
  expect_equal(apa_table(ci, digits = 10)$table$Lower, ci$table$Lower, tolerance = 1e-9)
  p <- plot_data(plot(res, type = "gpcm_locations", draw = FALSE))
  expect_identical(p$table, ci$table)
  expect_match(p$xlab, "fixed N(0,1)", fixed = TRUE)
  expect_match(p$subtitle, "Uncentered locations")
  expect_match(p$subtitle, "coverage unqualified")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(res, type = "gpcm_locations")
    expect_identical(gg$labels$x, p$xlab)
    expect_identical(plot_data(gg)$table, ci$table)
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
  report <- mfrm_report(res)
  expect_match(report$markdown, "3 experimental location/contrast intervals are available", fixed = TRUE)
  expect_match(report$markdown, "Location differences do not imply uniform")
  expect_identical(report$tables$gpcm_locations_checks, ci$checks)
  folder <- withr::local_tempdir()
  exported <- export_mfrm_results(res, output_dir = folder, include = c("tables", "report", "replay"), acknowledge_sensitive = TRUE)
  csv <- exported$written_files$Path[exported$written_files$Component == "table_gpcm_locations_intervals"]
  expect_length(csv, 1L)
  tab <- read.csv(csv, check.names = FALSE)
  expect_equal(tab$Lower, ci$table$Lower, tolerance = 1e-12)
  expect_identical(tab$CIEligible, ci$table$CIEligible)
  expect_identical(tab$ScaleReference, ci$table$ScaleReference)
  replay <- list.files(folder, pattern = "_replay[.]R$", full.names = FALSE)
  e <- new.env(parent = globalenv())
  withr::with_dir(folder, sys.source(replay, envir = e))
  expect_identical(e$res, res)
  expect_identical(e$report$markdown, report$markdown)
  bad <- fit; bad$opt$par[1] <- bad$opt$par[1] + .1
  expect_error(mfrm_results(bad, intervals = ci), "must match")
})
