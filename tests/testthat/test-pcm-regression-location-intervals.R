# Saved synthetic fit and independently integrated joint information; no refits.
pcm_regression_location <- local({
  result <- NULL
  function() {
    if (is.null(result)) {
      fixture <- readRDS(test_path("fixtures", "native-location-pcm-regression.rds"))
      result <<- list(fixture = fixture,
        ci = suppressWarnings(mfrm_facet_intervals(fixture$fit, "Rater")))
    }
    result
  }
})

test_that("regression location intervals retain all population cross-information", {
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
    simulate_mfrm_data = function(...) stop("No simulation"), .package = "mfrmr")
  x <- pcm_regression_location(); ci <- x$ci; fit <- x$fixture$fit
  expect_true(all(ci$checks$Passed))
  expect_true(all(ci$table$CIEligible))
  expect_identical(ci$fit, fit)
  expect_equal(ncol(ci$target_jacobian), 13L)
  expected <- ci$target_jacobian %*% x$fixture$reference_covariance %*% t(ci$target_jacobian)
  expect_equal(ci$covariance, expected, ignore_attr = TRUE, tolerance = 1e-5)
  expect_equal(ci$table$Lower, ci$table$Estimate - qnorm(.975)*ci$table$SE)
  expect_identical(ci$settings$population, "conditional_normal_regression")
  expect_identical(ci$population$Term, fit$population$design_columns)
  expect_equal(ci$population$Coefficient, unname(fit$population$coefficients))
  expect_true(all(ci$population$ResidualVariance == fit$population$sigma2))
  expect_false("Mean" %in% names(ci$population))
  # The full inverse differs from holding the population coordinates fixed.
  V <- ci$parameter_covariance
  H <- solve(V)
  conditional <- solve(H[1:10, 1:10])
  expect_gt(max(abs(V[1:10, 1:10]-conditional)), 1e-5)
})

test_that("regression source inconsistencies fail before information evaluation", {
  fit <- pcm_regression_location()$fixture$fit
  local_mocked_bindings(compute_mml_parameter_covariance = function(...) stop("No information for invalid source"),
    .package = "mfrmr")
  bad <- fit; bad$config$population_spec$design_matrix[, 2] <- 1
  expect_error(mfrm_facet_intervals(bad, "Rater"), "Population coding")
  bad <- fit; bad$population$coefficients[2] <- bad$population$coefficients[2] + .1
  expect_error(mfrm_facet_intervals(bad, "Rater"), "Saved population")
  expect_error(mfrm_facet_intervals(fit, "Rater", method = "sandwich"), "sandwich intervals are unavailable")
  gpcm <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))$fit
  pop <- gpcm$config$population_spec
  pop$design_matrix <- cbind(pop$design_matrix, X=seq_len(nrow(pop$design_matrix)))
  pop$design_columns <- colnames(pop$design_matrix)
  gpcm$config$population_spec <- pop
  expect_error(mfrm_facet_intervals(gpcm, "Rater"), "intercept-only normal GPCM")
})

test_that("saved regression intervals reach tables and plots without recomputation", {
  x <- pcm_regression_location(); ci <- unserialize(serialize(x$ci, NULL))
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
    compute_mml_parameter_covariance = function(...) stop("No covariance on replay"), .package = "mfrmr")
  res <- mfrm_results(ci$fit, intervals=list(locations=ci), include=c("fit", "plots"), compute="never")
  expect_identical(res$tables$facet_locations_intervals$Lower, ci$table$Lower)
  expect_identical(res$tables$facet_locations_population, ci$population)
  expect_identical(mfrm_report(res)$tables$facet_locations_population, ci$population)
  expect_identical(plot_data(plot(res, type="facet_locations", draw=FALSE))$table, ci$table)
  bad <- ci$fit
  bad$config$population_spec$design_matrix[1, 2] <- bad$config$population_spec$design_matrix[1, 2]+.01
  expect_error(mfrm_results(bad, intervals=ci, compute="never"), "must match")
})
