# Algebra/source-contract tests, not newly optimized or coverage-qualified fits.
fixed_sd_interval_source <- function(anchor_origin = FALSE, sd = 1.4) {
  fit <- readRDS(test_path("fixtures", "native-location-pcm-regression.rds"))$fit
  old <- expand_params(fit$opt$par, build_param_sizes(fit$config), fit$config)
  pop <- fit$population
  pop$fixed_sd <- sd
  pop$sigma2 <- sd^2
  if (anchor_origin) {
    pop$formula <- ~ 1
    pop$terms <- stats::terms(~ 1)
    pop$design_matrix <- matrix(1, nrow(pop$design_matrix), 1L,
      dimnames = list(rownames(pop$design_matrix), "(Intercept)"))
    pop$design_columns <- "(Intercept)"
    pop$coefficients <- c("(Intercept)" = .4)
    fit$config$noncenter_facet <- "Rater"
    levels <- fit$prep$levels$Rater
    anchors <- setNames(c(-.3, .3), levels[1:2])
    fit$config$facet_specs$Rater <- build_facet_constraint(levels,
      anchors = anchors, centered = FALSE,
      centering_scope = fit$config$facet_specs$Rater$centering_scope %||% "remainder")
    old$facets$Rater <- setNames(c(-.3, .3, rep(.5, length(levels) - 2L)), levels)
  }
  fit$population <- pop
  fit$config$population_spec <- compact_population_spec(pop, fit$prep$levels$Person)
  old$population$coefficients <- pop$coefficients
  old$population$sigma2 <- sd^2
  old$population$log_sigma2 <- log(sd^2)
  fit$opt$par <- collapse_expanded_params(old, fit$config)
  fit$facets$others <- build_other_facet_table(fit$config, fit$prep, old)
  fit$steps <- build_step_table(fit$config, fit$prep, old)
  fit
}

fixed_sd_interval_information <- function(fit, eligible = TRUE) {
  p <- length(fit$opt$par)
  V <- diag(.04, p) + outer(seq_len(p), seq_len(p)) * 1e-5
  list(covariance = if (eligible) V else NULL,
    information = list(status = "ok"),
    checks = data.frame(Check = c("Identification and boundary", "Category support",
      "Joint information", "Local solution", "Full gradient", "Stationarity and inverse",
      "Comparison information", "Quadrature comparison"), Passed = eligible,
      Detail = "Algebra-only injected qualification"),
    numerical_checks = data.frame(CurvatureScaledGradient = 0, InverseResidual = 0,
      ComparisonInverseResidual = 0, QuadratureScoreShift = 0, QuadratureCovarianceChange = 0),
    check = list(eligible = eligible, review = "Algebra-only injected qualification"))
}

test_that("fixed population variance has no free covariance or inferential interval", {
  for (anchored in c(FALSE, TRUE)) {
    fit <- fixed_sd_interval_source(anchored)
    inference <- fixed_sd_interval_information(fit)
    local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
      mfrm_native_location_inference = function(...) inference,
      make_mfrm_direct_evaluator = function(...) stop("Fixed SD has no zero-variance alternative"),
      .package = "mfrmr")
    ci <- suppressWarnings(mfrm_population_intervals(fit))
    beta <- build_param_slices(build_param_sizes(fit$config))$beta
    expect_equal(unname(vcov(ci)), unname(inference$covariance[beta, beta, drop = FALSE]))
    expect_identical(rownames(vcov(ci)), paste0("beta:", fit$population$design_columns))
    expect_false("log_residual_variance" %in% colnames(vcov(ci)))
    fixed <- subset(ci$table, Type == "residual_variance")
    expect_equal(fixed$Estimate, 1.4^2)
    expect_equal(fixed$WorkingEstimate, log(1.4^2))
    expect_identical(fixed$Status, "fixed")
    expect_identical(fixed$NaturalScaleStatus, "finite")
    expect_false(fixed$CIEligible)
    expect_true(all(is.na(fixed[c("SE", "Lower", "Upper", "WorkingSE", "WorkingLower", "WorkingUpper")])))
    expect_true(all(ci$table$CIEligible[ci$table$Type == "coefficient"]))
    L <- matrix(2, 1L, length(beta), dimnames = list("double_sum", fit$population$design_columns))
    contrast <- suppressWarnings(mfrm_population_intervals(fit, contrasts = L))
    expect_equal(contrast$contrast_covariance, L %*% vcov(ci) %*% t(L), ignore_attr = TRUE)
    expect_equal(vcov(contrast), vcov(ci))
    expect_identical(tail(contrast$table$Status, 1L), "fixed")
    expect_identical(mfrm_native_location_boundary_candidate(fit, inference), inference)
  }
})

test_that("anchored free-mean facet intervals preserve the origin and fixed SD scaling", {
  fit <- fixed_sd_interval_source(TRUE)
  inference <- fixed_sd_interval_information(fit)
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
    mfrm_native_location_inference = function(...) inference,
    make_mfrm_direct_evaluator = function(...) stop("No boundary evaluation"), .package = "mfrmr")
  native <- suppressWarnings(mfrm_facet_intervals(fit, "Rater"))
  expect_identical(native$settings$location_reference, "direct_anchor_origin")
  expect_identical(native$table$Status[1:2], c("fixed", "fixed"))
  expect_true(all(!native$table$CIEligible[1:2]))
  expect_true(all(is.na(native$table$Lower[1:2])))
  expect_equal(native$table$Estimate[1:2], c(-.3, .3))
  expect_equal(native$table$Estimate[-(1:2)], rep(.5, nrow(native$table) - 2L))
  expect_true(any(grepl("source calibration is not propagated", native$cautions)))
  C <- matrix(0, 1L, nrow(native$table), dimnames = list("new_minus_old", native$table$Target))
  C[1, c(1, ncol(C))] <- c(-1, 1)
  raw <- suppressWarnings(mfrm_facet_intervals(fit, "Rater", C))
  standardized <- suppressWarnings(mfrm_facet_intervals(fit, "Rater", C, scale = "standardized"))
  expect_equal(standardized$table$Estimate, raw$table$Estimate / 1.4)
  expect_equal(standardized$table$SE, raw$table$SE / 1.4)
  expect_equal(standardized$target_jacobian, raw$target_jacobian / 1.4)
  expect_false(standardized$settings$population_sd_estimated)
  expect_identical(standardized$settings$population_sd, 1.4)
  expect_false(any(standardized$table$CoverageQualified))
  direct <- suppressWarnings(mfrm_population_intervals(fit))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(native, path)
  local_mocked_bindings(mfrm_native_location_inference = function(...) stop("No repeated information"),
    .package = "mfrmr")
  reused <- suppressWarnings(mfrm_population_intervals(readRDS(path)))
  expect_identical(reused$table, direct$table)
  expect_identical(vcov(reused), vcov(direct))
  expect_true(reused$settings$inference_reused)
})

test_that("fixed SD source contracts reject changed constraints and variance metadata", {
  fit <- fixed_sd_interval_source(TRUE)
  expect_silent(mfrm_native_location_source(fit))
  source <- mfrm_native_location_identity(fit)
  changed <- fit; changed$population$fixed_sd <- 1
  expect_error(mfrm_native_location_source(changed), "fixed population SD")
  expect_false(identical(mfrm_native_location_identity(changed), source))
  changed <- fit; changed$config$population_spec$sigma2 <- 1
  expect_error(mfrm_native_location_source(changed), "fixed population SD")
  changed <- fit; changed$config$facet_specs$Rater$centered <- TRUE
  expect_error(mfrm_native_location_source(changed), "constraints")
  changed <- fit; changed$config$facet_specs$Rater$anchors[] <- NA_real_
  expect_error(mfrm_native_location_source(changed), "direct anchor")
  changed <- fit; changed$config$facet_specs$Rater$groups[1] <- "g"
  expect_error(mfrm_native_location_source(changed), "without groups")
  changed <- fit; changed$facets$others$Estimate[1] <- 100
  expect_error(mfrm_native_location_source(changed), "Saved locations")
  changed <- fit; changed$config$model <- "GPCM"
  expect_error(mfrm_native_location_source(changed), "RSM/PCM")
})

test_that("failed fixed SD inference retains coefficients and the supplied variance", {
  fit <- fixed_sd_interval_source(TRUE)
  local_mocked_bindings(mfrm_native_location_inference = function(...) fixed_sd_interval_information(fit, FALSE),
    .package = "mfrmr")
  ci <- suppressWarnings(mfrm_population_intervals(fit))
  expect_true(all(is.finite(ci$table$Estimate)))
  expect_true(all(!ci$table$CIEligible))
  expect_true(all(is.na(ci$table$SE)))
  expect_identical(ci$table$Status, c("numerical_review_failed", "fixed"))
  expect_true(all(is.na(vcov(ci))))
})
