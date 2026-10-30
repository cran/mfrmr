# Reuse retained fits. These checks must never optimize or simulate new data.
native_location_fixture <- local({
  saved <- list()
  function(name = c("separate", "paired", "single", "shared", "regression", "pcm", "fixed")) {
    name <- match.arg(name)
    if (!is.null(saved[[name]])) return(saved[[name]])
    if (name == "separate") {
      fit <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))$fit
    } else if (name == "fixed") {
      fit <- readRDS(test_path("fixtures", "native-location-fixed-gpcm.rds"))
    } else if (name == "pcm") {
      path <- test_path("..", "..", "validation-results", "gpcm-separated-owners-20260926", "fit", "pcm.rds")
      skip_if_not(file.exists(path), "Repository-only retained intercept-population PCM fit")
      fit <- readRDS(path)$value
    } else if (name == "shared") {
      path <- test_path("..", "..", "validation-results", "gpcm-probability-refit-20260925", "refit-788.rds")
      skip_if_not(file.exists(path), "Repository-only retained shared-owner fit")
      fit <- readRDS(path)$fit
    } else {
      path <- test_path("..", "..", "inst", "validation", "population-full-information-evidence-0.2.4.rds")
      skip_if_not(file.exists(path), "Repository-only independent population information archive")
      key <- switch(name, paired = "paired-61", single = "single-31", regression = "PCM-61")
      fit <- readRDS(path)$results[[key]]$result$fit
    }
    saved[[name]] <<- fit
    fit
  }
})

# Share interval qualification across tests; the finer-grid failure test also
# evaluates the same retained separate-owner Hessian once. Replay tests below
# replace fitting/covariance entry points with errors.
native_location_result <- local({
  saved <- list()
  function(name = "separate") {
    if (is.null(saved[[name]])) saved[[name]] <<- suppressWarnings(
      mfrm_facet_intervals(native_location_fixture(name), "Rater"))
    saved[[name]]
  }
})

test_that("population-SD contrasts propagate the full derivative and cross-covariance", {
  bases <- lapply(c("separate", "paired", "shared", "pcm", "fixed"), native_location_result)
  for (base in bases) {
    fit <- base$fit
    n <- nrow(base$table)
    C <- rbind(Difference = c(1, -1, rep(0, n - 2)),
      Reverse = c(-1, 1, rep(0, n - 2)), Double = c(2, -2, rep(0, n - 2)))
    colnames(C) <- base$table$Target
    V <- base$parameter_covariance
    local_mocked_bindings(mfrm_native_location_inference = function(...) list(
      covariance = V, information = list(status = "ok"), checks = base$checks,
      numerical_checks = base$numerical_checks,
      check = list(eligible = TRUE, review = "Retained native qualification")), .package = "mfrmr")
    ci <- suppressWarnings(mfrm_facet_intervals(fit, "Rater", C[, n:1],
      scale = "standardized", level = .9))
    sizes <- build_param_sizes(fit$config)
    # Independently differentiate the nonlinear expanded-parameter target.
    value <- function(p) {
      par <- expand_params(p, sizes, fit$config)
      drop(C %*% par$facets$Rater) / sqrt(par$population$sigma2 %||% 1)
    }
    h <- 1e-5
    J <- vapply(seq_along(fit$opt$par), function(j) {
      delta <- numeric(length(fit$opt$par)); delta[j] <- h
      (value(fit$opt$par + delta) - value(fit$opt$par - delta)) / (2 * h)
    }, numeric(nrow(C)))
    expect_equal(unname(ci$target_jacobian), unname(J), tolerance = 1e-8)
    expect_equal(ci$covariance, J %*% V %*% t(J), ignore_attr = TRUE, tolerance = 1e-8)
    expect_equal(ci$table$Estimate, unname(value(fit$opt$par)))
    expect_equal(ci$table$NativeEstimate, drop(C %*% base$table$Estimate), ignore_attr = TRUE)
    expect_equal(ci$table$Lower[1], -ci$table$Upper[2])
    expect_equal(ci$table$SE[3], 2 * ci$table$SE[1])
    expect_true(all(ci$table$CIEligible))
    expect_false(any(ci$table$CoverageQualified))
    expect_true(all(grepl("Population-SD-standardized", ci$table$InferenceReview)))
    expect_false(any(grepl("not standardized", ci$cautions)))
    expect_identical(ci$fit, fit)
    if (fit$config$model %in% c("RSM", "PCM")) {
      pop_native <- suppressWarnings(mfrm_population_intervals(base))
      pop_std <- suppressWarnings(mfrm_population_intervals(ci))
      expect_equal(pop_std$table, pop_native$table)
      expect_false(any(grepl("Population-SD-standardized contrasts", pop_std$cautions)))
    }
    var_d <- drop(C[1, , drop = FALSE] %*% base$covariance %*% C[1, ])
    sd <- ci$settings$population_sd
    d <- ci$table$NativeEstimate[1]
    if (isTRUE(fit$config$population_spec$active)) {
      at <- build_param_slices(sizes)$log_sigma2
      cov_d_lambda <- drop(C[1, , drop = FALSE] %*% base$target_jacobian %*% V[, at])
      expected <- (var_d + d^2 * V[at, at] / 4 - d * cov_d_lambda) / sd^2
      expect_equal(ci$table$SE[1]^2, expected, ignore_attr = TRUE)
      expect_gt(abs(expected - var_d / sd^2), 1e-7)
      expect_gt(abs(d * cov_d_lambda), 1e-7)
    } else {
      expect_identical(sd, 1)
      expect_equal(ci$table$SE[1]^2, unname(var_d))
    }
  }
})

test_that("zero estimated contrasts remain uncertain and overflow retains only points", {
  base <- native_location_result("separate"); b <- base$table$Estimate
  C <- rbind(ZeroEstimate = c(b[2] - b[3], b[3] - b[1], b[1] - b[2]),
    Overflow = c(1e200, -1e200, 0))
  colnames(C) <- base$table$Target
  local_mocked_bindings(mfrm_native_location_inference = function(...) list(
    covariance = base$parameter_covariance, information = list(status = "ok"),
    check = list(eligible = TRUE, review = "Retained native qualification")), .package = "mfrmr")
  ci <- suppressWarnings(mfrm_facet_intervals(base$fit, "Rater", C, scale = "standardized"))
  expect_equal(ci$table$Estimate[1], 0, tolerance = 1e-12)
  expect_gt(ci$table$SE[1], 0)
  expect_true(ci$table$CIEligible[1])
  expect_false(ci$table$CIEligible[2])
  expect_true(is.finite(ci$table$Estimate[2]))
  expect_true(all(is.na(ci$table[2, c("SE", "Lower", "Upper", "ModelSE", "ModelLower", "ModelUpper")])))
})

test_that("standardization rejects ambiguous targets and preserves numerical refusal", {
  fit <- native_location_fixture()
  C <- matrix(c(1, -1, 0), 1, dimnames = list("Difference", fit$config$facet_levels$Rater))
  local_mocked_bindings(compute_mml_parameter_covariance = function(...) stop("Unexpected information"),
    mfrm_native_location_inference = function(...) list(covariance = NULL,
      information = list(status = "unavailable"), check = list(eligible = FALSE,
        review = "Retained numerical refusal")), .package = "mfrmr")
  expect_error(mfrm_facet_intervals(fit, "Rater", scale = "standardized"), "explicit zero-sum")
  expect_error(mfrm_facet_intervals(fit, "Rater", C, scale = "standardized", method = "sandwich"), "method")
  expect_error(mfrm_facet_intervals(fit, "Rater", C + 1, scale = "standardized"), "zero-sum")
  expect_error(mfrm_facet_intervals(native_location_fixture("regression"), "Rater", C,
    scale = "standardized"), "residual SD")
  failed <- suppressWarnings(mfrm_facet_intervals(fit, "Rater", C, scale = "standardized"))
  expect_true(is.finite(failed$table$Estimate))
  expect_false(failed$table$CIEligible)
  expect_identical(failed$table$Status, "numerical_review_failed")
  expect_true(all(is.na(failed$table[c("SE", "Lower", "Upper", "ModelSE", "ModelLower", "ModelUpper")])))
  expect_match(failed$table$InferenceReview, "numerical refusal")
})

test_that("standardized intervals survive public calculation, saved reports and exports", {
  fit <- native_location_fixture()
  C <- matrix(c(1, -1, 0), 1, dimnames = list("Difference", fit$config$facet_levels$Rater))
  # One end-to-end call retains the real numerical admission checks.
  ci <- suppressWarnings(mfrm_facet_intervals(fit, "Rater", C, scale = "standardized"))
  expect_true(ci$table$CIEligible)
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
    compute_mml_parameter_covariance = function(...) stop("No information during replay"), .package = "mfrmr")
  folder <- tempfile("standardized-contrast-"); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  saveRDS(ci, file.path(folder, "ci.rds")); saved <- readRDS(file.path(folder, "ci.rds"))
  res <- mfrm_results(fit, intervals = list(standardized = saved), include = c("fit", "plots"), compute = "never")
  key <- "gpcm_standardized_intervals"
  expect_identical(res$tables[[key]]$Estimate, ci$table$Estimate)
  expect_identical(mfrm_report(res)$tables[[key]], res$tables[[key]])
  p <- plot_data(plot(saved, draw = FALSE))
  expect_identical(p$xlab, "Contrast (population SD units)")
  expect_match(p$subtitle, "estimated population SD uncertainty")
  expect_false(grepl("not divided", p$subtitle))
  expect_identical(p$table, ci$table)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(saved)
    expect_null(gg$labels$title); expect_null(gg$labels$subtitle); expect_null(gg$labels$caption)
    expect_identical(plot_data(gg)$table, ci$table)
  }
  suppressWarnings(export_mfrm_results(res, output_dir = file.path(folder, "export"),
    include = c("tables", "report", "replay"), preset = "starter", acknowledge_sensitive = TRUE))
  files <- list.files(file.path(folder, "export"), recursive = TRUE, full.names = TRUE)
  replay <- readRDS(files[grepl("\\.rds$", files)][1])
  expect_identical(replay$gpcm_inference$standardized$table, ci$table)
  expect_identical(replay$gpcm_inference$standardized$target_jacobian, ci$target_jacobian)
  fit$population$sigma2 <- 9
  expect_error(mfrm_results(fit, intervals = ci, compute = "never"), "must match")
})

test_that("native locations retain full nuisance covariance at fixed fitted points", {
  skip_if_private_validation_missing(c(
    "population-full-information-evidence-0.2.4.rds"
  ))
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
    simulate_mfrm_data = function(...) stop("No simulation"),
    mfrm_mml_person_scores_numeric = function(...) stop("No empirical Person-score rank gate"),
    .package = "mfrmr")
  for (name in c("separate", "paired", "shared", "pcm", "fixed")) {
    fit <- native_location_fixture(name); before <- fit
    ci <- native_location_result(name)
    expect_true(all(ci$checks$Passed))
    expect_true(all(ci$table$CIEligible))
    expect_identical(ci$settings$procedure, "mml_native_location_model_v1")
    expect_identical(ci$settings$scale, "native_ability")
    expect_identical(ci$settings$comparison_quad_points, 2 * ci$settings$quad_points - 1)
    expect_identical(ci$fit, before)
    expect_identical(fit, before)
    expect_identical(ci$fit$readiness, fit$readiness)
    expect_equal(ci$table$Lower, ci$table$Estimate - qnorm(.975) * ci$table$SE)
    sizes <- build_param_sizes(fit$config)
    slice <- build_param_slices(sizes)$Rater
    n <- length(fit$prep$levels$Rater)
    J <- matrix(0, n, length(fit$opt$par))
    J[, slice] <- rbind(diag(n - 1L), rep(-1, n - 1L))
    expect_equal(unname(ci$target_jacobian), J)
    expect_equal(ci$covariance, J %*% ci$parameter_covariance %*% t(J), ignore_attr = TRUE)
    expect_match(paste(ci$cautions, collapse = " "), "not standardized")
    expect_null(ci$person_scores)
    expect_null(ci$score_rank)
  }
  # Independently differentiated continuous-information reference, not a
  # second call to the package Hessian used to construct the interval.
  archive <- readRDS(test_path("..", "..", "inst", "validation",
    "population-full-information-evidence-0.2.4.rds"))$results[["paired-61"]]$result
  ci <- native_location_result("paired")
  J <- ci$target_jacobian
  expect_equal(ci$covariance, J %*% archive$reference$covariance %*% t(J),
    ignore_attr = TRUE, tolerance = 1e-5)
  expect_gt(abs(ci$parameter_covariance[1, 1] - 1 / archive$reference$hessian[1, 1]), 1e-5)
})

test_that("fixed-normal native locations recognize the public fit metadata", {
  fit <- native_location_fixture("fixed")
  expect_false(fit$population$active)
  expect_identical(fit$config$posterior_basis, "legacy_mml")
  expect_identical(fit$config$gpcm_mml_identification, "fixed_standard_normal")
  expect_true(all(native_location_result("fixed")$table$CIEligible))
  local_mocked_bindings(compute_mml_parameter_covariance = function(...) stop("No information for inconsistent coding"),
    .package = "mfrmr")
  for (change in list(
    function(f) { f$config$posterior_basis <- "fixed_standard_normal"; f },
    function(f) { f$config$gpcm_mml_identification <- "free_population"; f },
    function(f) { f$population$active <- TRUE; f })) {
    expect_error(mfrm_facet_intervals(change(fit), "Rater"), "Population coding")
  }
})

test_that("fixed and redundant native contrasts remain pointwise targets", {
  base <- native_location_result("separate")
  V <- base$parameter_covariance
  local_mocked_bindings(mfrm_native_location_inference = function(...) list(
    covariance = V, information = list(status = "ok"), checks = base$checks,
    check = list(eligible = TRUE, review = "Saved positive-control qualification")), .package = "mfrmr")
  C <- rbind(Difference = c(1, -1, 0), Reverse = c(-1, 1, 0), Sum = c(1, 1, 1))
  colnames(C) <- base$table$Target
  ci <- suppressWarnings(mfrm_facet_intervals(base$fit, "Rater", C[, 3:1], level = .9))
  expect_equal(ci$table$Estimate, as.numeric(C %*% base$table$Estimate))
  expect_equal(ci$covariance, C %*% base$covariance %*% t(C), ignore_attr = TRUE)
  expect_equal(ci$table$Lower[1], -ci$table$Upper[2])
  expect_identical(ci$table$Status, c("available", "available", "fixed"))
  expect_identical(ci$table$CIEligible, c(TRUE, TRUE, FALSE))
  expect_equal(ci$table$SE[3], 0)
  expect_true(is.na(ci$table$Lower[3]))
  expect_true(all(ci$covariance[3, ] == 0))
})

test_that("scope and altered metadata fail before information evaluation", {
  fit <- native_location_fixture()
  local_mocked_bindings(compute_mml_parameter_covariance = function(...) stop("Unexpected information evaluation"),
    .package = "mfrmr")
  expect_error(mfrm_facet_intervals(fit, "Rater", method = "sandwich"), "sandwich intervals are unavailable")
  expect_silent(mfrm_native_location_source(native_location_fixture("regression")))
  changes <- list(
    locations = function(f) { f$facets$others$Estimate[1] <- 999; f },
    steps = function(f) { f$steps$Estimate[1] <- 999; f },
    slopes = function(f) { f$slopes$Estimate[1] <- 999; f },
    population = function(f) { f$population$sigma2 <- 999; f },
    coding = function(f) { f$config$population_spec$person_lookup[1] <- NA_integer_; f },
    integration = function(f) { f$config$estimation_control$mml_integration <- "unknown"; f },
    constraints = function(f) { f$config$facet_specs$Rater$centered <- FALSE; f },
    weights = function(f) { f$prep$data$Weight[1] <- 2; f })
  patterns <- c("Saved locations", "Saved steps", "Saved slopes", "Saved population", "Population coding", "integration metadata", "sum-to-zero", "unit weights")
  for (i in seq_along(changes)) expect_error(
    mfrm_facet_intervals(changes[[i]](fit), "Rater"), patterns[i])
})

test_that("known nonidentification keeps points without differentiating", {
  fit <- native_location_fixture("single")
  local_mocked_bindings(compute_mml_parameter_covariance = function(...) stop("A known failure must stop first"),
    .package = "mfrmr")
  ci <- suppressWarnings(mfrm_facet_intervals(fit, "Rater"))
  expect_true(all(is.finite(ci$table$Estimate)))
  expect_true(all(!ci$table$CIEligible))
  expect_true(all(is.na(ci$table$SE) & is.na(ci$table$Lower) & is.na(ci$table$ModelLower)))
  expect_identical(ci$table$Status, rep("numerical_review_failed", nrow(ci$table)))
  expect_match(ci$table$InferenceReview[1], "identification failure")
  expect_null(ci$parameter_covariance)
})

test_that("finer-grid failures cannot leak model bounds", {
  fit <- native_location_fixture()
  ci <- native_location_result()
  information <- compute_mml_parameter_covariance(fit)
  expect_length(information$solution_information$gradient, length(fit$opt$par))
  q <- fit$config$estimation_control$quad_points
  mode <- "unavailable"
  local_mocked_bindings(compute_mml_parameter_covariance = function(res, retained_matrices = 0L) {
    expect_identical(res$opt$par, fit$opt$par)
    if (res$config$estimation_control$quad_points == q) return(information)
    expect_identical(res$config$estimation_control$quad_points, 2 * q - 1)
    expect_equal(retained_matrices, 2L)
    if (mode == "unavailable") return(list(status = "unavailable", detail = "Injected finer-grid failure"))
    if (mode == "nonfinite") {
      information$solution_information$gradient[1] <- Inf
    } else {
      information$cov <- 1.05 * information$cov
      information$hessian <- information$hessian / 1.05
    }
    information
  }, .package = "mfrmr")
  failed <- suppressWarnings(mfrm_facet_intervals(fit, "Rater"))
  expect_equal(failed$table$Estimate, ci$table$Estimate)
  expect_true(all(is.na(failed$table$Lower) & is.na(failed$table$ModelLower)))
  expect_match(failed$table$InferenceReview[1], "Injected finer-grid failure")
  expect_false(tail(failed$checks$Passed, 1))
  for (next_mode in c("nonfinite", "changed")) {
    mode <- next_mode
    failed <- suppressWarnings(mfrm_facet_intervals(fit, "Rater"))
    expect_equal(failed$table$Estimate, ci$table$Estimate)
    expect_true(all(!failed$table$CIEligible & is.na(failed$table$SE)))
    expect_false(tail(failed$checks$Passed, 1))
  }
  expect_gt(failed$numerical_checks$QuadratureCovarianceChange, .01)
})

test_that("retained information matrices count against the same workspace budget", {
  fit <- native_location_fixture()
  p <- length(fit$opt$par)
  # Ten dense matrices do not fit; reject before the Hessian is requested.
  withr::local_options(mfrmr.max_information_bytes = 9 * 8 * p^2)
  local_mocked_bindings(make_param_cache = function(...) stop("No evaluator allocation expected"),
    .package = "mfrmr")
  info <- compute_mml_parameter_covariance(fit, retained_matrices = 2L)
  expect_identical(info$status, "unavailable")
  expect_null(info$cov)
  expect_match(info$detail, "No matrix was allocated")
})

test_that("native interval reports and exports replay without covariance calculations", {
  # Resolve all retained results before forbidding numerical work during replay.
  inputs <- lapply(c("separate", "paired", "shared", "pcm", "fixed"), native_location_result)
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
    compute_mml_parameter_covariance = function(...) stop("No covariance evaluation"), .package = "mfrmr")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  for (ci in inputs) {
    gpcm <- identical(ci$fit$config$model, "GPCM")
    key <- if (gpcm) "gpcm_locations" else "facet_locations"
    table_key <- paste0(key, "_intervals")
    saveRDS(ci, path); ci <- readRDS(path)
    res <- mfrm_results(ci$fit, intervals = list(locations = ci), include = c("fit", "plots"), compute = "never")
    expect_identical(res$tables[[table_key]]$Lower, ci$table$Lower)
    p <- plot_data(plot(res, type = key, draw = FALSE))
    expect_identical(p$table, ci$table)
    expect_match(p$subtitle, "Native locations: asymptotic Wald intervals")
    expect_match(p$xlab, "native ability units")
    expect_identical(mfrm_report(res)$tables[[table_key]], res$tables[[table_key]])
    for (change in list(
      function(f) { f$config$facet_specs$Rater$centered <- FALSE; f },
      function(f) { f$config$estimation_control$quad_points <- 121L; f },
      function(f) { f$population$sigma2 <- 9; f })) {
      expect_error(mfrm_results(change(ci$fit), intervals = ci, compute = "never"), "must match")
    }
    folder <- tempfile("native-location-export-")
    on.exit(unlink(folder, recursive = TRUE), add = TRUE)
    suppressWarnings(export_mfrm_results(res, output_dir = folder,
      include = c("tables", "report", "replay"), preset = "starter", acknowledge_sensitive = TRUE))
    files <- list.files(folder, recursive = TRUE, full.names = TRUE)
    replay <- readRDS(files[grepl("\\.rds$", files)][1])
    saved <- if (gpcm) replay$gpcm_inference$locations else replay$facet_intervals$locations
    expect_identical(saved$table, ci$table)
    expect_identical(saved$source, ci$source)
  }
})
