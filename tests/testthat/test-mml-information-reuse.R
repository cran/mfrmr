test_that("scoring and population intervals share only the identical coarse information", {
  original <- compute_mml_parameter_covariance
  calls <- integer()
  local_mocked_bindings(compute_mml_parameter_covariance = function(res, retained_matrices = 0L, ...) {
    calls <<- c(calls, res$config$estimation_control$quad_points)
    original(res, retained_matrices, ...)
  }, fit_mfrm = function(...) stop("No refitting"), .package = "mfrmr")
  for (name in c("native-location-pcm-regression.rds", "native-location-rsm-zero-boundary.rds")) {
    fit <- readRDS(test_path("fixtures", name))$fit
    q <- fit$config$estimation_control$quad_points
    L <- matrix(c(0, 2), 1L, dimnames = list("difference", names(fit$population$coefficients)))
    calls <- integer()
    score <- prediction_source_scoring_readiness(fit)
    ci <- suppressWarnings(mfrm_population_intervals(fit, contrasts = L))
    expect_identical(calls, c(q, q, 2*q-1))
    calls <- integer()
    shared <- mfrm_mml_covariance_reuse(fit)
    expect_identical(prediction_source_scoring_readiness(fit, shared), score)
    reused <- suppressWarnings(mfrm_population_intervals_impl(fit, .95, L, shared))
    expect_identical(reused, ci)
    expect_identical(calls, c(q, 2*q-1))
    path <- tempfile(); withr::defer(unlink(path))
    saveRDS(reused, path)
    expect_identical(readRDS(path), ci)
  }
})

test_that("information reuse rejects changed inputs and honors the current memory budget", {
  fit <- list(opt = list(par = c(1, 2)), prep = list(data = 1:3),
    config = list(quad_points = 61, design = diag(2)))
  calls <- 0L
  local_mocked_bindings(compute_mml_parameter_covariance = function(res) {
    calls <<- calls+1L
    list(status = "ok", value = calls)
  }, .package = "mfrmr")
  shared <- mfrm_mml_covariance_reuse(fit)
  expect_identical(shared(fit), shared(fit)); expect_identical(calls, 1L)
  for (field in c("opt", "prep", "config")) {
    changed <- fit; changed[[field]] <- NULL
    expect_error(shared(changed), "unchanged source fit")
  }
  withr::local_options(mfrmr.max_information_bytes = 1)
  expect_identical(shared(fit)$value, 2L)
  expect_identical(shared(fit)$value, 2L)
})

test_that("information failures and warnings are never hidden by reuse", {
  fit <- list(); mode <- "error"; calls <- 0L
  local_mocked_bindings(compute_mml_parameter_covariance = function(res) {
    calls <<- calls+1L
    if (mode == "error") stop("evaluation failed")
    if (mode == "warning") warning("evaluation warning")
    list(status = if (mode == "refusal") "unavailable" else "ok")
  }, .package = "mfrmr")
  shared <- mfrm_mml_covariance_reuse(fit)
  for (i in 1:2) expect_error(shared(fit), "evaluation failed")
  mode <- "refusal"
  for (i in 1:2) expect_identical(shared(fit)$status, "unavailable")
  mode <- "warning"
  for (i in 1:2) expect_warning(shared(fit), "evaluation warning")
  expect_identical(calls, 6L)
})
