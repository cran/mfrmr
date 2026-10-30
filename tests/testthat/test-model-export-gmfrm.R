test_that("two-family exports retain native probabilities and constrained priors", {
  fit <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))$fit
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
    mfrm_gpcm_product_inference = function(...) stop("No inference"), .package = "mfrmr")
  for (binary in c(FALSE, TRUE)) {
    f <- fit
    if (binary) {
      f$config$n_cat <- 2L
      f$config$rating_max <- f$config$rating_min + 1L
      f$config$score_map <- f$config$score_map[1:2, ]
      f$config$step_specs <- mfrmr:::mfrmr_default_step_specs(f$config)
      f$prep$data$score_k <- f$prep$data$score_k %% 2L
    }
    f$config$facet_signs[1] <- 1
    sizes <- mfrmr:::build_param_sizes(f$config)
    par <- sin(seq_len(sum(unlist(sizes)))) * .6
    params <- mfrmr:::expand_params(par, sizes, f$config)
    b <- mfrmr:::mfrm_model_export_data(f, 1.3, .7, .4)
    s <- b$stan_data
    free <- drop(crossprod(s$location_basis, unlist(params$facets, use.names = FALSE)))
    log_slopes <- unlist(lapply(params$slope_components, `[[`, "LogSlope"), use.names = FALSE)
    slope_free <- drop(crossprod(s$slope_basis, log_slopes))
    expect_equal(drop(s$slope_basis %*% slope_free), log_slopes, tolerance = 1e-12)
    expect_equal(diag(s$slope_basis %*% diag(as.vector(s$log_slope_sd)^2) %*% t(s$slope_basis)),
      rep(.4^2, s$A), tolerance = 1e-12)
    expect_equal(sum(log_slopes[b$slope_levels$ScaleReference == "geometric_mean_one"]), 0,
      tolerance = 1e-12)
    theta <- seq(-2, 2, length.out = s$J)
    slope <- exp(drop(s$slope_X %*% slope_free))
    eta <- theta[s$person] + drop(s$X %*% free)
    logits <- (outer(eta, 0:(s$C - 1L)) -
      t(apply(cbind(0, params$steps_mat[s$owner, , drop = FALSE]), 1L, cumsum))) * slope
    probs <- exp(logits - apply(logits, 1L, max))
    probs <- probs / rowSums(probs)
    idx <- mfrmr:::build_indices(f$prep, step_facet = f$config$step_facet,
      slope_facet = f$config$slope_facet, gpcm_spec = f$config$gpcm_spec)
    eta_native <- mfrmr:::compute_eta(idx, params, f$config, theta_override = theta)
    native <- mfrmr:::compute_response_probability_bundle(f$config, idx, params, eta_native)
    expect_equal(unname(probs), unname(native$probs), tolerance = 1e-12)
    expect_identical(b$slope_levels$Owner, rep(f$config$slope_facet,
      lengths(f$config$gpcm_spec$component_levels)))
    expect_equal(s$Q, if (binary) 0L else f$config$n_cat - 2L)
  }
})

test_that("GMFRM Stan bundles protect files and reject unsupported routes before writing", {
  fit <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))$fit
  path <- tempfile(); on.exit(unlink(path, recursive = TRUE), add = TRUE)
  expect_error(export_mfrm_model(1, path), "must be an ordinary mfrm_fit")
  expect_false(dir.exists(path))
  expect_error(export_mfrm_model(fit, path, "brms"), "Stan|stan")
  expect_false(dir.exists(path))
  for (value in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(export_mfrm_model(fit, path, log_slope_prior_sd = value), "positive")
    expect_false(dir.exists(path))
  }
  bad <- fit; bad$config$gpcm_spec$component_levels[[1]] <- rev(bad$config$gpcm_spec$component_levels[[1]])
  expect_error(export_mfrm_model(bad, path), "disagree")
  expect_false(dir.exists(path))
  exported <- export_mfrm_model(fit, path, acknowledge_sensitive = TRUE)
  expect_true(all(file.exists(exported$files)))
  expect_true("slope-levels.csv" %in% names(exported$files))
  expect_identical(read.csv(exported$files[["slope-levels.csv"]]),
    readRDS(exported$files[["model-data.rds"]])$slope_levels)
  expect_silent(parse(exported$files[["run.R"]]))
  runner <- new.env(parent = globalenv())
  # Sourcing helpers must not compile, sample or require external software.
  expect_silent(sys.source(exported$files[["run.R"]], envir = runner))
  expect_true(is.function(runner$summarize_mfrm_stan))
  expect_true(is.function(runner$run_mfrm_stan))
  expect_error(runner$run_mfrm_stan("--unknown"), "Usage:")
  hashes <- tools::md5sum(exported$files)
  expect_error(export_mfrm_model(fit, path, acknowledge_sensitive = TRUE), "already exist")
  expect_identical(tools::md5sum(exported$files), hashes)
  expect_silent(export_mfrm_model(fit, path, overwrite = TRUE, acknowledge_sensitive = TRUE))
})
