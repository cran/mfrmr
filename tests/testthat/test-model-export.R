model_export_fixture <- local({
  cache <- new.env(parent = emptyenv())
  function(model = "RSM") {
    if (!exists(model, cache, inherits = FALSE)) {
      d <- load_mfrmr_data("example_core")
      d <- d[d$Person %in% unique(d$Person)[1:12], ]
      fit <- suppressWarnings(fit_mfrm(
        d, "Person", c("Rater", "Criterion"), "Score", model = model,
        method = "MML", step_facet = if (model == "PCM") "Criterion" else NULL,
        quad_points = 15))
      assign(model, fit, cache)
    }
    get(model, cache, inherits = FALSE)
  }
})

test_that("exported coordinates retain the native response likelihood and priors", {
  for (model in c("RSM", "PCM")) {
    fit <- model_export_fixture(model)
    # Exercise signed and uncentered effects without changing the response data.
    fit$config$facet_signs[["Rater"]] <- 1
    fit$config$facet_specs$Rater$centered <- FALSE
    fit$config$facet_specs$Rater$n_params <- length(fit$config$facet_levels$Rater)
    sizes <- mfrmr:::build_param_sizes(fit$config)
    par <- seq(-0.8, 0.9, length.out = sum(unlist(sizes)))
    # build_param_sizes contains only parameter block lengths.
    params <- mfrmr:::expand_params(par, sizes, fit$config)
    b <- mfrmr:::mfrm_model_export_data(fit, 1.7, 0.6)
    s <- b$stan_data
    locations <- unlist(params$facets, use.names = FALSE)
    free <- drop(crossprod(s$location_basis, locations))
    expect_equal(drop(s$location_basis %*% free), locations, tolerance = 1e-12)
    expect_equal(diag(s$location_basis %*% diag(s$location_sd^2) %*%
                        t(s$location_basis)), rep(1.7^2, s$L), tolerance = 1e-12)
    expect_equal(unname(diag(s$step_basis %*% t(s$step_basis))) * s$step_sd^2,
                 rep(0.6^2, s$C - 1L), tolerance = 1e-12)
    steps <- if (model == "RSM") matrix(params$steps, 1L) else params$steps_mat
    expect_equal(steps %*% s$step_basis %*% t(s$step_basis), unname(steps),
                 ignore_attr = TRUE, tolerance = 1e-12)
    theta <- seq(-2, 1.5, length.out = s$J)
    eta <- theta[s$person] + drop(s$X %*% free)
    logits <- outer(eta, 0:(s$C - 1L)) -
      t(apply(cbind(0, steps[s$owner, , drop = FALSE]), 1L, cumsum))
    probs <- exp(logits - apply(logits, 1L, max))
    probs <- probs / rowSums(probs)
    idx <- mfrmr:::build_indices(fit$prep, step_facet = fit$config$step_facet)
    native_eta <- mfrmr:::compute_eta(idx, params, fit$config, theta_override = theta)
    native <- mfrmr:::compute_response_probability_bundle(fit$config, idx, params, native_eta)
    expect_equal(unname(probs), unname(native$probs), tolerance = 1e-12)
  }
})

test_that("bundles are self-contained, protect files and retain category labels", {
  fit <- model_export_fixture()
  fit$config$score_map$OriginalScore <- c(10, 20, 40, 80)
  for (format in c("stan", "brms")) {
    path <- tempfile()
    on.exit(unlink(path, recursive = TRUE), add = TRUE)
    expect_warning(out <- export_mfrm_model(fit, path, format), "original level labels")
    expect_true(all(file.exists(out$files)))
    b <- readRDS(out$files[["model-data.rds"]])
    expect_equal(b$categories$Score, c(10, 20, 40, 80))
    expect_equal(as.integer(b$stan_data$y), fit$prep$data$score_k + 1L)
    expect_equal(b$levels$Level[b$levels$Facet == "Person"], fit$config$person_levels)
    expect_silent(parse(out$files[["run.R"]]))
    if (format == "brms") expect_silent(parse(out$files[["model.R"]]))
    old <- tools::md5sum(out$files)
    expect_error(export_mfrm_model(fit, path, format), "already exist")
    expect_identical(tools::md5sum(out$files), old)
    expect_silent(export_mfrm_model(fit, path, format, overwrite = TRUE,
                                    acknowledge_sensitive = TRUE))
  }
})

test_that("unsupported models and priors fail before creating output", {
  fit <- model_export_fixture()
  mutations <- list(
    function(f) { f$config$method <- "JML"; f },
    function(f) { f$config$model <- "GPCM"; f },
    function(f) { f$config$population_spec$active <- TRUE; f },
    function(f) { f$config$interaction_specs <- list(a = 1); f },
    function(f) { f$config$dummy_facets <- "Rater"; f },
    function(f) { f$config$facet_shrinkage <- "empirical_bayes"; f },
    function(f) { f$config$facet_specs$Rater$anchors[1] <- 0; f },
    function(f) { f$config$theta_spec$anchors[1] <- 0; f },
    function(f) { f$config$facet_specs$Rater$groups[1] <- "a"; f },
    function(f) { f$prep$data$Weight[1] <- 2; f },
    function(f) { class(f) <- c("mfrm_testlet", class(f)); f })
  for (modify in mutations) {
    path <- tempfile()
    expect_error(export_mfrm_model(modify(fit), path))
    expect_false(file.exists(path))
  }
  for (sd in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(export_mfrm_model(fit, tempfile(), location_prior_sd = sd), "positive")
    expect_error(export_mfrm_model(fit, tempfile(), step_prior_sd = sd), "positive")
  }
})

test_that("binary ladders and centered singleton facets use zero-dimensional bases", {
  fit <- model_export_fixture()
  fit$config$n_cat <- 2L
  fit$config$rating_max <- fit$config$rating_min + 1L
  fit$config$score_map <- fit$config$score_map[1:2, ]
  fit$config$step_specs <- mfrmr:::mfrmr_default_step_specs(fit$config)
  fit$prep$data$score_k <- fit$prep$data$score_k %% 2L
  for (f in fit$config$facet_names) {
    level <- fit$config$facet_levels[[f]][1L]
    fit$config$facet_levels[[f]] <- level
    fit$config$facet_specs[[f]] <- mfrmr:::build_facet_constraint(level)
    fit$prep$data[[f]] <- level
  }
  b <- mfrmr:::mfrm_model_export_data(fit, 1, 1)
  expect_identical(b$stan_data$Q, 0L)
  expect_identical(b$stan_data$P, 0L)
  expect_equal(dim(b$stan_data$X), c(nrow(fit$prep$data), 0))
  expect_equal(dim(b$stan_data$step_basis), c(1, 0))
})
