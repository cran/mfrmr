residual_bootstrap_fixture <- local({
  cache <- new.env(parent = emptyenv())
  function(model = "RSM") {
    if (!exists(model, cache, inherits = FALSE)) {
      data <- with_preserved_rng_seed(932, {
        d <- expand.grid(Person = paste0("p", 1:80),
          Rater = c("A", "B"), Criterion = c("X", "Y", "Z"))
        theta <- rnorm(80)
        eta <- theta[match(d$Person, unique(d$Person))] +
          rep(c(-.4, .4), each = 80, length.out = nrow(d))
        # Independent adjacent-category recurrence for the generating RSM.
        logp <- cbind(0, eta + .5, 2 * eta)
        p <- exp(logp) / rowSums(exp(logp))
        d$Score <- 1 + rowSums(runif(nrow(d)) > t(apply(p, 1, cumsum)))
        names(d)[1] <- "Candidate ID"
        d
      })
      cache[[model]] <- fit_mfrm(data, "Candidate ID", c("Rater", "Criterion"),
        "Score", model = model, step_facet = if (model == "PCM") "Criterion",
        method = "MML", quad_points = 21, maxit = 300, attach_diagnostics = FALSE)
    }
    cache[[model]]
  }
})

test_that("RSM and PCM references refit once per replicate across all scopes", {
  for (model in c("RSM", "PCM")) {
    fit <- residual_bootstrap_fixture(model)
    before <- serialize(fit, NULL)
    set.seed(701)
    rng <- .Random.seed
    p <- analyze_residual_pca(fit, mode = "both", parallel = TRUE,
      parallel_method = "model_bootstrap", parallel_reps = 3, seed = 193)
    expect_identical(.Random.seed, rng)
    expect_identical(serialize(fit, NULL), before)
    expect_equal(nrow(p$bootstrap_trials), 9L)
    expect_equal(length(unique(p$bootstrap_trials$Seed)), 3L)
    expect_true(all(p$bootstrap_trials$RefitAvailable))
    expect_true(all(p$parallel_status$ParallelAvailable))
    expect_equal(p$parallel_status$SuccessfulParallelReps, rep(3, 3))
    expect_equal(p$overall$parallel$draws |> dim(), c(3, 6))
    expect_equal(p$overall_table$ParallelCutoff,
      apply(p$overall$parallel$draws, 2, quantile, probs = .95, type = 8, names = FALSE))
    expect_identical(p$bootstrap_settings$Model, model)
    expect_identical(p$SupportsFormalInference, FALSE)
    full <- analyze_residual_pca(diagnose_mfrm(fit, diagnostic_mode = "both"), mode = "both")
    expect_identical(p$overall$residual_matrix, full$overall$residual_matrix)
    expect_identical(p$overall_table$Eigenvalue, full$overall_table$Eigenvalue)
    expect_identical(p$by_facet_table$Eigenvalue, full$by_facet_table$Eigenvalue)

    # Reproduce a saved replicate, including nondefault score origins and IDs.
    data <- mfrm_gpcm_bootstrap_generate(fit, p$bootstrap_trials$Seed[1])
    expect_identical(names(data), c("Candidate ID", "Rater", "Criterion", "Score"))
    expect_true(all(data$Score %in% 1:3))
    expect_identical(data$`Candidate ID`, as.character(fit$prep$data$Person))
    replicate <- mfrm_gpcm_bootstrap_refit(fit, data)
    replay <- analyze_residual_pca(replicate)
    expect_equal(replay$overall_table$Eigenvalue, p$overall$parallel$draws[1, ])

    saved <- tempfile(fileext = ".rds")
    saveRDS(p, saved)
    p <- readRDS(saved)
    unlink(saved)
    expect_output(print(summary(p)), "regenerates ratings and refits")
    plotted <- plot_residual_pca(p, plot_type = "parallel_scree", draw = FALSE)
    expect_match(plotted$data$subtitle, "model refitting")
    expect_match(plotted$data$subtitle, "no adjustment for scanning")
    expect_output(print(summary(p)), "does not control the chance of any false flag")
    expect_true(any(grepl("refitted-model", plotted$data$legend$label)))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      converted <- as_ggplot(plotted)
      expect_s3_class(converted, "ggplot")
      expect_equal(ggplot2::ggplot_build(converted)$data[[1]]$y,
        c(p$overall_table$Eigenvalue, p$overall_table$ParallelCutoff, p$overall_table$ParallelMean))
    }
    excess <- plot_residual_pca(p, plot_type = "parallel_excess", draw = FALSE)
    expect_match(excess$data$reference_lines$label, "Refitted-model")
  }
})

test_that("unresolved refits are recorded without successful-subset cutoffs", {
  fit <- residual_bootstrap_fixture()
  calls <- 0L
  local_mocked_bindings(mfrm_gpcm_bootstrap_refit = function(source, data) {
    calls <<- calls + 1L
    if (calls == 2L) stop("Deliberate convergence failure")
    if (calls == 3L) warning("Retained numerical caution")
    source
  })
  p <- analyze_residual_pca(fit, mode = "both", parallel = TRUE,
    parallel_method = "model_bootstrap", parallel_reps = 3, seed = 193)
  expect_identical(calls, 3L)
  expect_false(any(p$parallel_status$ParallelAvailable))
  expect_equal(p$parallel_status$SuccessfulParallelReps, rep(2L, 3))
  expect_true(all(is.na(p$overall$parallel$draws[2, ])))
  expect_match(p$bootstrap_trials$Error[4:6], "Deliberate convergence failure")
  expect_match(p$bootstrap_trials$Warning[7:9], "Retained numerical caution")
  expect_match(p$parallel_status$Error, "reference distribution")
  expect_equal(nrow(p$parallel_overall_table), 0L)
})

test_that("the fitted-model route rejects unsupported estimation targets", {
  fit <- residual_bootstrap_fixture()
  run <- function(x) analyze_residual_pca(x, parallel = TRUE,
    parallel_method = "model_bootstrap", parallel_reps = 1)
  expect_error(run(list(obs = data.frame())), "original fit_mfrm")
  bad <- fit; bad$config$method <- "JML"
  expect_error(run(bad), "RSM/PCM fitted by MML")
  bad <- fit; bad$config$model <- "GPCM"
  expect_error(run(bad), "RSM/PCM fitted by MML")
  for (change in c("weights", "anchors", "population", "quadrature", "shrinkage")) {
    bad <- fit
    if (change == "weights") bad$prep$data$Weight[1] <- 2
    if (change == "anchors") bad$config$replay_inputs$anchors <- data.frame(Value = 0)
    if (change == "population") bad$config$population_spec$active <- TRUE
    if (change == "quadrature") bad$config$estimation_control$mml_integration <- "adaptive"
    if (change == "shrinkage") bad$config$facet_shrinkage <- "normal"
    expect_error(run(bad), "fixed standard-normal population")
  }
  bad <- fit; bad$opt$par[1] <- Inf
  expect_error(run(bad), "numerical convergence")
})

test_that("a PCA failure withholds its scope without erasing valid facet scopes", {
  fit <- residual_bootstrap_fixture()
  original <- compute_pca_overall
  calls <- 0L
  local_mocked_bindings(mfrm_gpcm_bootstrap_refit = function(source, data) source,
    calc_marginal_fit_bundle = function(...) stop("Unused marginal diagnostics requested"),
    compute_pca_overall = function(...) {
      calls <<- calls + 1L
      out <- original(...)
      if (calls == 2L) out$error <- "Deliberate undefined residual correlation"
      out
    })
  p <- analyze_residual_pca(fit, mode = "both", parallel = TRUE,
    parallel_method = "model_bootstrap", parallel_reps = 2, seed = 193)
  expect_equal(p$parallel_status$ParallelAvailable, c(FALSE, TRUE, TRUE))
  expect_true(all(p$bootstrap_trials$RefitAvailable))
  expect_match(p$bootstrap_trials$Error[1], "undefined residual correlation")
  expect_true(all(is.na(p$overall$parallel$draws[1, ])))
})

test_that("the generator shares a Person ability and follows adjacent-category odds", {
  fit <- residual_bootstrap_fixture()
  params <- expand_params(fit$opt$par, build_param_sizes(fit$config), fit$config)
  expected <- with_preserved_rng_seed(57, {
    theta <- rnorm(length(fit$prep$levels$Person))
    d <- fit$prep$data
    eta <- theta[match(d$Person, fit$prep$levels$Person)]
    for (facet in fit$config$facet_names) {
      eta <- eta - params$facets[[facet]][match(d[[facet]], fit$prep$levels[[facet]])]
    }
    # Literal category odds, without the production probability evaluator.
    odds <- cbind(1, exp(eta - params$steps[1]),
      exp(2 * eta - sum(params$steps)))
    p <- odds / rowSums(odds)
    u <- runif(nrow(d))
    1 + (u > p[, 1]) + (u > p[, 1] + p[, 2])
  })
  expect_equal(mfrm_gpcm_bootstrap_generate(fit, 57)$Score, unname(expected))
})

test_that("unavailable observed PCA is distinguished from failed refitting", {
  fit <- residual_bootstrap_fixture()
  original <- compute_pca_overall
  calls <- 0L
  local_mocked_bindings(mfrm_gpcm_bootstrap_refit = function(source, data) source,
    compute_pca_overall = function(...) {
      calls <<- calls + 1L
      out <- original(...)
      if (calls == 1L) out$error <- "Observed pairwise matrix is indefinite."
      out
    })
  p <- analyze_residual_pca(fit, mode = "both", parallel = TRUE,
    parallel_method = "model_bootstrap", parallel_reps = 2, seed = 193)
  expect_true(all(p$bootstrap_trials$RefitAvailable))
  expect_equal(p$parallel_status$ParallelAvailable, c(FALSE, TRUE, TRUE))
  expect_match(p$parallel_status$Error[1], "observed residual PCA is unavailable")
  expect_false(grepl("successful replicates", p$parallel_status$Error[1]))
})

test_that("an unobserved score category stays in the declared model on replay", {
  fit <- residual_bootstrap_fixture()
  data <- mfrm_gpcm_bootstrap_generate(fit, 20)
  data$Score[data$Score == 3] <- 2
  refit <- suppressWarnings(mfrm_gpcm_bootstrap_refit(fit, data))
  expect_identical(refit$prep$score_map, fit$prep$score_map)
  expect_identical(build_param_sizes(refit$config), build_param_sizes(fit$config))
  expect_match(mfrmr_get_readiness_record(refit)$fit$ReasonCodes,
    "declared_category_unobserved")
})

test_that("sparse assignment rows remain fixed instead of being completed", {
  fit <- residual_bootstrap_fixture()
  data <- mfrm_gpcm_bootstrap_generate(fit, 881)
  data <- data[-seq(1, nrow(data), by = 7), ]
  sparse <- fit_mfrm(data, "Candidate ID", c("Rater", "Criterion"), "Score",
    method = "MML", quad_points = 21, maxit = 300, attach_diagnostics = FALSE)
  generated <- mfrm_gpcm_bootstrap_generate(sparse, 84)
  expect_equal(nrow(generated), nrow(data))
  expect_identical(generated[c("Candidate ID", "Rater", "Criterion")],
    data[c("Candidate ID", "Rater", "Criterion")] |> `rownames<-`(NULL))
  p <- analyze_residual_pca(sparse, parallel = TRUE,
    parallel_method = "model_bootstrap", parallel_reps = 2, seed = 193)
  expect_true(anyNA(p$overall$residual_matrix))
  expect_equal(nrow(p$bootstrap_trials), 2L)
  expect_true(all(p$bootstrap_trials$RefitAvailable))
  expect_false(any(grepl("observation pattern", p$bootstrap_trials$Error)))
})
