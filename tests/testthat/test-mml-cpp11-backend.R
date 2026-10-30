make_mml_backend_fixture <- function(model = c("RSM", "PCM")) {
  model <- match.arg(model)
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  prep <- mfrmr:::prepare_mfrm_data(
    d,
    person_col = "Person",
    facet_cols = c("Rater", "Task", "Criterion"),
    score_col = "Score"
  )
  roles <- mfrmr:::resolve_step_and_slope_facets(
    model = model,
    step_facet = if (identical(model, "PCM")) "Criterion" else NULL,
    slope_facet = NULL,
    facet_names = prep$facet_names
  )
  signs <- mfrmr:::build_facet_signs(prep$facet_names)
  idx <- mfrmr:::build_indices(
    prep,
    step_facet = roles$step_facet,
    slope_facet = roles$slope_facet
  )
  cfg <- mfrmr:::build_estimation_config(
    prep = prep,
    model = model,
    method = "MML",
    step_facet = roles$step_facet,
    slope_facet = roles$slope_facet,
    weight_col = NULL,
    facet_signs = signs$signs,
    positive_facets = signs$positive_facets,
    noncenter_facet = "Person",
    dummy_facets = character(0),
    anchor_df = NULL,
    group_anchor_df = NULL,
    population = NULL
  )
  start <- mfrmr:::build_initial_param_vector(cfg$config, cfg$sizes)
  params <- mfrmr:::expand_params(start, cfg$sizes, cfg$config)

  list(
    idx = idx,
    config = cfg$config,
    quad = mfrmr:::gauss_hermite_normal(5L),
    params = params,
    base_eta = mfrmr:::compute_base_eta(idx, params, cfg$config)
  )
}

test_that("cpp11 backend availability and gating are stable", {
  expect_true(mfrmr:::mfrm_cpp11_backend_available())
  expect_true(mfrmr:::mfrm_use_cpp11_backend(list(model = "RSM")))
  expect_true(mfrmr:::mfrm_use_cpp11_backend(list(model = "PCM")))
  expect_false(mfrmr:::mfrm_use_cpp11_backend(list(model = "GPCM")))
  expect_false(mfrmr:::mfrm_use_cpp11_backend(
    list(model = "PCM"),
    include_linear_part = TRUE
  ))

  old_opt <- options(mfrmr.use_cpp11_backend = TRUE)
  on.exit(options(old_opt), add = TRUE)
  expect_true(mfrmr:::mfrm_use_cpp11_backend(list(model = "RSM")))
  expect_true(mfrmr:::mfrm_use_cpp11_backend(list(model = "PCM")))
})

test_that("cpp11 MML kernel matches the pure-R reference for RSM and PCM", {
  skip_if_not(mfrmr:::mfrm_cpp11_backend_available())

  for (model in c("RSM", "PCM")) {
    fixture <- make_mml_backend_fixture(model)

    logprob_r <- mfrmr:::mfrm_mml_logprob_bundle_r(
      idx = fixture$idx,
      config = fixture$config,
      quad = fixture$quad,
      params = fixture$params,
      base_eta = fixture$base_eta,
      include_probs = TRUE
    )
    logprob_cpp <- mfrmr:::mfrm_mml_logprob_bundle_cpp11(
      idx = fixture$idx,
      config = fixture$config,
      quad = fixture$quad,
      params = fixture$params,
      base_eta = fixture$base_eta,
      include_probs = TRUE
    )

    expect_equal(logprob_cpp$log_prob_mat, logprob_r$log_prob_mat, tolerance = 1e-12, info = model)
    expect_equal(logprob_cpp$quad_basis$nodes, logprob_r$quad_basis$nodes, tolerance = 1e-12, info = model)
    expect_identical(logprob_cpp$person_int, logprob_r$person_int)

    for (q in seq_along(logprob_r$prob_list)) {
      expect_equal(
        unname(logprob_cpp$prob_list[[q]]),
        unname(logprob_r$prob_list[[q]]),
        tolerance = 1e-12,
        info = paste(model, q)
      )
    }

    posterior_r <- mfrmr:::mfrm_mml_posterior_bundle(logprob_r)
    posterior_cpp <- mfrmr:::mfrm_mml_posterior_bundle(logprob_cpp)
    expect_equal(
      posterior_cpp$person_bundle$log_marginal,
      posterior_r$person_bundle$log_marginal,
      tolerance = 1e-12,
      info = model
    )
    expect_equal(
      posterior_cpp$obs_posterior,
      posterior_r$obs_posterior,
      tolerance = 1e-12,
      info = model
    )

    expected_r <- mfrmr:::mfrm_mml_expected_category_bundle_r(
      logprob_bundle = logprob_r,
      posterior_bundle = posterior_r,
      include_p_geq = TRUE
    )
    expected_cpp <- mfrmr:::mfrm_mml_expected_category_bundle(
      logprob_bundle = logprob_cpp,
      posterior_bundle = posterior_cpp,
      include_p_geq = TRUE
    )

    expect_equal(
      unname(expected_cpp$posterior_prob),
      unname(expected_r$posterior_prob),
      tolerance = 1e-12,
      info = model
    )
    expect_equal(expected_cpp$expected_k, expected_r$expected_k, tolerance = 1e-12, info = model)
    expect_equal(expected_cpp$var_k, expected_r$var_k, tolerance = 1e-12, info = model)
    expect_equal(
      unname(expected_cpp$p_geq),
      unname(expected_r$p_geq),
      tolerance = 1e-12,
      info = model
    )
  }
})

test_that("automatic MML bundle dispatch uses cpp11 backend for RSM and PCM", {
  skip_if_not(mfrmr:::mfrm_cpp11_backend_available())
  old_opt <- options(mfrmr.use_cpp11_backend = TRUE)
  on.exit(options(old_opt), add = TRUE)

  for (model in c("RSM", "PCM")) {
    fixture <- make_mml_backend_fixture(model)
    auto_bundle <- mfrmr:::mfrm_mml_logprob_bundle(
      idx = fixture$idx,
      config = fixture$config,
      quad = fixture$quad,
      params = fixture$params,
      base_eta = fixture$base_eta,
      include_probs = TRUE
    )
    cpp_bundle <- mfrmr:::mfrm_mml_logprob_bundle_cpp11(
      idx = fixture$idx,
      config = fixture$config,
      quad = fixture$quad,
      params = fixture$params,
      base_eta = fixture$base_eta,
      include_probs = TRUE
    )
    expect_equal(auto_bundle$log_prob_mat, cpp_bundle$log_prob_mat, tolerance = 1e-12, info = model)
  }
})

test_that("shared-grid native probabilities preserve observation-specific results", {
  skip_if_not(mfrmr:::mfrm_cpp11_backend_available())

  for (model in c("RSM", "PCM")) {
    kernel <- if (model == "RSM") {
      mfrmr:::mfrm_cpp_rsm_logprob_bundle
    } else {
      mfrmr:::mfrm_cpp_pcm_logprob_bundle
    }
    for (condition in c("repeated", "unique", "person_nodes", "nonfinite")) {
      n <- 64L
      args <- list(
        score_k = as.integer(rep(0:4, length.out = n)),
        person_int = as.integer(rep(4:1, length.out = n)),
        base_eta = rep(c(-0.5, 0, -0, 0.75), length.out = n),
        person_nodes = matrix(rep(c(-2, 0, 2), each = 4L), nrow = 4L),
        weight = rep(c(0, 0.25, 1, 2), length.out = n)
      )
      if (model == "RSM") {
        args$step_cum <- c(0, -0.2, 0.1, 0.4, 0.8)
      } else {
        args$step_cum_obs <- matrix(rep(c(0, -0.2, 0.1, 0.4, 0.8), each = n), nrow = n)
        args$step_cum_obs[seq(3L, n, by = 3L), 2L] <- 0.125
      }
      if (condition == "unique") args$base_eta <- seq(-1, 1, length.out = n)
      if (condition == "person_nodes") args$person_nodes[2L, ] <- c(-1.5, 0.25, 2.5)
      if (condition == "nonfinite") args$base_eta[1L] <- NA_real_

      for (include_probs in c(FALSE, TRUE)) {
        args$include_probs <- include_probs
        observed <- do.call(kernel, args)
        # Single observation calls use the original uncached loop. They also
        # preserve scores, weights, PCM step ownership and Person-specific nodes.
        singles <- lapply(seq_len(n), function(i) {
          one <- args
          for (field in c("score_k", "person_int", "base_eta", "weight")) {
            one[[field]] <- one[[field]][i]
          }
          if (model == "PCM") one$step_cum_obs <- one$step_cum_obs[i, , drop = FALSE]
          do.call(kernel, one)
        })
        expected <- list(log_prob_mat = do.call(rbind, lapply(singles, `[[`, "log_prob_mat")))
        if (include_probs) {
          expected$prob_list <- lapply(seq_len(ncol(args$person_nodes)), function(q) {
            do.call(rbind, lapply(singles, function(x) x$prob_list[[q]]))
          })
        }
        expect_identical(observed, expected, info = paste(model, condition, include_probs))
      }
    }
  }
})


test_that("native cumulative category probabilities preserve right-to-left sums", {
  skip_if_not(mfrmr:::mfrm_cpp11_backend_available())
  values <- c(0, -0, .25, .Machine$double.xmin, 1e300, -1e300)
  for (shape in list(c(0L, 1L), c(0L, 5L), c(1L, 1L),
                     c(1L, 5L), c(257L, 5L), c(11L, 65L))) {
    probabilities <- matrix(rep(values, length.out = prod(shape)), shape[1L], shape[2L])
    expect_identical(mfrmr:::mfrm_cpp_compute_p_geq(probabilities),
                     mfrmr:::compute_P_geq_r(probabilities))
  }
  # The native zero-category contract predates this loop optimization; the
  # independent R implementation does not accept this input shape.
  expect_identical(mfrmr:::mfrm_cpp_compute_p_geq(matrix(numeric(), 3L, 0L)),
                   matrix(numeric(), 3L, 0L))
})

test_that("native cumulative probabilities preserve nonfinite positions", {
  skip_if_not(mfrmr:::mfrm_cpp11_backend_available())
  values <- c(0, -0, .25, .Machine$double.xmin, 1e300, -1e300,
              NA_real_, NaN, Inf, -Inf)
  for (shape in list(c(257L, 5L), c(11L, 65L))) {
    probabilities <- matrix(rep(values, length.out = prod(shape)), shape[1L], shape[2L])
    actual <- mfrmr:::mfrm_cpp_compute_p_geq(probabilities)
    expected <- mfrmr:::compute_P_geq_r(probabilities)
    # R's ?NA permits either NA or NaN when both participate in arithmetic;
    # compiler and platform choices can differ between R and native loops.
    expect_identical(is.na(actual), is.na(expected))
    expect_identical(actual[!is.na(actual)], expected[!is.na(expected)])
    # The last column is copied without arithmetic, including NA versus NaN.
    expect_identical(actual[, ncol(actual)], probabilities[, ncol(probabilities)])
  }
})
