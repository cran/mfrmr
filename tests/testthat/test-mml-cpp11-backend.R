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
  expect_false(mfrmr:::mfrm_use_cpp11_backend(list(model = "RSM")))
  expect_false(mfrmr:::mfrm_use_cpp11_backend(list(model = "PCM")))
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
