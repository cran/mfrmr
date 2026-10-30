test_that("large finite category ties retain normalized probabilities and log likelihoods", {
  for (big in c(1e12, 1e16)) {
    eta <- rep(big / 2, 2L)
    cumulative <- c(0, -big / 2, 0)
    expected <- matrix(c(0, .5, .5), 2L, 3L, byrow = TRUE)
    owner <- c(1L, 1L); score <- c(1L, 2L); weight <- c(.5, 2)
    expect_equal(category_prob_rsm(eta, cumulative), expected)
    expect_equal(category_prob_pcm(eta, matrix(cumulative, 1L), owner), expected)
    expect_equal(category_prob_gpcm(eta, matrix(cumulative, 1L), owner, 1), expected)
    expect_equal(loglik_rsm(eta, score, cumulative, weight), -sum(weight) * log(2))
    expect_equal(loglik_pcm(eta, score, matrix(cumulative, 1L), owner, weight), -sum(weight) * log(2))
    expect_equal(loglik_gpcm(eta, score, matrix(cumulative, 1L), owner, 1,
      weight = weight), -sum(weight) * log(2))
    # Finite tail log probabilities survive even when their probabilities underflow.
    expect_equal(loglik_rsm(eta[1], 0L, cumulative) / big, -1)
    limit <- mfrmr_jml_gpcm_joint_rate_boundary_limit(matrix(c(0, big, big), 1L),
      matrix(c(0, 1, 1), 1L), 1L, 1L, 1, 1, integer(0), 1L, integer(0), 1e-10)
    expect_true(limit$valid)
    expect_equal(limit$log_likelihood, -log(2))
    idx <- list(person = c(1L, 1L), score_k = score, step_idx = owner,
      slope_idx = owner, weight = weight)
    params <- list(steps = c(-big / 2, big / 2),
      steps_mat = matrix(c(-big / 2, big / 2), 1L), slopes = 1)
    for (model in c("RSM", "PCM", "GPCM")) {
      cfg <- list(model = model, n_person = 1L, n_cat = 3L)
      for (keep_probs in c(FALSE, TRUE)) {
        r <- mfrm_mml_logprob_bundle_r(idx, cfg, list(nodes = 0, weights = 1),
          params, eta, include_probs = keep_probs)
        expect_equal(as.numeric(r$log_prob_mat), -weight * log(2))
        if (keep_probs) expect_equal(r$prob_list[[1L]], expected)
        if (model != "GPCM") {
          cpp <- mfrm_mml_logprob_bundle_cpp11(idx, cfg, list(nodes = 0, weights = 1),
            params, eta, include_probs = keep_probs)
          expect_equal(cpp, r, tolerance = 1e-12)
        }
      }
      adaptive <- mfrmr_adaptive_person_kernel(idx, cfg, params, eta, 0, 1)
      at <- adaptive(0)
      expect_equal(at$log_likelihood, -sum(weight) * log(2))
      expect_equal(at$prob_list[[1L]], expected)
      expect_equal(adaptive(0, include_moments = FALSE)$log_likelihood, at$log_likelihood)
    }
  }
})

test_that("posterior normalization preserves ties at large negative log joints", {
  # Person rows are deliberately unsorted; two independent events for Person 2.
  lp <- rbind(c(-5e15, -5e15), c(-1e16, -1e16), c(-5e15, -5e15))
  x <- mfrm_mml_person_bundle(lp, c(2L, 1L, 2L),
    list(log_weights = matrix(0, 2L, 2L)), include_posterior = TRUE)
  expect_equal(unname(x$posterior), matrix(.5, 2L, 2L))
  expect_equal(unname(x$log_posterior), matrix(-log(2), 2L, 2L))
  expect_equal(x$person_ids, c(2L, 1L))
  expect_equal(unname(rowSums(x$posterior)), c(1, 1))
  # A representable difference must also survive the common large offset.
  tilted <- mfrm_mml_person_bundle(matrix(c(-1e16, -1e16 - 4), 1L), 1L,
    list(log_weights = matrix(0, 1L, 2L)), include_posterior = TRUE)
  expect_equal(as.numeric(tilted$posterior), c(plogis(4), plogis(-4)), tolerance = 1e-12)
  weighted <- mfrm_mml_person_bundle(matrix(-1e16, 1L, 2L), 1L,
    list(log_weights = matrix(log(c(.25, .75)), 1L)), include_posterior = TRUE)
  expect_equal(as.numeric(weighted$posterior), c(.25, .75), tolerance = 1e-12)
})

test_that("adaptive posterior review preserves the prior under a constant likelihood", {
  # The normalizing constant can be arbitrarily small without changing the
  # posterior when the likelihood has no ability dependence.
  testthat::local_mocked_bindings(mfrmr_adaptive_person_kernel = function(...)
    function(z, ...) list(log_likelihood = rep(-1e16, length(z)),
      score = -z, information = rep(1, length(z))), .package = "mfrmr")
  out <- mfrmr_adaptive_quadrature_review(list(person = 1L, score_k = 0L),
    list(model = "RSM", n_cat = 2L), list(steps = 0),
    gauss_hermite_normal(7L), "P", c(7L, 15L), base_eta = 0)
  expect_true(all(out$Status == "computed"))
  expect_equal(out$FixedEAP, c(0, 0), tolerance = 1e-12)
  expect_equal(out$AdaptiveEAP, c(0, 0), tolerance = 1e-12)
  expect_equal(out$FixedPosteriorSD, c(1, 1), tolerance = 1e-12)
  expect_equal(out$AdaptivePosteriorSD, c(1, 1), tolerance = 1e-12)
})

test_that("cached and uncached JML and MML retain the same large-logit objective", {
  d <- data.frame(Person = rep(paste0("p", 1:6), each = 2L),
    Rater = rep(c("A", "B"), 6L), Score = rep(c(1L, 0L, 2L, 1L), 3L))
  prep <- prepare_mfrm_data(d, "Person", "Rater", "Score", keep_original = TRUE)
  for (method in c("JML", "MML")) for (model in c("RSM", "PCM", "GPCM")) {
    built <- build_estimation_config(prep, model, method,
      if (model != "RSM") "Rater", if (model == "GPCM") "Rater",
      weight_col = NULL, facet_signs = c(Rater = -1), positive_facets = NULL,
      noncenter_facet = "Person", dummy_facets = NULL, anchor_df = NULL, group_anchor_df = NULL)
    cfg <- built$config; sizes <- built$sizes
    idx <- build_indices(prep, cfg$step_facet, cfg$slope_facet)
    par <- numeric(sum(unlist(sizes))); slices <- build_param_slices(sizes)
    par[slices$Rater] <- -5e15; par[slices$steps] <- -5e15
    quad <- gauss_hermite_normal(1L)
    reference_value <- if (method == "JML") mfrm_loglik_jml(par, idx, cfg, sizes) else
      mfrm_loglik_mml(par, idx, cfg, sizes, quad)
    reference_gradient <- if (method == "JML") mfrm_grad_jml(par, idx, cfg, sizes) else
      mfrm_grad_mml(par, idx, cfg, sizes, quad)
    expect_equal(reference_value, nrow(d) * log(2))
    expect_true(all(is.finite(reference_gradient)))
    for (reuse in c(FALSE, TRUE)) {
      e <- make_mfrm_direct_evaluator(method,
        make_param_cache(sizes, cfg, idx, is_mml = method == "MML"),
        idx, cfg, sizes, quad, reuse_probability_workspace = reuse)
      expect_equal(e$value(par), reference_value, tolerance = 1e-12)
      expect_equal(e$gradient(par), reference_gradient, tolerance = 1e-12)
      expect_identical(e$value(par), e$value(par))
      expect_identical(e$gradient(par), e$gradient(par))
    }
  }
})

test_that("optimizer readiness requires every terminal gradient coordinate", {
  for (basis in c("optimizer_gradient", "relative_loglik", "marginal_score_per_person")) {
    for (g in list(c(0, NA_real_), c(0, Inf), c(NA_real_, NaN), numeric(0), NULL, 0)) {
      result <- build_optimizer_diagnostics(list(convergence = 0L, par = c(0, 0)),
        gradient = g, convergence_basis = basis, gradient_tolerance = 1e-4)
      expect_false(identical(result$ConvergenceSeverity, "pass"))
      expect_true(is.na(result$TerminalGradientSupNorm))
      expect_true(is.na(result$TerminalGradientRMS))
    }
    good <- build_optimizer_diagnostics(list(convergence = 0L, par = c(0, 0)),
      gradient = c(1e-8, -1e-8), convergence_basis = basis, gradient_tolerance = 1e-4)
    expect_identical(good$ConvergenceSeverity, "pass")
    fixed <- build_optimizer_diagnostics(list(convergence = 0L, par = numeric(0)),
      gradient = numeric(0), convergence_basis = basis, gradient_tolerance = 1e-4)
    expect_identical(fixed$ConvergenceSeverity, "pass")
    expect_equal(fixed$TerminalGradientSupNorm, 0)
  }
  for (scale in c(1e-200, 1e200)) {
    metrics <- compute_gradient_metrics(c(scale, -scale, 0))
    expect_equal(metrics$TerminalGradientRMS / scale, sqrt(2/3), tolerance = 1e-12)
    expect_equal(metrics$TerminalGradientSupNorm / scale, 1)
  }
})

test_that("a terminal derivative error cannot qualify an optimizer code-zero result", {
  testthat::local_mocked_bindings(
    make_param_cache = function(...) list(),
    make_mfrm_direct_evaluator = function(...) list(value = function(p) sum(p^2),
      gradient = function(p) stop("terminal derivative unavailable"), diagnostics = function() list()),
    optim = function(par, fn, ...) list(par = par, value = fn(par), convergence = 0L,
      counts = c("function" = 1L, "gradient" = 1L)), .package = "mfrmr")
  for (method in c("JML", "MML")) {
    expect_warning(opt <- run_mfrm_direct_optimization(0, method, list(),
      list(model = "RSM"), list(), quad_points = 3L, maxit = 10L, reltol = 1e-6),
      "complete finite terminal gradient is unavailable")
    expect_identical(opt$convergence, 0L)
    expect_identical(opt$optimizer_diagnostics$ConvergenceReason, "code_zero_unavailable_gradient")
    expect_identical(mfrmr_readiness_numerical_component(opt)$State, "review")
    expect_false(any(opt$optimizer_polish$Stages$InferenceReady))
  }
})

test_that("batched adaptive nodes preserve scalar probabilities and posterior derivatives", {
  withr::local_options(mfrmr.use_cpp11_backend = FALSE)
  for (model in c("RSM", "PCM", "GPCM")) for (n in c(1L, 4L)) {
    idx <- list(person = rep(1L, n), score_k = (seq_len(n) - 1L) %% 4L,
      weight = c(2, 0, .25, 1)[seq_len(n)], step_idx = c(2L, 1L, 2L, 1L)[seq_len(n)],
      slope_idx = c(1L, 2L, 2L, 1L)[seq_len(n)])
    config <- list(model = model, n_cat = 4L, n_person = 1L)
    params <- list(steps = c(-.8, 0, .8),
      steps_mat = rbind(c(-1.2, .3, .9), c(-.4, -.1, .5)), slopes = c(.6, 1.5))
    evaluate <- mfrmr:::mfrmr_adaptive_person_kernel(
      idx, config, params, c(-1, 0, 1, .3)[seq_len(n)], -.8, .4)
    z <- c(25, -.3, -25, 0, 1.2)
    batched <- evaluate(z, include_linear_part = TRUE)
    scalar <- lapply(z, evaluate, include_linear_part = TRUE)
    for (field in c("log_likelihood", "score", "information")) {
      expect_identical(batched[[field]], vapply(scalar, `[[`, 0, field))
    }
    expect_identical(batched$prob_list, lapply(scalar, function(x) x$prob_list[[1L]]))
    expect_identical(batched$linear_part_list, lapply(scalar, function(x) x$linear_part_list[[1L]]))
    # Differentiate the log posterior, independently of the category-moment formulas.
    log_posterior <- function(x) evaluate(x, include_moments = FALSE)$log_likelihood +
      dnorm(x, log = TRUE)
    h <- 1e-3
    z <- c(-1.2, .3, 2)
    score <- (log_posterior(z - 2*h) - 8*log_posterior(z - h) +
      8*log_posterior(z + h) - log_posterior(z + 2*h)) / (12*h)
    information <- (log_posterior(z + 2*h) - 16*log_posterior(z + h) +
      30*log_posterior(z) - 16*log_posterior(z - h) + log_posterior(z - 2*h)) / (12*h^2)
    expect_lt(max(abs(evaluate(z)$score - score)), 1e-8)
    expect_lt(max(abs(evaluate(z)$information - information)), 1e-7)
  }
})
