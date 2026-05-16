test_that("step parameter counts match the sum-zero identification", {
  d <- mfrmr:::sample_mfrm_data(seed = 420)
  fit_args <- list(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "MML",
    quad_points = 5,
    maxit = 12
  )

  fit_rsm <- suppressMessages(suppressWarnings(
    do.call(fit_mfrm, c(fit_args, list(model = "RSM")))
  ))
  fit_pcm <- suppressMessages(suppressWarnings(
    do.call(fit_mfrm, c(fit_args, list(model = "PCM", step_facet = "Criterion")))
  ))
  fit_gpcm <- suppressMessages(suppressWarnings(
    do.call(fit_mfrm, c(fit_args, list(
      model = "GPCM",
      step_facet = "Criterion",
      slope_facet = "Criterion"
    )))
  ))

  n_steps <- fit_rsm$config$n_cat - 1L
  criterion_n <- length(fit_pcm$config$facet_levels$Criterion)

  sizes_rsm <- mfrmr:::build_param_sizes(fit_rsm$config)
  sizes_pcm <- mfrmr:::build_param_sizes(fit_pcm$config)
  sizes_gpcm <- mfrmr:::build_param_sizes(fit_gpcm$config)

  expect_equal(sizes_rsm$steps, max(n_steps - 1L, 0L))
  expect_equal(sizes_pcm$steps, criterion_n * max(n_steps - 1L, 0L))
  expect_equal(sizes_gpcm$steps, criterion_n * max(n_steps - 1L, 0L))
  expect_equal(sizes_gpcm$log_slopes, criterion_n - 1L)

  implied_k <- function(fit) {
    as.numeric((fit$summary$AIC[1] + 2 * fit$summary$LogLik[1]) / 2)
  }
  expect_equal(implied_k(fit_rsm), sum(unlist(sizes_rsm)), tolerance = 1e-8)
  expect_equal(implied_k(fit_pcm), sum(unlist(sizes_pcm)), tolerance = 1e-8)
  expect_equal(implied_k(fit_gpcm), sum(unlist(sizes_gpcm)), tolerance = 1e-8)

  expect_equal(sum(fit_rsm$steps$Estimate), 0, tolerance = 1e-8)
  step_sums <- tapply(fit_pcm$steps$Estimate, fit_pcm$steps$StepFacet, sum)
  expect_true(all(abs(step_sums) < 1e-8))
  gpcm_step_sums <- tapply(fit_gpcm$steps$Estimate, fit_gpcm$steps$StepFacet, sum)
  expect_true(all(abs(gpcm_step_sums) < 1e-8))
  expect_equal(prod(fit_gpcm$slopes$Estimate), 1, tolerance = 1e-8)
})

test_that("GPCM diagnostics expose joint-covariance SEs for steps and slopes", {
  d <- mfrmr:::sample_mfrm_data(seed = 422)
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(
      d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "MML",
      model = "GPCM",
      step_facet = "Criterion",
      slope_facet = "Criterion",
      quad_points = 5,
      maxit = 12
    )
  ))

  dx <- suppressMessages(suppressWarnings(diagnose_mfrm(fit, residual_pca = "none")))
  unc <- dx$parameter_uncertainty

  expect_true(is.list(unc))
  expect_true(all(c("steps", "slopes", "summary") %in% names(unc)))
  expect_true(all(c("SE", "CI_Lower", "CI_Upper", "SE_Method", "SE_Status") %in% names(unc$steps)))
  expect_true(all(c("LogSE", "SE", "LogCI_Lower", "LogCI_Upper", "CI_Lower", "CI_Upper") %in% names(unc$slopes)))
  expect_true(any(is.finite(unc$steps$SE)))
  expect_true(any(is.finite(unc$slopes$SE)))
  expect_true(all(unc$slopes$CI_Lower > 0, na.rm = TRUE))
  expect_true(unc$status %in% c("ok", "regularized"))

  attached <- mfrmr:::attach_diagnostics_to_fit(fit)
  expect_true("SE" %in% names(attached$steps))
  expect_true("SE" %in% names(attached$slopes))
  expect_true("LogSE" %in% names(attached$slopes))
  expect_true(any(is.finite(attached$steps$SE)))
  expect_true(any(is.finite(attached$slopes$SE)))
})

test_that("identified GPCM MML gradient agrees with finite differences", {
  d <- mfrmr:::sample_mfrm_data(seed = 421)
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(
      d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "MML",
      model = "GPCM",
      step_facet = "Criterion",
      slope_facet = "Criterion",
      quad_points = 5,
      maxit = 4
    )
  ))

  config <- fit$config
  sizes <- mfrmr:::build_param_sizes(config)
  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = config$step_facet,
    slope_facet = config$slope_facet,
    interaction_specs = config$interaction_specs
  )
  quad <- mfrmr:::gauss_hermite_normal(5)
  par <- fit$opt$par

  fn <- function(p) mfrmr:::mfrm_loglik_mml(p, idx, config, sizes, quad)
  analytic <- mfrmr:::mfrm_grad_mml(par, idx, config, sizes, quad)
  eps <- 1e-5
  numeric_grad <- vapply(seq_along(par), function(i) {
    p_hi <- p_lo <- par
    p_hi[i] <- p_hi[i] + eps
    p_lo[i] <- p_lo[i] - eps
    (fn(p_hi) - fn(p_lo)) / (2 * eps)
  }, numeric(1))

  expect_equal(analytic, numeric_grad, tolerance = 1e-4)
})
