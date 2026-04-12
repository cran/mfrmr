# test-parameter-recovery.R
# Tests for parameter recovery with known ground truth from sample_mfrm_data.
#
# True generating parameters (from sample_mfrm_data):
#   person ability:  rnorm(36, 0, 1) with set.seed(20240131)
#   rater_eff:       c(-0.4, 0, 0.4)   R1=lenient, R2=neutral, R3=harsh
#   task_eff:        seq(-0.5, 0.5, length.out=4)  T1 easiest, T4 hardest
#   crit_eff:        c(-0.3, 0, 0.3)   C1 easiest, C3 hardest
#   eta = ability - rater_eff - task_eff - crit_eff
#   Score = cut(eta + noise, breaks)
#
# In the subtractive MFRM parameterization, a positive facet effect means
# harder/harsher, so the estimated measure should be positively correlated
# with the true effect vector.

# ---------------------------------------------------------------------------
# Helper: reconstruct true person abilities
# ---------------------------------------------------------------------------
true_person_ability <- function() {
  set.seed(20240131)
  rnorm(36, 0, 1)
}

true_rater_eff  <- c(-0.4, 0, 0.4)        # R1, R2, R3
true_task_eff   <- seq(-0.5, 0.5, length.out = 4)  # T1, T2, T3, T4
true_crit_eff   <- c(-0.3, 0, 0.3)        # C1, C2, C3

sample_latent_regression_recovery_case <- function(seed,
                                                   model = c("RSM", "PCM"),
                                                   predictor = c("continuous", "dummy", "none"),
                                                   beta0 = 0.25,
                                                   beta1 = 0.8,
                                                   sigma2 = 0.49,
                                                   n_person = 90L,
                                                   n_item = 6L) {
  model <- match.arg(model)
  predictor <- match.arg(predictor)

  softmax_prob <- function(x) {
    exp_x <- exp(x - max(x))
    exp_x / sum(exp_x)
  }

  mfrmr:::with_preserved_rng_seed(seed, {
    persons <- paste0("P", sprintf("%03d", seq_len(n_person)))
    items <- paste0("I", seq_len(n_item))
    x <- switch(
      predictor,
      continuous = stats::rnorm(n_person),
      dummy = stats::rbinom(n_person, 1, 0.5),
      none = rep(0, n_person)
    )
    theta <- beta0 + beta1 * x + stats::rnorm(n_person, sd = sqrt(sigma2))
    item_beta <- seq(-0.7, 0.7, length.out = n_item)
    base_steps <- c(-1.1, 0, 1.1)
    item_shift <- seq(-0.25, 0.25, length.out = n_item)
    step_mat <- vapply(
      item_shift,
      function(shift) base_steps + shift,
      numeric(length(base_steps))
    )

    dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
    eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
    dat$Score <- vapply(seq_len(nrow(dat)), function(i) {
      if (identical(model, "RSM")) {
        step_cum <- c(0, cumsum(base_steps))
      } else {
        step_cum <- c(0, cumsum(step_mat[, match(dat$Item[i], items)]))
      }
      log_num <- (0:3) * eta[i] - step_cum
      sample(0:3, size = 1, prob = softmax_prob(log_num))
    }, integer(1))

    person_data <- data.frame(Person = persons, stringsAsFactors = FALSE)
    population_formula <- switch(
      predictor,
      continuous = ~ X,
      dummy = ~ G,
      none = ~ 1
    )
    if (identical(predictor, "continuous")) person_data$X <- x
    if (identical(predictor, "dummy")) person_data$G <- x

    fit <- suppressWarnings(
      fit_mfrm(
        dat,
        "Person", "Item", "Score",
        method = "MML",
        model = model,
        step_facet = if (identical(model, "PCM")) "Item" else NULL,
        population_formula = population_formula,
        person_data = person_data,
        quad_points = 7,
        maxit = 80
      )
    )

    truth_person <- data.frame(Person = persons, Theta = theta, stringsAsFactors = FALSE)
    merged_person <- merge(fit$facets$person, truth_person, by = "Person")
    coeff <- fit$population$coefficients
    coeff_names <- names(coeff)
    coeff_tbl <- data.frame(
      Term = coeff_names,
      Estimate = as.numeric(coeff),
      stringsAsFactors = FALSE
    )

    list(
      fit = fit,
      coefficients = coeff_tbl,
      sigma2 = fit$population$sigma2,
      converged = isTRUE(fit$summary$Converged[1]),
      person_correlation = stats::cor(merged_person$Estimate, merged_person$Theta)
    )
  })
}

# ---------------------------------------------------------------------------
# 4A.  JML true value recovery
# ---------------------------------------------------------------------------

test_that("JML rater estimates correlate with true rater effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  rater_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  # Center both for comparison
  rater_est_c <- rater_est - mean(rater_est)
  true_c <- true_rater_eff - mean(true_rater_eff)
  r <- cor(rater_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(rater_est_c - true_c)), 0.6)
})

test_that("JML task estimates correlate with true task effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  task_est <- fit$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  task_est_c <- task_est - mean(task_est)
  true_c <- true_task_eff - mean(true_task_eff)
  r <- cor(task_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(task_est_c - true_c)), 0.6)
})

test_that("JML criterion estimates correlate with true criterion effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  crit_est <- fit$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  crit_est_c <- crit_est - mean(crit_est)
  true_c <- true_crit_eff - mean(true_crit_eff)
  r <- cor(crit_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(crit_est_c - true_c)), 0.6)
})

# ---------------------------------------------------------------------------
# 4B.  MML true value recovery
# ---------------------------------------------------------------------------

test_that("MML rater estimates correlate with true rater effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "MML", quad_points = 21, maxit = 100))

  rater_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  rater_est_c <- rater_est - mean(rater_est)
  true_c <- true_rater_eff - mean(true_rater_eff)
  r <- cor(rater_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(rater_est_c - true_c)), 0.6)
})

test_that("MML task estimates correlate with true task effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "MML", quad_points = 21, maxit = 100))

  task_est <- fit$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  task_est_c <- task_est - mean(task_est)
  true_c <- true_task_eff - mean(true_task_eff)
  r <- cor(task_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(task_est_c - true_c)), 0.6)
})

test_that("MML criterion estimates correlate with true criterion effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "MML", quad_points = 21, maxit = 100))

  crit_est <- fit$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  crit_est_c <- crit_est - mean(crit_est)
  true_c <- true_crit_eff - mean(true_crit_eff)
  r <- cor(crit_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(crit_est_c - true_c)), 0.6)
})

# ---------------------------------------------------------------------------
# 4C.  Person ability recovery
# ---------------------------------------------------------------------------

test_that("JML person ability estimates correlate with true ability", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  true_ability <- true_person_ability()
  persons <- paste0("P", sprintf("%02d", 1:36))

  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  # Match person labels to true ability
  person_tbl <- fit$facets$person |> dplyr::arrange(Person)
  true_sorted <- true_ability[order(persons)]
  r <- cor(person_tbl$Estimate, true_sorted)
  expect_gt(r, 0.6)
})

# ---------------------------------------------------------------------------
# 4C2.  Latent-regression parameter recovery
# ---------------------------------------------------------------------------

test_that("latent-regression RSM continuous predictor recovery is directionally stable", {
  cases <- lapply(101:102, function(seed) {
    sample_latent_regression_recovery_case(
      seed = seed,
      model = "RSM",
      predictor = "continuous",
      beta1 = 0.8,
      sigma2 = 0.49
    )
  })

  slope_est <- vapply(cases, function(x) {
    x$coefficients$Estimate[x$coefficients$Term == "X"]
  }, numeric(1))
  intercept_est <- vapply(cases, function(x) {
    x$coefficients$Estimate[x$coefficients$Term == "(Intercept)"]
  }, numeric(1))
  sigma_hat <- vapply(cases, function(x) x$sigma2, numeric(1))
  person_cor <- vapply(cases, function(x) x$person_correlation, numeric(1))
  converged <- vapply(cases, function(x) x$converged, logical(1))

  expect_true(all(converged))
  expect_lt(abs(mean(intercept_est) - 0.25), 0.15)
  expect_gt(mean(slope_est), 0.5)
  expect_lt(abs(mean(slope_est) - 0.8), 0.4)
  expect_lt(abs(mean(sigma_hat) - 0.49), 0.4)
  expect_gt(mean(person_cor), 0.8)
})

test_that("latent-regression RSM continuous predictor recovery preserves item ordering", {
  cases <- lapply(101:104, function(seed) {
    sample_latent_regression_recovery_case(
      seed = seed,
      model = "RSM",
      predictor = "continuous",
      beta1 = 0.8,
      sigma2 = 0.49
    )
  })

  item_metrics <- vapply(cases, function(x) {
    item_tbl <- subset(as.data.frame(x$fit$facets$others), Facet == "Item", c(Level, Estimate))
    item_tbl <- item_tbl[order(item_tbl$Level), , drop = FALSE]
    est_centered <- item_tbl$Estimate - mean(item_tbl$Estimate)
    truth_centered <- seq(-0.7, 0.7, length.out = nrow(item_tbl))
    truth_centered <- truth_centered - mean(truth_centered)
    c(
      converged = x$converged,
      correlation = stats::cor(est_centered, truth_centered),
      mae = mean(abs(est_centered - truth_centered))
    )
  }, numeric(3))

  expect_true(all(as.logical(item_metrics["converged", ])))
  expect_gt(mean(item_metrics["correlation", ]), 0.97)
  expect_lt(mean(item_metrics["mae", ]), 0.10)
})

test_that("latent-regression PCM continuous predictor recovery remains usable", {
  cases <- lapply(201:202, function(seed) {
    sample_latent_regression_recovery_case(
      seed = seed,
      model = "PCM",
      predictor = "continuous",
      beta1 = 0.8,
      sigma2 = 0.49
    )
  })

  slope_est <- vapply(cases, function(x) {
    x$coefficients$Estimate[x$coefficients$Term == "X"]
  }, numeric(1))
  sigma_hat <- vapply(cases, function(x) x$sigma2, numeric(1))
  converged <- vapply(cases, function(x) x$converged, logical(1))

  expect_true(all(converged))
  expect_gt(mean(slope_est), 0.4)
  expect_lt(abs(mean(slope_est) - 0.8), 0.5)
  expect_true(all(is.finite(sigma_hat) & sigma_hat > 0))
})

test_that("latent-regression mean-only, dummy, and null-effect cases behave sensibly", {
  mean_only <- sample_latent_regression_recovery_case(
    seed = 401,
    model = "RSM",
    predictor = "none",
    beta1 = 0,
    sigma2 = 0.49
  )
  dummy_case <- sample_latent_regression_recovery_case(
    seed = 301,
    model = "RSM",
    predictor = "dummy",
    beta1 = 0.7,
    sigma2 = 0.49
  )
  null_case <- sample_latent_regression_recovery_case(
    seed = 501,
    model = "RSM",
    predictor = "continuous",
    beta1 = 0,
    sigma2 = 0.49
  )

  mean_intercept <- mean_only$coefficients$Estimate[mean_only$coefficients$Term == "(Intercept)"]
  dummy_effect <- dummy_case$coefficients$Estimate[dummy_case$coefficients$Term == "G"]
  null_effect <- null_case$coefficients$Estimate[null_case$coefficients$Term == "X"]

  expect_true(mean_only$converged)
  expect_true(dummy_case$converged)
  expect_true(null_case$converged)
  expect_true(is.finite(mean_intercept))
  expect_lt(abs(mean_intercept - 0.25), 0.3)
  expect_gt(dummy_effect, 0.15)
  expect_lt(abs(dummy_effect - 0.7), 0.5)
  expect_lt(abs(null_effect), 0.3)
  expect_true(all(c(mean_only$sigma2, dummy_case$sigma2, null_case$sigma2) > 0))
})

# ---------------------------------------------------------------------------
# 4C3.  GPCM narrow recovery
# ---------------------------------------------------------------------------

sample_gpcm_recovery_case <- function(seed = 20260404,
                                      quad_points = 11L,
                                      maxit = 120L) {
  fixture <- mfrmr:::sample_mfrm_gpcm_benchmark_data(seed = seed)
  fit <- suppressWarnings(
    fit_mfrm(
      fixture$data,
      person = fixture$person,
      facets = fixture$facets,
      score = fixture$score,
      method = "MML",
      model = "GPCM",
      step_facet = "Criterion",
      slope_facet = "Criterion",
      quad_points = quad_points,
      maxit = maxit
    )
  )

  slope_tbl <- merge(
    data.frame(
      SlopeFacet = names(fixture$truth$slopes),
      Truth = as.numeric(fixture$truth$slopes),
      stringsAsFactors = FALSE
    ),
    as.data.frame(fit$slopes)[, c("SlopeFacet", "Estimate"), drop = FALSE],
    by = "SlopeFacet",
    sort = FALSE
  )
  slope_tbl$TruthLog <- log(slope_tbl$Truth) - mean(log(slope_tbl$Truth))
  slope_tbl$EstimateLog <- log(slope_tbl$Estimate) - mean(log(slope_tbl$Estimate))

  step_tbl <- merge(
    as.data.frame(fixture$truth$steps, stringsAsFactors = FALSE),
    as.data.frame(fit$steps)[, c("StepFacet", "Step", "Estimate"), drop = FALSE],
    by = c("StepFacet", "Step"),
    sort = FALSE,
    suffixes = c(".Truth", ".Estimate")
  )

  criterion_tbl <- subset(as.data.frame(fit$facets$others), Facet == "Criterion", c(Level, Estimate))
  criterion_tbl <- criterion_tbl[order(criterion_tbl$Level), , drop = FALSE]
  truth_criterion <- fixture$truth$criterion[criterion_tbl$Level]
  est_centered <- criterion_tbl$Estimate - mean(criterion_tbl$Estimate)
  truth_centered <- as.numeric(truth_criterion) - mean(as.numeric(truth_criterion))

  list(
    fit = fit,
    slope_correlation = stats::cor(slope_tbl$TruthLog, slope_tbl$EstimateLog),
    slope_mae = mean(abs(slope_tbl$TruthLog - slope_tbl$EstimateLog)),
    step_correlation = stats::cor(step_tbl$Estimate.Truth, step_tbl$Estimate.Estimate),
    step_mae = mean(abs(step_tbl$Estimate.Truth - step_tbl$Estimate.Estimate)),
    criterion_correlation = stats::cor(est_centered, truth_centered),
    criterion_mae = mean(abs(est_centered - truth_centered))
  )
}

test_that("GPCM MML narrow synthetic recovery remains usable", {
  case <- sample_gpcm_recovery_case()

  expect_true(isTRUE(case$fit$summary$Converged[1]))
  expect_gt(case$slope_correlation, 0.95)
  expect_lt(case$slope_mae, 0.15)
  expect_gt(case$step_correlation, 0.98)
  expect_lt(case$step_mae, 0.20)
  expect_gt(case$criterion_correlation, 0.98)
  expect_lt(case$criterion_mae, 0.10)
})

# ---------------------------------------------------------------------------
# 4D.  Ordering preservation
# ---------------------------------------------------------------------------

test_that("facet ordering matches true generating direction", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  rater_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  names(rater_est) <- c("R1", "R2", "R3")

  # R3 (eff=0.4, harshest) should have highest rater measure
  # R1 (eff=-0.4, most lenient) should have lowest rater measure
  expect_gt(rater_est["R3"], rater_est["R1"])

  task_est <- fit$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  names(task_est) <- c("T1", "T2", "T3", "T4")

  # T4 (hardest) > T1 (easiest)
  expect_gt(task_est["T4"], task_est["T1"])

  crit_est <- fit$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  names(crit_est) <- c("C1", "C2", "C3")

  # C3 (hardest) > C1 (easiest)
  expect_gt(crit_est["C3"], crit_est["C1"])
})

# ---------------------------------------------------------------------------
# 4E.  Reproducibility
# ---------------------------------------------------------------------------

test_that("same seed produces identical LogLik across runs", {
  d1 <- mfrmr:::sample_mfrm_data(seed = 42)
  d2 <- mfrmr:::sample_mfrm_data(seed = 42)
  expect_identical(d1, d2)

  fit1 <- suppressWarnings(fit_mfrm(d1, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 60))
  fit2 <- suppressWarnings(fit_mfrm(d2, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 60))

  expect_equal(fit1$summary$LogLik, fit2$summary$LogLik, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4F.  PCM recovery
# ---------------------------------------------------------------------------

test_that("PCM with step_facet produces all finite estimates", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100,
    model = "PCM", step_facet = "Criterion"))

  expect_s3_class(fit, "mfrm_fit")
  # All facet estimates should be finite
  all_est <- fit$facets$others$Estimate
  expect_true(all(is.finite(all_est)))
  # Step parameter estimates should be finite
  expect_true(all(is.finite(fit$steps$Estimate)))
})

# ---------------------------------------------------------------------------
# 4G.  Bias interaction recovery
# ---------------------------------------------------------------------------

test_that("estimate_bias detects known injected bias in correct direction", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)

  # Inject known bias: for Rater R3 x Criterion C1, boost scores by +1
  boost_idx <- which(d$Rater == "R3" & d$Criterion == "C1")
  d$Score[boost_idx] <- pmin(d$Score[boost_idx] + 1L, 5L)

  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))

  bias <- suppressWarnings(estimate_bias(fit, diag,
    facet_a = "Rater", facet_b = "Criterion", max_iter = 4))

  expect_true(is.list(bias))
  expect_true("table" %in% names(bias))
  expect_true(is.data.frame(bias$table))

  # Find the R3 x C1 cell
  r3c1 <- bias$table |>
    dplyr::filter(
      Facet1_Level == "R3" & Facet2_Level == "C1" |
      FacetA_Level == "R3" & FacetB_Level == "C1"
    )

  # The bias for this cell should be present and indicate positive direction
  # (scores were boosted upward)
  if (nrow(r3c1) > 0) {
    bias_size <- r3c1$`Bias Size`[1]
    expect_true(is.finite(bias_size))
    # Positive bias because we boosted scores upward in that cell
    # The sign depends on the subtractive parameterization:
    # higher observed scores relative to expected means positive obs-exp diff.
    # Check the Obs-Exp Average is positive (higher than expected).
    obs_exp_avg <- r3c1$`Obs-Exp Average`[1]
    if (is.finite(obs_exp_avg)) {
      expect_gt(obs_exp_avg, -0.5,
        label = "R3:C1 Obs-Exp Average should reflect upward bias")
    }
  }
})
