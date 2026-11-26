# Unit and regression tests for Non-Bayesian Laplace-EB Markov Rater-Severity Drift MFRM

simulate_drift_panel <- function(seed = 20261122L, J_per_wave = 25L, K = 4L,
                                 with_facets = TRUE, bridge_n = 0L, D = 1.0) {
  set.seed(seed)
  waves <- paste0("W", 1:4)
  raters <- paste0("R", 1:6)
  Tn <- length(waves)
  R <- length(raters)
  J <- J_per_wave * Tn
  persons <- paste0("P", sprintf("%03d", seq_len(J)))
  theta <- numeric(J)
  for (w in seq_len(Tn)) {
    idx_w <- ((w - 1L) * J_per_wave + 1L):(w * J_per_wave)
    th_w <- stats::rnorm(J_per_wave, 0, 1)
    theta[idx_w] <- th_w - mean(th_w)
  }
  names(theta) <- persons

  beta_true <- rbind(
    R1 = c(-0.40, -0.39, -0.41, -0.40),
    R2 = c( 0.10,  0.11,  0.09,  0.10),
    R3 = c(-0.15, -0.14, -0.16, -0.15),
    R4 = c( 0.25,  0.24,  0.26,  0.25),
    R5 = c(-0.65, -0.15,  0.35,  0.85),
    R6 = c( 0.55,  0.55, -0.55, -0.55)
  )
  colnames(beta_true) <- waves

  steps_true <- seq(-0.7, 0.7, length.out = K - 1L)
  steps_true <- steps_true - mean(steps_true)
  crit_true <- c(C1 = -0.45, C2 = 0.05, C3 = 0.40)

  rows <- list()
  idx <- 0L
  for (w in seq_len(Tn)) {
    p_wave <- persons[((w - 1L) * J_per_wave + 1L):(w * J_per_wave)]
    if (bridge_n > 0L && w > 1L) {
      p_wave <- unique(c(persons[seq_len(bridge_n)], p_wave))
    }
    for (p in p_wave) {
      # Assign 4 raters per person-wave in balanced cyclic blocks
      p_num <- match(p, persons)
      r_sel <- ((p_num + w + 0:3) %% R) + 1L
      for (r in r_sel) {
        idx <- idx + 1L
        cr <- names(crit_true)[((idx - 1L) %% length(crit_true)) + 1L]
        f_off <- if (with_facets) crit_true[[cr]] else 0
        eta <- theta[[p]] - beta_true[r, w] - f_off
        g_vec <- c(0, cumsum(steps_true))
        lp <- D * ((0:(K - 1L)) * eta - g_vec)
        pr <- exp(lp - max(lp))
        pr <- pr / sum(pr)
        sc <- sample.int(K, size = 1L, prob = pr)
        rows[[idx]] <- data.frame(
          Person = p,
          Rater = raters[r],
          Wave = waves[w],
          Criterion = cr,
          Score = sc,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(
    data = do.call(rbind, rows),
    theta = theta,
    beta_true = beta_true,
    steps_true = steps_true,
    crit_true = crit_true
  )
}

test_that("analytical gradient and Louis observed Hessian match finite differences across facets, PCM/RSM, and K=2", {
  skip_on_cran()
  sim <- simulate_drift_panel(seed = 101L, J_per_wave = 8L, K = 4L, with_facets = TRUE, D = 1.7)
  sim$data$Task <- ifelse(seq_len(nrow(sim$data)) %% 2L == 0L, "T1", "T2")

  dat <- mfrmr:::.drift_prepare_input(
    sim$data, person = "Person", rater = "Rater", time = "Wave", score = "Score",
    facets = c("Criterion", "Task"), score_levels = 1:4, missing = "fail"
  )
  quad <- mfrmr:::.drift_quadrature_normal(51L, mean = 0, sd = 1)

  # 1) Dynamic + PCM (rater steps) + two fixed facets
  lay_pcm <- mfrmr:::.drift_layout(dat, time_mode = "dynamic", step_mode = "rater")
  set.seed(11)
  xi0 <- stats::rnorm(lay_pcm$n_free, 0, 0.25)
  ev_pcm <- mfrmr:::.drift_mml_eval(xi0, dat, quad, lay_pcm, D = 1.7, need_hessian = TRUE)

  fn_ll <- function(x) mfrmr:::.drift_mml_eval(x, dat, quad, lay_pcm, D = 1.7, need_hessian = FALSE)$ll
  gr_ll <- function(x) mfrmr:::.drift_mml_eval(x, dat, quad, lay_pcm, D = 1.7, need_hessian = FALSE)$grad

  eps <- 1e-6
  central_hess <- function(gr_fun, x) {
    H <- matrix(0, length(x), length(x))
    for (j in seq_along(x)) {
      xp <- x; xm <- x
      xp[j] <- xp[j] + eps; xm[j] <- xm[j] - eps
      H[, j] <- -(gr_fun(xp) - gr_fun(xm)) / (2 * eps)
    }
    H
  }
  num_grad <- vapply(seq_along(xi0), function(i) {
    xp <- xi0; xm <- xi0
    xp[i] <- xp[i] + eps; xm[i] <- xm[i] - eps
    (fn_ll(xp) - fn_ll(xm)) / (2 * eps)
  }, numeric(1))
  num_hess <- central_hess(gr_ll, xi0)

  expect_lt(max(abs(ev_pcm$grad - num_grad)), 1e-7)
  expect_lt(max(abs(ev_pcm$H_mml - num_hess)), 1e-7)

  # 2) Static + RSM (shared steps) + two fixed facets
  lay_rsm <- mfrmr:::.drift_layout(dat, time_mode = "static", step_mode = "shared")
  set.seed(22)
  xi1 <- stats::rnorm(lay_rsm$n_free, 0, 0.25)
  ev_rsm <- mfrmr:::.drift_mml_eval(xi1, dat, quad, lay_rsm, D = 1.0, need_hessian = TRUE)
  fn_rsm <- function(x) mfrmr:::.drift_mml_eval(x, dat, quad, lay_rsm, D = 1.0, need_hessian = FALSE)$ll
  gr_rsm <- function(x) mfrmr:::.drift_mml_eval(x, dat, quad, lay_rsm, D = 1.0, need_hessian = FALSE)$grad
  num_grad_rsm <- vapply(seq_along(xi1), function(i) {
    xp <- xi1; xm <- xi1
    xp[i] <- xp[i] + eps; xm[i] <- xm[i] - eps
    (fn_rsm(xp) - fn_rsm(xm)) / (2 * eps)
  }, numeric(1))
  num_hess_rsm <- central_hess(gr_rsm, xi1)
  expect_lt(max(abs(ev_rsm$grad - num_grad_rsm)), 1e-7)
  expect_lt(max(abs(ev_rsm$H_mml - num_hess_rsm)), 1e-7)

  # 3) Dichotomous scale (K = 2, 0:1)
  bin_df <- sim$data
  bin_df$Score <- as.integer(bin_df$Score >= 3L)
  dat_bin <- mfrmr:::.drift_prepare_input(
    bin_df, person = "Person", rater = "Rater", time = "Wave", score = "Score",
    facets = "Criterion", score_levels = 0:1, missing = "fail"
  )
  lay_bin <- mfrmr:::.drift_layout(dat_bin, time_mode = "dynamic", step_mode = "shared")
  set.seed(33)
  xi2 <- stats::rnorm(lay_bin$n_free, 0, 0.25)
  ev_bin <- mfrmr:::.drift_mml_eval(xi2, dat_bin, quad, lay_bin, D = 1.0, need_hessian = TRUE)
  gr_bin <- function(x) mfrmr:::.drift_mml_eval(x, dat_bin, quad, lay_bin, D = 1.0, need_hessian = FALSE)$grad
  num_hess_bin <- central_hess(gr_bin, xi2)
  expect_lt(max(abs(ev_bin$H_mml - num_hess_bin)), 1e-7)
})

test_that("fit_mfrm_drift recovers heterogeneous drifting raters, fixed facets, and Kass-Steffey SEs", {
  skip_on_cran()
  sim <- simulate_drift_panel(seed = 2026L, J_per_wave = 35L, K = 4L, with_facets = TRUE, bridge_n = 4L, D = 1.7)

  fit_het <- fit_mfrm_drift(
    sim$data,
    person = "Person",
    rater = "Rater",
    time = "Wave",
    score = "Score",
    facets = "Criterion",
    model = "RSM",
    drift_mode = "markov_hetero",
    logistic_scale = 1.7,
    quad_points = 61L
  )
  fit_lap <- fit_mfrm_drift(
    sim$data,
    person = "Person",
    rater = "Rater",
    time = "Wave",
    score = "Score",
    facets = "Criterion",
    model = "RSM",
    drift_mode = "markov_hetero",
    eb_method = "full_laplace",
    logistic_scale = 1.7,
    quad_points = 61L
  )
  fit_stat <- fit_mfrm_drift(
    sim$data,
    person = "Person",
    rater = "Rater",
    time = "Wave",
    score = "Score",
    facets = "Criterion",
    model = "RSM",
    drift_mode = "static",
    logistic_scale = 1.7,
    quad_points = 61L
  )
  fit_homo <- fit_mfrm_drift(
    sim$data,
    person = "Person",
    rater = "Rater",
    time = "Wave",
    score = "Score",
    facets = "Criterion",
    model = "RSM",
    drift_mode = "markov_homo",
    logistic_scale = 1.7,
    quad_points = 61L
  )

  expect_s3_class(fit_het, "mfrm_drift")
  expect_true(fit_het$checks$NumericalReady)
  expect_true(fit_het$checks$InformationPositive)
  expect_true(fit_het$checks$QuadratureStable)
  expect_equal(fit_het$checks$BridgePersons, 4L)
  expect_equal(fit_het$checks$WaveComponents, 1L)
  expect_false(fit_het$checks$StationaryCohortRequired)

  # Heterogeneous drift SD separates drifting raters (R5, R6) from stable raters (R1..R4)
  rd <- fit_het$rater_drift_table
  drifting_sd <- rd$DriftSD_EAP[rd$Rater %in% c("R5", "R6")]
  stable_sd <- rd$DriftSD_EAP[rd$Rater %in% c("R1", "R2", "R3", "R4")]
  expect_gt(min(drifting_sd), max(stable_sd))
  expect_gt(mean(drifting_sd), mean(stable_sd) * 1.5)
  expect_true(all(rd$Flagged[rd$Rater %in% c("R5", "R6")]))

  # Dynamic model beats static model on Laplace marginal posterior & Effective AIC
  expect_gt(fit_het$model_summary$LaplaceLogPost, fit_stat$model_summary$LaplaceLogPost + 8)
  expect_lt(fit_het$model_summary$EffectiveAIC, fit_stat$model_summary$EffectiveAIC - 15)

  # Homogeneous Markov estimates a single shared DriftSD across all raters
  expect_equal(length(unique(round(fit_homo$rater_drift_table$LogSigma, 10))), 1L)

  # EM-Laplace and Full-Laplace agree closely
  expect_gt(stats::cor(coef(fit_het, "severity"), coef(fit_lap, "severity")), 0.999)

  # Kass-Steffey TotalSE >= ConditionalSE everywhere
  expect_true(all(fit_het$severity_table$TotalSE + 1e-12 >= fit_het$severity_table$ConditionalSE))
  expect_equal(fit_het$severity_table$ChangeSE[fit_het$severity_table$WaveIndex == 1L], rep(0, 6L))

  # Sum-to-zero fixed-facet recovery
  ft <- fit_het$facet_table
  expect_equal(sum(ft$Estimate), 0, tolerance = 1e-10)
  expect_gt(stats::cor(ft$Estimate, as.numeric(sim$crit_true[ft$Level])), 0.95)
})

test_that("Markov RW1 exactly interpolates unobserved rater-wave cells and flags disjoint cohorts", {
  skip_on_cran()
  sim <- simulate_drift_panel(seed = 303L, J_per_wave = 18L, K = 3L, with_facets = FALSE, bridge_n = 0L)
  # Drop all Wave W2 ratings for Rater R2
  sub_df <- sim$data[!(sim$data$Rater == "R2" & sim$data$Wave == "W2"), , drop = FALSE]

  fit <- fit_mfrm_drift(
    sub_df,
    person = "Person",
    rater = "Rater",
    time = "Wave",
    score = "Score",
    model = "PCM",
    drift_mode = "markov_hetero",
    quad_points = 61L
  )

  expect_equal(fit$checks$BridgePersons, 0L)
  expect_equal(fit$checks$WaveComponents, 4L)
  expect_true(fit$checks$StationaryCohortRequired)
  expect_equal(fit$checks$UnobservedRaterWaveCells, 1L)

  st_r2 <- fit$severity_table[fit$severity_table$Rater == "R2", , drop = FALSE]
  expect_identical(st_r2$Observed, c(TRUE, FALSE, TRUE, TRUE))
  # Exact GMRF harmonic interpolation: with zero data at t=2, beta_{r,2} = 0.5 * (beta_{r,1} + beta_{r,3})
  expect_equal(st_r2$Estimate[2L], 0.5 * (st_r2$Estimate[1L] + st_r2$Estimate[3L]), tolerance = 1e-8)
  # Interpolated wave has strictly larger TotalSE than observed adjacent waves
  expect_gt(st_r2$TotalSE[2L], min(st_r2$TotalSE[c(1L, 3L)]))
})

test_that("S3 extractors, reference-wave FairM prediction, RDS serialization, and CVD-safe plots obey contracts", {
  skip_on_cran()
  sim <- simulate_drift_panel(seed = 404L, J_per_wave = 15L, K = 4L, with_facets = TRUE, bridge_n = 2L)
  fit <- fit_mfrm_drift(
    sim$data,
    person = "Person",
    rater = "Rater",
    time = "Wave",
    score = "Score",
    facets = "Criterion",
    model = "RSM",
    quad_points = 51L
  )

  # Print & Summary
  out_pr <- capture.output(print(fit))
  expect_match(paste(out_pr, collapse = "\n"), "Markov Rater-Severity Drift MFRM")
  s_obj <- summary(fit)
  expect_s3_class(s_obj, "summary.mfrm_drift")
  out_sm <- capture.output(print(s_obj))
  expect_match(paste(out_sm, collapse = "\n"), "Rater Drift Summary")

  # coef, vcov, confint, as.data.frame
  expect_length(coef(fit, "severity"), 6L * 4L)
  expect_length(coef(fit, "drift_sd"), 6L)
  expect_length(coef(fit, "contrasts"), 6L * 3L)
  expect_length(coef(fit, "steps"), 3L)
  expect_length(coef(fit, "facets"), 3L)
  expect_equal(nrow(vcov(fit, "total")), length(coef(fit, "all")))
  expect_equal(dim(vcov(fit, "hyper")), c(6L, 6L))

  ci_sd <- confint(fit, parm = "drift_sd")
  expect_true(all(ci_sd$Lower > 0 & ci_sd$Upper > ci_sd$Lower))
  ci_con <- confint(fit, parm = "contrasts")
  expect_equal(nrow(ci_con), 18L)
  expect_equal(nrow(as.data.frame(fit, component = "persons")), nrow(fit$person_table))

  # Reference-wave FairM sensitivity & newdata prediction
  f_w1 <- predict(fit, type = "fair_score", ref_time = "W1")
  f_w4 <- predict(fit, type = "fair_score", ref_time = "W4", ref_raters = c("R5", "R6"))
  expect_equal(nrow(f_w1), nrow(fit$person_table))
  expect_false(isTRUE(all.equal(f_w1$FairM, f_w4$FairM)))

  nd <- data.frame(
    Person = fit$person_table$Person[1:3],
    Rater = c("R1", "R5", "R6"),
    Wave = c("W1", "W4", "W4"),
    Criterion = c("C1", "C2", "C3"),
    stringsAsFactors = FALSE
  )
  pred_exp <- predict(fit, newdata = nd, type = "expected")
  pred_prb <- predict(fit, newdata = nd, type = "probabilities")
  expect_equal(nrow(pred_exp), 3L)
  expect_true(all(pred_exp$ExpectedScore >= 1 & pred_exp$ExpectedScore <= 4))
  expect_equal(rowSums(pred_prb), rep(1, 3L), tolerance = 1e-10)

  # Boundary check on diagnose_mfrm
  expect_error(diagnose_mfrm(fit), "dynamic `mfrm_drift` fits")

  # RDS round-trip
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(fit, tf)
  fit_rds <- readRDS(tf)
  expect_equal(coef(fit_rds, "severity"), coef(fit, "severity"))
})


test_that("a small live Markov calibration preserves prediction and saved-result contracts", {
  sim <- simulate_drift_panel(seed = 404L, J_per_wave = 8L, K = 3L,
    with_facets = FALSE, bridge_n = 2L)
  fit <- fit_mfrm_drift(sim$data, person = "Person", rater = "Rater",
    time = "Wave", score = "Score", model = "RSM", quad_points = 31L)
  expect_s3_class(fit, "mfrm_drift")
  expect_equal(fit$checks$BridgePersons, 2L)
  expect_equal(fit$checks$WaveComponents, 1L)
  expect_false(fit$checks$StationaryCohortRequired)
  expect_true(fit$checks$NumericalReady)
  expect_length(coef(fit, "severity"), 24L)
  nd <- sim$data[1:3, ]
  probabilities <- predict(fit, newdata = nd, type = "probabilities")
  expect_equal(unname(rowSums(probabilities)), rep(1, 3L), tolerance = 1e-10)
  expected <- predict(fit, newdata = nd, type = "expected")$ExpectedScore
  expect_equal(expected, as.numeric(probabilities %*% (1:3)), tolerance = 1e-10)
  saved <- unserialize(serialize(fit, NULL))
  expect_identical(coef(saved), coef(fit))
  expect_identical(predict(saved, newdata = nd, type = "probabilities"), probabilities)
})
