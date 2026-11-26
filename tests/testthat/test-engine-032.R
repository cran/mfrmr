test_that("check_mfrm_connectivity and pre-estimation fit_mfrm warning detect disconnected designs", {
  ratings <- load_mfrmr_data("example_core")
  expect_true("example_core" %in% load_mfrmr_data("list"))

  conn <- check_mfrm_connectivity(ratings, "Person", c("Rater", "Criterion"), "Score")
  expect_true(conn$connected)
  expect_equal(conn$n_subsets, 1L)
  expect_equal(conn$n_obs, nrow(ratings))

  d_disc <- data.frame(
    Person = c("P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4"),
    Rater = c("R1", "R2", "R1", "R2", "R3", "R4", "R3", "R4"),
    Criterion = c("C1", "C2", "C1", "C2", "C3", "C4", "C3", "C4"),
    Score = c(1, 2, 3, 2, 2, 4, 1, 3),
    stringsAsFactors = FALSE
  )
  conn2 <- check_mfrm_connectivity(d_disc, "Person", c("Rater", "Criterion"), "Score")
  expect_false(conn2$connected)
  expect_equal(conn2$n_subsets, 2L)

  expect_warning(
    fit_mfrm(d_disc, "Person", c("Rater", "Criterion"), "Score"),
    "Disconnected rating design"
  )
  expect_no_warning(
    fit_mfrm(d_disc, "Person", c("Rater", "Criterion"), "Score", warn_connectivity = FALSE)
  )
})

test_that("JML diagnostics retain fitted estimates and all estimators support the standard S3 API", {
  ratings <- load_mfrmr_data("example_core")

  # 1. JML RSM reporting preserves the calibration
  f_ref_jml <- fit_mfrm(
    ratings, "Person", c("Rater", "Criterion"), "Score",
    method = "JML", model = "RSM", rating_min = 1, rating_max = 4
  )
  dg_ref_jml <- diagnose_mfrm(f_ref_jml)
  ft_ref_jml <- fair_average_table(f_ref_jml, diagnostics = dg_ref_jml)
  f_jml <- f_ref_jml
  expect_true(f_jml$converged)
  non_ext <- which(f_jml$person$ExtremeFlag == 0L)
  ref_pm <- dg_ref_jml$measures[dg_ref_jml$measures$Facet == "Person", ]
  ref_pm <- ref_pm[match(f_jml$person$Person, ref_pm$Level), ]
  expect_equal(unname(f_jml$person$Estimate[non_ext]), unname(ref_pm$Estimate[non_ext]), tolerance = 1e-5)
  expect_equal(unname(f_jml$person$ModelSE[non_ext]), unname(ref_pm$ModelSE[non_ext]), tolerance = 1e-5)
  expect_equal(unname(f_jml$person$RealSE[non_ext]), unname(ref_pm$RealSE[non_ext]), tolerance = 1e-5)
  ref_pf <- ft_ref_jml$raw_by_facet$Person
  ref_pf <- ref_pf[match(f_jml$person$Person, ref_pf$Level), ]
  expect_equal(unname(f_jml$person$FairM[non_ext]), unname(ref_pf$FairM[non_ext]), tolerance = 1e-5)

  # 2. Jeffreys-penalized JML + Block-Schur SEs
  f_jef <- fit_mfrm(
    ratings, "Person", c("Rater", "Criterion"), "Score",
    method = "JML_jeffreys", model = "RSM", rating_min = 1, rating_max = 4
  )
  expect_true(f_jef$converged)
  expect_true(all(is.finite(f_jef$person$ModelSE) & f_jef$person$ModelSE > 0))
  expect_true(all(is.finite(f_jef$person$RealSE) & f_jef$person$RealSE >= f_jef$person$ModelSE - 1e-12))

  # 3. MML + Standard S3 API (predict with newdata, coef, vcov, confint, logLik, AIC, BIC, nobs)
  f_mml <- fit_mfrm(
    ratings, "Person", c("Rater", "Criterion"), "Score",
    method = "MML", model = "RSM", rating_min = 1, rating_max = 4,
    quad_points = 31L
  )
  expect_true(f_mml$converged)
  expect_equal(nobs(f_mml), nrow(ratings))
  expect_true(is.finite(as.numeric(logLik(f_mml))))
  expect_true(is.finite(AIC(f_mml)))
  expect_true(is.finite(BIC(f_mml)))

  expect_true(length(coef(f_mml, type = "facets")) > 0L)
  expect_true(length(coef(f_mml, type = "steps")) == 3L)
  V_des <- vcov(f_mml, component = "design", se_type = "sandwich")
  expect_true(is.matrix(V_des) && nrow(V_des) == ncol(V_des))
  ci_sw <- confint(f_mml, component = "facets", se_type = "sandwich")
  expect_equal(ncol(ci_sw), 2L)

  # General predict() on arbitrary newdata conditions
  nd <- expand.grid(
    theta = c(-1.5, 0.0, 1.5),
    Rater = unique(ratings$Rater)[1:2],
    Criterion = unique(ratings$Criterion)[1:2],
    stringsAsFactors = FALSE
  )
  pr_exp <- predict(f_mml, newdata = nd, type = "response", se.fit = TRUE)
  expect_length(pr_exp$fit, nrow(nd))
  expect_true(all(pr_exp$fit >= 1.0 & pr_exp$fit <= 4.0))
  expect_true(all(pr_exp$se.fit >= 0.0))

  pr_mat <- predict(f_mml, newdata = nd, type = "probs")
  expect_equal(dim(pr_mat), c(nrow(nd), 4L))
  expect_equal(unname(rowSums(pr_mat)), rep(1.0, nrow(nd)), tolerance = 1e-10)

  # Score WLE for brand-new examinee responses via predict(type = "wle")
  nd_wle <- data.frame(
    Person = c("NewP1", "NewP1", "NewP2", "NewP2"),
    Rater = rep(unique(ratings$Rater)[1:2], 2L),
    Criterion = rep(unique(ratings$Criterion)[1:2], 2L),
    Score = c(1L, 2L, 4L, 4L),
    stringsAsFactors = FALSE
  )
  wle_new <- predict(f_mml, newdata = nd_wle, type = "wle")
  expect_equal(nrow(wle_new), 2L)
  expect_gt(wle_new$Estimate[wle_new$Person == "NewP2"], wle_new$Estimate[wle_new$Person == "NewP1"])
  labelled <- transform(nd_wle, Score = factor(Score, levels = 4:1))
  expect_equal(predict(f_mml, labelled, type = "wle"), wle_new)
  for (bad_score in c(0, 99)) {
    invalid <- nd_wle
    invalid$Score[2L] <- bad_score
    expect_error(predict(f_mml, invalid, type = "wle"), "row 2.*outside")
  }
  for (bad_score in c(1.25, Inf, NA_real_)) {
    invalid <- nd_wle
    invalid$Score[2L] <- bad_score
    expect_error(predict(f_mml, invalid, type = "wle"), "finite, integer-valued")
  }
  expect_error(predict(f_mml, transform(nd_wle, Person = " "), type = "wle"), "Person ID")
  expect_error(predict(f_mml, nd_wle, type = "eap"), "EAP scoring of new ratings is unavailable")
  expect_equal(nrow(predict(f_mml, type = "eap")), f_mml$n_person)
})

test_that("step_constraints, declarative category_support, and explanatory RSM work with LRT and profile CIs", {
  skip_on_cran()
  ratings <- load_mfrmr_data("example_core")
  crit_lvls <- sort(unique(as.character(ratings$Criterion)))

  # 1. Step constraints across CML, JML_jeffreys, and MML (with estimated population) + profile CI
  decl <- data.frame(Level = crit_lvls[1:2], From = c(1, 1), To = c(2, 2), Group = c("sh1", "sh1"))
  for (m in c("CML", "JML_jeffreys", "MML")) {
    fc <- fit_mfrm(
      ratings, "Person", c("Rater", "Criterion"), "Score",
      model = "PCM", step_facet = "Criterion", method = m,
      population_formula = if (m == "MML") ~ 1 else NULL,
      step_constraints = decl
    )
    expect_true(fc$converged)
    s1 <- fc$steps[fc$steps$StepIndex == 1L & fc$steps$ScaleLevel %in% crit_lvls[1:2], ]
    expect_equal(s1$Estimate[1L], s1$Estimate[2L], tolerance = 1e-12)
    expect_equal(s1$ModelSE[1L], s1$ModelSE[2L], tolerance = 1e-12)
  }

  f_cml_rsm <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
                        model = "RSM", method = "CML")
  f_cml_hyb <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
                        model = "PCM", step_facet = "Criterion", method = "CML",
                        step_constraints = decl)
  f_cml_pcm <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
                        model = "PCM", step_facet = "Criterion", method = "CML")
  an_cml <- anova(f_cml_rsm, f_cml_hyb, f_cml_pcm)
  expect_identical(an_cml$comparison_basis$lrt_status, "computed")
  expect_equal(nrow(an_cml$lrt), 2L)

  ci_prof <- confint(f_cml_rsm, parm = 1L, se_type = "profile")
  expect_equal(dim(ci_prof), c(1L, 2L))
  expect_lt(ci_prof[1L, 1L], coef(f_cml_rsm, type = "free")[1L])
  expect_gt(ci_prof[1L, 2L], coef(f_cml_rsm, type = "free")[1L])

  # 2. Declarative category_support (mixed binary 1..2 and polytomous 1..4 criteria)
  d_mix <- ratings
  is_bin <- d_mix$Criterion == crit_lvls[1L]
  d_mix$Score[is_bin] <- ifelse(d_mix$Score[is_bin] >= 3L, 2L, 1L)
  cs_df <- data.frame(
    Level = crit_lvls,
    Min = rep(1L, length(crit_lvls)),
    Max = c(2L, rep(4L, length(crit_lvls) - 1L)),
    stringsAsFactors = FALSE
  )
  f_mix <- fit_mfrm(
    d_mix, "Person", c("Rater", "Criterion"), "Score",
    model = "PCM", step_facet = "Criterion", method = "MML",
    rating_min = 1, rating_max = 4, category_support = cs_df
  )
  expect_true(f_mix$converged)
  pr_mix <- predict(
    f_mix,
    newdata = data.frame(theta = c(-1, 0, 1), Rater = d_mix$Rater[1L], Criterion = crit_lvls[1L]),
    type = "probs"
  )
  # Categories "3" and "4" are structurally unsupported on binary Criterion 1 -> exact 0!
  expect_equal(unname(pr_mix[, c("3", "4")]), matrix(0.0, nrow = 3L, ncol = 2L))
  expect_equal(unname(rowSums(pr_mix[, c("1", "2")])), rep(1.0, 3L), tolerance = 1e-12)

  # Observed rating outside category_support must fail clearly
  d_bad <- d_mix
  d_bad$Score[which(is_bin)[1L]] <- 3L
  expect_error(
    fit_mfrm(d_bad, "Person", c("Rater", "Criterion"), "Score",
             model = "PCM", step_facet = "Criterion", category_support = cs_df),
    "violates declared `category_support`"
  )
  expect_error(predict(f_mix, d_bad[which(is_bin)[1L], ], type = "wle"), "category support")
  expect_equal(nrow(predict(f_mix, d_mix[which(is_bin)[1L], ], type = "wle")), 1L)

  # 3. Explanatory RSM (facet_formulas + facet_data + unseen level prediction + nested LRT)
  set.seed(42)
  sim_grid <- expand.grid(Person = paste0("P", 1:30), Rater = paste0("R", 1:3), Task = paste0("T", 1:6),
                          stringsAsFactors = FALSE)
  task_cov <- data.frame(
    Task = paste0("T", 1:6),
    Complexity = c(-1.5, -0.9, -0.3, 0.3, 0.9, 1.5),
    stringsAsFactors = FALSE
  )
  th_true <- stats::rnorm(30, 0, 1.0)[match(sim_grid$Person, paste0("P", 1:30))]
  r_true <- c(-0.4, 0.0, 0.4)[match(sim_grid$Rater, paste0("R", 1:3))]
  b_task <- 0.8 * task_cov$Complexity[match(sim_grid$Task, task_cov$Task)]
  eta_s <- th_true - r_true - b_task
  steps_s <- c(0, -0.7, 0.0, 0.7)
  sim_grid$Score <- vapply(eta_s, function(e) {
    lp <- (0:3) * e - cumsum(steps_s)
    sample(1:4, 1L, prob = exp(lp - max(lp)))
  }, integer(1))

  f_exp <- fit_mfrm(
    sim_grid, "Person", c("Rater", "Task"), "Score",
    model = "RSM", method = "MML",
    facet_formulas = list(Task = ~ Complexity),
    facet_data = list(Task = task_cov)
  )
  f_free <- fit_mfrm(
    sim_grid, "Person", c("Rater", "Task"), "Score",
    model = "RSM", method = "MML"
  )
  expect_true(f_exp$converged)
  expect_gt(f_exp$explanatory$Estimate[1L], 0.3)
  cmp_exp <- anova(f_exp, f_free)
  expect_identical(cmp_exp$comparison_basis$lrt_status, "computed")
  expect_equal(cmp_exp$lrt$df, 4L) # (6 - 1) - 1 = 4 df

  # Predict for a brand-new unseen task T_new using its covariate value!
  pr_unseen <- predict(
    f_exp,
    newdata = data.frame(theta = 0.0, Rater = "R2", Task = c("T_easy", "T_hard")),
    facet_data = list(Task = data.frame(Task = c("T_easy", "T_hard"), Complexity = c(-2.0, 2.0))),
    type = "response"
  )
  expect_gt(pr_unseen[1L], pr_unseen[2L])
})

test_that("FACETS and classical RajuSA-Welch degrees of freedom match hand calculations", {
  # Exact numerical verification of FACETS 4.5.0-4.5.3 RajuSA-Welch df (von Davier 2025)
  # v1 = 1 (se1 = 1), v2 = 4 (se2 = 2), n1 = 3, n2 = 10:
  # num = 25; den = 1/4 + 16/11 = 75/44 -> (num/den)*(11/15) = (44/3)*(11/15) = 484/45 = 10.7555556 -> floor = 10
  rw_facets <- mfrmr:::mfrm_cpp_raju_welch_df(1.0, 3.0, 2.0, 10.0, TRUE)
  rw_class <- mfrmr:::mfrm_cpp_raju_welch_df(1.0, 3.0, 2.0, 10.0, FALSE)
  expect_equal(rw_facets$df_unrounded, 484 / 45, tolerance = 1e-10)
  expect_equal(rw_facets$df, 10)
  # Classical Welch-Satterthwaite: 25 / (1/2 + 16/9) = 25 / (41/18) = 450/41 = 10.97561
  expect_equal(rw_class$df, 450 / 41, tolerance = 1e-10)
})

test_that("bias and DFF methods and all plot palettes work together", {
  skip_on_cran()
  ratings <- load_mfrmr_data("example_core")
  fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score", method = "MML")

  b_fac <- estimate_bias(fit, interaction_facets = "all", df_method = "facets")
  b_cls <- estimate_bias(fit, interaction_facets = "all", df_method = "classical")
  b_jnt <- estimate_bias(fit, interaction_facets = c("Rater", "Criterion"), method = "joint")
  expect_gt(nrow(b_fac$pairwise), 0L)
  expect_true(all(b_fac$pairwise$df == floor(b_fac$pairwise$df)))
  expect_false(identical(b_fac$pairwise$df, b_cls$pairwise$df))
  expect_equal(nrow(b_jnt$lrt), 1L)

  # Fast parallel DFF across methods
  for (dm in c("mml_lrt", "cml_lrt", "raju_welch", "mh", "logistic")) {
    dff_res <- analyze_dff(fit, facet = "Rater", group = "Group", data = ratings,
                           method = dm, min_obs = 3L, threads = 1L)
    expect_s3_class(dff_res, "mfrm_dff")
    expect_gt(nrow(dff_res$table), 0L)
  }

  # Unified interrater_agreement()
  ag <- interrater_agreement(fit)
  expect_s3_class(ag, "mfrm_agreement")
  expect_equal(nrow(ag$summary), 7L)
  expect_true(all(is.finite(ag$summary$Estimate)))

  # Unified CVD-safe ggplot2 plot() methods across okabe_ito, viridis, monochrome
  gt <- mfrm_g_study(ratings, scores = "Score", person = "Person",
                     facets = c("Rater", "Criterion"), method = "urgenova")
  info <- compute_information(fit)
  for (pal in c("okabe_ito", "viridis", "monochrome")) {
    expect_s3_class(plot(fit, type = "wright", palette = pal), "ggplot")
    expect_s3_class(plot(fit, type = "pathway", palette = pal), "ggplot")
    expect_s3_class(plot(fit, type = "steps", palette = pal), "ggplot")
    expect_s3_class(plot(gt, type = "dstudy", palette = pal), "ggplot")
    expect_s3_class(plot(b_fac, type = "pairwise", palette = pal), "ggplot")
    expect_s3_class(plot(info, type = "information", palette = pal), "ggplot")
  }
  mono <- mfrmr:::.mfrm_cvd_palette(4L, "monochrome")
  rgb_m <- grDevices::col2rgb(mono)
  expect_equal(rgb_m["red", ], rgb_m["green", ])
  expect_equal(rgb_m["green", ], rgb_m["blue", ])
})

test_that("unified default fit_mfrm(), 1/2-family GPCM, nested/multivariate G-Study, and diagnostic plots work", {
  skip_on_cran()
  ratings <- load_mfrmr_data("example_core")
  d0 <- ratings
  d0$Score <- d0$Score - 1L

  # Default fit_mfrm() routes to C++ mfrm_fit and supports group_anchors & shrinkage
  ga <- data.frame(Facet = "Rater", Level = c("R01", "R02"), Group = "G1", GroupValue = 0.20)
  f_ga <- fit_mfrm(d0, "Person", c("Rater", "Criterion"), "Score", group_anchors = ga)
  expect_s3_class(f_ga, "mfrm_fit")
  r_est <- f_ga$facets$Estimate[f_ga$facets$Facet == "Rater" & f_ga$facets$Level %in% c("R01", "R02")]
  expect_equal(mean(r_est), 0.20, tolerance = 1e-6)

  f_eb <- fit_mfrm(d0, "Person", c("Rater", "Criterion"), "Score",
                   control = list(facet_shrinkage = "empirical_bayes"))
  f_lap <- fit_mfrm(d0, "Person", c("Rater", "Criterion"), "Score",
                    control = list(facet_shrinkage = "laplace", facet_prior_sd = 1.0))
  expect_true(all(is.finite(f_eb$facets$Estimate)))
  expect_true(all(is.finite(f_lap$facets$Estimate)))

  # Downstream S3 summaries on mfrm_fit
  dg <- diagnose_mfrm(f_ga)
  expect_s3_class(summary(dg), "summary.mfrm_diagnostics")
  expect_s3_class(compute_person_fit_indices(dg), "mfrm_person_fit_indices")
  expect_s3_class(summary(estimate_bias(f_ga)), "summary.mfrm_bias")

  # 1-family & 2-family GPCM MML + JML + S3 methods
  f_rsm <- fit_mfrm(d0, "Person", c("Rater", "Criterion"), "Score", model = "RSM", quad_points = 21L)
  f_pcm <- fit_mfrm(d0, "Person", c("Rater", "Criterion"), "Score",
                    model = "PCM", step_facet = "Criterion", quad_points = 21L)
  f_gpcm1 <- fit_mfrm(d0, "Person", c("Rater", "Criterion"), "Score",
                      model = "GPCM", step_facet = "Criterion", slope_facet = "Criterion",
                      quad_points = 31L)
  expect_equal(attr(logLik(f_gpcm1), "df"), 19L)
  expect_equal(round(as.numeric(logLik(f_gpcm1)), 3), -892.610, tolerance = 1e-3)
  expect_equal(exp(mean(log(f_gpcm1$slopes$Estimate))), 1.0, tolerance = 1e-6)

  f_gpcm1_std <- fit_mfrm(d0, "Person", c("Rater", "Criterion"), "Score",
                          model = "GPCM", step_facet = "Criterion", slope_facet = "Criterion",
                          spec = list(gpcm_mml_identification = "fixed_standard_normal"),
                          quad_points = 21L)
  f_gpcm2 <- fit_mfrm(d0, "Person", c("Rater", "Criterion"), "Score",
                      model = "GPCM", step_facet = "Criterion",
                      slope_facet = c("Criterion", "Rater"),
                      spec = list(slope_reference = "Criterion"),
                      quad_points = 21L)
  expect_true(f_gpcm2$converged)
  expect_equal(nrow(f_gpcm2$slopes), 8L)

  an_chain <- anova(f_rsm, f_pcm, f_gpcm1_std, f_gpcm2)
  expect_identical(an_chain$comparison_basis$lrt_status, "computed")
  expect_equal(nrow(an_chain$lrt), 3L)

  expect_equal(length(coef(f_gpcm1, type = "slopes")), 4L)
  expect_equal(dim(vcov(f_gpcm1, component = "slopes")), c(4L, 4L))
  ci_sl <- confint(f_gpcm1, component = "slopes")
  expect_equal(nrow(ci_sl), 4L)
  wle_g <- predict(f_gpcm1, newdata = d0[1:16, ], type = "wle")
  expect_true(all(is.finite(wle_g$Estimate)))
  expect_s3_class(compute_information(f_gpcm1), "mfrm_information")

  f_gpcm_ao <- fit_mfrm(d0, "Person", c("Rater", "Criterion"), "Score",
                        model = "GPCM", step_facet = "Criterion", slope_facet = "Criterion",
                        spec = list(slope_type = "ability_only"),
                        control = list(slope_prior_sd = 0.5), quad_points = 21L)
  expect_true(f_gpcm_ao$converged)

  f_gpcm_jml <- fit_mfrm(d0, "Person", c("Rater", "Criterion"), "Score",
                         model = "GPCM", step_facet = "Criterion", slope_facet = "Criterion",
                         method = "JML", control = list(slope_prior_sd = 0.5))
  expect_true(f_gpcm_jml$converged)

  # Specialized diagnostic ggplot2 plot types
  for (tp in c("category_curves", "person_fit", "residual_qq", "slopes", "qc")) {
    expect_s3_class(plot(f_gpcm1, type = tp, palette = "okabe_ito"), "ggplot")
  }

  # Nested & multivariate weighted fast G-Study + confint()
  dn <- expand.grid(Person = paste0("P", 1:16), Task = paste0("T", 1:3), RInT = 1:2, stringsAsFactors = FALSE)
  dn$Rater <- paste0(dn$Task, "_R", dn$RInT)
  p_num <- as.numeric(factor(dn$Person)) / 8
  t_num <- as.numeric(factor(dn$Task)) / 5
  dn$S1 <- p_num + t_num + 0.2 * dn$RInT + ((seq_len(nrow(dn)) %% 5) - 2) * 0.15
  dn$S2 <- 0.8 * p_num - 0.3 * t_num + ((seq_len(nrow(dn)) %% 7) - 3) * 0.12
  dn$Comp <- 0.6 * dn$S1 + 0.4 * dn$S2

  gt_mv <- mfrm_g_study(dn, scores = c("S1", "S2"), person = "Person",
                        facets = c("Rater", "Task"), nesting = c(Rater = "Task"),
                        weights = c(0.6, 0.4), method = "urgenova")
  gt_uni <- mfrm_g_study(dn, scores = "Comp", person = "Person",
                         facets = c("Rater", "Task"), nesting = c(Rater = "Task"),
                         method = "urgenova")
  expect_equal(gt_mv$variance_components$RawEstimate,
               gt_uni$variance_components$RawEstimate, tolerance = 1e-10)
  ci_gt <- confint(gt_mv)
  expect_true(all(c("G", "Phi") %in% rownames(ci_gt)))
  expect_true(all(is.finite(ci_gt[c("G", "Phi"), ])))

  # Multi-facet unequal-level EB shrinkage + ANOVA G/D-Study ggplot2 methods
  d_op <- load_mfrmr_data("example_operational")
  f_op <- fit_mfrm(d_op, "Person", c("Rater", "Criterion"), "Score")
  f_op_sh <- apply_empirical_bayes_shrinkage(f_op, shrink_person = TRUE)
  expect_equal(nrow(f_op_sh$facets), 9L)
  expect_true(all(is.finite(f_op_sh$facets$ShrunkEstimate)))
  expect_true(all(is.finite(f_op_sh$person$ShrunkEstimate)))


})

test_that("shared-design batch fits retain each replication's ratings", {
  d <- load_mfrmr_data("example_core")
  d2 <- d
  d2$Score <- 5L - d$Score
  facets <- c("Rater", "Criterion")
  scores <- cbind(d$Score, d2$Score)
  shared <- mfrm_batch_fit(d, "Person", facets, "Score", score_matrix = scores)
  separate <- mfrm_batch_fit(list(d, d2), "Person", facets, "Score")

  for (i in seq_along(shared)) {
    copy <- unserialize(serialize(shared[[i]], NULL))
    expect_equal(as.numeric(logLik(copy)), as.numeric(logLik(separate[[i]])))
    expect_equal(copy$encoding$retained_data$Score, scores[, i])
    expect_equal(residuals(copy, type = "response"), residuals(separate[[i]], type = "response"))
  }

  malformed <- scores
  malformed[1L, 2L] <- 1.5
  expect_error(mfrm_batch_fit(d, "Person", facets, "Score", score_matrix = malformed),
               "finite integer-valued ratings", fixed = TRUE)
  malformed[1L, 2L] <- 5L
  expect_error(mfrm_batch_fit(d, "Person", facets, "Score", score_matrix = malformed,
                              rating_min = 1L, rating_max = 4L), "rubric bounds")
  d$Person[1L] <- NA
  expect_error(mfrm_batch_fit(d, "Person", facets, "Score", score_matrix = scores),
               "complete Person and facet identifiers", fixed = TRUE)
})
test_that("factor score labels cannot become internal codes in unified rating inputs", {
  ratings <- load_mfrmr_data("example_core")
  d <- transform(ratings, Score = factor(Score, levels = 4:1))
  expect_error(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score"), "Convert factor labels")
  expect_error(mfrm_analysis(d, "Person", c("Rater", "Criterion"), "Score"), "Convert factor labels")
  expect_error(mfrm_batch_fit(list(d), "Person", c("Rater", "Criterion"), "Score"), "Convert factor labels")
  fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score")
  fit$encoding$retained_data$Score <- d$Score
  expect_true("factor_score_input" %in% summary(fit)$readiness$reason_codes)
  expect_error(confint(fit, component = "facets"), "factor_score_input")
})

test_that("free MML population means are checked jointly with constrained facet locations", {
  d <- load_mfrmr_data("example_core")
  p <- match(d$Person, unique(d$Person))
  r <- match(d$Rater, unique(d$Rater))
  d$Role <- ifelse(d$Rater == "R04", "Teacher", ifelse(p %% 3L + 1L == r, "Self", "Peer"))
  facets <- c("Rater", "Criterion", "Role")
  fit <- fit_mfrm(d, "Person", facets, "Score", method = "MML")
  expect_true(summary(fit)$readiness$ready)
  for (model in c("RSM", "PCM")) {
    expect_error(fit_mfrm(d, "Person", facets, "Score", model = model, step_facet = "Criterion",
                          population_formula = ~ 1, warn_connectivity = FALSE),
                 "population and facet location design is rank-deficient")
  }
  expect_error(mfrm_batch_fit(list(d), "Person", facets, "Score", population_formula = ~ 1),
               "population and facet location design is rank-deficient")
  full_rank <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", population_formula = ~ 1)
  expect_true(summary(full_rank)$readiness$ready)
  fit$encoding$pop_mode <- 1L
  fit$encoding$pop_X_mat <- matrix(1, fit$n_person, 1)
  fit$population_location_rank <- NULL
  expect_true("population_location_rank_deficient" %in% summary(fit)$readiness$reason_codes)
  expect_error(confint(fit, component = "facets"), "population_location_rank_deficient")
})


test_that("live GMFRM estimation supports one and two discrimination facets", {
  ratings <- load_mfrmr_data("example_core")
  for (slopes in list("Criterion", c("Criterion", "Rater"))) {
    fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
      model = "GMFRM", method = "MML", step_facet = "Criterion",
      slope_facet = slopes, spec = list(slope_reference = "Criterion"), quad_points = 31L)
    expect_true(summary(fit)$readiness$ready)
    expect_equal(nrow(fit$slopes), 4L * length(slopes))
    pr <- predict(fit, type = "probs")
    expect_equal(unname(rowSums(pr)), rep(1, nrow(ratings)), tolerance = 1e-10)
    expect_true(all(is.finite(confint(fit, component = "slopes"))))
    reference <- fit$slopes$Estimate[fit$slopes$Facet == "Criterion"]
    expect_equal(exp(mean(log(reference))), 1, tolerance = 1e-6)
  }
})
