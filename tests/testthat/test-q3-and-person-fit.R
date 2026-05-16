# Tests for the 0.1.6 helpers introduced for local-dependence and
# person-fit reporting: q3_statistic(), compute_person_fit_indices(),
# and mfrm_generalizability().

local({
  .toy <<- load_mfrmr_data("example_core")
  .fit <<- make_toy_fit()
  .diag <<- make_toy_diagnostics(.fit)
})

# --- q3_statistic ---------------------------------------------------------

test_that("q3_statistic returns the documented shape", {
  q3 <- q3_statistic(.fit, diagnostics = .diag)
  expect_s3_class(q3, "mfrm_q3")
  expect_true(all(c("Level1", "Level2", "Q3", "N", "AbsQ3",
                    "YenFlag", "MaraisFlag", "RelativeFlag",
                    "Interpretation") %in% names(q3$pairs)))
  expect_named(q3$thresholds, c("yen", "marais", "relative_offset"))
})

test_that("q3_statistic respects custom thresholds", {
  q3_strict <- q3_statistic(.fit, diagnostics = .diag,
                             yen_threshold = 0.05,
                             marais_threshold = 0.10)
  expect_gte(q3_strict$summary$YenFlagged, 0L)
  q3_lax <- q3_statistic(.fit, diagnostics = .diag,
                         yen_threshold = 0.99,
                         marais_threshold = 0.99)
  expect_equal(q3_lax$summary$YenFlagged, 0L)
})

test_that("q3_statistic rejects unknown facet", {
  expect_error(
    q3_statistic(.fit, diagnostics = .diag, facet = "NotAFacet"),
    "must be one of"
  )
})

# --- compute_person_fit_indices ------------------------------------------

test_that("compute_person_fit_indices returns one row per person", {
  pf <- compute_person_fit_indices(.diag, fit = .fit)
  expect_true(is.data.frame(pf))
  expect_s3_class(pf, "mfrm_person_fit_indices")
  expect_true(all(c("Person", "N", "LogLik", "lz", "lz_star",
                    "lz_star_status", "lz_star_c", "lz_star_variance",
                    "ReportIndex", "ReportValue", "ReportFlagLevel",
                    "ReportFlag", "ReviewStatus", "ReviewReason",
                    "ReportCaveat")
                  %in% names(pf)))
  expect_true(all(pf$lz_star_status %in% c(
    "computed_jml_conditional_calibration", "insufficient_information"
  )))
  expect_true(any(is.finite(pf$lz_star)))
  expect_true(all(pf$ReportIndex %in% c("lz_star", "lz", "none")))
  expect_true(all(pf$ReviewStatus %in% c(
    "review_1pct", "review_5pct", "not_flagged", "not_available"
  )))
  # ECI4 was removed in 0.2.0 — it was a misnamed duplicate of the
  # Outfit ZSTD (linear Smith form), not Tatsuoka & Tatsuoka (1983).
  expect_false("ECI4" %in% names(pf))
  expect_equal(length(unique(pf$Person)), nrow(pf))

  spf <- summary(pf, top_n = 5)
  expect_s3_class(spf, "summary.mfrm_person_fit_indices")
  expect_true(all(c("overview", "status_summary", "report_index_summary",
                    "lz_star_status_summary", "top_review", "caveats",
                    "thresholds", "reporting_map") %in% names(spf)))
  expect_true(is.data.frame(spf$top_review))
  expect_equal(spf$overview$Persons, nrow(pf))
  expect_true(any(spf$report_index_summary$Value == "lz_star"))
  printed <- capture.output(print(spf))
  expect_true(any(grepl("Person-Fit Summary", printed, fixed = TRUE)))
})

test_that("compute_person_fit_indices works without fit (lz_star fit-required)", {
  pf <- compute_person_fit_indices(.diag, fit = NULL)
  expect_true(all(is.na(pf$lz_star)))
  expect_true(all(pf$lz_star_status == "fit_required"))
})

test_that("lz uses true Drasgow polytomous form via PrObserved", {
  # Closed-form check: build a tiny synthetic obs table directly with
  # PrObserved / ItemEntropy / ItemVarLogP populated and verify the
  # output matches the manual computation.
  P_item1 <- c(0.1, 0.2, 0.4, 0.3)
  P_item2 <- c(0.05, 0.15, 0.5, 0.3)
  log_p_item1 <- log(P_item1)
  log_p_item2 <- log(P_item2)

  ent1 <- sum(P_item1 * log_p_item1)
  ent2 <- sum(P_item2 * log_p_item2)
  var1 <- sum(P_item1 * log_p_item1^2) - ent1^2
  var2 <- sum(P_item2 * log_p_item2^2) - ent2^2

  fake_obs <- data.frame(
    Person = c("p1", "p1"),
    Observed = c(2, 3),
    Expected = c(2, 3),
    Residual = c(0, 0),
    PrObserved = c(P_item1[3], P_item2[4]),
    ItemEntropy = c(ent1, ent2),
    ItemVarLogP = c(var1, var2),
    stringsAsFactors = FALSE
  )
  fake_diag <- list(obs = fake_obs)

  pf <- compute_person_fit_indices(fake_diag, fit = NULL)
  expected_loglik <- log(P_item1[3]) + log(P_item2[4])
  expected_e_logp <- ent1 + ent2
  expected_var_logp <- var1 + var2
  expected_lz <- (expected_loglik - expected_e_logp) / sqrt(expected_var_logp)

  expect_equal(pf$LogLik, expected_loglik, tolerance = 1e-12)
  expect_equal(pf$lz, expected_lz, tolerance = 1e-12)
})

test_that("lz_star uses Snijders weight projection for JML-style estimates", {
  P_item1 <- c(0.1, 0.2, 0.4, 0.3)
  P_item2 <- c(0.05, 0.15, 0.5, 0.3)
  slope_item1 <- 0.75
  slope_item2 <- 1.40
  k_vals <- 0:3
  log_p_item1 <- log(P_item1)
  log_p_item2 <- log(P_item2)
  expected1 <- sum(P_item1 * k_vals)
  expected2 <- sum(P_item2 * k_vals)
  r1 <- k_vals - expected1
  r2 <- k_vals - expected2

  ent1 <- sum(P_item1 * log_p_item1)
  ent2 <- sum(P_item2 * log_p_item2)
  var1 <- sum(P_item1 * log_p_item1^2) - ent1^2
  var2 <- sum(P_item2 * log_p_item2^2) - ent2^2
  cov1 <- slope_item1 * sum(P_item1 * log_p_item1 * r1)
  cov2 <- slope_item2 * sum(P_item2 * log_p_item2 * r2)
  info1 <- slope_item1^2 * sum(P_item1 * r1^2)
  info2 <- slope_item2^2 * sum(P_item2 * r2^2)

  fake_obs <- data.frame(
    Person = c("p1", "p1"),
    Observed = c(2, 3),
    Expected = c(expected1, expected2),
    Residual = c(2 - expected1, 3 - expected2),
    PrObserved = c(P_item1[3], P_item2[4]),
    ItemEntropy = c(ent1, ent2),
    ItemVarLogP = c(var1, var2),
    ItemLogPScoreCov = c(cov1, cov2),
    ScoreInformation = c(info1, info2),
    ObservedScoreDerivative = c(slope_item1 * r1[3], slope_item2 * r2[4]),
    stringsAsFactors = FALSE
  )
  fake_diag <- list(obs = fake_obs)
  fake_fit <- structure(
    list(config = list(method = "JMLE"),
         summary = data.frame(Method = "JML")),
    class = "mfrm_fit"
  )

  c_n <- (cov1 + cov2) / (info1 + info2)
  corrected_var <- (var1 + var2) - (cov1 + cov2)^2 / (info1 + info2)
  centered_loglik <- (log(P_item1[3]) - ent1) + (log(P_item2[4]) - ent2)
  score_sum <- slope_item1 * r1[3] + slope_item2 * r2[4]
  expected_lz_star <- (centered_loglik - c_n * score_sum) / sqrt(corrected_var)

  pf <- compute_person_fit_indices(fake_diag, fit = fake_fit)
  expect_equal(pf$lz_star, expected_lz_star, tolerance = 1e-12)
  expect_equal(pf$lz_star_c, c_n, tolerance = 1e-12)
  expect_equal(pf$lz_star_variance, corrected_var, tolerance = 1e-12)
  expect_identical(pf$lz_star_status, "computed_jml_conditional_calibration")
  expect_identical(pf$ReportIndex, "lz_star")
  expect_equal(pf$ReportValue, expected_lz_star, tolerance = 1e-12)
})

test_that("lz_star is not applied to MML/EAP person scores", {
  mml_fit <- make_toy_fit(method = "MML", maxit = 10)
  mml_diag <- make_toy_diagnostics(mml_fit)
  pf <- compute_person_fit_indices(mml_diag, fit = mml_fit)
  expect_true(all(is.na(pf$lz_star)))
  expect_true(all(pf$lz_star_status == "not_applicable_eap"))
  expect_true(all(pf$ReportIndex %in% c("lz", "none")))
  expect_true(all(grepl("not_applicable_eap", pf$ReportCaveat, fixed = TRUE)))
})

test_that("lz_star falls back with explicit status when Snijders information is degenerate", {
  P_item <- c(0.2, 0.5, 0.3)
  log_p <- log(P_item)
  ent <- sum(P_item * log_p)
  var_logp <- sum(P_item * log_p^2) - ent^2
  fake_obs <- data.frame(
    Person = "p1",
    Observed = 1,
    Expected = sum(P_item * 0:2),
    Residual = 0,
    PrObserved = P_item[2],
    ItemEntropy = ent,
    ItemVarLogP = var_logp,
    ItemLogPScoreCov = 0,
    ScoreInformation = 0,
    ObservedScoreDerivative = 0,
    stringsAsFactors = FALSE
  )
  fake_fit <- structure(
    list(config = list(method = "JML"),
         summary = data.frame(Method = "JML")),
    class = "mfrm_fit"
  )
  pf <- compute_person_fit_indices(list(obs = fake_obs), fit = fake_fit)
  expect_true(is.na(pf$lz_star))
  expect_identical(pf$lz_star_status, "insufficient_information")
  expect_identical(pf$ReportIndex, "lz")
  expect_true(grepl("insufficient_information", pf$ReportCaveat, fixed = TRUE))
})

# --- mfrm_generalizability -----------------------------------------------

test_that("mfrm_generalizability returns variance components and G/Phi", {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    skip("lme4 (Suggests) not installed.")
  }
  gt <- mfrm_generalizability(.fit)
  expect_s3_class(gt, "mfrm_generalizability")
  expect_true(all(c("Source", "Variance", "ProportionVariance")
                  %in% names(gt$variance_components)))
  expect_true(all(c("G", "Phi") %in% names(gt$coefficients)))
})

test_that("mfrm_d_study projects G and Phi across planned counts", {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    skip("lme4 (Suggests) not installed.")
  }
  gt <- mfrm_generalizability(.fit)
  ds <- mfrm_d_study(
    gt,
    data.frame(Rater = c(2, 4), Criterion = c(2, 4)),
    residual_scaling = "sensitivity"
  )

  expect_s3_class(ds, "mfrm_d_study")
  expect_true(all(c(
    "Scenario", "n_Rater", "n_Criterion", "ResidualScaling",
    "ResidualDivisor", "G", "Phi", "GStatus", "PhiStatus"
  ) %in% names(ds)))
  expect_equal(nrow(ds), 6L)
  expect_true(all(c("highest_order", "single_condition", "none") %in% ds$ResidualScaling))
  expect_true(all(is.na(ds$G) | (ds$G >= 0 & ds$G <= 1)))
  expect_true(all(is.na(ds$Phi) | (ds$Phi >= 0 & ds$Phi <= 1)))

  highest <- ds[ds$ResidualScaling == "highest_order", , drop = FALSE]
  if (all(is.finite(highest$G))) {
    expect_gte(highest$G[2], highest$G[1] - 1e-8)
  }

  p <- plot(ds, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
  expect_identical(p$data$plot, "coefficients")
  expect_true(all(c("table", "series", "x_var") %in% names(p$data)))
  expect_true(is.list(plot_data(ds)))

  ds_grid <- mfrm_d_study(
    gt,
    expand.grid(Rater = 2:4, Criterion = 2:4),
    residual_scaling = "sensitivity"
  )
  p_panel <- plot(
    ds_grid,
    x_var = "n_Rater",
    group_var = "n_Criterion",
    panel_grid = c("Metric", "ResidualScaling"),
    draw = FALSE
  )
  expect_identical(p_panel$data$group_var, "n_Criterion")
  expect_identical(p_panel$data$panel_grid, c("Metric", "ResidualScaling"))
  expect_true(all(c("MetricFamily", "Series", "Panel") %in% names(p_panel$data$series)))
  expect_true(all(p_panel$data$series$MetricFamily == "G-theory"))

  p_heat <- plot(
    ds_grid,
    type = "heatmap",
    x_var = "n_Rater",
    y_var = "n_Criterion",
    metric = "Phi",
    draw = FALSE
  )
  expect_identical(p_heat$data$plot, "heatmap")
  expect_identical(p_heat$data$metric, "Phi")
  expect_identical(p_heat$data$panel_by, "ResidualScaling")
  expect_s3_class(p_heat$data$surface, "data.frame")

  p_error_heat <- plot(
    ds_grid,
    type = "heatmap",
    x_var = "n_Rater",
    y_var = "n_Criterion",
    metric = "AbsoluteErrorVariance",
    draw = FALSE
  )
  expect_identical(p_error_heat$data$metric, "AbsoluteErrorVariance")

  p_surface <- plot(
    ds_grid,
    type = "surface3d",
    x_var = "n_Rater",
    y_var = "n_Criterion",
    metric = "G",
    panel_by = "ResidualScaling",
    draw = FALSE
  )
  expect_identical(p_surface$data$plot, "surface3d")
  expect_identical(p_surface$data$metric, "G")

  default_ds <- mfrm_d_study(gt)
  expect_equal(nrow(default_ds), 1L)
  expect_identical(default_ds$ResidualScaling, "highest_order")
})
