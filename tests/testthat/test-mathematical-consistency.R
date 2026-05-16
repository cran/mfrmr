manual_category_moments <- function(probs) {
  k_vals <- 0:(ncol(probs) - 1L)
  expected <- as.vector(probs %*% k_vals)
  variance <- as.vector(probs %*% (k_vals^2)) - expected^2
  diff <- sweep(
    matrix(k_vals, nrow = nrow(probs), ncol = length(k_vals), byrow = TRUE),
    1L,
    expected,
    FUN = "-"
  )
  list(
    expected = expected,
    variance = variance,
    fourth = as.vector(rowSums(probs * diff^4))
  )
}

aggregate_curve_moments <- function(probabilities) {
  prob_df <- as.data.frame(probabilities, stringsAsFactors = FALSE)
  prob_df$CategoryScore <- suppressWarnings(as.numeric(as.character(prob_df$Category)))
  prob_df$ProbabilitySum <- prob_df$Probability
  prob_df$ExpectedManual <- prob_df$Probability * prob_df$CategoryScore
  prob_df$SecondMomentManual <- prob_df$Probability * prob_df$CategoryScore^2
  prob_df$CategoryInformationSum <- prob_df$CategoryInformation

  out <- stats::aggregate(
    cbind(
      ProbabilitySum,
      ExpectedManual,
      SecondMomentManual,
      CategoryInformationSum
    ) ~ CurveGroup + Theta + Slope + Model,
    data = prob_df,
    FUN = sum
  )
  out$ScoreVarianceManual <- out$SecondMomentManual - out$ExpectedManual^2
  out
}

make_consistency_gpcm_fit <- function(maxit = 8) {
  toy <- load_mfrmr_data("example_core")
  suppressMessages(suppressWarnings(
    fit_mfrm(
      toy, "Person", c("Rater", "Criterion"), "Score",
      method = "JML", model = "GPCM", step_facet = "Criterion", maxit = maxit
    )
  ))
}

make_consistency_pcm_fit <- function(maxit = 10) {
  toy <- load_mfrmr_data("example_core")
  suppressMessages(suppressWarnings(
    fit_mfrm(
      toy, "Person", c("Rater", "Criterion"), "Score",
      method = "JML", model = "PCM", step_facet = "Criterion", maxit = maxit
    )
  ))
}

expect_curve_moment_identities <- function(fit, theta_points = 13) {
  curves <- category_curves_report(fit, theta_points = theta_points, digits = 12)
  moments <- aggregate_curve_moments(curves$probabilities)
  expected <- as.data.frame(curves$expected_ogive, stringsAsFactors = FALSE)
  joined <- merge(
    expected,
    moments,
    by = c("CurveGroup", "Theta", "Slope", "Model"),
    sort = FALSE
  )

  expect_equal(joined$ProbabilitySum, rep(1, nrow(joined)), tolerance = 1e-10)
  expect_equal(joined$ExpectedScore, joined$ExpectedManual, tolerance = 1e-10)
  expect_equal(joined$ScoreVariance, joined$ScoreVarianceManual, tolerance = 1e-10)
  expect_equal(
    joined$Information,
    (joined$Slope^2) * joined$ScoreVariance,
    tolerance = 1e-10
  )
  expect_equal(
    joined$Information,
    joined$CategoryInformationSum,
    tolerance = 1e-10
  )
}

expect_plot_curves_match_report <- function(fit, theta_points = 13) {
  curves <- category_curves_report(fit, theta_points = theta_points, digits = 12)
  ccc <- plot(fit, type = "ccc", draw = FALSE, theta_points = theta_points)
  pathway <- plot(
    fit,
    type = "pathway",
    draw = FALSE,
    theta_points = theta_points,
    include_fit_measures = FALSE
  )

  ccc_joined <- merge(
    as.data.frame(ccc$data$probabilities, stringsAsFactors = FALSE),
    as.data.frame(curves$probabilities, stringsAsFactors = FALSE),
    by = c("CurveGroup", "Theta", "Category", "Slope", "Model"),
    suffixes = c(".plot", ".report"),
    sort = FALSE
  )
  pathway_joined <- merge(
    as.data.frame(pathway$data$expected, stringsAsFactors = FALSE),
    as.data.frame(curves$expected_ogive, stringsAsFactors = FALSE),
    by = c("CurveGroup", "Theta", "Slope", "Model"),
    suffixes = c(".plot", ".report"),
    sort = FALSE
  )

  expect_equal(
    ccc_joined$Probability.plot,
    ccc_joined$Probability.report,
    tolerance = 1e-10
  )
  expect_equal(
    ccc_joined$CategoryInformation.plot,
    ccc_joined$CategoryInformation.report,
    tolerance = 1e-10
  )
  expect_equal(
    pathway_joined$ExpectedScore.plot,
    pathway_joined$ExpectedScore.report,
    tolerance = 1e-10
  )
  expect_equal(
    pathway_joined$Information.plot,
    pathway_joined$Information.report,
    tolerance = 1e-10
  )
}

expect_iif_aggregates_to_tif <- function(fit, theta_points = 13) {
  info <- compute_information(fit, theta_points = theta_points)
  iif_df <- as.data.frame(info$iif, stringsAsFactors = FALSE)
  tif_df <- as.data.frame(info$tif, stringsAsFactors = FALSE)

  for (facet in unique(iif_df$Facet)) {
    facet_iif <- iif_df[iif_df$Facet == facet, , drop = FALSE]
    aggregated <- stats::aggregate(
      Information ~ Theta,
      data = facet_iif,
      FUN = sum
    )
    joined <- merge(tif_df, aggregated, by = "Theta", suffixes = c(".tif", ".iif"))
    expect_equal(
      joined$Information.iif,
      joined$Information.tif,
      tolerance = 1e-10,
      info = paste("facet", facet)
    )
  }
}

test_that("GPCM response-probability bundle obeys moment identities", {
  steps_mat <- matrix(
    c(-0.4, 0.1, 0.6,
      -0.2, 0.5, 0.4,
       0.1, 0.2, 0.7),
    nrow = 3,
    byrow = TRUE
  )
  eta <- c(-1.3, -0.8, -0.2, 0.2, 0.7, 1.1)
  step_idx <- c(1L, 2L, 3L, 1L, 2L, 3L)
  slope_idx <- c(1L, 2L, 3L, 1L, 2L, 3L)
  slopes <- c(0.7, 1.0, 1.4)

  bundle <- mfrmr:::compute_response_probability_bundle(
    config = list(model = "GPCM", n_cat = 4L),
    idx = list(step_idx = step_idx, slope_idx = slope_idx),
    params = list(steps_mat = steps_mat, slopes = slopes),
    eta = eta
  )
  manual <- manual_category_moments(bundle$probs)

  expect_equal(rowSums(bundle$probs), rep(1, length(eta)), tolerance = 1e-12)
  expect_true(all(bundle$probs >= 0))
  expect_equal(bundle$expected_k, manual$expected, tolerance = 1e-12)
  expect_equal(bundle$var_k, manual$variance, tolerance = 1e-12)
  expect_equal(bundle$fourth_central_moment, manual$fourth, tolerance = 1e-12)
  expect_equal(bundle$slope_obs, slopes[slope_idx], tolerance = 1e-12)
  expect_equal(
    bundle$score_information,
    (slopes[slope_idx]^2) * manual$variance,
    tolerance = 1e-12
  )
})

test_that("category-curve reports conserve probabilities, moments, and information", {
  expect_curve_moment_identities(make_toy_fit(maxit = 12, model = "RSM"))
  expect_curve_moment_identities(make_consistency_pcm_fit(), theta_points = 9)
  expect_curve_moment_identities(make_consistency_gpcm_fit(), theta_points = 9)
})

test_that("draw-free category plots expose the same curve data as reports", {
  expect_plot_curves_match_report(make_toy_fit(maxit = 12, model = "RSM"))
  expect_plot_curves_match_report(make_consistency_pcm_fit(), theta_points = 9)
  expect_plot_curves_match_report(make_consistency_gpcm_fit(), theta_points = 9)
})

test_that("design information and conditional SEM obey reciprocal exposure scaling", {
  fit <- make_toy_fit(maxit = 12, model = "RSM")
  fit_dup <- fit
  fit_dup$prep$data <- rbind(fit$prep$data, fit$prep$data)

  info <- compute_information(fit, theta_points = 13)
  info_dup <- compute_information(fit_dup, theta_points = 13)
  sem_plot <- plot_information(info, type = "sem", draw = FALSE)
  positive_info <- info$tif$Information > 0 & is.finite(info$tif$Information)

  expect_equal(
    info$tif$SE[positive_info],
    1 / sqrt(info$tif$Information[positive_info]),
    tolerance = 1e-12
  )
  expect_equal(info$conditional_sem$ConditionalSEM, info$tif$SE)
  expect_equal(info_dup$tif$Information, 2 * info$tif$Information)
  expect_equal(
    info_dup$tif$SE[positive_info],
    info$tif$SE[positive_info] / sqrt(2),
    tolerance = 1e-12
  )
  expect_equal(
    sem_plot$data$conditional_sem$ConditionalSEM,
    info$tif$SE,
    tolerance = 1e-12
  )
})

test_that("facet-level information contributions aggregate to the total curve", {
  expect_iif_aggregates_to_tif(make_toy_fit(maxit = 12, model = "RSM"))
  expect_iif_aggregates_to_tif(make_consistency_pcm_fit(), theta_points = 9)
  expect_iif_aggregates_to_tif(make_consistency_gpcm_fit(), theta_points = 9)
})
test_that("compute_P_geq uses stable R fallback when cpp11 backend option is disabled", {
  probs <- matrix(
    c(0.20, 0.30, 0.50,
      0.10, 0.40, 0.50),
    nrow = 2,
    byrow = TRUE
  )

  old <- getOption("mfrmr.use_cpp11_backend")
  on.exit(options(mfrmr.use_cpp11_backend = old), add = TRUE)
  options(mfrmr.use_cpp11_backend = FALSE)

  expect_equal(
    mfrmr:::compute_P_geq(probs),
    mfrmr:::compute_P_geq_r(probs)
  )
})
