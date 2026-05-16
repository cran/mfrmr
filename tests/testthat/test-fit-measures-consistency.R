manual_fit_band <- function(value, lower, upper, zstd_cut, type) {
  value <- suppressWarnings(as.numeric(value))
  out <- rep("not_available", length(value))
  ok <- is.finite(value)
  if (identical(type, "zstd")) {
    out[ok & value > zstd_cut] <- "underfit"
    out[ok & value < -zstd_cut] <- "overfit"
    out[ok & value <= zstd_cut & value >= -zstd_cut] <- "within_band"
  } else {
    out[ok & value > upper] <- "underfit"
    out[ok & value < lower] <- "overfit"
    out[ok & value <= upper & value >= lower] <- "within_band"
  }
  out
}

manual_fit_status <- function(tbl, lower, upper, zstd_cut) {
  infit_band <- manual_fit_band(tbl$Infit, lower, upper, zstd_cut, "mnsq")
  outfit_band <- manual_fit_band(tbl$Outfit, lower, upper, zstd_cut, "mnsq")
  infit_z_band <- manual_fit_band(tbl$InfitZSTD, lower, upper, zstd_cut, "zstd")
  outfit_z_band <- manual_fit_band(tbl$OutfitZSTD, lower, upper, zstd_cut, "zstd")
  underfit <- infit_band == "underfit" | outfit_band == "underfit" |
    infit_z_band == "underfit" | outfit_z_band == "underfit"
  overfit <- infit_band == "overfit" | outfit_band == "overfit" |
    infit_z_band == "overfit" | outfit_z_band == "overfit"
  available <- is.finite(tbl$Infit) | is.finite(tbl$Outfit) |
    is.finite(tbl$InfitZSTD) | is.finite(tbl$OutfitZSTD)
  ifelse(
    !available,
    "not_available",
    ifelse(underfit & overfit, "mixed",
           ifelse(underfit, "underfit",
                  ifelse(overfit, "overfit", "within_band")))
  )
}

manual_profile_counts <- function(tbl, facet_value, lower, upper, zstd_cut) {
  status <- manual_fit_status(tbl, lower, upper, zstd_cut)
  available <- status != "not_available"
  denom <- sum(available, na.rm = TRUE)
  safe_rate <- function(n) if (denom > 0L) n / denom else NA_real_
  under_n <- sum(status == "underfit", na.rm = TRUE)
  over_n <- sum(status == "overfit", na.rm = TRUE)
  mixed_n <- sum(status == "mixed", na.rm = TRUE)
  within_n <- sum(status == "within_band", na.rm = TRUE)
  not_available_n <- sum(status == "not_available", na.rm = TRUE)
  data.frame(
    Facet = facet_value,
    Rows = nrow(tbl),
    AvailableRows = denom,
    UnderfitRows = under_n,
    OverfitRows = over_n,
    MixedRows = mixed_n,
    WithinBandRows = within_n,
    NotAvailableRows = not_available_n,
    UnderfitRate = safe_rate(under_n),
    OverfitRate = safe_rate(over_n),
    MixedRate = safe_rate(mixed_n),
    AnyFlagRate = safe_rate(under_n + over_n + mixed_n),
    stringsAsFactors = FALSE
  )
}

test_that("fit_measures_table preserves df/ZSTD and CI formulas", {
  fit <- make_toy_fit(maxit = 12)
  diag <- diagnose_mfrm(fit, residual_pca = "none", fit_df_method = "both")
  fm <- fit_measures_table(
    fit,
    diagnostics = diag,
    include_person = FALSE,
    fit_df_method = "both",
    threshold_profiles = "all",
    ci_level = 0.90,
    top_n = Inf
  )
  tbl <- as.data.frame(fm$table, stringsAsFactors = FALSE)
  z_ci <- stats::qnorm(1 - (1 - 0.90) / 2)
  ci_ok <- is.finite(tbl$Measure) & is.finite(tbl$SE) & tbl$SE >= 0

  expect_equal(
    tbl$InfitZSTD_ENGINE,
    mfrmr:::zstd_from_mnsq(tbl$Infit, tbl$DF_Infit_ENGINE),
    tolerance = 1e-12
  )
  expect_equal(
    tbl$OutfitZSTD_ENGINE,
    mfrmr:::zstd_from_mnsq(tbl$Outfit, tbl$DF_Outfit_ENGINE),
    tolerance = 1e-12
  )
  expect_equal(
    tbl$InfitZSTD_FACETS,
    mfrmr:::zstd_from_mnsq_facets(tbl$Infit, tbl$DF_Infit_FACETS, cap = 9),
    tolerance = 1e-12
  )
  expect_equal(
    tbl$OutfitZSTD_FACETS,
    mfrmr:::zstd_from_mnsq_facets(tbl$Outfit, tbl$DF_Outfit_FACETS, cap = 9),
    tolerance = 1e-12
  )
  expect_equal(
    tbl$CI_Lower[ci_ok],
    tbl$Measure[ci_ok] - z_ci * tbl$SE[ci_ok],
    tolerance = 1e-12
  )
  expect_equal(
    tbl$CI_Upper[ci_ok],
    tbl$Measure[ci_ok] + z_ci * tbl$SE[ci_ok],
    tolerance = 1e-12
  )

  fm_facets <- fit_measures_table(
    fit,
    diagnostics = diagnose_mfrm(fit, residual_pca = "none"),
    include_person = FALSE,
    fit_df_method = "facets",
    top_n = Inf
  )
  facets_tbl <- as.data.frame(fm_facets$table, stringsAsFactors = FALSE)
  expect_true(all(facets_tbl$FitDfMethod == "facets_wright_masters"))
  expect_equal(facets_tbl$DF_Infit, facets_tbl$DF_Infit_FACETS, tolerance = 1e-12)
  expect_equal(facets_tbl$DF_Outfit, facets_tbl$DF_Outfit_FACETS, tolerance = 1e-12)
  expect_equal(facets_tbl$InfitZSTD, facets_tbl$InfitZSTD_FACETS, tolerance = 1e-12)
  expect_equal(facets_tbl$OutfitZSTD, facets_tbl$OutfitZSTD_FACETS, tolerance = 1e-12)
})

test_that("fit status, subsets, and threshold-profile rates are self-consistent", {
  fit <- make_toy_fit(maxit = 12)
  diag <- diagnose_mfrm(fit, residual_pca = "none", fit_df_method = "both")
  fm <- fit_measures_table(
    fit,
    diagnostics = diag,
    include_person = FALSE,
    fit_df_method = "both",
    threshold_profiles = "all",
    top_n = Inf
  )
  tbl <- as.data.frame(fm$table, stringsAsFactors = FALSE)
  lower <- fm$settings$lower
  upper <- fm$settings$upper
  zstd_cut <- fm$settings$zstd_cut
  manual_status <- manual_fit_status(tbl, lower, upper, zstd_cut)

  expect_equal(tbl$FitStatus, manual_status)
  expect_equal(fm$summary$Rows, nrow(tbl))
  expect_equal(fm$summary$UnderfitRows, sum(tbl$FitStatus == "underfit"))
  expect_equal(fm$summary$OverfitRows, sum(tbl$FitStatus == "overfit"))
  expect_equal(fm$summary$MixedRows, sum(tbl$FitStatus == "mixed"))
  expect_equal(fm$summary$WithinBandRows, sum(tbl$FitStatus == "within_band"))
  expect_equal(fm$summary$NotAvailableRows, sum(tbl$FitStatus == "not_available"))
  expect_equal(nrow(fm$underfit), sum(tbl$FitStatus == "underfit"))
  expect_equal(nrow(fm$overfit), sum(tbl$FitStatus == "overfit"))
  expect_equal(nrow(fm$mixed), sum(tbl$FitStatus == "mixed"))

  profile_tbl <- as.data.frame(fm$threshold_profiles, stringsAsFactors = FALSE)
  manual_profiles <- do.call(rbind, lapply(seq_len(nrow(profile_tbl)), function(i) {
    profile <- profile_tbl[i, , drop = FALSE]
    facets <- sort(unique(as.character(tbl$Facet)))
    counts <- rbind(
      manual_profile_counts(tbl, "All facets", profile$Lower, profile$Upper, zstd_cut),
      do.call(rbind, lapply(facets, function(facet) {
        manual_profile_counts(
          tbl[tbl$Facet == facet, , drop = FALSE],
          facet,
          profile$Lower,
          profile$Upper,
          zstd_cut
        )
      }))
    )
    cbind(
      profile[rep(1L, nrow(counts)), , drop = FALSE],
      ZSTDCut = zstd_cut,
      counts,
      row.names = NULL
    )
  }))
  profile_joined <- merge(
    as.data.frame(fm$profile_summary, stringsAsFactors = FALSE),
    manual_profiles,
    by = c("Profile", "Facet"),
    suffixes = c(".reported", ".manual"),
    sort = FALSE
  )
  count_cols <- c(
    "Rows", "AvailableRows", "UnderfitRows", "OverfitRows",
    "MixedRows", "WithinBandRows", "NotAvailableRows"
  )
  rate_cols <- c("UnderfitRate", "OverfitRate", "MixedRate", "AnyFlagRate")
  for (col in count_cols) {
    expect_equal(profile_joined[[paste0(col, ".reported")]],
                 profile_joined[[paste0(col, ".manual")]])
  }
  for (col in rate_cols) {
    expect_equal(profile_joined[[paste0(col, ".reported")]],
                 profile_joined[[paste0(col, ".manual")]],
                 tolerance = 1e-12)
  }
})

test_that("df-sensitivity classifications and summaries follow documented rules", {
  synth <- data.frame(
    Facet = "Rater",
    Level = c("flag", "large", "df", "small", "same", "missing"),
    Infit = c(1.1, 1.1, 1.1, 1.1, 1.1, NA),
    Outfit = c(1.1, 1.1, 1.1, 1.1, 1.1, NA),
    DF_Infit_ENGINE = c(100, 100, 100, 100, 100, NA),
    DF_Infit_FACETS = c(100, 100, 80, 100, 100, NA),
    DF_Outfit_ENGINE = c(100, 100, 100, 100, 100, NA),
    DF_Outfit_FACETS = c(100, 100, 80, 100, 100, NA),
    InfitZSTD_ENGINE = c(1.90, 0.00, 0.20, 0.10, 0.10, NA),
    InfitZSTD_FACETS = c(2.10, 0.70, 0.23, 0.20, 0.11, NA),
    OutfitZSTD_ENGINE = c(0.00, 0.00, 0.20, 0.10, 0.10, NA),
    OutfitZSTD_FACETS = c(0.00, 0.70, 0.23, 0.20, 0.11, NA),
    stringsAsFactors = FALSE
  )
  sens <- mfrmr:::build_fit_measure_df_sensitivity(
    synth,
    zstd_cut = 2,
    df_zstd_tolerance = 0.05,
    df_zstd_large_shift = 0.5,
    df_ratio_tolerance = 0.05
  )
  status <- stats::setNames(sens$DfSensitivityStatus, sens$Level)

  expect_equal(unname(status["flag"]), "flag_changed_by_df")
  expect_equal(unname(status["large"]), "large_zstd_shift")
  expect_equal(unname(status["df"]), "df_convention_difference")
  expect_equal(unname(status["small"]), "small_zstd_shift")
  expect_equal(unname(status["same"]), "same_or_rounding")
  expect_equal(unname(status["missing"]), "not_available")
  expect_true(sens$FlagChangedByDf[sens$Level == "flag"])

  summary_tbl <- mfrmr:::summarize_fit_measure_df_sensitivity(sens)
  expect_equal(summary_tbl$Rows, nrow(sens))
  expect_equal(summary_tbl$ComparedRows, 5L)
  expect_equal(summary_tbl$FlagChangedByDfRows, 1L)
  expect_equal(summary_tbl$LargeZSTDShiftRows, 1L)
  expect_equal(summary_tbl$DfConventionDifferenceRows, 1L)
  expect_equal(summary_tbl$SmallZSTDShiftRows, 1L)
  expect_equal(summary_tbl$SameOrRoundingRows, 1L)
})

test_that("fit-measure plot data recomputes CI and exposes df-sensitivity thresholds", {
  fit <- make_toy_fit(maxit = 12)
  fm <- fit_measures_table(
    fit,
    include_person = FALSE,
    fit_df_method = "both",
    df_zstd_tolerance = 0.01,
    df_zstd_large_shift = 0.25,
    top_n = Inf
  )
  ci_plot <- plot(fm, type = "measure_ci", ci_level = 0.80, draw = FALSE)
  ci_tbl <- as.data.frame(ci_plot$data$table, stringsAsFactors = FALSE)
  z_ci <- stats::qnorm(1 - (1 - 0.80) / 2)
  ci_ok <- is.finite(ci_tbl$Measure) & is.finite(ci_tbl$SE) & ci_tbl$SE >= 0

  expect_equal(ci_plot$data$ci_level, 0.80)
  expect_equal(ci_tbl$CI_Level, rep(0.80, nrow(ci_tbl)))
  expect_equal(
    ci_tbl$CI_Lower[ci_ok],
    ci_tbl$Measure[ci_ok] - z_ci * ci_tbl$SE[ci_ok],
    tolerance = 1e-12
  )
  expect_equal(
    ci_tbl$CI_Upper[ci_ok],
    ci_tbl$Measure[ci_ok] + z_ci * ci_tbl$SE[ci_ok],
    tolerance = 1e-12
  )

  df_plot <- plot(fm, type = "df_sensitivity", top_n = 5, draw = FALSE)
  df_tbl <- as.data.frame(df_plot$data$table, stringsAsFactors = FALSE)
  expect_equal(unname(df_plot$data$thresholds["tolerance"]), fm$settings$df_zstd_tolerance)
  expect_equal(unname(df_plot$data$thresholds["large_shift"]), fm$settings$df_zstd_large_shift)
  expect_lte(nrow(df_tbl), 5L)
  expect_true(all(is.finite(df_tbl$MaxAbsZSTDDiff_FACETS_vs_ENGINE)))
  expect_equal(
    df_tbl$MaxAbsZSTDDiff_FACETS_vs_ENGINE,
    sort(df_tbl$MaxAbsZSTDDiff_FACETS_vs_ENGINE, decreasing = TRUE),
    tolerance = 1e-12
  )
})
