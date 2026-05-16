test_that("separation, reliability, and strata match hand calculations", {
  measures <- tibble::tibble(
    Facet = rep("Rater", 3),
    Level = paste0("R", 1:3),
    Estimate = c(-1, 0, 1),
    ModelSE = rep(0.2, 3),
    RealSE = rep(0.3, 3),
    Infit = c(1.0, 1.4, 0.8),
    Outfit = c(1.1, 1.2, 0.9),
    PrecisionTier = rep("model_based", 3),
    Converged = rep(TRUE, 3)
  )

  out <- mfrmr:::calc_reliability(measures)
  row <- out[out$Facet == "Rater", , drop = FALSE]

  observed_var <- stats::var(measures$Estimate)
  model_error_var <- mean(measures$ModelSE^2)
  model_true_var <- observed_var - model_error_var
  model_rmse <- sqrt(model_error_var)
  model_true_sd <- sqrt(model_true_var)
  model_separation <- model_true_sd / model_rmse
  model_reliability <- model_true_var / observed_var
  model_strata <- (4 * model_separation + 1) / 3

  real_error_var <- mean(measures$RealSE^2)
  real_true_var <- observed_var - real_error_var
  real_rmse <- sqrt(real_error_var)
  real_true_sd <- sqrt(real_true_var)
  real_separation <- real_true_sd / real_rmse
  real_reliability <- real_true_var / observed_var
  real_strata <- (4 * real_separation + 1) / 3

  expect_equal(row$ObservedVariance, observed_var)
  expect_equal(row$ModelErrorVariance, model_error_var)
  expect_equal(row$ModelTrueVariance, model_true_var)
  expect_equal(row$Separation, model_separation)
  expect_equal(row$Reliability, model_reliability)
  expect_equal(row$Strata, model_strata)
  expect_equal(row$ModelRMSE, model_rmse)
  expect_equal(row$ModelTrueSD, model_true_sd)
  expect_equal(row$RealErrorVariance, real_error_var)
  expect_equal(row$RealTrueVariance, real_true_var)
  expect_equal(row$RealSeparation, real_separation)
  expect_equal(row$RealReliability, real_reliability)
  expect_equal(row$RealStrata, real_strata)
  expect_equal(row$MeanInfit, mean(measures$Infit))
  expect_equal(row$MeanOutfit, mean(measures$Outfit))
  expect_identical(row$ReliabilityUse, "primary_reporting")
})

test_that("facet precision summary includes sample and population bases", {
  measures <- tibble::tibble(
    Facet = rep("Criterion", 3),
    Level = paste0("C", 1:3),
    Estimate = c(-1, 0, 1),
    ModelSE = rep(0.2, 3),
    RealSE = rep(0.3, 3),
    Infit = c(1.0, 1.1, 0.9),
    Outfit = c(1.2, 1.0, 0.8)
  )

  out <- mfrmr:::build_facet_precision_summary(measures)
  sample_model <- out[out$Facet == "Criterion" &
                        out$DistributionBasis == "sample" &
                        out$SEMode == "model", , drop = FALSE]
  population_model <- out[out$Facet == "Criterion" &
                            out$DistributionBasis == "population" &
                            out$SEMode == "model", , drop = FALSE]

  expect_equal(nrow(out), 4L)
  expect_equal(sample_model$ObservedVariance, stats::var(measures$Estimate))
  expect_equal(
    population_model$ObservedVariance,
    mean((measures$Estimate - mean(measures$Estimate))^2)
  )
  expect_lt(population_model$Reliability, sample_model$Reliability)
})

test_that("overall and facet Infit/Outfit match weighted hand calculations", {
  obs <- tibble::tibble(
    Person = paste0("P", 1:4),
    Rater = c("R1", "R1", "R2", "R2"),
    StdSq = c(1, 4, 9, 16),
    Var = c(2, 1, 3, 2),
    Weight = c(1, 2, 1, 3)
  )

  overall <- mfrmr:::calc_overall_fit(obs)
  expected_infit <- sum(obs$StdSq * obs$Var * obs$Weight) /
    sum(obs$Var * obs$Weight)
  expected_outfit <- sum(obs$StdSq * obs$Weight) / sum(obs$Weight)

  expect_equal(overall$Infit, expected_infit)
  expect_equal(overall$Outfit, expected_outfit)
  expect_equal(overall$DF_Infit, sum(obs$Var * obs$Weight))
  expect_equal(overall$DF_Outfit, sum(obs$Weight))

  by_facet <- mfrmr:::calc_facet_fit(obs, facet_cols = "Rater")
  r1 <- by_facet[by_facet$Level == "R1", , drop = FALSE]
  r1_idx <- obs$Rater == "R1"
  expect_equal(
    r1$Infit,
    sum(obs$StdSq[r1_idx] * obs$Var[r1_idx] * obs$Weight[r1_idx]) /
      sum(obs$Var[r1_idx] * obs$Weight[r1_idx])
  )
  expect_equal(
    r1$Outfit,
    sum(obs$StdSq[r1_idx] * obs$Weight[r1_idx]) /
      sum(obs$Weight[r1_idx])
  )
})

test_that("FACETS-style fit df uses fourth-moment Wright-Masters formula", {
  obs <- tibble::tibble(
    Person = paste0("P", 1:4),
    Rater = c("R1", "R1", "R2", "R2"),
    StdSq = c(1, 4, 9, 16),
    Var = c(2, 1, 3, 2),
    FourthCentralMoment = c(10, 3, 15, 8),
    Weight = c(1, 2, 1, 3)
  )

  sum_w <- sum(obs$Weight)
  sum_var_w <- sum(obs$Var * obs$Weight)
  denom_infit <- sum(obs$Weight * (obs$FourthCentralMoment - obs$Var^2))
  denom_outfit <- sum(obs$Weight * (obs$FourthCentralMoment / obs$Var^2 - 1))
  expected_df_infit <- 2 * sum_var_w^2 / denom_infit
  expected_df_outfit <- 2 * sum_w^2 / denom_outfit

  overall <- mfrmr:::calc_overall_fit(obs, fit_df_method = "both")
  expect_equal(overall$DF_Infit, sum_var_w)
  expect_equal(overall$DF_Outfit, sum_w)
  expect_equal(overall$FitDfMethod, "engine_primary_facets_available")
  expect_equal(overall$FitZSTDCap, 9)
  expect_equal(overall$DF_Infit_FACETS, expected_df_infit)
  expect_equal(overall$DF_Outfit_FACETS, expected_df_outfit)
  expect_equal(
    overall$InfitZSTD_FACETS,
    mfrmr:::zstd_from_mnsq_facets(overall$Infit, expected_df_infit, cap = 9)
  )
  expect_equal(
    overall$OutfitZSTD_FACETS,
    mfrmr:::zstd_from_mnsq_facets(overall$Outfit, expected_df_outfit, cap = 9)
  )

  facets_primary <- mfrmr:::calc_overall_fit(obs, fit_df_method = "facets")
  expect_equal(facets_primary$DF_Infit, expected_df_infit)
  expect_equal(facets_primary$DF_Outfit, expected_df_outfit)
  expect_equal(facets_primary$FitDfMethod, "facets_wright_masters")
  expect_equal(facets_primary$FitZSTDCap, 9)
  expect_equal(facets_primary$InfitZSTD, overall$InfitZSTD_FACETS)
  expect_equal(facets_primary$OutfitZSTD, overall$OutfitZSTD_FACETS)
  expect_true(all(c("DF_Infit_ENGINE", "DF_Outfit_ENGINE") %in%
                    names(facets_primary)))

  engine_primary <- mfrmr:::calc_overall_fit(obs)
  expect_false(any(grepl("_FACETS$|FitDfMethod|FitZSTD", names(engine_primary))))

  by_facet <- mfrmr:::calc_facet_fit(obs, facet_cols = "Rater",
                                     fit_df_method = "both")
  r1 <- by_facet[by_facet$Level == "R1", , drop = FALSE]
  r1_idx <- obs$Rater == "R1"
  r1_sum_w <- sum(obs$Weight[r1_idx])
  r1_sum_var_w <- sum(obs$Var[r1_idx] * obs$Weight[r1_idx])
  r1_denom_infit <- sum(obs$Weight[r1_idx] *
                          (obs$FourthCentralMoment[r1_idx] - obs$Var[r1_idx]^2))
  r1_denom_outfit <- sum(obs$Weight[r1_idx] *
                           (obs$FourthCentralMoment[r1_idx] / obs$Var[r1_idx]^2 - 1))
  expect_equal(r1$DF_Infit_FACETS, 2 * r1_sum_var_w^2 / r1_denom_infit)
  expect_equal(r1$DF_Outfit_FACETS, 2 * r1_sum_w^2 / r1_denom_outfit)
})
