diagnostic_response_cases <- function() {
  out <- readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits
  out$GPCM <- readRDS(test_path("fixtures", "native-location-fixed-gpcm.rds"))
  out
}

test_that("diagnostic response paths match Persons by identity", {
  for (f in diagnostic_response_cases()) {
    g <- f
    g$facets$person <- g$facets$person[nrow(g$facets$person):1, ]
    bias <- list(table = data.frame(FacetA = "Rater", FacetB = "Criterion",
      FacetA_Level = as.character(f$prep$data$Rater[1]),
      FacetB_Level = as.character(f$prep$data$Criterion[1]),
      check.names = FALSE, "Bias Size" = .2))
    paths <- list(compute_obs_table, compute_prob_matrix,
      mfrm_residual_table_for_export,
      function(x) mfrm_residual_table_for_export(x, include_probabilities = TRUE),
      calc_expected_category_counts,
      function(x) compute_obs_table_with_bias(x, bias),
      function(x) compute_prob_matrix_with_bias(x, bias),
      function(x) compute_expected_score_vector_from_par(x, x$opt$par))
    for (path in paths) expect_identical(path(g), path(f))
    obs <- compute_obs_table(f)
    # Agreement with independently checked primary extractors in the interior.
    expect_equal(obs$Expected, as.numeric(fitted(f)), tolerance = 1e-12)
    expect_equal(obs$StdResidual, as.numeric(residuals(f, "pearson")), tolerance = 1e-12)
    expect_equal(obs$ConditionalVariance, attr(fitted(f), "diagnostics")$Variance,
      tolerance = 1e-12)
    expect_true(all(obs$PearsonAvailable))
    expect_identical(compute_obs_table_with_bias(f, NULL)$ResponseStatus, obs$ResponseStatus)
  }
})

test_that("missing or ambiguous Person estimates cannot change response assignment", {
  f <- diagnostic_response_cases()[[1]]
  reference <- compute_obs_table(f)
  person <- f$facets$person$Person[1]
  rows <- f$prep$data$Person == person
  for (remove in c(FALSE, TRUE)) {
    g <- f
    if (remove) g$facets$person <- g$facets$person[-1, ] else
      g$facets$person$Estimate[1] <- NA_real_
    obs <- compute_obs_table(g)
    expect_true(all(is.na(obs$Expected[rows])))
    expect_identical(obs$Expected[!rows], reference$Expected[!rows])
    expect_true(all(!obs$PearsonAvailable[rows]))
    expect_true(all(obs$ResponseStatus[rows] == "unavailable"))
    exported <- mfrm_residual_table_for_export(g, include_probabilities = TRUE)
    probabilities <- exported[, startsWith(names(exported), "PrCategory_"), drop = FALSE]
    expect_true(all(is.na(as.matrix(probabilities[rows, , drop = FALSE]))))
    expect_identical(mfrm_residual_table_for_export(g)$ResponseReason, obs$ResponseReason)
  }
  for (id in list(NA_character_, f$facets$person$Person[1])) {
    g <- f; g$facets$person$Person[2] <- id
    expect_error(compute_obs_table(g), "distinct, nonmissing identifiers")
    expect_error(compute_prob_matrix(g), "distinct, nonmissing identifiers")
  }
})

test_that("legacy cutoff preserves true variance and states its diagnostic basis", {
  f <- diagnostic_response_cases()[[1]]
  f$facets$person$Estimate[] <- -30
  obs <- compute_obs_table(f)
  v <- attr(fitted(f), "diagnostics")$Variance
  expect_equal(obs$ConditionalVariance, v)
  expect_true(all(v > 0 & v < 1e-10))
  expect_true(all(is.finite(obs$Expected) & is.finite(obs$Residual)))
  expect_true(all(is.na(obs$Var) & is.na(obs$StdResidual)))
  expect_true(all(obs$ResponseStatus == "small_variance"))
  expect_true(all(grepl("legacy diagnostic cutoff", obs$ResponseReason)))
  expect_true(all(is.finite(residuals(f, "pearson"))))
  expect_true(all(obs$PersonMeasureBasis == "saved EAP"))
  f$facets$person$Estimate[] <- -1000
  zero <- compute_obs_table(f)
  expect_true(all(zero$ConditionalVariance == 0))
  expect_true(all(zero$ResponseStatus == "zero_variance"))
  expect_true(all(grepl("undefined", zero$ResponseReason)))
})

test_that("JML boundary diagnostic proxies are explicit alongside primary estimates", {
  pattern <- rbind(c(0, 0), c(1, 1), c(0, 1), c(1, 0))
  x <- pattern[rep(seq_len(4L), c(2, 2, 4, 6)), ]
  data <- data.frame(Person = rep(seq_len(nrow(x)), each = 2),
    Rater = rep(c("A", "B"), nrow(x)), Score = as.vector(t(x)))
  f <- fit_mfrm(data, "Person", "Rater", "Score", method = "JML",
    maxit = 300, category_policy = "preserve")
  obs <- compute_obs_table(f)
  boundary <- is.infinite(obs$PrimaryPersonMeasure)
  expect_true(any(boundary))
  expect_true(all(is.finite(obs$PersonMeasure[boundary])))
  expect_true(all(obs$PersonMeasureBasis == "JML optimizer"))
  expect_true(all(grepl("differ from the primary", obs$ResponseReason[boundary])))
  expect_true(all(attr(fitted(f), "diagnostics")$Variance[boundary] == 0))
  expect_true(all(obs$ConditionalVariance[boundary] > 0))
  expect_identical(mfrm_residual_table_for_export(f)$PrimaryPersonMeasure, obs$PrimaryPersonMeasure)
  # The public diagnostic bundle retains the row and aggregate explanations.
  diag <- diagnose_mfrm(f, diagnostic_mode = "legacy")
  expect_identical(diag$obs, obs)
  expect_true(all(c("FitStatus", "FitReason", "FitUnavailableRows") %in% names(diag$fit)))
  expect_identical(diag$overall_fit, calc_overall_fit(obs, fit_df_method = "facets"))
})

test_that("aggregate fit never treats unavailable residuals as zero misfit", {
  obs <- data.frame(Rater = c("A", "A", "B", "B"),
    StdSq = c(1, 4, 9, 16), Var = c(2, 1, 3, 2),
    FourthCentralMoment = c(10, 3, 15, 8), Weight = c(1, 2, 1, 3))
  complete <- calc_overall_fit(obs)
  expect_equal(complete$Outfit, (1 + 8 + 9 + 48) / 7)
  expect_equal(complete$FitAvailableRows, 4L)
  for (column in c("StdSq", "Var", "Weight")) {
    bad <- obs; bad[[column]][1] <- NA_real_
    for (method in c("engine", "facets", "both")) {
      overall <- calc_overall_fit(bad, fit_df_method = method)
      facet <- calc_facet_fit(bad, "Rater", fit_df_method = method)
      expect_true(is.na(overall$Infit) && is.na(overall$Outfit))
      expect_true(is.na(overall$DF_Infit) && is.na(overall$DF_Outfit))
      expect_equal(overall$FitUnavailableRows, 1L)
      expect_equal(overall$FitAvailableRows, 3L)
      expect_equal(overall$FitStatus, "unavailable")
      expect_match(overall$FitReason, "no partial-row")
      expect_true(all(is.na(facet$Outfit[facet$Level == "A"])))
      expect_equal(facet$Outfit[facet$Level == "B"], (9 + 48) / 4)
      expect_equal(facet$FitUnavailableRows, c(1L, 0L))
    }
  }
  obs$FourthCentralMoment[1] <- NA_real_
  out <- calc_overall_fit(obs, fit_df_method = "both")
  expect_equal(out$Infit, complete$Infit)
  expect_true(is.na(out$DF_Infit_FACETS) && is.na(out$DF_Outfit_FACETS))
  obs$StdSq[] <- NA_real_
  expect_true(is.na(calc_overall_fit(obs)$Outfit))
  expect_equal(calc_overall_fit(obs)$FitUnavailableRows, 4L)
  empty <- calc_overall_fit(obs[FALSE, ])
  expect_true(is.na(empty$Outfit))
  expect_equal(empty$FitResponseRows, 0L)
})

test_that("public exports preserve explanations and small positive variance", {
  f <- diagnostic_response_cases()[[1]]
  f$facets$person$Estimate[] <- -30
  obs <- compute_obs_table(f)
  d <- structure(list(calculation_version = 2L, obs = obs, facet_names = f$config$facet_names,
    fit_readiness = as.data.frame(mfrmr_get_readiness_record(f)$fit),
    parameter_uncertainty = list(source = mfrm_native_location_identity(f))),
    class = "mfrm_diagnostics")
  columns <- c("PrimaryPersonMeasure", "PersonMeasureBasis", "ConditionalVariance",
    "PearsonAvailable", "ResponseStatus", "ResponseReason")
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  residual <- write_mfrm_residual_file(f, d, path, digits = 0)
  score <- facets_output_file_bundle(f, d, include = "score", digits = 0,
    score_se_method = "none")
  for (table in list(residual$table, score$scorefile, read.csv(path))) {
    expect_true(all(columns %in% names(table)))
    expect_true(all(table$ConditionalVariance > 0 & table$ConditionalVariance < 1e-10))
    expect_equal(table$ConditionalVariance / obs$ConditionalVariance,
      rep(1, nrow(obs)), tolerance = 1e-12)
    expect_identical(table$ResponseStatus, obs$ResponseStatus)
    expect_identical(table$ResponseReason, obs$ResponseReason)
  }
})
