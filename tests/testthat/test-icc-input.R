skip_if_not_installed("lme4")

test_that("ICC ratios and availability survive score-unit changes without rounding away variances", {
  data <- load_mfrmr_data("example_core")
  original <- compute_facet_icc(data, c("Rater", "Criterion"), "Score", "Person")
  for (scale in c(1e-4, 1e4)) {
    changed <- data
    changed$Score <- data$Score * scale
    out <- compute_facet_icc(changed, c("Rater", "Criterion"), "Score", "Person")
    expect_true(all(is.finite(out$ICC)))
    expect_equal(out$ICC, original$ICC)
    expect_equal(out$Variance / scale^2, original$Variance, tolerance = 1e-5)
    expect_identical(out$Interpretation, original$Interpretation)
    deff <- compute_facet_design_effect(changed, c("Rater", "Criterion"), out)
    expect_identical(deff, compute_facet_design_effect(data, c("Rater", "Criterion"), original))
    if (scale < 1) {
      expect_true(all(out$Variance > 0 & out$Variance < 1e-6))
      expect_output(print(out), "[1-9][0-9.]*e-[0-9]+")
      expect_output(summary(out), "[1-9][0-9.]*e-[0-9]+")
      csv <- tempfile(fileext = ".csv")
      on.exit(unlink(csv), add = TRUE)
      utils::write.csv(out, csv, row.names = FALSE)
      expect_equal(utils::read.csv(csv)$Variance, out$Variance)
    }
  }
  # A short seeded run verifies scale invariance of the interval path, not coverage.
  base_boot <- compute_facet_icc(data, c("Rater", "Criterion"), "Score", "Person",
                                ci_method = "boot", ci_boot_reps = 4, ci_boot_seed = 2026)
  data$Score <- data$Score * 1e-4
  small_boot <- compute_facet_icc(data, c("Rater", "Criterion"), "Score", "Person",
                                 ci_method = "boot", ci_boot_reps = 4, ci_boot_seed = 2026)
  expect_identical(small_boot$ICC_CI_Status, base_boot$ICC_CI_Status)
  expect_equal(attr(small_boot, "icc_ci")$bootstrap$draws,
               attr(base_boot, "icc_ci")$bootstrap$draws, tolerance = 1e-5)
  expect_equal(small_boot$ICC_CI_Lower, base_boot$ICC_CI_Lower, tolerance = 1e-4)
  expect_equal(small_boot$ICC_CI_Upper, base_boot$ICC_CI_Upper, tolerance = 1e-4)
})

test_that("constant retained scores remain unidentified at zero and nonzero values", {
  testthat::local_mocked_bindings(.mfrm_boot_mer = function(...) stop("must not run"), .package = "mfrmr")
  data <- load_mfrmr_data("example_core")
  for (value in c(0, 5, 1e-6)) {
    data$Score <- value
    expect_message(out <- compute_facet_icc(data, c("Rater", "Criterion"), "Score", "Person",
      ci_method = "boot", ci_boot_reps = 2), "ICCs are undefined")
    expect_true(all(is.na(out$ICC) & is.na(out$Variance)))
    expect_true(all(out$Interpretation == "Non-identifiable"))
    expect_true(all(out$ICC_CI_Status == "Undefined ICC" & out$ICC_CI_NReps == 0))
    expect_true(all(is.na(out$ICC_CI_Lower) & is.na(out$ICC_CI_Upper)))
    deff <- compute_facet_design_effect(data, "Rater", out)
    expect_true(is.na(deff$DesignEffect) && is.na(deff$EffectiveN))
    expect_output(print(deff), "per-facet approximation")
    expect_output(summary(deff), "descriptive row count")
  }
})

test_that("ICC preserves numeric score labels and literal column names", {
  data <- load_mfrmr_data("example_core")
  data$Score <- as.numeric(data$Score)^2 + 10
  baseline <- compute_facet_icc(data, c("Rater", "Criterion"), "Score", "Person")
  data$Score <- factor(data$Score, levels = rev(sort(unique(data$Score))))
  labelled <- compute_facet_icc(data, c("Rater", "Criterion"), "Score", "Person")
  expect_identical(labelled$Variance, baseline$Variance)
  expect_identical(labelled$ICC, baseline$ICC)
  data$Score <- as.character(data$Score)
  names(data)[names(data) == "Rater"] <- "Rater group"
  names(data)[names(data) == "Score"] <- "Rating score"
  named <- compute_facet_icc(data, c("Rater group", "Criterion"), "Rating score", "Person")
  expect_identical(named$Variance, baseline$Variance)
  expect_true("Rater group" %in% named$Facet)
  deff <- compute_facet_design_effect(data, "Rater group", named)
  expect_equal(deff$AvgClusterSize, nrow(data) / length(unique(data$`Rater group`)))
})

test_that("ICC omission and design effects use exactly the same retained rows", {
  data <- load_mfrmr_data("example_core")
  data$Rater <- as.character(data$Rater)
  baseline <- compute_facet_icc(data, c("Rater", "Criterion"), "Score", "Person")
  n <- nrow(data)
  extra <- data[1:3, ]
  extra$Score[c(1, 3)] <- NA_real_
  extra$Rater[1:2] <- c("Only in an excluded row", NA_character_)
  extra$Person[3] <- NA
  input <- rbind(data, extra)
  input$Unselected <- NA_real_
  expect_error(compute_facet_icc(input, c("Rater", "Criterion"), "Score", "Person"),
               "3 row.*missing")
  # Global na.action must not silently change which observations enter the fit.
  old_options <- options(na.action = "na.omit")
  on.exit(options(old_options), add = TRUE)
  out <- compute_facet_icc(input, c("Rater", "Criterion"), "Score", "Person", missing = "omit")
  expect_identical(out$Variance, baseline$Variance)
  expect_identical(out$ICC, baseline$ICC)
  usage <- attr(out, "data_usage")
  expect_identical(usage$missing, "omit")
  expect_equal(usage$counts, c(InputRows = n + 3L, UsedRows = n, ExcludedRows = 3L))
  expect_equal(usage$excluded_rows, n + 1:3)
  expect_equal(usage$missing_cells, data.frame(
    InputRow = n + c(1L, 3L, 3L, 2L), Column = c("Score", "Score", "Person", "Rater")))
  expect_equal(usage$observed_levels[["Rater"]], length(unique(data$Rater)))
  expect_equal(out$ExcludedRows, rep(3L, nrow(out)))
  expect_output(print(out), "3 excluded")
  expect_output(summary(out), "no missing values were imputed")

  deff <- compute_facet_design_effect(input, c("Rater", "Criterion"), out)
  k <- vapply(data[c("Rater", "Criterion")], function(x) length(unique(x)), integer(1))
  rho <- out$ICC[match(names(k), out$Facet)]
  expect_equal(deff$AvgClusterSize, unname(round(n / k, 3)))
  expect_equal(deff$DesignEffect, unname(round(1 + (n / k - 1) * rho, 3)))
  expect_equal(deff$EffectiveN, unname(round(n / (1 + (n / k - 1) * rho), 1)))
  expect_identical(attr(deff, "data_usage"), usage)
  expect_equal(deff$UsedRows, rep(n, 2))
  expect_output(print(deff), "3 excluded")
  expect_output(summary(deff), "3 excluded")
  # Saved ICC counts, not a later copy of the input, define the DEFF sample.
  expect_identical(compute_facet_design_effect(data[1, ], c("Rater", "Criterion"), out), deff)
  direct <- compute_facet_design_effect(input, c("Rater", "Criterion"),
                                        score = "Score", person = "Person", missing = "omit")
  expect_identical(direct, deff)
  expect_error(compute_facet_design_effect(input, "Rater", score = "Score", person = "Person"),
               "missing.*omit")
  expect_error(compute_facet_design_effect(input, "Unselected", out), "grouping column")
  expect_error(compute_facet_design_effect(input, "Residual", out), "grouping column")
  old <- out
  attr(old, "data_usage") <- NULL
  expect_output(print(old), "Row accounting unavailable")
  expect_error(compute_facet_design_effect(input, "Rater", old), "Rerun compute_facet_icc")

  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv))
  utils::write.csv(deff, csv, row.names = FALSE)
  expect_equal(utils::read.csv(csv)$ExcludedRows, c(3L, 3L))
})

test_that("malformed ICC inputs cannot be hidden by omission", {
  data <- load_mfrmr_data("example_core")
  for (value in c("not recorded", " ", "Inf", "NaN")) {
    invalid <- data
    invalid$Score <- as.character(invalid$Score)
    invalid$Score[1] <- value
    expect_error(compute_facet_icc(invalid, "Rater", "Score", missing = "omit"),
                 "nonnumeric or infinite")
  }
  for (value in list(Inf, -Inf, as.Date("2026-01-01"), TRUE, 1i, list(1))) {
    invalid <- data
    invalid$Score <- rep(value, nrow(data))
    expect_error(compute_facet_icc(invalid, "Rater", "Score", missing = "omit"),
                 "score column")
  }
  invalid <- data
  invalid$Score <- NA_real_
  expect_error(compute_facet_icc(invalid, "Rater", "Score", missing = "omit"), "No complete rows")
  invalid <- data
  invalid$Rater <- as.numeric(factor(invalid$Rater))
  invalid$Rater[1] <- NaN
  expect_error(compute_facet_icc(invalid, "Rater", "Score"), "1 row.*missing")
  invalid$Rater[1] <- Inf
  expect_error(compute_facet_icc(invalid, "Rater", "Score", missing = "omit"), "finite numeric")
  invalid$Rater <- as.character(invalid$Rater)
  invalid$Rater[1] <- " "
  expect_error(compute_facet_icc(invalid, "Rater", "Score", missing = "omit"), "blank labels")
  expect_error(compute_facet_icc(data, "Rater", "Score", reml = "TRUE"), "TRUE or FALSE")
  expect_error(compute_facet_icc(data, "Rater", "Rater"), "Grouping columns")
  expect_error(compute_facet_icc(data, c("Rater", "Rater"), "Score"), "distinct facet")
})

test_that("hierarchical entry points retain omission accounting and propagate input errors", {
  data <- load_mfrmr_data("example_core")
  data$Score[1] <- NA_real_
  data$Rater[2] <- NA
  data$Person[3] <- NA
  expect_error(analyze_hierarchical_structure(data, c("Rater", "Criterion"), igraph_layout = FALSE),
               "3 row.*missing")
  out <- analyze_hierarchical_structure(data, c("Rater", "Criterion"),
                                        igraph_layout = FALSE, missing = "omit")
  expect_identical(attr(out$icc, "data_usage"), attr(out$design_effect, "data_usage"))
  expect_equal(out$design_effect$UsedRows, rep(nrow(data) - 3L, 2))
  expect_output(print(out), "3 excluded")
  expect_output(summary(out), "3 excluded")
  expect_equal(sum(out$crosstabs[[1]]$N), nrow(data))
  expect_s3_class(analyze_hierarchical_structure(data, c("Rater", "Criterion"),
                   compute_icc = FALSE, igraph_layout = FALSE), "mfrm_hierarchical_structure")
  # This wrapper reads prep only; no MFRM estimation is needed for this case.
  fit <- structure(list(prep = list(data = data[-(1:3), ], facet_names = c("Person", "Rater", "Criterion"))),
                   class = "mfrm_fit")
  stored <- analyze_hierarchical_structure(fit, igraph_layout = FALSE)
  expect_identical(attr(stored$icc, "data_usage")$source, "Stored fitted rows")
  expect_identical(attr(stored$icc, "data_usage"), attr(stored$design_effect, "data_usage"))
  expect_output(print(stored), "earlier MFRM filtering is not included")
})
