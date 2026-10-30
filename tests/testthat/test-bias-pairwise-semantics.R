pairwise_reference_bundle <- function(sign = -1, counts = c(20, 20), rows = c(5, 5)) {
  b <- data.frame(Sq = 1:2, check.names = FALSE,
    `Observd Score` = 10, `Expctd Score` = 10,
    `Observd Count` = counts, `Obs-Exp Average` = 0,
    `Bias Size` = c(1, -1) * 2.2 * sqrt(2) * .2 / 2,
    `S.E.` = .2, t = 0, `d.f.` = counts - 1, `Prob.` = 1,
    Infit = 1, Outfit = 1, ObsN = rows,
    FacetA = "Reviewer", FacetB = "Task", FacetA_Level = "R1",
    FacetB_Level = c("A", "B"), FacetA_Index = 1, FacetB_Index = 1:2,
    FacetA_Measure = .7, FacetA_SE = NA_real_, FacetB_Measure = 0, FacetB_SE = .3)
  list(table = b, summary = data.frame(), chi_sq = data.frame(),
    facet_a = "Reviewer", facet_b = "Task",
    orientation_review = data.frame(Facet = c("Reviewer", "Task"), Sign = c(sign, -1)))
}

test_that("local contrasts solve the oriented predictor for either target sign", {
  for (sign in c(-1, 1)) {
    b <- pairwise_reference_bundle(sign)
    for (branch in c("facets", "original")) {
      fixed <- build_fixed_reports(b, branch = branch)
      p <- fixed$pairwise_table
      expect_equal(sign * p$`Local Measure1`, sign * .7 + b$table$`Bias Size`[1])
      expect_equal(sign * p$`Local Measure2`, sign * .7 + b$table$`Bias Size`[2])
      expect_equal(p$Contrast, sign * 2.2 * sqrt(2) * .2)
      expect_equal(p$SE1, .2)
      expect_equal(p$SE2, .2)
      expect_equal(p$SE, sqrt(.08))
      expect_equal(p$t, sign * 2.2)
      expect_equal(p$`d.f.`, 38)
      expect_equal(p$`Prob.`, 2 * pt(-2.2, 38))
      expect_identical(p$CalculationVersion, 3L)
      expect_identical(p$ScreenStatus, "available_conditional")
      expect_equal(plot(fixed, draw = FALSE)$data$table$Contrast, p$Contrast)
    }
    report <- bias_pairwise_report(b)
    expect_equal(report$table$Contrast, p$Contrast)
    expect_equal(report$summary$ScreenAvailable, 1)
  }
})

test_that("Welch screening is invariant to replication encoding", {
  compressed <- bias_pairwise_report(pairwise_reference_bundle())
  expanded <- bias_pairwise_report(pairwise_reference_bundle(rows = c(20, 20)))
  cols <- c("Contrast", "SE", "t", "d.f.", "Prob.", "Flag")
  expect_equal(compressed$table[cols], expanded$table[cols])
  expect_equal(compressed$table$Count1, 20)
  expect_equal(compressed$table$ObsN1, 5)
  expect_equal(expanded$table$ObsN1, 20)
  expect_lt(compressed$table$`Prob.`, .05)
  # The df reference remains defined for positive fractional replication counts > 1.
  fractional <- bias_pairwise_report(pairwise_reference_bundle(counts = c(2.5, 4.5)),
    df_method = "classical")$table
  expect_equal(fractional$`d.f.`, .08^2 / (.04^2 / 1.5 + .04^2 / 3.5))
})

test_that("unavailable constituents never become an available selected-subset screen", {
  for (bad in c(0, 1, NA_real_)) {
    b <- pairwise_reference_bundle(counts = c(bad, 5))
    p <- bias_pairwise_report(b)
    expect_true(is.finite(p$table$Contrast))
    expect_true(is.na(p$table$`d.f.`))
    expect_true(is.na(p$table$`Prob.`))
    expect_true(is.na(p$table$Flag))
    expect_equal(p$summary$ScreenUnavailable, 1)
    expect_true(is.na(p$summary$Flagged))
  }
  expect_true(is.na(welch_satterthwaite_df(c(.04, .04), c(0, 4))))
  expect_true(is.na(welch_satterthwaite_df(c(NA_real_, .04), c(4, 4))))
  expect_true(is.na(welch_satterthwaite_df(c(.04, .04), 4)))
  expect_equal(welch_satterthwaite_df(c(.04, .09), c(4, 8)),
    .13^2 / (.04^2 / 4 + .09^2 / 8))
  b <- pairwise_reference_bundle()
  b$table$`S.E.`[1] <- NA_real_
  p <- bias_pairwise_report(b)
  expect_true(is.na(p$table$SE))
  expect_true(is.na(p$table$`d.f.`))
  expect_true(is.na(p$table$Flag))
})

test_that("orientation comes from saved metadata and legacy reports require recreation", {
  b <- pairwise_reference_bundle()
  b$orientation_review <- NULL
  missing <- bias_pairwise_report(b)
  expect_true(is.na(missing$table$Contrast))
  expect_true(is.na(missing$table$Flag))
  expect_match(missing$table$ScreenReason, "orientation")
  expect_true(is.na(missing$summary$Flagged))
  fixed <- build_fixed_reports(pairwise_reference_bundle())
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  saveRDS(fixed, saved)
  expect_equal(plot(readRDS(saved), draw = FALSE)$data$table$Contrast,
    fixed$pairwise_table$Contrast)
  old <- fixed
  old$pairwise_table$CalculationVersion <- NULL
  expect_error(summary(old), "Recreate build_fixed_reports", fixed = TRUE)
  expect_error(plot(old, draw = FALSE), "Recreate build_fixed_reports", fixed = TRUE)
  old <- bias_pairwise_report(pairwise_reference_bundle())
  old$table$CalculationVersion <- NULL
  expect_error(summary(old), "Recreate build_fixed_reports", fixed = TRUE)
})

test_that("pairwise summaries retain the complete family before display truncation", {
  b <- pairwise_reference_bundle()
  third <- b$table[2, ]
  third$FacetB_Level <- "C"
  third$FacetB_Index <- 3
  third$`S.E.` <- NA_real_
  b$table <- rbind(b$table, third)
  out <- bias_pairwise_report(b, top_n = 1)
  expect_equal(nrow(out$table), 1)
  expect_equal(out$summary$Contrasts, 3)
  expect_equal(out$summary$ScreenAvailable, 1)
  expect_equal(out$summary$ScreenUnavailable, 2)
  all_rows <- bias_pairwise_report(b, top_n = 10)
  expect_equal(nrow(all_rows$table), 3)
  expect_equal(out$summary, all_rows$summary)
})

test_that("the full-precision screening boundary and target-facet reversal are explicit", {
  for (t in c(2 - 1e-6, 2, 2 + 1e-6)) {
    b <- pairwise_reference_bundle(sign = 1, counts = c(5, 5))
    b$table$FacetA_Measure <- 0
    b$table$`S.E.` <- c(.6, .8)
    b$table$`Bias Size` <- c(t, 0)
    out <- bias_pairwise_report(b)
    expect_equal(out$table$t, t)
    expect_identical(out$table$Flag, t >= 2)
  }
  b <- pairwise_reference_bundle()
  for (suffix in c("", "_Level", "_Index", "_Measure", "_SE")) {
    a <- paste0("FacetA", suffix); z <- paste0("FacetB", suffix)
    tmp <- b$table[[a]]; b$table[[a]] <- b$table[[z]]; b$table[[z]] <- tmp
  }
  b$facet_a <- "Task"; b$facet_b <- "Reviewer"
  expect_lt(bias_pairwise_report(b, target_facet = "Reviewer")$table$Contrast, 0)
  b <- pairwise_reference_bundle(sign = 1)
  b$table$FacetA_Measure <- 1e16
  b$table$`Bias Size` <- c(.1, -.1)
  expect_equal(bias_pairwise_report(b)$table$Contrast, .2)
  for (scale in c(1e-200, 1e200)) {
    b$table$FacetA_Measure <- 0
    b$table$`Bias Size` <- c(scale, -scale)
    b$table$`S.E.` <- scale
    p <- bias_pairwise_report(b)$table
    expect_equal(p$t, sqrt(2))
    expect_equal(p$`d.f.`, 38)
  }
})

test_that("retained calibration reproduces identical cell and pairwise precision under replication", {
  f <- readRDS(test_path("fixtures", "anchored-facet-information.rds"))$RSM$fit
  measures <- dplyr::bind_rows(
    dplyr::transmute(f$facets$person, Facet = "Person", Level = as.character(Person), Estimate),
    dplyr::select(f$facets$others, Facet, Level, Estimate))
  measures$SE <- NA_real_
  weighted <- f; weighted$prep$data$Weight <- 4
  expanded <- f
  expanded$prep$data <- expanded$prep$data[rep(seq_len(nrow(expanded$prep$data)), each = 4), ]
  expanded$prep$data$Weight <- 1
  result <- lapply(list(weighted, expanded), function(fit) {
    b <- estimate_bias_interaction(fit, list(obs = compute_obs_table(fit), measures = measures),
      facet_a = "Rater", facet_b = "Criterion")
    build_fixed_reports(b)$pairwise_table
  })
  cols <- c("Local Measure1", "Local Measure2", "SE1", "SE2", "Contrast", "SE", "t", "d.f.", "Prob.")
  expect_equal(result[[1]][cols], result[[2]][cols], tolerance = 1e-8)
  expect_equal(result[[1]]$Count1, result[[2]]$Count1)
  expect_equal(result[[1]]$ObsN1 * 4, result[[2]]$ObsN1)
  expect_true(all(result[[1]]$TargetSign == -1))
})
