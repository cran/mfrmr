test_that("planned trials, dependent targets and unknown family outcomes stay distinct", {
  roster <- expand.grid(Condition = "Null", Replicate = 1:4, Target = c("R1", "R2"))
  roster$Affected <- FALSE
  results <- roster[c("Condition", "Replicate", "Target")]
  results$Flag <- c(FALSE, TRUE, FALSE, NA, FALSE, NA, FALSE, NA)
  x <- mfrm_screening_performance(roster, results, "Prespecified test rule")
  family <- x$by_family[x$by_family$Targets > 0, ]
  expect_equal(family$Planned, 4)
  expect_equal(family$Available, 3)
  expect_equal(family$CompleteScreens, 2)
  expect_equal(family$Positive, 1)
  expect_equal(family$Rate, 1 / 3)
  expect_equal(family$AllTrialsLower, .25)
  expect_equal(family$AllTrialsUpper, .5)
  expect_identical(x$family_outcomes$Flag, c(FALSE, TRUE, FALSE, NA))
  expect_equal(x$by_target$Available, c(3, 2))
  expect_true(is.na(x$by_family$Rate[x$by_family$Targets == 0]))
  omitted <- mfrm_screening_performance(roster, results[results$Replicate != 4, ], "Same rule")
  expect_equal(omitted$by_family, x$by_family)
  expect_equal(sum(!omitted$outcomes$Reported), 2)
  empty <- mfrm_screening_performance(roster, results[FALSE, ], "No results returned")
  expect_true(all(empty$by_target$Unavailable == 4))
  expect_true(all(is.na(empty$by_target$Rate)))
  expect_equal(empty$by_target$AllTrialsLower, c(0, 0))
  expect_equal(empty$by_target$AllTrialsUpper, c(1, 1))
  saved <- tempfile(); saveRDS(x, saved)
  expect_identical(readRDS(saved), x)
  expect_equal(summary(x)$by_target, x$by_target)
})

test_that("known truth and exact binomial bounds preserve their statistical meaning", {
  roster <- expand.grid(Condition = "Test", Replicate = 1:100, Target = c("A", "B"))
  roster$Affected <- roster$Target == "B"
  results <- roster[c("Condition", "Replicate", "Target")]
  results$Flag <- roster$Affected
  x <- mfrm_screening_performance(roster, results, "Declared alternative at B")
  expect_equal(x$by_target$Rate, c(0, 1))
  expect_equal(x$by_target$MCUpper[1], 1 - .025^(1 / 100))
  expect_equal(x$by_target$MCLower[2], .025^(1 / 100))
  expect_equal(x$by_family$Rate, c(0, 1))
  expect_equal(x$by_target$MCSE, c(0, 0))
  expect_gt(x$by_target$MCUpper[1], 0)
  expect_lt(x$by_target$MCLower[2], 1)
})

test_that("identity checks prevent silent target changes and dropped planned trials", {
  roster <- data.frame(Condition = c("A\rB", "A"), Replicate = c("C", "B\rC"),
    Target = "D", Affected = FALSE)
  results <- roster[2:1, 1:3]; results$Flag <- c(FALSE, TRUE)
  x <- mfrm_screening_performance(roster, results, "Identity check")
  expect_identical(x$outcomes$Flag, c(TRUE, FALSE))
  expect_error(mfrm_screening_performance(roster, rbind(results, results[1, ]), "Rule"), "Duplicate")
  bad <- results; bad$Target[1] <- "Extra"
  expect_error(mfrm_screening_performance(roster, bad, "Rule"), "outside")
  bad <- roster; bad$Affected[1] <- NA
  expect_error(mfrm_screening_performance(bad, results, "Rule"), "missing truth")
  bad <- results; bad$Flag <- c(0, 1)
  expect_error(mfrm_screening_performance(roster, bad, "Rule"), "logical Flag")
  bad <- results; bad$Flag <- matrix(FALSE, nrow = 2, ncol = 2)
  expect_error(mfrm_screening_performance(roster, bad, "Rule"), "logical Flag")
  expect_error(mfrm_screening_performance(roster, results, ""), "prespecified")
  expect_error(mfrm_screening_performance(roster, results, "Rule", level = 1), "level")
  full <- expand.grid(Condition = "C", Replicate = 1:2, Target = c("A", "B")); full$Affected <- FALSE
  expect_error(mfrm_screening_performance(full[-1, ], results, "Rule"), "same planned targets")
  full$Affected[1] <- TRUE
  expect_error(mfrm_screening_performance(full, results, "Rule"), "truth labels")
})

test_that("performance graphics preserve unavailable outcomes without opening hidden devices", {
  roster <- expand.grid(Condition = c("Known", "Unavailable"), Replicate = 1:4, Target = "R1")
  roster$Affected <- FALSE
  results <- roster[1:3]; results$Flag <- ifelse(results$Condition == "Known", FALSE, NA)
  x <- mfrm_screening_performance(roster, results, "Rule")
  before <- grDevices::dev.cur()
  z <- plot(x, draw = FALSE)
  expect_identical(grDevices::dev.cur(), before)
  expect_equal(plot_data(z)$table$Available, c(4, 0))
  expect_true(is.na(plot_data(z)$table$Rate[2]))
  expect_equal(nrow(plot_data(plot(x, condition = "Known", scope = "target", draw = FALSE))$table), 1)
  expect_error(plot(x, metric = "detection", draw = FALSE), "No planned targets")
  expect_error(plot(x, condition = "Missing", draw = FALSE), "condition labels")
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_error(as_ggplot(z), "screening performance")
  grDevices::pdf(tempfile(fileext = ".pdf"))
  before <- graphics::par("mar"); plot(x)
  expect_equal(graphics::par("mar"), before)
  grDevices::dev.off()
})

test_that("unavailable bias statistics cannot become negative screens or lower false-flag rates", {
  local_mocked_bindings(estimate_bias = function(...) list(table = data.frame(
    FacetA_Level = c("R01", "R01", "R02", "R02"),
    FacetB_Level = c("C01", "C02", "C01", "C02"),
    `Prob.` = c(.001, NA, .8, NA), t = c(3, NA, 0, NA),
    `Bias Size` = c(.4, NA, .1, NA), check.names = FALSE)))
  x <- suppressWarnings(evaluate_mfrm_signal_detection(n_person = 12, n_rater = 2,
    n_criterion = 2, raters_per_person = 2, reps = 1, maxit = 30, seed = 772))
  expect_true(is.na(x$results$BiasDetected))
  expect_equal(x$results$BiasScreenFalsePositiveRate, .5)
  expect_equal(x$results$BiasNonTargetPlanned, 3)
  expect_equal(x$results$BiasNonTargetAvailable, 2)
  s <- summary(x)$detection_summary
  expect_true(is.na(s$BiasScreenRate))
  expect_equal(s$BiasScreenPlanned, 1)
  expect_equal(s$BiasScreenAvailable, 0)
  expect_equal(s$BiasScreenAllTrialsLower, 0)
  expect_equal(s$BiasScreenAllTrialsUpper, 1)
  expect_true(all(is.na(plot(x, signal = "bias", draw = FALSE)$data$y)))
  old <- x; old$results$BiasDetected <- FALSE
  old$results$BiasNonTargetAvailable <- NULL
  expect_true(is.na(summary(old)$detection_summary$BiasScreenRate))
  expect_true(is.na(summary(old)$detection_summary$BiasScreenFalsePositiveRate))
  expect_true(any(grepl("Rerun the simulation", summary(old)$notes)))
  expect_true(is.na(plot(old, signal = "bias", metric = "false_positive", draw = FALSE)$data$y))
  old$results$DIFDetected <- TRUE
  expect_true(is.na(plot(old, signal = "dif", draw = FALSE)$data$y))
  expect_identical(signal_eval_threshold(c(.01, NA, .8), c(NA, 3, .1), .05, 2), c(NA, NA, FALSE))
})

test_that("descriptive DIF classifications and saved missing contrasts remain unavailable", {
  local_mocked_bindings(analyze_dff = function(...) list(dif_table = data.frame(
    Level = c("C01", "C02"), Group1 = "A", Group2 = "B", Contrast = c(NA, .8),
    p_value = .01, p_adjusted = .02, ETS = NA_character_,
    Classification = "Descriptive", ClassificationSystem = "descriptive")))
  x <- suppressWarnings(evaluate_mfrm_signal_detection(n_person = 12, n_rater = 2,
    n_criterion = 2, raters_per_person = 2, dif_method = "refit", dif_level = "C02",
    reps = 1, maxit = 30, seed = 773))
  expect_true(x$results$DIFDetected)
  expect_true(is.na(x$results$DIFClassDetected))
  expect_true(is.na(x$results$DIFFalsePositiveRate))
  expect_equal(x$results$DIFNonTargetAvailable, 0)
  # No additional model fits: three saved independent-row placeholders test aggregation.
  x$results <- x$results[rep(1, 3), ]; x$results$rep <- 1:3
  x$results$DIFDetected <- c(TRUE, FALSE, FALSE)
  x$results$DIFClassDetected <- FALSE
  expect_equal(summary(x, digits = 0)$detection_summary$DIFPower, 1 / 3)
  expect_true(is.na(summary(x)$detection_summary$DIFClassificationPower))
  old <- x; old$results$DIFContrast <- NA_real_
  expect_true(is.na(summary(old)$detection_summary$DIFPower))
  expect_true(all(is.na(plot(old, draw = FALSE)$data$y)))
  expect_error(evaluate_mfrm_signal_detection(dif_p_cut = NA_real_), "dif_p_cut")
  expect_error(evaluate_mfrm_signal_detection(bias_p_cut = 2), "bias_p_cut")
  expect_error(evaluate_mfrm_signal_detection(bias_abs_t = -1), "bias_abs_t")
})
