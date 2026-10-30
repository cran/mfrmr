test_that("fit screens keep unknown classifications and exact denominators", {
  tbl <- data.frame(Infit = c(2, 2, 1, 1, NA, Inf, -1),
                    Outfit = c(1, NA, 1, NA, NA, 1, 1),
                    InfitZSTD = c(3, 3, 1, 1, NA, Inf, 0),
                    OutfitZSTD = c(1, NA, 1, NA, NA, 1, 0))
  screen <- fit_screening_summary(tbl)
  expect_identical(screen$flags$MeanSquare, c(TRUE, TRUE, FALSE, NA, NA, NA, NA))
  expect_equal(screen$counts$Classified, c(3, 4, 4))
  expect_equal(screen$counts$Unclassified, c(4, 3, 3))
  expect_equal(screen$counts$IncompleteStatistics, c(5, 4, 4))
  expect_equal(screen$counts$Flagged, c(2, 2, 2))
  expect_true(all(is.na(screen$counts$FlagRate)))
  expect_match(fit_screening_note(screen$counts[1, ]), "2 of 3 classified")
  none <- fit_screening_summary(data.frame(Infit = NA, Outfit = NA))
  expect_true(all(is.na(none$counts$Flagged)))
  expect_match(fit_screening_note(none$counts[1, ]), "no elements could be classified")
  expect_equal(fit_screening_summary(tbl[1, ])$counts$FlagRate, c(1, 1, 1))
  expect_equal(fit_screening_summary(data.frame())$counts$Elements, rep(0, 3))
})

test_that("fit rankings never mix ZSTD and mean-square deviation units", {
  tbl <- data.frame(Facet = "Rater", Level = c("A", "B"),
                    Infit = c(1, 4), Outfit = c(1, 1),
                    InfitZSTD = c(8, NA), OutfitZSTD = c(8, NA))
  text <- summarize_top_misfit_levels(tbl)
  expect_match(text, "Rater:A (|ZSTD| = 8.00)", fixed = TRUE)
  expect_false(grepl("Rater:B|MnSq", text))
  tbl$InfitZSTD[] <- NA_real_
  tbl$OutfitZSTD[] <- NA_real_
  expect_match(summarize_top_misfit_levels(tbl), "Rater:B (|MnSq - 1| = 3.00)", fixed = TRUE)
})

test_that("reports and diagnostic summaries retain incomplete fit screens", {
  fit <- make_toy_fit()
  diag <- make_toy_diagnostics(fit)
  for (slot in c("fit", "measures")) {
    for (name in c("Infit", "Outfit", "InfitZSTD", "OutfitZSTD")) {
      diag[[slot]][[name]][] <- NA_real_
    }
  }
  diag$overall_fit$Infit <- diag$overall_fit$Outfit <- NA_real_
  contract <- build_apa_reporting_contract(fit, diag)
  expect_true(is.na(contract$summaries$misfit_n))
  expect_false(grepl("CIEligible|CI_Lower|CI_Upper", paste(contract$report_text, collapse = " ")))
  expect_match(paste(contract$report_text, collapse = " "), "Overall mean-square screening was unavailable")
  expect_false(grepl("Overall mean-square fit was outside", paste(contract$report_text, collapse = " ")))
  s <- summary(diag)
  expect_true(all(is.na(s$fit_screening$Flagged)))
  expect_equal(nrow(s$top_fit), 0)
  expect_true(any(grepl("no elements could be classified", s$key_warnings)))
  printed <- paste(capture.output(print(s)), collapse = " ")
  expect_match(printed, "Element-fit screening")
  expect_match(printed, "no elements could be classified")
  expect_false(grepl("follow_up_needed", printed))
  bundle <- build_summary_table_bundle(s)
  expect_equal(bundle$tables$fit_screening, as.data.frame(s$fit_screening))
  s$fit_screening <- NULL
  expect_error(print(s), "Recreate this diagnostic summary")

  for (slot in c("fit", "measures")) {
    diag[[slot]]$Infit[1:3] <- c(2, 2, 1)
    diag[[slot]]$Outfit[1:3] <- c(NA, NA, 1)
  }
  expect_match(paste(build_visual_summary_map(fit, diag)$fit_diagnostics, collapse = " "),
               "2 of 3 classified elements flagged")
  expect_match(paste(build_visual_warning_map(fit, diag)$fit_diagnostics, collapse = " "),
               "unclassified")
  s <- summary(diag)
  expect_equal(nrow(s$misfit_flagged), 2)
  expect_equal(s$fit_screening$Flagged[1], 2)
  expect_true(is.na(s$fit_screening$FlagRate[1]))
  diag$fit <- data.frame()
  empty <- summary(diag)
  expect_true(all(is.na(empty$fit_screening$Flagged)))
  expect_true(any(grepl("no elements could be classified", empty$key_warnings)))
})
