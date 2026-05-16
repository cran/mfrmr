# Reporting blocks added in 0.1.6: facets_chisq / interrater summary /
# MnSq misfit auto-flag in summary.mfrm_diagnostics; targeting in
# summary.mfrm_fit. These tests pin the contract surfaced by the new
# print blocks so we can detect accidental shape changes without
# snapshotting the literal console output (which would tie us to a
# specific digit layout).

# Toy fit + diagnostics shared by every test_that() block in this
# file. `local_toy_fit()` (helper-fixtures.R) caches the JML fit and
# the legacy diagnostics so re-running the file does not pay the
# fitting cost more than once.
.toy <- load_mfrmr_data("example_core")
.fit <- make_toy_fit()
.diag <- make_toy_diagnostics(.fit)

# --- summary.mfrm_diagnostics ---------------------------------------------

test_that("summary(diag) surfaces facets_chisq", {
  s <- summary(.diag)
  expect_s3_class(s, "summary.mfrm_diagnostics")
  expect_true(!is.null(s$facets_chisq))
  expect_gt(nrow(s$facets_chisq), 0L)
  expect_true(all(c("Facet", "FixedChiSq", "FixedDF", "FixedProb") %in%
                    names(s$facets_chisq)))
})

test_that("summary(diag) surfaces inter-rater agreement summary", {
  s <- summary(.diag)
  expect_true(!is.null(s$interrater))
  expect_gt(nrow(s$interrater), 0L)
  ir_cols <- names(s$interrater)
  expect_true(all(c("Raters", "Pairs", "ExactAgreement", "MeanCorr") %in% ir_cols))
})

test_that("summary(diag) carries the MnSq misfit threshold pair", {
  s <- summary(.diag)
  expect_named(s$misfit_thresholds, c("lower", "upper"))
  expect_equal(unname(s$misfit_thresholds), c(0.5, 1.5))
})

test_that("summary(diag) auto-flag names worst MnSq element", {
  # Force a misfit-flagged diagnostic by injecting an Outfit > 1.5 row.
  poisoned <- .diag
  poisoned$fit$Outfit[1] <- 2.4
  poisoned$fit$Infit[1] <- 1.7
  poisoned$fit$OutfitZSTD[1] <- 3.1
  poisoned$fit$InfitZSTD[1] <- 2.6
  s <- summary(poisoned)
  joined <- paste(s$key_warnings, collapse = " | ")
  expect_true(grepl("MnSq misfit", joined))
  # Specific element-level naming: "Facet:Level (Infit=..., Outfit=...)"
  expect_true(grepl("Outfit=2\\.4", joined))
})

test_that("print(summary(diag)) emits the new blocks without error", {
  s <- summary(.diag)
  out <- utils::capture.output(print(s))
  joined <- paste(out, collapse = "\n")
  expect_true(grepl("Facet variability", joined))
  expect_true(grepl("Inter-rater agreement summary", joined))
})

# --- summary.mfrm_fit -----------------------------------------------------

test_that("summary(fit) surfaces a targeting block", {
  s <- summary(.fit)
  expect_s3_class(s, "summary.mfrm_fit")
  expect_true(!is.null(s$targeting))
  expect_gt(nrow(s$targeting), 0L)
  expect_true(all(c("Facet", "PersonMean", "FacetMean", "Targeting",
                     "PersonSD", "FacetSD", "SpreadRatio") %in%
                    names(s$targeting)))
})

test_that("summary(fit) targeting matches Person mean (sum-to-zero ID)", {
  s <- summary(.fit)
  expect_equal(
    as.numeric(s$targeting$Targeting),
    rep(as.numeric(s$person_overview$Mean[1]), nrow(s$targeting)) -
      as.numeric(s$targeting$FacetMean),
    tolerance = 1e-10
  )
})

test_that("print(summary(fit)) shows the targeting block", {
  s <- summary(.fit)
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("Targeting", out)))
})
