# Tests for fit_mfrm(attach_diagnostics = TRUE) — the opt-in merge of
# per-level SE / Infit / Outfit / PtMeaCorr onto fit$facets$others
# added in 0.1.6. Verifies the documented column contract on the
# `others` table and the PtMeaCorr column name.

local({
  .toy <<- load_mfrmr_data("example_core")
})

test_that("attach_diagnostics = FALSE preserves the minimal others table", {
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15)
  ))
  expect_false(isTRUE(fit$config$attached_diagnostics))
  expect_setequal(names(fit$facets$others),
                  c("Facet", "Level", "Estimate"))
})

test_that("attach_diagnostics = TRUE merges SE / Infit / Outfit / PtMeaCorr", {
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15, attach_diagnostics = TRUE)
  ))
  expect_true(isTRUE(fit$config$attached_diagnostics))
  expected <- c("Facet", "Level", "Estimate", "SE", "Infit",
                "Outfit", "PtMeaCorr")
  expect_true(all(expected %in% names(fit$facets$others)))
  # Values should be finite for at least one row of each metric.
  expect_true(any(is.finite(fit$facets$others$SE)))
  expect_true(any(is.finite(fit$facets$others$Infit)))
  expect_true(any(is.finite(fit$facets$others$Outfit)))
  expect_true(any(is.finite(fit$facets$others$PtMeaCorr)))
  # Estimate column is untouched by the merge.
  fit_bare <- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15)
  ))
  key_bare <- paste(fit_bare$facets$others$Facet,
                    fit_bare$facets$others$Level, sep = "|")
  key_full <- paste(fit$facets$others$Facet,
                    fit$facets$others$Level, sep = "|")
  match_idx <- match(key_bare, key_full)
  expect_equal(unname(fit$facets$others$Estimate[match_idx]),
               unname(fit_bare$facets$others$Estimate),
               tolerance = 1e-10)
})

test_that("attach_diagnostics records the attached columns in config", {
  fit <- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15, attach_diagnostics = TRUE)
  ))
  expect_true("attached_diagnostics_cols" %in% names(fit$config))
  expect_true(all(c("SE", "Infit", "Outfit") %in%
                    fit$config$attached_diagnostics_cols))
})

test_that("attach_diagnostics rejects non-logical inputs", {
  expect_error(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5, attach_diagnostics = "yes"),
    "attach_diagnostics"
  )
  expect_error(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5, attach_diagnostics = NA),
    "attach_diagnostics"
  )
  expect_error(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 5,
             attach_diagnostics = c(TRUE, FALSE)),
    "attach_diagnostics"
  )
})

test_that("attach_diagnostics_to_fit() is idempotent when called twice", {
  fit1 <- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 15, attach_diagnostics = TRUE)
  ))
  # Call internal helper a second time; column count should not grow.
  attach_fn <- getFromNamespace("attach_diagnostics_to_fit", "mfrmr")
  fit2 <- attach_fn(fit1)
  expect_equal(ncol(fit1$facets$others), ncol(fit2$facets$others))
  expect_setequal(names(fit1$facets$others), names(fit2$facets$others))
})

test_that("attach_diagnostics on MML fit also populates the canonical columns", {
  fit_mml <- suppressMessages(suppressWarnings(
    fit_mfrm(.toy, "Person", c("Rater", "Criterion"), "Score",
             method = "MML", quad_points = 7, maxit = 50,
             attach_diagnostics = TRUE)
  ))
  # MML fits have SE = ModelSE from observed information; must be finite.
  expect_true(any(is.finite(fit_mml$facets$others$SE)))
  expect_true(all(c("SE", "Infit", "Outfit") %in%
                    names(fit_mml$facets$others)))
})
