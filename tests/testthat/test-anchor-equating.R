# Tests for Phase 4: Anchoring & Equating Workflow

# ---------- shared fixtures (computed once) ----------
d1   <- load_mfrmr_data("study1")
d2   <- load_mfrmr_data("study2")
fit1 <- fit_mfrm(d1, person = "Person", facets = c("Rater", "Criterion"),
                 score = "Score", method = "JML")
fit2 <- fit_mfrm(d2, person = "Person", facets = c("Rater", "Criterion"),
                 score = "Score", method = "JML")

# ================================================================
# anchor_to_baseline
# ================================================================

test_that("anchor_to_baseline returns correct class and structure", {
  res <- suppressWarnings(
    anchor_to_baseline(d2, fit1, person = "Person",
                       facets = c("Rater", "Criterion"),
                       score = "Score")
  )

  expect_s3_class(res, "mfrm_anchored_fit")
  expect_true(is.list(res))
  expect_named(res, c("fit", "diagnostics", "baseline_anchors", "drift"),
               ignore.order = TRUE)

  # fit is an mfrm_fit
  expect_s3_class(res$fit, "mfrm_fit")

  # baseline_anchors is a tibble with expected columns
  expect_true(is.data.frame(res$baseline_anchors))
  expect_true(all(c("Facet", "Level", "Anchor") %in% names(res$baseline_anchors)))
  expect_true(nrow(res$baseline_anchors) > 0)

  # drift is a tibble with expected columns
  expect_true(is.data.frame(res$drift))
  drift_cols <- c("Facet", "Level", "Baseline", "New", "Drift",
                  "SE_Baseline", "SE_New", "SE_Diff", "Drift_SE_Ratio", "Flag")
  expect_true(all(drift_cols %in% names(res$drift)))
})

test_that("anchor_to_baseline self-anchoring yields near-zero drift", {
  # Anchor fit1 data to fit1 itself -> drift should be ~0
  res <- anchor_to_baseline(d1, fit1, person = "Person",
                            facets = c("Rater", "Criterion"),
                            score = "Score")

  expect_s3_class(res, "mfrm_anchored_fit")

  # All drifts should be very small (< 0.1 logits)
  if (nrow(res$drift) > 0) {
    expect_true(all(abs(res$drift$Drift) < 0.1),
                info = "Self-anchored drift should be near zero")
  }
})

test_that("anchor_to_baseline rejects non-mfrm_fit input", {
  expect_error(
    anchor_to_baseline(data.frame(), list(x = 1), "P", "F", "S"),
    "mfrm_fit"
  )
})

test_that("fit_mfrm surfaces malformed anchor schemas instead of silently dropping them", {
  toy <- load_mfrmr_data("example_core")
  bad_anchors <- data.frame(
    WrongFacet = "Rater",
    WrongLevel = "R1",
    WrongValue = 0,
    stringsAsFactors = FALSE
  )

  expect_warning(
    withCallingHandlers(
      fit_mfrm(
        toy,
        person = "Person",
        facets = c("Rater", "Criterion"),
        score = "Score",
        method = "JML",
        maxit = 15,
        anchors = bad_anchors,
        anchor_policy = "warn"
      ),
      warning = function(w) {
        if (grepl("Optimizer did not fully converge", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    "anchor_schema_mismatch"
  )

  expect_error(
    fit_mfrm(
      toy,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      method = "JML",
      maxit = 15,
      anchors = bad_anchors,
      anchor_policy = "error"
    ),
    "anchor_schema_mismatch"
  )
})

test_that("anchor_to_baseline S3 methods produce output", {
  res <- suppressWarnings(
    anchor_to_baseline(d2, fit1, person = "Person",
                       facets = c("Rater", "Criterion"),
                       score = "Score")
  )

  # summary returns expected class
  s <- summary(res)
  expect_s3_class(s, "summary.mfrm_anchored_fit")
  expect_true(is.numeric(s$n_anchored))
  expect_true(is.numeric(s$n_common))
  expect_true(is.numeric(s$n_flagged))

  # print methods produce output without error
  expect_output(print(res), "Anchored Fit Summary")
  expect_output(print(s), "Anchored Fit Summary")
})

# ================================================================
# detect_anchor_drift
# ================================================================

test_that("detect_anchor_drift returns correct class and structure", {
  drift <- detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2))

  expect_s3_class(drift, "mfrm_anchor_drift")
  expect_named(drift, c("drift_table", "summary", "common_elements", "common_by_facet", "config"),
               ignore.order = TRUE)

  # drift_table is a tibble with expected columns
  expect_true(is.data.frame(drift$drift_table))
  dt_cols <- c("Facet", "Level", "Reference", "Wave",
               "Ref_Est", "Wave_Est", "LinkOffset", "Drift", "SE_Ref", "SE_Wave", "SE",
               "Drift_SE_Ratio", "LinkSupportAdequate", "Flag")
  expect_true(all(dt_cols %in% names(drift$drift_table)))

  # common_elements has expected columns
  expect_true(is.data.frame(drift$common_elements))
  expect_true(all(c("Wave1", "Wave2", "N_Common") %in% names(drift$common_elements)))
  expect_true(is.data.frame(drift$common_by_facet))
  expect_true(all(c("Reference", "Wave", "Facet", "N_Common", "N_Retained",
                    "GuidelineMinCommon", "LinkSupportAdequate") %in% names(drift$common_by_facet)))

  # config preserves settings
  expect_equal(drift$config$reference, "Wave1")
  expect_equal(drift$config$method, "screened_common_element_alignment")
  expect_equal(drift$config$intended_use, "review_screen")
  expect_equal(drift$config$drift_threshold, 0.5)
  expect_equal(drift$config$min_common_per_facet, 5L)
  expect_equal(drift$config$waves, c("Wave1", "Wave2"))
})

test_that("detect_anchor_drift finds common elements", {
  drift <- detect_anchor_drift(list(W1 = fit1, W2 = fit2))

  # Should have at least some common elements
  expect_true(nrow(drift$common_elements) > 0)
  expect_true(all(drift$common_elements$N_Common >= 0))
})

test_that("detect_anchor_drift uses aligned drift and combined standard errors", {
  drift <- detect_anchor_drift(list(W1 = fit1, W2 = fit2))

  if (nrow(drift$drift_table) > 0) {
    expected_se <- sqrt(drift$drift_table$SE_Ref^2 + drift$drift_table$SE_Wave^2)
    expect_equal(drift$drift_table$SE, expected_se, tolerance = 1e-8)
    expect_equal(
      drift$drift_table$Drift_SE_Ratio,
      abs(drift$drift_table$Drift) / drift$drift_table$SE,
      tolerance = 1e-8
    )
  }
})

test_that("detect_anchor_drift warns when retained link support is thin", {
  d_small1 <- simulate_mfrm_data(
    n_person = 20,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    seed = 901
  )
  d_small2 <- simulate_mfrm_data(
    n_person = 20,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    seed = 902
  )
  fit_small1 <- suppressWarnings(
    fit_mfrm(d_small1, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 10)
  )
  fit_small2 <- suppressWarnings(
    fit_mfrm(d_small2, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 10)
  )

  expect_warning(
    drift <- detect_anchor_drift(list(W1 = fit_small1, W2 = fit_small2), facets = "Rater"),
    "Thin linking support"
  )

  expect_true(any(!drift$common_by_facet$LinkSupportAdequate))
})

test_that("detect_anchor_drift flagging logic works", {
  # Use a very small threshold to trigger flags
  drift <- detect_anchor_drift(list(W1 = fit1, W2 = fit2),
                               drift_threshold = 0.01,
                               flag_se_ratio = 0.01)

  # With such small thresholds, most elements should be flagged
  if (nrow(drift$drift_table) > 0) {
    expect_true(is.logical(drift$drift_table$Flag))
  }

  # Use a very large threshold to suppress flags
  drift_lax <- detect_anchor_drift(list(W1 = fit1, W2 = fit2),
                                   drift_threshold = 100,
                                   flag_se_ratio = 100)
  if (nrow(drift_lax$drift_table) > 0) {
    expect_equal(sum(drift_lax$drift_table$Flag), 0)
  }
})

test_that("detect_anchor_drift rejects invalid input", {
  expect_error(detect_anchor_drift(list()), "length")
  expect_error(detect_anchor_drift(list(a = 1, b = 2)), "mfrm_fit")
})

test_that("detect_anchor_drift S3 methods produce output", {
  drift <- detect_anchor_drift(list(W1 = fit1, W2 = fit2))

  s <- summary(drift)
  expect_s3_class(s, "summary.mfrm_anchor_drift")
  expect_true(is.numeric(s$n_comparisons))
  expect_true(is.numeric(s$n_flagged))

  expect_output(print(drift), "Anchor Drift Screen")
  expect_output(print(s), "Anchor Drift Screen")
})

# ================================================================
# build_equating_chain
# ================================================================

test_that("build_equating_chain returns correct class and structure", {
  chain <- build_equating_chain(list(Form1 = fit1, Form2 = fit2))

  expect_s3_class(chain, "mfrm_equating_chain")
  expect_named(chain, c("links", "cumulative", "element_detail", "common_by_facet", "config"),
               ignore.order = TRUE)

  # links is a tibble with expected columns
  expect_true(is.data.frame(chain$links))
  link_cols <- c("Link", "From", "To", "N_Common", "N_Retained",
                 "Min_Common_Per_Facet", "Min_Retained_Per_Facet",
                 "Offset_Prelim", "Offset", "Offset_SD", "Max_Residual",
                 "LinkSupportAdequate",
                 "Offset_Method")
  expect_true(all(link_cols %in% names(chain$links)))
  expect_equal(nrow(chain$links), 1)  # 2 fits -> 1 link

  # cumulative has one row per wave
  expect_true(is.data.frame(chain$cumulative))
  expect_equal(nrow(chain$cumulative), 2)
  expect_true(all(c("Wave", "Cumulative_Offset") %in% names(chain$cumulative)))

  # First wave offset is always 0
  expect_equal(chain$cumulative$Cumulative_Offset[1], 0)
  expect_true(is.data.frame(chain$common_by_facet))
  expect_equal(chain$config$method, "screened_common_element_alignment")
  expect_equal(chain$config$intended_use, "screened_linking_aid")
})

test_that("build_equating_chain with 3 fits produces 2 links", {
  # Use fit1 three times (artificial but tests chain logic)
  chain <- build_equating_chain(list(A = fit1, B = fit2, C = fit1))

  expect_equal(nrow(chain$links), 2)
  expect_equal(nrow(chain$cumulative), 3)
  expect_equal(chain$cumulative$Wave, c("A", "B", "C"))

  # Cumulative offset of first wave is 0
  expect_equal(chain$cumulative$Cumulative_Offset[1], 0)
})

test_that("build_equating_chain uses inverse-variance weighted offsets", {
  chain <- build_equating_chain(list(F1 = fit1, F2 = fit2))
  detail <- chain$element_detail

  if (nrow(detail) > 0) {
    w <- 1 / (detail$SE_From^2 + detail$SE_To^2)
    keep <- is.finite(w) & detail$Retained
    expected_offset <- stats::weighted.mean(detail$Diff[keep], w = w[keep])
    expect_equal(chain$links$Offset[1], expected_offset, tolerance = 1e-8)
  }
})

test_that("build_equating_chain warns when retained link support is thin", {
  d_small1 <- simulate_mfrm_data(
    n_person = 20,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    seed = 903
  )
  d_small2 <- simulate_mfrm_data(
    n_person = 20,
    n_rater = 4,
    n_criterion = 3,
    raters_per_person = 2,
    seed = 904
  )
  fit_small1 <- suppressWarnings(
    fit_mfrm(d_small1, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 10)
  )
  fit_small2 <- suppressWarnings(
    fit_mfrm(d_small2, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 10)
  )

  expect_warning(
    chain <- build_equating_chain(list(F1 = fit_small1, F2 = fit_small2), anchor_facets = "Rater"),
    "Thin linking support"
  )

  expect_true(any(!chain$links$LinkSupportAdequate))
  expect_true(any(chain$common_by_facet$N_Retained < chain$config$min_common_per_facet))
})

test_that("build_equating_chain rejects invalid input", {
  expect_error(build_equating_chain(list()), "length")
  expect_error(build_equating_chain(list(a = 1, b = 2)), "mfrm_fit")
})

test_that("build_equating_chain S3 methods produce output", {
  chain <- build_equating_chain(list(F1 = fit1, F2 = fit2))

  s <- summary(chain)
  expect_s3_class(s, "summary.mfrm_equating_chain")
  expect_true(is.numeric(s$n_flagged))

  expect_output(print(chain), "Screened Linking Chain")
  expect_output(print(s), "Screened Linking Chain")
})

# ================================================================
# plot_anchor_drift
# ================================================================

# Precompute drift and chain objects for all plot tests
drift_obj <- detect_anchor_drift(list(W1 = fit1, W2 = fit2))
chain_obj <- build_equating_chain(list(F1 = fit1, F2 = fit2))

test_that("plot_anchor_drift drift type returns data with draw=FALSE", {
  result <- plot_anchor_drift(drift_obj, type = "drift", draw = FALSE)

  expect_s3_class(result, "mfrm_plot_data")
  expect_identical(result$data$plot, "drift")
  expect_true(is.data.frame(result$data$table))
  expect_true(nrow(result$data$table) > 0)
  expect_true(all(c("title", "subtitle", "legend", "reference_lines") %in% names(result$data)))
})

test_that("plot_anchor_drift heatmap type returns data with draw=FALSE", {
  result <- plot_anchor_drift(drift_obj, type = "heatmap", draw = FALSE)

  expect_s3_class(result, "mfrm_plot_data")
  expect_identical(result$data$plot, "heatmap")
  expect_true(is.matrix(result$data$matrix))
})

test_that("plot_anchor_drift chain type returns data with draw=FALSE", {
  result <- plot_anchor_drift(chain_obj, type = "chain", draw = FALSE)

  expect_s3_class(result, "mfrm_plot_data")
  expect_identical(result$data$plot, "chain")
  expect_true(is.data.frame(result$data$table))
  expect_true(all(c("Wave", "Cumulative_Offset") %in% names(result$data$table)))
})

test_that("plot_anchor_drift drift type draws without error", {
  pdf(NULL)  # suppress graphical output
  on.exit(dev.off(), add = TRUE)

  expect_no_error(plot_anchor_drift(drift_obj, type = "drift"))
})

test_that("plot_anchor_drift chain type draws without error", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_no_error(plot_anchor_drift(chain_obj, type = "chain"))
})

test_that("plot_anchor_drift heatmap type draws without error", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_no_error(plot_anchor_drift(drift_obj, type = "heatmap"))
})

test_that("plot_anchor_drift accepts publication preset", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_no_error(plot_anchor_drift(drift_obj, type = "drift", preset = "publication"))
  expect_no_error(plot_anchor_drift(chain_obj, type = "chain", preset = "publication"))
})

test_that("plot_anchor_drift rejects unsupported type/class combo", {
  # chain object with drift type should error
  expect_error(plot_anchor_drift(chain_obj, type = "drift"),
               "Unsupported")
})

test_that("plot_anchor_drift facet filter works", {
  result <- plot_anchor_drift(drift_obj, type = "drift", facet = "Rater",
                              draw = FALSE)

  if (inherits(result, "mfrm_plot_data") && nrow(result$data$table) > 0) {
    expect_true(all(result$data$table$Facet == "Rater"))
  }
})
