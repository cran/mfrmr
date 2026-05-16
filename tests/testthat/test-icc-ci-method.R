# Tests for compute_facet_icc(ci_method = ...) and its passthrough in
# analyze_hierarchical_structure(). Targets the new CI paths added in
# 0.1.6 (profile + parametric bootstrap) and the canonical ci_method
# argument rename (0.1.6 deprecation of icc_ci_method).

skip_if_not_installed("lme4")

local({
  .toy <<- load_mfrmr_data("example_core")
})

test_that("compute_facet_icc(ci_method = 'none') is backward-compatible", {
  icc <- compute_facet_icc(.toy, facets = c("Rater", "Criterion"),
                           score = "Score", person = "Person")
  expect_s3_class(icc, "mfrm_facet_icc")
  expect_true(all(c("Facet", "Variance", "ICC", "Interpretation",
                    "InterpretationScale") %in% names(icc)))
  # CI columns exist but are NA.
  expect_true(all(c("ICC_CI_Lower", "ICC_CI_Upper",
                    "ICC_CI_Method") %in% names(icc)))
  expect_true(all(is.na(icc$ICC_CI_Lower)))
  expect_true(all(is.na(icc$ICC_CI_Upper)))
  expect_true(all(icc$ICC_CI_Method == "none"))
})

test_that("ci_method = 'profile' populates CI bounds that bracket the estimate", {
  icc <- compute_facet_icc(.toy, facets = c("Rater", "Criterion"),
                           score = "Score", person = "Person",
                           ci_method = "profile", ci_level = 0.95)
  expect_true(all(icc$ICC_CI_Method == "profile"))
  valid <- is.finite(icc$ICC_CI_Lower) & is.finite(icc$ICC_CI_Upper)
  expect_true(any(valid))
  # The point estimate should fall (approximately) inside its own CI
  # for valid rows.
  expect_true(all(
    icc$ICC_CI_Lower[valid] - 1e-4 <= icc$ICC[valid] &
      icc$ICC[valid] <= icc$ICC_CI_Upper[valid] + 1e-4
  ))
})

test_that("ci_method = 'boot' reports replicate count and quantile CI", {
  icc <- compute_facet_icc(.toy, facets = c("Rater", "Criterion"),
                           score = "Score", person = "Person",
                           ci_method = "boot",
                           ci_boot_reps = 50L, ci_boot_seed = 123L)
  expect_true(all(icc$ICC_CI_Method == "boot"))
  expect_true("ICC_CI_NReps" %in% names(icc))
  expect_true(all(icc$ICC_CI_NReps <= 50L & icc$ICC_CI_NReps >= 1L))
  # Bootstrap CI widths should be positive when valid.
  valid <- is.finite(icc$ICC_CI_Lower) & is.finite(icc$ICC_CI_Upper)
  expect_true(any(valid))
  expect_true(all(icc$ICC_CI_Upper[valid] - icc$ICC_CI_Lower[valid] >= 0))
})

test_that("ci_method is validated", {
  expect_error(
    compute_facet_icc(.toy, facets = c("Rater", "Criterion"),
                      score = "Score", person = "Person",
                      ci_method = "nonsense"),
    "arg"
  )
})

test_that("ci_level bounds are validated", {
  expect_error(
    compute_facet_icc(.toy, facets = c("Rater", "Criterion"),
                      score = "Score", person = "Person",
                      ci_method = "profile", ci_level = 1.2),
    "ci_level"
  )
  expect_error(
    compute_facet_icc(.toy, facets = c("Rater", "Criterion"),
                      score = "Score", person = "Person",
                      ci_method = "profile", ci_level = 0),
    "ci_level"
  )
})

test_that("ci_boot_seed makes bootstrap CI reproducible", {
  icc_a <- compute_facet_icc(.toy, facets = c("Rater", "Criterion"),
                              score = "Score", person = "Person",
                              ci_method = "boot",
                              ci_boot_reps = 40L, ci_boot_seed = 2026L)
  icc_b <- compute_facet_icc(.toy, facets = c("Rater", "Criterion"),
                              score = "Score", person = "Person",
                              ci_method = "boot",
                              ci_boot_reps = 40L, ci_boot_seed = 2026L)
  expect_equal(icc_a$ICC_CI_Lower, icc_b$ICC_CI_Lower, tolerance = 1e-8)
  expect_equal(icc_a$ICC_CI_Upper, icc_b$ICC_CI_Upper, tolerance = 1e-8)
})

test_that("analyze_hierarchical_structure passes ci_method through to ICC", {
  h <- suppressMessages(suppressWarnings(analyze_hierarchical_structure(
    .toy, facets = c("Rater", "Criterion"),
    person = "Person", score = "Score",
    ci_method = "profile"
  )))
  expect_true("ICC_CI_Method" %in% names(h$icc))
  expect_true(all(h$icc$ICC_CI_Method == "profile"))
  expect_true(any(is.finite(h$icc$ICC_CI_Lower)))
})

test_that("deprecated icc_ci_method still works with a warning", {
  suppressMessages(suppressWarnings({
    w <- testthat::capture_warnings(
      h <- analyze_hierarchical_structure(
        .toy, facets = c("Rater", "Criterion"),
        person = "Person", score = "Score",
        icc_ci_method = "profile"
      )
    )
  }))
  expect_true(any(grepl("icc_ci_method", paste(w, collapse = " "))))
  expect_true(any(is.finite(h$icc$ICC_CI_Lower)))
})

test_that("plot.mfrm_hierarchical_structure(type = 'icc') renders CI whiskers", {
  h <- suppressMessages(suppressWarnings(analyze_hierarchical_structure(
    .toy, facets = c("Rater", "Criterion"),
    person = "Person", score = "Score",
    ci_method = "profile"
  )))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(plot(h, type = "icc")))
})
