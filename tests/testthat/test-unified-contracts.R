test_that("estimation has one class and rejects removed controls and old saved schemas", {
  d <- load_mfrmr_data("example_core")
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", method = "JML")
  expect_identical(class(f), "mfrm_fit")
  expect_identical(class(diagnose_mfrm(f)), "mfrm_diagnostics")
  expect_false("fast" %in% names(formals(fit_mfrm)))
  expect_false(any(grepl("mfrm_fast_", names(asNamespace("mfrmr")), fixed = TRUE)))
  expect_error(fit_mfrm(d, fast = TRUE), "Unknown argument.*fast")
  expect_error(fit_mfrm(d, control = list(fast = FALSE)), "Unsupported.*fast")
  for (option in c("jml_correction_order", "fit_scope", "mml_engine", "attach_diagnostics")) {
    expect_error(do.call(fit_mfrm, c(list(data = d), setNames(list(TRUE), option))), "Unknown argument")
  }
  old <- f; old$schema_version <- NULL
  for (fun in list(summary, coef, vcov, confint, fitted, predict, plot, diagnose_mfrm))
    expect_error(fun(old), "Refit older saved objects")
  restored <- unserialize(serialize(f, NULL))
  expect_equal(predict(restored, type = "probs"), predict(f, type = "probs"))
  plan <- unique(d[c("Rater", "Criterion")])
  fair <- fair_average_table(f, reference_plan = plan)
  ci <- confint(fair)
  expect_equal(ci$Estimate, fair$reference_scores$Estimate)
  expect_true(all(ci$Lower <= ci$Estimate & ci$Estimate <= ci$Upper))
  expect_true(all(ci$Lower >= 1 & ci$Upper <= 4))
  expect_true(all(grepl("conditional", ci$Uncertainty)))
})

test_that("all G-study estimators use one contract while retaining statistical assumptions", {
  d <- load_mfrmr_data("example_core")
  # A balanced subset is sufficient for the three estimator/S3 contracts.
  d <- d[d$Person %in% unique(d$Person)[1:8] &
           d$Rater %in% unique(d$Rater)[1:3] &
           d$Criterion %in% unique(d$Criterion)[1:2], ]
  results <- lapply(c("anova", "minque0", "urgenova"), function(method)
    mfrm_g_study(d, score = "Score", person = "Person", facets = c("Rater", "Criterion"), method = method))
  for (g in results) {
    expect_identical(class(g), "mfrm_gtheory")
    expect_identical(g$schema_version, 1L)
    expect_true(all(c("components", "variance_components", "design", "method_code", "g_study_summary", "d_study") %in% names(g)))
    expect_s3_class(tidy(g), "data.frame")
    expect_s3_class(plot(g), "ggplot")
    expect_equal(summary(g)$variance_components, g$variance_components)
  }
  expect_equal(results[[1]]$variance_components$RawEstimate,
    results[[3]]$variance_components$RawEstimate, tolerance = 1e-7)
  expect_error(confint(results[[1]]), "assumption")
  ci <- confint(results[[1]], assumption = "normal", method = "mls")
  expect_true(all(is.finite(ci$Upper)))
  expect_error(confint(results[[3]], method = "mls"), "must be empty")
  incomplete <- mfrm_g_study(d[-1, ], score = "Score", person = "Person",
    facets = c("Rater", "Criterion"), method = "minque0")
  expect_equal(nrow(incomplete$d_study), 0L)
  expect_error(mfrm_d_study(incomplete), "explicit.*design_grid")
})

test_that("automatic quadrature preserves refined results and withholds unstable inference", {
  d <- load_mfrmr_data("example_core")
  expect_warning(coarse <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
    control = list(mml_integration = "auto", quad_points = 15L)), "quadrature orders")
  expect_equal(coarse$quadrature_check$orders, c(15L, 31L))
  expect_equal(coarse$control$quad_points, 31L)
  expect_false(summary(coarse)$readiness$ready)
  expect_true("quadrature_sensitivity_not_stable" %in% summary(coarse)$readiness$reason_codes)
  expect_error(confint(coarse), "quadrature_sensitivity")
})

test_that("fine automatic quadrature is stable and survives serialization", {
  skip_on_cran()
  d <- load_mfrmr_data("example_core")
  fine <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
    control = list(mml_integration = "auto", quad_points = 61L))
  expect_equal(fine$quadrature_check$orders, c(61L, 123L))
  expect_true(fine$quadrature_check$stable)
  expect_true(summary(fine)$readiness$ready)
  expect_identical(unserialize(serialize(fine, NULL))$quadrature_check, fine$quadrature_check)
})

test_that("data-review fit recipes use current estimation arguments", {
  d <- load_mfrmr_data("example_core")
  levels <- sort(unique(as.character(d$Rater)))
  W <- matrix(c(1, -1, rep(0, length(levels) - 2L)), nrow = 1L,
    dimnames = list("first_minus_second", levels))
  review <- describe_mfrm_data(d, "Person", c("Rater", "Criterion", "Group"), "Score",
    category_policy = "preserve", include_agreement = FALSE,
    comparisons = list(facet = "Rater", group = "Group", groups = c("A", "B"), contrasts = W))
  recipe <- review$comparison_review$fit_arguments
  f <- do.call(fit_mfrm, c(list(data = d), recipe))
  expect_identical(class(f), "mfrm_fit")
  expect_equal(nrow(predict(f, type = "probs")), nrow(d))
  expect_true(length(f$encoding$interaction_meta) > 0L)
})

test_that("anchor displacement summaries follow the unified calibration contract", {
  d <- load_mfrmr_data("example_core")
  f <- read_example_bundle()$fit
  # Identical estimates isolate table matching and interval arithmetic, not wave independence.
  x <- detect_anchor_drift(list(A = f, B = f), include_person = TRUE)
  s <- summary(unserialize(serialize(x, NULL)))
  expect_equal(s$n_comparisons, nrow(f$facets) + nrow(f$person))
  expect_equal(s$table$Drift, rep(0, s$n_comparisons))
  expect_equal(s$table$SE, sqrt(2) * c(f$facets$ModelSE, f$person$ModelSE))
  expect_equal(s$n_flagged_available, 0L)
  expect_equal(nrow(s$flagged), 0L)
  empty <- detect_anchor_drift(list(f, f), facets = "absent")
  expect_identical(summary(empty)$n_comparisons, 0L)
  expect_identical(summary(empty)$n_flagged, 0L)
  expect_error(detect_anchor_drift(list(f, f), reference = 1.2), "reference")
  expect_error(detect_anchor_drift(list(A = f, B = f), reference = "missing"), "reference")
  old <- f; old$schema_version <- NULL
  expect_error(detect_anchor_drift(list(f, old)), "Refit older saved objects")
  held <- f; held$converged <- FALSE
  unavailable <- detect_anchor_drift(list(f, held))
  expect_true(all(is.na(unavailable$table$Flagged)))
  expect_true(all(unavailable$table$Status == "inference_unavailable"))
  expect_equal(summary(unavailable)$n_unavailable, nrow(f$facets))
})

test_that("retired screening controls cannot silently change a requested analysis", {
  d <- load_mfrmr_data("example_core")
  f <- read_example_bundle()$fit
  expect_error(estimate_bias(f, max_iter = 2), "unused argument")
  expect_error(estimate_bias(f, diagnostics = diagnose_mfrm(f)), "unused argument")
  expect_error(analyze_dff(f, "Rater", "Group", method = "refit"), "arg")
  expect_error(analyze_dff(f, "Rater", "Group", comparison_review = list()), "unused argument")
  expect_error(compare_mfrm(f, f, response_diagnostics = list()), "current fit_mfrm")
  old <- f; old$schema_version <- NULL
  expect_error(compare_mfrm(f, old), "Refit older saved objects")
  expect_no_warning(cmp <- compare_mfrm(f, f, nested = TRUE, warn_constraints = FALSE))
  expect_identical(cmp$comparison_basis$lrt_status, "ineligible")
  expect_null(cmp$lrt)
})
