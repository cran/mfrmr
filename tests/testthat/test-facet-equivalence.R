local({
  toy <- load_mfrmr_data("example_core")
  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "MML", maxit = 150)
  diagnostics <- diagnose_mfrm(fit, residual_pca = "none")

  test_that("equivalence uses joint covariance with an analytic normal reference", {
    # Four centered independent N(0, v) errors have V = v * (I - 11'/4).
    # Thus each pair has SE sqrt(2v), even though marginal SE = sqrt(3v/4).
    # The old diagonal-only calculation incorrectly declared the first pair
    # equivalent at bound 0.3; the correct one-sided p is about 0.0668.
    reference <- fit
    rater_rows <- reference$facets$others$Facet == "Rater"
    reference$facets$others$Estimate[rater_rows] <- c(0, 0, 0.2, -0.2)
    v <- 0.02
    covariance <- compute_mml_parameter_covariance(fit)
    slice <- covariance$param_slices$Rater
    covariance$cov[slice, slice] <- v * (diag(3) - 1 / 4)
    local_mocked_bindings(compute_mml_parameter_covariance = function(...) covariance)
    eq <- analyze_facet_equivalence(reference, facet = "Rater", equivalence_bound = 0.3)
    expect_equal(eq$pairwise$SE_Diff, rep(0.2, 6), tolerance = 1e-12)
    expect_equal(eq$pairwise$P_TOST[1], pnorm(-1.5), tolerance = 1e-12)
    expect_false(any(eq$pairwise$Equivalent))
    expect_equal(eq$chi_square$FixedChiSq, 4, tolerance = 1e-12)
    expect_equal(eq$chi_square$FixedDF, 3)
    expect_equal(eq$chi_square$FixedProb, pchisq(4, 3, lower.tail = FALSE))
    expect_equal(eq$rope$DeviationSE, rep(sqrt(0.015), 4))
    expect_true(is.na(eq$summary$BF01))
    expect_match(eq$summary$BF01Label, "requires fitted likelihood")
    expect_equal(eq$summary$MultiplicityAdjustment, "none")
  })

  test_that("eligible MML output and supplied diagnostics preserve the same basis", {
    eq <- analyze_facet_equivalence(fit, facet = "Rater")
    expect_s3_class(eq, "mfrm_facet_equivalence")
    expect_named(eq, c("summary", "chi_square", "pairwise", "rope", "forest", "settings"))
    expect_equal(eq, analyze_facet_equivalence(fit, diagnostics, facet = "Rater"))
    expect_equal(eq$summary, summary(eq)$summary)
    expect_true(eq$summary$InferenceReady)
    expect_true(all(eq$pairwise$Equivalent ==
      (eq$pairwise$CI90_Lower > -0.5 & eq$pairwise$CI90_Upper < 0.5)))
    forest <- plot_facet_equivalence(fit, diagnostics, facet = "Rater", draw = FALSE)
    rope <- plot(eq, type = "rope", draw = FALSE)
    expect_equal(forest$data, eq$forest)
    expect_equal(rope$data, eq$forest)
    expect_match(forest$title, "deviation from facet mean")
    printed <- paste(capture.output(print(summary(eq))), collapse = "\n")
    expect_match(printed, "two one-sided normal-reference tests")
    expect_false(grepl("mml_observed_information_contrasts|pairwise_tost_summary|InferenceReady|Decision", printed))
    expect_true(forest$inference_ready)
    expect_equal(forest$covariance_basis, eq$summary$CovarianceBasis)
    reordered <- diagnostics
    reordered$measures <- reordered$measures[nrow(reordered$measures):1, ]
    expect_equal(eq, analyze_facet_equivalence(fit, reordered, facet = "Rater"))
  })

  test_that("grand-mean proximity includes mean uncertainty with a fixed anchor", {
    reference <- fit
    labels <- reference$config$facet_specs$Rater$levels
    reference$config$facet_specs$Rater <- build_facet_constraint(
      labels, anchors = stats::setNames(0, labels[1]), centered = FALSE
    )
    reference$facets$others$Estimate[reference$facets$others$Facet == "Rater"] <-
      c(0, 0.1, 0.2, 0.3)
    covariance <- compute_mml_parameter_covariance(fit)
    slice <- covariance$param_slices$Rater
    covariance$cov[slice, slice] <- diag(0.02, 3)
    local_mocked_bindings(compute_mml_parameter_covariance = function(...) covariance)
    eq <- analyze_facet_equivalence(reference, facet = "Rater")
    # Independent errors (0, e2, e3, e4), Var(ej)=v: mean variance is 3v/16;
    # a free level minus the mean has variance (9+1+1)v/16.
    expect_equal(eq$summary$GrandMean, 0.15)
    expect_equal(eq$rope$SE, c(0, rep(sqrt(0.02), 3)))
    expect_identical(eq$rope$Fixed, c(TRUE, FALSE, FALSE, FALSE))
    expect_true(is.na(eq$rope$CI_Lower[1]) && is.na(eq$rope$CI_Upper[1]))
    expect_equal(eq$rope$DeviationSE, sqrt(c(3, 11, 11, 11) * 0.02 / 16))
    expect_equal(eq$pairwise$SE_Diff, sqrt(c(rep(0.02, 3), rep(0.04, 3))))
    expect_equal(nrow(eq$pairwise), 6)
    expect_equal(eq$rope$ROPEPct[1], 100 * (
      pnorm(0.5, -0.15, sqrt(3 * 0.02 / 16)) -
        pnorm(-0.5, -0.15, sqrt(3 * 0.02 / 16))
    ))
  })

  test_that("actual anchored MML fits agree with and without supplied diagnostics", {
    anchored <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
      method = "MML", noncenter_facet = "Rater",
      anchors = data.frame(Facet = "Rater", Level = "R01", Anchor = .1))
    dx <- diagnose_mfrm(anchored)
    eq <- analyze_facet_equivalence(anchored, facet = "Rater")
    expect_identical(eq, analyze_facet_equivalence(anchored, dx, facet = "Rater"))
    fixed <- eq$forest$Element == "R01"
    expect_true(eq$forest$Fixed[fixed])
    expect_equal(eq$forest$Measure[fixed], .1)
    expect_equal(eq$forest$SE[fixed], 0)
    expect_true(is.na(eq$forest$CI_Lower[fixed]) && is.na(eq$forest$CI_Upper[fixed]))
    expect_gt(eq$forest$DeviationSE[fixed], 0)
    expect_true(all(is.finite(eq$pairwise$SE_Diff)))
    expect_identical(plot(eq, draw = FALSE)$data, eq$forest)
    bad <- dx; bad$measures$Fixed[bad$measures$Facet == "Rater"] <- FALSE
    expect_error(analyze_facet_equivalence(anchored, bad, facet = "Rater"), "fixed-value labels")
  })

  test_that("constrained contrasts are rejected before covariance roundoff can admit them", {
    reference <- fit
    labels <- reference$config$facet_specs$Rater$levels
    reference$config$facet_specs$Rater <- build_facet_constraint(
      labels, anchors = stats::setNames(0.25, labels[1]), centered = TRUE
    )
    reference$facets$others$Estimate[reference$facets$others$Facet == "Rater"] <-
      c(0.25, 0.1, 0.2, -0.3)
    # One fixed level plus three centered free levels gives only two free
    # directions for three joint contrasts, even with nonsingular free V.
    covariance <- compute_mml_parameter_covariance(fit)
    covariance$param_slices$Rater <- 1:2
    covariance$cov <- diag(c(0.02, 0.03))
    local_mocked_bindings(compute_mml_parameter_covariance = function(...) covariance)
    expect_error(analyze_facet_equivalence(reference, facet = "Rater"),
                 "independent under the model constraints")
    expect_error(plot_facet_equivalence(reference, facet = "Rater", draw = FALSE),
                 "independent under the model constraints")
  })

  test_that("ineligible fits and covariance cannot produce equivalence decisions", {
    for (method in c("MML", "JML")) {
      for (model in c("RSM", "PCM", "GPCM")) {
        candidate <- if (method == "MML" && model == "RSM") fit else fit_mfrm(
          toy, "Person", c("Rater", "Criterion"), "Score",
          method = method, model = model, maxit = 150,
          step_facet = if (model == "RSM") NULL else "Criterion",
          slope_facet = if (model == "GPCM") "Criterion" else NULL
        )
        if (method == "MML" && model != "GPCM") {
          expect_true(analyze_facet_equivalence(candidate)$summary$InferenceReady)
        } else {
          for (bound in c(0.5, 5)) {
            expect_error(analyze_facet_equivalence(candidate, equivalence_bound = bound),
                         "requires an inference-ready MML fit")
          }
          expect_error(plot_facet_equivalence(candidate, draw = FALSE),
                       "requires an inference-ready MML fit")
        }
      }
    }
    unready <- fit
    unready$readiness$fit$InferenceReady <- FALSE
    expect_error(analyze_facet_equivalence(unready, diagnostics), "inference-ready MML")
    covariance <- compute_mml_parameter_covariance(fit)
    local_mocked_bindings(compute_mml_parameter_covariance = function(...) covariance)
    for (status in c("regularized", "fallback", "not_applicable")) {
      covariance$status <- status
      expect_error(analyze_facet_equivalence(fit), "unregularized MML")
    }
    covariance$status <- "ok"
    covariance$cov[,] <- 0
    expect_error(analyze_facet_equivalence(fit), "positive-definite covariance")
  })

  test_that("ROPE mass is symmetric and retains narrow normal bands", {
    reference <- fit
    reference$facets$others$Estimate[reference$facets$others$Facet == "Rater"] <- c(-10, 10, 0, 0)
    covariance <- compute_mml_parameter_covariance(fit)
    slice <- covariance$param_slices$Rater
    covariance$cov[slice, slice] <- (diag(3) - 1 / 4) / .75
    local_mocked_bindings(compute_mml_parameter_covariance = function(...) covariance)
    eq <- analyze_facet_equivalence(reference, facet = "Rater", equivalence_bound = .1)
    expected <- 100 * integrate(dnorm, 9.9, 10.1, abs.tol = 0, rel.tol = 1e-10)$value
    expect_gt(eq$rope$ROPEPct[1], 0)
    expect_equal(eq$rope$ROPEPct[1:2] / expected, c(1, 1), tolerance = 1e-12)
    narrow <- analyze_facet_equivalence(reference, facet = "Rater", equivalence_bound = 1e-20)
    expected <- 100 * 2e-20 * dnorm(c(10, 10, 0, 0))
    expect_equal(narrow$rope$ROPEPct / expected, rep(1, 4), tolerance = 1e-12)
  })

  test_that("unresolved display intervals retain SEs and independent TOST decisions", {
    covariance <- compute_mml_parameter_covariance(fit)
    local_mocked_bindings(compute_mml_parameter_covariance = function(...) covariance)
    ordinary <- analyze_facet_equivalence(fit, facet = "Rater")
    narrow <- analyze_facet_equivalence(fit, facet = "Rater", ci_level = 1e-20)
    expect_identical(narrow$pairwise, ordinary$pairwise)
    expect_identical(narrow$rope$DeviationSE, ordinary$rope$DeviationSE)
    expect_true(all(is.na(narrow$rope$DeviationCI_Lower)))
    expect_true(all(is.na(narrow$rope$DeviationCI_Upper)))
    expect_true(all(narrow$rope$DeviationCIReason == "numerical_resolution"))
    expect_true(all(narrow$rope$CIReason == "numerical_resolution"))
    expect_true(all(narrow$rope$ROPEStatus == "unavailable"))
    expect_identical(narrow$summary$Decision, ordinary$summary$Decision)
    expect_match(paste(capture.output(print(summary(narrow))), collapse = " "), "numeric precision")
    payload <- plot(narrow, draw = FALSE)
    expect_identical(payload$data, narrow$forest)
    expect_match(payload$note, "numeric precision")
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
    expect_no_warning(plot(narrow))
  })

  test_that("mismatched or incomplete diagnostics cannot override the fit", {
    changed <- diagnostics
    rows <- which(changed$measures$Facet == "Rater")
    changed$measures$Estimate[rows[1]] <- changed$measures$Estimate[rows[1]] + 0.1
    expect_error(analyze_facet_equivalence(fit, changed), "fit/diagnostics mismatch")
    changed <- diagnostics
    changed$measures$SE[rows] <- 1e-6
    expect_error(analyze_facet_equivalence(fit, changed), "matching model-based SEs")
    changed <- diagnostics
    changed$measures <- changed$measures[-rows[1], ]
    expect_error(analyze_facet_equivalence(fit, changed), "all levels")
    changed <- diagnostics
    changed$measures$SupportsFormalInference[rows[1]] <- FALSE
    expect_error(analyze_facet_equivalence(fit, changed), "ordinary-inference eligibility")
    changed <- diagnostics
    changed$precision_profile$SupportsFormalInference <- FALSE
    expect_error(analyze_facet_equivalence(fit, changed), "do not support ordinary inference")
    changed <- diagnostics
    changed$measures$SupportsFormalInference <- NULL
    expect_error(analyze_facet_equivalence(fit, changed), "ordinary-inference eligibility")
    expect_error(analyze_facet_equivalence(fit, list(measures = diagnostics$measures)),
                 "mfrm_diagnostics object")
  })

  test_that("legacy and ineligible bundles are rejected by summary print and plots", {
    eq <- analyze_facet_equivalence(fit)
    for (kind in c("legacy", "ineligible", "missing_rank", "insufficient_rank", "old_readiness")) {
      changed <- eq
      if (kind == "legacy") changed$settings$covariance_basis <- NULL
      if (kind == "ineligible") changed$summary$InferenceReady <- FALSE
      if (kind == "missing_rank") changed$settings$contrast_rank <- NULL
      if (kind == "insufficient_rank") changed$settings$contrast_rank <- 1L
      if (kind == "old_readiness") changed$settings$readiness_contract_version <- NULL
      expect_error(summary(changed), "Recompute it")
      expect_error(print(changed), "Recompute it")
      expect_error(plot(changed, draw = FALSE), "Recompute it")
      expect_error(plot_facet_equivalence(changed, type = "rope", draw = FALSE), "Recompute it")
    }
  })
})
