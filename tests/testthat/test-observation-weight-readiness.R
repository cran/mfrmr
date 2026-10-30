test_that("observation weights restrict inference without changing unit-weight estimation", {
  toy <- load_mfrmr_data("example_core")
  fit_case <- function(data, model, weight = NULL) fit_mfrm(
    data, "Person", c("Rater", "Criterion"), "Score", weight = weight,
    model = model, step_facet = if (model == "PCM") "Criterion" else NULL,
    method = "MML", maxit = 150
  )
  for (model in c("RSM", "PCM")) {
    baseline <- fit_case(toy, model)
    unit <- toy
    unit$Weight <- 1
    explicit <- fit_case(unit, model, "Weight")
    expect_true(mfrm_inference_ready(baseline))
    expect_true(mfrm_inference_ready(explicit))
    expect_equal(explicit$opt$par, baseline$opt$par, tolerance = 1e-12)
    expect_equal(compute_mml_parameter_covariance(explicit)$cov,
                 compute_mml_parameter_covariance(baseline)$cov, tolerance = 1e-12)
    diag_unit <- diagnose_mfrm(explicit, residual_pca = "none")
    expect_true(diag_unit$precision_profile$SupportsFormalInference)
    expect_true(analyze_facet_equivalence(explicit, diag_unit)$summary$InferenceReady)

    for (policy in c("nonunit_row_varying", "nonunit_constant_within_person")) {
      weighted <- toy
      index <- if (policy == "nonunit_row_varying") toy$Rater else toy$Person
      weighted$Weight <- c(.5, 1, 1.5)[(match(index, unique(index)) - 1L) %% 3L + 1L]
      fit <- fit_case(weighted, model, "Weight")
      expect_identical(mfrm_ic_weight_policy(fit$prep, fit$config), policy)
      expect_identical(fit$readiness$fit$InputState, "review")
      expect_match(fit$readiness$fit$ReasonCodes,
                   "nonunit_observation_weights_inference_unvalidated")
      expect_false(mfrm_inference_ready(fit))
      expect_false(mfrm_inference_ready(fit$summary))
      diag <- diagnose_mfrm(fit, residual_pca = "none")
      expect_false(diag$precision_profile$SupportsFormalInference)
      expect_false(any(diag$measures$SupportsFormalInference))
      expect_false(any(diag$measures$CIEligible))
      expect_false(any(diag$parameter_uncertainty$steps$CIEligible))
      expect_false(any(diag$reliability$SupportsFormalInference))
      expect_true(all(is.finite(diag$measures$SE)))
      expect_true(all(is.finite(diag$parameter_uncertainty$steps$SE)))
      expect_match(diag$precision_profile$RecommendedUse, "diagnostic only", fixed = TRUE)
      attached <- attach_diagnostics_to_fit(fit)
      expect_false(any(attached$facets$others$SupportsFormalInference))
      expect_false(any(attached$facets$person$CIEligible))
      expect_false(any(attached$steps$CIEligible))
      expect_identical(summary(fit, diagnostics = diag)$decision$FormalInference, "No")
      expect_identical(summary(diag)$decision$FormalInference, "No")
      expect_match(summary(fit, diagnostics = diag)$decision$Why, "non-unit observation weights")
      expect_error(analyze_facet_equivalence(fit), "inference-ready MML")
      diag$precision_profile$SupportsFormalInference <- TRUE
      diag$measures$SupportsFormalInference <- TRUE
      expect_error(analyze_facet_equivalence(fit, diag), "inference-ready MML")
      expect_warning(plotted <- plot(fit, type = "wright", draw = FALSE), "Fit=review")
      expect_identical(plotted$data$interpretation_status, "review_only")
    }

    # Before the weight gate, saved v3 fits and detached diagnostics could say
    # TRUE. Unit objects of that version also require a current fit/re-audit.
    old <- explicit
    old$readiness$fit$ReadinessContractVersion <- "mfrmr-readiness-0.2.3-v3"
    old$summary$ReadinessContractVersion <- "mfrmr-readiness-0.2.3-v3"
    old <- unserialize(serialize(old, NULL))
    expect_false(mfrm_inference_ready(old))
    expect_false(mfrm_inference_ready(old$summary))
    expect_identical(summary(old)$decision$FormalInference, "No")
    expect_error(analyze_facet_equivalence(old, diag_unit), "inference-ready MML")
    diag_unit$fit_readiness$ReadinessContractVersion <- "mfrmr-readiness-0.2.3-v3"
    expect_identical(summary(diag_unit)$decision$FormalInference, "No")
    expect_false(summary(diag_unit)$precision_profile$SupportsFormalInference)
    expect_match(paste(capture.output(print(diag_unit)), collapse = "\n"), "legacy")
    expect_error(mfrm_results(explicit, diagnostics = diag_unit, compute = "never"),
                 "current readiness contract")
  }
})

test_that("the shared input audit fails closed for invalid or non-unit weights", {
  prep <- list(data = data.frame(Person = rep(c("P1", "P2"), each = 2), Weight = 1))
  for (w in list(rep(1, 4), rep(2, 4), c(1, 1, .5, .5), c(.5, 1.5, .5, 1.5),
                 c(1, NA, 1, 1), rep(0, 4), numeric(), NULL)) {
    candidate <- prep
    if (length(w) == 0L && !is.null(w)) candidate$data <- candidate$data[FALSE, ]
    candidate$data$Weight <- w
    audit <- mfrmr_readiness_input_component(candidate, list(), list(weight_col = "Weight"))
    expect_identical(audit$State, if (identical(w, rep(1, 4))) "pass" else
      if (is.null(w) || !length(w) || anyNA(w) || any(w <= 0)) "blocked" else "review")
  }
})
