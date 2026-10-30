local({
  toy <- load_mfrmr_data("example_core")
  fits <- diagnostics <- list()
  for (method in c("MML", "JML")) {
    for (spec in c("RSM", "PCM", "GPCM_Criterion", "GPCM_Rater")) {
      model <- sub("_.*", "", spec)
      owner <- if (spec == "GPCM_Rater") "Rater" else "Criterion"
      key <- paste(method, spec, sep = "_")
      fits[[key]] <- fit_mfrm(
        toy, "Person", c("Rater", "Criterion"), "Score",
        method = method, model = model, maxit = 150,
        step_facet = if (model == "RSM") NULL else owner,
        slope_facet = if (model == "GPCM") owner else NULL
      )
      diagnostics[[key]] <- diagnose_mfrm(fits[[key]], residual_pca = "none")
    }
  }

  test_that("matched saved diagnostics retain estimator-specific precision", {
    for (key in names(fits)) {
      fit <- fits[[key]]
      file <- tempfile(fileext = ".rds")
      saveRDS(diagnostics[[key]], file)
      dx <- readRDS(file)
      unlink(file)
      formal <- key %in% c("MML_RSM", "MML_PCM")
      brief <- summary(fit, diagnostics = dx)
      precision <- precision_review_report(fit, dx)
      expect_identical(brief$decision$FormalInference, if (formal) "Yes" else "No", info = key)
      expect_identical(precision$profile, as.data.frame(dx$precision_profile), info = key)
      if (startsWith(key, "JML")) {
        expect_identical(precision$profile$PrecisionTier, "exploratory", info = key)
        expect_false(any(dx$measures$SupportsFormalInference), info = key)
      }
      if (grepl("GPCM", key)) {
        expect_false(any(fit$slopes$SEEligible), info = key)
        expect_false(any(fit$slopes$CIEligible), info = key)
      }
      apa <- apa_table(fit, which = "measures", diagnostics = dx,
                       caption = "My table", note = "My note", digits = 12)
      expect_equal(unname(apa$table$Estimate), unname(dx$measures$Estimate),
                   tolerance = 1e-11, info = key)
      expect_identical(apa$table$CIEligible, dx$measures$CIEligible, info = key)
      fair <- fair_average_table(fit, diagnostics = dx, reference = "zero", fair_se = TRUE)
      expect_false(any(fair$stacked$FairCIEligible), info = key)
      expect_equal(fair$raw_by_facet,
                   fair_average_table(fit, reference = "zero", fair_se = TRUE)$raw_by_facet,
                   info = key)
      result <- mfrm_results(fit, diagnostics = dx,
                             include = c("fit", "diagnostics", "precision"))
      report <- mfrm_report(result)
      expect_identical(report$precision_evidence_summary$SupportsFormalInference,
                       formal, info = key)
    }
  })

  test_that("toy diagnostics follow the supplied fit on every call", {
    for (key in c("JML_RSM", "MML_RSM", "JML_RSM")) {
      dx <- make_toy_diagnostics(fits[[key]])
      expect_equal(dx$measures, diagnostics[[key]]$measures, info = key)
      expect_equal(dx$fit_readiness, diagnostics[[key]]$fit_readiness, info = key)
    }
    changed <- fits$JML_RSM
    changed$readiness$fit$InferenceReady <- FALSE
    changed$readiness$fit$FitReadiness <- "blocked"
    dx <- make_toy_diagnostics(changed)
    expect_false(dx$fit_readiness$InferenceReady)
    expect_identical(dx$fit_readiness$FitReadiness, "blocked")
  })

  routes <- list(
    summary_fit = function(f, d) summary(f, diagnostics = d),
    summary_facets = function(f, d) summary(f, profile = "facets", diagnostics = d),
    precision = function(f, d) precision_review_report(f, d),
    checklist = function(f, d) reporting_checklist(f, d),
    results = function(f, d) mfrm_results(f, diagnostics = d),
    apa = function(f, d) apa_table(f, which = "measures", diagnostics = d),
    apa_custom = function(f, d) apa_table(f, which = "measures", diagnostics = d,
                                         caption = "My table", note = "My note"),
    apa_facets = function(f, d) apa_table(f, which = "measures", diagnostics = d,
                                         branch = "facets"),
    apa_report = function(f, d) build_apa_outputs(f, d),
    wright = function(f, d) plot(f, diagnostics = d, draw = FALSE),
    wright_no_ci = function(f, d) plot(f, diagnostics = d, show_ci = FALSE, draw = FALSE),
    wright_unified = function(f, d) plot_wright_unified(f, diagnostics = d, draw = FALSE),
    fit_pathway = function(f, d) plot(f, type = "fit_pathway", diagnostics = d, draw = FALSE),
    visual_summary = function(f, d) build_visual_summaries(f, d),
    fair_table = function(f, d) fair_average_table(f, diagnostics = d),
    fair_plot = function(f, d) plot_fair_average(f, diagnostics = d, show_ci = TRUE, draw = FALSE),
    fair_plot_no_ci = function(f, d) plot_fair_average(f, diagnostics = d, show_ci = FALSE, draw = FALSE),
    manifest = function(f, d) build_mfrm_manifest(f, diagnostics = d),
    replay = function(f, d) build_mfrm_replay_script(f, diagnostics = d)
  )
  test_that("all connected reporting routes refuse mixed estimator or model diagnostics", {
    for (key in setdiff(names(fits), "MML_RSM")) {
      for (route in names(routes)) {
        expect_error(routes[[route]](fits[[key]], diagnostics$MML_RSM),
                     "fit/diagnostics mismatch", info = paste(key, route))
      }
    }
    # Refuse the opposite direction too; a conservative label cannot excuse
    # reporting estimates and standard errors from a different analysis.
    expect_error(summary(fits$MML_RSM, diagnostics = diagnostics$JML_RSM),
                 "fit/diagnostics mismatch")
  })

  test_that("missing locations or stale readiness cannot bypass source matching", {
    incomplete <- diagnostics$MML_RSM
    incomplete$measures <- NULL
    expect_error(mfrm_results(fits$JML_RSM, diagnostics = incomplete),
                 "fit/diagnostics mismatch")
    stale <- diagnostics$JML_RSM
    stale$fit_readiness$InferenceReady <- !stale$fit_readiness$InferenceReady
    expect_error(precision_review_report(fits$JML_RSM, stale),
                 "fit/diagnostics mismatch")
    legacy <- diagnostics$JML_RSM
    legacy$fit_readiness <- NULL
    for (route in names(routes)) {
      expect_error(routes[[route]](fits$JML_RSM, legacy),
                   "current readiness contract", info = route)
    }
  })
})
