owner_summary_fixture <- function() {
  d <- expand.grid(Person = paste0("p", 1:4),
    Task = c("same:1", "same,2", "third"), Rater = c("same:1", "same,2"))
  d$Score <- rep(0:2, length.out = nrow(d))
  owners <- c("Rubric: facet", "評価 者")
  names(d)[2:3] <- owners
  p <- mfrm_gmfrm_problem(d, 2L, gauss_hermite_normal(9L), slope_facets = owners)
  par <- p$start
  par[build_param_slices(p$common$sizes)$log_slopes] <- log(c(2, 3, 2, 8))
  table <- build_slope_table(p$common$config, p$common$prep,
    expand_params(par, p$common$sizes, p$common$config))
  readiness <- mfrmr_readiness_gpcm_slope_parameters(p$common$config, table)
  list(problem = p, owners = owners, par = par,
    slopes = apply_mfrm_slope_readiness(table, list(parameters = readiness)))
}

test_that("slope summaries preserve owner-specific counts, scales and unavailable estimates", {
  x <- owner_summary_fixture(); config <- x$problem$common$config
  s <- mfrm_fit_slope_overview(config, x$slopes)
  expect_identical(s$SlopeOwner, x$owners)
  expect_identical(s$Slopes, c(3L, 2L))
  # Declared log slopes: (log(2), log(3), -log(6)) and (log(2), log(8)).
  expect_equal(s$OptimizerGeometricMean, c(1, 4), tolerance = 1e-12)
  expect_equal(s$OptimizerMin, c(1/6, 2), tolerance = 1e-12)
  expect_equal(s$OptimizerMax, c(3, 8), tolerance = 1e-12)
  expect_identical(s$ScaleReference, c("geometric_mean_one", "fixed_standard_normal"))
  expect_true(all(is.na(s$GeometricMean)))
  expect_true(all(s$ValueBasis == "optimizer_trace"))
  expect_true(all(!s$PrimaryReady & s$CIEligible == 0L))
  expect_identical(mfrm_fit_slope_overview(config, x$slopes[c(5,2,4,1,3), ]), s)

  # A controlled readiness-table change in one owner must not change the other.
  first <- x$slopes$SlopeOwner == x$owners[1]
  x$slopes$PrimaryEstimate[first] <- x$slopes$Estimate[first]
  x$slopes$ParameterStatus[first] <- "finite"
  changed <- mfrm_fit_slope_overview(config, x$slopes)
  expect_identical(changed$PrimaryReady, c(TRUE, FALSE))
  expect_identical(changed[2, ], s[2, ])
  expect_error(mfrm_fit_slope_overview(config, x$slopes[-1, ]), "every fitted owner")
  expect_error(mfrm_fit_slope_overview(config, rbind(x$slopes, x$slopes[1, ])), "unique within")
  expect_error(mfrm_fit_slope_overview(config, x$slopes[, -1]), "SlopeOwner")
  config$slope_facet <- rev(config$slope_facet)
  expect_error(mfrm_fit_slope_overview(config, x$slopes), "owner order")
})

test_that("shared summary and printing retain both references without advertising inference", {
  x <- owner_summary_fixture(); p <- x$problem
  # Exercise the shared consumer contract, not a new public model constructor.
  input <- list(config = p$common$config, prep = p$common$prep, slopes = x$slopes,
    summary = data.frame(Model = "GPCM", Method = "MML", N = 24L, Persons = 4L,
      Facets = 2L, Categories = 3L, MMLEngineUsed = "em", MMLEngineRequested = "em",
      LogLik = p$marginal(x$par)$logLik, Converged = FALSE,
      ConvergenceSeverity = "review", ConvergenceStatus = "iteration_limit",
      TerminalGradientSupNorm = 1, ConvergenceDetail = "Iteration limit reached"))
  contract <- mfrm_fit_scale_contract(input)
  expect_identical(contract$SlopeBasis, "owner_specific_product_discrimination")
  expect_identical(contract$GpcmModelFamily, "two_owner_product_slope_gpcm")
  expect_identical(contract$GpcmSlopeComposition, "two_owner_product_first_gm1")
  expect_identical(contract$PopulationSD, 1)
  expect_identical(contract$FixedLatentSDSlopeField, "OptimizerEstimate")
  expect_no_warning(s <- mfrm_fit_summary_core(input))
  expect_identical(s$slope_overview, mfrm_fit_slope_overview(input$config, input$slopes))
  expect_identical(s$settings_overview$SlopeFacet, x$owners)
  expect_true(all(s$settings_overview$StepFacet == x$owners[2]))
  printed <- gsub("[[:space:]]+", " ", paste(capture.output(print(s)), collapse = " "))
  expect_match(printed, "Within-facet geometric mean 1", fixed = TRUE)
  expect_match(printed, "Fixed standard-normal ability scale", fixed = TRUE)
  expect_match(printed, "second family free", fixed = TRUE)
  expect_match(printed, "Formal inference: No", fixed = TRUE)
  expect_match(printed, "experimental component-slope intervals", fixed = TRUE)
  expect_false(grepl("plot(fit", printed, fixed = TRUE))
  expect_false(grepl("free_given_first_family_and_population", printed, fixed = TRUE))
  printed_fit <- gsub("[[:space:]]+", " ",
    paste(capture.output(print(structure(input, class = "mfrm_fit"))), collapse = " "))
  expect_match(printed_fit, "product of facet slopes", fixed = TRUE)
  expect_false(grepl("relative slopes with geometric mean 1", printed_fit, fixed = TRUE))
  path <- tempfile(fileext = ".rds"); withr::defer(unlink(path))
  saveRDS(s, path)
  expect_identical(readRDS(path), s)
  input$population <- list(active = TRUE, sigma2 = 2)
  expect_error(mfrm_fit_scale_contract(input), "fixed population scale")
})

test_that("internal EM saves shared summaries without promoting convergence or inference", {
  x <- owner_summary_fixture()
  fit <- mfrm_gmfrm_em(x$problem, start = x$par, maxit = 1L, score_tol = 1e-12)
  expect_false(inherits(fit, "mfrm_fit"))
  expect_identical(fit$slope_overview$SlopeOwner, x$owners)
  expect_equal(fit$slope_overview$OptimizerGeometricMean[1], 1, tolerance = 1e-12)
  expect_true(all(!fit$slope_overview$PrimaryReady & fit$slope_overview$CIEligible == 0L))
  expect_false(fit$converged)
  expect_identical(fit$reason, "iteration_limit")
  expect_identical(fit$scale_contract, mfrm_fit_scale_contract(list(config = x$problem$common$config)))
})

test_that("saved public single-owner summaries retain their existing values and scale", {
  fit <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))$fit
  s <- summary(fit, compute = "never")
  expect_equal(nrow(s$slope_overview), 1L)
  expect_identical(s$slope_overview$SlopeOwner, fit$config$slope_facet)
  expect_equal(s$slope_overview$OptimizerMin, min(fit$slopes$OptimizerEstimate))
  expect_equal(s$slope_overview$OptimizerMax, max(fit$slopes$OptimizerEstimate))
  expect_equal(s$slope_overview$OptimizerGeometricMean, 1, tolerance = 1e-12)
  expect_false("ScaleReference" %in% names(s$slope_overview))
  expect_identical(mfrm_fit_scale_contract(fit)$SlopeBasis, "geometric_mean_one_relative_discrimination")
  printed <- gsub("[[:space:]]+", " ", paste(capture.output(print(fit)), collapse = " "))
  expect_match(printed, "relative slopes with geometric mean 1", fixed = TRUE)
})

test_that("an actual EM result retains likelihood, posterior, stopping rule and owner identity", {
  x <- owner_summary_fixture(); p <- x$problem
  em <- mfrm_gmfrm_em(p, start = x$par, maxit = 1L, score_tol = 1e-12)
  fit <- mfrm_gmfrm_fit_result(p, em)
  expect_s3_class(fit, "mfrm_fit")
  expect_identical(fit$opt$par, em$par)
  expect_equal(fit$summary$LogLik, em$logLik, tolerance = 1e-12)
  expect_equal(fit$summary$Npar, length(em$par))
  expect_equal(fit$summary$Persons, p$n_person)
  expect_equal(fit$summary$N, nrow(p$specification$data))
  expect_identical(fit$summary$MMLEngineUsed, "em")
  expect_identical(fit$summary$IterationsBasis, "em_iterations")
  expect_identical(fit$summary$Iterations, 1L)
  expect_identical(fit$summary$ConvergenceBasis, "marginal_score_per_person")
  expect_equal(fit$summary$TerminalGradientSupNorm, em$max_score)
  expect_equal(fit$summary$GradientReviewTolerance, 1e-12)
  expect_identical(fit$opt$em_trace, em$trace)
  expect_false(fit$summary$Converged)
  expect_false(fit$summary$InferenceReady)
  expect_false(any(fit$summary$ICEligible, fit$summary$ICSelectable))
  expect_true(all(is.na(fit$summary[
    c("AIC", "BIC", "SABIC", "LegacyAIC", "LegacyBIC")
  ])))
  expect_identical(fit$summary$ICStatus, "unsupported_product_slopes")
  # Check posterior moments directly against the retained E-step, not another
  # fitted-object consumer, and retain calibration-conditional uncertainty.
  expected <- drop(em$posterior %*% p$nodes)
  sd <- sqrt(rowSums(em$posterior * sweep(matrix(p$nodes, p$n_person,
    length(p$nodes), byrow = TRUE), 1L, expected, "-")^2))
  expect_equal(fit$facets$person$Estimate, unname(expected), tolerance = 1e-12)
  expect_equal(fit$facets$person$PosteriorSD, unname(sd), tolerance = 1e-12)
  expect_true(all(is.na(fit$slopes$PrimaryEstimate)))

  expect_no_warning(s <- summary(fit, compute = "never"))
  expect_identical(s$slope_overview$SlopeOwner, x$owners)
  expect_identical(s$section_status$Status[s$section_status$Section == "fit_summary"], "ok")
  expect_false(any(s$required_visual$Available | s$required_visual$Required))
  expect_true(all(is.na(s$required_visual$Route)))
  printed <- gsub("[[:space:]]+", " ",
    paste(c(capture.output(print(fit)), capture.output(print(s))), collapse = " "))
  expect_match(printed, "experimental component-slope intervals", fixed = TRUE)
  expect_false(grepl("compare_mfrm\\(|plot\\(fit", printed))
  expect_error(summary(fit, profile = "facets", compute = "never"), "only the fit summary")
  expect_error(summary(fit, profile = "reporting", compute = "never"), "only the fit summary")
  newdata <- data.frame(Theta = c(-1, 0, 1))
  for (owner in x$owners) newdata[[owner]] <- p$levels[[owner]][1]
  evaluator <- mfrm_gpcm_response_evaluator(fit$prep, fit$config, newdata)
  expect_equal(evaluator$evaluate(fit$opt$par)$probabilities,
    p$response(em$par, newdata)$probabilities, tolerance = 1e-12)
  path <- tempfile(fileext = ".rds"); withr::defer(unlink(path))
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_identical(summary(restored, compute = "never"), s)
  expect_identical(restored$gmfrm$specification, p$specification)
  bad <- em; bad$converged <- TRUE
  expect_error(mfrm_gmfrm_fit_result(p, bad), "convergence flag")
  bad <- em; bad$specification$slope_facets <- rev(x$owners)
  expect_error(mfrm_gmfrm_fit_result(p, bad), "specification must match")
  bad <- em; bad$logLik <- em$logLik + 1
  expect_error(mfrm_gmfrm_fit_result(p, bad), "saved EM likelihood")
  shifted <- p; shifted$specification$quadrature$nodes <- p$nodes + .1
  bad <- em; bad$specification <- shifted$specification
  expect_error(mfrm_gmfrm_fit_result(shifted, bad), "standard-normal Gauss-Hermite")
})

test_that("EM score convergence does not imply inference or information-criterion eligibility", {
  x <- owner_summary_fixture(); p <- x$problem
  # Deliberately loose tolerance tests the separation of numerical stopping
  # from statistical readiness. This is not a well-estimated model example.
  tol <- 2 * max(abs(p$marginal(p$start)$gradient))
  em <- mfrm_gmfrm_em(p, score_tol = tol)
  expect_true(em$converged)
  fit <- mfrm_gmfrm_fit_result(p, em)
  expect_true(fit$summary$Converged)
  expect_identical(fit$summary$ConvergenceReason, "marginal_score_tolerance_met")
  expect_identical(fit$summary$Iterations, 0L)
  expect_identical(fit$readiness$fit$NumericalState, "ready")
  expect_false(fit$readiness$fit$InferenceReady)
  expect_false(fit$summary$ICEligible)
  expect_true(all(is.na(fit$slopes$PrimaryEstimate)))
  expect_error(build_optimizer_diagnostics(list(convergence=0L), gradient=0,
    convergence_basis="marginal_score_per_person"), "stopping tolerance")
})
