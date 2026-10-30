test_that("two-family adaptive fitting preserves its objective and output scope", {
  data <- expand.grid(Candidate = paste0("p", 1:12),
    A = c("same", "a2", "a3"), B = c("same", "b2"))
  data <- data[!(data$A == "a3" & data$B == "b2"), ]
  data$Rating <- rep(c(0, 2, 1, 1, 2, 0, 2, 1, 0, 1), length.out = nrow(data))
  names(data)[2:3] <- c("First owner", "Second owner")
  args <- list(data = data, person = "Candidate", facets = names(data)[3:2],
    score = "Rating", model = "GPCM", method = "MML",
    slope_facet = names(data)[2:3], step_facet = names(data)[3],
    noncenter_facet = names(data)[3], gpcm_mml_identification = "fixed_standard_normal",
    mml_engine = "direct", mml_integration = "adaptive", quad_points = 7L, maxit = 1L)
  fit <- suppressWarnings(do.call(fit_mfrm, args))
  expect_s3_class(fit, "mfrm_fit")
  expect_identical(fit$summary$MMLEngineUsed, "direct")
  expect_identical(fit$summary$MMLIntegration, "adaptive")
  expect_identical(fit$summary$ConvergenceBasis, "optimizer_gradient")
  expect_false(fit$summary$Converged)
  expect_false(fit$summary$InferenceReady)
  expect_false(fit$summary$ICEligible)
  expect_null(fit$opt$em_trace)
  expect_null(fit$config$replay_inputs$em_score_tol)
  expect_identical(fit$config$replay_inputs$reltol, 1e-9)
  expect_identical(fit$config$replay_inputs$gpcm_mml_start, "neutral_em")
  expect_identical(fit$opt$mml_initialization$table$Start, c("neutral", "em"))
  expect_equal(fit$opt$value, min(fit$opt$mml_initialization$table$FinalNLL))
  expect_true(all(summary(fit)$settings_overview$GpcmMmlInitialization == "neutral_em"))
  setup <- do.call(mfrm_gmfrm_problem, fit$gmfrm$specification)$common
  review <- mfrmr_adaptive_quadrature_review(setup$idx, fit$config,
    expand_params(fit$opt$par, setup$sizes, fit$config),
    gauss_hermite_normal(7L), fit$prep$levels$Person, 7L)
  expect_equal(fit$summary$LogLik, sum(review$AdaptiveLogMarginal), tolerance = 1e-12)
  ci <- confint(fit, parm = "slopes")
  expect_true(all(is.na(ci)))
  expect_match(attr(ci, "diagnostics")$InferenceReview[1], "Direct adaptive MML must pass")
  expect_output(print(fit), "adaptive Gauss-Hermite")
  expect_identical(summary(fit)$slope_overview$SlopeOwner, names(data)[2:3])
  grid <- data.frame(Theta = c(-1, 0, 1))
  grid[[names(data)[2]]] <- "a2"; grid[[names(data)[3]]] <- "b2"
  curves <- mfrm_curve_intervals(fit, grid)
  expect_s3_class(as_ggplot(curves), "ggplot")
  results <- mfrm_results(fit, intervals = list(curves = curves, slopes = ci))
  report <- mfrm_report(results)
  expect_identical(results$tables$gpcm_fitted_locations$FittedLocation, unname(fit$facets$others$Estimate))
  expect_identical(results$tables$gpcm_fitted_steps$FittedStep, unname(fit$steps$Estimate))
  expect_identical(report$tables$fitted_locations, results$tables$gpcm_fitted_locations)
  expect_identical(report$tables$fitted_steps, results$tables$gpcm_fitted_steps)
  expect_identical(results$tables$gpcm_initialization, fit$opt$mml_initialization$table)
  expect_identical(report$tables$gpcm_initialization, results$tables$gpcm_initialization)
  expect_match(report$tables$model_settings$Value[2], "direct adaptive")
  expect_match(report$markdown, "gpcm initialization")
  expect_match(report$markdown, "Direct adaptive MML did not pass")
  expect_false(grepl("EM stopped|EM stopping", report$markdown))
  expect_match(report$markdown, "0 experimental slope intervals are available")
  unresolved <- results
  unresolved$fit$summary$Converged <- TRUE
  unresolved$fit$opt$optimizer_diagnostics$ConvergenceSeverity <- "review"
  expect_match(mfrm_report(unresolved)$markdown, "did not pass its numerical")
  # A component profile must reject the same unconverged source.
  expect_error(confint(fit, parm = "slopes", method = "profile", slope = c(`First owner` = "a2")),
    "Direct adaptive MML must pass")
  expect_error(mfrm_response_diagnostics(fit), "Resolve the two-family convergence")
  expect_error(predict_mfrm_units(fit, data), "Resolve the two-family convergence")
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
  saveRDS(fit, file); restored <- readRDS(file)
  expect_identical(restored, fit)
  replay <- restored$config$replay_inputs
  replay$package_version <- NULL; replay$data <- data
  expect_identical(suppressWarnings(do.call(fit_mfrm, replay))$summary, fit$summary)
  comparison <- suppressWarnings(mml_quadrature_sensitivity(fit, data,
    quad_points = c(7, 9), adaptive_quad_points = c(7, 9)))
  expect_true(all(vapply(comparison$fits, function(x)
    identical(x$config$estimation_control$mml_integration, "adaptive"), TRUE)))
  expect_true(all(vapply(comparison$intervals, function(x) all(is.na(x)), TRUE)))
  expect_true(all(comparison$quadrature_review$Status == "computed"))
  expect_true(all(vapply(comparison$fits, function(x)
    identical(x$config$replay_inputs$gpcm_mml_start, "neutral_em"), TRUE)))
  legacy <- fit
  legacy$config$replay_inputs$gpcm_mml_start <- NULL
  legacy$config$estimation_control$gpcm_mml_start <- NULL
  legacy$opt$mml_initialization <- NULL
  expect_identical(mfrmr_gqs_refit_arguments(legacy, data, 9L)$gpcm_mml_start, "neutral")
  expect_null(mfrm_results(legacy, include = "fit", compute = "never")$tables$gpcm_initialization)
  changed <- fit; changed$config$replay_inputs$gpcm_mml_start <- "neutral"
  expect_error(mfrmr_gqs_refit_arguments(changed, data, 9L), "initialization and replay")
  neutral <- args; neutral$gpcm_mml_start <- "neutral"
  old <- suppressWarnings(do.call(fit_mfrm, neutral))
  expect_identical(old$opt$par, fit$opt$mml_initialization$attempts$neutral$value$par)
  expect_identical(old$opt$mml_initialization$table$Start, "neutral")
  old_review <- mfrmr_adaptive_quadrature_review(setup$idx, old$config,
    expand_params(old$opt$par, setup$sizes, old$config),
    gauss_hermite_normal(7L), old$prep$levels$Person, 7L)
  expect_gt(abs(sum(old_review$FixedLogMarginal) - old$summary$LogLik), 1e-7)
  invalid <- args; invalid$gpcm_mml_start <- "truth"
  expect_error(do.call(fit_mfrm, invalid), "arg.*should be one of")
  invalid$mml_integration <- "fixed"; invalid$mml_engine <- "em"
  expect_error(do.call(fit_mfrm, invalid), "only to two-family adaptive")
  args$em_score_tol <- 1e-6
  expect_error(do.call(fit_mfrm, args), "not supported by the two-family route")
  args$em_score_tol <- NULL; args$reltol <- 0
  expect_error(do.call(fit_mfrm, args), "finite positive")
})

test_that("adaptive two-family intervals retain the full covariance and output identity", {
  saved <- readRDS(test_path("fixtures", "gmfrm-joint-information.rds"))
  d <- saved$data
  levels(d$Task) <- c("same", "t2", "t3")
  levels(d$Rater) <- c("same", "r2", "r3")
  names(d) <- c("Candidate", "First owner", "Second owner", "Rating")
  fit <- fit_mfrm(d, person = "Candidate", facets = names(d)[3:2], score = "Rating",
    model = "GPCM", method = "MML", slope_facet = names(d)[2:3],
    step_facet = names(d)[3], noncenter_facet = names(d)[3],
    gpcm_mml_identification = "fixed_standard_normal", mml_engine = "direct",
    mml_integration = "adaptive", quad_points = 31L, maxit = 500L, reltol = 1e-10)
  expect_true(fit$summary$Converged)
  before <- fit
  expect_warning(ci <- confint(fit, parm = "slopes"), "Experimental two-family")
  expect_identical(fit, before)
  tab <- attr(ci, "diagnostics")
  expect_true(all(tab$CIEligible))
  expect_identical(tab$SlopeOwner, rep(names(d)[2:3], each = 3))
  expect_identical(tab$SlopeOwner[tab$SlopeLevel == "same"], names(d)[2:3])
  expect_true(all(is.na(tab$PValue)))
  expect_identical(attr(ci, "settings")$integration, "adaptive")
  expect_identical(attr(ci, "settings")$engine, "direct")
  expect_identical(attr(ci, "numerical_checks")$Integration, "adaptive")
  expect_true(all(attr(ci, "checks")$Passed))
  setup <- do.call(mfrm_gmfrm_problem, fit$gmfrm$specification)$common
  parameters <- expand_params(fit$opt$par, setup$sizes, fit$config)
  basis <- mfrmr_adaptive_quadrature_basis(setup$idx, fit$config, parameters,
    gauss_hermite_normal(61L), compute_base_eta(setup$idx, parameters, fit$config))
  reference <- gmfrm_louis_reference(fit$gmfrm$specification, fit$opt$par,
    basis$nodes, basis$log_weights)
  jac <- matrix(0, 6, 13)
  jac[1:3, 9:10] <- rbind(c(1, 0), c(0, 1), c(-1, -1))
  jac[4:6, 11:13] <- diag(3)
  covariance <- jac %*% solve(reference$information) %*% t(jac)
  expect_equal(unname(attr(ci, "covariance")), unname(covariance), tolerance = 1e-5)
  expect_equal(tab$LogSE, sqrt(diag(covariance)), tolerance = 1e-5)
  expect_gt(max(abs(covariance[1:3, 4:6])), 1e-5)
  expect_equal(colSums(covariance[1:3, ]), rep(0, 6), tolerance = 1e-10)
  expect_equal(mfrm_mml_person_scores_numeric(fit), reference$scores,
    tolerance = 1e-6, ignore_attr = TRUE)
  selected <- which(as.character(d$Candidate) %in% unique(as.character(d$Candidate))[c(1, 20, 160)])
  response <- mfrm_response_diagnostics(fit, rows = selected, group_by = names(d)[3], quad_points = 31L)
  expect_true(all(response$rows$Status == "available_conditional"))
  expect_match(response$settings$integration, "Adaptive quadrature under N\\(0,1\\)")
  expect_false(response$settings$calibration_uncertainty)
  expect_identical(names(response$source_data), names(d))
  scoring_data <- d[d$Candidate == d$Candidate[1], ]
  scoring_data$Candidate <- "NEW"
  scored <- predict_mfrm_units(fit, scoring_data, n_draws = 3, seed = 21)
  expect_identical(scored$settings$source_scoring_status, "conditional")
  expect_identical(scored$settings$scoring_algorithm, "adaptive_quadrature_eap_v2")
  expect_true(all(scored$estimates$ScoreIntegrationReady))
  expect_identical(scored$settings$two_family_calibration$slope_facets, names(d)[2:3])
  literal_data <- scoring_data; literal_data$Score <- literal_data$Rating
  expected_score <- gmfrm_scoring_reference(fit, literal_data)
  expect_equal(unname(unlist(scored$estimates[c("Estimate", "SD", "Lower", "Upper")])),
    unname(expected_score), tolerance = 1e-7)
  expect_s3_class(summary(scored), "summary.mfrm_unit_prediction")
  calibration <- freeze_mfrm_calibration(validate_mfrm_calibration(extract_mfrm_calibration(fit)))
  portable <- score_mfrm_calibration(calibration, scoring_data)
  expect_identical(calibration$header$schema_version, 6L)
  expect_identical(portable$settings$scoring_algorithm, "adaptive_quadrature_eap_v2")
  expect_equal(unname(as.matrix(portable$estimates[c("Estimate", "SD", "Lower", "Upper")])),
    unname(as.matrix(scored$estimates[c("Estimate", "SD", "Lower", "Upper")])), tolerance = 1e-10)
  expect_s3_class(summary(portable), "summary.mfrm_calibration_score")
  rows <- selected[d$Candidate[selected] == d$Candidate[selected[1]]]
  literal <- gmfrm_response_reference(fit, rows)
  at <- match(rows, selected)
  expect_equal(unname(response$probabilities[at, ]), literal$probabilities, tolerance = 1e-8)
  expect_equal(response$rows$ExpectedScore[at], literal$mean, tolerance = 1e-8)
  expect_equal(response$rows$PredictiveVariance[at], literal$variance, tolerance = 1e-8)
  expect_gt(max(literal$variance - literal$conditional_variance), .001)
  expect_gt(max(abs(literal$plugin - literal$probabilities)), .001)
  subset <- rev(selected[c(2, 8, 14)])
  fewer <- mfrm_response_diagnostics(fit, rows = subset, quad_points = 31L)
  expect_identical(fewer$rows$InputRow, subset)
  expect_equal(unname(fewer$probabilities), unname(response$probabilities[match(subset, selected), ]))
  # The original source is required, not a fabricated shortened conditioning record.
  incomplete <- fit
  incomplete$gmfrm$specification$data <- d[-selected[1], ]
  expect_error(mfrm_response_diagnostics(incomplete, rows = selected[2]), "must match")
  unready <- fit; unready$opt$optimizer_diagnostics$ConvergenceSeverity <- "review"
  expect_error(mfrm_response_diagnostics(unready, rows = 1), "convergence")
  changed <- fit; changed$opt$value <- changed$opt$value + 1
  expect_error(mfrm_response_diagnostics(changed, rows = 1), "adaptive likelihood")
  changed <- fit; changed$opt$par[1] <- changed$opt$par[1] + .1
  expect_error(mfrm_response_diagnostics(changed, rows = 1), "adaptive likelihood")
  changed <- fit; changed$opt$optimizer_diagnostics$GradientReviewTolerance <- 1
  expect_error(mfrm_response_diagnostics(changed, rows = 1), "adaptive likelihood")
  changed <- fit; changed$config$estimation_control$mml_integration <- "fixed"
  expect_error(mfrm_response_diagnostics(changed, rows = 1), "original fixed-N")
  low <- mfrm_response_diagnostics(fit, rows = selected, quad_points = 7L, group_by = names(d)[3])
  expect_true(any(low$rows$Status == "unavailable"))
  expect_equal(sum(low$measures$Observed), length(selected))
  unresolved <- low$measures$Available < low$measures$Observed
  expect_true(all(is.na(low$measures$Infit[unresolved])))
  expect_true(all(is.na(low$measures$Outfit[unresolved])))
  results <- mfrm_results(fit, intervals = list(slopes = ci), response_diagnostics = response)
  expect_match(mfrm_report(results)$markdown, "6 experimental slope intervals are available")
  expect_identical(results$tables$gpcm_slopes_settings$Integration, "adaptive")
  expect_s3_class(as_ggplot(ci), "ggplot")
  folder <- tempfile("gmfrm-adaptive-ci-"); withr::defer(unlink(folder, recursive = TRUE))
  exported <- export_mfrm_results(results, output_dir = folder, preset = "starter",
    acknowledge_sensitive = TRUE, plot_width = 1200, plot_height = 900)
  expect_equal(nrow(exported$plot_errors), 0L)
  restored <- readRDS(list.files(folder, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)[1])
  expect_identical(restored$gpcm_inference$slopes, ci)
  expect_identical(restored$response_diagnostics, response)
  expect_identical(restored$fit$opt$mml_initialization, fit$opt$mml_initialization)
  expect_identical(restored$tables$gpcm_initialization, results$tables$gpcm_initialization)
  expect_identical(restored$tables$response_residuals, response$rows)
  expect_identical(summary(restored$response_diagnostics), summary(response))
  expect_identical(mfrm_report(restored)$markdown, mfrm_report(results)$markdown)
  bad <- fit; bad$opt$value <- bad$opt$value + 1
  unavailable <- confint(bad, parm = "slopes")
  expect_true(all(is.na(unavailable)))
  expect_match(attr(unavailable, "diagnostics")$InferenceReview[1], "likelihood must match")
  expect_error(confint(bad, parm = "slopes", method = "profile", slope = c(`First owner` = "same")),
    "likelihood must match")
  expect_false(fit$summary$InferenceReady)
  expect_false(fit$summary$ICEligible)
  args <- fit$config$replay_inputs
  args$package_version <- NULL; args$data <- d; args$quad_points <- 3L
  coarse <- do.call(fit_mfrm, args)
  expect_true(coarse$summary$Converged)
  missing <- confint(coarse, parm = "slopes")
  expect_true(all(is.na(missing)))
  expect_identical(with(attr(missing, "checks"), Check[!Passed]), "Quadrature sensitivity")
  # Saved output must remain usable without reoptimizing or integrating.
  local_mocked_bindings(mfrmr_make_adaptive_mml_evaluator = function(...) stop("unexpected gradient"),
    mfrmr_adaptive_quadrature_basis = function(...) stop("unexpected integration"),
    mfrm_gpcm_response_evaluator = function(...) stop("unexpected probability"), .package = "mfrmr")
  attached <- mfrm_results(fit, intervals = list(slopes = ci), response_diagnostics = response)
  expect_identical(mfrm_report(attached)$markdown, mfrm_report(results)$markdown)
  expect_match(mfrm_report(attached)$markdown, "Adaptive quadrature")
  for (style in c("paired", "scatter")) {
    plotted <- plot(attached, type = "response_diagnostics", style = style, draw = FALSE)
    expect_identical(plot_data(plotted)$settings, response$settings)
    expect_match(plot_data(plotted)$caption, "no reference cutoffs")
    expect_s3_class(as_ggplot(plotted), "ggplot")
  }
})
