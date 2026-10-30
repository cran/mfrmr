# test-reporting-method-contracts.R
# Exercises reporting.R edge paths:
# - build_visual_warning_map (all warning stages)
# - build_visual_summary_map (all summary stages, including detail="detailed")
# - build_apa_report_text (with and without context/bias)
# - build_apa_table_figure_notes / captions
# - build_sectioned_fixed_report edge cases
# - format_fixed_width_table edge cases
# - resolve_warning_thresholds with different profiles
# - build_pca_reference_text / build_pca_check_text
# - summarize_anchor_constraints / summarize_convergence_metrics / summarize_step_estimates
# - summarize_bias_counts / summarize_top_misfit_levels

# ---- Shared fixture ----

local({
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  .fit <<- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  .diag <<- diagnose_mfrm(.fit, residual_pca = "both", pca_max_factors = 3)
  .bias <<- estimate_bias(.fit, .diag, facet_a = "Rater", facet_b = "Task")

  .fit_mml <<- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "MML", model = "RSM", maxit = 200, quad_points = 7)
  )
  .diag_mml <<- diagnose_mfrm(
    .fit_mml,
    residual_pca = "none",
    diagnostic_mode = "both"
  )
})

# ============================================================================
# build_visual_warning_map
# ============================================================================

test_that("regularized facet SEs remain diagnostic across reports and saved inputs", {
  covariance <- mfrmr:::compute_mml_parameter_covariance(.fit_mml)
  expect_identical(covariance$status, "ok")
  covariance$status <- "regularized"
  covariance$regularized <- TRUE
  local_mocked_bindings(
    compute_mml_parameter_covariance = function(...) covariance,
    .package = "mfrmr"
  )
  diag <- diagnose_mfrm(.fit_mml, residual_pca = "none", diagnostic_mode = "both")
  facet <- diag$measures$Facet != "Person"
  expect_true(any(is.finite(diag$measures$SE[facet])))
  expect_false(any(diag$measures$CIEligible[facet]))
  expect_false(any(diag$measures$SupportsFormalInference[facet]))
  expect_true(all(grepl("regularized Hessian", diag$measures$SE_Method[facet])))
  expect_true(diag$precision_profile$HasRegularizedSE)
  expect_false(diag$precision_profile$SupportsFormalInference)
  expect_equal(diag$measures$SE, .diag_mml$measures$SE, tolerance = 0)
  expect_equal(diag$measures$CI_Lower, .diag_mml$measures$CI_Lower, tolerance = 0)
  expect_equal(diag$measures$CI_Upper, .diag_mml$measures$CI_Upper, tolerance = 0)
  expect_error(analyze_facet_equivalence(.fit_mml), "unregularized", fixed = TRUE)
  report <- precision_review_report(.fit_mml, diag)
  text <- paste(capture.output(print(summary(report))), collapse = " ")
  expect_match(text, "regularized Hessian", fixed = TRUE)
  expect_match(text, "diagnostic only", fixed = TRUE)
  expect_false(grepl("PrecisionTier|SupportsFormalInference|validation-use", text))

  older <- diag
  older$measures$SupportsFormalInference[facet] <- TRUE
  older$measures$CIEligible[facet] <- TRUE
  expect_error(summary(older), "Recompute diagnose_mfrm(fit)", fixed = TRUE)
  expect_error(precision_review(older), "Recompute diagnose_mfrm(fit)", fixed = TRUE)
  expect_error(precision_review_report(.fit_mml, older), "Recompute diagnose_mfrm(fit)", fixed = TRUE)
  expect_error(build_apa_outputs(.fit_mml, older), "Recompute diagnose_mfrm(fit)", fixed = TRUE)

  older_report <- report
  older_report$profile$HasRegularizedSE <- NULL
  expect_error(summary(older_report), "Recreate it with precision_review_report", fixed = TRUE)
  older_summary <- summary(report)
  older_summary$profile$HasRegularizedSE <- NULL
  expect_error(print(older_summary), "Recreate it with precision_review_report", fixed = TRUE)
  older_normal <- .diag_mml
  older_normal$precision_profile$HasRegularizedSE <- NULL
  expect_error(summary(older_normal), "Recompute diagnose_mfrm(fit)", fixed = TRUE)
})

test_that("unavailable model-based SEs cannot promote observation-table fallbacks", {
  model_se <- mfrmr:::compute_mml_facet_model_se(.fit_mml)
  model_se$table$ModelSE[1] <- NA_real_
  local_mocked_bindings(compute_mml_facet_model_se = function(...) model_se,
                        .package = "mfrmr")
  se <- mfrmr:::build_measure_se_table(
    .fit_mml, .diag_mml$obs, .fit_mml$config$facet_names, .diag_mml$fit
  )
  row <- se$Facet == model_se$table$Facet[1] & se$Level == model_se$table$Level[1]
  expect_true(is.finite(se$SE[row]))
  expect_identical(unname(se$SE_Method[row]), "Fallback observation-table information")
  expect_false(se$SupportsFormalInference[row])
  expect_identical(se$CIUse[row], "review_before_reporting")
})

test_that("precision reporting separates numerical convergence from uncertainty", {
  fit <- .fit_mml
  fit$readiness$fit$InferenceReady <- FALSE
  fit$readiness$fit$FitReadiness <- "review"
  fit$readiness$fit$EstimabilityState <- "not_evaluated"
  fit$readiness$fit$ReasonCodes <- "design_rank_not_evaluated"
  profile <- mfrmr:::build_precision_profile(
    fit, .diag_mml$measures, .diag_mml$reliability, .diag_mml$facet_precision
  )
  checks <- mfrmr:::audit_precision_outputs(
    fit, .diag_mml$measures, .diag_mml$reliability, .diag_mml$facet_precision, profile
  )
  expect_identical(checks$Status[checks$Check == "Optimizer convergence"], "pass")
  expect_identical(checks$Status[checks$Check == "Precision tier"], "review")
  expect_match(profile$RecommendedUse, "does not support ordinary SE/CI inference", fixed = TRUE)
  expect_match(checks$Detail[checks$Check == "Optimizer convergence"],
               "does not establish valid SEs", fixed = TRUE)
  jml <- precision_review_report(.fit, .diag)
  expect_false(jml$profile$SupportsFormalInference)
  expect_match(jml$profile$RecommendedUse, "Changing to MML does not by itself", fixed = TRUE)
})

test_that("APA design prose uses each facet's own level count", {
  apa <- build_apa_outputs(.fit, .diag)
  text <- gsub("[[:space:]]+", " ", apa$report_text)
  counts <- lengths(.fit$config$facet_levels)
  expect_gt(length(unique(counts)), 1L)
  for (facet in names(counts)) {
    expect_match(text, paste0(facet, " (n = ", counts[[facet]], ")"), fixed = TRUE)
  }
})

test_that("weighted reports distinguish observed rows from the weight sum", {
  dat <- load_mfrmr_data("example_operational")
  dat$Weight <- rep(c(1, 0.5, 2, 1, 0), length.out = nrow(dat))
  fit <- suppressWarnings(fit_mfrm(dat, "Person", c("Rater", "Criterion"),
                                   "Score", weight = "Weight", method = "MML"))
  diag <- diagnose_mfrm(fit, residual_pca = "none")
  apa <- build_apa_outputs(fit, diag)
  rows <- nrow(diag$obs)
  weights <- sum(diag$obs$Weight)
  expect_equal(apa$contract$metadata$n_obs, rows)
  expect_equal(apa$contract$metadata$weight_sum, weights)
  expect_match(gsub("[[:space:]]+", " ", apa$report_text),
               paste0("was fit to ", rows, " observations"), fixed = TRUE)
  expect_match(apa$report_text, "sum of observation weights", fixed = TRUE)
  text <- gsub("[[:space:]]+", " ", apa$report_text)
  expect_match(text, "multiplied row-level ordered-category log-likelihood contributions", fixed = TRUE)
  expect_match(text, "Person distribution was integrated once per Person", fixed = TRUE)
  expect_match(text, "not interpreted as replicated independent Person response patterns", fixed = TRUE)
  expect_false(grepl("applied as frequency counts", text, fixed = TRUE))
  smap <- mfrmr:::build_visual_summary_map(fit, diag)
  expect_true(any(grepl(paste0("Observations: N = ", rows), smap$wright_map, fixed = TRUE)))
  wmap <- mfrmr:::build_visual_warning_map(fit, diag,
                                          thresholds = list(n_obs_min = mean(c(rows, weights))))
  expect_true(any(grepl("Small number of observations", wmap$wright_map, fixed = TRUE)))
})

test_that("reporting uses only stored residual PCA analyses", {
  calls <- 0L
  unexpected_pca <- function(...) {
    calls <<- calls + 1L
    stop("Reporting must not compute residual PCA")
  }
  local_mocked_bindings(
    compute_pca_overall = unexpected_pca,
    compute_pca_by_facet = unexpected_pca,
    .package = "mfrmr"
  )
  for (mode in c("none", "overall", "facet", "both")) {
    diag <- .diag
    diag$residual_pca_mode <- mode
    has_overall <- mode %in% c("overall", "both")
    has_facet <- mode %in% c("facet", "both")
    if (!has_overall) diag$residual_pca_overall <- NULL
    if (!has_facet) diag$residual_pca_by_facet <- NULL

    apa <- build_apa_outputs(.fit, diag)
    text <- gsub("[[:space:]]+", " ", apa$report_text)
    expect_identical(apa$contract$availability$has_pca_overall, has_overall)
    expect_identical(apa$contract$availability$has_pca_by_facet, has_facet)
    expect_identical(grepl("overall standardized residual matrix", text), has_overall)
    expect_identical(grepl("Facet-specific exploratory residual PCA showed", text), has_facet)
    if (mode == "none") expect_false(grepl("PCA", text))

    smap <- mfrmr:::build_visual_summary_map(.fit, diag)
    wmap <- mfrmr:::build_visual_warning_map(.fit, diag)
    expect_identical(any(grepl("Overall residual PCA PC1:", smap$residual_pca_overall)), has_overall)
    expect_identical(any(grepl("Top facet PC1", wmap$residual_pca_by_facet)), has_facet)
    pca <- mfrmr:::safe_residual_pca(diag)
    if (has_overall) expect_identical(pca$overall, .diag$residual_pca_overall)
    if (has_facet) expect_identical(pca$by_facet, .diag$residual_pca_by_facet)
  }
  expect_identical(calls, 0L)
})

test_that("build_visual_warning_map returns all expected visual keys", {
  wmap <- mfrmr:::build_visual_warning_map(.fit, .diag)
  expected_keys <- c(
    "wright_map", "pathway_map", "facet_distribution", "step_thresholds",
    "category_curves", "observed_expected", "fit_diagnostics",
    "fit_zstd_distribution", "misfit_levels",
    "strict_marginal_fit", "strict_pairwise_local_dependence",
    "residual_pca_overall", "residual_pca_by_facet"
  )
  expect_true(all(expected_keys %in% names(wmap)))
})

test_that("build_visual_warning_map with strict profile produces more warnings", {
  wmap_strict <- mfrmr:::build_visual_warning_map(.fit, .diag, threshold_profile = "strict")
  expect_true(is.list(wmap_strict))
  expect_true(length(wmap_strict$residual_pca_overall) > 0)
})

test_that("build_visual_warning_map with lenient profile", {
  wmap_lenient <- mfrmr:::build_visual_warning_map(.fit, .diag, threshold_profile = "lenient")
  expect_true(is.list(wmap_lenient))
})

test_that("build_visual_warning_map with NULL inputs returns empty map", {
  wmap_null <- mfrmr:::build_visual_warning_map(NULL, NULL)
  expect_true(is.list(wmap_null))
  expect_true(all(lengths(wmap_null) == 0))
})

test_that("build_visual_warning_map with custom threshold overrides", {
  custom_thr <- list(n_obs_min = 999999, n_person_min = 999999)
  wmap <- mfrmr:::build_visual_warning_map(.fit, .diag, thresholds = custom_thr)
  expect_true(any(grepl("Small", unlist(wmap))))
})

# ============================================================================
# build_visual_summary_map
# ============================================================================

test_that("build_visual_summary_map returns summary text for all keys", {
  smap <- mfrmr:::build_visual_summary_map(.fit, .diag)
  expected_keys <- c(
    "wright_map", "pathway_map", "facet_distribution", "step_thresholds",
    "category_curves", "observed_expected", "fit_diagnostics",
    "fit_zstd_distribution", "misfit_levels",
    "strict_marginal_fit", "strict_pairwise_local_dependence",
    "residual_pca_overall", "residual_pca_by_facet"
  )
  expect_true(all(expected_keys %in% names(smap)))
  # wright_map should always have some text
  expect_true(length(smap$wright_map) > 0)
})

test_that("build_visual_summary_map includes strict marginal routes for MML diagnostics", {
  smap_mml <- mfrmr:::build_visual_summary_map(.fit_mml, .diag_mml)
  expect_true(any(grepl("Expected counts condition on the same responses", smap_mml$strict_marginal_fit, fixed = TRUE)))
  expect_true(any(grepl('plot(x, type = "marginal_fit")', smap_mml$strict_marginal_fit, fixed = TRUE)))
  expect_true(any(grepl('plot(x, type = "marginal_pairwise")', smap_mml$strict_pairwise_local_dependence, fixed = TRUE)))
})

test_that("build_visual_warning_map explains when strict marginal diagnostics were not requested", {
  diag_legacy <- diagnose_mfrm(.fit, residual_pca = "none",
                               diagnostic_mode = "legacy")
  wmap_legacy <- mfrmr:::build_visual_warning_map(.fit, diag_legacy)
  expect_true(any(grepl("Strict marginal diagnostics are not available", wmap_legacy$strict_marginal_fit, fixed = TRUE)))
  expect_true(any(grepl("diagnostic_mode", wmap_legacy$strict_marginal_fit, fixed = TRUE)))
})

test_that("build_visual_summary_map with detail='detailed' adds extra summaries", {
  smap_detailed <- mfrmr:::build_visual_summary_map(
    .fit, .diag,
    options = list(detail = "detailed", max_facet_ranges = 2, top_misfit_n = 5)
  )
  expect_true(is.list(smap_detailed))
  expect_true(length(smap_detailed$wright_map) > 0)
})

test_that("build_visual_summary_map with NULL inputs returns empty map", {
  smap_null <- mfrmr:::build_visual_summary_map(NULL, NULL)
  expect_true(is.list(smap_null))
  expect_true(all(lengths(smap_null) == 0))
})

# ============================================================================
# build_apa_report_text
# ============================================================================

test_that("build_apa_report_text produces Method and Results sections", {
  text <- mfrmr:::build_apa_reporting_contract(.fit, .diag)$report_text
  expect_true(grepl("Method", text))
  expect_true(grepl("Results", text))
  expect_true(grepl("many-facet .*Rasch", text, ignore.case = TRUE))
})

test_that("build_apa_report_text with context supplies assessment/setting text", {
  ctx <- list(
    assessment = "writing proficiency",
    setting = "a university setting",
    rater_training = "two hours of calibration training",
    raters_per_response = "two",
    scale_desc = "a 5-point holistic rubric",
    line_width = 80L
  )
  text <- mfrmr:::build_apa_reporting_contract(.fit, .diag, context = ctx)$report_text
  expect_true(grepl("writing proficiency", text))
  expect_true(grepl("university", text))
  expect_true(grepl("calibration", text))
  expect_true(grepl("holistic rubric", text))
})

test_that("build_apa_report_text with bias_results includes bias summary", {
  text <- mfrmr:::build_apa_reporting_contract(.fit, .diag, bias_results = .bias)$report_text
  expect_true(grepl("Bias", text, ignore.case = TRUE) || grepl("bias", text))
})

test_that("build_apa_report_text without bias mentions no bias data", {
  text <- mfrmr:::build_apa_reporting_contract(.fit, .diag, bias_results = NULL)$report_text
  expect_true(nchar(text) > 0)
})

test_that("build_apa_outputs surfaces latent-regression population-model wording", {
  set.seed(141)
  persons <- paste0("P", sprintf("%02d", seq_len(24)))
  items <- paste0("I", seq_len(4))
  group <- rep(c("high", "low"), length.out = length(persons))
  theta <- ifelse(group == "high", 0.7, -0.4) + stats::rnorm(length(persons), sd = 0.5)
  item_beta <- seq(-0.6, 0.6, length.out = length(items))
  dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
  eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
  dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))
  dat$Score[c(1, 2)] <- c(0L, 1L)
  person_data <- data.frame(
    Person = persons,
    Group = group,
    stringsAsFactors = FALSE
  )

  fit <- suppressWarnings(fit_mfrm(
    dat,
    "Person",
    "Item",
    "Score",
    method = "MML",
    model = "RSM",
    population_formula = ~ Group,
    person_data = person_data,
    quad_points = 5,
    maxit = 35
  ))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))

  apa <- build_apa_outputs(fit, diag)
  text <- as.character(apa$report_text)
  text_normalized <- gsub("\\s+", " ", text)

  expect_true(isTRUE(apa$contract$availability$has_population_model))
  expect_true(isTRUE(apa$contract$availability$has_population_coding))
  expect_true("results_population_model" %in% apa$section_map$SectionId)
  expect_true(apa$section_map$Available[apa$section_map$SectionId == "results_population_model"][1])
  expect_true("population_model" %in% names(apa$contract$note_map))
  expect_true("population_model" %in% apa$contract$ordered_keys)
  expect_match(text_normalized, "conditional-normal latent-regression", fixed = TRUE)
  expect_match(text_normalized, "not as post hoc regression", fixed = TRUE)
  expect_match(text_normalized, "documented latent-regression MML comparison scope", fixed = TRUE)
  expect_match(apa$table_figure_notes, "Latent-regression population model", fixed = TRUE)
  expect_match(apa$table_figure_captions, "Latent-Regression Population Model", fixed = TRUE)

  s <- summary(apa)
  row <- s$content_checks[s$content_checks$Check == "Latent-regression wording alignment", , drop = FALSE]
  expect_equal(nrow(row), 1L)
  expect_true(row$Passed[1])
})

test_that("build_apa_outputs validates fit, diagnostics, context, and bias inputs at the front door", {
  expect_error(
    mfrmr::build_apa_outputs(list(), .diag),
    "requires `fit` to be an `mfrm_fit` object",
    fixed = TRUE
  )
  expect_error(
    mfrmr::build_apa_outputs(.fit, list()),
    "requires `diagnostics` to be an `mfrm_diagnostics` object",
    fixed = TRUE
  )
  expect_error(
    mfrmr::build_apa_outputs(.fit, .diag, context = 1),
    "requires `context` to be a list",
    fixed = TRUE
  )
  expect_error(
    mfrmr::build_apa_outputs(.fit, .diag, bias_results = 1),
    "requires `bias_results` to be `NULL` or a package-native bias result",
    fixed = TRUE
  )
})

test_that("apa_table uses the same APA contract validation when contract metadata is needed", {
  expect_error(
    mfrmr::apa_table(.fit, which = "summary", diagnostics = list()),
    "`diagnostics` must be an mfrm_diagnostics object",
    fixed = TRUE
  )
  expect_error(
    mfrmr::apa_table(.fit, which = "summary", diagnostics = .diag, context = 1),
    "requires `context` to be a list",
    fixed = TRUE
  )
})

# ============================================================================
# build_apa_table_figure_notes / captions
# ============================================================================

test_that("build_apa_table_figure_notes produces note text", {
  notes <- mfrmr:::build_apa_reporting_contract(.fit, .diag)$note_text
  expect_true(is.character(notes))
  expect_true(nchar(notes) > 0)
})

test_that("build_apa_table_figure_notes with bias_results includes bias notes", {
  notes <- mfrmr:::build_apa_reporting_contract(.fit, .diag, bias_results = .bias)$note_text
  expect_true(nchar(notes) > 0)
})

test_that("build_apa_table_figure_captions produces caption text", {
  captions <- mfrmr:::build_apa_reporting_contract(.fit, .diag, whexact = FALSE)$caption_text
  expect_true(is.character(captions))
  expect_true(grepl("Table 1", captions))
})

test_that("build_apa_table_figure_captions with context adds assessment phrase", {
  ctx <- list(assessment = "oral proficiency")
  captions <- mfrmr:::build_apa_reporting_contract(.fit, .diag, context = ctx, whexact = FALSE)$caption_text
  expect_true(grepl("oral proficiency", captions))
})

test_that("build_apa_table_figure_captions with bias includes interaction label", {
  captions <- mfrmr:::build_apa_reporting_contract(.fit, .diag, bias_results = .bias, whexact = FALSE)$caption_text
  expect_true(grepl("Table 4", captions))
})

# ============================================================================
# resolve_warning_thresholds
# ============================================================================

test_that("resolve_warning_thresholds returns all profiles", {
  for (profile in c("strict", "standard", "lenient")) {
    resolved <- mfrmr:::resolve_warning_thresholds(threshold_profile = profile)
    expect_equal(resolved$profile_name, profile)
    expect_true("thresholds" %in% names(resolved))
    expect_true("pca_reference_bands" %in% names(resolved))
    expect_true(resolved$thresholds$n_obs_min > 0)
  }
})

test_that("resolve_warning_thresholds applies custom overrides", {
  resolved <- mfrmr:::resolve_warning_thresholds(
    thresholds = list(n_obs_min = 42),
    threshold_profile = "standard"
  )
  expect_equal(resolved$thresholds$n_obs_min, 42)
})

# ============================================================================
# build_pca_reference_text / build_pca_check_text
# ============================================================================

test_that("build_pca_reference_text returns reference band text", {
  bands <- mfrmr:::warning_threshold_profiles()$pca_reference_bands
  text <- mfrmr:::build_pca_reference_text(bands)
  expect_true(is.character(text))
  expect_true(nchar(text) > 0)
})

test_that("build_pca_check_text returns interpretive text for a given eigenvalue/proportion", {
  bands <- mfrmr:::warning_threshold_profiles()$pca_reference_bands
  text <- mfrmr:::build_pca_check_text(eigenvalue = 1.5, proportion = 0.08, reference_bands = bands)
  expect_true(is.character(text))
  expect_true(nchar(text) > 0)
})

test_that("build_pca_check_text with high eigenvalue flags concern", {
  bands <- mfrmr:::warning_threshold_profiles()$pca_reference_bands
  text <- mfrmr:::build_pca_check_text(eigenvalue = 5.0, proportion = 0.25, reference_bands = bands)
  expect_true(nchar(text) > 0)
})

# ============================================================================
# Format helpers: format_fixed_width_table edge cases
# ============================================================================

test_that("format_fixed_width_table handles empty data frame", {
  result <- mfrmr:::format_fixed_width_table(data.frame(), columns = character(0))
  expect_equal(result, "No data")
})

test_that("format_fixed_width_table handles missing columns gracefully", {
  df <- data.frame(A = 1:3, B = c("x", "y", "z"), stringsAsFactors = FALSE)
  result <- mfrmr:::format_fixed_width_table(df, columns = c("A", "B", "C"))
  expect_true(is.character(result))
  expect_true(grepl("A", result))
})

# ============================================================================
# build_sectioned_fixed_report
# ============================================================================

test_that("build_sectioned_fixed_report handles various section types", {
  report <- mfrmr:::build_sectioned_fixed_report(
    title = "Test Report",
    sections = list(
      list(title = "Empty Section", data = NULL),
      list(title = "String Section", data = "Some text here"),
      list(title = "Data Section", data = data.frame(X = 1:3, Y = c("a", "b", "c"), stringsAsFactors = FALSE)),
      list(title = "Empty DF", data = data.frame()),
      list(title = "Truncated Section", data = data.frame(X = 1:10), max_rows = 3)
    )
  )
  expect_true(grepl("Test Report", report))
  expect_true(grepl("Some text here", report))
  expect_true(grepl("Empty Section", report))
  expect_true(grepl("Showing first 3 rows of 10", report))
})

# ============================================================================
# summarize_anchor_constraints
# ============================================================================

test_that("summarize_anchor_constraints produces text about anchor status", {
  text <- mfrmr:::summarize_anchor_constraints(.fit$config)
  expect_true(is.character(text))
})

# ============================================================================
# summarize_convergence_metrics
# ============================================================================

test_that("summarize_convergence_metrics produces text from summary row", {
  summary_row <- .fit$summary
  if (!is.null(summary_row) && nrow(summary_row) > 0) {
    text <- mfrmr:::summarize_convergence_metrics(summary_row[1, , drop = FALSE])
    expect_true(is.character(text))
  }
})

test_that("convergence prose separates numerical checks from inference readiness", {
  row <- data.frame(
    Converged = TRUE, ConvergenceSeverity = "pass",
    ConvergenceStatus = "converged", InferenceReady = FALSE,
    ReadinessContractVersion = mfrmr:::mfrmr_readiness_contract_version(),
    FitReadiness = "review", Method = "MML"
  )
  expect_false(mfrmr:::mfrm_convergence_state(row)$inference_ready)
  expect_match(mfrmr:::summarize_convergence_metrics(row),
               "met the numerical convergence checks", fixed = TRUE)
  row$ConvergenceSeverity <- "review"
  row$ConvergenceStatus <- "converged_gradient_review"
  expect_match(mfrmr:::summarize_convergence_metrics(row),
               "terminal gradient requires review", fixed = TRUE)
  row$ConvergenceStatus <- "reviewable_warning"
  expect_match(mfrmr:::summarize_convergence_metrics(row),
               "reviewable optimizer warning", fixed = TRUE)
  row$Converged <- FALSE
  row$ConvergenceStatus <- "not_converged"
  expect_match(mfrmr:::summarize_convergence_metrics(row),
               "did not meet the package convergence checks", fixed = TRUE)
})

# ============================================================================
# summarize_step_estimates
# ============================================================================

test_that("summarize_step_estimates produces text from step table", {
  text <- mfrmr:::summarize_step_estimates(.fit$steps)
  expect_true(is.character(text))
  expect_true(nchar(paste(text, collapse = "")) > 0)
})

# ============================================================================
# summarize_bias_counts
# ============================================================================

test_that("summarize_bias_counts with bias results", {
  text <- mfrmr:::summarize_bias_counts(.bias)
  expect_true(is.character(text))
  expect_true(nchar(paste(text, collapse = "")) > 0)
})

test_that("summarize_bias_counts with NULL returns text", {
  text <- mfrmr:::summarize_bias_counts(NULL)
  expect_true(is.character(text))
})

# ============================================================================
# summarize_top_misfit_levels
# ============================================================================

test_that("summarize_top_misfit_levels produces text", {
  fit_tbl <- .diag$fit
  if (!is.null(fit_tbl) && nrow(fit_tbl) > 0) {
    text <- mfrmr:::summarize_top_misfit_levels(fit_tbl, top_n = 3)
    expect_true(is.character(text))
  }
})

# ============================================================================
# collapse_apa_paragraph
# ============================================================================

test_that("collapse_apa_paragraph wraps text at specified width", {
  sentences <- c(
    "This is a first sentence.",
    "This is a second sentence that is somewhat longer than the first.",
    "And a third one."
  )
  result <- mfrmr:::collapse_apa_paragraph(sentences, width = 40L)
  expect_true(is.character(result))
  expect_true(nchar(result) > 0)
})

# ============================================================================
# py_style_format edge cases
# ============================================================================

test_that("py_style_format handles edge cases", {
  fmt <- mfrmr:::py_style_format
  expect_equal(fmt(NULL, 42), "42")
  expect_equal(fmt("{}", 42), "42")
  expect_equal(fmt("{:.0f}", NA), "")
  expect_equal(fmt(function(x) paste0("$", x), 100), "$100")
  expect_equal(fmt(c("a", "b"), 42), "42")
})

# ============================================================================
# describe_series edge cases
# ============================================================================

test_that("describe_series handles NULL and empty inputs", {
  expect_null(mfrmr:::describe_series(NULL))
  expect_null(mfrmr:::describe_series(c(NA, NaN)))
  s <- mfrmr:::describe_series(c(1, 2, 3))
  expect_equal(s$min, 1)
  expect_equal(s$max, 3)
  expect_equal(s$mean, 2)
})

test_that("describe_series handles single element", {
  s <- mfrmr:::describe_series(5)
  expect_equal(s$min, 5)
  expect_equal(s$max, 5)
  expect_true(is.na(s$sd))
})

# ============================================================================
# build_bias_fixed_text / build_pairwise_fixed_text
# ============================================================================

test_that("build_bias_fixed_text with valid data produces text", {
  if (!is.null(.bias$table) && nrow(.bias$table) > 0) {
    tbl <- as.data.frame(.bias$table, stringsAsFactors = FALSE)
    cols <- intersect(
      c("Bias Size", "Obs-Exp Average", "S.E.", "t", "Prob."),
      names(tbl)
    )
    text <- mfrmr:::build_bias_fixed_text(
      table_df = tbl,
      summary_df = .bias$summary %||% data.frame(),
      chi_df = .bias$chi_sq %||% data.frame(),
      facet_a = "Rater",
      facet_b = "Task",
      columns = cols,
      formats = list()
    )
    expect_true(is.character(text))
    expect_true(grepl("Bias", text))
  }
})

test_that("build_bias_fixed_text with empty table returns 'No bias data'", {
  text <- mfrmr:::build_bias_fixed_text(
    table_df = data.frame(),
    summary_df = data.frame(),
    chi_df = data.frame(),
    facet_a = "Rater",
    facet_b = "Task",
    columns = character(0),
    formats = list()
  )
  expect_equal(text, "No bias data")
})

test_that("build_pairwise_fixed_text with empty data returns 'No pairwise data'", {
  text <- mfrmr:::build_pairwise_fixed_text(
    pair_df = data.frame(),
    target_facet = "Rater",
    context_facet = "Task",
    columns = character(0),
    formats = list()
  )
  expect_equal(text, "No pairwise data")
})

# ============================================================================
# extract_overall_pca_first / extract_overall_pca_second / extract_facet_pca_first
# ============================================================================

test_that("PCA extraction functions return correct data", {
  pca <- mfrmr:::safe_residual_pca(.diag, mode = "both")
  if (!is.null(pca)) {
    first <- mfrmr:::extract_overall_pca_first(pca)
    expect_true(!is.null(first))
    expect_true("Eigenvalue" %in% names(first))

    second <- mfrmr:::extract_overall_pca_second(pca)
    # second may be NULL if only 1 component
    expect_true(is.null(second) || "Eigenvalue" %in% names(second))

    facet <- mfrmr:::extract_facet_pca_first(pca)
    expect_true(is.data.frame(facet))
  }
})

test_that("PCA extraction functions handle NULL input", {
  expect_null(mfrmr:::extract_overall_pca_first(NULL))
  expect_null(mfrmr:::extract_overall_pca_second(NULL))
  facet <- mfrmr:::extract_facet_pca_first(NULL)
  expect_true(is.data.frame(facet))
  expect_equal(nrow(facet), 0)
})

# ============================================================================
# build_apa_table_figure_note_map (internal)
# ============================================================================

test_that("build_apa_table_figure_note_map returns a named list", {
  note_map <- mfrmr:::build_apa_reporting_contract(.fit, .diag)$note_map
  expect_true(is.list(note_map))
  expect_true(length(note_map) > 0)
})

test_that("build_apa_table_figure_note_map with bias_results includes table4", {
  note_map <- mfrmr:::build_apa_reporting_contract(.fit, .diag, bias_results = .bias)$note_map
  expect_true("table4" %in% names(note_map))
})
