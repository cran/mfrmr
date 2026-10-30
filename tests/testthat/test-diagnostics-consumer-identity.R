diagnostic_consumer_fixture <- local({
  cases <- readRDS(test_path("fixtures", "mml-comparison-intervals.rds"))
  fit <- cases[["RSM-regression"]]$fit
  list(fit = fit, other = cases$RSM$fit,
       diagnostics = suppressWarnings(diagnose_mfrm(
         fit, diagnostic_mode = "legacy", residual_pca = "none")))
})

test_that("direct reports reject another fit's diagnostics with the same labels", {
  x <- diagnostic_consumer_fixture
  expect_identical(x$fit$facets$others[c("Facet", "Level")],
                   x$other$facets$others[c("Facet", "Level")])
  local_mocked_bindings(diagnose_mfrm = function(...) stop("Unexpected recomputation"),
                       .package = "mfrmr")
  consumers <- list(facet_quality_dashboard, facet_statistics_report, fit_measures_table,
                    facets_chisq_table, displacement_table, measurable_summary_table,
                    facet_small_sample_review,
                    function(fit, diagnostics) estimate_bias(fit, diagnostics,
                      facet_a = "Rater", facet_b = "Criterion", max_iter = 1),
                    function(fit, diagnostics) estimate_all_bias(fit, diagnostics,
                      pairs = list(c("Rater", "Criterion")), max_iter = 1),
                    function(fit, diagnostics) facets_output_file_bundle(fit,
                      diagnostics, include = "score", score_se_method = "none"))
  for (consume in consumers) {
    expect_error(consume(x$other, diagnostics = x$diagnostics), "fit/diagnostics mismatch")
    expect_no_error(consume(x$fit, diagnostics = x$diagnostics))
  }
  # Validation precedes the optional refresh of FACETS-style degrees of freedom.
  legacy <- x$diagnostics
  legacy$measures[c("DF_Infit_FACETS", "DF_Outfit_FACETS",
                    "InfitZSTD_FACETS", "OutfitZSTD_FACETS")] <- NULL
  expect_error(fit_measures_table(x$other, diagnostics = legacy,
                                 fit_df_method = "both"), "fit/diagnostics mismatch")
  expect_error(facets_fit_review(x$other, diagnostics = legacy), "fit/diagnostics mismatch")
  expect_error(analyze_dff(x$other, x$diagnostics, facet = "Rater", group = "Group",
                         method = "refit"), "fit/diagnostics mismatch")
  # Supplied secondary tables must not bypass the primary diagnostic check.
  expect_error(build_misfit_casebook(x$other, x$diagnostics,
    unexpected = structure(list(), class = "mfrm_unexpected"),
    displacement = structure(list(), class = "mfrm_displacement")), "fit/diagnostics mismatch")
  expect_error(facets_output_contract_review(x$other, x$diagnostics), "fit/diagnostics mismatch")
  for (type in c("qc", "residual_qq", "guttman", "reliability", "fit_pathway")) {
    expect_error(plot(x$other, type = type, diagnostics = x$diagnostics, draw = FALSE),
                 "fit/diagnostics mismatch")
  }
})

test_that("CSV export validates diagnostics before creating any output", {
  x <- diagnostic_consumer_fixture
  directory <- tempfile("diagnostic-identity-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  expect_error(export_mfrm_tables(x$other, diagnostics = x$diagnostics,
    output_dir = directory, acknowledge_sensitive = TRUE), "fit/diagnostics mismatch")
  expect_false(dir.exists(directory))

  expect_error(export_mfrm_bundle(x$other, diagnostics = x$diagnostics,
    output_dir = directory, include = "core_tables", acknowledge_sensitive = TRUE),
    "fit/diagnostics mismatch")
  expect_false(dir.exists(directory))
  expect_error(facets_output_file_bundle(x$other, diagnostics = x$diagnostics,
    output_dir = directory, write_files = TRUE), "fit/diagnostics mismatch")
  expect_false(dir.exists(directory))
  for (writer in list(write_mfrm_residual_file, write_mfrm_subset_file)) {
    path <- tempfile(fileext = ".csv")
    expect_error(writer(x$other, path = path, diagnostics = x$diagnostics),
                 "fit/diagnostics mismatch")
    expect_false(file.exists(path))
    expect_no_error(writer(x$fit, path = path, diagnostics = x$diagnostics))
    expect_true(file.exists(path))
    unlink(path)
  }

  files <- export_mfrm_tables(x$fit, diagnostics = x$diagnostics,
    output_dir = directory, tables = c("facets", "measures"), acknowledge_sensitive = TRUE)
  expect_true(all(file.exists(files$Path)))
  facets <- read.csv(files$Path[files$Table == "facets"])
  source <- x$diagnostics$measures
  key <- function(x) paste(x$Facet, x$Level, sep = "\r")
  idx <- match(key(facets), key(source))
  expect_equal(facets$Estimate, unname(source$Estimate[idx]), tolerance = 1e-12)
  expect_equal(facets$Infit, unname(source$Infit[idx]), tolerance = 1e-12)
})

test_that("diagnostic identity includes calibration beyond location estimates", {
  x <- diagnostic_consumer_fixture
  local_mocked_bindings(diagnose_mfrm = function(...) stop("Unexpected recomputation"),
                       .package = "mfrmr")
  variants <- list(
    function(f) { f$steps$Estimate[1] <- f$steps$Estimate[1] + .25; f },
    function(f) { f$opt$par[1] <- f$opt$par[1] + .25; f },
    function(f) { f$population$sigma2 <- f$population$sigma2 * 2; f },
    function(f) { f$config$estimation_control$quad_points <- f$config$estimation_control$quad_points + 2L; f },
    function(f) { f$config$population_spec$design_matrix[1, 2] <- 99; f }
  )
  for (change in variants) {
    expect_error(mfrm_results_validate_diagnostics_identity(change(x$fit), x$diagnostics),
                 "fit/diagnostics mismatch")
  }
  for (change in list(
    function(d) {
      names(d$parameter_uncertainty$source$constraints)[1] <- "different_constraint"
      d
    },
    function(d) { d$measures <- rbind(d$measures, d$measures[1, ]); d },
    function(d) { d$measures$Level[1] <- NA_character_; d },
    function(d) { d$obs$Weight <- NULL; d },
    function(d) { d$obs <- d$obs[nrow(d$obs):1, ]; d }
  )) {
    expect_error(mfrm_results_validate_diagnostics_identity(x$fit, change(x$diagnostics)),
                 "fit/diagnostics mismatch")
  }
  old <- x$diagnostics
  old$parameter_uncertainty$source <- NULL
  expect_error(fit_measures_table(x$fit, diagnostics = old), "no model refit is needed")
  expect_no_error(fit_measures_table(old))
})

test_that("matching saved diagnostics remain usable and preserve interval eligibility", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(diagnostic_consumer_fixture, path)
  x <- readRDS(path)
  local_mocked_bindings(diagnose_mfrm = function(...) stop("Unexpected recomputation"),
                       .package = "mfrmr")
  dashboard <- facet_quality_dashboard(x$fit, diagnostics = x$diagnostics)
  source <- x$diagnostics$measures
  selected <- source$Facet == dashboard$facet
  expect_equal(dashboard$detail$Estimate, unname(source$Estimate[selected]))
  measures <- fit_measures_table(x$fit, diagnostics = x$diagnostics)
  idx <- match(paste(measures$table$Facet, measures$table$Level),
               paste(source$Facet, source$Level))
  expect_identical(measures$table$CIEligible, source$CIEligible[idx])
  expect_no_error(facet_statistics_report(x$fit, diagnostics = x$diagnostics))
})
