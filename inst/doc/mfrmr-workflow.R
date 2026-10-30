## ----include = FALSE--------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
vignette_artifact <- function(name) {
  path <- system.file(
    "extdata", "vignette-artifacts", name,
    package = "mfrmr",
    mustWork = TRUE
  )
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character(),
    colClasses = "character"
  )
}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = !is_cran_check
)

## ----quick-start, fig.alt = "Wright map of fitted person abilities, rater severities, criterion difficulties and shared category steps for the synthetic ratings. All locations use the same logit scale."----
# # Load the package
# library(mfrmr)
# 
# # Load example ratings and look at the first six rows
# ratings <- load_mfrmr_data("example_operational")
# head(ratings)
# 
# # Fit the model
# fit <- fit_mfrm(
#   data = ratings,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   rating_min = 1, rating_max = 4,
#   category_policy = "preserve",
#   method = "MML",
#   model = "RSM"
# )
# 
# # Plot the results (Wright map)
# plot(fit)
# 
# # Save the summary, then display its tables
# results <- summary(fit)
# results$person_overview # One row summarizing person ability estimates
# results$facet_overview  # One row per facet: number of levels, mean, SD, range
# 
# # Check the interpretation status and recommended next step
# results$decision

## ----quick-start-estimates--------------------------------------------------
# estimates <- as.data.frame(fit)
# head(subset(estimates, Facet == "Person")) # First six persons
# subset(estimates, Facet == "Rater")       # All raters
# subset(estimates, Facet == "Criterion")   # All criteria

## ----quick-start-diagnostics------------------------------------------------
# diagnostics <- diagnose_mfrm(fit)
# diagnostic_summary <- summary(diagnostics)
# diagnostic_summary$decision

## ----csv-practice-file------------------------------------------------------
# library(mfrmr)
# csv_path <- tempfile(fileext = ".csv")
# write.csv(load_mfrmr_data("example_operational"), csv_path,
#           row.names = FALSE, na = "")

## ----csv-select-file, eval = FALSE------------------------------------------
# library(mfrmr)
# csv_path <- file.choose()

## ----csv-read---------------------------------------------------------------
# ratings <- read.csv(
#   csv_path,
#   colClasses = "character",
#   na.strings = "",
#   check.names = FALSE,
#   fileEncoding = "UTF-8-BOM"
# )
# # Treat the documented missing-score marker only in the score column
# ratings <- recode_missing_codes(ratings, columns = "Score", codes = "NA")
# head(ratings)
# names(ratings)
# table(ratings$Score, useNA = "ifany")

## ----csv-data-review--------------------------------------------------------
# data_review <- describe_mfrm_data(
#   data = ratings,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   rating_min = 1,
#   rating_max = 4,
#   category_policy = "preserve"
# )
# data_review$row_retention
# data_review$missing_by_column
# data_review$score_distribution
# data_review$design_connectivity

## ----csv-fit----------------------------------------------------------------
# csv_fit <- fit_mfrm(
#   data = ratings,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   rating_min = 1,
#   rating_max = 4,
#   category_policy = "preserve",
#   method = "MML",
#   model = "RSM"
# )
# csv_diagnostics <- diagnose_mfrm(csv_fit)
# csv_results <- summary(csv_fit, diagnostics = csv_diagnostics)
# csv_results$decision

## ----csv-results, fig.alt = "Wright map for the ratings read from the example CSV. Compare person, facet and step locations on the fitted scale; the display does not establish that model assumptions hold."----
# plot(csv_fit)
# csv_estimates <- as.data.frame(csv_fit)
# head(subset(csv_estimates, Facet == "Person"))
# subset(csv_estimates, Facet == "Rater")
# subset(csv_estimates, Facet == "Criterion")

## ----csv-missing-codes, eval = FALSE----------------------------------------
# ratings <- recode_missing_codes(
#   ratings, columns = "Score", codes = c("99", ".")
# )
# table(ratings$Score, useNA = "ifany")

## ----csv-wide-to-long-------------------------------------------------------
# wide <- data.frame(
#   Person = c("001", "001"),
#   Rater = c("R1", "R2"),
#   Content = c(3, 4),
#   Style = c(2, 3)
# )
# ratings_long <- tidyr::pivot_longer(
#   wide,
#   cols = c("Content", "Style"),
#   names_to = "Criterion",
#   values_to = "Score"
# )
# ratings_long

## ----reporting-input--------------------------------------------------------
# analysis_data <- ratings
# fit_toy <- fit
# diag_toy <- diagnostics
# data("mfrmr_example_operational_design", package = "mfrmr")
# analysis_roster <- mfrmr_example_operational_design

## ----reporting-csv-input, eval = FALSE--------------------------------------
# analysis_data <- ratings
# fit_toy <- csv_fit
# diag_toy <- csv_diagnostics
# analysis_roster <- NULL # Replace with your own planned roster if available

## ----ratings-setup----------------------------------------------------------
# data_review_toy <- describe_mfrm_data(
#   data = analysis_data,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   rating_min = 1,
#   rating_max = 4,
#   category_policy = "preserve",
#   expected_design = analysis_roster
# )
# data_summary_toy <- summary(data_review_toy)
# data_summary_toy$structural_missingness
# data_summary_toy$design_connectivity
# 
# fit_summary_toy <- summary(fit_toy)
# fit_summary_toy$overview[, c("Model", "Method", "N", "Persons", "Converged")]
# fit_summary_toy$decision

## ----ratings-setup-artifacts, echo = FALSE, eval = is_cran_check------------
vignette_artifact("workflow_fit_decision.csv")
vignette_artifact("workflow_fit_overview.csv")

## ----precision-and-fit, fig.alt = "Fitted person, facet and step locations with the requested uncertainty display. Interpret precision using the accompanying table and its stated interval basis."----
# diagnostic_summary_toy <- summary(diag_toy)
# diagnostic_summary_toy$decision
# diagnostic_summary_toy$overall_fit
# 
# precision_toy <- precision_review_report(fit_toy, diagnostics = diag_toy)
# precision_toy$profile[, c(
#   "PrecisionTier", "SupportsFormalInference", "PersonSEBasis", "NonPersonSEBasis"
# )]
# precision_toy$checks
# 
# # View the fitted locations with their uncertainty
# plot(fit_toy, diagnostics = diag_toy, show_ci = TRUE)

## ----diagnostic-overview-artifact, echo = FALSE, eval = is_cran_check-------
vignette_artifact("workflow_diagnostic_overview.csv")

## ----category-review, fig.alt = "Category-use and threshold summaries for the fitted rubric. Sparse categories and transition ordering should be reviewed together with the rating design."----
# scale_toy <- rating_scale_table(fit_toy, diagnostics = diag_toy)
# scale_toy$category_table[, c("Category", "Count", "AvgMeasure", "ExpectedMeasure", "Infit", "Outfit")]
# scale_toy$model_context
# scale_toy$owner_category_table[, c(
#   "StepFacet", "CategoryLabel", "RawCount", "Count", "Percent", "UsageStatus"
# )]
# scale_toy$threshold_table[, c("StepFacet", "Step", "Estimate", "SE", "CIUse")]
# plot(scale_toy)
# scale_toy$notes

## ----diagnostics-reporting, fig.alt = "Unusual rating observations selected by the stated residual and probability thresholds. Flags support review of the underlying records, not automatic removal."----
# unexpected_toy <- unexpected_response_table(
#   fit_toy,
#   diagnostics = diag_toy,
#   abs_z_min = 1.5,
#   prob_max = 0.4,
#   top_n = 10
# )
# unexpected_toy$summary
# unexpected_toy$table[, c("Person", "Rater", "Criterion", "Observed", "Expected", "StdResidual")]
# plot(unexpected_toy)

## ----diagnostic-bases, fig.alt = "Marginal fit display based on posterior-averaged response probabilities. Its probability basis differs from conditional residual diagnostics."----
# diagnostic_summary_toy$diagnostic_basis[, c("DiagnosticPath", "Status", "ReportingUse")]
# plot(fit_toy, type = "marginal_fit", diagnostics = diag_toy)

## ----rater-feedback, fig.alt = "Rater severity estimates with 95% normal intervals on the fitted logit scale. Positive estimates indicate stricter ratings relative to the fitted reference."----
# plot(fit_toy, type = "rater_severity", diagnostics = diag_toy, facet = "Rater", show_bands = FALSE)
# rater_review <- facet_quality_dashboard(
#   fit_toy, diagnostics = diag_toy, facet = "Rater"
# )
# knitr::kable(
#   rater_review$detail[, c("Level", "N", "Estimate", "SE", "Infit", "Outfit")],
#   digits = 3, row.names = FALSE, caption = "Rater estimates and fit"
# )
# knitr::kable(
#   rater_review$detail[, c("Level", "MissingMetrics", "SeverityFlag", "MisfitFlag")],
#   row.names = FALSE, caption = "Availability and screening flags"
# )
# rater_review$settings
# writeLines(strwrap(rater_review$notes, width = 72))

## ----primary-route----------------------------------------------------------
# res_toy <- mfrm_results(fit_toy, diagnostics = diag_toy, include = "publication")
# report_toy <- mfrm_report(res_toy, style = "apa")
# report_toy$first_screen[, c("Area", "Status", "MainIssue", "NextAction")]
# 
# checklist_toy <- reporting_checklist(fit_toy, diagnostics = diag_toy)
# subset(checklist_toy$checklist, !DraftReady,
#        c("Section", "Item", "NextAction"))

## ----report-overview-artifact, echo = FALSE, eval = is_cran_check-----------
vignette_artifact("workflow_next_actions.csv")

## ----report-table-----------------------------------------------------------
# measurements_toy <- fit_measures_table(fit_toy, diagnostics = diag_toy)
# rater_table_toy <- subset(
#   measurements_toy$table, Facet == "Rater",
#   c("Level", "N", "Measure", "SE", "CI_Lower", "CI_Upper", "Infit", "Outfit")
# )
# apa_table(rater_table_toy, digits = 3,
#           caption = "Rater severity and fit",
#           note = "N counts rating rows. Rater effects are centered within facet and measured in logits; higher values indicate stricter ratings. Intervals are 95% normal approximations from MML observed-information SEs. Infit and Outfit have reference value 1.")

## ----report-export----------------------------------------------------------
# # A new temporary directory keeps this synthetic example repeatable
# export_dir <- tempfile("mfrmr-workflow-export-")
# export_toy <- export_mfrm_results(
#   res_toy,
#   output_dir = export_dir,
#   include = c("default", "report"),
#   acknowledge_sensitive = TRUE
# )
# export_preview <- head(export_toy$written_files[, c("Component", "Format", "Path")])
# export_preview$Path <- basename(export_preview$Path)
# knitr::kable(export_preview, row.names = FALSE,
#              caption = "First six exported files (filenames only)")

## ----report-reload----------------------------------------------------------
# results_file <- subset(export_toy$written_files, Component == "results_rds")$Path
# restored_results <- readRDS(results_file)
# summary(restored_results)$decision

## ----export-artifact, echo = FALSE, eval = is_cran_check--------------------
vignette_artifact("workflow_export_files.csv")

## ----model-export-example---------------------------------------------------
# # Reuse the fit above; no additional estimation or backend is needed to export.
# stan_model <- export_mfrm_model(
#   fit_toy, tempfile("stan-model-"), format = "stan",
#   acknowledge_sensitive = TRUE
# )
# brms_model <- export_mfrm_model(
#   fit_toy, tempfile("brms-model-"), format = "brms",
#   acknowledge_sensitive = TRUE
# )
# names(stan_model$files)
# model_data <- readRDS(stan_model$files[["model-data.rds"]])
# model_data$categories
# model_data$priors

## ----saved-posterior-return-------------------------------------------------
# posterior <- read_mfrm_posterior(system.file(
#   "extdata", "posterior-example", package = "mfrmr"))
# summary(posterior)
# posterior_plot <- plot(posterior, quantity = "location", draw = FALSE)
# posterior_report <- mfrm_report(posterior, style = "apa")

## ----simulation-file-workflow, eval=FALSE-----------------------------------
# dir.create("results", showWarnings = FALSE)
# recovery <- simulate_mfrm_recovery(
#   n_person = 50, n_rater = 4, n_criterion = 4,
#   reps = 100, seed = 20261021, file = "results/recovery.rds"
# )
# summary(recovery)$rep_overview  # Retains failed attempts, too.

## ----observed-writing-example, eval = FALSE---------------------------------
# library(mfrmr)
# data("data.ratings1", package = "sirt")
# wide <- data.ratings1
# persons <- match(wide$idstud, unique(wide$idstud))
# raters <- match(wide$rater, unique(wide$rater))
# writing <- do.call(rbind, lapply(paste0("k", 1:5), function(criterion) {
#   data.frame(
#     Person = sprintf("P%03d", persons), Rater = paste0("R", raters),
#     Criterion = criterion, Score = wide[[criterion]]
#   )
# }))
# review <- describe_mfrm_data(
#   writing, person = "Person", facets = c("Rater", "Criterion"), score = "Score",
#   rating_min = 0, rating_max = 3
# )
# review$row_retention
# review$design_connectivity
# 
# writing_fit <- fit_mfrm(
#   writing, "Person", c("Rater", "Criterion"), "Score",
#   method = "MML", model = "PCM", step_facet = "Criterion",
#   rating_min = 0, rating_max = 3, quad_points = 61, maxit = 400
# )
# writing_diagnostics <- diagnose_mfrm(writing_fit, residual_pca = "none")
# feedback <- facet_quality_dashboard(
#   writing_fit, diagnostics = writing_diagnostics, facet = "Rater"
# )
# feedback$detail[, c("Level", "N", "Estimate", "SE", "Infit", "Outfit", "MisfitFlag")]
# feedback$settings
# plot(writing_fit, type = "rater_severity", diagnostics = writing_diagnostics, facet = "Rater", show_bands = FALSE)
# saveRDS(feedback, "writing-feedback.rds")
# archive <- export_mfrm_bundle(
#   writing_fit, diagnostics = writing_diagnostics, facet = "Rater",
#   output_dir = "writing-analysis", data = writing, acknowledge_sensitive = TRUE
# )

