## ----include = FALSE----------------------------------------------------------
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

## ----load-data----------------------------------------------------------------
# library(mfrmr)
# 
# list_mfrmr_data(details = TRUE)[, c("Key", "PrimaryUse", "Design", "CountBasis")]
# 
# data("ej2021_study1", package = "mfrmr")
# head(ej2021_study1)
# 
# study1_alt <- load_mfrmr_data("study1")
# identical(names(ej2021_study1), names(study1_alt))

## ----toy-setup----------------------------------------------------------------
# data("mfrmr_example_operational", package = "mfrmr")
# data("mfrmr_example_operational_design", package = "mfrmr")
# toy <- mfrmr_example_operational
# 
# data_review_toy <- describe_mfrm_data(
#   data = toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   rating_min = 1,
#   rating_max = 4,
#   expected_design = mfrmr_example_operational_design
# )
# data_summary_toy <- summary(data_review_toy)
# data_summary_toy$structural_missingness
# data_summary_toy$design_connectivity
# 
# fit_toy <- fit_mfrm(
#   data = toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM"
# )
# diag_toy <- diagnose_mfrm(
#   fit_toy,
#   residual_pca = "none",
#   diagnostic_mode = "both",
#   fit_df_method = "both"
# )
# 
# # Fast fit-only summary: this does not compute diagnostics.
# fit_summary_toy <- summary(fit_toy, profile = "fit", detail = "brief")
# 
# # Comprehensive review, reusing diagnostics already computed above.
# facets_summary_toy <- summary(
#   fit_toy,
#   profile = "facets",
#   detail = "brief",
#   diagnostics = diag_toy
# )
# res_toy <- facets_summary_toy$results
# 
# fit_summary_toy$overview
# fit_summary_toy$readiness
# fit_summary_toy$data_review
# summary(diag_toy)$overview
# facets_summary_toy
# 
# # Required first fitted-scale figure.
# plot(res_toy, type = "wright", preset = "publication", show_ci = TRUE, top_n = Inf)
# 
# # Optional FACETS-style ruler. Replace these examples with the study rubric.
# rubric_labels <- c(
#   "1" = "Level 1",
#   "2" = "Level 2",
#   "3" = "Level 3",
#   "4" = "Level 4"
# )
# plot(
#   res_toy,
#   type = "wright",
#   renderer = "facets",
#   category_labels = rubric_labels,
#   show_ci = FALSE,
#   preset = "publication"
# )
# # Setting show_ci = TRUE on this FACETS-style ruler is available as a
# # deliberate hybrid: FACETS ruler grammar plus mfrmr uncertainty intervals.
# 
# # Optional follow-up: Infit on x, measure on y; persons are explicit opt-in.
# plot(
#   res_toy,
#   type = "fit_pathway",
#   fit_stat = "Infit",
#   include_person = TRUE,
#   top_n_person = 12,
#   person_labels = "none",
#   facet_labels = "flagged",
#   preset = "publication"
# )

## ----summary-profile-controls, eval = FALSE-----------------------------------
# reporting_summary_toy <- summary(
#   fit_toy,
#   profile = "reporting",
#   diagnostics = diag_toy
# )
# 
# availability_only <- summary(
#   fit_toy,
#   profile = "facets",
#   compute = "never"
# )
# availability_only$section_status

## ----toy-setup-artifacts, echo = FALSE, eval = is_cran_check------------------
vignette_artifact("workflow_fit_overview.csv")
vignette_artifact("workflow_diagnostic_overview.csv")
vignette_artifact("workflow_plot_components.csv")

## ----primary-route------------------------------------------------------------
# report_toy <- mfrm_report(res_toy, style = "qc")
# 
# summary(res_toy)$next_actions
# summary(report_toy)$overview
# 
# # This is a controlled analysis archive, not a deidentified shareable export.
# export_dir <- file.path(tempdir(), "mfrmr-workflow-export")
# export_toy <- export_mfrm_results(
#   res_toy,
#   output_dir = export_dir,
#   include = c("default", "report"),
#   overwrite = TRUE,
#   acknowledge_sensitive = TRUE
# )
# head(export_toy$written_files)

## ----primary-route-artifacts, echo = FALSE, eval = is_cran_check--------------
vignette_artifact("workflow_next_actions.csv")
vignette_artifact("workflow_report_overview.csv")
vignette_artifact("workflow_export_files.csv")

## ----diagnostics-reporting----------------------------------------------------
# t4_toy <- unexpected_response_table(
#   fit_toy,
#   diagnostics = diag_toy,
#   abs_z_min = 1.5,
#   prob_max = 0.4,
#   top_n = 10
# )
# t12_toy <- fair_average_table(fit_toy, diagnostics = diag_toy)
# t13_toy <- bias_interaction_report(
#   estimate_bias(fit_toy, diag_toy,
#                 facet_a = "Rater", facet_b = "Criterion",
#                 max_iter = 2),
#   top_n = 10
# )
# 
# class(summary(t4_toy))
# class(summary(t12_toy))
# class(summary(t13_toy))
# 
# names(plot(t4_toy, draw = FALSE))
# names(plot(t12_toy, draw = FALSE))
# names(plot(t13_toy, draw = FALSE))
# 
# chk_toy <- reporting_checklist(fit_toy, diagnostics = diag_toy)
# subset(
#   chk_toy$checklist,
#   Section == "Visual Displays",
#   c("Item", "DraftReady", "NextAction")
# )

## ----diagnostics-reporting-artifacts, echo = FALSE, eval = is_cran_check------
vignette_artifact("workflow_summary_classes.csv")
vignette_artifact("workflow_plot_object_components.csv")
vignette_artifact("workflow_visual_checklist.csv")

## ----fit-full-----------------------------------------------------------------
# fit <- fit_mfrm(
#   data = ej2021_study1,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   quad_points = 7
# )
# 
# diag <- diagnose_mfrm(
#   fit,
#   residual_pca = "none",
#   diagnostic_mode = "both",
#   fit_df_method = "both"
# )
# 
# summary(fit, profile = "fit", detail = "brief")
# summary(diag)
# 
# # Keep the final figure flow explicit: fit -> Wright map -> follow-up plots.
# s <- summary(fit, profile = "facets", diagnostics = diag)
# res <- s$results
# s
# plot(res, type = "wright", preset = "publication", show_ci = TRUE, top_n = Inf)
# 
# # Optional closest FACETS-style asterisk ruler (without mfrmr CI overlays).
# plot(res, type = "wright", renderer = "facets",
#      category_labels = rubric_labels, show_ci = FALSE,
#      preset = "publication")
# 
# plot(
#   res,
#   type = "fit_pathway",
#   fit_stat = "Infit",
#   include_person = TRUE,
#   top_n_person = 12,
#   person_labels = "none",
#   facet_labels = "flagged",
#   preset = "publication"
# )

## ----fit-full-pca-------------------------------------------------------------
# diag_pca <- diagnose_mfrm(
#   fit,
#   residual_pca = "both",
#   pca_max_factors = 6
# )
# 
# summary(diag_pca)

## ----strict-rsm-pcm-----------------------------------------------------------
# fit_rsm_strict <- fit_mfrm(
#   data = toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   quad_points = 7,
#   maxit = 30
# )
# diag_rsm_strict <- diagnose_mfrm(
#   fit_rsm_strict,
#   diagnostic_mode = "both",
#   residual_pca = "none"
# )
# 
# fit_pcm_strict <- fit_mfrm(
#   data = toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "PCM",
#   step_facet = "Criterion",
#   quad_points = 7,
#   maxit = 30
# )
# diag_pcm_strict <- diagnose_mfrm(
#   fit_pcm_strict,
#   diagnostic_mode = "both",
#   residual_pca = "none"
# )
# 
# summary(diag_rsm_strict)$diagnostic_basis[, c("DiagnosticPath", "Status", "Basis")]
# summary(diag_pcm_strict)$diagnostic_basis[, c("DiagnosticPath", "Status", "Basis")]

## ----strict-screening---------------------------------------------------------
# screen_rsm <- evaluate_mfrm_diagnostic_screening(
#   design = list(person = 18, rater = 3, criterion = 3, assignment = 3),
#   reps = 1,
#   scenarios = c("well_specified", "local_dependence"),
#   model = "RSM",
#   maxit = 30,
#   quad_points = 7,
#   seed = 123
# )
# screen_pcm <- evaluate_mfrm_diagnostic_screening(
#   design = list(person = 18, rater = 3, criterion = 3, assignment = 3),
#   reps = 1,
#   scenarios = c("well_specified", "step_structure_misspecification"),
#   model = "PCM",
#   maxit = 30,
#   quad_points = 7,
#   seed = 123
# )
# 
# screen_rsm$performance_summary[, c("Scenario", "EvaluationUse", "LegacyAnyFlagRate", "StrictAnyFlagRate")]
# screen_pcm$performance_summary[, c("Scenario", "EvaluationUse", "LegacySensitivityProxy", "StrictSensitivityProxy", "DeltaStrictMinusLegacyFlagRate")]

## ----strict-reporting-route---------------------------------------------------
# chk_rsm_strict <- reporting_checklist(fit_rsm_strict, diagnostics = diag_rsm_strict)
# subset(
#   chk_rsm_strict$checklist,
#   Section == "Visual Displays" &
#     Item %in% c("QC / facet dashboard", "Strict marginal visuals", "Precision / information curves"),
#   c("Item", "Available", "DraftReady", "NextAction")
# )

## ----residual-pca-------------------------------------------------------------
# pca <- analyze_residual_pca(diag_pca, mode = "both")
# plot_residual_pca(pca, mode = "overall", plot_type = "scree")

## ----bias-apa-----------------------------------------------------------------
# data("mfrmr_example_bias", package = "mfrmr")
# bias_df <- mfrmr_example_bias
# fit_bias <- fit_mfrm(
#   bias_df,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   quad_points = 7
# )
# diag_bias <- diagnose_mfrm(fit_bias, residual_pca = "none")
# bias <- estimate_bias(fit_bias, diag_bias, facet_a = "Rater", facet_b = "Criterion")
# fixed <- build_fixed_reports(bias)
# apa <- build_apa_outputs(fit_bias, diag_bias, bias_results = bias)
# 
# mfrm_threshold_profiles()
# vis <- build_visual_summaries(fit_bias, diag_bias, threshold_profile = "standard")
# vis$warning_map$residual_pca_overall

## ----reporting-api------------------------------------------------------------
# spec <- specifications_report(fit, title = "Study run")
# data_qc <- data_quality_report(
#   fit,
#   data = ej2021_study1,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score"
# )
# iter <- estimation_iteration_report(fit, max_iter = 8)
# subset_rep <- subset_connectivity_report(fit, diagnostics = diag)
# facet_stats <- facet_statistics_report(fit, diagnostics = diag)
# cat_structure <- category_structure_report(fit, diagnostics = diag)
# cat_curves <- category_curves_report(fit, theta_points = 101)
# bias_rep <- bias_interaction_report(bias, top_n = 20)
# plot_bias_interaction(bias_rep, plot = "scatter")

## ----gd-study-reading-order, eval=FALSE---------------------------------------
# if (requireNamespace("lme4", quietly = TRUE)) {
#   gt <- mfrm_generalizability(fit)
#   gt$coefficients[, c("G", "Phi", "GStatus", "PhiStatus",
#                       "IdentificationStatus")]
# 
#   ds <- mfrm_d_study(
#     gt,
#     data.frame(Rater = c(2, 3, 4), Criterion = 4),
#     residual_scaling = "sensitivity"
#   )
#   ds[, c("n_Rater", "n_Criterion", "ResidualScaling",
#          "G", "Phi", "GStatus", "PhiStatus", "IdentificationStatus")]
# }

## ----design-prediction--------------------------------------------------------
# sim_spec <- build_mfrm_sim_spec(
#   n_person = 30,
#   n_rater = 4,
#   n_criterion = 4,
#   raters_per_person = 2,
#   assignment = "rotating"
# )
# 
# recovery <- suppressWarnings(
#   evaluate_mfrm_recovery(
#     sim_spec = sim_spec,
#     reps = 2,
#     maxit = 30,
#     include_diagnostics = TRUE,
#     diagnostic_fit_df_method = "both",
#     seed = 2
#   )
# )
# 
# summary(recovery)$recovery_summary[, c("ParameterType", "Facet", "RMSE", "Bias")]
# plot(recovery, type = "summary", metric = "rmse", draw = FALSE)$data$plot_table
# 
# recovery_review <- assess_mfrm_recovery(
#   recovery,
#   min_reps = 2,
#   min_se_available = NULL,
#   max_mcse_rmse_ratio = NULL,
#   max_rmse = c(facet = 1, step = 1, default = 1.5),
#   max_abs_bias = c(default = 0.75)
# )
# 
# summary(recovery_review)$checklist[, c("Section", "Item", "Status")]
# summary(recovery_review)$reading_order
# summary(recovery_review)$condition_reporting_notes
# summary(recovery_review)$condition_review
# summary(recovery_review)$diagnostic_reporting_notes
# summary(recovery_review)$diagnostic_review
# 
# status_plot <- plot(recovery_review, type = "status", draw = FALSE)
# status_plot$data$section_status
# status_plot$data$reading_order
# 
# metric_plot <- plot(recovery_review, type = "metrics", metric = "rmse", draw = FALSE)
# metric_plot$data$plot_table
# metric_plot$data$guidance
# 
# recovery_bundle <- build_summary_table_bundle(
#   recovery_review,
#   appendix_preset = "recommended"
# )
# recovery_bundle$table_index[, c("Table", "Rows", "Role")]
# 
# pred_pop <- predict_mfrm_population(
#   sim_spec = sim_spec,
#   reps = 2,
#   maxit = 30,
#   seed = 1
# )
# 
# summary(pred_pop)$forecast[, c("Facet", "MeanSeparation", "McseSeparation")]
# 
# keep_people <- unique(toy$Person)[1:18]
# toy_mml <- suppressWarnings(
#   fit_mfrm(
#     toy[toy$Person %in% keep_people, , drop = FALSE],
#     person = "Person",
#     facets = c("Rater", "Criterion"),
#     score = "Score",
#     method = "MML",
#     quad_points = 5,
#     maxit = 30
#   )
# )
# 
# new_units <- data.frame(
#   Person = c("NEW01", "NEW01"),
#   Rater = unique(toy$Rater)[1],
#   Criterion = unique(toy$Criterion)[1:2],
#   Score = c(2, 3)
# )
# 
# pred_units <- predict_mfrm_units(toy_mml, new_units, n_draws = 0)
# pv_units <- sample_mfrm_plausible_values(toy_mml, new_units, n_draws = 2, seed = 1)
# 
# summary(pred_units)$estimates[, c("Person", "Estimate", "Lower", "Upper")]
# summary(pv_units)$draw_summary[, c("Person", "Draws", "MeanValue")]

## ----recovery-export, eval = FALSE--------------------------------------------
# export_summary_appendix(
#   list(recovery = recovery, recovery_review = recovery_review),
#   output_dir = tempdir(),
#   prefix = "mfrmr_recovery_appendix",
#   preset = "recommended",
#   include_html = FALSE,
#   overwrite = TRUE
# )

