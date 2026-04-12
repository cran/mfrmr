## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load-data----------------------------------------------------------------
library(mfrmr)

list_mfrmr_data()

data("ej2021_study1", package = "mfrmr")
head(ej2021_study1)

study1_alt <- load_mfrmr_data("study1")
identical(names(ej2021_study1), names(study1_alt))

## ----toy-setup----------------------------------------------------------------
data("mfrmr_example_core", package = "mfrmr")
toy <- mfrmr_example_core

fit_toy <- fit_mfrm(
  data = toy,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "JML",
  model = "RSM",
  maxit = 15
)
diag_toy <- diagnose_mfrm(fit_toy, residual_pca = "none")

summary(fit_toy)$overview
summary(diag_toy)$overview
names(plot(fit_toy, draw = FALSE))

## ----diagnostics-reporting----------------------------------------------------
t4_toy <- unexpected_response_table(
  fit_toy,
  diagnostics = diag_toy,
  abs_z_min = 1.5,
  prob_max = 0.4,
  top_n = 10
)
t12_toy <- fair_average_table(fit_toy, diagnostics = diag_toy)
t13_toy <- bias_interaction_report(
  estimate_bias(fit_toy, diag_toy,
                facet_a = "Rater", facet_b = "Criterion",
                max_iter = 2),
  top_n = 10
)

class(summary(t4_toy))
class(summary(t12_toy))
class(summary(t13_toy))

names(plot(t4_toy, draw = FALSE))
names(plot(t12_toy, draw = FALSE))
names(plot(t13_toy, draw = FALSE))

chk_toy <- reporting_checklist(fit_toy, diagnostics = diag_toy)
subset(
  chk_toy$checklist,
  Section == "Visual Displays",
  c("Item", "DraftReady", "NextAction")
)

## ----fit-full-----------------------------------------------------------------
fit <- fit_mfrm(
  data = ej2021_study1,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",
  model = "RSM",
  quad_points = 7
)

diag <- diagnose_mfrm(
  fit,
  residual_pca = "none"
)

summary(fit)
summary(diag)

## ----fit-full-pca-------------------------------------------------------------
diag_pca <- diagnose_mfrm(
  fit,
  residual_pca = "both",
  pca_max_factors = 6
)

summary(diag_pca)

## ----strict-rsm-pcm-----------------------------------------------------------
fit_rsm_strict <- fit_mfrm(
  data = toy,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",
  model = "RSM",
  quad_points = 7,
  maxit = 10
)
diag_rsm_strict <- diagnose_mfrm(
  fit_rsm_strict,
  diagnostic_mode = "both",
  residual_pca = "none"
)

fit_pcm_strict <- fit_mfrm(
  data = toy,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",
  model = "PCM",
  step_facet = "Criterion",
  quad_points = 7,
  maxit = 10
)
diag_pcm_strict <- diagnose_mfrm(
  fit_pcm_strict,
  diagnostic_mode = "both",
  residual_pca = "none"
)

summary(diag_rsm_strict)$diagnostic_basis[, c("DiagnosticPath", "Status", "Basis")]
summary(diag_pcm_strict)$diagnostic_basis[, c("DiagnosticPath", "Status", "Basis")]

## ----strict-screening---------------------------------------------------------
screen_rsm <- evaluate_mfrm_diagnostic_screening(
  design = list(person = 18, rater = 3, criterion = 3, assignment = 3),
  reps = 1,
  scenarios = c("well_specified", "local_dependence"),
  model = "RSM",
  maxit = 8,
  quad_points = 7,
  seed = 123
)
screen_pcm <- evaluate_mfrm_diagnostic_screening(
  design = list(person = 18, rater = 3, criterion = 3, assignment = 3),
  reps = 1,
  scenarios = c("well_specified", "step_structure_misspecification"),
  model = "PCM",
  maxit = 8,
  quad_points = 7,
  seed = 123
)

screen_rsm$performance_summary[, c("Scenario", "EvaluationUse", "LegacyAnyFlagRate", "StrictAnyFlagRate")]
screen_pcm$performance_summary[, c("Scenario", "EvaluationUse", "LegacySensitivityProxy", "StrictSensitivityProxy", "DeltaStrictMinusLegacyFlagRate")]

## ----strict-reporting-route---------------------------------------------------
chk_rsm_strict <- reporting_checklist(fit_rsm_strict, diagnostics = diag_rsm_strict)
subset(
  chk_rsm_strict$checklist,
  Section == "Visual Displays" &
    Item %in% c("QC / facet dashboard", "Strict marginal visuals", "Precision / information curves"),
  c("Item", "Available", "DraftReady", "NextAction")
)

## ----residual-pca-------------------------------------------------------------
pca <- analyze_residual_pca(diag_pca, mode = "both")
plot_residual_pca(pca, mode = "overall", plot_type = "scree")

## ----bias-apa-----------------------------------------------------------------
data("mfrmr_example_bias", package = "mfrmr")
bias_df <- mfrmr_example_bias
fit_bias <- fit_mfrm(
  bias_df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",
  model = "RSM",
  quad_points = 7
)
diag_bias <- diagnose_mfrm(fit_bias, residual_pca = "none")
bias <- estimate_bias(fit_bias, diag_bias, facet_a = "Rater", facet_b = "Criterion")
fixed <- build_fixed_reports(bias)
apa <- build_apa_outputs(fit_bias, diag_bias, bias_results = bias)

mfrm_threshold_profiles()
vis <- build_visual_summaries(fit_bias, diag_bias, threshold_profile = "standard")
vis$warning_map$residual_pca_overall

## ----reporting-api------------------------------------------------------------
spec <- specifications_report(fit, title = "Study run")
data_qc <- data_quality_report(
  fit,
  data = ej2021_study1,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score"
)
iter <- estimation_iteration_report(fit, max_iter = 8)
subset_rep <- subset_connectivity_report(fit, diagnostics = diag)
facet_stats <- facet_statistics_report(fit, diagnostics = diag)
cat_structure <- category_structure_report(fit, diagnostics = diag)
cat_curves <- category_curves_report(fit, theta_points = 101)
bias_rep <- bias_interaction_report(bias, top_n = 20)
plot_bias_interaction(bias_rep, plot = "scatter")

## ----design-prediction--------------------------------------------------------
sim_spec <- build_mfrm_sim_spec(
  n_person = 30,
  n_rater = 4,
  n_criterion = 4,
  raters_per_person = 2,
  assignment = "rotating"
)

pred_pop <- predict_mfrm_population(
  sim_spec = sim_spec,
  reps = 2,
  maxit = 10,
  seed = 1
)

summary(pred_pop)$forecast[, c("Facet", "MeanSeparation", "McseSeparation")]

keep_people <- unique(toy$Person)[1:18]
toy_mml <- suppressWarnings(
  fit_mfrm(
    toy[toy$Person %in% keep_people, , drop = FALSE],
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    quad_points = 5,
    maxit = 15
  )
)

new_units <- data.frame(
  Person = c("NEW01", "NEW01"),
  Rater = unique(toy$Rater)[1],
  Criterion = unique(toy$Criterion)[1:2],
  Score = c(2, 3)
)

pred_units <- predict_mfrm_units(toy_mml, new_units, n_draws = 0)
pv_units <- sample_mfrm_plausible_values(toy_mml, new_units, n_draws = 2, seed = 1)

summary(pred_units)$estimates[, c("Person", "Estimate", "Lower", "Upper")]
summary(pv_units)$draw_summary[, c("Person", "Draws", "MeanValue")]

