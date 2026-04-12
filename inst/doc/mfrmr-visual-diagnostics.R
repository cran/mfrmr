## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(mfrmr)

toy <- load_mfrmr_data("example_core")

fit <- fit_mfrm(
  toy,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "JML",
  model = "RSM",
  maxit = 20
)

diag <- diagnose_mfrm(fit, residual_pca = "none")
checklist <- reporting_checklist(fit, diagnostics = diag)
subset(
  checklist$checklist,
  Section == "Visual Displays",
  c("Item", "Available", "NextAction")
)

## ----wright-------------------------------------------------------------------
plot(fit, type = "wright", preset = "publication", show_ci = TRUE)

## ----pathway------------------------------------------------------------------
plot(fit, type = "pathway", preset = "publication")

## ----unexpected---------------------------------------------------------------
plot_unexpected(
  fit,
  diagnostics = diag,
  abs_z_min = 1.5,
  prob_max = 0.4,
  plot_type = "scatter",
  preset = "publication"
)

## ----displacement-------------------------------------------------------------
plot_displacement(
  fit,
  diagnostics = diag,
  anchored_only = FALSE,
  plot_type = "lollipop",
  preset = "publication"
)

## ----strict-marginal----------------------------------------------------------
fit_strict <- fit_mfrm(
  toy,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",
  model = "RSM",
  quad_points = 7,
  maxit = 40
)

diag_strict <- diagnose_mfrm(
  fit_strict,
  residual_pca = "none",
  diagnostic_mode = "both"
)

strict_checklist <- reporting_checklist(fit_strict, diagnostics = diag_strict)
subset(
  strict_checklist$checklist,
  Section == "Visual Displays" &
    Item %in% c("QC / facet dashboard", "Strict marginal visuals"),
  c("Item", "Available", "NextAction")
)

plot_marginal_fit(
  diag_strict,
  top_n = 12,
  preset = "publication"
)

## ----linking------------------------------------------------------------------
sc <- subset_connectivity_report(fit, diagnostics = diag)
plot(sc, type = "design_matrix", preset = "publication")

## ----eval = FALSE-------------------------------------------------------------
# drift <- detect_anchor_drift(current_fit, baseline = baseline_anchors)
# plot_anchor_drift(drift, type = "heatmap", preset = "publication")

## ----residual-pca-------------------------------------------------------------
diag_pca <- diagnose_mfrm(fit, residual_pca = "both", pca_max_factors = 4)
pca <- analyze_residual_pca(diag_pca, mode = "both")
plot_residual_pca(pca, mode = "overall", plot_type = "scree", preset = "publication")

## ----bias---------------------------------------------------------------------
bias_df <- load_mfrmr_data("example_bias")

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

plot_bias_interaction(
  bias,
  plot = "facet_profile",
  preset = "publication"
)

