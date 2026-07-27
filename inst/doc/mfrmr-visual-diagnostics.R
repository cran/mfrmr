## ----include = FALSE----------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = !is_cran_check
)

## ----setup--------------------------------------------------------------------
# library(mfrmr)
# 
# toy <- load_mfrmr_data("example_operational")
# 
# fit <- fit_mfrm(
#   toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM"
# )
# 
# diag <- diagnose_mfrm(fit, residual_pca = "none")
# checklist <- reporting_checklist(fit, diagnostics = diag)
# subset(
#   checklist$checklist,
#   Section == "Visual Displays",
#   c("Item", "Available", "NextAction")
# )

## ----wright-------------------------------------------------------------------
# plot(fit, type = "wright", preset = "publication", show_ci = TRUE)

## ----wright-facets-style------------------------------------------------------
# plot(
#   fit,
#   type = "wright",
#   renderer = "facets",
#   category_labels = c(
#     `1` = "Beginning", `2` = "Developing", `3` = "Secure", `4` = "Advanced"
#   )
# )

## ----pathway------------------------------------------------------------------
# plot(fit, type = "pathway", preset = "publication")

## ----fit-pathway--------------------------------------------------------------
# plot(
#   fit,
#   type = "fit_pathway",
#   diagnostics = diag,
#   fit_stat = "Infit",
#   fit_scale = "mnsq",
#   include_person = TRUE,
#   show_ci = TRUE,
#   preset = "publication"
# )

## ----unexpected---------------------------------------------------------------
# plot_unexpected(
#   fit,
#   diagnostics = diag,
#   abs_z_min = 1.5,
#   prob_max = 0.4,
#   plot_type = "scatter",
#   preset = "publication"
# )

## ----displacement-------------------------------------------------------------
# plot_displacement(
#   fit,
#   diagnostics = diag,
#   anchored_only = FALSE,
#   plot_type = "lollipop",
#   preset = "publication"
# )

## ----strict-marginal----------------------------------------------------------
# fit_strict <- fit_mfrm(
#   toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   quad_points = 7,
#   maxit = 40
# )
# 
# diag_strict <- diagnose_mfrm(
#   fit_strict,
#   residual_pca = "none",
#   diagnostic_mode = "both"
# )
# 
# strict_checklist <- reporting_checklist(fit_strict, diagnostics = diag_strict)
# subset(
#   strict_checklist$checklist,
#   Section == "Visual Displays" &
#     Item %in% c("QC / facet dashboard", "Strict marginal visuals"),
#   c("Item", "Available", "NextAction")
# )
# 
# plot_marginal_fit(
#   diag_strict,
#   top_n = 12,
#   preset = "publication"
# )

## ----linking------------------------------------------------------------------
# sc <- subset_connectivity_report(fit, diagnostics = diag)
# plot(sc, type = "design_matrix", preset = "publication")

## ----eval = FALSE-------------------------------------------------------------
# drift <- detect_anchor_drift(current_fit, baseline = baseline_anchors)
# plot_anchor_drift(drift, type = "heatmap", preset = "publication")

## ----residual-pca-------------------------------------------------------------
# diag_pca <- diagnose_mfrm(fit, residual_pca = "both", pca_max_factors = 4)
# pca <- analyze_residual_pca(diag_pca, mode = "both")
# plot_residual_pca(pca, mode = "overall", plot_type = "scree", preset = "publication")

## ----bias---------------------------------------------------------------------
# bias_df <- load_mfrmr_data("example_bias")
# 
# fit_bias <- fit_mfrm(
#   bias_df,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   quad_points = 7
# )
# 
# diag_bias <- diagnose_mfrm(fit_bias, residual_pca = "none")
# bias <- estimate_bias(fit_bias, diag_bias, facet_a = "Rater", facet_b = "Criterion")
# 
# plot_bias_interaction(
#   bias,
#   plot = "facet_profile",
#   preset = "publication"
# )

## ----custom-plot-data---------------------------------------------------------
# plot(fit, type = "wright", preset = "monochrome")
# 
# wright_payload <- plot(fit, type = "wright", draw = FALSE, preset = "publication")
# plot_data_components(wright_payload)
# 
# locations <- plot_data(wright_payload, component = "locations")
# head(locations)
# 
# pathway_long <- plot_data(
#   fit,
#   type = "pathway",
#   component = "pathway_long",
#   preset = "publication"
# )
# head(pathway_long[, c("Layer", "CurveGroup", "Theta", "Value")])

## ----custom-guidance----------------------------------------------------------
# names(wright_payload$data)
# wright_payload$data$reference_lines

## ----response-time-review-----------------------------------------------------
# toy_rt <- toy
# toy_rt$ResponseTime <- 12 + (seq_len(nrow(toy_rt)) %% 7) +
#   as.numeric(toy_rt$Score)
# toy_rt$ResponseTime[1] <- 2
# toy_rt$ResponseTime[2] <- 38
# 
# rt <- response_time_review(
#   toy_rt,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   time = "ResponseTime",
#   rapid_quantile = 0.10,
#   slow_quantile = 0.90
# )
# 
# summary(rt)
# plot_response_time_review(rt, type = "distribution", preset = "publication")
# plot_response_time_review(rt, type = "person", preset = "publication")

## ----shrinkage-funnel---------------------------------------------------------
# fit_eb <- apply_empirical_bayes_shrinkage(fit)
# 
# shrink <- plot_shrinkage_funnel(
#   fit_eb,
#   show_ci = TRUE,
#   ci_level = 0.95,
#   preset = "publication",
#   draw = FALSE
# )
# 
# head(shrink$data$table[, c(
#   "Facet", "Level", "RawEstimate", "RawCI_Lower", "RawCI_Upper",
#   "ShrunkEstimate", "ShrunkCI_Lower", "ShrunkCI_Upper",
#   "ShrinkageFactor"
# )])
# 
# plot_shrinkage_funnel(
#   fit_eb,
#   show_ci = TRUE,
#   ci_level = 0.95,
#   preset = "publication"
# )

