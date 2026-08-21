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
# bias_df <- load_mfrmr_data("example_bias")
# 
# # The vignette uses compact quadrature so optional local execution stays fast.
# # For final DFF or linking evidence, refit with the package default or a higher
# # quadrature setting and record that setting in the analysis log.
# fit <- fit_mfrm(
#   bias_df,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   quad_points = 31
# )
# 
# diag <- diagnose_mfrm(fit, residual_pca = "none")

## ----connectivity-------------------------------------------------------------
# sc <- subset_connectivity_report(fit, diagnostics = diag)
# 
# sc$summary[, c("Subset", "Observations", "ObservationPercent")]
# plot(sc, type = "design_matrix", preset = "publication")

## ----anchors------------------------------------------------------------------
# anchors <- make_anchor_table(fit, facets = "Criterion")
# head(anchors)

## ----dff-residual-------------------------------------------------------------
# dff_resid <- analyze_dff(
#   fit,
#   diag,
#   facet = "Criterion",
#   group = "Group",
#   data = bias_df,
#   method = "residual"
# )
# 
# dff_resid$summary
# head(
#   dff_resid$dif_table[, c("Level", "Group1", "Group2", "Classification", "ClassificationSystem")],
#   8
# )
# plot_dif_heatmap(dff_resid)

## ----dff-refit----------------------------------------------------------------
# dff_refit <- analyze_dff(
#   fit,
#   diag,
#   facet = "Criterion",
#   group = "Group",
#   data = bias_df,
#   method = "refit"
# )
# 
# dff_refit$summary
# head(
#   dff_refit$dif_table[, c("Level", "Group1", "Group2", "Classification", "ContrastComparable")],
#   8
# )

## ----dff-follow-up------------------------------------------------------------
# dit <- dif_interaction_table(
#   fit,
#   diag,
#   facet = "Criterion",
#   group = "Group",
#   data = bias_df
# )
# 
# head(dit$table)
# 
# dr <- dif_report(dff_resid)
# cat(dr$narrative)

## ----model-estimated-interaction, eval = FALSE--------------------------------
# fit_add <- fit_mfrm(
#   bias_df,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   quad_points = 7
# )
# 
# fit_interaction <- fit_mfrm(
#   bias_df,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   facet_interactions = "Rater:Criterion",
#   quad_points = 31
# )
# 
# interaction_effect_table(fit_interaction)
# compare_mfrm(Additive = fit_add, Interaction = fit_interaction, nested = TRUE)

## ----drift-route, eval = FALSE------------------------------------------------
# declared_common_facets <- c("Criterion")
# 
# fit1 <- fit_mfrm(
#   wave1, "Person", c("Rater", "Criterion"), "Score",
#   method = "MML"
# )
# 
# anchored <- anchor_to_baseline(
#   wave2,
#   fit1,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   anchor_facets = declared_common_facets
# )
# 
# fit2 <- fit_mfrm(
#   wave2, "Person", c("Rater", "Criterion"), "Score",
#   method = "MML"
# )
# drift <- detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2))
# plot_anchor_drift(drift, type = "drift", preset = "publication")

