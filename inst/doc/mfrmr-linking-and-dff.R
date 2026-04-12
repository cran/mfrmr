## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(mfrmr)

bias_df <- load_mfrmr_data("example_bias")

fit <- fit_mfrm(
  bias_df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",
  model = "RSM",
  quad_points = 7
)

diag <- diagnose_mfrm(fit, residual_pca = "none")

## ----connectivity-------------------------------------------------------------
sc <- subset_connectivity_report(fit, diagnostics = diag)

sc$summary[, c("Subset", "Observations", "ObservationPercent")]
plot(sc, type = "design_matrix", preset = "publication")

## ----anchors------------------------------------------------------------------
anchors <- make_anchor_table(fit, facets = "Criterion")
head(anchors)

## ----dff-residual-------------------------------------------------------------
dff_resid <- analyze_dff(
  fit,
  diag,
  facet = "Criterion",
  group = "Group",
  data = bias_df,
  method = "residual"
)

dff_resid$summary
head(
  dff_resid$dif_table[, c("Level", "Group1", "Group2", "Classification", "ClassificationSystem")],
  8
)
plot_dif_heatmap(dff_resid)

## ----dff-refit----------------------------------------------------------------
dff_refit <- analyze_dff(
  fit,
  diag,
  facet = "Criterion",
  group = "Group",
  data = bias_df,
  method = "refit"
)

dff_refit$summary
head(
  dff_refit$dif_table[, c("Level", "Group1", "Group2", "Classification", "ContrastComparable")],
  8
)

## ----dff-follow-up------------------------------------------------------------
dit <- dif_interaction_table(
  fit,
  diag,
  facet = "Criterion",
  group = "Group",
  data = bias_df
)

head(dit$table)

dr <- dif_report(dff_resid)
cat(dr$narrative)

## ----drift-route, eval = FALSE------------------------------------------------
# d1 <- load_mfrmr_data("study1")
# d2 <- load_mfrmr_data("study2")
# 
# fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score",
#                  method = "JML", maxit = 25)
# fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score",
#                  method = "JML", maxit = 25)
# 
# anchored <- anchor_to_baseline(
#   d2,
#   fit1,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score"
# )
# 
# drift <- detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2))
# plot_anchor_drift(drift, type = "drift", preset = "publication")

