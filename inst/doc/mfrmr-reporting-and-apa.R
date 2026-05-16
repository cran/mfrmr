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
# toy <- load_mfrmr_data("example_core")
# 
# # The vignette uses compact quadrature so optional local execution stays fast.
# # For final manuscript reporting, refit with the package default or a higher
# # quadrature setting and record that setting in the analysis log.
# fit <- fit_mfrm(
#   toy,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   method = "MML",
#   model = "RSM",
#   quad_points = 7
# )
# 
# diag <- diagnose_mfrm(fit, residual_pca = "none")

## ----checklist----------------------------------------------------------------
# chk <- reporting_checklist(fit, diagnostics = diag)
# 
# head(
#   chk$checklist[, c("Section", "Item", "DraftReady", "Priority", "NextAction")],
#   10
# )

## ----precision----------------------------------------------------------------
# prec <- precision_review_report(fit, diagnostics = diag)
# 
# prec$profile
# prec$checks

## ----apa----------------------------------------------------------------------
# apa <- build_apa_outputs(
#   fit,
#   diagnostics = diag,
#   context = list(
#     assessment = "Writing assessment",
#     setting = "Local scoring study",
#     scale_desc = "0-4 rubric scale",
#     rater_facet = "Rater"
#   )
# )
# 
# cat(apa$report_text)

## ----section-map--------------------------------------------------------------
# apa$section_map[, c("SectionId", "Heading", "Available")]

## ----apa-tables---------------------------------------------------------------
# tbl_summary <- apa_table(fit, which = "summary")
# tbl_reliability <- apa_table(fit, which = "reliability", diagnostics = diag)
# 
# tbl_summary$caption
# tbl_reliability$note

## ----visuals------------------------------------------------------------------
# vis <- build_visual_summaries(
#   fit,
#   diagnostics = diag,
#   threshold_profile = "standard"
# )
# 
# names(vis)
# names(vis$warning_map)

## ----bias-screen--------------------------------------------------------------
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
# apa_bias <- build_apa_outputs(fit_bias, diagnostics = diag_bias, bias_results = bias)
# 
# apa_bias$section_map[, c("SectionId", "Available", "Heading")]

