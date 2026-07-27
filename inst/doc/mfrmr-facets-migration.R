## ----include = FALSE----------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = !is_cran_check
)

## ----facets-mode--------------------------------------------------------------
# library(mfrmr)
# data("mfrmr_example_operational", package = "mfrmr")
# 
# run <- run_mfrm_facets(
#   data = mfrmr_example_operational,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   model = "RSM",
#   method = "JML"
# )
# 
# names(run)

## ----facets-mode-summary------------------------------------------------------
# jml_status <- summary(run$fit, profile = "fit", detail = "brief")
# jml_status$overview[, c(
#   "Model", "Method", "Converged", "InferenceReady",
#   "ConvergenceSeverity"
# )]
# jml_status$readiness
# head(run$fair_average)

