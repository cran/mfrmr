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
# data("ej2021_study1", package = "mfrmr")
# 
# run <- run_mfrm_facets(
#   data = ej2021_study1,
#   person = "Person",
#   facets = c("Rater", "Criterion"),
#   score = "Score",
#   model = "RSM",
#   method = "JML"
# )
# 
# names(run)

## ----facets-mode-summary------------------------------------------------------
# summary(run$fit)
# head(run$fair_average)

