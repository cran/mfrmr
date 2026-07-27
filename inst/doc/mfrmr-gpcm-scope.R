## ----include = FALSE----------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = !is_cran_check
)

## ----capability-supported-----------------------------------------------------
# library(mfrmr)
# gpcm_capability_matrix("supported")[, c("Area", "Status")]

## ----capability-with-caveat---------------------------------------------------
# gpcm_capability_matrix("supported_with_caveat")[, c("Area", "Status")]

## ----capability-blocked-------------------------------------------------------
# gpcm_capability_matrix("blocked")[, c("Area", "Status", "RecommendedRoute")]

## ----capability-deferred------------------------------------------------------
# gpcm_capability_matrix("deferred")[, c("Area", "Status", "Boundary", "RecommendedRoute")]

