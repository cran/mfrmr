## ----include = FALSE----------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = !is_cran_check
)

## ----adaptive-order-review----------------------------------------------------
# library(mfrmr)
# toy <- load_mfrmr_data("example_operational")
# fit_adaptive <- fit_mfrm(
#   data = toy, person = "Person", facets = c("Rater", "Criterion"),
#   score = "Score", method = "MML", model = "RSM",
#   mml_integration = "adaptive", quad_points = 31
# )
# quadrature_review <- mml_quadrature_sensitivity(
#   fit_adaptive, toy, quad_points = c(31, 61)
# )
# summary(quadrature_review)

