## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5, fig.height = 4.2)

## -----------------------------------------------------------------------------
library(mfrmr)

ratings <- load_mfrmr_data("example_core")

gt <- mfrm_g_study(
  ratings,
  scores = "Score",
  person = "Person",
  facets = c("Rater", "Criterion"),
  method = "urgenova",
  d_study_grid = expand.grid(Rater = c(1, 2, 3, 4), Criterion = c(2, 3, 4))
)
print(gt)
confint(gt)
plans <- mfrm_d_study(gt)
plans
equal_budget <- plans[plans$TotalRatingsPerPerson == 8, ]
equal_budget[, c("Rater", "Criterion", "TotalRatingsPerPerson", "G", "Phi", "AbsSEM")]
plot(gt, type = "dstudy", palette = "okabe_ito")

## -----------------------------------------------------------------------------
nested_ratings <- load_mfrmr_data("example_gtheory_nested")
head(nested_ratings, 4)

gt_nested <- mfrm_g_study(
  nested_ratings,
  scores = c("Content", "Organization"),
  weights = c(0.6, 0.4),
  person = "Person",
  facets = c("Rater", "Task"),
  nesting = c(Rater = "Task"),
  method = "urgenova"
)
print(gt_nested)
gt_nested$univariate
confint(gt_nested, component = "coefficients")

