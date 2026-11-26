## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5, fig.height = 4.2)

## -----------------------------------------------------------------------------
library(mfrmr)
ratings <- load_mfrmr_data("example_core")
wf <- mfrm_analysis(ratings, "Person", c("Rater", "Criterion"), "Score")
glance(wf)[c("n_obs", "connected", "converged", "inference_status")]
plot(wf, type = "wright")
plan <- data.frame(Rater = c("R01", "R02"),
                   Criterion = c("Content", "Organization"))
report_file <- render_mfrm_report(
  wf, file = tempfile(fileext = ".html"), format = "html", reference_plan = plan
)
feedback <- render_mfrm_report(wf, rater = "R01")

## -----------------------------------------------------------------------------
head(ratings, 4)

conn <- check_mfrm_connectivity(
  ratings,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score"
)
conn$connected
conn$n_subsets

## -----------------------------------------------------------------------------
fit_rsm <- fit_mfrm(
  ratings,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  model = "RSM",
  method = "MML",
  person_score = "WLE"
)
print(fit_rsm)
fit_rsm$facets[, c("Facet", "Level", "Estimate", "ModelSE", "RealSE", "FairM", "InfitMnSq")]

fit_pcm <- fit_mfrm(
  ratings,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  model = "PCM",
  step_facet = "Criterion",
  method = "MML"
)

fit_gpcm <- fit_mfrm(
  ratings,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  model = "GPCM",
  step_facet = "Criterion",
  slope_facet = "Criterion",
  spec = list(gpcm_mml_identification = "fixed_standard_normal"),
  control = list(quad_points = 31L)
)
fit_gpcm$slopes[, c("Facet", "Level", "Estimate", "ModelSE")]
comparison <- anova(fit_rsm, fit_pcm, fit_gpcm)
comparison$lrt

## -----------------------------------------------------------------------------
predict(
  fit_rsm,
  newdata = data.frame(theta = c(-1, 0, 1), Rater = "R01", Criterion = "Content"),
  type = "response"
)
confint(fit_rsm, component = "facets")

plot(fit_rsm, type = "wright", palette = "okabe_ito")
plot(fit_gpcm, type = "slopes", palette = "monochrome")

## -----------------------------------------------------------------------------
set.seed(42)
rep_scores <- replicate(5L, sample(1:4, nrow(ratings), replace = TRUE))
batch_summary <- mfrm_batch_fit(
  ratings,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  score_matrix = rep_scores,
  method = "CML",
  model = "RSM",
  rating_min = 1,
  rating_max = 4,
  threads = 1L,
  return_fits = FALSE
)
batch_summary

