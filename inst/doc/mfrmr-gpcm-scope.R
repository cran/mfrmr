## ----include = FALSE----------------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = !is_cran_check
)

## ----gpcm-interval-route-guide------------------------------------------------
# interval_routes <- mfrmr::mfrmr_interval_guide("gpcm")
# interval_routes[interval_routes$Route %in% c(
#   "Fixed-facet model/sandwich intervals", "GPCM relative-slope intervals",
#   "GMFRM component-slope intervals", "GPCM slope profile intervals",
#   "Corrected-JML local root uncertainty"
# ), c("Route", "PrimaryHelper", "GPCMStatus")]

## ----corrected-jml-portable---------------------------------------------------
# # draft <- extract_mfrm_calibration(fit, scoring_quad_points = 101)
# # calibration <- freeze_mfrm_calibration(validate_mfrm_calibration(draft))
# # save_mfrm_calibration(calibration, "rubric-calibration.rds")
# # calibration <- load_mfrm_calibration("rubric-calibration.rds")
# # scores <- score_mfrm_calibration(calibration, new_ratings)
# # summary(scores)
# # plot(scores, type = "interval")
# # plot(scores, type = "precision")
# # as_ggplot(scores)

## ----two-family-workflow, eval=FALSE------------------------------------------
# fit_joint <- fit_mfrm(
#   ratings, person = "Candidate", facets = c("Task", "Assessor"), score = "Rating",
#   model = "GPCM", method = "MML",
#   slope_facet = c("Task", "Assessor"), step_facet = "Assessor",
#   noncenter_facet = "Assessor", gpcm_mml_identification = "fixed_standard_normal",
#   mml_engine = "em", rating_min = 0, rating_max = 3,
#   category_policy = "preserve", quad_points = 31, maxit = 500,
#   em_score_tol = 1e-6, progress = TRUE
# )
# summary(fit_joint)
# context <- expand.grid(
#   Theta = seq(-3, 3, length.out = 61),
#   Task = unique(ratings$Task)[1], Assessor = unique(ratings$Assessor)[1]
# )
# curves <- mfrm_curve_intervals(fit_joint, context)
# plot(curves) # provisional fitted probabilities; no confidence intervals
# ci_joint <- confint(fit_joint, parm = "slopes") # experimental component-slope intervals
# plot(ci_joint)
# apa_table(ci_joint, which = "checks")
# # For a figure without headings or notes:
# plot(curves, title = NULL, subtitle = NULL, caption = NULL)
# res <- mfrm_results(fit_joint, intervals = list(slopes = ci_joint, curves = curves))
# mfrm_report(res)
# # export_mfrm_results(res, output_dir = "joint_results")

## ----joint-quadrature, eval=FALSE---------------------------------------------
# # fit_joint above used 31 points. Include the original point count.
# integration <- mml_quadrature_sensitivity(
#   fit_joint, ratings, quad_points = c(31, 61),
#   adaptive_quad_points = c(31, 61)
# )
# summary(integration)
# integration$quadrature_review # fixed-calibration integration checks per Person
# integration$slopes # estimates, raw diagnostic SEs, interval availability
# ci_finer <- integration$intervals$q61 # already computed; no repeated fit
# apa_table(ci_finer, which = "numerical_checks")
# plot(ci_finer) # retains missing intervals and cautions
# # integration$fits$q61 is the corresponding fitted model.

## ----two-family-stan-export, eval=FALSE---------------------------------------
# stan_bundle <- export_mfrm_model(
#   fit_joint, "gmfrm-stan", format = "stan", log_slope_prior_sd = 1
# )
# model_data <- readRDS(stan_bundle$files[["model-data.rds"]])
# model_data$priors
# model_data$slope_levels
# # Review README.txt and run.R, then run Rscript run.R from the bundle directory.

## ----two-family-new-person, eval=FALSE----------------------------------------
# # new_ratings has the original Person, facet and score columns.
# scores <- predict_mfrm_units(fit_joint, new_ratings, scoring_quad_points = 121)
# summary(scores)
# saveRDS(scores, "two-family-new-person-scores.rds")

## ----two-family-score-results, eval=FALSE-------------------------------------
# res <- mfrm_results(fit_joint, scores = scores)
# res$tables$person_scores
# res$tables$scoring_integration
# report <- mfrm_report(res)

## ----two-family-portable-scoring, eval=FALSE----------------------------------
# draft <- extract_mfrm_calibration(fit_joint, scoring_quad_points = 121)
# calibration <- freeze_mfrm_calibration(validate_mfrm_calibration(draft))
# save_mfrm_calibration(calibration, "two-family-calibration.rds")
# 
# # In a later R session: the original fit and training data are not needed.
# calibration <- load_mfrm_calibration("two-family-calibration.rds")
# scores <- score_mfrm_calibration(calibration, new_ratings)
# summary(scores)

## ----gpcm-interval-example----------------------------------------------------
# library(mfrmr)
# ratings <- simulate_mfrm_data(
#   n_person = 80, n_rater = 3, n_criterion = 3, score_levels = 3,
#   model = "GPCM", step_facet = "Criterion", slope_facet = "Criterion",
#   slopes = c(C01 = 0.8, C02 = 1, C03 = 1.25), seed = 924090
# )
# fit <- fit_mfrm(
#   ratings, "Person", c("Rater", "Criterion"), "Score",
#   model = "GPCM", method = "MML", step_facet = "Criterion",
#   slope_facet = "Criterion", quad_points = 31, maxit = 400, reltol = 1e-10
# )

## ----gpcm-slope-comparisons---------------------------------------------------
# confint(fit, parm = "slopes", scale = "standardized")
# slope_levels <- fit$slopes$SlopeFacet
# comparison <- matrix(0, 1, length(slope_levels),
#                      dimnames = list("first versus second", slope_levels))
# comparison[1, 1:2] <- c(1, -1)
# ratio <- confint(fit, parm = "slopes", contrasts = comparison)
# difference <- confint(fit, parm = "slopes", contrasts = comparison,
#                       contrast_scale = "difference")
# attr(ratio, "diagnostics")[, c("Estimate", "CI_Lower", "CI_Upper", "PValue")]

## ----gpcm-sandwich-example----------------------------------------------------
# confint(fit, parm = "slopes", method = "sandwich", simultaneous = "bonferroni")
# # If Persons belong to independent schools, supply a complete Person/Cluster
# # mapping and use clusters = schools. Do not cluster separate rows of one Person.

## ----slope-profile, eval=FALSE------------------------------------------------
# # Choose a level from fit$slopes$SlopeFacet before looking for a favorable result.
# ci <- confint(fit, parm = "slopes", method = "profile", slope = "R01")
# print(ci)
# plot(ci, type = "profile")
# # Hide text when composing a larger figure:
# as_ggplot(ci, type = "profile", title = NULL, subtitle = NULL, caption = NULL)
# apa_table(ci, which = "profile_endpoints")
# apa_table(ci, which = "profile_checks", digits = 8)
# 
# # Keep the calculated profile; reporting does not fit the model again.
# res <- mfrm_results(fit, intervals = list(profile = ci), compute = "never")
# report <- mfrm_report(res)
# saveRDS(ci, "one-slope-profile.rds")

## ----gpcm-probability-intervals, fig.alt="Three fitted category probability curves for one rater and criterion, with pointwise calibration-uncertainty ribbons at fixed ability values. Colors and line types distinguish categories; these are not Person-score intervals."----
# grid <- expand.grid(Theta = seq(-3, 3, length.out = 41),
#                     Rater = "R01", Criterion = "C01")
# curves <- mfrm_curve_intervals(fit, grid)
# plot(curves, title = NULL, subtitle = NULL)

## ----gpcm-information-intervals, fig.alt="Fitted information per rating for the same rater and criterion across ability values, with a pointwise calibration-uncertainty ribbon. The curve is not total information across all ratings."----
# information <- mfrm_curve_intervals(fit, grid, type = "information")
# plot(information)

## ----gpcm-saved-inference, fig.alt="Population-SD-standardized discrimination estimates with explicitly selected model-based pointwise 95 percent intervals. The plot concerns discrimination, not rater severity or quality."----
# ci <- confint(fit, parm = "slopes", scale = "standardized")
# results <- mfrm_results(fit, intervals = list(slopes = ci, curves = curves),
#                         include = c("fit", "plots"), compute = "never")
# plot(results, type = "gpcm_slopes", title = NULL)
# apa_table(ci)
# head(plot_data(results, type = "gpcm_curves", component = "table"))

## ----empirical-writing-review, eval=FALSE-------------------------------------
# library(mfrmr)
# data("data.ratings1", package = "sirt")
# writing_ratings <- do.call(rbind, lapply(paste0("k", 1:5), function(k) {
#   data.frame(
#     Person = as.character(data.ratings1$idstud),
#     Rater = as.character(data.ratings1$rater),
#     Criterion = k, Score = data.ratings1[[k]]
#   )
# }))
# writing_review <- describe_mfrm_data(
#   writing_ratings, person = "Person", facets = c("Criterion", "Rater"),
#   score = "Score", rating_min = 0, rating_max = 3,
#   category_policy = "preserve", rater_facet = "Rater"
# )
# writing_review$design_connectivity
# writing_review$structural_missingness$summary
# # Count distinct raters for each student, independently of criterion rows.
# writing_pairs <- unique(writing_ratings[c("Person", "Rater")])
# table(table(writing_pairs$Person))
# xtabs(~ Rater + Score, writing_ratings)

## ----separate-owners, eval=FALSE----------------------------------------------
# fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#   model = "GPCM", method = "MML",
#   step_facet = "Rater", slope_facet = "Criterion")
# slopes <- confint(fit, parm = "slopes")
# grid <- expand.grid(Theta = seq(-2, 2, length.out = 21),
#                     Rater = "R01", Criterion = "C01") # use your fitted labels
# curves <- mfrm_curve_intervals(fit, grid)
# plot(curves, title = NULL, subtitle = NULL)
# results <- mfrm_results(fit, intervals = list(slopes = slopes, curves = curves),
#                         include = c("fit", "plots"), compute = "never")
# saveRDS(results, "criterion-slopes-rater-steps.rds")
# mfrm_report(readRDS("criterion-slopes-rater-steps.rds"))

## ----capability-supported-----------------------------------------------------
# library(mfrmr)
# gpcm_capability_matrix("supported")[, c("Area", "Status")]

## ----capability-with-caveat---------------------------------------------------
# gpcm_capability_matrix("supported_with_caveat")[, c("Area", "Status")]

## ----capability-blocked-------------------------------------------------------
# gpcm_capability_matrix("blocked")[, c("Area", "Status", "RecommendedRoute")]

## ----capability-deferred------------------------------------------------------
# gpcm_capability_matrix("deferred")[, c("Area", "Status", "Boundary", "RecommendedRoute")]

