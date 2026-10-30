## ----setup, include=FALSE-----------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
rtmb_available <- !is_cran_check && requireNamespace("RTMB", quietly = TRUE) && packageVersion("RTMB") >= "2.0"
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", dev.args = list(bg = "white"),
                     eval = !is_cran_check && rtmb_available)
library(mfrmr)

## ----data---------------------------------------------------------------------
# ratings <- load_mfrmr_data("example_core")

## ----fit-recipe, eval=FALSE---------------------------------------------------
# random_fit <- fit_mfrm_random_rater(
#   ratings, person = "Person", rater = "Rater", score = "Score",
#   facets = "Criterion", score_levels = 1:4, quad_points = 121
# )

## ----fit----------------------------------------------------------------------
# saved_example <- readRDS(system.file("examples", "extended-models.rds",
#                                       package = "mfrmr", mustWork = TRUE))$random_rater
# random_fit <- saved_example$fit
# random_fit$checks
# summary(random_fit)$data_usage

## ----population---------------------------------------------------------------
# random_fit$calibration[c("person_sd", "person_variance", "rater_sd")]

## ----calibration-intervals----------------------------------------------------
# confint(random_fit, parm = "calibration", level = .95)
# summary(random_fit, calibration_intervals = "normal", level = .95)$calibration

## ----quadrature-refit, eval=FALSE---------------------------------------------
# refined <- fit_mfrm_random_rater(
#   ratings, "Person", "Rater", "Score", facets = "Criterion",
#   score_levels = 1:4, person_sd = NULL, quad_points = 121
# )
# refined$checks

## ----severity, fig.width=7, fig.height=5, fig.alt="Conditional-mode severity estimates for four observed raters relative to the assumed population mean. Positive values indicate stricter ratings. The default shows points without interval bounds or quality classifications."----
# random_fit$raters
# plot(random_fit)

## ----explicit-normal----------------------------------------------------------
# normal_approximation <- confint(random_fit, parm = "raters")
# normal_approximation
# # Optional display: plot(random_fit, intervals = "normal")

## ----severity-precision, fig.width=7, fig.height=5, fig.alt="Four observed rater estimates plotted against explicitly requested approximate interval width. Saved plot data provide the exact estimates and bounds; nominal coverage is not established."----
# plot(random_fit, style = "precision", intervals = "normal", point_size = 3,
#      title = "Rater severity and approximate interval width")

## ----severity-distribution, fig.width=7, fig.height=5, fig.alt="Empirical cumulative distribution of four observed rater estimates. The shaded step curve describes point estimates, not a latent population distribution."----
# plot(random_fit, style = "distribution", reference = NULL,
#      title = "Distribution of the observed rater estimates")

## ----severity-custom, fig.width=7, fig.height=4, fig.alt="Monochrome rater point estimates in increasing estimate order, with the title and caption hidden. Text alternatives and the full table follow the figure."----
# figure <- plot(random_fit, draw = FALSE, sort = "estimate", palette = "mono",
#                show_title = FALSE, show_notes = FALSE, text_scale = 1.15)
# if (requireNamespace("ggplot2", quietly = TRUE)) print(as_ggplot(figure))
# # Text alternatives and the full table can accompany the exported image:
# plot_data(figure)$alt_text
# plot_data(figure)$table

## ----ordinary-reference-------------------------------------------------------
# ordinary_fit <- fit_mfrm(
#   ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score",
#   model = "RSM", method = "MML", rating_min = 1, rating_max = 4,
#   keep_original = TRUE, population_formula = ~1,
#   person_data = unique(ratings["Person"]), quad_points = 121
# )

## ----model-comparison, fig.width=7, fig.height=5, fig.alt="Fixed versus shared random rater effects, centered over the same observed panel. Points show descriptive changes relative to the equality line; no difference intervals or model preference are shown."----
# model_comparison <- compare_mfrm(
#   ordinary_fit, random_fit, labels = c("Fixed raters", "Shared random raters")
# )
# model_comparison$checks
# model_comparison$models
# subset(model_comparison$effects, Facet == "Rater")
# plot(model_comparison, facet = "Rater", style = "paired")

## ----profile------------------------------------------------------------------
# sd_interval <- confint(random_fit, parm = "rater_sd")
# sd_interval[, c("Lower", "Upper"), drop = FALSE]
# head(attr(sd_interval, "profile"))

## ----bootstrap-recipe, eval=FALSE---------------------------------------------
# # A small run checks the workflow; 19 draws cannot stabilize 2.5% tails.
# bootstrap_intervals <- mfrm_random_rater_intervals(
#   random_fit, nsim = 19, seed = 923701
# )

## ----bootstrap-intervals------------------------------------------------------
# bootstrap_intervals <- saved_example$intervals
# stopifnot(identical(bootstrap_intervals$source, random_fit))
# bootstrap_summary <- summary(bootstrap_intervals)
# bootstrap_summary$trials
# bootstrap_summary$availability
# bootstrap_summary$intervals[, c("Lower", "Upper"), drop = FALSE]
# confint(bootstrap_intervals, method = "error")[, c("Lower", "Upper"), drop = FALSE]
# # Reuse the same roots at a different level, with no further fitting:
# confint(bootstrap_intervals, level = .90)[, c("Lower", "Upper"), drop = FALSE]

## ----bootstrap-plot, fig.width=7, fig.height=5, fig.alt="Explicit bootstrap intervals for observed rater effects. Arrowheads mark unbounded limits caused by unresolved studentized errors; all planned refits remain in the adjacent trial table."----
# plot(bootstrap_intervals)
# bootstrap_intervals$trials[
#   !bootstrap_intervals$trials$FitReady |
#     bootstrap_intervals$trials$EstimatedBoundary, ]

## ----person-selection---------------------------------------------------------
# # Keep the complete source roster; request only these output rows.
# selected_persons <- as.character(unique(ratings$Person)[1:2])

## ----person-scoring-recipe, eval=FALSE----------------------------------------
# person_scores <- score_mfrm_random_rater(random_fit, persons = selected_persons)

## ----person-scoring-----------------------------------------------------------
# person_scores <- saved_example$scores
# stopifnot(identical(person_scores$table$Person, selected_persons))
# person_scores$table
# summary(person_scores)$data_usage

## ----person-score-plot, fig.width=8, fig.height=5, fig.alt="Conditional Person EAP estimates and equal-tail intervals using all shared-rater evidence in the scoring roster. Calibration is held fixed; these are not intervals for Person differences."----
# plot(person_scores)

## ----predictions--------------------------------------------------------------
# observed <- data.frame(
#   Rater = "R01", Criterion = "Content", Ability = c(-1, 0, 1)
# )
# p_observed <- predict(random_fit, observed, ability = "Ability", rater = "observed")
# replacement <- transform(observed, Rater = "New rater")
# p_new <- predict(random_fit, replacement, ability = "Ability", rater = "new")
# p_observed$expected_scores
# p_new$expected_scores
# p_new$probabilities

## ----response-diagnostics, fig.width=8, fig.height=5, fig.alt="Paired descriptive Infit and Outfit by rater for the selected Person's ratings, conditioned on the full roster. No fit thresholds or rater-quality flags are shown."----
# selected_rows <- which(ratings$Person == ratings$Person[1])
# response_review <- mfrm_response_diagnostics(random_fit,
#   rows = selected_rows, group_by = "Rater")
# response_review$measures
# plot(response_review, style = "paired")

## ----predictive-comparison, fig.width=8, fig.height=5, fig.alt="Differences in descriptive posterior predictive Infit for matched rater summaries. Zero means agreement between models, not acceptable model fit."----
# ordinary_response_review <- mfrm_response_diagnostics(ordinary_fit,
#   rows = selected_rows, group_by = "Rater")
# model_comparison <- compare_mfrm(ordinary_fit, random_fit,
#   labels = c("Ordinary RSM", "Shared-rater RSM"),
#   response_diagnostics = list(ordinary_response_review, response_review))
# model_comparison$responses$measures
# plot(model_comparison, metric = "infit", style = "difference")

## ----predictive-probabilities, fig.width=9, fig.height=5.5, fig.alt="Paired category probabilities from the two models, one panel per score category. The dashed diagonal marks equal probabilities; this monochrome display has no quality flags."----
# as_ggplot(model_comparison, metric = "probability", show_labels = FALSE,
#   palette = "mono", show_title = FALSE, show_notes = FALSE)

## ----person-comparison, fig.width=8, fig.height=5, fig.alt="Changes in conditional Person EAPs after centering each model at its fitted population mean. There are no model-difference confidence intervals or significance decisions."----
# source_scores <- person_scores # Already scored from the complete roster
# ordinary_person_scores <- score_mfrm_persons(ordinary_fit,
#   persons = source_scores$table$Person)
# model_comparison <- compare_mfrm(ordinary_fit, random_fit,
#   response_diagnostics = list(ordinary_response_review, response_review),
#   person_scores = list(ordinary_person_scores, source_scores))
# model_comparison$persons$table
# plot(model_comparison, metric = "person", style = "difference", show_labels = TRUE)

## ----model-wright, fig.width=10, fig.height=7, fig.alt="Model-aware Wright map of selected Person EAPs, conditional facet reference locations and adjacent category boundaries. Only Person points have conditional interval whiskers."----
# map_results <- mfrm_results(random_fit, scores = source_scores,
#   response_diagnostics = response_review, comparison = model_comparison)
# plot(map_results, type = "wright", show_labels = TRUE)

## ----model-fit-pathway, fig.width=8, fig.height=5.5, fig.alt="Rater reference locations against descriptive posterior Infit for the selected ratings. Rater labels identify points; the monochrome plot has no acceptable-fit bands or reference value of one."----
# as_ggplot(plot(map_results, type = "fit_pathway", facet = "Rater",
#   palette = "mono", show_title = FALSE, draw = FALSE))

## ----report-------------------------------------------------------------------
# res <- mfrm_results(random_fit, predictions = p_new, scores = person_scores,
#                     intervals = bootstrap_intervals, comparison = model_comparison,
#                     diagnostics = response_review)
# res$tables$interval_basis
# res$tables$person_scores
# res$tables$bootstrap_availability
# res$tables$comparison_effects
# report <- mfrm_report(res)
# summary(report, view = "reader")
# archive_dir <- tempfile("random-rater-report-")
# archive <- export_mfrm_results(res, archive_dir, preset = "starter",
#   acknowledge_sensitive = TRUE) # Synthetic example; real exports retain IDs
# stopifnot(nrow(archive$plot_errors) == 0)
# restored_results <- readRDS(file.path(archive_dir, "mfrmr_results_results.rds"))
# stopifnot(identical(restored_results$tables$bootstrap_trials,
#                     bootstrap_intervals$trials))

## ----save---------------------------------------------------------------------
# path <- tempfile(fileext = ".rds")
# saveRDS(list(ratings = ratings, fit = random_fit, population_interval = sd_interval,
#              person_scores = person_scores, rater_intervals = bootstrap_intervals,
#              observed = p_observed, replacement = p_new), path)
# saved <- readRDS(path)
# stopifnot(identical(predict(saved$fit, replacement, "Ability", rater = "new"), p_new))
# plot_data(plot(saved$fit, draw = FALSE))$table

