## ----setup, include=FALSE---------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(eval = !is_cran_check, collapse = TRUE, comment = "#>", cache = FALSE,
                     dev.args = list(bg = "white"))
library(mfrmr)

## ----data-------------------------------------------------------------------
# ratings <- load_mfrmr_data("example_core")
# head(ratings[c("Person", "Rater", "Criterion", "Score")], 8)

## ----fit-recipe, eval=FALSE-------------------------------------------------
# testlet_fit <- fit_mfrm_testlet(
#   ratings, person = "Person", score = "Score", testlet = "Rater",
#   facets = c("Rater", "Criterion"), score_levels = 1:4,
#   quad_points = 121
# )

## ----fit--------------------------------------------------------------------
# saved_example <- readRDS(system.file("examples", "extended-models.rds",
#                                       package = "mfrmr", mustWork = TRUE))$testlet
# testlet_fit <- saved_example$fit
# testlet_fit$checks
# summary(testlet_fit)$data_usage
# testlet_fit$calibration$variance
# testlet_fit$calibration$person_variance
# stopifnot(testlet_fit$checks$NumericalReady,
#           testlet_fit$checks$InformationPositive)

## ----feedback, fig.width=7, fig.height=5, fig.alt="Fixed-rater severity estimates with explicitly requested 95 percent observed-information normal bounds. These approximations have unresolved finite-sample coverage and are not rater-quality classifications."----
# subset(testlet_fit$calibration_table, Facet == "Rater")
# # Bounds require an explicit choice of approximation and nominal level.
# confint(testlet_fit, parm = "calibration", level = .95)
# plot(testlet_fit, facet = "Rater", intervals = "normal", level = .95)

## ----scoring, fig.width=7, fig.height=5, fig.alt="Four Person EAP estimates with conditional equal-tail intervals under the testlet model. Fitted calibration, including ability and local variances, is held fixed."----
# first <- ratings[ratings$Person %in% unique(ratings$Person)[1:4], ]
# scores <- predict(testlet_fit, first)
# scores$table
# plot(scores)
# stopifnot(all(scores$table$Status == "available_conditional"))

## ----missing, fig.width=7, fig.height=5, fig.alt="Person estimates after omitting missing assigned scores. An open circle identifies the entirely missing Person whose result comes only from the fitted population; other points use their remaining observed ratings."----
# incomplete <- first
# last_person <- tail(unique(incomplete$Person), 1)
# incomplete$Score[incomplete$Person == last_person] <- NA_real_
# incomplete$Score[1] <- NA_real_
# review <- predict(testlet_fit, incomplete, missing = "omit")
# review$table
# review$data_usage
# plot(review)
# stopifnot(sum(review$table$Status == "prior_only") == 1)

## ----saved------------------------------------------------------------------
# path <- tempfile(fileext = ".rds")
# saveRDS(testlet_fit, path)
# restored <- readRDS(path)
# stopifnot(identical(restored$parameters, testlet_fit$parameters))
# saved_plot <- plot(scores, draw = FALSE)
# plot_data(saved_plot)$table

## ----shared-accessors-------------------------------------------------------
# plot_data_components(scores)
# score_table <- plot_data(scores, component = "table")
# plot_data(testlet_fit, component = "table", facet = "Rater")
# stopifnot(identical(score_table, scores$table))

## ----shared-graphics, fig.width=7, fig.height=5, fig.alt="The same incomplete-roster Person estimates and conditional intervals rendered with ggplot2. The open point marks a prior-only result, not measured average performance."----
# if (requireNamespace("ggplot2", quietly = TRUE)) {
#   figure <- as_ggplot(review)
#   print(figure + ggplot2::theme_minimal(base_size = 12))
#   # For a file: ggplot2::ggsave("person-scores.png", figure,
#   #                            width = 7, height = 5, dpi = 300)
# }

## ----table-export-----------------------------------------------------------
# csv_path <- tempfile(fileext = ".csv")
# utils::write.csv(review$table, csv_path, row.names = FALSE, na = "NA")
# score_path <- tempfile(fileext = ".rds")
# saveRDS(review, score_path)
# stopifnot(identical(readRDS(score_path), review))

## ----ordinary-reference-----------------------------------------------------
# ordinary_fit <- fit_mfrm(
#   ratings, person = "Person", facets = c("Rater", "Criterion"), score = "Score",
#   model = "RSM", method = "MML", rating_min = 1, rating_max = 4,
#   keep_original = TRUE, population_formula = ~1,
#   person_data = unique(ratings["Person"]), quad_points = 121
# )

## ----model-comparison, fig.width=7, fig.height=5, fig.alt="Differences between centered rater effects from the testlet and ordinary RSM, plotted against their means. The horizontal zero line is descriptive, with no difference intervals or automatic decision threshold."----
# model_comparison <- compare_mfrm(
#   ordinary_fit, testlet_fit, labels = c("Ordinary RSM", "Testlet RSM")
# )
# model_comparison$checks
# model_comparison$effects
# plot(model_comparison, facet = "Rater", style = "difference")

## ----response-diagnostics, fig.width=8, fig.height=5, fig.alt="Paired descriptive posterior predictive Infit and Outfit for each rater; symbols distinguish the indices and no reference cutoffs are drawn."----
# response_review <- mfrm_response_diagnostics(testlet_fit, group_by = "Rater")
# response_review$measures
# plot(response_review, style = "paired")

## ----predictive-comparison, fig.width=8, fig.height=5, fig.alt="Differences in descriptive posterior predictive Infit for matched rater summaries. Zero means agreement between models, not acceptable model fit."----
# ordinary_response_review <- mfrm_response_diagnostics(ordinary_fit,
#   group_by = "Rater")
# model_comparison <- compare_mfrm(ordinary_fit, testlet_fit,
#   labels = c("Ordinary RSM", "Testlet RSM"),
#   response_diagnostics = list(ordinary_response_review, response_review))
# model_comparison$responses$measures
# plot(model_comparison, metric = "infit", style = "difference")

## ----predictive-probabilities, fig.width=9, fig.height=5.5, fig.alt="Paired category probabilities from the two models, one panel per score category. The dashed diagonal marks equal probabilities; this monochrome display has no quality flags."----
# as_ggplot(model_comparison, metric = "probability", show_labels = FALSE,
#   palette = "mono", show_title = FALSE, show_notes = FALSE)

## ----person-comparison, fig.width=8, fig.height=5, fig.alt="Changes in conditional Person EAPs after centering each model at its fitted population mean. There are no model-difference confidence intervals or significance decisions."----
# source_scores <- score_mfrm_persons(testlet_fit, persons = scores$table$Person)
# ordinary_person_scores <- score_mfrm_persons(ordinary_fit,
#   persons = source_scores$table$Person)
# model_comparison <- compare_mfrm(ordinary_fit, testlet_fit,
#   response_diagnostics = list(ordinary_response_review, response_review),
#   person_scores = list(ordinary_person_scores, source_scores))
# model_comparison$persons$table
# plot(model_comparison, metric = "person", style = "difference", show_labels = TRUE)

## ----model-wright, fig.width=10, fig.height=7, fig.alt="Model-aware Wright map of selected Person EAPs, conditional facet reference locations and adjacent category boundaries. Only Person points have conditional interval whiskers."----
# map_results <- mfrm_results(testlet_fit, scores = source_scores,
#   response_diagnostics = response_review, comparison = model_comparison)
# plot(map_results, type = "wright", show_labels = TRUE)

## ----model-fit-pathway, fig.width=8, fig.height=5.5, fig.alt="Rater reference locations against descriptive posterior Infit for the selected ratings. Rater labels identify points; the monochrome plot has no acceptable-fit bands or reference value of one."----
# as_ggplot(plot(map_results, type = "fit_pathway", facet = "Rater",
#   palette = "mono", show_title = FALSE, draw = FALSE))

## ----report-----------------------------------------------------------------
# res <- mfrm_results(testlet_fit, scores = source_scores, comparison = model_comparison,
#                     diagnostics = response_review)
# res$tables$interval_basis
# res$tables$scoring_status
# report <- mfrm_report(res)
# summary(report, view = "reader")
# archive_dir <- tempfile("testlet-report-")
# archive <- export_mfrm_results(res, archive_dir, preset = "starter",
#   acknowledge_sensitive = TRUE) # Synthetic example; real exports retain IDs
# stopifnot(nrow(archive$plot_errors) == 0)
# restored_results <- readRDS(file.path(archive_dir, "mfrmr_results_results.rds"))
# stopifnot(identical(restored_results$tables$person_scores, source_scores$table))

