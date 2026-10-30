## ----setup, include=FALSE-----------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(eval = !is_cran_check, collapse = TRUE, comment = "#>")

## ----fit----------------------------------------------------------------------
# library(mfrmr)
# ratings <- load_mfrmr_data("example_core")
# fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#                 model = "RSM", method = "MML")
# ordinary <- mfrm_facet_intervals(fit, "Rater")
# summary(ordinary)

## ----contrast-----------------------------------------------------------------
# contrast <- matrix(c(1, -1, 0, 0), nrow = 1,
#   dimnames = list("R01 minus R02", c("R01", "R02", "R03", "R04")))
# difference <- mfrm_facet_intervals(fit, "Rater", contrasts = contrast,
#                                   method = "sandwich")
# summary(difference)

## ----contrast-plot, fig.width=7, fig.height=4.2, fig.alt="The same estimated rater contrast with model-based and one-way Person-sandwich normal intervals, offset vertically for comparison. Zero denotes no contrast, not a practical-importance or rater-quality threshold."----
# plot(difference)

## ----cml-contrast-------------------------------------------------------------
# if (requireNamespace("lpSolve", quietly = TRUE)) {
#   cml_ratings <- load_mfrmr_data("example_core")
#   cml_fit <- fit_mfrm(cml_ratings, "Person", c("Rater", "Criterion"), "Score",
#     method = "CML", model = "RSM", rating_min = 1, rating_max = 4,
#     category_policy = "preserve")
#   cml_contrast <- matrix(c(1, -1, 0, 0), nrow = 1,
#     dimnames = list("R01 minus R02", c("R01", "R02", "R03", "R04")))
#   cml_interval <- mfrm_facet_intervals(cml_fit, "Rater", cml_contrast)
#   print(summary(cml_interval))
#   print(summary(cml_fit)$steps)
#   cml_fit$checks
# }

## ----cml-contrast-plot, fig.width=7, fig.height=4.2, fig.alt="A CML estimate of the R01 minus R02 severity contrast with a pointwise 95 percent conditional-information Wald interval on the native logit scale."----
# if (requireNamespace("lpSolve", quietly = TRUE)) plot(cml_interval, comparison = FALSE)

## ----individual-sheet---------------------------------------------------------
# diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
# feedback_results <- mfrm_results(fit, diagnostics = diagnostics,
#   intervals = list(raters = ordinary), compute = "never")
# sheet <- mfrm_report(feedback_results, style = "rater", facet = "Rater",
#   rater = "R01", interval = "raters", output = "html", max_cases = 5)
# sheet$report$tables$severity
# sheet$report$tables$uncertainty
# # Open sheet$path in a browser. Review before copying the HTML file for sharing.

## ----feedback-save-reopen-----------------------------------------------------
# analysis_file <- tempfile(fileext = ".rds")
# saveRDS(feedback_results, analysis_file)
# saved_results <- readRDS(analysis_file)
# 
# # Reuse the saved diagnostics and intervals for this recipient.
# saved_sheet <- mfrm_report(saved_results, style = "rater", facet = "Rater",
#   rater = "R01", interval = "raters", output = "html", max_cases = 5)
# recipient_file <- tempfile(fileext = ".html")
# stopifnot(file.copy(saved_sheet$path, recipient_file, overwrite = FALSE))

## ----feedback-table-files-----------------------------------------------------
# sheet_tables <- mfrm_report(saved_results, style = "rater", facet = "Rater",
#   rater = "R01", interval = "raters", output = "tables", max_cases = 5)
# sheet_tables$interpretation_notes
# table_dir <- tempfile("rater-tables-")
# dir.create(table_dir)
# for (name in names(sheet_tables)) {
#   if (ncol(sheet_tables[[name]]) > 0) {
#     write.csv(sheet_tables[[name]], file.path(table_dir, paste0(name, ".csv")),
#       row.names = FALSE)
#   }
# }

## ----saved-intervals----------------------------------------------------------
# results <- mfrm_results(fit, intervals = list(difference = difference),
#                         include = c("fit", "plots"), compute = "never")
# apa_table(difference)
# plot_data(results, type = "facet_difference", component = "table")
# report <- mfrm_report(results)
# report$tables$facet_difference_intervals

## ----saved-interval-plot, fig.width=7, fig.height=4.2, fig.alt="A saved rater-severity difference displayed with model-based and sandwich pointwise normal intervals. Vertical offsets and different line types distinguish the methods without relying only on color. The title and subtitle have been explicitly omitted."----
# as_ggplot(results, type = "facet_difference", title = NULL, subtitle = NULL,
#           preset = "monochrome")

## ----interval-archive---------------------------------------------------------
# archive <- export_mfrm_results(feedback_results,
#   output_dir = tempfile("rater-analysis-"),
#   include = c("tables", "report", "replay", "manifest"),
#   acknowledge_sensitive = TRUE)
# # Count files by format; inspect archive$written_files for the full index.
# table(archive$written_files$Format)
# archive$plot_errors
# saved_file <- archive$written_files$Path[
#   archive$written_files$Component == "results_rds"]
# reopened <- readRDS(saved_file)
# stopifnot(identical(reopened$facet_intervals, feedback_results$facet_intervals))

## ----pcm-saved-intervals------------------------------------------------------
# fit_pcm <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#                     model = "PCM", step_facet = "Criterion", method = "MML")
# pcm_difference <- mfrm_facet_intervals(fit_pcm, "Rater", contrasts = contrast,
#                                       method = "sandwich", level = 0.90)
# pcm_results <- mfrm_results(fit_pcm, intervals = list(difference = pcm_difference),
#                             include = c("fit", "plots"), compute = "never")
# apa_table(pcm_difference)

## ----clusters-----------------------------------------------------------------
# persons <- sort(unique(ratings$Person))
# schools <- data.frame(Person = persons,
#   Cluster = rep(sprintf("School%02d", 1:16), each = 3))
# school_difference <- mfrm_facet_intervals(
#   fit, "Rater", contrasts = contrast, method = "sandwich", clusters = schools
# )
# school_difference$settings[c("persons", "clusters", "cluster_score_rank")]
# summary(school_difference)

## ----unavailable--------------------------------------------------------------
# two_groups <- schools
# two_groups$Cluster <- rep(c("A", "B"), each = 24)
# unavailable <- mfrm_facet_intervals(
#   fit, "Rater", contrasts = contrast, method = "sandwich", clusters = two_groups
# )
# summary(unavailable)[, c("Target", "Estimate", "Lower", "Upper", "Status")]

## ----save---------------------------------------------------------------------
# saved <- tempfile(fileext = ".rds")
# saveRDS(difference, saved)
# restored <- readRDS(saved)
# summary(restored)
# plotted <- plot(restored, draw = FALSE)
# plot_data(plotted)$table

