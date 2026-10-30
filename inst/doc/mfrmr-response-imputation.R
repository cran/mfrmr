## ----setup, include=FALSE---------------------------------------------------
is_cran_check <- !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
knitr::opts_chunk$set(eval = !is_cran_check, collapse = TRUE, comment = "#>")

## ----roster-----------------------------------------------------------------
# library(mfrmr)
# example <- readRDS(system.file("examples", "response-imputation.rds", package = "mfrmr"))
# ratings <- example$ratings
# missing_assigned <- example$missing
# unassigned <- example$unassigned
# held_out_scores <- example$held_out
# impute_ids <- ratings$Event[missing_assigned]
# m <- ncol(example$scores)
# model <- example$imputation_model
# with(ratings, table(Assigned, MissingScore = is.na(Score)))

## ----sampling-diagnostics---------------------------------------------------
# diagnostics <- model$diagnostics
# c(MaxRhat = max(diagnostics$rhat),
#   MinBulkESS = min(diagnostics$ess_bulk),
#   MinTailESS = min(diagnostics$ess_tail))
# model$sampler
# head(model$draw_ids)

## ----regenerate, eval=FALSE-------------------------------------------------
# source(system.file("examples", "response-imputation.R", package = "mfrmr"))
# example <- response_imputation_example("mi-example-output")

## ----completions------------------------------------------------------------
# completed <- lapply(seq_len(m), function(i) {
#   result <- ratings
#   result$Score[missing_assigned] <- example$scores[, i]
#   result
# })
# review <- review_mfrm_imputations(
#   ratings, completed, person = "Person", facets = c("Rater", "Criterion"),
#   score = "Score", event_id = "Event", impute_ids = impute_ids, categories = 1:4,
#   assigned = "Assigned", imputation_model = model
# )
# review
# subset(summary(review), Facet == "Rater")
# stopifnot(all(vapply(review$completed, function(d) {
#   all(is.na(d$Score[unassigned])) &&
#     identical(d$Score[!is.na(ratings$Score)], ratings$Score[!is.na(ratings$Score)])
# }, logical(1))))

## ----imputation-diagnostics, fig.width=6.5, fig.height=8, fig.alt="Category probabilities for thirty missing assigned ratings; each cell shows its percentage."----
# probabilities <- example$probability_mean
# n <- nrow(probabilities)
# image(1:4, seq_len(n), t(probabilities[n:1, , drop = FALSE]),
#       zlim = c(0, 1), col = colorRampPalette(c("#ffffd9", "#7fcdbb", "#41b6c4"))(64), axes = FALSE,
#       xlab = "Score", ylab = "", main = "")
# axis(1, at = 1:4)
# axis(2, at = seq_len(n), labels = rev(impute_ids), las = 1, cex.axis = .7)
# for (j in 1:4) text(j, seq_len(n),
#   labels = sprintf("%.0f%%", probabilities[n:1, j] * 100), cex = .65,
#   col = "black")

## ----held-out-check---------------------------------------------------------
# prediction_check <- data.frame(
#   Event = impute_ids, HeldOutScore = held_out_scores,
#   PredictedMean = drop(probabilities %*% (1:4)),
#   Brier = rowSums((probabilities - outer(held_out_scores, 1:4, "=="))^2)
# )
# head(prediction_check)
# c(HeldOutMean = mean(held_out_scores),
#   PredictedMean = mean(prediction_check$PredictedMean),
#   Brier = mean(prediction_check$Brier))

## ----analyses---------------------------------------------------------------
# analyses <- fit_mfrm_imputed(review, model = "RSM", quad_points = 61)
# summary(analyses)[, c("Imputation", "Status", "Error", "Warnings")]

## ----pooled-result----------------------------------------------------------
# levels <- analyses$fits[[1]]$config$facet_specs$Rater$levels
# contrast <- matrix(0, 1, length(levels),
#                    dimnames = list("R04 minus R01", levels))
# contrast[1, c("R04", "R01")] <- c(1, -1)
# pooled <- pool_mfrm_imputed(analyses, facet = "Rater", contrasts = contrast)
# summary(pooled)

## ----pooled-plot, fig.width=7, fig.height=3, fig.alt="Pooled R04-minus-R01 severity contrast with a pointwise Rubin t interval incorporating within- and between-imputation covariance. The result is conditional on the supplied imputations and analysis assumptions."----
# plot(pooled)

## ----pooled-ggplot, fig.width=7, fig.height=4, fig.alt="The saved R04-minus-R01 contrast and its pointwise multiple-imputation t interval, rendered with ggplot. Changing the title does not recalculate the interval or change its confidence level."----
# if (requireNamespace("ggplot2", quietly = TRUE)) {
#   pooled_figure <- as_ggplot(pooled, title = NULL, preset = "monochrome")
#   print(pooled_figure)
#   plot_data(pooled_figure)$table[, c("Target", "Estimate", "Lower", "Upper", "DF")]
# }

## ----observed-comparison----------------------------------------------------
# observed_ratings <- subset(ratings, Assigned & !is.na(Score))
# observed_fit <- fit_mfrm(observed_ratings, "Person", c("Rater", "Criterion"),
#   "Score", model = "RSM", method = "MML", rating_min = 1, rating_max = 4,
#   quad_points = 121, keep_original = TRUE, attach_diagnostics = FALSE)
# observed_interval <- mfrm_facet_intervals(observed_fit, "Rater", contrasts = contrast)
# comparison <- rbind(
#   data.frame(Method = "Observed-score MML",
#              summary(observed_interval)[c("Estimate", "SE", "Lower", "Upper")]),
#   example$reference,
#   data.frame(Method = "Joint-model MI",
#              summary(pooled)[c("Estimate", "SE", "Lower", "Upper")])
# )
# comparison

## ----sensitivity------------------------------------------------------------
# lower_scores <- lapply(completed, function(d) {
#   d$Score[missing_assigned] <- pmax(1L, d$Score[missing_assigned] - 1L)
#   d
# })
# lower_review <- review_mfrm_imputations(
#   ratings, lower_scores, person = "Person", facets = c("Rater", "Criterion"),
#   score = "Score", event_id = "Event", impute_ids = impute_ids, categories = 1:4,
#   assigned = "Assigned", imputation_model = list(
#     baseline = model, assumption = "Missing scores one category lower, floored at 1")
# )
# lower_analyses <- fit_mfrm_imputed(lower_review, model = "RSM", quad_points = 61)
# lower_pooled <- pool_mfrm_imputed(lower_analyses, facet = "Rater", contrasts = contrast)
# comparison <- rbind(
#   cbind(Assumption = "Joint RSM baseline", summary(pooled)),
#   cbind(Assumption = "Missing scores one category lower", summary(lower_pooled))
# )
# comparison[, c("Assumption", "Estimate", "SE", "Lower", "Upper", "MonteCarloSE")]
# paired_change <- drop(lower_pooled$estimates - pooled$estimates)
# c(ContrastChange = mean(paired_change),
#   MonteCarloSEOfChange = sd(paired_change) / sqrt(m))

## ----save, eval=FALSE-------------------------------------------------------
# saveRDS(pooled, "rater-imputation-result.rds")
# saved <- readRDS("rater-imputation-result.rds")
# summary(saved)
# write.csv(summary(saved), "rater-imputation-table.csv", row.names = FALSE)
# saved$analyses$imputations$imputation_model$diagnostics
# saveRDS(list(baseline = pooled, lower_scores = lower_pooled),
#         "rater-imputation-sensitivity.rds")

