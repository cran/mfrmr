## ----setup, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", dev.args = list(bg = "white"))

## ----ratings----------------------------------------------------------------
library(mfrmr)

data("mfrmr_example_agreement", package = "mfrmr")
ratings <- mfrmr_example_agreement
head(ratings)
with(ratings, table(Rater, Score, useNA = "ifany"))

## ----check-keys-------------------------------------------------------------
key <- c("Person", "Task", "Rater")
stopifnot(anyDuplicated(ratings[key]) == 0L)

## ----pairwise-agreement-----------------------------------------------------
qwk <- pairwise_kappa_table(
  ratings, rater = "Rater", unit = c("Person", "Task"), score = "Score",
  categories = 0:3, weights = "quadratic"
)
qwk$table[, c("Rater1", "Rater2", "N", "Kappa", "Status")]

## ----pair-audit-------------------------------------------------------------
pair_audit <- cbind(qwk$cases, qwk$units[qwk$cases$Unit, , drop = FALSE])
subset(pair_audit, !Included,
       select = c(Person, Task, Rater1, Rater2,
                  Assigned1, Assigned2, Missing1, Missing2, Included))

## ----align-sources----------------------------------------------------------
scored <- subset(ratings, Rater == "R2", select = c(Person, Task, Score))
reference <- subset(ratings, Rater == "R1", select = c(Person, Task, Score))
names(reference)[3] <- "Reference"
stopifnot(anyDuplicated(scored[c("Person", "Task")]) == 0L,
          anyDuplicated(reference[c("Person", "Task")]) == 0L)
aligned <- merge(scored, reference, by = c("Person", "Task"), all = TRUE)

comparison <- score_comparison_table(
  scores = aligned$Score,
  reference = aligned$Reference,
  categories = 0:3
)
comparison$table[, c("N", "MeanScore", "MeanReference", "QWK", "Correlation", "SMD")]
mean(aligned$Score == aligned$Reference) # Exact agreement for this complete pair

## ----comparison-status------------------------------------------------------
comparison$criteria[, c("Metric", "Value", "Threshold", "Status", "Reason")]
comparison_cases <- cbind(aligned[c("Person", "Task")], comparison$cases)
comparison_cases

## ----second-reference-------------------------------------------------------
reference2 <- subset(ratings, Rater == "R3", select = c(Person, Task, Score))
names(reference2)[3] <- "Reference2"
aligned3 <- merge(aligned, reference2, by = c("Person", "Task"), all = TRUE)
with_baseline <- score_comparison_table(
  scores = aligned3$Score, reference = aligned3$Reference,
  reference2 = aligned3$Reference2, categories = 0:3
)
with_baseline$table[, c("NInput", "N", "NExcluded", "QWK", "ReferenceQWK", "QWKDegradation")]
cbind(aligned3[c("Person", "Task")], with_baseline$cases)

## ----panel-agreement--------------------------------------------------------
alpha <- krippendorff_alpha_table(
  ratings, rater = "Rater", unit = c("Person", "Task"), score = "Score",
  level = "ordinal", categories = 0:3
)
ac2 <- gwet_ac2_table(
  ratings, rater = "Rater", unit = c("Person", "Task"), score = "Score",
  categories = 0:3, weights = "quadratic"
)
alpha$table[, c("Alpha", "Level", "Status", "PairableUnits", "PairableRatings")]
ac2$table[, c("AC", "Coefficient", "Status", "PairableUnits", "NonmissingRatings")]
alpha$unit_summary

## ----score-count-figure, fig.width=7, fig.height=4.5, fig.alt="R1 versus R2 score counts. Four responses occur at each of (0, 1), (1, 2) and (2, 3). All lie one category above the equal-score diagonal; none has identical scores."----
can_plot <- requireNamespace("ggplot2", quietly = TRUE)
figures <- list()
selected_pair <- subset(qwk$table, Rater1 == "R1" & Rater2 == "R2")
score_counts <- subset(qwk$counts, Pair == selected_pair$Pair)
category_index <- qwk$categories$CategoryIndex
category_labels <- qwk$categories$Category

if (can_plot) {
  figures$score_counts <- ggplot2::ggplot(
    score_counts, ggplot2::aes(x = Category1, y = Category2)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                        colour = "grey45") +
    ggplot2::geom_point(ggplot2::aes(size = N), shape = 21,
                        fill = "grey85", colour = "grey30") +
    ggplot2::geom_text(ggplot2::aes(label = N), size = 3.5) +
    ggplot2::scale_size_area(max_size = 14, guide = "none") +
    ggplot2::scale_x_continuous(breaks = category_index, labels = category_labels) +
    ggplot2::scale_y_continuous(breaks = category_index, labels = category_labels) +
    ggplot2::coord_equal(xlim = range(category_index) + c(-0.4, 0.4),
                        ylim = range(category_index) + c(-0.4, 0.4)) +
    ggplot2::labs(x = "R1 score", y = "R2 score", title = NULL, caption = NULL) +
    ggplot2::theme_classic(base_size = 12)
  print(figures$score_counts)
}

## ----pairwise-qwk-figure, fig.width=7, fig.height=3.5, fig.alt="Pairwise quadratic weighted kappa: R1/R2, 0.571 on 12 responses; R1/R3, 0.833 on 10; R2/R3, 0.455 on 10. Points have no uncertainty bars. The vertical dotted line marks zero, not a pass threshold."----
pair_plot_data <- qwk$table
pair_plot_data$PairLabel <- with(pair_plot_data,
  paste0(Rater1, " / ", Rater2, " (N = ", N, ")",
         ifelse(Status == "available", "", paste0("; ", Status))))
pair_plot_data$PairLabel <- factor(pair_plot_data$PairLabel,
                                 levels = rev(pair_plot_data$PairLabel))

if (can_plot) {
  figures$pairwise_qwk <- ggplot2::ggplot(
    pair_plot_data, ggplot2::aes(x = Kappa, y = PairLabel)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey45") +
    ggplot2::geom_point(data = subset(pair_plot_data,
      Status == "available" & is.finite(Kappa)), size = 3) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::expand_limits(x = c(0, 1)) +
    ggplot2::labs(x = "Quadratic weighted kappa", y = NULL,
                  title = NULL, caption = NULL) +
    ggplot2::theme_classic(base_size = 12)
  print(figures$pairwise_qwk)
}

## ----figure-notes-----------------------------------------------------------
figure_titles <- c(score_counts = "Matched score counts for R1 and R2",
                  pairwise_qwk = "Pairwise quadratic weighted agreement")
figure_notes <- list(
  score_counts = c(
    "Circle area and printed counts describe common complete responses; the diagonal marks equal scores.",
    "Blank cells are zero retained counts, not missing-rating indicators. Inspect qwk$cases for exclusions.",
    "The reference is not assumed correct; observed score differences are not model-adjusted severity."
  ),
  pairwise_qwk = c(qwk$notes,
    "N is pair-specific. The zero line is a chance reference, not an acceptance threshold; no intervals are shown.")
)
figure_titles
figure_notes

## ----save-and-reopen--------------------------------------------------------
analysis <- list(
  ratings = ratings, qwk = qwk, alpha = alpha, ac2 = ac2,
  aligned = aligned, comparison = comparison,
  aligned3 = aligned3, with_baseline = with_baseline,
  figures = figures, figure_titles = figure_titles, figure_notes = figure_notes,
  figure_data = list(score_counts = score_counts, pairwise_qwk = pair_plot_data),
  session = sessionInfo()
)
output_dir <- tempfile("rating-agreement-")
dir.create(output_dir)
saveRDS(analysis, file.path(output_dir, "analysis.rds"))
write.csv(qwk$table, file.path(output_dir, "pairwise-qwk.csv"), row.names = FALSE)
write.csv(comparison$table, file.path(output_dir, "score-comparison.csv"), row.names = FALSE)
write.csv(comparison_cases, file.path(output_dir, "matched-cases.csv"), row.names = FALSE)
if (can_plot) {
  ggplot2::ggsave(file.path(output_dir, "score-counts.png"), figures$score_counts,
                 width = 7, height = 4.5, units = "in", dpi = 150)
  ggplot2::ggsave(file.path(output_dir, "pairwise-qwk.png"), figures$pairwise_qwk,
                 width = 7, height = 3.5, units = "in", dpi = 150)
}

restored <- readRDS(file.path(output_dir, "analysis.rds"))
restored$comparison$table[, c("N", "QWK", "Correlation", "SMD")]
restored$qwk$notes
# With ggplot2 installed: print(restored$figures$pairwise_qwk)

