.tidymfrm_plan_data <- function(gt) {
  ds <- gt$d_study
  allocation <- vapply(seq_len(nrow(ds)), function(i) {
    paste(paste(gt$facets, unlist(ds[i, gt$facets], use.names = FALSE)), collapse = " x ")
  }, character(1))
  do.call(rbind, lapply(c("G", "Phi"), function(metric) {
    label <- if (metric == "G") "G (relative)" else "Phi (absolute)"
    data.frame(PlanRow = seq_len(nrow(ds)), Ratings = ds$TotalRatingsPerPerson,
               Allocation = allocation, Metric = label, Coefficient = ds[[metric]],
               Target = rep(gt$targets[[metric]], nrow(ds)),
               Reason = ds[[paste0(metric, "Status")]])
  }))
}

.tidymfrm_plot_frontier <- function(gt, palette, base_size) {
  df <- .tidymfrm_plan_data(gt)
  available <- df[is.finite(df$Coefficient), , drop = FALSE]
  available$Attainment <- ifelse(available$Coefficient >= available$Target, "Target met", "Below target")
  frontier <- do.call(rbind, lapply(split(available, available$Metric), function(sub) {
    ordered <- sub[order(sub$Ratings, -sub$Coefficient), , drop = FALSE]
    best <- ordered[!duplicated(ordered$Ratings), , drop = FALSE]
    best[c(TRUE, diff(cummax(best$Coefficient)) > 0), , drop = FALSE]
  }))
  first_met <- available[available$Attainment == "Target met", , drop = FALSE]
  first_met <- first_met[order(first_met$Ratings, -first_met$Coefficient), , drop = FALSE]
  first_met <- first_met[!duplicated(first_met$Metric), , drop = FALSE]

  accent <- .tidymfrm_palette(1L, palette)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Ratings, y = .data$Coefficient)) +
    ggplot2::geom_blank(ggplot2::aes(y = .data$Target)) +
    ggplot2::geom_hline(data = unique(df[c("Metric", "Target")]), ggplot2::aes(yintercept = .data$Target),
                        colour = "#666666", linetype = "dotted") +
    ggplot2::facet_wrap(~ Metric) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    ggplot2::labs(x = "Ratings per person", y = "Projected coefficient") +
    .tidymfrm_theme(base_size)
  if (nrow(available)) {
    p <- p + ggplot2::geom_point(data = available, ggplot2::aes(shape = .data$Attainment),
                                colour = "#666666", size = 2.2) +
      ggplot2::scale_shape_manual(values = c("Below target" = 1, "Target met" = 16), name = NULL)
  }
  empty_panels <- setdiff(df$Metric, available$Metric)
  if (length(empty_panels)) {
    empty <- data.frame(Metric = empty_panels, Ratings = stats::median(df$Ratings), Coefficient = 0.45)
    p <- p + ggplot2::geom_text(data = empty, label = "No projections available", colour = "#666666")
  }
  if (NROW(frontier)) {
    p <- p + ggplot2::geom_step(data = frontier, colour = accent, linewidth = .9, na.rm = TRUE) +
      ggplot2::geom_point(data = frontier, colour = accent, shape = 18, size = 2.6)
  }
  if (nrow(first_met)) {
    p <- p + ggplot2::geom_label(data = first_met, ggplot2::aes(label = .data$Allocation),
                               vjust = ifelse(first_met$Coefficient > .9, 1.6, -0.6),
                               colour = accent, size = 2.8)
  }
  unavailable <- unique(df[!is.finite(df$Coefficient), c("Metric", "Reason")])
  attr(p, "mfrm_notes") <- c(
    "Each circle is a candidate allocation. The step line marks improvements in the best projection among plans at or below that rating count.",
    "Labels identify the best projection at the first evaluated rating count that meets each target. Differences are point projections.",
    paste("Unplotted projections:", nrow(df) - nrow(available)),
    if (nrow(unavailable)) paste(unavailable$Metric, unavailable$Reason, sep = " - "),
    if (length(gt$nesting)) paste(names(gt$nesting), "counts are per", gt$nesting),
    if (length(gt$scores) > 1L) paste("Composite weights", paste(gt$scores, signif(gt$weights, 3), collapse = ", "))
  )
  p
}

.tidymfrm_plot_gstudy <- function(gt, palette, base_size) {
  composite <- length(gt$scores) > 1L
  combined <- data.frame(Score = if (composite) "Composite" else gt$scores[1L],
                         gt$g_study_summary[c("G", "Phi", "GStatus", "PhiStatus")])
  sub <- gt$univariate
  if (composite) {
    sub$Score <- paste("Subscore", sub$Score)
    combined <- rbind(combined, sub[names(combined)])
  }
  df <- rbind(
    data.frame(Score = combined$Score, Metric = "G (relative)", Coefficient = combined$G, Reason = combined$GStatus),
    data.frame(Score = combined$Score, Metric = "Phi (absolute)", Coefficient = combined$Phi, Reason = combined$PhiStatus)
  )
  df$Score <- factor(df$Score, levels = rev(combined$Score))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Coefficient, y = .data$Score,
                                       colour = .data$Metric, shape = .data$Metric)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = .45), size = 3, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = .tidymfrm_palette(2L, palette), name = NULL) +
    ggplot2::scale_shape_manual(values = c(16, 17), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    ggplot2::labs(x = "Coefficient at observed assignment counts", y = NULL) +
    .tidymfrm_theme(base_size)
  unavailable <- df[!is.finite(df$Coefficient), , drop = FALSE]
  attr(p, "mfrm_notes") <- c(
    "Each point uses the observed harmonic assignment counts for the same G-Study. G addresses relative standing and Phi a fixed standard.",
    if (composite) paste("Composite weights", paste(gt$scores, signif(gt$weights, 3), collapse = ", ")),
    if (composite) "Composite coefficients use cross-score covariances. Evaluate each subscore before reporting it separately.",
    if (nrow(unavailable)) paste(unavailable$Score, unavailable$Metric, unavailable$Reason, sep = " - ")
  )
  p
}
