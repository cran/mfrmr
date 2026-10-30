#' Plot exploratory groups and imputation sensitivity
#'
#' Inspect within-sample separation, individual feature profiles, or pairwise
#' co-membership across imputations using existing clustering results.
#'
#' @param x An object returned by [mfrm_cluster_pam()], [mfrm_cluster_kmeans()] or
#'   [mfrm_cluster_imputed()], as appropriate.
#' @param type For a single partition, `"silhouette"` or `"profile"`.
#' @param feature A single selected feature name, required for `type = "profile"`.
#'   Numeric features show means and medians in their original units;
#'   categorical features show within-group proportions in the original level
#'   order (including unused factor levels).
#' @param ids For imputation heatmaps, an optional character vector of distinct
#'   entity IDs in the desired display order. Selection affects only the view,
#'   not clustering or the denominator. By default all IDs appear in input order.
#' @param labels Whether to display entity IDs on silhouettes or heatmaps.
#'   The default displays them for at most 50 entities. No entities are sampled
#'   when labels are hidden. Profile plots always label groups and levels.
#' @param draw Draw the plot when `TRUE`; `FALSE` only returns plotted values.
#' @param preset Plot style: `"standard"`, `"publication"`, `"compact"`, or
#'   `"monochrome"`.
#' @param ... Reserved for future use; additional arguments are rejected.
#' @return Invisibly, an `mfrm_plot_data` object. Its `data` contains the plotted
#'   table or matrix, title, subtitle, legend, and excluded IDs. Heatmaps also
#'   retain the displayed IDs and the number of imputations. Use [plot_data()]
#'   to extract this payload for custom graphics. Profile plots also retain
#'   `observations`: original-unit entity values and group/exclusion status,
#'   suitable for custom distribution plots; see [plot_data()].
#'   [as_ggplot()] converts imputation co-membership heatmaps using the saved
#'   matrix, ID order, label choice, colours and imputation count. The default
#'   and `component = "matrix"` retain the complete view and its fixed
#'   zero-to-one scale. Unavailable cells have both grey fill and crosses,
#'   distinguishing them from zero even in monochrome. Metadata, including
#'   excluded IDs, remain available with [plot_data()]. No values are
#'   recomputed or renormalized when IDs are selected. Use
#'   `ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)` to hide
#'   annotations while retaining the source data.
#'   Silhouette conversion retains negative widths, saved order and the overall
#'   mean reference, with group labels independent of colour. Numeric profiles
#'   retain original-unit means/medians and counts; circles and triangles are
#'   offset vertically to show coincident values without adding intervals.
#'   Categorical profiles retain the original category order, unused levels,
#'   group counts and a fixed zero-to-one proportion scale. Default and
#'   `component = "table"` preserve the complete selected view; categorical
#'   profiles also accept `component = "matrix"`. These conversions do not
#'   recluster, select groups, or estimate uncertainty.
#' @details
#' No model or clustering is refitted. Silhouette widths describe separation
#' in the fitted sample, not stability or probabilities. The dashed line is
#' the overall mean silhouette; excluded entities have no silhouette.
#'
#' Feature profiles describe one partition; numeric summaries have no confidence
#' intervals. For an imputed result, inspect a completion with
#' `plot(x$analyses[[1]], type = "profile", feature = "ExperienceYears")`.
#' Cluster numbers must not be averaged across imputations.
#'
#' Imputation heatmaps show the fraction of all supplied imputations in which
#' each pair belongs to the same group, on a fixed zero-to-one scale. Grey cells
#' are unavailable pairs involving excluded entities, not zero co-membership.
#' These fractions describe sensitivity to imputations under fixed settings,
#' not membership probabilities, sampling stability, or a consensus partition.
#' Rows and columns follow input order or explicit `ids`; no hierarchical
#' clustering is performed. PAM is nonhierarchical and these plots do not
#' provide a dendrogram. Use [mfrm_cluster_hierarchical()] for a separate
#' hierarchical analysis and its dendrogram.
#' @seealso [mfrm_cluster_pam()], [mfrm_cluster_imputed()], [mfrm_cluster_compare()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_clusters", run.dontrun = TRUE)
#' \dontrun{
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   raters <- data.frame(Rater = paste0("R", 1:6),
#'     ExperienceYears = c(1, 2, 3, 12, 13, 14),
#'     Specialty = factor(rep(c("Language", "Science"), each = 3)))
#'   groups <- mfrm_cluster_pam(mfrm_features(raters, "Rater",
#'     c("ExperienceYears", "Specialty")), k = 2)
#'   plot(groups)
#'   plot(groups, type = "profile", feature = "ExperienceYears")
#'   plot(groups, type = "profile", feature = "Specialty")
#'   values <- plot_data(plot(groups, draw = FALSE))
#'   values$table
#' }
#' }
#' @name plot.mfrm_clusters
NULL

#' @rdname plot.mfrm_clusters
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_clusters <- function(x, type = c("silhouette", "profile"),
                               feature = NULL, labels = NULL, draw = TRUE,
                               preset = "standard", ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  rlang::check_dots_empty()
  type <- match.arg(type)
  check_cluster_plot_flags(draw, labels)
  style <- resolve_plot_preset(preset)
  excluded <- x$membership$ID[is.na(x$membership$Cluster)]
  method <- if (is.null(x$settings$linkage)) x$settings$method else
    paste(x$settings$linkage, "linkage")
  subtitle <- sprintf("%s / %s | Included: %d | Excluded: %d",
    method, x$settings$distance, x$settings$included, length(excluded))
  if (!is.null(x$settings$space)) {
    subtitle <- paste0(subtitle, " | ", x$settings$space,
      if (!is.null(x$settings$components)) paste0(": ", x$settings$components) else "",
      if (isTRUE(x$settings$scale)) " | Sample-SD scaling" else " | Original-unit scaling")
  }
  if (type == "silhouette") {
    if (!is.null(feature)) stop("`feature` is only used for type = 'profile'.", call. = FALSE)
    tab <- x$membership[!is.na(x$membership$Cluster), , drop = FALSE]
    if (!nrow(tab) || any(!is.finite(tab$Silhouette))) {
      stop("Silhouettes are unavailable; fit with silhouette = TRUE or use a profile plot.", call. = FALSE)
    }
    tab <- tab[order(tab$Cluster, -tab$Silhouette, tab$ID), , drop = FALSE]
    show_labels <- labels %||% (nrow(tab) <= 50L)
    average <- mean(tab$Silhouette)
    out <- new_mfrm_plot_data("cluster_silhouette", list(table = tab,
      cluster_summary = x$cluster_summary, excluded_ids = excluded,
      labels = show_labels, title = "Silhouette widths", subtitle = subtitle,
      caption = sprintf("Dashed: overall mean %.3f | In-sample separation", average),
      reference_lines = new_reference_lines("x", average, "Overall mean", "dashed", "mean"),
      preset = style$name))
    if (!draw) return(invisible(out))
    apply_plot_preset(style)
    old <- graphics::par(mar = c(6, 7, 5, if (show_labels) 6 else 2))
    on.exit(graphics::par(old), add = TRUE)
    pos <- seq_len(nrow(tab))
    plot(NA, xlim = c(-1, 1), ylim = c(nrow(tab) + 0.5, 0.5),
      xlab = "Silhouette width", ylab = "", yaxt = "n", main = "")
    abline(v = 0, col = style$grid)
    rect(0, pos - 0.4, tab$Silhouette, pos + 0.4, border = NA,
      col = ifelse(tab$Cluster %% 2L == 1L, style$accent_primary, style$accent_secondary))
    abline(v = average, lty = 2, col = style$foreground)
    centers <- tapply(pos, tab$Cluster, mean)
    axis(2, at = centers, labels = paste("Group", names(centers)), las = 1)
    if (show_labels) axis(4, at = pos, labels = tab$ID, las = 1, tick = FALSE)
    return(invisible(out))
  }
  if (!is.character(feature) || length(feature) != 1L || is.na(feature) ||
      !feature %in% x$feature_data$features) {
    stop("`feature` must name one selected feature for type = 'profile'.", call. = FALSE)
  }
  numeric <- is.numeric(x$feature_data$data[[feature]])
  tab <- x$profiles[[if (numeric) "numeric" else "categorical"]]
  tab <- tab[tab$Feature == feature, , drop = FALSE]
  title <- paste("Feature profile:", feature)
  groups <- x$cluster_summary$Cluster
  group_labels <- sprintf("Group %d (n = %d)", groups, x$cluster_summary$N)
  payload <- list(table = tab, excluded_ids = excluded, feature = feature,
    title = title, subtitle = subtitle, preset = style$name)
  original <- x$feature_data
  ids <- original$row_summary$ID
  membership <- x$membership$Cluster[match(ids, x$membership$ID)]
  reasons <- original$missing[original$missing$Feature == feature, , drop = FALSE]
  payload$observations <- data.frame(ID = ids, Feature = feature,
    Value = original$data[[feature]], Cluster = membership,
    Included = !is.na(membership),
    MissingReason = reasons$Reason[match(ids, reasons$ID)], check.names = FALSE)
  if (numeric) {
    payload$legend <- new_plot_legend(c("Mean", "Median"), rep("summary", 2),
      rep("pch", 2), c("16", "17"))
  } else {
    levels <- unique(tab$Level)
    mat <- matrix(0, length(groups), length(levels), dimnames = list(group_labels, levels))
    mat[cbind(match(tab$Cluster, groups), match(tab$Level, levels))] <- tab$Proportion
    payload$matrix <- mat
    payload$legend <- cluster_proportion_legend(style)
  }
  payload$caption <- if (numeric) "Mean: circles | Median: triangles. Descriptive summaries without confidence intervals." else "Within-group proportion"
  out <- new_mfrm_plot_data("cluster_profile", payload)
  if (!draw) return(invisible(out))
  apply_plot_preset(style)
  if (!numeric) {
    draw_cluster_proportions(mat, title, subtitle, style, TRUE, "Within-group proportion")
  } else {
    old <- graphics::par(mar = c(6, 9, 5, 2))
    on.exit(graphics::par(old), add = TRUE)
    pos <- seq_len(nrow(tab))
    plot(tab$Mean, pos, xlim = range(c(tab$Mean, tab$Median)),
      ylim = c(nrow(tab) + 0.5, 0.5), pch = 16, col = style$accent_primary,
      xlab = paste(feature, "(original units)"), ylab = "", yaxt = "n", main = "")
    points(tab$Median, pos, pch = 17, col = style$accent_secondary)
    axis(2, at = pos, labels = group_labels, las = 1)
  }
  invisible(out)
}

#' @rdname plot.mfrm_clusters
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_imputed_clusters <- function(x, ids = NULL, labels = NULL, draw = TRUE,
                                       preset = "standard", ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  rlang::check_dots_empty()
  check_cluster_plot_flags(draw, labels)
  style <- resolve_plot_preset(preset)
  all_ids <- rownames(x$co_membership)
  if (is.null(ids)) ids <- all_ids
  if (!is.character(ids) || !length(ids) || anyNA(ids) || anyDuplicated(ids) ||
      !all(ids %in% all_ids)) {
    stop("`ids` must be distinct entity IDs present in the result.", call. = FALSE)
  }
  mat <- x$co_membership[ids, ids, drop = FALSE]
  excluded <- all_ids[is.na(diag(x$co_membership))]
  show_labels <- labels %||% (length(ids) <= 50L)
  method <- x$settings$method %||% "PAM"
  linkage <- x$settings$linkage
  title <- paste("Co-membership:", if (is.null(linkage)) method else paste(linkage, "linkage"))
  subtitle <- sprintf("Imputations: %d | Displayed: %d of %d | Excluded in view: %d",
    x$settings$imputations, length(ids), length(all_ids), sum(ids %in% excluded))
  out <- new_mfrm_plot_data("cluster_co_membership", list(matrix = mat,
    ids = ids, excluded_ids = excluded, imputations = x$settings$imputations,
    labels = show_labels, method = method, linkage = linkage, title = title, subtitle = subtitle,
    caption = "Fraction of all supplied imputations | Grey: unavailable",
    legend = cluster_proportion_legend(style), preset = style$name))
  if (!draw) return(invisible(out))
  apply_plot_preset(style)
  draw_cluster_proportions(mat, title, subtitle, style, show_labels,
    "Fraction of imputations | Grey: unavailable")
  invisible(out)
}

check_cluster_plot_flags <- function(draw, labels) {
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw) ||
      (!is.null(labels) && (!is.logical(labels) || length(labels) != 1L || is.na(labels)))) {
    stop("`draw` and non-NULL `labels` must each be TRUE or FALSE.", call. = FALSE)
  }
}

cluster_proportion_legend <- function(style) {
  new_plot_legend(c(as.character(seq(0, 1, 0.25)), "Unavailable"),
    c(rep("proportion", 5), "missing"), rep("fill", 6),
    c(grDevices::colorRampPalette(c("white", style$accent_primary))(5), "grey70"))
}

draw_cluster_proportions <- function(mat, title, subtitle, style, labels, caption) {
  old <- graphics::par(mar = c(if (labels) 7 else 4, if (labels) 9 else 4, 5, 7))
  on.exit(graphics::par(old), add = TRUE)
  nr <- nrow(mat)
  nc <- ncol(mat)
  plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5),
    xaxs = "i", yaxs = "i", xlab = "", ylab = "", axes = FALSE, main = "")
  rect(0.5, 0.5, nc + 0.5, nr + 0.5, col = "grey70", border = NA)
  image(seq(0.5, nc + 0.5), seq(0.5, nr + 0.5),
    t(mat[nr:1, , drop = FALSE]), add = TRUE, useRaster = TRUE,
    col = grDevices::colorRampPalette(c("white", style$accent_primary))(100),
    breaks = seq(0, 1, length.out = 101))
  if (labels) {
    axis(1, at = seq_len(nc), labels = colnames(mat), las = 2, tick = FALSE)
    axis(2, at = seq_len(nr), labels = rev(rownames(mat)), las = 1, tick = FALSE)
  }
  graphics::box()
  key <- cluster_proportion_legend(style)
  legend(nc + 0.5 + 0.03 * nc, nr + 0.5, legend = key$label,
    fill = key$value, bty = "n", xpd = NA, cex = 0.8)
}

# Preserve the complete co-membership view, including unavailable pairs.
.mfrmr_gg_co_membership <- function(x) {
  .require_mfrmr_ggplot2()
  d <- x$data
  mat <- d$matrix
  ids <- d$ids
  if (!is.matrix(mat) || !is.numeric(mat) || is.complex(mat) || !length(mat) ||
      nrow(mat) != ncol(mat) || !identical(rownames(mat), ids) ||
      !identical(colnames(mat), ids) || anyNA(ids) || anyDuplicated(ids) ||
      any(!is.na(mat) & (!is.finite(mat) | mat < 0 | mat > 1))) {
    stop("Co-membership data require a square zero-to-one matrix with matching saved IDs; unavailable cells must be NA. Recreate with plot(groups, draw = FALSE).", call. = FALSE)
  }
  if (length(d$imputations) != 1L || !is.numeric(d$imputations) ||
      !is.finite(d$imputations) || d$imputations < 1 || d$imputations != floor(d$imputations)) {
    stop("Co-membership data lack the saved number of imputations.", call. = FALSE)
  }
  key <- d$legend
  if (!is.data.frame(key) || !all(c("role", "value") %in% names(key)) || anyNA(key$role) ||
      sum(key$role == "proportion") != 5L || sum(key$role == "missing") != 1L ||
      anyNA(key$value)) {
    stop("Co-membership data lack the saved proportion and unavailable-cell colours. Recreate with plot(groups, draw = FALSE).", call. = FALSE)
  }
  n <- nrow(mat)
  tab <- expand.grid(Row = seq_len(n), Column = seq_len(n))
  tab$Y <- n + 1L - tab$Row
  tab$Fraction <- as.vector(mat)
  unavailable <- tab[is.na(tab$Fraction), , drop = FALSE]
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$Column, y = .data$Y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$Fraction)) +
    ggplot2::scale_fill_gradientn(colours = key$value[key$role == "proportion"],
      values = seq(0, 1, length.out = 5), limits = c(0, 1), breaks = seq(0, 1, 0.25),
      na.value = key$value[key$role == "missing"], name = "Fraction of\nimputations") +
    ggplot2::scale_x_continuous(breaks = if (isTRUE(d$labels)) seq_len(n) else NULL,
      labels = if (isTRUE(d$labels)) ids else NULL, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = if (isTRUE(d$labels)) seq_len(n) else NULL,
      labels = if (isTRUE(d$labels)) rev(ids) else NULL, expand = c(0, 0)) +
    ggplot2::coord_equal() + .mfrmr_gg_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (nrow(unavailable)) {
    p <- p + ggplot2::geom_point(data = unavailable,
      ggplot2::aes(shape = "Unavailable"), colour = "black", size = 2) +
      ggplot2::scale_shape_manual(values = c(Unavailable = 4), name = NULL)
  }
  d$caption <- paste("Fractions use all supplied imputations; not membership probabilities or sampling stability.",
    if (nrow(unavailable)) "Crosses mark unavailable pairs; zero means never in the same group.")
  p <- .mfrmr_gg_labs(p, d, x = NULL, y = NULL) +
    ggplot2::labs(alt = "Pairwise co-membership fractions on a fixed zero-to-one scale, in the saved ID order from top to bottom and left to right. Crosses mark unavailable cells. Values describe sensitivity across the supplied imputations, not membership probabilities, sampling stability or a consensus partition.")
  attr(p, "mfrmr_plot_data") <- x
  p
}

.mfrmr_gg_cluster_summary <- function(x) {
  .require_mfrmr_ggplot2()
  d <- x$data
  tab <- d$table
  silhouette <- identical(x$name, "cluster_silhouette")
  categorical <- !silhouette && !is.null(d$matrix)
  required <- if (silhouette) c("ID", "Cluster", "Silhouette") else if (categorical)
    c("Cluster", "Feature", "Level", "Proportion") else c("Cluster", "Mean", "Median", "N")
  if (!is.data.frame(tab) || !nrow(tab) || !all(required %in% names(tab))) {
    stop("Cluster plot data lack the saved summary table. Recreate them with plot(groups, draw = FALSE).", call. = FALSE)
  }
  style <- resolve_plot_preset(d$preset %||% "standard")
  if (silhouette) {
    refs <- normalize_reference_lines(d$reference_lines)
    if (!any(refs$axis == "x" & refs$role == "mean" & is.finite(refs$value), na.rm = TRUE)) {
      stop("Silhouette data lack the saved overall-mean reference. Recreate with plot(groups, draw = FALSE).", call. = FALSE)
    }
    if (any(!is.finite(tab$Silhouette) | abs(tab$Silhouette) > 1) || anyNA(tab$Cluster)) {
      stop("Saved silhouette widths must be finite and between -1 and 1.", call. = FALSE)
    }
    tab$Y <- nrow(tab) + 1L - seq_len(nrow(tab))
    tab$Colour <- ifelse(tab$Cluster %% 2L == 1L, style$accent_primary, style$accent_secondary)
    centers <- tapply(tab$Y, tab$Cluster, mean)
    p <- ggplot2::ggplot(tab) +
      ggplot2::geom_rect(ggplot2::aes(xmin = pmin(0, .data$Silhouette),
        xmax = pmax(0, .data$Silhouette), ymin = .data$Y - 0.4,
        ymax = .data$Y + 0.4, fill = .data$Colour)) +
      ggplot2::scale_fill_identity() +
      ggplot2::geom_vline(xintercept = 0, colour = style$grid) +
      ggplot2::scale_x_continuous(limits = c(-1, 1)) +
      ggplot2::scale_y_continuous(breaks = unname(centers), labels = paste("Group", names(centers)),
        sec.axis = ggplot2::dup_axis(name = NULL,
          breaks = if (isTRUE(d$labels)) rev(tab$Y) else NULL,
          labels = if (isTRUE(d$labels)) rev(tab$ID) else NULL))
    p <- .mfrmr_gg_add_references(p, d$reference_lines)
    d$caption <- "Dashed line: saved overall mean silhouette. In-sample separation, not sampling stability or rater quality."
    xlab <- "Silhouette width"
    alt <- "Silhouette widths in the saved within-group order on a fixed minus-one-to-one scale. Negative widths extend left of zero. Group labels do not rely on colour; the dashed line shows the saved overall mean. Excluded entities have no widths."
  } else if (categorical) {
    mat <- d$matrix
    key <- d$legend
    if (!is.matrix(mat) || !is.numeric(mat) || !length(mat) ||
        any(!is.finite(mat) | mat < 0 | mat > 1) ||
        is.null(rownames(mat)) || is.null(colnames(mat)) ||
        !is.data.frame(key) || !all(c("role", "value") %in% names(key)) ||
        anyNA(key$role) || sum(key$role == "proportion") != 5L) {
      stop("Categorical profiles need the saved zero-to-one matrix, group/category labels and proportion colours.", call. = FALSE)
    }
    grid <- expand.grid(Row = seq_len(nrow(mat)), Column = seq_len(ncol(mat)))
    grid$Y <- nrow(mat) + 1L - grid$Row
    grid$Proportion <- as.vector(mat)
    p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data$Column, y = .data$Y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$Proportion), colour = style$grid, linewidth = 0.35) +
      ggplot2::scale_fill_gradientn(colours = key$value[key$role == "proportion"],
        values = seq(0, 1, length.out = 5), limits = c(0, 1), breaks = seq(0, 1, .25),
        name = "Within-group\nproportion") +
      ggplot2::scale_x_continuous(breaks = seq_len(ncol(mat)), labels = colnames(mat)) +
      ggplot2::scale_y_continuous(breaks = seq_len(nrow(mat)), labels = rev(rownames(mat))) +
      ggplot2::coord_equal()
    d$caption <- "Within-group proportions in original category order, including unused levels. Descriptive summaries, not effects or confidence intervals."
    xlab <- d$feature
    alt <- "Categorical feature profiles by saved group and original category order, including unused levels. Cell values are within-group proportions on a fixed zero-to-one scale. Group labels retain sample sizes; these summaries do not show uncertainty or causal effects."
  } else {
    if (any(!is.finite(tab$Mean)) || any(!is.finite(tab$Median)) || anyNA(tab$Cluster) ||
        any(!is.finite(tab$N) | tab$N <= 0)) {
      stop("Numeric profiles need finite saved means, medians and positive counts.", call. = FALSE)
    }
    positions <- nrow(tab) + 1L - seq_len(nrow(tab))
    points <- rbind(data.frame(Y = positions + .06, Value = tab$Mean, Statistic = "Mean"),
      data.frame(Y = positions - .06, Value = tab$Median, Statistic = "Median"))
    p <- ggplot2::ggplot(points, ggplot2::aes(x = .data$Value, y = .data$Y,
      colour = .data$Statistic, shape = .data$Statistic)) +
      ggplot2::geom_point(size = 2.8) +
      ggplot2::scale_colour_manual(values = c(Mean = style$accent_primary, Median = style$accent_secondary)) +
      ggplot2::scale_shape_manual(values = c(Mean = 16, Median = 17)) +
      ggplot2::scale_y_continuous(breaks = rev(positions), limits = c(.5, nrow(tab) + .5),
        labels = rev(sprintf("Group %d (n = %d)", tab$Cluster, tab$N))) +
      ggplot2::labs(colour = "Summary", shape = "Summary")
    d$caption <- "Means and medians in original units, offset vertically for visibility. Descriptive summaries without confidence intervals."
    xlab <- paste(d$feature, "(original units)")
    alt <- "Numeric feature means and medians by group, in original units. Circles and triangles distinguish summaries even in monochrome; slight vertical offsets separate coincident values. Group labels retain sample sizes. No intervals or rater-quality judgments are shown."
  }
  p <- .mfrmr_gg_labs(p + .mfrmr_gg_theme(), d, x = xlab, y = NULL) + ggplot2::labs(alt = alt)
  if (categorical) p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  attr(p, "mfrmr_plot_data") <- x
  p
}
