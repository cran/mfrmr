#' Hierarchical grouping of external features
#'
#' Build an agglomerative hierarchy from weighted Gower dissimilarities, cut
#' it into a chosen number of groups, and inspect the retained dendrogram.
#'
#' @param x For clustering, a table reviewed with [mfrm_features()]. For plotting,
#'   an `mfrm_hierarchical_clusters` result.
#' @inheritParams mfrm_cluster
#' @param linkage `"average"` (default, UPGMA) or `"complete"`. Average linkage
#'   uses the mean dissimilarity over all cross-group pairs; complete linkage
#'   uses their maximum. Neither is selected automatically from the data.
#' @param type `"dendrogram"` (default), `"silhouette"`, or `"profile"`.
#' @inheritParams plot.mfrm_clusters
#' @param labels Whether to display entity IDs on dendrograms and silhouettes.
#'   The default labels at most 50 included entities. Hiding labels does not
#'   remove entities. Profile plots always label groups and levels.
#' @param ... Reserved for future use; additional arguments are rejected.
#' @return An `mfrm_hierarchical_clusters` object inheriting from `mfrm_clusters`.
#'   It retains the `hclust` object in `tree`, ID-aligned `membership` (ID,
#'   Cluster, Silhouette), `cluster_summary`, `profiles`, `feature_data`, and
#'   `settings` including the linkage. Hierarchical clustering does not select
#'   medoids; `medoids` is `NULL` and membership has no `Medoid` column.
#'   `summary()` returns the size/silhouette table. `plot()` invisibly returns
#'   `mfrm_plot_data`; dendrogram data include the tree, leaf order, memberships,
#'   group count, and excluded IDs.
#' @inheritSection mfrm_cluster Reading clustering metrics
#' @details
#' Uses [stats::hclust()] followed by [stats::cutree()] at the requested `k`.
#' Feature types, scaling, weights, missingness handling, silhouette definition,
#' and the 5,000-included-entity limit follow [mfrm_cluster_pam()]. No distance
#' transformation or automatic sampling is performed. This limit is not a
#' memory or speed guarantee.
#'
#' Only average and complete linkage are supported. Ward's minimum-variance
#' criterion requires an appropriate Euclidean geometry; this mixed-feature
#' Gower interface does not define a Euclidean conversion or a Ward analysis.
#' Clustering MFRM bias estimates is a separate methodological question from
#' grouping the external attributes accepted here. Measurement uncertainty is
#' not propagated.
#'
#' The tree is fitted independently of PAM. Its merge heights describe the
#' selected linkage on Gower dissimilarities, not significance or branch support.
#' Tied distances can produce alternative hierarchies and input order can affect
#' their resolution. A cut at `k` uses merge order even when heights tie, so a
#' horizontal height threshold need not uniquely identify that cut.
#'
#' The default plot draws this stored tree, with boxes marking the stored `k`
#' groups. Labels are shown for at most 50 included entities by default; hiding
#' labels does not sample or remove entities. Excluded entities have no leaves
#' but remain in the result and plot data. Silhouette and feature-profile views
#' reuse [plot.mfrm_clusters()]. Plots do not refit or choose groups.
#' [as_ggplot()] converts the stored dendrogram without refitting. The default
#' and `component = "tree"` retain the full tree and group boxes. Leaf order,
#' heights, label settings and excluded IDs are preserved. Dashed boxes differ
#' from the solid tree branches even in monochrome. Box widths and
#' physical text sizes can differ from base graphics. Tied heights retain the
#' merge-order partition; boxes do not imply a unique horizontal height cut.
#' Use `ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)` to remove
#' headings and annotations, and [plot_data()] to inspect the retained source
#' evidence. Silhouette/profile conversions use the same dedicated summary
#' renderers as [plot.mfrm_clusters()].
#'
#' Use [mfrm_cluster_compare()] to compare this partition with PAM or another
#' linkage on the same data. For multiple imputations, use
#' `mfrm_cluster_imputed(..., method = "hierarchical", linkage = "average")`.
#' Each completion retains its own tree; no pooled tree or branch support is
#' estimated. See `vignette("mfrmr-external-features", package = "mfrmr")`.
#' @references Murtagh, F. and Legendre, P. (2014). Ward's Hierarchical
#'   Agglomerative Clustering Method: Which Algorithms Implement Ward's
#'   Criterion? Journal of Classification, 31, 274--295.
#'   \doi{10.1007/s00357-014-9161-z}.
#' @seealso [mfrm_cluster_pam()], [mfrm_cluster_imputed()], [mfrm_cluster_compare()]
#' @examples
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # Fictional rater backgrounds; R7 has unrecorded experience.
#'   raters <- data.frame(Rater = paste0("R", 1:7),
#'     ExperienceYears = c(1, 2, 3, 12, 13, 14, NA),
#'     Specialty = factor(c(rep("Language", 3), rep("Science", 4))))
#'   features <- mfrm_features(raters, "Rater", c("ExperienceYears", "Specialty"))
#'   hierarchy <- mfrm_cluster_hierarchical(features, k = 2, missing = "omit")
#'   plot(hierarchy)
#'   plot(hierarchy, type = "silhouette")
#'   plot(hierarchy, type = "profile", feature = "ExperienceYears")
#'   comparison <- mfrm_cluster_compare(list(
#'     PAM = mfrm_cluster_pam(features, k = 2, missing = "omit"),
#'     Average = hierarchy,
#'     Complete = mfrm_cluster_hierarchical(features, k = 2,
#'       linkage = "complete", missing = "omit")))
#'   comparison$analysis_summary
#'   summary(comparison)
#' }
#' @export
mfrm_cluster_hierarchical <- function(x, k, weights = NULL,
                                      missing = c("error", "omit"),
                                      linkage = c("average", "complete")) {
  linkage <- match.arg(linkage)
  cluster_external_features(x, k, weights, match.arg(missing), linkage)
}

#' @rdname mfrm_cluster_hierarchical
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_hierarchical_clusters <- function(x, type = c("dendrogram", "silhouette", "profile"),
                                            feature = NULL, labels = NULL, draw = TRUE,
                                            preset = "standard", ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  rlang::check_dots_empty()
  type <- match.arg(type)
  if (type != "dendrogram") {
    return(invisible(plot.mfrm_clusters(x, type, feature, labels, draw, preset)))
  }
  if (!is.null(feature)) stop("`feature` is only used for type = 'profile'.", call. = FALSE)
  check_cluster_plot_flags(draw, labels)
  style <- resolve_plot_preset(preset)
  tree <- x$tree
  show_labels <- labels %||% (length(tree$order) <= 50L)
  excluded <- x$membership$ID[is.na(x$membership$Cluster)]
  leaf_order <- tree$labels[tree$order]
  title <- "Hierarchical feature groups"
  subtitle <- sprintf("%s linkage / Gower | Included: %d | Excluded: %d",
    switch(x$settings$linkage, average = "Average", complete = "Complete"),
    x$settings$included, length(excluded))
  out <- new_mfrm_plot_data("cluster_dendrogram", list(tree = tree,
    leaf_order = leaf_order, table = x$membership[match(leaf_order, x$membership$ID), , drop = FALSE],
    k = x$settings$k, excluded_ids = excluded, labels = show_labels,
    title = title, subtitle = subtitle, preset = style$name,
    caption = sprintf("Boxes: %d groups | Heights do not measure branch support", x$settings$k)))
  if (!draw) return(invisible(out))
  apply_plot_preset(style)
  old <- graphics::par(mar = c(5, 5, 5, 2), xpd = FALSE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot(tree, labels = if (show_labels) tree$labels else FALSE,
    hang = -1, main = "", sub = "", xlab = "", ylab = "Gower dissimilarity", axes = FALSE,
    cex = style$axis_cex)
  graphics::axis(2, at = pretty(c(0, tree$height)))
  # Keep outer boxes inside hclust's fixed plotting window, also for small trees.
  ends <- cumsum(rle(out$data$table$Cluster)$lengths)
  usr <- graphics::par("usr")
  padding <- 0.01 * diff(usr[1:2])
  cut_height <- mean(tree$height[length(tree$order) - x$settings$k + 0:1])
  graphics::rect(pmax(c(0, utils::head(ends, -1)) + 0.66, usr[1] + padding),
    usr[3] + 0.005 * diff(usr[3:4]), pmin(ends + 0.33, usr[2] - padding),
    cut_height, border = style$accent_primary)
  invisible(out)
}

# Lay out the saved binary merges; no distances or groups are re-estimated.
.mfrmr_gg_dendrogram <- function(x) {
  .require_mfrmr_ggplot2()
  d <- x$data
  tree <- d$tree
  n <- length(tree$order)
  if (!inherits(tree, "hclust") || n < 2L ||
      !identical(sort(tree$order), seq_len(n)) ||
      !is.matrix(tree$merge) || !identical(dim(tree$merge), c(n - 1L, 2L)) ||
      length(tree$height) != n - 1L || any(!is.finite(tree$height)) ||
      any(tree$height < 0) || length(tree$labels) != n ||
      !identical(d$leaf_order, tree$labels[tree$order]) ||
      !is.data.frame(d$table) || !all(c("ID", "Cluster") %in% names(d$table)) ||
      !identical(d$table$ID, d$leaf_order) || anyNA(d$table$Cluster)) {
    stop("Dendrogram data lack a matching saved tree, leaf order or membership table. Recreate them with plot(hierarchy, draw = FALSE).", call. = FALSE)
  }
  runs <- rle(d$table$Cluster)
  if (length(d$k) != 1L || !is.numeric(d$k) || is.na(d$k) ||
      d$k < 2L || d$k >= n || length(runs$lengths) != d$k ||
      anyDuplicated(runs$values)) {
    stop("Dendrogram groups must match the saved contiguous leaf groups and group count.", call. = FALSE)
  }
  leaf_x <- match(seq_len(n), tree$order)
  node_x <- numeric(n - 1L)
  segments <- matrix(NA_real_, 3L * (n - 1L), 4L,
    dimnames = list(NULL, c("X", "Y", "Xend", "Yend")))
  for (i in seq_len(n - 1L)) {
    children <- tree$merge[i, ]
    if (anyNA(children) || any(children != as.integer(children)) ||
        any(children == 0L | children < -n | children >= i)) {
      stop("The saved dendrogram contains invalid merge indices.", call. = FALSE)
    }
    xx <- yy <- numeric(2)
    for (j in 1:2) {
      child <- children[j]
      xx[j] <- if (child < 0) leaf_x[-child] else node_x[child]
      yy[j] <- if (child < 0) 0 else tree$height[child]
    }
    node_x[i] <- mean(xx)
    segments[3L * (i - 1L) + 1:3, ] <- rbind(
      c(xx[1], yy[1], xx[1], tree$height[i]),
      c(xx[2], yy[2], xx[2], tree$height[i]),
      c(xx[1], tree$height[i], xx[2], tree$height[i]))
  }
  ends <- cumsum(runs$lengths)
  bottom <- -0.04 * if (max(tree$height) > 0) max(tree$height) else 1
  boxes <- data.frame(Left = c(0, utils::head(ends, -1)) + 0.65,
    Right = ends + 0.35, Bottom = bottom,
    Top = mean(tree$height[n - d$k + 0:1]))
  style <- resolve_plot_preset(d$preset %||% "standard")
  p <- ggplot2::ggplot(as.data.frame(segments)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$X, y = .data$Y,
      xend = .data$Xend, yend = .data$Yend), colour = style$foreground, linewidth = 0.45) +
    ggplot2::geom_rect(data = boxes, ggplot2::aes(xmin = .data$Left,
      xmax = .data$Right, ymin = .data$Bottom, ymax = .data$Top),
      fill = NA, colour = style$accent_primary, linewidth = 0.55, linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = if (isTRUE(d$labels)) seq_len(n) else NULL,
      labels = if (isTRUE(d$labels)) d$leaf_order else NULL) +
    ggplot2::scale_y_continuous(breaks = function(limits) {
      ticks <- pretty(c(0, max(tree$height)))
      ticks[ticks >= 0 & ticks <= limits[2]]
    }) +
    .mfrmr_gg_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  caption <- sprintf("Dashed boxes: %d groups | Heights do not measure branch support", d$k)
  if (tree$height[n - d$k] == tree$height[n - d$k + 1L]) caption <- paste0(caption,
    "\nTied heights: boxes follow saved merge-order groups, not a unique height cut.")
  d$caption <- caption
  p <- .mfrmr_gg_labs(p, d, x = NULL, y = "Gower dissimilarity") +
    ggplot2::labs(
      alt = "Hierarchy of external features, using saved merge heights and leaf order. Dashed boxes mark the requested groups, including merge-order cuts at tied heights. Excluded entities have no leaves; heights do not measure significance or branch support.")
  attr(p, "mfrmr_plot_data") <- x
  p
}
