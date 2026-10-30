#' @rdname mfrm_pca
#' @details `plot()` method: Plot numeric feature PCA and optional exploratory groups
#'
#' Review explained variation, original-feature coefficients or entity scores
#' from a saved PCA without refitting it.
#' @param type `"scree"` (default), `"scores"`, or `"loadings"`.
#' @param groups Optional [mfrm_cluster_kmeans()], [mfrm_cluster_pam()] or
#'   [mfrm_cluster_hierarchical()] result used to color a scores plot. Its IDs,
#'   included entities and shared PCA features must match. Labels only annotate
#'   the view; they do not refit PCA or imply separation on every component.
#'   Groups use both colour and point shape, including monochrome output.
#'   Shapes repeat after six groups; inspect the ID-aligned table for crowded views.
#' @inheritParams plot.mfrm_clusters
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object with the exact plotted table,
#'   selected components, group colour/shape encoding, excluded IDs and axis meanings. Scree data retain the
#'   full variance table, including components not used for clustering.
#'   [plot_data()] extracts these values. [as_ggplot()] converts all three views
#'   using their saved axes and encodings; `component = "table"` retains the
#'   complete view. Excluded IDs and transformation metadata remain attached.
#'   Use `ggplot2::labs(title = NULL, subtitle = NULL)` to hide headings in
#'   the returned ggplot. Labels follow the saved `labels` setting; no entities
#'   are sampled when labels are hidden. Physical text and point sizes can
#'   differ between base graphics and ggplot. Converted scores use equal axis
#'   units and equal displayed spans to avoid a narrow panel when the selected
#'   components explain very different amounts of variation.
#' @details For `plot()`: Scores and loading signs are arbitrary. Loadings are eigenvector
#'   coefficients in the centered, weighted and optionally standardized space;
#'   they are not original-unit correlations. A two-component scores view omits
#'   other directions, so apparent overlap or separation is only a projection.
#'   No confidence region, group validity or rater-quality judgment is implied.
#' @seealso [mfrm_pca()], [mfrm_cluster_kmeans()]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot.mfrm_pca", run.dontrun = TRUE)
#' \dontrun{
#' attributes <- data.frame(ID = letters[1:6],
#'   ExperienceYears = c(1, 2, 4, 8, 10, 12),
#'   WorkshopHours = c(8, 16, 12, 24, 16, 32))
#' result <- mfrm_pca(mfrm_features(attributes, "ID", names(attributes)[-1]))
#' plot(result)
#' plot(result, type = "scores")
#' plot(result, type = "loadings", components = 1)
#' plot_data(plot(result, type = "scores", draw = FALSE))$table
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_pca <- function(x, type = c("scree", "scores", "loadings"),
                          components = NULL, groups = NULL, labels = NULL,
                          draw = TRUE, preset = "standard", ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  rlang::check_dots_empty()
  type <- match.arg(type)
  check_cluster_plot_flags(draw, labels)
  style <- resolve_plot_preset(preset)
  if (type != "scores" && !is.null(groups)) stop("`groups` is only used for scores plots.", call. = FALSE)
  excluded <- x$scores$ID[is.na(x$scores$PC1)]
  subtitle <- sprintf("Included: %d | Excluded: %d | Retained components: %d",
    x$settings$included, x$settings$excluded, x$settings$components)
  subtitle <- paste0(subtitle, if (isTRUE(x$settings$scale)) " | Sample-SD scaling" else " | Original-unit scaling")
  if (type == "scree") {
    if (!is.null(components)) stop("Scree plots show every fitted component; do not select components here.", call. = FALSE)
    tab <- x$variance
    title <- "Variation in numeric external features"
    xlab <- "Principal component"
    ylab <- "Explained variance (%)"
  } else {
    if (is.null(components)) components <- if (type == "scores") c(1L, 2L) else 1L
    if (!is.numeric(components) || is.complex(components) ||
        length(components) != (if (type == "scores") 2L else 1L) ||
        any(!is.finite(components)) || any(components != floor(components)) ||
        anyDuplicated(components) || any(components < 1 | components > x$settings$components)) {
      stop("Choose ", if (type == "scores") "two distinct" else "one",
        " retained component number(s).", call. = FALSE)
    }
    axes <- paste0("PC", components)
    if (type == "loadings") {
      tab <- data.frame(Feature = rownames(x$loadings), Loading = x$loadings[, axes])
      title <- paste("Feature coefficients:", axes)
      xlab <- "Coefficient in the transformed feature space"
      ylab <- ""
    } else {
      tab <- x$scores[!x$scores$ID %in% excluded, c("ID", axes), drop = FALSE]
      tab$Cluster <- NA_integer_
      if (!is.null(groups)) {
        if (!inherits(groups, "mfrm_clusters") ||
            !setequal(groups$membership$ID, x$scores$ID) ||
            !setequal(groups$membership$ID[!is.na(groups$membership$Cluster)], tab$ID) ||
            !all(x$feature_data$features %in% groups$feature_data$features)) {
          stop("Groups must retain the PCA IDs, included entities and selected numeric features.", call. = FALSE)
        }
        original <- x$feature_data
        compared <- groups$feature_data
        index <- match(original$row_summary$ID, compared$row_summary$ID)
        same <- isTRUE(all.equal(lapply(original$data[original$features], unname),
          lapply(compared$data[index, original$features, drop = FALSE], unname), tolerance = 0))
        if (!same) stop("Groups and PCA must use the same original feature values.", call. = FALSE)
        tab$Cluster <- groups$membership$Cluster[match(tab$ID, groups$membership$ID)]
      }
      percent <- 100 * x$variance$Proportion[components]
      title <- "Principal component scores"
      xlab <- sprintf("%s (%.1f%% of transformed variance)", axes[1], percent[1])
      ylab <- sprintf("%s (%.1f%% of transformed variance)", axes[2], percent[2])
    }
  }
  encoding <- NULL
  if (type == "scores") {
    levels <- sort(unique(tab$Cluster[!is.na(tab$Cluster)]))
    encoding <- data.frame(Cluster = levels,
      Colour = unname(.plot_series_colors(levels, style$name)),
      Shape = rep(c(16, 17, 15, 18, 1, 2), length.out = length(levels)))
  }
  show_labels <- labels %||% (nrow(tab) <= 50L)
  out <- new_mfrm_plot_data(paste0("feature_pca_", type), list(table = tab,
    components = components, encoding = encoding, transformation = x$transformation,
    weights = x$settings$weights, scale = x$settings$scale,
    excluded_ids = excluded, title = title,
    subtitle = subtitle, xlab = xlab, ylab = ylab, labels = show_labels,
    preset = style$name))
  if (!draw) return(invisible(out))
  apply_plot_preset(style)
  old <- graphics::par(mar = c(6, if (type == "loadings") 10 else 5, 5, 2))
  on.exit(graphics::par(old), add = TRUE)
  if (type == "scree") {
    plot(seq_len(nrow(tab)), 100 * tab$Proportion, type = "l",
      col = style$foreground, xlab = xlab, ylab = ylab, main = "", xaxt = "n")
    axis(1, at = seq_len(nrow(tab)))
    points(seq_len(nrow(tab)), 100 * tab$Proportion,
      pch = ifelse(tab$Retained, 16, 1),
      col = ifelse(tab$Retained, style$accent_primary, "grey65"))
    if (any(!tab$Retained)) legend("topright",
      legend = c("Retained", "Not retained"), pch = c(16, 1),
      col = c(style$accent_primary, "grey65"), bty = "n")
  } else if (type == "loadings") {
    plot(tab$Loading, seq_len(nrow(tab)), xlim = c(-1, 1), yaxt = "n",
      pch = 16, col = style$accent_primary, xlab = xlab, ylab = ylab, main = "")
    abline(v = 0, col = style$grid)
    axis(2, at = seq_len(nrow(tab)), labels = tab$Feature, las = 1)
  } else {
    color <- if (nrow(encoding)) encoding$Colour[match(tab$Cluster, encoding$Cluster)] else style$accent_primary
    shape <- if (nrow(encoding)) encoding$Shape[match(tab$Cluster, encoding$Cluster)] else 16
    plot(tab[[axes[1]]], tab[[axes[2]]], pch = shape, col = color,
      xlab = xlab, ylab = ylab, main = "", asp = 1)
    if (show_labels) text(tab[[axes[1]]], tab[[axes[2]]], labels = tab$ID, pos = 3, cex = 0.7)
    if (nrow(encoding)) legend("topright", legend = paste("Group", encoding$Cluster),
      col = encoding$Colour, pch = encoding$Shape, bty = "n")
  }
  invisible(out)
}

# Use the saved PCA axes and group encoding, never generic numeric-column order.
.mfrmr_gg_feature_pca <- function(x) {
  .require_mfrmr_ggplot2()
  d <- x$data
  tab <- d$table
  type <- sub("^feature_pca_", "", x$name)
  axes <- paste0("PC", d$components)
  required <- switch(type, scree = c("Proportion", "Retained"),
    scores = c("ID", axes, "Cluster"), loadings = c("Feature", "Loading"))
  if (!is.data.frame(tab) || !nrow(tab) || !all(required %in% names(tab)) ||
      (type == "scores" && (length(axes) != 2L || anyDuplicated(axes)))) {
    stop("PCA plot data lack the saved table or selected axes. Recreate them with plot(pca, type = \"", type, "\", draw = FALSE).", call. = FALSE)
  }
  style <- resolve_plot_preset(d$preset %||% "standard")
  if (type == "scree") {
    tab$Position <- seq_len(nrow(tab))
    tab$Retention <- factor(ifelse(tab$Retained, "Retained", "Not retained"),
      levels = c("Retained", "Not retained"))
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$Position, y = 100 * .data$Proportion)) +
      (if (nrow(tab) > 1L) ggplot2::geom_line(colour = style$foreground)) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$Retention, shape = .data$Retention), size = 2.6) +
      ggplot2::scale_colour_manual(values = c(Retained = style$accent_primary, "Not retained" = "grey65"), drop = FALSE) +
      ggplot2::scale_shape_manual(values = c(Retained = 16, "Not retained" = 1), drop = FALSE) +
      ggplot2::scale_x_continuous(breaks = tab$Position) +
      ggplot2::labs(colour = "Components", shape = "Components")
    if (all(tab$Retained)) p <- p + ggplot2::guides(colour = "none", shape = "none")
    alt <- "Explained variance for every fitted external-feature component. Filled points mark retained components; open points mark components not retained. Explained variance does not validate groups."
  } else if (type == "loadings") {
    tab$Position <- seq_len(nrow(tab))
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$Loading, y = .data$Position)) +
      ggplot2::geom_vline(xintercept = 0, colour = style$grid) +
      ggplot2::geom_point(colour = style$accent_primary, size = 2.6) +
      ggplot2::scale_x_continuous(limits = c(-1, 1)) +
      ggplot2::scale_y_continuous(breaks = tab$Position, labels = tab$Feature)
    alt <- "Coefficients for the selected principal component in the transformed external-feature space. Signs are arbitrary; coefficients are not original-unit correlations. The vertical reference marks zero."
  } else {
    encoding <- d$encoding
    grouped <- any(!is.na(tab$Cluster))
    if (grouped && (!is.data.frame(encoding) ||
        !all(c("Cluster", "Colour", "Shape") %in% names(encoding)) ||
        anyDuplicated(encoding$Cluster) || anyNA(tab$Cluster) ||
        !setequal(tab$Cluster, encoding$Cluster))) {
      stop("PCA plot data lack matching saved group colours and shapes. Recreate them with plot(pca, groups = groups, draw = FALSE).", call. = FALSE)
    }
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data[[axes[1]]], y = .data[[axes[2]]]))
    if (grouped) {
      p <- p + ggplot2::geom_point(ggplot2::aes(colour = factor(.data$Cluster),
        shape = factor(.data$Cluster)), size = 2.6) +
        ggplot2::scale_colour_manual(values = stats::setNames(encoding$Colour, encoding$Cluster),
          breaks = as.character(encoding$Cluster)) +
        ggplot2::scale_shape_manual(values = stats::setNames(encoding$Shape, encoding$Cluster),
          breaks = as.character(encoding$Cluster)) +
        ggplot2::labs(colour = "Group", shape = "Group")
    } else p <- p + ggplot2::geom_point(colour = style$accent_primary, shape = 16, size = 2.6)
    if (isTRUE(d$labels)) p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$ID),
      vjust = -0.7, size = 3, colour = style$foreground)
    # Equal units with equal displayed spans avoid narrow panels when one PC
    # explains much less variance. Limits change; stored scores do not.
    ranges <- lapply(tab[axes], range)
    span <- max(vapply(ranges, diff, numeric(1)))
    limits <- lapply(ranges, function(r) mean(r) + c(-0.5, 0.5) * span)
    p <- p + ggplot2::coord_fixed(xlim = limits[[1]], ylim = limits[[2]])
    alt <- paste("Scores on the two selected external-feature principal components with equal axis units.",
      if (grouped) "Saved groups use both colour and shape; shapes repeat after six groups." else "No groups are assigned in this view.",
      "Only included entities are plotted. This projection is not latent ability, rater quality or evidence of group validity.")
  }
  p <- .mfrmr_gg_labs(p + .mfrmr_gg_theme(), d, x = d$xlab, y = d$ylab) +
    ggplot2::labs(alt = alt)
  attr(p, "mfrmr_plot_data") <- x
  p
}
