.mfrm_network_display <- function(x, draw = TRUE, layout = NULL,
                                   display_min_weight = 0, graphical = list(),
                                   main = NULL, preset = .mfrm_default_plot_preset(),
                                   palette = NULL, label_angle = 45, metric = NULL, top_n = 20) {
  style <- resolve_plot_preset(preset)
  design <- inherits(x, c("mfrm_network_analysis", "mfrm_subset_connectivity"))
  halo <- inherits(x, "mfrm_halo_network")
  source <- if (inherits(x, "mfrm_network_analysis")) x$source_connectivity else x
  nodes <- as.data.frame(source$nodes)
  edges <- as.data.frame(source$edges)
  if (!design) {
    validate_descriptive_bundle(x, if (halo) "mfrm_halo_network" else "mfrm_rater_network")
    nodes <- as.data.frame(x$node_metrics)
    edges <- as.data.frame(x$edge_metrics)
  }
  id_col <- if (design || halo) "Node" else "Rater"
  if (!id_col %in% names(nodes) || !nrow(nodes) ||
      !all(c("From", "To", "Weight") %in% names(edges))) {
    stop("No network node/edge tables are available.", call. = FALSE)
  }
  ids <- as.character(nodes[[id_col]])
  if (anyNA(ids) || anyDuplicated(ids)) stop("Network IDs must be unique and nonmissing.", call. = FALSE)
  if (!is.numeric(display_min_weight) || length(display_min_weight) != 1L ||
      !is.finite(display_min_weight) || display_min_weight < 0) {
    stop("`display_min_weight` must be a nonnegative number.", call. = FALSE)
  }
  if (!is.list(graphical) || (length(graphical) &&
      (is.null(names(graphical)) || any(!nzchar(names(graphical))) ||
       anyDuplicated(names(graphical)) || any(names(graphical) %in% c("x", "layout"))))) {
    stop("`graphical` must be a named list of plot.igraph arguments other than x/layout.", call. = FALSE)
  }
  shown <- is.finite(edges$Weight) & abs(edges$Weight) >= display_min_weight
  displayed <- edges[shown, , drop = FALSE]
  directed <- !design && !halo && isTRUE(x$settings$directed)
  need_graph <- isTRUE(draw) || is.function(layout)
  if (need_graph) {
    if (!requireNamespace("igraph", quietly = TRUE)) stop("Drawing requires igraph.", call. = FALSE)
    graph <- igraph::graph_from_data_frame(displayed[, c("From", "To", "Weight")],
      directed = directed, vertices = data.frame(name = ids))
  }
  if (is.function(layout)) layout <- layout(graph)
  if (is.null(layout)) {
    angle <- seq(0, 2 * pi, length.out = length(ids) + 1L)[seq_along(ids)]
    layout <- cbind(cos(angle), sin(angle))
    rownames(layout) <- ids
  }
  if (!is.matrix(layout) || !is.numeric(layout) ||
      !identical(dim(layout), c(length(ids), 2L)) || any(!is.finite(layout))) {
    stop("`layout` must be a finite two-column matrix with one row per node.", call. = FALSE)
  }
  if (!is.null(rownames(layout))) {
    if (anyDuplicated(rownames(layout)) || !setequal(rownames(layout), ids))
      stop("Named layout rows must match all network node IDs.", call. = FALSE)
    layout <- layout[match(ids, rownames(layout)), , drop = FALSE]
  }
  group <- if (design) nodes$Facet else if (halo) nodes$Rater else rep("Rater", nrow(nodes))
  groups <- unique(as.character(group))
  colours <- stats::setNames(grDevices::hcl.colors(max(3L, length(groups)), "Dark 3")[seq_along(groups)], groups)
  if (identical(style$name, "monochrome")) colours[] <- grDevices::gray.colors(length(groups), start = .3, end = .8)
  if (!is.null(palette)) {
    if (!is.null(names(palette)) && all(groups %in% names(palette))) {
      colours <- palette[groups]
    } else if (!design && !halo) {
      pal <- resolve_palette(palette, defaults = c(agreement = style$accent_primary))
      colours[] <- pal["agreement"]
    }
  }
  width <- if (nrow(displayed)) 0.5 + 4 * abs(displayed$Weight) /
    max(1e-12, abs(displayed$Weight)) else numeric()
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    signed <- if ("SignedWeight" %in% names(displayed)) displayed$SignedWeight else displayed$Weight
    defaults <- list(x = graph, layout = layout, main = "",
      vertex.label = if (length(ids) > 80L) NA_character_ else if (design) nodes$Level else ids,
      vertex.color = unname(colours[as.character(group)]),
      vertex.size = if (design) 6 else 12, vertex.label.cex = 0.85,
      vertex.label.family = "sans", xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
      vertex.label.color = style$foreground, vertex.frame.color = style$foreground,
      edge.width = width, edge.color = ifelse(signed < 0, "#2166AC", "#777777"),
      edge.lty = ifelse(signed < 0, 2L, 1L),
      edge.arrow.size = if (directed) 0.4 else 0, edge.curved = if (directed) 0.15 else 0)
    defaults[names(graphical)] <- graphical
    do.call(.mfrm_plot_igraph, defaults)
  }
  notes <- c(x$caveats$Message,
    if (halo) c(
      "These are observed-score profile correlations, not model residual correlations or severity estimates. Graph centrality is not rater centrality (middle-category use).",
      "Missing or filtered edges do not mean zero correlation or absence of halo."),
    "Display filtering changes this picture only; computed metrics and source edges are unchanged.",
    "Node position and size do not estimate ability or severity unless explicitly mapped to a declared variable.")
  invisible(new_mfrm_plot_data(if (design) "subset_connectivity" else if (halo) "halo_network" else "rater_network",
    list(plot = "network", nodes = nodes, edges = displayed, source_edges = edges,
      coordinates = data.frame(Node = ids, X = layout[, 1], Y = layout[, 2]),
      directed = directed, mode = x$settings$mode, display_min_weight = display_min_weight,
      settings = x$settings, summary = x$summary, caveats = x$caveats,
      title = main %||% "Rating network", notes = notes,
      legend = new_plot_legend(groups, "node_group", "node_color", unname(colours)),
      reference_lines = new_reference_lines(), preset = style$name)))
}
