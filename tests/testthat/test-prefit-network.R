test_that("raw networks distinguish planned, observed and retained scope", {
  skip_if_not_installed("igraph")
  d <- data.frame(P = c("a", "a", "b", "c"),
    R = c("x", "x", "y", "z"), S = c(1, 2, 3, NA), W = c(2, 7, 1, 1))
  net <- mfrm_network_analysis(data = d, person = "P", facets = "R", score = "S", weight = "W")
  expect_equal(net$summary$Components, 4)
  expect_equal(sort(net$edge_metrics$Weight), c(1, 2))
  expect_equal(sum(net$node_metrics$Degree == 0), 2)
  expect_equal(net$source_connectivity$row_usage$Included, c(TRUE, TRUE, TRUE, FALSE))
  expect_null(net$source_connectivity$identification_audit)
  planned <- mfrm_network_analysis(data = d, person = "P", facets = "R", basis = "planned")
  expect_equal(planned$summary$Components, 3)
  selected <- mfrm_network_analysis(data = d, person = "P", facets = "R", score = "S", top_n_subsets = 1)
  expect_false(selected$summary$ScopeComplete)
  expect_true(selected$source_connectivity$settings$is_disjoint)
  expect_equal(selected$summary$SourceSubsets, 4)
  review <- build_mfrm_network_review(data = d, person = "P", facets = "R", score = "S")
  expect_equal(review$network_summary, net$summary)
  expect_s3_class(build_summary_table_bundle(review), "mfrm_summary_table_bundle")
  d$Criterion <- "shared"
  full <- mfrm_network_analysis(data = d, person = "P", facets = c("R", "Criterion"), basis = "planned")
  expect_true(full$summary$Connected)
  expect_match(full$source_connectivity$notes, "confounded")
})

test_that("raw connectivity preserves labels, isolated nodes and empty inputs", {
  skip_if_not_installed("igraph")
  d <- data.frame(P = c("a:b", "c"), check.names = FALSE)
  d[["P:a"]] <- c("b", "d")
  d$Score <- c(NA_real_, NA_real_)
  n <- mfrm_network_analysis(data = d, person = "P", facets = "P:a", score = "Score")
  expect_equal(n$summary$Nodes, 4)
  expect_equal(n$summary$Edges, 0)
  expect_equal(n$summary$Components, 4)
  expect_equal(anyDuplicated(n$node_metrics$Node), 0)
  expect_equal(sum(n$node_metrics$Degree), 0)
  expect_true(all(is.na(n$node_metrics$Closeness)))
  expect_no_warning(z <- subset_connectivity_report(data = d[FALSE, ], person = "P", facets = "P:a", score = "Score"))
  expect_equal(nrow(z$nodes), 0)
  expect_error(subset_connectivity_report(data = d, person = "P", facets = "P", score = "Score"), "distinct")
  expect_error(subset_connectivity_report(data = d, person = "P", facets = "P:a"), "score")
  expect_error(subset_connectivity_report(data = d, person = "P", facets = "P:a", score = "Score", basis = "planned"), "omit")
})

test_that("network display selection keeps calculation and layout identity", {
  skip_if_not_installed("igraph")
  d <- data.frame(P = c("a", "a", "b"), R = c("x", "y", "y"))
  n <- mfrm_network_analysis(data = d, person = "P", facets = "R", basis = "planned")
  before <- serialize(n, NULL)
  p <- plot(n, type = "network", draw = FALSE, display_min_weight = 2)
  expect_equal(nrow(p$data$edges), 0)
  expect_equal(nrow(p$data$source_edges), 3)
  expect_equal(nrow(p$data$nodes), 4)
  expect_error(as_ggplot(p), "network conversion")
  coords <- as.matrix(p$data$coordinates[, c("X", "Y")])
  rownames(coords) <- p$data$coordinates$Node
  q <- plot(n, type = "network", draw = FALSE, layout = coords[4:1, ])
  expect_equal(p$data$coordinates, q$data$coordinates)
  expect_identical(serialize(n, NULL), before)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(n, type = "network", layout = coords,
    graphical = list(vertex.label = NA, vertex.size = 5, edge.curved = .1)))
  expect_error(plot(n, type = "network", draw = FALSE, layout = coords[-1, ]), "one row")
  expect_error(plot(n, type = "network", draw = FALSE, graphical = list(layout = coords)), "other than")
})

test_that("legacy subset diagnostics preserve facet identities with colons", {
  # The retired edge builder remains an independent oracle for the public graph.
  build_subset_connectivity_edges <- function(fit, nodes_tbl, keep_subsets) {
    data <- as.data.frame(fit$prep$data %||% data.frame(), stringsAsFactors = FALSE)
    source_cols <- fit$config$source_columns %||% list()
    person_col <- if ("Person" %in% names(data)) "Person" else as.character(source_cols$person)
    facet_cols <- as.character(fit$prep$facet_names %||% fit$config$facet_names %||% source_cols$facets %||% character(0))
    cols <- c(person_col, facet_cols)
    cols <- cols[nzchar(cols)]
    if (!is.data.frame(data) || nrow(data) == 0L || length(cols) < 2L ||
        !all(cols %in% names(data)) || !is.data.frame(nodes_tbl) ||
        nrow(nodes_tbl) == 0L || !"Node" %in% names(nodes_tbl)) {
      return(data.frame())
    }
    node_lookup <- nodes_tbl[, intersect(c("Node", "Subset", "Facet", "Level"), names(nodes_tbl)), drop = FALSE]
    node_lookup$Node <- as.character(node_lookup$Node)
    if ("Subset" %in% names(node_lookup)) {
      node_lookup$Subset <- suppressWarnings(as.integer(node_lookup$Subset))
      node_lookup <- node_lookup[node_lookup$Subset %in% keep_subsets, , drop = FALSE]
    }
    if (nrow(node_lookup) == 0L) return(data.frame())
  
    long <- do.call(rbind, lapply(cols, function(col) {
      data.frame(
        Row = seq_len(nrow(data)),
        Facet = col,
        Level = as.character(data[[col]]),
        stringsAsFactors = FALSE
      )
    }))
    long <- long[!is.na(long$Level) & nzchar(long$Level), , drop = FALSE]
    long <- merge(
      long,
      node_lookup[, c("Facet", "Level", "Node", "Subset"), drop = FALSE],
      by = c("Facet", "Level"),
      all.x = FALSE,
      all.y = FALSE
    )
    if (nrow(long) == 0L) return(data.frame())
    split_rows <- split(long, long$Row)
    edge_list <- lapply(split_rows, function(row_tbl) {
      row_tbl <- row_tbl[!duplicated(row_tbl$Node), , drop = FALSE]
      row_tbl <- row_tbl[order(row_tbl$Node), , drop = FALSE]
      if (nrow(row_tbl) < 2L) return(NULL)
      cmb <- utils::combn(seq_len(nrow(row_tbl)), 2L)
      from <- row_tbl[cmb[1, ], , drop = FALSE]
      to <- row_tbl[cmb[2, ], , drop = FALSE]
      data.frame(
        Subset = from$Subset,
        From = from$Node,
        To = to$Node,
        FromFacet = from$Facet,
        FromLevel = from$Level,
        ToFacet = to$Facet,
        ToLevel = to$Level,
        stringsAsFactors = FALSE
      )
    })
    edges <- do.call(rbind, edge_list[!vapply(edge_list, is.null, logical(1))])
    if (!is.data.frame(edges) || nrow(edges) == 0L) return(data.frame())
    edges |>
      dplyr::group_by(
        .data$Subset, .data$From, .data$To, .data$FromFacet, .data$FromLevel,
        .data$ToFacet, .data$ToLevel
      ) |>
      dplyr::summarise(Weight = dplyr::n(), .groups = "drop") |>
      dplyr::arrange(.data$Subset, dplyr::desc(.data$Weight), .data$From, .data$To) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  d <- data.frame(Person = c("a:b", "c"), check.names = FALSE)
  d[["Person:a"]] <- c("b", "d")
  s <- calc_subsets(d, names(d))
  expect_equal(nrow(s$nodes), 4)
  expect_equal(anyDuplicated(s$nodes$Node), 0)
  expect_equal(nrow(s$summary), 2)
  expect_equal(s$nodes$Facet, rep(c("Person", "Person:a"), each = 2))
  expect_equal(s$nodes$Level, c("a:b", "c", "b", "d"))
  fit <- list(prep = list(data = d, facet_names = "Person:a"),
    config = list(source_columns = list(person = "OriginalID", facets = "Person:a")))
  edges <- build_subset_connectivity_edges(fit, s$nodes, s$summary$Subset)
  expect_equal(nrow(edges), 2)
  expect_true(all(edges$From != edges$To))
  expect_equal(sort(c(edges$From, edges$To)), sort(s$nodes$Node))
  expect_equal(sort(edges$Weight), c(1, 1))
  raw <- subset_connectivity_report(data = d, person = "Person", facets = "Person:a", basis = "planned")
  expect_equal(raw$edges[, c("From", "To", "Weight")], edges[, c("From", "To", "Weight")])
})
