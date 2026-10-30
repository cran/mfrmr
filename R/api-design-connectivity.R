.mfrm_design_connectivity <- function(data, person, facets, score = NULL,
                                      basis = "observed", weight = NULL,
                                      top_n_subsets = NULL, min_observations = 0) {
  cols <- c(person, facets)
  if (!is.data.frame(data) || anyDuplicated(names(data)) ||
      !is.character(person) || length(person) != 1L ||
      !is.character(facets) || !length(facets) || anyNA(cols) ||
      anyDuplicated(cols) || any(!nzchar(cols)) || !all(cols %in% names(data))) {
    stop("Supply `data`, one `person` column and distinct `facets` columns.", call. = FALSE)
  }
  if (any(!vapply(data[cols], function(x) is.atomic(x) && is.null(dim(x)), logical(1))))
    stop("Identifier columns must be atomic vectors, not lists or matrices.", call. = FALSE)
  if (!is.numeric(min_observations) || length(min_observations) != 1L ||
      !is.finite(min_observations) || min_observations < 0 ||
      (!is.null(top_n_subsets) && (!is.numeric(top_n_subsets) ||
       length(top_n_subsets) != 1L || !is.finite(top_n_subsets) ||
       top_n_subsets < 1 || top_n_subsets != round(top_n_subsets)))) {
    stop("Use nonnegative `min_observations` and a positive integer `top_n_subsets`.", call. = FALSE)
  }
  ids <- lapply(data[cols], as.character)
  valid_ids <- vapply(ids, function(x) !is.na(x) & nzchar(trimws(x)), logical(nrow(data)))
  valid_ids <- matrix(valid_ids, nrow = nrow(data), ncol = length(cols))
  used <- rowSums(valid_ids) == length(cols)
  reason <- ifelse(used, "Included", "Missing identifier")
  if (basis == "observed") {
    if (!is.character(score) || length(score) != 1L || is.na(score) ||
        !score %in% names(data) || !is.numeric(data[[score]]) ||
        any(is.infinite(data[[score]]))) {
      stop("Observed networks require a numeric `score` column; recode missing markers to NA.", call. = FALSE)
    }
    absent <- is.na(data[[score]])
    reason[used & absent] <- "Missing score"
    used <- used & !absent
    if (!is.null(weight)) {
      if (!is.character(weight) || length(weight) != 1L || is.na(weight) ||
          !weight %in% names(data) || !is.numeric(data[[weight]]) ||
          any(is.infinite(data[[weight]])) || any(data[[weight]] < 0, na.rm = TRUE)) {
        stop("`weight` must name a nonnegative finite numeric column (NA allowed).", call. = FALSE)
      }
      active <- !is.na(data[[weight]]) & data[[weight]] > 0
      reason[used & !active] <- "Missing or zero weight"
      used <- used & active
    }
  } else if (!is.null(score) || !is.null(weight)) {
    stop("Planned networks use assignment rows only; omit `score` and `weight`.", call. = FALSE)
  }
  # Keep known levels even when all their assignment rows are unscored.
  nodes <- do.call(rbind, lapply(seq_along(cols), function(j) {
    levels <- unique(ids[[j]][valid_ids[, j]])
    data.frame(Facet = rep(cols[j], length(levels)), Level = levels)
  }))
  nodes$Node <- make.unique(paste(nodes$Facet, nodes$Level, sep = ":"))
  nodes <- nodes[, c("Node", "Facet", "Level")]
  index <- lapply(seq_along(cols), function(j) {
    at <- which(nodes$Facet == cols[j])
    at[match(ids[[j]], nodes$Level[at])]
  })
  uf <- make_union_find(nodes$Node)
  pair_cols <- utils::combn(seq_along(cols), 2L)
  pairs <- lapply(seq_len(ncol(pair_cols)), function(j) {
    a <- index[[pair_cols[1L, j]]][used]
    b <- index[[pair_cols[2L, j]]][used]
    data.frame(From = nodes$Node[a], To = nodes$Node[b],
      FromFacet = nodes$Facet[a], FromLevel = nodes$Level[a],
      ToFacet = nodes$Facet[b], ToLevel = nodes$Level[b])
  })
  edges <- dplyr::count(do.call(rbind, pairs), .data$From, .data$To,
    .data$FromFacet, .data$FromLevel, .data$ToFacet, .data$ToLevel, name = "Weight")
  reverse <- edges$From > edges$To
  for (suffix in c("", "Facet", "Level")) {
    from <- paste0("From", suffix); to <- paste0("To", suffix)
    old <- edges[[from]][reverse]
    edges[[from]][reverse] <- edges[[to]][reverse]
    edges[[to]][reverse] <- old
  }
  for (i in seq_len(nrow(edges))) uf$union(edges$From[i], edges$To[i])
  roots <- vapply(nodes$Node, uf$find, character(1))
  nodes$Component <- unname(roots)
  nodes$Subset <- match(roots, unique(roots))
  edges$Subset <- nodes$Subset[match(edges$From, nodes$Node)]
  source_n <- length(unique(roots))
  summary <- data.frame(Subset = seq_len(source_n),
    Observations = tabulate(nodes$Subset[index[[1L]][used]], nbins = source_n))
  counts <- table(factor(nodes$Subset, levels = seq_len(source_n)),
                  factor(nodes$Facet, levels = cols))
  # Facet names can coincide with metadata names; listings retain all counts.
  for (col in setdiff(cols, names(summary))) summary[[col]] <- as.integer(counts[, col])
  summary$ObservationPercent <- if (sum(used)) 100 * summary$Observations / sum(used) else rep(NA_real_, source_n)
  summary <- summary[order(-summary$Observations, summary$Subset), , drop = FALSE]
  summary <- summary[summary$Observations >= min_observations, , drop = FALSE]
  if (!is.null(top_n_subsets)) summary <- utils::head(summary, top_n_subsets)
  nodes <- nodes[nodes$Subset %in% summary$Subset, , drop = FALSE]
  edges <- edges[edges$Subset %in% summary$Subset, , drop = FALSE]
  edges <- edges[order(edges$Subset, -edges$Weight, edges$From, edges$To), , drop = FALSE]
  rownames(edges) <- NULL
  rownames(summary) <- NULL
  listing <- nodes |>
    dplyr::group_by(.data$Subset, .data$Facet) |>
    dplyr::summarise(LevelsN = dplyr::n(), Levels = paste(sort(.data$Level), collapse = ", "), .groups = "drop") |>
    dplyr::left_join(summary[, c("Subset", "Observations", "ObservationPercent")], by = "Subset") |>
    dplyr::group_by(.data$Facet) |>
    dplyr::mutate(Ruler = vapply(round(20 * .data$LevelsN / max(1, .data$LevelsN)),
      function(n) paste0("[", strrep("=", n), strrep(".", 20 - n), "]"), character(1))) |>
    dplyr::ungroup()
  list(summary = summary, nodes = nodes, edges = as.data.frame(edges),
    listing = as.data.frame(listing),
    row_usage = data.frame(Row = seq_len(nrow(data)), Included = used, Reason = reason),
    settings = list(top_n_subsets = top_n_subsets %||% NA_integer_,
      min_observations = min_observations, source_subsets = source_n,
      retained_subsets = nrow(summary), is_disjoint = source_n > 1L,
      basis = basis, person = person, facets = facets, score = score, weight = weight,
      input_rows = nrow(data), included_rows = sum(used)))
}
