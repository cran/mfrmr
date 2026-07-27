calc_subsets_reference <- function(obs_df, facet_cols) {
  if (is.null(obs_df) || nrow(obs_df) == 0L || length(facet_cols) == 0L) {
    return(list(summary = tibble::tibble(), nodes = tibble::tibble()))
  }
  df <- obs_df[, facet_cols, drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) == 0L) {
    return(list(summary = tibble::tibble(), nodes = tibble::tibble()))
  }

  nodes <- unique(unlist(lapply(facet_cols, function(facet) {
    paste0(facet, ":", as.character(df[[facet]]))
  })))
  if (length(nodes) == 0L) {
    return(list(summary = tibble::tibble(), nodes = tibble::tibble()))
  }
  uf <- mfrmr:::make_union_find(nodes)

  for (i in seq_len(nrow(df))) {
    row_vals <- unlist(df[i, facet_cols, drop = FALSE], use.names = FALSE)
    row_nodes <- paste0(facet_cols, ":", as.character(row_vals))
    row_nodes <- row_nodes[!is.na(row_nodes)]
    if (length(row_nodes) < 2L) next
    base <- row_nodes[1L]
    for (node in row_nodes[-1L]) {
      uf$union(base, node)
    }
  }

  comp_ids <- vapply(nodes, uf$find, character(1))
  comp_levels <- unique(comp_ids)
  comp_index <- stats::setNames(seq_along(comp_levels), comp_levels)

  node_tbl <- tibble::tibble(
    Node = nodes,
    Component = comp_ids,
    Subset = comp_index[comp_ids],
    Facet = stringr::str_extract(nodes, "^[^:]+"),
    Level = stringr::str_replace(nodes, "^[^:]+:", "")
  )

  facet_counts <- node_tbl |>
    dplyr::group_by(.data$Subset, .data$Facet) |>
    dplyr::summarize(Levels = dplyr::n_distinct(.data$Level), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = "Facet",
      values_from = "Levels",
      values_fill = 0
    )

  row_subset <- vapply(seq_len(nrow(df)), function(i) {
    node <- paste0(
      facet_cols[1L], ":", as.character(df[[facet_cols[1L]]][i])
    )
    comp_index[uf$find(node)]
  }, integer(1))
  obs_counts <- tibble::tibble(Subset = row_subset) |>
    dplyr::count(.data$Subset, name = "Observations")

  summary_tbl <- facet_counts |>
    dplyr::left_join(obs_counts, by = "Subset") |>
    dplyr::mutate(Observations = tidyr::replace_na(.data$Observations, 0)) |>
    dplyr::arrange(dplyr::desc(.data$Observations))

  list(summary = summary_tbl, nodes = node_tbl)
}

expect_calc_subsets_parity <- function(data, facets) {
  expected <- calc_subsets_reference(data, facets)
  observed <- mfrmr:::calc_subsets(data, facets)
  testthat::expect_identical(observed$nodes, expected$nodes)
  testthat::expect_identical(observed$summary, expected$summary)
}

test_that("calc_subsets optimized edge traversal preserves reference output", {
  data <- data.frame(
    Person = c("P2", "P1", "P2", "P3", "P4", "P3", "P2", NA),
    Rater = c("R2", "R1", "R2", "R3", "R4", "R3", "R1", "R1"),
    Criterion = c("C1", "C1", "C2", "C3", "C4", "C3", "C1", "C1"),
    stringsAsFactors = FALSE
  )
  data <- rbind(data, data[c(1L, 3L, 3L, 6L), , drop = FALSE])

  expect_calc_subsets_parity(data, c("Person", "Rater", "Criterion"))
  expect_calc_subsets_parity(data, c("Rater", "Person", "Criterion"))
})

test_that("calc_subsets parity covers disconnected and single-facet designs", {
  disconnected <- data.frame(
    Person = c("P1", "P2", "P3", "P4", "P1", "P3"),
    Rater = c("R1", "R1", "R2", "R2", "R1", "R2"),
    Criterion = c("C1", "C1", "C2", "C2", "C1", "C2"),
    stringsAsFactors = FALSE
  )

  expect_calc_subsets_parity(
    disconnected,
    c("Person", "Rater", "Criterion")
  )
  expect_calc_subsets_parity(disconnected, "Person")
})

test_that("calc_subsets parity holds for repeated crossed observations", {
  crossed <- expand.grid(
    Person = sprintf("P%03d", seq_len(40L)),
    Rater = sprintf("R%02d", seq_len(6L)),
    Criterion = sprintf("C%02d", seq_len(5L)),
    stringsAsFactors = FALSE
  )
  crossed <- crossed[c(seq_len(nrow(crossed)), 1:300, 1:300), , drop = FALSE]

  expect_calc_subsets_parity(
    crossed,
    c("Person", "Rater", "Criterion")
  )
})
