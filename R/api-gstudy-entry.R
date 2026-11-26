mfrm_d_study <- function(x, design_grid = NULL, ...) {
    UseMethod("mfrm_d_study")
}

mfrm_d_study.default <- function(x, design_grid = NULL, ...) {
    stop("`x` must be a G-study from mfrm_g_study(), or a current mfrm_fit object.",
        call. = FALSE)
}

mfrm_d_study_count_columns <- function(grid, facets, columns) {
    if (!is.data.frame(grid) && (!is.list(grid) || is.object(grid)))
        stop("`design_grid` must be a data frame or named list of planned counts.", call. = FALSE)
    grid <- as.data.frame(grid, stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(grid))
        stop("`design_grid` must contain at least one design row.", call. = FALSE)
    if (anyNA(names(grid)) || anyDuplicated(names(grid)) || any(!nzchar(trimws(names(grid)))))
        stop("`design_grid` must have unique, nonmissing, nonblank column names.", call. = FALSE)
    candidates <- lapply(seq_along(facets), function(i) intersect(c(facets[i], paste0("n_", facets[i]),
        columns[i]), names(grid)))
    if (any(lengths(candidates) == 0L))
        stop("`design_grid` is missing count column(s) for random facet(s): ", paste(facets[lengths(candidates) ==
            0L], collapse = ", "), ". Use the facet name, `n_<facet>`, or its saved count-column name.",
            call. = FALSE)
    if (ncol(grid) != length(facets))
        stop("`design_grid` must contain exactly one count column per random facet; remove extra columns or duplicate aliases.",
            call. = FALSE)
    selected <- rep(NA_character_, length(facets))
    for (i in seq_along(facets)) {
        pending <- which(is.na(selected))
        available <- lapply(candidates[pending], setdiff, selected[!is.na(selected)])
        forced <- which(lengths(available) == 1L)
        if (any(lengths(available) == 0L) || !length(forced))
            stop("`design_grid` count columns must correspond one-to-one to the random facets; use the complete set of facet names or `n_<facet>` names.",
                call. = FALSE)
        j <- forced[1L]
        selected[pending[j]] <- available[[j]]
    }
    grid <- grid[selected]
    names(grid) <- columns
    grid
}
