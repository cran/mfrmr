#' @rdname mfrm_d_study
#' @export
mfrm_d_study <- function(x, design_grid = NULL, ...) {
  UseMethod("mfrm_d_study")
}

#' @rdname mfrm_d_study
#' @export
mfrm_d_study.default <- function(x, design_grid = NULL, ...) {
  stop("`x` must be a G-study from mfrm_g_study(), a main-effects result from mfrm_generalizability(), or an mfrm_fit object.", call. = FALSE)
}

#' @rdname mfrm_d_study_main_effects
#' @export
mfrm_d_study.mfrm_fit <- function(x, design_grid = NULL,
    object_facet = "Person", random_facets = NULL,
    residual_scaling = c("highest_order", "single_condition", "none", "sensitivity"), ...) {
  residual_scaling <- match.arg(residual_scaling)
  g <- mfrm_generalizability(x, object_facet = object_facet,
    random_facets = random_facets, ...)
  mfrm_d_study(g, design_grid, residual_scaling = residual_scaling)
}

#' @rdname mfrm_d_study_main_effects
#' @param object A main-effects D-study result.
#' @export
summary.mfrm_d_study <- function(object, ...) {
  rlang::check_dots_empty()
  validate_gtheory_output(object)
  as.data.frame(object)
}

mfrm_d_study_count_columns <- function(grid, facets, columns) {
  if (!is.data.frame(grid) && (!is.list(grid) || is.object(grid)))
    stop("`design_grid` must be a data frame or named list of planned counts.", call. = FALSE)
  grid <- as.data.frame(grid, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(grid)) stop("`design_grid` must contain at least one design row.", call. = FALSE)
  if (anyNA(names(grid)) || anyDuplicated(names(grid)) || any(!nzchar(trimws(names(grid)))))
    stop("`design_grid` must have unique, nonmissing, nonblank column names.", call. = FALSE)
  candidates <- lapply(seq_along(facets), function(i)
    intersect(c(facets[i], paste0("n_", facets[i]), columns[i]), names(grid)))
  if (any(lengths(candidates) == 0L))
    stop("`design_grid` is missing count column(s) for random facet(s): ",
      paste(facets[lengths(candidates) == 0L], collapse = ", "),
      ". Use the facet name, `n_<facet>`, or its saved count-column name.", call. = FALSE)
  if (ncol(grid) != length(facets))
    stop("`design_grid` must contain exactly one count column per random facet; remove extra columns or duplicate aliases.", call. = FALSE)
  selected <- rep(NA_character_, length(facets))
  # Resolve aliases jointly: a literal n_Rater facet must not reuse Rater's count.
  for (i in seq_along(facets)) {
    pending <- which(is.na(selected))
    available <- lapply(candidates[pending], setdiff, selected[!is.na(selected)])
    forced <- which(lengths(available) == 1L)
    if (any(lengths(available) == 0L) || !length(forced))
      stop("`design_grid` count columns must correspond one-to-one to the random facets; use the complete set of facet names or `n_<facet>` names.", call. = FALSE)
    j <- forced[1L]
    selected[pending[j]] <- available[[j]]
  }
  grid <- grid[selected]
  names(grid) <- columns
  grid
}
