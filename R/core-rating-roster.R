# One matching contract for observed-score agreement coefficients.
mfrm_rating_roster <- function(data, rater, unit, score) {
  if (!is.data.frame(data) || !nrow(data) || anyNA(names(data)) ||
      anyDuplicated(names(data)) || any(!nzchar(names(data)))) stop(
    "`data` must be a nonempty data frame with unique nonblank column names.", call. = FALSE)
  for (nm in c("rater", "unit", "score")) {
    v <- get(nm)
    if (!is.character(v) || !length(v) || anyNA(v) || any(!nzchar(v)) ||
        (nm != "unit" && length(v) != 1L)) stop(
      "`", nm, "` must name ", if (nm == "unit") "one or more columns." else "one column.",
      call. = FALSE)
  }
  selected <- c(rater, unit, score)
  if (anyDuplicated(selected) || !all(selected %in% names(data))) stop(
    "Rater, unit and score must name distinct existing columns.", call. = FALSE)
  ids <- as.data.frame(data[c(rater, unit)])
  for (nm in names(ids)) {
    v <- ids[[nm]]
    if (!is.null(dim(v)) || !(is.factor(v) || (!is.object(v) &&
        (is.character(v) || is.logical(v) || (is.numeric(v) && !is.complex(v))))) ||
        anyNA(v) || any(!nzchar(trimws(as.character(v)))) ||
        (is.numeric(v) && any(!is.finite(v)))) stop(
      "Identifiers must be nonmissing, finite and nonblank: ", nm, ".", call. = FALSE)
  }
  if (anyDuplicated(ids)) stop(
    "Duplicate rater-unit cells are not supported; define the unit explicitly instead of averaging repetitions.",
    call. = FALSE)
  raters <- unique(as.character(ids[[rater]]))
  if (length(raters) != length(unique(ids[[rater]]))) stop(
    "Distinct rater identifiers have identical text labels; supply unique character identifiers.",
    call. = FALSE)
  if (length(raters) < 2L) stop("At least two observed raters are required.", call. = FALSE)
  grouped <- dplyr::group_by(ids[unit], dplyr::across(dplyr::everything()))
  units <- as.data.frame(dplyr::group_keys(grouped))
  rownames(units) <- NULL
  list(ids = ids, raters = raters, unit_id = dplyr::group_indices(grouped), units = units)
}
