mfrm_validate_cut_scores <- function(cut_scores) {
    if (is.null(cut_scores))
        return(NULL)
    if (!is.numeric(cut_scores) || is.complex(cut_scores) || is.object(cut_scores) || !is.null(dim(cut_scores)) ||
        !length(cut_scores) || any(!is.finite(cut_scores)) || any(diff(cut_scores) <= 0))
        stop("`cut_scores` must be a finite, strictly increasing numeric vector.", call. = FALSE)
    labels <- names(cut_scores)
    if (!is.null(labels) && (anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels)))
        stop("Names of `cut_scores` must be distinct, nonempty labels.", call. = FALSE)
    cut_scores
}
