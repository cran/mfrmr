.mfrmr_empty_table <- function(types, tbl = FALSE) {
    x <- as.data.frame(lapply(types, vector, length = 0L), check.names = FALSE, stringsAsFactors = FALSE)
    if (tbl)
        tibble::as_tibble(x)
    else x
}


`%||%` <- function(x, y) if (is.null(x)) y else x
