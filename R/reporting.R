to_float <- function(value) {
    out <- suppressWarnings(as.numeric(value))
    ifelse(length(out) == 0, NA_real_, out)
}

fmt_count <- function(value) {
    val <- to_float(value)
    if (!is.finite(val))
        return("NA")
    if (abs(val - round(val)) < 1e-06)
        return(as.character(as.integer(round(val))))
    sprintf("%.0f", val)
}
