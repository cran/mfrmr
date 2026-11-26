.mfrm_gstudy_interval_context <- function(ci, g) {
    attr(ci, "source_summary") <- g[c("design", "components", "mean_products", "degrees_of_freedom")]
    needs_data <- any(grepl("REML|^Residual-contrast", ci$Method))
    attr(ci, "source_requires_data") <- needs_data
    if (needs_data)
        attr(ci, "source_data") <- g$data
    ci
}
