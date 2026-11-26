.mfrmr_summary_print_sections <- list(print.mfrm_data_description_01 = rbind(c(field = "design_connectivity",
    heading = "\nPerson-facet connectivity\n", cast = "frame", rounded = "no", guard = "present"), c("score_distribution",
    "\nScore distribution\n", "raw", "no", "present")), print.summary.mfrm_data_description_01 = rbind(c(field = "overview",
    heading = "\nOverview\n", cast = "frame", rounded = "yes", guard = "present"), c("missing", "\nMissing by column\n",
    "frame", "yes", "present"), c("structural_missingness", "\nPlanned assignment coverage\n", "frame",
    "yes", "present"), c("design_connectivity", "\nPerson-facet connectivity\n", "frame", "yes", "present"),
    c("linkage_summary", "\nFacet linkage support\n", "frame", "yes", "present"), c("score_distribution",
        "\nScore distribution\n", "frame", "yes", "present"), c("facet_overview", "\nFacet coverage\n",
        "frame", "yes", "present")), print.summary.mfrm_data_description_02 = rbind(c(field = "reporting_map",
    heading = "\nPaper reporting map\n", cast = "frame", rounded = "no", guard = "present")), print.mfrm_anchor_review_01 = rbind(c(field = "facet_summary",
    heading = "\nFacet summary\n", cast = "raw", rounded = "no", guard = "present")))

print_summary_sections <- function(x, sections, digits = 3L) {
    value <- function() do.call("$", list(x, field))
    for (i in seq_len(nrow(sections))) {
        field <- sections[i, "field"]
        guard <- sections[i, "guard"]
        available <- switch(guard, present = !is.null(value()) && nrow(value()) > 0L, fallback = nrow(value() %||%
            data.frame()) > 0L, required = nrow(value()) > 0L)
        if (available) {
            cat(sections[i, "heading"])
            table <- value()
            table <- switch(sections[i, "cast"], raw = table, frame = as.data.frame(table), frame_false = as.data.frame(table,
                stringsAsFactors = FALSE))
            if (sections[i, "rounded"] == "yes")
                table <- round_numeric_df(table, digits = digits)
            print(table, row.names = FALSE)
        }
    }
    invisible(NULL)
}
