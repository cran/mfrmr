
print.mfrm_bundle <- function(x, ...) {
    cat(sprintf("<%s>\n", class(x)[1L]))
    if (is.data.frame(x$table))
        print(utils::head(x$table, 20L), row.names = FALSE)
    else if (is.data.frame(x$stacked))
        print(utils::head(x$stacked, 20L), row.names = FALSE)
    invisible(x)
}


summary.mfrm_bundle <- function(object, ...) object

print_bundle_section <- function(title, table, digits = 3, round_numeric = TRUE) {
    if (is.null(table) || !is.data.frame(table) || nrow(table) == 0)
        return(invisible(NULL))
    cat("\n", title, "\n", sep = "")
    if (isTRUE(round_numeric)) {
        print(round_numeric_df(as.data.frame(table), digits = digits), row.names = FALSE)
    }
    else {
        print(as.data.frame(table), row.names = FALSE)
    }
    invisible(NULL)
}

summary.mfrm_bias <- function(object, digits = 3, top_n = 10, p_cut = 0.05, ...) {
    if (!is.list(object) || is.null(object$table) || nrow(object$table) == 0) {
        stop("`object` must be non-empty output from estimate_bias().")
    }
    digits <- max(0L, as.integer(digits))
    top_n <- max(1L, as.integer(top_n))
    p_cut <- max(0, min(1, as.numeric(p_cut[1L])))
    bt <- object$table
    b_vals <- if ("BiasSize" %in% names(bt))
        bt$BiasSize
    else bt$Estimate
    p_vals <- bt$pValue
    p_adj <- if ("pAdj" %in% names(bt))
        bt$pAdj
    else stats::p.adjust(p_vals, method = "holm")
    pair_lbl <- if ("Pair" %in% names(bt))
        paste(unique(bt$Pair), collapse = ", ")
    else paste(unique(bt$Interaction), collapse = ", ")
    overview <- tibble::tibble(FacetPair = pair_lbl, InteractionOrder = 2L, InteractionMode = object$method %||%
        "local", Cells = nrow(bt), Points = sum(is.finite(b_vals)), ScreenAvailable = sum(is.finite(p_vals)),
        ScreenUnavailable = sum(!is.finite(p_vals)), MeanAbsBias = mean(abs(b_vals), na.rm = TRUE), MaxAbsBias = max(abs(b_vals),
            na.rm = TRUE), ScreenPositive = if (all(is.na(p_vals)))
            NA_integer_
        else sum(p_vals <= p_cut, na.rm = TRUE), ScreeningCut = p_cut, BonferroniScreenPositive = if (all(is.na(p_vals)))
            NA_integer_
        else sum(stats::p.adjust(p_vals, method = "bonferroni") <= p_cut, na.rm = TRUE), HolmScreenPositive = if (all(is.na(p_adj)))
            NA_integer_
        else sum(p_adj <= p_cut, na.rm = TRUE))
    ord_idx <- order(-abs(b_vals))
    top_rows <- tibble::as_tibble(bt[utils::head(ord_idx, top_n), , drop = FALSE])
    out <- list(overview = overview, summary = object$summary %||% tibble::tibble(), pairwise = utils::head(object$pairwise %||%
        data.frame(), top_n), lrt = object$lrt, chi_sq = tibble::tibble(), final_iteration = tibble::tibble(),
        top_rows = top_rows, notes = c(if (!is.null(object$test_label)) paste("Pairwise contrast test:", object$test_label) else "Joint interaction parameterization.",
            if (!is.null(object$readiness)) paste("Numerical inference:", object$readiness$status, paste(object$readiness$reason_codes,
                collapse = "; "))), digits = digits)
    class(out) <- "summary.mfrm_bias"
    return(out)
}

print.summary.mfrm_bias <- function(x, ...) {
    digits <- as.integer(x$digits %||% 3L)
    if (!is.finite(digits))
        digits <- 3L
    cat("Many-Facet Measurement Bias Summary\n")
    if (!is.null(x$overview) && nrow(x$overview) > 0) {
        ov <- round_numeric_df(as.data.frame(x$overview), digits = digits)[1, , drop = FALSE]
        cat(sprintf("  Interaction facets: %s | Cells: %s\n", ov$FacetPair, ov$Cells))
        if ("InteractionOrder" %in% names(ov) && "InteractionMode" %in% names(ov)) {
            cat(sprintf("  Order: %s | Mode: %s\n", ov$InteractionOrder, ov$InteractionMode))
        }
        cat(sprintf("  Mean |Bias|: %s | Max |Bias|: %s | Screen-positive (p <= %.3f): %s\n", ov$MeanAbsBias,
            ov$MaxAbsBias, as.numeric(ov$ScreeningCut), ov$ScreenPositive))
        if ("LRScreenPositive" %in% names(ov)) {
            cat(sprintf("  GPCM profile-LR screen-positive (p <= %.3f): %s\n", as.numeric(ov$ScreeningCut),
                ov$LRScreenPositive))
        }
        if (all(c("BonferroniScreenPositive", "HolmScreenPositive") %in% names(ov))) {
            cat(sprintf("  Bonferroni screen-positive: %s | Holm screen-positive: %s (cut = %.3f, m = %s)\n",
                ov$BonferroniScreenPositive, ov$HolmScreenPositive, as.numeric(ov$ScreeningCut), ov$Cells))
        }
    }
    print_bundle_section("Largest cell effects", x$top_rows, digits)
    print_bundle_section("Pairwise contrasts", x$pairwise, digits)
    print_bundle_section("Joint interaction likelihood-ratio test", x$lrt, digits)
    if (length(x$notes) > 0) {
        cat("\nNotes\n")
        for (line in x$notes) cat(" - ", line, "\n", sep = "")
    }
    invisible(x)
}

round_numeric_df <- function(df, digits = 3L) {
    if (!is.data.frame(df) || nrow(df) == 0)
        return(df)
    out <- df
    numeric_cols <- vapply(out, is.numeric, logical(1))
    out[numeric_cols] <- lapply(out[numeric_cols], round, digits = digits)
    out
}

clean_summary_lines <- function(lines, max_n = Inf) {
    if (length(lines) == 0)
        return(character(0))
    lines <- trimws(as.character(lines))
    lines <- lines[!is.na(lines) & nzchar(lines)]
    lines <- unique(lines)
    if (is.finite(max_n) && length(lines) > max_n) {
        lines <- lines[seq_len(max_n)]
    }
    lines
}

mfrm_console_width <- function() {
    width <- getOption("width", 80L)
    if (length(width) != 1L || !is.numeric(width) || is.na(width) || !is.finite(width)) {
        width <- 80L
    }
    min(100L, max(50L, as.integer(width)))
}


print_bullet_section <- function(title, lines, prefix = " - ") {
    lines <- clean_summary_lines(lines)
    if (length(lines) == 0)
        return(invisible(NULL))
    cat("\n", title, "\n", sep = "")
    console_width <- mfrm_console_width()
    content_width <- max(30L, console_width - nchar(prefix))
    continuation <- paste(rep(" ", nchar(prefix)), collapse = "")
    for (line in lines) {
        wrapped <- strwrap(line, width = content_width)
        if (length(wrapped) == 0L)
            next
        cat(prefix, wrapped[1], "\n", sep = "")
        if (length(wrapped) > 1L) {
            for (part in wrapped[-1L]) cat(continuation, part, "\n", sep = "")
        }
    }
    invisible(NULL)
}

summary_note_display_table <- function(notes, preparation = FALSE) {
    notes <- as.data.frame(notes %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(notes) == 0L)
        return(notes)
    keep <- intersect(if (preparation)
        c("Stage", "Severity", "Condition", "Count", "Affected", "RecommendedAction")
    else c("Area", "Severity", "Condition", "Categories", "Message", "RecommendedAction"), names(notes))
    if (length(keep) == 0L)
        return(notes)
    notes[, keep, drop = FALSE]
}

print_note_section <- function(notes, title = "Caveats", preparation = FALSE) {
    notes <- summary_note_display_table(notes, preparation = preparation)
    if (nrow(notes) == 0L)
        return(invisible(NULL))
    cat("\n", title, "\n", sep = "")
    print(notes, row.names = FALSE)
    invisible(NULL)
}

print.mfrm_fair_average <- function(x, ..., digits = 3, n = 20L) {
    if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 1 || n != as.integer(n))
        stop("n must be a positive integer.", call. = FALSE)
    cat("Fair Average scores\n")
    if (!is.null(x$reference_scores)) {
        print_bundle_section("Common reference plan", utils::head(x$reference_scores, n), digits)
        cat("Intervals, when present, condition on the calibration and reference plan.\n")
    }
    print_bundle_section("Observed-design Fair(M) and measures", utils::head(x$stacked, n), digits)
    invisible(x)
}

print.mfrm_bias <- function(x, ...) {
    print(summary(x, ...))
    invisible(x)
}
