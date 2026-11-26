summary.mfrm_multivariate_d_study <- function(object, intervals = FALSE, parm = c("G", "Phi", "RelativeSEM",
    "AbsoluteSEM"), ...) {
    if (!identical(.mfrm_mvds_policy(object), "raw"))
        .mfrm_mvds_replay(object)
    if (isFALSE(intervals)) {
        rlang::check_dots_empty()
        if (!missing(parm))
            stop("`parm` applies only when `intervals = TRUE`.", call. = FALSE)
        return(object$coefficients)
    }
    if (isTRUE(intervals)) {
        if (!is.character(parm) || !length(parm) || !is.null(dim(parm)) || anyNA(parm) || anyDuplicated(parm) ||
            any(!parm %in% c("G", "Phi", "RelativeSEM", "AbsoluteSEM")))
            stop("`parm` must contain distinct G, Phi, RelativeSEM or AbsoluteSEM metrics.", call. = FALSE)
        intervals <- setNames(lapply(parm, function(metric) confint.mfrm_multivariate_d_study(object,
            parm = metric, ...)), parm)
    }
    else {
        rlang::check_dots_empty()
        if (!missing(parm))
            stop("Saved intervals retain their own metrics and settings.", call. = FALSE)
        if (inherits(intervals, "mfrm_dstudy_intervals"))
            intervals <- list(intervals)
    }
    if (!is.list(intervals) || !length(intervals) || !all(vapply(intervals, inherits, logical(1), "mfrm_dstudy_intervals")))
        stop("`intervals` must be FALSE, TRUE, a confint(d) result or a list of those results.", call. = FALSE)
    for (ci in intervals) .mfrm_mvds_summary_check(ci, object)
    labels <- names(intervals)
    if (is.null(labels))
        labels <- rep("", length(intervals))
    unnamed <- is.na(labels) | !nzchar(labels)
    labels[unnamed] <- vapply(intervals[unnamed], function(ci) ci$Metric[1L], character(1))
    names(intervals) <- make.unique(labels)
    columns <- unique(unlist(lapply(intervals, names), use.names = FALSE))
    rows <- lapply(seq_along(intervals), function(i) {
        ci <- intervals[[i]]
        row <- as.data.frame(setNames(lapply(ci, identity), names(ci)))
        for (name in setdiff(columns, names(row))) row[[name]] <- NA
        cbind(Interval = names(intervals)[i], row[columns])
    })
    table <- do.call(rbind, rows)
    rownames(table) <- NULL
    notes <- unique(c("Future plans are complete balanced ratings for every person and score; source counts and missingness remain separate. Weights are used as supplied.",
        "G/Phi are dimensionless; SEM is in score units. For plan differences, positive G/Phi and negative SEM favor the comparison plan, not necessarily its cost-effectiveness.",
        "No new simultaneous or selection-adjusted guarantee follows from combining results. Preserve method notes, unavailable states and any separate REML or residual-contrast IntervalCenter.",
        unlist(lapply(intervals, function(ci) attr(ci, "notes")), use.names = FALSE)))
    structure(list(table = table, intervals = intervals, source_design = object$gstudy$design, source_rows = object$gstudy$data_usage$counts,
        design_grid = object$design_grid, weights = object$weights, component_diagnostics = object$component_diagnostics,
        notes = notes), class = "summary.mfrm_multivariate_d_study")
}

.mfrm_mvds_interval_source <- function(ci, object) {
    .mfrm_mvds_interval_policy(ci, object)
    required <- c("Scenario", "Kind", "Score", "Metric", "Estimate", "PointStatus", "Lower", "Upper",
        "Level", "MethodCode", "Method", "Status", "Available")
    if (!all(required %in% names(ci)) || !nrow(ci))
        stop("Incomplete interval result; recreate it with confint(d).", call. = FALSE)
    source <- object$gstudy[c("design", "components", "mean_products", "degrees_of_freedom")]
    if (is.null(attr(ci, "source_summary")) || !isTRUE(all.equal(attr(ci, "source_summary"), source,
        tolerance = 0)))
        stop("Interval source summaries do not match; recreate intervals with confint(d) for this D-study.",
            call. = FALSE)
    if (!isTRUE(all.equal(attr(ci, "source_rows"), object$gstudy$data_usage$counts, tolerance = 0)))
        stop("Interval source-row accounting does not match this D-study.", call. = FALSE)
    needs_data <- any(c("reml_profile", "residual_contrast", "spectral_f") %in% ci$MethodCode) || ("delta" %in%
        ci$MethodCode && identical(object$gstudy$design$method, "Multivariate MINQUE(0)"))
    if (needs_data && (is.null(attr(ci, "source_data")) || !isTRUE(all.equal(attr(ci, "source_data"),
        object$gstudy$data, tolerance = 0))))
        stop("Interval source ratings do not match this D-study.", call. = FALSE)
}

.mfrm_mvds_summary_check <- function(ci, object) {
    .mfrm_mvds_interval_source(ci, object)
    metric <- unique(ci$Metric)
    kind <- unique(ci$Kind)
    score <- unique(ci$Score)
    base <- sub(" difference$", "", metric)
    if (length(metric) != 1L || anyNA(metric) || !base %in% c("G", "Phi", "RelativeSEM", "AbsoluteSEM") ||
        length(kind) != 1L || anyNA(kind) || !kind %in% c("Score", "Composite") || length(score) != 1L ||
        anyNA(score) || !isTRUE(all.equal(attr(ci, "design_grid"), object$design_grid, tolerance = 0)) ||
        !identical(attr(ci, "assumption"), "normal"))
        stop("Intervals must retain their metric, score, assumptions and matching plans.", call. = FALSE)
    tab <- object$coefficients
    tab <- tab[tab$Kind == kind & tab$Score == score, , drop = FALSE]
    if (nrow(tab) != nrow(object$design_grid))
        stop("Interval score/composite is absent from this D-study.", call. = FALSE)
    w <- if (kind == "Score")
        setNames(as.numeric(object$gstudy$design$scores == score), object$gstudy$design$scores)
    else if (is.matrix(object$weights))
        object$weights[, score]
    else object$weights
    if (!isTRUE(all.equal(attr(ci, "weights"), w, tolerance = 0)))
        stop("Interval weights do not match the selected score/composite.", call. = FALSE)
    ref <- attr(ci, "reference")$scenario
    plans <- seq_len(nrow(object$design_grid))
    difference <- metric != base
    if (difference) {
        if (length(ref) != 1L || is.na(ref) || !ref %in% plans || !all(c("Reference", "Value", "ReferenceValue") %in%
            names(ci)) || !all(ci$Reference == ref))
            stop("Difference intervals must retain a valid reference plan.", call. = FALSE)
        plans <- setdiff(plans, ref)
    }
    else if (!is.null(ref))
        stop("Individual-plan intervals must not have a reference plan.", call. = FALSE)
    if (anyNA(ci[c("Scenario", "MethodCode", "Level")]) || any(!is.finite(ci$Level) | ci$Level <= 0 |
        ci$Level >= 1) || !all(ci$Scenario %in% plans) || any(vapply(split(ci$Scenario, interaction(ci$MethodCode,
        ci$Level, drop = TRUE)), function(s) anyDuplicated(s) > 0L || !setequal(s, plans), logical(1))))
        stop("Retain one interval row per plan and method/level, including unavailable rows.", call. = FALSE)
    counts <- names(object$design_grid)
    if (!all(counts %in% names(ci)) || !isTRUE(all.equal(unname(as.matrix(ci[counts])), unname(as.matrix(object$design_grid[ci$Scenario,
        , drop = FALSE])), tolerance = 0)))
        stop("Interval facet counts differ from the D-study plans.", call. = FALSE)
    idx <- match(ci$Scenario, tab$Scenario)
    values <- tab[[base]][idx]
    if (difference && (!isTRUE(all.equal(ci$Value, values, tolerance = 0)) || !isTRUE(all.equal(ci$ReferenceValue,
        rep(tab[[base]][ref], nrow(ci)), tolerance = 0))))
        stop("Interval comparison values differ from the stored D-study.", call. = FALSE)
    expected <- if (difference)
        values - tab[[base]][ref]
    else values
    status <- if (difference)
        ifelse(is.finite(expected), "Available", "Point difference unavailable")
    else tab[[paste0(base, "Status")]][idx]
    if (!isTRUE(all.equal(ci$Estimate, expected, tolerance = 0)) || !identical(ci$PointStatus, status))
        stop("Interval point estimates/statuses differ from the stored D-study.", call. = FALSE)
    available <- is.finite(ci$Lower) & is.finite(ci$Upper)
    lower <- ci$Lower[available]
    upper <- ci$Upper[available]
    bounded <- ci$MethodCode[available] != "delta" | !base %in% c("G", "Phi")
    if (!identical(ci$Available, available) || any(lower > upper) || (!difference && any(bounded & lower <
        0)) || (base %in% c("G", "Phi") && any(bounded & (upper > 1 | lower < -as.integer(difference)))))
        stop("Inconsistent interval endpoints or availability; recreate intervals with confint(d).",
            call. = FALSE)
    invisible(NULL)
}

print.summary.mfrm_multivariate_d_study <- function(x, digits = NULL, row.names = FALSE, ...) {
    cat("Assessment plans with source-study uncertainty\n")
    .mfrm_mvgt_print_source(x$source_design, x$source_rows)
    cat("Future complete balanced plans; counts apply to every person and score.\n")
    .mfrm_mvgt_print_future_counts(x$source_design)
    cat("Intervals assume independent normal random effects; the assumption is not tested.\n")
    cat("Nominal confidence level(s):", paste(unique(x$table$Level), collapse = ", "), "\n")
    if (!is.null(x$source_design$score_convention))
        cat(x$source_design$score_convention, "\n")
    if (!is.null(x$weights)) {
        cat("Composite weights (not normalized):\n")
        print(x$weights, digits = digits)
    }
    columns <- c("Scenario", names(x$design_grid), "Kind", "Score", "Metric", "Reference", "Estimate",
        "ComponentTreatment", "IntervalEstimate", "SE", "IntervalSE", "IntervalCenter", "Lower", "Upper",
        "Level", "MethodCode")
    print(x$table[intersect(columns, names(x$table))], digits = digits, row.names = row.names, ...)
    cat("Methods and qualifications:\n")
    print(unique(x$table[intersect(c("MethodCode", "Method", "Validation"), names(x$table))]), row.names = row.names,
        ...)
    review <- is.na(x$table$PointStatus) | x$table$PointStatus != "Available" | is.na(x$table$Available) |
        !x$table$Available
    if (!is.null(x$table$FormulaComplete))
        review <- review | (!is.na(x$table$FormulaComplete) & !x$table$FormulaComplete)
    if (any(review)) {
        cat("Rows needing review (including unavailable or conservatively completed intervals):\n")
        columns <- c("Interval", "Scenario", "Score", "Metric", "Level", "PointStatus", "Status", "Available",
            "FormulaComplete")
        print(x$table[review, intersect(columns, names(x$table)), drop = FALSE], digits = digits, row.names = row.names,
            ...)
    }
    if (any(!x$component_diagnostics$PositiveSemidefinite))
        cat("Warning: non-PSD covariance components; inspect component_diagnostics.\n")
    cat("G/Phi are dimensionless; SEM is in score units. Available describes calculation, not coverage.\n")
    cat("For differences, positive G/Phi and negative SEM favor the comparison plan.\n")
    cat("Method assumptions and limitations: $notes; complete saved results: $intervals.\n")
    invisible(x)
}
