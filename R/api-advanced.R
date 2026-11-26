


analyze_dff <- function(fit, facet, group, data = NULL, focal = NULL,
    method = c("mml_lrt", "cml_lrt", "raju_welch", "mh", "logistic"),
    df_method = c("facets", "classical"), min_obs = 10, p_adjust = "holm", threads = 1L) {
    .mfrm_validate_fit(fit)
    .mfrm_analyze_dff(fit = fit, facet = facet, group = group, data = data, focal = focal,
        method = match.arg(method), df_method = match.arg(df_method), min_obs = min_obs,
        p_adjust = p_adjust, threads = threads)
}


summary.mfrm_dff <- function(object, ...) {
    tbl <- object$table
    flagged <- if ("Flagged" %in% names(tbl))
        tbl[!is.na(tbl$Flagged) & tbl$Flagged, , drop = FALSE]
    else tbl[FALSE, , drop = FALSE]
    ov <- data.frame(Method = as.character(object$method %||% "mml_lrt"), Facet = as.character(object$facet %||%
        ""), Group = as.character(object$group %||% ""), Reference = as.character(object$reference %||% ""),
        Focal = as.character(object$focal %||% ""), ElementsTested = sum(is.finite(tbl$pValue)), ElementsUnavailable = sum(!is.finite(tbl$pValue)),
        Flagged = if (all(is.na(tbl$Flagged)))
            NA_integer_
        else nrow(flagged), MaxAbsContrast = if (nrow(tbl) > 0L && any(is.finite(tbl$Contrast)))
            max(abs(tbl$Contrast), na.rm = TRUE)
        else NA_real_, stringsAsFactors = FALSE)
    return(structure(list(overview = ov, summary = ov, dif_table = tbl, table = tbl, flagged = flagged),
        class = "summary.mfrm_dff"))
}


print.summary.mfrm_dff <- function(x, ...) {
    ov <- x$overview
    cat(sprintf("Differential Facet Functioning Summary [%s | Facet: %s | Group: %s (%s vs %s)]\n", ov$Method,
        ov$Facet, ov$Group, ov$Reference, ov$Focal))
    cat(sprintf("  Elements tested: %d | Unavailable: %d | Flagged (pAdj < 0.05): %d | Max |Contrast|: %.4f\n\n",
        ov$ElementsTested, ov$ElementsUnavailable %||% 0L, ov$Flagged, ov$MaxAbsContrast))
    if (nrow(x$table) > 0L) {
        disp <- x$table
        nc <- vapply(disp, is.numeric, logical(1))
        disp[nc] <- lapply(disp[nc], round, 4)
        print(disp, row.names = FALSE)
    }
    return(invisible(x))
}




compute_information <- function(fit, theta_range = c(-6, 6), theta_points = 201L, rating_plan = NULL, exposure_basis = c("observed_total",
    "per_person_average")) {
    .mfrm_validate_fit(fit)
    basis <- if (missing(exposure_basis))
        "unique_contexts"
    else match.arg(exposure_basis)
    return(.mfrm_information_fit(fit, theta_range, theta_points, rating_plan, basis))
}

anchor_to_baseline <- function(new_data, baseline_fit, person, facets, score, anchor_facets = NULL, include_person = FALSE,
    weight = NULL, model = NULL, method = NULL, anchor_policy = "warn", mml_integration = NULL, centering_scope = NULL,
    ..., drift_reference = NULL) {
    .mfrm_validate_fit(baseline_fit)
    anchor_policy <- match.arg(anchor_policy, c("warn", "error", "silent"))
    model <- toupper(model %||% baseline_fit$model)
    method <- method %||% baseline_fit$method
    anchor_tbl <- as.data.frame(make_anchor_table(baseline_fit, facets = anchor_facets, include_person = include_person))
    if (!nrow(anchor_tbl))
        stop("No anchors could be extracted from the baseline.", call. = FALSE)
    args <- c(list(data = new_data, person = person, facets = facets, score = score, weight = weight, model = model,
        method = method), list(...))
    if (is.null(args$step_facet) && identical(model, baseline_fit$model)) {
        args$step_facet <- if (identical(model, "GPCM"))
            baseline_fit$encoding$step_facet %||% "Common"
        else baseline_fit$encoding$step_facet
    }
    if (is.null(args$slope_facet) && identical(model, "GPCM")) {
        args$slope_facet <- baseline_fit$encoding$slope_facets
        args$spec$slope_reference <- args$spec$slope_reference %||% baseline_fit$encoding$slope_reference
        args$spec$slope_type <- args$spec$slope_type %||% baseline_fit$encoding$slope_type
    }
    args$rating_min <- args$rating_min %||% baseline_fit$rating_min
    args$rating_max <- args$rating_max %||% baseline_fit$rating_max
    if (!is.null(centering_scope))
        args$spec$centering_scope <- centering_scope
    args$spec$centering_scope <- args$spec$centering_scope %||% "remainder"
    if (!is.null(mml_integration))
        args$control$mml_integration <- mml_integration
    cols <- c(person, facets, score, weight)
    if (any(!cols %in% names(new_data)))
        stop("Missing rating columns in new_data.", call. = FALSE)
    keep <- is.finite(as.numeric(new_data[[score]])) & stats::complete.cases(new_data[cols])
    if (!is.null(weight))
        keep <- keep & is.finite(new_data[[weight]]) & new_data[[weight]] > 0
    present <- vapply(seq_len(nrow(anchor_tbl)), function(i) {
        f <- anchor_tbl$Facet[i]
        f %in% facets && anchor_tbl$Level[i] %in% as.character(new_data[[f]][keep])
    }, logical(1))
    if (any(!present)) {
        msg <- paste("Baseline anchor elements absent or unsupported in the retained ratings:", paste(paste(anchor_tbl$Facet[!present],
            anchor_tbl$Level[!present], sep = "/"), collapse = ", "))
        if (anchor_policy == "error")
            stop(msg, call. = FALSE)
        if (anchor_policy == "warn")
            warning(msg, call. = FALSE)
    }
    anchors_used <- anchor_tbl[present, , drop = FALSE]
    if (!nrow(anchors_used))
        stop("No baseline anchors are available in the retained ratings.", call. = FALSE)
    new_fit <- do.call(fit_mfrm, c(args, list(anchors = anchors_used)))
    drift_df <- data.frame(Facet = anchors_used$Facet, Level = anchors_used$Level, Baseline = anchors_used$Anchor,
        NewFree = NA_real_, SE_New = NA_real_, Drift = NA_real_, SE_Diff = NA_real_, Drift_SE_Ratio = NA_real_,
        Flag = NA, Status = "independent_reference_required", stringsAsFactors = FALSE)
    if (!is.null(drift_reference)) {
        if (!is.data.frame(drift_reference) || !all(c("Facet", "Level") %in% names(drift_reference)) || !nrow(drift_reference) ||
            anyNA(drift_reference[c("Facet", "Level")]) || anyDuplicated(drift_reference[c("Facet", "Level")])) {
            stop("Supply unique Facet and Level rows in drift_reference.", call. = FALSE)
        }
        ref <- merge(drift_reference[c("Facet", "Level")], anchors_used[c("Facet", "Level", "Anchor")], by = c("Facet",
            "Level"))
        if (nrow(ref) != nrow(drift_reference))
            stop("Drift reference elements require retained baseline anchors.", call. = FALSE)
        if (!identical(model, baseline_fit$model) || !identical(method, baseline_fit$method) || !model %in%
            c("RSM", "PCM") || length(baseline_fit$encoding$explanatory_meta) || length(baseline_fit$encoding$interaction_meta) ||
            isTRUE(baseline_fit$encoding$pop_mode == 1L)) {
            stop("Drift comparisons with reference anchors require the same additive RSM/PCM estimator and step structure.",
                call. = FALSE)
        }
        candidates <- match(paste(drift_df$Facet, drift_df$Level, sep = "/"), paste(ref$Facet, ref$Level,
            sep = "/"))
        drift_df$Status[!is.na(candidates)] <- "reference_anchor"
        compare <- is.na(candidates) & drift_df$Facet %in% ref$Facet
        args$warn_connectivity <- FALSE
        args$spec$noncenter_facet <- union(args$spec$noncenter_facet %||% "Person", ref$Facet)
        reference_fit <- do.call(fit_mfrm, c(args, list(anchors = ref)))
        if (any(compare) && .mfrm_readiness(reference_fit)$ready && .mfrm_readiness(baseline_fit)$ready &&
            identical(reference_fit$encoding$step_facet, baseline_fit$encoding$step_facet) && identical(reference_fit$encoding$facet_signs,
            baseline_fit$encoding$facet_signs) && all(c(reference_fit$connectivity$n_subsets, baseline_fit$connectivity$n_subsets) ==
            1L)) {
            idx <- match(paste(drift_df$Facet, drift_df$Level, sep = "/"), paste(reference_fit$facets$Facet,
                reference_fit$facets$Level, sep = "/"))
            drift_df$NewFree[compare] <- reference_fit$facets$Estimate[idx[compare]]
            drift_df$SE_New[compare] <- reference_fit$facets$ModelSE[idx[compare]]
            drift_df$Drift[compare] <- drift_df$NewFree[compare] - drift_df$Baseline[compare]
            drift_df$SE_Diff[compare] <- drift_df$SE_New[compare]
            available <- compare & is.finite(drift_df$Drift) & is.finite(drift_df$SE_Diff) & drift_df$SE_Diff >
                0
            drift_df$Drift_SE_Ratio[available] <- abs(drift_df$Drift[available])/drift_df$SE_Diff[available]
            drift_df$Flag[available] <- abs(drift_df$Drift[available]) > 0.5 | drift_df$Drift_SE_Ratio[available] >
                2
            drift_df$Status[available] <- "conditional_on_baseline"
            drift_df$Status[compare & !available] <- "comparison_unavailable"
        }
        else {
            drift_df$Status[compare] <- "comparison_unavailable"
        }
    }
    return(structure(list(fit = new_fit, diagnostics = diagnose_mfrm(new_fit), baseline_anchors = anchor_tbl,
        drift = drift_df), class = c("mfrm_anchored_fit", "list")))
}

summary.mfrm_anchored_fit <- function(object, ...) {
    drift <- object$drift
    n_anchored <- nrow(object$baseline_anchors)
    n_common <- nrow(drift)
    n_available <- sum(!is.na(drift$Flag))
    n_flagged <- if (n_available) sum(drift$Flag, na.rm = TRUE) else NA_integer_
    .mfrm_validate_fit(object$fit)
    out <- list(n_anchored = n_anchored, n_common = n_common, n_available = n_available, n_flagged = n_flagged, drift_summary = if (n_common > 0L) {
        dplyr::summarise(dplyr::group_by(drift, .data$Facet), N = dplyr::n(), N_Available = sum(is.finite(.data$Drift)),
            Mean_Drift = if (any(is.finite(.data$Drift))) mean(abs(.data$Drift), na.rm = TRUE) else NA_real_,
            Max_Drift = if (any(is.finite(.data$Drift))) max(abs(.data$Drift), na.rm = TRUE) else NA_real_, N_Flagged = if (all(is.na(.data$Flag))) NA_integer_ else sum(.data$Flag, na.rm = TRUE), .groups = "drop")
    } else {
        tibble::tibble()
    }, flagged = dplyr::filter(drift, .data$Flag),
        converged = isTRUE(object$fit$converged),
        optimizer_code_zero = isTRUE(object$fit$converged),
        convergence_status = if (isTRUE(object$fit$converged)) "converged" else "not_converged")
    class(out) <- "summary.mfrm_anchored_fit"
    out
}

detect_anchor_drift <- function(fits, facets = NULL, drift_threshold = 0.5, flag_se_ratio = 2, reference = 1L,
    include_person = FALSE) {
    if (!is.list(fits) || length(fits) < 2L)
        stop("Supply at least two current fit_mfrm() results.", call. = FALSE)
    invisible(lapply(fits, .mfrm_validate_fit))
    for (threshold in list(drift_threshold, flag_se_ratio)) {
        if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold < 0)
            stop("Drift thresholds must be finite nonnegative numbers.", call. = FALSE)
    }
    if (!is.logical(include_person) || length(include_person) != 1L || is.na(include_person))
        stop("include_person must be TRUE or FALSE.", call. = FALSE)
    if (is.null(names(fits)))
        names(fits) <- paste0("Wave", seq_along(fits))
    if (anyNA(names(fits)) || any(!nzchar(names(fits))) || anyDuplicated(names(fits)))
        stop("Wave names must be nonempty and unique.", call. = FALSE)
    ref_idx <- if (is.character(reference))
        match(reference, names(fits))
    else reference
    if (!is.numeric(ref_idx) || length(ref_idx) != 1L || !is.finite(ref_idx) ||
        ref_idx != as.integer(ref_idx) || ref_idx < 1L || ref_idx > length(fits))
        stop("reference must identify one supplied wave by name or integer index.", call. = FALSE)
    element_table <- function(fit) {
        out <- fit$facets[, c("Facet", "Level", "Estimate", "ModelSE"), drop = FALSE]
        if (include_person) {
            out <- rbind(out, data.frame(Facet = "Person", Level = fit$person$Person,
                Estimate = fit$person$Estimate, ModelSE = fit$person$ModelSE))
        }
        out
    }
    ref_f <- fits[[ref_idx]]
    ref_tbl <- element_table(ref_f)
    if (!is.null(facets))
        ref_tbl <- ref_tbl[ref_tbl$Facet %in% facets, , drop = FALSE]
    rows <- list()
    for (w in setdiff(seq_along(fits), ref_idx)) {
        wf <- element_table(fits[[w]])
        m <- match(paste(ref_tbl$Facet, ref_tbl$Level, sep = "\r"), paste(wf$Facet, wf$Level, sep = "\r"))
        ok <- which(!is.na(m))
        if (length(ok) == 0L)
            next
        d_val <- wf$Estimate[m[ok]] - ref_tbl$Estimate[ok]
        se_d <- sqrt(ref_tbl$ModelSE[ok]^2 + wf$ModelSE[m[ok]]^2)
        if (!.mfrm_readiness(ref_f)$ready || !.mfrm_readiness(fits[[w]])$ready)
            se_d[] <- NA_real_
        z_val <- ifelse(se_d > 0, d_val/se_d, NA_real_)
        available <- is.finite(z_val)
        rows[[length(rows) + 1L]] <- data.frame(Reference = names(fits)[ref_idx], Wave = names(fits)[w], Facet = ref_tbl$Facet[ok],
            Level = ref_tbl$Level[ok], EstRef = ref_tbl$Estimate[ok], EstWave = wf$Estimate[m[ok]], Drift = d_val,
            SE = se_d, z = z_val, pValue = 2 * stats::pnorm(-abs(z_val)),
            Flagged = ifelse(available, abs(d_val) > drift_threshold | abs(z_val) > flag_se_ratio, NA),
            Status = ifelse(available, "independent_waves_common_scale_assumed", "inference_unavailable"),
            stringsAsFactors = FALSE)
    }
    tbl <- if (length(rows) > 0L)
        do.call(rbind, rows)
    else data.frame(Reference = character(), Wave = character(), Facet = character(),
        Level = character(), EstRef = numeric(), EstWave = numeric(), Drift = numeric(),
        SE = numeric(), z = numeric(), pValue = numeric(), Flagged = logical(), Status = character())
    return(structure(list(table = tbl, reference = names(fits)[ref_idx]), class = c("mfrm_anchor_drift", "mfrm_comparison")))
}

summary.mfrm_anchor_drift <- function(object, ...) {
    dt <- object$table
    out <- list(reference = object$reference, n_comparisons = nrow(dt),
        n_flagged = if (anyNA(dt$Flagged)) NA_integer_ else sum(dt$Flagged),
        n_flagged_available = sum(dt$Flagged, na.rm = TRUE), n_unavailable = sum(is.na(dt$Flagged)),
        table = dt, flagged = dt[which(dt$Flagged), , drop = FALSE])
    class(out) <- "summary.mfrm_anchor_drift"
    out
}
