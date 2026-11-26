compute_person_fit_indices <- function(diagnostics, fit = NULL) {
    if (inherits(diagnostics, "mfrm_fit") || inherits(diagnostics, "mfrm_diagnostics")) {
        source_fit <- if (inherits(diagnostics, "mfrm_fit"))
            diagnostics
        else diagnostics$fit
        .mfrm_validate_fit(source_fit)
        if (!is.null(fit) && !identical(fit, source_fit))
            stop("`fit` must match the source calibration in `diagnostics`.", call. = FALSE)
        fit <- source_fit
        diagnostics <- list(obs = .mfrm_person_fit_moments(fit))
    }
    if (!is.list(diagnostics) || is.null(diagnostics$obs))
        stop("Supply a fitted calibration or its diagnostics.", call. = FALSE)
    obs <- diagnostics$obs
    obs <- as.data.frame(obs, stringsAsFactors = FALSE)
    if (!nrow(obs) || anyNA(obs$Person) || any(!nzchar(as.character(obs$Person)))) {
        stop("Person-fit calculations require observations with non-missing Person identifiers.", call. = FALSE)
    }
    needed <- c("Person", "PrObserved", "ItemEntropy", "ItemVarLogP")
    missing_cols <- setdiff(needed, names(obs))
    if (length(missing_cols) > 0L) {
        stop("`diagnostics$obs` is missing required columns: ", paste(missing_cols, collapse = ", "), ". This typically means the diagnostics bundle was generated ",
            "by an older mfrmr version; re-run diagnose_mfrm(fit, ...) ", "to populate the per-observation probability columns.",
            call. = FALSE)
    }
    obs$Person <- as.character(obs$Person)
    obs$PrObserved <- suppressWarnings(as.numeric(obs$PrObserved))
    obs$ItemEntropy <- suppressWarnings(as.numeric(obs$ItemEntropy))
    obs$ItemVarLogP <- suppressWarnings(as.numeric(obs$ItemVarLogP))
    valid_probability <- is.finite(obs$PrObserved) & obs$PrObserved > 0 & obs$PrObserved <= 1
    log_p <- rep(NA_real_, nrow(obs))
    log_p[valid_probability] <- log(obs$PrObserved[valid_probability])
    agg <- by(data.frame(log_p = log_p, e_logp = obs$ItemEntropy, var_logp = obs$ItemVarLogP, stringsAsFactors = FALSE),
        obs$Person, function(d) {
            ok <- is.finite(d$log_p) & is.finite(d$e_logp) & is.finite(d$var_logp) & d$var_logp >= 0
            total <- nrow(d)
            n <- sum(ok)
            if (n < total) {
                return(c(TotalN = total, N = n, UnavailableN = total - n, LogLik = NA_real_, lz = NA_real_))
            }
            ll <- sum(d$log_p)
            e_ll <- sum(d$e_logp)
            var_ll <- sum(d$var_logp)
            lz_val <- if (is.finite(var_ll) && var_ll > 0)
                (ll - e_ll)/sqrt(var_ll)
            else NA_real_
            c(TotalN = total, N = n, UnavailableN = 0, LogLik = ll, lz = lz_val)
        })
    out <- do.call(rbind, lapply(agg, function(x) as.data.frame(t(x))))
    out$Person <- rownames(out)
    rownames(out) <- NULL
    out <- out[, c("Person", "TotalN", "N", "UnavailableN", "LogLik", "lz"), drop = FALSE]
    out$N <- as.integer(out$N)
    out$TotalN <- as.integer(out$TotalN)
    out$UnavailableN <- as.integer(out$UnavailableN)
    out$lz_status <- ifelse(out$UnavailableN > 0L, "incomplete_observation_terms", ifelse(is.finite(out$lz),
        "computed", "insufficient_variance"))
    numerical_ready <- is.null(fit) || isTRUE(fit$converged)
    lz_star_tbl <- compute_snijders_lz_star(obs, out$Person, fit = if (numerical_ready)
        fit
    else NULL)
    if (!numerical_ready)
        lz_star_tbl$lz_star_status[] <- "fit_not_converged"
    lz_idx <- match(out$Person, lz_star_tbl$Person)
    out$lz_star <- lz_star_tbl$lz_star[lz_idx]
    out$lz_star_status <- lz_star_tbl$lz_star_status[lz_idx]
    out$lz_star_c <- lz_star_tbl$lz_star_c[lz_idx]
    out$lz_star_variance <- lz_star_tbl$lz_star_variance[lz_idx]
    out <- add_person_fit_reporting_columns(out)
    out <- out[, c("Person", "TotalN", "N", "UnavailableN", "LogLik", "lz", "lz_status", "lz_star", "lz_star_status",
        "lz_star_c", "lz_star_variance", "lz_flag_5pct", "lz_flag_1pct", "lz_star_flag_5pct", "lz_star_flag_1pct",
        "ReportIndex", "ReportValue", "ReportFlagLevel", "ReportFlag", "ReviewStatus", "ReviewReason", "ReportCaveat"),
        drop = FALSE]
    out$InfitMnSq <- out$OutfitMnSq <- NA_real_
    if (!is.null(fit)) {
        ids <- match(out$Person, fit$person$Person)
        out$InfitMnSq <- fit$person$InfitMnSq[ids]
        out$OutfitMnSq <- fit$person$OutfitMnSq[ids]
    }
    out$ObsCount <- out$N
    out$L0 <- out$LogLik
    out$Flagged <- out$ReportFlag
    attr(out, "person_fit_thresholds") <- .mfrmr_const_person_fit_threshold_table
    attr(out, "person_fit_calculation_version") <- 2L
    class(out) <- c("mfrm_person_fit_indices", class(out))
    out
}

add_person_fit_reporting_columns <- function(out) {
    z_5pct <- stats::qnorm(0.975)
    z_1pct <- stats::qnorm(0.995)
    out$lz_flag_5pct <- ifelse(is.finite(out$lz), abs(out$lz) > z_5pct, NA)
    out$lz_flag_1pct <- ifelse(is.finite(out$lz), abs(out$lz) > z_1pct, NA)
    out$lz_star_flag_5pct <- ifelse(is.finite(out$lz_star), abs(out$lz_star) > z_5pct, NA)
    out$lz_star_flag_1pct <- ifelse(is.finite(out$lz_star), abs(out$lz_star) > z_1pct, NA)
    snijders_ready <- is.finite(out$lz_star) & identical_length_status(out$lz_star_status, "computed_jml_conditional_calibration")
    lz_ready <- !snijders_ready & is.finite(out$lz)
    out$ReportIndex <- ifelse(snijders_ready, "lz_star", ifelse(lz_ready, "lz", "none"))
    out$ReportValue <- ifelse(snijders_ready, out$lz_star, ifelse(lz_ready, out$lz, NA_real_))
    abs_report <- abs(out$ReportValue)
    out$ReportFlagLevel <- ifelse(!is.finite(abs_report), "not_available", ifelse(abs_report > z_1pct,
        "1pct", ifelse(abs_report > z_5pct, "5pct", "none")))
    out$ReportFlag <- ifelse(is.finite(out$ReportValue), out$ReportFlagLevel %in% c("5pct", "1pct"),
        NA)
    out$ReviewStatus <- ifelse(out$ReportFlagLevel == "1pct", "review_1pct", ifelse(out$ReportFlagLevel ==
        "5pct", "review_5pct", ifelse(out$ReportFlagLevel == "none", "not_flagged", "not_available")))
    threshold_reason <- ifelse(out$ReportFlagLevel == "1pct", paste0(ifelse(out$ReportIndex == "lz_star",
        "lz*", "lz"), " exceeds the absolute normal-reference cutoff 2.58."), ifelse(out$ReportFlagLevel ==
        "5pct", paste0(ifelse(out$ReportIndex == "lz_star", "lz*", "lz"), " exceeds the absolute normal-reference cutoff 1.96."),
        ifelse(out$ReportFlagLevel == "none", "No report-level flag under the practical two-sided thresholds.",
            "No finite person-fit statistic is available for report-level flagging.")))
    out$ReviewReason <- threshold_reason
    status <- as.character(out$lz_star_status)
    reason <- person_fit_correction_note(status)
    out$ReportCaveat <- ifelse(snijders_ready, "lz* conditions on the fitted non-person calibration; calibration uncertainty is omitted.",
        paste(ifelse(lz_ready, "lz is an uncorrected, unweighted response-pattern screen.", "No person-fit index is available."),
            reason))
    incomplete <- out$lz_status == "incomplete_observation_terms"
    out$ReportCaveat[incomplete] <- paste("Response probability or log-likelihood moment information is incomplete.",
        "A statistic is not calculated on a selected subset of responses.")
    out
}

person_fit_correction_note <- function(status) {
    notes <- c(computed_jml_conditional_calibration = "The ability correction conditions on the fitted calibration.",
        fit_required = "Supply the matching fit to assess the ability correction.", not_applicable_eap = "The implemented correction does not apply to MML/EAP person scores.",
        not_applicable_person_scoring = "The implemented correction requires unpenalized JML person scores; it does not cover WLE or adjusted/penalized scores.",
        nonunit_weights = "The implemented ability correction does not cover non-unit observation weights.",
        nonfinite_person_estimate = "The person estimate is not finite, so the ability correction is unavailable.",
        fit_not_converged = "Numerical convergence is not established; the ability correction is unavailable.",
        incomplete_observation_terms = "The ability correction requires complete probability and derivative information for all responses.",
        diagnostics_missing_snijders_terms = "Recompute diagnostics to supply the required score derivatives.",
        insufficient_information = "There is insufficient information for the ability correction.")
    out <- unname(notes[as.character(status)])
    out[is.na(out)] <- "The ability correction is unavailable for this result."
    out
}

identical_length_status <- function(x, value) {
    x <- as.character(x)
    !is.na(x) & x == value
}

person_fit_count_table <- function(x, col, label = col) {
    if (!col %in% names(x)) {
        return(data.frame(Variable = label, Value = character(), Rows = integer(), Proportion = numeric(),
            stringsAsFactors = FALSE))
    }
    vals <- as.character(x[[col]])
    vals[is.na(vals) | !nzchar(vals)] <- "missing"
    tab <- sort(table(vals), decreasing = TRUE)
    total <- sum(tab)
    data.frame(Variable = label, Value = names(tab), Rows = as.integer(tab), Proportion = if (total >
        0)
        as.numeric(tab)/total
    else NA_real_, stringsAsFactors = FALSE)
}

person_fit_review_priority <- function(status) {
    status <- as.character(status)
    out <- rep(0L, length(status))
    out[status == "review_5pct"] <- 1L
    out[status == "review_1pct"] <- 2L
    out[status == "not_available"] <- -1L
    out
}

person_fit_validate_table <- function(object) {
    if (!is.data.frame(object)) {
        stop("`object` must be a data frame returned by `compute_person_fit_indices()`.", call. = FALSE)
    }
    if (!identical(attr(object, "person_fit_calculation_version", exact = TRUE), 2L)) {
        stop("Recreate this person-fit result with compute_person_fit_indices(diagnostics, fit) to retain unavailable responses and flags.",
            call. = FALSE)
    }
    needed <- c("Person", "TotalN", "N", "UnavailableN", "LogLik", "lz", "lz_star", "lz_star_status",
        "ReportIndex", "ReportValue", "ReportFlagLevel", "ReportFlag", "ReviewStatus", "ReviewReason",
        "ReportCaveat")
    missing_cols <- setdiff(needed, names(object))
    if (length(missing_cols) > 0L) {
        stop("`object` is missing required person-fit columns: ", paste(missing_cols, collapse = ", "),
            ". Rebuild it with `compute_person_fit_indices()`.", call. = FALSE)
    }
    as.data.frame(object, stringsAsFactors = FALSE)
}

summary.mfrm_person_fit_indices <- function(object, digits = 3, top_n = 10, include_person = FALSE, ...) {
    tbl <- person_fit_validate_table(object)
    digits <- max(0L, as.integer(digits))
    top_n <- max(1L, as.integer(top_n))
    report_flag <- as.logical(tbl$ReportFlag)
    report_flag[is.na(report_flag)] <- FALSE
    status <- as.character(tbl$ReviewStatus)
    report_index <- as.character(tbl$ReportIndex)
    lz_star_status <- as.character(tbl$lz_star_status)
    report_value <- suppressWarnings(as.numeric(tbl$ReportValue))
    overview <- data.frame(Persons = nrow(tbl), ReportableRows = sum(is.finite(report_value)), ReportFlaggedRows = if (any(is.finite(report_value)))
        sum(report_flag)
    else NA_integer_, Review1PctRows = sum(status == "review_1pct", na.rm = TRUE), Review5PctRows = sum(status ==
        "review_5pct", na.rm = TRUE), NotFlaggedRows = sum(status == "not_flagged", na.rm = TRUE), NotAvailableRows = sum(status ==
        "not_available", na.rm = TRUE), SnijdersRows = sum(report_index == "lz_star", na.rm = TRUE),
        LzFallbackRows = sum(report_index == "lz", na.rm = TRUE), MissingReportIndexRows = sum(report_index ==
            "none", na.rm = TRUE), stringsAsFactors = FALSE)
    overview$FlagRate <- if (overview$ReportableRows > 0) {
        overview$ReportFlaggedRows/overview$ReportableRows
    }
    else {
        NA_real_
    }
    status_summary <- person_fit_count_table(tbl, "ReviewStatus", "ReviewStatus")
    report_index_summary <- person_fit_count_table(tbl, "ReportIndex", "ReportIndex")
    lz_star_status_summary <- person_fit_count_table(tbl, "lz_star_status", "lz_star_status")
    review_tbl <- tbl
    review_tbl$.priority <- person_fit_review_priority(review_tbl$ReviewStatus)
    review_tbl$.abs_report <- abs(suppressWarnings(as.numeric(review_tbl$ReportValue)))
    review_tbl$.abs_report[!is.finite(review_tbl$.abs_report)] <- -Inf
    review_tbl <- review_tbl[order(-review_tbl$.priority, -review_tbl$.abs_report, as.character(review_tbl$Person)),
        , drop = FALSE]
    review_tbl <- review_tbl[review_tbl$.priority > 0 | as.logical(review_tbl$ReportFlag %in% TRUE),
        , drop = FALSE]
    if (nrow(review_tbl) == 0L) {
        review_tbl <- tbl[order(-abs(report_value), as.character(tbl$Person)), , drop = FALSE]
    }
    keep_cols <- intersect(c("Person", "N", "ReportIndex", "ReportValue", "ReportFlagLevel", "ReviewStatus",
        "ReviewReason", "lz", "lz_star", "lz_star_status", "ReportCaveat"), names(review_tbl))
    top_review <- utils::head(review_tbl[, keep_cols, drop = FALSE], n = top_n)
    top_review[] <- lapply(top_review, function(x) {
        if (is.numeric(x))
            round(x, digits = digits)
        else x
    })
    caveats <- unique(tbl[, intersect(c("ReportIndex", "lz_star_status", "ReportCaveat"), names(tbl)),
        drop = FALSE])
    caveats <- caveats[order(as.character(caveats$ReportIndex), as.character(caveats$lz_star_status)),
        , drop = FALSE]
    thresholds <- attr(object, "person_fit_thresholds", exact = TRUE)
    if (!is.data.frame(thresholds) || nrow(thresholds) == 0L) {
        thresholds <- .mfrmr_const_person_fit_threshold_table
    }
    reporting_map <- data.frame(Area = c("First read", "Flagged person rows", "Response-level follow-up",
        "Plot/data handoff"), CoveredHere = c("yes", "yes", "no", "partial"), CompanionOutput = c("overview / status_summary / lz_star_status_summary",
        "top_review", "augment(workflow) / render_mfrm_report(workflow, rater = ID)", "plot(fit, type = \"person_fit\")"),
        stringsAsFactors = FALSE)
    notes <- c("The reported index uses lz* when its conditional ability correction is available; otherwise it uses uncorrected lz where possible.",
        "The two-sided normal cutoffs are screening references; neither individual nor multiple-person false-positive rates are calibrated for these fitted scores.",
        "Person-fit flags are screening evidence. Review response-level evidence before making substantive claims about a person.")
    out <- list(overview = overview, status_summary = status_summary, report_index_summary = report_index_summary,
        lz_star_status_summary = lz_star_status_summary, top_review = top_review, caveats = caveats,
        thresholds = thresholds, reporting_map = reporting_map, notes = notes, digits = digits, include_person = isTRUE(include_person),
        calculation_version = 2L)
    class(out) <- "summary.mfrm_person_fit_indices"
    out
}

print.summary.mfrm_person_fit_indices <- function(x, ...) {
    if (!identical(x$calculation_version, 2L)) {
        stop("Recreate this summary from compute_person_fit_indices(diagnostics, fit).", call. = FALSE)
    }
    cat("Person-Fit Summary\n")
    overview <- x$overview[, c("Persons", "ReportableRows", "NotAvailableRows", "ReportFlaggedRows",
        "FlagRate"), drop = FALSE]
    names(overview) <- c("Persons", "Available", "Unavailable", "Flagged", "Flag rate")
    print_bundle_section("Overview", overview, digits = x$digits)
    cat(sprintf("Ability-corrected lz*: %d; uncorrected lz: %d.\n", x$overview$SnijdersRows, x$overview$LzFallbackRows))
    if (isTRUE(x$include_person) && nrow(x$top_review)) {
        rows <- x$top_review[, c("Person", "N", "ReportIndex", "ReportValue", "ReviewReason"), drop = FALSE]
        rows$ReportIndex <- ifelse(rows$ReportIndex == "lz_star", "lz*", ifelse(rows$ReportIndex == "lz",
            "lz", "Unavailable"))
        names(rows) <- c("Person", "Responses", "Index", "Value", "Review")
        print_bundle_section("Person review", rows, digits = x$digits)
    }
    else {
        cat("Person identifiers are omitted; use include_person = TRUE to display the review rows.\n")
    }
    print_bullet_section("Interpretation", unique(c(x$caveats$ReportCaveat, x$notes)))
    invisible(x)
}

print.mfrm_person_fit_indices <- function(x, ...) {
    if (all(c("Person", "N", "LogLik", "lz_star", "ReviewStatus", "ReportValue", "ReportCaveat") %in%
        names(x))) {
        print(summary(x, ...))
    }
    else {
        print(as.data.frame(x), ...)
    }
    invisible(x)
}

compute_snijders_lz_star <- function(obs, persons, fit = NULL) {
    empty <- data.frame(Person = as.character(persons), lz_star = NA_real_, lz_star_status = "fit_required",
        lz_star_c = NA_real_, lz_star_variance = NA_real_, stringsAsFactors = FALSE)
    if (is.null(fit)) {
        return(empty)
    }
    if (!inherits(fit, "mfrm_fit")) {
        empty$lz_star_status <- "invalid_fit"
        return(empty)
    }
    .mfrm_validate_fit(fit)
    method <- fit$score_method
    if (identical(method, "JMLE")) method <- "JML"
    if (!identical(method, "JML")) {
        empty$lz_star_status <- if (method %in% c("MML", "EAP")) "not_applicable_eap" else
            "not_applicable_person_scoring"
        return(empty)
    }
    weights <- get_weights(obs)
    if (any(!is.finite(weights) | weights != 1)) {
        empty$lz_star_status <- "nonunit_weights"
        return(empty)
    }
    needed <- c("ItemLogPScoreCov", "ScoreInformation")
    missing_cols <- setdiff(needed, names(obs))
    if (length(missing_cols) > 0L) {
        empty$lz_star_status <- "diagnostics_missing_snijders_terms"
        return(empty)
    }
    obs$ItemLogPScoreCov <- suppressWarnings(as.numeric(obs$ItemLogPScoreCov))
    obs$ScoreInformation <- suppressWarnings(as.numeric(obs$ScoreInformation))
    if ("ObservedScoreDerivative" %in% names(obs)) {
        obs$ObservedScoreDerivative <- suppressWarnings(as.numeric(obs$ObservedScoreDerivative))
    }
    else if (all(c("ScoreSlope", "Observed", "Expected") %in% names(obs))) {
        obs$ObservedScoreDerivative <- suppressWarnings(as.numeric(obs$ScoreSlope) * (as.numeric(obs$Observed) -
            as.numeric(obs$Expected)))
    }
    else {
        empty$lz_star_status <- "diagnostics_missing_snijders_terms"
        return(empty)
    }
    valid_probability <- is.finite(obs$PrObserved) & obs$PrObserved > 0 & obs$PrObserved <= 1
    log_p <- rep(NA_real_, nrow(obs))
    log_p[valid_probability] <- log(obs$PrObserved[valid_probability])
    obs$.WCentered <- log_p - obs$ItemEntropy
    obs$Person <- as.character(obs$Person)
    rows <- split(obs, obs$Person)
    out <- empty
    for (i in seq_along(persons)) {
        person <- as.character(persons[i])
        d <- rows[[person]]
        if (is.null(d) || nrow(d) == 0L) {
            out$lz_star_status[i] <- "insufficient_information"
            next
        }
        ok <- is.finite(d$.WCentered) & is.finite(d$ItemVarLogP) & d$ItemVarLogP >= 0 & is.finite(d$ItemLogPScoreCov) &
            is.finite(d$ScoreInformation) & d$ScoreInformation >= 0 & is.finite(d$ObservedScoreDerivative)
        if (any(!ok)) {
            out$lz_star_status[i] <- "incomplete_observation_terms"
            next
        }
        if ("PersonMeasure" %in% names(d) && any(!is.finite(d$PersonMeasure))) {
            out$lz_star_status[i] <- "nonfinite_person_estimate"
            next
        }
        info_total <- sum(d$ScoreInformation)
        cov_total <- sum(d$ItemLogPScoreCov)
        var_logp_total <- sum(d$ItemVarLogP)
        if (!is.finite(info_total) || info_total <= 0 || !is.finite(var_logp_total) || var_logp_total <=
            0) {
            out$lz_star_status[i] <- "insufficient_information"
            next
        }
        c_n <- cov_total/info_total
        corrected_var <- var_logp_total - (cov_total^2/info_total)
        if (!is.finite(corrected_var) || corrected_var <= sqrt(.Machine$double.eps)) {
            out$lz_star_status[i] <- "insufficient_information"
            out$lz_star_c[i] <- c_n
            out$lz_star_variance[i] <- corrected_var
            next
        }
        centered_loglik <- sum(d$.WCentered)
        score_sum <- sum(d$ObservedScoreDerivative)
        out$lz_star[i] <- (centered_loglik - c_n * score_sum)/sqrt(corrected_var)
        out$lz_star_status[i] <- "computed_jml_conditional_calibration"
        out$lz_star_c[i] <- c_n
        out$lz_star_variance[i] <- corrected_var
    }
    out
}

# Both engines return a data frame; keep the earlier $table accessor for callers.
`$.mfrm_person_fit_indices` <- function(x, name) {
    if (identical(name, "table") && is.data.frame(x)) return(as.data.frame(x))
    .subset2(x, name)
}

.mfrm_person_fit_moments <- function(fit) {
    enc <- fit$encoding
    probabilities <- stats::predict(fit, type = "probs")
    categories <- seq_len(ncol(probabilities)) - 1L
    logp <- log(probabilities)
    logp[probabilities == 0] <- 0
    expected_logp <- rowSums(probabilities * logp)
    expected_score <- drop(probabilities %*% categories)
    derivative <- fit$alpha_obs * sweep(matrix(categories, nrow(probabilities),
        length(categories), byrow = TRUE), 1L, expected_score, "-")
    centered_logp <- logp - expected_logp
    data.frame(Person = enc$person_levels[enc$person0 + 1L],
        PrObserved = probabilities[cbind(seq_len(enc$n_obs), enc$score0 + 1L)],
        ItemEntropy = expected_logp,
        ItemVarLogP = rowSums(probabilities * centered_logp^2),
        ItemLogPScoreCov = rowSums(probabilities * centered_logp * derivative),
        ScoreInformation = rowSums(probabilities * derivative^2),
        ObservedScoreDerivative = fit$alpha_obs * (enc$score0 - expected_score),
        PersonMeasure = fit$person$Estimate[enc$person0 + 1L], Weight = enc$weight)
}
