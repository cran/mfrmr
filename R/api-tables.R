


fair_average_table <- function(fit, diagnostics = NULL, facets = NULL, umean = 0, uscale = 1, fair_se = FALSE,
    ci_level = 0.95, reference_plan = NULL) {
    .mfrm_validate_fit(fit)
    if (!is.null(diagnostics) && !identical(diagnostics$fit, fit))
        stop("Diagnostics must match the source calibration.", call. = FALSE)
    if (isTRUE(fair_se) && is.null(reference_plan))
        stop("Fair-average intervals are not implemented without a reference_plan; measure SEs are not fair-score SEs.",
            call. = FALSE)
    if (isTRUE(fair_se) && (length(ci_level) != 1L || !is.numeric(ci_level) || !is.finite(ci_level) || ci_level <=
        0 || ci_level >= 1))
        stop("ci_level must be between 0 and 1.", call. = FALSE)
    p_raw <- data.frame(Facet = "Person", Level = fit$person$Person, ObsCount = fit$person$ObsCount, ObsAvg = fit$person$ObsAvg,
        ExpAvg = fit$person$ExpAvg, FairM = fit$person$FairM, FairZ = fit$person$FairZ, Measure = umean +
            uscale * fit$person$Estimate, SE = uscale * fit$person$ModelSE, stringsAsFactors = FALSE)
    by_f <- list(Person = p_raw)
    for (fnm in unique(fit$facets$Facet)) {
        sub <- fit$facets[fit$facets$Facet == fnm, ]
        by_f[[fnm]] <- data.frame(Facet = fnm, Level = sub$Level, ObsCount = sub$ObsCount, ObsAvg = sub$ObsAvg,
            ExpAvg = sub$ExpAvg, FairM = sub$FairM, FairZ = sub$FairZ, Measure = umean + uscale * sub$Estimate,
            SE = uscale * sub$ModelSE, stringsAsFactors = FALSE)
    }
    if (!is.null(facets))
        by_f <- by_f[intersect(as.character(facets), names(by_f))]
    stacked <- do.call(rbind, by_f)
    rownames(stacked) <- NULL
    bundle <- list(stacked = stacked, raw_by_facet = by_f)
    if (!is.null(reference_plan)) {
        plan <- mfrm_fair_reference_plan(fit, reference_plan)
        bundle$reference_scores <- mfrm_fair_reference_scores(fit, plan, if (isTRUE(fair_se))
            ci_level)
        bundle$reference_plan <- plan
        bundle$reference_fit <- fit
        bundle$settings <- list(reference_target = "common_plan_expected_average", score_method = fit$score_method)
    }
    return(structure(bundle, class = c("mfrm_fair_average", "mfrm_bundle")))
}

score_category_support_profile <- function(prep = NULL, score_distribution = NULL) {
    support <- integer(0)
    observed <- integer(0)
    unused <- integer(0)
    if (!is.null(prep)) {
        rating_min <- suppressWarnings(as.integer(prep$rating_min %||% NA_integer_))
        rating_max <- suppressWarnings(as.integer(prep$rating_max %||% NA_integer_))
        if (is.finite(rating_min) && is.finite(rating_max) && rating_min <= rating_max) {
            support <- seq(rating_min, rating_max)
        }
        prep_data <- as.data.frame(prep$data %||% data.frame(), stringsAsFactors = FALSE)
        if (nrow(prep_data) > 0 && "Score" %in% names(prep_data)) {
            observed <- sort(unique(suppressWarnings(as.integer(prep_data$Score))))
            observed <- observed[is.finite(observed)]
        }
        unused <- sort(unique(suppressWarnings(as.integer(prep$unused_score_categories %||% integer(0)))))
        unused <- unused[is.finite(unused)]
    }
    score_distribution <- as.data.frame(score_distribution %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(score_distribution) > 0) {
        category_col <- if ("Category" %in% names(score_distribution)) {
            "Category"
        }
        else if ("Score" %in% names(score_distribution)) {
            "Score"
        }
        else {
            NA_character_
        }
        if (!is.na(category_col)) {
            dist_categories <- suppressWarnings(as.integer(score_distribution[[category_col]]))
            dist_categories <- dist_categories[is.finite(dist_categories)]
            support <- sort(unique(c(support, dist_categories)))
            raw_n <- suppressWarnings(as.numeric(score_distribution$RawN %||% score_distribution$Count %||%
                NA_real_))
            weighted_n <- suppressWarnings(as.numeric(score_distribution$WeightedN %||% raw_n))
            zero_mask <- (is.finite(raw_n) & raw_n <= 0) | (is.finite(weighted_n) & weighted_n <= 0)
            if (length(zero_mask) == nrow(score_distribution)) {
                zero_categories <- suppressWarnings(as.integer(score_distribution[[category_col]][zero_mask]))
                zero_categories <- zero_categories[is.finite(zero_categories)]
                unused <- sort(unique(c(unused, zero_categories)))
                positive_categories <- suppressWarnings(as.integer(score_distribution[[category_col]][!zero_mask]))
                positive_categories <- positive_categories[is.finite(positive_categories)]
                observed <- sort(unique(c(observed, positive_categories)))
            }
        }
    }
    if (length(support) == 0L) {
        support <- sort(unique(c(observed, unused)))
    }
    support <- support[is.finite(support)]
    if (length(support) == 0L) {
        return(tibble::tibble(Category = integer(0), ZeroCount = logical(0), UnusedCategoryType = character(0),
            WeaklyIdentified = logical(0), CategoryCaveat = character(0)))
    }
    if (length(observed) == 0L && length(unused) > 0L) {
        observed <- setdiff(support, unused)
    }
    observed_min <- suppressWarnings(min(observed, na.rm = TRUE))
    observed_max <- suppressWarnings(max(observed, na.rm = TRUE))
    zero_count <- support %in% unused
    internal <- if (is.finite(observed_min) && is.finite(observed_max)) {
        support > observed_min & support < observed_max
    }
    else {
        rep(FALSE, length(support))
    }
    unused_type <- ifelse(zero_count & internal, "internal", ifelse(zero_count, "boundary", "none"))
    caveat <- dplyr::case_when(unused_type == "internal" ~ "Zero-count intermediate category; adjacent thresholds are weakly identified.",
        unused_type == "boundary" ~ "Zero-count boundary category; document the retained support and avoid overinterpreting adjacent thresholds.",
        TRUE ~ "")
    tibble::tibble(Category = support, ZeroCount = zero_count, UnusedCategoryType = unused_type, WeaklyIdentified = zero_count,
        CategoryCaveat = caveat)
}

empty_mfrm_caveats <- .mfrmr_empty_table(c(Area = "character", Severity = "character", Condition = "character",
    Categories = "character", CategoryType = "character", Message = "character", RecommendedAction = "character",
    Details = "character"))

collect_mfrm_caveats <- function(prep = NULL, score_distribution = NULL, include_recode = TRUE,
    context = c("fit", "data")) {
    context <- match.arg(context)
    profile <- score_category_support_profile(prep = prep, score_distribution = score_distribution)
    out <- empty_mfrm_caveats
    support_phrase <- if (identical(context, "data"))
        "prepared score support"
    else "fitted score support"
    add_caveat <- function(severity, condition, categories, category_type, message, action) {
        out <<- rbind(out, data.frame(Area = "score_categories", Severity = severity, Condition = condition,
            Categories = paste(as.character(categories), collapse = ", "), CategoryType = category_type,
            Message = message, RecommendedAction = action, Details = "", stringsAsFactors = FALSE))
        invisible(NULL)
    }
    zero_rows <- if (nrow(profile) > 0 && "ZeroCount" %in% names(profile)) {
        profile[as.logical(profile$ZeroCount), , drop = FALSE]
    }
    else {
        profile[0, , drop = FALSE]
    }
    if (nrow(zero_rows) > 0) {
        internal <- zero_rows$Category[zero_rows$UnusedCategoryType == "internal"]
        boundary <- zero_rows$Category[zero_rows$UnusedCategoryType == "boundary"]
        if (length(internal) > 0L) {
            internal_text <- paste(as.character(internal), collapse = ", ")
            add_caveat(severity = "warning", condition = "zero_count_intermediate_score_category", categories = internal,
                category_type = "internal", message = paste0("Unused intermediate score categories retained in the ",
                  support_phrase, ": ", internal_text, ". Under the current sum-zero step parameterization, a zero-count internal category creates an unsupported adjacent-step recession direction in its fitted ladder. New fits stop before optimization; legacy or descriptive objects must not present the adjacent steps as finite identified estimates."),
                action = "Collect observations in the internal category, revise the score support or response model, or separate the affected ladder before fitting.")
        }
        if (length(boundary) > 0L) {
            boundary_text <- paste(as.character(boundary), collapse = ", ")
            add_caveat(severity = "review", condition = "zero_count_boundary_score_category", categories = boundary,
                category_type = "boundary", message = paste0("Unused boundary score categories retained in the ",
                  support_phrase, ": ", boundary_text, ". Document the zero-count category and avoid overinterpreting adjacent threshold estimates."),
                action = "Document the intended score support and verify `rating_min` / `rating_max` before reporting category functioning.")
        }
    }
    if (isTRUE(include_recode) && !is.null(prep)) {
        score_map <- as.data.frame(prep$score_map %||% data.frame(), stringsAsFactors = FALSE)
        if (nrow(score_map) > 0L && all(c("OriginalScore", "InternalScore") %in% names(score_map)) &&
            any(as.character(score_map$OriginalScore) != as.character(score_map$InternalScore))) {
            add_caveat(severity = "info", condition = "score_categories_recoded", categories = sort(unique(as.character(score_map$OriginalScore))),
                category_type = "recoded", message = "Observed score categories were internally recoded; inspect `describe_mfrm_data(...)$score_support$score_map` before interpreting category labels.",
                action = "Use `describe_mfrm_data(...)$score_support$score_map` when connecting internal category estimates back to original score labels.")
        }
        if (nrow(score_map) > 0L && "OriginalScore" %in% names(score_map)) {
            original_scores <- sort(unique(suppressWarnings(as.integer(score_map$OriginalScore))))
            original_scores <- original_scores[is.finite(original_scores)]
            if (length(original_scores) >= 2L) {
                original_support <- seq(min(original_scores), max(original_scores))
                missing_original <- setdiff(original_support, original_scores)
                retained_zero <- suppressWarnings(as.integer(zero_rows$Category %||% integer(0)))
                retained_zero <- retained_zero[is.finite(retained_zero)]
                missing_not_retained <- setdiff(missing_original, retained_zero)
                if (length(missing_not_retained) > 0L) {
                  add_caveat(severity = "warning", condition = "original_score_gap_before_recoding",
                    categories = missing_not_retained, category_type = "original_gap", message = paste0("Original score labels skip category ",
                      paste(missing_not_retained, collapse = ", "), " before the fitted internal scale is formed. If the intended scale includes ",
                      if (length(missing_not_retained) == 1L)
                        "this category"
                      else "these categories", ", refit with `rating_min`, `rating_max`, and `keep_original = TRUE` so zero-frequency categories remain visible in category-functioning output."),
                    action = "Check the declared rating scale; retain intended zero-count intermediate categories when threshold functioning is part of the analysis.")
                }
            }
        }
    }
    out
}

estimate_bias <- function(fit, facet_a = NULL, facet_b = NULL, interaction_facets = NULL,
    method = c("local", "joint"), df_method = c("facets", "classical"), p_adjust = "holm") {
    .mfrm_validate_fit(fit)
    .mfrm_estimate_bias_cpp(fit = fit, facet_a = facet_a, facet_b = facet_b,
        interaction_facets = interaction_facets, method = match.arg(method),
        df_method = match.arg(df_method), p_adjust = p_adjust)
}
