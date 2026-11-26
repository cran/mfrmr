normalize_anchor_df <- function(df) {
    if (is.null(df) || nrow(df) == 0) {
        return(tibble::tibble(Facet = character(0), Level = character(0), Anchor = numeric(0)))
    }
    nm <- tolower(names(df))
    facet_col <- which(nm %in% c("facet", "facets"))
    level_col <- which(nm %in% c("level", "element", "label"))
    anchor_col <- which(nm %in% c("anchor", "value", "measure"))
    if (length(facet_col) == 0 || length(level_col) == 0 || length(anchor_col) == 0) {
        return(tibble::tibble(Facet = character(0), Level = character(0), Anchor = numeric(0)))
    }
    dplyr::filter(tibble::tibble(Facet = as.character(df[[facet_col[1]]]), Level = as.character(df[[level_col[1]]]),
        Anchor = suppressWarnings(as.numeric(df[[anchor_col[1]]]))), !is.na(Facet), !is.na(Level))
}

normalize_group_anchor_df <- function(df) {
    if (is.null(df) || nrow(df) == 0) {
        return(tibble::tibble(Facet = character(0), Level = character(0), Group = character(0), GroupValue = numeric(0)))
    }
    nm <- tolower(names(df))
    facet_col <- which(nm %in% c("facet", "facets"))
    level_col <- which(nm %in% c("level", "element", "label"))
    group_col <- which(nm %in% c("group", "subset"))
    value_col <- which(nm %in% c("groupvalue", "value", "anchor"))
    if (length(facet_col) == 0 || length(level_col) == 0 || length(group_col) == 0 || length(value_col) ==
        0) {
        return(tibble::tibble(Facet = character(0), Level = character(0), Group = character(0), GroupValue = numeric(0)))
    }
    dplyr::filter(tibble::tibble(Facet = as.character(df[[facet_col[1]]]), Level = as.character(df[[level_col[1]]]),
        Group = as.character(df[[group_col[1]]]), GroupValue = suppressWarnings(as.numeric(df[[value_col[1]]]))),
        !is.na(Facet), !is.na(Level), !is.na(Group))
}

collect_anchor_levels <- function(prep) {
    facets <- c("Person", prep$facet_names)
    rows <- lapply(facets, function(facet) {
        lv <- prep$levels[[facet]]
        if (is.null(lv) || length(lv) == 0)
            return(tibble::tibble())
        tibble::tibble(Facet = as.character(facet), Level = as.character(lv))
    })
    dplyr::bind_rows(rows)
}

safe_join_key <- function(facet, level) {
    paste0(as.character(facet), "\r", as.character(level))
}

build_anchor_issue_counts <- function(issue_tables) {
    issue_names <- names(issue_tables)
    tibble::tibble(Issue = issue_names, N = vapply(issue_tables, nrow, integer(1)))
}

build_anchor_recommendations <- function(facet_summary, issue_counts, design_checks = NULL, min_common_anchors = 5L,
    min_obs_per_element = 30, min_obs_per_category = 10, noncenter_facet = "Person", dummy_facets = character(0),
    centering_scope = "remainder") {
    rec <- character(0)
    if (!is.null(issue_counts) && nrow(issue_counts) > 0) {
        n_overlap <- issue_counts$N[issue_counts$Issue == "overlap_anchor_group"]
        n_missing <- issue_counts$N[issue_counts$Issue == "missing_group_values"]
        n_dup_anchor <- issue_counts$N[issue_counts$Issue == "duplicate_anchors"]
        n_dup_group <- issue_counts$N[issue_counts$Issue == "duplicate_group_assignments"]
        n_group_conf <- issue_counts$N[issue_counts$Issue == "group_value_conflicts"]
        if (length(n_overlap) > 0 && n_overlap > 0) {
            rec <- c(rec, "Levels listed in both direct- and group-anchor tables retain both constraints; verify that the fixed value and group mean are jointly compatible.")
        }
        if (length(n_missing) > 0 && n_missing > 0) {
            rec <- c(rec, "Some group anchors had missing GroupValue; default 0 was applied using the legacy-compatible group-centering rule.")
        }
        if (length(n_dup_anchor) > 0 && n_dup_anchor > 0) {
            rec <- c(rec, "Duplicate anchors were detected; the last row per Facet-Level was retained.")
        }
        if (length(n_dup_group) > 0 && n_dup_group > 0) {
            rec <- c(rec, "A Facet-Level appeared in multiple groups; the last row per Facet-Level was retained.")
        }
        if (length(n_group_conf) > 0 && n_group_conf > 0) {
            rec <- c(rec, "Conflicting GroupValue settings were detected within the same Facet-Group; the most recent finite value was retained.")
        }
    }
    if (!is.null(facet_summary) && nrow(facet_summary) > 0) {
        link_tbl <- dplyr::filter(facet_summary, Facet != "Person", AnchoredLevels > 0, AnchoredLevels < min_common_anchors)
        if (nrow(link_tbl) > 0) {
            rec <- c(rec, paste0("Package count screen: consider >= ", min_common_anchors, " directly anchored levels per facet. This does not verify cross-run identity or empirical connectedness. Low-count facets: ",
                paste(link_tbl$Facet, collapse = ", "), "."))
        }
        fixed_tbl <- dplyr::filter(facet_summary, FreeLevels <= 0, !Facet %in% dummy_facets)
        if (nrow(fixed_tbl) > 0) {
            rec <- c(rec, paste0("All levels are covered by direct or group constraints in: ", paste(fixed_tbl$Facet,
                collapse = ", "), ". Group members can still vary subject to their mean constraint; coverage is not a count of free parameters."))
        }
    }
    if (!is.null(design_checks) && is.list(design_checks)) {
        if (!is.null(design_checks$low_observation_levels) && nrow(design_checks$low_observation_levels) >
            0) {
            low_facets <- unique(design_checks$low_observation_levels$Facet)
            rec <- c(rec, paste0("Package observation-count screen: about ", fmt_count(min_obs_per_element),
                " observations per element are used for review. Low-observation facets: ", paste(low_facets,
                  collapse = ", "), "."))
        }
        if (!is.null(design_checks$low_categories) && nrow(design_checks$low_categories) > 0) {
            cats <- paste(design_checks$low_categories$Category, collapse = ", ")
            rec <- c(rec, paste0("Package observation-count screen: about ", fmt_count(min_obs_per_category),
                " observations per rating category are used for review. Low categories: ", cats, "."))
        }
    }
    rec <- c(rec, "Supply anchors in internal logit coordinates, not Umean/Uscale display units. Use compatible score coding, orientation and reporting scales across linked analyses.",
        paste0("Current noncenter facet is '", noncenter_facet, "'. ", if (identical(centering_scope,
            "facet")) {
            "Other facets sum to zero across all levels, including anchors and group totals. The ungrouped remainder balances those contributions; a single remaining level is fixed by that balance."
        } else if (identical(centering_scope, "anchor_origin")) {
            "Direct or group anchors set their facet's origin, replacing its zero-sum constraint. Other facets are centered except the declared unanchored noncenter facet. This convention is available for ordinary additive RSM/PCM JML."
        } else {
            "Legacy remainder centering: in other facets, levels outside direct and group constraints sum to zero; a single remaining level is fixed at zero."
        }))
    unique(rec)
}

.mfrmr_empty_anchor_audit <- local({
    anchors <- tibble::tibble(Facet = character(), Level = character(), Anchor = numeric())
    groups <- tibble::tibble(Facet = character(), Level = character(), Group = character(), GroupValue = numeric())
    list(anchors = anchors, group_anchors = tibble::tibble(Facet = character(), Level = character(), Group = character(),
        GroupValue = logical()), issues = list(anchor_schema_mismatch = tibble::tibble(), group_anchor_schema_mismatch = tibble::tibble(),
        unknown_anchor_facets = anchors, unknown_anchor_levels = anchors, invalid_anchor_values = anchors,
        duplicate_anchors = tibble::tibble(Facet = character(), Level = character(), Rows = integer(), DistinctValues = integer(),
            Values = character()), unknown_group_facets = groups, unknown_group_levels = groups, invalid_group_labels = groups,
        duplicate_group_assignments = tibble::tibble(Facet = character(), Level = character(), Rows = integer(),
            DistinctGroups = integer(), Groups = character()), missing_group_values = tibble::tibble(Facet = character(),
            Group = character()), group_value_conflicts = tibble::tibble(Facet = character(), Group = character(),
            DistinctFiniteValues = integer(), FiniteValues = character()), overlap_anchor_group = tibble::tibble(Facet = character(),
            Level = character())))
})

.mfrmr_anchor_observation_counts <- function(value, weight) {
    levels <- sort(unique(value), method = "radix", na.last = TRUE)
    index <- match(value, levels)
    groups <- factor(index, levels = seq_along(levels))
    list(level = levels, raw = tabulate(index, nbins = length(levels)), weighted = unname(vapply(split(weight,
        groups), sum, numeric(1), na.rm = TRUE)))
}

audit_anchor_tables <- function(prep, anchor_df = NULL, group_anchor_df = NULL, min_common_anchors = 5L,
    min_obs_per_element = 30, min_obs_per_category = 10, noncenter_facet = "Person", dummy_facets = character(0),
    centering_scope = "remainder") {
    all_facets <- c("Person", prep$facet_names)
    level_df <- collect_anchor_levels(prep)
    no_anchors <- (is.null(anchor_df) || nrow(anchor_df) == 0L) && (is.null(group_anchor_df) || nrow(group_anchor_df) ==
        0L)
    if (no_anchors) {
        anchors_clean <- .mfrmr_empty_anchor_audit$anchors
        groups_clean <- .mfrmr_empty_anchor_audit$group_anchors
        issues <- .mfrmr_empty_anchor_audit$issues
        facet_summary <- dplyr::mutate(dplyr::summarize(dplyr::group_by(level_df, Facet), Levels = dplyr::n_distinct(Level), .groups = "drop"),
            AnchoredLevels = 0L, GroupedLevels = 0L, GroupCount = 0L, ConstrainedLevels = 0L, OverlapLevels = 0L)
    }
    else {
        valid_keys <- safe_join_key(level_df$Facet, level_df$Level)
        anchor_schema_issue <- if (!is.null(anchor_df) && nrow(anchor_df) > 0) {
            nm <- tolower(names(anchor_df))
            has_facet <- any(nm %in% c("facet", "facets"))
            has_level <- any(nm %in% c("level", "element", "label"))
            has_anchor <- any(nm %in% c("anchor", "value", "measure"))
            if (!(has_facet && has_level && has_anchor)) {
                tibble::tibble(Columns = paste(names(anchor_df), collapse = ", "), Required = "facet + level + anchor/value/measure")
            }
            else {
                tibble::tibble()
            }
        }
        else {
            tibble::tibble()
        }
        group_schema_issue <- if (!is.null(group_anchor_df) && nrow(group_anchor_df) > 0) {
            nm <- tolower(names(group_anchor_df))
            has_facet <- any(nm %in% c("facet", "facets"))
            has_level <- any(nm %in% c("level", "element", "label"))
            has_group <- any(nm %in% c("group", "subset"))
            has_value <- any(nm %in% c("groupvalue", "value", "anchor"))
            if (!(has_facet && has_level && has_group && has_value)) {
                tibble::tibble(Columns = paste(names(group_anchor_df), collapse = ", "), Required = "facet + level + group/subset + groupvalue/value/anchor")
            }
            else {
                tibble::tibble()
            }
        }
        else {
            tibble::tibble()
        }
        anchor_in <- dplyr::mutate(normalize_anchor_df(anchor_df), .Row = dplyr::row_number(), Facet = trimws(as.character(Facet)),
            Level = trimws(as.character(Level)), .Key = safe_join_key(Facet, Level), .ValidFacet = Facet %in%
                all_facets, .ValidLevel = .Key %in% valid_keys, .ValidValue = is.finite(Anchor))
        group_in <- dplyr::mutate(normalize_group_anchor_df(group_anchor_df), .Row = dplyr::row_number(), Facet = trimws(as.character(Facet)),
            Level = trimws(as.character(Level)), Group = trimws(as.character(Group)), .Key = safe_join_key(Facet,
                Level), .ValidFacet = Facet %in% all_facets, .ValidLevel = .Key %in% valid_keys, .ValidGroup = nzchar(Group),
            .FiniteGroupValue = is.finite(GroupValue))
        issues <- list(anchor_schema_mismatch = anchor_schema_issue, group_anchor_schema_mismatch = group_schema_issue,
            unknown_anchor_facets = dplyr::select(dplyr::filter(anchor_in, !.ValidFacet), Facet, Level, Anchor), unknown_anchor_levels = dplyr::select(dplyr::filter(anchor_in,
                .ValidFacet, !.ValidLevel), Facet, Level, Anchor), invalid_anchor_values = dplyr::select(dplyr::filter(anchor_in,
                .ValidFacet, .ValidLevel, !.ValidValue), Facet, Level, Anchor), duplicate_anchors = dplyr::filter(dplyr::summarize(dplyr::group_by(dplyr::filter(anchor_in,
                .ValidFacet, .ValidLevel, .ValidValue), Facet, Level), Rows = dplyr::n(), DistinctValues = dplyr::n_distinct(Anchor),
                Values = paste(unique(round(Anchor, 6)), collapse = ", "), .groups = "drop"), Rows >
                1 | DistinctValues > 1), unknown_group_facets = dplyr::select(dplyr::filter(group_in, !.ValidFacet),
                Facet, Level, Group, GroupValue), unknown_group_levels = dplyr::select(dplyr::filter(group_in, .ValidFacet,
                !.ValidLevel), Facet, Level, Group, GroupValue), invalid_group_labels = dplyr::select(dplyr::filter(group_in,
                .ValidFacet, .ValidLevel, !.ValidGroup), Facet, Level, Group, GroupValue), duplicate_group_assignments = dplyr::filter(dplyr::summarize(dplyr::group_by(dplyr::filter(group_in,
                .ValidFacet, .ValidLevel, .ValidGroup), Facet, Level), Rows = dplyr::n(), DistinctGroups = dplyr::n_distinct(Group),
                Groups = paste(unique(Group), collapse = ", "), .groups = "drop"), Rows > 1 | DistinctGroups >
                1))
        anchors_clean <- dplyr::select(dplyr::ungroup(dplyr::slice_tail(dplyr::group_by(dplyr::arrange(dplyr::filter(anchor_in, .ValidFacet, .ValidLevel,
            .ValidValue), .Row), Facet, Level), n = 1)), Facet, Level, Anchor)
        groups_clean <- dplyr::select(dplyr::ungroup(dplyr::slice_tail(dplyr::group_by(dplyr::arrange(dplyr::filter(group_in, .ValidFacet, .ValidLevel,
            .ValidGroup), .Row), Facet, Level), n = 1)), Facet, Level, Group, GroupValue, .FiniteGroupValue)
        group_value_tbl <- dplyr::summarize(dplyr::group_by(dplyr::arrange(groups_clean, Facet, Group), Facet, Group), .NFinite = sum(.FiniteGroupValue,
            na.rm = TRUE), ChosenGroupValue = if (any(.FiniteGroupValue))
            dplyr::last(GroupValue[.FiniteGroupValue])
        else 0, DistinctFiniteValues = dplyr::n_distinct(GroupValue[.FiniteGroupValue]), FiniteValues = paste(unique(round(GroupValue[.FiniteGroupValue],
            6)), collapse = ", "), .groups = "drop")
        issues$missing_group_values <- dplyr::select(dplyr::filter(group_value_tbl, .NFinite == 0), Facet, Group)
        issues$group_value_conflicts <- dplyr::select(dplyr::filter(group_value_tbl, DistinctFiniteValues > 1), Facet,
            Group, DistinctFiniteValues, FiniteValues)
        groups_clean <- dplyr::mutate(dplyr::rename(dplyr::left_join(dplyr::select(groups_clean, Facet, Level, Group), dplyr::select(group_value_tbl,
            Facet, Group, ChosenGroupValue), by = c("Facet", "Group")), GroupValue = ChosenGroupValue),
            GroupValue = ifelse(is.finite(GroupValue), GroupValue, 0))
        overlap_tbl <- dplyr::inner_join(dplyr::select(anchors_clean, Facet, Level), dplyr::select(groups_clean, Facet, Level),
            by = c("Facet", "Level"))
        issues$overlap_anchor_group <- overlap_tbl
        constrained_counts <- dplyr::summarize(dplyr::group_by(dplyr::distinct(dplyr::bind_rows(dplyr::select(anchors_clean, Facet, Level),
            dplyr::select(groups_clean, Facet, Level)), Facet, Level), Facet), ConstrainedLevels = dplyr::n_distinct(Level),
            .groups = "drop")
        facet_counts <- dplyr::summarize(dplyr::group_by(level_df, Facet), Levels = dplyr::n_distinct(Level), .groups = "drop")
        anchor_counts <- dplyr::summarize(dplyr::group_by(anchors_clean, Facet), AnchoredLevels = dplyr::n_distinct(Level),
            .groups = "drop")
        group_counts <- dplyr::summarize(dplyr::group_by(groups_clean, Facet), GroupedLevels = dplyr::n_distinct(Level), GroupCount = dplyr::n_distinct(Group),
            .groups = "drop")
        overlap_counts <- dplyr::summarize(dplyr::group_by(overlap_tbl, Facet), OverlapLevels = dplyr::n_distinct(Level),
            .groups = "drop")
        facet_summary <- dplyr::left_join(dplyr::left_join(dplyr::left_join(dplyr::left_join(facet_counts, anchor_counts, by = "Facet"),
            group_counts, by = "Facet"), constrained_counts, by = "Facet"), overlap_counts, by = "Facet")
    }
    facet_summary <- dplyr::arrange(dplyr::mutate(facet_summary, AnchoredLevels = tidyr::replace_na(AnchoredLevels,
        0L), GroupedLevels = tidyr::replace_na(GroupedLevels, 0L), GroupCount = tidyr::replace_na(GroupCount,
        0L), ConstrainedLevels = tidyr::replace_na(ConstrainedLevels, 0L), OverlapLevels = tidyr::replace_na(OverlapLevels,
        0L), FreeLevels = pmax(Levels - ConstrainedLevels, 0L), Noncenter = Facet == noncenter_facet,
        DummyFacet = Facet %in% dummy_facets, Centered = !Noncenter & !DummyFacet & !(identical(centering_scope,
            "anchor_origin") & ConstrainedLevels > 0L), OriginBasis = dplyr::case_when(DummyFacet ~ "fixed_zero",
            identical(centering_scope, "anchor_origin") & AnchoredLevels > 0L & GroupedLevels > 0L ~
                "direct_and_group_anchors", identical(centering_scope, "anchor_origin") & AnchoredLevels >
                0L ~ "direct_anchors", identical(centering_scope, "anchor_origin") & GroupedLevels >
                0L ~ "group_anchors", Noncenter & AnchoredLevels > 0L ~ "direct_anchors", Noncenter &
                GroupedLevels > 0L ~ "group_anchors", Noncenter ~ "estimated_origin", identical(centering_scope,
                "remainder") ~ "remainder_zero_sum", TRUE ~ "zero_sum")), match(Facet, all_facets))
    design_df <- dplyr::mutate(prep$data, Person = as.character(Person), Score = as.numeric(Score), Weight = as.numeric(Weight),
        Weight = ifelse(is.finite(Weight) & Weight > 0, Weight, 0))
    level_obs <- dplyr::bind_rows(lapply(all_facets, function(facet) {
        if (!facet %in% names(design_df))
            return(tibble::tibble())
        counts <- .mfrmr_anchor_observation_counts(as.character(design_df[[facet]]), design_df$Weight)
        tibble::tibble(Facet = rep(facet, length(counts$level)), Level = counts$level, RawN = counts$raw, WeightedN = counts$weighted)
    }))
    level_obs_summary <- if (nrow(level_obs) == 0) {
        tibble::tibble()
    }
    else {
        dplyr::summarize(dplyr::group_by(level_obs, Facet), Levels = dplyr::n(), MinObsPerLevel = min(WeightedN, na.rm = TRUE),
            MedianObsPerLevel = stats::median(WeightedN, na.rm = TRUE), RecommendedMinObs = as.numeric(min_obs_per_element),
            PassMinObs = all(WeightedN >= min_obs_per_element), .groups = "drop")
    }
    low_observation_levels <- if (nrow(level_obs) == 0) {
        tibble::tibble()
    }
    else {
        dplyr::arrange(dplyr::filter(level_obs, WeightedN < min_obs_per_element), Facet, WeightedN, Level)
    }
    counts <- .mfrmr_anchor_observation_counts(design_df$Score, design_df$Weight)
    category_counts <- tibble::tibble(Category = counts$level, RawN = counts$raw, WeightedN = counts$weighted,
        RecommendedMinObs = rep(as.numeric(min_obs_per_category), length(counts$level)), PassMinObs = counts$weighted >=
            min_obs_per_category)
    low_categories <- dplyr::filter(category_counts, !PassMinObs)
    design_checks <- list(level_observation_summary = level_obs_summary, low_observation_levels = low_observation_levels,
        category_counts = category_counts, low_categories = low_categories)
    issue_counts <- build_anchor_issue_counts(issues)
    rec <- build_anchor_recommendations(facet_summary = facet_summary, issue_counts = issue_counts, design_checks = design_checks,
        min_common_anchors = min_common_anchors, min_obs_per_element = min_obs_per_element, min_obs_per_category = min_obs_per_category,
        noncenter_facet = noncenter_facet, dummy_facets = dummy_facets, centering_scope = centering_scope)
    thresholds <- list(min_common_anchors = as.integer(min_common_anchors), min_obs_per_element = as.numeric(min_obs_per_element),
        min_obs_per_category = as.numeric(min_obs_per_category))
    list(anchors = anchors_clean, group_anchors = groups_clean, centering_scope = centering_scope, facet_summary = facet_summary,
        design_checks = design_checks, thresholds = thresholds, issues = issues, issue_counts = issue_counts,
        recommendations = rec)
}
