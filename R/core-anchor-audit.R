# ==============================================================================
# Anchor table audit and table-input helpers
# ==============================================================================
#
# Internal helpers for reading anchor / group-anchor tables, normalizing
# their structure, and generating the connectivity / overlap audit
# bundle consumed by `fit_mfrm()`. Split out of `mfrm_core.R` for
# 0.1.6 so the anchor-table layer lives in a single file. All
# functions here are internal (no @export); they are called from the
# preparation phase of `mfrm_estimate()` before the optimizer runs.

read_flexible_table <- function(text_value, file_input) {
  if (!is.null(file_input) && !is.null(file_input$datapath)) {
    sep <- ifelse(grepl("\\.tsv$|\\.txt$", file_input$name, ignore.case = TRUE), "\t", ",")
    return(read.csv(file_input$datapath,
                    sep = sep,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE))
  }
  if (is.null(text_value)) return(tibble())
  text_value <- trimws(text_value)
  if (!nzchar(text_value)) return(tibble())
  sep <- if (grepl("\t", text_value)) {
    "\t"
  } else if (grepl(";", text_value)) {
    ";"
  } else {
    ","
  }
  read.csv(text = text_value,
           sep = sep,
           header = TRUE,
           stringsAsFactors = FALSE,
           check.names = FALSE)
}

normalize_anchor_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(Facet = character(0), Level = character(0), Anchor = numeric(0)))
  }
  nm <- tolower(names(df))
  facet_col <- which(nm %in% c("facet", "facets"))
  level_col <- which(nm %in% c("level", "element", "label"))
  anchor_col <- which(nm %in% c("anchor", "value", "measure"))
  if (length(facet_col) == 0 || length(level_col) == 0 || length(anchor_col) == 0) {
    return(tibble(Facet = character(0), Level = character(0), Anchor = numeric(0)))
  }
  tibble(
    Facet = as.character(df[[facet_col[1]]]),
    Level = as.character(df[[level_col[1]]]),
    Anchor = suppressWarnings(as.numeric(df[[anchor_col[1]]]))
  ) |>
    filter(!is.na(Facet), !is.na(Level))
}

normalize_group_anchor_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(Facet = character(0), Level = character(0), Group = character(0), GroupValue = numeric(0)))
  }
  nm <- tolower(names(df))
  facet_col <- which(nm %in% c("facet", "facets"))
  level_col <- which(nm %in% c("level", "element", "label"))
  group_col <- which(nm %in% c("group", "subset"))
  value_col <- which(nm %in% c("groupvalue", "value", "anchor"))
  if (length(facet_col) == 0 || length(level_col) == 0 || length(group_col) == 0 || length(value_col) == 0) {
    return(tibble(Facet = character(0), Level = character(0), Group = character(0), GroupValue = numeric(0)))
  }
  tibble(
    Facet = as.character(df[[facet_col[1]]]),
    Level = as.character(df[[level_col[1]]]),
    Group = as.character(df[[group_col[1]]]),
    GroupValue = suppressWarnings(as.numeric(df[[value_col[1]]]))
  ) |>
    filter(!is.na(Facet), !is.na(Level), !is.na(Group))
}

collect_anchor_levels <- function(prep) {
  facets <- c("Person", prep$facet_names)
  rows <- lapply(facets, function(facet) {
    lv <- prep$levels[[facet]]
    if (is.null(lv) || length(lv) == 0) return(tibble())
    tibble(
      Facet = as.character(facet),
      Level = as.character(lv)
    )
  })
  bind_rows(rows)
}

safe_join_key <- function(facet, level) {
  paste0(as.character(facet), "\r", as.character(level))
}

build_anchor_issue_counts <- function(issue_tables) {
  issue_names <- names(issue_tables)
  tibble(
    Issue = issue_names,
    N = vapply(issue_tables, nrow, integer(1))
  )
}

build_anchor_recommendations <- function(facet_summary,
                                         issue_counts,
                                         design_checks = NULL,
                                         min_common_anchors = 5L,
                                         min_obs_per_element = 30,
                                         min_obs_per_category = 10,
                                         noncenter_facet = "Person",
                                         dummy_facets = character(0)) {
  rec <- character(0)

  if (!is.null(issue_counts) && nrow(issue_counts) > 0) {
    n_overlap <- issue_counts$N[issue_counts$Issue == "overlap_anchor_group"]
    n_missing <- issue_counts$N[issue_counts$Issue == "missing_group_values"]
    n_dup_anchor <- issue_counts$N[issue_counts$Issue == "duplicate_anchors"]
    n_dup_group <- issue_counts$N[issue_counts$Issue == "duplicate_group_assignments"]
    n_group_conf <- issue_counts$N[issue_counts$Issue == "group_value_conflicts"]

    if (length(n_overlap) > 0 && n_overlap > 0) {
      rec <- c(rec, "Levels listed in both anchor and group-anchor tables are directly anchored (fixed anchors take precedence).")
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
    link_tbl <- facet_summary |>
      filter(Facet != "Person", AnchoredLevels > 0, AnchoredLevels < min_common_anchors)
    if (nrow(link_tbl) > 0) {
      rec <- c(
        rec,
        paste0(
          "FACETS linking guideline: consider >= ", min_common_anchors,
          " common anchor levels per linking facet. Low-anchor facets: ",
          paste(link_tbl$Facet, collapse = ", "), "."
        )
      )
    }

    fixed_tbl <- facet_summary |>
      filter(FreeLevels <= 0, !Facet %in% dummy_facets)
    if (nrow(fixed_tbl) > 0) {
      rec <- c(
        rec,
        paste0(
          "Some facets are fully constrained (no free levels): ",
          paste(fixed_tbl$Facet, collapse = ", "),
          ". Verify this is intentional."
        )
      )
    }
  }

  if (!is.null(design_checks) && is.list(design_checks)) {
    if (!is.null(design_checks$low_observation_levels) && nrow(design_checks$low_observation_levels) > 0) {
      low_facets <- unique(design_checks$low_observation_levels$Facet)
      rec <- c(
        rec,
        paste0(
          "Linacre guideline: about ", fmt_count(min_obs_per_element),
          " observations per element are desirable. Low-observation facets: ",
          paste(low_facets, collapse = ", "), "."
        )
      )
    }
    if (!is.null(design_checks$low_categories) && nrow(design_checks$low_categories) > 0) {
      cats <- paste(design_checks$low_categories$Category, collapse = ", ")
      rec <- c(
        rec,
        paste0(
          "Linacre guideline: about ", fmt_count(min_obs_per_category),
          " observations per rating category are desirable. Low categories: ", cats, "."
        )
      )
    }
  }

  rec <- c(
    rec,
    "For linked analyses, keep Umean/Uscale from the source calibration so reporting origin and scaling stay consistent.",
    paste0("Current noncenter facet is '", noncenter_facet, "'. Other facets are centered unless constrained by anchors/group anchors.")
  )

  unique(rec)
}

audit_anchor_tables <- function(prep,
                                anchor_df = NULL,
                                group_anchor_df = NULL,
                                min_common_anchors = 5L,
                                min_obs_per_element = 30,
                                min_obs_per_category = 10,
                                noncenter_facet = "Person",
                                dummy_facets = character(0)) {
  all_facets <- c("Person", prep$facet_names)
  level_df <- collect_anchor_levels(prep)
  valid_keys <- safe_join_key(level_df$Facet, level_df$Level)
  anchor_schema_issue <- if (!is.null(anchor_df) && nrow(anchor_df) > 0) {
    nm <- tolower(names(anchor_df))
    has_facet <- any(nm %in% c("facet", "facets"))
    has_level <- any(nm %in% c("level", "element", "label"))
    has_anchor <- any(nm %in% c("anchor", "value", "measure"))
    if (!(has_facet && has_level && has_anchor)) {
      tibble::tibble(
        Columns = paste(names(anchor_df), collapse = ", "),
        Required = "facet + level + anchor/value/measure"
      )
    } else {
      tibble::tibble()
    }
  } else {
    tibble::tibble()
  }
  group_schema_issue <- if (!is.null(group_anchor_df) && nrow(group_anchor_df) > 0) {
    nm <- tolower(names(group_anchor_df))
    has_facet <- any(nm %in% c("facet", "facets"))
    has_level <- any(nm %in% c("level", "element", "label"))
    has_group <- any(nm %in% c("group", "subset"))
    has_value <- any(nm %in% c("groupvalue", "value", "anchor"))
    if (!(has_facet && has_level && has_group && has_value)) {
      tibble::tibble(
        Columns = paste(names(group_anchor_df), collapse = ", "),
        Required = "facet + level + group/subset + groupvalue/value/anchor"
      )
    } else {
      tibble::tibble()
    }
  } else {
    tibble::tibble()
  }

  anchor_in <- normalize_anchor_df(anchor_df) |>
    mutate(
      .Row = row_number(),
      Facet = trimws(as.character(Facet)),
      Level = trimws(as.character(Level)),
      .Key = safe_join_key(Facet, Level),
      .ValidFacet = Facet %in% all_facets,
      .ValidLevel = .Key %in% valid_keys,
      .ValidValue = is.finite(Anchor)
    )

  group_in <- normalize_group_anchor_df(group_anchor_df) |>
    mutate(
      .Row = row_number(),
      Facet = trimws(as.character(Facet)),
      Level = trimws(as.character(Level)),
      Group = trimws(as.character(Group)),
      .Key = safe_join_key(Facet, Level),
      .ValidFacet = Facet %in% all_facets,
      .ValidLevel = .Key %in% valid_keys,
      .ValidGroup = nzchar(Group),
      .FiniteGroupValue = is.finite(GroupValue)
    )

  issues <- list(
    anchor_schema_mismatch = anchor_schema_issue,
    group_anchor_schema_mismatch = group_schema_issue,
    unknown_anchor_facets = anchor_in |>
      filter(!.ValidFacet) |>
      select(Facet, Level, Anchor),
    unknown_anchor_levels = anchor_in |>
      filter(.ValidFacet, !.ValidLevel) |>
      select(Facet, Level, Anchor),
    invalid_anchor_values = anchor_in |>
      filter(.ValidFacet, .ValidLevel, !.ValidValue) |>
      select(Facet, Level, Anchor),
    duplicate_anchors = anchor_in |>
      filter(.ValidFacet, .ValidLevel, .ValidValue) |>
      group_by(Facet, Level) |>
      summarize(
        Rows = n(),
        DistinctValues = n_distinct(Anchor),
        Values = paste(unique(round(Anchor, 6)), collapse = ", "),
        .groups = "drop"
      ) |>
      filter(Rows > 1 | DistinctValues > 1),
    unknown_group_facets = group_in |>
      filter(!.ValidFacet) |>
      select(Facet, Level, Group, GroupValue),
    unknown_group_levels = group_in |>
      filter(.ValidFacet, !.ValidLevel) |>
      select(Facet, Level, Group, GroupValue),
    invalid_group_labels = group_in |>
      filter(.ValidFacet, .ValidLevel, !.ValidGroup) |>
      select(Facet, Level, Group, GroupValue),
    duplicate_group_assignments = group_in |>
      filter(.ValidFacet, .ValidLevel, .ValidGroup) |>
      group_by(Facet, Level) |>
      summarize(
        Rows = n(),
        DistinctGroups = n_distinct(Group),
        Groups = paste(unique(Group), collapse = ", "),
        .groups = "drop"
      ) |>
      filter(Rows > 1 | DistinctGroups > 1)
  )

  anchors_clean <- anchor_in |>
    filter(.ValidFacet, .ValidLevel, .ValidValue) |>
    arrange(.Row) |>
    group_by(Facet, Level) |>
    slice_tail(n = 1) |>
    ungroup() |>
    select(Facet, Level, Anchor)

  groups_clean <- group_in |>
    filter(.ValidFacet, .ValidLevel, .ValidGroup) |>
    arrange(.Row) |>
    group_by(Facet, Level) |>
    slice_tail(n = 1) |>
    ungroup() |>
    select(Facet, Level, Group, GroupValue, .FiniteGroupValue)

  group_value_tbl <- groups_clean |>
    arrange(Facet, Group) |>
    group_by(Facet, Group) |>
    summarize(
      .NFinite = sum(.FiniteGroupValue, na.rm = TRUE),
      ChosenGroupValue = if (any(.FiniteGroupValue)) dplyr::last(GroupValue[.FiniteGroupValue]) else 0,
      DistinctFiniteValues = n_distinct(GroupValue[.FiniteGroupValue]),
      FiniteValues = paste(unique(round(GroupValue[.FiniteGroupValue], 6)), collapse = ", "),
      .groups = "drop"
    )

  issues$missing_group_values <- group_value_tbl |>
    filter(.NFinite == 0) |>
    select(Facet, Group)

  issues$group_value_conflicts <- group_value_tbl |>
    filter(DistinctFiniteValues > 1) |>
    select(Facet, Group, DistinctFiniteValues, FiniteValues)

  groups_clean <- groups_clean |>
    select(Facet, Level, Group) |>
    left_join(group_value_tbl |> select(Facet, Group, ChosenGroupValue), by = c("Facet", "Group")) |>
    rename(GroupValue = ChosenGroupValue) |>
    mutate(GroupValue = ifelse(is.finite(GroupValue), GroupValue, 0))

  overlap_tbl <- inner_join(
    anchors_clean |> select(Facet, Level),
    groups_clean |> select(Facet, Level),
    by = c("Facet", "Level")
  )
  issues$overlap_anchor_group <- overlap_tbl

  constrained_counts <- bind_rows(
    anchors_clean |> select(Facet, Level),
    groups_clean |> select(Facet, Level)
  ) |>
    distinct(Facet, Level) |>
    group_by(Facet) |>
    summarize(ConstrainedLevels = n_distinct(Level), .groups = "drop")

  facet_counts <- level_df |>
    group_by(Facet) |>
    summarize(Levels = n_distinct(Level), .groups = "drop")

  anchor_counts <- anchors_clean |>
    group_by(Facet) |>
    summarize(AnchoredLevels = n_distinct(Level), .groups = "drop")

  group_counts <- groups_clean |>
    group_by(Facet) |>
    summarize(GroupedLevels = n_distinct(Level), GroupCount = n_distinct(Group), .groups = "drop")

  overlap_counts <- overlap_tbl |>
    group_by(Facet) |>
    summarize(OverlapLevels = n_distinct(Level), .groups = "drop")

  facet_summary <- facet_counts |>
    left_join(anchor_counts, by = "Facet") |>
    left_join(group_counts, by = "Facet") |>
    left_join(constrained_counts, by = "Facet") |>
    left_join(overlap_counts, by = "Facet") |>
    mutate(
      AnchoredLevels = tidyr::replace_na(AnchoredLevels, 0L),
      GroupedLevels = tidyr::replace_na(GroupedLevels, 0L),
      GroupCount = tidyr::replace_na(GroupCount, 0L),
      ConstrainedLevels = tidyr::replace_na(ConstrainedLevels, 0L),
      OverlapLevels = tidyr::replace_na(OverlapLevels, 0L),
      FreeLevels = pmax(Levels - ConstrainedLevels, 0L),
      Noncenter = Facet == noncenter_facet,
      DummyFacet = Facet %in% dummy_facets
    ) |>
    arrange(match(Facet, all_facets))

  design_df <- prep$data |>
    mutate(
      Person = as.character(Person),
      Score = as.numeric(Score),
      Weight = as.numeric(Weight),
      Weight = ifelse(is.finite(Weight) & Weight > 0, Weight, 0)
    )

  level_obs <- bind_rows(lapply(all_facets, function(facet) {
    if (!facet %in% names(design_df)) return(tibble())
    design_df |>
      mutate(Level = as.character(.data[[facet]])) |>
      group_by(Facet = facet, Level) |>
      summarize(
        RawN = n(),
        WeightedN = sum(Weight, na.rm = TRUE),
        .groups = "drop"
      )
  }))

  level_obs_summary <- if (nrow(level_obs) == 0) {
    tibble()
  } else {
    level_obs |>
      group_by(Facet) |>
      summarize(
        Levels = n(),
        MinObsPerLevel = min(WeightedN, na.rm = TRUE),
        MedianObsPerLevel = stats::median(WeightedN, na.rm = TRUE),
        RecommendedMinObs = as.numeric(min_obs_per_element),
        PassMinObs = all(WeightedN >= min_obs_per_element),
        .groups = "drop"
      )
  }

  low_observation_levels <- if (nrow(level_obs) == 0) {
    tibble()
  } else {
    level_obs |>
      filter(WeightedN < min_obs_per_element) |>
      arrange(Facet, WeightedN, Level)
  }

  category_counts <- design_df |>
    group_by(Category = Score) |>
    summarize(
      RawN = n(),
      WeightedN = sum(Weight, na.rm = TRUE),
      RecommendedMinObs = as.numeric(min_obs_per_category),
      PassMinObs = WeightedN >= min_obs_per_category,
      .groups = "drop"
    ) |>
    arrange(Category)

  low_categories <- category_counts |>
    filter(!PassMinObs)

  design_checks <- list(
    level_observation_summary = level_obs_summary,
    low_observation_levels = low_observation_levels,
    category_counts = category_counts,
    low_categories = low_categories
  )

  issue_counts <- build_anchor_issue_counts(issues)
  rec <- build_anchor_recommendations(
    facet_summary = facet_summary,
    issue_counts = issue_counts,
    design_checks = design_checks,
    min_common_anchors = min_common_anchors,
    min_obs_per_element = min_obs_per_element,
    min_obs_per_category = min_obs_per_category,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets
  )

  thresholds <- list(
    min_common_anchors = as.integer(min_common_anchors),
    min_obs_per_element = as.numeric(min_obs_per_element),
    min_obs_per_category = as.numeric(min_obs_per_category)
  )

  list(
    anchors = anchors_clean,
    group_anchors = groups_clean,
    facet_summary = facet_summary,
    design_checks = design_checks,
    thresholds = thresholds,
    issues = issues,
    issue_counts = issue_counts,
    recommendations = rec
  )
}
