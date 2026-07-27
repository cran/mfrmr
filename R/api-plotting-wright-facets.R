# FACETS Table 6-style Wright-map data and base-R renderer.
#
# The helpers in this file intentionally describe a visual layout.  They do
# not change, or claim numerical equivalence with, the calibration engine used
# to produce the measures in an mfrm_fit object.

match_wright_style <- function(x, renderer = NULL) {
  if (!is.null(renderer)) {
    renderer_value <- tolower(as.character(renderer[1]))
    renderer_choices <- c("native", "facets")
    if (!(renderer_value %in% renderer_choices)) {
      stop(
        sprintf(
          "`renderer = \"%s\"` is not recognized. Choose %s.",
          renderer_value,
          paste(shQuote(renderer_choices), collapse = " or ")
        ),
        call. = FALSE
      )
    }
    return(if (identical(renderer_value, "facets")) "facets_style" else "native")
  }
  value <- tolower(as.character(x[1] %||% "native"))
  choices <- c("native", "facets_style")
  if (!(value %in% choices)) {
    stop(
      sprintf(
        "`wright_style = \"%s\"` is not recognized. Choose %s.",
        value,
        paste(shQuote(choices), collapse = " or ")
      ),
      call. = FALSE
    )
  }
  value
}

# Half-point thresholds are expected-score crossings at the midpoint between
# adjacent retained score categories.  The status columns deliberately match
# the category-curve/report contract so a renderer never turns an out-of-range
# or ambiguous crossing into an apparently exact line.
build_half_point_threshold_table <- function(expected_df, categories) {
  expected_df <- as.data.frame(expected_df %||% data.frame(), stringsAsFactors = FALSE)
  cats_num <- suppressWarnings(as.numeric(as.character(categories)))
  if (nrow(expected_df) == 0L || length(cats_num) < 2L || !all(is.finite(cats_num))) {
    return(data.frame())
  }
  needed <- c("CurveGroup", "Theta", "ExpectedScore")
  if (!all(needed %in% names(expected_df))) return(data.frame())
  cats_num <- sort(cats_num)
  rows <- list()
  idx <- 1L
  for (group in unique(as.character(expected_df$CurveGroup))) {
    line <- expected_df[as.character(expected_df$CurveGroup) == group, , drop = FALSE]
    line <- line[order(suppressWarnings(as.numeric(line$Theta))), , drop = FALSE]
    theta <- suppressWarnings(as.numeric(line$Theta))
    expected <- suppressWarnings(as.numeric(line$ExpectedScore))
    ok <- is.finite(theta) & is.finite(expected)
    theta <- theta[ok]
    expected <- expected[ok]
    for (j in seq_len(length(cats_num) - 1L)) {
      target <- (cats_num[j] + cats_num[j + 1L]) / 2
      threshold <- NA_real_
      in_range <- FALSE
      crossing_count <- 0L
      if (length(theta) > 1L) {
        centered <- expected - target
        exact <- which(abs(centered) <= sqrt(.Machine$double.eps))
        if (length(exact) > 0L) {
          threshold <- theta[exact[1]]
          in_range <- TRUE
          crossing_count <- length(exact)
        } else {
          crossing <- which(centered[-length(centered)] * centered[-1L] <= 0)
          crossing_count <- length(crossing)
          if (length(crossing) > 0L) {
            i <- crossing[1]
            y1 <- expected[i]
            y2 <- expected[i + 1L]
            x1 <- theta[i]
            x2 <- theta[i + 1L]
            threshold <- if (abs(y2 - y1) > sqrt(.Machine$double.eps)) {
              x1 + (target - y1) * (x2 - x1) / (y2 - y1)
            } else {
              mean(c(x1, x2))
            }
            in_range <- TRUE
          }
        }
      }
      rows[[idx]] <- data.frame(
        CurveGroup = group,
        PairOrder = as.integer(j),
        LowerCategory = as.character(cats_num[j]),
        UpperCategory = as.character(cats_num[j + 1L]),
        TargetScore = target,
        HalfPointThreshold = threshold,
        InThetaRange = in_range,
        CrossingCount = as.integer(crossing_count),
        CrossingStatus = dplyr::case_when(
          !in_range ~ "outside_theta_range",
          crossing_count == 1L ~ "in_range",
          crossing_count > 1L ~ "multiple_crossings",
          TRUE ~ "review"
        ),
        CrossingLabel = paste0("E[X] = ", format(target)),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  if (length(rows) == 0L) return(data.frame())
  dplyr::bind_rows(rows) |>
    dplyr::arrange(.data$CurveGroup, .data$PairOrder)
}

validate_wright_range <- function(x) {
  if (is.null(x)) return(NULL)
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 2L || !all(is.finite(x)) || x[1] >= x[2]) {
    stop("`wright_range` must be NULL or a finite increasing length-2 numeric vector.", call. = FALSE)
  }
  x
}

wright_score_labels <- function(fit, categories, category_labels = NULL) {
  categories <- as.character(categories)
  score_map <- as.data.frame(fit$prep$score_map %||% data.frame(), stringsAsFactors = FALSE)
  original <- categories
  if (nrow(score_map) > 0L &&
      all(c("OriginalScore", "InternalScore") %in% names(score_map))) {
    idx <- match(categories, as.character(score_map$InternalScore))
    mapped <- as.character(score_map$OriginalScore[idx])
    original[!is.na(mapped)] <- mapped[!is.na(mapped)]
  }

  rubric <- rep(NA_character_, length(categories))
  if (!is.null(category_labels)) {
    if (is.data.frame(category_labels)) {
      label_tbl <- as.data.frame(category_labels, stringsAsFactors = FALSE)
      if (!all(c("Score", "Label") %in% names(label_tbl))) {
        stop("A data-frame `category_labels` must contain `Score` and `Label` columns.", call. = FALSE)
      }
      label_keys <- as.character(label_tbl$Score)
      if (anyNA(label_keys) || any(!nzchar(label_keys))) {
        stop("Every `Score` key in a data-frame `category_labels` must be non-missing and non-empty.", call. = FALSE)
      }
      if (anyDuplicated(label_keys)) {
        duplicated_keys <- unique(label_keys[duplicated(label_keys)])
        stop(
          "`category_labels$Score` contains duplicate key(s): ",
          paste(duplicated_keys, collapse = ", "),
          ". Supply one label per original score.",
          call. = FALSE
        )
      }
      missing_keys <- setdiff(unique(original), label_keys)
      if (length(missing_keys) > 0L) {
        stop(
          "`category_labels` does not cover retained original score key(s): ",
          paste(missing_keys, collapse = ", "),
          ". Key labels by the original scores in `fit$prep$score_map$OriginalScore`.",
          call. = FALSE
        )
      }
      rubric <- as.character(label_tbl$Label[match(original, label_keys)])
    } else {
      labels <- as.character(category_labels)
      label_names <- names(category_labels)
      if (!is.null(label_names) && any(nzchar(label_names))) {
        if (anyNA(label_names) || any(!nzchar(label_names))) {
          stop(
            "A named `category_labels` vector must name every label with its original score key.",
            call. = FALSE
          )
        }
        if (anyDuplicated(label_names)) {
          duplicated_keys <- unique(label_names[duplicated(label_names)])
          stop(
            "A named `category_labels` vector contains duplicate key(s): ",
            paste(duplicated_keys, collapse = ", "),
            ". Supply one label per original score.",
            call. = FALSE
          )
        }
        missing_keys <- setdiff(unique(original), label_names)
        if (length(missing_keys) > 0L) {
          internal_hint <- if (
            !identical(original, categories) &&
              all(categories %in% label_names)
          ) {
            paste0(
              " The supplied names appear to use internal score keys (",
              paste(categories, collapse = ", "),
              "); use the original keys shown in `fit$prep$score_map` instead."
            )
          } else {
            ""
          }
          stop(
            "Named `category_labels` does not cover retained original score key(s): ",
            paste(missing_keys, collapse = ", "),
            ".",
            internal_hint,
            call. = FALSE
          )
        }
        rubric <- labels[match(original, label_names)]
      } else if (length(labels) == length(categories)) {
        rubric <- labels
      } else {
        stop(
          "An unnamed `category_labels` vector must have one value per retained score category.",
          call. = FALSE
        )
      }
    }
  }
  rubric[is.na(rubric)] <- ""
  display <- ifelse(nzchar(rubric), paste0(original, " ", rubric), original)

  tibble::tibble(
    InternalScore = categories,
    OriginalScore = original,
    RubricLabel = rubric,
    DisplayLabel = display
  )
}

build_wright_score_rulers <- function(fit,
                                      theta_range,
                                      category_labels = NULL,
                                      include_steps = TRUE) {
  empty_steps <- tibble::tibble(
    CurveGroup = character(), Step = character(), StepIndex = integer(),
    Estimate = numeric(), LowerScore = character(), UpperScore = character(),
    LowerLabel = character(), UpperLabel = character(),
    TransitionLabel = character()
  )
  empty_half <- tibble::tibble(
    CurveGroup = character(), LowerScore = character(), UpperScore = character(),
    LowerLabel = character(), UpperLabel = character(),
    MeanHalfScore = numeric(), Estimate = numeric(), TransitionLabel = character(),
    PairOrder = integer(), InThetaRange = logical(), CrossingCount = integer(),
    CrossingStatus = character(), CrossingLabel = character()
  )
  if (!isTRUE(include_steps)) {
    return(list(
      category_labels = tibble::tibble(),
      step_ruler = empty_steps,
      score_transitions = empty_half
    ))
  }

  curve_spec <- tryCatch(build_step_curve_spec(fit), error = function(e) NULL)
  if (is.null(curve_spec)) {
    return(list(
      category_labels = tibble::tibble(),
      step_ruler = empty_steps,
      score_transitions = empty_half
    ))
  }
  label_tbl <- wright_score_labels(
    fit,
    categories = curve_spec$categories,
    category_labels = category_labels
  )
  step_points <- tibble::as_tibble(curve_spec$step_points)
  step_ruler <- if (nrow(step_points) > 0L) {
    step_points |>
      dplyr::mutate(
        StepIndex = as.integer(.data$StepIndex),
        Estimate = as.numeric(.data$Threshold),
        LowerScore = label_tbl$OriginalScore[pmax(1L, .data$StepIndex)],
        UpperScore = label_tbl$OriginalScore[pmin(nrow(label_tbl), .data$StepIndex + 1L)],
        LowerLabel = label_tbl$DisplayLabel[pmax(1L, .data$StepIndex)],
        UpperLabel = label_tbl$DisplayLabel[pmin(nrow(label_tbl), .data$StepIndex + 1L)],
        TransitionLabel = paste(.data$LowerLabel, .data$UpperLabel, sep = " -> ")
      ) |>
      dplyr::select(
        "CurveGroup", "Step", "StepIndex", "Estimate", "LowerScore", "UpperScore",
        "LowerLabel", "UpperLabel", "TransitionLabel"
      )
  } else {
    empty_steps
  }

  theta_range <- suppressWarnings(as.numeric(theta_range))
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = 1601L)
  curve_tbl <- build_curve_tables(curve_spec, theta_grid)
  exp_df <- as.data.frame(curve_tbl$expected, stringsAsFactors = FALSE)
  half_points <- build_half_point_threshold_table(exp_df, curve_spec$categories)
  score_transitions <- if (nrow(half_points) > 0L) {
    tibble::as_tibble(half_points) |>
      dplyr::mutate(
        LowerScore = label_tbl$OriginalScore[.data$PairOrder],
        UpperScore = label_tbl$OriginalScore[.data$PairOrder + 1L],
        LowerLabel = label_tbl$DisplayLabel[.data$PairOrder],
        UpperLabel = label_tbl$DisplayLabel[.data$PairOrder + 1L],
        MeanHalfScore = .data$TargetScore,
        Estimate = .data$HalfPointThreshold,
        TransitionLabel = paste(.data$LowerLabel, .data$UpperLabel, sep = " -> ")
      ) |>
      dplyr::select(
        "CurveGroup", "LowerScore", "UpperScore", "LowerLabel", "UpperLabel",
        "MeanHalfScore", "Estimate", "TransitionLabel", "PairOrder",
        "InThetaRange", "CrossingCount", "CrossingStatus", "CrossingLabel"
      )
  } else {
    empty_half
  }

  list(
    category_labels = label_tbl,
    step_ruler = tibble::as_tibble(step_ruler),
    score_transitions = tibble::as_tibble(score_transitions)
  )
}

nearest_wright_row <- function(x, lower, upper, rows_per_logit) {
  step <- 1 / rows_per_logit
  pmin(upper, pmax(lower, round((x - lower) / step) * step + lower))
}

build_wright_display_scale <- function(fit,
                                       person_values,
                                       location_tbl,
                                       category_labels = NULL,
                                       rows_per_logit = 2L,
                                       wright_range = NULL,
                                       include_steps = TRUE) {
  person_values <- suppressWarnings(as.numeric(person_values))
  location_tbl <- as.data.frame(location_tbl, stringsAsFactors = FALSE)
  location_values <- suppressWarnings(as.numeric(location_tbl$Estimate))
  base_values <- c(person_values, location_values)
  base_values <- base_values[is.finite(base_values)]
  if (length(base_values) == 0L) base_values <- c(-1, 1)

  boundary_levels <- as.data.frame(
    fit$data_review$boundary_levels %||% data.frame(),
    stringsAsFactors = FALSE
  )
  boundary_key <- character()
  if (nrow(boundary_levels) > 0L &&
      all(c("Facet", "Level") %in% names(boundary_levels))) {
    boundary_key <- unique(paste(
      as.character(boundary_levels$Facet),
      as.character(boundary_levels$Level),
      sep = "\r"
    ))
  }

  auto_range_policy <- "all_fitted_locations"
  range_basis <- base_values
  if (is.null(wright_range) && length(boundary_key) > 0L &&
      all(c("PlotType", "Group", "Label") %in% names(location_tbl))) {
    location_key <- paste(
      as.character(location_tbl$Group),
      as.character(location_tbl$Label),
      sep = "\r"
    )
    boundary_match <- location_tbl$PlotType == "Facet level" &
      location_key %in% boundary_key
    robust_basis <- c(person_values, location_values[!boundary_match])
    robust_basis <- robust_basis[is.finite(robust_basis)]
    if (any(boundary_match) && length(robust_basis) >= 2L) {
      range_basis <- robust_basis
      auto_range_policy <- "boundary_levels_at_ends"
    }
  }

  provisional <- c(floor(min(range_basis)), ceiling(max(range_basis)))
  if (diff(provisional) <= 0) provisional <- mean(provisional) + c(-1, 1)
  score_rulers <- build_wright_score_rulers(
    fit,
    theta_range = provisional + c(-6, 6),
    category_labels = category_labels,
    include_steps = include_steps
  )
  transition_values <- c(
    suppressWarnings(as.numeric(score_rulers$step_ruler$Estimate)),
    suppressWarnings(as.numeric(score_rulers$score_transitions$Estimate))
  )
  transition_values <- transition_values[is.finite(transition_values)]
  all_values <- c(range_basis, transition_values)
  if (is.null(wright_range)) {
    lower <- floor(min(all_values) * rows_per_logit) / rows_per_logit
    upper <- ceiling(max(all_values) * rows_per_logit) / rows_per_logit
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
      center <- mean(all_values, na.rm = TRUE)
      lower <- center - 1
      upper <- center + 1
    }
  } else {
    lower <- wright_range[1]
    upper <- wright_range[2]
  }

  list(
    lower = lower,
    upper = upper,
    auto_range_policy = auto_range_policy,
    boundary_key = boundary_key,
    score_rulers = score_rulers
  )
}

add_wright_ci_display_metadata <- function(tbl, lower, upper) {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  n <- nrow(tbl)
  numeric_column <- function(name) {
    if (name %in% names(tbl)) {
      suppressWarnings(as.numeric(tbl[[name]]))
    } else {
      rep(NA_real_, n)
    }
  }
  logical_column <- function(name) {
    if (name %in% names(tbl)) {
      value <- as.logical(tbl[[name]])
      value[is.na(value)] <- FALSE
      value
    } else {
      rep(FALSE, n)
    }
  }

  tbl$OriginalCI_Lower <- numeric_column("CI_Lower")
  tbl$OriginalCI_Upper <- numeric_column("CI_Upper")
  has_ci <- is.finite(tbl$OriginalCI_Lower) & is.finite(tbl$OriginalCI_Upper)
  tbl$DisplayCI_Lower <- ifelse(
    has_ci,
    pmin(upper, pmax(lower, tbl$OriginalCI_Lower)),
    NA_real_
  )
  tbl$DisplayCI_Upper <- ifelse(
    has_ci,
    pmin(upper, pmax(lower, tbl$OriginalCI_Upper)),
    NA_real_
  )
  tbl$CIClippedLower <- has_ci & tbl$OriginalCI_Lower < lower
  tbl$CIClippedUpper <- has_ci & tbl$OriginalCI_Upper > upper
  tbl$CIClipped <- tbl$CIClippedLower | tbl$CIClippedUpper

  boundary <- logical_column("BoundarySeparated")
  below <- logical_column("BelowRange")
  above <- logical_column("AboveRange")
  tbl$BoundaryEnd <- ifelse(
    boundary & below,
    "lower",
    ifelse(boundary & above, "upper", "none")
  )
  tbl$CISuppressed <- boundary & has_ci &
    (below | above | tbl$CIClipped)
  tbl$CIDisplayStatus <- ifelse(
    !has_ci,
    "not_available",
    ifelse(
      tbl$CISuppressed,
      "boundary_endpoint_marker",
      ifelse(tbl$CIClipped, "clipped_with_endpoint_marker", "full_interval")
    )
  )
  tbl
}

build_facets_style_wright_data <- function(fit,
                                           plot_data,
                                           category_labels = NULL,
                                           rows_per_logit = 2L,
                                           wright_range = NULL,
                                           extreme_placement = c("ends", "estimate"),
                                           persons_per_star = NULL,
                                           include_steps = TRUE) {
  rows_per_logit <- suppressWarnings(as.integer(rows_per_logit[1]))
  if (!is.finite(rows_per_logit) || rows_per_logit < 1L || rows_per_logit > 20L) {
    stop("`rows_per_logit` must be an integer from 1 through 20.", call. = FALSE)
  }
  wright_range <- validate_wright_range(wright_range)
  extreme_placement <- match.arg(tolower(as.character(extreme_placement[1])), c("ends", "estimate"))

  person_values <- suppressWarnings(as.numeric(plot_data$person$Estimate))
  location_tbl <- as.data.frame(plot_data$locations, stringsAsFactors = FALSE)
  display_scale <- build_wright_display_scale(
    fit = fit,
    person_values = person_values,
    location_tbl = location_tbl,
    category_labels = category_labels,
    rows_per_logit = rows_per_logit,
    wright_range = wright_range,
    include_steps = include_steps
  )
  lower <- display_scale$lower
  upper <- display_scale$upper
  auto_range_policy <- display_scale$auto_range_policy
  boundary_key <- display_scale$boundary_key
  score_rulers <- display_scale$score_rulers
  row_step <- 1 / rows_per_logit
  ruler_values <- seq(lower, upper + row_step / 10, by = row_step)
  ruler_values <- ruler_values[ruler_values <= upper + sqrt(.Machine$double.eps)]
  ruler_rows <- tibble::tibble(
    Row = seq_along(ruler_values),
    Logit = ruler_values,
    Major = abs(ruler_values - round(ruler_values)) < sqrt(.Machine$double.eps),
    Label = ifelse(
      abs(ruler_values - round(ruler_values)) < sqrt(.Machine$double.eps),
      format(round(ruler_values), trim = TRUE),
      format(round(ruler_values, 3), trim = TRUE)
    )
  )

  person <- as.data.frame(plot_data$person, stringsAsFactors = FALSE)
  extreme <- if ("Extreme" %in% names(person)) tolower(as.character(person$Extreme)) else rep("none", nrow(person))
  extreme[is.na(extreme) | !nzchar(extreme)] <- "none"
  person$OriginalEstimate <- suppressWarnings(as.numeric(person$Estimate))
  person$DisplayEstimate <- person$OriginalEstimate
  if (identical(extreme_placement, "ends")) {
    person$DisplayEstimate[extreme == "low"] <- lower
    person$DisplayEstimate[extreme == "high"] <- upper
  }
  person$BelowRange <- person$DisplayEstimate < lower
  person$AboveRange <- person$DisplayEstimate > upper
  person$DisplayEstimate <- pmin(upper, pmax(lower, person$DisplayEstimate))
  person$RulerValue <- nearest_wright_row(person$DisplayEstimate, lower, upper, rows_per_logit)
  person$ExtremePlacement <- extreme

  person_frequency <- person |>
    dplyr::filter(is.finite(.data$RulerValue)) |>
    dplyr::count(.data$RulerValue, name = "Count") |>
    dplyr::right_join(
      ruler_rows |>
        dplyr::transmute(RulerValue = .data$Logit),
      by = "RulerValue"
    ) |>
    dplyr::mutate(Count = dplyr::coalesce(.data$Count, 0L)) |>
    dplyr::arrange(.data$RulerValue)
  if (is.null(persons_per_star)) {
    persons_per_star <- max(1L, ceiling(max(person_frequency$Count, na.rm = TRUE) / 20))
  } else {
    persons_per_star <- suppressWarnings(as.numeric(persons_per_star[1]))
    if (!is.finite(persons_per_star) || persons_per_star <= 0) {
      stop("`persons_per_star` must be NULL or one positive number.", call. = FALSE)
    }
  }
  person_frequency$StarCount <- ifelse(
    person_frequency$Count > 0,
    pmax(1L, ceiling(person_frequency$Count / persons_per_star)),
    0L
  )
  person_frequency$Stars <- vapply(
    person_frequency$StarCount,
    function(n) if (n > 0L) paste(rep("*", n), collapse = "") else "",
    character(1)
  )

  facet_ruler <- as.data.frame(
    plot_data$locations[plot_data$locations$PlotType == "Facet level", , drop = FALSE],
    stringsAsFactors = FALSE
  )
  signs <- fit$config$facet_signs %||%
    stats::setNames(rep(-1, length(unique(facet_ruler$Group))), unique(facet_ruler$Group))
  sign_value <- suppressWarnings(as.numeric(signs[as.character(facet_ruler$Group)]))
  sign_value[!is.finite(sign_value)] <- -1
  facet_ruler$Facet <- as.character(facet_ruler$Group)
  facet_ruler$Level <- as.character(facet_ruler$Label)
  facet_ruler$BoundarySeparated <- paste(
    facet_ruler$Facet,
    facet_ruler$Level,
    sep = "\r"
  ) %in% boundary_key
  facet_ruler$Sign <- sign_value
  facet_ruler$Orientation <- ifelse(sign_value >= 0, "positive", "negative")
  facet_ruler$Header <- paste0(ifelse(sign_value >= 0, "+", "-"), facet_ruler$Facet)
  facet_ruler$OriginalEstimate <- suppressWarnings(as.numeric(facet_ruler$Estimate))
  facet_ruler$BelowRange <- facet_ruler$OriginalEstimate < lower
  facet_ruler$AboveRange <- facet_ruler$OriginalEstimate > upper
  facet_ruler$DisplayEstimate <- pmin(upper, pmax(lower, facet_ruler$OriginalEstimate))
  facet_ruler$RulerValue <- nearest_wright_row(facet_ruler$DisplayEstimate, lower, upper, rows_per_logit)
  facet_ruler$DisplayLabel <- ifelse(
    facet_ruler$BelowRange | facet_ruler$AboveRange,
    paste0("(", facet_ruler$Level, ")"),
    facet_ruler$Level
  )
  facet_ruler <- add_wright_ci_display_metadata(
    facet_ruler,
    lower = lower,
    upper = upper
  )
  facet_cells <- tibble::as_tibble(facet_ruler) |>
    dplyr::group_by(.data$Facet, .data$Header, .data$RulerValue) |>
    dplyr::summarise(
      CellLabel = paste(.data$DisplayLabel, collapse = ", "),
      N = dplyr::n(),
      .groups = "drop"
    )

  add_ruler_position <- function(tbl) {
    if (nrow(tbl) == 0L) return(tbl)
    tbl |>
      dplyr::mutate(
        InRange = is.finite(.data$Estimate) & .data$Estimate >= lower & .data$Estimate <= upper,
        DisplayEstimate = pmin(upper, pmax(lower, .data$Estimate)),
        RulerValue = nearest_wright_row(.data$DisplayEstimate, lower, upper, rows_per_logit),
        # FACETS-style frequency cells use the discrete ruler row, while
        # threshold lines retain their exact fitted logit coordinate.
        DrawValue = .data$DisplayEstimate
      )
  }
  step_ruler <- add_ruler_position(score_rulers$step_ruler)
  if (nrow(step_ruler) > 0L) {
    step_ruler$DrawLabel <- paste0(
      step_ruler$TransitionLabel,
      " [",
      formatC(step_ruler$Estimate, digits = 2L, format = "f"),
      "]"
    )
  }
  score_transitions <- add_ruler_position(score_rulers$score_transitions)
  step_headers <- if (nrow(step_ruler) > 0L) {
    tibble::tibble(
      CurveGroup = unique(as.character(step_ruler$CurveGroup)),
      Header = paste0("Scale:", unique(as.character(step_ruler$CurveGroup)))
    )
  } else {
    tibble::tibble(CurveGroup = character(), Header = character())
  }

  settings <- tibble::tibble(
    Renderer = "facets",
    WrightStyle = "facets_style",
    VisualCorrespondence = "FACETS Table 6-style layout; not a claim of FACETS numerical equivalence",
    LowerLogit = lower,
    UpperLogit = upper,
    RowsPerLogit = rows_per_logit,
    ExtremePlacement = extreme_placement,
    PersonsPerStar = persons_per_star,
    StarsPerPerson = 1 / persons_per_star,
    PersonN = nrow(person)
  )
  settings$AutoRangePolicy <- auto_range_policy
  settings$BoundaryLevelsAtEnds <- sum(facet_ruler$BoundaryEnd != "none")
  settings$CIClippedCount <- sum(facet_ruler$CIClipped)
  settings$BoundaryCIEndpointCount <- sum(facet_ruler$CISuppressed)
  settings$CIDisplayPolicy <- paste0(
    "Exact intervals remain in CI_Lower/CI_Upper; clipped intervals use ",
    "endpoint markers and boundary-separated intervals may be omitted from ",
    "the ruler."
  )
  headers <- dplyr::bind_rows(
    tibble::tibble(ColumnType = "person", Group = "Person", Header = "+Person"),
    unique(facet_ruler[, c("Facet", "Header"), drop = FALSE]) |>
      dplyr::transmute(ColumnType = "facet", Group = .data$Facet, Header = .data$Header),
    step_headers |>
      dplyr::transmute(ColumnType = "step", Group = .data$CurveGroup, Header = .data$Header)
  )

  list(
    ruler_rows = ruler_rows,
    person_frequency = tibble::as_tibble(person_frequency),
    person_placements = tibble::as_tibble(person),
    facet_ruler = tibble::as_tibble(facet_ruler),
    facet_cells = tibble::as_tibble(facet_cells),
    step_ruler = tibble::as_tibble(step_ruler),
    score_transitions = tibble::as_tibble(score_transitions),
    category_labels = tibble::as_tibble(score_rulers$category_labels),
    headers = tibble::as_tibble(headers),
    settings = settings
  )
}

draw_wright_facets_style <- function(plot_data,
                                     title = NULL,
                                     palette = NULL,
                                     show_ci = FALSE,
                                     ci_level = 0.95) {
  facets_data <- plot_data$facets_style
  if (is.null(facets_data) || !is.list(facets_data)) {
    stop("FACETS-style Wright-map data are missing.", call. = FALSE)
  }
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      facet_level = "#1b9e77",
      step_threshold = "#d95f02",
      person_hist = "gray25",
      person_star = "gray20",
      grid = "#e5e7eb",
      range = "#94a3b8",
      iqr = "#334155"
    )
  )
  headers <- as.data.frame(facets_data$headers, stringsAsFactors = FALSE)
  if (nrow(headers) == 0L) stop("No Wright-map columns are available.", call. = FALSE)
  n_col <- nrow(headers)
  settings <- facets_data$settings
  yr <- c(settings$LowerLogit[1], settings$UpperLogit[1])
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  auto_range_policy <- as.character(
    settings$AutoRangePolicy[1L] %||% "all_fitted_locations"
  )
  ci_clipped_count <- suppressWarnings(as.integer(
    settings$CIClippedCount[1L] %||% 0L
  ))
  boundary_ci_count <- suppressWarnings(as.integer(
    settings$BoundaryCIEndpointCount[1L] %||% 0L
  ))
  extra_footer <- isTRUE(show_ci) ||
    identical(auto_range_policy, "boundary_levels_at_ends")
  graphics::par(
    mar = c(if (extra_footer) 7 else 5.4, 4.8, 5.4, 2.2),
    mgp = c(2.6, 0.8, 0),
    xpd = FALSE
  )
  graphics::plot(
    NA_real_, NA_real_,
    xlim = c(0.45, n_col + 0.55),
    ylim = yr,
    xaxt = "n", yaxt = "n", xlab = "", ylab = "Logit ruler",
    main = title %||% "FACETS-style Wright map"
  )
  rows <- facets_data$ruler_rows
  major <- !is.na(rows$Major) & rows$Major
  if (nrow(rows) <= 30L) {
    axis_at <- rows$Logit
    axis_labels <- rows$Label
  } else {
    axis_at <- pretty(yr, n = 10L)
    axis_at <- axis_at[axis_at >= yr[1] & axis_at <= yr[2]]
    axis_labels <- format(axis_at, trim = TRUE)
  }
  grid_at <- if (nrow(rows) <= 30L) rows$Logit else axis_at
  graphics::abline(h = grid_at, col = grDevices::adjustcolor(pal["grid"], alpha.f = 0.38), lty = 1)
  if (nrow(rows) <= 30L) {
    graphics::abline(h = rows$Logit[major], col = grDevices::adjustcolor(pal["range"], alpha.f = 0.75), lwd = 1.1)
  }
  graphics::abline(v = seq(0.5, n_col + 0.5, by = 1), col = grDevices::adjustcolor(pal["grid"], alpha.f = 0.9))
  graphics::axis(2, at = axis_at, labels = axis_labels, las = 1, cex.axis = 0.78)
  graphics::axis(4, at = axis_at, labels = axis_labels, las = 1, cex.axis = 0.78)
  graphics::axis(3, at = seq_len(n_col), labels = headers$Header, tick = FALSE, line = 0.25, cex.axis = 0.78)

  person_x <- which(headers$ColumnType == "person")[1]
  person_freq <- facets_data$person_frequency
  if (is.finite(person_x) && nrow(person_freq) > 0L) {
    graphics::text(
      x = person_x,
      y = person_freq$RulerValue,
      labels = person_freq$Stars,
      cex = 0.78,
      family = "mono",
      col = pal["person_star"]
    )
  }

  facet_cells <- facets_data$facet_cells
  facet_ruler <- facets_data$facet_ruler
  facet_headers <- headers[headers$ColumnType == "facet", , drop = FALSE]
  spread_label_y <- function(y, lower, upper, min_separation) {
    y <- pmin(upper, pmax(lower, as.numeric(y)))
    finite <- which(is.finite(y))
    if (length(finite) <= 1L) return(y)
    ord <- finite[order(y[finite])]
    adjusted <- y
    for (j in 2:length(ord)) {
      adjusted[ord[j]] <- max(
        adjusted[ord[j]],
        adjusted[ord[j - 1L]] + min_separation
      )
    }
    overflow <- adjusted[ord[length(ord)]] - upper
    if (is.finite(overflow) && overflow > 0) adjusted[ord] <- adjusted[ord] - overflow
    underflow <- lower - adjusted[ord[1L]]
    if (is.finite(underflow) && underflow > 0) adjusted[ord] <- adjusted[ord] + underflow
    adjusted
  }
  for (i in seq_len(nrow(facet_headers))) {
    x_pos <- match(facet_headers$Header[i], headers$Header)
    cell <- facet_cells[facet_cells$Facet == facet_headers$Group[i], , drop = FALSE]
    if (nrow(cell) > 0L && !isTRUE(show_ci)) {
      cell_labels <- vapply(
        as.character(cell$CellLabel),
        function(value) paste(strwrap(value, width = 22L), collapse = "\n"),
        character(1)
      )
      graphics::text(
        x = x_pos,
        y = cell$RulerValue,
        labels = cell_labels,
        cex = 0.68,
        col = pal["facet_level"]
      )
    }
    if (isTRUE(show_ci) && all(c("SE", "Estimate") %in% names(facet_ruler))) {
      sub <- facet_ruler[facet_ruler$Facet == facet_headers$Group[i], , drop = FALSE]
      raw_label_y <- suppressWarnings(as.numeric(sub$DisplayEstimate))
      label_y <- spread_label_y(
        raw_label_y,
        lower = yr[1],
        upper = yr[2],
        min_separation = 0.04 * diff(yr)
      )
      moved <- is.finite(raw_label_y) & is.finite(label_y) &
        abs(raw_label_y - label_y) > sqrt(.Machine$double.eps)
      if (any(moved)) {
        graphics::segments(
          x0 = x_pos - 0.04,
          y0 = raw_label_y[moved],
          x1 = x_pos + 0.04,
          y1 = label_y[moved],
          col = grDevices::adjustcolor(pal["facet_level"], alpha.f = 0.45),
          lwd = 0.8
        )
      }
      graphics::text(
        x = x_pos,
        y = label_y,
        labels = as.character(sub$DisplayLabel),
        cex = 0.64,
        col = pal["facet_level"]
      )
      ci_ok <- is.finite(sub$OriginalCI_Lower) &
        is.finite(sub$OriginalCI_Upper)
      if (any(ci_ok)) {
        ci_x_all <- if (nrow(sub) <= 1L) {
          rep(x_pos - 0.36, nrow(sub))
        } else {
          seq(x_pos - 0.42, x_pos - 0.30, length.out = nrow(sub))
        }
        regular <- ci_ok & !sub$CISuppressed
        if (any(regular)) {
          ci_x <- ci_x_all[regular]
          lo <- sub$DisplayCI_Lower[regular]
          hi <- sub$DisplayCI_Upper[regular]
          graphics::segments(
            x0 = ci_x, y0 = lo,
            x1 = ci_x, y1 = hi,
            col = grDevices::adjustcolor(pal["facet_level"], alpha.f = 0.65), lwd = 1
          )
          graphics::segments(
            x0 = ci_x - 0.035, y0 = lo,
            x1 = ci_x + 0.035, y1 = lo,
            col = grDevices::adjustcolor(pal["facet_level"], alpha.f = 0.65)
          )
          graphics::segments(
            x0 = ci_x - 0.035, y0 = hi,
            x1 = ci_x + 0.035, y1 = hi,
            col = grDevices::adjustcolor(pal["facet_level"], alpha.f = 0.65)
          )
          if (any(sub$CIClippedLower[regular])) {
            lower_marker <- regular & sub$CIClippedLower
            graphics::points(
              x = ci_x_all[lower_marker], y = yr[1],
              pch = 25, bg = "white", col = pal["facet_level"], cex = 0.62
            )
          }
          if (any(sub$CIClippedUpper[regular])) {
            upper_marker <- regular & sub$CIClippedUpper
            graphics::points(
              x = ci_x_all[upper_marker], y = yr[2],
              pch = 24, bg = "white", col = pal["facet_level"], cex = 0.62
            )
          }
          graphics::points(
            x = ci_x,
            y = sub$DisplayEstimate[regular],
            pch = 16,
            cex = 0.42,
            col = pal["facet_level"]
          )
        }
        boundary_interval <- ci_ok & sub$CISuppressed
        if (any(boundary_interval)) {
          boundary_pch <- ifelse(
            sub$BoundaryEnd[boundary_interval] == "lower",
            25,
            ifelse(sub$BoundaryEnd[boundary_interval] == "upper", 24, 23)
          )
          graphics::points(
            x = ci_x_all[boundary_interval],
            y = sub$DisplayEstimate[boundary_interval],
            pch = boundary_pch,
            bg = "white",
            col = pal["facet_level"],
            cex = 0.78,
            lwd = 1.1
          )
        }
      }
    }
  }

  step_headers <- headers[headers$ColumnType == "step", , drop = FALSE]
  step_ruler <- facets_data$step_ruler
  half_ruler <- facets_data$score_transitions
  for (i in seq_len(nrow(step_headers))) {
    x_pos <- match(step_headers$Header[i], headers$Header)
    step_sub <- step_ruler[
      step_ruler$CurveGroup == step_headers$Group[i] & step_ruler$InRange,
      , drop = FALSE
    ]
    if (nrow(step_sub) > 0L) {
      step_y <- suppressWarnings(as.numeric(step_sub$DrawValue))
      graphics::segments(
        x0 = x_pos - 0.34, y0 = step_y,
        x1 = x_pos + 0.34, y1 = step_y,
        col = pal["step_threshold"], lwd = 1.5
      )
      label_pos <- rep(3L, nrow(step_sub))
      step_order <- order(step_y)
      if (length(step_order) > 1L) {
        close_pair <- diff(step_y[step_order]) < 0.09 * diff(yr)
        cluster <- cumsum(c(TRUE, !close_pair))
        for (cluster_id in unique(cluster)) {
          members <- step_order[cluster == cluster_id]
          if (length(members) > 1L) {
            label_pos[members] <- rep(c(1L, 3L), length.out = length(members))
          }
        }
      }
      draw_labels <- vapply(
        as.character(step_sub$DrawLabel),
        function(value) paste(strwrap(value, width = 20L), collapse = "\n"),
        character(1)
      )
      graphics::text(
        x = x_pos,
        y = step_y,
        labels = draw_labels,
        pos = label_pos, offset = 0.16, cex = 0.54,
        col = pal["step_threshold"]
      )
    }
    half_sub <- half_ruler[
      half_ruler$CurveGroup == step_headers$Group[i] & half_ruler$InRange,
      , drop = FALSE
    ]
    if (nrow(half_sub) > 0L) {
      graphics::segments(
        x0 = x_pos - 0.22, y0 = half_sub$DrawValue,
        x1 = x_pos + 0.22, y1 = half_sub$DrawValue,
        col = grDevices::adjustcolor(pal["step_threshold"], alpha.f = 0.55),
        lty = 3, lwd = 1
      )
    }
  }

  graphics::mtext(
    sprintf("* = %s person(s); rows/logit = %d", format(settings$PersonsPerStar[1], trim = TRUE), settings$RowsPerLogit[1]),
    side = 1, line = 2.7, adj = 0, cex = 0.75
  )
  graphics::mtext(
    "FACETS Table 6-style visual layout; estimates remain mfrmr estimates (not numerical equivalence).",
    side = 1, line = 3.7, adj = 0, cex = 0.67, col = "gray35"
  )
  if (isTRUE(show_ci)) {
    graphics::mtext(
      sprintf(
        "Whiskers show %g%% mfrmr confidence intervals; triangles mark bounds beyond the displayed ruler.",
        round(100 * ci_level)
      ),
      side = 1, line = 4.6, adj = 0, cex = 0.64,
      col = pal["facet_level"]
    )
  }
  if (identical(auto_range_policy, "boundary_levels_at_ends")) {
    graphics::mtext(
      if (isTRUE(show_ci) && boundary_ci_count > 0L) {
        paste0(
          "Boundary-separated levels use end triangles and their intervals are omitted from the ruler; ",
          "inspect OriginalEstimate and CI_Lower/CI_Upper."
        )
      } else {
        "Boundary-separated levels are placed at ruler ends; inspect OriginalEstimate for their untruncated values."
      },
      side = 1, line = if (isTRUE(show_ci)) 5.4 else 4.6,
      adj = 0, cex = 0.58, col = "#9A3412"
    )
  } else if (isTRUE(show_ci) && ci_clipped_count > 0L) {
    graphics::mtext(
      "CI endpoint triangles identify clipping; exact bounds remain in CI_Lower/CI_Upper.",
      side = 1, line = 5.4, adj = 0, cex = 0.6, col = "#9A3412"
    )
  }
  invisible(NULL)
}
