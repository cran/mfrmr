# ==============================================================================
# Data preparation, indexing, and small string / formatting utilities
# ==============================================================================
#
# Core helpers for `fit_mfrm()`-style pipelines. Split out of
# `mfrm_core.R` so the data-validation, indexing, and
# print / formatting utilities live in a single browseable file.
# Functions are internal (no @export); they are called directly by
# `mfrm_estimate()` and the surrounding orchestration helpers.

# ---- data preparation ----
mfrm_single_person_caution <- function(prep) {
  if (!identical(as.integer(prep$n_person), 1L)) return(character())
  paste("Only one Person remains. Calibration depends on that response vector and the assumed within-Person conditional independence.",
    "Any returned model-information interval is an unqualified normal approximation; its finite-sample coverage is not established.",
    "One Person cannot supply independent-Person sandwich uncertainty or empirical evidence of between-Person variation.")
}

mfrm_preparation_note <- function(stage,
                                  condition,
                                  severity = "info",
                                  count = NA_integer_,
                                  affected = "",
                                  message = "",
                                  action = "") {
  data.frame(
    Stage = as.character(stage %||% ""),
    Condition = as.character(condition %||% ""),
    Severity = as.character(severity %||% "info"),
    Count = suppressWarnings(as.integer(count %||% NA_integer_)),
    Affected = paste(as.character(affected %||% ""), collapse = ", "),
    Message = as.character(message %||% ""),
    RecommendedAction = as.character(action %||% ""),
    stringsAsFactors = FALSE
  )
}

bind_mfrm_preparation_notes <- function(notes) {
  notes <- notes[!vapply(notes, is.null, logical(1))]
  if (length(notes) == 0L) {
    return(data.frame(
      Stage = character(0),
      Condition = character(0),
      Severity = character(0),
      Count = integer(0),
      Affected = character(0),
      Message = character(0),
      RecommendedAction = character(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, notes)
}

mfrm_validate_facet_names <- function(facets) {
  reserved <- intersect(facets, c("Person", "Score", "Weight", "score_k",
    "theta", "steps", "interactions", "log_slopes", "beta", "log_sigma2"))
  if (length(reserved)) {
    stop("Facet columns cannot use reserved internal names: ",
      paste(reserved, collapse = ", "),
      ". Rename these facet columns before analysis; affected saved fits require refitting from the original data.",
      call. = FALSE)
  }
  invisible(NULL)
}

prepare_mfrm_data <- function(data, person_col, facet_cols, score_col,
                              rating_min = NULL, rating_max = NULL,
                              weight_col = NULL, keep_original = FALSE,
                              missing_codes = NULL,
                              repeated_ratings = c("review", "independent"),
                              allow_single_person = FALSE,
                              allow_constant_scores = FALSE) {
  repeated_ratings <- match.arg(repeated_ratings)
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  # Model columns, not inherited dplyr grouping, define the analysis rows.
  data <- dplyr::ungroup(data)
  preparation_notes <- list()
  emit_preparation_messages <- isTRUE(getOption("mfrmr.show_preparation_messages", FALSE)) &&
    !isTRUE(getOption("mfrmr._preparation_messages_announced", FALSE))
  add_preparation_note <- function(stage,
                                   condition,
                                   severity = "info",
                                   count = NA_integer_,
                                   affected = "",
                                   message = "",
                                   action = "") {
    preparation_notes[[length(preparation_notes) + 1L]] <<-
      mfrm_preparation_note(
        stage = stage,
        condition = condition,
        severity = severity,
        count = count,
        affected = affected,
        message = message,
        action = action
      )
    invisible(NULL)
  }

  required <- c(person_col, facet_cols, score_col)
  if (!is.null(weight_col)) {
    required <- c(required, weight_col)
  }
  if (length(unique(required)) != length(required)) {
    dup_names <- required[duplicated(required)]
    stop("The 'person', 'score', and 'facets' arguments must name distinct columns, ",
         "but duplicates were found: ", paste(dup_names, collapse = ", "), ". ",
         "Remove or rename the duplicated references.", call. = FALSE)
  }
  mfrm_validate_facet_names(facet_cols)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "), ". ",
         "fit_mfrm() expects long-format data with one person column, one score column, ",
         "and one or more facet columns. Available columns: ",
         paste(names(data), collapse = ", "), ". ",
         "Check spelling of person/facets/score arguments or reshape the data to long format.",
         call. = FALSE)
  }

  # Optional pre-processing step that converts
  # FACETS / SPSS / SAS sentinel values to NA so the drop_na() downstream
  # behaves intuitively. Accepts TRUE / "default" for the conventional
  # set, or a character vector of custom codes. Recoding is restricted
  # to person / facets / score (weight retains its original meaning).
  missing_audit <- NULL
  if (!is.null(missing_codes) && !isFALSE(missing_codes)) {
    default_missing_codes <- isTRUE(missing_codes) ||
      (is.character(missing_codes) &&
       length(missing_codes) == 1L &&
       identical(tolower(missing_codes), "default"))
    codes <- if (default_missing_codes) {
      c("99", "999", "-1", "N", "NA", "n/a", ".", "")
    } else {
      as.character(missing_codes)
    }
    recode_cols <- c(person_col, facet_cols, score_col)
    # The conventional set contains short tokens such as "N" and numeric
    # strings such as "99" that can be legitimate person or facet labels.
    # Applying it indiscriminately to identifier columns can silently delete a
    # complete person or facet level. Therefore the convenience policy acts on
    # the score column only. An explicit character vector remains an explicit
    # request to apply those codes across all selected model columns.
    applied_cols <- if (default_missing_codes) score_col else recode_cols
    data <- recode_missing_codes(
      data,
      columns = applied_cols,
      codes = codes,
      verbose = FALSE
    )
    applied_audit <- attr(data, "mfrm_missing_recoding")
    missing_audit <- data.frame(
      Column = recode_cols,
      Replaced = 0L,
      Scope = ifelse(
        recode_cols %in% applied_cols,
        if (default_missing_codes) "default_score_only" else "explicit_all_model_columns",
        "identifier_preserved_by_default"
      ),
      stringsAsFactors = FALSE
    )
    if (is.data.frame(applied_audit) && nrow(applied_audit) > 0L) {
      matched <- match(applied_audit$Column, missing_audit$Column)
      keep <- !is.na(matched)
      missing_audit$Replaced[matched[keep]] <- as.integer(applied_audit$Replaced[keep])
    }
    attr(data, "mfrm_missing_recoding") <- missing_audit
    replaced_n <- if (is.data.frame(missing_audit) &&
                      "Replaced" %in% names(missing_audit)) {
      sum(suppressWarnings(as.integer(missing_audit$Replaced)), na.rm = TRUE)
    } else {
      0L
    }
    if (replaced_n > 0L) {
      replaced_cols <- missing_audit$Column[
        suppressWarnings(as.integer(missing_audit$Replaced)) > 0L
      ]
      add_preparation_note(
        stage = "missing_code_recode",
        condition = "declared_missing_codes_replaced",
        severity = "review",
        count = replaced_n,
        affected = paste(replaced_cols, collapse = ", "),
        message = paste0(
          "Replaced ", replaced_n,
          " declared missing-code value(s) with NA before row filtering."
        ),
        action = "Inspect `fit$prep$missing_recoding` (or `describe_mfrm_data(...)$missing_recoding`) before reporting exclusions."
      )
    }
  }
  if (any(duplicated(names(data)))) {
    dupes <- unique(names(data)[duplicated(names(data))])
    if (any(required %in% dupes)) {
      stop("Selected columns include duplicate names in the data: ",
           paste(intersect(required, dupes), collapse = ", "), ". ",
           "Rename columns so each name is unique.", call. = FALSE)
    }
  }
  if (length(facet_cols) == 0) {
    stop("No facet columns were specified. ",
         "Supply at least one column name via 'facets' from the long-format rating table ",
         "(e.g., facets = c('Rater', 'Task')).", call. = FALSE)
  }

  cols <- c(person_col, facet_cols, score_col)
  if (!is.null(weight_col)) {
    cols <- c(cols, weight_col)
  }
  df <- data |>
    select(all_of(cols)) |>
    rename(
      Person = all_of(person_col),
      Score = all_of(score_col)
    )
  if (!is.null(weight_col)) {
    df <- df |> rename(Weight = all_of(weight_col))
  }

  raw_score <- as.character(df$Score)
  raw_weight <- if ("Weight" %in% names(df)) as.character(df$Weight) else NULL
  input_rows <- nrow(df)

  if (is.complex(df$Score) || !is.null(dim(df$Score))) {
    stop("`Score` must be a vector of real ordered integer category codes.", call. = FALSE)
  }
  score_num <- if (is.numeric(df$Score)) as.numeric(df$Score) else
    suppressWarnings(as.numeric(raw_score))
  # Practical integer tolerance. sqrt(.Machine$double.eps) (~1.5e-8) was
  # too strict: CSV round-trip artifacts like "1.0000001" (diff 1e-7
  # from an integer) were rejected as fractional. 1e-6 still catches
  # real fractional scores like 1.5 or 2.75 while accepting float
  # representation noise for integer codes.
  score_tol <- 1e-6
  nonfinite_score <- is.nan(score_num) |
    (!is.na(score_num) & !is.finite(score_num))
  bad_score <- is.na(score_num) & !is.nan(score_num) &
    !is.na(raw_score) & nzchar(trimws(raw_score))

  # If essentially every non-empty value is non-numeric (e.g. "low", "medium",
  # "high" text labels), the later `drop_na` would silently remove every row
  # and surface as an unhelpful "No valid observations" error. Fail loudly up
  # front with a targeted message instead.
  non_empty <- !is.na(raw_score) & nzchar(trimws(raw_score))
  n_non_empty <- sum(non_empty)
  if (n_non_empty > 0L && sum(bad_score) == n_non_empty) {
    examples <- utils::head(unique(raw_score[bad_score]), 5L)
    stop(
      "`Score` column appears to contain text labels (e.g. ",
      paste(shQuote(examples), collapse = ", "),
      ") rather than ordered integer category codes. Recode to integers ",
      "(for example 0/1 or 1:5) before calling fit_mfrm().",
      call. = FALSE
    )
  }

  if (any(nonfinite_score)) {
    examples <- utils::head(unique(raw_score[nonfinite_score]), 5L)
    stop(
      "`Score` must contain finite ordered integer category codes. ",
      "Non-finite value(s) were found: ",
      paste(examples, collapse = ", "), ". ",
      "Recode or remove these rows before fitting.",
      call. = FALSE
    )
  }

  if (any(bad_score)) {
    bad_score_n <- sum(bad_score)
    bad_score_msg <- paste0(
      "`Score` contained ", bad_score_n,
      " non-numeric value(s); affected row(s) will be removed before estimation."
    )
    add_preparation_note(
      stage = "score_coercion",
      condition = "non_numeric_score",
      severity = "warning",
      count = bad_score_n,
      affected = score_col,
      message = bad_score_msg,
      action = "Recode score values to ordered integer categories before fitting."
    )
    warning(
      bad_score_msg,
      call. = FALSE
    )
  }
  fractional_score <- is.finite(score_num) &
    (abs(score_num - round(score_num)) > score_tol)
  if (any(fractional_score)) {
    fractional_examples <- unique(raw_score[fractional_score])
    fractional_examples <- utils::head(fractional_examples, n = 5L)
    stop(
      "`Score` must contain ordered integer category codes (for example 0/1, 1/2, or 1:5). ",
      "Fractional value(s) were found: ", paste(fractional_examples, collapse = ", "), ". ",
      "Recode the score column explicitly before fitting.",
      call. = FALSE
    )
  }
  if (any(is.finite(score_num) & abs(round(score_num)) > .Machine$integer.max)) {
    stop("`Score` category codes exceed R's integer range. Recode categories explicitly before fitting.", call. = FALSE)
  }
  raw_person_id <- as.character(df$Person)
  raw_facet_id <- lapply(facet_cols, function(f) as.character(df[[f]]))
  names(raw_facet_id) <- facet_cols
  df <- df |>
    mutate(
      Person = enc2utf8(trimws(as.character(Person))),
      across(all_of(facet_cols), ~ enc2utf8(trimws(as.character(.x)))),
      Score = score_num
    )
  blank_person <- !is.na(df$Person) & !nzchar(df$Person)
  blank_facets <- vapply(facet_cols, function(facet) {
    sum(!is.na(df[[facet]]) & !nzchar(df[[facet]]))
  }, integer(1))
  if (any(blank_person) || any(blank_facets > 0L)) {
    blank_facet_names <- names(blank_facets)[blank_facets > 0L]
    blank_facet_labels <- if (length(blank_facet_names) > 0L) {
      paste0(blank_facet_names, " (", blank_facets[blank_facet_names], ")")
    } else {
      character(0)
    }
    affected <- c(
      if (any(blank_person)) paste0(person_col, " (", sum(blank_person), ")"),
      blank_facet_labels
    )
    stop(
      "Blank person or facet identifier(s) were found after trimming whitespace: ",
      paste(affected, collapse = ", "), ". ",
      "Fill, remove, or explicitly recode these identifiers before fitting; ",
      "blank labels cannot define measurement levels.",
      call. = FALSE
    )
  }
  # Detect Person / facet IDs that gained / lost surrounding whitespace
  # in the trim and warn so users do not silently end up with a "P01"
  # vs " P01 " split.
  trimmed_person_diff_n <- sum(raw_person_id != df$Person, na.rm = TRUE)
  if (trimmed_person_diff_n > 0L) {
    trim_msg <- paste0(
      "Trimmed leading/trailing whitespace from `", person_col, "` ",
      "in ", trimmed_person_diff_n, " row(s). Affected IDs were treated ",
      "as the trimmed value; pre-clean them upstream if you need the ",
      "original spelling."
    )
    add_preparation_note(
      stage = "id_normalization",
      condition = "trimmed_person_ids",
      severity = "info",
      count = trimmed_person_diff_n,
      affected = person_col,
      message = trim_msg,
      action = "Pre-clean ID whitespace upstream if the original spelling must be preserved."
    )
    if (isTRUE(emit_preparation_messages)) {
      message(trim_msg)
    }
  }
  for (facet in facet_cols) {
    trimmed_facet_diff_n <- sum(raw_facet_id[[facet]] != df[[facet]], na.rm = TRUE)
    if (trimmed_facet_diff_n > 0L) {
      add_preparation_note(
        stage = "id_normalization",
        condition = "trimmed_facet_ids",
        severity = "info",
        count = trimmed_facet_diff_n,
        affected = facet,
        message = paste0(
          "Trimmed leading/trailing whitespace from `", facet, "` in ",
          trimmed_facet_diff_n, " row(s)."
        ),
        action = "Pre-clean facet labels upstream if whitespace is meaningful."
      )
    }
  }
  if (!"Weight" %in% names(df)) {
    df <- df |> mutate(Weight = 1)
  } else {
    if (is.complex(df$Weight) || !is.null(dim(df$Weight))) {
      stop("`Weight` must be a vector of real numeric values.", call. = FALSE)
    }
    weight_num <- if (is.numeric(df$Weight)) as.numeric(df$Weight) else
      suppressWarnings(as.numeric(raw_weight))
    nonfinite_weight <- is.nan(weight_num) |
      (!is.na(weight_num) & !is.finite(weight_num))
    if (any(nonfinite_weight)) {
      examples <- utils::head(unique(raw_weight[nonfinite_weight]), 5L)
      stop(
        "`Weight` must contain finite numeric values. Non-finite value(s) ",
        "were found: ", paste(examples, collapse = ", "), ". ",
        "Use finite weights or remove the weight column; non-positive finite ",
        "weights are excluded during data preparation.",
        call. = FALSE
      )
    }
    bad_weight <- is.na(weight_num) & !is.nan(weight_num) &
      !is.na(raw_weight) & nzchar(trimws(raw_weight))
    if (any(bad_weight)) {
      bad_weight_n <- sum(bad_weight)
      bad_weight_msg <- paste0(
        "`Weight` contained ", bad_weight_n,
        " non-numeric value(s); affected row(s) will be removed before estimation."
      )
      add_preparation_note(
        stage = "weight_coercion",
        condition = "non_numeric_weight",
        severity = "warning",
        count = bad_weight_n,
        affected = weight_col,
        message = bad_weight_msg,
        action = "Recode weights to positive numeric values or remove the weight column."
      )
      warning(
        bad_weight_msg,
        call. = FALSE
      )
    }
    df <- df |> mutate(Weight = weight_num)
  }

  rows_before_drop <- nrow(df)
  omitted <- !stats::complete.cases(df) | (!is.na(df$Weight) & df$Weight <= 0)
  omitted_data <- as.data.frame(df[omitted, , drop = FALSE])
  omitted_input_rows <- which(omitted)
  df <- df |>
    tidyr::drop_na() |>
    filter(Weight > 0)
  rows_dropped <- rows_before_drop - nrow(df)
  row_retention <- data.frame(
    Stage = c("input_selected_columns", "after_missing_and_weight_filter"),
    Rows = c(input_rows, nrow(df)),
    DroppedRows = c(0L, rows_dropped),
    DroppedReason = c("", "missing values or non-positive weights"),
    stringsAsFactors = FALSE
  )
  if (rows_dropped > 0L) {
    drop_msg <- paste0(
      "Dropped ", rows_dropped, " row(s) with missing values or non-positive ",
      "weights before estimation. Pass `missing_codes = ...` to recode ",
      "user-specified missing markers, or pre-process upstream if you need ",
      "to keep the row."
    )
    add_preparation_note(
      stage = "row_filter",
      condition = "missing_or_nonpositive_weight_rows_dropped",
      severity = "review",
      count = rows_dropped,
      affected = paste(c(required), collapse = ", "),
      message = drop_msg,
      action = "Inspect `fit$prep$row_retention` and use `missing_codes = ...` or upstream cleaning if dropped rows were unintended."
    )
    if (isTRUE(emit_preparation_messages)) {
      message(drop_msg)
    }
  }

  if (nrow(df) == 0) {
    stop("No valid observations remain after removing missing values and ",
         "zero-weight rows. Check that person, facet, score, and weight columns ",
         "contain valid (non-NA, non-empty) data.", call. = FALSE)
  }

  # Repeated cell keys do not establish dependence. Keep the default review
  # unless the analyst explicitly declares conditionally independent events.
  if (nrow(df) > 0L) {
    key_cols <- c("Person", facet_cols)
    dup_mask <- duplicated(df[, key_cols, drop = FALSE]) |
                duplicated(df[, key_cols, drop = FALSE], fromLast = TRUE)
    n_dup <- sum(dup_mask)
    if (n_dup > 0L && identical(repeated_ratings, "independent")) {
      add_preparation_note(
        stage = "design_check",
        condition = "independent_repeated_ratings_declared",
        severity = "info", count = n_dup,
        affected = paste(key_cols, collapse = ", "),
        message = "Repeated Person-by-facet ratings are retained as separate conditionally independent events, as explicitly declared by the analyst.",
        action = "This declaration does not verify independence or distinguish accidental copies. Use a dependence model when ratings remain associated after conditioning on the model."
      )
    } else if (n_dup > 0L) {
      n_dup_cells <- nrow(unique(df[dup_mask, key_cols, drop = FALSE]))
      duplicate_msg <- paste0(
        "Detected ", n_dup, " row(s) across ", n_dup_cells,
        " duplicated Person x ",
        "(", paste(facet_cols, collapse = ", "), ") combination. MFRM ",
        "requires review of repeated ratings; correct accidental copies, model ",
        "dependence, or explicitly declare repeated_ratings = 'independent' ",
        "when conditional independence is justified. Continuing ",
        "with the rows as supplied."
      )
      add_preparation_note(
        stage = "design_check",
        condition = "duplicate_person_facet_cells",
        severity = "warning",
        count = n_dup,
        affected = paste(c("Person", facet_cols), collapse = ", "),
        message = duplicate_msg,
        action = "Review repeated events: correct accidental copies, model dependence, or explicitly declare conditional independence."
      )
      duplicate_warning_announced <- isTRUE(
        getOption("mfrmr._duplicate_warning_announced", FALSE)
      )
      if (!duplicate_warning_announced) {
        warning(duplicate_msg, call. = FALSE)
        if (!is.null(getOption("mfrmr._duplicate_warning_announced"))) {
          options(mfrmr._duplicate_warning_announced = TRUE)
        }
      }
    }
  }

  rounded_n <- sum(df$Score != round(df$Score))
  if (rounded_n > 0L) {
    add_preparation_note(
      stage = "score_coercion", condition = "score_roundoff_normalized",
      severity = "info", count = rounded_n, affected = score_col,
      message = "Normalized retained scores within 1e-6 of an integer to the nearest category code.",
      action = "Inspect the source scores; genuine fractional responses require explicit recoding."
    )
  }
  df <- df |>
    mutate(Score = as.integer(round(Score)))
  if (!is.finite(sum(df$Weight))) {
    stop("The sum of retained weights exceeds finite numeric range. Choose a representable weighting scheme explicitly; weights are not rescaled automatically.", call. = FALSE)
  }

  observed_score_values <- sort(unique(df$Score))

  if (length(unique(df$Score)) < 2 && !isTRUE(allow_constant_scores)) {
    stop("Only one score category found in the data (Score = ",
         unique(df$Score), "). ",
         "This calibration interface currently requires at least two distinct response categories.", call. = FALSE)
  }

  rating_min_supplied <- !is.null(rating_min)
  rating_max_supplied <- !is.null(rating_max)
  if (is.null(rating_min)) rating_min <- min(df$Score, na.rm = TRUE)
  if (is.null(rating_max)) rating_max <- max(df$Score, na.rm = TRUE)
  rating_min_source <- if (rating_min_supplied) "declared" else "observed"
  rating_max_source <- if (rating_max_supplied) "declared" else "observed"
  rating_range_source <- if (rating_min_supplied && rating_max_supplied) {
    "declared"
  } else if (!rating_min_supplied && !rating_max_supplied) {
    "observed"
  } else {
    "partly_declared"
  }
  inferred_rating_message <- paste0(
    "Rating range inferred from observed scores: ",
    paste(
      c(
        if (!rating_min_supplied) paste0("rating_min = ", rating_min),
        if (!rating_max_supplied) paste0("rating_max = ", rating_max)
      ),
      collapse = ", "
    ),
    ". Supply `rating_min`/`rating_max` explicitly if the declared scale ",
    "differs from the observed range."
  )
  # Keep routine fits quiet by default: large simulation/recovery workflows can
  # call this helper hundreds of times. The declared-vs-observed provenance is
  # stored in the returned prep object and surfaced by summary/data-description
  # outputs. Users who want the previous message stream can opt in with
  # options(mfrmr.show_inferred_rating_range = TRUE). fit_mfrm() sets the
  # session-scoped private option `mfrmr._rating_range_announced` so that an
  # opt-in message still appears at most once per top-level fit.
  if (!rating_min_supplied || !rating_max_supplied) {
    already_announced <- isTRUE(getOption("mfrmr._rating_range_announced"))
    show_inferred_message <- isTRUE(getOption("mfrmr.show_inferred_rating_range", FALSE))
    if (show_inferred_message && !already_announced) {
      message(inferred_rating_message)
      # Flip the flag if it exists (fit_mfrm set it to FALSE up front). When
      # called outside fit_mfrm, the option is NULL and we leave it as-is so
      # subsequent standalone calls still announce when the opt-in is active.
      if (!is.null(getOption("mfrmr._rating_range_announced"))) {
        options(mfrmr._rating_range_announced = TRUE)
      }
    }
  }
  if (!is.numeric(rating_min) || is.complex(rating_min) || length(rating_min) != 1L ||
      !is.finite(rating_min) || abs(rating_min - round(rating_min)) > score_tol ||
      abs(round(rating_min)) > .Machine$integer.max) {
    stop("`rating_min` must be a single finite integer category value.", call. = FALSE)
  }
  if (!is.numeric(rating_max) || is.complex(rating_max) || length(rating_max) != 1L ||
      !is.finite(rating_max) || abs(rating_max - round(rating_max)) > score_tol ||
      abs(round(rating_max)) > .Machine$integer.max) {
    stop("`rating_max` must be a single finite integer category value.", call. = FALSE)
  }
  rating_min <- as.integer(round(rating_min))
  rating_max <- as.integer(round(rating_max))
  if (rating_max <= rating_min) {
    stop("`rating_max` must be larger than `rating_min`.", call. = FALSE)
  }
  if (as.double(rating_max) - rating_min + 1 > .Machine$integer.max) {
    stop("The declared category range exceeds R's integer indexing capacity. Recode categories explicitly.", call. = FALSE)
  }
  explicit_rating_range <- rating_min_supplied || rating_max_supplied
  expected_vals <- seq(rating_min, rating_max)
  out_of_range <- observed_score_values[observed_score_values < rating_min | observed_score_values > rating_max]
  if (length(out_of_range) > 0L) {
    stop(
      "Observed `Score` categories fall outside the supplied rating range: ",
      paste(out_of_range, collapse = ", "),
      ". Adjust `rating_min`/`rating_max` or recode the score column before fitting.",
      call. = FALSE
    )
  }

  preserve_score_support <- isTRUE(keep_original)
  if (!isTRUE(keep_original)) {
    score_vals <- sort(unique(df$Score))
    observed_contiguous <- identical(score_vals, seq(min(score_vals), max(score_vals)))
    boundary_only_gap <- isTRUE(explicit_rating_range) &&
      observed_contiguous &&
      all(score_vals %in% expected_vals)
    if (!identical(score_vals, expected_vals) && !isTRUE(boundary_only_gap)) {
      recoded_vals <- seq(rating_min, rating_min + length(score_vals) - 1L)
      suspected_missing_codes <- intersect(score_vals, c(-1L, 99L, 999L))
      sentinel_note <- if (length(suspected_missing_codes) > 0L &&
                           (is.null(missing_codes) || isFALSE(missing_codes))) {
        paste0(
          " Value(s) ", paste(suspected_missing_codes, collapse = ", "),
          " are commonly used as missing-value codes in FACETS, SPSS, or SAS exports. ",
          "If they mean missing data here, refit with `missing_codes = TRUE` or an explicit code vector."
        )
      } else {
        ""
      }
      recode_msg <- paste0(
        "Observed `Score` categories were non-consecutive (",
        paste(score_vals, collapse = ", "),
        ") and were recoded internally to a contiguous scale (",
        paste(recoded_vals, collapse = ", "),
        ") because `keep_original = FALSE`. Inspect `fit$prep$score_map` ",
        "before interpreting category labels, or set `keep_original = TRUE` ",
        "to preserve the declared score support.", sentinel_note
      )
      add_preparation_note(
        stage = "score_support",
        condition = "score_categories_recoded",
        severity = "warning",
        count = length(score_vals),
        affected = score_col,
        message = recode_msg,
        action = if (length(suspected_missing_codes) > 0L &&
                     (is.null(missing_codes) || isFALSE(missing_codes))) {
          "Confirm whether common sentinel values represent missing data; then set `missing_codes`, `keep_original`, and the rating range explicitly."
        } else {
          "Inspect `fit$prep$score_map` and document the mapping, or preserve the intended support explicitly."
        }
      )
      score_recode_announced <- isTRUE(getOption("mfrmr._score_recode_announced"))
      if (!score_recode_announced) {
        warning(recode_msg, call. = FALSE)
        if (!is.null(getOption("mfrmr._score_recode_announced"))) {
          options(mfrmr._score_recode_announced = TRUE)
        }
      }
      df <- df |>
        mutate(Score = match(Score, score_vals) + rating_min - 1)
      rating_max <- rating_min + length(score_vals) - 1
      expected_vals <- seq(rating_min, rating_max)
    } else if (isTRUE(boundary_only_gap)) {
      preserve_score_support <- TRUE
    }
  }

  if (isTRUE(preserve_score_support)) {
    score_map <- tibble(
      OriginalScore = seq(rating_min, rating_max),
      InternalScore = seq(rating_min, rating_max)
    )
  } else {
    score_map <- tibble(
      OriginalScore = observed_score_values,
      InternalScore = seq(rating_min, rating_max)
    )
  }

  df <- df |>
    mutate(score_k = Score - rating_min)
  unused_score_categories <- setdiff(seq(rating_min, rating_max), sort(unique(df$Score)))

  # Guard against silently dropped rows. `score_k < 0` or `score_k >= n_cat`
  # produces a 0 or out-of-range matrix index later; R's `m[cbind(i, 0)]`
  # returns zero-length silently, so such rows would not contribute to
  # the likelihood while the reported n_obs stays unchanged.
  n_cat <- rating_max - rating_min + 1L
  bad_score_k <- !is.finite(df$score_k) | df$score_k < 0L | df$score_k >= n_cat
  if (any(bad_score_k)) {
    bad_scores <- sort(unique(df$Score[bad_score_k]))
    stop(
      "`Score` contains value(s) outside the declared category range [",
      rating_min, ", ", rating_max, "]: ",
      paste(utils::head(bad_scores, 10), collapse = ", "),
      if (length(bad_scores) > 10L) ", ..." else "",
      ". Check `rating_min`/`rating_max` or recode the score column before fitting.",
      call. = FALSE
    )
  }

  df <- df |>
    mutate(
      Person = factor(Person),
      across(all_of(facet_cols), ~ factor(.x))
    )

  facet_names <- facet_cols
  facet_levels <- lapply(facet_names, function(f) levels(df[[f]]))
  names(facet_levels) <- facet_names

  # Single-Person support is admitted by the caller; it is not a universal
  # identification theorem. Unreviewed fitting routes retain their input limit.
  n_person <- length(levels(df$Person))
  if (n_person < 2L && !isTRUE(allow_single_person)) {
    stop(
      "This fitting route requires at least two persons (got ", n_person, "). ",
      "One-Person calibration is supported only for full RSM/PCM JML, CML or fixed-population MML. ",
      "The selected route retains a software support limit, not a general mathematical impossibility claim.",
      call. = FALSE
    )
  }
  # A count alone is not an invalid-input condition. Interval outputs also
  # carry this caution, separately from their numerical eligibility checks.
  if (n_person == 1L) add_preparation_note(
    stage = "design_check", condition = "single_person_calibration",
    severity = "info", count = 1L, affected = levels(df$Person),
    message = mfrm_single_person_caution(list(n_person = n_person)),
    action = "Review the independent rating events, identification and boundary checks; do not interpret numerical availability as population-level validation.")
  # Record facets that have only a single observed level. They are structurally
  # identified at 0 by the sum-to-zero constraint and do not contribute
  # measurement information; users often supply them by mistake.
  single_level <- facet_names[vapply(facet_levels, length, integer(1)) <= 1L]
  if (length(single_level) > 0L) {
    single_level_msg <- paste0(
      "Facet(s) with only a single observed level: ",
      paste(shQuote(single_level), collapse = ", "),
      ". They will be fixed at 0 by the sum-to-zero constraint and cannot ",
      "inform the fit."
    )
    add_preparation_note(
      stage = "design_check",
      condition = "single_level_facet",
      severity = "review",
      count = length(single_level),
      affected = paste(single_level, collapse = ", "),
      message = single_level_msg,
      action = "Remove the facet from the model or collect at least two observed levels if the facet should be estimated."
    )
    if (isTRUE(emit_preparation_messages)) {
      message(single_level_msg)
    }
  }
  preparation_notes <- bind_mfrm_preparation_notes(preparation_notes)
  if (isTRUE(emit_preparation_messages) &&
      !is.null(getOption("mfrmr._preparation_messages_announced"))) {
    options(mfrmr._preparation_messages_announced = TRUE)
  }

  list(
    data = df,
    n_obs = nrow(df),
    repeated_ratings = repeated_ratings,
    weighted_n = sum(df$Weight, na.rm = TRUE),
    n_person = length(levels(df$Person)),
    rating_min = rating_min,
    rating_max = rating_max,
    rating_range_source = rating_range_source,
    rating_min_source = rating_min_source,
    rating_max_source = rating_max_source,
    rating_range_message = if (!rating_min_supplied || !rating_max_supplied) {
      inferred_rating_message
    } else {
      ""
    },
    score_map = score_map,
    unused_score_categories = unused_score_categories,
    facet_names = facet_names,
    levels = c(list(Person = levels(df$Person)), facet_levels),
    weight_col = if (!is.null(weight_col)) weight_col else NULL,
    keep_original = isTRUE(keep_original),
    row_retention = row_retention,
    omitted_data = omitted_data,
    omitted_input_rows = omitted_input_rows,
    preparation_notes = preparation_notes,
    source_columns = list(
      person = person_col,
      facets = facet_cols,
      score = score_col,
      weight = if (!is.null(weight_col)) weight_col else NULL
    ),
    missing_recoding = missing_audit
  )
}

build_indices <- function(prep, step_facet = NULL, slope_facet = NULL,
                          interaction_specs = NULL, gpcm_spec = NULL) {
  df <- prep$data
  facets_idx <- lapply(prep$facet_names, function(f) as.integer(df[[f]]))
  names(facets_idx) <- prep$facet_names
  step_idx <- if (!is.null(step_facet)) {
    as.integer(df[[step_facet]])
  } else {
    NULL
  }
  slope_idx <- if (length(slope_facet) > 1L) {
    if (is.null(gpcm_spec$cells) || !identical(slope_facet, gpcm_spec$slope_facet)) {
      stop("Multiple slope families require their fitted crossing specification.")
    }
    if (!identical(prep$levels[slope_facet], gpcm_spec$component_levels)) {
      stop("Slope indices require the fitted facet level ordering.")
    }
    key <- function(x) do.call(paste, c(unname(lapply(x, as.integer)), sep = ":"))
    index <- match(key(df[slope_facet]), key(gpcm_spec$cells))
    if (anyNA(index)) stop("Slope indices require crossings in the supplied specification.")
    index
  } else if (!is.null(slope_facet)) {
    as.integer(df[[slope_facet]])
  } else {
    NULL
  }
  # Pre-split observation indices by criterion for PCM (avoids repeated which())
  criterion_splits <- if (!is.null(step_idx)) {
    split(seq_len(nrow(df)), step_idx)
  } else {
    NULL
  }
  slope_splits <- if (!is.null(slope_idx)) {
    split(seq_len(nrow(df)), slope_idx)
  } else {
    NULL
  }
  interaction_idx <- if (length(interaction_specs %||% list()) > 0L) {
    out <- lapply(interaction_specs, function(spec) {
      a_idx <- facets_idx[[spec$facet_a]]
      b_idx <- facets_idx[[spec$facet_b]]
      as.integer(a_idx + (b_idx - 1L) * spec$n_a)
    })
    names(out) <- names(interaction_specs)
    out
  } else {
    list()
  }
  list(
    person = as.integer(df$Person),
    facets = facets_idx,
    interactions = interaction_idx,
    step_idx = step_idx,
    slope_idx = slope_idx,
    criterion_splits = criterion_splits,
    slope_splits = slope_splits,
    score_k = as.integer(df$score_k),
    weight = suppressWarnings(as.numeric(df$Weight))
  )
}

sample_mfrm_data <- function(seed = 20240131) {
  with_preserved_rng_seed(seed, {
    persons <- paste0("P", sprintf("%02d", 1:36))
    raters <- paste0("R", 1:3)
    tasks <- paste0("T", 1:4)
    criteria <- paste0("C", 1:3)
    df <- expand_grid(
      Person = persons,
      Rater = raters,
      Task = tasks,
      Criterion = criteria
    )
    ability <- rnorm(length(persons), 0, 1)
    rater_eff <- c(-0.4, 0, 0.4)
    task_eff <- seq(-0.5, 0.5, length.out = length(tasks))
    crit_eff <- c(-0.3, 0, 0.3)
    eta <- ability[match(df$Person, persons)] -
      rater_eff[match(df$Rater, raters)] -
      task_eff[match(df$Task, tasks)] -
      crit_eff[match(df$Criterion, criteria)]
    raw <- eta + rnorm(nrow(df), 0, 0.6)
    score <- as.integer(cut(
      raw,
      breaks = c(-Inf, -1.0, -0.3, 0.3, 1.0, Inf),
      labels = 1:5
    ))
    df$Score <- score
    df
  })
}

with_preserved_rng_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed[1]))
  force(expr)
}

guess_col <- function(cols, patterns, fallback = 1) {
  if (length(cols) == 0) return(character(0))
  hit <- which(stringr::str_detect(tolower(cols), paste(patterns, collapse = "|")))
  if (length(hit) > 0) return(cols[hit[1]])
  cols[min(fallback, length(cols))]
}
