# Match a requested response multiset, never a many-to-many group lookup.
.dff_match_observations <- function(fit, data, group, normalized = FALSE) {
  obs <- compute_obs_table(fit)
  facets <- fit$prep$facet_names
  source <- fit$prep$source_columns %||% fit$config$source_columns
  person <- if (normalized) "Person" else source$person %||% fit$config$person_col %||% "Person"
  score <- if (normalized) "Score" else source$score %||% fit$config$score_col %||% "Score"
  weight <- if (normalized) "Weight" else source$weight %||% fit$prep$weight_col
  required <- unique(c(person, facets, score, weight, group))
  if (anyDuplicated(names(data)) || !all(required %in% names(data))) stop(
    "Residual DFF data must contain unique columns for the fitted Person, all facets, score, weight (if used), and group.",
    call. = FALSE)
  data <- as.data.frame(data)
  missing_codes <- fit$config$replay_inputs$missing_codes
  if (!normalized && !is.null(missing_codes) && !isFALSE(missing_codes)) {
    defaults <- isTRUE(missing_codes) || identical(tolower(missing_codes), "default")
    data <- recode_missing_codes(data, columns = if (defaults) score else c(person, facets, score),
      codes = if (defaults) c("99", "999", "-1", "N", "NA", "n/a", ".", "") else missing_codes,
      verbose = FALSE)
  }
  ids <- c("Person", facets)
  observed <- as.data.frame(obs[c(ids, "Score", "Weight")])
  requested <- data[c(person, facets, score)]
  names(requested) <- c(ids, "Score")
  requested$Weight <- if (is.null(weight)) rep(1, nrow(data)) else data[[weight]]
  for (nm in ids) {
    observed[[nm]] <- enc2utf8(trimws(as.character(observed[[nm]])))
    requested[[nm]] <- enc2utf8(trimws(as.character(requested[[nm]])))
    requested[[nm]][!is.na(requested[[nm]]) & !nzchar(requested[[nm]])] <- NA_character_
  }
  if (!normalized && is.data.frame(fit$prep$score_map)) {
    map <- fit$prep$score_map
    observed$Score <- map$OriginalScore[match(observed$Score, map$InternalScore)]
  }
  y <- suppressWarnings(as.numeric(as.character(requested$Score)))
  valid_score <- is.finite(y) & abs(y - round(y)) <= 1e-6 &
    y <= .Machine$integer.max & y >= -.Machine$integer.max
  # Replay preparation's accepted integer conversion, including its tolerance.
  requested$Score <- NA_integer_
  requested$Score[valid_score] <- as.integer(y[valid_score])
  requested$Weight <- suppressWarnings(as.numeric(as.character(requested$Weight)))
  valid <- stats::complete.cases(requested) & is.finite(requested$Weight) & requested$Weight > 0
  gv <- trimws(as.character(data[[group]]))
  valid_group <- !is.na(gv) & nzchar(gv)
  audit <- data.frame(InputRow = seq_len(nrow(data)), FitRow = NA_integer_,
    GroupValue = gv, Included = FALSE,
    Reason = ifelse(!valid_group, "missing_group", ifelse(!valid, "invalid_or_missing_response", "not_in_fitted_responses")))
  take <- which(valid & valid_group)
  # Canonical allocation of exchangeable copies keeps floating-point sums stable.
  take <- take[order(gv[take], method = "radix")]
  # Synthetic key names avoid collisions with arbitrary user facet names.
  names(observed) <- names(requested) <- paste0("k", seq_len(ncol(observed)))
  keys <- names(observed)
  numbered <- function(x) {
    dplyr::group_by(x, dplyr::across(dplyr::everything())) |>
      dplyr::mutate(.occurrence = dplyr::row_number(), .copies = dplyr::n()) |>
      dplyr::ungroup()
  }
  reference <- numbered(observed)
  reference$.fit_row <- seq_len(nrow(obs))
  request <- numbered(requested[take, , drop = FALSE])
  request$.input_row <- take
  counts <- dplyr::distinct(reference[c(keys, ".copies")])
  check <- dplyr::left_join(request, counts, by = keys, suffix = c(".input", ".fit"))
  if (any(check$.copies.input > check$.copies.fit, na.rm = TRUE)) stop(
    "Residual DFF data request more copies of a response than the fitted data contain; check duplicated rows and group assignments.",
    call. = FALSE)
  # Equal-key fitted rows must be exchangeable for these residual summaries.
  moment_keys <- cbind(observed, Expected = obs$Expected, Var = obs$Var)
  moments <- dplyr::distinct(moment_keys)
  if (anyDuplicated(moments[keys])) stop(
    "Repeated fitted responses have different model moments; event-specific matching is required.", call. = FALSE)
  matched <- dplyr::left_join(request, reference[c(keys, ".occurrence", ".fit_row")],
    by = c(keys, ".occurrence"))
  audit$FitRow[take] <- matched$.fit_row
  audit$Included <- !is.na(audit$FitRow)
  audit$Reason[audit$Included] <- "included"
  unmatched <- sum(valid & valid_group & !audit$Included)
  if (unmatched) message("Excluded ", unmatched,
    " row(s) not found in the fitted responses; inspect `$cases` before comparing methods.")
  used <- which(audit$Included)
  if (!length(used)) stop("No requested response rows match the fitted data.", call. = FALSE)
  used <- used[order(audit$FitRow[used])]
  work <- obs[audit$FitRow[used], , drop = FALSE]
  for (nm in ids) work[[nm]] <- as.character(work[[nm]])
  work$.group_var <- audit$GroupValue[used]
  list(data = work, cases = audit)
}
