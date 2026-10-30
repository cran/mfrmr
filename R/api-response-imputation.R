#' Review multiple imputations of scores on assigned ratings
#'
#' Check several completed versions of a rating table before analyzing them.
#' Supply these completed data from an imputation model; this function checks
#' that observed scores, rating assignments and identifiers are preserved.
#' No person-by-facet grid is constructed. This function validates supplied
#' imputations; it does not choose or fit an imputation model.
#' Start with `review_mfrm_imputations()`, then use [fit_mfrm_imputed()] and
#' [pool_mfrm_imputed()] for eligible fixed-facet analyses. The older name
#' `mfrm_response_imputations()` is retained with its original `impute` argument.
#'
#' @param data Original long-format rating roster, including missing scores.
#' @param completed A list of at least two completed data frames, or a `mids`
#'   object from [mice::mice()] fitted to the long-format roster. Each completion
#'   must contain all original columns and rows. Row order may differ; rows
#'   are matched by `event_id`.
#' @param person,score,event_id Column names. `event_id` uniquely identifies a
#'   rating event, including repeated ratings of the same person and facets.
#' @param facets Nonempty character vector of facet column names.
#' @param impute_ids Character vector of event IDs explicitly selecting missing
#'   scores on assigned ratings. Observed scores cannot be selected.
#'   These are values from the `event_id` column, not a column name or a logical
#'   switch. For example, `c("E2", "E7")` selects those two rating events.
#' @param impute Compatibility name for `impute_ids`, used only by
#'   `mfrm_response_imputations()`. Supply one selection, not both argument names.
#' @param categories The full intended contiguous integer category vector,
#'   for example `0:4`. It is preserved across all completed analyses.
#' @param assigned Optional name of a complete logical column: `TRUE` denotes
#'   an assigned rating and `FALSE` an unassigned combination. Without this
#'   column every supplied row is declared assigned. Unassigned rows must have
#'   missing scores in the original and every completion.
#' @param imputation_model For a list of completions, the saved model or a
#'   nonempty list containing its specification, settings and diagnostics.
#'   Required so that the provenance is retained. For a `mids` input, that
#'   object is retained automatically; omit this argument.
#' @param missing How to handle assigned missing scores not selected by
#'   `impute_ids`: `"error"` (default) or explicit `"omit"`. Omitted events remain
#'   in the roster and every completion, with missing scores, and are excluded
#'   from each analysis. This choice is not a correction for nonresponse.
#'
#' @return An `mfrm_response_imputations` object containing the original
#'   `data`, aligned `completed` data sets, `imputation_model`, an `events`
#'   table with assignment/observation/imputation/omission status, a `support`
#'   table counting original observed and imputed scores by person/facet level,
#'   and `settings`. No imputations or failed analyses are silently discarded.
#'
#' @details The score must contain numeric integer category labels (numeric
#'   vectors, or character/factor labels such as `"0"`, `"1"`). Recoded sentinel
#'   missing values must already be `NA`. Identifiers, assignment indicators
#'   and observed values must not change. Only selected scores may be filled;
#'   missing auxiliary predictors may also be completed. For a `mids` input,
#'   its original data and score `where` selection must agree with this review.
#' Numeric identifiers must remain distinct when converted to text. Numeric
#' observed predictors are compared exactly, without rounding to text, and
#' must remain numeric in completions; integer/double storage may differ.
#'
#' A sparse assignment and a missing assigned score are different events.
#' Neither absent roster rows nor explicit unassigned rows are imputed.
#' An entirely imputed person or facet level is visible in `support`; its
#' analysis depends on the imputation model, not on observed ratings for that
#' level. Imputed links do not establish empirical connectedness.
#'
#' Proper multiple imputation must include uncertainty about missing values
#' and imputation parameters, reflect the ordinal score support and the
#' person/facet dependence, and be compatible with the intended analysis.
#' Include relevant observed predictors of nonresponse. An MAR analysis
#' requires an adequate conditional model; informative missingness beyond
#' observed predictors requires sensitivity analysis. Passing these software
#' checks does not establish MAR, model adequacy or interval coverage.
#'
#' For a wide-format or multilevel imputation model, reshape each completion
#' back to the original long roster using event identities and supply the
#' resulting list with the saved model. Do not average completed scores.
#' See `vignette("mfrmr-response-imputation", package = "mfrmr")` for a
#' joint RSM example with supplied posterior predictive completions, shared
#' Person draws, sampling diagnostics, direct observed-score inference and a
#' separate lower-score sensitivity analysis. The accompanying R/Stan script
#' regenerates that example; it is not a general imputation engine.
#' Under the same score model and ignorable missingness, observed-score MML
#' can directly estimate the fixed-facet target without completing scores.
#' MI does not create additional observed information.
#'
#' @seealso [fit_mfrm_imputed()], [pool_mfrm_imputed()], [mfrm_cluster_imputed()]
#' @examples
#' # Small supplied completions to illustrate input checks, not a fitted imputer.
#' ratings <- data.frame(Event = paste0("E", 1:4), Person = c("P1", "P1", "P2", "P2"),
#'   Rater = c("A", "B", "A", "B"), Score = c(0, NA, 1, 2))
#' first <- second <- ratings
#' first$Score[2] <- 1
#' second$Score[2] <- 2
#' reviewed <- review_mfrm_imputations(ratings, list(first, second),
#'   person = "Person", facets = "Rater", score = "Score", event_id = "Event",
#'   impute_ids = "E2", categories = 0:2,
#'   imputation_model = list(method = "Illustrative supplied completions"))
#' reviewed$events  # Only E2 is imputed; the three observed scores are retained.
#' summary(reviewed)  # Original observed support and imputation counts per level.
#' @name mfrm_response_imputations
#' @rdname mfrm_response_imputations
#' @export
review_mfrm_imputations <- function(data, completed, person, facets, score,
                                     event_id, impute_ids, categories,
                                     assigned = NULL, imputation_model = NULL,
                                     missing = c("error", "omit")) {
  missing <- match.arg(missing)
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)) ||
      anyNA(names(data)) || any(!nzchar(names(data)))) {
    stop("`data` must be a nonempty data frame with distinct column names.", call. = FALSE)
  }
  scalar_name <- function(x) is.character(x) && length(x) == 1L &&
    !is.na(x) && nzchar(x)
  if (!all(vapply(list(person, score, event_id), scalar_name, logical(1))) ||
      !is.character(facets) || !length(facets) || anyNA(facets) ||
      any(!nzchar(facets)) || (!is.null(assigned) && !scalar_name(assigned))) {
    stop("Supply column names for `person`, `facets`, `score`, `event_id` and optional `assigned`.", call. = FALSE)
  }
  columns <- c(person, facets, score, event_id, assigned)
  if (anyDuplicated(columns) || !all(columns %in% names(data))) {
    stop("Model, event and assignment columns must be distinct and present in `data`.", call. = FALSE)
  }
  if (!is.numeric(categories) || is.complex(categories) || length(categories) < 2L ||
      anyNA(categories) || any(!is.finite(categories)) ||
      any(categories != floor(categories)) || any(diff(categories) != 1) ||
      any(abs(categories) > .Machine$integer.max)) {
    stop("`categories` must be the full contiguous increasing integer score scale.", call. = FALSE)
  }
  data <- as.data.frame(data)
  ids <- as.character(data[[event_id]])
  id_columns <- c(event_id, person, facets)
  if (any(vapply(data[id_columns], function(z) anyNA(z) ||
      any(!nzchar(trimws(as.character(z)))), logical(1))) || anyDuplicated(ids)) {
    stop("Event IDs must be unique; event, person and facet IDs must be nonmissing and nonempty.", call. = FALSE)
  }
  for (nm in id_columns) {
    z <- data[[nm]]
    if (!(is.numeric(z) || is.character(z) || is.factor(z)) || is.complex(z) ||
        !is.null(dim(z)) || (is.numeric(z) && any(!is.finite(z))) ||
        length(unique(z)) != length(unique(as.character(z)))) {
      stop("Identifiers must be finite and retain distinct text labels: ", nm,
        ". Supply unambiguous character IDs.", call. = FALSE)
    }
  }
  is_assigned <- if (is.null(assigned)) rep(TRUE, nrow(data)) else data[[assigned]]
  if (!is.logical(is_assigned) || anyNA(is_assigned)) {
    stop("The `assigned` column must be complete and logical.", call. = FALSE)
  }
  original_score <- response_imputation_scores(data[[score]], categories)
  if (any(!is_assigned & !is.na(original_score))) {
    stop("Unassigned rows must have missing scores.", call. = FALSE)
  }
  if (!is.character(impute_ids) || !length(impute_ids) || anyNA(impute_ids) ||
      anyDuplicated(impute_ids) || !all(impute_ids %in% ids)) {
    stop("`impute_ids` must explicitly list distinct event IDs with missing assigned scores.", call. = FALSE)
  }
  selected <- ids %in% impute_ids
  if (any(selected & (!is_assigned | !is.na(original_score)))) {
    stop("Only missing scores on assigned ratings may be selected by `impute_ids`.", call. = FALSE)
  }
  omitted <- is_assigned & is.na(original_score) & !selected
  if (any(omitted) && missing == "error") {
    stop("Some assigned missing scores are not selected; review `impute_ids` or use missing = 'omit' explicitly.", call. = FALSE)
  }
  align <- function(z, original = FALSE) {
    if (!is.data.frame(z) || nrow(z) != nrow(data) || anyDuplicated(names(z)) ||
        !setequal(names(z), names(data)) || anyNA(z[[event_id]]) ||
        anyDuplicated(as.character(z[[event_id]])) ||
        !setequal(as.character(z[[event_id]]), ids)) {
      stop("Each data set must retain the original columns and all event IDs exactly once.", call. = FALSE)
    }
    z <- as.data.frame(z[match(ids, as.character(z[[event_id]])), names(data), drop = FALSE])
    rownames(z) <- NULL
    for (nm in setdiff(names(data), score)) {
      immutable <- nm %in% c(id_columns, assigned)
      rows <- if (immutable || original) rep(TRUE, nrow(data)) else !is.na(data[[nm]])
      before <- data[[nm]][rows]; after <- z[[nm]][rows]
      same <- if (is.numeric(before) || is.numeric(after)) {
        is.numeric(before) && is.numeric(after) &&
          isTRUE(all.equal(before, after, tolerance = 0, check.attributes = FALSE))
      } else identical(as.character(after), as.character(before))
      if (!same) {
        stop("Identifiers, assignment and observed predictor values must remain unchanged: ", nm, ".", call. = FALSE)
      }
    }
    values <- response_imputation_scores(z[[score]], categories)
    unchanged <- if (original) rep(TRUE, nrow(data)) else !selected
    if (!identical(values[unchanged], original_score[unchanged])) {
      stop("Observed and unselected scores must remain unchanged.", call. = FALSE)
    }
    if (!original && anyNA(values[selected])) {
      stop("Every selected missing score must be completed.", call. = FALSE)
    }
    z
  }
  if (inherits(completed, "mids")) {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("Install the optional 'mice' package to read mids completions.", call. = FALSE)
    }
    if (!is.null(imputation_model)) {
      stop("For a mids input, omit `imputation_model`; the mids object is retained.", call. = FALSE)
    }
    if (!is.numeric(completed$m) || length(completed$m) != 1L ||
        !is.finite(completed$m) || completed$m < 2 || completed$m != floor(completed$m)) {
      stop("A mids object must retain at least two imputations.", call. = FALSE)
    }
    align(completed$data, original = TRUE)
    where <- completed$where
    if (!is.matrix(where) || !is.logical(where) || anyNA(where) ||
        !identical(dim(where), dim(completed$data)) ||
        !identical(colnames(where), names(completed$data))) {
      stop("`completed$where` must record the mice imputation selection.", call. = FALSE)
    }
    where <- where[match(ids, as.character(completed$data[[event_id]])), , drop = FALSE]
    if (!identical(unname(where[, score]), selected) || any(where[, c(id_columns, assigned), drop = FALSE])) {
      stop("mice `where` must select exactly `impute_ids` scores and no identifier or assignment cells.", call. = FALSE)
    }
    imputation_model <- completed
    completed <- lapply(seq_len(completed$m), function(i) mice::complete(completed, action = i))
  }
  if (!is.list(completed) || is.data.frame(completed) || length(completed) < 2L) {
    stop("Supply at least two completed data frames or a mids object.", call. = FALSE)
  }
  if (!is.list(imputation_model) || !length(imputation_model)) {
    stop("Retain the saved imputation model or its specification and diagnostics in `imputation_model`.", call. = FALSE)
  }
  completed <- lapply(seq_along(completed), function(i) tryCatch(align(completed[[i]]),
    error = function(e) stop("Imputation ", i, ": ", conditionMessage(e), call. = FALSE)))
  events <- data.frame(ID = ids, Assigned = is_assigned,
    Observed = !is.na(original_score), Imputed = selected, Omitted = omitted,
    stringsAsFactors = FALSE)
  support <- do.call(rbind, lapply(c(person, facets), function(nm) {
    values <- as.character(data[[nm]])
    counts <- rowsum(1L * cbind(Assigned = is_assigned, Observed = events$Observed,
      Imputed = selected, Omitted = omitted), values, reorder = FALSE)
    data.frame(Facet = nm, Level = rownames(counts), counts, row.names = NULL)
  }))
  out <- list(data = data, completed = completed, imputation_model = imputation_model,
    events = events, support = support,
    settings = list(person = person, facets = facets, score = score, event_id = event_id,
      categories = as.integer(categories), assigned = assigned, missing = missing,
      imputations = length(completed)))
  class(out) <- "mfrm_response_imputations"
  out
}

#' @rdname mfrm_response_imputations
#' @export
mfrm_response_imputations <- function(data, completed, person, facets, score,
                                     event_id, impute, categories,
                                     assigned = NULL, imputation_model = NULL,
                                     missing = c("error", "omit")) {
  review_mfrm_imputations(data, completed, person, facets, score, event_id,
    impute_ids = impute, categories = categories, assigned = assigned,
    imputation_model = imputation_model, missing = match.arg(missing))
}

response_imputation_scores <- function(x, categories) {
  if (!(is.numeric(x) || is.character(x) || is.factor(x)) || is.complex(x)) {
    stop("Scores must use numeric integer category labels or NA.", call. = FALSE)
  }
  value <- suppressWarnings(as.numeric(as.character(x)))
  observed <- !is.na(x)
  if (any(!is.finite(value[observed])) || any(!value[observed] %in% categories)) {
    stop("All observed and imputed scores must belong to `categories`; no rounding or recoding is performed.", call. = FALSE)
  }
  value
}

#' @rdname mfrm_response_imputations
#' @param x,object An object returned by [review_mfrm_imputations()] or its
#'   compatibility wrapper [mfrm_response_imputations()].
#' @param ... Unused for print and summary methods.
#' @export
print.mfrm_response_imputations <- function(x, ...) {
  cat("Multiple imputations of assigned rating scores\n")
  cat(length(x$completed), "imputations;", sum(x$events$Imputed), "selected missing scores;",
      sum(!x$events$Assigned), "unassigned rows;", sum(x$events$Omitted), "omitted assigned rows\n")
  cat("Review the imputation model and original observed support before analysis.\n")
  invisible(x)
}

#' @rdname mfrm_response_imputations
#' @export
summary.mfrm_response_imputations <- function(object, ...) object$support
