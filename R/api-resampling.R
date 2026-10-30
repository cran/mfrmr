#' Build an observed-data resampling specification
#'
#' @param data A long-format observed MFRM data set with distinct column names.
#'   The `.mfrm_*` metadata columns listed below are reserved and must not
#'   already be present. Numeric identifiers must retain distinct text labels.
#' @param person Person/respondent identifier column.
#' @param facets Non-person facet columns used by the target MFRM fit.
#' @param score Ordered score column.
#' @param strata Optional person-level stratification columns, for example a
#'   `Region` or L1 group column. Each person must have at most one unique
#'   stratum combination.
#' @param preserve_facets Optional facet columns whose level coverage should be
#'   reviewed and, when possible, topped up after the stratified person draw.
#'   A common choice is the rater facet.
#' @param design Resampling design. `"stratified_subsample"` samples persons
#'   without replacement inside each stratum. `"stratified_bootstrap"` samples
#'   persons with replacement inside each stratum and re-keys duplicate person
#'   instances in the returned data.
#' @param reps Number of resampling replicates to draw.
#' @param sample_fraction Fraction of persons to draw within each stratum when
#'   `sample_n = NULL`. Targets are rounded upward and respect
#'   `min_per_stratum`. The default is `0.5` for both designs; use `1` for a
#'   same-size bootstrap.
#' @param sample_n Optional target number of persons to draw per stratum. Supply
#'   either one scalar used for every stratum, or a named numeric vector whose
#'   names match `spec$target_plan$Stratum`. An unnamed vector of length equal
#'   to the number of strata follows that table's order. Targets are rounded
#'   upward and capped at available Persons when sampling without replacement.
#' @param replace Optional logical override for replacement. By default,
#'   replacement is `FALSE` for `"stratified_subsample"` and `TRUE` for
#'   `"stratified_bootstrap"`.
#' @param seed Optional seed used by [draw_mfrm_resamples()].
#' @param min_per_stratum Minimum target persons per represented stratum.
#' @param topup_preserve_facets Logical; if `TRUE`, add extra person clusters
#'   when possible to recover missing levels of `preserve_facets`. This changes
#'   the sampling distribution and can increase stratum sizes; it is not an
#'   ordinary bootstrap. Set `FALSE` for ordinary cluster-bootstrap draws.
#'
#' @details
#' This helper defines a resampling design for observed-data stability checks.
#' It is intentionally separate from [build_mfrm_sim_spec()] and
#' [simulate_mfrm_recovery()]. The full-data estimates used with these draws
#' are reference estimates, not known truth, so downstream summaries should be
#' described as estimation stability, reproducibility, or agreement with a
#' full-data reference rather than strict parameter recovery.
#'
#' The design is person-clustered: all rows for a selected person are kept
#' together. For bootstrap draws, duplicated person clusters are re-keyed in the
#' returned data while the original identifier is retained in
#' `.mfrm_original_person`.
#'
#' @section Sampling assumptions and interpretation:
#' For ordinary stratified Person-cluster bootstrap draws, set
#' `design = "stratified_bootstrap"`, `sample_fraction = 1`,
#' `sample_n = NULL`, and `topup_preserve_facets = FALSE`. The unit is a
#' Person, not a rating row; unequal-sized clusters retain all their rows on
#' every draw. This targets variation across independent Persons within strata,
#' conditional on the observed non-Person facet levels and assignment patterns.
#' It does not reproduce sampling new raters or tasks, or dependence across
#' Persons induced by shared random effects or peer-review networks. Crossed
#' random effects require a resampling scheme matched to that sampling target.
#'
#' A reduced sample size or coverage top-up instead defines a stability
#' experiment. Top-up chooses an additional Person carrying a currently absent
#' preserved level, and can also bring other missing levels. It changes
#' inclusion probabilities and stratum totals. Use `sampling_scheme`, the
#' returned terminology and base/added/final counts to report what was drawn.
#' These helpers do not fit models or construct confidence intervals; the
#' presence of the word bootstrap does not qualify downstream interval coverage.
#' Reference estimates from the full sample still have estimation uncertainty.
#'
#' @section Labels and saved specifications:
#' Original Person keys and data values are retained. Stratum labels escape
#' percent signs and vertical bars, and distinguish missing values, empty
#' strings and literal `"(missing)"`/`"(empty)"` labels. Inspect `target_plan`
#' before supplying named `sample_n` values. Data column names
#' `.mfrm_original_person`, `.mfrm_draw_unit`, `.mfrm_stratum`, `.mfrm_topup`,
#' and `.mfrm_topup_reason` are reserved for returned metadata.
#'
#' `draw_mfrm_resamples()` rebuilds derived tables from the specification's
#' declared inputs, including saved older specifications. Corrected labels and
#' top-up rules can change old draws even with the same seed; retain old results
#' separately. `saveRDS(spec, ...)` preserves the input data and seed;
#' `saveRDS(draws, ...)` preserves the actual samples, per-replicate seeds and
#' manifests. Use an explicit seed to reproduce draws in a matching R/RNG
#' environment. Manifest-only results retain the specification and input data
#' but omit the materialized samples.
#'
#' @return An object of class `mfrm_resampling_spec` containing the input data,
#'   columns, target plan, person/stratum rosters, `sampling_scheme` and
#'   interpretation in `terminology`.
#' @references
#' Owen, A. B. and Eckles, D. (2012). Bootstrapping data arrays of arbitrary
#' order. \emph{The Annals of Applied Statistics}, 6(3), 895--927.
#' \doi{10.1214/12-AOAS547}. This discusses crossed-factor resampling;
#' it is not the method implemented by this Person-cluster helper.
#' @seealso [draw_mfrm_resamples()], [build_mfrm_sim_spec()]
#' @examples
#' # Bundled synthetic ratings: inspect stability in half-samples by Group.
#' ratings <- load_mfrmr_data("example_operational")
#' spec <- build_mfrm_resampling_spec(
#'   ratings, person = "Person", facets = c("Rater", "Criterion"),
#'   score = "Score", strata = "Group", preserve_facets = "Rater",
#'   reps = 2, sample_fraction = 0.5, seed = 99
#' )
#' draws <- draw_mfrm_resamples(spec)
#' summary(draws)$overview
#' draws$stratum_manifest  # Base, top-up and final Person-cluster counts.
#' @export
build_mfrm_resampling_spec <- function(data,
                                       person,
                                       facets,
                                       score,
                                       strata = NULL,
                                       preserve_facets = NULL,
                                       design = c("stratified_subsample", "stratified_bootstrap"),
                                       reps = 20,
                                       sample_fraction = 0.5,
                                       sample_n = NULL,
                                       replace = NULL,
                                       seed = NULL,
                                       min_per_stratum = 1,
                                       topup_preserve_facets = TRUE) {
  design <- match.arg(tolower(as.character(design[1])),
                      c("stratified_subsample", "stratified_bootstrap"))
  if (!is.data.frame(data) || nrow(data) == 0L || anyDuplicated(names(data)) ||
      anyNA(names(data)) || any(!nzchar(names(data)))) {
    stop("`data` must be a non-empty data.frame with unique, nonempty column names.", call. = FALSE)
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  reserved <- c(".mfrm_original_person", ".mfrm_draw_unit", ".mfrm_stratum",
    ".mfrm_topup", ".mfrm_topup_reason")
  if (any(reserved %in% names(data))) stop(
    "`data` contains reserved resampling columns: ",
    paste(intersect(reserved, names(data)), collapse = ", "),
    ". Rename them or supply the original input columns.", call. = FALSE)
  person <- resampling_scalar_name(person, "person")
  score <- resampling_scalar_name(score, "score")
  facets <- resampling_name_vector(facets, "facets", min_len = 1L)
  strata <- resampling_name_vector(strata, "strata", min_len = 0L)
  preserve_facets <- resampling_name_vector(preserve_facets, "preserve_facets", min_len = 0L)

  required <- unique(c(person, facets, score, strata, preserve_facets))
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop("`data` is missing required columns: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  for (name in unique(c(person, facets, strata, preserve_facets))) {
    value <- data[[name]]
    if (!(is.character(value) || is.factor(value) || is.logical(value) ||
          (is.numeric(value) && !is.object(value) && !is.complex(value))) ||
        !is.null(dim(value))) stop(
      "Resampling identifiers/strata must be character, factor, logical or numeric vectors: ",
      name, ".", call. = FALSE)
    if (is.numeric(value) && (any(!is.finite(value) & !is.na(value)) ||
        length(unique(value[!is.na(value)])) != length(unique(as.character(value[!is.na(value)]))))) stop(
      "Numeric identifiers/strata must be finite and retain distinct text labels: ",
      name, ". Supply distinct character labels.", call. = FALSE)
  }
  if (any(is.na(data[[person]]) | !nzchar(as.character(data[[person]])))) {
    stop("`data` contains missing or empty person identifiers.", call. = FALSE)
  }

  reps <- as.integer(reps[1])
  if (!is.finite(reps) || reps < 1L) {
    stop("`reps` must be a positive integer.", call. = FALSE)
  }
  min_per_stratum <- as.integer(min_per_stratum[1])
  if (!is.finite(min_per_stratum) || min_per_stratum < 1L) {
    stop("`min_per_stratum` must be a positive integer.", call. = FALSE)
  }
  if (is.null(replace)) {
    replace <- identical(design, "stratified_bootstrap")
  } else {
    replace <- isTRUE(replace[1])
  }
  if (identical(design, "stratified_bootstrap") && !isTRUE(replace)) {
    stop("`design = \"stratified_bootstrap\"` requires `replace = TRUE`.", call. = FALSE)
  }
  if (!is.null(sample_fraction)) {
    sample_fraction <- suppressWarnings(as.numeric(sample_fraction[1]))
    if (!is.finite(sample_fraction) || sample_fraction <= 0 || sample_fraction > 1) {
      stop("`sample_fraction` must be a finite number in (0, 1].", call. = FALSE)
    }
  }
  if (is.null(sample_n) && is.null(sample_fraction)) {
    stop("Supply either `sample_fraction` or `sample_n`.", call. = FALSE)
  }

  person_table <- resampling_person_table(data, person = person, strata = strata)
  stratum_overview <- resampling_stratum_overview(
    data = data,
    person = person,
    person_table = person_table
  )
  target_plan <- resampling_target_plan(
    stratum_overview = stratum_overview,
    sample_fraction = sample_fraction,
    sample_n = sample_n,
    min_per_stratum = min_per_stratum,
    replace = replace
  )
  preserve_overview <- resampling_preserve_overview(
    data = data,
    preserve_facets = preserve_facets
  )
  sampling_scheme <- if (!replace) "stratified_person_subsample" else if (
    all(target_plan$TargetPersons == target_plan$AvailablePersons)) {
    "stratified_person_bootstrap"
  } else "stratified_person_resample_with_replacement"
  if (isTRUE(topup_preserve_facets) && length(preserve_facets))
    sampling_scheme <- paste0("coverage_augmented_", sampling_scheme)

  spec <- list(
    data = data,
    columns = list(
      person = person,
      facets = facets,
      score = score,
      strata = strata,
      preserve_facets = preserve_facets
    ),
    design = design,
    sampling_scheme = sampling_scheme,
    reps = reps,
    sample_fraction = sample_fraction,
    sample_n = sample_n,
    replace = replace,
    seed = seed,
    min_per_stratum = min_per_stratum,
    topup_preserve_facets = isTRUE(topup_preserve_facets),
    person_table = person_table,
    stratum_overview = stratum_overview,
    target_plan = target_plan,
    preserve_overview = preserve_overview,
    terminology = c(
      "Observed-data resampling checks stability and reproducibility against a full-data reference.",
      "Full-data estimates are reference estimates, not known true parameters.",
      "Use parameter-recovery wording only for simulations with known generated truth.",
      "The resampling unit is the entire Person cluster; inference requires independent Persons within strata, conditional on the observed facet levels and assignment patterns.",
      "Person resampling does not reproduce sampling new raters/tasks or dependence shared across Persons, including peer-review networks.",
      "A same-size stratified Person bootstrap uses replacement, the original number of Persons in each stratum and no coverage top-up. Other sample sizes and coverage top-up are stability designs, not automatic bootstrap confidence-interval procedures.",
      paste0("Sampling scheme: ", sampling_scheme, ". Top-up can change stratum sizes and inclusion probabilities; inspect base and added cluster counts.")
    )
  )
  class(spec) <- "mfrm_resampling_spec"
  spec
}

#' Draw observed-data MFRM resamples
#'
#' @param spec Output from [build_mfrm_resampling_spec()].
#' @param keep_data Logical; if `TRUE`, return a list of replicate data frames.
#'   If `FALSE`, return manifest tables only.
#'
#' @details See [build_mfrm_resampling_spec()] for sampling assumptions,
#'   top-up, label handling and saved-specification compatibility. Sampling
#'   keeps each selected Person's entire rating cluster and re-keys every
#'   instance when drawing with replacement. Other facet labels are retained.
#'
#' @return An object of class `mfrm_resamples` with:
#' - `samples`: replicate data frames with original Person, draw-unit, stratum
#'   and top-up metadata, or an empty list when `keep_data = FALSE`.
#' - `manifest`: replicate seeds, sampling scheme, base/added/final cluster
#'   counts, original Person counts, rating-row counts and coverage checks.
#' - `stratum_manifest`: `BasePersonClusters` before top-up,
#'   `TopupPersonClusters`, final `SelectedPersonClusters` and
#'   `SelectedOriginalPersons`, and final `Rows` for each stratum.
#' - `preserve_manifest`: observed coverage of requested facet levels.
#' - `spec` and `terminology`: inputs, sampling assumptions and interpretation.
#' `summary()` retains full manifests as well as compact tables. For older
#' results with top-up, pre-top-up stratum counts are retained as base counts;
#' unrecorded final cluster counts are `NA`, with a note, while final row counts
#' remain available. This does not redraw old samples.
#' @seealso [build_mfrm_resampling_spec()]
#' @examples
#' # Same-size Person-cluster bootstrap within the observed Group strata.
#' ratings <- load_mfrmr_data("example_operational")
#' spec <- build_mfrm_resampling_spec(
#'   ratings, person = "Person", facets = c("Rater", "Criterion"),
#'   score = "Score", strata = "Group", reps = 2,
#'   design = "stratified_bootstrap", sample_fraction = 1,
#'   topup_preserve_facets = FALSE, seed = 99
#' )
#' draws <- draw_mfrm_resamples(spec, keep_data = FALSE)
#' summary(draws)$overview
#' @export
draw_mfrm_resamples <- function(spec, keep_data = TRUE) {
  if (!inherits(spec, "mfrm_resampling_spec")) {
    stop("`spec` must be output from build_mfrm_resampling_spec().", call. = FALSE)
  }
  # Rebuild derived rosters so old saved specs cannot retain merged identities.
  controls <- setdiff(names(formals(build_mfrm_resampling_spec)), c("data", names(spec$columns)))
  spec <- do.call(build_mfrm_resampling_spec,
    c(list(data = spec$data), spec$columns, spec[intersect(controls, names(spec))]))
  keep_data <- isTRUE(keep_data)

  with_preserved_rng_seed(spec$seed, {
    rep_seeds <- sample.int(.Machine$integer.max, size = spec$reps)
  })

  samples <- vector("list", spec$reps)
  manifest <- vector("list", spec$reps)
  stratum_manifest <- vector("list", spec$reps)
  preserve_manifest <- vector("list", spec$reps)

  for (rep_id in seq_len(spec$reps)) {
    draw <- with_preserved_rng_seed(rep_seeds[rep_id], {
      resampling_draw_one(spec, rep_id = rep_id, rep_seed = rep_seeds[rep_id])
    })
    if (keep_data) {
      samples[[rep_id]] <- draw$data
      names(samples)[rep_id] <- draw$resample_id
    }
    manifest[[rep_id]] <- draw$manifest
    stratum_manifest[[rep_id]] <- draw$stratum_manifest
    preserve_manifest[[rep_id]] <- draw$preserve_manifest
  }

  out <- list(
    spec = spec,
    samples = if (keep_data) samples else list(),
    manifest = dplyr::bind_rows(manifest),
    stratum_manifest = dplyr::bind_rows(stratum_manifest),
    preserve_manifest = dplyr::bind_rows(preserve_manifest),
    terminology = spec$terminology
  )
  class(out) <- "mfrm_resamples"
  out
}

#' @export
summary.mfrm_resampling_spec <- function(object, ...) {
  overview <- data.frame(
    Design = object$design,
    SamplingScheme = object$sampling_scheme %||% "legacy_scheme_unreviewed",
    Reps = object$reps,
    Replace = isTRUE(object$replace),
    Persons = length(unique(object$person_table[[object$columns$person]])),
    Rows = nrow(object$data),
    Strata = nrow(object$stratum_overview),
    PreserveFacets = paste(object$columns$preserve_facets, collapse = ", "),
    stringsAsFactors = FALSE
  )
  out <- list(
    overview = overview,
    columns = data.frame(
      Role = c("person", "score", rep("facet", length(object$columns$facets)),
               rep("stratum", length(object$columns$strata)),
               rep("preserve_facet", length(object$columns$preserve_facets))),
      Column = c(object$columns$person, object$columns$score, object$columns$facets,
                 object$columns$strata, object$columns$preserve_facets),
      stringsAsFactors = FALSE
    ),
    stratum_overview = object$stratum_overview,
    target_plan = object$target_plan,
    preserve_overview = object$preserve_overview,
    terminology = object$terminology
  )
  class(out) <- "summary.mfrm_resampling_spec"
  out
}

#' @export
print.mfrm_resampling_spec <- function(x, ...) {
  s <- summary(x)
  cat("MFRM Observed-Data Resampling Specification\n")
  print(s$overview, row.names = FALSE)
  if (nrow(s$target_plan) > 0L) {
    cat("\nStratum target plan\n")
    print(s$target_plan, row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.summary.mfrm_resampling_spec <- function(x, ...) {
  cat("MFRM Observed-Data Resampling Specification Summary\n")
  print(x$overview, row.names = FALSE)
  if (nrow(x$stratum_overview) > 0L) {
    cat("\nStrata\n")
    print(x$stratum_overview, row.names = FALSE)
  }
  if (nrow(x$preserve_overview) > 0L) {
    cat("\nPreserve-facet coverage basis\n")
    print(x$preserve_overview, row.names = FALSE)
  }
  cat("\nTerminology\n")
  for (line in x$terminology) cat(" - ", line, "\n", sep = "")
  invisible(x)
}

#' @export
summary.mfrm_resamples <- function(object, top_n = 10, ...) {
  manifest <- as.data.frame(object$manifest %||% data.frame(), stringsAsFactors = FALSE)
  stratum_manifest <- as.data.frame(object$stratum_manifest %||% data.frame(), stringsAsFactors = FALSE)
  preserve_manifest <- as.data.frame(object$preserve_manifest %||% data.frame(), stringsAsFactors = FALSE)
  terminology <- object$terminology
  if (nrow(stratum_manifest) && !"BasePersonClusters" %in% names(stratum_manifest)) {
    stratum_manifest$BasePersonClusters <- stratum_manifest$SelectedPersonClusters
    affected <- stratum_manifest$Replicate %in% manifest$Replicate[manifest$TopupPersonClusters > 0]
    stratum_manifest$TopupPersonClusters <- ifelse(affected, NA_integer_, 0L)
    stratum_manifest$SelectedPersonClusters[affected] <- NA_integer_
    stratum_manifest$SelectedOriginalPersons[affected] <- NA_integer_
    if (any(affected)) terminology <- c(terminology,
      "Legacy stratum counts describe the base draw only. Final cluster counts by stratum are unavailable where top-up was used; original base counts and final row counts are retained.")
  }

  coverage_rate <- if (nrow(manifest) == 0L) NA_real_ else {
    mean(manifest$PreserveCoverageComplete %in% TRUE, na.rm = TRUE)
  }
  topup_reps <- if (nrow(manifest) == 0L) NA_integer_ else {
    sum(suppressWarnings(as.integer(manifest$TopupPersonClusters)) > 0L, na.rm = TRUE)
  }
  gap_reps <- if (nrow(manifest) == 0L) NA_integer_ else {
    sum(!manifest$FallbackStatus %in% c("ok", "topup_used"), na.rm = TRUE)
  }
  overview <- data.frame(
    Design = as.character(object$spec$design),
    SamplingScheme = object$spec$sampling_scheme %||% "legacy_scheme_unreviewed",
    Reps = nrow(manifest),
    SamplesReturned = length(object$samples),
    Replace = isTRUE(object$spec$replace),
    PreserveCoverageCompleteRate = coverage_rate,
    TopupReps = topup_reps,
    GapOrFallbackReps = gap_reps,
    stringsAsFactors = FALSE
  )

  stratum_summary <- if (nrow(stratum_manifest) == 0L) {
    data.frame()
  } else {
    stratum_manifest |>
      dplyr::group_by(.data$Stratum) |>
      dplyr::summarise(
        AvailablePersons = dplyr::first(.data$AvailablePersons),
        MeanTargetPersons = mean(.data$TargetPersons, na.rm = TRUE),
        MeanSelectedPersonClusters = mean(.data$SelectedPersonClusters),
        MeanBasePersonClusters = if ("BasePersonClusters" %in% names(stratum_manifest))
          mean(.data$BasePersonClusters) else NA_real_,
        MeanTopupPersonClusters = if ("TopupPersonClusters" %in% names(stratum_manifest))
          mean(.data$TopupPersonClusters) else NA_real_,
        MeanRows = mean(.data$Rows, na.rm = TRUE),
        Reps = dplyr::n(),
        .groups = "drop"
      ) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  preserve_summary <- if (nrow(preserve_manifest) == 0L) {
    data.frame()
  } else {
    preserve_manifest |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(
        AvailableLevels = dplyr::first(.data$AvailableLevels),
        MeanCoveredLevels = mean(.data$CoveredLevels, na.rm = TRUE),
        CompleteCoverageRate = mean(.data$MissingLevels == 0L, na.rm = TRUE),
        Reps = dplyr::n(),
        .groups = "drop"
      ) |>
      as.data.frame(stringsAsFactors = FALSE)
  }

  out <- list(
    overview = overview,
    manifest = utils::head(manifest, n = top_n),
    stratum_summary = stratum_summary,
    preserve_summary = preserve_summary,
    terminology = terminology,
    full_manifest = manifest,
    full_stratum_manifest = stratum_manifest,
    full_preserve_manifest = preserve_manifest
  )
  class(out) <- "summary.mfrm_resamples"
  out
}

#' @rawNamespace S3method(print,mfrm_resamples,print.mfrm_bundle)
NULL

#' @export
print.summary.mfrm_resamples <- function(x, ...) {
  cat("MFRM Observed-Data Resampling Draws Summary\n")
  print(x$overview, row.names = FALSE)
  if (nrow(x$stratum_summary) > 0L) {
    cat("\nStratum summary\n")
    print(x$stratum_summary, row.names = FALSE)
  }
  if (nrow(x$preserve_summary) > 0L) {
    cat("\nPreserve-facet summary\n")
    print(x$preserve_summary, row.names = FALSE)
  }
  cat("\nTerminology\n")
  for (line in x$terminology) cat(" - ", line, "\n", sep = "")
  invisible(x)
}

resampling_scalar_name <- function(x, arg_name) {
  value <- as.character(x[1] %||% NA_character_)
  if (is.na(value) || !nzchar(value)) {
    stop("`", arg_name, "` must be a non-empty column name.", call. = FALSE)
  }
  value
}

resampling_name_vector <- function(x, arg_name, min_len = 1L) {
  if (is.null(x)) {
    out <- character(0)
  } else {
    out <- unique(as.character(x))
    out <- out[!is.na(out) & nzchar(out)]
  }
  if (length(out) < min_len) {
    stop("`", arg_name, "` must contain at least ", min_len,
         " non-empty column name", ifelse(min_len == 1L, "", "s"), ".",
         call. = FALSE)
  }
  out
}

resampling_missing_label <- function(x) {
  missing <- is.na(x)
  x <- as.character(x)
  # Escape separators and sentinel literals before constructing tuple labels.
  x <- gsub("%", "%25", x, fixed = TRUE)
  x <- gsub("|", "%7C", x, fixed = TRUE)
  x[!missing & x == "(missing)"] <- "%28missing%29"
  x[!missing & x == "(empty)"] <- "%28empty%29"
  x[!missing & x == ""] <- "(empty)"
  x[missing] <- "(missing)"
  x
}

resampling_person_table <- function(data, person, strata) {
  person_values <- as.character(data[[person]])
  if (length(strata) == 0L) {
    tbl <- tibble::tibble(.mfrm_person = unique(person_values))
    names(tbl)[names(tbl) == ".mfrm_person"] <- person
    tbl$.mfrm_stratum <- "all"
    return(tbl)
  }

  tbl <- data[, unique(c(person, strata)), drop = FALSE]
  tbl[[person]] <- person_values
  # The Person key remains unchanged; only stratum display labels are escaped.
  stratum_values <- lapply(data[, strata, drop = FALSE], resampling_missing_label)
  tbl$.mfrm_stratum <- do.call(paste, c(unname(stratum_values), sep = " | "))
  tbl <- tbl[!duplicated(tbl[, c(person, ".mfrm_stratum"), drop = FALSE]), , drop = FALSE]
  if (anyDuplicated(tbl[[person]])) {
    bad <- unique(tbl[[person]][duplicated(tbl[[person]])])
    stop(
      "Each person must map to at most one stratum combination. Conflicting persons: ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5L) ", ..." else "",
      ".",
      call. = FALSE
    )
  }
  tibble::as_tibble(tbl)
}

resampling_stratum_overview <- function(data, person, person_table) {
  counts <- data.frame(
    Stratum = person_table$.mfrm_stratum,
    Rows = tabulate(match(as.character(data[[person]]), person_table[[person]]),
      nbins = nrow(person_table)),
    stringsAsFactors = FALSE
  )
  out <- counts |>
    dplyr::group_by(.data$Stratum) |>
    dplyr::summarise(
      Persons = dplyr::n(),
      Rows = sum(.data$Rows, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$Stratum)
  as.data.frame(out, stringsAsFactors = FALSE)
}

resampling_target_plan <- function(stratum_overview,
                                   sample_fraction,
                                   sample_n,
                                   min_per_stratum,
                                   replace) {
  strata <- as.character(stratum_overview$Stratum)
  available <- as.integer(stratum_overview$Persons)
  if (is.null(sample_n)) {
    target <- ceiling(available * sample_fraction)
  } else {
    sample_n_num <- suppressWarnings(as.numeric(sample_n))
    if (any(!is.finite(sample_n_num) | sample_n_num < 1L)) {
      stop("`sample_n` must contain positive finite values.", call. = FALSE)
    }
    if (!is.null(names(sample_n)) && any(nzchar(names(sample_n)))) {
      nm <- names(sample_n)
      if (!all(strata %in% nm)) {
        stop("Named `sample_n` must include every stratum label: ",
             paste(strata, collapse = ", "), ".", call. = FALSE)
      }
      target <- sample_n_num[match(strata, nm)]
    } else if (length(sample_n_num) == 1L) {
      target <- rep(sample_n_num, length(strata))
    } else if (length(sample_n_num) == length(strata)) {
      target <- sample_n_num
    } else {
      stop("`sample_n` must be length 1, length equal to strata, or named by strata.",
           call. = FALSE)
    }
    target <- ceiling(target)
  }
  target <- pmax(as.integer(target), min_per_stratum)
  status <- rep("ok", length(target))
  if (!isTRUE(replace)) {
    capped <- target > available
    status[capped] <- "capped_at_available_persons"
    target <- pmin(target, available)
  }
  data.frame(
    Stratum = strata,
    AvailablePersons = available,
    TargetPersons = as.integer(target),
    Replace = isTRUE(replace),
    TargetStatus = status,
    stringsAsFactors = FALSE
  )
}

resampling_preserve_overview <- function(data, preserve_facets) {
  if (length(preserve_facets) == 0L) {
    return(data.frame())
  }
  rows <- lapply(preserve_facets, function(facet) {
    levels <- sort(unique(resampling_missing_label(data[[facet]])))
    data.frame(
      Facet = facet,
      AvailableLevels = length(levels),
      Levels = paste(levels, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

resampling_draw_one <- function(spec, rep_id, rep_seed) {
  person_col <- spec$columns$person
  data <- spec$data
  target_plan <- spec$target_plan
  person_table <- spec$person_table

  selected <- vector("list", nrow(target_plan))
  stratum_rows <- vector("list", nrow(target_plan))
  for (i in seq_len(nrow(target_plan))) {
    stratum <- target_plan$Stratum[i]
    pool <- person_table[person_table$.mfrm_stratum == stratum, , drop = FALSE]
    target <- as.integer(target_plan$TargetPersons[i])
    pick_idx <- sample.int(nrow(pool), size = target, replace = isTRUE(spec$replace))
    unit_tbl <- pool[pick_idx, , drop = FALSE]
    unit_tbl$.mfrm_topup <- FALSE
    unit_tbl$.mfrm_topup_reason <- NA_character_
    selected[[i]] <- unit_tbl
    stratum_rows[[i]] <- data.frame(
      Stratum = stratum,
      AvailablePersons = as.integer(nrow(pool)),
      TargetPersons = target,
      BasePersonClusters = as.integer(nrow(unit_tbl)),
      Rows = NA_integer_,
      TargetStatus = target_plan$TargetStatus[i],
      stringsAsFactors = FALSE
    )
  }
  selected_units <- dplyr::bind_rows(selected)
  selected_units <- resampling_topup_units(spec, selected_units)
  selected_units$.mfrm_draw_unit <- seq_len(nrow(selected_units))

  sample_data <- resampling_materialize_sample(spec, selected_units)
  stratum_manifest <- dplyr::bind_rows(stratum_rows)
  unit_by_stratum <- selected_units |>
    dplyr::group_by(.data$.mfrm_stratum) |>
    dplyr::summarise(SelectedPersonClusters = dplyr::n(),
      SelectedOriginalPersons = dplyr::n_distinct(.data[[person_col]]),
      TopupPersonClusters = sum(.data$.mfrm_topup %in% TRUE), .groups = "drop")
  stratum_manifest <- dplyr::left_join(stratum_manifest, unit_by_stratum,
    by = c("Stratum" = ".mfrm_stratum"))
  row_by_stratum <- sample_data |>
    dplyr::count(.data$.mfrm_stratum, name = "Rows")
  stratum_manifest <- stratum_manifest |>
    dplyr::left_join(row_by_stratum, by = c("Stratum" = ".mfrm_stratum"), suffix = c("", ".sample")) |>
    dplyr::mutate(Rows = dplyr::coalesce(.data$Rows.sample, .data$Rows)) |>
    dplyr::select(-dplyr::any_of("Rows.sample"))

  preserve_manifest <- resampling_preserve_manifest(spec, sample_data)
  missing_strata <- setdiff(spec$stratum_overview$Stratum, unique(sample_data$.mfrm_stratum))
  preserve_complete <- if (nrow(preserve_manifest) == 0L) TRUE else {
    all(preserve_manifest$MissingLevels == 0L)
  }
  status_parts <- unique(c(
    if (length(missing_strata) > 0L) "stratum_gap" else character(0),
    if (!isTRUE(preserve_complete)) "preserve_facet_gap" else character(0),
    if (any(selected_units$.mfrm_topup %in% TRUE)) "topup_used" else character(0),
    unique(stratum_manifest$TargetStatus[stratum_manifest$TargetStatus != "ok"])
  ))
  fallback_status <- if (length(status_parts) == 0L) "ok" else paste(status_parts, collapse = ";")
  resample_id <- paste0("resample_", sprintf("%03d", rep_id))

  manifest <- data.frame(
    ResampleID = resample_id,
    Replicate = as.integer(rep_id),
    Design = spec$design,
    SamplingScheme = spec$sampling_scheme,
    Replace = isTRUE(spec$replace),
    RepSeed = as.integer(rep_seed),
    PersonClusters = as.integer(nrow(selected_units)),
    BasePersonClusters = sum(!selected_units$.mfrm_topup),
    OriginalPersons = length(unique(selected_units[[person_col]])),
    Rows = nrow(sample_data),
    StrataRepresented = length(unique(sample_data$.mfrm_stratum)),
    StrataTotal = nrow(spec$stratum_overview),
    MissingStrata = paste(missing_strata, collapse = ", "),
    PreserveCoverageComplete = isTRUE(preserve_complete),
    MissingPreserveLevels = resampling_missing_preserve_label(preserve_manifest),
    TopupPersonClusters = sum(selected_units$.mfrm_topup %in% TRUE),
    FallbackStatus = fallback_status,
    stringsAsFactors = FALSE
  )
  stratum_manifest$ResampleID <- resample_id
  stratum_manifest$Replicate <- as.integer(rep_id)
  if (nrow(preserve_manifest) > 0L) {
    preserve_manifest$ResampleID <- resample_id
    preserve_manifest$Replicate <- as.integer(rep_id)
  } else {
    preserve_manifest$ResampleID <- character(0)
    preserve_manifest$Replicate <- integer(0)
  }

  list(
    resample_id = resample_id,
    data = sample_data,
    manifest = manifest,
    stratum_manifest = stratum_manifest[, c("ResampleID", "Replicate", setdiff(names(stratum_manifest), c("ResampleID", "Replicate")))],
    preserve_manifest = preserve_manifest[, c("ResampleID", "Replicate", setdiff(names(preserve_manifest), c("ResampleID", "Replicate")))]
  )
}

resampling_topup_units <- function(spec, selected_units) {
  if (!isTRUE(spec$topup_preserve_facets) || length(spec$columns$preserve_facets) == 0L) {
    return(selected_units)
  }
  person_col <- spec$columns$person
  data <- spec$data
  selected_people <- as.character(selected_units[[person_col]])
  selected_rows <- data[as.character(data[[person_col]]) %in% selected_people, , drop = FALSE]

  additions <- list()
  for (facet in spec$columns$preserve_facets) {
    all_levels <- sort(unique(resampling_missing_label(data[[facet]])))
    selected_levels <- sort(unique(resampling_missing_label(selected_rows[[facet]])))
    missing_levels <- setdiff(all_levels, selected_levels)
    if (length(missing_levels) == 0L) next

    for (level in missing_levels) {
      if (level %in% resampling_missing_label(selected_rows[[facet]])) next
      candidate_people <- unique(as.character(data[[person_col]][resampling_missing_label(data[[facet]]) == level]))
      if (!isTRUE(spec$replace)) {
        candidate_people <- setdiff(candidate_people, selected_people)
      }
      if (length(candidate_people) == 0L) next
      chosen <- sample(candidate_people, size = 1L)
      add <- spec$person_table[as.character(spec$person_table[[person_col]]) == chosen, , drop = FALSE][1, , drop = FALSE]
      add$.mfrm_topup <- TRUE
      add$.mfrm_topup_reason <- paste0(facet, "=", level)
      additions[[length(additions) + 1L]] <- add
      selected_people <- c(selected_people, chosen)
      selected_rows <- data[as.character(data[[person_col]]) %in% selected_people, , drop = FALSE]
    }
  }
  if (length(additions) == 0L) {
    return(selected_units)
  }
  dplyr::bind_rows(selected_units, dplyr::bind_rows(additions))
}

resampling_materialize_sample <- function(spec, selected_units) {
  person_col <- spec$columns$person
  data <- spec$data
  rows <- vector("list", nrow(selected_units))
  for (i in seq_len(nrow(selected_units))) {
    original_person <- as.character(selected_units[[person_col]][i])
    sub <- data[as.character(data[[person_col]]) == original_person, , drop = FALSE]
    sub$.mfrm_original_person <- original_person
    sub$.mfrm_draw_unit <- as.integer(selected_units$.mfrm_draw_unit[i])
    sub$.mfrm_stratum <- as.character(selected_units$.mfrm_stratum[i])
    sub$.mfrm_topup <- selected_units$.mfrm_topup[i] %in% TRUE
    sub$.mfrm_topup_reason <- as.character(selected_units$.mfrm_topup_reason[i] %||% NA_character_)
    if (isTRUE(spec$replace)) {
      sub[[person_col]] <- paste0(original_person, "__b", sprintf("%03d", selected_units$.mfrm_draw_unit[i]))
    }
    rows[[i]] <- sub
  }
  dplyr::bind_rows(rows)
}

resampling_preserve_manifest <- function(spec, sample_data) {
  if (length(spec$columns$preserve_facets) == 0L) {
    return(data.frame(
      Facet = character(),
      AvailableLevels = integer(),
      CoveredLevels = integer(),
      MissingLevels = integer(),
      MissingLevelNames = character(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(spec$columns$preserve_facets, function(facet) {
    all_levels <- sort(unique(resampling_missing_label(spec$data[[facet]])))
    covered <- sort(unique(resampling_missing_label(sample_data[[facet]])))
    missing <- setdiff(all_levels, covered)
    data.frame(
      Facet = facet,
      AvailableLevels = length(all_levels),
      CoveredLevels = length(intersect(all_levels, covered)),
      MissingLevels = length(missing),
      MissingLevelNames = paste(missing, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

resampling_missing_preserve_label <- function(preserve_manifest) {
  if (!is.data.frame(preserve_manifest) || nrow(preserve_manifest) == 0L) {
    return("")
  }
  bad <- preserve_manifest[preserve_manifest$MissingLevels > 0L, , drop = FALSE]
  if (nrow(bad) == 0L) {
    return("")
  }
  paste(paste0(bad$Facet, ":", bad$MissingLevelNames), collapse = "; ")
}
