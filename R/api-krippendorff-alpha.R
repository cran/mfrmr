#' Compute Krippendorff's alpha across raters
#'
#' Estimate observed-score agreement across two or more raters, allowing
#' unequal numbers of ratings per unit. This is Krippendorff's agreement
#' coefficient, not Cronbach's internal-consistency alpha or MFRM separation
#' reliability. No measurement-model fitting is performed.
#'
#' @param data Long-format data with at most one row per rater and unit.
#' @param rater Name of the rater identifier column.
#' @param unit Names of columns jointly identifying a rated unit, for example
#'   `c("Person", "Task")`. They match observations, not variance components.
#' @param score Name of the score column. Missing scores are `NA`.
#'   Nominal/ordinal scores may be numeric, character, factor or logical;
#'   factors compare by labels. Interval/ratio scores must be finite numeric
#'   values, and ratio scores must be nonnegative with a meaningful zero.
#' @param level Measurement level: `"nominal"` (default), `"ordinal"`,
#'   `"interval"` or `"ratio"`. Choose on substantive grounds; numeric
#'   category codes alone do not establish an interval or ratio scale.
#' @param categories For ordinal scores, the required distinct category values
#'   in their substantive order. Factor levels and numeric spacing never
#'   silently set this order. For nominal scores, optional permitted values;
#'   if `NULL`, observed nonmissing values are used. Values must have the same
#'   numeric/character/logical representation as `score` (factors use labels).
#'   Unused categories are allowed. Leave `NULL` for interval/ratio scores.
#'
#' @details
#' The coefficient is \eqn{1-D_o/D_e}. Only units with at least two nonmissing
#' scores contribute to either disagreement. A unit with m ratings contributes
#' its ordered rating pairs divided by m-1 to the coincidence counts, hence
#' contributes m pairable ratings overall. Units are not equally weighted
#' when m differs. The expected disagreement draws two distinct observations
#' without replacement from the pooled pairable ratings, with denominator
#' n(n-1), where n is the number of pairable ratings. This is not Cohen's kappa
#' chance model or MFRM model-implied expected agreement.
#'
#' The squared difference function is 0/1 mismatch for nominal scores,
#' frequency-based rank distance for ordinal scores, squared numeric difference
#' for interval scores, and \eqn{((c-k)/(c+k))^2} for ratio scores (0 when both
#' values are zero). Ordinal positions use cumulative pairable-category counts
#' minus half the category count, divided by n. The squared difference between
#' these positions is proportional to Krippendorff's ordinal distance. It is
#' not quadratic weighting of category numbers or equally spaced category ranks.
#' Unused ranks and singleton-only categories do not alter pairable distances.
#'
#' Interval distances are calculated after subtracting `score_offset`, dividing
#' by `score_half_range`, and dividing by 2, using only pairable values to set
#' the offset and half-range. Constant pairable values map to zero. This avoids
#' squaring very large or tiny raw scores. The common normalization leaves alpha
#' unchanged; reported disagreements are on this normalized scale, not original
#' squared score units. Ordinal positions are similarly normalized by n.
#' Nominal and ratio disagreements use their stated distance directly.
#' Computation avoids a dense category-by-category matrix. Ratio computation
#' still requires quadratic time in the number of distinct pairable values.
#'
#' No rows are imputed and repeated rater-unit cells are rejected, including
#' repetitions with missing scores. All input rows remain in `cases`.
#' `unit_summary` distinguishes missing scores on assigned rows from absent
#' assignments relative to the observed rater roster. Units with no input rows
#' at all are unknown to the function. Missingness and selective assignment
#' may affect representativeness; numerical handling does not correct that bias.
#'
#' Negative estimates are retained. With no pairable ratings or zero expected
#' disagreement, alpha is unavailable with a reason rather than set to 0 or 1.
#' A finite coefficient describes this pooled rating sample, not each rater's
#' accuracy against a reference or a guarantee of interchangeable ratings.
#' There are no confidence intervals, tests, quality cutoffs or bootstrap
#' claims. Dependencies among several tasks from one person require attention
#' in any subsequent sampling inference. Circular, bipolar and custom metrics
#' are not implemented.
#'
#' @return A plain list with `table` (alpha, observed/expected disagreement,
#'   availability and sample counts), `counts` (category values, nonmissing
#'   and pairable counts, plus ordinal positions), `cases` (input rows and
#'   exclusion reasons), `unit_summary` (assignment counts and unit disagreements),
#'   `units` (unique input unit keys), `settings` and `notes`. `Unit` indexes
#'   rows of `units`. Save with [base::saveRDS()]; export data frames with
#'   [utils::write.csv()]. Reopening the saved list performs no calculation.
#'
#' @references
#' Krippendorff, K. (2011; literature updated 2013).
#' *Computing Krippendorff's Alpha-Reliability*. Sections A-E, pp. 1-10.
#' \url{https://www.asc.upenn.edu/krippendorffs-alpha-reliability}
#' @seealso [pairwise_icc_table()], [interrater_agreement_table()],
#'   [score_comparison_table()]
#' @examples
#' data("mfrmr_example_agreement", package = "mfrmr")
#' ratings <- mfrmr_example_agreement
#' agreement <- krippendorff_alpha_table(ratings, rater = "Rater",
#'   unit = c("Person", "Task"), score = "Score", level = "ordinal", categories = 0:3)
#' agreement$table
#' agreement$unit_summary # Each response needs at least two observed ratings.
#' agreement$notes
#' # Ordinal alpha uses category ordering; it is not quadratic weighted kappa.
#' @export
krippendorff_alpha_table <- function(data, rater, unit, score,
                                     level = c("nominal", "ordinal", "interval", "ratio"),
                                     categories = NULL) {
  level <- match.arg(level)
  roster <- mfrm_rating_roster(data, rater, unit, score)
  y <- data[[score]]
  if (is.factor(y)) y <- as.character(y)
  valid_values <- function(v) !is.object(v) && is.null(dim(v)) &&
    (is.character(v) || is.logical(v) || (is.numeric(v) && !is.complex(v))) &&
    (!is.numeric(v) || all(is.na(v) | is.finite(v))) &&
    (!is.character(v) || all(is.na(v) | nzchar(trimws(v))))
  if (!valid_values(y)) stop("Scores must be finite numbers or nonblank category labels, with NA for missing values.", call. = FALSE)
  numeric_level <- level %in% c("interval", "ratio")
  if (numeric_level && (!is.numeric(y) || (level == "ratio" && any(y < 0, na.rm = TRUE)))) stop(
    "Interval scores must be numeric; ratio scores must also be nonnegative.", call. = FALSE)
  if (numeric_level && !is.null(categories)) stop(
    "Leave `categories` NULL for interval or ratio scores.", call. = FALSE)
  if (level == "ordinal" && is.null(categories)) stop(
    "Supply `categories` in the substantive ordinal order.", call. = FALSE)
  if (is.factor(categories)) categories <- as.character(categories)
  if (!is.null(categories) && (!valid_values(categories) || !length(categories) ||
      anyNA(categories) || anyDuplicated(categories) ||
      is.numeric(categories) != is.numeric(y) || is.character(categories) != is.character(y))) stop(
    "`categories` must be distinct nonmissing values of the same representation as scores.", call. = FALSE)
  if (is.null(categories)) categories <- if (numeric_level) sort(unique(y[!is.na(y)])) else unique(y[!is.na(y)])
  code <- match(y, categories)
  if (any(!is.na(y) & is.na(code))) stop("A score is outside the supplied categories.", call. = FALSE)
  u <- roster$unit_id; nu <- nrow(roster$units); k <- length(categories)
  observed <- tabulate(u[!is.na(y)], nbins = nu)
  assigned <- tabulate(u, nbins = nu)
  included <- !is.na(y) & observed[u] >= 2L
  counts <- tabulate(code[included], nbins = k)
  n <- sum(counts)
  positions <- rep(NA_real_, k)
  offset <- half_range <- NA_real_
  if (level == "ordinal" && n > 0) positions <- (cumsum(counts) - counts / 2) / n
  if (level == "interval" && n > 0) {
    retained <- y[included]
    offset <- min(retained) / 2 + max(retained) / 2
    half_range <- max(abs(retained - offset))
    take <- counts > 0
    positions[take] <- if (half_range > 0) (categories[take] - offset) / half_range / 2 else 0
  }
  # Mean distance over distinct ordered observations; no dense category matrix.
  disagreement <- function(indices) {
    m <- length(indices)
    if (m < 2L) return(NA_real_)
    if (level %in% c("ordinal", "interval")) {
      z <- positions[indices]
      return(2 * sum((z - mean(z))^2) / (m - 1))
    }
    take <- sort(unique(indices))
    f <- tabulate(match(indices, take), nbins = length(take))
    if (level == "nominal") return((m^2 - sum(f^2)) / (m * (m - 1)))
    v <- categories[take]
    total <- 0
    if (length(v) > 1L) for (j in 2:length(v)) {
      before <- seq_len(j - 1L)
      ratio <- v[before] / v[j]
      total <- total + 2 * f[j] * sum(f[before] * ((1 - ratio) / (1 + ratio))^2)
    }
    total / (m * (m - 1))
  }
  unit_disagreement <- rep(NA_real_, nu)
  groups <- split(code[included], u[included])
  if (length(groups)) unit_disagreement[as.integer(names(groups))] <- vapply(groups, disagreement, numeric(1))
  do <- de <- alpha <- NA_real_
  status <- "no_pairable_ratings"
  if (n > 0) {
    pairable <- observed >= 2L
    do <- sum(observed[pairable] * unit_disagreement[pairable]) / n
    de <- disagreement(code[included])
    status <- if (!is.finite(do) || !is.finite(de)) "nonfinite_calculation" else
      if (de <= 0) "zero_expected_disagreement" else "available"
    if (status == "available") alpha <- 1 - do / de
    if (!is.finite(alpha) && status == "available") {
      alpha <- NA_real_; status <- "nonfinite_calculation"
    }
  }
  list(table = data.frame(Alpha = alpha, Level = level, Status = status,
      ObservedDisagreement = do, ExpectedDisagreement = de,
      Raters = length(roster$raters), Units = nu, PairableUnits = sum(observed >= 2L),
      PairableRatings = n, NonmissingRatings = sum(observed), MissingScores = sum(is.na(y)),
      SingletonUnits = sum(observed == 1L), UnratedUnits = sum(observed == 0L)),
    counts = data.frame(Category = categories, NonmissingRatings = tabulate(code, nbins = k),
      PairableRatings = counts, OrdinalPosition = if (level == "ordinal") positions else rep(NA_real_, k)),
    cases = data.frame(Row = seq_along(y), Rater = as.character(roster$ids[[rater]]),
      Unit = u, Score = y, MissingScore = is.na(y), Included = included,
      Reason = ifelse(is.na(y), "missing_score", ifelse(included, "included", "fewer_than_two_ratings"))),
    unit_summary = data.frame(Unit = seq_len(nu), AssignedRows = assigned,
      NonmissingRatings = observed, MissingScores = assigned - observed,
      AbsentAssignments = length(roster$raters) - assigned, Included = observed >= 2L,
      ObservedDisagreement = unit_disagreement), units = roster$units,
    settings = list(rater = rater, unit = unit, score = score, level = level,
      categories = categories, score_offset = offset, score_half_range = half_range,
      missing = "units_with_two_or_more_nonmissing_ratings", uncertainty = "not_computed"),
    notes = c(
      "Krippendorff's alpha measures pooled observed-score agreement, not Cronbach's alpha, MFRM separation reliability or criterion-referenced rater accuracy.",
      "Only pairable ratings contribute to observed and expected disagreement; unequal rating counts weight units differently.",
      "Missingness and selective assignment may affect representativeness; no imputation or automatic correction is made.",
      "No confidence intervals, tests or quality cutoffs are supplied; negative estimates are retained.",
      if (level == "ordinal") "Ordinal distances use pairable-category frequencies, not equally spaced ranks or QWK weights." else
        if (level == "interval") "Disagreements use a common centered/scaled score metric, not original squared score units." else
          paste0("Disagreements use the ", level, " distance function.")))
}
