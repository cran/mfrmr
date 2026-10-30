#' Compute observed-score ICCs for each rater pair
#'
#' Compare each pair on its common, complete units using a two-way ANOVA ICC.
#' No MFRM fit, mixed-model estimation or extra package is required.
#'
#' @param data Long-format data with at most one row per rater and unit.
#' @param rater Name of the rater identifier column.
#' @param unit Names of columns jointly identifying a rated unit, for example
#'   `c("Person", "Task", "Criterion")`. These columns match observations;
#'   they are not separately estimated variance components.
#' @param score Name of a numeric score column. Scores must share a meaningful
#'   interval scale. Numeric ordinal category codes do not establish that
#'   assumption. Factors/characters are not silently converted to scores.
#' @param model `"mixed"` (default) describes the selected raters. `"random"`
#'   assumes raters were sampled from a population to which results may
#'   generalize. Two-way point formulas are identical; their interpretation
#'   differs. Choosing `"random"` does not establish random sampling.
#' @param definition `"agreement"` (default) includes systematic rater offsets;
#'   `"consistency"` removes their contribution.
#' @param type `"single"` (default) concerns one rating; `"average"` concerns
#'   the unweighted mean of the two ratings, not a mean over the full panel.
#'
#' @details
#' All observed rater identifiers are paired. Each pair uses only units with
#' two nonmissing scores. Other raters' missing scores do not exclude a unit
#' for that pair. Missing-score rows and absent assignments are distinguished
#' in `cases`. Duplicate rater-unit cells (including missing-score duplicates)
#' are rejected rather than averaged. Identifiers must be nonmissing and
#' nonblank. There is no imputation or observation weighting. Pairs with fewer
#' than two complete units remain in `table` with an unavailable ICC.
#'
#' Let R, C and E be the unit, rater and residual ANOVA mean squares on a
#' complete n-by-2 matrix. Following Koo and Li (2016, Table 3), the numerator
#' is R-E. The single-rating denominator is R+E for consistency, or
#' R+E+2(C-E)/n for absolute agreement. The average-rating denominator is R
#' for consistency, or R+(C-E)/n for agreement. `Form` records ICC(C,1),
#' ICC(A,1), ICC(C,2) or ICC(A,2); interpret these together with `Model`.
#' ICC(A,2) denotes an average of two ratings, not Shrout--Fleiss ICC(2,1).
#' One-way ICCs are not computed.
#'
#' Negative estimates are retained. A nonpositive denominator makes `ICC`
#' unavailable and sets `Status = "nonpositive_denominator"`; `RawICC`
#' retains the algebraic ratio when finite so it cannot masquerade as a valid
#' reliability estimate. An entirely constant matrix is undefined, not 1.
#' Internally scores are shifted and scaled to avoid overflow/underflow.
#' `anova` records scaled mean squares, the offset/scale, numerator and
#' denominator, so the calculation can be reconstructed. These scaled mean
#' squares are not variances in the original score units.
#'
#' Results describe each pair's retained unit mix. Different pairs can have
#' different units and score ranges; their ICCs cannot automatically rank
#' rater quality and should not be averaged into a panel ICC. Restricted unit
#' variability can lower ICC without larger absolute disagreement. Inference
#' requires a suitable sampling design and assumptions about errors and
#' dependence. Treating several tasks from one person as separate units does
#' not make them independent. No confidence intervals, hypothesis tests,
#' quality bands or rater-training decisions are supplied. This is distinct
#' from the component variance shares in [compute_facet_icc()], MFRM separation
#' reliability and the fitted expected agreement in [interrater_agreement_table()].
#'
#' @return A plain list containing `table` (one row per rater pair), `anova`
#'   (scaled ANOVA ingredients), `cases` (one row per pair and union of assigned
#'   units, with input-row references and inclusion flags), `units` (the unique
#'   input unit keys in input-column order), `settings` and `notes`. `Unit` in
#'   `cases` is a row index into `units`. Units assigned to neither member of
#'   a pair are outside that pair's audit. `Missing1`/`Missing2` indicate
#'   explicit missing-score rows; `Assigned1`/`Assigned2` also distinguish
#'   absent assignments. Save with [base::saveRDS()] and export data frames
#'   with [utils::write.csv()]. Reading saved output performs no calculation.
#'
#' @references
#' Koo, T. K., & Li, M. Y. (2016). A guideline of selecting and reporting
#' intraclass correlation coefficients for reliability research.
#' *Journal of Chiropractic Medicine, 15*(2), 155-163.
#' \doi{10.1016/j.jcm.2016.02.012}. Table 3; see also the 2017 erratum,
#' \doi{10.1016/j.jcm.2017.10.001}, correcting the one-way ICC(1,1) formula.
#' @seealso [interrater_agreement_table()], [compute_facet_icc()],
#'   [krippendorff_alpha_table()], [score_comparison_table()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("pairwise_icc_table", run.dontrun = TRUE)
#' \dontrun{
#' data("mfrmr_example_agreement", package = "mfrmr")
#' ratings <- mfrmr_example_agreement
#' absolute <- pairwise_icc_table(ratings, rater = "Rater",
#'   unit = c("Person", "Task"), score = "Score", definition = "agreement")
#' consistent <- pairwise_icc_table(ratings, rater = "Rater",
#'   unit = c("Person", "Task"), score = "Score", definition = "consistency")
#' absolute$table
#' consistent$table
#' # R1/R2 have perfect consistency, but R2's one-category offset reduces agreement.
#' # These ICCs treat 0--3 as numeric scores; justify equal score intervals.
#' }
#' @export
pairwise_icc_table <- function(data, rater, unit, score,
                               model = c("mixed", "random"),
                               definition = c("agreement", "consistency"),
                               type = c("single", "average")) {
  model <- match.arg(model); definition <- match.arg(definition); type <- match.arg(type)
  roster <- mfrm_rating_roster(data, rater, unit, score)
  ids <- roster$ids; raters <- roster$raters
  unit_id <- roster$unit_id; units <- roster$units
  y <- data[[score]]
  if (!is.numeric(y) || is.complex(y) || is.object(y) || !is.null(dim(y)) ||
      any(!is.na(y) & !is.finite(y))) stop(
    "Scores must be finite numeric values or NA; category labels are not converted.", call. = FALSE)
  pairs <- utils::combn(raters, 2, simplify = FALSE)
  table_rows <- anova_rows <- case_rows <- vector("list", length(pairs))
  for (i in seq_along(pairs)) {
    p <- pairs[[i]]
    a <- which(as.character(ids[[rater]]) == p[1])
    b <- which(as.character(ids[[rater]]) == p[2])
    scope <- sort(unique(c(unit_id[a], unit_id[b])))
    a <- a[match(scope, unit_id[a])]; b <- b[match(scope, unit_id[b])]
    complete <- !is.na(a) & !is.na(b) & !is.na(y[a]) & !is.na(y[b])
    n <- sum(complete)
    case_rows[[i]] <- data.frame(Pair = i, Rater1 = p[1], Rater2 = p[2], Unit = scope,
      Row1 = a, Row2 = b, Score1 = y[a], Score2 = y[b],
      Assigned1 = !is.na(a), Assigned2 = !is.na(b),
      Missing1 = !is.na(a) & is.na(y[a]), Missing2 = !is.na(b) & is.na(y[b]),
      Included = complete)
    msr <- msc <- mse <- offset <- scale <- numerator <- denominator <- raw <- icc <- NA_real_
    status <- "fewer_than_two_complete_units"
    if (n >= 2L) {
      z <- cbind(y[a[complete]], y[b[complete]])
      offset <- min(z) / 2 + max(z) / 2
      z <- z - offset
      scale <- max(abs(z))
      if (scale > 0) z <- z / scale
      row_mean <- rowMeans(z); col_mean <- colMeans(z); grand <- mean(z)
      msr <- 2 * sum((row_mean - grand)^2) / (n - 1)
      msc <- n * sum((col_mean - grand)^2)
      residual <- sweep(sweep(z, 1, row_mean), 2, col_mean) + grand
      mse <- sum(residual^2) / (n - 1)
      numerator <- msr - mse
      denominator <- if (type == "single") msr + mse else msr
      if (definition == "agreement") denominator <- denominator +
        (if (type == "single") 2 else 1) * (msc - mse) / n
      if (denominator != 0) raw <- numerator / denominator
      status <- if (!is.finite(denominator) || !is.finite(numerator)) "nonfinite_calculation" else
        if (denominator <= 0) "nonpositive_denominator" else "available"
      if (status == "available" && is.finite(raw)) icc <- raw
      if (status == "available" && !is.finite(raw)) status <- "nonfinite_calculation"
      if (!is.finite(raw)) raw <- NA_real_
    }
    table_rows[[i]] <- data.frame(Pair = i, Rater1 = p[1], Rater2 = p[2],
      N = n, CandidateUnits = length(scope), ExcludedUnits = length(scope) - n,
      ICC = icc, RawICC = raw, Status = status, Model = paste0("two_way_", model),
      Definition = definition, Type = type,
      Form = paste0("ICC(", if (definition == "agreement") "A" else "C", ",",
        if (type == "single") "1" else "2", ")"))
    anova_rows[[i]] <- data.frame(Pair = i, Rater1 = p[1], Rater2 = p[2], N = n,
      ScaledMSUnit = msr, ScaledMSRater = msc, ScaledMSError = mse,
      ScoreOffset = offset, ScoreScale = scale, Numerator = numerator, Denominator = denominator,
      DFUnit = if (n >= 2L) n - 1L else NA_integer_,
      DFRater = if (n >= 2L) 1L else NA_integer_,
      DFError = if (n >= 2L) n - 1L else NA_integer_)
  }
  list(table = do.call(rbind, table_rows), anova = do.call(rbind, anova_rows),
    cases = do.call(rbind, case_rows), units = units,
    settings = list(rater = rater, unit = unit, score = score, model = model,
      definition = definition, type = type, raters_per_pair = 2L,
      missing = "pairwise_complete_units", uncertainty = "not_computed"),
    notes = c(
      "Two-way observed-score ICCs describe each pair's retained unit mix; they are not MFRM separation reliability or facet variance shares.",
      "Different pair cohorts and unit variability affect comparisons; do not average these values into a panel ICC.",
      "No intervals, tests, quality bands or independence/representativeness guarantees are supplied.",
      if (model == "mixed") "Results concern the selected raters." else
        "Generalization requires a justified random-rater sampling design; selecting this model does not establish it.",
      if (type == "average") "The target is the unweighted mean of these two ratings, not the whole panel." else
        "The target is a single rating, not the mean of two ratings."))
}
