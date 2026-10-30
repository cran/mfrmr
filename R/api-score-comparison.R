#' Compare aligned scores against reference ratings
#'
#' Summarize agreement, association and mean differences for two score sources,
#' with an optional second-reference baseline. Sources can be human raters,
#' scoring systems or other scoring procedures using the same rubric.
#' No measurement-model fit is needed for the score comparisons.
#'
#' @param scores Numeric scores from the source being evaluated, before rounding.
#' @param reference Numeric reference ratings on the declared categories.
#'   A reference is a comparison source, not necessarily an error-free criterion.
#' @param reference2 Optional second reference rating of each response.
#'   Its agreement with `reference` supplies the baseline for degradation.
#' @param categories All permitted numeric score categories, in increasing
#'   order, including categories not observed in this sample.
#' @param score_categories Scores expressed on `categories`. Defaults to
#'   `scores`. For continuous scores, supply an explicitly categorized vector.
#'   No rounding, clipping, category inference or averaging is performed.
#' @param smd SMD denominator: `"pooled"` uses the square root of the mean of
#'   the two sample variances; `"reference"` uses the reference sample SD.
#'   Both SMDs are returned. SDs use the sample divisor `N - 1`.
#' @param screen `"none"` (default) applies no automatic thresholds.
#'   `"williamson2012"` selects the automated-scoring examples discussed by
#'   Williamson, Xi, and Breyer (2012): QWK and correlation >= .70,
#'   agreement degradation <= .10, and absolute SMD <= .15.
#' @param qwk_min,correlation_min Optional minimum QWK and Pearson correlation.
#' @param max_degradation Optional maximum baseline-minus-score agreement,
#'   applied separately to QWK and correlation. Negative degradation indicates
#'   higher agreement of `scores` with `reference` than the baseline.
#' @param max_abs_smd Optional maximum absolute selected SMD.
#'   Non-NULL threshold arguments override the selected screen. With
#'   `screen = "none"`, only explicitly supplied thresholds are applied.
#' @param fit Optional native additive RSM/PCM [fit_mfrm()] result containing
#'   the selected scorer as a rater. This function does not fit a model.
#' @param diagnostics Optional matching [diagnose_mfrm()] result, supplying
#'   saved Infit/Outfit. No diagnostics are computed here.
#' @param rater_facet,rater With `fit`, the explicit facet name and level
#'   identifying the scorer being evaluated.
#'
#' @details
#' Vectors must already be aligned: one entry per same response in the same
#' order, on one common rubric/score scale. Match Person, task, criterion and
#' occasion before calling. Repeated ratings are not collapsed. All statistics
#' use the same complete cases across supplied vectors, including
#' `score_categories` and `reference2`. Every row's inclusion and missing-input
#' flags are retained in `cases`; missingness may affect representativeness.
#'
#' QWK is \eqn{1-D_o/D_e}, where \eqn{D_o} is observed mean squared category
#' disagreement and \eqn{D_e} uses independent draws from empirical marginals.
#' Distances use numeric category values, not their ranks. The denominator is
#' not MFRM model-implied expected agreement. Constant identical ratings have
#' undefined QWK. Correlation uses unrounded scores. SMD is signed `scores`
#' minus `reference`; numerical screens use its absolute value.
#'
#' `ReferenceQWK` and `ReferenceCorrelation` compare `reference2` with
#' `reference`; degradation subtracts the score-reference statistic from that
#' baseline. These columns are unavailable without a second reference.
#' `criteria` retains all five statistics, even when no thresholds are requested.
#' Its status is `"not_requested"` for a finite statistic without a threshold,
#' `"unavailable"` for an undefined statistic, or `"within_threshold"` /
#' `"review"` when an inclusive threshold is applied. No combined approval,
#' significance test or sampling uncertainty is inferred.
#'
#' The optional Williamson screen concerns joint human/automated scoring;
#' it is not a universal standard for all score comparisons. It evaluates
#' unrounded statistics and does not reproduce the paper's unspecified rounding
#' convention or operational exceptions. Its overall SMD example is .15; the
#' subgroup discussion uses .10. The source's p. 8 prose uses human-SD SMD,
#' whereas p. 12 note 4 uses the pooled formula. The default follows the formula
#' and retains `SMDReference` for sensitivity review. Evaluation data/tasks
#' should be separate from model development as required by the intended use.
#' Construct relevance, reference quality, external validity, generalizability,
#' subgroup differences and consequences require additional evidence.
#'
#' Optional MFRM columns describe **all retained source ratings for the selected
#' rater**, which need not be the comparison cohort. Source counts, readiness,
#' model settings and scope are retained separately. Severity is in logits,
#' oriented so positive values imply lower expected ratings; its zero is the
#' fitted model reference. Severity and Infit/Outfit do not measure criterion
#' accuracy and receive no automatic quality threshold.
#'
#' @return A list with `table` (one comparison row), `criteria` (five statistics,
#'   optional thresholds, availability and reasons), `cases` (row disposition),
#'   `counts` (`scores_reference` and `reference2_reference` category tables),
#'   `settings`, `mfrm` (optional saved source metadata), and `notes`.
#'   Save with [base::saveRDS()] and reopen with [base::readRDS()]. Export the
#'   data frames with [utils::write.csv()]. No ICC, alpha or AC2 is computed.
#'
#' @references
#' Williamson, D. M., Xi, X., & Breyer, F. J. (2012). A framework for
#' evaluation and use of automated scoring. *Educational Measurement: Issues
#' and Practice, 31*(1), 2-13. \doi{10.1111/j.1745-3992.2011.00223.x}
#' Especially pp. 7-10 and p. 12, notes 1 and 4.
#' @seealso [interrater_agreement_table()], [mfrm_report()],
#'   [automated_scoring_table()] for the earlier automated-scoring interface.
#' @examples
#' comparison <- score_comparison_table(
#'   scores = c(0, 1, 2, 2, 3, 3), reference = c(0, 1, 1, 2, 3, 3),
#'   categories = 0:3
#' )
#' comparison$table
#' # The automated-scoring screen is an explicit application choice.
#' auto_comparison <- score_comparison_table(
#'   scores = c(0.2, 0.9, 1.3, 2.1, 2.7, 2.8),
#'   reference = c(0, 1, 1, 2, 3, 3), reference2 = c(0, 1, 2, 2, 2, 3),
#'   categories = 0:3, score_categories = c(0, 1, 1, 2, 3, 3),
#'   screen = "williamson2012"
#' )
#' auto_comparison$criteria
#' @export
score_comparison_table <- function(scores, reference, reference2 = NULL,
                                   categories, score_categories = scores,
                                   smd = c("pooled", "reference"),
                                   screen = c("none", "williamson2012"),
                                   qwk_min = NULL, correlation_min = NULL,
                                   max_degradation = NULL, max_abs_smd = NULL, fit = NULL,
                                   diagnostics = NULL, rater_facet = NULL,
                                   rater = NULL) {
  smd <- match.arg(smd)
  screen <- match.arg(screen)
  vectors <- list(Score = scores, Reference = reference,
                  ScoreCategory = score_categories)
  if (!is.null(reference2)) vectors$Reference2 <- reference2
  for (nm in names(vectors)) {
    value <- vectors[[nm]]
    if (!is.numeric(value) || !is.null(dim(value)) ||
        any(!is.na(value) & !is.finite(value))) {
      stop(nm, " must be a numeric vector containing finite scores or NA.", call. = FALSE)
    }
  }
  if (length(unique(lengths(vectors))) != 1L) {
    stop("All rating vectors must have the same length and response order.", call. = FALSE)
  }
  if (!is.numeric(categories) || !is.null(dim(categories)) || length(categories) < 2L ||
      any(!is.finite(categories)) || any(diff(categories) <= 0) ||
      !is.finite(diff(range(categories)))) {
    stop("`categories` must contain at least two finite, strictly increasing score values.", call. = FALSE)
  }
  for (nm in intersect(c("Reference", "Reference2", "ScoreCategory"), names(vectors))) {
    value <- vectors[[nm]]
    if (any(!is.na(value) & !value %in% categories)) {
      stop(nm, " contains scores outside `categories`; supply category scores explicitly.", call. = FALSE)
    }
  }
  thresholds <- list(qwk_min = qwk_min, correlation_min = correlation_min,
                     max_degradation = max_degradation, max_abs_smd = max_abs_smd)
  if (screen == "williamson2012") {
    defaults <- list(qwk_min = .70, correlation_min = .70,
                     max_degradation = .10, max_abs_smd = .15)
    for (nm in names(thresholds)) {
      if (is.null(thresholds[[nm]])) thresholds[[nm]] <- defaults[[nm]]
    }
  }
  for (nm in names(thresholds)) {
    value <- thresholds[[nm]]
    if (is.null(value)) next
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        (nm %in% c("qwk_min", "correlation_min") && abs(value) > 1) ||
        (nm %in% c("max_degradation", "max_abs_smd") && value < 0)) {
      stop("Invalid numerical threshold: ", nm, ".", call. = FALSE)
    }
  }
  missing <- as.data.frame(lapply(vectors, is.na))
  included <- rowSums(missing) == 0L
  names(missing) <- paste0("Missing", names(missing))
  cases <- cbind(data.frame(Row = seq_along(scores), Included = included), missing)
  a <- scores[included]; h <- reference[included]
  ac <- score_categories[included]
  hr <- if (!is.null(reference2)) reference2[included] else NULL
  n <- length(a)
  ah <- mfrm_score_qwk(ac, h, categories)
  hh <- if (!is.null(hr)) mfrm_score_qwk(hr, h, categories) else
    list(value = NA_real_, reason = "Second reference ratings were not supplied.", counts = NULL)
  correlation <- function(x, y) {
    if (length(x) < 2L) return(NA_real_)
    sx <- stats::sd(x); sy <- stats::sd(y)
    if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(NA_real_)
    as.numeric(stats::cor(x, y))
  }
  r <- correlation(a, h)
  rh <- if (!is.null(hr)) correlation(hr, h) else NA_real_
  mean_a <- if (n) mean(a) else NA_real_
  mean_h <- if (n) mean(h) else NA_real_
  sd_a <- if (n > 1L) stats::sd(a) else NA_real_
  sd_h <- if (n > 1L) stats::sd(h) else NA_real_
  # Scaling first avoids squaring a very large but finite SD.
  scale <- max(sd_a, sd_h)
  pooled <- if (is.finite(scale) && scale > 0)
    scale * sqrt(((sd_a / scale)^2 + (sd_h / scale)^2) / 2) else NA_real_
  divide <- function(d) {
    value <- (mean_a - mean_h) / d
    if (is.finite(d) && d > 0 && is.finite(value)) value else NA_real_
  }
  smd_pooled <- divide(pooled); smd_reference <- divide(sd_h)
  smd_value <- if (smd == "pooled") smd_pooled else smd_reference
  tab <- data.frame(NInput = length(scores), N = n, NExcluded = sum(!included),
    Reference2Supplied = !is.null(reference2),
    MeanScore = mean_a, MeanReference = mean_h, SDScore = sd_a, SDReference = sd_h,
    QWK = ah$value, ReferenceQWK = hh$value, QWKDegradation = hh$value - ah$value,
    Correlation = r, ReferenceCorrelation = rh, CorrelationDegradation = rh - r,
    SMD = smd_value, SMDBasis = smd, SMDPooled = smd_pooled, SMDReference = smd_reference,
    MFRMSeverity = NA_real_, MFRMInfit = NA_real_, MFRMOutfit = NA_real_,
    MFRMRatingRows = NA_integer_, MFRMWeightSum = NA_real_,
    MFRMFacet = NA_character_, MFRMRater = NA_character_,
    MFRMModel = NA_character_, MFRMMethod = NA_character_,
    MFRMStatus = "not_supplied", MFRMScope = "not_supplied")
  cor_reason <- "Correlation requires at least two paired scores and positive finite SDs."
  qwk_gap_reason <- paste(c(ah$reason, hh$reason)[nzchar(c(ah$reason, hh$reason))], collapse = " ")
  criteria <- data.frame(
    Metric = c("QWK", "Correlation", "QWKDegradation", "CorrelationDegradation", "AbsoluteSMD"),
    Value = c(ah$value, r, tab$QWKDegradation, tab$CorrelationDegradation, abs(smd_value)),
    Operator = c(">=", ">=", "<=", "<=", "<="),
    Threshold = unname(vapply(thresholds, function(x) x %||% NA_real_, numeric(1))[
      c("qwk_min", "correlation_min", "max_degradation", "max_degradation", "max_abs_smd")]),
    Reason = c(ah$reason, cor_reason, qwk_gap_reason,
      if (is.null(hr)) hh$reason else cor_reason,
      "SMD requires at least two paired scores and a positive finite selected denominator."))
  criteria$Available <- is.finite(criteria$Value)
  criteria$N <- n
  criteria$Basis <- c("categorized_scores", "unrounded_scores",
    "categorized_scores", "unrounded_scores", paste0("unrounded_scores_smd_", smd))
  within <- ifelse(criteria$Operator == ">=", criteria$Value >= criteria$Threshold,
                    criteria$Value <= criteria$Threshold)
  criteria$Status <- ifelse(!criteria$Available, "unavailable",
    ifelse(is.na(criteria$Threshold), "not_requested",
           ifelse(within, "within_threshold", "review")))
  criteria$Reason[criteria$Available] <- ""
  criteria$Reason[criteria$Status == "not_requested"] <- "No threshold requested."
  mfrm <- NULL
  if (!is.null(fit)) {
    if (!inherits(fit, "mfrm_fit") || !isTRUE(fit$config$model %in% c("RSM", "PCM"))) {
      stop("Optional MFRM columns require a native additive RSM/PCM fit.", call. = FALSE)
    }
    # Reuse the source, direction, diagnostic identity and eligibility contract
    # of individual-rater reporting; do not invoke fitting or diagnostics.
    sheet <- mfrm_rater_feedback(list(fit = fit, diagnostics = diagnostics),
      facet = rater_facet, rater = rater, audience = "researcher",
      label = "Selected scorer", max_cases = 0L, interval = NULL)
    tab$MFRMSeverity <- sheet$tables$severity$Severity
    saved_fit <- sheet$tables$fit
    if (nrow(saved_fit)) {
      tab$MFRMInfit <- saved_fit$Value[match("Infit", saved_fit$Statistic)]
      tab$MFRMOutfit <- saved_fit$Value[match("Outfit", saved_fit$Statistic)]
    }
    tab$MFRMRatingRows <- sheet$tables$exposure$RatingRows
    tab$MFRMWeightSum <- sheet$tables$exposure$WeightSum
    tab$MFRMFacet <- rater_facet
    tab$MFRMRater <- rater
    tab$MFRMModel <- fit$config$model
    tab$MFRMMethod <- fit$config$method
    tab$MFRMStatus <- if (isTRUE(.mfrm_fit_plot_readiness(fit)$ready)) "descriptive" else "source_needs_review"
    tab$MFRMScope <- "all_retained_source_ratings_for_selected_rater"
    mfrm <- list(facet = rater_facet, rater = rater,
      config = fit$config, review = sheet$review,
      notes = sheet$notes[c("severity", "fit", "exposure")])
  } else if (!is.null(diagnostics) || !is.null(rater_facet) || !is.null(rater)) {
    stop("MFRM diagnostics and rater identifiers require `fit`.", call. = FALSE)
  }
  list(table = tab, criteria = criteria, cases = cases,
    counts = list(scores_reference = ah$counts, reference2_reference = hh$counts),
    settings = list(categories = categories, smd = smd, thresholds = thresholds,
      screen = if (screen == "none" && any(!is.na(criteria$Threshold))) "custom" else screen,
      threshold_rounding = "none", missing = "common_complete_cases",
      qwk_distance = "squared_numeric_category_distance",
      reference = "Williamson, Xi & Breyer (2012), pp. 7-10; p. 12 notes 1 and 4",
      doi = "10.1111/j.1745-3992.2011.00223.x"),
    mfrm = mfrm,
    notes = c(
      "Descriptive point comparisons, without sampling uncertainty or a combined approval decision. Evaluation data and tasks must match the intended use and be separate from model development as appropriate.",
      "All score statistics use common complete cases. Score vectors must already align the same responses; repeated observations are not averaged.",
      "QWK uses categorized scores; correlation and SMD use the supplied unrounded scores. Requested thresholds evaluate unrounded statistics with inclusive limits.",
      "SMD is scores minus reference. The pooled default follows note 4; the reference-SD alternative follows the discrepant prose on p. 8. Both values are retained.",
      "No universal thresholds apply to score comparisons. The optional Williamson screen concerns joint human/automated scoring. Construct relevance, reference quality, external validity, generalizability, subgroups and consequences require separate evidence.",
      "Optional MFRM values describe the selected rater's whole fitted cohort, whose identity with the comparison sample is not established. Severity and fit are distinct from agreement and accuracy."))
}

mfrm_score_qwk <- function(x, y, categories) {
  # Match numeric values before formatting labels, which otherwise can collapse
  # nearby distinct numeric categories into one factor level.
  counts <- table(factor(match(x, categories), levels = seq_along(categories)),
                  factor(match(y, categories), levels = seq_along(categories)),
                  dnn = c("Scorer", "Reference"))
  dimnames(counts) <- list(Scorer = format(categories, digits = 17, trim = TRUE),
                           Reference = format(categories, digits = 17, trim = TRUE))
  n <- length(x)
  value <- NA_real_
  reason <- "At least two complete responses are required."
  if (n >= 2L) {
    # Normalize distances to avoid unnecessary overflow. Scaling cancels in
    # the observed/expected disagreement ratio.
    used <- rowSums(counts) + colSums(counts) > 0
    active_counts <- counts[used, used, drop = FALSE]
    points <- categories[used]
    if (length(points) == 1L) return(list(value = NA_real_,
      reason = "Chance disagreement is zero; QWK is undefined for identical constant ratings.",
      counts = counts))
    distance <- mfrm_agreement_distance(as.numeric(points), "quadratic")
    observed <- sum(distance * (active_counts / n))
    expected <- sum(distance * outer(rowSums(active_counts) / n, colSums(active_counts) / n))
    if (is.finite(expected) && expected > 0) {
      value <- 1 - observed / expected
      reason <- ""
    } else reason <- "Chance disagreement is zero; QWK is undefined for identical constant ratings."
  }
  list(value = value, reason = reason, counts = counts)
}
