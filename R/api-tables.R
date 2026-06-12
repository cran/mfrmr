#' Build an inter-rater agreement report
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param rater_facet Name of the rater facet. If `NULL`, inferred from facet names.
#' @param context_facets Optional context facets used to match observations for
#'   agreement. If `NULL`, all remaining facets (including `Person`) are used.
#' @param exact_warn Warning threshold for exact agreement.
#' @param corr_warn Warning threshold for pairwise correlation.
#' @param include_precision If `TRUE`, append rater severity spread indices from
#'   the facet precision summary when available.
#' @param top_n Optional maximum number of pair rows to keep.
#'
#' @details
#' This helper computes pairwise rater agreement on matched contexts
#' and returns both a pair-level table and a one-row summary. The output is
#' package-native and does not require knowledge of legacy report numbering.
#'
#' @section Interpreting output:
#' - `summary`: overall agreement level, number/share of flagged pairs.
#' - `pairs`: pairwise exact agreement, correlation, and direction/size gaps.
#' - `settings`: applied facet matching and warning thresholds.
#'
#' Pairs flagged by both low exact agreement and low correlation generally
#' deserve highest calibration priority.
#'
#' @section Typical workflow:
#' 1. Run with explicit `rater_facet` (and `context_facets` if needed).
#' 2. Review `summary(ir)` and top flagged rows in `ir$pairs`.
#' 3. Visualize with [plot_interrater_agreement()].
#'
#' @section Output columns:
#' The `pairs` data.frame contains:
#' \describe{
#'   \item{Rater1, Rater2}{Rater pair identifiers.}
#'   \item{N}{Number of matched-context observations for this pair.}
#'   \item{Exact}{Proportion of exact score agreements.}
#'   \item{ExpectedExact}{Expected exact agreement under chance.}
#'   \item{Adjacent}{Proportion of adjacent (+/- 1 category) agreements.}
#'   \item{MeanDiff}{Signed mean score difference (Rater1 - Rater2).}
#'   \item{MAD}{Mean absolute score difference.}
#'   \item{Corr}{Pearson correlation between paired scores.}
#'   \item{Flag}{Logical; `TRUE` when Exact < `exact_warn` or Corr < `corr_warn`.}
#'   \item{OpportunityCount, ExactCount, ExpectedExactCount, AdjacentCount}{Raw
#'     counts behind the agreement proportions.}
#' }
#'
#' The `summary` data.frame contains:
#' \describe{
#'   \item{RaterFacet}{Name of the rater facet analyzed.}
#'   \item{TotalPairs}{Number of rater pairs evaluated.}
#'   \item{ExactAgreement}{Mean exact agreement across all pairs.}
#'   \item{AgreementMinusExpected}{Observed exact agreement minus expected exact
#'     agreement.}
#'   \item{MeanCorr}{Mean pairwise correlation.}
#'   \item{FlaggedPairs, FlaggedShare}{Count and proportion of flagged pairs.}
#'   \item{RaterSeparation, RaterReliability}{Severity-spread indices for the
#'     rater facet, reported separately from agreement.}
#' }
#'
#' @return A named list with:
#' - `summary`: one-row inter-rater summary
#' - `pairs`: pair-level agreement table
#' - `settings`: applied options and thresholds
#'
#' @seealso [diagnose_mfrm()], [facets_chisq_table()], [plot_interrater_agreement()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' ir <- interrater_agreement_table(fit, rater_facet = "Rater")
#' # One-row overview: ExactAgreement, ExpectedExactAgreement, MeanCorr,
#' # RaterSeparation, and RaterReliability are the headline reportable
#' # statistics.
#' ir$summary
#' # Per-pair detail (Rater1 vs Rater2 with Exact, Adjacent, Corr, MAD).
#' head(ir$pairs)
#' p_ir <- plot(ir, draw = FALSE)
#' p_ir$data$plot
#' @export
interrater_agreement_table <- function(fit,
                                       diagnostics = NULL,
                                       rater_facet = NULL,
                                       context_facets = NULL,
                                       exact_warn = 0.50,
                                       corr_warn = 0.30,
                                       include_precision = TRUE,
                                       top_n = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  known_facets <- c("Person", fit$config$facet_names)
  if (is.null(rater_facet) || !nzchar(as.character(rater_facet[1]))) {
    if (!is.null(diagnostics$interrater$summary) &&
        nrow(diagnostics$interrater$summary) > 0 &&
        "RaterFacet" %in% names(diagnostics$interrater$summary)) {
      rater_facet <- as.character(diagnostics$interrater$summary$RaterFacet[1])
    } else {
      rater_facet <- infer_default_rater_facet(fit$config$facet_names)
    }
  } else {
    rater_facet <- as.character(rater_facet[1])
  }
  if (is.null(rater_facet) || !rater_facet %in% known_facets) {
    stop("`rater_facet` must match one of: ", paste(known_facets, collapse = ", "))
  }
  if (identical(rater_facet, "Person")) {
    stop("`rater_facet = 'Person'` is not supported. Use a non-person facet.")
  }

  if (is.null(context_facets)) {
    facet_cols <- known_facets
  } else {
    context_facets <- unique(as.character(context_facets))
    unknown <- setdiff(context_facets, known_facets)
    if (length(unknown) > 0) {
      stop("Unknown `context_facets`: ", paste(unknown, collapse = ", "))
    }
    context_facets <- setdiff(context_facets, rater_facet)
    if (length(context_facets) == 0) {
      stop("`context_facets` must include at least one facet different from `rater_facet`.")
    }
    facet_cols <- c(rater_facet, context_facets)
  }

  agreement <- calc_interrater_agreement(
    obs_df = diagnostics$obs,
    facet_cols = facet_cols,
    rater_facet = rater_facet,
    res = fit
  )

  pairs <- as.data.frame(agreement$pairs, stringsAsFactors = FALSE)
  flagged_n <- 0L
  if (nrow(pairs) > 0) {
    pairs <- pairs |>
      mutate(
        Pair = paste(.data$Rater1, .data$Rater2, sep = " | "),
        ExactGap = ifelse(is.finite(.data$ExpectedExact), .data$Exact - .data$ExpectedExact, NA_real_),
        LowExactFlag = is.finite(.data$Exact) & .data$Exact < exact_warn,
        LowCorrFlag = is.finite(.data$Corr) & .data$Corr < corr_warn,
        Flag = .data$LowExactFlag | .data$LowCorrFlag
      ) |>
      arrange(desc(.data$Flag), .data$Exact, .data$Corr)
    flagged_n <- sum(pairs$Flag, na.rm = TRUE)
    if (!is.null(top_n)) {
      pairs <- pairs |>
        slice_head(n = max(1L, as.integer(top_n)))
    }
  }

  agreement <- if (isTRUE(include_precision)) {
    augment_interrater_with_precision(
      agreement,
      reliability_tbl = diagnostics$reliability,
      rater_facet = rater_facet
    )
  } else {
    agreement
  }

  summary_tbl <- as.data.frame(agreement$summary, stringsAsFactors = FALSE)
  if (nrow(summary_tbl) > 0) {
    summary_tbl$FlaggedPairs <- flagged_n
    summary_tbl$FlaggedShare <- ifelse(summary_tbl$Pairs > 0, flagged_n / summary_tbl$Pairs, NA_real_)
  }

  out <- list(
    summary = summary_tbl,
    pairs = pairs,
    settings = list(
      rater_facet = rater_facet,
      context_facets = setdiff(facet_cols, rater_facet),
      exact_warn = exact_warn,
      corr_warn = corr_warn,
      include_precision = include_precision
    )
  )
  as_mfrm_bundle(out, "mfrm_interrater")
}

#' Build facet variability diagnostics with fixed/random reference tests
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param fixed_p_max Warning cutoff for fixed-effect chi-square p-values.
#' @param random_p_max Warning cutoff for random-effect chi-square p-values.
#' @param top_n Optional maximum number of facet rows to keep.
#'
#' @details
#' This helper summarizes facet-level variability with fixed and random
#' chi-square indices for spread and heterogeneity checks.
#'
#' @section Interpreting output:
#' - `table`: facet-level fixed/random chi-square and p-value flags.
#' - `summary`: number of significant facets and overall magnitude indicators.
#' - `thresholds`: p-value criteria used for flagging.
#'
#' Use this table together with inter-rater and displacement diagnostics to
#' distinguish global facet effects from local anomalies.
#'
#' @section Typical workflow:
#' 1. Run `facets_chisq_table(fit, ...)`.
#' 2. Inspect `summary(chi)` then facet rows in `chi$table`.
#' 3. Visualize with [plot_facets_chisq()].
#'
#' @section Output columns:
#' The `table` data.frame contains:
#' \describe{
#'   \item{Facet}{Facet name.}
#'   \item{Levels}{Number of estimated levels in this facet.}
#'   \item{MeanMeasure, SD}{Mean and standard deviation of level measures.}
#'   \item{FixedChiSq, FixedDF, FixedProb}{Fixed-effect chi-square test
#'     (null hypothesis: all levels equal). Significant result means the
#'     facet elements differ more than measurement error alone.}
#'   \item{RandomChiSq, RandomDF, RandomProb, RandomVar}{Random-effect test
#'     (null hypothesis: variation equals that of a random sample from a
#'     single population). Significant result suggests systematic
#'     heterogeneity beyond sampling variation.}
#'   \item{FixedFlag, RandomFlag}{Logical flags for significance.}
#' }
#'
#' @return A named list with:
#' - `table`: facet-level chi-square diagnostics
#' - `summary`: one-row summary
#' - `thresholds`: applied p-value thresholds
#'
#' @seealso [diagnose_mfrm()], [interrater_agreement_table()], [plot_facets_chisq()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' chi <- facets_chisq_table(fit)
#' summary(chi)
#' p_chi <- plot(chi, draw = FALSE)
#' p_chi$data$plot
#' @export
facets_chisq_table <- function(fit,
                               diagnostics = NULL,
                               fixed_p_max = 0.05,
                               random_p_max = 0.05,
                               top_n = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$measures) || nrow(diagnostics$measures) == 0) {
    stop("`diagnostics$measures` is empty. Run diagnose_mfrm() first.")
  }

  tbl <- calc_facets_chisq(diagnostics$measures)
  if (nrow(tbl) > 0) {
    tbl <- tbl |>
      mutate(
        FixedFlag = is.finite(.data$FixedProb) & .data$FixedProb < fixed_p_max,
        RandomFlag = is.finite(.data$RandomProb) & .data$RandomProb < random_p_max
      ) |>
      arrange(desc(.data$FixedChiSq))
    if (!is.null(top_n)) {
      tbl <- tbl |>
        slice_head(n = max(1L, as.integer(top_n)))
    }
  }

  summary_tbl <- if (nrow(tbl) == 0) {
    data.frame()
  } else {
    safe_max <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else max(x)
    }
    safe_mean <- function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else mean(x)
    }
    data.frame(
      Facets = nrow(tbl),
      FixedSignificant = sum(tbl$FixedFlag, na.rm = TRUE),
      RandomSignificant = sum(tbl$RandomFlag, na.rm = TRUE),
      MeanRandomVar = safe_mean(tbl$RandomVar),
      MaxFixedChiSq = safe_max(tbl$FixedChiSq),
      MaxRandomChiSq = safe_max(tbl$RandomChiSq),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    table = as.data.frame(tbl, stringsAsFactors = FALSE),
    summary = summary_tbl,
    thresholds = list(
      fixed_p_max = fixed_p_max,
      random_p_max = random_p_max
    )
  )
  as_mfrm_bundle(out, "mfrm_facets_chisq")
}

#' Build an unexpected-response screening report
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param abs_z_min Absolute standardized-residual cutoff.
#' @param prob_max Maximum observed-category probability cutoff.
#' @param top_n Maximum number of rows to return.
#' @param rule Flagging rule: `"either"` (default) or `"both"`.
#'
#' @details
#' A response is flagged as unexpected when:
#' - `rule = "either"`: `|StdResidual| >= abs_z_min` OR `ObsProb <= prob_max`
#' - `rule = "both"`: both conditions must be met.
#'
#' The table includes row-level observed/expected values, residuals,
#' observed-category probability, most-likely category, and a composite
#' severity score for sorting.
#'
#' @section Interpreting output:
#' - `summary`: prevalence of unexpected responses under current thresholds.
#' - `table`: ranked row-level diagnostics for case review.
#' - `thresholds`: active cutoffs and flagging rule.
#'
#' Compare results across `rule = "either"` and `rule = "both"` to assess how
#' conservative your screening should be.
#'
#' @section Typical workflow:
#' 1. Start with `rule = "either"` for broad screening.
#' 2. Re-run with `rule = "both"` for strict subset.
#' 3. Inspect top rows and visualize with [plot_unexpected()].
#'
#' @section Further guidance:
#' For a plot-selection guide and a longer walkthrough, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @section Output columns:
#' The `table` data.frame contains:
#' \describe{
#'   \item{Row}{Original row index in the prepared data.}
#'   \item{Person}{Person identifier (plus one column per facet).}
#'   \item{Score}{Observed score category.}
#'   \item{Observed, Expected}{Observed and model-expected score values.}
#'   \item{Residual, StdResidual}{Raw and standardized residuals.}
#'   \item{ObsProb}{Probability of the observed category under the model.}
#'   \item{MostLikely, MostLikelyProb}{Most probable category and its
#'     probability.}
#'   \item{Severity}{Composite severity index (higher = more unexpected).}
#'   \item{Direction}{"Higher than expected" or "Lower than expected".}
#'   \item{FlagLowProbability, FlagLargeResidual}{Logical flags for each
#'     criterion.}
#' }
#'
#' The `summary` data.frame contains:
#' \describe{
#'   \item{TotalObservations}{Total observations analyzed.}
#'   \item{UnexpectedN, UnexpectedPercent}{Count and share of flagged rows.}
#'   \item{AbsZThreshold, ProbThreshold}{Applied cutoff values.}
#'   \item{Rule}{"either" or "both".}
#' }
#'
#' @return A named list with:
#' - `table`: flagged response rows
#' - `summary`: one-row overview
#' - `thresholds`: applied thresholds
#'
#' @seealso [diagnose_mfrm()], [displacement_table()], [fair_average_table()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy_full <- load_mfrmr_data("example_core")
#' toy_people <- unique(toy_full$Person)[1:12]
#' toy <- toy_full[toy_full$Person %in% toy_people, , drop = FALSE]
#' fit <- suppressWarnings(
#'   fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' )
#' t4 <- unexpected_response_table(fit, abs_z_min = 1.5, prob_max = 0.4, top_n = 5)
#' summary(t4)
#' p_t4 <- plot(t4, draw = FALSE)
#' p_t4$data$plot
#' @export
unexpected_response_table <- function(fit,
                                      diagnostics = NULL,
                                      abs_z_min = 2,
                                      prob_max = 0.30,
                                      top_n = 100,
                                      rule = c("either", "both")) {
  rule <- match.arg(tolower(rule), c("either", "both"))
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  tbl <- calc_unexpected_response_table(
    obs_df = diagnostics$obs,
    probs = compute_prob_matrix(fit),
    facet_names = fit$config$facet_names,
    rating_min = fit$prep$rating_min,
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    top_n = top_n,
    rule = rule
  )
  summary_tbl <- summarize_unexpected_response_table(
    unexpected_tbl = tbl,
    total_observations = nrow(diagnostics$obs),
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    rule = rule
  )

  out <- list(
    table = tbl,
    summary = summary_tbl,
    thresholds = list(
      abs_z_min = abs_z_min,
      prob_max = prob_max,
      rule = rule
    )
  )
  as_mfrm_bundle(out, "mfrm_unexpected")
}

#' Build an adjusted-score reference table bundle
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param facets Optional subset of facets.
#' @param totalscore Include all observations for score totals (`TRUE`) or apply
#'   legacy extreme-row exclusion (`FALSE`).
#' @param umean Additive score-to-report origin shift.
#' @param uscale Multiplicative score-to-report scale.
#' @param udecimals Rounding digits used in formatted output.
#' @param reference Which adjusted-score reference to keep in formatted outputs:
#'   `"both"` (default), `"mean"`, or `"zero"`.
#' @param label_style Column-label style for formatted outputs:
#'   `"both"` (default), `"native"`, or `"legacy"`.
#' @param omit_unobserved If `TRUE`, remove unobserved levels.
#' @param xtreme Extreme-score adjustment amount.
#' @param fair_se Logical. When `TRUE` and `fit` is an MML bounded-`GPCM`
#'   fit, add structural delta-method standard errors and confidence limits
#'   for `Fair(M)` / `AdjustedAverage` and `Fair(Z)` /
#'   `StandardizedAdjustedAverage`. Person rows remain `NA` because MML
#'   person EAP estimates are not part of the structural Hessian. For
#'   `RSM`, `PCM`, and `JML` fits this option leaves fair-average SE columns
#'   unavailable.
#' @param ci_level Confidence level used when `fair_se = TRUE`; default
#'   `0.95`.
#'
#' @details
#' This function wraps the package's adjusted-score calculations and returns
#' both facet-wise and stacked tables. Historical display columns such as
#' `Fair(M) Average` and `Fair(Z) Average` are retained for compatibility, and
#' package-native aliases such as `AdjustedAverage`,
#' `StandardizedAdjustedAverage`, `ModelBasedSE`, and `FitAdjustedSE` are
#' appended to the formatted outputs.
#'
#' For the Rasch-family `RSM` / `PCM` branch, these tables follow the
#' standard FACETS Linacre construction: fair averages are
#' Rasch-measure-to-score transformations evaluated in a standardized
#' mean/zero-facet environment.
#'
#' Bounded `GPCM` fits are supported under a slope-aware
#' element-conditional construction. For each slope-facet element
#' \eqn{j^\star} the per-row fair-average is the GPCM expected score
#' \deqn{\mathrm{FA}_{p, j^\star} = \sum_k k \cdot P_{GPCM}(X = k \mid \theta_p, a_{j^\star}, \boldsymbol{\delta}_{j^\star})}
#' computed at that element's own discrimination \eqn{a_{j^\star}}
#' and threshold structure. Rows for non-slope facets (Person, Rater,
#' \ldots) use the geometric-mean-one slope by the GPCM
#' identification convention, so those rows remain continuous with
#' the standard PCM Linacre fair-average and reduce to it exactly
#' when all slopes equal one.
#' This is an identification-based reporting convention for the package's
#' bounded `GPCM` route, not a unique free-discrimination score-side analogue
#' to FACETS fair averages. Do not report it as FACETS score-side equivalence
#' or as an operational scoring rule unless that convention is substantively
#' justified.
#'
#' Standard errors on the fair-average value itself are opt-in for MML
#' bounded `GPCM` fits via `fair_se = TRUE`. The original `SE`,
#' `Model S.E.`, `ModelBasedSE`, `Real S.E.`, and `FitAdjustedSE` columns
#' retain the same meaning as for PCM (scaled facet-measure SEs); fair-average
#' uncertainty is reported under distinct columns such as `Fair(M) S.E.`,
#' `Fair(M) CI Lower`, and `AdjustedAverageSE`.
#'
#' @section Interpreting output:
#' - `stacked`: cross-facet table for global comparison.
#' - `by_facet`: per-facet formatted tables for reporting.
#' - `raw_by_facet`: unformatted values for custom analyses/plots.
#' - `settings`: scoring-transformation and filtering options used.
#'
#' Larger observed-vs-fair gaps can indicate systematic scoring tendencies by
#' specific facet levels.
#'
#' @section Typical workflow:
#' 1. Run `fair_average_table(fit, ...)`.
#' 2. Inspect `summary(t12)` and `t12$stacked`.
#' 3. Visualize with [plot_fair_average()].
#'
#' @section Output columns:
#' The `stacked` data.frame contains:
#' \describe{
#'   \item{Facet}{Facet name for this row.}
#'   \item{Level}{Element label within the facet.}
#'   \item{Obsvd Average}{Observed raw-score average.}
#'   \item{Fair(M) Average}{Model-adjusted reference average on the reported score scale.}
#'   \item{Fair(Z) Average}{Standardized adjusted reference average.}
#'   \item{ObservedAverage, AdjustedAverage, StandardizedAdjustedAverage}{Package-native aliases for the three average columns above.}
#'   \item{AdjustedAverageSE, AdjustedAverageCI_Lower, AdjustedAverageCI_Upper}{Optional structural delta-method uncertainty for `AdjustedAverage` when `fair_se = TRUE` and available.}
#'   \item{StandardizedAdjustedAverageSE, StandardizedAdjustedAverageCI_Lower, StandardizedAdjustedAverageCI_Upper}{Optional structural delta-method uncertainty for `StandardizedAdjustedAverage` when `fair_se = TRUE` and available.}
#'   \item{Measure}{Estimated logit measure for this level.}
#'   \item{SE}{Compatibility alias for the model-based standard error.}
#'   \item{ModelBasedSE, FitAdjustedSE}{Package-native aliases for `Model S.E.` and `Real S.E.`.}
#'   \item{Infit MnSq, Outfit MnSq}{Fit statistics for this level.}
#' }
#'
#' @section Standard-error caveat (read before quoting CIs):
#' The `SE`, `Model S.E.`, `ModelBasedSE`, `Real S.E.`, and `FitAdjustedSE`
#' columns in this table are the **measure-level** standard errors of the
#' underlying facet element (the same SE that would appear in
#' `summary(fit)$facets`), rescaled by the fair-average score scale factor
#' so the units line up with the reported `Fair(M) Average` / `Fair(Z) Average`
#' columns. They are **not** delta-method standard errors of the fair-average
#' values themselves. When `fair_se = TRUE`, the distinct `Fair(M) S.E.` /
#' `Fair(Z) S.E.` columns are computed by
#' propagating the joint covariance of the relevant facet element, the
#' threshold parameters, and the slope parameters through the gradient of
#' \eqn{\mathrm{E}[X \mid \theta_p, j^\star]}. This is a structural
#' covariance calculation: MML person EAP estimates are conditioned on rather
#' than included in the Hessian, so person rows receive unavailable fair-average
#' SEs. **Do not use the measure-level `SE` / `Model S.E.` columns as
#' \eqn{\pm 1.96 \cdot \mathrm{SE}} confidence-interval bounds on the
#' fair-average value.**
#'
#' @return A named list with:
#' - `by_facet`: named list of formatted data.frames
#' - `stacked`: one stacked data.frame across facets
#' - `raw_by_facet`: unformatted component tables
#' - `settings`: resolved options
#'
#' @seealso [diagnose_mfrm()], [unexpected_response_table()], [displacement_table()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' t12 <- fair_average_table(fit, udecimals = 2)
#' t12_native <- fair_average_table(fit, reference = "mean", label_style = "native")
#' summary(t12)
#' p_t12 <- plot(t12, draw = FALSE)
#' p_t12$data$plot
#' }
#'
#' @section References:
#' - Linacre, J. M. (1989). *Many-Facet Rasch Measurement*. MESA Press.
#' - Linacre, J. M. (1994). *Many-facet Rasch Measurement* (2nd ed.).
#'   MESA Press.
#' - Linacre, J. M. (2026). *A user's guide to FACETS, version 4.5.0*.
#'   Winsteps.com. <https://www.winsteps.com/facets.htm>
#'   (FACETS Table 12 corresponds to the fair-average
#'   construction implemented here for `RSM` / `PCM` fits; the
#'   slope-aware element-conditional construction for bounded `GPCM`
#'   is documented in this help page.)
#' - Andrich, D. (1978). A rating formulation for ordered response
#'   categories. *Psychometrika, 43*(4), 561-573.
#'   \doi{10.1007/BF02293814}
#' - Masters, G. N. (1982). A Rasch model for partial credit scoring.
#'   *Psychometrika, 47*(2), 149-174. \doi{10.1007/BF02296272}
#' - Muraki, E. (1992). A generalized partial credit model:
#'   Application of an EM algorithm. *Applied Psychological
#'   Measurement, 16*(2), 159-176. (Cited for the bounded `GPCM`
#'   slope-aware extension.)
#' @export
fair_average_table <- function(fit,
                               diagnostics = NULL,
                               facets = NULL,
                               totalscore = TRUE,
                               umean = 0,
                               uscale = 1,
                               udecimals = 2,
                               reference = c("both", "mean", "zero"),
                               label_style = c("both", "native", "legacy"),
                               omit_unobserved = FALSE,
                               xtreme = 0,
                               fair_se = FALSE,
                               ci_level = 0.95) {
  reference <- match.arg(reference)
  label_style <- match.arg(label_style)
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  fit_model <- as.character(fit$config$model %||% fit$summary$Model[1] %||% NA_character_)
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || is.null(diagnostics$measures)) {
    stop("`diagnostics` must include both `obs` and `measures`.")
  }

  bundle <- calc_fair_average_bundle(
    res = fit,
    diagnostics = diagnostics,
    facets = facets,
    totalscore = totalscore,
    umean = umean,
    uscale = uscale,
    udecimals = udecimals,
    reference = reference,
    label_style = label_style,
    omit_unobserved = omit_unobserved,
    xtreme = xtreme,
    fair_se = fair_se,
    ci_level = ci_level
  )
  bundle$settings <- list(
    facets = if (is.null(facets)) NULL else as.character(facets),
    totalscore = totalscore,
    umean = umean,
    uscale = uscale,
    udecimals = udecimals,
    reference = reference,
    label_style = label_style,
    omit_unobserved = omit_unobserved,
    xtreme = xtreme,
    fair_se = isTRUE(fair_se),
    ci_level = ci_level,
    model = fit_model,
    method = if (identical(fit_model, "GPCM")) "GPCM-slope-aware" else "PCM/RSM"
  )
  if (identical(fit_model, "GPCM")) {
    bundle$caveat <- paste0(
      "GPCM fair-averages use the slope-aware element-conditional ",
      "construction: each slope-facet element row uses that element's ",
	      "own discrimination, while non-slope-facet rows (persons, raters, ",
	      "...) use the geometric-mean-one slope from the GPCM identification ",
	      "convention. This is an identification-based reporting convention, ",
	      "not FACETS score-side equivalence. Standard errors on the ",
	      "fair-average value are opt-in: ",
      "the SE / Model S.E. / Real S.E. columns are scaled ",
      "facet-measure SEs, not fair-average SEs. Use `fair_se = TRUE` ",
      "to request structural delta-method fair-average SEs for non-person ",
      "rows when the MML observed-information Hessian is available. ",
      "See `gpcm_capability_matrix()` for the current support contract."
    )
  }
  as_mfrm_bundle(bundle, "mfrm_fair_average")
}

#' Compute displacement diagnostics for facet levels
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param facets Optional subset of facets.
#' @param anchored_only If `TRUE`, keep only directly/group anchored levels.
#' @param abs_displacement_warn Absolute displacement warning threshold.
#' @param abs_t_warn Absolute displacement t-value warning threshold.
#' @param top_n Optional maximum number of rows to keep after sorting.
#'
#' @details
#' Displacement is computed as a one-step Newton update:
#' `sum(residual) / sum(information)` for each facet level.
#' This approximates how much a level would move if constraints were relaxed.
#'
#' @section Interpreting output:
#' - `table`: level-wise displacement and flag indicators.
#' - `summary`: count/share of flagged levels.
#' - `thresholds`: displacement and t-value cutoffs.
#'
#' Large absolute displacement in anchored levels suggests potential instability
#' in anchor assumptions.
#'
#' @section Typical workflow:
#' 1. Run `displacement_table(fit, anchored_only = TRUE)` for anchor checks.
#' 2. Inspect `summary(disp)` then detailed rows.
#' 3. Visualize with [plot_displacement()].
#'
#' @section Output columns:
#' The `table` data.frame contains:
#' \describe{
#'   \item{Facet, Level}{Facet name and element label.}
#'   \item{Displacement}{One-step Newton displacement estimate (logits).}
#'   \item{DisplacementSE}{Standard error of the displacement.}
#'   \item{DisplacementT}{Displacement / SE ratio.}
#'   \item{Estimate, SE}{Current measure estimate and its standard error.}
#'   \item{N}{Number of observations involving this level.}
#'   \item{AnchorValue, AnchorStatus, AnchorType}{Anchor metadata.}
#'   \item{Flag}{Logical; `TRUE` when displacement exceeds thresholds.}
#' }
#'
#' @return A named list with:
#' - `table`: displacement diagnostics by level
#' - `summary`: one-row summary
#' - `thresholds`: applied thresholds
#'
#' @seealso [diagnose_mfrm()], [unexpected_response_table()], [fair_average_table()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' disp <- displacement_table(fit, anchored_only = FALSE)
#' summary(disp)
#' p_disp <- plot(disp, draw = FALSE)
#' p_disp$data$plot
#' @export
displacement_table <- function(fit,
                               diagnostics = NULL,
                               facets = NULL,
                               anchored_only = FALSE,
                               abs_displacement_warn = 0.5,
                               abs_t_warn = 2,
                               top_n = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  tbl <- calc_displacement_table(
    obs_df = diagnostics$obs,
    res = fit,
    measures = diagnostics$measures,
    abs_displacement_warn = abs_displacement_warn,
    abs_t_warn = abs_t_warn
  )
  if (!is.null(facets) && nrow(tbl) > 0) {
    tbl <- tbl |>
      filter(.data$Facet %in% as.character(facets))
  }
  if (isTRUE(anchored_only) && nrow(tbl) > 0) {
    tbl <- tbl |>
      filter(.data$AnchorType %in% c("Anchor", "Group"))
  }
  if (!is.null(top_n) && nrow(tbl) > 0) {
    top_n <- max(1L, as.integer(top_n))
    tbl <- tbl |>
      slice_head(n = top_n)
  }

  summary_tbl <- summarize_displacement_table(
    displacement_tbl = tbl,
    abs_displacement_warn = abs_displacement_warn,
    abs_t_warn = abs_t_warn
  )

  out <- list(
    table = tbl,
    summary = summary_tbl,
    thresholds = list(
      abs_displacement_warn = abs_displacement_warn,
      abs_t_warn = abs_t_warn
    )
  )
  as_mfrm_bundle(out, "mfrm_displacement")
}

#' Build a measurable-data summary
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#'
#' @details
#' This helper consolidates measurable-data diagnostics into a dedicated
#' report bundle: run-level summary, facet coverage, category usage, and
#' subset (connected-component) information.
#'
#' `summary(t5)` is supported through `summary()`.
#' `plot(t5)` is dispatched through `plot()` for class
#' `mfrm_measurable` (`type = "facet_coverage"`, `"category_counts"`,
#' `"subset_observations"`).
#'
#' @section Interpreting output:
#' - `summary`: overall measurable design status.
#' - `facet_coverage`: spread/precision by facet.
#' - `category_stats`: category usage and fit context.
#' - `subsets`: connectivity diagnostics (fragmented subsets reduce comparability).
#'
#' @section Typical workflow:
#' 1. Run `measurable_summary_table(fit)`.
#' 2. Check `summary(t5)` for subset/connectivity warnings.
#' 3. Use `plot(t5, ...)` to inspect facet/category/subset views.
#'
#' @section Further guidance:
#' For a plot-selection guide and a longer walkthrough, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @section Output columns:
#' The `summary` data.frame (one row) contains:
#' \describe{
#'   \item{Observations, TotalWeight}{Total observations and summed weight.}
#'   \item{Persons, Facets, Categories}{Design dimensions.}
#'   \item{ConnectedSubsets}{Number of connected subsets.}
#'   \item{LargestSubsetObs, LargestSubsetPct}{Largest subset coverage.}
#' }
#'
#' The `facet_coverage` data.frame contains:
#' \describe{
#'   \item{Facet}{Facet name.}
#'   \item{Levels}{Number of estimated levels.}
#'   \item{MeanSE}{Mean standard error across levels.}
#'   \item{MeanInfit, MeanOutfit}{Mean fit statistics across levels.}
#'   \item{MinEstimate, MaxEstimate}{Measure range for this facet.}
#' }
#'
#' The `category_stats` data.frame contains:
#' \describe{
#'   \item{Category}{Score category value.}
#'   \item{Count, Percent}{Observed count and percentage.}
#'   \item{Infit, Outfit, InfitZSTD, OutfitZSTD}{Category-level fit.}
#'   \item{ExpectedCount, DiffCount, LowCount}{Expected-observed comparison
#'     and low-count flag.}
#' }
#'
#' @return A named list with:
#' - `summary`: one-row measurable-data summary
#' - `facet_coverage`: per-facet coverage summary
#' - `category_stats`: category-level usage/fit summary
#' - `subsets`: subset summary table (when available)
#'
#' @seealso [diagnose_mfrm()], [rating_scale_table()], [describe_mfrm_data()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' t5 <- measurable_summary_table(fit)
#' summary(t5)
#' p_t5 <- plot(t5, draw = FALSE)
#' p_t5$data$plot
#' @export
measurable_summary_table <- function(fit, diagnostics = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  obs_df <- diagnostics$obs
  total_weight <- if ("Weight" %in% names(obs_df)) sum(obs_df$Weight, na.rm = TRUE) else nrow(obs_df)
  subset_tbl <- if (!is.null(diagnostics$subsets$summary)) as.data.frame(diagnostics$subsets$summary) else data.frame()
  subset_n <- if (nrow(subset_tbl) > 0 && "Subset" %in% names(subset_tbl)) nrow(subset_tbl) else NA_integer_
  largest_subset_obs <- if (nrow(subset_tbl) > 0 && "Observations" %in% names(subset_tbl)) max(subset_tbl$Observations, na.rm = TRUE) else NA_real_

  summary_tbl <- data.frame(
    Observations = nrow(obs_df),
    TotalWeight = total_weight,
    Persons = length(fit$prep$levels$Person),
    Facets = length(fit$config$facet_names),
    Categories = fit$prep$rating_max - fit$prep$rating_min + 1,
    ConnectedSubsets = subset_n,
    LargestSubsetObs = largest_subset_obs,
    LargestSubsetPct = ifelse(is.finite(largest_subset_obs) && nrow(obs_df) > 0, 100 * largest_subset_obs / nrow(obs_df), NA_real_),
    stringsAsFactors = FALSE
  )

  facet_coverage <- if (!is.null(diagnostics$measures) && nrow(diagnostics$measures) > 0) {
    as.data.frame(diagnostics$measures) |>
      group_by(.data$Facet) |>
      summarize(
        Levels = n(),
        MeanSE = mean(.data$SE, na.rm = TRUE),
        MeanInfit = mean(.data$Infit, na.rm = TRUE),
        MeanOutfit = mean(.data$Outfit, na.rm = TRUE),
        MinEstimate = min(.data$Estimate, na.rm = TRUE),
        MaxEstimate = max(.data$Estimate, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame()
  }

  category_stats <- as.data.frame(calc_category_stats(obs_df, res = fit, whexact = FALSE), stringsAsFactors = FALSE)

  out <- list(
    summary = summary_tbl,
    facet_coverage = facet_coverage,
    category_stats = category_stats,
    subsets = subset_tbl
  )
  as_mfrm_bundle(out, "mfrm_measurable")
}

#' Build a rating-scale diagnostics report
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param whexact Use exact ZSTD transformation for category fit.
#' @param drop_unused If `TRUE`, remove categories with zero count from the
#'   displayed category table; `summary` and `caveats` still retain the omitted
#'   score-support warning.
#'
#' @details
#' This helper provides category usage/fit statistics and threshold summaries
#' for reviewing score-category functioning.
#' The category usage portion is a global observed-score screen. In PCM fits
#' with a `step_facet`, threshold diagnostics should be interpreted within each
#' `StepFacet` rather than as one pooled whole-scale verdict.
#'
#' Typical checks:
#' - sparse category usage (`Count`, `ExpectedCount`)
#' - category fit (`Infit`, `Outfit`, `ZStd`)
#' - threshold ordering within each `StepFacet`
#'   (`threshold_table$Estimate`, `GapFromPrev`)
#'
#' @section Interpreting output:
#' Start with `summary`:
#' - `UsedCategories` close to total `Categories` suggests that most score
#'   categories are represented in the observed data.
#' - very small `MinCategoryCount` indicates potential instability.
#' - `ThresholdMonotonic = FALSE` indicates disordered thresholds within at
#'   least one threshold set. In PCM fits, inspect `threshold_table` by
#'   `StepFacet` before drawing scale-wide conclusions.
#'
#' Then inspect:
#' - `category_table` for global category-level misfit/sparsity.
#' - `threshold_table` for adjacent-step gaps and ordering within each
#'   `StepFacet`.
#'
#' @section Typical workflow:
#' 1. Fit model: [fit_mfrm()].
#' 2. Build diagnostics: [diagnose_mfrm()].
#' 3. Run `rating_scale_table()` and review `summary()`.
#' 4. Use `plot()` to visualize category profile quickly.
#'
#' @section Further guidance:
#' For a plot-selection guide and a longer walkthrough, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @section Output columns:
#' The `category_table` data.frame contains:
#' \describe{
#'   \item{Category}{Score category value.}
#'   \item{Count, Percent}{Observed count and percentage of total.}
#'   \item{AvgPersonMeasure}{Mean person measure for respondents in this
#'     category.}
#'   \item{Infit, Outfit}{Category-level fit statistics.}
#'   \item{InfitZSTD, OutfitZSTD}{Standardized fit values.}
#'   \item{ExpectedCount, DiffCount}{Expected count and observed-expected
#'     difference.}
#'   \item{LowCount}{Logical; `TRUE` if count is below minimum threshold.}
#'   \item{InfitFlag, OutfitFlag, ZSTDFlag}{Fit-based warning flags.}
#'   \item{ZeroCount, UnusedCategoryType, WeaklyIdentified, CategoryCaveat}{
#'     Structured score-support caveats for retained zero-count categories.}
#' }
#'
#' The `threshold_table` data.frame contains:
#' \describe{
#'   \item{Step}{Step label (e.g., "1-2", "2-3").}
#'   \item{Estimate}{Estimated threshold/step difficulty (logits).}
#'   \item{StepFacet}{Threshold family identifier when the fit uses facet-specific
#'     threshold sets.}
#'   \item{GapFromPrev}{Difference from the previous threshold within the same
#'     `StepFacet` when thresholds are facet-specific. Gaps below
#'     1.4 logits may indicate category underuse; gaps above 5.0 may
#'     indicate wide unused regions (Linacre, 2002).}
#'   \item{ThresholdMonotonic}{Logical flag repeated within each threshold set.
#'     For PCM fits, read this within `StepFacet`, not as a pooled item-bank
#'     verdict.}
#'   \item{LowerCategory, UpperCategory, WeaklyIdentified, ThresholdCaveat}{
#'     Adjacent score-category support metadata. Thresholds adjacent to retained
#'     zero-count categories are flagged for cautious interpretation.}
#' }
#'
#' @return A named list with:
#' - `category_table`: category-level counts, expected counts, fit, and ZSTD
#' - `threshold_table`: model step/threshold estimates
#' - `summary`: one-row summary (usage and threshold monotonicity)
#' - `caveats`: structured score-support warning/review rows
#' - `diagnostic_mode`: character scalar carried from
#'   `diagnostics$diagnostic_mode` (`"legacy"`, `"both"`, or
#'   `"marginal_fit"`); used by downstream reporting helpers to
#'   pick the correct expected-count basis
#' - `marginal_fit`: list bundle from `diagnostics$marginal_fit` when
#'   strict marginal fit was computed, otherwise `NULL`. Carries
#'   the raw OverallRMSD / OverallMaxAbsStdResidual / per-cell
#'   tables that feed the `MarginalOverallRMSD` columns in
#'   `summary`.
#'
#' @seealso [diagnose_mfrm()], [measurable_summary_table()], [plot.mfrm_fit()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' t8 <- rating_scale_table(fit)
#' summary(t8)
#' summary(t8)$summary
#' p_t8 <- plot(t8, draw = FALSE)
#' p_t8$data$plot
#'
#' @section References:
#' - Andrich, D. (1978). *A rating formulation for ordered response
#'   categories*. Psychometrika, 43(4), 561-573.
#'   \doi{10.1007/BF02293814}
#' - Masters, G. N. (1982). *A Rasch model for partial credit scoring*.
#'   Psychometrika, 47(2), 149-174. \doi{10.1007/BF02296272}
#' - Linacre, J. M. (2002). What do Infit and Outfit, mean-square and
#'   standardized mean? *Rasch Measurement Transactions, 16*(2), 878.
#'   (Source for the 0.5-1.5 mean-square acceptance band and the
#'   threshold-gap heuristics used in `summary(t8)$summary`.)
#' - Wind, S. A. (2023). *Detecting rating scale malfunctioning with the
#'   partial credit model and generalized partial credit model*.
#'   Educational and Psychological Measurement, 83(5), 953-983.
#'   \doi{10.1177/00131644221116292} (Recent simulation evidence on
#'   PCM- and GPCM-based rating-scale diagnostics; useful for
#'   interpreting the `summary(t8)$summary` flags in the bounded
#'   `GPCM` route.)
#' @export
rating_scale_table <- function(fit,
                               diagnostics = NULL,
                               whexact = FALSE,
                               drop_unused = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }

  cat_tbl <- as.data.frame(calc_category_stats(diagnostics$obs, res = fit, whexact = whexact), stringsAsFactors = FALSE)
  cat_tbl <- augment_category_table_with_marginal_fit(cat_tbl, diagnostics)
  cat_tbl <- annotate_score_category_caveats(cat_tbl, prep = fit$prep)
  if (isTRUE(drop_unused) && nrow(cat_tbl) > 0 && "Count" %in% names(cat_tbl)) {
    cat_tbl <- cat_tbl[cat_tbl$Count > 0, , drop = FALSE]
  }

  step_tbl <- as.data.frame(fit$steps, stringsAsFactors = FALSE)
  if (nrow(step_tbl) > 0 && all(c("Step", "Estimate") %in% names(step_tbl))) {
    monotonic_flag <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      if (length(x) < 2) {
        return(NA)
      }
      all(diff(x) >= -sqrt(.Machine$double.eps))
    }
    if ("StepFacet" %in% names(step_tbl)) {
      ord <- order(as.character(step_tbl$StepFacet), step_index_from_label(step_tbl$Step))
    } else {
      ord <- order(step_index_from_label(step_tbl$Step))
    }
    step_tbl <- step_tbl[ord, , drop = FALSE]
    est <- suppressWarnings(as.numeric(step_tbl$Estimate))
    if ("StepFacet" %in% names(step_tbl)) {
      groups <- as.character(step_tbl$StepFacet)
      step_tbl$GapFromPrev <- stats::ave(est, groups, FUN = function(x) c(NA_real_, diff(x)))
      step_tbl$ThresholdMonotonic <- stats::ave(
        est,
        groups,
        FUN = function(x) rep(monotonic_flag(x), length(x))
      )
    } else {
      step_tbl$GapFromPrev <- c(NA_real_, diff(est))
      step_tbl$ThresholdMonotonic <- rep(monotonic_flag(est), nrow(step_tbl))
    }
    step_tbl <- annotate_threshold_caveats(step_tbl, prep = fit$prep)
  }

  threshold_monotonic <- if (nrow(step_tbl) > 1 && "Estimate" %in% names(step_tbl)) {
    if ("StepFacet" %in% names(step_tbl)) {
      group_flags <- vapply(split(step_tbl$Estimate, step_tbl$StepFacet), monotonic_flag, logical(1))
      if (any(is.na(group_flags))) NA else all(group_flags)
    } else {
      monotonic_flag(step_tbl$Estimate)
    }
  } else {
    NA
  }

  marginal_fit_available <- has_marginal_fit_bundle(diagnostics)
  marginal_summary <- as.data.frame(diagnostics$marginal_fit$summary %||% data.frame(), stringsAsFactors = FALSE)
  marginal_flagged_categories <- if ("MarginalFitFlag" %in% names(cat_tbl)) {
    sum(as.logical(cat_tbl$MarginalFitFlag), na.rm = TRUE)
  } else {
    NA_integer_
  }
  caveats <- collect_mfrm_caveats(fit = fit)
  prep_unused <- suppressWarnings(as.numeric(fit$prep$unused_score_categories %||% numeric(0)))
  prep_unused <- prep_unused[is.finite(prep_unused)]
  unused_score_categories <- if (length(prep_unused) > 0L) {
    paste(as.character(prep_unused), collapse = ", ")
  } else if ("ZeroCount" %in% names(cat_tbl) && "Category" %in% names(cat_tbl)) {
    paste(as.character(cat_tbl$Category[as.logical(cat_tbl$ZeroCount)]), collapse = ", ")
  } else {
    ""
  }
  weak_thresholds <- if ("WeaklyIdentified" %in% names(step_tbl)) {
    sum(as.logical(step_tbl$WeaklyIdentified), na.rm = TRUE)
  } else {
    NA_integer_
  }

  summary_tbl <- data.frame(
    Categories = nrow(cat_tbl),
    UsedCategories = if ("Count" %in% names(cat_tbl)) sum(cat_tbl$Count > 0, na.rm = TRUE) else NA_integer_,
    UnusedScoreCategories = unused_score_categories,
    WeaklyIdentifiedThresholds = weak_thresholds,
    MinCategoryCount = if ("Count" %in% names(cat_tbl) && nrow(cat_tbl) > 0) min(cat_tbl$Count, na.rm = TRUE) else NA_real_,
    MaxCategoryCount = if ("Count" %in% names(cat_tbl) && nrow(cat_tbl) > 0) max(cat_tbl$Count, na.rm = TRUE) else NA_real_,
    MeanCategoryInfit = if ("Infit" %in% names(cat_tbl)) mean(cat_tbl$Infit, na.rm = TRUE) else NA_real_,
    MeanCategoryOutfit = if ("Outfit" %in% names(cat_tbl)) mean(cat_tbl$Outfit, na.rm = TRUE) else NA_real_,
    ThresholdMonotonic = threshold_monotonic,
    DiagnosticMode = as.character(diagnostics$diagnostic_mode %||% "legacy"),
    ExpectedCountBasis = if (marginal_fit_available) {
      "legacy_plugin + latent_integrated_first_order_counts"
    } else {
      "legacy_plugin"
    },
    MarginalFitAvailable = marginal_fit_available,
    MarginalOverallRMSD = if (marginal_fit_available) marginal_summary$OverallRMSD[1] %||% NA_real_ else NA_real_,
    MarginalMaxAbsStdResidual = if (marginal_fit_available) marginal_summary$OverallMaxAbsStdResidual[1] %||% NA_real_ else NA_real_,
    MarginalFlaggedCategories = marginal_flagged_categories,
    stringsAsFactors = FALSE
  )

  out <- list(
    category_table = cat_tbl,
    threshold_table = step_tbl,
    summary = summary_tbl,
    caveats = caveats,
    diagnostic_mode = as.character(diagnostics$diagnostic_mode %||% "legacy"),
    marginal_fit = diagnostics$marginal_fit %||% NULL
  )
  as_mfrm_bundle(out, "mfrm_rating_scale")
}

#' Build a bias-cell count report
#'
#' @param bias_results Output from [estimate_bias()].
#' @param min_count_warn Minimum count threshold for flagging sparse bias cells.
#' @param branch Output branch:
#'   `"facets"` keeps legacy manual-aligned naming, `"original"` returns compact QC-oriented names.
#' @param fit Optional [fit_mfrm()] result used to attach run context metadata.
#'
#' @details
#' This helper summarizes how many observations contribute to each
#' bias-cell estimate and flags sparse cells.
#'
#' Branch behavior:
#' - `"facets"`: keeps legacy manual-aligned column labels (`Sq`,
#'   `Observd Count`, `Obs-Exp Average`, `Model S.E.`) for side-by-side
#'   comparison with external workflows.
#' - `"original"`: keeps compact field names (`Count`, `BiasSize`, `SE`) for
#'   custom QC workflows and scripting.
#'
#' @section Interpreting output:
#' - `table`: cell-level contribution counts and low-count flags.
#' - `by_facet`: sparse-cell structure by each interaction facet.
#' - `summary`: overall low-count prevalence.
#' - `fit_overview`: optional run context (when `fit` is supplied).
#'
#' Low-count cells should be interpreted cautiously because bias-size estimates
#' can become unstable with sparse support.
#'
#' @section Typical workflow:
#' 1. Estimate bias with [estimate_bias()].
#' 2. Build `bias_count_table(...)` in desired branch.
#' 3. Review low-count flags before interpreting bias magnitudes.
#'
#' @section Further guidance:
#' For a plot-selection guide and a longer walkthrough, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @section Output columns:
#' The `table` data.frame contains, in the legacy-compatible branch:
#' \describe{
#'   \item{FacetA, FacetB}{Interaction facet level identifiers; generic
#'     names for the two interaction facets.}
#'   \item{Sq}{Sequential row number.}
#'   \item{Observd Count}{Number of observations for this cell.}
#'   \item{Obs-Exp Average}{Observed minus expected average for this cell.}
#'   \item{Model S.E.}{Standard error of the bias estimate.}
#'   \item{Infit, Outfit}{Fit statistics for this cell.}
#'   \item{LowCountFlag}{Logical; `TRUE` when count < `min_count_warn`.}
#' }
#'
#' The `summary` data.frame contains:
#' \describe{
#'   \item{InteractionFacets}{Names of the interaction facets.}
#'   \item{Cells, TotalCount}{Number of cells and total observations.}
#'   \item{LowCountCells, LowCountPercent}{Number and share of low-count
#'     cells.}
#' }
#'
#' @return A named list with:
#' - `table`: cell-level counts with low-count flags
#' - `by_facet`: named list of counts aggregated by each interaction facet
#' - `by_facet_a`, `by_facet_b`: first two facet summaries (legacy compatibility)
#' - `summary`: one-row summary
#' - `thresholds`: applied thresholds
#' - `branch`, `style`: output branch metadata
#' - `fit_overview`: optional one-row fit metadata when `fit` is supplied
#'
#' @seealso [estimate_bias()], [unexpected_after_bias_table()], [build_fixed_reports()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' t11 <- bias_count_table(bias)
#' t11_facets <- bias_count_table(bias, branch = "facets", fit = fit)
#' summary(t11)
#' p <- plot(t11, draw = FALSE)
#' p2 <- plot(t11, type = "lowcount_by_facet", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     t11,
#'     type = "cell_counts",
#'     draw = TRUE,
#'     main = "Bias Cell Counts (Customized)",
#'     palette = c(count = "#2b8cbe", low = "#cb181d"),
#'     label_angle = 45
#'   )
#' }
#' @export
bias_count_table <- function(bias_results,
                             min_count_warn = 10,
                             branch = c("original", "facets"),
                             fit = NULL) {
  branch <- match.arg(tolower(as.character(branch[1])), c("original", "facets"))

  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    stop("`bias_results` must be output from estimate_bias() with non-empty `table`.")
  }
  spec <- extract_bias_facet_spec(bias_results)
  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)
  if (is.null(spec) || length(spec$facets) < 2) {
    stop("`bias_results$table` does not include recognizable interaction facet columns.")
  }
  req_cols <- c(spec$level_cols, "Observd Count", "Bias Size", "S.E.")
  if (!all(req_cols %in% names(tbl))) {
    stop("`bias_results$table` does not include required columns.")
  }

  min_count_warn <- max(0, as.numeric(min_count_warn))
  level_tbl <- tbl[, spec$level_cols, drop = FALSE]
  level_tbl[] <- lapply(level_tbl, as.character)
  names(level_tbl) <- spec$facets

  cell_tbl <- dplyr::bind_cols(
    level_tbl,
    tbl |>
      dplyr::transmute(
        Count = suppressWarnings(as.numeric(.data$`Observd Count`)),
        BiasSize = suppressWarnings(as.numeric(.data$`Bias Size`)),
        SE = suppressWarnings(as.numeric(.data$`S.E.`)),
        Infit = suppressWarnings(as.numeric(.data$Infit)),
        Outfit = suppressWarnings(as.numeric(.data$Outfit))
      )
  ) |>
    dplyr::mutate(
      LowCountFlag = is.finite(.data$Count) & .data$Count < min_count_warn
    ) |>
    dplyr::arrange(.data$Count)

  by_facet <- stats::setNames(vector("list", length(spec$facets)), spec$facets)
  for (facet in spec$facets) {
    by_facet[[facet]] <- cell_tbl |>
      dplyr::group_by(Level = .data[[facet]]) |>
      dplyr::summarize(
        Cells = dplyr::n(),
        TotalCount = sum(.data$Count, na.rm = TRUE),
        MeanCount = mean(.data$Count, na.rm = TRUE),
        MinCount = min(.data$Count, na.rm = TRUE),
        LowCountCells = sum(.data$LowCountFlag, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(Facet = facet, .before = 1)
  }

  summary_tbl <- data.frame(
    InteractionFacets = paste(spec$facets, collapse = " x "),
    InteractionOrder = spec$interaction_order,
    InteractionMode = spec$interaction_mode,
    Branch = branch,
    Style = ifelse(branch == "facets", "facets_manual", "original"),
    FacetA = spec$facets[1],
    FacetB = spec$facets[2],
    Cells = nrow(cell_tbl),
    TotalCount = sum(cell_tbl$Count, na.rm = TRUE),
    MeanCount = mean(cell_tbl$Count, na.rm = TRUE),
    MedianCount = stats::median(cell_tbl$Count, na.rm = TRUE),
    MinCount = min(cell_tbl$Count, na.rm = TRUE),
    MaxCount = max(cell_tbl$Count, na.rm = TRUE),
    LowCountCells = sum(cell_tbl$LowCountFlag, na.rm = TRUE),
    LowCountPercent = ifelse(nrow(cell_tbl) > 0, 100 * sum(cell_tbl$LowCountFlag, na.rm = TRUE) / nrow(cell_tbl), NA_real_),
    stringsAsFactors = FALSE
  )

  table_out <- cell_tbl
  if (branch == "facets") {
    table_out <- cell_tbl |>
      dplyr::mutate(
        Sq = dplyr::row_number(),
        `Observd Count` = .data$Count,
        `Obs-Exp Average` = .data$BiasSize,
        `Model S.E.` = .data$SE
      ) |>
      dplyr::select(
        dplyr::all_of(spec$facets),
        "Sq",
        "Observd Count",
        "Obs-Exp Average",
        "Model S.E.",
        "Infit",
        "Outfit",
        "LowCountFlag"
      )
  }

  fit_overview <- data.frame()
  if (!is.null(fit) && inherits(fit, "mfrm_fit") &&
      is.data.frame(fit$summary) && nrow(fit$summary) > 0) {
    fit_overview <- fit$summary[1, , drop = FALSE]
    fit_overview$InteractionFacets <- paste(spec$facets, collapse = " x ")
    fit_overview$Branch <- branch
    fit_overview <- fit_overview[, c("InteractionFacets", "Branch", setdiff(names(fit_overview), c("InteractionFacets", "Branch"))), drop = FALSE]
  }

  out <- list(
    table = table_out,
    by_facet = by_facet,
    by_facet_a = by_facet[[spec$facets[1]]],
    by_facet_b = by_facet[[spec$facets[2]]],
    summary = summary_tbl,
    thresholds = list(min_count_warn = min_count_warn),
    branch = branch,
    style = ifelse(branch == "facets", "facets_manual", "original"),
    fit_overview = fit_overview
  )
  out <- as_mfrm_bundle(out, "mfrm_bias_count")
  class(out) <- unique(c(paste0("mfrm_bias_count_", branch), class(out)))
  out
}

#' Build an unexpected-after-adjustment screening report
#'
#' @param fit Output from [fit_mfrm()].
#' @param bias_results Output from [estimate_bias()].
#' @param diagnostics Optional output from [diagnose_mfrm()] for baseline comparison.
#' @param abs_z_min Absolute standardized-residual cutoff.
#' @param prob_max Maximum observed-category probability cutoff.
#' @param top_n Maximum number of rows to return.
#' @param rule Flagging rule: `"either"` or `"both"`.
#'
#' @details
#' This helper recomputes expected values and residuals after interaction
#' adjustments from [estimate_bias()] have been introduced.
#'
#' `summary(t10)` is supported through `summary()`.
#' `plot(t10)` is dispatched through `plot()` for class
#' `mfrm_unexpected_after_bias` (`type = "scatter"`, `"severity"`,
#' `"comparison"`).
#'
#' @section Interpreting output:
#' - `summary`: before/after unexpected counts and reduction metrics.
#' - `table`: residual unexpected responses after bias adjustment.
#' - `thresholds`: screening settings used in this comparison.
#'
#' Large reductions indicate bias terms explain part of prior unexpectedness;
#' persistent unexpected rows indicate remaining model-data mismatch.
#'
#' @section Typical workflow:
#' 1. Run [unexpected_response_table()] as baseline.
#' 2. Estimate bias via [estimate_bias()].
#' 3. Run `unexpected_after_bias_table(...)` and compare reductions.
#'
#' @section Further guidance:
#' For a plot-selection guide and a longer walkthrough, see
#' [mfrmr_visual_diagnostics] and
#' `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
#'
#' @section Output columns:
#' The `table` data.frame has the same structure as
#' [unexpected_response_table()] output, with an additional
#' `BiasAdjustment` column showing the bias correction applied to each
#' observation's expected value.
#'
#' The `summary` data.frame contains:
#' \describe{
#'   \item{TotalObservations}{Total observations analyzed.}
#'   \item{BaselineUnexpectedN}{Unexpected count before bias adjustment.}
#'   \item{AfterBiasUnexpectedN}{Unexpected count after adjustment.}
#'   \item{ReducedBy, ReducedPercent}{Reduction in unexpected count.}
#' }
#'
#' @return A named list with:
#' - `table`: unexpected responses after bias adjustment
#' - `summary`: one-row summary (includes baseline-vs-after counts)
#' - `thresholds`: applied thresholds
#' - `facets`: analyzed bias facet pair
#'
#' @seealso [estimate_bias()], [unexpected_response_table()], [bias_count_table()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' t10 <- unexpected_after_bias_table(fit, bias, diagnostics = diag, top_n = 20)
#' summary(t10)
#' p_t10 <- plot(t10, draw = FALSE)
#' p_t10$data$plot
#' @export
unexpected_after_bias_table <- function(fit,
                                        bias_results,
                                        diagnostics = NULL,
                                        abs_z_min = 2,
                                        prob_max = 0.30,
                                        top_n = 100,
                                        rule = c("either", "both")) {
  rule <- match.arg(tolower(rule), c("either", "both"))
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  stop_if_gpcm_out_of_scope(fit, "unexpected_after_bias_table()")
  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    stop("`bias_results` must be output from estimate_bias() with non-empty `table`.")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }

  obs_adj <- compute_obs_table_with_bias(fit, bias_results = bias_results)
  probs_adj <- compute_prob_matrix_with_bias(fit, bias_results = bias_results)

  tbl <- calc_unexpected_response_table(
    obs_df = obs_adj,
    probs = probs_adj,
    facet_names = fit$config$facet_names,
    rating_min = fit$prep$rating_min,
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    top_n = top_n,
    rule = rule
  )
  if (nrow(tbl) > 0 && "Row" %in% names(tbl) && "BiasAdjustment" %in% names(obs_adj)) {
    adj_df <- data.frame(Row = seq_len(nrow(obs_adj)), BiasAdjustment = obs_adj$BiasAdjustment, stringsAsFactors = FALSE)
    tbl <- tbl |>
      left_join(adj_df, by = "Row")
  }

  summary_tbl <- summarize_unexpected_response_table(
    unexpected_tbl = tbl,
    total_observations = nrow(obs_adj),
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    rule = rule
  )

  baseline <- unexpected_response_table(
    fit = fit,
    diagnostics = diagnostics,
    abs_z_min = abs_z_min,
    prob_max = prob_max,
    top_n = max(top_n, nrow(obs_adj)),
    rule = rule
  )
  baseline_n <- if (is.null(baseline$table)) NA_integer_ else nrow(baseline$table)
  after_n <- nrow(tbl)
  summary_tbl <- summary_tbl |>
    mutate(
      BaselineUnexpectedN = baseline_n,
      AfterBiasUnexpectedN = after_n,
      ReducedBy = ifelse(is.finite(baseline_n), baseline_n - after_n, NA_real_),
      ReducedPercent = ifelse(is.finite(baseline_n) && baseline_n > 0, 100 * (baseline_n - after_n) / baseline_n, NA_real_)
    )

  out <- list(
    table = tbl,
    summary = summary_tbl,
    thresholds = list(
      abs_z_min = abs_z_min,
      prob_max = prob_max,
      rule = rule
    ),
    facets = list(
      facet_a = bias_results$facet_a,
      facet_b = bias_results$facet_b,
      interaction_facets = bias_results$interaction_facets %||% c(bias_results$facet_a, bias_results$facet_b),
      interaction_order = bias_results$interaction_order %||% 2L,
      interaction_mode = bias_results$interaction_mode %||% "pairwise"
    )
  )
  as_mfrm_bundle(out, "mfrm_unexpected_after_bias")
}

resolve_table2_source_columns <- function(fit,
                                          person = NULL,
                                          facets = NULL,
                                          score = NULL,
                                          weight = NULL) {
  source <- fit$config$source_columns
  resolved_person <- if (!is.null(person)) as.character(person[1]) else as.character(source$person %||% "Person")
  resolved_facets <- if (!is.null(facets)) as.character(facets) else as.character(source$facets %||% fit$config$facet_names)
  resolved_score <- if (!is.null(score)) as.character(score[1]) else as.character(source$score %||% "Score")
  resolved_weight <- if (!is.null(weight)) as.character(weight[1]) else as.character(source$weight %||% NA_character_)
  if (length(resolved_weight) == 0 || is.na(resolved_weight) || !nzchar(resolved_weight)) {
    resolved_weight <- NA_character_
  }
  list(
    person = resolved_person,
    facets = resolved_facets,
    score = resolved_score,
    weight = resolved_weight
  )
}

compute_iteration_state <- function(par, idx, prep, config, sizes, quad_points = 15L) {
  params <- expand_params(par, sizes, config)
  if (config$method == "MML") {
    quad <- gauss_hermite_normal(quad_points)
    theta_tbl <- compute_person_eap(idx, config, params, quad)
    theta_diag <- suppressWarnings(as.numeric(theta_tbl$Estimate))
    if (length(theta_diag) != config$n_person || any(!is.finite(theta_diag))) {
      theta_diag <- rep(0, config$n_person)
    }
  } else {
    theta_diag <- suppressWarnings(as.numeric(params$theta))
  }

  eta <- compute_eta(idx, params, config, theta_override = theta_diag)
  if (config$model == "RSM") {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
  } else if (identical(config$model, "GPCM")) {
    slope_idx <- idx$slope_idx %||% idx$step_idx
    if (is.null(slope_idx)) {
      stop("GPCM iteration state requires slope-facet indices.", call. = FALSE)
    }
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_gpcm(
      eta = eta,
      step_cum_mat = step_cum_mat,
      criterion_idx = idx$step_idx,
      slopes = params$slopes,
      slope_idx = slope_idx
    )
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(eta, step_cum_mat, idx$step_idx)
  }
  k_vals <- 0:(ncol(probs) - 1)
  expected_k <- as.vector(probs %*% k_vals)
  expected <- prep$rating_min + expected_k
  observed <- as.numeric(prep$data$Score)
  weight <- if (!is.null(idx$weight)) as.numeric(idx$weight) else rep(1, length(observed))

  obs_df <- prep$data |>
    mutate(
      Observed = observed,
      Expected = expected,
      Residual = .data$Observed - .data$Expected,
      .Weight = weight
    )

  facet_cols <- c("Person", config$facet_names)
  max_elem_resid <- 0
  for (facet in facet_cols) {
    if (!facet %in% names(obs_df)) next
    sub <- obs_df |>
      group_by(.data[[facet]]) |>
      summarize(ResidualSum = sum(.data$Residual * .data$.Weight, na.rm = TRUE), .groups = "drop")
    if (nrow(sub) == 0) next
    idx_max <- which.max(abs(sub$ResidualSum))
    v <- sub$ResidualSum[idx_max]
    if (is.finite(v) && abs(v) > abs(max_elem_resid)) max_elem_resid <- v
  }
  score_range <- prep$rating_max - prep$rating_min
  max_elem_resid_pct <- ifelse(score_range > 0, 100 * max_elem_resid / score_range, NA_real_)

  obs_k <- observed - prep$rating_min
  n_cat <- ncol(probs)
  obs_count <- rep(0, n_cat)
  valid_k <- is.finite(obs_k) & obs_k >= 0 & obs_k < n_cat
  if (any(valid_k)) {
    obs_idx <- as.integer(obs_k[valid_k]) + 1L
    obs_count <- tapply(weight[valid_k], obs_idx, sum)
    obs_count <- {
      out <- rep(0, n_cat)
      out[as.integer(names(obs_count))] <- as.numeric(obs_count)
      out
    }
  }
  exp_count <- colSums(probs * weight, na.rm = TRUE)
  cat_diff <- obs_count - exp_count
  if (length(cat_diff) == 0 || all(!is.finite(cat_diff))) {
    max_cat_resid <- NA_real_
  } else {
    max_cat_resid <- cat_diff[which.max(abs(cat_diff))]
  }

  facet_params <- unlist(params$facets, use.names = FALSE)
  element_vec <- c(theta_diag, facet_params)
  step_vec <- if (config$model == "RSM") {
    as.numeric(params$steps)
  } else if (identical(config$model, "GPCM")) {
    c(as.numeric(as.vector(params$steps_mat)), as.numeric(params$log_slopes))
  } else {
    as.numeric(as.vector(params$steps_mat))
  }

  list(
    max_score_resid_elements = as.numeric(max_elem_resid),
    max_score_resid_pct = as.numeric(max_elem_resid_pct),
    max_score_resid_categories = as.numeric(max_cat_resid),
    element_vec = element_vec,
    step_vec = step_vec
  )
}

signal_legacy_name_deprecation <- function(old_name,
                                           new_name,
                                           suppress_if_called_from = NULL,
                                           when = "0.1.0") {
  if (isTRUE(getOption("mfrmr.suppress_legacy_name_warning", FALSE))) {
    return(invisible(NULL))
  }

  caller_call <- tryCatch(sys.call(-2), error = function(e) NULL)
  caller <- ""
  if (!is.null(caller_call) && length(caller_call) > 0) {
    caller <- as.character(caller_call[[1]])
    if (length(caller) == 0 || !is.finite(nchar(caller[1])) || is.na(caller[1])) {
      caller <- ""
    } else {
      caller <- caller[1]
    }
  }
  if (!is.null(suppress_if_called_from) &&
      is.character(suppress_if_called_from) &&
      length(suppress_if_called_from) > 0 &&
      nzchar(caller) &&
      caller %in% suppress_if_called_from) {
    return(invisible(NULL))
  }

  dep_env <- tryCatch(rlang::caller_env(2), error = function(e) rlang::caller_env())
  user_env <- tryCatch(rlang::caller_env(3), error = function(e) dep_env)
  lifecycle::deprecate_soft(
    when = when,
    what = paste0(old_name, "()"),
    with = paste0(new_name, "()"),
    id = paste0("mfrmr_", old_name, "_legacy_name"),
    env = dep_env,
    user_env = user_env
  )
  invisible(NULL)
}

with_legacy_name_warning_suppressed <- function(expr) {
  old_opt <- options(mfrmr.suppress_legacy_name_warning = TRUE)
  on.exit(options(old_opt), add = TRUE)
  eval.parent(substitute(expr))
}

as_mfrm_bundle <- function(x, class_name) {
  if (!is.list(x)) return(x)
  class(x) <- unique(c(class_name, "mfrm_bundle", class(x)))
  x
}

#' Build a legacy-compatible Table 1 specification summary
#'
#' @param fit Output from [fit_mfrm()].
#' @param title Optional analysis title.
#' @param data_file Optional data-file label (for reporting only).
#' @param output_file Optional output-file label (for reporting only).
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#'
#' @details
#' The legacy-compatible Table 1 layout groups model settings by function
#' (data, output, convergence).
#' This helper assembles those settings from a fitted object.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [specifications_report()].
#'
#' @section Output columns:
#' The `header` data.frame contains:
#' \describe{
#'   \item{Engine}{Estimation engine identifier (always `"mfrmr"`).}
#'   \item{Model}{`"RSM"` or `"PCM"`.}
#'   \item{Method}{`"JMLE"` or `"MML"`.}
#' }
#'
#' The `facet_labels` data.frame contains:
#' \describe{
#'   \item{Facet}{Facet column name.}
#'   \item{Elements}{Number of levels in this facet.}
#' }
#'
#' @return A named list with:
#' - `header`: one-row run header
#' - `data_spec`: key-value table for data/model settings
#' - `facet_labels`: facet names with level counts
#' - `output_spec`: key-value table for reporting-related defaults
#' - `convergence_control`: key-value table for optimizer controls/results
#' - `anchor_summary`: anchor/group-anchor summary by facet
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#'
#' @seealso [fit_mfrm()], [data_quality_report()], [estimation_iteration_report()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' t1 <- specifications_report(fit, title = "Toy run")
#' @keywords internal
#' @noRd
table1_specifications <- function(fit,
                                  title = NULL,
                                  data_file = NULL,
                                  output_file = NULL,
                                  include_fixed = FALSE) {
  signal_legacy_name_deprecation(
    old_name = "table1_specifications",
    new_name = "specifications_report",
    suppress_if_called_from = "specifications_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }

  cfg <- fit$config
  prep <- fit$prep
  ov <- as.data.frame(fit$summary, stringsAsFactors = FALSE)
  if (nrow(ov) == 0) stop("`fit$summary` is empty.")

  header <- data.frame(
    Engine = paste0("mfrmr ", as.character(utils::packageVersion("mfrmr"))),
    Title = ifelse(is.null(title), "", as.character(title[1])),
    DataFile = ifelse(is.null(data_file), "", as.character(data_file[1])),
    OutputFile = ifelse(is.null(output_file), "", as.character(output_file[1])),
    Model = as.character(ov$Model[1]),
    Method = as.character(ov$Method[1]),
    stringsAsFactors = FALSE
  )

  all_facets <- c("Person", cfg$facet_names)
  facet_labels <- data.frame(
    FacetIndex = seq_along(all_facets),
    Facet = all_facets,
    Elements = vapply(prep$levels[all_facets], length, integer(1)),
    stringsAsFactors = FALSE
  )

  est_ctl <- cfg$estimation_control %||% list()
  positive <- cfg$positive_facets %||% character(0)
  dummy <- cfg$dummy_facets %||% character(0)
  data_spec <- data.frame(
    Setting = c(
      "Facets",
      "Persons",
      "Categories",
      "RatingMin",
      "RatingMax",
      "NonCenteredFacet",
      "PositiveFacets",
      "DummyFacets",
      "StepFacet",
      "WeightColumn"
    ),
    Value = c(
      as.character(length(cfg$facet_names)),
      as.character(cfg$n_person),
      as.character(cfg$n_cat),
      as.character(prep$rating_min),
      as.character(prep$rating_max),
      as.character(cfg$noncenter_facet %||% ""),
      ifelse(length(positive) == 0, "", paste(positive, collapse = ", ")),
      ifelse(length(dummy) == 0, "", paste(dummy, collapse = ", ")),
      as.character(cfg$step_facet %||% ""),
      as.character(cfg$weight_col %||% "")
    ),
    stringsAsFactors = FALSE
  )

  output_spec <- data.frame(
    Setting = c(
      "UnexpectedAbsZThreshold",
      "UnexpectedProbThreshold",
      "DisplacementWarnLogit",
      "DisplacementWarnT",
      "FairScoreDefault"
    ),
    Value = c("2", "0.30", "0.5", "2", "Mean"),
    stringsAsFactors = FALSE
  )

  convergence_control <- data.frame(
    Setting = c(
      "MaxIterations",
      "RelativeTolerance",
      "QuadPoints",
      "Converged",
      "FunctionEvaluations",
      "OptimizerMessage"
    ),
    Value = c(
      as.character(est_ctl$maxit %||% NA_integer_),
      as.character(est_ctl$reltol %||% NA_real_),
      as.character(est_ctl$quad_points %||% NA_integer_),
      as.character(isTRUE(ov$Converged[1])),
      as.character(fit$opt$counts[["function"]] %||% NA_integer_),
      as.character(fit$opt$message %||% "")
    ),
    stringsAsFactors = FALSE
  )

  anchor_summary <- if (!is.null(cfg$anchor_summary)) {
    as.data.frame(cfg$anchor_summary, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }

  out <- list(
    header = header,
    data_spec = data_spec,
    facet_labels = facet_labels,
    output_spec = output_spec,
    convergence_control = convergence_control,
    anchor_summary = anchor_summary
  )

  if (isTRUE(include_fixed)) {
    out$fixed <- build_sectioned_fixed_report(
      title = "Legacy-compatible Table 1 Specification Summary",
      sections = list(
        list(title = "Header", data = header),
        list(title = "Data specification", data = data_spec),
        list(title = "Facet labels", data = facet_labels),
        list(title = "Output specification", data = output_spec),
        list(title = "Convergence control", data = convergence_control),
        list(title = "Anchor summary", data = anchor_summary)
      ),
      max_col_width = 24
    )
  }
  out
}

#' Build a legacy-compatible Table 2 data summary report
#'
#' @param fit Output from [fit_mfrm()].
#' @param data Optional raw data frame used for additional row-level review.
#' @param person Optional person column name in `data`.
#' @param facets Optional facet column names in `data`.
#' @param score Optional score column name in `data`.
#' @param weight Optional weight column name in `data`.
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#'
#' @details
#' When `data` is supplied, this function performs row-level validity checks
#' (missing identifiers, missing score, non-positive weight, out-of-range score)
#' and reports dropped rows in a legacy-compatible Table 2 layout.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [data_quality_report()].
#'
#' @section Output columns:
#' The `summary` data.frame contains:
#' \describe{
#'   \item{TotalLinesInData}{Total rows in the input data.}
#'   \item{ValidResponsesUsedForEstimation}{Rows used after filtering.}
#'   \item{MissingScoreRows, MissingFacetRows, MissingPersonRows}{Counts
#'     of rows excluded due to missing values.}
#' }
#'
#' The `row_review` data.frame contains:
#' \describe{
#'   \item{Status}{Row-status category (e.g., "Valid", "MissingScore").}
#'   \item{N}{Number of rows in this category.}
#' }
#'
#' @return A named list with:
#' - `summary`: one-row data summary
#' - `quality_overview`: compact area-level QC status summary
#' - `quality_flags`: prioritized user-facing QC flags and next actions
#' - `model_match`: model-match counts
#' - `row_review`: counts by row-status category
#' - `unknown_elements`: elements in `data` not present in fitted levels
#' - `category_counts`: observed category counts over the retained score support
#' - `score_support_review`: zero-count category and weak-threshold caveats
#' - `category_usage_by_facet`: facet-level category usage over the retained
#'   score support
#' - `category_usage_summary`: per-facet-level zero/sparse category summary
#' - `facet_response_patterns`: facet-level response-pattern summaries,
#'   including single-category and dominant-category use
#' - `caveats`: user-facing score-support warnings
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#'
#' @seealso [fit_mfrm()], [specifications_report()], [describe_mfrm_data()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' t2 <- data_quality_report(
#'   fit, data = toy, person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score"
#' )
#' @keywords internal
#' @noRd
table2_data_summary <- function(fit,
                                data = NULL,
                                person = NULL,
                                facets = NULL,
                                score = NULL,
                                weight = NULL,
                                min_category_count = 10,
                                dominant_category_cutoff = 0.95,
                                include_fixed = FALSE) {
  signal_legacy_name_deprecation(
    old_name = "table2_data_summary",
    new_name = "data_quality_report",
    suppress_if_called_from = "data_quality_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  min_category_count <- validate_data_quality_min_category_count(min_category_count)
  dominant_category_cutoff <- validate_data_quality_dominant_category_cutoff(dominant_category_cutoff)

  prep <- fit$prep
  cfg <- fit$config
  src <- resolve_table2_source_columns(
    fit = fit,
    person = person,
    facets = facets,
    score = score,
    weight = weight
  )

  fitted_df <- prep$data
  valid_used <- nrow(fitted_df)
  observed_category_counts <- fitted_df |>
    group_by(.data$Score) |>
    summarize(
      Count = n(),
      WeightedCount = sum(.data$Weight, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(.data$Score)
  category_counts <- build_data_quality_category_counts(
    prep = prep,
    observed_counts = observed_category_counts
  )
  score_support_review <- build_data_quality_score_support_review(category_counts)
  category_usage_by_facet <- build_data_quality_facet_category_usage(
    prep = prep,
    category_counts = category_counts,
    min_category_count = min_category_count
  )
  category_usage_summary <- build_data_quality_facet_category_summary(category_usage_by_facet)
  facet_response_patterns <- build_data_quality_facet_response_patterns(
    category_usage_by_facet = category_usage_by_facet,
    dominant_category_cutoff = dominant_category_cutoff
  )
  caveats <- collect_mfrm_caveats(
    fit = fit,
    score_distribution = category_counts,
    include_recode = TRUE
  )
  score_map <- as.data.frame(prep$score_map %||% data.frame(), stringsAsFactors = FALSE)

  if (is.null(data)) {
    summary_tbl <- data.frame(
      TotalLinesInData = NA_integer_,
      TotalDataLines = NA_integer_,
      TotalNonBlankResponsesFound = NA_integer_,
      ValidResponsesUsedForEstimation = valid_used,
      ZeroCountScoreCategories = sum(category_counts$ZeroCount %in% TRUE, na.rm = TRUE),
      IntermediateZeroCountScoreCategories = sum(category_counts$UnusedCategoryType == "internal", na.rm = TRUE),
      FacetLevelsWithZeroCategories = sum(category_usage_summary$ZeroCategories > 0, na.rm = TRUE),
      FacetLevelsWithIntermediateZeroCategories = sum(category_usage_summary$IntermediateZeroCategories > 0, na.rm = TRUE),
      FacetLevelsWithSparseCategories = sum(category_usage_summary$SparseCategories > 0, na.rm = TRUE),
      FacetLevelsWithSingleCategoryUse = sum(facet_response_patterns$SingleCategoryUse %in% TRUE, na.rm = TRUE),
      FacetLevelsWithDominantCategoryUse = sum(facet_response_patterns$DominantCategoryUse %in% TRUE, na.rm = TRUE),
      FacetLevelsWithBoundaryOnlyUse = sum(facet_response_patterns$BoundaryOnlyUse %in% TRUE, na.rm = TRUE),
      ScoreSupportCaveats = nrow(caveats),
      stringsAsFactors = FALSE
    )
    model_match <- data.frame(
      Model = paste(cfg$model, cfg$method, sep = " / "),
      MatchedResponses = valid_used,
      stringsAsFactors = FALSE
    )
    quality_flags <- build_data_quality_quality_flags(
      summary_tbl = summary_tbl,
      row_review = data.frame(),
      unknown_elements = data.frame(),
      category_counts = category_counts,
      category_usage_summary = category_usage_summary,
      facet_response_patterns = facet_response_patterns,
      caveats = caveats
    )
    quality_overview <- build_data_quality_quality_overview(
      summary_tbl = summary_tbl,
      quality_flags = quality_flags,
      row_review = data.frame(),
      unknown_elements = data.frame(),
      category_counts = category_counts,
      category_usage_summary = category_usage_summary,
      facet_response_patterns = facet_response_patterns,
      caveats = caveats
    )
    out <- list(
      summary = summary_tbl,
      quality_overview = quality_overview,
      quality_flags = quality_flags,
      model_match = model_match,
      row_review = data.frame(),
      unknown_elements = data.frame(),
      category_counts = as.data.frame(category_counts, stringsAsFactors = FALSE),
      score_support_review = score_support_review,
      category_usage_by_facet = category_usage_by_facet,
      category_usage_summary = category_usage_summary,
      facet_response_patterns = facet_response_patterns,
      caveats = caveats,
      score_map = score_map,
      settings = list(
        min_category_count = min_category_count,
        dominant_category_cutoff = dominant_category_cutoff
      )
    )
    if (isTRUE(include_fixed)) {
      out$fixed <- build_sectioned_fixed_report(
        title = "Legacy-compatible Table 2 Data Summary",
        sections = list(
          list(title = "Summary", data = summary_tbl),
          list(title = "Quality overview", data = quality_overview),
          list(title = "Quality flags", data = quality_flags),
          list(title = "Model match", data = model_match),
          list(title = "Row review", data = data.frame()),
          list(title = "Unknown elements", data = data.frame()),
          list(title = "Category counts", data = as.data.frame(category_counts, stringsAsFactors = FALSE)),
          list(title = "Score support review", data = score_support_review),
          list(title = "Facet category usage summary", data = category_usage_summary),
          list(title = "Facet response patterns", data = facet_response_patterns),
          list(title = "Caveats", data = caveats)
        ),
        max_col_width = 24
      )
    }
    return(out)
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }
  required <- c(src$person, src$facets, src$score)
  if (!all(required %in% names(data))) {
    stop("`data` is missing required columns: ", paste(setdiff(required, names(data)), collapse = ", "))
  }
  if (!is.na(src$weight) && !src$weight %in% names(data)) {
    stop("`weight` column not found in `data`.")
  }

  raw <- data
  n_total <- nrow(raw)
  person_ok <- !is.na(raw[[src$person]])
  facet_ok <- if (length(src$facets) == 0) rep(TRUE, n_total) else {
    apply(!is.na(raw[, src$facets, drop = FALSE]), 1, all)
  }
  score_num <- suppressWarnings(as.numeric(raw[[src$score]]))
  score_nonblank <- !is.na(score_num)
  if (!is.na(src$weight)) {
    w_num <- suppressWarnings(as.numeric(raw[[src$weight]]))
    weight_ok <- is.finite(w_num) & w_num > 0
  } else {
    weight_ok <- rep(TRUE, n_total)
  }
  valid_score_labels <- if (nrow(score_map) > 0L && "OriginalScore" %in% names(score_map)) {
    sort(unique(suppressWarnings(as.numeric(score_map$OriginalScore))))
  } else {
    numeric(0)
  }
  valid_score_labels <- valid_score_labels[is.finite(valid_score_labels)]
  range_ok <- if (length(valid_score_labels) > 0L) {
    score_nonblank & score_num %in% valid_score_labels
  } else {
    score_nonblank & score_num >= prep$rating_min & score_num <= prep$rating_max
  }

  row_status <- rep("valid", n_total)
  row_status[!person_ok] <- "missing_person"
  row_status[person_ok & !facet_ok] <- "missing_facet"
  row_status[person_ok & facet_ok & !score_nonblank] <- "missing_score"
  row_status[person_ok & facet_ok & score_nonblank & !weight_ok] <- "invalid_weight"
  row_status[person_ok & facet_ok & score_nonblank & weight_ok & !range_ok] <- "score_out_of_range"

  row_review <- data.frame(Status = row_status, stringsAsFactors = FALSE) |>
    count(.data$Status, name = "N") |>
    arrange(desc(.data$N), .data$Status)

  model_match <- data.frame(
    Model = paste(cfg$model, cfg$method, sep = " / "),
    MatchedResponses = sum(person_ok & facet_ok & score_nonblank, na.rm = TRUE),
    ValidResponses = sum(row_status == "valid", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  unknown_rows <- list()
  check_facets <- c("Person", cfg$facet_names)
  raw_names <- c(src$person, src$facets)
  names(raw_names) <- check_facets
  for (facet in check_facets) {
    col <- raw_names[[facet]]
    if (is.null(col) || !col %in% names(raw)) next
    known <- as.character(prep$levels[[facet]])
    vals <- as.character(raw[[col]])
    bad <- vals[!is.na(vals) & !(vals %in% known)]
    if (length(bad) == 0) next
    unknown_rows[[facet]] <- data.frame(
      Facet = facet,
      Level = sort(unique(bad)),
      stringsAsFactors = FALSE
    )
  }
  unknown_tbl <- if (length(unknown_rows) == 0) data.frame() else bind_rows(unknown_rows)

  summary_tbl <- data.frame(
    TotalLinesInData = n_total,
    TotalDataLines = n_total,
    TotalNonBlankResponsesFound = sum(score_nonblank, na.rm = TRUE),
    MissingScoreRows = sum(row_status == "missing_score", na.rm = TRUE),
    MissingFacetRows = sum(row_status == "missing_facet", na.rm = TRUE),
    MissingPersonRows = sum(row_status == "missing_person", na.rm = TRUE),
    InvalidWeightRows = sum(row_status == "invalid_weight", na.rm = TRUE),
    OutOfRangeScoreRows = sum(row_status == "score_out_of_range", na.rm = TRUE),
    ValidResponsesUsedForEstimation = valid_used,
    ZeroCountScoreCategories = sum(category_counts$ZeroCount %in% TRUE, na.rm = TRUE),
    IntermediateZeroCountScoreCategories = sum(category_counts$UnusedCategoryType == "internal", na.rm = TRUE),
    FacetLevelsWithZeroCategories = sum(category_usage_summary$ZeroCategories > 0, na.rm = TRUE),
    FacetLevelsWithIntermediateZeroCategories = sum(category_usage_summary$IntermediateZeroCategories > 0, na.rm = TRUE),
    FacetLevelsWithSparseCategories = sum(category_usage_summary$SparseCategories > 0, na.rm = TRUE),
    FacetLevelsWithSingleCategoryUse = sum(facet_response_patterns$SingleCategoryUse %in% TRUE, na.rm = TRUE),
    FacetLevelsWithDominantCategoryUse = sum(facet_response_patterns$DominantCategoryUse %in% TRUE, na.rm = TRUE),
    FacetLevelsWithBoundaryOnlyUse = sum(facet_response_patterns$BoundaryOnlyUse %in% TRUE, na.rm = TRUE),
    ScoreSupportCaveats = nrow(caveats),
    stringsAsFactors = FALSE
  )

  quality_flags <- build_data_quality_quality_flags(
    summary_tbl = summary_tbl,
    row_review = row_review,
    unknown_elements = unknown_tbl,
    category_counts = category_counts,
    category_usage_summary = category_usage_summary,
    facet_response_patterns = facet_response_patterns,
    caveats = caveats
  )
  quality_overview <- build_data_quality_quality_overview(
    summary_tbl = summary_tbl,
    quality_flags = quality_flags,
    row_review = row_review,
    unknown_elements = unknown_tbl,
    category_counts = category_counts,
    category_usage_summary = category_usage_summary,
    facet_response_patterns = facet_response_patterns,
    caveats = caveats
  )

  out <- list(
    summary = summary_tbl,
    quality_overview = quality_overview,
    quality_flags = quality_flags,
    model_match = model_match,
    row_review = row_review,
    unknown_elements = unknown_tbl,
    category_counts = as.data.frame(category_counts, stringsAsFactors = FALSE),
    score_support_review = score_support_review,
    category_usage_by_facet = category_usage_by_facet,
    category_usage_summary = category_usage_summary,
    facet_response_patterns = facet_response_patterns,
    caveats = caveats,
    score_map = score_map,
    settings = list(
      min_category_count = min_category_count,
      dominant_category_cutoff = dominant_category_cutoff
    )
  )
  if (isTRUE(include_fixed)) {
    out$fixed <- build_sectioned_fixed_report(
      title = "Legacy-compatible Table 2 Data Summary",
      sections = list(
        list(title = "Summary", data = summary_tbl),
        list(title = "Quality overview", data = quality_overview),
        list(title = "Quality flags", data = quality_flags),
        list(title = "Model match", data = model_match),
        list(title = "Row review", data = row_review),
        list(title = "Unknown elements", data = unknown_tbl),
        list(title = "Category counts", data = as.data.frame(category_counts, stringsAsFactors = FALSE)),
        list(title = "Score support review", data = score_support_review),
        list(title = "Facet category usage summary", data = category_usage_summary),
        list(title = "Facet response patterns", data = facet_response_patterns),
        list(title = "Caveats", data = caveats)
      ),
      max_col_width = 24
    )
  }
  out
}

#' Build a legacy-compatible Table 3 iteration report
#'
#' @param fit Output from [fit_mfrm()].
#' @param max_iter Maximum replay iterations (excluding optional initial row).
#' @param reltol Stopping tolerance for replayed max-logit change.
#' @param include_prox If `TRUE`, include an initial pseudo-row labeled `PROX`.
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#'
#' @details
#' The legacy-compatible Table 3 layout prints per-iteration score residual and
#' logit-change diagnostics.
#' The underlying optimizer in this package does not expose the exact internal
#' per-iteration path, so this function reconstructs an approximation by
#' repeatedly running one-iteration updates from the current parameter vector.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [estimation_iteration_report()].
#'
#' @section Output columns:
#' The `table` data.frame contains:
#' \describe{
#'   \item{Method}{Estimation method label.}
#'   \item{Iteration}{Iteration number.}
#'   \item{MaxLogitChangeElements}{Largest parameter change across
#'     elements in this iteration.}
#'   \item{MaxLogitChangeSteps}{Largest step-parameter change.}
#'   \item{Objective}{Objective function value (negative log-likelihood).}
#' }
#'
#' The `summary` data.frame contains:
#' \describe{
#'   \item{FinalConverged}{Logical; `TRUE` if the optimizer converged.}
#'   \item{FinalIterations}{Number of iterations used.}
#' }
#'
#' @return A named list with:
#' - `table`: iteration rows (residual and logit-change metrics)
#' - `summary`: one-row convergence summary
#' - `settings`: replay settings used
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#'
#' @seealso [fit_mfrm()], [specifications_report()], [data_quality_report()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' t3 <- estimation_iteration_report(fit, max_iter = 5)
#' @keywords internal
#' @noRd
table3_iteration_report <- function(fit,
                                    max_iter = 20,
                                    reltol = NULL,
                                    include_prox = TRUE,
                                    include_fixed = FALSE) {
  signal_legacy_name_deprecation(
    old_name = "table3_iteration_report",
    new_name = "estimation_iteration_report",
    suppress_if_called_from = "estimation_iteration_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  stop_if_gpcm_out_of_scope(fit, "estimation_iteration_report()")
  cfg <- fit$config
  prep <- fit$prep
  sizes <- build_param_sizes(cfg)
  idx <- build_indices(prep, step_facet = cfg$step_facet,
                       slope_facet = cfg$slope_facet,
                       interaction_specs = cfg$interaction_specs)
  est_ctl <- cfg$estimation_control %||% list()
  if (is.null(reltol) || !is.finite(reltol)) reltol <- as.numeric(est_ctl$reltol %||% 1e-6)
  quad_points <- as.integer(est_ctl$quad_points %||% 15L)
  max_iter <- max(1L, as.integer(max_iter))

  current_par <- build_initial_param_vector(cfg, sizes)
  prev_state <- compute_iteration_state(
    par = current_par,
    idx = idx,
    prep = prep,
    config = cfg,
    sizes = sizes,
    quad_points = quad_points
  )

  rows <- list()
  row_id <- 1L
  if (isTRUE(include_prox)) {
    rows[[row_id]] <- tibble(
      Method = "PROX",
      Iteration = 1L,
      MaxScoreResidualElements = prev_state$max_score_resid_elements,
      MaxScoreResidualPercent = prev_state$max_score_resid_pct,
      MaxScoreResidualCategories = prev_state$max_score_resid_categories,
      MaxLogitChangeElements = NA_real_,
      MaxLogitChangeSteps = NA_real_,
      Objective = NA_real_
    )
    row_id <- row_id + 1L
  }

  for (it in seq_len(max_iter)) {
    opt_step <- run_mfrm_optimization(
      start = current_par,
      method = cfg$method,
      idx = idx,
      config = cfg,
      sizes = sizes,
      quad_points = quad_points,
      maxit = 1L,
      reltol = reltol,
      suppress_convergence_warning = TRUE
    )
    current_par <- opt_step$par
    state <- compute_iteration_state(
      par = current_par,
      idx = idx,
      prep = prep,
      config = cfg,
      sizes = sizes,
      quad_points = quad_points
    )

    elem_change <- if (length(state$element_vec) == length(prev_state$element_vec) && length(state$element_vec) > 0) {
      max(abs(state$element_vec - prev_state$element_vec), na.rm = TRUE)
    } else {
      NA_real_
    }
    step_change <- if (length(state$step_vec) == length(prev_state$step_vec) && length(state$step_vec) > 0) {
      max(abs(state$step_vec - prev_state$step_vec), na.rm = TRUE)
    } else {
      NA_real_
    }

    rows[[row_id]] <- tibble(
      Method = cfg$method,
      Iteration = if (isTRUE(include_prox)) it + 1L else it,
      MaxScoreResidualElements = state$max_score_resid_elements,
      MaxScoreResidualPercent = state$max_score_resid_pct,
      MaxScoreResidualCategories = state$max_score_resid_categories,
      MaxLogitChangeElements = elem_change,
      MaxLogitChangeSteps = step_change,
      Objective = -opt_step$value
    )
    row_id <- row_id + 1L
    prev_state <- state

    change_vec <- c(elem_change, step_change)
    change_vec <- change_vec[is.finite(change_vec)]
    max_change <- if (length(change_vec) == 0) NA_real_ else max(change_vec)
    if (is.finite(max_change) && max_change < reltol) break
  }

  tbl <- bind_rows(rows)
  subset_tbl <- calc_subsets(compute_obs_table(fit), c("Person", cfg$facet_names))$summary
  connected <- if (!is.null(subset_tbl) && nrow(subset_tbl) > 0) nrow(subset_tbl) == 1 else NA

  summary_tbl <- data.frame(
    FinalConverged = isTRUE(fit$summary$Converged[1]),
    FinalIterations = as.integer(fit$summary$Iterations[1]),
    ReplayRows = nrow(tbl),
    ConnectedSubset = connected,
    stringsAsFactors = FALSE
  )

  out <- list(
    table = as.data.frame(tbl, stringsAsFactors = FALSE),
    summary = summary_tbl,
    settings = list(
      max_iter = max_iter,
      reltol = reltol,
      include_prox = include_prox,
      quad_points = quad_points,
      include_fixed = isTRUE(include_fixed)
    )
  )
  if (isTRUE(include_fixed)) {
    out$fixed <- build_sectioned_fixed_report(
      title = "Legacy-compatible Table 3 Iteration Report",
      sections = list(
        list(title = "Iteration rows", data = as.data.frame(tbl, stringsAsFactors = FALSE), max_rows = 200L),
        list(title = "Summary", data = summary_tbl),
        list(title = "Settings", data = as.data.frame(out$settings, stringsAsFactors = FALSE))
      ),
      max_col_width = 24
    )
  }
  out
}

#' Build a legacy-compatible Table 6.0.0 subset/disjoint-element listing
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param top_n_subsets Optional maximum number of subset rows to keep.
#' @param min_observations Minimum observations required to keep a subset row.
#'
#' @details
#' The legacy-compatible Table 6.0.0 layout reports disjoint subsets when the
#' design is not fully connected. This helper exposes the same design-check
#' idea with:
#' 1) subset-level counts and percentages, and
#' 2) per-subset/per-facet element listings.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [subset_connectivity_report()].
#'
#' @section Output columns:
#' The `summary` data.frame contains:
#' \describe{
#'   \item{Subset}{Subset index.}
#'   \item{Observations}{Number of observations in this subset.}
#'   \item{ObservationPercent}{Percentage of total observations.}
#' }
#'
#' The `listing` data.frame contains:
#' \describe{
#'   \item{Subset}{Subset index.}
#'   \item{Facet}{Facet name.}
#'   \item{LevelsN}{Number of levels in this facet within this subset.}
#'   \item{Levels}{Comma-separated level labels.}
#'   \item{Ruler}{ASCII ruler for visual comparison.}
#' }
#'
#' @return A named list with:
#' - `summary`: subset-level counts (including `ObservationPercent`)
#' - `listing`: facet-level element listing by subset
#' - `nodes`: node-level table from subset detection
#' - `settings`: applied filters and flags
#'
#' @seealso [diagnose_mfrm()], [measurable_summary_table()], [data_quality_report()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' t6 <- subset_connectivity_report(fit)
#' @keywords internal
#' @noRd
build_subset_connectivity_edges <- function(fit, nodes_tbl, keep_subsets) {
  data <- as.data.frame(fit$prep$data %||% data.frame(), stringsAsFactors = FALSE)
  source_cols <- fit$config$source_columns %||% list()
  person_col <- as.character(source_cols$person %||% "Person")
  facet_cols <- as.character(source_cols$facets %||% fit$config$facet_names %||% character(0))
  cols <- c(person_col, facet_cols)
  cols <- cols[nzchar(cols)]
  if (!is.data.frame(data) || nrow(data) == 0L || length(cols) < 2L ||
      !all(cols %in% names(data)) || !is.data.frame(nodes_tbl) ||
      nrow(nodes_tbl) == 0L || !"Node" %in% names(nodes_tbl)) {
    return(data.frame())
  }
  node_lookup <- nodes_tbl[, intersect(c("Node", "Subset", "Facet", "Level"), names(nodes_tbl)), drop = FALSE]
  node_lookup$Node <- as.character(node_lookup$Node)
  if ("Subset" %in% names(node_lookup)) {
    node_lookup$Subset <- suppressWarnings(as.integer(node_lookup$Subset))
    node_lookup <- node_lookup[node_lookup$Subset %in% keep_subsets, , drop = FALSE]
  }
  if (nrow(node_lookup) == 0L) return(data.frame())

  long <- do.call(rbind, lapply(cols, function(col) {
    data.frame(
      Row = seq_len(nrow(data)),
      Facet = col,
      Level = as.character(data[[col]]),
      stringsAsFactors = FALSE
    )
  }))
  long <- long[!is.na(long$Level) & nzchar(long$Level), , drop = FALSE]
  long$Node <- paste(long$Facet, long$Level, sep = ":")
  long <- merge(
    long,
    node_lookup[, c("Node", "Subset"), drop = FALSE],
    by = "Node",
    all.x = FALSE,
    all.y = FALSE
  )
  if (nrow(long) == 0L) return(data.frame())
  split_rows <- split(long, long$Row)
  edge_list <- lapply(split_rows, function(row_tbl) {
    row_tbl <- row_tbl[!duplicated(row_tbl$Node), , drop = FALSE]
    if (nrow(row_tbl) < 2L) return(NULL)
    cmb <- utils::combn(seq_len(nrow(row_tbl)), 2L)
    from <- row_tbl[cmb[1, ], , drop = FALSE]
    to <- row_tbl[cmb[2, ], , drop = FALSE]
    data.frame(
      Subset = from$Subset,
      From = from$Node,
      To = to$Node,
      FromFacet = from$Facet,
      FromLevel = from$Level,
      ToFacet = to$Facet,
      ToLevel = to$Level,
      stringsAsFactors = FALSE
    )
  })
  edges <- do.call(rbind, edge_list[!vapply(edge_list, is.null, logical(1))])
  if (!is.data.frame(edges) || nrow(edges) == 0L) return(data.frame())
  edges |>
    dplyr::group_by(
      .data$Subset, .data$From, .data$To, .data$FromFacet, .data$FromLevel,
      .data$ToFacet, .data$ToLevel
    ) |>
    dplyr::summarise(Weight = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(.data$Subset, dplyr::desc(.data$Weight), .data$From, .data$To) |>
    as.data.frame(stringsAsFactors = FALSE)
}

table6_subsets_listing <- function(fit,
                                   diagnostics = NULL,
                                   top_n_subsets = NULL,
                                   min_observations = 0) {
  signal_legacy_name_deprecation(
    old_name = "table6_subsets_listing",
    new_name = "subset_connectivity_report",
    suppress_if_called_from = "subset_connectivity_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }

  summary_tbl <- if (!is.null(diagnostics$subsets$summary)) {
    as.data.frame(diagnostics$subsets$summary, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  nodes_tbl <- if (!is.null(diagnostics$subsets$nodes)) {
    as.data.frame(diagnostics$subsets$nodes, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }

  if (nrow(summary_tbl) == 0 || nrow(nodes_tbl) == 0) {
    return(list(
      summary = summary_tbl,
      listing = data.frame(),
      nodes = nodes_tbl,
      edges = data.frame(),
      settings = list(
        top_n_subsets = if (is.null(top_n_subsets)) NA_integer_ else max(1L, as.integer(top_n_subsets)),
        min_observations = as.numeric(min_observations),
        is_disjoint = FALSE
      )
    ))
  }

  summary_tbl$Subset <- suppressWarnings(as.integer(summary_tbl$Subset))
  summary_tbl$Observations <- suppressWarnings(as.numeric(summary_tbl$Observations))
  total_obs <- sum(summary_tbl$Observations, na.rm = TRUE)
  summary_tbl$ObservationPercent <- ifelse(
    is.finite(total_obs) && total_obs > 0,
    100 * summary_tbl$Observations / total_obs,
    NA_real_
  )
  summary_tbl <- summary_tbl |>
    dplyr::arrange(dplyr::desc(.data$Observations), .data$Subset)

  min_observations <- as.numeric(min_observations)
  if (is.finite(min_observations) && min_observations > 0) {
    summary_tbl <- summary_tbl |>
      dplyr::filter(.data$Observations >= min_observations)
  }
  if (!is.null(top_n_subsets)) {
    summary_tbl <- summary_tbl |>
      dplyr::slice_head(n = max(1L, as.integer(top_n_subsets)))
  }

  keep_subsets <- unique(summary_tbl$Subset)
  nodes_tbl <- nodes_tbl |>
    dplyr::mutate(Subset = suppressWarnings(as.integer(.data$Subset))) |>
    dplyr::filter(.data$Subset %in% keep_subsets)

  listing_tbl <- if (nrow(nodes_tbl) == 0) {
    data.frame()
  } else {
    nodes_tbl |>
      dplyr::group_by(.data$Subset, .data$Facet) |>
      dplyr::summarise(
        LevelsN = dplyr::n_distinct(.data$Level),
        Levels = paste(sort(unique(as.character(.data$Level))), collapse = ", "),
        .groups = "drop"
      ) |>
      dplyr::left_join(
        summary_tbl |>
          dplyr::select("Subset", "Observations", "ObservationPercent"),
        by = "Subset"
      ) |>
      dplyr::group_by(.data$Facet) |>
      dplyr::mutate(MaxLevelsN = max(.data$LevelsN, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        Ruler = purrr::map2_chr(
          .data$LevelsN,
          .data$MaxLevelsN,
          function(levels_n, max_levels_n) {
            width <- 20L
            if (!is.finite(max_levels_n) || max_levels_n <= 0 || !is.finite(levels_n)) {
              return("")
            }
            fill <- as.integer(round(width * levels_n / max_levels_n))
            fill <- min(width, max(1L, fill))
            paste0("[", strrep("=", fill), strrep(".", width - fill), "]")
          }
        )
      ) |>
      dplyr::select(-"MaxLevelsN") |>
      dplyr::arrange(.data$Subset, .data$Facet)
  }

  list(
    summary = as.data.frame(summary_tbl, stringsAsFactors = FALSE),
    listing = as.data.frame(listing_tbl, stringsAsFactors = FALSE),
    nodes = as.data.frame(nodes_tbl, stringsAsFactors = FALSE),
    edges = build_subset_connectivity_edges(fit, as.data.frame(nodes_tbl, stringsAsFactors = FALSE), keep_subsets),
    settings = list(
      top_n_subsets = if (is.null(top_n_subsets)) NA_integer_ else max(1L, as.integer(top_n_subsets)),
      min_observations = min_observations,
      is_disjoint = nrow(summary_tbl) > 1
    )
  )
}

#' Build a legacy-compatible Table 6.2 facet-statistics graphic summary
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param metrics Numeric columns in `diagnostics$measures` to summarize.
#' @param ruler_width Width of the fixed-width ruler used for `M/S/Q/X` marks.
#'
#' @details
#' The legacy-compatible Table 6.2 layout describes each facet with a compact
#' graphical ruler where:
#' - `M` marks the mean
#' - `S` marks one SD from the mean
#' - `Q` marks two SD from the mean
#' - `X` marks three SD from the mean
#'
#' Rulers for a given metric share the same min/max scale across facets.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [facet_statistics_report()].
#'
#' @section Output columns:
#' The `table` data.frame contains:
#' \describe{
#'   \item{Metric}{Statistic name (e.g., "Measure", "Infit", "Outfit").}
#'   \item{Facet}{Facet name.}
#'   \item{Levels}{Number of levels.}
#'   \item{Mean, SD, Min, Max}{Summary statistics for this metric-facet
#'     combination.}
#'   \item{Ruler}{ASCII ruler string showing the distribution of level
#'     values. Markers: M = mean, S = +/- 1 SD, Q = +/- 2 SD.}
#' }
#'
#' @return A named list with:
#' - `table`: facet-by-metric summary with ruler strings
#' - `ranges`: global metric ranges used to draw rulers
#' - `settings`: applied metrics and ruler settings
#'
#' @seealso [diagnose_mfrm()], [summary.mfrm_fit()], [plot_facets_chisq()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' t62 <- facet_statistics_report(fit)
#' @keywords internal
#' @noRd
table6_2_facet_statistics <- function(fit,
                                      diagnostics = NULL,
                                      metrics = c("Estimate", "SE", "Infit", "Outfit"),
                                      ruler_width = 41) {
  signal_legacy_name_deprecation(
    old_name = "table6_2_facet_statistics",
    new_name = "facet_statistics_report",
    suppress_if_called_from = "facet_statistics_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (is.null(diagnostics$measures) || nrow(diagnostics$measures) == 0) {
    stop("`diagnostics$measures` is empty. Run diagnose_mfrm() first.")
  }

  measure_tbl <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
  if (!"Facet" %in% names(measure_tbl)) {
    stop("`diagnostics$measures` must include a `Facet` column.")
  }

  metrics <- unique(as.character(metrics))
  available_metrics <- metrics[metrics %in% names(measure_tbl)]
  if (length(available_metrics) == 0) {
    stop("None of `metrics` were found in `diagnostics$measures`.")
  }
  numeric_metrics <- available_metrics[
    vapply(measure_tbl[available_metrics], is.numeric, logical(1))
  ]
  if (length(numeric_metrics) == 0) {
    stop("Selected `metrics` must be numeric columns in `diagnostics$measures`.")
  }

  ruler_width <- max(21L, as.integer(ruler_width))

  metric_ranges <- purrr::map_dfr(numeric_metrics, function(metric) {
    vals <- suppressWarnings(as.numeric(measure_tbl[[metric]]))
    tibble::tibble(
      Metric = metric,
      GlobalMin = if (any(is.finite(vals))) min(vals, na.rm = TRUE) else NA_real_,
      GlobalMax = if (any(is.finite(vals))) max(vals, na.rm = TRUE) else NA_real_
    )
  })

  make_ruler <- function(mean_val, sd_val, global_min, global_max, width) {
    chars <- rep(".", width)
    priority <- c("." = 0L, "X" = 1L, "Q" = 2L, "S" = 3L, "M" = 4L)

    scale_pos <- function(x) {
      if (!is.finite(x) || !is.finite(global_min) || !is.finite(global_max)) {
        return(NA_integer_)
      }
      if (global_max <= global_min) {
        return(as.integer((width + 1L) / 2L))
      }
      pos <- 1L + as.integer(round((width - 1L) * (x - global_min) / (global_max - global_min)))
      as.integer(min(width, max(1L, pos)))
    }

    place_marker <- function(x, marker) {
      pos <- scale_pos(x)
      if (!is.finite(pos)) return(invisible(NULL))
      current <- chars[pos]
      if (priority[[marker]] >= priority[[current]]) {
        chars[pos] <<- marker
      }
      invisible(NULL)
    }

    if (is.finite(mean_val) && is.finite(sd_val) && sd_val >= 0) {
      for (k in c(-3, 3)) place_marker(mean_val + k * sd_val, "X")
      for (k in c(-2, 2)) place_marker(mean_val + k * sd_val, "Q")
      for (k in c(-1, 1)) place_marker(mean_val + k * sd_val, "S")
      place_marker(mean_val, "M")
    } else if (is.finite(mean_val)) {
      place_marker(mean_val, "M")
    }

    paste0("[", paste(chars, collapse = ""), "]")
  }

  out_tbl <- purrr::map_dfr(numeric_metrics, function(metric) {
    gmin <- metric_ranges$GlobalMin[match(metric, metric_ranges$Metric)]
    gmax <- metric_ranges$GlobalMax[match(metric, metric_ranges$Metric)]

    measure_tbl |>
      dplyr::group_by(.data$Facet) |>
      dplyr::summarise(
        Levels = dplyr::n(),
        Mean = mean(.data[[metric]], na.rm = TRUE),
        SD = stats::sd(.data[[metric]], na.rm = TRUE),
        Min = min(.data[[metric]], na.rm = TRUE),
        Max = max(.data[[metric]], na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        Metric = metric,
        GlobalMin = gmin,
        GlobalMax = gmax,
        Ruler = purrr::pmap_chr(
          list(.data$Mean, .data$SD, .data$GlobalMin, .data$GlobalMax),
          function(mean_val, sd_val, global_min, global_max) {
            make_ruler(
              mean_val = mean_val,
              sd_val = sd_val,
              global_min = global_min,
              global_max = global_max,
              width = ruler_width
            )
          }
        ),
        .before = 1
      )
  }) |>
    dplyr::arrange(.data$Metric, .data$Facet)

  list(
    table = as.data.frame(out_tbl, stringsAsFactors = FALSE),
    ranges = as.data.frame(metric_ranges, stringsAsFactors = FALSE),
    settings = list(
      metrics = numeric_metrics,
      ruler_width = ruler_width,
      marker_legend = c(M = "mean", S = "+/-1 SD", Q = "+/-2 SD", X = "+/-3 SD")
    )
  )
}

closest_theta_for_target <- function(theta, y, target) {
  if (length(theta) == 0 || length(y) == 0 || !is.finite(target)) return(NA_real_)
  ok <- is.finite(theta) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  theta <- theta[ok]
  y <- y[ok]
  theta[which.min(abs(y - target))]
}

has_marginal_fit_bundle <- function(diagnostics) {
  is.list(diagnostics) &&
    is.list(diagnostics$marginal_fit) &&
    isTRUE(diagnostics$marginal_fit$available)
}

validate_data_quality_min_category_count <- function(min_category_count) {
  min_category_count <- suppressWarnings(as.numeric(min_category_count[1]))
  if (!is.finite(min_category_count) || min_category_count < 0) {
    stop("`min_category_count` must be a non-negative finite number.", call. = FALSE)
  }
  min_category_count
}

validate_data_quality_dominant_category_cutoff <- function(dominant_category_cutoff) {
  dominant_category_cutoff <- suppressWarnings(as.numeric(dominant_category_cutoff[1]))
  if (!is.finite(dominant_category_cutoff) ||
      dominant_category_cutoff <= 0 ||
      dominant_category_cutoff > 1) {
    stop("`dominant_category_cutoff` must be a finite number in (0, 1].", call. = FALSE)
  }
  dominant_category_cutoff
}

score_category_support_profile <- function(prep = NULL, score_distribution = NULL) {
  support <- integer(0)
  observed <- integer(0)
  unused <- integer(0)

  if (!is.null(prep)) {
    rating_min <- suppressWarnings(as.integer(prep$rating_min %||% NA_integer_))
    rating_max <- suppressWarnings(as.integer(prep$rating_max %||% NA_integer_))
    if (is.finite(rating_min) && is.finite(rating_max) && rating_min <= rating_max) {
      support <- seq(rating_min, rating_max)
    }
    prep_data <- as.data.frame(prep$data %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(prep_data) > 0 && "Score" %in% names(prep_data)) {
      observed <- sort(unique(suppressWarnings(as.integer(prep_data$Score))))
      observed <- observed[is.finite(observed)]
    }
    unused <- sort(unique(suppressWarnings(as.integer(prep$unused_score_categories %||% integer(0)))))
    unused <- unused[is.finite(unused)]
  }

  score_distribution <- as.data.frame(score_distribution %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(score_distribution) > 0) {
    category_col <- if ("Category" %in% names(score_distribution)) {
      "Category"
    } else if ("Score" %in% names(score_distribution)) {
      "Score"
    } else {
      NA_character_
    }
    if (!is.na(category_col)) {
      dist_categories <- suppressWarnings(as.integer(score_distribution[[category_col]]))
      dist_categories <- dist_categories[is.finite(dist_categories)]
      support <- sort(unique(c(support, dist_categories)))

      raw_n <- suppressWarnings(as.numeric(score_distribution$RawN %||% score_distribution$Count %||% NA_real_))
      weighted_n <- suppressWarnings(as.numeric(score_distribution$WeightedN %||% raw_n))
      zero_mask <- (is.finite(raw_n) & raw_n <= 0) | (is.finite(weighted_n) & weighted_n <= 0)
      if (length(zero_mask) == nrow(score_distribution)) {
        zero_categories <- suppressWarnings(as.integer(score_distribution[[category_col]][zero_mask]))
        zero_categories <- zero_categories[is.finite(zero_categories)]
        unused <- sort(unique(c(unused, zero_categories)))

        positive_categories <- suppressWarnings(as.integer(score_distribution[[category_col]][!zero_mask]))
        positive_categories <- positive_categories[is.finite(positive_categories)]
        observed <- sort(unique(c(observed, positive_categories)))
      }
    }
  }

  if (length(support) == 0L) {
    support <- sort(unique(c(observed, unused)))
  }
  support <- support[is.finite(support)]
  if (length(support) == 0L) {
    return(tibble::tibble(
      Category = integer(0),
      ZeroCount = logical(0),
      UnusedCategoryType = character(0),
      WeaklyIdentified = logical(0),
      CategoryCaveat = character(0)
    ))
  }
  if (length(observed) == 0L && length(unused) > 0L) {
    observed <- setdiff(support, unused)
  }

  observed_min <- suppressWarnings(min(observed, na.rm = TRUE))
  observed_max <- suppressWarnings(max(observed, na.rm = TRUE))
  zero_count <- support %in% unused
  internal <- if (is.finite(observed_min) && is.finite(observed_max)) {
    support > observed_min & support < observed_max
  } else {
    rep(FALSE, length(support))
  }
  unused_type <- ifelse(zero_count & internal, "internal", ifelse(zero_count, "boundary", "none"))
  caveat <- dplyr::case_when(
    unused_type == "internal" ~ "Zero-count intermediate category; adjacent thresholds are weakly identified.",
    unused_type == "boundary" ~ "Zero-count boundary category; document the retained support and avoid overinterpreting adjacent thresholds.",
    TRUE ~ ""
  )

  tibble::tibble(
    Category = support,
    ZeroCount = zero_count,
    UnusedCategoryType = unused_type,
    WeaklyIdentified = zero_count,
    CategoryCaveat = caveat
  )
}

build_data_quality_category_counts <- function(prep, observed_counts = NULL) {
  observed_counts <- as.data.frame(observed_counts %||% data.frame(), stringsAsFactors = FALSE)
  rating_min <- suppressWarnings(as.integer(prep$rating_min %||% NA_integer_))
  rating_max <- suppressWarnings(as.integer(prep$rating_max %||% NA_integer_))
  support <- if (is.finite(rating_min) && is.finite(rating_max) && rating_min <= rating_max) {
    seq(rating_min, rating_max)
  } else if (nrow(observed_counts) > 0 && "Score" %in% names(observed_counts)) {
    sort(unique(suppressWarnings(as.integer(observed_counts$Score))))
  } else {
    integer(0)
  }
  support <- support[is.finite(support)]
  if (length(support) == 0L) {
    return(data.frame(
      Score = integer(0),
      Count = integer(0),
      WeightedCount = numeric(0),
      Percent = numeric(0),
      ZeroCount = logical(0),
      UnusedCategoryType = character(0),
      WeaklyIdentified = logical(0),
      SupportRole = character(0),
      CategoryCaveat = character(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- data.frame(Score = support, stringsAsFactors = FALSE)
  if (nrow(observed_counts) > 0 && all(c("Score", "Count") %in% names(observed_counts))) {
    observed_counts$Score <- suppressWarnings(as.integer(observed_counts$Score))
    out <- merge(out, observed_counts, by = "Score", all.x = TRUE, sort = FALSE)
    out <- out[match(support, out$Score), , drop = FALSE]
  } else {
    out$Count <- NA_integer_
    out$WeightedCount <- NA_real_
  }
  if (!"WeightedCount" %in% names(out)) out$WeightedCount <- out$Count
  out$Count <- suppressWarnings(as.integer(out$Count))
  out$WeightedCount <- suppressWarnings(as.numeric(out$WeightedCount))
  out$Count[is.na(out$Count)] <- 0L
  out$WeightedCount[is.na(out$WeightedCount)] <- 0
  total_weight <- sum(out$WeightedCount, na.rm = TRUE)
  out$Percent <- if (total_weight > 0) 100 * out$WeightedCount / total_weight else NA_real_

  profile <- as.data.frame(score_category_support_profile(prep = prep, score_distribution = out), stringsAsFactors = FALSE)
  if (nrow(profile) > 0 && "Category" %in% names(profile)) {
    names(profile)[names(profile) == "Category"] <- "Score"
    out <- merge(out, profile, by = "Score", all.x = TRUE, sort = FALSE)
    out <- out[match(support, out$Score), , drop = FALSE]
  }
  out$ZeroCount <- as.logical(out$ZeroCount %||% (out$Count == 0L))
  out$ZeroCount[is.na(out$ZeroCount)] <- out$Count[is.na(out$ZeroCount)] == 0L
  out$UnusedCategoryType <- as.character(out$UnusedCategoryType %||% "none")
  out$UnusedCategoryType[is.na(out$UnusedCategoryType) | !nzchar(out$UnusedCategoryType)] <- "none"
  out$WeaklyIdentified <- as.logical(out$WeaklyIdentified %||% out$ZeroCount)
  out$WeaklyIdentified[is.na(out$WeaklyIdentified)] <- out$ZeroCount[is.na(out$WeaklyIdentified)]
  out$CategoryCaveat <- as.character(out$CategoryCaveat %||% "")
  out$CategoryCaveat[is.na(out$CategoryCaveat)] <- ""
  out$SupportRole <- ifelse(
    out$ZeroCount & out$UnusedCategoryType == "internal",
    "zero_count_intermediate",
    ifelse(
      out$ZeroCount,
      "zero_count_boundary",
      "observed"
    )
  )
  row.names(out) <- NULL
  out
}

build_data_quality_score_support_review <- function(category_counts) {
  category_counts <- as.data.frame(category_counts %||% data.frame(), stringsAsFactors = FALSE)
  keep <- intersect(
    c("Score", "Count", "WeightedCount", "Percent", "ZeroCount",
      "UnusedCategoryType", "WeaklyIdentified", "SupportRole", "CategoryCaveat"),
    names(category_counts)
  )
  out <- category_counts[, keep, drop = FALSE]
  if (nrow(out) > 0 && "ZeroCount" %in% names(out)) {
    out <- out[order(!as.logical(out$ZeroCount), out$Score), , drop = FALSE]
  }
  row.names(out) <- NULL
  out
}

build_data_quality_facet_category_usage <- function(prep,
                                                    category_counts,
                                                    min_category_count = 10) {
  df <- as.data.frame(prep$data %||% data.frame(), stringsAsFactors = FALSE)
  facet_names <- as.character(prep$facet_names %||% character(0))
  category_counts <- as.data.frame(category_counts %||% data.frame(), stringsAsFactors = FALSE)
  support <- suppressWarnings(as.integer(category_counts$Score %||% integer(0)))
  support <- support[is.finite(support)]
  if (length(support) == 0L) {
    rating_min <- suppressWarnings(as.integer(prep$rating_min %||% NA_integer_))
    rating_max <- suppressWarnings(as.integer(prep$rating_max %||% NA_integer_))
    if (is.finite(rating_min) && is.finite(rating_max) && rating_min <= rating_max) {
      support <- seq(rating_min, rating_max)
    }
  }
  if (nrow(df) == 0L || length(facet_names) == 0L || length(support) == 0L) {
    return(data.frame(
      Facet = character(0),
      Level = character(0),
      Score = integer(0),
      Count = integer(0),
      WeightedCount = numeric(0),
      LevelTotalCount = integer(0),
      LevelTotalWeightedCount = numeric(0),
      PercentWithinLevel = numeric(0),
      ZeroCount = logical(0),
      SparseCount = logical(0),
      CategoryPosition = character(0),
      UsageStatus = character(0),
      ReviewStatus = character(0),
      CategoryCaveat = character(0),
      stringsAsFactors = FALSE
    ))
  }

  min_score <- min(support, na.rm = TRUE)
  max_score <- max(support, na.rm = TRUE)
  chunks <- vector("list", length(facet_names))
  for (i in seq_along(facet_names)) {
    facet <- facet_names[i]
    if (!facet %in% names(df)) next
    levels_f <- as.character(prep$levels[[facet]] %||% sort(unique(as.character(df[[facet]]))))
    levels_f <- levels_f[nzchar(levels_f)]
    if (length(levels_f) == 0L) next
    grid <- expand.grid(
      Level = levels_f,
      Score = support,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    counts <- df |>
      dplyr::group_by(Level = as.character(.data[[facet]]), Score = .data$Score) |>
      dplyr::summarise(
        Count = dplyr::n(),
        WeightedCount = sum(.data$Weight, na.rm = TRUE),
        .groups = "drop"
      ) |>
      as.data.frame(stringsAsFactors = FALSE)
    counts$Score <- suppressWarnings(as.integer(counts$Score))
    out <- merge(grid, counts, by = c("Level", "Score"), all.x = TRUE, sort = FALSE)
    out <- out[order(match(out$Level, levels_f), match(out$Score, support)), , drop = FALSE]
    out$Count <- suppressWarnings(as.integer(out$Count))
    out$WeightedCount <- suppressWarnings(as.numeric(out$WeightedCount))
    out$Count[is.na(out$Count)] <- 0L
    out$WeightedCount[is.na(out$WeightedCount)] <- 0
    totals <- out |>
      dplyr::group_by(.data$Level) |>
      dplyr::summarise(
        LevelTotalCount = sum(.data$Count, na.rm = TRUE),
        LevelTotalWeightedCount = sum(.data$WeightedCount, na.rm = TRUE),
        .groups = "drop"
      ) |>
      as.data.frame(stringsAsFactors = FALSE)
    out <- merge(out, totals, by = "Level", all.x = TRUE, sort = FALSE)
    out <- out[order(match(out$Level, levels_f), match(out$Score, support)), , drop = FALSE]
    out$PercentWithinLevel <- ifelse(
      out$LevelTotalWeightedCount > 0,
      100 * out$WeightedCount / out$LevelTotalWeightedCount,
      NA_real_
    )
    out$ZeroCount <- out$Count <= 0L
    out$SparseCount <- out$Count > 0L & (
      out$Count < min_category_count | out$WeightedCount < min_category_count
    )
    out$CategoryPosition <- ifelse(
      out$Score > min_score & out$Score < max_score,
      "intermediate",
      "boundary"
    )
    out$UsageStatus <- ifelse(
      out$ZeroCount,
      "zero",
      ifelse(out$SparseCount, "sparse", "adequate")
    )
    out$ReviewStatus <- ifelse(
      out$ZeroCount & out$CategoryPosition == "intermediate",
      "warning",
      ifelse(out$ZeroCount | out$SparseCount, "review", "ok")
    )
    out$CategoryCaveat <- dplyr::case_when(
      out$ZeroCount & out$CategoryPosition == "intermediate" ~
        "This facet level has zero observations in an intermediate score category; adjacent thresholds may be weakly supported for this level.",
      out$ZeroCount ~
        "This facet level has zero observations in a boundary score category; document the limited category use.",
      out$SparseCount ~
        paste0("This facet level has fewer than ", min_category_count,
               " observations in this score category; treat category-functioning evidence as sparse."),
      TRUE ~ ""
    )
    out$Facet <- facet
    out <- out[, c(
      "Facet", "Level", "Score", "Count", "WeightedCount",
      "LevelTotalCount", "LevelTotalWeightedCount", "PercentWithinLevel",
      "ZeroCount", "SparseCount", "CategoryPosition", "UsageStatus",
      "ReviewStatus", "CategoryCaveat"
    ), drop = FALSE]
    chunks[[i]] <- out
  }
  out <- do.call(rbind, chunks[!vapply(chunks, is.null, logical(1))])
  if (is.null(out)) out <- data.frame()
  row.names(out) <- NULL
  out
}

build_data_quality_facet_category_summary <- function(category_usage_by_facet) {
  tbl <- as.data.frame(category_usage_by_facet %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L || !all(c("Facet", "Level", "Score") %in% names(tbl))) {
    return(data.frame(
      Facet = character(0),
      Level = character(0),
      Categories = integer(0),
      ObservedCategories = integer(0),
      ZeroCategories = integer(0),
      IntermediateZeroCategories = integer(0),
      SparseCategories = integer(0),
      IssueCategories = integer(0),
      CategoryCoverageRate = numeric(0),
      MinCount = integer(0),
      MinWeightedCount = numeric(0),
      ReviewStatus = character(0),
      stringsAsFactors = FALSE
    ))
  }
  out <- tbl |>
    dplyr::group_by(.data$Facet, .data$Level) |>
    dplyr::summarise(
      Categories = dplyr::n_distinct(.data$Score),
      ObservedCategories = sum(!.data$ZeroCount, na.rm = TRUE),
      ZeroCategories = sum(.data$ZeroCount, na.rm = TRUE),
      IntermediateZeroCategories = sum(.data$ZeroCount & .data$CategoryPosition == "intermediate", na.rm = TRUE),
      SparseCategories = sum(.data$SparseCount, na.rm = TRUE),
      IssueCategories = sum(.data$ZeroCount | .data$SparseCount, na.rm = TRUE),
      MinCount = suppressWarnings(min(.data$Count, na.rm = TRUE)),
      MinWeightedCount = suppressWarnings(min(.data$WeightedCount, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      CategoryCoverageRate = ifelse(.data$Categories > 0, .data$ObservedCategories / .data$Categories, NA_real_),
      ReviewStatus = dplyr::case_when(
        .data$IntermediateZeroCategories > 0 ~ "warning",
        .data$ZeroCategories > 0 | .data$SparseCategories > 0 ~ "review",
        TRUE ~ "ok"
      )
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$IntermediateZeroCategories),
      dplyr::desc(.data$ZeroCategories),
      dplyr::desc(.data$SparseCategories),
      .data$Facet,
      .data$Level
    ) |>
    as.data.frame(stringsAsFactors = FALSE)
  row.names(out) <- NULL
  out
}

empty_data_quality_facet_response_patterns <- function() {
  data.frame(
    Facet = character(0),
    Level = character(0),
    Responses = integer(0),
    WeightedResponses = numeric(0),
    Categories = integer(0),
    ObservedCategories = integer(0),
    CategoryCoverageRate = numeric(0),
    DominantScore = integer(0),
    DominantCount = integer(0),
    DominantWeightedCount = numeric(0),
    DominantPercent = numeric(0),
    SingleCategoryUse = logical(0),
    DominantCategoryUse = logical(0),
    BoundaryOnlyUse = logical(0),
    PatternStatus = character(0),
    PatternLabel = character(0),
    PatternCaveat = character(0),
    stringsAsFactors = FALSE
  )
}

build_data_quality_facet_response_patterns <- function(category_usage_by_facet,
                                                       dominant_category_cutoff = 0.95) {
  tbl <- as.data.frame(category_usage_by_facet %||% data.frame(), stringsAsFactors = FALSE)
  dominant_category_cutoff <- suppressWarnings(as.numeric(dominant_category_cutoff[1]))
  if (!is.finite(dominant_category_cutoff) ||
      dominant_category_cutoff <= 0 ||
      dominant_category_cutoff > 1) {
    dominant_category_cutoff <- 0.95
  }
  required <- c("Facet", "Level", "Score", "Count", "WeightedCount", "CategoryPosition")
  if (nrow(tbl) == 0L || !all(required %in% names(tbl))) {
    return(empty_data_quality_facet_response_patterns())
  }
  tbl$Count <- suppressWarnings(as.numeric(tbl$Count))
  tbl$WeightedCount <- suppressWarnings(as.numeric(tbl$WeightedCount))
  tbl$ScoreNumeric <- suppressWarnings(as.numeric(tbl$Score))
  chunks <- split(tbl, paste(tbl$Facet, tbl$Level, sep = "\r"))
  rows <- lapply(chunks, function(x) {
    x <- x[order(x$ScoreNumeric, as.character(x$Score), na.last = TRUE), , drop = FALSE]
    responses <- sum(x$Count, na.rm = TRUE)
    weighted_responses <- sum(x$WeightedCount, na.rm = TRUE)
    observed <- x[x$Count > 0, , drop = FALSE]
    observed_categories <- nrow(observed)
    dominant_idx <- if (nrow(observed) > 0L) {
      which.max(observed$WeightedCount)
    } else {
      NA_integer_
    }
    dominant_score <- if (is.finite(dominant_idx)) observed$Score[dominant_idx] else NA
    dominant_count <- if (is.finite(dominant_idx)) observed$Count[dominant_idx] else NA_real_
    dominant_weighted <- if (is.finite(dominant_idx)) observed$WeightedCount[dominant_idx] else NA_real_
    dominant_percent <- if (is.finite(weighted_responses) && weighted_responses > 0) {
      dominant_weighted / weighted_responses
    } else if (is.finite(responses) && responses > 0) {
      dominant_count / responses
    } else {
      NA_real_
    }
    intermediate_used <- any(
      observed$CategoryPosition == "intermediate" & observed$Count > 0,
      na.rm = TRUE
    )
    single_category <- is.finite(responses) && responses > 0 && observed_categories <= 1L
    dominant_category <- is.finite(dominant_percent) &&
      dominant_percent >= dominant_category_cutoff &&
      !single_category
    boundary_only <- is.finite(responses) && responses > 0 &&
      observed_categories > 0L &&
      !intermediate_used &&
      any(x$CategoryPosition == "intermediate", na.rm = TRUE) &&
      !single_category
    status <- if (single_category) {
      "high"
    } else if (dominant_category || boundary_only) {
      "review"
    } else {
      "ok"
    }
    label <- if (single_category) {
      "single_category"
    } else if (dominant_category) {
      "dominant_category"
    } else if (boundary_only) {
      "boundary_only"
    } else {
      "mixed_category_use"
    }
    caveat <- dplyr::case_when(
      single_category ~
        paste0("This facet level used only score ", dominant_score,
               "; local severity/leniency and fit evidence should be reviewed before substantive interpretation."),
      dominant_category ~
        paste0("This facet level assigned score ", dominant_score,
               " to at least ", round(100 * dominant_category_cutoff),
               "% of responses; review for restricted category use."),
      boundary_only ~
        "This facet level used only boundary categories and no intermediate score categories; review local scale use.",
      TRUE ~ ""
    )
    data.frame(
      Facet = as.character(x$Facet[1]),
      Level = as.character(x$Level[1]),
      Responses = as.integer(responses),
      WeightedResponses = weighted_responses,
      Categories = dplyr::n_distinct(x$Score),
      ObservedCategories = as.integer(observed_categories),
      CategoryCoverageRate = ifelse(dplyr::n_distinct(x$Score) > 0,
                                    observed_categories / dplyr::n_distinct(x$Score),
                                    NA_real_),
      DominantScore = suppressWarnings(as.integer(dominant_score)),
      DominantCount = suppressWarnings(as.integer(dominant_count)),
      DominantWeightedCount = dominant_weighted,
      DominantPercent = dominant_percent,
      SingleCategoryUse = single_category,
      DominantCategoryUse = dominant_category,
      BoundaryOnlyUse = boundary_only,
      PatternStatus = status,
      PatternLabel = label,
      PatternCaveat = caveat,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(empty_data_quality_facet_response_patterns())
  out <- out[order(
    match(out$PatternStatus, c("high", "review", "ok")),
    -out$DominantPercent,
    out$Facet,
    out$Level
  ), , drop = FALSE]
  row.names(out) <- NULL
  out
}

empty_data_quality_quality_flags <- function() {
  data.frame(
    Area = character(0),
    Severity = character(0),
    Flag = character(0),
    Count = integer(0),
    Unit = character(0),
    PercentOfData = numeric(0),
    Action = character(0),
    stringsAsFactors = FALSE
  )
}

build_data_quality_quality_flags <- function(summary_tbl,
                                             row_review = NULL,
                                             unknown_elements = NULL,
                                             category_counts = NULL,
                                             category_usage_summary = NULL,
                                             facet_response_patterns = NULL,
                                             caveats = NULL) {
  summary_tbl <- as.data.frame(summary_tbl %||% data.frame(), stringsAsFactors = FALSE)
  row_review <- as.data.frame(row_review %||% data.frame(), stringsAsFactors = FALSE)
  unknown_elements <- as.data.frame(unknown_elements %||% data.frame(), stringsAsFactors = FALSE)
  category_counts <- as.data.frame(category_counts %||% data.frame(), stringsAsFactors = FALSE)
  category_usage_summary <- as.data.frame(category_usage_summary %||% data.frame(), stringsAsFactors = FALSE)
  facet_response_patterns <- as.data.frame(facet_response_patterns %||% data.frame(), stringsAsFactors = FALSE)
  caveats <- as.data.frame(caveats %||% data.frame(), stringsAsFactors = FALSE)

  rows <- list()
  total_data <- NA_real_
  if (nrow(summary_tbl) > 0L) {
    for (nm in c("TotalDataLines", "TotalLinesInData")) {
      if (nm %in% names(summary_tbl)) {
        total_data <- suppressWarnings(as.numeric(summary_tbl[[nm]][1]))
        if (is.finite(total_data) && total_data > 0) break
      }
    }
  }
  percent_of_data <- function(n) {
    if (!is.finite(total_data) || total_data <= 0) return(NA_real_)
    100 * n / total_data
  }
  add_flag <- function(area, severity, flag, count, unit, action, denominator = "data") {
    count <- suppressWarnings(as.numeric(count))
    if (!is.finite(count) || count <= 0) return(invisible(NULL))
    rows[[length(rows) + 1L]] <<- data.frame(
      Area = area,
      Severity = severity,
      Flag = flag,
      Count = as.integer(count),
      Unit = unit,
      PercentOfData = if (identical(denominator, "data")) percent_of_data(count) else NA_real_,
      Action = action,
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  if (nrow(summary_tbl) > 0L) {
    get_summary_count <- function(name) {
      if (!name %in% names(summary_tbl)) return(0)
      value <- suppressWarnings(as.numeric(summary_tbl[[name]][1]))
      ifelse(is.finite(value), value, 0)
    }
    add_flag(
      "Rows", "review", "Rows have missing scores",
      get_summary_count("MissingScoreRows"), "rows",
      "Check missing-code handling or resolve missing scores before final estimation."
    )
    add_flag(
      "Rows", "review", "Rows have missing facet values",
      get_summary_count("MissingFacetRows"), "rows",
      "Check facet identifiers and recode intended missing markers explicitly."
    )
    add_flag(
      "Rows", "review", "Rows have missing person identifiers",
      get_summary_count("MissingPersonRows"), "rows",
      "Resolve person identifiers before fitting or exclude the rows intentionally."
    )
    add_flag(
      "Rows", "review", "Rows have non-positive or invalid weights",
      get_summary_count("InvalidWeightRows"), "rows",
      "Inspect the weight column and use positive finite weights."
    )
    add_flag(
      "Rows", "review", "Rows have scores outside the fitted scale",
      get_summary_count("OutOfRangeScoreRows"), "rows",
      "Check the declared rating scale and the original score labels."
    )
    add_flag(
      "Score support", "high", "Intermediate score categories have zero observations",
      get_summary_count("IntermediateZeroCountScoreCategories"), "categories",
      "Treat adjacent thresholds as weakly supported and document the score-support gap.",
      denominator = "categories"
    )
    add_flag(
      "Score support", "review", "Declared score categories have zero observations",
      get_summary_count("ZeroCountScoreCategories"), "categories",
      "Inspect category use before interpreting threshold or category-functioning results.",
      denominator = "categories"
    )
    add_flag(
      "Facet category use", "high", "Facet levels have intermediate zero-category use",
      get_summary_count("FacetLevelsWithIntermediateZeroCategories"), "facet levels",
      "Review category use by facet level before making rater or criterion conclusions.",
      denominator = "facet_levels"
    )
    add_flag(
      "Facet category use", "review", "Facet levels have zero-category use",
      get_summary_count("FacetLevelsWithZeroCategories"), "facet levels",
      "Inspect the affected facet levels and document limited local category support.",
      denominator = "facet_levels"
    )
    add_flag(
      "Facet category use", "review", "Facet levels have sparse category use",
      get_summary_count("FacetLevelsWithSparseCategories"), "facet levels",
      "Interpret category-functioning evidence for these levels as sparse.",
      denominator = "facet_levels"
    )
  }

  has_row_summary_counts <- nrow(summary_tbl) > 0L &&
    any(c("MissingScoreRows", "MissingFacetRows", "MissingPersonRows",
          "InvalidWeightRows", "OutOfRangeScoreRows") %in% names(summary_tbl))
  if (!has_row_summary_counts &&
      nrow(row_review) > 0L && all(c("Status", "N") %in% names(row_review))) {
    non_valid <- row_review[!tolower(as.character(row_review$Status)) %in% "valid", , drop = FALSE]
    total_non_valid <- sum(suppressWarnings(as.numeric(non_valid$N)), na.rm = TRUE)
    add_flag(
      "Rows", "review", "Rows were excluded before estimation",
      total_non_valid, "rows",
      "Use `row_review` to identify whether exclusion comes from missingness, weights, or scale range."
    )
  }
  if (nrow(unknown_elements) > 0L) {
    add_flag(
      "Design match", "review", "Raw data contain levels absent from the fitted design",
      nrow(unknown_elements), "levels",
      "Check identifiers, facet labels, and whether the model should be refit with the intended design.",
      denominator = "levels"
    )
  }
  has_facet_usage_summary_counts <- nrow(summary_tbl) > 0L &&
    any(c("FacetLevelsWithZeroCategories",
          "FacetLevelsWithIntermediateZeroCategories",
          "FacetLevelsWithSparseCategories") %in% names(summary_tbl))
  if (!has_facet_usage_summary_counts &&
      nrow(category_usage_summary) > 0L && "ReviewStatus" %in% names(category_usage_summary)) {
    warning_levels <- sum(category_usage_summary$ReviewStatus %in% "warning", na.rm = TRUE)
    review_levels <- sum(category_usage_summary$ReviewStatus %in% "review", na.rm = TRUE)
    add_flag(
      "Facet category use", "high", "Facet levels need category-use review",
      warning_levels, "facet levels",
      "Start with levels marked `warning` in `category_usage_summary`.",
      denominator = "facet_levels"
    )
    add_flag(
      "Facet category use", "review", "Facet levels have category-use cautions",
      review_levels, "facet levels",
      "Review zero or sparse category use in `category_usage_by_facet`.",
      denominator = "facet_levels"
    )
  }
  if (nrow(facet_response_patterns) > 0L && "PatternStatus" %in% names(facet_response_patterns)) {
    single_levels <- sum(facet_response_patterns$SingleCategoryUse %in% TRUE, na.rm = TRUE)
    dominant_levels <- sum(facet_response_patterns$DominantCategoryUse %in% TRUE, na.rm = TRUE)
    boundary_only_levels <- sum(facet_response_patterns$BoundaryOnlyUse %in% TRUE, na.rm = TRUE)
    add_flag(
      "Facet response patterns", "high", "Facet levels use only one score category",
      single_levels, "facet levels",
      "Inspect `facet_response_patterns`; for raters, single-category use may indicate unusable or non-discriminating scoring.",
      denominator = "facet_levels"
    )
    add_flag(
      "Facet response patterns", "review", "Facet levels have dominant score-category use",
      dominant_levels, "facet levels",
      "Inspect dominant-category rates before interpreting local fit, severity, or category-functioning evidence.",
      denominator = "facet_levels"
    )
    add_flag(
      "Facet response patterns", "review", "Facet levels use only boundary score categories",
      boundary_only_levels, "facet levels",
      "Inspect whether the facet level skipped intermediate categories in a way that weakens local scale interpretation.",
      denominator = "facet_levels"
    )
  }

  if (nrow(caveats) > 0L && "Condition" %in% names(caveats)) {
    recoded <- sum(caveats$Condition %in% "score_categories_recoded", na.rm = TRUE)
    original_gap <- sum(caveats$Condition %in% "original_score_gap_before_recoding", na.rm = TRUE)
    add_flag(
      "Score support", "review", "Original score labels were recoded",
      recoded, "conditions",
      "Check `score_map` before reporting category or threshold labels.",
      denominator = "conditions"
    )
    add_flag(
      "Score support", "high", "Original score sequence had gaps before recoding",
      original_gap, "conditions",
      "Report the original-label gap and avoid treating the recoded scale as fully observed.",
      denominator = "conditions"
    )
  }

  if (length(rows) == 0L) return(empty_data_quality_quality_flags())
  out <- do.call(rbind, rows)
  severity_order <- c(high = 1L, review = 2L, notice = 3L)
  out$SeverityRank <- unname(severity_order[tolower(as.character(out$Severity))])
  out$SeverityRank[is.na(out$SeverityRank)] <- 99L
  out <- out[order(out$SeverityRank, out$Area, out$Flag), , drop = FALSE]
  out$SeverityRank <- NULL
  row.names(out) <- NULL
  out
}

empty_data_quality_quality_overview <- function() {
  data.frame(
    Area = character(0),
    Status = character(0),
    Count = integer(0),
    Unit = character(0),
    PercentOfData = numeric(0),
    Message = character(0),
    NextStep = character(0),
    stringsAsFactors = FALSE
  )
}

build_data_quality_quality_overview <- function(summary_tbl,
                                                quality_flags = NULL,
                                                row_review = NULL,
                                                unknown_elements = NULL,
                                                category_counts = NULL,
                                                category_usage_summary = NULL,
                                                facet_response_patterns = NULL,
                                                caveats = NULL) {
  summary_tbl <- as.data.frame(summary_tbl %||% data.frame(), stringsAsFactors = FALSE)
  quality_flags <- as.data.frame(quality_flags %||% data.frame(), stringsAsFactors = FALSE)
  row_review <- as.data.frame(row_review %||% data.frame(), stringsAsFactors = FALSE)
  unknown_elements <- as.data.frame(unknown_elements %||% data.frame(), stringsAsFactors = FALSE)
  category_counts <- as.data.frame(category_counts %||% data.frame(), stringsAsFactors = FALSE)
  category_usage_summary <- as.data.frame(category_usage_summary %||% data.frame(), stringsAsFactors = FALSE)
  facet_response_patterns <- as.data.frame(facet_response_patterns %||% data.frame(), stringsAsFactors = FALSE)
  caveats <- as.data.frame(caveats %||% data.frame(), stringsAsFactors = FALSE)

  total_data <- NA_real_
  if (nrow(summary_tbl) > 0L) {
    for (nm in c("TotalDataLines", "TotalLinesInData")) {
      if (!nm %in% names(summary_tbl)) next
      total_data <- suppressWarnings(as.numeric(summary_tbl[[nm]][1]))
      if (is.finite(total_data) && total_data > 0) break
    }
  }
  get_summary_count <- function(name) {
    if (!name %in% names(summary_tbl)) return(0)
    value <- suppressWarnings(as.numeric(summary_tbl[[name]][1]))
    ifelse(is.finite(value), value, 0)
  }
  percent_of_data <- function(n) {
    if (!is.finite(total_data) || total_data <= 0) return(NA_real_)
    100 * n / total_data
  }
  caveat_conditions <- if ("Condition" %in% names(caveats)) {
    as.character(caveats$Condition)
  } else {
    character(0)
  }

  row_issue_count <- sum(
    get_summary_count("MissingScoreRows"),
    get_summary_count("MissingFacetRows"),
    get_summary_count("MissingPersonRows"),
    get_summary_count("InvalidWeightRows"),
    get_summary_count("OutOfRangeScoreRows"),
    na.rm = TRUE
  )
  if (row_issue_count == 0 && nrow(row_review) > 0L && all(c("Status", "N") %in% names(row_review))) {
    non_valid <- row_review[!tolower(as.character(row_review$Status)) %in% "valid", , drop = FALSE]
    row_issue_count <- sum(suppressWarnings(as.numeric(non_valid$N)), na.rm = TRUE)
  }

  zero_categories <- get_summary_count("ZeroCountScoreCategories")
  intermediate_zero_categories <- get_summary_count("IntermediateZeroCountScoreCategories")
  recode_conditions <- sum(caveat_conditions %in% "score_categories_recoded", na.rm = TRUE)
  original_gap_conditions <- sum(caveat_conditions %in% "original_score_gap_before_recoding", na.rm = TRUE)
  score_support_count <- if (zero_categories > 0) zero_categories else recode_conditions + original_gap_conditions

  facet_zero_levels <- get_summary_count("FacetLevelsWithZeroCategories")
  facet_intermediate_zero_levels <- get_summary_count("FacetLevelsWithIntermediateZeroCategories")
  facet_sparse_levels <- get_summary_count("FacetLevelsWithSparseCategories")
  facet_issue_count <- max(facet_zero_levels, facet_intermediate_zero_levels, facet_sparse_levels, na.rm = TRUE)
  if (!is.finite(facet_issue_count)) facet_issue_count <- 0
  if (facet_issue_count == 0 && nrow(category_usage_summary) > 0L && "ReviewStatus" %in% names(category_usage_summary)) {
    facet_issue_count <- sum(!tolower(as.character(category_usage_summary$ReviewStatus)) %in% "ok", na.rm = TRUE)
  }

  single_category_levels <- get_summary_count("FacetLevelsWithSingleCategoryUse")
  dominant_category_levels <- get_summary_count("FacetLevelsWithDominantCategoryUse")
  boundary_only_levels <- get_summary_count("FacetLevelsWithBoundaryOnlyUse")
  response_pattern_count <- max(single_category_levels, dominant_category_levels, boundary_only_levels, na.rm = TRUE)
  if (!is.finite(response_pattern_count)) response_pattern_count <- 0
  if (response_pattern_count == 0 &&
      nrow(facet_response_patterns) > 0L &&
      "PatternStatus" %in% names(facet_response_patterns)) {
    response_pattern_count <- sum(
      !tolower(as.character(facet_response_patterns$PatternStatus)) %in% "ok",
      na.rm = TRUE
    )
  }

  unknown_count <- nrow(unknown_elements)

  rows <- list(
    data.frame(
      Area = "Rows",
      Status = if (row_issue_count > 0) "review" else "ok",
      Count = as.integer(row_issue_count),
      Unit = "rows",
      PercentOfData = percent_of_data(row_issue_count),
      Message = if (row_issue_count > 0) {
        "Some rows were excluded or need row-level review."
      } else {
        "No row-level missingness, invalid weight, or out-of-range score issue was found."
      },
      NextStep = if (row_issue_count > 0) {
        "Inspect `row_review` and resolve missingness, weights, or scale-range issues."
      } else {
        "Continue to score-support and facet-level category-use checks."
      },
      stringsAsFactors = FALSE
    ),
    data.frame(
      Area = "Score support",
      Status = if (intermediate_zero_categories > 0 || original_gap_conditions > 0) {
        "high"
      } else if (zero_categories > 0 || recode_conditions > 0) {
        "review"
      } else {
        "ok"
      },
      Count = as.integer(score_support_count),
      Unit = if (zero_categories > 0) "categories" else "conditions",
      PercentOfData = NA_real_,
      Message = if (intermediate_zero_categories > 0) {
        "At least one intermediate score category has zero observations."
      } else if (original_gap_conditions > 0) {
        "Original score labels had gaps before internal recoding."
      } else if (zero_categories > 0) {
        "At least one declared score category has zero observations."
      } else if (recode_conditions > 0) {
        "Original score labels were internally recoded."
      } else {
        "No score-support gap was found over the fitted score scale."
      },
      NextStep = if (intermediate_zero_categories > 0 || zero_categories > 0) {
        "Inspect `score_support_review`; document weak threshold support where needed."
      } else if (original_gap_conditions > 0 || recode_conditions > 0) {
        "Inspect `score_map` before reporting category or threshold labels."
      } else {
        "Continue to facet-level category-use checks."
      },
      stringsAsFactors = FALSE
    ),
    data.frame(
      Area = "Facet category use",
      Status = if (facet_intermediate_zero_levels > 0) {
        "high"
      } else if (facet_zero_levels > 0 || facet_sparse_levels > 0 || facet_issue_count > 0) {
        "review"
      } else {
        "ok"
      },
      Count = as.integer(facet_issue_count),
      Unit = "facet levels",
      PercentOfData = NA_real_,
      Message = if (facet_intermediate_zero_levels > 0) {
        "Some facet levels do not use one or more intermediate score categories."
      } else if (facet_zero_levels > 0 || facet_sparse_levels > 0 || facet_issue_count > 0) {
        "Some facet levels have zero or sparse local category use."
      } else {
        "No facet-level zero or sparse category-use issue was found."
      },
      NextStep = if (facet_issue_count > 0) {
        "Inspect `category_usage_summary` and `category_usage_by_facet` before local facet conclusions."
      } else {
        "Continue to design-match checks."
      },
      stringsAsFactors = FALSE
    ),
    data.frame(
      Area = "Facet response patterns",
      Status = if (single_category_levels > 0) {
        "high"
      } else if (dominant_category_levels > 0 || boundary_only_levels > 0 || response_pattern_count > 0) {
        "review"
      } else {
        "ok"
      },
      Count = as.integer(response_pattern_count),
      Unit = "facet levels",
      PercentOfData = NA_real_,
      Message = if (single_category_levels > 0) {
        "Some facet levels used only one score category."
      } else if (dominant_category_levels > 0) {
        "Some facet levels assigned one score category to most responses."
      } else if (boundary_only_levels > 0 || response_pattern_count > 0) {
        "Some facet levels show restricted local response patterns."
      } else {
        "No single-category or dominant-category facet response pattern was found."
      },
      NextStep = if (response_pattern_count > 0) {
        "Inspect `facet_response_patterns`; for raters, restricted response patterns may indicate scoring-use problems."
      } else {
        "Continue to design-match checks."
      },
      stringsAsFactors = FALSE
    ),
    data.frame(
      Area = "Design match",
      Status = if (unknown_count > 0) "review" else "ok",
      Count = as.integer(unknown_count),
      Unit = "levels",
      PercentOfData = NA_real_,
      Message = if (unknown_count > 0) {
        "Raw data include levels that are absent from the fitted design."
      } else {
        "No raw-data level outside the fitted design was found."
      },
      NextStep = if (unknown_count > 0) {
        "Inspect `unknown_elements`; check labels or refit with the intended design."
      } else {
        "Proceed with estimation diagnostics and reporting checks."
      },
      stringsAsFactors = FALSE
    )
  )

  out <- do.call(rbind, rows)
  if (nrow(quality_flags) > 0L && "Severity" %in% names(quality_flags)) {
    out$QualityFlags <- vapply(out$Area, function(area) {
      sum(quality_flags$Area %in% area, na.rm = TRUE)
    }, integer(1))
    out$HighSeverityFlags <- vapply(out$Area, function(area) {
      sum(quality_flags$Area %in% area &
            tolower(as.character(quality_flags$Severity)) %in% "high", na.rm = TRUE)
    }, integer(1))
  } else {
    out$QualityFlags <- 0L
    out$HighSeverityFlags <- 0L
  }
  status_rank <- c(high = 1L, review = 2L, ok = 3L)
  out$StatusRank <- unname(status_rank[tolower(as.character(out$Status))])
  out$StatusRank[is.na(out$StatusRank)] <- 99L
  out <- out[order(out$StatusRank, out$Area), , drop = FALSE]
  out$StatusRank <- NULL
  row.names(out) <- NULL
  out
}

empty_mfrm_caveats <- function() {
  data.frame(
    Area = character(0),
    Severity = character(0),
    Condition = character(0),
    Categories = character(0),
    CategoryType = character(0),
    Message = character(0),
    RecommendedAction = character(0),
    Details = character(0),
    stringsAsFactors = FALSE
  )
}

collect_mfrm_caveats <- function(fit = NULL,
                                 prep = NULL,
                                 score_distribution = NULL,
                                 include_recode = TRUE,
                                 context = c("fit", "data")) {
  context <- match.arg(context)
  if (!is.null(fit) && is.null(prep)) {
    prep <- fit$prep %||% NULL
  }
  profile <- score_category_support_profile(prep = prep, score_distribution = score_distribution)
  out <- empty_mfrm_caveats()
  support_phrase <- if (identical(context, "data")) "prepared score support" else "fitted score support"

  add_caveat <- function(severity, condition, categories, category_type, message, action) {
    out <<- rbind(
      out,
      data.frame(
        Area = "score_categories",
        Severity = severity,
        Condition = condition,
        Categories = paste(as.character(categories), collapse = ", "),
        CategoryType = category_type,
        Message = message,
        RecommendedAction = action,
        Details = "",
        stringsAsFactors = FALSE
      )
    )
    invisible(NULL)
  }

  zero_rows <- if (nrow(profile) > 0 && "ZeroCount" %in% names(profile)) {
    profile[as.logical(profile$ZeroCount), , drop = FALSE]
  } else {
    profile[0, , drop = FALSE]
  }
  if (nrow(zero_rows) > 0) {
    internal <- zero_rows$Category[zero_rows$UnusedCategoryType == "internal"]
    boundary <- zero_rows$Category[zero_rows$UnusedCategoryType == "boundary"]
    if (length(internal) > 0L) {
      internal_text <- paste(as.character(internal), collapse = ", ")
      add_caveat(
        severity = "warning",
        condition = "zero_count_intermediate_score_category",
        categories = internal,
        category_type = "internal",
        message = paste0(
          "Unused intermediate score categories retained in the ", support_phrase, ": ",
          internal_text,
          ". Adjacent threshold estimates are weakly identified; review `rating_scale_table()` / `category_structure_report()` and consider category collapsing before treating the thresholds as stable."
        ),
        action = "Review adjacent thresholds and category curves; collapse categories or collect additional data when threshold stability is required."
      )
    }
    if (length(boundary) > 0L) {
      boundary_text <- paste(as.character(boundary), collapse = ", ")
      add_caveat(
        severity = "review",
        condition = "zero_count_boundary_score_category",
        categories = boundary,
        category_type = "boundary",
        message = paste0(
          "Unused boundary score categories retained in the ", support_phrase, ": ",
          boundary_text,
          ". Document the zero-count category and avoid overinterpreting adjacent threshold estimates."
        ),
        action = "Document the intended score support and verify `rating_min` / `rating_max` before reporting category functioning."
      )
    }
  }

  if (isTRUE(include_recode) && !is.null(prep)) {
    score_map <- as.data.frame(prep$score_map %||% data.frame(), stringsAsFactors = FALSE)
    if (nrow(score_map) > 0L &&
        all(c("OriginalScore", "InternalScore") %in% names(score_map)) &&
        any(as.character(score_map$OriginalScore) != as.character(score_map$InternalScore))) {
      add_caveat(
        severity = "info",
        condition = "score_categories_recoded",
        categories = sort(unique(as.character(score_map$OriginalScore))),
        category_type = "recoded",
        message = "Observed score categories were internally recoded; inspect `fit$prep$score_map` before interpreting category labels.",
        action = "Use `fit$prep$score_map` when connecting internal category estimates back to original score labels."
      )
    }
    if (nrow(score_map) > 0L && "OriginalScore" %in% names(score_map)) {
      original_scores <- sort(unique(suppressWarnings(as.integer(score_map$OriginalScore))))
      original_scores <- original_scores[is.finite(original_scores)]
      if (length(original_scores) >= 2L) {
        original_support <- seq(min(original_scores), max(original_scores))
        missing_original <- setdiff(original_support, original_scores)
        retained_zero <- suppressWarnings(as.integer(zero_rows$Category %||% integer(0)))
        retained_zero <- retained_zero[is.finite(retained_zero)]
        missing_not_retained <- setdiff(missing_original, retained_zero)
        if (length(missing_not_retained) > 0L) {
          add_caveat(
            severity = "warning",
            condition = "original_score_gap_before_recoding",
            categories = missing_not_retained,
            category_type = "original_gap",
            message = paste0(
              "Original score labels skip category ",
              paste(missing_not_retained, collapse = ", "),
              " before the fitted internal scale is formed. If the intended scale includes ",
              if (length(missing_not_retained) == 1L) "this category" else "these categories",
              ", refit with `rating_min`, `rating_max`, and `keep_original = TRUE` so zero-frequency categories remain visible in category-functioning output."
            ),
            action = "Check the declared rating scale; retain intended zero-count intermediate categories when threshold functioning is part of the analysis."
          )
        }
      }
    }
  }

  out
}

annotate_score_category_caveats <- function(category_table, prep = NULL, score_distribution = NULL) {
  category_table <- as.data.frame(category_table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(category_table) == 0 || !"Category" %in% names(category_table)) {
    return(category_table)
  }
  profile <- score_category_support_profile(prep = prep, score_distribution = score_distribution)
  if (nrow(profile) == 0) {
    category_table$ZeroCount <- FALSE
    category_table$UnusedCategoryType <- "none"
    category_table$WeaklyIdentified <- FALSE
    category_table$CategoryCaveat <- ""
    return(category_table)
  }

  category_table |>
    dplyr::mutate(.CategoryKey = as.character(.data$Category)) |>
    dplyr::left_join(
      profile |>
        dplyr::mutate(.CategoryKey = as.character(.data$Category)) |>
        dplyr::select(
          ".CategoryKey",
          "ZeroCount",
          "UnusedCategoryType",
          "WeaklyIdentified",
          "CategoryCaveat"
        ),
      by = ".CategoryKey"
    ) |>
    dplyr::mutate(
      ZeroCount = dplyr::coalesce(.data$ZeroCount, FALSE),
      UnusedCategoryType = dplyr::coalesce(.data$UnusedCategoryType, "none"),
      WeaklyIdentified = dplyr::coalesce(.data$WeaklyIdentified, FALSE),
      CategoryCaveat = dplyr::coalesce(.data$CategoryCaveat, "")
    ) |>
    dplyr::select(-dplyr::all_of(".CategoryKey")) |>
    as.data.frame(stringsAsFactors = FALSE)
}

step_category_bounds <- function(step, prep = NULL) {
  step <- as.character(step)
  parts <- strsplit(step, "-", fixed = TRUE)
  lower <- rep(NA_integer_, length(parts))
  upper <- rep(NA_integer_, length(parts))
  for (i in seq_along(parts)) {
    if (length(parts[[i]]) >= 2L) {
      lower[i] <- suppressWarnings(as.integer(parts[[i]][1]))
      upper[i] <- suppressWarnings(as.integer(parts[[i]][2]))
    }
  }
  missing_bounds <- !is.finite(lower) | !is.finite(upper)
  if (any(missing_bounds) && !is.null(prep)) {
    rating_min <- suppressWarnings(as.integer(prep$rating_min %||% NA_integer_))
    rating_max <- suppressWarnings(as.integer(prep$rating_max %||% NA_integer_))
    if (is.finite(rating_min) && is.finite(rating_max) && rating_min <= rating_max) {
      support <- seq(rating_min, rating_max)
      step_index <- step_index_from_label(step)
      valid <- missing_bounds &
        is.finite(step_index) &
        step_index >= 1L &
        step_index < length(support)
      lower[valid] <- support[step_index[valid]]
      upper[valid] <- support[step_index[valid] + 1L]
    }
  }
  data.frame(LowerCategory = lower, UpperCategory = upper)
}

annotate_threshold_caveats <- function(threshold_table, prep = NULL) {
  threshold_table <- as.data.frame(threshold_table %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(threshold_table) == 0) {
    return(threshold_table)
  }
  profile <- score_category_support_profile(prep = prep)
  if (nrow(profile) == 0 || !"Step" %in% names(threshold_table)) {
    threshold_table$WeaklyIdentified <- FALSE
    threshold_table$ThresholdCaveat <- ""
    return(threshold_table)
  }

  if (!all(c("LowerCategory", "UpperCategory") %in% names(threshold_table))) {
    bounds <- step_category_bounds(threshold_table$Step, prep = prep)
    threshold_table$LowerCategory <- bounds$LowerCategory
    threshold_table$UpperCategory <- bounds$UpperCategory
  }

  profile_key <- profile |>
    dplyr::mutate(.CategoryKey = as.character(.data$Category)) |>
    dplyr::select(
      ".CategoryKey",
      "ZeroCount",
      "UnusedCategoryType",
      "CategoryCaveat"
    )

  annotated <- threshold_table |>
    dplyr::mutate(
      .LowerCategoryKey = as.character(.data$LowerCategory),
      .UpperCategoryKey = as.character(.data$UpperCategory)
    ) |>
    dplyr::left_join(
      profile_key |>
        dplyr::rename(
          LowerZeroCount = "ZeroCount",
          LowerUnusedCategoryType = "UnusedCategoryType",
          LowerCategoryCaveat = "CategoryCaveat"
        ),
      by = c(".LowerCategoryKey" = ".CategoryKey")
    ) |>
    dplyr::left_join(
      profile_key |>
        dplyr::rename(
          UpperZeroCount = "ZeroCount",
          UpperUnusedCategoryType = "UnusedCategoryType",
          UpperCategoryCaveat = "CategoryCaveat"
        ),
      by = c(".UpperCategoryKey" = ".CategoryKey")
    ) |>
    dplyr::mutate(
      LowerZeroCount = dplyr::coalesce(.data$LowerZeroCount, FALSE),
      UpperZeroCount = dplyr::coalesce(.data$UpperZeroCount, FALSE),
      LowerUnusedCategoryType = dplyr::coalesce(.data$LowerUnusedCategoryType, "none"),
      UpperUnusedCategoryType = dplyr::coalesce(.data$UpperUnusedCategoryType, "none"),
      WeaklyIdentified = .data$LowerZeroCount | .data$UpperZeroCount,
      ThresholdCaveat = dplyr::case_when(
        .data$WeaklyIdentified & (
          .data$LowerUnusedCategoryType == "internal" |
            .data$UpperUnusedCategoryType == "internal"
        ) ~ "Adjacent to a zero-count intermediate category; threshold estimate is weakly identified.",
        .data$WeaklyIdentified ~ "Adjacent to a zero-count boundary category; document support before interpreting this threshold.",
        TRUE ~ ""
      )
    ) |>
    dplyr::select(
      -dplyr::all_of(c(
        ".LowerCategoryKey",
        ".UpperCategoryKey",
        "LowerZeroCount",
        "UpperZeroCount",
        "LowerUnusedCategoryType",
        "UpperUnusedCategoryType",
        "LowerCategoryCaveat",
        "UpperCategoryCaveat"
      ))
    )

  as.data.frame(annotated, stringsAsFactors = FALSE)
}

augment_category_table_with_marginal_fit <- function(category_table, diagnostics) {
  if (!has_marginal_fit_bundle(diagnostics)) {
    return(category_table)
  }

  marginal_cells <- as.data.frame(
    diagnostics$marginal_fit$overall$cell_stats %||% data.frame(),
    stringsAsFactors = FALSE
  )
  if (nrow(marginal_cells) == 0 || !"Category" %in% names(marginal_cells)) {
    return(category_table)
  }

  keep <- intersect(
    c(
      "Category",
      "ObservedCount",
      "ExpectedCount",
      "ResidualCount",
      "ObservedProp",
      "ExpectedProp",
      "PropDiff",
      "StdResidual",
      "FlaggedAbsZ"
    ),
    names(marginal_cells)
  )
  marginal_cells <- marginal_cells[, keep, drop = FALSE]

  rename_map <- c(
    ObservedCount = "MarginalObservedCount",
    ExpectedCount = "MarginalExpectedCount",
    ResidualCount = "MarginalResidualCount",
    ObservedProp = "MarginalObservedProp",
    ExpectedProp = "MarginalExpectedProp",
    PropDiff = "MarginalPropDiff",
    StdResidual = "MarginalStdResidual",
    FlaggedAbsZ = "MarginalFitFlag"
  )
  for (nm in intersect(names(rename_map), names(marginal_cells))) {
    names(marginal_cells)[names(marginal_cells) == nm] <- rename_map[[nm]]
  }

  category_table |>
    dplyr::mutate(.CategoryKey = as.character(.data$Category)) |>
    dplyr::left_join(
      marginal_cells |>
        dplyr::mutate(.CategoryKey = as.character(.data$Category)) |>
        dplyr::select(-dplyr::all_of("Category")),
      by = ".CategoryKey"
    ) |>
    dplyr::select(-dplyr::all_of(".CategoryKey")) |>
    as.data.frame(stringsAsFactors = FALSE)
}

#' Build a legacy-compatible Table 8 bar-chart style scale-structure export
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()].
#' @param theta_range Theta/logit range used to derive mode/mean transition points.
#' @param theta_points Number of grid points used for transition-point search.
#' @param drop_unused If `TRUE`, remove zero-count categories from `category_table`.
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#' @param fixed_max_rows Maximum rows per section in the fixed-width text output.
#'
#' @details
#' The legacy-compatible Table 8 bar-chart output describes category structure with
#' mode / median / mean landmarks along the latent scale.
#' This helper returns those landmarks as numeric tables:
#' - mode peaks and mode transition points
#' - median thresholds (step calibrations)
#' - mean half-score transition points
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [category_structure_report()].
#'
#' @return A named list with:
#' - `category_table`: observed/expected category counts and fit
#' - `mode_peaks`: peak theta/probability by group and category
#' - `mode_boundaries`: theta points where modal category changes
#' - `median_thresholds`: threshold table (step-based), including
#'   weak-identification caveats for thresholds adjacent to retained
#'   zero-count categories
#' - `mean_halfscore_points`: theta points where expected score crosses half-scores
#' - `caveats`: structured score-support warning/review rows
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#' - `settings`: applied options
#'
#' @seealso [rating_scale_table()], [category_curves_report()], [plot.mfrm_fit()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' t8b <- category_structure_report(fit)
#' @keywords internal
#' @noRd
table8_barchart_export <- function(fit,
                                   diagnostics = NULL,
                                   theta_range = c(-6, 6),
                                   theta_points = 241,
                                   drop_unused = FALSE,
                                   include_fixed = FALSE,
                                   fixed_max_rows = 200) {
  signal_legacy_name_deprecation(
    old_name = "table8_barchart_export",
    new_name = "category_structure_report",
    suppress_if_called_from = "category_structure_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  fit_model <- as.character(fit$config$model %||% fit$summary$Model[1] %||% NA_character_)
  if (is.null(diagnostics) && !identical(fit_model, "GPCM")) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  obs_tbl <- if (!is.null(diagnostics$obs) && nrow(diagnostics$obs) > 0) {
    diagnostics$obs
  } else if (identical(fit_model, "GPCM")) {
    compute_obs_table(fit)
  } else {
    NULL
  }
  if (is.null(obs_tbl) || nrow(obs_tbl) == 0) {
    stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
  }
  theta_points <- max(51L, as.integer(theta_points))
  theta_range <- as.numeric(theta_range)
  if (length(theta_range) != 2 || !all(is.finite(theta_range)) || theta_range[1] >= theta_range[2]) {
    stop("`theta_range` must be a numeric length-2 vector with increasing values.")
  }

  category_table <- as.data.frame(calc_category_stats(obs_tbl, res = fit, whexact = FALSE), stringsAsFactors = FALSE)
  category_table <- augment_category_table_with_marginal_fit(category_table, diagnostics)
  category_table <- annotate_score_category_caveats(category_table, prep = fit$prep)
  if (isTRUE(drop_unused) && nrow(category_table) > 0 && "Count" %in% names(category_table)) {
    category_table <- category_table[category_table$Count > 0, , drop = FALSE]
  }

  curve_spec <- build_step_curve_spec(fit)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  curve_tbl <- build_curve_tables(curve_spec, theta_grid)
  prob_df <- as.data.frame(curve_tbl$probabilities, stringsAsFactors = FALSE)
  exp_df <- as.data.frame(curve_tbl$expected, stringsAsFactors = FALSE)

  mode_peaks <- prob_df |>
    dplyr::group_by(.data$CurveGroup, .data$Category) |>
    dplyr::arrange(dplyr::desc(.data$Probability), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::transmute(
      CurveGroup = .data$CurveGroup,
      Category = .data$Category,
      PeakTheta = .data$Theta,
      PeakProbability = .data$Probability
    ) |>
    dplyr::ungroup() |>
    as.data.frame(stringsAsFactors = FALSE)

  mode_boundaries <- prob_df |>
    dplyr::group_by(.data$CurveGroup, .data$Theta) |>
    dplyr::arrange(dplyr::desc(.data$Probability), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::transmute(
      CurveGroup = .data$CurveGroup,
      Theta = .data$Theta,
      ModalCategory = .data$Category
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$CurveGroup) |>
    dplyr::arrange(.data$Theta, .by_group = TRUE) |>
    dplyr::mutate(
      PrevCategory = dplyr::lag(.data$ModalCategory),
      PrevTheta = dplyr::lag(.data$Theta)
    ) |>
    dplyr::filter(
      !is.na(.data$PrevCategory),
      .data$ModalCategory != .data$PrevCategory
    ) |>
    dplyr::transmute(
      CurveGroup = .data$CurveGroup,
      LowerCategory = .data$PrevCategory,
      UpperCategory = .data$ModalCategory,
      ModeBoundaryTheta = (.data$PrevTheta + .data$Theta) / 2
    ) |>
    dplyr::ungroup() |>
    as.data.frame(stringsAsFactors = FALSE)

  median_thresholds <- curve_spec$step_points |>
    dplyr::mutate(
      StepIndex = as.integer(.data$StepIndex),
      LowerCategory = as.character(curve_spec$categories[.data$StepIndex]),
      UpperCategory = as.character(curve_spec$categories[.data$StepIndex + 1]),
      MedianThreshold = .data$Threshold
    ) |>
    dplyr::select("CurveGroup", "Step", "StepIndex", "LowerCategory", "UpperCategory", "MedianThreshold") |>
    annotate_threshold_caveats(prep = fit$prep) |>
    as.data.frame(stringsAsFactors = FALSE)

  cat_values <- suppressWarnings(as.numeric(curve_spec$categories))
  if (!all(is.finite(cat_values))) {
    cat_values <- seq_along(curve_spec$categories) - 1
  }
  mean_halfscore_points <- purrr::map_dfr(unique(exp_df$CurveGroup), function(g) {
    sub <- exp_df[exp_df$CurveGroup == g, , drop = FALSE]
    if (nrow(sub) == 0) return(tibble::tibble())
    mids <- (head(cat_values, -1) + tail(cat_values, -1)) / 2
    purrr::map_dfr(seq_along(mids), function(i) {
      tibble::tibble(
        CurveGroup = g,
        LowerCategory = as.character(curve_spec$categories[i]),
        UpperCategory = as.character(curve_spec$categories[i + 1]),
        MeanHalfScore = mids[i],
        MeanBoundaryTheta = closest_theta_for_target(sub$Theta, sub$ExpectedScore, mids[i])
      )
    })
  }) |>
    as.data.frame(stringsAsFactors = FALSE)

  out <- list(
    category_table = category_table,
    mode_peaks = mode_peaks,
    mode_boundaries = mode_boundaries,
    median_thresholds = median_thresholds,
    mean_halfscore_points = mean_halfscore_points,
    caveats = collect_mfrm_caveats(fit = fit),
    diagnostic_mode = as.character(diagnostics$diagnostic_mode %||% "legacy"),
    marginal_fit = diagnostics$marginal_fit %||% NULL,
    settings = list(
      theta_range = theta_range,
      theta_points = theta_points,
      drop_unused = isTRUE(drop_unused),
      include_fixed = isTRUE(include_fixed),
      fixed_max_rows = max(10L, as.integer(fixed_max_rows))
    )
  )
  if (isTRUE(include_fixed)) {
    fixed_max_rows <- max(10L, as.integer(fixed_max_rows))
    out$fixed <- build_sectioned_fixed_report(
      title = "Legacy-compatible Table 8 Bar-chart Export",
      sections = list(
        list(title = "Category table", data = category_table, max_rows = fixed_max_rows),
        list(title = "Mode peaks", data = mode_peaks, max_rows = fixed_max_rows),
        list(title = "Mode boundaries", data = mode_boundaries, max_rows = fixed_max_rows),
        list(title = "Median thresholds", data = median_thresholds, max_rows = fixed_max_rows),
        list(title = "Mean half-score transition points", data = mean_halfscore_points, max_rows = fixed_max_rows)
      ),
      max_col_width = 24
    )
  }
  out
}

#' Build a legacy-compatible graphfile-style Table 8 curves export
#'
#' @param fit Output from [fit_mfrm()].
#' @param theta_range Theta/logit range for curve coordinates.
#' @param theta_points Number of points on the theta grid.
#' @param digits Rounding digits for numeric graph output.
#' @param include_fixed If `TRUE`, include a legacy-compatible fixed-width text
#'   block.
#' @param fixed_max_rows Maximum rows shown in the fixed-width graph table.
#'
#' @details
#' The legacy-compatible `Graphfile=` output for Table 8 contains curve
#' coordinates:
#' scale number, measure, expected score, expected category, and
#' category probabilities by measure. This helper returns the same
#' information as data frames ready for CSV export.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [category_curves_report()].
#'
#' @return A named list with:
#' - `graphfile`: wide table with `Scale, Measure, Expected, ExpCat, Prob:*`
#' - `graphfile_syntactic`: same content with syntactic probability column names
#' - `probabilities`: long probability table (`Theta`, `Category`, `CurveGroup`)
#' - `cumulative_probabilities`: long cumulative probability table
#' - `cumulative_boundaries`: approximate 0.5 cumulative-probability thresholds
#' - `category_information`: long category-specific information contribution table
#' - `expected_ogive`: long expected-score table
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#' - `settings`: applied options
#'
#' @seealso [category_structure_report()], [rating_scale_table()], [plot.mfrm_fit()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' t8c <- category_curves_report(fit, theta_points = 101)
#' @keywords internal
#' @noRd
table8_curves_export <- function(fit,
                                 theta_range = c(-6, 6),
                                 theta_points = 241,
                                 digits = 4,
                                 include_fixed = FALSE,
                                 fixed_max_rows = 400) {
  signal_legacy_name_deprecation(
    old_name = "table8_curves_export",
    new_name = "category_curves_report",
    suppress_if_called_from = c("category_curves_report", "facets_output_file_bundle")
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  theta_points <- max(51L, as.integer(theta_points))
  theta_range <- as.numeric(theta_range)
  if (length(theta_range) != 2 || !all(is.finite(theta_range)) || theta_range[1] >= theta_range[2]) {
    stop("`theta_range` must be a numeric length-2 vector with increasing values.")
  }
  digits <- max(0L, as.integer(digits))

  curve_spec <- build_step_curve_spec(fit)
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  curve_tbl <- build_curve_tables(curve_spec, theta_grid)
  prob_df <- as.data.frame(curve_tbl$probabilities, stringsAsFactors = FALSE)
  exp_df <- as.data.frame(curve_tbl$expected, stringsAsFactors = FALSE)
  cumulative_tbls <- build_cumulative_category_tables(prob_df, curve_spec$categories)

  groups <- unique(as.character(exp_df$CurveGroup))
  scale_map <- setNames(seq_along(groups), groups)

  prob_wide <- prob_df |>
    dplyr::mutate(
      Category = as.character(.data$Category),
      ProbCol = paste0("Prob:", .data$Category)
    ) |>
    dplyr::select("CurveGroup", "Theta", "ProbCol", "Probability") |>
    tidyr::pivot_wider(names_from = "ProbCol", values_from = "Probability")

  graph_tbl <- exp_df |>
    dplyr::left_join(prob_wide, by = c("CurveGroup", "Theta")) |>
    dplyr::mutate(
      Scale = as.integer(scale_map[as.character(.data$CurveGroup)]),
      Measure = .data$Theta
    )

  prob_cols <- grep("^Prob:", names(graph_tbl), value = TRUE)
  if (length(prob_cols) == 0) {
    stop("No probability columns were generated for Table 8 curves export.")
  }
  category_values <- suppressWarnings(as.numeric(sub("^Prob:", "", prob_cols)))
  if (!all(is.finite(category_values))) {
    category_values <- seq_along(prob_cols) - 1
  }
  expected_vec <- suppressWarnings(as.numeric(graph_tbl$ExpectedScore))
  nearest_idx <- vapply(expected_vec, function(ev) {
    if (!is.finite(ev)) return(NA_integer_)
    as.integer(which.min(abs(ev - category_values)))
  }, integer(1))
  exp_cat <- ifelse(is.na(nearest_idx), NA_real_, category_values[nearest_idx])

  graph_tbl <- graph_tbl |>
    dplyr::mutate(ExpCat = exp_cat) |>
    dplyr::select("Scale", "CurveGroup", "Measure", "ExpectedScore", "ExpCat", dplyr::all_of(prob_cols)) |>
    dplyr::arrange(.data$Scale, .data$Measure)

  graph_tbl_syntactic <- graph_tbl
  names(graph_tbl_syntactic) <- gsub("^Prob:", "Prob_", names(graph_tbl_syntactic))
  names(graph_tbl_syntactic) <- sub("^ExpectedScore$", "Expected", names(graph_tbl_syntactic))

  graph_tbl_facets <- graph_tbl
  names(graph_tbl_facets)[names(graph_tbl_facets) == "ExpectedScore"] <- "Expected"
  graph_tbl_facets <- as.data.frame(graph_tbl_facets, check.names = FALSE, stringsAsFactors = FALSE)

  round_numeric <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)
    out <- df
    num_cols <- vapply(out, is.numeric, logical(1))
    out[num_cols] <- lapply(out[num_cols], round, digits = digits)
    out
  }

  out <- list(
    graphfile = round_numeric(graph_tbl_facets),
    graphfile_syntactic = round_numeric(graph_tbl_syntactic),
    probabilities = round_numeric(as.data.frame(prob_df, stringsAsFactors = FALSE)),
    cumulative_probabilities = round_numeric(cumulative_tbls$probabilities),
    cumulative_boundaries = round_numeric(cumulative_tbls$boundaries),
    category_information = round_numeric(
      as.data.frame(
        prob_df |>
          dplyr::select(
            "CurveGroup", "Theta", "Category", "Probability",
            "ExpectedScore", "ScoreVariance", "Information",
            "CategoryInformation", "CategoryInformationShare",
            "Slope", "Model"
          ),
        stringsAsFactors = FALSE
      )
    ),
    expected_ogive = round_numeric(as.data.frame(exp_df, stringsAsFactors = FALSE)),
    settings = list(
      theta_range = theta_range,
      theta_points = theta_points,
      digits = digits,
      include_fixed = isTRUE(include_fixed),
      fixed_max_rows = max(25L, as.integer(fixed_max_rows)),
      scales = as.data.frame(
        data.frame(
          Scale = as.integer(scale_map),
          CurveGroup = names(scale_map),
          stringsAsFactors = FALSE
        )
      )
    )
  )
  if (isTRUE(include_fixed)) {
    fixed_max_rows <- max(25L, as.integer(fixed_max_rows))
    out$fixed <- build_sectioned_fixed_report(
      title = "Legacy-compatible Table 8 Curves Graphfile",
      sections = list(
        list(title = "Graphfile wide table", data = out$graphfile, max_rows = fixed_max_rows),
        list(title = "Scale map", data = out$settings$scales),
        list(title = "Expected ogive", data = out$expected_ogive, max_rows = fixed_max_rows)
      ),
      max_col_width = 24
    )
  }
  out
}

build_cumulative_category_tables <- function(prob_df, categories) {
  prob_df <- as.data.frame(prob_df %||% data.frame(), stringsAsFactors = FALSE)
  categories_chr <- as.character(categories)
  if (nrow(prob_df) == 0L || length(categories_chr) == 0L) {
    return(list(probabilities = data.frame(), boundaries = data.frame()))
  }
  needed <- c("CurveGroup", "Theta", "Category", "Probability")
  if (!all(needed %in% names(prob_df))) {
    return(list(probabilities = data.frame(), boundaries = data.frame()))
  }

  row_list <- list()
  idx <- 1L
  groups <- unique(as.character(prob_df$CurveGroup))
  for (g in groups) {
    group_df <- prob_df[as.character(prob_df$CurveGroup) == g, , drop = FALSE]
    theta_vals <- sort(unique(suppressWarnings(as.numeric(group_df$Theta))))
    theta_vals <- theta_vals[is.finite(theta_vals)]
    for (theta in theta_vals) {
      sub <- group_df[suppressWarnings(as.numeric(group_df$Theta)) == theta, , drop = FALSE]
      ord <- match(as.character(sub$Category), categories_chr)
      sub <- sub[order(ord), , drop = FALSE]
      ord <- ord[order(ord)]
      keep <- !is.na(ord)
      sub <- sub[keep, , drop = FALSE]
      ord <- ord[keep]
      if (nrow(sub) == 0L) next
      probs <- suppressWarnings(as.numeric(sub$Probability))
      probs[!is.finite(probs)] <- NA_real_
      below <- cumsum(probs)
      above <- rev(cumsum(rev(probs)))
      for (j in seq_len(nrow(sub))) {
        category <- as.character(sub$Category[j])
        model <- as.character(sub$Model[j] %||% NA_character_)
        slope <- suppressWarnings(as.numeric(sub$Slope[j] %||% NA_real_))
        row_list[[idx]] <- data.frame(
          CurveGroup = g,
          Theta = theta,
          Direction = "at_or_below",
          BoundaryCategory = category,
          BoundaryOrder = as.integer(ord[j]),
          CategorySet = paste0("<= ", category),
          CumulativeProbability = below[j],
          Model = model,
          Slope = slope,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
        row_list[[idx]] <- data.frame(
          CurveGroup = g,
          Theta = theta,
          Direction = "at_or_above",
          BoundaryCategory = category,
          BoundaryOrder = as.integer(ord[j]),
          CategorySet = paste0(">= ", category),
          CumulativeProbability = above[j],
          Model = model,
          Slope = slope,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  cumulative <- if (length(row_list) == 0L) {
    data.frame()
  } else {
    dplyr::bind_rows(row_list) |>
      dplyr::arrange(.data$CurveGroup, .data$Direction, .data$BoundaryOrder, .data$Theta)
  }
  boundaries <- build_cumulative_boundary_table(cumulative, categories_chr)
  list(
    probabilities = as.data.frame(cumulative, stringsAsFactors = FALSE),
    boundaries = as.data.frame(boundaries, stringsAsFactors = FALSE)
  )
}

build_cumulative_boundary_table <- function(cumulative, categories_chr) {
  cumulative <- as.data.frame(cumulative %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(cumulative) == 0L || !"Direction" %in% names(cumulative)) {
    return(data.frame())
  }
  cumulative <- cumulative[cumulative$Direction == "at_or_below", , drop = FALSE]
  if (nrow(cumulative) == 0L) {
    return(data.frame())
  }
  max_order <- length(categories_chr)
  rows <- list()
  idx <- 1L
  groups <- unique(as.character(cumulative$CurveGroup))
  for (g in groups) {
    group_df <- cumulative[as.character(cumulative$CurveGroup) == g, , drop = FALSE]
    boundary_orders <- sort(unique(suppressWarnings(as.integer(group_df$BoundaryOrder))))
    boundary_orders <- boundary_orders[is.finite(boundary_orders) & boundary_orders < max_order]
    for (j in boundary_orders) {
      line <- group_df[suppressWarnings(as.integer(group_df$BoundaryOrder)) == j, , drop = FALSE]
      line <- line[order(suppressWarnings(as.numeric(line$Theta))), , drop = FALSE]
      theta <- suppressWarnings(as.numeric(line$Theta))
      y <- suppressWarnings(as.numeric(line$CumulativeProbability))
      ok <- is.finite(theta) & is.finite(y)
      theta <- theta[ok]
      y <- y[ok]
      threshold <- NA_real_
      in_range <- FALSE
      crossing_count <- 0L
      if (length(theta) > 0L) {
        centered <- y - 0.5
        exact <- which(abs(centered) <= sqrt(.Machine$double.eps))
        if (length(exact) > 0L) {
          threshold <- theta[exact[1]]
          in_range <- TRUE
          crossing_count <- length(exact)
        } else if (length(theta) > 1L) {
          crossing <- which(centered[-length(centered)] * centered[-1] <= 0)
          crossing_count <- length(crossing)
          if (length(crossing) > 0L) {
            i <- crossing[1]
            y1 <- y[i]
            y2 <- y[i + 1L]
            x1 <- theta[i]
            x2 <- theta[i + 1L]
            threshold <- if (is.finite(y1) && is.finite(y2) && abs(y2 - y1) > sqrt(.Machine$double.eps)) {
              x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)
            } else {
              mean(c(x1, x2))
            }
            in_range <- TRUE
          }
        }
      }
      rows[[idx]] <- data.frame(
        CurveGroup = g,
        BoundaryOrder = as.integer(j),
        LowerOrEqualCategory = categories_chr[j],
        AboveCategory = categories_chr[j + 1L],
        ThresholdCategory = categories_chr[j + 1L],
        CumulativeDirection = "at_or_below",
        TargetProbability = 0.5,
        ThurstonianThreshold = threshold,
        InThetaRange = in_range,
        CrossingCount = as.integer(crossing_count),
        BoundaryStatus = dplyr::case_when(
          !in_range ~ "outside_theta_range",
          crossing_count == 1L ~ "in_range",
          crossing_count > 1L ~ "multiple_crossings",
          TRUE ~ "review"
        ),
        BoundaryLabel = paste0("P(X <= ", categories_chr[j], ") = 0.5"),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  if (length(rows) == 0L) {
    return(data.frame())
  }
  dplyr::bind_rows(rows) |>
    dplyr::arrange(.data$CurveGroup, .data$BoundaryOrder)
}

#' Build a legacy-compatible output-file bundle (`GRAPH=` / `SCORE=`)
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()] (used for score file).
#' @param include Output components to include: `"graph"` and/or `"score"`.
#' @param theta_range Theta/logit range for graph coordinates.
#' @param theta_points Number of points on the theta grid for graph coordinates.
#' @param digits Rounding digits for numeric fields.
#' @param score_se_method For bounded `GPCM` scorefile exports, which
#'   observation-level score uncertainty columns to compute. `"both"`
#'   (default) includes native structural expected-score SEs and score-side
#'   delta-method SEs; `"native"` includes only the structural expected-score
#'   route; `"score_side"` includes only the score-side delta route; `"none"`
#'   records explicit `not_requested` status columns.
#' @param include_fixed If `TRUE`, include fixed-width text mirrors of output tables.
#' @param fixed_max_rows Maximum rows shown in fixed-width text blocks.
#' @param write_files If `TRUE`, write selected outputs to files in `output_dir`.
#' @param output_dir Output directory used when `write_files = TRUE`.
#' @param file_prefix Prefix used for output file names.
#' @param overwrite If `FALSE`, existing output files are not overwritten.
#'
#' @details
#' Legacy-compatible output files often include:
#' - graph coordinates for Table 8 curves (`GRAPH=` / `Graphfile=`), and
#' - observation-level modeled score lines (`SCORE=`-style inspection).
#'
#' This helper returns both as data frames and can optionally write
#' CSV/fixed-width text files to disk.
#'
#' `summary(out)` is supported through `summary()`.
#' `plot(out)` is dispatched through `plot()` for class
#' `mfrm_output_bundle` (`type = "graph_expected"`, `"score_residuals"`,
#' `"obs_probability"`, `"score_se"`).
#'
#' @section Interpreting output:
#' - `graphfile`: legacy-compatible wide curve coordinates (human-readable labels).
#' - `graphfile_syntactic`: same curves with syntactic column names for programmatic use.
#' - `scorefile`: observation-level observed/expected/residual diagnostics.
#' - `written_files`: traceability record of files produced when `write_files = TRUE`.
#'
#' For reproducible pipelines, prefer `graphfile_syntactic` and keep
#' `written_files` in run logs.
#'
#' @section Preferred route for new analyses:
#' For new scripts, prefer [category_curves_report()] or
#' [category_structure_report()] for scale outputs, then use
#' [export_mfrm_bundle()] for file handoff. Use
#' `facets_output_file_bundle()` only when a legacy-compatible graphfile or
#' scorefile contract is required.
#'
#' @section Bounded GPCM boundary:
#' For bounded `GPCM`, graph output and package-native scorefile output are
#' available with caveats. `include = "score"` returns observation-level fitted
#' expected score, residual, standardized residual, observed-category
#' probability, GPCM slope fields, and native structural delta-method
#' expected-score uncertainty and/or score-side delta-method SEs when the
#' required MML diagnostics are available. Use `score_se_method` to choose
#' `"both"` (default), `"native"`, `"score_side"`, or `"none"`.
#' The scorefile also carries explicit score-side caveat columns. It is not a
#' FACETS score-side equivalence file, does not export FACETS-equivalent
#' score-side standard errors, and does not establish an operational
#' score-scale decision. Use [gpcm_score_side_contract()] and
#' [gpcm_capability_matrix()] for the current scope.
#'
#' @section Typical workflow:
#' 1. Fit and diagnose model.
#' 2. Generate bundle with `include = c("graph", "score")`.
#' 3. Validate with `summary(out)` / `plot(out)`.
#' 4. Export with `write_files = TRUE` for reporting handoff.
#'
#' @return A named list including:
#' - `graphfile` / `graphfile_syntactic` when `"graph"` is requested
#' - `scorefile` when `"score"` is requested
#' - `graphfile_fixed` / `scorefile_fixed` when `include_fixed = TRUE`
#' - `written_files` when `write_files = TRUE`
#' - `settings`: applied options
#'
#' @seealso [category_curves_report()], [diagnose_mfrm()], [unexpected_response_table()],
#'   [export_mfrm_bundle()], [mfrmr_reports_and_tables],
#'   [mfrmr_compatibility_layer]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' out <- facets_output_file_bundle(fit, diagnostics = diagnose_mfrm(fit, residual_pca = "none"))
#' summary(out)
#' p_out <- plot(out, draw = FALSE)
#' p_out$data$plot
#' @export
facets_output_file_bundle <- function(fit,
                                      diagnostics = NULL,
                                      include = c("graph", "score"),
                                      theta_range = c(-6, 6),
                                      theta_points = 241,
                                      digits = 4,
                                      score_se_method = c("both", "native", "score_side", "none"),
                                      include_fixed = FALSE,
                                      fixed_max_rows = 400,
                                      write_files = FALSE,
                                      output_dir = NULL,
                                      file_prefix = "mfrmr_output",
                                      overwrite = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  include <- unique(tolower(as.character(include)))
  allowed <- c("graph", "score")
  bad <- setdiff(include, allowed)
  if (length(bad) > 0) {
    stop("Unsupported `include` values: ", paste(bad, collapse = ", "), ". Use: ", paste(allowed, collapse = ", "))
  }
  if (length(include) == 0) {
    stop("`include` must contain at least one of: graph, score.")
  }
  if (is.null(score_se_method)) {
    score_se_method <- "both"
  }
  score_se_method <- match.arg(
    tolower(as.character(score_se_method[1])),
    choices = c("both", "native", "score_side", "none")
  )
  fit_model <- as.character(fit$config$model %||% fit$summary$Model[1] %||% NA_character_)
  gpcm_score_side <- identical(fit_model, "GPCM") && "score" %in% include
  if (isTRUE(gpcm_score_side)) {
    stop_if_gpcm_out_of_scope(
      fit,
      "facets_output_file_bundle(include = \"score\")",
      supported = paste(
        "package-native bounded-GPCM scorefile export with explicit caveat",
        "columns; full FACETS output-contract review remains blocked"
      )
    )
  }
  digits <- max(0L, as.integer(digits))
  include_fixed <- isTRUE(include_fixed)
  fixed_max_rows <- max(25L, as.integer(fixed_max_rows))
  write_files <- isTRUE(write_files)
  overwrite <- isTRUE(overwrite)
  file_prefix <- as.character(file_prefix[1] %||% "mfrmr_output")
  if (!nzchar(file_prefix)) file_prefix <- "mfrmr_output"
  output_dir <- if (is.null(output_dir)) NULL else as.character(output_dir[1])

  if (write_files) {
    if (is.null(output_dir) || !nzchar(output_dir)) {
      stop("`output_dir` must be supplied when `write_files = TRUE`.")
    }
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(output_dir)) {
      stop("Failed to create `output_dir`: ", output_dir)
    }
  }

  written_files <- data.frame(
    Component = character(0),
    Format = character(0),
    Path = character(0),
    stringsAsFactors = FALSE
  )
  add_written <- function(component, format, path) {
    written_files <<- rbind(
      written_files,
      data.frame(
        Component = as.character(component),
        Format = as.character(format),
        Path = as.character(path),
        stringsAsFactors = FALSE
      )
    )
    invisible(NULL)
  }
  write_csv_if_needed <- function(df, filename, component) {
    if (!write_files) return(invisible(NULL))
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("Output file already exists: ", path, ". Set `overwrite = TRUE` to replace it.")
    }
    utils::write.csv(df, file = path, row.names = FALSE, na = "")
    add_written(component = component, format = "csv", path = path)
    invisible(NULL)
  }
  write_txt_if_needed <- function(text, filename, component) {
    if (!write_files) return(invisible(NULL))
    path <- file.path(output_dir, filename)
    if (file.exists(path) && !overwrite) {
      stop("Output file already exists: ", path, ". Set `overwrite = TRUE` to replace it.")
    }
    writeLines(as.character(text), con = path, useBytes = TRUE)
    add_written(component = component, format = "txt", path = path)
    invisible(NULL)
  }

  out <- list()
  if ("graph" %in% include) {
    graph <- with_legacy_name_warning_suppressed(
      table8_curves_export(
        fit = fit,
        theta_range = theta_range,
        theta_points = theta_points,
        digits = digits,
        include_fixed = include_fixed,
        fixed_max_rows = fixed_max_rows
      )
    )
    out$graphfile <- graph$graphfile
    out$graphfile_syntactic <- graph$graphfile_syntactic
    if (include_fixed) {
      out$graphfile_fixed <- graph$fixed
    }
    write_csv_if_needed(out$graphfile, paste0(file_prefix, "_graphfile.csv"), "graphfile")
    write_csv_if_needed(out$graphfile_syntactic, paste0(file_prefix, "_graphfile_syntactic.csv"), "graphfile_syntactic")
    if (include_fixed && !is.null(out$graphfile_fixed)) {
      write_txt_if_needed(out$graphfile_fixed, paste0(file_prefix, "_graphfile_fixed.txt"), "graphfile_fixed")
    }
  }

  if ("score" %in% include) {
    if (is.null(diagnostics)) {
      diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
    }
    if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
      stop("`diagnostics$obs` is empty. Run diagnose_mfrm() first.")
    }

    obs <- as.data.frame(diagnostics$obs, stringsAsFactors = FALSE)
    keep_cols <- c(
      "Person", fit$config$facet_names, "Observed", "Expected", "Residual",
      "StdResidual", "Var", "Weight", "score_k", "PersonMeasure",
      "ScoreSlope", "ScoreInformation", "ObservedScoreDerivative", "PrObserved"
    )
    keep_cols <- keep_cols[keep_cols %in% names(obs)]
    scorefile <- obs[, keep_cols, drop = FALSE]

    probs <- compute_prob_matrix(fit)
    if (!is.null(probs) && nrow(probs) == nrow(obs) && "score_k" %in% names(obs)) {
      k_idx <- suppressWarnings(as.integer(obs$score_k)) + 1L
      obs_prob <- vapply(seq_len(nrow(obs)), function(i) {
        k <- k_idx[i]
        if (is.finite(k) && k >= 1L && k <= ncol(probs)) {
          return(as.numeric(probs[i, k]))
        }
        NA_real_
      }, numeric(1))
      scorefile$ObsProb <- obs_prob
    } else if ("PrObserved" %in% names(scorefile) && !"ObsProb" %in% names(scorefile)) {
      scorefile$ObsProb <- scorefile$PrObserved
    }
    if (isTRUE(gpcm_score_side)) {
      if (score_se_method %in% c("both", "native")) {
        scorefile <- add_gpcm_scorefile_delta_se(scorefile, fit)
      } else {
        scorefile <- add_gpcm_scorefile_native_se_not_requested(scorefile)
      }
      if (score_se_method %in% c("both", "score_side")) {
        scorefile <- add_gpcm_score_side_delta_se(scorefile, fit, diagnostics = diagnostics)
      } else {
        scorefile <- add_gpcm_score_side_se_not_requested(scorefile)
      }
      scorefile$ScoreSideModel <- "bounded_GPCM"
      scorefile$ScoreSideStatus <- "supported_with_caveat"
      scorefile$ScoreSideEstimand <- "fitted_bounded_gpcm_expected_score"
      scorefile$ScoreSideCaveat <- paste(
        "Package-native bounded-GPCM scorefile: slope-aware expected score,",
        "residual, standardized residual, observed-category probability,",
        "score slope, native expected-score uncertainty, and score-side",
        "delta-method SEs are exported for sensitivity review when requested",
        "and available. This is not FACETS score-side equivalence, does not",
        "export FACETS-equivalent score-side standard errors, and is not an",
        "operational scoring decision."
      )
    }
    out$scorefile <- round_numeric_df(scorefile, digits = digits)
    if (isTRUE(gpcm_score_side)) {
      out$gpcm_score_side_contract <- gpcm_score_side_contract()
      out$gpcm_boundary <- gpcm_capability_boundary_table(
        fit,
        helper = "facets_output_file_bundle(include = \"score\")",
        extra_areas = c("FACETS output-contract score-side review")
      )
    }
    if (include_fixed) {
      out$scorefile_fixed <- build_sectioned_fixed_report(
        title = "Legacy-compatible SCORE Output",
        sections = list(
          list(title = "Score file", data = out$scorefile, max_rows = fixed_max_rows)
        ),
        max_col_width = 24
      )
    }
    write_csv_if_needed(out$scorefile, paste0(file_prefix, "_scorefile.csv"), "scorefile")
    if (isTRUE(gpcm_score_side)) {
      write_csv_if_needed(
        out$gpcm_score_side_contract,
        paste0(file_prefix, "_gpcm_score_side_contract.csv"),
        "gpcm_score_side_contract"
      )
      write_csv_if_needed(
        out$gpcm_boundary,
        paste0(file_prefix, "_gpcm_boundary.csv"),
        "gpcm_boundary"
      )
    }
    if (include_fixed && !is.null(out$scorefile_fixed)) {
      write_txt_if_needed(out$scorefile_fixed, paste0(file_prefix, "_scorefile_fixed.txt"), "scorefile_fixed")
    }
  }

  if (write_files) {
    out$written_files <- written_files
  }

  out$settings <- list(
    include = include,
    theta_range = as.numeric(theta_range),
    theta_points = as.integer(theta_points),
    digits = digits,
    include_fixed = include_fixed,
    fixed_max_rows = fixed_max_rows,
    write_files = write_files,
    output_dir = output_dir,
    file_prefix = file_prefix,
    overwrite = overwrite,
    score_se_method = score_se_method,
    score_side_status = if (isTRUE(gpcm_score_side)) "bounded_gpcm_supported_with_caveat" else "standard"
  )
  as_mfrm_bundle(out, "mfrm_output_bundle")
}

mfrm_resolve_file_format <- function(format = NULL, path = NULL) {
  if (is.null(format)) {
    ext <- tolower(tools::file_ext(as.character(path[1] %||% "")))
    if (ext %in% c("csv", "tsv")) {
      return(ext)
    }
    return("csv")
  }
  format <- tolower(as.character(format[1] %||% "csv"))
  match.arg(format, c("csv", "tsv"))
}

mfrm_validate_output_path <- function(path) {
  path <- path.expand(as.character(path[1] %||% ""))
  if (is.na(path) || !nzchar(path)) {
    stop("`path` must be a non-empty output file path.")
  }
  path
}

mfrm_ensure_output_parent <- function(path) {
  path <- mfrm_validate_output_path(path)
  parent <- dirname(path)
  if (!is.na(parent) && nzchar(parent) && !dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(parent)) {
    stop("Failed to create output directory: ", parent)
  }
  invisible(path)
}

mfrm_write_delimited_table <- function(df, path, format, overwrite = FALSE) {
  path <- mfrm_validate_output_path(path)
  format <- mfrm_resolve_file_format(format, path)
  mfrm_ensure_output_parent(path)
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Output file already exists: ", path, ". Set `overwrite = TRUE` to replace it.")
  }
  if (identical(format, "csv")) {
    utils::write.csv(df, file = path, row.names = FALSE, na = "")
  } else {
    utils::write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

mfrm_assert_output_writable <- function(path, overwrite = FALSE) {
  path <- mfrm_validate_output_path(path)
  mfrm_ensure_output_parent(path)
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Output file already exists: ", path, ". Set `overwrite = TRUE` to replace it.")
  }
  invisible(path)
}

mfrm_residual_table_for_export <- function(fit, diagnostics = NULL, include_probabilities = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  if (!is.null(diagnostics) && !is.null(diagnostics$obs)) {
    obs <- as.data.frame(diagnostics$obs, stringsAsFactors = FALSE)
  } else {
    obs <- as.data.frame(compute_obs_table(fit), stringsAsFactors = FALSE)
  }
  if (nrow(obs) == 0) {
    stop("No observation-level residual rows are available.")
  }

  keep_cols <- c(
    "Person", fit$config$facet_names, "Weight", "Score", "Observed",
    "Expected", "Residual", "StdResidual", "StdSq", "Var",
    "ScoreInformation", "ScoreSlope", "PrObserved", "ItemEntropy",
    "ItemVarLogP", "PersonMeasure"
  )
  keep_cols <- keep_cols[keep_cols %in% names(obs)]
  tbl <- obs[, keep_cols, drop = FALSE]

  if (isTRUE(include_probabilities)) {
    probs <- compute_prob_matrix(fit)
    if (!is.null(probs) && nrow(probs) == nrow(obs)) {
      categories <- fit$prep$rating_min:fit$prep$rating_max
      prob_tbl <- as.data.frame(probs, stringsAsFactors = FALSE)
      names(prob_tbl) <- paste0("PrCategory_", categories[seq_len(ncol(prob_tbl))])
      tbl <- cbind(tbl, prob_tbl)
    }
  }
  as.data.frame(tbl, stringsAsFactors = FALSE)
}

mfrm_subset_tables_for_export <- function(fit, diagnostics = NULL) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  subset_tbls <- NULL
  if (!is.null(diagnostics) && !is.null(diagnostics$subsets)) {
    subset_tbls <- diagnostics$subsets
  }
  if (is.null(subset_tbls) ||
      is.null(subset_tbls$summary) ||
      is.null(subset_tbls$nodes)) {
    obs <- if (!is.null(diagnostics) && !is.null(diagnostics$obs)) {
      as.data.frame(diagnostics$obs, stringsAsFactors = FALSE)
    } else {
      as.data.frame(compute_obs_table(fit), stringsAsFactors = FALSE)
    }
    subset_tbls <- calc_subsets(obs, c("Person", fit$config$facet_names))
  }
  list(
    summary = as.data.frame(subset_tbls$summary %||% data.frame(), stringsAsFactors = FALSE),
    nodes = as.data.frame(subset_tbls$nodes %||% data.frame(), stringsAsFactors = FALSE)
  )
}

mfrm_default_node_path <- function(path, format) {
  format <- mfrm_resolve_file_format(format, path)
  stem <- tools::file_path_sans_ext(path)
  paste0(stem, "_nodes.", format)
}

#' Write a standalone residual file
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. Supplying it
#'   avoids recomputing observation diagnostics.
#' @param path Output file path.
#' @param format File format: `"csv"` or `"tsv"`. If omitted, inferred from
#'   `path` when the extension is `.csv` or `.tsv`, otherwise `"csv"`.
#' @param digits Rounding digits for numeric columns.
#' @param overwrite If `FALSE`, existing files are not overwritten.
#' @param include_probabilities If `TRUE`, append model probabilities for all
#'   response categories as `PrCategory_*` columns.
#'
#' @details
#' The exported table is observation-level and model-native. It includes the
#' observed score, expected score, residual, standardized residual, variance,
#' score information, observed-category probability, and modeled person
#' measure when those quantities are available.
#'
#' This writer is separate from [facets_output_file_bundle()] because it is a
#' direct analysis handoff rather than a legacy graph/score bundle.
#'
#' @return A bundle with `table`, `summary`, `written_files`, and `settings`.
#' @seealso [diagnose_mfrm()], [facets_output_file_bundle()],
#'   [write_mfrm_subset_file()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' path <- tempfile(fileext = ".csv")
#' out <- write_mfrm_residual_file(fit, diag, path, overwrite = TRUE)
#' out$written_files
#' }
#' @export
write_mfrm_residual_file <- function(fit,
                                     diagnostics = NULL,
                                     path,
                                     format = c("csv", "tsv"),
                                     digits = 4,
                                     overwrite = FALSE,
                                     include_probabilities = FALSE) {
  if (missing(path)) {
    stop("`path` must be supplied.")
  }
  path <- mfrm_validate_output_path(path)
  format <- mfrm_resolve_file_format(if (missing(format)) NULL else format, path)
  digits <- max(0L, as.integer(digits))
  tbl <- round_numeric_df(
    mfrm_residual_table_for_export(
      fit = fit,
      diagnostics = diagnostics,
      include_probabilities = include_probabilities
    ),
    digits = digits
  )
  written_path <- mfrm_write_delimited_table(tbl, path, format, overwrite = overwrite)
  written <- data.frame(
    Component = "residual_file",
    Format = format,
    Path = written_path,
    stringsAsFactors = FALSE
  )
  out <- list(
    table = tbl,
    summary = data.frame(
      Component = "residual_file",
      Rows = nrow(tbl),
      Columns = ncol(tbl),
      Format = format,
      Path = written_path,
      stringsAsFactors = FALSE
    ),
    written_files = written,
    settings = list(
      format = format,
      path = written_path,
      digits = digits,
      overwrite = isTRUE(overwrite),
      include_probabilities = isTRUE(include_probabilities)
    )
  )
  as_mfrm_bundle(out, "mfrm_residual_file")
}

#' Write standalone subset-connectivity files
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. Supplying it
#'   avoids recomputing subset connectivity.
#' @param path Output file path for the subset summary table.
#' @param node_path Optional output file path for the node-level subset table.
#'   When `NULL` and `include_nodes = TRUE`, a sibling file ending in
#'   `_nodes.csv` or `_nodes.tsv` is created.
#' @param format File format: `"csv"` or `"tsv"`. If omitted, inferred from
#'   `path` when the extension is `.csv` or `.tsv`, otherwise `"csv"`.
#' @param digits Rounding digits for numeric columns.
#' @param overwrite If `FALSE`, existing files are not overwritten.
#' @param include_nodes If `TRUE`, also write the node-level facet/level to
#'   subset membership table.
#'
#' @details
#' Subsets are connected components in the observation design graph. The graph
#' links `Person` and all modeled facet levels that co-occur in an observation.
#' Multiple subsets mean the scale is not fully connected unless external
#' anchoring or a deliberate separate-calibration design justifies it.
#'
#' @return A bundle with `table`, `nodes`, `summary`, `written_files`, and
#'   `settings`.
#' @seealso [diagnose_mfrm()], [subset_connectivity_report()],
#'   [write_mfrm_residual_file()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' path <- tempfile(fileext = ".csv")
#' out <- write_mfrm_subset_file(fit, diag, path, overwrite = TRUE)
#' out$written_files
#' }
#' @export
write_mfrm_subset_file <- function(fit,
                                   diagnostics = NULL,
                                   path,
                                   node_path = NULL,
                                   format = c("csv", "tsv"),
                                   digits = 4,
                                   overwrite = FALSE,
                                   include_nodes = TRUE) {
  if (missing(path)) {
    stop("`path` must be supplied.")
  }
  path <- mfrm_validate_output_path(path)
  format <- mfrm_resolve_file_format(if (missing(format)) NULL else format, path)
  digits <- max(0L, as.integer(digits))
  subset_tbls <- mfrm_subset_tables_for_export(fit = fit, diagnostics = diagnostics)
  tbl <- round_numeric_df(subset_tbls$summary, digits = digits)
  nodes <- round_numeric_df(subset_tbls$nodes, digits = digits)

  resolved_node_path <- NA_character_
  if (isTRUE(include_nodes)) {
    if (is.null(node_path)) {
      node_path <- mfrm_default_node_path(path, format)
    }
    node_path <- mfrm_validate_output_path(node_path)
    if (identical(
      normalizePath(path, winslash = "/", mustWork = FALSE),
      normalizePath(node_path, winslash = "/", mustWork = FALSE)
    )) {
      stop("`node_path` must differ from `path` when `include_nodes = TRUE`.")
    }
    mfrm_assert_output_writable(node_path, overwrite = overwrite)
  }
  mfrm_assert_output_writable(path, overwrite = overwrite)

  summary_path <- mfrm_write_delimited_table(tbl, path, format, overwrite = overwrite)

  written <- data.frame(
    Component = "subset_summary",
    Format = format,
    Path = summary_path,
    stringsAsFactors = FALSE
  )
  if (isTRUE(include_nodes)) {
    resolved_node_path <- mfrm_write_delimited_table(nodes, node_path, format, overwrite = overwrite)
    written <- rbind(
      written,
      data.frame(
        Component = "subset_nodes",
        Format = format,
        Path = resolved_node_path,
        stringsAsFactors = FALSE
      )
    )
  }

  out <- list(
    table = tbl,
    nodes = nodes,
    summary = data.frame(
      Component = c("subset_summary", if (isTRUE(include_nodes)) "subset_nodes" else character(0)),
      Rows = c(nrow(tbl), if (isTRUE(include_nodes)) nrow(nodes) else integer(0)),
      Columns = c(ncol(tbl), if (isTRUE(include_nodes)) ncol(nodes) else integer(0)),
      Format = format,
      Path = c(summary_path, if (isTRUE(include_nodes)) resolved_node_path else character(0)),
      stringsAsFactors = FALSE
    ),
    written_files = written,
    settings = list(
      format = format,
      path = summary_path,
      node_path = if (isTRUE(include_nodes)) resolved_node_path else NA_character_,
      digits = digits,
      overwrite = isTRUE(overwrite),
      include_nodes = isTRUE(include_nodes)
    )
  )
  as_mfrm_bundle(out, "mfrm_subset_file")
}

extract_pca_eigenvalues <- function(pca_bundle) {
  if (is.null(pca_bundle)) return(numeric(0))

  eig <- numeric(0)
  if (!is.null(pca_bundle$pca) && "values" %in% names(pca_bundle$pca)) {
    eig <- suppressWarnings(as.numeric(pca_bundle$pca$values))
  }
  if (length(eig) == 0 && !is.null(pca_bundle$cor_matrix)) {
    eig <- tryCatch(
      suppressWarnings(as.numeric(eigen(pca_bundle$cor_matrix, symmetric = TRUE, only.values = TRUE)$values)),
      error = function(e) numeric(0)
    )
  }

  eig[is.finite(eig)]
}

build_pca_variance_table <- function(pca_bundle, facet = NULL) {
  eig <- extract_pca_eigenvalues(pca_bundle)
  if (length(eig) == 0) return(data.frame())

  total <- sum(eig, na.rm = TRUE)
  prop <- if (is.finite(total) && total > 0) eig / total else rep(NA_real_, length(eig))
  out <- data.frame(
    Component = seq_along(eig),
    Eigenvalue = eig,
    Proportion = prop,
    Cumulative = cumsum(prop),
    stringsAsFactors = FALSE
  )
  parallel_tbl <- build_pca_parallel_table(pca_bundle)
  if (nrow(parallel_tbl) > 0 && "Component" %in% names(parallel_tbl)) {
    parallel_cols <- setdiff(names(parallel_tbl), c("Eigenvalue", "Facet"))
    out <- dplyr::left_join(
      out,
      parallel_tbl[, parallel_cols, drop = FALSE],
      by = "Component"
    )
  }
  if (!is.null(facet)) out$Facet <- facet
  out
}

build_pca_parallel_table <- function(pca_bundle, facet = NULL) {
  if (is.null(pca_bundle) || is.null(pca_bundle$parallel) || is.null(pca_bundle$parallel$table)) {
    return(data.frame())
  }
  out <- as.data.frame(pca_bundle$parallel$table, stringsAsFactors = FALSE)
  if (nrow(out) == 0) return(out)
  if (!is.null(facet)) out$Facet <- facet
  out
}

build_residual_parallel_settings <- function(enabled,
                                             reps,
                                             quantile,
                                             method,
                                             seed) {
  data.frame(
    Enabled = isTRUE(enabled),
    Method = if (isTRUE(enabled)) as.character(method) else NA_character_,
    Reps = if (isTRUE(enabled)) as.integer(reps) else NA_integer_,
    Quantile = if (isTRUE(enabled)) as.numeric(quantile) else NA_real_,
    Seed = if (isTRUE(enabled) && !is.null(seed)) as.integer(seed) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

validate_residual_parallel_args <- function(parallel,
                                            parallel_reps,
                                            parallel_quantile,
                                            parallel_method,
                                            seed) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.")
  }
  parallel_reps <- as.integer(parallel_reps)
  if (length(parallel_reps) != 1L || is.na(parallel_reps) || parallel_reps < 1L) {
    stop("`parallel_reps` must be a positive integer.")
  }
  parallel_quantile <- as.numeric(parallel_quantile)
  if (length(parallel_quantile) != 1L || !is.finite(parallel_quantile) ||
      parallel_quantile <= 0 || parallel_quantile >= 1) {
    stop("`parallel_quantile` must be a number between 0 and 1.")
  }
  parallel_method <- match.arg(tolower(as.character(parallel_method[1])), "residual_permutation")
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be NULL or a single integer-like value.")
    }
  }
  list(
    parallel = parallel,
    parallel_reps = parallel_reps,
    parallel_quantile = parallel_quantile,
    parallel_method = parallel_method,
    seed = seed
  )
}

seed_for_residual_parallel <- function(seed, index = 0L) {
  if (is.null(seed)) return(NULL)
  seed <- as.integer(seed)
  index <- as.integer(index)
  if (index <= 0L) return(seed)
  as.integer(((seed + index - 1L) %% .Machine$integer.max) + 1L)
}

residual_parallel_eigenvalues <- function(residual_matrix) {
  mat <- as.matrix(residual_matrix)
  storage.mode(mat) <- "double"
  if (nrow(mat) < 2L || ncol(mat) < 2L) return(numeric(0))
  finite_cols <- colSums(is.finite(mat)) >= 2L
  mat <- mat[, finite_cols, drop = FALSE]
  if (nrow(mat) < 2L || ncol(mat) < 2L) return(numeric(0))

  cor_mat <- tryCatch(
    suppressWarnings(stats::cor(mat, use = "pairwise.complete.obs")),
    error = function(e) NULL
  )
  if (is.null(cor_mat) || nrow(cor_mat) < 2L || ncol(cor_mat) < 2L) {
    return(numeric(0))
  }
  cor_mat[!is.finite(cor_mat)] <- 0
  diag(cor_mat) <- 1
  cor_mat <- ensure_positive_definite(cor_mat)
  eig <- tryCatch(
    suppressWarnings(as.numeric(eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values)),
    error = function(e) numeric(0)
  )
  eig[is.finite(eig)]
}

permute_residual_columns <- function(residual_matrix) {
  out <- as.matrix(residual_matrix)
  storage.mode(out) <- "double"
  for (j in seq_len(ncol(out))) {
    idx <- which(is.finite(out[, j]))
    if (length(idx) > 1L) {
      out[idx, j] <- sample(out[idx, j], size = length(idx), replace = FALSE)
    }
  }
  out
}

compute_residual_parallel_analysis <- function(residual_matrix,
                                               observed_eigenvalues = NULL,
                                               reps = 200L,
                                               quantile = 0.95,
                                               method = "residual_permutation",
                                               seed = NULL) {
  settings <- list(
    method = method,
    reps = as.integer(reps),
    quantile = as.numeric(quantile),
    seed = if (is.null(seed)) NA_integer_ else as.integer(seed)
  )
  empty <- function(error = NULL, warning = NULL) {
    list(
      table = data.frame(),
      successful_reps = 0L,
      error = error,
      warning = warning,
      settings = settings
    )
  }
  if (is.null(residual_matrix)) {
    return(empty("Residual matrix is unavailable."))
  }

  mat <- tryCatch({
    tmp <- as.matrix(residual_matrix)
    storage.mode(tmp) <- "double"
    tmp
  }, error = function(e) NULL)
  if (is.null(mat) || nrow(mat) < 2L || ncol(mat) < 2L) {
    return(empty("Residual matrix must have at least two rows and two columns."))
  }

  keep_cols <- colSums(is.finite(mat)) >= 2L
  mat <- mat[, keep_cols, drop = FALSE]
  if (nrow(mat) < 2L || ncol(mat) < 2L) {
    return(empty("Residual matrix has fewer than two usable columns."))
  }

  observed <- if (is.null(observed_eigenvalues)) residual_parallel_eigenvalues(mat) else as.numeric(observed_eigenvalues)
  observed <- observed[is.finite(observed)]
  if (length(observed) == 0) {
    return(empty("Observed residual eigenvalues are unavailable."))
  }

  if (!is.null(seed)) {
    seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (seed_exists) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
      if (seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  eig_list <- vector("list", as.integer(reps))
  ok <- logical(as.integer(reps))
  for (i in seq_len(as.integer(reps))) {
    eig <- residual_parallel_eigenvalues(permute_residual_columns(mat))
    if (length(eig) > 0L) {
      eig_list[[i]] <- eig
      ok[i] <- TRUE
    }
  }

  eig_list <- eig_list[ok]
  successful <- length(eig_list)
  if (successful == 0L) {
    return(empty("No valid permuted residual correlation matrices were produced."))
  }

  n_components <- length(observed)
  null_mat <- matrix(NA_real_, nrow = successful, ncol = n_components)
  for (i in seq_len(successful)) {
    k <- min(length(eig_list[[i]]), n_components)
    null_mat[i, seq_len(k)] <- eig_list[[i]][seq_len(k)]
  }

  null_mean <- colMeans(null_mat, na.rm = TRUE)
  null_sd <- apply(null_mat, 2L, stats::sd, na.rm = TRUE)
  null_cutoff <- apply(
    null_mat,
    2L,
    stats::quantile,
    probs = as.numeric(quantile),
    na.rm = TRUE,
    names = FALSE,
    type = 8L
  )
  excess <- observed - null_cutoff
  out <- data.frame(
    Component = seq_len(n_components),
    Eigenvalue = observed,
    ParallelMean = null_mean,
    ParallelSD = null_sd,
    ParallelCutoff = null_cutoff,
    ParallelQuantile = rep(as.numeric(quantile), n_components),
    ExcessOverParallelCutoff = excess,
    ExceedsParallelCutoff = is.finite(excess) & excess > 0,
    ParallelReps = rep(as.integer(reps), n_components),
    SuccessfulParallelReps = rep(as.integer(successful), n_components),
    ParallelMethod = rep(method, n_components),
    stringsAsFactors = FALSE
  )
  list(
    table = out,
    successful_reps = as.integer(successful),
    error = NULL,
    warning = if (successful < as.integer(reps)) {
      paste0("Only ", successful, " of ", as.integer(reps), " residual permutations produced valid eigenvalues.")
    } else {
      NULL
    },
    settings = settings
  )
}

attach_residual_parallel <- function(pca_bundle,
                                     reps,
                                     quantile,
                                     method,
                                     seed = NULL) {
  if (!is.list(pca_bundle)) return(pca_bundle)
  bundle_error <- as.character(pca_bundle$error %||% "")
  if (nzchar(bundle_error)) {
    pca_bundle$parallel <- list(
      table = data.frame(),
      successful_reps = 0L,
      error = paste0("PCA result unavailable: ", bundle_error),
      warning = NULL,
      settings = list(
        method = method,
        reps = as.integer(reps),
        quantile = as.numeric(quantile),
        seed = if (is.null(seed)) NA_integer_ else as.integer(seed)
      )
    )
    return(pca_bundle)
  }
  pca_bundle$parallel <- compute_residual_parallel_analysis(
    residual_matrix = pca_bundle$residual_matrix,
    observed_eigenvalues = extract_pca_eigenvalues(pca_bundle),
    reps = reps,
    quantile = quantile,
    method = method,
    seed = seed
  )
  pca_bundle
}

build_residual_parallel_status <- function(out_overall, out_by_facet = list()) {
  one <- function(bundle, scope, facet = NA_character_) {
    pa <- if (is.list(bundle)) bundle$parallel else NULL
    data.frame(
      Scope = scope,
      Facet = facet,
      ParallelAvailable = is.list(pa) && !is.null(pa$table) && nrow(pa$table) > 0L,
      SuccessfulParallelReps = if (is.list(pa)) as.integer(pa$successful_reps %||% 0L) else NA_integer_,
      Error = if (is.list(pa)) as.character(pa$error %||% "") else "",
      Warning = if (is.list(pa)) as.character(pa$warning %||% "") else "",
      stringsAsFactors = FALSE
    )
  }
  rows <- list()
  if (is.list(out_overall)) rows[[length(rows) + 1L]] <- one(out_overall, "overall")
  if (length(out_by_facet) > 0L) {
    rows <- c(rows, lapply(names(out_by_facet), function(f) {
      one(out_by_facet[[f]], "facet", f)
    }))
  }
  if (length(rows) == 0L) {
    return(data.frame(
      Scope = character(0),
      Facet = character(0),
      ParallelAvailable = logical(0),
      SuccessfulParallelReps = integer(0),
      Error = character(0),
      Warning = character(0),
      stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(rows)
}

infer_facet_names <- function(diagnostics) {
  if (!is.null(diagnostics$facet_names) && length(diagnostics$facet_names) > 0) {
    return(unique(as.character(diagnostics$facet_names)))
  }

  if (!is.null(diagnostics$measures) && "Facet" %in% names(diagnostics$measures)) {
    f <- unique(as.character(diagnostics$measures$Facet))
    f <- setdiff(f, "Person")
    if (length(f) > 0) return(f)
  }

  if (!is.null(diagnostics$obs)) {
    nm <- names(diagnostics$obs)
    skip <- c(
      "Person", "Score", "Weight", "score_k", "PersonMeasure",
      "Observed", "Expected", "Var", "Residual", "StdResidual", "StdSq"
    )
    f <- setdiff(nm, skip)
    if (length(f) > 0) return(f)
  }

  character(0)
}

#' Run exploratory residual PCA summaries
#'
#' Legacy-compatible residual diagnostics can be inspected in two ways:
#' 1) overall residual PCA on the person x combined-facet matrix
#' 2) facet-specific residual PCA on person x facet-level matrices
#'
#' @param diagnostics Output from [diagnose_mfrm()] or [fit_mfrm()].
#' @param mode `"overall"`, `"facet"`, or `"both"`.
#' @param facets Optional subset of facets for facet-specific PCA.
#' @param pca_max_factors Maximum number of retained components.
#' @param parallel Logical; if `TRUE`, add residual-permutation parallel
#'   analysis to the PCA tables.
#' @param parallel_reps Number of residual permutations used when
#'   `parallel = TRUE`.
#' @param parallel_quantile Upper null quantile used as the exploratory
#'   comparison cutoff. The default (`0.95`) follows the common parallel
#'   analysis convention.
#' @param parallel_method Parallel-analysis null method. Currently
#'   `"residual_permutation"` is implemented: standardized residuals are
#'   permuted within each residual column, preserving each column's residual
#'   distribution and missingness pattern while breaking residual association.
#' @param seed Optional integer seed for reproducible residual permutations.
#'
#' @details
#' The function works on standardized residual structures derived from
#' [diagnose_mfrm()]. When a fitted object from [fit_mfrm()] is supplied,
#' diagnostics are computed internally.
#'
#' Conceptually, this follows the Rasch residual-PCA tradition of examining
#' structure in model residuals after the primary Rasch dimension has been
#' extracted. In `mfrmr`, however, the implementation is an **exploratory
#' many-facet adaptation**: it works on standardized residual matrices built as
#' person x combined-facet or person x facet-level layouts, rather than
#' reproducing FACETS/Winsteps residual-contrast tables one-to-one.
#'
#' Residual PCA should therefore be reported as residual-structure evidence,
#' not as a formal proof of unidimensionality. It also should not be described
#' as DIMTEST or UNIDIM: those essential-unidimensionality tests require a
#' separate item-response-layer definition that is not uniquely determined by a
#' many-facet long data set. In applied MFRM reporting, residual PCA is best
#' triangulated with global residual fit, element fit, and Q3-style
#' local-dependence screens.
#'
#' Output tables use:
#' - `Component`: principal-component index (1, 2, ...)
#' - `Eigenvalue`: eigenvalue for each component
#' - `Proportion`: component variance proportion
#' - `Cumulative`: cumulative variance proportion
#'
#' When `parallel = TRUE`, the variance tables additionally include
#' data-driven null summaries:
#' - `ParallelMean`: mean permuted-residual eigenvalue
#' - `ParallelCutoff`: `parallel_quantile` cutoff of permuted eigenvalues
#' - `ExcessOverParallelCutoff`: observed eigenvalue minus the cutoff
#' - `ExceedsParallelCutoff`: whether the observed eigenvalue exceeds the
#'   permutation cutoff
#'
#' The default `parallel_reps = 200` is intended as a practical review setting.
#' For stable final reporting of the 95% cutoff, use a larger value when the
#' residual matrix size makes that computationally reasonable.
#'
#' For `mode = "facet"` or `"both"`, `by_facet_table` additionally includes
#' a `Facet` column.
#'
#' `summary(pca)` is supported through `summary()`.
#' `plot(pca)` is dispatched through `plot()` for class
#' `mfrm_residual_pca`. Available types include `"overall_scree"`,
#' `"facet_scree"`, `"overall_parallel_scree"`,
#' `"facet_parallel_scree"`, `"overall_parallel_excess"`,
#' `"facet_parallel_excess"`, `"overall_loadings"`, and
#' `"facet_loadings"`.
#'
#' @section Interpreting output:
#' Use `overall_table` first:
#' - early components with noticeably larger eigenvalues or proportions
#'   suggest stronger residual structure that may deserve follow-up. Small
#'   early components can be described as evidence consistent with the specified
#'   one-dimensional facet structure only when fit and local-dependence screens
#'   tell the same story.
#'
#' Then inspect `by_facet_table`:
#' - helps localize which facet contributes most to residual structure.
#'
#' Finally, inspect loadings via [plot_residual_pca()] to identify which
#' variables/elements drive each component.
#'
#' @section References:
#' The residual-PCA idea follows the Rasch residual-structure literature,
#' especially Linacre's discussions of principal components of Rasch residuals.
#' The current `mfrmr` implementation should be interpreted as an exploratory
#' extension for many-facet workflows rather than as a direct reproduction of a
#' single FACETS/Winsteps output table.
#'
#' The optional parallel analysis follows Horn's data-driven eigenvalue
#' comparison logic and later recommendations to compare observed eigenvalues
#' with high quantiles of an empirical null distribution. Because `mfrmr`
#' applies it to standardized Rasch-family residual matrices, the null
#' distribution is generated by within-column residual permutation rather than
#' by simulating raw item scores.
#'
#' - Horn, J. L. (1965). A rationale and test for the number of factors in
#'   factor analysis. *Psychometrika*, 30, 179-185.
#' - Glorfeld, L. W. (1995). An improvement on Horn's parallel analysis
#'   methodology for selecting the correct number of factors to retain.
#'   *Educational and Psychological Measurement*, 55, 377-393.
#' - Hayton, J. C., Allen, D. G., & Scarpello, V. (2004). Factor retention
#'   decisions in exploratory factor analysis: A tutorial on parallel analysis.
#'   *Organizational Research Methods*, 7, 191-205.
#' - Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality assessment
#'   of ordered polytomous items with parallel analysis. *Psychological
#'   Methods*, 16, 209-220.
#' - Linacre, J. M. (1998). *Structure in Rasch residuals: Why principal
#'   components analysis (PCA)?* Rasch Measurement Transactions, 12(2), 636.
#' - Linacre, J. M. (1998). *Detecting multidimensionality: Which residual
#'   data-type works best?* Journal of Outcome Measurement, 2(3), 266-283.
#' - Eckes, T. (2005). Examining rater effects in TestDaF writing and speaking
#'   performance assessments: A many-facet Rasch analysis. *Language Assessment
#'   Quarterly*, 2(3), 197-221.
#' - Yamashita, T. (2024). An application of many-facet Rasch measurement to
#'   evaluate automated essay scoring: A case of ChatGPT-4.0. *Research Methods
#'   in Applied Linguistics*, 3(3), 100133.
#' - Uto, M. (2021). A multidimensional generalized many-facet Rasch model for
#'   rubric-based performance assessment. *Behaviormetrika*, 48(2), 425-457.
#' - Aryadoust, V., Ng, L. Y., & Sayama, H. (2021). A comprehensive review of
#'   Rasch measurement in language assessment: Recommendations and guidelines
#'   for research. *Language Testing*, 38(1), 6-40.
#' - Tseng, W.-T. (2016). Measuring English vocabulary size via computerized
#'   adaptive testing. *Computers & Education*, 97, 69-85.
#'
#' @section Typical workflow:
#' 1. Fit model and run [diagnose_mfrm()] with `residual_pca = "none"` or `"both"`.
#' 2. Call `analyze_residual_pca(..., mode = "both")`.
#' 3. Review `summary(pca)`, then plot scree/loadings.
#' 4. Cross-check with fit/misfit diagnostics before conclusions.
#'
#' @return
#' A named list with:
#' - `mode`: resolved mode used for computation
#' - `facet_names`: facets analyzed
#' - `overall`: overall PCA bundle (or `NULL`)
#' - `by_facet`: named list of facet PCA bundles
#' - `overall_table`: variance table for overall PCA
#' - `by_facet_table`: stacked variance table across facets
#' - `parallel_settings`, `parallel_overall_table`,
#'   `parallel_by_facet_table`, and `parallel_status`: returned for every call;
#'   the parallel tables are populated when `parallel = TRUE`
#' - `errors`: named list of any per-facet PCA errors that were
#'   caught and turned into `NA_real_` rows in the variance tables
#'   (e.g., `psych::principal()` failure on a near-singular residual
#'   matrix). The list is empty when every facet PCA succeeded.
#' - `warnings`: named list of non-fatal PCA warnings captured from the
#'   underlying PCA engine. These indicate exploratory boundary conditions,
#'   not confirmatory evidence.
#'
#' @seealso [diagnose_mfrm()], [plot_residual_pca()], [mfrmr_visual_diagnostics]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' pca <- analyze_residual_pca(diag, mode = "both")
#' pca2 <- analyze_residual_pca(fit, mode = "both")
#' summary(pca)
#' p <- plot_residual_pca(pca, mode = "overall", plot_type = "scree", draw = FALSE)
#' p$data$plot
#' head(p$data)
#' pca_pa <- analyze_residual_pca(diag, mode = "overall", parallel = TRUE, parallel_reps = 10)
#' head(pca_pa$overall_table)
#' head(pca$overall_table)
#' }
#' @export
analyze_residual_pca <- function(diagnostics,
                                 mode = c("overall", "facet", "both"),
                                 facets = NULL,
                                 pca_max_factors = 10L,
                                 parallel = FALSE,
                                 parallel_reps = 200L,
                                 parallel_quantile = 0.95,
                                 parallel_method = c("residual_permutation"),
                                 seed = NULL) {
  mode <- match.arg(tolower(mode), c("overall", "facet", "both"))
  parallel_args <- validate_residual_parallel_args(
    parallel = parallel,
    parallel_reps = parallel_reps,
    parallel_quantile = parallel_quantile,
    parallel_method = parallel_method,
    seed = seed
  )
  parallel <- parallel_args$parallel
  parallel_reps <- parallel_args$parallel_reps
  parallel_quantile <- parallel_args$parallel_quantile
  parallel_method <- parallel_args$parallel_method
  seed <- parallel_args$seed
  # `"auto"` defers the cap to the per-matrix rank (min(10, ncol-1,
  # nrow-1) inside compute_pca_*). Previously a non-integer value was
  # silently coerced to NA, which broke psych::principal() downstream.
  pca_max_factors <- .resolve_pca_max_factors(pca_max_factors)

  if (inherits(diagnostics, "mfrm_fit")) {
    diagnostics <- diagnose_mfrm(
      diagnostics,
      residual_pca = "none",
      pca_max_factors = pca_max_factors
    )
  }

  if (!is.list(diagnostics)) {
    stop("`diagnostics` must be output from diagnose_mfrm() or fit_mfrm().")
  }

  if (is.null(diagnostics$obs) || nrow(diagnostics$obs) == 0) {
    stop("diagnostics$obs is empty. Run diagnose_mfrm() first.")
  }

  facet_names <- infer_facet_names(diagnostics)
  if (!is.null(facets)) {
    facet_names <- intersect(facet_names, as.character(facets))
    if (length(facet_names) == 0) {
      stop("No matching facets found in diagnostics for `facets`.")
    }
  }

  out_overall <- NULL
  out_by_facet <- list()

  if (mode %in% c("overall", "both")) {
    can_reuse <- is.null(facets) && !is.null(diagnostics$residual_pca_overall)
    out_overall <- if (can_reuse) {
      diagnostics$residual_pca_overall
    } else {
      compute_pca_overall(
        obs_df = diagnostics$obs,
        facet_names = facet_names,
        max_factors = pca_max_factors
      )
    }
  }

  if (mode %in% c("facet", "both")) {
    can_reuse <- is.null(facets) && !is.null(diagnostics$residual_pca_by_facet)
    out_by_facet <- if (can_reuse) {
      diagnostics$residual_pca_by_facet
    } else {
      compute_pca_by_facet(
        obs_df = diagnostics$obs,
        facet_names = facet_names,
        max_factors = pca_max_factors
      )
    }

    if (!is.null(facets) && length(out_by_facet) > 0) {
      out_by_facet <- out_by_facet[intersect(names(out_by_facet), facet_names)]
    }
  }

  if (isTRUE(parallel)) {
    if (is.list(out_overall)) {
      out_overall <- attach_residual_parallel(
        out_overall,
        reps = parallel_reps,
        quantile = parallel_quantile,
        method = parallel_method,
        seed = seed_for_residual_parallel(seed, 0L)
      )
    }
    if (length(out_by_facet) > 0L) {
      out_by_facet <- stats::setNames(
        lapply(seq_along(out_by_facet), function(i) {
          attach_residual_parallel(
            out_by_facet[[i]],
            reps = parallel_reps,
            quantile = parallel_quantile,
            method = parallel_method,
            seed = seed_for_residual_parallel(seed, i)
          )
        }),
        names(out_by_facet)
      )
    }
  }

  overall_table <- build_pca_variance_table(out_overall)
  by_facet_table <- if (length(out_by_facet) == 0) {
    data.frame()
  } else {
    tbls <- lapply(names(out_by_facet), function(f) {
      build_pca_variance_table(out_by_facet[[f]], facet = f)
    })
    tbls <- tbls[vapply(tbls, nrow, integer(1)) > 0]
    if (length(tbls) == 0) data.frame() else dplyr::bind_rows(tbls)
  }
  pca_errors <- list(
    overall = if (is.list(out_overall)) as.character(out_overall$error %||% NA_character_) else NA_character_,
    by_facet = if (length(out_by_facet) == 0) {
      data.frame(Facet = character(0), Error = character(0), stringsAsFactors = FALSE)
    } else {
      err_tbl <- data.frame(
        Facet = names(out_by_facet),
        Error = vapply(out_by_facet, function(x) as.character(x$error %||% ""), character(1)),
        stringsAsFactors = FALSE
      )
      err_tbl[nzchar(err_tbl$Error), , drop = FALSE]
    }
  )
  pca_warnings <- list(
    overall = if (is.list(out_overall)) as.character(out_overall$warning %||% NA_character_) else NA_character_,
    by_facet = if (length(out_by_facet) == 0) {
      data.frame(Facet = character(0), Warning = character(0), stringsAsFactors = FALSE)
    } else {
      warn_tbl <- data.frame(
        Facet = names(out_by_facet),
        Warning = vapply(out_by_facet, function(x) as.character(x$warning %||% ""), character(1)),
        stringsAsFactors = FALSE
      )
      warn_tbl[nzchar(warn_tbl$Warning), , drop = FALSE]
    }
  )

  out <- list(
    mode = mode,
    facet_names = facet_names,
    overall = out_overall,
    by_facet = out_by_facet,
    overall_table = overall_table,
    by_facet_table = by_facet_table,
    parallel_settings = build_residual_parallel_settings(
      enabled = parallel,
      reps = parallel_reps,
      quantile = parallel_quantile,
      method = parallel_method,
      seed = seed
    ),
    parallel_overall_table = build_pca_parallel_table(out_overall),
    parallel_by_facet_table = if (length(out_by_facet) == 0L) {
      data.frame()
    } else {
      parallel_tbls <- lapply(names(out_by_facet), function(f) {
        build_pca_parallel_table(out_by_facet[[f]], facet = f)
      })
      parallel_tbls <- parallel_tbls[vapply(parallel_tbls, nrow, integer(1)) > 0]
      if (length(parallel_tbls) == 0L) data.frame() else dplyr::bind_rows(parallel_tbls)
    },
    parallel_status = build_residual_parallel_status(out_overall, out_by_facet),
    errors = pca_errors,
    warnings = pca_warnings
  )
  as_mfrm_bundle(out, "mfrm_residual_pca")
}

# Resolve pca_max_factors = "auto" to NA_integer_ so downstream
# compute_pca_*() defers the cap to min(10, ncol - 1, nrow - 1) for
# each residual matrix. Any other value is coerced to integer.
.resolve_pca_max_factors <- function(x) {
  if (is.character(x) && length(x) == 1L && tolower(x) == "auto") {
    return(NA_integer_)
  }
  as.integer(x)
}

resolve_pca_input <- function(x) {
  if (is.null(x)) stop("Input cannot be NULL.")
  if (!is.null(x$overall_table) || !is.null(x$by_facet_table)) return(x)
  if (inherits(x, "mfrm_fit")) return(analyze_residual_pca(x, mode = "both"))
  if (!is.null(x$obs)) return(analyze_residual_pca(x, mode = "both"))
  stop("Input must be fit from fit_mfrm(), diagnostics from diagnose_mfrm(), or result from analyze_residual_pca().")
}

extract_loading_table <- function(pca_bundle, component = 1L, top_n = 20L) {
  if (is.null(pca_bundle) || is.null(pca_bundle$pca) || is.null(pca_bundle$pca$loadings)) {
    return(data.frame())
  }

  loads <- tryCatch(as.matrix(unclass(pca_bundle$pca$loadings)), error = function(e) NULL)
  if (is.null(loads) || nrow(loads) == 0) return(data.frame())
  if (component > ncol(loads) || component < 1) return(data.frame())

  vals <- suppressWarnings(as.numeric(loads[, component]))
  vars <- rownames(loads)
  if (is.null(vars)) vars <- paste0("V", seq_along(vals))

  ok <- is.finite(vals)
  if (!any(ok)) return(data.frame())

  tbl <- data.frame(
    Variable = vars[ok],
    Loading = vals[ok],
    stringsAsFactors = FALSE
  )
  tbl <- tbl[order(abs(tbl$Loading), decreasing = TRUE), , drop = FALSE]
  top_n <- max(1L, as.integer(top_n))
  head(tbl, n = min(nrow(tbl), top_n))
}

#' Visualize residual PCA results
#'
#' @param x Output from [analyze_residual_pca()], [diagnose_mfrm()], or [fit_mfrm()].
#' @param mode `"overall"` or `"facet"`.
#' @param facet Facet name for `mode = "facet"`.
#' @param plot_type `"scree"`, `"parallel_scree"`,
#'   `"parallel_excess"`, or `"loadings"`.
#' @param component Component index for loadings plot.
#' @param top_n Maximum number of variables shown in loadings plot.
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`, or `"monochrome"`).
#' @param draw If `TRUE`, draws the plot using base graphics.
#'
#' @details
#' `x` can be either:
#' - output of [analyze_residual_pca()], or
#' - a diagnostics object from [diagnose_mfrm()] (PCA is computed internally), or
#' - a fitted object from [fit_mfrm()] (diagnostics and PCA are computed internally).
#'
#' Plot types:
#' - `"scree"`: component vs eigenvalue line plot
#' - `"parallel_scree"`: observed eigenvalues with residual-permutation
#'   parallel-analysis mean and upper cutoff
#' - `"parallel_excess"`: observed eigenvalue minus the parallel-analysis
#'   cutoff by component
#' - `"loadings"`: horizontal bar chart of top absolute loadings
#'
#' For `mode = "facet"` and `facet = NULL`, the first available facet is used.
#'
#' @section Interpreting output:
#' - `plot_type = "scree"`: look for dominant early components relative
#'   to later components and the unit-eigenvalue reference line. Treat
#'   this as exploratory residual-structure screening, not a standalone
#'   unidimensionality test or a DIMTEST/UNIDIM substitute.
#' - `plot_type = "parallel_scree"` or `"parallel_excess"`: use only after
#'   running [analyze_residual_pca()] with `parallel = TRUE`. Components
#'   above the residual-permutation cutoff are candidates for follow-up,
#'   not proof of multidimensionality.
#' - `plot_type = "loadings"`: identifies variables/elements driving each
#'   component; inspect both sign and absolute magnitude.
#'
#' Facet mode (`mode = "facet"`) helps localize residual structure to a
#' specific facet after global PCA review.
#'
#' @section Typical workflow:
#' 1. Run [diagnose_mfrm()] with `residual_pca = "overall"` or `"both"`.
#' 2. Build PCA object via [analyze_residual_pca()] (or pass diagnostics directly).
#' 3. Use scree plot first, then loadings plot for targeted interpretation.
#'
#' @return
#' A named list of plotting data (class `mfrm_plot_data`) with:
#' - `plot`: `"scree"`, `"parallel_scree"`, `"parallel_excess"`, or
#'   `"loadings"`
#' - `mode`: `"overall"` or `"facet"`
#' - `facet`: facet name (or `NULL`)
#' - `title`: plot title text
#' - `data`: underlying table used for plotting
#'
#' @seealso [analyze_residual_pca()], [diagnose_mfrm()]
#' @examples
#' toy_full <- load_mfrmr_data("example_core")
#' toy_people <- unique(toy_full$Person)[1:24]
#' toy <- toy_full[match(toy_full$Person, toy_people, nomatch = 0L) > 0L, , drop = FALSE]
#' fit <- suppressWarnings(
#'   fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' )
#' diag <- diagnose_mfrm(fit, residual_pca = "overall")
#' pca <- analyze_residual_pca(diag, mode = "overall")
#' plt <- plot_residual_pca(pca, mode = "overall", plot_type = "scree", draw = FALSE)
#' head(plt$data)
#' \donttest{
#' pca_pa <- analyze_residual_pca(diag, mode = "overall", parallel = TRUE, parallel_reps = 10)
#' pa <- plot_residual_pca(pca_pa, mode = "overall", plot_type = "parallel_scree", draw = FALSE)
#' head(pa$data)
#' }
#' plt_load <- plot_residual_pca(
#'   pca, mode = "overall", plot_type = "loadings", component = 1, draw = FALSE
#' )
#' head(plt_load$data)
#' if (interactive()) {
#'   plot_residual_pca(pca, mode = "overall", plot_type = "scree", preset = "publication")
#' }
#' @export
plot_residual_pca <- function(x,
                              mode = c("overall", "facet"),
                              facet = NULL,
                              plot_type = c("scree", "parallel_scree", "parallel_excess", "loadings"),
                              component = 1L,
                              top_n = 20L,
                              preset = c("standard", "publication", "compact", "monochrome"),
                              draw = TRUE) {
  mode <- match.arg(tolower(mode), c("overall", "facet"))
  plot_type <- match.arg(tolower(plot_type), c("scree", "parallel_scree", "parallel_excess", "loadings"))
  style <- resolve_plot_preset(preset)
  pca_obj <- resolve_pca_input(x)

  pca_bundle <- NULL
  title_suffix <- ""

  if (mode == "overall") {
    pca_bundle <- pca_obj$overall
    title_suffix <- "Overall Residual PCA"
  } else {
    if (is.null(pca_obj$by_facet) || length(pca_obj$by_facet) == 0) {
      stop("No facet-level PCA results available.")
    }
    if (is.null(facet)) facet <- names(pca_obj$by_facet)[1]
    if (!facet %in% names(pca_obj$by_facet)) {
      stop("Requested facet not found in PCA results.")
    }
    pca_bundle <- pca_obj$by_facet[[facet]]
    title_suffix <- paste0("Residual PCA - ", facet)
  }

  if (plot_type %in% c("scree", "parallel_scree", "parallel_excess")) {
    tbl <- build_pca_variance_table(pca_bundle)
    if (nrow(tbl) == 0) stop("No eigenvalues available for scree plot.")
    has_parallel <- all(c("ParallelMean", "ParallelCutoff", "ExcessOverParallelCutoff",
                          "ExceedsParallelCutoff") %in% names(tbl))
    if (plot_type %in% c("parallel_scree", "parallel_excess") && !has_parallel) {
      stop("Parallel-analysis results are unavailable. Run analyze_residual_pca(..., parallel = TRUE) first.")
    }
    title <- paste0(title_suffix, " (Scree)")
    subtitle <- if (mode == "overall") {
      "Variance explained by residual components"
    } else {
      paste0("Facet-specific scree profile: ", facet)
    }
    # Rasch-conventional secondary-dimension reference bands on the residual
    # eigenvalue scale (see Linacre, 2026, A User's Guide to Winsteps 5.11.0):
    # 1.4 critical minimum, 2.0 noticeable, 3.0 strong second dimension.
    rasch_refs <- c(1, 1.4, 2, 3)
    rasch_ref_labels <- c(
      "Unit-eigenvalue",
      "Critical minimum (1.4)",
      "Noticeable second dim (2.0)",
      "Strong second dim (3.0)"
    )

    if (plot_type == "parallel_scree") {
      title <- paste0(title_suffix, " (Parallel Scree)")
      q_percent <- formatC(100 * unique(tbl$ParallelQuantile)[1], format = "fg", digits = 4)
      q_label <- paste0(q_percent, "% residual-permutation cutoff")
      subtitle <- if (mode == "overall") {
        "Observed residual eigenvalues compared with the permutation null"
      } else {
        paste0("Facet-specific permutation comparison: ", facet)
      }
      if (isTRUE(draw)) {
        apply_plot_preset(style)
        yr <- range(c(tbl$Eigenvalue, tbl$ParallelMean, tbl$ParallelCutoff, rasch_refs), finite = TRUE)
        graphics::plot(
          x = tbl$Component,
          y = tbl$Eigenvalue,
          type = "b",
          pch = 16,
          col = style$accent_primary,
          xlab = "Component",
          ylab = "Eigenvalue",
          ylim = yr,
          main = title
        )
        graphics::abline(h = pretty(graphics::par("usr")[3:4], n = 5), col = style$grid, lty = 1)
        graphics::abline(v = pretty(graphics::par("usr")[1:2], n = 5), col = style$grid, lty = 1)
        graphics::lines(
          x = tbl$Component,
          y = tbl$ParallelCutoff,
          type = "b",
          pch = 17,
          lty = 2,
          col = style$accent_secondary
        )
        graphics::lines(
          x = tbl$Component,
          y = tbl$ParallelMean,
          type = "l",
          lty = 3,
          col = style$accent_tertiary
        )
        graphics::abline(h = rasch_refs, lty = c(2, 3, 3, 3), col = style$neutral)
      }
      out <- new_mfrm_plot_data(
        "residual_pca",
        list(
          plot = "parallel_scree",
          mode = mode,
          facet = if (mode == "facet") facet else NULL,
          title = title,
          subtitle = subtitle,
          legend = new_plot_legend(
            label = c("Observed residual eigenvalues", q_label, "Parallel mean", rasch_ref_labels),
            role = c("component", "parallel_cutoff", "parallel_mean", rep("reference", length(rasch_refs))),
            aesthetic = c("line-point", "line-point", "line", rep("line", length(rasch_refs))),
            value = c(style$accent_primary, style$accent_secondary, style$accent_tertiary,
                      rep(style$neutral, length(rasch_refs)))
          ),
          reference_lines = new_reference_lines("h", rasch_refs, rasch_ref_labels,
                                                c("dashed", "dotted", "dotted", "dotted"),
                                                rep("reference", length(rasch_refs))),
          data = tbl,
          preset = style$name
        )
      )
      return(invisible(out))
    }

    if (plot_type == "parallel_excess") {
      title <- paste0(title_suffix, " (Parallel Excess)")
      subtitle <- "Observed residual eigenvalue minus permutation cutoff"
      if (isTRUE(draw)) {
        apply_plot_preset(style)
        cols <- ifelse(tbl$ExceedsParallelCutoff, style$accent_secondary, style$neutral)
        graphics::barplot(
          height = tbl$ExcessOverParallelCutoff,
          names.arg = tbl$Component,
          col = cols,
          border = style$background,
          xlab = "Component",
          ylab = "Eigenvalue minus cutoff",
          main = title
        )
        graphics::abline(h = 0, lty = 2, col = style$neutral)
      }
      out <- new_mfrm_plot_data(
        "residual_pca",
        list(
          plot = "parallel_excess",
          mode = mode,
          facet = if (mode == "facet") facet else NULL,
          title = title,
          subtitle = subtitle,
          legend = new_plot_legend(
            label = c("Above cutoff", "At/below cutoff"),
            role = c("parallel_excess", "parallel_excess"),
            aesthetic = c("bar", "bar"),
            value = c(style$accent_secondary, style$neutral)
          ),
          reference_lines = new_reference_lines("h", 0, "Permutation cutoff", "dashed", "reference"),
          data = tbl,
          preset = style$name
        )
      )
      return(invisible(out))
    }

    if (isTRUE(draw)) {
      apply_plot_preset(style)
      graphics::plot(
        x = tbl$Component,
        y = tbl$Eigenvalue,
        type = "b",
        pch = 16,
        col = style$accent_primary,
        xlab = "Component",
        ylab = "Eigenvalue",
        main = title
      )
      graphics::abline(h = pretty(graphics::par("usr")[3:4], n = 5), col = style$grid, lty = 1)
      graphics::abline(v = pretty(graphics::par("usr")[1:2], n = 5), col = style$grid, lty = 1)
      graphics::abline(h = rasch_refs, lty = c(2, 3, 3, 3), col = style$neutral)
    }

    out <- new_mfrm_plot_data(
      "residual_pca",
      list(
        plot = "scree",
        mode = mode,
        facet = if (mode == "facet") facet else NULL,
        title = title,
        subtitle = subtitle,
        legend = new_plot_legend(
          label = c("Residual eigenvalues", rasch_ref_labels),
          role = c("component", rep("reference", length(rasch_refs))),
          aesthetic = c("line-point", rep("line", length(rasch_refs))),
          value = c(style$accent_primary, rep(style$neutral, length(rasch_refs)))
        ),
        reference_lines = new_reference_lines("h", rasch_refs, rasch_ref_labels,
                                              c("dashed", "dotted", "dotted", "dotted"),
                                              rep("reference", length(rasch_refs))),
        data = tbl,
        preset = style$name
      )
    )
    return(invisible(out))
  }

  load_tbl <- extract_loading_table(
    pca_bundle = pca_bundle,
    component = as.integer(component),
    top_n = as.integer(top_n)
  )
  if (nrow(load_tbl) == 0) stop("No loadings available for the requested component.")

  load_tbl <- load_tbl[order(load_tbl$Loading), , drop = FALSE]
  title <- paste0(title_suffix, " (Loadings: PC", as.integer(component), ")")
  subtitle <- paste0("Top ", min(nrow(load_tbl), as.integer(top_n)), " absolute loadings")

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    cols <- ifelse(load_tbl$Loading >= 0, style$accent_tertiary, style$accent_secondary)
    graphics::barplot(
      height = load_tbl$Loading,
      names.arg = load_tbl$Variable,
      horiz = TRUE,
      las = 1,
      col = cols,
      border = style$background,
      xlab = "Loading",
      main = title
    )
    graphics::abline(v = 0, lty = 2, col = style$neutral)
  }

  out <- new_mfrm_plot_data(
    "residual_pca",
    list(
      plot = "loadings",
      mode = mode,
      facet = if (mode == "facet") facet else NULL,
      title = title,
      subtitle = subtitle,
      legend = new_plot_legend(
        label = c("Positive loadings", "Negative loadings"),
        role = c("loading", "loading"),
        aesthetic = c("bar", "bar"),
        value = c(style$accent_tertiary, style$accent_secondary)
      ),
      reference_lines = new_reference_lines("v", 0, "Zero-loading reference", "dashed", "reference"),
      component = as.integer(component),
      data = load_tbl,
      preset = style$name
    )
  )
  invisible(out)
}

#' Estimate bias and interaction screening terms
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param facet_a First facet name. Provide together with `facet_b` for the
#'   classic pairwise 2-way interaction. Ignored when `interaction_facets`
#'   is supplied.
#' @param facet_b Second facet name. See `facet_a`.
#' @param interaction_facets Character vector of two or more facets to model as
#'   one interaction effect. When supplied, this takes precedence over
#'   `facet_a`/`facet_b`. Use this form (rather than `facet_a`/`facet_b`)
#'   whenever you want 3+ way interactions, since `facet_a/facet_b` is
#'   restricted to the pairwise case.
#' @param max_abs Bound for absolute bias size.
#' @param omit_extreme Omit extreme-only elements.
#' @param max_iter Iteration cap.
#' @param tol Convergence tolerance.
#'
#' @details
#' **Bias (interaction) in MFRM** refers to a systematic departure from
#' the additive model: a specific rater-criterion (or higher-order)
#' combination produces scores that are consistently higher or lower than
#' predicted by the main effects alone.  For example, Rater A might be
#' unexpectedly harsh on Criterion 2 despite being lenient overall.
#'
#' Mathematically, the bias term \eqn{b_{jc}} for rater \eqn{j} on
#' criterion \eqn{c} modifies the linear predictor:
#'
#' \deqn{\eta_{njc} = \theta_n - \delta_j - \beta_c - b_{jc}}
#'
#' For `RSM` / `PCM`, the function estimates \eqn{b_{jc}} from the residuals
#' of the fitted additive model using an iterative recalibration screen aligned
#' with the many-facet bias literature (Myford & Wolfe, 2003, 2004):
#'
#' \deqn{b_{jc} = \frac{\sum_n (X_{njc} - E_{njc})}
#'                     {\sum_n \mathrm{Var}_{njc}}}
#'
#' Each iteration updates expected scores using the current bias estimates,
#' then re-computes the bias. Convergence is reached when the maximum absolute
#' change in bias estimates falls below `tol`. For bounded `GPCM`, the same
#' additive-bias idea is evaluated with the slope-aware GPCM kernel and
#' conditional profile-likelihood follow-up columns; those quantities remain
#' screening evidence because theta, facet, step, and slope estimates are held
#' fixed.
#'
#' - For two-way mode, use `facet_a` and `facet_b` (or `interaction_facets`
#'   with length 2).
#' - For higher-order mode, provide `interaction_facets` with length >= 3.
#'
#' @section What this screening means:
#' `estimate_bias()` summarizes interaction departures from the additive MFRM.
#' It is best read as a targeted screening tool for potentially noteworthy
#' cells or facet combinations that may merit substantive review.
#'
#' @section What this screening does not justify:
#' - `t` and `Prob.` are screening metrics, not formal inferential quantities.
#' - A flagged interaction cell is not, by itself, proof of rater bias or
#'   construct-irrelevant variance.
#' - Non-flagged cells should not be over-read as evidence that interaction
#'   effects are absent.
#'
#' @section Interpreting output:
#' Use `summary` for global magnitude, then inspect `table` for cell-level
#' interaction effects.
#'
#' Prioritize rows with:
#' - larger `|Bias Size|` (effect on logit scale; \eqn{> 0.5} logits is
#'   typically noteworthy, \eqn{> 1.0} is large)
#' - larger `|t|` among the screening metrics (\eqn{|t| \ge 2} suggests a
#'   screen-positive interaction cell)
#' - smaller `Prob.` among the screening metrics
#'
#' A positive `Obs-Exp Average` means the cell produced *higher* scores
#' than the additive model predicts (unexpected leniency); negative
#' means unexpected harshness.
#'
#' `iteration` helps verify whether iterative recalibration stabilized.
#' If the maximum change on the final iteration is still above `tol`,
#' consider increasing `max_iter`.
#'
#' @section Typical workflow:
#' 1. Fit and diagnose model.
#' 2. Run `estimate_bias(...)` for target interaction facets.
#' 3. Review `summary(bias)` and `bias$table`.
#' 4. Visualize/report via [plot_bias_interaction()] and [build_fixed_reports()].
#'
#' @section Interpreting key output columns:
#' In `bias$table`, the most-used columns are:
#' - `Bias Size`: estimated interaction effect \eqn{b_{jc}} (logit scale)
#' - `t` and `Prob.`: screening metrics, not formal inferential quantities
#' - `Obs-Exp Average`: direction and practical size of observed-vs-expected
#'   gap on the raw-score metric
#' - for bounded `GPCM`, `LR ChiSq`, `LR Prob.`, and `Profile CI Lower` /
#'   `Profile CI Upper`: conditional profile-likelihood checks for a single
#'   additive bias shift, holding the fitted person, facet, step, and slope
#'   estimates fixed
#'
#' The `chi_sq` element provides a fixed-effect heterogeneity screen across all
#' interaction cells.
#'
#' @section Recommended next step:
#' Use [plot_bias_interaction()] to inspect the flagged cells visually, then
#' integrate the result with DFF, linking, or substantive scoring review before
#' making formal claims about fairness or invariance.
#'
#' @return
#' An object of class `mfrm_bias` with:
#' - `table`: interaction rows with effect size, SE, screening t/p metadata,
#'   reporting-use flags, fit columns, and bounded-`GPCM`
#'   profile-likelihood columns when available
#' - `summary`: compact summary statistics
#' - `chi_sq`: fixed-effect chi-square style screening summary
#' - `facet_a`, `facet_b`: first two analyzed facet names (legacy compatibility)
#' - `interaction_facets`, `interaction_order`, `interaction_mode`: full
#'   interaction metadata
#' - `iteration`: iteration history/metadata
#' - `orientation_review`: facet-orientation sign-consistency review table
#' - `mixed_sign`: logical flag indicating whether bias-size signs flip
#'   across facets in a way that complicates direction interpretation
#' - `direction_note`: one-line interpretive note describing the
#'   dominant bias direction (empty when not applicable)
#' - `recommended_action`: one-line recommended-action label routing
#'   the user to the appropriate follow-up helper
#' - `inference_tier`: summary label indicating that the bias rows are
#'   intended for screening and follow-up review in this release
#' - `optimization_failures`: per-cell record of any inner-loop
#'   optimizer failures encountered while estimating the bias
#'   parameters; empty when every cell converged cleanly
#'
#' @seealso [build_fixed_reports()], [build_apa_outputs()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' s_bias <- summary(bias)
#' s_bias$overview
#' # Look for: `MaxAbsBias` < ~0.5 logits and `Significant = 0` mean
#' #   no cell exceeded the screen. The `BonferroniSignificant` /
#' #   `HolmSignificant` columns count cells that survive multiple-
#' #   testing correction; both being 0 is a stronger "no bias"
#' #   signal than the raw screen-positive count alone.
#' s_bias$top_rows
#' # Look for: rows with `|t|` > 2 and |Bias Size| > 0.5 logits warrant
#' #   review (large effect AND statistically reliable). Rows with only
#' #   one of those triggered are usually small-cell artefacts.
#' p_bias <- plot_bias_interaction(bias, draw = FALSE)
#' p_bias$data$plot
#'
#' @section References:
#' - Linacre, J. M. (1989). *Many-Facet Rasch Measurement*. MESA Press.
#'   (FACETS Table 13 corresponds to the bias / interaction
#'   estimation that this helper implements.)
#' - Eckes, T. (2005). Examining rater effects in TestDaF writing
#'   and speaking performance assessments: A many-facet Rasch
#'   analysis. *Language Assessment Quarterly, 2*(3), 197-221.
#' - Eckes, T. (2015). *Introduction to many-facet Rasch
#'   measurement: Analyzing and evaluating rater-mediated
#'   assessments* (2nd ed.). Peter Lang.
#' - Myford, C. M., & Wolfe, E. W. (2003). Detecting and measuring
#'   rater effects using many-facet Rasch measurement: Part I.
#'   *Journal of Applied Measurement, 4*(4), 386-422.
#' - Myford, C. M., & Wolfe, E. W. (2004). Detecting and measuring
#'   rater effects using many-facet Rasch measurement: Part II.
#'   *Journal of Applied Measurement, 5*(2), 189-227.
#' @export
estimate_bias <- function(fit,
                          diagnostics,
                          facet_a = NULL,
                          facet_b = NULL,
                          interaction_facets = NULL,
                          max_abs = 10,
                          omit_extreme = TRUE,
                          max_iter = 4,
                          tol = 1e-3) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm(). ",
         "Got: ", paste(class(fit), collapse = "/"), ".", call. = FALSE)
  }
  fit_model <- as.character(fit$config$model %||% fit$summary$Model[1] %||% NA_character_)
  if (missing(diagnostics)) {
    stop("`diagnostics` is required. Call diagnose_mfrm(fit) first and ",
         "pass the result as the second argument: ",
         "estimate_bias(fit, diagnose_mfrm(fit), facet_a = ..., facet_b = ...).",
         call. = FALSE)
  }
  if (!is.list(diagnostics) || is.null(diagnostics$obs)) {
    stop("`diagnostics` must be the output of diagnose_mfrm(). ",
         "Run: diagnostics <- diagnose_mfrm(fit)", call. = FALSE)
  }
  if (is.null(interaction_facets)) {
    if (is.null(facet_a) || is.null(facet_b)) {
      stop("Provide either `interaction_facets` (length >= 2) ",
           "or both `facet_a` and `facet_b`.", call. = FALSE)
    }
    interaction_facets <- c(as.character(facet_a[1]), as.character(facet_b[1]))
  }
  interaction_facets <- as.character(interaction_facets)
  interaction_facets <- interaction_facets[!is.na(interaction_facets) & nzchar(interaction_facets)]
  interaction_facets <- unique(interaction_facets)
  if (length(interaction_facets) < 2) {
    stop("`interaction_facets` must contain at least two facet names.")
  }
  if (is.null(facet_a) && length(interaction_facets) >= 1) facet_a <- interaction_facets[1]
  if (is.null(facet_b) && length(interaction_facets) >= 2) facet_b <- interaction_facets[2]

  # Validate that every requested facet label actually names a facet in the
  # fit. Without this check, a typo (e.g. "Raters" with trailing s) used to
  # fall through as an empty list, which silently masked the error.
  known_facets <- as.character(fit$config$facet_names %||% character())
  unknown <- setdiff(as.character(interaction_facets), known_facets)
  if (length(unknown) > 0L) {
    stop(
      "`facet_a` / `facet_b` / `interaction_facets` refer to facets ",
      "that are not part of this fit. Unknown: ",
      paste(shQuote(unknown), collapse = ", "),
      ". Available facets: ",
      paste(shQuote(known_facets), collapse = ", "), ".",
      call. = FALSE
    )
  }

  out <- estimate_bias_interaction(
    res = fit,
    diagnostics = diagnostics,
    facet_a = facet_a,
    facet_b = facet_b,
    interaction_facets = interaction_facets,
    max_abs = max_abs,
    omit_extreme = omit_extreme,
    max_iter = max_iter,
    tol = tol
  )
  if (is.list(out) && length(out) > 0) {
    class(out) <- c("mfrm_bias", class(out))
    if (identical(fit_model, "GPCM")) {
      out$method <- "GPCM-slope-aware"
      out$caveat <- paste0(
        "GPCM bias estimates use the slope-aware GPCM kernel: the bias ",
        "parameter is the additive shift on the linear predictor that ",
        "maximises the GPCM log-likelihood for the interaction cell. ",
        "For GPCM rows, `LR ChiSq`, `LR Prob.`, and profile-CI columns ",
        "compare that fitted shift with zero by conditional profile ",
        "likelihood. All reported bias quantities still hold theta, steps, ",
        "slopes, and other facet estimates fixed, so use them for screening ",
        "and follow-up review rather than as standalone fairness claims."
      )
    }
  }
  out
}

#' Build legacy-compatible fixed-width text reports
#'
#' @param bias_results Output from [estimate_bias()].
#' @param target_facet Optional target facet for pairwise contrast table.
#' @param branch Output branch:
#'   `"facets"` keeps the legacy-compatible fixed-width layout;
#'   `"original"` returns compact sectioned fixed-width text for report drafts.
#'
#' @details
#' This function generates plain-text, fixed-width output intended to be read in
#' console/log environments or exported into text reports.
#'
#' The pairwise section (Table 14 style) is only generated for 2-way bias runs.
#' For higher-order interactions (`interaction_facets` length >= 3), the function
#' returns the bias table text and a note explaining why pairwise contrasts were
#' skipped.
#'
#' @section Interpreting output:
#' - `bias_fixed`: fixed-width table of interaction effects.
#' - `pairwise_fixed`: pairwise contrast text (2-way only).
#' - `pairwise_table`: structured contrast table.
#' - `interaction_label`: facets used for the bias run.
#'
#' @section Typical workflow:
#' 1. Run [estimate_bias()].
#' 2. Build text bundle with `build_fixed_reports(...)`.
#' 3. Use `summary()`/`plot()` for quick checks, then export text blocks.
#'
#' @section Preferred route for new analyses:
#' For new reporting workflows, prefer [bias_interaction_report()] and
#' [build_apa_outputs()]. Use `build_fixed_reports()` when a fixed-width text
#' artifact is specifically required for a compatibility handoff.
#'
#' @return
#' A named list with class `mfrm_fixed_reports` (and a branch-specific
#' subclass `mfrm_fixed_reports_<branch>`):
#' - `bias_fixed`: fixed-width interaction table text
#' - `pairwise_fixed`: fixed-width pairwise contrast text
#' - `pairwise_table`: underlying pairwise data.frame
#' - `branch`: character scalar `"original"` or `"facets"` echoing
#'   which fixed-width style was rendered
#' - `style`: character scalar carrying the resolved style preset
#'   used when building the text artifact
#' - `interaction_label`: human-readable label for the interaction
#'   that drove the bias run (`"Rater x Criterion"`-style); `NA`
#'   when no bias rows are available
#' - `target_facet`: character scalar identifying which facet was
#'   used as the target facet for pairwise contrasts; `NA` when no
#'   pairwise contrasts were requested or available
#'
#' @seealso [estimate_bias()], [build_apa_outputs()], [bias_interaction_report()],
#'   [mfrmr_reports_and_tables], [mfrmr_compatibility_layer]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' fixed <- build_fixed_reports(bias)
#' fixed_original <- build_fixed_reports(bias, branch = "original")
#' summary(fixed)
#' p <- plot(fixed, draw = FALSE)
#' p2 <- plot(fixed, type = "pvalue", draw = FALSE)
#' if (interactive()) {
#'   plot(
#'     fixed,
#'     type = "contrast",
#'     draw = TRUE,
#'     main = "Pairwise Contrasts (Customized)",
#'     palette = c(pos = "#1b9e77", neg = "#d95f02"),
#'     label_angle = 45
#'   )
#' }
#' }
#' @export
build_fixed_reports <- function(bias_results,
                                target_facet = NULL,
                                branch = c("facets", "original")) {
  branch <- match.arg(tolower(as.character(branch[1])), c("facets", "original"))
  style <- ifelse(branch == "facets", "facets_manual", "original")

  make_empty_bundle <- function(msg) {
    out <- list(
      bias_fixed = msg,
      pairwise_fixed = "No pairwise data",
      pairwise_table = tibble::tibble(),
      branch = branch,
      style = style,
      interaction_label = NA_character_,
      target_facet = as.character(target_facet %||% NA_character_)
    )
    out <- as_mfrm_bundle(out, "mfrm_fixed_reports")
    class(out) <- unique(c(paste0("mfrm_fixed_reports_", branch), class(out)))
    out
  }

  if (!is.null(bias_results) &&
      !inherits(bias_results, "mfrm_bias") &&
      !(is.list(bias_results) && !is.data.frame(bias_results))) {
    stop(
      "`bias_results` must be NULL, output from estimate_bias(), or a list-like bias bundle with a `table` component.",
      call. = FALSE
    )
  }

  if (is.null(bias_results) || is.null(bias_results$table) || nrow(bias_results$table) == 0) {
    return(make_empty_bundle("No bias data"))
  }

  if (!is.data.frame(bias_results$table)) {
    stop("`bias_results$table` must be a data.frame-like bias table.", call. = FALSE)
  }

  spec <- extract_bias_facet_spec(bias_results)
  if (is.null(spec) || length(spec$facets) < 2) {
    stop(
      "`bias_results` must come from estimate_bias() or another package-native bias helper with recognizable interaction facet columns.",
      call. = FALSE
    )
  }

  facets <- spec$facets
  if (!is.null(target_facet)) {
    target_facet <- as.character(target_facet[1] %||% NA_character_)
    if (!is.na(target_facet) && nzchar(target_facet) && !target_facet %in% facets) {
      stop(
        "`target_facet` must be one of the interaction facets in `bias_results`: ",
        paste(facets, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  interaction_label <- paste(facets, collapse = " x ")
  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)

  core_cols <- c(
    "Sq", "Observd Score", "Expctd Score", "Observd Count", "Obs-Exp Average",
    "Bias Size", "S.E.", "t", "d.f.", "Prob.", "Infit", "Outfit"
  )
  detail_cols <- c(spec$index_cols, spec$level_cols, spec$measure_cols)
  keep_cols <- c(core_cols, detail_cols)
  keep_cols <- keep_cols[keep_cols %in% names(tbl)]
  tbl_display <- tbl |>
    dplyr::select(dplyr::all_of(keep_cols))

  if ("S.E." %in% names(tbl_display)) {
    tbl_display <- dplyr::rename(tbl_display, `Model S.E.` = `S.E.`)
  }
  if ("Infit" %in% names(tbl_display)) {
    tbl_display <- dplyr::rename(tbl_display, `Infit MnSq` = Infit)
  }
  if ("Outfit" %in% names(tbl_display)) {
    tbl_display <- dplyr::rename(tbl_display, `Outfit MnSq` = Outfit)
  }

  for (i in seq_along(facets)) {
    facet_i <- facets[i]
    idx_col <- spec$index_cols[i]
    lvl_col <- spec$level_cols[i]
    meas_col <- spec$measure_cols[i]
    if (idx_col %in% names(tbl_display)) {
      names(tbl_display)[names(tbl_display) == idx_col] <- paste0(facet_i, " N")
    }
    if (lvl_col %in% names(tbl_display)) {
      names(tbl_display)[names(tbl_display) == lvl_col] <- facet_i
    }
    if (meas_col %in% names(tbl_display)) {
      names(tbl_display)[names(tbl_display) == meas_col] <- paste0(facet_i, " measr")
    }
  }

  bias_cols <- names(tbl_display)
  bias_formats <- list(
    Sq = "{}",
    `Observd Score` = "{:.2f}",
    `Expctd Score` = "{:.2f}",
    `Observd Count` = "{:.0f}",
    `Obs-Exp Average` = "{:.2f}",
    `Bias Size` = "{:.2f}",
    `Model S.E.` = "{:.2f}",
    t = "{:.2f}",
    `d.f.` = "{:.0f}",
    `Prob.` = "{:.4f}",
    `Infit MnSq` = "{:.2f}",
    `Outfit MnSq` = "{:.2f}"
  )
  for (facet_i in facets) {
    bias_formats[[paste0(facet_i, " N")]] <- "{:.0f}"
    bias_formats[[paste0(facet_i, " measr")]] <- "{:.2f}"
  }

  if (branch == "facets") {
    bias_fixed <- build_bias_fixed_text(
      table_df = tbl_display,
      summary_df = bias_results$summary,
      chi_df = bias_results$chi_sq,
      facet_a = facets[1],
      facet_b = if (length(facets) >= 2) facets[2] else "",
      interaction_label = interaction_label,
      columns = bias_cols,
      formats = bias_formats
    )
  } else {
    bias_fixed <- build_sectioned_fixed_report(
      title = paste0("Bias interaction summary: ", interaction_label),
      sections = list(
        list(
          title = "Top interaction rows",
          data = tbl_display,
          columns = bias_cols,
          formats = bias_formats,
          max_rows = 40L
        ),
        list(
          title = "Summary",
          data = as.data.frame(bias_results$summary %||% data.frame(), stringsAsFactors = FALSE)
        ),
        list(
          title = "Chi-square",
          data = as.data.frame(bias_results$chi_sq %||% data.frame(), stringsAsFactors = FALSE)
        )
      ),
      max_col_width = 18,
      min_col_width = 6
    )
  }

  if (length(facets) != 2) {
    pairwise_tbl <- tibble::tibble()
    pairwise_fixed <- paste0(
      "Legacy-compatible pairwise contrasts are available only for 2-way interactions.\n",
      "Current interaction: ", interaction_label, " (order ", length(facets), ")."
    )
  } else {
    facet_a <- facets[1]
    facet_b <- facets[2]
    if (is.null(target_facet)) target_facet <- facet_a
    context_facet <- ifelse(target_facet == facet_a, facet_b, facet_a)
    pairwise_tbl <- calc_bias_pairwise(bias_results$table, target_facet, context_facet)

    pairwise_fixed <- if (nrow(pairwise_tbl) == 0) {
      "No pairwise data"
    } else {
      pair_display <- pairwise_tbl |>
        dplyr::select(
          `Target N`,
          Target,
          `Target Measure`,
          `Target S.E.`,
          `Context1 N`,
          Context1,
          `Local Measure1`,
          SE1,
          `Obs-Exp Avg1`,
          Count1,
          `Context2 N`,
          Context2,
          `Local Measure2`,
          SE2,
          `Obs-Exp Avg2`,
          Count2,
          Contrast,
          SE,
          t,
          `d.f.`,
          `Prob.`
        ) |>
        dplyr::rename(
          `Target Measr` = `Target Measure`,
          `Context1 Measr` = `Local Measure1`,
          `Context1 S.E.` = SE1,
          `Context2 Measr` = `Local Measure2`,
          `Context2 S.E.` = SE2
        )

      pair_cols <- c(
        "Target N", "Target", "Target Measr", "Target S.E.",
        "Context1 N", "Context1", "Context1 Measr", "Context1 S.E.",
        "Obs-Exp Avg1", "Count1", "Context2 N", "Context2", "Context2 Measr", "Context2 S.E.",
        "Obs-Exp Avg2", "Count2", "Contrast", "SE", "t", "d.f.", "Prob."
      )

      pair_formats <- list(
        `Target N` = "{:.0f}",
        `Target Measr` = "{:.2f}",
        `Target S.E.` = "{:.2f}",
        `Context1 N` = "{:.0f}",
        `Context1 Measr` = "{:.2f}",
        `Context1 S.E.` = "{:.2f}",
        `Obs-Exp Avg1` = "{:.2f}",
        Count1 = "{:.0f}",
        `Context2 N` = "{:.0f}",
        `Context2 Measr` = "{:.2f}",
        `Context2 S.E.` = "{:.2f}",
        `Obs-Exp Avg2` = "{:.2f}",
        Count2 = "{:.0f}",
        Contrast = "{:.2f}",
        SE = "{:.2f}",
        t = "{:.2f}",
        `d.f.` = "{:.0f}",
        `Prob.` = "{:.4f}"
      )

      if (branch == "facets") {
        build_pairwise_fixed_text(
          pair_df = pair_display,
          target_facet = target_facet,
          context_facet = context_facet,
          columns = pair_cols,
          formats = pair_formats
        )
      } else {
        build_sectioned_fixed_report(
          title = paste0("Pairwise contrast summary: ", target_facet, " within ", context_facet),
          sections = list(
            list(
              title = "Pairwise rows",
              data = pair_display,
              columns = pair_cols,
              formats = pair_formats,
              max_rows = 40L
            )
          ),
          max_col_width = 18,
          min_col_width = 6
        )
      }
    }
  }

  out <- list(
    bias_fixed = bias_fixed,
    pairwise_fixed = pairwise_fixed,
    pairwise_table = pairwise_tbl,
    branch = branch,
    style = style,
    interaction_label = interaction_label,
    target_facet = as.character(target_facet %||% NA_character_)
  )
  out <- as_mfrm_bundle(out, "mfrm_fixed_reports")
  class(out) <- unique(c(paste0("mfrm_fixed_reports_", branch), class(out)))
  out
}

normalize_bias_plot_input <- function(x,
                                      diagnostics = NULL,
                                      facet_a = NULL,
                                      facet_b = NULL,
                                      interaction_facets = NULL,
                                      max_abs = 10,
                                      omit_extreme = TRUE,
                                      max_iter = 4,
                                      tol = 1e-3) {
  if (is.list(x) && !is.null(x$table) && !is.null(x$summary) && !is.null(x$chi_sq)) {
    return(x)
  }
  if (inherits(x, "mfrm_fit")) {
    if (is.null(interaction_facets)) {
      if (is.null(facet_a) || is.null(facet_b)) {
        stop("When `x` is mfrm_fit, provide `interaction_facets` or both `facet_a` and `facet_b`.")
      }
      interaction_facets <- c(as.character(facet_a[1]), as.character(facet_b[1]))
    }
    interaction_facets <- as.character(interaction_facets)
    interaction_facets <- interaction_facets[!is.na(interaction_facets) & nzchar(interaction_facets)]
    interaction_facets <- unique(interaction_facets)
    if (length(interaction_facets) < 2) {
      stop("`interaction_facets` must contain at least two facet names.")
    }
    if (is.null(facet_a)) facet_a <- interaction_facets[1]
    if (is.null(facet_b) && length(interaction_facets) >= 2) facet_b <- interaction_facets[2]
    if (is.null(facet_b)) {
      stop("`interaction_facets` must contain at least two facet names.")
    }
    if (is.null(diagnostics)) {
      diagnostics <- diagnose_mfrm(x, residual_pca = "none")
    }
    return(estimate_bias(
      fit = x,
      diagnostics = diagnostics,
      facet_a = facet_a,
      facet_b = facet_b,
      interaction_facets = interaction_facets,
      max_abs = max_abs,
      omit_extreme = omit_extreme,
      max_iter = max_iter,
      tol = tol
    ))
  }
  stop("`x` must be output from estimate_bias() or an mfrm_fit object.")
}

#' Build a legacy-compatible Table 13 bias-plot export bundle
#'
#' @param x Output from [estimate_bias()] or [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()] (used when `x` is fit).
#' @param facet_a First facet name (required when `x` is fit).
#' @param facet_b Second facet name (required when `x` is fit).
#' @param interaction_facets Character vector of two or more facets (required
#'   when `x` is fit and higher-order interaction output is needed).
#' @param max_abs Bound for absolute bias size when estimating from fit.
#' @param omit_extreme Omit extreme-only elements when estimating from fit.
#' @param max_iter Iteration cap for bias estimation when `x` is fit.
#' @param tol Convergence tolerance for bias estimation when `x` is fit.
#' @param top_n Maximum number of ranked rows to keep.
#' @param abs_t_warn Warning cutoff for absolute t statistics.
#' @param abs_bias_warn Warning cutoff for absolute bias size.
#' @param p_max Warning cutoff for p-values.
#' @param sort_by Ranking key: `"abs_t"`, `"abs_bias"`, or `"prob"`.
#'
#' @details
#' The legacy-compatible Table 13 layout is often inspected graphically
#' (bias size, observed-minus-expected average, significance). This helper
#' prepares a plotting-ready bundle with:
#' - ranked table for lollipop/strip displays
#' - scatter table (`Obs-Exp Average` vs `Bias Size`)
#' - summary and threshold metadata
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [bias_interaction_report()].
#'
#' @return A named list with:
#' - `ranked_table`: top-ranked bias rows with flags
#' - `scatter_data`: bias scatter data with flags
#' - `facet_profile`: per-facet level profile (`MeanAbsBias`, `FlagRate`)
#' - `plot_long`: unified long-format plotting table spanning scatter,
#'   ranked, profile, and heatmap views
#' - `plot_annotations`, `flag_summary`, `plot_settings`: reusable threshold,
#'   screening-trigger, and rendering metadata
#' - `summary`: one-row overview
#' - `thresholds`: applied cutoffs
#' - `facet_a`, `facet_b`: first two analyzed facet names
#' - `interaction_facets`, `interaction_order`, `interaction_mode`: full
#'   interaction metadata
#'
#' @seealso [estimate_bias()], [plot_bias_interaction()], [build_fixed_reports()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' t13 <- bias_interaction_report(bias, top_n = 10)
#' @keywords internal
#' @noRd
table13_bias_plot_export <- function(x,
                                     diagnostics = NULL,
                                     facet_a = NULL,
                                     facet_b = NULL,
                                     interaction_facets = NULL,
                                     max_abs = 10,
                                     omit_extreme = TRUE,
                                     max_iter = 4,
                                     tol = 1e-3,
                                     top_n = 50,
                                     abs_t_warn = 2,
                                     abs_bias_warn = 0.5,
                                     p_max = 0.05,
                                     sort_by = c("abs_t", "abs_bias", "prob")) {
  signal_legacy_name_deprecation(
    old_name = "table13_bias_plot_export",
    new_name = "bias_interaction_report",
    suppress_if_called_from = c("bias_interaction_report", "plot_table13_bias", "plot_bias_interaction")
  )
  sort_by <- match.arg(tolower(sort_by), c("abs_t", "abs_bias", "prob"))
  top_n <- max(1L, as.integer(top_n))
  abs_t_warn <- abs(as.numeric(abs_t_warn))
  abs_bias_warn <- abs(as.numeric(abs_bias_warn))
  p_max <- max(0, min(1, as.numeric(p_max)))

  bias_results <- normalize_bias_plot_input(
    x = x,
    diagnostics = diagnostics,
    facet_a = facet_a,
    facet_b = facet_b,
    interaction_facets = interaction_facets,
    max_abs = max_abs,
    omit_extreme = omit_extreme,
    max_iter = max_iter,
    tol = tol
  )
  spec <- extract_bias_facet_spec(bias_results)
  if (is.null(spec) || length(spec$facets) < 2) {
    stop("`bias_results$table` does not include recognizable interaction facet columns.")
  }
  interaction_facets <- spec$facets
  interaction_order <- spec$interaction_order
  interaction_mode <- spec$interaction_mode
  facet_a <- interaction_facets[1]
  facet_b <- interaction_facets[2]

  tbl <- as.data.frame(bias_results$table, stringsAsFactors = FALSE)
  req <- c(spec$level_cols, "Obs-Exp Average", "Bias Size", "t", "Prob.")
  if (!all(req %in% names(tbl))) {
    stop("`bias_results$table` does not include required Table 13 columns.")
  }

  level_df <- tbl[, spec$level_cols, drop = FALSE]
  level_df[] <- lapply(level_df, as.character)
  names(level_df) <- paste0("Level", seq_along(spec$level_cols))

  tbl2 <- dplyr::bind_cols(
    data.frame(
      InteractionFacets = paste(interaction_facets, collapse = " x "),
      InteractionOrder = interaction_order,
      InteractionMode = interaction_mode,
      FacetA = facet_a,
      FacetB = facet_b,
      stringsAsFactors = FALSE
    )[rep(1, nrow(tbl)), , drop = FALSE],
    level_df,
    tbl |>
      dplyr::transmute(
        ObsExpAverage = suppressWarnings(as.numeric(.data$`Obs-Exp Average`)),
        BiasSize = suppressWarnings(as.numeric(.data$`Bias Size`)),
        SE = if ("S.E." %in% names(tbl)) suppressWarnings(as.numeric(.data$`S.E.`)) else NA_real_,
        t = suppressWarnings(as.numeric(.data$t)),
        Prob = suppressWarnings(as.numeric(.data$`Prob.`)),
        ObservedCount = if ("Observd Count" %in% names(tbl)) suppressWarnings(as.numeric(.data$`Observd Count`)) else NA_real_,
        LRChiSq = if ("LR ChiSq" %in% names(tbl)) suppressWarnings(as.numeric(.data$`LR ChiSq`)) else NA_real_,
        LRDF = if ("LR d.f." %in% names(tbl)) suppressWarnings(as.numeric(.data$`LR d.f.`)) else NA_real_,
        LRProb = if ("LR Prob." %in% names(tbl)) suppressWarnings(as.numeric(.data$`LR Prob.`)) else NA_real_,
        ProfileCILower = if ("Profile CI Lower" %in% names(tbl)) suppressWarnings(as.numeric(.data$`Profile CI Lower`)) else NA_real_,
        ProfileCIUpper = if ("Profile CI Upper" %in% names(tbl)) suppressWarnings(as.numeric(.data$`Profile CI Upper`)) else NA_real_,
        ProfileCILevel = if ("Profile CI Level" %in% names(tbl)) suppressWarnings(as.numeric(.data$`Profile CI Level`)) else NA_real_,
        ProfileCIStatus = if ("Profile CI Status" %in% names(tbl)) as.character(.data$`Profile CI Status`) else NA_character_,
        LikelihoodBasis = if ("Likelihood Basis" %in% names(tbl)) as.character(.data$`Likelihood Basis`) else NA_character_
      )
  ) |>
    dplyr::mutate(
      Pair = do.call(paste, c(level_df, sep = " | ")),
      AbsT = abs(.data$t),
      AbsBias = abs(.data$BiasSize),
      TFlag = is.finite(.data$AbsT) & .data$AbsT >= abs_t_warn,
      BiasFlag = is.finite(.data$AbsBias) & .data$AbsBias >= abs_bias_warn,
      PFlag = is.finite(.data$Prob) & .data$Prob <= p_max,
      Flag = .data$TFlag | .data$BiasFlag | .data$PFlag
    )

  for (i in seq_along(interaction_facets)) {
    tbl2[[paste0("Facet", i)]] <- interaction_facets[i]
    tbl2[[paste0("Facet", i, "_Level")]] <- level_df[[i]]
  }
  # Keep legacy aliases for 2-way compatibility.
  tbl2$FacetA_Level <- tbl2$Facet1_Level
  tbl2$FacetB_Level <- tbl2$Facet2_Level

  ord <- switch(
    sort_by,
    abs_t = order(tbl2$AbsT, decreasing = TRUE, na.last = NA),
    abs_bias = order(tbl2$AbsBias, decreasing = TRUE, na.last = NA),
    prob = order(tbl2$Prob, decreasing = FALSE, na.last = NA)
  )
  ranked <- if (length(ord) == 0) tbl2[0, , drop = FALSE] else tbl2[ord, , drop = FALSE]
  if (nrow(ranked) > top_n) ranked <- ranked[seq_len(top_n), , drop = FALSE]

  scatter_cols <- c(
    "InteractionFacets", "InteractionOrder", "InteractionMode",
    "FacetA", "FacetA_Level", "FacetB", "FacetB_Level",
    paste0("Facet", seq_along(interaction_facets)),
    paste0("Facet", seq_along(interaction_facets), "_Level"),
    "Pair", "ObsExpAverage", "BiasSize", "SE", "t", "Prob",
    "LRChiSq", "LRDF", "LRProb", "ProfileCILower", "ProfileCIUpper",
    "ProfileCILevel", "ProfileCIStatus", "LikelihoodBasis",
    "ObservedCount", "Flag", "TFlag", "BiasFlag", "PFlag"
  )
  scatter_cols <- unique(scatter_cols)
  scatter <- tbl2 |>
    dplyr::select(dplyr::all_of(scatter_cols))

  profile_rows <- lapply(seq_along(interaction_facets), function(i) {
    facet_i <- interaction_facets[i]
    level_col <- paste0("Facet", i, "_Level")
    tbl2 |>
      dplyr::group_by(Level = .data[[level_col]]) |>
      dplyr::summarize(
        Cells = dplyr::n(),
        MeanAbsBias = mean(.data$AbsBias, na.rm = TRUE),
        MeanAbsT = mean(.data$AbsT, na.rm = TRUE),
        Flagged = sum(.data$Flag, na.rm = TRUE),
        FlagRate = ifelse(dplyr::n() > 0, 100 * sum(.data$Flag, na.rm = TRUE) / dplyr::n(), NA_real_),
        .groups = "drop"
      ) |>
      dplyr::mutate(Facet = facet_i, .before = 1)
  })
  facet_profile <- dplyr::bind_rows(profile_rows)

  n_flagged <- sum(tbl2$Flag, na.rm = TRUE)
  n_computed <- sum(is.finite(tbl2$AbsT) | is.finite(tbl2$AbsBias), na.rm = TRUE)
  flag_status <- if (nrow(tbl2) == 0L) {
    "No cells constructed for this facet pair."
  } else if (n_computed == 0L) {
    paste0(
      "No flag statistics were computable (all t/bias columns returned NA); ",
      "inspect `bias_results$optimization_failures` before interpreting."
    )
  } else if (n_flagged == 0L) {
    sprintf(
      "No cells crossed screening thresholds (|t| >= %s, |Bias| >= %s, p <= %s); top-ranked cells still appear in ranked_table for review.",
      format(abs_t_warn), format(abs_bias_warn), format(p_max)
    )
  } else {
    sprintf("%d of %d cell(s) flagged.", n_flagged, nrow(tbl2))
  }

  summary_tbl <- data.frame(
    InteractionFacets = paste(interaction_facets, collapse = " x "),
    InteractionOrder = interaction_order,
    InteractionMode = interaction_mode,
    FacetA = facet_a,
    FacetB = facet_b,
    Cells = nrow(tbl2),
    Flagged = n_flagged,
    FlaggedPercent = ifelse(nrow(tbl2) > 0, 100 * n_flagged / nrow(tbl2), NA_real_),
    MeanAbsT = mean(tbl2$AbsT, na.rm = TRUE),
    MeanAbsBias = mean(tbl2$AbsBias, na.rm = TRUE),
    FlagStatus = flag_status,
    stringsAsFactors = FALSE
  )

  list(
    ranked_table = as.data.frame(ranked, stringsAsFactors = FALSE),
    scatter_data = as.data.frame(scatter, stringsAsFactors = FALSE),
    facet_profile = as.data.frame(facet_profile, stringsAsFactors = FALSE),
    summary = summary_tbl,
    thresholds = list(
      abs_t_warn = abs_t_warn,
      abs_bias_warn = abs_bias_warn,
      p_max = p_max,
      sort_by = sort_by,
      top_n = top_n
    ),
    facet_a = facet_a,
    facet_b = facet_b,
    interaction_facets = interaction_facets,
    interaction_order = interaction_order,
    interaction_mode = interaction_mode,
    orientation_review = as.data.frame(bias_results$orientation_review %||% data.frame(), stringsAsFactors = FALSE),
    mixed_sign = isTRUE(bias_results$mixed_sign),
    direction_note = as.character(bias_results$direction_note %||% ""),
    recommended_action = as.character(bias_results$recommended_action %||% "")
  )
}

#' Plot interaction-bias diagnostics using base R
#'
#' @param x Output from [bias_interaction_report()], [estimate_bias()], or [fit_mfrm()].
#' @param plot Plot type: `"scatter"`, `"ranked"`, `"heatmap"`,
#'   `"abs_t_hist"`, or `"facet_profile"`.
#' @param diagnostics Optional output from [diagnose_mfrm()] when `x` is fit.
#' @param facet_a First facet name (required when `x` is fit).
#' @param facet_b Second facet name (required when `x` is fit).
#' @param interaction_facets Character vector of two or more facets (required
#'   when `x` is fit and higher-order interaction output is needed).
#' @param top_n Maximum number of ranked rows to show.
#' @param abs_t_warn Warning cutoff for absolute t.
#' @param abs_bias_warn Warning cutoff for absolute bias size.
#' @param p_max Warning cutoff for p-values.
#' @param sort_by Ranking key for `"ranked"` plot.
#' @param main Optional plot title override.
#' @param palette Optional named color overrides (`normal`, `flag`, `hist`,
#'   `profile`).
#' @param label_angle Label angle hint for ranked/profile labels.
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`, or `"monochrome"`).
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [plot_bias_interaction()].
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [bias_interaction_report()], [estimate_bias()], [plot_displacement()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 30)
#' p13 <- plot_bias_interaction(
#'   fit,
#'   diagnostics = diagnose_mfrm(fit, residual_pca = "none"),
#'   facet_a = "Rater",
#'   facet_b = "Criterion",
#'   draw = FALSE
#' )
#' @keywords internal
#' @noRd
build_bias_heatmap_data <- function(scatter) {
  empty_matrix <- matrix(numeric(0), nrow = 0L, ncol = 0L)
  reference_notes <- data.frame(
    Reference = c(
      "Linacre (1989); FACETS Table 13 / Table 14 interaction logic",
      "FACETS Tutorial 3: bias / interaction report",
      "mfrmr estimate_bias() screening label"
    ),
    Role = c(
      "Local behavior is contrasted with general facet behavior before substantive interpretation.",
      "Table 13 highlights conspicuous interaction cells; plots are follow-up aids, not replacements for the table.",
      "Heatmap cells inherit the conditional screening interpretation of estimate_bias()."
    ),
    stringsAsFactors = FALSE
  )
  interpretation_guide <- data.frame(
    Topic = c("Cell", "Sign", "Flag outline", "Higher-order interactions", "Use of results"),
    Guidance = c(
      "Each cell represents a Facet1 level by Facet2 level interaction row from the bias table.",
      "Positive bias means higher observed ratings than the additive model expected; negative bias means lower observed ratings.",
      "Flagged cells crossed at least one current screening rule for |t|, |bias|, or tail probability.",
      "When the source bias table contains more than two facets, the heatmap keeps the largest absolute bias row for each first-two-facet cell and reports the number of collapsed rows.",
      "Use the heatmap for pattern recognition and follow-up; it does not upgrade screening statistics to confirmatory inference."
    ),
    stringsAsFactors = FALSE
  )

  scatter <- as.data.frame(scatter %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(scatter) == 0L ||
      !all(c("Facet1_Level", "Facet2_Level", "BiasSize") %in% names(scatter))) {
    return(list(
      cells = data.frame(),
      matrix = empty_matrix,
      flag_matrix = matrix(FALSE, nrow = 0L, ncol = 0L),
      count_matrix = empty_matrix,
      x_levels = character(0),
      y_levels = character(0),
      collapsed = FALSE,
      interpretation_guide = interpretation_guide,
      reference_notes = reference_notes
    ))
  }

  sc <- scatter[is.finite(suppressWarnings(as.numeric(scatter$BiasSize))), , drop = FALSE]
  if (nrow(sc) == 0L) {
    return(list(
      cells = data.frame(),
      matrix = empty_matrix,
      flag_matrix = matrix(FALSE, nrow = 0L, ncol = 0L),
      count_matrix = empty_matrix,
      x_levels = character(0),
      y_levels = character(0),
      collapsed = FALSE,
      interpretation_guide = interpretation_guide,
      reference_notes = reference_notes
    ))
  }

  sc$HeatmapY <- as.character(sc$Facet1_Level)
  sc$HeatmapX <- as.character(sc$Facet2_Level)
  sc$BiasSize <- suppressWarnings(as.numeric(sc$BiasSize))
  sc$AbsBias <- abs(sc$BiasSize)
  sc$Flag <- as.logical(sc$Flag %||% FALSE)
  sc$TFlag <- as.logical(sc$TFlag %||% FALSE)
  sc$BiasFlag <- as.logical(sc$BiasFlag %||% FALSE)
  sc$PFlag <- as.logical(sc$PFlag %||% FALSE)
  sc$ObservedCount <- suppressWarnings(as.numeric(sc$ObservedCount %||% NA_real_))
  sc$SE <- suppressWarnings(as.numeric(sc$SE %||% NA_real_))
  sc$t <- suppressWarnings(as.numeric(sc$t %||% NA_real_))
  sc$Prob <- suppressWarnings(as.numeric(sc$Prob %||% NA_real_))

  cells <- sc |>
    dplyr::group_by(.data$HeatmapY, .data$HeatmapX) |>
    dplyr::arrange(dplyr::desc(.data$AbsBias), .by_group = TRUE) |>
    dplyr::summarise(
      Facet1_Level = dplyr::first(.data$Facet1_Level),
      Facet2_Level = dplyr::first(.data$Facet2_Level),
      Pair = dplyr::first(.data$Pair),
      BiasSize = dplyr::first(.data$BiasSize),
      AbsBias = dplyr::first(.data$AbsBias),
      SE = dplyr::first(.data$SE),
      t = dplyr::first(.data$t),
      Prob = dplyr::first(.data$Prob),
      ObservedCount = if (any(is.finite(.data$ObservedCount))) {
        sum(.data$ObservedCount[is.finite(.data$ObservedCount)], na.rm = TRUE)
      } else {
        NA_real_
      },
      Flag = any(.data$Flag, na.rm = TRUE),
      TFlag = any(.data$TFlag, na.rm = TRUE),
      BiasFlag = any(.data$BiasFlag, na.rm = TRUE),
      PFlag = any(.data$PFlag, na.rm = TRUE),
      CollapsedRows = dplyr::n(),
      RepresentativeRule = "largest_abs_bias",
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$HeatmapY, .data$HeatmapX) |>
    as.data.frame(stringsAsFactors = FALSE)

  y_levels <- sort(unique(as.character(cells$HeatmapY)))
  x_levels <- sort(unique(as.character(cells$HeatmapX)))
  mat_bias <- matrix(NA_real_, nrow = length(y_levels), ncol = length(x_levels),
                     dimnames = list(y_levels, x_levels))
  mat_flag <- matrix(FALSE, nrow = length(y_levels), ncol = length(x_levels),
                     dimnames = list(y_levels, x_levels))
  mat_count <- matrix(NA_real_, nrow = length(y_levels), ncol = length(x_levels),
                      dimnames = list(y_levels, x_levels))

  for (k in seq_len(nrow(cells))) {
    i <- match(as.character(cells$HeatmapY[k]), y_levels)
    j <- match(as.character(cells$HeatmapX[k]), x_levels)
    if (is.finite(i) && is.finite(j)) {
      mat_bias[i, j] <- as.numeric(cells$BiasSize[k])
      mat_flag[i, j] <- isTRUE(as.logical(cells$Flag[k]))
      mat_count[i, j] <- as.numeric(cells$ObservedCount[k])
    }
  }

  list(
    cells = cells,
    matrix = mat_bias,
    flag_matrix = mat_flag,
    count_matrix = mat_count,
    x_levels = x_levels,
    y_levels = y_levels,
    collapsed = any(cells$CollapsedRows > 1L),
    interpretation_guide = interpretation_guide,
    reference_notes = reference_notes
  )
}

build_bias_plot_long <- function(ranked,
                                 scatter,
                                 profile,
                                 heatmap_cells,
                                 heatmap_levels) {
  ranked <- as.data.frame(ranked %||% data.frame(), stringsAsFactors = FALSE)
  scatter <- as.data.frame(scatter %||% data.frame(), stringsAsFactors = FALSE)
  profile <- as.data.frame(profile %||% data.frame(), stringsAsFactors = FALSE)
  heatmap_cells <- as.data.frame(heatmap_cells %||% data.frame(), stringsAsFactors = FALSE)
  x_levels <- as.character(heatmap_levels$x %||% character())
  y_levels <- as.character(heatmap_levels$y %||% character())

  common_cols <- c(
    "Layer", "PlotType", "Series", "Pair", "Facet", "Level",
    "FacetA_Level", "FacetB_Level", "X", "Y", "XLabel", "YLabel",
    "Value", "ValueName", "Rank", "Flag", "TFlag", "BiasFlag", "PFlag",
    "ObservedCount", "DisplayedByDefault"
  )
  empty <- data.frame(
    Layer = character(),
    PlotType = character(),
    Series = character(),
    Pair = character(),
    Facet = character(),
    Level = character(),
    FacetA_Level = character(),
    FacetB_Level = character(),
    X = numeric(),
    Y = numeric(),
    XLabel = character(),
    YLabel = character(),
    Value = numeric(),
    ValueName = character(),
    Rank = integer(),
    Flag = logical(),
    TFlag = logical(),
    BiasFlag = logical(),
    PFlag = logical(),
    ObservedCount = numeric(),
    DisplayedByDefault = logical(),
    stringsAsFactors = FALSE
  )
  normalize <- function(df) {
    for (nm in common_cols) {
      if (!nm %in% names(df)) df[[nm]] <- empty[[nm]]
    }
    df[, common_cols, drop = FALSE]
  }

  scatter_long <- if (nrow(scatter) > 0L) {
    data.frame(
      Layer = "scatter",
      PlotType = "scatter",
      Series = as.character(scatter$Pair %||% seq_len(nrow(scatter))),
      Pair = as.character(scatter$Pair %||% NA_character_),
      Facet = NA_character_,
      Level = NA_character_,
      FacetA_Level = as.character(scatter$FacetA_Level %||% NA_character_),
      FacetB_Level = as.character(scatter$FacetB_Level %||% NA_character_),
      X = suppressWarnings(as.numeric(scatter$ObsExpAverage %||% NA_real_)),
      Y = suppressWarnings(as.numeric(scatter$BiasSize %||% NA_real_)),
      XLabel = "Obs-Exp Average",
      YLabel = "Bias Size",
      Value = suppressWarnings(as.numeric(scatter$BiasSize %||% NA_real_)),
      ValueName = "BiasSize",
      Rank = NA_integer_,
      Flag = as.logical(scatter$Flag %||% FALSE),
      TFlag = as.logical(scatter$TFlag %||% FALSE),
      BiasFlag = as.logical(scatter$BiasFlag %||% FALSE),
      PFlag = as.logical(scatter$PFlag %||% FALSE),
      ObservedCount = suppressWarnings(as.numeric(scatter$ObservedCount %||% NA_real_)),
      DisplayedByDefault = TRUE,
      stringsAsFactors = FALSE
    )
  } else empty

  ranked_long <- if (nrow(ranked) > 0L) {
    rank_idx <- seq_len(nrow(ranked))
    data.frame(
      Layer = "ranked",
      PlotType = "ranked",
      Series = as.character(ranked$Pair %||% rank_idx),
      Pair = as.character(ranked$Pair %||% NA_character_),
      Facet = NA_character_,
      Level = NA_character_,
      FacetA_Level = as.character(ranked$FacetA_Level %||% NA_character_),
      FacetB_Level = as.character(ranked$FacetB_Level %||% NA_character_),
      X = rank_idx,
      Y = suppressWarnings(as.numeric(ranked$BiasSize %||% NA_real_)),
      XLabel = as.character(rank_idx),
      YLabel = as.character(ranked$Pair %||% NA_character_),
      Value = suppressWarnings(as.numeric(ranked$BiasSize %||% NA_real_)),
      ValueName = "BiasSize",
      Rank = rank_idx,
      Flag = as.logical(ranked$Flag %||% FALSE),
      TFlag = as.logical(ranked$TFlag %||% FALSE),
      BiasFlag = as.logical(ranked$BiasFlag %||% FALSE),
      PFlag = as.logical(ranked$PFlag %||% FALSE),
      ObservedCount = suppressWarnings(as.numeric(ranked$ObservedCount %||% NA_real_)),
      DisplayedByDefault = TRUE,
      stringsAsFactors = FALSE
    )
  } else empty

  profile_long <- if (nrow(profile) > 0L &&
      all(c("Facet", "Level", "MeanAbsBias", "FlagRate") %in% names(profile))) {
    data.frame(
      Layer = "facet_profile",
      PlotType = "facet_profile",
      Series = paste0(as.character(profile$Facet), ": ", as.character(profile$Level)),
      Pair = NA_character_,
      Facet = as.character(profile$Facet),
      Level = as.character(profile$Level),
      FacetA_Level = NA_character_,
      FacetB_Level = NA_character_,
      X = suppressWarnings(as.numeric(profile$MeanAbsBias)),
      Y = suppressWarnings(as.numeric(profile$FlagRate)),
      XLabel = "Mean |Bias Size|",
      YLabel = "Flag Rate",
      Value = suppressWarnings(as.numeric(profile$MeanAbsBias)),
      ValueName = "MeanAbsBias",
      Rank = NA_integer_,
      Flag = suppressWarnings(as.numeric(profile$FlagRate)) > 0,
      TFlag = NA,
      BiasFlag = NA,
      PFlag = NA,
      ObservedCount = suppressWarnings(as.numeric(profile$Cells %||% NA_real_)),
      DisplayedByDefault = TRUE,
      stringsAsFactors = FALSE
    )
  } else empty

  heatmap_long <- if (nrow(heatmap_cells) > 0L &&
      all(c("Facet1_Level", "Facet2_Level", "BiasSize") %in% names(heatmap_cells))) {
    data.frame(
      Layer = "heatmap_cell",
      PlotType = "heatmap",
      Series = as.character(heatmap_cells$Pair %||% seq_len(nrow(heatmap_cells))),
      Pair = as.character(heatmap_cells$Pair %||% NA_character_),
      Facet = NA_character_,
      Level = NA_character_,
      FacetA_Level = as.character(heatmap_cells$Facet1_Level),
      FacetB_Level = as.character(heatmap_cells$Facet2_Level),
      X = match(as.character(heatmap_cells$Facet2_Level), x_levels),
      Y = match(as.character(heatmap_cells$Facet1_Level), y_levels),
      XLabel = as.character(heatmap_cells$Facet2_Level),
      YLabel = as.character(heatmap_cells$Facet1_Level),
      Value = suppressWarnings(as.numeric(heatmap_cells$BiasSize)),
      ValueName = "BiasSize",
      Rank = NA_integer_,
      Flag = as.logical(heatmap_cells$Flag %||% FALSE),
      TFlag = as.logical(heatmap_cells$TFlag %||% FALSE),
      BiasFlag = as.logical(heatmap_cells$BiasFlag %||% FALSE),
      PFlag = as.logical(heatmap_cells$PFlag %||% FALSE),
      ObservedCount = suppressWarnings(as.numeric(heatmap_cells$ObservedCount %||% NA_real_)),
      DisplayedByDefault = TRUE,
      stringsAsFactors = FALSE
    )
  } else empty

  out <- rbind(
    normalize(scatter_long),
    normalize(ranked_long),
    normalize(profile_long),
    normalize(heatmap_long)
  )
  rownames(out) <- NULL
  out
}

build_bias_plot_annotations <- function(thresholds) {
  data.frame(
    AnnotationType = c(
      "bias_threshold_negative",
      "bias_center",
      "bias_threshold_positive",
      "abs_t_threshold",
      "p_threshold"
    ),
    Axis = c("y", "y", "y", "x", "p"),
    Value = c(
      -abs(suppressWarnings(as.numeric(thresholds$abs_bias_warn %||% NA_real_))),
      0,
      abs(suppressWarnings(as.numeric(thresholds$abs_bias_warn %||% NA_real_))),
      abs(suppressWarnings(as.numeric(thresholds$abs_t_warn %||% NA_real_))),
      suppressWarnings(as.numeric(thresholds$p_max %||% NA_real_))
    ),
    Label = c(
      "Negative bias review threshold",
      "Centered bias reference",
      "Positive bias review threshold",
      "Absolute t screening threshold",
      "Tail-area screening threshold"
    ),
    stringsAsFactors = FALSE
  )
}

build_bias_flag_summary <- function(scatter) {
  scatter <- as.data.frame(scatter %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(scatter) == 0L) {
    return(data.frame(
      Trigger = c("any_flag", "t", "bias_size", "probability"),
      Rows = integer(4),
      Percent = rep(NA_real_, 4),
      stringsAsFactors = FALSE
    ))
  }
  n <- nrow(scatter)
  counts <- c(
    any_flag = sum(as.logical(scatter$Flag %||% FALSE), na.rm = TRUE),
    t = sum(as.logical(scatter$TFlag %||% FALSE), na.rm = TRUE),
    bias_size = sum(as.logical(scatter$BiasFlag %||% FALSE), na.rm = TRUE),
    probability = sum(as.logical(scatter$PFlag %||% FALSE), na.rm = TRUE)
  )
  data.frame(
    Trigger = names(counts),
    Rows = as.integer(counts),
    Percent = if (n > 0L) 100 * as.numeric(counts) / n else NA_real_,
    stringsAsFactors = FALSE
  )
}

plot_table13_bias <- function(x,
                              plot = c("scatter", "ranked", "heatmap", "abs_t_hist", "facet_profile"),
                              diagnostics = NULL,
                              facet_a = NULL,
                              facet_b = NULL,
                              interaction_facets = NULL,
                              top_n = 40,
                              abs_t_warn = 2,
                              abs_bias_warn = 0.5,
                              p_max = 0.05,
                              sort_by = c("abs_t", "abs_bias", "prob"),
                              show_ci = FALSE,
                              ci_level = 0.95,
                              main = NULL,
                              palette = NULL,
                              label_angle = 45,
                              preset = c("standard", "publication", "compact", "monochrome"),
                              draw = TRUE) {
  signal_legacy_name_deprecation(
    old_name = "plot_table13_bias",
    new_name = "plot_bias_interaction",
    suppress_if_called_from = "plot_bias_interaction"
  )
  plot <- match.arg(tolower(plot), c("scatter", "ranked", "heatmap", "abs_t_hist", "facet_profile"))
  sort_by <- match.arg(tolower(sort_by), c("abs_t", "abs_bias", "prob"))
  top_n <- max(1L, as.integer(top_n))
  label_angle <- suppressWarnings(as.numeric(label_angle[1]))
  if (!is.finite(label_angle)) label_angle <- 45
  las_rank <- if (label_angle >= 45) 2 else 1
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number in (0, 1).", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  pal <- resolve_palette(
    palette = palette,
    defaults = c(
      normal = style$accent_primary,
      flag = style$fail,
      hist = style$fill_soft,
      profile = style$accent_tertiary
    )
  )

  bundle <- if (is.list(x) && all(c("ranked_table", "scatter_data", "summary", "thresholds") %in% names(x))) {
    x
  } else {
    with_legacy_name_warning_suppressed(
      table13_bias_plot_export(
        x = x,
        diagnostics = diagnostics,
        facet_a = facet_a,
        facet_b = facet_b,
        interaction_facets = interaction_facets,
        top_n = top_n,
        abs_t_warn = abs_t_warn,
        abs_bias_warn = abs_bias_warn,
        p_max = p_max,
        sort_by = sort_by
      )
    )
  }

  ranked <- as.data.frame(bundle$ranked_table, stringsAsFactors = FALSE)
  scatter <- as.data.frame(bundle$scatter_data, stringsAsFactors = FALSE)
  profile <- as.data.frame(bundle$facet_profile %||% data.frame(), stringsAsFactors = FALSE)
  thr <- bundle$thresholds

  # Add CI bounds for Bias Size when requested.  Bounded GPCM rows carry
  # conditional profile-likelihood limits; otherwise we use the
  # conditional plug-in SE column.
  ci_requested <- isTRUE(show_ci) && plot %in% c("scatter", "ranked")
  profile_ci_available <- ci_requested &&
    all(c("ProfileCILower", "ProfileCIUpper") %in% names(ranked)) &&
    any(is.finite(suppressWarnings(as.numeric(ranked$ProfileCILower))) &
          is.finite(suppressWarnings(as.numeric(ranked$ProfileCIUpper))))
  se_ci_available <- ci_requested &&
    !profile_ci_available &&
    "SE" %in% names(ranked) &&
    any(is.finite(suppressWarnings(as.numeric(ranked$SE))))
  ci_available <- profile_ci_available || se_ci_available
  if (ci_available) {
    if (profile_ci_available) {
      ranked$CI_Lower <- suppressWarnings(as.numeric(ranked$ProfileCILower))
      ranked$CI_Upper <- suppressWarnings(as.numeric(ranked$ProfileCIUpper))
      ranked$CI_Level <- suppressWarnings(as.numeric(ranked$ProfileCILevel %||% 0.95))
      ranked$CI_Method <- "conditional profile likelihood"
      if (all(c("ProfileCILower", "ProfileCIUpper") %in% names(scatter))) {
        scatter$CI_Lower <- suppressWarnings(as.numeric(scatter$ProfileCILower))
        scatter$CI_Upper <- suppressWarnings(as.numeric(scatter$ProfileCIUpper))
        scatter$CI_Level <- suppressWarnings(as.numeric(scatter$ProfileCILevel %||% 0.95))
        scatter$CI_Method <- "conditional profile likelihood"
      }
    } else {
      z_ci <- stats::qnorm(1 - (1 - ci_level) / 2)
      bias_r <- suppressWarnings(as.numeric(ranked$BiasSize))
      se_r   <- suppressWarnings(as.numeric(ranked$SE))
      ranked$CI_Lower <- bias_r - z_ci * se_r
      ranked$CI_Upper <- bias_r + z_ci * se_r
      ranked$CI_Level <- ci_level
      ranked$CI_Method <- "conditional plug-in SE"
      if ("SE" %in% names(scatter)) {
        bias_s <- suppressWarnings(as.numeric(scatter$BiasSize))
        se_s   <- suppressWarnings(as.numeric(scatter$SE))
        scatter$CI_Lower <- bias_s - z_ci * se_s
        scatter$CI_Upper <- bias_s + z_ci * se_s
        scatter$CI_Level <- ci_level
        scatter$CI_Method <- "conditional plug-in SE"
      }
    }
    if (profile_ci_available && !all(c("CI_Lower", "CI_Upper") %in% names(scatter)) &&
        "SE" %in% names(scatter)) {
      z_ci <- stats::qnorm(1 - (1 - ci_level) / 2)
      bias_s <- suppressWarnings(as.numeric(scatter$BiasSize))
      se_s   <- suppressWarnings(as.numeric(scatter$SE))
      scatter$CI_Lower <- bias_s - z_ci * se_s
      scatter$CI_Upper <- bias_s + z_ci * se_s
      scatter$CI_Level <- ci_level
      scatter$CI_Method <- "conditional plug-in SE"
    }
  }
  plot_title <- switch(
    plot,
    scatter = "Bias interaction scatter",
    ranked = "Ranked bias interaction size",
    heatmap = "Bias interaction heatmap",
    abs_t_hist = "Screening |t| distribution",
    facet_profile = "Facet-level bias profile"
  )
  if (!is.null(main)) plot_title <- as.character(main[1])
  plot_subtitle <- paste0(
    "Interaction facets: ",
    paste(bundle$summary$InteractionFacets[1] %||% paste(interaction_facets %||% c(facet_a, facet_b), collapse = " x ")),
    "; screening thresholds |t| >= ", format(abs_t_warn),
    ", |bias| >= ", format(abs_bias_warn)
  )
  plot_legend <- switch(
    plot,
    scatter = new_plot_legend(
      label = c("Within review band", "Screen-positive cell"),
      role = c("status", "status"),
      aesthetic = c("point", "point"),
      value = c(pal["normal"], pal["flag"])
    ),
    ranked = new_plot_legend(
      label = c("Within review band", "Screen-positive pair"),
      role = c("status", "status"),
      aesthetic = c("point", "point"),
      value = c(pal["normal"], pal["flag"])
    ),
    heatmap = new_plot_legend(
      label = c("Negative bias", "Centred", "Positive bias", "Flagged outline"),
      role = c("magnitude", "magnitude", "magnitude", "alert"),
      aesthetic = c("fill", "fill", "fill", "border"),
      value = c("#2166AC", "#FFFFFF", "#B2182B", pal["flag"])
    ),
    abs_t_hist = new_plot_legend(
      label = "Absolute screening t",
      role = "screening",
      aesthetic = "histogram",
      value = pal["hist"]
    ),
    facet_profile = new_plot_legend(
      label = c("No flagged cells", "Contains flagged cells"),
      role = c("status", "status"),
      aesthetic = c("point", "point"),
      value = c(pal["profile"], pal["flag"])
    )
  )
  plot_reference <- switch(
    plot,
    scatter = new_reference_lines(
      axis = c("h", "h", "h", "v"),
      value = c(-thr$abs_bias_warn, 0, thr$abs_bias_warn, 0),
      label = c("Bias review threshold", "Centered bias reference", "Bias review threshold", "Residual balance reference"),
      linetype = c("dashed", "solid", "dashed", "dashed"),
      role = c("threshold", "reference", "threshold", "reference")
    ),
    ranked = new_reference_lines(
      axis = c("v", "v", "v"),
      value = c(-thr$abs_bias_warn, 0, thr$abs_bias_warn),
      label = c("Bias review threshold", "Centered bias reference", "Bias review threshold"),
      linetype = c("dashed", "solid", "dashed"),
      role = c("threshold", "reference", "threshold")
    ),
    heatmap = new_reference_lines(
      axis = character(0), value = numeric(0),
      label = character(0), linetype = character(0),
      role = character(0)
    ),
    abs_t_hist = new_reference_lines("v", thr$abs_t_warn, "Screening |t| threshold", "dashed", "threshold"),
    facet_profile = new_reference_lines("v", thr$abs_bias_warn, "Mean |bias| review threshold", "dashed", "threshold")
  )
  heatmap_data <- build_bias_heatmap_data(scatter)
  plot_long <- build_bias_plot_long(
    ranked = ranked,
    scatter = scatter,
    profile = profile,
    heatmap_cells = heatmap_data$cells,
    heatmap_levels = list(x = heatmap_data$x_levels, y = heatmap_data$y_levels)
  )
  plot_annotations <- build_bias_plot_annotations(thr)
  flag_summary <- build_bias_flag_summary(scatter)
  plot_settings <- data.frame(
    Plot = plot,
    TopN = top_n,
    SortBy = sort_by,
    AbsTWarn = suppressWarnings(as.numeric(thr$abs_t_warn %||% abs_t_warn)),
    AbsBiasWarn = suppressWarnings(as.numeric(thr$abs_bias_warn %||% abs_bias_warn)),
    PMax = suppressWarnings(as.numeric(thr$p_max %||% p_max)),
    CIRequested = isTRUE(show_ci),
    CIAvailable = isTRUE(ci_available),
    CILevel = if (isTRUE(show_ci)) ci_level else NA_real_,
    Preset = style$name,
    stringsAsFactors = FALSE
  )

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    if (plot == "scatter") {
      scatter_plot <- scatter[is.finite(scatter$ObsExpAverage) & is.finite(scatter$BiasSize), , drop = FALSE]
      if (nrow(scatter_plot) == 0) {
        graphics::plot.new()
        graphics::title(main = plot_title)
        graphics::text(0.5, 0.5, "No data")
      } else {
        col <- ifelse(as.logical(scatter_plot$Flag), pal["flag"], pal["normal"])
        # Widen ylim when CI whiskers are drawn so the plot fits them.
        ylim <- if (ci_available && all(c("CI_Lower", "CI_Upper") %in% names(scatter_plot))) {
          range(c(scatter_plot$BiasSize, scatter_plot$CI_Lower,
                  scatter_plot$CI_Upper, -thr$abs_bias_warn,
                  thr$abs_bias_warn), finite = TRUE)
        } else {
          NULL
        }
        graphics::plot(
          x = scatter_plot$ObsExpAverage,
          y = scatter_plot$BiasSize,
          pch = 16,
          col = col,
          xlab = "Obs-Exp Average",
          ylab = "Bias Size (logits)",
          main = plot_title,
          ylim = ylim
        )
        if (ci_available && all(c("CI_Lower", "CI_Upper") %in% names(scatter_plot))) {
          valid <- is.finite(scatter_plot$CI_Lower) & is.finite(scatter_plot$CI_Upper)
          if (any(valid)) {
            graphics::segments(
              x0 = scatter_plot$ObsExpAverage[valid],
              y0 = scatter_plot$CI_Lower[valid],
              x1 = scatter_plot$ObsExpAverage[valid],
              y1 = scatter_plot$CI_Upper[valid],
              col = col[valid], lwd = 1
            )
          }
        }
        graphics::abline(h = pretty(graphics::par("usr")[3:4], n = 5), col = style$grid, lty = 1)
        graphics::abline(v = pretty(graphics::par("usr")[1:2], n = 5), col = style$grid, lty = 1)
        graphics::abline(h = c(-thr$abs_bias_warn, 0, thr$abs_bias_warn), lty = c(2, 1, 2), col = c(style$neutral, style$axis, style$neutral))
        graphics::abline(v = 0, lty = 2, col = style$neutral)
      }
    } else if (plot == "ranked") {
      ranked_plot <- ranked[is.finite(suppressWarnings(as.numeric(ranked$BiasSize))), , drop = FALSE]
      if (nrow(ranked_plot) == 0) {
        graphics::plot.new()
        graphics::title(main = plot_title)
        graphics::text(0.5, 0.5, "No data")
      } else {
        sub <- ranked_plot[seq_len(min(nrow(ranked_plot), top_n)), , drop = FALSE]
        y <- rev(seq_len(nrow(sub)))
        vals <- rev(suppressWarnings(as.numeric(sub$BiasSize)))
        lbl <- truncate_axis_label(rev(as.character(sub$Pair)), width = 28L)
        col <- ifelse(rev(as.logical(sub$Flag)), pal["flag"], pal["normal"])
        # CI whiskers need their own x-range so the plot fits them.
        have_sub_ci <- ci_available && all(c("CI_Lower", "CI_Upper") %in% names(sub))
        ci_lo <- if (have_sub_ci) rev(suppressWarnings(as.numeric(sub$CI_Lower))) else NULL
        ci_hi <- if (have_sub_ci) rev(suppressWarnings(as.numeric(sub$CI_Upper))) else NULL
        xlim <- if (have_sub_ci) {
          range(c(vals, ci_lo, ci_hi, -thr$abs_bias_warn, thr$abs_bias_warn),
                finite = TRUE)
        } else {
          NULL
        }
        graphics::plot(
          x = vals,
          y = y,
          type = "n",
          yaxt = "n",
          ylab = "",
          xlab = "Bias Size (logits)",
          main = plot_title,
          xlim = xlim
        )
        graphics::abline(v = pretty(graphics::par("usr")[1:2], n = 5), col = style$grid, lty = 1)
        graphics::segments(0, y, vals, y, col = style$neutral)
        if (have_sub_ci) {
          valid <- is.finite(ci_lo) & is.finite(ci_hi)
          if (any(valid)) {
            graphics::segments(
              x0 = ci_lo[valid], y0 = y[valid],
              x1 = ci_hi[valid], y1 = y[valid],
              col = col[valid], lwd = 2
            )
          }
        }
        graphics::points(vals, y, pch = 16, col = col)
        graphics::axis(side = 2, at = y, labels = lbl, las = las_rank, cex.axis = 0.75)
        graphics::abline(v = c(-thr$abs_bias_warn, 0, thr$abs_bias_warn), lty = c(2, 1, 2), col = c(style$neutral, style$axis, style$neutral))
      }
    } else if (plot == "heatmap") {
      # Cell-color heatmap of bias size for the two interaction facets.
      # Diverging palette, centred at zero. Flagged cells get a heavy
      # outline so they remain visible against the colour fill.
      mat_bias <- heatmap_data$matrix
      mat_flag <- heatmap_data$flag_matrix
      a_lvls <- heatmap_data$y_levels
      b_lvls <- heatmap_data$x_levels
      if (length(mat_bias) == 0L || length(a_lvls) == 0L ||
          length(b_lvls) == 0L) {
        graphics::plot.new()
        graphics::title(main = plot_title)
        graphics::text(0.5, 0.5, "No data")
      } else {
        # Diverging palette around zero with symmetric range.
        max_abs <- max(abs(mat_bias), na.rm = TRUE)
        if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
        cols <- grDevices::colorRampPalette(
          c("#2166AC", "#FFFFFF", "#B2182B")
        )(101)
        old_par <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(old_par), add = TRUE)
        graphics::par(mar = c(max(5, label_angle / 9 + 4),
                              max(5, max(nchar(a_lvls)) * 0.55 + 2),
                              3, 1))
        graphics::image(
          x = seq_len(length(b_lvls)),
          y = seq_len(length(a_lvls)),
          z = t(mat_bias),
          col = cols,
          zlim = c(-max_abs, max_abs),
          xaxt = "n", yaxt = "n",
          xlab = bundle$summary$InteractionFacets[1] %||% "",
          ylab = "",
          main = plot_title
        )
        graphics::axis(1, at = seq_len(length(b_lvls)), labels = b_lvls,
                       las = if (label_angle >= 45) 2 else 1,
                       cex.axis = 0.8)
        graphics::axis(2, at = seq_len(length(a_lvls)), labels = a_lvls,
                       las = 1, cex.axis = 0.8)
        # Flagged cells: thick outline so they stand out.
        for (i in seq_len(nrow(mat_flag))) {
          for (j in seq_len(ncol(mat_flag))) {
            if (isTRUE(mat_flag[i, j])) {
              graphics::rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5,
                             border = pal["flag"], lwd = 2)
            }
            v <- mat_bias[i, j]
            if (is.finite(v)) {
              graphics::text(j, i, sprintf("%.2f", v),
                             cex = 0.7,
                             col = if (abs(v) > max_abs * 0.6) "white" else "black")
            }
          }
        }
      }
    } else if (plot == "abs_t_hist") {
      tvals <- abs(suppressWarnings(as.numeric(scatter$t)))
      tvals <- tvals[is.finite(tvals)]
      if (length(tvals) == 0) {
        graphics::plot.new()
        graphics::title(main = plot_title)
        graphics::text(0.5, 0.5, "No data")
      } else {
        graphics::hist(
          x = tvals,
          breaks = "FD",
          col = pal["hist"],
          border = style$background,
          xlab = "|t|",
          main = plot_title
        )
        graphics::abline(v = pretty(graphics::par("usr")[1:2], n = 5), col = style$grid, lty = 1)
        graphics::abline(v = thr$abs_t_warn, lty = 2, col = style$neutral)
      }
    } else {
      if (nrow(profile) == 0 || !all(c("Facet", "Level", "MeanAbsBias", "FlagRate") %in% names(profile))) {
        graphics::plot.new()
        graphics::title(main = plot_title)
        graphics::text(0.5, 0.5, "No data")
      } else {
        prof <- profile |>
          dplyr::mutate(
            MeanAbsBias = suppressWarnings(as.numeric(.data$MeanAbsBias)),
            FlagRate = suppressWarnings(as.numeric(.data$FlagRate))
          ) |>
          dplyr::filter(is.finite(.data$MeanAbsBias)) |>
          dplyr::arrange(.data$Facet, dplyr::desc(.data$MeanAbsBias))
        if (nrow(prof) == 0) {
          graphics::plot.new()
          graphics::title(main = plot_title)
          graphics::text(0.5, 0.5, "No data")
        } else {
          lbl <- truncate_axis_label(paste0(prof$Facet, ": ", prof$Level), width = 34L)
          cols <- ifelse(is.finite(prof$FlagRate) & prof$FlagRate > 0, pal["flag"], pal["profile"])
          graphics::dotchart(
            x = prof$MeanAbsBias,
            labels = lbl,
            pch = 16,
            col = cols,
            cex = 0.8,
            cex.axis = 0.75,
            xlab = "Mean |Bias Size| (logits)",
            main = plot_title
          )
          graphics::abline(v = pretty(graphics::par("usr")[1:2], n = 5), col = style$grid, lty = 1)
          graphics::abline(v = thr$abs_bias_warn, lty = 2, col = style$neutral)
        }
      }
    }
  }

  out <- new_mfrm_plot_data(
    "table13_bias",
    list(
      plot = plot,
      ranked_table = ranked,
      scatter_data = scatter,
      facet_profile = profile,
      plot_long = plot_long,
      plot_annotations = plot_annotations,
      flag_summary = flag_summary,
      plot_settings = plot_settings,
      heatmap_cells = heatmap_data$cells,
      heatmap_matrix = heatmap_data$matrix,
      heatmap_flag_matrix = heatmap_data$flag_matrix,
      heatmap_count_matrix = heatmap_data$count_matrix,
      heatmap_levels = list(
        x = heatmap_data$x_levels,
        y = heatmap_data$y_levels,
        collapsed = heatmap_data$collapsed
      ),
      interpretation_guide = heatmap_data$interpretation_guide,
      reference_notes = heatmap_data$reference_notes,
      summary = bundle$summary,
      thresholds = thr,
      title = plot_title,
      subtitle = plot_subtitle,
      legend = plot_legend,
      reference_lines = plot_reference,
      preset = style$name
    )
  )
  invisible(out)
}

# Human-friendly API aliases (preferred names)
