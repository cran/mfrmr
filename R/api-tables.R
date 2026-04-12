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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' ir <- interrater_agreement_table(fit, rater_facet = "Rater")
#' summary(ir)$summary
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
#'   fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 10)
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
#'
#' @details
#' This function wraps the package's adjusted-score calculations and returns
#' both facet-wise and stacked tables. Historical display columns such as
#' `Fair(M) Average` and `Fair(Z) Average` are retained for compatibility, and
#' package-native aliases such as `AdjustedAverage`,
#' `StandardizedAdjustedAverage`, `ModelBasedSE`, and `FitAdjustedSE` are
#' appended to the formatted outputs.
#'
#' In the current release, these tables are source-backed only for the
#' Rasch-family `RSM` / `PCM` branch. FACETS documents fair averages as
#' Rasch-measure-to-score transformations evaluated in a standardized
#' mean/zero-facet environment. The bounded `GPCM` branch already has a
#' generalized ordered-category probability kernel, but this package has not
#' yet validated a slope-aware analogue of that fair-average score contract.
#' `fair_average_table()` therefore stops for `GPCM` fits instead of silently
#' reusing the Rasch-only calculation.
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
#'   \item{Measure}{Estimated logit measure for this level.}
#'   \item{SE}{Compatibility alias for the model-based standard error.}
#'   \item{ModelBasedSE, FitAdjustedSE}{Package-native aliases for `Model S.E.` and `Real S.E.`.}
#'   \item{Infit MnSq, Outfit MnSq}{Fit statistics for this level.}
#' }
#'
#' @return A named list with:
#' - `by_facet`: named list of formatted data.frames
#' - `stacked`: one stacked data.frame across facets
#' - `raw_by_facet`: unformatted internal tables
#' - `settings`: resolved options
#'
#' @seealso [diagnose_mfrm()], [unexpected_response_table()], [displacement_table()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t12 <- fair_average_table(fit, udecimals = 2)
#' t12_native <- fair_average_table(fit, reference = "mean", label_style = "native")
#' summary(t12)
#' p_t12 <- plot(t12, draw = FALSE)
#' p_t12$data$plot
#' }
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
                               xtreme = 0) {
  reference <- match.arg(reference)
  label_style <- match.arg(label_style)
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }
  fit_model <- as.character(fit$config$model %||% fit$summary$Model[1] %||% NA_character_)
  if (identical(fit_model, "GPCM")) {
    stop(
      "`fair_average_table()` is not yet validated for `GPCM` fits. ",
      gpcm_fair_average_rationale(),
      call. = FALSE
    )
  }
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
    xtreme = xtreme
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
    xtreme = xtreme
  )
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
#'
#' @seealso [diagnose_mfrm()], [measurable_summary_table()], [plot.mfrm_fit()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t8 <- rating_scale_table(fit)
#' summary(t8)
#' summary(t8)$summary
#' p_t8 <- plot(t8, draw = FALSE)
#' p_t8$data$plot
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
#'   \item{FacetA, FacetB}{Interaction facet level identifiers; placeholder
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
#' @param data Optional raw data frame used for additional row-level audit.
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
#' The `row_audit` data.frame contains:
#' \describe{
#'   \item{Status}{Row-status category (e.g., "Valid", "MissingScore").}
#'   \item{N}{Number of rows in this category.}
#' }
#'
#' @return A named list with:
#' - `summary`: one-row data summary
#' - `model_match`: model-match counts
#' - `row_audit`: counts by row-status category
#' - `unknown_elements`: elements in `data` not present in fitted levels
#' - `category_counts`: observed category counts in fitted data
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#'
#' @seealso [fit_mfrm()], [specifications_report()], [describe_mfrm_data()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
                                include_fixed = FALSE) {
  signal_legacy_name_deprecation(
    old_name = "table2_data_summary",
    new_name = "data_quality_report",
    suppress_if_called_from = "data_quality_report"
  )
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().")
  }

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
  category_counts <- fitted_df |>
    group_by(.data$Score) |>
    summarize(
      Count = n(),
      WeightedCount = sum(.data$Weight, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(.data$Score)

  if (is.null(data)) {
    summary_tbl <- data.frame(
      TotalLinesInData = NA_integer_,
      TotalDataLines = NA_integer_,
      TotalNonBlankResponsesFound = NA_integer_,
      ValidResponsesUsedForEstimation = valid_used,
      stringsAsFactors = FALSE
    )
    model_match <- data.frame(
      Model = paste(cfg$model, cfg$method, sep = " / "),
      MatchedResponses = valid_used,
      stringsAsFactors = FALSE
    )
    out <- list(
      summary = summary_tbl,
      model_match = model_match,
      row_audit = data.frame(),
      unknown_elements = data.frame(),
      category_counts = as.data.frame(category_counts, stringsAsFactors = FALSE)
    )
    if (isTRUE(include_fixed)) {
      out$fixed <- build_sectioned_fixed_report(
        title = "Legacy-compatible Table 2 Data Summary",
        sections = list(
          list(title = "Summary", data = summary_tbl),
          list(title = "Model match", data = model_match),
          list(title = "Row audit", data = data.frame()),
          list(title = "Unknown elements", data = data.frame()),
          list(title = "Category counts", data = as.data.frame(category_counts, stringsAsFactors = FALSE))
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
  range_ok <- score_nonblank & score_num >= prep$rating_min & score_num <= prep$rating_max

  row_status <- rep("valid", n_total)
  row_status[!person_ok] <- "missing_person"
  row_status[person_ok & !facet_ok] <- "missing_facet"
  row_status[person_ok & facet_ok & !score_nonblank] <- "missing_score"
  row_status[person_ok & facet_ok & score_nonblank & !weight_ok] <- "invalid_weight"
  row_status[person_ok & facet_ok & score_nonblank & weight_ok & !range_ok] <- "score_out_of_range"

  row_audit <- data.frame(Status = row_status, stringsAsFactors = FALSE) |>
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
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary_tbl,
    model_match = model_match,
    row_audit = row_audit,
    unknown_elements = unknown_tbl,
    category_counts = as.data.frame(category_counts, stringsAsFactors = FALSE)
  )
  if (isTRUE(include_fixed)) {
    out$fixed <- build_sectioned_fixed_report(
      title = "Legacy-compatible Table 2 Data Summary",
      sections = list(
        list(title = "Summary", data = summary_tbl),
        list(title = "Model match", data = model_match),
        list(title = "Row audit", data = row_audit),
        list(title = "Unknown elements", data = unknown_tbl),
        list(title = "Category counts", data = as.data.frame(category_counts, stringsAsFactors = FALSE))
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
  idx <- build_indices(prep, step_facet = cfg$step_facet)
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' t6 <- subset_connectivity_report(fit)
#' @keywords internal
#' @noRd
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
#' - `expected_ogive`: long expected-score table
#' - `fixed`: fixed-width report text (when `include_fixed = TRUE`)
#' - `settings`: applied options
#'
#' @seealso [category_structure_report()], [rating_scale_table()], [plot.mfrm_fit()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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

#' Build a legacy-compatible output-file bundle (`GRAPH=` / `SCORE=`)
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()] (used for score file).
#' @param include Output components to include: `"graph"` and/or `"score"`.
#' @param theta_range Theta/logit range for graph coordinates.
#' @param theta_points Number of points on the theta grid for graph coordinates.
#' @param digits Rounding digits for numeric fields.
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
#' `"obs_probability"`).
#'
#' @section Interpreting output:
#' - `graphfile`: legacy-compatible wide curve coordinates (human-readable labels).
#' - `graphfile_syntactic`: same curves with syntactic column names for programmatic use.
#' - `scorefile`: observation-level observed/expected/residual diagnostics.
#' - `written_files`: audit trail of files produced when `write_files = TRUE`.
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
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
  fit_model <- as.character(fit$config$model %||% fit$summary$Model[1] %||% NA_character_)
  if (identical(fit_model, "GPCM") && "score" %in% include) {
    stop_if_gpcm_out_of_scope(
      fit,
      "facets_output_file_bundle(include = \"score\")",
      supported = "fitting, core summary output, fixed-calibration posterior scoring, compute_information(), pathway/CCC plotting, category curve/structure reports, and graph-only output bundles"
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
      "StdResidual", "Var", "Weight", "score_k", "PersonMeasure"
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
    }
    out$scorefile <- round_numeric_df(scorefile, digits = digits)
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
    overwrite = overwrite
  )
  as_mfrm_bundle(out, "mfrm_output_bundle")
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
  if (!is.null(facet)) out$Facet <- facet
  out
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
#' Output tables use:
#' - `Component`: principal-component index (1, 2, ...)
#' - `Eigenvalue`: eigenvalue for each component
#' - `Proportion`: component variance proportion
#' - `Cumulative`: cumulative variance proportion
#'
#' For `mode = "facet"` or `"both"`, `by_facet_table` additionally includes
#' a `Facet` column.
#'
#' `summary(pca)` is supported through `summary()`.
#' `plot(pca)` is dispatched through `plot()` for class
#' `mfrm_residual_pca`. Available types include `"overall_scree"`,
#' `"facet_scree"`, `"overall_loadings"`, and `"facet_loadings"`.
#'
#' @section Interpreting output:
#' Use `overall_table` first:
#' - early components with noticeably larger eigenvalues or proportions
#'   suggest stronger residual structure that may deserve follow-up.
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
#' - Linacre, J. M. (1998). *Structure in Rasch residuals: Why principal
#'   components analysis (PCA)?* Rasch Measurement Transactions, 12(2), 636.
#' - Linacre, J. M. (1998). *Detecting multidimensionality: Which residual
#'   data-type works best?* Journal of Outcome Measurement, 2(3), 266-283.
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
#'
#' @seealso [diagnose_mfrm()], [plot_residual_pca()], [mfrmr_visual_diagnostics]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "both")
#' pca <- analyze_residual_pca(diag, mode = "both")
#' pca2 <- analyze_residual_pca(fit, mode = "both")
#' summary(pca)
#' p <- plot_residual_pca(pca, mode = "overall", plot_type = "scree", draw = FALSE)
#' p$data$plot
#' head(p$data)
#' head(pca$overall_table)
#' }
#' @export
analyze_residual_pca <- function(diagnostics,
                                 mode = c("overall", "facet", "both"),
                                 facets = NULL,
                                 pca_max_factors = 10L) {
  mode <- match.arg(tolower(mode), c("overall", "facet", "both"))

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

  out <- list(
    mode = mode,
    facet_names = facet_names,
    overall = out_overall,
    by_facet = out_by_facet,
    overall_table = overall_table,
    by_facet_table = by_facet_table,
    errors = pca_errors
  )
  as_mfrm_bundle(out, "mfrm_residual_pca")
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
#' @param plot_type `"scree"` or `"loadings"`.
#' @param component Component index for loadings plot.
#' @param top_n Maximum number of variables shown in loadings plot.
#' @param preset Visual preset (`"standard"`, `"publication"`, or `"compact"`).
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
#' - `"loadings"`: horizontal bar chart of top absolute loadings
#'
#' For `mode = "facet"` and `facet = NULL`, the first available facet is used.
#'
#' @section Interpreting output:
#' - `plot_type = "scree"`: look for dominant early components relative
#'   to later components and the unit-eigenvalue reference line. Treat
#'   this as exploratory residual-structure screening, not a standalone
#'   unidimensionality test.
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
#' - `plot`: `"scree"` or `"loadings"`
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
#'   fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 15)
#' )
#' diag <- diagnose_mfrm(fit, residual_pca = "overall")
#' pca <- analyze_residual_pca(diag, mode = "overall")
#' plt <- plot_residual_pca(pca, mode = "overall", plot_type = "scree", draw = FALSE)
#' head(plt$data)
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
                              plot_type = c("scree", "loadings"),
                              component = 1L,
                              top_n = 20L,
                              preset = c("standard", "publication", "compact"),
                              draw = TRUE) {
  mode <- match.arg(tolower(mode), c("overall", "facet"))
  plot_type <- match.arg(tolower(plot_type), c("scree", "loadings"))
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

  if (plot_type == "scree") {
    tbl <- build_pca_variance_table(pca_bundle)
    if (nrow(tbl) == 0) stop("No eigenvalues available for scree plot.")
    title <- paste0(title_suffix, " (Scree)")
    subtitle <- if (mode == "overall") {
      "Variance explained by residual components"
    } else {
      paste0("Facet-specific scree profile: ", facet)
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
      graphics::abline(h = 1, lty = 2, col = style$neutral)
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
          label = c("Residual eigenvalues", "Unit-eigenvalue reference"),
          role = c("component", "reference"),
          aesthetic = c("line-point", "line"),
          value = c(style$accent_primary, style$neutral)
        ),
        reference_lines = new_reference_lines("h", 1, "Unit-eigenvalue reference", "dashed", "reference"),
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

#' Estimate legacy-compatible bias/interaction terms iteratively
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param facet_a First facet name.
#' @param facet_b Second facet name.
#' @param interaction_facets Character vector of two or more facets to model as
#'   one interaction effect. When supplied, this takes precedence over
#'   `facet_a`/`facet_b`.
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
#' The function estimates \eqn{b_{jc}} from the residuals of the fitted
#' (additive) model using iterative recalibration in a legacy-compatible style
#' (Myford & Wolfe, 2003, 2004):
#'
#' \deqn{b_{jc} = \frac{\sum_n (X_{njc} - E_{njc})}
#'                     {\sum_n \mathrm{Var}_{njc}}}
#'
#' Each iteration updates expected scores using the current bias
#' estimates, then re-computes the bias.  Convergence is reached when
#' the maximum absolute change in bias estimates falls below `tol`.
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
#'   reporting-use flags, and fit columns
#' - `summary`: compact summary statistics
#' - `chi_sq`: fixed-effect chi-square style screening summary
#' - `facet_a`, `facet_b`: first two analyzed facet names (legacy compatibility)
#' - `interaction_facets`, `interaction_order`, `interaction_mode`: full
#'   interaction metadata
#' - `iteration`: iteration history/metadata
#'
#' @seealso [build_fixed_reports()], [build_apa_outputs()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
#' summary(bias)
#' p_bias <- plot_bias_interaction(bias, draw = FALSE)
#' p_bias$data$plot
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
  stop_if_gpcm_out_of_scope(fit, "estimate_bias()")
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
  }
  out
}

#' Build legacy-compatible fixed-width text reports
#'
#' @param bias_results Output from [estimate_bias()].
#' @param target_facet Optional target facet for pairwise contrast table.
#' @param branch Output branch:
#'   `"facets"` keeps the legacy-compatible fixed-width layout;
#'   `"original"` returns compact sectioned fixed-width text for internal reporting.
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
#' - `pairwise_table`: machine-readable contrast table.
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
#' A named list with:
#' - `bias_fixed`: fixed-width interaction table text
#' - `pairwise_fixed`: fixed-width pairwise contrast text
#' - `pairwise_table`: underlying pairwise data.frame
#'
#' @seealso [estimate_bias()], [build_apa_outputs()], [bias_interaction_report()],
#'   [mfrmr_reports_and_tables], [mfrmr_compatibility_layer]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
#' - `summary`: one-row overview
#' - `thresholds`: applied cutoffs
#' - `facet_a`, `facet_b`: first two analyzed facet names
#' - `interaction_facets`, `interaction_order`, `interaction_mode`: full
#'   interaction metadata
#'
#' @seealso [estimate_bias()], [plot_bias_interaction()], [build_fixed_reports()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
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
        ObservedCount = if ("Observd Count" %in% names(tbl)) suppressWarnings(as.numeric(.data$`Observd Count`)) else NA_real_
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

  summary_tbl <- data.frame(
    InteractionFacets = paste(interaction_facets, collapse = " x "),
    InteractionOrder = interaction_order,
    InteractionMode = interaction_mode,
    FacetA = facet_a,
    FacetB = facet_b,
    Cells = nrow(tbl2),
    Flagged = sum(tbl2$Flag, na.rm = TRUE),
    FlaggedPercent = ifelse(nrow(tbl2) > 0, 100 * sum(tbl2$Flag, na.rm = TRUE) / nrow(tbl2), NA_real_),
    MeanAbsT = mean(tbl2$AbsT, na.rm = TRUE),
    MeanAbsBias = mean(tbl2$AbsBias, na.rm = TRUE),
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
    orientation_audit = as.data.frame(bias_results$orientation_audit %||% data.frame(), stringsAsFactors = FALSE),
    mixed_sign = isTRUE(bias_results$mixed_sign),
    direction_note = as.character(bias_results$direction_note %||% ""),
    recommended_action = as.character(bias_results$recommended_action %||% "")
  )
}

#' Plot interaction-bias diagnostics using base R
#'
#' @param x Output from [bias_interaction_report()], [estimate_bias()], or [fit_mfrm()].
#' @param plot Plot type: `"scatter"`, `"ranked"`, `"abs_t_hist"`, or
#'   `"facet_profile"`.
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
#' @param preset Visual preset (`"standard"`, `"publication"`, or `"compact"`).
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @section Lifecycle:
#' Soft-deprecated. Prefer [plot_bias_interaction()].
#'
#' @return A plotting-data object of class `mfrm_plot_data`.
#' @seealso [bias_interaction_report()], [estimate_bias()], [plot_displacement()]
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
#' p13 <- plot_bias_interaction(
#'   fit,
#'   diagnostics = diagnose_mfrm(fit, residual_pca = "none"),
#'   facet_a = "Rater",
#'   facet_b = "Criterion",
#'   draw = FALSE
#' )
#' @keywords internal
#' @noRd
plot_table13_bias <- function(x,
                              plot = c("scatter", "ranked", "abs_t_hist", "facet_profile"),
                              diagnostics = NULL,
                              facet_a = NULL,
                              facet_b = NULL,
                              interaction_facets = NULL,
                              top_n = 40,
                              abs_t_warn = 2,
                              abs_bias_warn = 0.5,
                              p_max = 0.05,
                              sort_by = c("abs_t", "abs_bias", "prob"),
                              main = NULL,
                              palette = NULL,
                              label_angle = 45,
                              preset = c("standard", "publication", "compact"),
                              draw = TRUE) {
  signal_legacy_name_deprecation(
    old_name = "plot_table13_bias",
    new_name = "plot_bias_interaction",
    suppress_if_called_from = "plot_bias_interaction"
  )
  plot <- match.arg(tolower(plot), c("scatter", "ranked", "abs_t_hist", "facet_profile"))
  sort_by <- match.arg(tolower(sort_by), c("abs_t", "abs_bias", "prob"))
  top_n <- max(1L, as.integer(top_n))
  label_angle <- suppressWarnings(as.numeric(label_angle[1]))
  if (!is.finite(label_angle)) label_angle <- 45
  las_rank <- if (label_angle >= 45) 2 else 1
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
  plot_title <- switch(
    plot,
    scatter = "Bias interaction scatter",
    ranked = "Ranked bias interaction size",
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
    abs_t_hist = new_reference_lines("v", thr$abs_t_warn, "Screening |t| threshold", "dashed", "threshold"),
    facet_profile = new_reference_lines("v", thr$abs_bias_warn, "Mean |bias| review threshold", "dashed", "threshold")
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
        graphics::plot(
          x = scatter_plot$ObsExpAverage,
          y = scatter_plot$BiasSize,
          pch = 16,
          col = col,
          xlab = "Obs-Exp Average",
          ylab = "Bias Size (logits)",
          main = plot_title
        )
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
        graphics::plot(
          x = vals,
          y = y,
          type = "n",
          yaxt = "n",
          ylab = "",
          xlab = "Bias Size (logits)",
          main = plot_title
        )
        graphics::abline(v = pretty(graphics::par("usr")[1:2], n = 5), col = style$grid, lty = 1)
        graphics::segments(0, y, vals, y, col = style$neutral)
        graphics::points(vals, y, pch = 16, col = col)
        graphics::axis(side = 2, at = y, labels = lbl, las = las_rank, cex.axis = 0.75)
        graphics::abline(v = c(-thr$abs_bias_warn, 0, thr$abs_bias_warn), lty = c(2, 1, 2), col = c(style$neutral, style$axis, style$neutral))
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
