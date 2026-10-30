#' Describe rater centrality using RSM residuals and rater-specific PCM steps
#'
#' Place indirect centrality indices alongside their model basis, exposure and
#' availability. These statistics describe fitted ratings; they do not classify
#' raters or estimate a centrality parameter.
#'
#' @param fit Optional native additive RSM fit from [fit_mfrm()], for residual
#'   versus expected-score correlations and ordinary Infit/Outfit.
#' @param threshold_fit Optional native additive PCM fit with `step_facet`
#'   equal to `rater_facet`, for the spread of rater-specific step estimates.
#'   At least one of `fit` and `threshold_fit` must be supplied.
#' @param rater_facet Name of the fitted non-Person rater facet.
#'
#' @details
#' Eckes and Jin (2021, pp. 68-70, 77, 83-85) distinguish RSM residual-based
#' indices from threshold SD under a rater-specific PCM. Accordingly, this
#' function never assigns a shared RSM threshold SD to individual raters or
#' substitutes a criterion-specific PCM for a rater-specific PCM. Supplying
#' only one model leaves the other indices unavailable with an explicit reason.
#' Both models must use the same retained response roster, including row
#' multiplicities, category mapping, facet signs, centering and anchors. Row
#' order may differ. Each model's estimator and readiness are retained.
#'
#' Current support is ordinary JML or MML, unit observation weights, additive
#' RSM/PCM and no facet shrinkage. Imported fits, interactions, generalized
#' slope models, corrected JML, testlet and shared-random-rater models require
#' their own definitions. There is no model fitting, quadrature or covariance
#' estimation here. RSM response expectations are evaluated at saved point
#' estimates using the same probability calculation as [diagnose_mfrm()].
#' MML uses saved Person EAP scores: these are plug-in expectations, not
#' posterior-integrated or posterior-predictive expectations. Estimation and
#' assignment effects can affect these indices even without true centrality.
#'
#' `ResidualExpectedCorrelation` is the unweighted Pearson correlation of raw
#' residuals (observed minus expected score) with expected scores within a rater.
#' More negative correlations are consistent with scores compressed toward the
#' middle. Standardized residuals are not used. A correlation requires at least
#' two complete rows with nonconstant residuals and expected scores. Missing or
#' nonfinite expectations, or a nonfinite primary Person/facet estimate involved
#' in a row, make that rater's residual indices unavailable; incomplete subsets
#' are not silently analyzed. Infit and Outfit additionally require strictly
#' positive finite response variances in every row. Infit is
#' `sum(residual^2) / sum(variance)`; Outfit is `mean(residual^2 / variance)`.
#'
#' `ThresholdSD` is the sample SD across the complete set of a rater's `m`
#' adjacent-category step estimates (divisor `m - 1`, explicitly chosen here;
#' the paper does not specify the divisor). `ThresholdVariance` is its square.
#' Larger spreads are consistent with centrality and smaller spreads with
#' extremity. Steps are centered offsets in logits, not cumulative thresholds;
#' adding the same rater severity/location to each step leaves their SD
#' unchanged. At least two finite steps are required: a binary rubric cannot
#' supply this spread. Steps are kept in category order. Disorder is reported,
#' not repaired by sorting estimates; spread alone does not establish proper
#' category functioning. `ObservedSD` describes the fitted internal score
#' coding, with original-category labels retained in `categories`.
#'
#' Infit/Outfit are supporting fit descriptions. Low Infit can also reflect
#' halo, restricted score ranges or other mechanisms. No default .75, -.30 or
#' 1.40 cutoff, combined flag, significance test or confidence interval is
#' supplied. Threshold uncertainty, sparse exposure, different ability mixes,
#' criterion mixtures and dependence between ratings affect interpretation.
#' This function does not implement the Bayesian FM-SC centrality parameter
#' omega, its uncertainty or ability adjustments, and does not reproduce FACETS
#' numerical output merely by using related formulas. The package's `ej2021`
#' example data are synthetic, not the original TestDaF samples.
#'
#' In the direct FM-SC model, rater severity shifts category locations while
#' a positive rater-specific omega multiplies category step offsets. Omega
#' therefore changes the response model and its ability estimates; it is not
#' the SD of separately fitted PCM steps or a GPCM discrimination coefficient.
#' Jin and Wang (2018) introduce this model; Eckes and Jin (2021, p. 71) give
#' its RSM formulation. They are references for the distinction, not additional
#' estimation methods available through this function.
#'
#' @return A list containing `table` (one row per rater), `thresholds`
#'   (category-ordered PCM step estimates), `categories` (counts on the retained
#'   roster), `sources` (model settings, estimates and readiness), `settings`
#'   and `notes`. Each index has an availability reason; source readiness is
#'   separate from numerical availability. Save with [base::saveRDS()] and
#'   reopen with [base::readRDS()]; export the data frames with
#'   [utils::write.csv()]. Reopening the result does not recalculate indices.
#'   This plain list has no dedicated `plot()` or automatic aggregate-report
#'   section. For a custom figure, retain each index's reason and source-status
#'   columns; do not mix threshold SD, residual correlation and fit on one scale.
#'   Graph centrality from [rater_halo_network_analysis()] has a different meaning.
#'
#' @references
#' Eckes, T., & Jin, K.-Y. (2021). Measuring rater centrality effects in writing
#' assessment: A Bayesian facets modeling approach. *Psychological Test and
#' Assessment Modeling, 63*(1), 65-94. Definitions and actual analyses:
#' pp. 68-70 and 77; comparisons and limitations: pp. 83-88.
#' [Full text](https://www.psychologie-aktuell.com/fileadmin/download/ptam/1-2021/Seiten_aus_PTAM_2021-1_ebook_4.pdf).
#'
#' Jin, K.-Y., & Wang, W.-C. (2018). A new facets model for rater's
#' centrality/extremity response style. *Journal of Educational Measurement,
#' 55*(4), 543-563. \doi{10.1111/jedm.12191}. Direct model, not implemented here.
#' @seealso [score_comparison_table()], [interrater_agreement_table()],
#'   [diagnose_mfrm()], [mfrm_report()], [mfrmr_methods_reference]
#' @examples
#' ratings <- load_mfrmr_data("example_core")
#' rsm <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", model = "RSM", maxit = 100)
#' pcm <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", model = "PCM", step_facet = "Rater",
#'                 maxit = 100)
#' centrality <- rater_centrality_table(rsm, pcm, rater_facet = "Rater")
#' centrality$table
#' centrality$thresholds
#' # Review both SourceStatus columns before interpreting finite indices.
#' @export
rater_centrality_table <- function(fit = NULL, threshold_fit = NULL,
                                   rater_facet = "Rater") {
  if (is.null(fit) && is.null(threshold_fit)) stop(
    "Supply an RSM `fit`, a rater-specific PCM `threshold_fit`, or both.", call. = FALSE)
  if (!is.character(rater_facet) || length(rater_facet) != 1L ||
      is.na(rater_facet) || !nzchar(rater_facet)) stop(
    "`rater_facet` must be one nonempty facet name.", call. = FALSE)
  if (!is.null(fit)) mfrm_centrality_source_check(fit, "RSM", rater_facet)
  if (!is.null(threshold_fit)) mfrm_centrality_source_check(threshold_fit, "PCM", rater_facet)
  if (!is.null(fit) && !is.null(threshold_fit)) {
    a <- normalize_compare_signature(fit)
    b <- normalize_compare_signature(threshold_fit)
    fields <- c("facets", "rating_min", "rating_max", "score_map",
                "noncenter_facet", "dummy_facets", "anchors", "group_anchors")
    if (!all(vapply(fields, function(nm) same_signature_component(a[[nm]], b[[nm]]), logical(1))) ||
        !identical(fit$config$facet_signs[sort(fit$config$facet_names)],
                   threshold_fit$config$facet_signs[sort(threshold_fit$config$facet_names)]) ||
        !identical(mfrm_centrality_roster(fit), mfrm_centrality_roster(threshold_fit))) stop(
      "RSM and PCM must use the same retained response roster, category mapping and facet constraints.",
      call. = FALSE)
  }
  source <- fit %||% threshold_fit
  data <- source$prep$data
  raters <- as.character(source$prep$levels[[rater_facet]])
  status <- function(x) {
    if (is.null(x)) "not_supplied" else if (isTRUE(.mfrm_fit_plot_readiness(x)$ready))
      "descriptive" else "source_needs_review"
  }
  obs <- if (!is.null(fit)) compute_obs_table(fit) else NULL
  if (!is.null(obs)) {
    # Boundary Person values may be displayed as infinite although the
    # optimizer retains a finite proxy. Do not use that proxy as a valid row.
    valid_primary <- is.finite(fit$facets$person$Estimate[
      match(as.character(obs$Person), as.character(fit$facets$person$Person))])
    for (facet in fit$config$facet_names) {
      coef <- fit$facets$others[fit$facets$others$Facet == facet, , drop = FALSE]
      valid_primary <- valid_primary & is.finite(coef$Estimate[
        match(as.character(obs[[facet]]), as.character(coef$Level))])
    }
    obs$.CentralityAvailable <- valid_primary & is.finite(obs$Observed) &
      is.finite(obs$Expected) & is.finite(obs$Residual)
    fit_stats <- calc_facet_fit(obs, rater_facet)
  }
  steps <- if (!is.null(threshold_fit)) as.data.frame(threshold_fit$steps) else NULL
  categories <- source$prep$score_map
  category_rows <- vector("list", length(raters))
  rows <- lapply(seq_along(raters), function(i) {
    rater <- raters[i]
    take <- which(as.character(data[[rater_facet]]) == rater)
    score <- data$Score[take]
    row <- data.frame(Rater = rater, RatingRows = length(take),
      Persons = length(unique(data$Person[take])), ObservedSD = if (length(take) >= 2L) stats::sd(score) else NA_real_,
      ResidualExpectedCorrelation = NA_real_, Infit = NA_real_, Outfit = NA_real_,
      ResidualRows = 0L, ResidualUnavailableRows = length(take),
      ResidualReason = "RSM fit not supplied.", FitReason = "RSM fit not supplied.",
      ThresholdSD = NA_real_, ThresholdVariance = NA_real_, ThresholdCount = 0L,
      ThresholdsOrdered = NA, ThresholdReason = "Rater-specific PCM not supplied.",
      ResidualMethod = if (!is.null(fit)) fit$config$method else NA_character_,
      ThresholdMethod = if (!is.null(threshold_fit)) threshold_fit$config$method else NA_character_,
      ResidualSourceStatus = status(fit), ThresholdSourceStatus = status(threshold_fit))
    if (!is.null(obs)) {
      at <- which(as.character(obs[[rater_facet]]) == rater)
      r <- obs[at, , drop = FALSE]
      row$ResidualRows <- sum(r$.CentralityAvailable)
      row$ResidualUnavailableRows <- nrow(r) - row$ResidualRows
      row$ResidualReason <- if (nrow(r) < 2L) "Fewer than two retained ratings." else
        if (any(!r$.CentralityAvailable)) "Incomplete or nonfinite primary estimates/residuals." else
          if (stats::sd(r$Expected) == 0 || stats::sd(r$Residual) == 0)
            "Expected scores or raw residuals are constant." else "available"
      if (row$ResidualReason == "available") {
        row$ResidualExpectedCorrelation <- stats::cor(r$Residual, r$Expected)
        if (!is.finite(row$ResidualExpectedCorrelation)) {
          row$ResidualExpectedCorrelation <- NA_real_
          row$ResidualReason <- "Correlation is numerically unavailable."
        }
      }
      row$FitReason <- if (!nrow(r) || any(!r$.CentralityAvailable) ||
        any(!is.finite(r$Var) | r$Var <= 0)) "Incomplete residuals or nonpositive response variance." else "available"
      if (row$FitReason == "available") {
        j <- match(rater, as.character(fit_stats$Level))
        row$Infit <- fit_stats$Infit[j]; row$Outfit <- fit_stats$Outfit[j]
        if (!is.finite(row$Infit) || !is.finite(row$Outfit)) {
          row$Infit <- row$Outfit <- NA_real_
          row$FitReason <- "Fit statistics are numerically unavailable."
        }
      }
    }
    if (!is.null(steps)) {
      j <- which(as.character(steps$StepFacet) == rater)
      s <- steps$Estimate[j]
      coef <- threshold_fit$facets$others
      primary <- coef$Estimate[coef$Facet == rater_facet & as.character(coef$Level) == rater]
      row$ThresholdCount <- length(s)
      row$ThresholdReason <- if (length(s) < 2L) "Fewer than two thresholds (binary or shorter rubric)." else
        if (any(!is.finite(s)) || length(primary) != 1L || !is.finite(primary))
          "Incomplete or nonfinite primary rater/threshold estimates." else "available"
      if (row$ThresholdReason == "available") {
        row$ThresholdSD <- stats::sd(s)
        row$ThresholdVariance <- stats::var(s)
        row$ThresholdsOrdered <- all(diff(s) > 0)
        if (!is.finite(row$ThresholdSD) || !is.finite(row$ThresholdVariance)) {
          row$ThresholdSD <- row$ThresholdVariance <- NA_real_
          row$ThresholdReason <- "Threshold spread is numerically unavailable."
        }
      }
    }
    counts <- vapply(categories$InternalScore, function(k) sum(score == k), integer(1))
    category_rows[[i]] <<- data.frame(Rater = rater,
      ModelScore = categories$InternalScore, OriginalScore = categories$OriginalScore,
      Ratings = counts, Proportion = if (length(take)) counts / length(take) else NA_real_)
    row
  })
  if (!is.null(steps)) {
    names(steps)[names(steps) == "StepFacet"] <- "Rater"
    steps$SDDivisor <- threshold_fit$config$n_cat - 2L
  } else steps <- data.frame(Rater = character(), Step = character(), Estimate = numeric(), SDDivisor = integer())
  source_record <- function(x) if (is.null(x)) NULL else list(
    config = x$config, steps = x$steps, facets = x$facets,
    readiness = mfrmr_get_readiness_record(x), source_status = status(x))
  list(table = do.call(rbind, rows), thresholds = steps,
    categories = do.call(rbind, category_rows),
    sources = list(rsm = source_record(fit), pcm = source_record(threshold_fit)),
    settings = list(rater_facet = rater_facet, weights = "unit", sd_divisor = "count_minus_one",
      residual_basis = "RSM_point_score_plugin_raw_residuals", thresholds_basis = "rater_specific_PCM_steps",
      matched_roster = !is.null(fit) && !is.null(threshold_fit),
      reference = "Eckes & Jin (2021), pp. 68-70, 77, 83-88"),
    notes = c(
      "Indirect descriptive centrality indices, not an estimated FM-SC omega, diagnostic classification or adjustment of ability estimates.",
      "RSM residual-expected correlations use raw residuals and point-score plug-in expectations. MML uses saved Person EAP scores, not posterior integration.",
      "Larger PCM threshold SDs and more negative RSM residual-expected correlations can suggest centrality. Exposure, ability/criterion mixtures, threshold error and rating dependence affect interpretation.",
      "Infit/Outfit are supporting descriptions; low fit does not uniquely identify centrality. No automatic cutoffs, tests or confidence intervals are supplied.",
      "Threshold SD uses m-1 across m complete steps. Shared RSM and non-rater PCM thresholds cannot provide this index. Disordered steps are retained and marked.",
      "Finite values from a source needing review remain provisional. Numerical availability does not qualify diagnostic accuracy or establish equivalence with FACETS."))
}

mfrm_centrality_source_check <- function(fit, model, facet) {
  cfg <- fit$config
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      !identical(cfg$model, model) || !isTRUE(cfg$method %in% c("JML", "MML")) ||
      length(cfg$interaction_specs) || length(fit$interactions$specs) ||
      mfrm_has_jml_adjustment(fit) || !identical(cfg$facet_shrinkage %||% "none", "none")) stop(
    "Centrality indices require a native additive ", model, " JML/MML fit without shrinkage or interactions.", call. = FALSE)
  if (!facet %in% cfg$facet_names || identical(facet, "Person")) stop(
    "`rater_facet` must identify a fitted non-Person facet.", call. = FALSE)
  if (model == "PCM" && !identical(cfg$step_facet, facet)) stop(
    "`threshold_fit` must have PCM steps owned by `rater_facet`.", call. = FALSE)
  d <- fit$prep$data
  if (!is.data.frame(d) || !nrow(d) || !"Weight" %in% names(d) ||
      any(!is.finite(d$Weight) | d$Weight != 1)) stop(
    "Centrality indices currently require retained ratings with unit weights.", call. = FALSE)
  # Compare by the model's category order, not sorted step estimates.
  params <- expand_params(fit$opt$par, build_param_sizes(cfg), cfg)
  expected <- build_step_table(cfg, fit$prep, params)
  columns <- intersect(c("StepFacet", "Step", "Estimate"), names(expected))
  if (!is.data.frame(fit$steps) || !all(columns %in% names(fit$steps)) ||
      !isTRUE(all.equal(as.data.frame(fit$steps[columns]), as.data.frame(expected[columns]),
                         tolerance = 1e-10, check.attributes = FALSE))) stop(
    "Saved thresholds must match the fitted parameters and category order.", call. = FALSE)
  invisible(TRUE)
}

mfrm_centrality_roster <- function(fit) {
  columns <- c("Person", sort(fit$config$facet_names), "Score", "score_k", "Weight")
  d <- as.data.frame(fit$prep$data[columns])
  d[] <- lapply(d, function(x) if (is.factor(x)) as.character(x) else x)
  d <- d[do.call(order, c(d, list(na.last = TRUE, method = "radix"))), , drop = FALSE]
  rownames(d) <- NULL
  d
}
