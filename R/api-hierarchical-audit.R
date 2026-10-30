# ==============================================================================
# Hierarchical structure and small-sample review
# ==============================================================================
#
# Background and literature:
#
# - Linacre (2026), A User's Guide to FACETS, notes that rater estimates are
#   "more sensitive to link reductions" than examinee or task estimates.
#   Rasch sample-size guidelines (Linacre, 1994, 2021) recommend:
#     * >= 10 observations per category for stable scale reporting
#     * ~30 persons for +-1.0 logit stability at 95% CI
#     * ~100 persons for +-0.5 logit stability at 95% CI
#     * 250+ for high-stakes or published item calibrations
#   mfrmr applies these numerical bands to facet elements as well, since a
#   facet element with < 10 ratings is essentially a pilot-data level.
#
# - Myford & Wolfe (2004) Part II classified rater effects as
#   severity/leniency, central tendency, randomness (inaccuracy), halo,
#   and differential severity/leniency. This module adds only the review layer
#   needed to screen adequacy; bias screening for central tendency / halo
#   remains out of the current fit_mfrm() surface.
#
# - McEwen (2018) decomposed incomplete rating designs along four
#   attributes: rater coverage, repetition size, design structure, and
#   rater order. The `analyze_hierarchical_structure()` cross-tab and
#   nesting reports follow the first three.
#
# - Koo & Li (2016) provide the now-standard ICC interpretation cutoffs
#   used below for rater-mediated assessment:
#     < 0.5 Poor, 0.5-0.75 Moderate, 0.75-0.9 Good, > 0.9 Excellent.
#
# - The design-effect formula is Kish (1965): Deff = 1 + (m - 1) * rho,
#   where m is the average cluster size (ratings per facet element) and
#   rho is the intra-class correlation.
#
# - FACETS itself does not surface ICC or a Kish design effect; it
#   reports rater separation/reliability on the Rasch metric. Because
#   FACETS advises that rater reliability "near 0.0 is preferred"
#   (Linacre's winsteps help, `reliability.htm`), mfrmr reports ICC
#   here as a complementary descriptive summary rather than a
#   replacement for facet separation/reliability.
# ==============================================================================


# ---- internal helpers -----------------------------------------------------

#' @keywords internal
#' @noRd
.ha_entropy <- function(x, weights = rep(1, length(x))) {
  # Shannon entropy H(X) in nats using observed frequencies. NA-safe:
  # NA levels are dropped before counting.
  keep <- !is.na(x) & weights > 0
  x <- x[keep]
  weights <- weights[keep]
  if (length(x) == 0L) return(0)
  weights <- weights / max(weights)
  freq <- as.numeric(rowsum(weights, match(x, unique(x)), reorder = FALSE))
  p <- freq / sum(freq)
  p <- p[p > 0]
  -sum(p * log(p))
}

#' @keywords internal
#' @noRd
.ha_conditional_entropy <- function(y, x, weights = rep(1, length(x))) {
  # H(Y | X) = sum_x p(x) * H(Y | X = x)
  keep <- !is.na(x) & !is.na(y) & weights > 0
  y <- y[keep]; x <- x[keep]
  weights <- weights[keep]
  if (length(x) == 0L) return(0)
  weights <- weights / max(weights)
  groups <- split(seq_along(y), match(x, unique(x)))
  sum(vapply(groups, function(i) sum(weights[i]) * .ha_entropy(y[i], weights[i]), numeric(1))) /
    sum(weights)
}

#' @keywords internal
#' @noRd
.ha_nesting_classify <- function(idx) {
  if (!is.finite(idx)) return(NA_character_)
  if (idx >= 0.99) "Fully nested"
  else if (idx >= 0.95) "Near-perfectly nested"
  else if (idx >= 0.50) "Partially nested"
  else "Crossed"
}

#' @keywords internal
#' @noRd
.ha_sample_classify <- function(n, thresholds) {
  if (!is.finite(n)) return(NA_character_)
  sparse <- as.numeric(thresholds[["sparse"]] %||% 10)
  marginal <- as.numeric(thresholds[["marginal"]] %||% 30)
  standard <- as.numeric(thresholds[["standard"]] %||% 50)
  if (n < sparse) "sparse"
  else if (n < marginal) "marginal"
  else if (n < standard) "standard"
  else "strong"
}


# ---- 1. detect_facet_nesting ----------------------------------------------

#' Detect nesting structure between facets
#'
#' Classifies both directions of each pair of facets (optionally including `Person`)
#' as crossed, partially nested, near-perfectly nested, or fully nested,
#' based on a conditional-entropy index:
#' \deqn{\text{nesting\_index}(A \to B) = 1 - H(B \mid A) / H(B).}
#' An index near 1 means that knowing the level of `A` essentially
#' determines the level of `B` (A is nested in B).
#'
#' This is a pure descriptive review of the observed design. It does not
#' affect estimation; fit_mfrm() continues to treat all facets as fixed
#' effects.
#'
#' @param data Data frame in long format (one row per rating).
#' @param facets Character vector of facet column names.
#' @param person Optional name of the person column (adds Person to the
#'   nesting matrix if supplied).
#' @param weight_col Optional name of a numeric column of finite, nonnegative
#'   weights. Zero-weight rows are excluded; at least one weight must be positive.
#'   Weighted frequencies are used directly, without rounding or replicating rows.
#'
#' @section Classification bands:
#' - `"Fully nested"`: nesting index >= 0.99.
#' - `"Near-perfectly nested"`: 0.95 <= index < 0.99.
#' - `"Partially nested"`: 0.50 <= index < 0.95.
#' - `"Crossed"`: index < 0.50.
#'
#' The direction column records which facet is nested in which, or
#' `"crossed"` when neither direction is above 0.95.
#'
#' @section Missing values and weights:
#' Each pair uses rows observed on both columns and having positive weight.
#' Both marginal entropies and both conditional entropies use those same rows.
#' `UsedRows` and `ExcludedRows` report this pair-specific scope; they count
#' original rows, not the sum of weights. Different pairs may use different rows.
#' Multiplying all weights by a positive constant leaves the indices unchanged.
#' Weights are descriptive frequency weights, not a survey-variance correction.
#' The index is undefined (`NA`) when its denominator entropy is zero, including
#' pairs with no usable rows. Missing-value exclusion does not establish that
#' the observed pattern represents the missing ratings.
#'
#' @section Interpreting output:
#' A `Direction` value of `"Rater nested in Region"` means that every
#' rater appears in exactly one region (or very close to it). For
#' additive fixed-effects MFRM, this is a concern: the severity of a
#' rater is confounded with region-level variance that the model cannot
#' partition. Consider reporting the nesting direction explicitly and,
#' when relevant, refitting without the nested facet or moving to a
#' hierarchical estimation tool (e.g. `lme4::lmer`, `brms`, `TAM`) to
#' separate the variance components.
#'
#' `Direction = "crossed"` is the most common reading when both nesting
#' indices are below 0.5; the two facets largely co-occur at multiple
#' combinations, which is the setting Linacre (1989) assumed.
#'
#' @section Typical workflow:
#' 1. Call `detect_facet_nesting(data, facets)` before fitting.
#' 2. If any pair is flagged as nested or partially nested, review the
#'    numeric index and the `LevelsA`/`LevelsB` counts.
#' 3. For downstream reporting, use [analyze_hierarchical_structure()]
#'    to bundle this output with ICC and design-effect summaries, which
#'    [reporting_checklist()] can record as available report material.
#'
#' @return A list of class `mfrm_facet_nesting` with:
#' - `pairwise_table`: one row per facet pair, reporting both directions with
#'   `NestingIndex_AinB`, `NestingIndex_BinA`, classification strings,
#'   `Direction`, `UsedRows` and `ExcludedRows`.
#' - `summary`: a one-line summary table with facet counts and whether
#'   any non-crossed structure was detected.
#' - `facets`: the facet vector that was reviewed.
#'
#' @seealso [facet_small_sample_review()],
#'   [analyze_hierarchical_structure()], [compute_facet_icc()],
#'   [compute_facet_design_effect()], [fit_mfrm()] (see "Fixed effects
#'   assumption" in its details).
#'
#' @references
#' Shannon, C. E. (1948). A mathematical theory of communication.
#' *Bell System Technical Journal, 27*, 379-423, 623-656.
#' <https://web.mit.edu/6.976/www/handout/shannon.pdf>
#'
#' McEwen, M. R. (2018). *The effects of incomplete rating designs on
#' results from many-facets-Rasch model analyses* (Doctoral thesis,
#' Brigham Young University). <https://scholarsarchive.byu.edu/etd/6689/>
#'
#' Linacre, J. M. (1989). *Many-facet Rasch measurement*. MESA Press.
#'
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' nesting <- detect_facet_nesting(toy, c("Rater", "Criterion"))
#' summary(nesting)
#'
#' # Synthetic example: raters fully nested within regions.
#' d <- data.frame(
#'   Person = rep(paste0("P", formatC(1:20, width = 2, flag = "0")),
#'                each = 6),
#'   Rater  = rep(paste0("R", 1:6), 20),
#'   Region = rep(rep(c("A", "A", "B", "B", "C", "C"), 20)),
#'   Score  = sample(0:4, 120, replace = TRUE),
#'   stringsAsFactors = FALSE
#' )
#' nest <- detect_facet_nesting(d, c("Rater", "Region"))
#' nest$pairwise_table[, c("FacetA", "FacetB",
#'                         "NestingIndex_AinB", "Direction")]
#' @export
detect_facet_nesting <- function(data, facets, person = NULL,
                                 weight_col = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  facets <- as.character(facets)
  missing_cols <- setdiff(c(facets, person), names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }

  all_cols <- unique(c(person, facets))
  w <- rep(1, nrow(data))
  if (!is.null(weight_col)) {
    if (!is.character(weight_col) || length(weight_col) != 1L ||
        is.na(weight_col) || !weight_col %in% names(data)) {
      stop("`weight_col` must name a column in `data`.", call. = FALSE)
    }
    w <- data[[weight_col]]
    if (!is.numeric(w) || is.complex(w) || !is.null(dim(w)) ||
        any(!is.finite(w)) || any(w < 0) || !any(w > 0)) {
      stop("Weights must be finite, nonnegative numeric values with at least one positive value.", call. = FALSE)
    }
  }
  n_cols <- length(all_cols)
  if (n_cols < 2L) {
    return(structure(
      list(
        pairwise_table = data.frame(),
        summary = data.frame(
          NFacets = n_cols,
          AnyNested = FALSE,
          Note = "At least two facets required for a nesting review.",
          stringsAsFactors = FALSE
        ),
        facets = all_cols
      ),
      class = "mfrm_facet_nesting"
    ))
  }

  pairs <- list()
  idx <- 1L
  for (i in seq_len(n_cols - 1L)) {
    for (j in seq(i + 1L, n_cols)) {
      a <- all_cols[i]; b <- all_cols[j]
      keep <- !is.na(data[[a]]) & !is.na(data[[b]]) & w > 0
      av <- data[[a]][keep]; bv <- data[[b]][keep]; weights <- w[keep]
      ha <- .ha_entropy(av, weights)
      hb <- .ha_entropy(bv, weights)
      hb_given_a <- .ha_conditional_entropy(bv, av, weights)
      ha_given_b <- .ha_conditional_entropy(av, bv, weights)
      # Clip only arithmetic roundoff: conditional entropy cannot exceed its marginal.
      nesting_a_in_b <- if (hb > 0) max(0, min(1, 1 - hb_given_a / hb)) else NA_real_
      nesting_b_in_a <- if (ha > 0) max(0, min(1, 1 - ha_given_b / ha)) else NA_real_
      cls_a_in_b <- .ha_nesting_classify(nesting_a_in_b)
      cls_b_in_a <- .ha_nesting_classify(nesting_b_in_a)
      direction <- if (is.na(nesting_a_in_b) || is.na(nesting_b_in_a)) {
        "insufficient_data"
      } else if (nesting_a_in_b >= 0.95 && nesting_b_in_a >= 0.95) {
        "isomorphic"
      } else if (nesting_a_in_b >= 0.95) {
        paste0(a, " nested in ", b)
      } else if (nesting_b_in_a >= 0.95) {
        paste0(b, " nested in ", a)
      } else {
        "crossed"
      }
      pairs[[idx]] <- data.frame(
        FacetA = a,
        FacetB = b,
        LevelsA = dplyr::n_distinct(av),
        LevelsB = dplyr::n_distinct(bv),
        UsedRows = sum(keep),
        ExcludedRows = sum(!keep),
        NestingIndex_AinB = round(nesting_a_in_b, 4),
        NestingIndex_BinA = round(nesting_b_in_a, 4),
        ClassificationAinB = cls_a_in_b,
        ClassificationBinA = cls_b_in_a,
        Direction = direction,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  pairwise <- do.call(rbind, pairs)
  any_nested <- any(
    pairwise$ClassificationAinB %in% c("Fully nested", "Near-perfectly nested") |
    pairwise$ClassificationBinA %in% c("Fully nested", "Near-perfectly nested")
  )
  summary_tbl <- data.frame(
    NFacets = n_cols,
    NPairs = nrow(pairwise),
    AnyNested = any_nested,
    FullyNestedPairs = sum(
      pairwise$ClassificationAinB == "Fully nested" |
        pairwise$ClassificationBinA == "Fully nested",
      na.rm = TRUE
    ),
    CrossedPairs = sum(pairwise$Direction == "crossed", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      pairwise_table = pairwise,
      summary = summary_tbl,
      facets = all_cols
    ),
    class = "mfrm_facet_nesting"
  )
}


# ---- 2. facet_small_sample_review ------------------------------------------

#' Review rating counts for each facet level
#'
#' Reports per-level observation counts, SE, and fit statistics for every
#' level of every facet in a fitted MFRM model, and classifies each level
#' as `"sparse"`, `"marginal"`, `"standard"`, or `"strong"` using
#' configurable count bands. These labels identify levels to inspect;
#' they do not establish measurement precision or sample adequacy.
#'
#' `N` counts retained rating rows for a level in `fit$prep$data`, without
#' summing observation weights. For a rater it counts ratings, not distinct
#' examinees; for a person it counts that person's ratings. Repeated or
#' dependent ratings need not supply independent information. Review the
#' counts alongside targeting, connectedness, category use, the estimation
#' method, and available uncertainty estimates before drawing conclusions.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output from the same `fit`.
#'   When supplied, per-level `Infit` and `Outfit` are added. `Estimate` and
#'   `SE` come from the fitted facet tables; unavailable SEs remain missing.
#' @param thresholds Named numeric vector of count bands. Defaults are
#'   `c(sparse = 10, marginal = 30, standard = 50)`. Supply exactly these
#'   three names with finite, positive values satisfying
#'   `sparse < marginal < standard`; their order in the vector is immaterial.
#'   The defaults are package screening choices, not precision thresholds.
#'
#' @section Interpreting output:
#' With the default thresholds, the four labels mean only `N < 10`,
#' `10 <= N < 30`, `30 <= N < 50`, and `N >= 50`, respectively.
#' `WorstCategory` is the lowest count band represented within a facet.
#' Neither `"standard"` nor `"strong"` guarantees narrow intervals, correct
#' coverage, good targeting, or suitability for publication. Do not pool or
#' discard levels solely because of a count label.
#'
#' Linacre (1994) relates Rasch item/person precision to sample size under
#' targeting and model assumptions. Those guidelines do not establish a
#' precision guarantee for these raw many-facet rating counts. Jones and
#' Wind (2018) study repeated ratings in incomplete designs and highlight
#' examinee precision. Integrating over a population distribution in MML
#' does not remove uncertainty in individual person scores. Interpret `SE`
#' according to the fitted estimator; this count review adds no interval
#' qualification or hierarchical model.
#'
#' @section Typical workflow:
#' 1. Fit with `fit_mfrm()`; optionally also produce `diagnostics`
#'    with `diagnose_mfrm()` if you want per-level Infit/Outfit.
#' 2. Call `facet_small_sample_review(fit, diagnostics)`.
#' 3. Read the `facet_summary` first: it highlights the worst level
#'    per facet. The `summary` table gives counts in each band.
#' 4. If any facet is flagged as sparse or marginal, discuss it in the
#'    Methods section; [build_apa_outputs()] already adds a sentence
#'    about the band when `fit$summary$FacetSampleSizeFlag` is set.
#'
#' @return A list of class `mfrm_facet_sample_review` with:
#' - `table`: one row per `(Facet, Level)` with `N`, `Estimate`, `SE`,
#'   `Infit`, `Outfit`, and `SampleCategory`.
#' - `summary`: counts of levels in each sample-size category, by facet.
#' - `facet_summary`: count range and lowest count band per facet.
#' - `thresholds`: the applied count bands.
#' - `notes`: count definition and interpretation limits, also shown by
#'   `print()` and `summary()`.
#'
#' @seealso [detect_facet_nesting()], [analyze_hierarchical_structure()],
#'   [compute_facet_icc()], [compute_facet_design_effect()],
#'   [reporting_checklist()].
#'
#' @references
#' Linacre, J. M. (1994). Sample size and item calibration stability.
#' *Rasch Measurement Transactions, 7*(4), 328.
#' \url{https://rasch.org/rmt/rmt74m.htm}
#'
#' Jones, E., & Wind, S. A. (2018). Using repeated ratings to improve
#' measurement precision in incomplete rating designs. *Journal of
#' Applied Measurement, 19*(2), 148-161.
#' PMID: 29894984.
#'
#' @examples
#' toy <- load_mfrmr_data("example_operational")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", quad_points = 7, maxit = 30)
#' review <- facet_small_sample_review(fit)
#' summary(review)
#'
#' # Custom thresholds (e.g. a stricter protocol).
#' strict <- facet_small_sample_review(
#'   fit,
#'   thresholds = c(sparse = 15, marginal = 40, standard = 100)
#' )
#' strict$facet_summary
#' @name facet_small_sample_review
#' @export
#' @order 1
facet_small_sample_review <- function(fit, diagnostics = NULL,
                                      thresholds = c(sparse = 10,
                                                     marginal = 30,
                                                     standard = 50)) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (!is.null(diagnostics)) {
    mfrm_results_validate_diagnostics_identity(fit, diagnostics, helper = "facet_small_sample_review()")
  }
  bands <- c("sparse", "marginal", "standard")
  if (!is.numeric(thresholds) || is.complex(thresholds) ||
      !is.null(dim(thresholds)) || length(thresholds) != 3L ||
      !setequal(names(thresholds), bands) || any(!is.finite(thresholds)) ||
      any(thresholds <= 0) || any(diff(thresholds[bands]) <= 0)) {
    stop("`thresholds` must contain exactly the finite, positive numeric values ",
         "`sparse < marginal < standard`, named accordingly.", call. = FALSE)
  }
  thresholds <- thresholds[bands]

  prep <- fit$prep
  data <- prep$data
  facets <- c("Person", prep$facet_names)
  rows <- list()

  # Build per-level count and estimate tables.
  # fit$facets$person has columns Person / Estimate (no Level/Facet); other
  # facets live in fit$facets$others with Facet / Level / Estimate.
  other_tbl <- if (!is.null(fit$facets$others)) fit$facets$others else NULL

  extract_estimates <- function(facet) {
    if (facet == "Person") {
      p <- fit$facets$person
      if (is.null(p) || nrow(p) == 0L) return(NULL)
      id_col <- if ("Level" %in% names(p)) "Level" else if ("Person" %in% names(p)) "Person" else names(p)[1]
      data.frame(
        Level = as.character(p[[id_col]]),
        Estimate = suppressWarnings(as.numeric(p$Estimate %||% NA_real_)),
        SE = suppressWarnings(as.numeric(p$SE %||% p$ModelSE %||% NA_real_)),
        stringsAsFactors = FALSE
      )
    } else if (!is.null(other_tbl) && nrow(other_tbl) > 0L) {
      tbl <- other_tbl[as.character(other_tbl$Facet) == facet, , drop = FALSE]
      if (nrow(tbl) == 0L) return(NULL)
      data.frame(
        Level = as.character(tbl$Level),
        Estimate = suppressWarnings(as.numeric(tbl$Estimate %||% NA_real_)),
        SE = suppressWarnings(as.numeric(tbl$SE %||% tbl$ModelSE %||% NA_real_)),
        stringsAsFactors = FALSE
      )
    } else NULL
  }

  for (facet in facets) {
    obs_counts <- data |>
      dplyr::count(Level = as.character(.data[[facet]]), name = "N")
    est_tbl <- extract_estimates(facet)
    merged <- if (!is.null(est_tbl)) {
      dplyr::left_join(obs_counts, est_tbl, by = "Level")
    } else {
      dplyr::mutate(obs_counts, Estimate = NA_real_, SE = NA_real_)
    }

    if (!is.null(diagnostics) && !is.null(diagnostics$measures) &&
        "Facet" %in% names(diagnostics$measures)) {
      diag_rows <- diagnostics$measures[
        as.character(diagnostics$measures$Facet) == facet, , drop = FALSE]
      if (nrow(diag_rows) > 0L) {
        diag_tbl <- data.frame(
          Level = as.character(diag_rows$Level),
          Infit = suppressWarnings(as.numeric(diag_rows$Infit %||% NA_real_)),
          Outfit = suppressWarnings(as.numeric(diag_rows$Outfit %||% NA_real_)),
          stringsAsFactors = FALSE
        )
        merged <- dplyr::left_join(merged, diag_tbl, by = "Level")
      } else {
        merged <- dplyr::mutate(merged, Infit = NA_real_, Outfit = NA_real_)
      }
    } else {
      merged <- dplyr::mutate(merged, Infit = NA_real_, Outfit = NA_real_)
    }

    merged$Facet <- facet
    merged$SampleCategory <- vapply(merged$N, .ha_sample_classify,
                                    character(1), thresholds = thresholds)
    rows[[facet]] <- merged
  }
  audit_tbl <- dplyr::bind_rows(rows) |>
    dplyr::select(Facet, Level, N, Estimate, SE, Infit, Outfit,
                  SampleCategory)

  category_order <- c("sparse", "marginal", "standard", "strong")
  level_count_by_facet <- audit_tbl |>
    dplyr::count(Facet, SampleCategory) |>
    tidyr::pivot_wider(
      names_from = SampleCategory, values_from = n, values_fill = 0L
    )
  for (nm in category_order) {
    if (!nm %in% names(level_count_by_facet)) {
      level_count_by_facet[[nm]] <- 0L
    }
  }
  level_count_by_facet <- level_count_by_facet[, c("Facet", category_order)]

  facet_summary <- audit_tbl |>
    dplyr::group_by(Facet) |>
    dplyr::summarize(
      Levels = dplyr::n(),
      MinN = min(N, na.rm = TRUE),
      MedianN = stats::median(N, na.rm = TRUE),
      MaxN = max(N, na.rm = TRUE),
      WorstCategory = category_order[
        min(match(SampleCategory, category_order), na.rm = TRUE)
      ],
      .groups = "drop"
    )

  structure(
    list(
      table = as.data.frame(audit_tbl, stringsAsFactors = FALSE),
      summary = as.data.frame(level_count_by_facet, stringsAsFactors = FALSE),
      facet_summary = as.data.frame(facet_summary, stringsAsFactors = FALSE),
      thresholds = as.list(thresholds),
      notes = c(
        "N counts retained rating rows per level, not distinct persons, weight sums, or effective sample sizes.",
        "Count bands are screening labels; they do not establish precision, interval coverage, or targeting."
      )
    ),
    class = "mfrm_facet_sample_review"
  )
}

# ---- 3. compute_facet_icc / design effect ---------------------------------

#' Compute intra-class correlations for each facet
#'
#' Fits a random-effects variance-components model
#' `Score ~ 1 + (1 | Person) + (1 | Facet1) + (1 | Facet2) + ...`
#' using `lme4::lmer` (in `Suggests`) and returns the proportion of
#' observed score variance attributable to each facet. This is a
#' descriptive summary complementary to the Rasch-metric rater
#' separation/reliability reported elsewhere.
#'
#' @param data Data frame in long format.
#' @param facets Character vector of facet column names.
#' @param score Name of the score column.
#' @param person Optional person column. If supplied it is added as a
#'   separate random intercept so Person-level variance is partitioned
#'   out.
#' @param reml Logical; whether to fit with REML. Default `TRUE`.
#' @param ci_method Confidence-interval method for the ICC column.
#'   One of `"none"` (default, point estimate only) or `"boot"`
#'   (parametric percentile bootstrap via [lme4::bootMer()]). Each
#'   simulated data set is refitted and its full variance decomposition
#'   used to calculate the ICC ratios. The former `"profile"` method
#'   is no longer supported; see Updating saved intervals below.
#' @param ci_level Confidence level when `ci_method != "none"`; default
#'   `0.95`. Koo & Li (2016) recommend banding the CI rather than the
#'   point estimate when classifying reliability as Poor / Moderate /
#'   Good / Excellent.
#' @param ci_boot_reps Number of bootstrap replicates used when
#'   `ci_method = "boot"`. An integer of at least 2; default `1000`.
#'   Small counts give imprecise tail quantiles; choose enough replicates
#'   for the precision required in the application.
#' @param ci_boot_seed Optional integer seed for the bootstrap path
#'   (between 0 and `.Machine$integer.max`). `NULL` uses the current
#'   random-number state. Bootstrap simulation advances that state.
#' @param ci_boot_parallel Parallelisation strategy for the
#'   parametric-bootstrap CI path, passed through to
#'   [lme4::bootMer()]: `"no"` (default), `"multicore"` (POSIX
#'   `mclapply`), or `"snow"` (PSOCK cluster). `"multicore"` does
#'   nothing on Windows and falls back to serial; use `"snow"` there.
#' @param ci_boot_ncpus Number of CPUs to use for the parallel
#'   bootstrap path (ignored when `ci_boot_parallel = "no"`). A positive
#'   integer. Interactive progress is available for serial execution.
#' @param missing How to handle missing scores or selected grouping values:
#'   `"error"` (default) or explicit complete-case omission with `"omit"`.
#'   Numeric character/factor score labels retain their numeric values.
#'   Nonnumeric scores, infinite values, and blank grouping labels are refused.
#'
#' @section Rows used:
#' Missingness is checked only in the score, facets, and optional person column.
#' `InputRows`, `UsedRows`, and `ExcludedRows` are included in the table.
#' `attr(x, "data_usage")` retains these counts, excluded input row positions,
#' missing columns per row, and observed grouping-level counts after omission.
#' [compute_facet_design_effect()] uses these retained counts for its sample
#' sizes. Omission does not impute scores or correct missing-data bias.
#'
#' @section Score units and zero variation:
#' A small positive variance is not treated as zero using a fixed cutoff.
#' Multiplying scores by a nonzero constant leaves the variance shares
#' unchanged, up to fitting precision, while variances change by its square.
#' Variance estimates are retained without decimal rounding. If all retained
#' scores are equal, or the fitted total variance is not positive and finite,
#' variances and ICCs are unavailable (`NA`); numerical fitting residue is not
#' interpreted as observed variation. A constant-response bootstrap refit is
#' also unavailable and withholds the interval.
#'
#' @section Interpreting output:
#' The `Interpretation` column uses **two scales** so the same numeric
#' ICC reads correctly for each facet role:
#'
#' - For the `person` facet, higher ICC = better. Koo & Li (2016, p. 161)
#'   bands are applied: `< 0.5` Poor, `[0.5, 0.75]` Moderate,
#'   `(0.75, 0.9]` Good, `> 0.9` Excellent. The strict `>` boundary at
#'   0.9 follows Koo & Li's wording "values greater than 0.90 indicate
#'   excellent reliability" (so an ICC of exactly 0.9 reads as Good).
#' - For non-person facets (Rater, Criterion, Task, Region, ...) the
#'   same numeric value is a **variance share**: how much of the total
#'   observed score variance sits at that facet. The bands used here
#'   are different (`Trivial share` < 0.05, `Small share` < 0.15,
#'   `Moderate share` < 0.30, `Large share` >= 0.30), and a large
#'   rater share is generally *bad* news (raters disagree about
#'   averages), not good news.
#'
#' The `InterpretationScale` column explicitly records which scale
#' applies to each row, so downstream reporting does not confuse the
#' two. FACETS (Linacre, 2026) reports rater separation/reliability on
#' the Rasch metric instead of an ICC; mfrmr surfaces both, with the
#' Rasch-metric version in `diagnostics$reliability` and this
#' variance-share view here.
#' For agreement or consistency between specific pairs on common units, use
#' [pairwise_icc_table()] instead. A small rater variance share here does not
#' establish pairwise agreement; residual disagreement can still be large.
#'
#' Set `ci_method = "boot"` to request intervals alongside the point
#' estimates. The `Interpretation` column still uses point estimates.
#' The bootstrap simulates Gaussian random effects and errors from the
#' fitted model; it does not correct model misspecification or missing-data
#' bias. Percentile coverage can be unreliable near zero variance components
#' or with few grouping levels.
#'
#' Intervals are withheld if the original fit has convergence problems or
#' warnings, or if any requested bootstrap refit fails, has convergence
#' problems or warnings, or produces an undefined ICC. Finite draws are
#' retained but never silently selected to calculate an interval. Singular
#' fits (zero random-effect components) are recorded separately and retained
#' when they converge; they are not automatically treated as failures.
#' Inspect `ICC_CI_Status` and `attr(x, "icc_ci")` before reporting intervals.
#'
#' @section Updating saved intervals:
#' The former `ci_method = "profile"` transformed separate standard-deviation
#' intervals while holding other variance components fixed. These are not
#' profile-likelihood intervals for the ICC ratio and should not be reported
#' as ICC confidence intervals. Requests now stop with an explanation.
#' Rerun [compute_facet_icc()] or [analyze_hierarchical_structure()] with the
#' original data and settings, choosing `ci_method = "boot"` explicitly if
#' intervals are needed. Saved bootstrap results from earlier versions must
#' also be rerun to obtain complete failure accounting. Printing, summarizing,
#' or plotting old interval results cannot correct their calculations.
#'
#' @section Typical workflow:
#' 1. Fit the MFRM model with `fit_mfrm()` for the Rasch-metric
#'    separation/reliability.
#' 2. Call `compute_facet_icc(data, facets, score, person)` to get the
#'    complementary variance-share summary.
#' 3. Feed into [compute_facet_design_effect()] to convert ICCs and
#'    average cluster sizes into descriptive, per-facet design-effect
#'    approximations. These do not estimate the precision of the full design.
#'
#' @return A data.frame of class `mfrm_facet_icc` with one row per
#'   variance component (including a `"Residual"` row) and columns:
#' - `Facet`: the grouping factor name (or `"Residual"`).
#' - `Variance`: unrounded variance estimate (REML by default, ML if
#'   `reml = FALSE`); `NA` when the variance shares are undefined.
#' - `ICC`: variance share (`Variance / sum(Variance)`), in `[0, 1]`.
#' - `Interpretation`: band label according to the facet's scale.
#' - `InterpretationScale`: `"Koo-Li reliability"` for the person
#'   facet, `"Variance share"` for others.
#' - `ICC_CI_Lower` / `ICC_CI_Upper` / `ICC_CI_Level` / `ICC_CI_Method`:
#'   CI bounds (unavailable bounds are `NA`), requested level, and method.
#' - `ICC_CI_Status`: whether intervals are available and, otherwise, why.
#' - `ICC_CI_NRequested` / `ICC_CI_NReps` / `ICC_CI_NUnavailable`: requested
#'   bootstrap count, number of converged refits with all ICCs finite, and
#'   number without such a result (absent for `"none"`). Counts are `NA` when
#'   the bootstrap aborts without returning its draws. Warnings can withhold
#'   intervals even when all draws are finite and all refits converge.
#'
#' The `icc_ci` attribute retains the original fit's convergence and singularity
#' diagnostics and warnings. When bootstrap results are returned, its `bootstrap`
#' entry contains every draw, per-refit convergence and singularity indicators,
#' the number of refit errors, and lme4's message/warning/error tables.
#'
#' @seealso [compute_facet_design_effect()], [pairwise_icc_table()],
#'   [analyze_hierarchical_structure()], [detect_facet_nesting()],
#'   [facet_small_sample_review()].
#'
#' @concept confidence intervals
#' @concept hierarchical structure
#' @concept ICC
#'
#' @references
#' Koo, T. K., & Li, M. Y. (2016). A guideline of selecting and
#' reporting intraclass correlation coefficients for reliability
#' research. *Journal of Chiropractic Medicine, 15*(2), 155-163.
#'
#' Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting
#' linear mixed-effects models using lme4. *Journal of Statistical
#' Software, 67*(1), 1-48.
#'
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("compute_facet_icc", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   icc <- compute_facet_icc(toy, facets = c("Rater", "Criterion"),
#'                            score = "Score", person = "Person")
#'   print(icc)
#'   # Look for:
#'   # - Person ICC reads as Koo & Li (2016) reliability: < 0.5 poor,
#'   #   0.5-0.75 moderate, 0.75-0.9 good, > 0.9 excellent.
#'   # - Rater / Criterion ICC reads as variance share, NOT reliability;
#'   #   small values indicate limited systematic facet differences, and
#'   #   shares > 0.10 hint at meaningful systematic facet differences.
#'   # - `Interpretation` summarises the variance-share band the helper
#'   #   has assigned to each row.
#' }
#' }
#' @export
compute_facet_icc <- function(data, facets, score,
                              person = NULL, reml = TRUE,
                              ci_method = c("none", "boot"),
                              ci_level = 0.95,
                              ci_boot_reps = 1000L,
                              ci_boot_seed = NULL,
                              ci_boot_parallel = c("no", "multicore", "snow"),
                              ci_boot_ncpus = 1L,
                              missing = c("error", "omit")) {
  missing <- match.arg(missing)
  ci_method <- .icc_ci_method(ci_method)
  ci_boot_parallel <- match.arg(ci_boot_parallel)
  if (!is.numeric(ci_level) || is.complex(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (ci_method == "boot") {
    for (arg in c("ci_boot_reps", "ci_boot_ncpus", "ci_boot_seed")) {
      value <- get(arg)
      if (arg == "ci_boot_seed" && is.null(value)) next
      minimum <- switch(arg, ci_boot_reps = 2, ci_boot_ncpus = 1, 0)
      if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
          !is.finite(value) || value < minimum ||
          value > .Machine$integer.max || value != floor(value)) {
        stop("`", arg, "` must be a single integer between ", minimum,
             " and .Machine$integer.max.", call. = FALSE)
      }
    }
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    message("`compute_facet_icc()` requires the `lme4` package ",
            "(in Suggests). Install it and retry.")
    return(structure(
      data.frame(
        Facet = character(0), Variance = numeric(0),
        ICC = numeric(0), Interpretation = character(0),
        stringsAsFactors = FALSE
      ),
      class = c("mfrm_facet_icc", "data.frame")
    ))
  }
  if (!is.data.frame(data) || nrow(data) == 0L || anyNA(names(data)) ||
      anyDuplicated(names(data)) || any(!nzchar(names(data)))) {
    stop("`data` must be a nonempty data frame with unique, nonmissing column names.", call. = FALSE)
  }
  if (!is.character(score) || length(score) != 1L || is.na(score) || !nzchar(score) ||
      !is.character(facets) || !length(facets) || anyNA(facets) ||
      anyDuplicated(facets) || any(!nzchar(facets)) ||
      (!is.null(person) && (!is.character(person) || length(person) != 1L ||
                           is.na(person) || !nzchar(person)))) {
    stop("Supply a score column name, distinct facet names, and an optional person column name.", call. = FALSE)
  }
  re_terms <- unique(c(person, facets))
  if (any(re_terms %in% c(score, "Residual"))) {
    stop("Grouping columns cannot be the score column or named `Residual`.", call. = FALSE)
  }
  if (!is.logical(reml) || length(reml) != 1L || is.na(reml)) {
    stop("`reml` must be TRUE or FALSE.", call. = FALSE)
  }
  needed_cols <- c(score, re_terms)
  missing_cols <- setdiff(needed_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  data <- as.data.frame(data)
  for (col in re_terms) {
    value <- data[[col]]
    if (!is.null(dim(value)) || !(is.factor(value) || (!is.object(value) &&
        (is.character(value) || is.logical(value) || (is.numeric(value) && !is.complex(value))))) ||
        (is.numeric(value) && any(is.infinite(value)))) {
      stop("Grouping column '", col, "' must contain finite numeric, character, factor, or logical labels, or NA.", call. = FALSE)
    }
    if (any(!is.na(value) & !nzchar(trimws(as.character(value))))) {
      stop("Grouping column '", col, "' contains blank labels; use NA for missing values.", call. = FALSE)
    }
    if (!is.factor(value)) {
      labels <- as.character(value)
      labels[is.na(value)] <- NA_character_
      data[[col]] <- as.factor(labels)
    }
  }
  value <- data[[score]]
  if (!is.null(dim(value)) || !(is.factor(value) || (!is.object(value) &&
      (is.character(value) || (is.numeric(value) && !is.complex(value)))))) {
    stop("The score column must be numeric or numeric character/factor labels, with NA for missing values.", call. = FALSE)
  }
  data[[score]] <- suppressWarnings(as.numeric(if (is.factor(value)) as.character(value) else value))
  if (any(!is.na(value) & !is.finite(data[[score]]))) {
    stop("The score column contains nonnumeric or infinite values; correct them explicitly and use NA for missing values.", call. = FALSE)
  }
  absent <- is.na(data[, needed_cols, drop = FALSE])
  keep <- rowSums(absent) == 0L
  cells <- which(absent, arr.ind = TRUE)
  data_usage <- list(source = "Supplied data", missing = missing,
    counts = c(InputRows = nrow(data), UsedRows = sum(keep), ExcludedRows = sum(!keep)),
    excluded_rows = which(!keep),
    missing_cells = data.frame(InputRow = cells[, 1L], Column = needed_cols[cells[, 2L]],
                               row.names = NULL))
  if (any(!keep) && missing == "error") {
    stop(sum(!keep), " row(s) have missing scores or selected grouping values. Review the data or choose missing = 'omit' explicitly.", call. = FALSE)
  }
  if (!any(keep)) stop("No complete rows remain for the ICC model.", call. = FALSE)
  data <- data[keep, , drop = FALSE]
  data_usage$observed_levels <- vapply(data[re_terms], function(x) length(unique(x)), integer(1))

  quoted_names <- vapply(c(score, re_terms), function(name) {
    deparse1(as.name(name), backtick = TRUE)
  }, character(1))
  formula <- stats::as.formula(paste0(
    quoted_names[1L], " ~ 1 + ",
    paste0("(1 | ", quoted_names[-1L], ")", collapse = " + ")
  ))
  # Retain fit warnings separately from convergence codes and boundary messages.
  lmer_warnings <- character(0)
  fit <- tryCatch(
    withCallingHandlers(
      .mfrm_icc_lmer(formula, data = data, REML = reml, na.action = stats::na.fail),
      warning = function(w) {
        lmer_warnings <<- c(lmer_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    stop("lme4::lmer failed: ", conditionMessage(fit), call. = FALSE)
  }
  if (length(lmer_warnings) > 0L) {
    message("compute_facet_icc(): lme4 reported ",
            length(lmer_warnings), " fit warning(s); ",
            "the ICC table is returned but results may be unreliable. ",
            "First message: ", lmer_warnings[1])
  }

  vc <- as.data.frame(lme4::VarCorr(fit))
  vc <- vc[is.na(vc$var2), c("grp", "vcov")]
  total_var <- sum(vc$vcov)
  # Constant responses can leave positive fitting residue. Detect them on the
  # observed scale without rejecting legitimate small-unit variance estimates.
  undefined_icc <- all(data[[score]] == data[[score]][1L]) ||
    !is.finite(total_var) || total_var <= 0
  icc_vec <- if (!undefined_icc) {
    vc$vcov / total_var
  } else {
    rep(NA_real_, length(vc$vcov))
  }
  if (undefined_icc) {
    message("compute_facet_icc(): ICCs are undefined: scores are constant or ",
            "the fitted total variance is not positive and finite. ",
            "Returning NA variances and ICCs.")
  }

  # Variance-share labels. For a _person_ facet in rater-mediated data, this
  # corresponds to the Koo & Li (2016) reliability interpretation
  # (higher = better). For rater, criterion, or similar non-person facets
  # the ICC is the variance share, not reliability, and the convention is
  # the opposite direction: small shares are desirable (raters agree).
  person_label <- person %||% "Person"
  var_share_band <- function(i) {
    if (!is.finite(i)) return(NA_character_)
    if (i < 0.05) "Trivial share"
    else if (i < 0.15) "Small share"
    else if (i < 0.30) "Moderate share"
    else "Large share"
  }
  koo_li_band <- function(i) {
    if (!is.finite(i)) return(NA_character_)
    # Koo & Li (2016, p. 161): "values greater than 0.90 indicate excellent
    # reliability." Strict > at 0.9 places ICC = 0.9 in Good, not Excellent.
    if (i < 0.5) "Poor"
    else if (i < 0.75) "Moderate"
    else if (i <= 0.9) "Good"
    else "Excellent"
  }
  interpret <- vapply(seq_along(icc_vec), function(k) {
    if (!is.finite(icc_vec[k])) {
      return("Non-identifiable")
    }
    grp <- as.character(vc$grp[k])
    if (identical(grp, person_label)) koo_li_band(icc_vec[k])
    else var_share_band(icc_vec[k])
  }, character(1))

  out <- data.frame(
    Facet = vc$grp,
    Variance = if (undefined_icc) NA_real_ else vc$vcov,
    ICC = round(icc_vec, 4),
    Interpretation = interpret,
    InterpretationScale = ifelse(
      as.character(vc$grp) == person_label,
      "Koo-Li reliability",
      "Variance share"
    ),
    stringsAsFactors = FALSE
  )

  out$ICC_CI_Lower <- NA_real_
  out$ICC_CI_Upper <- NA_real_
  out$ICC_CI_Level <- ci_level
  out$ICC_CI_Method <- ci_method
  out$ICC_CI_Status <- "Not requested"
  fit_diagnostics <- .icc_fit_diagnostics(fit)
  fit_diagnostics$warnings <- lmer_warnings
  ci_details <- list(calculation_version = 2L, fit = fit_diagnostics)
  if (ci_method == "boot") {
    out$ICC_CI_NRequested <- as.integer(ci_boot_reps)
    out$ICC_CI_NReps <- 0L
    out$ICC_CI_NUnavailable <- as.integer(ci_boot_reps)
    if (!all(is.finite(icc_vec))) {
      out$ICC_CI_Status <- "Undefined ICC"
    } else if (!fit_diagnostics$converged || length(lmer_warnings) > 0L) {
      out$ICC_CI_Status <- "Original fit requires review"
    } else {
      ci_result <- tryCatch(
        .compute_icc_ci(
          fit = fit, vc_grp = as.character(vc$grp), ci_level = ci_level,
          boot_reps = ci_boot_reps, boot_seed = ci_boot_seed,
          boot_parallel = ci_boot_parallel, boot_ncpus = ci_boot_ncpus
        ),
        error = function(e) e
      )
      if (inherits(ci_result, "error")) {
        out$ICC_CI_Status <- "Bootstrap failed"
        out$ICC_CI_NReps <- out$ICC_CI_NUnavailable <- NA_integer_
        ci_details$error <- conditionMessage(ci_result)
        message("compute_facet_icc(): bootstrap failed: ", ci_details$error,
                ". Returning point estimates without intervals.")
      } else {
        out$ICC_CI_Lower <- round(ci_result$lower, 4)
        out$ICC_CI_Upper <- round(ci_result$upper, 4)
        out$ICC_CI_Status <- ci_result$status
        out$ICC_CI_NReps <- ci_result$n_reps
        out$ICC_CI_NUnavailable <- ci_boot_reps - ci_result$n_reps
        ci_details$bootstrap <- ci_result$diagnostics
      }
    }
  }
  for (name in names(data_usage$counts)) out[[name]] <- data_usage$counts[[name]]
  structure(out, class = c("mfrm_facet_icc", "data.frame"),
            icc_ci = ci_details, data_usage = data_usage)
}

.icc_ci_method <- function(method) {
  if (identical(method, "profile")) {
    stop('`ci_method = "profile"` has been withdrawn: transforming separate ',
         'variance-component bounds does not give an ICC confidence interval. ',
         'Use "none", or explicitly choose "boot" for a parametric bootstrap.',
         call. = FALSE)
  }
  match.arg(method, c("none", "boot"))
}

.icc_fit_diagnostics <- function(fit) {
  conv <- fit@optinfo$conv
  codes <- c(conv$opt, conv$lme4$code)
  list(converged = all(is.finite(codes)) && all(codes == 0),
       singular = lme4::isSingular(fit), convergence = conv)
}

# Each row contains the jointly refitted variance shares plus fit indicators.
# Boundary components remain in the distribution; failed refits do not disappear.
.compute_icc_ci <- function(fit, vc_grp, ci_level,
                            boot_reps = 1000L, boot_seed = NULL,
                            boot_parallel = "no", boot_ncpus = 1L) {
  n_grp <- length(vc_grp)
  icc_of <- function(fit_b) {
    vc <- as.data.frame(lme4::VarCorr(fit_b))
    vc <- vc[is.na(vc$var2), c("grp", "vcov")]
    tot <- sum(vc$vcov)
    response <- lme4::getME(fit_b, "y")
    share <- if (all(response == response[1L]) || !is.finite(tot) || tot <= 0) {
      rep(NA_real_, n_grp)
    } else {
      as.numeric(stats::setNames(vc$vcov / tot, as.character(vc$grp))[vc_grp])
    }
    conv <- fit_b@optinfo$conv
    codes <- c(conv$opt, conv$lme4$code)
    c(share, converged = as.numeric(all(is.finite(codes)) && all(codes == 0)),
      singular = as.numeric(lme4::isSingular(fit_b)))
  }
  bootstrap_warnings <- character()
  cl <- NULL
  if (boot_parallel == "snow" && boot_ncpus > 1L) {
    cl <- parallel::makePSOCKcluster(boot_ncpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, loadNamespace("lme4"))
  }
  b <- withCallingHandlers(
    .mfrm_boot_mer(fit, FUN = icc_of, nsim = boot_reps, seed = boot_seed,
                  type = "parametric", use.u = FALSE,
                  parallel = boot_parallel, ncpus = boot_ncpus, cl = cl,
                  .progress = if (interactive() && boot_parallel == "no") "txt" else "none"),
    warning = function(w) {
      bootstrap_warnings <<- c(bootstrap_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  draws <- b$t[, seq_len(n_grp), drop = FALSE]
  colnames(draws) <- vc_grp
  converged <- b$t[, n_grp + 1L] == 1
  singular <- b$t[, n_grp + 2L] == 1
  finite <- apply(is.finite(draws) & draws >= 0 & draws <= 1, 1L, all)
  usable <- finite & !is.na(converged) & converged
  n_reps <- sum(usable)
  all_messages <- attr(b, "boot.all.msgs")
  has_warnings <- length(bootstrap_warnings) > 0L ||
    sum(all_messages$`factory-warning`) > 0L
  n_errors <- attr(b, "bootFail") %||% 0L
  status <- if (n_reps != boot_reps || n_errors > 0L) {
    "Incomplete bootstrap"
  } else if (has_warnings) {
    "Bootstrap warnings require review"
  } else "Available"
  lo <- hi <- rep(NA_real_, n_grp)
  if (status == "Available") {
    bounds <- apply(draws, 2L, stats::quantile,
                    probs = c((1 - ci_level) / 2, (1 + ci_level) / 2),
                    names = FALSE)
    lo <- bounds[1L, ]
    hi <- bounds[2L, ]
  }
  list(lower = unname(lo), upper = unname(hi), n_reps = n_reps, status = status,
       diagnostics = list(draws = draws, converged = converged, singular = singular,
                          usable = usable, n_errors = n_errors,
                          failure_messages = attr(b, "boot.fail.msgs"),
                          messages = all_messages, warnings = bootstrap_warnings,
                          seed = boot_seed, parallel = boot_parallel, ncpus = boot_ncpus))
}

# Saved results cannot acquire corrected numerical intervals merely by printing.
.check_icc_intervals <- function(x) {
  methods <- x$ICC_CI_Method
  if (any(methods == "profile", na.rm = TRUE) ||
      (any(methods == "boot", na.rm = TRUE) &&
       !identical(attr(x, "icc_ci")$calculation_version, 2L))) {
    stop("These saved ICC intervals need to be recomputed. Rerun ",
         "compute_facet_icc() or analyze_hierarchical_structure() with the ",
         'original data and settings, choosing ci_method = "boot" explicitly ',
         'for intervals or "none" for point estimates.', call. = FALSE)
  }
  invisible(x)
}

#' Compute descriptive design-effect approximations for each facet
#'
#' Combines per-facet average cluster size with ICC estimates to return
#' the Kish-style approximation `Deff = 1 + (m - 1) * rho`, where `m`
#' is the average number of observations per facet element and `rho` is
#' the ICC variance share. Each facet is evaluated separately.
#'
#' @param data Data frame in long format, used to fit the ICC model when
#'   `icc_table` is `NULL`. With a supplied ICC table, sample sizes come from
#'   its retained row accounting, not from `data`.
#' @param facets Character vector of facet column names.
#' @param icc_table Output from [compute_facet_icc()] (optional; will be
#'   computed on the fly when `NULL`). Must retain its `data_usage` attribute.
#'   Rerun older saved ICC results from their original data and settings before
#'   calculating design effects; their analysis sample cannot be reconstructed
#'   from the ICC table alone.
#' @param score Score column name; required when `icc_table` is `NULL`.
#' @param person Person column; passed through to compute_facet_icc().
#' @param missing Missing-value policy passed to [compute_facet_icc()] when
#'   `icc_table` is `NULL`. A supplied table retains its original policy.
#'
#' @section Interpreting output:
#' The formula describes the variance inflation of an unweighted mean under
#' a single clustering factor with equal cluster sizes, independent clusters,
#' and common within-cluster correlation. This helper substitutes the average
#' cluster size and one fitted facet variance share. It does not calculate
#' the variance of a specified estimator under the full sampling design.
#'
#' - `Deff = 1` means this approximation adds no inflation for that facet;
#'   it does not establish independence of observations or adequate precision.
#' - `Deff > 1` signals potential clustering influence under this approximation.
#'   `EffectiveN = UsedRows / Deff` is a descriptive equivalent row count,
#'   not a count of independent Persons or an assurance of matching precision.
#' - Unequal cluster sizes, crossed or nested dependencies, sampling weights,
#'   and finite-population corrections are not accounted for. Per-facet values
#'   must not be added or multiplied to obtain an overall design effect.
#' - Reported `ICC` is pulled from `icc_table$ICC` (the variance share);
#'   interpretation is the same as in [compute_facet_icc()].
#'
#' @section Typical workflow:
#' 1. Run [compute_facet_icc()] to get the variance-component shares.
#' 2. Feed the result and the data into
#'    `compute_facet_design_effect(data, facets, icc_table = icc)`.
#' 3. Use these values to flag facets for design review. For standard errors,
#'    sample-size planning, or comparisons of precision, use an estimator
#'    and variance calculation that represent the actual design.
#'
#' @return A data.frame of class `mfrm_facet_design_effect` with columns
#'   `Facet`, `AvgClusterSize`, `ICC`, `DesignEffect`, `EffectiveN`, `InputRows`,
#'   `UsedRows`, and `ExcludedRows`. Its `data_usage` attribute is retained
#'   from the ICC result. Cluster sizes and effective sample sizes use only
#'   rows included in that model.
#'
#' @seealso [compute_facet_icc()], [analyze_hierarchical_structure()].
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. New York: Wiley.
#'
#' Park, I., & Lee, H. (2004). Design effects for the weighted mean and total
#' estimators under complex survey sampling. *Survey Methodology, 30*(2),
#' 183-193. \url{https://www150.statcan.gc.ca/n1/pub/12-001-x/2004002/article/7751-eng.pdf}
#'
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("compute_facet_design_effect", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   icc <- compute_facet_icc(toy, facets = c("Rater", "Criterion"),
#'                            score = "Score", person = "Person")
#'   deff <- compute_facet_design_effect(toy,
#'                                       facets = c("Rater", "Criterion"),
#'                                       icc_table = icc)
#'   print(deff)
#'   # Review clustering influence; EffectiveN is a descriptive row count.
#' }
#' }
#' @export
compute_facet_design_effect <- function(data, facets, icc_table = NULL,
                                        score = NULL, person = NULL,
                                        missing = c("error", "omit")) {
  missing <- match.arg(missing)
  if (!is.character(facets) || !length(facets) || anyNA(facets) ||
      anyDuplicated(facets) || any(!nzchar(facets))) {
    stop("`facets` must contain distinct, nonmissing column names.", call. = FALSE)
  }
  if (is.null(icc_table)) {
    if (is.null(score)) {
      stop("Supply either `icc_table` or `score`.", call. = FALSE)
    }
    icc_table <- compute_facet_icc(data, facets = facets,
                                   score = score, person = person, missing = missing)
  }
  usage <- attr(icc_table, "data_usage")
  if (is.null(usage$counts) || is.null(usage$observed_levels)) {
    stop("ICC row accounting is unavailable. Rerun compute_facet_icc() with the original data and settings before calculating design effects.", call. = FALSE)
  }
  if (any(!facets %in% names(usage$observed_levels)) || any(!facets %in% icc_table$Facet)) {
    stop("Each requested facet must be a grouping column in the supplied ICC result.", call. = FALSE)
  }
  total_n <- usage$counts[["UsedRows"]]
  out_rows <- lapply(facets, function(f) {
    k <- usage$observed_levels[[f]]
    avg_m <- if (k > 0) total_n / k else NA_real_
    rho <- suppressWarnings(as.numeric(
      icc_table$ICC[match(f, icc_table$Facet)]
    ))
    deff <- if (is.finite(rho) && is.finite(avg_m)) {
      1 + (avg_m - 1) * rho
    } else NA_real_
    eff_n <- if (is.finite(deff) && deff > 0) total_n / deff else NA_real_
    data.frame(
      Facet = f,
      AvgClusterSize = round(avg_m, 3),
      ICC = round(rho, 4),
      DesignEffect = round(deff, 3),
      EffectiveN = round(eff_n, 1),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out_rows)
  for (name in names(usage$counts)) out[[name]] <- usage$counts[[name]]
  structure(out, class = c("mfrm_facet_design_effect", "data.frame"), data_usage = usage)
}


# ---- 4. analyze_hierarchical_structure ------------------------------------

# Only reviewed columns matter; row order and unused factor levels do not.
.ha_data_identity <- function(data, facets, person, score = NULL) {
  columns <- unique(c(person, facets, score))
  values <- lapply(data[columns], function(x) {
    labels <- as.character(x)
    labels[is.na(x)] <- NA_character_
    labels
  })
  if (!is.null(score)) {
    value <- data[[score]]
    values[[score]] <- as.numeric(if (is.factor(value)) as.character(value) else value)
  }
  rows <- do.call(order, c(unname(values), list(na.last = TRUE, method = "radix")))
  values <- lapply(values, function(x) unname(x[rows]))
  list(version = 1L, facets = facets, person = person, score = score,
       rows = nrow(data), fingerprint = mfrm_checkpoint_fingerprint(values, canonical = TRUE))
}

.ha_validate_data_identity <- function(fit, review) {
  if (is.null(review)) return(invisible(TRUE))
  if (!inherits(review, "mfrm_hierarchical_structure")) {
    stop("`hierarchical_structure` must come from analyze_hierarchical_structure().", call. = FALSE)
  }
  source <- review$data_identity
  if (!is.list(source) || !identical(source$version, 1L) || is.null(source$fingerprint)) {
    stop("`hierarchical_structure` lacks its reviewed-data identity. Recreate it with ",
         "analyze_hierarchical_structure() from the saved fit or matching retained rows ",
         "and original review settings; no MFRM refit is needed.", call. = FALSE)
  }
  original <- fit$config$source_columns %||% fit$prep$source_columns %||% list()
  person_names <- unique(c("Person", original$person))
  score_names <- unique(c("Score", original$score))
  mapping <- stats::setNames(fit$prep$facet_names, fit$prep$facet_names)
  mapping[person_names] <- "Person"
  scope_matches <- identical(source$facets, review$facets) &&
    is.character(source$facets) && length(source$facets) >= 2L &&
    all(source$facets %in% names(mapping)) &&
    (is.null(source$person) || (length(source$person) == 1L && source$person %in% person_names)) &&
    (is.null(source$score) || (length(source$score) == 1L && source$score %in% score_names))
  mapping[score_names] <- "Score"
  columns <- unique(c(source$person, source$facets, source$score))
  matches <- FALSE
  if (isTRUE(scope_matches) && all(unname(mapping[columns]) %in% names(fit$prep$data))) {
    data <- as.data.frame(fit$prep$data[unname(mapping[columns])])
    names(data) <- columns
    matches <- identical(source, .ha_data_identity(data, source$facets, source$person, source$score))
  }
  if (!matches) {
    stop("fit/hierarchical_structure mismatch: the reviewed columns or rating rows ",
         "do not match this fit's retained data. Recreate the review from this fit ",
         "or matching retained rows. Reviews of other data remain usable on their own.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Analyze the hierarchical structure of a rating design
#'
#' One-stop review that combines the nesting, cross-tabulation, ICC, and
#' design-effect reports into a single object. Its summary feeds into
#' [reporting_checklist()] after checking the reviewed data against the fit.
#'
#' @param data Data frame in long format, or an `mfrm_fit` (its
#'   `prep$data` is used).
#' @param facets Character vector of facet column names. When `data` is
#'   an `mfrm_fit`, defaults to `fit$prep$facet_names`.
#' @param person Person column name. Defaults to `"Person"`.
#' @param score Score column name. Defaults to `"Score"`.
#' @param compute_icc Logical; if `TRUE` and `lme4` is available, adds
#'   ICC and design-effect tables.
#' @param ci_method ICC confidence-interval method passed through to
#'   [compute_facet_icc()]. One of `"none"` (default, point estimate
#'   only) or `"boot"`. The former `"profile"` method is refused because it
#'   did not calculate an ICC profile-likelihood interval. Deprecated alias:
#'   `icc_ci_method` (kept for backward compatibility, emits a
#'   lifecycle warning).
#' @param ci_level Confidence level when `ci_method != "none"`;
#'   default `0.95`. Deprecated alias: `icc_ci_level`.
#' @param ci_boot_reps Number of bootstrap replicates when
#'   `ci_method = "boot"`. Default `1000`. Deprecated alias:
#'   `icc_ci_boot_reps`.
#' @param ci_boot_seed Optional RNG seed for reproducible bootstrap
#'   CIs. Deprecated alias: `icc_ci_boot_seed`.
#' @param igraph_layout Logical; if `TRUE` and `igraph` is available,
#'   adds a connectivity component summary using a bipartite graph over
#'   person x facet levels.
#' @param icc_ci_method,icc_ci_level,icc_ci_boot_reps,icc_ci_boot_seed
#'   Deprecated compatibility spellings of the `ci_*` arguments above.
#'   Supplying a non-`NULL` value routes through
#'   [lifecycle::deprecate_warn()] and overrides the canonical
#'   `ci_*` argument.
#' @param missing Missing-value policy for the ICC and design-effect tables,
#'   passed to [compute_facet_icc()]. Invalid input or model errors stop the
#'   requested ICC analysis. Nesting, cross-tabulation, and connectivity tables
#'   describe the supplied design, including rows excluded from the ICC model.
#'   When `data` is a fit, ICC row counts start from its stored fitted rows;
#'   they cannot recover exclusions made before fitting the MFRM.
#'
#' @section Interpreting output:
#' - `nesting`: a
#'   [detect_facet_nesting()] object with every facet pair classified
#'   as Crossed / Partially / Near-perfectly / Fully nested.
#' - `crosstabs`: list of `(LevelA, LevelB, N)` long-format tables,
#'   one per facet pair. Plot via `plot(x, type = "crosstab",
#'   pair = "FacetA__FacetB")`.
#' - `icc`: per-facet variance shares. See
#'   [compute_facet_icc()] for the two-scale interpretation.
#' - `design_effect`: per-facet design-effect approximations and descriptive
#'   equivalent row counts; not precision estimates for the full design.
#' - `connectivity`: number of bipartite components linking
#'   Person x facet levels. A single component is required for a
#'   common measurement scale; multiple components indicate a
#'   disconnected design.
#'
#' @section Typical workflow:
#' 1. Optional: fit the MFRM with `fit_mfrm()`.
#' 2. Call `analyze_hierarchical_structure(fit)` (or on the raw data).
#' 3. Read `summary(x)` for the condensed view.
#' 4. Pass it as `hierarchical_structure` to [reporting_checklist()] and
#'    save the review with `saveRDS()` when its full tables are needed later.
#'    [build_mfrm_manifest()] records a separate fit-level rating-count
#'    review; it does not accept or archive this structural-review object.
#'
#' @section Saved-result reuse:
#' A reporting checklist checks the reviewed column roles and retained rating
#' rows, including scores when an ICC table was computed. It accepts a subset
#' of fitted facets and reports that scope explicitly. Reviews from raw data
#' also work when the selected columns and rows match, using either canonical
#' `Person`/`Score` names or the fit's original source-column names.
#' Row order, unused factor levels, unused columns and MFRM estimates do not
#' affect this match. Nesting and row counts are unweighted; when ICC is absent,
#' the match does not depend on scores. These rules do not establish that a
#' partial review covers all facets or that the design supports valid inference.
#' A review of additional rows, other strata or columns absent from the fit
#' remains usable on its own, but cannot mark this fit's checklist item as
#' available. Recreate older reviews without `data_identity` from the saved
#' fit or matching retained rows and original options; no MFRM refit is needed.
#' See [mfrmr_compatibility_layer].
#'
#' @return A list of class `mfrm_hierarchical_structure` with:
#' - `nesting`: output of [detect_facet_nesting()].
#' - `crosstabs`: list of pairwise observation-count data.frames (long
#'   format, suitable for heatmap plotting).
#' - `icc`: output of [compute_facet_icc()] when requested.
#' - `design_effect`: output of [compute_facet_design_effect()] when
#'   requested.
#' - `connectivity`: named list with bipartite-graph component summary
#'   when `igraph` is available.
#' - `summary`: one-row summary used by downstream reporting helpers.
#' - `facets`: character vector of facet names that were reviewed
#'   (echoed for downstream reporting helpers that need to label rows
#'   by review scope).
#' - `data_identity`: reviewed column roles, input row count and a compact
#'   fingerprint for reuse checks. The optional person column also contributes
#'   to nesting pair counts; `summary$NFacets` counts the requested `facets`.
#'
#' `plot(x, type = "icc")` draws variance shares and any available intervals
#' without an embedded title or footnote. Save its returned object as `p` and
#' use `plot_data(p, "notes")` to read the interval method or the reason bounds
#' are unavailable. These notes survive `saveRDS()`/`readRDS()`; they do not
#' change the underlying ICC estimates.
#'
#' @seealso [detect_facet_nesting()], [facet_small_sample_review()],
#'   [compute_facet_icc()], [compute_facet_design_effect()],
#'   [reporting_checklist()], [build_mfrm_manifest()], [fit_mfrm()].
#'
#' @concept confidence intervals
#' @concept hierarchical structure
#' @concept reporting workflow
#'
#' @references
#' McEwen, M. R. (2018). *The effects of incomplete rating designs on
#' results from many-facets-Rasch model analyses* (Doctoral thesis,
#' Brigham Young University). <https://scholarsarchive.byu.edu/etd/6689/>
#'
#' Linacre, J. M. (2026). *A User's Guide to FACETS, Version 4.5.0*.
#' Winsteps.com.
#'
#' Kish, L. (1965). *Survey Sampling*. New York: Wiley.
#'
#' Koo, T. K., & Li, M. Y. (2016). A guideline of selecting and
#' reporting intraclass correlation coefficients for reliability
#' research. *Journal of Chiropractic Medicine, 15*(2), 155-163.
#'
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("analyze_hierarchical_structure", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' hs <- analyze_hierarchical_structure(toy,
#'                                      facets = c("Rater", "Criterion"),
#'                                      compute_icc = FALSE,
#'                                      igraph_layout = FALSE)
#' summary(hs)
#'
#' # Full review when lme4 and igraph are available.
#' if (requireNamespace("lme4", quietly = TRUE) &&
#'     requireNamespace("igraph", quietly = TRUE)) {
#'   hs_full <- analyze_hierarchical_structure(toy,
#'                                             facets = c("Rater", "Criterion"))
#'   summary(hs_full)
#'   plot(hs_full, type = "icc")
#' }
#' }
#' @export
analyze_hierarchical_structure <- function(data,
                                           facets = NULL,
                                           person = "Person",
                                           score = "Score",
                                           compute_icc = TRUE,
                                           ci_method = c("none", "boot"),
                                           ci_level = 0.95,
                                           ci_boot_reps = 1000L,
                                           ci_boot_seed = NULL,
                                           igraph_layout = TRUE,
                                           icc_ci_method = NULL,
                                           icc_ci_level = NULL,
                                           icc_ci_boot_reps = NULL,
                                           icc_ci_boot_seed = NULL,
                                           missing = c("error", "omit")) {
  missing <- match.arg(missing)
  data_source <- if (inherits(data, "mfrm_fit")) "Stored fitted rows" else "Supplied data"
  # Deprecated `icc_ci_*` spellings route through lifecycle and
  # override the canonical `ci_*` values when supplied. This unifies
  # the API with compute_facet_icc() while preserving compatibility.
  if (!is.null(icc_ci_method)) {
    lifecycle::deprecate_warn(
      when = "0.1.6",
      what = "analyze_hierarchical_structure(icc_ci_method = )",
      with = "analyze_hierarchical_structure(ci_method = )"
    )
    ci_method <- icc_ci_method
  }
  if (!is.null(icc_ci_level)) {
    lifecycle::deprecate_warn(
      when = "0.1.6",
      what = "analyze_hierarchical_structure(icc_ci_level = )",
      with = "analyze_hierarchical_structure(ci_level = )"
    )
    ci_level <- icc_ci_level
  }
  if (!is.null(icc_ci_boot_reps)) {
    lifecycle::deprecate_warn(
      when = "0.1.6",
      what = "analyze_hierarchical_structure(icc_ci_boot_reps = )",
      with = "analyze_hierarchical_structure(ci_boot_reps = )"
    )
    ci_boot_reps <- icc_ci_boot_reps
  }
  if (!is.null(icc_ci_boot_seed)) {
    lifecycle::deprecate_warn(
      when = "0.1.6",
      what = "analyze_hierarchical_structure(icc_ci_boot_seed = )",
      with = "analyze_hierarchical_structure(ci_boot_seed = )"
    )
    ci_boot_seed <- icc_ci_boot_seed
  }
  ci_method <- .icc_ci_method(ci_method)
  if (inherits(data, "mfrm_fit")) {
    fit_ref <- data
    if (is.null(facets)) facets <- fit_ref$prep$facet_names
    data <- fit_ref$prep$data
  }
  facets <- as.character(facets)
  if (length(facets) < 2L) {
    stop("`analyze_hierarchical_structure()` needs at least two facets.",
         call. = FALSE)
  }

  # 1. Nesting
  nesting <- detect_facet_nesting(data, facets, person = person)

  # 2. Cross-tabulations (long format)
  crosstabs <- list()
  pair_idx <- 1L
  for (i in seq_len(length(facets) - 1L)) {
    for (j in seq(i + 1L, length(facets))) {
      a <- facets[i]; b <- facets[j]
      ctab <- data |>
        dplyr::count(LevelA = as.character(.data[[a]]),
                     LevelB = as.character(.data[[b]]), name = "N") |>
        dplyr::mutate(FacetA = a, FacetB = b)
      crosstabs[[paste(a, b, sep = "__")]] <- as.data.frame(ctab,
                                                            stringsAsFactors = FALSE)
      pair_idx <- pair_idx + 1L
    }
  }

  # 3. ICC and design effect
  icc_tbl <- NULL
  deff_tbl <- NULL
  icc_available <- isTRUE(compute_icc) &&
    requireNamespace("lme4", quietly = TRUE) &&
    !is.null(score)
  if (icc_available) {
    icc_tbl <- compute_facet_icc(data, facets = facets, score = score,
                        person = person,
                        ci_method = ci_method,
                        ci_level = ci_level,
                        ci_boot_reps = ci_boot_reps,
                        ci_boot_seed = ci_boot_seed, missing = missing)
    attr(icc_tbl, "data_usage")$source <- data_source
    if (!is.null(icc_tbl) && nrow(icc_tbl) > 0) {
      deff_tbl <- compute_facet_design_effect(data, facets = facets,
                                    icc_table = icc_tbl,
                                    score = score, person = person)
    }
  }

  # 4. Connectivity via bipartite graph
  connectivity <- NULL
  if (isTRUE(igraph_layout) &&
      requireNamespace("igraph", quietly = TRUE) &&
      !is.null(person) && person %in% names(data)) {
    edges <- dplyr::distinct(
      data[, c(person, facets), drop = FALSE]
    )
    el <- do.call(rbind, lapply(facets, function(f) {
      data.frame(
        from = paste0("P:", as.character(edges[[person]])),
        to = paste0(f, ":", as.character(edges[[f]])),
        stringsAsFactors = FALSE
      )
    }))
    el <- dplyr::distinct(el)
    g <- igraph::graph_from_data_frame(el, directed = FALSE)
    comps <- igraph::components(g)
    connectivity <- list(
      n_components = as.integer(comps$no),
      largest_component_size = as.integer(max(comps$csize)),
      component_sizes = as.integer(comps$csize),
      isolates = sum(comps$csize == 1L)
    )
  }

  # 5. Summary
  any_sparse <- if (!is.null(icc_tbl) && nrow(icc_tbl) > 0) {
    any(icc_tbl$ICC > 0.10, na.rm = TRUE)
  } else NA
  summary_tbl <- data.frame(
    NFacets = length(facets),
    NestedPairs = if (!is.null(nesting$summary$FullyNestedPairs)) {
      nesting$summary$FullyNestedPairs
    } else 0L,
    CrossedPairs = if (!is.null(nesting$summary$CrossedPairs)) {
      nesting$summary$CrossedPairs
    } else 0L,
    ICCAvailable = !is.null(icc_tbl) && nrow(icc_tbl) > 0,
    ConnectivityComponents = if (!is.null(connectivity)) {
      connectivity$n_components
    } else NA_integer_,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      nesting = nesting,
      crosstabs = crosstabs,
      icc = icc_tbl,
      design_effect = deff_tbl,
      connectivity = connectivity,
      summary = summary_tbl,
      facets = facets,
      data_identity = .ha_data_identity(data, facets, person,
                                       score = if (icc_available) score else NULL)
    ),
    class = "mfrm_hierarchical_structure"
  )
}


# ---- S3 methods -----------------------------------------------------------

#' @export
print.mfrm_facet_nesting <- function(x, ...) {
  cat("mfrm_facet_nesting\n")
  cat("  Facets reviewed:", paste(x$facets, collapse = ", "), "\n")
  if (nrow(x$pairwise_table) > 0) {
    cat("  Pairs:", nrow(x$pairwise_table), "\n")
    cat("  Any nested pair:",
        isTRUE(x$summary$AnyNested), "\n")
  }
  cat("Use `summary(x)` for the full pairwise table.\n")
  invisible(x)
}

#' @export
summary.mfrm_facet_nesting <- function(object, ...) {
  cat("mfrm_facet_nesting\n\n")
  cat("Summary:\n")
  print(object$summary, row.names = FALSE)
  if (nrow(object$pairwise_table) > 0) {
    cat("\nPairwise nesting:\n")
    display_cols <- c("FacetA", "FacetB", "LevelsA", "LevelsB",
                      "NestingIndex_AinB", "NestingIndex_BinA",
                      "Direction")
    print(object$pairwise_table[, display_cols, drop = FALSE],
          row.names = FALSE)
  }
  invisible(object)
}

#' @rdname facet_small_sample_review
#' @details `plot()` method: Plot a facet sample-size review
#'
#' Per-level observation counts rendered as a horizontal bar chart
#' coloured by the count band assigned in
#' [facet_small_sample_review()]. Vertical dashed lines mark the
#' sparse / marginal / standard thresholds so reviewers see where
#' every facet level sits relative to the selected screening thresholds.
#' Counts and colours do not establish measurement precision.
#' @param x An `mfrm_facet_sample_review` object.
#' @param top_n Optional integer; trim the y-axis to the `top_n`
#'   smallest level counts per facet. `NULL` (default) keeps all.
#' @param preset One of `"standard"`, `"publication"`, `"compact"`, `"monochrome"`.
#' @param ... Reserved.
#' @return For `plot()`: Invisibly, the data.frame used for the plot, with the review's
#'   interpretation limits in its `notes` attribute.
#' @seealso [facet_small_sample_review()].
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_facet_sample_review <- function(x, top_n = NULL,
                                          preset = c("standard",
                                                     "publication",
                                                     "compact",
                                                     "monochrome"),
                                          ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  style <- resolve_plot_preset(preset)
  tbl <- as.data.frame(x$table, stringsAsFactors = FALSE)
  attr(tbl, "notes") <- x$notes
  if (is.null(tbl) || nrow(tbl) == 0L) {
    plot.new()
    text(0.5, 0.5, "Review table empty.")
    return(invisible(tbl))
  }
  band_colors <- c(
    sparse   = "#D73027",
    marginal = "#FDAE61",
    standard = "#66BD63",
    strong   = "#1A9850"
  )
  tbl$BarColor <- band_colors[tbl$SampleCategory]
  tbl$BarColor[is.na(tbl$BarColor)] <- style$neutral
  if (!is.null(top_n) && is.finite(top_n) && top_n > 0) {
    tbl <- tbl |>
      dplyr::group_by(Facet) |>
      dplyr::slice_min(order_by = .data$N, n = as.integer(top_n),
                       with_ties = FALSE) |>
      dplyr::ungroup() |>
      as.data.frame(stringsAsFactors = FALSE)
  }
  tbl <- tbl[order(tbl$Facet, tbl$N), , drop = FALSE]
  labels <- paste0(tbl$Facet, " / ", tbl$Level)

  thr <- x$thresholds %||% list(sparse = 10, marginal = 30, standard = 50)
  xmax <- max(c(tbl$N, as.numeric(thr$standard), 1), na.rm = TRUE) * 1.1

  old_par <- graphics::par()["mar"]
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4, max(8, min(18, max(nchar(labels)) * 0.55)), 3, 1))
  barplot(
    height = tbl$N,
    names.arg = labels,
    horiz = TRUE,
    las = 1,
    xlim = c(0, xmax),
    col = tbl$BarColor,
    border = NA,
    xlab = "Retained rating rows per level",
    main = "",
    cex.names = 0.8
  )
  abline(
    v = c(as.numeric(thr$sparse),
          as.numeric(thr$marginal),
          as.numeric(thr$standard)),
    lty = 2, col = style$neutral
  )
  legend(
    "bottomright", bty = "n", cex = 0.8,
    legend = c(sprintf("< %s sparse", thr$sparse),
               sprintf("< %s marginal", thr$marginal),
               sprintf("< %s standard", thr$standard),
               "strong"),
    fill = unname(band_colors)
  )
  attr(tbl, "notes") <- x$notes
  invisible(tbl)
}

#' Plot the pairwise nesting index matrix
#'
#' Renders the directed nesting index
#' \eqn{1 - H(B \mid A)/H(B)} as a heatmap between facet pairs,
#' highlighting fully nested relationships close to 1. Colour scale
#' runs from 0 (crossed, white / cold) to 1 (fully nested, dark).
#'
#' @param x An `mfrm_facet_nesting` object.
#' @param preset Plot preset.
#' @param ... Reserved.
#' @return Invisibly, the matrix rendered.
#' @seealso [detect_facet_nesting()],
#'   [analyze_hierarchical_structure()].
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @export
plot.mfrm_facet_nesting <- function(x,
                                    preset = c("standard",
                                               "publication",
                                               "compact",
                                               "monochrome"),
                                    ...) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  style <- resolve_plot_preset(preset)
  pair <- as.data.frame(x$pairwise_table, stringsAsFactors = FALSE)
  if (is.null(pair) || nrow(pair) == 0L) {
    plot.new()
    text(0.5, 0.5, "At least two facets required.")
    return(invisible(NULL))
  }
  facets <- x$facets
  n <- length(facets)
  m <- matrix(NA_real_, nrow = n, ncol = n,
              dimnames = list(facets, facets))
  for (i in seq_len(nrow(pair))) {
    a <- pair$FacetA[i]; b <- pair$FacetB[i]
    # nesting index A in B at m[a, b]
    m[a, b] <- pair$NestingIndex_AinB[i]
    m[b, a] <- pair$NestingIndex_BinA[i]
  }
  diag(m) <- 1

  cols <- grDevices::hcl.colors(20, palette = "Blues 3", rev = TRUE)
  old_par <- graphics::par()["mar"]
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(5, 5, 3, 2))
  image(
    x = seq_len(n), y = seq_len(n),
    z = m,
    col = cols,
    zlim = c(0, 1),
    xaxt = "n", yaxt = "n",
    xlab = "Nested in (column facet)",
    ylab = "Nested facet (row)",
    main = ""
  )
  axis(1, at = seq_len(n), labels = facets, las = 2, cex.axis = 0.8)
  axis(2, at = seq_len(n), labels = facets, las = 1, cex.axis = 0.8)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (is.finite(m[i, j])) {
      text(j, i, sprintf("%.2f", m[i, j]), cex = 0.75,
                     col = if (m[i, j] > 0.6) "white" else "black")
    }
  }
  invisible(m)
}

#' @export
print.mfrm_facet_sample_review <- function(x, ...) {
  cat("mfrm_facet_sample_review\n")
  if (length(x$notes)) cat(paste(x$notes, collapse = "\n"), "\n")
  cat("  Thresholds (sparse / marginal / standard):",
      paste(unlist(x$thresholds), collapse = " / "), "\n")
  cat("  Facets:", nrow(x$facet_summary), "\n")
  cat("  Sparse levels total:",
      sum(x$summary$sparse, na.rm = TRUE), "\n")
  cat("Use `summary(x)` for the detailed breakdown.\n")
  invisible(x)
}

#' @export
summary.mfrm_facet_sample_review <- function(object, ...) {
  cat("mfrm_facet_sample_review\n\n")
  if (length(object$notes)) cat(paste(object$notes, collapse = "\n"), "\n\n")
  cat("Per-facet summary:\n")
  print(object$facet_summary, row.names = FALSE)
  cat("\nSample-size category counts by facet:\n")
  print(object$summary, row.names = FALSE)
  sparse_rows <- object$table[object$table$SampleCategory == "sparse", ,
                              drop = FALSE]
  if (nrow(sparse_rows) > 0) {
    cat("\nSparse levels (n <", object$thresholds$sparse, "):\n")
    print(sparse_rows[, c("Facet", "Level", "N",
                          "Estimate", "SE", "SampleCategory")],
          row.names = FALSE)
  }
  invisible(object)
}

.print_icc_data_usage <- function(x) {
  usage <- attr(x, "data_usage", exact = TRUE)
  if (is.null(usage)) {
    cat("  Row accounting unavailable; rerun the ICC analysis to obtain it.\n")
  } else {
    cat(sprintf("  ICC rows: %d input, %d used, %d excluded.\n",
                usage$counts[["InputRows"]], usage$counts[["UsedRows"]],
                usage$counts[["ExcludedRows"]]))
    if (identical(usage$source, "Stored fitted rows")) {
      cat("  Counts start from stored fitted rows; earlier MFRM filtering is not included.\n")
    }
    if (usage$counts[["ExcludedRows"]] > 0L) {
      cat("  Incomplete rows were explicitly omitted; no missing values were imputed.\n")
    }
  }
}

#' @export
print.mfrm_facet_icc <- function(x, ...) {
  .check_icc_intervals(x)
  cat("mfrm_facet_icc\n")
  .print_icc_data_usage(x)
  if (nrow(x) == 0L) {
    cat("  (empty; lme4 unavailable or fit failed)\n")
  } else {
    print.data.frame(x, row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.mfrm_facet_icc <- function(object, ...) {
  .check_icc_intervals(object)
  # Condensed view that separates the two interpretation scales so
  # readers don't conflate person reliability with non-person variance
  # share; see `compute_facet_icc()` "Interpreting output".
  cat("Facet ICC summary (mfrmr)\n")
  .print_icc_data_usage(object)
  if (!is.data.frame(object) || nrow(object) == 0L) {
    cat("  (empty; lme4 unavailable or fit failed)\n")
    return(invisible(object))
  }
  scales <- if ("InterpretationScale" %in% names(object)) {
    split(object, object$InterpretationScale)
  } else {
    list(`(unscaled)` = object)
  }
  for (nm in names(scales)) {
    cat(sprintf("  -- %s --\n", nm))
    print.data.frame(scales[[nm]], row.names = FALSE)
  }
  invisible(object)
}

#' @export
print.mfrm_facet_design_effect <- function(x, ...) {
  cat("mfrm_facet_design_effect (per-facet approximation)\n")
  cat("  EffectiveN is a descriptive row count, not full-design precision.\n")
  .print_icc_data_usage(x)
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#' @export
summary.mfrm_facet_design_effect <- function(object, ...) {
  cat("Per-facet design-effect approximations (mfrmr)\n")
  cat("  EffectiveN is a descriptive row count, not full-design precision.\n")
  .print_icc_data_usage(object)
  if (!is.data.frame(object) || nrow(object) == 0L) {
    cat("  (empty)\n")
    return(invisible(object))
  }
  worst <- which.max(suppressWarnings(as.numeric(object$DesignEffect)))
  cat("  Largest DesignEffect: ",
      if (length(worst) == 1L) {
        sprintf("%s = %.2f (EffectiveN = %.1f)",
                object$Facet[worst],
                as.numeric(object$DesignEffect[worst]),
                as.numeric(object$EffectiveN[worst]))
      } else "NA", "\n", sep = "")
  print.data.frame(object, row.names = FALSE)
  invisible(object)
}

#' @export
print.mfrm_hierarchical_structure <- function(x, ...) {
  cat("mfrm_hierarchical_structure\n")
  cat("  Facets:", paste(x$facets, collapse = ", "), "\n")
  cat("  Nested pairs:",
      x$summary$NestedPairs %||% 0L, "\n")
  cat("  Crossed pairs:",
      x$summary$CrossedPairs %||% 0L, "\n")
  if (isTRUE(x$summary$ICCAvailable)) {
    cat("  ICC table: available (", nrow(x$icc), " facets)\n", sep = "")
    .print_icc_data_usage(x$icc)
  } else {
    cat("  ICC table: unavailable (install `lme4` or set compute_icc = FALSE)\n")
  }
  if (!is.null(x$connectivity)) {
    cat("  Connectivity components:",
        x$connectivity$n_components, "\n")
  }
  cat("Use `summary(x)` for the full report.\n")
  invisible(x)
}

#' @export
summary.mfrm_hierarchical_structure <- function(object, ...) {
  if (!is.null(object$icc)) .check_icc_intervals(object$icc)
  cat("mfrm_hierarchical_structure\n\n")
  cat("Summary:\n")
  print(object$summary, row.names = FALSE)

  cat("\nNesting review:\n")
  print(object$nesting$pairwise_table[,
    c("FacetA", "FacetB", "NestingIndex_AinB",
      "NestingIndex_BinA", "Direction"),
    drop = FALSE], row.names = FALSE)

  if (!is.null(object$icc) && nrow(object$icc) > 0) {
    cat("\nICC (lme4 variance-components):\n")
    print(object$icc, row.names = FALSE)
  }
  if (!is.null(object$design_effect) && nrow(object$design_effect) > 0) {
    cat("\nDesign-effect approximations:\n")
    print(object$design_effect, row.names = FALSE)
  }
  if (!is.null(object$connectivity)) {
    cat("\nBipartite connectivity (via igraph):\n")
    cat("  Components:", object$connectivity$n_components,
        "\n  Largest component:", object$connectivity$largest_component_size,
        "\n  Isolates:", object$connectivity$isolates, "\n")
  }
  invisible(object)
}

#' @export
plot.mfrm_hierarchical_structure <- function(x, type = c("crosstab", "icc"),
                                             pair = NULL, ...) {
  type <- match.arg(type)
  # A saved return value can be plotted again with a different view.
  attr(x, "mfrmr_notes") <- NULL
  if (type == "crosstab") {
    if (length(x$crosstabs) == 0L) {
      stop("No cross-tabulations available.", call. = FALSE)
    }
    pair_name <- if (!is.null(pair)) {
      if (!pair %in% names(x$crosstabs)) {
        stop("Unknown pair: ", pair,
             ". Available: ", paste(names(x$crosstabs), collapse = ", "),
             call. = FALSE)
      }
      pair
    } else {
      names(x$crosstabs)[1L]
    }
    tbl <- x$crosstabs[[pair_name]]
    mat <- tidyr::pivot_wider(tbl[, c("LevelA", "LevelB", "N")],
                              names_from = LevelB, values_from = N,
                              values_fill = 0L)
    rnames <- as.character(mat$LevelA)
    mat_num <- as.matrix(mat[, -1L, drop = FALSE])
    rownames(mat_num) <- rnames
    image(
      x = seq_len(nrow(mat_num)),
      y = seq_len(ncol(mat_num)),
      z = mat_num,
      xaxt = "n", yaxt = "n",
      xlab = tbl$FacetA[1L],
      ylab = tbl$FacetB[1L],
      main = "",
      col = grDevices::hcl.colors(20, palette = "YlGnBu", rev = TRUE)
    )
    axis(1L, at = seq_len(nrow(mat_num)),
                   labels = rownames(mat_num), las = 2L, cex.axis = 0.7)
    axis(2L, at = seq_len(ncol(mat_num)),
                   labels = colnames(mat_num), las = 2L, cex.axis = 0.7)
  } else if (type == "icc") {
    if (is.null(x$icc) || nrow(x$icc) == 0L) {
      stop("No ICC table available; re-run with compute_icc = TRUE and lme4 installed.",
           call. = FALSE)
    }
    icc_tbl <- x$icc
    .check_icc_intervals(icc_tbl)
    has_ci <- all(c("ICC_CI_Lower", "ICC_CI_Upper") %in% names(icc_tbl)) &&
      any(is.finite(icc_tbl$ICC_CI_Lower) & is.finite(icc_tbl$ICC_CI_Upper))
    ci_level <- if (has_ci && "ICC_CI_Level" %in% names(icc_tbl)) {
      suppressWarnings(as.numeric(icc_tbl$ICC_CI_Level[1L]))
    } else NA_real_
    y_max <- max(
      1,
      max(c(icc_tbl$ICC, icc_tbl$ICC_CI_Upper), na.rm = TRUE) * 1.1
    )
    main_txt <- "Facet ICC (variance component share)"
    if (has_ci && is.finite(ci_level)) {
      main_txt <- sprintf("%s\n%g%% parametric bootstrap CI",
                          main_txt, 100 * ci_level)
    }
    missing_ci <- !is.finite(icc_tbl$ICC_CI_Lower) | !is.finite(icc_tbl$ICC_CI_Upper)
    if (any(missing_ci & icc_tbl$ICC_CI_Method == "boot", na.rm = TRUE)) {
      main_txt <- paste0(main_txt, "\nIntervals unavailable: ",
                         paste(unique(icc_tbl$ICC_CI_Status[missing_ci]), collapse = "; "))
    }
    attr(x, "mfrmr_notes") <- data.frame(
      Type = "description", Text = main_txt, stringsAsFactors = FALSE)
    mids <- barplot(
      height = icc_tbl$ICC,
      names.arg = icc_tbl$Facet,
      main = "",
      ylab = "ICC",
      ylim = c(0, y_max)
    )
    if (has_ci) {
      valid <- is.finite(icc_tbl$ICC_CI_Lower) & is.finite(icc_tbl$ICC_CI_Upper)
      nonzero <- valid & icc_tbl$ICC_CI_Upper > icc_tbl$ICC_CI_Lower
      if (any(nonzero)) {
        arrows(
          x0 = mids[nonzero], y0 = icc_tbl$ICC_CI_Lower[nonzero],
          x1 = mids[nonzero], y1 = icc_tbl$ICC_CI_Upper[nonzero],
          angle = 90, code = 3, length = 0.05, col = "black", lwd = 1.5
        )
      }
      flat <- valid & icc_tbl$ICC_CI_Upper == icc_tbl$ICC_CI_Lower
      if (any(flat)) {
        segments(mids[flat] - 0.05, icc_tbl$ICC_CI_Lower[flat],
                           mids[flat] + 0.05, icc_tbl$ICC_CI_Upper[flat], lwd = 1.5)
      }
    }
  }
  invisible(x)
}
