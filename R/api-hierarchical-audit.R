# ==============================================================================
# Hierarchical structure and small-sample review (added in 0.1.6)
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
#   and differential severity/leniency. 0.1.6 adds only the review layer
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
.ha_entropy <- function(x) {
  # Shannon entropy H(X) in nats using observed frequencies. NA-safe:
  # NA levels are dropped before counting.
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(0)
  freq <- tabulate(match(x, unique(x)))
  p <- freq / sum(freq)
  p <- p[p > 0]
  -sum(p * log(p))
}

#' @keywords internal
#' @noRd
.ha_conditional_entropy <- function(y, x) {
  # H(Y | X) = sum_x p(x) * H(Y | X = x)
  keep <- !is.na(x) & !is.na(y)
  y <- y[keep]; x <- x[keep]
  if (length(x) == 0L) return(0)
  split_y <- split(y, x)
  weights <- lengths(split_y) / length(x)
  sum(weights * vapply(split_y, .ha_entropy, numeric(1)))
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

#' @keywords internal
#' @noRd
.ha_extract_fit_data <- function(fit) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  prep <- fit$prep
  data <- prep$data
  list(
    data = data,
    person_col = "Person",
    facets = prep$facet_names,
    score_col = "Score"
  )
}


# ---- 1. detect_facet_nesting ----------------------------------------------

#' Detect nesting structure between facets
#'
#' Classifies every ordered pair of facets (optionally including `Person`)
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
#' @param weight_col Optional name of a weight column; if supplied, rows
#'   are replicated proportionally when counting element co-occurrences.
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
#'    [build_mfrm_manifest()] then records for reproducibility.
#'
#' @return A list of class `mfrm_facet_nesting` with:
#' - `pairwise_table`: one row per ordered facet pair with
#'   `NestingIndex_AinB`, `NestingIndex_BinA`, classification strings,
#'   and `Direction`.
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

  all_cols <- c(if (!is.null(person)) person else NULL, facets)
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

  # Apply weight column by replicating rows when present (so entropy
  # reflects weighted frequencies).
  if (!is.null(weight_col) && weight_col %in% names(data)) {
    w <- suppressWarnings(as.numeric(data[[weight_col]]))
    w <- ifelse(is.finite(w) & w > 0, pmax(1L, round(w)), 1L)
    data <- data[rep(seq_len(nrow(data)), times = w), , drop = FALSE]
  }

  pairs <- list()
  idx <- 1L
  for (i in seq_len(n_cols - 1L)) {
    for (j in seq(i + 1L, n_cols)) {
      a <- all_cols[i]; b <- all_cols[j]
      ha <- .ha_entropy(data[[a]])
      hb <- .ha_entropy(data[[b]])
      hb_given_a <- .ha_conditional_entropy(data[[b]], data[[a]])
      ha_given_b <- .ha_conditional_entropy(data[[a]], data[[b]])
      nesting_a_in_b <- if (hb > 0) 1 - hb_given_a / hb else NA_real_
      nesting_b_in_a <- if (ha > 0) 1 - ha_given_b / ha else NA_real_
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
        LevelsA = dplyr::n_distinct(data[[a]], na.rm = TRUE),
        LevelsB = dplyr::n_distinct(data[[b]], na.rm = TRUE),
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

#' Review per-facet-level sample adequacy
#'
#' Reports per-level observation counts, SE, and fit statistics for every
#' level of every facet in a fitted MFRM model, and classifies each level
#' as `"sparse"`, `"marginal"`, `"standard"`, or `"strong"` against the
#' Linacre sample-size bands.
#'
#' In mfrmr every facet is a fixed effect (see `?fit_mfrm`, "Fixed
#' effects assumption"), so a level with very few ratings contributes an
#' estimate with wide SE but no shrinkage toward the facet mean. This
#' helper surfaces those levels up front so users can decide whether to
#' drop them, pool them, or move to a hierarchical model outside mfrmr.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output. When supplied,
#'   per-level `Infit`, `Outfit`, and `ModelSE` are added to the report.
#' @param thresholds Named numeric vector of count bands. Defaults are
#'   `c(sparse = 10, marginal = 30, standard = 50)`. These are adapted
#'   from Linacre (1994): the 30-level band preserves Linacre's
#'   approximately `+-1.0 logit at 95% CI` line, while the `sparse < 10`
#'   floor and the `standard = 50` watermark are mfrmr-specific screening
#'   choices below Linacre's 30-examinee minimum and between Linacre's
#'   30 and 100 thresholds.
#'
#' @section Interpreting output:
#' - `"sparse"` (n < 10): level-level estimate is unstable; SE will be
#'   wide; consider combining with adjacent levels or treating as
#'   exploratory only.
#' - `"marginal"` (10 <= n < 30): below Linacre (1994) 95% CI
#'   +-1.0 logit threshold; usable as screening only.
#' - `"standard"` (30 <= n < 50): meets baseline stability; reasonable
#'   for publication if fit statistics are acceptable.
#' - `"strong"` (n >= 50): well-targeted; facet estimate is robust.
#'
#' Because mfrmr has no shrinkage by default, sparse and marginal levels
#' do not "borrow strength" from other levels. Jones and Wind (2018)
#' report that rater estimates are particularly sensitive to thin
#' linking; the `Facet = "Person"` row is usually less of a concern
#' because the person prior integrates out the uncertainty.
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
#' - `facet_summary`: smallest observed level count per facet.
#' - `thresholds`: the applied count bands.
#'
#' @seealso [detect_facet_nesting()], [analyze_hierarchical_structure()],
#'   [compute_facet_icc()], [compute_facet_design_effect()],
#'   [reporting_checklist()].
#'
#' @references
#' Linacre, J. M. (2026). *A User's Guide to FACETS, Version 4.5.0*.
#' Winsteps.com. <https://www.winsteps.com/facets.htm>
#'
#' Linacre, J. M. (1994). Sample size and item calibration stability.
#' *Rasch Measurement Transactions, 7*(4), 328.
#' <https://www.rasch.org/rmt/rmt74m.htm>
#'
#' Jones, E., & Wind, S. A. (2018). Using repeated ratings to improve
#' measurement precision in incomplete rating designs. *Journal of
#' Applied Measurement, 19*(2), 148-161.
#'
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
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
facet_small_sample_review <- function(fit, diagnostics = NULL,
                                      thresholds = c(sparse = 10,
                                                     marginal = 30,
                                                     standard = 50)) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (is.null(names(thresholds)) ||
      !all(c("sparse", "marginal", "standard") %in% names(thresholds))) {
    stop("`thresholds` must be a named numeric vector containing ",
         "`sparse`, `marginal`, and `standard`.", call. = FALSE)
  }

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
        max(match(SampleCategory, category_order), na.rm = TRUE)
      ],
      .groups = "drop"
    )

  structure(
    list(
      table = as.data.frame(audit_tbl, stringsAsFactors = FALSE),
      summary = as.data.frame(level_count_by_facet, stringsAsFactors = FALSE),
      facet_summary = as.data.frame(facet_summary, stringsAsFactors = FALSE),
      thresholds = as.list(thresholds)
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
#'   One of `"none"` (default, point estimate only), `"profile"`
#'   (a **first-order approximation**: marginal likelihood-profile
#'   bounds for each variance-component SD via
#'   [lme4::confint.merMod()] with `method = "profile"`, squared to
#'   variances, then plugged into the ICC ratio while holding the
#'   other components at their point estimate; fast and deterministic,
#'   the default recommendation for reporting), or `"boot"`
#'   (parametric bootstrap via [lme4::bootMer()]; slower but robust
#'   to non-normal ICC sampling distributions because each bootstrap
#'   replicate resamples the full variance decomposition jointly).
#' @param ci_level Confidence level when `ci_method != "none"`; default
#'   `0.95`. Koo & Li (2016) recommend banding the CI rather than the
#'   point estimate when classifying reliability as Poor / Moderate /
#'   Good / Excellent.
#' @param ci_boot_reps Number of bootstrap replicates used when
#'   `ci_method = "boot"`. Default `1000`.
#' @param ci_boot_seed Optional integer seed for the bootstrap path
#'   (`NULL` leaves the RNG state untouched).
#' @param ci_boot_parallel Parallelisation strategy for the
#'   parametric-bootstrap CI path, passed through to
#'   [lme4::bootMer()]: `"no"` (default), `"multicore"` (POSIX
#'   `mclapply`), or `"snow"` (PSOCK cluster). `"multicore"` does
#'   nothing on Windows and falls back to serial; in that case use
#'   `"snow"` with [parallel::makeCluster()] in scope.
#' @param ci_boot_ncpus Number of CPUs to use for the parallel
#'   bootstrap path (ignored when `ci_boot_parallel = "no"`). The
#'   per-replicate progress bar is suppressed under parallel
#'   execution because worker processes cannot push updates to the
#'   parent's cli console.
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
#'
#' Note: Koo & Li (2016) recommend applying the reliability bands to
#' the **95% confidence interval** of the ICC rather than to the point
#' estimate alone. Set `ci_method = "profile"` (default `"none"`) to
#' obtain likelihood-profile CI bounds alongside the point estimate,
#' or `ci_method = "boot"` for a parametric bootstrap with
#' `ci_boot_reps` replicates. The returned data frame gains
#' `ICC_CI_Lower` / `ICC_CI_Upper` columns so downstream reporting can
#' apply the band to the CI rather than the point estimate. The
#' `Interpretation` column still uses the point estimate so
#' callers who want CI-aware banding can implement it externally from
#' the supplied bounds.
#'
#' @section Typical workflow:
#' 1. Fit the MFRM model with `fit_mfrm()` for the Rasch-metric
#'    separation/reliability.
#' 2. Call `compute_facet_icc(data, facets, score, person)` to get the
#'    complementary variance-share summary.
#' 3. Feed into [compute_facet_design_effect()] to convert ICCs and
#'    average cluster sizes into Kish (1965) design effects.
#'
#' @return A data.frame of class `mfrm_facet_icc` with one row per
#'   variance component (including a `"Residual"` row) and columns:
#' - `Facet`: the grouping factor name (or `"Residual"`).
#' - `Variance`: REML variance estimate.
#' - `ICC`: variance share (`Variance / sum(Variance)`), in `[0, 1]`.
#' - `Interpretation`: band label according to the facet's scale.
#' - `InterpretationScale`: `"Koo-Li reliability"` for the person
#'   facet, `"Variance share"` for others.
#' - `ICC_CI_Lower` / `ICC_CI_Upper` / `ICC_CI_Level` / `ICC_CI_Method`:
#'   CI bounds, level, and method (populated when `ci_method != "none"`;
#'   `NA_real_` otherwise).
#' - `ICC_CI_NReps`: bootstrap replicate count when
#'   `ci_method = "boot"` (absent otherwise).
#'
#' @seealso [compute_facet_design_effect()],
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
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   icc <- compute_facet_icc(toy, facets = c("Rater", "Criterion"),
#'                            score = "Score", person = "Person")
#'   print(icc)
#'   # Look for:
#'   # - Person ICC reads as Koo & Li (2016) reliability: < 0.5 poor,
#'   #   0.5-0.75 moderate, 0.75-0.9 good, > 0.9 excellent.
#'   # - Rater / Criterion ICC reads as variance share, NOT reliability;
#'   #   here SMALL values are desirable (raters / items agree), and
#'   #   shares > 0.10 hint at meaningful systematic facet differences.
#'   # - `Interpretation` summarises the variance-share band the helper
#'   #   has assigned to each row.
#' }
#' }
#' @export
compute_facet_icc <- function(data, facets, score,
                              person = NULL, reml = TRUE,
                              ci_method = c("none", "profile", "boot"),
                              ci_level = 0.95,
                              ci_boot_reps = 1000L,
                              ci_boot_seed = NULL,
                              ci_boot_parallel = c("no", "multicore", "snow"),
                              ci_boot_ncpus = 1L) {
  ci_method <- match.arg(ci_method)
  ci_boot_parallel <- match.arg(ci_boot_parallel)
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number in (0, 1).", call. = FALSE)
  }
  ci_boot_reps <- max(1L, as.integer(ci_boot_reps))
  ci_boot_ncpus <- max(1L, as.integer(ci_boot_ncpus))
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
  missing_cols <- setdiff(c(score, facets, person), names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  # Coerce to factor so lmer treats them as grouping factors.
  for (col in c(facets, person)) {
    if (!is.null(col) && !is.factor(data[[col]])) {
      data[[col]] <- as.factor(as.character(data[[col]]))
    }
  }
  data[[score]] <- suppressWarnings(as.numeric(data[[score]]))
  data <- data[is.finite(data[[score]]), , drop = FALSE]
  if (nrow(data) == 0L) {
    stop("No finite rows in `score`.", call. = FALSE)
  }

  re_terms <- c(if (!is.null(person)) person else NULL, facets)
  formula <- stats::as.formula(paste0(
    score, " ~ 1 + ",
    paste0("(1 | ", re_terms, ")", collapse = " + ")
  ))
  # Capture lme4 convergence warnings so users know when the variance
  # decomposition rests on a near-singular fit. Silently refitting under
  # suppressWarnings() hid genuine problems (e.g. all-same scores).
  lmer_warnings <- character(0)
  fit <- tryCatch(
    withCallingHandlers(
      lme4::lmer(formula, data = data, REML = isTRUE(reml)),
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
            length(lmer_warnings), " convergence warning(s); ",
            "the ICC table is returned but results may be unreliable. ",
            "First message: ", lmer_warnings[1])
  }

  vc <- as.data.frame(lme4::VarCorr(fit))
  vc <- vc[is.na(vc$var2), c("grp", "vcov")]
  total_var <- sum(vc$vcov, na.rm = TRUE)
  # Treat "effectively zero" variance as non-identifiable; a positive
  # tolerance (sqrt(.Machine$double.eps)) prevents the "Variance 0 but
  # ICC 0.27" artifact that arises when lme4 returns boundary-singular
  # components on the order of 1e-30.
  zero_tol <- sqrt(.Machine$double.eps)
  is_singular <- is.finite(total_var) && total_var <= zero_tol
  icc_vec <- if (is.finite(total_var) && total_var > zero_tol) {
    vc$vcov / total_var
  } else {
    rep(NA_real_, length(vc$vcov))
  }
  if (isTRUE(is_singular)) {
    message("compute_facet_icc(): total variance is numerically zero ",
            "(non-identifiable); returning NA ICCs. ",
            "This usually means the score column has no within-facet spread.")
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
      return(if (isTRUE(is_singular)) "Non-identifiable" else NA_character_)
    }
    grp <- as.character(vc$grp[k])
    if (identical(grp, person_label)) koo_li_band(icc_vec[k])
    else var_share_band(icc_vec[k])
  }, character(1))

  out <- data.frame(
    Facet = vc$grp,
    Variance = round(vc$vcov, 6),
    ICC = round(icc_vec, 4),
    Interpretation = interpret,
    InterpretationScale = ifelse(
      as.character(vc$grp) == person_label,
      "Koo-Li reliability",
      "Variance share"
    ),
    stringsAsFactors = FALSE
  )

  # Optional ICC confidence intervals. Koo & Li (2016) recommend
  # applying the reliability bands to the 95% CI rather than the point
  # estimate; this block adds CI columns so callers can implement that
  # recommendation. Two methods are supported:
  #   * "profile": likelihood-profile bounds on the standard-deviation
  #     components via lme4::confint(method = "profile"), then
  #     transformed to ICC share (Variance_j / sum of squared bounds).
  #     Fast and deterministic, but may fail on singular fits.
  #   * "boot": parametric bootstrap via lme4::bootMer, draws
  #     `ci_boot_reps` simulated datasets from the fitted model and
  #     refits to build the empirical CI distribution. Slow, but
  #     robust to non-normality of the ICC sampling distribution.
  out$ICC_CI_Lower <- NA_real_
  out$ICC_CI_Upper <- NA_real_
  out$ICC_CI_Level <- ci_level
  out$ICC_CI_Method <- ci_method
  if (ci_method != "none" && !is_singular && is.finite(total_var)) {
    ci_result <- tryCatch(
      .compute_icc_ci(
        fit = fit, vc_grp = as.character(vc$grp),
        method = ci_method, ci_level = ci_level,
        boot_reps = ci_boot_reps, boot_seed = ci_boot_seed,
        boot_parallel = ci_boot_parallel,
        boot_ncpus = ci_boot_ncpus
      ),
      error = function(e) {
        message("compute_facet_icc(): CI computation (",
                ci_method, ") failed: ", conditionMessage(e),
                ". Returning point estimates only.")
        NULL
      }
    )
    if (!is.null(ci_result)) {
      out$ICC_CI_Lower <- round(ci_result$lower, 4)
      out$ICC_CI_Upper <- round(ci_result$upper, 4)
      if ("n_reps" %in% names(ci_result)) {
        out$ICC_CI_NReps <- ci_result$n_reps
      }
    }
  }
  structure(out, class = c("mfrm_facet_icc", "data.frame"))
}

#' @keywords internal
#' @noRd
# Map lme4::confint()-style row names to VarCorr component positions.
#
# lme4 uses two row-name conventions for the SD components returned
# by stats::confint(, method = "profile") / "Wald" / "boot":
#
#   (a) Terse form: ".sig01", ".sig02", ..., ".sigNN" for each random-
#       effect SD (in VarCorr formula order), plus ".sigma" for the
#       residual SD. This is the default return from
#       stats::confint(fit, method = "profile") on older / current
#       lme4 releases and for simple random-intercept-only models.
#
#   (b) Verbose form: "sd_(Intercept)|<group>" for each random-effect
#       SD plus "sigma" for the residual SD. Some lme4 branches and
#       broom.mixed post-processing surface this form.
#
# This helper returns an integer vector of length `length(vc_grp)`,
# where each entry is the row index in `row_names` corresponding to
# that VarCorr component's SD, or NA_integer_ when no row matches.
# Extracting it lets us add format-variant regression tests
# (see test-lme4-confint-helper.R) instead of relying on the inline
# grep pattern.
.lme4_confint_components <- function(row_names, vc_grp) {
  row_names <- as.character(row_names)
  vc_grp <- as.character(vc_grp)
  n_grp <- length(vc_grp)
  ordered <- rep(NA_integer_, n_grp)
  if (n_grp == 0L || length(row_names) == 0L) return(ordered)

  # Terse form first: ".sig01"..".sigNN" + ".sigma".
  terse_sig <- grep("^\\.sig[0-9]+$", row_names)
  terse_sigma <- which(row_names == ".sigma")
  if (length(terse_sig) > 0L || length(terse_sigma) > 0L) {
    re_positions <- which(vc_grp != "Residual")
    resid_position <- which(vc_grp == "Residual")
    for (k in seq_along(re_positions)) {
      if (k <= length(terse_sig)) ordered[re_positions[k]] <- terse_sig[k]
    }
    if (length(resid_position) == 1L && length(terse_sigma) == 1L) {
      ordered[resid_position] <- terse_sigma
    }
    return(ordered)
  }

  # Verbose form: "sd_<anything>|<group>" + "sigma". We use a literal
  # suffix match via endsWith() so group names containing regex
  # metacharacters (e.g. "Rater[1]") do not need manual escaping.
  for (i in seq_len(n_grp)) {
    grp <- vc_grp[i]
    if (identical(grp, "Residual")) {
      idx <- which(row_names == "sigma")
    } else {
      idx <- which(endsWith(row_names, paste0("|", grp)))
    }
    if (length(idx) == 1L) ordered[i] <- idx
  }
  ordered
}

# Internal helper: compute ICC confidence intervals by profile likelihood
# or parametric bootstrap, aligned with the component order returned by
# lme4::VarCorr.
.compute_icc_ci <- function(fit, vc_grp, method, ci_level,
                            boot_reps = 1000L, boot_seed = NULL,
                            boot_parallel = "no", boot_ncpus = 1L) {
  alpha <- 1 - ci_level
  n_grp <- length(vc_grp)
  if (identical(method, "profile")) {
    # Profile CIs for each SD component plus residual sigma. Row
    # labels depend on lme4 version (see `.lme4_confint_components`
    # for the supported formats).
    ci <- suppressWarnings(suppressMessages(
      stats::confint(fit, level = ci_level, method = "profile")
    ))
    ordered_rows <- .lme4_confint_components(rownames(ci), vc_grp)
    comp_sd_lo <- rep(NA_real_, n_grp)
    comp_sd_hi <- rep(NA_real_, n_grp)
    for (i in seq_len(n_grp)) {
      if (!is.finite(ordered_rows[i]) || ordered_rows[i] < 1L) next
      comp_sd_lo[i] <- ci[ordered_rows[i], 1L]
      comp_sd_hi[i] <- ci[ordered_rows[i], 2L]
    }
    # Approximate ICC CI by holding other components at their point
    # estimate while the focal component moves across its sd CI.
    vc_point <- as.data.frame(lme4::VarCorr(fit))
    vc_point <- vc_point[is.na(vc_point$var2), c("grp", "vcov")]
    point_var <- stats::setNames(vc_point$vcov, as.character(vc_point$grp))[vc_grp]
    lo <- rep(NA_real_, n_grp)
    hi <- rep(NA_real_, n_grp)
    for (i in seq_len(n_grp)) {
      if (!is.finite(comp_sd_lo[i]) || !is.finite(comp_sd_hi[i])) next
      var_lo <- comp_sd_lo[i]^2
      var_hi <- comp_sd_hi[i]^2
      others <- sum(point_var[-i], na.rm = TRUE)
      total_lo <- var_lo + others
      total_hi <- var_hi + others
      lo[i] <- if (total_lo > 0) var_lo / total_lo else NA_real_
      hi[i] <- if (total_hi > 0) var_hi / total_hi else NA_real_
    }
    return(list(lower = lo, upper = hi))
  }
  # Parametric bootstrap path.
  if (!is.null(boot_seed) && is.finite(boot_seed)) {
    set.seed(as.integer(boot_seed))
  }
  # Show an interactive progress bar for parametric bootstrap CIs because
  # `lme4::bootMer` can take many seconds even at default `nsim = 1000`.
  # cli::cli_progress_bar is silent in non-interactive contexts (e.g.
  # R CMD check) and respects options(cli.progress_show_after).
  # The bar is disabled when bootMer runs in parallel because worker
  # processes hold their own copy of `progress_id` / `reps_done`, so
  # in-process updates do not reach the parent.
  reps_done <- 0L
  total_reps <- as.integer(boot_reps)
  progress_id <- NULL
  use_progress <- total_reps > 1L && identical(boot_parallel, "no")
  if (use_progress) {
    progress_id <- cli::cli_progress_bar(
      name = "compute_facet_icc(boot)",
      total = total_reps,
      format = "{cli::pb_spin} bootstrap ICC: {cli::pb_current}/{cli::pb_total} [{cli::pb_elapsed}]",
      clear = TRUE,
      .envir = parent.frame()
    )
    on.exit(cli::cli_progress_done(id = progress_id), add = TRUE)
  }
  icc_of <- function(fit_b) {
    vc <- as.data.frame(lme4::VarCorr(fit_b))
    vc <- vc[is.na(vc$var2), c("grp", "vcov")]
    tot <- sum(vc$vcov, na.rm = TRUE)
    out <- if (!is.finite(tot) || tot <= 0) {
      rep(NA_real_, n_grp)
    } else {
      share <- stats::setNames(vc$vcov / tot, as.character(vc$grp))
      as.numeric(share[vc_grp])
    }
    if (!is.null(progress_id)) {
      reps_done <<- reps_done + 1L
      cli::cli_progress_update(id = progress_id, set = reps_done)
    }
    out
  }
  b <- suppressWarnings(suppressMessages(
    lme4::bootMer(fit, FUN = icc_of, nsim = boot_reps,
                  type = "parametric",
                  parallel = boot_parallel, ncpus = boot_ncpus,
                  use.u = FALSE)
  ))
  t_mat <- b$t
  lo <- apply(t_mat, 2L, stats::quantile, probs = alpha / 2,
              na.rm = TRUE, names = FALSE)
  hi <- apply(t_mat, 2L, stats::quantile, probs = 1 - alpha / 2,
              na.rm = TRUE, names = FALSE)
  list(lower = unname(lo), upper = unname(hi),
       n_reps = sum(stats::complete.cases(t_mat)))
}

#' Compute Kish design effects for each facet
#'
#' Combines per-facet average cluster size with ICC estimates to return
#' the Kish (1965) design effect `Deff = 1 + (m - 1) * rho`, where `m`
#' is the average number of observations per facet element and `rho` is
#' the ICC.
#'
#' @param data Data frame in long format.
#' @param facets Character vector of facet column names.
#' @param icc_table Output from [compute_facet_icc()] (optional; will be
#'   computed on the fly when `NULL`).
#' @param score Score column name; required when `icc_table` is `NULL`.
#' @param person Person column; passed through to compute_facet_icc().
#'
#' @section Interpreting output:
#' - `Deff = 1`: facet behaves like simple random sampling; no
#'   clustering-induced variance inflation.
#' - `Deff > 1`: variance of the mean estimate is inflated by a factor
#'   of `Deff` relative to SRS. `EffectiveN = N / Deff` is the sample
#'   size one would need under SRS to achieve the same precision. For
#'   rater-mediated designs, `Deff` well above 1 on the Rater facet
#'   means rater-level clustering is noticeable; consider whether
#'   rater generalisation is warranted.
#' - Reported `ICC` is pulled from `icc_table$ICC` (the variance share);
#'   interpretation is the same as in [compute_facet_icc()].
#'
#' @section Typical workflow:
#' 1. Run [compute_facet_icc()] to get the variance-component shares.
#' 2. Feed the result and the data into
#'    `compute_facet_design_effect(data, facets, icc_table = icc)`.
#' 3. Use `Deff` as part of the Methods discussion when generalising
#'    over raters or sites. Large `Deff` values argue for reporting
#'    robust SEs or moving to a hierarchical model.
#'
#' @return A data.frame of class `mfrm_facet_design_effect` with columns
#'   `Facet`, `AvgClusterSize`, `ICC`, `DesignEffect`, and `EffectiveN`.
#'
#' @seealso [compute_facet_icc()], [analyze_hierarchical_structure()].
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. New York: Wiley.
#'
#' Park, I., & Lee, H. (2001). The design effect: Do we know all about
#' it? In *Proceedings of the American Statistical Association, Survey
#' Research Methods Section* (pp. 143-148).
#'
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   icc <- compute_facet_icc(toy, facets = c("Rater", "Criterion"),
#'                            score = "Score", person = "Person")
#'   deff <- compute_facet_design_effect(toy,
#'                                       facets = c("Rater", "Criterion"),
#'                                       icc_table = icc)
#'   print(deff)
#'   # Large DesignEffect -> modest EffectiveN relative to raw N.
#' }
#' }
#' @export
compute_facet_design_effect <- function(data, facets, icc_table = NULL,
                                        score = NULL, person = NULL) {
  facets <- as.character(facets)
  if (is.null(icc_table)) {
    if (is.null(score)) {
      stop("Supply either `icc_table` or `score`.", call. = FALSE)
    }
    icc_table <- compute_facet_icc(data, facets = facets,
                                   score = score, person = person)
  }
  total_n <- nrow(data)
  out_rows <- lapply(facets, function(f) {
    if (!f %in% names(data)) {
      return(data.frame(
        Facet = f, AvgClusterSize = NA_real_, ICC = NA_real_,
        DesignEffect = NA_real_, EffectiveN = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    k <- length(unique(stats::na.omit(data[[f]])))
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
  structure(do.call(rbind, out_rows),
            class = c("mfrm_facet_design_effect", "data.frame"))
}


# ---- 4. analyze_hierarchical_structure ------------------------------------

#' Analyze the hierarchical structure of a rating design
#'
#' One-stop review that combines the nesting, cross-tabulation, ICC, and
#' design-effect reports into a single object. Designed to be reused by
#' the publication-workflow surface: its summary feeds into
#' `reporting_checklist()`, and its tables are picked up by
#' `build_mfrm_manifest()` for reproducibility bundles.
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
#'   only), `"profile"`, or `"boot"`. Deprecated alias:
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
#'   Deprecated spellings of the `ci_*` arguments above, retained for
#'   one release. Supplying a non-`NULL` value routes through
#'   [lifecycle::deprecate_warn()] and overrides the canonical
#'   `ci_*` argument.
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
#' - `design_effect`: Kish (1965) `Deff` and `EffectiveN`.
#' - `connectivity`: number of bipartite components linking
#'   Person x facet levels. A single component is required for a
#'   common measurement scale; multiple components indicate a
#'   disconnected design.
#'
#' @section Typical workflow:
#' 1. Optional: fit the MFRM with `fit_mfrm()`.
#' 2. Call `analyze_hierarchical_structure(fit)` (or on the raw data).
#' 3. Read `summary(x)` for the condensed view.
#' 4. Feed the object to [reporting_checklist()] and
#'    [build_mfrm_manifest()] to record the review in publication
#'    bundles. `build_apa_outputs()` uses the fit-level
#'    `FacetSampleSizeFlag` to add a Methods sentence automatically.
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
#' Winsteps.com. <https://www.winsteps.com/facets.htm>
#'
#' Kish, L. (1965). *Survey Sampling*. New York: Wiley.
#'
#' Koo, T. K., & Li, M. Y. (2016). A guideline of selecting and
#' reporting intraclass correlation coefficients for reliability
#' research. *Journal of Chiropractic Medicine, 15*(2), 155-163.
#'
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' hs <- analyze_hierarchical_structure(toy,
#'                                      facets = c("Rater", "Criterion"),
#'                                      compute_icc = FALSE,
#'                                      igraph_layout = FALSE)
#' summary(hs)
#'
#' \donttest{
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
                                           ci_method = c("none", "profile", "boot"),
                                           ci_level = 0.95,
                                           ci_boot_reps = 1000L,
                                           ci_boot_seed = NULL,
                                           igraph_layout = TRUE,
                                           icc_ci_method = NULL,
                                           icc_ci_level = NULL,
                                           icc_ci_boot_reps = NULL,
                                           icc_ci_boot_seed = NULL) {
  # Deprecated `icc_ci_*` spellings route through lifecycle and
  # override the canonical `ci_*` values when supplied. This unifies
  # the API with compute_facet_icc() while preserving one release of
  # backward compatibility.
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
  ci_method <- match.arg(ci_method)
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
    !is.null(score) && score %in% names(data)
  if (icc_available) {
    icc_tbl <- tryCatch(
      compute_facet_icc(data, facets = facets, score = score,
                        person = person,
                        ci_method = ci_method,
                        ci_level = ci_level,
                        ci_boot_reps = ci_boot_reps,
                        ci_boot_seed = ci_boot_seed),
      error = function(e) {
        message("ICC computation failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(icc_tbl) && nrow(icc_tbl) > 0) {
      deff_tbl <- tryCatch(
        compute_facet_design_effect(data, facets = facets,
                                    icc_table = icc_tbl,
                                    score = score, person = person),
        error = function(e) NULL
      )
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
      facets = facets
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

#' Plot a facet sample-size review
#'
#' Per-level observation counts rendered as a horizontal bar chart
#' coloured by the Linacre sample-size band assigned in
#' [facet_small_sample_review()]. Vertical dashed lines mark the
#' sparse / marginal / standard thresholds so reviewers see where
#' every facet level sits relative to the Linacre (1994) guidance.
#'
#' @param x An `mfrm_facet_sample_review` object.
#' @param top_n Optional integer; trim the y-axis to the `top_n`
#'   smallest level counts per facet. `NULL` (default) keeps all.
#' @param preset One of `"standard"`, `"publication"`, `"compact"`, `"monochrome"`.
#' @param ... Reserved.
#' @return Invisibly, the data.frame used for the plot.
#' @seealso [facet_small_sample_review()].
#' @export
plot.mfrm_facet_sample_review <- function(x, top_n = NULL,
                                          preset = c("standard",
                                                     "publication",
                                                     "compact",
                                                     "monochrome"),
                                          ...) {
  style <- resolve_plot_preset(preset)
  tbl <- as.data.frame(x$table, stringsAsFactors = FALSE)
  if (is.null(tbl) || nrow(tbl) == 0L) {
    graphics::plot.new()
    graphics::title(main = "Facet sample-size review")
    graphics::text(0.5, 0.5, "Review table empty.")
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

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4, max(8, min(18, max(nchar(labels)) * 0.55)), 3, 1))
  graphics::barplot(
    height = tbl$N,
    names.arg = labels,
    horiz = TRUE,
    las = 1,
    xlim = c(0, xmax),
    col = tbl$BarColor,
    border = NA,
    xlab = "Observations per level",
    main = "Facet sample-size review (Linacre bands)",
    cex.names = 0.8
  )
  graphics::abline(
    v = c(as.numeric(thr$sparse),
          as.numeric(thr$marginal),
          as.numeric(thr$standard)),
    lty = 2, col = style$neutral
  )
  graphics::legend(
    "bottomright", bty = "n", cex = 0.8,
    legend = c(sprintf("< %s sparse", thr$sparse),
               sprintf("< %s marginal", thr$marginal),
               sprintf("< %s standard", thr$standard),
               "strong"),
    fill = unname(band_colors)
  )
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
#' @export
plot.mfrm_facet_nesting <- function(x,
                                    preset = c("standard",
                                               "publication",
                                               "compact",
                                               "monochrome"),
                                    ...) {
  style <- resolve_plot_preset(preset)
  pair <- as.data.frame(x$pairwise_table, stringsAsFactors = FALSE)
  if (is.null(pair) || nrow(pair) == 0L) {
    graphics::plot.new()
    graphics::title(main = "Facet nesting (pairwise)")
    graphics::text(0.5, 0.5, "At least two facets required.")
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
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(5, 5, 3, 2))
  graphics::image(
    x = seq_len(n), y = seq_len(n),
    z = m,
    col = cols,
    zlim = c(0, 1),
    xaxt = "n", yaxt = "n",
    xlab = "Nested in (column facet)",
    ylab = "Nested facet (row)",
    main = "Pairwise nesting index"
  )
  graphics::axis(1, at = seq_len(n), labels = facets, las = 2, cex.axis = 0.8)
  graphics::axis(2, at = seq_len(n), labels = facets, las = 1, cex.axis = 0.8)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (is.finite(m[i, j])) {
      graphics::text(j, i, sprintf("%.2f", m[i, j]), cex = 0.75,
                     col = if (m[i, j] > 0.6) "white" else "black")
    }
  }
  invisible(m)
}

#' @export
print.mfrm_facet_sample_review <- function(x, ...) {
  cat("mfrm_facet_sample_review\n")
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

#' @export
print.mfrm_facet_icc <- function(x, ...) {
  cat("mfrm_facet_icc\n")
  if (nrow(x) == 0L) {
    cat("  (empty; lme4 unavailable or fit failed)\n")
  } else {
    print.data.frame(x, row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.mfrm_facet_icc <- function(object, ...) {
  # Condensed view that separates the two interpretation scales so
  # readers don't conflate person reliability with non-person variance
  # share; see `compute_facet_icc()` "Interpreting output".
  cat("Facet ICC summary (mfrmr)\n")
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
  cat("mfrm_facet_design_effect (Kish, 1965)\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#' @export
summary.mfrm_facet_design_effect <- function(object, ...) {
  cat("Kish design-effect summary (mfrmr)\n")
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
    cat("\nDesign effects (Kish):\n")
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
    graphics::image(
      x = seq_len(nrow(mat_num)),
      y = seq_len(ncol(mat_num)),
      z = mat_num,
      xaxt = "n", yaxt = "n",
      xlab = tbl$FacetA[1L],
      ylab = tbl$FacetB[1L],
      main = paste0("Cross-tabulation: ", pair_name),
      col = grDevices::hcl.colors(20, palette = "YlGnBu", rev = TRUE)
    )
    graphics::axis(1L, at = seq_len(nrow(mat_num)),
                   labels = rownames(mat_num), las = 2L, cex.axis = 0.7)
    graphics::axis(2L, at = seq_len(ncol(mat_num)),
                   labels = colnames(mat_num), las = 2L, cex.axis = 0.7)
  } else if (type == "icc") {
    if (is.null(x$icc) || nrow(x$icc) == 0L) {
      stop("No ICC table available; re-run with compute_icc = TRUE and lme4 installed.",
           call. = FALSE)
    }
    icc_tbl <- x$icc
    has_ci <- all(c("ICC_CI_Lower", "ICC_CI_Upper") %in% names(icc_tbl)) &&
      any(is.finite(icc_tbl$ICC_CI_Lower) | is.finite(icc_tbl$ICC_CI_Upper))
    ci_level <- if (has_ci && "ICC_CI_Level" %in% names(icc_tbl)) {
      suppressWarnings(as.numeric(icc_tbl$ICC_CI_Level[1L]))
    } else NA_real_
    ci_method <- if (has_ci && "ICC_CI_Method" %in% names(icc_tbl)) {
      as.character(icc_tbl$ICC_CI_Method[1L])
    } else NA_character_
    y_max <- max(
      1,
      max(c(icc_tbl$ICC, icc_tbl$ICC_CI_Upper), na.rm = TRUE) * 1.1
    )
    main_txt <- "Facet ICC (variance component share)"
    if (has_ci && is.finite(ci_level)) {
      main_txt <- sprintf("%s\n%g%% CI via %s",
                          main_txt, round(100 * ci_level), ci_method)
    }
    mids <- graphics::barplot(
      height = icc_tbl$ICC,
      names.arg = icc_tbl$Facet,
      main = main_txt,
      ylab = "ICC",
      ylim = c(0, y_max)
    )
    graphics::abline(h = c(0.5, 0.75, 0.9), lty = 2,
                     col = c("grey70", "grey50", "grey30"))
    if (has_ci) {
      valid <- is.finite(icc_tbl$ICC_CI_Lower) & is.finite(icc_tbl$ICC_CI_Upper)
      if (any(valid)) {
        graphics::arrows(
          x0 = mids[valid], y0 = icc_tbl$ICC_CI_Lower[valid],
          x1 = mids[valid], y1 = icc_tbl$ICC_CI_Upper[valid],
          angle = 90, code = 3, length = 0.05, col = "black", lwd = 1.5
        )
      }
    }
  }
  invisible(x)
}
