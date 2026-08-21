# ==============================================================================
# mfrmr standardized/aggregated-residual Q3-style screen
# ==============================================================================
#
# `q3_statistic()` computes a Q3-style index of local response dependence
# between facet-level pairs from the residuals stored on a diagnostics
# bundle. Q3 was originally defined by Yen (1984, eqs. 7-8, p. 127) as the
# Pearson correlation of *raw* residuals (d_{ik} = u_{ik} - P_hat). mfrmr
# uses *standardized* residuals because that is what `diagnose_mfrm()`
# stores; this is a deliberate departure from Yen's original definition
# and is documented in the help page.
#
# Threshold attribution (descriptive convention vs primary source):
# - Yen (1984) did not propose a fixed 0.20 cutoff. The package's legacy rule
#   applies strict `|Q3-style| > threshold`; it is an uncalibrated screening
#   heuristic for this standardized/aggregated-residual index.
# - The |Q3| > 0.30 "strict" cutoff is a community convention summarized
#   in Marais (2013, p. 121) as "often considered" - not Marais's own
#   proposal. Marais's actual recommendation is to interpret Q3 *relative*
#   to the average pairwise correlation in the same dataset.
# - Christensen, Makransky, & Horton (2017) demonstrate empirically that
#   *no single critical value* is appropriate across designs: percentiles
#   of the Q3 null distribution depend on number of items, sample size, and
#   number of categories. They recommend a parametric bootstrap. Their
#   `Q3_*` statistic compares the *maximum* Q3 to the average, not every
#   pair to the average. mfrmr's `relative_offset` flag is a screening
#   approximation, not a re-implementation of Q3_*.
#
# This helper exposes the same numerical surface that
# `plot_local_dependence_heatmap()` draws, but returns a tidy
# data.frame with thresholds and interpretation labels for use in
# manuscript tables.

#' mfrmr standardized/aggregated-residual Q3-style screen
#'
#' Computes mfrmr's Q3-style screening index: the Pearson correlation of
#' **standardized residuals**, after mean-aggregation to one residual per
#' Person-by-level cell, between every pair of levels of a chosen facet.
#' Large absolute values identify pairs whose residuals show stronger positive
#' co-movement or negative counter-movement than the main-effects model
#' expects. This is not the raw-residual Yen (1984) Q3 statistic.
#'
#' @section Statistic definition and distinction from Yen (1984):
#' The name **Q3-style** is deliberate. This implementation differs from
#' Yen's (1984) original Q3 definition in two respects that together affect
#' threshold interpretation.
#'
#' **(1) Standardized vs raw residuals.** Yen (1984, eqs. 7-8, p. 127)
#' defines `Q3 = cor(d_i, d_j)` where `d_{ik} = u_{ik} - P_hat_{ik}` is
#' the **raw** residual. mfrmr uses **standardized** residuals
#' `Z = (u - P_hat) / sqrt(Var(u))` because that is what
#' `diagnose_mfrm()` stores. Because model-based residual scale factors can
#' vary across observations, the resulting pair correlations and their null
#' distribution differ from raw-residual Q3.
#'
#' **(2) Mean-aggregation.** When the facet being paired (e.g. `Rater`)
#' has multiple residual rows per (Person, Level) cell because of
#' additional facets in the design (e.g. multiple `Criterion` rows per
#' Person-Rater cell), the standardized residuals are first
#' **mean-aggregated to one value per (Person, Level) cell**, and the
#' Pearson correlation is taken over those mean-aggregated residuals.
#' Yen's original formulation takes the correlation directly over
#' per-(Person, Item) residuals, without aggregation. Aggregation changes the
#' analysis unit to persons. Each correlation is based on persons observed at
#' both levels, so its magnitude and null distribution are not directly
#' comparable with raw-row Q3.
#'
#' For both reasons, treat the values returned here as a **screening
#' summary** rather than a direct substitute for the published Q3
#' thresholds. A formal raw-residual Q3 procedure would require a separately
#' implemented and validated design-specific bootstrap; mfrmr 0.2.3 does not
#' provide that procedure.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output. Computed
#'   on demand when omitted.
#' @param facet Facet whose levels are paired (default `"Rater"`).
#' @param min_pairs Minimum number of persons with finite aggregated residuals
#'   at both levels required to retain a pair; a single integer of at least
#'   three. Pairs below the threshold drop out
#'   of the table (mirrors [plot_local_dependence_heatmap()]).
#' @param yen_threshold Legacy-named argument retained for compatibility. It
#'   applies the strict absolute heuristic `|Q3-style| > threshold` (default
#'   `0.20`). Yen (1984) did not propose this as a general cutoff, and it is not
#'   calibrated for mfrmr's standardized/aggregated-residual index.
#' @param marais_threshold Stricter community-convention threshold
#'   (default `0.30`). Marais (2013, p. 121) reports this as a value
#'   "often considered" in the literature, not as her own
#'   recommendation; her actual recommendation is the relative
#'   comparison implemented by `relative_offset`.
#' @param relative_offset Screening offset for the relative-flag rule
#'   `|Q3 - mean(Q3)| > relative_offset` (default `0.20`). This is a
#'   simplified screening approximation of the relative comparison
#'   advocated by Marais (2013) and operationalized by Christensen
#'   et al. (2017) as `Q3_* = Q3_max - mean(Q3)`. Christensen et al.
#'   (2017) demonstrate empirically that no single critical value is
#'   appropriate across designs and recommend a parametric bootstrap;
#'   the fixed `0.20` here is a screening default, not a substitute
#'   for that bootstrap.
#'
#' @return An object of class `mfrm_q3` containing:
#' \describe{
#'   \item{`pairs`}{A data frame with one row per facet-level pair
#'     and columns `Level1`, `Level2`, `Q3`, `N`, `AbsQ3`,
#'     `YenFlag`, `MaraisFlag`, `RelativeFlag`, and a textual
#'     `Interpretation` summarising which thresholds were exceeded. `N` is the
#'     number of persons with finite aggregated residuals at both levels.
#'     `YenFlag` is a compatibility name and does not imply that Yen (1984)
#'     proposed the fixed default.}
#'   \item{`summary`}{One-row tibble with `MeanQ3`, `MaxAbsQ3`,
#'     and the three flagged-pair counts.}
#'   \item{`thresholds`}{The thresholds used, for reproducibility. The `yen`
#'     name is retained for compatibility.}
#'   \item{`facet`}{The facet whose levels were paired.}
#' }
#'
#' @section References:
#' - Yen, W. M. (1984). Effects of local item dependence on the fit
#'   and equating performance of the three-parameter logistic model.
#'   *Applied Psychological Measurement, 8*(2), 125-145.
#'   \doi{10.1177/014662168400800201}
#' - Chen, W.-H., & Thissen, D. (1997). Local dependence indexes for
#'   item pairs using item response theory. *Journal of Educational
#'   and Behavioral Statistics, 22*(3), 265-289.
#'   \doi{10.3102/10769986022003265}
#' - Marais, I. (2013). Local dependence. In K. B. Christensen,
#'   S. Kreiner, & M. Mesbah (Eds.), *Rasch models in health*
#'   (pp. 111-130). London: ISTE / Wiley.
#' - Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical
#'   values for Yen's Q3: Identification of local dependence in the
#'   Rasch model using residual correlations.
#'   *Applied Psychological Measurement, 41*(3), 178-194.
#'   \doi{10.1177/0146621616677520}
#' - Nason, K., & DeMars, C. (2025). Another look at Yen's Q3: Is .2 an
#'   appropriate cut-off? *Journal of Educational Measurement, 62*(2),
#'   345-359. \doi{10.1111/jedm.12432}
#'
#' @seealso [plot_local_dependence_heatmap()], [diagnose_mfrm()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' q3 <- q3_statistic(fit)
#' q3$summary
#' # Look for: use MaxAbsQ3 and the pair table to rank follow-up. The default
#' #   0.20/0.30 absolute rules are uncalibrated heuristics for this Q3-style
#' #   index, not standalone local-independence tests.
#' head(q3$pairs)
#' }
#' @export
q3_statistic <- function(fit,
                         diagnostics = NULL,
                         facet = "Rater",
                         min_pairs = 5L,
                         yen_threshold = 0.20,
                         marais_threshold = 0.30,
                         relative_offset = 0.20) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  facet <- as.character(facet[1])
  facet_names <- as.character(fit$config$facet_names %||% character(0))
  if (!facet %in% facet_names) {
    stop("`facet` must be one of: ", paste(facet_names, collapse = ", "), ".",
         call. = FALSE)
  }
  if (!is.numeric(min_pairs) || length(min_pairs) != 1L ||
      !is.finite(min_pairs) || min_pairs < 3 || min_pairs != floor(min_pairs)) {
    stop("`min_pairs` must be a single integer of at least 3.", call. = FALSE)
  }
  min_pairs <- as.integer(min_pairs)
  valid_positive_scalar <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  }
  if (!valid_positive_scalar(yen_threshold) ||
      !valid_positive_scalar(marais_threshold) ||
      !valid_positive_scalar(relative_offset)) {
    stop("`yen_threshold`, `marais_threshold`, and `relative_offset` ",
         "must be finite positive numeric scalars.", call. = FALSE)
  }
  yen_threshold <- as.numeric(yen_threshold)
  marais_threshold <- as.numeric(marais_threshold)
  relative_offset <- as.numeric(relative_offset)
  if (is.null(diagnostics)) {
    diagnostics <- suppressMessages(
      diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "legacy")
    )
  }

  # Reuse the heatmap helper's plot data; it already implements the
  # standardized-residual pivot + pairwise Pearson correlation that
  # Q3 is defined as.
  h <- plot_local_dependence_heatmap(
    fit = fit,
    diagnostics = diagnostics,
    facet = facet,
    min_pairs = min_pairs,
    draw = FALSE
  )

  pairs_df <- as.data.frame(h$data$pairs %||% data.frame(),
                             stringsAsFactors = FALSE)
  if (nrow(pairs_df) == 0L) {
    return(structure(
      list(
        pairs = data.frame(Level1 = character(0), Level2 = character(0),
                            Q3 = numeric(0), N = integer(0),
                            AbsQ3 = numeric(0), YenFlag = logical(0),
                            MaraisFlag = logical(0),
                            RelativeFlag = logical(0),
                            Interpretation = character(0),
                            stringsAsFactors = FALSE),
        summary = data.frame(MeanQ3 = NA_real_, MaxAbsQ3 = NA_real_,
                              YenFlagged = 0L, MaraisFlagged = 0L,
                              RelativeFlagged = 0L,
                              stringsAsFactors = FALSE),
        thresholds = c(yen = yen_threshold,
                        marais = marais_threshold,
                        relative_offset = relative_offset),
        facet = facet
      ),
      class = c("mfrm_q3", "list")
    ))
  }

  pairs_df$Q3 <- suppressWarnings(as.numeric(pairs_df$ResidualCor))
  pairs_df$AbsQ3 <- abs(pairs_df$Q3)
  mean_q3 <- mean(pairs_df$Q3, na.rm = TRUE)

  pairs_df$YenFlag <- is.finite(pairs_df$AbsQ3) &
                       pairs_df$AbsQ3 > yen_threshold
  pairs_df$MaraisFlag <- is.finite(pairs_df$AbsQ3) &
                         pairs_df$AbsQ3 > marais_threshold
  pairs_df$RelativeFlag <- is.finite(pairs_df$Q3) &
                           abs(pairs_df$Q3 - mean_q3) > relative_offset

  pairs_df$Interpretation <- vapply(seq_len(nrow(pairs_df)), function(i) {
    flags <- character(0)
    if (isTRUE(pairs_df$MaraisFlag[i])) {
      flags <- c(flags, sprintf(
        "stricter heuristic (|Q3-style| > %.2f)",
        marais_threshold
      ))
    }
    if (isTRUE(pairs_df$YenFlag[i]) && !isTRUE(pairs_df$MaraisFlag[i])) {
      flags <- c(flags, sprintf(
        "fixed absolute rule (|Q3-style| > %.2f)",
        yen_threshold
      ))
    }
    if (isTRUE(pairs_df$RelativeFlag[i])) {
      flags <- c(flags, sprintf(
        "relative-to-mean (offset > %.2f)",
        relative_offset
      ))
    }
    if (length(flags) == 0L) "OK" else paste(flags, collapse = ", ")
  }, character(1))

  pairs_df <- pairs_df[order(pairs_df$AbsQ3, decreasing = TRUE), , drop = FALSE]

  summary_df <- data.frame(
    MeanQ3 = mean_q3,
    MaxAbsQ3 = max(pairs_df$AbsQ3, na.rm = TRUE),
    YenFlagged = sum(pairs_df$YenFlag, na.rm = TRUE),
    MaraisFlagged = sum(pairs_df$MaraisFlag, na.rm = TRUE),
    RelativeFlagged = sum(pairs_df$RelativeFlag, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  out <- list(
    pairs = pairs_df[, c("Level1", "Level2", "Q3", "N", "AbsQ3",
                          "YenFlag", "MaraisFlag", "RelativeFlag",
                          "Interpretation"), drop = FALSE],
    summary = summary_df,
    thresholds = c(yen = yen_threshold,
                   marais = marais_threshold,
                   relative_offset = relative_offset),
    facet = facet
  )
  class(out) <- c("mfrm_q3", "list")
  out
}

#' @export
print.mfrm_q3 <- function(x, ...) {
  cat("mfrmr standardized/aggregated-residual Q3-style screen\n")
  cat(sprintf("  Facet: %s\n", x$facet))
  cat(sprintf("  Mean Q3-style correlation: %.3f | Max absolute: %.3f\n",
              x$summary$MeanQ3, x$summary$MaxAbsQ3))
  cat(sprintf(paste0(
    "  Flagged pairs: fixed absolute rule %d / stricter absolute rule %d / ",
    "relative-to-mean rule %d (of %d)\n"
  ),
              x$summary$YenFlagged, x$summary$MaraisFlagged,
              x$summary$RelativeFlagged, nrow(x$pairs)))
  cat(sprintf(paste0(
    "  Heuristic thresholds: absolute %.2f / stricter %.2f / ",
    "relative offset %.2f\n"
  ),
              as.numeric(x$thresholds["yen"]),
              as.numeric(x$thresholds["marais"]),
              as.numeric(x$thresholds["relative_offset"])))
  cat("  Fixed absolute rules are uncalibrated heuristics for this index.\n")
  cat("  Interpretation: screening evidence, not a standalone local-independence test.\n")
  invisible(x)
}
