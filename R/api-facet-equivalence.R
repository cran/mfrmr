#' Analyze practical equivalence within a facet
#'
#' @param fit Output from [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()]. When `NULL`,
#'   diagnostics are computed with `residual_pca = "none"`.
#' @param facet Character scalar naming the non-person facet to evaluate. If
#'   `NULL`, the function prefers a rater-like facet and otherwise uses the
#'   first model facet.
#' @param equivalence_bound Practical-equivalence bound in logits. Default
#'   `0.5`.
#' @param conf_level Confidence level used for the forest-style interval view.
#'   Default `0.95`.
#'
#' @details
#' This function tests whether facet elements (e.g., raters) are similar
#' enough to be treated as practically interchangeable, rather than merely
#' testing whether they differ significantly.  This is the key distinction
#' from a standard chi-square heterogeneity test: absence of evidence
#' for difference is not evidence of equivalence.
#'
#' The function uses existing facet estimates and their standard errors
#' from `diagnostics$measures`; no re-estimation is performed.
#'
#' The bundle combines four complementary views:
#'
#' 1. **Fixed chi-square test**: tests \eqn{H_0}: all element measures
#'    are equal.  A non-significant result is *necessary but not
#'    sufficient* for interchangeability. It is reported as context, not
#'    as direct evidence of equivalence.
#'
#' 2. **Pairwise TOST (Two One-Sided Tests)**: for each pair of
#'    elements, tests whether the difference falls within
#'    \eqn{\pm}`equivalence_bound`.  The TOST procedure (Schuirmann,
#'    1987) rejects the null hypothesis of *non-equivalence* when both
#'    one-sided tests are significant at level \eqn{\alpha}.  A pair is
#'    declared "Equivalent" when the TOST p-value < 0.05.
#'
#' 3. **BIC-based Bayes-factor heuristic**: an approximate screening
#'    tool (not full Bayesian inference) that compares the evidence for
#'    a common-facet model (all elements equal) against a heterogeneity
#'    model (elements differ).  Values > 3 favour the common-facet
#'    model; < 1/3 favour heterogeneity.
#'
#' 4. **ROPE-style grand-mean proximity**: the proportion of each
#'    element's normal-approximation confidence distribution that falls
#'    within \eqn{\pm}`equivalence_bound` of the weighted grand mean.
#'    This is a descriptive proximity summary, not a Bayesian ROPE
#'    decision rule around a prespecified null value.
#'
#' **Choosing `equivalence_bound`**: the default of 0.5 logits is a
#' moderate criterion.  For high-stakes certification, 0.3 logits may
#' be appropriate; for exploratory or low-stakes contexts, 1.0 logits
#' may suffice.  The bound should reflect the smallest difference that
#' would be practically meaningful in your application.
#'
#' @section What this analysis means:
#' `analyze_facet_equivalence()` is a practical-interchangeability screen. It
#' asks whether facet levels are close enough, under a user-defined logit
#' bound, to be treated as practically similar for the current use case.
#'
#' @section What this analysis does not justify:
#' - A non-significant chi-square result is not evidence of equivalence.
#' - Forest/ROPE displays are descriptive and do not replace the pairwise TOST
#'   decision rule.
#' - The BIC-based Bayes-factor summary is a heuristic screen, not a full
#'   Bayesian equivalence analysis.
#'
#' @section Interpreting output:
#' Start with `summary$Decision`, which is a conservative summary of the
#' pairwise TOST results. Then use the remaining tables as context:
#' - `chi_square`: is there broad heterogeneity in the facet?
#' - `pairwise`: which specific pairs meet the practical-equivalence bound?
#' - `rope` / `forest`: how close is each level to the facet grand mean?
#'
#' Smaller `equivalence_bound` values make the criterion stricter. If the
#' decision is `"partial_pairwise_equivalence"`, that means some pairwise
#' contrasts satisfy the practical-equivalence bound but not all of them do.
#'
#' @section Decision rule:
#' The final `Decision` is a pairwise TOST summary rather than a global
#' equivalence proof. If all pairwise contrasts satisfy the practical-
#' equivalence bound, the facet is labeled `"all_pairs_equivalent"`. If at
#' least one, but not all, pairwise contrasts are equivalent, the facet is
#' labeled `"partial_pairwise_equivalence"`. If no pairwise contrasts meet the
#' practical-equivalence bound, the facet is labeled
#' `"no_pairwise_equivalence_established"`. The chi-square, Bayes-factor, and
#' grand-mean proximity summaries are reported as descriptive context.
#'
#' @section How to read the main outputs:
#' - `summary`: one-row pairwise-TOST decision summary and aggregate context.
#' - `pairwise`: pair-level TOST detail; use this for the primary inferential
#'   read.
#' - `chi_square`: broad heterogeneity screen.
#' - `rope` / `forest`: level-wise proximity to the weighted grand mean.
#'
#' @section Recommended next step:
#' If the result is borderline or high-stakes, re-run the analysis with a
#' tighter or looser `equivalence_bound`, then inspect `pairwise` and
#' [plot_facet_equivalence()] before deciding how strongly to claim
#' interchangeability.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#' 2. Run `analyze_facet_equivalence()` for the facet you want to screen.
#' 3. Read `summary` and `chi_square` first.
#' 4. Use [plot_facet_equivalence()] to inspect which levels drive the result.
#'
#' @section Output:
#' The returned bundle has class `mfrm_facet_equivalence` and includes:
#' - `summary`: one-row overview with convergent decision
#' - `chi_square`: fixed chi-square / separation summary
#' - `pairwise`: pairwise TOST detail table
#' - `rope`: element-wise ROPE probabilities around the weighted grand mean
#' - `forest`: element-wise estimate, confidence interval, and ROPE status
#' - `settings`: applied facet and threshold settings
#'
#' @return A named list with class `mfrm_facet_equivalence`.
#' @seealso [facets_chisq_table()], [fair_average_table()], [plot_facet_equivalence()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 25)
#' eq <- analyze_facet_equivalence(fit, facet = "Rater")
#' eq$summary[, c("Facet", "Elements", "Decision", "MeanROPE")]
#' head(eq$pairwise[, c("ElementA", "ElementB", "Equivalent")])
#' @export
analyze_facet_equivalence <- function(fit,
                                      diagnostics = NULL,
                                      facet = NULL,
                                      equivalence_bound = 0.5,
                                      conf_level = 0.95) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (!is.numeric(equivalence_bound) || length(equivalence_bound) != 1 || !is.finite(equivalence_bound) ||
      equivalence_bound <= 0) {
    stop("`equivalence_bound` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number between 0 and 1.", call. = FALSE)
  }

  if (is.null(diagnostics)) {
    diagnostics <- diagnose_mfrm(fit, residual_pca = "none")
  }
  if (!is.list(diagnostics) || is.null(diagnostics$measures)) {
    stop("`diagnostics` must include `measures`, typically from diagnose_mfrm().", call. = FALSE)
  }

  measures <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
  if (nrow(measures) == 0) {
    stop("`diagnostics$measures` is empty.", call. = FALSE)
  }
  if (!all(c("Facet", "Estimate", "SE") %in% names(measures))) {
    stop("`diagnostics$measures` must contain Facet, Estimate, and SE columns.", call. = FALSE)
  }

  facet_names <- as.character(fit$config$facet_names %||% character(0))
  if (length(facet_names) == 0 && "Facet" %in% names(measures)) {
    facet_names <- unique(as.character(measures$Facet))
  }
  if (length(facet_names) == 0) {
    stop("No non-person facets are available for equivalence analysis.", call. = FALSE)
  }

  if (is.null(facet) || !nzchar(as.character(facet[1]))) {
    rater_guess <- infer_default_rater_facet(facet_names)
    facet <- if (!is.null(rater_guess) && nzchar(rater_guess)) rater_guess else facet_names[1]
  } else {
    facet <- as.character(facet[1])
  }
  if (!facet %in% facet_names) {
    stop("`facet` must be one of: ", paste(facet_names, collapse = ", "), ".", call. = FALSE)
  }

  facet_df <- measures[as.character(measures$Facet) == facet, , drop = FALSE]
  if (nrow(facet_df) < 2) {
    stop("Facet '", facet, "' has fewer than 2 estimated levels.", call. = FALSE)
  }

  label_col <- if ("Level" %in% names(facet_df)) {
    "Level"
  } else if ("Element" %in% names(facet_df)) {
    "Element"
  } else {
    NULL
  }
  if (is.null(label_col)) {
    facet_df$Level <- paste0(facet, "_", seq_len(nrow(facet_df)))
    label_col <- "Level"
  }

  facet_df$Estimate <- suppressWarnings(as.numeric(facet_df$Estimate))
  facet_df$SE <- suppressWarnings(as.numeric(facet_df$SE))
  facet_df <- facet_df[is.finite(facet_df$Estimate) & is.finite(facet_df$SE) & facet_df$SE > 0, , drop = FALSE]
  if (nrow(facet_df) < 2) {
    stop("Facet '", facet, "' does not have at least 2 levels with finite Estimate and SE values.", call. = FALSE)
  }

  labels <- as.character(facet_df[[label_col]])
  est <- as.numeric(facet_df$Estimate)
  se <- as.numeric(facet_df$SE)
  weights <- 1 / (se ^ 2)
  grand_mean <- stats::weighted.mean(est, w = weights)

  n_elem <- length(est)
  df_chi <- n_elem - 1
  chi2_val <- sum(weights * (est - grand_mean) ^ 2)
  p_chi <- if (df_chi > 0) stats::pchisq(chi2_val, df = df_chi, lower.tail = FALSE) else NA_real_
  sep_sd <- if (n_elem > 1) stats::sd(est) else NA_real_
  rmse <- sqrt(mean(se ^ 2))
  true_sd <- sqrt(max((sep_sd %||% NA_real_) ^ 2 - rmse ^ 2, 0))
  separation <- if (is.finite(rmse) && rmse > 0) true_sd / rmse else NA_real_
  reliability <- if (is.finite(separation)) separation ^ 2 / (1 + separation ^ 2) else NA_real_

  n_obs <- suppressWarnings(as.numeric(fit$prep$n_obs %||% NA_real_))
  bic_diff <- if (is.finite(n_obs) && n_obs > 1 && is.finite(chi2_val) && is.finite(df_chi)) {
    chi2_val - df_chi * log(n_obs)
  } else {
    NA_real_
  }
  bf01 <- if (is.finite(bic_diff)) exp(max(min(-bic_diff / 2, 700), -700)) else NA_real_
  bf_label <- classify_equivalence_bf(bf01)

  z_ci <- stats::qnorm(1 - (1 - conf_level) / 2)
  z_tost <- stats::qnorm(0.95)
  pair_idx <- utils::combn(seq_len(n_elem), 2)
  pairwise_tbl <- if (is.null(dim(pair_idx))) {
    data.frame()
  } else {
    pairwise_rows <- lapply(seq_len(ncol(pair_idx)), function(k) {
      i <- pair_idx[1, k]
      j <- pair_idx[2, k]
      diff <- est[i] - est[j]
      se_diff <- sqrt(se[i] ^ 2 + se[j] ^ 2)
      if (!is.finite(se_diff) || se_diff <= 0) return(NULL)
      z_lower <- (diff + equivalence_bound) / se_diff
      z_upper <- (diff - equivalence_bound) / se_diff
      p_lower <- stats::pnorm(z_lower, lower.tail = FALSE)
      p_upper <- stats::pnorm(z_upper, lower.tail = TRUE)
      p_tost <- max(p_lower, p_upper)
      data.frame(
        ElementA = labels[i],
        ElementB = labels[j],
        Diff = diff,
        SE_Diff = se_diff,
        CI90_Lower = diff - z_tost * se_diff,
        CI90_Upper = diff + z_tost * se_diff,
        P_Lower = p_lower,
        P_Upper = p_upper,
        P_TOST = p_tost,
        Equivalent = is.finite(p_tost) & p_tost < 0.05,
        stringsAsFactors = FALSE
      )
    })
    pairwise_tbl <- do.call(rbind, pairwise_rows)
    if (is.null(pairwise_tbl)) pairwise_tbl <- data.frame()
    pairwise_tbl
  }
  n_pairs <- nrow(pairwise_tbl)
  n_equiv <- if (n_pairs > 0) sum(pairwise_tbl$Equivalent, na.rm = TRUE) else 0L

  rope_tbl <- data.frame(
    Element = labels,
    Measure = est,
    Deviation = est - grand_mean,
    SE = se,
    CI_Lower = est - z_ci * se,
    CI_Upper = est + z_ci * se,
    stringsAsFactors = FALSE
  )
  rope_tbl$ROPEPct <- 100 * (
    stats::pnorm(equivalence_bound, mean = rope_tbl$Deviation, sd = rope_tbl$SE) -
      stats::pnorm(-equivalence_bound, mean = rope_tbl$Deviation, sd = rope_tbl$SE)
  )
  lower_bound <- grand_mean - equivalence_bound
  upper_bound <- grand_mean + equivalence_bound
  rope_tbl$ROPEStatus <- ifelse(
    rope_tbl$CI_Lower >= lower_bound & rope_tbl$CI_Upper <= upper_bound,
    "inside",
    ifelse(
      rope_tbl$CI_Lower > upper_bound | rope_tbl$CI_Upper < lower_bound,
      "outside",
      "overlap"
    )
  )

  mean_rope <- if (nrow(rope_tbl) > 0) mean(rope_tbl$ROPEPct, na.rm = TRUE) else NA_real_
  all_pairs_equivalent <- n_pairs > 0L && n_equiv == n_pairs
  any_pair_equivalent <- n_pairs > 0L && n_equiv > 0L
  decision <- if (all_pairs_equivalent) {
    "all_pairs_equivalent"
  } else if (any_pair_equivalent) {
    "partial_pairwise_equivalence"
  } else {
    "no_pairwise_equivalence_established"
  }

  chi_square_tbl <- data.frame(
    Facet = facet,
    Elements = n_elem,
    GrandMean = grand_mean,
    FixedChiSq = chi2_val,
    FixedDF = df_chi,
    FixedProb = p_chi,
    Separation = separation,
    Reliability = reliability,
    stringsAsFactors = FALSE
  )

  summary_tbl <- data.frame(
    Facet = facet,
    Elements = n_elem,
    EquivalenceBound = equivalence_bound,
    GrandMean = grand_mean,
    FixedProb = p_chi,
    PairwiseComparisons = n_pairs,
    PairwiseEquivalent = n_equiv,
    PairwiseEquivalentPct = if (n_pairs > 0) 100 * n_equiv / n_pairs else NA_real_,
    BF01 = bf01,
    BF01Label = bf_label,
    MeanROPE = mean_rope,
    AllPairsEquivalent = all_pairs_equivalent,
    AnyPairEquivalent = any_pair_equivalent,
    PairwiseDecisionBasis = "pairwise_tost_summary",
    Decision = decision,
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary_tbl,
    chi_square = chi_square_tbl,
    pairwise = as.data.frame(pairwise_tbl, stringsAsFactors = FALSE),
    rope = as.data.frame(rope_tbl, stringsAsFactors = FALSE),
    forest = as.data.frame(rope_tbl, stringsAsFactors = FALSE),
    settings = list(
      facet = facet,
      equivalence_bound = equivalence_bound,
      conf_level = conf_level
    )
  )
  as_mfrm_bundle(out, "mfrm_facet_equivalence")
}

classify_equivalence_bf <- function(bf01) {
  if (!is.finite(bf01)) return("Not available")
  if (bf01 > 100) return("Extreme evidence for common-facet model")
  if (bf01 > 30) return("Very strong evidence for common-facet model")
  if (bf01 > 10) return("Strong evidence for common-facet model")
  if (bf01 > 3) return("Moderate evidence for common-facet model")
  if (bf01 > 1) return("Anecdotal evidence for common-facet model")
  if (bf01 > (1 / 3)) return("Anecdotal evidence for heterogeneity")
  if (bf01 > (1 / 10)) return("Moderate evidence for heterogeneity")
  "Strong evidence for heterogeneity"
}

#' Plot facet-equivalence results
#'
#' @param x Output from [analyze_facet_equivalence()] or [fit_mfrm()].
#' @param diagnostics Optional output from [diagnose_mfrm()] when `x` is an
#'   `mfrm_fit` object.
#' @param facet Facet to analyze when `x` is an `mfrm_fit` object.
#' @param type Plot type: `"forest"` (default) or `"rope"`.
#' @param draw If `TRUE` (default), draw the plot. If `FALSE`, return the
#'   prepared plotting data.
#' @param ... Additional graphical arguments passed to base plotting functions.
#'
#' @details
#' `plot_facet_equivalence()` is a visual companion to
#' [analyze_facet_equivalence()]. It does not recompute the equivalence
#' analysis; it only reshapes and displays the returned results.
#'
#' @section Plot types:
#' - `"forest"` places each level on the logit scale with its confidence
#'   interval and shades the practical-equivalence region around the weighted
#'   grand mean.
#' - `"rope"` shows the percentage of each level's uncertainty mass that falls
#'   inside the ROPE.
#'
#' @section Interpreting output:
#' In the **forest plot**, the shaded band marks the ROPE
#' (\eqn{\pm}`equivalence_bound` around the weighted grand mean).
#' Levels whose entire confidence interval lies inside this band are
#' close to the facet grand mean under this descriptive screen. Levels whose
#' interval extends outside the band are more displaced from the facet average.
#' Overlapping intervals between two elements suggest they are not
#' reliably separable, but overlap alone does not establish formal
#' equivalence---use the TOST results for that.
#'
#' In the **ROPE bar chart**, each bar shows the proportion of the
#' element's normal-approximation distribution that falls inside the
#' ROPE-style grand-mean proximity.  Values > 95\% indicate that most of
#' the element's normal-approximation uncertainty falls near the facet
#' average; 50--95\% is indeterminate; < 50\% suggests the element is
#' meaningfully displaced from that average.
#'
#' @section Typical workflow:
#' 1. Run [analyze_facet_equivalence()].
#' 2. Start with `type = "forest"` to see the facet on the logit scale.
#' 3. Switch to `type = "rope"` when you want a ranking of levels by
#'    grand-mean proximity.
#'
#' @return Invisibly returns the plotting data. If `draw = FALSE`, the plotting
#'   data are returned without drawing.
#' @seealso [analyze_facet_equivalence()]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 25)
#' eq <- analyze_facet_equivalence(fit, facet = "Rater")
#' pdat <- plot_facet_equivalence(eq, type = "forest", draw = FALSE)
#' c(pdat$facet, pdat$type)
#' @export
plot_facet_equivalence <- function(x,
                                   diagnostics = NULL,
                                   facet = NULL,
                                   type = c("forest", "rope"),
                                   draw = TRUE,
                                   ...) {
  type <- match.arg(type)
  if (inherits(x, "mfrm_fit")) {
    x <- analyze_facet_equivalence(
      fit = x,
      diagnostics = diagnostics,
      facet = facet
    )
  }
  if (!inherits(x, "mfrm_facet_equivalence")) {
    stop("`x` must be output from analyze_facet_equivalence() or fit_mfrm().", call. = FALSE)
  }

  forest_df <- as.data.frame(x$forest %||% data.frame(), stringsAsFactors = FALSE)
  settings <- x$settings %||% list()
  summary_tbl <- as.data.frame(x$summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(forest_df) == 0) {
    stop("No equivalence rows are available for plotting.", call. = FALSE)
  }

  out <- list(
    data = forest_df,
    facet = as.character(settings$facet %||% summary_tbl$Facet[1] %||% ""),
    grand_mean = suppressWarnings(as.numeric(summary_tbl$GrandMean[1] %||% NA_real_)),
    equivalence_bound = suppressWarnings(as.numeric(settings$equivalence_bound %||% summary_tbl$EquivalenceBound[1] %||% NA_real_)),
    type = type
  )
  if (!isTRUE(draw)) {
    return(out)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  dots <- list(...)

  if (identical(type, "forest")) {
    ord <- order(forest_df$Measure, decreasing = FALSE, na.last = TRUE)
    forest_df <- forest_df[ord, , drop = FALSE]
    ypos <- seq_len(nrow(forest_df))
    cols <- ifelse(
      forest_df$ROPEStatus == "inside",
      "#2E8B57",
      ifelse(forest_df$ROPEStatus == "outside", "#C0392B", "#D68910")
    )
    xlim <- range(c(forest_df$CI_Lower, forest_df$CI_Upper,
                    out$grand_mean - out$equivalence_bound,
                    out$grand_mean + out$equivalence_bound), finite = TRUE)
    do.call(graphics::plot, c(list(
      x = forest_df$Measure,
      y = ypos,
      xlim = xlim,
      yaxt = "n",
      ylab = "",
      xlab = "Measure (logits)",
      main = paste0(out$facet, ": facet equivalence"),
      pch = 19,
      col = cols
    ), dots))
    graphics::axis(2, at = ypos, labels = forest_df$Element, las = 2)
    graphics::rect(
      xleft = out$grand_mean - out$equivalence_bound,
      ybottom = 0.5,
      xright = out$grand_mean + out$equivalence_bound,
      ytop = nrow(forest_df) + 0.5,
      border = NA,
      col = grDevices::adjustcolor("#2E8B57", alpha.f = 0.12)
    )
    graphics::abline(v = out$grand_mean, lty = 2, col = "gray40")
    graphics::segments(forest_df$CI_Lower, ypos, forest_df$CI_Upper, ypos, col = cols, lwd = 2)
    graphics::points(forest_df$Measure, ypos, pch = 19, col = cols)
  } else {
    ord <- order(forest_df$ROPEPct, decreasing = TRUE, na.last = TRUE)
    forest_df <- forest_df[ord, , drop = FALSE]
    cols <- ifelse(
      forest_df$ROPEPct >= 95,
      "#2E8B57",
      ifelse(forest_df$ROPEPct < 50, "#C0392B", "#D68910")
    )
    mids <- graphics::barplot(
      height = forest_df$ROPEPct,
      names.arg = forest_df$Element,
      las = 2,
      ylim = c(0, 100),
      col = cols,
      ylab = "% in ROPE",
      main = paste0(out$facet, ": ROPE probabilities"),
      ...
    )
    graphics::abline(h = c(50, 80, 95), lty = c(3, 3, 2), col = c("gray60", "gray60", "gray40"))
    out$bar_midpoints <- mids
  }

  invisible(out)
}

#' @export
plot.mfrm_facet_equivalence <- function(x, ...) {
  plot_facet_equivalence(x, ...)
}
