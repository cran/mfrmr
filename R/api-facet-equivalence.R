#' Analyze practical equivalence within a facet
#'
#' @param fit Output from [fit_mfrm()]. Requires an inference-ready MML fit
#'   with unregularized observed-information covariance and estimable contrasts.
#' @param diagnostics Optional matching output from [diagnose_mfrm()]. Supplied
#'   diagnostics must retain all target-facet levels, model-based SEs and
#'   ordinary-inference eligibility. Estimates and covariance are always taken
#'   from `fit`; supplying diagnostics cannot override its restrictions.
#'   Fixed levels instead require their fixed-value labels and absent sampling SEs.
#' @param facet Character scalar naming a non-person facet. When `NULL`, a
#'   rater-like facet is preferred, otherwise the first model facet is used.
#' @param equivalence_bound Positive practical-equivalence bound in logits.
#'   The default `0.5` is not a universal threshold. Choose the smallest
#'   practically meaningful difference for the intended use before inspecting
#'   the equivalence results.
#' @param ci_level Confidence level for level and grand-mean-deviation intervals
#'   (default `0.95`). Pairwise TOST always uses alpha 0.05 and 90% intervals.
#' @param conf_level Deprecated alias for `ci_level`; when supplied it takes
#'   precedence and emits a lifecycle deprecation warning.
#'
#' @details
#' Pair differences use the full constrained covariance from the MML observed
#' information: \eqn{\mathrm{Var}(A-B) = \mathrm{Var}(A) + \mathrm{Var}(B) -
#' 2\mathrm{Cov}(A,B)}. No model is refitted. JML, inference-ineligible fits,
#' missing or regularized covariance, and singular contrast covariance stop
#' with an error. Levels with fixed or unavailable contrasts are not silently
#' dropped. Known anchors are treated as fixed; their uncertainty is excluded.
#' Non-unit observation weights are inference-ineligible, including weights
#' normalized to mean one. Older bundles must retain the current readiness
#' contract as well as the covariance basis before they can be displayed.
#' A numerically verified inverse of ill-conditioned information is retained
#' with a warning, not silently treated as strong evidence. When applicable,
#' `cautions`, `information_review` and table column `InferenceCaution` preserve
#' that warning through summaries and plots. Numerical verification does not
#' establish the accuracy of the normal or chi-square approximations.
#'
#' The heterogeneity table uses a joint Wald chi-square test of equality of
#' the facet levels. Non-significant heterogeneity is neither necessary nor
#' sufficient for practical equivalence. `FixedChiSq`, `FixedDF`, and
#' `FixedProb` retain their column names but use this joint contrast test.
#' Separation and reliability remain descriptive summaries.
#'
#' `GrandMean` is the equally weighted mean of the facet estimates. For each
#' deviation from that mean, uncertainty includes the covariance with the
#' estimated mean. `ROPEPct` is the mass of its normal confidence distribution
#' inside the practical bound; it is descriptive, not a Bayesian posterior
#' probability or a separate equivalence decision.
#'
#' The former BIC/Bayes-factor heuristic is unavailable: a Wald statistic is
#' not a fitted likelihood comparison. For compatibility, `BF01` is `NA` and
#' `BF01Label` states why no value is supplied.
#'
#' @section What this analysis means:
#' The analysis asks whether differences between facet levels fall within a
#' prespecified practical bound under the fitted model. These are asymptotic
#' normal-approximation tests; numerical eligibility does not establish
#' finite-sample coverage or adequacy of the rating design.
#'
#' @section What this analysis does not justify:
#' A non-significant difference is not evidence of equivalence. Pairwise
#' conclusions are unadjusted for multiplicity: selecting some positive pairs
#' does not provide family-wise error control for that selected set. Estimated
#' linking or anchor uncertainty, population transport, and model
#' misspecification require separate evaluation.
#'
#' @section Decision rule:
#' Each pair uses two one-sided normal tests at alpha 0.05. `Equivalent` is
#' true when both tests reject non-equivalence, equivalently when its 90%
#' interval lies strictly inside the bound. `Decision` summarizes all pairs
#' as `"all_pairs_equivalent"`, `"partial_pairwise_equivalence"`, or
#' `"no_pairwise_equivalence_established"`. No pair may be omitted from the
#' all-pairs summary. Heterogeneity and ROPE summaries do not change this rule.
#' The all-pairs conclusion is an intersection-union decision: if any pair is
#' truly outside the bound, declaring all pairs equivalent requires falsely
#' passing that pair's test. Its asymptotic type-I error is therefore at most
#' alpha when each pair's approximation is valid. This does not provide
#' simultaneous confidence intervals or error control for selected pairs.
#'
#' @section Interpreting output:
#' Start with `summary$Decision` and examine the corresponding `pairwise`
#' differences, SEs and intervals. A negative result can reflect imprecision
#' or a material difference. `chi_square` addresses exact equality, while
#' `rope` and `forest` describe proximity to the facet mean.
#'
#' @section How to read the main outputs:
#' - `summary`: pairwise decision, covariance basis and multiplicity convention.
#' - `pairwise`: differences, covariance-aware SEs, 90% intervals and TOST tests.
#' - `chi_square`: joint Wald heterogeneity test and descriptive separation.
#' - `rope` / `forest`: `Measure`, marginal `SE` and `CI_Lower`/`CI_Upper`, plus
#'   `Deviation`, `DeviationSE` and `DeviationCI_Lower`/`DeviationCI_Upper` for
#'   the difference from the estimated facet mean. `Fixed` marks an anchored
#'   or constraint-fixed location: its conditional marginal `SE` is zero and
#'   its marginal interval is absent. Its deviation from an estimated mean,
#'   or contrast with an estimated level, can still have nonzero uncertainty.
#'   The forest plot displays these deviations, not the anchored locations.
#'   No uncertainty in supplied anchor values is propagated.
#' - `CIReason`, `DeviationCIReason`, and pairwise `CI90Reason` distinguish
#'   `"available"`, `"fixed"`, `"numerical_resolution"` (the requested width
#'   rounds away), and `"outside_numeric_range"`. Unavailable endpoints are
#'   `NA`; estimates and SEs remain. `ROPEStatus = "unavailable"` means the
#'   deviation interval could not be displayed, not that it is inside the bound.
#'   The fixed-alpha TOST p-values do not depend on the display `ci_level`.
#'
#' @section Recommended next step:
#' Review numerical integration and the model's uncertainty assumptions before
#' interpreting a borderline result. Sensitivity to another practical bound
#' should be reported transparently, without selecting a bound to obtain a
#' desired decision.
#'
#' @section Typical workflow:
#' 1. Fit and review an MML model with [fit_mfrm()].
#' 2. Prespecify the practical bound and run `analyze_facet_equivalence()`.
#' 3. Read `summary` and `pairwise`.
#' 4. Use [plot(x, type = "facet_equivalence")][plot.mfrm_fit] for descriptive grand-mean proximity.
#'
#' @section Output:
#' A bundle with `summary`, `chi_square`, `pairwise`, `rope`, `forest`, and
#' `settings`. Older bundles without the current inference/covariance basis
#' must be recomputed before using `summary()`, `print()`, or plotting.
#'
#' @return A named list with class `mfrm_facet_equivalence`.
#' @seealso [facets_chisq_table()], [fair_average_table()], [plot(x, type = "facet_equivalence")][plot.mfrm_fit]
#' @concept confidence intervals
#' @concept facet equivalence
#' @concept reporting workflow
#' @references
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#' procedure and the power approach for assessing the equivalence of
#' average bioavailability. *Journal of Pharmacokinetics and
#' Biopharmaceutics, 15*(6), 657-680.
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", quad_points = 31, maxit = 150)
#' eq <- analyze_facet_equivalence(fit, facet = "Rater")
#' eq$summary[, c("Facet", "Elements", "Decision", "MeanROPE")]
#' head(eq$pairwise[, c("ElementA", "ElementB", "Equivalent")])
#' @export
#' @order 1
analyze_facet_equivalence <- function(fit,
                                      diagnostics = NULL,
                                      facet = NULL,
                                      equivalence_bound = 0.5,
                                      ci_level = 0.95,
                                      conf_level = NULL) {
  stop_if_product_slopes(fit, "analyze_facet_equivalence()")
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (!is.numeric(equivalence_bound) || length(equivalence_bound) != 1 || !is.finite(equivalence_bound) ||
      equivalence_bound <= 0) {
    stop("`equivalence_bound` must be a positive finite number.", call. = FALSE)
  }
  # `conf_level` is the deprecated spelling; the canonical name
  # elsewhere in mfrmr is `ci_level`. When both are supplied we honor
  # `conf_level` and route the notification through
  # lifecycle so users can control verbosity with
  # options(lifecycle_verbosity = "..."). `conf_level` will be
  # deprecated; prefer `ci_level` in new code.
  if (!is.null(conf_level)) {
    lifecycle::deprecate_warn(
      when = "0.1.6",
      what = "analyze_facet_equivalence(conf_level = )",
      with = "analyze_facet_equivalence(ci_level = )"
    )
    ci_level <- conf_level
  }
  conf_level <- ci_level
  if (!is.numeric(conf_level) || length(conf_level) != 1 || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`ci_level` must be a single number between 0 and 1.", call. = FALSE)
  }

  if (!identical(fit$config$method, "MML") || !mfrm_inference_ready(fit)) {
    stop("Facet equivalence requires an inference-ready MML fit; ordinary inference is unavailable for this fit.",
         call. = FALSE)
  }
  if (!is.null(diagnostics)) {
    mfrm_results_validate_diagnostics_identity(
      fit, diagnostics, helper = "analyze_facet_equivalence()"
    )
    if (!isTRUE(diagnostics$precision_profile$SupportsFormalInference)) {
      stop("Supplied diagnostics do not support ordinary inference.", call. = FALSE)
    }
  }

  facet_names <- as.character(fit$config$facet_names)
  if (is.null(facet)) {
    facet <- infer_default_rater_facet(facet_names) %||% facet_names[1]
  }
  if (length(facet) != 1L || is.na(facet) || !facet %in% facet_names) {
    stop("`facet` must be one of: ", paste(facet_names, collapse = ", "), ".", call. = FALSE)
  }
  spec <- fit$config$facet_specs[[facet]]
  labels <- as.character(spec$levels)
  n_elem <- length(labels)
  if (n_elem < 2L) {
    stop("Facet '", facet, "' has fewer than 2 estimated levels.", call. = FALSE)
  }
  facet_df <- fit$facets$others[fit$facets$others$Facet == facet, , drop = FALSE]
  idx <- match(labels, as.character(facet_df$Level))
  est <- as.numeric(facet_df$Estimate[idx])
  if (anyNA(idx) || anyDuplicated(facet_df$Level) || !all(is.finite(est))) {
    stop("Complete finite facet estimates are required; no levels may be dropped.", call. = FALSE)
  }

  # Always derive uncertainty from this fit, never from supplied SE columns.
  covariance <- compute_mml_parameter_covariance(fit)
  if (!identical(covariance$status, "ok") || is.null(covariance$cov)) {
    stop("Facet equivalence requires unregularized MML observed-information covariance.", call. = FALSE)
  }
  slice <- covariance$param_slices[[facet]]
  jac <- constraint_jacobian(spec)
  fixed <- if (ncol(jac)) rowSums(abs(jac)) == 0 else rep(TRUE, n_elem)
  facet_cov <- symmetrize_matrix(
    jac %*% covariance$cov[slice, slice, drop = FALSE] %*% t(jac)
  )
  contrasts <- cbind(diag(n_elem - 1L), -1)
  contrast_rank <- qr(contrasts %*% jac)$rank
  if (contrast_rank < n_elem - 1L) {
    stop("All facet contrasts must be independent under the model constraints; fixed or unidentified contrasts cannot be tested.",
         call. = FALSE)
  }
  contrast_cov <- contrasts %*% facet_cov %*% t(contrasts)
  contrast_chol <- tryCatch(chol(contrast_cov), error = function(e) NULL)
  if (!all(is.finite(facet_cov)) || is.null(contrast_chol)) {
    stop("All facet contrasts must have positive-definite covariance; fixed or unidentified contrasts cannot be tested.",
         call. = FALSE)
  }
  se <- covariance_diag_se(facet_cov)
  if (!all(is.finite(se))) {
    stop("Facet covariance has invalid marginal variances.", call. = FALSE)
  }
  if (!is.null(diagnostics)) {
    rows <- diagnostics$measures[diagnostics$measures$Facet == facet, , drop = FALSE]
    row_idx <- match(labels, as.character(rows$Level))
    if (!all(c("Level", "SE", "SupportsFormalInference") %in% names(rows)) ||
        nrow(rows) != n_elem || anyNA(row_idx) || anyDuplicated(rows$Level) ||
        !isTRUE(all(rows$SupportsFormalInference[row_idx[!fixed]])) ||
        !isTRUE(all.equal(as.numeric(rows$SE[row_idx[!fixed]]), se[!fixed], tolerance = 1e-8,
                          check.attributes = FALSE)) ||
        (any(fixed) && !("Fixed" %in% names(rows) &&
          all(rows$Fixed[row_idx[fixed]] %in% TRUE) && all(is.na(rows$SE[row_idx[fixed]]))))) {
      stop("Supplied facet diagnostics must contain all levels with matching model-based SEs and ordinary-inference eligibility for estimated levels, and fixed-value labels for constants.",
           call. = FALSE)
    }
  }

  grand_mean <- mean(est)
  centering <- diag(n_elem) - 1 / n_elem
  deviation_se <- covariance_diag_se(centering %*% facet_cov %*% t(centering))
  df_chi <- n_elem - 1L
  standardized <- forwardsolve(t(contrast_chol), contrasts %*% est)
  chi2_val <- sum(standardized ^ 2)
  p_chi <- stats::pchisq(chi2_val, df = df_chi, lower.tail = FALSE)
  sep_sd <- stats::sd(est)
  rmse <- sqrt(mean(se ^ 2))
  true_sd <- sqrt(max(sep_sd ^ 2 - rmse ^ 2, 0))
  separation <- if (rmse > 0) true_sd / rmse else NA_real_
  reliability <- separation ^ 2 / (1 + separation ^ 2)

  z_ci <- stats::qnorm((1 - conf_level) / 2, lower.tail = FALSE)
  z_tost <- stats::qnorm(0.95)
  interval <- function(estimate, se, z, fixed = rep(FALSE, length(estimate))) {
    lower <- estimate - z * se
    upper <- estimate + z * se
    reason <- rep("available", length(estimate))
    reason[!is.finite(lower) | !is.finite(upper)] <- "outside_numeric_range"
    unresolved <- reason == "available" & !(lower < estimate & upper > estimate)
    reason[unresolved] <- "numerical_resolution"
    reason[fixed] <- "fixed"
    lower[reason != "available"] <- upper[reason != "available"] <- NA_real_
    list(lower = lower, upper = upper, reason = reason)
  }
  pair_idx <- utils::combn(seq_len(n_elem), 2)
  pairwise_tbl <- if (is.null(dim(pair_idx))) {
    data.frame()
  } else {
    pairwise_rows <- lapply(seq_len(ncol(pair_idx)), function(k) {
      i <- pair_idx[1, k]
      j <- pair_idx[2, k]
      diff <- est[i] - est[j]
      se_diff <- sqrt(facet_cov[i, i] + facet_cov[j, j] - 2 * facet_cov[i, j])
      z_lower <- (diff + equivalence_bound) / se_diff
      z_upper <- (diff - equivalence_bound) / se_diff
      p_lower <- stats::pnorm(z_lower, lower.tail = FALSE)
      p_upper <- stats::pnorm(z_upper, lower.tail = TRUE)
      p_tost <- max(p_lower, p_upper)
      ci <- interval(diff, se_diff, z_tost)
      data.frame(
        ElementA = labels[i],
        ElementB = labels[j],
        Diff = diff,
        SE_Diff = se_diff,
        CI90_Lower = ci$lower,
        CI90_Upper = ci$upper,
        CI90Reason = ci$reason,
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

  marginal_ci <- interval(est, se, z_ci, fixed)
  deviation_ci <- interval(est - grand_mean, deviation_se, z_ci)
  rope_tbl <- data.frame(
    Element = labels,
    Fixed = fixed,
    Measure = est,
    Deviation = est - grand_mean,
    SE = se,
    DeviationSE = deviation_se,
    DeviationCI_Lower = deviation_ci$lower,
    DeviationCI_Upper = deviation_ci$upper,
    DeviationCIReason = deviation_ci$reason,
    CI_Lower = marginal_ci$lower,
    CI_Upper = marginal_ci$upper,
    CIReason = marginal_ci$reason,
    stringsAsFactors = FALSE
  )
  rope_tbl$ROPEPct <- 100 * vapply(seq_len(n_elem), function(i) {
    # Reflection avoids subtracting two lower-tail probabilities near one.
    center <- abs(rope_tbl$Deviation[i]); sd <- rope_tbl$DeviationSE[i]
    upper <- stats::pnorm(equivalence_bound, center, sd)
    mass <- upper - stats::pnorm(-equivalence_bound, center, sd)
    if (mass < upper * 1e-6) {
      # Integrate on a fixed domain when a narrow band loses CDF precision.
      half_width <- equivalence_bound / sd
      mass <- half_width * stats::integrate(function(t)
        stats::dnorm(half_width * t - center / sd), -1, 1,
        rel.tol = 1e-10, abs.tol = 0)$value
    }
    mass
  }, numeric(1))
  rope_tbl$ROPEStatus <- ifelse(
    rope_tbl$DeviationCI_Lower >= -equivalence_bound & rope_tbl$DeviationCI_Upper <= equivalence_bound,
    "inside",
    ifelse(
      rope_tbl$DeviationCI_Lower > equivalence_bound | rope_tbl$DeviationCI_Upper < -equivalence_bound,
      "outside",
      "overlap"
    )
  )
  rope_tbl$ROPEStatus[rope_tbl$DeviationCIReason != "available"] <- "unavailable"

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
    TestBasis = "Wald test of joint facet contrasts",
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
    BF01 = NA_real_,
    BF01Label = "Unavailable: requires fitted likelihood comparison",
    MeanROPE = mean_rope,
    AllPairsEquivalent = all_pairs_equivalent,
    AnyPairEquivalent = any_pair_equivalent,
    PairwiseDecisionBasis = "pairwise_tost_summary",
    InferenceReady = TRUE,
    CovarianceBasis = "mml_observed_information_contrasts",
    MultiplicityAdjustment = "none",
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
      ci_level = conf_level,
      mean_basis = "equal_weight_facet_mean",
      contrast_rank = contrast_rank,
      readiness_contract_version = mfrmr_readiness_contract_version(),
      covariance_basis = "mml_observed_information_contrasts"
    )
  )
  if (any(c(pairwise_tbl$CI90Reason, rope_tbl$CIReason,
            rope_tbl$DeviationCIReason) %in% c("numerical_resolution", "outside_numeric_range"))) {
    out$notes <- "Some interval endpoints are unavailable at the requested numeric precision; inspect CIReason, DeviationCIReason and CI90Reason. Estimates, SEs and the fixed-alpha TOST remain separate from display intervals."
  }
  cautions <- mfrm_mml_information_caution(covariance)
  if (length(cautions)) {
    out$cautions <- cautions
    out$information_review <- covariance$solution_information$inverse_review
    for (name in c("summary", "chi_square", "pairwise", "rope", "forest")) {
      out[[name]]$InferenceCaution <- paste(cautions, collapse = " ")
    }
    warning(paste(cautions, collapse = " "), call. = FALSE)
  }
  as_mfrm_bundle(out, "mfrm_facet_equivalence")
}

validate_facet_equivalence_bundle <- function(x) {
  if (!isTRUE(x$summary$InferenceReady) ||
      !identical(x$settings$readiness_contract_version,
                 mfrmr_readiness_contract_version()) ||
      !identical(x$settings$contrast_rank, x$summary$Elements - 1L) ||
      !identical(x$settings$covariance_basis, "mml_observed_information_contrasts") ||
      !identical(x$summary$CovarianceBasis, "mml_observed_information_contrasts")) {
    stop("This equivalence bundle lacks the current inference and covariance basis. Recompute it with analyze_facet_equivalence() from an eligible MML fit.",
         call. = FALSE)
  }
  invisible(x)
}

#' Plot facet-equivalence results
#'
#' @param x Output from [analyze_facet_equivalence()] or an eligible MML
#'   [fit_mfrm()] object. Legacy equivalence bundles must be recomputed.
#' @param diagnostics Optional matching output from [diagnose_mfrm()] when
#'   `x` is an `mfrm_fit` object.
#' @param facet Facet to analyze when `x` is an `mfrm_fit` object.
#' @param type Plot type: `"forest"` (default) or `"rope"`.
#' @param draw If `TRUE` (default), draw the plot. If `FALSE`, return the
#'   prepared plotting data.
#' @param ... Additional graphical arguments passed to base plotting functions.
#'
#' @details
#' Fit inputs use the same eligibility checks as [analyze_facet_equivalence()].
#' Bundle inputs display the already calculated results. Both routes require
#' the current inference and covariance basis, including when `draw = FALSE`.
#'
#' @section Plot types:
#' - `"forest"` shows each level's deviation from the equally weighted facet
#'   mean, with covariance-aware deviation intervals and the practical region
#'   around zero. The raw marginal measure intervals remain in the data table.
#' - `"rope"` shows the normal confidence-distribution mass within that region.
#'
#' @section Interpreting output:
#' Both plots describe grand-mean proximity. Colors in the forest plot indicate
#' whether the deviation interval is inside, outside, or overlaps the practical
#' region. Neither plot establishes pairwise equivalence or a Bayesian
#' probability. Read the pairwise TOST results for pair-specific conclusions.
#'
#' @section Typical workflow:
#' 1. Run [analyze_facet_equivalence()] with a prespecified practical bound.
#' 2. Use `type = "forest"` to inspect deviations and their uncertainty.
#' 3. Use `type = "rope"` for a descriptive proximity view.
#'
#' @return Invisibly returns the plotting data and inference/covariance basis.
#'   With `draw = FALSE`, returns the data without drawing.
#' @seealso [analyze_facet_equivalence()]
#' @concept confidence intervals
#' @concept facet equivalence
#' @concept visual diagnostics
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot_facet_equivalence", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", quad_points = 31, maxit = 150)
#' eq <- analyze_facet_equivalence(fit, facet = "Rater")
#' pdat <- plot(eq, type = "facet_equivalence", view = "forest", draw = FALSE)
#' c(pdat$facet, pdat$type)
#' }
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
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

  validate_facet_equivalence_bundle(x)

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
    type = type,
    title = paste0(as.character(settings$facet), ": deviation from facet mean"),
    note = paste(c("Mean-deviation intervals do not establish pairwise equivalence.",
      x$notes, x$cautions), collapse = " "),
    cautions = x$cautions,
    information_review = x$information_review,
    inference_ready = TRUE,
    covariance_basis = settings$covariance_basis
  )
  if (!isTRUE(draw)) {
    return(out)
  }

  dots <- list(...)

  if (identical(type, "forest")) {
    ord <- order(forest_df$Deviation, decreasing = FALSE, na.last = TRUE)
    forest_df <- forest_df[ord, , drop = FALSE]
    ypos <- seq_len(nrow(forest_df))
    cols <- ifelse(
      forest_df$ROPEStatus == "inside",
      "#2E8B57",
      ifelse(forest_df$ROPEStatus == "outside", "#C0392B",
        ifelse(forest_df$ROPEStatus == "unavailable", "grey50", "#D68910"))
    )
    xlim <- range(c(forest_df$Deviation, forest_df$DeviationCI_Lower, forest_df$DeviationCI_Upper,
                    -out$equivalence_bound, out$equivalence_bound), finite = TRUE)
    do.call(graphics::plot, c(list(
      x = forest_df$Deviation,
      y = ypos,
      xlim = xlim,
      yaxt = "n",
      ylab = "",
      xlab = "Deviation from facet mean (logits)",
      main = "",
      pch = 19,
      col = cols
    ), dots))
    graphics::axis(2, at = ypos, labels = forest_df$Element, las = 2)
    graphics::rect(
      xleft = -out$equivalence_bound,
      ybottom = 0.5,
      xright = out$equivalence_bound,
      ytop = nrow(forest_df) + 0.5,
      border = NA,
      col = grDevices::adjustcolor("#2E8B57", alpha.f = 0.12)
    )
    graphics::abline(v = 0, lty = 2, col = "gray40")
    graphics::segments(forest_df$DeviationCI_Lower, ypos, forest_df$DeviationCI_Upper, ypos, col = cols, lwd = 2)
    graphics::points(forest_df$Deviation, ypos, pch = 19, col = cols)
  } else {
    ord <- order(forest_df$ROPEPct, decreasing = TRUE, na.last = TRUE)
    forest_df <- forest_df[ord, , drop = FALSE]
    cols <- "#4477AA"
    mids <- graphics::barplot(
      height = forest_df$ROPEPct,
      names.arg = forest_df$Element,
      las = 2,
      ylim = c(0, 100),
      col = cols,
      ylab = "Normal reference mass within bound (%)",
      main = "",
      ...
    )
    out$bar_midpoints <- mids
  }

  invisible(out)
}

#' @rdname analyze_facet_equivalence
#' @details `plot()` method: Plot a saved facet-equivalence analysis
#' @param x A result from [analyze_facet_equivalence()].
#' @param ... Use `type = "forest"` (default) for deviations and intervals or
#'   `type = "rope"` for normal confidence-distribution mass. Other controls
#'   are described in [plot_facet_equivalence()].
#' @details For `plot()`: When an interval cannot be represented at the requested numeric precision,
#' the forest retains its estimate as a grey point without a whisker. Inspect
#' `data$DeviationCIReason` and `note` in the returned payload for the reason.
#' This does not change the pairwise TOST decision stored in the analysis.
#' @return For `plot()`: The existing equivalence plot payload and inference restrictions.
#'   Grand-mean proximity does not establish pairwise equivalence.
#' @seealso [mfrmr_plot_migration]
#' @export
plot.mfrm_facet_equivalence <- function(x, ...) {
  args <- list(...)
  if (identical(args$type, "facet_equivalence")) {
    args$type <- NULL
    args <- .mfrm_plot_rename(args, "view", "type")
  }
  .mfrm_plot_call(plot_facet_equivalence, x, args)
}
