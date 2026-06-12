# ==============================================================================
# Empirical-Bayes and Laplace shrinkage for small-N facets (added in 0.1.6)
# ==============================================================================
#
# mfrmr estimates every non-person facet as a fixed effect with a sum-to-zero
# identification constraint (see `?fit_mfrm` "Fixed effects assumption").
# That keeps the model close to the Linacre (1989) specification, but when a
# facet has only a handful of levels (e.g. 3 raters in a language-testing
# study) the resulting severity estimates carry wide standard errors and
# are not "shrunk" toward the facet mean.
#
# The helpers below add an optional, post-hoc partial-pooling layer:
#
# - `.compute_facet_shrinkage()` implements the James-Stein / empirical-
#   Bayes formula (Efron & Morris, 1973; see references below).
# - `apply_empirical_bayes_shrinkage()` is the public wrapper that takes an
#   existing `mfrm_fit` and returns a fit augmented with shrunk columns and
#   a `shrinkage_report` table.
# - `fit_mfrm(..., facet_shrinkage = "empirical_bayes")` wires the same
#   helper into the primary entry point so users can opt in up front.
#
# The default `facet_shrinkage = "none"` keeps the 0.1.5 / 0.1.6 behaviour
# unchanged; every new column is additive.
#
# References:
# - Efron, B., & Morris, C. (1973). Combining possibly related estimation
#   problems. Journal of the Royal Statistical Society: Series B, 35(3),
#   379-402.
# - Efron, B. (2021). Empirical Bayes: Concepts and methods. Stanford
#   University technical report. <https://efron.ckirby.su.domains/papers/
#   2021EB-concepts-methods.pdf>
# - Morris, C. N. (1983). Parametric empirical Bayes inference: Theory and
#   applications. Journal of the American Statistical Association, 78(381),
#   47-55.
# - Wright, B. D. (1998). Estimating Rasch measures for extreme scores.
#   Rasch Measurement Transactions, 12(2), 632-633.
#   <https://www.rasch.org/rmt/rmt122h.htm>
# - Jones, E., & Wind, S. A. (2018). Using repeated ratings to improve
#   measurement precision in incomplete rating designs. Journal of Applied
#   Measurement, 19(2), 148-161.
# ==============================================================================


#' @keywords internal
#' @noRd
.compute_facet_shrinkage <- function(estimates, ses,
                                     method = c("empirical_bayes"),
                                     prior_sd = NULL,
                                     min_levels = 3L) {
  # Core math:
  #   tau^2      = max(0, sum(estimates^2)/K - mean(se^2))  (method-of-moments,
  #                                                          under the sum-to-zero
  #                                                          identification the
  #                                                          facet mean is exactly
  #                                                          zero -- no DF lost to
  #                                                          mean estimation)
  #   B_j        = se_j^2 / (tau^2 + se_j^2)                (shrinkage factor)
  #   est^EB_j   = (1 - B_j) * est_j                         (SZ -> shrink to 0)
  #   se^EB_j    = sqrt((1 - B_j) * se_j^2)                  (naive Morris 1983
  #                                                          posterior SE;
  #                                                          tau^2 treated as known)
  #
  # Returns a list with parallel shrunk vectors so downstream callers can
  # merge them back onto the facet table by row order.
  method <- match.arg(method)

  estimates <- as.numeric(estimates)
  ses <- as.numeric(ses)
  if (length(estimates) != length(ses)) {
    stop("`estimates` and `ses` must have the same length.", call. = FALSE)
  }

  K <- length(estimates)
  out_template <- list(
    shrunk_estimates = rep(NA_real_, K),
    shrunk_ses = rep(NA_real_, K),
    shrinkage_factors = rep(NA_real_, K),
    tau2 = NA_real_,
    mean_se2 = NA_real_,
    n_levels = K,
    n_levels_used = 0L,
    method = method,
    prior_sd_source = "empirical",
    note = NA_character_
  )

  valid <- is.finite(estimates) & is.finite(ses) & ses > 0
  K_use <- sum(valid)
  if (K_use < min_levels) {
    out_template$note <- paste0(
      "Fewer than ", min_levels, " valid levels; shrinkage not applied."
    )
    # Pass through unchanged.
    out_template$shrunk_estimates <- estimates
    out_template$shrunk_ses <- ses
    out_template$shrinkage_factors <- rep(0, K)
    return(out_template)
  }

  mean_se2 <- mean(ses[valid]^2)
  # Population variance of point estimates around zero (sum-to-zero gives
  # the population mean = 0). Using the raw second moment makes the
  # estimator consistent under the exchangeable prior assumption.
  raw_var <- sum(estimates[valid]^2) / K_use

  # Pull prior variance either from the supplied value or from the
  # method-of-moments estimator. The typical choice is MoM.
  if (!is.null(prior_sd)) {
    tau2 <- as.numeric(prior_sd)^2
    out_template$prior_sd_source <- "user"
  } else {
    tau2 <- max(0, raw_var - mean_se2)
    out_template$prior_sd_source <- "empirical"
  }

  # Compute per-level shrinkage. Levels with non-finite SEs pass through
  # with factor 0 (no shrinkage applied).
  B <- rep(0, K)
  shrunk_est <- estimates
  shrunk_se <- ses
  if (tau2 > 0) {
    B[valid] <- ses[valid]^2 / (tau2 + ses[valid]^2)
    shrunk_est[valid] <- (1 - B[valid]) * estimates[valid]
    shrunk_se[valid] <- sqrt((1 - B[valid]) * ses[valid]^2)
  } else {
    # tau^2 collapsed to zero: signal that no between-level variance is
    # discernible above measurement error. Shrinkage factor = 1 collapses
    # every level to 0; that is the correct EB answer but can look
    # surprising, so we also record a note for downstream reporting.
    B[valid] <- 1
    shrunk_est[valid] <- 0
    shrunk_se[valid] <- sqrt(pmax(0, 0)) # zero variance after complete pooling
    out_template$note <- paste0(
      "Between-level variance below measurement error; ",
      "all levels collapsed to the facet mean."
    )
  }

  out_template$shrunk_estimates <- shrunk_est
  out_template$shrunk_ses <- shrunk_se
  out_template$shrinkage_factors <- B
  out_template$tau2 <- tau2
  out_template$mean_se2 <- mean_se2
  out_template$n_levels_used <- K_use
  out_template
}


#' @keywords internal
#' @noRd
.apply_shrinkage_to_fit <- function(fit,
                                    method = c("empirical_bayes", "laplace"),
                                    facet_prior_sd = NULL,
                                    shrink_person = FALSE) {
  method <- match.arg(method)
  if (!inherits(fit, "mfrm_fit")) {
    stop("`.apply_shrinkage_to_fit()` requires an mfrm_fit.", call. = FALSE)
  }

  others <- fit$facets$others
  if (is.null(others) || nrow(others) == 0L) {
    fit$shrinkage_report <- data.frame(
      Facet = character(0), NLevels = integer(0), NLevelsUsed = integer(0),
      Tau2 = numeric(0), MeanSE2 = numeric(0), MeanShrinkage = numeric(0),
      EffectiveDF = numeric(0), Method = character(0),
      PriorSource = character(0), Note = character(0),
      stringsAsFactors = FALSE
    )
    fit$config$facet_shrinkage <- as.character(method)
    fit$config$facet_prior_sd <- facet_prior_sd
    return(fit)
  }

  others <- as.data.frame(others, stringsAsFactors = FALSE)
  others$ShrunkEstimate <- NA_real_
  others$ShrunkSE <- NA_real_
  others$ShrinkageFactor <- NA_real_

  # Per-level SEs live on `diagnostics$measures`, not on
  # `fit$facets$others`. Pull them once; callers that want to avoid the
  # diagnose_mfrm() cost can pre-compute and stash SE columns on
  # `fit$facets$others` (e.g. via `ModelSE`) before invoking shrinkage.
  se_col <- if ("ModelSE" %in% names(others)) {
    "ModelSE"
  } else if ("SE" %in% names(others)) {
    "SE"
  } else {
    NA_character_
  }
  if (is.na(se_col)) {
    measures <- tryCatch(
      suppressMessages(suppressWarnings(
        diagnose_mfrm(fit, residual_pca = "none")
      ))$measures,
      error = function(e) NULL
    )
    if (is.null(measures) ||
        !all(c("Facet", "Level", "ModelSE") %in% names(measures))) {
      stop("Shrinkage requires per-level standard errors. ",
           "Run diagnose_mfrm(fit) and stash SE on fit$facets$others ",
           "before calling, or use the default fit_mfrm(..., ",
           "facet_shrinkage = ...) path.", call. = FALSE)
    }
    measures <- as.data.frame(measures, stringsAsFactors = FALSE)
    measures$Level <- as.character(measures$Level)
    others$Level <- as.character(others$Level)
    others <- merge(
      others,
      measures[, c("Facet", "Level", "ModelSE")],
      by = c("Facet", "Level"),
      all.x = TRUE,
      sort = FALSE
    )
    se_col <- "ModelSE"
  }

  facets <- unique(as.character(others$Facet))
  report_rows <- list()
  for (f in facets) {
    idx <- which(as.character(others$Facet) == f)
    if (length(idx) == 0L) next
    res <- .compute_facet_shrinkage(
      estimates = others$Estimate[idx],
      ses = others[[se_col]][idx],
      method = if (identical(method, "laplace")) "empirical_bayes" else method,
      prior_sd = facet_prior_sd,
      min_levels = 3L
    )
    others$ShrunkEstimate[idx] <- res$shrunk_estimates
    others$ShrunkSE[idx] <- res$shrunk_ses
    others$ShrinkageFactor[idx] <- res$shrinkage_factors
    effective_df <- res$n_levels_used - sum(res$shrinkage_factors[
      is.finite(res$shrinkage_factors)
    ])
    report_rows[[f]] <- data.frame(
      Facet = f,
      NLevels = res$n_levels,
      NLevelsUsed = res$n_levels_used,
      Tau2 = res$tau2,
      MeanSE2 = res$mean_se2,
      MeanShrinkage = if (res$n_levels_used > 0) {
        mean(res$shrinkage_factors[is.finite(res$shrinkage_factors)])
      } else NA_real_,
      EffectiveDF = effective_df,
      Method = res$method,
      PriorSource = res$prior_sd_source,
      Note = res$note %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }

  fit$facets$others <- others
  fit$shrinkage_report <- do.call(rbind, report_rows)
  rownames(fit$shrinkage_report) <- NULL

  # Record settings so downstream consumers (manifest, replay, checklist,
  # APA narrative) can report them consistently.
  fit$config$facet_shrinkage <- as.character(method)
  fit$config$facet_prior_sd <- facet_prior_sd

  # Optional: shrink Person estimates too. JML exposes each theta as a
  # fixed effect so EB has real bite; MML already integrates the prior
  # so shrinkage is usually redundant but kept behind an explicit opt-in.
  if (isTRUE(shrink_person)) {
    person <- fit$facets$person
    if (!is.null(person) && "Estimate" %in% names(person) &&
        "SE" %in% names(person)) {
      res_p <- .compute_facet_shrinkage(
        estimates = person$Estimate,
        ses = person$SE,
        method = "empirical_bayes",
        prior_sd = NULL,
        min_levels = 3L
      )
      person$ShrunkEstimate <- res_p$shrunk_estimates
      person$ShrunkSE <- res_p$shrunk_ses
      person$ShrinkageFactor <- res_p$shrinkage_factors
      fit$facets$person <- person
      person_report <- data.frame(
        Facet = "Person",
        NLevels = res_p$n_levels,
        NLevelsUsed = res_p$n_levels_used,
        Tau2 = res_p$tau2,
        MeanSE2 = res_p$mean_se2,
        MeanShrinkage = if (res_p$n_levels_used > 0) {
          mean(res_p$shrinkage_factors[is.finite(res_p$shrinkage_factors)])
        } else NA_real_,
        EffectiveDF = res_p$n_levels_used -
          sum(res_p$shrinkage_factors[is.finite(res_p$shrinkage_factors)]),
        Method = res_p$method,
        PriorSource = res_p$prior_sd_source,
        Note = res_p$note %||% NA_character_,
        stringsAsFactors = FALSE
      )
      fit$shrinkage_report <- rbind(fit$shrinkage_report, person_report)
      rownames(fit$shrinkage_report) <- NULL
    }
  }

  # Summary-level surfacing. Downstream tables and narrative read these.
  if (!is.null(fit$summary) && nrow(fit$summary) > 0L) {
    fit$summary$FacetShrinkage <- as.character(method)
    fit$summary$FacetShrinkageTau2Mean <- if (nrow(fit$shrinkage_report) > 0L) {
      mean(fit$shrinkage_report$Tau2, na.rm = TRUE)
    } else NA_real_
  }

  fit
}


#' Apply empirical-Bayes shrinkage to fitted non-person facet estimates
#'
#' Post-hoc shrinkage helper that augments an `mfrm_fit` with James-Stein
#' / empirical-Bayes shrunk estimates for each non-person facet. The
#' shrinkage variance \eqn{\hat{\tau}^2} is estimated by method of
#' moments from the facet-level point estimates and their standard
#' errors:
#' \deqn{\hat{\tau}^2 = \max\!\left(0,
#'   \frac{1}{K}\sum_{j=1}^{K}\hat{\delta}_j^{2} -
#'   \overline{\mathrm{SE}^2}\right),}
#' where the first term is the population variance of the facet point
#' estimates around their *known* mean of zero (the mfrmr sum-to-zero
#' identification pins the facet mean exactly at 0, so no degree of
#' freedom is consumed by mean estimation). The shrinkage factor is
#' \eqn{B_j = \mathrm{SE}_j^2 / (\hat{\tau}^2 + \mathrm{SE}_j^2)}, and
#' the shrunk point / standard error are
#' \eqn{\hat{\delta}_j^{EB} = (1 - B_j)\hat{\delta}_j} and
#' \eqn{\mathrm{SE}_j^{EB} = \sqrt{(1 - B_j)\mathrm{SE}_j^2}}.
#' The posterior SE form treats \eqn{\hat{\tau}^2} as known; it omits
#' the Morris (1983, eqs. 4.1-4.2, p. 51) confidence-interval correction
#' \eqn{v \cdot \hat{\delta}_j^{2}} with
#' \eqn{v = 2 B_j^2 / (K - r - 2)}, where \eqn{r} is the number of
#' regression coefficients used to model the prior mean (under mfrmr's
#' sum-to-zero pinning, \eqn{r = 0}, so the divisor is \eqn{K - 2}).
#' This correction adds variance proportional to the squared deviation
#' \eqn{\hat{\delta}_j^{2}}, accounting for uncertainty in
#' \eqn{\hat{\tau}^2}. Under the equal-variance assumption
#' \eqn{\hat{\delta}_j^{2} \approx \hat{\tau}^2}, the omitted variance is
#' on the order of \eqn{2 / (K - 2)} times the reported posterior
#' variance \eqn{V(1 - B_j)}, so the true SE is approximately
#' \eqn{\sqrt{1 + 2/(K - 2)}} times the reported `ShrunkSE`. Magnitudes:
#' SE understated by ~73\% at \eqn{K = 3}, ~29\% at \eqn{K = 5}, ~15\%
#' at \eqn{K = 8}, ~7\% at \eqn{K = 15}. For a small-K facet, treat
#' `ShrunkSE` as a lower bound rather than a calibrated posterior SE.
#'
#' `fit$facets$others` gains `ShrunkEstimate`, `ShrunkSE`, and
#' `ShrinkageFactor` columns, and `fit$shrinkage_report` records the
#' per-facet \eqn{\hat{\tau}^2}, mean shrinkage, and effective degrees
#' of freedom (\eqn{\mathrm{EffectiveDF}_f = \sum_j (1 - B_j)}, which
#' matches the "effective number of parameters" defined by
#' Efron & Morris, 1973). The original `Estimate` / `SE` columns are
#' preserved.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()] with a non-empty
#'   `facets$others` table.
#' @param facet_prior_sd Optional numeric scalar. When supplied, the
#'   shrinkage variance is fixed at `facet_prior_sd^2` instead of being
#'   estimated from the data. Useful when a prior is elicited from
#'   expert knowledge or a previous fit.
#' @param shrink_person Logical. When `TRUE`, the same empirical-Bayes
#'   shrinkage is also applied to `fit$facets$person`. Default `FALSE`,
#'   since MML person estimates already reflect a N(0, sigma^2) prior.
#'
#' @return The same `mfrm_fit`, with augmented columns and a new
#'   `shrinkage_report` list entry, and with
#'   `fit$config$facet_shrinkage` set to `"empirical_bayes"`.
#'
#' @section Typical workflow:
#' 1. Fit the model as usual with `fit_mfrm()`.
#' 2. Call `apply_empirical_bayes_shrinkage(fit)` when small-N facets
#'    are present (see [facet_small_sample_review()]).
#' 3. Report both the original and shrunk estimates in the manuscript,
#'    citing Efron & Morris (1973). `build_apa_outputs()` will add the
#'    sentence automatically when `fit$config$facet_shrinkage` is set.
#'
#' @seealso [fit_mfrm()] (which accepts `facet_shrinkage` directly),
#'   [facet_small_sample_review()], [compute_facet_icc()].
#'
#' @references
#' Efron, B., & Morris, C. (1973). Combining possibly related
#' estimation problems. *Journal of the Royal Statistical Society:
#' Series B, 35*(3), 379-402.
#'
#' Efron, B. (2021). *Empirical Bayes: Concepts and methods*
#' (Technical report). Department of Statistics, Stanford University.
#' <https://efron.ckirby.su.domains/papers/2021EB-concepts-methods.pdf>
#'
#' Morris, C. N. (1983). Parametric empirical Bayes inference: Theory
#' and applications. *Journal of the American Statistical Association,
#' 78*(381), 47-55.
#'
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' fit_eb <- apply_empirical_bayes_shrinkage(fit)
#' fit_eb$shrinkage_report
#' # Look for:
#' # - `Tau2` is the estimated between-level prior variance per facet.
#' #   `Tau2 = 0` means the data did not justify any pooling and the
#' #   shrunken estimates equal the raw estimates (`MeanShrinkage = 0`).
#' # - `MeanShrinkage` near 0 = little movement, near 1 = heavy pooling
#' #   toward 0. Small-N facets typically pull values further than
#' #   well-identified ones.
#' # - `EffectiveDF` is the implied "effective number of parameters"
#' #   (Efron & Morris 1973); EffectiveDF much smaller than the row
#' #   count of the facet means most levels were pooled together.
#' head(fit_eb$facets$others[, c("Facet", "Level", "Estimate",
#'                                "ShrunkEstimate", "ShrinkageFactor")])
#' # Look for: rows where `ShrinkageFactor` is large (close to 1) had
#' #   their estimates pulled most strongly toward the facet mean (0).
#' @export
apply_empirical_bayes_shrinkage <- function(fit,
                                            facet_prior_sd = NULL,
                                            shrink_person = FALSE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit from fit_mfrm().", call. = FALSE)
  }
  .apply_shrinkage_to_fit(
    fit = fit,
    method = "empirical_bayes",
    facet_prior_sd = facet_prior_sd,
    shrink_person = shrink_person
  )
}


#' @keywords internal
#' @noRd
.build_shrinkage_plot_data <- function(fit) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`.build_shrinkage_plot_data()` requires an mfrm_fit.", call. = FALSE)
  }
  others <- fit$facets$others %||% data.frame()
  report <- fit$shrinkage_report
  mode <- as.character(fit$config$facet_shrinkage %||% "none")
  if (identical(mode, "none") || is.null(others) || nrow(others) == 0L ||
      !"ShrunkEstimate" %in% names(others)) {
    return(list(
      table = data.frame(
        Facet = character(0), Level = character(0),
        Estimate = numeric(0), SE = numeric(0),
        ShrunkEstimate = numeric(0), ShrunkSE = numeric(0),
        ShrinkageFactor = numeric(0),
        stringsAsFactors = FALSE
      ),
      report = report,
      mode = mode
    ))
  }
  se_col <- if ("ModelSE" %in% names(others)) "ModelSE" else "SE"
  tbl <- data.frame(
    Facet = as.character(others$Facet),
    Level = as.character(others$Level),
    Estimate = suppressWarnings(as.numeric(others$Estimate)),
    SE = if (se_col %in% names(others)) {
      suppressWarnings(as.numeric(others[[se_col]]))
    } else {
      rep(NA_real_, nrow(others))
    },
    ShrunkEstimate = suppressWarnings(as.numeric(others$ShrunkEstimate)),
    ShrunkSE = suppressWarnings(as.numeric(others$ShrunkSE)),
    ShrinkageFactor = suppressWarnings(as.numeric(others$ShrinkageFactor)),
    stringsAsFactors = FALSE
  )
  list(table = tbl, report = report, mode = mode)
}

#' @keywords internal
#' @noRd
.validate_shrinkage_ci_level <- function(ci_level) {
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number in (0, 1).", call. = FALSE)
  }
  as.numeric(ci_level)
}

#' @keywords internal
#' @noRd
.add_shrinkage_ci_columns <- function(tbl,
                                      ci_level = 0.95,
                                      raw_estimate_col = "Estimate",
                                      raw_se_col = "SE",
                                      raw_prefix = "",
                                      shrunk_estimate_col = "ShrunkEstimate",
                                      shrunk_se_col = "ShrunkSE",
                                      shrunk_prefix = "Shrunk") {
  ci_level <- .validate_shrinkage_ci_level(ci_level)
  z_ci <- stats::qnorm(1 - (1 - ci_level) / 2)
  raw_lower <- if (nzchar(raw_prefix)) paste0(raw_prefix, "CI_Lower") else "CI_Lower"
  raw_upper <- if (nzchar(raw_prefix)) paste0(raw_prefix, "CI_Upper") else "CI_Upper"
  shrunk_lower <- paste0(shrunk_prefix, "CI_Lower")
  shrunk_upper <- paste0(shrunk_prefix, "CI_Upper")

  raw_est <- if (raw_estimate_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[raw_estimate_col]])) else rep(NA_real_, nrow(tbl))
  raw_se <- if (raw_se_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[raw_se_col]])) else rep(NA_real_, nrow(tbl))
  shrunk_est <- if (shrunk_estimate_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[shrunk_estimate_col]])) else rep(NA_real_, nrow(tbl))
  shrunk_se <- if (shrunk_se_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[shrunk_se_col]])) else rep(NA_real_, nrow(tbl))

  tbl[[raw_lower]] <- ifelse(is.finite(raw_est) & is.finite(raw_se),
                             raw_est - z_ci * raw_se, NA_real_)
  tbl[[raw_upper]] <- ifelse(is.finite(raw_est) & is.finite(raw_se),
                             raw_est + z_ci * raw_se, NA_real_)
  tbl[[shrunk_lower]] <- ifelse(is.finite(shrunk_est) & is.finite(shrunk_se),
                                shrunk_est - z_ci * shrunk_se, NA_real_)
  tbl[[shrunk_upper]] <- ifelse(is.finite(shrunk_est) & is.finite(shrunk_se),
                                shrunk_est + z_ci * shrunk_se, NA_real_)
  tbl$CI_Level <- ci_level
  tbl
}

#' @keywords internal
#' @noRd
.draw_shrinkage_plot <- function(data_list, style, title,
                                  show_ci = FALSE,
                                  ci_level = 0.95) {
  tbl <- data_list$table
  if (nrow(tbl) == 0L) {
    graphics::plot.new()
    graphics::title(main = title %||% "Shrinkage plot")
    graphics::text(0.5, 0.5,
                   "No shrinkage applied.\nFit with facet_shrinkage = 'empirical_bayes'.")
    return(invisible(NULL))
  }

  # Order: group by facet, then by magnitude of shrinkage (largest at
  # bottom so the visual attention falls on most-shrunken levels).
  tbl <- tbl[order(tbl$Facet, -abs(tbl$ShrinkageFactor)), , drop = FALSE]
  tbl$Row <- seq_len(nrow(tbl))
  labels <- paste0(tbl$Facet, " / ", tbl$Level)

  # x-range with margin for error bars.
  x_vals <- c(tbl$Estimate, tbl$ShrunkEstimate)
  se_vals <- c(tbl$SE, tbl$ShrunkSE)
  x_vals <- x_vals[is.finite(x_vals)]
  if (isTRUE(show_ci) && any(is.finite(se_vals))) {
    ci_level <- .validate_shrinkage_ci_level(ci_level)
    z_ci <- stats::qnorm(1 - (1 - ci_level) / 2)
    x_vals <- c(x_vals, tbl$Estimate - z_ci * tbl$SE,
                tbl$Estimate + z_ci * tbl$SE,
                tbl$ShrunkEstimate - z_ci * tbl$ShrunkSE,
                tbl$ShrunkEstimate + z_ci * tbl$ShrunkSE)
    x_vals <- x_vals[is.finite(x_vals)]
  }
  if (length(x_vals) == 0L) x_vals <- c(-1, 1)
  xr <- range(x_vals, na.rm = TRUE)
  if (diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
  xlim <- xr + c(-0.1, 0.1) * diff(xr)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4, max(7, min(14, max(nchar(labels)) * 0.5)), 3, 1))
  graphics::plot(
    x = tbl$Estimate,
    y = tbl$Row,
    xlim = xlim,
    ylim = c(nrow(tbl) + 0.5, 0.5),
    yaxt = "n",
    xlab = "logit",
    ylab = "",
    pch = 16,
    col = style$accent_primary,
    main = title,
    cex = 1.1
  )
  graphics::abline(v = 0, lty = 2, col = style$neutral)
  graphics::axis(2, at = seq_len(nrow(tbl)), labels = labels,
                 las = 1, cex.axis = 0.8)

  # Arrows from original to shrunk (shrinkage direction).
  for (i in seq_len(nrow(tbl))) {
    if (is.finite(tbl$Estimate[i]) && is.finite(tbl$ShrunkEstimate[i]) &&
        abs(tbl$Estimate[i] - tbl$ShrunkEstimate[i]) > 1e-6) {
      graphics::arrows(
        x0 = tbl$Estimate[i], y0 = i,
        x1 = tbl$ShrunkEstimate[i], y1 = i,
        length = 0.06, angle = 20,
        col = style$neutral, lwd = 1
      )
    }
  }

  # Shrunk estimate on top (open circles) so both are visible.
  graphics::points(
    x = tbl$ShrunkEstimate, y = tbl$Row,
    pch = 21, bg = "white",
    col = style$accent_tertiary, cex = 1.1
  )

  # Optional CI error bars for both.
  if (isTRUE(show_ci)) {
    ci_level <- .validate_shrinkage_ci_level(ci_level)
    z_ci <- stats::qnorm(1 - (1 - ci_level) / 2)
    for (i in seq_len(nrow(tbl))) {
      if (is.finite(tbl$Estimate[i]) && is.finite(tbl$SE[i])) {
        graphics::segments(
          x0 = tbl$Estimate[i] - z_ci * tbl$SE[i], y0 = i - 0.12,
          x1 = tbl$Estimate[i] + z_ci * tbl$SE[i], y1 = i - 0.12,
          col = style$accent_primary
        )
      }
      if (is.finite(tbl$ShrunkEstimate[i]) && is.finite(tbl$ShrunkSE[i])) {
        graphics::segments(
          x0 = tbl$ShrunkEstimate[i] - z_ci * tbl$ShrunkSE[i], y0 = i + 0.12,
          x1 = tbl$ShrunkEstimate[i] + z_ci * tbl$ShrunkSE[i], y1 = i + 0.12,
          col = style$accent_tertiary
        )
      }
    }
  }

  # Horizontal band separators between facets.
  facet_boundaries <- cumsum(table(tbl$Facet)[unique(tbl$Facet)])
  for (b in utils::head(facet_boundaries, -1L)) {
    graphics::abline(h = b + 0.5, lty = 3, col = style$grid)
  }

  legend_labels <- c("Original", "Shrunk", "Shrinkage direction")
  legend_pch <- c(16, 21, NA)
  legend_lty <- c(NA, NA, 1)
  legend_col <- c(style$accent_primary, style$accent_tertiary, style$neutral)
  legend_bg <- c(NA, "white", NA)
  if (isTRUE(show_ci)) {
    legend_labels <- c(legend_labels, sprintf("%g%% CI", round(100 * ci_level)))
    legend_pch <- c(legend_pch, NA)
    legend_lty <- c(legend_lty, 1)
    legend_col <- c(legend_col, style$accent_primary)
    legend_bg <- c(legend_bg, NA)
  }
  graphics::legend(
    "topright",
    legend = legend_labels,
    pch = legend_pch,
    lty = legend_lty,
    col = legend_col,
    pt.bg = legend_bg,
    bg = "white",
    cex = 0.85
  )
  invisible(NULL)
}

#' Extract the shrinkage report from a fitted mfrm_fit
#'
#' Lightweight accessor that returns the per-facet empirical-Bayes
#' shrinkage table stored on a fit when `facet_shrinkage != "none"`.
#' Returns `NULL` (with a message) when no shrinkage has been applied
#' so callers can probe without error.
#'
#' @param fit An `mfrm_fit` object.
#' @return A data.frame with one row per facet (and optionally
#'   `"Person"`) or `NULL` when shrinkage has not been applied.
#' @seealso [apply_empirical_bayes_shrinkage()], [fit_mfrm()].
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30,
#'                 facet_shrinkage = "empirical_bayes")
#' shrinkage_report(fit)
#' @export
shrinkage_report <- function(fit) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit from fit_mfrm().", call. = FALSE)
  }
  rep_tbl <- fit$shrinkage_report
  if (is.null(rep_tbl) || !is.data.frame(rep_tbl) || nrow(rep_tbl) == 0L) {
    message("No shrinkage applied; ",
            "call fit_mfrm(..., facet_shrinkage = 'empirical_bayes') or ",
            "apply_empirical_bayes_shrinkage(fit) first.")
    return(NULL)
  }
  rep_tbl
}
