# ==============================================================================
# Empirical-Bayes and Laplace shrinkage for small-N facets
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
# - `.compute_facet_shrinkage()` implements a zero-centered normal-model
#   moment adjustment, not a James-Stein risk guarantee.
# - `apply_empirical_bayes_shrinkage()` is the public wrapper that takes an
#   existing `mfrm_fit` and returns a fit augmented with shrunk columns and
#   a `shrinkage_report` table.
# - `fit_mfrm(..., facet_shrinkage = "empirical_bayes")` wires the same
#   helper into the primary entry point so users can opt in up front.
#
# The default `facet_shrinkage = "none"` keeps the unshrunk behavior
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


shrinkage_interpretation <- function() {
  "Descriptive zero-centered adjustment; plug-in SEs/bands omit prior-variance uncertainty and cross-level covariance. Zero SE after full pooling is not perfect precision."
}

validate_shrinkage_output <- function(fit) {
  report <- fit$shrinkage_report
  if (is.data.frame(report) && nrow(report) > 0L &&
      (!all(c("Interpretation", "CalculationVersion") %in% names(report)) ||
       !all(report$CalculationVersion %in% 2L))) {
    stop("Reapply apply_empirical_bayes_shrinkage() with the original prior settings to refresh this saved report, then regenerate plots and exports; no MFRM refit is needed. Use shrink_person = FALSE unless the fit uses ordinary JML.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.compute_facet_shrinkage <- function(estimates, ses,
                                     method = c("empirical_bayes"),
                                     prior_sd = NULL,
                                     min_levels = 3L) {
  # Zero-centered moment adjustment; tau2 is held fixed in the plug-in SE.
  # Identification/anchor constraints and sampling covariance are not modeled.
  # tau2 = max(0, mean(estimates^2) - mean(se^2)) on complete pairs.
  # B = se^2 / (tau2 + se^2); adjusted estimate = (1-B) * estimate.
  # Plug-in SE = sqrt((1-B) * se^2).
  #
  # Returns a list with parallel shrunk vectors so downstream callers can
  # merge them back onto the facet table by row order.
  method <- match.arg(method)
  if (!is.null(prior_sd) && (!is.numeric(prior_sd) || length(prior_sd) != 1L ||
      !is.finite(prior_sd) || prior_sd < 0 || !is.finite(prior_sd^2) ||
      (prior_sd > 0 && prior_sd^2 == 0))) {
    stop("`facet_prior_sd` must be NULL or a non-negative finite scalar whose square is representable (positive when the SD is positive).", call. = FALSE)
  }

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
    retained_weights = rep(NA_real_, K),
    tau2 = NA_real_,
    mean_se2 = NA_real_,
    n_levels = K,
    n_levels_used = 0L,
    method = method,
    prior_sd_source = if (is.null(prior_sd)) "empirical" else "user",
    note = NA_character_
  )

  valid <- is.finite(estimates) & is.finite(ses) & ses > 0
  K_use <- sum(valid)
  out_template$n_levels_used <- K_use
  if (K_use < min_levels) {
    out_template$note <- paste0(
      "Fewer than ", min_levels, " valid levels; shrinkage not applied."
    )
    # Pass through unchanged.
    out_template$shrunk_estimates <- estimates
    out_template$shrunk_ses <- ses
    out_template$shrinkage_factors <- ifelse(valid, 0, NA_real_)
    out_template$retained_weights <- ifelse(valid, 1, NA_real_)
    return(out_template)
  }

  # Scale before squaring/summing: a finite mean square need not have a
  # representable sum of squares, or even representable individual squares.
  mean_square <- function(x) {
    scale <- max(abs(x))
    if (scale == 0) return(0)
    value <- (scale * sqrt(mean((x / scale)^2)))^2
    if (value == 0) NA_real_ else value
  }
  mean_se2 <- mean_square(ses[valid])
  # A zero-centered moment heuristic on complete estimate/SE pairs.
  # Identification constraints do not supply a known population mean or
  # independent sampling errors; their covariance is not modeled here.
  raw_var <- mean_square(estimates[valid])

  # Pull prior variance either from the supplied value or from the
  # method-of-moments estimator. The typical choice is MoM.
  if (!is.null(prior_sd)) {
    tau2 <- as.numeric(prior_sd)^2
    out_template$prior_sd_source <- "user"
  } else {
    tau2 <- max(0, raw_var - mean_se2)
    out_template$prior_sd_source <- "empirical"
  }

  # Invalid pairs pass through with unavailable factors.
  B <- rep(NA_real_, K)
  retained <- rep(NA_real_, K)
  shrunk_est <- estimates
  shrunk_se <- ses
  if (!is.finite(tau2)) {
    tau2 <- NA_real_
    out_template$note <- "Moment variance is outside numerical range; shrinkage not applied. Original estimates and SEs are retained."
  } else if (tau2 > 0) {
    tau <- sqrt(tau2)
    se <- ses[valid]
    ratio <- pmin(tau, se) / pmax(tau, se)
    denominator <- sqrt(1 + ratio^2)
    retained_root <- ifelse(se <= tau, 1 / denominator, ratio / denominator)
    B[valid] <- ifelse(se <= tau, ratio^2 / (1 + ratio^2), 1 / (1 + ratio^2))
    retained[valid] <- retained_root^2
    # Do not subtract a factor rounded to one, or multiply a tiny variance
    # weight by an overflowing squared SE.
    shrunk_est[valid] <- (estimates[valid] * retained_root) * retained_root
    shrunk_se[valid] <- pmin(tau, se) / denominator
  } else {
    # Full pooling under the selected zero-variance normal approximation.
    B[valid] <- 1
    retained[valid] <- 0
    shrunk_est[valid] <- 0
    shrunk_se[valid] <- 0
    out_template$note <- paste0(
      "The selected prior variance is zero; eligible levels are fully pooled at zero. ",
      "Zero plug-in SE does not establish perfect precision."
    )
  }
  if (!is.finite(mean_se2)) {
    mean_se2 <- NA_real_
    out_template$note <- paste(stats::na.omit(c(out_template$note,
      "Mean squared SE is outside numerical range and is unavailable.")), collapse = " ")
  }

  out_template$shrunk_estimates <- shrunk_est
  out_template$shrunk_ses <- shrunk_se
  out_template$shrinkage_factors <- B
  out_template$retained_weights <- retained
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

  if (!is.logical(shrink_person) || length(shrink_person) != 1L || is.na(shrink_person)) {
    stop("`shrink_person` must be a single logical value.", call. = FALSE)
  }
  if (isTRUE(shrink_person) &&
      (!identical(fit$config$method, "JML") || mfrm_has_jml_adjustment(fit))) {
    stop("`shrink_person = TRUE` requires ordinary JML estimates with sampling SEs. EAP posterior SDs are not sampling SEs. Use shrink_person = FALSE to adjust non-Person facets while retaining the original Person estimates.", call. = FALSE)
  }
  fit$config$facet_shrinkage <- as.character(method)
  fit$config$facet_prior_sd <- facet_prior_sd
  fit$config$shrinkage_settings <- list(
    facet_prior_sd = facet_prior_sd, shrink_person = shrink_person,
    applied_after_fit = FALSE
  )
  # Reapplication always starts from the original estimates and SEs. In
  # particular, turning Person shrinkage off must not retain old adjustments.
  person <- fit$facets$person
  if (is.data.frame(person)) {
    person[intersect(c("ShrunkEstimate", "ShrunkSE", "ShrinkageFactor"), names(person))] <- NULL
    fit$facets$person <- person
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
    effective_df <- if (any(is.finite(res$retained_weights))) {
      sum(res$retained_weights[is.finite(res$retained_weights)])
    } else NA_real_
    report_rows[[f]] <- data.frame(
      Facet = f,
      NLevels = res$n_levels,
      NLevelsUsed = res$n_levels_used,
      Tau2 = res$tau2,
      MeanSE2 = res$mean_se2,
      MeanShrinkage = if (any(is.finite(res$shrinkage_factors))) {
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

  # Optional descriptive adjustment of ordinary-JML Person estimates.
  # Posterior SDs are not observation-error SEs for this normal model.
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
        MeanShrinkage = if (any(is.finite(res_p$shrinkage_factors))) {
          mean(res_p$shrinkage_factors[is.finite(res_p$shrinkage_factors)])
        } else NA_real_,
        EffectiveDF = if (any(is.finite(res_p$retained_weights))) {
          sum(res_p$retained_weights[is.finite(res_p$retained_weights)])
        } else NA_real_,
        Method = res_p$method,
        PriorSource = res_p$prior_sd_source,
        Note = res_p$note %||% NA_character_,
        stringsAsFactors = FALSE
      )
      fit$shrinkage_report <- rbind(fit$shrinkage_report, person_report)
      rownames(fit$shrinkage_report) <- NULL
    }
  }

  fit$shrinkage_report$SupportsFormalInference <- FALSE
  fit$shrinkage_report$Interpretation <- shrinkage_interpretation()
  fit$shrinkage_report$CalculationVersion <- 2L

  # Summary-level surfacing. Downstream tables and narrative read these.
  if (!is.null(fit$summary) && nrow(fit$summary) > 0L) {
    fit$summary$FacetShrinkage <- as.character(method)
    fit$summary$FacetShrinkageTau2Mean <- if (nrow(fit$shrinkage_report) > 0L) {
      tau <- fit$shrinkage_report$Tau2
      if (any(is.finite(tau))) mean(tau[is.finite(tau)]) else NA_real_
    } else NA_real_
  }

  fit
}


#' Apply empirical-Bayes shrinkage to fitted non-person facet estimates
#'
#' Post-hoc normal-model shrinkage helper that augments an `mfrm_fit`
#' with descriptive empirical-Bayes adjustments for each non-person facet. The
#' shrinkage variance \eqn{\hat{\tau}^2} is estimated by method of
#' moments from the facet-level point estimates and their standard
#' errors:
#' \deqn{\hat{\tau}^2 = \max\!\left(0,
#'   \frac{1}{K}\sum_{j=1}^{K}\hat{\delta}_j^{2} -
#'   \overline{\mathrm{SE}^2}\right),}
#' where zero is the chosen shrinkage target and the moment is calculated on
#' complete finite estimate/positive-SE pairs. A sum-to-zero identification
#' constraint does not establish a known population mean or independent errors.
#' The shrinkage factor is
#' \eqn{B_j = \mathrm{SE}_j^2 / (\hat{\tau}^2 + \mathrm{SE}_j^2)}, with
#' \eqn{\hat{\delta}_j^{EB} = (1-B_j)\hat{\delta}_j} and
#' \eqn{\mathrm{ShrunkSE}_j = \sqrt{(1-B_j)\mathrm{SE}_j^2}}.
#'
#' `ShrunkSE` is a plug-in normal-model quantity conditional on the selected
#' prior variance. It omits uncertainty in that variance, cross-level covariance
#' and the fitted identification/anchor constraints. It is neither a calibrated
#' SE nor a guaranteed lower bound. Unequal shrinkage can break the original
#' sum-to-zero constraint. Original estimates, anchors and predictions are
#' unchanged; this helper does not refit a hierarchical MFRM.
#' The normal-normal formula is exact for an independent Gaussian observation
#' with a known Gaussian prior variance. Applying it to jointly calibrated
#' MFRM estimates is a descriptive approximation; a James-Stein improvement
#' in risk is not established for this adjustment.
#'
#' A zero estimated prior variance fully pools eligible estimates at zero and
#' gives zero plug-in SE; this does not establish perfect precision. Optional
#' plot whiskers are descriptive normal bands. Their legacy `CI` column names
#' do not confer confidence-interval coverage; `SupportsFormalInference` is
#' false. Invalid estimate/SE pairs retain their original values with unavailable
#' shrinkage factors. With fewer than three valid pairs no shrinkage is applied.
#' Unrepresentable moment variances also leave the original values unchanged,
#' with unavailable factors and an explanatory `Note` in the report.
#'
#' The report records complete-pair counts and mean shrinkage over those pairs.
#' `EffectiveDF` sums retained weights over eligible levels, conditional on the
#' chosen prior variance; it is not model degrees of freedom for testing or IC.
#' Reapply this helper to existing fits to refresh old reports; no MFRM refit
#' is needed. Regenerate previously saved plots and exports as well. Updated
#' reports carry `CalculationVersion = 2`: positive prior variances no longer
#' produce zero plug-in SE merely because a shrinkage factor rounded to one.
#' Reapplication replaces the previous adjustment; `shrink_person = FALSE`
#' removes any previous Person shrinkage columns. The original fitting inputs
#' remain recorded separately from the latest post-fit shrinkage settings.
#' [build_mfrm_replay_script()] replays the adjustment after fitting, preserving
#' whether diagnostic SEs were attached by [fit_mfrm()]. Manually edited SEs
#' or other table edits require their own reproducible editing steps.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()] with a non-empty
#'   `facets$others` table.
#' @param facet_prior_sd Optional numeric scalar. When supplied, the
#'   shrinkage variance is fixed at `facet_prior_sd^2` instead of being
#'   estimated from the data. Its square must be representable, and positive
#'   when the supplied SD is positive. Useful when a prior is elicited from
#'   expert knowledge or a previous fit.
#' @param shrink_person Logical. When `TRUE`, the same empirical-Bayes
#'   adjustment is also applied to ordinary-JML Person estimates with
#'   available sampling SEs. Default `FALSE`. Other estimators are refused:
#'   MML EAP estimates already use a population prior, and their posterior
#'   SDs are not sampling SEs for a second normal-model update.
#'
#' @return The same `mfrm_fit`, with augmented columns and a new
#'   `shrinkage_report` list entry, and with
#'   `fit$config$facet_shrinkage` set to `"empirical_bayes"`.
#'
#' @section Typical workflow:
#' 1. Fit the model as usual with `fit_mfrm()`.
#' 2. Call `apply_empirical_bayes_shrinkage(fit)` when small-N facets
#'    are present and zero-centered pooling is substantively defensible
#'    (see [facet_small_sample_review()]); small counts alone do not require it.
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
#' # Saved synthetic JML/RSM analysis of example_core.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' fit_eb <- apply_empirical_bayes_shrinkage(fit)
#' fit_eb$shrinkage_report
#' # Look for:
#' # - `Tau2` is the estimated between-level prior variance per facet.
#' #   `Tau2 = 0` fully pools eligible estimates at zero
#' #   (`MeanShrinkage = 1`). Zero plug-in SE is not perfect precision.
#' # - `MeanShrinkage` near 0 = little movement, near 1 = heavy pooling
#' #   toward the chosen zero target; inspect the SEs and prior variance.
#' # - `EffectiveDF` sums retained weights over eligible levels, conditional
#' #   on the chosen prior variance; it is not degrees of freedom for testing.
#' head(fit_eb$facets$others[, c("Facet", "Level", "Estimate",
#'                                "ShrunkEstimate", "ShrinkageFactor")])
#' # Look for: rows where `ShrinkageFactor` is large (close to 1) had
#' #   their estimates pulled most strongly toward the facet mean (0).
#' @export
apply_empirical_bayes_shrinkage <- function(fit,
                                            facet_prior_sd = NULL,
                                            shrink_person = FALSE) {
  stop_if_product_slopes(fit, "apply_empirical_bayes_shrinkage()")
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit from fit_mfrm().", call. = FALSE)
  }
  fit <- .apply_shrinkage_to_fit(
    fit = fit,
    method = "empirical_bayes",
    facet_prior_sd = facet_prior_sd,
    shrink_person = shrink_person
  )
  fit$config$shrinkage_settings$applied_after_fit <- TRUE
  fit
}


#' @keywords internal
#' @noRd
.build_shrinkage_plot_data <- function(fit) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`.build_shrinkage_plot_data()` requires an mfrm_fit.", call. = FALSE)
  }
  validate_shrinkage_output(fit)
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
  tbl$SupportsFormalInference <- rep(FALSE, nrow(tbl))
  tbl$Interpretation <- rep(shrinkage_interpretation(), nrow(tbl))
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
  z_ci <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
  raw_lower <- if (nzchar(raw_prefix)) paste0(raw_prefix, "CI_Lower") else "CI_Lower"
  raw_upper <- if (nzchar(raw_prefix)) paste0(raw_prefix, "CI_Upper") else "CI_Upper"
  shrunk_lower <- paste0(shrunk_prefix, "CI_Lower")
  shrunk_upper <- paste0(shrunk_prefix, "CI_Upper")

  raw_est <- if (raw_estimate_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[raw_estimate_col]])) else rep(NA_real_, nrow(tbl))
  raw_se <- if (raw_se_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[raw_se_col]])) else rep(NA_real_, nrow(tbl))
  shrunk_est <- if (shrunk_estimate_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[shrunk_estimate_col]])) else rep(NA_real_, nrow(tbl))
  shrunk_se <- if (shrunk_se_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[shrunk_se_col]])) else rep(NA_real_, nrow(tbl))

  tbl[[raw_lower]] <- ifelse(is.finite(raw_est) & is.finite(raw_se) & raw_se >= 0,
                             raw_est - z_ci * raw_se, NA_real_)
  tbl[[raw_upper]] <- ifelse(is.finite(raw_est) & is.finite(raw_se) & raw_se >= 0,
                             raw_est + z_ci * raw_se, NA_real_)
  tbl[[shrunk_lower]] <- ifelse(is.finite(shrunk_est) & is.finite(shrunk_se) & shrunk_se >= 0,
                                shrunk_est - z_ci * shrunk_se, NA_real_)
  tbl[[shrunk_upper]] <- ifelse(is.finite(shrunk_est) & is.finite(shrunk_se) & shrunk_se >= 0,
                                shrunk_est + z_ci * shrunk_se, NA_real_)
  tbl$CI_Level <- ci_level
  tbl$SupportsFormalInference <- rep(FALSE, nrow(tbl))
  tbl$IntervalInterpretation <- rep(shrinkage_interpretation(), nrow(tbl))
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
    z_ci <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
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

  old_par <- graphics::par()["mar"]
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
    main = "",
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
    z_ci <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
    for (i in seq_len(nrow(tbl))) {
      if (is.finite(tbl$Estimate[i]) && is.finite(tbl$SE[i]) && tbl$SE[i] >= 0) {
        graphics::segments(
          x0 = tbl$Estimate[i] - z_ci * tbl$SE[i], y0 = i - 0.12,
          x1 = tbl$Estimate[i] + z_ci * tbl$SE[i], y1 = i - 0.12,
          col = style$accent_primary
        )
      }
      if (is.finite(tbl$ShrunkEstimate[i]) && is.finite(tbl$ShrunkSE[i]) && tbl$ShrunkSE[i] >= 0) {
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
    legend_labels <- c(legend_labels, sprintf("%g%% plug-in band", round(100 * ci_level)))
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
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("shrinkage_report", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 300,
#'                 facet_shrinkage = "empirical_bayes")
#' shrinkage_report(fit)
#' }
#' @export
shrinkage_report <- function(fit) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit from fit_mfrm().", call. = FALSE)
  }
  validate_shrinkage_output(fit)
  rep_tbl <- fit$shrinkage_report
  if (is.null(rep_tbl) || !is.data.frame(rep_tbl) || nrow(rep_tbl) == 0L) {
    message("No shrinkage applied; ",
            "call fit_mfrm(..., facet_shrinkage = 'empirical_bayes') or ",
            "apply_empirical_bayes_shrinkage(fit) first.")
    return(NULL)
  }
  rep_tbl
}
