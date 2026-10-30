# ==============================================================================
# Person fit indices: lz (Drasgow et al. 1985)
# ==============================================================================
#
# `compute_person_fit_indices()` extends the Infit / Outfit / ZSTD
# columns that `diagnose_mfrm()` already returns with a person-level
# standardized log-likelihood statistic:
#
# - lz (Drasgow, Levine & Williams, 1985): standardized log-likelihood
#   under the fitted model. Computed in its proper polytomous form
#   using the per-observation category probability `Pr(X = x | theta)`
#   from the model (not the Gaussian-residual approximation used in
#   earlier mfrmr releases). The centering term is the per-item
#   negentropy E[log P_k] = sum_k P_k log P_k and the variance term is
#   sum_k P_k (log P_k)^2 - (sum_k P_k log P_k)^2, summed across the
#   person's observations. Asymptotically standard normal under the
#   conditional-independence assumption when person ability is known.
# - lz_star: Snijders (2001) correction for estimated person parameters,
#   computed for JML/fixed-effect person estimates conditional on the fitted
#   non-person calibration. The polytomous extension uses
#   w_tilde_k = log(P_k) - c * d log(P_k) / d theta. For MML EAP person
#   scores, the column is returned as NA because EAP is not the ML/MAP/WLE
#   estimating-equation setting covered by Snijders' correction.
#
# Note: earlier mfrmr versions reported an `ECI4` column whose
# implementation was the standardized chi-square (sum StdSq - n) /
# sqrt(2 n) rather than the Tatsuoka & Tatsuoka (1983) extended-caution
# index. That column was a duplicate of the existing `OutfitZSTD`
# statistic in `diagnose_mfrm()$measures` under the linear (Smith
# whexact) approximation and has been removed in 0.2.0; use
# `OutfitZSTD` instead.

#' Person fit indices: lz and Snijders-corrected lz*
#'
#' Computes person-level fit statistics for an MFRM bundle, extending
#' the Infit / Outfit / ZSTD columns that `diagnose_mfrm()$measures`
#' already exposes with the standardized log-likelihood `lz` and, when
#' justified by the person-estimation method, Snijders' `lz*`.
#'
#' @param diagnostics Output from [diagnose_mfrm()].
#' @param fit Optional `mfrm_fit` from [fit_mfrm()]. Required to decide
#'   whether the person estimates are JML/fixed-effect estimates for which
#'   the Snijders (2001) correction is computed. MML/EAP person scores
#'   return `NA` for `lz_star` with an explanatory status. Supplied fit and
#'   diagnostics must identify the same analysis.
#'
#' @return A data frame of class `mfrm_person_fit_indices` with one row per
#' Person and columns:
#' \describe{
#'   \item{`Person`}{Person ID.}
#'   \item{`TotalN`, `N`, `UnavailableN`}{Total responses, responses with
#'     valid probability and log-likelihood moments, and unavailable responses.
#'     An incomplete response set does not produce a person-fit statistic.}
#'   \item{`LogLik`}{Sum of log P(X = x | theta) under the fitted
#'     model. Computed from the per-observation category probability
#'     `PrObserved` (the model probability of the observed category),
#'     not from a Gaussian residual approximation.}
#'   \item{`lz`}{Drasgow et al. (1985) standardized log-likelihood,
#'     in its polytomous form, treating each retained response once regardless
#'     of calibration weights. `lz_status` records computational availability.}
#'   \item{`lz_star`}{Snijders-corrected `lz*` when the source fit used
#'     JML/fixed-effect person estimates, conditioning on the fitted
#'     non-person calibration, and the diagnostics include the required
#'     derivative terms. Numerical convergence, a finite person estimate and
#'     unit observation weights are required; otherwise `NA`.}
#'   \item{`lz_star_status`}{Status string for `lz_star`, such as
#'     `"computed_jml_conditional_calibration"`, `"fit_required"`,
#'     `"not_applicable_eap"`, or `"insufficient_information"`.}
#'   \item{`lz_star_c`}{Estimated Snijders projection coefficient `c_n`
#'     for each person, when available.}
#'   \item{`lz_star_variance`}{Corrected variance denominator used for
#'     `lz_star`, when available.}
#'   \item{`lz_flag_5pct`, `lz_flag_1pct`}{Logical flags for practical
#'     two-sided `lz` thresholds of `|z| > 1.96` and `|z| > 2.58`.}
#'   \item{`lz_star_flag_5pct`, `lz_star_flag_1pct`}{The same flags for
#'     `lz_star`. All index flags, including `ReportFlag`, remain `NA` when
#'     the corresponding statistic is unavailable.}
#'   \item{`ReportIndex`, `ReportValue`, `ReportFlagLevel`,
#'     `ReportFlag`, `ReviewStatus`, `ReviewReason`, `ReportCaveat`}{Compact
#'     reporting columns. `ReportIndex` prefers `lz_star` when the Snijders
#'     correction was computed; otherwise it falls back to `lz` with an
#'     explicit caveat.}
#' }
#'
#' The `5pct` and `1pct` flag names refer to standard-normal reference cutoffs:
#' |lz| > 1.96 and |lz| > 2.58. Treat them as screening thresholds, not a
#' guarantee of those false-positive rates for fitted person scores.
#' Read `ReportCaveat` and `lz_star_status` with every flagged result. When
#' `lz_star_status == "computed_jml_conditional_calibration"`, `lz_star`
#' applies Snijders' estimated-ability correction for JML person estimates,
#' conditional on the fitted non-person parameters. This does not propagate
#' non-person calibration uncertainty. For MML/EAP person scores, use `lz`
#' with its documented caveat rather than treating EAP scores as if they
#' satisfied the Snijders estimating equation. Positive and negative flags
#' indicate more and less predictable patterns, respectively; neither alone
#' establishes invalid responses or justifies exclusion. The implemented
#' correction does not cover non-unit observation weights, even if normalized
#' to mean one. No multiple-person error-rate control is provided.
#'
#' Small positive probabilities are used without a lower floor. Zero, invalid
#' or missing observed probabilities make the person statistic unavailable.
#' Missing derivative information withholds the correction for the whole person.
#' With the matching fit supplied, older probability-moment diagnostics are
#' refreshed from that fit without re-estimation. Otherwise recompute
#' diagnostics first. Older saved person-fit tables/summaries must be rebuilt.
#'
#' Note: this implementation reads the model category probabilities
#' directly from the diagnostics bundle. Earlier mfrmr releases used
#' a Gaussian-residual approximation
#' \eqn{\log P(X = x) \approx -\tfrac{1}{2}(R^2/V) - \tfrac{1}{2}\log(2\pi V)}
#' as a stand-in for \eqn{\log P}, which overstated the per-item
#' variance of \eqn{\log P} for polytomous items, shrinking the
#' reported `lz` toward zero. Numerical `lz` values are therefore
#' not directly comparable across mfrmr releases; treat the values
#' returned here as the polytomous statistic and re-evaluate any
#' historical `|lz| > 1.96` flagging that was based on the earlier
#' approximation.
#'
#' @section References:
#' - Drasgow, F., Levine, M. V., & Williams, E. A. (1985).
#'   Appropriateness measurement with polychotomous item response
#'   models and standardized indices.
#'   *British Journal of Mathematical and Statistical Psychology, 38*(1), 67-86.
#' - Snijders, T. A. B. (2001). Asymptotic null distribution of
#'   person fit statistics with estimated person parameter.
#'   *Psychometrika, 66*(3), 331-342.
#' - Magis, D., Raiche, G., & Beland, S. (2012). A didactic
#'   presentation of Snijders's lz* index of person fit with emphasis
#'   on response model selection and ability estimation.
#'   *Journal of Educational and Behavioral Statistics, 37*(1), 57-81.
#' - Sinharay, S. (2016). Asymptotically correct standardization of
#'   person-fit statistics beyond dichotomous items.
#'   *Psychometrika, 81*(4), 992-1013.
#'
#' @seealso [diagnose_mfrm()]
#' @examples
#' library(mfrmr)
#'
#' # Saved synthetic MML/RSM analysis of example_operational.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "operational-mml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#'
#' # Use the diagnostics saved with this fit
#' diagnostics <- example_fit$diagnostics
#'
#' # Screen for unusual person response patterns
#' person_fit <- compute_person_fit_indices(diagnostics, fit = fit)
#' head(person_fit[, c("Person", "N", "ReportIndex", "ReportValue", "ReviewStatus")])
#'
#' # Inspect flagged persons and the limitations of their reported index
#' subset(person_fit, ReportFlag,
#'        c("Person", "ReportValue", "ReviewReason", "ReportCaveat"))
#' # With MML/EAP scores, lz_star is unavailable; lz is an uncorrected screening index
#' # A flag alone does not establish invalid responses or justify excluding a person
#' @export
#' @order 1
compute_person_fit_indices <- function(diagnostics, fit = NULL) {
  if (is.null(diagnostics) || !is.list(diagnostics) ||
      is.null(diagnostics$obs)) {
    stop("`diagnostics` must be a non-empty `mfrm_diagnostics` bundle.",
         call. = FALSE)
  }
  if (!is.null(fit)) {
    mfrm_results_validate_diagnostics_identity(fit, diagnostics, helper = "compute_person_fit_indices()")
  }
  obs <- diagnostics$obs
  if (inherits(diagnostics, "mfrm_diagnostics") &&
      !identical(attr(obs, "person_fit_moments_version", exact = TRUE), 2L)) {
    if (is.null(fit)) {
      stop("Recreate these diagnostics with diagnose_mfrm(fit, ...) or supply the matching fit to refresh person-fit probabilities.", call. = FALSE)
    }
    obs <- compute_obs_table(fit)
  }
  obs <- as.data.frame(obs, stringsAsFactors = FALSE)
  if (!nrow(obs) || anyNA(obs$Person) || any(!nzchar(as.character(obs$Person)))) {
    stop("Person-fit calculations require observations with non-missing Person identifiers.", call. = FALSE)
  }
  needed <- c("Person", "PrObserved", "ItemEntropy", "ItemVarLogP")
  missing_cols <- setdiff(needed, names(obs))
  if (length(missing_cols) > 0L) {
    stop("`diagnostics$obs` is missing required columns: ",
         paste(missing_cols, collapse = ", "),
         ". This typically means the diagnostics bundle was generated ",
         "by an older mfrmr version; re-run diagnose_mfrm(fit, ...) ",
         "to populate the per-observation probability columns.",
         call. = FALSE)
  }
  obs$Person <- as.character(obs$Person)
  obs$PrObserved <- suppressWarnings(as.numeric(obs$PrObserved))
  obs$ItemEntropy <- suppressWarnings(as.numeric(obs$ItemEntropy))
  obs$ItemVarLogP <- suppressWarnings(as.numeric(obs$ItemVarLogP))

  # True Drasgow et al. (1985) lz computed from the polytomous category
  # probabilities. log_p is the model log-probability of the observed
  # category; e_logp is the per-item negentropy E[log P]; var_logp is
  # the per-item Var(log P). All three terms are summed across each
  # person's observations to form lz = (l - E[l]) / sqrt(Var[l]).
  valid_probability <- is.finite(obs$PrObserved) & obs$PrObserved > 0 & obs$PrObserved <= 1
  log_p <- rep(NA_real_, nrow(obs))
  log_p[valid_probability] <- log(obs$PrObserved[valid_probability])

  agg <- by(
    data.frame(log_p = log_p,
               e_logp = obs$ItemEntropy,
               var_logp = obs$ItemVarLogP,
               stringsAsFactors = FALSE),
    obs$Person,
    function(d) {
      ok <- is.finite(d$log_p) & is.finite(d$e_logp) &
            is.finite(d$var_logp) & d$var_logp >= 0
      total <- nrow(d)
      n <- sum(ok)
      if (n < total) {
        return(c(TotalN = total, N = n, UnavailableN = total - n, LogLik = NA_real_, lz = NA_real_))
      }
      ll <- sum(d$log_p)
      e_ll <- sum(d$e_logp)
      var_ll <- sum(d$var_logp)
      lz_val <- if (is.finite(var_ll) && var_ll > 0) (ll - e_ll) / sqrt(var_ll) else NA_real_
      c(TotalN = total, N = n, UnavailableN = 0, LogLik = ll, lz = lz_val)
    }
  )
  out <- do.call(rbind, lapply(agg, function(x) as.data.frame(t(x))))
  out$Person <- rownames(out)
  rownames(out) <- NULL
  out <- out[, c("Person", "TotalN", "N", "UnavailableN", "LogLik", "lz"), drop = FALSE]
  out$N <- as.integer(out$N)
  out$TotalN <- as.integer(out$TotalN)
  out$UnavailableN <- as.integer(out$UnavailableN)
  out$lz_status <- ifelse(out$UnavailableN > 0L, "incomplete_observation_terms",
                          ifelse(is.finite(out$lz), "computed", "insufficient_variance"))

  numerical_ready <- is.null(fit) || !fit$config$method %in% c("JML", "JMLE") ||
    identical(as.character(mfrmr_get_readiness_record(fit)$fit$NumericalState[1]), "ready")
  lz_star_tbl <- compute_snijders_lz_star(obs, out$Person, fit = if (numerical_ready) fit else NULL)
  if (!numerical_ready) lz_star_tbl$lz_star_status[] <- "fit_not_converged"
  lz_idx <- match(out$Person, lz_star_tbl$Person)
  out$lz_star <- lz_star_tbl$lz_star[lz_idx]
  out$lz_star_status <- lz_star_tbl$lz_star_status[lz_idx]
  out$lz_star_c <- lz_star_tbl$lz_star_c[lz_idx]
  out$lz_star_variance <- lz_star_tbl$lz_star_variance[lz_idx]
  out <- add_person_fit_reporting_columns(out)
  out <- out[, c("Person", "TotalN", "N", "UnavailableN", "LogLik", "lz", "lz_status", "lz_star", "lz_star_status",
                 "lz_star_c", "lz_star_variance",
                 "lz_flag_5pct", "lz_flag_1pct",
                 "lz_star_flag_5pct", "lz_star_flag_1pct",
                 "ReportIndex", "ReportValue", "ReportFlagLevel", "ReportFlag",
                 "ReviewStatus", "ReviewReason", "ReportCaveat"),
             drop = FALSE]
  attr(out, "person_fit_thresholds") <- person_fit_threshold_table()
  attr(out, "person_fit_calculation_version") <- 2L
  class(out) <- c("mfrm_person_fit_indices", class(out))
  out
}

person_fit_threshold_table <- function() {
  data.frame(
    Threshold = c("5pct", "1pct"),
    TwoSidedAlpha = c(0.05, 0.01),
    AbsZ = c(stats::qnorm(0.975), stats::qnorm(0.995)),
    Rule = c("|z| > 1.96", "|z| > 2.58"),
    stringsAsFactors = FALSE
  )
}

add_person_fit_reporting_columns <- function(out) {
  z_5pct <- stats::qnorm(0.975)
  z_1pct <- stats::qnorm(0.995)

  out$lz_flag_5pct <- ifelse(is.finite(out$lz), abs(out$lz) > z_5pct, NA)
  out$lz_flag_1pct <- ifelse(is.finite(out$lz), abs(out$lz) > z_1pct, NA)
  out$lz_star_flag_5pct <- ifelse(is.finite(out$lz_star), abs(out$lz_star) > z_5pct, NA)
  out$lz_star_flag_1pct <- ifelse(is.finite(out$lz_star), abs(out$lz_star) > z_1pct, NA)

  snijders_ready <- is.finite(out$lz_star) &
    identical_length_status(out$lz_star_status, "computed_jml_conditional_calibration")
  lz_ready <- !snijders_ready & is.finite(out$lz)

  out$ReportIndex <- ifelse(snijders_ready, "lz_star",
                            ifelse(lz_ready, "lz", "none"))
  out$ReportValue <- ifelse(snijders_ready, out$lz_star,
                            ifelse(lz_ready, out$lz, NA_real_))
  abs_report <- abs(out$ReportValue)
  out$ReportFlagLevel <- ifelse(
    !is.finite(abs_report), "not_available",
    ifelse(abs_report > z_1pct, "1pct",
           ifelse(abs_report > z_5pct, "5pct", "none"))
  )
  out$ReportFlag <- ifelse(is.finite(out$ReportValue), out$ReportFlagLevel %in% c("5pct", "1pct"), NA)
  out$ReviewStatus <- ifelse(
    out$ReportFlagLevel == "1pct", "review_1pct",
    ifelse(out$ReportFlagLevel == "5pct", "review_5pct",
           ifelse(out$ReportFlagLevel == "none", "not_flagged",
                  "not_available"))
  )

  threshold_reason <- ifelse(
    out$ReportFlagLevel == "1pct",
    paste0(ifelse(out$ReportIndex == "lz_star", "lz*", "lz"), " exceeds the absolute normal-reference cutoff 2.58."),
    ifelse(
      out$ReportFlagLevel == "5pct",
      paste0(ifelse(out$ReportIndex == "lz_star", "lz*", "lz"), " exceeds the absolute normal-reference cutoff 1.96."),
      ifelse(
        out$ReportFlagLevel == "none",
        "No report-level flag under the practical two-sided thresholds.",
        "No finite person-fit statistic is available for report-level flagging."
      )
    )
  )
  out$ReviewReason <- threshold_reason

  status <- as.character(out$lz_star_status)
  reason <- person_fit_correction_note(status)
  out$ReportCaveat <- ifelse(snijders_ready,
    "lz* conditions on the fitted non-person calibration; calibration uncertainty is omitted.",
    paste(ifelse(lz_ready, "lz is an uncorrected, unweighted response-pattern screen.",
      "No person-fit index is available."), reason))
  incomplete <- out$lz_status == "incomplete_observation_terms"
  out$ReportCaveat[incomplete] <- paste(
    "Response probability or log-likelihood moment information is incomplete.",
    "A statistic is not calculated on a selected subset of responses."
  )
  out
}

person_fit_correction_note <- function(status) {
  notes <- c(
    computed_jml_conditional_calibration = "The ability correction conditions on the fitted calibration.",
    fit_required = "Supply the matching fit to assess the ability correction.",
    not_applicable_eap = "The implemented correction does not apply to MML/EAP person scores.",
    nonunit_weights = "The implemented ability correction does not cover non-unit observation weights.",
    nonfinite_person_estimate = "The person estimate is not finite, so the ability correction is unavailable.",
    fit_not_converged = "Numerical convergence is not established; the ability correction is unavailable.",
    incomplete_observation_terms = "The ability correction requires complete probability and derivative information for all responses.",
    diagnostics_missing_snijders_terms = "Recompute diagnostics to supply the required score derivatives.",
    insufficient_information = "There is insufficient information for the ability correction."
  )
  out <- unname(notes[as.character(status)])
  out[is.na(out)] <- "The ability correction is unavailable for this result."
  out
}

identical_length_status <- function(x, value) {
  x <- as.character(x)
  !is.na(x) & x == value
}

person_fit_count_table <- function(x, col, label = col) {
  if (!col %in% names(x)) {
    return(data.frame(
      Variable = label,
      Value = character(),
      Rows = integer(),
      Proportion = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  vals <- as.character(x[[col]])
  vals[is.na(vals) | !nzchar(vals)] <- "missing"
  tab <- sort(table(vals), decreasing = TRUE)
  total <- sum(tab)
  data.frame(
    Variable = label,
    Value = names(tab),
    Rows = as.integer(tab),
    Proportion = if (total > 0) as.numeric(tab) / total else NA_real_,
    stringsAsFactors = FALSE
  )
}

person_fit_review_priority <- function(status) {
  status <- as.character(status)
  out <- rep(0L, length(status))
  out[status == "review_5pct"] <- 1L
  out[status == "review_1pct"] <- 2L
  out[status == "not_available"] <- -1L
  out
}

person_fit_validate_table <- function(object) {
  if (!is.data.frame(object)) {
    stop("`object` must be a data frame returned by `compute_person_fit_indices()`.",
         call. = FALSE)
  }
  if (!identical(attr(object, "person_fit_calculation_version", exact = TRUE), 2L)) {
    stop("Recreate this person-fit result with compute_person_fit_indices(diagnostics, fit) to retain unavailable responses and flags.", call. = FALSE)
  }
  needed <- c(
    "Person", "TotalN", "N", "UnavailableN", "LogLik", "lz", "lz_star", "lz_star_status",
    "ReportIndex", "ReportValue", "ReportFlagLevel", "ReportFlag",
    "ReviewStatus", "ReviewReason", "ReportCaveat"
  )
  missing_cols <- setdiff(needed, names(object))
  if (length(missing_cols) > 0L) {
    stop("`object` is missing required person-fit columns: ",
         paste(missing_cols, collapse = ", "),
         ". Rebuild it with `compute_person_fit_indices()`.",
         call. = FALSE)
  }
  as.data.frame(object, stringsAsFactors = FALSE)
}

#' @rdname compute_person_fit_indices
#' @details `summary()` method: Summarize person-fit indices
#' @details `summary()` method: `summary()` for [compute_person_fit_indices()] output gives a compact,
#' report-ready reading order: first the number of persons and flags, then
#' status counts, then the highest-priority review rows. The summary keeps
#' `lz_star` availability visible so users do not silently treat uncorrected
#' `lz` as Snijders-corrected output.
#' @param object Output from [compute_person_fit_indices()].
#' @param digits Number of digits used when printing numeric columns.
#' @param top_n Number of review rows retained in `top_review`.
#' @param include_person If `TRUE`, print person identifiers in the review
#'   preview. The default retains the rows in the returned object but suppresses
#'   identifiers from console output.
#' @param ... Unused.
#' @return For `summary()`: An object of class `summary.mfrm_person_fit_indices` with:
#' - `overview`
#' - `status_summary`
#' - `report_index_summary`
#' - `lz_star_status_summary`
#' - `top_review`
#' - `caveats`
#' - `thresholds`
#' - `reporting_map`
#' - `notes`
#' @export
summary.mfrm_person_fit_indices <- function(object,
                                            digits = 3,
                                            top_n = 10,
                                            include_person = FALSE,
                                            ...) {
  tbl <- person_fit_validate_table(object)
  digits <- max(0L, as.integer(digits))
  top_n <- max(1L, as.integer(top_n))

  report_flag <- as.logical(tbl$ReportFlag)
  report_flag[is.na(report_flag)] <- FALSE
  status <- as.character(tbl$ReviewStatus)
  report_index <- as.character(tbl$ReportIndex)
  lz_star_status <- as.character(tbl$lz_star_status)
  report_value <- suppressWarnings(as.numeric(tbl$ReportValue))

  overview <- data.frame(
    Persons = nrow(tbl),
    ReportableRows = sum(is.finite(report_value)),
    ReportFlaggedRows = if (any(is.finite(report_value))) sum(report_flag) else NA_integer_,
    Review1PctRows = sum(status == "review_1pct", na.rm = TRUE),
    Review5PctRows = sum(status == "review_5pct", na.rm = TRUE),
    NotFlaggedRows = sum(status == "not_flagged", na.rm = TRUE),
    NotAvailableRows = sum(status == "not_available", na.rm = TRUE),
    SnijdersRows = sum(report_index == "lz_star", na.rm = TRUE),
    LzFallbackRows = sum(report_index == "lz", na.rm = TRUE),
    MissingReportIndexRows = sum(report_index == "none", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  overview$FlagRate <- if (overview$ReportableRows > 0) {
    overview$ReportFlaggedRows / overview$ReportableRows
  } else {
    NA_real_
  }

  status_summary <- person_fit_count_table(tbl, "ReviewStatus", "ReviewStatus")
  report_index_summary <- person_fit_count_table(tbl, "ReportIndex", "ReportIndex")
  lz_star_status_summary <- person_fit_count_table(tbl, "lz_star_status", "lz_star_status")

  review_tbl <- tbl
  review_tbl$.priority <- person_fit_review_priority(review_tbl$ReviewStatus)
  review_tbl$.abs_report <- abs(suppressWarnings(as.numeric(review_tbl$ReportValue)))
  review_tbl$.abs_report[!is.finite(review_tbl$.abs_report)] <- -Inf
  review_tbl <- review_tbl[order(-review_tbl$.priority, -review_tbl$.abs_report,
                                 as.character(review_tbl$Person)), ,
                           drop = FALSE]
  review_tbl <- review_tbl[review_tbl$.priority > 0 |
                             as.logical(review_tbl$ReportFlag %in% TRUE), ,
                           drop = FALSE]
  if (nrow(review_tbl) == 0L) {
    review_tbl <- tbl[order(-abs(report_value), as.character(tbl$Person)), ,
                      drop = FALSE]
  }
  keep_cols <- intersect(
    c("Person", "N", "ReportIndex", "ReportValue", "ReportFlagLevel",
      "ReviewStatus", "ReviewReason", "lz", "lz_star", "lz_star_status",
      "ReportCaveat"),
    names(review_tbl)
  )
  top_review <- utils::head(review_tbl[, keep_cols, drop = FALSE], n = top_n)
  top_review[] <- lapply(top_review, function(x) {
    if (is.numeric(x)) round(x, digits = digits) else x
  })

  caveats <- unique(tbl[, intersect(c("ReportIndex", "lz_star_status", "ReportCaveat"),
                                    names(tbl)),
                        drop = FALSE])
  caveats <- caveats[order(as.character(caveats$ReportIndex),
                           as.character(caveats$lz_star_status)), ,
                     drop = FALSE]

  thresholds <- attr(object, "person_fit_thresholds", exact = TRUE)
  if (!is.data.frame(thresholds) || nrow(thresholds) == 0L) {
    thresholds <- person_fit_threshold_table()
  }

  reporting_map <- data.frame(
    Area = c(
      "First read",
      "Flagged person rows",
      "Response-level follow-up",
      "Plot/data handoff"
    ),
    CoveredHere = c("yes", "yes", "no", "partial"),
    CompanionOutput = c(
      "overview / status_summary / lz_star_status_summary",
      "top_review",
      "build_misfit_casebook() / unexpected_response_table()",
      "plot(..., type = \"person_fit\", fit_index = \"loglik\", draw = FALSE) / plot_data()"
    ),
    stringsAsFactors = FALSE
  )

  notes <- c(
    "The reported index uses lz* when its conditional ability correction is available; otherwise it uses uncorrected lz where possible.",
    "The two-sided normal cutoffs are screening references; neither individual nor multiple-person false-positive rates are calibrated for these fitted scores.",
    "Person-fit flags are screening evidence. Review response-level evidence before making substantive claims about a person."
  )

  out <- list(
    overview = overview,
    status_summary = status_summary,
    report_index_summary = report_index_summary,
    lz_star_status_summary = lz_star_status_summary,
    top_review = top_review,
    caveats = caveats,
    thresholds = thresholds,
    reporting_map = reporting_map,
    notes = notes,
    digits = digits,
    include_person = isTRUE(include_person),
    calculation_version = 2L
  )
  class(out) <- "summary.mfrm_person_fit_indices"
  out
}

#' @export
print.summary.mfrm_person_fit_indices <- function(x, ...) {
  if (!identical(x$calculation_version, 2L)) {
    stop("Recreate this summary from compute_person_fit_indices(diagnostics, fit).", call. = FALSE)
  }
  cat("Person-Fit Summary\n")
  overview <- x$overview[, c("Persons", "ReportableRows", "NotAvailableRows", "ReportFlaggedRows", "FlagRate"), drop = FALSE]
  names(overview) <- c("Persons", "Available", "Unavailable", "Flagged", "Flag rate")
  print_bundle_section("Overview", overview, digits = x$digits)
  cat(sprintf("Ability-corrected lz*: %d; uncorrected lz: %d.\n",
              x$overview$SnijdersRows, x$overview$LzFallbackRows))
  if (isTRUE(x$include_person) && nrow(x$top_review)) {
    rows <- x$top_review[, c("Person", "N", "ReportIndex", "ReportValue", "ReviewReason"), drop = FALSE]
    rows$ReportIndex <- ifelse(rows$ReportIndex == "lz_star", "lz*",
                              ifelse(rows$ReportIndex == "lz", "lz", "Unavailable"))
    names(rows) <- c("Person", "Responses", "Index", "Value", "Review")
    print_bundle_section("Person review", rows, digits = x$digits)
  } else {
    cat("Person identifiers are omitted; use include_person = TRUE to display the review rows.\n")
  }
  print_bullet_section("Interpretation", unique(c(x$caveats$ReportCaveat, x$notes)))
  invisible(x)
}

#' @export
print.mfrm_person_fit_indices <- function(x, ...) {
  if (all(c("Person", "N", "LogLik", "lz_star", "ReviewStatus", "ReportValue", "ReportCaveat") %in% names(x))) {
    print(summary(x, ...))
  } else {
    print(as.data.frame(x), ...)
  }
  invisible(x)
}

compute_snijders_lz_star <- function(obs, persons, fit = NULL) {
  empty <- data.frame(
    Person = as.character(persons),
    lz_star = NA_real_,
    lz_star_status = "fit_required",
    lz_star_c = NA_real_,
    lz_star_variance = NA_real_,
    stringsAsFactors = FALSE
  )
  if (is.null(fit)) {
    return(empty)
  }
  if (!inherits(fit, "mfrm_fit")) {
    empty$lz_star_status <- "invalid_fit"
    return(empty)
  }

  method <- toupper(as.character(
    fit$config$method %||% fit$summary$Method[1] %||% NA_character_
  )[1])
  if (identical(method, "JMLE")) method <- "JML"
  if (!identical(method, "JML")) {
    empty$lz_star_status <- "not_applicable_eap"
    return(empty)
  }

  weights <- get_weights(obs)
  if (any(!is.finite(weights) | weights != 1)) {
    empty$lz_star_status <- "nonunit_weights"
    return(empty)
  }

  needed <- c("ItemLogPScoreCov", "ScoreInformation")
  missing_cols <- setdiff(needed, names(obs))
  if (length(missing_cols) > 0L) {
    empty$lz_star_status <- "diagnostics_missing_snijders_terms"
    return(empty)
  }

  obs$ItemLogPScoreCov <- suppressWarnings(as.numeric(obs$ItemLogPScoreCov))
  obs$ScoreInformation <- suppressWarnings(as.numeric(obs$ScoreInformation))
  if ("ObservedScoreDerivative" %in% names(obs)) {
    obs$ObservedScoreDerivative <- suppressWarnings(as.numeric(obs$ObservedScoreDerivative))
  } else if (all(c("ScoreSlope", "Observed", "Expected") %in% names(obs))) {
    obs$ObservedScoreDerivative <- suppressWarnings(
      as.numeric(obs$ScoreSlope) * (as.numeric(obs$Observed) - as.numeric(obs$Expected))
    )
  } else {
    empty$lz_star_status <- "diagnostics_missing_snijders_terms"
    return(empty)
  }

  valid_probability <- is.finite(obs$PrObserved) & obs$PrObserved > 0 & obs$PrObserved <= 1
  log_p <- rep(NA_real_, nrow(obs))
  log_p[valid_probability] <- log(obs$PrObserved[valid_probability])
  obs$.WCentered <- log_p - obs$ItemEntropy
  obs$Person <- as.character(obs$Person)

  rows <- split(obs, obs$Person)
  out <- empty
  for (i in seq_along(persons)) {
    person <- as.character(persons[i])
    d <- rows[[person]]
    if (is.null(d) || nrow(d) == 0L) {
      out$lz_star_status[i] <- "insufficient_information"
      next
    }
    ok <- is.finite(d$.WCentered) &
      is.finite(d$ItemVarLogP) & d$ItemVarLogP >= 0 &
      is.finite(d$ItemLogPScoreCov) &
      is.finite(d$ScoreInformation) & d$ScoreInformation >= 0 &
      is.finite(d$ObservedScoreDerivative)
    if (any(!ok)) {
      out$lz_star_status[i] <- "incomplete_observation_terms"
      next
    }

    if ("PersonMeasure" %in% names(d) && any(!is.finite(d$PersonMeasure))) {
      out$lz_star_status[i] <- "nonfinite_person_estimate"
      next
    }
    info_total <- sum(d$ScoreInformation)
    cov_total <- sum(d$ItemLogPScoreCov)
    var_logp_total <- sum(d$ItemVarLogP)
    if (!is.finite(info_total) || info_total <= 0 ||
        !is.finite(var_logp_total) || var_logp_total <= 0) {
      out$lz_star_status[i] <- "insufficient_information"
      next
    }

    c_n <- cov_total / info_total
    corrected_var <- var_logp_total - (cov_total^2 / info_total)
    if (!is.finite(corrected_var) || corrected_var <= sqrt(.Machine$double.eps)) {
      out$lz_star_status[i] <- "insufficient_information"
      out$lz_star_c[i] <- c_n
      out$lz_star_variance[i] <- corrected_var
      next
    }

    centered_loglik <- sum(d$.WCentered)
    score_sum <- sum(d$ObservedScoreDerivative)
    out$lz_star[i] <- (centered_loglik - c_n * score_sum) / sqrt(corrected_var)
    out$lz_star_status[i] <- "computed_jml_conditional_calibration"
    out$lz_star_c[i] <- c_n
    out$lz_star_variance[i] <- corrected_var
  }
  out
}
