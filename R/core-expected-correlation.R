#' @section Expected point-measure correlation under ordinary JML:
#' For ordinary additive RSM/PCM JML without group-mean anchors,
#' `measures$PTMEAExpected` adds a conditional reference beside observed
#' `PTMEA`. Let `x` be the oriented generating measure, `e` the expected score,
#' `v` its conditional model variance, and `N` the sum of frequency weights
#' within the element. The reference is
#' \deqn{\frac{\operatorname{Cov}_w(x,e)}{
#' \sqrt{\operatorname{Var}_w(x)\{\operatorname{Var}_w(e)+(N-2)\bar v_w/N\}}}.}
#' Weighted moments use the same observations and normalized frequency weights.
#' `PTMEAExpectedAsymptotic` omits the `(N-2)/N` adjustment. The correction is
#' Linacre's published finite-count reporting approximation; it is not an exact
#' expectation of the sample correlation, a significance test or a confidence
#' interval. Calibration is held fixed. Different administered PCM scales are
#' retained when calculating each observation's expectation and variance.
#'
#' `PTMEAExpectedN`, `PTMEAExpectedFactor`, `PTMEAExpectedStatus` and
#' `PTMEAExpectedReason` identify the calculation and unavailable states.
#' At `N = 2` the variance term vanishes; this does not establish perfect
#' measurement or fit. Below `N = 2`, the corrected reference is unavailable
#' (FACETS 4.5.0 can instead display a boundary value). The asymptotic reference
#' remains available if its moments exist. Arbitrarily normalizing weights
#' changes the corrected reference: use replication semantics for comparison.
#' Constant predictors, invalid moments, or observations involving unanchored
#' extreme main effects withhold the complete element's expected reference;
#' they are not silently dropped. Fixed finite anchors remain usable.
#'
#' FACETS 4.5.0 fixed-calibration PCM checks cover missing assignments and
#' unit, integer, fractional and unequal replication weights. These checks
#' establish a reporting convention within that scope, not general diagnostic
#' calibration. Observed `PTMEA` retains its existing definition. MML, corrected
#' JML, GPCM, fitted interactions and group-mean anchors do not receive the new
#' expected-correlation columns. Recreate older diagnostics from the saved fit
#' to obtain them; old saved results are not relabeled or silently recalculated.
#' See Linacre (2008), *The Expected Value of a Point-Biserial (or Similar)
#' Correlation*, Rasch Measurement Transactions, 22:1, 1154, including its later
#' finite-count note: \url{https://www.rasch.org/rmt/rmt221e.htm}.
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("diagnose_mfrm", run.dontrun = TRUE)
#' \dontrun{
#' # Each row is one bundled example rating; use ordinary JML for this reference
#' ratings <- load_mfrmr_data("example_core")
#' jml_fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'                      method = "JML", model = "RSM")
#' jml_diagnostics <- diagnose_mfrm(jml_fit, residual_pca = "none",
#'                                  diagnostic_mode = "legacy")
#' head(jml_diagnostics$measures[, c("Facet", "Level", "PTMEA",
#'      "PTMEAExpected", "PTMEAExpectedN", "PTMEAExpectedStatus")])
#' # Compare observed PTMEA with its conditional reference; neither is a test
#' }
#' @name diagnose_mfrm
#' @rdname diagnose_mfrm
NULL

# Moment ratio with a declared frequency count, evaluated without raw-square
# cancellation. Missing rows invalidate the complete scope, not just a subset.
.mfrm_expected_correlation_moments <- function(x, e, v, w, eligible = rep(TRUE, length(x))) {
  n <- sum(w)
  out <- data.frame(PTMEAExpected = NA_real_, PTMEAExpectedAsymptotic = NA_real_,
    PTMEAExpectedN = n, PTMEAExpectedFactor = if (is.finite(n) && n >= 2) (n-2)/n else NA_real_,
    PTMEAExpectedStatus = "unavailable", PTMEAExpectedReason = "Incomplete conditional moments or invalid weights.")
  if (!length(x) || any(!is.finite(w) | w < 0) || !is.finite(n) || n <= 0) return(out)
  keep <- w > 0
  x <- x[keep]; e <- e[keep]; v <- v[keep]; w <- w[keep]; eligible <- eligible[keep]
  if (any(!is.finite(x) | !is.finite(e) | !is.finite(v) | v < 0)) return(out)
  if (anyNA(eligible) || !all(eligible)) {
    out$PTMEAExpectedReason <- "The scope includes unanchored extreme main effects or nonfinite primary Person estimates."
    return(out)
  }
  w <- w/n
  x <- x-x[1L]; e <- e-e[1L]
  x <- x-sum(w*x); e <- e-sum(w*e)
  vx <- sum(w*x*x); ve <- sum(w*e*e); cv <- sum(w*x*e); vm <- sum(w*v)
  if (!all(is.finite(c(vx,ve,cv,vm))) || vx <= 0) {
    out$PTMEAExpectedReason <- "The oriented predictor is constant or its moments are nonfinite."
    return(out)
  }
  ratio <- function(factor) {
    denom <- sqrt(vx)*sqrt(ve+factor*vm)
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    rho <- cv/denom
    if (!is.finite(rho) || abs(rho) > 1+1e-10) return(NA_real_)
    max(-1,min(1,rho))
  }
  out$PTMEAExpectedAsymptotic <- ratio(1)
  if (n < 2) {
    out$PTMEAExpectedStatus <- "insufficient_frequency"
    out$PTMEAExpectedReason <- "Frequency count below two; the finite-count correction is unavailable."
    return(out)
  }
  out$PTMEAExpected <- ratio(out$PTMEAExpectedFactor)
  if (is.finite(out$PTMEAExpected)) {
    out$PTMEAExpectedStatus <- "available_conditional"
    out$PTMEAExpectedReason <- "Conditional finite-count reporting approximation; not an exact sample-correlation expectation or inferential test."
  } else {
    out$PTMEAExpectedReason <- "The corrected denominator is zero or nonfinite."
  }
  out
}

.mfrm_expected_ptmea <- function(res, obs, facet_cols) {
  eligible <- rep(TRUE,nrow(obs))
  extreme <- .mfrm_unanchored_extremes(res,obs)
  for (f in names(extreme)) eligible <- eligible & !as.character(obs[[f]]) %in% extreme[[f]]
  primary <- res$facets$person$Estimate[match(as.character(obs$Person), as.character(res$facets$person$Person))]
  eligible <- eligible & is.finite(primary)
  w <- get_weights(obs)
  variance <- obs[["ConditionalVariance"]] %||% obs$Var
  purrr::map_dfr(facet_cols,function(f) {
    parts <- split(seq_len(nrow(obs)),as.character(obs[[f]]))
    purrr::map_dfr(names(parts),function(level) {
      rows <- parts[[level]]
      dplyr::bind_cols(tibble(Facet=f,Level=level),
        .mfrm_expected_correlation_moments(obs$MeasureSum[rows], obs$Expected[rows],
          variance[rows], w[rows], eligible[rows]))
    })
  })
}
