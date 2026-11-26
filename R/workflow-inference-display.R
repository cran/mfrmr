.tidymfrm_intervals <- function(fit, component, level, real = TRUE) {
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("`conf.level` must be between 0 and 1.", call. = FALSE)
  n <- NROW(fit[[component]])
  if (!summary(fit)$readiness$ready) return(matrix(NA_real_, n, 2L))
  stats::confint(fit, component = component, level = level,
                 se_type = if (real) "real" else "model")
}

.tidymfrm_cut_decisions <- function(fit, cuts, readiness = summary(fit)$readiness) {
  if (!is.numeric(cuts) || any(!is.finite(cuts)))
    stop("`cut_scores` must contain finite numeric scores.", call. = FALSE)
  pt <- fit$person
  enc <- fit$encoding
  sc <- fit$step_cumulative
  # Fair(M) neutralizes facet locations, uses the geometric reference slope,
  # and averages expected scores over each person's weighted step-group mix.
  base <- sum(vapply(enc$facets, function(fc) {
    tb <- if (fc %in% enc$base_facets) fit$facets[fit$facets$Facet == fc, ] else
      fit$interactions[fit$interactions$Interaction == fc, ]
    values <- tb$Estimate[is.finite(tb$Estimate)]
    enc$facet_signs[match(fc, enc$facets)] * if (length(values)) mean(values) else 0
  }, numeric(1)))
  valid_alpha <- is.finite(fit$alpha_obs) & fit$alpha_obs > 0 & enc$weight > 0
  alpha <- if (any(valid_alpha)) exp(stats::weighted.mean(log(fit$alpha_obs[valid_alpha]), enc$weight[valid_alpha])) else 1
  mixes <- lapply(seq_len(nrow(pt)), function(i) {
    rows <- which(enc$person0 == i - 1L & enc$weight > 0)
    w <- numeric(nrow(sc))
    for (g in seq_len(nrow(sc))) w[g] <- sum(enc$weight[rows][enc$step_group0[rows] == g - 1L])
    w / sum(w)
  })
  expected <- function(theta, mix) {
    active <- which(mix > 0)
    sum(vapply(active, function(g) {
      k <- seq.int(enc$step_cat_min[g], enc$step_cat_max[g])
      cum <- sc[g, k + 1L]
      lin <- if (enc$slope_formula == 1L) alpha * theta + base else alpha * (theta + base)
      lp <- (k - min(k)) * lin - if (enc$slope_formula == 1L) cum else alpha * cum
      p <- exp(lp - max(lp)); p <- p / sum(p)
      mix[g] * (fit$rating_min + sum(k * p))
    }, numeric(1)))
  }
  tibble::as_tibble(do.call(rbind, lapply(sort(unique(cuts)), function(ct) {
    threshold <- vapply(mixes, function(mix) {
      if (any(!is.finite(mix))) return(NA_real_)
      active <- which(mix > 0)
      bounds <- fit$rating_min + c(sum(mix[active] * enc$step_cat_min[active]),
                                   sum(mix[active] * enc$step_cat_max[active]))
      if (ct <= bounds[1L] || ct >= bounds[2L]) return(NA_real_)
      width <- 8
      while (width < 1e6 && (expected(-width, mix) > ct || expected(width, mix) < ct)) width <- width * 2
      tryCatch(stats::uniroot(function(th) expected(th, mix) - ct,
                              c(-width, width), tol = 1e-9)$root,
               error = function(e) NA_real_)
    }, numeric(1))
    probability <- function(se) {
      ok <- readiness$ready & is.finite(threshold) & is.finite(pt$Estimate) & is.finite(se) & se > 0
      p <- rep(NA_real_, nrow(pt))
      p[ok] <- stats::pnorm((pt$Estimate[ok] - threshold[ok]) / se[ok])
      p
    }
    pm <- probability(pt$ModelSE); pr <- probability(pt$RealSE)
    count <- function(p) if (any(is.finite(p))) sum(p > .05 & p < .95, na.rm = TRUE) else NA_integer_
    finite_cut <- threshold[is.finite(threshold)]
    pass <- pt$FairM >= ct
    why <- unique(c(readiness$reason_codes,
      if (any(!is.finite(threshold))) "cut_outside_open_fair_score_range",
      if (readiness$ready && any(!is.finite(pm) | !is.finite(pr))) "some_person_intervals_unavailable"))
    data.frame(CutFairM = ct,
      CutLogit = if (length(finite_cut) > 0L && diff(range(finite_cut)) < 1e-8) finite_cut[1L] else NA_real_,
      CutLogitMin = if (length(finite_cut)) min(finite_cut) else NA_real_,
      CutLogitMax = if (length(finite_cut)) max(finite_cut) else NA_real_,
      PassCount = if (any(!is.na(pass))) sum(pass, na.rm = TRUE) else NA_integer_,
      PassAvailable = sum(!is.na(pass)),
      PassRatePct = if (any(!is.na(pass))) 100 * mean(pass, na.rm = TRUE) else NA_real_,
      Borderline90_ModelSE = count(pm), Borderline90_RealSE = count(pr),
      Available_ModelSE = sum(is.finite(pm)), Available_RealSE = sum(is.finite(pr)),
      Unavailable_ModelSE = sum(!is.finite(pm)), Unavailable_RealSE = sum(!is.finite(pr)),
      InferenceStatus = if (!length(why)) "ready" else "review_required",
      InferenceReason = paste(why, collapse = "; "), stringsAsFactors = FALSE)
  })))
}

.tidymfrm_model_label <- function(fit, fit_summary) {
  fit_summary$model_label %||% .mfrm_model_label(fit$model, unique(fit$slopes$Facet))
}

.tidymfrm_linkage_notes <- function(wf) {
  x <- wf$connectivity$rater_linkage
  if (is.null(x)) return(character())
  result <- sprintf("Person-by-%s graph: %d connected component(s).", x$facet, x$n_components)
  if (x$n_components > 1L) {
    result <- c(result, paste("No persons are shared across these components.",
                              "Cross-component comparisons depend on calibration assumptions or anchors."))
    if (!is.null(wf$fit) && identical(wf$fit$method, "MML") &&
        identical(wf$fit$population$mode, "fixed")) {
      result <- c(result, "This MML calibration uses a common fixed ability distribution.")
    }
  }
  if (length(x$bridge_persons)) {
    result <- c(result, paste(
      "Removing any one of these linking persons separates the observed rater panels:",
      paste(x$bridge_persons, collapse = ", ")
    ))
  }
  result
}

.tidymfrm_observed_scores <- function(wf, raters) {
  enc <- wf$fit$encoding
  score <- enc$score0 + wf$fit$rating_min
  categories <- seq.int(wf$fit$rating_min, wf$fit$rating_max)
  rater <- enc$retained_data[[wf$spec$rater_facet]]
  counts <- table(factor(rater, levels = raters$Level), factor(score, levels = categories))
  totals <- rowSums(counts)
  raters$AtMinimum <- as.integer(counts[, 1L])
  raters$AtMaximum <- as.integer(counts[, ncol(counts)])
  range_end <- function(i, end) {
    used <- categories[counts[i, ] > 0]
    if (length(used)) end(used) else NA_real_
  }
  raters$ObservedMin <- vapply(seq_len(nrow(counts)), range_end, numeric(1), end = min)
  raters$ObservedMax <- vapply(seq_len(nrow(counts)), range_end, numeric(1), end = max)
  raters$ScorePattern <- ifelse(totals > 0 & raters$AtMinimum == totals, "All ratings at minimum",
                                ifelse(totals > 0 & raters$AtMaximum == totals, "All ratings at maximum", ""))
  long <- as.data.frame(counts, responseName = "Count", stringsAsFactors = FALSE)
  names(long)[1:2] <- c("Level", "Category")
  long$Category <- as.numeric(long$Category)
  list(raters = raters, counts = long)
}

.tidymfrm_diagnostic_fit_basis <- function(fit) {
  fit[c("model", "method", "score_method", "loglik", "xi")]
}

.tidymfrm_check_diagnostics <- function(wf) {
  diagnostics <- wf$diagnostics
  fit <- wf$fit
  if (is.null(diagnostics) || is.null(fit)) return(wf)

  matches <- if (!is.null(diagnostics$fit_basis)) {
    identical(diagnostics$fit_basis, .tidymfrm_diagnostic_fit_basis(fit))
  } else {
    # Older saved workflows have no fit identity but retain fitted residuals.
    observations <- diagnostics$observations
    isTRUE(all.equal(observations$Expected, as.numeric(stats::fitted(fit, type = "response")),
                     tolerance = 1e-10)) &&
      isTRUE(all.equal(observations$StdResidual, as.numeric(stats::residuals(fit, type = "pearson")),
                       tolerance = 1e-10))
  }
  if (!matches) wf[c("diagnostics", "method_comparison")] <- list(NULL)
  wf
}

# Recheck saved workflows at every reporting entry without refitting the model.
.tidymfrm_refresh <- function(wf) {
  wf <- .tidymfrm_check_diagnostics(wf)
  old_linkage <- wf$connectivity$rater_linkage
  if (!"rater_linkage" %in% names(wf$connectivity) ||
      (!is.null(old_linkage) && !"robustness" %in% names(old_linkage))) {
    sp <- wf$spec
    conn <- check_mfrm_connectivity(wf$data, sp$person, sp$facets, sp$score,
      rating_min = sp$rating_min, rating_max = sp$rating_max, rater_facet = sp$rater_facet, network = TRUE)
    wf$connectivity["rater_linkage"] <- list(conn$rater_linkage)
  }
  if (!is.null(wf$fit$encoding) && !isTRUE(wf$connectivity_calibrated)) {
    wf$connectivity <- c(.mfrm_connectivity_summary(wf$fit$encoding), wf$connectivity["rater_linkage"])
    wf$connectivity_calibrated <- TRUE
  }
  bootstrap <- wf[["bootstrap"]]
  if (!is.null(bootstrap) && !("BootstrapAvailable" %in% names(bootstrap))) {
    wf["bootstrap"] <- list(NULL)
    wf$bootstrap_reason <- paste(
      "The saved bootstrap predates model-specification and refit checks.",
      "Refit with the current workflow for bootstrap intervals."
    )
  }
  if (!is.null(wf$gtheory)) {
    gt <- .mfrm_refresh_gstudy(wf$gtheory)
    gt$d_study$MeetsTargetG <- gt$d_study$G >= gt$targets[["G"]]
    gt$d_study$MeetsTargetPhi <- gt$d_study$Phi >= gt$targets[["Phi"]]
    for (key in c("g", "phi")) {
      col <- if (key == "g") "G" else "Phi"
      ok <- gt$d_study[[paste0("MeetsTarget", col)]] %in% TRUE
      rows <- gt$d_study[ok, , drop = FALSE]
      gt[[paste0("recommended_for_", key)]] <- rows[utils::head(order(rows$TotalRatingsPerPerson, -rows[[col]]), 1L), , drop = FALSE]
    }
    wf$gtheory <- gt
  }
  if (is.null(wf$fit) || is.null(wf$diagnostics)) return(wf)
  r <- summary(wf$fit)$readiness
  wf$diagnostics$readiness <- r
  rt <- wf$diagnostics$raters
  if (is.null(wf$diagnostics$category_counts)) {
    scores <- .tidymfrm_observed_scores(wf, rt)
    rt <- scores$raters
    wf$diagnostics$category_counts <- scores$counts
  }
  rt$InferenceStatus <- r$status
  rt$InferenceReason <- paste(r$reason_codes, collapse = "; ")
  if (!r$ready) {
    rt$DiagnosticFlag <- "Review required"
    rt$RecommendedAction <- paste("Review calibration:", rt$InferenceReason)
    if (!is.null(wf$diagnostics$dff)) {
      tb <- wf$diagnostics$dff
      for (nm in intersect(c("pValue", "pAdj", "Flagged"), names(tb))) tb[[nm]] <- NA
      tb$InferenceStatus <- r$status; tb$InferenceReason <- rt$InferenceReason[1L]
      wf$diagnostics$dff <- tb
    }
  } else {
    endpoint <- nzchar(rt$ScorePattern)
    rt$RecommendedAction[endpoint] <- "Check the raw score records and rubric use for the observed endpoint concentration."
  }
  wf$diagnostics$raters <- rt
  if (!is.null(wf$diagnostics$decisions)) wf$diagnostics$decisions <-
    .tidymfrm_cut_decisions(wf$fit, wf$diagnostics$decisions$CutFairM, r)
  wf
}

.tidymfrm_flag_count <- function(rt) {
  ok <- !is.na(rt$DiagnosticFlag) & rt$DiagnosticFlag != "Review required"
  if (any(ok)) sum(rt$DiagnosticFlag[ok] != "Acceptable") else NA_integer_
}
