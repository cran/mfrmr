# ==============================================================================
# Tidy Extractors: tidy(), glance(), augment()
# ==============================================================================

































generics::tidy





















generics::glance




















generics::augment



tidy.tidymfrm <- function(x,
                          component = c("raters", "facets", "persons", "steps", "slopes",
                                        "g_study", "d_study", "dff", "drift"),
                          conf.level = 0.95,
                          use_real_se = TRUE,
                          ...) {
  component <- match.arg(component)
  x <- .tidymfrm_refresh(x)

  if (component == "raters") {
    if (is.null(x$diagnostics)) stop("No rater diagnostics found; run `step_mfrm_diagnose()`.", call. = FALSE)
    rt <- x$diagnostics$raters
    ci <- .tidymfrm_intervals(x$fit, "facets", conf.level, use_real_se)
    ci <- ci[match(paste(rt$Facet, rt$Level), paste(x$fit$facets$Facet, x$fit$facets$Level)), , drop = FALSE]
    return(tibble::tibble(
      facet = rt$Facet,
      term = rt$Level,
      estimate = rt$Estimate,
      std.error = rt$ModelSE,
      real.se = rt$RealSE,
      conf.low = ci[, 1L], conf.high = ci[, 2L],
      inference_status = rt$InferenceStatus, inference_reason = rt$InferenceReason,
      fair_m = rt$FairM,
      infit_mnsq = rt$InfitMnSq,
      infit_zstd = rt$InfitZSTD,
      outfit_mnsq = rt$OutfitMnSq,
      outfit_zstd = rt$OutfitZSTD,
      observed_min = rt$ObservedMin,
      observed_max = rt$ObservedMax,
      at_minimum = rt$AtMinimum,
      at_maximum = rt$AtMaximum,
      score_pattern = rt$ScorePattern,
      flag = rt$DiagnosticFlag,
      action = rt$RecommendedAction
    ))
  }
  if (component %in% c("facets", "persons", "steps", "slopes")) {
    return(tidy.mfrm_fit(x$fit, component = component, conf.level = conf.level, use_real_se = use_real_se))
  }
  if (component == "g_study") {
    return(if (is.null(x$gtheory)) tibble::tibble() else tidy(x$gtheory, conf.level = conf.level))
  }
  if (component == "d_study") {
    return(x$gtheory$d_study)
  }
  if (component == "dff") {
    if (is.null(x$diagnostics$dff)) return(tibble::tibble())
    return(x$diagnostics$dff)
  }
  if (component == "drift") {
    if (is.null(x$drift)) return(tibble::tibble())
    return(tibble::as_tibble(as.data.frame(x$drift, component = "severity")))
  }
}



tidy.mfrm_fit <- function(x, component = c("facets", "persons", "steps", "slopes"), conf.level = 0.95, use_real_se = TRUE,
    ...) {
    .mfrm_validate_fit(x)
    component <- match.arg(component)
    ci <- .tidymfrm_intervals(x, if (component == "persons")
        "person"
    else component, conf.level, use_real_se)
    readiness <- summary(x)$readiness
    status <- readiness$status
    reason <- paste(readiness$reason_codes, collapse = "; ")
    if (component == "slopes") {
        sl <- x$slopes
        if (is.null(sl) || nrow(sl) == 0L)
            return(tibble::tibble())
        return(tibble::tibble(facet = sl$Facet, term = sl$Level, estimate = sl$Estimate, log.estimate = sl$LogEstimate,
            std.error = sl$ModelSE, sandwich.se = sl$SandwichSE, conf.low = ci[, 1L], conf.high = ci[, 2L],
            inference_status = status, inference_reason = reason, formula = sl$Formula))
    }
    if (component == "facets") {
        fc <- x$facets
        return(tibble::tibble(facet = fc$Facet, term = fc$Level, estimate = fc$Estimate, std.error = fc$ModelSE,
            real.se = fc$RealSE, sandwich.se = fc$SandwichSE, conf.low = ci[, 1L], conf.high = ci[, 2L], inference_status = status,
            inference_reason = reason, fair_m = fc$FairM, infit_mnsq = fc$InfitMnSq, infit_zstd = fc$InfitZSTD,
            outfit_mnsq = fc$OutfitMnSq, outfit_zstd = fc$OutfitZSTD))
    }
    if (component == "persons") {
        pt <- x$person
        return(tibble::tibble(term = pt$Person, estimate = pt$Estimate, std.error = pt$ModelSE, real.se = pt$RealSE,
            wle = pt$WLE, wle.se = pt$WLE_SE, conf.low = ci[, 1L], conf.high = ci[, 2L], inference_status = status,
            inference_reason = reason, fair_m = pt$FairM, infit_mnsq = pt$InfitMnSq, infit_zstd = pt$InfitZSTD,
            extreme_flag = pt$ExtremeFlag))
    }
    st <- x$steps
    tibble::tibble(scale = st$ScaleLevel, step = st$Step, estimate = st$Estimate, std.error = st$ModelSE,
        conf.low = ci[, 1L], conf.high = ci[, 2L], inference_status = status, inference_reason = reason)
}



tidy.mfrm_gtheory <- function(x, component = c("g_study", "d_study"), conf.level = 0.95, ...) {
  .mfrm_validate_gstudy(x)
  if (x$method_code != "urgenova") {
    component <- match.arg(component)
    if (component == "d_study") return(tibble::as_tibble(mfrm_d_study(x)$coefficients))
    rows <- lapply(names(x$components), function(source) {
      m <- x$components[[source]]
      ij <- which(upper.tri(m, diag = TRUE), arr.ind = TRUE)
      data.frame(Component = source, Score1 = rownames(m)[ij[,1L]],
        Score2 = colnames(m)[ij[,2L]], Estimate = m[ij],
        interval_status = "not_requested", stringsAsFactors = FALSE)
    })
    return(tibble::as_tibble(do.call(rbind, rows)))
  }
  component <- match.arg(component)
  if (component == "d_study") return(tibble::as_tibble(mfrm_d_study(x)))
  ci <- stats::confint(x, level = conf.level, component = "variance_components")
  out <- tibble::as_tibble(x$variance_components)
  out$conf.low <- ci[, 1L]
  out$conf.high <- ci[, 2L]
  out$interval_status <- unname(attr(ci, "status"))
  out
}



glance.tidymfrm <- function(x, ...) {
  x <- .tidymfrm_refresh(x)
  fit_summary <- if (!is.null(x$fit)) summary(x$fit)
  readiness <- fit_summary$readiness
  tibble::tibble(
    model = if (!is.null(x$fit)) x$fit$model else NA_character_,
    method = if (!is.null(x$fit)) x$fit$method else NA_character_,
    score_method = if (!is.null(x$fit)) x$fit$score_method else NA_character_,
    n_obs = x$connectivity$n_obs,
    n_person = x$connectivity$n_person,
    connected = x$connectivity$connected,
    n_subsets = x$connectivity$n_subsets,
    converged = if (!is.null(x$fit)) x$fit$converged else NA,
    inference_ready = if (!is.null(readiness)) readiness$ready else NA,
    inference_status = if (!is.null(readiness)) readiness$status else "not_assessed",
    logLik = if (!is.null(x$fit)) x$fit$loglik else NA_real_,
    eap_reliability = if (!is.null(x$fit)) x$fit$population$eap_reliability else NA_real_,
    g_coef = if (!is.null(x$gtheory)) x$gtheory$g_study_summary$G else NA_real_,
    phi_coef = if (!is.null(x$gtheory)) x$gtheory$g_study_summary$Phi else NA_real_,
    flagged_raters = if (!is.null(x$diagnostics)) .tidymfrm_flag_count(x$diagnostics$raters) else NA_integer_,
    unavailable_raters = if (!is.null(x$diagnostics)) sum(x$diagnostics$raters$DiagnosticFlag == "Review required") else NA_integer_,
    unexpected_ratings = if (!is.null(x$diagnostics)) nrow(x$diagnostics$unexpected) else NA_integer_,
    elapsed_sec = x$elapsed_sec,
    model_label = if (!is.null(fit_summary)) .tidymfrm_model_label(x$fit, fit_summary) else NA_character_
  )
}



glance.mfrm_fit <- function(x, ...) {
    .mfrm_validate_fit(x)
    fit_summary <- summary(x)
    readiness <- fit_summary$readiness
    tibble::tibble(model = x$model, method = x$method, score_method = x$score_method, n_obs = x$n_obs, n_person = x$n_person,
        connected = x$connectivity$connected, n_subsets = x$connectivity$n_subsets, converged = x$converged,
        inference_ready = readiness$ready, inference_status = readiness$status, logLik = x$loglik, eap_reliability = x$population$eap_reliability,
        elapsed_sec = x$elapsed_sec, model_label = .tidymfrm_model_label(x, fit_summary))
}



augment.tidymfrm <- function(x, ...) {
  x <- .tidymfrm_check_diagnostics(x)
  if (!is.null(x$diagnostics$observations)) {
    return(x$diagnostics$observations)
  }
  if (is.null(x$fit)) stop("Run `step_mfrm_measure()` before `augment()`.", call. = FALSE)
  .tidymfrm_observation_residuals(x$data, x$spec$person, x$spec$facets, x$spec$score, x$fit)
}
