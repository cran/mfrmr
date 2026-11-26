# ==============================================================================
# Unified End-to-End & Pipeable Workflow: Connectivity -> MFRM -> G/D -> Audit
# Calls ONLY exported public functions of `mfrmr` (zero triple-colon calls)
# ==============================================================================

.tidymfrm_infer_rater_facet <- function(facets) {
  hit <- grep("rater|judge|examiner|marker|grader|reviewer", facets, ignore.case = TRUE, value = TRUE)
  if (length(hit) >= 1L) hit[1L] else facets[1L]
}

.tidymfrm_observation_residuals <- function(data, person, facets, score, fit) {
  raw_y <- as.numeric(data[[score]])
  keep <- is.finite(raw_y) & !is.na(data[[person]])
  for (fc in facets) keep <- keep & !is.na(data[[fc]])
  d <- data[keep, , drop = FALSE]

  exp_y <- as.numeric(stats::fitted(fit, type = "response"))
  std_res <- as.numeric(stats::residuals(fit, type = "pearson"))
  probs <- stats::fitted(fit, type = "probs")
  y_idx <- as.integer(round(as.numeric(d[[score]]) - fit$rating_min)) + 1L
  prob_obs <- probs[cbind(seq_len(nrow(d)), y_idx)]
  obs_orig <- as.numeric(d[[score]])
  resid <- obs_orig - exp_y
  mod_sd <- ifelse(is.finite(std_res) & abs(std_res) > 1e-12, abs(resid / std_res), NA_real_)

  out <- tibble::as_tibble(d[, c(person, facets, score), drop = FALSE])
  out$.row_id <- which(keep)
  out$Observed <- obs_orig
  out$Expected <- exp_y
  out$Residual <- resid
  out$ModelSD <- mod_sd
  out$StdResidual <- std_res
  out$ObsProb <- prob_obs
  out$UnexpectedFlag <- is.finite(std_res) & (abs(std_res) >= 2.0)
  out
}


























mfrm_workflow <- function(data,
                          person,
                          facets,
                          score,
                          rater_facet = NULL,
                          occasion = NULL,
                          rating_min = NULL,
                          rating_max = NULL,
                          fit = NULL,
                          expected_design = NULL,
                          rater_panels = NULL) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("`data` must be a non-empty data.frame.", call. = FALSE)
  }
  if (is.null(rater_facet)) rater_facet <- .tidymfrm_infer_rater_facet(facets)
  if (!is.character(rater_facet) || length(rater_facet) != 1L || is.na(rater_facet) ||
      !rater_facet %in% facets) {
    stop("`rater_facet` must name one of the columns in `facets`.", call. = FALSE)
  }
  anchor_drift <- NULL
  if (inherits(fit, "mfrm_anchored_fit")) {
    anchor_drift <- fit$drift
    fit <- fit$fit
  }
  if (!is.null(fit)) {
    .mfrm_validate_fit(fit)
    enc <- fit$encoding
    if (!identical(c(person, facets, score), c(enc$person_col, enc$base_facets, enc$score_col))) {
      stop("The fitted calibration must use the same person, facet, and score columns.", call. = FALSE)
    }
    if ((!is.null(rating_min) && rating_min != fit$rating_min) ||
        (!is.null(rating_max) && rating_max != fit$rating_max)) {
      stop("Rating bounds must match the fitted calibration.", call. = FALSE)
    }
    columns <- c(person, facets, score)
    complete <- is.finite(as.numeric(data[[score]])) & stats::complete.cases(data[columns])
    same <- vapply(columns, function(column) {
      identical(as.character(data[[column]][complete]), as.character(enc$retained_data[[column]]))
    }, logical(1))
    if (!all(same)) stop("Supply the same retained ratings in their original order.", call. = FALSE)
    rating_min <- fit$rating_min
    rating_max <- fit$rating_max
  }
  conn <- if (is.null(fit)) {
    check_mfrm_connectivity(
      data = data, person = person, facets = facets, score = score,
      rating_min = rating_min, rating_max = rating_max, rater_facet = rater_facet,
      network = TRUE, rater_panels = rater_panels
    )
  } else {
    c(.mfrm_connectivity_summary(enc), list(
      rater_linkage = .mfrm_rater_linkage(enc, rater_facet, network = TRUE, rater_panels = rater_panels)
    ))
  }
  roster <- if (!is.null(expected_design)) {
    summary(describe_mfrm_data(data, person, facets, score,
      rating_min = rating_min, rating_max = rating_max, include_agreement = FALSE,
      expected_design = expected_design))$structural_missingness
  }
  structure(
    list(
      data = data,
      spec = list(
        person = person,
        facets = facets,
        score = score,
        rater_facet = rater_facet,
        occasion = occasion,
        rating_min = rating_min,
        rating_max = rating_max
      ),
      connectivity = conn,
      connectivity_calibrated = !is.null(fit),
      planned_coverage = roster,
      fit = fit,
      anchor_drift = anchor_drift,
      drift = NULL,
      method_comparison = NULL,
      bootstrap = NULL,
      gtheory = NULL,
      diagnostics = NULL,
      elapsed_sec = 0.0
    ),
    class = "tidymfrm"
  )
}

















































step_mfrm_measure <- function(wf,
                              model = c("RSM", "PCM", "GPCM"),
                              method = c("MML", "JML", "JML_jeffreys", "CML"),
                              step_facet = NULL,
                              slope_facet = NULL,
                              slope_reference = NULL,
                              slope_type = c("full_predictor", "ability_only"),
                              slope_prior_sd = NULL,
                              occasion = NULL,
                              person_score = c("auto", "EAP", "WLE"),
                              compare_methods = NULL,
                              population_formula = NULL,
                              person_data = NULL,
                              step_constraints = NULL,
                              category_support = NULL,
                              facet_formulas = NULL,
                              facet_data = NULL,
                              quad_points = 31L,
                              mml_integration = c("fixed", "adaptive"),
                              xtreme = 0.3,
                              warn_connectivity = TRUE,
                              bootstrap = 0L,
                              threads = 1L,
                              gpcm_mml_identification = NULL) {
  if (!inherits(wf, "tidymfrm")) stop("`wf` must be a `tidymfrm` object.", call. = FALSE)
  t0 <- proc.time()[[3L]]
  model <- .mfrm_model_code(model)
  method <- match.arg(method[1L], c("MML", "JML", "JML_jeffreys", "CML"))
  slope_type <- match.arg(slope_type)
  person_score <- match.arg(person_score)
  mml_integration <- match.arg(mml_integration)
  sp <- wf$spec
  occ_col <- occasion %||% sp$occasion

  model_spec <- list(slope_reference = slope_reference, slope_type = slope_type)
  if (!is.null(gpcm_mml_identification)) model_spec$gpcm_mml_identification <- gpcm_mml_identification
  fit_with_method <- function(estimator, warn) fit_mfrm(
    data = wf$data,
    person = sp$person,
    facets = sp$facets,
    score = sp$score,
    model = model,
    method = estimator,
    step_facet = step_facet,
    slope_facet = slope_facet,
    person_score = person_score,
    rating_min = sp$rating_min,
    rating_max = sp$rating_max,
    population_formula = population_formula,
    person_data = person_data,
    step_constraints = step_constraints,
    category_support = category_support,
    facet_formulas = facet_formulas,
    facet_data = facet_data,
    warn_connectivity = warn,
    spec = model_spec,
    control = list(
      quad_points = quad_points,
      mml_integration = mml_integration,
      xtreme = if (identical(estimator, "JML")) xtreme else 0.0,
      slope_prior_sd = slope_prior_sd
    )
  )
  fit <- fit_with_method(method, warn_connectivity)
  wf$fit <- fit
  wf$connectivity <- c(.mfrm_connectivity_summary(fit$encoding), wf$connectivity["rater_linkage"])
  wf$connectivity_calibrated <- TRUE
  wf[c("anchor_drift", "drift", "method_comparison", "diagnostics")] <- list(NULL)
  wf$spec$model <- model
  wf$spec$method <- method
  wf$spec$step_facet <- step_facet
  wf$spec$slope_facet <- slope_facet
  wf$spec$occasion <- occ_col
  wf$spec$quad_points <- quad_points

  if (!is.null(occ_col) && nzchar(occ_col)) {
    if (!occ_col %in% names(wf$data)) {
      stop("`occasion` column '", occ_col, "' not found in `data`.", call. = FALSE)
    }
    other_f <- setdiff(sp$facets, c(sp$rater_facet, occ_col))
    wf$drift <- fit_mfrm_drift(
      data = wf$data,
      person = sp$person,
      rater = sp$rater_facet,
      time = occ_col,
      score = sp$score,
      facets = other_f,
      model = if (identical(model, "RSM")) "RSM" else "PCM",
      quad_points = 81L,
      missing = "omit"
    )
  }

  if (!is.null(compare_methods) && length(compare_methods) > 0L) {
    all_m <- unique(c(method, compare_methods))
    cmp_rows <- list()
    for (m in all_m) {
      fm <- if (identical(m, method)) {
        fit
      } else {
        fit_with_method(m, FALSE)
      }
      ok_p <- is.finite(fit$person$Estimate) & is.finite(fm$person$Estimate)
      ok_f <- is.finite(fit$person$FairM) & is.finite(fm$person$FairM)
      rf0 <- fit$facets[fit$facets$Facet == sp$rater_facet, , drop = FALSE]
      rfm <- fm$facets[fm$facets$Facet == sp$rater_facet, , drop = FALSE]
      cmp_rows[[length(cmp_rows) + 1L]] <- data.frame(
        Method = m,
        ScoreMethod = fm$score_method,
        Converged = fm$converged,
        InferenceStatus = summary(fm)$readiness$status,
        PersonCor = if (sum(ok_p) >= 2L) stats::cor(fit$person$Estimate[ok_p], fm$person$Estimate[ok_p]) else NA_real_,
        PersonMAD = if (sum(ok_p) >= 1L) stats::median(abs(fit$person$Estimate[ok_p] - fm$person$Estimate[ok_p])) else NA_real_,
        FairM_MaxDiff = if (sum(ok_f) >= 1L) max(abs(fit$person$FairM[ok_f] - fm$person$FairM[ok_f])) else NA_real_,
        RaterMaxDiff = max(abs(rf0$Estimate - rfm$Estimate), na.rm = TRUE),
        ElapsedSec = fm$elapsed_sec,
        stringsAsFactors = FALSE
      )
    }
    wf$method_comparison <- tibble::as_tibble(do.call(rbind, cmp_rows))
  }

  wf["bootstrap"] <- list(NULL)
  wf$bootstrap_reason <- NULL
  if (is.numeric(bootstrap) && bootstrap >= 2L) {
    B <- as.integer(bootstrap)[1L]
    readiness <- summary(fit)$readiness
    if (!isTRUE(readiness$ready)) {
      wf$bootstrap_reason <- paste("The source calibration requires review:",
                                   paste(readiness$reason_codes, collapse = ", "))
    } else {
      probs <- stats::fitted(fit, type = "probs")
      if (any(!is.finite(probs)) || any(probs < 0) || any(rowSums(probs) <= 0)) {
        wf$bootstrap_reason <- "Some retained ratings have unavailable fitted category probabilities."
      } else {
        cats_int <- as.integer(colnames(probs))
        smat <- t(vapply(seq_len(nrow(probs)), function(i) {
          sample(cats_int, size = B, replace = TRUE, prob = probs[i, ])
        }, integer(B)))
        bfits <- mfrm_batch_fit(
          fit$encoding$retained_data, person = sp$person, facets = sp$facets, score = sp$score,
          score_matrix = smat, model = model, method = method, step_facet = step_facet,
          slope_facet = slope_facet, slope_reference = slope_reference, slope_type = slope_type,
          gpcm_mml_identification = fit$encoding$gpcm_mml_identification,
          rating_min = fit$rating_min, rating_max = fit$rating_max, person_score = person_score,
          population_formula = population_formula, person_data = person_data,
          step_constraints = step_constraints, category_support = category_support,
          facet_formulas = facet_formulas, facet_data = facet_data, slope_prior_sd = slope_prior_sd,
          quad_points = quad_points, mml_integration = mml_integration,
          xtreme = if (identical(method, "JML")) xtreme else 0.0, threads = threads
        )
        n_levels <- nrow(fit$facets)
        estimates <- matrix(vapply(bfits, function(bf) {
          if (isTRUE(summary(bf)$readiness$ready)) bf$facets$Estimate else rep(NA_real_, n_levels)
        }, numeric(n_levels)), nrow = n_levels)
        available <- rowSums(is.finite(estimates))
        intervals <- vapply(seq_len(n_levels), function(i) {
          values <- estimates[i, is.finite(estimates[i, ])]
          if (length(values) < 2L) return(rep(NA_real_, 3L))
          c(stats::sd(values), stats::quantile(values, probs = c(0.025, 0.975), names = FALSE))
        }, numeric(3L))
        wf$bootstrap <- tibble::tibble(
          Facet = fit$facets$Facet,
          Level = fit$facets$Level,
          Estimate = fit$facets$Estimate,
          ModelSE = fit$facets$ModelSE,
          BootstrapSE = intervals[1L, ],
          BootLower95 = intervals[2L, ],
          BootUpper95 = intervals[3L, ],
          BootstrapRequested = B,
          BootstrapAvailable = available
        )
        if (any(available < 2L)) {
          wf$bootstrap_reason <- "Some facet levels have fewer than two usable bootstrap refits."
        }
      }
    }
  }

  wf$elapsed_sec <- wf$elapsed_sec + (proc.time()[[3L]] - t0)
  wf
}































step_mfrm_gtheory <- function(wf,
                              g_facets = NULL,
                              nesting = NULL,
                              g_scores = NULL,
                              g_weights = NULL,
                              d_study_grid = NULL,
                              target_g = 0.80,
                              target_phi = 0.80,
                              g_data = NULL) {
  if (!inherits(wf, "tidymfrm")) stop("`wf` must be a `tidymfrm` object.", call. = FALSE)
  for (target in list(target_g, target_phi)) {
    if (!is.numeric(target) || length(target) != 1L || !is.finite(target) || target < 0 || target > 1) {
      stop("`target_g` and `target_phi` must each be a number between 0 and 1, such as 0.80.", call. = FALSE)
    }
  }
  t0 <- proc.time()[[3L]]
  sp <- wf$spec
  use_facets <- g_facets %||% sp$facets
  use_scores <- g_scores %||% sp$score

  gt <- mfrm_g_study(
    data = g_data %||% wf$data,
    scores = use_scores,
    person = sp$person,
    facets = use_facets,
    nesting = nesting,
    weights = g_weights,
    method = "urgenova",
    missing = "omit",
    d_study_grid = d_study_grid
  )
  gt$input_rows <- nrow(g_data %||% wf$data)
  gt$input_source <- if (is.null(g_data)) "ratings" else "supplied G-Study data"
  ds <- gt$d_study
  ds$MeetsTargetG <- ds$G >= target_g
  ds$MeetsTargetPhi <- ds$Phi >= target_phi

  rec_phi <- ds[ds$MeetsTargetPhi %in% TRUE, , drop = FALSE]
  if (nrow(rec_phi) > 0L) {
    rec_phi <- rec_phi[order(rec_phi$TotalRatingsPerPerson, -rec_phi$Phi), , drop = FALSE][1L, , drop = FALSE]
  }
  rec_g <- ds[ds$MeetsTargetG %in% TRUE, , drop = FALSE]
  if (nrow(rec_g) > 0L) {
    rec_g <- rec_g[order(rec_g$TotalRatingsPerPerson, -rec_g$G), , drop = FALSE][1L, , drop = FALSE]
  }

  gt$d_study <- tibble::as_tibble(ds)
  gt$variance_components <- tibble::as_tibble(gt$variance_components)
  gt$recommended_for_g <- rec_g
  gt$recommended_for_phi <- rec_phi
  gt$targets <- c(G = target_g, Phi = target_phi)

  wf$gtheory <- gt
  wf$gtheory_reason <- NULL
  wf$elapsed_sec <- wf$elapsed_sec + (proc.time()[[3L]] - t0)
  wf
}









































step_mfrm_diagnose <- function(wf,
                               mnsq_range = c(0.70, 1.30),
                               zstd_limit = 2.0,
                               bias_df_method = c("facets", "classical"),
                               dff_group = NULL,
                               dff_method = c("mml_lrt", "cml_lrt", "raju_welch", "mh", "logistic"),
                               cut_scores = NULL,
                               threads = 1L) {
  if (!inherits(wf, "tidymfrm") || is.null(wf$fit)) {
    stop("Run `step_mfrm_measure()` before `step_mfrm_diagnose()`.", call. = FALSE)
  }
  t0 <- proc.time()[[3L]]
  bias_df_method <- match.arg(bias_df_method)
  dff_method <- match.arg(dff_method)
  sp <- wf$spec
  fit <- wf$fit

  rf <- fit$facets[fit$facets$Facet == sp$rater_facet, , drop = FALSE]
  sev_sd <- if (nrow(rf) >= 2L) stats::sd(rf$Estimate, na.rm = TRUE) else 0.0

  flag <- rep("Acceptable", nrow(rf))
  action <- rep("Maintain current rubric calibration.", nrow(rf))
  for (i in seq_len(nrow(rf))) {
    inf_m <- rf$InfitMnSq[i]
    inf_z <- rf$InfitZSTD[i]
    est <- rf$Estimate[i]
    if (is.finite(inf_m) && inf_m > mnsq_range[2L] && (!is.finite(inf_z) || abs(inf_z) >= zstd_limit)) {
      flag[i] <- "Misfit (Unpredictable)"
      action[i] <- "Inspect unexpected rating rows for rubric drift or idiosyncratic criteria."
    } else if (is.finite(inf_m) && inf_m < mnsq_range[1L] && (!is.finite(inf_z) || abs(inf_z) >= zstd_limit)) {
      flag[i] <- "Muted / Centrality"
      action[i] <- "Review potential central-tendency or halo restriction across criteria."
    } else if (is.finite(est) && sev_sd > 0 && abs(est) > max(1.0, 2.0 * sev_sd)) {
      flag[i] <- if (est > 0) "Markedly Severe" else "Markedly Lenient"
      action[i] <- "Model adjusts additive severity via Fair(M); discuss anchor exemplars in training."
    }
  }
  rf$DiagnosticFlag <- flag
  rf$RecommendedAction <- action
  raters_tbl <- tibble::as_tibble(rf)

  obs_res <- .tidymfrm_observation_residuals(wf$data, sp$person, sp$facets, sp$score, fit)
  unexp <- obs_res[obs_res$UnexpectedFlag, , drop = FALSE]
  if (nrow(unexp) > 0L) {
    unexp <- unexp[order(-abs(unexp$StdResidual)), , drop = FALSE]
  }

  bias_res <- NULL
  bias_reason <- NULL
  if (length(sp$facets) >= 2L) {
    bias_attempt <- tryCatch(
      estimate_bias(fit, interaction_facets = "all", df_method = bias_df_method),
      error = identity
    )
    if (inherits(bias_attempt, "error")) bias_reason <- conditionMessage(bias_attempt) else bias_res <- bias_attempt
  }

  dff_tbl <- NULL
  if (!is.null(dff_group)) {
    dff_obj <- analyze_dff(
      fit = fit,
      facet = sp$rater_facet,
      group = dff_group,
      data = fit$encoding$retained_data,
      method = dff_method,
      df_method = bias_df_method,
      min_obs = 2L,
      threads = threads
    )
    dff_tbl <- tibble::as_tibble(dff_obj$table)
  }

  decisions_tbl <- if (length(cut_scores)) .tidymfrm_cut_decisions(fit, cut_scores) else NULL

  st <- fit$steps
  disordered <- FALSE
  if (nrow(st) >= 2L) {
    for (lvl in unique(st$ScaleLevel)) {
      sub_s <- st$Estimate[st$ScaleLevel == lvl & is.finite(st$Estimate)]
      if (length(sub_s) >= 2L && any(diff(sub_s) < 0, na.rm = TRUE)) disordered <- TRUE
    }
  }

  wf$diagnostics <- list(
    fit_basis = .tidymfrm_diagnostic_fit_basis(fit),
    report_statistics = .mfrm_diagnostic_report_statistics(wf),
    raters = raters_tbl,
    observations = obs_res,
    unexpected = unexp,
    bias = bias_res,
    bias_reason = bias_reason,
    dff = dff_tbl,
    decisions = decisions_tbl,
    step_disordered = disordered,
    mnsq_range = mnsq_range,
    zstd_limit = zstd_limit
  )
  wf$elapsed_sec <- wf$elapsed_sec + (proc.time()[[3L]] - t0)
  .tidymfrm_refresh(wf)
}














































mfrm_analysis <- function(data,
                     person,
                     facets,
                     score,
                     model = c("RSM", "PCM", "GPCM"),
                     method = c("MML", "JML", "JML_jeffreys", "CML"),
                     step_facet = NULL,
                     slope_facet = NULL,
                     slope_reference = NULL,
                     slope_type = c("full_predictor", "ability_only"),
                     slope_prior_sd = NULL,
                     occasion = NULL,
                     person_score = c("auto", "EAP", "WLE"),
                     rater_facet = NULL,
                     rating_min = NULL,
                     rating_max = NULL,
                     compare_methods = NULL,
                     population_formula = NULL,
                     person_data = NULL,
                     step_constraints = NULL,
                     category_support = NULL,
                     facet_formulas = NULL,
                     facet_data = NULL,
                     quad_points = 31L,
                     mml_integration = c("fixed", "adaptive"),
                     xtreme = 0.3,
                     warn_connectivity = TRUE,
                     nesting = NULL,
                     g_scores = NULL,
                     g_weights = NULL,
                     d_study_grid = NULL,
                     target_g = 0.80,
                     target_phi = 0.80,
                     mnsq_range = c(0.70, 1.30),
                     zstd_limit = 2.0,
                     bias_df_method = c("facets", "classical"),
                     dff_group = NULL,
                     dff_method = c("mml_lrt", "cml_lrt", "raju_welch", "mh", "logistic"),
                     cut_scores = NULL,
                     bootstrap = 0L,
                     threads = 1L,
                     gpcm_mml_identification = NULL,
                     g_facets = NULL,
                     expected_design = NULL,
                     rater_panels = NULL,
                     g_data = NULL) {
  wf <- mfrm_workflow(
    data = data,
    person = person,
    facets = facets,
    score = score,
    rater_facet = rater_facet,
    occasion = occasion,
    rating_min = rating_min,
    rating_max = rating_max,
    expected_design = expected_design,
    rater_panels = rater_panels
  )
  wf <- step_mfrm_measure(
    wf,
    model = model,
    method = method,
    step_facet = step_facet,
    slope_facet = slope_facet,
    slope_reference = slope_reference,
    slope_type = slope_type,
    slope_prior_sd = slope_prior_sd,
    occasion = occasion,
    person_score = person_score,
    compare_methods = compare_methods,
    population_formula = population_formula,
    person_data = person_data,
    step_constraints = step_constraints,
    category_support = category_support,
    facet_formulas = facet_formulas,
    facet_data = facet_data,
    quad_points = quad_points,
    mml_integration = mml_integration,
    xtreme = xtreme,
    warn_connectivity = warn_connectivity,
    bootstrap = bootstrap,
    threads = threads,
    gpcm_mml_identification = gpcm_mml_identification
  )
  wf <- tryCatch(
    step_mfrm_gtheory(
      wf,
      g_facets = g_facets,
      nesting = nesting,
      g_scores = g_scores,
      g_weights = g_weights,
      d_study_grid = d_study_grid,
      target_g = target_g,
      target_phi = target_phi,
      g_data = g_data
    ),
    mfrm_gstudy_unavailable = function(e) {
      wf$gtheory_reason <- conditionMessage(e)
      warning("G-Study unavailable: ", wf$gtheory_reason, call. = FALSE)
      wf
    }
  )
  wf <- step_mfrm_diagnose(
    wf,
    mnsq_range = mnsq_range,
    zstd_limit = zstd_limit,
    bias_df_method = bias_df_method,
    dff_group = dff_group,
    dff_method = dff_method,
    cut_scores = cut_scores,
    threads = threads
  )
  wf
}



tidymfrm <- mfrm_analysis


predict.tidymfrm <- function(object, newdata = NULL,
                             type = c("response", "probs", "link", "fair", "wle", "eap"), ...) {
  if (is.null(object$fit)) stop("Run `step_mfrm_measure()` before `predict()`.", call. = FALSE)
  stats::predict(object$fit, newdata = newdata, type = match.arg(type), ...)
}


coef.tidymfrm <- function(object, ...) {
  if (is.null(object$fit)) stop("Run `step_mfrm_measure()` before `coef()`.", call. = FALSE)
  stats::coef(object$fit, ...)
}


vcov.tidymfrm <- function(object, ...) {
  if (is.null(object$fit)) stop("Run `step_mfrm_measure()` before `vcov()`.", call. = FALSE)
  stats::vcov(object$fit, ...)
}


confint.tidymfrm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(object$fit)) stop("Run `step_mfrm_measure()` before `confint()`.", call. = FALSE)
  stats::confint(object$fit, parm = parm, level = level, ...)
}


anova.tidymfrm <- function(object, ...) {
  workflows <- c(list(object), list(...))
  if (length(workflows) < 2L || !all(vapply(workflows, inherits, logical(1), "tidymfrm"))) {
    stop("Supply at least two calibrated tidymfrm workflows to anova().", call. = FALSE)
  }
  fits <- lapply(workflows, `[[`, "fit")
  if (any(vapply(fits, is.null, logical(1)))) {
    stop("Run step_mfrm_measure() for every workflow before anova().", call. = FALSE)
  }
  labels <- make.unique(vapply(fits, `[[`, character(1), "model"))
  names(fits) <- labels
  fit_env <- list2env(fits, parent = environment())
  comparison <- do.call(stats::anova, lapply(labels, as.name), envir = fit_env)
  comparison$table$InferenceStatus <- vapply(fits, function(fit) summary(fit)$readiness$status, character(1))
  comparison
}


print.tidymfrm <- function(x, ...) {
  x <- .tidymfrm_refresh(x)
  cat("== mfrmr: End-to-End Many-Facet Assessment Workflow ==================\n")
  cat(sprintf(
    "* Design & Connectivity : %d ratings | %d persons | Subsets = %d (%s)\n",
    x$connectivity$n_obs, x$connectivity$n_person, x$connectivity$n_subsets,
    if (x$connectivity$connected) "MULTI-FACET GRAPH CONNECTED" else "DISCONNECTED - review links!"
  ))
  cat(paste0("  ", .tidymfrm_linkage_notes(x), "\n"), sep = "")
  if (!is.null(x$fit)) {
    fit_summary <- summary(x$fit)
    cat(sprintf(
      "* Pillar 1 [Measurement]: %s / %s (Score: %s) | Converged = %s | logLik = %.2f\n",
      .tidymfrm_model_label(x$fit, fit_summary), x$fit$method, x$fit$score_method, x$fit$converged, x$fit$loglik
    ))
    readiness <- fit_summary$readiness
    cat("  Inference:", readiness$status, paste(readiness$reason_codes, collapse = ", "), "\n")
    if (is.finite(x$fit$population$eap_reliability)) {
      cat(sprintf("  Latent EAP Reliability = %.4f (Pop SD = %.3f logits)\n",
                  x$fit$population$eap_reliability, x$fit$population$sigma))
    }
  }
  if (!is.null(x$gtheory)) {
    gs <- x$gtheory$g_study_summary
    cat(sprintf(
      "* Pillar 2 [G / D-Study]: Observed G = %.4f (RelSEM = %.3f) | Phi = %.4f (AbsSEM = %.3f)\n",
      gs$G, gs$RelSEM, gs$Phi, gs$AbsSEM
    ))
  } else if (!is.null(x$gtheory_reason)) {
    cat("* Pillar 2 [G / D-Study]: Unavailable.", x$gtheory_reason, "\n")
  }
  if (!is.null(x$diagnostics)) {
    rt <- x$diagnostics$raters
    n_flg <- .tidymfrm_flag_count(rt)
    cat("  Rater screens unavailable:", sum(rt$DiagnosticFlag == "Review required"), "\n")
    cat(sprintf(
      "* Pillar 3 [Rater Audit]: %d/%d raters flagged | %d unexpected ratings (|z| >= 2) | Steps %s\n",
      n_flg, nrow(rt), nrow(x$diagnostics$unexpected),
      if (x$diagnostics$step_disordered) "DISORDERED" else "ordered"
    ))
  }
  cat(sprintf("* Total Runtime         : %.4f sec\n", x$elapsed_sec))
  invisible(x)
}


summary.tidymfrm <- function(object, ..., reference_plan = NULL, fair_level = 0.95) {
  object <- .tidymfrm_refresh(object)
  fit_summary <- if (!is.null(object$fit)) summary(object$fit)
  list(
    readiness = fit_summary$readiness,
    metadata = .mfrm_report_metadata(object),
    report_tables = .mfrm_report_tables(object, reference_plan, fair_level),
    model_label = if (!is.null(fit_summary)) .tidymfrm_model_label(object$fit, fit_summary),
    connectivity = object$connectivity,
    planned_coverage = object$planned_coverage,
    measurement_facets = if (!is.null(object$fit)) tibble::as_tibble(object$fit$facets) else NULL,
    measurement_steps = if (!is.null(object$fit)) tibble::as_tibble(object$fit$steps) else NULL,
    measurement_slopes = if (!is.null(object$fit) && !is.null(object$fit$slopes)) tibble::as_tibble(object$fit$slopes) else NULL,
    drift = if (!is.null(object$drift)) tibble::as_tibble(as.data.frame(object$drift, component = "severity")) else NULL,
    explanatory = if (!is.null(object$fit) && !is.null(object$fit$explanatory)) tibble::as_tibble(object$fit$explanatory) else NULL,
    method_comparison = object$method_comparison,
    bootstrap = object$bootstrap,
    bootstrap_reason = object$bootstrap_reason,
    variance_components = if (!is.null(object$gtheory)) object$gtheory$variance_components else NULL,
    gtheory_reason = object$gtheory_reason,
    d_study = if (!is.null(object$gtheory)) object$gtheory$d_study else NULL,
    rater_diagnostics = if (!is.null(object$diagnostics)) object$diagnostics$raters else NULL,
    rater_category_counts = object$diagnostics$category_counts,
    pairwise_bias = if (!is.null(object$diagnostics$bias)) tibble::as_tibble(object$diagnostics$bias$pairwise) else NULL,
    unexpected_top10 = if (!is.null(object$diagnostics)) utils::head(object$diagnostics$unexpected, 10L) else NULL,
    dff = if (!is.null(object$diagnostics)) object$diagnostics$dff else NULL,
    decisions = if (!is.null(object$diagnostics)) object$diagnostics$decisions else NULL
  )
}
