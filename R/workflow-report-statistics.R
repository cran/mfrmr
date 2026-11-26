.mfrm_report_metadata <- function(wf) {
  fit <- wf$fit
  if (is.null(fit)) return(list(notes = character()))
  enc <- fit$encoding
  control <- fit$control
  bias <- switch(fit$method, JML_jeffreys = "Jeffreys/Firth penalized joint likelihood",
                 JML = "No explicit calibration bias correction", "Not a joint calibration")
  method_note <- switch(fit$method,
    JML = "Joint person and facet estimation can expand the measure scale with few ratings per person. Extreme-score handling affects person measures.",
    JML_jeffreys = "Jeffreys-penalized joint calibration reduces first-order bias under its model assumptions.",
    MML = "Marginal likelihood integrates over the specified ability population. EAP scores shrink toward that fitted population, conditional on the calibration.",
    CML = "Conditional likelihood conditions on person totals to remove ability from facet calibration.")
  fitting <- data.frame(
    Model = summary(fit)$model_label, Estimator = fit$method, PersonScoring = fit$score_method,
    CalibrationBiasReduction = bias,
    FacetL2Penalty = control$facet_l2_pen %||% NA_real_,
    SlopeL2Penalty = control$slope_l2_pen %||% NA_real_,
    Population = if (fit$method == "MML") fit$population$mode %||% NA_character_ else NA_character_,
    QuadraturePoints = if (fit$method == "MML") control$quad_points %||% NA_integer_ else NA_integer_,
    AdaptiveQuadrature = if (fit$method == "MML") control$adaptive_quad %||% NA else NA,
    StepFacet = enc$step_facet %||% "Common",
    SlopePlacement = if (identical(fit$model, "GPCM")) enc$slope_type else NA_character_,
    Centering = if (is.null(enc$centering_scope)) NA_character_ else
      if (identical(enc$centering_scope, 1L)) "facet" else "remainder",
    RatingRows = fit$n_obs, OmittedRatingRows = nrow(wf$data) - fit$n_obs,
    Converged = isTRUE(fit$converged), InferenceStatus = .mfrm_readiness(fit)$status
  )
  list(calibration = fitting,
       notes = c(paste("Estimator:", fit$method, "| Person scoring:", fit$score_method, "|", bias),
                 method_note,
                 "Optimizer convergence and inference readiness are reported separately. Person scoring and calibration bias correction refer to different parameters."))
}

.mfrm_displacement_report <- function(fit) {
  enc <- fit$encoding
  anchors <- lapply(seq_along(enc$base_facets), function(j) {
    idx <- which(enc$anchor_fixed[[j]] == 1L)
    data.frame(Facet = rep(enc$base_facets[j], length(idx)),
               Level = enc$facet_levels[[j]][idx], Type = rep("Fixed element", length(idx)),
               Target = enc$anchor_val[[j]][idx])
  })
  specs <- lapply(seq_along(anchors), function(j) as.list(anchors[[j]]$Level))
  anchors <- do.call(rbind, anchors)
  levels_by_row <- unlist(specs, recursive = FALSE)
  if (NROW(enc$group_anchors)) {
    groups <- enc$group_anchors
    col <- function(choices) intersect(choices, names(groups))[1L]
    group_facet <- as.character(groups[[col(c("Facet", "facet"))]])
    group_name <- as.character(groups[[col(c("Group", "group"))]])
    group_level <- as.character(groups[[col(c("Level", "level", "Element", "element"))]])
    group_value <- groups[[col(c("GroupValue", "Value", "Anchor", "value", "anchor"))]]
    for (rows in split(seq_len(nrow(groups)), .mfrm_mvgt_group(lapply(list(group_facet, group_name), factor)))) {
      j <- match(group_facet[rows[1L]], enc$base_facets)
      if (is.na(j) || is.null(enc$group_anchor_meta[[enc$base_facets[j]]])) next
      present <- intersect(group_level[rows], enc$facet_levels[[j]])
      if (!length(present)) next
      anchors <- rbind(anchors, data.frame(Facet = enc$base_facets[j], Level = group_name[rows[1L]],
                                          Type = "Group mean", Target = as.numeric(group_value[rows[1L]])))
      levels_by_row <- c(levels_by_row, list(present))
    }
  }
  if (!nrow(anchors)) return(NULL)
  ready <- .mfrm_readiness(fit)$ready
  p <- if (ready) predict(fit, type = "probs")
  cats <- if (ready) as.numeric(colnames(p))
  expected <- if (ready) drop(p %*% cats)
  variance <- if (ready) pmax(0, drop(p %*% cats^2) - expected^2)
  anchors$Rows <- integer(nrow(anchors))
  anchors$Displacement <- anchors$ConditionalSE <- NA_real_
  anchors$Status <- if (ready) "no_usable_conditional_information" else "calibration_unavailable"
  for (row in seq_len(nrow(anchors))) {
    j <- match(anchors$Facet[row], enc$base_facets)
    r <- which(enc$facet_levels[[j]][enc$facet_mat0[, j] + 1L] %in% levels_by_row[[row]])
    anchors$Rows[row] <- length(r)
    # This is a local diagnostic direction with all other measures held fixed.
    if (!ready || !length(r)) next
    if (anchors$Type[row] == "Group mean" && any(enc$anchor_fixed[[j]][match(levels_by_row[[row]], enc$facet_levels[[j]])] == 1L)) {
      anchors$Status[row] <- "group_direction_overlaps_fixed_elements"
      next
    }
    if (length(enc$interaction_meta) || length(enc$explanatory_meta)) {
      anchors$Status[row] <- "conditional_direction_not_supported"
      next
    }
    alpha <- if (identical(enc$slope_formula, 0L)) fit$alpha_obs[r] else rep(1, length(r))
    covariate <- enc$facet_signs[j] * alpha
    info <- sum(enc$weight[r] * covariate^2 * variance[r])
    if (!is.finite(info) || info <= 0) next
    anchors$Displacement[row] <- sum(enc$weight[r] * covariate *
      (enc$score0[r] + fit$rating_min - expected[r])) / info
    anchors$ConditionalSE[row] <- 1 / sqrt(info)
    anchors$Status[row] <- "conditional_one_step"
  }
  anchors
}

.mfrm_diagnostic_report_statistics <- function(wf) {
  fit <- wf$fit
  if (is.null(fit)) return(NULL)
  sp <- wf$spec
  context <- unique(c(setdiff(sp$facets, sp$rater_facet), sp$occasion))
  agreement <- tryCatch(
    interrater_agreement(fit, rater = sp$rater_facet, context_facets = context, model_expected = TRUE),
    error = identity
  )
  checked <- .tidymfrm_check_diagnostics(wf)
  reliability <- checked$diagnostics$report_statistics$reliability
  if (is.null(reliability)) reliability <- diagnose_mfrm(fit)$reliability
  if (NROW(reliability))
    reliability$HomogeneityBasis <- "ModelSE diagonal; no covariance adjustment"
  list(reliability = reliability,
       agreement = if (!inherits(agreement, "error")) agreement,
       agreement_reason = if (inherits(agreement, "error")) conditionMessage(agreement),
       displacement = .mfrm_displacement_report(fit))
}

.mfrm_multivariate_report <- function(g) {
  comp <- g$vcov_components
  if (!length(comp)) return(list())
  min_eigen <- function(x) if (all(is.finite(x)))
    min(eigen(x, symmetric = TRUE, only.values = TRUE)$values) else NA_real_
  valid <- function(x) is.finite(min_eigen(x)) && min_eigen(x) >= -1e-8 * max(abs(x))
  covariance <- do.call(rbind, lapply(names(comp), function(nm) {
    x <- comp[[nm]]
    idx <- which(lower.tri(x, diag = TRUE), arr.ind = TRUE)
    data.frame(Component = nm, ScoreA = rownames(x)[idx[, 1L]], ScoreB = colnames(x)[idx[, 2L]],
               Estimate = x[idx], MinEigenvalue = min_eigen(x),
               CovarianceStatus = if (valid(x)) "positive_semidefinite" else "invalid_covariance")
  }))
  observed_n <- g$obs_n_per_facet
  scales <- vapply(g$non_person_subsets, function(i) 1 / prod(observed_n[i]), numeric(1))
  w <- g$weights
  effective <- do.call(rbind, lapply(c("Universe score", "Relative error", "Absolute error"), function(target) {
    i <- switch(target, "Universe score" = 1L, "Relative error" = which(g$variance_components$ContainsPerson)[-1L],
                "Absolute error" = seq_along(comp)[-1L])
    v <- Reduce(`+`, Map(function(x, s) x * s, comp[i], scales[i]), init = matrix(0, length(w), length(w)))
    denominator <- drop(crossprod(w, v %*% w))
    status <- if (any(g$variance_components$Truncated[i])) "scalar_component_truncated" else
      if (!valid(v)) "non_psd_target_covariance" else
      if (!is.finite(denominator) || denominator <= 0) "zero_or_undefined_target_variance" else "available"
    weight <- if (identical(status, "available")) w * drop(v %*% w) / denominator else rep(NA_real_, length(w))
    data.frame(Target = target, Score = g$scores, NominalWeight = w,
               EffectiveWeight = weight, TargetVariance = denominator, Status = status)
  }))
  list("Subscore variance-covariance components" = covariance,
       "Subscore generalizability" = g$univariate,
       "Effective weights at observed rating counts" = effective)
}

.mfrm_report_tables <- function(wf, reference_plan = NULL, fair_level = 0.95) {
  fit <- wf$fit
  dg <- wf$diagnostics
  stat <- dg$report_statistics
  if (is.null(stat) && !is.null(dg)) stat <- .mfrm_diagnostic_report_statistics(wf)
  extreme <- if (!is.null(fit)) {
    flag <- fit$person$ExtremeFlag
    n <- sum(!is.na(flag))
    count <- c(sum(flag == -1L, na.rm = TRUE), sum(flag == 1L, na.rm = TRUE))
    data.frame(Endpoint = c("Minimum", "Maximum"), Persons = count, TotalPersons = n,
               Proportion = if (n) count / n else NA_real_)
  }
  network <- wf$connectivity$rater_linkage$robustness
  fair <- if (!is.null(reference_plan)) {
    if (is.null(fit)) stop("Attach a calibration before supplying reference_plan.", call. = FALSE)
    fair_average_table(fit, reference_plan = reference_plan, fair_se = TRUE, ci_level = fair_level)
  }
  c(list("Separation, strata, and fixed-effect homogeneity" = stat$reliability,
         "Observed inter-rater agreement" = stat$agreement$summary,
         "Expected agreement and model-based kappa" = stat$agreement$model_expected,
         "Agreement availability" = if (length(stat$agreement_reason)) data.frame(Reason = stat$agreement_reason),
         "Extreme person scores" = extreme,
         "Anchor displacement" = stat$displacement,
         "Linked benchmark comparison" = wf$anchor_drift,
         "Rating network robustness" = network$overview,
         "Panel effective resistance" = network$panel_resistance,
         "Observed-design generalizability" = wf$gtheory$g_study_summary,
         "Fair Average reference plan" = fair$reference_plan,
         "Common-plan Fair Average" = fair$reference_scores),
    .mfrm_multivariate_report(wf$gtheory))
}
