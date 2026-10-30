# Internal two-owner representation. Effective slopes belong to observed
# crossings; component slopes retain their own facet and level labels. Only
# the first family has geometric mean one (the population is fixed N(0, 1)).
mfrm_has_product_slopes <- function(fit) {
  is.list(fit) && is.list(fit$config) && !is.null(fit$config$gpcm_spec$component_levels)
}

stop_if_product_slopes <- function(fit, helper) {
  stop_if_jml_adjustment(fit, helper)
  if (mfrm_has_product_slopes(fit)) {
    stop(new_gpcm_scope_error(paste0("`", helper, "` is not available for two slope families. ",
      "Use summary(fit) and mfrm_curve_intervals(fit, newdata) for fitted curves without intervals; ",
      "category_curves_report(fit, newdata = contexts) supplies full-context curves and cumulative thresholds without intervals. ",
      "confint(fit, parm = 'slopes') separately checks experimental component log-Wald intervals. ",
      "mfrm_response_diagnostics(fit) supplies descriptive posterior residuals with the fitted integration method. ",
      "An explicit method = 'profile' call can separately review one named component slope. ",
      "saved curves can be passed to mfrm_results(fit, include = c('fit', 'plots'), compute = 'never', intervals = curves)."),
      helper = helper, status = "blocked", area = "Two-family GPCM outputs"))
  }
  invisible(NULL)
}

build_gpcm_product_slope_spec <- function(prep, slope_facets, step_facet) {
  if (length(slope_facets) != 2L || anyDuplicated(slope_facets) ||
      !all(slope_facets %in% prep$facet_names)) {
    stop("Product slopes require two distinct observed facets.", call. = FALSE)
  }
  levels <- prep$levels[slope_facets]
  n <- lengths(levels)
  if (any(n < 2L)) stop("Product slopes require at least two levels per facet.")
  ids <- as.data.frame(unname(lapply(prep$data[slope_facets], as.integer)), check.names = FALSE)
  names(ids) <- slope_facets
  cells <- unique(ids)
  cells <- cells[order(cells[[1]], cells[[2]]), , drop = FALSE]
  rownames(cells) <- NULL
  nr <- nrow(cells)
  first <- which(cells[[1]] < n[1])
  reference <- which(cells[[1]] == n[1])
  design <- Matrix::sparseMatrix(
    i = c(first, rep(reference, n[1] - 1L), seq_len(nr)),
    j = c(cells[[1]][first], rep(seq_len(n[1] - 1L), each = length(reference)),
          n[1] - 1L + cells[[2]]),
    x = c(rep(1, length(first)), rep(-1, length(reference) * (n[1] - 1L)), rep(1, nr)),
    dims = c(nr, sum(n) - 1L))
  if (nrow(design) < ncol(design) ||
      as.integer(Matrix::rankMatrix(design, method = "qr")) < ncol(design)) {
    stop("Facet crossings do not identify both slope families.", call. = FALSE)
  }
  list(active = TRUE, slope_facet = slope_facets, step_facet = step_facet,
       component_levels = levels, cells = cells, log_slope_design = design,
       n_params = ncol(design), identification = "first_family_sum_zero_log_slopes",
       scale_reference = "fixed_standard_normal",
       reduction_reference = "PCM when all slopes equal 1")
}

gpcm_log_slope_design <- function(spec) {
  spec$log_slope_design %||% sum_zero_jacobian(length(spec$levels))
}

# Component identities and identification meanings used by tables and targets.
mfrm_gpcm_product_metadata <- function(spec) {
  n <- lengths(spec$component_levels)
  data.frame(Owner = rep(names(n), n),
    Level = unlist(spec$component_levels, use.names = FALSE),
    ScaleReference = rep(c("geometric_mean_one", "fixed_standard_normal"), n),
    Identification = rep(c("sum_to_zero_log_slopes",
      "free_given_first_family_and_population"), n))
}

# One row per declared owner; references describe the model, not fit readiness.
mfrm_gpcm_slope_roles <- function(config) {
  owners <- as.character(config$slope_facet %||% NA_character_)
  spec <- config$gpcm_spec
  if (is.null(spec$component_levels)) {
    if (length(owners) != 1L) {
      stop("Multiple slope owners require their fitted product specification.", call. = FALSE)
    }
    return(data.frame(SlopeOwner = owners, ScaleReference = "geometric_mean_one",
      Identification = "sum_to_zero_log_slopes"))
  }
  if (!identical(owners, names(spec$component_levels)) || length(owners) != 2L ||
      !identical(spec$identification, "first_family_sum_zero_log_slopes") ||
      !identical(spec$scale_reference, "fixed_standard_normal") ||
      !identical(config$method, "MML") || isTRUE(config$population_spec$active)) {
    stop("Slope summaries require the fitted owner order and scale constraints.", call. = FALSE)
  }
  roles <- unique(mfrm_gpcm_product_metadata(spec)[c("Owner", "ScaleReference", "Identification")])
  names(roles)[1] <- "SlopeOwner"
  rownames(roles) <- NULL
  roles
}

# Match within owners, without delimiter-based composite keys or assuming row order.
# Legacy single-owner tables do not carry SlopeOwner; they remain supported.
mfrm_match_slope_table <- function(owners, levels, table) {
  index <- rep(NA_integer_, length(levels))
  if (!length(levels) || !nrow(table)) return(index)
  if (is.null(table[["SlopeOwner"]]) && length(unique(owners)) != 1L) {
    stop("Multiple slope owners require an explicit SlopeOwner column.", call. = FALSE)
  }
  for (owner in unique(owners)) {
    rows <- if (is.null(table[["SlopeOwner"]])) seq_len(nrow(table)) else which(table[["SlopeOwner"]] == owner)
    if (anyDuplicated(table$SlopeFacet[rows])) {
      stop("Slope levels must be unique within each owner.", call. = FALSE)
    }
    selected <- which(owners == owner)
    index[selected] <- rows[match(levels[selected], table$SlopeFacet[rows])]
  }
  index
}

# Targets for the existing joint-information/delta-method calculation. These
# matrices describe uncertainty in estimated parameters, not G-study score
# covariance components. No interval eligibility is established here.
mfrm_gpcm_product_slope_targets <- function(par, config) {
  spec <- config$gpcm_spec
  if (is.null(spec$log_slope_design) ||
      !identical(config$model, "GPCM") || !identical(config$method, "MML") ||
      isTRUE(config$population_spec$active)) {
    stop("Product-slope targets currently require the fixed-standard-normal MML model.")
  }
  sizes <- build_param_sizes(config)
  if (length(par) != sum(unlist(sizes)) || any(!is.finite(par))) {
    stop("Product-slope targets require the complete finite free-parameter vector.")
  }
  invisible(expand_params(par, sizes, config))
  n <- lengths(spec$component_levels)
  component_design <- as.matrix(Matrix::bdiag(
    sum_zero_jacobian(n[1]), diag(n[2])))
  metadata <- mfrm_gpcm_product_metadata(spec)
  contexts <- as.data.frame(unname(Map(function(ids, levels) levels[ids],
                                spec$cells, spec$component_levels)), check.names = FALSE)
  # data.frame list-name expansion can transliterate names in a C locale.
  names(contexts) <- names(spec$component_levels)
  make_target <- function(design, labels, table, description) {
    jac <- matrix(0, nrow(design), length(par))
    jac[, build_param_slices(sizes)[["log_slopes"]]] <- as.matrix(design)
    value <- as.numeric(jac %*% par)
    list(value = value, estimate = exp(value), jacobian = jac,
         labels = labels, log_scale = TRUE, metadata = table,
         scale = "fixed_standard_normal", target = description)
  }
  list(components = make_target(component_design,
         paste(metadata$Owner, metadata$Level, sep = ":"), metadata,
         "Facet slope components on the declared fixed population scale"),
       effective = make_target(spec$log_slope_design,
         paste0("crossing:", seq_len(nrow(contexts))), contexts,
         "Effective slope products for observed crossings"))
}

project_gpcm_log_slope_gradient <- function(gradient, spec) {
  if (is.null(spec$log_slope_design)) return(project_sum_zero_gradient(gradient))
  as.numeric(Matrix::crossprod(spec$log_slope_design, gradient))
}

collapse_gpcm_log_slopes <- function(log_slopes, spec) {
  if (is.null(spec$log_slope_design)) return(log_slopes[seq_len(spec$n_params %||% 0L)])
  design <- spec$log_slope_design
  if (length(log_slopes) != nrow(design) || any(!is.finite(log_slopes))) {
    stop("Effective log slopes must match the fitted crossings.", call. = FALSE)
  }
  free <- as.numeric(Matrix::qr.coef(Matrix::qr(design), log_slopes))
  if (max(abs(as.numeric(design %*% free) - log_slopes)) >
      1e-10 * max(1, abs(log_slopes))) {
    stop("Effective log slopes do not follow the fitted product structure.", call. = FALSE)
  }
  free
}

# Called only after mfrm_gmfrm_problem() validates the limited model contract.
# Public dispatch preserves this model specification through its admitted outputs.
mfrm_gmfrm_common_setup <- function(data, max_score,
                                    slope_facets = c("Task", "Rater"),
                                    person = "Person", score = "Score") {
  data[c(person, slope_facets)] <- lapply(data[c(person, slope_facets)], as.character)
  prep <- prepare_mfrm_data(data, person, slope_facets, score,
                            rating_min = 0, rating_max = max_score, keep_original = TRUE)
  step_facet <- slope_facets[2L]
  built <- build_estimation_config(prep, "GPCM", "MML", step_facet, slope_facets[1L],
    weight_col = NULL, facet_signs = setNames(rep(-1, 2L), slope_facets),
    positive_facets = character(), noncenter_facet = step_facet,
    dummy_facets = character(), anchor_df = NULL, group_anchor_df = NULL)
  config <- built$config
  config$slope_facet <- slope_facets
  config$gpcm_spec <- build_gpcm_product_slope_spec(prep, slope_facets, step_facet)
  identity <- mfrmr_gpcm_model_identity("GPCM", config$gpcm_spec)
  config$gpcm_model_family <- identity$model_family
  config$gpcm_slope_composition <- identity$slope_composition
  config$gpcm_mml_identification <- "fixed_population_first_family_gm1"
  config$gpcm_common_discrimination <- "free_second_family"
  sizes <- build_param_sizes(config)
  idx <- build_indices(prep, step_facet, slope_facets, gpcm_spec = config$gpcm_spec)
  list(prep = prep, config = config, sizes = sizes, idx = idx,
       parameter_map = mfrmr_mml_optimizer_parameter_map(prep, idx, config, sizes))
}

# Qualification for a local Wald approximation, not a global boundary or
# sampling-coverage certificate. Reuse the full marginal information: frozen-Q
# curvature and the inverse slope block are not the covariance of joint estimates.
mfrm_gpcm_product_inference <- function(fit) {
  reused <- mfrmr_gmfrm_inference_reuse(fit)
  if (!is.null(reused)) return(reused)
  checks <- data.frame(Check = character(), Passed = logical(), Detail = character())
  out <- list(covariance = NULL, information = NULL, checks = checks,
    check = list(eligible = FALSE, review = "Two-family inference was not evaluated."))
  record <- function(name, passed, detail) {
    out$checks <<- rbind(out$checks, data.frame(Check=name, Passed=isTRUE(passed), Detail=detail))
    isTRUE(passed)
  }
  refuse <- function() {
    out$check$review <- paste(out$checks$Detail[!out$checks$Passed], collapse=" ")
    out
  }
  config <- fit$config
  adaptive <- mfrmr_adaptive_integration(config)
  if (!record("Scope", mfrm_has_product_slopes(fit) &&
      identical(config$method, "MML") && !isTRUE(config$population_spec$active) &&
      identical(config$posterior_basis, "fixed_standard_normal") &&
      identical(config$estimation_control$mml_engine_used, if (adaptive) "direct" else "em") &&
      is.null(config$weight_col) && all(fit$prep$data$Weight == 1) &&
      identical(config$noncenter_facet, config$slope_facet[2]),
      "Two-family intervals require unweighted fixed-N(0,1) MML with fixed-grid EM or adaptive direct integration.")) return(refuse())
  specification <- fit$gmfrm$specification
  reference <- tryCatch(do.call(mfrm_gmfrm_common_setup,
    specification[setdiff(names(specification),"quadrature")]), error=function(e) NULL)
  without_factor_cache <- function(spec) {
    # Sparse QR factors can differ across Matrix/SuiteSparse builds. The
    # design entries and all model metadata still require exact agreement.
    if (inherits(spec$log_slope_design, "dgCMatrix"))
      spec$log_slope_design@factors <- list()
    spec
  }
  if (!record("Model identity", !is.null(reference) &&
      identical(lapply(reference$prep$data, identity),lapply(fit$prep$data, identity)) &&
      identical(reference$config$facet_specs,config$facet_specs) &&
      identical(without_factor_cache(reference$config$gpcm_spec),
                without_factor_cache(config$gpcm_spec)) &&
      identical(reference$config$facet_signs,config$facet_signs) &&
      identical(reference$config$interaction_specs,config$interaction_specs) &&
      identical(reference$config$n_cat,config$n_cat),
      "The saved data, constraints and product specification must match the original fitted model.")) return(refuse())
  roles <- tryCatch(mfrm_gpcm_slope_roles(config), error=function(e) NULL)
  target <- tryCatch(mfrm_gpcm_product_slope_targets(fit$opt$par,config)$components,
    error=function(e) NULL)
  index <- if (!is.null(target)) tryCatch(mfrm_match_slope_table(target$metadata$Owner,
    target$metadata$Level, fit$slopes), error=function(e) integer()) else integer()
  if (!record("Slope identity", !is.null(roles) && !is.null(target) &&
      length(index) == nrow(fit$slopes) && !anyNA(index) &&
      isTRUE(all.equal(target$estimate, unname(fit$slopes$Estimate[index]), tolerance=1e-10)),
      "Slope owners, levels, scale constraints and retained estimates must agree.")) return(refuse())
  if (adaptive) {
    if (!record("Direct convergence", isTRUE(fit$summary$Converged) &&
        identical(fit$opt$convergence, 0L) &&
        identical(fit$opt$optimizer_diagnostics$ConvergenceBasis, "optimizer_gradient") &&
        identical(fit$opt$optimizer_diagnostics$ConvergenceSeverity, "pass"),
        "Direct adaptive MML must pass its optimizer and gradient convergence checks.")) return(refuse())
  } else if (!record("EM convergence", isTRUE(fit$summary$Converged) &&
      identical(fit$opt$optimizer_diagnostics$ConvergenceBasis, "marginal_score_per_person"),
      "EM must meet its per-Person marginal-score stopping rule.")) return(refuse())
  category <- audit_mfrm_category_support(fit$prep,config,build_param_sizes(config))
  if (!record("Category support", identical(category$readiness$CategoryState, "adequate"),
      "Category support must be adequate for a local Wald approximation; inspect the declared score ladder.")) return(refuse())
  info <- compute_mml_parameter_covariance(fit)
  out$information <- info
  if (!record("Joint information", identical(info$status,"ok") && !isTRUE(info$regularized) &&
      !is.null(info$cov) && all(is.finite(info$cov)),
      paste("An unregularized inverse of the full marginal information is required.",info$detail))) return(refuse())
  ev <- info$solution_information$evaluation_summary
  if (!record("Likelihood identity", nrow(ev)==1L && is.finite(fit$opt$value) &&
      is.finite(ev$ReevaluatedObjective) &&
      abs(ev$ReevaluatedObjective-fit$opt$value) <= 1e-10*max(1,abs(fit$opt$value)),
      "The retained likelihood must match reevaluation at the saved data and parameters.")) return(refuse())
  p <- ncol(info$cov); n <- length(fit$prep$levels$Person)
  budget <- getOption("mfrmr.max_information_bytes",256*1024^2)
  if (!record("Score workspace", 4*8*as.double(n)*p <= budget,
      "Person-score derivative workspace must fit within mfrmr.max_information_bytes.")) return(refuse())
  scores <- tryCatch(mfrm_mml_person_scores_numeric(fit),error=function(e) e)
  fine <- tryCatch(mfrm_mml_person_scores_numeric(fit,relative_step=5e-6),error=function(e) e)
  if (inherits(scores,"error") || inherits(fine,"error")) {
    record("Score derivatives",FALSE,"Person marginal-score derivatives could not be verified.")
    return(refuse())
  }
  column_scale <- sqrt(colSums(scores^2))
  scaled <- sweep(scores,2,pmax(column_scale,.Machine$double.eps),"/")
  scaled_fine <- sweep(fine,2,pmax(column_scale,.Machine$double.eps),"/")
  rank <- mfrmr_transformation_rank_ladder(scaled)
  fine_rank <- mfrmr_transformation_rank_ladder(scaled_fine)
  out$score_rank <- rank$rank_ladder
  score_rank <- if (rank$valid) min(rank$rank_ladder$Rank) else NA_integer_
  fine_score_rank <- if (fine_rank$valid) min(fine_rank$rank_ladder$Rank) else NA_integer_
  rank_deficient <- min(score_rank, fine_score_rank) < p
  derivative_error <- max(abs(scaled-scaled_fine))
  q <- config$estimation_control$quad_points
  out$numerical_checks <- data.frame(QuadraturePoints=q,ComparisonPoints=2L*q-1L,
    Integration=if (adaptive) "adaptive" else "fixed",
    FreeParameters=p,ScoreRank=score_rank,FineScoreRank=fine_score_rank,
    ScoreRankDeficient=rank_deficient,DerivativeDifference=derivative_error,
    MaximumMeanScore=NA_real_,CurvatureScaledGradient=NA_real_,InverseResidual=NA_real_,
    QuadratureScoreShift=NA_real_,QuadratureCovarianceChange=NA_real_)
  # Empirical score rank governs an OPG/sandwich covariance, not the inverse
  # observed Hessian used here. Retain it as a diagnostic, not an N > p veto.
  if (!record("Score derivative stability", rank$valid && fine_rank$valid &&
      is.finite(derivative_error) && derivative_error <= 1e-6,
      "Person marginal-score derivatives must be finite and agree within 1e-6 after column scaling at both derivative steps.")) return(refuse())
  gradient <- -colSums(scores)
  scaled_gradient <- sqrt(max(0,drop(crossprod(gradient,info$cov%*%gradient))))
  tol <- if (adaptive) 1e-6 else min(fit$gmfrm$controls$score_tol %||% 1e-6,1e-6)
  inverse_error <- norm(info$hessian%*%info$cov-diag(p),"I")
  out$numerical_checks$MaximumMeanScore <- max(abs(gradient))/n
  out$numerical_checks$CurvatureScaledGradient <- scaled_gradient
  out$numerical_checks$InverseResidual <- inverse_error
  out$numerical_checks$OptimizationCaution <- scaled_gradient > 1e-4
  if (!record("Stationary local maximum", is.finite(scaled_gradient) &&
      max(abs(gradient))/n <= tol && scaled_gradient <= .01 && inverse_error <= 1e-6,
      paste0("The fresh per-Person score must meet ", if (adaptive) "1e-6" else "min(em_score_tol, 1e-6)",
        ", standardized Newton displacement must be at most 0.01, and inverse residual at most 1e-6. Displacement above 1e-4 retains a warning."))) return(refuse())
  refined_fit <- fit; refined_fit$config$estimation_control$quad_points <- 2L*q-1L
  refined <- tryCatch(compute_mml_parameter_covariance(refined_fit),error=function(e) NULL)
  refined_scores <- tryCatch(mfrm_mml_person_scores_numeric(refined_fit),error=function(e) NULL)
  if (!record("Refined quadrature information", !is.null(refined) && identical(refined$status,"ok") &&
      !isTRUE(refined$regularized) && !is.null(refined_scores),
      "Higher-order quadrature must supply positive unregularized marginal information and finite scores.")) return(refuse())
  root <- chol(info$hessian)
  gradient_change <- -colSums(refined_scores)-gradient
  shift <- sqrt(max(0,drop(crossprod(gradient_change,info$cov%*%gradient_change))))
  change <- norm(root%*%(refined$cov-info$cov)%*%t(root),"2")
  out$numerical_checks$QuadratureScoreShift <- shift
  out$numerical_checks$QuadratureCovarianceChange <- change
  if (!record("Quadrature sensitivity", is.finite(shift) && is.finite(change) &&
      shift <= .01 && change <= .01,
      "The q versus 2q-1 comparison must be at most 0.01 in standardized score displacement and covariance change; otherwise refit with more quadrature points.")) return(refuse())
  out$covariance <- info$cov
  optimization_caution <- if (scaled_gradient > 1e-4) sprintf(
    "Small optimization residual: standardized Newton displacement is %.3g (above 1e-4, at most 0.01). The intervals retain this local numerical approximation.",
    scaled_gradient) else ""
  rank_caution <- if (rank_deficient) sprintf(
    "Observed Person-score rank is deficient (%d of %d free parameters). These model-based intervals use the observed Hessian, not an empirical score covariance; full-parameter sandwich inference is unavailable and finite-sample coverage is unqualified.",
    min(score_rank, fine_score_rank), p) else ""
  out$check <- list(eligible=TRUE,
    review="Experimental local Wald approximation from full joint marginal information; sampling coverage is not qualified.",
    caution=trimws(paste("Experimental two-family intervals: local numerical checks passed, but global identification, boundary absence and sampling coverage are not established.",
      optimization_caution,rank_caution,paste(mfrm_mml_information_caution(info),collapse=" "))))
  out
}
