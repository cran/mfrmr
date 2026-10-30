# At each reachable (rating, total) node, compare incoming paths with one
# canonical prefix. Each difference extends to an observed total; telescoping
# these differences spans every pair of conditional response patterns.
mfrmr_cml_support_design <- function(problem) {
  p <- problem$parameters
  if (as.double(p) * p > 2e6) return(NULL)
  rows <- list()
  for (g in problem$groups) {
    totals <- which(g$counts > 0) - 1L
    m <- length(g$features); K <- nrow(g$features[[1L]]) - 1L
    if (as.double(max(totals) + 1L) * p > 2e6) return(NULL)
    previous <- matrix(0, 1L, p)
    for (j in seq_len(m)) {
      current <- matrix(NA_real_, min(j * K, max(totals)) + 1L, p)
      for (s in seq_len(nrow(current)) - 1L) {
        if (!any(totals >= s & totals <= s + (m - j) * K)) next
        for (k in seq.int(max(0L, s - nrow(previous) + 1L), min(K, s))) {
          prefix <- previous[s - k + 1L, ]
          if (anyNA(prefix)) next
          candidate <- prefix + g$features[[j]][k + 1L, ]
          if (anyNA(current[s + 1L, ])) current[s + 1L, ] <- candidate else {
            delta <- candidate - current[s + 1L, ]
            if (any(delta != 0)) {
              if ((as.double(length(rows)) + 1) * p > 2e6) return(NULL)
              rows[[length(rows) + 1L]] <- delta
            }
          }
        }
      }
      previous <- current
    }
  }
  if (length(rows)) unique(do.call(rbind, rows)) else matrix(0, 1L, p)
}

mfrmr_fit_cml_comparisons <- function(prep, idx, cfg, sizes, review, a) {
  if ("duplicate_person_facet_cells" %in% prep$preparation_notes$Condition)
    stop("CML intervals require resolution of repeated rating cells; correct duplicates or explicitly declare repeated_ratings = 'independent' when justified.", call. = FALSE)
  table <- review$table
  table$EstimateStatus <- "unavailable"
  table$EstimateReason <- "The declared target is unidentified or its conditional support was not evaluated."
  table$SE <- table$Lower <- table$Upper <- NA_real_
  table$NumericalState <- table$BoundaryState <- "not_evaluated"
  table$QuadratureState <- "not_applicable"
  table$GradientMaxAbs <- NA_real_
  table$InferenceMethod <- "conditional_information_wald"
  table$ConfidenceLevel <- .95
  eligible <- table$ConditionalState %in% c("target_identified_only", "all_parameters_identified")
  covariance <- matrix(NA_real_, nrow(table), nrow(table),
    dimnames = list(table$Comparison, table$Comparison))
  problem <- mfrmr_cml_problem(prep, idx, cfg, sizes)
  checks <- list(boundary = NULL, numerical = NULL)
  quotient <- fitted <- NULL
  notes <- c("Only the declared identifiable interaction contrasts are reported; unidentified rows and covariance entries remain NA.",
    "Optimization removes conditional-support null directions and preserves the response law given each Person total and roster; internal coordinates are not identified cell effects.",
    "Pointwise 95% Wald intervals use the full conditional information on the identified subspace, without regularization. They require the correctly specified response model, independent Persons, a stable fixed-dimensional identified subspace and information growing in every direction of that subspace; finite-sample coverage is not guaranteed.",
    "More Persons alone does not establish the information-growth condition. Response-dependent assignment, stopping or missingness requires separate justification that the conditional response law is preserved.",
    "A recession direction anywhere in the conditioned model withholds all estimates; extended-MLE inference is not implemented.",
    "Person-constant main effects disappear. CML estimates neither Person abilities nor an ability population SD, so these contrasts are on the original model scale.")
  finish <- function() {
    out <- tibble::as_tibble(table)
    attr(out, "comparison_review") <- review
    attr(out, "fit_checks") <- checks
    attr(out, "notes") <- notes
    structure(list(comparisons = out, covariance = covariance, comparison_review = review,
      checks = checks, representation = quotient, optimization = fitted,
      person_totals = problem$person_totals, prep = prep, config = cfg,
      controls = a[c("maxit", "reltol", "optimizer")], notes = notes),
      class = "mfrm_comparison_fit")
  }
  if (!any(eligible) || any(table$ConditionalState == "not_evaluated")) return(finish())
  support <- mfrmr_cml_support_design(problem)
  quotient <- mfrmr_comparison_quotient(list(design = support, map = problem$map), sizes, table$Rank[1L])
  T <- quotient$transform
  reduced <- problem
  reduced$parameters <- ncol(T)
  reduced$observed <- drop(crossprod(T, problem$observed))
  reduced$groups <- lapply(problem$groups, function(g) {
    g$features <- lapply(g$features, function(x) x %*% T)
    g
  })
  zero <- mfrmr_cml_evaluate(problem, numeric(problem$parameters))
  checks$boundary <- tryCatch(mfrmr_cml_boundary(problem, zero$gradient),
    error = function(e) list(state = "not_evaluated", reason = conditionMessage(e)))
  table$BoundaryState <- checks$boundary$state
  if (checks$boundary$state != "no_recession_detected") {
    table$EstimateReason[eligible] <- "CML recession was detected or could not be checked; no finite comparison estimate or interval delivered."
    return(finish())
  }
  initial <- build_initial_param_vector(cfg, sizes)
  start <- drop(crossprod(quotient$scaled_basis,
    quotient$scales * initial[quotient$linear_indices]))
  fitted <- tryCatch(mfrmr_cml_optimize(reduced, start, a), error = identity)
  if (inherits(fitted, "error")) {
    checks$numerical <- conditionMessage(fitted)
    fitted <- NULL
    table$NumericalState <- "review"
    table$EstimateReason[eligible] <- checks$numerical
    return(finish())
  }
  checks$numerical <- fitted$checks
  C <- review$target_map[eligible, , drop = FALSE] %*% quotient$linear_basis
  estimates <- drop(C %*% fitted$opt$par)
  wald <- mfrmr_cml_wald(estimates, C, fitted$covariance, .95, fitted$checks$information_reason)
  for (column in c("Estimate", "SE", "Lower", "Upper")) table[[column]][eligible] <- wald$table[[column]]
  covariance[eligible, eligible] <- wald$covariance
  information_available <- identical(fitted$checks$information_status, "available")
  table$NumericalState <- if (information_available) "ready" else "review"
  table$GradientMaxAbs <- fitted$checks$gradient_max
  table$FormalInferenceEligible[eligible] <- wald$table$CIEligible
  table$EstimateStatus[eligible] <- if (information_available) "conditional_wald" else "conditional_point"
  table$EstimateReason[eligible] <- if (information_available)
    "Identified conditional contrast; finite-MLE, stationary-point and unregularized information checks passed." else
    wald$table$InferenceReview
  if (!information_available) notes <- c(notes, paste("Point estimates retained; uncertainty unavailable:",
    fitted$checks$information_reason))
  finish()
}

# Reuse qualified target covariance; changing the confidence level never refits.
mfrmr_cml_comparison_intervals <- function(fit, facet, contrasts, method, clusters,
                                         adjust, level, scale, comparison_review) {
  if (!identical(fit$config$method, "CML"))
    stop("Comparison-only JML/MML fits remain descriptive; only CML supplies comparison intervals.", call. = FALSE)
  if (scale != "native")
    stop("CML does not estimate a population SD; standardized contrasts require an MML fit that estimates that SD.", call. = FALSE)
  if (method != "model" || !is.null(clusters) || !identical(adjust, FALSE))
    stop("CML comparison intervals require method = 'model', clusters = NULL and adjust = FALSE.", call. = FALSE)
  if (!is.numeric(level) || is.complex(level) || is.object(level) || !is.null(dim(level)) ||
      length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("Supply 0 < level < 1.", call. = FALSE)
  review <- fit$comparison_review
  if (!is.null(contrasts) || !is.null(comparison_review) ||
      (!is.null(facet) && !identical(facet, review$specification$facet)))
    stop("Use the CML comparison fit's saved targets; omit contrasts, comparison_review and any different facet.", call. = FALSE)
  level <- unname(level)
  tab <- mfrmr_comparison_basis(as.data.frame(fit$comparisons))
  eligible <- tab$FormalInferenceEligible %in% TRUE
  tab$Target <- tab$Comparison
  tab$Status <- ifelse(eligible, "available", "unavailable")
  tab$CIEligible <- eligible
  tab$InferenceReview <- tab$EstimateReason
  tab$CoverageQualified <- FALSE
  tab$ConfidenceLevel <- level
  if (level != .95) {
    margin <- stats::qnorm((1 - level) / 2, lower.tail = FALSE) * tab$SE[eligible]
    tab$Lower[eligible] <- tab$Estimate[eligible] - margin
    tab$Upper[eligible] <- tab$Estimate[eligible] + margin
  }
  tab$ModelSE <- tab$SE; tab$ModelLower <- tab$Lower; tab$ModelUpper <- tab$Upper
  notes <- sub("Pointwise 95% Wald", paste0("Pointwise ", format(100 * level, trim = TRUE), "% Wald"),
    fit$notes, fixed = TRUE)
  attr(tab, "notes") <- notes
  structure(list(table = tab, covariance = fit$covariance, model_covariance = fit$covariance,
    parameter_covariance = fit$optimization$covariance,
    model_parameter_covariance = fit$optimization$covariance,
    contrasts = review$weights, comparison_review = review, fit_checks = fit$checks,
    settings = list(facet = paste(review$specification$facet, review$specification$group, sep = ":"),
      method = "model", level = level, adjust = FALSE, adjustment_factor = 1,
      persons = fit$prep$n_person, clusters = NA_integer_, cluster_score_rank = NA_integer_,
      free_parameters = if (is.null(fit$representation)) NA_integer_ else ncol(fit$representation$transform),
      reference = "pointwise normal", quad_points = NA_integer_,
      procedure = "cml_conditional_information", scale = "native_ability",
      target = "facet_by_origin_interaction_contrast",
      parameterization = "identified_conditional_subspace"),
    fit = fit, cautions = notes), class = "mfrm_facet_intervals")
}
