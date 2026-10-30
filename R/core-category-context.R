# Explicit rating contexts use the same free-coordinate predictor as curve inference.
category_context_curve_spec <- function(fit, newdata) {
  stop_if_jml_adjustment(fit, "category_curves_report")
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      !identical(fit$config$model, "GPCM") || !identical(fit$config$method, "MML")) {
    stop("Explicit `newdata` requires a native GPCM MML fit.", call. = FALSE)
  }
  facets <- fit$config$facet_names
  if (!is.data.frame(newdata) || !nrow(newdata) || anyDuplicated(names(newdata)) ||
      !all(facets %in% names(newdata)) || "Theta" %in% names(newdata)) {
    stop("`newdata` must contain every non-Person facet and no Theta column.", call. = FALSE)
  }
  newdata <- as.data.frame(newdata[facets], stringsAsFactors = FALSE)
  query <- newdata
  query$Theta <- 0
  evaluator <- mfrm_gpcm_response_evaluator(fit$prep, fit$config, query)
  value <- evaluator$evaluate(fit$opt$par)
  labels <- make.unique(vapply(seq_len(nrow(newdata)), function(i) {
    paste(paste(facets, vapply(newdata[i, , drop = FALSE], as.character, character(1)),
      sep = " = "), collapse = "; ")
  }, character(1)))
  profiles <- data.frame(CurveGroup = labels, InputRow = seq_len(nrow(newdata)),
    StepOwner = fit$config$step_facet,
    StepLevel = as.character(newdata[[fit$config$step_facet]]),
    SlopeOwner = paste(fit$config$slope_facet, collapse = ", "),
    EffectiveSlope = value$effective_slope, PredictorOffset = value$predictor_offset,
    ObservedContext = evaluator$contexts$ObservedContext)
  groups <- setNames(lapply(seq_len(nrow(newdata)), function(i) {
    list(name = labels[i], step_cum = value$step_cum[i, ],
      slope = value$effective_slope[i], predictor_offset = value$predictor_offset[i])
  }), labels)
  list(model = "GPCM", categories = seq_len(fit$config$n_cat) - 1L + fit$prep$rating_min,
    groups = groups, profiles = profiles, newdata = newdata,
    curve_basis = "fitted_rating_context", predictor_offset = 0,
    basis_description = paste("Full fitted predictor at the supplied rating contexts, including facet locations",
      "and interactions; two-family discriminations are multiplied. Theta is on the native fitted scale."))
}

# Compare unnormalized lower/upper masses. Pair exponential terms before
# subtraction so tiny positive slopes do not round every probability to 1/K.
category_cumulative_mass_balance <- function(theta, group, boundary) {
  k <- seq_along(group$step_cum) - 1L
  eta <- theta + group$predictor_offset
  log_mass <- group$slope * (k * eta - group$step_cum)
  if (any(!is.finite(log_mass))) return(NA_real_)
  shift <- max(log_mass)
  n <- length(k)
  paired <- seq_len(min(boundary, n - boundary))
  lower <- paired
  upper <- boundary + paired
  delta <- group$slope * ((k[lower] - k[upper]) * eta -
    (group$step_cum[lower] - group$step_cum[upper]))
  pairs <- exp(pmax(log_mass[lower], log_mass[upper]) - shift) *
    sign(delta) * (-expm1(-abs(delta)))
  remaining_lower <- setdiff(seq_len(boundary), lower)
  remaining_upper <- setdiff(seq.int(boundary + 1L, n), upper)
  sum(pairs) + sum(exp(log_mass[remaining_lower] - shift)) -
    sum(exp(log_mass[remaining_upper] - shift))
}

# Positive GPCM slopes make each cumulative boundary unique. Never substitute
# a search endpoint when the target probability is outside the requested range.
category_context_boundaries <- function(boundaries, spec, theta_range, fit) {
  out <- boundaries
  out$AndrichThreshold <- out$AdjacentThreshold <- out$PredictorOffset <-
    out$EffectiveSlope <- out$RootResidual <- out$LowerCumulativeProbability <-
    out$UpperCumulativeProbability <- NA_real_
  out$BoundaryDetail <- ""
  out$RootMethod <- "paired_mass_uniroot"
  out$RootTolerance <- 1e-10
  out$ProbabilityTolerance <- 1e-8
  for (i in seq_len(nrow(out))) {
    grp <- spec$groups[[out$CurveGroup[i]]]
    j <- out$BoundaryOrder[i]
    binary <- length(grp$step_cum) == 2L
    out$RootMethod[i] <- if (binary) "analytic_binary" else "paired_mass_uniroot"
    cumulative <- function(theta) {
      p <- category_prob_gpcm(theta + grp$predictor_offset,
        matrix(grp$step_cum, nrow = 1L), 1L, grp$slope, 1L)
      sum(p[1, seq_len(j)])
    }
    ends <- vapply(theta_range, cumulative, numeric(1))
    out$LowerCumulativeProbability[i] <- ends[1]
    out$UpperCumulativeProbability[i] <- ends[2]
    out$PredictorOffset[i] <- grp$predictor_offset
    out$EffectiveSlope[i] <- grp$slope
    out$AndrichThreshold[i] <- diff(grp$step_cum)[j]
    out$AdjacentThreshold[i] <- out$AndrichThreshold[i] - grp$predictor_offset
    out$ThurstonianThreshold[i] <- NA_real_
    out$InThetaRange[i] <- FALSE
    out$CrossingCount[i] <- 0L
    out$BoundaryStatus[i] <- "outside_theta_range"
    # For a binary response the common positive slope cancels analytically.
    # Do not let a subnormal slope underflow the local-precision check.
    balance <- if (binary) {
      function(theta) diff(grp$step_cum) - (theta + grp$predictor_offset)
    } else {
      function(theta) category_cumulative_mass_balance(theta, grp, j)
    }
    bracket <- vapply(theta_range, balance, numeric(1))
    if (any(!is.finite(bracket))) {
      out$BoundaryStatus[i] <- "numerical_failure"
      out$BoundaryDetail[i] <- "Cumulative mass comparison is nonfinite at the search endpoints."
      next
    }
    if (bracket[1] < 0 || bracket[2] > 0) next
    if (binary) {
      root <- out$AdjacentThreshold[i]
    } else {
      root <- tryCatch(stats::uniroot(balance, interval = theta_range,
        f.lower = bracket[1], f.upper = bracket[2],
        tol = out$RootTolerance[i], check.conv = TRUE)$root,
        error = function(e) e)
    }
    if (inherits(root, "error")) {
      out$BoundaryStatus[i] <- "numerical_failure"
      out$BoundaryDetail[i] <- conditionMessage(root)
      next
    }
    neighborhood <- root + c(-1, 1) * out$RootTolerance[i]
    signs <- vapply(neighborhood, balance, numeric(1))
    if (!is.finite(root) || root < theta_range[1] || root > theta_range[2] ||
        neighborhood[1] >= root || neighborhood[2] <= root ||
        any(!is.finite(signs)) || signs[1] <= 0 || signs[2] >= 0) {
      out$BoundaryStatus[i] <- "numerical_failure"
      out$BoundaryDetail[i] <- "Root neighborhood did not establish the requested ability tolerance."
      next
    }
    out$RootResidual[i] <- cumulative(root) - .5
    if (!is.finite(out$RootResidual[i]) || abs(out$RootResidual[i]) > out$ProbabilityTolerance[i]) {
      out$BoundaryStatus[i] <- "numerical_failure"
      out$BoundaryDetail[i] <- "Root did not meet the cumulative-probability tolerance."
      next
    }
    out$ThurstonianThreshold[i] <- root
    out$InThetaRange[i] <- TRUE
    out$CrossingCount[i] <- 1L
    out$BoundaryStatus[i] <- "in_range"
  }
  out$LowerCategoryLabel <- category_original_labels(as.numeric(out$LowerOrEqualCategory), fit)
  out$UpperCategoryLabel <- category_original_labels(as.numeric(out$AboveCategory), fit)
  out$BoundaryLabel <- paste0("P(X <= ", out$LowerCategoryLabel, ") = 0.5")
  out$CurveBasis <- spec$curve_basis
  out$SE <- out$CI_Lower <- out$CI_Upper <- NA_real_
  out$CIEligible <- FALSE
  out$UncertaintyStatus <- if (mfrm_has_product_slopes(fit)) "unavailable_two_family" else "not_computed"
  out
}
