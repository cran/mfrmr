# Linearize a pooled available-row mean over independent simulation trials.
recovery_trial_metric <- function(x, rep, trials, root_mean_square = FALSE) {
  ok <- is.finite(x) & !is.na(rep)
  x <- x[ok]
  rep <- rep[ok]
  if (!length(x)) return(c(value = NA_real_, mcse = NA_real_, reps = 0))
  scale <- max(abs(x))
  y <- if (scale > 0) x / scale else x
  if (root_mean_square) y <- y^2
  index <- match(rep, trials)
  if (anyNA(index)) stop("Recovery rows refer to unrecorded trials.", call. = FALSE)
  n <- tabulate(index, nbins = length(trials))
  sums <- numeric(length(trials))
  totals <- rowsum(y, index, reorder = FALSE)
  sums[as.integer(rownames(totals))] <- totals[, 1L]
  mu <- sum(sums) / sum(n)
  available <- sum(n > 0)
  mcse <- NA_real_
  if (available > 1L) {
    influence <- sums - mu * n
    mcse <- sqrt(sum(influence^2) * length(trials) / (length(trials) - 1)) / sum(n)
    if (root_mean_square) mcse <- if (mu > 0) mcse / (2 * sqrt(mu)) else NA_real_
  }
  c(value = scale * if (root_mean_square) sqrt(mu) else mu,
    mcse = scale * mcse, reps = available)
}

recovery_identified_truth <- function(truth, model, step_facet, population = FALSE, levels = NULL) {
  steps <- recovery_normalize_step_table(truth$step_table)
  step_mean <- tapply(steps$Estimate, steps$StepFacet, mean)
  common_shift <- 0
  if (model %in% c("PCM", "GPCM") && step_facet %in% names(truth$facets)) {
    facet <- truth$facets[[step_facet]]
    shift <- if ("Common" %in% names(step_mean)) rep(step_mean[["Common"]], length(facet)) else
      as.numeric(step_mean[names(facet)])
    truth$facets[[step_facet]] <- facet + shift
  } else if (identical(model, "RSM") && length(step_mean) == 1L) {
    common_shift <- unname(step_mean[1])
  }
  if (population && "(Intercept)" %in% names(truth$population$coefficients)) {
    facet_means <- vapply(names(truth$facets), function(facet) {
      x <- truth$facets[[facet]]
      if (!is.null(levels[[facet]])) x <- x[levels[[facet]]]
      mean(x)
    }, numeric(1))
    truth$population$coefficients["(Intercept)"] <- truth$population$coefficients["(Intercept)"] -
      sum(facet_means) - common_shift
  }
  truth
}

recovery_target_counts <- function(truth, rep, model, step_facet, include_person) {
  rows <- list()
  add <- function(type, facet, scale, n) {
    rows[[length(rows) + 1L]] <<- tibble::tibble(rep = rep, ParameterType = type,
      Facet = facet, ComparisonScale = scale, ExpectedRows = as.integer(n))
  }
  if (include_person && length(truth$person)) add("person", "Person", "logit", length(truth$person))
  for (facet in names(truth$facets)) add("facet", facet, "logit", length(truth$facets[[facet]]))
  steps <- recovery_normalize_step_table(truth$step_table)
  if (nrow(steps)) {
    if (model %in% c("PCM", "GPCM") && all(steps$StepFacet == "Common")) {
      for (level in names(truth$facets[[step_facet]])) add("step", level, "logit", nrow(steps))
    } else for (level in unique(steps$StepFacet)) add("step", level, "logit", sum(steps$StepFacet == level))
  }
  if (is.data.frame(truth$slope_table) && nrow(truth$slope_table))
    add("slope", "SlopeFacet", "log_slope", nrow(truth$slope_table))
  if (length(truth$population$coefficients))
    add("population", "population", "coefficient", length(truth$population$coefficients))
  if (length(truth$population$sigma2)) add("population", "population", "variance", 1L)
  dplyr::bind_rows(rows)
}

recovery_interval_rows <- function(rows, fit, truth) {
  rows$RawSE <- rows$SE
  rows$SE <- rows$CoverageEstimate <- rows$CoverageTruth <- rows$Lower95 <- rows$Upper95 <- NA_real_
  rows$CIEligible <- FALSE
  rows$CIReason <- "interval_not_implemented_for_parameter_group"
  rows$CoverageBasis <- "centered_facet_model_interval"
  rows$CIReason[rows$ParameterType == "person"] <- "person_score_uncertainty_is_not_a_centered_structural_interval"
  facets <- unique(rows$Facet[rows$ParameterType == "facet"])
  for (facet in facets) {
    index <- which(rows$ParameterType == "facet" & rows$Facet == facet)
    if (!as.character(fit$config$method %||% "unknown")[1] %in% c("MML", "CML")) {
      rows$CIReason[index] <- "jml_point_estimation_not_structural_truth_inference"
      next
    }
    levels <- as.character(fit$facets$others$Level[fit$facets$others$Facet == facet])
    if (!setequal(rows$Level[index], levels) ||
        !setequal(rows$Level[index], names(truth$facets[[facet]])) ||
        any(!is.finite(rows$Estimate[index])) ||
        any(!is.finite(rows$Truth[index]))) {
      rows$CIReason[index] <- "incomplete_alignment_group"
      next
    }
    m <- length(index)
    contrast <- diag(m) - 1 / m
    dimnames(contrast) <- list(rows$Level[index], rows$Level[index])
    ci <- tryCatch(mfrm_facet_intervals(fit, facet, contrasts = contrast), error = function(e) e)
    rows$CoverageTruth[index] <- rows$Truth[index] - mean(rows$Truth[index])
    if (inherits(ci, "error")) {
      rows$CIReason[index] <- paste0("interval_unavailable: ", conditionMessage(ci))
      next
    }
    tab <- ci$table[match(rows$Level[index], ci$table$Target), , drop = FALSE]
    centered <- rows$Estimate[index] - mean(rows$Estimate[index])
    matched <- is.finite(tab$Estimate) & abs(tab$Estimate - centered) <=
      sqrt(.Machine$double.eps) * pmax(1, abs(centered))
    se_available <- matched & is.finite(tab$SE) & tab$SE > 0
    eligible <- tab$CIEligible %in% TRUE & se_available &
      is.finite(tab$Lower) & is.finite(tab$Upper)
    rows$CoverageEstimate[index] <- tab$Estimate
    rows$SE[index] <- ifelse(se_available, tab$SE, NA_real_)
    rows$Lower95[index] <- ifelse(eligible, tab$Lower, NA_real_)
    rows$Upper95[index] <- ifelse(eligible, tab$Upper, NA_real_)
    rows$CIEligible[index] <- eligible
    rows$CIReason[index] <- ifelse(eligible, "available", paste(tab$Status, tab$InferenceReview, sep = ": "))
    rows$CIReason[index[!matched]] <- "interval_point_does_not_match_recovery_target"
  }
  recovery_review_rows(rows)
}

recovery_review_rows <- function(rows) {
  rows <- tibble::as_tibble(rows %||% tibble::tibble())
  if (!nrow(rows)) return(data.frame(rep = integer(), ParameterType = character(),
    Facet = character(), ComparisonScale = character(), Truth = numeric(),
    EstimateAligned = numeric(), ErrorAligned = numeric(), ErrorRaw = numeric(),
    SE = numeric(), Covered95 = logical(), RecoveryComparable = logical(),
    RecoveryBasis = character(), CoverageBasis = character(), CIReason = character()))
  required <- c("CIEligible", "CoverageBasis", "CoverageTruth", "Lower95", "Upper95", "CIReason")
  if (!all(required %in% names(rows))) {
    rows$RawSE <- rows[["RawSE"]] %||% rows$SE
    rows$SE <- NA_real_
    rows$CIEligible <- FALSE
    rows$CIReason <- "legacy_interval_basis_unverified"
    rows$CoverageBasis <- "unverified_legacy"
    rows$CoverageTruth <- rows$Lower95 <- rows$Upper95 <- NA_real_
  }
  eligible <- rows$CIEligible %in% TRUE & is.finite(rows$CoverageTruth) &
    is.finite(rows$Lower95) & is.finite(rows$Upper95)
  rows$Covered95 <- ifelse(eligible, rows$CoverageTruth >= rows$Lower95 &
    rows$CoverageTruth <= rows$Upper95, NA)
  rows
}
