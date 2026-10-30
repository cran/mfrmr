# Fixed-plan bounds for the true expected score, not Wald limits for JML.
# An observed plan cell gets its plan mass divided by its independent response
# count. Unobserved cells retain their entire known score range. Hoeffding's
# inequality is then applied once to the signed event sum for each target.
mfrm_fair_bounded_support <- function(fit) {
  prep <- fit$prep
  if (!identical(prep$rating_min_source, "declared") ||
      !identical(prep$rating_max_source, "declared"))
    stop(paste("Bounded-response intervals require declared rating_min and rating_max from the rubric;",
      "refit with these limits and category_policy = 'preserve'. Observed extrema are not known support."), call. = FALSE)
  low <- prep$rating_min; high <- prep$rating_max
  if (length(low) != 1L || length(high) != 1L || !all(is.finite(c(low, high))) ||
      high <= low || low != floor(low) || high != floor(high))
    stop("Bounded-response intervals require finite integer rubric limits.", call. = FALSE)
  map <- prep$score_map
  if (!is.data.frame(map) || !all(c("OriginalScore", "InternalScore") %in% names(map)) ||
      nrow(map) != high - low + 1 ||
      !identical(as.numeric(map$OriginalScore), as.numeric(seq(low, high))) ||
      !identical(as.numeric(map$InternalScore), as.numeric(seq(low, high))))
    stop(paste("Bounded-response intervals require unchanged full rubric coding;",
      "refit with the declared rating_min, rating_max and category_policy = 'preserve'.",
      "Data-dependent category collapse does not define a fixed score target."), call. = FALSE)
  list(range = c(lower = as.numeric(low), upper = as.numeric(high)), provenance = "declared_rubric",
    coding = "identity")
}

mfrm_fair_bounded_weights <- function(data, plan, facets, persons) {
  keys <- function(d) do.call(paste, c(lapply(facets, function(f)
    match(as.character(d[[f]]), unique(as.character(plan[[f]])))), sep = ":"))
  cell <- match(keys(data), keys(plan))
  pi <- match(as.character(data$Person), persons)
  valid <- !is.na(pi) & !is.na(cell)
  rows <- which(valid)
  # Integer coordinate keys avoid collisions from arbitrary user level labels.
  group <- paste(pi[valid], cell[valid], sep = ":")
  counts <- table(group)
  w <- plan$NormalizedWeight[cell[valid]] / as.numeric(counts[group])
  B <- Matrix::sparseMatrix(i = rows, j = pi[valid], x = w,
    dims = c(nrow(data), length(persons)))
  seen <- lapply(seq_along(persons), function(j) sort(unique(cell[valid & pi == j])))
  observed <- vapply(seen, function(c) sum(plan$NormalizedWeight[c]), 0)
  missing <- vapply(seen, function(c) sum(plan$NormalizedWeight[setdiff(seq_len(nrow(plan)), c)]), 0)
  list(weights = B, support = data.frame(Person = persons,
    ObservedPlanWeight = observed, UnobservedPlanWeight = missing,
    RelevantResponses = tabulate(pi[valid], nbins = length(persons)),
    CompletePlan = lengths(seen) == nrow(plan), stringsAsFactors = FALSE))
}

mfrm_fair_bounded_result <- function(fit, plan, parm, L, estimate, level, has_contrasts) {
  support <- mfrm_fair_bounded_support(fit)
  d <- fit$prep$data
  if (any(d$Weight != 1))
    stop("Bounded-response intervals require unit weights; replicated rows must represent independent responses.", call. = FALSE)
  low <- unname(support$range["lower"]); high <- unname(support$range["upper"])
  if (any(!is.finite(d$Score)) || any(d$Score < low | d$Score > high))
    stop("Bounded-response intervals require finite scores within their declared range.", call. = FALSE)
  design <- mfrm_fair_bounded_weights(d, plan, fit$config$facet_names, parm)
  U <- design$weights %*% Matrix::t(Matrix::Matrix(L, sparse = TRUE))
  observed <- as.numeric(Matrix::crossprod(U, d$Score))
  unobserved <- sweep(L, 2L, design$support$UnobservedPlanWeight, "*")
  missing_low <- rowSums(pmax(unobserved, 0) * low + pmin(unobserved, 0) * high)
  missing_high <- rowSums(pmax(unobserved, 0) * high + pmin(unobserved, 0) * low)
  support_low <- rowSums(pmax(L, 0) * low + pmin(L, 0) * high)
  support_high <- rowSums(pmax(L, 0) * high + pmin(L, 0) * low)
  # Normalize before squaring: small/large but finite contrast multipliers must
  # not turn a nonzero radius into zero or infinity through squared weights.
  scale <- as.numeric(Matrix::colSums(abs(U)))
  unit <- methods::as(U, "dgCMatrix")
  if (length(unit@x)) unit@x <- unit@x / rep.int(scale, diff(unit@p))
  norm <- scale * sqrt(as.numeric(Matrix::colSums(unit^2)))
  # log(2) - log1p(-level) avoids loss from forming 2/(1-level).
  radius <- norm * (high - low) * sqrt((log(2) - log1p(-level)) / 2)
  lower <- pmax(support_low, observed + missing_low - radius)
  upper <- pmin(support_high, observed + missing_high + radius)
  if (any(!is.finite(c(observed, radius, lower, upper, support_low, support_high))) ||
      any(lower > upper)) stop("Bounded-response target arithmetic is nonfinite or inconsistent.", call. = FALSE)
  complete <- vapply(seq_len(nrow(L)), function(i)
    all(design$support$CompletePlan[L[i, ] != 0]), TRUE)
  informative <- (upper - lower) < (support_high - support_low) * (1 - 1e-12)
  design_estimate <- ifelse(complete, observed, NA_real_)
  review <- paste("Pointwise Hoeffding bound for a fixed plan and contrast under independent bounded responses",
    "with stable within-cell means. Unobserved cells retain their full score range.",
    "Fitted points are unchanged; these bounds do not correct their bias or establish fairness.")
  out <- data.frame(Target = rownames(L), Estimate = unname(estimate), SE = NA_real_,
    CI_Lower = unname(lower), CI_Upper = unname(upper), Level = level,
    Available = TRUE, Informative = informative,
    Status = ifelse(informative, "conservative_bound", "uninformative_bound"),
    Planned = NA_integer_, Resolved = NA_integer_, Unresolved = NA_integer_,
    ExpectedTailDraws = NA_real_, CIEligible = TRUE,
    Uncertainty = "independent_bounded_responses", InferenceReview = review,
    Reason = ifelse(!informative, "full_target_support", ifelse(complete,
      "complete_plan_support", "unobserved_plan_cells_bounded")),
    DesignEstimate = unname(design_estimate), ObservedContribution = unname(observed),
    SamplingRadius = unname(radius), UnobservedLower = unname(missing_low),
    UnobservedUpper = unname(missing_high), SupportLower = unname(support_low),
    SupportUpper = unname(support_high), row.names = NULL)
  attr(out, "reference_plan") <- plan
  attr(out, "support") <- design$support
  attr(out, "person_event_weights") <- design$weights
  attr(out, "source_identity") <- list(scoring = mfrm_fair_fit_identity(fit))
  attr(out, "settings") <- list(parm = parm, contrasts = L, level = level,
    estimator = if (mfrm_is_jeffreys_fit(fit)) "jeffreys_canonical_joint" else "ordinary_jml",
    method = "hoeffding_bounded_responses", source_events = "scoring_only",
    score_support = support,
    anchor_uncertainty = "anchors_not_used_by_bound",
    coverage_basis = "finite_sample_inequality_under_declared_assumptions",
    multiple_targets = "pointwise_not_simultaneous", gradient_tol = NA_real_,
    target = if (has_contrasts) "common_plan_expected_average_contrast" else
      "common_plan_expected_average")
  class(out) <- c("mfrm_fair_intervals", "data.frame")
  out
}
