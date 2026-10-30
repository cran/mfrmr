#' Evaluate classification accuracy, consistency and conditional precision
#'
#' Evaluate a fixed future rating plan under an ordinary RSM/PCM calibration.
#' Lee-style recursion gives the exact finite summed-score distribution;
#' the Rudner comparison uses a normal approximation to an ML ability estimate.
#'
#' Both methods are returned. For an operational decision based on the declared
#' category-score sum, report the `lee_sum` row. For a decision based on an
#' ML ability estimate, `rudner_normal` is an approximate comparison, subject
#' to the limitations below. It does not evaluate an EAP/MAP decision rule.
#' The rows need not describe the same decision; compare their accuracy as
#' competing calculations only after matching the actual operational rules.
#' `Status = "available"` means the calculation succeeded, not that the normal
#' approximation or the fitted model is accurate for the intended population.
#'
#' @param fit A native ordinary RSM/PCM `fit_mfrm()` result (MML or JML),
#'   without interactions, random effects or a JML correction.
#' @param rating_plan Data frame with one row per future rating and exactly
#'   the calibrated facet identifier columns. All levels must be known.
#'   Duplicate complete facet cells, missing identifiers, score columns and
#'   observation weights are not supported. Each row contributes one rating.
#' @param cut_scores Finite, strictly increasing true-ability category cuts
#'   on the calibration scale. Optional names label the cuts.
#' @param score_cuts Nondecreasing integer cuts on the sum of **model category
#'   indices**, each rating scored `0, ..., K-1`. One cut per `cut_scores`.
#'   Values from zero through the maximum possible sum plus one are allowed.
#'   Repeated cuts retain empty categories; zero and maximum-plus-one retain
#'   unreachable edge categories. These are explicit operational decisions,
#'   not automatically rounded expected scores or sums of original labels.
#' @param theta Finite numeric vector of true-ability points to evaluate.
#'   Supply a grid around the cuts to inspect conditional precision there.
#'   Fitted EAPs are not automatically substituted for true abilities.
#' @param theta_weights Optional nonnegative population point/quadrature
#'   weights corresponding to `theta`, with at least one positive entry.
#'   Weights are normalized and retained. Density values alone are generally
#'   insufficient on an unequally spaced grid. No default population is used.
#' @param operational_cuts Strictly increasing ability-estimate decision cuts
#'   for the Rudner approximation, one per true cut. Defaults to `cut_scores`.
#' @param facets Optional character vector mapping calibrated facet names to
#'   `rating_plan` columns, as in [predict_mfrm_units()].
#' @param readiness_policy `"error"` requires a ready source calibration.
#'   `"review"` permits a nonready but valid source and retains its assessment
#'   in `source`; it does not establish statistical adequacy. Invalid retained
#'   parameters cannot be overridden.
#'
#' @details
#' Categories are numbered from 1. Values equal to a cut enter the higher
#' category: lower endpoints are inclusive, upper endpoints exclusive. The
#' top category includes the maximum possible sum. True categories use
#' `cut_scores`; operational categories use `score_cuts` for `lee_sum` and
#' `operational_cuts` for `rudner_normal`. These can describe different rules.
#' Expected sums at ability cuts are reported for inspection, never rounded
#' into a decision rule. In an ordinary fixed-plan Rasch model the ML score
#' equation is `S = E(S | theta)`. A matched sum cut therefore requires the
#' actual ML decision mapping, including ties and extreme-score policy.
#' EAP/MAP decisions do not automatically share that mapping.
#'
#' With conditional operational category probabilities `p[h]` and true
#' category `H(theta)`, accuracy is `p[H(theta)]`; consistency is `sum(p^2)`.
#' Consistency concerns two conditionally independent repetitions with the
#' **same** calibrated plan. Overclassification and underclassification sum
#' probabilities above and below the true category. Population versions are
#' weighted averages over the supplied true-ability distribution. The joint
#' matrix contains population joint probabilities, not row-normalized error
#' rates. Accuracy need not exceed consistency in a multicategory decision.
#' With two categories, marginal overclassification is the joint probability
#' of a true fail and an observed pass; underclassification is the joint
#' probability of a true pass and an observed fail. Divide each by its true
#' category's population probability for the error rate conditional on that
#' category. A zero denominator makes that conditional rate undefined.
#' Conditioning instead on the observed pass/fail category uses its column
#' sum in the joint matrix and answers how often that reported decision is wrong.
#'
#' Information sums conditional category-score variances over the plan.
#' `SEM = 1/sqrt(Information)` is a conditional **ability-scale** SEM with
#' calibration fixed, not a posterior SD, sum-score SEM or an SE propagating
#' facet estimation uncertainty. `rudner_normal` assumes an ML estimator's
#' sampling distribution is normal with mean `theta` and this SEM. Short
#' plans, extremes, estimator bias and model misfit can invalidate that
#' approximation. Zero numerical information gives infinite SEM and makes
#' Rudner results unavailable; Lee probabilities can remain available.
#' Arithmetic overflow or unresolved probability mass causes an error.
#' Narrow normal-probability bands use direct integration when subtraction
#' of cumulative probabilities would lose precision.
#' Finite precision does not guarantee relative accuracy in vanishing tails.
#'
#' Local independence is assumed even when rows share a task or rater.
#' Calibration, cuts and the plan are held fixed. There is no uncertainty
#' interval, quality threshold, automatic cut selection, random-rater/task
#' generalization, retake policy, or propagation of calibration/population
#' estimation error. Model-based accuracy is not observed reference agreement.
#' This is also distinct from posterior ability-band probabilities in
#' [predict_mfrm_units()]. Passing calibration checks does not validate
#' the future plan or the normal approximation.
#'
#' @return A plain list. `table` contains conditional indices, expected sum,
#'   information, SEM and availability for each point/method. `SourceStatus`
#'   in `table`, `probabilities`, `scores`, `cuts`, `marginal` and
#'   `classification_matrix` preserves the separate calibration assessment,
#'   including in CSV exports. `Status`, where present, describes numerical
#'   availability; an available row with `SourceStatus = "review_only"` still
#'   needs source review. Inspect `source$reason_codes` for the reasons.
#'   `probabilities`
#'   contains operational category probabilities; `scores` contains the exact
#'   summed-score probability mass. `cuts` reports both cut scales and SEM at
#'   true and operational ability cuts. `marginal` and `classification_matrix`
#'   are empty unless population weights are supplied. Any unavailable point
#'   with positive weight makes its method's marginal results unavailable;
#'   such points are never dropped or renormalized away. `population` retains
#'   input points and raw/normalized weights. `rating_plan`, `calibration`
#'   (offsets, cumulative thresholds and category map), `source`, `settings`
#'   and `notes` retain the calculation's conditions. A source review status
#'   applies to the entire result, including numerically available rows.
#'   Save with [base::saveRDS()] or export tables with [utils::write.csv()].
#'   A table CSV alone does not retain the plan, cut definitions, population
#'   distribution and notes; keep the full result for reproducible reuse.
#'   Reading saved output does not refit or recalculate anything.
#'
#' @references
#' Rudner, L. M. (2001). Computing the expected proportions of misclassified
#' examinees. *Practical Assessment, Research, and Evaluation, 7*, Article 14.
#' \doi{10.7275/an9m-2035}.
#'
#' Rudner, L. M. (2005). Expected classification accuracy.
#' *Practical Assessment, Research, and Evaluation, 10*, Article 13.
#' \doi{10.7275/56a5-6b14}.
#'
#' Lee, W.-C. (2010). Classification consistency and accuracy for complex
#' assessments using item response theory. *Journal of Educational
#' Measurement, 47*(1), 1-17. \doi{10.1111/j.1745-3984.2009.00096.x}.
#'
#' Wyse, A. E., & Hao, S. (2012). An evaluation of item response theory
#' classification accuracy and consistency indices. *Applied Psychological
#' Measurement, 36*(7), 602-624. \doi{10.1177/0146621612451522}.
#'
#' Lathrop, Q. N., & Cheng, Y. (2013). Two approaches to estimation of
#' classification accuracy rate under item response theory. *Applied
#' Psychological Measurement, 37*(3), 226-241. \doi{10.1177/0146621612471888}.
#' @seealso [predict_mfrm_units()], [score_comparison_table()]
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("classification_table", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_operational")
#' fit <- fit_mfrm(ratings, person = "Person", facets = c("Rater", "Criterion"),
#'   score = "Score", method = "MML", model = "RSM",
#'   rating_min = 1, rating_max = 4, category_policy = "preserve")
#' # The planned roster includes the intended ratings even if scores are missing.
#' data("mfrmr_example_operational_design", package = "mfrmr")
#' rating_plan <- subset(mfrmr_example_operational_design, Person == "P001",
#'   select = c("Rater", "Criterion"))
#' rating_plan # Six ratings: two raters, each scoring three criteria.
#' # Illustrative rules: ability at least 0 logits; original-score sum at least 18.
#' # score_cuts uses 0:3 category indices, not the original 1:4 score labels:
#' # subtract 1 for each of the six ratings, so 18 - 6 = 12.
#' # These are different decision rules, not automatically equivalent cutoffs.
#' result <- classification_table(fit, rating_plan,
#'   cut_scores = 0, score_cuts = 12, theta = c(-0.5, 0, 0.5))
#' result$table
#' result$cuts
#' result$settings$calibration # These probabilities hold calibration fixed.
#' # No population average is inferred from this three-point inspection grid.
#' }
#' @export
classification_table <- function(fit, rating_plan, cut_scores, score_cuts, theta,
                                 theta_weights = NULL, operational_cuts = cut_scores,
                                 facets = NULL, readiness_policy = c("error", "review")) {
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      mfrm_extended_fit(fit) || mfrm_has_jml_adjustment(fit) ||
      !fit$config$model %in% c("RSM", "PCM") ||
      !fit$config$method %in% c("MML", "JML") ||
      length(fit$config$interaction_specs))
    stop("`fit` must be a native ordinary RSM/PCM MML/JML calibration without interactions or random effects.", call. = FALSE)
  readiness_policy <- match.arg(readiness_policy)
  cut_scores <- mfrm_validate_cut_scores(cut_scores)
  if (is.null(cut_scores)) stop("Supply at least one true ability cut.", call. = FALSE)
  operational_cuts <- mfrm_validate_cut_scores(operational_cuts)
  if (length(operational_cuts) != length(cut_scores))
    stop("Supply one `operational_cuts` value per true ability cut.", call. = FALSE)
  finite_vector <- function(x) is.numeric(x) && !is.complex(x) && !is.object(x) &&
    is.null(dim(x)) && length(x) > 0L && all(is.finite(x))
  if (!finite_vector(theta)) stop("`theta` must be a nonempty finite numeric vector.", call. = FALSE)
  if (!is.null(theta_weights) && (!finite_vector(theta_weights) ||
      length(theta_weights) != length(theta) || any(theta_weights < 0) || !any(theta_weights > 0)))
    stop("`theta_weights` must match theta and be finite, nonnegative, with positive total mass.", call. = FALSE)
  if (!is.null(facets) && (anyNA(facets) || anyNA(names(facets))))
    stop("`facets` must not contain missing names.", call. = FALSE)
  mapping <- resolve_prediction_facets(fit, facets)
  if (!is.data.frame(rating_plan) || !nrow(rating_plan) || anyDuplicated(names(rating_plan)) ||
      anyDuplicated(unname(mapping)) || !setequal(names(rating_plan), unname(mapping)))
    stop("`rating_plan` must contain exactly the mapped facet columns and at least one rating.", call. = FALSE)
  plan <- as.data.frame(rating_plan[unname(mapping)], stringsAsFactors = FALSE)
  names(plan) <- names(mapping)
  for (f in names(plan)) {
    ids <- plan[[f]]
    if (!(is.character(ids) || is.factor(ids) || is.numeric(ids)) ||
        !is.null(dim(ids)) || anyNA(ids) || any(!nzchar(trimws(as.character(ids)))))
      stop("Rating-plan identifiers must be nonmissing, nonblank vectors.", call. = FALSE)
    plan[[f]] <- as.character(ids)
    if (any(!plan[[f]] %in% fit$config$facet_levels[[f]]))
      stop("Unknown rating-plan level in facet `", f, "`.", call. = FALSE)
  }
  if (anyDuplicated(plan)) stop("Duplicate complete facet cells in `rating_plan` are not supported.", call. = FALSE)
  max_sum <- nrow(plan) * (fit$config$n_cat - 1)
  if (!finite_vector(score_cuts) || length(score_cuts) != length(cut_scores) ||
      any(score_cuts != floor(score_cuts)) || any(diff(score_cuts) < 0) ||
      any(score_cuts < 0 | score_cuts > max_sum + 1))
    stop("`score_cuts` must be nondecreasing integers from 0 through maximum sum + 1, one per true cut.", call. = FALSE)
  source <- prediction_source_scoring_readiness(fit)
  if (identical(source$status, "unavailable"))
    stop("The retained calibration is invalid; review cannot enable classification.", call. = FALSE)
  if (!isTRUE(source$ready) && readiness_policy == "error")
    stop("The source calibration is not ready: ", paste(source$reason_codes, collapse = ", "),
         ". Use readiness_policy = 'review' only for explicit review.", call. = FALSE)
  cfg <- fit$config
  params <- expand_params(fit$opt$par, build_param_sizes(cfg), cfg)
  idx <- list(score_k = rep(0L, nrow(plan)), facets = lapply(names(plan), function(f)
    match(plan[[f]], cfg$facet_levels[[f]])), interactions = list())
  names(idx$facets) <- names(plan)
  offsets <- compute_base_eta(idx, params, cfg)
  steps <- if (cfg$model == "RSM")
    matrix(params$steps, nrow(plan), cfg$n_cat - 1L, byrow = TRUE) else
    params$steps_mat[idx$facets[[cfg$step_facet]], , drop = FALSE]
  cumulative <- t(vapply(seq_len(nrow(plan)), function(j) c(0, cumsum(steps[j, ])),
                         numeric(cfg$n_cat)))
  n <- length(theta); categories <- seq_len(length(cut_scores) + 1L)
  true <- findInterval(theta, cut_scores) + 1L
  methods <- c("lee_sum", "rudner_normal")
  tables <- probabilities <- scores <- vector("list", n)
  for (i in seq_len(n)) {
    point <- mfrm_classification_point(theta[i], offsets, cumulative, score_cuts, operational_cuts)
    p <- rbind(point$lee, point$rudner)
    status <- c("available", if (all(is.finite(point$rudner))) "available" else "zero_information")
    tables[[i]] <- data.frame(Point = i, Theta = theta[i], TrueCategory = true[i], Method = methods,
      Accuracy = p[, true[i]], Consistency = rowSums(p^2),
      Overclassification = rowSums(p[, categories > true[i], drop = FALSE]),
      Underclassification = rowSums(p[, categories < true[i], drop = FALSE]),
      ExpectedSum = point$expected, Information = point$information, SEM = point$sem, Status = status,
      row.names = NULL)
    # An unavailable method cannot report a misleading zero for an empty error tail.
    tables[[i]][status != "available", c("Accuracy", "Consistency", "Overclassification", "Underclassification")] <- NA_real_
    probabilities[[i]] <- data.frame(Point = i, Theta = theta[i], Method = rep(methods, each = length(categories)),
      Category = rep(categories, 2), Probability = as.vector(t(p)), Status = rep(status, each = length(categories)))
    scores[[i]] <- data.frame(Point = i, Theta = theta[i], SumScore = seq_along(point$pmf) - 1L,
                              Probability = point$pmf)
  }
  table <- do.call(rbind, tables)
  probabilities <- do.call(rbind, probabilities)
  cut_moments <- function(x) t(vapply(x, function(t) {
    z <- mfrm_classification_response(t, offsets, cumulative)
    c(z$expected, z$information, z$sem)
  }, numeric(3)))
  at_true <- cut_moments(cut_scores); at_operational <- cut_moments(operational_cuts)
  cuts <- data.frame(Cut = names(cut_scores) %||% paste0("Cut", seq_along(cut_scores)),
    TrueCut = unname(cut_scores), OperationalAbilityCut = unname(operational_cuts),
    ScoreCut = unname(score_cuts), ExpectedSumAtTrueCut = at_true[, 1],
    InformationAtTrueCut = at_true[, 2], SEMAtTrueCut = at_true[, 3],
    ExpectedSumAtOperationalCut = at_operational[, 1],
    InformationAtOperationalCut = at_operational[, 2], SEMAtOperationalCut = at_operational[, 3])
  population <- data.frame(Point = seq_len(n), Theta = unname(theta),
    Weight = if (is.null(theta_weights)) NA_real_ else unname(theta_weights), NormalizedWeight = NA_real_)
  marginal <- table[FALSE, c("Method", "Accuracy", "Consistency", "Overclassification", "Underclassification", "Status")]
  joint <- data.frame(Method = character(), TrueCategory = integer(), ObservedCategory = integer(),
                       Probability = numeric(), Status = character())
  if (!is.null(theta_weights)) {
    w <- theta_weights / max(theta_weights); w <- w / sum(w)
    population$NormalizedWeight <- unname(w)
    marginal <- do.call(rbind, lapply(methods, function(method) {
      d <- table[table$Method == method & theta_weights[table$Point] > 0, , drop = FALSE]
      ok <- all(d$Status == "available")
      out <- data.frame(Method = method)
      for (col in c("Accuracy", "Consistency", "Overclassification", "Underclassification"))
        out[[col]] <- if (ok) sum(d[[col]] * w[d$Point]) else NA_real_
      out$Status <- if (ok) "available" else "positive_weight_point_unavailable"
      out
    }))
    joint <- do.call(rbind, lapply(methods, function(method) {
      d <- probabilities[probabilities$Method == method & theta_weights[probabilities$Point] > 0, ]
      out <- expand.grid(TrueCategory = categories, ObservedCategory = categories)
      ok <- all(d$Status == "available")
      out$Probability <- vapply(seq_len(nrow(out)), function(j) {
        if (!ok) return(NA_real_)
        take <- true[d$Point] == out$TrueCategory[j] & d$Category == out$ObservedCategory[j]
        sum(d$Probability[take] * w[d$Point[take]])
      }, 0)
      data.frame(Method = method, out, Status = if (ok) "available" else "positive_weight_point_unavailable")
    }))
  }
  table$SourceStatus <- rep(source$status, nrow(table))
  probabilities$SourceStatus <- rep(source$status, nrow(probabilities))
  scores <- do.call(rbind, scores)
  scores$SourceStatus <- rep(source$status, nrow(scores))
  cuts$SourceStatus <- rep(source$status, nrow(cuts))
  marginal$SourceStatus <- rep(source$status, nrow(marginal))
  joint$SourceStatus <- rep(source$status, nrow(joint))
  score_map <- as.data.frame(fit$prep$score_map)
  score_map$CategoryScore <- score_map$InternalScore - cfg$rating_min
  list(table = table, probabilities = probabilities, scores = scores, cuts = cuts,
    marginal = marginal, classification_matrix = joint, population = population,
    rating_plan = plan, calibration = list(model = cfg$model, method = cfg$method,
      offsets = offsets, cumulative_steps = cumulative, score_map = score_map), source = source,
    settings = list(facets = mapping, readiness_policy = readiness_policy,
      category_rule = "lower_inclusive_upper_exclusive", calibration = "fixed_point_estimates",
      repetition = "conditionally_independent_identical_plan", maximum_sum = max_sum),
    notes = c("Lee uses declared category-sum cuts; Rudner approximates ML ability decisions using declared ability cuts.",
      "SEM is inverse-information precision on the ability scale, conditional on fixed calibration and local independence.",
      "Source review status applies to all rows; numerical availability does not establish source or approximation adequacy.",
      "The classification matrix contains joint population probabilities; its true-category row sums are the population proportions, not necessarily one.",
      "Marginal overclassification and underclassification average over the full supplied population. For two categories, divide each joint error probability by its true-category population proportion for a conditional false-pass or false-fail rate; a zero denominator is undefined.",
      "Conditioning on an observed decision instead uses that observed category's column sum in the joint matrix and answers a different question.",
      "Population indices require supplied true-ability weights; calibration, cut and population estimation uncertainty is excluded."))
}

# Shared probability kernel; central moments avoid E(Y^2)-E(Y)^2 cancellation.
mfrm_classification_response <- function(theta, offsets, cumulative) {
  p <- mfrm_jml_probability_bundle(theta + offsets, rep(0L, length(offsets)), "PCM",
    cumulative, criterion_idx = seq_along(offsets))$probs
  if (any(!is.finite(p)) || any(p < 0) || any(abs(rowSums(p) - 1) > 1e-10))
    stop("Classification response probabilities are numerically unresolved.", call. = FALSE)
  k <- seq_len(ncol(p)) - 1L
  expected <- as.vector(p %*% k)
  centered <- matrix(k, nrow(p), length(k), byrow = TRUE) - expected
  information <- sum(p * centered^2)
  list(probs = p, expected = sum(expected), information = information, sem = 1 / sqrt(information))
}

mfrm_classification_point <- function(theta, offsets, cumulative, score_cuts, operational_cuts) {
  x <- mfrm_classification_response(theta, offsets, cumulative)
  pmf <- 1
  for (j in seq_len(nrow(x$probs))) {
    next_mass <- numeric(length(pmf) + ncol(x$probs) - 1L)
    for (k in seq_len(ncol(x$probs))) {
      cells <- seq_along(pmf) + k - 1L
      next_mass[cells] <- next_mass[cells] + pmf * x$probs[j, k]
    }
    pmf <- next_mass
  }
  if (any(!is.finite(pmf)) || any(pmf < 0) || abs(sum(pmf) - 1) > 1e-10)
    stop("Classification summed-score probability mass is numerically unresolved.", call. = FALSE)
  pmf <- pmf / sum(pmf)
  bands <- findInterval(seq_along(pmf) - 1L, score_cuts) + 1L
  lee <- vapply(seq_len(length(score_cuts) + 1L), function(h) sum(pmf[bands == h]), 0)
  rudner <- rep(NA_real_, length(lee))
  if (is.finite(x$sem)) {
    cuts <- c(-Inf, operational_cuts, Inf)
    edges <- (cuts - theta) / x$sem
    rudner <- vapply(seq_along(lee), function(h) {
      lo <- edges[h]; hi <- edges[h + 1L]
      tails <- if (lo >= 0) stats::pnorm(c(lo, hi), lower.tail = FALSE) else
        stats::pnorm(c(hi, lo))
      value <- tails[1L] - tails[2L]
      width <- (cuts[h + 1L] - cuts[h]) / x$sem
      if (is.finite(lo) && is.finite(width) && width > 0 && value < 1e-4*tails[1L]) {
        # Integrate over [0,1], retaining width before subtracting theta and
        # avoiding both CDF cancellation and underflow in the integrand.
        peak <- stats::dnorm(max(lo, min(0, lo + width)), log = TRUE)
        if (is.finite(peak)) {
          mass <- stats::integrate(function(t) exp(stats::dnorm(lo + t*width, log = TRUE) - peak),
            0, 1, rel.tol = 1e-10, abs.tol = 0)$value
          value <- exp(log(width) + peak + log(mass))
        }
      }
      value
    }, 0)
    if (any(!is.finite(rudner)) || any(rudner < 0) || abs(sum(rudner) - 1) > 1e-10)
      stop("Classification normal probabilities are numerically unresolved.", call. = FALSE)
  }
  c(x[c("expected", "information", "sem")], list(pmf = pmf, lee = lee, rudner = rudner))
}
