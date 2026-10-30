# Common-plan expected scores and the fixed-Person calibration/refit pipeline.
# These helpers deliberately do not use posterior EAPs or reporting replacements.
mfrm_fair_reference_fit <- function(fit) {
  if ((!inherits(fit, "mfrm_fit") && !mfrm_is_jeffreys_fit(fit)) || inherits(fit, "mfrm_imported_fit") ||
      mfrm_extended_fit(fit) || mfrm_has_jml_adjustment(fit) ||
      !identical(fit$config$method, "JML") ||
      !fit$config$model %in% c("RSM", "PCM") ||
      length(fit$config$interaction_specs) ||
      !identical(fit$config$facet_shrinkage %||% "none", "none"))
    stop("Common-plan scores require native additive RSM/PCM ordinary or explicit Jeffreys JML without shrinkage.", call. = FALSE)
}

mfrm_fair_fit_identity <- function(fit) {
  if (!mfrm_is_jeffreys_fit(fit)) return(prediction_fit_identity(fit))
  mfrm_checkpoint_fingerprint(list(estimator = fit$config$estimator, class = class(fit),
    config = fit$config, prep = fit$prep, opt = fit$opt, facets = fit$facets,
    steps = fit$steps, probabilities = fit$probabilities), canonical = TRUE)
}

mfrm_fair_reference_plan <- function(fit, plan) {
  mfrm_fair_reference_fit(fit)
  facets <- fit$config$facet_names
  if ("Weight" %in% facets)
    stop("The reference plan reserves the column name 'Weight'.", call. = FALSE)
  if (!is.data.frame(plan) || !nrow(plan) || anyDuplicated(names(plan)) ||
      !setequal(setdiff(names(plan), "Weight"), facets))
    stop("`reference_plan` must contain every fitted facet and optional Weight only.", call. = FALSE)
  out <- as.data.frame(plan[facets], stringsAsFactors = FALSE)
  for (f in facets) {
    x <- plan[[f]]
    if (!(is.character(x) || is.factor(x) || is.numeric(x)) ||
        !is.null(dim(x)) || anyNA(x) || any(!as.character(x) %in% fit$config$facet_levels[[f]]))
      stop("Unknown or invalid reference-plan level in facet `", f, "`.", call. = FALSE)
    out[[f]] <- as.character(x)
  }
  if (anyDuplicated(out)) stop("Use weights for duplicate reference-plan cells.", call. = FALSE)
  w <- plan$Weight %||% rep(1, nrow(out))
  if (!is.numeric(w) || !is.null(dim(w)) || any(!is.finite(w)) || any(w <= 0))
    stop("Reference-plan weights must be finite and positive.", call. = FALSE)
  out$Weight <- w
  out$NormalizedWeight <- (w / max(w)) / sum(w / max(w))
  out
}

mfrm_fair_reference_kernel <- function(fit, plan, par = fit$opt$par) {
  cfg <- fit$config
  pa <- expand_params(par, build_param_sizes(cfg), cfg)
  idx <- list(score_k = rep(0L, nrow(plan)), facets = lapply(cfg$facet_names,
    function(f) match(plan[[f]], cfg$facet_levels[[f]])), interactions = list())
  names(idx$facets) <- cfg$facet_names
  if (anyNA(unlist(idx$facets))) stop("Reference plan contains unknown fitted levels.", call. = FALSE)
  offsets <- compute_base_eta(idx, pa, cfg)
  steps <- if (cfg$model == "RSM")
    matrix(pa$steps, nrow(plan), cfg$n_cat - 1L, byrow = TRUE) else
    pa$steps_mat[idx$facets[[cfg$step_facet]], , drop = FALSE]
  cumulative <- t(vapply(seq_len(nrow(plan)), function(j) c(0, cumsum(steps[j, ])),
                         numeric(cfg$n_cat)))
  list(offsets = offsets, cumulative = cumulative)
}

mfrm_fair_reference_scores <- function(fit, plan) {
  kernel <- mfrm_fair_reference_kernel(fit, plan)
  p <- fit$facets$person
  theta <- p$PrimaryEstimate %||% p$Estimate
  score <- vapply(theta, function(t) {
    if (is.na(t)) return(NA_real_)
    if (is.infinite(t)) return(if (t > 0) fit$prep$rating_max else fit$prep$rating_min)
    pr <- mfrm_jml_probability_bundle(t + kernel$offsets, rep(0L, nrow(plan)),
      "PCM", kernel$cumulative, criterion_idx = seq_len(nrow(plan)))$probs
    sum(plan$NormalizedWeight * drop(pr %*% (seq_len(ncol(pr)) - 1L))) + fit$prep$rating_min
  }, 0)
  joint <- fit$config$boundary_audit$joint_additive
  structural <- joint$target_status
  invalid <- !mfrm_is_jeffreys_fit(fit) && (!isTRUE(joint$complete) ||
    any(structural$ParameterClass != "person" &
      structural$CandidateStatus %in% c("unbounded_low", "unbounded_high", "unbounded_both")))
  if (invalid) score[] <- NA_real_
  data.frame(Person = as.character(p$Person), Estimate = unname(score),
    Status = ifelse(!is.finite(score), "unavailable_calibration",
      ifelse(is.finite(theta), "point_estimate", "boundary_limit")),
    Target = "common_plan_expected_average", Units = "internal_score",
    Uncertainty = "not_computed", stringsAsFactors = FALSE)
}

# Check the actual ordinary-JML equations, not just optim's return code.
mfrm_fair_refit_check <- function(fit, gradient_tol) {
  mfrm_fair_reference_fit(fit)
  cfg <- fit$config
  if (!identical(as.integer(fit$opt$convergence), 0L)) stop("Optimizer did not converge.", call. = FALSE)
  audit <- cfg$estimability_audit
  if (!isTRUE(audit$complete) || length(audit$zero_coordinates) || nrow(audit$null_blocks))
    stop("Additive identification is unresolved or deficient.", call. = FALSE)
  if (mfrm_is_jeffreys_fit(fit)) {
    d <- mfrm_fair_jeffreys_design(fit)
    objective <- mfrm_fair_jeffreys_objective(d$category, d$offset, d$y)
    z <- objective$evaluate(fit$opt$par, TRUE)
    grad <- if (length(z$gradient)) max(abs(z$gradient)) else 0
    if (!is.finite(grad) || grad > gradient_tol ||
        any(z$probabilities <= 0 | z$probabilities >= 1))
      stop("Jeffreys penalized stationarity/interior check failed.", call. = FALSE)
    return(list(max_gradient = grad, nll = z$nll, penalized_objective = z$value,
      information_logdet = z$logdet, solve_residual = z$solve_residual,
      estimator = "jeffreys_canonical_joint"))
  }
  joint <- cfg$boundary_audit$joint_additive
  if (!isTRUE(joint$complete)) stop("Joint boundary audit is incomplete.", call. = FALSE)
  tab <- joint$target_status
  if (any(tab$ParameterClass != "person" &
          tab$CandidateStatus %in% c("unbounded_low", "unbounded_high", "unbounded_both")))
    stop("Structural calibration has a certified recession direction.", call. = FALSE)
  idx <- build_indices(fit$prep, cfg$step_facet, cfg$slope_facet, cfg$interaction_specs)
  sizes <- build_param_sizes(cfg)
  g <- mfrm_grad_jml(fit$opt$par, idx, cfg, sizes)
  grad <- if (length(g)) max(abs(g)) else 0
  nll <- mfrm_loglik_jml(fit$opt$par, idx, cfg, sizes)
  if (!is.finite(nll) || !is.finite(grad) || grad > gradient_tol)
    stop("JML stationarity check failed (maximum absolute gradient = ", signif(grad, 5), ").", call. = FALSE)
  list(max_gradient = grad, nll = nll)
}

mfrm_fair_anchor_values <- function(fit, anchors) {
  person <- fit$facets$person
  other <- fit$facets$others
  tab <- rbind(data.frame(Facet = "Person", Level = person$Person,
    Estimate = person$PrimaryEstimate %||% person$Estimate),
    as.data.frame(other[c("Facet", "Level", "Estimate")]))
  values <- vapply(seq_len(nrow(anchors)), function(i) {
    j <- which(tab$Facet == anchors$Facet[i] & tab$Level == anchors$Level[i])
    if (length(j) == 1L) tab$Estimate[j] else NA_real_
  }, 0)
  if (any(!is.finite(values))) stop("An estimated anchor is missing or nonfinite in the calibration fit.", call. = FALSE)
  values
}

# The explicit row maps identify response EVENTS, not merely matching labels.
# Matching cell labels cannot tell whether two rows are the same occasion.
mfrm_fair_response_map <- function(source, generator, map, label) {
  n <- nrow(source$prep$data)
  if (!is.numeric(map) || length(map) != n || any(!is.finite(map)) ||
      any(map != floor(map)) || any(map < 1 | map > nrow(generator$prep$data)) || anyDuplicated(map))
    stop(label, " must identify one distinct generating row per prepared source row.", call. = FALSE)
  columns <- c("Person", source$config$facet_names)
  a <- source$prep$data[columns]; b <- generator$prep$data[map, columns, drop = FALSE]
  if (!all(vapply(columns, function(f) identical(as.character(a[[f]]), as.character(b[[f]])), TRUE)) ||
      !identical(as.numeric(source$prep$data$Score), as.numeric(generator$prep$data$Score[map])))
    stop(label, " does not match the prepared Person, facet and observed-score rows.", call. = FALSE)
  as.integer(map)
}

mfrm_fair_refit <- function(source, scores, anchors = source$config$replay_inputs$anchors) {
  args <- source$config$replay_inputs
  dat <- as.data.frame(source$prep$data[c("Person", source$config$facet_names)])
  # Keep factor level order: it defines the saved affine parameter coordinates.
  names(dat)[1] <- args$person
  dat[[args$score]] <- scores
  if (!is.null(args$weight)) dat[[args$weight]] <- 1
  args$package_version <- NULL; args$data <- dat; args$anchors <- anchors
  args$attach_diagnostics <- FALSE
  # Pin the complete source category support when a draw omits a category.
  args$rating_min <- source$prep$rating_min; args$rating_max <- source$prep$rating_max
  fit <- do.call(fit_mfrm, args)
  if (!identical(source$prep$levels, fit$prep$levels) ||
      !identical(build_param_sizes(source$config), build_param_sizes(fit$config)) ||
      !isTRUE(all.equal(source$prep$score_map, fit$prep$score_map, check.attributes = FALSE)))
    stop("Refitting changed parameter identities or category support.", call. = FALSE)
  fit
}

mfrm_fair_generator_probs <- function(fit) {
  cfg <- fit$config; pa <- expand_params(fit$opt$par, build_param_sizes(cfg), cfg)
  p <- fit$facets$person
  theta <- p$PrimaryEstimate %||% p$Estimate
  if (any(!is.finite(theta)))
    stop("The generating fit needs finite primary estimates for every Person.", call. = FALSE)
  idx <- build_indices(fit$prep, cfg$step_facet, cfg$slope_facet, cfg$interaction_specs)
  eta <- theta[match(as.character(fit$prep$data$Person), p$Person)] + compute_base_eta(idx, pa, cfg)
  if (cfg$model == "RSM") category_prob_rsm(eta, c(0, cumsum(pa$steps))) else
    category_prob_pcm(eta, t(apply(pa$steps_mat, 1L, function(s) c(0, cumsum(s)))), idx$step_idx)
}

mfrm_fair_bootstrap_bounds <- function(estimate, truth, draws, level) {
  roots <- sweep(draws, 2L, truth, "-")
  missing <- !is.finite(roots); low <- high <- roots
  low[missing] <- -Inf; high[missing] <- Inf
  tail <- (1 - level) / 2
  lower <- estimate - apply(high, 2L, stats::quantile, probs = 1 - tail, type = 1, names = FALSE)
  upper <- estimate - apply(low, 2L, stats::quantile, probs = tail, type = 1, names = FALSE)
  list(lower = lower, upper = upper, available = colSums(!missing), roots = roots)
}

# Binomial ranks enclose each fixed generating-distribution quantile. Allocate
# four one-sided errors to cover both basic endpoints for ONE target jointly.
# F(q-) <= p <= F(q) also makes the bounds conservative for ties/discrete roots.
# Missing roots stay unknown: monotone order-statistic envelopes cover every
# completion, without assuming that numerical failures are missing at random.
mfrm_fair_bootstrap_mc <- function(estimate, roots, level, mc_level) {
  n <- nrow(roots)
  p <- c((1 + level) / 2, (1 - level) / 2)
  error <- (1 - mc_level) / 4
  lower_rank <- stats::qbinom(error, n, p)
  upper_rank <- stats::qbinom(1 - error, n, p) + 1
  targets <- colnames(roots) %||% as.character(seq_len(ncol(roots)))
  rows <- lapply(seq_len(ncol(roots)), function(j) {
    x <- roots[, j]; missing <- !is.finite(x)
    lo <- hi <- x; lo[missing] <- -Inf; hi[missing] <- Inf
    lo <- c(-Inf, sort(lo), Inf); hi <- c(-Inf, sort(hi), Inf)
    data.frame(Target = targets[j], Endpoint = c("lower", "upper"),
      Probability = p, RankLower = lower_rank, RankUpper = upper_rank,
      Lower = estimate[j] - hi[upper_rank + 1],
      Upper = estimate[j] - lo[lower_rank + 1],
      MCLevel = mc_level, Planned = n, Unresolved = sum(missing),
      RankLimited = lower_rank == 0 | upper_rank == n + 1,
      row.names = NULL)
  })
  do.call(rbind, rows)
}
