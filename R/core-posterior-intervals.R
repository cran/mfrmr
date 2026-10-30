# Continuous equal-tail intervals, conditional on the supplied calibration/prior.
# Quadrature point masses approximate moments, but their step CDF need not
# approximate posterior quantiles accurately even when those moments converge.
mfrmr_person_posterior_intervals <- function(idx, config, params, base_eta,
                                             person_labels, population_spec = NULL,
                                             interval_level = 0.95, cut_scores = NULL) {
  basis <- resolve_person_quadrature_basis(list(nodes = c(0,1), weights = c(1,1)),
    population_spec %||% materialize_population_spec(config, params), person_count = length(person_labels))
  mu <- if (isTRUE(basis$transformed)) basis$mu else rep(0, length(person_labels))
  sigma <- if (isTRUE(basis$transformed)) basis$sigma else 1
  weights <- idx$weight %||% rep(1, length(idx$score_k))
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("Posterior intervals require finite non-negative response weights.", call. = FALSE)
  }
  alpha <- (1 - interval_level) / 2
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 0.5) {
    stop("`interval_level` must be strictly between zero and one.", call. = FALSE)
  }
  ids <- unique(idx$person)
  groups <- split(seq_along(idx$person), factor(idx$person, levels = ids))
  result <- lapply(seq_along(groups), function(i) {
    rows <- groups[[i]]; person <- ids[i]
    local_idx <- list(person = rep(1L,length(rows)), score_k = idx$score_k[rows],
      weight = weights[rows], step_idx = idx$step_idx[rows], slope_idx = idx$slope_idx[rows])
    evaluate <- mfrmr_adaptive_person_kernel(local_idx, config, params, base_eta[rows], mu[person], sigma)
    location <- mfrmr_adaptive_posterior_location(evaluate)
    mode <- location$mode; scale <- location$scale
    peak_likelihood <- evaluate(mode)$log_likelihood
    peak_prior <- stats::dnorm(mode, log = TRUE)
    log_density <- function(u) {
      z <- mode + scale * u
      # Remove the common likelihood offset before adding the prior density.
      (evaluate(z, include_moments = FALSE)$log_likelihood - peak_likelihood) +
        (stats::dnorm(z, log = TRUE) - peak_prior)
    }
    integral <- function(lower, upper) {
      if (lower == upper) return(0)
      stats::integrate(function(u) exp(log_density(u)), lower, upper,
        rel.tol = 1e-10, abs.tol = 0, subdivisions = 1000L)$value
    }
    bounded <- FALSE
    for (limit in c(8,16,32,64,128,256)) {
      left <- integral(-limit,0); right <- integral(0,limit); mass <- left + right
      edge <- mode + scale * c(-limit,limit)
      score <- scale * evaluate(edge)$score
      if (!is.finite(mass) || mass <= 0 || score[1L] <= 0 || score[2L] >= 0) next
      # Strict concavity gives tangent-exponential upper bounds on omitted tails.
      log_tail <- log_density(c(-limit,limit)) - log(abs(score))
      if (all(log_tail - log(mass) < log(min(1e-12,alpha*1e-10)/2))) {
        bounded <- TRUE; break
      }
    }
    if (!bounded) stop("Could not bound the omitted posterior mass for an interval.", call. = FALSE)
    lower_tail <- function(u) (if (u <= 0) integral(-limit,u) else left + integral(0,u)) / mass
    upper_tail <- function(u) (if (u >= 0) integral(u,limit) else right + integral(u,0)) / mass
    # Solve the upper survival probability directly to avoid 1-alpha cancellation.
    lower <- stats::uniroot(function(u) lower_tail(u)-alpha, c(-limit,limit),
      tol=1e-10, check.conv=TRUE)$root
    upper <- stats::uniroot(function(u) upper_tail(u)-alpha, c(-limit,limit),
      tol=1e-10, check.conv=TRUE)$root
    answer <- mu[person] + sigma * (mode + scale * c(Lower=lower,Upper=upper))
    if (!is.null(cut_scores)) {
      cuts <- ((unname(cut_scores)-mu[person])/sigma-mode)/scale
      answer <- c(answer, mfrm_cut_probabilities(integral, cuts, mass, -limit, limit))
    }
    answer
  })
  do.call(rbind,result)
}
