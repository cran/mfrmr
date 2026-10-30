# Target-specific restricted likelihood; this does not change MINQUE/ANOVA points.
.mfrm_mvgt_profile_minimize <- function(starts, objective, lower = rep(0, length(starts[[1L]]))) {
  attempt <- function(expr) {
    messages <- character()
    ans <- tryCatch(withCallingHandlers(expr, warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) list(par = starts[[1L]], objective = Inf,
      convergence = 1L, message = conditionMessage(e)))
    ans$warnings <- messages
    ans
  }
  fits <- lapply(starts, function(start) attempt(stats::nlminb(start, objective,
    lower = lower, control = list(iter.max = 1000L, eval.max = 2000L, rel.tol = 1e-11))))
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "objective"))]]
  if (is.finite(best$objective)) {
    second <- attempt({
      fit <- minqa::bobyqa(best$par, objective, lower = lower,
        control = list(maxfun = 10000L, rhobeg = .1*max(1, abs(best$par)), rhoend = 1e-7))
      list(par = fit$par, objective = fit$fval, convergence = fit$ierr,
        message = fit$msg, evaluations = fit$feval, optimizer = "bobyqa")
    })
    fits[[length(fits)+1L]] <- second
    if (second$convergence == 0 && second$objective <= best$objective+1e-7) best <- second
  }
  list(best = best, attempts = fits)
}

.mfrm_mvgt_profile_gradient <- function(par, objective, lower = rep(0, length(par)), refinements = 0L) {
  # Richardson checks truncation error, including a one-sided derivative at zero.
  value <- objective(par)
  derivatives <- vapply(seq_along(par), function(j) {
    h <- 1e-5*max(1, abs(par[j]))
    central <- par[j]-h >= lower[j]
    derivative <- function(step) {
      a <- b <- par
      a[j] <- a[j]+step
      if (central) {
        b[j] <- b[j]-step
        (objective(a)-objective(b))/(2*step)
      } else {
        b[j] <- b[j]+2*step
        (-3*value+4*objective(a)-objective(b))/(2*step)
      }
    }
    coarse <- derivative(h)
    best <- c(gradient = NA_real_, error = Inf)
    for (refinement in 0:refinements) {
      fine <- derivative(h/2)
      error <- abs(fine-coarse)/3
      if (is.finite(error) && error < best["error"])
        best <- c(gradient = (4*fine-coarse)/3, error = error)
      if (is.finite(error) && error < .0005) break
      coarse <- fine; h <- h/2
    }
    best
  }, numeric(2))
  gradient <- derivatives[1L, ]
  kkt <- max(ifelse(par > lower+1e-8, abs(gradient), pmax(-gradient, 0)))
  error <- max(derivatives[2L, ])
  list(gradient = gradient, error = error, kkt = kkt,
    valid = is.finite(kkt) && is.finite(error) && kkt < .002 && error < .0005)
}

.mfrm_mvgt_profile_objective <- function(y, groups) {
  d <- data.frame(response = y)
  terms <- paste0("group", seq_len(length(groups)-1L))
  for (j in seq_along(terms)) d[[terms[j]]] <- factor(groups[[j]])
  build <- function() {
    parsed <- lme4::lFormula(stats::reformulate(paste0("(1|", terms, ")"),
      response = "response"), data = d, REML = TRUE)
    order <- match(terms, names(parsed$reTrms$cnms))
    if (anyNA(order) || length(order) != length(parsed$reTrms$theta))
      stop("Unexpected random-effect representation in the REML calculation.", call. = FALSE)
    list(deviance = do.call(lme4::mkLmerDevfun, parsed), order = order)
  }
  state <- build()
  function(v, details = FALSE) {
    # lme4 reorders grouping terms. Variance coordinates permit boundary checks.
    # Invalid trials must not poison the mutable lme4 factorization state.
    unavailable <- if (details) list(objective = Inf, sigma2 = NA_real_) else Inf
    if (length(v) != length(terms) || any(!is.finite(v)) || any(v < 0))
      return(unavailable)
    if (is.null(state)) state <<- build()
    theta <- numeric(length(v)); theta[state$order] <- sqrt(v)
    suspect <- FALSE
    value <- tryCatch(withCallingHandlers(state$deviance(theta),
      warning = function(w) suspect <<- TRUE), error = function(e) {
        state <<- NULL
        warning("REML trial rejected: ", conditionMessage(e), call. = FALSE)
        Inf
      })
    if (suspect || !is.finite(value)) {
      state <<- NULL
      return(unavailable)
    }
    if (!details) return(value)
    environment <- environment(state$deviance)
    sigma2 <- (environment$pp$sqrL(1)+environment$resp$wrss())/(length(y)-1L)
    if (!is.finite(sigma2) || sigma2 <= 0) {
      state <<- NULL
      return(unavailable)
    }
    list(objective = value, sigma2 = sigma2)
  }
}

.mfrm_mvgt_profile_source <- function(g, weights) {
  for (package in c("lme4", "minqa")) if (!requireNamespace(package, quietly = TRUE))
    stop("`method = 'reml_profile'` requires the optional ", package, " package.", call. = FALSE)
  ids <- c(g$design$person, unname(g$design$facets))
  scores <- g$design$scores
  if (!is.data.frame(g$data) || !all(c(ids, scores) %in% names(g$data)) ||
      nrow(g$data) != g$data_usage$counts[["UsedRows"]])
    stop("REML profiles require the original used ratings in the G-study; recreate it from the source data.", call. = FALSE)
  groups <- lapply(g$data[ids], factor)
  subsets <- if (length(ids) == 2L) list(1L, 2L, 1:2) else
    list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
  if (!is.null(g$design$nesting)) {
    labels <- c("Person", names(g$design$facets))
    child <- match(names(g$design$nesting), labels)
    parent <- match(unname(g$design$nesting), labels)
    subsets <- list(1L, parent, c(child, parent), c(1L, parent), 1:3)
  }
  groupings <- lapply(subsets, function(s) .mfrm_mvgt_group(groups[s]))
  names(groupings) <- names(g$components)
  active <- which(weights[scores] != 0)
  y <- as.matrix(g$data[scores[active]])
  if (any(!is.finite(y))) stop("Invalid stored source ratings.", call. = FALSE)
  unit <- apply(abs(y), 2L, max); unit[unit == 0] <- 1
  log_units <- log(unit)+log(abs(weights[scores][active]))
  log_unit <- max(log_units)
  w <- sign(weights[scores][active])*exp(log_units-log_unit)
  # Normalize weighted score contributions together, before squaring their SD.
  y <- sweep(y, 2L, unit, "/")
  y <- drop(sweep(y, 2L, colMeans(y), "-") %*% w)
  y_unit <- max(abs(y))
  scale <- if (y_unit == 0) 0 else y_unit*sqrt(sum((y/y_unit)^2)/(length(y)-1L))
  estimate <- vapply(g$components, function(m) {
    m <- sweep(sweep(m[active, active, drop = FALSE], 1L, unit, "/"), 2L, unit, "/")
    as.numeric(crossprod(w, m %*% w))
  }, numeric(1))
  ans <- list(weights = weights, estimate = estimate, groups = groupings,
    raw_estimate = sign(estimate)*exp(log(abs(estimate))+2*log_unit),
    y = y, score_scale = exp(log_unit+log(scale)),
    valid = FALSE, status = "Degenerate selected score")
  if (!is.finite(scale) || scale <= 0) return(ans)
  y <- y/scale; estimate <- (estimate/scale)/scale
  objective <- .mfrm_mvgt_profile_objective(y, groupings)
  k <- length(groupings)
  initial <- if (is.finite(estimate[k]) && estimate[k] > 0)
    pmax(estimate[-k], 0)/estimate[k] else rep(1, k-1L)
  fit <- .mfrm_mvgt_profile_minimize(list(initial, rep(1, k-1L), rep(.1, k-1L)), objective)
  check <- if (is.finite(fit$best$objective))
    .mfrm_mvgt_profile_gradient(fit$best$par, objective) else list(valid = FALSE)
  ans$y <- y; ans$estimate <- estimate; ans$objective <- objective
  ans$fit <- fit$best; ans$attempts <- fit$attempts; ans$check <- check
  ans$valid <- fit$best$convergence == 0 && check$valid
  ans$status <- if (ans$valid) "Numerically checked REML source" else "Unresolved REML source optimum"
  ans$relative_variances <- c(fit$best$par, 1)
  ans$sigma2 <- objective(fit$best$par, details = TRUE)$sigma2
  ans$components <- ans$sigma2*ans$relative_variances
  ans
}

# A cutoff must be distinguishable in a difference of restricted deviances.
.mfrm_mvgt_profile_cutoff <- function(level, objective) {
  critical <- stats::qchisq(level, 1)
  tolerance <- min(1e-5, 1e-5*critical)
  resolution <- 64*.Machine$double.eps*max(1, abs(objective))
  list(critical = critical, tolerance = tolerance, resolution = resolution,
    valid = is.finite(critical) && critical > 0 && tolerance >= resolution)
}

# Keep normalized profile evidence separate from loss in reported score units.
.mfrm_mvgt_profile_restore <- function(ans, restore) {
  raw <- c(ans$center, ans$bounds, ans$candidate_bounds)
  values <- restore(raw)
  lost <- is.finite(raw) & (!is.finite(values) | (raw != 0 & values == 0))
  values[lost] <- NA_real_
  status <- ans$status
  available <- all(is.finite(ans$bounds))
  identical_errors <- identical(ans$status, "Identical future errors")
  if (available && (any(lost[2:3]) ||
      (!identical_errors && values[2L] >= values[3L]))) {
    values[2:3] <- NA_real_
    status <- "Profile interval unrepresentable in reported units"
    available <- FALSE
  }
  list(values = values, status = status, available = available)
}

.mfrm_mvgt_profile_at <- function(source, rho, error) {
  expand <- function(par) {
    v <- c(0, par, 1); v[1L] <- rho/(1-rho)*sum(error*v)
    v[-length(v)]
  }
  objective <- function(par) source$objective(expand(par))
  k <- length(error)
  fit <- .mfrm_mvgt_profile_minimize(list(source$fit$par[-1L],
    rep(1, k-2L), rep(.1, k-2L)), objective)
  check <- if (is.finite(fit$best$objective))
    .mfrm_mvgt_profile_gradient(fit$best$par, objective) else list(valid = FALSE)
  lr <- fit$best$objective-source$fit$objective
  list(target = rho, lr = lr, fit = fit$best, attempts = fit$attempts, check = check,
    relative_variances = c(expand(fit$best$par), 1),
    valid = source$valid && fit$best$convergence == 0 && check$valid &&
      is.finite(lr) && lr >= -1e-6)
}

.mfrm_mvgt_profile_set <- function(source, error, level) {
  ans <- list(bounds = c(NA_real_, NA_real_), candidate_bounds = c(NA_real_, NA_real_),
    center = NA_real_, status = source$status, error_coefficients = error)
  if (is.null(source$relative_variances)) return(ans)
  v <- source$relative_variances
  center <- v[1L]/(v[1L]+sum(error*v))
  ans$center <- center
  if (!source$valid) return(ans)
  if (!is.finite(center) || center < 0 || center >= 1) {
    ans$status <- if (is.finite(center) && center == 1)
      "REML coefficient center rounds to one at numeric resolution; profile search unavailable" else
      "REML coefficient center outside the finite [0,1) profiling range"
    return(ans)
  }
  cutoff <- .mfrm_mvgt_profile_cutoff(level, source$fit$objective)
  ans$cutoff <- cutoff; ans$critical <- critical <- cutoff$critical
  if (!cutoff$valid) {
    ans$status <- "Likelihood-ratio cutoff below numeric resolution"
    return(ans)
  }
  cache <- new.env(parent = emptyenv())
  evaluate <- function(rho) {
    key <- sprintf("%.17g", rho)
    if (!exists(key, cache, inherits = FALSE))
      assign(key, .mfrm_mvgt_profile_at(source, rho, error), cache)
    get(key, cache, inherits = FALSE)
  }
  statistic <- function(rho) {
    p <- evaluate(rho)
    if (p$valid) p$lr-critical else NA_real_
  }
  grid <- sort(unique(c(seq(0, .95, by = .05), .975, .99, .999, center)))
  value <- vapply(grid, statistic, numeric(1))
  if (is.finite(tail(value, 1)) && tail(value, 1) <= 0) {
    last <- tail(grid, 1)
    for (rho in 1-(1-last)*c(.1, .01, .001)) {
      if (rho >= 1) break
      grid <- c(grid, rho); value <- c(value, statistic(rho))
      if (is.finite(tail(value, 1)) && tail(value, 1) > 0) break
    }
  }
  crossings <- which(is.finite(head(value, -1)) & is.finite(tail(value, -1)) &
    head(value, -1)*tail(value, -1) < 0)
  roots <- grid[which(is.finite(value) & value == 0)]
  failures <- character()
  for (j in crossings) {
    root <- tryCatch(stats::uniroot(statistic, grid[c(j, j+1L)],
      tol = min(1e-8, (grid[j+1L]-grid[j])*1e-6), check.conv = TRUE)$root, error = function(e) {
        failures <<- c(failures, conditionMessage(e)); NA_real_
      })
    if (is.finite(root)) {
      roots <- c(roots, root)
      residual <- statistic(root)
      if (!is.finite(residual) || abs(residual) > cutoff$tolerance)
        failures <- c(failures, "A numerical endpoint did not recover the likelihood-ratio cutoff.")
    }
  }
  roots <- sort(unique(roots))
  breaks <- unique(c(0, roots, 1))
  sets <- matrix(numeric(), 0L, 2L, dimnames = list(NULL, c("Lower", "Upper")))
  for (j in seq_len(length(breaks)-1L)) {
    midpoint <- statistic(mean(breaks[c(j, j+1L)]))
    if (is.finite(midpoint) && midpoint <= 0) sets <- rbind(sets, breaks[c(j, j+1L)])
  }
  profiles <- as.list(cache)
  scan <- do.call(rbind, lapply(profiles, function(p) data.frame(Target = p$target,
    LR = p$lr, KKT = if (is.null(p$check$kkt)) NA_real_ else p$check$kkt,
    GradientError = if (is.null(p$check$error)) NA_real_ else p$check$error, Valid = p$valid)))
  scan <- scan[order(scan$Target), ]; rownames(scan) <- NULL
  ans$critical <- critical; ans$sets <- sets; ans$roots <- roots
  ans$scan <- scan; ans$profiles <- profiles; ans$failures <- failures
  if (nrow(sets)) ans$candidate_bounds <- c(min(sets[, 1L]), max(sets[, 2L]))
  resolved <- !length(failures) && all(scan$Valid) && nrow(sets) > 0L &&
    is.finite(tail(value, 1)) && tail(value, 1) > 0 &&
    ans$candidate_bounds[1L] <= center && ans$candidate_bounds[2L] >= center &&
    ans$candidate_bounds[1L] < ans$candidate_bounds[2L]
  ans$status <- "Unresolved numerical profile; inspect candidate bounds and diagnostics"
  if (resolved) {
    ans$bounds <- ans$candidate_bounds
    ans$status <- if (nrow(sets) == 1L) "Computed profile approximation" else
      "Envelope of disconnected numerical profile sets"
  }
  ans
}

.mfrm_mvgt_profile_error <- function(d, i, metric) {
  g <- d$gstudy
  counts <- as.numeric(d$design_grid[i, ])
  if (!is.null(g$design$nesting)) {
    child <- counts[match(g$design$count_columns[[names(g$design$nesting)]], names(d$design_grid))]
    parent <- counts[match(g$design$count_columns[[unname(g$design$nesting)]], names(d$design_grid))]
    error <- c(0, 0, 0, 1/parent, 1/(child*parent))
    extra <- c(0, 1/parent, 1/(child*parent), 0, 0)
  } else if (length(counts) == 1L) {
    error <- c(0, 0, 1/counts); extra <- c(0, 1/counts, 0)
  } else {
    error <- c(0, 0, 0, 1/counts[1L], 1/counts[2L], 0, 1/prod(counts))
    extra <- c(0, 1/counts[1L], 1/counts[2L], 0, 0, 1/prod(counts), 0)
  }
  if (metric %in% c("Phi", "AbsoluteSEM")) error <- error+extra
  stats::setNames(error, names(g$components))
}

.mfrm_mvgt_profile_table <- function(d, tab, weights, metric, level) {
  g <- d$gstudy
  sem <- metric %in% c("RelativeSEM", "AbsoluteSEM")
  source <- .mfrm_mvgt_profile_source(g, weights)
  if (sem && (!is.finite(source$score_scale) || source$score_scale <= 0) && source$valid)
    stop("SEM units exceed numeric range; rescale the scores or weights.", call. = FALSE)
  rows <- diagnostics <- vector("list", nrow(tab))
  for (i in seq_len(nrow(tab))) {
    error <- .mfrm_mvgt_profile_error(d, i, metric)
    duplicate <- which(vapply(diagnostics[seq_len(i-1L)], function(x) {
      old <- x$error_coefficients
      if (sem) identical(unname(old/tail(old, 1)), unname(error/tail(error, 1))) else
        identical(old, error)
    }, logical(1)))
    ans <- if (length(duplicate)) diagnostics[[duplicate[1L]]] else
      if (sem) .mfrm_mvgt_variance_set(source, error, level) else
        .mfrm_mvgt_profile_set(source, error, level)
    if (sem && length(duplicate)) {
      factor <- unname(tail(error, 1)/tail(ans$error_coefficients, 1))
      for (name in c("bounds", "candidate_bounds", "center", "anchor"))
        if (!is.null(ans[[name]])) ans[[name]] <- ans[[name]]*factor
      if (!is.null(ans$scan)) ans$scan$Target <- ans$scan$Target*factor
      ans$profiles <- lapply(ans$profiles, function(p) {
        p$target <- p$target*factor; p$weights <- error; p
      })
      ans$error_coefficients <- error
    }
    restore <- if (sem) function(x) sqrt(x)*source$score_scale else identity
    reported <- .mfrm_mvgt_profile_restore(ans, restore)
    values <- reported$values
    row <- tab[i, c("Scenario", names(d$design_grid), "Kind", "Score"), drop = FALSE]
    row$Metric <- metric; row$Estimate <- tab[[metric]][i]
    row$PointStatus <- tab[[paste0(metric, "Status")]][i]
    row$IntervalCenter <- values[1L]
    row$Lower <- values[2L]; row$Upper <- values[3L]
    row$CandidateLower <- values[4L]; row$CandidateUpper <- values[5L]
    row$Level <- level; row$Method <- "REML target profile (approximate)"
    row$Status <- reported$status
    row$Validation <- "Experimental chi-square calibration; boundary and finite-sample coverage not guaranteed"
    rows[[i]] <- row; diagnostics[[i]] <- ans
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  # Closures hold live lme4 reference objects; retain only portable evidence.
  source$objective <- NULL
  attr(out, "reml_source") <- source
  attr(out, "source_data") <- g$data
  attr(out, "weights") <- weights; attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- g$design[c("facets", "counts", "method", "balanced", "complete", "score_convention")]
  attr(out, "source_rows") <- g$data_usage$counts
  attr(out, "diagnostics") <- diagnostics
  attr(out, "notes") <- c(
    "Explicit experimental REML target-profile intervals under independent normal random effects and positive residual variance. Missing source cells must be noninformative conditional on the declared fixed assignment. Future plans remain complete balanced plans.",
    "Estimate and PointStatus retain the original ANOVA/MINQUE D-study, including unavailable points. IntervalCenter is the separately estimated REML coefficient. Both concern the same population target; the source point can lie outside the interval. Plots retain the source points, not IntervalCenter.",
    "The chi-square(1) cutoff is an asymptotic approximation, not an exact boundary mixture or finite-sample coverage guarantee. Small checks gave some below-nominal inclusion rates, with substantial Monte Carlo uncertainty; this method is not the automatic default and is not distribution robust.",
    "The numerical inversion checks a grid and brackets cutoff crossings, retaining accepted pieces and their envelope. Multiple starts, BOBYQA and variance-coordinate gradient checks do not prove global optima or exclude narrower undetected pieces. Available describes numerical completion only.",
    "Failed source/profile checks or an unresolved upper tail leave Lower/Upper unavailable. A coefficient center rounded to one also withholds the profile interval with a numeric-resolution status; it does not establish perfect reliability or zero uncertainty. CandidateLower/Upper and diagnostics retain exploratory endpoints; they are not accepted confidence bounds. Uncertainty is never replaced by zero or an automatic [0,1].",
    "Intervals are pointwise for a fixed score/composite and future plan. They do not supply component, covariance-element, plan-difference or simultaneous intervals. Increasing future counts adds no source information.",
    "Raw used ratings are needed once for REML fitting. Save the complete result with saveRDS(); plotting and reloading use saved intervals without refitting. Inspect reml_source and row-aligned diagnostics for objective values, optimizer attempts, gradients and cutoff scans.",
    "A likelihood-ratio cutoff below numerical resolution, or a nonzero interval lost when restoring score units, leaves bounds unavailable with a specific Status. Estimates and representable REML centers are retained; diagnostics$cutoff records the criterion tolerance and resolution.")
  if (sem) {
    notes <- attr(out, "notes")
    notes[2L] <- "Estimate and PointStatus retain the original ANOVA/MINQUE SEM, including unavailable points. IntervalCenter is the separate REML SEM in score units; plots retain the original points. Measurement SEM and its sampling uncertainty are different quantities."
    notes[4L] <- "For SEM, invert the error-variance profile and take square roots. Search outward from the REML center to the first checked cutoff crossing in each direction. This gives a numerical interval around that center; it does not enumerate disconnected remote accepted regions or certify global nuisance optima. All required evaluations and endpoints must pass the same numerical criteria. Available describes numerical completion, not coverage."
    notes[5L] <- "An unresolved source, required profile evaluation, cutoff root or tail leaves Lower/Upper unavailable; resolved CandidateLower/Upper are retained separately. A search limit is not an infinite confidence bound."
    attr(out, "notes") <- c(notes,
      "In a limited normal-model SEM check with 40 samples per condition, true-value inclusion was 92.5-100% for positive-component populations and 82.5-100% at checked boundaries, counting failures as misses. The 82.5% result (MCSE 6.0 percentage points) had 97.5% numerical availability. Numerical completion is therefore not a nominal-coverage guarantee; these checks concern true-value acceptance rather than full endpoint coverage.",
      "SEM diagnostics retain error-variance targets in the standardized scalar-score units used for fitting. Multiply their square roots by reml_source$score_scale to obtain SEM units. Zero nonresidual variances are permitted; residual variance is strictly positive. The constrained calculation uses log1p variance ratios to preserve exact zeros.")
  }
  .mfrm_mvgt_interval_result(out, d, source, "reml_profile")
}

# Proportional future errors share one ratio (G/Phi) or error SD (SEM).
# Keep the observed assignment intact; only the future parent count becomes one.
.mfrm_mvgt_profile_difference <- function(d, tab, weights, metric, level, reference) {
  g <- d$gstudy
  count_column <- if (is.null(g$design$nesting)) names(d$design_grid)[1L] else
    g$design$count_columns[[unname(g$design$nesting)]]
  counts <- d$design_grid[[count_column]]
  base_grid <- d$design_grid[1L, , drop = FALSE]
  base_grid[[count_column]] <- 1
  base_d <- mfrm_d_study(g, base_grid, d$weights)
  base_tab <- base_d$coefficients
  base_tab <- base_tab[base_tab$Kind == tab$Kind[1L] & base_tab$Score == tab$Score[1L], , drop = FALSE]
  base <- .mfrm_mvgt_profile_table(base_d, base_tab, weights, metric, level)
  sem <- metric %in% c("RelativeSEM", "AbsoluteSEM")
  a <- counts[reference]
  factor <- (a-counts)/(sqrt(a)+sqrt(counts))/sqrt(a)/sqrt(counts)
  image <- function(bounds) {
    if (any(!is.finite(bounds))) return(list(lower = rep(NA_real_, length(counts)),
      upper = rep(NA_real_, length(counts)), status = rep(base$Status, length(counts))))
    if (!sem) return(.mfrm_mvgt_count_difference_range(
      stats::qlogis(bounds), counts, a, "Available"))
    lo <- factor*bounds[1L]; hi <- factor*bounds[2L]
    list(lower = pmin(lo, hi), upper = pmax(lo, hi),
      status = ifelse(counts == a, "Identical counts", "Available"))
  }
  bounds <- image(c(base$Lower, base$Upper))
  candidates <- image(c(base$CandidateLower, base$CandidateUpper))
  center <- image(rep(base$IntervalCenter, 2L))$lower
  out <- tab[c("Scenario", names(d$design_grid), "Kind", "Score")]
  out$Reference <- as.integer(reference); out$ReferenceCount <- a
  out$ReferenceValue <- tab[[metric]][reference]; out$Value <- tab[[metric]]
  out$Metric <- paste(metric, "difference")
  out$Estimate <- out$Value-out$ReferenceValue
  out$PointStatus <- ifelse(is.finite(out$Estimate), "Available", "Point difference unavailable")
  out$IntervalCenter <- center
  out$Lower <- bounds$lower; out$Upper <- bounds$upper
  out$CandidateLower <- candidates$lower; out$CandidateUpper <- candidates$upper
  out$Level <- level; out$Method <- "REML profile-set image (approximate)"
  out$Status <- if (base$Available) paste("REML profile image:", bounds$status) else base$Status
  out$Validation <- base$Validation
  out <- out[-reference, , drop = FALSE]; rownames(out) <- NULL
  attr(out, "base_interval") <- base
  for (name in c("reml_source", "source_data", "weights", "assumption", "source_design", "source_rows"))
    attr(out, name) <- attr(base, name)
  attr(out, "reference") <- list(scenario = as.integer(reference), count = a,
    value = tab[[metric]][reference], point_status = tab[[paste0(metric, "Status")]][reference])
  attr(out, "diagnostics") <- lapply(seq_along(counts)[-reference], function(i)
    list(base_scenario = 1L, count = counts[i], reference_count = a,
      sem_multiplier = if (sem) factor[i] else NULL))
  notes <- attr(base, "notes")
  notes[grepl("^Intervals are pointwise", notes)] <-
    "All returned comparisons map one common base interval. Whenever it contains the population target, all mapped intervals contain their differences jointly for this fixed score/composite. The base method's approximate calibration and numerical limitations remain; this is not an exact or distribution-robust coverage guarantee."
  notes[grepl("^Estimate and PointStatus", notes)] <-
    "Estimate and PointStatus retain the original ANOVA/MINQUE plan difference, including unavailable points. IntervalCenter is the separately estimated REML difference. Plots retain source differences, not REML centers."
  notes[grepl("^Raw used ratings", notes)] <-
    "Raw used ratings are needed once for REML fitting. Save the complete result with saveRDS(); plotting and reloading do not refit. Inspect base_interval for the shared profile scans, objective values, optimizer attempts and gradient checks."
  attr(out, "notes") <- c(
    if (sem) "SEM differences multiply one shared base SEM interval by 1/sqrt(comparison count) - 1/sqrt(reference count), reversing bounds for negative multipliers. Negative differences favor the comparison plan." else
      "G/Phi differences map the entire shared base interval, including interior extrema. They do not subtract marginal confidence limits. Positive differences favor the comparison plan.",
    "base_interval retains the canonical future count of one, its accepted/candidate bounds and fitting diagnostics. In nested plans it retains the common future child count; the observed incomplete or unequal assignment is never averaged or replaced. ReferenceCount counts parents. Only proportional future error vectors are supported.",
    "Estimate is the original ANOVA/MINQUE plan difference and IntervalCenter is the separate REML difference. A missing source point stays missing. Failed base intervals stay unavailable even for identical plans; candidate images are exploratory, not accepted bounds.",
    notes)
  .mfrm_mvgt_interval_result(out, d, attr(base, "reml_source"), "reml_profile")
}

.mfrm_mvgt_profile_bind <- function(results) {
  columns <- unique(unlist(lapply(results, names), use.names = FALSE))
  rows <- lapply(results, function(x) {
    x <- as.data.frame(x)
    for (name in setdiff(columns, names(x))) x[[name]] <- NA
    x[columns]
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  for (result in results) for (name in setdiff(names(attributes(result)), c("names", "row.names", "class")))
    if (is.null(attr(out, name))) attr(out, name) <- attr(result, name)
  attr(out, "diagnostics") <- unlist(lapply(results, function(x) {
    diagnostics <- attr(x, "diagnostics")
    if (is.null(diagnostics)) rep(list(NULL), nrow(x)) else diagnostics
  }), recursive = FALSE)
  attr(out, "notes") <- unique(unlist(lapply(results, function(x) attr(x, "notes")), use.names = FALSE))
  # Methods can have different base intervals/centers; retain each full result.
  attr(out, "method_results") <- stats::setNames(results,
    vapply(results, function(x) unique(x$MethodCode), character(1)))
  class(out) <- c("mfrm_dstudy_intervals", "data.frame")
  out
}
