# Direct REML profiles for comparison-minus-reference G, Phi and SEM.
# The chi-square cutoff is experimental, especially at singular zero contrasts.
.mfrm_mvgt_contrast_value <- function(v, a, b) {
  u <- v[1]
  ea <- sum(a * v)
  eb <- sum(b * v)
  u / (u + ea) * (ea - eb) / (u + eb)
}

.mfrm_mvgt_contrast_gradient <- function(v, a, b) {
  u <- v[1]
  ea <- sum(a * v)
  eb <- sum(b * v)
  ans <- u * a / (u + ea)^2 - u * b / (u + eb)^2
  ans[1] <- eb / (u + eb)^2 - ea / (u + ea)^2
  ans
}

.mfrm_mvgt_contrast_person_roots <- function(ea, eb, delta) {
  stopifnot(is.finite(ea), is.finite(eb), ea > 0, eb > 0, is.finite(delta), delta != 0)
  limit <- tanh((log(ea) - log(eb)) / 4)
  if (sign(delta) != sign(limit) || abs(delta) > abs(limit) * (1 + 32 * .Machine$double.eps)) {
    return(numeric())
  }
  # log roots are symmetric about log(sqrt(ea*eb)); avoid quadratic cancellation.
  x <- max(0, (abs(limit) - abs(delta)) / abs(delta))
  distance <- 2 * asinh(sqrt(x) * cosh((log(ea) - log(eb)) / 4))
  exp((log(ea) + log(eb)) / 2 + c(-distance, distance))
}

.mfrm_mvgt_contrast_range <- function(a, b) {
  active <- a + b > 0
  ratio <- a[active] / b[active]
  c(min(0, tanh(log(min(ratio)) / 4)), max(0, tanh(log(max(ratio)) / 4)))
}

.mfrm_mvgt_contrast_zero_possible <- function(a, b) {
  h <- b - a
  tail(h, 1) == 0 || (tail(h, 1) > 0 && any(head(h, -1) < 0)) ||
    (tail(h, 1) < 0 && any(head(h, -1) > 0))
}

.mfrm_mvgt_contrast_bound_attained <- function(a, b, target) {
  target == 0 || abs(tanh(log(tail(a, 1) / tail(b, 1)) / 4) - target) < 32 * .Machine$double.eps
}

.mfrm_mvgt_contrast_starts <- function(source, a, b, target) {
  k <- length(a)
  starts <- list()
  oriented_a <- if (target > 0) a else b
  oriented_b <- if (target > 0) b else a
  h <- oriented_a - ((1 + abs(target)) / (1 - abs(target)))^2 * oriented_b
  pivot <- which.max(h[-k])
  limits <- .mfrm_mvgt_contrast_range(a, b)
  at_boundary <- any(abs(target - limits) < 32 * .Machine$double.eps) &&
    .mfrm_mvgt_contrast_bound_attained(a, b, target)
  if (at_boundary) {
    ratio <- oriented_a / oriented_b
    residual_ratio <- tail(ratio, 1)
    v <- source$relative_variances
    v[which(is.finite(ratio) & abs(ratio / residual_ratio - 1) > 64 * .Machine$double.eps)] <- 0
    roots <- .mfrm_mvgt_contrast_person_roots(sum(a * v), sum(b * v), target)
    return(lapply(roots, function(u) {
      v[1] <- u
      log1p(v[-k])
    }))
  }
  if (h[pivot] <= 0 && h[k] < 0) {
    return(starts)
  }
  for (original in list(source$relative_variances, c(1, rep(.1, k - 2), 1), rep(1, k))) {
    v <- original
    if (sum(h * v) <= 0) {
      if (h[pivot] > 0) {
        v[pivot] <- v[pivot] + (-sum(h * v) + .1) / h[pivot]
      } else {
        # A feasible residual-dominated ratio is reached by reducing other errors.
        negative <- which(h < 0)
        v[negative] <- v[negative] * h[k] / (-2 * sum(h[negative] * v[negative]))
      }
    }
    roots <- .mfrm_mvgt_contrast_person_roots(sum(a * v), sum(b * v), target)
    for (u in roots) {
      v[1] <- u
      starts[[length(starts) + 1L]] <- log1p(v[-k])
    }
  }
  starts
}

.mfrm_mvgt_contrast_kkt <- function(par, objective, constraint_jacobian, lower = rep(0, length(par))) {
  check <- .mfrm_mvgt_profile_gradient(par, objective, lower, refinements = 4L)
  gradient <- check$gradient
  h <- constraint_jacobian(par)
  free <- par > lower + 1e-7
  divisor <- sum(h[free]^2)
  if (divisor > 1e-20) {
    lambda <- -sum(h[free] * gradient[free]) / divisor
  } else {
    # If the equality fixes a boundary component, its multiplier is not unique.
    positive <- !free & h > 0
    negative <- !free & h < 0
    lo <- max(c(-Inf, -gradient[positive] / h[positive]))
    hi <- min(c(Inf, -gradient[negative] / h[negative]))
    lambda <- if (lo <= hi) max(lo, min(hi, 0)) else NA_real_
  }
  stationarity <- gradient + lambda * h
  kkt <- max(ifelse(free, abs(stationarity), pmax(-stationarity, 0)))
  list(
    kkt = kkt, error = check$error, multiplier = lambda, gradient = gradient,
    constraint_gradient = h, stationarity = stationarity,
    valid = is.finite(kkt) && kkt < .002 && check$error < .0005
  )
}

.mfrm_mvgt_contrast_profile_at <- function(source, a, b, target, algorithm = "NLOPT_LD_SLSQP", zero_person_branch = TRUE,
                                           starts_extra = list(source$relative_variances), previous = NULL) {
  if (!is.null(previous) && is.null(previous$attempts)) previous <- NULL
  if (target != 0 && any(abs(target - .mfrm_mvgt_contrast_range(a, b)) < 32 * .Machine$double.eps) &&
    .mfrm_mvgt_contrast_bound_attained(a, b, target)) {
    return(.mfrm_mvgt_contrast_boundary_profile_at(source, a, b, target))
  }
  k <- length(a)
  decode <- function(x) c(expm1(x), 1)
  objective <- function(x) source$objective(expm1(x))
  gradient <- function(x) .mfrm_mvgt_profile_gradient(x, objective)$gradient
  if (target == 0) {
    # Equal-error branch; the U=0 branch is fitted independently below.
    h <- (a - b) / max(abs(a - b))
    constraint <- function(x) sum(h * decode(x))
    jacobian <- function(x) h[-k] * exp(x)
    starts <- lapply(list(source$fit$par, rep(.1, k - 1), rep(1, k - 1)), log1p)
    if (!.mfrm_mvgt_contrast_zero_possible(a, b)) starts <- list()
  } else {
    constraint <- function(x) .mfrm_mvgt_contrast_value(decode(x), a, b) - target
    jacobian <- function(x) .mfrm_mvgt_contrast_gradient(decode(x), a, b)[-k] * exp(x)
    starts <- if (is.null(previous)) .mfrm_mvgt_contrast_starts(source, a, b, target) else list()
    starts <- c(starts, lapply(starts_extra, function(v) log1p(v[-k])))
  }
  fits <- lapply(starts, function(start) {
    messages <- character()
    fit <- tryCatch(withCallingHandlers(nloptr::nloptr(
      x0 = start, eval_f = objective,
      eval_grad_f = if (algorithm == "NLOPT_LD_SLSQP") gradient else NULL,
      lb = rep(0, k - 1), eval_g_eq = constraint,
      eval_jac_g_eq = if (algorithm == "NLOPT_LD_SLSQP") jacobian else NULL,
      opts = list(
        algorithm = algorithm, xtol_rel = 1e-9, ftol_abs = 1e-9,
        tol_constraints_eq = 1e-10, maxeval = 3000, print_level = 0
      )
    ), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    if (inherits(fit, "error")) {
      return(list(
        objective = Inf, par = start,
        relative_variances = decode(start), convergence = NA_integer_, message = conditionMessage(fit),
        constraint = constraint(start), check = list(kkt = Inf, error = Inf, valid = FALSE),
        branch = "failed optimization attempt", warnings = unique(messages), valid = FALSE
      ))
    }
    v <- decode(fit$solution)
    check <- .mfrm_mvgt_contrast_kkt(fit$solution, objective, jacobian)
    residual <- constraint(fit$solution)
    list(
      objective = fit$objective, par = fit$solution, relative_variances = v,
      convergence = fit$status, message = fit$message, iterations = fit$iterations,
      constraint = residual, check = check, branch = if (target == 0) "equal errors" else if (v[1] <= sqrt(sum(a * v) * sum(b * v))) "lower Person root" else "upper Person root",
      warnings = unique(messages), valid = fit$status %in% 1:4 && abs(residual) < 1e-8 && check$valid
    )
  })
  if (!is.null(previous)) fits <- c(previous$attempts, fits)
  if (target == 0 && zero_person_branch) {
    p <- .mfrm_mvgt_variance_at(source, 0, c(1, rep(0, k - 1)))
    fits[[length(fits) + 1L]] <- list(
      objective = p$fit$objective, par = p$fit$par,
      relative_variances = p$relative_variances, convergence = p$fit$convergence,
      check = p$check, branch = "zero Person variance", constraint = 0, valid = p$valid
    )
  }
  values <- vapply(fits, `[[`, numeric(1), "objective")
  if (!length(values)) stop("No feasible starting point for the target difference.")
  accepted <- which(vapply(fits, `[[`, logical(1), "valid"))
  best <- fits[[if (length(accepted)) accepted[which.min(values[accepted])] else which.min(values)]]
  gap <- best$objective - min(values)
  lr <- best$objective - source$fit$objective
  list(
    target = target, best = best, attempts = fits, lr = lr,
    unchecked_objective_gap = gap,
    valid = source$valid && best$valid && is.finite(gap) && gap <= 1e-6 && is.finite(lr) && lr >= -1e-6,
    algorithm = algorithm, a = a, b = b
  )
}

.mfrm_mvgt_contrast_boundary_profile_at <- function(source, a, b, target) {
  k <- length(a)
  ratio <- a / b
  residual_ratio <- tail(ratio, 1)
  forbidden <- which(is.finite(ratio) & abs(ratio / residual_ratio - 1) > 64 * .Machine$double.eps)
  free <- setdiff(seq_len(k - 1), c(1, forbidden))
  decode <- function(x) {
    v <- numeric(k)
    v[k] <- 1
    v[free] <- expm1(x)
    v[1] <- sqrt(sum(a * v) * sum(b * v))
    v
  }
  objective <- function(x) source$objective(head(decode(x), -1))
  if (length(free)) {
    fitted <- .mfrm_mvgt_profile_minimize(lapply(list(
      source$fit$par[free],
      rep(.1, length(free)), rep(1, length(free))
    ), log1p), objective)
    fit <- fitted$best
    check <- .mfrm_mvgt_profile_gradient(fit$par, objective, refinements = 4L)
  } else {
    fit <- list(par = numeric(), objective = objective(numeric()), convergence = 0L)
    fitted <- list(attempts = list(fit))
    check <- list(kkt = 0, error = 0, valid = TRUE)
  }
  residual <- .mfrm_mvgt_contrast_value(decode(fit$par), a, b) - target
  best <- list(
    objective = fit$objective, par = fit$par, relative_variances = decode(fit$par),
    convergence = fit$convergence, check = check, branch = "Attained coefficient boundary",
    constraint = residual, valid = fit$convergence == 0 && check$valid && abs(residual) < 1e-10
  )
  lr <- best$objective - source$fit$objective
  list(
    target = target, best = best, attempts = fitted$attempts, lr = lr,
    valid = source$valid && best$valid && is.finite(lr) && lr >= -1e-6, a = a, b = b, free = free, forbidden = forbidden,
    boundary_reduction = TRUE
  )
}


.mfrm_mvgt_sem_contrast_value <- function(v, a, b) {
  ea <- sum(a * v)
  eb <- sum(b * v)
  (eb - ea) / (sqrt(eb) + sqrt(ea))
}

.mfrm_mvgt_sem_contrast_starts <- function(source, a, b, target) {
  k <- length(a)
  h <- sign(target) * (b - a)
  pivot <- which.max(h[-k])
  if (h[pivot] <= 0 && h[k] <= 0) {
    return(list())
  }
  starts <- list()
  for (v in list(source$relative_variances, c(1, rep(.1, k - 2), 1), rep(1, k))) {
    if (sum(h * v) <= 0) {
      if (h[pivot] > 0) {
        v[pivot] <- v[pivot] + (.1 - sum(h * v)) / h[pivot]
      } else {
        negative <- which(h < 0)
        v[negative] <- v[negative] * h[k] / (-2 * sum(h[negative] * v[negative]))
      }
    }
    difference <- .mfrm_mvgt_sem_contrast_value(v, a, b)
    logscale <- 2 * (log(abs(target)) - log(abs(difference)))
    starts[[length(starts) + 1]] <- c(log1p(v[-k]), logscale)
  }
  starts
}

.mfrm_mvgt_sem_contrast_objective <- function(source, x) {
  k <- length(x)
  v <- expm1(x[-k])
  if (any(!is.finite(v)) || !is.finite(x[k])) {
    return(Inf)
  }
  at <- source$objective(v, details = TRUE)
  if (!is.finite(at$objective)) {
    return(Inf)
  }
  logr <- x[k] - log(at$sigma2)
  at$objective + (length(source$y) - 1) * (logr + expm1(-logr))
}

.mfrm_mvgt_sem_contrast_profile_at <- function(
  source, a, b, target, algorithm = "NLOPT_LD_SLSQP",
  starts_extra = list(list(relative_variances = source$relative_variances, sigma2 = source$sigma2)),
  previous = NULL
) {
  if (!is.null(previous) && is.null(previous$attempts)) previous <- NULL
  k <- length(a)
  h <- b - a
  zero_possible <- .mfrm_mvgt_contrast_zero_possible(a, b)
  if (target == 0) {
    if (!zero_possible) {
      return(list(
        target = target, valid = TRUE, infeasible = TRUE, lr = Inf,
        status = "Zero contrast is not attained with positive residual variance"
      ))
    }
    p <- .mfrm_mvgt_contrast_profile_at(source, a, b, 0, algorithm, zero_person_branch = FALSE)
    p$best$sigma2 <- source$objective(p$best$relative_variances[-k], details = TRUE)$sigma2
    p$infeasible <- FALSE
    return(p)
  }
  starts <- .mfrm_mvgt_sem_contrast_starts(source, a, b, target)
  if (!length(starts)) {
    return(list(
      target = target, valid = TRUE, infeasible = TRUE, lr = Inf,
      status = "Contrast has the wrong sign for these error weights"
    ))
  }
  if (!is.null(previous)) starts <- list()
  starts <- c(starts, lapply(starts_extra, function(p) c(log1p(head(p$relative_variances, -1)), log(p$sigma2))))
  objective <- function(x) .mfrm_mvgt_sem_contrast_objective(source, x)
  lower <- c(rep(0, k - 1), -Inf)
  gradient <- function(x) .mfrm_mvgt_profile_gradient(x, objective, lower)$gradient
  constraint <- function(x) exp(x[k] / 2) * .mfrm_mvgt_sem_contrast_value(c(expm1(x[-k]), 1), a, b) - target
  jacobian <- function(x) {
    v <- c(expm1(x[-k]), 1)
    ea <- sum(a * v)
    eb <- sum(b * v)
    scale <- exp(x[k] / 2)
    c(
      scale * (b[-k] / (2 * sqrt(eb)) - a[-k] / (2 * sqrt(ea))) * exp(x[-k]),
      scale * .mfrm_mvgt_sem_contrast_value(v, a, b) / 2
    )
  }
  attempts <- lapply(starts, function(start) {
    messages <- character()
    fit <- tryCatch(withCallingHandlers(nloptr::nloptr(
      x0 = start, eval_f = objective,
      eval_grad_f = if (algorithm == "NLOPT_LD_SLSQP") gradient else NULL, lb = lower,
      eval_g_eq = constraint, eval_jac_g_eq = if (algorithm == "NLOPT_LD_SLSQP") jacobian else NULL,
      opts = list(
        algorithm = algorithm, xtol_rel = 1e-9, ftol_abs = 1e-9,
        tol_constraints_eq = 1e-10, maxeval = 3000, print_level = 0
      )
    ), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    if (inherits(fit, "error")) {
      return(list(
        objective = Inf, par = start,
        relative_variances = c(expm1(start[-k]), 1), sigma2 = exp(start[k]),
        convergence = NA_integer_, message = conditionMessage(fit), constraint = constraint(start),
        check = list(kkt = Inf, error = Inf, valid = FALSE), warnings = unique(messages), valid = FALSE
      ))
    }
    check <- .mfrm_mvgt_contrast_kkt(fit$solution, objective, jacobian, lower)
    residual <- constraint(fit$solution)
    list(
      objective = fit$objective, par = fit$solution,
      relative_variances = c(expm1(fit$solution[-k]), 1), sigma2 = exp(fit$solution[k]),
      convergence = fit$status, message = fit$message, iterations = fit$iterations,
      constraint = residual, check = check, warnings = unique(messages),
      valid = fit$status %in% 1:4 && abs(residual) < 1e-8 && check$valid
    )
  })
  if (!is.null(previous)) attempts <- c(previous$attempts, attempts)
  values <- vapply(attempts, `[[`, numeric(1), "objective")
  accepted <- which(vapply(attempts, `[[`, logical(1), "valid"))
  best <- attempts[[if (length(accepted)) accepted[which.min(values[accepted])] else which.min(values)]]
  gap <- best$objective - min(values)
  lr <- best$objective - source$fit$objective
  list(
    target = target, best = best, attempts = attempts, lr = lr,
    unchecked_objective_gap = gap, infeasible = FALSE,
    valid = source$valid && best$valid && is.finite(gap) && gap <= 1e-6 && is.finite(lr) && lr >= -1e-6,
    algorithm = algorithm, a = a, b = b
  )
}


# Finite scans retain detected pieces; they do not prove global or remote-set exclusion.
.mfrm_mvgt_contrast_profile_set <- function(source, a, b, level = .95, sem = FALSE, initial = list()) {
  critical <- stats::qchisq(level, 1)
  unavailable <- list(
    bounds = c(NA_real_, NA_real_), candidate_bounds = c(NA_real_, NA_real_),
    valid = FALSE, center = NA_real_, status = source$status, a = a, b = b, critical = critical, sem = sem
  )
  if (is.null(source$relative_variances)) {
    return(unavailable)
  }
  center <- if (sem) {
    sqrt(source$sigma2) * .mfrm_mvgt_sem_contrast_value(source$relative_variances, a, b)
  } else {
    .mfrm_mvgt_contrast_value(source$relative_variances, a, b)
  }
  unavailable$center <- center
  if (!source$valid || !is.finite(center)) {
    return(unavailable)
  }
  if (all(a == b)) {
    return(list(
      bounds = c(0, 0), candidate_bounds = c(0, 0), valid = TRUE,
      center = 0, status = "Identical future errors", sets = matrix(c(0, 0), 1, 2),
      a = a, b = b, critical = critical, sem = sem
    ))
  }
  cutoff <- .mfrm_mvgt_profile_cutoff(level, source$fit$objective)
  unavailable$cutoff <- cutoff
  if (!cutoff$valid) {
    unavailable$status <- "Likelihood-ratio cutoff below numeric resolution"
    return(unavailable)
  }
  if (sem) {
    h <- b - a
    bounds <- c(if (any(h < 0)) -Inf else 0, if (any(h > 0)) Inf else 0)
    attained <- rep(.mfrm_mvgt_contrast_zero_possible(a, b), 2)
    anchor <- sqrt(sum(a * source$components)) / 2
    grid <- sort(unique(c(center + anchor * c(-2, -1, -.5, -.25, 0, .25, .5, 1, 2), 0)))
    grid <- grid[grid >= bounds[1] & grid <= bounds[2]]
  } else {
    bounds <- .mfrm_mvgt_contrast_range(a, b)
    anchor <- diff(bounds)
    attained <- vapply(bounds, function(x) .mfrm_mvgt_contrast_bound_attained(a, b, x), logical(1))
    grid <- sort(unique(c(bounds[1] + diff(bounds) * seq(.01, .99, length.out = 17), center, 0)))
    grid <- grid[grid >= bounds[1] & grid <= bounds[2]]
    grid <- sort(unique(c(grid, bounds[attained])))
  }
  # The grid explicitly includes zero; remove floating duplicates of that node.
  grid[abs(grid) < 64 * .Machine$double.eps * anchor] <- 0
  grid <- sort(unique(grid))
  cache <- new.env(parent = emptyenv())
  visited <- new.env(parent = emptyenv())
  failures <- character()
  extras <- if (sem) {
    list(list(relative_variances = source$relative_variances, sigma2 = source$sigma2))
  } else {
    list(source$relative_variances)
  }
  for (p in initial) assign(sprintf("%.17g", p$target), p, cache)
  evaluate <- function(target) {
    key <- sprintf("%.17g", target)
    old <- if (exists(key, cache, inherits = FALSE)) get(key, cache, inherits = FALSE) else NULL
    if (is.null(old) || (target != 0 && !exists(key, visited, inherits = FALSE))) {
      p <- tryCatch(if (sem) {
        .mfrm_mvgt_sem_contrast_profile_at(source, a, b, target, starts_extra = extras, previous = old)
      } else {
        .mfrm_mvgt_contrast_profile_at(source, a, b, target, starts_extra = extras, previous = old)
      }, error = function(e) {
        list(target = target, valid = FALSE, lr = NA_real_, status = conditionMessage(e))
      })
      assign(key, p, cache)
    }
    assign(key, TRUE, visited)
    get(key, cache, inherits = FALSE)
  }
  statistic <- function(target) {
    p <- evaluate(target)
    if (p$valid) p$lr - critical else NA_real_
  }
  null <- evaluate(0)
  if (isTRUE(null$valid) && !isTRUE(null$infeasible)) {
    if (sem) {
      extras <- c(extras, list(null$best))
    } else {
      candidates <- Filter(function(p) isTRUE(p$valid) && identical(p$branch, "equal errors"), null$attempts)
      if (length(candidates)) {
        best <- candidates[[which.min(vapply(candidates, `[[`, numeric(1), "objective"))]]
        extras <- c(extras, list(best$relative_variances))
      }
    }
  }
  value <- vapply(grid, statistic, numeric(1))
  tail_resolved <- rep(FALSE, 2)
  for (side in 1:2) {
    index <- if (side == 1) 1L else length(grid)
    edge <- grid[index]
    v <- value[index]
    if (!is.na(v) && v > 0) {
      tail_resolved[side] <- TRUE
      next
    }
    if (is.finite(bounds[side]) && edge == bounds[side] && !is.na(v)) {
      tail_resolved[side] <- attained[side]
      next
    }
    for (j in 1:12) {
      target <- if (is.finite(bounds[side])) {
        bounds[side] + (edge - bounds[side]) * 10^(-j)
      } else {
        edge + if (side == 1) -anchor * 2^j else anchor * 2^j
      }
      if (target == bounds[side] && !attained[side]) break
      v <- statistic(target)
      grid <- c(grid, target)
      value <- c(value, v)
      if (is.na(v)) break
      if (v > 0) {
        tail_resolved[side] <- TRUE
        break
      }
    }
  }
  ord <- order(grid)
  grid <- grid[ord]
  value <- value[ord]
  crossings <- which(!is.na(head(value, -1)) & !is.na(tail(value, -1)) &
    sign(head(value, -1)) != sign(tail(value, -1)))
  roots <- grid[which(!is.na(value) & value == 0)]
  for (j in crossings) {
    pair <- grid[c(j, j + 1)]
    root <- tryCatch(stats::uniroot(statistic, pair,
      tol = min(1e-9, abs(diff(pair)) * 1e-7),
      check.conv = TRUE
    )$root, error = function(e) {
      failures <<- c(failures, conditionMessage(e))
      NA_real_
    })
    if (is.finite(root)) {
      residual <- statistic(root)
      if (!is.finite(residual) || abs(residual) > cutoff$tolerance) {
        failures <- c(failures, "Endpoint did not recover the likelihood-ratio cutoff")
      } else {
        roots <- c(roots, root)
      }
    }
  }
  roots <- sort(unique(roots))
  breaks <- sort(unique(c(min(grid), roots, max(grid))))
  sets <- matrix(numeric(), 0, 2, dimnames = list(NULL, c("Lower", "Upper")))
  if (length(breaks) > 1) {
    for (j in seq_len(length(breaks) - 1)) {
      at <- statistic(mean(breaks[c(j, j + 1)]))
      if (is.finite(at) && at <= 0) sets <- rbind(sets, breaks[c(j, j + 1)])
    }
  }
  profiles <- as.list(cache)[ls(visited)]
  scan <- do.call(rbind, lapply(profiles, function(p) {
    data.frame(
      Target = p$target, LR = p$lr,
      Valid = p$valid, Infeasible = isTRUE(p$infeasible)
    )
  }))
  scan <- scan[order(scan$Target), ]
  rownames(scan) <- NULL
  candidate <- if (nrow(sets)) c(min(sets[, 1]), max(sets[, 2])) else c(NA_real_, NA_real_)
  valid <- source$valid && all(scan$Valid) && all(tail_resolved) && !length(failures) &&
    all(is.finite(candidate)) && candidate[1] <= center && candidate[2] >= center &&
    candidate[1] < candidate[2]
  list(
    bounds = if (valid) candidate else c(NA_real_, NA_real_), candidate_bounds = candidate,
    valid = valid, center = center, sets = sets, roots = roots, scan = scan, profiles = profiles,
    target_bounds = bounds, attained = attained, tail_resolved = tail_resolved, failures = failures,
    critical = critical, cutoff = cutoff, sem = sem, a = a, b = b,
    status = if (!valid) {
      "Unresolved difference profile; inspect candidate bounds and diagnostics"
    } else if (nrow(sets) > 1) {
      "Envelope of disconnected numerical profile sets"
    } else {
      "Computed difference profile approximation"
    }
  )
}

.mfrm_mvgt_contrast_table <- function(d, tab, weights, metric, level, reference) {
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    stop("Nonproportional REML comparisons require the optional nloptr package.", call. = FALSE)
  }
  g <- d$gstudy
  source <- .mfrm_mvgt_profile_source(g, weights)
  sem <- metric %in% c("RelativeSEM", "AbsoluteSEM")
  scale <- if (sem) source$score_scale else 1
  if (source$valid && (!is.finite(scale) || scale <= 0)) {
    stop("Interval units exceed numeric range; rescale the scores or weights.", call. = FALSE)
  }
  a <- .mfrm_mvgt_profile_error(d, reference, metric)
  indices <- seq_len(nrow(tab))[-reference]
  rows <- diagnostics <- vector("list", length(indices))
  for (j in seq_along(indices)) {
    i <- indices[j]
    b <- .mfrm_mvgt_profile_error(d, i, metric)
    duplicate <- which(vapply(diagnostics[seq_len(j - 1L)], function(x) {
      identical(x$b, b)
    }, logical(1)))
    ans <- if (length(duplicate)) {
      diagnostics[[duplicate[1L]]]
    } else {
      .mfrm_mvgt_contrast_profile_set(source, a, b, level, sem)
    }
    row <- tab[i, c("Scenario", names(d$design_grid), "Kind", "Score"), drop = FALSE]
    row$Reference <- as.integer(reference)
    row$ReferenceValue <- tab[[metric]][reference]
    row$Value <- tab[[metric]][i]
    row$Metric <- paste(metric, "difference")
    row$Estimate <- row$Value - row$ReferenceValue
    row$PointStatus <- if (is.finite(row$Estimate)) "Available" else "Point difference unavailable"
    reported <- .mfrm_mvgt_profile_restore(ans, function(x) x*scale)
    values <- reported$values
    row$IntervalCenter <- values[1L]
    row$Lower <- values[2L]; row$Upper <- values[3L]
    row$CandidateLower <- values[4L]; row$CandidateUpper <- values[5L]
    row$Level <- level
    row$Method <- "REML difference profile (approximate)"
    row$Status <- reported$status
    row$Validation <- "Experimental chi-square calibration; singularities, boundaries and finite-sample coverage not established"
    rows[[j]] <- row
    diagnostics[[j]] <- ans
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  source$objective <- NULL
  attr(out, "reml_source") <- source
  attr(out, "source_data") <- g$data
  attr(out, "weights") <- weights
  attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- g$design[c("facets", "counts", "method", "balanced", "complete", "score_convention")]
  attr(out, "source_rows") <- g$data_usage$counts
  attr(out, "reference") <- list(
    scenario = as.integer(reference),
    counts = d$design_grid[reference, , drop = FALSE], value = tab[[metric]][reference],
    point_status = tab[[paste0(metric, "Status")]][reference]
  )
  attr(out, "diagnostics") <- diagnostics
  attr(out, "notes") <- c(
    "Differences are comparison minus reference. Positive G/Phi differences and negative SEM differences favor the comparison plan. The two plans share the same source likelihood; marginal confidence limits are not subtracted.",
    "Explicit experimental REML difference profiles assume independent normal random effects and strictly positive residual variance. Missing source cells must be noninformative conditional on the declared fixed assignment. Future plans are complete balanced designs; both facet counts may change.",
    "Estimate and PointStatus retain the original ANOVA/MINQUE difference. IntervalCenter is the separately estimated REML difference; it does not replace an unavailable original point. Plots show original differences, which can lie outside the separately calculated intervals.",
    "The chi-square(1) cutoff is uncalibrated at general boundaries and singularities. At zero G/Phi difference the null includes both zero Person variance and equal error variances; their intersection is singular. Finite-sample coverage for these comparisons is not established, so this is not an automatic default or a recommended substitute for a qualified procedure.",
    "Multiple starts include both positive Person roots and the equal-error null continuation. Independent constraint and stationarity checks, a grid and checked cutoff roots support numerical completion; they do not certify global optima or exclude narrow undetected accepted pieces. Lower/Upper enclose detected pieces, retained in diagnostics$sets. Comparisons are pointwise, not simultaneous.",
    "An unresolved source, required evaluation, cutoff crossing or tail leaves Lower/Upper unavailable. CandidateLower/Upper retain exploratory endpoints, not accepted confidence bounds. A finite search limit is not an infinite endpoint. Available describes numerical completion, not coverage.",
    "Save the complete result with saveRDS(). Plotting and reloading do not refit. reml_source and row-aligned diagnostics retain optimizer attempts, gradients, constraint residuals and likelihood scans. SEM diagnostics use standardized score units; multiply targets, bounds and sets by reml_source$score_scale to recover score units.",
    "A likelihood-ratio cutoff below numerical resolution, or a nonzero interval lost when restoring score units, leaves bounds unavailable with a specific Status. Estimates and representable REML centers are retained; diagnostics$cutoff records the criterion tolerance and resolution."
  )
  .mfrm_mvgt_interval_result(out, d, source, "reml_profile")
}
