# Invert signed MLS null intervals for balanced crossed G/Phi.
# No finite-sample coverage claim follows from computing all polynomial roots.
.mfrm_mvgt_mls_matrix <- function(df, signs, tail) {
  pos <- which(signs > 0)
  neg <- which(signs < 0)
  k <- length(df)
  G <- 1 - df / stats::qchisq(tail, df, lower.tail = FALSE)
  H <- df / stats::qchisq(tail, df) - 1
  K <- matrix(0, k, k)
  diag(K)[pos] <- G[pos]^2
  diag(K)[neg] <- H[neg]^2
  for (i in pos) {
    for (j in neg) {
      f <- stats::qf(tail, df[i], df[j], lower.tail = FALSE)
      K[i, j] <- K[j, i] <- (f - 2 + 1 / f - G[i]^2 * f - H[j]^2 / f) / 2
    }
  }
  if (length(pos) > 1 && length(neg)) {
    for (q in seq_len(length(pos) - 1)) {
      for (r in seq.int(q + 1, length(pos))) {
        i <- pos[q]
        j <- pos[r]
        d <- df[i] + df[j]
        g <- 1 - d / stats::qchisq(tail, d, lower.tail = FALSE)
        K[i, j] <- K[j, i] <- (g^2 * d^2 / (df[i] * df[j]) - G[i]^2 * df[i] / df[j] -
          G[j]^2 * df[j] / df[i]) / (2 * (length(pos) - 1))
      }
    }
  }
  K * outer(signs, signs)
}

.mfrm_mvgt_mls_prepare <- function(map, error, level = .95) {
  a <- unname(map$A[1, ])
  b <- a + drop(error %*% map$A)
  coefficient_scale <- max(abs(c(a, b)))
  a <- a / coefficient_scale
  b <- b / coefficient_scale
  crossing <- a / b
  crossing[!is.finite(crossing)] <- NA_real_
  breaks <- sort(unique(c(0, 1, crossing[crossing > 0 & crossing < 1 & !is.na(crossing)])))
  matrices <- lapply(seq_len(length(breaks) - 1), function(i) {
    signs <- sign(a - mean(breaks[i + 0:1]) * b)
    list(
      lower = .mfrm_mvgt_mls_matrix(map$df, signs, (1 - level) / 2),
      upper = .mfrm_mvgt_mls_matrix(map$df, -signs, (1 - level) / 2)
    )
  })
  list(
    a = a, b = b, df = map$df, breaks = breaks, crossing = crossing,
    matrices = matrices, level = level, error = error
  )
}

.mfrm_mvgt_mls_poly_roots <- function(p, lo, hi) {
  while (length(p) > 1 && tail(p, 1) == 0) p <- head(p, -1)
  if (length(p) == 1) {
    return(numeric())
  }
  r <- polyroot(p / max(abs(p)))
  r <- Re(r[abs(Im(r)) <= 1e-9 * pmax(1, abs(Re(r)))])
  r <- r[is.finite(r) & r > lo & r < hi]
  if (length(r) && any(vapply(r, function(x) {
    abs(sum(p * x^(seq_along(p) - 1))) > 1e-9 * sum(abs(p))
  }, logical(1)))) {
    stop("Unresolved MLS polynomial roots", call. = FALSE)
  }
  r
}

.mfrm_mvgt_mls_null_interval <- function(ms, problem, rho, snap = FALSE) {
  h <- problem$a - rho * problem$b
  # Only known sign-change coordinates are snapped, not arbitrary targets.
  if (snap) h[which(!is.na(problem$crossing) & problem$crossing == rho)] <- 0
  point <- sum(ms * h)
  point + c(
    -.mfrm_mvgt_mls_width(ms, problem$df, h, (1 - problem$level) / 2),
    .mfrm_mvgt_mls_width(ms, problem$df, -h, (1 - problem$level) / 2)
  )
}

.mfrm_mvgt_mls_invert <- function(ms, problem) {
  stopifnot(all(is.finite(ms)), all(ms > 0))
  ms[problem$a == 0 & problem$b == 0] <- 0
  ms <- ms / max(ms)
  point <- sum(ms * problem$a) / sum(ms * problem$b)
  polynomial <- function(K, C) {
    Q <- crossprod(C, K %*% C)
    c(Q[1, 1], Q[1, 2] + Q[2, 1], Q[2, 2])
  }
  pieces <- list()
  nodes <- problem$breaks
  unresolved <- list()
  polys <- list()
  for (i in seq_along(problem$matrices)) {
    lo <- problem$breaks[i]
    hi <- problem$breaks[i + 1]
    # Form polynomials near the relevant endpoint, before squaring. Expanding
    # in rho loses narrow intervals near one through cancellation of O(1) terms.
    origin <- if (is.finite(point) && abs(point - hi) < abs(point - lo)) hi else lo
    C <- cbind(ms * (problem$a - origin * problem$b), -ms * problem$b)
    step <- max(abs(C[, 1])) / max(abs(C[, 2]))
    if (!is.finite(step) || step == 0) step <- hi - lo
    C[, 2] <- C[, 2] * step
    C <- C / max(abs(C))
    mu <- colSums(C)
    mu2 <- c(mu[1]^2, 2 * mu[1] * mu[2], mu[2]^2)
    vl <- polynomial(problem$matrices[[i]]$lower, C)
    vu <- polynomial(problem$matrices[[i]]$upper, C)
    roots <- unlist(lapply(list(mu, vl, vu, mu2 - vl, mu2 - vu),
      .mfrm_mvgt_mls_poly_roots,
      lo = (lo - origin) / step, hi = (hi - origin) / step
    ))
    coordinates <- origin + step * roots
    if (any(coordinates <= lo | coordinates >= hi)) {
      stop("MLS root separation is not representable on the coefficient scale", call. = FALSE)
    }
    cuts <- sort(unique(c(lo, hi, coordinates)))
    polys[[i]] <- list(lower = vl, upper = vu, origin = origin, step = step, cuts = cuts)
    for (j in seq_len(length(cuts) - 1)) {
      at <- mean(cuts[j + 0:1])
      ci <- .mfrm_mvgt_mls_null_interval(ms, problem, at)
      if (!all(is.finite(ci))) {
        unresolved[[length(unresolved) + 1L]] <- cuts[j + 0:1]
        next
      }
      if (ci[1] <= 0 && ci[2] >= 0) {
        pieces[[length(pieces) + 1L]] <- cuts[j + 0:1]
      }
    }
    nodes <- c(nodes, cuts)
  }
  nodes <- sort(unique(nodes))
  node_ci <- t(vapply(nodes, function(r) .mfrm_mvgt_mls_null_interval(ms, problem, r, snap = TRUE), numeric(2)))
  node_scale <- vapply(nodes, function(r) sum(abs(ms * (problem$a - r * problem$b))), numeric(1))
  node_tol <- 256 * .Machine$double.eps * pmax(node_scale, rowSums(abs(node_ci)))
  node_ok <- apply(is.finite(node_ci), 1, all)
  accepted <- node_ok & node_ci[, 1] <= node_tol & node_ci[, 2] >= -node_tol
  endpoints <- c(unlist(pieces), nodes[accepted])
  valid <- !length(unresolved) && all(node_ok)
  list(
    valid = valid, empty = valid && !length(endpoints),
    bounds = if (valid && length(endpoints)) range(endpoints) else c(NA_real_, NA_real_),
    pieces = if (length(pieces)) do.call(rbind, pieces) else matrix(numeric(), 0, 2),
    nodes = nodes, node_accepted = accepted, node_ci = node_ci,
    unresolved = unresolved, polynomials = polys
  )
}

.mfrm_mvgt_mls_complete <- function(fit, ms, problem) {
  if (fit$valid) {
    return(list(
      available = !fit$empty, bounds = fit$bounds,
      formula_complete = TRUE, undefined = matrix(numeric(), 0, 2), status = "Computed MLS set"
    ))
  }
  failed <- list(
    available = FALSE, bounds = c(NA_real_, NA_real_), formula_complete = FALSE,
    undefined = matrix(numeric(), 0, 2), status = "Unresolved numerical inversion"
  )
  if (!length(fit$polynomials) || !length(fit$unresolved)) {
    return(failed)
  }
  undefined <- do.call(rbind, fit$unresolved)
  if (any(!is.finite(unlist(fit$polynomials)))) {
    return(failed)
  }
  # Certify the reason with the original scalar formula, not a matrix label.
  ms[problem$a == 0 & problem$b == 0] <- 0
  for (r in rowMeans(undefined)) {
    ci <- .mfrm_mvgt_mls_null_interval(ms / max(ms), problem, r)
    if (!anyNA(ci) || any(is.infinite(ci))) {
      return(failed)
    }
  }
  bad_node <- which(!apply(is.finite(fit$node_ci), 1, all))
  if (length(bad_node)) {
    undefined <- rbind(
      undefined,
      cbind(fit$nodes[bad_node], fit$nodes[bad_node])
    )
  }
  allpoints <- c(fit$pieces, fit$nodes[fit$node_accepted], undefined)
  bounds <- range(allpoints)
  stopifnot(all(is.finite(bounds)), bounds[1] >= 0, bounds[2] <= 1)
  list(
    available = TRUE, bounds = bounds, formula_complete = FALSE, undefined = undefined,
    status = "Conservative enclosure includes undefined MLS widths"
  )
}

.mfrm_mvgt_mls_coefficient_table <- function(d, tab, p, metric, level) {
  rows <- diagnostics <- vector("list", nrow(tab))
  for (i in seq_len(nrow(tab))) {
    error <- .mfrm_mvgt_profile_error(d, i, metric)
    duplicate <- which(vapply(diagnostics[seq_len(i - 1L)], function(x) {
      identical(x$error, error)
    }, logical(1)))
    ans <- if (length(duplicate)) {
      diagnostics[[duplicate[1L]]]
    } else {
      tryCatch(
        {
          if (any(!is.finite(p$ms) | p$ms <= 0)) {
            stop("Positive mean squares required", call. = FALSE)
          }
          problem <- .mfrm_mvgt_mls_prepare(p$design, error, level)
          fit <- .mfrm_mvgt_mls_invert(p$ms, problem)
          complete <- .mfrm_mvgt_mls_complete(fit, p$ms, problem)
          complete$numerical_complete <- isTRUE(fit$valid) || isTRUE(complete$available)
          complete$empty <- isTRUE(fit$empty)
          if (complete$empty) complete$status <- "Empty MLS confidence set"
          complete$fit <- fit
          complete$problem <- problem
          complete
        },
        error = function(e) {
          list(
            available = FALSE, bounds = c(NA_real_, NA_real_),
            formula_complete = FALSE, numerical_complete = FALSE, empty = FALSE,
            status = conditionMessage(e)
          )
        }
      )
    }
    ans$error <- error
    row <- tab[i, c("Scenario", names(d$design_grid), "Kind", "Score"), drop = FALSE]
    row$Metric <- metric
    row$Estimate <- tab[[metric]][i]
    row$PointStatus <- tab[[paste0(metric, "Status")]][i]
    row$Lower <- ans$bounds[1L]
    row$Upper <- ans$bounds[2L]
    row$Level <- level
    row$Method <- "MLS null inversion (approximate)"
    row$NumericalComplete <- ans$numerical_complete
    row$FormulaComplete <- ans$formula_complete
    row$Empty <- ans$empty
    row$Status <- ans$status
    row$Validation <- "Normal-model MLS approximation; not exact or distribution robust"
    rows[[i]] <- row
    diagnostics[[i]] <- ans
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "weights") <- p$weights
  attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- d$gstudy$design[c("facets", "counts", "method", "balanced", "complete", "score_convention")]
  attr(out, "source_rows") <- d$gstudy$data_usage$counts
  attr(out, "mean_squares") <- list(
    normalized = p$ms, scale = p$restore,
    df = p$design$df, component_map = p$design$A
  )
  attr(out, "diagnostics") <- diagnostics
  attr(out, "notes") <- c(
    "Approximate intervals for individual G/Phi targets from a balanced Person-by-two-crossed-random-facet source under independent normal effects. Select a prespecified score or composite. These are not nested, incomplete-source or plan-difference intervals.",
    "For each candidate coefficient, invert the unrestricted Ting MLS interval for its signed linear null in the expected mean squares. Source degrees of freedom determine the asymmetric chi-square/F corrections. Original ANOVA estimates and point failures are preserved; no component clipping or fitted-nuisance bootstrap is used.",
    "The returned bounds enclose all retained pieces. FormulaComplete is FALSE if a squared MLS width is negative and therefore undefined for some candidates. Those candidates are conservatively retained, with their ranges in diagnostics; the width is not replaced by zero. NumericalComplete concerns the inversion, while Available concerns returned bounds. Neither certifies coverage.",
    "When root separation cannot be represented on the coefficient scale, bounds are unavailable with an explicit Status and NumericalComplete=FALSE; the point estimate remains. A genuine parameter-space boundary is not itself a failure.",
    "Small sources can give very wide or full [0,1] intervals, including with large proposed future counts. More future ratings do not add source information. Limited normal-model evaluations support this approximation but do not establish uniform finite-sample coverage; see the help for results and scope.",
    "The reported sampling evaluations concern 95% intervals. At very low requested confidence levels, squared chi-square corrections can lose their direction; the two-term F limiting identity is not guaranteed there.",
    "Intervals are pointwise for one prespecified score and plan. Curve ribbons add no simultaneous coverage, and interval overlap does not test a plan difference. An empty set is distinct from zero uncertainty or a failed numerical calculation. Save the complete result with saveRDS() to retain notes and diagnostics without recalculation."
  )
  .mfrm_mvgt_interval_result(out, d, p, "mls")
}
