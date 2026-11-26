.mfrm_mvgt_two_facet_intervals <- function(d, tab, p, metric, level, methods, reference = NULL) {
    setup <- if ("mean_square_region" %in% methods)
        .mfrm_mvgt_ms_region_setup(p$design, level)
    else NULL
    sem <- metric %in% c("RelativeSEM", "AbsoluteSEM")
    nesting <- d$gstudy$design$nesting
    errors <- lapply(seq_len(nrow(tab)), function(i) unname(.mfrm_mvgt_profile_error(d, i, metric)))
    rows <- diagnostics <- list()
    for (method in methods) for (i in seq_len(nrow(tab))) {
        error <- errors[[i]]
        ans <- if (method == "approximate_f")
            .mfrm_mvgt_approximate_f(p$ms, p$design, error, level)
        else .mfrm_mvgt_ms_region(p$ms, setup, error, ratio = !sem)
        row <- tab[i, c("Scenario", names(d$design_grid), "Kind", "Score"), drop = FALSE]
        row$Metric <- metric
        row$Estimate <- tab[[metric]][i]
        row$PointStatus <- tab[[paste0(metric, "Status")]][i]
        bounds <- if (sem)
            sqrt(ans$bounds) * sqrt(p$restore)
        else ans$bounds
        row$Lower <- bounds[1L]
        row$Upper <- bounds[2L]
        row$Level <- level
        row$MethodCode <- method
        row$Method <- if (method == "approximate_f")
            "Approximate F"
        else "Mean-square region"
        row$Status <- ans$status
        if (sem && all(is.finite(ans$bounds)) && (any(!is.finite(bounds)) || any(bounds == 0 & ans$bounds >
            0) || (ans$bounds[1L] < ans$bounds[2L] && bounds[1L] == bounds[2L]))) {
            row$Lower <- row$Upper <- NA_real_
            row$Status <- "Unrepresentable SEM endpoints"
        }
        rows[[length(rows) + 1L]] <- row
        diagnostics[[length(diagnostics) + 1L]] <- ans
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    attr(out, "weights") <- p$weights
    attr(out, "assumption") <- "normal"
    attr(out, "source_design") <- d$gstudy$design[c("facets", "counts", "method", "balanced", "complete",
        "score_convention")]
    attr(out, "source_rows") <- d$gstudy$data_usage$counts
    attr(out, "diagnostics") <- diagnostics
    attr(out, "mean_squares") <- list(normalized = p$ms, scale = p$restore, df = p$design$df, component_map = p$design$A)
    attr(out, "notes") <- c("Experimental methods for a complete balanced Person-by-two-random-facet ANOVA source, one selected score or prespecified composite, and independent normal effects. Finite-sample evaluation is limited; these are not distribution-robust intervals.",
        "More future facet levels do not add source information. Methods use the same source and point estimate; interval overlap and method agreement do not test plan differences or establish coverage.",
        if ("approximate_f" %in% methods) paste("Approximate F freezes Satterthwaite denominator degrees of freedom at the unrestricted ANOVA coefficient and inverts F inequalities within [0,1].",
            "Wong and McGraw advise against negative component estimates; calculable results are retained here as an experimental extension with NegativeComponent flagged. Invalid moments or empty intersections remain unavailable."),
        if ("mean_square_region" %in% methods) paste("Mean-square region projects one independent chi-square simultaneous region onto nonnegative component vectors. Coefficients use linear-fractional optimization; SEM uses linear error-variance projections followed by square roots.",
            "Under the normal model with positive mean-square expectations, its coverage is at least the requested level, jointly over future plans for this fixed score. It can be very conservative; coefficient intervals can be almost [0,1]. An empty feasible region or failed numerical check remains unavailable."),
        "Negative source component estimates are preserved. PointStatus and Status describe different outputs. Available records numerical availability, not validation; inspect Validation and diagnostics, including unavailable rows.",
        "Numerically unresolved F quantiles, mean-square region limits or SEM transformations leave missing bounds with an explicit Status. Point estimates remain unchanged. These checks concern arithmetic precision, not a minimum sample size or model adequacy.",
        "Save the full result with saveRDS(). Use as.data.frame(result) for interval values and attr(result, \"notes\") for interpretation notes.")
    if (!is.null(nesting) && "approximate_f" %in% methods)
        attr(out, "notes") <- c(attr(out, "notes"), "Nested approximate F is derived from the expected mean squares and population target. Printed absolute-agreement formulas in Wong/McGraw (1999, p.278) omit terms required by that map; the help gives the algebra and limiting-case checks. This is not a literal reproduction of those printed formulas.")
    if (!is.null(reference)) {
        attr(out, "individual_intervals") <- out[c("Scenario", "Estimate", "PointStatus", "Lower", "Upper",
            "Status")]
        attr(out, "reference") <- list(scenario = as.integer(reference), counts = d$design_grid[reference,
            , drop = FALSE], value = out$Estimate[reference], point_status = out$PointStatus[reference],
            interval = c(out$Lower[reference], out$Upper[reference]), status = out$Status[reference],
            diagnostics = diagnostics[[reference]])
        lower <- out$Lower - out$Upper[reference]
        upper <- out$Upper - out$Lower[reference]
        feasible <- if (!sem)
            t(vapply(errors, function(b) .mfrm_mvgt_contrast_range(errors[[reference]], b), numeric(2)))
        available <- is.finite(lower) & is.finite(upper)
        status <- out$Status
        if (!all(is.finite(c(out$Lower[reference], out$Upper[reference]))))
            status[] <- paste("Reference interval unavailable:", out$Status[reference])
        linear <- vector("list", length(errors))
        universe <- NULL
        if (any(available & seq_along(errors) != reference)) {
            for (i in which(available)) linear[[i]] <- if (sem)
                diagnostics[[i]]
            else .mfrm_mvgt_ms_region(p$ms, setup, errors[[i]], ratio = FALSE)
            if (!sem)
                universe <- .mfrm_mvgt_ms_region(p$ms, setup, c(1, rep(0, length(p$ms) - 1L)), ratio = FALSE)
        }
        for (i in which(available)) {
            delta_error <- errors[[i]] - errors[[reference]]
            lower_sign <- if (sem)
                all(delta_error >= 0)
            else all(delta_error <= 0)
            upper_sign <- if (sem)
                all(delta_error <= 0)
            else all(delta_error >= 0)
            if (lower_sign)
                lower[i] <- max(0, lower[i])
            if (upper_sign)
                upper[i] <- min(0, upper[i])
            if (!sem) {
                lower[i] <- max(lower[i], feasible[i, 1L])
                upper[i] <- min(upper[i], feasible[i, 2L])
            }
            if (i != reference && any(delta_error != 0)) {
                delta <- .mfrm_mvgt_ms_region(p$ms, setup, if (sem)
                  delta_error
                else -delta_error, ratio = FALSE)
                shared <- .mfrm_mvgt_shared_difference(linear[[reference]]$bounds, linear[[i]]$bounds,
                  delta$bounds, if (!sem)
                    universe$bounds)
                shared$linear_projections <- list(reference = linear[[reference]], comparison = linear[[i]],
                  difference = delta, universe = universe)
                candidate <- shared$bounds * if (sem)
                  sqrt(p$restore)
                else 1
                if (all(is.finite(shared$bounds)) && (any(!is.finite(candidate)) || any(candidate ==
                  0 & shared$bounds != 0) || (diff(shared$bounds) > 0 && diff(candidate) == 0))) {
                  candidate[] <- NA_real_
                  shared$status <- "Unrepresentable difference endpoints"
                }
                shared$bounds <- candidate
                shared$previous_bounds <- c(lower[i], upper[i])
                shared$used <- all(is.finite(candidate))
                if (shared$used) {
                  lower[i] <- max(lower[i], candidate[1L])
                  upper[i] <- min(upper[i], candidate[2L])
                }
                diagnostics[[i]]$shared_difference <- shared
            }
            status[i] <- if (all(delta_error == 0))
                "Identical error weights"
            else "Conservative joint-region enclosure"
            if (lower[i] > upper[i]) {
                lower[i] <- upper[i] <- NA_real_
                status[i] <- "Inconsistent numerical projections"
            }
        }
        out$Reference <- as.integer(reference)
        out$ReferenceValue <- out$Estimate[reference]
        out$Value <- out$Estimate
        out$Estimate <- out$Value - out$ReferenceValue
        out$PointStatus <- ifelse(is.finite(out$Estimate), "Available", "Point difference unavailable")
        out$Metric <- paste(metric, "difference")
        out$Lower <- lower
        out$Upper <- upper
        out$Status <- status
        if (!sem) {
            out$FeasibleLower <- feasible[, 1L]
            out$FeasibleUpper <- feasible[, 2L]
        }
        out$Method <- "Mean-square region difference enclosure"
        out <- out[-reference, , drop = FALSE]
        rownames(out) <- NULL
        attr(out, "diagnostics") <- diagnostics[-reference]
        attr(out, "notes") <- c(attr(out, "notes"), "Differences are comparison minus reference. Individual projections share one simultaneous mean-square region. Their endpoint differences form a conservative enclosure, not the tight image of the difference and not subtraction of independent or marginal confidence intervals.",
            "Shared error-variance differences and bounds on their nonnegative multiplying factor can further narrow this enclosure on the same region. Each comparison's shared_difference diagnostic records the candidate and linear projections; an unavailable refinement retains the preceding enclosure. This is not a tight nonlinear projection or a change in confidence level.",
            "Known componentwise ordering of error weights restricts the sign; identical error weights give [0,0] only when the shared region projections are available. Negative SEM differences and positive G/Phi differences favor the comparison plan.",
            if (!sem) "G/Phi difference enclosures are also intersected with the analytic population range FeasibleLower/Upper implied by the two error-weight vectors and nonnegative variance components. This removes impossible values without adding source information or altering the common-event coverage argument. Original ANOVA points are not clipped to this range.",
            "Under the positive-mean-square normal model the same source-region event protects the returned comparisons jointly for this fixed score/composite. Enclosures can be wide; individual_intervals and reference retain projections and diagnostics. No paired-delta SE is implied.")
    }
    .mfrm_mvgt_interval_result(out, d, p, out$MethodCode)
}

.mfrm_mvgt_shared_difference <- function(a, b, delta, u = NULL) {
    ans <- list(bounds = c(NA_real_, NA_real_), status = "Shared linear projection unavailable")
    if (any(!is.finite(c(a, b, delta, u))))
        return(ans)
    ans$status <- "Unrepresentable shared difference factor"
    scale <- max(a, b, u)
    if (scale <= 0 || any(c(a, b) <= 0) || any(vapply(list(a, b, delta, u), function(z) length(z) > 1L &&
        diff(z) < 0, logical(1))) || (!is.null(u) && any(u < 0)))
        return(ans)
    original <- c(a, b, delta, u)
    a <- a/scale
    b <- b/scale
    delta <- delta/scale
    u <- u/scale
    if (any(c(a, b, delta, u) == 0 & original != 0))
        return(ans)
    if (length(u)) {
        f <- function(x, aa, bb) (x/(x + aa))/(x + bb)
        peak <- max(u[1L], min(u[2L], sqrt(a[1L]) * sqrt(b[1L])))
        factor <- c(min(f(u, a[2L], b[2L])), f(peak, a[1L], b[1L]))
    }
    else factor <- c(1/(sqrt(a[2L]) + sqrt(b[2L])), 1/(sqrt(a[1L]) + sqrt(b[1L])))
    ans$factor <- factor
    ans$scale <- scale
    if (length(u) && any(factor == 0 & u > 0))
        return(ans)
    products <- outer(delta, factor)
    if (any(!is.finite(c(factor, products))) || any(products == 0 & outer(delta != 0, factor != 0, "&")))
        return(ans)
    bounds <- range(products) * if (length(u))
        1
    else sqrt(scale)
    if (any(!is.finite(bounds)) || any(bounds == 0 & range(products) != 0) || (diff(range(products)) >
        0 && diff(bounds) == 0))
        return(ans)
    ans$bounds <- bounds
    ans$status <- "Computed shared-difference enclosure"
    ans
}

.mfrm_mvgt_approximate_f <- function(ms, map, error, level) {
    ans <- list(bounds = c(NA_real_, NA_real_), status = "Invalid F moments", raw_roots = c(NA_real_,
        NA_real_), denominator_df = NA_real_)
    if (max(ms) <= 0)
        return(ans)
    ms <- ms/max(ms)
    theta <- drop(map$A %*% ms)
    rho <- theta[1L]/(theta[1L] + sum(error * theta))
    a <- map$A[1L, ]
    delta <- drop(error %*% map$A)
    b0 <- -a[-1L]/a[1L]
    b1 <- (a[-1L] + delta[-1L])/a[1L]
    wm <- (b0 + rho * b1) * ms[-1L]
    scale <- max(abs(wm))
    if (!is.finite(rho) || !is.finite(scale) || scale <= 0)
        return(ans)
    z <- wm/scale
    total <- sum(z)
    nu <- total^2/sum(z^2/map$df[-1L])
    ans$denominator_df <- nu
    ans$unrestricted_estimate <- rho
    if (!is.finite(nu) || nu <= 0 || total <= 0)
        return(ans)
    tail <- (1 - level)/2
    q <- c(stats::qf(tail, map$df[1L], nu, lower.tail = FALSE), stats::qf(tail, map$df[1L], nu))
    if (any(!is.finite(q) | q <= 0)) {
        ans$status <- "Unrepresentable F quantile"
        return(ans)
    }
    ans$f_quantiles <- q
    ans$quantile_log_gap <- log(q[1L]) - log(q[2L])
    ans$quantile_resolution <- 64 * .Machine$double.eps * max(1, abs(log(q)))
    if (ans$quantile_log_gap <= ans$quantile_resolution) {
        ans$status <- "F quantile interval below numeric resolution"
        return(ans)
    }
    c0 <- sum(b0 * ms[-1L])
    c1 <- sum(b1 * ms[-1L])
    divisor <- ms[1L] + q * c1
    ans$raw_roots <- (ms[1L] - q * c0)/divisor
    inequalities <- rbind(c(c0, c1), c(q[1L] * c0 - ms[1L], divisor[1L]), c(ms[1L] - q[2L] * c0, -divisor[2L]))
    if (any(!is.finite(inequalities))) {
        ans$status <- "Unrepresentable F inequalities"
        return(ans)
    }
    bounds <- c(0, 1)
    for (i in seq_len(nrow(inequalities))) {
        intercept <- inequalities[i, 1L]
        slope <- inequalities[i, 2L]
        if (slope > 0)
            bounds[1L] <- max(bounds[1L], -intercept/slope)
        else if (slope < 0)
            bounds[2L] <- min(bounds[2L], -intercept/slope)
        else if (intercept < 0)
            bounds <- c(1, 0)
    }
    ans$status <- "Empty parameter-space interval"
    if (bounds[1L] <= bounds[2L] && c0 + c1 * mean(bounds) > 0) {
        if (bounds[1L] == bounds[2L] && all(is.finite(ans$raw_roots)) && ans$raw_roots[1L] == ans$raw_roots[2L]) {
            ans$status <- "F interval below numeric resolution"
            return(ans)
        }
        ans$bounds <- bounds
        ans$status <- "Computed approximation"
    }
    ans
}

.mfrm_mvgt_ms_region_setup <- function(map, level) {
    L <- solve(map$A)
    units <- diag(L)
    tail <- -expm1(log(level)/length(map$df))/2
    list(units = units, B = sweep(L, 2L, units, "/"), low = map$df/stats::qchisq(tail, map$df, lower.tail = FALSE),
        high = map$df/stats::qchisq(tail, map$df), tail = tail)
}

.mfrm_mvgt_ms_region <- function(ms, setup, error, ratio = TRUE) {
    ans <- list(bounds = c(NA_real_, NA_real_), status = "Degenerate mean squares", solver_status = c(NA_integer_,
        NA_integer_), residual = c(NA_real_, NA_real_))
    if (any(ms <= 0))
        return(ans)
    ms_scale <- max(ms)
    ms <- ms/ms_scale
    lo <- ms * setup$low
    hi <- ms * setup$high
    ans$mean_square_bounds <- cbind(Lower = lo, Upper = hi)
    ans$tail <- setup$tail
    if (any(!is.finite(c(lo, hi)))) {
        ans$status <- "Unrepresentable chi-square bounds"
        return(ans)
    }
    ans$relative_width <- (hi - lo)/hi
    ans$region_resolution <- 64 * .Machine$double.eps
    if (any(lo <= 0 | ans$relative_width <= ans$region_resolution)) {
        ans$status <- "Mean-square region below numeric resolution"
        return(ans)
    }
    if (ratio) {
        person <- c(1, rep(0, length(ms) - 1L))
        u <- person/setup$units
        v <- (person + error)/setup$units
        scale <- max(v)
        u <- u/scale
        v <- v/scale
        C <- rbind(cbind(setup$B, -hi), cbind(-setup$B, lo), c(v, 0))
        rhs <- c(rep(0, 2L * length(ms)), 1)
        directions <- c(rep("<=", 2L * length(ms)), "=")
        objective <- c(u, 0)
    }
    else {
        objective <- error/setup$units
        scale <- max(abs(objective))
        if (scale == 0)
            scale <- 1
        objective <- objective/scale
        C <- rbind(setup$B, -setup$B)
        rhs <- c(hi, -lo)
        directions <- rep("<=", length(rhs))
    }
    row_scale <- pmax(1, apply(abs(C), 1L, max), abs(rhs))
    C <- C/row_scale
    rhs <- rhs/row_scale
    fits <- lapply(c("min", "max"), function(direction) lpSolve::lp(direction, objective, C, directions,
        rhs))
    ans$solver_status <- vapply(fits, function(f) as.integer(f$status), integer(1))
    if (any(ans$solver_status != 0L)) {
        ans$status <- if (any(ans$solver_status == 2L))
            "Empty mean-square region"
        else "Linear-program solver failure"
        return(ans)
    }
    bounds <- vapply(fits, function(f) f$objval, numeric(1))
    ans$residual <- vapply(fits, function(f) {
        got <- drop(C %*% f$solution)
        max(0, -f$solution, (got - rhs)[directions == "<="], abs((got - rhs)[directions == "="]), abs(sum(objective *
            f$solution) - f$objval))
    }, numeric(1))
    tol <- 1e-07
    nonnegative <- ratio || all(error >= 0)
    if (any(!is.finite(c(bounds, ans$residual))) || any(ans$residual > tol) || (nonnegative && bounds[1L] <
        -tol) || (ratio && bounds[2L] > 1 + tol) || bounds[1L] > bounds[2L] || (ratio && any(vapply(fits,
        function(f) tail(f$solution, 1L) <= 0, logical(1))))) {
        ans$status <- "Linear-program numerical check failed"
        return(ans)
    }
    ans$bounds <- if (ratio)
        pmin(1, pmax(0, bounds))
    else (if (nonnegative)
        pmax(0, bounds)
    else bounds) * scale * ms_scale
    if (any(!is.finite(ans$bounds)) || (!ratio && any(ans$bounds == 0 & bounds != 0 & (!nonnegative |
        bounds > 0)))) {
        ans$bounds[] <- NA_real_
        ans$status <- "Unrepresentable linear projection endpoints"
        return(ans)
    }
    ans$status <- "Computed conservative interval"
    ans
}
