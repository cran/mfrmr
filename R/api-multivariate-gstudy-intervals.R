.mfrm_gstudy_confint_moments <- function(object, parm = NULL, level = 0.95, score = NULL, weights = NULL,
    assumption, method = "auto", draws = 20000L, seed = 20261016L, ...) {
    rlang::check_dots_empty()
    if (missing(assumption) || !identical(assumption, "normal"))
        stop("This method requires `assumption = 'normal'`; it does not test the random-effects assumption.",
            call. = FALSE)
    if (!is.numeric(level) || is.complex(level) || is.object(level) || !is.null(dim(level)) || length(level) !=
        1L || !is.finite(level) || level <= 0 || level >= 1)
        stop("`level` must be a number strictly between zero and one.", call. = FALSE)
    if (!is.character(method) || length(method) != 1L || is.na(method))
        stop("`method` must name one interval method.", call. = FALSE)
    if (method == "residual_contrast") {
        if (!missing(draws) || !missing(seed))
            stop("`draws` and `seed` apply only to `method = 'gci'`.", call. = FALSE)
        return(.mfrm_gstudy_interval_context(.mfrm_mvgt_residual_intervals(object, parm, score, weights,
            level), object))
    }
    if (is.character(score) && length(score) == 2L && is.null(weights)) {
        if (!method %in% c("gci", "regression_region", "reml_profile"))
            stop("For two-score covariance intervals, explicitly select `method = 'gci', 'regression_region', 'reml_profile' or residual-only 'residual_contrast'.",
                call. = FALSE)
        if (method != "gci" && (!missing(draws) || !missing(seed)))
            stop("`draws` and `seed` apply only to `method = 'gci'`.", call. = FALSE)
        if (method == "reml_profile")
            return(.mfrm_gstudy_interval_context(.mfrm_mvgt_covariance_profile(object, parm, score, level),
                object))
        out <- .mfrm_mvgt_covariance_intervals(object, score, level, method, draws, seed)
        return(.mfrm_gstudy_interval_context(out[.mfrm_mvgt_component_selection(parm, out$Source), ,
            drop = FALSE], object))
    }
    if (!method %in% c("auto", "mls", "reml_profile") || !missing(draws) || !missing(seed))
        stop("Single-score or composite component intervals use `method = 'auto', 'mls', 'reml_profile' or residual-only 'residual_contrast', without draws or seed.",
            call. = FALSE)
    if (method == "reml_profile")
        return(.mfrm_gstudy_interval_context(.mfrm_mvgt_component_profile(object, parm, score, weights,
            level), object))
    projection <- .mfrm_mvgt_anova_projection(object, score, weights)
    design <- projection$design
    sources <- rownames(design$A)
    weights <- projection$weights
    score <- projection$score
    kind <- projection$kind
    ms <- projection$ms
    estimate <- projection$estimate
    restore <- projection$restore
    selected <- .mfrm_mvgt_component_selection(parm, sources)
    intervals <- lapply(seq_along(sources), function(i) .mfrm_mvgt_linear_variance_interval(ms, design$df,
        design$A[i, ], level))
    se <- vapply(intervals, `[[`, numeric(1), "se") * restore
    raw <- t(vapply(intervals, `[[`, numeric(2), "raw")) * restore
    status <- vapply(intervals, `[[`, character(1), "status")
    for (i in which(status %in% c("Available", "Includes zero"))) {
        original <- intervals[[i]]$raw
        if (any(!is.finite(raw[i, ])) || any(raw[i, ] == 0 & original > 0)) {
            status[i] <- "Unrepresentable component endpoints"
        }
        else if (raw[i, 1L] == raw[i, 2L] && original[1L] < original[2L]) {
            status[i] <- "Component interval below numeric resolution"
        }
    }
    estimate <- estimate * restore
    if (any(!is.finite(c(estimate, se))))
        stop("Estimates or SEs exceed numeric range; rescale the scores or weights.", call. = FALSE)
    lower <- pmax(raw[, 1L], 0)
    upper <- raw[, 2L]
    unavailable <- !status %in% c("Available", "Includes zero")
    lower[unavailable] <- upper[unavailable] <- NA_real_
    out <- data.frame(Source = sources, Score = score, Estimate = unname(estimate), SE = unname(se),
        Lower = lower, Upper = upper, RawLower = raw[, 1L], RawUpper = raw[, 2L], Level = level, Method = vapply(intervals,
            `[[`, character(1), "method"), Status = status)
    out <- out[selected, , drop = FALSE]
    rownames(out) <- NULL
    attr(out, "weights") <- weights
    attr(out, "kind") <- kind
    attr(out, "assumption") <- assumption
    attr(out, "source_design") <- object$design[c("facets", "counts", "nesting", "child_counts", "balanced",
        "complete", "method", "score_convention")]
    attr(out, "source_rows") <- object$data_usage$counts
    attr(out, "notes") <- c("Pointwise normal-random-effects intervals for a prespecified score or fixed composite; MLS intervals are approximate. Not simultaneous or distribution robust.",
        "Estimates and SEs have squared score units. SE is sampling uncertainty in a component, not measurement SEM or uncertainty in G/Phi.",
        "Raw estimates and endpoints are retained. Lower/Upper intersect the interval with nonnegative variances; consult Status for empty or unavailable intervals.",
        "Positive intervals whose width or endpoints cannot be represented remain unavailable with an explicit Status; a rounded zero width is not exact population knowledge.",
        "SE uses the df + 2 moment correction. Weights are used as supplied; saveRDS() preserves assumptions and source context.")
    .mfrm_gstudy_interval_context(out, object)
}

.mfrm_mvgt_anova_projection <- function(object, score = NULL, weights = NULL) {
    design <- .mfrm_mvgt_anova_summary(object)
    sources <- rownames(design$A)
    scores <- object$design$scores
    selection <- .mfrm_mvgt_interval_score(scores, score, weights)
    weights <- selection$weights
    score <- selection$score
    kind <- selection$kind
    wscale <- max(abs(weights))
    w <- weights/wscale
    active <- which(w != 0)
    mscale <- max(abs(unlist(lapply(object$mean_products, function(m) m[active, active, drop = FALSE]))))
    unit <- if (mscale == 0)
        1
    else mscale
    restore <- unit * wscale * wscale
    if (!is.finite(restore) || restore == 0)
        stop("Projected variances exceed numeric range; rescale the scores or weights.", call. = FALSE)
    project <- function(m) as.numeric(crossprod(w[active], (m[active, active, drop = FALSE]/unit) %*%
        w[active]))
    ms <- vapply(object$mean_products, project, numeric(1))
    estimate <- vapply(object$components, project, numeric(1))
    tol <- 64 * .Machine$double.eps * length(scores)^2
    if (any(!is.finite(c(ms, estimate))) || any(ms < -tol))
        stop("Stored mean products give invalid mean squares for the selected score.", call. = FALSE)
    ms <- pmax(ms, 0)
    mapped <- drop(design$A %*% ms)
    if (any(abs(estimate - mapped) > tol * pmax(1, rowSums(abs(design$A)))))
        stop("Stored components disagree with their ANOVA mean squares; recreate the G-study.", call. = FALSE)
    list(design = design, weights = weights, score = score, kind = kind, ms = ms, estimate = estimate,
        mapped = mapped, restore = restore)
}

.mfrm_mvgt_interval_score <- function(scores, score, weights) {
    if (!is.null(score) && !is.null(weights))
        stop("Choose either `score` or `weights`, not both.", call. = FALSE)
    kind <- if (is.null(weights))
        "Score"
    else "Composite"
    if (is.null(weights)) {
        if (is.null(score))
            score <- scores[1L]
        if (!is.character(score) || length(score) != 1L || is.na(score) || !score %in% scores)
            stop("`score` must name one original score.", call. = FALSE)
        weights <- setNames(as.numeric(scores == score), scores)
    }
    else {
        if (!is.numeric(weights) || is.object(weights) || is.complex(weights) || !is.null(dim(weights)) ||
            length(weights) != length(scores) || is.null(names(weights)) || anyNA(names(weights)) ||
            anyDuplicated(names(weights)) || !setequal(names(weights), scores) || any(!is.finite(weights)) ||
            !any(weights != 0))
            stop("`weights` must name every score once, with finite values and at least one nonzero weight.",
                call. = FALSE)
        weights <- weights[scores]
        score <- "Composite"
    }
    list(weights = weights, score = score, kind = kind)
}

.mfrm_mvgt_component_selection <- function(parm, sources) {
    selected <- seq_along(sources)
    if (!is.null(parm)) {
        if (is.character(parm))
            selected <- match(parm, sources)
        else {
            if (!is.numeric(parm) || is.complex(parm) || is.object(parm) || !is.null(dim(parm)) || any(!is.finite(parm)) ||
                any(parm != floor(parm)))
                stop("`parm` must contain component names or integer positions.", call. = FALSE)
            selected <- parm
        }
        if (!length(selected) || anyNA(selected) || anyDuplicated(selected) || any(selected < 1 | selected >
            length(sources)))
            stop("`parm` must select distinct existing components.", call. = FALSE)
    }
    selected
}

.mfrm_mvgt_anova_summary <- function(object) {
    if (!inherits(object, "mfrm_gtheory") || !identical(object$design$method, "Balanced multivariate ANOVA") ||
        !isTRUE(object$design$balanced) || !isTRUE(object$design$complete))
        stop("MLS component intervals require a complete balanced ANOVA G-study; for MINQUE(0), explicitly select 'reml_profile' or residual-only 'residual_contrast' after reviewing its limitations.",
            call. = FALSE)
    design <- .mfrm_mvgt_anova_interval_design(object$design)
    sources <- rownames(design$A)
    scores <- object$design$scores
    valid_matrices <- function(m) is.list(m) && identical(names(m), sources) && all(vapply(m, function(a) is.matrix(a) &&
        is.numeric(a) && !is.complex(a) && identical(dimnames(a), list(scores, scores)) && all(is.finite(a)) &&
        identical(a, t(a)), logical(1)))
    if (!is.character(scores) || !length(scores) || anyNA(scores) || any(!nzchar(scores)) || anyDuplicated(scores) ||
        !valid_matrices(object$components) || !valid_matrices(object$mean_products) || !identical(object$degrees_of_freedom$Source,
        sources) || !isTRUE(all.equal(object$degrees_of_freedom$DF, unname(design$df), tolerance = 0)))
        stop("Stored ANOVA matrices or degrees of freedom are incomplete or inconsistent with the design; recreate the G-study.",
            call. = FALSE)
    design
}

.mfrm_mvgt_anova_interval_design <- function(design) {
    n <- design$counts
    factors <- names(n)
    if (!is.numeric(n) || is.complex(n) || !is.null(dim(n)) || !length(n) || any(!is.finite(n)) || any(n <
        2 | n != floor(n)) || !identical(factors[1L], "Person") || !.mfrm_mvgt_valid_facets(factors[-1L]))
        stop("The G-study factor counts or identities are incomplete or altered.", call. = FALSE)
    nesting <- .mfrm_mvgt_nesting(design$nesting, factors[-1L])
    sources <- .mfrm_mvgt_sources(factors, nesting)
    if (!is.null(nesting)) {
        P <- n[1L]
        T <- n[unname(nesting)]
        children <- design$child_counts
        if (!is.numeric(children) || length(children) != T || any(!is.finite(children)) || any(children <
            2 | children != floor(children)) || length(unique(children)) != 1L || sum(children) != n[names(nesting)])
            stop("The balanced nested child counts are incomplete or inconsistent.", call. = FALSE)
        R <- children[1L]
        A <- rbind(c(1, 0, 0, -1, 0)/(R * T), c(0, 1, -1, -1, 1)/(P * R), c(0, 0, 1, 0, -1)/P, c(0, 0,
            0, 1, -1)/R, c(0, 0, 0, 0, 1))
        df <- c(P - 1, T - 1, T * (R - 1), (P - 1) * (T - 1), T * (P - 1) * (R - 1))
    }
    else {
        subsets <- if (length(n) == 2L)
            list(1L, 2L, 1:2)
        else list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
        L <- matrix(0, length(subsets), length(subsets))
        for (i in seq_along(subsets)) for (j in seq_along(subsets)) {
            if (all(subsets[[i]] %in% subsets[[j]]))
                L[i, j] <- prod(n[setdiff(seq_along(n), subsets[[j]])])
        }
        A <- solve(L)
        df <- vapply(subsets, function(s) prod(n[s] - 1), numeric(1))
    }
    dimnames(A) <- list(sources, sources)
    list(A = A, df = unname(df))
}

.mfrm_mvgt_mls_width <- function(ms, df, a, tail) {
    p <- which(a > 0)
    n <- which(a < 0)
    z <- abs(a) * ms
    scale <- max(z)
    if (scale == 0)
        return(0)
    z <- z/scale
    G <- 1 - df/stats::qchisq(tail, df, lower.tail = FALSE)
    H <- df/stats::qchisq(tail, df) - 1
    terms <- c((G[p] * z[p])^2, (H[n] * z[n])^2)
    for (i in p) for (j in n) {
        f <- stats::qf(tail, df[i], df[j], lower.tail = FALSE)
        cross <- f - 2 + 1/f - G[i]^2 * f - H[j]^2/f
        terms <- c(terms, cross * z[i] * z[j])
    }
    if (length(p) > 1L && length(n)) {
        for (j in seq.int(2L, length(p))) for (i in seq_len(j - 1L)) {
            s <- p[i]
            t <- p[j]
            total <- df[s] + df[t]
            joint <- 1 - total/stats::qchisq(tail, total, lower.tail = FALSE)
            cross <- (joint^2 * total^2/(df[s] * df[t]) - G[s]^2 * df[s]/df[t] - G[t]^2 * df[t]/df[s])/(length(p) -
                1L)
            terms <- c(terms, cross * z[s] * z[t])
        }
    }
    value <- sum(terms)
    if (!is.finite(value) || value < -64 * .Machine$double.eps * sum(abs(terms)))
        return(NA_real_)
    sqrt(max(0, value)) * scale
}
