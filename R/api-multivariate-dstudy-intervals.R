confint.mfrm_multivariate_d_study <- function(object, parm = "G", level = 0.95, score = NULL, composite = NULL,
    assumption, draws = 100000L, seed = 1L, reference = NULL, method = "auto", ...) {
    rlang::check_dots_empty()
    if (missing(assumption) || !identical(assumption, "normal"))
        stop("This interval requires `assumption = 'normal'`; the assumption is not tested.", call. = FALSE)
    if (!is.character(parm) || length(parm) != 1L || is.na(parm) || !parm %in% c("G", "Phi", "RelativeSEM",
        "AbsoluteSEM"))
        stop("Select one of 'G', 'Phi', 'RelativeSEM' or 'AbsoluteSEM'.", call. = FALSE)
    sem <- parm %in% c("RelativeSEM", "AbsoluteSEM")
    if (!is.numeric(level) || is.object(level) || is.complex(level) || !is.null(dim(level)) || length(level) !=
        1L || !is.finite(level) || level <= 0 || level >= 1)
        stop("`level` must be strictly between zero and one.", call. = FALSE)
    if (!inherits(object, "mfrm_multivariate_d_study"))
        stop("`object` must be a multivariate D-study result.", call. = FALSE)
    if (!identical(.mfrm_mvds_policy(object), "raw")) {
        d <- .mfrm_mvds_replay(object)
        raw <- mfrm_d_study(d$gstudy, d$design_grid, d$weights)
        ci <- confint.mfrm_multivariate_d_study(raw, parm = parm, level = level, score = score, composite = composite,
            assumption = assumption, draws = draws, seed = seed, reference = reference, method = method)
        return(.mfrm_mvds_adjusted_intervals(ci, d))
    }
    g <- object$gstudy
    d <- mfrm_d_study(g, object$design_grid, object$weights)
    anova <- identical(g$design$method, "Balanced multivariate ANOVA") && isTRUE(g$design$balanced) &&
        isTRUE(g$design$complete)
    nested <- !is.null(g$design$nesting)
    child <- if (nested)
        d$design_grid[[g$design$count_columns[[names(g$design$nesting)]]]]
    facets <- length(g$design$counts) - 1L
    one_target <- facets == 1L || (nested && all(child == g$design$child_counts[1L]))
    if (!is.character(method) || !length(method) || anyNA(method) || anyDuplicated(method) || !is.null(dim(method)))
        stop("`method` must contain distinct method names.", call. = FALSE)
    if (identical(method, "auto")) {
        if (!anova)
            stop("Automatic intervals require a balanced ANOVA source. For a MINQUE source, explicitly select 'delta' for G/Phi, 'reml_profile', 'spectral_f' for one-facet G, or 'residual_contrast' for one-facet RelativeSEM, after reviewing the method's limitations.",
                call. = FALSE)
        if (!sem && !one_target && !is.null(reference))
            stop("Automatic coefficient plan-difference intervals require one random measurement facet or unchanged nested child count per parent. Otherwise explicitly select 'mean_square_region' or 'reml_profile' after reviewing their limitations.",
                call. = FALSE)
        if (!sem && !one_target && nested)
            stop("Automatic nested coefficient intervals require unchanged child count per parent. Otherwise explicitly select 'delta', 'approximate_f', 'mean_square_region' or 'reml_profile' after reviewing their limitations.",
                call. = FALSE)
        method <- if (sem || !one_target)
            "mls"
        else if (parm == "G")
            "exact_f"
        else "gci"
    }
    supported <- if (!anova)
        "reml_profile"
    else if (sem)
        c("mls", if (facets == 2L) "mean_square_region", "reml_profile")
    else c(if (one_target) if (parm == "G") "exact_f" else "gci", if (facets == 2L) c("approximate_f",
        "mean_square_region"), if (facets == 2L && !nested) c("mls", if (parm == "G") "score"), "reml_profile")
    if (facets == 1L && parm == "G")
        supported <- c(supported, "spectral_f")
    if (facets == 1L && parm == "RelativeSEM")
        supported <- c(supported, "residual_contrast")
    if (!sem)
        supported <- c(supported, "delta")
    if (any(!method %in% supported))
        stop("Unsupported method for this design and metric. Available: ", paste(supported, collapse = ", "),
            ".", call. = FALSE)
    if (length(method) > 1L && ((sem && !any(c("reml_profile", "residual_contrast") %in% method)) ||
        (!any(c("reml_profile", "spectral_f", "delta") %in% method) && any(method %in% c("exact_f", "gci")))))
        stop("Select one method, or compare 'approximate_f' and 'mean_square_region' for coefficients.",
            call. = FALSE)
    if (!is.null(reference)) {
        if ("delta" %in% method)
            stop("Delta intervals for differences from a reference plan are not supported here. Omit `reference` for marginal delta intervals.",
                call. = FALSE)
        if (!sem && "mls" %in% method)
            stop("Coefficient MLS `reference` intervals are not supported; marginal null sets do not give a plan-difference interval.",
                call. = FALSE)
        if ("score" %in% method)
            stop("Score `reference` intervals are not supported; marginal score sets do not give a plan-difference interval.",
                call. = FALSE)
        if ("approximate_f" %in% method)
            stop("Approximate-F `reference` intervals are not supported; select 'mean_square_region' for a conservative joint-region difference enclosure.",
                call. = FALSE)
        if (identical(method, "mls") && facets != 1L && (!nested || length(unique(child)) != 1L))
            stop("MLS `reference` comparisons require one facet or nested plans with the same future child count per parent across plans. Otherwise select 'mean_square_region'.",
                call. = FALSE)
    }
    if (!is.null(reference)) {
        if (!is.numeric(reference) || is.object(reference) || is.complex(reference) || !is.null(dim(reference)) ||
            length(reference) != 1L || !is.finite(reference) || reference != floor(reference) || nrow(d$design_grid) <
            2L || reference < 1L || reference > nrow(d$design_grid))
            stop("`reference` must name a valid row in a design grid with at least two plans.", call. = FALSE)
    }
    target <- .mfrm_mvds_interval_target(d, score, composite)
    score <- target$score
    composite <- target$composite
    w <- target$weights
    tab <- target$table
    if (length(method) > 1L && any(c("reml_profile", "residual_contrast", "spectral_f", "delta") %in%
        method))
        return(.mfrm_mvgt_profile_bind(lapply(method, function(code) confint.mfrm_multivariate_d_study(d,
            parm = parm, level = level, score = score, composite = composite, assumption = assumption,
            draws = draws, seed = seed, reference = reference, method = code))))
    if (identical(method, "delta"))
        return(.mfrm_mvgt_delta_table(d, tab, w, parm, level))
    if (identical(method, "spectral_f"))
        return(.mfrm_mvgt_spectral_table(d, tab, w, level, reference))
    if (identical(method, "residual_contrast"))
        return(.mfrm_mvgt_residual_sem(d, tab, w, level, reference))
    if ("reml_profile" %in% method) {
        if (!is.null(reference)) {
            if (facets == 1L || (nested && length(unique(child)) == 1L))
                return(.mfrm_mvgt_profile_difference(d, tab, w, parm, level, reference))
            return(.mfrm_mvgt_contrast_table(d, tab, w, parm, level, reference))
        }
        return(.mfrm_mvgt_profile_table(d, tab, w, parm, level))
    }
    p <- .mfrm_mvgt_anova_projection(g, weights = w)
    if (!sem && "mls" %in% method) {
        if (length(method) > 1L)
            return(.mfrm_mvgt_profile_bind(lapply(method, function(code) confint.mfrm_multivariate_d_study(d,
                parm = parm, level = level, score = score, composite = composite, assumption = assumption,
                method = code))))
        return(.mfrm_mvgt_mls_coefficient_table(d, tab, p, parm, level))
    }
    if ("score" %in% method) {
        if (length(method) > 1L)
            return(.mfrm_mvgt_profile_bind(lapply(method, function(code) confint.mfrm_multivariate_d_study(d,
                parm = parm, level = level, score = score, composite = composite, assumption = assumption,
                draws = draws, seed = seed, method = code))))
        return(.mfrm_mvgt_score_table(d, tab, p, level))
    }
    if (identical(method, "mls"))
        return(.mfrm_mvgt_sem_table(d, tab, p, parm, level, reference))
    if (all(method %in% c("approximate_f", "mean_square_region")))
        return(.mfrm_mvgt_two_facet_intervals(d, tab, p, parm, level, method, reference))
    one <- .mfrm_mvgt_one_facet_target(d, p)
    if (parm == "Phi")
        return(.mfrm_mvgt_interval_result(.mfrm_mvgt_phi_table(d, tab, p, level, draws, seed, reference,
            one), d, p, method))
    interval <- .mfrm_mvgt_g_f_interval(one$ms[1L], one$ms[3L], one$df[c(1L, 3L)], one$observed, one$future,
        level)
    out <- tab[c("Scenario", names(d$design_grid), "Kind", "Score")]
    out$Metric <- "G"
    out$Estimate <- tab$G
    out$PointStatus <- tab$GStatus
    out$Lower <- interval$lower
    out$Upper <- interval$upper
    out$RawRatioLower <- exp(interval$log_bounds[1L])
    out$RawRatioUpper <- exp(interval$log_bounds[2L])
    out$Level <- level
    out$Method <- "Exact F under normal effects"
    out$Status <- interval$projection_status
    rownames(out) <- NULL
    attr(out, "weights") <- p$weights
    attr(out, "assumption") <- assumption
    attr(out, "source_design") <- g$design[c("facets", "counts", "method", "balanced", "complete", "score_convention")]
    attr(out, "source_rows") <- g$data_usage$counts
    attr(out, "pivot") <- list(log_ratio = interval$log_ratio, log_bounds = interval$log_bounds, df = one$df[c(1L,
        3L)])
    attr(out, "notes") <- c("Exact F coverage requires the balanced independent normal-effects model and positive residual variance. It is not distribution robust.",
        "Intervals describe population G for each future count, using uncertainty and degrees of freedom from the source study. They are not Person-rank, pass/fail, Phi or plan-difference intervals.",
        "All counts for this score use the same ratio interval. Do not treat rows as independent or infer a difference from interval overlap.",
        "Intersect the mean-square expectation-ratio interval with R >= 1 before transforming to G. An empty intersection is not zero uncertainty. RawRatioLower/Upper are ratio endpoints, not G endpoints.",
        "PointStatus and Status describe separate outputs. A positive-width ratio interval that rounds to identical G endpoints is unavailable on that scale; retained log-ratio bounds may still resolve plan differences. If the base log-ratio width itself rounds away, neither interval is available. SaveRDS() retains fixed weights, assumptions and source context.")
    if (!is.null(reference)) {
        counts <- one$future
        bounds <- .mfrm_mvgt_g_difference_range(interval, one$observed, counts, counts[reference])
        out$Reference <- as.integer(reference)
        out$ReferenceCount <- counts[reference]
        out$ReferenceValue <- tab$G[reference]
        out$Value <- tab$G
        out$Metric <- "G difference"
        out$Estimate <- out$Value - out$ReferenceValue
        out$PointStatus <- ifelse(is.finite(out$Estimate), "Available", "Point difference unavailable")
        out$Lower <- bounds$lower
        out$Upper <- bounds$upper
        out$Status <- bounds$status
        out$Method <- "F confidence-set image under normal effects"
        out <- out[-reference, , drop = FALSE]
        rownames(out) <- NULL
        attr(out, "reference") <- list(scenario = as.integer(reference), count = counts[reference], value = tab$G[reference],
            point_status = tab$GStatus[reference])
        attr(out, "notes") <- c("G differences map the entire nonnegative F confidence set, including any interior extremum. They do not subtract marginal endpoints or use a delta approximation.",
            "Coverage is at least the requested level under balanced independent normal effects with positive residual variance; it can be conservative. The same event protects the returned comparisons jointly for this one fixed score/composite.",
            "Positive differences favor the comparison plan; the sign from adding levels follows from the model. The interval quantifies size, not cost-effectiveness, pass/fail accuracy or changes in Phi.",
            "RawRatioLower/Upper retain the base ratio interval, not difference endpoints. PointStatus and Status describe separate outputs. An empty base set or unresolved base log-ratio width remains unavailable, including for repeated counts.",
            "SaveRDS() retains the reference, fixed weights, assumptions and source context. Existing D-study plots do not display these difference intervals.")
    }
    .mfrm_mvgt_interval_result(out, d, p, method)
}

.mfrm_mvgt_one_facet_target <- function(d, p) {
    g <- d$gstudy
    if (is.null(g$design$nesting))
        return(list(ms = p$ms, df = p$design$df, observed = unname(g$design$counts[2L]), future = d$design_grid[[1L]]))
    parent <- unname(g$design$nesting)
    list(ms = p$ms[c(1L, 2L, 4L)], df = p$design$df[c(1L, 2L, 4L)], observed = unname(g$design$counts[[parent]]),
        future = d$design_grid[[g$design$count_columns[[parent]]]])
}

.mfrm_mvgt_interval_result <- function(out, d, p, method) {
    out$MethodCode <- rep(method, length.out = nrow(out))
    out$Available <- is.finite(out$Lower) & is.finite(out$Upper)
    out$NegativeComponent <- any(p$estimate < 0)
    if (is.null(out$Validation))
        out$Validation <- ifelse(out$MethodCode == "approximate_f", "Experimental approximate F; finite-sample undercoverage documented",
            ifelse(out$MethodCode == "mean_square_region", "Experimental numerical projection; conservative joint coverage under the normal model",
                ifelse(out$MethodCode == "gci", "Generalized interval; frequentist coverage not exact",
                  "Normal-model F confidence set")))
    attr(out, "design_grid") <- d$design_grid
    attr(out, "source_summary") <- d$gstudy[c("design", "components", "mean_products", "degrees_of_freedom")]
    nesting <- d$gstudy$design$nesting
    if (!is.null(nesting)) {
        g <- d$gstudy
        child <- d$design_grid[[g$design$count_columns[[names(nesting)]]]]
        parent <- d$design_grid[[g$design$count_columns[[unname(nesting)]]]]
        attr(out, "source_design") <- g$design[c("facets", "counts", "nesting", "child_counts", "method",
            "balanced", "complete", "score_convention")]
        attr(out, "nesting_context") <- list(child = names(nesting), parent = unname(nesting), source_children_per_parent = if (length(unique(g$design$child_counts)) ==
            1L) unname(g$design$child_counts[1L]) else g$design$child_counts, future = data.frame(Scenario = seq_along(child),
            ChildrenPerParent = child, Parents = parent))
        if (!is.null(attr(out, "reference"))) {
            ref <- attr(out, "reference")
            if (!is.null(ref$count))
                ref$count_facet <- unname(nesting)
            ref$children_per_parent <- child[ref$scenario]
            attr(out, "reference") <- ref
        }
        note <- if (all(method == "reml_profile") && !is.null(attr(out, "base_interval")))
            "Nested REML comparisons keep a common future child count per parent and vary the parent count. The source may have different or unequal child counts; its full observed assignment remains in the likelihood."
        else if (all(method == "reml_profile"))
            "Nested counts mean children per parent and numbers of parents; both may change. REML uses the five grouping kernels of the observed nested assignment, including interactions, without requiring balanced mean-square pivots."
        else if (all(method == "mls"))
            "Nested counts mean children per parent and numbers of parents. For SEM comparisons, every plan must keep the same future child count; the common base SEM then scales with the inverse square root of the parent count. Exact versus MLS coverage follows the reported Method and Validation, including when the future child count differs from the source."
        else if (all(method %in% c("exact_f", "gci")))
            "The nested child count per parent equals the source count in every plan. Averaging within each Person-by-parent cell gives the one-facet sampling law with effective parent variance v_parent + v_child/R and error variance v_P:parent + v_E/R. This retains both random facets; it does not make a facet fixed. The original source and child counts are retained. ReferenceCount, when present, counts parents."
        else if (all(method == "delta"))
            "Nested counts mean children per parent and numbers of parents; both may change. The delta projection retains the five original ANOVA/MINQUE components and their sampling covariance; no unchanged-child reduction is used."
        else "Nested counts mean children per parent and numbers of parents; both may change. The five-component mean-square map retains the full random nested design. Explicit experimental methods do not use the unchanged-child one-facet reduction."
        attr(out, "notes") <- c(note, attr(out, "notes"))
    }
    class(out) <- c("mfrm_dstudy_intervals", "data.frame")
    out
}

.mfrm_mvgt_g_difference_range <- function(interval, observed, counts, reference) {
    x <- pmax(interval$log_bounds, 0)
    log_t <- x + log(-expm1(-x)) - log(observed)
    .mfrm_mvgt_count_difference_range(log_t, counts, reference, interval$status)
}

.mfrm_mvgt_count_difference_range <- function(log_t, counts, reference, base_status) {
    lower <- upper <- rep(NA_real_, length(counts))
    status <- rep(base_status, length(counts))
    if (base_status %in% c("Available", "Includes zero", "Upper endpoint at zero")) {
        softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
        for (i in seq_along(counts)) {
            b <- counts[i]
            a <- reference
            values <- c(0, 0)
            finite <- is.finite(log_t)
            if (a != b)
                values[finite] <- sign(b - a) * exp(log(abs(b - a)) + log_t[finite] - softplus(log(a) +
                  log_t[finite]) - softplus(log(b) + log_t[finite]))
            lost <- a != b && any(values[finite] == 0)
            stationary <- -(log(a) + log(b))/2
            if (stationary >= log_t[1L] && stationary <= log_t[2L])
                values <- c(values, sign(b - a) * tanh(log1p(abs(b - a)/min(a, b))/4))
            lower[i] <- min(values)
            upper[i] <- max(values)
            status[i] <- if (a == b)
                "Identical counts"
            else if (base_status == "Upper endpoint at zero")
                base_status
            else if (lower[i] <= 0 && upper[i] >= 0)
                "Includes zero"
            else "Available"
            if (lost || (a != b && log_t[1L] < log_t[2L] && lower[i] == upper[i])) {
                lower[i] <- upper[i] <- NA_real_
                status[i] <- "Difference interval below numeric resolution"
            }
        }
    }
    list(lower = lower, upper = upper, status = status)
}

.mfrm_mvgt_phi_table <- function(d, tab, p, level, draws, seed, reference = NULL, one = .mfrm_mvgt_one_facet_target(d,
    p)) {
    integer_scalar <- function(x) is.numeric(x) && !is.object(x) && !is.complex(x) && is.null(dim(x)) &&
        length(x) == 1L && is.finite(x) && x == floor(x)
    tail <- (1 - level)/2
    if (!integer_scalar(draws) || draws < 1000 || draws > .Machine$integer.max || draws * tail < 100)
        stop("`draws` must be an integer >= 1000 with at least 100 expected draws in each tail.", call. = FALSE)
    if (!integer_scalar(seed) || seed < 0 || seed > .Machine$integer.max)
        stop("`seed` must be a nonnegative integer.", call. = FALSE)
    n <- unname(d$gstudy$design$counts[1L])
    k <- one$observed
    m <- one$future
    transform <- function(rho) rho/(rho + (1 - rho)/m)
    bounds <- rep(NA_real_, 2L)
    mc_bounds <- rep(NA_real_, 4L)
    zero <- NA_real_
    status <- "Degenerate mean squares"
    if (all(one$ms > 0)) {
        z <- with_preserved_rng_seed(seed, .mfrm_mvgt_phi_pivots(one$ms, n, k, draws))
        status <- "Numerical interval failure"
        if (all(is.finite(z) & z >= 0 & z <= 1)) {
            z <- sort(z)
            bounds <- unname(stats::quantile(z, c(tail, 1 - tail), type = 7))
            band_tail <- if (is.null(reference))
                0.025
            else 0.0125
            rank <- as.integer(c(stats::qbinom(band_tail, draws, tail), stats::qbinom(1 - band_tail,
                draws, tail) + 1, stats::qbinom(band_tail, draws, 1 - tail), stats::qbinom(1 - band_tail,
                draws, 1 - tail) + 1))
            mc_bounds <- c(0, z, 1)[rank + 1L]
            zero <- mean(z == 0)
            status <- if (bounds[2L] == 0)
                "Upper endpoint at zero"
            else if (bounds[1L] == bounds[2L])
                "Pivotal interval below numeric resolution"
            else if (bounds[1L] == 0)
                "Includes zero"
            else "Available"
        }
    }
    out <- tab[c("Scenario", names(d$design_grid), "Kind", "Score")]
    out$Metric <- "Phi"
    out$Estimate <- tab$Phi
    out$PointStatus <- tab$PhiStatus
    out$Lower <- transform(bounds[1L])
    out$Upper <- transform(bounds[2L])
    out$RawRatioLower <- out$RawRatioUpper <- NA_real_
    out$Level <- level
    out$Method <- "Generalized CI under normal effects"
    out$Status <- status
    collapsed <- is.finite(out$Lower) & is.finite(out$Upper) & out$Lower == out$Upper & (bounds[2L] >
        bounds[1L] | bounds[2L] > 0)
    out$Lower[collapsed] <- out$Upper[collapsed] <- NA_real_
    out$Status[collapsed] <- "Phi interval below numeric resolution"
    rownames(out) <- NULL
    attr(out, "weights") <- p$weights
    attr(out, "assumption") <- "normal"
    attr(out, "source_design") <- d$gstudy$design[c("facets", "counts", "method", "balanced", "complete",
        "score_convention")]
    attr(out, "source_rows") <- d$gstudy$data_usage$counts
    attr(out, "pivot") <- list(single_rating_bounds = bounds, df = one$df)
    attr(out, "monte_carlo") <- list(draws = draws, draws_used = if (all(one$ms > 0)) draws else 0L,
        seed = seed, rng_kind = RNGkind(), quantile_type = 7L, zero_proportion = zero, band_level = 0.95,
        endpoint_bands = data.frame(Scenario = out$Scenario, LowerMCLower = transform(mc_bounds[1L]),
            LowerMCUpper = transform(mc_bounds[2L]), UpperMCLower = transform(mc_bounds[3L]), UpperMCUpper = transform(mc_bounds[4L])))
    attr(out, "monte_carlo")$endpoint_bands[collapsed, -1L] <- NA_real_
    attr(out, "notes") <- c("Phi includes the random facet main effect in absolute error. Independent normal effects and a complete balanced one-facet source are required.",
        "This generalized confidence interval does not guarantee exact frequentist coverage. More pivotal draws reduce numerical error, not model uncertainty or coverage error.",
        "The GCI truncates its pivotal numerator at zero, not the source variance estimates. PointStatus and Status describe separate outputs; a zero endpoint is not proof of zero uncertainty.",
        "The monte_carlo attribute contains separate 95% numerical bands for the endpoints. Increase draws if these bands affect the intended decision; save the whole result with saveRDS().",
        "All future counts transform the same single-rating interval and use source degrees of freedom. Rows are dependent; overlap does not test a plan difference. These are not pass/fail accuracy intervals.",
        "Widths lost to floating-point resolution leave endpoints and their numerical bands unavailable. A base pivotal interval with lost resolution cannot supply difference bounds, including for identical plans. Diagnostic single-rating bounds remain in pivot; the sampled zero atom retains its separate status.")
    if (!is.null(reference)) {
        image <- function(rho, state = status) .mfrm_mvgt_count_difference_range(stats::qlogis(rho),
            m, m[reference], state)
        mapped <- image(bounds)
        outer <- image(mc_bounds[c(1L, 4L)])
        inner <- if (all(is.finite(mc_bounds)) && mc_bounds[2L] <= mc_bounds[3L])
            image(mc_bounds[c(2L, 3L)])
        else list(lower = outer$upper, upper = outer$lower)
        out$Reference <- as.integer(reference)
        out$ReferenceCount <- m[reference]
        out$ReferenceValue <- tab$Phi[reference]
        out$Value <- tab$Phi
        out$Metric <- "Phi difference"
        out$Estimate <- out$Value - out$ReferenceValue
        out$PointStatus <- ifelse(is.finite(out$Estimate), "Available", "Point difference unavailable")
        out$Lower <- mapped$lower
        out$Upper <- mapped$upper
        out$Status <- mapped$status
        out$Method <- "Generalized confidence-set image under normal effects"
        out <- out[-reference, , drop = FALSE]
        rownames(out) <- NULL
        attr(out, "reference") <- list(scenario = as.integer(reference), count = m[reference], value = tab$Phi[reference],
            point_status = tab$PhiStatus[reference])
        mc <- attr(out, "monte_carlo")
        mc$base_status <- status
        mc$base_endpoint_bands <- matrix(mc_bounds, nrow = 2L, byrow = TRUE, dimnames = list(c("Lower",
            "Upper"), c("MCLower", "MCUpper")))
        mc$band_scope <- "At least 95% jointly over both numerical endpoints and returned comparisons; not coefficient coverage"
        mc$endpoint_bands <- data.frame(Scenario = out$Scenario, LowerMCLower = outer$lower[-reference],
            LowerMCUpper = inner$lower[-reference], UpperMCLower = inner$upper[-reference], UpperMCUpper = outer$upper[-reference])
        attr(out, "monte_carlo") <- mc
        attr(out, "notes") <- c("Phi differences map the whole single-rating generalized interval, including any interior extremum. They do not subtract marginal endpoints or take quantiles of transformed differences.",
            "Whenever the base generalized set contains the true ratio, all mapped comparisons contain their true differences. Approximate frequentist coverage of the base method is retained; there is no exact nominal-coverage guarantee.",
            "Positive differences favor the comparison plan. Adding levels already determines the direction under this model; the interval describes size, not cost-effectiveness or pass/fail accuracy.",
            "The monte_carlo attribute gives at least 95% joint numerical bands from Bonferroni-protected base endpoints and their full set images. These describe quantile error, not sampling coverage. More draws do not add source information.",
            "PointStatus and Status remain separate. The base pivot bounds, zero proportion and reference context are retained; a sampled upper endpoint at zero is not proof of zero uncertainty. Base widths lost to numeric resolution leave all difference intervals unavailable. Save the full result with saveRDS().")
    }
    out
}

.mfrm_mvgt_phi_pivots <- function(ms, n, k, draws) {
    facet <- (k - 1) * ms[2L]/stats::rchisq(draws, k - 1)
    person <- (n - 1) * ms[1L]/stats::rchisq(draws, n - 1)
    residual <- (n - 1) * (k - 1) * ms[3L]/stats::rchisq(draws, (n - 1) * (k - 1))
    pmax(0, person - residual)/(person + (k/n) * facet + (k - 1 - k/n) * residual)
}

.mfrm_mvgt_g_f_interval <- function(person_ms, residual_ms, df, observed, future, level) {
    lower <- upper <- rep(NA_real_, length(future))
    log_bounds <- rep(NA_real_, 2L)
    log_ratio <- NA_real_
    status <- "Zero residual mean square"
    if (residual_ms > 0) {
        log_ratio <- log(person_ms) - log(residual_ms)
        tail <- (1 - level)/2
        log_bounds <- log_ratio - log(c(stats::qf(tail, df[1L], df[2L], lower.tail = FALSE), stats::qf(tail,
            df[1L], df[2L])))
        status <- "Numerical interval failure"
        if (!anyNA(log_bounds) && log_bounds[1L] <= log_bounds[2L]) {
            status <- "Empty parameter-space interval"
            if (is.finite(log_bounds[2L]) && log_bounds[2L] >= 0 && log_bounds[1L] == log_bounds[2L]) {
                status <- "F ratio interval below numeric resolution"
            }
            else if (log_bounds[2L] >= 0) {
                bounded <- pmax(log_bounds, 0)
                transform <- function(x) {
                  u <- -expm1(-x)
                  u/(u + (observed/future) * exp(-x))
                }
                lower <- transform(bounded[1L])
                upper <- transform(bounded[2L])
                status <- if (bounded[1L] == 0)
                  "Includes zero"
                else "Available"
            }
        }
    }
    projection_status <- rep(status, length(future))
    collapsed <- is.finite(lower) & is.finite(upper) & lower == upper & log_bounds[2L] > max(0, log_bounds[1L])
    lower[collapsed] <- upper[collapsed] <- NA_real_
    projection_status[collapsed] <- "G interval below numeric resolution"
    list(lower = lower, upper = upper, log_ratio = log_ratio, log_bounds = log_bounds, status = status,
        projection_status = projection_status)
}
