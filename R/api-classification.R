classification_table <- function(fit, rating_plan = NULL, cut_scores = 0, score_cuts = NULL, theta = seq(-3,
    3, by = 0.5), theta_weights = NULL, operational_cuts = cut_scores, facets = NULL, readiness_policy = c("error",
    "review")) {
    .mfrm_validate_fit(fit)
    if (!fit$model %in% c("RSM", "PCM") || !fit$method %in% c("MML", "JML") || identical(fit$control$method_code,
        1L) || (!is.null(fit$interactions) && nrow(fit$interactions) > 0L)) {
        stop("`fit` must be a native ordinary RSM/PCM MML/JML calibration without interactions or random effects.",
            call. = FALSE)
    }
    readiness_policy <- match.arg(readiness_policy)
    cut_scores <- mfrm_validate_cut_scores(cut_scores)
    if (is.null(cut_scores))
        stop("Supply at least one true ability cut.", call. = FALSE)
    operational_cuts <- mfrm_validate_cut_scores(operational_cuts)
    if (length(operational_cuts) != length(cut_scores)) {
        stop("Supply one `operational_cuts` value per true ability cut.", call. = FALSE)
    }
    finite_vector <- function(x) is.numeric(x) && !is.complex(x) && !is.object(x) && is.null(dim(x)) && length(x) >
        0L && all(is.finite(x))
    if (!finite_vector(theta))
        stop("`theta` must be a nonempty finite numeric vector.", call. = FALSE)
    if (!is.null(theta_weights) && (!finite_vector(theta_weights) || length(theta_weights) != length(theta) ||
        any(theta_weights < 0) || !any(theta_weights > 0))) {
        stop("`theta_weights` must match theta and be finite, nonnegative, with positive total mass.", call. = FALSE)
    }
    enc <- fit$encoding
    base_f <- enc$base_facets
    if (is.null(facets)) {
        mapping <- stats::setNames(base_f, base_f)
    }
    else {
        if (anyNA(facets) || anyNA(names(facets)))
            stop("`facets` must not contain missing names.", call. = FALSE)
        if (is.null(names(facets))) {
            if (length(facets) != length(base_f))
                stop("Unnamed `facets` must match model facets length.", call. = FALSE)
            mapping <- stats::setNames(as.character(facets), base_f)
        }
        else {
            if (!setequal(names(facets), base_f))
                stop("Named `facets` must match model facets: ", paste(base_f, collapse = ", "), ".", call. = FALSE)
            mapping <- as.character(facets[base_f])
            names(mapping) <- base_f
        }
    }
    if (is.null(rating_plan)) {
        rating_plan <- unique(enc$retained_data[, base_f, drop = FALSE])
        rownames(rating_plan) <- NULL
        names(rating_plan) <- unname(mapping)
    }
    if (!is.data.frame(rating_plan) || !nrow(rating_plan) || anyDuplicated(names(rating_plan)) || anyDuplicated(unname(mapping)) ||
        !setequal(names(rating_plan), unname(mapping))) {
        stop("`rating_plan` must contain exactly the mapped facet columns and at least one rating.", call. = FALSE)
    }
    plan <- as.data.frame(rating_plan[unname(mapping)], stringsAsFactors = FALSE)
    names(plan) <- names(mapping)
    for (fi in seq_along(base_f)) {
        f <- base_f[fi]
        ids <- plan[[f]]
        if (!(is.character(ids) || is.factor(ids) || is.numeric(ids)) || !is.null(dim(ids)) || anyNA(ids) ||
            any(!nzchar(trimws(as.character(ids))))) {
            stop("Rating-plan identifiers must be nonmissing, nonblank vectors.", call. = FALSE)
        }
        plan[[f]] <- as.character(ids)
        if (any(!plan[[f]] %in% enc$facet_levels[[fi]])) {
            stop("Unknown rating-plan level in facet `", f, "`.", call. = FALSE)
        }
    }
    if (anyDuplicated(plan))
        stop("Duplicate complete facet cells in `rating_plan` are not supported.", call. = FALSE)
    max_sum <- nrow(plan) * (enc$n_cat - 1L)
    if (is.null(score_cuts)) {
        score_cuts <- pmin(max_sum, pmax(1L, round(seq(1, length(cut_scores))/(length(cut_scores) + 1) * max_sum)))
    }
    if (!finite_vector(score_cuts) || length(score_cuts) != length(cut_scores) || any(score_cuts != floor(score_cuts)) ||
        any(diff(score_cuts) < 0) || any(score_cuts < 0 | score_cuts > max_sum + 1)) {
        stop("`score_cuts` must be nondecreasing integers from 0 through maximum sum + 1, one per true cut.",
            call. = FALSE)
    }
    source <- .mfrm_readiness(fit)
    if (identical(fit$method, "JML") && any(fit$person$ExtremeFlag != 0L)) {
        source$reason_codes <- c(source$reason_codes, "boundary_not_scoring_ready:has_exclusions")
        source$ready <- FALSE
        source$status <- "review_required"
    }
    if (!isTRUE(source$ready) && readiness_policy == "error") {
        stop("The source calibration is not ready: ", paste(source$reason_codes, collapse = ", "), ". Use readiness_policy = 'review' only for explicit review.",
            call. = FALSE)
    }
    offsets <- numeric(nrow(plan))
    for (fi in seq_along(base_f)) {
        fnm <- base_f[fi]
        f_sub <- fit$facets[fit$facets$Facet == fnm, ]
        b_vec <- f_sub$Estimate[match(plan[[fnm]], f_sub$Level)]
        offsets <- offsets + enc$facet_signs[fi] * b_vec
    }
    sg_idx <- if (identical(fit$model, "PCM"))
        match(plan[[enc$step_facet]], enc$step_owners)
    else rep(1L, nrow(plan))
    cumulative <- fit$step_cumulative[sg_idx, , drop = FALSE]
    score_map <- data.frame(ObservedScore = enc$rating_min:enc$rating_max, InternalScore = enc$rating_min:enc$rating_max,
        CategoryScore = seq.int(0L, enc$n_cat - 1L))
    cfg_model <- fit$model
    cfg_method <- fit$method
    n <- length(theta)
    categories <- seq_len(length(cut_scores) + 1L)
    true <- findInterval(theta, cut_scores) + 1L
    methods <- c("lee_sum", "rudner_normal")
    tables <- probabilities <- scores <- vector("list", n)
    for (i in seq_len(n)) {
        point <- mfrm_classification_point(theta[i], offsets, cumulative, score_cuts, operational_cuts)
        p <- rbind(point$lee, point$rudner)
        status <- c("available", if (all(is.finite(point$rudner))) "available" else "zero_information")
        tables[[i]] <- data.frame(Point = i, Theta = theta[i], TrueCategory = true[i], Method = methods, Accuracy = p[,
            true[i]], Consistency = rowSums(p^2), Overclassification = rowSums(p[, categories > true[i], drop = FALSE]),
            Underclassification = rowSums(p[, categories < true[i], drop = FALSE]), ExpectedSum = point$expected,
            Information = point$information, SEM = point$sem, Status = status, row.names = NULL)
        tables[[i]][status != "available", c("Accuracy", "Consistency", "Overclassification", "Underclassification")] <- NA_real_
        probabilities[[i]] <- data.frame(Point = i, Theta = theta[i], Method = rep(methods, each = length(categories)),
            Category = rep(categories, 2), Probability = as.vector(t(p)), Status = rep(status, each = length(categories)))
        scores[[i]] <- data.frame(Point = i, Theta = theta[i], SumScore = seq_along(point$pmf) - 1L, Probability = point$pmf)
    }
    table <- do.call(rbind, tables)
    probabilities <- do.call(rbind, probabilities)
    cut_moments <- function(x) t(vapply(x, function(t) {
        z <- mfrm_classification_response(t, offsets, cumulative)
        c(z$expected, z$information, z$sem)
    }, numeric(3)))
    at_true <- cut_moments(cut_scores)
    at_operational <- cut_moments(operational_cuts)
    cuts <- data.frame(Cut = names(cut_scores) %||% paste0("Cut", seq_along(cut_scores)), TrueCut = unname(cut_scores),
        OperationalAbilityCut = unname(operational_cuts), ScoreCut = unname(score_cuts), ExpectedSumAtTrueCut = at_true[,
            1], InformationAtTrueCut = at_true[, 2], SEMAtTrueCut = at_true[, 3], ExpectedSumAtOperationalCut = at_operational[,
            1], InformationAtOperationalCut = at_operational[, 2], SEMAtOperationalCut = at_operational[,
            3])
    population <- data.frame(Point = seq_len(n), Theta = unname(theta), Weight = if (is.null(theta_weights))
        NA_real_
    else unname(theta_weights), NormalizedWeight = NA_real_)
    marginal <- table[FALSE, c("Method", "Accuracy", "Consistency", "Overclassification", "Underclassification",
        "Status")]
    joint <- data.frame(Method = character(), TrueCategory = integer(), ObservedCategory = integer(), Probability = numeric(),
        Status = character())
    if (!is.null(theta_weights)) {
        w <- theta_weights/max(theta_weights)
        w <- w/sum(w)
        population$NormalizedWeight <- unname(w)
        marginal <- do.call(rbind, lapply(methods, function(method) {
            d <- table[table$Method == method & theta_weights[table$Point] > 0, , drop = FALSE]
            ok <- all(d$Status == "available")
            out <- data.frame(Method = method)
            for (col in c("Accuracy", "Consistency", "Overclassification", "Underclassification")) out[[col]] <- if (ok)
                sum(d[[col]] * w[d$Point])
            else NA_real_
            out$Status <- if (ok)
                "available"
            else "positive_weight_point_unavailable"
            out
        }))
        joint <- do.call(rbind, lapply(methods, function(method) {
            d <- probabilities[probabilities$Method == method & theta_weights[probabilities$Point] > 0, ]
            out <- expand.grid(TrueCategory = categories, ObservedCategory = categories)
            ok <- all(d$Status == "available")
            out$Probability <- vapply(seq_len(nrow(out)), function(j) {
                if (!ok)
                  return(NA_real_)
                take <- true[d$Point] == out$TrueCategory[j] & d$Category == out$ObservedCategory[j]
                sum(d$Probability[take] * w[d$Point[take]])
            }, 0)
            data.frame(Method = method, out, Status = if (ok)
                "available"
            else "positive_weight_point_unavailable")
        }))
    }
    table$SourceStatus <- rep(source$status, nrow(table))
    probabilities$SourceStatus <- rep(source$status, nrow(probabilities))
    scores <- do.call(rbind, scores)
    scores$SourceStatus <- rep(source$status, nrow(scores))
    cuts$SourceStatus <- rep(source$status, nrow(cuts))
    marginal$SourceStatus <- rep(source$status, nrow(marginal))
    joint$SourceStatus <- rep(source$status, nrow(joint))
    list(table = table, probabilities = probabilities, scores = scores, cuts = cuts, marginal = marginal,
        classification_matrix = joint, population = population, rating_plan = plan, calibration = list(model = cfg_model,
            method = cfg_method, offsets = offsets, cumulative_steps = cumulative, score_map = score_map),
        source = source, settings = list(facets = mapping, readiness_policy = readiness_policy, category_rule = "lower_inclusive_upper_exclusive",
            calibration = "fixed_point_estimates", repetition = "conditionally_independent_identical_plan",
            maximum_sum = max_sum), notes = c("Lee uses declared category-sum cuts; Rudner approximates ML ability decisions using declared ability cuts.",
            "SEM is inverse-information precision on the ability scale, conditional on fixed calibration and local independence.",
            "Source review status applies to all rows; numerical availability does not establish source or approximation adequacy.",
            "The classification matrix contains joint population probabilities; its true-category row sums are the population proportions, not necessarily one.",
            "Marginal overclassification and underclassification average over the full supplied population. For two categories, divide each joint error probability by its true-category population proportion for a conditional false-pass or false-fail rate; a zero denominator is undefined.",
            "Conditioning on an observed decision instead uses that observed category's column sum in the joint matrix and answers a different question.",
            "Population indices require supplied true-ability weights; calibration, cut and population estimation uncertainty is excluded."))
}

mfrm_classification_response <- function(theta, offsets, cumulative) {
    p <- mfrm_jml_probability_bundle(theta + offsets, rep(0L, length(offsets)), "PCM", cumulative, criterion_idx = seq_along(offsets))$probs
    if (any(!is.finite(p)) || any(p < 0) || any(abs(rowSums(p) - 1) > 1e-10))
        stop("Classification response probabilities are numerically unresolved.", call. = FALSE)
    k <- seq_len(ncol(p)) - 1L
    expected <- as.vector(p %*% k)
    centered <- matrix(k, nrow(p), length(k), byrow = TRUE) - expected
    information <- sum(p * centered^2)
    list(probs = p, expected = sum(expected), information = information, sem = 1/sqrt(information))
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
    pmf <- pmf/sum(pmf)
    bands <- findInterval(seq_along(pmf) - 1L, score_cuts) + 1L
    lee <- vapply(seq_len(length(score_cuts) + 1L), function(h) sum(pmf[bands == h]), 0)
    rudner <- rep(NA_real_, length(lee))
    if (is.finite(x$sem)) {
        cuts <- c(-Inf, operational_cuts, Inf)
        edges <- (cuts - theta)/x$sem
        rudner <- vapply(seq_along(lee), function(h) {
            lo <- edges[h]
            hi <- edges[h + 1L]
            tails <- if (lo >= 0)
                stats::pnorm(c(lo, hi), lower.tail = FALSE)
            else stats::pnorm(c(hi, lo))
            value <- tails[1L] - tails[2L]
            width <- (cuts[h + 1L] - cuts[h])/x$sem
            if (is.finite(lo) && is.finite(width) && width > 0 && value < 1e-04 * tails[1L]) {
                peak <- stats::dnorm(max(lo, min(0, lo + width)), log = TRUE)
                if (is.finite(peak)) {
                  mass <- stats::integrate(function(t) exp(stats::dnorm(lo + t * width, log = TRUE) -
                    peak), 0, 1, rel.tol = 1e-10, abs.tol = 0)$value
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
