.mfrm_gstudy_vcov_moments <- function(object, score = NULL, weights = NULL, assumption, ...) {
    rlang::check_dots_empty()
    if (missing(assumption) || !identical(assumption, "normal"))
        stop("Sampling covariance requires `assumption = 'normal'`; the assumption is not tested.", call. = FALSE)
    g <- object
    scores <- g$design$scores
    sources <- names(g$components)
    if (!inherits(g, "mfrm_gtheory") || !is.character(scores) || !length(scores) || anyNA(scores) ||
        anyDuplicated(scores) || !.mfrm_mvgt_valid_facets(names(g$design$facets)) || !identical(sources,
        .mfrm_mvgt_sources(c("Person", names(g$design$facets)), g$design$nesting)) || !all(vapply(g$components,
        function(m) is.matrix(m) && is.numeric(m) && !is.complex(m) && identical(dimnames(m), list(scores,
            scores)) && all(is.finite(m)) && identical(m, t(m)), logical(1))))
        stop("The G-study components or identities are incomplete or altered; recreate the result.",
            call. = FALSE)
    pair <- is.character(score) && length(score) == 2L
    if (pair) {
        if (!is.null(weights) || anyNA(score) || anyDuplicated(score) || !all(score %in% scores))
            stop("Select two distinct original scores, without weights, for covariance-component SEs.",
                call. = FALSE)
    }
    else {
        selection <- .mfrm_mvgt_interval_score(scores, score, weights)
        score <- selection$score
        weights <- selection$weights
    }
    k <- length(sources)
    if (identical(g$design$method, "Balanced multivariate ANOVA")) {
        design <- .mfrm_mvgt_anova_summary(g)
        for (m in g$mean_products) {
            unit <- max(abs(m))
            if (unit > 0 && min(eigen(m/unit, symmetric = TRUE, only.values = TRUE)$values) < -64 * .Machine$double.eps *
                nrow(m))
                stop("Stored ANOVA mean products must be positive semidefinite; recreate the G-study.",
                  call. = FALSE)
        }
        mapped <- design$A %*% do.call(rbind, lapply(g$mean_products, as.vector))
        original <- do.call(rbind, lapply(g$components, as.vector))
        if (!isTRUE(all.equal(unname(mapped), unname(original), tolerance = 1e-10)))
            stop("Stored components disagree with their mean products; recreate the G-study.", call. = FALSE)
    }
    else if (identical(g$design$method, "Multivariate MINQUE(0)")) {
        groups <- .mfrm_mvgt_moment_groups(g)
        design <- if (isTRUE(g$design$balanced) && isTRUE(g$design$complete))
            .mfrm_mvgt_anova_interval_design(g$design)
        else NULL
    }
    else stop("Unsupported G-study estimator.", call. = FALSE)
    if (!is.null(design)) {
        A <- design$A
        L <- solve(A)
        left <- vapply(seq_len(k), function(j) as.vector(tcrossprod(A[, j])), numeric(k * k))
        right <- t(vapply(seq_len(k), function(j) as.vector(tcrossprod(L[j, ])), numeric(k * k)))
        P <- sweep(left, 2L, 2/design$df, "*") %*% right
    }
    else P <- .mfrm_mvgt_moment_polynomial(groups, g$estimation$kernel_gram)
    map <- .mfrm_mvgt_moment_map(P, pair)
    if (pair) {
        v1 <- vapply(g$components, function(m) m[score[1L], score[1L]], numeric(1))
        v2 <- vapply(g$components, function(m) m[score[2L], score[2L]], numeric(1))
        estimate <- vapply(g$components, function(m) m[score[1L], score[2L]], numeric(1))
        unit1 <- max(abs(v1))
        unit2 <- max(abs(v2))
        if (unit1 == 0)
            unit1 <- 1
        if (unit2 == 0)
            unit2 <- 1
        unit <- sqrt(unit1) * sqrt(unit2)
        result <- .mfrm_mvgt_moment_correct(map, estimate/unit, v1/unit1, v2/unit2)
    }
    else {
        ws <- max(abs(weights))
        w <- weights/ws
        active <- which(w != 0)
        ms <- max(abs(unlist(lapply(g$components, function(m) m[active, active, drop = FALSE]))))
        if (ms == 0)
            ms <- 1
        theta <- vapply(g$components, function(m) drop(crossprod(w[active], (m[active, active, drop = FALSE]/ms) %*%
            w[active])), numeric(1))
        unit <- ms * ws * ws
        estimate <- theta * unit
        result <- .mfrm_mvgt_moment_correct(map, theta)
        if (!is.null(design)) {
            means <- if (identical(g$design$method, "Balanced multivariate ANOVA"))
                vapply(g$mean_products, function(m) drop(crossprod(w[active], (m[active, active, drop = FALSE]/ms) %*%
                  w[active])), numeric(1))
            else drop(L %*% theta)
            result$covariance <- tcrossprod(sweep(A, 2L, means * sqrt(2/(design$df + 2)), "*"))
            result$se <- sqrt(diag(result$covariance))
            result$estimable[] <- TRUE
        }
    }
    if (!is.finite(unit) || unit <= 0 || any(!is.finite(estimate)))
        stop("Sampling-covariance units exceed numeric range; rescale scores or weights.", call. = FALSE)
    out <- (result$covariance * unit) * unit
    if (any(is.infinite(out)) || any(!is.na(out) & result$covariance != 0 & out == 0))
        stop("Sampling covariance exceeds numeric range; rescale scores or weights.", call. = FALSE)
    dimnames(out) <- list(sources, sources)
    variance <- diag(out)
    se <- result$se * unit
    status <- ifelse(!diag(result$estimable), "Sampling variance not identifiable", ifelse(variance <
        0, "Negative sampling-variance estimate", ifelse(variance == 0, "Zero sampling-variance estimate",
        "Available")))
    attr(out, "summary") <- data.frame(Source = sources, Estimate = unname(estimate), SamplingVariance = unname(variance),
        SE = unname(se), SEStatus = unname(status))
    eigenvalues <- if (all(is.finite(result$covariance)))
        eigen(result$covariance, symmetric = TRUE, only.values = TRUE)$values
    else NA_real_
    attr(out, "diagnostics") <- list(moment_rank = map$rank, moment_dimension = map$dimension, singular_values = map$singular_values,
        estimable = result$estimable, reconstruction_error = map$error, positive_semidefinite = if (all(is.finite(eigenvalues))) min(eigenvalues) >=
            -sqrt(.Machine$double.eps) * max(abs(eigenvalues)) else NA, normalized_covariance = result$covariance,
        component_scale = unit)
    attr(out, "score") <- score
    attr(out, "weights") <- weights
    attr(out, "assumption") <- "normal"
    attr(out, "source_design") <- g$design
    attr(out, "source_rows") <- g$data_usage$counts
    attr(out, "notes") <- c("Sampling covariance of original ANOVA/MINQUE component estimates, conditional on the declared assignment under independent normal random effects. This is not measurement SEM, covariance between individual ratings or the SE of a separate REML estimator.",
        "Gaussian product-moment correction gives unbiased identifiable sampling-covariance entries under the model. Taking square roots is not an unbiasedness guarantee for SEs and does not justify normal-reference confidence intervals.",
        "Negative sampling variances are retained, with unavailable SEs. Unidentified entries are NA. An indefinite matrix is not repaired; inspect summary and diagnostics before propagating it to another target.",
        "Scores and composite weights are prespecified; omission must preserve the model conditional on assignment. No robustness, informative-missingness correction or general confidence-interval coverage is claimed. Save the complete result with saveRDS().")
    out
}

.mfrm_mvgt_moment_groups <- function(g) {
    ids <- c(g$design$person, unname(g$design$facets))
    data <- g$data
    if (anyDuplicated(ids) || !is.data.frame(data) || !all(ids %in% names(data)) || anyNA(data[ids]) ||
        nrow(data) != g$design$rows || anyDuplicated(data[ids]))
        stop("MINQUE sampling covariance needs the retained observed assignments; recreate the G-study.",
            call. = FALSE)
    factors <- lapply(data[ids], factor)
    actual_counts <- vapply(factors, nlevels, integer(1))
    subsets <- if (length(ids) == 2L)
        list(1L, 2L, 1:2)
    else list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
    if (!is.null(g$design$nesting)) {
        labels <- c("Person", names(g$design$facets))
        child <- match(names(g$design$nesting), labels)
        parent <- match(unname(g$design$nesting), labels)
        subsets <- list(1L, parent, c(child, parent), c(1L, parent), 1:3)
        children <- vapply(split(factors[[child]], factors[[parent]]), function(x) length(unique(x)),
            integer(1))
        actual_counts[child] <- sum(children)
    }
    potential <- if (is.null(g$design$nesting))
        prod(as.double(actual_counts))
    else actual_counts[1L] * actual_counts[child]
    complete <- unname(nrow(data) == potential)
    balanced <- complete && (is.null(g$design$nesting) || length(unique(children)) == 1L)
    if (!identical(as.numeric(g$design$counts), as.numeric(actual_counts)) || !identical(isTRUE(g$design$complete),
        complete) || !identical(isTRUE(g$design$balanced), balanced) || (!is.null(g$design$nesting) &&
        !identical(as.numeric(g$design$child_counts), as.numeric(children))))
        stop("Stored design counts or completeness do not match the observed assignments.", call. = FALSE)
    groups <- lapply(subsets, function(s) .mfrm_mvgt_group(factors[s]))
    k <- length(groups)
    N <- nrow(data)
    size <- lapply(groups, function(g) as.double(tabulate(g)[g]))
    total <- vapply(size, sum, numeric(1))
    S <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in seq_len(i)) {
        counts <- tabulate(.mfrm_mvgt_group(list(groups[[i]], groups[[j]])))
        S[i, j] <- S[j, i] <- sum(as.double(counts)^2) - 2 * sum(size[[i]] * size[[j]])/N + total[i] *
            total[j]/N^2
    }
    gram <- g$estimation$kernel_gram
    products <- g$estimation$quadratic_products
    if (!is.matrix(gram) || !identical(dim(gram), c(k, k)) || any(!is.finite(gram)) || !isTRUE(all.equal(unname(gram),
        S, tolerance = 1e-10)) || !is.list(products) || length(products) != k)
        stop("Stored MINQUE equations do not match the observed assignment; recreate the G-study.", call. = FALSE)
    rhs <- do.call(rbind, lapply(products, as.vector))
    theta <- do.call(rbind, lapply(g$components, as.vector))
    if (!isTRUE(all.equal(unname(gram %*% theta), unname(rhs), tolerance = 1e-10)))
        stop("Stored MINQUE components disagree with their quadratic products; recreate the G-study.",
            call. = FALSE)
    groups
}

.mfrm_mvgt_moment_polynomial <- function(groups, gram, block = 32L) {
    N <- length(groups[[1L]])
    k <- length(groups)
    scale <- sqrt(diag(gram))
    inverse <- solve(gram/outer(scale, scale))/outer(scale, scale)
    apply_K <- function(X, g) rowsum(X, g)[g, , drop = FALSE]
    apply_A <- function(X) {
        X <- sweep(X, 2, colMeans(X), "-")
        B <- vapply(groups, function(g) {
            x <- apply_K(X, g)
            as.vector(sweep(x, 2, colMeans(x), "-"))
        }, numeric(length(X)))
        B %*% t(inverse)
    }
    trace <- matrix(0, k * k, k * k)
    block <- max(1L, min(block, floor(1e+06/(N * k * k))))
    for (first in seq.int(1L, N, by = block)) {
        columns <- seq.int(first, min(N, first + block - 1L))
        E <- matrix(0, N, length(columns))
        E[cbind(columns, seq_along(columns))] <- 1
        right <- vapply(groups, function(g) as.vector(apply_A(apply_K(E, g))), numeric(length(E) * k))
        A <- matrix(apply_A(E), nrow = N)
        left <- vapply(groups, function(g) as.vector(apply_K(A, g)), numeric(length(E) * k))
        trace <- trace + crossprod(matrix(right, nrow = length(E)), matrix(left, nrow = length(E)))
    }
    P <- 2 * matrix(aperm(array(trace, rep(k, 4)), c(1, 3, 2, 4)), k * k, k * k)
    transposed <- as.vector(t(matrix(seq_len(k * k), k)))
    P <- (P + P[transposed, ])/2
    P <- (P + P[, transposed])/2
    P
}

.mfrm_mvgt_moment_symmetric <- function(P) {
    k <- as.integer(sqrt(nrow(P)))
    ij <- which(upper.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE)
    index <- ij[, 1] + (ij[, 2] - 1L) * k
    E <- matrix(0, k * k, nrow(ij))
    for (j in seq_len(nrow(ij))) {
        x <- matrix(0, k, k)
        x[ij[j, 1], ij[j, 2]] <- x[ij[j, 2], ij[j, 1]] <- 1
        E[, j] <- as.vector(x)
    }
    list(k = k, index = index, expand = E, L = P[index, , drop = FALSE] %*% E, diagonal = which(ij[,
        1] == ij[, 2]))
}

.mfrm_mvgt_moment_map <- function(P, pair = FALSE) {
    s <- .mfrm_mvgt_moment_symmetric(P)
    L <- s$L
    I <- diag(nrow(L))
    M <- if (pair)
        rbind(cbind(I, L), cbind(L/2, I + L/2))
    else I + L
    target <- if (pair)
        cbind(L/2, L/2)
    else L
    singular <- svd(M)
    keep <- singular$d > sqrt(.Machine$double.eps) * max(singular$d)
    inverse <- sweep(singular$v[, keep, drop = FALSE], 2, singular$d[keep], "/") %*% t(singular$u[, keep,
        drop = FALSE])
    weights <- target %*% inverse
    error <- sqrt(rowSums((weights %*% M - target)^2))/pmax(sqrt(rowSums(target^2)), .Machine$double.eps)
    list(symmetric = s, weights = weights, estimable = error < 1e-08, error = error, rank = sum(keep),
        dimension = nrow(M), singular_values = singular$d, moment = M, target = target, pair = pair)
}

.mfrm_mvgt_moment_correct <- function(map, theta, variance1 = NULL, variance2 = NULL) {
    s <- map$symmetric
    obs <- as.vector(tcrossprod(theta))[s$index]
    if (map$pair) {
        B <- (outer(variance1, variance2) + outer(variance2, variance1))/2
        obs <- c(B[s$index], obs)
    }
    corrected <- drop(map$weights %*% obs)
    raw <- matrix(s$expand %*% corrected, s$k)
    available <- matrix(s$expand %*% as.numeric(map$estimable), s$k) > 0
    value <- raw
    value[!available] <- NA_real_
    list(covariance = value, raw_covariance = raw, estimable = available, se = ifelse(is.finite(diag(value)) &
        diag(value) >= 0, sqrt(pmax(diag(value), 0)), NA_real_))
}
