.mfrm_mvgt_component_covariance <- function(g, theta) {
    n <- g$design$counts
    data <- g$data
    ids <- c(g$design$person, if (is.null(g$design$facets)) c(g$design$rater, g$design$task) else unname(g$design$facets))
    if (!is.data.frame(data) || length(ids) != 3L || !all(ids %in% names(data)) || nrow(data) != g$design$rows ||
        anyNA(data[ids]) || anyDuplicated(data[ids]) || !identical(unname(vapply(data[ids], function(x) length(unique(x)),
        integer(1))), unname(n))) {
        stop("The G-study must retain its analyzed data and matching facet identities; recreate the G-study.",
            call. = FALSE)
    }
    subsets <- list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
    if (isTRUE(g$design$complete)) {
        if (nrow(data) != prod(as.double(n)))
            stop("The stored complete design is inconsistent.", call. = FALSE)
        L <- matrix(0, 7, 7)
        for (i in 1:7) for (j in 1:7) if (all(subsets[[i]] %in% subsets[[j]])) {
            L[i, j] <- prod(n[setdiff(1:3, subsets[[j]])])
        }
        inverse <- solve(L)
        df <- vapply(subsets, function(s) prod(n[s] - 1), numeric(1))
        C <- inverse %*% diag(2 * drop(L %*% theta)^2/df) %*% t(inverse)
    }
    else {
        id_groups <- lapply(data[ids], factor)
        groups <- lapply(subsets, function(s) .mfrm_mvgt_group(id_groups[s]))
        S <- g$estimation$kernel_gram
        if (!is.matrix(S) || !identical(dim(S), c(7L, 7L)) || any(!is.finite(S))) {
            stop("The incomplete G-study needs its MINQUE covariance-component equations; recreate it.",
                call. = FALSE)
        }
        scale <- sqrt(diag(S))
        inverse <- solve(S/outer(scale, scale))/outer(scale, scale)
        apply_V <- function(X) Reduce(`+`, lapply(1:7, function(j) theta[j] * rowsum(X, groups[[j]])[groups[[j]],
            , drop = FALSE]))
        apply_A <- function(X) {
            X <- sweep(X, 2L, colMeans(X), "-")
            B <- vapply(groups, function(group) {
                KX <- rowsum(X, group)[group, , drop = FALSE]
                as.vector(sweep(KX, 2L, colMeans(KX), "-"))
            }, numeric(length(X)))
            B %*% t(inverse)
        }
        N <- nrow(data)
        block <- max(1L, min(32L, floor(1e+06/N)))
        C <- matrix(0, 7, 7)
        for (first in seq.int(1L, N, by = block)) {
            columns <- seq.int(first, min(N, first + block - 1L))
            E <- matrix(0, N, length(columns))
            E[cbind(columns, seq_along(columns))] <- 1
            AV <- apply_A(apply_V(E))
            VA <- apply_V(matrix(apply_A(E), nrow = N))
            C <- C + 2 * crossprod(AV, matrix(VA, nrow = length(E)))
        }
    }
    if (any(!is.finite(C)))
        stop("Sampling covariance exceeds numeric range; review the design and rescale scores.", call. = FALSE)
    (C + t(C))/2
}

summary.mfrm_multivariate_d_comparison <- function(object, ..., planning = FALSE) {
    rlang::check_dots_empty()
    if (!is.logical(planning) || length(planning) != 1L || is.na(planning))
        stop("`planning` must be TRUE or FALSE.", call. = FALSE)
    if (planning) {
        if (is.null(object$planning))
            stop("This result has no stored cost/target review. Use `planning = FALSE` to return the available plan comparisons.", call. = FALSE)
        review <- object$planning
        plans <- review$plans
        counts <- object$design_grid[plans$Scenario, , drop = FALSE]
        names(counts) <- tail(make.unique(c(names(plans), names(counts)), sep = "_count"), ncol(counts))
        review$plans[names(counts)] <- counts
        return(review)
    }
    object$comparisons
}

print.mfrm_multivariate_d_comparison <- function(x, digits = NULL, row.names = FALSE, ...) {
    cat("D-study differences from reference scenario", x$reference, "\n")
    .mfrm_mvgt_print_source(x$source_design, x$source_rows)
    cat(x$kind, ": ", x$score, "\n", sep = "")
    if (identical(x$kind, "Composite")) {
        cat("Composite weights (not normalized):\n")
        print(x$weights, digits = digits)
    }
    cat("Future complete crossed plans; counts apply to every person and score.\n")
    if (!is.null(x$source_design$score_convention))
        cat(x$source_design$score_convention, "\n")
    print(x$design_grid, digits = digits, row.names = row.names, ...)
    print(x$comparisons, digits = digits, row.names = row.names, ...)
    cat(format(100 * x$level), "% approximate pointwise intervals; normal random effects required.\n",
        sep = "")
    cat("Positive G/Phi or negative SEM differences favor the comparison plan.\n")
    cat("An interval containing zero does not establish equivalence.\n")
    if (any(!x$component_diagnostics$PositiveSemidefinite) || !x$sampling_covariance_psd) {
        cat("Non-PSD covariance estimates: inspect component and sampling-covariance diagnostics.\n")
    }
    if (!is.null(x$planning)) {
        cat("\nResource feasibility and point targets (all plans retained):\n")
        print(summary(x, planning = TRUE)$plans, digits = digits, row.names = row.names, ...)
        print(x$planning$costs, digits = digits, row.names = row.names, ...)
        print(x$planning$differences[c("Scenario", "Metric", "TargetStatus", "Direction")], row.names = row.names)
        cat(paste(x$planning$notes, collapse = "\n"), "\n")
    }
    invisible(x)
}
