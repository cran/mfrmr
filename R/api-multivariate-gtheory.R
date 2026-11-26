mfrm_g_study <- function(data, scores = NULL, person = "Person", rater = "Rater", task = "Task", method = c("anova",
    "minque0", "urgenova"), missing = c("error", "omit"), facets = NULL, nesting = NULL, score = NULL, weights = NULL,
    d_study_grid = NULL, truncate_negative = TRUE) {
    if (is.null(scores) && !is.null(score))
        scores <- score
    method <- match.arg(method)
    missing <- match.arg(missing)
    if (!is.null(weights) && !identical(method, "urgenova")) {
        stop("`weights` in mfrm_g_study() requires method = 'urgenova'. For ANOVA or MINQUE(0), fit the separate scores first and supply composite weights to mfrm_d_study().",
            call. = FALSE)
    }
    if (identical(method, "urgenova")) {
        f_use <- if (!is.null(facets))
            unname(as.character(facets))
        else c(rater, task)[!vapply(list(rater, task), is.null, logical(1))]
        required <- c(person, f_use, scores)
        if (identical(missing, "error") && all(required %in% names(data)) && anyNA(data[required]))
            stop("Missing values are not allowed with missing = 'error'.", call. = FALSE)
        return(.mfrm_gstudy_henderson(data = data, person = person, facets = f_use, scores = scores, nesting = nesting,
            weights = weights, d_study_grid = d_study_grid, truncate_negative = truncate_negative))
    }
    valid_names <- function(x) is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(trimws(x))) &&
        !anyDuplicated(x)
    if (!is.data.frame(data) || nrow(data) == 0L || !valid_names(names(data))) {
        stop("`data` must be a nonempty data frame with unique, nonblank column names.", call. = FALSE)
    }
    if (is.null(facets)) {
        facet_args <- list(Rater = rater, Task = task)
        facet_args <- facet_args[!vapply(facet_args, is.null, logical(1))]
        facet_names <- names(facet_args)
        count_columns <- c(Rater = "Raters", Task = "Tasks")[facet_names]
    }
    else {
        if (!missing(rater) || !missing(task)) {
            stop("Use either `facets` or `rater`/`task`, not both.", call. = FALSE)
        }
        if (!is.character(facets) || !is.null(dim(facets)) || !valid_names(facets)) {
            stop("`facets` must specify one or two distinct identifier columns.", call. = FALSE)
        }
        facet_names <- if (is.null(names(facets)))
            facets
        else names(facets)
        if (!.mfrm_mvgt_valid_facets(facet_names)) {
            stop("Supply one or two unique, nonblank facet labels without ':' or reserved result names (such as Person, Residual, Score, or Scenario).",
                call. = FALSE)
        }
        facet_args <- as.list(unname(facets))
        count_columns <- setNames(facet_names, facet_names)
        rater <- task <- NULL
    }
    one_facet <- length(facet_args) == 1L
    nesting <- .mfrm_mvgt_nesting(nesting, facet_names)
    id_args <- c(list(person), facet_args)
    ids <- unlist(id_args, use.names = FALSE)
    if (!length(facet_args) || any(lengths(id_args) != 1L) || !all(vapply(id_args, is.character, logical(1))) ||
        !valid_names(ids) || !valid_names(scores) || any(scores %in% ids) || !all(c(ids, scores) %in% names(data))) {
        stop("Supply distinct person, facet, and score column names present in `data` (at least one facet).",
            call. = FALSE)
    }
    data <- as.data.frame(data[c(ids, scores)])
    groups <- lapply(data[ids], function(x) {
        if (!is.null(dim(x)) || !(is.factor(x) || (!is.object(x) && (is.character(x) || is.logical(x) || (is.numeric(x) &&
            !is.complex(x))))) || (is.numeric(x) && any(!is.na(x) & !is.finite(x))) || any(!is.na(x) & !nzchar(trimws(as.character(x))))) {
            stop("Factor identifiers must be finite, nonmissing, nonblank labels.", call. = FALSE)
        }
        factor(x)
    })
    if (!all(vapply(data[scores], function(x) is.numeric(x) && !is.object(x) && !is.complex(x) && is.null(dim(x)) &&
        all(is.na(x) | is.finite(x)), logical(1)))) {
        stop("All selected scores must be finite numeric values or NA; correct invalid scores explicitly.",
            call. = FALSE)
    }
    absent <- is.na(data)
    keep <- rowSums(absent) == 0L
    cells <- which(absent, arr.ind = TRUE)
    data_usage <- list(source = "data", missing = missing, counts = c(InputRows = nrow(data), UsedRows = sum(keep),
        ExcludedRows = sum(!keep)), excluded_rows = which(!keep), missing_cells = data.frame(InputRow = cells[,
        1L], Column = names(data)[cells[, 2L]], row.names = NULL))
    if (any(!keep) && missing == "error") {
        stop("Scores must be finite numeric values and identifiers nonmissing, nonblank; choose missing = 'omit' explicitly to exclude incomplete rows.",
            call. = FALSE)
    }
    if (!any(keep))
        stop("No complete rows remain for the G-study.", call. = FALSE)
    identified <- rowSums(absent[, ids, drop = FALSE]) == 0L
    if (anyDuplicated(data[identified, ids, drop = FALSE])) {
        stop("Duplicate cells in the selected design are not supported; do not average replicates silently.",
            call. = FALSE)
    }
    data <- data[keep, , drop = FALSE]
    groups <- lapply(groups, function(g) droplevels(g[keep]))
    factor_names <- c("Person", facet_names)
    counts <- setNames(vapply(groups, nlevels, integer(1)), factor_names)
    if (any(counts < 2L))
        stop("Each included factor needs at least two observed levels.", call. = FALSE)
    child_counts <- NULL
    potential_cells <- prod(as.double(counts))
    if (!is.null(nesting)) {
        child <- match(names(nesting), factor_names)
        parent <- match(unname(nesting), factor_names)
        pairs <- unique(data.frame(Parent = groups[[parent]], Child = groups[[child]]))
        child_counts <- setNames(as.integer(table(pairs$Parent)), levels(groups[[parent]]))
        counts[child] <- nrow(pairs)
        potential_cells <- as.double(counts[1L]) * sum(child_counts)
    }
    complete <- nrow(data) == potential_cells
    balanced <- complete && (is.null(nesting) || length(unique(child_counts)) == 1L)
    if (!balanced && method == "anova") {
        stop("ANOVA requires a complete, balanced design; consider method = 'minque0' for an incomplete or unequal observed design.",
            call. = FALSE)
    }
    if (!is.null(nesting) && method == "anova" && any(child_counts < 2L)) {
        stop("Cannot separate parent and nested-child components with one child per parent.", call. = FALSE)
    }
    y <- as.matrix(data[scores])
    y <- sweep(y, 2L, colMeans(y), "-")
    if (any(!is.finite(y)))
        stop("Score centering overflowed; rescale the scores.", call. = FALSE)
    score_scale <- sqrt(colSums(y^2)/(nrow(y) - 1))
    if (any(!is.finite(score_scale)) || any(score_scale == 0 & colSums(abs(y)) > 0)) {
        stop("Score variation exceeds numeric range; rescale the scores.", call. = FALSE)
    }
    names(score_scale) <- scores
    subsets <- if (one_facet)
        list(1L, 2L, 1:2)
    else list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
    if (!is.null(nesting))
        subsets <- list(1L, parent, c(child, parent), c(1L, parent), 1:3)
    sources <- .mfrm_mvgt_sources(factor_names, nesting)
    groupings <- setNames(lapply(subsets, function(s) .mfrm_mvgt_group(groups[s])), sources)
    estimation <- NULL
    if (method == "minque0") {
        fitted <- .mfrm_mvgt_minque0(y, groupings)
        components <- fitted$components
        estimation <- fitted$estimation
        mean_products <- degrees_of_freedom <- NULL
    }
    else {
        effects <- mean_products <- setNames(vector("list", length(sources)), sources)
        df <- setNames(numeric(length(sources)), sources)
        for (i in seq_along(subsets)) {
            subset <- subsets[[i]]
            group <- groupings[[i]]
            effect <- (rowsum(y, group)/tabulate(group))[group, , drop = FALSE]
            for (j in seq_len(i - 1L)) {
                if (all(subsets[[j]] %in% subset))
                  effect <- effect - effects[[j]]
            }
            effects[[i]] <- effect
            df[i] <- if (is.null(nesting))
                prod(counts[subset] - 1)
            else {
                P <- counts[1L]
                T <- counts[parent]
                R <- child_counts[1L]
                c(P - 1, T - 1, T * (R - 1), (P - 1) * (T - 1), T * (P - 1) * (R - 1))[i]
            }
            mean_products[[i]] <- crossprod(effect)/df[i]
        }
        m <- mean_products
        if (!is.null(nesting)) {
            components <- list((m[[1L]] - m[[4L]])/(R * T), (m[[2L]] - m[[3L]] - m[[4L]] + m[[5L]])/(P * R),
                (m[[3L]] - m[[5L]])/P, (m[[4L]] - m[[5L]])/R, m[[5L]])
        }
        else if (one_facet) {
            components <- list((m[[1L]] - m[[3L]])/counts[2L], (m[[2L]] - m[[3L]])/counts[1L], m[[3L]])
        }
        else {
            n1 <- counts[2L]
            n2 <- counts[3L]
            components <- list((m[[1L]] - m[[4L]] - m[[5L]] + m[[7L]])/(n1 * n2), (m[[2L]] - m[[4L]] - m[[6L]] +
                m[[7L]])/(counts[1L] * n2), (m[[3L]] - m[[5L]] - m[[6L]] + m[[7L]])/(counts[1L] * n1), (m[[4L]] -
                m[[7L]])/n2, (m[[5L]] - m[[7L]])/n1, (m[[6L]] - m[[7L]])/counts[1L], m[[7L]])
        }
        names(components) <- sources
        degrees_of_freedom <- data.frame(Source = sources, DF = unname(df))
    }
    if (any(!is.finite(unlist(components)))) {
        stop("Covariance estimation exceeded numeric range; rescale the scores.", call. = FALSE)
    }
    diagnostics <- .mfrm_mvgt_diagnostics(components, score_scale)
    .mfrm_new_gstudy(structure(list(components = components, mean_products = mean_products, degrees_of_freedom = degrees_of_freedom,
        estimation = estimation, data_usage = data_usage, component_diagnostics = diagnostics, score_scale = score_scale,
        design = list(person = person, rater = rater, task = task, scores = scores, facets = setNames(ids[-1L],
            facet_names), count_columns = count_columns, counts = counts, levels = lapply(groups, levels),
            rows = nrow(data), nesting = nesting, child_counts = child_counts, balanced = balanced, complete = complete,
            potential_cells = potential_cells, observed_fraction = nrow(data)/potential_cells, method = if (method ==
                "anova") "Balanced multivariate ANOVA" else "Multivariate MINQUE(0)", score_convention = if (!is.null(nesting)) paste0("Means over random ",
                unname(nesting), " conditions and their distinct random ", names(nesting), " conditions; ",
                count_columns[[names(nesting)]], " is the count per ", unname(nesting)) else if (is.null(facets)) paste("Means over common random",
                paste(tolower(count_columns), collapse = " and ")) else paste("Means over common random conditions:",
                paste(facet_names, collapse = ", ")), calculation_version = if (!is.null(nesting)) 3L else if (method ==
                "minque0") 2L else 1L), data = data), class = "mfrm_gtheory"), method, data)
}

.mfrm_mvgt_valid_facets <- function(x) {
    is.character(x) && length(x) %in% 1:2 && !anyNA(x) && all(nzchar(trimws(x))) && !anyDuplicated(x) &&
        !any(grepl(":", x, fixed = TRUE)) && !any(x %in% c("Person", "Residual", "Scenario", "Kind",
        "Score", "UniverseVariance", "RelativeErrorVariance", "AbsoluteErrorVariance", "G", "Phi", "RelativeSEM",
        "GStatus", "PhiStatus", "RelativeSEMStatus", "AbsoluteSEMStatus", "ComponentPSD", "AbsoluteSEM",
        "Status", "Metric", "X", "Group", "Value", "Panel"))
}

.mfrm_mvgt_nesting <- function(nesting, facets) {
    if (is.null(nesting))
        return(NULL)
    if (length(facets) != 2L || !is.character(nesting) || !is.null(dim(nesting)) || length(nesting) !=
        1L || is.na(nesting) || is.null(names(nesting)) || is.na(names(nesting)) || !names(nesting) %in%
        facets || !nesting %in% facets || names(nesting) == unname(nesting)) {
        stop("`nesting` must name one child facet and its other parent facet, for example c(Rater = 'Task'); nesting within Person is not supported.",
            call. = FALSE)
    }
    nesting
}

.mfrm_mvgt_sources <- function(factors, nesting = NULL) {
    if (!is.null(nesting))
        return(c("Person", unname(nesting), paste0(names(nesting), "(", nesting, ")"), paste0("Person:",
            nesting), "Residual"))
    c(factors, if (length(factors) == 3L) utils::combn(factors, 2L, paste, collapse = ":"), "Residual")
}

.mfrm_mvgt_group <- function(groups) {
    as.integer(interaction(lapply(groups, as.integer), drop = TRUE))
}

.mfrm_gt_kernel_gram <- function(groupings, centered = TRUE) {
    n <- as.double(length(groupings[[1L]]))
    k <- length(groupings)
    sizes <- lapply(groupings, tabulate)
    row_sizes <- Map(function(g, size) as.double(size[g]), groupings, sizes)
    totals <- vapply(row_sizes, sum, numeric(1))
    gram <- matrix(0, k, k, dimnames = list(names(groupings), names(groupings)))
    for (i in seq_len(k)) for (j in seq_len(i)) {
        joint_sizes <- tabulate(.mfrm_mvgt_group(list(groupings[[i]], groupings[[j]])))
        value <- sum(as.double(joint_sizes)^2)
        if (centered)
            value <- value - 2 * sum(row_sizes[[i]] * row_sizes[[j]])/n + totals[i] * totals[j]/n^2
        gram[i, j] <- gram[j, i] <- value
    }
    gram
}

.mfrm_mvgt_minque0 <- function(y, groupings) {
    gram <- .mfrm_gt_kernel_gram(groupings)
    k <- length(groupings)
    sizes <- lapply(groupings, tabulate)
    if (any(!is.finite(gram)) || any(diag(gram) <= 0)) {
        stop("The observed design cannot support the requested covariance components.", call. = FALSE)
    }
    scale <- sqrt(diag(gram))
    normalized <- gram/outer(scale, scale)
    eigenvalues <- eigen(normalized, symmetric = TRUE, only.values = TRUE)$values
    tolerance <- sqrt(.Machine$double.eps)
    rank <- sum(eigenvalues > tolerance * max(eigenvalues))
    if (rank < k) {
        stop("Cannot separate covariance components in this observed design (numerical rank ", rank,
            " of ", k, "). Review facet overlap and confounding; no estimates were returned.", call. = FALSE)
    }
    products <- lapply(groupings, function(g) crossprod(rowsum(y, g)))
    rhs <- do.call(rbind, lapply(products, as.vector))
    estimates <- solve(normalized, rhs/scale)/scale
    components <- setNames(lapply(seq_len(k), function(i) {
        matrix(estimates[i, ], ncol(y), dimnames = list(colnames(y), colnames(y)))
    }), names(groupings))
    list(components = components, estimation = list(kernel_gram = gram, quadratic_products = products,
        scaled_eigenvalues = eigenvalues, rank = rank, relative_tolerance = tolerance, condition_number = max(eigenvalues)/min(eigenvalues),
        component_support = data.frame(Source = names(groupings), Groups = lengths(sizes), MinimumRows = vapply(sizes,
            min, integer(1)), MaximumRows = vapply(sizes, max, integer(1)), row.names = NULL)))
}

.mfrm_mvgt_diagnostics <- function(components, score_scale) {
    scale <- score_scale
    scale[scale == 0] <- 1
    do.call(rbind, lapply(names(components), function(source) {
        standardized <- sweep(sweep(components[[source]], 1L, scale, "/"), 2L, scale, "/")
        values <- eigen(standardized, symmetric = TRUE, only.values = TRUE)$values
        tolerance <- sqrt(.Machine$double.eps) * max(1, abs(values))
        data.frame(Source = source, MinimumScaledEigenvalue = min(values), Tolerance = tolerance, PositiveSemidefinite = min(values) >=
            -tolerance, Rank = sum(values > tolerance), Dimension = length(values), row.names = NULL)
    }))
}

.mfrm_d_study_moments <- function(x, design_grid = NULL, weights = NULL, negative_components = c("raw",
    "genova_algorithm", "genova_ems"), ...) {
    rlang::check_dots_empty()
    negative_components <- match.arg(negative_components)
    if (!inherits(x, "mfrm_gtheory") || !(identical(x$design$calculation_version, 1L) ||
        identical(x$design$calculation_version, 2L) || identical(x$design$calculation_version, 3L))) {
        stop("`x` must be a current mfrm_gtheory result.", call. = FALSE)
    }
    scores <- x$design$scores
    counts <- x$design$counts
    facet_names <- names(counts)[-1L]
    count_columns <- x$design$count_columns
    if (is.null(count_columns) && all(facet_names %in% c("Rater", "Task"))) {
        count_columns <- c(Rater = "Raters", Task = "Tasks")[facet_names]
    }
    if (!is.numeric(counts) || !identical(names(counts)[1L], "Person") || !.mfrm_mvgt_valid_facets(facet_names) ||
        any(!is.finite(counts)) || any(counts < 2 | counts != floor(counts)) || !identical(names(count_columns),
        facet_names) || !.mfrm_mvgt_valid_facets(unname(count_columns))) {
        stop("The G-study factor or count identities are incomplete or altered; recreate the result.",
            call. = FALSE)
    }
    one_facet <- length(facet_names) == 1L
    nesting <- .mfrm_mvgt_nesting(x$design$nesting, facet_names)
    if (identical(x$design$calculation_version, 3L) != !is.null(nesting)) {
        stop("The G-study nesting specification is incomplete or altered; recreate the result.", call. = FALSE)
    }
    if (!is.null(nesting)) {
        child_counts <- x$design$child_counts
        if (!is.numeric(child_counts) || length(child_counts) != counts[[unname(nesting)]] || any(!is.finite(child_counts)) ||
            any(child_counts < 1 | child_counts != floor(child_counts)) || sum(child_counts) != counts[[names(nesting)]]) {
            stop("The G-study child counts are incomplete or altered; recreate the result.", call. = FALSE)
        }
    }
    sources <- .mfrm_mvgt_sources(names(counts), nesting)
    if (!is.character(scores) || !length(scores) || anyNA(scores) || anyDuplicated(scores) || !identical(names(x$components),
        sources) || !is.numeric(x$score_scale) || !identical(names(x$score_scale), scores) || any(!is.finite(x$score_scale)) ||
        any(x$score_scale < 0) || !all(vapply(x$components, function(a) is.matrix(a) && is.numeric(a) &&
        !is.complex(a) && identical(dim(a), rep(length(scores), 2L)) && identical(dimnames(a), list(scores,
        scores)) && all(is.finite(a)) && identical(a, t(a)), logical(1)))) {
        stop("The G-study covariance matrices or score identities are incomplete or altered; recreate the result.",
            call. = FALSE)
    }
    diagnostics <- .mfrm_mvgt_diagnostics(x$components, x$score_scale)
    grid_names <- unname(count_columns)
    if (is.null(design_grid)) {
        if (identical(x$design$complete, FALSE) || (!is.null(nesting) && length(unique(child_counts)) !=
            1L)) {
            stop("An incomplete or unequal G-study requires an explicit `design_grid` for a future complete balanced design; observed pool counts are not planned replication.",
                call. = FALSE)
        }
        future_counts <- setNames(counts[-1L], grid_names)
        if (!is.null(nesting))
            future_counts[count_columns[[names(nesting)]]] <- child_counts[1L]
        design_grid <- as.data.frame(as.list(future_counts), check.names = FALSE)
    }
    design_grid <- mfrm_d_study_count_columns(design_grid, facet_names, grid_names)
    if (!is.data.frame(design_grid) || !nrow(design_grid) || !setequal(names(design_grid), grid_names) ||
        ncol(design_grid) != length(grid_names) || !all(vapply(design_grid, function(n) is.numeric(n) &&
        !is.object(n) && is.null(dim(n)) && !is.complex(n) && all(is.finite(n)) && all(n >= 1 & n <=
        .Machine$integer.max & n == floor(n)), logical(1)))) {
        stop(paste0("`design_grid` must contain only positive integer ", paste(grid_names, collapse = " and "),
            " columns, one scenario per row."), call. = FALSE)
    }
    design_grid <- as.data.frame(design_grid[grid_names])
    rownames(design_grid) <- NULL
    vectors <- diag(length(scores))
    colnames(vectors) <- scores
    kind <- rep("Score", length(scores))
    if (!is.null(weights)) {
        if (is.matrix(weights)) {
            if (!is.numeric(weights) || is.object(weights) || is.complex(weights) || nrow(weights) !=
                length(scores) || !ncol(weights) || is.null(rownames(weights)) || anyNA(rownames(weights)) ||
                anyDuplicated(rownames(weights)) || !setequal(rownames(weights), scores) || is.null(colnames(weights)) ||
                anyNA(colnames(weights)) || anyDuplicated(colnames(weights)) || any(!nzchar(trimws(colnames(weights)))) ||
                any(!is.finite(weights)) || any(colSums(weights != 0) == 0)) {
                stop("A `weights` matrix must name every score once in its rows and have unique, nonempty composite names in its columns, with finite values and at least one nonzero weight per column.",
                  call. = FALSE)
            }
            weights <- weights[scores, , drop = FALSE]
            vectors <- cbind(vectors, weights)
            kind <- c(kind, rep("Composite", ncol(weights)))
        }
        else {
            if (!is.numeric(weights) || is.object(weights) || is.complex(weights) || !is.null(dim(weights)) ||
                length(weights) != length(scores) || is.null(names(weights)) || anyNA(names(weights)) ||
                anyDuplicated(names(weights)) || !setequal(names(weights), scores) || any(!is.finite(weights)) ||
                !any(weights != 0)) {
                stop("`weights` must name every score once with finite values and at least one nonzero value.",
                  call. = FALSE)
            }
            weights <- weights[scores]
            vectors <- cbind(vectors, Composite = weights)
            kind <- c(kind, "Composite")
        }
    }
    projection_scale <- apply(abs(vectors), 2L, max)
    if (any(!is.finite(projection_scale^2)) || any(projection_scale^2 == 0)) {
        stop("Composite variances exceed numeric range; rescale the weights.", call. = FALSE)
    }
    vectors <- sweep(vectors, 2L, projection_scale, "/")
    c <- .mfrm_mvds_components(x, negative_components)
    effective_diagnostics <- .mfrm_mvgt_diagnostics(c, x$score_scale)
    covariances <- vector("list", nrow(design_grid))
    tables <- vector("list", nrow(design_grid))
    for (i in seq_len(nrow(design_grid))) {
        n1 <- as.double(design_grid[[1L]][i])
        if (!is.null(nesting)) {
            n_child <- as.double(design_grid[[count_columns[[names(nesting)]]]][i])
            n_parent <- as.double(design_grid[[count_columns[[unname(nesting)]]]][i])
            relative <- c[[4L]]/n_parent + c[[5L]]/(n_child * n_parent)
            absolute <- relative + c[[2L]]/n_parent + c[[3L]]/(n_child * n_parent)
        }
        else if (one_facet) {
            relative <- c[[3L]]/n1
            absolute <- relative + c[[2L]]/n1
        }
        else {
            n2 <- as.double(design_grid[[2L]][i])
            relative <- c[[4L]]/n1 + c[[5L]]/n2 + c[[7L]]/(n1 * n2)
            absolute <- relative + c[[2L]]/n1 + c[[3L]]/n2 + c[[6L]]/(n1 * n2)
        }
        covariances[[i]] <- list(Universe = c$Person, RelativeError = relative, AbsoluteError = absolute)
        unit_variances <- matrix(vapply(covariances[[i]], function(a) colSums(vectors * (a %*% vectors)),
            numeric(ncol(vectors))), nrow = ncol(vectors), ncol = 3L)
        variances <- unit_variances * projection_scale^2
        if (any(!is.finite(variances)) || any(variances == 0 & unit_variances != 0)) {
            stop("D-study projection exceeded numeric range; rescale scores or weights.", call. = FALSE)
        }
        u <- unit_variances[, 1L]
        r <- unit_variances[, 2L]
        a <- unit_variances[, 3L]
        metrics <- c("G", "Phi", "RelativeSEM", "AbsoluteSEM")
        values <- matrix(NA_real_, length(u), 4L, dimnames = list(NULL, metrics))
        metric_status <- matrix("Available", length(u), 4L, dimnames = list(NULL, paste0(metrics, "Status")))
        for (j in 1:2) {
            error <- if (j == 1L)
                r
            else a
            metric_status[error < 0, c(j, j + 2L)] <- "Negative projected error variance"
            metric_status[u < 0, j] <- "Negative universe variance"
            metric_status[u == 0 & error == 0, j] <- "Zero total variance"
            ok <- metric_status[, j] == "Available"
            scale <- pmax(u[ok], error[ok])
            values[ok, j] <- (u[ok]/scale)/(u[ok]/scale + error[ok]/scale)
            ok_sem <- metric_status[, j + 2L] == "Available"
            values[ok_sem, j + 2L] <- sqrt(error[ok_sem]) * projection_scale[ok_sem]
        }
        available_count <- rowSums(metric_status == "Available")
        status <- ifelse(available_count == 4L, "Available", ifelse(available_count == 0L, "Unavailable",
            "Partially available"))
        tables[[i]] <- data.frame(Scenario = i, design_grid[i, , drop = FALSE], Kind = kind, Score = colnames(vectors),
            UniverseVariance = variances[, 1L], RelativeErrorVariance = variances[, 2L], AbsoluteErrorVariance = variances[,
                3L], values, Status = status, metric_status, ComponentPSD = all(diagnostics$PositiveSemidefinite),
            ComponentTreatment = negative_components, row.names = NULL, check.names = FALSE)
    }
    structure(list(coefficients = do.call(rbind, tables), covariances = covariances, design_grid = design_grid,
        weights = weights, component_diagnostics = diagnostics, negative_components = negative_components,
        effective_components = c, effective_component_diagnostics = effective_diagnostics, gstudy = x),
        class = "mfrm_multivariate_d_study")
}

.mfrm_gstudy_print_moments <- function(x, digits = NULL, row.names = FALSE, ...) {
    cat("Generalizability (G-Study): multivariate observed scores\n")
    label <- if (identical(x$design$complete, FALSE))
        "Incomplete observed design:"
    else if (!is.null(x$design$nesting))
        "Complete nested design:"
    else "Complete crossed design:"
    counts <- x$design$counts
    units <- c(Person = "persons", Rater = "raters", Task = "tasks")[names(counts)]
    units[is.na(units)] <- paste(names(counts)[is.na(units)], "levels")
    cat(label, paste(counts, units, collapse = ", "), "\n")
    if (!is.null(x$design$nesting)) {
        cat(names(x$design$nesting), "identities are local to", unname(x$design$nesting), "; the observed child count is the total across parents.\n")
        cat("Children per parent:", paste(names(x$design$child_counts), x$design$child_counts, sep = " = ",
            collapse = ", "), "\n")
    }
    cat(x$design$method, "; highest-order interaction and residual are combined.\n", sep = "")
    if (!is.null(x$data_usage)) {
        cat("Rows:", x$data_usage$counts[["UsedRows"]], "used of", x$data_usage$counts[["InputRows"]],
            ";", x$data_usage$counts[["ExcludedRows"]], "explicitly omitted.\n")
    }
    if (!is.null(x$estimation)) {
        cat("Covariance-component rank:", x$estimation$rank, "; scaled moment-system condition number:",
            format(x$estimation$condition_number, digits = 4), "\n")
    }
    print(x$component_diagnostics, digits = digits, row.names = row.names, ...)
    cat("Raw covariance components are retained; this table shows point estimates.\n")
    if (identical(x$design$method, "Balanced multivariate ANOVA"))
        cat("Use confint(x, assumption = 'normal') for eligible score-component SEs and intervals.\n")
    else cat("Use vcov(x, assumption = 'normal') for component sampling SEs; confint(x, method = 'reml_profile', assumption = 'normal') provides explicit experimental intervals.\n")
    invisible(x)
}


print.mfrm_multivariate_d_study <- function(x, digits = NULL, row.names = FALSE, ...) {
    if (!identical(.mfrm_mvds_policy(x), "raw"))
        .mfrm_mvds_replay(x)
    cat("Multivariate observed-score D-study\n")
    .mfrm_mvgt_print_source(x$gstudy$design, x$gstudy$data_usage$counts)
    cat("Future complete balanced plans; counts apply to every person and score.\n")
    .mfrm_mvgt_print_future_counts(x$gstudy$design)
    convention <- x$gstudy$design$score_convention
    if (is.null(convention))
        convention <- paste("Means over common random", paste(tolower(names(x$design_grid)), collapse = " and "))
    cat(paste0(convention, "; weights are used as supplied.\n"))
    cat(.mfrm_mvds_component_note(x), "\n")
    if (!is.null(x$weights)) {
        cat("Composite weights (not normalized):\n")
        print(if (is.matrix(x$weights))
            x$weights
        else matrix(x$weights, ncol = 1L, dimnames = list(names(x$weights), "Composite")), digits = digits)
    }
    if (any(!x$component_diagnostics$PositiveSemidefinite)) {
        cat("Warning: non-PSD covariance components in the original G-study; inspect component_diagnostics and the selected treatment.\n")
    }
    print(x$coefficients, digits = digits, row.names = row.names, ...)
    cat("Point projections conditional on estimated covariance components; no confidence intervals.\n")
    invisible(x)
}

.mfrm_mvgt_print_source <- function(design, rows) {
    if (is.null(design)) {
        cat("Source design not recorded in this saved comparison; consult the original G-study.\n")
        return(invisible(NULL))
    }
    cat("Source G-study:", design$method, "\n")
    cat("Observed level counts (not future replication):", paste(names(design$counts), design$counts,
        sep = " = ", collapse = ", "), "\n")
    if (!is.null(design$nesting))
        cat("The observed child count is the total across parents.\n")
    if (!is.null(rows)) {
        cat("Source rows:", rows[["UsedRows"]], "used of", rows[["InputRows"]], ";", rows[["ExcludedRows"]],
            "explicitly omitted.\n")
    }
    else if (!is.null(design$rows))
        cat("Source rows analyzed:", design$rows, "\n")
    if (identical(design$complete, FALSE) || identical(design$balanced, FALSE)) {
        cat("The observed source is incomplete or unequal; future counts do not reproduce its observed assignments.\n")
    }
    invisible(NULL)
}

.mfrm_mvgt_print_future_counts <- function(design) {
    if (!is.null(design$nesting)) {
        child <- names(design$nesting)
        cat("Future", design$count_columns[[child]], "counts levels per", unname(design$nesting), "; not the total across parents.\n")
    }
    invisible(NULL)
}
