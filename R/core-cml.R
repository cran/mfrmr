

mfrmr_cml_problem <- function(prep, idx, cfg, sizes) {
    adjacent <- mfrmr_estimability_adjacent_design(prep, idx, cfg, sizes, include_person = FALSE, include_population_beta = FALSE)
    p <- sum(unlist(sizes))
    K <- cfg$n_cat - 1L
    n <- nrow(prep$data)
    if (!p || !setequal(adjacent$map$OptimizerIndex, seq_len(p)))
        stop("CML requires free structural coordinates covered by the native model map.", call. = FALSE)
    order <- order(adjacent$map$OptimizerIndex)
    A <- adjacent$design[, order, drop = FALSE]
    keys <- do.call(paste, c(idx$facets, sep = ":"))
    cells <- unique(keys)
    features <- lapply(match(cells, keys), function(i) {
        increments <- as.matrix(A[i + (seq_len(K) - 1L) * n, , drop = FALSE])
        out <- matrix(0, K + 1L, p)
        for (k in seq_len(K)) out[k + 1L, ] <- out[k, ] + increments[k, ]
        out
    })
    cell <- match(keys, cells)
    persons <- split(seq_len(n), idx$person)
    groups <- list()
    group_index <- new.env(hash = TRUE, parent = emptyenv())
    observed <- numeric(p)
    totals <- vapply(persons, function(at) sum(idx$score_k[at]), numeric(1))
    maximum <- lengths(persons) * K
    extreme <- totals == 0 | totals == maximum
    for (i in which(!extreme)) {
        at <- persons[[i]]
        roster <- sort(cell[at])
        key <- paste(roster, collapse = ":")
        g <- group_index[[key]]
        if (is.null(g)) {
            g <- length(groups) + 1L
            group_index[[key]] <- g
            groups[[g]] <- list(features = features[roster], counts = numeric(length(at) * K + 1L))
        }
        groups[[g]]$counts[totals[i] + 1L] <- groups[[g]]$counts[totals[i] + 1L] + 1L
        for (j in at) observed <- observed + features[[cell[j]]][idx$score_k[j] + 1L, ]
    }
    list(groups = groups, observed = observed, parameters = p, map = adjacent$map[order, , drop = FALSE],
        person_totals = data.frame(Person = prep$levels$Person[as.integer(names(persons))], Total = totals,
            Maximum = maximum, Extreme = extreme, row.names = NULL))
}
